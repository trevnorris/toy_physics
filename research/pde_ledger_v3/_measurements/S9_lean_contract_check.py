#!/usr/bin/env python3
"""Sequential, bounded S9 verification and isolated mathematical mutations.

One Lean worker, 4096 MiB Lean allocator limit, low scheduling priority. Logs,
sources and incremental records survive an error. A timeout/resource failure is
an instrument failure, never a successfully rejected mutation. No CAS jobs.
"""
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

BASE = Path(__file__).resolve().parents[1]
LEAN = BASE / 'lean'
SCRATCH = LEAN / 's9/_scratch/closure'
REPORT = BASE / '_measurements/S9_lean_contract_checks.json'
ORDER = ('Action', 'PlaneWave', 'Analytic', 'Variation', 'FiniteAction',
         'Spectrum', 'Certificate', 'VariationalCertificate', 'Madelung')
CASES = (
    ('action_shear_sign', 'Action',
     'rho / 2 * normSq (J 0) - mu / 2 * normSq (jetCurl J)',
     'rho / 2 * normSq (J 0) + mu / 2 * normSq (jetCurl J)',
     {'lagrangian_variation'}),
    ('relative_action_expansion_sign', 'Variation',
     'relativeAction rho mu u h s =\n      s * (∫ x, linearDensity rho mu (fieldJet u x) (fieldJet h x)) +\n      s ^ 2 * (∫ x, lagrangian rho mu (fieldJet h x))',
     'relativeAction rho mu u h s =\n      -(s * (∫ x, linearDensity rho mu (fieldJet u x) (fieldJet h x))) -\n      s ^ 2 * (∫ x, lagrangian rho mu (fieldJet h x))',
     {'relativeAction_expansion'}),
    ('phase_wrong_coordinate', 'Madelung',
     '(hbar / mass) * coordDeriv i.succ theta x',
     '(hbar / mass) * coordDeriv 0 theta x', {'velocity_phasePerturbation'}),
    ('phase_sine_sign', 'Madelung',
     '(-(hbar / mass) * A) • k',
     '((hbar / mass) * A) • k', {'linearVelocity_eq'}),
)
EXAMPLES = (
    ('transverse_count', '''theorem domain_control :
    Module.finrank ℝ (transverseSpace (![0, 0, 1] : Vec)) = REL := by
  have hk : (![0, 0, 1] : Vec) ≠ 0 := by
    intro h
    have h2 := congrFun h 2
    norm_num [Matrix.cons_val_two] at h2
  have hd := transverseSpace_finrank hk
  norm_num [hd]
''', '2', '3'),
    ('longitudinal_count', '''theorem domain_control :
    Module.finrank ℝ (longitudinalSpace (![0, 0, 1] : Vec)) = REL := by
  have hk : (![0, 0, 1] : Vec) ≠ 0 := by
    intro h
    have h2 := congrFun h 2
    norm_num [Matrix.cons_val_two] at h2
  have hd := longitudinalSpace_finrank hk
  norm_num [hd]
''', '1', '2'),
    ('zero_wavevector_domain', '''theorem domain_control :
    ModalStationary 1 1 0 (0 : Vec) (![0, 0, 1] : Vec) ∧
      (![0, 0, 1] : Vec) REL longitudinalSpace (0 : Vec) := by
  constructor
  · rw [modal_stationary_iff]
    ext i
    simp [modalOperator, normSq, dot]
  · norm_num [longitudinalSpace]
''', '∉', '∈'),
    ('longitudinal_not_transverse', '''theorem domain_control :
    Madelung.cosAmplitude 1 1 1 (![0, 0, 1] : Vec) REL
      transverseSpace (![0, 0, 1] : Vec) := by
  norm_num [Madelung.cosAmplitude, transverseSpace, dotLinear, dot, Matrix.cons_val_two]
''', '∉', '∈'),
    ('nonzero_prefactor_required_for_range', '''theorem domain_control :
    Madelung.cosAmplitude 0 1 1 (![0, 0, 1] : Vec) REL 0 := by
  norm_num [Madelung.cosAmplitude]
''', '=', '≠'),
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--reuse-build', action='store_true',
                        help='Reuse successful builds only when all canonical source hashes match.')
    args = parser.parse_args()
    SCRATCH.mkdir(parents=True, exist_ok=True)
    sources = sorted((LEAN / 's9/S9Pilot').glob('*.lean')) + [LEAN / 's9/S9Pilot.lean']
    before = {str(p.relative_to(BASE)): digest(p) for p in sources}
    prior = json.loads(REPORT.read_text()) if args.reuse_build and REPORT.exists() else {}
    reusable = {r['name']: r for r in prior.get('checks', []) if r['outcome'] == 'PASS'} if (
        prior.get('canonical_sha256_before') == before) else {}
    report = {'status': 'RUNNING', 'started_utc': datetime.now(timezone.utc).isoformat(),
              'instrument_sha256': digest(Path(__file__)),
              'resources': {'concurrent_lean_processes': 1, 'threads': 1,
                            'lean_allocator_limit_mib': 4096, 'process_timeout_seconds': 300,
                            'limits_note': 'Lean allocator limit is not an OS RSS cap; one process at a time.'},
              'canonical_sha256_before': before, 'checks': []}

    def save():
        temporary = REPORT.with_suffix('.new')
        temporary.write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
        temporary.replace(REPORT)

    def run(name, path, expected, required=(), output_module=None):
        rel = str(path.relative_to(LEAN))
        command = ['nice', '-n', '10', 'lake', 'env', 'lean', '-j1', '-M4096',
                   '-DwarningAsError=true']
        if output_module:
            command += ['-o', str(LEAN / '.lake/build/lib/lean' / output_module)]
        command.append(rel)
        if output_module and name in reusable:
            old = reusable[name]
            if old['source_sha256'] == digest(path) and old['command'] == command:
                record = {**old, 'reused_build_from_instrument_sha256': prior['instrument_sha256']}
                report['checks'].append(record)
                save()
                return record
        start = time.monotonic()
        log = SCRATCH / (name + '.log')
        with log.open('w') as stream:
            result = subprocess.run(command, cwd=LEAN, stdout=stream, stderr=subprocess.STDOUT,
                                    env={**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'},
                                    timeout=300)
        output = log.read_text()
        source = path.read_text()
        decls = [(source[:m.start()].count('\n') + 1, m[1]) for m in
                 re.finditer(r'^\s*(?:theorem|def|abbrev)\s+(\w+)', source, re.M)]
        errors = list(re.finditer(re.escape(rel) + r':(\d+):\d+: error:', output))
        diagnostics = []
        for i, match in enumerate(errors):
            line = int(match[1])
            decl = next((n for pos, n in reversed(decls) if pos <= line), '')
            block = output[match.start():errors[i+1].start() if i+1 < len(errors) else len(output)]
            diagnostics.append({'declaration': decl, 'line': line, 'message': block})
        invalid = re.search(r'unknown (?:module|namespace|identifier)|unexpected token|'
                            r'maximum (?:recursion|heartbeats)|excessive memory|out of memory|'
                            r'PANIC|failed to create|No such file', output, re.I)
        mathematical = r'unsolved goals|Type mismatch|failed to close|proved that the proposition|Tactic `rewrite` failed'
        intended = [d for d in diagnostics if d['declaration'] in required and
                    re.search(mathematical, d['message'], re.I)]
        valid = (result.returncode == 0 and not re.search(r'\b(?:error|warning):', output)) if expected == 'PASS' else (
            result.returncode == 1 and bool(intended) and not invalid)
        record = {'name': name, 'expected': expected, 'outcome': expected if valid else 'UNEXPECTED',
                  'command': command, 'source_sha256': digest(path), 'exit_status': result.returncode,
                  'wall_seconds': round(time.monotonic() - start, 3),
                  'required_failure_in': sorted(required), 'diagnostics': diagnostics,
                  'output': output, 'output_sha256': digest(log)}
        report['checks'].append(record)
        save()
        if not valid and expected == 'PASS':
            raise RuntimeError(f'{name}: unexpected result; inspect {log}')
        return record

    save()
    try:
        # Rebuild only the small S9 library sequentially; reuse pinned Mathlib/Physlib.
        for module in ORDER:
            run('build_' + module, LEAN / 's9/S9Pilot' / (module + '.lean'), 'PASS',
                output_module='S9Pilot/' + module + '.olean')
        audit = run('axiom_audit', LEAN / 's9/S9Pilot.lean', 'PASS', output_module='S9Pilot.olean')
        axioms = re.findall(r"depends on axioms: \[([^\]]*)\]", audit['output'])
        assert len(axioms) == 32, ('axiom audit count', len(axioms))
        assert all(set(a.replace(' ', '').split(',')) <= {'propext', 'Classical.choice', 'Quot.sound', ''}
                   for a in axioms)
        assert not any(re.search(r'^\s*(?:axiom\b|.*\b(?:sorry|admit)\b)', p.read_text(), re.M) for p in sources)
        report['axiom_audit_count'] = len(axioms)
        for name, module, old, new, required in CASES:
            canonical = LEAN / 's9/S9Pilot' / (module + '.lean')
            source = canonical.read_text()
            assert source.count(old) == 1, (name, source.count(old))
            mutant = SCRATCH / (name + '.lean')
            mutant.write_text(source.replace(old, new))
            record = run(name, mutant, 'REJECTED', required)
            record.update(canonical_source=str(canonical.relative_to(BASE)),
                          canonical_sha256=digest(canonical), replacement={'old': old, 'new': new})
            save()
        for name, statement, correct, wrong in EXAMPLES:
            for label, relation, expected in [('positive', correct, 'PASS'), ('mutant', wrong, 'REJECTED')]:
                path = SCRATCH / (name + '_' + label + '.lean')
                path.write_text('import S9Pilot.Madelung\nopen S9Pilot\nnoncomputable section\n' + statement.replace('REL', relation))
                record = run(name + '_' + label, path, expected, {'domain_control'})
                record['source'] = path.read_text()
                save()
        after = {str(p.relative_to(BASE)): digest(p) for p in sources}
        assert before == after, 'canonical source changed during verification'
        assert all(c['outcome'] != 'UNEXPECTED' for c in report['checks']), 'unexpected mutation outcomes; inspect recorded diagnostics'
        report.update(status='PASS', canonical_sha256_after=after)
    except Exception as error:
        report.update(status='ERROR', error=repr(error))
        raise
    finally:
        report['finished_utc'] = datetime.now(timezone.utc).isoformat()
        save()


if __name__ == '__main__':
    main()
