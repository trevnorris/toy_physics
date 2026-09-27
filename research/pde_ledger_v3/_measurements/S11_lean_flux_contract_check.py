#!/usr/bin/env python3
"""Sequential isolated F1–F4 verification; invoke inside the host resource guard.

No reuse is needed for these four small modules. Existing proof objects are never
on LEAN_PATH. The pinned external cache is reused and explicitly recorded.
"""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import sys
import time

BASE = Path(__file__).resolve().parents[1]
LEAN = BASE / 'lean'
sys.path.insert(0, str(LEAN))
import verify

ROOT = 'S11ScatteringFlux'
M = BASE / '_measurements'


def control_sources():
    cases = [
        ('interference', 'flux coherent both', '4', '2', 'coherent_flux'),
        ('conjugation', 'flux oneCurrent phaseAmplitude', '1', '-1', 'phase_flux'),
        ('basis_metric', 'flux (pullback oneCurrent doubleBasis) realAmplitude', '4', '1', 'basis_flux'),
        ('orientation', 'outgoing 3 1', '-2', '4', 'oriented_outgoing_witness'),
        ('signed_incident', 'incident .left (-3)', '-3', '3', 'signed_incident_witness'),
        ('zero_denominator', 'fraction 7 0', 'none', 'some 0', 'zero_denominator_witness'),
        ('normalization', 'fraction 2 4', 'some (1 / 2)', 'some 2', 'positive_fraction_witness'),
        ('negative_domain', 'fraction 2 (-4)', 'some (-1 / 2)', 'some (1 / 2)', 'negative_fraction_witness'),
        ('no_unconditional_unit_bound', 'fraction 4 2', 'some 2', 'some 1', 'fraction_above_one_witness'),
        ('no_unconditional_conservation', '(1 : ℝ) / 4 + 1 / 4', '(1 / 2)', '1', ''),
        ('current_reality_premise', '(pair (fun _ _ : Fin 1 => Complex.I) realAmplitude realAmplitude).im',
         '1', '0', 'nonhermitian_imaginary_witness'),
    ]
    prefix = 'import S11ScatteringFlux\nopen S11ScatteringFlux Matrix\n'
    result = []
    for name, lhs, good, bad, lemma in cases:
        for suffix, rhs, expected in [('positive', good, 'PASS'), ('mutant', bad, 'REJECTED')]:
            # norm_num does not run the default constructor simprocs. Supply
            # proved Option facts, retaining the literal-False acceptance rule.
            rules = lemma
            if name == 'zero_denominator':
                rules += ', Ne.symm (Option.some_ne_none (0 : ℝ))'
            elif lhs.startswith('fraction '):
                rules += ', Option.some.injEq'
            source = prefix + f'theorem contract_control : {lhs} = {rhs} := by\n  norm_num only [{rules}]\n'
            result.append({'name': name + '_' + suffix, 'source': source, 'expected': expected,
                           'required_failure_in': 'contract_control' if expected == 'REJECTED' else []})
    extra = {
        'empty_positive': 'example (J : Matrix (Fin 0) (Fin 0) ℂ) (x : Fin 0 → ℂ) : flux J x = 0 := empty_flux J x',
        'indefinite_null_positive': 'example : flux (!![1,0;0,-1] : Matrix (Fin 2) (Fin 2) ℂ) both = 0 ∧ both ≠ 0 := null_flux_nonzero_amplitude',
        'balance_positive': 'example : (1 : ℝ) / 4 + 1 / 4 = 1 - 2 / 4 := conditional_balance (by norm_num) (by norm_num)',
        'zero_defect_positive': 'example {c s d : ℝ} (h : c + s + 0 = d) (hd : d ≠ 0) : c / d + s / d = 1 := by\n  simpa using conditional_balance h hd',
    }
    result += [{'name': n, 'source': prefix + s + '\n', 'expected': 'PASS'} for n, s in extra.items()]
    return result


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--run-directory', type=Path, required=True)
    args = parser.parse_args()
    run = args.run_directory.resolve()
    run.relative_to(BASE / '_scratch/S11_lean_flux')
    run.mkdir(parents=True, exist_ok=False)
    logs, snapshot, objects = run / 'logs', run / 'snapshot', run / 'objects'
    for p in (logs, snapshot, objects):
        p.mkdir()
    result = {'status': 'RUNNING', 'checks': [], 'resources': {'lean_workers': 1,
              'lean_allocator_mib': 2048, 'per_process_seconds': 180,
              'outer_guard': '2 GiB/no swap/one CPU; receipt in host guard directory'},
              'limits': ['No canonical source-replacement mutations.',
                         'Fresh local objects; pinned external transitive cache baseline.',
                         'No physical scattering/current-conservation certificate.']}
    report = M / 'S11_lean_flux_contract_checks.json'
    def save():
        text = json.dumps(result, indent=2) + '\n'
        (run / 'report.json').write_text(text)
        temp = report.with_suffix('.new')
        temp.write_text(text)
        temp.replace(report)
    save()
    try:
        env, executable, dependencies = verify.external_environment()
        order = verify.build_order([ROOT])
        files = {verify.module_path(m) for m in order}
        files |= {Path('lean') / n for n in ['lean-toolchain', 'lake-manifest.json', 'lakefile.toml', 'verify.py']}
        files |= {Path('_measurements') / n for n in ['S11_lean_flux_contract_check.py', 'S11_lean_flux_source_check.py']}
        direct = {}
        for module in order:
            source = (BASE / verify.module_path(module)).read_text()
            for child in verify.IMPORT.findall(source):
                if verify.LOCAL.fullmatch(child):
                    continue
                stem = child.replace('.', '/')
                p = next((Path(d) / (stem + '.olean') for d in env['LEAN_PATH'].split(os.pathsep)
                          if (Path(d) / (stem + '.olean')).is_file()), None)
                if p is None:
                    raise ValueError('Missing external import: ' + child)
                src = LEAN / '.lake/packages/mathlib' / (stem + '.lean')
                direct[child] = {'object_path': str(p), 'object_sha256': verify.sha(p),
                                 'source_path': str(src), 'source_sha256': verify.sha(src)}
        for rel in files:
            target = snapshot / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(BASE / rel, target)
        before = {str(p): verify.sha(snapshot / p) for p in files}
        result.update(source_sha256=before, dependencies=dependencies, direct_imports=direct,
                      order=order, run_directory=str(run))
        env['LEAN_PATH'] = str(objects) + os.pathsep + env['LEAN_PATH']
        env.pop('PYTHONOPTIMIZE', None)
        env['PYTHONDONTWRITEBYTECODE'] = '1'
        save()
        def execute(name, cmd, source=None, expected='PASS', required=None, obj=None):
            log = logs / (name + '.log')
            started = time.monotonic()
            code = verify.run_process(cmd, snapshot / 'lean', log, 180, env)
            output = log.read_text()
            valid = (code == 0 and not re.search(r'\b(?:warning|error)(?:\([^)]*\))?:', output)) if expected == 'PASS' else verify.adjudicate(code, output, source, required)
            row = {'name': name, 'command': cmd, 'expected': expected, 'exit_status': code,
                   'passed': valid, 'seconds': round(time.monotonic() - started, 3),
                   'log': str(log), 'log_sha256': verify.sha(log), 'output': output}
            if source is not None:
                row.update(source=source, source_sha256=__import__('hashlib').sha256(source.encode()).hexdigest(),
                           required_failure_in=required or [])
            if obj is not None:
                row.update(object_path=str(obj), object_sha256=verify.sha(obj) if obj.exists() else None)
                valid = row['passed'] = valid and obj.is_file()
            result['checks'].append(row)
            save()
            if not valid:
                raise RuntimeError('Verification failed: ' + name)
            return output
        # Execute only the selected original native algebra, never production.
        execute('native', [sys.executable, str(M / 'S11_lean_flux_source_check.py')])
        native_path = M / 'S11_lean_flux_source_checks.json'
        native = verify.load(native_path)
        if native['status'] != 'PASS':
            raise ValueError('Native check did not pass')
        result['native_report_sha256'] = verify.sha(native_path)
        command = [executable, '-j1', '-M2048', '-DwarningAsError=true']
        audits = 0
        for module in order:
            path = snapshot / verify.module_path(module)
            obj = objects / (module.replace('.', '/') + '.olean')
            obj.parent.mkdir(parents=True, exist_ok=True)
            out = execute('build_' + module, command + ['-o', str(obj), str(path)], obj=obj)
            count = verify.audit_axioms(out, path.read_text())
            audits += count
            result['checks'][-1]['axiom_audits'] = count
        for c in control_sources():
            path = snapshot / 'lean' / (c['name'] + '.lean')
            path.write_text(c['source'])
            execute(c['name'], command + [str(path)], c['source'], c['expected'], c.get('required_failure_in'))
        if any(verify.sha(BASE / p) != h or verify.sha(snapshot / p) != h for p, h in before.items()):
            raise ValueError('Local source drift')
        for item in direct.values():
            if verify.sha(item['object_path']) != item['object_sha256'] or verify.sha(item['source_path']) != item['source_sha256']:
                raise ValueError('External dependency drift')
        _, _, after = verify.external_environment()
        if after != dependencies or verify.sha(native_path) != result['native_report_sha256']:
            raise ValueError('Pin/native report drift')
        if any(verify.sha(BASE / p) != h for p, h in native['source_sha256'].items()):
            raise ValueError('Native source drift')
        result.update(status='PASS', axiom_audits=audits,
                      objects={str(p.relative_to(objects)): verify.sha(p) for p in objects.rglob('*.olean')})
    except BaseException as error:
        result.update(status='ERROR', error=repr(error))
        (run / 'issues.log').write_text(repr(error) + '\n')
        raise
    finally:
        save()


if __name__ == '__main__':
    main()
