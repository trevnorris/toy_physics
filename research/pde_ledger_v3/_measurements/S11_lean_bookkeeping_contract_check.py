#!/usr/bin/env python3
"""Sequential isolated P1–P4 verification; invoke inside the host resource guard.

No reuse is needed for these seven local modules including the unchanged Current import. Existing proof objects are never
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

ROOT = 'S11ScatteringBookkeeping'
M = BASE / '_measurements'


def control_sources():
    second='q2 (scalarB 1) (scalarB 4) (scalarB 5) (scalarA 1) (scalarA 2) (scalarA 3)'
    cases=[
      ('path_ratio','rectangle (scalarA 1) (scalarA 2) (scalarA 3) (scalarA 4) 1 2 0','17','25','path_witness'),
      ('delta_zero_jet','delta (scalarA 2) (scalarA 3) (scalarA 4) 1 2 0','16','14','delta_witness'),
      ('first_interference','q1 (scalarB 1) (scalarB 4) (scalarB 5) (scalarA 1) (scalarA 2) (scalarA 3)','8','4','first_witness'),
      ('current_variation',second,'31','10','second_witness'),
      ('second_amplitude',second,'31','25','second_witness'),
      ('baseline_interference',second,'31','9','second_witness'),
      ('higher_terms','q3 (scalarB 1) (scalarB 4) (scalarB 5) (scalarA 1) (scalarA 2) (scalarA 3)','72','0','higher_witness'),
      ('induced_leading','q2 (scalarB 1) (scalarB 4) (scalarB 5) 0 (scalarA 2) (scalarA 3)','4','31','induced_witness'),
      ('subtracted_is_not_induced','flux (scalarB 1) (scalarA 1 + scalarA 1) - flux (scalarB 1) (scalarA 1)','3','1','subtracted_witness'),
      ('omitted_parent_second','q2 (scalarB 1) 0 0 (scalarA 1) 0 (scalarA 3)','6','0','parent_second_witness'),
      ('epsilon_squared','flux (scalarB 1) ((2 : ℂ) • scalarA 1)','4','2','epsilon_witness'),
      ('leading_normalization','c0 (2 : ℝ) 2','1','2','c0'),
      ('incident_first','c1 (2 : ℝ) 3 2 1','1','3/2','c1, c0'),
      ('incident_second','c2 (2 : ℝ) 3 5 2 1 2','1','2','c2, c1, c0'),
      ('negative_leading','c0 (2 : ℝ) (-2)','-1','1','c0'),
      ('zero_model_two_constants','((0 : ℝ) * 0 = 0 ∧ (0 : ℝ) * 1 = 0)','True','False','true_and, eq_iff_iff, iff_self, true_iff_false'),
    ]
    prefix='import S11ScatteringBookkeeping\nopen S11ScatteringBookkeeping S11ScatteringFlux Matrix\n'
    out=[]
    for name,lhs,good,bad,lemma in cases:
      for suffix,rhs,expected in [('positive',good,'PASS'),('mutant',bad,'REJECTED')]:
        source=prefix+f'theorem contract_control : {lhs} = {rhs} := by\n  norm_num only [{lemma}]\n'
        out.append({'name':name+'_'+suffix,'source':source,'expected':expected,'required_failure_in':'contract_control' if expected=='REJECTED' else []})
    extra={
      'actual_quotient_equations_positive':'example : (2:ℂ) * c0 2 2 = 2 ∧ (2:ℂ)*c1 2 3 2 1 + 1*c0 2 2 = 3 ∧ (2:ℂ)*c2 2 3 5 2 1 2 + 1*c1 2 3 2 1 + 2*c0 2 2 = 5 := quotient_equations 2 3 5 2 1 2 (by norm_num)',
      'zero_obstruction_positive':'example (u : ℝ) : (0 : ℝ) * u ≠ 1 := zero_leading_obstruction 1 (by norm_num) u',
      'zero_nonunique_positive':'example : (0 : ℝ) * 0 = 0 ∧ (0 : ℝ) * 1 = 0 ∧ (0 : ℝ) ≠ 1 := zero_model_two_solutions',
      'empty_amplitude_positive':'example (B : Matrix (Fin 0) (Fin 0) ℂ) (a : Fin 0 → ℂ) : flux B a = 0 := by simp [flux, pair, dotProduct]',
    }
    return out+[{'name':name,'source':prefix+src+'\n','expected':'PASS'} for name,src in extra.items()]


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--run-directory', type=Path, required=True)
    parser.add_argument('--native-evidence', type=Path)
    args = parser.parse_args()
    run = args.run_directory.resolve()
    run.relative_to(BASE / '_scratch/S11_lean_bookkeeping')
    run.mkdir(parents=True, exist_ok=False)
    logs, snapshot, objects = run / 'logs', run / 'snapshot', run / 'objects'
    for p in (logs, snapshot, objects):
        p.mkdir()
    result = {'status': 'RUNNING', 'checks': [], 'resources': {'lean_workers': 1,
              'lean_allocator_mib': 2048, 'per_process_seconds': 180,
              'outer_guard': '2 GiB/no swap/one CPU; receipt in host guard directory'},
              'limits': ['No canonical source-replacement mutations.',
                         'Fresh local objects; pinned external transitive cache baseline.',
                         'No physical scattering or parent-theory Taylor certificate.']}
    report = M / 'S11_lean_bookkeeping_contract_checks.json'
    def save():
        text = json.dumps(result, indent=2) + '\n'
        (run / 'report.json').write_text(text)
        temp = report.with_suffix('.new')
        temp.write_text(text)
        temp.replace(report)
    save()
    try:
        preserve = verify.load(M / 'S11_lean_bookkeeping_preserved_inputs.json')
        def preserved():
            for key in ['source_sha256','shared_olean_sha256']:
                for rel,h in preserve[key].items():
                    if verify.sha(BASE/rel)!=h:raise ValueError('Historical preservation failure: '+rel)
        preserved()
        env, executable, dependencies = verify.external_environment()
        order = verify.build_order([ROOT])
        files = {verify.module_path(m) for m in order}
        files |= {Path('lean') / n for n in ['lean-toolchain', 'lake-manifest.json', 'lakefile.toml', 'verify.py']}
        files |= {Path('_measurements') / n for n in ['S11_lean_bookkeeping_contract_check.py', 'S11_lean_bookkeeping_source_check.py']}
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
        # Reuse only inspected, hash-matching native evidence; never infer a pass.
        native_path = M / 'S11_lean_bookkeeping_source_checks.json'
        if args.native_evidence:
            if args.native_evidence.resolve() != native_path.resolve():
                raise ValueError('Unexpected native evidence path')
            native = verify.load(native_path)
            if native['status'] != 'PASS' or not all(c['passed'] for c in native['checks']):
                raise ValueError('Native evidence not passing')
            if verify.sha(M / 'S11_lean_bookkeeping_source_check.py') != native['instrument_sha256']:
                raise ValueError('Native instrument drift')
            if any(verify.sha(BASE/p)!=h for p,h in native['source_sha256'].items()):
                raise ValueError('Native source drift')
            result['checks'].append({'name':'native','expected':'PASS','passed':True,
                'mode':'revalidated_recorded_native_evidence','report_sha256':verify.sha(native_path),
                'original_run':'verification_run1; not rerun in this attempt'})
        else:
            execute('native', [sys.executable, str(M / 'S11_lean_bookkeeping_source_check.py')])
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
        preserved()
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
