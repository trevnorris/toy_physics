#!/usr/bin/env python3
"""Portable, sequential replay of completed ledger Lean contracts.

Historical reports supply test *inputs*, never cached pass results or required
local object hashes. Each run gets fresh local objects and immutable input
copies under _scratch. External package caches must be installed by setup.sh.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time

LEAN = Path(__file__).resolve().parent
BASE = LEAN.parent
# Deliberately explicit: an unfinished library is not a completed contract.
CONTRACTS = {
    's9': ('S9_lean', ['S9Pilot']),
    's10': ('S10_lean', ['S10Pilot', 'S10Controls', 'S10Anisotropic', 'S10Audit']),
    'homogeneous': ('S11_lean', ['S11Homogeneous']),
    'd2': ('S11_lean_invariant', ['S11Invariants']),
    'd2-dynamics': ('S11_lean_dynamics', ['S11OddDynamics']),
    'd3': ('S11_lean_d3', ['S11D3Invariants']),
    'd3-bulk': ('S11_lean_d3_bulk', ['S11D3Bulk']),
    'd4': ('S11_lean_d4', ['S11D4Invariants']),
    'd4-odd': ('S11_lean_d4_odd', ['S11D4Odd']),
    'analytic-error': ('S11_lean_analytic_error', ['S11AnalyticError']),
    'variable': ('S11_lean_variable', ['S11VariableCoefficients']),
    'poles': ('S11_lean_pole', ['S11NonlinearPole']),
    'd4-bulk': ('S11_lean_d4_bulk', ['S11D4Bulk']),
    'd5': ('S11_lean_d5', ['S11D5Invariants']),
    'd5-bulk': ('S11_lean_d5_bulk', ['S11D5Bulk']),
}
IMPORT = re.compile(r'^import\s+([\w.]+)\s*$', re.M)
LOCAL = re.compile(r'S(?:9|10|11)\w*(?:\.\w+)*\Z')
ALLOWED_AXIOMS = {'propext', 'Classical.choice', 'Quot.sound'}
BAD_INSTRUMENT = re.compile(
    r'unknown (?:module|namespace|identifier)|unexpected token|maximum (?:recursion|heartbeats)|'
    r'excessive memory|out of memory|PANIC|No such file|failed to create|interrupted|'
    r'failed to synthesize|declaration uses .sorry.|unknown constant', re.I)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def inside(root, relative):
    """Accept normalized in-tree references, never an escaping report path."""
    path = (root / relative).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError(f'Path escapes snapshot: {relative}')
    return path


def module_path(module):
    if not LOCAL.fullmatch(module):
        raise ValueError(f'Not a ledger module: {module}')
    step = 's11' if module.startswith('S11') else 's10' if module.startswith('S10') else 's9'
    return Path('lean') / step / (module.replace('.', '/') + '.lean')


def build_order(roots):
    ordered, visiting, seen = [], set(), set()
    def visit(module):
        if module in seen:
            return
        if module in visiting:
            raise ValueError(f'Import cycle: {module}')
        visiting.add(module)
        for child in IMPORT.findall((BASE / module_path(module)).read_text()):
            if LOCAL.fullmatch(child):
                visit(child)
        visiting.remove(module)
        seen.add(module)
        ordered.append(module)
    for root in roots:
        visit(root)
    return ordered


def control_source(record):
    if 'source' in record:
        source = record['source']
    else:
        source = inside(BASE, record['canonical_source']).read_text()
        patch = record['replacement']
        if source.count(patch['old']) != 1:
            raise ValueError(f"Mutation no longer applies exactly once: {record['name']}")
        source = source.replace(patch['old'], patch['new'], 1)
    if hashlib.sha256(source.encode()).hexdigest() != record['source_sha256']:
        raise ValueError(f"Control source differs from reviewed input: {record['name']}")
    return source


def specifications(names):
    specs = {}
    for name in names:
        prefix, roots = CONTRACTS[name]
        path = BASE / '_measurements' / (prefix + '_contract_checks.json')
        report = load(path)
        if report['status'] != 'PASS':
            raise ValueError(f'Contract does not have passing evidence: {path}')
        controls = []
        for record in report['checks']:
            # Canonical builds are recreated from the current import graph.
            # S10 also has explicit canonical positive rechecks; retain those.
            if 'source' in record or 'replacement' in record:
                source = control_source(record)
            elif name == 's10':
                source = inside(LEAN, record['command'][-1]).read_text()
                if hashlib.sha256(source.encode()).hexdigest() != record['source_sha256']:
                    raise ValueError(f"Canonical positive changed: {record['name']}")
            else:
                continue
            controls.append({
                'name': record['name'], 'source': source,
                'expected': record.get('expected', record.get('expected_outcome')),
                'required': record.get('required_failure_in', []),
                'source_mutation': 'replacement' in record,
            })
        specs[name] = {'prefix': prefix, 'roots': roots, 'controls': controls,
                       'report': str(path.relative_to(BASE)), 'report_sha256': sha(path)}
    return specs


def audit_axioms(output, source):
    expected = len(re.findall(r'^#print axioms\s', source, re.M))
    lists = re.findall(r'depends on axioms:\s*\[([^\]]*)\]', output)
    empty = len(re.findall(r'does not depend on any axioms', output))
    if len(lists) + empty != expected:
        raise ValueError(f'Axiom audit incomplete: expected {expected}, got {len(lists) + empty}')
    for values in lists:
        axioms = {re.sub(r'\s+', '', a) for a in values.split(',')} - {''}
        if not axioms <= ALLOWED_AXIOMS:
            raise ValueError(f'Unapproved axioms: {axioms - ALLOWED_AXIOMS}')
    return expected


def adjudicate(code, output, source, required, source_mutation=False):
    """A failed compiler invocation is not automatically a rejected theorem."""
    if code != 1 or BAD_INSTRUMENT.search(output):
        return False
    required = [required] if isinstance(required, str) else required
    declarations = [(source[:m.start()].count('\n') + 1, m[1]) for m in re.finditer(
        r'^\s*(?:@\[[^\n]*\]\s*)?(?:theorem|def|abbrev)\s+(\w+)', source, re.M)]
    errors = list(re.finditer(r'^.*?:(\d+):\d+: error(?:\([^)]*\))?:', output, re.M))
    intended = []
    for i, error in enumerate(errors):
        declaration = next((n for line, n in reversed(declarations) if line <= int(error[1])), '')
        message = output[error.end():errors[i+1].start() if i+1 < len(errors) else len(output)]
        if declaration in required:
            intended.append(message)
    if required == ['contract_control'] and not source_mutation:
        return (len(errors) == len(intended) == 1 and intended[0].count('⊢ False') == 1
                and not re.search(r'\bwarning(?:\([^)]*\))?:', output))
    # Older source mutants can also break downstream rewrites. Only a named
    # false identity/type obligation counts, never a later rewrite or warning.
    if any(re.search(r'unsolved goals|type mismatch|failed to close|proved that the proposition',
                    message, re.I) for message in intended):
        return True
    # This reviewed S9 source mutant replaces a spatial derivative by time.
    # It reaches a false displayed identity at a rewrite, not `unsolved goals`.
    # Whitelist its exact frozen source and goal; do not accept arbitrary
    # rewrite failures. Its separate unused-variable warning is not evidence.
    known_s9 = '6acfb14822a7f43336acb63c87bd38f9e8aee8beb766f0b74aad9a8f8070e5ed'
    return (source_mutation and hashlib.sha256(source.encode()).hexdigest() == known_s9
            and any('Tactic `rewrite` failed' in m and '⊢' in m
                    and 'waveCovector omega k 0' in m and 'k i)' in m for m in intended))


def run_process(command, cwd, log, timeout, env):
    with log.open('w') as stream:
        process = subprocess.Popen(command, cwd=cwd, env=env, stdout=stream,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            return process.wait(timeout=timeout)
        except BaseException:
            # Includes cancellation and timeout. Children cannot outlive a
            # launcher which exits before them on SIGTERM.
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                pass
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
            raise


def external_environment():
    """Use pinned external caches only; exclude all original local objects."""
    env = json.loads(subprocess.check_output(
        ['lake', 'env', sys.executable, '-c', 'import json,os; print(json.dumps(dict(os.environ)))'],
        cwd=LEAN, text=True, timeout=60))
    packages, paths = [], []
    for package in load(LEAN / 'lake-manifest.json')['packages']:
        path = LEAN / '.lake/packages' / package['name'].strip('«»')
        head = subprocess.check_output(['git', '-C', str(path), 'rev-parse', 'HEAD'], text=True).strip()
        dirty = subprocess.check_output(
            ['git', '-C', str(path), 'status', '--porcelain', '--untracked-files=no'], text=True)
        if head != package['rev'] or dirty:
            raise ValueError(f'Package pin/source mismatch: {package["name"]}')
        packages.append({'name': package['name'], 'commit': head})
        paths.append(str(path / '.lake/build/lib/lean'))
    # Lean itself supplies its standard library path. External native library
    # paths from `lake env` remain intact, but no local ledger olean is used.
    env['LEAN_PATH'] = os.pathsep.join(paths)
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', LEAN_NUM_THREADS='1')
    executable = subprocess.check_output(['elan', 'which', 'lean'], cwd=LEAN, text=True).strip()
    version = subprocess.check_output([executable, '--version'], text=True).strip()
    selected = (LEAN / 'lean-toolchain').read_text().strip().split(':v')[-1]
    if f'version {selected},' not in version and f'version {selected} ' not in version:
        raise ValueError(f'Wrong Lean version: {version}')
    return env, executable, {'packages': packages, 'lean_version': version}


def native_prefixes(names):
    prefixes = [CONTRACTS[n][0] for n in names if n != 's10']
    if 'variable' in names:
        prefixes = ['S11_lean_d3_bulk', 'S11_lean_d4_odd'] + prefixes
    return list(dict.fromkeys(prefixes))


def native_inputs(prefixes):
    """Bounded data/helper copies; none of the production modules are imported."""
    files = set()
    def add(rel):
        path = inside(BASE, rel)
        files.add(path.relative_to(BASE))
    for prefix in prefixes:
        add(f'_measurements/{prefix}_source_check.py')
        report = load(BASE / f'_measurements/{prefix}_source_checks.json')
        sources = report.get('source_sha256', {})
        if isinstance(sources, dict):
            for rel in sources:
                add(rel)
        else:  # Original D2 report predates the path-keyed schema.
            add('scripts/S11_stray_longitudinal_sympy_audit.py')
            add('mathematica/S11_stray_longitudinal_mathematica_audit.wl')
        if prefix == 'S11_lean_pole':
            add('_measurements/S11_lean_pole_preserved_inputs.json')
    return files


def main(argv=None):
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('contracts', nargs='*', help='all or names shown by --list')
    parser.add_argument('--list', action='store_true')
    parser.add_argument('--plan', action='store_true', help='validate inputs and show work; no compiler')
    parser.add_argument('--doctor', action='store_true', help='check installed versions, pins and direct imports')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--build-only', action='store_true', help='proof builds/axiom audits only')
    mode.add_argument('--native-only', action='store_true', help='compact SymPy source checks only')
    parser.add_argument('--timeout', type=int, default=600, help='seconds per child process')
    args = parser.parse_args(argv)
    if args.list:
        print('\n'.join(f'{name}: {", ".join(roots)}' for name, (_, roots) in CONTRACTS.items()))
        return 0
    names = list(CONTRACTS) if args.contracts == ['all'] or args.doctor else args.contracts
    if not names or any(n not in CONTRACTS for n in names) or args.timeout <= 0:
        parser.error('Choose all or contract names from --list; timeout must be positive.')
    names = list(dict.fromkeys(names))
    specs = specifications(names)
    order = build_order([root for s in specs.values() for root in s['roots']])
    prefixes = [] if args.build_only else native_prefixes(names)
    if args.native_only and not prefixes:
        parser.error('S10 has no compact native replay in this runner; select its full proof/control suite.')
    files = {module_path(m) for m in order}
    files |= {Path('lean') / p for p in ['lean-toolchain', 'lakefile.toml', 'lake-manifest.json',
                                       'verify.py', 'requirements.txt']}
    files |= {Path(s['report']) for s in specs.values()}
    files |= native_inputs(prefixes)
    generators = []
    for dim in ['d3', 'd4', 'd5']:
        if any(m.startswith(f'S11{dim.upper()}Invariants') for m in order) and not args.native_only:
            generators.append(Path('_measurements') / f'S11_lean_{dim}_generate.py')
    files |= set(generators)
    # Generators check helper/audit files too, all included in their import graph.
    missing = [str(p) for p in files if not (BASE / p).is_file()]
    if missing:
        raise ValueError(f'Missing input files (restore checkout/annex content): {missing}')
    plan = {'contracts': names, 'local_modules': len(order),
            'control_executions': sum(len(s['controls']) for s in specs.values()),
            'native_checks': prefixes, 'generators': list(map(str, generators)),
            'mode': 'native' if args.native_only else 'build' if args.build_only else 'full'}
    if args.plan:
        print(json.dumps(plan, indent=2))
        return 0
    env, executable, dependencies = (dict(os.environ), None, {}) if args.native_only else external_environment()
    if not args.build_only or generators:
        import sympy
        import mpmath
        if (sympy.__version__, mpmath.__version__) != ('1.14.0', '1.3.0'):
            raise ValueError('Use the pinned .venv from setup.sh (SymPy 1.14.0, mpmath 1.3.0).')
    # Check direct cache imports before creating a long-running proof job.
    direct = {}
    if not args.native_only:
        for module in order:
            for child in IMPORT.findall((BASE / module_path(module)).read_text()):
                if LOCAL.fullmatch(child):
                    continue
                stem = child.replace('.', '/')
                candidates = [Path(p) / (stem + '.olean') for p in env['LEAN_PATH'].split(os.pathsep)]
                obj = next((p for p in candidates if p.is_file()), None)
                if obj is None:
                    raise ValueError(f'Missing external import {child}; run bash setup.sh.')
                direct[child] = {'olean_sha256': sha(obj), 'path': str(obj)}
                source = obj.parents[len(Path(stem).parts) + 3] / (stem + '.lean')
                if source.is_file():
                    direct[child]['source_sha256'] = sha(source)
        dependencies['direct_imports'] = direct
    if args.doctor:
        print(f'Ready: pinned Lean, {len(dependencies["packages"])} packages, '
              f'{len(direct)} direct external imports, Python/SymPy and all completed control inputs.')
        return 0
    runs = BASE / '_scratch/lean_portable'
    runs.mkdir(parents=True, exist_ok=True)
    run = Path(tempfile.mkdtemp(prefix=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ_'), dir=runs))
    snapshot = run / 'snapshot'
    for rel in sorted(files):
        target = snapshot / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(BASE / rel, target)
    before = {str(rel): sha(snapshot / rel) for rel in sorted(files)}
    local_lean = snapshot / 'lean'
    objects = run / 'objects'
    objects.mkdir()
    env['LEAN_PATH'] = str(objects) + (os.pathsep + env['LEAN_PATH'] if 'LEAN_PATH' in env else '')
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    # Existing native/generator instruments use load-bearing Python asserts.
    # Do not let a user's optimization setting silently disable those checks.
    env.pop('PYTHONOPTIMIZE', None)
    logs = run / 'logs'
    logs.mkdir()
    result = {'status': 'RUNNING', 'plan': plan, 'started_utc': datetime.now(timezone.utc).isoformat(),
              'input_sha256': before, 'dependencies': dependencies, 'checks': [],
              'limits': ['Fresh local proof/control execution, not a new independent fidelity review.',
                         'Pinned external caches are reused; no historical local object hashes are required.',
                         'Native checks are compact selected-helper checks, not CAS production reruns.'],
              'resources': {'workers': 1, 'lean_threads': 1, 'allocator_mib': 4096,
                            'process_timeout_seconds': args.timeout}}
    report_path = run / 'report.json'
    def save():
        temporary = report_path.with_suffix('.new')
        temporary.write_text(json.dumps(result, indent=2) + '\n')
        temporary.replace(report_path)
    def execute(name, command, cwd, check):
        log = logs / (name + '.log')
        started = time.monotonic()
        code = run_process(command, cwd, log, args.timeout, env)
        output = log.read_text()
        valid = check(code, output)
        result['checks'].append({'name': name, 'passed': valid, 'exit_status': code,
                                 'seconds': round(time.monotonic() - started, 3),
                                 'command': command, 'log': str(log.relative_to(run)),
                                 'log_sha256': sha(log)})
        save()
        if not valid:
            raise RuntimeError(f'{name} failed; see {log}')
        return output
    def positive(code, output):
        return code == 0 and not re.search(r'\b(?:error|warning)(?:\([^)]*\))?:', output)
    save()
    print(f'Running silently; report: {report_path}', flush=True)
    try:
        for generator in generators:
            execute(generator.stem, [sys.executable, str(snapshot / generator), '--check'], snapshot, positive)
        if not args.native_only:
            lean_command = [executable, '-j1', '-M4096', '-DwarningAsError=true']
            for module in order:
                rel = module_path(module).relative_to('lean')
                source = (local_lean / rel).read_text()
                obj = objects / (module.replace('.', '/') + '.olean')
                obj.parent.mkdir(parents=True, exist_ok=True)
                output = execute('build_' + module, lean_command + ['-o', str(obj), str(rel)], local_lean, positive)
                audit_count = audit_axioms(output, source)
                if not obj.is_file():
                    raise RuntimeError(f'Compiler did not create {obj}')
                result['checks'][-1].update(axiom_audits=audit_count, olean_sha256=sha(obj))
                save()
            if not args.build_only:
                for name, spec in specs.items():
                    for control in spec['controls']:
                        path = local_lean / '_controls' / name / (control['name'] + '.lean')
                        path.parent.mkdir(parents=True, exist_ok=True)
                        path.write_text(control['source'])
                        check = positive if control['expected'] == 'PASS' else (
                            lambda code, out, c=control: adjudicate(code, out, c['source'], c['required'], c['source_mutation']))
                        execute(name + '_' + control['name'], lean_command + [str(path.relative_to(local_lean))], local_lean, check)
        for prefix in prefixes:
            instrument = snapshot / '_measurements' / (prefix + '_source_check.py')
            execute(prefix + '_native', [sys.executable, str(instrument)], snapshot, positive)
            native = load(snapshot / '_measurements' / (prefix + '_source_checks.json'))
            if native['status'] != 'PASS':
                raise RuntimeError(f'Native source check failed: {prefix}')
        # Native output records are expected to change only in the snapshot.
        outputs = {f'_measurements/{p}_source_checks.json' for p in prefixes}
        for rel, digest in before.items():
            if rel not in outputs and sha(snapshot / rel) != digest:
                raise RuntimeError(f'Input changed during execution: {rel}')
            if sha(BASE / rel) != digest:
                raise RuntimeError(f'Live source changed during execution: {rel}; rerun from a stable checkout.')
        for item in direct.values():
            if sha(item['path']) != item['olean_sha256']:
                raise RuntimeError(f'External cache changed during execution: {item["path"]}')
        result['objects'] = {str(p.relative_to(objects)): sha(p) for p in objects.rglob('*.olean')}
        result['status'] = 'PASS'
    except BaseException as error:
        result.update(status='FAILED', error=f'{type(error).__name__}: {error}')
        raise
    finally:
        result['finished_utc'] = datetime.now(timezone.utc).isoformat()
        save()
    print(f'PASS: {report_path}', flush=True)
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as error:
        print(f'Verification stopped: {error}', file=sys.stderr)
        sys.exit(1)
