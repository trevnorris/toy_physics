#!/usr/bin/env python3
"""Fixed N6 covariance knives and shared, payload-only harness infrastructure.

Manifest: K_junk (module ACTUAL_JUNK), K_circular (run amplitude block),
K_rank (prolonged_phi.image_of images assignment). Each has the literal FORM,
identity and same-site coefficient x2 from the fixed four-engine directive.
Only production run/emit/PIT payloads are observed. No payload classification.
"""
from __future__ import annotations
import argparse
import ast
import contextlib
import difflib
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent.parent
PREFIX = 'S11CC2_'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def dumps(value):
    return json.dumps(value, ensure_ascii=True, separators=(',', ':'), sort_keys=False)


def knife(name, scope, old, form, twice, dead):
    return dict(name=name, scope=scope, old=old, form=form, twice=twice, dead=dead)


CIRCULAR = '''    mu_pred = predicted_amplitude(mu_e, phi)
    mu_actual, mu_baseline = actual_amplitudes(a, b, inputs, alpha, rho,
                                               ACTUAL_A_RHO, ACTUAL_JUNK)'''
RANK = '''            images[atom] = b.total_derivative(image_of(parent), direction,
                                             background_depth=3)'''
CERTIFIED = [('S11CC2_N6COV_' + s, None) for s in (
    'R_COV', 'SOURCE_ACTUAL', 'SOURCE_PREDICTED', 'R_COV_INCREMENT',
    'R_COV_CONTROL_DELTA', 'FROZEN_PHI', 'PHI_DOMAIN_CENSUS')]

def covkeys(*names):
    return [('S11CC2_N6COV_' + s, None) for s in names]

KNIVES = [
    knife('K_junk', 'ACTUAL_JUNK', 'ACTUAL_JUNK = sp.Integer(0)',
          'ACTUAL_JUNK = sp.Integer(1)', 'ACTUAL_JUNK = sp.Integer(2)',
          covkeys('SOURCE_PREDICTED', 'FROZEN_PHI', 'PHI_DOMAIN_CENSUS')),
    knife('K_circular', 'run', CIRCULAR,
          '''    mu_actual, mu_baseline = actual_amplitudes(a, b, inputs, alpha, rho,
                                               ACTUAL_A_RHO, ACTUAL_JUNK)
    mu_pred = mu_actual''',
          CIRCULAR.replace('mu_pred = predicted_amplitude(mu_e, phi)',
                           'mu_pred = [2 * term for term in predicted_amplitude(mu_e, phi)]'),
          covkeys('SOURCE_ACTUAL', 'FROZEN_PHI', 'PHI_DOMAIN_CENSUS')),
    knife('K_rank', 'prolonged_phi', RANK,
          '''            images[atom] = (image_of(parent) if len(paths[atom][1]) >= 2 else
                            b.total_derivative(image_of(parent), direction,
                                               background_depth=3))''',
          '''            images[atom] = ((2 if len(paths[atom][1]) >= 2 else 1) *
                            b.total_derivative(image_of(parent), direction,
                                               background_depth=3))''', covkeys('SOURCE_ACTUAL')),
]
CONFIG = dict(engine='S11c_c2_N6_covariance_sympy.py', seed=110603,
              certified=CERTIFIED, knives=KNIVES, extra=[])


def patch_source(source, scope, old, new):
    tree = ast.parse(source)
    nodes = [node for node in tree.body if
             (isinstance(node, ast.FunctionDef) and node.name == scope) or
             (isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == scope for t in node.targets))]
    if len(nodes) != 1:
        raise ValueError('scope_count')
    node = nodes[0]
    lines = source.splitlines(keepends=True)
    before = ''.join(lines[node.lineno - 1:node.end_lineno])
    # K_operand_swap is restricted further to the single ds assignment in run.
    if old == 'ms[s].get(w, sp.S.Zero)':
        assignments = [x for x in ast.walk(node) if isinstance(x, ast.Assign)
                       and any(isinstance(t, ast.Name) and t.id == 'ds' for t in x.targets)]
        if len(assignments) != 1:
            raise ValueError('ds_count')
        node = assignments[0]
        before = ''.join(lines[node.lineno - 1:node.end_lineno])
    if before.count(old) != 1:
        raise ValueError('fragment_count')
    after = before.replace(old, new, 1)
    result = ''.join(lines[:node.lineno - 1]) + after + ''.join(lines[node.end_lineno:])
    ast.parse(result)
    return result


def variants(config):
    yield 'CANONICAL', None, None, 'CANONICAL'
    yield 'UNABLATED_COPY', None, None, 'UNABLATED_COPY'
    for k in config['knives']:
        for suffix, field, cls in [('FORM', 'form', 'FORM'),
                                   ('IDENTITY', 'old', 'NO-OP / IDENTITY'),
                                   ('X2', 'twice', 'COEFFICIENT x2')]:
            yield k['name'] + '/' + suffix, k, k[field], cls
    for k in config.get('extra', []):
        yield k['name'], k, k['form'], 'COEFFICIENT'


def selected(config):
    return list(dict.fromkeys(tuple(k) for k in config['certified'] +
                             [key for knife_ in config['knives'] for key in knife_['dead']]))


def scalars(data):
    # No recursive table traversal and no strings other than literal rationals.
    if not isinstance(data, dict):
        return {'value': data} if type(data) is int else {}
    return {k: v for k, v in data.items() if type(v) is int or
            (isinstance(v, str) and re.fullmatch(r'-?\d+/[1-9]\d*', v))}


def worker(config, path, tree, destination):
    """Fresh worker: ordinary sibling-first imports, engine-owned emit bindings."""
    sys.path.insert(0, str(tree / 'scripts'))
    os.environ['S11CB_PROJECTION_WORKERS'] = '1'
    spool = tree / 'engine.stdout'
    with spool.open('w') as stream, contextlib.redirect_stdout(stream):
        spec = importlib.util.spec_from_file_location(path.stem, path)
        engine = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = engine
        spec.loader.exec_module(engine)
        args = argparse.Namespace(anchoring='LAB_HELD', density='RHOBR_CONSTANT',
                                  draws=8, seed=config['seed'])
        if 'covariance' in path.name:
            n = engine.n
            old_n, old_r = n.emit, engine.r.emit
            n.emit = engine.r.emit = engine.emit
            try:
                engine.run(args)
            finally:
                n.emit, engine.r.emit = old_n, old_r
        elif 'reconcile' in path.name:
            n = engine.n
            old_n = n.emit
            n.emit = engine.emit
            try:
                engine.run(args)
            finally:
                n.emit = old_n
        else:
            n = engine
            engine.run(args)
    objects, metadata = {}, []
    wanted = selected(config)
    for line in spool.open():
        if not line.startswith('{'):
            continue
        item = json.loads(line)
        key = (item.get('object'), item.get('probe', item.get('route')))
        data = item.get('data')
        if key in wanted:
            if isinstance(data, dict) and 'numerator_denominator' in data:
                if dumps(key) in objects:
                    raise ValueError('duplicate_pit_tag')
                # Hash the complete PARSED emitted PIT data, before projection.
                full_sha = n.sha(data)
                columns = data['columns']
                bitmap = data['nonzero_modular_numerator']
                if len(columns) != len(bitmap) or len({dumps(k) for k in columns}) != len(columns):
                    raise ValueError('pit_schema')
                objects[dumps(key)] = dict(sha256=full_sha, columns=columns,
                    nonzero_modular_numerator=''.join('1' if b else '0' for b in bitmap),
                    nonzero_count=sum(bitmap), column_count=len(columns))
            elif key[0].endswith(('FROZEN_PHI', 'PHI_DOMAIN_CENSUS')):
                objects[dumps(key)] = dict(sha256=n.sha(data), scalars=scalars(data))
        elif key[0] and key[0].endswith(('PIT_PROVENANCE', 'PIT_REJECTIONS', 'PRIMES')):
            meta = dict(object=key[0], sha256=n.sha(data), scalars=scalars(data))
            if key[0].endswith('PRIMES'):
                meta['primes'] = data['values']
            metadata.append(meta)
        del item, data
    missing = [key for key in wanted if dumps(key) not in objects]
    imports = {name: dict(path=str(Path(module.__file__).resolve()),
                         sha256=digest(Path(module.__file__).read_bytes()))
               for name, module in list(sys.modules.items())
               if getattr(module, '__file__', None) and
               (name.startswith('S11') or name == 'ledger_fold')}
    result = dict(objects=objects, metadata=metadata, missing=missing, imports=imports,
                  python=sys.version.split()[0], sympy=engine.sp.__version__)
    destination.write_text(dumps(result))
    spool.unlink()


class Transcript:
    def __init__(self, path):
        self.path, self.lines, self.schemas = path, [], {}

    def emit(self, tag, value):
        line = tag + ' ' + dumps(value) + '\n'
        self.lines.append(line)
        print(line, end='', flush=True)

    def project(self, record):
        if record is None:
            return {'MISSING': 1}
        result = dict(record)
        columns = result.pop('columns', None)
        if columns is not None:
            ident = digest(dumps(columns).encode())
            if ident not in self.schemas:
                self.schemas[ident] = columns
                self.emit('COLUMN_KEYS', dict(id=ident, columns=columns))
            result['columns'] = {'ref': ident}
        return result

    def triple(self, label, key, baseline, corrupted, role):
        b, c = baseline.get(dumps(key)), corrupted.get(dumps(key))
        pb, pc = self.project(b), self.project(c)
        delta = compact_delta(b, c)
        self.emit(role, dict(label=label, object=key[0], probe=key[1],
                             baseline=pb, corrupted=pc, diff=delta))

    def finish(self, status):
        body = ''.join(self.lines).encode()
        self.emit('TRANSCRIPT_GUARD', dict(bytes_before_guard=len(body),
                  sha256_before_guard=digest(body), exit_code=status))
        self.path.write_text('```text\n' + ''.join(self.lines) + '```\n')


def compact_delta(b, c):
    if b is None or c is None:
        return dict(baseline='MISSING' if b is None else 'PRESENT',
                    corrupted='MISSING' if c is None else 'PRESENT')
    result = {'digest_equal': int(b['sha256'] == c['sha256'])}
    if 'columns' in b and 'columns' in c:
        bm = dict(zip(map(dumps, b['columns']), map(int, b['nonzero_modular_numerator'])))
        cm = dict(zip(map(dumps, c['columns']), map(int, c['nonzero_modular_numerator'])))
        shared = [k for k in bm if k in cm]
        result.update(nonzero_count=c['nonzero_count'] - b['nonzero_count'],
                      shared_key_count=len(shared),
                      bitmap_delta=[[json.loads(k), cm[k] - bm[k]] for k in shared if cm[k] != bm[k]],
                      MISSING=[dict(key=json.loads(k), side='corrupted') for k in bm if k not in cm] +
                              [dict(key=json.loads(k), side='baseline') for k in cm if k not in bm])
    else:
        bm, cm = b.get('scalars', {}), c.get('scalars', {})
        from fractions import Fraction
        result['scalar_delta'] = {k: str(Fraction(cm[k]) - Fraction(bm[k])) for k in bm if k in cm}
        result['MISSING'] = [dict(key=k, side='corrupted') for k in bm if k not in cm] + [dict(key=k, side='baseline') for k in cm if k not in bm]
    return result


def copy_tree(tree):
    target = tree / 'scripts'
    target.mkdir()
    # All local Python siblings, including the ledger exports, resolve in this tree.
    for path in (ROOT / 'scripts').glob('*.py'):
        shutil.copy2(path, target / path.name)
    (tree / '_measurements').mkdir()
    spec = ROOT / '_measurements/S11c_c2_N6_route2_spec_astra.md'
    shutil.copy2(spec, tree / '_measurements' / spec.name)


def source_diff(before, after, name):
    return ''.join(difflib.unified_diff(before.splitlines(True), after.splitlines(True),
                                       fromfile=name + ':before', tofile=name + ':after'))


def bounded_diagnostics(tree, exit_code):
    result = {'exit_code': exit_code}
    for name in ('engine.stdout', 'worker.stdout', 'worker.stderr'):
        p = tree / name
        if p.exists():
            data = p.read_bytes()
            result[name] = dict(bytes=len(data), sha256=digest(data))
            # A stderr tail is hashed, never echoed: exceptions can embed objects.
            if name.endswith('stderr'):
                result[name]['tail_bytes'] = min(len(data), 2048)
                result[name]['tail_sha256'] = digest(data[-2048:])
    return result


def run_harness(config, script):
    engine = ROOT / 'scripts' / config['engine']
    original = engine.read_text()
    source_sha = digest(engine.read_bytes())
    transcript = Transcript(ROOT / '_measurements' / (script.stem + '.md'))
    # Validate all exact sites and complete ASTs before starting any engine.
    for _, k, replacement, _ in variants(config):
        if k:
            patch_source(original, k['scope'], k['old'], replacement)
    import sympy as sp
    transcript.emit('MANIFEST', dict(engine=str(engine), case=['LAB_HELD', 'RHOBR_CONSTANT'],
        invocation=shlex.join([sys.executable, str(script)]), cwd=os.getcwd(),
        pinned_args=['--anchoring', 'LAB_HELD', '--density', 'RHOBR_CONSTANT', '--draws', '8', '--seed', str(config['seed'])],
        production_sha256=source_sha, harness_sha256=digest(script.read_bytes()),
        infrastructure_sha256=digest(Path(__file__).read_bytes()), python=sys.version.split()[0], sympy=sp.__version__,
        knives=config['knives'], additional_coefficients=config.get('extra', []),
        printed_objects=selected(config), digest_input='n.sha(parsed emitted PIT data), before projection',
        bitmap_delta='changed shared keys only; unlisted shared deltas are 0; absence is MISSING'))
    baseline, status = None, 0
    for label, k, replacement, classification in variants(config):
        with tempfile.TemporaryDirectory(prefix='s11cc2_n6_') as directory:
            tree = Path(directory)
            copy_tree(tree)
            source = original if k is None else patch_source(original, k['scope'], k['old'], replacement)
            temporary_engine = tree / 'scripts' / engine.name
            temporary_engine.write_text(source)
            path = engine if label == 'CANONICAL' else temporary_engine
            dest = tree / 'compact.json'
            command = [sys.executable, '-B', str(script), '--worker', str(path), str(tree), str(dest)]
            patch = source_diff(original, source, engine.name)
            transcript.emit('RUN_MANIFEST', dict(label=label, classification=classification,
                scope=k['scope'] if k else None, old=k['old'] if k else None,
                replacement=replacement, exact_fragment_count=1 if k else 0,
                production_sha256=source_sha, executed_sha256=digest(path.read_bytes()),
                source_path=str(path), sibling_path=str(tree / 'scripts'), invocation=shlex.join(command),
                knife_diff=patch, knife_diff_sha256=digest(patch.encode())))
            with (tree / 'worker.stdout').open('wb') as out, (tree / 'worker.stderr').open('wb') as err:
                proc = subprocess.run(command, stdout=out, stderr=err,
                    env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONHASHSEED='0', S11CB_PROJECTION_WORKERS='1'))
            result = json.loads(dest.read_text()) if dest.exists() else dict(objects={}, missing=selected(config))
            objects = result['objects']
            if baseline is None:
                baseline = objects
            for key in config['certified']:
                role = 'NATIVE_CONTROL' if 'CONTROL_INDEPENDENCE_' in key[0] else 'CERTIFIED'
                transcript.triple(label, key, baseline, objects, role)
            for key in (k['dead'] if k else list(dict.fromkeys(key for kn in config['knives'] for key in kn['dead']))):
                transcript.triple(label, key, baseline, objects, 'DEAD')
            transcript.emit('ENGINE_METADATA', dict(label=label, values=result.get('metadata', [])))
            transcript.emit('SOURCE_PROVENANCE', dict(label=label, imports=result.get('imports', {})))
            transcript.emit('INTEGRITY_GUARDS', dict(label=label, parseable=1,
                missing=result['missing'], serialization=int(dest.exists()),
                production_immutable=int(digest(engine.read_bytes()) == source_sha),
                **bounded_diagnostics(tree, proc.returncode)))
            if proc.returncode or result['missing'] or digest(engine.read_bytes()) != source_sha:
                status = 2
                break
    transcript.finish(status)
    return status


def main(config=CONFIG, script=Path(__file__).resolve()):
    if len(sys.argv) > 1 and sys.argv[1] == '--worker':
        try:
            worker(config, *(Path(p) for p in sys.argv[2:]))
        except BaseException as exc:
            # Bounded, non-object-bearing diagnostics even for engine exceptions.
            print(dumps({'worker_exception_type': type(exc).__name__}), file=sys.stderr)
            return 2
        return 0
    return run_harness(config, script)

if __name__ == '__main__':
    sys.exit(main())
