#!/usr/bin/env python3
"""Blind WL N6 fixed manifest: K_carrier, K_junk, K_split_route.
All labels: MATERIAL_ADVECTED/RHO4_CONSTANT, four valid draws/prime/cell.
The two budget edits are BUDGET_RESTRICTION / NON-KNIFE. Production source
is never executed unrestricted. Exact full emitted assignment text is hashed
in this driver before compact projection; raw output exists only in /tmp.
"""
from __future__ import annotations
import ast
from collections import OrderedDict
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
import S11c_c2_N6_covariance_ablation_harness as h

ROOT = Path(__file__).resolve().parent.parent
ENGINE = ROOT / 'mathematica/S11c_c2_N6_mathematica_audit.wl'
HARNESS = ROOT / 'mathematica/S11c_c2_N6_ablation_harness.wl'
CASE = ['MATERIAL_ADVECTED', 'RHO4_CONSTANT']
ITERATOR = '{case, Tuples[{{"LAB_HELD", "MATERIAL_ADVECTED"}, {"RHO4_CONSTANT", "RHOBR_CONSTANT"}}]}'
DRAWS = 'drawCount = If[0 < bound < 1, Max[8, Ceiling[Log[2^-80/Max[1, unionCount]]/Log[bound]]], 8];'
BUDGET = [('driver_Do', ITERATOR, '{case, {{"MATERIAL_ADVECTED", "RHO4_CONSTANT"}}}'),
          ('probeCase', DRAWS, 'drawCount = 4;')]


def rc(*names):
    return [('WL_S11CC2_N6RC_' + name, None) for name in names]


def cov(*names):
    return [('WL_S11CC2_N6COV_' + name, None) for name in names]

SOURCE = 'sourceChannel = joinFaceMaps[buildContraction[cm[#], es[#] - ms[#], response[#], #, False] &];'
CERTIFIED = rc('R_N6', 'SPLIT_CHECK') + cov('R_COV') + [
    ('WL_S11CC2_N6_' + x + '_GUARD_RESIDUAL', None) for x in ('SLOT', 'CLOSURE')]
KNIVES = [
    h.knife('K_carrier', 'materialNormalKnife', 'materialNormalKnife = 0;',
      'materialNormalKnife = 1;', 'materialNormalKnife = 2;',
      rc('CARRIER_EULERIAN', 'SOURCE_EULERIAN', 'SOURCE_MATERIAL', 'SOURCE_BRIDGE_RESIDUAL') +
      cov('SOURCE_ACTUAL', 'SOURCE_PREDICTED', 'R_COV', 'FROZEN_PHI')),
    h.knife('K_junk', 'actualJunkCoefficient', 'actualJunkCoefficient = 0;',
      'actualJunkCoefficient = 1;', 'actualJunkCoefficient = 2;',
      rc('CARRIER_EULERIAN', 'CARRIER_MATERIAL', 'CARRIER_BRIDGE_RESIDUAL') +
      cov('SOURCE_PREDICTED', 'FROZEN_PHI', 'PHI_DOMAIN_CENSUS')),
    h.knife('K_split_route', 'buildCase', SOURCE, SOURCE.replace(', False]', ', True]'),
      SOURCE.replace('es[#] - ms[#]', '2 (es[#] - ms[#])'),
      rc('CARRIER_EULERIAN', 'CARRIER_MATERIAL', 'CARRIER_BRIDGE_RESIDUAL', 'SOURCE_EULERIAN', 'SOURCE_MATERIAL') + cov('R_COV')),
]
CONFIG = dict(certified=CERTIFIED, knives=KNIVES)


def split_top(text, separator=','):
    """Split held WL text without evaluating any symbol, arithmetic or head.

    Brackets, strings, escaped quotes, associations and nested comments are
    balanced. Returned spans preserve exact emitted/source spelling and order.
    """
    start = i = 0
    stack = []
    quote = False
    comment = 0
    while i < len(text):
        pair = text[i:i+2]
        char = text[i]
        if comment:
            if pair == '(*': comment += 1; i += 2; continue
            if pair == '*)': comment -= 1; i += 2; continue
            i += 1; continue
        if quote:
            if char == '\\': i += 2; continue
            if char == '"': quote = False
            i += 1; continue
        if pair == '(*': comment = 1; i += 2; continue
        if char == '"': quote = True; i += 1; continue
        if pair == '<|': stack.append('|>'); i += 2; continue
        if pair == '|>':
            if not stack or stack.pop() != '|>': raise ValueError('association_balance')
            i += 2; continue
        if char in '[{(': stack.append({'[':']', '{':'}', '(' : ')'}[char])
        elif char in ']})':
            if not stack or stack.pop() != char: raise ValueError('bracket_balance')
        elif not stack and text.startswith(separator, i):
            yield text[start:i]
            i += len(separator); start = i; continue
        i += 1
    if stack or quote or comment: raise ValueError('incomplete_wl_text')
    yield text[start:]


def association(text):
    text = text.strip()
    if not (text.startswith('<|') and text.endswith('|>')):
        raise ValueError('association_schema')
    result = OrderedDict()
    for item in split_top(text[2:-2]):
        if not item.strip(): continue
        parts = list(split_top(item, '->'))
        if len(parts) != 2: raise ValueError('rule_schema')
        key, value = (x.strip() for x in parts)
        if key in result: raise ValueError('duplicate_wl_key')
        result[key] = value
    return result


def literal(text):
    # Only the numeric/string/list subset used by emitted indices and sparse
    # storage is decoded. No eval, symbolic fallback, or circuit reconstruction.
    return ast.literal_eval(text.replace('{', '[').replace('}', ']'))


def sparse_row_counts(text):
    if not text.startswith('SparseArray[') or not text.endswith(']'):
        raise ValueError('sparse_head')
    args = list(split_top(text[len('SparseArray['):-1]))
    if len(args) != 4 or args[0].strip() != 'Automatic':
        raise ValueError('sparse_storage')
    shape, default, storage = literal(args[1]), literal(args[2]), literal(args[3])
    if len(shape) != 2 or default != 0 or storage[0] != 1:
        raise ValueError('sparse_dimensions')
    pointers, column_indices = storage[1]
    values = storage[2]
    if len(pointers) != shape[0] + 1 or pointers[-1] != len(values) or len(column_indices) != len(values):
        raise ValueError('sparse_index')
    counts = [sum(v != 0 for v in values[pointers[i]:pointers[i+1]]) for i in range(shape[0])]
    return counts, shape


def patch_source(source, scope, old, replacement):
    # Top-level semicolons identify complete held assignments/Do definitions.
    statements = list(split_top(source, ';'))
    candidates = []
    offset = 0
    for statement in statements:
        clean = re.sub(r'\(\*.*?\*\)', '', statement, flags=re.S).strip()
        match = (clean.startswith('Do[') and 'buildCase[case]' in clean) if scope == 'driver_Do' else (
            bool(re.match(re.escape(scope) + r'(?:\[|\s*=)', clean)))
        if match: candidates.append((offset, statement))
        offset += len(statement) + 1
    if len(candidates) != 1: raise ValueError('wl_scope_count')
    start, fragment = candidates[0]
    # Include the terminal semicolon, when it belongs to the literal patch.
    if source[start+len(fragment):start+len(fragment)+1] == ';': fragment += ';'
    if fragment.count(old) != 1: raise ValueError('wl_fragment_count')
    patched = source[:start] + fragment.replace(old, replacement, 1) + source[start+len(fragment):]
    list(split_top(patched, ';'))
    return patched


def compact_output(spool, wanted):
    selected = {key[0] for key in wanted}
    raw = {}
    for line in spool.open():
        tag, sep, payload = line.rstrip('\n').partition(' = ')
        if tag in selected or tag in ('WL_S11CC2_N6_LOCAL_PROBE', 'WL_S11CC2_N6_LOCAL_CONSTRUCTION'):
            if tag in raw: raise ValueError('duplicate_wl_tag')
            # Exact assignment bytes excluding only the line delimiter.
            raw[tag] = (h.digest(line.rstrip('\n').encode()), payload)
    metadata_tag = 'WL_S11CC2_N6_LOCAL_PROBE'
    construction_tag = 'WL_S11CC2_N6_LOCAL_CONSTRUCTION'
    if metadata_tag not in raw or construction_tag not in raw:
        return dict(objects={}, metadata=[], missing=wanted)
    construction = association(raw[construction_tag][1])
    catalogs, metadata = {}, []
    for case, case_text in association(raw[metadata_tag][1]).items():
        values = association(case_text)
        sample_index = literal(values['"SAMPLE_INDEX"'])
        groups = OrderedDict()
        for row, (cell, prime, draw) in enumerate(sample_index):
            groups.setdefault((cell, prime), []).append(row)
        catalog = OrderedDict()
        for (cell, prime), rows in groups.items():
            key = '{' + case + ', "PROBE", ' + str(cell) + ', ' + str(prime) + '}'
            valid, rejected, leaves = literal(construction[key])
            if valid != len(rows): raise ValueError('sample_index_count')
            catalog[key] = dict(rows=rows, valid=valid, rejected=rejected, circuit_leaves=leaves)
        catalogs[case] = (sample_index, catalog)
        bounds = []
        for item in split_top(values['"BOUNDS"'][1:-1]):
            entry = association(item)
            bounds.append({k.strip('"'): v for k, v in entry.items()
                           if re.fullmatch(r'-?\d+(?:/[1-9]\d*)?', v)})
        metadata.append(dict(case=literal(case), sha256=raw[metadata_tag][0],
            seeds=literal(values['"SEEDS"']), bounds=bounds,
            family_cardinality=int(values['"FAMILY_CARDINALITY"'])))
    objects = {}
    for tag, _ in wanted:
        if tag not in raw: continue
        sha, payload = raw.pop(tag)
        case_map = association(payload)
        fingerprint = OrderedDict()
        component_keys = []
        scalar_fields = OrderedDict()
        for case, case_text in case_map.items():
            entries = association(case_text)
            if tag.endswith(('FROZEN_PHI', 'PHI_DOMAIN_CENSUS')):
                scalar_fields[case] = {k.strip('"'): v for k, v in entries.items()
                                      if re.fullmatch(r'-?\d+(?:/[1-9]\d*)?', v)}
                continue
            index, catalog = catalogs[case]
            counts = [0] * len(index)
            column_count = 0
            for key, component in entries.items():
                fields = association(component)
                row_counts, shape = sparse_row_counts(fields['"PROBE_NUMERATORS"'])
                if shape[0] != len(index): raise ValueError('probe_row_count')
                component_keys.append([literal(case), key, fields['"COMPONENT_AXES"']])
                counts = [a+b for a,b in zip(counts, row_counts)]
                column_count += shape[1]
            for key, record in catalog.items():
                fingerprint[key] = dict(valid=record['valid'], rejected=record['rejected'],
                    nonzero_count=sum(counts[row] for row in record['rows']),
                    circuit_leaves=record['circuit_leaves'],
                    numerator_count=record['valid'] * column_count)
        objects[h.dumps((tag, None))] = dict(sha256=sha, fingerprint=fingerprint,
            component_keys=component_keys, scalars=scalar_fields)
        del payload, case_map
    raw.clear()
    return dict(objects=objects, metadata=metadata,
                missing=[key for key in wanted if h.dumps(key) not in objects])


def project(transcript, value):
    if value is None: return {'MISSING': 1}
    result = dict(value)
    keys = result.pop('component_keys')
    if keys:
        ident = h.digest(h.dumps(keys).encode())
        if ident not in transcript.schemas:
            transcript.schemas[ident] = keys
            transcript.emit('COMPONENT_KEYS', dict(id=ident, keys=keys))
        result['component_keys'] = {'ref': ident}
    fingerprint = result.pop('fingerprint')
    if fingerprint:
        ident = h.digest(h.dumps(fingerprint).encode())
        if ident not in transcript.schemas:
            transcript.schemas[ident] = fingerprint
            transcript.emit('PIT_FINGERPRINT', dict(id=ident, fingerprint=fingerprint))
        result['fingerprint'] = {'ref': ident}
    return result


def delta(b, c):
    if b is None or c is None:
        return h.compact_delta(b, c)
    result = dict(digest_equal=int(b['sha256'] == c['sha256']), fingerprint_delta={}, MISSING=[])
    for key, values in b['fingerprint'].items():
        if key not in c['fingerprint']:
            result['MISSING'].append(dict(key=key, side='corrupted')); continue
        moved = {field: c['fingerprint'][key][field] - value for field, value in values.items()
                 if c['fingerprint'][key][field] != value}
        if moved: result['fingerprint_delta'][key] = moved
    result['MISSING'].extend(dict(key=key, side='baseline') for key in c['fingerprint'] if key not in b['fingerprint'])
    bm, cm = {h.dumps(k) for k in b['component_keys']}, {h.dumps(k) for k in c['component_keys']}
    result['MISSING'].extend(dict(component=json.loads(k), side='corrupted') for k in sorted(bm-cm))
    result['MISSING'].extend(dict(component=json.loads(k), side='baseline') for k in sorted(cm-bm))
    from fractions import Fraction
    result['scalar_delta'] = {}
    for case, fields in b['scalars'].items():
        if case not in c['scalars']:
            result['MISSING'].append(dict(key=case, side='corrupted')); continue
        result['scalar_delta'][case] = {k: str(Fraction(c['scalars'][case][k])-Fraction(v))
                                       for k,v in fields.items() if k in c['scalars'][case]}
        result['MISSING'].extend(dict(key=[case,k], side='corrupted') for k in fields if k not in c['scalars'][case])
        result['MISSING'].extend(dict(key=[case,k], side='baseline') for k in c['scalars'][case] if k not in fields)
    result['MISSING'].extend(dict(key=case, side='baseline') for case in c['scalars'] if case not in b['scalars'])
    return result


def main():
    script = Path(__file__).resolve()
    transcript = h.Transcript(ROOT / '_measurements/S11c_c2_N6_ablation_harness_wl.md')
    original = ENGINE.read_text()
    prod_sha = h.digest(ENGINE.read_bytes())
    budget_source = original
    budget_patches = []
    for scope, old, new in BUDGET:
        next_source = patch_source(budget_source, scope, old, new)
        diff = h.source_diff(budget_source, next_source, ENGINE.name)
        budget_patches.append(dict(classification='BUDGET_RESTRICTION / NON-KNIFE', scope=scope,
            old=old, replacement=new, exact_fragment_count=1, diff=diff, diff_sha256=h.digest(diff.encode()),
            before_sha256=h.digest(budget_source.encode()), after_sha256=h.digest(next_source.encode())))
        budget_source = next_source
    for _, k, replacement, _ in h.variants(CONFIG):
        if k: patch_source(budget_source, k['scope'], k['old'], replacement)
    transcript.emit('MANIFEST', dict(engine=str(ENGINE), case=CASE, draws_per_prime_cell=4,
        invocation=shlex.join([sys.executable, str(script)]), cwd=os.getcwd(),
        production_sha256=prod_sha, harness_sha256=h.digest(HARNESS.read_bytes()), driver_sha256=h.digest(script.read_bytes()),
        python=sys.version.split()[0], wolfram_version='each label: HELD_SOURCE_GUARDS.wolfram_version',
        knives=KNIVES, budget_restrictions=budget_patches, printed_objects=h.selected(CONFIG),
        digest_input='exact full emitted WL_S11CC2_* = payload text, without newline, before compaction',
        delta_schema='changed shared fields only; unlisted shared deltas are 0; absence is MISSING'))
    baseline, status = None, 0
    for label, k, replacement, classification in h.variants(CONFIG):
        with tempfile.TemporaryDirectory(prefix='s11cc2_n6_wl_') as directory:
            tree = Path(directory)
            source = budget_source if k is None else patch_source(budget_source, k['scope'], k['old'], replacement)
            prod_path, budget_path, exec_path = (tree / name for name in ('production.wl', 'budget.wl', ENGINE.name))
            prod_path.write_text(original); budget_path.write_text(budget_source); exec_path.write_text(source)
            guard = tree / 'held_guards.json'
            plan = dict(production=str(prod_path), budget=str(budget_path), executed=str(exec_path), guard=str(guard),
                        knife_count=int(k is not None), scope=k['scope'] if k else 'NONE',
                        old=k['old'] if k else 'Null', replacement=replacement if k else 'Null')
            plan_path = tree / 'plan.json'; plan_path.write_text(h.dumps(plan))
            command = ['timeout', '--kill-after=5', '600', shutil.which('wolframscript') or 'wolframscript', '-file', str(HARNESS)]
            knife_diff = h.source_diff(budget_source, source, ENGINE.name)
            transcript.emit('RUN_MANIFEST', dict(label=label, classification=classification, scope=plan['scope'],
                old=plan['old'], replacement=plan['replacement'], production_sha256=prod_sha,
                budget_sha256=h.digest(budget_source.encode()), executed_sha256=h.digest(source.encode()),
                source_path=str(exec_path), production_path=str(ENGINE), budget_patches=budget_patches,
                knife_diff=knife_diff, knife_diff_sha256=h.digest(knife_diff.encode()),
                invocation='S11CC2_HARNESS_PLAN=' + shlex.quote(str(plan_path)) + ' ' + shlex.join(command)))
            with (tree / 'engine.stdout').open('wb') as out, (tree / 'worker.stderr').open('wb') as err:
                proc = subprocess.run(command, stdout=out, stderr=err, env=dict(os.environ, S11CC2_HARNESS_PLAN=str(plan_path)))
            parse_error = None
            try:
                result = compact_output(tree / 'engine.stdout', h.selected(CONFIG))
            except Exception as exc:
                parse_error = type(exc).__name__
                result = dict(objects={}, metadata=[], missing=h.selected(CONFIG))
            if baseline is None: baseline = result['objects']
            objects = result['objects']
            sets = [('CERTIFIED', CERTIFIED), ('DEAD', k['dead'] if k else list(dict.fromkeys(key for kn in KNIVES for key in kn['dead'])))]
            for role, keys in sets:
                for key in keys:
                    b, c = baseline.get(h.dumps(key)), objects.get(h.dumps(key))
                    transcript.emit(role, dict(label=label, object=key[0], baseline=project(transcript,b),
                                              corrupted=project(transcript,c), diff=delta(b,c)))
            transcript.emit('ENGINE_METADATA', dict(label=label, values=result['metadata']))
            transcript.emit('HELD_SOURCE_GUARDS', dict(label=label, values=json.loads(guard.read_text()) if guard.exists() else {'MISSING': 1}))
            transcript.emit('INTEGRITY_GUARDS', dict(label=label, missing=result['missing'], parser_error_type=parse_error,
                production_immutable=int(h.digest(ENGINE.read_bytes()) == prod_sha), **h.bounded_diagnostics(tree, proc.returncode)))
            if proc.returncode or result['missing'] or not guard.exists() or h.digest(ENGINE.read_bytes()) != prod_sha:
                status = 2; break
    transcript.finish(status)
    return status

if __name__ == '__main__':
    sys.exit(main())
