#!/usr/bin/env python3
"""Mechanical retrieval for the O2 record; never imports an engine or CAS.

Only reads filed text/JSON, retrieves stored fields, counts literal matches,
and reports the shape of named stored collections. JSON formatting changes
whitespace only. No differences, residuals, predicates of physics, or verdicts
are computed. All comparison values below were printed by the accepted
comparator. --write writes only the record's measurements file.
"""
from __future__ import annotations

import argparse
import contextlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
V3 = ROOT / 'research/pde_ledger_v3'
SELF = 'research/pde_ledger_v3/scripts/O2_record_measurements.py'
OUTPUT = V3 / 'steps/_measurements/O2_record_measurements.md'
STREAM = V3 / 'scripts/out/O2_cross_engine_comparator.out'
REGISTER = 'research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md'
BASELINE = 'c9db665d'

SOURCES = [
    'directives/O2_steady_brane_balance_scoping.md',
    'directives/O2_premise_decision_list.md',
    'directives/O2_input_contract.md',
    'directives/O2_SHARED_PHYSICS.md',
    'directives/_measurements/O2_production_runs.md',
    'directives/_measurements/O2_build_r3_review_disposition.md',
    'directives/_measurements/O2_comparator_production_run.md',
    'directives/_measurements/O2_comparator_build_r5_review_disposition.md',
    'directives/_measurements/O2_record_directive_review_disposition.md',
]


def dump(obj):
    print(json.dumps(obj, ensure_ascii=False, separators=(',', ':')))


def rows():
    with STREAM.open() as source:
        for number, line in enumerate(source, 1):
            yield number, json.loads(line), line


def fields(obj, keys):
    return {key: obj[key] for key in keys if key in obj}


def comparison_nodes(obj, path=()):
    """Walk named stored result collections, without evaluating operands."""
    yield path, obj
    for key, child in obj.get('children', {}).items():
        yield from comparison_nodes(child, (*path, key))
    for index, child in enumerate(obj.get('operands', [])):
        yield from comparison_nodes(child, (*path, 'operands', index))


def lookup(name):
    if name == 'existence':
        for rel in (
            'scripts/out/O2_live_balance_sympy_audit.out',
            'mathematica/out/O2_live_balance_mathematica_audit.out',
            'scripts/out/O2_cross_engine_comparator.out',
            'scripts/out/O2_cross_engine_comparator_accounting.jsonl',
            'scripts/out/O2_cross_engine_comparator_catalog.json',
        ):
            dump({'path': rel, 'exists': (V3 / rel).exists()})
    elif name == 'sources':
        for rel in SOURCES:
            print('\nSOURCE ' + rel)
            for number, line in enumerate((V3 / rel).read_text().splitlines(), 1):
                if rel.endswith('O2_build_r3_review_disposition.md') and not (
                        13 <= number <= 26 or 45 <= number <= 57):
                    continue
                if rel.endswith('O2_comparator_build_r5_review_disposition.md') and not (
                        number <= 12 or 62 <= number <= 92):
                    continue
                print(f'{number}: {line}')
        rel = 'scripts/O2_cross_engine_comparator.py'
        print('\nSOURCE ' + rel + ' (role-binding label provenance)')
        for number, line in enumerate((V3 / rel).read_text().splitlines(), 1):
            if 1100 <= number <= 1129:
                print(f'{number}: {line}')
    elif name == 'catalog':
        # Complete stored catalog, not a newly inferred join map.
        print((V3 / 'scripts/out/O2_cross_engine_comparator_catalog.json').read_text(), end='')
    elif name == 'accounting':
        with (V3 / 'scripts/out/O2_cross_engine_comparator_accounting.jsonl').open() as source:
            for line in source:
                if json.loads(line).get('kind') == 'stream_object_count':
                    print(line, end='')
    elif name == 'scope':
        for number, obj, line in rows():
            if obj['kind'] == 'comparison_scope':
                print(f'STREAM_LINE {number}')
                print(line, end='')
    elif name == 'index':
        for number, obj, _ in rows():
            if obj['kind'] in ('accounting', 'balance_operands'):
                dump({'stream_line': number, **fields(obj, (
                    'kind', 'row', 'py_path', 'wl_path', 'accounting',
                    'parsed_leaves', 'accounted_leaves', 'compared_leaves',
                    'unaccounted_leaves'))})
    elif name == 'differences':
        # Includes empty arrays. Nothing is selected on a computed equality,
        # nonemptiness, zero, tolerance, or acceptance predicate.
        for number, obj, _ in rows():
            if obj['kind'] == 'structure':
                dump({'stream_line': number, **fields(obj, (
                    'kind', 'row', 'serialization_differences', 'serialization_scope'))})
                for action in obj['action_comparison']:
                    dump({'stream_line': number, 'row': obj['row'], **fields(action, (
                        'role', 'unpaired_reason', 'differences'))})
            elif obj['kind'] == 'residual':
                for path, result in comparison_nodes(obj['comparison']):
                    if 'structural_delta' in result:
                        dump({'stream_line': number, 'row': obj['row'],
                              'comparison_path': path, **fields(result, (
                                  'structural_delta', 'structural_delta_scope'))})
            elif obj['kind'] == 'balance_comparison':
                dump({'stream_line': number, **fields(obj, (
                    'kind', 'row', 'component')),
                    'entry_differences': obj['comparison']['entry_differences']})
    elif name == 'inventory-limits':
        for number, obj, _ in rows():
            if obj['kind'] == 'balance_entries':
                print(f"STREAM_LINE {number} ROW {obj['row']} COMPONENT {obj['component']}")
                for engine, entries in obj['entries'].items():
                    dump({'stored_object': 'entries/' + engine, 'length': len(entries)})
                    for entry in entries:
                        dump({'engine': engine, **fields(entry, (
                            'role', 'orientation', 'open_free'))})
    elif name == 'closed-results':
        # All stored outcome/reason/residual leaves, not only successful rows.
        for number, obj, _ in rows():
            if obj['kind'] in ('residual', 'balance_comparison'):
                print(f"STREAM_LINE {number} ROW {obj['row']} COMPONENT {obj.get('component', [])}")
                if obj['kind'] == 'balance_comparison':
                    dump({'comparison_path': ['closed'],
                          **obj['comparison']['closed_residual']})
                else:
                    for path, result in comparison_nodes(obj['comparison']):
                        dump({'comparison_path': path, **fields(result, (
                            'outcome', 'reason', 'residual', 'head',
                            'py_present', 'wl_present',
                            'point_policy', 'points')),
                            **({'children_keys': list(result['children'])}
                               if 'children' in result else {})})
    elif name == 'unjoined':
        for number, obj, _ in rows():
            if obj['kind'] == 'unjoined':
                dump({'stream_line': number, **fields(obj, (
                    'kind', 'accounting', 'engine', 'path', 'reason',
                    'parsed_leaves', 'compared_leaves'))})
    elif name == 'shape':
        catalog = json.loads((V3 / 'scripts/out/O2_cross_engine_comparator_catalog.json').read_text())
        for key in ('joins', 'actions', 'balances', 'names'):
            dump({'stored_object': 'catalog/' + key, 'length': len(catalog[key])})
        kinds = ('comparison_scope', 'operands', 'structure', 'residual', 'accounting',
                     'balance_operands', 'balance_entries', 'balance_comparison',
                     'unjoined', 'stream_object_count')
        reasons = (
            'held, OPEN or unsupported applied head',
            'named operand has no supplied scalar value',
            'container structure differs', 'nested sibling absent',
            'text, native boolean or name is not a subtractable operand',
            'bindings alone are not value evidence')
        residual_tokens = ('"outcome":"exact"', '"outcome":"not_formed"',
                           '"residual":"0"', '"residual":',
                           *(json.dumps('reason') + ':' + json.dumps(reason)
                             for reason in reasons))
        counts = {kind: 0 for kind in kinds}
        residual_counts = {token: 0 for token in residual_tokens}
        paired = empty_paired = 0
        empty_delta = '{"head":[],"live_arguments":[],"named_OPEN_operands":[],"orientation":[],"role":[]}'
        for _, obj, line in rows():
            for kind in kinds:
                counts[kind] += line.count('"kind":"' + kind + '"')
            if obj['kind'] == 'residual':
                for token in residual_tokens:
                    residual_counts[token] += line.count(token)
            elif obj['kind'] == 'structure':
                for action in obj['action_comparison']:
                    stored = json.dumps(action, ensure_ascii=False, separators=(',', ':'))
                    # Literal null is the comparator's filed paired-role marker.
                    if '"unpaired_reason":null' in stored:
                        paired += stored.count('"unpaired_reason":null')
                        empty_paired += stored.count('"differences":' + empty_delta)
        for kind, count in counts.items():
            dump({'stored_object': 'stream', 'literal': '"kind":"' + kind + '"', 'count': count})
        for token, count in residual_counts.items():
            dump({'stored_object': 'stream/residual rows', 'literal': token, 'count': count})
        dump({'stored_object': 'structure/action_comparison',
              'literal': '"unpaired_reason":null', 'count': paired})
        dump({'stored_object': 'paired structure/action_comparison',
              'literal': '"differences":' + empty_delta, 'count': empty_paired})
    elif name == 'register':
        before = subprocess.run(
            ['git', 'show', BASELINE + ':' + REGISTER], cwd=ROOT,
            check=True, text=True, capture_output=True).stdout
        after = (ROOT / REGISTER).read_text()
        for label, content in ((BASELINE, before), ('working-tree', after)):
            print('REGISTER ' + label)
            dump({'literal': '### R-', 'count': content.count('### R-')})
            selected = False
            for line in content.splitlines():
                if line.startswith('### '):
                    selected = any(line.startswith('### ' + entry + ' —') for entry in (
                        'R-S1-02', 'R-S1-03', 'R-S8-04', 'R-S8-05', 'R-S8-06',
                        'R-S12-01', 'R-S12-02', 'R-O2-01'))
                elif line.startswith('## '):
                    selected = False
                if selected:
                    print(line)
            for heading, ending in (
                    ('## Entry schema', '\n---'),
                    ('### O2 pass', '\n### Method for population passes'),
                    ('**Directive/register disagreement', '\n#### Inputs'),
                    ('### Method for population passes', '\n## When this is done')):
                start = content.find(heading)
                if start >= 0:
                    end = content.find(ending, start + len(heading))
                    print(content[start:end if end >= 0 else len(content)])
    else:
        raise ValueError(name)


SECTIONS = (
    ('M0 Existence', 'existence'), ('M1 Governing source text and label provenance', 'sources'),
    ('M2 Printed comparison scope', 'scope'), ('M3 Stored catalog', 'catalog'),
    ('M4 Row classes, joins and literal counts', 'shape'), ('M5 Joined-row paths and accounting', 'index'),
    ('M6 Every printed difference', 'differences'),
    ('M7 Oriented balance entries', 'inventory-limits'),
    ('M8 Closed results and reasons for unformed residuals', 'closed-results'),
    ('M9 Every unjoined path and its printed reason', 'unjoined'),
    ('M10 Filed accounting', 'accounting'), ('M11 Register before and after', 'register'),
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lookup', choices=[name for _, name in SECTIONS])
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        with OUTPUT.open('w') as target:
            target.write('# O2 record measurements\n\n'
                         'Generator: `scripts/O2_record_measurements.py`. Retrieval only; '
                         'no engine/comparator rerun, CAS, normalization or reconciliation.\n\n'
                         f'Regenerate from repository root: `python3 {SELF} --write`. '
                         'Cited source text and acceptance qualifications are retained literally. '
                         'JSON projections preserve stored values; whitespace alone is reformatted. '
                         'Line numbers refer to the filed comparison stream.\n\n')
            for title, name in SECTIONS:
                target.write(f'## {title}\n\n```text\n$ python3 {SELF} --lookup {name}\n')
                with contextlib.redirect_stdout(target):
                    lookup(name)
                target.write('\n```\n\n')
    elif args.lookup:
        lookup(args.lookup)
    else:
        parser.error('use --lookup or --write')


if __name__ == '__main__':
    main()
