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
    'scripts/O2_live_balance_sympy_audit.py',
    'mathematica/O2_live_balance_mathematica_audit.wl',
    'scripts/O2_cross_engine_comparator.py',
]


def dump(obj):
    print(json.dumps(obj, ensure_ascii=False, separators=(',', ':')))


def rows():
    with STREAM.open() as source:
        for number, line in enumerate(source, 1):
            yield number, json.loads(line), line


def fields(obj, keys):
    return {key: obj[key] for key in keys if key in obj}


def retrieve_comparison(obj, path=()):
    """Retrieve all leaf results, including every stored structural delta."""
    if 'children' in obj:
        dump({'comparison_path': path, **fields(obj, ('outcome',))})
        for key, child in obj['children'].items():
            retrieve_comparison(child, (*path, key))
    else:
        dump({'comparison_path': path, **obj})


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
                print(f'{number}: {line}')
    elif name == 'catalog':
        # Complete stored catalog, not a newly inferred join map.
        print((V3 / 'scripts/out/O2_cross_engine_comparator_catalog.json').read_text(), end='')
    elif name == 'accounting':
        print((V3 / 'scripts/out/O2_cross_engine_comparator_accounting.jsonl').read_text(), end='')
    elif name == 'scope':
        for number, obj, line in rows():
            if obj['kind'] == 'comparison_scope':
                print(f'STREAM_LINE {number}')
                print(line, end='')
    elif name == 'index':
        for number, obj, _ in rows():
            dump({'stream_line': number, **fields(obj, (
                'kind', 'row', 'component', 'engine', 'tag', 'path',
                'py_path', 'wl_path', 'accounting', 'reason',
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
                print(f"STREAM_LINE {number} ROW {obj['row']}")
                retrieve_comparison(obj['comparison'])
            elif obj['kind'] == 'balance_comparison':
                dump({'stream_line': number, **obj})
            elif obj['kind'] == 'unjoined':
                dump({'stream_line': number, **fields(obj, (
                    'kind', 'accounting', 'engine', 'path', 'reason',
                    'parsed_leaves', 'compared_leaves'))})
    elif name == 'inventory-limits':
        for number, obj, _ in rows():
            if obj['kind'] == 'structure':
                for action in obj['action_comparison']:
                    print(f"STREAM_LINE {number} ROW {obj['row']} ROLE {action['role']}")
                    dump(fields(action, ('py', 'wl', 'occurrences')))
            elif obj['kind'] == 'balance_entries':
                print(f"STREAM_LINE {number} ROW {obj['row']} COMPONENT {obj['component']}")
                for engine, entries in obj['entries'].items():
                    for entry in entries:
                        dump({'engine': engine, **fields(entry, (
                            'role', 'orientation', 'open_free', 'inventories', 'limit'))})
    elif name == 'closed-results':
        # All stored outcome/reason/residual leaves, not only successful rows.
        def retrieve(obj, path=()):
            dump({'comparison_path': path, **fields(obj, (
                'outcome', 'reason', 'residual', 'py_value', 'wl_value',
                'point_policy', 'points', 'closed_residual'))})
            for key, child in obj.get('children', {}).items():
                retrieve(child, (*path, key))
            for index, child in enumerate(obj.get('operands', [])):
                retrieve(child, (*path, 'operands', index))
        for number, obj, _ in rows():
            if obj['kind'] in ('residual', 'balance_comparison'):
                print(f"STREAM_LINE {number} ROW {obj['row']} COMPONENT {obj.get('component', [])}")
                retrieve(obj['comparison'])
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
        for kind in ('comparison_scope', 'operands', 'structure', 'residual', 'accounting',
                     'balance_operands', 'balance_entries', 'balance_comparison',
                     'unjoined', 'stream_object_count'):
            with STREAM.open() as source:
                count = sum(line.count('"kind":"' + kind + '"') for line in source)
            dump({'literal': '"kind":"' + kind + '"', 'count': count})
    elif name == 'register':
        before = subprocess.run(
            ['git', 'show', BASELINE + ':' + REGISTER], cwd=ROOT,
            check=True, text=True, capture_output=True).stdout
        after = (ROOT / REGISTER).read_text()
        for label, content in ((BASELINE, before), ('working-tree', after)):
            print('REGISTER ' + label)
            dump({'literal': '### R-', 'count': content.count('### R-')})
            for line in content.splitlines():
                if line.startswith('### R-') or line.startswith('- **source**'):
                    print(line)
        print('REGISTER_BASELINE_FULL_TEXT')
        print(before, end='')
        print('\nREGISTER_WORKING_TREE_FULL_TEXT')
        print(after, end='')
    else:
        raise ValueError(name)


SECTIONS = (
    ('M0 Existence', 'existence'), ('M1 Source text and engine/comparator sources', 'sources'),
    ('M2 Printed comparison scope', 'scope'), ('M3 Stored catalog', 'catalog'),
    ('M4 Row classes, joins and literal counts', 'shape'), ('M5 Complete row index', 'index'),
    ('M6 Every printed difference and residual', 'differences'),
    ('M7 Occurrence inventories and their printed limits', 'inventory-limits'),
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
                         'The source headers and acceptance qualifications are retained literally. '
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
