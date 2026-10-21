#!/usr/bin/env python3
"""Retrieve the stored S12 evidence; never import or execute a scientific engine.

Only verbatim text/JSON field retrieval, named-object shape and literal-match
counts are performed. No expression is evaluated, normalized or compared by
this generator. All residuals and comparison outcomes were already printed by
the accepted comparator. Run from any directory; --section prints one lookup.
Default output is the S12 record's measurements file. No scratch from another
job is read. The pre-pass register is pinned to the record-author baseline.
"""
from __future__ import annotations

import argparse
from collections import Counter
import io
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[3]
V3 = ROOT / 'research/pde_ledger_v3'
BASELINE = '560c6c51'
REPAIR_BASELINE = '8608bc71'
SCRIPT = 'research/pde_ledger_v3/scripts/S12_record_measurements.py'
OUTPUT = V3 / 'steps/_measurements/S12_record_measurements.md'
COMPARATOR = V3 / 'scripts/out/S12_cross_engine_comparator.out'


def dump(value):
    return json.dumps(value, ensure_ascii=False, separators=(',', ':'))


def source(out, path, first=1, last=None, revision=None):
    text = ((ROOT / path).read_text() if revision is None else subprocess.check_output(
        ['git', 'show', f'{revision}:{path}'], cwd=ROOT, text=True))
    lines = text.splitlines()
    print(f'FILE {path}' + (f' AT {revision}' if revision else ''), file=out)
    for i, line in enumerate(lines, 1):
        if i >= first and (last is None or i <= last):
            print(f'{i}: {line}', file=out)


def transcript(out, folder, dimensions=False):
    path = V3 / folder / 'out/S12_order_source.out'
    print('FILE ' + str(path.relative_to(ROOT)), file=out)
    # Continuation lines belong to their original tag, including Matrix rows.
    take = False
    for i, line in enumerate(path.read_text().splitlines(), 1):
        m = re.match(r'^([A-Za-z][A-Za-z0-9_/]*): ', line)
        if m:
            tag = m.group(1)
            take = ('dimension' in tag) == dimensions
            if tag == 'WORKING_ON' or tag.endswith('/comparison'):
                take = False
        if take:
            print(f'{i}: {line}', file=out)


def walk(value, path=()):
    """Enumerate stored JSON nodes, without interpreting their contents."""
    yield path, value
    if isinstance(value, dict):
        for key, child in value.items():
            yield from walk(child, path + (key,))
    elif isinstance(value, list):
        for i, child in enumerate(value):
            yield from walk(child, path + (i,))


def literal_counts(value, key):
    return dict(Counter(dump(node[key]) for _, node in walk(value)
                        if isinstance(node, dict) and key in node))


def comparison_projection(value):
    """Retrieve every nonsample field; name sample shape and literal counts."""
    if isinstance(value, dict):
        result = {}
        for key, child in value.items():
            if key == 'samples':
                result['stored_samples_count'] = len(child)
                result['stored_sample_residual_literal_counts'] = literal_counts(child, 'residual')
            else:
                result[key] = comparison_projection(child)
        return result
    if isinstance(value, list):
        return [comparison_projection(child) for child in value]
    return value


TITLES = {
    'M01': 'Step contract, premises, model point and retained OPEN',
    'M02': 'Acceptance dispositions and instrument notes',
    'M03': 'Both transcripts: nondimensional objects, verbatim',
    'M04': 'Both constructions: source definitions, verbatim',
    'M05': 'Comparator configuration, printed summary and record shape',
    'M06': 'Every accounting row, including all unmatched occurrences and leaf paths',
    'M07': 'Every comparison and correspondence, sample shape and literal counts',
    'M08': 'Every non-"0" residual node and every false text comparison, verbatim',
    'M09': 'Every explicit coverage gap and unpaired component, verbatim',
    'M10': 'Both transcripts: dimensional objects, verbatim',
    'M11': 'Separate inventories and their source rows',
    'M12': 'Register schema, population method, baseline and final census',
    'M13': 'Identified prior art and historical second-route disagreement',
    'M14': 'Filed build/comparator run records and comparator contract',
    'M15': 'Final S12 register pass and status lines',
    'M16': 'Literal residual/outcome counts across stored comparisons',
    'M17': 'Round-1 repair: owner handoffs, coefficient diagnostics and unit attribution',
}


def retrieve(section, rows):
    out = io.StringIO()
    if section == 'M01':
        source(out, 'research/pde_ledger_v3/directives/S12_decision_list.md')
    elif section == 'M02':
        for name in ['S12_build_review_r0_disposition.md', 'S12_build_review_r1_disposition.md',
                     'S12_build_review_r2_disposition.md',
                     'S12_comparator_build_review_r0_disposition.md',
                     'S12_comparator_build_review_r1_disposition.md',
                     'S12_comparator_build_review_r2_disposition.md',
                     'light_sector_closeout_review_r2_disposition.md',
                     'S12_record_decision_list_review_p1_disposition.md']:
            source(out, 'research/pde_ledger_v3/directives/_measurements/' + name)
    elif section in ('M03', 'M10'):
        for folder in ['scripts', 'mathematica']:
            transcript(out, folder, dimensions=section == 'M10')
    elif section == 'M04':
        for path in ['scripts/S12_order_source.py', 'mathematica/S12_order_source.wl']:
            source(out, 'research/pde_ledger_v3/' + path)
    elif section == 'M05':
        for r in rows:
            if r['kind'] in ['configuration', 'summary']:
                print(dump(r), file=out)
        print('stored record kind cardinalities ' + dump(dict(Counter(r['kind'] for r in rows))), file=out)
        print('WORKING_ON occurrence cardinalities ' + dump(dict(Counter(
            r['engine'] for r in rows if r['kind'] == 'accounting' and r['tag'] == 'WORKING_ON'))), file=out)
    elif section == 'M06':
        for r in rows:
            if r['kind'] == 'accounting':
                # Operands/parses are already preserved in the accepted output.
                print(dump({k: v for k, v in r.items() if k not in ['original_operand', 'parse']}), file=out)
    elif section == 'M07':
        for r in rows:
            if r['kind'] not in ['operands', 'accounting', 'residual_before_domain',
                                 'test_functions', 'probe_domain', 'resources', 'configuration', 'summary']:
                print(dump(comparison_projection(r)), file=out)
    elif section == 'M08':
        for r in rows:
            if r['kind'] != 'comparison':
                continue
            for path, node in walk(r['comparison']):
                if not isinstance(node, dict):
                    continue
                if 'as_computed_residual' in node and node['as_computed_residual'] != '0':
                    print(dump({'tag': r['tag'], 'path': path,
                                'stored_node': comparison_projection(node)}), file=out)
                    # Literal mismatch with the string "0", not a numerical test.
                    for sample in node.get('samples', []):
                        if sample.get('residual') != '0':
                            print(dump({'tag': r['tag'], 'path': path, 'stored_sample': sample}), file=out)
                if 'text_difference' in node and node['text_difference'].get('equal_text') is False:
                    print(dump({'tag': r['tag'], 'path': path, 'stored_node': node}), file=out)
    elif section == 'M09':
        for r in rows:
            if r['kind'] != 'comparison':
                continue
            for path, node in walk(r['comparison']):
                if isinstance(node, dict):
                    for key in ['not_compared', 'not_formed', 'unpaired_components',
                                'unresolved_coverage', 'py_coverage', 'wl_coverage', 'role', 'normalization']:
                        if key in node:
                            print(dump({'tag': r['tag'], 'path': path + (key,), 'stored_value': node[key]}), file=out)
    elif section == 'M11':
        for path, first, last in [
            ('research/pde_ledger_v3/V3_STEP_PLAN.md', 579, 637),
            ('research/pde_ledger_v2/notes/stage045_nonvariational_block_prep.md', 14, 60),
            ('docs/two_throat_simulation_handoff_spec.md', 338, 342),
            ('research/pde_ledger_v3/steps/O2_steady_brane_balance.md', 74, 91),
            ('research/pde_ledger_v3/steps/O2_steady_brane_balance.md', 389, 413),
            ('research/pde_ledger_v3/directives/O2_premise_decision_list.md', 28, 38),
            ('research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md', 129, 140),
            ('research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md', 236, 286),
            ('research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md', 578, 612),
            ('research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md', 72, 92),
            ('docs/oriented_throat_return_source_correction.md', 1, 28)]:
            source(out, path, first, last)
    elif section == 'M12':
        path = 'research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md'
        before = subprocess.check_output(['git', 'show', f'{BASELINE}:{path}'], cwd=ROOT, text=True)
        reviewed = subprocess.check_output(['git', 'show', f'{REPAIR_BASELINE}:{path}'], cwd=ROOT, text=True)
        after = (ROOT / path).read_text()
        for name, text in [('before ' + BASELINE, before), ('before repair ' + REPAIR_BASELINE, reviewed),
                           ('after working tree', after)]:
            headings = [line for line in text.splitlines() if line.startswith('### R-')]
            print(name + ': entry headings cardinality ' + str(len(headings)), file=out)
            print('\n'.join(headings), file=out)
        source(out, path, 42, 69, revision=BASELINE)
        source(out, path, 825, 840, revision=BASELINE)
        lines = after.splitlines()
        for heading in ['### Method for population passes', '### R-S1-01', '### R-S1-02',
                        '### R-S12-01', '### R-S12-02',
                        '### R-O2-01', '### R-S21-01', '### R-S22-01']:
            start = next(i for i, line in enumerate(lines) if line.startswith(heading))
            stop = next((i for i in range(start + 1, len(lines)) if lines[i].startswith('### ')), len(lines))
            source(out, path, start + 1, stop)
        # Both versions of every source-bearing entry, retrieved rather than
        # inferring a changed-source predicate from the text.
        for heading in ['### R-S1-01', '### R-S1-02', '### R-S12-02', '### R-S21-01']:
            original_lines = before.splitlines()
            start = next(i for i, line in enumerate(original_lines) if line.startswith(heading))
            stop = next((i for i in range(start + 1, len(original_lines))
                         if original_lines[i].startswith('### ')), len(original_lines))
            source(out, path, start + 1, stop, revision=BASELINE)
    elif section == 'M13':
        path = 'research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md'
        # These source blocks precede the inserted S12 pass, so line numbers are stable.
        source(out, path, 884, 920, revision=BASELINE)
        source(out, path, 1020, 1047, revision=BASELINE)
        source(out, path, 1126, 1176, revision=BASELINE)
        source(out, 'research/pde_ledger_v3/directives/O2_steady_brane_balance_scoping.md', 182, 208)
    elif section == 'M14':
        for path in ['directives/S12_comparator_decision_list.md',
                     'directives/_measurements/S12_build_run_record.md',
                     'directives/_measurements/S12_comparator_run_record.md',
                     'directives/_measurements/S12_comparator_builder_report.md',
                     'directives/_measurements/S12_comparator_production_run.md']:
            source(out, 'research/pde_ledger_v3/' + path)
        source(out, 'research/pde_ledger_v3/scripts/S12_cross_engine_comparator.py', 308, 325)
        source(out, 'research/pde_ledger_v3/scripts/S12_cross_engine_comparator.py', 820, 884)
        source(out, 'research/pde_ledger_v3/scripts/S12_cross_engine_comparator.py', 948, 1005)
        source(out, 'research/pde_ledger_v3/scripts/S12_cross_engine_comparator.py', 1040, 1186)
        source(out, 'research/pde_ledger_v3/scripts/S12_cross_engine_comparator.py', 1220, 1310)
    elif section == 'M15':
        path = 'research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md'
        lines = (ROOT / path).read_text().splitlines()
        source(out, path, 1, 22)
        start = next((i for i, l in enumerate(lines) if l.startswith('### S12 pass')), None)
        if start is not None:
            stop = next((i for i in range(start + 1, len(lines)) if lines[i].startswith('### ')), len(lines))
            source(out, path, start + 1, stop)
        source(out, 'research/pde_ledger_v3/steps/S9b_flow_oracle_comparison.md', 1, 8)
        source(out, 'STATUS.md')
    elif section == 'M16':
        comparisons = [r['comparison'] for r in rows if r['kind'] == 'comparison']
        for key in ['as_computed_residual', 'on_supplied_domain_residual', 'residual',
                    'equal_text', 'not_compared', 'not_formed', 'unresolved_coverage',
                    'exact_unresolved', 'py_coverage', 'wl_coverage']:
            print(key + ' ' + dump(literal_counts(comparisons, key)), file=out)
        for key in ['face_order', 'f_range', 'G_range']:
            print('probe_domain/' + key + ' ' + dump(dict(Counter(
                r['supplied_domain_values'][key] for r in rows if r['kind'] == 'probe_domain'))), file=out)
    elif section == 'M17':
        source(out, 'research/pde_ledger_v3/directives/_measurements/S12_record_r0_review_disposition.md')
        source(out, 'research/pde_ledger_v3/steps/S9b_flow_oracle_comparison.md', 151, 178)
        source(out, 'research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md', 14, 33)
        source(out, 'research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md', 821, 826, revision=REPAIR_BASELINE)
        source(out, 'research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md', 129, 140)
        source(out, 'research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md', 236, 248)
        source(out, 'research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md', 276, 282)
        source(out, 'research/pde_ledger_v3/directives/S12_comparator_decision_list.md', 51, 59)
        source(out, 'research/pde_ledger_v3/directives/S12_comparator_decision_list.md', 118, 119)
        tags = [f'{package}/face/{face}/coefficient_operand'
                for package in ['material', 'general_order', 'conditional_fast_ordering']
                for face in [1, 2]]
        for row in rows:
            if row.get('tag') not in tags:
                continue
            if row['kind'] == 'operands':
                print(dump({'kind': row['kind'], 'tag': row['tag'],
                            'py': {key: row['py'][key] for key in ['tag', 'line', 'raw']},
                            'wl': {key: row['wl'][key] for key in ['tag', 'line', 'raw']}}), file=out)
            elif row['kind'] == 'comparison':
                print(dump(comparison_projection(row)), file=out)
        source(out, 'research/pde_ledger_v3/scripts/S12_order_source.py', 219, 244)
        source(out, 'research/pde_ledger_v3/mathematica/S12_order_source.wl', 25, 42)
        source(out, 'research/pde_ledger_v3/scripts/out/S12_order_source.out', 35, 123)
        for folder in ['scripts', 'mathematica']:
            path = V3 / folder / 'out/S12_order_source.out'
            print(str(path.relative_to(ROOT)) + ' primitive_function_dimensions tag literal count ' +
                  str(sum(line.startswith('primitive_function_dimensions: ')
                          for line in path.read_text().splitlines())), file=out)
        for row in rows:
            if row['kind'] == 'accounting' and row['tag'] == 'primitive_function_dimensions':
                print(dump({key: value for key, value in row.items()
                            if key not in ['parse', 'original_operand']}), file=out)
    else:
        raise ValueError(section)
    return out.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--section', choices=TITLES)
    args = parser.parse_args()
    rows = [json.loads(line) for line in COMPARATOR.read_text().splitlines()]
    if args.section:
        print(retrieve(args.section, rows), end='')
        return
    output = ['# S12 record measurements\n\n',
              'Generated by `' + SCRIPT + '`. Run `python3 ' + SCRIPT + '` to regenerate.\n\n',
              'Commands below run from the repository root. Their literal stdout is retained. ',
              'JSON values are retrieved from the accepted production output, never recomputed. ',
              'A literal other than `"0"` can be the Boolean residual `"False"`; no numerical ',
              'or logical equivalence is inferred by this generator. Sample counts count stored ',
              'nodes, not independent tests. Unmatched operands remain in the original annex output.\n\n']
    for section, title in TITLES.items():
        output += [f'## {section} — {title}\n\n',
                   f'```text\n$ python3 {SCRIPT} --section {section}\n',
                   retrieve(section, rows), '```\n\n']
    OUTPUT.write_text(''.join(output))


if __name__ == '__main__':
    main()
