#!/usr/bin/env python3
"""Retrieve the filed S14a part 1 evidence; perform no CAS or new physics check.

Only text/JSON retrieval, named stored-object shapes and literal-match counts
are used. No engine/comparator imports, evaluation, simplification or reruns.
Working files, if any, belong under _scratch/s14a/rec/work/.
"""
from pathlib import Path
from collections import Counter
import argparse
import io
import json
import subprocess

ROOT = Path(__file__).resolve().parents[3]
V3 = 'research/pde_ledger_v3/'
SCRIPT = V3 + 'scripts/S14a_1_record_measurements.py'
BASE = 'f5d283c1'
OUTPUT = ROOT / V3 / 'steps/_measurements/S14a_1_record_measurements.md'
COMPARATOR = ROOT / V3 / 'scripts/out/S14a_1_cross_engine_comparator.out'
TITLES = {
    'M01': 'Contracts, dated premises, model point and plan',
    'M02': 'Acceptance dispositions and carried limitations',
    'M03': 'Both accepted transcripts, verbatim with line numbers',
    'M04': 'Accepted engine constructions and comparator scope/accounting source',
    'M05': 'Printed configuration, tables, summary and record-kind counts',
    'M06': 'Every tag and binding accounting row, including all unpaired components',
    'M07': 'Every comparison, with stored sample shapes and literal residual counts',
    'M08': 'Every nonliteral-zero/False exact residual and every unequal text pair',
    'M09': 'Every explicit comparison coverage field and endpoint ledger',
    'M10': 'Literal outcome counts, note differences and named binding references',
    'M11': 'Separate inventory source rows and earlier projected material balance',
    'M12': 'Register schema, method, baseline/final census and source-bearing entries',
    'M13': 'Second-route prior-art sources and inherited handoffs',
    'M14': 'Run records, builder report and comparator decision list',
    'M15': 'S12 acceptance, J1/J2 sources, baseline and final affected text',
    'M16': 'Final S14a register pass, status and exact allowed ancillary diff',
    'M17': 'Round-0 repair sources and stored sample cardinalities by form',
}


def dump(value):
    return json.dumps(value, ensure_ascii=False, separators=(',', ':'))


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT, text=True)


def source(out, path, first=1, last=None, revision=None):
    value = git('show', revision + ':' + path) if revision else (ROOT / path).read_text()
    print('FILE ' + path + (' AT ' + revision if revision else ''), file=out)
    for number, line in enumerate(value.splitlines(), 1):
        if number >= first and (last is None or number <= last):
            print(f'{number}: {line}', file=out)


def block(out, path, heading, revision=None):
    value = git('show', revision + ':' + path) if revision else (ROOT / path).read_text()
    lines = value.splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith(heading))
    end = next((i for i in range(start + 1, len(lines)) if lines[i].startswith('### ')), len(lines))
    source(out, path, start + 1, end, revision)


def walk(value, path=()):
    yield path, value
    if isinstance(value, dict):
        for key, child in value.items():
            yield from walk(child, path + (key,))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from walk(child, path + (index,))


def literal_counts(value, key):
    return dict(Counter(dump(node[key]) for _, node in walk(value)
                        if isinstance(node, dict) and key in node))


def sample_projection(value):
    """Keep every nonsample field; count literal stored sample outcomes only."""
    if isinstance(value, dict):
        result = {}
        for key, child in value.items():
            if key == 'samples':
                result['stored_samples_count'] = len(child)
                for field in ['residual', 'value_coverage', 'exclusive',
                              'not_evaluated', 'unresolved_coverage']:
                    result['stored_sample_' + field + '_literal_counts'] = literal_counts(child, field)
            else:
                result[key] = sample_projection(child)
        return result
    if isinstance(value, list):
        return [sample_projection(child) for child in value]
    return value


def retrieve(section, rows):
    out = io.StringIO()
    if section == 'M01':
        for name in ['S14a_1_record_decision_list.md', 'S14a_1_decision_list.md']:
            source(out, V3 + 'directives/' + name)
        source(out, V3 + 'directives/S12_decision_list.md')
        source(out, V3 + 'V3_STEP_PLAN.md', 579, 646)
        source(out, V3 + 'V3_STEP_PLAN.md', 652, 668)
        source(out, V3 + 'steps/S9b_light_bending.md', 74, 84)
    elif section == 'M02':
        names = [f'S14a_1_build_review_r{i}_disposition.md' for i in range(3)]
        names += [f'S14a_1_comparator_build_review_r{i}_disposition.md' for i in range(2)]
        names += [f'S14a_1_decision_list_review_p{i}_disposition.md' for i in range(1, 4)]
        names += ['S14a_1_record_decision_list_review_p1_disposition.md', 'S12_record_r3_review_disposition.md']
        for name in names:
            source(out, V3 + 'directives/_measurements/' + name)
        print(git('show', '-s', '--format=%H %s', 'ac909aaf', 'e9eed614', '94fddf2f', 'de81b9fc', BASE), end='', file=out)
    elif section == 'M03':
        for folder in ['scripts', 'mathematica']:
            source(out, V3 + folder + '/out/S14a_1_projection.out')
    elif section == 'M04':
        for path in ['scripts/S14a_1_projection.py', 'mathematica/S14a_1_projection.wl']:
            source(out, V3 + path)
        for first, last in [(1, 9), (280, 440), (1198, 1410), (1420, 1595)]:
            source(out, V3 + 'scripts/S14a_1_cross_engine_comparator.py', first, last)
    elif section == 'M05':
        for row in rows:
            if row['kind'] in ['configuration', 'summary', 'resources', 'test_functions', 'probe_domain']:
                print(dump(row), file=out)
        print('stored record-kind cardinalities ' + dump(dict(Counter(row['kind'] for row in rows))), file=out)
    elif section == 'M06':
        for row in rows:
            if row['kind'] in ['tag_accounting', 'binding_accounting', 'progress_marker']:
                # Original operands/parses remain in the annex and M03.
                print(dump({k: v for k, v in row.items() if k not in ['raw', 'parse']}), file=out)
    elif section == 'M07':
        for row in rows:
            if row['kind'] in ['comparison', 'residual_before_domain']:
                print(dump(sample_projection(row)), file=out)
    elif section == 'M08':
        for row in rows:
            if row['kind'] != 'comparison':
                continue
            for path, node in walk(row['comparison']):
                if not isinstance(node, dict):
                    continue
                if 'as_computed_residual' in node and node['as_computed_residual'] not in ['0', 'False']:
                    print(dump({'binding': row['binding'], 'path': path, 'stored_node': sample_projection(node)}), file=out)
                    for sample in node.get('samples', []):
                        if sample.get('residual') not in ['0', 'False']:
                            print(dump({'binding': row['binding'], 'path': path, 'stored_sample': sample}), file=out)
                if 'text_difference' in node and node['text_difference'].get('equal_text') is False:
                    print(dump({'binding': row['binding'], 'path': path, 'stored_node': node}), file=out)
    elif section == 'M09':
        fields = ['not_compared', 'not_formed', 'unpaired_components', 'unresolved_coverage',
                  'exact_unresolved', 'py_coverage', 'wl_coverage', 'closed_endpoint_integrals',
                  'closed_measure_operands', 'domain_context_census', 'not_evaluated', 'value_coverage']
        for row in rows:
            if row['kind'] != 'comparison':
                continue
            for path, node in walk(row['comparison']):
                if isinstance(node, dict):
                    for key in fields:
                        if key in node:
                            print(dump({'binding': row['binding'], 'path': path + (key,), 'stored_value': node[key]}), file=out)
    elif section == 'M10':
        comparisons = [row['comparison'] for row in rows if row['kind'] == 'comparison']
        for key in ['as_computed_residual', 'on_supplied_domain_residual', 'residual', 'equal_text',
                    'not_compared', 'not_formed', 'unresolved_coverage', 'exact_unresolved',
                    'not_evaluated', 'value_coverage', 'exclusive']:
            print(key + ' ' + dump(literal_counts(comparisons, key)), file=out)
        for row in rows[0]['binding_table']:
            print(dump({key: row[key] for key in ['key', 'object', 'py', 'wl', 'condition_context']}), file=out)
        for key in ['face_order', 'shell_order', 'radius_positive']:
            print('probe_domain/' + key + ' ' + dump(dict(Counter(row[key] for row in rows if row['kind'] == 'probe_domain'))), file=out)
        print('tag accounting classes ' + dump(dict(Counter(row['accounting'] for row in rows if row['kind'] == 'tag_accounting'))), file=out)
        print('unjoined reason basis literal counts ' + dump(dict(Counter(row['reason_basis'] for row in rows if row['kind'] == 'tag_accounting' and row['accounting'] == 'unjoined'))), file=out)
        print('unjoined full-reason literal counts ' + dump(dict(Counter(row['reason'] for row in rows if row['kind'] == 'tag_accounting' and row['accounting'] == 'unjoined'))), file=out)
    elif section == 'M11':
        for path, first, last in [
            ('software/stage1_solver/reports/pathA_21_emergent_G_mass_bridge.md', 14, 22),
            ('software/stage1_solver/reports/pathA_21_emergent_G_mass_bridge.md', 62, 71),
            ('software/stage1_solver/reports/pathA_21_emergent_G_mass_bridge.md', 97, 104),
            (V3 + 'DEFECT_REGISTER.md', 167, 174),
            ('research/pde_ledger_v2/notes/stages/ledger_stage005_sound_speed_light_ratio.md', 216, 220),
            ('research/pde_ledger_v2/notes/parameter_register.md', 147, 154),
            ('research/pde_ledger_v2/notes/parameter_register.md', 235, 235),
            ('research/pde_ledger_v2/notes/stages/ledger_stage004_gnls_action_dimensional_foundation.md', 139, 150),
            ('research/pde_ledger_v2/notes/stages/ledger_stage002_matter_stress_force_assembly.md', 51, 55),
            ('research/pde_ledger_v2/notes/stages/ledger_stage002_matter_stress_force_assembly.md', 109, 113),
            ('research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md', 78, 87),
            (V3 + 'directives/O2_SHARED_PHYSICS.md', 480, 494),
            (V3 + 'V3_STEP_PLAN.md', 619, 641),
            (V3 + 'V3_STEP_PLAN.md', 663, 665)]:
            source(out, path, first, last)
    elif section == 'M12':
        path = V3 + 'SUBSTRATE_REQUIREMENTS.md'
        for revision, value in [(BASE, git('show', BASE + ':' + path)), ('working tree', (ROOT / path).read_text())]:
            headings = [line for line in value.splitlines() if line.startswith('### R-')]
            print(revision + ': entry headings cardinality ' + str(len(headings)), file=out)
            print(revision + ': literal **status** OPEN count ' + str(value.count('**status** OPEN')), file=out)
            print('\n'.join(headings), file=out)
        source(out, path, 1, 115, BASE)
        block(out, path, '### Method for population passes', BASE)
        for heading in ['### R-S1-01', '### R-S1-02', '### R-S12-02', '### R-S21-01']:
            for revision in [BASE, None]:
                block(out, path, heading, revision)
        for heading in ['### R-S12-01', '### R-S22-01', '### R-O2-01']:
            block(out, path, heading, BASE)
    elif section == 'M13':
        source(out, V3 + 'directives/O2_steady_brane_balance_scoping.md', 182, 208)
        for heading in ['### O2 pass', '### S9b pass', '### Light-closing pass']:
            block(out, V3 + 'SUBSTRATE_REQUIREMENTS.md', heading, BASE)
        source(out, V3 + 'steps/S9b_flow_oracle_comparison.md', 151, 178)
        source(out, V3 + 'steps/S12_order_source.md', 360, 439, BASE)
    elif section == 'M14':
        for path in ['directives/S14a_1_comparator_decision_list.md',
                     'directives/_measurements/S14a_1_build_runs.md',
                     'directives/_measurements/S14a_1_comparator_run_record.md',
                     'directives/_measurements/S14a_1_comparator_builder_report.md',
                     'directives/_measurements/S14a_1_comparator_production_run.md']:
            source(out, V3 + path)
    elif section == 'M15':
        source(out, V3 + 'directives/_measurements/S12_record_r3_review_disposition.md')
        source(out, V3 + 'steps/_measurements/S12_record_measurements.md', 1, 5)
        # Retrieve the section holding the historical round-3 count, verbatim.
        path = V3 + 'steps/_measurements/S12_record_measurements.md'
        lines = (ROOT / path).read_text().splitlines()
        first = next(i for i, line in enumerate(lines) if line.startswith('## M12 '))
        last = next(i for i in range(first + 1, len(lines)) if lines[i].startswith('## '))
        source(out, path, first + 1, last)
        for revision in [BASE, None]:
            source(out, V3 + 'steps/S12_order_source.md', revision=revision)
            block(out, V3 + 'SUBSTRATE_REQUIREMENTS.md', '### S12 pass', revision)
        source(out, V3 + 'steps/S9b_flow_oracle_comparison.md', 1, 8)
    elif section == 'M16':
        block(out, V3 + 'SUBSTRATE_REQUIREMENTS.md', '### S14a part 1 pass')
        source(out, V3 + 'SUBSTRATE_REQUIREMENTS.md', 1, 28)
        source(out, 'STATUS.md')
        paths = [V3 + 'SUBSTRATE_REQUIREMENTS.md', V3 + 'steps/S12_order_source.md',
                 V3 + 'directives/_measurements/S14a_1_build_runs.md', 'STATUS.md']
        print('$ git diff ' + BASE + ' --name-only', file=out)
        print(git('diff', BASE, '--name-only'), end='', file=out)
        print('The pre-existing .claude/skills/prep-for-compact/SKILL.md modification is outside this work; its contents are not retrieved.', file=out)
        print('$ git diff ' + BASE + ' -- ' + ' '.join(paths), file=out)
        print(git('diff', BASE, '--', *paths), end='', file=out)
    elif section == 'M17':
        source(out, V3 + 'directives/_measurements/S14a_1_record_r0_review_disposition.md')
        print('$ git show --format=%h --name-only 25cb3a86 -- research/pde_ledger_v3/directives/_measurements/S14a_1_record_r0_review_disposition.md research/pde_ledger_v3/directives/_measurements/S14a_1_record_r0_review_lookups.md', file=out)
        print(git('show', '--format=%h', '--name-only', '25cb3a86', '--',
                  V3 + 'directives/_measurements/S14a_1_record_r0_review_disposition.md',
                  V3 + 'directives/_measurements/S14a_1_record_r0_review_lookups.md'), end='', file=out)
        source(out, V3 + 'directives/_measurements/S14a_1_comparator_build_review_r0_disposition.md', 54, 67)
        source(out, V3 + 'directives/_measurements/S14a_1_comparator_build_review_r1_disposition.md', 90, 100)
        for row in rows:
            if row['kind'] != 'comparison' or row['binding'] not in [
                    'b50', 'b72', 'b128', 'b148', 'b171', 'b173', 'b174']:
                continue
            all_samples = []
            for index, part in enumerate(row['comparison']['components']):
                samples = part['samples']
                all_samples.extend(samples)
                print(dump({'binding': row['binding'], 'component_index': index,
                            'part': part['part'], 'stored_samples_count': len(samples),
                            'stored_sample_form_literal_counts': dict(Counter(sample['form'] for sample in samples)),
                            'stored_sample_residual_literal_counts_by_form': {
                                form: literal_counts([sample for sample in samples if sample['form'] == form], 'residual')
                                for form in ['affine', 'quadratic', 'cubic']}}), file=out)
            print(dump({'binding': row['binding'], 'all_components_stored_samples_count': len(all_samples),
                        'all_components_stored_sample_form_literal_counts': dict(Counter(
                            sample['form'] for sample in all_samples))}), file=out)
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
    result = ['# S14a part 1 record measurements\n\n',
              'Generated by `' + SCRIPT + '`. Commands run from the repository root.\n\n',
              'This generator retrieves filed evidence only. JSON sample counts are literal ',
              'counts of stored nodes, not new evaluations or independent tests. An exact string ',
              'other than `0` or `False` is retained without a new equivalence decision. ',
              'M03 holds both complete transcripts; M05 the full reference tables; M06 every ',
              'unjoined tag and unpaired component; M07 every compared part; M08 all printed ',
              'exact/text differences and nonzero sample outcomes; M09 every explicit coverage field. ',
              'Full original JSON operands and sample values remain in the annex comparator output.\n\n']
    for section, title in TITLES.items():
        result += [f'## {section} — {title}\n\n',
                   f'```text\n$ python3 {SCRIPT} --section {section}\n', retrieve(section, rows), '```\n\n']
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(''.join(result))


if __name__ == '__main__':
    main()
