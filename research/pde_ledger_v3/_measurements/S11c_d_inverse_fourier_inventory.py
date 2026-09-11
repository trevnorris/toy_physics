#!/usr/bin/env python3
"""Inventory carrier/source, inverse, branch and action reconstruction records."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from ledger_fold import _restore


def association(value):
    return {str(k): v for k, v in value}


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('transcript', type=Path)
    parser.add_argument('--focused', action='store_true')
    args = parser.parse_args()
    records = {}
    tags = Counter()
    selected = ('ELEMENT_CENSUS_', 'PROFILE_CARRIER_REDUCTION_', 'CARRIER_',
                'INVERSE_FOURIER_', 'RECONSTRUCTION_', 'BRANCH_RECONSTRUCTION_',
                'ACTION_INTEGRAL_REDUCTION_', 'BUILD_INPUT_DIGESTS', 'PROCESS_COMPLETION',
                'METADATA_CARRIER_', 'METADATA_INVERSE_FOURIER_', 'METADATA_RECONSTRUCTION_',
                'METADATA_BRANCH_RECONSTRUCTION_', 'OUTSTANDING_CONSTRUCTIONS',
                'REDUCED_DIMENSION_', 'INVERSE_CHECK_DIMENSION_')
    for line in args.transcript.open():
        tag, sep, body = line.partition(': ')
        if not sep:
            continue
        name = tag.removeprefix('PY_S11CD_')
        tags[name] += 1
        if name.startswith(selected):
            records[name] = body.strip()
    gaps = []
    expected, actual = defaultdict(set), defaultdict(set)
    slots = defaultdict(set)
    expected_integrals, expected_branches = {}, {}
    action_counts, branch_counts = Counter(), Counter()
    literal = defaultdict(lambda: {'objects': 0, 'entries': 0, 'nonzero': []})
    projections = defaultdict(lambda: {'objects': 0, 'entries': 0, 'maximum_absolute_projection': 0.})
    zero_units = []
    carrier_arguments = Counter()
    inverse_units = Counter()
    for name, raw in records.items():
        if name.startswith('ELEMENT_CENSUS_'):
            suffix = name.removeprefix('ELEMENT_CENSUS_')
            split = next((suffix.split('_'+slot, 1)[0], slot) for slot in
                         ('VALUE', 'MULTIGRADE', 'DIMENSION_L_T_M', 'COMPUTED_BRANCH_BINDINGS',
                          'FOURIER_PROFILE_BINDINGS') if suffix.endswith('_'+slot))
            case, slot = split
            slots[case].add(slot)
            value = association(_restore(raw))
            expected[case].update(sp.srepr(v) for v in value['hat_arguments'])
            if slot == 'VALUE':
                expected_integrals[case] = int(value['distinct_integrals'])
            if slot == 'COMPUTED_BRANCH_BINDINGS':
                expected_branches[case] = len(value['branch_arguments'])
        if name.startswith('CARRIER_INVERSE_ARGUMENT_OPERANDS_'):
            suffix = name.removeprefix('CARRIER_INVERSE_ARGUMENT_OPERANDS_')
            case = suffix.rsplit('_', 1)[0]
            node = _restore(raw)[0]
            actual[case].add(sp.srepr(node))
            carrier_arguments[node.func.__name__] += 1
            for prefix in ('CARRIER_INVERSE_SOURCE_OPERANDS_', 'CARRIER_INVERSE_REDUCED_OPERANDS_',
                           'CARRIER_INVERSE_VALUES_', 'CARRIER_INVERSE_RESIDUAL_', 'CARRIER_SOURCE_IMAGE_RESIDUAL_'):
                if prefix+suffix not in records:
                    gaps.append(('missing_carrier_record', prefix+suffix))
        for prefix in ('CARRIER_INVERSE_RESIDUAL_', 'CARRIER_SOURCE_IMAGE_RESIDUAL_',
                       'CARRIER_INVERSE_REMAINDER_CENSUS_',
                       'INVERSE_FOURIER_WEAK_KERNEL_CHARACTER_RESIDUAL_', 'BRANCH_RECONSTRUCTION_RESIDUAL_'):
            if name.startswith(prefix) or name == prefix.rstrip('_'):
                value = _restore(raw)
                leaves = list(value) if isinstance(value, sp.Tuple) else [value]
                entry = literal[prefix]
                entry['objects'] += 1
                entry['entries'] += len(leaves)
                if any(v != 0 for v in leaves):
                    entry['nonzero'].append(name)
        for prefix in ('RECONSTRUCTION_ROW_RESIDUAL_', 'RECONSTRUCTION_INTEGRAL_RESIDUAL_'):
            if name.startswith(prefix):
                packet = association(_restore(raw))
                values = [abs(complex(v)) for _, _, samples in packet['OBJECT_SHA_AND_NUMERIC_PIT'] for v in samples]
                entry = projections[prefix]
                entry['objects'] += 1
                entry['entries'] += len(values)
                entry['maximum_absolute_projection'] = max([entry['maximum_absolute_projection'], *values])
        if name.startswith(('CARRIER_', 'INVERSE_FOURIER_', 'RECONSTRUCTION_', 'BRANCH_RECONSTRUCTION_')):
            if 'METADATA_'+name not in records:
                gaps.append(('missing_metadata', name))
        if name.startswith(('METADATA_CARRIER_', 'METADATA_INVERSE_FOURIER_')):
            for path, entry in _restore(raw):
                dims = association(entry)['DIMENSION_L_T_M']
                if str(dims) == 'ZERO_MAP':
                    zero_units.append((name, str(path)))
                if name.startswith('METADATA_CARRIER_INVERSE_RESIDUAL_'):
                    inverse_units[str(tuple(dims))] += 1
        if name.startswith('ACTION_INTEGRAL_REDUCTION_'):
            suffix = name.removeprefix('ACTION_INTEGRAL_REDUCTION_')
            action_counts[suffix.rsplit('_', 1)[0]] += 1
            if 'RECONSTRUCTION_INTEGRAL_OPERANDS_'+suffix not in records:
                gaps.append(('missing_action_inverse', suffix))
        if name.startswith('BRANCH_RECONSTRUCTION_OPERANDS_'):
            suffix = name.removeprefix('BRANCH_RECONSTRUCTION_OPERANDS_')
            branch_counts[suffix.rsplit('_', 1)[0]] += 1
    if not args.focused:
        wanted = {row+'_'+alpha+'_'+rho for row in ('s11cc2ClosedSlabOperator', 's11cc2ClosedCouplingKernel')
                  for alpha in ('LAB_HELD', 'MATERIAL_ADVECTED') for rho in ('RHO4_CONSTANT', 'RHOBR_CONSTANT')}
        if set(expected) != wanted:
            gaps.append(('row_case_census', sorted(wanted-set(expected)), sorted(set(expected)-wanted)))
        for case in expected.keys() | actual.keys():
            if expected[case] != actual[case]:
                gaps.append(('carrier_census', case, sorted(expected[case]-actual[case]),
                             sorted(actual[case]-expected[case])))
            if len(slots[case]) != 5:
                gaps.append(('slot_census', case, sorted(slots[case])))
            if expected_integrals.get(case) != action_counts[case]:
                gaps.append(('action_measure_census', case, expected_integrals.get(case), action_counts[case]))
            if expected_branches.get(case) != branch_counts[case]:
                gaps.append(('branch_census', case, expected_branches.get(case), branch_counts[case]))
    duplicates = [name for name, count in tags.items() if count != 1]
    dimension_records = {name: raw for name, raw in records.items()
                         if name.startswith(('REDUCED_DIMENSION_', 'INVERSE_CHECK_DIMENSION_'))}
    output = {'transcript': str(args.transcript.resolve()),
              'bytes': args.transcript.stat().st_size,
              'sha256': hashlib.sha256(args.transcript.read_bytes()).hexdigest(),
              'unique_tags': len(tags), 'duplicate_tags': duplicates,
              'process_completion': records.get('PROCESS_COMPLETION'),
              'carrier_counts_by_case': {k: len(v) for k, v in actual.items()},
              'source_census_counts_by_case': {k: len(v) for k, v in expected.items()},
              'five_slot_counts': {k: len(v) for k, v in slots.items()},
              'action_integral_counts': dict(action_counts), 'branch_counts': dict(branch_counts),
              'carrier_types': dict(carrier_arguments),
              'literal_residuals': dict(literal), 'residual_projections': dict(projections),
              'inverse_residual_dimensions': dict(inverse_units), 'zero_map_metadata': zero_units,
              'coverage_gaps': gaps, 'dimension_records': dimension_records,
              'outstanding': records.get('OUTSTANDING_CONSTRUCTIONS')}
    print(json.dumps(output, indent=2))
    raise SystemExit(bool(gaps or duplicates or zero_units or not output['process_completion']))


if __name__ == '__main__':
    run()
