#!/usr/bin/env python3
"""Index existing JSON/source bytes only; never restore or execute scientific objects."""
import ast
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path('/var/projects/toy_physics')
M = ROOT / 'research/pde_ledger_v3/_measurements'
SC = ROOT / '_scratch/s11c/s11c-defect-source-composition-20261002/diagnostic-01/complete'
OUTPUT = M / 'S11c_d_defect_packet_pressure_transport_source_map.json'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read(path):
    return json.loads(path.read_text())


def receipt(path):
    data = path.read_bytes()
    return {'path': str(path), 'sha256': digest(data), 'bytes': len(data)}


def main():
    if OUTPUT.exists():
        raise FileExistsError('Preserve existing preparation; do not overwrite it.')
    prior_path = M / 'S11c_d_defect_packet_pressure_unit_source_map.json'
    prior = read(prior_path)
    files = {}

    def add(path):
        key = str(path.relative_to(ROOT))
        files.setdefault(key, receipt(path))
        return key

    # The existing map is an index, not an independent scientific certificate.
    # Verify its byte receipts without executing any encoded expression.
    for r in prior['records'].values():
        got = receipt(Path(r['path']))
        assert (got['sha256'], got['bytes']) == (r['sha256'], r['bytes'])
    selected_path = Path(prior['records']['preflight-input/selected/pressure-addresses.json']['path'])
    selected = read(selected_path)['selected']
    selected_key = add(selected_path)
    assert len(selected) == 544 and len({a['addressId'] for a in selected}) == 544
    transforms_key = add(SC / 'whole-coefficient-transforms.json')
    transforms = read(ROOT / transforms_key)

    profile_records = {}
    for p in sorted(SC.glob('profile-map-*.json')):
        profile_records[add(p)] = read(p)
    assert len(profile_records) == 392

    grade_groups = {}
    names = ['plus-source', 'minus-source'] + [
        'THETA_BALANCE-' + slot + face
        for slot in ('delta_p_', 'd_w_delta_p_') for face in ('plus', 'minus')]
    for name in names:
        # Full denominator certificates include saved operands and their returns;
        # a quotient-remainder zero by itself is not a regularity certificate.
        paths = [SC / (name + '-operands.json'), SC / (name + '-split.json')]
        for suffix in ('regular-denominator*', 'saved-zero-*', 'quotient-*'):
            paths.extend(sorted(SC.glob(name + '-' + suffix + '.json')))
        grade_groups[name] = [add(p) for p in dict.fromkeys(paths)]
        assert any('regular-denominator-input' in p for p in grade_groups[name])
    source_tables = {}
    for face in ('plus', 'minus'):
        for grade in ('00', '10', '01', '11'):
            p = SC / (face + '-source-jets-' + grade + '.json')
            key = add(p)
            source_tables[face, grade] = (key, read(p))
            for suffix in ('input', 'raw', 'return'):
                add(SC / (face + '-source-jet-reconstruction-' + grade + '-' + suffix + '.json'))

    source_locations, consumer_locations, address_locations = {}, {}, []
    for i, a in enumerate(selected):
        face = a['face']
        sg = ''.join(map(str, a['sourceGrade']))
        source_key, source_table = source_tables[face, sg]
        matches = [j for j, item in enumerate(source_table['jets']) if
                   item['spec'] == a['jet'] and item['atom'] == a['sourceAtom'] and
                   item['coefficient'] == a['sourceOriginal'] and item['field'] == a['sourceField']]
        assert len(matches) == 1
        j = matches[0]
        source_id = face + '/' + sg + '/' + str(j)
        if source_id not in source_locations:
            pm = [p for p, record in profile_records.items() if
                  record['original'] == a['sourceOriginal'] and record['oneDimensional'] == a['sourceField']]
            assert pm
            source_locations[source_id] = {
                'record': source_key, 'jsonPointer': '/jets/' + str(j),
                'profileCandidateRecords': pm,
                'matchKind': 'Complete encoded JSON equality only; duplicates preserved.',
                'unitStatus': 'Required coefficient dimension is historical metadata, not independently verified here.'}
        slot = ('delta_p_' if a['slot'] == 'pressure' else 'd_w_delta_p_') + face
        name = 'THETA_BALANCE-' + slot
        cg = str(tuple(a['consumerGrade']))
        split_path = SC / (name + '-split.json')
        assert read(split_path)['retained'][cg] == a['consumerOriginal']
        consumer_id = slot + '/' + cg
        if consumer_id not in consumer_locations:
            pm = [p for p, record in profile_records.items() if record['oneDimensional'] == a['consumerField']]
            assert pm
            consumer_locations[consumer_id] = {
                'record': add(split_path), 'jsonPointer': '/retained/' + cg,
                'profileOutputCandidateRecords': pm,
                'matchKind': 'Encoded output equality only; epsilon division and input mapping remain runtime obligations.',
                'originalToProfileInputCertified': False}
        field_ids = []
        for role in ('source', 'consumer'):
            field_id = a[role + 'Transform']['coefficientId']
            assert transforms[field_id]['field'] == a[role + 'Field']
            assert digest(a[role + 'Field']['srepr'].encode()) == field_id
            field_ids.append(field_id)
        address_locations.append({'addressId': a['addressId'], 'selectedJsonPointer': '/selected/' + str(i),
                                  'source': source_id, 'consumer': consumer_id, 'fieldIds': field_ids})

    worker = M / 'S11c_d_defect_source_composition.py'
    source = worker.read_text()
    fragments = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.FunctionDef) and node.name in (
                'grade_split', 'quotient_recurrence', 'profile', 'wave', 'transform', 'join_source_input'):
            text = ast.get_source_segment(source, node)
            fragments.append({'name': node.name, 'line': node.lineno, 'endLine': node.end_lineno,
                              'source': add(worker), 'text': text, 'sha256': digest(text.encode()), 'executed': False})
    result = {
        'status': 'METADATA_ONLY_PRESSURE_TRANSPORT_EVIDENCE_INDEX_NOT_A_CERTIFICATE',
        'createdUtc': datetime.now(timezone.utc).isoformat(),
        'preparationScript': receipt(Path(__file__)), 'priorSourceMap': receipt(prior_path),
        'originalCompositionResultCommit': '7455ee78',
        'nativeUnitBuildReviewPendingSeparately': 'native-units-build-review',
        'records': files, 'selectedRecord': selected_key, 'transformRecord': transforms_key,
        'gradeGroups': grade_groups, 'sourceLocations': source_locations,
        'consumerLocations': consumer_locations, 'addressLocations': address_locations,
        'sourceFragments': fragments,
        'counts': {'priorReceiptsChecked': len(prior['records']), 'newIndexedFiles': len(files),
                   'profileRecords': len(profile_records), 'gradeGroups': len(grade_groups),
                   'sourceLocations': len(source_locations), 'consumerLocations': len(consumer_locations),
                   'addresses': len(address_locations), 'addressStatuses': dict(Counter(a['status'] for a in selected))},
        'obligations': [
            'Join accepted unbound native units to the SAME live-density/chemical/velocity/stage2/epsilon bindings.',
            'Transport each unit through independent dimensionless grades using complete saved rational quotient and denominator evidence; do not recompute grades.',
            'Join full native source-jet coefficient and profile inputs to saved fields. Candidate JSON matches are locators, not symbolic or dimensional proof.',
            'Join consumer epsilon removal to actual profile inputs; output equality alone is insufficient.',
            'Carry physical W and native L-scaled profile derivatives before numeric bindings; do not assign physical units to numeric constants.',
            'Retain units inferred from original gamma registries as inference-dependent. Zero coefficients have expected units only.',
            'Join every complete template, H/J/direct/contact/PV and Fourier measure with both normal signs and whole direct/iteration exactly once.',
            'Restore scientific operands only under the unchanged pooled guard. No completed functions or inventories may replay.'],
        'scientificRestoration': False, 'newAlgebraOrUnitsComputed': False,
        'priorFunctionsReplayed': False, 'newExternalSubmission': False, 'readyGate': False,
        'pressureValueComputed': False}
    with OUTPUT.open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps({'path': str(OUTPUT), 'bytes': OUTPUT.stat().st_size, 'counts': result['counts']}))


if __name__ == '__main__':
    main()
