#!/usr/bin/env python3
"""Locate completed A6–A8 consumers; inspect bytes/source, never CAS values."""
import argparse
import ast
import hashlib
import json
from pathlib import Path

M = Path(__file__).resolve().parent
LEDGER = M.parent
CASES = tuple(a + '__' + d for a in ('LAB_HELD', 'MATERIAL_ADVECTED')
              for d in ('RHO4_CONSTANT', 'RHOBR_CONSTANT'))
seen = {}


def pin(path):
    path = Path(path)
    digest = hashlib.sha256()
    size = 0
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
            size += len(block)
    value = {'path': str(path), 'resolvedPath': str(path.resolve(strict=True)),
             'bytes': size, 'sha256': digest.hexdigest()}
    seen[str(path)] = value
    return value


def checkpoint(stem):
    path = M / ('S11c_d_' + stem + '_checkpoint.json')
    record = pin(path)
    data = json.loads(path.read_text())
    return data, record


def artifact(cp, name):
    recorded = cp['artifacts'][name]
    value = pin(Path(cp['runDirectory']) / name)
    if (value['sha256'], value['bytes']) != (recorded['sha256'], recorded['bytes']):
        raise ValueError(('changed artifact bytes', name))
    return {**value, 'manifestKey': name, 'scientificPacketRestored': False}


def source(cp, filename, functions):
    path = M / filename
    value = pin(path)
    if value['sha256'] != cp['sourceFiles'][str(path.relative_to(LEDGER))]:
        raise ValueError(('changed consumed source', filename))
    text = path.read_text()
    nodes = {n.name: n for n in ast.parse(text).body if isinstance(n, ast.FunctionDef)}
    return {**value, 'functions': [
        {'name': name, 'line': nodes[name].lineno, 'endLine': nodes[name].end_lineno,
         'sha256': hashlib.sha256(ast.get_source_segment(text, nodes[name]).encode()).hexdigest()}
        for name in functions]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    ends, ep = checkpoint('remaining_case_end_sources')
    binding, bp = checkpoint('remaining_case_profile_bindings')
    response, rp = checkpoint('remaining_case_profile_response')
    coordinate, cp = checkpoint('remaining_case_coordinate_sources')
    profiles, pp = checkpoint('profile_form')
    original_profile = artifact(profiles, 'profile-form.pickle')
    profile_input = artifact(binding, 'accepted-profile/profile-form.pickle')
    response_input = artifact(response, 'profile-inputs/accepted-profile/profile-form.pickle')
    if len({p['sha256'] for p in (original_profile, profile_input, response_input)}) != 1:
        raise ValueError('profile source route bytes differ')
    uniform = []
    for case in CASES:
        packet = artifact(ends, 'cases/' + case + '/uniform-source.pickle')
        inventory = {end: ends['endInventory'][case + '__' + end]
                     for end in ('REFERENCE', 'LEFT', 'RIGHT')}
        uniform.append({'case': case, 'packet': packet,
                        'savedFieldsFromConstructor': ['records[REFERENCE|LEFT|RIGHT].coupling[TH|HT]',
                                                      'differences[LEFT|RIGHT][TH|HT]'],
                        'completedProducerInventory': inventory})
    profile_cases = []
    for case in CASES:
        if case == 'LAB_HELD__RHO4_CONSTANT':
            item = {'case': case, 'wholeAcceptedProfile': original_profile}
        else:
            item = {'case': case,
                    'binding': artifact(binding, 'cases/' + case + '/case-binding.pickle'),
                    'response': artifact(response, 'cases/' + case + '/continuum/profile-form.pickle'),
                    'savedFieldsFromConstructor': ['moments', 'baselineShapes', 'alteredShapes']}
        profile_cases.append(item)
    sources = [
        source(ends, 'S11c_d_uniform_source.py', ['construct', 'emit_result']),
        source(ends, 'S11c_d_remaining_case_end_sources.py', ['main']),
        source(binding, 'S11c_d_remaining_case_profile_bindings.py', ['load', 'main']),
        source(response, 'S11c_d_remaining_case_profile_response.py', ['load', 'compare_case']),
        source(coordinate, 'S11c_d_coordinate_source.py', ['first_jet_mutation', 'construct']),
        source(coordinate, 'S11c_d_remaining_case_coordinate_sources.py', ['native_record_body', 'construct']),
    ]
    output = {'status': 'COMPLETED_CONTROL_CONSUMER_NAVIGATION_NOT_PHYSICS_ACCEPTANCE',
              'checkpoints': [ep, bp, rp, cp, pp], 'sources': sources,
              'A6': {'cases': uniform, 'disposition': 'All four triplet packet routes located; no new uniform construction needed merely to obtain these operands.'},
              'A7': {'originalProfile': original_profile, 'bindingInput': profile_input,
                     'responseInput': response_input, 'cases': profile_cases,
                     'disposition': 'Same saved moment packet bytes; completed binding source checks actual original/altered profiles and endpoints before reuse, and response source carries those fields.'},
              'A8': {'disposition': 'Completed literal first-jet mutations feed coordinate record shape fields with original/address/unit and saved involution joins. This does not supply the separately uncomputed density-factor closed-operator insertion.',
                     'coordinateCheckpoint': cp},
              'scientificPacketRestorations': 0, 'scientificCalls': 0,
              'independentPhysicsReview': False}
    # File identity only. No completed scientific comparison or consumer executes.
    for path, old in tuple(seen.items()):
        if pin(path) != old:
            raise ValueError(('input changed during navigation', path))
    with args.output.open('x') as stream:
        json.dump(output, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'output': str(args.output), 'sha256': pin(args.output)['sha256'],
                      'uniformCases': len(uniform), 'profileCases': len(profile_cases),
                      'scientificCalls': 0}))


if __name__ == '__main__':
    main()
