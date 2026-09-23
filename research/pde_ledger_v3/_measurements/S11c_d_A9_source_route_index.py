#!/usr/bin/env python3
"""Index existing A9 inputs as bytes/source text; never import a science module.

This is navigation, not validation of mathematical values. Checkpoint hashes of
pickles are inherited metadata. Transcript ranges are freshly byte-checked.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path

LEDGER = Path(__file__).resolve().parents[1]
PROJECT = LEDGER.parents[1]
M = LEDGER / '_measurements'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def info(path):
    return {'path': str(path), 'resolvedPath': str(path.resolve(strict=True)),
            'sha256': digest(path), 'bytes': path.stat().st_size}


def checkpoint(stem, select):
    path = M / ('S11c_d_' + stem + '_checkpoint.json')
    data = json.loads(path.read_text())
    run = Path(data['runDirectory'])
    chosen = []
    for name, item in data['artifacts'].items():
        if not select(name):
            continue
        logical = run / name
        resolved = logical.resolve(strict=True)
        size = resolved.stat().st_size
        if size != item['bytes']:
            raise ValueError(('changed saved size', str(logical)))
        chosen.append({'manifestKey': name, 'logicalPath': str(logical),
                       'resolvedPath': str(resolved), 'bytes': size,
                       'manifestSha256': item['sha256'],
                       'payloadHashRecomputed': False, 'payloadRestored': False})
    if not chosen:
        raise ValueError(('empty selection', stem))
    return {'checkpoint': info(path), 'statusAsStored': data.get('status'),
            'runDirectory': str(run), 'routes': chosen}


def source(path, names):
    text = path.read_text()
    tree = ast.parse(text)
    found = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name in names:
            segment = ast.get_source_segment(text, node)
            found.append({'name': node.name, 'line': node.lineno,
                          'endLine': node.end_lineno,
                          'sourceSha256': hashlib.sha256(segment.encode()).hexdigest()})
    if {v['name'] for v in found} != set(names):
        raise ValueError(('source names not found', str(path), names))
    return {**info(path), 'definitions': found}


def transcript(stage, select):
    stem = ('S11c_c1_bulk_closure_sympy_audit' if stage == 'c1'
            else 'S11c_c2_selfenergy_fold_sympy_audit')
    path = LEDGER / 'scripts/out' / (stem + '.out')
    manifest_path = (PROJECT / '_scratch/s11c/s11c-thickness-coordinate-20260914'
                     / (stage + '_full') / 'manifest.json')
    manifest = json.loads(manifest_path.read_text())
    expected = manifest['artifacts']['full.out']
    selected = []
    h = hashlib.sha256()
    offset = 0
    # Stream printed records. Do not eval srepr, import the codec or build CAS values.
    with path.open('rb') as stream:
        for number, line in enumerate(stream, 1):
            h.update(line)
            tag = line[:256].partition(b': ')[0].decode('ascii', errors='replace')
            if select(tag):
                labels = ('DRIVEN_BOUNDARY_SYSTEM', 'OUTGOING_PHI_FLAT',
                          'OUTGOING_PHI_SCATTERED', 'CONTROL_SURFACE_FLUX',
                          'FARFIELD_LIMIT', 'DELTA_P_SOURCE', 'IDENTIFICATIONS',
                          'REFERENCE_TRACE_MAP', 'REFERENCE_PRESSURE', 'NORMAL_JET',
                          'DENSITY_BINDING', 'RESOLVENT_INVERSE_OPERAND')
                fields = {k: line.find(("Str('" + k + "')").encode()) for k in labels}
                selected.append({'tag': tag, 'line': number, 'byteOffset': offset,
                                 'bytesIncludingNewline': len(line),
                                 'sha256': hashlib.sha256(line).hexdigest(),
                                 'literalFieldOffsets': {k: v for k, v in fields.items() if v >= 0}})
            offset += len(line)
    if offset != expected['bytes'] or h.hexdigest() != expected['sha256']:
        raise ValueError(('transcript does not match completed manifest', stage))
    if not selected or manifest['exit_code'] != 0:
        raise ValueError(('no completed saved transcript route', stage))
    script = LEDGER / 'scripts' / (stem + '.py')
    key = str(script.relative_to(LEDGER))
    return {'path': str(path), 'resolvedPath': str(path.resolve(strict=True)),
            'sha256': h.hexdigest(), 'bytes': offset, 'producerManifest': info(manifest_path),
            'currentProducerSourceMatches': digest(script) == manifest['source_hashes_after'][key],
            'records': selected, 'CASPayloadParsedOrRecomputed': False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    checkpoints = [
        checkpoint('uniform_source', lambda n: n in ('uniform-source.pickle', 'common.pickle',
                    'reference.pickle', 'left.pickle', 'right.pickle')),
        checkpoint('remaining_case_uniform', lambda n: n == 'remaining-case-uniform.pickle'
                    or n == 'input-routes.json'),
        checkpoint('profile_form', lambda n: n in ('binding-checks.pickle', 'profile-form.pickle')),
        checkpoint('remaining_case_profile_response', lambda n: n == 'remaining-case-profile-form.pickle'),
        checkpoint('coordinate_source', lambda n: n in ('density-advection.pickle', 'coordinate-source.pickle')),
        checkpoint('remaining_case_coordinate_inputs', lambda n: n.startswith('cases/') and
                   n.endswith(('density-operands.pickle', 'density-checks.json', 'density-family-pair.pickle'))),
        checkpoint('remaining_case_coordinate_sources', lambda n: n.startswith('coordinate-cases/')
                   and n.endswith('/coordinate-source.pickle')),
        checkpoint('remaining_case_first_jet_sources', lambda n: n.startswith('cases/')
                   and n.endswith('/first-jet-sources.pickle')),
        checkpoint('reduced_action_source', lambda n: n in ('reduced-action.pickle', 'actions.pickle')),
    ]
    specifications = {
        'scripts/S11c_d_mixing_scattering_sympy_audit.py': ['ReducedPencil', 'ReducedActionAssembly',
            'EdgeReconstruction', 'save_reduced_action_cache'],
        '_measurements/S11c_d_uniform_source.py': ['construct', 'emit_result'],
        '_measurements/S11c_d_profile_form.py': ['moments'],
        '_measurements/S11c_d_coordinate_source.py': ['first_jet_mutation', 'construct'],
        '_measurements/S11c_d_remaining_case_coordinate_inputs.py': ['density_operands'],
        'scripts/S11c_c1_bulk_closure_sympy_audit.py': ['outgoing_farfield_poynting', 'response_operator_case'],
        'scripts/S11c_c2_selfenergy_fold_sympy_audit.py': ['reference_pressure_kernels', 'build_face', 'run'],
    }
    result = {
        'status': 'A9_SAVED_SOURCE_ROUTE_INDEX_NOT_SCIENTIFIC_CLEARANCE',
        'method': 'Read named JSON manifests, source AST and literal output bytes only; no project/scientific import, pickle restoration, CAS evaluation, numerical call or previous worker replay.',
        'checkpoints': checkpoints,
        'sources': [source(LEDGER / name, definitions) for name, definitions in specifications.items()],
        'transcripts': [
            transcript('c1', lambda t: t in ('PY_S11CC1_ENERGY_FACE_TRACTION_OPERAND',
                       'PY_S11CC1_ENERGY_BULK_FARFIELD_FLUX_OPERAND', 'PY_S11CC1_ENERGY_RESIDUAL')),
            transcript('c2', lambda t: t.startswith('PY_S11CC2_FOLD_SYMBOL_MAP_'))],
        'scienceCalls': 0, 'readerSource': info(Path(__file__).resolve()),
        'remaining': 'Restore only selected required saved values under the guard; resolve their exact case/grade/domain coverage. Navigation alone closes no physics control.'}
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'output': str(args.output), 'checkpoints': len(checkpoints),
                      'selectedRoutes': sum(len(v['routes']) for v in checkpoints),
                      'transcriptRecords': sum(len(v['records']) for v in result['transcripts']),
                      'sourceMatches': [v['currentProducerSourceMatches'] for v in result['transcripts']],
                      'scienceCalls': 0}))


if __name__ == '__main__':
    main()
