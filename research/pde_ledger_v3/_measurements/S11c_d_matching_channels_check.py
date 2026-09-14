#!/usr/bin/env python3
"""Construct and validate both-end matching channels from accepted sources."""
import argparse
import ast
import contextlib
import json
from pathlib import Path
import pickle
import resource
import shutil
import time
from types import SimpleNamespace

import numpy as np
import sympy as sp

import S11c_d_end_normalization_check as normalization
from S11c_d_end_normalization_check import ROOT, engine, digest, decoded_lines, _restore
from S11c_d_output_codec import restore_emission_index

PREFIX = 'MATCHING_CHANNELS_LAB_HELD_RHO4_CONSTANT'
SOURCES = (*normalization.SOURCES, '_measurements/S11c_d_matching_channels_check.py',
           '_measurements/S11c_d_variable_profile_matching_plan.md')


def save(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def source_context(end):
    path = ROOT/'_measurements'/('S11c_d_end_normalization_'+end.lower()+'_thickness_repair_checkpoint.json')
    checkpoint = json.loads(path.read_text())
    if checkpoint['end'] != end or checkpoint['unaccountedResidualNormsAboveDiagnosticThreshold']:
        raise ValueError('accepted matching-end normalization required')
    base = Path(checkpoint['runDirectory'])
    for name, item in checkpoint['artifacts'].items():
        if digest(base/name) != item['sha256']:
            raise ValueError(('normalization artifact changed', end, name))
    for name, item in checkpoint['publications'].items():
        if digest(ROOT/item['path']) != item['sha256']:
            raise ValueError(('normalization publication changed', end, name))
    for name, expected in checkpoint['sourceFiles'].items():
        if digest(base/'source'/name) != expected:
            raise ValueError(('normalization source snapshot changed', end, name))
        if name not in ('scripts/S11c_d_mixing_scattering_sympy_audit.py',
                        '_measurements/S11c_d_end_normalization_validate.py') and digest(ROOT/name) != expected:
            raise ValueError(('normalization helper changed', end, name))
    frozen = ast.parse((base/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py').read_text())
    current = ast.parse(Path(engine.__file__).read_text())
    current.body = [n for n in current.body if getattr(n, 'name', None) != 'TwoEndedMatchingChannels']
    if ast.dump(frozen) != ast.dump(current):
        raise ValueError('existing native engine changed beyond the new matching class')
    arguments = json.loads((base/'arguments.json').read_text())
    args = SimpleNamespace(**{k: Path(v) if v is not None and k != 'end' else v for k, v in arguments.items()})
    pairing, current, inputs, provenance = normalization.load_pairing(args)
    native, coverage, bindings, physical = normalization.load_native(args, pairing, inputs, provenance)
    modal, known = pickle.loads((base/'modal.pickle').read_bytes())
    adjoint, adjoint_known = pickle.loads((base/'adjoint.pickle').read_bytes())
    if modal['NATIVE_RECORDS'] != native or modal['NATIVE_COVERAGE'] != coverage:
        raise ValueError('complete native root/lift join')
    for a, b, c in zip(modal['RECORDS'], adjoint['RECORDS'], native):
        for key in ('ROOT_DISK_INDEX', 'NORMAL_LIFT_SIGN', 'NULLITY'):
            if a[key] != b[key] or a[key] != int(c[key]):
                raise ValueError(('modal/adjoint/native join', end, key))
    if len(modal['RECORDS']) != len(adjoint['RECORDS']) or len(modal['RECORDS']) != len(native):
        raise ValueError('complete candidate census')
    engine.PHYSICAL_METADATA.dimensions.known.update(known)
    engine.PHYSICAL_METADATA.dimensions.known.update(adjoint_known)
    builder = engine.ModalCurrentSubspaces(pairing, current, bindings)
    source = {'normalizationCheckpoint': str(path.relative_to(ROOT)),
        'normalizationCheckpointSha256': digest(path), 'inputSha256': provenance['inputSha256'],
        'producerManifestSha256': provenance['producerManifestSha256'], 'end': end,
        'modalSha256': digest(base/'modal.pickle'), 'adjointSha256': digest(base/'adjoint.pickle'),
        'unitFrame': inputs.frame, 'profileSpecification': inputs.specification['profiles'],
        'physicalPencilJoin': provenance['physicalPencilJoin']}
    return engine.TwoEndedMatchingChannels(builder, modal, adjoint, end), source, checkpoint


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-directory', type=Path, required=True)
    parser.add_argument('--publish', action='store_true')
    args = parser.parse_args()
    base = args.run_directory.resolve(); base.relative_to(ROOT.parents[1]/'_scratch/s11c')
    base.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    pins = {name: digest(ROOT/name) for name in SOURCES}
    for name in SOURCES:
        target = base/'source'/name; target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/name, target)
    contexts, results, sources, summaries = {}, {}, {}, {}
    with (base/'full.out').open('w') as transcript, contextlib.redirect_stdout(transcript):
        for end in ('LEFT', 'RIGHT'):
            context, source, checkpoint = source_context(end)
            contexts[end], sources[end] = context, source
            result = context.construct(); results[end] = result
            normalization.atomic(base/(end.lower()+'.pickle'), pickle.dumps(result, protocol=5))
            context.emit(result, source, PREFIX+'_'+end)
            actual = {(r['RECORD_INDEX'], r['BASIS_COLUMN'], r['DIRECTION']) for r in result['CHANNELS']}
            expected = {(r['INDEX'], column, direction) for r in checkpoint['orientationRecords']
                        if r['OUTWARD_FLUX_CLASSIFICATION_DEFINED']
                        for direction, key in (('INCOMING', 'INCOMING_BASIS_COLUMNS'), ('OUTGOING', 'OUTGOING_BASIS_COLUMNS'))
                        for column in r[key]}
            summaries[end] = {'counts': result['COUNTS'], 'orientationSetResidual': len(actual ^ expected),
                'residualNorms': {k: float(np.linalg.norm(v)) for k, v in result['RESIDUALS'].items()},
                'interRootCurrentNorm': float(np.linalg.norm(result['INTER_ROOT_CURRENT'])),
                'candidateRankResiduals': {key: sum(abs(r['INFO'][key]) for r in result['CANDIDATES']) for key in
                    ('RIGHT_RANK_MINUS_NULLITY', 'LEFT_RANK_MINUS_NULLITY', 'NATIVE_NULLITY_RESIDUAL')}}
            normalization.put(PREFIX+'_'+end, 'ORIENTATION_SET_RESIDUAL', len(actual ^ expected), context.builder.modes)
        names = [tag.removeprefix('PY_S11CD_') for tag in engine.EMISSION_LINES if not tag.startswith('PY_S11CD_METADATA_')]
        keys = {name: 's11cd'+''.join(v.title() for v in name.split('_')) for name in names}
        modes = contexts['RIGHT'].builder.modes
        normalization.put(PREFIX, 'WRITE_KEYS', keys, modes)
        index = engine.emission_index(engine.EMISSION_LINES)
        normalization.put(PREFIX, 'EMISSION_LINES', index, modes)
    # Replay the computed emitter after closing the transcript. This checks
    # every serialized object, fingerprint, unit and grade against its source.
    entries = {}
    for line in decoded_lines(base/'full.out'):
        tag, separator, payload = line.rstrip('\n').partition(': ')
        if not separator or tag in entries:
            raise ValueError('matching transcript framing/census')
        entries[tag] = _restore(payload)
    seen = set(); original = engine.emit
    def check_emit(name, value):
        tag = 'PY_S11CD_'+name
        if tag in seen or entries.get(tag) != engine.cas(value):
            raise ValueError(('matching computed/emitted join', tag))
        seen.add(tag)
    engine.emit = check_emit
    try:
        for end in ('LEFT', 'RIGHT'):
            contexts[end].emit(results[end], sources[end], PREFIX+'_'+end)
            normalization.put(PREFIX+'_'+end, 'ORIENTATION_SET_RESIDUAL', summaries[end]['orientationSetResidual'], modes)
        normalization.put(PREFIX, 'WRITE_KEYS', keys, modes)
        normalization.put(PREFIX, 'EMISSION_LINES', index, modes)
    finally:
        engine.emit = original
    if seen != set(entries) or len(keys) != len(set(keys.values())) or set(keys.values()) & set(engine.IMPORT_KEYS):
        raise ValueError('matching emission or write-key census')
    final = 'PY_S11CD_'+PREFIX+'_EMISSION_LINES'
    indexed = restore_emission_index({str(k): v for k, v in entries[final]}, list(entries)[:list(entries).index(final)])
    if set(indexed) != set(list(entries)[:list(entries).index(final)]):
        raise ValueError('matching emission index census')
    metadata_paths = 0
    for tag, body in entries.items():
        if tag.startswith('PY_S11CD_METADATA_'):
            for row in body:
                d = {str(k): v for k, v in row}
                if len(d['DIMENSION_L_T_M']) != 3 or any(v.free_symbols for v in d['DIMENSION_L_T_M']):
                    raise ValueError(('unresolved matching dimensions', tag))
                if not all(k in d for k in ('PATHS', 'MULTIGRADE', 'EPSILON_LAMBDA_SUPPORT')):
                    raise ValueError('missing matching grades')
                metadata_paths += len(d['PATHS'])
    record = {'runDirectory': str(base), 'sourceFiles': pins, 'sources': sources, 'ends': summaries,
        'tagCount': len(entries), 'metadataPaths': metadata_paths, 'writeKeyCount': len(keys),
        'wallSeconds': time.monotonic()-started, 'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'artifacts': {name: {'bytes': (base/name).stat().st_size, 'sha256': digest(base/name)}
                      for name in ('full.out', 'left.pickle', 'right.pickle')},
        'scope': 'One supplied case: both-end full candidate bases and cross-mode boundary currents. No interior solve, S-matrix, continuum re-expansion or global spectral coverage.'}
    save(base/'checks.json', record)
    if pins != {name: digest(ROOT/name) for name in SOURCES}:
        raise ValueError('matching sources changed during execution')
    if sources['LEFT']['inputSha256'] != sources['RIGHT']['inputSha256']:
        raise ValueError('different two-end profile/parameter input')
    if any(row['orientationSetResidual'] or any(row['candidateRankResiduals'].values()) or
           any(value > 1e-8 for value in row['residualNorms'].values()) for row in summaries.values()):
        raise ValueError('matching channel residual requires inspection; operands preserved')
    if args.publish:
        target = ROOT/'scripts/out/S11c_d_matching_channels.out'
        checkpoint = ROOT/'_measurements/S11c_d_matching_channels_checkpoint.json'
        if target.exists() or target.is_symlink() or checkpoint.exists():
            raise ValueError('matching publication already exists')
        temporary = target.with_name('.'+target.name+'.new')
        shutil.copyfile(base/'full.out', temporary)
        if digest(temporary) != record['artifacts']['full.out']['sha256']:
            raise ValueError('matching publication hash')
        temporary.replace(target)
        record['publication'] = {'path': str(target.relative_to(ROOT)), **record['artifacts']['full.out']}
        save(checkpoint, record)
    print(json.dumps({k: record[k] for k in ('ends', 'tagCount', 'metadataPaths', 'wallSeconds')}, indent=2))


if __name__ == '__main__':
    main()
