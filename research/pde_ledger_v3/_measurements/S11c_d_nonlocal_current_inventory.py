#!/usr/bin/env python3
"""Inventory the recorded current construction and its unresolved input join."""
import argparse
import hashlib
import json
from pathlib import Path
import pickle
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import S11c_d_mixing_scattering_sympy_audit as engine
from ledger_fold import _restore
from S11c_d_output_codec import decoded_lines


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    base = Path(manifest['run_directory'])
    for name, record in manifest['artifacts'].items():
        if sha(base/name) != record['sha256']:
            raise ValueError(('artifact digest mismatch', name))
    if manifest['source_hashes_before'] != manifest['source_hashes_after']:
        raise ValueError('source changed during run')
    b, c, a = pickle.loads((base/'objects.pickle').read_bytes())
    tags, gaps, nonfinite, zero_maps = set(), [], [], []
    dimension_constraints, unknown_dimensions = None, set()
    for line in decoded_lines(base/'full.out'):
        tag, _, body = line.partition(': ')
        if tag in tags:
            raise ValueError(('duplicate tag', tag))
        tags.add(tag)
        value = _restore(body)
        if value.has(sp.nan, sp.zoo, sp.oo, -sp.oo):
            nonfinite.append(tag)
        if '_METADATA_' in tag and "Str('ZERO_MAP')" in body:
            zero_maps.append(tag)
        if '_METADATA_' in tag:
            unknown_dimensions.update(str(s) for s in value.free_symbols if str(s).startswith('s11cdUnit'))
        if tag == 'PY_S11CD_NONLOCAL_CURRENT_PREFLIGHT_DIMENSION_CONSTRAINTS':
            dimension_constraints = sp.srepr(value)
    for tag in tags:
        if '_METADATA_' not in tag and tag.replace('PY_S11CD_', 'PY_S11CD_METADATA_', 1) not in tags:
            gaps.append(tag)
    residuals = {}
    for prefix, group in (('slab', b), ('acoustic', a)):
        for key, value in group.items():
            if key.endswith('_RESIDUAL'):
                values = [v for _, v in engine.leaves(engine.cas(value))]
                residuals[prefix+'_'+key] = {'scalarCount':len(values), 'nonzeroCount':sum(v != 0 for v in values)}
    residuals['faceEquations'] = {'scalarCount':sum(key.endswith('_RESIDUAL') for face in a['FACE_RECORDS'] for key in face),
                                'nonzeroCount':sum(v != 0 for face in a['FACE_RECORDS'] for key, v in face.items() if key.endswith('_RESIDUAL'))}
    # The k_W coefficient anchors the mechanical row orientation without
    # comparing or repairing the inherited shear sector.
    command = manifest['command']
    producer = json.loads(Path(command[command.index('--manifest')+1]).read_text())
    source_cache = Path(producer['run_directory'])/'symbols/REFERENCE_LAB_HELD_RHO4_CONSTANT.pickle'
    _, _, _, known, strong, _ = pickle.loads(source_cache.read_bytes())
    if sha(source_cache) != producer['artifacts']['symbols/REFERENCE_LAB_HELD_RHO4_CONSTANT.pickle']['sha256']:
        raise ValueError('source cache digest mismatch')
    symbols = {s.name:s for s in strong.free_symbols}
    energy_row = -b['CONSTRAINED_MECHANICAL_EULER_DERIVATIVES'][7]/b['HARMONIC_VARIATION_NORMALIZATION']
    eps = next(s for s in energy_row.free_symbols if s.name == 'epsilon_shape')
    energy_row = sp.expand(energy_row).coeff(eps, 2)
    kw = symbols['k_W']
    field = next(f for f in energy_row.atoms(sp.Function) if f.func.__name__ == 's11cdCurrentPlusE')
    source_anchor = sp.diff(strong[4, 4], kw)
    energy_anchor = sp.cancel(sp.diff(energy_row, kw)/field)
    anchor_residual = sp.cancel(source_anchor-energy_anchor)
    mass = a['REDUCED_MASS_FACE_JOIN_RESIDUAL']
    mechanical_sum = a['REDUCED_MECHANICAL_FACE_SUM']
    data = {'instrumentSha256':sha(Path(__file__)), 'producerManifestSha256':sha(args.manifest),
            'tagCount':len(tags), 'metadataGaps':gaps, 'nonfiniteTags':nonfinite, 'zeroMapMetadataTags':zero_maps,
            'dimensionConstraints':dimension_constraints, 'unknownDimensionSymbols':sorted(unknown_dimensions),
            'residuals':residuals,
            'mechanicalOrientationAnchor':{'source':sp.srepr(source_anchor), 'energy':sp.srepr(energy_anchor),
                'residual':sp.srepr(anchor_residual), 'dimensionLTM':[str(2*d) for d in known[symbols['W_0']]],
                'multigrade':[[0, 0, 0]], 'epsilonLambdaSupport':[[0, 0]]},
            'massJoin':sp.srepr(mass), 'mechanicalSum':sp.srepr(mechanical_sum),
            'diagonalDepthIntegral':sp.srepr(a['DIAGONAL_DEPTH_INTEGRAL']),
            'boundaryModulusLimit':sp.srepr(a['DEPTH_BOUNDARY_MODULUS_LIMIT'])}
    args.output.write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps(data, indent=2))


if __name__ == '__main__':
    run()
