#!/usr/bin/env python3
"""Check numerical parameter gaps against the accepted native end operators."""
import ast
import hashlib
import json
from pathlib import Path
import pickle
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import S11c_d_mixing_scattering_sympy_audit as engine


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    inventory_path = ROOT/'_measurements/S11c_d_variable_profile_parameter_inventory.json'
    inventory = json.loads(inventory_path.read_text())
    names = set(inventory['unboundGradientEnergyCoefficients'])
    source_path = ROOT/'_measurements/S11c_d_reduced_action_source_checkpoint.json'
    source = json.loads(source_path.read_text())
    base = ROOT.parents[1]/'_scratch/s11c/s11c-thickness-coordinate-20260914/d_full'
    manifest_path = base/'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    if digest(manifest_path) != source['sourceJoins']['producerManifestSha256']:
        raise ValueError('accepted native end producer identity')
    engine_name = str(engine.HERE.relative_to(ROOT))
    frozen_path = base/'source'/engine_name
    if digest(frozen_path) != manifest['source_hashes_after'][engine_name]:
        raise ValueError('native producer source snapshot changed')
    frozen, current = ast.parse(frozen_path.read_text()), ast.parse(engine.HERE.read_text())
    constructors = {}
    for name in ('EdgeReduction', 'ReducedPencil', 'ConstantEndPencil', 'dag_substitute', 'polynomial_terms'):
        a = next(n for n in frozen.body if getattr(n, 'name', None) == name)
        b = next(n for n in current.body if getattr(n, 'name', None) == name)
        constructors[name] = ast.dump(a) == ast.dump(b)
    matrices = {}
    for end in ('REFERENCE', 'LEFT', 'RIGHT'):
        name = 'symbols/'+end+'_LAB_HELD_RHO4_CONSTANT.pickle'
        path = base/name
        if digest(path) != manifest['artifacts'][name]['sha256']:
            raise ValueError(('native end cache changed', end))
        with path.open('rb') as stream:
            packet = pickle.load(stream)
        matrices[end] = {'cacheSha256': digest(path), 'matrices': {}}
        for label, value in (('WEAK', packet[0]), ('STRONG', packet[4])):
            symbols = engine.dag_free_symbols(value)
            present = sorted(names & {s.name for s in symbols})
            matrices[end]['matrices'][label] = {'shape': list(value.shape),
                'presentInventoriedCoefficients': present, 'presentCount': len(present)}
    record = {'case': 'LAB_HELD_RHO4_CONSTANT', 'coefficientCount': len(names),
        'inventorySha256': digest(inventory_path), 'sourceCheckpointSha256': digest(source_path),
        'nativeManifestSha256': digest(manifest_path), 'nativeConstructorAstJoins': constructors,
        'endOperatorDependencies': matrices,
        'scope': 'Exact symbolic dependency census of the accepted native strong and weak end matrices. The proposed numerical extension must preserve every existing parameter and profile.'}
    print(json.dumps(record, indent=2), flush=True)
    if not all(constructors.values()) or any(m['presentCount'] for e in matrices.values() for m in e['matrices'].values()):
        raise ValueError('end-parameter dependency requires inspection')


if __name__ == '__main__':
    main()
