#!/usr/bin/env python3
"""Join the approved development input to the unchanged endpoint input."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MEASUREMENTS = ROOT/'_measurements'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check():
    paths = {name: MEASUREMENTS/('S11c_d_'+suffix+'.json') for name, suffix in {
        'proposal': 'variable_profile_development_input_proposal',
        'input': 'variable_profile_development_input',
        'base': 'channel_preflight_input',
        'dependencies': 'variable_profile_parameter_checkpoint',
        'inventory': 'variable_profile_parameter_inventory',
        'matching': 'matching_channels_checkpoint'}.items()}
    values = {n: json.loads(p.read_text()) for n, p in paths.items()}
    proposal, selected, base = (values[n] for n in ('proposal', 'input', 'base'))
    names = sorted(values['inventory']['unboundGradientEnergyCoefficients'])
    if selected != proposal['input'] or digest(paths['base']) != proposal['baseInputSha256']:
        raise ValueError('approved proposal/base identity changed')
    if selected['profiles'] != base['profiles'] or selected['unit_frame'] != base['unit_frame']:
        raise ValueError('profile or reference-unit frame changed')
    differences = {k: str(Fraction(selected['parameters'][k])-Fraction(v))
                   for k, v in base['parameters'].items()}
    added = sorted(set(selected['parameters'])-set(base['parameters']))
    if any(Fraction(v) for v in differences.values()) or added != names:
        raise ValueError('input extension changed an existing parameter or inventory')
    if any(Fraction(selected['parameters'][n]) != Fraction(i, 101)
           for i, n in enumerate(names, 1)):
        raise ValueError('selected coefficients differ from the approved rule')
    dependencies = values['dependencies']
    if dependencies['inventorySha256'] != digest(paths['inventory']):
        raise ValueError('end-dependency inventory changed')
    if not all(dependencies['nativeConstructorAstJoins'].values()):
        raise ValueError('end-dependency constructor join')
    if any(m['presentCount'] or m['presentInventoriedCoefficients']
           for end in dependencies['endOperatorDependencies'].values()
           for m in end['matrices'].values()):
        raise ValueError('extended coefficients occur in accepted end pencils')
    endpoint_joins = {end: record['inputSha256'] == digest(paths['base'])
                      for end, record in values['matching']['sources'].items()}
    if not all(endpoint_joins.values()):
        raise ValueError('matching endpoint input identity')
    return {'status': 'USER_APPROVED_INPUT_JOINED',
        'authorization': 'Go ahead with the proposed values',
        'files': {n: {'path': str(p.relative_to(ROOT)), 'sha256': digest(p)} for n, p in paths.items()},
        'existingParameterDifferences': differences, 'addedCoefficientCount': len(added),
        'addedCoefficientUnits': proposal['addedCoefficientUnits'],
        'matchingBaseInputJoins': endpoint_joins,
        'inheritedEndOperatorDependencies': dependencies['endOperatorDependencies'],
        'scope': 'The selected extension preserves all old parameters, profiles and units. '
                 'All added coefficients are absent from the accepted strong/weak end pencils. '
                 'Existing end records retain their original input hashes; this projection joins their reuse.'}


if __name__ == '__main__':
    print(json.dumps(check(), indent=2))
