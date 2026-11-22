"""Compare the independent strip minima and first-order susceptibility with the primary readouts."""
import argparse
import json
from pathlib import Path


def read_json(path):
    return json.loads(path.read_text())


def compare(linear, finite, primary):
    threshold = primary['comparison']['critical_YH_over_sigma']
    linear_values = [
        ('weighted_second_moment', linear['weighted_moment'], primary['comparison']['flat_I2']),
        ('neutral_tangent', linear['B_neutral_native'], primary['native']['tangent']),
        ('unit_modulus_tangent', linear['B_neutral_native'] + linear['weighted_moment'],
         primary['registered_trial']['tangent']),
        ('critical_modulus', linear['threshold_Y'], threshold),
    ]
    for name, independent, reported in linear_values:
        yield dict(check=name, independent=independent, primary=reported,
                   error=independent-reported, tolerance=1.e-9)
    for modulus, branch in [(0., 'native'), (1., 'registered_trial')]:
        fine, coarse = finite[20, modulus], finite[10, modulus]
        fields = primary[branch]
        for row in fields['training'] + fields['held_out']:
            key = str(row['curvature'])
            m20 = fine['energy_difference_moment'][key]
            m10 = coarse['energy_difference_moment'][key]
            extrapolated = (4.*m20 - m10) / 3.
            expected = row['holding_moment']
            yield dict(check='signed_curvature_moment', modulus=modulus, curvature=row['curvature'],
                       independent=m20, extrapolated=extrapolated, primary=expected,
                       error=extrapolated-expected, tolerance=2.e-10)
        outer = fine['energy_difference_moment']['0.08']
        held = fine['energy_difference_moment']['0.04']
        observed_ratio = abs(outer / (2.*held) - 1.)
        reported_ratio = abs(fields['predictions'][1]/fields['held_out'][1]['holding_moment'] - 1.)
        yield dict(check='reserved_curvature_relative_prediction_error', modulus=modulus,
                   independent=observed_ratio, primary=reported_ratio,
                   error=observed_ratio-reported_ratio, tolerance=1.e-5)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ['flat_sensitivity_stdout', 'signed_moments_stdout', 'primary_json']:
        parser.add_argument(argument, type=Path)
    args = parser.parse_args()
    rows = [json.loads(line) for line in args.signed_moments_stdout.read_text().splitlines()]
    finite = {(row['elements_per_width'], row['Y']): row for row in rows}
    results = list(compare(read_json(args.flat_sensitivity_stdout), finite,
                           read_json(args.primary_json)))
    for finding in results:
        print(json.dumps(finding), flush=True)
    if any(abs(row['error']) > row['tolerance'] for row in results):
        raise ArithmeticError('Independent strip response disagrees with a primary readout')
