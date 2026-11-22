"""Compare the independent triangular work with the main mapped-grid stationary result."""
import argparse
import hashlib
import json
from pathlib import Path

CASES = ((1, 'broad'), (2, 'broad'), (2, 'mid'), (2, 'tight'))
SOURCES = [
    'software/stage1_solver/tools/mouth_live_waist_checks/shape_work.py',
    'software/stage1_solver/tools/mouth_live_waist_checks/closed_wave.py',
    'software/stage1_solver/tools/mouth_lip_checks/shape_check.py',
    'software/stage1_solver/tools/mouth_lip_checks/duct.py',
    'software/stage1_solver/tools/mouth_exchange_geometry.py',
    'software/stage1_solver/tools/mouth_conversion_shape_mode.py',
    'research/pde_ledger_v3/_measurements/mouth_conversion_shape.json',
    'research/pde_ledger_v3/_measurements/mouth_sealed_shape_relax.json',
]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_guarded(directory):
    guard = json.loads((directory / 'outcome.json').read_text())
    assert guard['limitsVerified'] and guard['exitCode'] == guard['childOutcome']['exitCode'] == 0
    assert guard['childOutcome']['guardReason'] is None
    capture = directory / 'stdout'
    return [json.loads(line) for line in capture.read_text().splitlines()], digest(capture)


def signed_slope(rows, direction, step, energy):
    energies = {r['offset']: r[energy] for r in rows
                if r['kind'] == 'signed' and r['coordinate'] == direction}
    return (energies[step] - energies[-step]) / (2 * step)


def compare(primary_dir, independent_dir, exterior_dir):
    result = {}
    for radius, name in CASES:
        case = f'r{radius}_{name}_72'
        primary, primary_hash = read_guarded(primary_dir / case)
        independent, independent_hash = read_guarded(independent_dir / f'fine{radius}{name}')
        # Two original exterior jobs did not start a child; use the admitted outcome.
        outer_case = 'first_r2mid70' if radius == 2 and name == 'mid' else f'r{radius}{name}70'
        exterior, exterior_hash = read_guarded(exterior_dir / outer_case)
        basis = independent[0]
        assert basis['data_sha256'] == digest(basis['source'])
        assert basis['radius'] == radius and basis['state'] == name
        data = {
            'input_order_sha256': basis['data_sha256'],
            'output_sha256': {'triangular': independent_hash, 'mapped': primary_hash, 'outer': exterior_hash},
            'closed_frequencies': basis['independent_clamped'],
            'nearest_closed_index': basis['nearest'],
            'open_frequency': basis['primary_pole'],
        }
        for work in independent[1:]:
            direction = work['direction']
            coord = {'waist': 'waist', 'axial': 'depth'}[direction]
            outer_axis = {'waist': 'radial', 'axial': 'axial'}[direction]
            outer_slope = next(r['slope'] for r in exterior
                               if r['axis'] == outer_axis and r.get('h') == .001)
            mapped = signed_slope(primary, coord, .001, 'total_inside_and_drives')
            mapped_double_step = signed_slope(primary, coord, .002, 'total_inside_and_drives')
            parts = {key: work[key] for key in (
                'fixed_action_closed_wave', 'resolved_core_order_at_transported_chi', 'linear_sealed_carrier')}
            independent_sum = sum(parts.values())
            assert min(mapped, mapped_double_step, independent_sum, outer_slope) > 0
            row = {
                'signed_step': .002,
                'independent_interior_terms': parts,
                'independent_interior_sum': independent_sum,
                'mapped_interior_0.001': mapped,
                'mapped_interior_0.002': mapped_double_step,
                'relative_interior_difference': independent_sum / mapped - 1,
                'mapped_exterior_order_0.001_n70': outer_slope,
                'mapped_total_0.001_exterior_n70': mapped + outer_slope,
            }
            if direction == 'axial':
                row['closed_wave_exact_axial'] = work['axial_wave_identity']
                row['closed_wave_axial_fd_minus_exact'] = parts['fixed_action_closed_wave'] - work['axial_wave_identity']
            data[direction] = row
        result[case] = data
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--primary', type=Path, required=True)
    parser.add_argument('--independent', type=Path, required=True)
    parser.add_argument('--exterior', type=Path, required=True)
    args = parser.parse_args()
    record = {
        'method': 'Independent physical P1 triangles with third-order Gaussian product quadrature. Closed shear cavity on saved live order, weak sealed-carrier control and scalar interior virtual work. Mapped-grid comparison uses reoptimized order, density, outgoing shear and resolved carrier. The exterior column is reduced from the main tagged outer half-front run at mesh 70.',
        'coordinates': 'waist means b(t) + 4t(1-t)q at fixed endpoints and length. Axial means log(depth). Positive energy derivatives drive a negative generalized force.',
        'limitations': 'This independently checks force signs at shared input states. It does not solve for those states independently. The nearest of three closed eigenfrequencies is a pole approximation. The deep R2 branch is the second closed pole. The wave and carrier models differ between the two discretizations, so a small shift in magnitudes is expected. A registered shear boundary at half order does not measure radiation through an unmodeled exterior shear tail.',
        'inputs': {str(p): digest(p) for p in SOURCES},
        'cases': compare(args.primary, args.independent, args.exterior),
    }
    print(json.dumps(record, indent=2))
