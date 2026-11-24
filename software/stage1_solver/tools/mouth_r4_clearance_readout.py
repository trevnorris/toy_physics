"""Reduce the two smooth-map quadratures and independently differenced reservoir work."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_cp10_spatial_crystal import MODES
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read


def verify_guard(folder):
    outcome = read(folder / 'outcome.json')
    limits = read(folder / 'limit-validation.json')
    if (outcome['exitCode'] != 0 or outcome['childOutcome']['exitCode'] != 0
            or not outcome['limitsVerified'] or not limits['verified']
            or int(limits['actual']['memory.swap.max']) != 0):
        raise ValueError(f'Invalid guarded measurement at {folder}')


def as_list(x):
    return np.asarray(x, dtype=float).tolist()


def refinement_orders(magnitudes):
    with np.errstate(divide='ignore',invalid='ignore'):
        rates = np.log2(magnitudes[:-1]/magnitudes[1:])
    return np.where(np.isfinite(rates), rates.astype(object), None).tolist()


def reduce_degree(shapes, depth, bath):
    if (shapes['checkpoint'] != CHECKPOINT_SHA or depth['checkpoint_sha256'] != CHECKPOINT_SHA
            or bath['checkpoint_sha256'] != CHECKPOINT_SHA
            or shapes['common_quadrature'] != depth['partition_sha256']
            or shapes['mask_sha256'] != depth['mask_sha256']
            or shapes['quadrature_degree'] != depth['degree']
            or shapes['force_modes'] != MODES[1:]):
        raise ValueError('Different control preparations or quadratures')
    if abs(depth['common_bath_depth_force'] - bath['analytic']) > 1e-12:
        raise ArithmeticError('Common reservoir force changed')
    levels = []
    for level in ('1', '2', '4'):
        s = shapes['results'][level]
        d = depth['rows'][level]
        # Both integrations independently give the same stored energy and explicit
        # order work. The latter uses the physical order knots, including their motion.
        region = np.array(s['synthetic_region_delta'])
        axes = np.array(d['error_by_region'])
        cross = np.max(np.abs(region[:, :2] - axes[:, [0, 2]]))
        if cross > 2e-12:
            raise ArithmeticError(f'Common-point partition mismatch {cross}')
        combined = np.c_[region[:, 0], axes[:, 1], region[:, 2:], axes[:, 2:]]
        analytic = np.c_[d['full_force'], s['analytic_total']]
        raw = np.r_[s['delta'][0], np.diff(d['full_force'])[0], s['delta'][2:]]
        if np.max(abs(raw - combined[:, :7].sum(axis=0))) > 2e-12:
            raise ArithmeticError('Regional force or energy missing from closure')
        other_force_error = analytic[:, None, 1:] - np.array(s['central_total'])
        depth_errors = np.column_stack([
            np.array(v['analytic_minus_fd_total']) + bath['rows'][str(v['h'])]['analytic_minus_fd']
            for v in d['fd']
        ])
        # Signed depth probes move piecewise-smooth knots. Keep every step visible
        # and use the smaller central difference for the six-force comparison.
        fd = np.c_[depth_errors[:, -1], other_force_error[:, -1, :]]
        rich = np.c_[depth_errors[:, -1], s['total_analytic_minus_richardson']]
        levels.append(dict(level=int(level), totals=as_list(combined.sum(axis=0)),
                           by_region=as_list(combined), cross_integral_difference=float(cross),
                           crystal_energy_banks=[v[0] for v in s['material_spatial']],
                           analytic_total_banks=as_list(analytic),
                           analytic_minus_fd_total_by_step=dict(depth=as_list(depth_errors),
                               other_five=as_list(other_force_error)),
                           analytic_minus_smaller_central=as_list(fd),
                           analytic_minus_depth_central_and_five_richardson=as_list(rich),
                           depth_one_sided_max_spread=float(np.max(np.asarray(d['directional_spread'])[:,:,1])),
                           common_points_in_shape_integral=s['points']))
    values = np.array([v['totals'] for v in levels])
    magnitudes = abs(values)
    by_region = abs(np.array([v['by_region'] for v in levels]))
    return dict(levels=levels,
                magnitude_decreases_both_refinements=np.all(np.diff(magnitudes,axis=0)<0,axis=0).tolist(),
                signs_consistent_both_refinements=np.all(np.sign(values[1:])==np.sign(values[:-1]),axis=0).tolist(),
                observed_orders=refinement_orders(magnitudes),
                observed_orders_by_region=refinement_orders(by_region))


def run(scratch, destination):
    area = scratch/'r4_clearance_smooth01'
    verify_guard(area/'bath_guard')
    bath_file = area/'bath.json'
    bath = read(bath_file)
    files = [bath_file]
    regions = Path(__file__).resolve().parents[3]/'research/pde_ledger_v3/_measurements/mouth_r4_clearance_regions.json'
    record = dict(checkpoint_sha256=CHECKPOINT_SHA, baseline_commit='344030b4',
                  difference='sampled spatial Q1 inverse minus exact smooth inverse, same physical quadrature',
                  columns=['stored_crystal_energy', *MODES, 'order_work', 'remaining_depth_work'],
                  regions=['smooth interior A', 'weld B', 'saved melting front and fold neighborhood C'],
                  region_record_sha256=digest(regions),
                  banks=['analytic inverse map', 'sampled spatial inverse map'],
                  signed_depth_steps=[4e-7,2e-7], signed_other_steps=[1e-4,5e-5],
                  bath_depth=bath, degrees={})
    for n in (6,8):
        verify_guard(area/f'shape_guard_{n}')
        verify_guard(area/f'depth_guard_{n}')
        shape_path = area/f'shape_{n}/result.json'
        depth_path = area/f'depth/{n}/result.json'
        files.extend([shape_path, depth_path])
        shape, depth = read(shape_path), read(depth_path)
        if shape['mask_sha256'] != read(regions)['masks_sha256']:
            raise ValueError('Unexpected weld or melt-front partition')
        record['degrees'][n] = reduce_degree(shape, depth, bath)
    record['inputs'] = {str(p.relative_to(scratch)): digest(p) for p in files}
    sources = ['mouth_r4_clearance_shapes.cpp','mouth_r4_clearance_shapes.py',
               'mouth_r4_smooth_points.hpp','mouth_r4_smooth_depth.cpp','mouth_r4_smooth_depth.py',
               'mouth_r4_clearance_bath.py',Path(__file__).name,'mouth_cp10_spatial_crystal.cpp',
               'mouth_free_shape_solid.cpp','mouth_cp10_field_shape.py']
    record['sources'] = {s: digest(Path(__file__).with_name(s)) for s in sources}
    with destination.open('x') as f:
        json.dump(record,f,indent=2,allow_nan=False)
        f.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('scratch',type=Path)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    run(args.scratch.resolve(),args.output)
