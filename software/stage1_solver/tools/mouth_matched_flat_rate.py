"""Compare stationary curved and restrained planar conversion under matched baths.

The normal coordinate points from the ordered brane to receiving bulk. Both
models use the same neutral phases, equation of state and passive mobility.
A planar translation restraint supplies the external momentum difference.
Its release gives an initial material crossing rate, not a later steady state.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson
from scipy.optimize import brentq

from mouth_free_front_wave_head import bath_head, emit
from mouth_held_order_current import measure


RECORD = (Path(__file__).resolve().parents[3] /
          'research/pde_ledger_v3/_measurements/mouth_coupled_conversion.json')
MOBILITY = 0.5
EXTENT = 25.0


def curved_controls():
    reference = json.loads(RECORD.read_text())['runs']['curved_pressure']
    sha = hashlib.sha256(reference['literal_stdout'].encode()).hexdigest()
    if (sha != reference['stdout_sha256'] or not reference['limits_verified']
            or reference['guard_exit'] != 0 or reference['worker_exit'] != 0):
        raise ValueError('Unverified curved reference')
    rows = [json.loads(line) for line in reference['literal_stdout'].splitlines()]
    controls = {row['trial_action_per_inner_area']: row for row in rows
                if row['test'] == 'free_spherical_front'
                and row['brane_pressure_excess'] == .0045}
    emit(test='reference', reference_sha256=sha, equilibrium_well_difference=0.,
         kappa=2., mobility=MOBILITY, density_brane=1.)
    for action in (0., .001):
        row = controls[action]
        j = row['current_per_inner_area']
        emit(test='adopted_curve', action=action, pressure=row['brane_pressure_excess'],
             current_per_front_area=row['current_at_front'],
             reservoir_energy_per_constituent=row['reservoir_power']/j,
             heat_per_constituent=row['conversion_loss']/j,
             reservoir_power_per_front_area=row['current_at_front']*row['reservoir_power']/j)
    return controls


def flat(pressure_drop, width=1., spacing=.125, tolerance=1e-9):
    return measure(-pressure_drop, MOBILITY, width, EXTENT, spacing, tolerance, return_fields=True)


def report_face(tag, pressure, fields):
    row, sol, fixture = fields
    j, strength = sol.p
    x = np.linspace(-EXTENT, EXTENT, 8001)
    f, df, q, dq, B, _ = sol.sol(x)
    n = f**2
    ddq = sol.sol(x, 1)[3]
    mu_free = 2*q*(1-q)*(1-2*q)-2*ddq
    integrated_loss = simpson(MOBILITY*(mu_free-strength*fixture(x))**2, x=x)
    half = brentq(lambda s: sol.sol(s)[2]-.5, -EXTENT, EXTENT)
    amp, _, chi, slope, _, _ = sol.sol(half)
    curvature = sol.sol(half, 1)[3]
    u = j/amp**2
    mu_free_half = 2*chi*(1-chi)*(1-2*chi)-2*curvature
    velocity_after_release = u + MOBILITY*mu_free_half/slope
    velocity_from_holder = MOBILITY*strength*fixture(half)/slope
    momentum_residual = row['normal_flux_difference']-row['fixture_on_medium']
    emit(test='supported_face', ensemble=tag, brane_pressure_excess=pressure,
         mobility=MOBILITY, width=row['fixture_width'], spacing=row['mesh_step'],
         rms_residual=row['rms_residual_max'], nodes=row['nodes'],
         bulk_density=n[-1], current_per_front_area=j, head=B[0]-B[-1],
         reservoir_power_per_front_area=row['bath_power'],
         conversion_power_per_front_area=integrated_loss,
         energy_per_constituent_from_conversion=integrated_loss/j,
         conservative_well_energy_per_constituent=0., holder_on_material=row['fixture_on_medium'],
         momentum_flux_jump=row['normal_flux_difference'],
         pressure_jump_no_support_if_j_relative_zero=n[-1]**5/5-n[0]**5/5,
         min_order=min(q), max_order=max(q), maximum_order_rise=max(np.diff(q)),
         end_order_stress_max=max(abs((dq*dq-(q*(1-q))**2)[[0,-1]])),
         initial_level_position=half, half_order_carrier_speed=u,
         initial_level_speed_on_release=velocity_after_release,
         release_speed_from_removed_restraint=velocity_from_holder,
         initial_relative_current_on_release=j-amp**2*velocity_after_release)
    if abs(momentum_residual)>1e-9 or abs(row['bath_power']-integrated_loss)>2e-11:
        raise ArithmeticError('Momentum or sustained power balance failed')
    if abs(velocity_after_release-velocity_from_holder)>2.e-9:
        raise ArithmeticError('Two half-order release speed evaluations disagree')
    if min(n)<=0 or np.max(np.diff(q))>1e-9:
        raise ArithmeticError('Invalid order or density')
    return j, B[0]-B[-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=('pressure', 'head', 'compact-holder'), required=True)
    args = parser.parse_args()
    controls = curved_controls()
    curved = controls[.001]
    pressure = curved['brane_pressure_excess']
    if args.case == 'compact-holder':
        report_face('common_static_pressure_compact_holder', pressure, flat(pressure, .8, .080, 2e-10))
        return
    if args.case == 'pressure':
        j_face, head_face = report_face('common_static_pressure', pressure, flat(pressure))
        for action in (0., .001):
            reference = controls[action]
            emit(test='comparison', ensemble='common_static_pressure', action=action,
                 j_face_over_front=j_face/reference['current_at_front'],
                 face_to_curved_energy_per_constituent=head_face/reference['head_brane_minus_bulk'])
        return

    # Recompute the flat pressure that reproduces the full specific-head drop
    # of the occupied spherical drain. The curved bath velocities differ from
    # the planar velocities even when the far static pressures match.
    occupied_head = curved['head_brane_minus_bulk']
    samples = {}
    def difference(trial):
        fields = flat(trial)
        samples[trial] = fields
        n_bulk = (1-5*trial)**.2
        j = fields[0]['current_to_bulk']
        return bath_head(1.,j)-bath_head(n_bulk,j)-occupied_head

    pressure_for_head = brentq(difference, .00455, .00475, xtol=2e-14)
    j_matched, solved_head = report_face('common_specific_head_drop', pressure_for_head,
                                         samples[pressure_for_head])
    emit(test='comparison', ensemble='common_specific_head_drop', action=.001,
         face_pressure_excess=pressure_for_head, curved_pressure_excess=pressure,
         j_face_over_front=j_matched/curved['current_at_front'],
         solved_head=solved_head, reference_head=occupied_head)
    if abs(solved_head-occupied_head)>1e-11:
        raise ArithmeticError('Matched head condition failed')


if __name__=='__main__':
    main()
