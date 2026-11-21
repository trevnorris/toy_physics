"""Force and material-tag audit on a cold, flat, partially softened neutral patch."""
import json
from itertools import combinations

import numpy as np
from scipy.integrate import solve_ivp

from neutral_patch_conversion import NEUTRAL
from neutral_patch_resting_order_hold import hessian, rest as equilibrate
from neutral_patch_sites import COUNT, Patch, emit
from neutral_patch_turning_neighbor import AREA_PER_SITE


ENTRY = 0.95
DRIVES = (0., .001, -.001, .0005, -.0005)
TIMES = np.array([0., .25, .5, .75, 1., 1.25, 1.5, 2.])


def pulse(t):
    return np.sin(np.pi*t)**2 if t < 1. else 0.


def upper_tags(patch, energy, distances, reference_depth):
    # Both rules use the same mutual nearest neighbours. Only the energy
    # denominator changes. The upper plane can move in either tangent direction.
    current_edges = patch.neighbors(energy, distances, .6) & patch.cross
    reference_edges = current_edges & (-energy > .6*reference_depth)
    return {name: np.isin(np.arange(COUNT, 2*COUNT), patch.j[edges])
            for name, edges in [('current_well', current_edges), ('entry_well', reference_edges)]}


def trial(patch, rest, reference_depth, amplitude):
    def observe(t, upper):
        points = rest.copy()
        points[:COUNT, 0] += amplitude*pulse(t)
        points[COUNT:] = upper.reshape(COUNT, 2)
        return points, patch.evaluate(points)

    def drift(t, upper):
        # Native unit-mobility overdamped sites at zero temperature.
        _, (_, gradient, _) = observe(t, upper)
        return -gradient[COUNT:].ravel()

    # End the driven integration exactly when the pulse ends. The free
    # segment tests the relaxation of these very same site identities.
    runs = []
    y = rest[COUNT:].ravel()
    for begin, end in [(0., 1.), (1., 2.)]:
        ts = TIMES[(TIMES > begin) & (TIMES <= end)]
        sol = solve_ivp(drift, (begin, end), y, t_eval=ts, method='DOP853',
                        rtol=2e-11, atol=2e-13)
        if not sol.success:
            raise ArithmeticError(sol.message)
        runs.extend(sol.y.T)
        y = sol.y[:, -1]
    states = [rest[COUNT:].ravel(), *runs]
    energy0, _, d0 = patch.evaluate(rest)
    cohorts = upper_tags(patch, energy0, d0, reference_depth)
    output = []
    for t, state in zip(TIMES, states, strict=True):
        points, (energy, gradient, distances) = observe(t, state)
        tags = upper_tags(patch, energy, distances, reference_depth)
        speed = -gradient[COUNT:, 0]
        force = float(sum(speed))  # unit mobility, summed material force
        np.testing.assert_allclose(force, gradient[:COUNT, 0].sum(), rtol=0, atol=2e-11)
        # Separate work readout on the same snapshot, using only cross bonds.
        e = []
        for h in (-.00001, .00001):
            shifted = points.copy()
            shifted[COUNT:, 0] += h
            e.append(patch.evaluate(shifted)[0][patch.cross].sum())
        work_force = (e[0]-e[1])/.00002
        np.testing.assert_allclose(work_force, force, rtol=1e-5, atol=2e-9)
        grouped = {}
        for name, mask in cohorts.items():
            grouped[name] = dict(count=int(mask.sum()), changed=int(sum(tags[name] != mask)),
                tagged_speed=float(speed[mask].mean()) if mask.any() else None,
                other_speed=float(speed[~mask].mean()) if (~mask).any() else None)
        output.append(dict(t=float(t), material_force=force, energy_check=work_force,
                           upper_mean_speed=float(speed.mean()), upper_site_speeds=speed.tolist(),
                           groups=grouped))
    return output


def main():
    frame = json.loads(NEUTRAL.read_text())['saved_measurements']['frames1321.json']['contents']['central']
    initial = np.array(frame['rested_xy']).ravel()
    # A clamped drive face, all coordinates on the receiving face free.
    mobile = np.eye(4*COUNT)[:, 2*COUNT:]
    area = COUNT * AREA_PER_SITE
    for chi in (ENTRY, .5, .25, .1):
        patch = Patch('central')
        baseline_depth = patch.depth.copy()
        patch.depth[patch.cross] *= chi/ENTRY
        equilibrium, residual = equilibrate(patch, initial, mobile)
        minimum = float(np.linalg.eigvalsh(mobile.T @ hessian(patch, equilibrium) @ mobile)[0])
        if minimum <= 0.:
            raise ArithmeticError('The stationary patch has an unstable unconstrained site mode')
        rest = equilibrium.reshape(2*COUNT, 2)
        all_runs = {amp: trial(patch, rest, baseline_depth, amp) for amp in DRIVES}
        snapshots = []
        for i, time in enumerate(TIMES):
            full, half = DRIVES[1], DRIVES[3]
            def odd(amp, name):
                return (all_runs[amp][i][name]-all_runs[-amp][i][name])/(2*amp)
            signed = {str(a): all_runs[a][i] for a in DRIVES}
            snapshot = dict(time=float(time), signed=signed,
                full_stress_per_peak_step=odd(full, 'material_force')/area,
                half_stress_per_peak_step=odd(half, 'material_force')/area)
            if chi == .25:
                # Same identities on both drives. Holding the uniform quarter
                # count fixes three labels per face, but does not specify which.
                v = (np.array(signed[str(full)]['upper_site_speeds']) -
                     np.array(signed[str(-full)]['upper_site_speeds']))/(2*full)
                contrast = []
                for ids in combinations(range(COUNT), COUNT//4):
                    ordered = np.isin(np.arange(COUNT), ids)
                    contrast.append(v[ordered].mean() - v[~ordered].mean())
                snapshot['count_only_relative_speed_per_peak_step'] = dict(
                    assignments=len(contrast), minimum=min(contrast), maximum=max(contrast))
            snapshots.append(snapshot)
        full = np.array([s['full_stress_per_peak_step'] for s in snapshots])
        half = np.array([s['half_stress_per_peak_step'] for s in snapshots])
        linear_error = float(np.linalg.norm(full-half)/np.linalg.norm(half))
        if linear_error >= 1e-3:
            raise ArithmeticError(f'Half-step material response differs by {linear_error}')
        emit(test='neutral_shear_tags', chi=chi, equilibrium_residual=residual, min_upper_stiffness=minimum,
             sites_per_plane=COUNT, nominal_planar_area=area, half_step_relative_error=linear_error,
             max_abs_zero_force=max(abs(r['material_force']) for r in all_runs[0.]),
             upper_rest_cohorts=all_runs[0.][0]['groups'], snapshots=snapshots)


if __name__=='__main__':
    main()
