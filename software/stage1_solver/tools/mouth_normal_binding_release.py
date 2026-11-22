"""Release radial and fourth-coordinate motion for two matched sideways-force laws."""
import argparse

import numpy as np
from scipy.integrate import solve_ivp

from mouth_calibrated_flare_coating import read_inputs
from mouth_middle_tag_grip import cubic_trace
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import BETWEEN, emit


def run(mesh, tolerance):
    trace, intake = cubic_trace(mesh)
    patch, _, _ = read_inputs(mesh)
    entry, downstream = 1.175, 1.65
    strength = 2*patch['shared_scale']*BETWEEN[0]
    alpha = BETWEEN[1]/PITCH
    wz = np.array([entry])
    radius, slope = trace.geometry.radius(wz)
    speed = trace.speed(wz).item()
    initial = np.array([radius.item(), entry, slope.item()*speed, speed])
    transit = trace.transit(entry, downstream, 12)
    for compensation in (0., 1.):
        def field(r, w):
            radius, slope = trace.geometry.radius(w)
            order, derivative = trace.sample(trace.chi, w)
            exp = np.exp(-alpha*(r-radius))
            shape = exp*(exp-2) + compensation
            radial = 2*alpha*exp*(1-exp)
            potential = strength*order*shape
            return potential, -strength*order*radial, strength*(order*radial*slope-derivative*shape)

        def rates(t, y):
            _, fr, fw = field(y[0:1], y[1:2])
            return [y[2], y[3], fr.item(), fw.item()]

        def turn(t, y):
            return y[3]

        def exit_slice(t, y):
            return y[1]-downstream

        turn.terminal = exit_slice.terminal = True
        turn.direction, exit_slice.direction = -1, 1
        e0 = field(initial[0:1], initial[1:2])[0].item() + .5*sum(initial[2:]**2)
        state = initial
        start = drift = 0.
        nodes = 0
        while True:
            sol = solve_ivp(rates, (start, start+transit), state, method='DOP853',
                            events=(turn, exit_slice), rtol=tolerance, atol=tolerance*1e-2,
                            max_step=.012)
            if not sol.success:
                raise ValueError(sol.message)
            energies = field(sol.y[0], sol.y[1])[0] + .5*(sol.y[2]**2+sol.y[3]**2)
            drift = max(drift, np.max(abs(energies-e0)).item())
            nodes += len(sol.t)-1
            state, start = sol.y[:, -1], sol.t[-1]
            if sol.status == 1:
                break
        event = 'first_axial_turn' if len(sol.t_events[0]) else 'crossed_downstream_slice'
        emit(test='release', law='fixed_continuum_reference' if not compensation else 'level_well_bottom',
             scale=strength, intake=intake, mesh=mesh, rtol=tolerance, initial_state=initial.tolist(),
             initial_energy=e0, event=event, radial_position=state[0], depth=state[1],
             radial_speed=state[2], axial_speed=state[3], at_time=start, accepted_steps=nodes,
             maximum_energy_drift=drift)
        if drift > 2e-7:
            raise ArithmeticError('Conservative trajectory failed its energy check')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh', type=int, required=True, choices=(26,34))
    p.add_argument('--rtol', type=float, required=True)
    args = p.parse_args()
    run(args.mesh, args.rtol)
