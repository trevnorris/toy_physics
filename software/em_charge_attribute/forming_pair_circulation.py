"""Orbital angular momentum of a neutral scalar wave in a prescribed splitting well."""

import json
import numpy as np


def separation(time, forming):
    tau = time / 4
    return (2 * tau**2 * (3 - 2*tau), 3*tau*(1-tau)) if forming else (0., 0.)


def half_reading(u, v, potential, x, w, spacing, side):
    ux, uw = np.gradient(u, spacing, edge_order=2)
    energy = (v*v + ux*ux + uw*uw + potential*u*u)/2
    px, pw = -v*ux, -v*uw
    weights = np.where(side*w > 0, 1., np.where(w == 0, .5, 0.)) * spacing**2
    def integral(field):
        return float(np.sum(weights * field))
    total = integral(energy)
    xc, wc = integral(x*energy)/total, integral(w*energy)/total
    spin = integral((x-xc)*pw - (w-wc)*px)
    edge = integral(energy*((abs(x) >= 8) | (abs(w) >= 8)))/total
    band = abs(abs(w) - 2)
    wells = integral(energy*np.where(band < 1, 1., np.where(band == 1, .5, 0.)))/total
    return dict(L=spin, E=total, x=xc, w=wc, edge=edge, well_band=wells)


def run(spacing, forming):
    axis = np.arange(-9., 9.+spacing/2, spacing)
    x, w = axis[:, None], axis[None, :]
    gauss = np.exp(-(x*x+w*w)/2)
    # Independent evolution of the finite background and the unit velocity seed.
    # Their superpositions solve both handed trials and the two controls exactly.
    u = np.stack([x*gauss, np.zeros_like(gauss)])
    v = np.stack([np.zeros_like(gauss), (w*w-.5)*gauss])
    for field in (u, v):
        field[:, [0, -1], :] = 0
        field[:, :, [0, -1]] = 0
    dt = spacing/8
    duration = 4
    steps = round(duration/dt)
    def potential(time):
        d, speed = separation(time, forming)
        return x*x + (w*w-d*d)**2/4, -d*speed*(w*w-d*d)
    def acceleration(field, pot):
        result = np.zeros_like(field)
        center = field[:, 1:-1, 1:-1]
        result[:, 1:-1, 1:-1] = (
            (field[:, 2:, 1:-1] + field[:, :-2, 1:-1] + field[:, 1:-1, 2:] +
             field[:, 1:-1, :-2] - 4*center)/spacing**2 - pot[1:-1, 1:-1]*center)
        return result
    def drive_power(slope):
        return float(np.sum(slope*(u[0]+.01*u[1])**2)*spacing**2/2)
    pot, slope = potential(0)
    acc = acceleration(u, pot)
    work = 0.
    previous_power = drive_power(slope)
    initial_energy = None
    for step in range(steps + 1):
        time = step*dt
        if step in (0, steps//2, steps):
            for seed in (.01, -.01, .005, 0.):
                wave, velocity = u[0]+seed*u[1], v[0]+seed*v[1]
                upper, lower = [half_reading(wave, velocity, pot, x, w, spacing, s)
                                for s in (1, -1)]
                result = dict(h=spacing, forming=forming, t=time, seed=seed,
                              up=upper, down=lower)
                if seed == .01:
                    total = upper["E"] + lower["E"]
                    if step == 0:
                        initial_energy = total
                    result["energy_work_residual_over_initial"] = (total-initial_energy-work)/initial_energy
                print(json.dumps(result, allow_nan=False), flush=True)
        if step == steps:
            break
        v += dt*acc/2
        u += dt*v
        pot, slope = potential(time+dt)
        acc = acceleration(u, pot)
        v += dt*acc/2
        power = drive_power(slope)
        work += dt*(previous_power+power)/2
        previous_power = power


if __name__ == "__main__":
    for h in (1/8, 1/16):
        for split in (False, True):
            run(h, split)
