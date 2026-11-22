"""Necessary energy for one freely moving site to reach a downstream slice.

Use the archived stationary mouth and the same trial smooth radial Morse well
as the supplied-flow experiment. Give the site independent inertia in w. The
analytic minimum is over the whole receiving slice, so the check allows a
site to leave the well, take another path, or carry sideways momentum.
"""
import argparse

import numpy as np
from scipy.optimize import minimize_scalar

from mouth_calibrated_flare_coating import read_inputs
from mouth_middle_tag_grip import cubic_trace
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import BETWEEN, emit, morse


def run(mesh):
    trace, intake = cubic_trace(mesh)
    patch, _, _ = read_inputs(mesh)
    scale = 2 * patch['shared_scale']
    depth, _, rest = BETWEEN
    entry = np.array([1.175])
    chi_in = trace.sample(trace.chi, entry)[0].item()
    speed_w = trace.speed(entry).item()
    slope = trace.geometry.radius(entry)[1].item()
    speed_r = slope * speed_w
    kinetic = 0.5 * (speed_w**2 + speed_r**2)
    entry_potential = -scale * depth * chi_in
    energy = kinetic + entry_potential
    emit(test='preparation', mesh=mesh, intake=intake, entry=entry.item(),
         chi_in=chi_in, radial_to_axial_speed=slope, speed_w=speed_w,
         speed_r=speed_r, kinetic=kinetic, potential=entry_potential,
         energy=energy, zero_work_min_chi=-energy/(scale * depth),
         units='specific energy in squared mouth radii per squared model time')
    assert 0 < chi_in <= 1 and abs(intake - .25) < 1e-8
    for label, w in [('entry', entry.item()), ('half_order', trace.crossing(.5)),
                     ('quarter_order', trace.crossing(.25)), ('gauge', 1.55), ('exit', 1.65)]:
        at = np.array([w])
        chi = trace.sample(trace.chi, at)[0].item()
        radius = trace.geometry.radius(at)[0].item()
        floor = -scale * depth * chi

        def potential(r):
            return scale * chi * morse(rest + (r - radius)/PITCH, *BETWEEN)[0]

        minimized = minimize_scalar(potential, bounds=(radius - .1*PITCH, radius + PITCH),
                                    method='bounded', options={'xatol': 1e-14})
        if not minimized.success or abs(minimized.fun - floor) > 3e-13:
            raise ArithmeticError('Slice-minimum check failed')
        rise = floor - entry_potential
        emit(test='slice', label=label, depth=w, chi=chi, radius=radius,
             analytic_minimum=floor, numerical_minimum=minimized.fun,
             potential_rise=rise, energy_shortfall=floor-energy,
             minimum_added_work=max(0., floor-energy),
             critical_entry_speed=(2*max(0., rise))**.5,
             # Anisotropy check at the same radial and axial launch speeds.
             minimum_axial_to_radial_mass_ratio=max(0., (2*rise-speed_r**2)/speed_w**2))
    level = -energy / (scale * depth)
    if 0 < level < chi_in:
        boundary = trace.crossing(level)
        if abs(trace.sample(trace.chi, np.array([boundary]))[0].item()-level) > 1e-11:
            raise ArithmeticError('Inconsistent energy boundary')
        emit(test='necessary_reach_boundary', depth=boundary, chi=level,
             scope='First forbidden lower-order slices beyond this crossing at zero external work. Actual sites can turn earlier.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=int, required=True, choices=(26,34))
    run(parser.parse_args().mesh)
