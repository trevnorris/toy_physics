#!/usr/bin/env python3
"""Compute the class-P bulk pullback, reference-face trace, and Fourier pairing.

Inputs are the supplied rest-frame acoustic equations and the class-P ansatz.
Every printed expression is derived from those symbolic inputs.
"""

import sympy as sp

x1, x2, x3, w, t = sp.symbols("x1 x2 x3 w t", real=True)
k2, omega, rho_m, c_s0 = sp.symbols(
    "k2 omega rho_m c_s0", nonzero=True, real=True
)
s = sp.symbols("s", nonzero=True, real=True)
W0, eta = sp.symbols("W0 eta", real=True)
Phi = sp.Function("Phi")(x1, w)
w1 = sp.Function("w1")(x1)
phase = sp.exp(sp.I * (k2 * x2 - omega * t))
phi = Phi * phase
coordinates = (x1, x2, x3, w)


def amplitude(expression):
    return sp.simplify(sp.expand(expression / phase))


velocity = [sp.diff(phi, coordinate) for coordinate in coordinates]
pressure = -rho_m * sp.diff(phi, t)
wave_equation = sp.diff(phi, t, 2) + c_s0**2 * sum(
    sp.diff(phi, coordinate, 2) for coordinate in coordinates
)
phi_ww = sp.solve(
    sp.Eq(sp.simplify(wave_equation / phase), 0), sp.diff(Phi, w, 2)
)[0]

velocity_amp = [amplitude(component) for component in velocity]
pressure_amp = amplitude(pressure)
velocity_normal_jet = [
    sp.simplify(amplitude(sp.diff(component, w)).subs(sp.diff(Phi, w, 2), phi_ww))
    for component in velocity
]
pressure_normal_jet = amplitude(sp.diff(pressure, w))

delta_rho_amplitude = sp.Function("delta_rho")(x1, w)
delta_rho = delta_rho_amplitude * phase
continuity = sp.diff(delta_rho, t) + rho_m * sum(
    sp.diff(velocity[index], coordinates[index]) for index in range(4)
)
density_amp = sp.solve(
    sp.Eq(sp.simplify(continuity / phase), 0), delta_rho_amplitude
)[0]
density_amp = sp.simplify(density_amp.subs(sp.diff(Phi, w, 2), phi_ww))
density_normal_jet = sp.simplify(sp.diff(density_amp, w))
density_time_amp = sp.simplify(-sp.I * omega * density_amp)
current_amp = [sp.simplify(rho_m * component) for component in velocity_amp]
current_normal_jet = [
    sp.simplify(rho_m * component) for component in velocity_normal_jet
]

print("PRESSURE_TRACE", pressure_amp)
print("VELOCITY_TRACE", velocity_amp)
print("PRESSURE_NORMAL_JET", pressure_normal_jet)
print("VELOCITY_NORMAL_JET", velocity_normal_jet)
print("DENSITY_TRACE", density_amp)
print("DENSITY_MINUS_PRESSURE_OVER_C2", sp.simplify(density_amp - pressure_amp / c_s0**2))
print("DENSITY_NORMAL_JET", density_normal_jet)
print("DENSITY_TIME", density_time_amp)
print("CURRENT_TRACE", current_amp)
print("CURRENT_NORMAL_JET", current_normal_jet)

# The engine's affine perturbation is based at the flat reference face.
Phi_function = sp.Function("Phi")
w_ref = s * W0 / 2
w_bg = s * W0 * (1 + eta * w1) / 2
exact_at_background = Phi_function(x1, w_bg)
affine_from_reference = Phi_function(x1, w_ref) + (w_bg - w_ref) * sp.diff(
    Phi_function(x1, w), w
).subs(w, w_ref)
affine_from_background = Phi_function(x1, w_bg) + (w_bg - w_ref) * sp.diff(
    Phi_function(x1, w), w
).subs(w, w_bg)


def first_eta(expression):
    return sp.simplify(sp.series(expression, eta, 0, 2).removeO().doit())


print(
    "TRACE_EXACT_MINUS_REFERENCE_AFFINE",
    first_eta(exact_at_background - affine_from_reference),
)
print(
    "TRACE_EXACT_MINUS_BACKGROUND_AFFINE",
    first_eta(exact_at_background - affine_from_background),
)

# Spatial Fourier pairing: trial and test phases are conjugate.
trial_amplitude, test_amplitude = sp.symbols("trial_amplitude test_amplitude")
period = 2 * sp.pi / k2
trial = trial_amplitude * sp.exp(sp.I * k2 * x2)
test_same = test_amplitude * sp.exp(sp.I * k2 * x2)
test_conjugate = test_amplitude * sp.exp(-sp.I * k2 * x2)
print(
    "PAIRING_SAME_PHASE",
    sp.simplify(sp.integrate(test_same * sp.diff(trial, x2), (x2, 0, period))),
)
print(
    "PAIRING_CONJUGATE_PHASE",
    sp.simplify(
        sp.integrate(test_conjugate * sp.diff(trial, x2), (x2, 0, period))
    ),
)
