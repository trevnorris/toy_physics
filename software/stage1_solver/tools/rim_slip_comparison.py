"""Compare normal rim slip for a held drain and a translating drain at 0PN."""

import numpy as np
import sympy as sp

from held_crossflow_budget import archived_packet
from mouth_feed_weak_check import check_file, dx4


def rim_projection():
    a, d, distance, source = sp.symbols("a d l Q", positive=True)
    k, tau, u, v = sp.symbols("k tau U V", real=True)
    m = sp.symbols("m", real=True)
    gamma = 2 * sp.pi * a**2 * k * tau * sp.integrate(m**2, (m, -1, 1))
    force = gamma * (u - v)
    gamma_hold = sp.diff(force, u)
    gamma_drift = -sp.diff(force, v)
    print("held coefficient, opposing-drift coefficient =", gamma_hold, gamma_drift, flush=True)
    if gamma_hold != gamma_drift or gamma_hold != 4 * sp.pi * a**2 * k * tau / 3:
        raise ArithmeticError("Normal-slip projection failed")

    # Companion at +d on the force axis. l grows from d-a to d+a while
    # m = cos(theta) decreases from +1 to -1. The projected material velocity
    # here is the undisturbed potential flow, a conditional test input.
    mu = (d**2 + a**2 - distance**2) / (2 * a * d)
    radial_velocity = source * (d * mu - a) / (4 * sp.pi * distance**3)
    integrand = sp.Rational(3, 2) * mu * radial_velocity * distance / (a * d)
    primitive = sp.integrate(sp.expand(integrand), distance)
    dipole = sp.factor(primitive.subs(distance, d + a) - primitive.subs(distance, d - a))
    uniform = source / (4 * sp.pi * d**2)
    print("point-drain velocity dipole, uniform reference, difference =",
          dipole, uniform, sp.factor(dipole - uniform), flush=True)
    if sp.factor(dipole - uniform) != 0:
        raise ArithmeticError("Companion dipole mismatch")


def carrier_frames():
    for speed in (0.01, 0.02):
        path = archived_packet("mouth_intake_redistribution.json", f"n32_restoring_slope_u{speed:g}.npz")
        held = check_file(path, "restoring_slope")
        with np.load(path, allow_pickle=False) as packet:
            face = packet["field"][-1]
            current = packet["current"]
            dx = float(packet["dx"])
            if float(packet["speed"]) != speed:
                raise ValueError("Incorrect feed phase")
        density = abs(face)**2
        lab_velocity = (face.conj() * dx4(face, dx)).imag / density
        rate = float(dx**3 * current.sum())
        moving_carry = dx**3 * np.sum(current * lab_velocity)
        # A uniform sideways boost leaves the density stress unchanged.
        held_flux = held[3] + held[4]
        moving_flux = moving_carry + held[4]
        boost_flux = speed * rate
        print(f"U={speed:g} V={-speed:g} J={rate:.12g} held_lab={held_flux:+.12g} "
              f"moving_lab={moving_flux:+.12g} difference={held_flux-moving_flux:+.12g} "
              f"U*J={boost_flux:+.12g}", flush=True)
        np.testing.assert_allclose([rate, held_flux - moving_flux], [0.25, boost_flux], rtol=0, atol=2e-15)


if __name__ == "__main__":
    rim_projection()
    carrier_frames()
