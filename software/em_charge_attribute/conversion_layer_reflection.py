"""Exact taper solution and scattering readout for two transverse inertia laws."""

import json
from math import pi

import numpy as np
import sympy as sp
from scipy.optimize import brentq
from scipy.special import jv


def check_taper(exponent):
    # (x u')' + H^2 x^m u = 0, from bulk edge x=0 to brane edge x=1.
    x, H = sp.symbols("x H", positive=True)
    p = sp.Rational(exponent + 1, 2)
    u = sp.besselj(0, H*x**p/p)
    equation = sp.simplify(sp.expand_func(sp.diff(x*sp.diff(u, x), x)
                                         + H**2*x**exponent*u))
    v = sp.symbols("v")
    phase = sp.series(sp.atan(sp.besselj(1, v).series(v, 0, 5).removeO()
                              / sp.besselj(0, v).series(v, 0, 5).removeO()), v, 0, 5)
    print(f"m={exponent} order-equation residual={equation}, half reflection phase={phase}",
          flush=True)
    if equation != 0:
        raise ArithmeticError("Taper solution fails the wave equation")


def layer_edge(k, width, exponent, loss=0.0):
    # Units c=mu0=rho0=1, zero common gap. Time convention exp(-i omega t).
    inner_k = np.sqrt(k*k + 1j*loss*k) if loss else k
    v = 2*inner_k*width/(exponent+1)
    # Bare displacement and traction / (mu0 k) on the ordered side of the ramp.
    return jv(0, v), (inner_k/k)*jv(1, v)


def rim_reflection(k, center, width, exponent, loss=0.0):
    u, t = layer_edge(k, width, exponent, loss)
    return np.exp(2j*k*(center-width/2))*(u+1j*t)/(u-1j*t)


def first_modes(center, width, exponent):
    def node(k):
        u, t = layer_edge(k, width, exponent)
        angle = k*(center-width/2)
        return u*np.cos(angle) - t*np.sin(angle)

    scan = np.linspace(1e-6, 8, 161)
    modes = [brentq(node, a, b, xtol=5e-15) for a, b in zip(scan, scan[1:])
             if node(a)*node(b) < 0][:2]
    if len(modes) != 2:
        raise ArithmeticError("Missing resonance")
    return np.array(modes)


def readings(exponent):
    base_center, base_width = 1., .10
    probes = first_modes(base_center, base_width, exponent)
    def inferred_length(roots):
        return np.sqrt(2)*pi/np.sqrt(roots[1]**2-roots[0]**2)
    base_length = inferred_length(probes)
    for name, center, width in (("base", 1., .10), ("shift", 1.05, .10), ("widen", 1., .20)):
        modes = first_modes(center, width, exponent) if name != "base" else probes
        roots = rim_reflection(modes, center, width, exponent)
        for loss in (0., .1):
            r = rim_reflection(probes, center, width, exponent, loss)
            ref = rim_reflection(probes, base_center, base_width, exponent, loss)
            phase = np.angle(r/ref)  # continuous branch for these small changes
            result = {"m": exponent, "geometry": name, "loss": loss, "power": list(abs(r)**2),
                      "relative_phase": list(phase), "phase_shift": list(phase/(2*probes))}
            if not loss:
                result.update(modes=list(modes), spacing_length=inferred_length(modes),
                              spacing_shift=inferred_length(modes)-base_length)
            print(json.dumps(result, allow_nan=False), flush=True)
            if np.max(abs(roots+1)) > 1e-11:
                raise ArithmeticError("Resonance does not satisfy the rim node")
            if (np.max(abs(r)) > 1+1e-12 or np.max(abs(phase)) >= pi/2
                    or (not loss and np.max(abs(abs(r)-1)) > 1e-12)):
                raise ArithmeticError("Nonpassive, unwrapped or nonconservative result")
            if name == "shift" and np.max(abs(phase/(2*probes)-(center-base_center))) > 1e-12:
                raise ArithmeticError("Translation phase check failed")


if __name__ == "__main__":
    for power_exponent in (0, 1):
        check_taper(power_exponent)
        readings(power_exponent)
