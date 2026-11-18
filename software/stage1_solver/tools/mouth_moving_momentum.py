"""Moving-mouth momentum balance and a two-trace diagnostic on the resting drain."""

import numpy as np
import sympy as sp

from held_crossflow_budget import archived_packet
from mouth_feed_weak_check import dx4


def moving_balance():
    t, z = sp.symbols("t z", real=True)
    position = sp.Function("Z")(t)
    delta_minus, delta_plus, speed = sp.symbols("p_minus p_plus Wn")
    flux_minus, flux_plus = sp.symbols("Pi_minus Pi_plus")
    step = sp.Heaviside(z - position)
    density = delta_minus + (delta_plus - delta_minus)*step
    flux = flux_minus + (flux_plus - flux_minus)*step
    singular = sp.expand(sp.diff(density, t) + sp.diff(flux, z))
    jump = singular.coeff(sp.DiracDelta(z - position)).subs(position.diff(t), speed)
    incoming = flux_minus - speed*delta_minus
    outgoing = flux_plus - speed*delta_plus
    print("bulk jump in the moving balance =", jump)
    if sp.expand(jump - outgoing + incoming) != 0:
        raise ArithmeticError("Moving normal transport failed")

    # For a carrier inventory with side edges, subtract its moving-frame
    # momentum m*V*N. Use number conservation for the same region.
    m, v, a, number = sp.symbols("m V a N")
    rates = sp.symbols("J_in J_out J_edge")
    momentum_rates = sp.symbols("T_in T_out T_edge")
    force = sp.Symbol("f")
    def net(readings):
        left, right, edge = readings
        return left - right - edge
    lab = force + net(momentum_rates)
    frame_rates = [p - m*v*j for p, j in zip(momentum_rates, rates)]
    remaining = sp.expand(lab - m*(a*number + v*net(rates)) - force - net(frame_rates))
    print("coordinate acceleration term in carrier-frame momentum rate =", remaining)
    if remaining != -m*number*a:
        raise ArithmeticError("Carrier frame or edge transport failed")


def trace_difference():
    x, w = sp.symbols("x w", real=True)
    beta, m, eps, wn = sp.symbols("beta m eps Wn", real=True)
    density, slope, shape, theta, normal_speed = [
        sp.Function(name)(x) for name in ("n", "b", "h", "theta", "vN")
    ]
    phase = theta + normal_speed*w  # velocity potential in this calculation
    def trace(b):
        n = density*sp.exp(b*w)
        j = n*(sp.diff(phase, w) - wn)
        stress = -beta*n*sp.diff(sp.log(n), x, w)
        return sp.Matrix([j, m*j*sp.diff(phase, x), stress]).subs(w, 0)
    difference = (trace(slope) - trace(slope + eps*shape)).applyfunc(sp.simplify)
    print("feed minus receiver [number, carriage, mixed stress] =", list(difference))
    if difference != sp.Matrix([0, 0, beta*density*eps*shape.diff(x)]):
        raise ArithmeticError("Moving trace subtraction failed")


def saved_aperture():
    path = archived_packet("mouth_feed_weak_run.json", "n32_restoring_slope_u0.npz")
    with np.load(path, allow_pickle=False) as data:
        face, j, mouth = data["field"][-1], data["current"], data["mouth"]
        dx = float(data["dx"])
        if float(data["speed"]) != 0 or face.shape != (32, 32, 32):
            raise ValueError("Expected the saved resting drain")
    coord = (np.arange(32) + 0.5 - 16)*dx
    x, y, z = np.meshgrid(coord, coord, coord, indexing="ij", sparse=True)
    r2 = x*x + y*y + z*z
    if not np.array_equal(r2 < 1, mouth) or np.any(j[~mouth]) or not np.all(j[mouth] > 0):
        raise ValueError("Unexpected aperture or current support")
    n = np.abs(face)**2
    vx = (face.conj()*dx4(face, dx)).imag/n
    b = -2*np.maximum(1-r2, 0)**3*(1-1/np.sqrt(n))
    # A local normal log-density slope variation, supported well inside
    # the mouth even after the four-point stencil. All other traces agree.
    h = x*np.maximum(1-4*r2, 0)**3
    def integral(field):
        return dx**3*np.sum(field[mouth])
    def total(slope, frame_speed):
        return integral(j*(vx-frame_speed) - n*dx4(slope, dx)/4)
    rate = integral(j)
    print(f"saved inward number rate = {rate:.16g}")
    if np.any(dx4(h, dx)[~mouth]) or not np.isclose(rate, 0.25, rtol=0, atol=1e-13):
        raise ArithmeticError("Aperture support or rate failed")
    expected = integral(n*dx4(h, dx))/4
    for eps in (0, 0.01, -0.01, 0.005):
        frame_speed = 0.02
        lab = total(b, 0) - total(b + eps*h, 0)
        boosted = total(b, frame_speed) - total(b + eps*h, frame_speed)
        one_side_shift = total(b, frame_speed) - total(b, 0)
        print(f"slope amplitude={eps:g} face difference={lab:.15g} "
              f"boosted={boosted:.15g} single-face frame shift={one_side_shift:.15g}")
        np.testing.assert_allclose([lab, boosted, one_side_shift],
                                   [eps*expected, eps*expected, -rate*frame_speed],
                                   rtol=1e-8, atol=2e-17)
    print(f"face difference per unit trial slope amplitude = {expected:.15g}", flush=True)
    if not np.isfinite(expected) or abs(expected) < 1e-13:
        raise ArithmeticError("Trial slope did not distinguish the two traces")


if __name__ == "__main__":
    moving_balance()
    trace_difference()
    saved_aperture()
