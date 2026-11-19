"""Cycle work by independently injected, equal-speed parcels at a moving loop.

Time dependence is Re[A exp(i omega t)]. Positive torque and work act on the
coordinate. In the shared-edge experiment either incident face samples one of
three cyclic placements of its marked contact; the face holonomy is conveyed
to that contact. All other faces keep their original held Hamiltonians.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

from mouth_conversion_shared_links import RATE, prepared
from mouth_conversion_twist import H0, RHO_UP, adj, contact, read, tr_product
from mouth_collar_order_preparation import calibration


DWELL = 1.
FREQUENCIES = (.25, 1., 2., 4., 8.)


def linear(p, rates, orientation, order):
    """Finite residence-age susceptibility per unit coordinate amplitude."""
    h, v = contact(p, 0.)
    w = np.zeros_like(h)
    for spin in (0, 1):
        w[:, 4 + spin, spin] = -h[:, 4 + spin, spin]
        w[:, spin, 4 + spin] = -h[:, spin, 4 + spin]
    energy, u = np.linalg.eigh(h)
    gap = energy[:, :, None] - energy[:, None, :]
    rho = adj(u) @ RHO_UP @ u
    # Coordinate derivatives include the sign of each face orientation.
    v = (adj(u) @ v @ u) * orientation[:, None, None]
    w = adj(u) @ w @ u

    def population(age):
        # Integral of the unperturbed density matrix between age 0 and age.
        phase = gap * age
        return rho * age * np.exp(-.5j * phase) * np.sinc(phase / (2*np.pi))

    def force(operator):
        return -rates @ tr_product(operator, population(DWELL))

    freq = np.asarray(FREQUENCIES)
    chi = np.full(len(freq), force(w), complex)
    nodes, weights = leggauss(order)
    for node, weight in zip(nodes, weights):
        delay = .5 * DWELL * (node + 1)
        past = population(DWELL - delay)
        comm = v @ past - past @ v
        propagated = np.exp(-1j * gap * delay) * comm
        kernel = rates @ tr_product(1j * propagated, v)
        chi += .5 * DWELL * weight * np.exp(-1j * freq * delay) * kernel
    return force(v), chi


def direct(p, rates, orientation, omega, amplitude, background_tag, phases=24):
    """Solve the time-dependent parcel equation, and book entry and exit drift."""
    theta = (np.arange(phases) * 2*np.pi/phases)[None, :]
    p = p[:, None, :, :]
    signed_rate = orientation[:, None]
    base, _ = contact(p, 0.)
    shape = (len(rates), phases, 6, 3)
    incoming = np.zeros(shape, complex)
    incoming[:, :, ::2, :] = np.eye(3)
    size = incoming.size
    y0 = np.r_[incoming.ravel(), np.zeros(2*len(rates)*phases)]

    def operators(age):
        angle = theta + omega * age
        face_bias = signed_rate * amplitude * np.cos(angle)
        h = np.broadcast_to(base, shape[:2] + (6, 6)).copy()
        dh = np.zeros_like(h)
        closing = -np.exp(-1j * face_bias)
        for spin in (0, 1):
            h[..., 4 + spin, spin] = closing
            h[..., spin, 4 + spin] = closing.conj()
            dh[..., 4 + spin, spin] = -1j * closing * signed_rate
            dh[..., spin, 4 + spin] = dh[..., 4 + spin, spin].conj()
        return angle, h, dh

    def rhs(age, y):
        z = y[:size].reshape(shape)
        angle, h, dh = operators(age)
        rho = z @ adj(z) / 3
        reaction = -tr_product(dh, rho)
        velocity = -omega * amplitude * np.sin(angle)
        return np.r_[(-1j * h @ z).ravel(), (velocity * reaction).ravel(),
                     (np.exp(-1j * angle) * reaction).ravel()]

    sol = solve_ivp(rhs, (0., DWELL), y0, method="DOP853", rtol=2e-12, atol=1e-14)
    if not sol.success:
        raise RuntimeError(sol.message)
    z = sol.y[:size, -1].reshape(shape)
    out = z @ adj(z) / 3
    work, harmonic = sol.y[size:, -1].reshape(2, len(rates), phases)
    entrance = -tr_product(operators(0.)[1]-H0, RHO_UP)
    leave = tr_product(operators(DWELL)[1]-H0, out)
    tag = tr_product(H0, out-RHO_UP)
    drift = entrance + leave
    residual = np.max(abs(tag + drift + work.real))
    power = rates @ work.mean(axis=1).real
    period = 2*np.pi/omega
    return dict(tag=background_tag+rates@tag.mean(axis=1),
                drift=-background_tag+rates@drift.mean(axis=1), power=power,
                work=period*power, chi=(2/amplitude)*(rates@harmonic.mean(axis=1)),
                residual=residual)


def models():
    triad = np.eye(3)
    for label, order in (("triad", triad), ("triad_mirror", triad[::-1])):
        yield label, order[None, :, :], np.array([RATE]), np.ones(1), 0.
    r = calibration()
    pairs, weights, neighbors, orientation, edge = prepared(2, r)
    idx = neighbors[edge]
    print(f"mesh level=2 edge={edge} adjacent={idx.tolist()} signs={orientation[edge].tolist()} "
          f"rates={(RATE*weights[idx]).tolist()} marks_per_face=3")
    for label, p in zip(("mesh", "mesh_mirror"), pairs):
        # Held spectators can convert internal energy at the expense of drift.
        spectator = np.ones(len(p), dtype=bool)
        spectator[idx] = False
        tag0, _, _, _ = read(p[spectator], 0.)
        background = RATE * weights[spectator] @ tag0
        yield (label, np.concatenate([np.roll(p[idx], k, axis=1) for k in range(3)]),
               np.tile(RATE*weights[idx]/3, 3), np.tile(orientation[edge], 3), background)


def main():
    for label, p, rates, sign, background in models():
        torque, chi = linear(p, rates, sign, 48)
        _, coarse = linear(p, rates, sign, 24)
        print(f"response model={label} tau0={torque:.14g} "
              f"quad_error={np.max(abs(chi-coarse)):.3g}")
        assert np.max(abs(chi-coarse)) < 2e-13
        for freq, response in zip(FREQUENCIES, chi):
            print(f"linear model={label} omega={freq:g} chi_re={response.real:.14g} "
                  f"chi_im={response.imag:.14g} power_over_A2={freq*response.imag/2:.14g} "
                  f"work_per_cycle_over_A2={np.pi*response.imag:.14g}", flush=True)
            if freq not in (1., 4.):
                continue
            errors = []
            for a in (.08, -.08, .04):
                exact = direct(p, rates, sign, freq, a, background)
                ratio = exact["work"]/a**2
                dchi = abs(response-exact["chi"])
                difference = ratio - np.pi*response.imag
                print(f"direct model={label} omega={freq:g} A={a:+g} "
                      f"mean_tag_rate={exact['tag']:.14g} mean_drift_rate={exact['drift']:.14g} "
                      f"power_to_coordinate={exact['power']:.14g} "
                      f"work_per_cycle={exact['work']:.14g} work_per_cycle_over_A2={ratio:.14g} "
                      f"chi_re={exact['chi'].real:.14g} chi_im={exact['chi'].imag:.14g} "
                      f"cohort_ledger_error={exact['residual']:.3g} "
                      f"work_coefficient_minus_linear={difference:.4g}", flush=True)
                assert exact["residual"] < 2e-11
                assert abs(exact["work"] - np.pi*a*a*exact["chi"].imag) < 2e-13
                errors.append(dchi)
            print(f"amplitude_check model={label} omega={freq:g} "
                  f"chi_err_A08={errors[0]:.5g} chi_err_minus_A08={errors[1]:.5g} "
                  f"chi_err_A04={errors[2]:.5g}", flush=True)
            assert errors[2] < .27*errors[0] and abs(errors[0]-errors[1]) < 2e-12


if __name__ == "__main__":
    main()
