"""Check order stress, energy, and curvature on the stationary rounded front."""

import sympy as sp


def compare(label, left, right):
    residual = sp.factor(left - right)
    print(f"{label} | {left} | {right} | difference {residual}", flush=True)
    if residual != 0:
        raise AssertionError(label)


def diffuse_layer():
    s = sp.symbols("s", real=True)
    chi = sp.Function("chi")(s)
    well = sp.Function("W")(chi)
    area = sp.Function("A", positive=True)(s)
    kappa, mobility, n = sp.symbols("kappa M n", positive=True)
    speed, source = sp.symbols("u_N S", real=True)
    grad = chi.diff(s)
    curvature = area.diff(s) / area
    energy = well + kappa * grad**2 / 2
    integrand = area * energy
    mu = (integrand.diff(chi) - integrand.diff(grad).diff(s)) / area
    normal_stress = energy - kappa * grad**2
    tangent_stress = energy
    divergence = normal_stress.diff(s) + curvature * (normal_stress - tangent_stress)
    compare("parcel force and normal-tube stress divergence", mu * grad, divergence)

    # Reduce on a locally planar kink; H may vary on its normal tube.
    kink_mu = sp.simplify(mu.subs(chi.diff(s, 2), well.diff(chi) / kappa))
    compare("kink force density", kink_mu * grad, -kappa * curvature * grad**2)
    compare("stationary passive crossing speed", -mobility * kink_mu / grad,
            mobility * kappa * curvature)

    # Separately close the total conversion rate with a possible external source.
    driven_mu = (source - speed * grad) / mobility
    parcel_power = speed * driven_mu * grad
    compare("parcel power from the order force", parcel_power,
            source * driven_mu - mobility * driven_mu**2)
    number_flux = n * speed
    compare("passive heat per volume divided by number flux",
            (-parcel_power / number_flux).subs(source, 0), speed * grad**2 / (mobility * n))


def stationary_graph():
    r = sp.symbols("r", positive=True)
    h = sp.Function("h")(r)
    sigma, n = sp.symbols("sigma n", positive=True)
    slope = h.diff(r)
    g = sp.sqrt(1 + slope**2)
    area_measure = 4 * sp.pi * r**2
    sheet_energy = area_measure * sigma * g
    edge_gradient = sheet_energy.diff(slope)
    curvature = sp.factor(-edge_gradient.diff(r) / (area_measure * sigma))
    compare("graph curvature from its area variation", curvature,
            -h.diff(r, 2) / g**3 - 2 * slope / (r * g))
    tangent_curvature = sp.factor(curvature.subs(h.diff(r, 2), -2 * slope / r))
    compare("rim curvature with zero linear Laplacian", tangent_curvature,
            -2 * slope**3 / (r * g**3))

    rate, trial_k = sp.Rational(1, 4), sp.Rational(2, 5)
    lip = {slope: -rate / (4 * sp.pi * trial_k), r: 1}
    tension = sp.Rational(1, 3)
    linear_load = tension * rate / trial_k
    exact_load = (-edge_gradient).subs(lip).subs(sigma, tension)
    print("R=1 projected linear load, geometric axial load, geometric rim curvature |",
          linear_load, float(exact_load), float(tangent_curvature.subs(lip)), flush=True)
    compare("exact axial load", exact_load, linear_load / g.subs(lip))

    q = sp.symbols("q", real=True)
    # q is number per projected area per time; u_N uses true surface area.
    sheet_work = (tension * q / trial_k) * (q / (n * g)) * g
    compare("positive work required per projected area in the trial", sheet_work,
            tension * q**2 / (trial_k * n))
    projected_volume = sp.integrate(area_measure, (r, 0, 1))
    print("unit-density Cauchy lower bound for trial power |",
          tension * rate**2 / (trial_k * projected_volume), flush=True)


if __name__ == "__main__":
    diffuse_layer()
    stationary_graph()
