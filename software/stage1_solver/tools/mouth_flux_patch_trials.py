"""Exact admissible collar trials with a dipolar electric-flux perturbation.

Use rho=r/R, Psi=sqrt(K/u) f/R and energy units 4*pi*K**2/(u*R).
The connection in the north chart is m*(1-z)+4*z*(1-z)*b(rho)
with angular index k, z=(1+cos(theta))/2 and m=n-C_h.  Each component
has its own phase exp(i*k*phi).  Profiles f are constant in radius.
Free material boundaries and a free inner flux profile are understood.
The exterior vacuum dipole minimizes the field energy for its trace.
"""

import itertools

import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize


z, rho = sp.symbols("z rho", positive=True)
aspect = sp.Rational(5, 4)
g = sp.Rational(3, 2)
weak = sp.Rational(6, 5)
gamma = sp.Rational(1, 5)  # gamma*u/K**2
area_weight = z*(1-z)


def integral(expr, var, lo, hi):
    return sp.factor(sp.integrate(sp.cancel(expr), (var, lo, hi)))


length = aspect - 1
volume = integral(rho**2, rho, 1, aspect)
energy_unit = g**2 * integral(rho**-2, rho, 1, aspect) / 2  # B3 in our units
v = sp.Matrix([1, rho**2-1])
trace = v.subs(rho, aspect)
radial = lambda t: integral(t, rho, 1, aspect)
mean_v = v.applyfunc(radial)
mean_vv = (v*v.T).applyfunc(radial)
# Angular averages of the two physical field components.  Relative to
# radial monopole n/(2*rho**2), these are 2*b*cos(theta)/rho**2 and
# -b'*sin(theta)/rho, giving zero integrated cross term with the monopole.
normal = integral((2*(2*z-1))**2, z, 0, 1)
tangential = integral(4*area_weight, z, 0, 1)
tail = aspect/rho
tail_cost = integral(normal*tail**2/rho**2 + tangential*tail.diff(rho)**2,
                     rho, aspect, sp.oo)
magnetic = gamma*((normal*v*v.T/rho**2 + tangential*v.diff(rho)*v.diff(rho).T
                   ).applyfunc(radial) + tail_cost*trace*trace.T)


def section_prefactor(m, k):
    return z**(sp.Rational(abs(k-m), 2))*(1-z)**(sp.Rational(abs(k), 2))


def exact_trial(n, gains, indices, coefficients, relaxed):
    density = 0
    angular = 0
    current = 0
    for m, k, gain, coeff in zip((n+2, n-2), indices, gains, coefficients):
        p = sum(c*(2*z-1)**j for j, c in enumerate(coeff))
        # Express f' via its logarithmic prefactor derivative to retain
        # polynomial integrands at the south- and north-chart poles.
        pref2 = section_prefactor(m, k)**2
        dp = p.diff(z) + p*(abs(k-m)/(2*z)-abs(k)/(2*(1-z)))
        f2 = pref2*p**2
        az = k-m*(1-z)
        density += f2
        angular += integral(area_weight*pref2*dp**2 +
                            (az**2/(4*area_weight)-gain)*f2, z, 0, 1)
        current += integral(az*f2, z, 0, 1)
    quadratic = magnetic + mean_vv*integral(4*area_weight*density, z, 0, 1)
    source = mean_v*current
    b = quadratic.inv()*source if relaxed else sp.zeros(2, 1)
    monopole = gamma*n**2/4
    field_patch = (b.T*magnetic*b)[0]
    matter = (length*angular + volume*integral(density**2, z, 0, 1)/2
              + (b.T*(quadratic-magnetic)*b)[0] - 2*(b.T*source)[0])
    return tuple(sp.factor(value/energy_unit) for value in
                 (monopole, field_patch, matter, monopole+field_patch+matter)), b


# Floating optimization only proposes rational section coefficients. Every
# energy subsequently reported as a witness is integrated and checked exactly.
nodes, weights = leggauss(20)
nodes = (1+nodes)/2
weights = weights/2
w = nodes*(1-nodes)
mv = np.array(mean_v, dtype=float).ravel()
mvv = np.array(mean_vv, dtype=float)
mag = np.array(magnetic, dtype=float)


def search(n, gains, relaxed):
    candidates = []
    for indices in itertools.product(*(range(min(0, m), max(0, m)+1)
                                       for m in (n+2, n-2))):
        basis, derivative, azimuth = [], [], []
        for m, k in zip((n+2, n-2), indices):
            a, b = abs(k-m)/2, abs(k)/2
            pref = nodes**a*(1-nodes)**b
            polys = np.array([np.ones_like(nodes), 2*nodes-1, (2*nodes-1)**2])
            dpolys = np.array([np.zeros_like(nodes), 2*np.ones_like(nodes), 4*(2*nodes-1)])
            basis.append(polys*pref)
            derivative.append(pref*(dpolys + polys*(a/nodes-b/(1-nodes))))
            azimuth.append(k-m*(1-nodes))
        basis, derivative, azimuth = map(np.array, (basis, derivative, azimuth))
        pot = azimuth**2/(4*w)-np.array(gains, dtype=float)[:, None]

        def cost_and_gradient(vec):
            coeff = vec.reshape(2, 3)
            f = np.einsum("hj,hji->hi", coeff, basis)
            df = np.einsum("hj,hji->hi", coeff, derivative)
            den = np.sum(f**2, axis=0)
            linear = mv*np.sum(weights*azimuth*f**2)
            hessian_half = mag + mvv*np.dot(weights, 4*w*den)
            dip = np.linalg.solve(hessian_half, linear) if relaxed else np.zeros(2)
            d1 = float(length)*w*df
            d0 = f*(float(length)*pot + float(volume)*den +
                    (dip @ mvv @ dip)*4*w - 2*(mv @ dip)*azimuth)
            gradient = 2*np.sum(weights*(d1[:, None, :]*derivative +
                                         d0[:, None, :]*basis), axis=2)
            energy = (float(gamma)*n**2/4 + float(length)*np.sum(weights*(w*df**2 + pot*f**2))
                      + float(volume)*np.dot(weights, den**2)/2 - dip @ linear)
            return energy, gradient.ravel()

        # Include an interior start as well as each single-helicity face.
        for start in ((1, 1), (1, 0), (0, 1)):
            initial = np.array([[s, 0, 0] for s in start], dtype=float)
            res = minimize(cost_and_gradient, initial.ravel(), method="BFGS", jac=True,
                           options={"gtol": 1e-10})
            candidates.append((res.fun, indices, res.x.reshape(2, 3)))
    numerical, indices, coeff = min(candidates, key=lambda row: row[0])
    rational = [[sp.Rational(round(p*1000), 1000) for p in row] for row in coeff]
    parts, dip = exact_trial(n, gains, indices, rational, relaxed)
    return indices, rational, parts, dip, abs(float(parts[-1]*energy_unit)-numerical)


def main():
    print("units (B3, gamma), radial v, mean(v), mean(vv):",
          energy_unit, gamma, list(v), list(mean_v), mean_vv)
    print("angular Maxwell coefficients (normal, tangential), exterior dipole cost, dipole matrix:",
          normal, tangential, tail_cost, magnetic)
    print("vacuum dipole Euler residual, shell net perturbation:",
          sp.factor(tangential*tail.diff(rho, 2) - normal*tail/rho**2),
          integral(2*(2*z-1), z, 0, 1))
    checks = []
    for n, gains in ((0, (g, g)), (0, (g, weak)), (-1, (g, g)),
                     (-1, (g, weak)), (1, (g, weak))):
        energies = []
        for relaxed in (False, True):
            indices, coeff, parts, dip, rounding = search(n, gains, relaxed)
            print(f"n={n} gains={gains} relaxed={relaxed} k={indices} coeff={coeff}")
            print("  dipole coefficients:", list(dip))
            print("  E/B3 monopole, patch, material, total (exact):", *parts)
            print("  E/B3 monopole, patch, material, total (decimal):",
                  *(f"{float(p):.10f}" for p in parts))
            energies.append(parts[-1])
            checks.extend((parts[1] >= 0, parts[0] >= 0, rounding < 1e-5))
        checks.append(energies[1] <= energies[0])
    if not all(checks):
        raise ArithmeticError("A reported trial value failed its conditioning or normalization check.")


if __name__ == "__main__":
    main()
