"""Measure a stationary 4D receiver dipole from radial velocities on two cuts.

Length unit a and velocity unit D(a)/(rho * 2 pi^2 * a^3). Both mouths use
the same units. In the uniform-discharge comparison the harmonic coefficients
are solved from the boundary velocities of both spheres at each separation.
"""
import numpy as np
import sympy as sp
from scipy.special import eval_gegenbauer, roots_chebyu


def symbolic_checks():
    r, inner, outer = sp.symbols("R Ri Ro", positive=True)
    z, g, b, higher = sp.symbols("z g B H", real=True)
    # Sphere mean of an axially symmetric function in four dimensions.
    # Include an odd l=3 harmonic which contaminates a two-pole reading.
    phi = (g*r + b/r**3)*z + higher*sp.chebyshevu(3, z)/r**5
    radial = sp.diff(phi, r)
    moment = sp.simplify(4*sp.integrate(
        2/sp.pi * sp.sqrt(1-z**2)*z*radial, (z, -1, 1)))
    amplitude = sp.factor((moment.subs(r, outer)-moment.subs(r, inner))/
                          (3*(inner**-4-outer**-4)))
    poles = sp.factor((radial.subs(z, 1)-radial.subs(z, -1))/2 - moment)
    print(f"angular moment = {moment}; two-radius extraction = {amplitude}; "
          f"pole-only excess = {poles}", flush=True)
    if amplitude != b or sp.simplify(moment-g+3*b/r**4) != 0 or poles == 0:
        raise ArithmeticError("angular extraction")


def normal_basis(center, radius, separation, mu, degree):
    """Outward derivative on a cut, for both sets of decaying 4D harmonics."""
    sine = np.sqrt(1-mu**2)
    y = radius*sine
    columns = []
    for source in (0.0, separation):
        dx = center + radius*mu - source
        dist = np.hypot(dx, y)
        cosine = dx/dist
        dr = (dx*mu + y*sine)/dist
        for ell in range(degree+1):
            angular = eval_gegenbauer(ell, 1, cosine)
            derivative = 2*eval_gegenbauer(ell-1, 2, cosine) if ell else 0.0
            normal = (-(ell+2)*dr*angular+(mu-cosine*dr)*derivative)/dist**(ell+3)
            columns.append(normal)
    return np.column_stack(columns)


def coefficients(separation, degree, uniform):
    mu, weights = roots_chebyu(3*degree+7)
    basis = np.vstack([normal_basis(c, 1.0, separation, mu, degree) for c in (0.0, separation)])
    coeff = np.zeros(basis.shape[1])
    monopoles = [0, degree+1]
    coeff[monopoles] = -0.5  # Integrated outward flux fixes each monopole.
    if uniform:
        free = np.ones(len(coeff), dtype=bool)
        free[monopoles] = False
        weights = np.sqrt(np.tile(weights, 2))
        target = weights*(1-basis@coeff)
        fit, _, rank, _ = np.linalg.lstsq(weights[:, None]*basis[:, free], target, rcond=None)
        if rank != len(fit):
            raise ArithmeticError("rank deficient boundary solve")
        coeff[free] = fit
    return coeff


def measure(separation, degree, uniform):
    coeff = coefficients(separation, degree, uniform)
    mu, weights = roots_chebyu(103)
    weights = weights*2/np.pi
    boundary_error = 0.0
    rate_error = 0.0
    for c in (0.0, separation):
        velocity = normal_basis(c, 1.0, separation, mu, degree)@coeff
        rate_error = max(rate_error, abs(weights@velocity-1))
        boundary_error = max(boundary_error, np.max(np.abs(velocity-1)))
    moments = []
    for r in (1.2, 1.6):
        velocity = normal_basis(0, r, separation, mu, degree)@coeff
        moments.append(4*weights@(mu*velocity))
    b = (moments[1]-moments[0])/(3*(1.2**-4-1.6**-4))
    fitted_dipole = 2*coeff[1] if degree else 0.0
    return b, moments, rate_error, boundary_error, abs(b-fitted_dipole)


def main():
    symbolic_checks()
    for uniform in (False, True):
        previous = None
        for degree in (14, 24) if uniform else (0,):
            results = []
            worst_rate = worst_bc = worst_extract = 0.0
            for h in (0.1, 0.05, 0.025):
                pair = []
                for d in (4-h, 4+h):
                    value, moments, rate, bc, extract = measure(d, degree, uniform)
                    print(f"uniform={uniform} L={degree} d={d:.3f} "
                          f"M(1.2)={moments[0]:+.12g} M(1.6)={moments[1]:+.12g} B={value:+.12g}")
                    worst_rate = max(worst_rate, rate)
                    worst_bc = max(worst_bc, bc)
                    worst_extract = max(worst_extract, extract)
                    pair.append(value)
                slope = (pair[1]-pair[0])/(2*h)
                results.append(slope)
                print(f"h={h:.3f}  secant={slope:+.12g}")
            print(f"max error in either normalized rate={worst_rate:.4g} "
                  f"max |vN-1| on either mouth={worst_bc:.4g} extraction error={worst_extract:.4g}")
            if worst_rate > 1e-11 or worst_extract > 1e-12:
                raise ArithmeticError("rate or angular projection")
            if uniform and degree == 24:
                refined = results[-1]+(results[-1]-results[-2])/3
                shift = max(abs(a-b) for a, b in zip(results, previous))
                print(f"uniform extrapolated derivative={refined:+.12g}; degree refinement={shift:.4g}")
                if worst_bc > 1e-8 or shift > 1e-9:
                    raise ArithmeticError("unresolved harmonic truncation")
            if not uniform and max(abs(x) for x in results) > 1e-11:
                raise ArithmeticError("transparent monopole control")
            previous = results


if __name__ == "__main__":
    main()
