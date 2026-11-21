"""Axisymmetric four-dimensional incompressible potential flow through a smooth neck.

The reflecting wall is the coordinate surface v=1/2, and the radius-one
parabolic cut touches it at u=-1. Coordinates use an unshifted neck plane w=0.
Number density and mass density are one. Integrate the Cartesian gradient
across the physical cut separately from the confocal flux normalization.
"""

from math import pi, sqrt

import sympy as sp


def symbolic_checks():
    u, v, t, v0 = sp.symbols("u v t v0", real=True)
    d = 1 - v0**2
    # Divide the graph equation by a, substitute both confocal coordinates
    # and clear the nonzero cap-slope denominator 2 t d.
    radius_sq = (1 + u**2) * (1 - v**2)
    lip_radius_sq = (1 + t**2) * d
    cap_eq = 2*t*d*(u*v/v0 + t) - (lip_radius_sq - radius_sq)
    # The derivative along the streamline is proportional to q.
    sign_factor = sp.diff(cap_eq, u)/2
    square = sp.factor(sign_factor**2 - (1-v**2)*cap_eq)
    print(f"crossing_factor_squared_on_cap={square}")
    assert sp.simplify(square - (v**2-v0**2)*(1-v**2 + t**2*d/v0**2)) == 0

    r = sp.symbols("r", real=True)
    h = -r**2/sp.sqrt(6)
    clearance = sp.factor(sp.Rational(1, 2) + 3*h**2 - r**2)
    print(f"wall_radius_squared_minus_cut_radius_squared={clearance}")
    assert clearance == (r-1)**2*(r+1)**2/2


A = sqrt(2/3)
TARGET_FLUX = 0.25
# Integral of sqrt(1-v*v) from v=1/2 to v=1.
ANGULAR_INTEGRAL = pi/6 - sqrt(3)/8
AMPLITUDE = TARGET_FLUX/(4*pi*A*A*ANGULAR_INTEGRAL)


def current(r, w):
    """Differentiate the inverse-coordinate phase at fixed laboratory r,w."""
    k = 1 - (r*r+w*w)/(A*A)
    gap = sqrt(k*k + 4*(w/A)**2)
    y = (k + gap)/2  # v squared, positive throughout the neck
    u = w/(A*sqrt(y))
    dy_dr = -2*r*y/(A*A*gap)
    dy_dw = 2*w*(1-y)/(A*A*gap)
    du_dr = -u*dy_dr/(2*y)
    du_dw = 1/(A*sqrt(y)) - u*dy_dw/(2*y)
    dphase_du = AMPLITUDE/(1+u*u)**1.5
    return dphase_du*du_dr, dphase_du*du_dw


def cap_current(r):
    jr, jw = current(r, -r*r/sqrt(6))
    slope = -2*r/sqrt(6)
    return jr, jw, jw-slope*jr


def integrate_cut(n):
    terms = []
    for i in range(n+1):
        r = i/n
        weight = 1 if i in (0, n) else 2 + 2*(i % 2)
        terms.append(weight * r*r * cap_current(r)[2])
    return (4*pi/(3*n))*sum(terms)


if __name__ == "__main__":
    symbolic_checks()
    print(f"amplitude={AMPLITUDE:.14g} phase_difference={2*AMPLITUDE:.14g}")
    for r in (0, 0.5, 0.9, 0.99, 0.999, 1):
        jr, jw, q = cap_current(r)
        print(f"r={r:g} j_r={jr:+.12g} j_w={jw:+.12g} q={q:+.12g}")
    for intervals in (256, 1024):
        rate = integrate_cut(intervals)
        print(f"cut_intervals={intervals} inward_rate={rate:.14g}")
    wall_slip = []
    for u in (-2, -1, 0, 1, 2):
        w = A*u/2
        r = sqrt(0.5+3*w*w)
        jr, jw = current(r, w)
        wall_slip.append(abs(r*jr - 3*w*jw)/sqrt(r*r+9*w*w))
    print(f"max_wall_normal_current={max(wall_slip):.12g}")
    assert abs(rate-TARGET_FLUX) < 1e-10
    assert abs(cap_current(1)[2]) < 1e-14
    assert max(wall_slip) < 1e-14
