"""Material sites constrained to a translating sphere, free to slide tangentially.

One isolated radial mouth. Uniform sites, no local replacement, independent
areal mass and dilation modulus. Axisymmetric dipole, harmonic spin response,
a shared bond between both materials, and the outgoing angular wave. Coordinates
in physical lengths; u is displacement along the tangential projection of z.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_pair_constituent_response import BOND, EXCHANGE, HBAR, PIN, analytic_site, site
from mouth_pair_translating_cuts import MODULUS, SPEED, hankel


AREA_MASS, AREA_MODULUS = 10., .06


def geometry(radius, eta, modulus):
    """Differentiate displaced marked latitude rings before reducing to a dipole."""
    v, weight = leggauss(80)
    theta = np.arccos(v)
    sin = np.sin(theta)
    # Integral divided by (2/3)*surface area.
    avg = lambda f: .75*(weight@f)
    step = radius*1.e-4
    def moved(u):
        angle = theta-u*sin/radius
        stretch = np.sin(angle)/sin*(1-u*v/radius)
        r = radius*np.stack([np.sin(angle), np.cos(angle)])
        turn = eta*(v-np.cos(angle))/(2*radius*sin)
        return r, stretch, turn
    plus, minus = moved(step), moved(-step)
    dr, dilation, turn = [(a-b)/(2*step) for a, b in zip(plus,minus)]
    overlap = AREA_MASS*avg(dr[1])
    mass = AREA_MASS*avg(np.sum(dr**2, axis=0))
    stiffness = modulus*avg(dilation**2)
    native = -avg(turn*sin)
    density_shape = np.max(abs(dilation+2*v/radius))
    print(f"geometry R={radius:g} eta={eta:+d} K_A={modulus:g} M_Xu={overlap:.10g} "
          f"M_uu={mass:.10g} K_uu={stiffness:.10g} b={native:.10g} "
          f"max_dilation_shape_error={density_shape:.3g}")
    assert density_shape<1e-8 and abs(np.dot(weight, dilation))<1e-11
    return overlap, mass, stiffness, native


def field(radius, omega):
    kR = omega*radius/SPEED
    return -(1+kR*hankel(1, kR, True)/hankel(1,kR))


def solve(radius, eta, ratio, modulus, moments):
    pol, gap = analytic_site(EXCHANGE)
    omega = ratio*gap/HBAR
    surface = MODULUS*field(radius,omega)/radius
    static = MODULUS/radius
    h = 4*radius**2*pol*PIN
    m = h*(HBAR/gap)**2
    b = eta/(2*radius**2)
    # Unknowns are exterior angle a, geometric-frame slip s, tangential shift u.
    matrix = np.array([
        [surface+BOND, -BOND, 0.],
        [-BOND, BOND+h-m*omega**2, m*b*omega**2],
        [0., m*b*omega**2, 2*modulus/radius**2-(AREA_MASS+m*b*b)*omega**2]])
    drive = np.array([b*(static+BOND), -BOND*b, AREA_MASS*omega**2])
    a, s, u = np.linalg.solve(matrix,drive)  # all three responses per unit X
    # Direct alternative: coordinates a, tagged-site angle U=s-b*u, and u.
    # The spin is driven by local pin rotation and the bond's return torque;
    # the material shift is driven by both its actual mass and contact force.
    overlap, mass, stiffness, b_site = moments
    polarization, spin_chi, _ = site(EXCHANGE, eta, ratio*gap)
    comp = spin_chi/(2*radius*polarization)**2
    pin = 4*radius**2*polarization*PIN
    full = np.array([
        [surface+BOND, -BOND, -BOND*b_site],
        [-comp*BOND, 1+comp*BOND, b_site*(comp*BOND+spin_chi*PIN/polarization)],
        [-b_site*BOND, b_site*(BOND+pin), stiffness-mass*omega**2+b_site**2*(BOND+pin)]])
    source = [b*(static+BOND), -comp*BOND*b, overlap*omega**2-b_site*BOND*b]
    direct, tag, shift = np.linalg.solve(full,source)
    error = max(abs(direct/a-1), abs(shift/u-1), abs(tag+b_site*shift-s)/(abs(a)+abs(s)))
    fixed_a, _ = np.linalg.solve(matrix[:2,:2],drive[:2])
    torque = BOND*(a-b-s)
    # Force applied per (2/3) area and per X. Terms in phase with X do no mean work.
    force = -b*(static*a+torque)-omega**2*AREA_MASS*(1.5+u)
    work = (1j*omega*force).real/2
    flux = -surface.imag*omega*abs(a)**2/2
    balance = abs(work/flux-1)
    print(f"mobile R={radius:g} eta={eta:+d} K_A={modulus:g} freq_over_gap={ratio:g} "
          f"u_over_X=({u.real:.9g},{u.imag:.9g}) "
          f"far_over_carried_abs={abs(a/fixed_a):.9g} far_over_carried_phase={np.angle(a/fixed_a):.9g} "
          f"direct_error={error:.3g} wave_over_drive={flux/work:.10g}")
    assert error<2e-6 and balance<1e-6 and work>0
    return a, u


def main():
    for radius in (5., 10.):
        moments = {eta: geometry(radius, eta, AREA_MODULUS) for eta in (1, -1)}
        for ratio in (.05, .30, .65, 1.15):
            plus = solve(radius, 1, ratio, AREA_MODULUS, moments[1])
            minus = solve(radius, -1, ratio, AREA_MODULUS, moments[-1])
            assert abs(plus[0]/minus[0]+1)<1e-10 and abs(plus[1]/minus[1]-1)<1e-10
        rigid = 1e6*AREA_MODULUS
        stiff = solve(radius,1,.65,rigid,geometry(radius,1,rigid))
        assert abs(stiff[1])<1e-5


if __name__=="__main__":
    main()
