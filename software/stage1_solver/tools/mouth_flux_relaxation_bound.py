"""Exact energy certificates for angular flux relaxation in the round collar.

Use R=K=u=1. The energy functional is divided by 4*pi in radial integrals,
so an angular integral is an average on the unit sphere. All sector bounds
below are then divided by the similarly normalized energy unit B.
"""

import sympy as sp


def main():
    Q = sp.Rational
    g, aspect, gamma = Q(3, 2), Q(5, 4), Q(1, 5)
    r, z = sp.symbols("r z", positive=True)
    length, volume, inv = [sp.integrate(f, (r, 1, aspect))
                           for f in (1, r**2, r**-2)]
    energy_unit = g**2 * inv / 2
    y = gamma / (4 * energy_unit)  # radial monopole b=n/2
    # Smooth degree-one uniform-connection trial, continued constantly in r.
    f = sp.sqrt(z) * (1 + Q(2, 5) * (1 - z))
    kinetic = z*(1-z)*sp.diff(f, z)**2 + (1-z)*f**2/(4*z)
    T, M, P4 = [sp.integrate(sp.cancel(v), (z, 0, 1))
                for v in (kinetic, f**2, f**4)]
    unit_upper = y - length**2*(g*M-T)**2/(2*volume*P4*energy_unit)
    neutral_gain = g - 1  # |n-C_h|/2 is one for both n=0 components
    wrong_gain = max(g-Q(3, 2), Q(4, 5)*g-Q(1, 2))
    print("B/(4pi), Coulomb/B, unit trial T M P4, E/B:",
          energy_unit, y, T, M, P4, unit_upper)

    # b=r^2 B_r-n/2 has mean zero, while P is the signed material density
    # from the helicity-dependent first-order (Pauli) gradient inequalities.
    # avg(b*P) = avg(b*(P-avg(P))) and Var(P) <= avg(S^2).
    stiffness = sp.Symbol("D", positive=True)
    b, p = sp.symbols("b p", real=True)
    square = stiffness*b**2/r**2 + b*p
    b_star = sp.solve(sp.diff(square, b), b)[0]
    remainder = sp.factor(square - square.subs(b, b_star))
    quartic_fraction = 1 + 2*square.subs(b, b_star)/(r*p)**2
    required_D = sp.solve(-neutral_gain**2/(g**2*quartic_fraction)-unit_upper,
                          stiffness)[0]
    print("curvature square remainder, quartic fraction, D at neutral crossing:",
          remainder, quartic_fraction, required_D)
    # Uniform n benefits from no angular flux stiffness; for zero-mean b,
    # the unit sphere Poincare gap of two applies to every mixture and r.
    chi_simple = (required_D - gamma)/2
    print("plain Maxwell quartic fraction, sufficient chi threshold (strict):",
          quartic_fraction.subs(stiffness, gamma), chi_simple)

    def sector_bounds(effective):
        retained = quartic_fraction.subs(stiffness, effective)
        neutral = -(neutral_gain/g)**2/retained
        opposite = y - (wrong_gain/g)**2/retained
        multiple = 4*y - 1  # rough volume bound retaining the full quartic
        return retained, neutral, opposite, multiple

    sample = sector_bounds(gamma + 2*2)
    print("chi=2 retained quartic, neutral / split opposite / multiple floors:",
          *sample)
    print("chi=2 margins above unit:", *(x-unit_upper for x in sample[1:]))

    # For multipole l of b, a coexact potential contributes |b'|²/[l(l+1)]
    # to the tangential Maxwell term. Beyond the collar complete the square
    # against b'+l*b/r to get the boundary cost b(aspect)^2/[(l+1)*aspect].
    ell = sp.symbols("ell", integer=True, positive=True)
    radial_mode, radial_derivative = sp.symbols("v dv", real=True)
    eigenvalue = ell*(ell+1)
    maxwell_mode = radial_derivative**2/eigenvalue + radial_mode**2/r**2
    boundary_derivative = (2*radial_mode*radial_derivative/r-radial_mode**2/r**2)/(ell+1)
    tail_residual = sp.factor(maxwell_mode + boundary_derivative -
                              (radial_derivative + ell*radial_mode/r)**2/eigenvalue)
    # Weighted triangle and Cauchy-Schwarz give X=int |v|²/r² <= c*Q,
    # Q=int |v'|²/eigenvalue + |v(aspect)|²/[(ell+1)*aspect]. For the
    # derivative term int (aspect-r)/r² <= int (aspect-r) in R=1 units.
    int_bound = sp.integrate(aspect-r, (r, 1, aspect))
    c = inv*aspect*(ell+1) + eigenvalue*int_bound
    D_dipole_Maxwell = gamma*(1+1/c.subs(ell, 1))
    chi_tail = max((required_D-D_dipole_Maxwell)/2, (required_D-gamma)/6)
    print("exterior Maxwell square residual, c(l), dipole Maxwell coefficient:",
          tail_residual, sp.factor(c), D_dipole_Maxwell)
    print("including exterior, sufficient chi threshold (strict):", chi_tail)

    # Base Maxwell law, high angular orders only. Radial ball projection to
    # S<=g is a contraction of the covariant gradients and lowers the local
    # potential where S>g. In angular Coulomb gauge use the uniform Laplacian
    # bound and the elementary inequality for |D_uniform psi-i a psi|².
    t = g/(g + gamma*eigenvalue/aspect**2)
    floor_high = -((neutral_gain+t)/g)**2
    order = 1
    while floor_high.subs(ell, order) <= unit_upper:
        order += 1
    print("Maxwell cutoff comparison (L, t, E0/B floor, gap to old unit):")
    for level in (order-1, order):
        floor = floor_high.subs(ell, level)
        print(level, t.subs(ell, level), floor, floor-unit_upper)

    assert g*M > T and y > 0 and unit_upper < 0
    assert sp.factor(remainder - stiffness*(b-b_star)**2/r**2) == 0
    assert quartic_fraction.subs(stiffness, gamma) < 0
    assert 0 < chi_tail < chi_simple < 2
    assert tail_residual == 0 and all(x > unit_upper for x in sample[1:])
    assert floor_high.subs(ell, order-1) <= unit_upper < floor_high.subs(ell, order)


if __name__ == "__main__":
    main()
