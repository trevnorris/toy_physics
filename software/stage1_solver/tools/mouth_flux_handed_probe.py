"""Symmetry test for the circular coherence gains on a neutral mouth patch.

The normalized e_h denotes a weak probe. G is the Hermitian quadratic gain
after subtraction of the gradient cost, in the circular convention used for
the mouth. The real axial spin L below is separately specified material data.
A w reflection with paired points held at the same brane coordinates reverses
the projected face normal. Two choices of the paired internal L are tested.
"""

import sympy as sp


def main():
    a, b, c, d = sp.symbols("a b c d", real=True)
    J = sp.Matrix([[0, -1], [1, 0]])
    F = sp.diag(1, -1)
    gain = sp.Matrix([[a, c-sp.I*d], [c+sp.I*d, b]])
    # Commuting with the generator gives invariance under every planar
    # rotation. Imposing the tangent reflection then tests the odd coefficient.
    commutator = gain*J-J*gain
    rotation_solution = sp.solve(list(commutator), (b, c), dict=True)[0]
    isotropic = gain.subs(rotation_solution)
    reflected = F*isotropic*F.T
    reflection_solution = sp.solve(list(reflected-isotropic), d, dict=True)[0]
    print("rotation commutator, rotation constraints:", commutator, rotation_solution)
    print("isotropic G, G under tangent reflection, fixed-state constraint:",
          isotropic, reflected, reflection_solution)

    gains = {}
    helicity_swaps = []
    for h in (1, -1):
        e = sp.Matrix([1, sp.I*h])/sp.sqrt(2)
        gains[h] = sp.expand((e.H*isotropic*e)[0])
        helicity_swaps.append(F*e == e.conjugate())
    gf, gw = sp.Rational(9, 5), sp.Rational(1, 10)
    contrast_fit = sp.solve([gains[1]-gf, gains[-1]-gw], (a, d), dict=True)[0]
    contrast = sp.factor((gains[1]-gains[-1])/(gains[1]+gains[-1]))
    print("g_+, g_-, reflected probe helicities swapped:", gains[1], gains[-1],
          *helicity_swaps)
    print("chosen favored plus, mean and odd coefficients, normalized contrast:",
          contrast_fit[a], contrast_fit[d], contrast.subs(contrast_fit))

    # The standing-wave energy at fixed total circular action does not bias
    # which imbalance is prepared. Identifying wave spin with a material gain
    # would require an independent response coupling.
    total, imbalance, omega = sp.symbols("I S omega", real=True)
    actions = ((total+imbalance)/2, (total-imbalance)/2)
    energy = sp.expand(omega*sum(actions))
    slope = sp.diff(energy, imbalance)
    print("degenerate wave at fixed total action, energy and imbalance slope:",
          energy, slope)

    # Isotropy permits the lab z axis to be chosen along the fixed axial L.
    # n_s=s*rhat; average via x=cos(theta), uniform on [-1,1].
    x, chi = sp.symbols("x chi", real=True)
    spin = sp.Symbol("L", positive=True)
    def average(f):
        return sp.integrate(f, (x, -1, 1))/2

    projection_checks = []
    for s in (1, -1):
        fixed_projection = s*spin*x
        extrema = (fixed_projection.subs(x, -s), fixed_projection.subs(x, s))
        print(f"s={s} fixed axial L, <L.n>, <(L.n)^2>, min, max:",
              average(fixed_projection), average(fixed_projection**2), *extrema)
        projection_checks.extend((average(fixed_projection) == 0,
                                  average(fixed_projection**2) == spin**2/3,
                                  extrema == (-spin, spin)))
    # At the same point use radial components for n and L. A real axial field
    # transforms as det(F3)*F3 under an ordinary brane tangent reflection.
    F3 = sp.diag(1, -1, 1)
    normal = sp.Matrix([0, 0, 1])
    axial = sp.Matrix([0, 0, chi])
    reflected_projection = (F3*normal).dot(F3.det()*F3*axial)
    print("L=chi*n spherical mean, mean square, min, max, tangent-reflected projection:",
          average(chi), average(chi**2), chi, chi, reflected_projection)

    # h is the initially favored helicity relative to n, with nonzero coupling
    # and nonzero chi. The energy certificate then selects n_flux=-h*s.
    s, h = sp.symbols("s h", real=True)
    for name, axial_sign in (("held", 1), ("reversed", -1)):
        # At a given paired brane point the projected normal reverses.
        # The relative axial sign is an additional preparation choice.
        projection_ratio = sp.cancel((-normal).dot(axial_sign*axial)/normal.dot(axial))
        new_h = projection_ratio*h
        new_flux = sp.expand(-new_h*(-s))
        print("w mirror axial texture", name, "(new chi, h, unit flux/original):",
              projection_ratio*chi, new_h, sp.cancel(new_flux/(-h*s)))
    print("independent +/-chi at fixed s (two selected unit fluxes):", -h*s, h*s)

    assert isotropic == a*sp.eye(2)+sp.I*d*J
    assert reflected == isotropic.subs(d, -d) and reflection_solution == {d: 0}
    assert all(helicity_swaps) and slope == 0
    assert all(projection_checks) and reflected_projection == -chi


if __name__ == "__main__":
    main()
