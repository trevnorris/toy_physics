"""Pair energy and radial force for a passive field coupled to mouth handedness.

Assume well-separated, fixed-size spherical mouths with mirrored one-mouth
states of equal self-energy, signed units n_i=-s_i*chi_i, and chi_i=+/-1.
The alignment charge is lambda*chi_i. Its real field has a positive quadratic
energy and is minimized on R^3 with zero field at infinity. Point-source self
energies are subtracted. A force is computed at fixed internal signs. A common
zero-potential charge reservoir is assumed only in comparisons between sectors.
"""

import sympy as sp


def main():
    d, mu, kappa, lam, C = sp.symbols("d mu kappa lambda C", positive=True)
    # Solve the spherical homogeneous field equation with the decaying mode,
    # and normalize via the flux through an infinitesimal surface. The field
    # source factor lambda*chi/kappa multiplies this unit Green function.
    coefficient = sp.Symbol("a", real=True)
    exterior_mode = coefficient*sp.exp(-mu*d)/d
    operator = -sp.diff(d**2*sp.diff(exterior_mode, d), d)/d**2 + mu**2*exterior_mode
    flux = sp.limit(-4*sp.pi*d**2*sp.diff(exterior_mode, d), d, 0, dir="+")
    normalization = sp.solve(flux-1, coefficient)[0]
    green = exterior_mode.subs(coefficient, normalization)
    print("radial operator residual, short-distance flux, normalized Green function:",
          sp.simplify(operator), flux, green)

    s1, s2, chi1, chi2 = sp.symbols("s1 s2 chi1 chi2", real=True)
    n1, n2 = -s1*chi1, -s2*chi2
    # On shell quadratic energy = (1/2)*source.field. The separation-independent
    # self terms in -(1/2)*source.field are removed once at each mouth.
    partner_fields = (lam*chi2*green/kappa, lam*chi1*green/kappa)
    mediated = -lam*(chi1*partner_fields[0]+chi2*partner_fields[1])/2
    pair = sp.factor(mediated+C*n1*n2/d)
    alignment_strength = lam**2/(4*sp.pi*kappa)
    A = sp.Symbol("A", positive=True)
    compact_pair = sp.factor(pair.subs(lam**2, 4*sp.pi*kappa*A))
    force = sp.factor(-sp.diff(compact_pair, d))
    print("alignment amplitude A, mediated cross energy:", alignment_strength, mediated)
    print("pair energy, outward pair force at fixed signs:", compact_pair, force)

    tests = []
    for match in (1, -1):
        aligned = {s1: 1, s2: match, chi1: 1, chi2: 1}
        opposite = {**aligned, chi2: -1}
        difference = sp.factor(compact_pair.subs(aligned)-compact_pair.subs(opposite))
        Fplus = sp.factor(force.subs(aligned))
        print(f"s1*s2={match} aligned chi (n1*n2, E_aligned-E_opposed, force):",
              (n1*n2).subs(aligned), difference, Fplus)
        print(f"s1*s2={match} opposed chi (n1*n2, force):",
              (n1*n2).subs(opposite), force.subs(opposite))
        tests.append(sp.simplify(force.subs(opposite)+Fplus) == 0)

    # Express observables on the like-orientation aligned branch using the
    # local energy advantage delta=A*exp(-mu*d)-C (positive when selected).
    delta = sp.Symbol("delta", positive=True)
    aligned_like = {s1: 1, s2: 1, chi1: 1, chi2: 1}
    chosen = {A: (C+delta)*sp.exp(mu*d)}
    Uchosen = sp.factor(compact_pair.subs(aligned_like).subs(chosen))
    Fchosen = sp.factor(force.subs(aligned_like).subs(chosen))
    marginal = sp.factor(Fchosen.subs(delta, 0))
    print("like orientations, aligned advantage delta>0 (U, F, F at equal energies):",
          Uchosen, Fchosen, marginal)

    # At a finite positive inverse range the extra force vanishes relative to
    # Coulomb at large distance. For zero inverse range both act as 1/d^2.
    print("mu=0 (d*U, d^2*F):", sp.factor((d*compact_pair).subs(mu, 0)),
          sp.factor((d**2*force).subs(mu, 0)))
    far_energy = sp.limit(d*compact_pair, d, sp.oo)
    far_force = sp.limit(d**2*force, d, sp.oo)
    print("fixed mu>0 and d to infinity (d*U, d^2*F):", far_energy, far_force)
    tests.extend((sp.simplify(operator) == 0, sp.simplify(flux.subs(coefficient, normalization)) == 1,
                  sp.simplify(Uchosen+delta/d) == 0,
                  sp.simplify(Fchosen+(delta+(C+delta)*mu*d)/d**2) == 0,
                  marginal.is_negative, sp.simplify(far_force-C*n1*n2) == 0,
                  sp.simplify(far_energy-C*n1*n2) == 0))
    assert all(tests)


if __name__ == "__main__":
    main()
