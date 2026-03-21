
import sympy as sp
import numpy as np

# ---------------------------------------------------------------------------
# 2PN Family-1 TF inertia scale from the parent n=5 GNLS / hydrodynamic PDE
# ---------------------------------------------------------------------------

def section(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def main() -> None:
    # ------------------------------------------------------------------
    # Symbols
    # ------------------------------------------------------------------
    nEOS = sp.symbols('nEOS', positive=True)
    a0, Lam = sp.symbols('a0 Lam', positive=True)
    Sigma, rhoGeom, betaGeom = sp.symbols('Sigma rhoGeom betaGeom', positive=True)
    mPsi, KEOS, OmegaIn = sp.symbols('mPsi KEOS OmegaIn', positive=True)
    s = sp.symbols('s', real=True)

    # ------------------------------------------------------------------
    # 1) Family-1 interior TF profile from the frozen n-polytrope
    #
    # Parent hydrostatic/GNLS closure inside the throat:
    #   h(rho) + (1/2) mPsi OmegaIn^2 w^2 = muTF
    # with
    #   h(rho) = n K rho^(n-1)/(n-1).
    #
    # On the "filled-to-endcap" branch, the TF support reaches |w| = L/2,
    # so muTF = (1/2) mPsi OmegaIn^2 (L/2)^2.
    # ------------------------------------------------------------------
    alpha = sp.simplify(1 / (nEOS - 1))
    u = sp.symbols('u', real=True)

    c0 = sp.simplify(
        sp.Rational(1, 2) * sp.integrate((1 - u**2)**alpha, (u, -1, 1))
    )
    c2 = sp.simplify(
        sp.Rational(1, 2) * sp.integrate(u**2 * (1 - u**2)**alpha, (u, -1, 1))
    )
    c2OverC0 = sp.simplify(sp.factor(sp.together(sp.simplify(c2 / c0))))
    c2OverC0Closed = sp.simplify(1 / (2 * alpha + 3))
    mLLGeneral = sp.simplify(sp.factor(c2 / (4 * c0)))
    mLLClosed = sp.simplify((nEOS - 1) / (4 * (3 * nEOS - 1)))

    L0 = Lam * a0
    rhoCenter = sp.simplify(
        (((nEOS - 1) * mPsi * OmegaIn**2 * L0**2) / (8 * nEOS * KEOS))**alpha
    )
    rhoEffTF = sp.simplify(c0 * rhoCenter)

    mHatGeneral = sp.Matrix([
        [sp.Rational(3, 5), 0],
        [0, sp.simplify(mLLGeneral)]
    ])

    # n = 5 specialization
    nFive = sp.Integer(5)
    c0n5 = sp.simplify(sp.factor(c0.subs(nEOS, nFive)))
    rhoCenterN5 = sp.simplify(sp.factor(rhoCenter.subs(nEOS, nFive)))
    rhoEffN5 = sp.simplify(sp.factor(rhoEffTF.subs(nEOS, nFive)))
    mHatN5 = sp.simplify(mHatGeneral.subs(nEOS, nFive))

    # ------------------------------------------------------------------
    # 2) Carry forward the same static geometry Hessian as the previous step
    # ------------------------------------------------------------------
    a, L = sp.symbols('a L', positive=True)
    V = sp.Rational(4, 3) * sp.pi * a**3 * L
    A = 4 * sp.pi * a**2 * L + sp.Rational(8, 3) * sp.pi * a**3

    sigma = Sigma / a0**3
    Pvac = rhoGeom * Sigma / a0**4
    kappab = betaGeom * Sigma / a0
    Egeom = sp.expand(Pvac * V + sigma * A + kappab * a**2 / L)

    H = sp.hessian(Egeom, (a, L))
    g = sp.Matrix([sp.diff(V, a), sp.diff(V, L)])

    subs0 = {a: a0, L: Lam * a0}
    H0 = sp.simplify(H.subs(subs0))
    V0 = sp.simplify(V.subs(subs0))
    g0 = sp.simplify(g.subs(subs0))

    hBar = sp.simplify(H0 * a0**2 / Sigma)
    gBar = sp.simplify(g0 / V0 * a0)

    # Static geometry channel (independent of inertia model)
    Delta0 = sp.simplify((gBar.T * hBar.inv() * gBar)[0] / Sigma)

    # Dynamic response on the Family-1 TF inertia branch
    Ytf = sp.simplify((gBar.T * (hBar - s * mHatN5).inv() * gBar)[0] / Sigma)
    Delta2tf = sp.simplify((gBar.T * hBar.inv() * mHatN5 * hBar.inv() * gBar)[0] / Sigma)
    lamEffTF = sp.simplify(Delta0 / Delta2tf)

    # ------------------------------------------------------------------
    # 3) EM-worked point
    # ------------------------------------------------------------------
    x01 = sp.Float('2.40482555769577276862163187933', 50)
    LamEM = sp.N(sp.sqrt(2) * sp.pi / x01, 50)
    rhoEx = sp.Rational(1, 10)
    betaEx = sp.Integer(12)
    target = sp.Rational(109, 280)

    DeltaUnit = sp.N(Delta0.subs({a0: 1, Lam: LamEM, rhoGeom: rhoEx, betaGeom: betaEx, Sigma: 1}), 60)
    SigmaStar = sp.N(DeltaUnit / target, 60)

    hNum = np.array(
        sp.N(hBar.subs({a0: 1, Lam: LamEM, rhoGeom: rhoEx, betaGeom: betaEx}), 60).tolist(),
        dtype=float
    )
    mNum = np.array(
        sp.N(mHatN5, 60).tolist(),
        dtype=float
    )
    gNum = np.array(
        sp.N(gBar.subs({a0: 1, Lam: LamEM}), 60).tolist(),
        dtype=float
    ).flatten()

    evals, evecs = np.linalg.eig(np.linalg.solve(mNum, hNum))
    order = np.argsort(evals)
    evals = np.real(evals[order])
    evecs = np.real(evecs[:, order])

    # Mass normalize modes
    for i in range(evecs.shape[1]):
        norm = np.sqrt(evecs[:, i].T @ mNum @ evecs[:, i])
        evecs[:, i] /= norm

    residues = (gNum @ evecs) ** 2 / (float(SigmaStar) * evals)
    residueFracs = residues / residues.sum()
    lamEffTFNum = float(sp.N(lamEffTF.subs({a0: 1, Lam: LamEM, rhoGeom: rhoEx, betaGeom: betaEx}), 40))

    # Physical scaling once rhoEffTF is inserted
    rhoEffWorked = sp.simplify(rhoEffN5.subs({Lam: LamEM, a0: 1}))
    omegaMinusSq = sp.simplify(evals[0] * Sigma / (rhoEffWorked * V0.subs({a0: 1, Lam: LamEM})))
    omegaPlusSq = sp.simplify(evals[1] * Sigma / (rhoEffWorked * V0.subs({a0: 1, Lam: LamEM})))
    omegaEffSq = sp.simplify(lamEffTFNum * Sigma / (rhoEffWorked * V0.subs({a0: 1, Lam: LamEM})))

    # Low-frequency Padé quality
    sgrid = np.linspace(0.0, 0.1 * evals[0], 400)
    exactGrid = sum(Ri / (1.0 - sgrid / li) for Ri, li in zip(residues, evals))
    padeGrid = float(target) / (1.0 - sgrid / lamEffTFNum)
    relErr = np.max(np.abs((padeGrid - exactGrid) / exactGrid))

    section("1) Family-1 TF bulk inertia from the parent n-polytrope")
    print("alpha = 1/(nEOS-1) =")
    sp.pprint(alpha)
    print("")
    print("c0(n) = (1/2) ∫_{-1}^1 (1-u^2)^alpha du =")
    sp.pprint(sp.factor(c0))
    print("")
    print("c2(n) = (1/2) ∫_{-1}^1 u^2 (1-u^2)^alpha du =")
    sp.pprint(sp.factor(c2))
    print("")
    print("Exact ratio c2/c0 =")
    sp.pprint(sp.factor(c2OverC0))
    print("Closed form =")
    sp.pprint(sp.factor(c2OverC0Closed))
    print("")
    print("General TF axial inertia coefficient mLL(n) = c2/(4 c0) =")
    sp.pprint(sp.factor(mLLGeneral))
    print("Closed form =")
    sp.pprint(sp.factor(mLLClosed))
    print("")
    print("Central TF density on the filled-to-endcap branch =")
    sp.pprint(sp.factor(rhoCenter))
    print("")
    print("Effective bulk inertia density rhoEffTF = c0 * rhoCenter =")
    sp.pprint(sp.factor(rhoEffTF))
    print("")
    print("General reduced inertia metric mHatGeneral =")
    sp.pprint(mHatGeneral)

    section("2) Frozen n = 5 specialization")
    print("c0(5) =")
    sp.pprint(sp.factor(c0n5))
    print("≈", sp.N(c0n5, 20))
    print("")
    print("rhoCenter(n=5) =")
    sp.pprint(sp.factor(rhoCenterN5))
    print("")
    print("rhoEffTF(n=5) =")
    sp.pprint(sp.factor(rhoEffN5))
    print("")
    print("mHatTF(n=5) =")
    sp.pprint(mHatN5)
    print("")
    print("Check: axial inertia renormalizes from the uniform 1/12 to")
    print("  ", sp.simplify(mHatN5[1, 1]))

    section("3) Geometry breathing response with TF inertia")
    print("Dimensionless static Hessian hBar =")
    sp.pprint(hBar)
    print("")
    print("Dimensionless volume-coupling vector gBar =")
    sp.pprint(gBar)
    print("")
    print("Static coefficient Delta0 (unchanged by inertia model) =")
    sp.pprint(sp.factor(Delta0))
    print("")
    print("Exact TF dynamic monopole susceptibility Ytf(s) =")
    sp.pprint(sp.factor(Ytf))
    print("")
    print("Low-frequency slope Delta2tf =")
    sp.pprint(sp.factor(Delta2tf))
    print("")
    print("One-pole Pade lambda_eff(TF) =")
    sp.pprint(sp.factor(lamEffTF))

    section("4) EM-worked point with TF inertia")
    print("Lam_EM =", LamEM)
    print("Sigma* =", SigmaStar)
    print("")
    print("mHatTF(worked point) =")
    print(mNum)
    print("")
    print("Dimensionless poles lambda_i of mHatTF^{-1} hBar =")
    print(evals)
    print("")
    print("Positive residues R_i =")
    print(residues)
    print("Residue fractions =")
    print(residueFracs)
    print("")
    print("Pade lambda_eff(TF) =", lamEffTFNum)
    print("Max relative error on 0 <= s <= 0.1 lambda_- =", relErr)
    print("")
    print("Worked-point rhoEffTF(Lam_EM, a0=1) =")
    sp.pprint(sp.factor(rhoEffWorked))
    print("")
    print("Physical pole scales after inserting rhoEffTF:")
    print("Omega_-^2 =")
    sp.pprint(sp.factor(omegaMinusSq))
    print("")
    print("Omega_+^2 =")
    sp.pprint(sp.factor(omegaPlusSq))
    print("")
    print("Omega_eff^2 =")
    sp.pprint(sp.factor(omegaEffSq))

    section("5) Interpretation")
    print("1) The Family-1 parent PDE fixes the free bulk inertia scale:")
    print("      rhoEff -> rhoEffTF(n=5, OmegaIn, KEOS, L0, mPsi).")
    print("")
    print("2) The same TF reduction fixes the axial breathing moment exactly:")
    print("      1/12  ->  1/14 for the frozen n = 5 EOS.")
    print("")
    print("3) With the carried-forward static Hessian, the monopole wall channel remains")
    print("   an exact two-pole Stieltjes response with positive residues, and the")
    print("   single-pole reduction stays highly accurate.")
    print("")
    print("4) What is still open is the genuine soft-wall / surface inertia completion,")
    print("   not the bulk inertia scale itself.")


if __name__ == "__main__":
    main()
