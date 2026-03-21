
import sympy as sp
import numpy as np

# ---------------------------------------------------------------------------
# 2PN geometry-Hessian / monopole-breathing closure
# ---------------------------------------------------------------------------

def section(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)

def main() -> None:
    # Symbols
    a0, Lam, Sigma, rho, beta = sp.symbols('a0 Lam Sigma rho beta', positive=True)
    a, L = sp.symbols('a L', positive=True)

    # 4D cylinder throat geometry: 3-ball cross section of radius a, length L
    V = sp.Rational(4, 3) * sp.pi * a**3 * L
    A = 4 * sp.pi * a**2 * L + sp.Rational(8, 3) * sp.pi * a**3

    # Minimal geometry energy completion:
    #   E_geom = P_vac V + sigma A + kappa_b a^2/L
    sigma = Sigma / a0**3
    Pvac = rho * Sigma / a0**4
    kappab = beta * Sigma / a0

    E = sp.expand(Pvac * V + sigma * A + kappab * a**2 / L)

    x = sp.Matrix([a, L])
    H = sp.hessian(E, (a, L))
    g = sp.Matrix([sp.diff(V, a), sp.diff(V, L)])

    subs0 = {a: a0, L: Lam * a0}
    H0 = sp.simplify(H.subs(subs0))
    g0 = sp.simplify(g.subs(subs0))
    V0 = sp.simplify(V.subs(subs0))

    # Natural geometry-side compressibility / projector coefficient
    DeltaGeom = sp.factor(sp.simplify((g0.T * H0.inv() * g0)[0] / V0**2))
    Hbar = sp.simplify(H0 / (Sigma / a0**2))
    detHbar = sp.factor(sp.det(Hbar))

    # Baseline no-go: set kappa_b -> 0
    HbarBase = sp.simplify(Hbar.subs({beta: 0}))
    detHbarBase = sp.factor(detHbar.subs({beta: 0}))

    # Exact positivity / positive-support thresholds
    betaStab = sp.simplify(sp.pi * Lam**3 * (rho + 2)**2 / (2 * Lam * rho + 3 * Lam + 2))
    betaDelta = sp.simplify(sp.pi * Lam * (2 * Lam * rho + 5 * Lam - 2) / 4)

    section("1) Exact geometry-side setup")
    print("V(a,L) =", V)
    print("A(a,L) =", A)
    print("E_geom =", E)
    print("")
    print("Reference-point Hessian H0 =")
    sp.pprint(H0)
    print("")
    print("Natural volume-work coupling vector g0 = grad V |_(a0, L0) =")
    sp.pprint(g0)
    print("")
    print("Natural dimensionless monopole projector coefficient:")
    print("Delta_geom = (grad V)^T H0^{-1} (grad V) / V0^2 =")
    sp.pprint(DeltaGeom)

    section("2) Baseline no-go: P_vac V + sigma A alone")
    print("Hbar(beta=0) =")
    sp.pprint(HbarBase)
    print("")
    print("det Hbar(beta=0) =")
    sp.pprint(detHbarBase)
    print("")
    print("So the baseline P_vac V + sigma A model is not a passive 2DOF geometry Hessian.")
    print("It cannot by itself generate the monopole breathing auxiliary.")

    section("3) Minimal curvature completion and exact support formula")
    print("Hbar = H0 / (Sigma / a0^2) =")
    sp.pprint(Hbar)
    print("")
    print("det Hbar =")
    sp.pprint(detHbar)
    print("")
    print("Positive-definite Hbar requires:")
    print("  beta > 0")
    print("  beta > beta_stab(Lam, rho) with")
    sp.pprint(betaStab)
    print("")
    print("Positive Delta_geom additionally requires:")
    print("  beta > beta_Delta(Lam, rho) with")
    sp.pprint(betaDelta)
    print("")
    print("So the minimal curvature-completed geometry model closes the monopole channel")
    print("whenever beta > max(beta_stab, beta_Delta).")

    # Concrete EM-cavity-limit example
    x01 = sp.Float('2.40482555769577276862163187933', 40)
    LamEM = sp.N(sp.sqrt(2) * sp.pi / x01, 40)
    rhoEx = sp.Rational(1, 10)
    betaEx = sp.Integer(12)
    target = sp.Rational(109, 280)

    DeltaUnit = sp.simplify(DeltaGeom.subs({a0: 1, Sigma: 1, Lam: LamEM, rho: rhoEx, beta: betaEx}))
    SigmaStar = sp.N(DeltaUnit / target, 30)
    DeltaMatched = sp.N(DeltaGeom.subs({a0: 1, Sigma: SigmaStar, Lam: LamEM, rho: rhoEx, beta: betaEx}), 30)

    Hnum = np.array(sp.N(Hbar.subs({a0: 1, Lam: LamEM, rho: rhoEx, beta: betaEx}) * SigmaStar, 50).tolist(), dtype=float)
    ghat = np.array(sp.N((g0 / V0).subs({a0: 1, Lam: LamEM}), 50).tolist(), dtype=float).flatten()

    evals, evecs = np.linalg.eigh(Hnum)
    proj = evecs.T @ ghat
    contrib = proj**2 / evals
    frac = contrib / contrib.sum()
    dominant = int(np.argmax(contrib))
    relerr_onepole = abs(float(target) - float(contrib[dominant])) / float(target)

    section("4) EM-branch worked example")
    print("x01    =", x01)
    print("Lam_EM =", LamEM)
    print("rho    =", rhoEx)
    print("beta   =", betaEx)
    print("")
    print("beta_stab(Lam_EM, rho)  =", sp.N(betaStab.subs({Lam: LamEM, rho: rhoEx}), 20))
    print("beta_Delta(Lam_EM, rho) =", sp.N(betaDelta.subs({Lam: LamEM, rho: rhoEx}), 20))
    print("")
    print("Sigma_* needed to match Delta_geom = 109/280:")
    print("Sigma_* =", SigmaStar)
    print("")
    print("Delta_geom(Sigma_*) =", DeltaMatched)
    print("Target 109/280      =", sp.N(target, 30))
    print("Residual            =", sp.N(DeltaMatched - target, 30))
    print("")
    print("Hessian eigenvalues at the worked point:")
    print(evals)
    print("")
    print("Coupling vector ghat = grad V / V0 =")
    print(ghat)
    print("")
    print("Mode-resolved contributions to Delta_geom:")
    print(contrib)
    print("Contribution fractions:")
    print(frac)
    print("")
    print("Dominant-mode relative one-pole error:")
    print(relerr_onepole)

    section("5) Interpretation")
    print("1) The previous global monopole projector is not ad hoc.")
    print("   It is the exact low-frequency compressibility obtained by integrating out")
    print("   the reduced geometry DOFs (a,L) with natural volume-work coupling.")
    print("")
    print("2) The baseline P_vac V + sigma A model fails exactly as the notes warned:")
    print("   without the curvature/bend completion, the geometry Hessian is not passive.")
    print("")
    print("3) In the EM-limit worked example, one eigenmode carries",
          100.0 * frac[dominant], "% of the total static monopole response,")
    print("   so the earlier single breathing auxiliary is a controlled reduction of the")
    print("   full 2DOF geometry sector rather than a separate new assumption.")

if __name__ == "__main__":
    main()
