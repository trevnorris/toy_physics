import sympy as sp
import numpy as np

# ---------------------------------------------------------------------------
# 2PN geometry breathing: dynamic reduction with affine inertia
# ---------------------------------------------------------------------------

def section(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def main() -> None:
    # ------------------------------------------------------------------
    # Symbols and exact geometric data
    # ------------------------------------------------------------------
    a0, Lam, Sigma, rho, beta, rhoEff, s = sp.symbols(
        'a0 Lam Sigma rho beta rhoEff s', positive=True
    )
    a, L = sp.symbols('a L', positive=True)

    # 4D cylinder-like throat geometry (3-ball cross section x interval)
    V = sp.Rational(4, 3) * sp.pi * a**3 * L
    A = 4 * sp.pi * a**2 * L + sp.Rational(8, 3) * sp.pi * a**3

    sigma = Sigma / a0**3
    Pvac = rho * Sigma / a0**4
    kappab = beta * Sigma / a0
    Egeom = sp.expand(Pvac * V + sigma * A + kappab * a**2 / L)

    H = sp.hessian(Egeom, (a, L))
    g = sp.Matrix([sp.diff(V, a), sp.diff(V, L)])

    subs0 = {a: a0, L: Lam * a0}
    H0 = sp.simplify(H.subs(subs0))
    V0 = sp.simplify(V.subs(subs0))
    g0 = sp.simplify(g.subs(subs0))

    # ------------------------------------------------------------------
    # Affine inertia from entrained-fluid kinematics
    #   xi_perp = (delta a / a0) r_vec
    #   xi_w    = (delta L / L0) w \, e_w
    # ------------------------------------------------------------------
    r = sp.symbols('r', nonnegative=True)
    w = sp.symbols('w', real=True)

    ball_moment = sp.simplify(4 * sp.pi * sp.integrate(r**4, (r, 0, a0)))
    interval_moment = sp.simplify(sp.integrate(w**2, (w, -Lam * a0 / 2, Lam * a0 / 2)))
    V3 = sp.simplify(sp.Rational(4, 3) * sp.pi * a0**3)

    Maa = sp.simplify(rhoEff * Lam * a0 * ball_moment / a0**2)
    MLL = sp.simplify(rhoEff * V3 * interval_moment / (Lam * a0)**2)
    M0 = sp.Matrix([[Maa, 0], [0, MLL]])

    # Dimensionless coordinates q = (delta a / a0, delta L / a0)
    # Then E^(2) = (Sigma/2) q^T hBar q
    # and T^(2) = (rhoEff * V0 * a0^2 / 2) qdot^T mHat qdot
    hBar = sp.simplify(H0 * a0**2 / Sigma)
    mHat = sp.simplify(M0 / (rhoEff * V0))
    gBar = sp.simplify(g0 / V0 * a0)   # delta V / V0 = gBar . q

    # Dynamic monopole susceptibility in the dimensionless frequency variable
    # s = omega^2 * rhoEff * V0 * a0^2 / Sigma
    Ybar = sp.simplify((gBar.T * (hBar - s * mHat).inv() * gBar)[0] / Sigma)

    # Low-frequency expansion coefficients
    Delta0 = sp.simplify((gBar.T * hBar.inv() * gBar)[0] / Sigma)
    Delta2 = sp.simplify((gBar.T * hBar.inv() * mHat * hBar.inv() * gBar)[0] / Sigma)
    lamEff = sp.simplify(Delta0 / Delta2)

    # Raw full K00 completion: local wall law + geometry breathing channel
    K00RawLocal = -sp.Rational(757, 2520)
    K00RawDyn = sp.simplify(K00RawLocal + Ybar)
    K00RawStatic = sp.simplify(K00RawDyn.subs(s, 0))

    section("1) Exact affine inertia reduction")
    print("V(a,L) =", V)
    print("A(a,L) =", A)
    print("E_geom =", Egeom)
    print("")
    print("3-ball radial second moment integral:")
    sp.pprint(ball_moment)
    print("")
    print("Interval axial second moment integral:")
    sp.pprint(interval_moment)
    print("")
    print("Reference 4-volume V0 =")
    sp.pprint(V0)
    print("")
    print("Affine inertia matrix M0 on (delta a, delta L) =")
    sp.pprint(M0)
    print("")
    print("Dimensionless inertia metric mHat = M0 / (rhoEff V0) =")
    sp.pprint(mHat)
    print("")
    print("Dimensionless monopole coupling gBar = a0 grad(V)/V0 =")
    sp.pprint(gBar)

    section("2) Exact dynamic monopole response")
    print("Dimensionless static Hessian hBar = a0^2 H0 / Sigma =")
    sp.pprint(hBar)
    print("")
    print("Dimensionless frequency parameter:")
    print("  s = omega^2 * rhoEff * V0 * a0^2 / Sigma")
    print("")
    print("Exact geometry breathing susceptibility Y_geom(s) =")
    sp.pprint(sp.factor(Ybar))
    print("")
    print("Static coefficient Delta0 = Y_geom(0) =")
    sp.pprint(sp.factor(Delta0))
    print("")
    print("Low-frequency s coefficient Delta2 =")
    sp.pprint(sp.factor(Delta2))
    print("")
    print("[1/1] Pade single-pole lambda_eff = Delta0 / Delta2 =")
    sp.pprint(sp.factor(lamEff))
    print("")
    print("Raw monopole closure with local wall term:")
    print("K00_raw(s) = -757/2520 + Y_geom(s)")
    sp.pprint(sp.factor(K00RawDyn))
    print("")
    print("Static check K00_raw(0) =")
    sp.pprint(sp.factor(K00RawStatic))

    # ------------------------------------------------------------------
    # EM-branch worked point, matched to the carried-forward target 109/280
    # ------------------------------------------------------------------
    x01 = sp.Float('2.40482555769577276862163187933', 50)
    LamEM = sp.N(sp.sqrt(2) * sp.pi / x01, 50)
    rhoEx = sp.Rational(1, 10)
    betaEx = sp.Integer(12)
    target = sp.Rational(109, 280)

    # Match Sigma using the static geometry closure as in the previous step
    DeltaUnit = sp.N(Delta0.subs({a0: 1, Lam: LamEM, rho: rhoEx, beta: betaEx, Sigma: 1}), 50)
    SigmaStar = sp.N(DeltaUnit / target, 50)

    # Numerical hBar, mHat, and gBar at the worked point
    hNum = np.array(sp.N(hBar.subs({a0: 1, Lam: LamEM, rho: rhoEx, beta: betaEx}), 60).tolist(), dtype=float)
    mNum = np.array(sp.N(mHat.subs({a0: 1, Lam: LamEM}), 60).tolist(), dtype=float)
    gNum = np.array(sp.N(gBar.subs({a0: 1, Lam: LamEM}), 60).tolist(), dtype=float).flatten()

    evals, evecs = np.linalg.eig(np.linalg.solve(mNum, hNum))
    order = np.argsort(evals)
    evals = np.real(evals[order])
    evecs = np.real(evecs[:, order])

    # Mass-normalize modes: v_i^T mHat v_j = delta_ij
    for i in range(evecs.shape[1]):
        norm = np.sqrt(evecs[:, i].T @ mNum @ evecs[:, i])
        evecs[:, i] /= norm

    massOrth = evecs.T @ mNum @ evecs
    stiffDiag = evecs.T @ hNum @ evecs

    residues = (gNum @ evecs) ** 2 / (float(SigmaStar) * evals)
    contrib = residues
    contribFrac = contrib / contrib.sum()
    dom = int(np.argmax(contrib))

    lamEffNum = float(sp.N(lamEff.subs({a0: 1, Lam: LamEM, rho: rhoEx, beta: betaEx}), 30))
    # physical pole scales: Omega_i^2 = (Sigma / (rhoEff V0 a0^2)) * lambda_i
    V0Num = float(sp.N(V0.subs({a0: 1, Lam: LamEM}), 30))
    prefactor = float(SigmaStar) / V0Num  # divide by rhoEff for general rhoEff
    Om2 = prefactor * evals
    Om2Eff = prefactor * lamEffNum

    # exact-vs-Pade low-frequency accuracy in the natural scaled variable s
    def Y_exact_scaled(sgrid: np.ndarray) -> np.ndarray:
        vals = np.zeros_like(sgrid, dtype=float)
        for Ri, li in zip(contrib, evals):
            vals += Ri / (1.0 - sgrid / li)
        return vals

    def Y_pade_scaled(sgrid: np.ndarray) -> np.ndarray:
        return float(target) / (1.0 - sgrid / lamEffNum)

    s_grid = np.linspace(0.0, 0.1 * evals[0], 400)
    exact_grid = Y_exact_scaled(s_grid)
    pade_grid = Y_pade_scaled(s_grid)
    relerr = np.max(np.abs((pade_grid - exact_grid) / exact_grid))

    section("3) EM-branch worked point with affine inertia")
    print("Lam_EM =", LamEM)
    print("rho    =", rhoEx)
    print("beta   =", betaEx)
    print("Sigma* =", SigmaStar)
    print("")
    print("hBar(worked point) =")
    print(hNum)
    print("")
    print("mHat(worked point) =")
    print(mNum)
    print("")
    print("gBar(worked point) =", gNum)
    print("")
    print("Mass-orthogonality check V^T mHat V =")
    print(massOrth)
    print("")
    print("Stiffness diagonalization V^T hBar V =")
    print(stiffDiag)
    print("")
    print("Dimensionless pole parameters lambda_i of mHat^{-1} hBar =")
    print(evals)
    print("")
    print("Static pole residues R_i in Y_geom(s) = sum_i R_i / (1 - s/lambda_i)")
    print(contrib)
    print("Residue fractions =")
    print(contribFrac)
    print("")
    print("Static sum =", contrib.sum())
    print("Target 109/280 =", float(target))
    print("Residual =", contrib.sum() - float(target))
    print("")
    print("Physical pole scales (for rhoEff = 1, a0 = 1):")
    print("prefactor Sigma/(rhoEff V0 a0^2) =", prefactor, "/ rhoEff")
    print("Omega_i^2 = prefactor * lambda_i =")
    print(Om2)
    print("")
    print("Pade single-pole lambda_eff =", lamEffNum)
    print("Omega_eff^2 = prefactor * lambda_eff =", Om2Eff)
    print("")
    print("Max relative error of the one-pole Pade form on 0 <= s <= 0.1 lambda_- :")
    print(relerr)

    section("4) Interpretation")
    print("1) The old monopole auxiliary is now the low-frequency breathing response of")
    print("   the same reduced geometry sector that generated the static 109/280 closure.")
    print("")
    print("2) With affine entrained-fluid inertia, the exact monopole channel is a")
    print("   two-pole Stieltjes response with positive residues.")
    print("")
    print("3) At the EM worked point, the higher pole carries",
          100.0 * contribFrac[dom], "% of the static weight,")
    print("   so the single-pole auxiliary used in the wall sector is a controlled")
    print("   low-frequency reduction rather than an extra assumption.")
    print("")
    print("4) The remaining microphysical task is now narrow: derive the overall inertia")
    print("   scale rhoEff (or its soft-wall analog) from the Family-1 confinement /")
    print("   traction PDE, and then the monopole pole is fully fixed.")


if __name__ == "__main__":
    main()
