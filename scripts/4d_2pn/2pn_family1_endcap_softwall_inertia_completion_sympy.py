import sympy as sp
import numpy as np
import mpmath as mp

# ---------------------------------------------------------------------------
# 2PN Family-1 endcap soft-wall inertia completion
# ---------------------------------------------------------------------------
# This step derives the filled-to-endcap soft-cap correction and then folds it
# into the carried-forward monopole breathing response.  The same script also
# shows the leading separated-order composite branch in which the new endcap
# result is combined with the previously derived sidewall branch.
# ---------------------------------------------------------------------------


def section(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


# ---------------------------------------------------------------------------
# Logistic wall profile and numerical helpers
# ---------------------------------------------------------------------------
mp.mp.dps = 80


def S(x: mp.mpf) -> mp.mpf:
    return mp.mpf('0.5') * (1 + mp.tanh(x))


def endcap_turning_point(alpha: mp.mpf, p: int) -> mp.mpf:
    # Solve x + alpha S(x)^p = 0 for the local reduced cap profile.
    f = lambda x: x + alpha * (S(x) ** p)
    a = -max(mp.mpf('20'), 10 * alpha + 5)
    b = mp.mpf('2')
    # Bracket scan for safety.
    xs = [a + (b - a) * i / 200 for i in range(201)]
    vals = [f(x) for x in xs]
    lo, hi = None, None
    for i in range(len(xs) - 1):
        if vals[i] == 0:
            return xs[i]
        if vals[i] * vals[i + 1] < 0:
            lo, hi = xs[i], xs[i + 1]
            break
    if lo is None:
        lo, hi = -alpha - 1, mp.mpf('0.5')
    return mp.findroot(f, (lo, hi))


# Local reduced profile for the filled-to-endcap branch.
# The thin-endcap scaling is W_cap / mu_c = 2 eps_z alpha_z, so the local
# profile is controlled by (-x - alpha_z S(x)^p)_+^(1/4).

def g_local(x: mp.mpf, alpha: mp.mpf, p: int) -> mp.mpf:
    arg = -x - alpha * (S(x) ** p)
    return arg ** mp.mpf('0.25') if arg > 0 else mp.mpf('0')


def g_sharp(x: mp.mpf) -> mp.mpf:
    return (-x) ** mp.mpf('0.25') if x < 0 else mp.mpf('0')


def nu_k(alpha: mp.mpf, p: int, k: int) -> mp.mpf:
    xt = endcap_turning_point(alpha, p)

    def integrand(x: mp.mpf) -> mp.mpf:
        return (x ** k) * (g_local(x, alpha, p) - g_sharp(x))

    return mp.quad(integrand, [-mp.inf, xt, 0, mp.inf])


# Full dimensionless cap profile in u = 2 w / L.
# The soft endcap enters as 2 eps_z alpha_z S((|u|-1)/eps_z)^p.

def f_full(u: mp.mpf, eps_z: mp.mpf, alpha_z: mp.mpf, p_z: int) -> mp.mpf:
    x = (abs(u) - 1) / eps_z
    arg = 1 - u * u - 2 * eps_z * alpha_z * (S(x) ** p_z)
    return arg ** mp.mpf('0.25') if arg > 0 else mp.mpf('0')


def c0_full(eps_z: mp.mpf, alpha_z: mp.mpf, p_z: int) -> mp.mpf:
    f = lambda u: f_full(u, eps_z, alpha_z, p_z)
    return mp.mpf('0.5') * mp.quad(f, [-2, -1, 0, 1, 2])


def c2_full(eps_z: mp.mpf, alpha_z: mp.mpf, p_z: int) -> mp.mpf:
    f = lambda u: u * u * f_full(u, eps_z, alpha_z, p_z)
    return mp.mpf('0.5') * mp.quad(f, [-2, -1, 0, 1, 2])


# ---------------------------------------------------------------------------
# Geometry / breathing response carried forward from earlier steps
# ---------------------------------------------------------------------------
def geometry_response_arrays() -> tuple[np.ndarray, np.ndarray, float, float, float]:
    a0, Lam, rho, beta, Sigma = sp.symbols('a0 Lam rho beta Sigma', positive=True)
    a, L = sp.symbols('a L', positive=True)

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
    hBar = sp.simplify(H0 * a0**2 / Sigma)
    gBar = sp.simplify(g0 / V0 * a0)
    Delta0 = sp.simplify((gBar.T * hBar.inv() * gBar)[0] / Sigma)

    x01 = sp.Float('2.40482555769577276862163187933', 50)
    LamEM = sp.N(sp.sqrt(2) * sp.pi / x01, 60)
    rhoEx = sp.Rational(1, 10)
    betaEx = sp.Integer(12)
    target = sp.Rational(109, 280)

    DeltaUnit = sp.N(Delta0.subs({a0: 1, Lam: LamEM, rho: rhoEx, beta: betaEx, Sigma: 1}), 80)
    SigmaStar = float(sp.N(DeltaUnit / target, 40))
    hNum = np.array(sp.N(hBar.subs({a0: 1, Lam: LamEM, rho: rhoEx, beta: betaEx}), 60).tolist(), dtype=float)
    gNum = np.array(sp.N(gBar.subs({a0: 1, Lam: LamEM}), 60).tolist(), dtype=float).flatten()
    V0Num = float(sp.N(V0.subs({a0: 1, Lam: LamEM}), 50))
    return hNum, gNum, SigmaStar, float(LamEM), V0Num


def pole_data(hNum: np.ndarray, gNum: np.ndarray, SigmaStar: float, mHat: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    evals, evecs = np.linalg.eig(np.linalg.solve(mHat, hNum))
    order = np.argsort(evals)
    evals = np.real(evals[order])
    evecs = np.real(evecs[:, order])
    for i in range(evecs.shape[1]):
        norm = np.sqrt(evecs[:, i].T @ mHat @ evecs[:, i])
        evecs[:, i] /= norm
    residues = (gNum @ evecs) ** 2 / (SigmaStar * evals)
    lam_eff = residues.sum() / np.sum(residues / evals)
    return evals, residues, lam_eff


def max_relerr(evals: np.ndarray, residues: np.ndarray, lam_eff: float) -> float:
    s_grid = np.linspace(0.0, 0.1 * evals[0], 400)
    exact = np.zeros_like(s_grid)
    for Ri, li in zip(residues, evals):
        exact += Ri / (1.0 - s_grid / li)
    pade = residues.sum() / (1.0 - s_grid / lam_eff)
    return float(np.max(np.abs((pade - exact) / exact)))


def main() -> None:
    # Exact TF baseline on the filled-to-endcap n = 5 branch.
    c0 = sp.sqrt(sp.pi) * sp.gamma(sp.Rational(5, 4)) / (2 * sp.gamma(sp.Rational(7, 4)))
    c2 = sp.Rational(2, 7) * c0
    c0_num = mp.mpf(str(sp.N(c0, 50)))
    c2_num = mp.mpf(str(sp.N(c2, 50)))

    section("1) Filled-to-endcap soft-cap scaling")
    eps_z, alpha_z = sp.symbols('eps_z alpha_z', positive=True)
    print("Baseline n = 5 filled-to-endcap TF weight:")
    print("  rho_0(u) \propto (1 - u^2)^(1/4),   u = 2 w / L")
    print("")
    print("Near the cap, u = 1 + eps_z x, so")
    print("  1 - u^2 = eps_z (-2 x - eps_z x^2).")
    print("")
    print("Therefore a genuine thin-endcap layer must scale as")
    print("  V_cap / mu_c = 2 eps_z alpha_z S(x)^p,")
    print("not O(1) in mu_c.  This is the key difference from the sidewall.")
    print("")
    print("Reduced local profile:")
    print("  g_{alpha,p}(x) = (-x - alpha S(x)^p)_+^(1/4)")
    print("")
    print("So the endcap correction is weaker than the sidewall by one extra")
    print("power of the vanishing TF profile at the cap, namely eps_z^(5/4).")

    section("2) Exact asymptotic form of the endcap correction")
    print("Define the defect moments")
    print("  nu_k(alpha,p) = Integral x^k [g_{alpha,p}(x) - (-x)_+^(1/4)] dx")
    print("")
    print("Then for n = 5:")
    print("  c0^cap = c0 + 2^(1/4) eps_z^(5/4) nu_0 + O(eps_z^(9/4))")
    print("  c2^cap = c2 + 2^(1/4) eps_z^(5/4) nu_0 + O(eps_z^(9/4))")
    print("")
    Acoeff = sp.simplify(sp.root(2, 4) / c0)
    Bcoeff = sp.simplify(5 * sp.root(2, 4) / (28 * c0))
    print("Therefore")
    print("  rho_eff^(TF+cap) / rho_eff^TF = 1 + A_cap nu_0 eps_z^(5/4) + ...")
    print("  M_LL^(TF+cap) = 1/14 + B_cap nu_0 eps_z^(5/4) + ...")
    print("")
    print("A_cap =")
    sp.pprint(Acoeff)
    print("")
    print("B_cap =")
    sp.pprint(Bcoeff)
    print("")
    print("Numerically:")
    print("  A_cap =", float(sp.N(Acoeff, 20)))
    print("  B_cap =", float(sp.N(Bcoeff, 20)))

    # Representative branch and direct verification.
    eps_rep = mp.mpf('0.05')
    alpha_rep = mp.mpf('1')
    p_rep = 2
    nu0_rep = nu_k(alpha_rep, p_rep, 0)
    xstar_rep = endcap_turning_point(alpha_rep, p_rep)

    dc_rep = (mp.power(2, mp.mpf('0.25'))) * (eps_rep ** mp.mpf('1.25')) * nu0_rep
    c0_as = c0_num + dc_rep
    c2_as = c2_num + dc_rep
    Mll_as = c2_as / (4 * c0_as)

    c0_ex = c0_full(eps_rep, alpha_rep, p_rep)
    c2_ex = c2_full(eps_rep, alpha_rep, p_rep)
    Mll_ex = c2_ex / (4 * c0_ex)

    section("3) Representative steep-cap branch and direct full-profile check")
    print("Representative endcap branch:")
    print("  eps_z =", eps_rep)
    print("  alpha_z =", alpha_rep)
    print("  p_z =", p_rep)
    print("")
    print("Local turning point x_* solving x + alpha S(x)^p = 0:")
    print("  x_* =", xstar_rep)
    print("")
    print("Defect moment nu_0 =", nu0_rep)
    print("")
    print("Asymptotic vs direct full-profile values:")
    print("  c0 asymptotic =", c0_as)
    print("  c0 direct     =", c0_ex)
    print("  relative error =", (c0_as - c0_ex) / c0_ex)
    print("")
    print("  c2 asymptotic =", c2_as)
    print("  c2 direct     =", c2_ex)
    print("  relative error =", (c2_as - c2_ex) / c2_ex)
    print("")
    print("  M_LL asymptotic =", Mll_as)
    print("  M_LL direct     =", Mll_ex)
    print("  relative error  =", (Mll_as - Mll_ex) / Mll_ex)
    print("")
    print("So the leading eps_z^(5/4) asymptotic is already very accurate at eps_z = 0.05.")

    # Carried-forward geometry response.
    hNum, gNum, SigmaStar, LamEM, V0Num = geometry_response_arrays()
    target = float(sp.N(sp.Rational(109, 280), 30))
    base_mHat = np.array([[3 / 5, 0], [0, 1 / 14]], dtype=float)
    base_evals, base_res, base_lam_eff = pole_data(hNum, gNum, SigmaStar, base_mHat)
    base_relerr = max_relerr(base_evals, base_res, base_lam_eff)

    Rcap = float(c0_ex / c0_num)
    mHat_cap = np.array([[3 / 5, 0], [0, float(Mll_ex)]], dtype=float)
    cap_evals, cap_res, cap_lam_eff = pole_data(hNum, gNum, SigmaStar, mHat_cap)
    cap_relerr = max_relerr(cap_evals, cap_res, cap_lam_eff)
    Om2_ratio_cap = (cap_evals / Rcap) / base_evals

    section("4) Dynamic monopole response with the endcap correction")
    print("Baseline TF branch:")
    print("  lambda_- =", base_evals[0])
    print("  lambda_+ =", base_evals[1])
    print("  residues =", base_res)
    print("  lambda_eff =", base_lam_eff)
    print("  max relative Pade error =", base_relerr)
    print("")
    print("Endcap-corrected branch:")
    print("  rho_eff factor R_cap = c0^cap / c0 =", Rcap)
    print("  mHat_cap =")
    print(mHat_cap)
    print("")
    print("  lambda_- =", cap_evals[0])
    print("  lambda_+ =", cap_evals[1])
    print("  residues =", cap_res)
    print("  residue fractions =", cap_res / cap_res.sum())
    print("  lambda_eff =", cap_lam_eff)
    print("  max relative Pade error =", cap_relerr)
    print("")
    print("Physical pole-squared ratios relative to the sharp-wall TF baseline:")
    print("  Omega_-^2 / Omega_-^2(sharp) =", Om2_ratio_cap[0])
    print("  Omega_+^2 / Omega_+^2(sharp) =", Om2_ratio_cap[1])
    print("")
    print("Static sum of residues =", cap_res.sum(), "(target 109/280 =", target, ")")

    # ------------------------------------------------------------------
    # Leading separated-order composite branch: combine the carried-forward
    # sidewall result with the new endcap result.
    # ------------------------------------------------------------------
    Rside = 0.9060975247692787
    Maa_side = 0.5623811549096673

    Rfull = Rside * Rcap
    mHat_full = np.array([[Maa_side, 0], [0, float(Mll_ex)]], dtype=float)
    full_evals, full_res, full_lam_eff = pole_data(hNum, gNum, SigmaStar, mHat_full)
    full_relerr = max_relerr(full_evals, full_res, full_lam_eff)
    Om2_ratio_full = (full_evals / Rfull) / base_evals

    section("5) Leading separated-order full-wall composite branch")
    print("Carried-forward sidewall branch (from the previous step):")
    print("  R_side  =", Rside)
    print("  M_aa    =", Maa_side)
    print("")
    print("New endcap branch:")
    print("  R_cap   =", Rcap)
    print("  M_LL    =", float(Mll_ex))
    print("")
    print("Leading separated-order composite:")
    print("  R_full = R_side * R_cap =", Rfull)
    print("  mHat_full =")
    print(mHat_full)
    print("")
    print("  lambda_- =", full_evals[0])
    print("  lambda_+ =", full_evals[1])
    print("  residues =", full_res)
    print("  residue fractions =", full_res / full_res.sum())
    print("  lambda_eff =", full_lam_eff)
    print("  max relative Pade error =", full_relerr)
    print("")
    print("Physical pole-squared ratios relative to the sharp-wall TF baseline:")
    print("  Omega_-^2 / Omega_-^2(sharp) =", Om2_ratio_full[0])
    print("  Omega_+^2 / Omega_+^2(sharp) =", Om2_ratio_full[1])
    print("")
    print("This is the first full wall-completed monopole breathing branch in the")
    print("current Family-1 program, at least to separated leading order in the")
    print("sidewall and endcap thickness parameters.")

    section("6) Interpretation")
    print("1) The endcap layer is parametrically weaker than the sidewall because the")
    print("   filled-to-endcap TF profile already vanishes at the cap.")
    print("")
    print("2) The correct cap scaling is W_cap / mu_c = O(eps_z), and the first")
    print("   nontrivial correction is O(eps_z^(5/4)) on the n = 5 branch.")
    print("")
    print("3) Even after adding the cap correction, the monopole wall channel remains")
    print("   an excellent positive two-pole Stieltjes response with a one-pole")
    print("   reduction that is accurate well below the lower pole.")
    print("")
    print("4) Combining the carried-forward sidewall branch with the new cap branch")
    print("   gives a near-final full-wall dynamic monopole response.  The remaining")
    print("   gap is the fully coupled sidewall-cap derivation beyond separated order.")


if __name__ == "__main__":
    main()
