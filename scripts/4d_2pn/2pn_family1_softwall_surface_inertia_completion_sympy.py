
import math
from functools import lru_cache
from typing import Dict, List, Tuple

import numpy as np
import sympy as sp
from scipy import integrate
from scipy.special import roots_legendre


def section(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


# -----------------------------------------------------------------------------
# Family-1 soft-wall radial boundary layer
# -----------------------------------------------------------------------------

def smooth_step(x: float) -> float:
    return 0.5 * (1.0 + math.tanh(x))


def xi_star(alpha: float, p: float) -> float:
    """
    Turning point of the n=5 TF wall profile:
        1 - alpha * S(x)^p = 0  =>  S(x) = alpha^(-1/p).
    """
    return math.atanh(2.0 * (alpha ** (-1.0 / p)) - 1.0)


def wall_profile(x: float, alpha: float, p: float) -> float:
    """
    Local n=5 TF wall profile at fixed axial slice:
        f_{alpha,p}(x) = (1 - alpha S(x)^p)_+^{1/4}.
    """
    val = 1.0 - alpha * (smooth_step(x) ** p)
    return max(val, 0.0) ** 0.25


def defect_integrand(x: float, alpha: float, p: float, k: int) -> float:
    """
    Boundary-layer defect relative to a sharp step.
    """
    inside_sharp = 1.0 if x < 0.0 else 0.0
    return (x ** k) * (wall_profile(x, alpha, p) - inside_sharp)


@lru_cache(maxsize=None)
def wall_moment(alpha: float, p: float, k: int) -> float:
    """
    m_k(alpha,p) = ∫ x^k [f_{alpha,p}(x) - Theta(-x)] dx
    """
    xs = xi_star(alpha, p)
    f = lambda xx: defect_integrand(xx, alpha, p, k)

    if xs >= 0.0:
        res1 = integrate.quad(f, -np.inf, 0.0, limit=300)[0]
        res2 = integrate.quad(f, 0.0, xs, limit=300)[0]
        return float(res1 + res2)

    res1 = integrate.quad(f, -np.inf, xs, limit=300)[0]
    res2 = integrate.quad(f, xs, 0.0, limit=300)[0]
    return float(res1 + res2)


def sharp_c0_n5(npts: int = 120) -> float:
    """
    c0(5) = (1/2) ∫_{-1}^1 (1-y^2)^(1/4) dy
    """
    nodes, weights = roots_legendre(npts)
    return float(0.5 * np.dot(weights, (1.0 - nodes ** 2) ** 0.25))


def averaged_wall_moment(alpha0: float, p: float, k: int, npts: int = 120) -> float:
    r"""
    \bar m_k(alpha0,p) =
      [1/(2 c0)] ∫_{-1}^1 (1-y^2)^(1/4) m_k(alpha0/(1-y^2),p) dy
    where alpha0 = V0 / mu_c is the central wall-to-core ratio.
    """
    nodes, weights = roots_legendre(npts)
    c0 = sharp_c0_n5(npts)

    vals = np.array(
        [wall_moment(float(alpha0 / (1.0 - y * y)), float(p), int(k)) for y in nodes],
        dtype=float,
    )
    numer = 0.5 * np.dot(weights, (1.0 - nodes ** 2) ** 0.25 * vals)
    return float(numer / c0)


def symbolic_thinwall_series() -> Dict[str, sp.Expr]:
    eps, mb0, mb1, mb2 = sp.symbols('eps mb0 mb1 mb2', real=True)
    J2 = sp.Rational(1, 3) + eps * mb0 + 2 * eps ** 2 * mb1 + eps ** 3 * mb2
    J4 = sp.Rational(1, 5) + eps * mb0 + 4 * eps ** 2 * mb1 + 6 * eps ** 3 * mb2

    mass_factor = sp.expand(3 * J2)
    Maa = sp.series(sp.simplify(J4 / J2), eps, 0, 4).removeO()
    return {
        "eps": eps,
        "mb0": mb0,
        "mb1": mb1,
        "mb2": mb2,
        "J2": sp.expand(J2),
        "J4": sp.expand(J4),
        "mass_factor": sp.expand(mass_factor),
        "Maa": sp.expand(Maa),
    }


def worked_geometry_hbar_gbar() -> Tuple[np.ndarray, np.ndarray, float, float]:
    """
    Rebuild the carried-forward geometry Hessian and coupling vector at the
    EM worked point used in the previous monopole geometry scripts.
    """
    a0, Lam, Sigma, rhoGeom, betaGeom = sp.symbols(
        'a0 Lam Sigma rhoGeom betaGeom', positive=True
    )
    a, L = sp.symbols('a L', positive=True)

    V = sp.Rational(4, 3) * sp.pi * a ** 3 * L
    A = 4 * sp.pi * a ** 2 * L + sp.Rational(8, 3) * sp.pi * a ** 3

    sigma = Sigma / a0 ** 3
    Pvac = rhoGeom * Sigma / a0 ** 4
    kappab = betaGeom * Sigma / a0

    Egeom = sp.expand(Pvac * V + sigma * A + kappab * a ** 2 / L)
    H = sp.hessian(Egeom, (a, L))
    g = sp.Matrix([sp.diff(V, a), sp.diff(V, L)])

    subs0 = {a: a0, L: Lam * a0}
    H0 = sp.simplify(H.subs(subs0))
    V0 = sp.simplify(V.subs(subs0))
    g0 = sp.simplify(g.subs(subs0))

    hBar = sp.simplify(H0 * a0 ** 2 / Sigma)
    gBar = sp.simplify(g0 / V0 * a0)

    x01 = sp.Float('2.40482555769577276862163187933', 50)
    LamEM = float(sp.N(sp.sqrt(2) * sp.pi / x01, 50))

    hNum = np.array(
        sp.N(
            hBar.subs({a0: 1, Lam: LamEM, rhoGeom: sp.Rational(1, 10), betaGeom: 12}),
            60,
        ).tolist(),
        dtype=float,
    )
    gNum = np.array(
        sp.N(gBar.subs({a0: 1, Lam: LamEM}), 60).tolist(),
        dtype=float,
    ).flatten()

    # Static geometry scale matching carried-forward Delta0 = 109/280
    delta_unit = float(
        sp.N(
            (
                gBar.T
                * hBar.inv()
                * gBar
            )[0].subs({a0: 1, Lam: LamEM, rhoGeom: sp.Rational(1, 10), betaGeom: 12}),
            60,
        )
    )
    sigma_star = delta_unit / float(sp.Rational(109, 280))
    return hNum, gNum, LamEM, sigma_star


def eig_residues_for_metric(
    Maa: float, Mll: float, hbar: np.ndarray, gbar: np.ndarray, sigma_star: float
) -> Tuple[np.ndarray, np.ndarray]:
    M = np.diag([Maa, Mll])
    evals, evecs = np.linalg.eig(np.linalg.solve(M, hbar))
    order = np.argsort(np.real(evals))
    evals = np.real(evals[order])
    evecs = np.real(evecs[:, order])

    for i in range(evecs.shape[1]):
        norm = math.sqrt(float(evecs[:, i].T @ M @ evecs[:, i]))
        evecs[:, i] /= norm

    residues = (gbar @ evecs) ** 2 / (sigma_star * evals)
    return evals, residues


def pade_pole(evals: np.ndarray, residues: np.ndarray) -> float:
    delta0 = float(np.sum(residues))
    delta2 = float(np.sum(residues / evals))
    return delta0 / delta2


def max_pade_relative_error(evals: np.ndarray, residues: np.ndarray) -> float:
    lam_eff = pade_pole(evals, residues)
    sgrid = np.linspace(0.0, 0.1 * float(evals[0]), 400)
    exact = np.zeros_like(sgrid)
    for R, lam in zip(residues, evals):
        exact += R / (1.0 - sgrid / lam)
    pade = float(np.sum(residues)) / (1.0 - sgrid / lam_eff)
    return float(np.max(np.abs((pade - exact) / exact)))


def leading_pole_shift_coeffs(hbar: np.ndarray, gbar: np.ndarray, sigma_star: float) -> np.ndarray:
    Mll = 1.0 / 14.0
    Maa0 = 3.0 / 5.0

    def evals_only(Maa: float) -> np.ndarray:
        M = np.diag([Maa, Mll])
        evals = np.linalg.eigvals(np.linalg.solve(M, hbar))
        evals = np.real(evals)
        evals.sort()
        return evals

    lam0 = evals_only(Maa0)
    h = 1.0e-7
    dlam = (evals_only(Maa0 + h) - evals_only(Maa0 - h)) / (2.0 * h)

    # Maa = 3/5 + (6/5) eps * \bar m0 + ...
    # Rmass = 1 + 3 eps * \bar m0 + ...
    coeff = (dlam / lam0) * (6.0 / 5.0) - 3.0
    return coeff


def evaluate_case(
    p: float,
    alpha0: float,
    eps_r: float,
    hbar: np.ndarray,
    gbar: np.ndarray,
    sigma_star: float,
    baseline_evals: np.ndarray,
) -> Dict[str, float]:
    mb0 = averaged_wall_moment(alpha0, p, 0)
    mb1 = averaged_wall_moment(alpha0, p, 1)
    mb2 = averaged_wall_moment(alpha0, p, 2)

    Rmass = (
        1.0
        + 3.0 * eps_r * mb0
        + 6.0 * eps_r ** 2 * mb1
        + 3.0 * eps_r ** 3 * mb2
    )

    Maa = (
        3.0 / 5.0
        + (6.0 / 5.0) * eps_r * mb0
        + eps_r ** 2 * ((42.0 / 5.0) * mb1 - (18.0 / 5.0) * mb0 ** 2)
        + eps_r ** 3 * (
            (54.0 / 5.0) * mb0 ** 3
            - (162.0 / 5.0) * mb0 * mb1
            + (81.0 / 5.0) * mb2
        )
    )

    evals, residues = eig_residues_for_metric(Maa, 1.0 / 14.0, hbar, gbar, sigma_star)
    lam_eff = pade_pole(evals, residues)
    max_err = max_pade_relative_error(evals, residues)
    omega_sq_ratios = (evals / baseline_evals) / Rmass

    # Local sharp-wall threshold at the throat center:
    # xi_star(0) = 0 <=> alpha0 = 2^p.
    alpha_threshold = 2.0 ** p
    xistar0 = xi_star(alpha0, p)

    return {
        "p": p,
        "alpha0": alpha0,
        "eps_r": eps_r,
        "alpha_threshold": alpha_threshold,
        "xi_star0": xistar0,
        "mb0": mb0,
        "mb1": mb1,
        "mb2": mb2,
        "Rmass": Rmass,
        "Maa": Maa,
        "lam_minus": float(evals[0]),
        "lam_plus": float(evals[1]),
        "R_minus": float(residues[0]),
        "R_plus": float(residues[1]),
        "lam_eff": float(lam_eff),
        "max_err": float(max_err),
        "omega_sq_ratio_minus": float(omega_sq_ratios[0]),
        "omega_sq_ratio_plus": float(omega_sq_ratios[1]),
    }


def main() -> None:
    series = symbolic_thinwall_series()

    hbar, gbar, LamEM, sigma_star = worked_geometry_hbar_gbar()
    baseline_evals, baseline_res = eig_residues_for_metric(
        3.0 / 5.0, 1.0 / 14.0, hbar, gbar, sigma_star
    )
    baseline_lam_eff = pade_pole(baseline_evals, baseline_res)
    baseline_err = max_pade_relative_error(baseline_evals, baseline_res)
    coeff = leading_pole_shift_coeffs(hbar, gbar, sigma_star)

    section("1) Universal Family-1 thin-wall formulas")
    print("For the radial n=5 TF wall profile")
    print("  f_{alpha,p}(xi) = (1 - alpha * S(xi)^p)_+^(1/4),   S(x) = (1+tanh x)/2")
    print("define defect moments m_k(alpha,p) = ∫ xi^k [f - Theta(-xi)] dxi.")
    print("")
    print("Boundary-layer moments J2 and J4 through O(eps_r^3):")
    print("J2 =")
    sp.pprint(series["J2"])
    print("")
    print("J4 =")
    sp.pprint(series["J4"])
    print("")
    print("Cross-sectional mass factor R_mass = 3 J2 =")
    sp.pprint(series["mass_factor"])
    print("")
    print("Radial inertia coefficient M_aa = J4/J2 through O(eps_r^3) =")
    sp.pprint(series["Maa"])
    print("")
    print("Leading-order lesson:")
    print("  both the mass-scale correction and the radial inertia shift are controlled")
    print("  by the same averaged wall moment \\bar m0 at O(eps_r).")

    section("2) Axially averaged wall moments on the filled-to-endcap n=5 TF branch")
    print("At axial coordinate y = 2 w/L, the local chemical potential scales as")
    print("  mu(y) = mu_c * (1-y^2),")
    print("so the local wall ratio becomes")
    print("  alpha(y) = alpha0 / (1-y^2).")
    print("")
    print("The averaged moments are")
    print("  \\bar m_k(alpha0,p) = [1/(2 c0)] ∫_{-1}^1 (1-y^2)^(1/4) m_k(alpha0/(1-y^2),p) dy,")
    print("with c0 = (1/2) ∫_{-1}^1 (1-y^2)^(1/4) dy.")
    print("")
    print("Central turning-point threshold:")
    print("  xi_*(0)=0  <=>  alpha0 = 2^p.")
    print("So alpha0 > 2^p means the wall support already turns off inside the nominal radius a0.")

    section("3) Worked geometry point carried forward from the monopole branch")
    print("Lam_EM =", LamEM)
    print("Sigma*  =", sigma_star)
    print("Baseline sharp-wall dynamic poles (TF bulk only) =")
    print("  lambda_- =", baseline_evals[0])
    print("  lambda_+ =", baseline_evals[1])
    print("Baseline residues =")
    print("  R_- =", baseline_res[0])
    print("  R_+ =", baseline_res[1])
    print("Baseline Pade pole =", baseline_lam_eff)
    print("Baseline max relative error on 0 <= s <= 0.1 lambda_- =", baseline_err)
    print("")
    print("Leading small-eps_r physical pole-shift coefficients at the EM worked point:")
    print("  Omega_-^2 / Omega_-^2(sharp) = 1 + c_- * eps_r * \\bar m0 + O(eps_r^2)")
    print("  Omega_+^2 / Omega_+^2(sharp) = 1 + c_+ * eps_r * \\bar m0 + O(eps_r^2)")
    print("with")
    print("  c_- =", coeff[0])
    print("  c_+ =", coeff[1])
    print("Since steep walls give \\bar m0 < 0 on the physical branch, the soft wall")
    print("raises the physical monopole poles.")

    cases = [
        (2.0, 10.0, 0.05),  # main steep Family-1 reference
        (2.0, 20.0, 0.05),  # stiffer wall contrast
        (4.0, 10.0, 0.05),  # sharper power contrast
    ]

    section("4) Representative steep-wall branches")
    case_data: List[Dict[str, float]] = []
    for p, alpha0, eps_r in cases:
        data = evaluate_case(p, alpha0, eps_r, hbar, gbar, sigma_star, baseline_evals)
        case_data.append(data)

        print(f"\nCase p={p:g}, alpha0={alpha0:g}, eps_r={eps_r:g}")
        print("  threshold 2^p =", data["alpha_threshold"])
        print("  xi_*(0)      =", data["xi_star0"])
        print("  averaged moments:")
        print("    mbar0 =", data["mb0"])
        print("    mbar1 =", data["mb1"])
        print("    mbar2 =", data["mb2"])
        print("  mass factor R_mass =", data["Rmass"])
        print("  radial metric M_aa =", data["Maa"])
        print("  dynamic poles:")
        print("    lambda_- =", data["lam_minus"])
        print("    lambda_+ =", data["lam_plus"])
        print("  residues:")
        print("    R_- =", data["R_minus"])
        print("    R_+ =", data["R_plus"])
        print("  physical pole-squared ratios vs sharp-wall TF baseline:")
        print("    Omega_-^2 ratio =", data["omega_sq_ratio_minus"])
        print("    Omega_+^2 ratio =", data["omega_sq_ratio_plus"])
        print("  Pade lambda_eff =", data["lam_eff"])
        print("  max relative error =", data["max_err"])

    section("5) Interpretation")
    main_case = case_data[0]
    print("For the representative steep Family-1 branch")
    print("  (p, alpha0, eps_r) = (2, 10, 0.05),")
    print("the actual soft-wall boundary layer gives:")
    print("  R_mass  =", main_case["Rmass"])
    print("  M_aa    =", main_case["Maa"])
    print("  lambda_- =", main_case["lam_minus"])
    print("  lambda_+ =", main_case["lam_plus"])
    print("  Omega_-^2 / Omega_-^2(sharp) =", main_case["omega_sq_ratio_minus"])
    print("  Omega_+^2 / Omega_+^2(sharp) =", main_case["omega_sq_ratio_plus"])
    print("")
    print("So the remaining dynamic 'surface inertia' is not a new phenomenological")
    print("chi-type term. It is a derived Family-1 wall-layer correction controlled by")
    print("the actual wall parameters (alpha0 = V0/mu_c, p, eps_r = d_r/a0).")
    print("")
    print("The static 109/280 monopole closure is unchanged, because only the inertia")
    print("metric is corrected here. The dynamic monopole response remains an exact")
    print("positive two-pole Stieltjes function, and its one-pole Padé reduction stays")
    print("at the ~10^-4 level of accuracy on the natural low-frequency band.")
    print("")
    print("What is still left if you want to push further:")
    print("  1) include the endcap soft-wall layer to correct M_LL beyond the TF core,")
    print("  2) tie the same wall profile to the earlier tangential traction/support law,")
    print("  3) and then fold the derived wall-corrected monopole channel back into the")
    print("     full 2PN throat-response / wake reconstruction.")


if __name__ == "__main__":
    main()
