
import math
from typing import Callable, Tuple

import mpmath as mp
import sympy as sp


def assert_close(actual: float, expected: float, tol: float = 1e-12) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{actual} !~= {expected} (tol={tol})")


def bisect_root(
    func: Callable[[mp.mpf], mp.mpf],
    lo: mp.mpf,
    hi: mp.mpf,
    tol: mp.mpf = mp.mpf("1e-30"),
    max_iter: int = 300,
) -> mp.mpf:
    f_lo = func(lo)
    f_hi = func(hi)
    if f_lo == 0:
        return lo
    if f_hi == 0:
        return hi
    if f_lo > 0 or f_hi < 0:
        raise ValueError(f"Invalid bracket: f(lo)={f_lo}, f(hi)={f_hi}")
    for _ in range(max_iter):
        mid = (lo + hi) / 2
        f_mid = func(mid)
        if abs(f_mid) < tol or abs(hi - lo) < tol:
            return mid
        if f_mid > 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def main() -> None:
    mp.mp.dps = 80
    sp.init_printing()
    print("=== Stage 232 SymPy audit: known 5PN data injection and current branch verdict ===")

    # ------------------------------------------------------------------
    # 1. Refreshed exact Lambda_EM geometry
    # ------------------------------------------------------------------
    x01 = mp.besseljzero(0, 1)
    Lambda_ell = 20 * mp.sqrt(2) * mp.pi / x01
    chi_s = Lambda_ell / 2
    eta = Lambda_ell
    kappa = mp.mpf(9) / 5 * Lambda_ell**2

    print("\nRefreshed Family-1 geometry:")
    print("x01 =", x01)
    print("Lambda_ell =", Lambda_ell)
    print("chi_s =", chi_s)
    print("eta =", eta)
    print("kappa =", kappa)

    assert_close(float(Lambda_ell), 36.94973154240256, tol=5e-14)
    assert_close(float(chi_s), 18.47486577120128, tol=5e-14)
    assert_close(float(kappa), 2457.508789900114, tol=5e-12)

    # ------------------------------------------------------------------
    # 2. Robin support equation and zeta_max ceiling
    # ------------------------------------------------------------------
    def f_robin(y: mp.mpf) -> mp.mpf:
        return y * mp.tan(y) - eta

    y = bisect_root(f_robin, mp.mpf("1e-30"), mp.pi / 2 - mp.mpf("1e-30"))
    A_K = (kappa + mp.pi**2 / 4) / (kappa + y**2)
    zeta_max = A_K * mp.pi**2 / 4

    print("\nRobin support branch:")
    print("y solving y tan y = eta:", y)
    print("A_K =", A_K)
    print("zeta_max =", zeta_max)

    assert_close(float(y), 1.5294278190457656, tol=5e-15)
    assert_close(float(A_K), 1.0000521380385143, tol=5e-15)
    assert_close(float(zeta_max), 2.4675297457259358, tol=5e-15)

    # ------------------------------------------------------------------
    # 3. Exact support-drop kernel and endpoint formulas
    # ------------------------------------------------------------------
    alpha = mp.sqrt(kappa)
    denom = alpha * mp.sinh(alpha) + eta * mp.cosh(alpha)

    def K_kernel(x: mp.mpf) -> mp.mpf:
        return (
            mp.cosh(alpha * x)
            + (eta / alpha) * mp.sinh(alpha * x)
            - mp.cosh(alpha * (1 - x))
        ) / denom

    def Sigma_Pe(x: mp.mpf, Pe: mp.mpf) -> mp.mpf:
        # Stable for very large Pe.
        return Pe * mp.e**(Pe * (x - 1)) / (1 - mp.e**(-Pe))

    def Delta_quad(Pe: mp.mpf) -> mp.mpf:
        integrand = lambda xx: K_kernel(xx) * Sigma_Pe(xx, Pe)
        return mp.quad(integrand, [0, 1])

    Delta_0 = eta * (mp.cosh(alpha) - 1) / (alpha**2 * denom)
    Delta_inf = (mp.cosh(alpha) + (eta / alpha) * mp.sinh(alpha) - 1) / denom

    print("\nSupport-drop endpoints:")
    print("Delta_0 =", Delta_0)
    print("Delta_inf =", Delta_inf)

    assert_close(float(Delta_0), 1.7377393923469950e-4, tol=5e-18)
    assert_close(float(Delta_inf), 2.0172162594593645e-2, tol=5e-16)

    # Stable closed form for Delta(Pe; kappa, eta), obtained by exact integration.
    def Delta_closed(Pe: mp.mpf) -> mp.mpf:
        p = mp.mpf(Pe)
        q = mp.e**(-p)
        Jc = mp.mpf("0.5") * (
            (mp.e**alpha - q) / (p + alpha) + (mp.e**(-alpha) - q) / (p - alpha)
        )
        Js = mp.mpf("0.5") * (
            (mp.e**alpha - q) / (p + alpha) - (mp.e**(-alpha) - q) / (p - alpha)
        )
        return p / (1 - q) * (
            ((1 - mp.cosh(alpha)) * Jc + (eta / alpha + mp.sinh(alpha)) * Js) / denom
        )

    # Cross-check the closed form against direct quadrature on representative points.
    for Pe_test in (mp.mpf("1.0"), mp.mpf("10.0"), mp.mpf("100.0")):
        direct = Delta_quad(Pe_test)
        closed = Delta_closed(Pe_test)
        if abs(direct - closed) > mp.mpf("1e-28"):
            raise AssertionError(f"Delta quadrature mismatch at Pe={Pe_test}: {direct} vs {closed}")

    # Check the endpoint theorem numerically.
    assert abs(Delta_closed(mp.mpf("1e-8")) - Delta_0) < mp.mpf("5e-9")
    assert Delta_0 < Delta_closed(mp.mpf("1000")) < Delta_inf

    print("Verified Delta(Pe) closed form against direct quadrature and endpoint bounds.")

    # ------------------------------------------------------------------
    # 4. Wall-depth extractions and figures of merit
    # ------------------------------------------------------------------
    Theta_w_chi = mp.mpf("4.06863235008162")
    Theta_w_J = mp.mpf("0.927552032539308")

    Xi_chi = 100 * Theta_w_chi * Lambda_ell**2
    Xi_J = 100 * Theta_w_J * Lambda_ell**2

    print("\nWall-depth extracted figures of merit:")
    print("Xi_chi =", Xi_chi)
    print("Xi_J =", Xi_J)

    assert_close(float(Xi_chi), 5.5548332017764099e5, tol=5e-8)
    assert_close(float(Xi_J), 1.2663707072528143e5, tol=5e-8)

    # ------------------------------------------------------------------
    # 5. Exact fixed-point roots Pe = Xi * Delta(Pe; kappa, eta)
    # ------------------------------------------------------------------
    def solve_branch_root(Xi: mp.mpf) -> Tuple[mp.mpf, mp.mpf, mp.mpf]:
        lo = Xi * Delta_0
        hi = Xi * Delta_inf

        def fixed_point_residual(Pe: mp.mpf) -> mp.mpf:
            return Pe - Xi * Delta_closed(Pe)

        root = bisect_root(fixed_point_residual, lo, hi, tol=mp.mpf("1e-40"), max_iter=400)
        return root, lo, hi

    Pe_chi, lo_chi, hi_chi = solve_branch_root(Xi_chi)
    Pe_J, lo_J, hi_J = solve_branch_root(Xi_J)

    print("\nFixed-point roots:")
    print("Pe_*^(chi) =", Pe_chi)
    print("Pe_*^(J) =", Pe_J)
    print("Bracket chi: ", lo_chi, "<= Pe_* <= ", hi_chi)
    print("Bracket J:   ", lo_J, "<= Pe_* <= ", hi_J)

    assert lo_chi <= Pe_chi <= hi_chi
    assert lo_J <= Pe_J <= hi_J
    assert abs(Pe_chi - Xi_chi * Delta_closed(Pe_chi)) < mp.mpf("1e-30")
    assert abs(Pe_J - Xi_J * Delta_closed(Pe_J)) < mp.mpf("1e-30")

    assert_close(float(Pe_chi), 11155.7265863205869, tol=5e-10)
    assert_close(float(Pe_J), 2504.9703142859238, tol=5e-10)

    # ------------------------------------------------------------------
    # 6. Overlap boost and physical support ratios
    # ------------------------------------------------------------------
    def Omega_Pe(Pe: mp.mpf) -> mp.mpf:
        p = mp.mpf(Pe)
        q = mp.e**(-p)
        return mp.pi * p * (2 * p + mp.pi * q) / ((4 * p**2 + mp.pi**2) * (1 - q))

    zeta_phys_chi = A_K * Omega_Pe(Pe_chi) ** 2
    zeta_phys_J = A_K * Omega_Pe(Pe_J) ** 2

    rho_alpha_max_chi = 1 + zeta_phys_chi
    rho_alpha_max_J = 1 + zeta_phys_J

    print("\nPhysical support/source ratios:")
    print("zeta_phys^(chi) =", zeta_phys_chi)
    print("rho_alpha,max^(chi) =", rho_alpha_max_chi)
    print("zeta_phys^(J) =", zeta_phys_J)
    print("rho_alpha,max^(J) =", rho_alpha_max_J)

    assert_close(float(zeta_phys_chi), 2.4675296478814376, tol=5e-15)
    assert_close(float(rho_alpha_max_chi), 3.4675296478814376, tol=5e-15)
    assert_close(float(zeta_phys_J), 2.4675278051675084, tol=5e-15)
    assert_close(float(rho_alpha_max_J), 3.4675278051675084, tol=5e-15)

    # ------------------------------------------------------------------
    # 7. Inject the known 5PN data into the same-charge audit chain
    # ------------------------------------------------------------------
    zeta_req = mp.mpf(1) / 3
    rho_alpha_req = mp.mpf(4) / 3

    margin_zeta_chi = zeta_phys_chi - zeta_req
    margin_zeta_J = zeta_phys_J - zeta_req
    margin_rho_chi = rho_alpha_max_chi - rho_alpha_req
    margin_rho_J = rho_alpha_max_J - rho_alpha_req

    ratio_zeta_chi = zeta_phys_chi / zeta_req
    ratio_zeta_J = zeta_phys_J / zeta_req
    ratio_rho_chi = rho_alpha_max_chi / rho_alpha_req
    ratio_rho_J = rho_alpha_max_J / rho_alpha_req

    gap_to_ceiling_chi = zeta_max - zeta_phys_chi
    gap_to_ceiling_J = zeta_max - zeta_phys_J

    print("\nInjected safety margins:")
    print("margin_zeta^(chi) =", margin_zeta_chi)
    print("margin_rho^(chi) =", margin_rho_chi)
    print("margin_zeta^(J) =", margin_zeta_J)
    print("margin_rho^(J) =", margin_rho_J)

    print("\nUseful ratios:")
    print("zeta_phys^(chi) / zeta_req =", ratio_zeta_chi)
    print("rho_alpha,max^(chi) / (4/3) =", ratio_rho_chi)
    print("zeta_phys^(J) / zeta_req =", ratio_zeta_J)
    print("rho_alpha,max^(J) / (4/3) =", ratio_rho_J)

    print("\nCeiling gaps:")
    print("zeta_max - zeta_phys^(chi) =", gap_to_ceiling_chi)
    print("zeta_max - zeta_phys^(J) =", gap_to_ceiling_J)

    assert margin_zeta_chi > 0 and margin_zeta_J > 0
    assert margin_rho_chi > 0 and margin_rho_J > 0
    assert gap_to_ceiling_chi > 0 and gap_to_ceiling_J > 0

    assert_close(float(margin_zeta_chi), 2.1341963145481043, tol=5e-15)
    assert_close(float(margin_zeta_J), 2.1341944718341751, tol=5e-15)
    assert_close(float(ratio_zeta_chi), 7.402588943644313, tol=5e-15)
    assert_close(float(ratio_zeta_J), 7.402583415502525, tol=5e-15)
    assert_close(float(ratio_rho_chi), 2.600647235911078, tol=5e-15)
    assert_close(float(ratio_rho_J), 2.600645853875631, tol=5e-15)

    # ------------------------------------------------------------------
    # 8. Carried exact same-charge verdict
    # ------------------------------------------------------------------
    print("\nVerdict:")
    print("1. Support/source is numerically located and strongly non-bottlenecked.")
    print("2. The earlier audit ordering is unchanged: dynamic resonance is still not the first kill condition.")
    print("3. The first unresolved gate remains the PDE-selected static orbit-lock / coherent placement packet.")
    print("4. The carried unresolved packet is (d ln R_tr, d ln R_target, d ln epsilon_eta), together with N_Q = 1.")

    print("\nAll Stage 232 checks passed.")


if __name__ == "__main__":
    main()
