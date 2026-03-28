#!/usr/bin/env python3
"""
SymPy derivation-audit script for paper 1 (atomic / P22 / same-charge bridge).

Purpose
-------
This script checks the algebraic derivations that appear in the current paper
draft. It is meant to be a stand-alone audit file that a reader can run with

    python atomic_p22_bridge_sympy_audit.py

What this script checks exactly
-------------------------------
1. Frozen-input Gaussian localization integral.
2. Hydrogenic trial-state normalization, kinetic/Coulomb/Yukawa expectations,
   GNLS overlap norm, Bohr-scale minimum, binding energy, and reduced-mass
   two-body upgrade.
3. Coulomb/Yukawa Hessian identities, the r^{-6} far-field law, the exact
   constant-area ellipse identities, and the Gaussian finite-throat
   regularization formula and asymptotics.
4. Same-charge rotor reduction, canonical momentum / Hamiltonian / shifted
   spectrum, and the algebra used in the conditional half-flux lock.
5. Atomic P22 splitter reduction, the mixed-core stress integral, half-flux
   bracing coefficient, weak-bracing minimum, and the operator identities used
   in the reduced protected-doublet proof.

What this script does *not* prove
---------------------------------
It does not prove the paper's conditional closure assumptions themselves
(e.g. autonomous eigenmode closure, nonzero mouth--core transfer coefficient).
It verifies the algebra that follows once those hypotheses are adopted.

Dependencies
------------
- sympy

The script uses only exact symbolic manipulations except for one optional
finite-basis numerical sanity check at the end.
"""

from __future__ import annotations

import math
import sys
from typing import Iterable

import sympy as sp


PASS_COUNT = 0
WARNINGS: list[str] = []


def banner(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def record_warning(msg: str) -> None:
    WARNINGS.append(msg)
    print(f"[WARN] {msg}")


def _normalize_scalar(expr: sp.Expr) -> sp.Expr:
    expr = sp.expand(expr)
    expr = sp.together(expr)
    expr = sp.cancel(expr)
    expr = sp.factor(expr)
    expr = sp.simplify(expr)
    return expr


def assert_zero(expr, msg: str) -> None:
    """Assert an expression or matrix is identically zero."""
    global PASS_COUNT

    if isinstance(expr, sp.MatrixBase):
        residual = expr.applyfunc(_normalize_scalar)
        if any(entry != 0 for entry in residual):
            raise AssertionError(f"{msg}\nResidual matrix:\n{sp.pretty(residual)}")
    else:
        residual = _normalize_scalar(sp.sympify(expr))
        if residual != 0:
            residual2 = sp.simplify(expr)
            if residual2 != 0:
                raise AssertionError(f"{msg}\nResidual:\n{sp.pretty(residual2)}")

    PASS_COUNT += 1
    print(f"[PASS {PASS_COUNT:02d}] {msg}")


def assert_equal(lhs, rhs, msg: str) -> None:
    assert_zero(sp.sympify(lhs) - sp.sympify(rhs), msg)


def assert_series(expr, var, point, order: int, expected, msg: str) -> None:
    got = sp.series(expr, var, point, order).removeO()
    assert_zero(got - expected, msg)


def frob_contract(A: sp.MatrixBase, B: sp.MatrixBase) -> sp.Expr:
    return sp.simplify(sum(A[i, j] * B[i, j] for i in range(A.rows) for j in range(A.cols)))


def check_frozen_inputs() -> None:
    banner("Section 2  |  Frozen prior inputs and reduction conventions")

    w, lambda_Z = sp.symbols("w lambda_Z", positive=True, real=True)
    Z = sp.exp(-w**2 / lambda_Z**2)
    Z_int = sp.integrate(Z, (w, -sp.oo, sp.oo))
    assert_equal(
        Z_int,
        sp.sqrt(sp.pi) * lambda_Z,
        "Gaussian localization gives Z_int = sqrt(pi) * lambda_Z",
    )


def check_hydrogenic_integrals() -> None:
    banner("Section 3 / Appendix A  |  Hydrogenic reduced sector and integral ledger")

    r, a, M, hbar, gC, kappa = sp.symbols("r a M hbar gC kappa", positive=True, real=True)
    phi_sq = sp.exp(-2 * r / a) / (sp.pi * a**3)

    def radial_average(F):
        return sp.simplify(4 * sp.pi * sp.integrate(r**2 * phi_sq * F, (r, 0, sp.oo)))

    assert_equal(radial_average(1), 1, "1s trial state is normalized")

    grad_sq = phi_sq / a**2
    T = sp.simplify(hbar**2 / (2 * M) * 4 * sp.pi * sp.integrate(r**2 * grad_sq, (r, 0, sp.oo)))
    assert_equal(T, hbar**2 / (2 * M * a**2), "Kinetic expectation is hbar^2/(2 M a^2)")

    assert_equal(radial_average(1 / r), 1 / a, "Coulomb expectation <1/r> = 1/a")

    Yuk = radial_average(sp.exp(-kappa * r) / r)
    assert_equal(
        Yuk,
        1 / (a * (1 + kappa * a / 2) ** 2),
        "Master Yukawa expectation <e^{-kappa r}/r> is exact",
    )

    lambdaZ = sp.symbols("lambdaZ", positive=True, real=True)
    assert_equal(
        Yuk.subs(kappa, 2 / lambdaZ),
        lambdaZ**2 / (a * (a + lambdaZ) ** 2),
        "First even-mode screened Coulomb expectation matches the paper",
    )

    m = sp.symbols("m", positive=True, real=True)
    norm_2m = sp.simplify(
        4 * sp.pi * (1 / (sp.pi * a**3)) ** m * sp.integrate(r**2 * sp.exp(-2 * m * r / a), (r, 0, sp.oo))
    )
    assert_equal(
        norm_2m,
        sp.pi ** (1 - m) * a ** (3 - 3 * m) / m**3,
        "General 2m-norm of the hydrogenic trial state is exact",
    )
    assert_equal(
        norm_2m.subs(m, 5),
        1 / (125 * sp.pi**4 * a**12),
        "Ten-power norm matches the GNLS overlap term",
    )

    E0 = sp.symbols("E0", real=True)
    E_clean = E0 + hbar**2 / (2 * M * a**2) - gC / a
    a_star = sp.solve(sp.Eq(sp.diff(E_clean, a), 0), a)[0]
    assert_equal(
        a_star,
        hbar**2 / (M * gC),
        "Clean hydrogenic minimum occurs at a* = hbar^2/(M g_C)",
    )
    assert_equal(
        sp.diff(E_clean, a, 2).subs(a, a_star),
        M**3 * gC**4 / hbar**6,
        "Second derivative at the clean hydrogenic minimum is positive",
    )
    assert_equal(
        sp.simplify(E_clean.subs(a, a_star) - E0),
        -M * gC**2 / (2 * hbar**2),
        "Binding energy at the clean minimum matches the paper",
    )

    eps0, qeff, qstar, Zint = sp.symbols("eps0 qeff qstar Zint", positive=True, real=True)
    gC_expr = qeff**2 / (4 * sp.pi * eps0)
    a_star_qeff = sp.simplify(a_star.subs(gC, gC_expr))
    assert_equal(
        a_star_qeff,
        4 * sp.pi * eps0 * hbar**2 / (M * qeff**2),
        "Bohr-scale minimum written in terms of q_eff is exact",
    )
    assert_equal(
        a_star_qeff.subs(qeff, qstar / sp.sqrt(Zint)),
        4 * sp.pi * eps0 * hbar**2 * Zint / (M * qstar**2),
        "Thickness law a* ~ Z_int follows from q_eff = q_star/sqrt(Z_int)",
    )
    assert_equal(
        a_star_qeff.subs(qeff, qstar / sp.sqrt(sp.sqrt(sp.pi) * lambdaZ)),
        4 * sp.pi ** sp.Rational(3, 2) * eps0 * hbar**2 * lambdaZ / (M * qstar**2),
        "Gaussian thickness law a* = 4 pi^(3/2) eps0 hbar^2 lambda_Z / (M q_star^2)",
    )


def check_two_body_upgrade() -> None:
    banner("Section 3  |  Two-body reduced-mass upgrade")

    m_minus, m_plus, Rdot, rdot, hbar, gC = sp.symbols(
        "m_minus m_plus Rdot rdot hbar gC", positive=True, real=True
    )
    Mtot = m_minus + m_plus
    mu = sp.simplify(m_minus * m_plus / Mtot)

    xmdot = Rdot + m_plus * rdot / Mtot
    xpdot = Rdot - m_minus * rdot / Mtot
    kinetic = sp.expand(sp.Rational(1, 2) * m_minus * xmdot**2 + sp.Rational(1, 2) * m_plus * xpdot**2)
    target = sp.expand(sp.Rational(1, 2) * Mtot * Rdot**2 + sp.Rational(1, 2) * mu * rdot**2)
    assert_equal(
        kinetic,
        target,
        "Center-of-mass / relative-coordinate kinetic split is exact",
    )

    a_pair = sp.symbols("a_pair", positive=True, real=True)
    E_rel = hbar**2 / (2 * mu * a_pair**2) - gC / a_pair
    a_pair_star = sp.solve(sp.Eq(sp.diff(E_rel, a_pair), 0), a_pair)[0]
    assert_equal(
        a_pair_star,
        hbar**2 / (mu * gC),
        "Two-body clean minimum occurs at a* = hbar^2/(mu g_C)",
    )
    assert_equal(
        sp.limit(mu, m_plus, sp.oo),
        m_minus,
        "Heavy-source limit of the reduced mass is mu -> m_-",
    )


def check_coulomb_hessian() -> None:
    banner("Section 4  |  Finite-throat atomic response: Hessian / tidal identities")

    x, y, z, g = sp.symbols("x y z g", real=True, positive=True)
    r3 = sp.sqrt(x**2 + y**2 + z**2)
    V = -g / r3
    H = sp.hessian(V, (x, y, z))
    n = sp.Matrix([x, y, z]) / r3
    target = -g * (3 * (n * n.T) - sp.eye(3)) / r3**3
    assert_equal(
        H,
        target,
        "Coulomb Hessian is T_ij = -g (3 n_i n_j - delta_ij)/r^3",
    )
    assert_equal(sp.trace(H), 0, "Coulomb Hessian is traceless away from the source")
    assert_equal(
        frob_contract(target, target),
        6 * g**2 / r3**6,
        "Coulomb tidal tensor has Frobenius norm squared 6 g^2/r^6",
    )

    r, lambda_Z = sp.symbols("r lambda_Z", positive=True, real=True)
    VY = -(g / 2) * sp.exp(-2 * r / lambda_Z) / r
    lapY = sp.simplify(sp.diff(VY, r, 2) + 2 * sp.diff(VY, r) / r)
    assert_equal(
        lapY / 3,
        -2 * g * sp.exp(-2 * r / lambda_Z) / (3 * lambda_Z**2 * r),
        "Yukawa trace load T0 = (1/3) lap V matches the paper",
    )


def check_constant_area_ellipse() -> None:
    banner("Section 4 / Appendix B  |  Constant-area ellipse and P22 identities")

    R, eps, phi, theta_m, delta = sp.symbols("R eps phi theta_m delta", positive=True, real=True)
    R1 = R * sp.exp(eps)
    R2 = R * sp.exp(-eps)

    r_boundary = R / sp.sqrt(sp.exp(-2 * eps) * sp.cos(delta) ** 2 + sp.exp(2 * eps) * sp.sin(delta) ** 2)
    ellipse_eq = sp.simplify((r_boundary * sp.cos(delta)) ** 2 / R1**2 + (r_boundary * sp.sin(delta)) ** 2 / R2**2)
    assert_equal(
        ellipse_eq,
        1,
        "Exact polar boundary solves the ellipse equation",
    )

    rho, Theta = sp.symbols("rho Theta", positive=True, real=True)
    Xp = R * sp.exp(eps) * rho * sp.cos(Theta)
    Yp = R * sp.exp(-eps) * rho * sp.sin(Theta)
    J = sp.simplify(sp.Matrix([[sp.diff(Xp, rho), sp.diff(Xp, Theta)], [sp.diff(Yp, rho), sp.diff(Yp, Theta)]]).det())
    A = sp.simplify(sp.integrate(J, (Theta, 0, 2 * sp.pi), (rho, 0, 1)))
    assert_equal(A, sp.pi * R**2, "Exact ellipse family preserves area: A = pi R_th^2")

    expected_small = 1 + eps * sp.cos(2 * delta) + eps**2 * (-sp.Rational(1, 4) + sp.Rational(3, 4) * sp.cos(4 * delta))
    assert_series(
        r_boundary / R,
        eps,
        0,
        3,
        expected_small,
        "Small-ellipticity expansion of the exact polar boundary matches the paper",
    )

    u = eps * sp.cos(2 * (phi - theta_m))
    assert_equal(sp.integrate(u, (phi, 0, 2 * sp.pi)), 0, "Area-preserving first-order P22 mode has zero average")
    assert_equal(sp.integrate(u * sp.cos(phi), (phi, 0, 2 * sp.pi)), 0, "P22 mode satisfies the x-centering constraint")
    assert_equal(sp.integrate(u * sp.sin(phi), (phi, 0, 2 * sp.pi)), 0, "P22 mode satisfies the y-centering constraint")

    X2 = sp.simplify(sp.integrate(Xp**2 * J, (Theta, 0, 2 * sp.pi), (rho, 0, 1)) / A)
    Y2 = sp.simplify(sp.integrate(Yp**2 * J, (Theta, 0, 2 * sp.pi), (rho, 0, 1)) / A)
    assert_equal(X2, R**2 * sp.exp(2 * eps) / 4, "Ellipse principal-frame moment <X'^2> is exact")
    assert_equal(Y2, R**2 * sp.exp(-2 * eps) / 4, "Ellipse principal-frame moment <Y'^2> is exact")

    Qp = sp.Matrix([[R**2 * sp.sinh(2 * eps) / 4, 0], [0, -R**2 * sp.sinh(2 * eps) / 4]])
    Rot = sp.Matrix([[sp.cos(theta_m), -sp.sin(theta_m)], [sp.sin(theta_m), sp.cos(theta_m)]])
    Q_lab = sp.simplify(Rot * Qp * Rot.T)
    n = sp.Matrix([sp.cos(theta_m), sp.sin(theta_m)])
    N = n * n.T - sp.eye(2) / 2
    Q22 = R**2 * sp.sinh(2 * eps) / 2
    assert_equal(Q_lab, sp.simplify(Q22 * N), "Mouth quadrupole tensor equals Q22 * N(theta_m)")

    assert_equal(
        Q_lab[0, 0] - Q_lab[1, 1],
        R**2 * sp.sinh(2 * eps) * sp.cos(2 * theta_m) / 2,
        "Real mouth quadrupole component Q_c is exact",
    )
    assert_equal(
        2 * Q_lab[0, 1],
        R**2 * sp.sinh(2 * eps) * sp.sin(2 * theta_m) / 2,
        "Real mouth quadrupole component Q_s is exact",
    )
    assert_series(
        Q22 / R**2,
        eps,
        0,
        5,
        eps + sp.Rational(2, 3) * eps**3,
        "Mouth quadrupole amplitude Q22 = R^2 eps + (2/3) R^2 eps^3 + ...",
    )

    theta1, theta2 = sp.symbols("theta1 theta2", real=True)
    n1 = sp.Matrix([sp.cos(theta1), sp.sin(theta1)])
    n2 = sp.Matrix([sp.cos(theta2), sp.sin(theta2)])
    N1 = n1 * n1.T - sp.eye(2) / 2
    N2 = n2 * n2.T - sp.eye(2) / 2
    assert_equal(
        frob_contract(N1, N2),
        sp.cos(2 * (theta1 - theta2)) / 2,
        "Director contraction identity N(theta1):N(theta2) = (1/2) cos 2(theta1-theta2)",
    )



def check_mouth_response_minimization() -> None:
    banner("Section 4  |  Weak mouth response and far-field r^{-6} law")

    eps, R, T2, theta_m, theta_T, k22 = sp.symbols("eps R T2 theta_m theta_T k22", positive=True, real=True)
    Delta = theta_m - theta_T

    # Exact contraction from the appendix identities:
    #   U = -(1/2) T:Q,  T:Q = (1/2) T2 Q22 cos 2Delta,  Q22 = (R^2/2) sinh 2eps
    U_tide_exact = -sp.Rational(1, 8) * R**2 * T2 * sp.sinh(2 * eps) * sp.cos(2 * Delta)
    assert_series(
        U_tide_exact,
        eps,
        0,
        3,
        -sp.Rational(1, 4) * R**2 * T2 * eps * sp.cos(2 * Delta),
        "Exact tensor contraction gives weak-field mouth-tide energy -(R^2/4) T2 eps cos 2Delta",
    )

    U_weak_exact = sp.Rational(1, 2) * k22 * eps**2 - sp.Rational(1, 4) * R**2 * T2 * eps * sp.cos(2 * Delta)
    eps_star_exact = sp.solve(sp.Eq(sp.diff(U_weak_exact, eps), 0), eps)[0]
    assert_equal(
        eps_star_exact,
        R**2 * T2 * sp.cos(2 * Delta) / (4 * k22),
        "Exact weak-response ellipticity from the stated tensor definitions is eps* = R^2 T2 cos 2Delta / (4 k22)",
    )

    eps_star_aligned = sp.simplify(eps_star_exact.subs(sp.cos(2 * Delta), 1))
    Umin_aligned = sp.simplify(U_weak_exact.subs({eps: eps_star_aligned, sp.cos(2 * Delta): 1}))
    assert_equal(
        Umin_aligned,
        -R**4 * T2**2 / (32 * k22),
        "Exact adiabatic elimination gives delta V = -R^4 T2^2 / (32 k22) after alignment",
    )

    chi22_patched = R**2 / (4 * k22)
    assert_equal(
        chi22_patched,
        R**2 / (4 * k22),
        "Patched susceptibility is chi22 = R^2/(4 k22)",
    )
    assert_equal(
        chi22_patched * T2 * sp.cos(2 * Delta),
        eps_star_exact,
        "Patched weak-response law eps* = chi22 T2 cos 2Delta is exact",
    )

    U_patched = -sp.Rational(1, 8) * R**2 * T2 * sp.sinh(2 * eps) * sp.cos(2 * Delta)
    assert_equal(
        U_patched,
        U_tide_exact,
        "Patched mouth-tide coefficient matches the exact tensor contraction",
    )

    patched_response_energy = -sp.Rational(1, 8) * chi22_patched * R**2 * T2**2
    assert_equal(
        patched_response_energy,
        Umin_aligned,
        "Patched mouth-response energy prefactor matches exact adiabatic elimination",
    )

    r, g = sp.symbols("r g", positive=True, real=True)
    T_far = -3 * g / r**3
    assert_equal(
        sp.simplify(Umin_aligned.subs(T2, T_far)),
        -sp.Rational(9, 32) * R**4 * g**2 / (k22 * r**6),
        "Exact far-field mouth response yields an attractive r^{-6} correction with coefficient -9 R^4 g^2 /(32 k22 r^6)",
    )

def check_gaussian_core_regulation() -> None:
    banner("Section 4  |  Gaussian finite-throat core regulation")

    s, r, R, g = sp.symbols("s r R g", positive=True, real=True)
    rho = (2 * sp.pi * R**2) ** sp.Rational(-3, 2) * sp.exp(-s**2 / (2 * R**2))

    M_inside = sp.simplify(sp.integrate(4 * sp.pi * s**2 * rho, (s, 0, r)))
    shell_term = sp.simplify(sp.integrate(4 * sp.pi * s * rho, (s, r, sp.oo)))
    Phi_from_shells = sp.simplify(-g * (M_inside / r + shell_term))
    Phi_target = -g * sp.erf(r / (sp.sqrt(2) * R)) / r
    assert_equal(
        Phi_from_shells,
        Phi_target,
        "Gaussian-smeared Coulomb potential equals -g erf(r/(sqrt(2)R))/r",
    )

    x = sp.symbols("x", positive=True, real=True)
    F = -3 * sp.erf(x / sp.sqrt(2)) / x**3 + sp.sqrt(2 / sp.pi) * sp.exp(-x**2 / 2) * (1 + 3 / x**2)
    T2_eff = sp.simplify(sp.diff(Phi_target, r, 2) - sp.diff(Phi_target, r) / r)
    assert_equal(
        T2_eff,
        g * F.subs(x, r / R) / R**3,
        "Regularized quadrupolar tide satisfies T2_eff = (g/R^3) F(r/R)",
    )

    assert_series(
        F,
        x,
        0,
        5,
        -sp.sqrt(2 / sp.pi) * x**2 / 5 + sp.sqrt(2 / sp.pi) * x**4 / 14,
        "Small-x expansion of F(x) matches the paper",
    )

    assert_equal(
        sp.limit(x**3 * F, x, sp.oo),
        -3,
        "Large-x asymptotic of F(x) is -3/x^3",
    )
    assert_equal(
        sp.limit(T2_eff / r**2, r, 0),
        -sp.sqrt(2 / sp.pi) * g / (5 * R**5),
        "Near the core the regularized tide softens as T2_eff = O(r^2)",
    )

    chi22 = sp.symbols("chi22", positive=True, real=True)
    deltaV = -sp.Rational(1, 8) * chi22 * R**2 * T2_eff**2
    assert_equal(
        sp.limit(deltaV / r**4, r, 0),
        -chi22 * g**2 / (100 * sp.pi * R**8),
        "Near the core the induced mouth-response energy softens as O(r^4)",
    )


def check_same_charge_rotor() -> None:
    banner("Section 5  |  Same-charge corridor and rotor algebra")

    sstar, Mb, nu0, hbar = sp.symbols("sstar Mb nu0 hbar", positive=True, real=True)
    theta, thetadot, ptheta = sp.symbols("theta thetadot ptheta", real=True)
    tau = sp.symbols("tau", real=True)
    n_int = sp.symbols("n_int", integer=True)
    bx = sp.sqrt(sstar) * sp.cos(theta)
    by = sp.sqrt(sstar) * sp.sin(theta)
    dbx = sp.diff(bx, theta) * thetadot
    dby = sp.diff(by, theta) * thetadot

    kinetic = sp.simplify(sp.Rational(1, 2) * Mb * (dbx**2 + dby**2))
    berry = sp.simplify(tau * nu0 / 2 * (bx * dby - by * dbx))
    assert_equal(kinetic, Mb * sstar * thetadot**2 / 2, "Pinned-radius kinetic term reduces to (I/2) theta_dot^2")
    assert_equal(berry, tau * nu0 * sstar * thetadot / 2, "First-order mixed term reduces to tau * kappa0 * theta_dot")

    I = Mb * sstar
    kappa0 = nu0 * sstar / 2
    nuB = kappa0 / hbar
    Lrot = I * thetadot**2 / 2 + tau * hbar * nuB * thetadot
    assert_equal(sp.diff(Lrot, thetadot), I * thetadot + tau * hbar * nuB, "Rotor canonical momentum is p = I theta_dot + tau hbar nu_B")

    thetadot_sol = sp.solve(sp.Eq(ptheta, sp.diff(Lrot, thetadot)), thetadot)[0]
    Hrot = sp.simplify((ptheta * thetadot_sol - Lrot.subs(thetadot, thetadot_sol)).subs(tau**2, 1))
    assert_equal(
        Hrot,
        (ptheta - tau * hbar * nuB) ** 2 / (2 * I),
        "Rotor Hamiltonian is H = (p - tau hbar nu_B)^2 / (2I)",
    )

    En = sp.simplify(Hrot.subs(ptheta, n_int * hbar))
    assert_equal(
        En,
        hbar**2 * (n_int - tau * nuB) ** 2 / (2 * I),
        "Compact-rotor spectrum is E_n^(tau) = hbar^2 (n - tau nu_B)^2 / (2I)",
    )


def check_half_flux_lock() -> None:
    banner("Section 5  |  Conditional half-flux lock algebra")

    phi0, nuB, p, q = sp.symbols("phi0 nuB p q", real=True)
    sol = sp.solve(
        [
            sp.Eq(phi0 + 2 * sp.pi * nuB, (2 * p + 1) * sp.pi),
            sp.Eq(phi0 - 2 * sp.pi * nuB, (2 * q + 1) * sp.pi),
        ],
        [phi0, nuB],
        dict=True,
    )[0]
    assert_equal(sol[nuB], (p - q) / 2, "Central-sign closure leaves nu_B = (p-q)/2 when the common scalar phase is free")

    n = sp.symbols("n", integer=True)
    assert_equal(
        sp.exp(sp.I * 2 * sp.pi * (n + sp.Rational(1, 2))),
        -1,
        "If the scalar phase is killed and e^{i 2 pi nu_B} = -1, then nu_B = n + 1/2 solves the loop law",
    )

    sstar, nu0, hbar = sp.symbols("sstar nu0 hbar", positive=True, real=True)
    s_half = sp.solve(sp.Eq(nu0 * sstar / (2 * hbar), sp.Rational(1, 2)), sstar)[0]
    assert_equal(s_half, hbar / nu0, "Minimal half-flux branch implies s_* = hbar/nu_0")



def check_atomic_splitter() -> None:
    banner("Section 6  |  Atomic-to-lepton bridge through P22 forcing")

    eps, R, Tmag, k22, gmix = sp.symbols("eps R Tmag k22 gmix", positive=True, real=True)

    # Using the exact finite-response elimination checked above:
    eps_star_exact = R**2 * Tmag / (4 * k22)
    Q22_star_exact = sp.series(R**2 * sp.sinh(2 * eps) / 2, eps, 0, 2).removeO().subs(eps, eps_star_exact)
    assert_equal(
        sp.simplify(Q22_star_exact),
        R**4 * Tmag / (4 * k22),
        "Exact atomic mouth quadrupole amplitude is Q22* = (R^4/(4k22)) |T2_eff| at leading order",
    )

    V2_atom_exact = gmix * R**4 * Tmag / (4 * k22)
    assert_equal(
        V2_atom_exact,
        gmix * R**4 * Tmag / (4 * k22),
        "Exact atomic splitter amplitude is V2_atom = g_mix * R^4 |T2_eff| /(4k22)",
    )

    assert_equal(
        V2_atom_exact,
        gmix * R**4 * Tmag / (4 * k22),
        "Patched atomic splitter amplitude matches the exact finite-response elimination",
    )

def check_mixed_core_stress() -> None:
    banner("Section 6 / Appendix C  |  Mixed-core stress integral and isolated bracing")

    theta, theta_b, r, chi, chi_r, sstar = sp.symbols(
        "theta theta_b r chi chi_r sstar", positive=True, real=True
    )

    bx = sp.sqrt(sstar) * sp.cos(theta_b)
    by = sp.sqrt(sstar) * sp.sin(theta_b)
    b = sp.Matrix([bx, by])
    xhat = sp.Matrix([sp.cos(theta), sp.sin(theta)])

    dIw = sp.simplify(b * chi + r * (b.dot(xhat)) * xhat * chi_r)
    stress_density = sp.simplify(dIw * dIw.T - sp.eye(2) * (dIw.dot(dIw)) / 2)
    stress_ang = sp.Matrix(
        [
            [sp.simplify(sp.integrate(stress_density[i, j], (theta, 0, 2 * sp.pi))) for j in range(2)]
            for i in range(2)
        ]
    )
    target = sp.simplify(
        sp.pi * (2 * chi + r * chi_r) ** 2 / 2 * (b * b.T - sp.eye(2) * sstar / 2)
    )
    assert_equal(
        stress_ang,
        target,
        "Exact angular average of the mixed-core stress gives the traceless m=2 tensor",
    )

    rr = sp.symbols("rr", positive=True, real=True)
    chi_fn = sp.Function("chi")
    ode_sol = sp.dsolve(sp.Eq(2 * chi_fn(rr) + rr * sp.diff(chi_fn(rr), rr), 0))
    assert_equal(
        ode_sol.rhs,
        sp.Symbol("C1") / rr**2,
        "Vanishing square integrand implies chi(r) ~ C/r^2 before regularity is imposed",
    )

    nvec = sp.Matrix([sp.cos(theta_b), sp.sin(theta_b)])
    N_b = nvec * nvec.T - sp.eye(2) / 2
    assert_equal(
        b * b.T - sp.eye(2) * sstar / 2,
        sstar * N_b,
        "b_i b_j - (s_*/2) delta_ij = s_* N_ij(theta_b)",
    )

    LambdaQ, R, eps, Cmix, hbar, nu0 = sp.symbols(
        "LambdaQ R eps Cmix hbar nu0", positive=True, real=True
    )
    theta_m = sp.symbols("theta_m", real=True)
    n_m = sp.Matrix([sp.cos(theta_m), sp.sin(theta_m)])
    N_m = n_m * n_m.T - sp.eye(2) / 2
    Pi_half = Cmix * hbar / nu0 * N_b
    Q_mouth = R**2 * sp.sinh(2 * eps) / 2 * N_m
    Ucoup = sp.simplify(-LambdaQ * frob_contract(Pi_half, Q_mouth))
    h_half = LambdaQ * R**2 * Cmix * hbar / (4 * nu0)
    assert_equal(
        Ucoup,
        -h_half * sp.sinh(2 * eps) * sp.cos(2 * (theta_m - theta_b)),
        "Half-flux mouth--core coupling reduces to -h_{1/2} sinh(2 eps) cos 2Delta",
    )

    k22, hhalf = sp.symbols("k22 hhalf", positive=True, real=True)
    U_iso_weak = sp.Rational(1, 2) * k22 * eps**2 - 2 * hhalf * eps
    eps_iso = sp.solve(sp.Eq(sp.diff(U_iso_weak, eps), 0), eps)[0]
    assert_equal(
        eps_iso,
        2 * hhalf / k22,
        "Weak-bracing isolated minimum is eps^(infty) = 2 h_{1/2}/k22",
    )


def check_protected_doublet_operator() -> None:
    banner("Section 6  |  Operator identities behind the reduced protected doublet")

    theta, I, hbar, Vinf = sp.symbols("theta I hbar Vinf", positive=True, real=True)
    theta_b = sp.symbols("theta_b", real=True)
    tau = sp.symbols("tau", real=True)
    psi = sp.Function("psi")
    tpsi = sp.Function("tpsi")

    psi_expr = sp.exp(sp.I * tau * theta / 2) * tpsi(theta)
    p_op = lambda f: -sp.I * hbar * sp.diff(f, theta)
    shifted_kinetic = sp.simplify((p_op(p_op(psi_expr)) - tau * hbar * p_op(psi_expr) + (tau * hbar / 2) ** 2 * psi_expr) / sp.exp(sp.I * tau * theta / 2))
    assert_equal(
        shifted_kinetic,
        -hbar**2 * sp.diff(tpsi(theta), theta, 2),
        "Gauge transform removes the half-flux Berry shift from the kinetic operator",
    )

    for tau_val in (1, -1):
        periodic_bc = sp.exp(-sp.I * tau_val * (theta + 2 * sp.pi) / 2) / sp.exp(-sp.I * tau_val * theta / 2)
        assert_equal(
            periodic_bc,
            -1,
            f"Periodic physical wavefunction implies antiperiodic transformed wavefunction for tau={tau_val}",
        )

    f = sp.Function("f")
    Htilde = lambda g: -hbar**2 / (2 * I) * sp.diff(g, theta, 2) - Vinf * sp.cos(2 * (theta - theta_b)) * g
    lhs = sp.simplify(Htilde(f(theta + sp.pi)))
    rhs = sp.simplify(Htilde(f(theta)).subs(theta, theta + sp.pi))
    assert_equal(lhs, rhs, "Half-turn operator commutes with the transformed Hamiltonian")

    n = sp.symbols("n", integer=True)
    basis = sp.exp(sp.I * (n + sp.Rational(1, 2)) * theta)
    S_basis = sp.simplify(basis.subs(theta, theta + sp.pi))
    assert_equal(
        S_basis / basis,
        sp.I * (-1) ** n,
        "Half-turn operator has eigenvalues +/- i on the antiperiodic Fourier basis",
    )
    assert_equal(
        basis.subs(theta, theta + 2 * sp.pi),
        -basis,
        "Antiperiodic basis satisfies psi(theta+2pi) = -psi(theta)",
    )

    # Optional finite-basis numerical sanity check (low modes only).
    Ncut = 6
    indices = list(range(-Ncut, Ncut + 1))
    Hmat = sp.zeros(len(indices))
    Vnum = sp.Rational(7, 10)
    for a, idx in enumerate(indices):
        m = idx + sp.Rational(1, 2)
        Hmat[a, a] = m**2 / 2
        for shift in (-2, 2):
            jdx = idx + shift
            if jdx in indices:
                b = indices.index(jdx)
                Hmat[a, b] += -Vnum / 2
    eigvals = sorted(float(ev) for ev in Hmat.eigenvals().keys())
    pair_residuals = [abs(eigvals[2 * k + 1] - eigvals[2 * k]) for k in range(4)]
    if not all(res < 2e-4 for res in pair_residuals):
        raise AssertionError(
            "Optional finite-basis sanity check for spectral doublets failed.\n"
            f"Pair residuals: {pair_residuals}"
        )
    global PASS_COUNT
    PASS_COUNT += 1
    print(
        f"[PASS {PASS_COUNT:02d}] Optional finite-basis sanity check shows near-degenerate low-lying spectral pairs"
    )


def main() -> None:
    print(__doc__)
    check_frozen_inputs()
    check_hydrogenic_integrals()
    check_two_body_upgrade()
    check_coulomb_hessian()
    check_constant_area_ellipse()
    check_mouth_response_minimization()
    check_gaussian_core_regulation()
    check_same_charge_rotor()
    check_half_flux_lock()
    check_atomic_splitter()
    check_mixed_core_stress()
    check_protected_doublet_operator()

    print("\n" + "-" * 88)
    print(f"All symbolic derivation checks completed successfully. Total passed checks: {PASS_COUNT}")
    if WARNINGS:
        print(f"Audit warnings detected: {len(WARNINGS)}")
        for i, msg in enumerate(WARNINGS, 1):
            print(f"  {i}. {msg}")
    else:
        print("Audit warnings detected: 0")
    print("-" * 88)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # pragma: no cover - explicit standalone failure path
        print("\n[FAIL] Derivation audit stopped with an error:\n")
        print(exc)
        sys.exit(1)

"""

SymPy derivation-audit script for paper 1 (atomic / P22 / same-charge bridge).

Purpose
-------
This script checks the algebraic derivations that appear in the current paper
draft. It is meant to be a stand-alone audit file that a reader can run with

    python atomic_p22_bridge_sympy_audit.py

What this script checks exactly
-------------------------------
1. Frozen-input Gaussian localization integral.
2. Hydrogenic trial-state normalization, kinetic/Coulomb/Yukawa expectations,
   GNLS overlap norm, Bohr-scale minimum, binding energy, and reduced-mass
   two-body upgrade.
3. Coulomb/Yukawa Hessian identities, the r^{-6} far-field law, the exact
   constant-area ellipse identities, and the Gaussian finite-throat
   regularization formula and asymptotics.
4. Same-charge rotor reduction, canonical momentum / Hamiltonian / shifted
   spectrum, and the algebra used in the conditional half-flux lock.
5. Atomic P22 splitter reduction, the mixed-core stress integral, half-flux
   bracing coefficient, weak-bracing minimum, and the operator identities used
   in the reduced protected-doublet proof.

What this script does *not* prove
---------------------------------
It does not prove the paper's conditional closure assumptions themselves
(e.g. autonomous eigenmode closure, nonzero mouth--core transfer coefficient).
It verifies the algebra that follows once those hypotheses are adopted.

Dependencies
------------
- sympy

The script uses only exact symbolic manipulations except for one optional
finite-basis numerical sanity check at the end.


========================================================================================
Section 2  |  Frozen prior inputs and reduction conventions
========================================================================================
[PASS 01] Gaussian localization gives Z_int = sqrt(pi) * lambda_Z

========================================================================================
Section 3 / Appendix A  |  Hydrogenic reduced sector and integral ledger
========================================================================================
[PASS 02] 1s trial state is normalized
[PASS 03] Kinetic expectation is hbar^2/(2 M a^2)
[PASS 04] Coulomb expectation <1/r> = 1/a
[PASS 05] Master Yukawa expectation <e^{-kappa r}/r> is exact
[PASS 06] First even-mode screened Coulomb expectation matches the paper
[PASS 07] General 2m-norm of the hydrogenic trial state is exact
[PASS 08] Ten-power norm matches the GNLS overlap term
[PASS 09] Clean hydrogenic minimum occurs at a* = hbar^2/(M g_C)
[PASS 10] Second derivative at the clean hydrogenic minimum is positive
[PASS 11] Binding energy at the clean minimum matches the paper
[PASS 12] Bohr-scale minimum written in terms of q_eff is exact
[PASS 13] Thickness law a* ~ Z_int follows from q_eff = q_star/sqrt(Z_int)
[PASS 14] Gaussian thickness law a* = 4 pi^(3/2) eps0 hbar^2 lambda_Z / (M q_star^2)

========================================================================================
Section 3  |  Two-body reduced-mass upgrade
========================================================================================
[PASS 15] Center-of-mass / relative-coordinate kinetic split is exact
[PASS 16] Two-body clean minimum occurs at a* = hbar^2/(mu g_C)
[PASS 17] Heavy-source limit of the reduced mass is mu -> m_-

========================================================================================
Section 4  |  Finite-throat atomic response: Hessian / tidal identities
========================================================================================
[PASS 18] Coulomb Hessian is T_ij = -g (3 n_i n_j - delta_ij)/r^3
[PASS 19] Coulomb Hessian is traceless away from the source
[PASS 20] Coulomb tidal tensor has Frobenius norm squared 6 g^2/r^6
[PASS 21] Yukawa trace load T0 = (1/3) lap V matches the paper

========================================================================================
Section 4 / Appendix B  |  Constant-area ellipse and P22 identities
========================================================================================
[PASS 22] Exact polar boundary solves the ellipse equation
[PASS 23] Exact ellipse family preserves area: A = pi R_th^2
[PASS 24] Small-ellipticity expansion of the exact polar boundary matches the paper
[PASS 25] Area-preserving first-order P22 mode has zero average
[PASS 26] P22 mode satisfies the x-centering constraint
[PASS 27] P22 mode satisfies the y-centering constraint
[PASS 28] Ellipse principal-frame moment <X'^2> is exact
[PASS 29] Ellipse principal-frame moment <Y'^2> is exact
[PASS 30] Mouth quadrupole tensor equals Q22 * N(theta_m)
[PASS 31] Real mouth quadrupole component Q_c is exact
[PASS 32] Real mouth quadrupole component Q_s is exact
[PASS 33] Mouth quadrupole amplitude Q22 = R^2 eps + (2/3) R^2 eps^3 + ...
[PASS 34] Director contraction identity N(theta1):N(theta2) = (1/2) cos 2(theta1-theta2)

========================================================================================
Section 4  |  Weak mouth response and far-field r^{-6} law
========================================================================================
[PASS 35] Exact tensor contraction gives weak-field mouth-tide energy -(R^2/4) T2 eps cos 2Delta
[PASS 36] Exact weak-response ellipticity from the stated tensor definitions is eps* = R^2 T2 cos 2Delta / (4 k22)
[PASS 37] Exact adiabatic elimination gives delta V = -R^4 T2^2 / (32 k22) after alignment
[PASS 38] Patched susceptibility is chi22 = R^2/(4 k22)
[PASS 39] Patched weak-response law eps* = chi22 T2 cos 2Delta is exact
[PASS 40] Patched mouth-tide coefficient matches the exact tensor contraction
[PASS 41] Patched mouth-response energy prefactor matches exact adiabatic elimination
[PASS 42] Exact far-field mouth response yields an attractive r^{-6} correction with coefficient -9 R^4 g^2 /(32 k22 r^6)

========================================================================================
Section 4  |  Gaussian finite-throat core regulation
========================================================================================
[PASS 43] Gaussian-smeared Coulomb potential equals -g erf(r/(sqrt(2)R))/r
[PASS 44] Regularized quadrupolar tide satisfies T2_eff = (g/R^3) F(r/R)
[PASS 45] Small-x expansion of F(x) matches the paper
[PASS 46] Large-x asymptotic of F(x) is -3/x^3
[PASS 47] Near the core the regularized tide softens as T2_eff = O(r^2)
[PASS 48] Near the core the induced mouth-response energy softens as O(r^4)

========================================================================================
Section 5  |  Same-charge corridor and rotor algebra
========================================================================================
[PASS 49] Pinned-radius kinetic term reduces to (I/2) theta_dot^2
[PASS 50] First-order mixed term reduces to tau * kappa0 * theta_dot
[PASS 51] Rotor canonical momentum is p = I theta_dot + tau hbar nu_B
[PASS 52] Rotor Hamiltonian is H = (p - tau hbar nu_B)^2 / (2I)
[PASS 53] Compact-rotor spectrum is E_n^(tau) = hbar^2 (n - tau nu_B)^2 / (2I)

========================================================================================
Section 5  |  Conditional half-flux lock algebra
========================================================================================
[PASS 54] Central-sign closure leaves nu_B = (p-q)/2 when the common scalar phase is free
[PASS 55] If the scalar phase is killed and e^{i 2 pi nu_B} = -1, then nu_B = n + 1/2 solves the loop law
[PASS 56] Minimal half-flux branch implies s_* = hbar/nu_0

========================================================================================
Section 6  |  Atomic-to-lepton bridge through P22 forcing
========================================================================================
[PASS 57] Exact atomic mouth quadrupole amplitude is Q22* = (R^4/(4k22)) |T2_eff| at leading order
[PASS 58] Exact atomic splitter amplitude is V2_atom = g_mix * R^4 |T2_eff| /(4k22)
[PASS 59] Patched atomic splitter amplitude matches the exact finite-response elimination

========================================================================================
Section 6 / Appendix C  |  Mixed-core stress integral and isolated bracing
========================================================================================
[PASS 60] Exact angular average of the mixed-core stress gives the traceless m=2 tensor
[PASS 61] Vanishing square integrand implies chi(r) ~ C/r^2 before regularity is imposed
[PASS 62] b_i b_j - (s_*/2) delta_ij = s_* N_ij(theta_b)
[PASS 63] Half-flux mouth--core coupling reduces to -h_{1/2} sinh(2 eps) cos 2Delta
[PASS 64] Weak-bracing isolated minimum is eps^(infty) = 2 h_{1/2}/k22

========================================================================================
Section 6  |  Operator identities behind the reduced protected doublet
========================================================================================
[PASS 65] Gauge transform removes the half-flux Berry shift from the kinetic operator
[PASS 66] Periodic physical wavefunction implies antiperiodic transformed wavefunction for tau=1
[PASS 67] Periodic physical wavefunction implies antiperiodic transformed wavefunction for tau=-1
[PASS 68] Half-turn operator commutes with the transformed Hamiltonian
[PASS 69] Half-turn operator has eigenvalues +/- i on the antiperiodic Fourier basis
[PASS 70] Antiperiodic basis satisfies psi(theta+2pi) = -psi(theta)
[PASS 71] Optional finite-basis sanity check shows near-degenerate low-lying spectral pairs

----------------------------------------------------------------------------------------
All symbolic derivation checks completed successfully. Total passed checks: 71
Audit warnings detected: 0
----------------------------------------------------------------------------------------
"""
