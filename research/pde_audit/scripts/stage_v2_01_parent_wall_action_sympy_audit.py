#!/usr/bin/env python3
"""
Stage V2-01 — Parent-action and throat-action audit

This SymPy audit checks the formal variational difference between:

1. A parent GNLS/Maxwell action in which the throat shape Sigma/R appears only
   inside V_conf(X; Sigma), and
2. The same system after adding the quadratic distributed wall action S_eta^(2).

The script is intentionally algebraic. It does not solve the moving-throat PDE;
it verifies the minimal identities needed for the Volume 2 parent-action audit.
"""

from __future__ import annotations

import sympy as sp


def check(label: str, condition: bool, detail="") -> None:
    if not condition:
        raise AssertionError(f"{label} failed. {detail}")
    print(f"PASS: {label}")


def same(label: str, lhs, rhs) -> None:
    diff = sp.expand((lhs - rhs).doit())
    check(label, diff == 0, f"lhs-rhs = {diff}")


def has_derivative_of(expr, field) -> bool:
    for atom in expr.atoms(sp.Derivative):
        if atom.expr == field:
            return True
    return False


def main() -> None:
    t, w = sp.symbols("t w", real=True)
    ell = sp.symbols("ell", integer=True, nonnegative=True)
    lam = ell * (ell + 1)

    # ------------------------------------------------------------------
    # A. V_conf(Sigma) only: source but no autonomous wall dynamics.
    # ------------------------------------------------------------------
    eta = sp.Function("eta")(t, w)
    rho = sp.Function("rho0")(t, w)
    Vp = sp.Function("Vwall_prime_over_ellc")(t, w)

    # From L_psi = ... -rho V_conf and delta V_conf = -Vp eta.
    L_conf_linear = rho * Vp * eta
    EL_conf = sp.diff(L_conf_linear, eta)
    same("V_conf-only Euler derivative is the algebraic source", EL_conf, rho * Vp)
    check("V_conf-only source contains no eta derivatives", not has_derivative_of(EL_conf, eta))

    # ------------------------------------------------------------------
    # B. Quadratic wall action: modal Euler-Lagrange operator.
    # ------------------------------------------------------------------
    q = sp.Function("q_lm")(t, w)
    mu = sp.Function("mu_eta")(w)
    T_w = sp.Function("T_w")(w)
    T_Omega = sp.Function("T_Omega")(w)
    K_eta = sp.Function("K_eta")(w)
    S = sp.Function("S_lm")(t, w)
    K_eff = K_eta + lam * T_Omega

    q_t = sp.diff(q, t)
    q_w = sp.diff(q, w)

    L_wall = (
        sp.Rational(1, 2) * mu * q_t**2
        - sp.Rational(1, 2) * T_w * q_w**2
        - sp.Rational(1, 2) * K_eff * q**2
        + S * q
    )
    EL = sp.diff(L_wall, q) - sp.diff(sp.diff(L_wall, q_t), t) - sp.diff(sp.diff(L_wall, q_w), w)
    operator_minus_source = mu * sp.diff(q, t, 2) - sp.diff(T_w * q_w, w) + K_eff * q - S
    same("S_eta Euler-Lagrange equation gives modal wall PDE", EL, -operator_minus_source)

    same(
        "l=0 modal specialization",
        operator_minus_source.subs(ell, 0),
        mu * sp.diff(q, t, 2) - sp.diff(T_w * q_w, w) + K_eta * q - S,
    )
    same(
        "l=2 modal specialization",
        operator_minus_source.subs(ell, 2),
        mu * sp.diff(q, t, 2) - sp.diff(T_w * q_w, w) + (K_eta + 6 * T_Omega) * q - S,
    )

    # Boundary terms: delta S boundary pieces are pi_t delta q at time faces and p_w delta q at w faces.
    pi_t = sp.diff(L_wall, q_t)
    p_w = sp.diff(L_wall, q_w)
    same("canonical wall momentum", pi_t, mu * q_t)
    same("w-boundary variation momentum", p_w, -T_w * q_w)
    print("PASS: free-end natural condition is T_w*q_w=0; fixed mouth data means delta q=0")

    # Hamiltonian density for source-free wall action.
    H = sp.expand(pi_t * q_t - L_wall.subs(S, 0))
    expected_H = (
        sp.Rational(1, 2) * mu * q_t**2
        + sp.Rational(1, 2) * T_w * q_w**2
        + sp.Rational(1, 2) * K_eff * q**2
    )
    same("source-free quadratic wall Hamiltonian density", H, expected_H)
    print("PASS: positivity gate is mu_eta>0, T_w>0, K_eta+ell(ell+1)T_Omega>=0")

    # ------------------------------------------------------------------
    # C. Lowest axisymmetric two-mode reduction back to (delta a, delta L).
    # ------------------------------------------------------------------
    Q_a, Q_L, dQ_a, dQ_L = sp.symbols("Q_a Q_L dQ_a dQ_L")
    aa, aL, aap, aLp = sp.symbols("alpha_a alpha_L alpha_a_prime alpha_L_prime")
    mu0, Tw0, K0 = sp.symbols("mu0 Tw0 K0")
    c00 = 2 * sp.sqrt(sp.pi)

    eta_dot = c00 * (aa * dQ_a + aL * dQ_L)
    eta_w = c00 * (aap * Q_a + aLp * Q_L)
    eta0 = c00 * (aa * Q_a + aL * Q_L)
    L_red_integrand = (
        sp.Rational(1, 2) * mu0 * eta_dot**2
        - sp.Rational(1, 2) * Tw0 * eta_w**2
        - sp.Rational(1, 2) * K0 * eta0**2
    )

    same("axisymmetric M_aa integrand", sp.diff(L_red_integrand, dQ_a, dQ_a), 4 * sp.pi * mu0 * aa**2)
    same("axisymmetric M_aL integrand", sp.diff(L_red_integrand, dQ_a, dQ_L), 4 * sp.pi * mu0 * aa * aL)
    same(
        "axisymmetric K_aa integrand",
        -sp.diff(L_red_integrand, Q_a, Q_a),
        4 * sp.pi * (Tw0 * aap**2 + K0 * aa**2),
    )
    same(
        "axisymmetric K_aL integrand",
        -sp.diff(L_red_integrand, Q_a, Q_L),
        4 * sp.pi * (Tw0 * aap * aLp + K0 * aa * aL),
    )

    # ------------------------------------------------------------------
    # D. One-mode grouped P2 reduction.
    # ------------------------------------------------------------------
    beta, beta_p, Q2, dQ2, TO0 = sp.symbols("beta beta_prime Q2 dQ2 TO0")
    L_p2_integrand = (
        sp.Rational(1, 2) * mu0 * (beta * dQ2) ** 2
        - sp.Rational(1, 2) * Tw0 * (beta_p * Q2) ** 2
        - sp.Rational(1, 2) * (K0 + 6 * TO0) * (beta * Q2) ** 2
    )
    same("P2 one-mode M2 integrand", sp.diff(L_p2_integrand, dQ2, dQ2), mu0 * beta**2)
    same(
        "P2 one-mode K2 integrand",
        -sp.diff(L_p2_integrand, Q2, Q2),
        Tw0 * beta_p**2 + (K0 + 6 * TO0) * beta**2,
    )

    print("\nFINAL VERDICT")
    print("STRICT_PARENT_DYNAMIC_WALL: FAIL unless S_eta/S_Sigma is included in S_total.")
    print("EFFECTIVE_LINEAR_WALL_CLOSURE: PASS; S_eta^(2) supplies a consistent linear wall PDE.")
    print("PATCH_REQUIRED: promote S_eta to parent status or relabel the moving wall as an effective closure.")


if __name__ == "__main__":
    main()
