#!/usr/bin/env python3
"""
Step 13 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Starts from the Step-12 direct microscopic monomials
       C_tr,* , C_nt,* , epsilon_eta
   and proves that equality of those invariants between any two positive
   microscopic states is governed by the *same* 3x8 matrix M_* that appeared in
   the infinitesimal compatibility ledger.
2. Solves the exact finite invariant-fibre equations and reconstructs the exact
   five-parameter multiplicative similarity orbit G_*.
3. Builds the tangent generator matrix G for that orbit and verifies
       M_* G = 0.
4. Builds a canonical gauge-fixed quotient section S with
       M_* S = I_3,
   so every microscopic finite log-ratio vector decomposes as
       Delta x = G g + S q,
   where g are five similarity-orbit coordinates and
       q = (Delta ln C_tr, Delta ln C_nt, Delta ln epsilon_eta)
   are the three exact quotient coordinates.
5. Rewrites the coherent weak-axisymmetric observables
       Theta_1, Xi_1, R_1 + Xi_1
   purely in the quotient coordinates.
6. Shows that the carried quartic anomaly target lives only in the
   (q_tr, q_nt)-plane, with q_eta orthogonal to the direct Xi_1 gate.
7. Identifies the simplest canonical quotient representative of the carried
   quartic anomaly target: on the tracking-rigid, dressing-rigid slice the
   entire missing layer is represented by a pure mu_W drift.

Interpretation
--------------
Step 12 located the live microscopic variables. Step 13 removes the remaining
redundancy. The anomaly problem is no longer an 8-dimensional microscopic drift
problem: it is exactly a 3-dimensional quotient problem, and the direct quartic
match itself lives in a 2-dimensional quotient plane.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
    if isinstance(expr, sp.MatrixBase):
        simplified = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(simplified)
        if any(entry != 0 for entry in simplified):
            raise AssertionError(f"{name} is not zero")
    else:
        simplified = sp.simplify(sp.expand(expr))
        print(f"{name} = {simplified}")
        if simplified != 0:
            raise AssertionError(f"{name} is not zero")


def main() -> None:
    banner("STEP 13 — EXACT SIMILARITY ORBIT AND QUOTIENT CLOSURE")

    # ------------------------------------------------------------------
    # I. Direct monomials and exact finite log-ratio equations
    # ------------------------------------------------------------------
    subbanner("XIII.1 — Exact finite invariant-fibre equations")

    chi_star, delta_star = sp.symbols("chi0_star deltaU_star", positive=True, real=True)
    E_star, F_star = sp.symbols("E_star F_star", real=True)

    sigma = sp.symbols("sigma", positive=True, real=True)
    lambda_W, c_etaU, gamma = sp.symbols(
        "lambda_W c_etaU gamma", positive=True, real=True
    )
    K_U, K_eta, K_W = sp.symbols("K_U K_eta K_W", positive=True, real=True)
    mu_W, T_U, L = sp.symbols("mu_W T_U L", positive=True, real=True)

    Delta_l, Delta_c, Delta_g = sp.symbols("Delta_lambda Delta_c Delta_gamma", real=True)
    Delta_U, Delta_eta, Delta_W = sp.symbols("Delta_U Delta_eta Delta_W", real=True)
    Delta_mu, Delta_T = sp.symbols("Delta_mu Delta_T", real=True)

    x = sp.Matrix([
        Delta_l,
        Delta_c,
        Delta_g,
        Delta_U,
        Delta_eta,
        Delta_W,
        Delta_mu,
        Delta_T,
    ])

    # Two positive microscopic states related by finite log-ratios.
    lam_t = lambda_W * sp.exp(Delta_l)
    c_t = c_etaU * sp.exp(Delta_c)
    g_t = gamma * sp.exp(Delta_g)
    KU_t = K_U * sp.exp(Delta_U)
    Keta_t = K_eta * sp.exp(Delta_eta)
    KW_t = K_W * sp.exp(Delta_W)
    mu_t = mu_W * sp.exp(Delta_mu)
    TU_t = T_U * sp.exp(Delta_T)

    Ctr = sp.simplify(
        (gamma * c_etaU / K_U) ** (1 + delta_star)
        * (sp.pi**2 * T_U / (L**2 * K_U)) ** (1 + chi_star)
    )
    Ctr_t = sp.simplify(
        (g_t * c_t / KU_t) ** (1 + delta_star)
        * (sp.pi**2 * TU_t / (L**2 * KU_t)) ** (1 + chi_star)
    )

    Cnt = sp.simplify(
        (lambda_W**2 * mu_W / (K_eta * K_W**2))
        * (gamma**2 * lambda_W**2 * sigma / (K_U * K_W)) ** E_star
        * (sp.pi**2 * T_U / (L**2 * K_U)) ** (-F_star)
    )
    Cnt_t = sp.simplify(
        (lam_t**2 * mu_t / (Keta_t * KW_t**2))
        * (g_t**2 * lam_t**2 * sigma / (KU_t * KW_t)) ** E_star
        * (sp.pi**2 * TU_t / (L**2 * KU_t)) ** (-F_star)
    )

    eps_eta = sp.simplify(c_etaU**2 / (K_U * K_eta))
    eps_eta_t = sp.simplify(c_t**2 / (KU_t * Keta_t))

    q_tr = sp.simplify(sp.log(Ctr_t / Ctr))
    q_nt = sp.simplify(sp.log(Cnt_t / Cnt))
    q_eta = sp.simplify(sp.log(eps_eta_t / eps_eta))

    q = sp.Matrix([q_tr, q_nt, q_eta])

    M = sp.Matrix(
        [
            [
                0,
                1 + delta_star,
                1 + delta_star,
                -(2 + chi_star + delta_star),
                0,
                0,
                0,
                1 + chi_star,
            ],
            [
                2 * (1 + E_star),
                0,
                2 * E_star,
                F_star - E_star,
                -1,
                -(2 + E_star),
                1,
                -F_star,
            ],
            [0, 2, 0, -1, -1, 0, 0, 0],
        ]
    )

    expect_zero("finite quotient drift - M_* Delta x", q - M * x)

    print("q = (Delta ln C_tr, Delta ln C_nt, Delta ln epsilon_eta)^T =")
    sp.pprint(q)
    print("M_* =")
    sp.pprint(M)

    # ------------------------------------------------------------------
    # II. Exact five-parameter similarity orbit and its tangent matrix
    # ------------------------------------------------------------------
    subbanner("XIII.2 — Exact five-parameter similarity orbit")

    alpha = sp.simplify((1 + delta_star) / (1 + chi_star))
    Delta_eta_orb = sp.simplify(2 * Delta_c - Delta_U)
    Delta_T_orb = sp.simplify(Delta_U - alpha * (Delta_g + Delta_c - Delta_U))
    Delta_mu_orb = sp.simplify(
        2 * Delta_c
        - Delta_U
        + 2 * Delta_W
        - 2 * Delta_l
        - E_star * (2 * Delta_g + 2 * Delta_l - Delta_U - Delta_W)
        - F_star * alpha * (Delta_g + Delta_c - Delta_U)
    )

    x_orb = sp.Matrix(
        [
            Delta_l,
            Delta_c,
            Delta_g,
            Delta_U,
            Delta_eta_orb,
            Delta_W,
            Delta_mu_orb,
            Delta_T_orb,
        ]
    )

    expect_zero("orbit invariants M_* x_orb", M * x_orb)

    print("Exact orbit solve:")
    print("  Delta_eta =")
    sp.pprint(Delta_eta_orb)
    print("  Delta_T =")
    sp.pprint(Delta_T_orb)
    print("  Delta_mu =")
    sp.pprint(Delta_mu_orb)

    # Tangent generator matrix for the five free orbit coordinates
    G = sp.Matrix(
        [
            [1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0],
            [0, 0, 0, 1, 0],
            [0, 2, 0, -1, 0],
            [0, 0, 0, 0, 1],
            [
                -2 * (1 + E_star),
                2 - F_star * alpha,
                -2 * E_star - F_star * alpha,
                -1 + E_star + F_star * alpha,
                2 + E_star,
            ],
            [0, -alpha, -alpha, 1 + alpha, 0],
        ]
    )

    expect_zero("M_* G", M * G)
    print("Tangent generator matrix G =")
    sp.pprint(G)

    # ------------------------------------------------------------------
    # III. Exact quotient section and right inverse
    # ------------------------------------------------------------------
    subbanner("XIII.3 — Canonical exact quotient section")

    qtr_sym, qnt_sym, qeta_sym = sp.symbols("q_tr q_nt q_eta", real=True)

    # General exact solution of M_* Delta x = q with the same five free orbit coordinates.
    Delta_eta_gen = sp.simplify(2 * Delta_c - Delta_U - qeta_sym)
    Delta_T_gen = sp.simplify(Delta_T_orb + qtr_sym / (1 + chi_star))
    Delta_mu_gen = sp.simplify(Delta_mu_orb + qnt_sym - qeta_sym + F_star * qtr_sym / (1 + chi_star))

    x_gen = sp.Matrix(
        [
            Delta_l,
            Delta_c,
            Delta_g,
            Delta_U,
            Delta_eta_gen,
            Delta_W,
            Delta_mu_gen,
            Delta_T_gen,
        ]
    )

    q_vec = sp.Matrix([qtr_sym, qnt_sym, qeta_sym])
    expect_zero("general exact solve M_* x_gen - q", M * x_gen - q_vec)

    # Canonical gauge-fixed section: set the five orbit coordinates to zero.
    S = sp.Matrix(
        [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, -1],
            [0, 0, 0],
            [F_star / (1 + chi_star), 1, -1],
            [1 / (1 + chi_star), 0, 0],
        ]
    )

    expect_zero("M_* S - I_3", M * S - sp.eye(3))

    print("Canonical right-inverse / quotient section S =")
    sp.pprint(S)

    # Verify exact finite section directly on the monomials.
    Keta_q = K_eta * sp.exp(-qeta_sym)
    TU_q = T_U * sp.exp(qtr_sym / (1 + chi_star))
    mu_q = mu_W * sp.exp(qnt_sym - qeta_sym + F_star * qtr_sym / (1 + chi_star))

    Ctr_q = sp.simplify(
        (gamma * c_etaU / K_U) ** (1 + delta_star)
        * (sp.pi**2 * TU_q / (L**2 * K_U)) ** (1 + chi_star)
    )
    Cnt_q = sp.simplify(
        (lambda_W**2 * mu_q / (Keta_q * K_W**2))
        * (gamma**2 * lambda_W**2 * sigma / (K_U * K_W)) ** E_star
        * (sp.pi**2 * TU_q / (L**2 * K_U)) ** (-F_star)
    )
    eps_q = sp.simplify(c_etaU**2 / (K_U * Keta_q))

    expect_zero("section Delta ln C_tr - q_tr", sp.log(Ctr_q / Ctr) - qtr_sym)
    expect_zero("section Delta ln C_nt - q_nt", sp.log(Cnt_q / Cnt) - qnt_sym)
    expect_zero("section Delta ln epsilon_eta - q_eta", sp.log(eps_q / eps_eta) - qeta_sym)

    print("Canonical finite section representative:")
    print("  K_eta -> exp(-q_eta) K_eta")
    print("  T_U   -> exp(q_tr/(1+chi_*)) T_U")
    print("  mu_W  -> exp(q_nt - q_eta + F_* q_tr/(1+chi_*)) mu_W")
    print("  all other five microscopic variables fixed")

    # ------------------------------------------------------------------
    # IV. Exact orbit + quotient decomposition
    # ------------------------------------------------------------------
    subbanner("XIII.4 — Exact orbit + quotient decomposition")

    g1, g2, g3, g4, g5 = sp.symbols("g_1 g_2 g_3 g_4 g_5", real=True)
    gvec = sp.Matrix([g1, g2, g3, g4, g5])
    x_decomp = sp.simplify(G * gvec + S * q_vec)

    expect_zero("M_* (G g + S q) - q", M * x_decomp - q_vec)
    expect_zero(
        "general solution - decomposition",
        x_gen.subs(
            {
                Delta_l: g1,
                Delta_c: g2,
                Delta_g: g3,
                Delta_U: g4,
                Delta_W: g5,
            }
        )
        - x_decomp,
    )

    print("Exact decomposition:")
    print("  Delta x = G g + S q")
    print("with orbit coordinates g and quotient coordinates")
    print("  q = (Delta ln C_tr, Delta ln C_nt, Delta ln epsilon_eta)^T")

    # ------------------------------------------------------------------
    # V. Observables and quartic anomaly in quotient coordinates
    # ------------------------------------------------------------------
    subbanner("XIII.5 — Coherent observables and the quartic anomaly in quotient coordinates")

    Ctr_coeff = sp.simplify(
        chi_star * delta_star / ((1 + chi_star) * (1 + delta_star) * (1 + chi_star + delta_star))
    )
    A_tr = sp.simplify(2 * chi_star / ((1 + chi_star) * (1 + delta_star)))
    eps_eta_star = sp.symbols("epsilon_eta_star", positive=True, real=True)

    Theta_q = sp.simplify(-Ctr_coeff * qtr_sym)
    Xi_q = sp.simplify(A_tr * qtr_sym + qnt_sym)
    RplusXi_q = sp.simplify(-eps_eta_star / (1 - eps_eta_star) * qeta_sym)

    print("Theta_1(q) =")
    sp.pprint(Theta_q)
    print("Xi_1(q) =")
    sp.pprint(Xi_q)
    print("R_1 + Xi_1 =")
    sp.pprint(RplusXi_q)

    print("Immediate consequence:")
    print("  The direct quartic anomaly gate lives only in the (q_tr, q_nt)-plane.")
    print("  q_eta is the third exact quotient coordinate, but it is orthogonal to the")
    print("  direct Xi_1 match at this order.")

    # ------------------------------------------------------------------
    # VI. Canonical g-2 matching slices
    # ------------------------------------------------------------------
    subbanner("XIII.6 — Canonical quotient representatives of the quartic anomaly target")

    Lambda1 = sp.symbols("Lambda_1", real=True)
    qeta_free, qtr_free = sp.symbols("qeta_free qtr_free", real=True)
    qnt_match = sp.simplify(Lambda1 - A_tr * qtr_free)

    x_match = sp.simplify(S * sp.Matrix([qtr_free, qnt_match, qeta_free]))

    print("General canonical section of the quartic anomaly target Xi_1 = Lambda_1:")
    sp.pprint(x_match)

    # Simplest direct slice: tracking-rigid and dressing-rigid.
    x_simple = sp.simplify(x_match.subs({qtr_free: 0, qeta_free: 0}))
    expect_zero("simple slice quotient image - (0,Lambda_1,0)", M * x_simple - sp.Matrix([0, Lambda1, 0]))

    print("Tracking-rigid, dressing-rigid canonical representative:")
    print("  q_tr = 0, q_eta = 0, q_nt = Lambda_1")
    print("  Delta x =")
    sp.pprint(x_simple)
    print("So in the chosen quotient gauge the carried quartic target is represented")
    print("entirely by a pure mu_W drift.")

    banner("STEP 13 COMPLETE")


if __name__ == "__main__":
    main()
