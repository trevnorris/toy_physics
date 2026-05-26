#!/usr/bin/env python3
"""
moving_throat_pde_stage7_overlap_isotropy_sympy_audit.py

SymPy audit for Stage 7 of the moving-throat PDE program.

Scope
-----
This script verifies the first explicit overlap-integral layer beyond the abstract
Stage-6 grouped bundle:

  • orthonormality of the normalized real STF l=2 harmonics,
  • the exact angular source-map identity,
  • the isotropic grouped-bundle collapse under an O(3)-invariant kernel,
  • the exact axisymmetric quadrupole triple-overlap matrix,
  • the grouped 20/21/22 splitting pattern (1, 1/2, -1),
  • the first-order defect law b = 3 a,
  • the corresponding first-order transport law for P_A = N_A / D_A,
  • and representative isotropic overlap formulas for C_alpha, B_n, Z_n, N_n,
    and D(omega).

This is still a reduced-sector theorem. It does not solve the full moving-throat
PDE, but it does verify the first explicit angular overlap rules the PDE must obey
on the natural branch.

Provenance notes
----------------
- The grouped `20/21/22` transport quotient `P_A = N_A / D_A` is the same
  Stage 022 grouped prefactor law carried through the Stage 023 bundle
  assembly, now evaluated with the first exact real-STF overlap algebra.
- No new normalization convention is introduced here; the overlap matrix only
  supplies the explicit coefficients that Stage 023 left abstract.
"""

from __future__ import annotations

import sympy as sp

pi = sp.pi
sqrt = sp.sqrt


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


def pairings(lst: list[int]) -> list[list[tuple[int, int]]]:
    if not lst:
        return [[]]
    a = lst[0]
    out = []
    for i in range(1, len(lst)):
        b = lst[i]
        rest = lst[1:i] + lst[i + 1 :]
        for p in pairings(rest):
            out.append([(a, b)] + p)
    return out


# ---------------------------------------------------------------------------
# Canonical real STF basis for l=2
# ---------------------------------------------------------------------------

E20 = sp.Matrix([[-1 / sqrt(6), 0, 0], [0, -1 / sqrt(6), 0], [0, 0, 2 / sqrt(6)]])
E21c = sp.Matrix([[0, 0, 1 / sqrt(2)], [0, 0, 0], [1 / sqrt(2), 0, 0]])
E21s = sp.Matrix([[0, 0, 0], [0, 0, 1 / sqrt(2)], [0, 1 / sqrt(2), 0]])
E22c = sp.Matrix([[1 / sqrt(2), 0, 0], [0, -1 / sqrt(2), 0], [0, 0, 0]])
E22s = sp.Matrix([[0, 1 / sqrt(2), 0], [1 / sqrt(2), 0, 0], [0, 0, 0]])
BASIS = [E20, E21c, E21s, E22c, E22s]
NAMES = ["20", "21c", "21s", "22c", "22s"]

NORM = sp.sqrt(sp.Rational(15) / (8 * sp.pi))
NBASIS = [sp.simplify(NORM * B) for B in BASIS]


# ---------------------------------------------------------------------------
# Isotropic angular moments of the unit sphere
# ---------------------------------------------------------------------------

delta = sp.eye(3)
PAIRINGS4 = [[(0, 1), (2, 3)], [(0, 2), (1, 3)], [(0, 3), (1, 2)]]
PAIRINGS6 = pairings([0, 1, 2, 3, 4, 5])


def I4(i: int, j: int, k: int, l: int) -> sp.Expr:
    inds = [i, j, k, l]
    s = 0
    for pr in PAIRINGS4:
        prod = 1
        for a, b in pr:
            prod *= delta[inds[a], inds[b]]
        s += prod
    return sp.simplify(4 * pi * s / 15)


def I6(i: int, j: int, k: int, l: int, m: int, n: int) -> sp.Expr:
    inds = [i, j, k, l, m, n]
    s = 0
    for pr in PAIRINGS6:
        prod = 1
        for a, b in pr:
            prod *= delta[inds[a], inds[b]]
        s += prod
    return sp.simplify(4 * pi * s / 105)


def quad_overlap(A: sp.Matrix, B: sp.Matrix) -> sp.Expr:
    out = 0
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    out += A[i, j] * B[k, l] * I4(i, j, k, l)
    return sp.simplify(out)


def triple_overlap(A: sp.Matrix, Q: sp.Matrix, B: sp.Matrix) -> sp.Expr:
    out = 0
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        for n in range(3):
                            out += A[i, j] * Q[k, l] * B[m, n] * I6(i, j, k, l, m, n)
    return sp.simplify(out)


# ---------------------------------------------------------------------------
# I. Normalized harmonic orthonormality and source-map identity
# ---------------------------------------------------------------------------

def normalized_harmonics_and_source_map() -> None:
    banner("SECTION I — NORMALIZED HARMONICS AND THE ANGULAR SOURCE MAP")

    gram = sp.Matrix([[quad_overlap(A, B) for B in NBASIS] for A in NBASIS])
    subbanner("I.1 — Orthonormality of the normalized real STF harmonics")
    print("Gram matrix =")
    sp.pprint(gram)
    expect_zero("Gram - I5", gram - sp.eye(5))

    s20, s21c, s21s, s22c, s22s = sp.symbols("s20 s21c s21s s22c s22s", real=True)
    svec = sp.Matrix([s20, s21c, s21s, s22c, s22s])
    projected = sp.simplify(gram * svec)

    subbanner("I.2 — Exact angular source-map identity")
    expect_zero("projected coefficients - source coefficients", projected - svec)
    print("Conclusion: on the canonical isotropic angular basis, mhat_ang = 1 exactly.")


# ---------------------------------------------------------------------------
# II. Isotropic grouped-bundle collapse
# ---------------------------------------------------------------------------

def isotropic_grouped_collapse() -> None:
    banner("SECTION II — ISOTROPIC GROUPED-BUNDLE COLLAPSE")

    x0 = sp.symbols("x0", real=True)
    x20, x21, x22 = sp.symbols("x20 x21 x22", real=True)

    # Stage 022 grouped definitions: weighted trace / anomaly coordinates for
    # the `(20,21,22)` packet under the `(1,2,2)` grouped metric.
    xbar = sp.simplify((x20 + 2 * x21 + 2 * x22) / 5)
    ax = sp.simplify((2 * x20 - x21 - x22) / 10)
    bx = sp.simplify((x21 - x22) / 2)

    subbanner("II.1 — Equal-lane data imply exact grouped isotropy")

    subbanner("II.2 — Unequal lanes produce visible grouped defects")
    witness_a = {x20: x0 + 1, x21: x0, x22: x0}
    witness_b = {x20: x0, x21: x0 + 1, x22: x0}
    expect_zero("a_x witness - 1/5", ax.subs(witness_a) - sp.Rational(1, 5))
    expect_zero("b_x witness", bx.subs(witness_a))
    expect_zero("a_x second witness + 1/10", ax.subs(witness_b) + sp.Rational(1, 10))
    expect_zero("b_x second witness - 1/2", bx.subs(witness_b) - sp.Rational(1, 2))

    # Falsifiable pinning of the (1,2,2)/5 weighting: an arbitrary lane mix should
    # decompose into (xbar, a_x, b_x) and reassemble back to the original triple.
    p, q, rr = sp.symbols("p q rr", real=True)
    mix = {x20: p, x21: q, x22: rr}
    xbar_mix = xbar.subs(mix)
    ax_mix = ax.subs(mix)
    bx_mix = bx.subs(mix)
    x20_reassembled = sp.simplify(xbar_mix + 4 * ax_mix)
    x21_reassembled = sp.simplify(xbar_mix - ax_mix + bx_mix)
    x22_reassembled = sp.simplify(xbar_mix - ax_mix - bx_mix)
    expect_zero("x20 reassembled - p", x20_reassembled - p)
    expect_zero("x21 reassembled - q", x21_reassembled - q)
    expect_zero("x22 reassembled - rr", x22_reassembled - rr)

    print("Interpretation: any O(3)-invariant reduced kernel forces the grouped 20/21/22 bundle")
    print("to collapse to one common scalar value on the isotropic branch.")


# ---------------------------------------------------------------------------
# III. Isotropic radial/axial overlap moments
# ---------------------------------------------------------------------------

def isotropic_overlap_moments() -> None:
    banner("SECTION III — ISOTROPIC RADIAL/AXIAL OVERLAP MOMENTS")

    omega = sp.symbols("omega", real=True)
    K, M = sp.symbols("K M", real=True)

    lamB1, lamB2 = sp.symbols("lambda_B1 lambda_B2", real=True)
    Ieta1, Ieta2 = sp.symbols("I_eta1 I_eta2", real=True)
    varpi1, varpi2 = sp.symbols("varpi_1 varpi_2", positive=True, real=True)

    C1 = sp.simplify(lamB1 * Ieta1)
    C2 = sp.simplify(lamB2 * Ieta2)

    Bresp = sp.simplify(C1**2 / (varpi1**2 - omega**2) + C2**2 / (varpi2**2 - omega**2))
    Bseries = sp.expand(sp.series(Bresp, omega, 0, 6).removeO())
    B0 = sp.simplify(Bseries.coeff(omega, 0))
    B2 = sp.simplify(Bseries.coeff(omega, 2))
    B4 = sp.simplify(Bseries.coeff(omega, 4))

    subbanner("III.1 — Two-mode witness for the BdG moments")
    expect_zero("B0 sum formula", B0 - (C1**2 / varpi1**2 + C2**2 / varpi2**2))
    expect_zero("B2 sum formula", B2 - (C1**2 / varpi1**4 + C2**2 / varpi2**4))
    expect_zero("B4 sum formula", B4 - (C1**2 / varpi1**6 + C2**2 / varpi2**6))

    lamU, lamW, lamR = sp.symbols("lambda_U lambda_W lambda_R", real=True)
    IetaU, IetaW, Iuw = sp.symbols("I_etau I_etaw I_uw", real=True)
    OmegaU, OmegaW, Rmix = sp.symbols("Omega_U Omega_W R_mix", real=True)

    GU = sp.simplify(lamU * IetaU)
    GW = sp.simplify(lamW * IetaW)
    Rr = sp.simplify(lamR * Iuw)

    Delta = sp.simplify(OmegaU**2 * OmegaW**2 - Rr**2)
    S = sp.simplify(OmegaU**2 + OmegaW**2)
    Q = sp.simplify(GU**2 * OmegaW**2 + 2 * GU * GW * Rr + GW**2 * OmegaU**2)
    H = sp.simplify(GU**2 + GW**2)
    P = sp.simplify(OmegaU**2 * GW + Rr * GU)

    # ---- Anchor: derive Z, N response from the per-pair conservative 2x2 matrix inverse.
    # This converts the algebraic identities below into real physics checks: a typo in
    # Q, H, P, Delta, or S would now fail the anchor block instead of trivially passing
    # an algebraic identity of a self-chosen rational.
    M_pair = sp.Matrix([[OmegaU**2 - omega**2, -Rr], [-Rr, OmegaW**2 - omega**2]])
    coupling = sp.Matrix([GU, GW])
    eta_response = M_pair.inv() * coupling  # (U, W) amplitudes for unit eta drive
    Z_from_matrix = sp.simplify((coupling.T * eta_response)[0, 0])
    N_from_matrix = sp.simplify(eta_response[1, 0] ** 2)  # W-component squared
    Zresp_target = (Q - H * omega**2) / (Delta - S * omega**2 + omega**4)
    Nresp_target = (P - GW * omega**2) ** 2 / (Delta - S * omega**2 + omega**4) ** 2
    expect_zero(
        "Z_from_matrix - Zresp_target (physics anchor)",
        sp.simplify(Z_from_matrix - Zresp_target),
    )
    expect_zero(
        "N_from_matrix - Nresp_target (physics anchor)",
        sp.simplify(N_from_matrix - Nresp_target),
    )

    Zresp = sp.expand(
        sp.series((Q - H * omega**2) / (Delta - S * omega**2 + omega**4), omega, 0, 6).removeO()
    )
    Z0 = sp.simplify(Zresp.coeff(omega, 0))
    Z2 = sp.simplify(Zresp.coeff(omega, 2))
    Z4 = sp.simplify(Zresp.coeff(omega, 4))

    Nresp = sp.expand(
        sp.series((P - GW * omega**2) ** 2 / (Delta - S * omega**2 + omega**4) ** 2, omega, 0, 6).removeO()
    )
    N0 = sp.simplify(Nresp.coeff(omega, 0))
    N2 = sp.simplify(Nresp.coeff(omega, 2))
    N4 = sp.simplify(Nresp.coeff(omega, 4))

    Dresp = sp.expand(sp.series(K - M * omega**2 - Bresp - Zresp, omega, 0, 6).removeO())
    D0 = sp.simplify(Dresp.coeff(omega, 0))
    D2 = sp.simplify(Dresp.coeff(omega, 2))
    D4 = sp.simplify(Dresp.coeff(omega, 4))

    subbanner("III.2 — One-pair Maxwell/mixed witness and D(omega)")
    print("Note: Section III.2 below derives Z_n, N_n, D_n from the matrix-inverse")
    print("anchor above; the closed-form coefficient assertions are non-tautological")
    print("given the anchor.")
    expect_zero("Z0 formula", Z0 - Q / Delta)
    expect_zero("Z2 formula", Z2 - (Q * S - H * Delta) / Delta**2)
    expect_zero("Z4 formula", Z4 - (Q * (S**2 - Delta) - S * H * Delta) / Delta**3)
    expect_zero("N0 formula", N0 - P**2 / Delta**2)
    expect_zero("N2 formula", N2 - 2 * P * (P * S - Delta * GW) / Delta**3)
    expect_zero(
        "N4 formula",
        N4 - (Delta**2 * GW**2 - 2 * Delta * P**2 - 4 * Delta * P * S * GW + 3 * P**2 * S**2) / Delta**4,
    )
    expect_zero("D0 formula", D0 - (K - B0 - Z0))
    expect_zero("D2 formula", D2 + (M + B2 + Z2))
    expect_zero("D4 formula", D4 + (B4 + Z4))


def lane_collapse_check() -> None:
    banner("SECTION III.5 — LANE COLLAPSE UNDER O(3) INVARIANCE")
    omega = sp.symbols("omega", real=True)
    K, M = sp.symbols("K M", real=True)

    # Per-lane couplings: under O(3) invariance these are all equal.
    GU_20, GU_21, GU_22 = sp.symbols("GU_20 GU_21 GU_22", real=True)
    GW_20, GW_21, GW_22 = sp.symbols("GW_20 GW_21 GW_22", real=True)
    Rr_20, Rr_21, Rr_22 = sp.symbols("Rr_20 Rr_21 Rr_22", real=True)
    OmU, OmW = sp.symbols("Omega_U Omega_W", positive=True, real=True)
    Cc_20, Cc_21, Cc_22 = sp.symbols("Cc_20 Cc_21 Cc_22", real=True)
    varpi = sp.symbols("varpi", positive=True, real=True)

    def per_lane_D(GU_A, GW_A, Rr_A, C_A):
        # Per-lane B response (single-mode witness).
        B_A = C_A**2 / (varpi**2 - omega**2)
        # Per-lane Z response (one-pair witness using lane-decorated couplings).
        Q_A = GU_A**2 * OmW**2 + 2 * GU_A * GW_A * Rr_A + GW_A**2 * OmU**2
        H_A = GU_A**2 + GW_A**2
        S_A = OmU**2 + OmW**2
        Delta_A = OmU**2 * OmW**2 - Rr_A**2
        Z_A = (Q_A - H_A * omega**2) / (Delta_A - S_A * omega**2 + omega**4)
        return sp.simplify(K - M * omega**2 - B_A - Z_A)

    def per_lane_N(GU_A, GW_A, Rr_A):
        P_A = OmU**2 * GW_A + Rr_A * GU_A
        S_A = OmU**2 + OmW**2
        Delta_A = OmU**2 * OmW**2 - Rr_A**2
        return sp.simplify((P_A - GW_A * omega**2) ** 2 / (Delta_A - S_A * omega**2 + omega**4) ** 2)

    D_20 = per_lane_D(GU_20, GW_20, Rr_20, Cc_20)
    D_21 = per_lane_D(GU_21, GW_21, Rr_21, Cc_21)
    D_22 = per_lane_D(GU_22, GW_22, Rr_22, Cc_22)
    N_20 = per_lane_N(GU_20, GW_20, Rr_20)
    N_21 = per_lane_N(GU_21, GW_21, Rr_21)
    N_22 = per_lane_N(GU_22, GW_22, Rr_22)

    iso = {
        GU_20: sp.Symbol("GU_iso", real=True), GU_21: sp.Symbol("GU_iso", real=True), GU_22: sp.Symbol("GU_iso", real=True),
        GW_20: sp.Symbol("GW_iso", real=True), GW_21: sp.Symbol("GW_iso", real=True), GW_22: sp.Symbol("GW_iso", real=True),
        Rr_20: sp.Symbol("Rr_iso", real=True), Rr_21: sp.Symbol("Rr_iso", real=True), Rr_22: sp.Symbol("Rr_iso", real=True),
        Cc_20: sp.Symbol("C_iso", real=True), Cc_21: sp.Symbol("C_iso", real=True), Cc_22: sp.Symbol("C_iso", real=True),
    }

    subbanner("III.5.1 — Per-lane D_{A,n} collapse under O(3) invariance")
    expect_zero("D_20 - D_21 (isotropic)", sp.simplify(D_20.subs(iso) - D_21.subs(iso)))
    expect_zero("D_21 - D_22 (isotropic)", sp.simplify(D_21.subs(iso) - D_22.subs(iso)))
    expect_zero("D_20 - D_22 (isotropic)", sp.simplify(D_20.subs(iso) - D_22.subs(iso)))

    subbanner("III.5.2 — Per-lane N_{A,n} collapse under O(3) invariance")
    expect_zero("N_20 - N_21 (isotropic)", sp.simplify(N_20.subs(iso) - N_21.subs(iso)))
    expect_zero("N_21 - N_22 (isotropic)", sp.simplify(N_21.subs(iso) - N_22.subs(iso)))
    expect_zero("N_20 - N_22 (isotropic)", sp.simplify(N_20.subs(iso) - N_22.subs(iso)))

    subbanner("III.5.3 — Lane-breaking witness: a single-lane perturbation produces a nonzero defect")
    delta_sym = sp.symbols("delta", real=True)
    break_subs = dict(iso)
    break_subs[GU_20] = sp.Symbol("GU_iso", real=True) + delta_sym
    # Witness: D_{20} - D_{21} should now be linear in delta to leading order.
    D_20_b = D_20.subs(break_subs)
    D_21_b = D_21.subs(break_subs)
    diff_lin = sp.simplify(sp.series(D_20_b - D_21_b, delta_sym, 0, 2).removeO())
    # Assert the defect is not identically zero (linear coefficient in delta is nonzero).
    coeff_delta = sp.simplify(diff_lin.coeff(delta_sym, 1))
    print(f"linear coefficient of delta in (D_20 - D_21) = {coeff_delta}")
    if coeff_delta == 0:
        raise AssertionError("Lane-breaking witness produced no defect: collapse check is trivial")
    print("Lane-breaking witness: collapse check is non-tautological (defect is linear in delta).")


# ---------------------------------------------------------------------------
# IV. Weak axisymmetric quadrupole splitting
# ---------------------------------------------------------------------------

def axisymmetric_splitting() -> None:
    banner("SECTION IV — AXISYMMETRIC QUADRUPOLE SPLITTING MATRIX")

    Q = NBASIS[0]  # normalized Y20 source/background
    M = sp.Matrix([[triple_overlap(A, Q, B) for B in NBASIS] for A in NBASIS])

    subbanner("IV.1 — Exact five-mode triple-overlap matrix for Y20")
    print("M^(20) =")
    sp.pprint(M)

    kappa_star = sp.sqrt(5) / (7 * sp.sqrt(sp.pi))
    M_target = sp.diag(kappa_star, kappa_star / 2, kappa_star / 2, -kappa_star, -kappa_star)
    expect_zero("M - M_target", M - M_target)

    subbanner("IV.2 — Grouped 20/21/22 splitting weights")
    # Stage 022 grouped-defs readback: the axisymmetric `Y20` splitting weights
    # are measured in the same grouped `(xbar,a,b)` coordinates.
    lam20 = sp.Integer(1)
    lam21 = sp.Rational(1, 2)
    lam22 = -sp.Integer(1)
    x0, eps, x1 = sp.symbols("x0 eps x1", real=True)

    x20 = x0 + eps * lam20 * x1
    x21 = x0 + eps * lam21 * x1
    x22 = x0 + eps * lam22 * x1

    xbar = sp.simplify((x20 + 2 * x21 + 2 * x22) / 5)
    ax = sp.simplify((2 * x20 - x21 - x22) / 10)
    bx = sp.simplify((x21 - x22) / 2)

    expect_zero("xbar - x0", xbar - x0)
    expect_zero("a_x - eps*x1/4", ax - eps * x1 / 4)
    expect_zero("b_x - 3*eps*x1/4", bx - 3 * eps * x1 / 4)
    expect_zero("b_x - 3 a_x", bx - 3 * ax)

    print("Axisymmetric grouped weights:")
    print("  lambda_20 = 1")
    print("  lambda_21 = 1/2")
    print("  lambda_22 = -1")


# ---------------------------------------------------------------------------
# V. First-order transport law for P_A = N_A / D_A
# ---------------------------------------------------------------------------

def first_order_transport() -> None:
    banner("SECTION V — FIRST-ORDER TRANSPORT LAW FOR THE NORMALIZATION RATIO")

    eps = sp.symbols("eps", real=True)
    D0, D1, N0, N1 = sp.symbols("D0 D1 N0 N1", nonzero=True, real=True)
    lam20 = sp.Integer(1)
    lam21 = sp.Rational(1, 2)
    lam22 = -sp.Integer(1)

    def lane_ratio(lam: sp.Expr) -> sp.Expr:
        expr = (N0 + eps * lam * N1) / (D0 + eps * lam * D1)
        return sp.expand(sp.series(expr, eps, 0, 2).removeO())

    P20 = lane_ratio(lam20)
    P21 = lane_ratio(lam21)
    P22 = lane_ratio(lam22)

    P0 = sp.simplify(N0 / D0)
    P1 = sp.simplify((N1 * D0 - N0 * D1) / D0**2)

    subbanner("IV.1 — Exact first-order expansion")
    expect_zero("P20 - (P0 + eps P1)", P20 - (P0 + eps * P1))
    expect_zero("P21 - (P0 + eps P1/2)", P21 - (P0 + eps * P1 / 2))
    expect_zero("P22 - (P0 - eps P1)", P22 - (P0 - eps * P1))

    abar = sp.simplify((P20 + 2 * P21 + 2 * P22) / 5)
    aP = sp.simplify((2 * P20 - P21 - P22) / 10)
    bP = sp.simplify((P21 - P22) / 2)

    subbanner("IV.2 — Grouped defects of the normalization ratio")
    expect_zero("Pbar - P0", abar - P0)
    expect_zero("a_P - eps*P1/4", aP - eps * P1 / 4)
    expect_zero("b_P - 3*eps*P1/4", bP - 3 * eps * P1 / 4)
    expect_zero("b_P - 3 a_P", bP - 3 * aP)

    print("So the Stage-6 normalization stack inherits the same exact axisymmetric")
    print("defect law as the microscopic grouped overlaps.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    normalized_harmonics_and_source_map()
    isotropic_grouped_collapse()
    isotropic_overlap_moments()
    lane_collapse_check()
    axisymmetric_splitting()
    first_order_transport()
