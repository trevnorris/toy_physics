#!/usr/bin/env python3
"""
same_charge_barrier_audit_stage008_microscopic_xi1_primitive_compiler_sympy_audit.py

Stage 008 — compute the weak-axisymmetric same-charge scalar Xi_1=P1/P0
microscopically from a primitive deformation of the explicit finite-throat one-port
branch, while preserving the conservative grouped response to first order.

What this script does
---------------------
1. Proves the exact arbitrary-base first-order formulas for
      u2^(1), u4^(1), Xi_1,
   and the exact compensation surface that keeps the conservative grouped response
   fixed on a one-pole base branch.
2. Specializes those formulas to the concrete Stage-006 compatibility point of the
   explicit isotropic finite-throat one-port model.
3. Builds the primitive weak-axisymmetric compiler from logarithmic microscopic
   slopes
      (x_K, x_M, x_lambdaB, x_varpi, x_lambdaU, x_lambdaW, x_lambdaR,
       x_OmegaU, x_OmegaW)
   to
      D01, D21, D41, N01, Xi_1.
4. Runs the first mechanism sieve:
      - wall-only family,
      - pure BdG family,
      - mixed-sector-only family.
5. Shows that wall-only and pure BdG families are killed by the conservative
   first-order compensation surface, while the mixed sector retains a nontrivial
   nullspace and therefore a live Xi_1 corridor.
6. Translates one explicit mixed-sector surviving family into the transported
   Stage-007 same-charge ceiling budgets.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import sympy as sp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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
        simplified = expr.applyfunc(lambda z: sp.factor(sp.simplify(z)))
        print(f"{name} =")
        sp.pprint(simplified)
        if any(entry != 0 for entry in simplified):
            raise AssertionError(f"{name} is not zero")
    else:
        simplified = sp.factor(sp.simplify(expr))
        print(f"{name} = {simplified}")
        if simplified != 0:
            raise AssertionError(f"{name} is not zero")


def fmt(x: float) -> str:
    return f"{x:.15f}"


def expect_close(name: str, expr: sp.Expr, tol: float = 1e-30) -> None:
    val = abs(complex(sp.N(expr, 80)))
    print(f"{name} ~= {val}")
    if val > tol:
        raise AssertionError(f"{name} is not within tolerance {tol}")


# ---------------------------------------------------------------------------
# Carried constants from earlier same-charge stages
# ---------------------------------------------------------------------------

KAPPA = 2 * sp.sqrt(2) / sp.pi
P_COMPAT_SAMPLE = sp.Float("0.002069792318062885")
P_BOTH_10 = sp.Float("0.0028313316855593175")
P_ONE_10 = sp.Float("0.0035965105896846573")
P_BOTH_30 = sp.Float("0.00817339430971383")
P_ONE_30 = sp.Float("0.0116633929790174")

SLOPE_NAMES = [
    "xK",
    "xM",
    "xLambdaB",
    "xVarpi",
    "xLambdaU",
    "xLambdaW",
    "xLambdaR",
    "xOmegaU",
    "xOmegaW",
]


def zero_dict() -> dict[str, sp.Expr]:
    return {name: sp.Integer(0) for name in SLOPE_NAMES}


@dataclass(frozen=True)
class PrimitiveParams:
    lamB: sp.Expr
    lamU: sp.Expr
    lamW: sp.Expr
    lamR: sp.Expr
    OmegaU: sp.Expr
    OmegaW: sp.Expr
    varpi: sp.Expr
    M: sp.Expr


@dataclass(frozen=True)
class BaseBundle:
    lamB: sp.Expr
    lamU: sp.Expr
    lamW: sp.Expr
    lamR: sp.Expr
    OmegaU: sp.Expr
    OmegaW: sp.Expr
    varpi: sp.Expr
    M: sp.Expr
    C: sp.Expr
    GU: sp.Expr
    GW: sp.Expr
    R: sp.Expr
    Delta: sp.Expr
    S2: sp.Expr
    H: sp.Expr
    Q: sp.Expr
    P: sp.Expr
    B0: sp.Expr
    B2: sp.Expr
    B4: sp.Expr
    Z0: sp.Expr
    Z2: sp.Expr
    Z4: sp.Expr
    N0: sp.Expr
    Kcompat: sp.Expr
    D0: sp.Expr
    D2: sp.Expr
    D4: sp.Expr
    u2: sp.Expr
    u4: sp.Expr


@dataclass(frozen=True)
class PrimitiveSlopeCompiler:
    D01: dict[str, sp.Expr]
    D21: dict[str, sp.Expr]
    D41: dict[str, sp.Expr]
    N01: dict[str, sp.Expr]
    Xi: dict[str, sp.Expr]
    eq1: dict[str, sp.Expr]
    eq2: dict[str, sp.Expr]


# ---------------------------------------------------------------------------
# Part I. Exact arbitrary-base first-order formulas
# ---------------------------------------------------------------------------

def verify_exact_generic_compensation_surface() -> None:
    banner("PART I — EXACT ARBITRARY-BASE FIRST-ORDER COMPENSATION SURFACE")

    D0, D2, D4 = sp.symbols("D0 D2 D4", positive=True, real=True)
    D01, D21, D41 = sp.symbols("D01 D21 D41", real=True)
    N0, N01 = sp.symbols("N0 N01", positive=True, real=True)

    u2 = sp.simplify(-D2 / D0)
    u4 = sp.simplify((D2**2 - D0 * D4) / D0**2)

    u2_1 = sp.simplify((-D0 * D21 + D2 * D01) / D0**2)
    u4_1 = sp.simplify(
        (-D0 * (D0 * D41 + D01 * D4 - 2 * D2 * D21) + 2 * D01 * (D0 * D4 - D2**2)) / D0**3
    )
    Xi = sp.simplify(N01 / N0 - D01 / D0)

    print("u2 =")
    sp.pprint(u2)
    print("u4 =")
    sp.pprint(u4)
    print("u2^(1) =")
    sp.pprint(u2_1)
    print("u4^(1) =")
    sp.pprint(u4_1)
    print("Xi_1 =")
    sp.pprint(Xi)

    D21_from_u2 = sp.simplify(sp.solve(sp.Eq(u2_1, 0), D21)[0])
    u4_1_reduced = sp.simplify(u4_1.subs(D21, D21_from_u2))
    D41_from_u4 = sp.simplify(sp.solve(sp.Eq(u4_1_reduced, 0), D41)[0])

    print("D21 required by u2^(1)=0 =")
    sp.pprint(D21_from_u2)
    print("u4^(1) after imposing u2^(1)=0 =")
    sp.pprint(u4_1_reduced)
    print("D41 required by u4^(1)=0 as well =")
    sp.pprint(D41_from_u4)

    expect_zero("D21 + u2 D01 on the compensation surface", sp.simplify(D21_from_u2 + u2 * D01))
    expect_zero("D41 - (D4/D0) D01 on the compensation surface", sp.simplify(D41_from_u4 - (D4 / D0) * D01))
    expect_zero("D4/D0 - (u2^2 - u4)", sp.simplify(D4 / D0 - (u2**2 - u4)))
    expect_zero(
        "one-pole specialization D41 + 3 u2^2 D01",
        sp.simplify(D41_from_u4.subs({D4: -3 * u2**2 * D0}) + 3 * u2**2 * D01),
    )


# ---------------------------------------------------------------------------
# Part II. Concrete base branch from Stage 006
# ---------------------------------------------------------------------------

def primitive_base_bundle(p: PrimitiveParams) -> BaseBundle:
    C = sp.simplify(KAPPA * p.lamB)
    GU = sp.simplify(p.lamU)
    GW = sp.simplify(KAPPA * p.lamW)
    R = sp.simplify(KAPPA * p.lamR)

    Delta = sp.simplify(p.OmegaU**2 * p.OmegaW**2 - R**2)
    S2 = sp.simplify(p.OmegaU**2 + p.OmegaW**2)
    H = sp.simplify(GU**2 + GW**2)
    Q = sp.simplify(GU**2 * p.OmegaW**2 + 2 * GU * GW * R + GW**2 * p.OmegaU**2)
    P = sp.simplify(p.OmegaU**2 * GW + R * GU)

    B0 = sp.simplify(C**2 / p.varpi**2)
    B2 = sp.simplify(C**2 / p.varpi**4)
    B4 = sp.simplify(C**2 / p.varpi**6)

    Z0 = sp.simplify(Q / Delta)
    Z2 = sp.simplify((Q * S2 - H * Delta) / Delta**2)
    Z4 = sp.simplify((Q * (S2**2 - Delta) - S2 * H * Delta) / Delta**3)

    N0 = sp.simplify(P**2 / Delta**2)

    Kcompat = sp.simplify(3 * (p.M + B2 + Z2) ** 2 / (B4 + Z4) + B0 + Z0)
    D0 = sp.simplify(Kcompat - B0 - Z0)
    D2 = sp.simplify(-(p.M + B2 + Z2))
    D4 = sp.simplify(-(B4 + Z4))

    u2 = sp.simplify(-D2 / D0)
    u4 = sp.simplify((D2**2 - D0 * D4) / D0**2)

    return BaseBundle(
        lamB=p.lamB,
        lamU=p.lamU,
        lamW=p.lamW,
        lamR=p.lamR,
        OmegaU=p.OmegaU,
        OmegaW=p.OmegaW,
        varpi=p.varpi,
        M=p.M,
        C=C,
        GU=GU,
        GW=GW,
        R=R,
        Delta=Delta,
        S2=S2,
        H=H,
        Q=Q,
        P=P,
        B0=B0,
        B2=B2,
        B4=B4,
        Z0=Z0,
        Z2=Z2,
        Z4=Z4,
        N0=N0,
        Kcompat=Kcompat,
        D0=D0,
        D2=D2,
        D4=D4,
        u2=u2,
        u4=u4,
    )


def print_concrete_base_branch() -> BaseBundle:
    banner("PART II — CONCRETE STAGE-006 COMPATIBILITY BRANCH")

    params = PrimitiveParams(
        lamB=sp.Rational(1, 2),
        lamU=sp.Rational(3, 10),
        lamW=sp.Rational(2, 5),
        lamR=sp.Rational(1, 4),
        OmegaU=sp.Integer(1),
        OmegaW=sp.Rational(7, 5),
        varpi=sp.Integer(2),
        M=sp.Integer(1),
    )
    b = primitive_base_bundle(params)

    print("kappa =")
    sp.pprint(sp.N(KAPPA, 20))
    print("(C, G_U, G_W, R) =")
    sp.pprint(sp.Matrix([sp.N(b.C, 20), sp.N(b.GU, 20), sp.N(b.GW, 20), sp.N(b.R, 20)]))
    print("(Delta, Q, P) =")
    sp.pprint(sp.Matrix([sp.N(b.Delta, 20), sp.N(b.Q, 20), sp.N(b.P, 20)]))
    print("(B0, B2, B4) =")
    sp.pprint(sp.Matrix([sp.N(b.B0, 20), sp.N(b.B2, 20), sp.N(b.B4, 20)]))
    print("(Z0, Z2, Z4) =")
    sp.pprint(sp.Matrix([sp.N(b.Z0, 20), sp.N(b.Z2, 20), sp.N(b.Z4, 20)]))
    print("N0 =")
    sp.pprint(sp.N(b.N0, 20))
    print("K_compat =")
    sp.pprint(sp.N(b.Kcompat, 20))
    print("(D0, D2, D4) =")
    sp.pprint(sp.Matrix([sp.N(b.D0, 20), sp.N(b.D2, 20), sp.N(b.D4, 20)]))
    print("(u2, u4, u4 - 4 u2^2) =")
    sp.pprint(sp.Matrix([sp.N(b.u2, 20), sp.N(b.u4, 20), sp.N(sp.simplify(b.u4 - 4 * b.u2**2), 20)]))

    expect_zero("compatibility-point one-pole defect", sp.simplify(b.u4 - 4 * b.u2**2))
    expect_close("compatibility-point P0_target - carried sample", sp.simplify(b.N0 / b.D0 - P_COMPAT_SAMPLE), tol=1e-15)

    return b


# ---------------------------------------------------------------------------
# Part III. Primitive logarithmic slope compiler
# ---------------------------------------------------------------------------

def primitive_slope_compiler(b: BaseBundle) -> PrimitiveSlopeCompiler:
    dB0 = zero_dict()
    dB2 = zero_dict()
    dB4 = zero_dict()
    dDelta = zero_dict()
    dS2 = zero_dict()
    dH = zero_dict()
    dQ = zero_dict()
    dP = zero_dict()

    # BdG sector
    dB0["xLambdaB"] = sp.simplify(2 * b.B0)
    dB0["xVarpi"] = sp.simplify(-2 * b.B0)
    dB2["xLambdaB"] = sp.simplify(2 * b.B2)
    dB2["xVarpi"] = sp.simplify(-4 * b.B2)
    dB4["xLambdaB"] = sp.simplify(2 * b.B4)
    dB4["xVarpi"] = sp.simplify(-6 * b.B4)

    # Conservative Maxwell/mixed sector
    A0 = sp.simplify(b.OmegaU**2 * b.OmegaW**2)
    dDelta["xOmegaU"] = sp.simplify(2 * A0)
    dDelta["xOmegaW"] = sp.simplify(2 * A0)
    dDelta["xLambdaR"] = sp.simplify(-2 * b.R**2)

    dS2["xOmegaU"] = sp.simplify(2 * b.OmegaU**2)
    dS2["xOmegaW"] = sp.simplify(2 * b.OmegaW**2)

    dH["xLambdaU"] = sp.simplify(2 * b.GU**2)
    dH["xLambdaW"] = sp.simplify(2 * b.GW**2)

    term1 = sp.simplify(b.GU**2 * b.OmegaW**2)
    term2 = sp.simplify(2 * b.GU * b.GW * b.R)
    term3 = sp.simplify(b.GW**2 * b.OmegaU**2)
    dQ["xLambdaU"] = sp.simplify(2 * term1 + term2)
    dQ["xOmegaW"] = sp.simplify(2 * term1)
    dQ["xLambdaW"] = sp.simplify(term2 + 2 * term3)
    dQ["xLambdaR"] = sp.simplify(term2)
    dQ["xOmegaU"] = sp.simplify(2 * term3)

    termP1 = sp.simplify(b.OmegaU**2 * b.GW)
    termP2 = sp.simplify(b.R * b.GU)
    dP["xOmegaU"] = sp.simplify(2 * termP1)
    dP["xLambdaW"] = sp.simplify(termP1)
    dP["xLambdaR"] = sp.simplify(termP2)
    dP["xLambdaU"] = sp.simplify(termP2)

    dZ0 = zero_dict()
    dZ2 = zero_dict()
    dZ4 = zero_dict()
    dN0 = zero_dict()
    N01 = zero_dict()
    D01 = zero_dict()
    D21 = zero_dict()
    D41 = zero_dict()
    Xi = zero_dict()
    eq1 = zero_dict()
    eq2 = zero_dict()

    for name in SLOPE_NAMES:
        dd = dDelta[name]
        ds = dS2[name]
        dh = dH[name]
        dq = dQ[name]
        dp = dP[name]

        dZ0[name] = sp.simplify((dq * b.Delta - b.Q * dd) / b.Delta**2)
        dZ2[name] = sp.simplify(
            (b.Delta * (-b.Delta * dh - b.H * dd + b.Q * ds + b.S2 * dq) + 2 * dd * (b.Delta * b.H - b.Q * b.S2)) / b.Delta**3
        )
        dZ4[name] = sp.simplify(
            -(
                b.Delta**2 * b.H * ds
                + b.Delta**2 * b.S2 * dh
                + b.Delta**2 * dq
                - 2 * b.Delta * b.H * b.S2 * dd
                - 2 * b.Delta * b.Q * b.S2 * ds
                - 2 * b.Delta * b.Q * dd
                - b.Delta * b.S2**2 * dq
                + 3 * b.Q * b.S2**2 * dd
            ) / b.Delta**4
        )
        dN0[name] = sp.simplify(2 * b.P * dp / b.Delta**2 - 2 * b.P**2 * dd / b.Delta**3)

        D01[name] = sp.simplify((b.Kcompat if name == "xK" else 0) - dB0[name] - dZ0[name])
        D21[name] = sp.simplify((-(b.M) if name == "xM" else 0) - dB2[name] - dZ2[name])
        D41[name] = sp.simplify(-dB4[name] - dZ4[name])
        N01[name] = sp.simplify(dN0[name])

        Xi[name] = sp.simplify(N01[name] / b.N0 - D01[name] / b.D0)
        eq1[name] = sp.simplify(D21[name] + b.u2 * D01[name])
        eq2[name] = sp.simplify(D41[name] - (b.D4 / b.D0) * D01[name])

    return PrimitiveSlopeCompiler(D01=D01, D21=D21, D41=D41, N01=N01, Xi=Xi, eq1=eq1, eq2=eq2)


def print_compiler(c: PrimitiveSlopeCompiler) -> None:
    banner("PART III — PRIMITIVE MICROSCOPE: D01, D21, D41, N01, Xi_1")

    def print_family(title: str, data: dict[str, sp.Expr]) -> None:
        print(title)
        for name in SLOPE_NAMES:
            val = sp.N(data[name], 18)
            if val != 0:
                print(f"  {name:9s} -> {val}")
        print()

    print_family("D01 coefficients", c.D01)
    print_family("D21 coefficients", c.D21)
    print_family("D41 coefficients", c.D41)
    print_family("N01 coefficients", c.N01)
    print_family("Xi_1 coefficients", c.Xi)
    print_family("eq1 := D21 + u2 D01 coefficients", c.eq1)
    print_family("eq2 := D41 - (D4/D0) D01 coefficients", c.eq2)


# ---------------------------------------------------------------------------
# Part IV. Mechanism sieve
# ---------------------------------------------------------------------------

def verify_wall_only_no_go_generic() -> None:
    banner("PART IV — WALL-ONLY FAMILY: EXACT GENERIC NO-GO")

    K, M = sp.symbols("K M", positive=True, real=True)
    u2 = sp.symbols("u2", positive=True, real=True)
    r = sp.symbols("r", real=True)  # r = D4 / D0
    xK, xM = sp.symbols("xK xM", real=True)

    eq1 = sp.simplify(-M * xM + u2 * K * xK)
    eq2 = sp.simplify(-r * K * xK)

    print("eq1 (wall-only) =")
    sp.pprint(eq1)
    print("eq2 (wall-only) =")
    sp.pprint(eq2)

    det_wall = sp.simplify(sp.Matrix([[u2 * K, -M], [-r * K, 0]]).det())
    print("Wall-only compensation determinant =")
    sp.pprint(det_wall)

    print("Conclusion: if r = D4/D0 != 0, then eq2 forces xK = 0, and then eq1 forces xM = 0.")


def bdg_only_matrix_and_det(b: BaseBundle, c: PrimitiveSlopeCompiler) -> tuple[sp.Matrix, sp.Expr]:
    cols = ["xLambdaB", "xVarpi"]
    mat = sp.Matrix([
        [sp.simplify(c.eq1[name] / 2) for name in cols],
        [sp.simplify(c.eq2[name] / 2) for name in cols],
    ])

    # Exact compact determinant in terms of B0,B2,B4,u2,D4/D0
    B0, B2, B4, u2, r = sp.symbols("B0 B2 B4 u2 r", real=True)
    generic = sp.Matrix([
        [-(B2 + u2 * B0), 2 * B2 + u2 * B0],
        [-(B4 - r * B0), 3 * B4 - r * B0],
    ])
    det_generic = sp.simplify(generic.det())
    det_sample = sp.simplify(det_generic.subs({B0: b.B0, B2: b.B2, B4: b.B4, u2: b.u2, r: b.D4 / b.D0}))
    return mat, det_sample


def print_bdg_only_no_go(b: BaseBundle, c: PrimitiveSlopeCompiler) -> None:
    banner("PART V — PURE BdG FAMILY: EXACT SAMPLE-POINT NO-GO")

    mat, det_sample = bdg_only_matrix_and_det(b, c)
    print("Pure-BdG compensation matrix on (x_lambdaB, x_varpi) =")
    sp.pprint(sp.N(mat, 20))
    print("Exact pure-BdG determinant formula =")
    B0, B2, B4, u2, r = sp.symbols("B0 B2 B4 u2 r", real=True)
    sp.pprint(sp.simplify(-B0 * B2 * r - 2 * B0 * B4 * u2 - B2 * B4))
    print("Sample-point determinant =")
    sp.pprint(sp.N(det_sample, 20))
    if sp.N(det_sample, 20) == 0:
        raise AssertionError("Pure-BdG determinant vanished unexpectedly on the sample point.")
    print("Conclusion: on the Stage-006 compatibility point, the pure-BdG family has only the trivial compensated solution.")


def mixed_only_nullspace(c: PrimitiveSlopeCompiler) -> tuple[sp.Matrix, list[sp.Matrix], list[sp.Expr]]:
    cols = ["xLambdaU", "xLambdaW", "xLambdaR", "xOmegaU", "xOmegaW"]
    mat = sp.Matrix([
        [sp.N(c.eq1[name], 30) for name in cols],
        [sp.N(c.eq2[name], 30) for name in cols],
    ])
    null = mat.nullspace()
    xis: list[sp.Expr] = []
    for vec in null:
        xi = sp.simplify(sum(sp.N(c.Xi[name], 30) * vec[i] for i, name in enumerate(cols)))
        xis.append(xi)
    return mat, null, xis


def print_mixed_only_survivor(c: PrimitiveSlopeCompiler) -> tuple[sp.Matrix, list[sp.Matrix], list[sp.Expr]]:
    banner("PART VI — MIXED-SECTOR-ONLY FAMILY: EXPLICIT SURVIVING CORRIDOR")

    cols = ["xLambdaU", "xLambdaW", "xLambdaR", "xOmegaU", "xOmegaW"]
    mat, null, xis = mixed_only_nullspace(c)

    print("Mixed-only compensation matrix on (x_lambdaU, x_lambdaW, x_lambdaR, x_OmegaU, x_OmegaW) =")
    sp.pprint(sp.N(mat, 20))
    print(f"rank = {mat.rank()} ; nullity = {len(null)}")

    for i, vec in enumerate(null, start=1):
        print(f"\nNull basis vector v{i} =")
        sp.pprint(sp.N(vec, 20))
        print(f"Xi_1(v{i}) = {sp.N(xis[i-1], 20)}")

    return mat, null, xis


# ---------------------------------------------------------------------------
# Part V. Stage-007 transported ceilings on the surviving mixed family
# ---------------------------------------------------------------------------

def print_transported_ceiling_for_first_mixed_family(xis: list[sp.Expr]) -> None:
    banner("PART VII — TRANSPORTED SAME-CHARGE CEILINGS ON THE FIRST MIXED FAMILY")

    sigma1 = float(abs(sp.N(xis[0], 30)))
    print(f"Choose the first mixed null direction v1 and write the microscopic amplitude as t.")
    print(f"Then Xi_1 = sigma_1 * t with")
    print(f"  sigma_1 = {sigma1:.18f}")
    print()
    print("Transported robust ceiling budgets from Stage 007 become")
    print("  |epsilon * t| <= budget / sigma_1")

    rows = [
        ("10% loss, both wall-like poles", float(P_BOTH_10) / float(P_COMPAT_SAMPLE) - 1.0),
        ("10% loss, nonempty wall-like corridor", float(P_ONE_10) / float(P_COMPAT_SAMPLE) - 1.0),
        ("30% loss, both wall-like poles", float(P_BOTH_30) / float(P_COMPAT_SAMPLE) - 1.0),
        ("30% loss, nonempty wall-like corridor", float(P_ONE_30) / float(P_COMPAT_SAMPLE) - 1.0),
    ]
    for name, xi_budget in rows:
        print(f"[{name}]  |epsilon * t| <= {xi_budget / sigma1:.18f}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    verify_exact_generic_compensation_surface()
    base = print_concrete_base_branch()
    compiler = primitive_slope_compiler(base)
    print_compiler(compiler)
    verify_wall_only_no_go_generic()
    print_bdg_only_no_go(base, compiler)
    _, _, xis = print_mixed_only_survivor(compiler)
    print_transported_ceiling_for_first_mixed_family(xis)

    banner("STAGE 008 LEDGER")
    print("1. The actual same-charge scalar Xi_1 = P1/P0 is now compiled microscopically")
    print("   from primitive weak-axisymmetric slopes of the explicit finite-throat one-port branch.")
    print("2. On any one-pole base branch, preserving the conservative grouped response to first order")
    print("   forces D21 = -u2 D01 and D41 = (D4/D0) D01.")
    print("3. Wall-only deformations die generically once D4/D0 != 0, and pure BdG deformations")
    print("   also die on the Stage-006 compatibility point.")
    print("4. The mixed sector retains a nontrivial compensated nullspace and therefore a live Xi_1 corridor.")
    print("5. The idea survives this stage, but only in a constrained mixed-sector corridor rather than")
    print("   as a pure wall or pure support/BdG mechanism.")


if __name__ == "__main__":
    main()
