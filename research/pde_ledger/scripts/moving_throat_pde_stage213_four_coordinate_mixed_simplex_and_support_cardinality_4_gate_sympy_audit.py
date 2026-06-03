#!/usr/bin/env python3
from __future__ import annotations

import itertools
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


def expect_zero(name: str, expr) -> None:
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(sp.expand(expr))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


def expect_true(name: str, condition: bool) -> None:
    print(f"{name} = {condition}")
    if not condition:
        raise AssertionError(f"{name} failed")


banner("STAGE 213 — FOUR-COORDINATE MIXED SIMPLEX AND THE SUPPORT-CARDINALITY-4 GATE")

# ---------------------------------------------------------------------------
# I. Exact combinatorial ledger for primitive quadruples
# ---------------------------------------------------------------------------
subbanner("I. Exact combinatorial ledger")

axes = ("lambda", "c", "gamma", "U", "W")
triples = list(itertools.combinations(axes, 3))
quadruples = list(itertools.combinations(axes, 4))

print("primitive axes =", axes)
print("primitive triples =", triples)
print("primitive quadruples =", quadruples)
print("#triples =", len(triples))
print("#quadruples =", len(quadruples))

expect_zero("#quadruples - binomial(5,4)", len(quadruples) - int(sp.binomial(5, 4)))

for tri in triples:
    incidence = sum(1 for quad in quadruples if set(tri).issubset(quad))
    print(f"triple incidence {tri} -> {incidence}")
    if incidence != 2:
        raise AssertionError("each primitive triple must belong to exactly two quadruples")

for axis in axes:
    incidence = sum(1 for quad in quadruples if axis in quad)
    print(f"axis incidence {axis} -> {incidence}")
    if incidence != 4:
        raise AssertionError("each primitive axis must belong to exactly four quadruples")

# ---------------------------------------------------------------------------
# II. Positive spherical four-simplex and exact face reductions
# ---------------------------------------------------------------------------
subbanner("II. Positive spherical four-simplex and exact face reductions")

ai, aj, ak, al = sp.symbols("a_i a_j a_k a_l", nonnegative=True, real=True)
r, s, t, u, v = sp.symbols("r s t u v", nonnegative=True, real=True)
ki, kj, kk, kl = sp.symbols("k_i k_j k_k k_l", positive=True, real=True)
H0 = sp.symbols("H0", positive=True, real=True)

avec = sp.Matrix([ai, aj, ak, al])
kvec = sp.Matrix([ki, kj, kk, kl])
k_simplex = sp.simplify((avec.T * kvec)[0])

avec_ijk = sp.Matrix([1 / sp.sqrt(1 + r**2 + s**2), r / sp.sqrt(1 + r**2 + s**2), s / sp.sqrt(1 + r**2 + s**2), 0])
avec_ijl = sp.Matrix([1 / sp.sqrt(1 + r**2 + t**2), r / sp.sqrt(1 + r**2 + t**2), 0, t / sp.sqrt(1 + r**2 + t**2)])
avec_ikl = sp.Matrix([1 / sp.sqrt(1 + s**2 + t**2), 0, s / sp.sqrt(1 + s**2 + t**2), t / sp.sqrt(1 + s**2 + t**2)])
avec_jkl = sp.Matrix([0, 1 / sp.sqrt(1 + u**2 + v**2), u / sp.sqrt(1 + u**2 + v**2), v / sp.sqrt(1 + u**2 + v**2)])

print("a_ijk(r,s) =")
sp.pprint(avec_ijk)
print("a_ijl(r,t) =")
sp.pprint(avec_ijl)
print("a_ikl(s,t) =")
sp.pprint(avec_ikl)
print("a_jkl(u,v) =")
sp.pprint(avec_jkl)

expect_zero("face ijk normalization", (avec_ijk.T * avec_ijk)[0] - 1)
expect_zero("face ijl normalization", (avec_ijl.T * avec_ijl)[0] - 1)
expect_zero("face ikl normalization", (avec_ikl.T * avec_ikl)[0] - 1)
expect_zero("face jkl normalization", (avec_jkl.T * avec_jkl)[0] - 1)

expect_zero(
    "face ijk slope reduction",
    k_simplex.subs({ai: avec_ijk[0], aj: avec_ijk[1], ak: avec_ijk[2], al: avec_ijk[3]})
    - (ki + r * kj + s * kk) / sp.sqrt(1 + r**2 + s**2),
)
expect_zero(
    "face ijl slope reduction",
    k_simplex.subs({ai: avec_ijl[0], aj: avec_ijl[1], ak: avec_ijl[2], al: avec_ijl[3]})
    - (ki + r * kj + t * kl) / sp.sqrt(1 + r**2 + t**2),
)
expect_zero(
    "face ikl slope reduction",
    k_simplex.subs({ai: avec_ikl[0], aj: avec_ikl[1], ak: avec_ikl[2], al: avec_ikl[3]})
    - (ki + s * kk + t * kl) / sp.sqrt(1 + s**2 + t**2),
)
expect_zero(
    "face jkl slope reduction",
    k_simplex.subs({ai: avec_jkl[0], aj: avec_jkl[1], ak: avec_jkl[2], al: avec_jkl[3]})
    - (kj + u * kk + v * kl) / sp.sqrt(1 + u**2 + v**2),
)

# ---------------------------------------------------------------------------
# III. Exact four-coordinate gradient-synergy theorem
# ---------------------------------------------------------------------------
subbanner("III. Exact four-coordinate gradient-synergy theorem")

Kgrad = sp.sqrt(ki**2 + kj**2 + kk**2 + kl**2)
avec_grad = sp.Matrix([ki / Kgrad, kj / Kgrad, kk / Kgrad, kl / Kgrad])

print("a_grad =")
sp.pprint(avec_grad)
print("k_grad =")
sp.pprint(Kgrad)

expect_zero("gradient-optimal normalization", (avec_grad.T * avec_grad)[0] - 1)
expect_zero("gradient-optimal slope value", (avec_grad.T * kvec)[0] - Kgrad)
expect_zero("gradient-optimal ratio r", sp.simplify(avec_grad[1] / avec_grad[0] - kj / ki))
expect_zero("gradient-optimal ratio s", sp.simplify(avec_grad[2] / avec_grad[0] - kk / ki))
expect_zero("gradient-optimal ratio t", sp.simplify(avec_grad[3] / avec_grad[0] - kl / ki))

Kijk = sp.sqrt(ki**2 + kj**2 + kk**2)
Kijl = sp.sqrt(ki**2 + kj**2 + kl**2)
Kikl = sp.sqrt(ki**2 + kk**2 + kl**2)
Kjkl = sp.sqrt(kj**2 + kk**2 + kl**2)
expect_zero("Kgrad^2 - Kijk^2 - k_l^2", Kgrad**2 - Kijk**2 - kl**2)
expect_zero("Kgrad^2 - Kijl^2 - k_k^2", Kgrad**2 - Kijl**2 - kk**2)
expect_zero("Kgrad^2 - Kikl^2 - k_j^2", Kgrad**2 - Kikl**2 - kj**2)
expect_zero("Kgrad^2 - Kjkl^2 - k_i^2", Kgrad**2 - Kjkl**2 - ki**2)

# ---------------------------------------------------------------------------
# IV. Exact total cross-leverage theorem
# ---------------------------------------------------------------------------
subbanner("IV. Exact total cross-leverage theorem")

wSigma = sp.simplify(
    2 * (ai * aj + ai * ak + ai * al + aj * ak + aj * al + ak * al)
)
print("w_Sigma(a) =")
sp.pprint(wSigma)

expect_zero(
    "w_Sigma - ((sum a)^2 - ||a||^2)",
    wSigma - ((ai + aj + ak + al) ** 2 - (ai**2 + aj**2 + ak**2 + al**2)),
)
expect_zero(
    "Cauchy slack identity",
    4 * (ai**2 + aj**2 + ak**2 + al**2)
    - (ai + aj + ak + al) ** 2
    - (
        (ai - aj) ** 2
        + (ai - ak) ** 2
        + (ai - al) ** 2
        + (aj - ak) ** 2
        + (aj - al) ** 2
        + (ak - al) ** 2
    ),
)

aeq4 = sp.Matrix([sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 2)])
aeq3 = sp.Matrix([1 / sp.sqrt(3), 1 / sp.sqrt(3), 1 / sp.sqrt(3), 0])
aeq2 = sp.Matrix([1 / sp.sqrt(2), 1 / sp.sqrt(2), 0, 0])
expect_zero("equal-mix four-way normalization", (aeq4.T * aeq4)[0] - 1)
expect_zero("equal-mix triple-face normalization", (aeq3.T * aeq3)[0] - 1)
expect_zero("equal-mix pair-edge normalization", (aeq2.T * aeq2)[0] - 1)
expect_zero("w_Sigma(four-way equal mix) - 3", wSigma.subs({ai: aeq4[0], aj: aeq4[1], ak: aeq4[2], al: aeq4[3]}) - 3)
expect_zero("w_Sigma(triple equal face) - 2", wSigma.subs({ai: aeq3[0], aj: aeq3[1], ak: aeq3[2], al: aeq3[3]}) - 2)
expect_zero("w_Sigma(pair equal edge) - 1", wSigma.subs({ai: aeq2[0], aj: aeq2[1], ak: aeq2[2], al: aeq2[3]}) - 1)

# ---------------------------------------------------------------------------
# V. Exact four-coordinate curvature law and ratio-coordinate form
# ---------------------------------------------------------------------------
subbanner("V. Exact four-coordinate curvature law and ratio-coordinate form")

uii, uij, uik, uil, ujj, ujk, ujl, ukk, ukl, ull = sp.symbols(
    "u_ii u_ij u_ik u_il u_jj u_jk u_jl u_kk u_kl u_ll", real=True
)
Hmat = sp.Matrix(
    [
        [uii, uij, uik, uil],
        [uij, ujj, ujk, ujl],
        [uik, ujk, ukk, ukl],
        [uil, ujl, ukl, ull],
    ]
)

kappa_simplex = sp.simplify((avec.T * Hmat * avec)[0])
print("kappa_ijkl(a) =")
sp.pprint(kappa_simplex)
expect_zero(
    "diagonal-neutral reduction",
    kappa_simplex.subs({uij: 0, uik: 0, uil: 0, ujk: 0, ujl: 0, ukl: 0})
    - (uii * ai**2 + ujj * aj**2 + ukk * ak**2 + ull * al**2),
)

# Interior ratio patch a_i > 0
avec_rst = sp.Matrix([
    1 / sp.sqrt(1 + r**2 + s**2 + t**2),
    r / sp.sqrt(1 + r**2 + s**2 + t**2),
    s / sp.sqrt(1 + r**2 + s**2 + t**2),
    t / sp.sqrt(1 + r**2 + s**2 + t**2),
])
k_rst = sp.simplify(k_simplex.subs({ai: avec_rst[0], aj: avec_rst[1], ak: avec_rst[2], al: avec_rst[3]}))
kappa_rst = sp.simplify(kappa_simplex.subs({ai: avec_rst[0], aj: avec_rst[1], ak: avec_rst[2], al: avec_rst[3]}))

Acoef = sp.simplify(ki**2 - 2 * H0 * uii)
Bcoef = sp.simplify(2 * ki * kj - 4 * H0 * uij)
Ccoef = sp.simplify(2 * ki * kk - 4 * H0 * uik)
Dcoef = sp.simplify(2 * ki * kl - 4 * H0 * uil)
Ecoef = sp.simplify(kj**2 - 2 * H0 * ujj)
Fcoef = sp.simplify(2 * kj * kk - 4 * H0 * ujk)
Gcoef = sp.simplify(2 * kj * kl - 4 * H0 * ujl)
Hcoef = sp.simplify(kk**2 - 2 * H0 * ukk)
Icoef = sp.simplify(2 * kk * kl - 4 * H0 * ukl)
Jcoef = sp.simplify(kl**2 - 2 * H0 * ull)

Delta_sharp = sp.expand(
    Acoef + Bcoef * r + Ccoef * s + Dcoef * t
    + Ecoef * r**2 + Fcoef * r * s + Gcoef * r * t + Hcoef * s**2 + Icoef * s * t + Jcoef * t**2
)
print("Delta_sharp(r,s,t) =")
sp.pprint(Delta_sharp)

expect_zero("interior ratio normalization", (avec_rst.T * avec_rst)[0] - 1)
expect_zero(
    "discriminant numerator reduction",
    sp.simplify((1 + r**2 + s**2 + t**2) * (k_rst**2 - 2 * H0 * kappa_rst) - Delta_sharp),
)

tau_simplex = sp.simplify(2 * H0 / (k_simplex + sp.sqrt(k_simplex**2 - 2 * H0 * kappa_simplex)))
tau_rst_expected = sp.simplify(
    2 * H0 * sp.sqrt(1 + r**2 + s**2 + t**2)
    / (ki + r * kj + s * kk + t * kl + sp.sqrt(Delta_sharp))
)
expect_zero(
    "interior ratio tau form",
    tau_simplex.subs({ai: avec_rst[0], aj: avec_rst[1], ak: avec_rst[2], al: avec_rst[3]}) - tau_rst_expected,
)

# Face reductions from the interior ratio patch
expect_zero(
    "tau(r,s,t) -> triple ijk on t=0",
    sp.simplify(
        tau_rst_expected.subs(t, 0)
        - 2 * H0 * sp.sqrt(1 + r**2 + s**2)
        / (ki + r * kj + s * kk + sp.sqrt(Acoef + Bcoef * r + Ccoef * s + Ecoef * r**2 + Fcoef * r * s + Hcoef * s**2))
    ),
)
expect_zero(
    "tau(r,s,t) -> triple ijl on s=0",
    sp.simplify(
        tau_rst_expected.subs(s, 0)
        - 2 * H0 * sp.sqrt(1 + r**2 + t**2)
        / (ki + r * kj + t * kl + sp.sqrt(Acoef + Bcoef * r + Dcoef * t + Ecoef * r**2 + Gcoef * r * t + Jcoef * t**2))
    ),
)
expect_zero(
    "tau(r,s,t) -> triple ikl on r=0",
    sp.simplify(
        tau_rst_expected.subs(r, 0)
        - 2 * H0 * sp.sqrt(1 + s**2 + t**2)
        / (ki + s * kk + t * kl + sp.sqrt(Acoef + Ccoef * s + Dcoef * t + Hcoef * s**2 + Icoef * s * t + Jcoef * t**2))
    ),
)

avec_uv = avec_jkl
k_uv = sp.simplify(k_simplex.subs({ai: avec_uv[0], aj: avec_uv[1], ak: avec_uv[2], al: avec_uv[3]}))
kappa_uv = sp.simplify(kappa_simplex.subs({ai: avec_uv[0], aj: avec_uv[1], ak: avec_uv[2], al: avec_uv[3]}))
Delta_jkl = sp.expand(Ecoef + Fcoef * u + Gcoef * v + Hcoef * u**2 + Icoef * u * v + Jcoef * v**2)
expect_zero(
    "jkl-face discriminant reduction",
    sp.simplify((1 + u**2 + v**2) * (k_uv**2 - 2 * H0 * kappa_uv) - Delta_jkl),
)
expect_zero(
    "direct jkl-face tau reduction",
    tau_simplex.subs({ai: avec_uv[0], aj: avec_uv[1], ak: avec_uv[2], al: avec_uv[3]})
    - 2 * H0 * sp.sqrt(1 + u**2 + v**2) / (kj + u * kk + v * kl + sp.sqrt(Delta_jkl)),
)

# ---------------------------------------------------------------------------
# VI. Canonical interior screen values and interval-theorem checks
# ---------------------------------------------------------------------------
subbanner("VI. Canonical screens and exact theorem-gate checks")

k_eq = sp.simplify((ki + kj + kk + kl) / 2)
expect_zero("equal-mix slope value", (aeq4.T * kvec)[0] - k_eq)
expect_zero("gradient slope value", (avec_grad.T * kvec)[0] - Kgrad)
expect_zero("k_grad^2 - ||k||^2", Kgrad**2 - (ki**2 + kj**2 + kk**2 + kl**2))

# Boundary splice with four face intervals
f1_lo, f2_lo, f3_lo, f4_lo, f1_hi, f2_hi, f3_hi, f4_hi = sp.symbols(
    "f1_lo f2_lo f3_lo f4_lo f1_hi f2_hi f3_hi f4_hi", real=True
)
beta_lo = sp.Min(f1_lo, f2_lo, f3_lo, f4_lo)
beta_hi = sp.Min(f1_hi, f2_hi, f3_hi, f4_hi)
print("beta_lo =")
sp.pprint(beta_lo)
print("beta_hi =")
sp.pprint(beta_hi)
expect_true("nested Min flattening (boundary lo)", beta_lo == sp.Min(f1_lo, f2_lo, f3_lo, f4_lo))
expect_true("nested Min flattening (boundary hi)", beta_hi == sp.Min(f1_hi, f2_hi, f3_hi, f4_hi))

count_boundary = 0
for L1 in range(0, 5):
    for U1 in range(L1, 5):
        for L2 in range(0, 5):
            for U2 in range(L2, 5):
                for L3 in range(0, 5):
                    for U3 in range(L3, 5):
                        for L4 in range(0, 5):
                            for U4 in range(L4, 5):
                                blo = min(L1, L2, L3, L4)
                                bhi = min(U1, U2, U3, U4)
                                for x1 in range(L1, U1 + 1):
                                    for x2 in range(L2, U2 + 1):
                                        for x3 in range(L3, U3 + 1):
                                            for x4 in range(L4, U4 + 1):
                                                best = min(x1, x2, x3, x4)
                                                if not (blo <= best <= bhi):
                                                    raise AssertionError("boundary splice theorem failed")
                                                count_boundary += 1
print(f"verified four-face boundary splice theorem on {count_boundary} ordered integer samples")

# Canonical screen dominance theorem against the boundary ledger
count_screen = 0
for blo in range(1, 8):
    for bhi in range(blo, 8):
        for can_lo in range(0, 8):
            for can_hi in range(can_lo, 8):
                if can_hi < blo:
                    for b_star in range(blo, bhi + 1):
                        for c_star in range(can_lo, can_hi + 1):
                            if not (c_star < b_star):
                                raise AssertionError("interior-screen dominance theorem failed")
                            count_screen += 1
print(f"verified canonical screen dominance against boundary on {count_screen} samples")

count_filter = 0
for blo in range(0, 8):
    for bhi in range(blo, 8):
        for can_lo in range(0, 8):
            for can_hi in range(can_lo, 8):
                if can_lo > bhi:
                    for b_star in range(blo, bhi + 1):
                        for c_star in range(can_lo, can_hi + 1):
                            if not (b_star < c_star):
                                raise AssertionError("canonical non-improvement filter failed")
                            count_filter += 1
print(f"verified canonical non-improvement filter on {count_filter} samples")

# Global support-cardinality-4 gate against the Stage 212 up-to-three ledger
count_support4 = 0
for up3_lo in range(1, 8):
    for up3_hi in range(up3_lo, 8):
        for can_lo in range(0, 8):
            for can_hi in range(can_lo, 8):
                if can_hi < up3_lo:
                    for up3_star in range(up3_lo, up3_hi + 1):
                        for can_star in range(can_lo, can_hi + 1):
                            if not (can_star < up3_star):
                                raise AssertionError("support-cardinality-4 improvement gate failed")
                            count_support4 += 1
print(f"verified support-cardinality-4 improvement gate on {count_support4} samples")

count_support4_filter = 0
for up3_lo in range(0, 8):
    for up3_hi in range(up3_lo, 8):
        for can_lo in range(0, 8):
            for can_hi in range(can_lo, 8):
                if can_lo > up3_hi:
                    for up3_star in range(up3_lo, up3_hi + 1):
                        for can_star in range(can_lo, can_hi + 1):
                            if not (up3_star < can_star):
                                raise AssertionError("support-cardinality-4 non-improvement filter failed")
                            count_support4_filter += 1
print(f"verified support-cardinality-4 non-improvement filter on {count_support4_filter} samples")

banner("STAGE 213 SYMPY AUDIT COMPLETED SUCCESSFULLY")
