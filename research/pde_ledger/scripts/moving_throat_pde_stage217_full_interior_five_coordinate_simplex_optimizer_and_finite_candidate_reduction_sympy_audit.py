#!/usr/bin/env python3
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


banner("STAGE 217 — UNIQUE FIVE-COORDINATE SIMPLEX INTERIOR OPTIMIZER")

# ---------------------------------------------------------------------------
# I. Unique five-coordinate simplex combinatorics
# ---------------------------------------------------------------------------
subbanner("I. Unique support-cardinality-5 simplex and face count")

I5 = ["lambda", "c", "gamma", "U", "W"]
print("Primitive five-coordinate simplex:")
print(I5)

faces = [
    tuple(x for j, x in enumerate(I5) if j != i)
    for i in range(len(I5))
]
print("Codimension-one quadruple faces:")
for face in faces:
    print(face)

if len(faces) != 5:
    raise AssertionError("There should be exactly five codimension-one faces")

# ---------------------------------------------------------------------------
# II. Exact interior chart, discriminant, and stationary numerators
# ---------------------------------------------------------------------------
subbanner("II. Exact interior chart and stationary numerator identities")

H0 = sp.symbols("H0", positive=True, real=True)
r, s, t, u, y = sp.symbols("r s t u y", positive=True, real=True)

kL, kc, kg, kU, kW = sp.symbols("k_lambda k_c k_gamma k_U k_W", positive=True, real=True)

A, B, C, D, E, F, G, H, I, J, K, L, M, N, O = sp.symbols(
    "A B C D E F G H I J K L M N O", real=True
)

S = sp.simplify(1 + r**2 + s**2 + t**2 + u**2)
Klin = sp.simplify(kL + r * kc + s * kg + t * kU + u * kW)
Delta = sp.expand(
    A + B * r + C * s + D * t + E * u
    + F * r**2 + G * r * s + H * r * t + I * r * u
    + J * s**2 + K * s * t + L * s * u
    + M * t**2 + N * t * u + O * u**2
)

Mr = sp.simplify(S * kc - r * Klin)
Ms = sp.simplify(S * kg - s * Klin)
Mt = sp.simplify(S * kU - t * Klin)
Mu = sp.simplify(S * kW - u * Klin)

Lr = sp.simplify(S * sp.diff(Delta, r) - 2 * r * Delta)
Ls = sp.simplify(S * sp.diff(Delta, s) - 2 * s * Delta)
Lt = sp.simplify(S * sp.diff(Delta, t) - 2 * t * Delta)
Lu = sp.simplify(S * sp.diff(Delta, u) - 2 * u * Delta)

Phi = (Klin + sp.sqrt(Delta)) / sp.sqrt(S)

expect_zero(
    "r-stationary numerator identity",
    2 * sp.sqrt(Delta) * S ** sp.Rational(3, 2) * sp.diff(Phi, r) - (2 * Mr * sp.sqrt(Delta) + Lr),
)
expect_zero(
    "s-stationary numerator identity",
    2 * sp.sqrt(Delta) * S ** sp.Rational(3, 2) * sp.diff(Phi, s) - (2 * Ms * sp.sqrt(Delta) + Ls),
)
expect_zero(
    "t-stationary numerator identity",
    2 * sp.sqrt(Delta) * S ** sp.Rational(3, 2) * sp.diff(Phi, t) - (2 * Mt * sp.sqrt(Delta) + Lt),
)
expect_zero(
    "u-stationary numerator identity",
    2 * sp.sqrt(Delta) * S ** sp.Rational(3, 2) * sp.diff(Phi, u) - (2 * Mu * sp.sqrt(Delta) + Lu),
)

# ---------------------------------------------------------------------------
# III. Lifted polynomial system and degree ledger
# ---------------------------------------------------------------------------
subbanner("III. Lifted polynomial stationary system and exact degree ledger")

Fr = sp.expand(2 * Mr * y + Lr)
Fs = sp.expand(2 * Ms * y + Ls)
Ft = sp.expand(2 * Mt * y + Lt)
Fu = sp.expand(2 * Mu * y + Lu)
FDelta = sp.expand(y**2 - Delta)

vars_lift = [r, s, t, u, y]
deg_Fr = sp.Poly(Fr, *vars_lift).total_degree()
deg_Fs = sp.Poly(Fs, *vars_lift).total_degree()
deg_Ft = sp.Poly(Ft, *vars_lift).total_degree()
deg_Fu = sp.Poly(Fu, *vars_lift).total_degree()
deg_FDelta = sp.Poly(FDelta, *vars_lift).total_degree()

print(f"deg(F_r) = {deg_Fr}")
print(f"deg(F_s) = {deg_Fs}")
print(f"deg(F_t) = {deg_Ft}")
print(f"deg(F_u) = {deg_Fu}")
print(f"deg(F_Delta) = {deg_FDelta}")

if (deg_Fr, deg_Fs, deg_Ft, deg_Fu, deg_FDelta) != (3, 3, 3, 3, 2):
    raise AssertionError("Unexpected lifted degree pattern")

bezout_lift = deg_Fr * deg_Fs * deg_Ft * deg_Fu * deg_FDelta
print(f"Lifted Bezout bound = {bezout_lift}")
if bezout_lift != 162:
    raise AssertionError("Unexpected lifted Bezout bound")

# ---------------------------------------------------------------------------
# IV. Direct square-root-free elimination
# ---------------------------------------------------------------------------
subbanner("IV. Direct square-root-free elimination and projected candidate bound")

Crs = sp.expand(Ms * Lr - Mr * Ls)
Crt = sp.expand(Mt * Lr - Mr * Lt)
Cru = sp.expand(Mu * Lr - Mr * Lu)
Sr = sp.expand(Lr**2 - 4 * Mr**2 * Delta)

vars_proj = [r, s, t, u]
deg_Crs = sp.Poly(Crs, *vars_proj).total_degree()
deg_Crt = sp.Poly(Crt, *vars_proj).total_degree()
deg_Cru = sp.Poly(Cru, *vars_proj).total_degree()
deg_Sr = sp.Poly(Sr, *vars_proj).total_degree()

print(f"deg(C_rs) = {deg_Crs}")
print(f"deg(C_rt) = {deg_Crt}")
print(f"deg(C_ru) = {deg_Cru}")
print(f"deg(S_r) = {deg_Sr}")

if (deg_Crs, deg_Crt, deg_Cru, deg_Sr) != (5, 5, 5, 6):
    raise AssertionError("Unexpected projected degree pattern")

bezout_proj = deg_Crs * deg_Crt * deg_Cru * deg_Sr
print(f"Projected one-chart Bezout bound = {bezout_proj}")
if bezout_proj != 750:
    raise AssertionError("Unexpected projected Bezout bound")

# ---------------------------------------------------------------------------
# V. Special reduction 1: diagonal-isotropic curvature
# ---------------------------------------------------------------------------
subbanner("V. Special reduction 1 — diagonal-isotropic curvature")

nu = sp.symbols("nu", real=True)
diag_subs = {
    A: kL**2 - 2 * H0 * nu,
    B: 2 * kL * kc,
    C: 2 * kL * kg,
    D: 2 * kL * kU,
    E: 2 * kL * kW,
    F: kc**2 - 2 * H0 * nu,
    G: 2 * kc * kg,
    H: 2 * kc * kU,
    I: 2 * kc * kW,
    J: kg**2 - 2 * H0 * nu,
    K: 2 * kg * kU,
    L: 2 * kg * kW,
    M: kU**2 - 2 * H0 * nu,
    N: 2 * kU * kW,
    O: kW**2 - 2 * H0 * nu,
}
grad_ratios = {
    r: kc / kL,
    s: kg / kL,
    t: kU / kL,
    u: kW / kL,
}

Mr_diag = sp.simplify(Mr.subs(diag_subs))
Ms_diag = sp.simplify(Ms.subs(diag_subs))
Mt_diag = sp.simplify(Mt.subs(diag_subs))
Mu_diag = sp.simplify(Mu.subs(diag_subs))
Delta_diag = sp.simplify(Delta.subs(diag_subs))
Lr_diag = sp.simplify(Lr.subs(diag_subs))
Ls_diag = sp.simplify(Ls.subs(diag_subs))
Lt_diag = sp.simplify(Lt.subs(diag_subs))
Lu_diag = sp.simplify(Lu.subs(diag_subs))
Klin_diag = Klin

expect_zero("L_r(diag) - 2 K_lin M_r", sp.simplify(Lr_diag - 2 * Klin_diag * Mr_diag))
expect_zero("L_s(diag) - 2 K_lin M_s", sp.simplify(Ls_diag - 2 * Klin_diag * Ms_diag))
expect_zero("L_t(diag) - 2 K_lin M_t", sp.simplify(Lt_diag - 2 * Klin_diag * Mt_diag))
expect_zero("L_u(diag) - 2 K_lin M_u", sp.simplify(Lu_diag - 2 * Klin_diag * Mu_diag))

expect_zero("M_r at gradient-optimal ratios", sp.simplify(Mr_diag.subs(grad_ratios)))
expect_zero("M_s at gradient-optimal ratios", sp.simplify(Ms_diag.subs(grad_ratios)))
expect_zero("M_t at gradient-optimal ratios", sp.simplify(Mt_diag.subs(grad_ratios)))
expect_zero("M_u at gradient-optimal ratios", sp.simplify(Mu_diag.subs(grad_ratios)))
expect_zero("L_r at gradient-optimal ratios", sp.simplify(Lr_diag.subs(grad_ratios)))
expect_zero("L_s at gradient-optimal ratios", sp.simplify(Ls_diag.subs(grad_ratios)))
expect_zero("L_t at gradient-optimal ratios", sp.simplify(Lt_diag.subs(grad_ratios)))
expect_zero("L_u at gradient-optimal ratios", sp.simplify(Lu_diag.subs(grad_ratios)))

# ---------------------------------------------------------------------------
# VI. Special reduction 2: full fivefold symmetry
# ---------------------------------------------------------------------------
subbanner("VI. Special reduction 2 — full fivefold symmetry")

k, nu_d, nu_o = sp.symbols("k nu_d nu_o", positive=True, real=True)
sym_subs = {
    kL: k,
    kc: k,
    kg: k,
    kU: k,
    kW: k,
    A: k**2 - 2 * H0 * nu_d,
    B: 2 * k**2 - 4 * H0 * nu_o,
    C: 2 * k**2 - 4 * H0 * nu_o,
    D: 2 * k**2 - 4 * H0 * nu_o,
    E: 2 * k**2 - 4 * H0 * nu_o,
    F: k**2 - 2 * H0 * nu_d,
    G: 2 * k**2 - 4 * H0 * nu_o,
    H: 2 * k**2 - 4 * H0 * nu_o,
    I: 2 * k**2 - 4 * H0 * nu_o,
    J: k**2 - 2 * H0 * nu_d,
    K: 2 * k**2 - 4 * H0 * nu_o,
    L: 2 * k**2 - 4 * H0 * nu_o,
    M: k**2 - 2 * H0 * nu_d,
    N: 2 * k**2 - 4 * H0 * nu_o,
    O: k**2 - 2 * H0 * nu_d,
}

bary_subs = {r: 1, s: 1, t: 1, u: 1}

expect_zero("M_r at equal-mix barycenter", sp.simplify(Mr.subs(sym_subs).subs(bary_subs)))
expect_zero("M_s at equal-mix barycenter", sp.simplify(Ms.subs(sym_subs).subs(bary_subs)))
expect_zero("M_t at equal-mix barycenter", sp.simplify(Mt.subs(sym_subs).subs(bary_subs)))
expect_zero("M_u at equal-mix barycenter", sp.simplify(Mu.subs(sym_subs).subs(bary_subs)))
expect_zero("L_r at equal-mix barycenter", sp.simplify(Lr.subs(sym_subs).subs(bary_subs)))
expect_zero("L_s at equal-mix barycenter", sp.simplify(Ls.subs(sym_subs).subs(bary_subs)))
expect_zero("L_t at equal-mix barycenter", sp.simplify(Lt.subs(sym_subs).subs(bary_subs)))
expect_zero("L_u at equal-mix barycenter", sp.simplify(Lu.subs(sym_subs).subs(bary_subs)))

subbanner("VII. Stage 217 audit verdict")
print("Verified:")
print("1. The unique support-cardinality-5 interior chart carries exact stationary numerators.")
print("2. The lifted stationary system has degree pattern (3,3,3,3,2) and Bezout bound 162.")
print("3. The direct projected one-chart elimination has degrees (5,5,5,6) and bound 750.")
print("4. The gradient-optimal interior ray is exact under diagonal-isotropic curvature.")
print("5. The equal-mix barycenter is exact under full fivefold symmetry.")
print("6. The unique five-coordinate interior problem is therefore a finite algebraic candidate problem.")
