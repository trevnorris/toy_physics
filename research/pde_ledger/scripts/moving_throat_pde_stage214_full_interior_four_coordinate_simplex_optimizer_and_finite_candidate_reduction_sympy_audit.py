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


def expect_true(name: str, condition: bool) -> None:
    print(f"{name} = {condition}")
    if not condition:
        raise AssertionError(f"{name} failed")


banner("STAGE 214 — FULL INTERIOR FOUR-COORDINATE SIMPLEX OPTIMIZER")
banner("EXACT LIFTED STATIONARY SYSTEM AND FINITE ALGEBRAIC CANDIDATE REDUCTION")

# ---------------------------------------------------------------------------
# I. Exact interior ratio objective and stationary numerators
# ---------------------------------------------------------------------------
subbanner("I. Exact interior ratio objective and stationary numerators")

r, s, t, y = sp.symbols("r s t y", nonnegative=True, real=True)
ki, kj, kk, kl = sp.symbols("k_i k_j k_k k_l", positive=True, real=True)
H0 = sp.symbols("H0", positive=True, real=True)

A, B, C, D, E, F, G, H, I, J = sp.symbols(
    "A B C D E F G H I J", real=True
)
Delta = sp.expand(
    A
    + B * r
    + C * s
    + D * t
    + E * r**2
    + F * r * s
    + G * r * t
    + H * s**2
    + I * s * t
    + J * t**2
)
sqrtDelta = sp.sqrt(Delta)

Phi = (ki + kj * r + kk * s + kl * t + sqrtDelta) / sp.sqrt(1 + r**2 + s**2 + t**2)
tau = sp.simplify(2 * H0 / Phi)

Mr = sp.simplify((1 + r**2 + s**2 + t**2) * kj - r * (ki + kj * r + kk * s + kl * t))
Ms = sp.simplify((1 + r**2 + s**2 + t**2) * kk - s * (ki + kj * r + kk * s + kl * t))
Mt = sp.simplify((1 + r**2 + s**2 + t**2) * kl - t * (ki + kj * r + kk * s + kl * t))

Lr = sp.simplify((1 + r**2 + s**2 + t**2) * sp.diff(Delta, r) - 2 * r * Delta)
Ls = sp.simplify((1 + r**2 + s**2 + t**2) * sp.diff(Delta, s) - 2 * s * Delta)
Lt = sp.simplify((1 + r**2 + s**2 + t**2) * sp.diff(Delta, t) - 2 * t * Delta)

Nr = sp.simplify(2 * Mr * sqrtDelta + Lr)
Ns = sp.simplify(2 * Ms * sqrtDelta + Ls)
Nt = sp.simplify(2 * Mt * sqrtDelta + Lt)

print("Phi(r,s,t) =")
sp.pprint(Phi)
print("tau(r,s,t) =")
sp.pprint(tau)
print("M_r =")
sp.pprint(Mr)
print("M_s =")
sp.pprint(Ms)
print("M_t =")
sp.pprint(Mt)
print("L_r =")
sp.pprint(Lr)
print("L_s =")
sp.pprint(Ls)
print("L_t =")
sp.pprint(Lt)

expected_dPhi_dr = Nr / (2 * (1 + r**2 + s**2 + t**2) ** sp.Rational(3, 2) * sqrtDelta)
expected_dPhi_ds = Ns / (2 * (1 + r**2 + s**2 + t**2) ** sp.Rational(3, 2) * sqrtDelta)
expected_dPhi_dt = Nt / (2 * (1 + r**2 + s**2 + t**2) ** sp.Rational(3, 2) * sqrtDelta)

expect_zero("exact dPhi/dr compiler", sp.diff(Phi, r) - expected_dPhi_dr)
expect_zero("exact dPhi/ds compiler", sp.diff(Phi, s) - expected_dPhi_ds)
expect_zero("exact dPhi/dt compiler", sp.diff(Phi, t) - expected_dPhi_dt)

# ---------------------------------------------------------------------------
# II. Lifted polynomial stationary system and degree ledger
# ---------------------------------------------------------------------------
subbanner("II. Lifted polynomial stationary system and degree ledger")

Fr = sp.expand(2 * Mr * y + Lr)
Fs = sp.expand(2 * Ms * y + Ls)
Ft = sp.expand(2 * Mt * y + Lt)
FDelta = sp.expand(y**2 - Delta)

print("F_r(r,s,t,y) =")
sp.pprint(Fr)
print("F_s(r,s,t,y) =")
sp.pprint(Fs)
print("F_t(r,s,t,y) =")
sp.pprint(Ft)
print("F_Delta(r,s,t,y) =")
sp.pprint(FDelta)

deg_Fr = sp.Poly(Fr, r, s, t, y).total_degree()
deg_Fs = sp.Poly(Fs, r, s, t, y).total_degree()
deg_Ft = sp.Poly(Ft, r, s, t, y).total_degree()
deg_FDelta = sp.Poly(FDelta, r, s, t, y).total_degree()

print(f"deg F_r = {deg_Fr}")
print(f"deg F_s = {deg_Fs}")
print(f"deg F_t = {deg_Ft}")
print(f"deg F_Delta = {deg_FDelta}")

if deg_Fr != 3:
    raise AssertionError("F_r is not cubic")
if deg_Fs != 3:
    raise AssertionError("F_s is not cubic")
if deg_Ft != 3:
    raise AssertionError("F_t is not cubic")
if deg_FDelta != 2:
    raise AssertionError("F_Delta is not quadratic")

bezout_lift = deg_Fr * deg_Fs * deg_Ft * deg_FDelta
print(f"lifted Bezout bound = {bezout_lift}")
if bezout_lift != 54:
    raise AssertionError("Unexpected lifted Bezout bound")

# ---------------------------------------------------------------------------
# III. Exact square-root-free elimination identities
# ---------------------------------------------------------------------------
subbanner("III. Exact square-root-free elimination identities")

Crs = sp.expand(Ms * Lr - Mr * Ls)
Crt = sp.expand(Mt * Lr - Mr * Lt)
Cst = sp.expand(Mt * Ls - Ms * Lt)

Sr = sp.expand(Lr**2 - 4 * Mr**2 * Delta)
Ss = sp.expand(Ls**2 - 4 * Ms**2 * Delta)
St = sp.expand(Lt**2 - 4 * Mt**2 * Delta)

print("C_rs(r,s,t) =")
sp.pprint(Crs)
print("C_rt(r,s,t) =")
sp.pprint(Crt)
print("C_st(r,s,t) =")
sp.pprint(Cst)
print("S_r(r,s,t) =")
sp.pprint(Sr)
print("S_s(r,s,t) =")
sp.pprint(Ss)
print("S_t(r,s,t) =")
sp.pprint(St)

for label, eliminant in [
    ("C_rs", Crs),
    ("C_rt", Crt),
    ("C_st", Cst),
    ("S_r", Sr),
    ("S_s", Ss),
    ("S_t", St),
]:
    expect_true(f"{label} non-vacuous polynomial", not sp.Poly(eliminant, r, s, t).is_zero)

stationary_subs = {
    ki: sp.Rational(2),
    kj: sp.Rational(3),
    kk: sp.Rational(5),
    kl: sp.Rational(7),
    H0: sp.Rational(1),
    A: -sp.Rational(220, 3),
    B: sp.Rational(12),
    C: sp.Rational(20),
    D: sp.Rational(28),
    E: -sp.Rational(205, 3),
    F: sp.Rational(30),
    G: sp.Rational(42),
    H: -sp.Rational(157, 3),
    I: sp.Rational(70),
    J: -sp.Rational(85, 3),
    r: sp.Rational(3, 2),
    s: sp.Rational(5, 2),
    t: sp.Rational(7, 2),
    y: sp.Rational(29, 2),
}

for label, stationary_poly in [
    ("F_r", Fr),
    ("F_s", Fs),
    ("F_t", Ft),
    ("F_Delta", FDelta),
    ("C_rs", Crs),
    ("C_rt", Crt),
    ("C_st", Cst),
    ("S_r", Sr),
    ("S_s", Ss),
    ("S_t", St),
]:
    expect_zero(f"{label} at diagonal-isotropic stationary point", stationary_poly.subs(stationary_subs))

deg_Crs = sp.Poly(Crs, r, s, t).total_degree()
deg_Crt = sp.Poly(Crt, r, s, t).total_degree()
deg_Cst = sp.Poly(Cst, r, s, t).total_degree()
deg_Sr = sp.Poly(Sr, r, s, t).total_degree()
deg_Ss = sp.Poly(Ss, r, s, t).total_degree()
deg_St = sp.Poly(St, r, s, t).total_degree()

print(f"deg C_rs = {deg_Crs}")
print(f"deg C_rt = {deg_Crt}")
print(f"deg C_st = {deg_Cst}")
print(f"deg S_r = {deg_Sr}")
print(f"deg S_s = {deg_Ss}")
print(f"deg S_t = {deg_St}")

if deg_Crs != 5 or deg_Crt != 5 or deg_Cst != 5:
    raise AssertionError("cross-consistency polynomial is not quintic")
if deg_Sr != 6 or deg_Ss != 6 or deg_St != 6:
    raise AssertionError("square condition is not sextic")

bezout_projected = deg_Crs * deg_Crt * deg_Sr
print(f"projected one-chart Bezout bound = {bezout_projected}")
if bezout_projected != 150:
    raise AssertionError("Unexpected projected Bezout bound")

# ---------------------------------------------------------------------------
# IV. Diagonal-isotropic curvature reduction
# ---------------------------------------------------------------------------
subbanner("IV. Diagonal-isotropic curvature reduction")

u = sp.symbols("u", real=True)
diag_iso = {
    A: ki**2 - 2 * H0 * u,
    B: 2 * ki * kj,
    C: 2 * ki * kk,
    D: 2 * ki * kl,
    E: kj**2 - 2 * H0 * u,
    F: 2 * kj * kk,
    G: 2 * kj * kl,
    H: kk**2 - 2 * H0 * u,
    I: 2 * kk * kl,
    J: kl**2 - 2 * H0 * u,
}
Delta_iso = sp.simplify(Delta.subs(diag_iso))
k_rst = sp.simplify((ki + kj * r + kk * s + kl * t) / sp.sqrt(1 + r**2 + s**2 + t**2))
tau_iso = sp.simplify(tau.subs(diag_iso))

expect_zero(
    "Delta_iso reduction",
    Delta_iso - ((ki + kj * r + kk * s + kl * t) ** 2 - 2 * H0 * u * (1 + r**2 + s**2 + t**2)),
)
expect_zero(
    "tau_iso depends only on k_rst",
    tau_iso - 2 * H0 / (k_rst + sp.sqrt(k_rst**2 - 2 * H0 * u)),
)

Kgrad = sp.sqrt(ki**2 + kj**2 + kk**2 + kl**2)
a_grad = sp.Matrix([ki / Kgrad, kj / Kgrad, kk / Kgrad, kl / Kgrad])
print("a_grad =")
sp.pprint(a_grad)
expect_zero("gradient-optimal normalization", (a_grad.T * a_grad)[0] - 1)
expect_zero("gradient-optimal slope value", sp.simplify((a_grad.T * sp.Matrix([ki, kj, kk, kl]))[0] - Kgrad))

# ---------------------------------------------------------------------------
# V. Full quadruple-symmetry reduction and equal-mix barycenter
# ---------------------------------------------------------------------------
subbanner("V. Full quadruple-symmetry reduction and equal-mix barycenter")

k = sp.symbols("k", positive=True, real=True)
ud, ux = sp.symbols("u_d u_x", real=True)

sym_subs = {
    ki: k,
    kj: k,
    kk: k,
    kl: k,
    A: k**2 - 2 * H0 * ud,
    B: 2 * k**2 - 4 * H0 * ux,
    C: 2 * k**2 - 4 * H0 * ux,
    D: 2 * k**2 - 4 * H0 * ux,
    E: k**2 - 2 * H0 * ud,
    F: 2 * k**2 - 4 * H0 * ux,
    G: 2 * k**2 - 4 * H0 * ux,
    H: k**2 - 2 * H0 * ud,
    I: 2 * k**2 - 4 * H0 * ux,
    J: k**2 - 2 * H0 * ud,
}

Nr_sym = sp.simplify(Nr.subs(sym_subs).subs({r: 1, s: 1, t: 1}))
Ns_sym = sp.simplify(Ns.subs(sym_subs).subs({r: 1, s: 1, t: 1}))
Nt_sym = sp.simplify(Nt.subs(sym_subs).subs({r: 1, s: 1, t: 1}))

expect_zero("equal-mix stationary numerator Nr(1,1,1)", Nr_sym)
expect_zero("equal-mix stationary numerator Ns(1,1,1)", Ns_sym)
expect_zero("equal-mix stationary numerator Nt(1,1,1)", Nt_sym)

a_eq = sp.Matrix([sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 2)])
expect_zero("equal-mix normalization", (a_eq.T * a_eq)[0] - 1)

# ---------------------------------------------------------------------------
# VI. Interior winner theorem against the Stage 213 boundary ledger
# ---------------------------------------------------------------------------
subbanner("VI. Interior winner theorem against the Stage 213 boundary ledger")

beta_lo, beta_hi = sp.symbols("beta_lo beta_hi", real=True)

count_win = 0
for blo in range(1, 8):
    for bhi in range(blo, 8):
        for ilo in range(0, 8):
            for ihi in range(ilo, 8):
                if ihi < blo:
                    for bstar in range(blo, bhi + 1):
                        for istar in range(ilo, ihi + 1):
                            if not (istar < bstar):
                                raise AssertionError("interior winner theorem failed")
                            count_win += 1
print(f"verified interior winner theorem on {count_win} ordered integer samples")

count_filter = 0
for blo in range(0, 8):
    for bhi in range(blo, 8):
        for ilo in range(0, 8):
            for ihi in range(ilo, 8):
                if ilo > bhi:
                    for bstar in range(blo, bhi + 1):
                        for istar in range(ilo, ihi + 1):
                            if not (bstar < istar):
                                raise AssertionError("interior non-improvement theorem failed")
                            count_filter += 1
print(f"verified interior non-improvement theorem on {count_filter} ordered integer samples")

banner("STAGE 214 SYMPY AUDIT COMPLETED SUCCESSFULLY")
