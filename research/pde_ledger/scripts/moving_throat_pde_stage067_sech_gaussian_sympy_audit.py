
#!/usr/bin/env python3
"""
moving_throat_pde_stage50_sech_gaussian_sympy_audit.py

SymPy-backed audit for the explicit sech–Gaussian source/support coherence benchmark.

What this script checks
-----------------------
1. Exact transverse norms for the sech and Gaussian profiles.
2. Exact algebraic form of the coherence factor in terms of the dimensionless overlap I(r).
3. Exact self-duality implication:
      I(r) = (r/sqrt(pi)) I(pi/r)  =>  C^2(r) = C^2(pi/r).
4. Exact stationary-point consequence at r = sqrt(pi).
5. Numerical benchmark values:
      r_* = sqrt(pi),
      C_res^2 = C^2(r_*),
      P_res = 1 / C_res^2.
6. Numerical monotonicity scan on the constructive branch.

This is a mixed symbolic/numeric referee audit:
- SymPy handles all exact norm/duality algebra.
- mpmath handles the non-elementary overlap integral numerics.
"""

from __future__ import annotations

import mpmath as mp
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


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 50 — EXACT SECH–GAUSSIAN COHERENCE BENCHMARK")

# ---------------------------------------------------------------------------
# 1. Exact norm algebra
# ---------------------------------------------------------------------------

subbanner("1. Exact norms and dimensionless coherence form")

wf, wg, r = sp.symbols("w_f w_g r", positive=True, real=True)
Nss = 2 * wf
Npp = wg * sp.sqrt(sp.pi / 2)

print("N_(sigma sigma) =", Nss)
print("N_(phi phi)     =", Npp)

# Derive the norms by direct integration to anchor the declared values.
y = sp.symbols("y", real=True)
Nss_integral = sp.integrate((sp.sech(y / wf) ** 2).rewrite(sp.cosh), (y, -sp.oo, sp.oo))
Npp_integral = sp.integrate(sp.exp(-2 * y ** 2 / wg ** 2), (y, -sp.oo, sp.oo))
print("integrate(sech(y/w_f)^2)        =", sp.simplify(Nss_integral))
print("integrate(exp(-2 y^2/w_g^2))    =", sp.simplify(Npp_integral))
expect_zero("N_(sigma sigma) integral - 2 w_f", Nss_integral - Nss)
expect_zero("N_(phi phi) integral - w_g sqrt(pi/2)", Npp_integral - Npp)

# Dimensionless overlap variable.
I = sp.Function("I")
C2 = sp.simplify(I(r)**2 / (r * sp.sqrt(2 * sp.pi)))
print("C^2(r) =", C2)

# ---------------------------------------------------------------------------
# 2. Exact duality algebra
# ---------------------------------------------------------------------------

subbanner("2. Exact duality implication")

duality_rhs = (r / sp.sqrt(sp.pi)) * I(sp.pi / r)
C2_dual = sp.simplify(duality_rhs**2 / (r * sp.sqrt(2 * sp.pi)))
C2_target = sp.simplify(I(sp.pi / r)**2 / ((sp.pi / r) * sp.sqrt(2 * sp.pi)))
# Algebraic implication only: substitutes I -> (r/sqrt(pi)) I(pi/r) into C^2(r) and
# checks it equals C^2(pi/r). Holds for ANY function I; the duality identity for the
# sech-Gaussian overlap is exercised numerically in section 5.
expect_zero("C^2(r) - C^2(pi/r) under duality", C2_dual - C2_target)

print("Therefore: if I(r) = (r/sqrt(pi)) I(pi/r), then C^2(r) = C^2(pi/r).")

# ---------------------------------------------------------------------------
# 3. Exact stationary-point consequence
# ---------------------------------------------------------------------------

subbanner("3. Exact stationary point at the self-dual ratio")

rstar = sp.sqrt(sp.pi)
Istar, Iprime_left, Iprime_dual = sp.symbols(
    "I_star Iprime_left Iprime_dual",
    real=True,
)
duality_tangent = sp.simplify(
    Iprime_left - (Istar / sp.sqrt(sp.pi) - sp.sqrt(sp.pi) * Iprime_dual / r)
)
duality_tangent_at_rstar = sp.simplify(
    duality_tangent.subs(r, rstar).subs(Iprime_dual, Iprime_left)
)
print("differentiated overlap duality at r_* =", duality_tangent_at_rstar)
# Tautological: the substitution Iprime_left -> Istar/(2*sqrt(pi)) is the solution of
# the preceding equation. This checks calculus, not the sech-Gaussian profile.
expect_zero(
    "self-dual overlap-slope relation",
    duality_tangent_at_rstar.subs(Iprime_left, Istar / (2 * sp.sqrt(sp.pi))),
)

dC2_selfdual = sp.simplify(
    (2 * Istar * Iprime_left * rstar - Istar**2)
    / (rstar**2 * sp.sqrt(2 * sp.pi))
)
# Tautological: dC2_selfdual is a hand-written derivative formula, then the slope value
# derived above is substituted back. Stationarity of a symmetric differentiable function
# at the symmetric point is a calculus identity, not specific to sech-Gaussian.
# The substantive stationary-point evidence is the numerical monotonicity scan below.
expect_zero(
    "stationary derivative of C^2 at the self-dual point",
    dC2_selfdual.subs(Iprime_left, Istar / (2 * sp.sqrt(sp.pi))),
)

delta_bad = sp.symbols("delta_bad", nonzero=True, real=True)
dC2_broken = sp.simplify(
    dC2_selfdual.subs(Iprime_left, Istar / (2 * sp.sqrt(sp.pi)) + delta_bad)
)
print("broken self-dual derivative =", dC2_broken)
if dC2_broken == 0:
    raise AssertionError("Expected a perturbed self-dual tangent to break stationarity.")

print("Hence the self-dual point r_* = sqrt(pi) is an exact stationary point.")

# ---------------------------------------------------------------------------
# 4. Numerical overlap benchmark
# ---------------------------------------------------------------------------

subbanner("4. Numerical evaluation of the non-elementary overlap")

mp.mp.dps = 60

def I_num(rr: mp.mpf) -> mp.mpf:
    f = lambda x: mp.sech(x) * mp.e**(-(x**2) / (rr**2))
    return mp.quad(f, [-mp.inf, mp.inf])

def C2_num(rr: mp.mpf) -> mp.mpf:
    val = I_num(rr)
    return (val**2) / (rr * mp.sqrt(2 * mp.pi))

rstar_num = mp.sqrt(mp.pi)
C2_star = C2_num(rstar_num)
Pres = 1 / C2_star

print("r_* =", rstar_num)
print("C_res^2 =", C2_star)
print("P_res   =", Pres)
print("1 - C_res^2 =", 1 - C2_star)

# ---------------------------------------------------------------------------
# 5. Numerical duality checks
# ---------------------------------------------------------------------------

subbanner("5. Numerical duality checks")

samples = [mp.mpf("0.75"), mp.mpf("1.0"), mp.mpf("1.2"), mp.mpf("1.5"), mp.mpf("2.0")]
for rr in samples:
    lhs = I_num(rr)
    rhs = (rr / mp.sqrt(mp.pi)) * I_num(mp.pi / rr)
    diff = abs(lhs - rhs)
    print(f"r = {rr}:  |I(r) - r/sqrt(pi) I(pi/r)| = {diff}")
    if diff > mp.mpf("1e-40"):
        raise AssertionError("Numerical duality check failed.")

# ---------------------------------------------------------------------------
# 6. Numerical monotonicity scan
# ---------------------------------------------------------------------------

subbanner("6. Numerical monotonicity scan")

left_grid = [mp.mpf("0.55"), mp.mpf("0.7"), mp.mpf("0.9"), mp.mpf("1.1"), mp.mpf("1.35"), mp.mpf("1.6"), rstar_num]
vals_left = [C2_num(x) for x in left_grid]
for x, y in zip(left_grid, vals_left):
    print(f"C^2({x}) = {y}")

# Strict increase up to the self-dual point in the sampled constructive branch.
for a, b in zip(vals_left[:-1], vals_left[1:]):
    if not (b > a):
        raise AssertionError("Constructive-branch monotonicity failed before the self-dual point.")

right_grid = [rstar_num, mp.mpf("2.0"), mp.mpf("2.5"), mp.mpf("3.0"), mp.mpf("4.0")]
vals_right = [C2_num(x) for x in right_grid]
for x, y in zip(right_grid, vals_right):
    print(f"C^2({x}) = {y}")

for a, b in zip(vals_right[:-1], vals_right[1:]):
    if not (b < a):
        raise AssertionError("Constructive-branch monotonicity failed after the self-dual point.")

banner("FINAL LEDGER")
print("Exact symbolic results:")
print("  N_(sigma sigma) = 2 w_f")
print("  N_(phi phi)     = w_g sqrt(pi/2)")
print("  C^2(r)          = I(r)^2 / [ r sqrt(2 pi) ]")
print("  duality         = I(r) = (r/sqrt(pi)) I(pi/r)")
print("  symmetry        = C^2(r) = C^2(pi/r)")
print("  stationary point= r_* = sqrt(pi)")
print()
print("Numerical benchmark:")
print(f"  C_res^2 = {C2_star}")
print(f"  P_res   = {Pres}")
print("Interpretation:")
print("  The explicit independent sech–Gaussian family nearly saturates the exact")
print("  matched-layer ideal, but it does not by itself prove support-threshold survival.")
