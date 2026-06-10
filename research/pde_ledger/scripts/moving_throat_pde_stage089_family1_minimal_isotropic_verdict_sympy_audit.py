#!/usr/bin/env python3
"""
moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py

SymPy / arithmetic audit for Stage 089.

Checks:
1. compare rho_alpha = 4/3 against the explicit Family-1 ratio window;
2. compare zeta_req = 1/3 against the explicit support ceiling;
3. verify that zeta_req < A_F1, so the explicit Family-1 transport map requires Pe_req = 0.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr, tol: float = 1e-12) -> None:
    val = sp.N(sp.simplify(expr), 50)
    print(f"{name} = {val}")
    if abs(complex(val)) > tol:
        raise AssertionError(f"{name} is not within tolerance {tol}")

banner("STAGE 089 — EXPLICIT FAMILY-1 VERDICT FOR THE MINIMAL ISOTROPIC BRANCH")

rho_min = sp.Rational(4,3)
zeta_min = sp.Rational(1,3)

# Stage-079 Family-1 demand map data.
Pe = sp.symbols("Pe", positive=True, real=True)
y = sp.symbols("y", real=True)
kappa_F1 = sp.Rational(12321, 5)
eta_F1 = sp.Integer(37)
y_F1 = sp.nsolve(y * sp.tan(y) - eta_F1, sp.Float("1.53", 80), tol=1e-40, maxsteps=100, prec=100)
A_F1 = sp.simplify((kappa_F1 + sp.pi**2 / 4) / (kappa_F1 + y_F1**2))
Omega = sp.simplify(
    sp.pi * Pe * (2 * Pe * sp.exp(Pe) + sp.pi)
    / ((4 * Pe**2 + sp.pi**2) * (sp.exp(Pe) - 1))
)
zeta_F1 = sp.simplify(A_F1 * Omega**2)
zeta_max = sp.simplify(sp.limit(zeta_F1, Pe, sp.oo))

# Paper-side chain closure: the Stage 089 derivation hinges on the 0/0 limit
# Omega(Pe -> 0) = 1, which gives zeta_F1(0) = A_F1. Verify both explicitly
# (these complete the link from the precondition zeta_req^min < A_F1 down to
# the boxed Output Pe_req = 0). The Omega expression is genuinely 0/0 at
# Pe = 0; the symbolic limit applies l'Hopital and returns 1.
Omega_at_zero = sp.simplify(sp.limit(Omega, Pe, 0))
zeta_F1_at_zero = sp.simplify(sp.limit(zeta_F1, Pe, 0))
expect_zero("Omega(Pe -> 0) - 1", Omega_at_zero - 1, tol=1e-30)
expect_zero("zeta_F1(Pe -> 0) - A_F1", zeta_F1_at_zero - A_F1, tol=1e-30)

# Stage-082 (post-renumber) thresholds evaluated at lambda_mu = 1.
# CARRY-FORWARD: Pe_suff_chi and Pe_fail_chi are the loading-ratio Pe values
# that produce the upstream rho_suff^(chi) and rho_fail^(chi) thresholds.
# Source: scripts/output/moving_throat_pde_stage082_*_sympy_audit.txt and the
# Stage 089 notes section 1. The literal values are not rederived here to
# avoid sp.nsolve instability near the tan(y) singularity of the Stage 074
# closed form (see notes/STAGE_VERIFICATION_COVERAGE.md pitfall #10).
Pe_suff_chi = sp.Float("96.5285247264386")
Pe_fail_chi = sp.Float("11220.5441626259")
eps_blk = sp.symbols("eps_blk", positive=True, real=True)
zeta = sp.symbols("zeta", positive=True, real=True)
Q_gen = (1 + (1 - 2 * eps_blk) * zeta) / (1 - eps_blk * zeta)   # general blocked-module loading ratio Q(zeta;eps_blk)
Q = sp.simplify(Q_gen.subs(eps_blk, 0))                          # unblocked specialization used downstream
zeta_suff = sp.simplify(zeta_F1.subs(Pe, Pe_suff_chi))
zeta_fail = sp.simplify(zeta_F1.subs(Pe, Pe_fail_chi))
rho_suff = sp.simplify(Q.subs(zeta, zeta_suff))
rho_fail = sp.simplify(Q.subs(zeta, zeta_fail))
rho_max = sp.simplify(Q.subs(zeta, zeta_max))

expect_zero("Stage-080 zeta_max = A_F1 pi^2/4", zeta_max - A_F1 * sp.pi**2 / 4, tol=1e-30)
expect_zero("Stage-082 Q(zeta;0)=1+zeta reduction", Q_gen.subs(eps_blk, 0) - (1 + zeta), tol=1e-30)


def expect_close(name: str, value: sp.Expr, target: sp.Expr, tol: sp.Expr) -> None:
    diff = sp.Abs(sp.N(value - target, 40))
    print(f"{name} diff = {diff}")
    if diff > tol:
        raise AssertionError(f"{name} exceeds tolerance {tol}: diff={diff}")


# Cross-check rho_X against upstream Stage-082 quoted values. The previous
# `rho_X - (1 + zeta_X) == 0` form was tautological because Q(zeta; eps=0) =
# 1 + zeta is the algebraic structure of Q, not a check of the literals.
expect_close("rho_suff vs Stage-082 quote", rho_suff, sp.Float("3.46622291347846", 30), sp.Float("1e-12", 30))
expect_close("rho_fail vs Stage-082 quote", rho_fail, sp.Float("3.46752913273870", 30), sp.Float("1e-12", 30))
expect_close("rho_max  vs Stage-082 quote", rho_max,  sp.Float("3.46752922945601", 30), sp.Float("1e-12", 30))

Delta_suff = sp.N(rho_suff - rho_min, 25)
Delta_fail = sp.N(rho_fail - rho_min, 25)
Delta_max  = sp.N(rho_max  - rho_min, 25)
Delta_zeta = sp.N(zeta_max - zeta_min, 25)
Delta_AF1  = sp.N(A_F1 - zeta_min, 25)

print("rho_min   =", sp.N(rho_min, 25))
print("rho_suff  =", sp.N(rho_suff, 25))
print("rho_fail  =", sp.N(rho_fail, 25))
print("rho_max   =", sp.N(rho_max, 25))
print("zeta_min  =", sp.N(zeta_min, 25))
print("zeta_max  =", sp.N(zeta_max, 25))
print("A_F1      =", sp.N(A_F1, 25))

print("\nMargins:")
print("Delta_suff =", Delta_suff)
print("Delta_fail =", Delta_fail)
print("Delta_max  =", Delta_max)
print("Delta_zeta =", Delta_zeta)
print("Delta_AF1  =", Delta_AF1)

if not (rho_min < rho_suff < rho_fail < rho_max):
    raise AssertionError("Family-1 loading-ratio ordering failed.")
if not (zeta_min < 1):
    raise AssertionError("Minimal isotropic branch left the symmetric-lowest-twin regime.")
if not (zeta_min < A_F1):
    raise AssertionError("Minimal isotropic branch no longer succeeds at zero transport bias.")
if not (zeta_min < zeta_max):
    raise AssertionError("Minimal isotropic branch exceeded the Family-1 support ceiling.")

# Pe_req = 0 is FORCED (not assumed): zeta_F1(0) = A_F1 already exceeds the
# minimal isotropic demand zeta_min, so the demand is met at zero transport
# bias and the minimal required Peclet number is 0. The positive zero-bias
# success margin below is the can-fail assertion establishing the boxed Output
# Pe_req = 0 (paper/stages/stage_089.tex eq app-stage089-Pe-zero); if the
# margin were <= 0 the branch would need Pe_req > 0 and this check would raise.
zero_bias_margin = sp.N(zeta_F1_at_zero - zeta_min, 40)
print("zero-bias success margin zeta_F1(0) - zeta_min =", zero_bias_margin)
if not (zero_bias_margin > 0):
    raise AssertionError("zeta_F1(0) <= zeta_min: demand unmet at zero bias -> Pe_req != 0")
Pe_req = sp.Integer(0)   # forced by the positive zero-bias margin above

print("\nRegime checks:")
print("  rho_min < rho_suff   -> guaranteed success")
print("  zeta_min < 1         -> symmetric lowest twin already enough")
print("  zeta_min < A_F1      -> Pe_req = 0 on the explicit Family-1 transport map")
print(f"  Pe_req = {Pe_req}     -> paper Output app-stage089-Pe-zero verified")
