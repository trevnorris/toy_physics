#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(t):
    line="="*88
    print("\n"+line); print(t); print(line)

def subbanner(t):
    line="-"*88
    print("\n"+line); print(t); print(line)

def expect_zero(name, expr):
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")

def expect_close(name, value, target, tol):
    res = abs(float(sp.N(value, 30) - sp.N(target, 30)))
    print(f"{name} residual = {res}")
    if res > tol:
        raise AssertionError(f"{name} off by {res} > tol {tol}")

Pi = sp.symbols("Pi", positive=True, real=True)
pi = sp.pi

banner("STAGE 142 — SELF-CONSISTENT MOUTH-BRANCH LAW")

# r_F1: Family-1 reduced mixed-core ratio. Carried forward from upstream
# (see notes/stages/moving_throat_pde_stage142_selfconsistent_mouth_branch.md
# section 1; original derivation is in the upstream "shell/mixed core" block
# referenced by paper/stages/stage_142.tex Inputs field).
r = sp.sqrt(sp.Integer(4107) - 100*pi**2)/(10*pi)
gPi = 2*Pi*(2*Pi*sp.exp(Pi)+pi)/((4*Pi**2+pi**2)*(sp.exp(Pi)-1))
# S_q(Pi) closed form: carried forward from the self-matched mouth-susceptibility
# closure (Stage 140 / Sigma_0 = (20/9) That_m^2). The closed form here is
# S(Pi, pi/2), evaluated at the fixed second argument pi/2.
Sq = Pi*(((pi/2)*sp.tanh(pi/2)) + Pi*(sp.exp(-Pi)*sp.sech(pi/2)-1))/((1-sp.exp(-Pi))*(((pi/2)**2)-Pi**2))
Rq = sp.simplify((gPi-r)**2/(1+r**2))
Sigma0 = Pi/(1-Rq*Sq)
That = sp.sqrt(sp.Rational(9,20)*Sigma0)

subbanner("Core-to-mouth reduction")
print("r_F1 =", r)
print("g_Pi =", gPi)
print("S_q(Pi) =", Sq)
print("R_q(Pi) =", Rq)
print("Sigma0(Pi) =", Sigma0)
print("That(Pi) =", That)

gminus = sp.simplify(r - sp.Rational(1,2)*sp.sqrt(1+r**2))
Rq_minus = sp.simplify(((gminus-r)**2)/(1+r**2))
expect_zero("R_q(g_minus)-1/4", Rq_minus - sp.Rational(1,4))

subbanner("Canonical compensation point")
Pi_star = sp.N(sp.nsolve(sp.Eq(gPi, gminus), 1.5), 30)
g_star = sp.N(gPi.subs(Pi, Pi_star), 30)
Rq_star = sp.N(Rq.subs(Pi, Pi_star), 30)
Sq_star = sp.N(Sq.subs(Pi, Pi_star), 30)
Sigma_star = sp.N(Sigma0.subs(Pi, Pi_star), 30)
That_star = sp.N(That.subs(Pi, Pi_star), 30)

print("g_minus^F1 =", sp.N(gminus, 30))
print("Pi_*       =", Pi_star)
print("g(Pi_*)    =", g_star)
print("R_q(Pi_*)  =", Rq_star)
print("S_q(Pi_*)  =", Sq_star)
print("Sigma0(Pi_*) =", Sigma_star)
print("That(Pi_*)   =", That_star)

if abs(float(g_star - sp.N(gminus,30))) > 1e-12:
    raise AssertionError("Pi_* does not solve the compensation equation accurately enough.")

# (R1, solver-consistency) R_q(Pi_*) = 1/4 holds IDENTICALLY at the solved
# point because Pi_* was found from gPi(Pi_*) = g_-, and R_q(g_-) = 1/4 for any
# r. So this is ONLY a redundant solver-residual check, NOT a paper-claim test.
Rq_star_residual = abs(float(Rq_star - sp.Rational(1,4)))
print(f"R_q(Pi_*) - 1/4 (solver-consistency) = {Rq_star - sp.Rational(1,4)}")
# 1e-15 tracks nsolve's actual precision here (residual ~1.945e-18); do NOT
# tighten to 1e-20 (that was too tight and was already loosened in a prior batch).
if Rq_star_residual > 1e-15:
    raise AssertionError(f"R_q(Pi_*) does not equal 1/4 at nsolve'd Pi_* (residual {Rq_star_residual}).")

# (R1, NON-tautological anchor) Evaluate R_q at STAGE 131's INDEPENDENTLY-DERIVED Pi_*
# (NOT 142's own nsolve output, NOT the self-solved point). Stage 131 found this parent
# mouth-threshold bias by a structurally different route (cleared-denominator FindRoot,
# batch-4-verified); 142's gPi crosses g_- there ONLY IF the hardcoded gPi closed form is
# right, so R_q lands on 1/4. A typo in gPi or a wrong r breaks this.
Pi_ext = sp.Float("1.50882951349315558300555075595", 30)  # Stage 131 Pi_* (independent)
Rq_at_ext = sp.N(Rq.subs(Pi, Pi_ext), 30)
Rq_ext_residual = abs(float(Rq_at_ext - sp.Rational(1,4)))
print(f"R_q(Pi_ext) - 1/4 (independent anchor) = {Rq_at_ext - sp.Rational(1,4)}")
if Rq_ext_residual > 1e-12:
    raise AssertionError(f"R_q at external Pi_*={Pi_ext} is not 1/4 (residual {Rq_ext_residual}); gPi or r is off.")

expect_close("g_-^{F1} value", gminus, sp.Float("0.7580350789446628269196808904", 30), 1e-25)
expect_close("Pi_* value",      Pi_star, sp.Float("1.5088295134931555274704351177", 30), 1e-12)
expect_close("S_q(Pi_*) value", Sq_star, sp.Float("0.6580759376054292719303153134", 30), 1e-12)
expect_close("Sigma_0(Pi_*) value", Sigma_star, sp.Float("1.8059411109563538072179672471", 30), 1e-12)
expect_close("That(Pi_*) value", That_star, sp.Float("0.9014840541742040227024016887", 30), 1e-12)

banner("STAGE 142 LEDGER")
print("Self-consistent Family-1 mouth branch:")
print("  Pi = Sigma0 * [1 - R_q(Pi) S_q(Pi)]")
print("  Sigma0(Pi) = Pi / (1 - R_q(Pi) S_q(Pi))")
print("  That(Pi)   = sqrt(9 Sigma0(Pi) / 20)")
print()
print("Canonical point:")
print(f"  Pi_*       = {Pi_star}")
print(f"  Sigma0_*   = {Sigma_star}")
print(f"  That_*     = {That_star}")
