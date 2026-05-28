#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(t):
    line="="*88
    print("\n"+line); print(t); print(line)

def subbanner(t):
    line="-"*88
    print("\n"+line); print(t); print(line)

Pi = sp.symbols("Pi", positive=True, real=True)
pi = sp.pi

banner("STAGE 144 — UNIQUE REGULAR CANONICAL MOUTH BRANCH")

r = sp.sqrt(sp.Integer(4107) - 100*pi**2)/(10*pi)
gminus = sp.simplify(r - sp.Rational(1,2)*sp.sqrt(1+r**2))
gplus  = sp.simplify(r + sp.Rational(1,2)*sp.sqrt(1+r**2))
gPi = 2*Pi*(2*Pi*sp.exp(Pi)+pi)/((4*Pi**2+pi**2)*(sp.exp(Pi)-1))
Sq = Pi*(((pi/2)*sp.tanh(pi/2)) + Pi*(sp.exp(-Pi)*sp.sech(pi/2)-1))/((1-sp.exp(-Pi))*(((pi/2)**2)-Pi**2))
Rq = (gPi-r)**2/(1+r**2)
Sigma0 = Pi/(1-Rq*Sq)
That = sp.sqrt(sp.Rational(9,20)*Sigma0)

subbanner("Compensated branches")
print("g_-^F1 =", sp.N(gminus, 30))
print("g_+^F1 =", sp.N(gplus, 30))
print("2/pi   =", sp.N(2/pi, 30))

Pi_star = sp.N(sp.nsolve(sp.Eq(gPi, gminus), 1.5), 30)
That_star = sp.N(That.subs(Pi, Pi_star), 30)

Pi_match = sp.N(sp.nsolve(sp.Eq(gPi, pi/4), 1.9), 30)
That_match = sp.N(That.subs(Pi, Pi_match), 30)

print("Pi_*          =", Pi_star)
print("That(Pi_*)    =", That_star)
Sigma0_star = sp.N(Sigma0.subs(Pi, Pi_star), 30)
Sigma0_match = sp.N(Sigma0.subs(Pi, Pi_match), 30)
print("Sigma0(Pi_*)  =", Sigma0_star)
print("Sigma0(Pi_match)=", Sigma0_match)
print("Pi_match      =", Pi_match)
print("That(Pi_match)=", That_match)

if not (float(Pi_star) > 0 and float(Pi_match) > float(Pi_star)):
    raise AssertionError("Unexpected ordering of canonical and matched-derivative points.")
# Numerical targets from notes/stages/moving_throat_pde_stage144_unique_regular_canonical_branch.md
gminus_N = float(sp.N(gminus, 30))
gplus_N  = float(sp.N(gplus, 30))
two_over_pi_N = float(sp.N(2/pi, 30))
if not (gplus_N > 1):
    raise AssertionError(f"upper branch must satisfy g_+^F1 > 1, got {gplus_N}")
if not (two_over_pi_N < gminus_N < 1):
    raise AssertionError(f"lower branch must satisfy 2/pi < g_-^F1 < 1, got {gminus_N}")
TOL = sp.Rational(1, 10**12)
if not (abs(sp.N(Pi_star - sp.Rational('1.50882951349316'), 30)) < TOL):
    raise AssertionError(f"Pi_* drift: got {Pi_star}, expected 1.50882951349316")
if not (abs(sp.N(That_star - sp.Rational('0.901484054174205'), 30)) < TOL):
    raise AssertionError(f"That(Pi_*) drift: got {That_star}, expected 0.901484054174205")
if not (abs(sp.N(Sigma0_star - sp.Rational('1.80594111095636'), 30)) < TOL):
    raise AssertionError(f"Sigma0(Pi_*) drift: got {Sigma0_star}, expected 1.80594111095636")
if not (abs(sp.N(Pi_match - sp.Rational('1.90848600654854'), 30)) < TOL):
    raise AssertionError(f"Pi_match drift: got {Pi_match}, expected 1.90848600654854")
if not (abs(sp.N(That_match - sp.Rational('1.01132972803599'), 30)) < TOL):
    raise AssertionError(f"That(Pi_match) drift: got {That_match}, expected 1.01132972803599")
print("PASS: numerical-target assertions (upper/lower branches, Pi_*, Sigma0_*, That_*, Pi_match, That_match)")

banner("STAGE 144 LEDGER")
print("Positive-source theorem gives 0 <= g <= 1, so the upper compensated branch is impossible.")
print("The exponential mouth family is monotone and spans (2/pi, 1).")
print("Since g_-^F1 lies strictly inside that interval, it is reached at one unique finite Pi_*.")
print("That point has finite moderate traction, whereas g=1 is only the singular Pi->infty limit.")
print()
print(f"Pi_*       = {Pi_star}")
print(f"That_*     = {That_star}")
print(f"Pi_match   = {Pi_match}")
print(f"That_match = {That_match}")
