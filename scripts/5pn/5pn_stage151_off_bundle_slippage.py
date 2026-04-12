#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage151_160_common import *

banner("STAGE 151 — EXACT OFF-BUNDLE SLIPPAGE DECOMPOSITION")

# Symbols from the Stage-147 normal coordinate.
dlnZq_rhow, dlncsw, dlncs, dlnTm, dlnvw0, dlna, dlnLW = sp.symbols(
    'dlnZq_rhow dlncsw dlncs dlnTm dlnvw0 dlna dlnLW', real=True
)
epsL, epsv, epsT = sp.symbols('eps_L eps_v eps_T', real=True)
gstar, rstar = sp.symbols('gstar rstar', positive=True, real=True)

A = sp.simplify(gstar + 1/(4*root1p(rstar)))
B = sp.simplify(1/(2*root1p(rstar)))
C = sp.simplify(2*gstar + 3/(4*root1p(rstar)))

subbanner("1. Carry-forward Stage-147 normal coordinate")
delta_perp = sp.simplify(
    A*dlnZq_rhow
    + 3*A*dlncsw
    + B*dlncs
    - gstar*dlnTm
    - (gstar + B)*dlnvw0
    - 2*A*dlna
    - C*dlnLW
)
print("delta_perp =")
sp.pprint(delta_perp)

subbanner("2. Exact lower-branch transport laws")
dlnvw0_br = sp.simplify(sp.Rational(1,2)*dlnZq_rhow + sp.Rational(3,2)*dlncsw + dlncs - sp.Rational(5,2)*dlna)
dlnTm_br  = sp.simplify(sp.Rational(1,2)*dlnZq_rhow + sp.Rational(3,2)*dlncsw - dlncs - sp.Rational(3,2)*dlna)
print("delta ln v_w0 (branch) =", dlnvw0_br)
print("delta ln T_m  (branch) =", dlnTm_br)

subbanner("3. Define the three scalar slippages")
print("eps_L =", epsL)
print("eps_v =", epsv)
print("eps_T =", epsT)

subs151 = {
    dlnLW: dlna + epsL,
    dlnvw0: dlnvw0_br + epsv,
    dlnTm: dlnTm_br + epsT,
}
delta_perp_slip = sp.simplify(delta_perp.subs(subs151))
target = sp.simplify(-gstar*epsT - (gstar + B)*epsv - C*epsL)

expect_zero("delta_perp with slippages - target", delta_perp_slip - target)
print("delta_perp after slippage substitution =")
sp.pprint(delta_perp_slip)

subbanner("4. Weighted scalar slippage combination")
eps_perp = sp.simplify(gstar*epsT + (gstar + B)*epsv + C*epsL)
expect_zero("delta_perp + eps_perp", delta_perp_slip + eps_perp)
print("eps_perp =", eps_perp)

subbanner("5. Mouth-bias transport with off-bundle slippage")
dSigma0, dS = sp.symbols('dSigma0 dS', real=True)
Sigma0_can_sym, S_can_sym = sp.symbols('Sigma0_can_sym S_can_sym', positive=True, real=True)
deltaPi_tan = sp.simplify((1 - S_can_sym/4)*dSigma0 - Sigma0_can_sym*dS/4)
deltaPi = sp.simplify(deltaPi_tan - Sigma0_can_sym*S_can_sym*eps_perp/root1p(rstar))
print("deltaPi_tan =")
sp.pprint(deltaPi_tan)
print("deltaPi =")
sp.pprint(deltaPi)

subbanner("6. Numerical Family-1 coefficients")
numsubs = {gstar: gF1, rstar: rF1, Sigma0_can_sym: Sigma0_can, S_can_sym: S_can}
print("A_* =", sp.N(A.subs(numsubs), 15))
print("B_* =", sp.N(B.subs(numsubs), 15))
print("C_* =", sp.N(C.subs(numsubs), 15))
print("delta_perp ≈")
sp.pprint(sp.N(target.subs(numsubs), 16))
print("deltaPi_tan coefficients:")
print("  dSigma0 =", sp.N((1 - S_can/4), 16))
print("  dS      =", sp.N(-Sigma0_can/4, 16))
print("  eps_perp coefficient =", sp.N(-(Sigma0_can*S_can/root1p(rF1)), 16))

banner("STAGE 151 LEDGER")
print("1. The full first-order off-bundle defect collapses to three scalar slippages:")
print("      eps_L, eps_v, eps_T.")
print("2. They enter the normal coordinate only through one weighted scalar")
print("      eps_perp = g_* eps_T + (g_*+B_*) eps_v + C_* eps_L.")
print("3. Exactly, delta_perp = - eps_perp.")
print("4. The same scalar controls the non-tangential mouth-bias defect.")
print("5. So the next microscopic theorem gate is to determine the actual sources of")
print("   eps_L, eps_v, eps_T rather than reopening the whole isotropic bundle.")
