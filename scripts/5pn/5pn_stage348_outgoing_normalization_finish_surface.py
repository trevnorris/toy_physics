#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage346_349_common import *

banner("STAGE 348 — OUTGOING NORMALIZATION FINISH SURFACE")

chiQ = sp.symbols("chi_Q", positive=True, real=True)
sigma_star = sp.symbols("sigma_star", positive=True, real=True)
Xi_slip, dPi_tan = sp.symbols("Xi_slip deltaPi_tan", real=True)

DeltaQ = sp.simplify(-sigma_star/(1 - sigma_star) * Xi_slip * dPi_tan)
NQ = sp.simplify(1/(1 + DeltaQ))
NQ_from_chiQ = sp.simplify(1/chiQ)

subbanner("I. Exact natural-source-map outgoing surface")
print("N_Q = 1 / chi_Q =")
sp.pprint(NQ_from_chiQ)
expect_zero("N_Q - 1 on canonical outgoing branch", (NQ_from_chiQ - 1).subs(chiQ, 1))

subbanner("II. Exact lower-parent compensation family at first order")
print("Delta_Q =")
sp.pprint(DeltaQ)
print("N_Q =")
sp.pprint(NQ)
expect_zero("N_Q - 1 when Xi_slip = 0", (NQ - 1).subs(Xi_slip, 0))
expect_zero("N_Q - 1 when deltaPi_tan = 0", (NQ - 1).subs(dPi_tan, 0))

print("\nInterpretation:")
print("  The exact outgoing finish surface is chi_Q = 1 on the natural source-map branch.")
print("  On the linearized lower-parent compensation family, Xi_slip * deltaPi_tan = 0")
print("  is sufficient to land on N_Q = 1.")
