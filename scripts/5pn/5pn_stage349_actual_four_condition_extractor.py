#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage346_349_common import *

banner("STAGE 349 — ACTUAL FOUR-CONDITION EXTRACTOR")

chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
epsW, epsEta = sp.symbols("epsilon_W epsilon_eta", positive=True, real=True)
ZW, Lambda = sp.symbols("Z_W Lambda", positive=True, real=True)
chiQ = sp.symbols("chi_Q", positive=True, real=True)

u_chi0, u_deltaU = sp.symbols("dlnchi_0 dlndelta_U", real=True)
u_epsW, u_epsEta = sp.symbols("dlnepsilon_W dlnepsilon_eta", real=True)
u_ZW, u_Lambda = sp.symbols("dlnZ_W dlnLambda", real=True)

eps = sp.simplify(epsW * (9*deltaU + 11) / (11*(1 + deltaU)))
Rtr = sp.simplify((1 + chi0/(1 + deltaU)) / (1 + chi0))
Rtarget = sp.simplify(Lambda * (1 - epsEta) * (1 - eps)**2 / (ZW * (1 + chi0)**2))
NQ = sp.simplify(1/chiQ)

dln_eps = sp.simplify(
    linearized_log(
        eps,
        {epsW: u_epsW, deltaU: u_deltaU},
    )
)
C1 = sp.simplify(
    linearized_log(
        Rtr,
        {chi0: u_chi0, deltaU: u_deltaU},
    )
)
C2 = sp.simplify(
    linearized_log(
        Rtarget,
        {
            chi0: u_chi0,
            deltaU: u_deltaU,
            epsW: u_epsW,
            epsEta: u_epsEta,
            ZW: u_ZW,
            Lambda: u_Lambda,
        },
    )
)
C3 = u_epsEta
C4 = sp.simplify(NQ - 1)

u_deltaU_req = sp.solve(sp.Eq(C1, 0), u_deltaU)[0]
u_Lambda_req = sp.solve(sp.Eq(C2.subs({u_deltaU: u_deltaU_req, u_epsEta: 0}), 0), u_Lambda)[0]

subbanner("I. Actual branch packet")
print("R_tr =")
sp.pprint(Rtr)
print("R_target =")
sp.pprint(Rtarget)
print("epsilon_eta =")
sp.pprint(epsEta)
print("N_Q =")
sp.pprint(NQ)

subbanner("II. Four finish-line conditionals")
print("C1 = d ln R_tr")
sp.pprint(C1)
print("C2 = d ln R_target")
sp.pprint(C2)
print("C3 = d ln epsilon_eta")
sp.pprint(C3)
print("C4 = N_Q - 1")
sp.pprint(C4)

subbanner("III. Exact landing compiler")
print("Required codrifts:")
sp.pprint(sp.Eq(u_deltaU, u_deltaU_req))
sp.pprint(sp.Eq(u_Lambda, u_Lambda_req))
print("and")
sp.pprint(sp.Eq(u_epsEta, 0))
print("plus outgoing normalization")
sp.pprint(sp.Eq(chiQ, 1))

expect_zero("C1 on actual orbit-lock surface", C1.subs(u_deltaU, u_deltaU_req))
expect_zero(
    "C2 on actual target-ratio surface",
    C2.subs({u_deltaU: u_deltaU_req, u_epsEta: 0, u_Lambda: u_Lambda_req}),
)
expect_zero("C3 on epsilon_eta-invariant surface", C3.subs(u_epsEta, 0))
expect_zero("C4 on canonical outgoing branch", C4.subs(chiQ, 1))

print("\nVerdict:")
print("  The current file stack fixes the actual coherent-branch values only symbolically.")
print("  What it gives exactly is the four-condition landing surface in the physical branch")
print("  variables. No numerical PDE-selected point is present yet in the notes/scripts.")
