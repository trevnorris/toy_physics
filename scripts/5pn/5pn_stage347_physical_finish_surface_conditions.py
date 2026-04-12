#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage346_349_common import *

banner("STAGE 347 — PHYSICAL FINISH-SURFACE CONDITIONS IN ACTUAL BRANCH VARIABLES")

chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
epsW, epsEta = sp.symbols("epsilon_W epsilon_eta", positive=True, real=True)
ZW, Lambda = sp.symbols("Z_W Lambda", positive=True, real=True)

u_chi0, u_deltaU = sp.symbols("dlnchi_0 dlndelta_U", real=True)
u_epsW, u_epsEta = sp.symbols("dlnepsilon_W dlnepsilon_eta", real=True)
u_ZW, u_Lambda = sp.symbols("dlnZ_W dlnLambda", real=True)

eps = sp.simplify(epsW * (9*deltaU + 11) / (11*(1 + deltaU)))
Rtr = sp.simplify((1 + chi0/(1 + deltaU)) / (1 + chi0))
Rtarget = sp.simplify(Lambda * (1 - epsEta) * (1 - eps)**2 / (ZW * (1 + chi0)**2))

dln_eps = sp.simplify(
    linearized_log(
        eps,
        {
            epsW: u_epsW,
            deltaU: u_deltaU,
        },
    )
)
dln_Rtr = sp.simplify(
    linearized_log(
        Rtr,
        {
            chi0: u_chi0,
            deltaU: u_deltaU,
        },
    )
)
dln_Rtarget = sp.simplify(
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

subbanner("I. Exact physical drift laws")
print("d ln epsilon =")
sp.pprint(dln_eps)
print("d ln R_tr =")
sp.pprint(dln_Rtr)
print("d ln R_target =")
sp.pprint(dln_Rtarget)
print("d ln epsilon_eta =")
sp.pprint(u_epsEta)

subbanner("II. The actual finish surfaces")
track_surface = sp.solve(sp.Eq(dln_Rtr, 0), u_deltaU)[0]
rtarget_surface = sp.solve(sp.Eq(dln_Rtarget.subs(u_epsEta, 0).subs(u_deltaU, track_surface), 0), u_Lambda)[0]

print("Orbit-lock surface from d ln R_tr = 0:")
sp.pprint(sp.Eq(u_deltaU, track_surface))
print("Target-ratio surface after imposing d ln epsilon_eta = 0 and the tracking surface:")
sp.pprint(sp.Eq(u_Lambda, rtarget_surface))

expect_zero("dlnR_tr on orbit-lock surface", dln_Rtr.subs(u_deltaU, track_surface))
expect_zero(
    "dlnR_target on physical finish surface",
    dln_Rtarget.subs({u_epsEta: 0, u_deltaU: track_surface, u_Lambda: rtarget_surface}),
)

print("\nInterpretation:")
print("  The first three finish-line conditionals are exact codrift surfaces in the actual")
print("  coherent branch variables (chi_0, delta_U, epsilon_W, epsilon_eta, Z_W, Lambda).")
print("  They are not support conditions and do not involve zeta.")
