
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage199_201_common import *

"""
Stage 199 — exact final branch residual packet from grouped-lane bundle data.

What this script does
---------------------
1. Takes the actual grouped-lane low-frequency bundle data
      (D_A0,D_A2,D_A4,N_A0,N_A2,N_A4) for A in {20,21,22}
   together with the source-map factor mhat_0.
2. Compiles them into the normalized grouped response moments
      u_2^(A), u_4^(A), P_0^(A), P_2^(A), P_4^(A).
3. Extracts the exact grouped isotropy defects for both the conservative response
   and the outgoing prefactor.
4. Isolates the two scalar endgame defects:
      Delta_pole = ubar_4 - 4 ubar_2^2
      Delta_norm = mhat_0^2 Pbar_0 - 54 G c_s^5 / (5 a^5 c^5).
5. Packages the minimal exact final branch residual packet that the completed PDE
   must drive to zero on the reduced GR-compatible branch.
"""

banner("STAGE 199 — EXACT FINAL BRANCH RESIDUAL PACKET")

G, c, cs, a, mhat0 = sp.symbols("G c c_s a mhat_0", positive=True, real=True)

D20_0, D21_0, D22_0 = sp.symbols("D20_0 D21_0 D22_0", positive=True, real=True)
D20_2, D21_2, D22_2 = sp.symbols("D20_2 D21_2 D22_2", real=True)
D20_4, D21_4, D22_4 = sp.symbols("D20_4 D21_4 D22_4", real=True)

N20_0, N21_0, N22_0 = sp.symbols("N20_0 N21_0 N22_0", positive=True, real=True)
N20_2, N21_2, N22_2 = sp.symbols("N20_2 N21_2 N22_2", real=True)
N20_4, N21_4, N22_4 = sp.symbols("N20_4 N21_4 N22_4", real=True)

subbanner("I. Exact grouped-lane normalized response moments")
mom20 = response_moments_from_D(D20_0, D20_2, D20_4)
mom21 = response_moments_from_D(D21_0, D21_2, D21_4)
mom22 = response_moments_from_D(D22_0, D22_2, D22_4)

print("u_2^(20), u_4^(20) =")
sp.pprint(sp.Matrix([mom20["u2"], mom20["u4"]]))
print("u_2^(21), u_4^(21) =")
sp.pprint(sp.Matrix([mom21["u2"], mom21["u4"]]))
print("u_2^(22), u_4^(22) =")
sp.pprint(sp.Matrix([mom22["u2"], mom22["u4"]]))

u2grp = grouped_trace_anomaly(mom20["u2"], mom21["u2"], mom22["u2"])
u4grp = grouped_trace_anomaly(mom20["u4"], mom21["u4"], mom22["u4"])

print("Grouped (ubar_2,a_2,b_2) =")
sp.pprint(sp.Matrix([u2grp["bar"], u2grp["a"], u2grp["b"]]))
print("Grouped (ubar_4,a_4,b_4) =")
sp.pprint(sp.Matrix([u4grp["bar"], u4grp["a"], u4grp["b"]]))

subbanner("II. Exact grouped-lane outgoing prefactor moments")
P20 = prefactor_moments(D20_0, D20_2, D20_4, N20_0, N20_2, N20_4)
P21 = prefactor_moments(D21_0, D21_2, D21_4, N21_0, N21_2, N21_4)
P22 = prefactor_moments(D22_0, D22_2, D22_4, N22_0, N22_2, N22_4)

print("P_0^(20), P_2^(20), P_4^(20) =")
sp.pprint(sp.Matrix([P20["P0"], P20["P2"], P20["P4"]]))
print("P_0^(21), P_2^(21), P_4^(21) =")
sp.pprint(sp.Matrix([P21["P0"], P21["P2"], P21["P4"]]))
print("P_0^(22), P_2^(22), P_4^(22) =")
sp.pprint(sp.Matrix([P22["P0"], P22["P2"], P22["P4"]]))

P0grp = grouped_trace_anomaly(P20["P0"], P21["P0"], P22["P0"])
P2grp = grouped_trace_anomaly(P20["P2"], P21["P2"], P22["P2"])
P4grp = grouped_trace_anomaly(P20["P4"], P21["P4"], P22["P4"])

print("Grouped prefactor P_0 packet =")
sp.pprint(sp.Matrix([P0grp["bar"], P0grp["a"], P0grp["b"]]))
print("Grouped prefactor P_2 packet =")
sp.pprint(sp.Matrix([P2grp["bar"], P2grp["a"], P2grp["b"]]))
print("Grouped prefactor P_4 packet =")
sp.pprint(sp.Matrix([P4grp["bar"], P4grp["a"], P4grp["b"]]))

subbanner("III. Exact isotropic-branch scalar defects")
D0, D2, D4, N0 = sp.symbols("D_0 D_2 D_4 N_0", positive=True, real=True)
Delta_pole_iso = pole_defect_from_D(D0, D2, D4)
P0_iso = sp.simplify(N0 / D0)
Delta_norm_iso = sp.simplify(mhat0**2 * P0_iso - 54*G*cs**5/(5*a**5*c**5))

print("Delta_pole^(iso) =")
sp.pprint(Delta_pole_iso)
print("Delta_norm^(iso) =")
sp.pprint(Delta_norm_iso)

# Check exact isotropic specialization of grouped defects.
subs_iso = {
    D20_0: D0, D21_0: D0, D22_0: D0,
    D20_2: D2, D21_2: D2, D22_2: D2,
    D20_4: D4, D21_4: D4, D22_4: D4,
    N20_0: N0, N21_0: N0, N22_0: N0,
    N20_2: sp.Symbol("N_2"), N21_2: sp.Symbol("N_2"), N22_2: sp.Symbol("N_2"),
    N20_4: sp.Symbol("N_4"), N21_4: sp.Symbol("N_4"), N22_4: sp.Symbol("N_4"),
}

expect_zero("a_2 on isotropic branch", u2grp["a"].subs(subs_iso))
expect_zero("b_2 on isotropic branch", u2grp["b"].subs(subs_iso))
expect_zero("a_4 on isotropic branch", u4grp["a"].subs(subs_iso))
expect_zero("b_4 on isotropic branch", u4grp["b"].subs(subs_iso))
expect_zero("a(P_0) on isotropic branch", P0grp["a"].subs(subs_iso))
expect_zero("b(P_0) on isotropic branch", P0grp["b"].subs(subs_iso))

expect_zero(
    "ubar_4 - 4 ubar_2^2 under isotropic specialization",
    sp.simplify((u4grp["bar"] - 4*u2grp["bar"]**2).subs(subs_iso) - Delta_pole_iso),
)

subbanner("IV. First-order transport of prefactor anisotropy")
dD20_0, dD21_0, dD22_0 = sp.symbols("dD20_0 dD21_0 dD22_0", real=True)
dN20_0, dN21_0, dN22_0 = sp.symbols("dN20_0 dN21_0 dN22_0", real=True)
eps = sp.symbols("eps", real=True)

P0_lin20 = sp.simplify(prefactor_moments(D0 + eps*dD20_0, D2, D4, N0 + eps*dN20_0, 0, 0)["P0"])
P0_lin21 = sp.simplify(prefactor_moments(D0 + eps*dD21_0, D2, D4, N0 + eps*dN21_0, 0, 0)["P0"])
P0_lin22 = sp.simplify(prefactor_moments(D0 + eps*dD22_0, D2, D4, N0 + eps*dN22_0, 0, 0)["P0"])

dP20 = sp.expand(sp.series(P0_lin20, eps, 0, 2).removeO().coeff(eps, 1))
dP21 = sp.expand(sp.series(P0_lin21, eps, 0, 2).removeO().coeff(eps, 1))
dP22 = sp.expand(sp.series(P0_lin22, eps, 0, 2).removeO().coeff(eps, 1))

dDgrp = grouped_trace_anomaly(dD20_0, dD21_0, dD22_0)
dNgrp = grouped_trace_anomaly(dN20_0, dN21_0, dN22_0)
dPgrp = grouped_trace_anomaly(dP20, dP21, dP22)
P0base = sp.simplify(N0 / D0)

expect_zero("linear a(P_0) transport",
            sp.simplify(dPgrp["a"] - (dNgrp["a"] - P0base*dDgrp["a"]) / D0))
expect_zero("linear b(P_0) transport",
            sp.simplify(dPgrp["b"] - (dNgrp["b"] - P0base*dDgrp["b"]) / D0))

subbanner("V. Exact final branch residual packet")
Delta_branch = sp.Matrix([
    u2grp["a"], u2grp["b"],
    u4grp["a"], u4grp["b"],
    P0grp["a"], P0grp["b"],
    sp.simplify(u4grp["bar"] - 4*u2grp["bar"]**2),
    sp.simplify(mhat0**2 * P0grp["bar"] - 54*G*cs**5/(5*a**5*c**5)),
])

print("Delta_branch =")
sp.pprint(Delta_branch)

banner("STAGE 199 LEDGER")
print("1. The full grouped-lane bundle is now compiled into one exact residual packet:")
print("      Delta_branch = (a_2,b_2,a_4,b_4,a(P_0),b(P_0),Delta_pole,Delta_norm).")
print("2. The isotropic grouped branch is exactly the zero set of the first six components.")
print("3. The one-pole conservative test is Delta_pole = 0.")
print("4. The remaining quadrupole-normalization test is Delta_norm = 0.")
print("5. So the completed PDE no longer has to answer a vague question — it must drive")
print("   a specific finite-dimensional residual packet to zero on the physical branch.")
