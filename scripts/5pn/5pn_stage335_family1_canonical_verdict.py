#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage333_335_common import *

"""
Stage 335 — explicit Family-1 verdict for the canonical isotropic branch.

What this script does
---------------------
1. Compares the canonical isotropic loading ratio rho_alpha = 4/3 to the explicit
   Family-1 success/failure/ceiling windows.
2. Verifies the large exact margins to the explicit thresholds.
3. Uses the exact zero-bias support trigger A_F1 to prove Pe_req = 0 on the
   canonical isotropic branch.
4. Records the resulting reduced finish-line statement:
   the explicit support/source side is not the active bottleneck anymore; the
   remaining reduced question is the actual rho_alpha selected by the
   passive/outgoing branch.
"""

banner("STAGE 335 — EXPLICIT FAMILY-1 VERDICT FOR THE CANONICAL ISOTROPIC BRANCH")

rho_min = sp.Rational(4, 3)
zeta_min = sp.Rational(1, 3)

rho_suff = sp.Float("3.46622291347846")
rho_fail = sp.Float("3.46752913273870")
rho_max  = sp.Float("3.46752922945601")
A_F1     = sp.Float("1.00005192880220")
zeta_max_F1 = sp.Float("2.46752922945601")

subbanner("I. Exact margins to the explicit Family-1 windows")
Delta_suff = sp.N(rho_suff - rho_min, 30)
Delta_fail = sp.N(rho_fail - rho_min, 30)
Delta_max  = sp.N(rho_max  - rho_min, 30)

print("rho_alpha^(min) = 4/3")
print("Delta_suff = rho_suff - 4/3 =")
sp.pprint(Delta_suff)
print("Delta_fail = rho_fail - 4/3 =")
sp.pprint(Delta_fail)
print("Delta_max  = rho_max  - 4/3 =")
sp.pprint(Delta_max)

if not (rho_min < rho_suff < rho_fail < rho_max):
    raise AssertionError("Family-1 window ordering failed")
print("Ordering check passed: rho_min < rho_suff < rho_fail < rho_max")

subbanner("II. Exact support-ratio and regime comparison")
zeta_margin = sp.N(zeta_max_F1 - zeta_min, 30)
print("zeta_req^(min) = 1/3")
print("zeta_max^(F1) - zeta_req^(min) =")
sp.pprint(zeta_margin)
if not (0 < zeta_min < 1):
    raise AssertionError("Canonical isotropic branch is not in the symmetric-lowest-twin window")
print("Regime check passed: 0 < zeta_req^(min) < 1, so the canonical isotropic branch")
print("lies strictly inside the symmetric-lowest-twin regime.")

subbanner("III. Zero-transport-bias result on the explicit Family-1 branch")
print("A_F1 =")
sp.pprint(A_F1)
if not (zeta_min < A_F1):
    raise AssertionError("Zero-bias Family-1 support trigger failed")
print("Since zeta_req^(min) = 1/3 < A_F1, the exact transport-bias requirement is")
print("Pe_req = 0 on the explicit Family-1 branch.")

subbanner("IV. Reduced finish-line statement")
print("The explicit Family-1 support/source side is comfortably compatible with the")
print("canonical isotropic passive/outgoing quadrupole precursor.")
print("So the remaining reduced question is not on the support/source side anymore.")
print("It is exactly the outgoing-branch loading ratio rho_alpha selected by the")
print("actual passive/outgoing moving-throat quadrupole branch.")
