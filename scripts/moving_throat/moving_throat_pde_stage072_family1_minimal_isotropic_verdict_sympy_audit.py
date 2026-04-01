#!/usr/bin/env python3
"""
moving_throat_pde_stage72_family1_minimal_isotropic_verdict_sympy_audit.py

SymPy / arithmetic audit for Stage 72.

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

banner("STAGE 72 — EXPLICIT FAMILY-1 VERDICT FOR THE MINIMAL ISOTROPIC BRANCH")

rho_min = sp.Rational(4,3)
zeta_min = sp.Rational(1,3)

rho_suff = sp.nsimplify("3.46622291347846")
rho_fail = sp.nsimplify("3.46752913273870")
rho_max  = sp.nsimplify("3.46752922945601")
zeta_max = sp.nsimplify("2.46752922945601")
A_F1     = sp.nsimplify("1.00005192880220")

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

print("\nRegime checks:")
print("  rho_min < rho_suff   -> guaranteed success")
print("  zeta_min < 1         -> symmetric lowest twin already enough")
print("  zeta_min < A_F1      -> Pe_req = 0 on the explicit Family-1 transport map")
