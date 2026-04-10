#!/usr/bin/env python3
"""
5pn_stage67_reduced_pde_writeup_skeleton.py

Stage 67 audit: full reduced moving-throat PDE write-up skeleton and final remaining theorem gap.
"""

from __future__ import annotations

import math
import mpmath as mp
import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.together(sp.expand(expr)))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 67 — REDUCED MOVING-THROAT PDE WRITE-UP SKELETON")

subbanner("67.1 — Hierarchy now in hand")
print("Exact parent bulk system:")
print("  i hbar D_t psi = [ -(hbar^2/2m) D_A D_A + V_conf + h(|psi|^2) ] psi")
print("  partial_M(Z F^(MN)) + (1/xi) partial^N(partial·A) = mu_0 J^N")
print()
print("Reduced support/source throat operator:")
print("  partial_t sigma + partial_s J = 0")
print("  J = -D_sigma partial_s sigma + v_sigma sigma")
print("  -T_X partial_s^2 phi + K_X phi = Lambda_phi sigma")
print("  with Pe_* = Xi Delta(Pe_*;kappa,eta)")
print()
print("Explicit lowest-lane support ratio:")
print("  zeta_phys(Pe,eta;kappa) = Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2)")
print()
print("Selected-branch quadrupole demand:")
print("  zeta_req(Pi_tr,C_mix,eps_blk) = (Pi_tr - C_mix) / [C_mix - eps_blk (2 C_mix - Pi_tr)]")

subbanner("67.2 — Full reduced theorem gate")
print("The full reduced moving-throat PDE is governed by")
print("  R_quad = zeta_req - zeta_phys(Pe_*).")
print("This is the cleanest current write-up form of the remaining theorem gate.")

subbanner("67.3 — Family-1 specialization")
print("On the explicit Family-1 branch:")
print("  kappa_F1 = 12321/5,  eta_F1 = 37,  Xi_F1 = 1369 Upsilon_w = 136900 Theta_w.")
print("The support/source side is therefore explicit; what remains open is the outgoing")
print("quadrupole-normalization branch itself.")

banner("STAGE 67 FINAL LEDGER")
print("Stage 67 is not a new isolated calculation. It is the write-up skeleton showing that")
print("the reduced moving-throat PDE program is fully organized, with one remaining theorem gap:")
print("the physical passive/outgoing quadrupole-normalization branch.")
