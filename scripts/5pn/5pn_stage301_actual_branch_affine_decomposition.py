#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage300_302_common import *

"""
Stage 301 — exact affine decomposition of an actual branch tangent.

This stage solves the Stage-299 monomial-drift equations for one convenient
split of the microscopic drifts, exhibits an explicit five-dimensional orbit basis
plus a three-dimensional normal basis, and proves that any actual microscopic drift
vector decomposes uniquely into
    tangent-to-similarity-orbit part + off-orbit defect coordinates.
"""

banner("STAGE 301 — ACTUAL BRANCH AFFINE DECOMPOSITION")

chi0s, deltaUs, epsWs, epss = sp.symbols("chi0_star deltaU_star epsilonW_star epsilon_star", positive=True, real=True)
Mstar, E_star, F_star = monomial_drift_matrix(chi0s, deltaUs, epsWs, epss)

# Drift coordinates.
dl_lam, dl_c, dl_gam, dl_KU, dl_Keta, dl_KW, dl_muW, dl_TU = sp.symbols(
    "dl_lambda dl_c dl_gamma dl_KU dl_Keta dl_KW dl_muW dl_TU", real=True
)
qtr, qnt, qeta = sp.symbols("q_tr q_nt q_eta", real=True)

subbanner("I. Exact affine solve for the dependent drifts")
dlKeta = sp.simplify(2 * dl_c - dl_KU - qeta)
dlTU = sp.simplify(
    dl_KU + (qtr - (1 + deltaUs) * (dl_gam + dl_c - dl_KU)) / (1 + chi0s)
)
dlMuW = sp.simplify(
    qnt - qeta - dl_KU + 2 * dl_KW + 2 * dl_c - 2 * dl_lam
    + E_star * (dl_KU + dl_KW - 2 * dl_gam - 2 * dl_lam)
    + F_star * ((1 + deltaUs) * (dl_KU - dl_c - dl_gam) + qtr) / (1 + chi0s)
)

print("dl_Keta =")
sp.pprint(dlKeta)
print("dl_TU =")
sp.pprint(dlTU)
print("dl_muW =")
sp.pprint(dlMuW)

dx_aff = sp.Matrix([dl_lam, dl_c, dl_gam, dl_KU, dlKeta, dl_KW, dlMuW, dlTU])
print("Affine branch-drift solve dx_aff =")
sp.pprint(dx_aff)
expect_zero("M_* dx_aff - q", sp.simplify(Mstar * dx_aff - sp.Matrix([qtr, qnt, qeta])))

subbanner("II. Convenient exact orbit basis (five free tangent directions)")
free_syms = [dl_lam, dl_c, dl_gam, dl_KU, dl_KW]
orbit_basis = {}
for sym in free_syms:
    subs = {
        dl_lam: 0, dl_c: 0, dl_gam: 0, dl_KU: 0, dl_KW: 0,
        qtr: 0, qnt: 0, qeta: 0,
    }
    subs[sym] = 1
    vec = sp.Matrix([
        subs[dl_lam],
        subs[dl_c],
        subs[dl_gam],
        subs[dl_KU],
        sp.simplify(dlKeta.subs(subs)),
        subs[dl_KW],
        sp.simplify(dlMuW.subs(subs)),
        sp.simplify(dlTU.subs(subs)),
    ])
    orbit_basis[str(sym)] = vec
    print(f"orbit basis vector for {sym} =")
    sp.pprint(vec)
    expect_zero(f"M_* orbit[{sym}]", sp.simplify(Mstar * vec))

subbanner("III. Exact normal basis (three off-orbit defect directions)")
normal_basis = {}
for name, vals in {
    "n_tr": {qtr: 1, qnt: 0, qeta: 0},
    "n_nt": {qtr: 0, qnt: 1, qeta: 0},
    "n_eta": {qtr: 0, qnt: 0, qeta: 1},
}.items():
    subs = {
        dl_lam: 0, dl_c: 0, dl_gam: 0, dl_KU: 0, dl_KW: 0,
        **vals,
    }
    vec = sp.Matrix([
        0,
        0,
        0,
        0,
        sp.simplify(dlKeta.subs(subs)),
        0,
        sp.simplify(dlMuW.subs(subs)),
        sp.simplify(dlTU.subs(subs)),
    ])
    normal_basis[name] = vec
    print(f"{name} =")
    sp.pprint(vec)

Nmat = sp.Matrix.hstack(normal_basis["n_tr"], normal_basis["n_nt"], normal_basis["n_eta"])
print("Normal-basis matrix N =")
sp.pprint(Nmat)
expect_zero("M_* N - I_3", sp.simplify(Mstar * Nmat - sp.eye(3)))

subbanner("IV. Exact affine decomposition theorem")
qvec = sp.Matrix([qtr, qnt, qeta])
orbit_vec = dx_aff.subs({qtr: 0, qnt: 0, qeta: 0})
reconstructed = sp.simplify(orbit_vec + Nmat * qvec)
print("orbit part dx_orbit =")
sp.pprint(orbit_vec)
print("reconstructed dx = dx_orbit + N q =")
sp.pprint(reconstructed)
expect_zero("dx_aff - (dx_orbit + N q)", sp.simplify(dx_aff - reconstructed))

subbanner("V. Interpretation")
print("Any actual microscopic grouped weak-axisymmetric drift vector now decomposes uniquely as:")
print("  dx_actual = dx_orbit + q_tr n_tr + q_nt n_nt + q_eta n_eta.")
print("Here dx_orbit is tangent to the exact similarity orbit, while")
print("(q_tr,q_nt,q_eta) are the three off-orbit defect coordinates.")
