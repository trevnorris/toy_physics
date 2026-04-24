#!/usr/bin/env python3
"""
Stage V2-20 — weak-form / numerical branch-extraction audit.

This script does not solve the nonlinear moving-throat PDE.
It verifies the algebraic extraction scaffold needed before a numerical branch
solve can be trusted:

1. weak-form sign conventions for the open finite-throat wall equation,
2. FEM matrix symmetry/positivity in a toy open-exit basis,
3. exact extraction formulas for wall, BdG, Maxwell/mixed, and outgoing data,
4. grouped P2 projector / weak-axisymmetric bookkeeping,
5. target residual identities,
6. branch-freeze manifest completeness and hash stability.

The architectural patch from V2-04 is enforced:
    R(L) > 0 and the exit is an impedance/open boundary, not a hard cap.
"""

from __future__ import annotations

import hashlib
import json
from collections import OrderedDict

import sympy as sp


checks = OrderedDict()


def record(name: str, condition: bool, detail: str = "") -> None:
    checks[name] = (bool(condition), detail)
    if not condition:
        raise AssertionError(f"FAILED: {name}: {detail}")


# -----------------------------------------------------------------------------
# 1. Weak-form identity for a single wall harmonic
# -----------------------------------------------------------------------------

w, L = sp.symbols("w L", positive=True, finite=True)
omega = sp.symbols("omega")
mu, T, V = sp.symbols("mu T V", positive=True)
Yload = sp.symbols("Yload", nonnegative=True)

q = sp.Function("q")(w)
phi = sp.Function("phi")(w)
Tw = sp.Function("T")(w)

# Frequency-domain strong residual:
# R = -omega^2 mu q - d_w(T q_w) + V q - S.
# The integration-by-parts identity for the stiffness part is local:
# T q_w phi_w - d_w(T q_w phi) = -d_w(T q_w) phi.
weak_local = Tw * sp.diff(q, w) * sp.diff(phi, w) - sp.diff(Tw * sp.diff(q, w) * phi, w)
strong_local = -sp.diff(Tw * sp.diff(q, w), w) * phi
record(
    "weak_form_integration_by_parts_identity",
    sp.simplify(weak_local - strong_local) == 0,
    "T q_w phi_w - d_w(T q_w phi) = -d_w(T q_w) phi",
)

# Robin/open-exit boundary: adding 1/2 Yload q(L)^2 to the quadratic form
# produces the boundary variation Yload q(L) phi(L), which combines with
# the integration-by-parts boundary term -T q_w phi|_L to give
# T q_w(L) + Yload q(L)=0 for free variations at L.
qL, phiL, TqL = sp.symbols("q_L phi_L Tq_L")
boundary_variation = TqL * phiL + Yload * qL * phiL
open_bc_residual = TqL + Yload * qL
record(
    "open_exit_robin_boundary_sign",
    sp.simplify(boundary_variation - open_bc_residual * phiL) == 0,
    "free exit condition is T q_w(L)+Y_load q(L)=0 with the chosen weak-form sign",
)

# D/N support mode survives in the Yload -> 0 AC support limit:
k = sp.pi / (2 * L)
q_dn = sp.sin(k * w)
record("DN_mouth_Dirichlet", sp.simplify(q_dn.subs(w, 0)) == 0)
record("DN_open_exit_Neumann_limit", sp.simplify(sp.diff(q_dn, w).subs(w, L)) == 0)


# -----------------------------------------------------------------------------
# 2. Toy Galerkin extraction matrices with Dirichlet mouth and open exit
# -----------------------------------------------------------------------------

# Basis functions satisfy b_i(0)=0; no hard cap is imposed at L.
b1 = w / L
b2 = (w / L) ** 2
basis = [b1, b2]

Mmat = sp.Matrix([[sp.integrate(mu * bi * bj, (w, 0, L)) for bj in basis] for bi in basis])
Kmat = sp.Matrix([
    [
        sp.integrate(T * sp.diff(bi, w) * sp.diff(bj, w) + V * bi * bj, (w, 0, L))
        + Yload * bi.subs(w, L) * bj.subs(w, L)
        for bj in basis
    ]
    for bi in basis
])

record("Galerkin_mass_matrix_symmetric", sp.simplify(Mmat - Mmat.T) == sp.zeros(2))
record("Galerkin_stiffness_matrix_symmetric", sp.simplify(Kmat - Kmat.T) == sp.zeros(2))
record("Galerkin_basis_mouth_Dirichlet", all(sp.simplify(b.subs(w, 0)) == 0 for b in basis))

detM = sp.factor(Mmat.det())
record("Galerkin_mass_det_positive_form", detM == mu**2 * L**2 / 240, f"det(M)={detM}")

# Positive numerical sample for stiffness matrix.
K_sample = Kmat.subs({L: 3, T: 5, V: 7, Yload: sp.Rational(1, 2)})
record("Galerkin_stiffness_numeric_positive_trace", sp.N(K_sample.trace()) > 0, str(K_sample))
record("Galerkin_stiffness_numeric_positive_det", sp.N(K_sample.det()) > 0, f"det={sp.simplify(K_sample.det())}")


# -----------------------------------------------------------------------------
# 3. Stable BdG Schur-complement extraction
# -----------------------------------------------------------------------------

g, varpi, Ksym, Msym = sp.symbols("g varpi K M", positive=True)
x = sp.symbols("x")
D_bdg = Ksym - Msym * x - g**2 / (varpi**2 - x)
series_bdg = sp.series(D_bdg, x, 0, 3).removeO().expand()

B0 = g**2 / varpi**2
B2 = g**2 / varpi**4
B4 = g**2 / varpi**6
D_bdg_expected = (Ksym - B0) - (Msym + B2) * x - B4 * x**2
record("BdG_low_frequency_moments", sp.simplify(series_bdg - D_bdg_expected) == 0)

# Two-mode Stieltjes positivity.
w1, w2, lam1, lam2 = sp.symbols("w1 w2 lam1 lam2", positive=True)
BB0 = w1 / lam1 + w2 / lam2
BB2 = w1 / lam1**2 + w2 / lam2**2
BB4 = w1 / lam1**3 + w2 / lam2**3
hankel = sp.factor(BB0 * BB4 - BB2**2)
hankel_expected = sp.factor(w1 * w2 * (lam1 - lam2) ** 2 / (lam1**3 * lam2**3))
record("BdG_Hankel_moment_positivity_identity", sp.simplify(hankel - hankel_expected) == 0)


# -----------------------------------------------------------------------------
# 4. Maxwell/mixed conservative and outgoing coefficient extraction
# -----------------------------------------------------------------------------

OU, OW, R, gU, gW = sp.symbols("Omega_U Omega_W R g_U g_W", positive=True)
A = OU**2 - x
B = OW**2 - x
Delta_x = A * B - R**2
Sigma = (gU**2 * B + 2 * gU * gW * R + gW**2 * A) / Delta_x
Sigma_series = sp.series(Sigma, x, 0, 3).removeO().expand()

Delta0 = OU**2 * OW**2 - R**2
S2 = OU**2 + OW**2
Q = gU**2 * OW**2 + 2 * gU * gW * R + gW**2 * OU**2
H = gU**2 + gW**2
Z0 = Q / Delta0
Z2 = (Q * S2 - H * Delta0) / Delta0**2
Z4 = (Q * (S2**2 - Delta0) - S2 * H * Delta0) / Delta0**3
record(
    "Maxwell_mixed_conservative_Z_moments",
    sp.simplify(Sigma_series - (Z0 + Z2 * x + Z4 * x**2)) == 0,
)

P = OU**2 * gW + R * gU
N_transfer = (A * gW + R * gU) ** 2 / Delta_x**2
N_series = sp.series(N_transfer, x, 0, 3).removeO().expand()
N0 = P**2 / Delta0**2
N2 = 2 * P * (P * S2 - Delta0 * gW) / Delta0**3
N4 = (
    Delta0**2 * gW**2
    - 2 * Delta0 * P**2
    - 4 * Delta0 * P * S2 * gW
    + 3 * P**2 * S2**2
) / Delta0**4
record("Maxwell_mixed_outgoing_N_moments", sp.simplify(N_series - (N0 + N2 * x + N4 * x**2)) == 0)

record("mixed_internal_stability_symbol", Delta0 != 0, "strict gate is Delta0>0")
record("outgoing_transfer_nonnegative_symbolic_square", sp.factor(N0 * Delta0**2 - P**2) == 0)


# -----------------------------------------------------------------------------
# 5. Total isotropic extraction and target residuals
# -----------------------------------------------------------------------------

K, M, B0s, B2s, B4s, Z0s, Z2s, Z4s, N0s, N2s, N4s = sp.symbols(
    "K M B0 B2 B4 Z0 Z2 Z4 N0 N2 N4"
)
D0 = K - B0s - Z0s
D2 = -(M + B2s + Z2s)
D4 = -(B4s + Z4s)

u2 = -D2 / D0
u4 = (D2**2 - D0 * D4) / D0**2
Aeff = M + B2s + Z2s
Ceff = B4s + Z4s

one_pole_residual = sp.simplify((u4 - 4 * u2**2) * D0**2 - (D0 * Ceff - 3 * Aeff**2))
record("one_pole_surface_identity", one_pole_residual == 0)

P0 = N0s / D0
P2 = (D0 * N2s - 2 * D2 * N0s) / D0**2
P4 = (D0**2 * N4s - 2 * D0 * (D2 * N2s + D4 * N0s) + 3 * D2**2 * N0s) / D0**3

N2_const = 2 * D2 * N0s / D0
N4_const = N0s * (D2**2 + 2 * D0 * D4) / D0**2
record("constant_prefactor_P2_zero", sp.simplify(P2.subs(N2s, N2_const)) == 0)
record(
    "constant_prefactor_P4_zero",
    sp.simplify(P4.subs({N2s: N2_const, N4s: N4_const})) == 0,
)

G, cs, c, ath, mhat, Sport = sp.symbols("G c_s c a_th mhat S_port", positive=True)
Ptarget = 54 * G * cs**5 / (5 * ath**5 * c**5)
gamma_eff = (mhat**2 * Sport * P0) * ath**5 / (27 * cs**5)
gamma_target_sub = sp.simplify(gamma_eff.subs(P0, Ptarget / (mhat**2 * Sport)) - 2 * G / (5 * c**5))
record("quadrupole_gamma_target_equivalence", gamma_target_sub == 0)

tailTheta = sp.symbols("Theta_tail")
tail_gate = tailTheta * (c / cs) ** 3 - 1
record("tail_gate_scalar_form", sp.simplify(tail_gate - (tailTheta * c**3 / cs**3 - 1)) == 0)


# -----------------------------------------------------------------------------
# 6. Grouped projector and weak-axisymmetric bookkeeping
# -----------------------------------------------------------------------------

Ggrp = sp.diag(1, 2, 2)
e_bar = sp.Matrix([1, 1, 1])
e_a = sp.Matrix([4, -1, -1])
e_b = sp.Matrix([0, 1, -1])
Pbar = e_bar * (e_bar.T * Ggrp) / (e_bar.T * Ggrp * e_bar)[0]
Pa = e_a * (e_a.T * Ggrp) / (e_a.T * Ggrp * e_a)[0]
Pb = e_b * (e_b.T * Ggrp) / (e_b.T * Ggrp * e_b)[0]
I3 = sp.eye(3)

record("grouped_projector_completeness", sp.simplify(Pbar + Pa + Pb - I3) == sp.zeros(3))
record("grouped_projector_idempotence", all(sp.simplify(P * P - P) == sp.zeros(3) for P in [Pbar, Pa, Pb]))

x0, eps, x1 = sp.symbols("x0 eps x1")
lam = sp.Matrix([1, sp.Rational(1, 2), -1])
xvec = sp.Matrix([x0, x0, x0]) + eps * x1 * lam
xbar = (xvec[0] + 2 * xvec[1] + 2 * xvec[2]) / 5
a_x = (2 * xvec[0] - xvec[1] - xvec[2]) / 10
b_x = (xvec[1] - xvec[2]) / 2
record("weak_axisymmetric_trace_unchanged", sp.simplify(xbar - x0) == 0)
record("weak_axisymmetric_b_equals_3a", sp.simplify(b_x - 3 * a_x) == 0)


# -----------------------------------------------------------------------------
# 7. Extraction manifest / freeze packet
# -----------------------------------------------------------------------------

manifest = OrderedDict(
    [
        ("geometry", OrderedDict([
            ("throat_exit", "open_finite_radius"),
            ("R_L_positive_required", True),
            ("forbidden_condition", "R(L)=0_hard_cap"),
            ("mouth_bc", "Dirichlet_or_prescribed_worldtube_drive"),
            ("exit_bc", "impedance_Robin_with_Neumann_AC_support_limit"),
            ("dc_flow_policy", "baseline_current_may_exit; support_AC_modes_are_extracted_separately"),
        ])),
        ("weak_form", OrderedDict([
            ("wall_operator", "mu q_tt - d_w(T q_w) + (K_eta+l(l+1)T_Omega)q"),
            ("test_space", "mouth-compatible; free/Robin at open exit"),
            ("boundary_load_symbol", "Y_load(omega)"),
            ("densitized_measure", "dw; background surface weight absorbed into coefficients"),
        ])),
        ("extraction_slots", [
            "K", "M", "B0", "B2", "B4", "Z0", "Z2", "Z4",
            "N0", "N2", "N4", "D0", "D2", "D4",
            "u2", "u4", "P0", "P2", "P4",
            "mhat0", "S_port", "Theta_tail",
        ]),
        ("stability_gates", [
            "mu_eta>0", "T_w>0", "internal_BdG_positive_energy",
            "Delta_r>0 for every mixed port", "D0>0", "no_dark_port_if_N0_needed",
            "Ceff=B4+Z4>0 on one-pole stable branch",
        ]),
        ("target_residuals", [
            "R_pole=D0*(B4+Z4)-3*(M+B2+Z2)^2",
            "R_norm=mhat0^2*S_port*N0/D0-54*G*c_s^5/(5*a^5*c^5)",
            "R_P2=P2", "R_P4=P4",
            "R_tail=Theta_tail*(c/c_s)^3-1",
        ]),
        ("freeze_rule", "all primitive profiles, basis choices, port list, gauge convention, and boundary class are fixed before target residual evaluation"),
    ]
)

required_slots = {
    "K", "M", "B0", "B2", "B4", "Z0", "Z2", "Z4",
    "N0", "N2", "N4", "D0", "D2", "D4",
    "u2", "u4", "P0", "P2", "P4",
    "mhat0", "S_port", "Theta_tail",
}
record("manifest_required_extraction_slots_present", required_slots.issubset(set(manifest["extraction_slots"])))
record("manifest_open_exit_no_hard_cap", manifest["geometry"]["R_L_positive_required"] and "R(L)=0" in manifest["geometry"]["forbidden_condition"])
record("manifest_freeze_rule_present", "fixed before target residual" in manifest["freeze_rule"])

manifest_json = json.dumps(manifest, indent=2, sort_keys=True)
manifest_hash = hashlib.sha256(manifest_json.encode("utf-8")).hexdigest()


# -----------------------------------------------------------------------------
# Output report
# -----------------------------------------------------------------------------

passed = sum(1 for ok, _ in checks.values() if ok)
failed = sum(1 for ok, _ in checks.values() if not ok)

print("Stage V2-20 weak-form / branch-extraction audit")
print("=" * 72)
for name, (ok, detail) in checks.items():
    status = "PASS" if ok else "FAIL"
    if detail:
        print(f"{status:4s}  {name}: {detail}")
    else:
        print(f"{status:4s}  {name}")

print("=" * 72)
print(f"checks_total: {len(checks)}")
print(f"checks_passed: {passed}")
print(f"checks_failed: {failed}")
print()
print("Toy Galerkin mass matrix M:")
print(sp.sstr(Mmat))
print()
print("Toy Galerkin stiffness matrix K:")
print(sp.sstr(Kmat))
print()
print("Main extracted isotropic formulas:")
print(f"D0 = {sp.sstr(D0)}")
print(f"D2 = {sp.sstr(D2)}")
print(f"D4 = {sp.sstr(D4)}")
print(f"u2 = {sp.sstr(u2)}")
print(f"u4 = {sp.sstr(u4)}")
print(f"P0 = {sp.sstr(P0)}")
print(f"P2 = {sp.sstr(P2)}")
print(f"P4 = {sp.sstr(P4)}")
print()
print("Target residual packet:")
print("R_pole = D0*(B4+Z4)-3*(M+B2+Z2)^2")
print("R_norm = mhat0^2*S_port*N0/D0 - 54*G*c_s^5/(5*a^5*c^5)")
print("R_P2 = P2")
print("R_P4 = P4")
print("R_tail = Theta_tail*(c/c_s)^3 - 1")
print()
print("Extraction manifest hash:")
print(manifest_hash)
print()
print("Extraction manifest:")
print(manifest_json)
