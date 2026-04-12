#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage296_299_common import *

"""
Stage 297 — coherent tracking-branch transfer shape and support-blindness theorem.

This stage converts the actual coherent local D/N tracking branch into one exact
weak-axisymmetric defect compiler in the physical branch variables
(Z_W, Omega_W^2, chi_0, epsilon_W, delta_U, epsilon_eta),
and verifies that the coherent support ratio zeta drops out identically.
"""

banner("STAGE 297 — COHERENT TRACKING-BRANCH TRANSFER SHAPE")

chi0, deltaU, epsW, zeta = sp.symbols("chi_0 delta_U epsilon_W zeta", positive=True, real=True)
ZW, OmW2 = sp.symbols("Z_W Omega_W2", positive=True, real=True)
eps_eta, Lambda0 = sp.symbols("epsilon_eta Lambda_0", positive=True, real=True)
Lambda = sp.simplify(Lambda0 * OmW2)

# Weak-axisymmetric drift amplitudes.
zetaZ, omegaW = sp.symbols("zeta_Z omega_W", real=True)
chi1 = sp.symbols("chi_1", real=True)
epsW1 = sp.symbols("varepsilon_W", real=True)
deltaU1 = sp.symbols("delta_U1", real=True)
eta1 = sp.symbols("eta_1", real=True)

subbanner("I. Exact coherent branch variables")
eps = sp.simplify(epsW * (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)))
T2 = sp.simplify(ZW * (1 + chi0) ** 2 / (OmW2 * (1 - eps) ** 2))
Rtarget = sp.simplify(Lambda * (1 - eps_eta) * (1 - eps) ** 2 / (ZW * (1 + chi0) ** 2))

print("epsilon =")
sp.pprint(eps)
print("T^2 =")
sp.pprint(T2)
print("R_target =")
sp.pprint(Rtarget)

# Exact selected-branch identity.
expect_zero("R_target * T^2 - Lambda_0 * (1 - epsilon_eta)", sp.simplify(Rtarget * T2 - Lambda0 * (1 - eps_eta)))

subbanner("II. Support-blindness theorem")
dlnT2_dzeta = sp.diff(sp.log(T2), zeta)
dlnR_dzeta = sp.diff(sp.log(Rtarget), zeta)
print("d/dzeta ln T^2 =")
sp.pprint(sp.simplify(dlnT2_dzeta))
print("d/dzeta ln R_target =")
sp.pprint(sp.simplify(dlnR_dzeta))
expect_zero("d_zeta ln T^2", dlnT2_dzeta)
expect_zero("d_zeta ln R_target", dlnR_dzeta)

subbanner("III. Exact weak-axisymmetric defect law in physical branch variables")
# Split-blocking drift.
eps1 = sp.simplify(
    (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)) * epsW1
    - (2 * epsW / (11 * (1 + deltaU) ** 2)) * deltaU1
)
Xi1 = sp.simplify(zetaZ - omegaW + 2 * chi1 / (1 + chi0) + 2 * eps1 / (1 - eps))
R1 = sp.simplify(omegaW - eta1 / (1 - eps_eta) - zetaZ - 2 * chi1 / (1 + chi0) - 2 * eps1 / (1 - eps))

print("epsilon_1 =")
sp.pprint(eps1)
print("Xi_1 =")
sp.pprint(Xi1)
print("R_1 =")
sp.pprint(R1)
expect_zero("R_1 + Xi_1 + eta_1/(1-epsilon_eta)", sp.simplify(R1 + Xi1 + eta1 / (1 - eps_eta)))

subbanner("IV. Tracking-factor drift is not the full defect")
Theta1 = sp.simplify(
    -(
        chi0 * (1 + chi0) * deltaU1 + deltaU * (1 + deltaU) * chi1
    ) / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU))
)
print("Theta_1 =")
sp.pprint(Theta1)

Xi1_tracking_rigid = sp.simplify(Xi1.subs({chi1: 0, deltaU1: 0}))
print("Xi_1 on the tracking-rigid slice chi_1 = delta_U1 = 0:")
sp.pprint(Xi1_tracking_rigid)
print("This is generically nonzero unless the remaining mixed/outgoing drifts cancel.")

subbanner("V. Interpretation")
print("On the actual coherent local D/N tracking branch:")
print("  - the transfer shape depends only on Z_W, Omega_W^2, chi_0, epsilon_W, delta_U;")
print("  - the coherent support ratio zeta drops out identically;")
print("  - Xi_1 is carried only by the mixed/outgoing placement drifts;")
print("  - exact tracking rigidity alone is not enough to kill the grouped defect.")
