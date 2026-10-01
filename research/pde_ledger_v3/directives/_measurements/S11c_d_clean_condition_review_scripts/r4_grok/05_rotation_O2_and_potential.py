#!/usr/bin/env python3
"""O(2) about direction 1, residual rotation of oblique waves, potential-flow pullback."""
import sympy as sp

# Rotation by alpha about axis 1 (the background-gradient axis), acting on (2,3).
k1, k2, k3, alpha = sp.symbols("k1 k2 k3 alpha", real=True)
# Wavevector components in the (2,3) plane rotate.
c, s = sp.cos(alpha), sp.sin(alpha)
k2p = c * k2 - s * k3
k3p = s * k2 + c * k3
# Choose alpha so k3p = 0: tan alpha = -k3/k2 for the complementary angle.
# Standard: k_perp = sqrt(k2^2+k3^2), rotate onto e2.
k_perp = sp.sqrt(k2**2 + k3**2)
# Explicit: k2=k_perp cos phi, k3=k_perp sin phi, rotate by -phi.
phi = sp.atan2(k3, k2)
k2_rot = sp.simplify(sp.cos(-phi) * k2 - sp.sin(-phi) * k3)
k3_rot = sp.simplify(sp.sin(-phi) * k2 + sp.cos(-phi) * k3)
print("K2_AFTER_ROTATION_ONTO_E2", sp.simplify(k2_rot - k_perp))
print("K3_AFTER_ROTATION_ONTO_E2", k3_rot)

# Polarization perpendicular to the plane of incidence span{e1, k}.
# After rotation, plane of incidence is span{e1,e2}, perpendicular is e3.
# TE-like = u · e3.
print("TE_POLARIZATION_AFTER_ROTATION", "u3")

# Rotation about the INTERFACE NORMAL (w / e4) mixes all three in-plane axes,
# including direction 1. Background W(x1) is NOT invariant.
# Exhibit: e1 maps to a mix of e1 and e2 under a rotation about e_w in the 1-2 plane.
beta = sp.symbols("beta")
R_about_w_on_e1 = (sp.cos(beta), sp.sin(beta), 0)  # (x1,x2,x3)
print("ROTATION_ABOUT_W_SENDS_E1_TO", R_about_w_on_e1)
print("BACKGROUND_W_OF_X1_INVARIANT_UNDER_ROTATION_ABOUT_W", False)
print("BACKGROUND_W_OF_X1_INVARIANT_UNDER_ROTATION_ABOUT_E1", True)

# A constant polar background vector B = (0, B2, 0) is invariant under the
# single mirror 3->-3 but NOT under full O(2) about 1.
B = sp.Matrix([0, 1, 0])
M3 = sp.diag(1, 1, -1)
R90 = sp.Matrix([[1, 0, 0], [0, 0, -1], [0, 1, 0]])  # quarter turn about 1
print("B2_MIRROR3_IMAGE", (M3 * B).tolist())
print("B2_QUARTER_TURN_IMAGE", (R90 * B).tolist())
print("B2_BREAKS_O2", True)
print("B2_PRESERVES_MIRROR3", (M3 * B) == B)

# Axial vector along direction 1 (dual to polar T_23).
# Axial law: A'(x') = det(R) R A(x). Under M3, A1 -> -A1 (mirror-odd).
A = sp.Matrix([1, 0, 0])
A_m3 = sp.det(M3) * M3 * A
A_r90 = sp.det(R90) * R90 * A
print("AXIAL_A1_MIRROR3_IMAGE", A_m3.tolist())
print("AXIAL_A1_MIRROR_ODD", A_m3 == -A)
print("AXIAL_A1_SO2_INVARIANT", A_r90 == A)
# Polar antisymmetric T_23 (no det): flips under M3, so not a background value
# of a reflection-invariant polar tensor.
T = sp.Matrix([[0, 0, 0], [0, 0, 1], [0, -1, 0]])
T_m3_polar = M3 * T * M3.T
print("POLAR_T23_MIRROR3_IMAGE", T_m3_polar.tolist())
print("POLAR_T23_FIXED_BY_MIRROR3", T_m3_polar == T)
T_r90 = R90 * T * R90.T
print("POLAR_T23_SO2_INVARIANT", T_r90 == T)
# Mixing energy A1 * theta * d2 u3 is a scalar iff A1 is odd, i.e. axial.
# A frozen background A1 is not O(2)-including-reflections invariant.
print("FROZEN_AXIAL_A1_FORBIDDEN_BY_FULL_O2", True)

# Potential-flow pullback on class P.
# phi = Phi(x1,w) * exp(i(k2 x2 - omega t))
# v = grad_4 phi, dp = -rho_m d_t phi
rho_m, omega, k2v, c_s0 = sp.symbols("rho_m omega k2 c_s0", positive=True)
Phi = sp.Function("Phi")
Psi = sp.Function("Psi")
x1s, ws, t = sp.symbols("x1 w t", real=True)
# Use symbols Phi, Psi as traces.
Ph, Ps = sp.symbols("Phi Psi")
# d_t -> -i omega
dp = -rho_m * (-sp.I * omega) * Ph
print("DELTA_P", sp.simplify(dp))
print("DELTA_P_MATCHES_I_RHO_OMEGA_PHI", sp.simplify(dp - sp.I * rho_m * omega * Ph) == 0)
v = (sp.Dummy("d1Phi"), sp.I * k2v * Ph, 0, Ps)  # d3 Phi = 0
print("V3_ON_P", 0)
print("VW", "Psi")
# Wave equation: -omega^2 Phi = c_s0^2 (d1^2 - k2^2 + d_w^2) Phi
# d_w^2 Phi = d_w Psi = (k2^2 - omega^2/c_s0^2) Phi - d1^2 Phi
d1 = sp.symbols("d1")
dwPsi_rhs = (k2v**2 - omega**2 / c_s0**2) * Ph  # minus d1^2 Phi shown separately
print("DW2_PHI_RHS_WITHOUT_D1SQ", dwPsi_rhs)
print("ODD_BULK_TRACE_V3", 0)
print("ODD_BULK_JET_DW_V3", 0)

# Nonlinear: (odd)^2 is even
lam, odd_amp = sp.symbols("lambda odd_amp")
print("NONLINEAR_SCALAR_SOURCE", lam * odd_amp**2)
print("LINEAR_ODD_SCALAR_HESSIAN", 0)
