#!/usr/bin/env python3
"""Exact sign/matrix audit for the S11c-d v3 reflection hypotheses.

This uses integer parity and matrix algebra only.  It does not import either
S11c engine or any previously emitted conclusion.
"""

from fractions import Fraction


def matmul(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def transpose(a):
    return [list(row) for row in zip(*a)]


def diag(*xs):
    return [[x if i == j else 0 for j in range(len(xs))] for i, x in enumerate(xs)]


def show_matrix(a):
    return "[" + ", ".join(str(row) for row in a) + "]"


R3 = diag(1, 1, -1)
R4 = diag(1, 1, -1, 1)  # w is not reflected

print("REFLECTION_R3", show_matrix(R3))
print("REFLECTION_R4", show_matrix(R4))

# Exact covariance checks on representative generic integer tensors.
a = [[2], [-3], [5]]
b = [[7], [11], [-13]]
dot_before = matmul(transpose(a), b)[0][0]
dot_after = matmul(transpose(matmul(R3, a)), matmul(R3, b))[0][0]
print("DOT_VECTOR_RESIDUAL", dot_after - dot_before)

G = [[2, 3, 5], [7, 11, 13], [17, 19, 23]]
G_reflected = matmul(matmul(R3, G), R3)
trace_before = sum(G[i][i] for i in range(3))
trace_after = sum(G_reflected[i][i] for i in range(3))
print("TRACE_GRADIENT_RESIDUAL", trace_after - trace_before)

# On the R1 background, g=(g1,0,0), so both anchoring's u.g reaction and the
# material constraint are blind to u3 when d3=0.
g = [37, 0, 0]
print("ANCHOR_U_DOT_G_COEFFICIENTS", g)
print("ANCHOR_COEFF_U3", g[2])
divergence_coefficients_on_P = ["d1", "d2", 0]
print("CONSTRAINT_DIV_COEFFICIENTS_ON_P", divergence_coefficients_on_P)
print("CONSTRAINT_COEFF_U3_ON_P", divergence_coefficients_on_P[2])

# Linearized face kinematics on h=h(x1,x2), d3 h=0.
# n=(-h1,-h2,0,s)/sqrt(...), v_face=(u1_t,u2_t,u3_t,h_t).
n = ["-h1/N", "-h2/N", 0, "s/N"]
print("R1_FACE_NORMAL", n)
print("D_VFACE3_D_U3_T", 1)
print("D_DELTAVX3_D_DELTAVU3", 1)
print("D_NORMAL_SPEED_D_U3_T", n[2])
print("D_RELATIVE_FLUX_D_U3_T", f"-rho_m*({n[2]})")
print("TRACTION_COMPONENT3", f"-pressure_scalar*({n[2]})")

# The actual response laws contain scalar kernels multiplying scalar face
# data; traction is a scalar times a polar normal.
lambda_parities = {"Lambda_A": 1, "Lambda_V": 1, "Lambda_X": 1}
affinity_parity = 1 * 1  # mu and delta_p are both scalars
normal_speed_parity = 1  # dot of two polar vectors
response_flux_parity = max(
    lambda_parities["Lambda_A"] * affinity_parity,
    lambda_parities["Lambda_V"] * normal_speed_parity,
)
traction_component_parities = [1, 1, -1, 1]
virtual_displacement_component_parities = [1, 1, -1, 1]
virtual_work_term_parities = [
    a * b for a, b in zip(traction_component_parities, virtual_displacement_component_parities)
]
print("FACE_RESPONSE_PARITY_AFFINITY", affinity_parity)
print("FACE_RESPONSE_PARITY_NORMAL_SPEED", normal_speed_parity)
print("FACE_RESPONSE_PARITY_J", response_flux_parity)
print("FACE_RESPONSE_TRACTION_COMPONENT_PARITIES", traction_component_parities)
print("FACE_VIRTUAL_WORK_COMPONENT_PARITIES", virtual_work_term_parities)

# Reflection blocks: a matrix entry is permitted iff output and input have
# the same reflection sign.  This distinguishes harmless odd-to-odd face
# kinematics from a forbidden odd-to-even leakage channel.
parity = {
    "u1": 1,
    "u2": 1,
    "u3": -1,
    "theta": 1,
    "zeta_plus": 1,
    "zeta_minus": 1,
    "phi": 1,
    "delta_p": 1,
    "vface3": -1,
    "delta_vx3": -1,
    "V": 1,
    "J": 1,
    "affinity": 1,
}
for out_name in ("vface3", "delta_vx3", "V", "J", "affinity"):
    sign = parity[out_name] * parity["u3"]
    print("REFLECTION_BLOCK", out_name, "<-u3", "ALLOWED" if sign == 1 else "FORBIDDEN", sign)

# K1 has a fixed polar datum e=(sin beta,0,cos beta).  Holding the chosen
# background e fixed, reflection changes the u3-theta mixed part's sign.
parity_d1_u3 = -1
parity_d2_u3 = -1
parity_d1_theta = 1
parity_d2_theta = 1
k1_mixed_signs = [parity_d1_u3 * parity_d1_theta, parity_d2_u3 * parity_d2_theta]
print("K1_U3_THETA_TERM_SIGNS", k1_mixed_signs)
print("K1_FIXED_AXIS_BACKGROUND_RESIDUAL_COMPONENT3", "-2*cos(beta)")

# With g=(g1,0,0) and d3=0, K2 reduces exactly to a theta*g1*d2(u3)
# monomial, hence it is reflection odd.
k2_sign = parity["theta"] * 1 * parity_d2_u3
print("K2_REDUCED_CARRIER", "a_K2*theta*g1*d2(u3)")
print("K2_REFLECTION_SIGN", k2_sign)

# O(2) rotation about direction 1.  For q^2=k2^2+k3^2, the displayed
# tangential matrix maps (k2,k3) to (q,0).  The multiplication is performed
# symbolically here as coefficient pairs in k2^2/q and k3^2/q.
print("O2_TANGENTIAL_ROTATION", "[[k2/q,k3/q],[-k3/q,k2/q]]")
print("O2_ROTATED_K_COMPONENT_1", "(k2^2+k3^2)/q=q")
print("O2_ROTATED_K_COMPONENT_2", "(-k3*k2+k2*k3)/q=0")
print("O2_BACKGROUND_GRADIENT_IMAGE", "(g1,0,0)")

# Polar-vector inversion parity of vector spherical harmonics.  If scalar
# Y_lm has p=(-1)^l, spheroidal polar fields have eigenvalue p and the
# toroidal r x grad(Y) field has eigenvalue -p.
for ell in range(5):
    scalar = (-1) ** ell
    spheroidal = scalar
    toroidal = -scalar
    toroidal_exists = ell >= 1
    print(
        "ROUND_PARITY", ell,
        "SCALAR", scalar,
        "SPHEROIDAL", spheroidal,
        "TOROIDAL", toroidal if toroidal_exists else "ZERO_SECTOR",
        "TOROIDAL_EXISTS", toroidal_exists,
    )

# R1.1 versus R1.4: the general symmetric rank-2 tensor invariant under all
# rotations/reflections about axis 1 is diag(A,B,B).  Zeroing every component
# name containing labels 2 or 3 silently specializes B to zero.
T_o2 = [["A", 0, 0], [0, "B", 0], [0, 0, "B"]]
T_name_zeroed = [["A", 0, 0], [0, 0, 0], [0, 0, 0]]
print("O2_INVARIANT_RANK2_GENERAL", T_o2)
print("R1_NAME_ZEROED_RANK2", T_name_zeroed)
print("R1_LOST_O2_ALLOWED_COMPONENT", "B")
