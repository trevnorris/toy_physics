#!/usr/bin/env python3
"""Independent symbolic checks for the S11c-c1 shared-physics review."""

import sympy as sp


def section(title):
    print(f"\n=== {title} ===")


section("1. wavy half-space Dirichlet-to-Neumann kernel")
# The boundary is y=f(x), with y increasing into the exterior.  D0 is the
# flat outgoing conormal DtN: D0(k)=i*q(k).  Expanding a prescribed boundary
# trace g gives G1=-D0 f D0 + f D0**2 - grad(f).grad.
q, qp, kdotkp, k2, kp2, kappa2 = sp.symbols(
    "q_k q_kprime k_dot_kprime k_squared kprime_squared kappa_squared"
)
Dk, Dkp = sp.I * q, sp.I * qp
K_from_terms = sp.expand(-Dk * Dkp + Dkp**2 + (kdotkp - kp2))
K_reduced = sp.factor(
    K_from_terms.subs(qp**2, kappa2 - kp2)
)
K_expected = q * qp + kdotkp - kappa2
rho_bulk, omega_sym, fhat = sp.symbols("rho_bulk omega f_hat")
Z1_kernel = sp.I * rho_bulk * omega_sym * fhat * K_expected / (q * qp)
K_diagonal = sp.simplify(
    K_expected.subs({qp: q, kdotkp: k2, q**2: kappa2 - k2})
)
# Simultaneous substitution is made explicitly because q**2 is nested after
# qp -> q.
K_diagonal = sp.simplify((q**2 + k2 - kappa2).subs(q**2, kappa2 - k2))
K_frozen_output = sp.factor(
    K_expected.subs(q, qp).subs(qp**2, kappa2 - kp2)
)
print("G1 physical-space = -D0(f D0 g) + f D0^2 g - grad(f).grad(g)")
print("G1 kernel before dispersion reduction =", K_from_terms)
print("G1 kernel after q(k')^2=kappa^2-|k'|^2 =", K_reduced)
print("kernel - expected two-momentum form =", sp.simplify(K_reduced - K_expected))
print("impedance Z1(k,k') =", Z1_kernel)
print("rigid-shift diagonal kernel =", K_diagonal)
print("kernel after illegal q(k)->q(k') freeze =", K_frozen_output)


section("2. outward face-normal parity")
s, Wx, Wy, Wz = sp.symbols("s W_x W_y W_z", real=True)
# h_s=s W/2 and F_s=w-h_s.  The orientation law requires n_s=s grad(F_s)/|grad(F_s)|.
grad_F = sp.Matrix([-s * Wx / 2, -s * Wy / 2, -s * Wz / 2, 1])
n_oriented_unnormalized = sp.simplify(s) * grad_F
print("s*grad(F_s) =", list(n_oriented_unnormalized))
print("first-order in-plane tilt =", list(n_oriented_unnormalized[:3, :]))
print("vertical component =", n_oriented_unnormalized[3])
print("tilt under s -> -s residual =", [
    sp.simplify(expr.subs(s, -s) - expr) for expr in n_oriented_unnormalized[:3, :]
])


section("3. flat B0c permeable closure")
Z, LA, LV, LX, rho_m, rho_br, V, mu_theta = sp.symbols(
    "Z Lambda_A Lambda_V Lambda_X rho_m rho_br V mu_theta", nonzero=True
)
p, J, vb = sp.symbols("delta_p J v_bulk")
mu_s = mu_theta / rho_br
affinity = mu_s - p / rho_m
eqs = [
    sp.Eq(J, LA * affinity + LV * V),
    sp.Eq(vb, V + J / rho_m),
    sp.Eq(p, Z * vb),
]
solution = sp.solve(eqs, (p, J, vb), dict=True)[0]
p_expected = sp.factor(
    Z * ((1 + LV / rho_m) * V + LA * mu_theta / (rho_m * rho_br))
    / (1 + LA * Z / rho_m**2)
)
J_expected = sp.factor((solution[J]))
traction_scalar = sp.factor(solution[p] + LX * (mu_s - solution[p] / rho_m))
print("delta_p =", sp.factor(solution[p]))
print("delta_p - [I+Lambda_A Z/rho_m^2]^-1 Z*drive =", sp.factor(solution[p] - p_expected))
print("J =", J_expected)
print("minus-normal traction scalar delta_p+Lambda_X*A =", traction_scalar)
print("d(delta_p)/d(Lambda_X) =", sp.diff(solution[p], LX))
print("d(J)/d(Lambda_X) =", sp.diff(solution[J], LX))
print("d(traction scalar)/d(Lambda_X) =", sp.factor(sp.diff(traction_scalar, LX)))


section("4. global normal scaling is secular")
x, z, eps = sp.symbols("x z epsilon", real=True)
f = sp.Function("f")(x)
u = sp.Function("u")(x, z)
# Under y=(1+eps*f(x))*z, to O(eps):
# partial_y=(1-eps*f) partial_z and partial_x|y=partial_x|z-eps*z*f_x partial_z.
laplacian_first_variation = sp.expand(
    -2 * f * sp.diff(u, z, 2)
    - 2 * z * sp.diff(f, x) * sp.diff(u, x, z)
    - z * sp.diff(f, x, 2) * sp.diff(u, z)
)
k, qwave = sp.symbols("k q", real=True)
plane = sp.exp(sp.I * (k * x + qwave * z))
plane_source = sp.simplify(laplacian_first_variation.subs(u, plane).doit() / plane)
phase_first_variation = sp.diff(
    sp.exp(sp.I * qwave * z * (1 + eps * f)), eps
).subs(eps, 0) / sp.exp(sp.I * qwave * z)
print("delta(Laplacian) u =", laplacian_first_variation)
print("delta(Laplacian) on exp(i k x+i q z), divided by wave =", plane_source)
print("delta outgoing phase / outgoing phase =", sp.simplify(phase_first_variation))
print("coefficients proportional to unbounded z =", [
    sp.factor(-2 * sp.diff(f, x)), sp.factor(-sp.diff(f, x, 2)), sp.I * qwave * f
])


section("5. convective normal-wavenumber correction near grazing")
omega, c, v0, q0, dq, kabs = sp.symbols(
    "omega c v0 q0 delta_q k_abs", nonzero=True
)
disp = (omega - v0 * (q0 + dq))**2 - c**2 * (kabs**2 + (q0 + dq)**2)
# Linearize in v0 and dq, treating dq=O(v0), then impose rest dispersion.
lin = sp.expand(disp).subs({v0**2: 0, dq**2: 0, v0 * dq: 0})
lin_on_shell = sp.expand(lin.subs(kabs**2, omega**2 / c**2 - q0**2))
dq_solution = sp.solve(sp.Eq(lin_on_shell, 0), dq)[0]
print("linearized convective dispersion on rest shell =", sp.factor(lin_on_shell))
print("delta_q =", sp.factor(dq_solution))
print("relative delta_q/q0 =", sp.factor(dq_solution / q0))
print("timescale ratio =", sp.Abs(q0 * v0 / omega))
print("normal-root relative-correction ratio =", sp.Abs(omega * v0 / (c**2 * q0)))
print("exact dispersion at rest-frame grazing =", sp.factor(
    ((omega - v0 * q)**2 - c**2 * (kabs**2 + q**2)).subs(kabs**2, omega**2 / c**2)
))


section("6. entrywise real part is not the Hermitian part")
A = sp.Matrix([[0, sp.I], [0, 0]])
H = sp.simplify((A + A.conjugate().T) / 2)
entrywise_real = A.applyfunc(sp.re)
vec = sp.Matrix([1, -sp.I])
quadratic_H = sp.simplify((vec.conjugate().T * H * vec)[0])
quadratic_A_real = sp.simplify(sp.re((vec.conjugate().T * A * vec)[0]))
print("A =", A)
print("entrywise Re(A) =", entrywise_real)
print("Hermitian part (A+A^dagger)/2 =", H)
print("x^dagger H x for x=(1,-i) =", quadratic_H)
print("Re(x^dagger A x) =", quadratic_A_real)
print("entrywise Re(A) - Hermitian part =", entrywise_real - H)


section("7. traction-derived slab work versus outgoing bulk flux")
pamp, Vamp = sp.symbols("p_amp V_amp")
traction_normal = -pamp
slab_work = sp.expand(traction_normal * Vamp)
outgoing_flux = sp.expand(pamp * Vamp)
wrong_traction_work = sp.expand((-traction_normal) * Vamp)
print("traction normal component for t=-p*n =", traction_normal)
print("slab work derived from t dot v_face =", slab_work)
print("outgoing bulk flux amplitude =", outgoing_flux)
print("correct energy-sum residual =", sp.expand(slab_work + outgoing_flux))
print("energy-sum residual after reversing traction =", sp.expand(wrong_traction_work + outgoing_flux))
print("typed +pV slab operand minus outgoing flux =", sp.expand(pamp * Vamp - outgoing_flux))


section("8. advected density is absent from mu_s at linear wave order")
rho0, density_jet, mu1, eps_wave, sigma = sp.symbols(
    "rho0 density_jet mu1 epsilon sigma", nonzero=True
)
rho_advected = rho0 - eps_wave * sigma * density_jet
mu_theta_linear = eps_wave * mu1
mu_s_advected = sp.series(mu_theta_linear / rho_advected, eps_wave, 0, 3).removeO()
mu_s_unadvected = mu_theta_linear / rho0
print("mu_s with advected rho, through epsilon^2 =", sp.expand(mu_s_advected))
print("mu_s with background rho held =", mu_s_unadvected)
print("difference =", sp.factor(mu_s_advected - mu_s_unadvected))
print("coefficient of epsilon^1 in difference =", sp.expand(mu_s_advected - mu_s_unadvected).coeff(eps_wave, 1))
print("coefficient of epsilon^2 in difference =", sp.expand(mu_s_advected - mu_s_unadvected).coeff(eps_wave, 2))


section("9. sign-definiteness is not decidable from a first-order truncation on a lossless nullspace")
a, b, eta = sp.symbols("a b eta", positive=True)
H_exact = sp.Matrix([[a, eta * b], [eta * b, eta**2 * b**2 / a]])
H_first = sp.Matrix([[a, eta * b], [eta * b, 0]])
print("exact positive-semidefinite completion =", H_exact)
print("det(exact completion) =", sp.factor(H_exact.det()))
print("first-shape-order truncation =", H_first)
print("det(first-order truncation) =", sp.factor(H_first.det()))
print("order of the negative determinant in eta =", sp.Poly(H_first.det(), eta).degree())
