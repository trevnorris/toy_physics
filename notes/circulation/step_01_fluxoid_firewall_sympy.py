#!/usr/bin/env python3
"""
Step 1 SymPy audit: gauge-invariant circulation / fluxoid constraint.

This script audits local gauge cancellation, single-valued-gauge loop
integrals, and the integer phase winding implied by single-valued psi.

Run:
    python step_01_fluxoid_firewall_sympy.py
"""
import sympy as sp

# Branch and physical symbols.
n = sp.symbols("n", integer=True)
hbar, m, e_star, Z_int = sp.symbols("hbar m e_star Z_int", positive=True)
eta_Q = sp.symbols("eta_Q")  # later restricted to +/-1
q_star = eta_Q * e_star
q_eff = q_star / sp.sqrt(Z_int)

# Gauge-covariant phase-gradient component along a loop.
dtheta, dchi, A_i = sp.symbols("dtheta dchi A_i")
P_i = dtheta - q_star / hbar * A_i

# Gauge transformation: theta -> theta + q_* chi/hbar, A_i -> A_i + partial_i chi.
P_i_prime = (dtheta + q_star / hbar * dchi) - q_star / hbar * (A_i + dchi)
gauge_invariance_residual = sp.simplify(P_i_prime - P_i)

# Explicit loop audit. Let t parameterize a small circle linking the mouth.
# A single-valued gauge function has zero winding around the loop.
t = sp.symbols("t", real=True)
Delta_theta = sp.symbols("Delta_theta", real=True)
theta0, chi0, chi_c, chi_s, chi_w = sp.symbols("theta0 chi0 chi_c chi_s chi_w", real=True)
chi_single = chi0 + chi_c * sp.cos(t) + chi_s * sp.sin(t)
single_valued_chi_winding = sp.simplify(
    sp.integrate(sp.diff(chi_single, t), (t, 0, 2 * sp.pi))
)

# psi = sqrt(rho) exp(i theta) is single-valued only when exp(i Delta_theta)=1.
# Equivalently, Re/Im give cos(Delta_theta)=1 and sin(Delta_theta)=0.
phase_real_solutions = sp.solveset(
    sp.Eq(sp.cos(Delta_theta), 1), Delta_theta, domain=sp.S.Reals
)
phase_imag_solutions = sp.solveset(
    sp.Eq(sp.sin(Delta_theta), 0), Delta_theta, domain=sp.S.Reals
)
phase_winding_solutions = sp.Intersection(phase_real_solutions, phase_imag_solutions)

# A pure-gauge A_t = d_t chi_single can appear in both theta and A_t. The
# eta_Q-dependent pieces cancel before the loop integral is evaluated. The
# endpoint phase jump Delta_theta is left symbolic until the single-valuedness
# result is applied.
A_t_single = sp.diff(chi_single, t)
theta_path = theta0 + Delta_theta * t / (2 * sp.pi) + q_star / hbar * (chi_single - chi_single.subs(t, 0))
P_t_single = sp.simplify(sp.diff(theta_path, t) - q_star / hbar * A_t_single)
fluxoid_integral_delta = sp.simplify(sp.integrate(P_t_single, (t, 0, 2 * sp.pi)))
fluxoid_integral = fluxoid_integral_delta.subs(Delta_theta, 2 * sp.pi * n)

# A multi-valued chi_w t would shift the winding. This is not an allowed small
# single-valued gauge transformation unless chi_w = 0.
large_gauge_winding = sp.simplify(
    sp.integrate(sp.diff(chi_w * t, t), (t, 0, 2 * sp.pi))
)

Gamma_n = sp.simplify(hbar / m * fluxoid_integral)
kappa = sp.simplify(2 * sp.pi * hbar / m)

assert gauge_invariance_residual == 0
assert single_valued_chi_winding == 0
assert phase_winding_solutions.contains(2 * sp.pi * n) is sp.S.true
assert phase_winding_solutions.contains(sp.pi) is sp.S.false
assert P_t_single == Delta_theta / (2 * sp.pi)
assert fluxoid_integral_delta == Delta_theta
assert fluxoid_integral == 2 * sp.pi * n
assert sp.simplify(Gamma_n - kappa * n) == 0

print("STEP 1: FLUXOID LOOP AUDIT")
print("Gauge-invariant loop integrand P_i = dtheta - q_*/hbar A_i")
print("Gauge residual P_i' - P_i =", gauge_invariance_residual)
print()
print("Single-valued chi(t) =", chi_single)
print("Integral d chi around loop =", single_valued_chi_winding)
print("psi single-valued equations: cos(Delta_theta)=1 and sin(Delta_theta)=0")
print("Allowed Delta_theta set =", phase_winding_solutions)
print("Pure-gauge loop integrand after cancellation P_t =", P_t_single)
print("Fluxoid integral before quantization =", fluxoid_integral_delta)
print("Fluxoid integral after Delta_theta=2*pi*n =", fluxoid_integral)
print("Hydrodynamic circulation Gamma_n =", Gamma_n)
print("kappa =", kappa)
print()
print("Large/non-single-valued chi_w t winding would be", large_gauge_winding)
print("q_* =", q_star)
print("q_eff bookkeeping label =", q_eff)
print("Verdict: single-valued psi gives integer phase winding; single-valued chi adds no winding.")
