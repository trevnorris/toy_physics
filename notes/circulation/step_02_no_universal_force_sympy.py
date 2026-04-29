#!/usr/bin/env python3
"""
Step 2 SymPy audit: the 3D fluxoid constraint alone does not determine a
radial force. A current/moment closure and an ensemble choice are additional
assumptions.

Run:
    python step_02_no_universal_force_sympy.py
"""
import sympy as sp

# Separation and a positive mutual-inductance prefactor, M(d)=B/d^3.
d, B = sp.symbols("d B", positive=True)
n1, n2 = sp.symbols("n1 n2", real=True)
sigma1, sigma2 = sp.symbols("sigma1 sigma2", real=True)  # orientation signs +/-1
alpha1, alpha2 = sp.symbols("alpha1 alpha2", real=True)  # plumbing/current coefficients

# Supply an explicit current-like closure and derive the leading coaxial force
# from M(d), rather than writing the force sign directly.
I1_closure = alpha1 * sigma1 * n1
I2_closure = alpha2 * sigma2 * n2
M_of_d = B / d**3
U_current_closure = -I1_closure * I2_closure * M_of_d
F_current_closure = sp.simplify(-sp.diff(U_current_closure, d))

# Two closures with identical n but opposite plumbing sign on mouth 2.
unit_case = {alpha1: 1, sigma1: 1, sigma2: 1, n1: 1, n2: 1}
F_same_plumbing = sp.simplify(F_current_closure.subs({**unit_case, alpha2: 1}))
F_flipped_plumbing = sp.simplify(F_current_closure.subs({**unit_case, alpha2: -1}))

# Fixed-current mechanics from a Legendre transform.
L1, L2, M_sym = sp.symbols("L1 L2 M", positive=True)
I1, I2, Phi1, Phi2 = sp.symbols("I1 I2 Phi1 Phi2", real=True)
Lmat = sp.Matrix([[L1, M_sym], [M_sym, L2]])
Ivec = sp.Matrix([I1, I2])
Phi_from_I = Lmat * Ivec

field_energy_I = sp.simplify(sp.Rational(1, 2) * (Ivec.T * Lmat * Ivec)[0])
fixed_current_potential = sp.simplify(field_energy_I - (Ivec.T * Phi_from_I)[0])

field_cross_coeff = sp.diff(field_energy_I, I1, I2)
current_potential_cross_coeff = sp.diff(fixed_current_potential, I1, I2)

field_energy_gradient_if_used_as_potential = sp.simplify(
    -sp.diff(field_energy_I.subs(M_sym, M_of_d), d)
)
fixed_current_force_from_legendre = sp.simplify(
    -sp.diff(fixed_current_potential.subs(M_sym, M_of_d), d)
)

# Fixed-flux energy is obtained from the inverse inductance matrix. This is a
# separate ensemble, not the same expression with a hand-flipped sign.
Phivec = sp.Matrix([Phi1, Phi2])
fixed_flux_energy = sp.simplify(sp.Rational(1, 2) * (Phivec.T * Lmat.inv() * Phivec)[0])
Delta_L = L1 * L2 - M_sym**2
fixed_flux_energy_expected = sp.simplify(
    (L2 * Phi1**2 - 2 * M_sym * Phi1 * Phi2 + L1 * Phi2**2) / (2 * Delta_L)
)
weak_fixed_flux_slope = sp.simplify(sp.diff(fixed_flux_energy, M_sym).subs(M_sym, 0))
weak_fixed_flux_force = sp.simplify(-weak_fixed_flux_slope * sp.diff(M_of_d, d))
fixed_flux_force_exact = sp.factor(
    -sp.diff(fixed_flux_energy.subs(M_sym, M_of_d), d)
)
fixed_flux_force_series = sp.series(fixed_flux_force_exact, B, 0, 3).removeO()
fixed_flux_force_series_expected = (
    -3 * B * Phi1 * Phi2 / (L1 * L2 * d**4)
    + 3 * B**2 * (L1 * Phi2**2 + L2 * Phi1**2) / (L1**2 * L2**2 * d**7)
)

assert F_current_closure == -3 * B * alpha1 * alpha2 * sigma1 * sigma2 * n1 * n2 / d**4
assert F_same_plumbing == -3 * B / d**4
assert F_flipped_plumbing == 3 * B / d**4
assert field_cross_coeff == M_sym
assert current_potential_cross_coeff == -M_sym
assert sp.simplify(field_energy_gradient_if_used_as_potential + fixed_current_force_from_legendre) == 0
assert sp.simplify(fixed_flux_energy - fixed_flux_energy_expected) == 0
assert weak_fixed_flux_slope == -Phi1 * Phi2 / (L1 * L2)
assert sp.simplify(fixed_flux_force_series - fixed_flux_force_series_expected) == 0

print("STEP 2: NO UNIVERSAL RADIAL FORCE FROM FLUXOID ALONE")
print("Closure currents:")
print("I1 =", I1_closure)
print("I2 =", I2_closure)
print("M(d) =", M_of_d)
print("Derived fixed-current potential U_I =", U_current_closure)
print("Derived force F_d = -dU_I/dd =", F_current_closure)
print("F_d < 0 means attraction.")
print()
print("Same n1=n2=+1 attracts with alpha1*alpha2>0:", F_same_plumbing)
print("Same n1=n2=+1 repels with alpha1*alpha2<0:", F_flipped_plumbing)
print()
print("Inductance matrix L =", Lmat)
print("Field energy at fixed current coordinates E_B =", field_energy_I)
print("Fluxes Phi = L I =", Phi_from_I)
print("Legendre fixed-current potential G_I = E_B - I.Phi =", fixed_current_potential)
print("cross coefficient in E_B =", field_cross_coeff)
print("cross coefficient in G_I =", current_potential_cross_coeff)
print("Force from G_I with M(d)=B/d^3 =", fixed_current_force_from_legendre)
print("Naively using E_B as a potential would give", field_energy_gradient_if_used_as_potential)
print()
print("Fixed-flux energy U_Phi = 1/2 Phi^T L^-1 Phi =", fixed_flux_energy)
print("weak-coupling dU_Phi/dM at M=0 =", weak_fixed_flux_slope)
print("weak-coupling fixed-flux force =", weak_fixed_flux_force)
print("exact fixed-flux force with M(d)=B/d^3 =", fixed_flux_force_exact)
print("fixed-flux force series through B^2 =", fixed_flux_force_series)
print("At weak coupling the fixed-flux leading sign matches fixed current after I ~= L^-1 Phi.")
print("Verdict: fluxoid n fixes holonomy; force sign requires a closure and stated ensemble.")
