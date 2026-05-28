
#!/usr/bin/env python3
"""
moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.py

SymPy audit for the outlet-consistent mouth closure:
    M_s = 4 Sigma_m,   M_q = -Sigma_m.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

Pi, Sigma_m = sp.symbols("Pi Sigma_m", positive=True, real=True)

def S(Pi, kappa):
    return sp.simplify(
        Pi * (kappa * sp.tanh(kappa) + Pi * (sp.exp(-Pi) / sp.cosh(kappa) - 1))
        / ((1 - sp.exp(-Pi)) * (kappa**2 - Pi**2))
    )

banner("STAGE 135 — OUTLET-CONSISTENT ONE-PARAMETER CLOSURE")
S_q = sp.simplify(S(Pi, sp.pi/2))

# Step 1: outlet substitution check.
# Build the generic mouth-gain law M_s + M_q*S_q with independent symbols, then
# substitute the outlet-consistent reduction M_s = 4*Sigma_m, M_q = -Sigma_m and
# verify it reduces to Sigma_m*(4 - S_q).
M_s_sym, M_q_sym = sp.symbols("M_s_sym M_q_sym", real=True)
generic_law = M_s_sym + M_q_sym * S_q
reduced_law = generic_law.subs({M_s_sym: 4*Sigma_m, M_q_sym: -Sigma_m})
expected_law = Sigma_m * (4 - S_q)
residual_sub = sp.simplify(reduced_law - expected_law)
print("M_s + M_q*S_q - Sigma_m*(4 - S_q) =", residual_sub)
assert residual_sub == 0, f"Outlet substitution failed: residual = {residual_sub}"

Pi_eq = sp.simplify(Sigma_m * (4 - S_q))
print("Pi =")
sp.pprint(Pi_eq)

Pi_star = sp.Float("1.50882951349316")
S_star = sp.N(S_q.subs(Pi, Pi_star), 30)
Sigma_var = sp.symbols("Sigma_var", positive=True, real=True)
Sigma_star = sp.N(sp.solve(sp.Eq(Pi_star, Sigma_var * (4 - S_star)), Sigma_var)[0], 30)
M_s_star = sp.N(4 * Sigma_star, 30)
M_q_star = sp.N(-Sigma_star, 30)

print("S_q(Pi_star) =", S_star)

# Step 2: assert the inequality 0 < S_q(Pi_*) < 1 (paper deliverable iii).
s_in_range = bool(0 < S_star < 1)
print("0 < S_q(Pi_star) < 1 ->", s_in_range)
assert s_in_range, f"S_q(Pi_*) out of range: S_star = {S_star}"

print("Sigma_m^* =", Sigma_star)
print("M_s^* =", M_s_star)
print("M_q^* =", M_q_star)

# Step 3: numerical anchors against notes-stated values (paper deliverables iv-vi).
Sigma_target = sp.Float("0.451485277739090", 30)
M_s_target = sp.Float("1.80594111095636", 30)
M_q_target = sp.Float("-0.451485277739090", 30)

assert abs(Sigma_star - Sigma_target) < sp.Float("1e-12", 30), \
    f"Sigma_m^* mismatch: {Sigma_star} vs {Sigma_target}"
assert abs(M_s_star - M_s_target) < sp.Float("1e-11", 30), \
    f"M_s^* mismatch: {M_s_star} vs {M_s_target}"
assert abs(M_q_star - M_q_target) < sp.Float("1e-12", 30), \
    f"M_q^* mismatch: {M_q_star} vs {M_q_target}"
print("Sigma_m^*, M_s^*, M_q^* anchored to notes values within tolerance.")

# Step 4: mixed-lane correction M_q^* * S_q(Pi_*) (notes line 125-127 ref value).
mixed_correction = sp.N(M_q_star * S_star, 30)
mixed_target = sp.Float("-0.297111597463199", 30)
print("M_q^* * S_q(Pi_*) =", mixed_correction)
assert abs(mixed_correction - mixed_target) < sp.Float("1e-11", 30), \
    f"mixed-lane correction mismatch: {mixed_correction} vs {mixed_target}"

# Numerical sanity probe (was the original tautological closure residual).
residual = sp.N(Pi_star - Sigma_star * (4 - S_star), 30)
print("Pi_star - Sigma_star*(4 - S_star) =", residual)
