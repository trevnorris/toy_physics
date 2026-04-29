#!/usr/bin/env python3
"""
Step 6 SymPy audit: minimal mixed-sector plumbing closure and sign inheritance.

The static transfer factor is used as a nonnegative magnitude factor. It does
not choose the circulation sign. Electric branch signs are intentionally not
included in this reduced current-closure audit.

Run:
    python step_06_mixed_sector_plumbing_closure_sympy.py
"""
import sympy as sp

# Reduced Maxwell/mixed one-port symbols.
OmegaU, OmegaW = sp.symbols("Omega_U Omega_W", positive=True)
Rmix, gU, gW = sp.symbols("R_mix g_U g_W", real=True)
Delta = OmegaU**2 * OmegaW**2 - Rmix**2
Transfer0 = sp.factor((OmegaU**2 * gW + Rmix * gU)**2 / Delta**2)

transfer_num, transfer_den = sp.together(Transfer0).as_numer_denom()
num_base, num_exp = sp.factor(transfer_num).as_base_exp()
den_base, den_exp = sp.factor(transfer_den).as_base_exp()
Transfer_base = sp.simplify((OmegaU**2 * gW + Rmix * gU) / Delta)
tau = sp.symbols("tau", real=True)

# Generic mixed/current closure after transfer. Lambda_A is real here: its sign
# is the unresolved PDE-side plumbing question, not a symbol assumption.
K, d = sp.symbols("K d", positive=True)
Lambda1, Lambda2 = sp.symbols("Lambda1 Lambda2", real=True)
n1, n2 = sp.symbols("n1 n2", real=True)
sigma1, sigma2 = sp.symbols("sigma1 sigma2", real=True)

F_mixed = sp.simplify(
    -K * Transfer0 * Lambda1 * Lambda2 * sigma1 * sigma2 * n1 * n2 / d**4
)

# Facing-mouth branch under identical passive mixed closure.
facing_subs = {sigma1: 1, sigma2: -1, n1: 1}
F_facing_opposite_local = sp.simplify(F_mixed.subs({**facing_subs, n2: -1}))
F_facing_same_local = sp.simplify(F_mixed.subs({**facing_subs, n2: 1}))

assert num_exp == 2
assert den_exp == 2
assert sp.simplify(transfer_num - (OmegaU**2 * gW + Rmix * gU)**2) == 0
assert sp.simplify(transfer_den - Delta**2) == 0
assert sp.simplify(Transfer0 - Transfer_base**2) == 0
assert (tau**2).is_nonnegative
assert F_facing_opposite_local == -K * Lambda1 * Lambda2 * Transfer0 / d**4
assert F_facing_same_local == K * Lambda1 * Lambda2 * Transfer0 / d**4

print("STEP 6: MIXED-SECTOR PLUMBING CLOSURE")
print("Stage-4-style static transfer factor:")
print("N_0 =", Transfer0)
print("Numerator base/exponent =", num_base, num_exp)
print("Denominator base/exponent =", den_base, den_exp)
print("N_0 square base =", Transfer_base)
print("Symbolic square certificate: tau real => tau^2 nonnegative is", (tau**2).is_nonnegative)
print("For real nonzero Delta, N_0 is nonnegative.")
print()
print("Minimal mixed/current closure force with transfer magnitude included:")
print("F_d =", F_mixed)
print()
print("Facing-mouth identical passive branch:")
print("  opposite local swirl: F_d =", F_facing_opposite_local)
print("  same local swirl:     F_d =", F_facing_same_local)
print()
print("Attraction for opposite local swirl requires Lambda1*Lambda2*N_0 > 0.")
print("Verdict: Lambda_A signs remain a plumbing-law assumption, not a SymPy assumption.")
