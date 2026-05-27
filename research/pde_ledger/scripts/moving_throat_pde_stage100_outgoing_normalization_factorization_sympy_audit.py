#!/usr/bin/env python3
"""
Stage 100 SymPy audit: outgoing-normalization factorization mhat_0^2 chi_Q N_Q = 1.

Carry-forward annotations (per [[batch-IV1-paper-alignment]] Cluster B direction
(c)): paper card `\\stagefield{Checks}` items
  (ii) "higher odd terms begin beyond the point-particle 2.5PN coefficient"
       — verified at stage 102 (higher_odd_irrelevance) via omega^7 series
       extension and the tauQ-derivative gate; this stage's series truncates
       at omega^5 by design.
  (iii) "outgoing l=2 DtN fingerprint against the normalized z=omega a/c_s
       expansion" — the chi_Q = 1 identification by DtN comparison lives at
       stage 097 (single_normalization_defect); this stage carries chi_Q as
       a free real symbol.
This stage owns Check (i) — the `mhat_0^2 chi_Q N_Q = 1` closure — which is
derived below by imposing the observable-side condition `mhat_0^2 * Gamma_5 =
Gamma_5_target` on the script-derived Gamma_5(K_0, chi_Q, Omega).
"""
from __future__ import annotations
import sympy as sp

G, c = sp.symbols('G c', positive=True, real=True)
Omega = sp.symbols('Omega', positive=True, real=True)
K0, mhat0 = sp.symbols('K0 mhat0', positive=True, real=True)
# chi_Q is a free real parameter at this stage. Its pinning to 1 by the l=2
# DtN fingerprint is upstream (stage 097). Declaring it positive would falsely
# constrain the factorization assertion below.
chiQ = sp.symbols('chiQ', real=True)
NQ = sp.symbols('N_Q', positive=True, real=True)
omega = sp.symbols('omega', real=True)

sigma_can = sp.Rational(9, 8) / Omega**5
Y = sp.Rational(3,4) + sp.Rational(1,4) / (1 - omega**2/Omega**2 - sp.I*chiQ*sigma_can*omega**5)
Yser = sp.expand(sp.series(Y, omega, 0, 6).removeO())

K2 = sp.simplify(K0 * Yser.coeff(omega, 2))
K4 = sp.simplify(K0 * Yser.coeff(omega, 4))
Gamma5 = sp.simplify(sp.im(Yser.coeff(omega, 5)) * K0)

K0_t = sp.simplify(64*G*Omega**5/(45*c**5))
K2_t = sp.simplify(K0_t/(4*Omega**2))
K4_t = sp.simplify(K0_t/(4*Omega**4))
Gamma5_t = sp.simplify(2*G/(5*c**5))
NQ_derived = sp.simplify(K0 / K0_t)

print('Yhat_Q^ret series =', Yser)
print('K2 =', K2)
print('K4 =', K4)
print('Gamma5 =', Gamma5)
print('NQ =', NQ_derived)
print('K2/K2_target - NQ =', sp.simplify(K2/K2_t - NQ_derived))
print('K4/K4_target - NQ =', sp.simplify(K4/K4_t - NQ_derived))
print('Gamma5/Gamma5_target - chiQ*NQ =', sp.simplify(Gamma5/Gamma5_t - chiQ*NQ_derived))
assert sp.simplify(K2/K2_t - NQ_derived) == 0
assert sp.simplify(K4/K4_t - NQ_derived) == 0
assert sp.simplify(Gamma5/Gamma5_t - chiQ*NQ_derived) == 0

# Substantive closure: impose the observable-side condition
#   mhat_0^2 * Gamma_5 = Gamma_5_target
# (the point-particle 2.5PN coefficient with source-map factor mhat_0^2 matches
# the GR target). Substituting the script-derived Gamma_5 = chi_Q N_Q Gamma_5_target,
# the residual factorizes as Gamma_5_target * (mhat_0^2 chi_Q N_Q - 1).
# This is non-tautological because if the series derivation of Gamma_5 had a
# wrong factor, the ratio below would not equal (mhat_0^2 chi_Q N_Q - 1).
closure_residual = sp.simplify(mhat0**2 * Gamma5 - Gamma5_t)
closure_ratio = sp.simplify(closure_residual / Gamma5_t)
closure_check = sp.simplify(closure_ratio - (mhat0**2 * chiQ * NQ_derived - 1))
print('closure_residual / Gamma5_target =', closure_ratio)
print('closure_ratio - (mhat0^2 chi_Q N_Q - 1) =', closure_check)
assert closure_check == 0

# Headline ledger entry: on the closure (mhat_0^2 Gamma_5 = Gamma_5_target),
# the factorization mhat_0^2 chi_Q N_Q = 1 is forced.
print('\nClosure ledger: mhat_0^2 * Gamma_5 = Gamma_5_target  <=>  mhat_0^2 chi_Q N_Q = 1')
print('STAGE 100 AUDIT PASSED')
