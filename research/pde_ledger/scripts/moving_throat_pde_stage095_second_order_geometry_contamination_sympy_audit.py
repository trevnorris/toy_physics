#!/usr/bin/env python3
"""
Stage 78 SymPy audit

Checks that the first nonzero geometry contamination produced by weak scalar/l=2
mixing enters at O(chi^2) and extracts the exact low-frequency coefficients.
"""

from __future__ import annotations
import sympy as sp

chi, M0, G0, G2, G4, w = sp.symbols('chi M0 G0 G2 G4 w', real=True)
G0 = sp.symbols('G0', nonzero=True, real=True)
Kstat, Kpole, OQ = sp.symbols('Kstat Kpole OQ', nonzero=True, real=True)

Dg = G0 + G2*w**2 + G4*w**4
corr = sp.expand(sp.series(-chi**2 * M0**2 / Dg, w, 0, 6).removeO())
print('Schur-complement correction =', corr)

K0corr = sp.simplify(corr.subs(w, 0))
K2corr = sp.simplify(sp.expand(corr).coeff(w, 2))
K4corr = sp.simplify(sp.expand(corr).coeff(w, 4))
print('K0corr =', K0corr)
print('K2corr =', K2corr)
print('K4corr =', K4corr)

assert sp.factor(K0corr).has(chi**2)
assert sp.factor(K2corr).has(chi**2)
assert sp.factor(K4corr).has(chi**2)

# Dimensionless contamination numbers from Stage 75.
eps2 = sp.simplify(OQ**2 * K2corr / Kpole)
eps4 = sp.simplify(OQ**4 * K4corr / Kpole)
print('eps2 =', eps2)
print('eps4 =', eps4)

# Pole fraction expansion.
cpole = sp.simplify((1 + eps4) / (4 * (1 + eps2)**2))
cpole_series = sp.expand(sp.series(cpole, chi, 0, 3).removeO())
print('c_pole =', cpole_series)

delta = sp.simplify(cpole_series - sp.Rational(1, 4))
print('delta c_pole =', delta)
assert sp.simplify(delta.subs(chi, 0)) == 0

# Verify first nonzero term is O(chi^2).
chi1 = sp.simplify(sp.diff(cpole, chi).subs(chi, 0))
chi2 = sp.simplify(sp.diff(cpole, chi, 2).subs(chi, 0))
print('d c_pole / dchi |0 =', chi1)
print('d^2 c_pole / dchi^2 |0 =', chi2)
assert chi1 == 0

print('\nStage 78 theorem verified: contamination begins at O(chi^2).')
