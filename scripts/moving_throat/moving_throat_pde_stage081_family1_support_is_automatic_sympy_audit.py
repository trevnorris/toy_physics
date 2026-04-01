#!/usr/bin/env python3
"""
Stage 81 SymPy audit: actual isotropic branch support demand is automatic for any
explicit family with zeta_max > 1 on the admissible blocked interval.
"""
from __future__ import annotations
import sympy as sp

eps, zmax = sp.symbols('eps zmax', positive=True, real=True)
rho = sp.Rational(4, 3)
zeta_req = sp.simplify((rho - 1) / (1 - eps * (2 - rho)))
assert sp.simplify(zeta_req - 1 / (3 - 2 * eps)) == 0

dz = sp.simplify(sp.diff(zeta_req, eps))
assert sp.simplify(dz - 2 / (3 - 2 * eps)**2) == 0

# Worst-case admissible blocked value at eps = 1/zmax.
zeta_edge = sp.simplify(zeta_req.subs(eps, 1 / zmax))
gap = sp.simplify(zmax - zeta_edge)
gap_factored = sp.factor(gap)
assert sp.simplify(gap_factored - 3 * zmax * (zmax - 1) / (3 * zmax - 2)) == 0

# Family-1 specialization.
zmax_F1 = sp.N('2.46752922945601')
zeta_edge_F1 = sp.N(zeta_edge.subs(zmax, zmax_F1), 30)
gap_F1 = sp.N(zmax_F1 - zeta_edge_F1, 30)
assert gap_F1 > 0

print('zeta_req(eps) =', zeta_req)
print('d zeta_req / d eps =', dz)
print('zeta_edge =', zeta_edge)
print('zmax - zeta_edge =', gap_factored)
print('Family-1 zeta_edge =', zeta_edge_F1)
print('Family-1 margin =', gap_F1)
print('\nSTAGE 81 AUDIT PASSED')
