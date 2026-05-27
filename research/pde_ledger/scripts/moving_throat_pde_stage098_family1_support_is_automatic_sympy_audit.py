#!/usr/bin/env python3
"""
Stage 098 SymPy audit: actual isotropic branch support demand is automatic for any
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
# zeta_max^(F1) is an external carry-forward; not derived in this stage.
# Source: notes/stages/moving_throat_pde_stage098_family1_support_is_automatic.md (Family-1 specialization).
zmax_F1 = sp.N('2.46752922945601')
zeta_edge_F1 = sp.N(zeta_edge.subs(zmax, zmax_F1), 30)
gap_F1 = sp.N(zmax_F1 - zeta_edge_F1, 30)
assert gap_F1 > 0
# Numeric pin (matches Mathematica's expectApprox targets):
zeta_edge_F1_target = sp.Float('0.456730991107963169017835980412', 30)
gap_F1_target = sp.Float('2.01079823834804688464927835412', 30)
assert abs(zeta_edge_F1 - zeta_edge_F1_target) < sp.Float('1e-15', 30)
assert abs(gap_F1 - gap_F1_target) < sp.Float('1e-15', 30)

print('zeta_req(eps) =', zeta_req)
print('d zeta_req / d eps =', dz)
print('zeta_edge =', zeta_edge)
print('zmax - zeta_edge =', gap_factored)
print('Family-1 zeta_edge =', zeta_edge_F1)
print('Family-1 margin =', gap_F1)
print('\nSTAGE 098 AUDIT PASSED')
