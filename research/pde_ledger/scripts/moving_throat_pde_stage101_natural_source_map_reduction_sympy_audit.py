#!/usr/bin/env python3
"""
Stage 101 SymPy audit: natural source-map reduction (mhat_0 = 1 branch).

Carry-forward annotations (per [[batch-IV1-paper-alignment]] Cluster B direction
(c)): paper card `\\stagefield{Checks}` items
  (2) "higher odd terms begin beyond the point-particle 2.5PN coefficient"
      — anchored at stage 102 (higher_odd_irrelevance); the omega^7 series
      extension and tauQ-derivative gates verify the first non-(omega^2,
      omega^5) term sits at omega^7, in both engines.
  (3) "outgoing l=2 DtN fingerprint vs normalized z=omega a/c_s expansion"
      — anchored at stage 097 (single_normalization_defect); the chi_Q = 1
      identification by DtN comparison lives there. Notes line 41-51 of this
      stage explicitly attributes the canonical compact branch's chi_Q = 1
      to that upstream stage.
This stage owns Check (1) — the factorization mhat_0^2 chi_Q N_Q = 1 keeping
source, conservative, and outgoing factors separate — verified below via
substantive input-equation anchors.
"""
from __future__ import annotations
import sympy as sp

mhat0, chiQ, NQ, DeltaQ = sp.symbols('mhat0 chiQ NQ DeltaQ', positive=True, real=True)

# exact factorized condition
sol_NQ = sp.solve(sp.Eq(mhat0**2 * chiQ * NQ, 1), NQ)[0]
print('NQ from exact factorized odd normalization =', sol_NQ)
print('point-particle natural branch mhat0->1 gives =', sp.simplify(sol_NQ.subs(mhat0, 1)))
print('canonical compact outgoing branch chiQ=1 gives =', sp.simplify(sol_NQ.subs({mhat0:1, chiQ:1})))

# Anchor the reductions to the INPUT factorization mhat0^2 * chiQ * NQ = 1.
# Substituting the proposed solved NQ on the natural source-map branch must zero the residual.
assert sp.simplify((mhat0**2 * chiQ * NQ - 1).subs({mhat0: 1, NQ: 1/chiQ})) == 0, \
    'point-particle natural branch reduction NQ = 1/chiQ failed against input factorization'
assert sp.simplify((mhat0**2 * chiQ * NQ - 1).subs({mhat0: 1, chiQ: 1, NQ: 1})) == 0, \
    'canonical compact outgoing branch NQ = 1 failed against input factorization'

expr_delta = sp.simplify((1/(1+DeltaQ)) - 1)
series_delta = sp.expand(sp.series(expr_delta, DeltaQ, 0, 3).removeO())
print('NQ - 1 in terms of DeltaQ =', expr_delta)
print('small-DeltaQ series =', series_delta)
print('check exact replacement chiQ=1+DeltaQ:', sp.simplify(sol_NQ.subs({mhat0:1, chiQ:1+DeltaQ}) - 1/(1+DeltaQ)))

# Anchor exact NQ = 1/(1 + DeltaQ) on the natural source-map branch to the input factorization,
# with chiQ = 1 + DeltaQ.
assert sp.simplify((mhat0**2 * chiQ * NQ - 1).subs({mhat0: 1, chiQ: 1 + DeltaQ, NQ: 1/(1 + DeltaQ)})) == 0, \
    'exact replacement NQ = 1/(1+DeltaQ) failed against input factorization'
# Confirm the small-DeltaQ linearization matches the paper's stated form NQ - 1 = -DeltaQ + O(DeltaQ^2).
# The script's series_delta is the order-2 truncation of -DeltaQ/(1+DeltaQ).
assert sp.expand(series_delta - (-DeltaQ + DeltaQ**2)) == 0, \
    'small-DeltaQ series does not match -DeltaQ + DeltaQ**2'

print('\nSTAGE 101 AUDIT PASSED')
