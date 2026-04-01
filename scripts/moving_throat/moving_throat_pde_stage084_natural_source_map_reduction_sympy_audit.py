#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

mhat0, chiQ, NQ, DeltaQ = sp.symbols('mhat0 chiQ NQ DeltaQ', positive=True, real=True)

# exact factorized condition
sol_NQ = sp.solve(sp.Eq(mhat0**2 * chiQ * NQ, 1), NQ)[0]
print('NQ from exact factorized odd normalization =', sol_NQ)
print('point-particle natural branch mhat0->1 gives =', sp.simplify(sol_NQ.subs(mhat0, 1)))
print('canonical compact outgoing branch chiQ=1 gives =', sp.simplify(sol_NQ.subs({mhat0:1, chiQ:1})))

expr_delta = sp.simplify((1/(1+DeltaQ)) - 1)
series_delta = sp.expand(sp.series(expr_delta, DeltaQ, 0, 3).removeO())
print('NQ - 1 in terms of DeltaQ =', expr_delta)
print('small-DeltaQ series =', series_delta)
print('check exact replacement chiQ=1+DeltaQ:', sp.simplify(sol_NQ.subs({mhat0:1, chiQ:1+DeltaQ}) - 1/(1+DeltaQ)))
