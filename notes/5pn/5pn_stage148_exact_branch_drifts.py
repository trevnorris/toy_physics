
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage141_150_common import *

banner("STAGE 148 — EXACT LOWER-BRANCH DRIFTS")

dlnZq, dlnrhow, dlncsw, dlncs, dlna, dlnLW, dlnvw0, dlnTm = sp.symbols(
    'dlnZq dlnrhow dlncsw dlncs dlna dlnLW dlnvw0 dlnTm', real=True
)
dlnTheta = sp.symbols('dlnTheta', real=True)

subbanner("1. Exact D/N geometric co-transport")
print("delta ln L_W = delta ln a")

subbanner("2. Exact background-flow and traction transport")
eq_r = sp.Eq(dlnZq + 2*dlncs + 3*dlncsw - dlnrhow - 2*dlnvw0 - 2*dlna - 3*dlnLW, 0)
eq_g = sp.Eq(dlnZq + 3*dlncsw - dlnrhow - dlnTm - dlnvw0 - 2*dlna - 2*dlnLW, 0)
sol = sp.solve([eq_r.subs(dlnLW, dlna), eq_g.subs(dlnLW, dlna)], [dlnvw0, dlnTm], dict=True)[0]
print("delta ln v_w0 =", sp.simplify(sol[dlnvw0]))
print("delta ln T_m  =", sp.simplify(sol[dlnTm]))
target_v = sp.Rational(1,2)*(dlnZq - dlnrhow) + sp.Rational(3,2)*dlncsw + dlncs - sp.Rational(5,2)*dlna
target_t = sp.Rational(1,2)*(dlnZq - dlnrhow) + sp.Rational(3,2)*dlncsw - dlncs - sp.Rational(3,2)*dlna
expect_zero("v_w0 transport", sol[dlnvw0] - target_v)
expect_zero("T_m transport", sol[dlnTm] - target_t)

subbanner("3. Product/ratio factorization")
ratio = sp.simplify(sol[dlnvw0] - sol[dlnTm])
product = sp.simplify(sol[dlnvw0] + sol[dlnTm])
expect_zero("ratio drift", ratio - (2*dlncs - dlna))
expect_zero("product drift", product - (dlnZq + 3*dlncsw - dlnrhow - 4*dlna))
print("delta ln(v_w0 / T_m) =", ratio)
print("delta ln(v_w0 T_m)   =", product)

subbanner("4. Frozen n=5 wall-EOS reduction")
# c_sw^2 proportional rho_w^4 => dln c_sw = 2 dln rho_w
v_n5 = sp.simplify(sol[dlnvw0].subs(dlncsw, 2*dlnrhow))
t_n5 = sp.simplify(sol[dlnTm].subs(dlncsw, 2*dlnrhow))
print("delta ln v_w0 (n=5) =", v_n5)
print("delta ln T_m  (n=5) =", t_n5)
expect_zero("n=5 ratio drift", sp.simplify(v_n5 - t_n5) - (2*dlncs - dlna))
expect_zero("n=5 product drift", sp.simplify(v_n5 + t_n5) - (dlnZq + 5*dlnrhow - 4*dlna))

subbanner("5. Collapse of the Stage-147 off-family defect")
r = sp.symbols('r', positive=True, real=True)
root = root1p(r)
gstar = lower_branch_g(r)
Astar = sp.simplify(gstar + 1 / (4 * root))
Bstar = sp.simplify(1 / (2 * root))
delta_perp = sp.simplify(
    Astar * (dlnZq - dlnrhow)
    + 3 * Astar * dlncsw
    + Bstar * dlncs
    - gstar * sol[dlnTm]
    - (gstar + Bstar) * sol[dlnvw0]
    - 2 * Astar * dlna
    - (2 * gstar + sp.Rational(3,1)/(4 * root)) * dlna
)
expect_zero("delta_perp on the exact lower branch", sp.expand(delta_perp.subs(dlnLW, dlna)))
