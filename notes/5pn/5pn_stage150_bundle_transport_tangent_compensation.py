
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage141_150_common import *

banner("STAGE 150 — BUNDLE TRANSPORT AND TANGENT-COMPENSATION")

dlnTheta, dlnKs, dlnKq, dlnP0 = sp.symbols('dlnTheta dlnKs dlnKq dlnP0', real=True)
dlnrhow, dlncsw, dlnell, dlna, dlnLW, dlncs, dlnZq, dlnTm, dlnvw0 = sp.symbols(
    'dlnrhow dlncsw dlnell dlna dlnLW dlncs dlnZq dlnTm dlnvw0', real=True
)
dlngs, dlngq, dlnlambda = sp.symbols('dlngs dlngq dlnlambda', real=True)
dlnrc, dlnr, dlng = sp.symbols('dlnrc dlnr dlng', real=True)

subbanner("1. Carry-forward inversion formulas")
subs149 = {
    dlnrhow: dlnTheta/2,
    dlna: dlnKs/2 - dlnTheta/4,
    dlncs: dlnKs/2 - dlnTheta/4 + dlnP0/5,
    dlnZq: dlnKq - 2*dlnP0/5,
}
dlncsw_expr = dlnTheta
dlnell_expr = -dlnTheta
dlnLW_expr = dlnKs/2 - dlnTheta/4

print("delta ln c_sw =", dlncsw_expr)
print("delta ln ell =", dlnell_expr)

subbanner("2. Exact bundle transport of T_m and v_w0")
# Stage-148 laws:
dlnvw0_expr = sp.simplify((sp.Rational(1,2)*(dlnZq - dlnrhow) + sp.Rational(3,2)*dlncsw + dlncs - sp.Rational(5,2)*dlna).subs(subs149).subs({dlncsw: dlncsw_expr}))
dlnTm_expr = sp.simplify((sp.Rational(1,2)*(dlnZq - dlnrhow) + sp.Rational(3,2)*dlncsw - dlncs - sp.Rational(3,2)*dlna).subs(subs149).subs({dlncsw: dlncsw_expr}))
expect_zero("v_w0 transport residual", dlnvw0_expr - (-sp.Rational(3,4)*dlnKs + sp.Rational(1,2)*dlnKq + sp.Rational(13,8)*dlnTheta))
expect_zero("T_m transport residual", dlnTm_expr - (-sp.Rational(5,4)*dlnKs + sp.Rational(1,2)*dlnKq + sp.Rational(15,8)*dlnTheta - sp.Rational(2,5)*dlnP0))
print("delta ln v_w0 =", dlnvw0_expr)
print("delta ln T_m  =", dlnTm_expr)

ratio_expr = sp.simplify(dlnvw0_expr - dlnTm_expr)
product_expr = sp.simplify(dlnvw0_expr + dlnTm_expr)
expect_zero("ratio drift", ratio_expr - (sp.Rational(1,2)*dlnKs - sp.Rational(1,4)*dlnTheta + sp.Rational(2,5)*dlnP0))
expect_zero("product drift", product_expr - (-2*dlnKs + dlnKq + sp.Rational(7,2)*dlnTheta - sp.Rational(2,5)*dlnP0))
print("delta ln(v_w0/T_m) =", ratio_expr)
print("delta ln(v_w0 T_m) =", product_expr)

subbanner("3. Exact bundle transport of g_s, g_q, lambda")
dlngs_expr = sp.simplify(dlnTm_expr + 2*subs149[dlna] + dlnell_expr)
dlngq_expr = sp.simplify(subs149[dlnZq] - sp.Rational(3,2)*dlnLW_expr)
dlnlambda_expr = sp.simplify(dlnvw0_expr + 2*subs149[dlna] + dlnell_expr + sp.Rational(1,2)*dlnLW_expr)
expect_zero("g_s transport residual", dlngs_expr - (-sp.Rational(1,4)*dlnKs + sp.Rational(1,2)*dlnKq + sp.Rational(3,8)*dlnTheta - sp.Rational(2,5)*dlnP0))
expect_zero("g_q transport residual", dlngq_expr - (-sp.Rational(3,4)*dlnKs + dlnKq + sp.Rational(3,8)*dlnTheta - sp.Rational(2,5)*dlnP0))
expect_zero("lambda transport residual", dlnlambda_expr - (sp.Rational(1,2)*dlnKs + sp.Rational(1,2)*dlnKq))
print("delta ln g_s =", dlngs_expr)
print("delta ln g_q =", dlngq_expr)
print("delta ln lambda =", dlnlambda_expr)

subbanner("4. Tangent-compensation theorem")
dlnrc_expr = sp.simplify(2*dlnlambda_expr - dlnKs - dlnKq)
dlnr_expr = sp.simplify(dlnlambda_expr - (dlnKs + dlnKq)/2)
dlng_expr = sp.simplify(dlngq_expr + dlnKs/2 - dlngs_expr - dlnKq/2)
expect_zero("delta ln r_c", dlnrc_expr)
expect_zero("delta ln mathfrak r", dlnr_expr)
expect_zero("delta ln mathfrak g", dlng_expr)

channel1 = sp.simplify(dlngq_expr + dlnKs - dlngs_expr - dlnlambda_expr)
channel2 = sp.simplify(dlnKs + dlnKq - 2*dlnlambda_expr)
expect_zero("channel1", channel1)
expect_zero("channel2", channel2)

r = sp.symbols('r', positive=True, real=True)
gstar = lower_branch_g(r)
delta_perp = sp.simplify(gstar*channel1 + channel2/(4*root1p(r)))
expect_zero("delta_perp", delta_perp)
print("Arbitrary first-order isotropic bundle drift is tangent to the exact parent compensation family.")
