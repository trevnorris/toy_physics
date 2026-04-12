
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage141_150_common import *

banner("STAGE 149 — BUNDLE INVERSION OF THE LAST FOUR DRIFTS")

dlnTheta, dlnKs, dlnKq, dlnP0 = sp.symbols('dlnTheta dlnKs dlnKq dlnP0', real=True)
dlnrhow, dlna, dlncs, dlnZq = sp.symbols('dlnrhow dlna dlncs dlnZq', real=True)
dlnN0, dlnD0 = sp.symbols('dlnN0 dlnD0', real=True)

subbanner("1. The four exact branch laws")
eqs = [
    sp.Eq(dlnTheta, 2*dlnrhow),
    sp.Eq(dlnKs, 2*dlna + dlnrhow),
    sp.Eq(dlnKq, dlnZq + 2*dlncs - 2*dlna),
    sp.Eq(dlnP0, 5*(dlncs - dlna)),
]
sol = sp.solve(eqs, [dlnrhow, dlna, dlncs, dlnZq], dict=True)[0]
print(sol)

expect_zero("delta ln rho_w", sol[dlnrhow] - dlnTheta/2)
expect_zero("delta ln a", sol[dlna] - (dlnKs/2 - dlnTheta/4))
expect_zero("delta ln c_s", sol[dlncs] - (dlnKs/2 - dlnTheta/4 + dlnP0/5))
expect_zero("delta ln Z_q", sol[dlnZq] - (dlnKq - 2*dlnP0/5))

subbanner("2. Full-bundle form using P0 = N0 / D0")
cs_bundle = sp.simplify(sol[dlncs].subs(dlnP0, dlnN0 - dlnD0))
zq_bundle = sp.simplify(sol[dlnZq].subs(dlnP0, dlnN0 - dlnD0))
print("delta ln c_s =", cs_bundle)
print("delta ln Z_q =", zq_bundle)

subbanner("3. Frozen-wall corollary")
Theta_chi = sp.Float('4.06863235008162')
lam_mu = sp.symbols('lambda_mu', positive=True, real=True)
rho_chi = sp.sqrt(Theta_chi / 25) / lam_mu
print("rho_w^(chi) =", sp.N(rho_chi.subs(lam_mu, 1), 16), "* lambda_mu^{-1}")
