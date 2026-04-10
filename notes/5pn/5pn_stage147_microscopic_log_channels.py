
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage141_150_common import *

banner("STAGE 147 — MICROSCOPIC LOG CHANNELS")

r = sp.symbols('r', positive=True, real=True)
root = root1p(r)
gstar = lower_branch_g(r)
Astar = sp.simplify(gstar + 1 / (4 * root))
Bstar = sp.simplify(1 / (2 * root))

dlnKs, dlnZq, dlnTm, dlnvw0, dlnJs, dlnIsq, dlnLW = sp.symbols(
    'dlnKs dlnZq dlnTm dlnvw0 dlnJs dlnIsq dlnLW', real=True
)
dlna, dlnell, dlncs, dlncsw, dlnrhow = sp.symbols(
    'dlna dlnell dlncs dlncsw dlnrhow', real=True
)

subbanner("1. Exact microscopic meaning of the two channels")
print("delta ln[(g_q K_s)/(g_s lambda)] = delta ln(g/r)")
print("delta ln[(K_s K_q)/lambda^2] = - delta ln r_c")

subbanner("2. Exact parent-action formulas")
channel1_general = dlnKs + dlnZq - dlnTm - dlnvw0 - dlnJs - dlnIsq - sp.Rational(3, 2) * dlnLW
channel2_general = dlnKs + dlnZq + 2 * dlncs - 2 * dlnvw0 - 2 * dlnIsq - 2 * dlnLW
print("channel1_general =", channel1_general)
print("channel2_general =", channel2_general)

subbanner("3. Uniform-overlap + D/N simplification")
# dln Js = 2 dln a + dln ell ; dln Isq = 2 dln a + dln ell + 1/2 dln L_W
channel1_dn = sp.expand(channel1_general.subs({
    dlnJs: 2 * dlna + dlnell,
    dlnIsq: 2 * dlna + dlnell + sp.Rational(1, 2) * dlnLW,
}))
channel2_dn = sp.expand(channel2_general.subs({
    dlnIsq: 2 * dlna + dlnell + sp.Rational(1, 2) * dlnLW,
}))
expect_zero("channel1 D/N form", channel1_dn - (dlnKs + dlnZq - dlnTm - dlnvw0 - 4*dlna - 2*dlnell - 2*dlnLW))
expect_zero("channel2 D/N form", channel2_dn - (dlnKs + dlnZq + 2*dlncs - 2*dlnvw0 - 4*dlna - 2*dlnell - 3*dlnLW))
print("channel1_dn =", channel1_dn)
print("channel2_dn =", channel2_dn)

subbanner("4. Healing-locked shell simplification")
# K_s ∝ a^2 c_sw / rho_w ; ell ∝ c_sw^{-1}
channel1_h = sp.expand((dlnZq + 3*dlncsw - dlnrhow - dlnTm - dlnvw0 - 2*dlna - 2*dlnLW))
channel2_h = sp.expand((dlnZq + 2*dlncs + 3*dlncsw - dlnrhow - 2*dlnvw0 - 2*dlna - 3*dlnLW))
print("channel1_healing =", channel1_h)
print("channel2_healing =", channel2_h)

subbanner("5. Explicit off-family scalar")
delta_perp = sp.simplify(gstar * channel1_h + channel2_h / (4 * root))
target = sp.expand(
    Astar * (dlnZq - dlnrhow)
    + 3 * Astar * dlncsw
    + Bstar * dlncs
    - gstar * dlnTm
    - (gstar + Bstar) * dlnvw0
    - 2 * Astar * dlna
    - (2 * gstar + sp.Rational(3,1)/(4 * root)) * dlnLW
)
expect_zero("delta_perp explicit channel form", delta_perp - target)
print("delta_perp =", target)
print("A_* =", Astar)
print("B_* =", Bstar)
print("A_* numeric =", sp.N(Astar.subs(r, rF1), 16))
print("B_* numeric =", sp.N(Bstar.subs(r, rF1), 16))
