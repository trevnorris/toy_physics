
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.factor(sp.together(sp.expand(expr))))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

"""
5pn_stage101_parent_action_core_extraction.py

Stage 101 audit: parent-action extraction of the concrete core parameters
(K_s, K_q, lambda, g_s, g_q).
"""

banner("STAGE 101 — PARENT-ACTION EXTRACTION OF CORE PARAMETERS")

a, ell, H_w, hbar, mpsi, rho_w = sp.symbols("a ell H_w hbar m_psi rho_w", positive=True, real=True)
c_sw, c_s = sp.symbols("c_s_w c_s", positive=True, real=True)
Z_q, mu0 = sp.symbols("Z_q mu_0", positive=True, real=True)
L_W, qstar, v_w0, T_m = sp.symbols("L_W q_* v_w0 T_m", positive=True, real=True)
z = sp.symbols("z", real=True)

# Canonical thin-wall moments
I_f = sp.Rational(1, 3)
I_g = sp.Rational(4, 15)
print("I_f =", I_f)
print("I_g =", I_g)

# 1) Shell stiffness
K_s = sp.simplify(4 * sp.pi * a**2 * (H_w * ell / 3 + hbar**2 / (15 * mpsi * rho_w * ell)))
print("K_s =")
sp.pprint(K_s)

ell_lock = sp.simplify(hbar / (2 * mpsi * c_sw))
K_s_locked = sp.simplify(K_s.subs(H_w, mpsi * c_sw**2 / rho_w).subs(ell, ell_lock))
K_s_locked_target = sp.simplify(3 * sp.pi * a**2 * hbar**2 / (5 * mpsi * rho_w * ell_lock))
expect_zero("healing-lock K_s identity", K_s_locked - K_s_locked_target)

# 2) Mixed side-channel stiffness
chi = sp.sqrt(2 / L_W) * sp.sin(sp.pi * z / (2 * L_W))
I_half = sp.simplify(sp.integrate(chi**2, (z, 0, L_W)))
I_half_prime = sp.simplify(sp.integrate(sp.diff(chi, z)**2, (z, 0, L_W)))
expect_zero("half-wave norm", I_half - 1)
expect_zero("half-wave derivative norm", I_half_prime - sp.pi**2 / (4 * L_W**2))

K_q = sp.simplify(Z_q / mu0 * sp.pi**2 * c_s**2 / (4 * L_W**2))
print("K_q =")
sp.pprint(K_q)

# 3) Shell/mixed hybridization
J_s = sp.simplify(4 * sp.pi * a**2 * ell * I_f)
I_q = sp.simplify(sp.integrate(chi, (z, 0, L_W)))
I_sq = sp.simplify(J_s * I_q)
expect_zero("J_s formula", J_s - 4 * sp.pi * a**2 * ell / 3)
expect_zero("I_q formula", I_q - 2 * sp.sqrt(2 * L_W) / sp.pi)
expect_zero("I_sq formula", I_sq - 8 * sp.sqrt(2) * a**2 * ell * sp.sqrt(L_W) / 3)

lam = sp.simplify(-qstar * v_w0 * I_sq)
print("lambda =")
sp.pprint(lam)

# 4) Mouth couplings
g_s = sp.simplify(T_m * J_s)
g_q = sp.simplify(Z_q / mu0 * sp.diff(chi, z).subs(z, 0))
expect_zero("g_s formula", g_s - T_m * 4 * sp.pi * a**2 * ell / 3)
expect_zero("g_q formula", g_q - Z_q / mu0 * sp.pi / (sp.sqrt(2) * L_W**sp.Rational(3, 2)))

banner("STAGE 101 FINAL LEDGER")
print("Explicit parent-action core parameters:")
print("  K_s =", K_s)
print("  K_q =", K_q)
print("  lambda =", lam)
print("  g_s =", g_s)
print("  g_q =", g_q)
print()
print("On the healing-locked shell branch:")
print("  ell = hbar/(2 m_psi c_{s,w})")
print("  K_s = 3 pi a^2 hbar^2 / (5 m_psi rho_w ell)")
print()
print("So the reduced Stage-97 two-channel core is no longer abstract; every entry has an explicit overlap formula.")
