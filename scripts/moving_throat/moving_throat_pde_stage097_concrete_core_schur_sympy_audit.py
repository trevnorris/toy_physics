
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr: sp.Expr) -> None:
    s = sp.simplify(sp.expand(expr))
    print(f"{name} = {s}")
    if s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 97 — CONCRETE TWO-CHANNEL CORE OUTLET MODEL")

K_s, K_q, lam = sp.symbols("K_s K_q lam", positive=True, real=True)
g_s, g_q = sp.symbols("g_s g_q", real=True)
kappa0, gamma0, z = sp.symbols("kappa0 gamma0 z", positive=True, real=True)
D = sp.symbols("D")

# Concrete linear core system in algebraic form:
# [[K_s, lam],[lam, -K_q D]] [s,q]^T = u [g_s,g_q]^T
M = sp.Matrix([[K_s, lam], [lam, -K_q * D]])
c = sp.Matrix([g_s, g_q])

delta_D = sp.apart((c.T * M.inv() * c)[0], D)
print("delta_Lambda(D) =")
sp.pprint(delta_D)

rho_c = g_s**2 / K_s
r_c = lam**2 / (K_s * K_q)
sigma_tilde = (K_s * g_q - lam * g_s)**2 / (K_s**2 * K_q)
sigma_c = sigma_tilde / (1 + r_c)
kappa_c = kappa0 / (1 + r_c)
gamma_c = gamma0 / (1 + r_c)

target_D = rho_c - sigma_tilde / (D + r_c)
expect_zero("Schur form identity", delta_D - target_D)

D_bare = 1 - kappa0 * z**2 - sp.I * gamma0 * z**5
delta_z = sp.simplify(delta_D.subs(D, D_bare))
target_z = sp.simplify(
    rho_c - sigma_c / (1 - kappa_c * z**2 - sp.I * gamma_c * z**5)
)
expect_zero("low-frequency normalized outlet identity", delta_z - target_z)

print("\nExact core-level identifications:")
print("rho_c   =", sp.simplify(rho_c))
print("sigma_c =", sp.simplify(sigma_c))
print("kappa_c =", sp.simplify(kappa_c))
print("gamma_c =", sp.simplify(gamma_c))
