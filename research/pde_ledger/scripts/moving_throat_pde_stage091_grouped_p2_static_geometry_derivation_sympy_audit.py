#!/usr/bin/env python3
"""
Stage 091 SymPy audit.

Derive the minimal isotropic conservative quadrupole module directly from a
static geometry contact term plus one isotropic grouped-P2 conservative pole.

Carry-forward annotations (paper card `\\stagefield{Checks}` items resolved
upstream, per [[batch-IV1-paper-alignment]] Cluster A direction (a)):

- "Check `l=0` and `l=2` orthogonality before applying the geometry firewall":
  exercised at stage 094 (Isotropic Geometry-Decoupling Theorem) via 15 angular
  integrals plus the Laplace eigenvalue check in both engines. This stage uses
  the scalar `Kcons(omega)` ansatz that the upstream orthogonality result
  authorises.
"""

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


omega = sp.symbols("omega", real=True)
Kgeom, Kpole, OmegaQ = sp.symbols("K_geom K_pole Omega_Q", positive=True, real=True)

banner("STAGE 091 — GROUPED-P2 + STATIC-GEOMETRY DERIVATION")

Kcons = sp.simplify(Kgeom + Kpole / (1 - omega**2 / OmegaQ**2))
series = sp.expand(sp.series(Kcons, omega, 0, 6).removeO())
K0 = sp.simplify(series.coeff(omega, 0))
K2 = sp.simplify(series.coeff(omega, 2))
K4 = sp.simplify(series.coeff(omega, 4))

print("K_Q^cons(omega) =", Kcons)
print("Series =", series)
print("K0 =", K0)
print("K2 =", K2)
print("K4 =", K4)

expect_zero("K0 - (Kgeom + Kpole)", K0 - (Kgeom + Kpole))
expect_zero("K2 - Kpole/OmegaQ^2", K2 - Kpole / OmegaQ**2)
expect_zero("K4 - Kpole/OmegaQ^4", K4 - Kpole / OmegaQ**4)

branch_identity = sp.simplify(K0 * K4 - 4 * K2**2)
print("Branch identity K0*K4 - 4*K2^2 =", branch_identity)

Kgeom_sol = sp.solve(sp.Eq(branch_identity, 0), Kgeom)[0]
print("K_geom forced by branch identity =", Kgeom_sol)
expect_zero("K_geom - 3*K_pole", Kgeom_sol - 3 * Kpole)

K0_on_branch = sp.simplify(K0.subs(Kgeom, Kgeom_sol))
print("K0 on branch =", K0_on_branch)
expect_zero("K0 - 4*K_pole on branch", K0_on_branch - 4 * Kpole)

Yhat = sp.simplify(Kcons.subs(Kgeom, Kgeom_sol) / K0_on_branch)
Yhat_expected = sp.simplify(sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - omega**2 / OmegaQ**2))
print("Normalized module on branch =", Yhat)
expect_zero("Yhat - [3/4 + 1/4/(1-omega^2/OmegaQ^2)]", Yhat - Yhat_expected)

rho_alpha = sp.simplify(K0_on_branch / Kgeom_sol)
zeta_req = sp.simplify((K0_on_branch - Kgeom_sol) / Kgeom_sol)
print("rho_alpha = alpha_req/alpha_mix =", rho_alpha)
print("zeta_req  = (alpha_req-alpha_mix)/alpha_mix =", zeta_req)
expect_zero("rho_alpha - 4/3", rho_alpha - sp.Rational(4, 3))
expect_zero("zeta_req - 1/3", zeta_req - sp.Rational(1, 3))

banner("FINAL LEDGER")
print("If geometry contributes only the static completion and grouped-P2 supplies one isotropic pole,")
print("then the minimal branch identity forces")
print("  K_pole = K0/4,")
print("  K_geom = 3K0/4,")
print("  Yhat_Q^cons = 3/4 + (1/4)/(1-omega^2/OmegaQ^2),")
print("  rho_alpha = 4/3,")
print("  zeta_req  = 1/3.")
