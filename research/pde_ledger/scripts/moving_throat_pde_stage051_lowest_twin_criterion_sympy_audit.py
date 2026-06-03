#!/usr/bin/env python3
"""
moving_throat_pde_stage34_lowest_twin_criterion_sympy_audit.py

SymPy audit for Stage 51 of the moving-throat PDE program.

Checks
------
1. Exact tracking-branch product Pi_tr = F_tr * G_tr.
2. Elimination of zeta_req in favor of the continuum product law.
3. Exact endpoint limits of Pi_tr on the stable branch.
4. Exact threshold formulas:
      - Lambda_twin,req,
      - M_mix^(twin,req),
      - Z_W^(twin,req).
5. Exact quadratic root xi_(2x) solving G_tr = 2 M_mix.

Provenance notes
----------------
- `F_tr` / `G_tr` are the Stage 045 tracking-branch transport factors, carried
  through the Stage 050/034 product law without renaming.
- The `Z_W` threshold is checked through the same Stage 047 coherent map, so
  the `M_mix <-> Z_W` bridge stays explicit instead of being silently replayed.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


pi = sp.pi

# Symbols
xi, delta, R = sp.symbols("xi delta R", positive=True, real=True)
Lambda, eps = sp.symbols("Lambda eps", positive=True, real=True)
eps_eta, chi0 = sp.symbols("eps_eta chi0", real=True)
Mmix, ZW = sp.symbols("M_mix Z_W", positive=True, real=True)
Pi = sp.symbols("Pi", positive=True, real=True)

# Tracking-branch functions
Gtr = sp.simplify(9 * xi * (xi + delta) / (9 * delta + (9 + 2 * R**2) * xi))
Ftr = sp.simplify(
    (9 * delta + (9 + 2 * R**2) * xi) ** 2
    * (9 * delta + (9 + 2 * R) * xi) ** 2
    / (81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + (9 + 2 * R**2) * xi**2) ** 2)
)
Pi_tr = sp.simplify(sp.factor(Ftr * Gtr))

banner("STAGE 51 — EXACT TRACKING-BRANCH PRODUCT")
print("G_tr(xi,delta;R) =")
sp.pprint(Gtr)
print("\nF_tr(xi,delta;R) =")
sp.pprint(Ftr)
print("\nPi_tr = F_tr * G_tr =")
sp.pprint(Pi_tr)

Pi_expected = sp.simplify(
    xi * (xi + delta) * (9 * delta + (9 + 2 * R) * xi) ** 2 * (9 * delta + (9 + 2 * R**2) * xi)
    / (9 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + (9 + 2 * R**2) * xi**2) ** 2)
)
expect_zero("Pi_tr - expected closed form", Pi_tr - Pi_expected)

banner("ENDPOINT LIMITS")
Pi0 = sp.simplify(sp.limit(Pi_tr, xi, 0, dir="+"))
print("Pi_tr(xi->0+) =", Pi0)
if Pi0 != 0:
    raise AssertionError("Pi_tr(xi->0+) should be zero.")
Pi1 = sp.limit(Pi_tr, xi, 1, dir="-")
print("Pi_tr(xi->1-) =", Pi1)
if Pi1 is not sp.oo:
    raise AssertionError("Pi_tr(xi->1-) should diverge to +oo.")

banner("ELIMINATION OF zeta_req")
Cmix = sp.simplify(8 * Lambda * (1 - eps) / pi**2)
Sreq = sp.simplify(Pi / Cmix)
zeta_req = sp.simplify((Sreq - 1) / (1 + eps * (Sreq - 2)))

print("C_mix = 8 Lambda (1-eps) / pi^2 =")
sp.pprint(Cmix)
print("\nS_req = Pi / C_mix =")
sp.pprint(Sreq)
print("\nzeta_req(Pi,C_mix,eps) =")
sp.pprint(zeta_req)

expect_zero("zeta_req at Pi = C_mix", zeta_req.subs(Pi, Cmix))
expect_zero("zeta_req at Pi = 2 C_mix minus 1", zeta_req.subs(Pi, 2 * Cmix) - 1)

# Exact equivalence zeta_req <= 1  <=>  Pi <= 2 Cmix on the support-needed branch
# Algebraically, zeta_req - 1 vanishes exactly at Pi = 2 Cmix.
expect_zero(
    "zeta_req - 1 at Pi = 2 C_mix",
    sp.factor((zeta_req - 1).subs(Pi, 2 * Cmix)),
)

# Exact substitution of the physical branch product
zeta_req_branch = sp.simplify(zeta_req.subs(Pi, Pi_tr))
print("\nPhysical-branch zeta_req obtained by Pi -> Pi_tr =")
sp.pprint(zeta_req_branch)

banner("EXACT THRESHOLD SCALES")
Lambda_twin_req = sp.simplify(pi**2 * Pi_tr / (16 * (1 - eps)))
Mmix_twin_req = sp.simplify(Gtr / 2)
ZW_twin_req = sp.simplify(pi**2 * (1 - eps_eta) * (1 - eps) * Gtr / (16 * (1 + chi0) ** 2))

print("Lambda_twin,req =")
sp.pprint(Lambda_twin_req)
print("\nM_mix^(twin,req) =")
sp.pprint(Mmix_twin_req)
print("\nZ_W^(twin,req) =")
sp.pprint(ZW_twin_req)

# Stage 047/030 coherent forward map: Z_W = pi^2 (1-eps_eta)(1-eps) M_mix / [8 (1+chi0)^2].
ZW_from_Mmix = sp.simplify(pi**2 * (1 - eps_eta) * (1 - eps) * Mmix / (8 * (1 + chi0) ** 2))
# Apply the forward map to the M_mix threshold M_mix = G_tr/2 and compare to ZW_twin_req.
ZW_threshold_via_map = sp.simplify(ZW_from_Mmix.subs(Mmix, Gtr / 2))
expect_zero(
    "Z_W^(twin,req) - forward-map(M_mix=G_tr/2)",
    ZW_twin_req - ZW_threshold_via_map,
)

banner("EXACT TWIN-SATURATION DEPTH")
xi_2x = sp.simplify(
    (
        2 * Mmix * (9 + 2 * R**2)
        - 9 * delta
        + sp.sqrt((2 * Mmix * (9 + 2 * R**2) - 9 * delta) ** 2 + 648 * Mmix * delta)
    )
    / 18
)
print("xi_(2x)(M_mix,delta;R) =")
sp.pprint(xi_2x)

expect_zero("G_tr(xi_(2x)) - 2 M_mix", sp.simplify(Gtr.subs(xi, xi_2x) - 2 * Mmix))

banner("STAGE 51 THEOREM LEDGER")
print("Exact lowest-twin sufficiency criterion:")
print("  Pi_tr(xi_req,delta;R_tr) <= 16 Lambda (1-eps) / pi^2")
print()
print("Equivalent exact thresholds:")
print("  Lambda_twin,req = pi^2 Pi_tr / [16(1-eps)]")
print("  M_mix^(twin,req) = G_tr/2")
print("  Z_W^(twin,req) = pi^2 (1-eps_eta)(1-eps) G_tr / [16(1+chi0)^2]")
print()
print("The exact twin-saturation depth at fixed M_mix is the closed root xi_(2x) above.")
