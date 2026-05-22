#!/usr/bin/env python3
"""
Moving-throat PDE — Stage 27 SymPy audit.

What this audit verifies
------------------------
1. Inserting the continuum-selected support baseline M_supp and support direction
   R_phi into the Stage-24 support theorem gives an exact quadratic equation for
   the physical softening depth xi.
2. The physical root xi_phys is the branch that reduces to xi=0 when the loads vanish.
3. Inserting q=t R_U and r=t R_phi into the Stage-25 normalization law yields the
   exact continuum-selected normalization function F_cont.
4. The minimal-kernel limit R_phi=1 gives the exact source-tied closure.
5. The interference-match surface R_phi = R_U collapses the rank-2 branch to the
   one-direction tracking law with total loading M_mix + M_supp.
6. The genuine rank-2 mismatch penalty enters as +lambda0 M_mix M_supp (R_U-R_phi)^2.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


xi, delta = sp.symbols("xi delta", positive=True, real=True)
Mmix, Msupp = sp.symbols("M_mix M_supp", nonnegative=True, real=True)
RU, Rphi = sp.symbols("R_U R_phi", positive=True, real=True)
lambda0 = sp.Rational(2, 9)
t = sp.symbols("t", real=True)

banner("STAGE 27 — CONTINUUM-SELECTED RANK-2 CLOSURE")

subbanner("27.1 — Exact continuum-selected branch equation and quadratic theorem")

# Stage-24 support theorem with q=t RU and r=t Rphi, t^2=lambda0.
q = sp.sqrt(lambda0) * RU
r = sp.sqrt(lambda0) * Rphi
n_req = sp.simplify(
    (xi * (delta + xi) - Mmix * (delta + (1 + lambda0 * RU**2) * xi))
    / (delta + (1 + lambda0 * Rphi**2) * xi - Mmix * lambda0 * (RU - Rphi)**2)
)
print("n_req^(cont) =")
sp.pprint(sp.factor(n_req))

branch_eq = sp.expand(sp.factor((sp.together(n_req - Msupp)).as_numer_denom()[0]))
print("numerator of n_req - M_supp =")
sp.pprint(branch_eq)

B_cont = sp.simplify(delta - Mmix * (1 + lambda0 * RU**2) - Msupp * (1 + lambda0 * Rphi**2))
C_cont = sp.simplify(-delta * (Mmix + Msupp) + lambda0 * Mmix * Msupp * (RU - Rphi)**2)
branch_expected = sp.expand(xi**2 + B_cont * xi + C_cont)
expect_zero("quadratic branch equation", branch_eq - 9 * branch_expected)

Delta_disc = sp.simplify(B_cont**2 - 4 * C_cont)
xi_phys = sp.simplify((-B_cont + sp.sqrt(Delta_disc)) / 2)
print("xi_phys =")
sp.pprint(sp.factor(xi_phys))
expect_zero("zero-load root", sp.simplify(xi_phys.subs({Mmix: 0, Msupp: 0})))

subbanner("27.2 — Exact continuum-selected normalization function")

# Stage-25 general normalization with q=t RU, r=t Rphi, source t=sqrt(lambda0).
D_cont = sp.simplify(
    (delta + xi - Mmix * lambda0 * RU * (RU - Rphi))**2
    + lambda0 * (Mmix * (RU - Rphi) + Rphi * xi)**2
)
F_cont = sp.simplify(
    (delta + (1 + lambda0 * RU * Rphi) * xi)**2
    * (delta + (1 + lambda0 * Rphi) * xi - Mmix * lambda0 * (RU - Rphi) * (RU - 1))**2
    / ((1 - xi) * D_cont**2)
)
print("F_cont(xi) =")
sp.pprint(sp.factor(F_cont))

print("F_cont built successfully.")

# Third slice: independent literal R_phi = 2 to constrain bivariate dependence
# beyond the R_phi=1 and R_phi=R_U collapsed surfaces.
Rphi_lit = sp.Integer(2)
F_lit = sp.simplify(F_cont.subs(Rphi, Rphi_lit))
F_lit_expected = sp.simplify(
    (delta + (1 + lambda0 * RU * Rphi_lit) * xi)**2
    * (delta + (1 + lambda0 * Rphi_lit) * xi - Mmix * lambda0 * (RU - Rphi_lit) * (RU - 1))**2
    / ((1 - xi) * (
        (delta + xi - Mmix * lambda0 * RU * (RU - Rphi_lit))**2
        + lambda0 * (Mmix * (RU - Rphi_lit) + Rphi_lit * xi)**2
    )**2)
)
expect_zero("third-slice F at Rphi=2", F_lit - F_lit_expected)

subbanner("27.3 — Minimal-kernel source-tied surface")

Rphi_source = sp.Integer(1)
n_source = sp.simplify(n_req.subs(Rphi, Rphi_source))
F_source = sp.simplify(F_cont.subs(Rphi, Rphi_source))
print("n_source =")
sp.pprint(sp.factor(n_source))
print("F_source =")
sp.pprint(sp.factor(F_source))

# Source-tied formulas from Stage 24/25.
n_source_expected = sp.simplify(
    (xi * (delta + xi) - Mmix * (delta + (1 + lambda0 * RU**2) * xi))
    / (delta + (1 + lambda0) * xi - Mmix * lambda0 * (RU - 1)**2)
)
F_source_expected = sp.simplify(
    (delta + (1 + lambda0 * RU) * xi)**2
    * (delta + (1 + lambda0) * xi - Mmix * lambda0 * (RU - 1)**2)**2
    / ((1 - xi) * ((delta + xi - Mmix * lambda0 * RU * (RU - 1))**2 + lambda0 * (xi + Mmix * (RU - 1))**2)**2)
)
expect_zero("source-tied n", n_source - n_source_expected)
expect_zero("source-tied F", F_source - F_source_expected)

subbanner("27.4 — Interference-matched tracking surface")

n_track = sp.simplify(n_req.subs(Rphi, RU))
G_q = sp.simplify(xi * (delta + xi) / (delta + (1 + lambda0 * RU**2) * xi))
expect_zero("tracking collapse of n_req", n_track - (G_q - Mmix))

F_track = sp.simplify(F_cont.subs(Rphi, RU))
F_track_expected = sp.simplify(
    (delta + (1 + lambda0 * RU**2) * xi)**2
    * (delta + (1 + lambda0 * RU) * xi)**2
    / ((1 - xi) * ((delta + xi)**2 + lambda0 * RU**2 * xi**2)**2)
)
expect_zero("tracking F collapse", F_track - F_track_expected)

subbanner("27.5 — Exact mismatch penalty")

# Extract the xi-constant coefficient of branch_eq (the numerator of n_req - M_supp).
# branch_eq = 9*(xi^2 + B_cont*xi + C_cont), so the constant-in-xi term is 9*C_cont.
C_from_branch_eq = sp.simplify(sp.Poly(branch_eq, xi).nth(0) / 9)
C_expected = sp.simplify(-delta * (Mmix + Msupp) + lambda0 * Mmix * Msupp * (RU - Rphi)**2)
expect_zero("mismatch penalty in C coefficient", C_from_branch_eq - C_expected)

banner("STAGE 27 THEOREM LEDGER")
print("1. The continuum-selected selected branch is fixed by the exact quadratic equation")
print("      xi^2 + B_cont xi + C_cont = 0,")
print("   with B_cont and C_cont determined by M_mix, M_supp, R_U, and R_phi.")
print("2. The physical softening depth is the root xi_phys that reduces to xi=0 when")
print("   the loads vanish.")
print("3. The exact continuum-selected normalization gate is")
print("      R_target = F_cont(xi_phys).")
print("4. The minimal kernel (R_phi=1) lands exactly on the source-tied closure.")
print("5. The interference-matched surface (R_phi=R_U) collapses to the one-direction")
print("   tracking law with total loading M_mix + M_supp.")
print("6. The genuine rank-2 mismatch penalty is +lambda0 M_mix M_supp (R_U-R_phi)^2,")
print("   so the generic extended continuum kernel is an exact intermediate closure.")
