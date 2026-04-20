#!/usr/bin/env python3
"""
Stage 100 SymPy audit.

Checks:
1. Positive scale/argument deformations preserve the canonical fingerprint only
   on the harmless beta = 1 branch.
2. Pure Robin loading preserves the canonical branch only trivially at rho_R = 0.
3. A standalone mixed pole cannot preserve the canonical even branch except by
   disappearing.
4. The hybrid outlet has a trivial cancellation branch and a nontrivial
   compensated Robin-mixed branch.
5. The concrete two-channel core balance surface reproduces the nontrivial
   compensated branch.
6. The D/N mixed-tube normalization realizes kappa_c = 1/3 and gamma_c = 1/9,
   while leaving the actual microscopic throat-core realization open.
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


banner("STAGE 100 — CONCRETE OUTLET-CORE STATUS")

z = sp.symbols("z", real=True)
S, beta = sp.symbols("S beta", nonzero=True, real=True)
rho = sp.symbols("rho", real=True)
sigma, kappa, gamma = sp.symbols("sigma kappa gamma", real=True)
Ks, Kq, lam = sp.symbols("K_s K_q lam", positive=True, real=True)
gs, gq = sp.symbols("g_s g_q", real=True)
kappa0, gamma0, a = sp.symbols("kappa0 gamma0 a", positive=True, real=True)
Lw = sp.symbols("L_W", positive=True, real=True)
I = sp.I

Lambda_out = -3 + z**2 / 3 + z**4 / 9 + I * z**5 / 9
Y_can = sp.expand(sp.series((-3) / Lambda_out, z, 0, 6).removeO())

banner("1. Harmless scale/argument class")
Y_arg = sp.expand(sp.series((-3 * S) / (S * Lambda_out.subs(z, beta * z)), z, 0, 6).removeO())
m2_arg = sp.simplify(Y_arg.coeff(z, 2))
m4_arg = sp.simplify(Y_arg.coeff(z, 4))
chi_arg = sp.simplify(sp.im(Y_arg.coeff(z, 5)) / sp.Rational(1, 27))
beta_solutions = sp.solve(
    [sp.Eq(m2_arg, sp.Rational(1, 9)), sp.Eq(m4_arg, sp.Rational(4, 81))],
    [beta],
    dict=True,
)
print("scale/argument solutions =", beta_solutions)
if {sol[beta] for sol in beta_solutions} != {sp.Integer(-1), sp.Integer(1)}:
    raise AssertionError("Expected beta = ±1 from canonical-even matching.")
expect_zero("positive harmless branch has beta = 1 and chi_Q = 1", chi_arg.subs(beta, 1) - 1)

banner("2. Pure Robin class")
Lambda_R = Lambda_out + rho
Y_R = sp.expand(sp.series((-3 + rho) / Lambda_R, z, 0, 6).removeO())
c2_R = sp.simplify(Y_R.coeff(z, 2))
c4_R = sp.simplify(Y_R.coeff(z, 4))
chi_R = sp.simplify(sp.im(Y_R.coeff(z, 5)) / sp.Rational(1, 27))
robin_even_solutions = sp.solve(
    [sp.Eq(c2_R, sp.Rational(1, 9)), sp.Eq(c4_R, sp.Rational(4, 81))],
    [rho],
    dict=True,
)
print("pure Robin canonical-even solutions =", robin_even_solutions)
if robin_even_solutions != [{rho: sp.Integer(0)}]:
    raise AssertionError("Pure Robin should preserve the canonical branch only trivially.")
expect_zero("pure Robin odd norm is trivial on rho = 0", chi_R.subs(rho, 0) - 1)

banner("3. Standalone mixed-pole class")
Lambda_mix = sp.expand(sp.series(Lambda_out - sigma / (1 - kappa * z**2 - I * gamma * z**5), z, 0, 6).removeO())
L0_mix = sp.simplify(Lambda_mix.coeff(z, 0))
L2_mix = sp.simplify(Lambda_mix.coeff(z, 2))
L4_mix = sp.simplify(Lambda_mix.coeff(z, 4))
L5_mix = sp.simplify(sp.im(Lambda_mix.coeff(z, 5)))
kappa_match = sp.solve(sp.Eq(-L2_mix / L0_mix, sp.Rational(1, 9)), kappa)[0]
sigma_match = sp.solve(
    sp.Eq((L2_mix**2 / L0_mix**2 - L4_mix / L0_mix).subs(kappa, kappa_match), sp.Rational(4, 81)),
    sigma,
)[0]
chi_mix = sp.simplify((-L5_mix / L0_mix) / sp.Rational(1, 27))
print("standalone mixed-pole kappa match =", kappa_match)
print("standalone mixed-pole sigma match =", sigma_match)
expect_zero("formal even-match forces kappa = -1/9", kappa_match + sp.Rational(1, 9))
expect_zero("standalone mixed pole disappears on the canonical branch", sigma_match)
expect_zero("odd norm is then trivial", chi_mix.subs(sigma, 0) - 1)

banner("4. Hybrid outlet class split")
Lambda_hyb = sp.expand(sp.series(Lambda_out + rho - sigma / (1 - kappa * z**2 - I * gamma * z**5), z, 0, 6).removeO())
L0_hyb = sp.simplify(Lambda_hyb.coeff(z, 0))
L2_hyb = sp.simplify(Lambda_hyb.coeff(z, 2))
L4_hyb = sp.simplify(Lambda_hyb.coeff(z, 4))
L5_hyb = sp.simplify(sp.im(Lambda_hyb.coeff(z, 5)))
hybrid_solutions = sp.solve(
    [
        sp.Eq(-L2_hyb / L0_hyb, sp.Rational(1, 9)),
        sp.Eq(L2_hyb**2 / L0_hyb**2 - L4_hyb / L0_hyb, sp.Rational(4, 81)),
    ],
    [rho, kappa],
    dict=True,
)
print("hybrid canonical-even branches =", hybrid_solutions)
branch_cancel = next(sol for sol in hybrid_solutions if sol[kappa] == 0)
branch_comp = next(sol for sol in hybrid_solutions if sol[kappa] == sp.Rational(1, 3))

chi_cancel = sp.simplify(((-L5_hyb / L0_hyb) / sp.Rational(1, 27)).subs(branch_cancel))
chi_comp = sp.simplify(((-L5_hyb / L0_hyb) / sp.Rational(1, 27)).subs(branch_comp))
expect_zero("hybrid cancellation branch odd norm", chi_cancel - (1 - 9 * sigma * gamma))
expect_zero(
    "hybrid cancellation branch is trivial when gamma = 0",
    sp.expand(sp.series(Lambda_hyb.subs(branch_cancel).subs(gamma, 0) - Lambda_out, z, 0, 6).removeO()),
)
expect_zero("compensated branch odd norm", chi_comp.subs(gamma, sp.Rational(1, 9)) - 1)
expect_zero(
    "compensated branch collapses to a pure scale deformation",
    sp.expand(
        sp.series(
            Lambda_hyb.subs(branch_comp).subs(gamma, sp.Rational(1, 9)) - (1 - sigma) * Lambda_out,
            z,
            0,
            6,
        ).removeO()
    ),
)

banner("5. Concrete core realization of the compensated class")
r_c = sp.simplify(lam**2 / (Ks * Kq))
rho_c = sp.simplify(gs**2 / Ks)
sigma_c = sp.simplify((Ks * gq - lam * gs) ** 2 / (Ks**2 * Kq * (1 + r_c)))
kappa_c = sp.simplify(kappa0 / (1 + r_c))
gamma_c = sp.simplify(gamma0 / (1 + r_c))
gq_solutions = sp.solve(sp.Eq(rho_c - 4 * sigma_c, 0), gq)
print("core-balance surface branches =", gq_solutions)
sigma_star = sp.simplify(gs**2 / (4 * Ks))
expect_zero(
    "both core-balance branches give the same sigma_*",
    sigma_c.subs(gq, gq_solutions[0]) - sigma_c.subs(gq, gq_solutions[1]),
)
expect_zero("core-balance sigma_* value", sigma_c.subs(gq, gq_solutions[0]) - sigma_star)

Lw_required = sp.solve(sp.Eq(4 * Lw**2 / (sp.pi**2 * a**2), (1 + r_c) / 3), Lw)[0]
expect_zero("D/N tube fixes kappa_c = 1/3", kappa_c.subs(kappa0, 4 * Lw_required**2 / (sp.pi**2 * a**2)) - sp.Rational(1, 3))
expect_zero("bare mixed normalization fixes gamma_c = 1/9", gamma_c.subs(gamma0, (1 + r_c) / 9) - sp.Rational(1, 9))

delta_core = sp.simplify(
    rho_c - sigma_c / (1 - kappa_c * z**2 - I * gamma_c * z**5)
).subs({
    gq: gq_solutions[0],
    kappa0: (1 + r_c) / 3,
    gamma0: (1 + r_c) / 9,
})
delta_core_expected = sp.simplify(4 * sigma_star - sigma_star / (1 - z**2 / 3 - I * z**5 / 9))
expect_zero("concrete core collapses to the compensated hybrid class", sp.expand(sp.series(delta_core - delta_core_expected, z, 0, 6).removeO()))

banner("6. Classification capstone")
classification_rows = [
    ("scale/argument", True, True, False, "harmless beta = 1 pure-scale branch"),
    ("pure Robin", False, False, False, "rho_R = 0 only"),
    ("standalone mixed pole", False, False, False, "sigma_W = 0 only (formal kappa = -1/9)"),
    ("hybrid cancellation", True, True, False, "gamma_W = 0 reduces to exact cancellation"),
    ("compensated Robin-mixed core realization", True, True, True, "balance surface + D/N tube normalization"),
]
for row in classification_rows:
    print(row)

nontrivial_survivors = [name for name, even_ok, odd_ok, nontrivial, _ in classification_rows if even_ok and odd_ok and nontrivial]
if nontrivial_survivors != ["compensated Robin-mixed core realization"]:
    raise AssertionError(f"Unexpected nontrivial survivors: {nontrivial_survivors}")

print("\nOpen microscopic question:")
print("  The explicit low-frequency classification is closed at the reduced-model level,")
print("  but the actual moving-throat core still has to realize the balance surface and")
print("  D/N tube normalization. This script does not assert that realization.")

