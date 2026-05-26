#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 41 SymPy audit.

Checks:
1. exact support-drop kernel and its positive derivative identity,
2. normalization of the constructive source family Sigma_Pe,
3. exact antiderivative formulas for Ic and Is,
4. exact uniform-source drop Delta_0,
5. exact sharp-bottom endpoint Delta_inf,
6. fixed-point branch bracket Pe_* in [Xi Delta_0, Xi Delta_inf].
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


banner("STAGE 41 — COUPLED SUPPORT/SOURCE OPERATOR")

x, s = sp.symbols("x s", real=True)
Pe, alpha, eta, Xi = sp.symbols("Pe alpha eta Xi", positive=True, real=True)

W = sp.simplify(alpha * sp.sinh(alpha) + eta * sp.cosh(alpha))
K = sp.simplify(
    (sp.cosh(alpha * x) + (eta / alpha) * sp.sinh(alpha * x) - sp.cosh(alpha * (1 - x))) / W
)
Kprime = sp.simplify(sp.diff(K, x))
print("K_(alpha,eta)(x) =", K)
print("dK/dx =", Kprime)
expect_zero(
    "Kprime identity",
    Kprime
    - (alpha * sp.sinh(alpha * x) + eta * sp.cosh(alpha * x) + alpha * sp.sinh(alpha * (1 - x))) / W,
)

# --- dK/dx > 0 numerator positivity sweep (notes section 4) ---
kprime_num = alpha * sp.sinh(alpha * x) + eta * sp.cosh(alpha * x) + alpha * sp.sinh(alpha * (1 - x))
for a_val in [sp.Rational(1, 10), sp.Integer(1), sp.Integer(3)]:
    for e_val in [sp.Rational(1, 10), sp.Integer(1), sp.Integer(10)]:
        for x_val in [sp.Rational(0), sp.Rational(1, 4), sp.Rational(1, 2), sp.Rational(3, 4), sp.Rational(1)]:
            val = float(sp.N(kprime_num.subs({alpha: a_val, eta: e_val, x: x_val})))
            if val <= 0:
                raise AssertionError(
                    f"kernel numerator non-positive at alpha={a_val}, eta={e_val}, x={x_val}: {val}"
                )
print("kernel numerator positivity sweep = PASS")

Sigma = sp.simplify(Pe * sp.exp(Pe * x) / (sp.exp(Pe) - 1))
print("Sigma_Pe(x) =", Sigma)
expect_zero("Sigma normalization", sp.integrate(Sigma, (x, 0, 1)) - 1)

# Exact auxiliary integrals for the closed-form support drop.
Fc = sp.exp(Pe * x) * (Pe * sp.cosh(alpha * x) - alpha * sp.sinh(alpha * x)) / (Pe**2 - alpha**2)
Fs = sp.exp(Pe * x) * (Pe * sp.sinh(alpha * x) - alpha * sp.cosh(alpha * x)) / (Pe**2 - alpha**2)
expect_zero("Ic antiderivative check", sp.diff(Fc, x) - sp.exp(Pe * x) * sp.cosh(alpha * x))
expect_zero("Is antiderivative check", sp.diff(Fs, x) - sp.exp(Pe * x) * sp.sinh(alpha * x))

Ic = sp.simplify(Fc.subs(x, 1) - Fc.subs(x, 0))
Is = sp.simplify(Fs.subs(x, 1) - Fs.subs(x, 0))
print("Ic(Pe,alpha) =", Ic)
print("Is(Pe,alpha) =", Is)

Delta = sp.simplify(
    Pe / (sp.exp(Pe) - 1) * ((1 - sp.cosh(alpha)) * Ic + (eta / alpha + sp.sinh(alpha)) * Is) / W
)
print("Delta(Pe;alpha,eta) =", Delta)

Delta0 = sp.simplify(sp.limit(Delta, Pe, 0))
Delta0_expected = sp.simplify(eta * (sp.cosh(alpha) - 1) / (alpha**2 * W))
print("Delta_0 =", Delta0)
expect_zero("Delta0 formula", Delta0 - Delta0_expected)
expect_zero("Delta0 integral identity", Delta0 - sp.integrate(K, (x, 0, 1)))

Delta_inf = sp.simplify(K.subs(x, 1))
Delta_inf_expected = sp.simplify((sp.cosh(alpha) + (eta / alpha) * sp.sinh(alpha) - 1) / W)
print("Delta_inf =", Delta_inf)
expect_zero("Delta_inf direct substitution (sanity)", Delta_inf - Delta_inf_expected)

# --- BVP independence check via numerical kernel integral ---
# Notes section 3: the support drop equals integral_0^1 K(x) Sigma_Pe(x) dx, where K is
# the support-drop kernel from the Green-function representation of -Phi'' + alpha^2 Phi
# = Sigma_Pe with Robin (x=0) / Neumann (x=1) BCs. The symbolic version of this integral
# (sp.integrate over [0,1] with symbolic alpha) does not terminate in sympy; we instead
# sample concrete (alpha, eta, Pe) values, integrate K*Sigma_Pe numerically, and compare
# to the closed-form Delta. This catches sign/BC errors in the kernel ansatz that would
# survive the existing Delta_0 = integral(K) check (which only exercises the Pe -> 0
# limit).
kernel_test_combos = [
    (sp.Rational(1, 2), sp.Integer(1), sp.Integer(1)),
    (sp.Integer(1), sp.Integer(1), sp.Integer(2)),
    (sp.Integer(3), sp.Rational(1, 10), sp.Rational(1, 2)),
    (sp.Integer(1), sp.Integer(10), sp.Integer(5)),
]
for a_val, e_val, p_val in kernel_test_combos:
    integrand_num = (K * Sigma).subs({alpha: a_val, eta: e_val, Pe: p_val})
    kernel_int = float(sp.N(sp.integrate(integrand_num, (x, 0, 1))))
    delta_val = float(sp.N(Delta.subs({alpha: a_val, eta: e_val, Pe: p_val})))
    if abs(kernel_int - delta_val) > 1e-10:
        raise AssertionError(
            f"kernel integral {kernel_int} != Delta {delta_val} at "
            f"alpha={a_val}, eta={e_val}, Pe={p_val}"
        )
print("Delta = integral(K * Sigma_Pe) numerical sweep = PASS")

# Fixed-point branch bracket.
Pe_lo = sp.simplify(Xi * Delta0)
Pe_hi = sp.simplify(Xi * Delta_inf)
print("Pe_lo = Xi Delta_0 =", Pe_lo)
print("Pe_hi = Xi Delta_inf =", Pe_hi)

# Bracket non-emptiness: Delta_inf >= Delta_0 for all alpha, eta > 0.
bracket_gap = sp.simplify(sp.together(Delta_inf - Delta0))
bracket_gap_expected = sp.simplify(
    ((alpha**2 - eta) * (sp.cosh(alpha) - 1) + alpha * eta * sp.sinh(alpha))
    / (alpha**2 * W)
)
expect_zero("bracket gap closed form", bracket_gap - bracket_gap_expected)
for a_val in [sp.Rational(1, 10), sp.Integer(1), sp.Integer(3)]:
    for e_val in [sp.Rational(1, 10), sp.Integer(1), sp.Integer(10)]:
        val = float(sp.N(bracket_gap.subs({alpha: a_val, eta: e_val})))
        if val <= 0:
            raise AssertionError(
                f"bracket gap non-positive at alpha={a_val}, eta={e_val}: {val}"
            )
print("bracket gap positivity sweep = PASS")

# --- Delta(Pe; alpha, eta) monotonicity sweep on the constructive branch ---
# Paper claim: Delta_0 <= Delta(Pe) <= Delta_inf for Pe >= 0, alpha, eta > 0.
sample_alpha = [sp.Rational(1, 10), sp.Integer(1), sp.Integer(3)]
sample_eta = [sp.Rational(1, 10), sp.Integer(1), sp.Integer(10)]
sample_Pe = [sp.Rational(1, 2), sp.Integer(1), sp.Integer(3), sp.Integer(10)]
for a_val in sample_alpha:
    for e_val in sample_eta:
        d0 = float(sp.N(Delta0.subs({alpha: a_val, eta: e_val})))
        dinf = float(sp.N(Delta_inf.subs({alpha: a_val, eta: e_val})))
        for p_val in sample_Pe:
            # Delta has a removable singularity at Pe == alpha (Pe^2 - alpha^2 = 0
            # in denominator with numerator -> 0); subs() doesn't take the limit.
            if p_val == a_val:
                continue
            d_val = float(sp.N(Delta.subs({alpha: a_val, eta: e_val, Pe: p_val})))
            if d_val < d0 - 1e-9:
                raise AssertionError(
                    f"Delta(Pe={p_val}) < Delta_0 at alpha={a_val}, eta={e_val}: "
                    f"Delta={d_val}, Delta_0={d0}"
                )
            if d_val > dinf + 1e-9:
                raise AssertionError(
                    f"Delta(Pe={p_val}) > Delta_inf at alpha={a_val}, eta={e_val}: "
                    f"Delta={d_val}, Delta_inf={dinf}"
                )
print("Delta(Pe) monotonicity sweep = PASS")

# --- IVT bracket-existence check: F(Xi*Delta_0) <= 0 and F(Xi*Delta_inf) >= 0 ---
# Paper claim (notes section 5): F(Pe) := Pe - Xi*Delta(Pe; alpha, eta) satisfies
# F(Xi*Delta_0) <= 0 and F(Xi*Delta_inf) >= 0, so a constructive root exists.
sample_Xi = [sp.Rational(1, 2), sp.Integer(1), sp.Integer(2)]
for a_val in sample_alpha:
    for e_val in sample_eta:
        d0 = float(sp.N(Delta0.subs({alpha: a_val, eta: e_val})))
        dinf = float(sp.N(Delta_inf.subs({alpha: a_val, eta: e_val})))
        for xi_val in sample_Xi:
            pe_lo_val = float(xi_val) * d0
            pe_hi_val = float(xi_val) * dinf
            # Skip if the sweep happens to land within numerical noise of Pe = alpha
            # (Delta has a removable 0/0 there; substitution doesn't take the limit).
            if abs(pe_lo_val - float(a_val)) < 1e-9 or abs(pe_hi_val - float(a_val)) < 1e-9:
                continue
            d_at_lo = float(sp.N(Delta.subs({alpha: a_val, eta: e_val, Pe: sp.nsimplify(pe_lo_val, rational=False)})))
            d_at_hi = float(sp.N(Delta.subs({alpha: a_val, eta: e_val, Pe: sp.nsimplify(pe_hi_val, rational=False)})))
            F_lo = pe_lo_val - float(xi_val) * d_at_lo
            F_hi = pe_hi_val - float(xi_val) * d_at_hi
            if F_lo > 1e-9:
                raise AssertionError(
                    f"F(Xi*Delta_0) > 0 at alpha={a_val}, eta={e_val}, Xi={xi_val}: F_lo={F_lo}"
                )
            if F_hi < -1e-9:
                raise AssertionError(
                    f"F(Xi*Delta_inf) < 0 at alpha={a_val}, eta={e_val}, Xi={xi_val}: F_hi={F_hi}"
                )
print("F-sign IVT bracket existence sweep = PASS")

# Delta_inf is the sharp-bottom (Pe -> oo) endpoint of Delta(Pe).
Delta_inf_limit = sp.simplify(sp.limit(Delta, Pe, sp.oo))
print("Delta(Pe -> oo) =", Delta_inf_limit)
expect_zero("Delta_inf as Pe -> oo limit", sp.simplify((Delta_inf_limit - Delta_inf_expected).rewrite(sp.exp)))

# Weak-coupling branch law.
Delta_series = sp.series(Delta, Pe, 0, 2).removeO()
print("Delta(Pe) small-Pe series =", Delta_series)
expect_zero("weak-coupling constant term", sp.expand(Delta_series).coeff(Pe, 0) - Delta0)
Pe1_coeff = sp.simplify(sp.expand(Delta_series).coeff(Pe, 1))
print("Delta(Pe) first-order coefficient =", Pe1_coeff)
Pe1_val = sp.N(Pe1_coeff.subs({alpha: sp.Integer(1), eta: sp.Integer(1)}))
if sp.Abs(Pe1_val) < sp.Rational(1, 10**8):
    raise AssertionError(
        f"weak-coupling first-order coefficient vanishes at alpha=eta=1: {Pe1_val}"
    )
print("weak-coupling first-order coefficient nonvanishing at alpha=eta=1: PASS")

# --- Weak-coupling branch law: Pe_*(Xi) = Xi*Delta_0 + O(Xi^2). ---
# By the implicit function theorem applied to F(Pe, Xi) := Pe - Xi*Delta(Pe; alpha, eta) = 0
# at (Pe, Xi) = (0, 0), we have dPe_*/dXi|_{Xi=0} = Delta(0; alpha, eta) = Delta_0.
F = Pe - Xi * Delta
dF_dPe = sp.diff(F, Pe)
dF_dXi = sp.diff(F, Xi)
# Evaluate at Pe = 0, Xi = 0. Take limits because Delta has a 0/0 form at Pe = 0.
dF_dPe_at_origin = sp.limit(dF_dPe.subs(Xi, 0), Pe, 0)
dF_dXi_at_origin = sp.limit(dF_dXi.subs(Xi, 0), Pe, 0)
# At Xi=0, F = Pe, so dF/dPe = 1 and dF/dXi = -Delta(0) = -Delta_0.
dPe_star_dXi_at_zero = sp.simplify(-dF_dXi_at_origin / dF_dPe_at_origin)
expect_zero("weak-coupling branch slope = Delta_0", dPe_star_dXi_at_zero - Delta0_expected)

print("\nStage 41 audit passed.")
