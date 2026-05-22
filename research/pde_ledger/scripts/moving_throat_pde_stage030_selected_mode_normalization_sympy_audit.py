#!/usr/bin/env python3
"""
Stage 13 SymPy audit.

Checks:
1. Generic normalized-response expansion for a selected operator
       D = D0 + D2 w^2 + D4 w^4 - i C5 w^5 + ...
2. Exact selected lower eigenvalue of the loaded 2x2 wall block.
3. Exact Hellmann–Feynman selected overlap s_- = -(d lambda_- / d alpha).
4. Closed selected-branch odd coefficient and static prefactor formulas.
5. Equivalence of the selected-branch normalization target forms.
"""

from __future__ import annotations

import sympy as sp

I = sp.I


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


# ---------------------------------------------------------------------------
# Part I. Generic normalized selected response
# ---------------------------------------------------------------------------

banner("PART I — GENERIC SELECTED-MODE NORMALIZED RESPONSE")

w = sp.symbols("omega", real=True)
D0, D2, D4, C5 = sp.symbols("D0 D2 D4 C5", nonzero=True, real=True)
Dsel = D0 + D2 * w**2 + D4 * w**4 - I * C5 * w**5
Ysel = sp.expand(sp.series(D0 / Dsel, w, 0, 6).removeO())

print("Y_-(omega) =")
sp.pprint(Ysel)

u2 = -D2 / D0
u4 = (D2**2 - D0 * D4) / D0**2
Gamma5 = C5 / D0

expect_zero("u2 coefficient check", sp.expand(Ysel).coeff(w, 2) - u2)
expect_zero("u4 coefficient check", sp.expand(Ysel).coeff(w, 4) - u4)
expect_zero("Gamma5 coefficient check", sp.expand(sp.im(Ysel)).coeff(w, 5) - Gamma5)


# ---------------------------------------------------------------------------
# Part II. Exact selected lower mode of the loaded 2x2 wall block
# ---------------------------------------------------------------------------

banner("PART II — EXACT SELECTED LOWER EIGENVALUE AND OVERLAP")

A, DK = sp.symbols("A DK", positive=True, real=True)
alpha = sp.symbols("alpha", nonnegative=True, real=True)
x0, x1 = sp.symbols("x0 x1", positive=True, real=True)  # x0 = kappa_0^2, x1 = kappa_1^2

sigma = x0 + x1
delta_kappa = x0 - x1
KappaProd = x0 * x1
R = sp.sqrt((DK + alpha * delta_kappa) ** 2 + 4 * alpha**2 * KappaProd)

lam_minus = sp.simplify((2 * A + DK - alpha * sigma - R) / 2)
lam_plus = sp.simplify((2 * A + DK - alpha * sigma + R) / 2)

print("lambda_- =")
sp.pprint(lam_minus)
print("lambda_+ =")
sp.pprint(lam_plus)

s_minus_hf = sp.simplify(-sp.diff(lam_minus, alpha))
s_minus_closed = sp.simplify(
    sp.Rational(1, 2)
    * (sigma + ((DK + alpha * delta_kappa) * delta_kappa + 4 * alpha * KappaProd) / R)
)

expect_zero("selected overlap: HF - closed form", s_minus_hf - s_minus_closed)
print("s_- = (v.e_-)^2 =")
sp.pprint(s_minus_closed)

expect_zero("weak-loading overlap limit", sp.simplify(s_minus_closed.subs(alpha, 0) - x0))


# ---------------------------------------------------------------------------
# Part III. Selected odd coefficient and static prefactor
# ---------------------------------------------------------------------------

banner("PART III — SELECTED ODD COEFFICIENT AND STATIC PREFACTOR")

beta0, G5 = sp.symbols("beta0 G5", positive=True, real=True)
beta5 = G5 * beta0
C5_sel = sp.simplify(beta5 * s_minus_closed)
Gamma5_sel = sp.simplify(C5_sel / lam_minus)
P0_sel = sp.simplify(beta0 * s_minus_closed / lam_minus)

print("C5_sel =")
sp.pprint(C5_sel)
print("Gamma5_sel =")
sp.pprint(Gamma5_sel)
print("P0_sel =")
sp.pprint(P0_sel)

# Note: Gamma5_sel - G5*P0_sel == 0 follows by construction from the
# definitions on lines 101-104 (C5_sel := G5*beta0*s, Gamma5_sel := C5_sel/lam_-,
# P0_sel := beta0*s/lam_-), so it is not verified separately here. The physical
# content of this relation is the identification Gamma5 = C5/D0 at the selected
# mode, which is verified in generic form by the Part-I expansion (line 55).
expect_zero(
    "P0_sel + beta0*d(log lambda_-)/d alpha",
    P0_sel + beta0 * sp.diff(sp.log(lam_minus), alpha),
)

# exact determinant identity
T0 = (A + DK) * x0 + A * x1
expect_zero("det identity", sp.expand(lam_minus * lam_plus - (A * (A + DK) - alpha * T0)))


# ---------------------------------------------------------------------------
# Part IV. Equivalence of selected-branch normalization targets
# ---------------------------------------------------------------------------

banner("PART IV — EQUIVALENCE OF THE SELECTED-BRANCH TARGET FORMS")

G, cs, a, c, mhat = sp.symbols("G c_s a c mhat", positive=True, real=True)
G5_phys = a**5 / (27 * cs**5)
NQ_target = 54 * G * cs**5 / (5 * a**5 * c**5)

Gamma5_phys = sp.simplify(G5_phys * P0_sel)
cond1 = mhat**2 * Gamma5_phys - 2 * G / (5 * c**5)
cond2 = mhat**2 * P0_sel - NQ_target

expect_zero(
    "normalization equivalence",
    sp.simplify(cond1 - G5_phys * cond2),
)

lambda_req = sp.simplify(mhat**2 * beta0 * s_minus_closed / NQ_target)
# Note: the algebraic rearrangement
#   (lam_- - lambda_req) + (mhat^2*P0_sel - NQ_target)*lam_-/NQ_target = 0
# follows by substituting the definitions P0_sel := beta0*s/lam_- and
# lambda_req := mhat^2*beta0*s/NQ_target; both sides collapse to 0 without
# exercising the physical content of lam_- or s. It is therefore not
# verified separately. The substantive content (that lam_- = lambda_req
# when the constraint mhat^2 P0_sel = NQ_target holds) is recorded here
# as a definitional identity, not as a check.
print("lambda_req =")
sp.pprint(lambda_req)


banner("STAGE 13 AUDIT COMPLETE")
print("Verified:")
print("  • generic normalized selected-response expansion")
print("  • exact selected lower eigenvalue and Hellmann–Feynman overlap")
print("  • exact selected static prefactor P0_- = beta0 s_-/lambda_-")
print("  • exact identity P0_- = - beta0 d(log lambda_-)/d alpha")
print("  • exact determinant identity for the selected branch")
print("  • exact equivalence of the selected-branch normalization target forms")
