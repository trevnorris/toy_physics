#!/usr/bin/env python3
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


def expect_zero(name: str, expr) -> None:
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(sp.expand(expr))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


banner("STAGE 205 — DIRECTIONAL HESSIAN AND QUADRATIC ROOT REFINEMENT")

# ---------------------------------------------------------------------------
# I. Exact directional Hessian identities in free log space
# ---------------------------------------------------------------------------
subbanner("I. Exact directional Hessian identities")
Phi0 = sp.symbols("Phi0", positive=True, real=True)
s1, s2, s3, s4, s5 = sp.symbols("s1 s2 s3 s4 s5", real=True)
g1, g2, g3, g4, g5 = sp.symbols("g1 g2 g3 g4 g5", real=True)

H11, H12, H13, H14, H15 = sp.symbols("H11 H12 H13 H14 H15", real=True)
H22, H23, H24, H25 = sp.symbols("H22 H23 H24 H25", real=True)
H33, H34, H35 = sp.symbols("H33 H34 H35", real=True)
H44, H45 = sp.symbols("H44 H45", real=True)
H55 = sp.symbols("H55", real=True)

svec = sp.Matrix([s1, s2, s3, s4, s5])
gvec = sp.Matrix([g1, g2, g3, g4, g5])
Hchi = sp.Matrix(
    [
        [H11, H12, H13, H14, H15],
        [H12, H22, H23, H24, H25],
        [H13, H23, H33, H34, H35],
        [H14, H24, H34, H44, H45],
        [H15, H25, H35, H45, H55],
    ]
)

Phi1 = sp.simplify((svec.T * gvec)[0])
Phi2 = sp.simplify((svec.T * Hchi * svec)[0])
L0 = sp.simplify(Phi1 / Phi0)
Hlog = sp.simplify(Hchi / Phi0 - (gvec * gvec.T) / Phi0**2)
L1 = sp.simplify((svec.T * Hlog * svec)[0])

print("Phi1 =")
sp.pprint(Phi1)
print("Phi2 =")
sp.pprint(Phi2)
print("L0 =")
sp.pprint(L0)
print("L1 =")
sp.pprint(L1)
expect_zero("L1 - (Phi2/Phi0 - Phi1^2/Phi0^2)", L1 - (Phi2 / Phi0 - Phi1**2 / Phi0**2))
expect_zero("Phi2 - Phi0*(L1 + L0^2)", Phi2 - Phi0 * (L1 + L0**2))

# ---------------------------------------------------------------------------
# II. Exact quadratic affine predictor (positive-slope continuation branch)
# ---------------------------------------------------------------------------
subbanner("II. Exact quadratic affine predictor")
Phi0a, Phi1a, Phi2a = sp.symbols("Phi0a Phi1a Phi2a", positive=True, real=True)
tau = sp.symbols("tau", real=True)
Delta_aff = sp.simplify(Phi1a**2 - 2 * Phi2a * (Phi0a - 1))
tau_aff = sp.simplify((1 - Phi0a) / Phi1a)
tau_quad = sp.simplify(2 * (1 - Phi0a) / (Phi1a + sp.sqrt(Delta_aff)))
res_aff = sp.simplify((Phi0a - 1) + Phi1a * tau_quad + sp.Rational(1, 2) * Phi2a * tau_quad**2)

print("Delta_aff =")
sp.pprint(Delta_aff)
print("tau_aff =")
sp.pprint(tau_aff)
print("tau_quad =")
sp.pprint(tau_quad)
expect_zero("quadratic affine residual", res_aff)
expect_zero("limit Phi2 -> 0 gives tau_aff", sp.simplify(sp.limit(tau_quad, Phi2a, 0) - tau_aff))

Phi0a_n, Phi1a_abs, Phi2a_n = sp.symbols("Phi0a_n Phi1a_abs Phi2a_n", positive=True, real=True)
Phi1a_n = -Phi1a_abs
Delta_aff_n = sp.simplify(Phi1a_n**2 - 2 * Phi2a_n * (Phi0a_n - 1))
tau_aff_n = sp.simplify((1 - Phi0a_n) / Phi1a_n)
tau_quad_n = sp.simplify(2 * (1 - Phi0a_n) / (Phi1a_n - sp.sqrt(Delta_aff_n)))
res_aff_n = sp.simplify(
    (Phi0a_n - 1) + Phi1a_n * tau_quad_n + sp.Rational(1, 2) * Phi2a_n * tau_quad_n**2
)
expect_zero("quadratic affine residual (negative slope)", res_aff_n)
expect_zero(
    "limit Phi2 -> 0 gives tau_aff (negative slope)",
    sp.simplify(sp.limit(tau_quad_n, Phi2a_n, 0) - tau_aff_n),
)

# ---------------------------------------------------------------------------
# III. Exact quadratic logarithmic predictor (positive-slope continuation branch)
# ---------------------------------------------------------------------------
subbanner("III. Exact quadratic logarithmic predictor")
Phi0l, L0l, L1l = sp.symbols("Phi0l L0l L1l", positive=True, real=True)
Delta_log = sp.simplify(L0l**2 - 2 * L1l * sp.log(Phi0l))
tau_log = sp.simplify(-sp.log(Phi0l) / L0l)
tau_log2 = sp.simplify(-2 * sp.log(Phi0l) / (L0l + sp.sqrt(Delta_log)))
res_log = sp.simplify(sp.log(Phi0l) + L0l * tau_log2 + sp.Rational(1, 2) * L1l * tau_log2**2)

print("Delta_log =")
sp.pprint(Delta_log)
print("tau_log =")
sp.pprint(tau_log)
print("tau_log2 =")
sp.pprint(tau_log2)
expect_zero("quadratic log residual", res_log)
expect_zero("limit L1 -> 0 gives tau_log", sp.simplify(sp.limit(tau_log2, L1l, 0) - tau_log))

Phi0l_n, L0l_abs, L1l_n = sp.symbols("Phi0l_n L0l_abs L1l_n", positive=True, real=True)
L0l_n = -L0l_abs
Delta_log_n = sp.simplify(L0l_n**2 - 2 * L1l_n * sp.log(Phi0l_n))
tau_log_n = sp.simplify(-sp.log(Phi0l_n) / L0l_n)
tau_log2_n = sp.simplify(-2 * sp.log(Phi0l_n) / (L0l_n - sp.sqrt(Delta_log_n)))
res_log_n = sp.simplify(
    sp.log(Phi0l_n) + L0l_n * tau_log2_n + sp.Rational(1, 2) * L1l_n * tau_log2_n**2
)
expect_zero("quadratic log residual (negative slope)", res_log_n)
expect_zero(
    "limit L1 -> 0 gives tau_log (negative slope)",
    sp.simplify(sp.limit(tau_log2_n, L1l_n, 0) - tau_log_n),
)

# ---------------------------------------------------------------------------
# IV. Turning-point / tangency formulas
# ---------------------------------------------------------------------------
subbanner("IV. Turning-point and tangency formulas")
Phi0t, Phi2t = sp.symbols("Phi0t Phi2t", real=True)
radicand_tp = sp.simplify(2 * (1 - Phi0t) / Phi2t)
criterion_tp = sp.simplify((1 - Phi0t) * Phi2t)
expect_zero(
    "turning-point radicand/product sign bridge",
    sp.simplify(radicand_tp * Phi2t**2 / 2 - criterion_tp),
)

apos, bpos = sp.symbols("apos bpos", positive=True, real=True)
positive_product_cases = [
    radicand_tp.subs({1 - Phi0t: apos, Phi2t: bpos}),
    radicand_tp.subs({1 - Phi0t: -apos, Phi2t: -bpos}),
]
negative_product_cases = [
    radicand_tp.subs({1 - Phi0t: apos, Phi2t: -bpos}),
    radicand_tp.subs({1 - Phi0t: -apos, Phi2t: bpos}),
]
print("turning-point radicand on (1-Phi0)*Phi2 > 0 cases =")
sp.pprint([sp.simplify(case) for case in positive_product_cases])
if not all(sp.ask(sp.Q.positive(sp.simplify(case))) for case in positive_product_cases):
    raise AssertionError("turning-point positive criterion did not imply a real radicand")
print("turning-point radicand on (1-Phi0)*Phi2 < 0 cases =")
sp.pprint([sp.simplify(case) for case in negative_product_cases])
if not all(sp.ask(sp.Q.negative(sp.simplify(case))) for case in negative_product_cases):
    raise AssertionError("turning-point negative criterion did not imply no real root")

tau_tp = sp.sqrt(radicand_tp)
res_tp_plus = sp.simplify((Phi0t - 1) + sp.Rational(1, 2) * Phi2t * tau_tp**2)
expect_zero("turning-point root residual on real-root branch", res_tp_plus)

Phi0g, Phi1g, Phi2g = sp.symbols("Phi0g Phi1g Phi2g", real=True)
quad_model_g = sp.simplify((Phi0g - 1) + Phi1g * tau + sp.Rational(1, 2) * Phi2g * tau**2)
Delta_tangent = sp.simplify(sp.Rational(1, 2) * Phi2g * tau**2)
print("tangency model Delta(tau) at Phi0=1, Phi1=0 =")
sp.pprint(Delta_tangent)
expect_zero(
    "tangency model from quadratic closure",
    sp.simplify(quad_model_g.subs({Phi0g: 1, Phi1g: 0}) - Delta_tangent),
)

# ---------------------------------------------------------------------------
# V. Small-defect expansions around the closure slice
# ---------------------------------------------------------------------------
subbanner("V. Curvature-corrected local expansions")
eps, L0e, L1e = sp.symbols("eps L0e L1e", positive=True, real=True)
Phi0e = sp.simplify(1 + eps)
Phi1e = sp.simplify(Phi0e * L0e)
Phi2e = sp.simplify(Phi0e * (L1e + L0e**2))

tau_aff_e = sp.simplify((1 - Phi0e) / Phi1e)
tau_quad_e = sp.simplify(2 * (1 - Phi0e) / (Phi1e + sp.sqrt(Phi1e**2 - 2 * Phi2e * (Phi0e - 1))))
tau_log_e = sp.simplify(-sp.log(Phi0e) / L0e)
tau_log2_e = sp.simplify(-2 * sp.log(Phi0e) / (L0e + sp.sqrt(L0e**2 - 2 * L1e * sp.log(Phi0e))))

corr_aff = sp.simplify(
    sp.series(tau_quad_e - tau_aff_e + Phi2e * tau_aff_e**2 / (2 * Phi1e), eps, 0, 3).removeO()
)
print("series[tau_quad - tau_aff + (Phi2/(2 Phi1)) tau_aff^2] =")
sp.pprint(corr_aff)
expect_zero("ordinary quadratic correction formula", corr_aff)

corr_log = sp.simplify(
    sp.series(tau_log2_e - tau_log_e + (L1e / (2 * L0e)) * tau_log_e**2, eps, 0, 3).removeO()
)
print("series[tau_log2 - tau_log + (L1/(2 L0)) tau_log^2] =")
sp.pprint(corr_log)
expect_zero("logarithmic quadratic correction formula", corr_log)

predictor_gap = sp.simplify(sp.series(tau_log2_e - tau_quad_e, eps, 0, 4).removeO())
print("series[tau_log2 - tau_quad] =")
sp.pprint(predictor_gap)
expect_zero(
    "quadratic predictors agree through O(eps^2)",
    sp.simplify(sp.series((tau_log2_e - tau_quad_e) / eps**3, eps, 0, 1).removeO())
    - (L0e**2 + 3 * L1e) / (6 * L0e**3),
)

banner("STAGE 205 SYMPY AUDIT PASSED")
