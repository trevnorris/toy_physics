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


banner("STAGE 193 — THREE-COORDINATE MIXED-SIMPLEX AND THE CANONICAL TRIPLE-SCREEN AUDIT")

# ---------------------------------------------------------------------------
# I. Positive spherical simplex and exact boundary reductions
# ---------------------------------------------------------------------------
subbanner("I. Positive spherical simplex and exact boundary reductions")

ai, aj, ak = sp.symbols("a_i a_j a_k", nonnegative=True, real=True)
r, s, nu = sp.symbols("r s nu", nonnegative=True, real=True)
ki, kj, kk = sp.symbols("k_i k_j k_k", positive=True, real=True)
H0 = sp.symbols("H0", positive=True, real=True)

avec = sp.Matrix([ai, aj, ak])
kvec = sp.Matrix([ki, kj, kk])
k_simplex = sp.simplify((avec.T * kvec)[0])

# Edge parameterizations
avec_ij = sp.Matrix([1 / sp.sqrt(1 + r**2), r / sp.sqrt(1 + r**2), 0])
avec_ik = sp.Matrix([1 / sp.sqrt(1 + s**2), 0, s / sp.sqrt(1 + s**2)])
avec_jk = sp.Matrix([0, 1 / sp.sqrt(1 + nu**2), nu / sp.sqrt(1 + nu**2)])

print("a_ij(r) =")
sp.pprint(avec_ij)
print("a_ik(s) =")
sp.pprint(avec_ik)
print("a_jk(nu) =")
sp.pprint(avec_jk)

expect_zero("edge ij normalization", (avec_ij.T * avec_ij)[0] - 1)
expect_zero("edge ik normalization", (avec_ik.T * avec_ik)[0] - 1)
expect_zero("edge jk normalization", (avec_jk.T * avec_jk)[0] - 1)

expect_zero(
    "edge ij slope reduction",
    k_simplex.subs({ai: avec_ij[0], aj: avec_ij[1], ak: avec_ij[2]}) - (ki + r * kj) / sp.sqrt(1 + r**2),
)
expect_zero(
    "edge ik slope reduction",
    k_simplex.subs({ai: avec_ik[0], aj: avec_ik[1], ak: avec_ik[2]}) - (ki + s * kk) / sp.sqrt(1 + s**2),
)
expect_zero(
    "edge jk slope reduction",
    k_simplex.subs({ai: avec_jk[0], aj: avec_jk[1], ak: avec_jk[2]}) - (kj + nu * kk) / sp.sqrt(1 + nu**2),
)

# ---------------------------------------------------------------------------
# II. Exact three-coordinate gradient-synergy theorem
# ---------------------------------------------------------------------------
subbanner("II. Exact three-coordinate gradient-synergy theorem")

Kgrad = sp.sqrt(ki**2 + kj**2 + kk**2)
avec_grad = sp.Matrix([ki / Kgrad, kj / Kgrad, kk / Kgrad])

print("a_grad =")
sp.pprint(avec_grad)
print("k_grad =")
sp.pprint(Kgrad)

expect_zero("gradient-optimal normalization", (avec_grad.T * avec_grad)[0] - 1)
expect_zero("gradient-optimal slope value", (avec_grad.T * kvec)[0] - Kgrad)
expect_zero("gradient-optimal ratio r", sp.simplify(avec_grad[1] / avec_grad[0] - kj / ki))
expect_zero("gradient-optimal ratio s", sp.simplify(avec_grad[2] / avec_grad[0] - kk / ki))

Kij = sp.sqrt(ki**2 + kj**2)
Kik = sp.sqrt(ki**2 + kk**2)
Kjk = sp.sqrt(kj**2 + kk**2)
expect_zero("Kgrad^2 - Kij^2 - k_k^2", Kgrad**2 - Kij**2 - kk**2)
expect_zero("Kgrad^2 - Kik^2 - k_j^2", Kgrad**2 - Kik**2 - kj**2)
expect_zero("Kgrad^2 - Kjk^2 - k_i^2", Kgrad**2 - Kjk**2 - ki**2)

# ---------------------------------------------------------------------------
# III. Exact total cross-leverage theorem
# ---------------------------------------------------------------------------
subbanner("III. Exact total cross-leverage theorem")

wSigma = sp.simplify(2 * (ai * aj + ai * ak + aj * ak))
print("w_Sigma(a) =")
sp.pprint(wSigma)

expect_zero(
    "w_Sigma - ((sum a)^2 - ||a||^2)",
    wSigma - ((ai + aj + ak) ** 2 - (ai**2 + aj**2 + ak**2)),
)
expect_zero(
    "Cauchy slack identity",
    3 * (ai**2 + aj**2 + ak**2)
    - (ai + aj + ak) ** 2
    - ((ai - aj) ** 2 + (ai - ak) ** 2 + (aj - ak) ** 2),
)

aeq = sp.Matrix([1 / sp.sqrt(3), 1 / sp.sqrt(3), 1 / sp.sqrt(3)])
aeq_pair = sp.Matrix([1 / sp.sqrt(2), 1 / sp.sqrt(2), 0])
expect_zero("equal-mix barycenter normalization", (aeq.T * aeq)[0] - 1)
expect_zero("pairwise equal-edge normalization", (aeq_pair.T * aeq_pair)[0] - 1)
expect_zero("w_Sigma(equal-mix barycenter) - 2", wSigma.subs({ai: aeq[0], aj: aeq[1], ak: aeq[2]}) - 2)
expect_zero("w_Sigma(pairwise equal edge) - 1", wSigma.subs({ai: aeq_pair[0], aj: aeq_pair[1], ak: aeq_pair[2]}) - 1)

# ---------------------------------------------------------------------------
# IV. Exact three-coordinate curvature law and edge reductions
# ---------------------------------------------------------------------------
subbanner("IV. Exact three-coordinate curvature law and edge reductions")

uii, uij, uik, ujj, ujk, ukk = sp.symbols(
    "u_ii u_ij u_ik u_jj u_jk u_kk", real=True
)
Hmat = sp.Matrix(
    [
        [uii, uij, uik],
        [uij, ujj, ujk],
        [uik, ujk, ukk],
    ]
)

kappa_simplex = sp.simplify((avec.T * Hmat * avec)[0])
print("kappa_ijk(a) =")
sp.pprint(kappa_simplex)
expect_zero(
    "edge ij curvature reduction",
    kappa_simplex.subs({ai: avec_ij[0], aj: avec_ij[1], ak: avec_ij[2]})
    - (uii + 2 * uij * r + ujj * r**2) / (1 + r**2),
)
expect_zero(
    "edge ik curvature reduction",
    kappa_simplex.subs({ai: avec_ik[0], aj: avec_ik[1], ak: avec_ik[2]})
    - (uii + 2 * uik * s + ukk * s**2) / (1 + s**2),
)
expect_zero(
    "edge jk curvature reduction",
    kappa_simplex.subs({ai: avec_jk[0], aj: avec_jk[1], ak: avec_jk[2]})
    - (ujj + 2 * ujk * nu + ukk * nu**2) / (1 + nu**2),
)
expect_zero(
    "diagonal-neutral reduction",
    kappa_simplex.subs({uij: 0, uik: 0, ujk: 0}) - (uii * ai**2 + ujj * aj**2 + ukk * ak**2),
)

# ---------------------------------------------------------------------------
# V. Exact fixed-simplex root map and interior ratio-coordinate form
# ---------------------------------------------------------------------------
subbanner("V. Exact fixed-simplex root map and ratio-coordinate form")

tau_simplex = sp.simplify(2 * H0 / (k_simplex + sp.sqrt(k_simplex**2 - 2 * H0 * kappa_simplex)))
print("tau_ijk(a) =")
sp.pprint(tau_simplex)

# Interior ratio patch a_i > 0
avec_rs = sp.Matrix([1 / sp.sqrt(1 + r**2 + s**2), r / sp.sqrt(1 + r**2 + s**2), s / sp.sqrt(1 + r**2 + s**2)])
k_rs = sp.simplify(k_simplex.subs({ai: avec_rs[0], aj: avec_rs[1], ak: avec_rs[2]}))
kappa_rs = sp.simplify(kappa_simplex.subs({ai: avec_rs[0], aj: avec_rs[1], ak: avec_rs[2]}))

A = sp.simplify(ki**2 - 2 * H0 * uii)
B = sp.simplify(2 * ki * kj - 4 * H0 * uij)
C = sp.simplify(2 * ki * kk - 4 * H0 * uik)
D = sp.simplify(kj**2 - 2 * H0 * ujj)
E = sp.simplify(2 * kj * kk - 4 * H0 * ujk)
F = sp.simplify(kk**2 - 2 * H0 * ukk)

Delta_sharp = sp.expand(A + B * r + C * s + D * r**2 + E * r * s + F * s**2)
print("Delta_sharp(r,s) =")
sp.pprint(Delta_sharp)

expect_zero(
    "interior ratio normalization",
    (avec_rs.T * avec_rs)[0] - 1,
)
expect_zero(
    "discriminant numerator reduction",
    sp.simplify((1 + r**2 + s**2) * (k_rs**2 - 2 * H0 * kappa_rs) - Delta_sharp),
)

tau_rs_expected = sp.simplify(
    2 * H0 * sp.sqrt(1 + r**2 + s**2)
    / (ki + r * kj + s * kk + sp.sqrt(Delta_sharp))
)
expect_zero(
    "interior ratio tau form",
    tau_simplex.subs({ai: avec_rs[0], aj: avec_rs[1], ak: avec_rs[2]}) - tau_rs_expected,
)

# Boundary reductions from the interior ratio patch
Aij = sp.simplify(ki**2 - 2 * H0 * uii)
Bij = sp.simplify(2 * ki * kj - 4 * H0 * uij)
Cij = sp.simplify(kj**2 - 2 * H0 * ujj)

Aik = sp.simplify(ki**2 - 2 * H0 * uii)
Bik = sp.simplify(2 * ki * kk - 4 * H0 * uik)
Cik = sp.simplify(kk**2 - 2 * H0 * ukk)

Ajk = sp.simplify(kj**2 - 2 * H0 * ujj)
Bjk = sp.simplify(2 * kj * kk - 4 * H0 * ujk)
Cjk = sp.simplify(kk**2 - 2 * H0 * ukk)

expect_zero(
    "tau(r,s) -> pairwise ij on s=0",
    sp.simplify(tau_rs_expected.subs(s, 0) - 2 * H0 * sp.sqrt(1 + r**2) / (ki + r * kj + sp.sqrt(Aij + Bij * r + Cij * r**2))),
)
expect_zero(
    "tau(r,s) -> pairwise ik on r=0",
    sp.simplify(tau_rs_expected.subs(r, 0) - 2 * H0 * sp.sqrt(1 + s**2) / (ki + s * kk + sp.sqrt(Aik + Bik * s + Cik * s**2))),
)
expect_zero(
    "direct jk edge tau reduction",
    tau_simplex.subs({ai: avec_jk[0], aj: avec_jk[1], ak: avec_jk[2]})
    - 2 * H0 * sp.sqrt(1 + nu**2) / (kj + nu * kk + sp.sqrt(Ajk + Bjk * nu + Cjk * nu**2)),
)

# ---------------------------------------------------------------------------
# VI. Canonical interior screen values
# ---------------------------------------------------------------------------
subbanner("VI. Canonical interior screen values")

k_eq = sp.simplify((ki + kj + kk) / sp.sqrt(3))
k_grad = sp.simplify(Kgrad)
print("k_eq =")
sp.pprint(k_eq)
print("k_grad =")
sp.pprint(k_grad)
expect_zero("equal-mix slope value", (aeq.T * kvec)[0] - k_eq)
expect_zero("gradient slope value", (avec_grad.T * kvec)[0] - k_grad)
expect_zero("k_grad^2 - ||k||^2", k_grad**2 - (ki**2 + kj**2 + kk**2))

banner("STAGE 193 SYMPY AUDIT COMPLETED SUCCESSFULLY")
