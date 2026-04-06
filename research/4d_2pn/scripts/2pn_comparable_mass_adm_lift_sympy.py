
#!/usr/bin/env python3
"""
2PN comparable-mass ADM lift prototype
--------------------------------------
Purpose
-------
Continue the 2PN program past the one-body DtN closure and determine the
minimal additional comparable-mass cross sector needed to match the standard
ADM 2PN Hamiltonian.

Method
------
1. Freeze the already-passed sectors:
   - Newtonian + full frozen 1PN Lagrangian,
   - DtN-corrected 2PN one-body/self sector,
   - quadratic local mass scaling with lambda_rho = 1/2.
2. Use the generic-frame ADM Hamiltonians H_1PN and H_2PN from Appendix C
   (Eqs. C.7-C.8) of the Schäfer--Jaranowski 2018 review as the target.
3. Use the perturbative Legendre transform for
       L = L0 + eps L1 + eps^2 L2
   with quadratic L0:
       H1 = -L1(v0),
       H2 = -L2(v0) + 1/2 A0^T M^{-1} A0,
   where A0 = (∂L1/∂v)|_{v=v0}, v0 = p/m.
4. Compute the H_2PN residual of the DtN-corrected self/static candidate.
5. Solve that residual against a compact invariant 2PN cross-sector basis.
6. Rebuild the full candidate and verify exact H_2PN ADM matching.

This is still a target-matching / bookkeeping result, not yet a constructive
wake derivation.
"""
import sympy as sp

def banner(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)

# ---------------------------------------------------------------------------
# Symbols and helper definitions
# ---------------------------------------------------------------------------

G, r = sp.symbols('G r', positive=True)
m1, m2 = sp.symbols('m1 m2', positive=True)

v1x, v1y, v1z, v2x, v2y, v2z = sp.symbols('v1x v1y v1z v2x v2y v2z')
p1x, p1y, p1z, p2x, p2y, p2z = sp.symbols('p1x p1y p1z p2x p2y p2z')

v1 = sp.Matrix([v1x, v1y, v1z])
v2 = sp.Matrix([v2x, v2y, v2z])
p1 = sp.Matrix([p1x, p1y, p1z])
p2 = sp.Matrix([p2x, p2y, p2z])

def dot(a, b):
    return sp.expand(sum(ai * bi for ai, bi in zip(a, b)))

# choose n along z for coefficient extraction
a = dot(v1, v1)
b = dot(v2, v2)
c = dot(v1, v2)
d = v1z
e = v2z

p1sq = dot(p1, p1)
p2sq = dot(p2, p2)
pp = dot(p1, p2)
np1 = p1z
np2 = p2z

# ---------------------------------------------------------------------------
# Frozen Lagrangian sectors
# ---------------------------------------------------------------------------

L0 = sp.Rational(1, 2) * m1 * a + sp.Rational(1, 2) * m2 * b + G * m1 * m2 / r

L1 = (
    (m1 * a**2 + m2 * b**2) / 8
    + (G * m1 * m2 / r) * (
        sp.Rational(3, 2) * (a + b)
        - sp.Rational(7, 2) * c
        - sp.Rational(1, 2) * d * e
    )
    - G**2 * m1 * m2 * (m1 + m2) / (2 * r**2)
)

# DtN-corrected self/static 2PN sector:
#  - free  v^6 coefficient     = 1/16
#  - self  U v^4 coefficient   = 7/8
#  - self  U^2 v^2 coefficient = 2
#  - static U^3 coefficient    = 1/4  (lambda_rho = 1/2)
L2_self_static = (
    (m1 * a**3 + m2 * b**3) / 16
    + (7 * G * m1 * m2 / (8 * r)) * (a**2 + b**2)
    + (2 * G**2 * m1 * m2 / r**2) * (m2 * a + m1 * b)
    + G**3 * m1 * m2 * (m1**2 + m2**2) / (4 * r**3)
)

# ---------------------------------------------------------------------------
# Perturbative Legendre transform at 2PN
# ---------------------------------------------------------------------------

v0subs = {
    v1x: p1x / m1, v1y: p1y / m1, v1z: p1z / m1,
    v2x: p2x / m2, v2y: p2y / m2, v2z: p2z / m2,
}
subs_p_to_mv = {
    p1x: m1 * v1x, p1y: m1 * v1y, p1z: m1 * v1z,
    p2x: m2 * v2x, p2y: m2 * v2y, p2z: m2 * v2z,
}

Minv = sp.diag(1 / m1, 1 / m1, 1 / m1, 1 / m2, 1 / m2, 1 / m2)
A0 = sp.Matrix([sp.diff(L1, vv).subs(v0subs) for vv in [v1x, v1y, v1z, v2x, v2y, v2z]])
quad = sp.expand((A0.T * Minv * A0)[0] / 2)

H0_model = sp.expand(p1sq / (2 * m1) + p2sq / (2 * m2) - G * m1 * m2 / r)
H1_model = sp.expand(-L1.subs(v0subs))
H2_model_no_cross = sp.expand(-L2_self_static.subs(v0subs) + quad)

# ---------------------------------------------------------------------------
# ADM Hamiltonian targets (generic frame)
# Appendix C, Eqs. (C.7)-(C.8), in the notation of the 2018 review
# ---------------------------------------------------------------------------

H1_target = sp.expand(
    -(p1sq**2) / (8 * m1**3) - (p2sq**2) / (8 * m2**3)
    + (G * m1 * m2 / (4 * r)) * (
        -6 * p1sq / m1**2
        - 6 * p2sq / m2**2
        + 14 * pp / (m1 * m2)
        + 2 * np1 * np2 / (m1 * m2)
    )
    + G**2 * m1 * m2 * (m1 + m2) / (2 * r**2)
)

H2_base = (
    p1sq**3 / (16 * m1**5)
    + G * m1 * m2 / (8 * r) * (
        5 * p1sq**2 / m1**4
        - sp.Rational(11, 2) * p1sq * p2sq / (m1**2 * m2**2)
        - pp**2 / (m1**2 * m2**2)
        + 5 * p1sq * np2**2 / (m1**2 * m2**2)
        - 6 * pp * np1 * np2 / (m1**2 * m2**2)
        - sp.Rational(3, 2) * np1**2 * np2**2 / (m1**2 * m2**2)
    )
    + G**2 * m1 * m2 / (4 * r**2) * (
        m2 * (10 * p1sq / m1**2 + 19 * p2sq / m2**2)
        - sp.Rational(1, 2) * (m1 + m2) * (27 * pp + 6 * np1 * np2) / (m1 * m2)
    )
    - G**3 * m1 * m2 * (m1**2 + 5 * m1 * m2 + m2**2) / (8 * r**3)
)
swap = {m1: m2, m2: m1, p1x: p2x, p1y: p2y, p1z: p2z, p2x: p1x, p2y: p1y, p2z: p1z}
H2_target = sp.expand(H2_base + H2_base.xreplace(swap))

# ---------------------------------------------------------------------------
# Basic target checks
# ---------------------------------------------------------------------------

banner("1) Frozen 1PN and DtN-corrected 2PN setup")

print("Model H1 - target H1 =")
print(sp.expand(H1_model - H1_target))

residual_no_cross = sp.expand(H2_model_no_cross - H2_target)

print("\nNo-cross 2PN comparable-mass residual (model - ADM target):")
print(residual_no_cross)

# ---------------------------------------------------------------------------
# Solve the missing comparable-mass cross sector
# ---------------------------------------------------------------------------

# Required additional Lagrangian block satisfies:
#    H2_new = H2_old - L2_cross(v0)
# so we need
#    L2_cross(v0) = residual_no_cross
# therefore
#    L2_cross(v) = residual_no_cross with p -> m v.
L2_cross_required = sp.expand(residual_no_cross.subs(subs_p_to_mv))

banner("2) Solve a compact invariant 2PN cross-sector basis")

quartic_basis = [
    m1 * m2 * c * (a + b),                     # q1
    m1 * m2 * d * e * (a + b),                 # q2
    m1 * m2 * a * b,                           # q3
    m1 * m2 * c**2,                            # q4
    m1 * m2 * (a * e**2 + b * d**2),           # q5
    m1 * m2 * c * d * e,                       # q6
    m1 * m2 * d**2 * e**2,                     # q7
]

quadratic_basis = [
    m1 * m2 * (m2 * a + m1 * b),               # t1
    m1 * m2 * (m1 * a + m2 * b),               # t2
    m1 * m2 * (m1 + m2) * c,                   # t3
    m1 * m2 * (m1 + m2) * d * e,               # t4
    m1 * m2 * (m2 * d**2 + m1 * e**2),         # t5
    m1 * m2 * (m1 * d**2 + m2 * e**2),         # t6
]

static_basis = [m1**2 * m2**2]                # s1

q_syms = sp.symbols('q1:8')
t_syms = sp.symbols('t1:7')
s_sym = sp.symbols('s1')

cross_ansatz = sp.expand(
    G / r * sum(ci * bi for ci, bi in zip(q_syms, quartic_basis))
    + G**2 / r**2 * sum(ci * bi for ci, bi in zip(t_syms, quadratic_basis))
    + G**3 / r**3 * s_sym * static_basis[0]
)

poly = sp.Poly(sp.expand(L2_cross_required - cross_ansatz), v1x, v1y, v1z, v2x, v2y, v2z)
equations = [sp.Eq(coeff, 0) for _, coeff in poly.terms()]
solutions = sp.solve(equations, list(q_syms) + list(t_syms) + [s_sym], dict=True)

if not solutions:
    raise RuntimeError("No coefficient solution found for the chosen cross basis.")

sol = solutions[0]

print("Unique solution in the chosen basis:")
for sym in list(q_syms) + list(t_syms) + [s_sym]:
    print(f"  {sym} = {sp.simplify(sol[sym])}")

# Build explicit cross block
L2_cross = sp.expand(cross_ansatz.subs(sol))

print("\nRequired added cross block L2_cross:")
print(L2_cross)

# ---------------------------------------------------------------------------
# Final exact match check
# ---------------------------------------------------------------------------

banner("3) Full 2PN candidate and exact ADM-H2PN match")

L2_full = sp.expand(L2_self_static + L2_cross)
H2_full = sp.expand(-L2_full.subs(v0subs) + quad)

print("Full H2 - target H2 =")
print(sp.expand(H2_full - H2_target))

# ---------------------------------------------------------------------------
# Display the clean invariant form of the added cross block and the total L2 block
# ---------------------------------------------------------------------------

banner("4) Clean invariant form")

L2_cross_clean = sp.expand(
    G * m1 * m2 / r * (
        - sp.Rational(7, 4) * c * (a + b)
        - sp.Rational(1, 4) * d * e * (a + b)
        + sp.Rational(11, 8) * a * b
        + sp.Rational(1, 4) * c**2
        - sp.Rational(5, 8) * (a * e**2 + b * d**2)
        + sp.Rational(3, 2) * c * d * e
        + sp.Rational(3, 8) * d**2 * e**2
    )
    + G**2 * m1 * m2 / r**2 * (
        sp.Rational(11, 8) * (m1 * a + m2 * b)
        - sp.Rational(15, 4) * (m1 + m2) * c
        + sp.Rational(15, 8) * (m1 * d**2 + m2 * e**2)
    )
    + sp.Rational(5, 4) * G**3 * m1**2 * m2**2 / r**3
)

print("Added cross block (clean form) =")
print(L2_cross_clean)

L2_total_clean = sp.expand(L2_self_static + L2_cross_clean)
print("\nTotal 2PN Lagrangian block L2_total =")
print(L2_total_clean)

# test-mass sanity
L2_cross_testmass = sp.simplify(L2_cross_clean.subs({v2x: 0, v2y: 0, v2z: 0}) / m1)
print("\nCross block per unit test mass with body 2 at rest, m1 -> 0:")
print(sp.limit(L2_cross_testmass, m1, 0))

banner("5) Takeaways")
print("- The DtN-corrected one-body closure removes the old U^2 v^2 test-mass mismatch.")
print("- The remaining comparable-mass problem is linear at the 2PN Legendre stage.")
print("- In the chosen 7+6+1 basis, the coefficient solve is unique.")
print("- Three quadratic basis coefficients vanish automatically: t1 = t4 = t5 = 0.")
print("- The added static cross term is +5/4 G^3 m1^2 m2^2 / r^3,")
print("  which upgrades the full Lagrangian static mass polynomial to")
print("      +(G^3 m1 m2 / 4 r^3) (m1^2 + 5 m1 m2 + m2^2).")
print("- With this added block, the full candidate reproduces the generic-frame ADM H2PN target exactly.")
