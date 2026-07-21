#!/usr/bin/env python3
"""Independent SymPy recomputation from the three committed text sources.

The two-body kernel is represented by an off-diagonal Green response G*m_gg.
Differentiating the self-subtracted on-shell functional at G=0 extracts the
coefficient A in U_12 = s1*s2*A/(4*pi*R), without importing a force sign.
"""

import sympy as sp


# Stable coupled scalar block and its reduced far-field response matrix.
b, k = sp.symbols("b k", positive=True)
c, z_b = sp.symbols("c z_b", real=True)
z_g = sp.symbols("z_g", positive=True)
D = b * k - c**2
kappa = D / b
Z = sp.Matrix([[1, 0], [-(c / b) * z_b, z_g]])
m = sp.simplify(Z.T * sp.diag(1 / b, 1 / kappa) * Z)
m_gg = sp.factor(m[1, 1])
det_m = sp.factor(m.det())

assert sp.simplify(m_gg - b * z_g**2 / D) == 0
assert sp.simplify(det_m - z_g**2 / D) == 0

# Make D>0 explicit to let SymPy certify m_gg>0.
Delta = sp.symbols("Delta", positive=True)
m_gg_stable = sp.simplify(m_gg.subs(k, (Delta + c**2) / b))
assert m_gg_stable == b * z_g**2 / Delta
assert sp.ask(sp.Q.positive(m_gg_stable)) is True

witness_subs = {
    b: sp.Integer(2),
    k: sp.Integer(1),
    c: sp.Rational(1, 2),
    z_b: sp.Integer(1),
    z_g: sp.Integer(1),
}
assert D.subs(witness_subs) == sp.Rational(7, 4)
assert m.subs(witness_subs) == sp.Matrix(
    [[sp.Rational(4, 7), sp.Rational(-2, 7)],
     [sp.Rational(-2, 7), sp.Rational(8, 7)]]
)

# K maps oriented source/flux amplitudes to one-body plus cross responses.
# G is a bookkeeping variable for 1/(4*pi*R); S is the positive self response.
G, S = sp.symbols("G S", positive=True)
phi, q, j, g = sp.symbols("phi q j g", positive=True)
lam = sp.symbols("lambda", real=True)
s1, s2 = sp.symbols("s1 s2", nonzero=True, real=True)
orient = sp.Matrix([s1, s2])
K = sp.Matrix([[S, m_gg * G], [m_gg * G, S]])


def pair_coefficient(functional):
    """Remove the G=0 self part and extract the leading s1*s2*G term."""
    interaction = sp.simplify(functional - functional.subs(G, 0))
    return sp.factor(sp.diff(interaction, G).subs(G, 0) / (s1 * s2))


# V: Q is the reaction required to keep the value fixed.  The correct member
# is E0-Q.phi; evaluating E0 alone is the bare/wrong mutation.
held_value = phi * orient
reaction = K.inv() * held_value
E0_V = sp.Rational(1, 2) * (reaction.T * K * reaction)[0]
F_V = sp.simplify(E0_V - (reaction.T * held_value)[0])
A_V = pair_coefficient(F_V)
A_V_bare = pair_coefficient(E0_V)

# M: q is held flux, while g remains a fixed source whose work must be
# subtracted.  With g=0 this conjugate is exactly the bare stored energy.
total_M = (q + g) * orient
field_M = K * total_M
E0_M = sp.Rational(1, 2) * (total_M.T * field_M)[0]
F_M = sp.simplify(E0_M - (g * orient).T.dot(field_M))
A_M = pair_coefficient(F_M)
A_M_wrong = pair_coefficient(E0_M)

# J: both j and g are fixed sources, so all their work is subtracted.
total_J = (j + g) * orient
field_J = K * total_J
E0_J = sp.Rational(1, 2) * (total_J.T * field_J)[0]
F_J = sp.simplify(E0_J - total_J.T.dot(field_J))
A_J = pair_coefficient(F_J)
A_J_bare = pair_coefficient(E0_J)

# MIXED: lambda=0 subtracts only committed-source work (M); lambda=1 also
# subtracts all q work (the J endpoint with j=q).
mixed_work = (g + lam * q) * orient
F_MIXED = sp.simplify(E0_M - mixed_work.T.dot(field_M))
A_MIXED = pair_coefficient(F_MIXED)

claimed_A_V = m_gg * phi**2 / S**2
claimed_A_M = m_gg * (q**2 - g**2)
claimed_A_J = -m_gg * (j + g)**2
claimed_A_MIXED = m_gg * ((1 - 2 * lam) * q**2 - 2 * lam * q * g - g**2)

assert sp.simplify(A_V - claimed_A_V) == 0
assert sp.simplify(A_M - claimed_A_M) == 0
assert sp.simplify(A_J - claimed_A_J) == 0
assert sp.simplify(A_MIXED - claimed_A_MIXED) == 0

# Mutation controls and the special pure-flux M equivalence.
assert sp.simplify(A_V_bare + A_V) == 0
assert sp.simplify(A_J_bare + A_J) == 0
assert sp.simplify(A_M.subs(g, 0) - A_M_wrong.subs(g, 0)) == 0
assert sp.simplify(A_M_wrong - A_M - 2 * m_gg * g * (q + g)) == 0

# Range of MIXED(lambda): it decreases strictly.  It has all three signs iff
# q>g; q=g has a null endpoint and negative interior; q<g is strictly negative.
mixed_factored = sp.factor(A_MIXED)
lambda_zero = sp.solve(sp.Eq(A_MIXED, 0), lam)[0]
assert sp.simplify(mixed_factored - m_gg * (g + q) * (-g - 2 * lam * q + q)) == 0
assert lambda_zero == (q - g) / (2 * q)
assert sp.simplify(sp.diff(A_MIXED, lam) + 2 * m_gg * q * (q + g)) == 0
assert sp.simplify(A_MIXED.subs(lam, 0) - A_M) == 0
assert sp.simplify(A_MIXED.subs(lam, 1) - A_J.subs(j, q)) == 0

# Three-dimensional harmonic Green behavior gives U~1/R and F_out~1/R^2.
R, A = sp.symbols("R A", positive=True)
U = s1 * s2 * A / (4 * sp.pi * R)
F_out = -sp.diff(U, R)
assert F_out == s1 * s2 * A / (4 * sp.pi * R**2)

verdict = [
    "m_gg: CONFIRMED (b*z_g^2/D)",
    "det m: CONFIRMED (z_g^2/D)",
    "D*: CONFIRMED (7/4; m*=(1/7)[[4,-2],[-2,8]])",
    "A_V sign: CONFIRMED (positive for nonzero held magnitude phi)",
    "A_M expression+sign: CONFIRMED (m_gg*(q^2-g^2); indefinite)",
    "A_J sign: CONFIRMED (-m_gg*(j+g)^2 < 0)",
    "A_MIXED range: CONFIRMED (all three signs iff q>g; q=g null/negative; q<g negative)",
    "1/R^2: CONFIRMED",
    "M-wrong-functional: CONFIRMED (pure-flux M=bare; coupled mutation must omit/mishandle held-h work subtraction)",
    "KEY finding: CONFIRMED",
]
assert len(verdict) <= 15
print("\n".join(verdict))
