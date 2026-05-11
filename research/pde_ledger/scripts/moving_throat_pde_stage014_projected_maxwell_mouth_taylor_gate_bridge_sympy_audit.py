#!/usr/bin/env python3
"""
Projected Maxwell → mouth-Taylor → 5PN gate bridge.

This script starts from the previously derived one-port projected-Maxwell primitive
bridge and adds a concrete near-mouth Taylor ansatz:

    X_proj(ell) = X_0 + ell * mu1 * X_0' + O(ell^2)

for the primitive one-port packet X in {Q, S2, H_port, Delta, P, G_W}.

It then pushes those slippages into the moving-throat 5PN bottlenecks

    Xi_load = P1/P0,
    K1      = D21 + D01/9,
    H_even  = D41 - (2/3) D21 - D01/27,

and the constant-prefactor transport slots delta P2, delta P4.
"""

import sympy as sp


def assert_zero(label, expr):
    residue = sp.factor(sp.together(sp.simplify(expr)))
    if residue != 0:
        raise AssertionError(f"{label} failed: {sp.sstr(residue)}")


def assert_nonzero(label, expr):
    residue = sp.factor(sp.together(sp.simplify(expr)))
    if residue == 0:
        raise AssertionError(f"{label} unexpectedly vanished")


# Base one-port variables
Q, S2, Hport, Delta, P, Gw = sp.symbols('Q S2 H Delta P Gw', nonzero=True)
D0, D2, D4, N0, Ptarget = sp.symbols('D0 D2 D4 N0 Ptarget', nonzero=True)
mu1 = sp.symbols('mu1')

# Primitive slippages
q1, s1, h1, d1, p1, g1 = sp.symbols('q1 s1 h1 d1 p1 g1')

# Mouth derivatives
Qx, Sx, Hx, Dx, Px, Gx = sp.symbols('Qx S2x Hx Deltax Px Gwx')

# Primitive one-port formulas
z0 = (Delta*q1 - Q*d1)/Delta**2
z2 = (-Delta**2*h1 + Delta*(Hport*d1 + Q*s1 + S2*q1) - 2*Q*S2*d1)/Delta**3
z4 = (-Delta**2*Hport*s1 - Delta**2*S2*h1 - Delta**2*q1 + 2*Delta*Hport*S2*d1 + 2*Delta*Q*S2*s1 + 2*Delta*Q*d1 + Delta*S2**2*q1 - 3*Q*S2**2*d1)/Delta**4

n0 = 2*P*(Delta*p1 - P*d1)/Delta**3
n2 = -(2*Delta**2*(Gw*p1 + P*g1) - 2*Delta*P*(2*Gw*d1 + P*s1 + 2*S2*p1) + 6*P**2*S2*d1)/Delta**4
n4 = 2*(Delta**3*Gw*g1 - Delta**2*Gw**2*d1 - 2*Delta**2*Gw*P*s1 - 2*Delta**2*Gw*S2*p1 - 2*Delta**2*P*S2*g1 - 2*Delta**2*P*p1 + 6*Delta*Gw*P*S2*d1 + 3*Delta*P**2*S2*s1 + 3*Delta*P**2*d1 + 3*Delta*P*S2**2*p1 - 6*P**2*S2**2*d1)/Delta**5

# Mouth-local one-sided Taylor ansatz
subs_der = {q1: mu1*Qx, s1: mu1*Sx, h1: mu1*Hx, d1: mu1*Dx, p1: mu1*Px, g1: mu1*Gx}

z0d = sp.simplify(z0.subs(subs_der)/mu1)
z2d = sp.simplify(z2.subs(subs_der)/mu1)
z4d = sp.simplify(z4.subs(subs_der)/mu1)
n0d = sp.simplify(n0.subs(subs_der)/mu1)
n2d = sp.simplify(n2.subs(subs_der)/mu1)
n4d = sp.simplify(n4.subs(subs_der)/mu1)

# 5PN bottlenecks
Xi = sp.simplify((2*p1/P - 2*d1/Delta + q1/(D0*Delta) - Q*d1/(D0*Delta**2)).subs(subs_der)/mu1)
K1 = sp.simplify((-(z2 + z0/sp.Integer(9))).subs(subs_der)/mu1)
He = sp.simplify((-(z4) + sp.Rational(2,3)*z2 - z0/sp.Integer(27)).subs(subs_der)/mu1)

S, T = sp.symbols('S T', nonzero=True)
Compat = sp.simplify(((n0/Ptarget) - 6*S*z2/T + 3*S**2*z4/T**2).subs(subs_der)/mu1)

deltaP2 = sp.simplify(((D0**2*n2 - 2*D0*D2*n0 + 2*D0*N0*z2 - 2*D2*N0*z0)/D0**3).subs(subs_der)/mu1)
deltaP4 = sp.simplify(((D0**3*n4 - 2*D0**2*D2*n2 - 2*D0**2*D4*n0 + 2*D0**2*N0*z4 + 3*D0*D2**2*n0 - 2*D0*D2*N0*z2 - 2*D0*D4*N0*z0 + 2*D2**2*N0*z0)/D0**4).subs(subs_der)/mu1)

# Mechanism sieve
qd_only = sp.solve([sp.Eq(K1.subs({Sx:0, Hx:0}), 0), sp.Eq(He.subs({Sx:0, Hx:0}), 0)], [Qx, Dx], dict=True)
sh_only = sp.solve([sp.Eq(K1.subs({Qx:0, Dx:0}), 0), sp.Eq(He.subs({Qx:0, Dx:0}), 0)], [Sx, Hx], dict=True)
comp_surface = sp.solve([sp.Eq(K1, 0), sp.Eq(He, 0)], [Hx, Sx], dict=True)[0]

qd_matrix = sp.Matrix([
    [sp.diff(K1.subs({Sx: 0, Hx: 0}), Qx), sp.diff(K1.subs({Sx: 0, Hx: 0}), Dx)],
    [sp.diff(He.subs({Sx: 0, Hx: 0}), Qx), sp.diff(He.subs({Sx: 0, Hx: 0}), Dx)],
])
sh_matrix = sp.Matrix([
    [sp.diff(K1.subs({Qx: 0, Dx: 0}), Sx), sp.diff(K1.subs({Qx: 0, Dx: 0}), Hx)],
    [sp.diff(He.subs({Qx: 0, Dx: 0}), Sx), sp.diff(He.subs({Qx: 0, Dx: 0}), Hx)],
])
Hx_den = sp.factor(sp.denom(sp.together(comp_surface[Hx])))
Sx_den = sp.factor(sp.denom(sp.together(comp_surface[Sx])))

coef_Xi_Px = sp.simplify(sp.diff(Xi, Px))
coef_dP2_Gx = sp.simplify(sp.diff(deltaP2, Gx))
coef_dP4_Gx = sp.simplify(sp.diff(deltaP4, Gx))

assert_zero("z0 derivative map", z0d - (Delta*Qx - Q*Dx)/Delta**2)
assert_zero("n0 derivative map", n0d - 2*P*(Delta*Px - P*Dx)/Delta**3)

for sym in (Sx, Hx, Gx):
    assert_zero(f"Xi_load independence from {sym}", sp.diff(Xi, sym))
for sym in (Px, Gx):
    assert_zero(f"K1 independence from {sym}", sp.diff(K1, sym))
    assert_zero(f"H_even independence from {sym}", sp.diff(He, sym))
assert_zero("d Xi_load / d Pprime", coef_Xi_Px - 2/P)
assert_zero("d deltaP2 / d G_W prime", coef_dP2_Gx + 2*P/(D0*Delta**2))
assert_nonzero("deltaP4 should depend on G_W prime", coef_dP4_Gx)
assert_nonzero("source/denominator sieve determinant", qd_matrix.det())
assert_nonzero("spectral sieve determinant", sh_matrix.det())
assert_zero("compensation Hport denominator", Hx_den - 9*Delta**2*(Delta*Hport - Q*S2))
assert_zero("compensation S2 denominator", Sx_den - 9*Delta*(Delta*Hport - Q*S2))
assert_nonzero("mutated compensation denominator should fail", Hx_den - 9*Delta**2*(Delta*Hport + Q*S2))

if qd_only != [{Qx: 0, Dx: 0}]:
    raise AssertionError(f"Unexpected pure source/denominator solve: {qd_only}")
if sh_only != [{Sx: 0, Hx: 0}]:
    raise AssertionError(f"Unexpected pure spectral solve: {sh_only}")
assert_zero("compensation K1", K1.subs(comp_surface))
assert_zero("compensation H_even", He.subs(comp_surface))
assert_zero("compensation denominator tracks Z2 slot", (Delta*Hport - Q*S2) + Delta**2*((Q*S2 - Hport*Delta)/Delta**2))

def p(expr):
    return sp.factor(sp.simplify(expr))

print("="*96)
print("1) Mouth-local Taylor ansatz")
print("="*96)
print("Take a normalized one-sided mouth kernel w(u), u >= 0, with first moment mu1 = ∫ u w(u) du.")
print("Then for any primitive one-port quantity X(s),")
print("  X_proj(ell) = X(0) + ell * mu1 * X'(0) + O(ell^2).")
print("So")
print("  q1 = mu1*Q'(0),   s1 = mu1*S2'(0),   h1 = mu1*H_port'(0),")
print("  d1 = mu1*Delta'(0), p1 = mu1*P'(0),  g1 = mu1*G_W'(0).")
print("For a symmetric interior kernel, mu1 = 0 and the whole first layer vanishes.\n")

print("="*96)
print("2) Primitive-to-bundle mouth-derivative map (divide out mu1)")
print("="*96)
print("z0/mu1 =")
print(" ", p(z0d))
print("z2/mu1 =")
print(" ", p(z2d))
print("z4/mu1 =")
print(" ", p(z4d))
print("n0/mu1 =")
print(" ", p(n0d))
print("n2/mu1 =")
print(" ", p(n2d))
print("n4/mu1 =")
print(" ", p(n4d))
print()

print("="*96)
print("3) Direct projected-Maxwell contribution to the 5PN bottlenecks")
print("="*96)
print("Xi_load^(proj)/mu1 =")
print(" ", p(Xi))
print("K1^(proj)/mu1 =")
print(" ", p(K1))
print("H_even^(proj)/mu1 =")
print(" ", p(He))
print()
print("Immediate dependency pattern:")
print("  - Xi_load depends only on (Q', Delta', P').")
print("  - K1 and H_even depend only on (Q', Delta', S2', H_port').")
print("  - G_W' does not enter Xi_load, K1, or H_even.")
print("  - G_W' first enters the constant-prefactor transport delta P2, delta P4.\n")

print("="*96)
print("4) Isotropic compatibility and constant-prefactor transport")
print("="*96)
print("delta Compat / mu1 =")
print(" ", p(Compat))
print("delta P2 / mu1 =")
print(" ", p(deltaP2))
print("delta P4 / mu1 =")
print(" ", p(deltaP4))
print()
print("Useful direct coefficients:")
print("  d(Xi_load)/dP'  =")
print("   ", p(coef_Xi_Px))
print("  d(delta P2)/dG_W' =")
print("   ", p(coef_dP2_Gx))
print("  d(delta P4)/dG_W' =")
print("   ", p(coef_dP4_Gx))
print()

print("="*96)
print("5) Mechanism sieve")
print("="*96)
print("Case A: source/denominator sector only, S2' = H_port' = 0")
print("  solve K1 = H_even = 0  ->")
print("  ", qd_only)
print("Case B: spectral sector only, Q' = Delta' = 0")
print("  solve K1 = H_even = 0  ->")
print("  ", sh_only)
print()
print("So neither a pure (Q', Delta') correction nor a pure (S2', H_port') correction")
print("can close the even gates nontrivially.")
print()

print("="*96)
print("6) Exact mixed compensation surface for the even gates")
print("="*96)
print("Solving K1 = H_even = 0 for (H_port', S2') gives")
print("H_port' =")
print(" ", p(comp_surface[Hx]))
print("S2' =")
print(" ", p(comp_surface[Sx]))
print()
print("The compensation denominator is proportional to")
print("  Delta*H_port - Q*S2.")
print("This is exactly the primitive combination behind the conservative Z2 slot:")
print("  Z2 = (Q*S2 - H_port*Delta)/Delta^2.")
print("So the even-gate compensation becomes singular on the surface Delta*H_port = Q*S2")
print("(equivalently Z2 = 0).")
print()

print("="*96)
print("7) Readout")
print("="*96)
print("The mouth-local projected-Maxwell packet factorizes into three jobs:")
print("  (i)  (Q', Delta', S2', H_port') repair the conservative even gates K1 and H_even,")
print("  (ii) P' directly tunes the weak-axisymmetric prefactor slope Xi_load,")
print("  (iii) G_W' first tunes the constant-prefactor transport through delta P2 and delta P4.")
print()
print("This is precisely the kind of multi-channel structure you would expect from a real")
print("near-throat missing piece, and it is much sharper than a single scalar correction.")
print("STATUS: PASS")
