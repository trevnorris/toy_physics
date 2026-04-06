#!/usr/bin/env python3
import sympy as sp
import mpmath as mp

mp.mp.dps = 50

X, Y, Z = sp.symbols('X Y Z')
B0, B1, A1, B2, A2, C2 = sp.symbols('B0 B1 A1 B2 A2 C2')
p_iso, q_ax = sp.symbols('p_iso q_ax')

def double_factorial(n: int) -> int:
    if n <= 0:
        return 1
    out = 1
    for k in range(n, 0, -2):
        out *= k
    return out

def sphere_avg_monomial(a: int, b: int, c: int) -> sp.Rational:
    # <x^(2a) y^(2b) z^(2c)> over the unit 2-sphere S^2
    num = double_factorial(2*a - 1) * double_factorial(2*b - 1) * double_factorial(2*c - 1)
    den = double_factorial(2*(a+b+c) + 1)
    return sp.Rational(num, den)

def sphere_avg_poly(poly: sp.Expr) -> sp.Expr:
    poly = sp.expand(poly)
    result = sp.Rational(0)
    for monom, coeff in sp.Poly(poly, X, Y, Z).terms():
        a, b, c = monom
        result += coeff * sphere_avg_monomial(a, b, c)
    return sp.simplify(result)

P2 = (3*Z - 1) / 2

weights = {
    "1perp": X,
    "10": Z,
    "20": ((3*Z - 1) / 2)**2,
    "21": 3*X*Z,
    "22": 3*X*Y,
}

cvals = {}
dvals = {}
for key, w in weights.items():
    den = sphere_avg_poly(w)
    cvals[key] = sp.simplify(sphere_avg_poly(w * P2) / den)
    dvals[key] = sp.simplify(sphere_avg_poly(w * P2**2) / den)

Ytarget = {
    "0": sp.Rational(45, 4),
    "1perp": sp.Rational(7, 2),
    "10": sp.Integer(4),
    "20": sp.Rational(9, 4),
    "21": sp.Rational(3, 2),
    "22": sp.Rational(3, 8),
}
Ztarget = {k: sp.simplify(1/v) for k, v in Ytarget.items()}

sol_odd = sp.solve([
    sp.Eq(B1 + A1*cvals["1perp"], Ztarget["1perp"]),
    sp.Eq(B1 + A1*cvals["10"], Ztarget["10"]),
], [B1, A1], dict=True)[0]

sol_even_first_partial = sp.solve([
    sp.Eq(B2 + A2*cvals["20"], Ztarget["20"]),
    sp.Eq(B2 + A2*cvals["21"], Ztarget["21"]),
], [B2, A2], dict=True)[0]
pred_Z22_first = sp.simplify((B2 + A2*cvals["22"]).subs(sol_even_first_partial))
sol_even_first_full = sp.solve([
    sp.Eq(B2 + A2*cvals["20"], Ztarget["20"]),
    sp.Eq(B2 + A2*cvals["21"], Ztarget["21"]),
    sp.Eq(B2 + A2*cvals["22"], Ztarget["22"]),
], [B2, A2], dict=True)

sol_even_second = sp.solve([
    sp.Eq(B2 + A2*cvals["20"] + C2*dvals["20"], Ztarget["20"]),
    sp.Eq(B2 + A2*cvals["21"] + C2*dvals["21"], Ztarget["21"]),
    sp.Eq(B2 + A2*cvals["22"] + C2*dvals["22"], Ztarget["22"]),
], [B2, A2, C2], dict=True)[0]

source_sol = sp.solve([
    sp.Eq((sp.Rational(2,1)/sp.sqrt(5))*(p_iso + q_ax/sp.Integer(3)), sp.Rational(4,1)/sp.sqrt(5)),
    sp.Eq(sp.Rational(2,3)*q_ax, sp.Rational(5,4)),
], [p_iso, q_ax], dict=True)[0]

B0_val = Ztarget["0"]

# 4D-ball DtN kernel
def Lambda4_ball(l: int, z: mp.mpf) -> mp.mpf:
    nu = l + 1
    Jnu = mp.besselj(nu, z)
    Jprime = 0.5*(mp.besselj(nu-1, z) - mp.besselj(nu+1, z))
    return -1 + z * Jprime / Jnu

def dLambda4_ball(l: int, z: mp.mpf) -> mp.mpf:
    h = mp.mpf('1e-8')
    return (Lambda4_ball(l, z + h) - Lambda4_ball(l, z - h)) / (2*h)

def find_root(l: int, B: mp.mpf, guess: float) -> mp.mpf:
    f = lambda zz: Lambda4_ball(l, zz) + B
    return mp.findroot(f, guess)

B0_num = mp.mpf(sp.N(B0_val, 50))
B1_base = mp.mpf(sp.N(sol_odd[B1], 50))
A1_val = mp.mpf(sp.N(sol_odd[A1], 50))
B2_base = mp.mpf(sp.N(sol_even_second[B2], 50))
A2_val = mp.mpf(sp.N(sol_even_second[A2], 50))
C2_val = mp.mpf(sp.N(sol_even_second[C2], 50))

# exact static channel impedances
B1_perp = mp.mpf(sp.N(Ztarget["1perp"], 50))
B1_0 = mp.mpf(sp.N(Ztarget["10"], 50))
B2_0 = mp.mpf(sp.N(Ztarget["20"], 50))
B2_1 = mp.mpf(sp.N(Ztarget["21"], 50))
B2_2 = mp.mpf(sp.N(Ztarget["22"], 50))

z0_exact = find_root(0, B0_num, 0.6)
z1_base = find_root(1, B1_base, 2.55)
z1_perp_exact = find_root(1, B1_perp, 2.56)
z1_0_exact = find_root(1, B1_0, 2.53)
z2_base = find_root(2, B2_base, 4.25)
z2_0_exact = find_root(2, B2_0, 3.90)
z2_1_exact = find_root(2, B2_1, 4.03)
z2_2_exact = find_root(2, B2_2, 4.82)

c1_perp_num = mp.mpf(sp.N(cvals["1perp"], 50))
c1_0_num = mp.mpf(sp.N(cvals["10"], 50))
c2_0_num = mp.mpf(sp.N(cvals["20"], 50))
c2_1_num = mp.mpf(sp.N(cvals["21"], 50))
c2_2_num = mp.mpf(sp.N(cvals["22"], 50))
d2_0_num = mp.mpf(sp.N(dvals["20"], 50))
d2_1_num = mp.mpf(sp.N(dvals["21"], 50))
d2_2_num = mp.mpf(sp.N(dvals["22"], 50))

z1_perp_pert = z1_base - (A1_val*c1_perp_num)/dLambda4_ball(1, z1_base)
z1_0_pert = z1_base - (A1_val*c1_0_num)/dLambda4_ball(1, z1_base)

z2_0_pert = z2_base - (A2_val*c2_0_num + C2_val*d2_0_num)/dLambda4_ball(2, z2_base)
z2_1_pert = z2_base - (A2_val*c2_1_num + C2_val*d2_1_num)/dLambda4_ball(2, z2_base)
z2_2_pert = z2_base - (A2_val*c2_2_num + C2_val*d2_2_num)/dLambda4_ball(2, z2_base)

z0_small = 4/mp.sqrt(45)

def fmt_q(x):
    return str(sp.simplify(x))

def fmt_n(x):
    return mp.nstr(x, 16)

lines = []
lines.append("2PN axisymmetric Robin-wall / P2-perturbation prelim")
lines.append("")
lines.append("Channel expectation values")
for key in ["1perp", "10", "20", "21", "22"]:
    lines.append(f"  <P2>_{key}   = {fmt_q(cvals[key])}")
for key in ["1perp", "10", "20", "21", "22"]:
    lines.append(f"  <P2^2>_{key} = {fmt_q(dvals[key])}")
lines.append("")
lines.append("Odd-sector exact first-order P2 wall fit")
lines.append(f"  B1 = {fmt_q(sol_odd[B1])}")
lines.append(f"  A1 = {fmt_q(sol_odd[A1])}")
lines.append(f"  Z_1perp = {fmt_q((B1 + A1*cvals['1perp']).subs(sol_odd))}  -> Y_1perp = {fmt_q(Ytarget['1perp'])}")
lines.append(f"  Z_10    = {fmt_q((B1 + A1*cvals['10']).subs(sol_odd))}  -> Y_10    = {fmt_q(Ytarget['10'])}")
lines.append("")
lines.append("Even-sector first-order P2 wall no-go")
lines.append(f"  partial fit from (20,21): B2 = {fmt_q(sol_even_first_partial[B2])}, A2 = {fmt_q(sol_even_first_partial[A2])}")
lines.append(f"  predicted Z_22(first-order) = {fmt_q(pred_Z22_first)} but target is {fmt_q(Ztarget['22'])}")
lines.append(f"  full first-order solve exists? {bool(sol_even_first_full)}")
lines.append("")
lines.append("Even-sector exact second-order P2^2 wall fit")
lines.append(f"  B2 = {fmt_q(sol_even_second[B2])}")
lines.append(f"  A2 = {fmt_q(sol_even_second[A2])}")
lines.append(f"  C2 = {fmt_q(sol_even_second[C2])}")
for key in ["20", "21", "22"]:
    expr = sp.simplify((B2 + A2*cvals[key] + C2*dvals[key]).subs(sol_even_second))
    lines.append(f"  Z_{key} = {fmt_q(expr)}  -> Y_{key} = {fmt_q(Ytarget[key])}")
lines.append("")
lines.append("Monopole stiffness")
lines.append(f"  B0 = {fmt_q(B0_val)}  -> Y_0 = {fmt_q(Ytarget['0'])}")
lines.append("")
lines.append("Axisymmetric source profile")
lines.append(f"  p_iso = {fmt_q(source_sol[p_iso])}")
lines.append(f"  q_ax  = {fmt_q(source_sol[q_ax])}")
lines.append("")
lines.append("4D-ball DtN pole roots z = a*Omega/c_s from Lambda_l(z) + B_channel = 0")
lines.append(f"  z0 exact       = {fmt_n(z0_exact)}")
lines.append(f"  z0 small-z est = {fmt_n(z0_small)}")
lines.append(f"  z1 base        = {fmt_n(z1_base)}")
lines.append(f"  z1_perp exact  = {fmt_n(z1_perp_exact)}")
lines.append(f"  z1_perp pert   = {fmt_n(z1_perp_pert)}")
lines.append(f"  z10 exact      = {fmt_n(z1_0_exact)}")
lines.append(f"  z10 pert       = {fmt_n(z1_0_pert)}")
lines.append(f"  z2 base        = {fmt_n(z2_base)}")
lines.append(f"  z20 exact      = {fmt_n(z2_0_exact)}")
lines.append(f"  z20 pert       = {fmt_n(z2_0_pert)}")
lines.append(f"  z21 exact      = {fmt_n(z2_1_exact)}")
lines.append(f"  z21 pert       = {fmt_n(z2_1_pert)}")
lines.append(f"  z22 exact      = {fmt_n(z2_2_exact)}")
lines.append(f"  z22 pert       = {fmt_n(z2_2_pert)}")
lines.append("")
lines.append("Interpretation")
lines.append("  * finite B0 gives a finite monopole pole (Omega0 != 0)")
lines.append("  * first-order P2 splitting is already enough for the odd dipole sector")
lines.append("  * the even quadrupole support pattern needs a second-order P2^2 wall term")
lines.append("  * this is a concrete axisymmetric throat-wall PDE scaffold, not just an algebraic fit")

print("\n".join(lines))
