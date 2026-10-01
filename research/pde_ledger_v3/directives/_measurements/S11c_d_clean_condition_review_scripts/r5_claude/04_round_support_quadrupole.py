#!/usr/bin/env python3
"""H-round on the throat the twist mode itself holds open (unit-sphere angular algebra, exact).

Sphere integrals use the exact monomial moments  int x^a y^b z^c dOmega
  = 2 G((a+1)/2) G((b+1)/2) G((c+1)/2) / G((a+b+c+3)/2)   (all exponents even; else 0).
l = 2 harmonics are the pure-l=2 polynomials (x+iy)^2, z(x+iy), 2z^2-x^2-y^2, z(x-iy), (x-iy)^2 (unnormalised).

(a) Coherent l=1 toroidal pattern u = Re[(a.T) e^{-iwt}], T_i = e_i x r, a = p + i q.
    <|u|^2> = (1/2) sum_ij B_ij T_i.T_j,  B = p p^T + q q^T = Re(a a^dagger).  Print its l=2 projections,
    their relation to the traceless part of B, det B, and controls.
(b) Thickness background W = W0(r)(1 + delta P2(z)); tangential gradient g_t (the u.g form of the material-advected
    anchoring).  Print C(m', m) = int (g_t . T_(1m')) conj(Y_2m) dOmega for T_z, T_x + iT_y, T_x - iT_y, with delta
    live and the round control delta = 0, and the l=0 / l=1 projections.  Print inversion and y-mirror parities.
"""
import sympy as sp

x, y, z = sp.symbols("x y z", real=True)
delta = sp.symbols("delta", real=True)
r = sp.Matrix([x, y, z])
E = [sp.Matrix([1, 0, 0]), sp.Matrix([0, 1, 0]), sp.Matrix([0, 0, 1])]
T = [e.cross(r) for e in E]

def mono(a, b, c):
    if a % 2 or b % 2 or c % 2:
        return sp.Integer(0)
    G = sp.gamma
    return sp.simplify(2 * G(sp.Rational(a + 1, 2)) * G(sp.Rational(b + 1, 2)) * G(sp.Rational(c + 1, 2))
                       / G(sp.Rational(a + b + c + 3, 2)))

def sph_int(f):
    poly = sp.Poly(sp.expand(f), x, y, z)
    return sp.simplify(sum(coef * mono(*mon) for mon, coef in poly.terms()))

Y2 = {2: (x + sp.I * y) ** 2, 1: z * (x + sp.I * y), 0: 2 * z**2 - x**2 - y**2,
      -1: z * (x - sp.I * y), -2: (x - sp.I * y) ** 2}
Y1 = {1: x + sp.I * y, 0: z, -1: x - sp.I * y}
cj = lambda e: sp.conjugate(e).subs({sp.conjugate(x): x, sp.conjugate(y): y, sp.conjugate(z): z})

# ---- (a)
p = sp.Matrix(sp.symbols("p1 p2 p3", real=True))
q = sp.Matrix(sp.symbols("q1 q2 q3", real=True))
B = p * p.T + q * q.T
def avg_u2(Bm):
    return sp.Rational(1, 2) * sum(Bm[i, j] * T[i].dot(T[j]) for i in range(3) for j in range(3))
def l2_proj(Bm):
    f = avg_u2(Bm)
    return {m: sp.factor(sph_int(f * cj(Y2[m]))) for m in Y2}
proj = l2_proj(B)
print("A_L2_PROJECTIONS_GENERAL", proj)
dev = B - B.trace() / 3 * sp.eye(3)
print("A_L2_PROJECTION_M0_OVER_TRACELESS_B33", sp.simplify(proj[0] / dev[2, 2]))
print("A_DET_B", sp.simplify(B.det()))
l1, l2_ = sp.symbols("lambda1 lambda2", nonnegative=True)
tt = l1 + l2_
fro = sp.expand((l1 - tt / 3) ** 2 + (l2_ - tt / 3) ** 2 + (tt / 3) ** 2)
print("A_TRACELESS_FROBENIUS_EIGENVALUES_l1_l2_0", fro, " minus (trace/3)^2 =", sp.factor(fro - tt**2 / 9))
print("CONTROL_ISOTROPIC_INCOHERENT_L2", l2_proj(sp.eye(3)))
print("CONTROL_LINEAR_TZ_L2", l2_proj(sp.diag(0, 0, 1)))
print("CIRCULAR_TX_PLUS_ITY_L2", l2_proj(sp.diag(1, 1, 0)))
print("A_L0_PROJECTION_GENERAL", sp.factor(sph_int(avg_u2(B))))

# ---- (b)
f = (3 * z**2 - (x**2 + y**2 + z**2)) / 2          # r^2 P2(cos theta)
gradf = sp.Matrix([sp.diff(f, v) for v in (x, y, z)])
g_t = delta * (gradf - r * r.dot(gradf))            # tangential part on the unit sphere (x^2+y^2+z^2 = 1 there)
patterns = {"T_z": T[2], "T_x_plus_iT_y": T[0] + sp.I * T[1], "T_x_minus_iT_y": T[0] - sp.I * T[1]}
for name, field in patterns.items():
    s = sp.expand(g_t.dot(field))
    row = {m: sp.factor(sph_int(s * cj(Y2[m]))) for m in Y2}
    print("B_COUPLING_L2", name, row)
    print("CONTROL_ROUND_DELTA0", name, {m: v.subs(delta, 0) for m, v in row.items()})
    print("B_COUPLING_L0_L1", name, sph_int(s), {m: sph_int(s * cj(Y1[m])) for m in Y1})

# ---- parities (polar vector: P[u](r) = M u(M^{-1} r); scalar: f(M^{-1} r))
def vec_par(u, M):
    img = M * u.subs({x: (M * r)[0], y: (M * r)[1], z: (M * r)[2]}, simultaneous=True)
    if sp.simplify(img - u) == sp.zeros(3, 1): return 1
    if sp.simplify(img + u) == sp.zeros(3, 1): return -1
    return "MIXED"
def sca_par(fn, M):
    img = fn.subs({x: (M * r)[0], y: (M * r)[1], z: (M * r)[2]}, simultaneous=True)
    if sp.simplify(img - fn) == 0: return 1
    if sp.simplify(img + fn) == 0: return -1
    return "MIXED"
Inv, My = -sp.eye(3), sp.diag(1, -1, 1)
print("PARITY_INVERSION", "T_z", vec_par(T[2], Inv), "T_x+iT_y", vec_par(T[0] + sp.I * T[1], Inv),
      "Y_21", sca_par(Y2[1], Inv), "Y_20", sca_par(Y2[0], Inv), "deformation_P2", sca_par(f, Inv))
print("PARITY_MIRROR_Y", "T_z", vec_par(T[2], My), "Y_20", sca_par(Y2[0], My),
      "deformation_P2", sca_par(f, My), "T_x+iT_y", vec_par(T[0] + sp.I * T[1], My), "Y_21", sca_par(Y2[1], My))
