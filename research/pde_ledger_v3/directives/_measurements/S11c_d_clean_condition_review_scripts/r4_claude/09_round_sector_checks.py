#!/usr/bin/env python3
"""Leg script 09: H-round bookkeeping, computed on explicit Cartesian fields.

(a) Inversion parity P u(x) = -u(-x) (vector), P f(x) = f(-x) (scalar) of the toroidal field
    T_lm = f(r) r x grad(r^l Y_lm) and the spheroidal fields f(r) grad(r^l Y_lm),
    f(r) x Y_lm-type radial field, and the scalar f(r) r^l Y_lm, for l = 0..3 (solid harmonics).
(b) Helicity H = int u.curl u for a pure toroidal field, for a circular (quadrature) toroidal
    pair at time t, and for a toroidal + c * spheroidal mixture (Gaussian radial weight).
(c) Spin-type angular momentum int u x d_t u for the quadrature toroidal pair.
(d) Swirl background v0 = Omega z_hat x r: its inversion parity, and the Coriolis image
    Omega z_hat x T_10 decomposed: its curl-free (spheroidal) part and its l=0 scalar content
    int (div of the image) * r^0 (projection onto the breathing scalar) with Gaussian weight.
Prints computed objects only.
"""
import sympy as sp

x, y, z, t, om, c, Om = sp.symbols("x y z t omega c Omega", real=True)
X = (x, y, z)
r2 = x**2 + y**2 + z**2
g = sp.exp(-r2)

def grad(f):
    return [sp.diff(f, q) for q in X]

def curl(u):
    return [sp.diff(u[2], y) - sp.diff(u[1], z), sp.diff(u[0], z) - sp.diff(u[2], x), sp.diff(u[1], x) - sp.diff(u[0], y)]

def div(u):
    return sum(sp.diff(u[i], X[i]) for i in range(3))

def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]

def inv_vec(u):
    return [sp.expand(-ui.subs({x: -x, y: -y, z: -z}, simultaneous=True)) for ui in u]

def inv_sca(f):
    return sp.expand(f.subs({x: -x, y: -y, z: -z}, simultaneous=True))

def parity_vec(u):
    u = [sp.expand(ui) for ui in u]
    if all(ui == 0 for ui in u):
        return "ZERO_FIELD"
    pu = inv_vec(u)
    if all(sp.simplify(pu[i] - u[i]) == 0 for i in range(3)):
        return 1
    if all(sp.simplify(pu[i] + u[i]) == 0 for i in range(3)):
        return -1
    return "MIXED"

def parity_sca(f):
    f = sp.expand(f)
    pf = inv_sca(f)
    return 1 if sp.simplify(pf - f) == 0 else (-1 if sp.simplify(pf + f) == 0 else "MIXED")

# real solid harmonics r^l Y_lm (one representative m per l plus m variety at l=2)
solid = {
    0: [sp.Integer(1)],
    1: [z, x],
    2: [2*z**2 - x**2 - y**2, x*z, x*y],
    3: [z*(2*z**2 - 3*x**2 - 3*y**2), x*(4*z**2 - x**2 - y**2)],
}
rvec = [x, y, z]
for l, harms in solid.items():
    for Y in harms:
        T = [g * comp for comp in cross(rvec, grad(Y))]
        S_grad = [g * comp for comp in grad(Y)]
        S_rad = [g * Y * comp for comp in rvec]
        print("ROUND l", l, "Y", Y, "TOROIDAL_PARITY", parity_vec(T), "DIV_T", sp.simplify(div(T)),
              "SPHEROIDAL_GRAD_PARITY", parity_vec(S_grad), "SPHEROIDAL_RADIAL_PARITY", parity_vec(S_rad),
              "SCALAR_PARITY", parity_sca(g * Y))

def integrate_all(expr):
    expr = sp.expand(expr)
    return sp.simplify(sp.integrate(expr, (x, -sp.oo, sp.oo), (y, -sp.oo, sp.oo), (z, -sp.oo, sp.oo)))

# (b),(c): l = 1 toroidal fields (rotations about x and y with Gaussian weight)
Tx = [g * comp for comp in cross([1, 0, 0], rvec)]   # = g*(0,-z,y)
Ty = [g * comp for comp in cross([0, 1, 0], rvec)]   # = g*(z,0,-x)
u_quad = [sp.cos(om*t)*Tx[i] + sp.sin(om*t)*Ty[i] for i in range(3)]
H_pure = integrate_all(sum(Tx[i]*curl(Tx)[i] for i in range(3)))
H_quad = integrate_all(sum(u_quad[i]*curl(u_quad)[i] for i in range(3)))
Sx = curl(Tx)  # l=1 poloidal (spheroidal, opposite parity), divergence-free
u_mix = [Tx[i] + c*Sx[i] for i in range(3)]
H_mix = integrate_all(sum(u_mix[i]*curl(u_mix)[i] for i in range(3)))
print("HELICITY_PURE_TOROIDAL", H_pure)
print("HELICITY_QUADRATURE_TOROIDAL_PAIR", H_quad)
print("HELICITY_TOROIDAL_PLUS_c_POLOIDAL", H_mix)
print("PARITY_OF_POLOIDAL_PART", parity_vec(Sx), "PARITY_OF_TOROIDAL_PART", parity_vec(Tx))
dut = [sp.diff(ui, t) for ui in u_quad]
Lspin = [integrate_all(comp) for comp in cross(u_quad, dut)]
print("SPIN_ANGULAR_MOMENTUM_QUADRATURE_PAIR", [sp.simplify(Li) for Li in Lspin])

# (d) swirl
v0 = cross([0, 0, Om], rvec)
print("SWIRL_INVERSION_PARITY", parity_vec(v0))
T10 = [g * comp for comp in cross([0, 0, 1], rvec)]   # toroidal l=1 m=0
cor = cross([0, 0, Om], T10)                           # Coriolis image Omega z x T10
print("CORIOLIS_IMAGE", [sp.factor(ci) for ci in cor], "PARITY", parity_vec(cor))
print("CORIOLIS_IMAGE_CURL", [sp.simplify(ci) for ci in curl(cor)])
print("CORIOLIS_IMAGE_DIV", sp.simplify(div(cor)))
print("CORIOLIS_IMAGE_BREATHING_PROJECTION int div*exp(-r^2)", integrate_all(div(cor) * g))
print("CORIOLIS_IMAGE_L2_PROJECTION int div*(2z^2-x^2-y^2)*exp(-r^2)", integrate_all(div(cor) * (2*z**2 - x**2 - y**2) * g))
print("CORIOLIS_IMAGE_ONTO_T10 int cor.T10", integrate_all(sum(cor[i]*T10[i] for i in range(3))))
