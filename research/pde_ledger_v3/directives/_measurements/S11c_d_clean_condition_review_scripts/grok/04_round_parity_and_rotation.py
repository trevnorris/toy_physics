#!/usr/bin/env python3
"""H-round inversion parity of toroidal vs spheroidal vs scalar,
and H-planar rotation of an oblique wavevector about the interface normal.
"""
import sympy as sp

print("=== H-PLANAR ROTATION ABOUT INTERFACE NORMAL y ===")
# k = (k_x, 0, k_z) in (x, y, z); rotation about y by alpha
kx, kz, alpha = sp.symbols("k_x k_z alpha", real=True)
# R_y(alpha) on (x,z): x' = x cos a + z sin a, z' = -x sin a + z cos a
# wavevector transforms as a covector the same way
kxp = kx * sp.cos(alpha) + kz * sp.sin(alpha)
kzp = -kx * sp.sin(alpha) + kz * sp.cos(alpha)
# Choose alpha so k_z' = 0: tan alpha = kz/kx  (if kx != 0), or alpha = pi/2 if kx=0
# residual k_z' after alpha = atan2(kz, kx)
a_star = sp.atan2(kz, kx)
kzp_star = sp.simplify(-kx * sp.sin(a_star) + kz * sp.cos(a_star))
kxp_star = sp.simplify(kx * sp.cos(a_star) + kz * sp.sin(a_star))
print(f"alpha_* = atan2(k_z, k_x)")
print(f"k_z'(alpha_*) = {kzp_star}")
print(f"k_x'(alpha_*) = {kxp_star}")
# Numerical witnesses
for kx_n, kz_n in [(3, 4), (0, 5), (-2, 1), (1, 0)]:
    a = sp.atan2(sp.Integer(kz_n), sp.Integer(kx_n))
    kz_p = float(-kx_n * sp.sin(a) + kz_n * sp.cos(a))
    kx_p = float(kx_n * sp.cos(a) + kz_n * sp.sin(a))
    mag = float(sp.sqrt(kx_n**2 + kz_n**2))
    print(f"WITNESS k=({kx_n},{kz_n}): k'=({kx_p:.12g}, {kz_p:.12g}) |k'|={abs(kx_p):.12g} |k|={mag:.12g}")

print("\nRotation maps each single Fourier mode to class P. A superposition of")
print("distinct incidence planes is a sum of such modes; each TE amplitude is")
print("classified separately. Rotation does not map a general field to one P field.")

print("\n=== H-ROUND INVERSION PARITY ===")
# Scalar Y_lm ~ r^l homogeneous harmonic, inversion r -> -r multiplies by (-1)^l
# Toroidal polar vector: u = r × ∇Ψ, Ψ scalar of order l
# Spheroidal: u = ∇(Φ) + ... or u = f(r) Y_lm \hat{r} + ...

x, y, z = sp.symbols("x y z", real=True)
r2 = x**2 + y**2 + z**2
# l=1, m=0 scalar: Psi = z  (Y_10 ~ z/r, times r for a linear field)
Psi_10 = z
# grad Psi
gPsi = (sp.diff(Psi_10, x), sp.diff(Psi_10, y), sp.diff(Psi_10, z))
# u_tor = r × grad Psi
u_tor = (
    y * gPsi[2] - z * gPsi[1],
    z * gPsi[0] - x * gPsi[2],
    x * gPsi[1] - y * gPsi[0],
)
u_tor = tuple(sp.simplify(c) for c in u_tor)
print(f"l=1 m=0 toroidal u = r × ∇z = {u_tor}")
# inversion I: (x,y,z) -> (-x,-y,-z). A polar vector transforms as
# u'(r') = - u(r)  when r' = -r, i.e. u'(-r) = - u(r)
# Evaluate u at -r and compare to -u(r)
u_at_minus = tuple(c.subs({x: -x, y: -y, z: -z}) for c in u_tor)
minus_u = tuple(-c for c in u_tor)
# For a field of inversion parity eta, u(-r) = eta * (polar transformation)
# polar vector of parity +1 (like a gradient of even scalar): u(-r) = -u(r)
# polar vector of parity -1: u(-r) = +u(r)
# Standard: inversion parity eta means u_i(-r) = eta * (+/- depending on polar)
# For polar vectors, the intrinsic inversion of components:
#   u_i(-x) = P * u_i(x)  with P = (-1)^{l+1} for toroidal, (-1)^l for spheroidal
# in the sense of the spherical-harmonic factor, AFTER the polar-vector minus.
#
# Direct: u_tor(r) = (y, -x, 0) for Psi=z? Let's see:
print(f"u_tor(r)      = {u_tor}")
print(f"u_tor(-r)     = {u_at_minus}")
print(f"-u_tor(r)     = {minus_u}")
print(f"u(-r) + u(r)  = {tuple(sp.simplify(a+b) for a,b in zip(u_at_minus, u_tor))}")
print(f"u(-r) - u(r)  = {tuple(sp.simplify(a-b) for a,b in zip(u_at_minus, u_tor))}")
print(f"u(-r) + (-u)  match polar-even? {u_at_minus == minus_u}")
print(f"u(-r) == u(r) axial-like? {u_at_minus == u_tor}")

# Scalar l=1: S(-r) = -S(r), parity (-1)^l = -1
S10 = z
print(f"\nscalar l=1: S(r)=z, S(-r)={S10.subs({z:-z})}, (-1)^l S = {-S10}")

# Spheroidal l=1: u = ∇S = (0,0,1), a polar vector
u_sph = (sp.diff(S10, x), sp.diff(S10, y), sp.diff(S10, z))
u_sph_minus = tuple(c.subs({x: -x, y: -y, z: -z}) for c in u_sph)
print(f"spheroidal l=1 u=∇z = {u_sph}")
print(f"u_sph(-r) = {u_sph_minus}")
print(f"-u_sph(r) = {tuple(-c for c in u_sph)}")
print(f"spheroidal: u(-r) == -u(r)? {u_sph_minus == tuple(-c for c in u_sph)}")

# l=2 scalar ~ x^2-y^2, parity (+1)
S20 = x**2 - y**2
print(f"\nscalar l=2: S(-r)={S20.subs({x:-x,y:-y})}, S={S20}, equal? {sp.expand(S20.subs({x:-x,y:-y})-S20)==0}")

# toroidal l=2: Psi = x^2 - y^2, u = r × ∇Psi
g2 = (sp.diff(S20, x), sp.diff(S20, y), sp.diff(S20, z))
u_tor2 = (
    sp.simplify(y * g2[2] - z * g2[1]),
    sp.simplify(z * g2[0] - x * g2[2]),
    sp.simplify(x * g2[1] - y * g2[0]),
)
u_tor2_minus = tuple(sp.simplify(c.subs({x: -x, y: -y, z: -z})) for c in u_tor2)
print(f"l=2 toroidal u = {u_tor2}")
print(f"u(-r) = {u_tor2_minus}")
print(f"-u(r) = {tuple(sp.simplify(-c) for c in u_tor2)}")
print(f"l=2 toroidal u(-r)==-u(r)? {u_tor2_minus == tuple(sp.simplify(-c) for c in u_tor2)}")
print(f"l=2 toroidal u(-r)==+u(r)? {u_tor2_minus == tuple(sp.simplify(c) for c in u_tor2)}")

print("\nConvention used by A: toroidal inversion parity (-1)^{l+1},")
print("scalar/spheroidal (-1)^l.")
print("Computed (polar-vector component test u_i(-r) vs ± u_i(r)):")
print("  l=1 toroidal u=(y,-x,0): u(-r)=(-y,x,0)= -u(r)  => extra minus,")
print("      matching polar-vector * (-1)^{l+1}=(-1)^2=+1 times the polar minus")
print("      gives u(-r)=-u(r). Scalar l=1 has S(-r)=-S(r)=(-1)^l S.")
print("  Opposite parities at the same l, so they do not mix under inversion-even ops.")

# Same-parity different-l: toroidal l and scalar l+1
print("\n=== SAME INVERSION PARITY, DIFFERENT l ===")
print("toroidal l parity (-1)^{l+1} equals scalar l' parity (-1)^{l'} when l'=l+1.")
print("SO(3) conservation of l still forbids that mixing on an O(3) background.")
print("A's 'therefore' at each (l,m) uses both inversion and rotation.")
