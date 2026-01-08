import sympy as sp

# ----------------------------------------------------------------------
# 1. Symbols and basic objects
# ----------------------------------------------------------------------
r, theta, w = sp.symbols('r theta w', real=True)
a, L = sp.symbols('a L', positive=True)         # throat radius and bulk depth scale
rho0, eps = sp.symbols('rho0 eps', real=True)   # density scale and anisotropy parameter
G = sp.symbols('G', positive=True)              # gravitational constant (for later)

phi = sp.symbols('phi', real=True)

# Legendre P2(cos theta) = (3 cos^2 theta - 1)/2
P2 = (3*sp.cos(theta)**2 - 1) / 2

# ----------------------------------------------------------------------
# 2. 4D density ansatz and dimensional reduction over w
# ----------------------------------------------------------------------
# 4D density: Gaussian in r and w, with small angular anisotropy eps * P2
rho4 = rho0 * sp.exp(-r**2/a**2) * sp.exp(-w**2/L**2) * (1 + eps * P2)

# Integrate over w to get effective 3D density on the brane
rho3 = sp.integrate(rho4, (w, -sp.oo, sp.oo))
rho3 = sp.simplify(rho3)

print("Effective 3D density rho3(r, theta) =")
sp.pprint(rho3)
print()

# ----------------------------------------------------------------------
# 3. Total mass M (integrate rho3 over 3D space)
# ----------------------------------------------------------------------
dV = r**2 * sp.sin(theta)  # volume element in spherical coords (r, theta, phi)

M = sp.integrate(
    rho3 * dV,
    (r, 0, sp.oo),
    (theta, 0, sp.pi),
    (phi, 0, 2*sp.pi)
)

M_simplified = sp.simplify(M)
print("Total mass M =")
sp.pprint(M_simplified)
print()

# ----------------------------------------------------------------------
# 4. Quadrupole moment Q_20
#    Convention: Q_20 ~ ∫ rho3 * r^2 * P2(cos theta) d^3x
#    (Up to overall numeric factors; we're mainly interested in scaling.)
# ----------------------------------------------------------------------
Q20 = sp.integrate(
    rho3 * r**2 * P2 * dV,
    (r, 0, sp.oo),
    (theta, 0, sp.pi),
    (phi, 0, 2*sp.pi)
)

Q20_simplified = sp.simplify(Q20)
print("Quadrupole-like moment Q_20 =")
sp.pprint(Q20_simplified)
print()

# ----------------------------------------------------------------------
# 5. Ratio Q_20 / M to see the scaling
# ----------------------------------------------------------------------
ratio = sp.simplify(Q20_simplified / M_simplified)
print("Ratio Q_20 / M =")
sp.pprint(ratio)
print()

# ----------------------------------------------------------------------
# 6. Optional: Far-field potential (schematic)
#    Phi(r, theta) ~ -G M / r - G Q_20 * P2(cos theta) / r^3
#    This is just to show structure; we don't actually integrate here.
# ----------------------------------------------------------------------
Phi_far = -G * M_simplified / r - G * Q20_simplified * P2 / r**3
Phi_far_simplified = sp.simplify(Phi_far)

print("Schematic far-field potential Phi_far(r, theta) (up to overall constants):")
sp.pprint(Phi_far_simplified)
print("Note: quadrupole correction scales like eps * (a^2 / r^3).")

