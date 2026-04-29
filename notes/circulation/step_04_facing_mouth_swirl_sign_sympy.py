#!/usr/bin/env python3
"""
Step 4 SymPy audit: translating global current-loop attraction into local swirl
labels for parallel-normal versus facing-mouth geometry.

The orientation signs are derived from oriented loop integrals, not assigned by
hand. Positive local swirl is the right-handed loop orientation about the local
mouth normal.

Run:
    python step_04_facing_mouth_swirl_sign_sympy.py
"""
import sympy as sp

K, d, R = sp.symbols("K d R", positive=True)
Gamma0 = sp.symbols("Gamma0", positive=True)
n1, n2 = sp.symbols("n1 n2", real=True)
phi = sp.symbols("phi", real=True)

ex = sp.Matrix([1, 0, 0])
ey = sp.Matrix([0, 1, 0])
ez = sp.Matrix([0, 0, 1])
dhat = ez


def oriented_loop(e_hat, f_hat):
    """Return circulation per unit local label and area vector."""
    r = R * (sp.cos(phi) * e_hat + sp.sin(phi) * f_hat)
    dl_dphi = sp.diff(r, phi)
    tangent = sp.simplify(dl_dphi / R)
    v_phi = sp.simplify(Gamma0 / (2 * sp.pi * R) * tangent)
    circulation_per_label = sp.simplify(
        sp.integrate(v_phi.dot(dl_dphi), (phi, 0, 2 * sp.pi))
    )
    area_vector = sp.simplify(
        sp.integrate(r.cross(dl_dphi) / 2, (phi, 0, 2 * sp.pi))
    )
    sigma = sp.simplify(area_vector.dot(dhat) / (sp.pi * R**2))
    return circulation_per_label, area_vector, sigma


# Local right-handed bases: e_hat x f_hat = local normal.
gamma_plus, area_plus, sigma_plus = oriented_loop(ex, ey)       # normal +z
gamma_minus, area_minus, sigma_minus = oriented_loop(ex, -ey)   # normal -z

# Leading fixed-current force. F_d < 0 means attraction.
sigma1, sigma2 = sp.symbols("sigma1 sigma2", real=True)
F = -K * sigma1 * sigma2 * n1 * n2 / d**4

# Parallel-local-normal geometry: both positive local swirls have +z area vector.
F_parallel_same = sp.simplify(F.subs({sigma1: sigma_plus, sigma2: sigma_plus, n1: 1, n2: 1}))
F_parallel_opposite = sp.simplify(F.subs({sigma1: sigma_plus, sigma2: sigma_plus, n1: 1, n2: -1}))

# Facing-mouth geometry: mouth 1 normal +z and mouth 2 normal -z.
F_facing_same_local = sp.simplify(F.subs({sigma1: sigma_plus, sigma2: sigma_minus, n1: 1, n2: 1}))
F_facing_opposite_local = sp.simplify(F.subs({sigma1: sigma_plus, sigma2: sigma_minus, n1: 1, n2: -1}))

assert gamma_plus == Gamma0
assert gamma_minus == Gamma0
assert area_plus == sp.pi * R**2 * ez
assert area_minus == -sp.pi * R**2 * ez
assert sigma_plus == 1
assert sigma_minus == -1
assert F_parallel_same == -K/d**4
assert F_parallel_opposite == K/d**4
assert F_facing_same_local == K/d**4
assert F_facing_opposite_local == -K/d**4

print("STEP 4: LOCAL SWIRL VERSUS GLOBAL CURRENT SIGN")
print("Positive local swirl about +z:")
print("  circulation integral ∮v_phi.dl per unit label =", gamma_plus)
print("  area/current-axis vector =", area_plus.T)
print("  sigma = area.dhat/(pi R^2) =", sigma_plus)
print("Positive local swirl about -z:")
print("  circulation integral ∮v_phi.dl per unit label =", gamma_minus)
print("  area/current-axis vector =", area_minus.T)
print("  sigma = area.dhat/(pi R^2) =", sigma_minus)
print()
print("Derived leading force: F_d =", F)
print("F_d < 0 means attraction.")
print()
print("Parallel-local-normal geometry:")
print("  same local n:      ", F_parallel_same, "  attraction")
print("  opposite local n:  ", F_parallel_opposite, "  repulsion")
print()
print("Facing-mouth geometry (local normals opposite):")
print("  same local n:      ", F_facing_same_local, "  repulsion")
print("  opposite local n:  ", F_facing_opposite_local, "  attraction")
print()
print("Verdict: the facing-mouth sigma product follows from oriented loop geometry.")
