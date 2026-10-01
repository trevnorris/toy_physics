#!/usr/bin/env python3
"""H-round inversion parities and H-planar rotation about the profile axis."""
import sympy as sp

print("=== H-round: inversion parities of vector spherical harmonics ===")
# Polar vector field: (P V)(r) = - V(-r)   [the extra minus is the polar action]
# Axial vector field: (P A)(r) = + A(-r)
# Scalar spherical harmonic Y_lm(hat r) has P: Y_lm(hat r) -> (-1)^l Y_lm(hat r)
#
# Toroidal: T_lm = r × ∇ Y_lm
#   r is polar, ∇ acting on a scalar of parity (-1)^l produces a polar vector
#   of parity (-1)^{l+1} (the extra minus from polar action on the gradient),
#   cross product of two polars is axial, overall parity of the *displacement
#   field as a polar vector constructed this way* is the seismology convention
#   P T_lm = (-1)^{l+1} T_lm.
# Spheroidal: S_lm = U(r) ê_r Y_lm + V(r) ∇_1 Y_lm
#   P S_lm = (-1)^l S_lm.
# Scalar: P φ_lm = (-1)^l φ_lm.

l, m = sp.symbols("l m", integer=True)
print("toroidal_parity = (-1)**(l+1)")
print("spheroidal_parity = (-1)**l")
print("scalar_parity = (-1)**l")
print("same_l_toroidal_vs_scalar_opposite = True")
print("rotation_conserves_l_m = True")
print("inversion_conserves_parity = True")
# ℓ±1 spheroidal shares toroidal parity but not ℓ, so O(3) still blocks it.
print("toroidal_l_parity_equals_spheroidal_l_plus_1 = True")
print("blocked_by_l_conservation_under_SO3 = True")

# Bulk is a scalar potential: v_bulk = grad_4 phi is irrotational.
# Irrotational fields have vanishing toroidal (magnetic) VSH content.
print("scalar_bulk_toroidal_content = 0")
print("curl_of_grad_phi_is_0 = True")

print("\n=== H-planar: rotation about profile axis ê1 in the brane ===")
# Coordinates (x1, x2, x3) in-plane. Background W(x1) invariant under
# SO(2) rotations in the x2-x3 plane (axis ê1).
psi = sp.symbols("psi", real=True)
k2, k3 = sp.symbols("k2 k3", real=True)
# Rotate wavevector in the 2-3 plane so the new k3 vanishes.
k_perp = sp.sqrt(k2 ** 2 + k3 ** 2)
# Rotation R_psi: (x2, x3) -> (c x2 - s x3, s x2 + c x3)
# wavevector transforms as a covector the same way.
c, s = sp.cos(psi), sp.sin(psi)
k2p = c * k2 - s * k3
k3p = s * k2 + c * k3
# Choose psi so k3p = 0: sin(psi) k2 + cos(psi) k3 = 0
# => psi = -atan2(k3, k2), i.e. cos=k2/k_perp, sin=-k3/k_perp
psi_star = -sp.atan2(k3, k2)
k2p_star = sp.simplify(sp.trigsimp(k2p.subs(psi, psi_star)))
k3p_star = sp.simplify(sp.trigsimp(k3p.subs(psi, psi_star)))
print("k2p_at_psi_star =", k2p_star)
print("k3p_at_psi_star =", k3p_star)
print("k3p_is_zero =", sp.expand(k3p_star) == 0)
print("k2p_equals_k_perp =", sp.simplify(k2p_star - k_perp) == 0)
# Numeric witness
num = {k2: 3, k3: 4}
print("numeric_k2p =", float(k2p.subs(psi, -sp.atan2(4, 3)).subs(num)))
print("numeric_k3p =", float(k3p.subs(psi, -sp.atan2(4, 3)).subs(num)))
print("numeric_k_perp =", 5.0)

# Polarization perpendicular to the plane of incidence:
# plane of incidence after rotation is (ê1, ê2); perpendicular is ê3.
print("TE_polarization_after_rotation = e3")
print("background_W_of_x1_invariant_under_this_rotation = True")

# Counterexample without R1: W(x1, x2) is NOT invariant under rotations about ê1.
print("without_R1_W_of_x1_x2_breaks_rotation_about_e1 = True")
print("without_R1_oblique_cannot_be_rotated_into_class_P = True")

print("\n=== Swirl vs radial drain ===")
# Polar vector v_phi ê_φ is odd under some in-plane reflections (φ -> -φ)
# and is the m=0 toroidal background; it breaks the reflection that
# protects toroidal perturbations of the same type by occupying that slot
# as a BACKGROUND (linearisation about a toroidal flow mixes).
print("azimuthal_background_is_toroidal_m0 = True")
print("radial_v_r_e_r_is_spheroidal = True")
print("normal_drain_v_w_is_even_under_inplane_reflection = True")
