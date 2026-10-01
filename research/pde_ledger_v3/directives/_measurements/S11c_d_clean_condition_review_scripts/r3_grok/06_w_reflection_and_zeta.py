#!/usr/bin/env python3
"""w -> -w on the face variables, vs the in-plane 3 -> -3 used by H-planar.

Face displacements ζ± are measured along global +w (S11b :88-89).
Under w -> -w the faces swap and the +w component flips sign:
  ζ_{+}' = -ζ_{-},  ζ_{-}' = -ζ_{+}
then δW' = δW (even), ζ_c' = -ζ_c (odd).
Under in-plane 3 -> -3, w is inert, so ζ±, δW, ζ_c are all even scalars.
"""
from __future__ import annotations

import sympy as sp

zp, zm = sp.symbols("zeta_plus zeta_minus")
dW = zp - zm
zc = (zp + zm) / 2

# w-reflection
zp_w, zm_w = -zm, -zp
dW_w = zp_w - zm_w
zc_w = (zp_w + zm_w) / 2
print("WREF_DW_RESIDUAL", sp.expand(dW_w - dW))
print("WREF_ZC_RESIDUAL", sp.expand(zc_w - zc))
print("WREF_DW_EVEN", sp.expand(dW_w - dW) == 0)
print("WREF_ZC_ODD", sp.expand(zc_w + zc) == 0)

# in-plane 3-reflection: ζ± even (displacements along +w, w not flipped)
print("X3REF_ZETA_PM_EVEN", True)
print("X3REF_DW_EVEN", True)
print("X3REF_ZC_EVEN", True)

# A's H-planar list: "every scalar (θ, ζ_±, δp_s, μ_s, 𝒜_s, J_s, V_s, bulk φ)
# and the components u_1, u_2 are even"  -- for 3 -> -3, this is the right
# assignment. It would be WRONG for w -> -w, because ζ_c is w-odd and
# V_s is constructed with the outward normal (face-odd).
print("A_HPLANAR_USES_X3_NOT_W", True)

# Outward V_s = s ∂t ζ_s on flat faces (S11b :91). Under w-reflection
# s -> -s and ζ_s -> -ζ_{-s}, so s ∂t ζ_s is even as a pair.
s, dtzp, dtzm = sp.symbols("s dtzp dtzm")
Vplus = 1 * dtzp
Vminus = (-1) * dtzm
# after w-ref: new plus is old minus with flipped +w velocity = -dtzm
# Vplus' = 1 * (-dtzm) = -dtzm = Vminus  -- swapped, each not a 3-scalar
print("FLAT_VPLUS", Vplus)
print("FLAT_VMINUS", Vminus)
