#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

L, Theta = sp.symbols('L Theta', positive=True, real=True)
Ks, Kq, lam, gs, gq = sp.symbols('K_s K_q lam g_s g_q', positive=True, real=True)
rc = sp.symbols('r_c', positive=True, real=True)

rho_c = sp.simplify(gs**2 / Ks)
sigma_c = sp.simplify((Ks*gq - lam*gs)**2 / (Ks*(Ks*Kq + lam**2)))

# --- F1 (R3): independent matrix-Schur reconstruction of rho_c, sigma_c. ---
# The PHYSICAL core block (notes stage097 :16-33; owner script stage114 :25-49) is the
# two-channel stiffness matrix acting on the internal core coordinates (s, q); rho_c and
# sigma_c are NOT inputs to it. Inverting it is an independent derivation primitive.
kappa0, gamma0 = sp.symbols('kappa0 gamma0', positive=True, real=True)  # bare mixed pair, reused in F2
D_sch = sp.symbols('D_sch', positive=True)
M_core = sp.Matrix([[Ks, lam], [lam, -Kq * D_sch]])   # physical core stiffness matrix
v_coup = sp.Matrix([gs, gq])                          # mouth coupling vector (g_s, g_q)
# Mouth feedback delta_Lambda_core = v^T M_core^{-1} v  (Schur elimination of (s, q)).
delta_Lambda_schur = sp.apart((v_coup.T * M_core.inv() * v_coup)[0], D_sch)
# rho_c is the D -> oo limit (mixed side-channel frozen out); sigma_c is the residual that
# the finite-D term removes at the static point D = 1 (D_W_bare(0) = 1).
rho_c_schur = sp.simplify(sp.limit(delta_Lambda_schur, D_sch, sp.oo))
sigma_c_schur = sp.simplify(rho_c_schur - delta_Lambda_schur.subs(D_sch, 1))
assert sp.simplify(rho_c - rho_c_schur) == 0, "rho_c does not match the M_core Schur residue"
assert sp.simplify(sigma_c - sigma_c_schur) == 0, "sigma_c does not match the M_core Schur residue"
print('rho_c, sigma_c reproduced from explicit two-channel core Schur complement (M_core).')

Ms = sp.simplify(L * rho_c / Theta)
Mq = sp.simplify(-L * sigma_c / Theta)

print('rho_c =', rho_c)
print('sigma_c =', sigma_c)
print('M_s =', Ms)
print('M_q =', Mq)

# Anchor M_s, M_q against the literal forms quoted in paper/stages/stage_137.tex.
# Ms and Mq above were built via the route L*rho_c/Theta and -L*sigma_c/Theta;
# the *_paper forms are constructed directly from the paper-card primitives.
# Any sign/factor error in rho_c, sigma_c, Ms, or Mq propagation will fail this.
Ms_paper = L * gs**2 / (Ks * Theta)
Mq_paper = -L * (Ks*gq - lam*gs)**2 / (Ks * (Ks*Kq + lam**2) * Theta)
assert sp.simplify(Ms - Ms_paper) == 0
assert sp.simplify(Mq - Mq_paper) == 0
print('M_s matches paper card closed form.')
print('M_q matches paper card closed form.')

# F2 (R1): full-susceptibility anchor against the matrix-Schur source (NOT X - X).
# The reduced envelope rho_c - sigma_c/(1 - kappa_c z^2 - I gamma_c z^5) must equal the
# matrix route v^T M_core^{-1} v evaluated on the bare denominator D_W_bare(z), using the
# Stage 97/114 coefficient relations (notes stage097 :55-89; stage114 :34-49). A wrong
# factor in the hand-typed sigma_c leaves a NONZERO residual against the matrix route.
z_var = sp.symbols('z_var', positive=True, real=True)
r_c = lam**2 / (Ks * Kq)
kappa_c = kappa0 / (1 + r_c)
gamma_c = gamma0 / (1 + r_c)
D_W_bare = 1 - kappa0*z_var**2 - sp.I*gamma0*z_var**5
# Independent source: the inverted physical matrix on the bare denominator.
delta_Lambda_matrix = sp.simplify(delta_Lambda_schur.subs(D_sch, D_W_bare))
# Reduced envelope built from the hand-assigned rho_c, sigma_c and Schur coefficient maps.
delta_Lambda_reduced = rho_c - sigma_c / (1 - kappa_c*z_var**2 - sp.I*gamma_c*z_var**5)
assert sp.simplify(delta_Lambda_matrix - delta_Lambda_reduced) == 0, (
    "reduced core susceptibility does not match the M_core Schur source on D_W_bare(z)"
)
print('Reduced core susceptibility matches the matrix-Schur source (full z dependence).')
# Static specialization (z -> 0) now tests DERIVED content, since both sides trace to M_core.
static_limit = sp.limit(delta_Lambda_matrix, z_var, 0)
assert sp.simplify(static_limit - (rho_c_schur - sigma_c_schur)) == 0, (
    "static core residue does not match rho_c_schur - sigma_c_schur from M_core"
)
print('Static core residue matches rho_c_schur - sigma_c_schur from M_core.')

# F3 (R2): outlet consistency with a NONZERO mixed channel (paper Checks item 1).
# Family-1 fixed-point law (notes stage137 :94-101): Pi = M_s + M_q * S_q(Pi). At a
# NONZERO susceptibility S_q the map value is Pi = M_s + M_q * S_q; the mixed contribution
# M_q * S_q must equal the matrix-Schur reconstruction -L*sigma_c_schur*S_q/Theta (M_q
# rebuilt from the inverted physical matrix of F1), pinning its sign AND its Schur factor.
Pi_var, Sq_var = sp.symbols('Pi_var S_q_var', real=True)
Pi_map = Ms + Mq * Sq_var                       # fixed-point map at a generic, NONZERO S_q
mixed_contribution = sp.simplify(Pi_map - Ms)   # isolates the M_q * S_q term (not deleted)
Mq_from_schur = -L * sigma_c_schur / Theta      # M_q rebuilt from the matrix-Schur sigma_c
assert sp.simplify(mixed_contribution - Mq_from_schur * Sq_var) == 0, (
    "M_q * S_q outlet term does not match -L*sigma_c_schur*S_q/Theta (sign/factor of M_q)"
)
print('Outlet mixed channel M_q*S_q matches the matrix-Schur reconstruction (S_q != 0).')
# Sanity: a flipped-sign M_q would NOT satisfy the above (guard is non-vacuous).
assert sp.simplify(mixed_contribution - (-Mq_from_schur) * Sq_var) != 0, (
    "outlet check is vacuous: +M_q and -M_q both pass"
)
print('Outlet consistency (paper Checks item 1) verified with nonzero S_q.')

# Check equivalence with r_c notation.
expr_rc = sp.simplify((Ks*gq - lam*gs)**2 / (Ks**2*Kq*(1 + lam**2/(Ks*Kq))))
print('sigma_c (r_c form) =', expr_rc)
assert sp.simplify(sigma_c - expr_rc) == 0

print('\nFinal explicit gain map verified.')
