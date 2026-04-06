
import sympy as s

# -----------------------------------------------------------------------------
# 2PN tensor-wake channel decomposition prototype (SymPy)
# -----------------------------------------------------------------------------
# Purpose:
#   1. Start from the solved 2PN comparable-mass ADM lift.
#   2. Show that the quartic cross sector is not an arbitrary 7-term fit:
#        - two terms are exactly a universal 1/2 leg-dressing of the frozen
#          1PN vector wake,
#        - the remaining five terms are exactly a rank-2 tensor-wake overlap
#          in a minimal projector-channel basis.
#   3. Show that the G^2/r^2 velocity sector likewise decomposes into:
#        - a purely parallel local-potential dressing of the 1PN vector wake,
#        - plus diagonal transverse / longitudinal tensor-potential channels.
#   4. Export a fully constructive pairwise 2PN cross module and its first
#      3-body predictions in a local-potential environment.
#
# Notes:
#   - This is still algebraic / channel-level, not yet a full inner PDE.
#   - But it sharpens the target for tensor_wake_2pn_rebuild.wl enormously:
#     the solved 2PN cross block already factorizes into symmetry-respecting
#     vector and rank-2 tensor channels with small rational strengths.
# -----------------------------------------------------------------------------

s.init_printing(use_unicode=True)

def section(name: str) -> None:
    print(f"\n=== {name} ===")

def show(label: str, expr) -> None:
    print(f"{label}:\n{s.simplify(expr)}\n")

# symbols
G, c = s.symbols('G c', positive=True)
mA, mB, mC = s.symbols('mA mB mC', positive=True)
rAB, rAC, rBC = s.symbols('rAB rAC rBC', positive=True)

# pair invariants for pair AB (n = n_AB)
vA2, vB2 = s.symbols('vA2 vB2')
vAB = s.symbols('vAB')      # v_A . v_B
vAn, vBn = s.symbols('vAn vBn')  # radial projections onto n_AB

# 1PN wake coefficients (frozen from vector_wake_rebuild.wl)
Cpar = -s.Rational(7, 2)
CL = -s.Rational(1, 2)

# solved 2PN cross coefficients (from the ADM lift)
q1 = -s.Rational(7, 4)
q2 = -s.Rational(1, 4)
q3 = s.Rational(11, 8)
q4 = s.Rational(1, 4)
q5 = -s.Rational(5, 8)
q6 = s.Rational(3, 2)
q7 = s.Rational(3, 8)

t2 = s.Rational(11, 8)
t3 = -s.Rational(15, 4)
t6 = s.Rational(15, 8)
s1 = s.Rational(5, 4)

section("TARGET / solved added 2PN cross block from the ADM lift")
L2_cross_target = s.expand(
    G * mA * mB / (c**4 * rAB) * (
        q1 * vAB * (vA2 + vB2)
        + q2 * vAn * vBn * (vA2 + vB2)
        + q3 * vA2 * vB2
        + q4 * vAB**2
        + q5 * (vA2 * vBn**2 + vB2 * vAn**2)
        + q6 * vAB * vAn * vBn
        + q7 * vAn**2 * vBn**2
    )
    + G**2 * mA * mB / (c**4 * rAB**2) * (
        t2 * (mA * vA2 + mB * vB2)
        + t3 * (mA + mB) * vAB
        + t6 * (mA * vAn**2 + mB * vBn**2)
    )
    + s1 * G**3 * mA**2 * mB**2 / (c**4 * rAB**3)
)
show("L2_cross_target", L2_cross_target)

section("VECTOR WAKE / universal 1/2 leg-dressing fixes q1 and q2 exactly")
sigma = s.Rational(1, 2)
L2_vec_kin = s.expand(
    G * mA * mB / (c**4 * rAB)
    * sigma * (vA2 + vB2) * (Cpar * vAB + CL * vAn * vBn)
)
show("L2_vec_kin", L2_vec_kin)
show("Recovered q1 from sigma*Cpar", s.simplify(sigma * Cpar))
show("Recovered q2 from sigma*CL", s.simplify(sigma * CL))

quartic_target = s.expand(
    G * mA * mB / (c**4 * rAB) * (
        q1 * vAB * (vA2 + vB2)
        + q2 * vAn * vBn * (vA2 + vB2)
        + q3 * vA2 * vB2
        + q4 * vAB**2
        + q5 * (vA2 * vBn**2 + vB2 * vAn**2)
        + q6 * vAB * vAn * vBn
        + q7 * vAn**2 * vBn**2
    )
)
quartic_residual = s.expand(quartic_target - L2_vec_kin)
show("quartic_residual after removing vector leg-dressing", quartic_residual)

section("TENSOR WAKE / define minimal projector-channel basis")
# transverse/longitudinal scalar channels
TA = vA2 - vAn**2
TB = vB2 - vBn**2
LA = vAn**2
LB = vBn**2

# transverse-shear channel contraction (2D transverse plane)
SAB = s.expand((vAB - vAn * vBn)**2 - s.Rational(1, 2) * TA * TB)

# mixed transverse-longitudinal tensor contraction
MAB = s.expand(2 * (vAB - vAn * vBn) * vAn * vBn)

# scalar-sector basis
TLmix = s.expand(TA * LB + TB * LA)

show("TA = transverse trace on A", TA)
show("TB = transverse trace on B", TB)
show("SAB = transverse shear contraction", SAB)
show("MAB = mixed transverse-longitudinal contraction", MAB)
show("TLmix = scalar transverse/longitudinal mixing", TLmix)
show("LA*LB = pure longitudinal scalar contraction", s.expand(LA * LB))

# solve residual quartic block in the tensor basis
kTT, kS, kM, kTL, kLL = s.symbols('kTT kS kM kTL kLL')
quartic_tensor_ansatz = s.expand(
    G * mA * mB / (c**4 * rAB) * (
        kTT * TA * TB
        + kS * SAB
        + kM * MAB
        + kTL * TLmix
        + kLL * LA * LB
    )
)

mons_quartic = [
    vA2 * vB2,
    vAB**2,
    vA2 * vBn**2,
    vB2 * vAn**2,
    vAB * vAn * vBn,
    vAn**2 * vBn**2,
]
quartic_sol = s.solve(
    [s.expand(quartic_tensor_ansatz - quartic_residual).coeff(m) for m in mons_quartic],
    [kTT, kS, kM, kTL, kLL],
    dict=True,
)
quartic_sol = quartic_sol[0]
show("Quartic tensor-channel solution", quartic_sol)

L2_tensor_quartic = s.expand(quartic_tensor_ansatz.subs(quartic_sol))
show("L2_tensor_quartic", L2_tensor_quartic)
show("Quartic tensor residual check", s.expand(L2_tensor_quartic - quartic_residual))

# rank/minimality sanity check
coeff_matrix = s.Matrix([
    [s.expand(TA * TB).coeff(m) for m in mons_quartic],
    [s.expand(SAB).coeff(m) for m in mons_quartic],
    [s.expand(MAB).coeff(m) for m in mons_quartic],
    [s.expand(TLmix).coeff(m) for m in mons_quartic],
    [s.expand(LA * LB).coeff(m) for m in mons_quartic],
])
show("Tensor basis rank", coeff_matrix.rank())

# scalar-sector response matrix and positivity
Kscalar = s.Matrix([
    [quartic_sol[kTT], quartic_sol[kTL]],
    [quartic_sol[kTL], quartic_sol[kLL]],
])
show("Scalar-sector response matrix Kscalar", Kscalar)
show("det(Kscalar)", s.factor(Kscalar.det()))
show("Eigenvalues(Kscalar)", Kscalar.eigenvals())

section("QUADRATIC G^2/r^2 / local-potential decomposition")
UA = G * mB / rAB
UB = G * mA / rAB

tauPar, betaT, betaL = s.symbols('tauPar betaT betaL')
L2_quad_ansatz = s.expand(
    G * mA * mB / (c**4 * rAB) * (
        tauPar * (UA + UB) * Cpar * vAB
        + betaT * (UA * TB + UB * TA)
        + betaL * (UA * LB + UB * LA)
    )
)

quad_target = s.expand(
    G**2 * mA * mB / (c**4 * rAB**2) * (
        t2 * (mA * vA2 + mB * vB2)
        + t3 * (mA + mB) * vAB
        + t6 * (mA * vAn**2 + mB * vBn**2)
    )
)
quad_rescaled = s.expand((L2_quad_ansatz - quad_target) * c**4 * rAB**2 / (G**2 * mA * mB))
quad_sol = s.solve(
    [
        s.expand(quad_rescaled).coeff(mA * vA2),
        s.expand(quad_rescaled).coeff(mB * vB2),
        s.expand(quad_rescaled).coeff(mA * vAB),
        s.expand(quad_rescaled).coeff(mB * vAB),
        s.expand(quad_rescaled).coeff(mA * vAn**2),
        s.expand(quad_rescaled).coeff(mB * vBn**2),
    ],
    [tauPar, betaT, betaL],
    dict=True,
)[0]
show("Quadratic local-potential solution", quad_sol)

L2_quad_constructive = s.expand(L2_quad_ansatz.subs(quad_sol))
show("L2_quad_constructive", L2_quad_constructive)
show("Quadratic residual check", s.expand(L2_quad_constructive - quad_target))

section("FULL CONSTRUCTIVE CROSS MODULE / exact reconstruction")
L2_static_cross = s.expand(s1 * G**3 * mA**2 * mB**2 / (c**4 * rAB**3))
L2_constructive = s.expand(L2_vec_kin + L2_tensor_quartic + L2_quad_constructive + L2_static_cross)

show("L2_constructive", L2_constructive)
show("Full constructive residual check", s.expand(L2_constructive - L2_cross_target))

section("LOCAL-POTENTIAL EXTENSION / first 3-body predictions for pair AB in a 3-body environment")
UA3 = G * (mB / rAB + mC / rAC)
UB3 = G * (mA / rAB + mC / rBC)

L2_AB_in_3body_env = s.expand(
    G * mA * mB / (c**4 * rAB) * (
        quad_sol[tauPar] * (UA3 + UB3) * Cpar * vAB
        + quad_sol[betaT] * (UA3 * TB + UB3 * TA)
        + quad_sol[betaL] * (UA3 * LB + UB3 * LA)
    )
)

uAB, uAC, uBC = s.symbols('uAB uAC uBC')
L2_AB_in_3body_env_sub = s.expand(
    L2_AB_in_3body_env.subs({1 / rAB: uAB, 1 / rAC: uAC, 1 / rBC: uBC})
)
show("Pair AB constructive module inside a 3-body environment", L2_AB_in_3body_env)

show(
    "Coefficient on G^2 m_A m_B m_C vAB /(c^4 rAB rAC)",
    s.simplify(L2_AB_in_3body_env_sub.coeff(G**2 * mA * mB * mC * vAB * uAB * uAC / c**4)),
)
show(
    "Coefficient on G^2 m_A m_B m_C vAB /(c^4 rAB rBC)",
    s.simplify(L2_AB_in_3body_env_sub.coeff(G**2 * mA * mB * mC * vAB * uAB * uBC / c**4)),
)
show(
    "Coefficient on G^2 m_A m_B m_C vB2 /(c^4 rAB rAC)",
    s.simplify(L2_AB_in_3body_env_sub.coeff(G**2 * mA * mB * mC * vB2 * uAB * uAC / c**4)),
)
show(
    "Coefficient on G^2 m_A m_B m_C vBn^2 /(c^4 rAB rAC)",
    s.simplify(L2_AB_in_3body_env_sub.coeff(G**2 * mA * mB * mC * vBn**2 * uAB * uAC / c**4)),
)
show(
    "Coefficient on G^2 m_A m_B m_C vA2 /(c^4 rAB rBC)",
    s.simplify(L2_AB_in_3body_env_sub.coeff(G**2 * mA * mB * mC * vA2 * uAB * uBC / c**4)),
)
show(
    "Coefficient on G^2 m_A m_B m_C vAn^2 /(c^4 rAB rBC)",
    s.simplify(L2_AB_in_3body_env_sub.coeff(G**2 * mA * mB * mC * vAn**2 * uAB * uBC / c**4)),
)

section("TAKEAWAYS")
print("1. The solved quartic 2PN cross block splits exactly into:")
print("      (i) a universal sigma = 1/2 dressing of the frozen 1PN vector wake, and")
print("      (ii) a rank-2 tensor-wake residual in a 5-dimensional projector basis.")
print("2. The quartic tensor residual coefficients are")
print("      kTT = 3/2, kS = 1/4, kM = 1, kTL = 3/4, kLL = 9/4.")
print("3. The scalar tensor subsector has response matrix")
print("      [[3/2, 3/4], [3/4, 9/4]],")
print("   whose determinant is positive and whose eigenvalues are both positive.")
print("4. The G^2/r^2 velocity block is reproduced exactly by:")
print("      - a purely parallel local-potential dressing of the 1PN vector wake with tau_parallel = 15/14,")
print("      - plus diagonal tensor-potential channels beta_T = 11/8 and beta_L = 13/4.")
print("5. The full constructive cross module matches the solved ADM-lift cross block exactly.")
print("6. The first 3-body coefficients induced by the local-potential pieces are already fixed for pair AB in a 3-body environment:")
print("      -15/4 on vAB/(rAB rAC) and vAB/(rAB rBC),")
print("      11/8 on vB^2/(rAB rAC) and vA^2/(rAB rBC),")
print("      15/8 on vBn^2/(rAB rAC) and vAn^2/(rAB rBC).")
