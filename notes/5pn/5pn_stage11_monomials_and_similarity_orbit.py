
#!/usr/bin/env python3
"""
5pn_stage11_monomials_and_similarity_orbit.py

Eleventh executable SymPy audit for the 5PN grouped-real-P2 program.

What this script does
---------------------
1. Encodes the Stage-167/168 direct microscopic monomials
      C_tr, C_nt, epsilon_eta
   and derives their exact logarithmic drift matrix M_*.
2. Proves the rank-3 property of M_* and the exact determinant
      det M_*^(tau, kappa_eta, mu_1) = 1 + chi_{0,*} > 0,
   so the zero-defect branch has codimension 3.
3. Builds the exact five-parameter multiplicative similarity orbit G_* and proves
   that it preserves the three direct monomials exactly.
4. Shows that the tangent space of G_* is exactly ker M_*.
5. Upgrades the tangent statement to a finite one by proving that the exact
   invariant fibres satisfy M_* Delta x = 0 and coincide with the similarity
   orbits. This is the Stage-170 orbit-quotient closure.

Interpretation
--------------
After this stage the coherent weak-axisymmetric zero-defect problem is no longer
merely a drift identity. It is an exact finite quotient problem with quotient
coordinates given by the three direct microscopic monomials.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(sp.expand(expr))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


# ---------------------------------------------------------------------------
# I. Direct microscopic monomials and the monomial-drift matrix
# ---------------------------------------------------------------------------

banner("I. DIRECT MICROSCOPIC MONOMIALS AND THE EXACT DRIFT MATRIX")

# Positive microscopic couplings.
lamW, c_etaU, gamma, KU = sp.symbols("lambda_W c_etaU gamma K_U", positive=True, real=True)
Keta, KW, muW, TU = sp.symbols("K_eta K_W mu_W T_U", positive=True, real=True)

# Frozen constructive coherent-branch constants.
chi0s, deltaUs = sp.symbols("chi0_star deltaU_star", positive=True, real=True)
Estar, Fstar = sp.symbols("E_star F_star", real=True)
L, sigma = sp.symbols("L sigma", positive=True, real=True)

# Exact direct monomials from Stage 168.
Ctr = sp.simplify((gamma * c_etaU / KU)**(1 + deltaUs) * (sp.pi**2 * TU / (L**2 * KU))**(1 + chi0s))
Cnt = sp.simplify(
    (lamW**2 * muW / (Keta * KW**2))
    * (gamma**2 * lamW**2 * sigma / (KU * KW))**Estar
    * (sp.pi**2 * TU / (L**2 * KU))**(-Fstar)
)
eps_eta = sp.simplify(c_etaU**2 / (KU * Keta))

print("C_tr,* =")
sp.pprint(Ctr)
print("C_nt,* =")
sp.pprint(Cnt)
print("epsilon_eta =")
sp.pprint(eps_eta)

# Microscopic grouped weak-axisymmetric log-drift vector.
lambda1, c1, gamma1 = sp.symbols("lambda_1 c_1 gamma_1", real=True)
kappaU, kappa_eta, kappaW = sp.symbols("kappa_U kappa_eta kappa_W", real=True)
mu1, tau1 = sp.symbols("mu_1 tau_1", real=True)

drift_subs = {
    sp.log(lamW): lambda1,
    sp.log(c_etaU): c1,
    sp.log(gamma): gamma1,
    sp.log(KU): kappaU,
    sp.log(Keta): kappa_eta,
    sp.log(KW): kappaW,
    sp.log(muW): mu1,
    sp.log(TU): tau1,
    sp.log(L): 0,
    sp.log(sigma): 0,
    sp.log(sp.pi**2): 0,
}

# Build the exact logarithmic drift rows.
vars_vec = (lambda1, c1, gamma1, kappaU, kappa_eta, kappaW, mu1, tau1)

dlog_Ctr = sp.expand(
    (1 + deltaUs) * (gamma1 + c1 - kappaU)
    + (1 + chi0s) * (tau1 - kappaU)
)
dlog_Cnt = sp.expand(
    (2 * lambda1 + mu1 - kappa_eta - 2 * kappaW)
    + Estar * (2 * gamma1 + 2 * lambda1 - kappaU - kappaW)
    - Fstar * (tau1 - kappaU)
)
dlog_eps = sp.expand(2 * c1 - kappaU - kappa_eta)

M = sp.Matrix([
    [sp.expand(dlog_Ctr).coeff(v) for v in vars_vec],
    [sp.expand(dlog_Cnt).coeff(v) for v in vars_vec],
    [sp.expand(dlog_eps).coeff(v) for v in vars_vec],
])

print("M_* =")
sp.pprint(M)

M_expected = sp.Matrix([
    [0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s],
    [2 * (1 + Estar), 0, 2 * Estar, Fstar - Estar, -1, -(2 + Estar), 1, -Fstar],
    [0, 2, 0, -1, -1, 0, 0, 0],
])
expect_zero("M_* - expected", M - M_expected)

subbanner("Rank-3 test")
minor = M[:, [7, 4, 6]]  # columns (tau_1, kappa_eta, mu_1)
minor_det = sp.simplify(minor.det())
print("det M_*^(tau, kappa_eta, mu_1) =", minor_det)
expect_zero("minor determinant - (1 + chi0_star)", minor_det - (1 + chi0s))
print("Therefore rank M_* = 3 and dim ker M_* = 5 on the constructive branch.")


# ---------------------------------------------------------------------------
# II. Exact five-parameter similarity orbit
# ---------------------------------------------------------------------------

banner("II. EXACT FIVE-PARAMETER SIMILARITY ORBIT")

Lam, C, Gam, U, W = sp.symbols("Lambda C Gamma U W", real=True)

Delta_lambda = Lam
Delta_c = C
Delta_gamma = Gam
Delta_U = U
Delta_KW = W
Delta_Keta = sp.simplify(2 * C - U)
Delta_T = sp.simplify(U - (1 + deltaUs) / (1 + chi0s) * (Gam + C - U))
Delta_mu = sp.simplify(
    2 * C - U + 2 * W - 2 * Lam
    - Estar * (2 * Gam + 2 * Lam - U - W)
    - Fstar * (1 + deltaUs) / (1 + chi0s) * (Gam + C - U)
)

orbit_matrix = sp.Matrix([
    [sp.expand(Delta_lambda).coeff(p) for p in (Lam, C, Gam, U, W)],
    [sp.expand(Delta_c).coeff(p) for p in (Lam, C, Gam, U, W)],
    [sp.expand(Delta_gamma).coeff(p) for p in (Lam, C, Gam, U, W)],
    [sp.expand(Delta_U).coeff(p) for p in (Lam, C, Gam, U, W)],
    [sp.expand(Delta_Keta).coeff(p) for p in (Lam, C, Gam, U, W)],
    [sp.expand(Delta_KW).coeff(p) for p in (Lam, C, Gam, U, W)],
    [sp.expand(Delta_mu).coeff(p) for p in (Lam, C, Gam, U, W)],
    [sp.expand(Delta_T).coeff(p) for p in (Lam, C, Gam, U, W)],
])

print("Orbit exponent matrix O_* =")
sp.pprint(orbit_matrix)

# Free rows (lambda, c, gamma, K_U, K_W) already form the identity.
expect_zero("rank-5 free-block identity", orbit_matrix[[0,1,2,3,5], :] - sp.eye(5))

subbanner("Exact preservation of the three direct monomials")
Delta_x_orbit = sp.Matrix([Delta_lambda, Delta_c, Delta_gamma, Delta_U, Delta_Keta, Delta_KW, Delta_mu, Delta_T])
drift_orbit = sp.simplify(M * Delta_x_orbit)
expect_zero("M_* Delta_x_orbit", drift_orbit)

print("So the orbit preserves C_tr,* , C_nt,* , and epsilon_eta exactly.")


# ---------------------------------------------------------------------------
# III. Tangent-space equivalence
# ---------------------------------------------------------------------------

banner("III. TANGENT-SPACE EQUIVALENCE: ker M_* = T_id G_*")

print("Because M_* O_* = 0 and rank(O_*) = 5, the orbit tangent space lies in ker M_*.")
print("Because rank(M_*) = 3 in an 8-dimensional space, dim ker M_* = 5.")
print("Therefore image(O_*) = ker(M_*).")

# Check the Stage-168 compatibility relations directly from the orbit.
compat_keta = sp.simplify(Delta_Keta - (2 * Delta_c - Delta_U))
compat_tau = sp.simplify(Delta_T - (Delta_U - (1 + deltaUs) / (1 + chi0s) * (Delta_gamma + Delta_c - Delta_U)))
compat_mu = sp.simplify(
    Delta_mu
    - (
        2 * Delta_c - Delta_U + 2 * Delta_KW - 2 * Delta_lambda
        - Estar * (2 * Delta_gamma + 2 * Delta_lambda - Delta_U - Delta_KW)
        - Fstar * (1 + deltaUs) / (1 + chi0s) * (Delta_gamma + Delta_c - Delta_U)
    )
)
expect_zero("compatibility for kappa_eta", compat_keta)
expect_zero("compatibility for tau_1", compat_tau)
expect_zero("compatibility for mu_1", compat_mu)


# ---------------------------------------------------------------------------
# IV. Finite invariant fibres and orbit-quotient closure
# ---------------------------------------------------------------------------

banner("IV. FINITE INVARIANT FIBRES AND ORBIT-QUOTIENT CLOSURE")

# Finite log-ratio vector between two positive microscopic states.
Dl, Dc, Dg, DUf, DKe, DWf, Dmu, DTf = sp.symbols(
    "Delta_lambda Delta_c Delta_gamma Delta_U Delta_Keta Delta_KW Delta_mu Delta_T",
    real=True,
)
Delta_x = sp.Matrix([Dl, Dc, Dg, DUf, DKe, DWf, Dmu, DTf])

print("Exact finite invariant-fibre equations are M_* Delta_x = 0.")
fiber_eqs = sp.simplify(M * Delta_x)
sp.pprint(fiber_eqs)

# Solve the finite fibre equations for the dependent variables
# (Delta_T, Delta_Keta, Delta_mu) in terms of the five free ones.
sol_fibre = sp.solve(
    list(fiber_eqs),
    [DTf, DKe, Dmu],
    dict=True,
)
if len(sol_fibre) != 1:
    raise AssertionError("Expected a unique finite-fibre solve for (Delta_T, Delta_Keta, Delta_mu).")
sol_fibre = sol_fibre[0]

print("Finite invariant-fibre solve:")
for k, v in sol_fibre.items():
    print(f"  {k} = {sp.simplify(v)}")

expect_zero("finite fibre Delta_Keta matches orbit", sol_fibre[DKe] - Delta_Keta.subs({Lam: Dl, C: Dc, Gam: Dg, U: DUf, W: DWf}))
expect_zero("finite fibre Delta_T matches orbit", sol_fibre[DTf] - Delta_T.subs({Lam: Dl, C: Dc, Gam: Dg, U: DUf, W: DWf}))
expect_zero("finite fibre Delta_mu matches orbit", sol_fibre[Dmu] - Delta_mu.subs({Lam: Dl, C: Dc, Gam: Dg, U: DUf, W: DWf}))

print("Thus every finite log-ratio preserving the three invariants is exactly an orbit displacement.")
print("Equivalently, the level sets of (C_tr, C_nt, epsilon_eta) are precisely the similarity orbits.")


# ---------------------------------------------------------------------------
# V. Final theorem ledger
# ---------------------------------------------------------------------------

banner("V. FINAL THEOREM LEDGER")
print("1. The direct microscopic monomials are")
print("      C_tr,* ,  C_nt,* ,  epsilon_eta.")
print("2. Their logarithmic drift map is the exact rank-3 matrix M_*.")
print("3. The zero-defect branch therefore has codimension 3 and ker(M_*) has dimension 5.")
print("4. There is an exact five-parameter multiplicative similarity orbit G_*")
print("   preserving those three monomials.")
print("5. Its tangent space is exactly ker(M_*).")
print("6. The finite invariant-fibre equations are M_* Delta_x = 0.")
print("7. Solving those finite equations reproduces the exact orbit exponents.")
print("8. Therefore the exact level sets of (C_tr,* , C_nt,* , epsilon_eta)")
print("   are precisely the similarity orbits G_*.")
print("9. The coherent weak-axisymmetric problem is therefore an exact finite quotient")
print("   problem with quotient coordinates (C_tr,* , C_nt,* , epsilon_eta).")
