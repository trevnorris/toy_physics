
#!/usr/bin/env python3
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

def expect_zero(name: str, expr):
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


"""
5pn_stage169_similarity_orbit_closure.py

Executable SymPy audit for the moving-throat Stage 169 theorem.

What this script does
---------------------
1. Builds the exact monomial-drift matrix M_* for
      (C_tr,* , C_nt,* , epsilon_eta).
2. Proves rank(M_*) = 3 on the constructive coherent branch and therefore
      dim ker M_* = 5.
3. Constructs the exact five-parameter multiplicative similarity orbit G_*.
4. Proves that G_* preserves the three direct monomials exactly and that
      ker M_* = T_id G_*.
"""

banner("STAGE 169 — EXACT MICROSCOPIC SIMILARITY ORBIT")

# ---------------------------------------------------------------------------
# I. Exact monomial-drift matrix
# ---------------------------------------------------------------------------

subbanner("I. Exact monomial-drift matrix M_*")

chi0s, deltaUs = sp.symbols("chi0_star deltaU_star", positive=True, real=True)
Estar, Fstar = sp.symbols("E_star F_star", real=True)

lambda1, c1, gamma1 = sp.symbols("lambda_1 c_1 gamma_1", real=True)
kappaU, kappa_eta, kappaW = sp.symbols("kappa_U kappa_eta kappa_W", real=True)
mu1, tau1 = sp.symbols("mu_1 tau_1", real=True)
vars_vec = (lambda1, c1, gamma1, kappaU, kappa_eta, kappaW, mu1, tau1)

dlog_Ctr = sp.expand((1 + deltaUs) * (gamma1 + c1 - kappaU) + (1 + chi0s) * (tau1 - kappaU))
dlog_Cnt = sp.expand((2 * lambda1 + mu1 - kappa_eta - 2 * kappaW)
                     + Estar * (2 * gamma1 + 2 * lambda1 - kappaU - kappaW)
                     - Fstar * (tau1 - kappaU))
dlog_eps = sp.expand(2 * c1 - kappaU - kappa_eta)

M = sp.Matrix([
    [sp.expand(dlog_Ctr).coeff(v) for v in vars_vec],
    [sp.expand(dlog_Cnt).coeff(v) for v in vars_vec],
    [sp.expand(dlog_eps).coeff(v) for v in vars_vec],
])

M_expected = sp.Matrix([
    [0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s],
    [2 * (1 + Estar), 0, 2 * Estar, Fstar - Estar, -1, -(2 + Estar), 1, -Fstar],
    [0, 2, 0, -1, -1, 0, 0, 0],
])
print("M_* =")
sp.pprint(M)
expect_zero("M_* - expected", M - M_expected)

minor = M[:, [7, 4, 6]]
minor_det = sp.simplify(minor.det())
print("det M_*^(tau_1, kappa_eta, mu_1) =", minor_det)
expect_zero("determinant - (1 + chi0_star)", minor_det - (1 + chi0s))
print("So rank M_* = 3 and dim ker M_* = 5 on the constructive branch.")

# ---------------------------------------------------------------------------
# II. Exact five-parameter similarity orbit
# ---------------------------------------------------------------------------

subbanner("II. Exact five-parameter multiplicative similarity orbit")

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

Delta_x_orbit = sp.Matrix([
    Delta_lambda, Delta_c, Delta_gamma, Delta_U,
    Delta_Keta, Delta_KW, Delta_mu, Delta_T
])

print("Delta_x_orbit =")
sp.pprint(Delta_x_orbit)

expect_zero("orbit preserves all three monomials", M * Delta_x_orbit)

# ---------------------------------------------------------------------------
# III. Tangent-space equivalence
# ---------------------------------------------------------------------------

subbanner("III. Tangent-space equivalence ker M_* = T_id G_*")

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

# Image(O_*) ⊂ ker M_*
expect_zero("M_* O_*", M * orbit_matrix)

# rank(O_*) = 5 from the free coordinates (lambda,c,gamma,K_U,K_W)
expect_zero("free-coordinate block is identity", orbit_matrix[[0,1,2,3,5], :] - sp.eye(5))

print("Because image(O_*) has dimension 5 and dim ker M_* = 5,")
print("we conclude image(O_*) = ker M_* = T_id G_*.")

banner("FINAL STAGE-169 LEDGER")
print("1. The three direct microscopic monomials define the exact rank-3 drift matrix M_*.")
print("2. On the constructive coherent branch, det M_*^(tau_1,kappa_eta,mu_1) = 1 + chi0_* > 0.")
print("3. Therefore dim ker M_* = 5.")
print("4. The exact five-parameter multiplicative similarity orbit G_* preserves the monomials exactly.")
print("5. Its tangent space is exactly ker M_*, so the zero-defect branch is a codimension-3 similarity manifold, not fine-tuning.")
