#!/usr/bin/env python3
"""
Stage 44 SymPy audit — microscopic gain thresholds and operator phase diagram.
"""
from __future__ import annotations
import sympy as sp


def banner(t: str):
    print("\n" + "="*88)
    print(t)
    print("="*88)


def expect_zero(name: str, expr):
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} != 0")


banner("STAGE 44 — MICROSCOPIC GAIN THRESHOLDS")

# symbols
chi, Lam, KX, TX, L = sp.symbols('chi_sigma Lambda_phi K_X T_X L', positive=True, real=True)
kappa, eta, Pe_req = sp.symbols('kappa eta Pe_req', positive=True, real=True)
alpha = sp.sqrt(kappa)

# endpoint functions from Stage 41
Delta0 = sp.simplify(eta * (sp.cosh(alpha) - 1) / (alpha**2 * (alpha*sp.sinh(alpha) + eta*sp.cosh(alpha))))
Delta_inf = sp.simplify((sp.cosh(alpha) + (eta/alpha)*sp.sinh(alpha) - 1) / (alpha*sp.sinh(alpha) + eta*sp.cosh(alpha)))
print("Delta_0 =", Delta0)
print("Delta_inf =", Delta_inf)

# microscopic Xi and gain
Xi_micro = sp.simplify(chi * Lam**2 * L**2 / TX)
G_def = sp.simplify(chi * Lam**2 / KX)
print("Xi_micro =", Xi_micro)
print("G_micro definition =", G_def)
expect_zero("Xi_micro - kappa G_micro",
            Xi_micro.subs(TX, KX*L**2/kappa) - kappa*G_def)

# threshold surfaces
Xi_fail = sp.simplify(Pe_req / Delta_inf)
Xi_suff = sp.simplify(Pe_req / Delta0)
G_fail = sp.simplify(Xi_fail / kappa)
G_suff = sp.simplify(Xi_suff / kappa)
print("G_fail =", G_fail)
print("G_suff =", G_suff)

chi_fail = sp.simplify(TX * Pe_req / (Lam**2 * L**2 * Delta_inf))
chi_suff = sp.simplify(TX * Pe_req / (Lam**2 * L**2 * Delta0))
Lam2_fail = sp.simplify(TX * Pe_req / (chi * L**2 * Delta_inf))
Lam2_suff = sp.simplify(TX * Pe_req / (chi * L**2 * Delta0))
print("chi_fail =", chi_fail)
print("chi_suff =", chi_suff)
print("Lambda^2_fail =", Lam2_fail)
print("Lambda^2_suff =", Lam2_suff)
expect_zero("chi_fail from G_fail",
            chi_fail.subs(TX, KX*L**2/kappa) - (KX/Lam**2) * G_fail)
expect_zero("chi_suff from G_suff",
            chi_suff.subs(TX, KX*L**2/kappa) - (KX/Lam**2) * G_suff)
expect_zero("Lambda^2_fail from G_fail",
            Lam2_fail.subs(TX, KX*L**2/kappa) - (KX/chi) * G_fail)
expect_zero("Lambda^2_suff from G_suff",
            Lam2_suff.subs(TX, KX*L**2/kappa) - (KX/chi) * G_suff)

# exact soft-support limits kappa -> 0+
Delta0_k0 = sp.simplify(sp.limit(Delta0, kappa, 0, dir='+'))
Delta_inf_k0 = sp.simplify(sp.limit(Delta_inf, kappa, 0, dir='+'))
print("lim_{kappa->0+} Delta_0 =", Delta0_k0)
print("lim_{kappa->0+} Delta_inf =", Delta_inf_k0)
expect_zero("Delta0 soft-support limit - 1/2", Delta0_k0 - sp.Rational(1,2))
expect_zero("Delta_inf soft-support limit - 1", Delta_inf_k0 - 1)
# leading threshold scalings
expect_zero("kappa*G_fail soft-support limit - Pe_req",
            sp.simplify(sp.limit(kappa*G_fail, kappa, 0, dir='+') - Pe_req))
expect_zero("kappa*G_suff soft-support limit - 2 Pe_req",
            sp.simplify(sp.limit(kappa*G_suff, kappa, 0, dir='+') - 2*Pe_req))

# exact eta -> infinity limit
Delta0_eta_inf = sp.simplify(sp.limit(Delta0, eta, sp.oo))
Delta_inf_eta_inf = sp.simplify(sp.limit(Delta_inf, eta, sp.oo))
print("lim_{eta->oo} Delta_0 =", Delta0_eta_inf)
print("lim_{eta->oo} Delta_inf =", Delta_inf_eta_inf)
expect_zero("Delta0 eta->oo formula",
            Delta0_eta_inf - (1 - sp.sech(alpha))/kappa)
expect_zero("Delta_inf eta->oo formula",
            Delta_inf_eta_inf - sp.tanh(alpha)/alpha)

G_fail_inf = sp.simplify(sp.limit(G_fail, eta, sp.oo))
G_suff_inf = sp.simplify(sp.limit(G_suff, eta, sp.oo))
print("G_fail^(inf) =", G_fail_inf)
print("G_suff^(inf) =", G_suff_inf)
expect_zero("G_fail^(inf) formula",
            G_fail_inf - Pe_req/(alpha*sp.tanh(alpha)))
expect_zero("G_suff^(inf) formula",
            G_suff_inf - Pe_req/(1 - sp.sech(alpha)))

# combined stiff-support asymptotics on the compliant branch
z = sp.symbols('z', positive=True, real=True)
G_fail_inf_z = sp.simplify(G_fail_inf.subs(alpha, z))
G_suff_inf_z = sp.simplify(G_suff_inf.subs(alpha, z))
fail_stiff = sp.simplify(sp.limit(z * G_fail_inf_z, z, sp.oo) - Pe_req)
suff_stiff = sp.simplify(sp.limit(G_suff_inf_z, z, sp.oo) - Pe_req)
expect_zero("stiff-support compliant-mouth limit: sqrt(kappa)*G_fail -> Pe_req", fail_stiff)
expect_zero("stiff-support compliant-mouth limit: G_suff -> Pe_req", suff_stiff)

banner("STAGE 44 THEOREM LEDGER")
print("1. The microscopic coupling is Xi_micro = kappa * G_micro with G_micro = chi_sigma Lambda_phi^2 / K_X.")
print("2. The exact operator phase diagram is set by G_fail = Pe_req/[kappa Delta_inf] and G_suff = Pe_req/[kappa Delta_0].")
print("3. Soft support is strongly disfavored: G_fail ~ Pe_req/kappa, G_suff ~ 2 Pe_req/kappa as kappa -> 0+.")
print("4. In the highly compliant-mouth branch, G_fail^(inf) = Pe_req/[sqrt(kappa) tanh(sqrt(kappa))] and G_suff^(inf) = Pe_req/[1-sech(sqrt(kappa))].")
print("5. For stiff support with a compliant mouth, G_fail ~ Pe_req/sqrt(kappa) while G_suff -> Pe_req.")
