#!/usr/bin/env python3

"""
moving_throat_pde_stage35_nontwin_asymmetry_threshold_sympy_audit.py

SymPy audit for Stage 35 of the moving-throat PDE program.

Checks
------
1. Exact branch-product regime classification in terms of C_mix and Pi_tr.
2. Exact required support ratio zeta_req(Pi_tr,C_mix,eps).
3. Exact monotonicity of zeta_req with respect to Pi_tr on the support-needed branch.
4. Exact excess beyond the symmetric twin.
5. Exact lowest-lane asymmetry thresholds:
      - overlap boost at fixed stiffness,
      - support softening at fixed overlap.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


# Symbols
Pi_tr, Cmix, eps = sp.symbols("Pi_tr C_mix eps", positive=True, real=True)
KW, Kphi0, Omega0 = sp.symbols("K_W K_phi0 Omega0", positive=True, real=True)

banner("STAGE 35 — EXACT BRANCH-PRODUCT REGIME CLASSIFICATION")

Sreq = sp.simplify(Pi_tr / Cmix)
zeta_req = sp.simplify((Sreq - 1) / (1 + eps * (Sreq - 2)))

print("S_req = Pi_tr / C_mix =")
sp.pprint(Sreq)
print("\nzeta_req(Pi_tr,C_mix,eps) =")
sp.pprint(zeta_req)

expect_zero("zeta_req at Pi_tr = C_mix", zeta_req.subs(Pi_tr, Cmix))
expect_zero("zeta_req at Pi_tr = 2 C_mix minus 1", zeta_req.subs(Pi_tr, 2 * Cmix) - 1)

dz_dPi = sp.simplify(sp.diff(zeta_req, Pi_tr))
print("\nd zeta_req / d Pi_tr =")
sp.pprint(sp.factor(dz_dPi))

# Verify the exact closed form used in the note.
dz_expected = sp.simplify(Cmix * (1 - eps) / (Cmix - eps * (2 * Cmix - Pi_tr)) ** 2)
expect_zero("dzeta_req/dPi - expected", dz_dPi - dz_expected)

banner("EXACT EXCESS BEYOND THE SYMMETRIC TWIN")
Delta_zeta = sp.simplify(zeta_req - 1)
print("Delta_zeta = zeta_req - 1 =")
sp.pprint(Delta_zeta)

Delta_expected = sp.simplify((1 - eps) * (Pi_tr - 2 * Cmix) / (Cmix - eps * (2 * Cmix - Pi_tr)))
expect_zero("Delta_zeta - expected", Delta_zeta - Delta_expected)

banner("GENERAL LOWEST-LANE ASYMMETRY THRESHOLDS")
zeta0_phys = sp.simplify(KW * Omega0**2 / Kphi0)
print("zeta_0^(phys) =")
sp.pprint(zeta0_phys)

# Exact threshold equations
Omega0_req_sq = sp.simplify(zeta_req * Kphi0 / KW)
Kphi0_req = sp.simplify(KW * Omega0**2 / zeta_req)

print("\nOmega_(0,req)^2 =")
sp.pprint(Omega0_req_sq)
print("\nK_(phi,0)^(req) =")
sp.pprint(Kphi0_req)

Omega0_sq_sym = sp.symbols("Omega0_sq", positive=True, real=True)
sol_Omega0_sq = sp.solve(
    (KW * Omega0_sq_sym / Kphi0) - zeta_req, Omega0_sq_sym
)
assert len(sol_Omega0_sq) == 1, "expected unique Omega0^2 solution"
expect_zero(
    "solve(zeta_phys = zeta_req) for Omega0^2 - expected",
    sp.simplify(sol_Omega0_sq[0] - Omega0_req_sq),
)

sol_Kphi0 = sp.solve(
    (KW * Omega0**2 / Kphi0) - zeta_req, Kphi0
)
assert len(sol_Kphi0) == 1, "expected unique Kphi0 solution"
expect_zero(
    "solve(zeta_phys = zeta_req) for Kphi0 - expected",
    sp.simplify(sol_Kphi0[0] - Kphi0_req),
)

banner("SYMMETRIC TWIN DIAGNOSTICS")
zeta_twin = sp.simplify(zeta0_phys.subs({Omega0: 1, Kphi0: KW}))
expect_zero("zeta_0^(twin) - 1", zeta_twin - 1)

Omega_req_equal_stiff = sp.simplify(sp.sqrt(zeta_req))
Kphi_req_equal_ov = sp.simplify(KW / zeta_req)
soft_frac = sp.simplify(1 - Kphi_req_equal_ov / KW)

print("At equal stiffness K_phi0 = K_W, required overlap boost is")
sp.pprint(Omega_req_equal_stiff)

print("\nAt equal overlap Omega0 = 1, required softened stiffness is")
sp.pprint(Kphi_req_equal_ov)

print("\nExact softening fraction at equal overlap =")
sp.pprint(soft_frac)

soft_expected = sp.simplify((1 - eps) * (Pi_tr - 2 * Cmix) / (Pi_tr - Cmix))
expect_zero("softening fraction - expected", soft_frac - soft_expected)

banner("STAGE 35 THEOREM LEDGER")
print("Exact regime split:")
print("  Pi_tr <= C_mix            : mixed-only branch already enough")
print("  C_mix < Pi_tr <= 2 C_mix  : symmetric lowest twin branch is enough")
print("  Pi_tr > 2 C_mix           : non-twin asymmetry is required")
print()
print("Exact support-ratio formula:")
print("  zeta_req = (Pi_tr - C_mix) / [ C_mix - eps (2 C_mix - Pi_tr) ]")
print()
print("Exact lowest-lane rescue thresholds:")
print("  Omega_0^2 >= zeta_req K_phi0 / K_W")
print("  K_phi0 <= K_W Omega_0^2 / zeta_req")
