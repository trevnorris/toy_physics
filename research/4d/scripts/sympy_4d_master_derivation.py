#!/usr/bin/env python3
"""
sympy_4d_master_derivation.py

A lightweight SymPy "symbolic closure" harness for the Paper-7 hard-mode 4D model.

What it does (symbolic, not numerical):
  1) Derives the gauged 4D GNLS equation from the frozen matter Lagrangian
     using C1 minimal coupling (covariant derivatives).
  2) Derives the 4+1D Maxwell + localization + gauge-fixing PDE bundle
     (in explicit component form) and verifies it matches the compact divergence form.
  3) Verifies two key diagnostic identities that the numerical code must respect:
       - Continuity (charge/mass) conservation for the matter sector
       - Moving-wall chain rule: d/dt V_conf(X; a(t), L(t)) = V_a a_dot + V_L L_dot
  4) Emits compact LaTeX for the equation bundle (optional).

This file is intended to be:
  - a regression-testable "authority" on signs and coefficients,
  - the Python-side counterpart to the WL/Mathematica harness mentioned in the roadmap.

No external inputs; everything is kept symbolic.

Usage:
  python sympy_4d_master_derivation.py --check
  python sympy_4d_master_derivation.py --latex

Notes:
  - We treat psi and psib as independent fields (standard in variational derivations).
    Interpret psib as psi* at the end.
  - Geometry DOFs a(t), L(t) enter only through V_conf (and could also enter Z, m_gamma, etc).
"""

from __future__ import annotations

import argparse
import sympy as sp
from sympy.calculus.euler import euler_equations


def _coords_5():
    # 4+1D spacetime coordinates: (t, x, y, z, w)
    t, x, y, z, w = sp.symbols("t x y z w", real=True)
    return t, x, y, z, w


def derive_gauged_gnls():
    """
    Derive the gauged 4D GNLS equation from the frozen matter Lagrangian density:

      L_psi = (iħ/2)(ψ* D_t ψ - ψ (D_t ψ)*)
              - (ħ²/2m) |D_i ψ|²
              - V_conf |ψ|²
              - U(|ψ|²)

    with frozen EOS: U(ρ)= (K/4) ρ^5  =>  U'(ρ)= (5K/4) ρ^4.
    """
    t, x, y, z, w = _coords_5()
    coords = [t, x, y, z, w]
    spatial = [x, y, z, w]

    # Parameters
    hbar, m, K, q = sp.symbols("hbar m K q", positive=True, real=True)

    # Geometry DOFs
    a = sp.Function("a")(t)
    Lgeo = sp.Function("L")(t)

    # Fields
    psi = sp.Function("psi")(t, x, y, z, w)
    psib = sp.Function("psib")(t, x, y, z, w)  # interpret as psi* at the end

    # Gauge potentials (C1 minimal coupling)
    A0 = sp.Function("A0")(t, x, y, z, w)
    Ax = sp.Function("Ax")(t, x, y, z, w)
    Ay = sp.Function("Ay")(t, x, y, z, w)
    Az = sp.Function("Az")(t, x, y, z, w)
    Aw = sp.Function("Aw")(t, x, y, z, w)
    Avec = [Ax, Ay, Az, Aw]

    # Confinement potential
    Vconf = sp.Function("Vconf")(x, y, z, w, a, Lgeo)

    rho = psib * psi

    # Covariant derivatives (matches the roadmap convention)
    Dt_psi = sp.diff(psi, t) + sp.I * q / hbar * A0 * psi
    Dt_psib = sp.diff(psib, t) - sp.I * q / hbar * A0 * psib

    Di_psi = [sp.diff(psi, xi) - sp.I * q / hbar * Avec[i] * psi for i, xi in enumerate(spatial)]
    Di_psib = [sp.diff(psib, xi) + sp.I * q / hbar * Avec[i] * psib for i, xi in enumerate(spatial)]

    # |D psi|^2 in the psi/psib-independent-field convention
    Dpsi_sq = sum(Di_psi[i] * Di_psib[i] for i in range(4))

    # Frozen EOS potential
    U = K / 4 * rho**5
    Uprime = sp.Rational(5, 4) * K * rho**4  # U'(rho) for U=(K/4) rho^5

    # Lagrangian density
    Lpsi = (
        sp.I * hbar / 2 * (psib * Dt_psi - psi * Dt_psib)
        - (hbar**2 / (2 * m)) * Dpsi_sq
        - Vconf * rho
        - U
    )

    # Euler–Lagrange equations (psi and psib treated as independent)
    eqs = euler_equations(Lpsi, (psi, psib), coords)

    # Identify the equation containing d/dt psi (this is the "psi EOM")
    eq_psi = next(eq for eq in eqs if eq.lhs.has(sp.Derivative(psi, t)))
    eq_psib = next(eq for eq in eqs if eq.lhs.has(sp.Derivative(psib, t)))

    # Compact expected form:
    #   iħ D_t ψ + (ħ²/2m) D_i D_i ψ - (Vconf + U'(ρ)) ψ = 0
    D2_psi = 0
    for i, xi in enumerate(spatial):
        Ai = Avec[i]
        D_i_psi = Di_psi[i]
        D2_psi += sp.diff(D_i_psi, xi) - sp.I * q / hbar * Ai * D_i_psi

    expected_psi = sp.I * hbar * Dt_psi + (hbar**2 / (2 * m)) * D2_psi - (Vconf + Uprime) * psi

    # A compact current consistent with the roadmap:
    #   J^0 = q rho   (or neutrality-modified later)
    #   J^i = (qħ/m) Im(psib * D_i psi) = (qħ/2im)(psib D_i psi - psi D_i psib)
    J0 = q * rho
    Ji = [(q * hbar / (2 * sp.I * m)) * (psib * Di_psi[i] - psi * Di_psib[i]) for i in range(4)]

    return {
        "coords": (t, x, y, z, w),
        "params": (hbar, m, K, q),
        "fields": {"psi": psi, "psib": psib, "A0": A0, "Avec": Avec},
        "geometry": {"a": a, "L": Lgeo},
        "Vconf": Vconf,
        "rho": rho,
        "U": U,
        "Uprime": Uprime,
        "Lpsi": Lpsi,
        "EL_eq_psi": eq_psi,
        "EL_eq_psib": eq_psib,
        "expected_compact_eq_psi": expected_psi,
        "current": {"J0": J0, "Ji": Ji},
    }


def check_gnls_compact_form(bundle):
    """Assert the derived Euler–Lagrange psi equation matches the compact expected form."""
    eq_psi = bundle["EL_eq_psi"]
    expected = bundle["expected_compact_eq_psi"]
    diff = sp.simplify(eq_psi.lhs - expected)
    assert diff == 0, f"GNLS compact-form check failed: residual={diff}"


def check_continuity(bundle):
    """
    Verify the continuity identity:
       d_t rho + div_4(J) = 0
    using the *ungauged* formula for simplicity, with the EOM substituted.

    (This is a diagnostic identity in the freeze sheet, and a regression check in numerics.)
    """
    t, x, y, z, w = bundle["coords"]
    hbar, m, K, q = bundle["params"]
    psi = bundle["fields"]["psi"]
    psib = bundle["fields"]["psib"]
    rho = bundle["rho"]
    Vconf = bundle["Vconf"]
    Uprime = bundle["Uprime"]

    # For the continuity check we set gauge fields to zero to keep expressions smaller.
    # (Gauge-invariant continuity can be added later by using covariant currents.)
    subs_A0A = {bundle["fields"]["A0"]: 0}
    subs_A0A.update({Ai: 0 for Ai in bundle["fields"]["Avec"]})

    # Ungauged Laplacians
    lap_psi = sp.diff(psi, x, 2) + sp.diff(psi, y, 2) + sp.diff(psi, z, 2) + sp.diff(psi, w, 2)
    lap_psib = sp.diff(psib, x, 2) + sp.diff(psib, y, 2) + sp.diff(psib, z, 2) + sp.diff(psib, w, 2)

    # GNLS EOMs in ungauged form:
    #  iħ ψ_t = -(ħ²/2m) ∇²ψ + (V + U') ψ
    psi_t_expr = (-sp.I / hbar) * (-(hbar**2 / (2 * m)) * lap_psi + (Vconf + Uprime) * psi)
    psib_t_expr = (sp.I / hbar) * (-(hbar**2 / (2 * m)) * lap_psib + (Vconf + Uprime) * psib)

    rho_t = sp.diff(rho, t).subs({sp.diff(psi, t): psi_t_expr, sp.diff(psib, t): psib_t_expr}).subs(subs_A0A)

    # Standard Schrödinger mass/number current (ungauged):
    Jx = (hbar / (2 * sp.I * m)) * (psib * sp.diff(psi, x) - psi * sp.diff(psib, x))
    Jy = (hbar / (2 * sp.I * m)) * (psib * sp.diff(psi, y) - psi * sp.diff(psib, y))
    Jz = (hbar / (2 * sp.I * m)) * (psib * sp.diff(psi, z) - psi * sp.diff(psib, z))
    Jw = (hbar / (2 * sp.I * m)) * (psib * sp.diff(psi, w) - psi * sp.diff(psib, w))

    divJ = sp.diff(Jx, x) + sp.diff(Jy, y) + sp.diff(Jz, z) + sp.diff(Jw, w)

    residual = sp.simplify(rho_t + divJ)
    assert residual == 0, f"Continuity check failed: residual={residual}"


def check_moving_wall_chain_rule(bundle):
    """
    Verify the identity (explicitly demanded in the hard-mode plan):

      ∂_t V_conf(X; a(t), L(t)) = (∂_a V_conf) a_dot + (∂_L V_conf) L_dot
    """
    t, x, y, z, w = bundle["coords"]
    a = bundle["geometry"]["a"]
    Lgeo = bundle["geometry"]["L"]
    Vconf = bundle["Vconf"]

    V_t = sp.diff(Vconf, t)
    chain = sp.diff(Vconf, a) * sp.diff(a, t) + sp.diff(Vconf, Lgeo) * sp.diff(Lgeo, t)

    residual = sp.simplify(V_t - chain)
    assert residual == 0, f"Chain rule check failed: residual={residual}"


def derive_maxwell_localized():
    """
    Derive Maxwell + localization + gauge fixing equations in 4+1D.

    We build the bulk gauge Lagrangian density:
      L = -(Z/(4μ0)) F_{MN}F^{MN}
          - (1/(2ξμ0)) (∂_M A^M)^2
          + (mγ^2/(2μ0)) A_M A^M
          - J^M A_M

    Then Euler–Lagrange gives (for each N):
      ∂_M ( Z F^{MN} ) + mγ^2 A^N + (1/ξ) ∂^N(∂·A) = μ0 J^N

    This matches the roadmap PDE form (with Z and/or mγ chosen/frozen later).
    """
    t, x, y, z, w = _coords_5()
    coords = [t, x, y, z, w]

    mu0, xi = sp.symbols("mu0 xi", positive=True, real=True)

    # Minkowski metric signature (-,+,+,+,+) for index raising (diag only)
    eta = [-1, 1, 1, 1, 1]

    # Gauge potential components A_M
    A0 = sp.Function("A0")(t, x, y, z, w)
    Ax = sp.Function("Ax")(t, x, y, z, w)
    Ay = sp.Function("Ay")(t, x, y, z, w)
    Az = sp.Function("Az")(t, x, y, z, w)
    Aw = sp.Function("Aw")(t, x, y, z, w)
    A = [A0, Ax, Ay, Az, Aw]

    # Localization controls
    Z = sp.Function("Z")(t, x, y, z, w)
    mg2 = sp.Function("mg2")(t, x, y, z, w)

    # External + matter current (kept general here)
    J0 = sp.Function("J0")(t, x, y, z, w)
    Jx = sp.Function("Jx")(t, x, y, z, w)
    Jy = sp.Function("Jy")(t, x, y, z, w)
    Jz = sp.Function("Jz")(t, x, y, z, w)
    Jw = sp.Function("Jw")(t, x, y, z, w)
    J = [J0, Jx, Jy, Jz, Jw]

    # Field strength tensor F_{MN} = ∂_M A_N - ∂_N A_M
    F = {(M, N): sp.diff(A[N], coords[M]) - sp.diff(A[M], coords[N]) for M in range(5) for N in range(5)}
    # Raised components: F^{MN} = η^{MM} η^{NN} F_{MN} for diagonal metric
    Fup = {(M, N): eta[M] * eta[N] * F[(M, N)] for M in range(5) for N in range(5)}

    # Contractions
    FF = sum(F[(M, N)] * Fup[(M, N)] for M in range(5) for N in range(5))  # F_{MN}F^{MN}
    divA = sum(sp.diff(eta[M] * A[M], coords[M]) for M in range(5))  # ∂_M A^M
    A2 = sum(eta[M] * A[M] ** 2 for M in range(5))  # A_M A^M

    L = -(Z / (4 * mu0)) * FF - (1 / (2 * xi * mu0)) * divA**2 + (mg2 / (2 * mu0)) * A2 - sum(J[M] * A[M] for M in range(5))

    eqs = euler_equations(L, tuple(A), tuple(coords))

    # Compact expected residual (multiplied by μ0):
    #   R^N = ∂_M(Z F^{MN}) + mg2 A^N + (1/ξ) ∂^N(divA) - μ0 J^N
    residuals = {}
    for N in range(5):
        exp_res = sum(sp.diff(Z * Fup[(M, N)], coords[M]) for M in range(5))
        exp_res += mg2 * (eta[N] * A[N])  # A^N
        exp_res += (1 / xi) * (eta[N] * sp.diff(divA, coords[N]))  # ∂^N(divA)
        exp_res -= mu0 * J[N]
        residuals[N] = sp.simplify(mu0 * eqs[N].lhs - exp_res)

    return {
        "coords": (t, x, y, z, w),
        "params": (mu0, xi),
        "fields": {"A": A, "Z": Z, "mg2": mg2, "J": J},
        "EL_eqs": eqs,
        "expected_residual_diffs": residuals,
        "divA": divA,
        "F": F,
        "Fup": Fup,
    }


def check_maxwell_bundle(bundle):
    diffs = bundle["expected_residual_diffs"]
    bad = {k: v for k, v in diffs.items() if v != 0}
    assert not bad, f"Maxwell bundle check failed for components: {bad}"


def _latex_block(title: str, expr) -> str:
    return "\n".join([f"% --- {title}", sp.latex(expr), ""])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="run assert-style symbolic checks")
    ap.add_argument("--latex", action="store_true", help="print compact LaTeX for the equation bundle")
    args = ap.parse_args()

    gnls = derive_gauged_gnls()
    maxw = derive_maxwell_localized()

    if args.check:
        check_gnls_compact_form(gnls)
        check_continuity(gnls)
        check_moving_wall_chain_rule(gnls)
        check_maxwell_bundle(maxw)
        print("All symbolic checks passed.")

    if args.latex:
        psi = gnls["fields"]["psi"]
        psib = gnls["fields"]["psib"]
        print("% === Gauged GNLS equation (psi) ===")
        print(sp.latex(sp.Eq(gnls["expected_compact_eq_psi"], 0)))
        print()

        print("% === Continuity ingredients (ungauged form shown) ===")
        t, x, y, z, w = gnls["coords"]
        hbar, m, K, q = gnls["params"]
        Jx = (hbar / (2 * sp.I * m)) * (psib * sp.diff(psi, x) - psi * sp.diff(psib, x))
        Jy = (hbar / (2 * sp.I * m)) * (psib * sp.diff(psi, y) - psi * sp.diff(psib, y))
        Jz = (hbar / (2 * sp.I * m)) * (psib * sp.diff(psi, z) - psi * sp.diff(psib, z))
        Jw = (hbar / (2 * sp.I * m)) * (psib * sp.diff(psi, w) - psi * sp.diff(psib, w))
        print(sp.latex(sp.Eq(sp.diff(gnls["rho"], t) + sp.diff(Jx, x) + sp.diff(Jy, y) + sp.diff(Jz, z) + sp.diff(Jw, w), 0)))
        print()

        print("% === Maxwell + localization + gauge fixing (compact residual form) ===")
        t2, x2, y2, z2, w2 = maxw["coords"]
        mu0, xi = maxw["params"]
        A = maxw["fields"]["A"]
        Z = maxw["fields"]["Z"]
        mg2 = maxw["fields"]["mg2"]
        J = maxw["fields"]["J"]
        eta = [-1, 1, 1, 1, 1]
        coords = [t2, x2, y2, z2, w2]
        divA = maxw["divA"]
        # Print componentwise equation skeleton
        for N, name in enumerate(["0", "x", "y", "z", "w"]):
            # residual: ∂_M(Z F^{MN}) + mg2 A^N + (1/ξ) ∂^N(divA) = μ0 J^N
            # We'll print as a symbolic equality for each component.
            Fup = maxw["Fup"]
            lhs = sum(sp.diff(Z * Fup[(M, N)], coords[M]) for M in range(5))
            lhs += mg2 * (eta[N] * A[N])
            lhs += (1 / xi) * (eta[N] * sp.diff(divA, coords[N]))
            eq = sp.Eq(lhs, mu0 * J[N])
            print(f"% N={name}")
            print(sp.latex(eq))
            print()

    if not (args.check or args.latex):
        # Default: print a small human-readable summary
        print("Derived bundles are constructed. Use --check to run symbolic regression checks.")
        print("Use --latex to print LaTeX for the compact equation bundle.")


if __name__ == "__main__":
    main()
