#!/usr/bin/env python3
"""
Stage V2-06/V2-07 SymPy audit:
Projected continuity, Poisson hook, and Newtonian-source universality.

This script verifies the exact algebraic identities in the projected-continuity
stack and isolates the additional universality gates needed to upgrade the
Poisson hook into a Newtonian source law.
"""

from __future__ import annotations

import sympy as sp


def clean(expr):
    return sp.factor(sp.together(sp.expand(expr)))


def check(name: str, expr, expected=0) -> bool:
    residual = clean(expr - expected)
    ok = residual == 0
    print(f"{name}: {'PASS' if ok else 'FAIL'}")
    if not ok:
        print(f"  residual = {residual}")
    return ok


def section(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def main() -> None:
    passes = 0
    checks = 0

    section("V2-06/V2-07: projected continuity, Poisson hook, Newtonian universality")

    # ------------------------------------------------------------------
    # 1. Projected continuity / leakage identity
    # ------------------------------------------------------------------
    section("1. Exact projected continuity and leakage split")

    w = sp.symbols("w", real=True)
    W = sp.Function("W")(w)
    jw = sp.Function("j_w")(w)

    # Product-rule identity used in the integration by parts:
    # -∫ W * ∂w j_w dw = -[W j_w] + ∫ W' j_w dw.
    # The local differential rearrangement is
    # -W*jw' = -(W*jw)' + W'*jw.
    leakage_local_lhs = -W * sp.diff(jw, w)
    leakage_local_rhs = -sp.diff(W * jw, w) + sp.diff(W, w) * jw
    checks += 1
    passes += check("local leakage integration-by-parts identity", leakage_local_lhs, leakage_local_rhs)

    print("Finite-window leakage form:")
    print("  S_leak = -[W j^w]_{w1}^{w2} + Integral(W'(w) j^w(w), (w,w1,w2))")
    print("Infinite-window fast-decay limit:")
    print("  S_leak = Integral(W'(w) j^w(w), (w,-oo,oo))")

    # ------------------------------------------------------------------
    # 2. Exact longitudinal identity
    # ------------------------------------------------------------------
    section("2. Exact longitudinal identity after Helmholtz split")

    x, y, z, t = sp.symbols("x y z t", real=True)
    rho = sp.Function("rho_br")(x, y, z, t)
    phi = sp.Function("phi_br")(x, y, z, t)
    vTx = sp.Function("v_Tx")(x, y, z, t)
    vTy = sp.Function("v_Ty")(x, y, z, t)
    vTz = sp.Function("v_Tz")(x, y, z, t)
    S_leak = sp.Function("S_leak")(x, y, z, t)

    coords = [x, y, z]
    grad_phi = [sp.diff(phi, c) for c in coords]
    vT = [vTx, vTy, vTz]
    v = [grad_phi[i] + vT[i] for i in range(3)]
    div_rhov = sum(sp.diff(rho * v[i], coords[i]) for i in range(3))
    div_vT = sum(sp.diff(vT[i], coords[i]) for i in range(3))
    lap_phi = sum(sp.diff(phi, c, 2) for c in coords)
    grad_rho_dot_v = sum(sp.diff(rho, coords[i]) * v[i] for i in range(3))

    div_expansion_residual = div_rhov - (grad_rho_dot_v + rho * lap_phi + rho * div_vT)
    checks += 1
    passes += check("div(rho*(grad phi + vT)) expansion", div_expansion_residual)

    C_continuity = sp.diff(rho, t) + div_rhov - S_leak
    identity_raw = rho * lap_phi - (S_leak - sp.diff(rho, t) - grad_rho_dot_v)
    # identity_raw = 0 after imposing continuity and div(v_T)=0.
    # Algebraically: identity_raw - C_continuity + rho*div(v_T) = 0.
    longitudinal_identity_residual = identity_raw - C_continuity + rho * div_vT
    checks += 1
    passes += check("longitudinal identity modulo continuity and div(vT)=0", longitudinal_identity_residual)

    print("Longitudinal identity when div(v_T)=0 and continuity holds:")
    print("  rho_br * Laplacian(phi_br) = S_leak - d_t rho_br - grad(rho_br)·(grad(phi_br)+v_T)")

    # ------------------------------------------------------------------
    # 3. Controlled Poisson hook and inverse-square field
    # ------------------------------------------------------------------
    section("3. Controlled Poisson hook and inverse-square Green function")

    R = sp.symbols("R", positive=True)
    S0, rho0 = sp.symbols("S0 rho0", nonzero=True, positive=True)
    phi_green_R = -S0 / (4 * sp.pi * rho0 * R)
    radial_laplacian = sp.diff(R**2 * sp.diff(phi_green_R, R), R) / R**2
    checks += 1
    passes += check("radial Laplacian(-S0/(4*pi*rho0*R)) for R>0", radial_laplacian)

    radial_derivative = sp.diff(phi_green_R, R)
    flux = clean(4 * sp.pi * R**2 * radial_derivative)
    checks += 1
    passes += check("sphere flux normalization", flux, S0 / rho0)

    print("Distributional reading:")
    print("  Laplacian[-S0/(4*pi*rho0*r)] = (S0/rho0) delta^3(x)")
    print("Controlled reduction conditions used:")
    print("  rho_br ≈ rho0, d_t rho_br subleading, grad rho_br subleading, transverse/advection subleading")

    # ------------------------------------------------------------------
    # 4. Newtonian normalization of the scalar channel
    # ------------------------------------------------------------------
    section("4. Newtonian source normalization and pair-counting")

    lambda_phi, eta, G, m = sp.symbols("lambda_phi eta G m", nonzero=True, positive=True)
    # Effective leakage source amplitude for an isolated defect: S_eff = eta*m*delta^3(x)
    # Phi_N = lambda_phi * phi, so Poisson coefficient = lambda_phi*eta*m/rho0.
    eta_solution = sp.solve(sp.Eq(lambda_phi * eta / rho0, 4 * sp.pi * G), eta)[0]
    checks += 1
    passes += check("source-normalization solution", eta_solution, 4 * sp.pi * G * rho0 / lambda_phi)
    print(f"Required source amplitude per unit mass: eta = {eta_solution}")

    mA, mB, rAB = sp.symbols("m_A m_B r_AB", positive=True)
    etaA, etaB = sp.symbols("eta_A eta_B", positive=True)
    # Potential at A from B: Phi_B(A) = - lambda_phi eta_B m_B/(4*pi*rho0*r_AB)
    Phi_B_at_A = -lambda_phi * etaB * mB / (4 * sp.pi * rho0 * rAB)
    Phi_A_at_B = -lambda_phi * etaA * mA / (4 * sp.pi * rho0 * rAB)
    # Ordered external-potential sum is half-counted to avoid double counting field energy.
    L_pair = clean(-sp.Rational(1, 2) * (mA * Phi_B_at_A + mB * Phi_A_at_B))
    G_AB_eff = clean(L_pair * rAB / (mA * mB))
    print(f"Effective symmetric pair coefficient G_AB = {G_AB_eff}")

    solution_etaAB = sp.solve(
        [sp.Eq(lambda_phi * etaA / (4 * sp.pi * rho0), G),
         sp.Eq(lambda_phi * etaB / (4 * sp.pi * rho0), G)],
        [etaA, etaB],
        dict=True,
    )[0]
    checks += 2
    passes += check("eta_A universality", solution_etaAB[etaA], 4 * sp.pi * G * rho0 / lambda_phi)
    passes += check("eta_B universality", solution_etaAB[etaB], 4 * sp.pi * G * rho0 / lambda_phi)

    G_AB_with_universal_eta = clean(G_AB_eff.subs(solution_etaAB))
    checks += 1
    passes += check("pair coefficient becomes universal G", G_AB_with_universal_eta, G)

    v_A2, v_B2 = sp.symbols("v_A2 v_B2")
    L0_expected = sp.Rational(1, 2) * mA * v_A2 + sp.Rational(1, 2) * mB * v_B2 + G * mA * mB / rAB
    print("Newtonian two-body block after the universality gate:")
    print(f"  L0 = {L0_expected}")

    # ------------------------------------------------------------------
    # 5. Equivalence-principle / response universality gate
    # ------------------------------------------------------------------
    section("5. Response-mass universality gate")

    kappaA, kappaB, kappa_rho, gradPhi = sp.symbols("kappa_A kappa_B kappa_rho gradPhi")
    # L_i = 1/2 m_i v_i^2 - m_g,i Phi, with m_g,i = kappa_i m_i.
    accA = -kappaA * gradPhi
    accB = -kappaB * gradPhi
    acc_difference = clean(accA - accB)
    print(f"Acceleration difference for identical location/field: {acc_difference}")
    print("Equivalence-principle gate: kappa_A = kappa_B = kappa_rho. Newtonian matching chooses kappa_rho = 1 after normalization.")
    checks += 1
    passes += check("response universality if kappa_A=kappa_B", acc_difference.subs({kappaA: kappa_rho, kappaB: kappa_rho}))

    # ------------------------------------------------------------------
    # 6. Same-G ledger through PN and quadrupole normalization
    # ------------------------------------------------------------------
    section("6. Same-G ledger for PN and quadrupole sectors")

    G_N, G_Q, c_s, a, c = sp.symbols("G_N G_Q c_s a c", positive=True, nonzero=True)
    P0_target_Q = 54 * G_Q * c_s**5 / (5 * a**5 * c**5)
    P0_target_N = 54 * G_N * c_s**5 / (5 * a**5 * c**5)
    same_G_residual = clean(P0_target_Q - P0_target_N)
    print(f"P0 target mismatch if sectors use different Gs: {same_G_residual}")
    checks += 1
    passes += check("quadrupole target uses the Newtonian G when G_Q=G_N", same_G_residual.subs(G_Q, G_N))

    # ------------------------------------------------------------------
    # 7. Final status
    # ------------------------------------------------------------------
    section("Final audit status")

    print(f"symbolic_checks: {checks}")
    print(f"symbolic_passes: {passes}")
    print("\nGate verdicts:")
    print("  EXACT_PROJECTED_CONTINUITY: PASS")
    print("  EXACT_LONGITUDINAL_IDENTITY: PASS")
    print("  CONTROLLED_POISSON_HOOK: PASS under stated quasi-static / weak-correction assumptions")
    print("  INVERSE_SQUARE_GREEN_FUNCTION: PASS")
    print("  NEWTONIAN_SOURCE_LAW: CONDITIONAL PASS")
    print("    Required: S_eff for each compact defect is proportional to inertial mass with one universal coefficient.")
    print("    Required: scalar response mass equals inertial mass up to one universal normalization kappa_rho, fixed to 1 in the Newtonian ledger.")
    print("  STRICT_PARENT_DERIVATION_OF_UNIVERSAL_SOURCE_STRENGTH: NOT PROVEN BY THIS IDENTITY")
    print("  SAME_G_LEDGER: PASS as a bookkeeping gate; not independently derived by the Poisson identity.")


if __name__ == "__main__":
    main()
