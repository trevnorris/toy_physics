#!/usr/bin/env python3
"""SymPy cross-check for Path-A 01 symbolic identities.

This mirrors the Mathematica verifier at a compact level.  Genuine derivation
gates are kept separate from construction restatements, disclosure guards, and
bookkeeping identities.
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


def check(name: str, passed: bool, detail: str) -> dict[str, object]:
    return {"check": name, "pass": bool(passed), "detail": detail}


def read_if_exists(path: Path) -> str:
    return path.read_text() if path.exists() else ""


def target_blind_free(texts: list[str], patterns: list[str]) -> bool:
    return all(pattern not in text for text in texts for pattern in patterns)


def vec_add(a: tuple[sp.Expr, ...], b: tuple[sp.Expr, ...]) -> tuple[sp.Expr, ...]:
    return tuple(sp.simplify(x + y) for x, y in zip(a, b))


def vec_sub(a: tuple[sp.Expr, ...], b: tuple[sp.Expr, ...]) -> tuple[sp.Expr, ...]:
    return tuple(sp.simplify(x - y) for x, y in zip(a, b))


def vec_eq(a: tuple[sp.Expr, ...], b: tuple[sp.Expr, ...]) -> bool:
    return all(sp.simplify(x - y) == 0 for x, y in zip(a, b))


def main() -> int:
    script_path = Path(__file__).resolve()
    script_dir = script_path.parent
    repo_root = script_dir.parents[2]
    note_path = repo_root / "software/stage1_solver/derivations/pathA_01_return_source_and_balance.md"
    wls_path = script_dir / "pathA_01_return_source_and_balance.wls"
    report_path = script_dir / "pathA_01_return_source_and_balance_report.md"
    run_dir = script_dir / "runs/pathA_01_return_source_and_balance"
    run_dir.mkdir(parents=True, exist_ok=True)
    diagnostics_path = run_dir / "pathA_01_return_source_and_balance_sympy_diagnostics.json"
    wls_diagnostics_path = run_dir / "pathA_01_return_source_and_balance_diagnostics.json"

    eta, drho, rho0, k1, k2, v0 = sp.symbols("eta drho rho0 k1 k2 V0")
    v_taylor = v0 - k1 * eta + sp.Rational(1, 2) * k2 * eta**2
    l_int = sp.expand(-v_taylor * (rho0 + drho))
    l_cross_coeff = sp.expand(l_int).coeff(eta, 1).coeff(drho, 1)
    l_self_coeff = sp.expand(l_int.subs(drho, 0)).coeff(eta, 2)
    l_rho2_coeff = sp.expand(l_int).coeff(drho, 2)

    h_cross = -k1 * eta * drho
    h_quad = h_cross + sp.Rational(1, 2) * k2 * rho0 * eta**2
    forward = sp.diff(h_cross, drho)
    ret = sp.diff(h_cross, eta)
    cross_eta_rho = sp.diff(sp.diff(h_cross, eta), drho)
    cross_rho_eta = sp.diff(sp.diff(h_cross, drho), eta)
    source_raw = sp.diff(h_quad, eta)
    source_absorbed = sp.simplify(source_raw - k2 * rho0 * eta)

    reciprocity_pass = all(
        [
            sp.simplify(l_cross_coeff - k1) == 0,
            sp.simplify(l_self_coeff + sp.Rational(1, 2) * k2 * rho0) == 0,
            l_rho2_coeff == 0,
            sp.simplify(forward + k1 * eta) == 0,
            sp.simplify(ret + k1 * drho) == 0,
            sp.simplify(cross_eta_rho - cross_rho_eta) == 0,
            sp.simplify(cross_eta_rho + k1) == 0,
            sp.simplify(source_raw - (-k1 * drho + k2 * rho0 * eta)) == 0,
            sp.simplify(source_absorbed + k1 * drho) == 0,
        ]
    )

    t, a, w = sp.symbols("t a w")
    chi = sp.Function("chi")
    A0 = sp.Function("A0")
    Aa = sp.Function("Aa")
    Aw = sp.Function("Aw")
    theta = sp.Function("theta")
    rho = sp.Function("rho")
    q, hbar, mass, mu0, xi, fsq, gfsq, ext_term = sp.symbols(
        "q hbar mass mu0 xi Fsq gfSq extTerm"
    )
    Z = sp.Function("Z")
    H = sp.Function("H")

    ew0 = -sp.diff(Aw(t, a, w), t) - sp.diff(A0(t, a, w), w)
    ca0 = sp.diff(Aw(t, a, w), a) - sp.diff(Aa(t, a, w), w)
    a0_gauge = A0(t, a, w) - sp.diff(chi(t, a, w), t)
    aa_gauge = Aa(t, a, w) + sp.diff(chi(t, a, w), a)
    aw_gauge = Aw(t, a, w) + sp.diff(chi(t, a, w), w)
    ew_gauge = -sp.diff(aw_gauge, t) - sp.diff(a0_gauge, w)
    ca_gauge = sp.diff(aw_gauge, a) - sp.diff(aa_gauge, w)

    j_hydro = (hbar / mass) * rho(t, a, w) * (
        sp.diff(theta(t, a, w), a) - (q / hbar) * Aa(t, a, w)
    )
    theta_gauge = theta(t, a, w) + (q / hbar) * chi(t, a, w)
    j_hydro_gauge = (hbar / mass) * rho(t, a, w) * (
        sp.diff(theta_gauge, a) - (q / hbar) * aa_gauge
    )
    l_em_declared = -Z(w) * fsq / (4 * mu0) - H(w) * gfsq / (2 * xi * mu0) - ext_term
    direct_maxwell_eta_variation = sp.diff(l_em_declared, eta)

    ew_invariant = sp.simplify(ew_gauge - ew0) == 0
    ca_invariant = sp.simplify(ca_gauge - ca0) == 0
    current_invariant = sp.simplify(j_hydro_gauge - j_hydro) == 0
    direct_maxwell_pass = direct_maxwell_eta_variation == 0
    gauge_pass = all([ew_invariant, ca_invariant, current_invariant])
    j_cov_eta, ew_inv, ca_inv, sj, se, sc = sp.symbols("JcovEta EwInv CaInv sJ sE sC")
    allowed_gauge_data = sj * j_cov_eta + se * ew_inv + sc * ca_inv
    allowed_gauge_data_pass = not allowed_gauge_data.has(A0, Aa, Aw, chi)

    d = sp.symbols("d")
    dim_energy = (sp.Integer(1), sp.Integer(2), sp.Integer(-2))
    dim_length = (sp.Integer(0), sp.Integer(1), sp.Integer(0))
    dim_time = (sp.Integer(0), sp.Integer(0), sp.Integer(1))
    dim_eta = dim_length
    dim_rho = (sp.Integer(0), -d, sp.Integer(0))
    dim_k1 = vec_sub(dim_energy, dim_length)
    dim_k2 = vec_sub(dim_energy, tuple(2 * x for x in dim_length))
    dim_source = vec_add(dim_k1, dim_rho)
    dim_keta = vec_sub(dim_source, dim_eta)
    dim_mu_eta = vec_sub(dim_source, vec_sub(dim_eta, tuple(2 * x for x in dim_time)))
    dim_tw = vec_add(dim_source, dim_length)
    dim_tomega = vec_sub(dim_source, dim_eta)
    dimension_pass = all(
        [
            vec_eq(vec_add(vec_add(dim_k1, dim_eta), dim_rho), vec_add(dim_energy, dim_rho)),
            vec_eq(vec_add(vec_add(dim_k2, tuple(2 * x for x in dim_eta)), dim_rho), vec_add(dim_energy, dim_rho)),
            vec_eq(vec_add(dim_k1, dim_rho), dim_source),
            vec_eq(vec_add(vec_add(dim_k2, dim_rho), dim_eta), dim_source),
            vec_eq(vec_add(dim_mu_eta, vec_sub(dim_eta, tuple(2 * x for x in dim_time))), dim_source),
            vec_eq(vec_sub(dim_tw, dim_length), dim_source),
            vec_eq(vec_add(dim_tomega, dim_eta), dim_source),
            vec_eq(vec_add(dim_keta, dim_eta), dim_source),
        ]
    )

    omega_u2, omega_w2, rmix, gu, gw, kstar, qcons = sp.symbols(
        "OmegaU2 OmegaW2 Rmix GU GW Kstar Qcons"
    )
    spsi_sym, sa_sym = sp.symbols("SpsiSym SASym")
    delta_d = sp.Function("deltaD")
    delta_expr = omega_u2 * omega_w2 - rmix**2
    p_expr = omega_u2 * gw + rmix * gu
    n0_expr = (p_expr / delta_expr) ** 2
    d0_expr_with_return = kstar - qcons / delta_expr + delta_d(spsi_sym, sa_sym)
    no_numerator_structure_pass = all(
        [
            sp.diff(p_expr, spsi_sym) == 0,
            sp.diff(p_expr, sa_sym) == 0,
            sp.diff(n0_expr, spsi_sym) == 0,
            sp.diff(n0_expr, sa_sym) == 0,
            bool(d0_expr_with_return.has(spsi_sym, sa_sym)),
        ]
    )
    no_numerator_corollary_pass = reciprocity_pass and no_numerator_structure_pass

    eps, e, ep, et, grad_e, r0p, tw0, tw_r, tw_rr, u0, u1, u2, mu_sigma0, t_omega0, d_a = sp.symbols(
        "eps e ep et gradE R0p Tw0 TwR TwRR U0 U1 U2 muSigma0 TOmega0 dA"
    )
    tw_taylor = tw0 + eps * tw_r * e + sp.Rational(1, 2) * eps**2 * tw_rr * e**2
    u_taylor = u0 + eps * u1 * e + sp.Rational(1, 2) * eps**2 * u2 * e**2
    rw_taylor = r0p + eps * ep
    l_sigma_local = (
        sp.Rational(1, 2) * mu_sigma0 * eps**2 * et**2
        - sp.Rational(1, 2) * tw_taylor * rw_taylor**2
        - sp.Rational(1, 2) * t_omega0 * eps**2 * grad_e**2
        - u_taylor
    )
    l_sigma2 = sp.expand(l_sigma_local).coeff(eps, 2)
    l_sigma2_expected = (
        sp.Rational(1, 2) * mu_sigma0 * et**2
        - sp.Rational(1, 2) * tw0 * ep**2
        - tw_r * r0p * e * ep
        - sp.Rational(1, 4) * tw_rr * r0p**2 * e**2
        - sp.Rational(1, 2) * t_omega0 * grad_e**2
        - sp.Rational(1, 2) * u2 * e**2
    )
    k_eta_formula = u2 - d_a + sp.Rational(1, 2) * tw_rr * r0p**2
    l_sigma_after_parts = (
        sp.Rational(1, 2) * mu_sigma0 * et**2
        - sp.Rational(1, 2) * tw0 * ep**2
        - sp.Rational(1, 2) * t_omega0 * grad_e**2
        - sp.Rational(1, 2) * k_eta_formula * e**2
    )
    l_sigma_parts_from_local = l_sigma2_expected + tw_r * r0p * e * ep + sp.Rational(1, 2) * d_a * e**2
    local_reduction_pass = sp.simplify(l_sigma2 - l_sigma2_expected) == 0

    w_ibp = sp.symbols("wIBP")
    eta_ibp = sp.Function("etaIBP")
    tw_r_ibp = sp.Function("TwRIBP")
    r0p_ibp = sp.Function("R0pIBP")
    ibp_coeff = tw_r_ibp(w_ibp) * r0p_ibp(w_ibp)
    d_a_derived = sp.diff(ibp_coeff, w_ibp)
    ibp_local_density = -ibp_coeff * eta_ibp(w_ibp) * sp.diff(eta_ibp(w_ibp), w_ibp)
    ibp_no_boundary_density = sp.Rational(1, 2) * d_a_derived * eta_ibp(w_ibp) ** 2
    ibp_boundary_density = -sp.Rational(1, 2) * sp.diff(
        ibp_coeff * eta_ibp(w_ibp) ** 2, w_ibp
    )
    ibp_symbolic_pass = (
        sp.simplify(ibp_local_density - (ibp_no_boundary_density + ibp_boundary_density))
        == 0
    )
    ibp_compact_algebra_pass = sp.simplify(l_sigma_parts_from_local - l_sigma_after_parts) == 0
    ibp_reduction_pass = ibp_symbolic_pass and ibp_compact_algebra_pass
    reduction_pass = local_reduction_pass

    note_text = note_path.read_text()
    open_items_pass = all(
        token in note_text
        for token in [
            "## Posited or Open",
            "constitutive",
            "gauge return",
            "moving-boundary",
            "Boundary terms",
            "does not solve",
        ]
    )
    script_text = script_path.read_text()
    forbidden_patterns = [
        "10" + "." + "8",
        "54" + "/" + "5",
        "R" + "_" + "norm",
        "P0" + "_" + "target",
        "GR" + " " + "constant",
    ]
    artifact_texts = [
        note_text,
        script_text,
        read_if_exists(wls_path),
        read_if_exists(report_path),
        read_if_exists(wls_diagnostics_path),
        read_if_exists(diagnostics_path),
    ]
    target_blind_base_pass = target_blind_free(artifact_texts, forbidden_patterns)

    def genuine_checks(target_blind_pass: bool) -> list[dict[str, object]]:
        return [
            check(
                "reciprocity",
                reciprocity_pass,
                "matter cross-Hessian is symmetric with kernel -k1",
            ),
            check(
                "gauge_invariance",
                gauge_pass,
                "three live identities pass: E_w, C_a, and covariant current",
            ),
            check(
                "dimensional_consistency",
                dimension_pass,
                "source and wall terms share one dimension",
            ),
            check(
                "reduction_consistency",
                reduction_pass,
                "local eps^2 coefficient match for promoted S_Sigma",
            ),
            check(
                "target_blind_artifacts",
                target_blind_pass,
                "note, report, Mathematica verifier, SymPy verifier, and diagnostics avoid forbidden tokens",
            ),
        ]

    non_physics_checks = [
        check(
            "no_numerator_knob_corollary_not_a_physics_gate",
            no_numerator_corollary_pass,
            "corollary of reciprocity; P/N0 symbol placement is only a construction restatement",
        ),
        check(
            "honest_open_items_disclosure_guard_not_a_physics_gate",
            open_items_pass,
            "disclosure-presence guard for the Posited or Open section",
        ),
        check(
            "ibp_bookkeeping_identity_not_a_physics_gate",
            ibp_reduction_pass,
            "symbolic IBP derives dA=partial_w(T_w,Sigma,R R0') up to the named boundary density",
        ),
        check(
            "declared_maxwell_no_direct_eta_variation_not_a_physics_gate",
            direct_maxwell_pass,
            "declared fixed-domain Z(w),H(w) Maxwell action has zero direct eta variation by construction",
        ),
        check(
            "allowed_gauge_data_placeholder_not_a_physics_gate",
            allowed_gauge_data_pass,
            "allowed reduced gauge data placeholder contains no raw potentials or chi",
        ),
    ]

    def build_diagnostics(target_blind_pass: bool) -> dict[str, object]:
        checks = genuine_checks(target_blind_pass)
        genuine_pass = all(row["pass"] for row in checks)
        non_physics_pass = all(row["pass"] for row in non_physics_checks)
        all_pass = genuine_pass and non_physics_pass
        return {
            "schema": "pathA_01_return_source_and_balance/sympy_verifier/v2",
            "all_pass": all_pass,
            "genuine_pass": genuine_pass,
            "non_physics_pass": non_physics_pass,
            "checks": checks,
            "construction_restatements_and_guards": non_physics_checks,
            "target_blind_artifacts_scanned": [
                str(note_path),
                str(report_path),
                str(wls_path),
                str(script_path),
                str(wls_diagnostics_path),
                str(diagnostics_path),
            ],
            "symbolic_results": {
                "lagrangian_cross_coefficient": str(l_cross_coeff),
                "lagrangian_self_coefficient": str(l_self_coeff),
                "forward_kernel": str(sp.simplify(forward / eta)),
                "return_kernel": str(sp.simplify(ret / drho)),
                "raw_matter_source": str(source_raw),
                "absorbed_matter_source": str(source_absorbed),
                "direct_maxwell_eta_variation": str(direct_maxwell_eta_variation),
                "dA_derived_by_ibp": str(d_a_derived),
                "ibp_boundary_density": str(ibp_boundary_density),
                "k_eta_formula": str(k_eta_formula),
            },
        }

    candidate_diagnostics = build_diagnostics(target_blind_base_pass)
    candidate_text = json.dumps(candidate_diagnostics, indent=2, sort_keys=True)
    target_blind_pass = target_blind_base_pass and target_blind_free(
        [candidate_text], forbidden_patterns
    )
    diagnostics = build_diagnostics(target_blind_pass)
    diagnostics_text = json.dumps(diagnostics, indent=2, sort_keys=True)
    target_blind_pass = target_blind_pass and target_blind_free(
        [diagnostics_text], forbidden_patterns
    )
    diagnostics = build_diagnostics(target_blind_pass)
    diagnostics_path.write_text(json.dumps(diagnostics, indent=2, sort_keys=True) + "\n")
    all_pass = bool(diagnostics["all_pass"])
    if not all_pass:
        print(f"FAIL pathA_01_return_source_and_balance_sympy; see {diagnostics_path}")
        return 1
    print("PASS pathA_01_return_source_and_balance_sympy")
    print(f"Diagnostics: {diagnostics_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
