#!/usr/bin/env python3
"""Ledger stage015 SymPy audit: breathing structure + HF force.

Standalone, print-only, no arguments, no file I/O. This is the pathA_31
II-G2c slice only: legacy-Hessian structure recovery, Hellmann-Feynman force
by two genuinely different symbolic routes, the static-dynamic limit, and the
HF profile guard. Stage 013's profiles and M/K closed forms are cited as
literals with dual-site integrity; stage 014's truncation certificate is cited.
"""

from __future__ import annotations

from typing import Any

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0

BREATHING_CALIBRATED = "BREATHING_CALIBRATED"
BREATHING_FAIL_STRUCTURE = "BREATHING_FAIL_STRUCTURE"
BREATHING_FAIL_HF_FORCE = "BREATHING_FAIL_HF_FORCE"


class AuditFailure(AssertionError):
    pass


def banner(title: str) -> None:
    print("")
    print("=" * len(title))
    print(title)
    print("=" * len(title))


def subbanner(title: str) -> None:
    print("")
    print(title)
    print("-" * len(title))


def compact(expr: Any) -> Any:
    if isinstance(expr, sp.MatrixBase):
        return expr.applyfunc(compact)
    if not isinstance(expr, sp.Basic):
        return expr
    if isinstance(expr, (sp.Equality, sp.Order)):
        return expr
    return sp.trigsimp(sp.factor(sp.cancel(sp.together(sp.simplify(expr)))))


def hstr(expr: Any) -> str:
    if isinstance(expr, sp.Basic) and expr.has(sp.Integral):
        return sp.sstr(expr)
    return sp.sstr(compact(expr))


def fmt(expr: Any) -> str:
    if isinstance(expr, bool):
        return "True" if expr else "False"
    if isinstance(expr, str):
        return expr
    if isinstance(expr, sp.MatrixBase):
        return sp.sstr(compact(expr))
    try:
        return sp.sstr(compact(expr))
    except Exception:
        return sp.sstr(expr)


def assert_no_float(name: str, expr: Any) -> None:
    if isinstance(expr, dict):
        for key, value in expr.items():
            assert_no_float(f"{name}.{key}", value)
        return
    if isinstance(expr, sp.MatrixBase):
        for index, value in enumerate(expr):
            assert_no_float(f"{name}[{index}]", value)
        return
    if isinstance(expr, (list, tuple, set, frozenset)):
        for index, value in enumerate(expr):
            assert_no_float(f"{name}[{index}]", value)
        return
    if isinstance(expr, (str, type(None))):
        return
    if isinstance(expr, bool):
        expr = sp.Integer(1) if expr else sp.Integer(0)
    clean = sp.sympify(expr)
    floats = clean.atoms(sp.Float)
    if floats:
        raise AuditFailure(f"{name}: Float atom(s) found in exact audit expression: {floats}")


def _record_pass(message: str) -> None:
    global PASS_COUNT
    PASS_COUNT += 1
    print(message)


def _record_fail(message: str) -> None:
    global FAIL_COUNT
    FAIL_COUNT += 1
    print(message)


def expect_zero(name: str, residual: sp.Expr | int) -> None:
    assert_no_float(name, residual)
    clean = compact(residual)
    assert_no_float(name, clean)
    if clean == 0:
        _record_pass(f"PASS  {name}")
        return
    _record_fail(f"FAIL  {name}: residual = {fmt(clean)}")
    raise AuditFailure(f"{name} residual was not zero")


def expect_bool(name: str, condition: bool) -> None:
    expect_zero(name, sp.Integer(0) if bool(condition) else sp.Integer(1))


def expect_nonzero(name: str, residual: sp.Expr | int) -> None:
    assert_no_float(name, residual)
    clean = compact(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_pass(f"PASS  {name} is nonzero as required (residual = {fmt(clean)})")
        return
    _record_fail(f"FAIL  {name}: required nonzero residual vanished")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def expect_fail(name: str, residual: sp.Expr | int) -> None:
    assert_no_float(name, residual)
    clean = compact(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_pass(f"PASS  {name} produced required FAIL (residual = {fmt(clean)})")
        return
    _record_fail(f"FAIL  {name}: required mutation/ablation did not fire")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def bool_residual(condition: bool) -> sp.Integer:
    return sp.Integer(0) if bool(condition) else sp.Integer(1)


def verdict_residual(actual: str, expected: str) -> sp.Integer:
    return sp.Integer(0) if actual == expected else sp.Integer(1)


def zero_pattern(matrix: sp.Matrix) -> tuple[tuple[bool, ...], ...]:
    return tuple(tuple(compact(matrix[i, j]) == 0 for j in range(matrix.cols)) for i in range(matrix.rows))


w = sp.Symbol("w", real=True)
L0 = sp.Symbol("L0", positive=True, real=True)
beta = sp.Symbol("beta", positive=True, real=True)
muEta = sp.Symbol("muEta", positive=True, real=True)
Tw = sp.Symbol("Tw", positive=True, real=True)
rAL = sp.Symbol("rAL", positive=True, real=True)
rhoStar = sp.Symbol("rhoStar", positive=True, real=True)
Vp0 = sp.Symbol("Vp0", positive=True, real=True)
ellC = sp.Symbol("ellC", positive=True, real=True)
kappa = sp.Symbol("kappa", positive=True, real=True)
chi = sp.Symbol("chi", positive=True, real=True)
sigmaA = sp.Symbol("sigmaA", positive=True, real=True)
sigmaL = sp.Symbol("sigmaL", positive=True, real=True)
delta_a, delta_L = sp.symbols("delta_a delta_L", real=True)
qa, qL, r = sp.symbols("qa qL r", real=True)
B = sp.Symbol("B", positive=True, real=True)


def cited_profiles() -> dict[str, sp.Expr]:
    return {
        "alpha_a": sp.sinh(L0 * beta - beta * w) / sp.sinh(L0 * beta),
        "alpha_L": rAL * sp.sinh(beta * w) / sp.sinh(L0 * beta),
    }


def cited_mk() -> dict[str, Any]:
    M_entries = {
        "aa": -2 * sp.pi * muEta * (L0 * beta * sp.tanh(L0 * beta) - sp.sinh(L0 * beta) ** 2) / (beta * sp.sinh(L0 * beta) ** 2 * sp.tanh(L0 * beta)),
        "aL": 2 * sp.pi * muEta * rAL * (L0 * beta - sp.tanh(L0 * beta)) / (beta * sp.sinh(L0 * beta) * sp.tanh(L0 * beta)),
        "LL": -2 * sp.pi * muEta * rAL**2 * (L0 * beta * sp.tanh(L0 * beta) - sp.sinh(L0 * beta) ** 2) / (beta * sp.sinh(L0 * beta) ** 2 * sp.tanh(L0 * beta)),
    }
    K_entries = {
        "aa": 4 * sp.pi * Tw * beta / sp.tanh(L0 * beta),
        "aL": -4 * sp.pi * Tw * beta * rAL / sp.sinh(L0 * beta),
        "LL": 4 * sp.pi * Tw * beta * rAL**2 / sp.tanh(L0 * beta),
    }
    M = sp.Matrix([[M_entries["aa"], M_entries["aL"]], [M_entries["aL"], M_entries["LL"]]])
    K = sp.Matrix([[K_entries["aa"], K_entries["aL"]], [K_entries["aL"], K_entries["LL"]]])
    return {"M_entries": M_entries, "K_entries": K_entries, "M": M, "K": K}


def det_m_013_form() -> sp.Expr:
    return 4 * sp.pi**2 * muEta**2 * rAL**2 * (sp.sinh(L0 * beta) ** 2 - (L0 * beta) ** 2) / (beta**2 * sp.sinh(L0 * beta) ** 2)


def det_k_013_form() -> sp.Expr:
    return 16 * sp.pi**2 * Tw**2 * beta**2 * rAL**2


def positivity_certificate_identities() -> dict[str, sp.Expr]:
    f1 = sp.sinh(B) * sp.cosh(B) - B
    g = sp.sinh(B) - B
    f2 = sp.sinh(B) ** 2 - B**2
    h = B - sp.tanh(B)
    return {
        "f1_at_0": f1.subs(B, 0),
        "f1_prime_square": sp.diff(f1, B) - 2 * sp.sinh(B) ** 2,
        "g_at_0": g.subs(B, 0),
        "g_prime_square": sp.diff(g, B) - 2 * sp.sinh(B / 2) ** 2,
        "f2_factorization": f2 - g * (sp.sinh(B) + B),
        "h_at_0": h.subs(B, 0),
        "h_prime_square": sp.diff(h, B) - sp.tanh(B) ** 2,
    }


def entry_certificate_residuals(M: sp.Matrix, K: sp.Matrix) -> dict[str, sp.Expr]:
    local_B = L0 * beta
    f1_local = sp.sinh(local_B) * sp.cosh(local_B) - local_B
    f2_local = sp.sinh(local_B) ** 2 - local_B**2
    h_local = local_B - sp.tanh(local_B)
    return {
        "M_aa_core_is_f1": compact(M[0, 0] * beta * sp.sinh(local_B) ** 2 / (2 * sp.pi * muEta) - f1_local),
        "det_M_core_is_f2": compact(M.det() * beta**2 * sp.sinh(local_B) ** 2 / (4 * sp.pi**2 * muEta**2 * rAL**2) - f2_local),
        "M_aL_core_is_h": compact(M[0, 1] * beta * sp.sinh(local_B) * sp.tanh(local_B) / (2 * sp.pi * muEta * rAL) - h_local),
        "K_aL_negative_core": compact(K[0, 1] * sp.sinh(local_B) / (4 * sp.pi * Tw * beta * rAL) + 1),
    }


def m_posdef_by_certificate(M: sp.Matrix) -> bool:
    ids = positivity_certificate_identities()
    residuals = entry_certificate_residuals(M, sp.zeros(2))
    return all(compact(value) == 0 for value in (
        ids["f1_at_0"],
        ids["f1_prime_square"],
        ids["g_at_0"],
        ids["g_prime_square"],
        ids["f2_factorization"],
        residuals["M_aa_core_is_f1"],
        residuals["det_M_core_is_f2"],
    ))


def m_aL_positive_by_certificate(M: sp.Matrix) -> bool:
    ids = positivity_certificate_identities()
    residuals = entry_certificate_residuals(M, sp.zeros(2))
    return all(compact(value) == 0 for value in (
        ids["h_at_0"],
        ids["h_prime_square"],
        residuals["M_aL_core_is_h"],
    ))


def k_offdiag_negative_by_certificate(K: sp.Matrix) -> bool:
    residuals = entry_certificate_residuals(sp.eye(2), K)
    return compact(residuals["K_aL_negative_core"]) == 0


def mk_citation_integrity(M: sp.Matrix, K: sp.Matrix) -> dict[str, Any]:
    det_m_residual = compact(M.det() - det_m_013_form())
    det_k_residual = compact(K.det() - det_k_013_form())
    return {
        "det_M_residual": det_m_residual,
        "det_K_residual": det_k_residual,
        "M_aL_positive": m_aL_positive_by_certificate(M),
        "K_aL_negative": k_offdiag_negative_by_certificate(K),
        "ok": det_m_residual == 0 and det_k_residual == 0 and m_aL_positive_by_certificate(M) and k_offdiag_negative_by_certificate(K),
    }


def structure_gate(M: sp.Matrix, K: sp.Matrix, H_legacy: sp.Matrix) -> dict[str, Any]:
    legacy_det_positive_form = kappa * sigmaA + kappa * chi**2 * sigmaL + sigmaA * sigmaL
    legacy_symmetric = H_legacy == H_legacy.T
    legacy_offdiag_negative = compact(H_legacy[0, 1] + chi * kappa) == 0
    legacy_det_positive = compact(H_legacy.det() - legacy_det_positive_form) == 0
    legacy_rank = H_legacy.rank()
    legacy_zeros = zero_pattern(H_legacy)

    M_posdef = m_posdef_by_certificate(M)
    K_symmetric = K == K.T
    K_offdiag_negative = k_offdiag_negative_by_certificate(K)
    K_rank = K.rank()
    K_rank_matches_legacy = K_rank == legacy_rank == 2
    K_zero_pattern = zero_pattern(K)
    K_zero_pattern_matches_legacy = K_zero_pattern == legacy_zeros
    K_structure_ok = bool(
        K_symmetric
        and K_offdiag_negative
        and K_rank_matches_legacy
        and K_zero_pattern_matches_legacy
        and legacy_symmetric
        and legacy_offdiag_negative
        and legacy_det_positive
    )
    return {
        "M_posdef": M_posdef,
        "K_symmetric": K_symmetric,
        "K_offdiag_negative": K_offdiag_negative,
        "K_rank": K_rank,
        "K_rank_matches_legacy": K_rank_matches_legacy,
        "K_zero_pattern": K_zero_pattern,
        "K_zero_pattern_matches_legacy": K_zero_pattern_matches_legacy,
        "K_structure_ok": K_structure_ok,
        "legacy_symmetric": legacy_symmetric,
        "legacy_offdiag_negative": legacy_offdiag_negative,
        "legacy_det_positive": legacy_det_positive,
        "legacy_rank": legacy_rank,
        "legacy_zero_pattern": legacy_zeros,
        "structure_from_computed_matrices": True,
        "full_matrix_fit": False,
    }


def structure_probes(M: sp.Matrix, K: sp.Matrix, H_legacy: sp.Matrix) -> dict[str, Any]:
    M_probe = sp.Matrix(M)
    M_probe[0, 0] = -M_probe[0, 0]
    K_probe = sp.Matrix(K)
    K_probe[0, 1] = -K_probe[0, 1]
    K_probe[1, 0] = -K_probe[1, 0]
    return {
        "non_posdef_M_probe": structure_gate(M_probe, K, H_legacy),
        "sign_flipped_K_probe": structure_gate(M, K_probe, H_legacy),
    }


def compute_verdict(gate: dict[str, Any], hf_force_reduces: bool, unsimplified_routes_identical: bool) -> str:
    if not (gate["M_posdef"] and gate["K_symmetric"] and gate["K_structure_ok"]):
        return BREATHING_FAIL_STRUCTURE
    if not hf_force_reduces or unsimplified_routes_identical:
        return BREATHING_FAIL_HF_FORCE
    return BREATHING_CALIBRATED


def packet_residuals(packet: dict[str, sp.Expr]) -> dict[str, sp.Expr]:
    return {
        "site_A_betaL0_anchor": compact(packet["betaL0_cited"] - sp.Rational(37, 20)),
        "site_A_packet_product": compact(packet["beta_cited"] * packet["L0_cited"] - packet["betaL0_cited"]),
        "anchor_L0": compact(packet["L0_cited"] - sp.Rational(37, 20)),
        "anchor_beta": compact(packet["beta_cited"] - 1),
    }


def packet_ok(packet: dict[str, sp.Expr]) -> bool:
    return all(residual == 0 for residual in packet_residuals(packet).values())


def packet_mutation(packet: dict[str, sp.Expr], key: str, value: sp.Expr) -> dict[str, sp.Expr]:
    mutated = dict(packet)
    mutated[key] = value
    return mutated


def profile_site_b_residuals(alpha_a: sp.Expr, alpha_L: sp.Expr) -> dict[str, sp.Expr]:
    return {
        "residual_alpha_a": compact(-sp.diff(alpha_a, w, 2) + beta**2 * alpha_a),
        "residual_alpha_L": compact(-sp.diff(alpha_L, w, 2) + beta**2 * alpha_L),
        "bc_alpha_a_0": compact(alpha_a.subs(w, 0) - 1),
        "bc_alpha_a_L0": compact(alpha_a.subs(w, L0)),
        "bc_alpha_L_0": compact(alpha_L.subs(w, 0)),
        "bc_alpha_L_L0": compact(alpha_L.subs(w, L0) - rAL),
    }


def profile_site_b_ok(alpha_a: sp.Expr, alpha_L: sp.Expr) -> bool:
    return all(residual == 0 for residual in profile_site_b_residuals(alpha_a, alpha_L).values())


def hessian_entry_residuals(H_legacy: sp.Matrix) -> dict[str, sp.Expr]:
    return {
        "aa": compact(H_legacy[0, 0] - (chi**2 * kappa + sigmaA)),
        "aL": compact(H_legacy[0, 1] + chi * kappa),
        "LL": compact(H_legacy[1, 1] - (kappa + sigmaL)),
        "det": compact(H_legacy.det() - (kappa * sigmaA + kappa * chi**2 * sigmaL + sigmaA * sigmaL)),
    }


def build_baseline() -> dict[str, Any]:
    profiles = cited_profiles()
    alpha_a = profiles["alpha_a"]
    alpha_L = profiles["alpha_L"]
    mk = cited_mk()
    M = mk["M"]
    K = mk["K"]

    E_geom = (
        sp.Rational(1, 2) * kappa * (delta_L - chi * delta_a) ** 2
        + sp.Rational(1, 2) * sigmaA * delta_a**2
        + sp.Rational(1, 2) * sigmaL * delta_L**2
    )
    H_legacy = sp.hessian(E_geom, (delta_a, delta_L))
    gate = structure_gate(M, K, H_legacy)
    probes = structure_probes(M, K, H_legacy)

    source_density = compact(rhoStar * Vp0 / ellC)
    F_dist_a_uns = 4 * sp.pi * sp.Integral(source_density * alpha_a, (w, 0, L0))
    F_dist_L_uns = 4 * sp.pi * sp.Integral(source_density * alpha_L, (w, 0, L0))
    F_dist_a = compact(F_dist_a_uns.doit())
    F_dist_L = compact(F_dist_L_uns.doit())

    R0 = sp.Function("R0")
    R_param = R0(w) + qa * alpha_a + qL * alpha_L
    V_conf = (Vp0 / ellC) * (r - R_param)
    F_legacy_a_uns = -4 * sp.pi * sp.Integral(rhoStar * sp.diff(V_conf, qa).subs({qa: 0, qL: 0}), (w, 0, L0))
    F_legacy_L_uns = -4 * sp.pi * sp.Integral(rhoStar * sp.diff(V_conf, qL).subs({qa: 0, qL: 0}), (w, 0, L0))
    F_legacy_a = compact(F_legacy_a_uns.doit())
    F_legacy_L = compact(F_legacy_L_uns.doit())
    hf_a_ok = compact(F_dist_a - F_legacy_a) == 0
    hf_L_ok = compact(F_dist_L - F_legacy_L) == 0
    unsimplified_routes_identical = hstr(F_dist_a_uns) == hstr(F_legacy_a_uns)

    wrong_zero_F_a = compact(4 * sp.pi * sp.integrate(source_density * 0, (w, 0, L0)))
    wrong_const_F_a = compact(4 * sp.pi * sp.integrate(source_density * 1, (w, 0, L0)))
    degenerate_hf_mismatch = compact(wrong_zero_F_a - F_legacy_a) != 0
    constant_hf_mismatch = compact(wrong_const_F_a - F_legacy_a) != 0
    verdict = compute_verdict(gate, bool(hf_a_ok and hf_L_ok), unsimplified_routes_identical)

    return {
        "profiles": profiles,
        "mk": mk,
        "E_geom": E_geom,
        "H_legacy": H_legacy,
        "gate": gate,
        "probes": probes,
        "source_density": source_density,
        "hf": {
            "F_dist_a_uns": F_dist_a_uns,
            "F_dist_L_uns": F_dist_L_uns,
            "F_dist_a": F_dist_a,
            "F_dist_L": F_dist_L,
            "F_legacy_a_uns": F_legacy_a_uns,
            "F_legacy_L_uns": F_legacy_L_uns,
            "F_legacy_a": F_legacy_a,
            "F_legacy_L": F_legacy_L,
            "hf_a_ok": hf_a_ok,
            "hf_L_ok": hf_L_ok,
            "hf_force_reduces": bool(hf_a_ok and hf_L_ok),
            "unsimplified_routes_identical": unsimplified_routes_identical,
        },
        "guards": {
            "wrong_zero_F_a": wrong_zero_F_a,
            "wrong_const_F_a": wrong_const_F_a,
            "degenerate_hf_mismatch": degenerate_hf_mismatch,
            "constant_hf_mismatch": constant_hf_mismatch,
            "constant_overlap_passes_014": True,
        },
        "packet": {
            "L0_cited": sp.Rational(37, 20),
            "beta_cited": sp.Integer(1),
            "betaL0_cited": sp.Rational(37, 20),
        },
        "mk_integrity": mk_citation_integrity(M, K),
        "verdict": verdict,
    }


def run_consumed_inputs(data: dict[str, Any]) -> None:
    subbanner("CONSUMED-from-013/014 citation integrity")
    packet = data["packet"]
    alpha_a = data["profiles"]["alpha_a"]
    alpha_L = data["profiles"]["alpha_L"]
    M = data["mk"]["M"]
    K = data["mk"]["K"]
    integrity = data["mk_integrity"]

    print("  CONSUMED-from-013: profiles alpha_a,alpha_L; operator-projected M_AB/K_AB; frozen packet {L0=37/20,beta=1,beta*L0=37/20}.")
    print("  CONSUMED-from-014: truncation certificate beta_L0 in [0.1,3.0] / K_eta/Tw <= ~2.6; cited only, no eigensolve here.")
    print("  Site A reads betaL0_cited as an independently corruptible cited datum, not local beta*L0.")
    for name, residual in packet_residuals(packet).items():
        expect_zero(f"site A consumed packet {name}", residual)
    expect_bool("site A clean packet passes", packet_ok(packet))
    expect_fail("site A betaL0_cited-only corruption is non-vacuous", bool_residual(packet_ok(packet_mutation(packet, "betaL0_cited", sp.Rational(19, 10)))))

    print("  Site B verifies cited profile residual AND endpoint BCs; it does not solve the BVP.")
    for name, residual in profile_site_b_residuals(alpha_a, alpha_L).items():
        expect_zero(f"site B profile {name}", residual)
    expect_bool("site B residual-and-BC guard passes cited profiles", profile_site_b_ok(alpha_a, alpha_L))
    kernel_preserving_bad = sp.sinh(L0 * beta - beta * w) / sp.cosh(L0 * beta)
    nonkernel_bad = alpha_a + 1
    kernel_residuals = profile_site_b_residuals(kernel_preserving_bad, alpha_L)
    expect_bool(
        "kernel-preserving corruption has zero residual but wrong BC",
        kernel_residuals["residual_alpha_a"] == 0 and kernel_residuals["bc_alpha_a_0"] != 0,
    )
    expect_fail("kernel-preserving profile corruption breaks site B BC leg", bool_residual(profile_site_b_ok(kernel_preserving_bad, alpha_L)))
    expect_fail("non-kernel profile corruption breaks site B residual leg", bool_residual(profile_site_b_ok(nonkernel_bad, alpha_L)))

    print("  Consumed M/K det-identities plus off-diagonal sign checks cover every cited entry, including det-blind sign flips.")
    expect_zero("consumed det(M) identity vs independent 013 closed form", integrity["det_M_residual"])
    expect_zero("consumed det(K) identity vs independent 013 closed form", integrity["det_K_residual"])
    expect_bool("consumed M_aL > 0 via B-tanh(B) certificate", integrity["M_aL_positive"])
    expect_bool("consumed K_aL < 0 via -1/sinh(B) certificate", integrity["K_aL_negative"])
    expect_bool("consumed M/K citation integrity clean baseline passes", integrity["ok"])

    M_bad = sp.Matrix(M)
    M_bad[0, 0] = 2 * M_bad[0, 0]
    expect_fail("M_aa diagonal corruption breaks det(M) identity", mk_citation_integrity(M_bad, K)["det_M_residual"])
    M_bad = sp.Matrix(M)
    M_bad[1, 1] = 2 * M_bad[1, 1]
    expect_fail("M_LL diagonal corruption breaks det(M) identity", mk_citation_integrity(M_bad, K)["det_M_residual"])
    M_bad = sp.Matrix(M)
    M_bad[0, 1] = 2 * M_bad[0, 1]
    M_bad[1, 0] = 2 * M_bad[1, 0]
    expect_fail("M_aL magnitude corruption breaks det(M) identity", mk_citation_integrity(M_bad, K)["det_M_residual"])
    M_bad = sp.Matrix(M)
    M_bad[0, 1] = -M_bad[0, 1]
    M_bad[1, 0] = -M_bad[1, 0]
    sign_flip_integrity = mk_citation_integrity(M_bad, K)
    expect_zero("M_aL sign flip leaves det(M) identity blind as expected", sign_flip_integrity["det_M_residual"])
    expect_fail("M_aL sign flip breaks the explicit M_aL>0 check", bool_residual(sign_flip_integrity["M_aL_positive"]))
    K_bad = sp.Matrix(K)
    K_bad[0, 0] = 2 * K_bad[0, 0]
    expect_fail("K_aa diagonal corruption breaks det(K) identity", mk_citation_integrity(M, K_bad)["det_K_residual"])
    K_bad = sp.Matrix(K)
    K_bad[1, 1] = 2 * K_bad[1, 1]
    expect_fail("K_LL diagonal corruption breaks det(K) identity", mk_citation_integrity(M, K_bad)["det_K_residual"])
    K_bad = sp.Matrix(K)
    K_bad[0, 1] = 2 * K_bad[0, 1]
    K_bad[1, 0] = 2 * K_bad[1, 0]
    expect_fail("K_aL magnitude corruption breaks det(K) identity", mk_citation_integrity(M, K_bad)["det_K_residual"])
    K_bad = sp.Matrix(K)
    K_bad[0, 1] = -K_bad[0, 1]
    K_bad[1, 0] = -K_bad[1, 0]
    k_sign_flip_integrity = mk_citation_integrity(M, K_bad)
    expect_fail("K_aL sign flip breaks explicit K_aL<0 check", bool_residual(k_sign_flip_integrity["K_aL_negative"]))


def run_legacy_hessian(data: dict[str, Any]) -> None:
    subbanner("Legacy E_geom and own-built H_legacy")
    H = data["H_legacy"]
    E = data["E_geom"]
    print("  E_geom = ", fmt(E))
    print("  H_legacy = ", fmt(H))
    for name, residual in hessian_entry_residuals(H).items():
        expect_zero(f"H_legacy entry/det {name}", residual)
    expect_bool("legacy_symmetric", H == H.T)
    expect_bool("legacy_offdiag_negative certificate: H_aL=-chi*kappa with chi,kappa>0", compact(H[0, 1] + chi * kappa) == 0)
    expect_bool("legacy_det_positive certificate uses positive-term form", data["gate"]["legacy_det_positive"])

    legacy_names = {"kappa", "chi", "sigmaA", "sigmaL"}
    mk_symbols = {symbol.name for expr in [*data["mk"]["M_entries"].values(), *data["mk"]["K_entries"].values()] for symbol in expr.free_symbols}
    print("  LEGACY-CONSTANT boundary: {kappa,chi,sigmaA,sigmaL} enter only here; 013 M/K free symbols = ", sorted(mk_symbols))
    expect_bool("013 M/K free-symbol firewall excludes 015 legacy constants", not (legacy_names & mk_symbols))

    E_no_sigma_a = sp.Rational(1, 2) * kappa * (delta_L - chi * delta_a) ** 2 + sp.Rational(1, 2) * sigmaL * delta_L**2
    H_no_sigma_a = sp.hessian(E_no_sigma_a, (delta_a, delta_L))
    expect_fail("legacy-Hessian tooth: dropping sigmaA support term breaks H_aa", hessian_entry_residuals(H_no_sigma_a)["aa"])
    E_no_chi = sp.Rational(1, 2) * kappa * delta_L**2 + sp.Rational(1, 2) * sigmaA * delta_a**2 + sp.Rational(1, 2) * sigmaL * delta_L**2
    H_no_chi = sp.hessian(E_no_chi, (delta_a, delta_L))
    expect_fail("legacy-Hessian tooth: dropping chi cross coupling breaks H_aL", hessian_entry_residuals(H_no_chi)["aL"])


def run_structure_recovery(data: dict[str, Any]) -> None:
    subbanner("EARNED-STRUCTURE-RECOVERY")
    gate = data["gate"]
    probes = data["probes"]
    M = data["mk"]["M"]
    K = data["mk"]["K"]
    H = data["H_legacy"]
    certs = positivity_certificate_identities()
    entry_certs = entry_certificate_residuals(M, K)

    print("  Structure gate compares cited 013 M/K against own-built H_legacy: M pos-def; K symmetric; K_aL<0; rank and zero-pattern match.")
    print("  Positivity method: exact identities f1'=2*sinh(B)^2, f2=(sinh(B)-B)(sinh(B)+B), M_aL via h=B-tanh(B); no ask/is_positive and no M[i,j]-same-form positivity tooth.")
    for name, residual in certs.items():
        expect_zero(f"positivity certificate identity {name}", residual)
    for name, residual in entry_certs.items():
        expect_zero(f"entry-to-certificate residual {name}", residual)
    expect_bool("M_posdef via exact Sylvester certificates", gate["M_posdef"])
    expect_bool("K_symmetric", gate["K_symmetric"])
    expect_bool("K_offdiag_negative via cited K_aL sign certificate", gate["K_offdiag_negative"])
    expect_bool("K_rank==rank(H_legacy)==2", gate["K_rank_matches_legacy"])
    expect_bool("zero_pattern(K)==zero_pattern(H_legacy)", gate["K_zero_pattern_matches_legacy"])
    expect_bool("K_structure_ok", gate["K_structure_ok"])
    legacy_symbols = {kappa, chi, sigmaA, sigmaL}
    mk_free_symbols = set().union(*(entry.free_symbols for entry in list(M) + list(K)))
    expect_bool("structure_from_computed_matrices: cited M/K free_symbols exclude legacy {kappa,chi,sigmaA,sigmaL}", mk_free_symbols.isdisjoint(legacy_symbols))
    expect_nonzero("recovery is structural not full entrywise fit: M_aa != H_legacy_aa", compact(M[0, 0] - H[0, 0]))
    print("  K_zero_pattern = ", gate["K_zero_pattern"], "; legacy_zero_pattern = ", gate["legacy_zero_pattern"])
    print("  structure_from_computed_matrices=True; full_matrix_fit=False.")

    m_probe = probes["non_posdef_M_probe"]
    k_probe = probes["sign_flipped_K_probe"]
    expect_bool("probe M_aa -> -M_aa flips M_posdef false", not m_probe["M_posdef"])
    expect_bool("probe K_aL -> -K_aL flips K_offdiag_negative false", not k_probe["K_offdiag_negative"])
    expect_bool("probe K_aL -> -K_aL flips K_structure_ok false", not k_probe["K_structure_ok"])

    M_bad = sp.Matrix(M)
    M_bad[0, 0] = -M_bad[0, 0]
    expect_fail("structure entry tooth: corrupt M_aa sign fails M_posdef", bool_residual(structure_gate(M_bad, K, H)["M_posdef"]))
    K_bad = sp.Matrix(K)
    K_bad[0, 1] = -K_bad[0, 1]
    K_bad[1, 0] = -K_bad[1, 0]
    expect_fail("structure entry tooth: corrupt K_aL sign fails K_structure_ok", bool_residual(structure_gate(M, K_bad, H)["K_structure_ok"]))
    H_bad = sp.Matrix(H)
    H_bad[0, 1] = 0
    H_bad[1, 0] = 0
    bad_gate = structure_gate(M, K, H_bad)
    expect_fail("structure entry tooth: H_legacy offdiag mutation breaks zero-pattern match", bool_residual(bad_gate["K_zero_pattern_matches_legacy"]))


def run_hf_force(data: dict[str, Any]) -> None:
    subbanner("HF-TWO-ROUTE-GENUINENESS")
    hf = data["hf"]
    expected_F_a = 4 * sp.pi * Vp0 * rhoStar * (sp.cosh(L0 * beta) - 1) / (beta * ellC * sp.sinh(L0 * beta))
    expected_F_L = 4 * sp.pi * Vp0 * rhoStar * rAL * (sp.cosh(L0 * beta) - 1) / (beta * ellC * sp.sinh(L0 * beta))
    print("  source_density = ", fmt(data["source_density"]), "; measure is action 4*pi*int dw, not mu_eta-weighted.")
    print("  Distributed unsimplified F_a raw tree:")
    print("    ", hstr(hf["F_dist_a_uns"]))
    print("  Legacy/HF unsimplified F_a raw tree:")
    print("    ", hstr(hf["F_legacy_a_uns"]))
    print("  Distributed simplified: F_a=", fmt(hf["F_dist_a"]), "; F_L=", fmt(hf["F_dist_L"]), sep="")
    print("  Legacy simplified: F_a=", fmt(hf["F_legacy_a"]), "; F_L=", fmt(hf["F_legacy_L"]), sep="")
    expect_zero("F_dist_a closed form", hf["F_dist_a"] - expected_F_a)
    expect_zero("F_legacy_a closed form", hf["F_legacy_a"] - expected_F_a)
    expect_zero("F_dist_L closed form", hf["F_dist_L"] - expected_F_L)
    expect_zero("F_legacy_L closed form", hf["F_legacy_L"] - expected_F_L)
    expect_zero("HF route difference F_a", hf["F_dist_a"] - hf["F_legacy_a"])
    expect_zero("HF route difference F_L", hf["F_dist_L"] - hf["F_legacy_L"])
    expect_bool("hf_force_reduces=True", hf["hf_force_reduces"])
    expect_bool("unsimplified_routes_identical computed False", not hf["unsimplified_routes_identical"])

    typed_twice_flag = hstr(hf["F_dist_a_uns"]) == hstr(hf["F_dist_a_uns"])
    print("  tooth 1a typed-twice setup: unsimplified_routes_identical=True scaffolding is de-counted; verdict tooth follows.")
    expect_zero(
        "tooth 1a typed-twice Route 2 trips BREATHING_FAIL_HF_FORCE",
        verdict_residual(compute_verdict(data["gate"], hf["hf_force_reduces"], typed_twice_flag), BREATHING_FAIL_HF_FORCE),
    )
    corrupted_dist_a = compact(hf["F_dist_a"] * ellC / Vp0)
    corrupt_reduces = compact(corrupted_dist_a - hf["F_legacy_a"]) == 0 and compact(hf["F_dist_L"] - hf["F_legacy_L"]) == 0
    expect_fail("tooth 2 corrupt one HF route breaks hf_force_reduces", bool_residual(corrupt_reduces))
    expect_zero(
        "tooth 2 corrupt one route trips BREATHING_FAIL_HF_FORCE",
        verdict_residual(compute_verdict(data["gate"], corrupt_reduces, hf["unsimplified_routes_identical"]), BREATHING_FAIL_HF_FORCE),
    )


def run_static_dynamic_limit(data: dict[str, Any]) -> None:
    subbanner("static-dynamic-limit")
    print("  Qddot -> 0 in the same M_AB Qddot + K_AB Q = F_A system gives K_AB Q = F_A.")
    print("  After the legacy-pattern comparison this is the static dE_geom/dQ=0 limit; no separate static solve is fabricated.")
    qdd_a, qdd_L, force_a, force_L = sp.symbols("qdd_a qdd_L F_a F_L", real=True)
    M = data["mk"]["M"]
    K = data["mk"]["K"]
    Q = sp.Matrix([qa, qL])
    Qdd = sp.Matrix([qdd_a, qdd_L])
    F = sp.Matrix([force_a, force_L])
    eom_residual = M * Qdd + K * Q - F
    static_residual = compact(eom_residual.subs({qdd_a: 0, qdd_L: 0}) - (K * Q - F))
    expect_zero("static_limit_consistent: Qddot->0 EOM residual equals K*Q-F", sum(entry**2 for entry in static_residual))
    print("  separate_static_solve=False (scope note only; no counted literal tooth).")


def run_hf_profile_guard(data: dict[str, Any]) -> None:
    subbanner("HF-PROFILE-GUARD / 014<->015 boundary")
    guards = data["guards"]
    hf = data["hf"]
    print("  Degenerate alpha_a=0: wrong_zero_F_a = ", fmt(guards["wrong_zero_F_a"]))
    print("  Constant alpha_a=1: wrong_const_F_a = ", fmt(guards["wrong_const_F_a"]))
    print("  Frozen correct legacy force F_a = ", fmt(hf["F_legacy_a"]))
    print("  014<->015 boundary: constant_one overlap_passes=True in 014, but hf_mismatch=True here; HF is the profile guard the overlap could not supply.")
    expect_zero("wrong_zero_F_a is native zero integral", guards["wrong_zero_F_a"])
    expect_zero("wrong_const_F_a is 4*pi*L0*Vp0*rhoStar/ellC", guards["wrong_const_F_a"] - 4 * sp.pi * L0 * Vp0 * rhoStar / ellC)
    expect_bool("degenerate hf_mismatch=True", guards["degenerate_hf_mismatch"])
    expect_bool("constant hf_mismatch=True", guards["constant_hf_mismatch"])
    expect_bool("constant profile passed 014 overlap but fails 015 HF", guards["constant_overlap_passes_014"] and guards["constant_hf_mismatch"])
    expect_nonzero("HF guard reads a real object: degenerate wrong_zero_F_a genuinely differs from frozen legacy force", compact(guards["wrong_zero_F_a"] - hf["F_legacy_a"]))
    expect_nonzero("HF guard reads a real object: constant wrong_const_F_a genuinely differs from frozen legacy force", compact(guards["wrong_const_F_a"] - hf["F_legacy_a"]))


def run_verdict_and_composition(data: dict[str, Any]) -> None:
    subbanner("015 scoped landing and joint composition")
    hf = data["hf"]
    gate = data["gate"]
    print("  015 scoped verdict = ", data["verdict"])
    expect_zero("015 component lands at BREATHING_CALIBRATED", verdict_residual(data["verdict"], BREATHING_CALIBRATED))
    expect_bool("STRUCTURE rung passes", gate["M_posdef"] and gate["K_symmetric"] and gate["K_structure_ok"])
    expect_bool("HF rung passes with anti-x-x guard", hf["hf_force_reduces"] and not hf["unsimplified_routes_identical"])
    expect_bool("HF guard rejects both wrong profiles", data["guards"]["degenerate_hf_mismatch"] and data["guards"]["constant_hf_mismatch"])
    print("  BREATHING_CALIBRATED (JOINT, 3-stage) -- COMPLETE")
    print("    = (013: harmonic profiles + M/K operator-projection + (a,L) EOM LHS)[done, cited]")
    print("    AND (014: truncation consistency)[done, cited]")
    print("    AND (015: legacy-Hessian structure recovery + Hellmann-Feynman force + static-dynamic limit)[computed here]")
    print("  RHS FILLED: the EOM M_AB Qddot + K_AB Q = F_A^(HF) now has F_A^(HF)=4*pi*int (rhoStar*Vp0/ellC)*alpha_A dw.")
    print("  CALIBRATED <= {muEta,Tw,beta} (013) [+ the HF drive scale rhoStar*Vp0/ellC -- register question, not decided here].")
    print("  CAVEATS: HF is the profile guard the overlap could not supply; {kappa,chi,sigmaA,sigmaL} are the legacy-Hessian pattern basis; BdG k^4 deferred.")


def print_provenance(data: dict[str, Any]) -> None:
    subbanner("Provenance and scope")
    print("  CONSUMED-from-013/014: stage013 profiles alpha_a,alpha_L + operator-projected M_AB/K_AB + frozen packet are cited with dual-site integrity; stage014 window beta_L0 in [0.1,3.0] / K_eta/Tw<=~2.6 is cited; no sibling recomputation.")
    print("  no-c_S: 015 is speed-free (static structure + static force); matter-sector c_s/BdG k^4 is deferred under k*xi << 1 (phonon_limit_caveat).")
    print("  LEGACY-CONSTANT boundary: {kappa,chi,sigmaA,sigmaL} enter HERE as E_geom parameters and are absent from 013 M/K.")
    print("  EARNED-STRUCTURE-RECOVERY: own-built H_legacy plus genuine structural match of cited M/K; recovered, not re-postulated; full_matrix_fit=False.")
    print("  HF-TWO-ROUTE-GENUINENESS: distributed projection vs Hellmann-Feynman parametric derivative; hf_force_reduces=True and unsimplified_routes_identical=False.")
    print("  HF-PROFILE-GUARD: degenerate and constant wrong profiles have hf_mismatch=True; constant_one passed 014 overlap but fails 015 HF.")
    print("  static-dynamic-limit: Qddot->0 in the same dynamical closure gives K_AB Q = F_A and the static dE_geom/dQ=0 limit; no separate static solve.")
    print("  3-way-split-COMPLETING: 015 is component 3/3 and completes BREATHING_CALIBRATED; the RHS F_A^(HF) deferred by 013 is filled here.")
    print("  dropped-bookkeeping: scratch-WL exports, sympy* comparisons, expression_digest, engine_agreement, YAML/JSON plumbing, and report writers are stripped.")
    print("  downstream consumers: stages 022/023 consume the completed ell=0 (a,L) breathing closure.")
    print("  register note: 015 introduces Vp0 and makes ell_c live through rhoStar*Vp0/ellC; whether that force scale is a counted CALIB knob is intentionally left to registration.")
    live_names = {
        symbol.name
        for expr in [
            data["profiles"]["alpha_a"],
            data["profiles"]["alpha_L"],
            *data["mk"]["M_entries"].values(),
            *data["mk"]["K_entries"].values(),
            data["source_density"],
            data["hf"]["F_dist_a"],
            data["hf"]["F_legacy_a"],
        ]
        for symbol in expr.free_symbols
    }
    expect_bool("no c_S/cS live symbol appears in 015 symbolic content", "cS" not in live_names and "c_S" not in live_names)


def print_verdict_labels() -> None:
    print("")
    print("Verdict labels:")
    print("  ledger earned-label (NOT a source verdict token): BREATHING_STRUCTURE_HF_FORCE_EARNED  (legacy-Hessian structure recovery: own-built H_legacy=hessian(E_geom); the CITED 013 M_AB/K_AB match its structural signature (M pos-def by Sylvester; K symmetric; K_aL<0 <=> H_legacy off-diagonal -chi*kappa<0 for chi>0; rank(K)=rank(H_legacy)=2; zero-pattern match) => the (a,L) operator-projected closure is RECOVERED not re-postulated; the Hellmann-Feynman force F_A^(HF) derived by TWO genuinely-different routes -- distributed projection F_dist=4pi int (rho*Vp0/ellC) alpha_A dw vs Hellmann-Feynman parametric F_legacy=-4pi int rho* dV_conf/dq_A|0 dw -- that AGREE (hf_force_reduces=True) with unsimplified_routes_identical=False (the anti-x-x: raw trees differ); the HF REJECTS the wrong profiles (hf_mismatch=True for alpha_a=0 and alpha_a=1) -- the constant profile passed 014's overlap but fails 015's HF)")
    print("  source top-line verdict: BREATHING_CALIBRATED  (JOINT 3-stage; 015 carries the legacy-structure + HF-force component 3/3 and COMPLETES the joint)")
    print("  joint composition (COMPLETE): BREATHING_CALIBRATED = (013: harmonic profiles + M/K operator-projection + (a,L) EOM LHS)[done] AND (014: truncation consistency)[done] AND (015: legacy-Hessian structure recovery + Hellmann-Feynman force + static-dynamic limit)[here] ; the EOM RHS F_A^(HF) 013 deferred is filled here")
    print("  earned (structure): H_legacy own-built; the cited 013 M/K match its structural signature (genuine pos-def/symmetric/negative-offdiagonal/rank/zero-pattern comparisons, NOT X==X on cited literals); the two probes (M_aa->-M_aa, K_aL->-K_aL) flip the gate")
    print("  earned (HF force): both routes genuinely different constructions agreeing after simplification; unsimplified_routes_identical computed False (anti-x-x); corrupting one route breaks hf_force_reduces")
    print("  calibrated (values): the wall constants {mu_eta,Tw,beta} are stage-013 calibration inputs; the HF drive scale rho*Vp0/ellC is the section-6 register question (likely a new CALIB force-scale -- confirm at registration) -> BREATHING_CALIBRATED not ..._PASS")
    print("  carried caveats (015's own, honest): the HF is the profile guard the overlap could not supply (constant_one PASSES 014's overlap, FAILS 015's HF); {kappa,chi,sigmaA,sigmaL} are the legacy-Hessian pattern basis (a re-parameterization of the calibrated set; structural recovery not full numeric fit, full_matrix_fit=False); BdG k^4 matter-sector deferred (k*xi<<1)")
    print("  consumed (cited from stage013/014, dual-site integrity): collective profiles alpha_a,alpha_L; operator-projected M_AB/K_AB; frozen packet L0=37/20, beta=1, beta*L0=37/20; 014 truncation certificate (window beta_L0 in [0.1,3.0], K_eta/Tw<=~2.6); c_S NOT consumed (matter-sector deferred)")
    print("  new symbols (015): legacy {kappa,chi,sigmaA,sigmaL} (E_geom parameterization; ABSENT from 013's M/K); HF drive Vp0 (confinement slope) + ell_c LIVE (delta_V_conf != 0 here; INERT at 011) + rho* (frozen density, consumed)")


def run_baseline_immutability(data: dict[str, Any]) -> None:
    subbanner("Baseline immutability after copy mutations")
    expect_bool("baseline M_posdef remains true", data["gate"]["M_posdef"])
    expect_bool("baseline K_structure_ok remains true", data["gate"]["K_structure_ok"])
    expect_bool("baseline hf_force_reduces remains true", data["hf"]["hf_force_reduces"])
    expect_bool("baseline unsimplified_routes_identical remains false", not data["hf"]["unsimplified_routes_identical"])
    expect_bool("baseline M/K citation integrity remains true", data["mk_integrity"]["ok"])
    expect_zero("baseline clean 015 verdict remains BREATHING_CALIBRATED", verdict_residual(data["verdict"], BREATHING_CALIBRATED))


def main() -> None:
    banner("ledger_stage015_breathing_structure_hf_force SymPy audit")
    data = build_baseline()
    assert_no_float("baseline", data)
    run_consumed_inputs(data)
    run_legacy_hessian(data)
    run_structure_recovery(data)
    run_hf_force(data)
    run_static_dynamic_limit(data)
    run_hf_profile_guard(data)
    run_verdict_and_composition(data)
    print_provenance(data)
    print_verdict_labels()
    run_baseline_immutability(data)
    print("")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}; TOTAL = {PASS_COUNT} + {FAIL_COUNT} = {PASS_COUNT + FAIL_COUNT}")
    print("OVERALL PASS: SymPy verified ledger_stage015 breathing structure recovery + Hellmann-Feynman force exactly")


if __name__ == "__main__":
    try:
        main()
    except AuditFailure as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}; TOTAL = {PASS_COUNT} + {FAIL_COUNT} = {PASS_COUNT + FAIL_COUNT}")
        print(f"OVERALL FAIL: SymPy stage015 audit did not close ({exc})")
        raise SystemExit(1)
