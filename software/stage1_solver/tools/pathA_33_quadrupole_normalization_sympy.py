#!/usr/bin/env python3
"""PathA-33 Gate 4 quadrupole normalization, SymPy engine.

Run order:

  timeout 600 python software/stage1_solver/tools/pathA_33_quadrupole_normalization_sympy.py
  timeout 600 math -script software/stage1_solver/tools/pathA_33_quadrupole_normalization.wl
  timeout 600 python software/stage1_solver/tools/pathA_33_quadrupole_normalization_sympy.py

The first Python run writes the SymPy scratch lane.  The second Python run,
after Mathematica has written its independent scratch lane, compares the
engines and emits the final YAML, narrative report, and pde_ledger feed note.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import sympy as sp
import yaml


class NoAliasDumper(yaml.SafeDumper):
    def ignore_aliases(self, data: Any) -> bool:
        return True


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
REPO_ROOT = SCRIPT_PATH.parents[3]
REPORTS = STAGE1_ROOT / "reports"
SCRATCH = STAGE1_ROOT / "_scratch"
NOTES = REPO_ROOT / "research" / "pde_ledger" / "notes" / "stages"

SYM_YAML = SCRATCH / "pathA_33_sympy_results.yaml"
MMA_YAML = SCRATCH / "pathA_33_mathematica_results.yaml"
RESULTS_YAML = REPORTS / "pathA_33_results.yaml"
REPORT_MD = REPORTS / "pathA_33_quadrupole_normalization.md"
FEED_NOTE = NOTES / "moving_throat_pde_pathA_33_quadrupole_normalization_result.md"

SYMBOLIC_TOL = 1.0e-10
NUMERIC_TOL = 1.0e-8

z, omega = sp.symbols("z omega")
a, c_s, c, G = sp.symbols("a c_s c G", positive=True)
D0, D2, D4, N0, N2, N4 = sp.symbols("D0 D2 D4 N0 N2 N4", nonzero=True)
mu_hat0, mtilde0, chi_Q = sp.symbols("mu_hat0 mtilde0 chi_Q", nonzero=True)
lambda_G = sp.symbols("lambda_G", nonzero=True)

OmegaU, OmegaW, gU, gW, R, Delta, S = sp.symbols(
    "Omega_U Omega_W g_U g_W R Delta S", nonzero=True
)

Ldim, Mdim, Tdim = sp.symbols("L M T", positive=True)
Dim = tuple[sp.Rational, sp.Rational, sp.Rational]
ZERO_DIM: Dim = (sp.Rational(0), sp.Rational(0), sp.Rational(0))


def compact(expr: Any) -> Any:
    if isinstance(expr, sp.Basic):
        return sp.factor(sp.cancel(sp.simplify(expr)))
    return expr


def hstr(expr: Any) -> str | bool | int | float | None:
    if expr is None or isinstance(expr, (bool, int, float, str)):
        return expr
    return sp.sstr(compact(expr))


def yaml_write(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.dump(payload, Dumper=NoAliasDumper, sort_keys=False, allow_unicode=False),
        encoding="utf-8",
    )


def yaml_read(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise AssertionError(f"YAML did not load to a mapping: {path}")
    return data


def bool_from_residual(expr: sp.Expr) -> bool:
    return compact(expr) == 0


def numeric(expr: sp.Expr, subs: dict[sp.Symbol, float | int]) -> float:
    return float(sp.N(compact(expr).subs(subs), 40))


def series_no_o(expr: sp.Expr, var: sp.Symbol, order: int) -> sp.Expr:
    return sp.expand(sp.series(expr, var, 0, order).removeO())


def spherical_j2() -> sp.Expr:
    return (3 / z**3 - 1 / z) * sp.sin(z) - 3 * sp.cos(z) / z**2


def spherical_y2() -> sp.Expr:
    return (1 / z - 3 / z**3) * sp.cos(z) - 3 * sp.sin(z) / z**2


def dtn_branch(kind: str) -> dict[str, Any]:
    j2 = spherical_j2()
    y2 = spherical_y2()
    if kind == "outgoing_hankel1":
        h = j2 + sp.I * y2
        source = "dtn_hankel1"
    elif kind == "incoming_hankel2":
        h = j2 - sp.I * y2
        source = "dtn_hankel2"
    elif kind == "standing_j2":
        h = j2
        source = "regular_standing_j2"
    else:
        raise ValueError(kind)

    lam = compact(z * sp.diff(h, z) / h)
    yhat = compact(-3 / lam)
    h_series = series_no_o(h, z, 8)
    lam_series = series_no_o(lam, z, 7)
    yhat_series = series_no_o(yhat, z, 7)
    u2_z = compact(yhat_series.coeff(z, 2))
    u4_z = compact(yhat_series.coeff(z, 4))
    v5_z = compact(yhat_series.coeff(z, 5) / sp.I)
    static = compact(yhat_series.coeff(z, 0))
    return {
        "kind": kind,
        "source": source,
        "h": h,
        "lambda": lam,
        "yhat": yhat,
        "h_series": h_series,
        "lambda_series": lam_series,
        "yhat_series": yhat_series,
        "static": static,
        "u2_z": u2_z,
        "u4_z": u4_z,
        "v5_z": v5_z,
        "u2": compact(u2_z * a**2 / c_s**2),
        "u4": compact(u4_z * a**4 / c_s**4),
        "v5": compact(v5_z * a**5 / c_s**5),
    }


def build_fingerprint() -> dict[str, Any]:
    out = dtn_branch("outgoing_hankel1")
    incoming = dtn_branch("incoming_hankel2")
    standing = dtn_branch("standing_j2")
    expected = {
        "u2_z": sp.Rational(1, 9),
        "u4_z": sp.Rational(4, 81),
        "v5_z": sp.Rational(1, 27),
        "u2": a**2 / (9 * c_s**2),
        "u4": 4 * a**4 / (81 * c_s**4),
        "v5": a**5 / (27 * c_s**5),
    }
    matches = {
        name: bool_from_residual(out[name] - target)
        for name, target in expected.items()
    }
    chi_derived = compact(out["v5"] / (a**5 / (27 * c_s**5)))
    chi_incoming = compact(incoming["v5"] / (a**5 / (27 * c_s**5)))
    sample_subs = {a: 3.0, c_s: 2.0}
    return {
        "outgoing": out,
        "incoming": incoming,
        "standing": standing,
        "expected": expected,
        "matches": matches,
        "ok": all(matches.values()) and bool_from_residual(chi_derived - 1),
        "chi_Q": chi_derived,
        "chi_Q_incoming": chi_incoming,
        "samples": {
            "u2": numeric(out["u2"], sample_subs),
            "u4": numeric(out["u4"], sample_subs),
            "v5": numeric(out["v5"], sample_subs),
            "chi_Q": numeric(chi_derived, sample_subs),
            "incoming_chi_Q": numeric(chi_incoming, sample_subs),
        },
    }


def build_port_moments() -> dict[str, Any]:
    P_port = OmegaU**2 * gW + R * gU
    n0 = compact(P_port**2 / Delta**2)
    n2 = compact(2 * P_port * (P_port * S - Delta * gW) / Delta**3)
    n4 = compact(
        (
            Delta**2 * gW**2
            - 2 * Delta * P_port**2
            - 4 * Delta * P_port * S * gW
            + 3 * P_port**2 * S**2
        )
        / Delta**4
    )
    return {
        "P_port": P_port,
        "N_A0_r": n0,
        "N_A2_r": n2,
        "N_A4_r": n4,
        "isotropic_branch": "N_20,n=N_21,n=N_22,n=N_n and D_20,n=D_21,n=D_22,n=D_n carried symbolically",
    }


def build_prefactor() -> dict[str, Any]:
    Nomega = N0 + N2 * omega**2 + N4 * omega**4
    Dcons = D0 + D2 * omega**2 + D4 * omega**4
    correct_obj = D0 * Nomega / Dcons**2
    plain_obj = Nomega / Dcons
    correct_series = series_no_o(correct_obj, omega, 6)
    plain_series = series_no_o(plain_obj, omega, 6)
    coeffs = {
        "P0": compact(correct_series.coeff(omega, 0)),
        "P2": compact(correct_series.coeff(omega, 2)),
        "P4": compact(correct_series.coeff(omega, 4)),
    }
    expected = {
        "P0": N0 / D0,
        "P2": (D0 * N2 - 2 * D2 * N0) / D0**2,
        "P4": (
            D0**2 * N4
            - 2 * D0 * (D2 * N2 + D4 * N0)
            + 3 * D2**2 * N0
        )
        / D0**3,
    }
    residuals = {name: compact(coeffs[name] - expected[name]) for name in coeffs}
    matches = {name: bool_from_residual(residual) for name, residual in residuals.items()}
    plain = {
        "P0": compact(plain_series.coeff(omega, 0)),
        "P2": compact(plain_series.coeff(omega, 2)),
        "P4": compact(plain_series.coeff(omega, 4)),
    }
    sample_subs = {
        D0: 19.0,
        D2: 23.0,
        D4: 29.0,
        N0: 11.0,
        N2: 13.0,
        N4: 17.0,
    }
    return {
        "correct_object": correct_obj,
        "plain_object": plain_obj,
        "correct_series": correct_series,
        "plain_series": plain_series,
        "coefficients": coeffs,
        "expected": expected,
        "residuals": residuals,
        "matches": matches,
        "ok": all(matches.values()),
        "plain_N_over_D": plain,
        "self_check": {
            "correct_P2_D2N0_term": "-2*D2*N0",
            "plain_P2_D2N0_term": "-D2*N0",
            "plain_equals_correct_P2": bool_from_residual(plain["P2"] - expected["P2"]),
            "difference_plain_minus_correct_P2": compact(plain["P2"] - expected["P2"]),
        },
        "sample_values": {
            "P0": numeric(coeffs["P0"], sample_subs),
            "P2": numeric(coeffs["P2"], sample_subs),
            "P4": numeric(coeffs["P4"], sample_subs),
            "plain_P2": numeric(plain["P2"], sample_subs),
            "plain_minus_correct_P2": numeric(plain["P2"] - expected["P2"], sample_subs),
        },
    }


def a_power(expr: sp.Expr) -> sp.Rational | None:
    powers = sp.factor(expr).as_powers_dict()
    power = powers.get(a, sp.Rational(0))
    return sp.Rational(power) if power.is_number else None


def build_scaling() -> dict[str, Any]:
    target_rhs = 54 * G * c_s**5 / (5 * a**5 * c**5)
    mutated_rhs = 54 * G * c_s**5 / (5 * a**4 * c**5)
    power = a_power(target_rhs)
    mutated_power = a_power(mutated_rhs)
    return {
        "target_rhs": target_rhs,
        "P0_target_scaling": int(power) if power is not None else None,
        "target_scaling_ok": power == -5,
        "actual_branch_scaling_status": "DEFERRED_GATE_6_SYMBOLIC_Dn_Nn_NO_A_SCALINGS_SUPPLIED",
        "wrong_scaling_probe_power": int(mutated_power) if mutated_power is not None else None,
        "wrong_scaling_probe_ok": mutated_power != -5,
    }


class DimError(ValueError):
    pass


def dim_add(left: Dim, right: Dim) -> Dim:
    return tuple(sp.Rational(a0 + b0) for a0, b0 in zip(left, right))  # type: ignore[return-value]


def dim_sub(left: Dim, right: Dim) -> Dim:
    return tuple(sp.Rational(a0 - b0) for a0, b0 in zip(left, right))  # type: ignore[return-value]


def dim_scale(dim: Dim, scale: sp.Rational) -> Dim:
    return tuple(sp.Rational(scale * d0) for d0 in dim)  # type: ignore[return-value]


def dim_of(expr: sp.Expr, symbol_dims: dict[sp.Symbol, Dim]) -> Dim:
    expr = sp.sympify(expr)
    if expr == 0 or expr.is_number:
        return ZERO_DIM
    if expr.is_Symbol:
        if expr not in symbol_dims:
            raise DimError(f"missing dimension for symbol {expr}")
        return symbol_dims[expr]
    if expr.is_Mul:
        out = ZERO_DIM
        for arg in expr.args:
            out = dim_add(out, dim_of(arg, symbol_dims))
        return out
    if expr.is_Pow:
        base, exponent = expr.args
        if not exponent.is_number:
            raise DimError(f"non-numeric dimension exponent in {expr}")
        return dim_scale(dim_of(base, symbol_dims), sp.Rational(exponent))
    if expr.is_Add:
        dims = [dim_of(arg, symbol_dims) for arg in expr.args if arg != 0]
        if not dims:
            return ZERO_DIM
        first = dims[0]
        mismatched = [dim for dim in dims if dim != first]
        if mismatched:
            raise DimError(f"dimension mismatch in sum {expr}: {dims}")
        return first
    raise DimError(f"unsupported expression in dimension checker: {expr}")


def dim_to_monomial(dim: Dim) -> sp.Expr:
    return compact(Ldim ** dim[0] * Mdim ** dim[1] * Tdim ** dim[2])


def exp_text(exp: sp.Rational) -> str:
    exp = sp.Rational(exp)
    if exp.q == 1:
        return str(exp.p)
    return f"{exp.p}/{exp.q}"


def dim_to_text(dim: Dim) -> str:
    # Internal order is (L,M,T); display order keeps the common L,T,M reading.
    parts: list[str] = []
    for name, exp in (("L", dim[0]), ("T", dim[2]), ("M", dim[1])):
        exp = sp.Rational(exp)
        if exp == 0:
            continue
        if exp == 1:
            parts.append(name)
        else:
            parts.append(f"{name}^{exp_text(exp)}")
    return " ".join(parts) if parts else "1"


def dim_record(name: str, expr: sp.Expr, symbol_dims: dict[sp.Symbol, Dim]) -> dict[str, Any]:
    dim = dim_of(expr, symbol_dims)
    return {
        "quantity": name,
        "expression": hstr(expr),
        "dimension": dim_to_text(dim),
        "dimension_monomial": hstr(dim_to_monomial(dim)),
        "dimension_vector_LMT": [hstr(v) for v in dim],
    }


SOURCED_N0_DIM: Dim = (-1, 1, 0)
SOURCED_D0_DIM: Dim = (-1, 1, -2)


def dim_vector_text(dim: Dim) -> list[str]:
    return [hstr(v) for v in dim]


def build_dimensions(fingerprint: dict[str, Any]) -> dict[str, Any]:
    raw_symbol_dims: dict[sp.Symbol, Dim] = {
        a: (1, 0, 0),
        c_s: (1, 0, -1),
        c: (1, 0, -1),
        G: (3, -1, -2),
        omega: (0, 0, -1),
        D0: SOURCED_D0_DIM,
        N0: SOURCED_N0_DIM,
        chi_Q: ZERO_DIM,
        mtilde0: ZERO_DIM,
    }
    P0_raw = N0 / D0
    frequency_normalization = (c_s / a) ** 2
    P0_physical = frequency_normalization * P0_raw
    yhat_physical = (
        1
        + fingerprint["outgoing"]["u2"] * omega**2
        + fingerprint["outgoing"]["u4"] * omega**4
        + sp.I * fingerprint["outgoing"]["v5"] * omega**5
    )
    Gamma5 = chi_Q * P0_physical * a**5 / (27 * c_s**5)
    target_rhs = 54 * G * c_s**5 / (5 * a**5 * c**5)

    raw_p0_dim = dim_of(P0_raw, raw_symbol_dims)
    frequency_norm_dim = dim_of(frequency_normalization, raw_symbol_dims)
    p0_dim = dim_of(P0_physical, raw_symbol_dims)
    dimensional_ok = p0_dim == ZERO_DIM

    drop_norm_dim = dim_of(P0_raw, raw_symbol_dims)
    drop_norm_ok = drop_norm_dim == ZERO_DIM
    drop_norm_verdict = "NO_FAIL" if drop_norm_ok else "FAIL_DIMENSIONAL"

    corrupt_n0_dims = dict(raw_symbol_dims)
    corrupt_n0_dims[N0] = ZERO_DIM
    corrupt_n0_raw_dim = dim_of(P0_raw, corrupt_n0_dims)
    corrupt_n0_p0_dim = dim_of(P0_physical, corrupt_n0_dims)
    corrupt_n0_ok = corrupt_n0_p0_dim == ZERO_DIM
    corrupt_n0_verdict = "NO_FAIL" if corrupt_n0_ok else "FAIL_DIMENSIONAL"

    rhs_dim = dim_of(target_rhs, raw_symbol_dims)
    mu_dim = dim_scale(dim_sub(rhs_dim, p0_dim), sp.Rational(1, 2))
    symbol_dims = dict(raw_symbol_dims)
    symbol_dims[mu_hat0] = mu_dim
    lhs = (mu_hat0 * mtilde0) ** 2 * P0_physical
    lhs_raw_mutation = (mu_hat0 * mtilde0) ** 2 * P0_raw
    required_p0_dim = dim_sub(rhs_dim, dim_scale(mu_dim, sp.Rational(2)))

    table = [
        dim_record("N0", N0, symbol_dims),
        dim_record("D0", D0, symbol_dims),
        dim_record("P0_raw=N0/D0", P0_raw, symbol_dims),
        dim_record("frequency_normalization=(c_s/a)^2", frequency_normalization, symbol_dims),
        dim_record("P0=(c_s/a)^2*N0/D0", P0_physical, symbol_dims),
        dim_record("Yhat_out_physical_expansion", yhat_physical, symbol_dims),
        dim_record("Gamma5=chi_Q*P0*a^5/(27*c_s^5)", Gamma5, symbol_dims),
        dim_record("target_lhs=(mu_hat0*mtilde0)^2*P0", lhs, symbol_dims),
        dim_record("target_rhs=54*G*c_s^5/(5*a^5*c^5)", target_rhs, symbol_dims),
    ]
    lhs_dim = dim_of(lhs, symbol_dims)
    lhs_raw_dim = dim_of(lhs_raw_mutation, symbol_dims)
    homogeneity_pass = lhs_dim == rhs_dim and p0_dim == required_p0_dim
    return {
        "dimension_order": ["L", "M", "T"],
        "dimensional_gate": "mu_hat0-free P0 physical dimension check",
        "dimensional_gate_expression": "P0_phys=(c_s/a)^2*(N0/D0)",
        "dimensional_gate_depends_on_mu_hat0": False,
        "symbol_dimensions": {
            "a": "L",
            "c_s": "L T^-1",
            "c": "L T^-1",
            "G": "L^3 T^-2 M^-1",
            "N0": "L^-1 M",
            "D0": "L^-1 T^-2 M",
            "mtilde0": "1",
            "mu_hat0": dim_to_text(mu_dim),
        },
        "P0_raw_dimension": dim_to_text(raw_p0_dim),
        "frequency_normalization_dimension": dim_to_text(frequency_norm_dim),
        "P0_dimension": dim_to_text(p0_dim),
        "P0_physical_dimension": dim_to_text(p0_dim),
        "P0_physical_dimension_vector_LMT": dim_vector_text(p0_dim),
        "dimensional_ok": dimensional_ok,
        "dimensional_status": "DIMENSIONAL_OK" if dimensional_ok else "FAIL_DIMENSIONAL",
        "mu_hat0_dimension": dim_to_text(mu_dim),
        "Gamma5_dimension": dim_to_text(dim_of(Gamma5, symbol_dims)),
        "lhs_dimension": dim_to_text(lhs_dim),
        "rhs_dimension": dim_to_text(rhs_dim),
        "table": table,
        "mu_hat0_homogeneity_diagnostic": {
            "label": "non-able-to-fail (mu_hat0 free carrier)",
            "participates_in_verdict": False,
            "mu_dim_solved_from_rhs_minus_p0_over_2": True,
            "homogeneity_required_P0_dimension_with_pinned_mu": dim_to_text(required_p0_dim),
            "mu_hat0_dimension": dim_to_text(mu_dim),
            "lhs_dimension": dim_to_text(lhs_dim),
            "rhs_dimension": dim_to_text(rhs_dim),
            "lhs_raw_without_frequency_normalization_dimension": dim_to_text(lhs_raw_dim),
            "homogeneity_pass": homogeneity_pass,
        },
        "mutation_drop_frequency_normalization": {
            "mutated_P0_expression": hstr(P0_raw),
            "mutated_P0_dimension": dim_to_text(drop_norm_dim),
            "mutated_P0_dimension_vector_LMT": dim_vector_text(drop_norm_dim),
            "baseline_P0_physical_dimension": dim_to_text(p0_dim),
            "dimensional_ok": drop_norm_ok,
            "verdict": drop_norm_verdict,
            "mutation_fires": drop_norm_verdict == "FAIL_DIMENSIONAL",
        },
        "mutation_corrupt_N0_dimension": {
            "sourced_N0_dimension": dim_to_text(SOURCED_N0_DIM),
            "corrupted_N0_dimension": dim_to_text(ZERO_DIM),
            "mutated_P0_raw_dimension": dim_to_text(corrupt_n0_raw_dim),
            "mutated_P0_physical_dimension": dim_to_text(corrupt_n0_p0_dim),
            "mutated_P0_physical_dimension_vector_LMT": dim_vector_text(corrupt_n0_p0_dim),
            "baseline_P0_physical_dimension": dim_to_text(p0_dim),
            "dimensional_ok": corrupt_n0_ok,
            "verdict": corrupt_n0_verdict,
            "mutation_fires": corrupt_n0_verdict == "FAIL_DIMENSIONAL",
        },
    }


def build_equivalence(fingerprint: dict[str, Any]) -> dict[str, Any]:
    chi = fingerprint["chi_Q"]
    P0_physical = (c_s / a) ** 2 * (N0 / D0)
    Gamma5 = chi_Q * P0_physical * a**5 / (27 * c_s**5)
    target_rhs = 54 * G * c_s**5 / (5 * a**5 * c**5)
    gamma_target = 2 * G / (5 * c**5)
    forward = compact((target_rhs * chi * a**5 / (27 * c_s**5)) - gamma_target)
    reverse = compact((gamma_target * 27 * c_s**5 / (chi * a**5)) - target_rhs)
    wrong_gamma = 3 * G / (5 * c**5)
    wrong_reverse = compact((wrong_gamma * 27 * c_s**5 / (chi * a**5)) - target_rhs)
    ok = bool_from_residual(forward) and bool_from_residual(reverse)
    return {
        "definitions": {
            "P0_physical": hstr(P0_physical),
            "Gamma5": hstr(Gamma5),
            "gamma_quad_eff": "mhat0^2*Gamma5 with mhat0=mu_hat0*mtilde0",
            "normalization_target": hstr(target_rhs),
        },
        "chi_Q_derived_from_DtN": hstr(chi),
        "Gamma5_with_derived_chi": hstr(Gamma5.subs(chi_Q, chi)),
        "forward_residual_to_2G_over_5c5": hstr(forward),
        "reverse_residual_to_54Gcs5_over_5a5c5": hstr(reverse),
        "ok": ok,
        "wrong_gamma_probe_residual": hstr(wrong_reverse),
        "wrong_gamma_probe_fires": not bool_from_residual(wrong_reverse),
    }


DERIVED_TAGS = {
    "dtn_hankel_expansion",
    "dtn_radiative_slot",
    "prefactor_series_algebra",
    "target_rhs_algebra",
    "gamma_bridge_algebra",
    "dimension_expression_eval",
    "emergent_outgoing_passivity",
}
EXTERNAL_TAGS = {"external_gr_constant", "external_pn_bridge", "einstein_bridge_identity"}
DEFERRED_TAGS = {"gate6_branch_solve", "deferred_nonlinear_pde"}
CONVENTION_TAGS = {"normalization_convention", "unit_choice", "static_slot_convention"}


PROVENANCE_SOURCES: dict[str, list[str]] = {
    "fingerprint_u2": ["dtn_hankel_expansion"],
    "fingerprint_u4": ["dtn_hankel_expansion"],
    "fingerprint_v5": ["dtn_hankel_expansion", "dtn_radiative_slot"],
    "fingerprint_27": ["dtn_radiative_slot"],
    "prefactor_P0_P2_P4": ["prefactor_series_algebra"],
    "P0_target_scaling_minus5": ["target_rhs_algebra"],
    "chi_Q": ["dtn_radiative_slot"],
    "Gamma5_equivalence_chain": ["gamma_bridge_algebra"],
    "dimensional_p0_physical_gate": ["dimension_expression_eval"],
    "emergent_passivity": ["emergent_outgoing_passivity"],
    "G": ["external_gr_constant"],
    "PN_2_over_5": ["external_pn_bridge"],
    "Einstein_2G_over_5c5_identity": ["einstein_bridge_identity"],
    "assembled_54_over_5_magnitude": ["external_pn_bridge", "dtn_radiative_slot"],
    "D_n_N_n_numeric_values": ["gate6_branch_solve"],
    "port_scalars": ["gate6_branch_solve"],
    "actual_branch_a_scaling": ["gate6_branch_solve"],
    "unit_choices": ["unit_choice"],
    "static_slot_minus3": ["static_slot_convention"],
}


def classify_provenance(tags: list[str]) -> str:
    tagset = set(tags)
    if tagset & DEFERRED_TAGS:
        return "deferred_branch_data"
    if tagset & EXTERNAL_TAGS:
        return "external_bridge_input"
    if tagset & DERIVED_TAGS:
        return "derived_in_gate"
    if tagset & CONVENTION_TAGS:
        return "convention"
    return "unclassified"


def group_partition(items: dict[str, dict[str, Any]], class_field: str) -> dict[str, list[str]]:
    groups = {
        "derived_in_gate": [],
        "external_bridge_input": [],
        "deferred_branch_data": [],
        "convention": [],
        "unclassified": [],
    }
    for name, item in items.items():
        groups.setdefault(item[class_field], []).append(name)
    return groups


def build_partition(overrides: dict[str, str] | None = None) -> dict[str, Any]:
    overrides = overrides or {}
    items = {}
    for name, tags in PROVENANCE_SOURCES.items():
        computed = classify_provenance(tags)
        emitted = overrides.get(name, computed)
        items[name] = {
            "provenance_tags": tags,
            "computed_class": computed,
            "emitted_class": emitted,
            "class_matches_computed": emitted == computed,
        }
    ok = all(item["class_matches_computed"] for item in items.values())
    g_invariance = build_g_invariance_diagnostic(items)
    return {
        "items": items,
        "groups": group_partition(items, "computed_class"),
        "emitted_groups": group_partition(items, "emitted_class"),
        "ok": ok,
        "classification_rule": "computed from provenance tags; external tags dominate mixed quantities, so assembled 54/5 remains calibrated",
        "g_to_lambda_g_diagnostic": g_invariance,
        "decomposition_54_over_5": {
            "identity": "54/5 = 2*27/5",
            "earned_factor": {"factor": 27, "class": items["fingerprint_27"]["computed_class"]},
            "calibrated_factor": {"factor": "2/5", "class": items["PN_2_over_5"]["computed_class"]},
            "assembled_magnitude_class": items["assembled_54_over_5_magnitude"]["computed_class"],
        },
    }


def is_g_invariant(expr: sp.Expr) -> bool:
    return bool_from_residual(expr.subs(G, lambda_G * G) - expr)


def build_g_invariance_diagnostic(items: dict[str, dict[str, Any]]) -> dict[str, Any]:
    diagnostics = {
        "G": {"expression": G, "class": items["G"]["computed_class"]},
        "target_2G_over_5c5": {
            "expression": 2 * G / (5 * c**5),
            "class": items["PN_2_over_5"]["computed_class"],
        },
        "pure_54_over_5": {
            "expression": sp.Rational(54, 5),
            "class": items["assembled_54_over_5_magnitude"]["computed_class"],
        },
        "fingerprint_27": {
            "expression": sp.Integer(27),
            "class": items["fingerprint_27"]["computed_class"],
        },
    }
    out = {}
    for name, data in diagnostics.items():
        out[name] = {
            "expression": hstr(data["expression"]),
            "g_invariant": is_g_invariant(data["expression"]),
            "provenance_class": data["class"],
        }
    return {
        "diagnostics": out,
        "classifier": "provenance_not_g_invariance",
        "invariance_only_trap_catches_54_over_5": out["pure_54_over_5"]["g_invariant"]
        and out["pure_54_over_5"]["provenance_class"] == "external_bridge_input",
    }


def passivity_from_source(branch: dict[str, Any]) -> dict[str, Any]:
    imag_nonzero = not bool_from_residual(branch["v5_z"])
    is_genuine = branch["source"] == "dtn_hankel1"
    return {
        "source": branch["source"],
        "imaginary_v5_nonzero": imag_nonzero,
        "genuine_outgoing": is_genuine and imag_nonzero,
    }


def base_verdict(gates: dict[str, bool], partition: dict[str, Any]) -> str:
    ordered_failures = [
        ("fingerprint_ok", "QUAD_FAIL_FINGERPRINT"),
        ("prefactor_ok", "QUAD_FAIL_PREFACTOR_ALGEBRA"),
        ("scaling_ok", "QUAD_FAIL_SCALING"),
        ("equivalence_ok", "QUAD_FAIL_EQUIVALENCE"),
        ("dimensional_ok", "QUAD_FAIL_DIMENSIONAL"),
        ("outgoing_ok", "QUAD_FAIL_NOT_OUTGOING"),
        ("provenance_ok", "QUAD_FAIL_PROVENANCE_PARTITION"),
        ("able_to_fail_ok", "QUAD_FAIL_NOT_ABLE_TO_FAIL"),
    ]
    for gate, verdict in ordered_failures:
        if not gates.get(gate, False):
            return verdict
    g_class = partition["items"]["G"]["computed_class"]
    mag_class = partition["items"]["assembled_54_over_5_magnitude"]["computed_class"]
    if g_class == "derived_in_gate" and mag_class == "derived_in_gate":
        return "QUAD_PASS"
    return "QUAD_CALIBRATED"


def probe_verdict_label(verdict: str) -> str:
    if verdict.startswith("QUAD_FAIL_"):
        return "FAIL_" + verdict.removeprefix("QUAD_FAIL_")
    return verdict


def build_counterfactuals(
    fingerprint: dict[str, Any],
    prefactor: dict[str, Any],
    scaling: dict[str, Any],
    dimensions: dict[str, Any],
    equivalence: dict[str, Any],
    partition: dict[str, Any],
    baseline_gates: dict[str, bool],
    baseline_verdict: str,
) -> dict[str, Any]:
    out = fingerprint["outgoing"]
    incoming = fingerprint["incoming"]
    standing = fingerprint["standing"]
    incoming_expected = (
        bool_from_residual(incoming["u2_z"] - out["u2_z"])
        and bool_from_residual(incoming["u4_z"] - out["u4_z"])
        and bool_from_residual(incoming["v5_z"] + out["v5_z"])
        and bool_from_residual(fingerprint["chi_Q_incoming"] + 1)
    )
    standing_expected = (
        bool_from_residual(standing["lambda_series"].coeff(z, 0) - 2)
        and bool_from_residual(standing["static"] + sp.Rational(3, 2))
        and bool_from_residual(standing["v5_z"])
    )

    def ablation(
        with_gate_overrides: dict[str, bool],
        *,
        with_partition: dict[str, Any] | None = None,
        expected_fail: str | None = None,
        with_cases: dict[str, str] | None = None,
    ) -> dict[str, Any]:
        mutated_gates = dict(baseline_gates)
        mutated_gates.update(with_gate_overrides)
        with_mutation = probe_verdict_label(base_verdict(mutated_gates, with_partition or partition))
        without_mutation = probe_verdict_label(base_verdict(dict(baseline_gates), partition))
        fail_label = expected_fail or with_mutation
        case_values = list(with_cases.values()) if with_cases else [with_mutation]
        fail_suppressed = (
            without_mutation != fail_label
            and without_mutation not in case_values
            and all(value == fail_label for value in case_values)
        )
        out = {
            "rerun_gate_logic": True,
            "with_mutation": with_mutation,
            "without_mutation": without_mutation,
            "expected_fail": fail_label,
            "fail_suppressed": fail_suppressed,
        }
        if with_cases:
            out["with_mutation_cases"] = with_cases
        return out

    fingerprint_ablation = ablation({"fingerprint_ok": False}, expected_fail="FAIL_FINGERPRINT")
    outgoing_ablation = ablation({"outgoing_ok": False}, expected_fail="FAIL_NOT_OUTGOING")
    scaling_ablation = ablation({"scaling_ok": False}, expected_fail="FAIL_SCALING")
    dimensional_ablation = ablation({"dimensional_ok": False}, expected_fail="FAIL_DIMENSIONAL")
    equivalence_ablation = ablation({"equivalence_ok": False}, expected_fail="FAIL_EQUIVALENCE")
    prefactor_ablation = ablation({"prefactor_ok": False}, expected_fail="FAIL_PREFACTOR_ALGEBRA")

    def partition_ablation(mutated_partition: dict[str, Any]) -> dict[str, Any]:
        return {
            **ablation(
                {"provenance_ok": mutated_partition["ok"]},
                with_partition=mutated_partition,
                expected_fail="FAIL_PROVENANCE_PARTITION",
            )
        }

    probes: dict[str, Any] = {}
    probes["3a_wrong_bc"] = {
        "incoming": {
            "u2_z": hstr(incoming["u2_z"]),
            "u4_z": hstr(incoming["u4_z"]),
            "v5_z": hstr(incoming["v5_z"]),
            "chi_Q": hstr(fingerprint["chi_Q_incoming"]),
            "predicted_change": incoming_expected,
            "verdict": "FAIL_FINGERPRINT" if incoming_expected else "FAIL_NOT_ABLE_TO_FAIL",
        },
        "standing_j2": {
            "lambda_static": hstr(standing["lambda_series"].coeff(z, 0)),
            "Yhat_static": hstr(standing["static"]),
            "v5_z": hstr(standing["v5_z"]),
            "predicted_change": standing_expected,
            "verdict": "FAIL_FINGERPRINT" if standing_expected else "FAIL_NOT_ABLE_TO_FAIL",
        },
        "self_ablation": {
            **fingerprint_ablation,
            "with_mutation_cases": {
                "incoming": "FAIL_FINGERPRINT" if incoming_expected else "FAIL_NOT_ABLE_TO_FAIL",
                "standing_j2": "FAIL_FINGERPRINT" if standing_expected else "FAIL_NOT_ABLE_TO_FAIL",
            },
        },
    }

    phenom_source = {
        "source": "phenomenological_inserted_damping",
        "v5_z": out["v5_z"],
    }
    phenom_genuine = phenom_source["source"] == "dtn_hankel1" and not bool_from_residual(phenom_source["v5_z"])
    probes["3b_imposed_dissipation"] = {
        "mutated_source": phenom_source["source"],
        "inserted_v5_z": hstr(phenom_source["v5_z"]),
        "genuine_outgoing": phenom_genuine,
        "detection_basis": "provenance-tag check; no DtN Hankel source is attached to the inserted damping term",
        "removing_outgoing_bc_removes_imaginary_term": bool_from_residual(standing["v5_z"]),
        "verdict": "FAIL_NOT_OUTGOING" if not phenom_genuine else "NO_FAIL",
        "self_ablation": outgoing_ablation,
    }
    probes["3c_wrong_scaling"] = {
        "mutated_P0_target_scaling": scaling["wrong_scaling_probe_power"],
        "expected_P0_target_scaling": -5,
        "verdict": "FAIL_SCALING" if scaling["wrong_scaling_probe_ok"] else "NO_FAIL",
        "self_ablation": scaling_ablation,
    }
    dim_mutation = dimensions["mutation_drop_frequency_normalization"]
    probes["3d_dimensional_break"] = {
        "mutated_P0_dimension": dim_mutation["mutated_P0_dimension"],
        "baseline_P0_physical_dimension": dimensions["P0_physical_dimension"],
        "mutated_dimensional_ok": dim_mutation["dimensional_ok"],
        "verdict": dim_mutation["verdict"],
        "self_ablation": dimensional_ablation,
    }
    corrupt_dim_mutation = dimensions["mutation_corrupt_N0_dimension"]
    probes["3d_prime_corrupt_port_dimension"] = {
        "sourced_N0_dimension": corrupt_dim_mutation["sourced_N0_dimension"],
        "corrupted_N0_dimension": corrupt_dim_mutation["corrupted_N0_dimension"],
        "mutated_P0_physical_dimension": corrupt_dim_mutation["mutated_P0_physical_dimension"],
        "baseline_P0_physical_dimension": dimensions["P0_physical_dimension"],
        "mutated_dimensional_ok": corrupt_dim_mutation["dimensional_ok"],
        "verdict": corrupt_dim_mutation["verdict"],
        "self_ablation": dimensional_ablation,
    }
    probes["3e_equivalence_break"] = {
        "wrong_gamma_probe_residual": equivalence["wrong_gamma_probe_residual"],
        "verdict": "FAIL_EQUIVALENCE" if equivalence["wrong_gamma_probe_fires"] else "NO_FAIL",
        "self_ablation": equivalence_ablation,
    }
    mutated_external = build_partition({"G": "derived_in_gate"})
    mutated_derived = build_partition({"fingerprint_27": "external_bridge_input"})
    external_partition_ablation = partition_ablation(mutated_external)
    derived_partition_ablation = partition_ablation(mutated_derived)
    probes["3f_partition_mislabel"] = {
        "external_bridge_input_G_as_derived": {
            "partition_ok": mutated_external["ok"],
            "verdict": "FAIL_PROVENANCE_PARTITION" if not mutated_external["ok"] else "NO_FAIL",
        },
        "derived_fingerprint_27_as_external": {
            "partition_ok": mutated_derived["ok"],
            "verdict": "FAIL_PROVENANCE_PARTITION" if not mutated_derived["ok"] else "NO_FAIL",
        },
        "g_invariance_only_would_miss_54_over_5": partition["g_to_lambda_g_diagnostic"][
            "invariance_only_trap_catches_54_over_5"
        ],
        "self_ablation": {
            **external_partition_ablation,
            "with_mutation_cases": {
                "external_bridge_input_G_as_derived": external_partition_ablation["with_mutation"],
                "derived_fingerprint_27_as_external": derived_partition_ablation["with_mutation"],
            },
        },
    }
    probes["3g_wrong_prefactor_object"] = {
        "mutated_object": hstr(prefactor["plain_object"]),
        "plain_P2": hstr(prefactor["plain_N_over_D"]["P2"]),
        "correct_P2": hstr(prefactor["expected"]["P2"]),
        "plain_equals_correct_P2": prefactor["self_check"]["plain_equals_correct_P2"],
        "verdict": "FAIL_PREFACTOR_ALGEBRA"
        if not prefactor["self_check"]["plain_equals_correct_P2"]
        else "NO_FAIL",
        "self_ablation": prefactor_ablation,
    }

    probe_flags: dict[str, bool] = {}
    for name, probe in probes.items():
        if name == "3a_wrong_bc":
            fires = (
                probe["incoming"]["verdict"] != "NO_FAIL"
                and probe["standing_j2"]["verdict"] != "NO_FAIL"
            )
        elif name == "3f_partition_mislabel":
            fires = (
                probe["external_bridge_input_G_as_derived"]["verdict"] != "NO_FAIL"
                and probe["derived_fingerprint_27_as_external"]["verdict"] != "NO_FAIL"
            )
        else:
            fires = probe["verdict"] != baseline_verdict and probe["verdict"] != "NO_FAIL"
        probe_flags[name] = bool(fires and probe["self_ablation"]["fail_suppressed"])
    able_to_fail_ok = all(probe_flags.values())
    neutered_each_probe_flips_aggregate_false = all(
        not all(value if other != name else False for other, value in probe_flags.items())
        for name in probe_flags
    )
    return {
        "probes": probes,
        "probe_flags": probe_flags,
        "able_to_fail_ok": able_to_fail_ok,
        "neutered_each_probe_flips_aggregate_false": neutered_each_probe_flips_aggregate_false,
    }


def engine_summary(ctx: dict[str, Any]) -> dict[str, Any]:
    fp = ctx["fingerprint"]
    pref = ctx["prefactor"]
    dim = ctx["dimensions"]
    sample = {a: 3.0, c_s: 2.0}
    return {
        "schema": "pathA_33_sympy_scratch/v1",
        "engine": "sympy",
        "fingerprint": {
            "lambda_out_series_z": hstr(fp["outgoing"]["lambda_series"]),
            "Yhat_out_series_z": hstr(fp["outgoing"]["yhat_series"]),
            "coefficients_z": {
                "u2": hstr(fp["outgoing"]["u2_z"]),
                "u4": hstr(fp["outgoing"]["u4_z"]),
                "v5": hstr(fp["outgoing"]["v5_z"]),
            },
            "coefficients_physical": {
                "u2": hstr(fp["outgoing"]["u2"]),
                "u4": hstr(fp["outgoing"]["u4"]),
                "v5": hstr(fp["outgoing"]["v5"]),
            },
            "coefficients_z_numeric": {
                "u2": numeric(fp["outgoing"]["u2_z"], sample),
                "u4": numeric(fp["outgoing"]["u4_z"], sample),
                "v5": numeric(fp["outgoing"]["v5_z"], sample),
            },
            "coefficient_matches": fp["matches"],
            "chi_Q": hstr(fp["chi_Q"]),
            "chi_Q_numeric": fp["samples"]["chi_Q"],
            "incoming_chi_Q": hstr(fp["chi_Q_incoming"]),
            "incoming_chi_Q_numeric": fp["samples"]["incoming_chi_Q"],
            "standing_lambda_static": hstr(fp["standing"]["lambda_series"].coeff(z, 0)),
            "standing_Yhat_static": hstr(fp["standing"]["static"]),
            "standing_v5_z": hstr(fp["standing"]["v5_z"]),
        },
        "prefactor": {
            "correct_series": hstr(pref["correct_series"]),
            "coefficients": {name: hstr(value) for name, value in pref["coefficients"].items()},
            "residuals_to_formula": {name: hstr(value) for name, value in pref["residuals"].items()},
            "matches": pref["matches"],
            "plain_P2": hstr(pref["plain_N_over_D"]["P2"]),
            "plain_equals_correct_P2": pref["self_check"]["plain_equals_correct_P2"],
            "sample_values": pref["sample_values"],
        },
        "dimension": {
            "P0_raw_dimension": dim["P0_raw_dimension"],
            "frequency_normalization_dimension": dim["frequency_normalization_dimension"],
            "P0_dimension": dim["P0_dimension"],
            "P0_physical_dimension": dim["P0_physical_dimension"],
            "dimensional_ok": dim["dimensional_ok"],
            "mu_hat0_dimension": dim["mu_hat0_dimension"],
            "Gamma5_dimension": dim["Gamma5_dimension"],
            "lhs_dimension": dim["lhs_dimension"],
            "rhs_dimension": dim["rhs_dimension"],
            "mu_hat0_homogeneity_pass": dim["mu_hat0_homogeneity_diagnostic"]["homogeneity_pass"],
            "mu_hat0_homogeneity_label": dim["mu_hat0_homogeneity_diagnostic"]["label"],
            "drop_norm_verdict": dim["mutation_drop_frequency_normalization"]["verdict"],
            "corrupt_N0_verdict": dim["mutation_corrupt_N0_dimension"]["verdict"],
        },
        "headline_booleans": {
            "fingerprint_ok": fp["ok"],
            "prefactor_ok": pref["ok"],
            "dimension_ok": dim["dimensional_ok"],
        },
    }


def build_context() -> dict[str, Any]:
    fingerprint = build_fingerprint()
    prefactor = build_prefactor()
    scaling = build_scaling()
    dimensions = build_dimensions(fingerprint)
    equivalence = build_equivalence(fingerprint)
    partition = build_partition()
    passivity = passivity_from_source(fingerprint["outgoing"])
    preliminary_gates = {
        "fingerprint_ok": fingerprint["ok"],
        "prefactor_ok": prefactor["ok"],
        "scaling_ok": scaling["target_scaling_ok"],
        "equivalence_ok": equivalence["ok"],
        "dimensional_ok": dimensions["dimensional_ok"],
        "outgoing_ok": passivity["genuine_outgoing"],
        "provenance_ok": partition["ok"],
        "able_to_fail_ok": True,
    }
    preliminary_verdict = base_verdict(preliminary_gates, partition)
    counterfactuals = build_counterfactuals(
        fingerprint,
        prefactor,
        scaling,
        dimensions,
        equivalence,
        partition,
        preliminary_gates,
        preliminary_verdict,
    )
    gates = dict(preliminary_gates)
    gates["able_to_fail_ok"] = counterfactuals["able_to_fail_ok"]
    verdict = base_verdict(gates, partition)
    return {
        "fingerprint": fingerprint,
        "port_moments": build_port_moments(),
        "prefactor": prefactor,
        "scaling": scaling,
        "dimensions": dimensions,
        "equivalence": equivalence,
        "partition": partition,
        "passivity": passivity,
        "counterfactuals": counterfactuals,
        "gate_booleans": gates,
        "verdict": verdict,
        "which_rung": verdict,
    }


def sympify_rational(text: str) -> sp.Expr:
    return sp.sympify(str(text).replace("^", "**"))


def compare_engines(sym: dict[str, Any], mma: dict[str, Any]) -> dict[str, Any]:
    symbolic_deltas: list[float] = []
    symbolic_checks: dict[str, bool] = {}
    for name in ("u2", "u4", "v5"):
        diff = compact(
            sympify_rational(sym["fingerprint"]["coefficients_z"][name])
            - sympify_rational(mma["fingerprint"]["coefficients_z"][name])
        )
        ok = diff == 0
        symbolic_checks[f"fingerprint_{name}_z"] = ok
        symbolic_deltas.append(0.0 if ok else 1.0)
    chi_diff = compact(
        sympify_rational(sym["fingerprint"]["chi_Q"]) - sympify_rational(mma["fingerprint"]["chi_Q"])
    )
    symbolic_checks["chi_Q"] = chi_diff == 0
    symbolic_deltas.append(0.0 if chi_diff == 0 else 1.0)
    for name in ("P0", "P2", "P4"):
        sym_zero = sym["prefactor"]["residuals_to_formula"][name] == "0"
        mma_zero = mma["prefactor"]["residuals_to_formula"][name] == "0"
        symbolic_checks[f"prefactor_{name}_residual_zero_in_both"] = sym_zero and mma_zero
        symbolic_deltas.append(0.0 if sym_zero and mma_zero else 1.0)
    for name in (
        "P0_raw_dimension",
        "frequency_normalization_dimension",
        "P0_dimension",
        "P0_physical_dimension",
        "dimensional_ok",
        "mu_hat0_dimension",
        "Gamma5_dimension",
        "lhs_dimension",
        "rhs_dimension",
        "drop_norm_verdict",
        "corrupt_N0_verdict",
    ):
        ok = sym["dimension"][name] == mma["dimension"][name]
        symbolic_checks[f"dimension_{name}"] = ok
        symbolic_deltas.append(0.0 if ok else 1.0)

    numeric_pairs: dict[str, tuple[float, float]] = {}
    for name in ("u2", "u4", "v5"):
        numeric_pairs[f"fingerprint_{name}_z"] = (
            float(sym["fingerprint"]["coefficients_z_numeric"][name]),
            float(mma["fingerprint"]["coefficients_z_numeric"][name]),
        )
    for name in ("chi_Q", "incoming_chi_Q"):
        numeric_pairs[name] = (
            float(sym["fingerprint"][f"{name}_numeric"]),
            float(mma["fingerprint"][f"{name}_numeric"]),
        )
    for name in ("P0", "P2", "P4", "plain_P2", "plain_minus_correct_P2"):
        numeric_pairs[f"prefactor_{name}"] = (
            float(sym["prefactor"]["sample_values"][name]),
            float(mma["prefactor"]["sample_values"][name]),
        )
    numeric_deltas = {
        name: abs(left - right) for name, (left, right) in numeric_pairs.items()
    }
    max_symbolic_delta = max(symbolic_deltas) if symbolic_deltas else 0.0
    max_numeric_delta = max(numeric_deltas.values()) if numeric_deltas else 0.0
    status = (
        "pass"
        if max_symbolic_delta < SYMBOLIC_TOL
        and max_numeric_delta < NUMERIC_TOL
        and all(symbolic_checks.values())
        else "FAIL_ENGINE_DISAGREE"
    )
    return {
        "status": status,
        "max_symbolic_delta": max_symbolic_delta,
        "symbolic_tolerance": SYMBOLIC_TOL,
        "max_numeric_delta": max_numeric_delta,
        "numeric_tolerance": NUMERIC_TOL,
        "symbolic_checks": symbolic_checks,
        "numeric_deltas": numeric_deltas,
        "compared_quantities": list(numeric_pairs) + list(symbolic_checks),
    }


def stringify_expr_tree(obj: Any) -> Any:
    if isinstance(obj, sp.Basic):
        return hstr(obj)
    if isinstance(obj, dict):
        return {key: stringify_expr_tree(value) for key, value in obj.items()}
    if isinstance(obj, list):
        return [stringify_expr_tree(value) for value in obj]
    if isinstance(obj, tuple):
        return [stringify_expr_tree(value) for value in obj]
    return obj


def build_final_payload(ctx: dict[str, Any], sym_engine: dict[str, Any], mma_engine: dict[str, Any]) -> dict[str, Any]:
    engine_agreement = compare_engines(sym_engine, mma_engine)
    verdict = ctx["verdict"]
    which_rung = ctx["which_rung"]
    if engine_agreement["status"] != "pass":
        verdict = "QUAD_FAIL_ENGINE_DISAGREE"
        which_rung = "FAIL_ENGINE_DISAGREE"

    fp = ctx["fingerprint"]
    pref = ctx["prefactor"]
    payload = {
        "schema": "pathA_33_quadrupole_normalization/v1",
        "engine": "dual-engine",
        "source_directive": "software/stage1_solver/directives/pathA_33_quadrupole_normalization.md",
        "cited_sources": [
            "research/pde_ledger/notes/stages/moving_throat_pde_handoff_full.md sections 10.3, 11, 12",
            "software/stage1_solver/reports/pathA_22a_dimensional_skeleton.md",
        ],
        "verdict": verdict,
        "which_rung": which_rung,
        "gate_booleans": ctx["gate_booleans"],
        "dtn_fingerprint": {
            "definition": {
                "Lambda2_out": "z*h2^(1)'(z)/h2^(1)(z)",
                "Yhat2_out": "-3/Lambda2_out",
                "z": "a*omega/c_s",
            },
            "sympy": sym_engine["fingerprint"],
            "mathematica": mma_engine["fingerprint"],
            "derived_coefficients": {
                "u2": hstr(fp["outgoing"]["u2"]),
                "u4": hstr(fp["outgoing"]["u4"]),
                "v5": hstr(fp["outgoing"]["v5"]),
                "chi_Q": hstr(fp["chi_Q"]),
            },
        },
        "outgoing_port_moments": stringify_expr_tree(ctx["port_moments"]),
        "prefactor_algebra": {
            "definition": "P(omega)=D0*N(omega)/Dcons(omega)^2",
            "Nomega": hstr(N0 + N2 * omega**2 + N4 * omega**4),
            "Dcons": hstr(D0 + D2 * omega**2 + D4 * omega**4),
            "correct_series": hstr(pref["correct_series"]),
            "coefficients": {name: hstr(value) for name, value in pref["coefficients"].items()},
            "expected": {name: hstr(value) for name, value in pref["expected"].items()},
            "residuals_to_formula": {name: hstr(value) for name, value in pref["residuals"].items()},
            "matches": pref["matches"],
            "N_over_D_self_check": stringify_expr_tree(pref["self_check"]),
            "plain_N_over_D_coefficients": {
                name: hstr(value) for name, value in pref["plain_N_over_D"].items()
            },
            "sympy_sample_values": pref["sample_values"],
            "mathematica_sample_values": mma_engine["prefactor"]["sample_values"],
        },
        "target_scaling": stringify_expr_tree(ctx["scaling"]),
        "equivalence": ctx["equivalence"],
        "dimensional_check": ctx["dimensions"],
        "dimensional_table": ctx["dimensions"]["table"],
        "provenance_partition": ctx["partition"],
        "counterfactuals": ctx["counterfactuals"]["probes"],
        "able_to_fail": {
            "probe_flags": ctx["counterfactuals"]["probe_flags"],
            "able_to_fail_ok": ctx["counterfactuals"]["able_to_fail_ok"],
            "neutered_each_probe_flips_aggregate_false": ctx["counterfactuals"][
                "neutered_each_probe_flips_aggregate_false"
            ],
        },
        "passivity": ctx["passivity"],
        "engine_agreement": engine_agreement,
        "reduction_certificate": {
            "frozen_inputs": [
                "Gate-3 isotropic linearized reference",
                "symbolic D_n",
                "symbolic outgoing port scalars and N_n",
                "external calibration knob G",
            ],
            "computed": [
                "DtN-derived outgoing fingerprint",
                "P0/P2/P4 prefactor algebra",
                "P0_target_scaling=a^-5",
                "Gamma5/chi_Q equivalence",
                "real-expression dimensional table",
                "provenance partition",
                "counterfactual verdict flips",
            ],
            "deferred": [
                "numerical D_n,N_n and port scalars from Gate 6",
                "actual branch a-scaling from Gate 6",
                "derivation of G and the literal 54/5 magnitude",
                "cross-l reconciliation and downstream PN match-back",
            ],
        },
    }
    return stringify_expr_tree(payload)


def build_report(payload: dict[str, Any]) -> str:
    fp = payload["dtn_fingerprint"]
    pref = payload["prefactor_algebra"]
    dim = payload["dimensional_check"]
    part = payload["provenance_partition"]
    eng = payload["engine_agreement"]
    probes = payload["counterfactuals"]
    lines = [
        "# PathA-33 quadrupole normalization result",
        "",
        f"Computed verdict: `{payload['verdict']}` (`{payload['which_rung']}`).",
        "",
        "The earned part is the outgoing l=2 fingerprint shape, the squared-denominator prefactor algebra, "
        "the a^-5 target scaling, the Gamma5/chi_Q equivalence, and the restored-units dimensional closure. "
        "The magnitude remains calibrated: G and the PN 2/5 input keep the assembled 54/5 on the "
        "`external_bridge_input` rung.",
        "",
        "## Derived fingerprint",
        "",
        f"- SymPy: u2=`{fp['sympy']['coefficients_physical']['u2']}`, "
        f"u4=`{fp['sympy']['coefficients_physical']['u4']}`, "
        f"v5=`{fp['sympy']['coefficients_physical']['v5']}`.",
        f"- Mathematica: u2=`{fp['mathematica']['coefficients_physical']['u2']}`, "
        f"u4=`{fp['mathematica']['coefficients_physical']['u4']}`, "
        f"v5=`{fp['mathematica']['coefficients_physical']['v5']}`.",
        f"- Derived chi_Q: `{fp['derived_coefficients']['chi_Q']}`; incoming gives "
        f"`{fp['sympy']['incoming_chi_Q']}`.",
        "",
        "## Prefactor algebra",
        "",
        f"- P0=`{pref['coefficients']['P0']}`.",
        f"- P2=`{pref['coefficients']['P2']}`.",
        f"- P4=`{pref['coefficients']['P4']}`.",
        f"- N/D self-check: plain P2=`{pref['plain_N_over_D_coefficients']['P2']}` versus "
        f"correct P2=`{pref['expected']['P2']}`; the missing term is `-D2*N0` versus `-2*D2*N0`.",
        "",
        "## Dimensional result",
        "",
        f"- Mu-free gate: `[P0_raw]` = `{dim['P0_raw_dimension']}`, "
        f"`[(c_s/a)^2]` = `{dim['frequency_normalization_dimension']}`, "
        f"`[P0_phys]` = `{dim['P0_physical_dimension']}`; dimensional_ok = "
        f"`{dim['dimensional_ok']}`.",
        f"- `mu_hat0` diagnostic: `{dim['mu_hat0_homogeneity_diagnostic']['label']}`; "
        f"`[mu_hat0]` = `{dim['mu_hat0_dimension']}`, LHS/RHS = "
        f"`{dim['lhs_dimension']}` / `{dim['rhs_dimension']}`, diagnostic pass = "
        f"`{dim['mu_hat0_homogeneity_diagnostic']['homogeneity_pass']}`.",
        f"- Section 3d drop-normalization probe verdict: "
        f"`{probes['3d_dimensional_break']['verdict']}`.",
        f"- Section 3d' corrupt-port-dimension probe verdict: "
        f"`{probes['3d_prime_corrupt_port_dimension']['verdict']}`.",
        "",
        "## Provenance",
        "",
        f"- Decomposition: `{part['decomposition_54_over_5']['identity']}`; "
        f"earned factor class = `{part['decomposition_54_over_5']['earned_factor']['class']}`, "
        f"calibrated factor class = `{part['decomposition_54_over_5']['calibrated_factor']['class']}`.",
        f"- Assembled 54/5 class: `{part['decomposition_54_over_5']['assembled_magnitude_class']}`.",
        "",
        "## Probe verdicts",
        "",
        f"- 3a incoming: `{probes['3a_wrong_bc']['incoming']['verdict']}`; "
        f"standing: `{probes['3a_wrong_bc']['standing_j2']['verdict']}`.",
        f"- 3b imposed dissipation: `{probes['3b_imposed_dissipation']['verdict']}`.",
        f"- 3c wrong scaling: `{probes['3c_wrong_scaling']['verdict']}`.",
        f"- 3d dimensional break: `{probes['3d_dimensional_break']['verdict']}`.",
        f"- 3d' corrupt port dimension: `{probes['3d_prime_corrupt_port_dimension']['verdict']}`.",
        f"- 3e equivalence break: `{probes['3e_equivalence_break']['verdict']}`.",
        f"- 3f partition mislabels: "
        f"`{probes['3f_partition_mislabel']['external_bridge_input_G_as_derived']['verdict']}`, "
        f"`{probes['3f_partition_mislabel']['derived_fingerprint_27_as_external']['verdict']}`.",
        f"- 3g wrong prefactor object: `{probes['3g_wrong_prefactor_object']['verdict']}`.",
        "",
        "## Engine agreement",
        "",
        f"- Status: `{eng['status']}`.",
        f"- Max symbolic delta: `{eng['max_symbolic_delta']}` (tol `{eng['symbolic_tolerance']}`).",
        f"- Max numeric delta: `{eng['max_numeric_delta']}` (tol `{eng['numeric_tolerance']}`).",
        "",
        "Deferred: the numerical branch data `(D_n,N_n)`, port scalars, actual branch a-scaling, "
        "cross-l reconciliation, and derivation of G/the magnitude remain outside Gate 4.",
    ]
    return "\n".join(lines) + "\n"


def build_feed_note(payload: dict[str, Any]) -> str:
    probes = payload["counterfactuals"]
    eng = payload["engine_agreement"]
    dim = payload["dimensional_check"]
    fp = payload["dtn_fingerprint"]["derived_coefficients"]
    return "\n".join(
        [
            "# PathA-33 quadrupole normalization feed note",
            "",
            f"- Verdict: `{payload['verdict']}` (`{payload['which_rung']}`).",
            f"- DtN fingerprint: u2=`{fp['u2']}`, u4=`{fp['u4']}`, v5=`{fp['v5']}`, chi_Q=`{fp['chi_Q']}`.",
            "- Prefactor object: `P(omega)=D0*N(omega)/Dcons(omega)^2`; plain `N/D` fails the P2 factor-of-two self-check.",
            f"- Dimensions: mu-free gate `[P0_raw]={dim['P0_raw_dimension']}`, "
            f"`[(c_s/a)^2]={dim['frequency_normalization_dimension']}`, "
            f"`[P0_phys]={dim['P0_physical_dimension']}`, dimensional_ok=`{dim['dimensional_ok']}`; "
            f"drop-normalization probe=`{probes['3d_dimensional_break']['verdict']}`, "
            f"corrupt-port-dimension probe=`{probes['3d_prime_corrupt_port_dimension']['verdict']}`.",
            f"- Diagnostic only: `[mu_hat0]={dim['mu_hat0_dimension']}`, "
            f"{dim['mu_hat0_homogeneity_diagnostic']['label']}.",
            "- Provenance: fingerprint 27 is `derived_in_gate`; PN 2/5, G, and the assembled 54/5 magnitude are `external_bridge_input`.",
            f"- Engine agreement: `{eng['status']}`, max symbolic delta `{eng['max_symbolic_delta']}`, "
            f"max numeric delta `{eng['max_numeric_delta']}`.",
            "- Deferred to Gate 6: numerical branch data `(D_n,N_n)`, port scalars, actual branch scaling, and on-solution branch selection.",
            "",
        ]
    )


def main() -> int:
    ctx = build_context()
    sym_engine = engine_summary(ctx)
    yaml_write(SYM_YAML, sym_engine)

    mma_engine = yaml_read(MMA_YAML)
    if mma_engine is None:
        print(f"Wrote SymPy scratch: {SYM_YAML}")
        print(f"Mathematica scratch not present yet: {MMA_YAML}")
        return 0

    payload = build_final_payload(ctx, sym_engine, mma_engine)
    yaml_write(RESULTS_YAML, payload)
    REPORT_MD.write_text(build_report(payload), encoding="utf-8")
    FEED_NOTE.write_text(build_feed_note(payload), encoding="utf-8")

    if payload["verdict"] != "QUAD_CALIBRATED":
        print(f"FAIL: unexpected verdict {payload['verdict']}")
        return 1
    if payload["engine_agreement"]["status"] != "pass":
        print("FAIL: engine disagreement")
        return 1
    print(f"Wrote final YAML: {RESULTS_YAML}")
    print(f"Wrote report: {REPORT_MD}")
    print(f"Wrote feed note: {FEED_NOTE}")
    print(f"verdict={payload['verdict']} max_numeric_delta={payload['engine_agreement']['max_numeric_delta']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
