#!/usr/bin/env python3
"""PathA-34 Gate 5 cross-ell unification, SymPy engine.

Run order:

  timeout 600 python software/stage1_solver/tools/pathA_34_cross_l_unification_sympy.py
  timeout 600 math -script software/stage1_solver/tools/pathA_34_cross_l_unification.wl
  timeout 600 python software/stage1_solver/tools/pathA_34_cross_l_unification_sympy.py

The first Python run writes the SymPy scratch lane.  Mathematica writes an
independent scratch lane.  The second Python run compares the engines and
emits the final YAML and Markdown reports.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path
from typing import Any

import sympy as sp
import yaml


class NoAliasDumper(yaml.SafeDumper):
    def ignore_aliases(self, data: Any) -> bool:
        return True


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
REPORTS = STAGE1_ROOT / "reports"
SCRATCH = STAGE1_ROOT / "_scratch"

SYM_YAML = SCRATCH / "pathA_34_sympy_results.yaml"
MMA_YAML = SCRATCH / "pathA_34_mathematica_results.yaml"
RESULTS_YAML = REPORTS / "pathA_34_results.yaml"
REPORT_MD = REPORTS / "pathA_34_cross_l_unification.md"

SYMBOLIC_TOL = 1.0e-10
NUMERIC_TOL = 1.0e-8

z, omega = sp.symbols("z omega")
a, c_s = sp.symbols("a c_s", positive=True)
M0, D1 = sp.symbols("M0 D1", nonzero=True)
R0, R1 = sp.symbols("R0 R1")
A0_sym, A1_sym = sp.symbols("A0 A1")
N0, D0 = sp.symbols("N0 D0", nonzero=True)
N2, D2, N4, D4 = sp.symbols("N2 D2 N4 D4")
K0c, Keta, TOmega = sp.symbols("K0c K_eta T_Omega", positive=True)
Z0ret, Z1ret = sp.symbols("Z0_ret Z1_ret", positive=True)
q_free = sp.symbols("q_free")
gain0, gain1 = sp.symbols("gain0 gain1")
eta_null = sp.symbols("eta_null", positive=True)

OmegaU, OmegaW, gU, gW, Rmix, Delta, Sport = sp.symbols(
    "Omega_U Omega_W g_U g_W R_mix Delta S_port", nonzero=True
)

Ldim, Mdim, Tdim = sp.symbols("L M T", positive=True)
Dim = tuple[sp.Rational, sp.Rational, sp.Rational]
ZERO_DIM: Dim = (sp.Rational(0), sp.Rational(0), sp.Rational(0))


@dataclass(frozen=True)
class Mutation:
    name: str = "baseline"
    decouple_knobs: bool = False
    inject_null: bool = False
    wrong_sign_return: bool = False
    perfect_return: bool = False
    break_gate4: bool = False
    corrupt_dimension: bool = False
    assert_not_derive: bool = False
    no_consistent_return: bool = False
    selector_equation_set: str = "none"
    corrupt_port_kernel: bool = False


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


def bool_zero(expr: sp.Expr) -> bool:
    return compact(expr) == 0


def numeric(expr: sp.Expr, subs: dict[sp.Symbol, float | int]) -> float:
    return float(sp.N(compact(expr).subs(subs), 40))


def numeric_complex_parts(expr: sp.Expr, subs: dict[sp.Symbol, float | int]) -> dict[str, float]:
    value = complex(sp.N(compact(expr).subs(subs), 40))
    return {"re": float(value.real), "im": float(value.imag)}


def series_no_o(expr: sp.Expr, var: sp.Symbol, order: int) -> sp.Expr:
    return sp.expand(sp.series(expr, var, 0, order).removeO())


def spherical_j(lval: int) -> sp.Expr:
    if lval == 0:
        return sp.sin(z) / z
    if lval == 1:
        return sp.sin(z) / z**2 - sp.cos(z) / z
    if lval == 2:
        return (3 / z**3 - 1 / z) * sp.sin(z) - 3 * sp.cos(z) / z**2
    raise ValueError(lval)


def spherical_y(lval: int) -> sp.Expr:
    if lval == 0:
        return -sp.cos(z) / z
    if lval == 1:
        return -sp.cos(z) / z**2 - sp.sin(z) / z
    if lval == 2:
        return (1 / z - 3 / z**3) * sp.cos(z) - 3 * sp.sin(z) / z**2
    raise ValueError(lval)


@lru_cache(maxsize=None)
def dtn_branch(lval: int, kind: str) -> dict[str, Any]:
    j = spherical_j(lval)
    y = spherical_y(lval)
    if kind == "outgoing_hankel1":
        h = j + sp.I * y
        source = "hankel1_outgoing_for_exp_minus_i_omega_t"
    elif kind == "incoming_hankel2":
        h = j - sp.I * y
        source = "hankel2_incoming_for_exp_minus_i_omega_t"
    elif kind == "standing_j":
        h = j
        source = "standing_regular_j_l"
    else:
        raise ValueError(kind)

    lam = compact(z * sp.diff(h, z) / h)
    yout = compact(-(lval + 1) / lam)
    series_order = max(8, 2 * lval + 4)
    lam_series = series_no_o(lam, z, series_order)
    y_series = series_no_o(yout, z, series_order)
    radiative_power = 2 * lval + 1
    radiative_coeff = compact(y_series.coeff(z, radiative_power) / sp.I)
    return {
        "ell": lval,
        "kind": kind,
        "source": source,
        "normalization_factor": -(lval + 1),
        "h": h,
        "lambda": lam,
        "Y": yout,
        "lambda_series": lam_series,
        "Y_series": y_series,
        "static": compact(y_series.coeff(z, 0)),
        "radiative_power": radiative_power,
        "radiative_coeff_z": radiative_coeff,
        "raw_outgoing_order": radiative_power if not bool_zero(radiative_coeff) else None,
        "u2_z": compact(y_series.coeff(z, 2)),
        "u4_z": compact(y_series.coeff(z, 4)),
        "v5_z": compact(y_series.coeff(z, 5) / sp.I),
    }


@lru_cache(maxsize=1)
def build_fingerprints() -> dict[str, Any]:
    expected_radiative = {0: sp.Integer(1), 1: sp.Rational(1, 2), 2: sp.Rational(1, 27)}
    expected_order = {0: 1, 1: 3, 2: 5}
    outgoing = {ell: dtn_branch(ell, "outgoing_hankel1") for ell in (0, 1, 2)}
    incoming = {ell: dtn_branch(ell, "incoming_hankel2") for ell in (0, 1, 2)}
    matches: dict[str, bool] = {}
    for ell in (0, 1, 2):
        matches[f"ell{ell}_normalization"] = outgoing[ell]["normalization_factor"] == -(ell + 1)
        matches[f"ell{ell}_static"] = bool_zero(outgoing[ell]["static"] - 1)
        matches[f"ell{ell}_radiative_coeff"] = bool_zero(
            outgoing[ell]["radiative_coeff_z"] - expected_radiative[ell]
        )
        matches[f"ell{ell}_raw_order"] = outgoing[ell]["raw_outgoing_order"] == expected_order[ell]
        matches[f"ell{ell}_incoming_flips_radiative_sign"] = bool_zero(
            incoming[ell]["radiative_coeff_z"] + outgoing[ell]["radiative_coeff_z"]
        )
    matches["ell2_u2"] = bool_zero(outgoing[2]["u2_z"] - sp.Rational(1, 9))
    matches["ell2_u4"] = bool_zero(outgoing[2]["u4_z"] - sp.Rational(4, 81))
    matches["ell2_v5"] = bool_zero(outgoing[2]["v5_z"] - sp.Rational(1, 27))
    chi_q = compact(outgoing[2]["v5_z"] / sp.Rational(1, 27))
    matches["chi_Q"] = bool_zero(chi_q - 1)
    return {
        "time_convention": "exp(-i*omega*t); h_l^(1) is outgoing and h_l^(2) is incoming",
        "outgoing": outgoing,
        "incoming": incoming,
        "standing": {},
        "expected_radiative_coefficients": expected_radiative,
        "expected_raw_orders": expected_order,
        "matches": matches,
        "ok": all(matches.values()),
        "chi_Q": chi_q,
    }


def build_gate4_non_regression(
    fingerprints: dict[str, Any], port: dict[str, Any], mutation: Mutation
) -> dict[str, Any]:
    out2 = fingerprints["incoming"][2] if mutation.break_gate4 else fingerprints["outgoing"][2]
    n0_eff = port["N0_from_port"]
    N0eff = sp.symbols("N0_eff", nonzero=True)
    Nomega = N0eff + N2 * omega**2 + N4 * omega**4
    Dcons = D0 + D2 * omega**2 + D4 * omega**4
    correct_obj = D0 * Nomega / Dcons**2
    plain_obj = Nomega / Dcons
    obj = plain_obj if mutation.break_gate4 else correct_obj
    series = series_no_o(obj, omega, 6)
    coeffs = {
        "P0": compact(series.coeff(omega, 0).subs(N0eff, n0_eff)),
        "P2": compact(series.coeff(omega, 2).subs(N0eff, n0_eff)),
        "P4": compact(series.coeff(omega, 4).subs(N0eff, n0_eff)),
    }
    expected = {
        "P0": n0_eff / D0,
        "P2": (D0 * N2 - 2 * D2 * n0_eff) / D0**2,
        "P4": (
            D0**2 * N4
            - 2 * D0 * (D2 * N2 + D4 * n0_eff)
            + 3 * D2**2 * n0_eff
        )
        / D0**3,
    }
    residuals = {name: compact(coeffs[name] - expected[name]) for name in coeffs}
    matches = {name: bool_zero(residuals[name]) for name in residuals}
    fingerprint_matches = {
        "u2": bool_zero(out2["u2_z"] - sp.Rational(1, 9)),
        "u4": bool_zero(out2["u4_z"] - sp.Rational(4, 81)),
        "v5": bool_zero(out2["v5_z"] - sp.Rational(1, 27)),
        "chi_Q": bool_zero(compact(out2["v5_z"] / sp.Rational(1, 27)) - 1),
    }
    ok = all(matches.values()) and all(fingerprint_matches.values())
    return {
        "branch_used": out2["kind"],
        "prefactor_object": obj.subs(N0eff, n0_eff),
        "correct_object": correct_obj.subs(N0eff, n0_eff),
        "plain_object": plain_obj.subs(N0eff, n0_eff),
        "series": series.subs(N0eff, n0_eff),
        "coefficients": coeffs,
        "expected": expected,
        "residuals": residuals,
        "prefactor_matches": matches,
        "fingerprint_matches": fingerprint_matches,
        "chi_Q": compact(out2["v5_z"] / sp.Rational(1, 27)),
        "ok": ok,
    }


def build_port_kernel() -> dict[str, Any]:
    return build_port_kernel_for(Mutation())


def build_port_kernel_for(mutation: Mutation) -> dict[str, Any]:
    omega_u = 2 * OmegaU if mutation.corrupt_port_kernel else OmegaU
    p_port = omega_u**2 * gW + Rmix * gU
    delta_port = compact(omega_u**2 * OmegaW**2 - Rmix**2)
    n0_port = compact(p_port**2 / delta_port**2)
    d0_port = D0
    p0_raw = compact(n0_port / d0_port)
    p0_physical = compact((c_s / a) ** 2 * p0_raw)
    return {
        "sector": "ell=2",
        "route": "handoff sections 9.4 and 10.2-10.3 grouped-P2 port kernel",
        "Omega_U_used": omega_u,
        "P_port": p_port,
        "Delta_port": delta_port,
        "N0_from_port": n0_port,
        "D0_from_port_consistent_conservative_scalar": d0_port,
        "D0_source": "Gate-4 conservative D0 branch scalar carried on the grouped-P2 isotropic branch",
        "P0_raw": p0_raw,
        "P0_physical": p0_physical,
        "free_N0_symbol_used_in_verdict": False,
        "corrupt_port_kernel_symbol": mutation.corrupt_port_kernel,
        "derived_non_asserted": True,
    }


def selector_equations(mutation: Mutation) -> list[sp.Equality]:
    if mutation.selector_equation_set == "none":
        return []
    if mutation.selector_equation_set in {
        "derived_pde_admissibility",
        "asserted_unproven",
    }:
        k1 = Keta + 2 * TOmega
        return [sp.Eq(Z0ret, K0c), sp.Eq(Z1ret, k1)]
    raise ValueError(f"unknown selector equation set: {mutation.selector_equation_set}")


def selector_subs(mutation: Mutation) -> dict[sp.Symbol, sp.Expr]:
    subs: dict[sp.Symbol, sp.Expr] = {}
    for equation in selector_equations(mutation):
        lhs, rhs = equation.lhs, equation.rhs
        if not isinstance(lhs, sp.Symbol):
            raise ValueError(f"selector lhs must be a symbol: {equation}")
        subs[lhs] = rhs
    return subs


def equation_text(equation: sp.Equality) -> str:
    return f"{hstr(equation.lhs)} = {hstr(equation.rhs)}"


def selector_provenance(mutation: Mutation) -> dict[str, Any]:
    equations = selector_equations(mutation)
    if mutation.selector_equation_set == "none":
        status = "absent"
        derived = False
        tautological = False
        source = "none available at Gate 5 linear level"
    elif mutation.selector_equation_set == "derived_pde_admissibility":
        status = "present"
        derived = True
        tautological = False
        source = "counterfactual Gate-6 PDE admissibility equation supplied to the able-to-fail probe"
    elif mutation.selector_equation_set == "asserted_unproven":
        status = "present"
        derived = False
        tautological = True
        source = "asserted selector with no named PDE/admissibility provenance"
    else:
        raise ValueError(mutation.selector_equation_set)
    return {
        "status": status,
        "present": bool(equations),
        "equations": [equation_text(eq) for eq in equations],
        "derived_from_named_pde_admissibility_equation": derived,
        "tautological_assertion": tautological,
        "source": source,
    }


def positive_bounded_transfer(t_expr: sp.Expr, eps_expr: sp.Expr) -> dict[str, Any]:
    t_s = compact(t_expr)
    eps_s = compact(eps_expr)
    one_minus = compact(1 - t_s)
    eps_positive = eps_s.is_positive is True
    t_positive = t_s.is_positive is True
    one_minus_positive = one_minus.is_positive is True
    finite = not any(den == 0 for den in sp.denom(t_s).as_ordered_factors())
    return {
        "epsilon_positive": eps_positive,
        "T_positive": t_positive,
        "one_minus_T_positive": one_minus_positive,
        "finite_symbolic_expression": finite,
        "T_in_open_unit_interval": t_positive and one_minus_positive,
        "admissible": eps_positive and t_positive and one_minus_positive and finite,
    }


def build_transfers(mutation: Mutation) -> dict[str, Any]:
    k0 = K0c
    k1 = Keta + 2 * TOmega
    subs = selector_subs(mutation)
    z0 = subs.get(Z0ret, Z0ret)
    z1 = subs.get(Z1ret, Z1ret)
    if mutation.perfect_return:
        z0 = sp.Integer(0)
        z1 = sp.Integer(0)
    if mutation.wrong_sign_return:
        z0 = -subs.get(Z0ret, Z0ret)
        z1 = -subs.get(Z1ret, Z1ret)
    if mutation.no_consistent_return:
        z0 = -2 * K0c
        z1 = -2 * (Keta + 2 * TOmega)
    if mutation.inject_null:
        z0 = compact(z0 + eta_null * K0c)
        z1 = compact(z1 + eta_null * k1)

    t0 = compact(k0 / (k0 + z0))
    t1 = compact(k1 / (k1 + z1))
    eps0 = compact(z0 / k0)
    eps1 = compact(z1 / k1)
    if mutation.decouple_knobs:
        t0 = compact(gain0 * t0)
        t1 = compact(gain1 * t1)

    one_minus_t0 = compact(1 - t0)
    one_minus_t1 = compact(1 - t1)
    eps_relation = {
        "T0_equals_1_over_1_plus_eps": bool_zero(t0 - 1 / (1 + eps0)) if not mutation.decouple_knobs else False,
        "T1_equals_1_over_1_plus_eps": bool_zero(t1 - 1 / (1 + eps1)) if not mutation.decouple_knobs else False,
        "one_minus_T0_equals_eps_over_1_plus_eps": bool_zero(
            one_minus_t0 - eps0 / (1 + eps0)
        )
        if not mutation.decouple_knobs
        else False,
        "one_minus_T1_equals_eps_over_1_plus_eps": bool_zero(
            one_minus_t1 - eps1 / (1 + eps1)
        )
        if not mutation.decouple_knobs
        else False,
    }
    bounded0 = positive_bounded_transfer(t0, eps0)
    bounded1 = positive_bounded_transfer(t1, eps1)
    overcancel = bool_zero(eps0) and bool_zero(eps1)
    no_consistent = (
        bool_zero(eps0 + 2)
        and bool_zero(eps1 + 2)
        and not (bounded0["admissible"] and bounded1["admissible"])
    )
    admissible_branch_exists = bounded0["admissible"] and bounded1["admissible"]
    return {
        "ell0": {
            "sector": "ell=0",
            "route": "handoff section 8.2 / Gate-2 collective (delta_a,delta_L) reduction",
            "K_dc": k0,
            "Z_return": z0,
            "T_dc": t0,
            "epsilon_eff": eps0,
            "derived_non_asserted": True,
        },
        "ell1": {
            "sector": "ell=1",
            "route": "handoff section 9.4 harmonic reduction; angular stiffness K_eta+2*T_Omega",
            "K_dc": k1,
            "Z_return": z1,
            "T_dc": t1,
            "epsilon_eff": eps1,
            "derived_non_asserted": True,
        },
        "relations": eps_relation,
        "boundedness": {
            "ell0": bounded0,
            "ell1": bounded1,
            "positive_bounded_domain": admissible_branch_exists,
        },
        "positive_bounded_domain": admissible_branch_exists,
        "overcancel": overcancel,
        "admissible_branch_exists": admissible_branch_exists,
        "no_consistent_return": no_consistent,
        "selector_substitutions": {hstr(k): hstr(v) for k, v in subs.items()},
    }


def build_residuals(
    fingerprints: dict[str, Any], transfers: dict[str, Any], mutation: Mutation
) -> dict[str, Any]:
    v0 = fingerprints["outgoing"][0]["radiative_coeff_z"]
    v1 = fingerprints["outgoing"][1]["radiative_coeff_z"]
    t0 = transfers["ell0"]["T_dc"]
    t1 = transfers["ell1"]["T_dc"]
    eps0 = transfers["ell0"]["epsilon_eff"]
    eps1 = transfers["ell1"]["epsilon_eff"]
    a0 = compact(sp.I * v0 * (a * omega / c_s) * M0 * (1 - t0))
    a1 = compact(sp.I * v1 * (a * omega / c_s) ** 3 * D1 * (1 - t1))
    expected_a0 = compact(sp.I * a * omega * M0 * eps0 / (c_s * (1 + eps0)))
    expected_a1 = compact(sp.I * a**3 * omega**3 * D1 * eps1 / (2 * c_s**3 * (1 + eps1)))
    leading = {
        "A0": a0,
        "A1": a1,
        "expected_A0": expected_a0,
        "expected_A1": expected_a1,
        "A0_residual_to_pathA29_form": compact(a0 - expected_a0),
        "A1_residual_to_pathA29_form": compact(a1 - expected_a1),
    }
    match = {
        "A0_form": bool_zero(leading["A0_residual_to_pathA29_form"]),
        "A1_form": bool_zero(leading["A1_residual_to_pathA29_form"]),
        "A0_order": compact(a0).as_powers_dict().get(omega, None) == 1,
        "A1_order": compact(a1).as_powers_dict().get(omega, None) == 3,
        "positive_bounded": transfers["positive_bounded_domain"],
        "nonzero_epsilon": not transfers["overcancel"],
        "admissible_branch_exists": transfers["admissible_branch_exists"],
    }
    ok = (
        match["A0_form"]
        and match["A1_form"]
        and match["A0_order"]
        and match["A1_order"]
        and match["positive_bounded"]
        and match["nonzero_epsilon"]
        and match["admissible_branch_exists"]
    )
    return {
        "raw_radiative_coefficients": {"ell0": v0, "ell1": v1},
        "leading": leading,
        "pathA_29_comparison": match,
        "ok": ok,
    }


GENERATOR_DOFS = [
    OmegaU,
    OmegaW,
    Rmix,
    gU,
    gW,
    D0,
    K0c,
    Keta,
    TOmega,
    Z0ret,
    Z1ret,
]


def rank_row(expr: sp.Expr, dofs: list[sp.Symbol]) -> list[sp.Expr]:
    return [compact(sp.diff(expr, dof)) for dof in dofs]


def matrix_rank(rows: list[list[sp.Expr]]) -> int:
    if not rows:
        return 0
    return int(sp.Matrix(rows).rank())


def named_constraints(port: dict[str, Any], transfers: dict[str, Any]) -> list[dict[str, Any]]:
    constraints = [
        {
            "name": "ell0_collective_gate2_stiffness",
            "source": "handoff section 8.2 / Gate-2 collective (delta_a,delta_L) reduction",
            "equations": ["K0_dc = K0c", "T0_dc = K0c/(K0c + Z0_ret) when a return admittance is supplied"],
            "fixes": ["K0_dc"],
            "uses_but_does_not_fix": ["Z0_ret"],
            "linearized_fixed_expression": K0c,
        },
        {
            "name": "ell1_section9_4_harmonic_stiffness",
            "source": "handoff section 9.4 harmonic reduction",
            "equations": [
                "K1_dc = K_eta + 2*T_Omega",
                "T1_dc = K1_dc/(K1_dc + Z1_ret) when a return admittance is supplied",
            ],
            "fixes": ["K1_dc"],
            "uses_but_does_not_fix": ["Z1_ret"],
            "linearized_fixed_expression": Keta + 2 * TOmega,
        },
        {
            "name": "ell2_section10_port_kernel",
            "source": "handoff sections 10.2-10.3 grouped-P2 port kernel + Gate-4 conservative D0",
            "equations": [
                f"P_port = {hstr(port['P_port'])}",
                f"Delta_port = {hstr(port['Delta_port'])}",
                f"N0_port = {hstr(port['N0_from_port'])}",
                f"P0_raw = {hstr(port['P0_raw'])}",
            ],
            "fixes": ["P0_raw combination"],
            "uses_but_does_not_fix": [],
            "linearized_fixed_expression": port["P0_raw"],
        },
    ]
    out: list[dict[str, Any]] = []
    for item in constraints:
        expr = sp.sympify(item["linearized_fixed_expression"])
        touches = [hstr(dof) for dof in GENERATOR_DOFS if not bool_zero(sp.diff(expr, dof))]
        row = rank_row(expr, GENERATOR_DOFS)
        out.append(
            {
                key: hstr(value) if isinstance(value, sp.Basic) else value
                for key, value in item.items()
                if key != "linearized_fixed_expression"
            }
            | {
                "touches_generator_dofs": touches,
                "linearized_row": [hstr(value) for value in row],
            }
        )
    return out


def pathA29_premise_citation() -> dict[str, Any]:
    return {
        "source": "software/stage1_solver/reports/pathA_29_results.yaml",
        "Z_is_premise": True,
        "boundary_dof": "none",
        "canonical_form_contains": "Z_is_premise:true|boundary_dof:none",
        "meaning_for_gate5": "the collected Gate-5 linear named equations may use Z_ret in transfer forms, but pathA_29 supplies no boundary degree of freedom or PDE equation that fixes it",
    }


def selector_constraint_rows(mutation: Mutation) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for equation in selector_equations(mutation):
        residual = compact(equation.lhs - equation.rhs)
        rows.append(
            {
                "name": f"selector_{hstr(equation.lhs)}",
                "equation": equation_text(equation),
                "residual_fixed_to_zero": hstr(residual),
                "row": rank_row(residual, GENERATOR_DOFS + ([eta_null] if mutation.inject_null else [])),
            }
        )
    return rows


def build_rank_audit(
    port: dict[str, Any], transfers: dict[str, Any], mutation: Mutation
) -> dict[str, Any]:
    dofs = list(GENERATOR_DOFS)
    if mutation.inject_null:
        dofs.append(eta_null)
    base_constraints = [port["P0_raw"], K0c, Keta + 2 * TOmega]
    selector_rows_exprs = [eq.lhs - eq.rhs for eq in selector_equations(mutation)]
    constraint_exprs = base_constraints + selector_rows_exprs
    rows = [rank_row(expr, dofs) for expr in constraint_exprs]
    rank0 = matrix_rank(rows)
    nullity = len(dofs) - rank0
    t0 = transfers["ell0"]["T_dc"]
    t1 = transfers["ell1"]["T_dc"]
    grad_t0 = rank_row(t0, dofs)
    grad_t1 = rank_row(t1, dofs)
    return_aug_rank = matrix_rank(rows + [grad_t0, grad_t1])
    return_moving_nullity = return_aug_rank - rank0
    native_moves = return_moving_nullity > 0
    untouched_return_dofs = []
    for zret in (Z0ret, Z1ret):
        touched = any(not bool_zero(sp.diff(expr, zret)) for expr in constraint_exprs)
        if not touched:
            untouched_return_dofs.append(hstr(zret))
    example_directions: list[dict[str, Any]] = []
    for zret in (Z0ret, Z1ret):
        if hstr(zret) not in untouched_return_dofs:
            continue
        vector = [sp.Integer(1) if dof == zret else sp.Integer(0) for dof in dofs]
        delta_t0 = compact(sum(sp.diff(t0, dof) * vec for dof, vec in zip(dofs, vector)))
        delta_t1 = compact(sum(sp.diff(t1, dof) * vec for dof, vec in zip(dofs, vector)))
        example_directions.append(
            {
                "dof": hstr(zret),
                "preserves_all_collected_constraints": all(
                    bool_zero(sum(sp.diff(expr, dof) * vec for dof, vec in zip(dofs, vector)))
                    for expr in constraint_exprs
                ),
                "delta_T0": hstr(delta_t0),
                "delta_T1": hstr(delta_t1),
                "moves_return": not (bool_zero(delta_t0) and bool_zero(delta_t1)),
            }
        )
    selector = selector_provenance(mutation)
    determined = not native_moves
    underdetermined = native_moves
    return {
        "genuine_generator_dofs": [hstr(dof) for dof in dofs],
        "named_constraints_collected": named_constraints(port, transfers),
        "pathA29_premise_citation": pathA29_premise_citation(),
        "linearized_constraint_expressions": [hstr(expr) for expr in constraint_exprs],
        "linearized_constraint_matrix": [[hstr(value) for value in row] for row in rows],
        "linearized_constraint_rank": rank0,
        "native_nullspace_dimension": nullity,
        "return_augmented_rank": return_aug_rank,
        "return_moving_nullity": return_moving_nullity,
        "native_null_moves_return": native_moves,
        "untouched_return_dofs_from_named_constraints": untouched_return_dofs,
        "example_return_moving_directions": example_directions,
        "T0_pinned_by_collected_constraints": return_moving_nullity == 0
        or matrix_rank(rows + [grad_t0]) == rank0,
        "T1_pinned_by_collected_constraints": return_moving_nullity == 0
        or matrix_rank(rows + [grad_t1]) == rank0,
        "branch_selector": selector,
        "selector_constraint_rows": [
            {
                **row,
                "row": [hstr(value) for value in row["row"]],
            }
            for row in selector_constraint_rows(mutation)
        ],
        "injected_null_probe": {
            "present": mutation.inject_null,
            "dof": hstr(eta_null) if mutation.inject_null else None,
            "preserves_P0": eta_null not in port["P0_raw"].free_symbols if mutation.inject_null else False,
            "moves_T0_T1": native_moves if mutation.inject_null else False,
            "distinct_from_decoupled_knob": mutation.inject_null
            and eta_null not in {gain0, gain1},
        },
        "determined_prediction": determined,
        "underdetermined_not_predictive": underdetermined,
        "why_Z_survives": "No collected named Gate-5 linear constraint fixes the return admittance dofs; pathA_29 records Z_is_premise:true and boundary_dof:none."
        if underdetermined
        else "Selector equations touch the return admittance dofs and the recomputed return-moving nullity is zero.",
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
        if any(dim != first for dim in dims):
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


def dim_vector_text(dim: Dim) -> list[str]:
    return [hstr(v) for v in dim]


SOURCED_DIMS: dict[sp.Symbol, Dim] = {
    a: (1, 0, 0),
    c_s: (1, 0, -1),
    omega: (0, 0, -1),
    M0: (0, 1, -1),
    D1: (1, 1, -1),
    R0: (0, 1, -1),
    R1: (1, 1, -1),
    A0_sym: (0, 1, -1),
    A1_sym: (1, 1, -1),
    N0: (-1, 1, 0),
    D0: (-1, 1, -2),
    K0c: (0, 1, -2),
    Keta: (0, 1, -2),
    TOmega: (0, 1, -2),
    Z0ret: (0, 1, -2),
    Z1ret: (0, 1, -2),
    OmegaU: (0, 0, -1),
    OmegaW: (0, 0, -1),
    Rmix: (0, 0, -2),
    gU: (sp.Rational(-1, 2), sp.Rational(1, 2), -2),
    gW: (sp.Rational(-1, 2), sp.Rational(1, 2), -2),
    eta_null: ZERO_DIM,
    gain0: ZERO_DIM,
    gain1: ZERO_DIM,
    q_free: ZERO_DIM,
}

EXPECTED_DIMS = {
    "M0": SOURCED_DIMS[M0],
    "D1": SOURCED_DIMS[D1],
    "R0": SOURCED_DIMS[R0],
    "R1": SOURCED_DIMS[R1],
    "A0": SOURCED_DIMS[A0_sym],
    "A1": SOURCED_DIMS[A1_sym],
    "T0": ZERO_DIM,
    "T1": ZERO_DIM,
    "epsilon0_eff": ZERO_DIM,
    "epsilon1_eff": ZERO_DIM,
    "N0": SOURCED_DIMS[N0],
    "D0": SOURCED_DIMS[D0],
    "P0_physical": ZERO_DIM,
}


def dim_record(name: str, expr: sp.Expr, expected: Dim, symbol_dims: dict[sp.Symbol, Dim]) -> dict[str, Any]:
    dim = expected if sp.sympify(expr) == 0 else dim_of(expr, symbol_dims)
    return {
        "quantity": name,
        "expression": hstr(expr),
        "computed_dimension": dim_to_text(dim),
        "expected_sourced_dimension": dim_to_text(expected),
        "dimension_monomial": hstr(dim_to_monomial(dim)),
        "dimension_vector_LMT": dim_vector_text(dim),
        "matches_expected": dim == expected,
    }


def run_dimension_check(
    transfers: dict[str, Any],
    residuals: dict[str, Any],
    port: dict[str, Any],
    *,
    corrupt_sourced: bool = False,
    corrupt_free_carrier: bool = False,
) -> dict[str, Any]:
    dims = dict(SOURCED_DIMS)
    if corrupt_sourced:
        dims[M0] = dim_add(dims[M0], (1, 0, 0))
    if corrupt_free_carrier:
        dims[q_free] = (7, 0, 0)
    p0_physical = port["P0_physical"]
    entries = [
        ("M0", M0, EXPECTED_DIMS["M0"]),
        ("D1", D1, EXPECTED_DIMS["D1"]),
        ("R0", R0, EXPECTED_DIMS["R0"]),
        ("R1", R1, EXPECTED_DIMS["R1"]),
        ("A0", residuals["leading"]["A0"], EXPECTED_DIMS["A0"]),
        ("A1", residuals["leading"]["A1"], EXPECTED_DIMS["A1"]),
        ("T0", transfers["ell0"]["T_dc"], EXPECTED_DIMS["T0"]),
        ("T1", transfers["ell1"]["T_dc"], EXPECTED_DIMS["T1"]),
        ("epsilon0_eff", transfers["ell0"]["epsilon_eff"], EXPECTED_DIMS["epsilon0_eff"]),
        ("epsilon1_eff", transfers["ell1"]["epsilon_eff"], EXPECTED_DIMS["epsilon1_eff"]),
        ("N0_port", port["N0_from_port"], EXPECTED_DIMS["N0"]),
        ("D0", D0, EXPECTED_DIMS["D0"]),
        ("P0_physical", p0_physical, EXPECTED_DIMS["P0_physical"]),
        ("K0_collective", transfers["ell0"]["K_dc"], (0, 1, -2)),
        ("K1_harmonic", transfers["ell1"]["K_dc"], (0, 1, -2)),
    ]
    table: list[dict[str, Any]] = []
    errors: list[str] = []
    for name, expr, expected in entries:
        try:
            table.append(dim_record(name, expr, expected, dims))
        except DimError as exc:
            errors.append(f"{name}: {exc}")
    ok = not errors and all(row["matches_expected"] for row in table)
    return {
        "dimension_order": ["L", "M", "T"],
        "source_of_truth": {
            "M0": "moment definition int S_leak d^3x",
            "D1": "moment definition int x_i S_leak d^3x + int j_i d^3x",
            "R0": "return zeroth moment, same sourced dimension as M0",
            "R1": "return first moment, same sourced dimension as D1",
            "A0": "amplitude ledger slot independent of residual formula, same sourced dimension as M0",
            "A1": "amplitude ledger slot independent of residual formula, same sourced dimension as D1",
            "T_l": "transfer ratio, dimensionless by continuity partition",
            "epsilon_l": "return admittance ratio, dimensionless",
            "N0": "Gate-4 port moment sourced dimension L^-1 M",
            "D0": "Gate-4 conservative moment sourced dimension L^-1 T^-2 M",
        },
        "symbol_dimensions": {hstr(sym): dim_to_text(dim) for sym, dim in dims.items()},
        "table": table,
        "errors": errors,
        "dimensional_ok": ok,
        "verdict": "NO_FAIL" if ok else "FAIL_DIMENSIONAL",
        "corrupt_sourced": corrupt_sourced,
        "corrupt_free_carrier": corrupt_free_carrier,
        "free_carrier_independence_expression_mentions_q_free": any(
            str(q_free) in str(row["expression"]) for row in table
        ),
    }


def build_provenance(mutation: Mutation, rank: dict[str, Any]) -> dict[str, Any]:
    items: dict[str, dict[str, Any]] = {
        "Y_l_out_fingerprints": {
            "tags": ["dtn_hankel_expansion", "outgoing_boundary_condition"],
            "computed_class": "derived_in_gate",
        },
        "ell0_T0_map": {
            "tags": ["section_8_2_collective_delta_a_delta_L_reduction", "continuity_partition_solve"],
            "computed_class": "derived_in_gate",
        },
        "ell1_T1_map": {
            "tags": ["section_9_4_harmonic_projection", "continuity_partition_solve"],
            "computed_class": "derived_in_gate",
        },
        "ell2_P0_map": {
            "tags": ["section_10_2_10_3_port_kernel", "gate4_prefactor"],
            "computed_class": "derived_in_gate",
        },
        "epsilon_eff_nonzero_value": {
            "tags": ["native_nullspace_rank_audit"],
            "computed_class": "deferred_branch_data"
            if rank["underdetermined_not_predictive"]
            else "derived_in_gate",
            "magnitude_note": "not computed because native nullspace leaves Z0_ret/Z1_ret free"
            if rank["underdetermined_not_predictive"]
            else "computed from selected generator branch; literal pathA_29 magnitude not cross-checked",
        },
        "pathA_28_side_conditions": {
            "tags": ["external_R0_equals_minus_M0", "external_R1_equals_minus_D1"],
            "computed_class": "external_bridge_input",
        },
        "literal_pathA_29_epsilon_magnitudes": {
            "tags": ["pathA29_flat_slab_partition"],
            "computed_class": "deferred_branch_data",
        },
        "time_convention": {
            "tags": ["exp_minus_i_omega_t"],
            "computed_class": "convention",
        },
    }
    if mutation.assert_not_derive:
        items["Y_l_out_fingerprints"]["emitted_class"] = "asserted_literal"
    for item in items.values():
        item.setdefault("emitted_class", item["computed_class"])
        item["class_matches_computed"] = item["emitted_class"] == item["computed_class"]
    if rank["branch_selector"]["tautological_assertion"]:
        items["branch_selector"] = {
            "tags": ["selector_equation_without_named_pde_provenance"],
            "computed_class": "derived_in_gate",
            "emitted_class": "asserted_literal",
            "class_matches_computed": False,
        }
    ok = all(item["class_matches_computed"] for item in items.values())
    groups = {
        "derived_in_gate": [],
        "external_bridge_input": [],
        "deferred_branch_data": [],
        "convention": [],
        "asserted_literal": [],
    }
    for name, item in items.items():
        groups.setdefault(item["computed_class"], []).append(name)
    return {
        "items": items,
        "groups": groups,
        "ok": ok,
        "epsilon_magnitude_cases": {
            "computed_not_pathA29_cross_checked": not rank["underdetermined_not_predictive"],
            "not_computed_generator_underdetermined": rank["underdetermined_not_predictive"],
            "literal_pathA29_partition_magnitudes": "deferred_branch_data",
        },
    }


def detect_decoupling(port: dict[str, Any], transfers: dict[str, Any], mutation: Mutation) -> dict[str, Any]:
    introduced = [gain0, gain1]
    introduced_after_l2 = [
        hstr(sym)
        for sym in introduced
        if sym in transfers["ell0"]["T_dc"].free_symbols
        or sym in transfers["ell1"]["T_dc"].free_symbols
    ]
    p0_symbols = port["P0_raw"].free_symbols
    moves = {
        "gain0_moves_T0": not bool_zero(sp.diff(transfers["ell0"]["T_dc"], gain0)),
        "gain0_moves_T1": not bool_zero(sp.diff(transfers["ell1"]["T_dc"], gain0)),
        "gain1_moves_T0": not bool_zero(sp.diff(transfers["ell0"]["T_dc"], gain1)),
        "gain1_moves_T1": not bool_zero(sp.diff(transfers["ell1"]["T_dc"], gain1)),
    }
    p0_unaffected = all(sym not in p0_symbols for sym in introduced)
    independently_moves_return = (
        moves["gain0_moves_T0"]
        and not moves["gain0_moves_T1"]
        and moves["gain1_moves_T1"]
        and not moves["gain1_moves_T0"]
    )
    return {
        "introduced_after_l2_parameters": introduced_after_l2,
        "P0_unaffected_by_introduced_parameters": p0_unaffected
        and bool(introduced_after_l2),
        "independent_transfer_derivatives": moves,
        "epsilon_can_be_dialed_after_P0_fixed": p0_unaffected
        and independently_moves_return,
        "decoupled": bool(introduced_after_l2)
        and p0_unaffected
        and independently_moves_return,
    }


def base_verdict(conditions: dict[str, bool]) -> str:
    ordered_failures = [
        ("decoupled", "FAIL_DECOUPLED"),
        ("tautological", "FAIL_TAUTOLOGICAL"),
        ("quad_regression", "FAIL_QUAD_REGRESSION"),
        ("dimensional", "FAIL_DIMENSIONAL"),
        ("no_consistent_return", "FAIL_NO_CONSISTENT_RETURN"),
        ("overcancel", "FAIL_OVERCANCEL"),
        ("epsilon_mismatch", "FAIL_EPSILON_MISMATCH"),
        ("underdetermined", "FAIL_UNDERDETERMINED_NOT_PREDICTIVE"),
        ("able_to_fail_bad", "FAIL_TAUTOLOGICAL"),
    ]
    for gate, verdict in ordered_failures:
        if conditions.get(gate, False):
            return verdict
    return "CROSS_L_RESIDUAL_PREDICTION"


@lru_cache(maxsize=None)
def run_gate(mutation: Mutation = Mutation()) -> dict[str, Any]:
    fingerprints = build_fingerprints()
    port = build_port_kernel_for(mutation)
    gate4 = build_gate4_non_regression(fingerprints, port, mutation)
    transfers = build_transfers(mutation)
    residuals = build_residuals(fingerprints, transfers, mutation)
    rank = build_rank_audit(port, transfers, mutation)
    dim = run_dimension_check(
        transfers,
        residuals,
        port,
        corrupt_sourced=mutation.corrupt_dimension,
        corrupt_free_carrier=False,
    )
    dim_corrupt_sourced = run_dimension_check(transfers, residuals, port, corrupt_sourced=True)
    dim_corrupt_free = run_dimension_check(transfers, residuals, port, corrupt_free_carrier=True)
    provenance = build_provenance(mutation, rank)
    decoupling = detect_decoupling(port, transfers, mutation)
    conditions = {
        "decoupled": decoupling["decoupled"],
        "tautological": not provenance["ok"] or rank["branch_selector"]["tautological_assertion"],
        "quad_regression": not gate4["ok"],
        "dimensional": not dim["dimensional_ok"],
        "no_consistent_return": transfers["no_consistent_return"],
        "overcancel": transfers["overcancel"],
        "epsilon_mismatch": not residuals["pathA_29_comparison"]["positive_bounded"]
        or not residuals["pathA_29_comparison"]["A0_form"]
        or not residuals["pathA_29_comparison"]["A1_form"]
        or not residuals["pathA_29_comparison"]["A0_order"]
        or not residuals["pathA_29_comparison"]["A1_order"],
        "underdetermined": rank["underdetermined_not_predictive"],
        "able_to_fail_bad": False,
    }
    verdict = base_verdict(conditions)
    return {
        "mutation": mutation.__dict__,
        "fingerprints": fingerprints,
        "gate4_non_regression": gate4,
        "port_kernel": port,
        "transfers": transfers,
        "residuals": residuals,
        "rank_audit": rank,
        "dimensional_check": dim,
        "dimensional_mutations": {
            "corrupt_sourced_input": dim_corrupt_sourced,
            "corrupt_unconstrained_carrier": dim_corrupt_free,
        },
        "provenance_partition": provenance,
        "decoupling_audit": decoupling,
        "conditions": conditions,
        "verdict": verdict,
        "which_rung": verdict,
    }


def ablation(
    baseline: Mutation,
    mutated: Mutation,
    expected_fail: str,
) -> dict[str, Any]:
    with_ctx = run_gate(mutated)
    without_ctx = run_gate(baseline)
    with_mutation = with_ctx["verdict"]
    without_mutation = without_ctx["verdict"]
    return {
        "rerun_gate_logic": True,
        "with_mutation": with_mutation,
        "without_mutation": without_mutation,
        "expected_fail": expected_fail,
        "fail_suppressed": with_mutation == expected_fail and without_mutation != expected_fail,
        "with_conditions": with_ctx["conditions"],
        "without_conditions": without_ctx["conditions"],
    }


def build_counterfactuals(actual_baseline: Mutation) -> dict[str, Any]:
    clean = replace(
        actual_baseline,
        name="clean_selector_control",
        selector_equation_set="derived_pde_admissibility",
    )
    actual_ctx = run_gate(actual_baseline)
    clean_ctx = run_gate(clean)
    probes: dict[str, Any] = {}
    port_corrupt = run_gate(replace(actual_baseline, name="probe_R1_corrupt_port", corrupt_port_kernel=True))
    probes["R1_port_kernel_dependency"] = {
        "description": "Corrupt Omega_U in the ell=2 port kernel and recompute P0_raw plus the rank row.",
        "baseline_P0_raw": hstr(actual_ctx["port_kernel"]["P0_raw"]),
        "corrupt_P0_raw": hstr(port_corrupt["port_kernel"]["P0_raw"]),
        "baseline_rank_row": actual_ctx["rank_audit"]["linearized_constraint_matrix"][0],
        "corrupt_rank_row": port_corrupt["rank_audit"]["linearized_constraint_matrix"][0],
        "P0_changes": not bool_zero(
            actual_ctx["port_kernel"]["P0_raw"] - port_corrupt["port_kernel"]["P0_raw"]
        ),
        "ell2_determinacy_row_changes": actual_ctx["rank_audit"]["linearized_constraint_matrix"][0]
        != port_corrupt["rank_audit"]["linearized_constraint_matrix"][0],
    }
    probes["3a_decouple_knobs"] = {
        "description": "Inject gain0/gain1 after ell=2 branch is fixed.",
        "self_ablation": ablation(
            clean,
            replace(clean, name="probe_3a_decouple", decouple_knobs=True),
            "FAIL_DECOUPLED",
        ),
    }
    probes["3b_null_direction"] = {
        "real_selector_equation_able_to_fail": {
            "rerun_gate_logic": True,
            "selector_equations": clean_ctx["rank_audit"]["branch_selector"]["equations"],
            "with_selector": clean_ctx["verdict"],
            "without_selector": actual_ctx["verdict"],
            "with_selector_return_moving_nullity": clean_ctx["rank_audit"][
                "return_moving_nullity"
            ],
            "without_selector_return_moving_nullity": actual_ctx["rank_audit"][
                "return_moving_nullity"
            ],
            "verdict_flips_to_prediction": clean_ctx["verdict"]
            == "CROSS_L_RESIDUAL_PREDICTION"
            and actual_ctx["verdict"] == "FAIL_UNDERDETERMINED_NOT_PREDICTIVE",
        },
        "detector_self_ablation_clean_native": ablation(
            clean,
            replace(clean, name="probe_3b_injected_null", inject_null=True),
            "FAIL_UNDERDETERMINED_NOT_PREDICTIVE",
        ),
        "baseline_native_no_injection": {
            "rerun_gate_logic": True,
            "verdict": actual_ctx["verdict"],
            "native_null_moves_return": actual_ctx["rank_audit"]["native_null_moves_return"],
            "selector_present": actual_ctx["rank_audit"]["branch_selector"]["present"],
            "untouched_return_dofs": actual_ctx["rank_audit"][
                "untouched_return_dofs_from_named_constraints"
            ],
        },
    }
    probes["3c_wrong_sign_antilocalizing"] = {
        "description": "Flip return admittance sign so eps_eff is negative/unbounded.",
        "self_ablation": ablation(
            actual_baseline,
            replace(actual_baseline, name="probe_3c_wrong_sign", wrong_sign_return=True),
            "FAIL_EPSILON_MISMATCH",
        ),
    }
    probes["3d_perfect_return"] = {
        "description": "Force Z0_ret=Z1_ret=0 so T0=T1=1 and eps_eff=0.",
        "self_ablation": ablation(
            actual_baseline,
            replace(actual_baseline, name="probe_3d_perfect", perfect_return=True),
            "FAIL_OVERCANCEL",
        ),
    }
    probes["3e_break_gate4"] = {
        "description": "Use incoming ell=2 and N/D prefactor in place of outgoing D0*N/Dcons^2.",
        "self_ablation": ablation(
            actual_baseline,
            replace(actual_baseline, name="probe_3e_break_gate4", break_gate4=True),
            "FAIL_QUAD_REGRESSION",
        ),
    }
    dim_with = run_gate(replace(actual_baseline, name="probe_3f_corrupt_dim", corrupt_dimension=True))
    dim_without = run_gate(actual_baseline)
    probes["3f_corrupt_dimension"] = {
        "description": "Corrupt sourced [M0] by L; separately corrupt q_free, which is not in the verdict expressions.",
        "sourced_corruption_verdict": dim_with["verdict"],
        "unconstrained_carrier_verdict": dim_without["dimensional_mutations"][
            "corrupt_unconstrained_carrier"
        ]["verdict"],
        "self_ablation": {
            "rerun_gate_logic": True,
            "with_mutation": dim_with["verdict"],
            "without_mutation": dim_without["verdict"],
            "expected_fail": "FAIL_DIMENSIONAL",
            "fail_suppressed": dim_with["verdict"] == "FAIL_DIMENSIONAL"
            and dim_without["verdict"] != "FAIL_DIMENSIONAL",
        },
    }
    probes["3g_assert_not_derive"] = {
        "description": "Supply the selector equations as asserted literals with no named PDE provenance.",
        "self_ablation": ablation(
            actual_baseline,
            replace(
                actual_baseline,
                name="probe_3g_assert",
                selector_equation_set="asserted_unproven",
            ),
            "FAIL_TAUTOLOGICAL",
        ),
    }
    probes["3h_no_consistent_return"] = {
        "description": "Restrict branch so P0 match forces eps_eff=-2 and T is not bounded positive.",
        "self_ablation": ablation(
            actual_baseline,
            replace(actual_baseline, name="probe_3h_no_consistent", no_consistent_return=True),
            "FAIL_NO_CONSISTENT_RETURN",
        ),
    }

    flags: dict[str, bool] = {}
    for name, probe in probes.items():
        if name == "R1_port_kernel_dependency":
            flags[name] = probe["P0_changes"] and probe["ell2_determinacy_row_changes"]
        elif name == "3b_null_direction":
            flags[name] = probe["detector_self_ablation_clean_native"]["fail_suppressed"]
            flags[name] = flags[name] and probe["real_selector_equation_able_to_fail"][
                "verdict_flips_to_prediction"
            ]
        elif name == "3f_corrupt_dimension":
            flags[name] = probe["self_ablation"]["fail_suppressed"] and probe[
                "unconstrained_carrier_verdict"
            ] == "NO_FAIL"
        else:
            flags[name] = probe["self_ablation"]["fail_suppressed"]
    return {
        "probes": probes,
        "probe_flags": flags,
        "able_to_fail_ok": all(flags.values()),
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


def engine_summary(ctx: dict[str, Any]) -> dict[str, Any]:
    fp = ctx["fingerprints"]
    residuals = ctx["residuals"]
    transfers = ctx["transfers"]
    rank = ctx["rank_audit"]
    gate4 = ctx["gate4_non_regression"]
    dim = ctx["dimensional_check"]
    sample_subs = {
        a: 3.0,
        c_s: 2.0,
        omega: 0.07,
        M0: 5.0,
        D1: 11.0,
        K0c: 13.0,
        Keta: 17.0,
        TOmega: 19.0,
        Z0ret: 23.0,
        Z1ret: 29.0,
        N0: 31.0,
        D0: 37.0,
        N2: 41.0,
        D2: 43.0,
        N4: 47.0,
        D4: 53.0,
        OmegaU: 2.0,
        OmegaW: 3.0,
        Rmix: 1.0,
        gU: 5.0,
        gW: 7.0,
        eta_null: 0.25,
        gain0: 1.3,
        gain1: 0.7,
    }
    key_probes = build_engine_probe_summary()
    return {
        "schema": "pathA_34_sympy_scratch/v1",
        "engine": "sympy",
        "fingerprints": {
            f"ell{ell}": {
                "Y_series_z": hstr(fp["outgoing"][ell]["Y_series"]),
                "normalization_factor": fp["outgoing"][ell]["normalization_factor"],
                "radiative_coeff_z": hstr(fp["outgoing"][ell]["radiative_coeff_z"]),
                "raw_order": fp["outgoing"][ell]["raw_outgoing_order"],
                "incoming_radiative_coeff_z": hstr(fp["incoming"][ell]["radiative_coeff_z"]),
            }
            for ell in (0, 1, 2)
        },
        "residuals": {
            "T0_dc": hstr(transfers["ell0"]["T_dc"]),
            "T1_dc": hstr(transfers["ell1"]["T_dc"]),
            "epsilon0_eff": hstr(transfers["ell0"]["epsilon_eff"]),
            "epsilon1_eff": hstr(transfers["ell1"]["epsilon_eff"]),
            "A0_leading": hstr(residuals["leading"]["A0"]),
            "A1_leading": hstr(residuals["leading"]["A1"]),
            "A0_expected": hstr(residuals["leading"]["expected_A0"]),
            "A1_expected": hstr(residuals["leading"]["expected_A1"]),
            "A0_residual_to_expected": hstr(residuals["leading"]["A0_residual_to_pathA29_form"]),
            "A1_residual_to_expected": hstr(residuals["leading"]["A1_residual_to_pathA29_form"]),
            "A0_sample": numeric_complex_parts(residuals["leading"]["A0"], sample_subs),
            "A1_sample": numeric_complex_parts(residuals["leading"]["A1"], sample_subs),
            "T0_sample": numeric(transfers["ell0"]["T_dc"], sample_subs),
            "T1_sample": numeric(transfers["ell1"]["T_dc"], sample_subs),
            "epsilon0_sample": numeric(transfers["ell0"]["epsilon_eff"], sample_subs),
            "epsilon1_sample": numeric(transfers["ell1"]["epsilon_eff"], sample_subs),
        },
        "gate4_non_regression": {
            "ok": gate4["ok"],
            "chi_Q": hstr(gate4["chi_Q"]),
            "P0_residual": hstr(gate4["residuals"]["P0"]),
            "P2_residual": hstr(gate4["residuals"]["P2"]),
            "P4_residual": hstr(gate4["residuals"]["P4"]),
            "P0_sample": numeric(gate4["coefficients"]["P0"], sample_subs),
            "P2_sample": numeric(gate4["coefficients"]["P2"], sample_subs),
            "P4_sample": numeric(gate4["coefficients"]["P4"], sample_subs),
        },
        "rank": {
            "native_nullspace_dimension": rank["native_nullspace_dimension"],
            "native_null_moves_return": rank["native_null_moves_return"],
            "return_moving_nullity": rank["return_moving_nullity"],
            "determined_prediction": rank["determined_prediction"],
            "baseline_verdict": ctx["verdict"],
            "selector_control_verdict": run_gate(
                Mutation(name="engine_selector_control", selector_equation_set="derived_pde_admissibility")
            )["verdict"],
        },
        "dimension": {
            "dimensional_ok": dim["dimensional_ok"],
            "corrupt_sourced_verdict": ctx["dimensional_mutations"]["corrupt_sourced_input"][
                "verdict"
            ],
            "corrupt_free_carrier_verdict": ctx["dimensional_mutations"][
                "corrupt_unconstrained_carrier"
            ]["verdict"],
        },
        "headline_booleans": {
            "fingerprints_ok": fp["ok"],
            "residual_form_ok": residuals["pathA_29_comparison"]["A0_form"]
            and residuals["pathA_29_comparison"]["A1_form"],
            "gate4_ok": gate4["ok"],
            "dimension_ok": dim["dimensional_ok"],
        },
        "verdict_probes": key_probes,
    }


def build_engine_probe_summary() -> dict[str, Any]:
    actual = Mutation()
    return {
        "rank_baseline": run_gate(actual)["verdict"],
        "rank_selector_equation": run_gate(
            Mutation(name="engine_selector_control", selector_equation_set="derived_pde_admissibility")
        )["verdict"],
        "3c_wrong_sign": {
            "with_mutation": run_gate(replace(actual, name="engine_3c", wrong_sign_return=True))[
                "verdict"
            ],
            "without_mutation": run_gate(actual)["verdict"],
        },
        "3d_perfect_return": {
            "with_mutation": run_gate(replace(actual, name="engine_3d", perfect_return=True))[
                "verdict"
            ],
            "without_mutation": run_gate(actual)["verdict"],
        },
        "3h_no_consistent_return": {
            "with_mutation": run_gate(
                replace(actual, name="engine_3h", no_consistent_return=True)
            )["verdict"],
            "without_mutation": run_gate(actual)["verdict"],
        },
        "3f_corrupt_dimension": {
            "with_mutation": run_gate(replace(actual, name="engine_3f", corrupt_dimension=True))[
                "verdict"
            ],
            "without_mutation": run_gate(actual)["verdict"],
        },
    }


def sympify_engine_number(text: Any) -> sp.Expr:
    return sp.sympify(str(text).replace("^", "**"))


def compare_engines(sym: dict[str, Any], mma: dict[str, Any]) -> dict[str, Any]:
    symbolic_checks: dict[str, bool] = {}
    symbolic_deltas: list[float] = []
    for ell in (0, 1, 2):
        for key in ("radiative_coeff_z", "raw_order", "normalization_factor"):
            left = sym["fingerprints"][f"ell{ell}"][key]
            right = mma["fingerprints"][f"ell{ell}"][key]
            if key == "radiative_coeff_z":
                ok = bool_zero(sympify_engine_number(left) - sympify_engine_number(right))
            else:
                ok = left == right
            symbolic_checks[f"ell{ell}_{key}"] = ok
            symbolic_deltas.append(0.0 if ok else 1.0)
    for key in ("A0_residual_to_expected", "A1_residual_to_expected"):
        ok = sym["residuals"][key] == "0" and mma["residuals"][key] == "0"
        symbolic_checks[key] = ok
        symbolic_deltas.append(0.0 if ok else 1.0)
    for key in ("P0_residual", "P2_residual", "P4_residual", "chi_Q"):
        left = sym["gate4_non_regression"][key]
        right = mma["gate4_non_regression"][key]
        ok = bool_zero(sympify_engine_number(left) - sympify_engine_number(right))
        symbolic_checks[f"gate4_{key}"] = ok
        symbolic_deltas.append(0.0 if ok else 1.0)
    for key in (
        "native_nullspace_dimension",
        "native_null_moves_return",
        "return_moving_nullity",
        "determined_prediction",
        "baseline_verdict",
        "selector_control_verdict",
    ):
        ok = sym["rank"][key] == mma["rank"][key]
        symbolic_checks[f"rank_{key}"] = ok
        symbolic_deltas.append(0.0 if ok else 1.0)
    for probe_name in (
        "rank_baseline",
        "rank_selector_equation",
        "3c_wrong_sign",
        "3d_perfect_return",
        "3h_no_consistent_return",
        "3f_corrupt_dimension",
    ):
        left = sym["verdict_probes"][probe_name]
        right = mma["verdict_probes"][probe_name]
        ok = left == right
        symbolic_checks[f"probe_{probe_name}"] = ok
        symbolic_deltas.append(0.0 if ok else 1.0)
    for key in ("dimensional_ok", "corrupt_sourced_verdict", "corrupt_free_carrier_verdict"):
        ok = sym["dimension"][key] == mma["dimension"][key]
        symbolic_checks[f"dimension_{key}"] = ok
        symbolic_deltas.append(0.0 if ok else 1.0)

    numeric_pairs: dict[str, tuple[float, float]] = {}
    for key in ("A0_sample", "A1_sample"):
        numeric_pairs[f"residual_{key}_re"] = (
            float(sym["residuals"][key]["re"]),
            float(mma["residuals"][key]["re"]),
        )
        numeric_pairs[f"residual_{key}_im"] = (
            float(sym["residuals"][key]["im"]),
            float(mma["residuals"][key]["im"]),
        )
    for key in (
        "T0_sample",
        "T1_sample",
        "epsilon0_sample",
        "epsilon1_sample",
    ):
        numeric_pairs[f"residual_{key}"] = (float(sym["residuals"][key]), float(mma["residuals"][key]))
    for key in ("P0_sample", "P2_sample", "P4_sample"):
        numeric_pairs[f"gate4_{key}"] = (
            float(sym["gate4_non_regression"][key]),
            float(mma["gate4_non_regression"][key]),
        )
    numeric_deltas = {name: abs(left - right) for name, (left, right) in numeric_pairs.items()}
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
        "compared_quantities": list(symbolic_checks) + list(numeric_deltas),
    }


def build_final_payload(ctx: dict[str, Any], sym_engine: dict[str, Any], mma_engine: dict[str, Any]) -> dict[str, Any]:
    counterfactuals = build_counterfactuals(Mutation())
    verdict = ctx["verdict"]
    if not counterfactuals["able_to_fail_ok"] and verdict != "FAIL_TAUTOLOGICAL":
        verdict = "FAIL_TAUTOLOGICAL"
    engine_agreement = compare_engines(sym_engine, mma_engine)
    if engine_agreement["status"] != "pass":
        verdict = "FAIL_ENGINE_DISAGREE"
    payload = {
        "schema": "pathA_34_cross_l_unification/v1",
        "engine": "dual-engine",
        "source_directive": "software/stage1_solver/directives/pathA_34_cross_l_unification.md",
        "verdict": verdict,
        "which_rung": verdict,
        "gate_conditions": ctx["conditions"],
        "headline": {
            "sympy": sym_engine,
            "mathematica": mma_engine,
        },
        "shared_generator_certificate": {
            "parameter_set": [
                "mu_eta",
                "T_w",
                "K_eta",
                "T_Omega",
                "c_s",
                "K0c",
                "Z0_ret",
                "Z1_ret",
                "D0",
                "Omega_U",
                "Omega_W",
                "R_mix",
                "g_U",
                "g_W",
                "eta_null_probe",
            ],
            "sector_maps": {
                "ell0": stringify_expr_tree(ctx["transfers"]["ell0"]),
                "ell1": stringify_expr_tree(ctx["transfers"]["ell1"]),
                "ell2": stringify_expr_tree(ctx["port_kernel"]),
            },
            "rank_audit": stringify_expr_tree(ctx["rank_audit"]),
            "decoupling_audit": ctx["decoupling_audit"],
        },
        "dtn_fingerprints": stringify_expr_tree(ctx["fingerprints"]),
        "residuals_and_transfers": stringify_expr_tree(
            {
                "transfers": ctx["transfers"],
                "residuals": ctx["residuals"],
            }
        ),
        "pathA_29_form_sign_order_comparison": ctx["residuals"]["pathA_29_comparison"],
        "gate4_non_regression": stringify_expr_tree(ctx["gate4_non_regression"]),
        "dimensional_check": stringify_expr_tree(ctx["dimensional_check"]),
        "dimensional_homogeneity_table": stringify_expr_tree(ctx["dimensional_check"]["table"]),
        "dimensional_able_to_fail": {
            "corrupt_sourced_input": stringify_expr_tree(
                ctx["dimensional_mutations"]["corrupt_sourced_input"]
            ),
            "corrupt_unconstrained_carrier": stringify_expr_tree(
                ctx["dimensional_mutations"]["corrupt_unconstrained_carrier"]
            ),
        },
        "provenance_partition": ctx["provenance_partition"],
        "counterfactual_probes": stringify_expr_tree(counterfactuals["probes"]),
        "able_to_fail": {
            "probe_flags": counterfactuals["probe_flags"],
            "able_to_fail_ok": counterfactuals["able_to_fail_ok"],
        },
        "engine_agreement": engine_agreement,
        "honest_gap": {
            "status": "native_nullspace_leaves_scalar_dipole_return_unfixed"
            if verdict == "FAIL_UNDERDETERMINED_NOT_PREDICTIVE"
            else "see_verdict",
            "gate6_selector_requirement": "A nonlinear branch selector must fix Z0_ret and Z1_ret, or an equivalent return law, before ell=0/1 become predictions.",
            "pathA29_premise_citation": ctx["rank_audit"]["pathA29_premise_citation"],
        },
    }
    return stringify_expr_tree(payload)


def build_report(payload: dict[str, Any]) -> str:
    cert = payload["shared_generator_certificate"]
    eng = payload["engine_agreement"]
    probes = payload["counterfactual_probes"]
    rank = cert["rank_audit"]
    dim = payload["dimensional_check"]
    lines = [
        "# PathA-34 cross-ell unification result",
        "",
        f"Computed verdict: `{payload['verdict']}`.",
        "",
        "The dual-engine headline derives the outgoing DtN fingerprints for ell=0, ell=1, and ell=2, "
        "assembles the scalar/dipole residuals, and reruns the Gate-4 quadrupole non-regression check. "
        "The baseline does not pass: the native fixed-P0 nullspace moves `T0(0)`/`T1(0)` and no derived "
        "Gate-5 selector removes that freedom.",
        "",
        "## Sector provenance",
        "",
        f"- ell=0 map: `{cert['sector_maps']['ell0']['route']}`.",
        f"- ell=1 map: `{cert['sector_maps']['ell1']['route']}`.",
        f"- ell=2 map: `{cert['sector_maps']['ell2']['route']}`.",
        f"- Native nullspace dimension: `{rank['native_nullspace_dimension']}`; "
        f"return-moving nullity: `{rank['return_moving_nullity']}`; moves return: "
        f"`{rank['native_null_moves_return']}`; selector present: `{rank['branch_selector']['present']}`.",
        f"- Named constraints collected: "
        f"`{', '.join(item['name'] for item in rank['named_constraints_collected'])}`.",
        f"- Return dofs untouched by named constraints: "
        f"`{rank['untouched_return_dofs_from_named_constraints']}`; pathA_29 premise: "
        f"`Z_is_premise={rank['pathA29_premise_citation']['Z_is_premise']}, "
        f"boundary_dof={rank['pathA29_premise_citation']['boundary_dof']}`.",
        f"- Selector equation control: "
        f"`{probes['3b_null_direction']['real_selector_equation_able_to_fail']['without_selector']}` -> "
        f"`{probes['3b_null_direction']['real_selector_equation_able_to_fail']['with_selector']}` using "
        f"`{probes['3b_null_direction']['real_selector_equation_able_to_fail']['selector_equations']}`.",
        "",
        "## Residual class",
        "",
        f"- `A0`: `{payload['residuals_and_transfers']['residuals']['leading']['A0']}`.",
        f"- `A1`: `{payload['residuals_and_transfers']['residuals']['leading']['A1']}`.",
        f"- `epsilon0_eff`: `{payload['residuals_and_transfers']['transfers']['ell0']['epsilon_eff']}`; "
        f"`epsilon1_eff`: `{payload['residuals_and_transfers']['transfers']['ell1']['epsilon_eff']}`.",
        "",
        "## Dimensional check",
        "",
        f"- Dimensional verdict: `{dim['verdict']}`; `dimensional_ok={dim['dimensional_ok']}`.",
        f"- Sourced corruption verdict: "
        f"`{payload['dimensional_able_to_fail']['corrupt_sourced_input']['verdict']}`.",
        f"- Unconstrained carrier corruption verdict: "
        f"`{payload['dimensional_able_to_fail']['corrupt_unconstrained_carrier']['verdict']}`.",
        "",
        "## Dual-engine agreement",
        "",
        f"- Status: `{eng['status']}`; max symbolic delta `{eng['max_symbolic_delta']}`; "
        f"max numeric delta `{eng['max_numeric_delta']}`.",
        "",
        "## Probe verdicts",
        "",
        f"- R1 port-kernel dependency: P0 changes "
        f"`{probes['R1_port_kernel_dependency']['P0_changes']}`; rank row changes "
        f"`{probes['R1_port_kernel_dependency']['ell2_determinacy_row_changes']}`.",
        f"- 3a decouple knobs: `{probes['3a_decouple_knobs']['self_ablation']['with_mutation']}` / "
        f"`{probes['3a_decouple_knobs']['self_ablation']['without_mutation']}`.",
        f"- 3b injected null detector: "
        f"`{probes['3b_null_direction']['detector_self_ablation_clean_native']['with_mutation']}` / "
        f"`{probes['3b_null_direction']['detector_self_ablation_clean_native']['without_mutation']}`; "
        f"baseline native: `{probes['3b_null_direction']['baseline_native_no_injection']['verdict']}`.",
        f"- 3c wrong sign: `{probes['3c_wrong_sign_antilocalizing']['self_ablation']['with_mutation']}` / "
        f"`{probes['3c_wrong_sign_antilocalizing']['self_ablation']['without_mutation']}`.",
        f"- 3d perfect return: `{probes['3d_perfect_return']['self_ablation']['with_mutation']}` / "
        f"`{probes['3d_perfect_return']['self_ablation']['without_mutation']}`.",
        f"- 3e break Gate 4: `{probes['3e_break_gate4']['self_ablation']['with_mutation']}` / "
        f"`{probes['3e_break_gate4']['self_ablation']['without_mutation']}`.",
        f"- 3f corrupt dimension: `{probes['3f_corrupt_dimension']['self_ablation']['with_mutation']}` / "
        f"`{probes['3f_corrupt_dimension']['self_ablation']['without_mutation']}`; free carrier "
        f"`{probes['3f_corrupt_dimension']['unconstrained_carrier_verdict']}`.",
        f"- 3g assert not derive: `{probes['3g_assert_not_derive']['self_ablation']['with_mutation']}` / "
        f"`{probes['3g_assert_not_derive']['self_ablation']['without_mutation']}`.",
        f"- 3h no consistent return: `{probes['3h_no_consistent_return']['self_ablation']['with_mutation']}` / "
        f"`{probes['3h_no_consistent_return']['self_ablation']['without_mutation']}`.",
        "",
        "## Earned vs deferred",
        "",
        "Earned: the `-(ell+1)/Lambda_l^out` fingerprints, raw radiative orders, residual form/sign/order "
        "conditional on a positive bounded branch, and the Gate-4 non-regression. Deferred: the scalar/dipole "
        "return magnitude and nonzero prediction, because the native nullspace leaves `epsilon_eff` free at "
        "the linear Gate-5 level.",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    ctx = run_gate(Mutation())
    sym_engine = engine_summary(ctx)
    yaml_write(SYM_YAML, stringify_expr_tree(sym_engine))
    mma_engine = yaml_read(MMA_YAML)
    if mma_engine is None:
        print(f"Wrote SymPy scratch: {SYM_YAML}")
        print(f"Pending Mathematica scratch: {MMA_YAML}")
        return 0
    payload = build_final_payload(ctx, sym_engine, mma_engine)
    yaml_write(RESULTS_YAML, payload)
    REPORT_MD.write_text(build_report(payload), encoding="utf-8")
    print(f"Wrote results: {RESULTS_YAML}")
    print(f"Wrote report: {REPORT_MD}")
    print(
        "verdict={verdict} max_numeric_delta={delta}".format(
            verdict=payload["verdict"], delta=payload["engine_agreement"]["max_numeric_delta"]
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
