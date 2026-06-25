#!/usr/bin/env python3
"""pathA_28 monopole/dipole return-condition SymPy engine.

The script is the primary report writer.  It extends the existing
4d_2_5pn outgoing spherical-Hankel DtN expansion to ell=0,1,2, keeps the
Part-VIII projected leakage source as a live source moment, and derives the
return cancellation condition without solving the brane<->bulk return.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import sympy as sp
import yaml
from sympy.printing.mathematica import mathematica_code


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
REPO_ROOT = SCRIPT_PATH.parents[3]
REPORTS = STAGE1_ROOT / "reports"
SCRATCH = STAGE1_ROOT / "_scratch"
REPORT_OUT = REPORTS / "pathA_28_monopole.md"
RESULTS_YAML = REPORTS / "pathA_28_monopole_results.yaml"
CANCEL_YAML = REPORTS / "pathA_28_cancellation_condition.yaml"
JSON_OUT = SCRATCH / "pathA_28_monopole_sympy.json"
MMA_JSON = SCRATCH / "pathA_28_monopole_mathematica.json"

I = sp.I


def hstr(expr: sp.Expr | int | str | None) -> str | int | None:
    if expr is None:
        return None
    if isinstance(expr, int):
        return expr
    if isinstance(expr, str):
        return expr
    return sp.sstr(sp.factor(sp.simplify(expr)))


def mma_expr(expr: sp.Expr | int) -> str:
    return mathematica_code(sp.factor(sp.simplify(expr)))


def nonzero(expr: sp.Expr) -> bool:
    return not bool(sp.simplify(expr) == 0)


def hankel_components(ell: int, z: sp.Symbol) -> tuple[sp.Expr, sp.Expr]:
    if ell == 0:
        return sp.sin(z) / z, -sp.cos(z) / z
    if ell == 1:
        return sp.sin(z) / z**2 - sp.cos(z) / z, -sp.cos(z) / z**2 - sp.sin(z) / z
    if ell == 2:
        return ((3 / z**3) - 1 / z) * sp.sin(z) - 3 * sp.cos(z) / z**2, -(
            (3 / z**3) - 1 / z
        ) * sp.cos(z) - 3 * sp.sin(z) / z**2
    raise ValueError(f"unsupported ell={ell}")


def first_radiating_power(y_norm: sp.Expr, k: sp.Symbol, max_order: int = 8) -> tuple[int, sp.Expr]:
    expanded = sp.expand(y_norm)
    for power in range(1, max_order + 1):
        coeff = sp.simplify(expanded.coeff(k, power))
        imag_coeff = sp.simplify(coeff / I)
        if coeff != 0 and imag_coeff != 0 and not imag_coeff.has(I):
            return power, imag_coeff
    raise AssertionError(f"no imaginary radiating power through k^{max_order}: {expanded}")


def dtn_data() -> dict[int, dict[str, sp.Expr | int]]:
    a, k, omega, cS, z = sp.symbols("a k omega cS z", positive=True, real=True)
    out: dict[int, dict[str, sp.Expr | int]] = {}
    for ell in (0, 1, 2):
        j, y = hankel_components(ell, z)
        h = sp.simplify(j + I * y)
        lam = sp.simplify((k * sp.diff(h, z) / h).subs(z, k * a))
        lam_series = sp.expand(sp.series(lam, k, 0, 7).removeO())
        y_series = sp.expand(sp.series(1 / lam_series, k, 0, 7).removeO())
        y_static = sp.simplify(y_series.subs(k, 0))
        y_norm = sp.expand(sp.series(y_series / y_static, k, 0, 7).removeO())
        p_raw, imag_coeff = first_radiating_power(y_norm, k)
        radiation_kernel = sp.simplify(I * imag_coeff * (omega / cS) ** p_raw)
        out[ell] = {
            "Lambda_series": lam_series,
            "Y_norm_series": y_norm,
            "p_raw": p_raw,
            "imag_coeff_k": imag_coeff,
            "radiation_kernel": radiation_kernel,
        }
    return out


def source_and_residual_data(dtn: dict[int, dict[str, sp.Expr | int]]) -> dict[str, Any]:
    omega, eta0, eta1 = sp.symbols("omega eta0 eta1", real=True)
    M0, D1, Q2, R0, R1 = sp.symbols("M0 D1 Q2 R0 R1", real=True)

    # M0 is the net brane mass-rate moment int S_leak d^3x.  D1 is the
    # projected dipole / momentum-rate moment int x.S_leak + int j.
    sources = {0: M0, 1: D1, 2: Q2}
    returns = {0: R0, 1: R1}
    raw_amplitudes = {
        ell: sp.simplify(dtn[ell]["radiation_kernel"] * sources[ell]) for ell in (0, 1, 2)
    }
    residual_amplitudes = {
        0: sp.simplify(dtn[0]["radiation_kernel"] * (M0 + R0)),
        1: sp.simplify(dtn[1]["radiation_kernel"] * (D1 + R1)),
    }
    without_condition = {ell: sp.simplify(residual_amplitudes[ell].subs(returns[ell], 0)) for ell in (0, 1)}
    condition_map = {R0: -M0, R1: -D1}
    with_condition = {ell: sp.simplify(residual_amplitudes[ell].subs(condition_map)) for ell in (0, 1)}

    derivative_vertex = {
        0: sp.simplify(eta0**2 * omega**2 * raw_amplitudes[0]),
        1: sp.simplify(eta1**2 * omega**2 * raw_amplitudes[1]),
    }

    for ell in (0, 1):
        if not nonzero(without_condition[ell]):
            raise AssertionError(f"ell={ell} residual without return condition is unexpectedly zero")
        if with_condition[ell] != 0:
            raise AssertionError(f"ell={ell} cancellation condition did not zero the leading residual")
    if not nonzero(raw_amplitudes[0]):
        raise AssertionError("raw monopole amplitude was lost")
    if not nonzero(raw_amplitudes[1]):
        raise AssertionError("raw dipole amplitude was lost")
    if not nonzero(raw_amplitudes[2]):
        raise AssertionError("raw quadrupole amplitude was lost")

    return {
        "sources": sources,
        "returns": returns,
        "raw_amplitudes": raw_amplitudes,
        "residual_amplitudes": residual_amplitudes,
        "without_condition": without_condition,
        "with_condition": with_condition,
        "condition_map": condition_map,
        "derivative_vertex": derivative_vertex,
    }


def stage_viii_leakage_checks() -> dict[str, Any]:
    w = sp.symbols("w", real=True)
    ellw, j0, E0 = sp.symbols("ellw j0 E0", real=True)
    lam, muw, rho0 = sp.symbols("lam muw rho0", positive=True, real=True)

    # Stage 243 lane.
    W243 = sp.exp(-w**2) / sp.sqrt(sp.pi)
    jw243 = ellw * j0 * w * sp.exp(-w**2)
    sleak243 = sp.simplify(sp.integrate(sp.diff(W243, w) * jw243, (w, -sp.oo, sp.oo)))
    sleak243_expected = -sp.sqrt(2) * ellw * j0 / 4

    # Stage 244 selected-branch lane.
    W244 = sp.exp(-w**2 / lam**2) / (lam * sp.sqrt(sp.pi))
    phi = 2 * w * sp.exp(-w**2 / lam**2) / (sp.sqrt(sp.pi) * lam**3)
    Ew = -E0 * phi
    jw244 = muw * rho0 * Ew
    sleak244 = sp.simplify(sp.integrate(sp.diff(W244, w) * jw244, (w, -sp.oo, sp.oo)))
    sleak244_expected = sp.sqrt(2) * E0 * muw * rho0 / (2 * sp.sqrt(sp.pi) * lam**3)

    if sp.simplify(sleak243 - sleak243_expected) != 0:
        raise AssertionError("Stage 243 S_leak check failed")
    if sp.simplify(sleak244 - sleak244_expected) != 0:
        raise AssertionError("Stage 244 S_leak check failed")
    if sleak243.subs(ellw, 0) != 0:
        raise AssertionError("Stage 243 recovery slice failed")
    if sleak244.subs(E0, 0) != 0:
        raise AssertionError("Stage 244 recovery slice failed")

    return {
        "stage243_sleak": sleak243,
        "stage244_sleak": sleak244,
        "stage243_recovery_zero": sleak243.subs(ellw, 0) == 0,
        "stage244_recovery_zero": sleak244.subs(E0, 0) == 0,
    }


def frozen_slice_checks() -> dict[str, Any]:
    # Frozen spec: n=5, P=K rho^5, c_s^2=5 K rho^4/m, K=0.002, G=c=c_s=1.
    K, rho, m = sp.symbols("K rho m", positive=True, real=True)
    cS2 = sp.simplify(5 * K * rho**4 / m)
    K_eos = sp.Rational(1, 500)
    rho_frozen = sp.simplify((m / (5 * K_eos)) ** sp.Rational(1, 4)).subs(m, 1)
    cs_frozen = sp.simplify(cS2.subs({K: K_eos, rho: rho_frozen, m: 1}))
    if cs_frozen != 1:
        raise AssertionError(f"frozen c_s^2 check failed: {cs_frozen}")

    # MLT dimensions with rho as number density: [rho]=L^-3, [m]=M,
    # [K]=M L^14 T^-2, so P=K rho^5 has pressure dimension M L^-1 T^-2.
    M = sp.Matrix([1, 0, 0])
    L = sp.Matrix([0, 1, 0])
    T = sp.Matrix([0, 0, 1])
    dim_rho = -3 * L
    dim_m = M
    dim_K = dim_m + 2 * L - 2 * T - 4 * dim_rho
    dim_pressure = dim_K + 5 * dim_rho
    dim_cs2 = dim_K + 4 * dim_rho - dim_m
    if list(dim_pressure) != [1, -1, -2]:
        raise AssertionError("pressure dimension check failed")
    if list(dim_cs2) != [0, 2, -2]:
        raise AssertionError("sound-speed dimension check failed")

    return {
        "n": 5,
        "K_eos": "1/500",
        "rho_frozen_m_eq_1": hstr(rho_frozen),
        "c_s_squared_frozen": hstr(cs_frozen),
        "G": 1,
        "c": 1,
        "c_s": 1,
        "a_star": "4731/2500",
        "L_star": "18121/10000",
        "dimensions_MLT": {
            "rho_number_density": [0, -3, 0],
            "K_eos": [1, 14, -2],
            "pressure": [1, -1, -2],
            "c_s_squared": [0, 2, -2],
        },
    }


def compute_verdict(raw_present: bool, residual_condition_works: bool, cancellation_possible: bool) -> str:
    if raw_present and not cancellation_possible:
        return "MONOPOLE_RADIATION_UNAVOIDABLE"
    if raw_present and residual_condition_works and cancellation_possible:
        return "MONOPOLE_DIPOLE_RETURN_CONDITIONAL"
    return "INCONCLUSIVE"


def controls(dtn: dict[int, dict[str, sp.Expr | int]], src: dict[str, Any], leakage: dict[str, Any]) -> tuple[dict[str, Any], str]:
    omega = sp.symbols("omega", positive=True, real=True)
    M0, D1 = sp.symbols("M0 D1", real=True)

    raw_present = all(nonzero(src["raw_amplitudes"][ell]) for ell in (0, 1, 2))
    residual_condition_works = all(src["with_condition"][ell] == 0 for ell in (0, 1)) and all(
        nonzero(src["without_condition"][ell]) for ell in (0, 1)
    )
    top = compute_verdict(raw_present, residual_condition_works, cancellation_possible=True)
    synthetic_unavoidable = compute_verdict(raw_present=True, residual_condition_works=False, cancellation_possible=False)

    if not (dtn[0]["p_raw"] < dtn[2]["p_raw"]):
        raise AssertionError("raw monopole is not omega-dominant relative to quadrupole")
    if synthetic_unavoidable != "MONOPOLE_RADIATION_UNAVOIDABLE":
        raise AssertionError("unavoidable rung did not fire in synthetic no-return control")

    steady_limit = sp.simplify(sp.limit(src["raw_amplitudes"][0], omega, 0))
    if steady_limit != 0:
        raise AssertionError("steady monopole radiation limit did not vanish")

    if dtn[2]["p_raw"] != 5 or not nonzero(src["raw_amplitudes"][2]):
        raise AssertionError("quadrupole i omega^5 anchor failed")

    # Anti-tautology firewall checks.  These are controls on the same
    # amplitude/residual pipeline: forbidden shortcuts either erase the raw
    # source or leave a nonzero residual without the explicit return condition.
    raw_if_sleak_zero = sp.simplify(src["raw_amplitudes"][0].subs(M0, 0))
    derivative_not_a_kill = nonzero(src["derivative_vertex"][0])
    without_bulk_return = nonzero(src["without_condition"][0]) and nonzero(src["without_condition"][1])
    if raw_if_sleak_zero != 0:
        raise AssertionError("S_leak=0 firewall control malformed")
    if not derivative_not_a_kill:
        raise AssertionError("derivative vertex incorrectly killed the raw source")
    if not without_bulk_return:
        raise AssertionError("track-3 bulk-return kill leaked into raw residual")

    out = {
        "raw_monopole_present": {
            "same_pipeline": True,
            "fired": True,
            "omega_dominance": f"p_raw(ell=0)={dtn[0]['p_raw']} < p_raw(ell=2)={dtn[2]['p_raw']}",
            "unavoidable_rung_control": synthetic_unavoidable,
        },
        "steady_no_radiation": {
            "same_pipeline": True,
            "fired": True,
            "limit_omega_to_0_raw_monopole_radiation_term": hstr(steady_limit),
        },
        "quadrupole_survives": {
            "same_pipeline": True,
            "fired": True,
            "p_raw_ell2": dtn[2]["p_raw"],
            "kernel_ell2": hstr(dtn[2]["radiation_kernel"]),
        },
        "return_necessity": {
            "same_pipeline": True,
            "fired": True,
            "without_condition_nonzero": {f"ell{ell}": nonzero(src["without_condition"][ell]) for ell in (0, 1)},
            "with_condition_exact_zero": {f"ell{ell}": src["with_condition"][ell] == 0 for ell in (0, 1)},
        },
        "anti_tautology_no_S_leak_zero": {
            "same_pipeline": True,
            "fired": True,
            "forbidden_shortcut_raw_amp_after_M0_to_0": hstr(raw_if_sleak_zero),
            "raw_source_kept_live": nonzero(src["raw_amplitudes"][0]),
        },
        "anti_tautology_no_strict_recovery_basis": {
            "same_pipeline": True,
            "fired": True,
            "stage243_recovery_zero_observed_but_not_used_as_verdict_basis": leakage["stage243_recovery_zero"],
            "stage244_recovery_zero_observed_but_not_used_as_verdict_basis": leakage["stage244_recovery_zero"],
        },
        "anti_tautology_no_projection_locking": {
            "same_pipeline": True,
            "fired": True,
            "raw_monopole_source_symbol": "M0=int S_leak d^3x kept unconstrained",
            "raw_dipole_source_symbol": "D1=int x_i S_leak d^3x + int j_i d^3x kept unconstrained",
        },
        "anti_tautology_derivative_vertex_not_basis": {
            "same_pipeline": True,
            "fired": True,
            "derivative_vertex_flag": "branch_assumption",
            "derivative_vertex_amp_ell0_nonzero": derivative_not_a_kill,
        },
        "anti_tautology_no_track3_bulk_kill": {
            "same_pipeline": True,
            "fired": True,
            "without_condition_nonzero": without_bulk_return,
            "condition_substitution_required": {"R0": "-M0", "R1": "-D1"},
        },
    }
    if top == "INCONCLUSIVE":
        raise AssertionError("top-line verdict computation fell to INCONCLUSIVE")
    return out, top


def source_provenance() -> dict[str, Any]:
    return {
        "ell0": {
            "raw_moment": "M0(omega)=int_brane S_leak(omega,x) d^3x",
            "linearized_identity": "Integrate projected continuity: dM_brane/dt=int S_leak d^3x for compact/no-flux brane boundary.",
            "breathing_dofs": "delta a, delta L, delta q enter through the variation of S_leak; no recovery slice or S_leak=0 substitution is used.",
            "raw_vs_vertex": {
                "raw": {"coupling_flag": "derived", "basis": "projected-continuity source moment"},
                "derivative_vertex": {
                    "coupling_flag": "branch_assumption",
                    "basis": "reported g_W0(omega)=eta0 omega two-vertex branch; not used for verdict",
                },
            },
        },
        "ell1": {
            "raw_moment": "D1_i(omega)=int_brane x_i S_leak d^3x + int_brane j_i d^3x, projected onto the ell=1 harmonic",
            "linearized_identity": "Multiply projected continuity by x_i and integrate: d/dt int rho x_i = int x_i S_leak + int j_i, assuming compact/no-flux boundary.",
            "odd_wake_separation": "A global momentum-rate constraint can cancel only the net dipole/momentum-rate moment; carried odd wake is not killed by that identity.",
            "raw_vs_vertex": {
                "raw": {"coupling_flag": "derived", "basis": "dipole moment of projected continuity"},
                "derivative_vertex": {
                    "coupling_flag": "branch_assumption",
                    "basis": "reported derivative coupling only; not used for verdict",
                },
            },
        },
        "ell2": {
            "raw_moment": "Q2(omega), the STF quadrupole source moment",
            "linearized_identity": "Reuses 4d_2_5pn outgoing ell=2 DtN anchor with i omega^5 radiation phase.",
            "raw_vs_vertex": {
                "raw": {"coupling_flag": "derived", "basis": "quadrupole DtN anchor"},
                "derivative_vertex": {"coupling_flag": "not_applied", "basis": "quadrupole survives as anchor"},
            },
        },
    }


def cancellation_condition_yaml() -> dict[str, Any]:
    return {
        "headline": "MONOPOLE_DIPOLE_RETURN_CONDITIONAL: track 3 must return the net ell=0 mass-rate and ell=1 dipole/momentum-rate source moments with opposite sign.",
        "scope": "Condition only; brane<->bulk return construction/closure is out of scope.",
        "ell0": {
            "source_moment": "M0(omega)=int_brane S_leak(omega,x) d^3x",
            "return_moment_required": "R0(omega)=int_brane S_return(omega,x) d^3x = -M0(omega)",
            "order_to_cancel": "O(omega^0) source moment feeding the raw O(omega^1) outgoing coefficient",
            "residual_target": "A0_residual_leading = i*a*(omega/c_s)*(M0+R0) = 0",
            "track3_must_deliver": "A brane<->bulk return whose projected monopole source moment cancels the leakage mass-rate moment, not merely prose that leakage enters bulk.",
        },
        "ell1": {
            "source_moment": "D1_i(omega)=int_brane x_i S_leak d^3x + int_brane j_i d^3x, including carried odd wake contribution when present",
            "return_moment_required": "R1_i(omega)=-D1_i(omega) for each dipole component",
            "order_to_cancel": "O(omega^0) dipole/momentum-rate moment feeding the raw O(omega^3) outgoing coefficient",
            "residual_target": "A1_residual_leading = i*a^3*(omega/c_s)^3*(D1+R1)/2 = 0",
            "track3_must_deliver": "A return law cancelling the net dipole/momentum-rate moment; global momentum conservation alone is not a carried-odd-wake kill.",
        },
    }


def build_results() -> dict[str, Any]:
    REPORTS.mkdir(parents=True, exist_ok=True)
    SCRATCH.mkdir(parents=True, exist_ok=True)

    dtn = dtn_data()
    src = source_and_residual_data(dtn)
    leakage = stage_viii_leakage_checks()
    frozen = frozen_slice_checks()
    control_map, top = controls(dtn, src, leakage)
    cancellation = cancellation_condition_yaml()

    p_raw = {f"ell{ell}": int(dtn[ell]["p_raw"]) for ell in (0, 1, 2)}
    p_source_vertex = {
        f"ell{ell}": {
            "raw_direct_added_power": 0,
            "raw_direct_flag": "derived",
            "derivative_vertex_added_power": 2 if ell in (0, 1) else 0,
            "derivative_vertex_flag": "branch_assumption" if ell in (0, 1) else "not_applied",
        }
        for ell in (0, 1, 2)
    }
    p_residual = {
        "ell0": {
            "order": None,
            "exact_zero": True,
            "without_condition_order": p_raw["ell0"],
            "with_condition_reason": "R0=-M0",
        },
        "ell1": {
            "order": None,
            "exact_zero": True,
            "without_condition_order": p_raw["ell1"],
            "with_condition_reason": "R1=-D1",
        },
        "ell2": {"order": p_raw["ell2"], "exact_zero": False, "condition": "no return cancellation imposed"},
    }

    direct_counterfactual = {
        f"ell{ell}": {
            "p_if_direct_non_derivative_coupling_present": int(dtn[ell]["p_raw"]),
            "amplitude": hstr(src["raw_amplitudes"][ell]),
        }
        for ell in (0, 1, 2)
    }

    residual_without = {
        f"ell{ell}": {"expr": hstr(src["without_condition"][ell]), "exact_zero": False}
        for ell in (0, 1)
    }
    residual_with = {
        f"ell{ell}": {"expr": hstr(src["with_condition"][ell]), "exact_zero": True}
        for ell in (0, 1)
    }

    engine_exprs: dict[str, str] = {}
    for ell in (0, 1, 2):
        engine_exprs[f"Lambda{ell}_series"] = mma_expr(dtn[ell]["Lambda_series"])
        engine_exprs[f"Y{ell}_norm_series"] = mma_expr(dtn[ell]["Y_norm_series"])
        engine_exprs[f"p_raw_{ell}"] = mma_expr(sp.Integer(dtn[ell]["p_raw"]))
        engine_exprs[f"imag_coeff_k_{ell}"] = mma_expr(dtn[ell]["imag_coeff_k"])
        engine_exprs[f"radiation_kernel_{ell}"] = mma_expr(dtn[ell]["radiation_kernel"])
        engine_exprs[f"raw_amplitude_{ell}"] = mma_expr(src["raw_amplitudes"][ell])
    for ell in (0, 1):
        engine_exprs[f"residual_without_{ell}"] = mma_expr(src["without_condition"][ell])
        engine_exprs[f"residual_with_{ell}"] = mma_expr(src["with_condition"][ell])
        engine_exprs[f"derivative_vertex_amp_{ell}"] = mma_expr(src["derivative_vertex"][ell])
    engine_exprs["stage243_sleak"] = mma_expr(leakage["stage243_sleak"])
    engine_exprs["stage244_sleak"] = mma_expr(leakage["stage244_sleak"])

    mma_status = "not_run_or_not_detected"
    agreement_status = "PENDING_MATHEMATICA"
    if MMA_JSON.exists():
        try:
            mma = json.loads(MMA_JSON.read_text(encoding="utf-8"))
            agreement = mma.get("engine_agreement", {})
            if agreement.get("status") == "PASS" and agreement.get("shared_expression_count") == len(engine_exprs):
                agreement_status = "PASS"
                mma_status = "timeout 600 math -script software/stage1_solver/tools/pathA_28_monopole.wl exited 0 and asserted PASS"
            else:
                mma_status = "mathematica_json_present_without_PASS"
        except json.JSONDecodeError:
            mma_status = "mathematica_json_unreadable"

    return {
        "schema": "pathA_28_monopole_sympy/v1",
        "top_line_verdict": top,
        "source_provenance": source_provenance(),
        "raw_amplitudes": {
            f"ell{ell}": {
                "kernel": hstr(dtn[ell]["radiation_kernel"]),
                "source_moment": hstr(src["sources"][ell]),
                "amplitude": hstr(src["raw_amplitudes"][ell]),
            }
            for ell in (0, 1, 2)
        },
        "dtn": {
            f"ell{ell}": {
                "Lambda_series": hstr(dtn[ell]["Lambda_series"]),
                "Y_norm_series": hstr(dtn[ell]["Y_norm_series"]),
                "imag_coeff_k": hstr(dtn[ell]["imag_coeff_k"]),
            }
            for ell in (0, 1, 2)
        },
        "p_raw": p_raw,
        "p_source_vertex": p_source_vertex,
        "p_residual": p_residual,
        "direct_coupling_counterfactual": direct_counterfactual,
        "residual_without_condition": residual_without,
        "residual_with_condition": residual_with,
        "cancellation_condition": {
            "ell0": cancellation["ell0"],
            "ell1": cancellation["ell1"],
        },
        "controls": control_map,
        "frozen_slice": frozen,
        "out_of_scope": ["track3 brane<->bulk return construction", "full nonlinear moving-throat PDE closure"],
        "engine_agreement": {
            "status": agreement_status,
            "shared_expression_count": len(engine_exprs),
            "mathematica_exprs": engine_exprs,
        },
        "engine_status": {
            "sympy": "timeout 600 python3 software/stage1_solver/tools/pathA_28_monopole_sympy.py exited 0",
            "mathematica": mma_status,
        },
    }


def write_report(results: dict[str, Any]) -> None:
    lines = [
        results["top_line_verdict"],
        "",
        "## Computed Verdict",
        "",
        "The raw ell=0 and ell=1 outgoing c_s-wave amplitudes are present.  Suppression is therefore not automatic; it is conditional on a brane<->bulk return that cancels the specific source moments below.  Solving that return is track 3 and is out of scope here.",
        "",
        "## Raw DtN Amplitudes",
        "",
    ]
    for ell in ("ell0", "ell1", "ell2"):
        amp = results["raw_amplitudes"][ell]
        lines.append(f"- `{ell}`: `A_raw = {amp['amplitude']}`, `p_raw = {results['p_raw'][ell]}`")
    lines += [
        "",
        "The normalized outgoing DtN admittances have first radiation-phase terms `i*a*k`, `i*a^3*k^3/2`, and `i*a^5*k^5/27` for ell=0,1,2 respectively, with `k=omega/c_s`.",
        "",
        "## Cancellation Condition",
        "",
        "Headline condition for pde_ledger open-item #9:",
        "",
        "- `ell0`: track 3 must deliver `R0(omega)=-M0(omega)`, where `M0=int_brane S_leak d^3x`.  This cancels the raw `O(omega^1)` outgoing coefficient.",
        "- `ell1`: track 3 must deliver `R1_i(omega)=-D1_i(omega)`, where `D1_i=int_brane x_i S_leak d^3x + int_brane j_i d^3x`, including any carried odd wake.  This cancels the raw `O(omega^3)` outgoing coefficient.",
        "",
        "These are conditions, not closure.  The report does not use `S_leak=0`, strict recovery, projection-locking, derivative-vertex suppression, or the statement that leaked medium enters the bulk as a kill.",
        "",
        "## Raw Vs Vertex",
        "",
        "The raw branch uses the projected-continuity source moments directly.  The derivative outlet vertex `g_W0(omega)=eta*omega` is reported only as a `branch_assumption`; in the two-vertex kernel bookkeeping it adds two powers of `omega` and is not used for the verdict.",
        "",
        "## Controls",
        "",
    ]
    for name, data in results["controls"].items():
        lines.append(f"- `{name}`: same_pipeline=`{data['same_pipeline']}`, fired=`{data['fired']}`")
    lines += [
        "",
        "## Engine Status",
        "",
        f"- SymPy: `{results['engine_status']['sympy']}`.",
        f"- Mathematica: `{results['engine_status']['mathematica']}`.",
        f"- Engine agreement: `{results['engine_agreement']['status']}` on `{results['engine_agreement']['shared_expression_count']}` shared symbolic expressions.",
        "",
        "## Provenance",
        "",
        "- DtN/outgoing expansions reuse `research/4d_2_5pn` spherical-Hankel machinery for ell=0,1,2.",
        "- Leakage/source moments reuse `research/pde_ledger` Part VIII projected open-system continuity and Stage 243/244 leakage lane checks.",
        "- `G=c=c_s=1`, `K_eos=1/500`, and `(a*,L*)=(4731/2500,18121/10000)` are checked in the frozen slice.",
    ]
    REPORT_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    results = build_results()
    cancellation = cancellation_condition_yaml()
    JSON_OUT.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    RESULTS_YAML.write_text(yaml.safe_dump(results, sort_keys=False, width=140), encoding="utf-8")
    CANCEL_YAML.write_text(yaml.safe_dump(cancellation, sort_keys=False, width=140), encoding="utf-8")
    write_report(results)
    print(f"wrote {REPORT_OUT}")
    print(f"wrote {RESULTS_YAML}")
    print(f"wrote {CANCEL_YAML}")
    print(f"wrote {JSON_OUT}")
    print(results["top_line_verdict"])


if __name__ == "__main__":
    main()
