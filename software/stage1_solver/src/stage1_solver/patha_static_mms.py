"""Chunk-1a manufactured-solution validation for the Path-A static balance."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
import subprocess
from typing import Any

import numpy as np
import sympy as sp
import torch

from .backend import configure_backend, tensor
from .boundaries import BoundaryCondition
from .config import BackendConfig, config_hash_from_dict
from .grid import WallGrid, WallGridSpec
from .mms import run_convergence_study
from .patha_static_balance import (
    SSigmaSpec,
    resolve_s_sigma,
    static_balance_operator,
    static_balance_terms,
)


@dataclass(frozen=True)
class PathAStaticBalanceMMSConfig:
    name: str = "pathA_chunk1a_static_balance_mms"
    w_min: float = 0.0
    w_max: float = 1.4
    grid_levels: tuple[int, ...] = (32, 64, 128, 256)
    s_sigma_spec: SSigmaSpec = field(default_factory=SSigmaSpec.patha_static_mms)
    min_observed_order: float = 1.85
    final_error_max: float = 5.0e-4
    min_nonzero_term_linf: float = 1.0e-5
    dual_engine_abs_tol: float = 2.0e-11
    run_root: str = "software/stage1_solver/runs/pathA_chunk1a_static_balance_mms"
    report_path: str = "software/stage1_solver/pathA_chunk1a_static_balance_mms_report.md"
    backend: BackendConfig = field(default_factory=BackendConfig)

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["s_sigma_spec"] = self.s_sigma_spec.to_dict()
        data["newton"] = {
            "nonlinear_solve": "not_used_for_chunk1a_mms",
            "operator_application": "manufactured_field_only",
        }
        return data

    def config_hash(self) -> str:
        return config_hash_from_dict(self.to_dict())


def _torch_from_numpy(values: np.ndarray, *, dtype: torch.dtype, device: str) -> torch.Tensor:
    return tensor(np.asarray(values), dtype=dtype, device=device)


def _stage_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _repo_root() -> Path:
    return _stage_root().parents[1]


def _check(name: str, passed: bool, detail: str) -> dict[str, Any]:
    return {"check": name, "pass": bool(passed), "detail": detail}


def _parameters(spec: SSigmaSpec) -> dict[str, float]:
    return dict(spec.parameters)


def _mms_sympy_expressions(spec: SSigmaSpec) -> dict[str, Any]:
    params = _parameters(spec)
    w, R = sp.symbols("w R", real=True)
    w_min = sp.Float(params["w_min"])
    w_max = sp.Float(params["w_max"])
    length = w_max - w_min
    x = (w - w_min) / length

    bump = 256 * x**4 * (1 - x) ** 4
    R0 = sp.Float(1.08) + sp.Float(0.12) * bump + sp.Float(0.025) * sp.sin(
        sp.pi * x
    ) ** 4
    T_w = (
        sp.Float(params["tw_base"])
        + sp.Float(params["tw_w_sin"]) * sp.sin(2 * sp.pi * x)
        + sp.Float(params["tw_r1"]) * R
        + sp.Float(params["tw_r2"]) * R**2
        + sp.Float(params["tw_rw"]) * R * sp.cos(sp.pi * x)
    )
    U = (
        sp.Float(params["u_base"])
        + sp.Float(params["u_r2"]) * R**2
        + sp.Float(params["u_r3"]) * R**3
        + sp.Float(params["u_rw"]) * R * sp.sin(sp.pi * x)
    )
    R0_w = sp.diff(R0, w)
    T_w_R = sp.diff(T_w, R)
    U_R = sp.diff(U, R)
    flux_divergence = -sp.diff(T_w.subs(R, R0) * R0_w, w)
    gradient_square = sp.Rational(1, 2) * T_w_R.subs(R, R0) * R0_w**2
    potential_gradient = U_R.subs(R, R0)
    source = flux_divergence + gradient_square + potential_gradient

    return {
        "w_symbol": w,
        "R_symbol": R,
        "R0_expr": R0,
        "T_w_expr": T_w,
        "T_w_R_expr": T_w_R,
        "U_expr": U,
        "U_R_expr": U_R,
        "flux_divergence_expr": flux_divergence,
        "gradient_square_expr": gradient_square,
        "potential_gradient_expr": potential_gradient,
        "source_expr": source,
        "R0": sp.lambdify(w, R0, "numpy"),
        "source": sp.lambdify(w, source, "numpy"),
        "flux_divergence": sp.lambdify(w, flux_divergence, "numpy"),
        "gradient_square": sp.lambdify(w, gradient_square, "numpy"),
        "potential_gradient": sp.lambdify(w, potential_gradient, "numpy"),
    }


def _mathematica_script_path() -> Path:
    return _stage_root() / "mathematica" / "pathA_02_chunk1a_static_balance_mms.wls"


def _mathematica_diagnostics_path() -> Path:
    return (
        _stage_root()
        / "mathematica"
        / "runs"
        / "pathA_02_chunk1a_static_balance_mms"
        / "pathA_02_chunk1a_static_balance_mms_diagnostics.json"
    )


def _run_mathematica_mms() -> dict[str, Any]:
    script = _mathematica_script_path()
    if not script.exists():
        raise FileNotFoundError(f"Mathematica verifier missing: {script}")
    completed = subprocess.run(
        ["timeout", "600", "wolframscript", "-script", str(script)],
        cwd=_repo_root(),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "Mathematica static-balance MMS verifier failed\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    diagnostics_path = _mathematica_diagnostics_path()
    return json.loads(diagnostics_path.read_text(encoding="utf-8"))


def _dual_engine_check(
    expressions: dict[str, Any],
    mathematica: dict[str, Any],
    *,
    abs_tol: float,
) -> dict[str, Any]:
    sample_w = np.asarray(mathematica["sample_w"], dtype=np.float64)
    sympy_source = np.asarray(expressions["source"](sample_w), dtype=np.float64)
    math_source = np.asarray(mathematica["source_samples"], dtype=np.float64)
    diff = np.abs(sympy_source - math_source)
    term_diffs: dict[str, float] = {}
    for key, expr_key in (
        ("flux_divergence", "flux_divergence"),
        ("gradient_square", "gradient_square"),
        ("potential_gradient", "potential_gradient"),
    ):
        sympy_term = np.asarray(expressions[expr_key](sample_w), dtype=np.float64)
        math_term = np.asarray(mathematica["term_samples"][key], dtype=np.float64)
        term_diffs[key] = float(np.max(np.abs(sympy_term - math_term)))
    max_abs_diff = float(np.max(diff))
    return {
        "sample_w": sample_w.tolist(),
        "sympy_source_samples": sympy_source.tolist(),
        "mathematica_source_samples": math_source.tolist(),
        "max_abs_diff": max_abs_diff,
        "term_max_abs_diff": term_diffs,
        "passed": max_abs_diff <= abs_tol and all(value <= abs_tol for value in term_diffs.values()),
        "diagnostics_path": str(_mathematica_diagnostics_path()),
    }


def _forbidden_target_tokens() -> list[str]:
    return [
        "".join(("R", "_", "norm")),
        "54" + "/" + "5",
        "10" + "." + "8",
        "P0" + "_" + "target",
        "GR" + "-" + "constant",
        "GR" + " " + "constant",
    ]


def _target_token_scan(paths: list[Path]) -> dict[str, Any]:
    hits: list[dict[str, str]] = []
    tokens = _forbidden_target_tokens()
    for path in paths:
        if not path.exists() or path.is_dir():
            continue
        text = path.read_text(encoding="utf-8")
        for token in tokens:
            if token in text:
                hits.append({"path": str(path), "token": token})
    return {"passed": not hits, "hits": hits, "path_count": len(paths)}


def _scan_paths(report_path: Path | None = None) -> list[Path]:
    stage_root = _stage_root()
    paths = [
        stage_root / "src" / "stage1_solver" / "patha_static_balance.py",
        stage_root / "src" / "stage1_solver" / "patha_static_mms.py",
        stage_root / "tests" / "test_patha_static_balance.py",
        _mathematica_script_path(),
    ]
    if report_path is not None:
        paths.append(report_path)
    return paths


def run_patha_static_balance_mms(
    config: PathAStaticBalanceMMSConfig | None = None,
    *,
    run_mathematica: bool = True,
    write_report: bool = True,
) -> dict[str, Any]:
    cfg = config or PathAStaticBalanceMMSConfig()
    dtype = configure_backend(cfg.backend)
    expressions = _mms_sympy_expressions(cfg.s_sigma_spec)
    provider = resolve_s_sigma(cfg.s_sigma_spec)
    full_config = cfg.to_dict()
    lower_value = float(expressions["R0"](cfg.w_min))
    upper_value = float(expressions["R0"](cfg.w_max))
    lower_bc = BoundaryCondition.dirichlet(lower_value)
    upper_bc = BoundaryCondition.dirichlet(upper_value)

    def build_level(nw: int):
        grid = WallGrid.create(
            WallGridSpec(w_min=cfg.w_min, w_max=cfg.w_max, nw=nw),
            dtype=dtype,
            device=cfg.backend.device,
        )
        return grid, f"nw_{nw}", grid.dw, grid.to_dict(), grid.cell_widths

    def evaluate_level(grid: WallGrid):
        w_centers = grid.w_centers.detach().cpu().numpy()
        radius_np = expressions["R0"](w_centers)
        source_np = expressions["source"](w_centers)
        radius = _torch_from_numpy(radius_np, dtype=dtype, device=cfg.backend.device)
        source = _torch_from_numpy(source_np, dtype=dtype, device=cfg.backend.device)
        terms = static_balance_terms(
            radius,
            grid,
            s_sigma=provider,
            source=source,
            lower_bc=lower_bc,
            upper_bc=upper_bc,
        )
        residual = static_balance_operator(
            radius,
            grid,
            s_sigma=provider,
            source=source,
            lower_bc=lower_bc,
            upper_bc=upper_bc,
        )
        diagnostics = {
            "residual_linf": float(torch.max(torch.abs(residual)).detach().cpu().item()),
            "flux_divergence_linf": float(
                torch.max(torch.abs(terms.flux_divergence)).detach().cpu().item()
            ),
            "gradient_square_linf": float(
                torch.max(torch.abs(terms.gradient_square)).detach().cpu().item()
            ),
            "potential_gradient_linf": float(
                torch.max(torch.abs(terms.potential_gradient)).detach().cpu().item()
            ),
            "lower_bc": lower_bc.to_dict(),
            "upper_bc": upper_bc.to_dict(),
        }
        return terms.lhs, source, diagnostics

    convergence = asdict(
        run_convergence_study(
            name=cfg.name,
            description=(
                "Path-A chunk-1a nonlinear static throat balance with R-dependent "
                "T_w, live gradient-square term, and nonzero U_R."
            ),
            continuum_source=(
                "-d_w(T_w(R0,w) d_w R0) + 0.5 T_w_R(R0,w)(d_w R0)^2 + U_R(R0,w)"
            ),
            manufactured_field=str(expressions["R0_expr"]),
            forcing_derivation=(
                "SymPy differentiates the continuum balance analytically; Mathematica "
                "exports an independent sample check of the same continuum expression."
            ),
            levels=cfg.grid_levels,
            build_level=build_level,
            evaluate_level=evaluate_level,
            config=full_config,
            run_root=cfg.run_root,
            min_observed_order=cfg.min_observed_order,
            final_error_max=cfg.final_error_max,
            config_hash=cfg.config_hash(),
        )
    )

    mathematica = _run_mathematica_mms() if run_mathematica else {}
    dual_engine = (
        _dual_engine_check(expressions, mathematica, abs_tol=cfg.dual_engine_abs_tol)
        if run_mathematica
        else {"passed": False, "max_abs_diff": None, "diagnostics_path": None}
    )
    final_row = convergence["rows"][-1]
    term_nonzero = {
        "flux_divergence": final_row["flux_divergence_linf"] >= cfg.min_nonzero_term_linf,
        "gradient_square": final_row["gradient_square_linf"] >= cfg.min_nonzero_term_linf,
        "potential_gradient": final_row["potential_gradient_linf"] >= cfg.min_nonzero_term_linf,
    }

    report_path = Path(cfg.report_path)
    target_scan = _target_token_scan(_scan_paths(report_path if report_path.exists() else None))
    genuine_checks = [
        _check(
            "dual_engine_forcing_agreement",
            bool(dual_engine["passed"]),
            f"SymPy and Mathematica source samples agree to {dual_engine['max_abs_diff']:.3e}.",
        ),
        _check(
            "second_order_static_balance_mms",
            bool(convergence["pass_checks"]["observed_order"])
            and bool(convergence["pass_checks"]["final_error"]),
            "Discrete LHS compared with independently derived continuum source on a refinement ladder.",
        ),
        _check(
            "flux_divergence_term_nonzero",
            term_nonzero["flux_divergence"],
            f"final L_inf={final_row['flux_divergence_linf']:.3e}.",
        ),
        _check(
            "gradient_square_term_nonzero",
            term_nonzero["gradient_square"],
            f"final L_inf={final_row['gradient_square_linf']:.3e}.",
        ),
        _check(
            "potential_gradient_term_nonzero",
            term_nonzero["potential_gradient"],
            f"final L_inf={final_row['potential_gradient_linf']:.3e}.",
        ),
        _check(
            "target_token_scan_clean",
            bool(target_scan["passed"]),
            f"scanned {target_scan['path_count']} chunk-1a code/report files.",
        ),
    ]
    restatement_checks = [
        _check(
            "spec_roundtrip_not_a_physics_gate",
            SSigmaSpec.from_dict(cfg.s_sigma_spec.to_dict()) == cfg.s_sigma_spec,
            "Serializable registry selector round-trips without callables.",
        ),
        _check(
            "dirichlet_end_values_not_a_physics_gate",
            lower_bc.kind == "dirichlet" and upper_bc.kind == "dirichlet",
            "Chunk-1a MMS prescribes both manufactured end values; open-exit choice is deferred.",
        ),
    ]

    pass_checks = dict(convergence["pass_checks"])
    pass_checks.update(
        {
            "dual_engine_forcing_agreement": bool(dual_engine["passed"]),
            "target_token_scan_clean": bool(target_scan["passed"]),
            **{f"{name}_nonzero": value for name, value in term_nonzero.items()},
        }
    )
    result = {
        "name": cfg.name,
        "config": full_config,
        "config_hash": cfg.config_hash(),
        "s_sigma_spec": cfg.s_sigma_spec.to_dict(),
        "s_sigma_digest": cfg.s_sigma_spec.digest(),
        "convergence": convergence,
        "dual_engine": dual_engine,
        "genuine_checks": genuine_checks,
        "restatement_checks_not_a_physics_gate": restatement_checks,
        "target_token_scan": target_scan,
        "pass_checks": pass_checks,
        "passed": all(pass_checks.values())
        and all(row["pass"] for row in genuine_checks)
        and all(row["pass"] for row in restatement_checks),
    }

    if write_report:
        write_patha_static_balance_report(result, report_path)
        target_scan = _target_token_scan(_scan_paths(report_path))
        result["target_token_scan"] = target_scan
        result["pass_checks"]["target_token_scan_clean"] = bool(target_scan["passed"])
        for row in result["genuine_checks"]:
            if row["check"] == "target_token_scan_clean":
                row["pass"] = bool(target_scan["passed"])
                row["detail"] = f"scanned {target_scan['path_count']} chunk-1a code/report files."
        result["passed"] = all(result["pass_checks"].values()) and all(
            row["pass"] for row in result["genuine_checks"]
        ) and all(row["pass"] for row in result["restatement_checks_not_a_physics_gate"])
        write_patha_static_balance_report(result, report_path)

    return result


def _fmt(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.6e}"
    return str(value)


def _table(headers: list[str], rows: list[dict[str, Any]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(_fmt(row.get(header)) for header in headers) + " |")
    return "\n".join(lines)


def write_patha_static_balance_report(results: dict[str, Any], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    convergence = results["convergence"]
    dual = results["dual_engine"]
    final_row = convergence["rows"][-1]
    lines: list[str] = []
    lines.append("# Path-A Chunk 1a Static Balance MMS Report")
    lines.append("")
    lines.append(f"Overall gate: {'PASS' if results['passed'] else 'FAIL'}")
    lines.append(f"Config hash: `{results['config_hash']}`")
    lines.append(f"S_Sigma spec digest: `{results['s_sigma_digest']}`")
    lines.append("")
    lines.append("## Built Surface")
    lines.append("")
    lines.append(
        "Implemented a serializable `S_Sigma` spec/registry and a separate nonlinear "
        "flat-`dw` FV static-balance operator for `R0(w)` with Dirichlet face values."
    )
    lines.append("")
    lines.append("## Genuine Gates")
    lines.append("")
    lines.append(_table(["check", "pass", "detail"], results["genuine_checks"]))
    lines.append("")
    lines.append("## Construction Restatements")
    lines.append("")
    lines.append(_table(["check", "pass", "detail"], results["restatement_checks_not_a_physics_gate"]))
    lines.append("")
    lines.append("## MMS Refinement")
    lines.append("")
    lines.append(convergence["description"])
    lines.append(f"Continuum source: `{convergence['continuum_source']}`")
    lines.append(f"Manufactured field: `{convergence['manufactured_field']}`")
    lines.append(f"Forcing derivation: {convergence['forcing_derivation']}")
    lines.append("")
    lines.append(
        _table(
            [
                "grid",
                "spacing",
                "error",
                "observed_order",
                "reference_norm",
                "residual_linf",
                "flux_divergence_linf",
                "gradient_square_linf",
                "potential_gradient_linf",
            ],
            convergence["rows"],
        )
    )
    lines.append("")
    lines.append(
        "Finest-grid summary: "
        f"error={final_row['error']:.6e}, observed_order={final_row['observed_order']:.6e}."
    )
    lines.append("")
    lines.append("## Dual Engine")
    lines.append("")
    lines.append(
        "SymPy and Mathematica independently evaluated the continuum manufactured "
        f"source samples; max absolute difference = {_fmt(dual['max_abs_diff'])}."
    )
    gradient_square_term_diff = dual.get("term_max_abs_diff", {}).get("gradient_square")
    lines.append(
        "Because the finest-grid gradient-square L_inf "
        f"({_fmt(final_row['gradient_square_linf'])}) is small relative to flux-divergence "
        f"L_inf ({_fmt(final_row['flux_divergence_linf'])}), that term's correctness is "
        "validated primarily by the dual-engine continuum-term agreement "
        f"(`term_max_abs_diff.gradient_square` = {_fmt(gradient_square_term_diff)}) and "
        "its live residual contribution, not by the global second-order residual metric alone."
    )
    lines.append(f"Mathematica diagnostics: `{dual['diagnostics_path']}`")
    lines.append("")
    lines.append("## Validation Scope")
    lines.append("")
    lines.append(
        "This MMS exercises only the `patha_static_mms_v1` provider; the "
        "`smooth_positive_placeholder_v1` family derivatives (`T_w_R`, `T_w_RR`, `U_R`, "
        "`U_RR`) are analytically asserted and hand-verified, with MMS validation of the "
        "actual solve family deferred to chunk 1b/1c."
    )
    lines.append("")
    lines.append("## Under-Specified Items")
    lines.append("")
    lines.append(
        "None for chunk 1a. The background subtraction, radial reduction, and open-exit "
        "BC choices remain deferred to chunk 1b and are not used by this MMS."
    )
    lines.append("")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main() -> int:
    result = run_patha_static_balance_mms()
    print(f"wrote Path-A chunk-1a report: {Path(result['config']['report_path'])}")
    if not result["passed"]:
        print("Path-A chunk-1a static-balance gate failed")
        return 1
    print("Path-A chunk-1a static-balance gate passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
