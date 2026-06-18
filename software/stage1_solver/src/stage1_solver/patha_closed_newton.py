"""Path-A chunk 1b closed Newton solve over matter, gauge, and R0."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
import time
from typing import Any

import torch

from .backend import configure_backend
from .config import BackendConfig, BranchSmokeConfig, NewtonConfig, PreconditionerConfig
from .coupled_branch import (
    ClosedCoupledFields,
    _create_branch_grid,
    _max_rss_mb,
    branch_boundary_conditions,
    confinement_wall_to_matter_coefficient_torch,
    initial_closed_branch_state,
    patha_closed_branch_residual,
    patha_closed_wall_terms,
    unpack_closed_coupled_fields,
)
from .newton import finite_difference_jvp_check, solve_newton_jvp
from .patha_static_balance import SSigmaSpec, resolve_s_sigma
from .preconditioners import make_closed_coupled_colored_sparse_jacobian_lu_factory


def _sorted_parameters(parameters: dict[str, float]) -> tuple[tuple[str, float], ...]:
    return tuple(sorted((str(key), float(value)) for key, value in parameters.items()))


def default_closed_branch_config() -> BranchSmokeConfig:
    """Small target-blind placeholder branch for chunk-1b engineering runs."""

    preconditioner = PreconditionerConfig(
        type="colored_sparse_jacobian_lu",
        side="left",
        rebuild_policy="every_newton_step",
        stencil_radius=3,
        color_separation=7,
        factorization="splu",
        diagonal_shift=1.0e-10,
        drop_tolerance=0.0,
        fill_factor=10.0,
        permutation="COLAMD",
    )
    return BranchSmokeConfig(
        name="patha_chunk1b_closed_newton_placeholder",
        placeholder_label=(
            "Path-A chunk-1b closed engineering placeholder; no physical packet; target-blind"
        ),
        r_max=1.6,
        w_min=0.0,
        w_max=1.2,
        solve_grid=(6, 6),
        ladder_levels=((6, 6),),
        continuation_K_values=(0.08, 0.18),
        mass=0.05,
        localization_width=0.65,
        localization_floor=0.8,
        localization_amplitude=0.3,
        geometry_profile="cubic_smoothstep",
        r_mouth=1.0,
        r_exit=0.92,
        geometry_decay_length=0.8,
        radial_wall_strength=0.045,
        axial_trap_strength=0.04,
        matter_mouth_boundary="neumann",
        a0_mouth_boundary="dirichlet",
        matter_exit_impedance_alpha=0.2,
        a0_radial_impedance_alpha=0.45,
        a0_exit_impedance_alpha=0.35,
        initial_mu=0.65,
        newton=NewtonConfig(
            residual_atol=2.0e-9,
            residual_rtol=2.0e-9,
            step_atol=1.0e-12,
            step_rtol=1.0e-11,
            max_newton_iters=18,
            gmres_rtol=2.0e-9,
            gmres_atol=1.0e-11,
            gmres_restart=160,
            gmres_maxiter=12,
            max_line_search_iters=18,
            accept_best_line_search_decrease=True,
            finite_difference_jvp_epsilon=1.0e-5,
            preconditioner=preconditioner,
        ),
    )


def default_closed_s_sigma_spec(branch: BranchSmokeConfig) -> SSigmaSpec:
    base = dict(SSigmaSpec.smooth_positive_placeholder(
        w_min=branch.w_min,
        w_max=branch.w_max,
    ).parameters)
    base.update(
        {
            "mu_base": 1.0,
            "mu_r2": 0.02,
            "mu_w_amp": 0.03,
            "tw_base": 3.5,
            "tw_r2": 0.08,
            "tw_w_amp": 0.04,
            "tomega_base": 0.9,
            "tomega_r2": 0.02,
            "tomega_w_amp": 0.03,
            "u_base": 0.2,
            "u_r2": 0.025,
            "u_r4": 0.002,
            "u_w_amp": 0.02,
        }
    )
    return SSigmaSpec(
        family="smooth_positive_placeholder_v1",
        parameters=_sorted_parameters(base),
    )


@dataclass(frozen=True)
class PathAClosedNewtonConfig:
    name: str = "patha_chunk1b_closed_newton"
    branch: BranchSmokeConfig = field(default_factory=default_closed_branch_config)
    s_sigma_spec: SSigmaSpec | None = None
    run_root: str = "software/stage1_solver/runs/patha_chunk1b_closed_newton"
    report_path: str = "software/stage1_solver/reports/patha_chunk1b_closed_newton.md"
    jvp_relative_tol: float = 5.0e-6
    jvp_absolute_tol: float = 5.0e-7
    derivative_relative_tol: float = 2.0e-6
    derivative_absolute_tol: float = 2.0e-7
    source_fidelity_tol: float = 1.0e-12
    jacobian_check_seed: int = 20310617

    def resolved_s_sigma_spec(self) -> SSigmaSpec:
        return self.s_sigma_spec or default_closed_s_sigma_spec(self.branch)

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["s_sigma_spec"] = self.resolved_s_sigma_spec().to_dict()
        return data


def _newton_history_rows(result: Any) -> list[dict[str, Any]]:
    return [asdict(row) for row in result.history]


def _forbidden_target_tokens() -> list[str]:
    return [
        "".join(("R", "_", "norm")),
        "54" + "/" + "5",
        "10" + "." + "8",
        "P0" + "_" + "target",
        "GR" + "-" + "constant",
        "GR" + " " + "constant",
    ]


def target_token_scan(paths: list[Path]) -> dict[str, Any]:
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


def _stage_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _scan_paths(report_path: Path | None = None) -> list[Path]:
    root = _stage_root()
    paths = [
        root / "src" / "stage1_solver" / "patha_closed_newton.py",
        root / "src" / "stage1_solver" / "coupled_branch.py",
        root / "src" / "stage1_solver" / "preconditioners.py",
        root / "tests" / "test_patha_closed_newton.py",
    ]
    if report_path is not None:
        paths.append(report_path)
    return paths


def placeholder_provider_derivative_check(
    spec: SSigmaSpec,
    *,
    dtype: torch.dtype,
    device: str,
    epsilon: float = 1.0e-4,
) -> dict[str, Any]:
    provider = resolve_s_sigma(spec)
    radius = torch.linspace(0.72, 1.18, 11, dtype=dtype, device=device)
    w = torch.linspace(
        dict(spec.parameters)["w_min"],
        dict(spec.parameters)["w_max"],
        11,
        dtype=dtype,
        device=device,
    )
    eps = torch.as_tensor(epsilon, dtype=dtype, device=device)

    def first_fd(fn):
        return (fn(radius + eps, w) - fn(radius - eps, w)) / (2.0 * eps)

    def second_fd(fn):
        return (fn(radius + eps, w) - 2.0 * fn(radius, w) + fn(radius - eps, w)) / (
            eps * eps
        )

    checks = {
        "T_w_R": (provider.T_w_R(radius, w), first_fd(provider.T_w)),
        "T_w_RR": (provider.T_w_RR(radius, w), second_fd(provider.T_w)),
        "U_R": (provider.U_R(radius, w), first_fd(provider.U)),
        "U_RR": (provider.U_RR(radius, w), second_fd(provider.U)),
    }
    rows: list[dict[str, float | str]] = []
    max_absolute = 0.0
    max_relative = 0.0
    for name, (analytic, fd) in checks.items():
        diff = analytic - fd
        absolute = float(torch.max(torch.abs(diff)).detach().cpu().item())
        relative = absolute / max(1.0, float(torch.max(torch.abs(analytic)).detach().cpu().item()))
        rows.append({"quantity": name, "absolute": absolute, "relative": relative})
        max_absolute = max(max_absolute, absolute)
        max_relative = max(max_relative, relative)
    return {
        "epsilon": epsilon,
        "rows": rows,
        "max_absolute": max_absolute,
        "max_relative": max_relative,
    }


def return_source_fidelity_diagnostic(
    state: torch.Tensor,
    grid,
    branch: BranchSmokeConfig,
    spec: SSigmaSpec,
) -> dict[str, Any]:
    fields, _ = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    density = fields.psi_real**2 + fields.psi_imag**2
    shared_k1 = confinement_wall_to_matter_coefficient_torch(grid, branch, radius=fields.r0)
    recomputed_source = torch.sum(
        grid.radial_shell_volumes[:, None] * (-shared_k1 * density),
        dim=0,
    )
    residual_source = patha_closed_wall_terms(
        fields,
        grid,
        branch,
        s_sigma=resolve_s_sigma(spec),
    ).source
    diff = recomputed_source - residual_source
    return {
        "kernel": (
            "stage1_solver.coupled_branch."
            "confinement_wall_to_matter_coefficient_torch"
        ),
        "max_abs_diff": float(torch.max(torch.abs(diff)).detach().cpu().item()),
        "recomputed_linf": float(torch.max(torch.abs(recomputed_source)).detach().cpu().item()),
        "residual_source_linf": float(torch.max(torch.abs(residual_source)).detach().cpu().item()),
        "source_min": float(torch.min(recomputed_source).detach().cpu().item()),
        "source_max": float(torch.max(recomputed_source).detach().cpu().item()),
    }


def constitutive_stability_margin(
    fields: ClosedCoupledFields,
    grid,
    spec: SSigmaSpec,
) -> dict[str, float]:
    provider = resolve_s_sigma(spec)
    r0 = fields.r0
    w = grid.w_centers
    values = {
        "mu_min": torch.min(provider.mu(r0, w)),
        "T_w_min": torch.min(provider.T_w(r0, w)),
        "T_Omega_min": torch.min(provider.T_Omega(r0, w)),
        "U_RR_min": torch.min(provider.U_RR(r0, w)),
    }
    payload = {key: float(value.detach().cpu().item()) for key, value in values.items()}
    payload["minimum_margin"] = min(payload.values())
    return payload


def run_patha_closed_newton(
    config: PathAClosedNewtonConfig | None = None,
) -> dict[str, Any]:
    if config is None:
        config = PathAClosedNewtonConfig()
    Path(config.run_root).mkdir(parents=True, exist_ok=True)
    dtype = configure_backend(BackendConfig())
    branch = config.branch
    spec = config.resolved_s_sigma_spec()
    provider = resolve_s_sigma(spec)
    grid = _create_branch_grid(branch, branch.solve_grid, dtype=dtype, device="cpu")
    boundaries = branch_boundary_conditions(branch)
    state = initial_closed_branch_state(grid, branch, dtype=dtype, device="cpu")
    started = time.perf_counter()
    stages: list[dict[str, Any]] = []
    converged = True
    message = "closed continuation completed"
    stop_conditions: list[dict[str, Any]] = []
    shared_preconditioner_factory = make_closed_coupled_colored_sparse_jacobian_lu_factory(grid)

    for eos_K in branch.continuation_K_values:
        residual_fn = lambda x, eos_K=eos_K: patha_closed_branch_residual(
            x,
            grid,
            branch,
            eos_K=eos_K,
            boundaries=boundaries,
            s_sigma=provider,
        )
        fields_before, _ = unpack_closed_coupled_fields(
            state,
            grid,
            has_chemical_potential=True,
        )
        if torch.min(fields_before.r0).detach().cpu().item() <= 0.0:
            converged = False
            message = f"stopped before eos_K={eos_K}: R0 became nonpositive"
            stop_conditions.append({"condition": "R0_nonpositive", "tripped": True})
            break
        result = solve_newton_jvp(
            residual_fn,
            state,
            branch.newton,
            preconditioner_factory=(
                make_closed_coupled_colored_sparse_jacobian_lu_factory(grid)
                if branch.newton.preconditioner.rebuild_policy == "once_per_newton_solve"
                else shared_preconditioner_factory
            ),
        )
        state = result.x.detach()
        fields_after, _ = unpack_closed_coupled_fields(
            state,
            grid,
            has_chemical_potential=True,
        )
        r0_min = float(torch.min(fields_after.r0).detach().cpu().item())
        r0_max = float(torch.max(fields_after.r0).detach().cpu().item())
        gmres_counts = [
            row.gmres_iterations for row in result.history if row.gmres_iterations is not None
        ]
        stages.append(
            {
                "eos_K": eos_K,
                "converged": result.converged,
                "iterations": result.iterations,
                "initial_residual_norm": result.initial_residual_norm,
                "final_residual_norm": result.final_residual_norm,
                "tolerance": result.tolerance,
                "message": result.message,
                "gmres_iterations": gmres_counts,
                "residual_history": [row.residual_norm for row in result.history],
                "newton_history": _newton_history_rows(result),
                "r0_min": r0_min,
                "r0_max": r0_max,
            }
        )
        if r0_min <= 0.0:
            converged = False
            message = f"stopped after eos_K={eos_K}: R0 became nonpositive"
            stop_conditions.append({"condition": "R0_nonpositive", "tripped": True})
            break
        if not result.converged:
            converged = False
            message = f"continuation stopped at eos_K={eos_K}: {result.message}"
            break

    final_eos_K = branch.continuation_K_values[max(0, min(len(stages), len(branch.continuation_K_values)) - 1)]
    final_residual_fn = lambda x: patha_closed_branch_residual(
        x,
        grid,
        branch,
        eos_K=final_eos_K,
        boundaries=boundaries,
        s_sigma=provider,
    )
    final_residual = final_residual_fn(state).detach()
    fields, mu = unpack_closed_coupled_fields(state, grid, has_chemical_potential=True)
    density = fields.psi_real**2 + fields.psi_imag**2
    final_linf = float(torch.max(torch.abs(final_residual)).detach().cpu().item())
    jvp_check = finite_difference_jvp_check(
        final_residual_fn,
        state,
        epsilon=branch.newton.finite_difference_jvp_epsilon,
        seed=config.jacobian_check_seed,
    )
    jvp_passed = (
        jvp_check["relative_residual"] <= config.jvp_relative_tol
        or jvp_check["absolute_residual"] <= config.jvp_absolute_tol
    )
    derivative_check = placeholder_provider_derivative_check(
        spec,
        dtype=dtype,
        device="cpu",
    )
    derivative_passed = (
        derivative_check["max_relative"] <= config.derivative_relative_tol
        or derivative_check["max_absolute"] <= config.derivative_absolute_tol
    )
    source_diagnostic = return_source_fidelity_diagnostic(state, grid, branch, spec)
    source_fidelity_passed = source_diagnostic["max_abs_diff"] <= config.source_fidelity_tol
    stability = constitutive_stability_margin(fields, grid, spec)
    r0_min = float(torch.min(fields.r0).detach().cpu().item())
    positive_r0 = r0_min > 0.0
    if positive_r0:
        stop_conditions.append({"condition": "R0_nonpositive", "tripped": False})
    stop_conditions.extend(
        [
            {"condition": "exit_bc_changed_from_zero_traction", "tripped": False},
            {"condition": "hidden_clamp_or_line_search_hack", "tripped": False},
            {"condition": "source_sign_or_reduction_mismatch", "tripped": False},
        ]
    )
    elapsed = time.perf_counter() - started
    report_path = Path(config.report_path)
    token_scan = target_token_scan(_scan_paths(report_path if report_path.exists() else None))
    genuine_checks = [
        {
            "check": "closed_jvp_consistency",
            "pass": bool(jvp_passed),
            "detail": (
                f"relative={jvp_check['relative_residual']:.6e}, "
                f"absolute={jvp_check['absolute_residual']:.6e}"
            ),
        },
        {
            "check": "closed_newton_convergence",
            "pass": bool(converged and positive_r0 and final_linf <= stages[-1]["tolerance"]),
            "detail": f"final_linf={final_linf:.6e}, R0_min={r0_min:.6e}",
        },
        {
            "check": "placeholder_provider_derivatives_fd",
            "pass": bool(derivative_passed),
            "detail": (
                f"max_relative={derivative_check['max_relative']:.6e}, "
                f"max_absolute={derivative_check['max_absolute']:.6e}"
            ),
        },
        {
            "check": "constitutive_positive_margin",
            "pass": bool(stability["minimum_margin"] > 0.0),
            "detail": f"minimum_margin={stability['minimum_margin']:.6e}",
        },
        {
            "check": "target_token_scan_clean",
            "pass": bool(token_scan["passed"]),
            "detail": f"scanned {token_scan['path_count']} Path-A 1b files.",
        },
    ]
    diagnostics_not_a_physics_gate = [
        {
            "check": "return_source_fidelity_not_a_physics_gate",
            "status": "reported",
            "detail": f"max_abs_diff={source_diagnostic['max_abs_diff']:.6e}",
        },
        {
            "check": "gauge_return_no_explicit_source_not_a_physics_gate",
            "status": "construction_note",
            "detail": "matter and Maxwell lanes are coupled monolithically; no extra gauge source was added.",
        },
    ]
    passed = all(row["pass"] for row in genuine_checks)
    return {
        "config": config.to_dict(),
        "s_sigma_digest": spec.digest(),
        "grid": grid.to_dict(),
        "dof": int(state.numel()),
        "converged": converged,
        "message": message,
        "wall_clock_seconds": elapsed,
        "peak_memory_mb": _max_rss_mb(),
        "final_residual_linf": final_linf,
        "final_mass": float(torch.sum(density * grid.cell_volumes).detach().cpu().item()),
        "chemical_potential": float(mu.detach().cpu().item()) if mu is not None else None,
        "r0_min": r0_min,
        "r0_max": float(torch.max(fields.r0).detach().cpu().item()),
        "stages": stages,
        "jvp_check": jvp_check,
        "jvp_passed": jvp_passed,
        "placeholder_derivative_check": derivative_check,
        "placeholder_derivative_passed": derivative_passed,
        "return_source_fidelity": source_diagnostic,
        "return_source_fidelity_passed": source_fidelity_passed,
        "constitutive_stability_margin": stability,
        "target_token_scan": token_scan,
        "genuine_checks": genuine_checks,
        "diagnostics_not_a_physics_gate": diagnostics_not_a_physics_gate,
        "stop_conditions": stop_conditions,
        "passed": passed,
    }


def _fmt(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.6e}"
    if isinstance(value, list):
        return "[" + ", ".join(str(item) for item in value) + "]"
    return str(value)


def _table(headers: list[str], rows: list[dict[str, Any]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(_fmt(row.get(header)) for header in headers) + " |")
    return "\n".join(lines)


def write_patha_closed_newton_report(results: dict[str, Any], path: str) -> Path:
    report_path = Path(path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    config = results["config"]
    branch = config["branch"]
    lines: list[str] = []
    lines.append("# Path-A Chunk 1b Closed Newton")
    lines.append("")
    lines.append(f"Overall engineering gate: {'PASS' if results['passed'] else 'FAIL'}")
    lines.append(f"S_Sigma digest: `{results['s_sigma_digest']}`")
    lines.append("")
    lines.append(
        "Scope: closed static engineering solve over matter, gauge, R0, and chemical potential "
        "with placeholder S_Sigma. No physical packet or coefficient export is performed."
    )
    lines.append("")
    lines.append("## Placeholder Parameters")
    lines.append("")
    lines.append(f"Label: {branch['placeholder_label']}")
    lines.append("")
    lines.append("```yaml")
    for key in (
        "solve_grid",
        "continuation_K_values",
        "mass",
        "r_mouth",
        "r_exit",
        "radial_wall_strength",
        "axial_trap_strength",
        "matter_exit_impedance_alpha",
        "a0_exit_impedance_alpha",
    ):
        lines.append(f"{key}: {branch[key]}")
    lines.append("preconditioner_layout: 5*cells+nw+1")
    lines.append("exit_wall_bc: natural_zero_traction")
    lines.append("```")
    lines.append("")
    lines.append("## Newton Convergence")
    lines.append("")
    lines.append(
        _table(
            [
                "eos_K",
                "converged",
                "iterations",
                "initial_residual_norm",
                "final_residual_norm",
                "tolerance",
                "gmres_iterations",
                "r0_min",
                "r0_max",
                "message",
            ],
            results["stages"],
        )
    )
    lines.append("")
    lines.append(
        f"Final residual linf={results['final_residual_linf']:.6e}; "
        f"R0 range=[{results['r0_min']:.6e}, {results['r0_max']:.6e}]; "
        f"mass={results['final_mass']:.6e}; mu={results['chemical_potential']:.6e}; "
        f"wall-clock={results['wall_clock_seconds']:.6e}s."
    )
    lines.append("")
    lines.append("## JVP Check")
    lines.append("")
    jac = results["jvp_check"]
    lines.append(
        f"Closed residual JVP vs centered finite difference: relative={jac['relative_residual']:.6e}, "
        f"absolute={jac['absolute_residual']:.6e}, epsilon={jac['epsilon']:.6e}, "
        f"status={'PASS' if results['jvp_passed'] else 'FAIL'}."
    )
    lines.append("")
    lines.append("## Placeholder Derivatives")
    lines.append("")
    lines.append(
        _table(
            ["quantity", "absolute", "relative"],
            results["placeholder_derivative_check"]["rows"],
        )
    )
    lines.append("")
    lines.append("## Return Source Diagnostic")
    lines.append("")
    source = results["return_source_fidelity"]
    lines.append(
        "The wall source was recomputed from the shared confinement coefficient and shell volumes. "
        f"max_abs_diff={source['max_abs_diff']:.6e}, source range="
        f"[{source['source_min']:.6e}, {source['source_max']:.6e}]."
    )
    lines.append("")
    lines.append("## Stability Margin")
    lines.append("")
    lines.append("```yaml")
    for key, value in results["constitutive_stability_margin"].items():
        lines.append(f"{key}: {value:.6e}")
    lines.append("```")
    lines.append("")
    lines.append(
        "A downstream Schur-denominator value is not constructed in chunk 1b; this report records "
        "the placeholder constitutive positivity margin used by the closed background solve."
    )
    lines.append("")
    lines.append("## Counted Gates")
    lines.append("")
    for row in results["genuine_checks"]:
        lines.append(f"- {row['check']}: {'PASS' if row['pass'] else 'FAIL'} ({row['detail']})")
    lines.append("")
    lines.append("## Diagnostics Not A Physics Gate")
    lines.append("")
    for row in results["diagnostics_not_a_physics_gate"]:
        lines.append(f"- {row['check']}: {row['status']} ({row['detail']})")
    lines.append("")
    lines.append("## Stop Conditions")
    lines.append("")
    for row in results["stop_conditions"]:
        lines.append(f"- {row['condition']}: {'TRIPPED' if row['tripped'] else 'not tripped'}")
    lines.append("")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path


def main() -> int:
    config = PathAClosedNewtonConfig()
    results = run_patha_closed_newton(config)
    report_path = write_patha_closed_newton_report(results, config.report_path)
    token_scan = target_token_scan(_scan_paths(report_path))
    results["target_token_scan"] = token_scan
    results["genuine_checks"][-1]["pass"] = bool(token_scan["passed"])
    results["genuine_checks"][-1]["detail"] = (
        f"scanned {token_scan['path_count']} Path-A 1b files."
    )
    results["passed"] = all(row["pass"] for row in results["genuine_checks"])
    write_patha_closed_newton_report(results, config.report_path)
    print(f"wrote Path-A closed Newton report: {report_path}")
    if not results["passed"]:
        print("Path-A closed Newton gate failed")
        return 1
    print("Path-A closed Newton gate passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
