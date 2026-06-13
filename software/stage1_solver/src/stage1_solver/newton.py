"""Damped Newton-GMRES using torch JVPs."""

from __future__ import annotations

from dataclasses import dataclass, field
import time
from typing import Any, Callable

import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres
import torch

from .backend import jvp
from .config import NewtonConfig


@dataclass
class NewtonIteration:
    iteration: int
    residual_norm: float
    step_norm: float | None = None
    line_search_alpha: float | None = None
    gmres_info: int | None = None
    gmres_iterations: int | None = None
    preconditioner_type: str = "none"
    preconditioner_setup_seconds: float | None = None
    preconditioner_info: dict[str, Any] = field(default_factory=dict)


@dataclass
class NewtonResult:
    x: torch.Tensor
    converged: bool
    iterations: int
    initial_residual_norm: float
    final_residual_norm: float
    tolerance: float
    history: list[NewtonIteration] = field(default_factory=list)
    message: str = ""


@dataclass(frozen=True)
class PreconditionerBuildContext:
    residual_fn: Callable[[torch.Tensor], torch.Tensor]
    x: torch.Tensor
    rhs: np.ndarray
    iteration: int
    config: NewtonConfig


@dataclass(frozen=True)
class BuiltPreconditioner:
    operator: LinearOperator | None
    metadata: dict[str, Any] = field(default_factory=dict)


PreconditionerFactory = Callable[[PreconditionerBuildContext], BuiltPreconditioner]


def _norm(values: torch.Tensor, kind: str) -> float:
    if kind == "linf":
        return float(torch.max(torch.abs(values)).detach().cpu().item())
    if kind == "l2":
        return float(torch.linalg.norm(values).detach().cpu().item())
    raise ValueError(f"Unsupported residual norm {kind!r}")


def solve_newton_jvp(
    residual_fn: Callable[[torch.Tensor], torch.Tensor],
    x0: torch.Tensor,
    config: NewtonConfig,
    preconditioner_factory: PreconditionerFactory | None = None,
) -> NewtonResult:
    if config.linear_solver != "gmres_jvp":
        raise ValueError(f"Unsupported linear solver {config.linear_solver!r}")
    if config.line_search != "armijo":
        raise ValueError(f"Unsupported line search {config.line_search!r}")

    x = x0.detach().clone()
    r = residual_fn(x).detach()
    initial_norm = _norm(r, config.residual_norm)
    tolerance = max(config.residual_atol, config.residual_rtol * initial_norm)
    history = [NewtonIteration(iteration=0, residual_norm=initial_norm)]

    if initial_norm <= tolerance:
        return NewtonResult(
            x=x,
            converged=True,
            iterations=0,
            initial_residual_norm=initial_norm,
            final_residual_norm=initial_norm,
            tolerance=tolerance,
            history=history,
            message="initial residual met tolerance",
        )

    for iteration in range(1, config.max_newton_iters + 1):
        x_for_jvp = x.detach().clone()
        rhs = -r.detach().cpu().numpy().astype(np.float64)
        dim = rhs.size
        gmres_iterations = 0
        preconditioner_operator: LinearOperator | None = None
        preconditioner_metadata: dict[str, Any] = {"type": config.preconditioner.type}
        preconditioner_setup_seconds: float | None = None

        def matvec(vector_np: np.ndarray) -> np.ndarray:
            direction = torch.as_tensor(vector_np, dtype=x.dtype, device=x.device)
            jv = jvp(residual_fn, x_for_jvp, direction)
            return jv.detach().cpu().numpy().astype(np.float64)

        def callback(_residual_norm: float) -> None:
            nonlocal gmres_iterations
            gmres_iterations += 1

        linear_op = LinearOperator((dim, dim), matvec=matvec, dtype=np.float64)
        if config.preconditioner.type != "none":
            if config.preconditioner.side != "left":
                raise ValueError(
                    f"Unsupported preconditioner side {config.preconditioner.side!r}"
                )
            if preconditioner_factory is None:
                raise ValueError(
                    f"Preconditioner {config.preconditioner.type!r} requires a factory"
                )
            preconditioner_started = time.perf_counter()
            built_preconditioner = preconditioner_factory(
                PreconditionerBuildContext(
                    residual_fn=residual_fn,
                    x=x_for_jvp,
                    rhs=rhs,
                    iteration=iteration,
                    config=config,
                )
            )
            preconditioner_setup_seconds = time.perf_counter() - preconditioner_started
            preconditioner_operator = built_preconditioner.operator
            preconditioner_metadata = dict(built_preconditioner.metadata)
            preconditioner_metadata.setdefault("type", config.preconditioner.type)
            preconditioner_metadata.setdefault("side", config.preconditioner.side)
        elif preconditioner_factory is not None:
            raise ValueError("A preconditioner factory was provided while config type is 'none'")
        step_np, gmres_info = gmres(
            linear_op,
            rhs,
            M=preconditioner_operator,
            rtol=config.gmres_rtol,
            atol=config.gmres_atol,
            restart=config.gmres_restart,
            maxiter=config.gmres_maxiter,
            callback=callback,
            callback_type="pr_norm",
        )
        if gmres_info != 0:
            return NewtonResult(
                x=x,
                converged=False,
                iterations=iteration - 1,
                initial_residual_norm=initial_norm,
                final_residual_norm=_norm(r, config.residual_norm),
                tolerance=tolerance,
                history=history,
                message=f"GMRES failed with info={gmres_info}",
            )

        step = torch.as_tensor(step_np, dtype=x.dtype, device=x.device)
        step_norm = _norm(step, "l2")
        x_norm = max(_norm(x, "l2"), 1.0)
        if step_norm <= config.step_atol + config.step_rtol * x_norm:
            r_now = residual_fn(x).detach()
            final_norm = _norm(r_now, config.residual_norm)
            return NewtonResult(
                x=x,
                converged=final_norm <= tolerance,
                iterations=iteration - 1,
                initial_residual_norm=initial_norm,
                final_residual_norm=final_norm,
                tolerance=tolerance,
                history=history,
                message="step tolerance reached",
            )

        old_norm = _norm(r, config.residual_norm)
        alpha = 1.0
        accepted = False
        best_x = x
        best_r = r
        best_norm = old_norm
        accepted_alpha = None
        for _ in range(config.max_line_search_iters):
            candidate_x = x + alpha * step
            candidate_r = residual_fn(candidate_x).detach()
            candidate_norm = _norm(candidate_r, config.residual_norm)
            if np.isfinite(candidate_norm) and candidate_norm < best_norm:
                best_x = candidate_x.detach()
                best_r = candidate_r.detach()
                best_norm = candidate_norm
                accepted_alpha = alpha
            armijo_bound = (1.0 - config.line_search_c1 * alpha) * old_norm
            if np.isfinite(candidate_norm) and (
                candidate_norm <= armijo_bound or candidate_norm <= tolerance
            ):
                x = candidate_x.detach()
                r = candidate_r.detach()
                accepted = True
                accepted_alpha = alpha
                break
            alpha *= config.line_search_shrink

        if not accepted:
            if config.accept_best_line_search_decrease and best_norm < old_norm:
                x = best_x
                r = best_r
                accepted = True
            else:
                return NewtonResult(
                    x=x,
                    converged=False,
                    iterations=iteration - 1,
                    initial_residual_norm=initial_norm,
                    final_residual_norm=old_norm,
                    tolerance=tolerance,
                    history=history,
                    message="line search failed to reduce residual",
                )

        new_norm = _norm(r, config.residual_norm)
        history.append(
            NewtonIteration(
                iteration=iteration,
                residual_norm=new_norm,
                step_norm=step_norm,
                line_search_alpha=accepted_alpha,
                gmres_info=gmres_info,
                gmres_iterations=gmres_iterations,
                preconditioner_type=preconditioner_metadata.get(
                    "type",
                    config.preconditioner.type,
                ),
                preconditioner_setup_seconds=preconditioner_setup_seconds,
                preconditioner_info=preconditioner_metadata,
            )
        )
        if new_norm <= tolerance:
            return NewtonResult(
                x=x,
                converged=True,
                iterations=iteration,
                initial_residual_norm=initial_norm,
                final_residual_norm=new_norm,
                tolerance=tolerance,
                history=history,
                message="residual tolerance reached",
            )

    final_norm = _norm(r, config.residual_norm)
    return NewtonResult(
        x=x,
        converged=False,
        iterations=config.max_newton_iters,
        initial_residual_norm=initial_norm,
        final_residual_norm=final_norm,
        tolerance=tolerance,
        history=history,
        message="maximum Newton iterations reached",
    )


def finite_difference_jvp_check(
    residual_fn: Callable[[torch.Tensor], torch.Tensor],
    x: torch.Tensor,
    *,
    epsilon: float,
    seed: int,
) -> dict[str, float]:
    generator = torch.Generator(device=x.device)
    generator.manual_seed(seed)
    direction = torch.randn(x.shape, dtype=x.dtype, device=x.device, generator=generator)
    direction = direction / torch.linalg.norm(direction)
    jv = jvp(residual_fn, x.detach(), direction).detach()
    fd = (
        residual_fn((x + epsilon * direction).detach()).detach()
        - residual_fn((x - epsilon * direction).detach()).detach()
    ) / (2.0 * epsilon)
    diff = jv - fd
    abs_residual = float(torch.linalg.norm(diff).detach().cpu().item())
    denom = max(1.0, float(torch.linalg.norm(jv).detach().cpu().item()))
    rel_residual = abs_residual / denom
    fd_norm = float(torch.linalg.norm(fd).detach().cpu().item())
    jvp_norm = float(torch.linalg.norm(jv).detach().cpu().item())
    return {
        "epsilon": float(epsilon),
        "absolute_residual": abs_residual,
        "relative_residual": rel_residual,
        "fd_norm": fd_norm,
        "jvp_norm": jvp_norm,
    }
