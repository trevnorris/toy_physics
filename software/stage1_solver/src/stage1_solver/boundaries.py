"""Boundary operators for cell-centered finite-volume grids."""

from __future__ import annotations

from dataclasses import dataclass

import torch


@dataclass(frozen=True)
class BoundaryCondition:
    kind: str
    alpha: float
    beta: float
    gamma: float

    @staticmethod
    def dirichlet(value: float) -> "BoundaryCondition":
        return BoundaryCondition(kind="dirichlet", alpha=1.0, beta=0.0, gamma=float(value))

    @staticmethod
    def robin(alpha: float, beta: float, gamma: float = 0.0) -> "BoundaryCondition":
        if beta == 0.0:
            raise ValueError("Use dirichlet() for beta=0 boundary conditions")
        return BoundaryCondition(
            kind="robin",
            alpha=float(alpha),
            beta=float(beta),
            gamma=float(gamma),
        )

    @staticmethod
    def neumann(derivative: float) -> "BoundaryCondition":
        return BoundaryCondition.robin(alpha=0.0, beta=1.0, gamma=float(derivative))

    def to_dict(self) -> dict[str, float | str]:
        return {
            "kind": self.kind,
            "alpha": self.alpha,
            "beta": self.beta,
            "gamma": self.gamma,
        }


def ghost_value(interior: torch.Tensor, spacing: float, bc: BoundaryCondition) -> torch.Tensor:
    """Return the ghost-cell value implied by alpha*u + beta*du_dn = gamma.

    The ghost is placed one cell spacing outside the boundary face and the
    interior value is one half spacing inside. For both left and right
    boundaries, ``du_dn`` is the outward normal derivative.
    """

    h = float(spacing)
    if bc.kind == "dirichlet":
        return 2.0 * bc.gamma - interior
    if bc.kind == "robin":
        denom = 0.5 * bc.alpha + bc.beta / h
        if denom == 0.0:
            raise ZeroDivisionError("Degenerate Robin boundary operator")
        numer = torch.as_tensor(bc.gamma, dtype=interior.dtype, device=interior.device)
        numer = numer - (0.5 * bc.alpha - bc.beta / h) * interior
        return numer / denom
    raise ValueError(f"Unsupported boundary kind {bc.kind!r}")
