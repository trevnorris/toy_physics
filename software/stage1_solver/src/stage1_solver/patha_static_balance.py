"""Path-A static throat-balance constitutive registry and FV operator."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping

import torch

from .boundaries import BoundaryCondition
from .config import stable_hash, stable_json
from .grid import WallGrid


ParameterItems = tuple[tuple[str, float], ...]


def _sorted_parameters(parameters: Mapping[str, float]) -> ParameterItems:
    return tuple(sorted((str(key), float(value)) for key, value in parameters.items()))


@dataclass(frozen=True)
class SSigmaSpec:
    """Serializable, hashable selector for a registered ``S_Sigma`` family."""

    family: str
    parameters: ParameterItems

    @staticmethod
    def smooth_positive_placeholder(*, w_min: float = 0.0, w_max: float = 1.0) -> "SSigmaSpec":
        """Smooth positive placeholder family for target-blind engineering runs."""

        return SSigmaSpec(
            family="smooth_positive_placeholder_v1",
            parameters=_sorted_parameters(
                {
                    "w_min": w_min,
                    "w_max": w_max,
                    "mu_base": 1.0,
                    "mu_r2": 0.04,
                    "mu_w_amp": 0.07,
                    "tw_base": 1.2,
                    "tw_r2": 0.06,
                    "tw_w_amp": 0.09,
                    "tomega_base": 0.85,
                    "tomega_r2": 0.03,
                    "tomega_w_amp": 0.04,
                    "u_base": 0.35,
                    "u_r2": 0.12,
                    "u_r4": 0.015,
                    "u_w_amp": 0.03,
                }
            ),
        )

    @staticmethod
    def patha_static_mms(*, w_min: float = 0.0, w_max: float = 1.4) -> "SSigmaSpec":
        """Manufactured family used only by the chunk-1a static-balance MMS."""

        return SSigmaSpec(
            family="patha_static_mms_v1",
            parameters=_sorted_parameters(
                {
                    "w_min": w_min,
                    "w_max": w_max,
                    "mu_base": 1.05,
                    "mu_r2": 0.025,
                    "mu_w_amp": 0.035,
                    "tw_base": 1.32,
                    "tw_w_sin": 0.18,
                    "tw_r1": 0.11,
                    "tw_r2": 0.04,
                    "tw_rw": 0.03,
                    "tomega_base": 0.82,
                    "tomega_r1": 0.04,
                    "tomega_w_cos": 0.025,
                    "u_base": 0.42,
                    "u_r2": 0.20,
                    "u_r3": 0.03,
                    "u_rw": 0.05,
                }
            ),
        )

    @staticmethod
    def homogeneous_isotropic_hooke(
        *,
        tau: float,
        a: float = 1.0,
        w_min: float = 0.0,
        w_max: float = 37.0 / 20.0,
    ) -> "SSigmaSpec":
        """Frozen Path-A GATE-A homogeneous isotropic Hookean wall."""

        return SSigmaSpec(
            family="homogeneous_isotropic_hooke_v1",
            parameters=_sorted_parameters(
                {
                    "tau": tau,
                    "a": a,
                    "w_min": w_min,
                    "w_max": w_max,
                }
            ),
        )

    @staticmethod
    def from_dict(data: Mapping[str, Any]) -> "SSigmaSpec":
        if set(data) != {"family", "parameters"}:
            raise ValueError("S_Sigma spec requires exactly family and parameters")
        raw_parameters = data["parameters"]
        if isinstance(raw_parameters, Mapping):
            parameters = _sorted_parameters(raw_parameters)
        else:
            parameters = _sorted_parameters(dict(raw_parameters))
        spec = SSigmaSpec(family=str(data["family"]), parameters=parameters)
        if spec.family not in _S_SIGMA_REGISTRY:
            raise ValueError(f"Unregistered S_Sigma family {spec.family!r}")
        return spec

    def to_dict(self) -> dict[str, Any]:
        return {"family": self.family, "parameters": dict(self.parameters)}

    def to_json(self) -> str:
        return stable_json(self.to_dict())

    def digest(self) -> str:
        return stable_hash(self.to_dict())


class SSigmaProvider:
    """Callable tensor provider resolved from a serializable ``SSigmaSpec``."""

    def __init__(self, spec: SSigmaSpec):
        self.spec = spec
        self._parameters = dict(spec.parameters)

    def _p(self, name: str) -> float:
        return self._parameters[name]

    def _x(self, w: torch.Tensor) -> torch.Tensor:
        length = self._p("w_max") - self._p("w_min")
        if length <= 0.0:
            raise ValueError("S_Sigma provider requires w_max > w_min")
        return (w - self._p("w_min")) / length

    def mu(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        raise NotImplementedError

    def T_w(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        raise NotImplementedError

    def T_w_R(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        raise NotImplementedError

    def T_w_RR(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        raise NotImplementedError

    def T_Omega(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        raise NotImplementedError

    def U(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        raise NotImplementedError

    def U_R(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        raise NotImplementedError

    def U_RR(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        raise NotImplementedError


class _SmoothPositivePlaceholderProvider(SSigmaProvider):
    def mu(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        x = self._x(w)
        return self._p("mu_base") + self._p("mu_r2") * R**2 + self._p("mu_w_amp") * (
            1.0 + torch.sin(torch.pi * x) ** 2
        )

    def T_w(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        x = self._x(w)
        return self._p("tw_base") + self._p("tw_r2") * R**2 + self._p("tw_w_amp") * (
            1.0 + torch.sin(2.0 * torch.pi * x) ** 2
        )

    def T_w_R(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        return 2.0 * self._p("tw_r2") * R + 0.0 * w

    def T_w_RR(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        return torch.zeros_like(R + w) + 2.0 * self._p("tw_r2")

    def T_Omega(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        x = self._x(w)
        return (
            self._p("tomega_base")
            + self._p("tomega_r2") * R**2
            + self._p("tomega_w_amp") * (1.0 + torch.cos(2.0 * torch.pi * x) ** 2)
        )

    def U(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        x = self._x(w)
        return (
            self._p("u_base")
            + self._p("u_r2") * R**2
            + self._p("u_r4") * R**4
            + self._p("u_w_amp") * (1.0 + torch.sin(torch.pi * x) ** 2)
        )

    def U_R(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        return 2.0 * self._p("u_r2") * R + 4.0 * self._p("u_r4") * R**3 + 0.0 * w

    def U_RR(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        return 2.0 * self._p("u_r2") + 12.0 * self._p("u_r4") * R**2 + 0.0 * w


class _PathAStaticMMSProvider(SSigmaProvider):
    def mu(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        x = self._x(w)
        return self._p("mu_base") + self._p("mu_r2") * R**2 + self._p("mu_w_amp") * (
            1.0 + torch.sin(torch.pi * x) ** 2
        )

    def T_w(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        x = self._x(w)
        return (
            self._p("tw_base")
            + self._p("tw_w_sin") * torch.sin(2.0 * torch.pi * x)
            + self._p("tw_r1") * R
            + self._p("tw_r2") * R**2
            + self._p("tw_rw") * R * torch.cos(torch.pi * x)
        )

    def T_w_R(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        x = self._x(w)
        return self._p("tw_r1") + 2.0 * self._p("tw_r2") * R + self._p(
            "tw_rw"
        ) * torch.cos(torch.pi * x)

    def T_w_RR(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        return torch.zeros_like(R + w) + 2.0 * self._p("tw_r2")

    def T_Omega(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        x = self._x(w)
        return (
            self._p("tomega_base")
            + self._p("tomega_r1") * R
            + self._p("tomega_w_cos") * (1.0 + torch.cos(2.0 * torch.pi * x))
        )

    def U(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        x = self._x(w)
        return (
            self._p("u_base")
            + self._p("u_r2") * R**2
            + self._p("u_r3") * R**3
            + self._p("u_rw") * R * torch.sin(torch.pi * x)
        )

    def U_R(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        x = self._x(w)
        return (
            2.0 * self._p("u_r2") * R
            + 3.0 * self._p("u_r3") * R**2
            + self._p("u_rw") * torch.sin(torch.pi * x)
        )

    def U_RR(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        return 2.0 * self._p("u_r2") + 6.0 * self._p("u_r3") * R + 0.0 * w


class _HomogeneousIsotropicHookeProvider(SSigmaProvider):
    def _tau_tensor(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        tau = self._p("tau")
        if tau <= 0.0:
            raise ValueError("homogeneous_isotropic_hooke_v1 requires tau > 0")
        if self._p("a") <= 0.0:
            raise ValueError("homogeneous_isotropic_hooke_v1 requires a > 0")
        if self._p("w_max") <= self._p("w_min"):
            raise ValueError("homogeneous_isotropic_hooke_v1 requires w_max > w_min")
        return torch.zeros_like(R + w) + tau

    def _tomega_tensor(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        return self._tau_tensor(R, w) / (self._p("a") ** 2)

    def mu(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        return self._tau_tensor(R, w)

    def T_w(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        return self._tau_tensor(R, w)

    def T_w_R(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        self._tau_tensor(R, w)
        return torch.zeros_like(R + w)

    def T_w_RR(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        self._tau_tensor(R, w)
        return torch.zeros_like(R + w)

    def T_Omega(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        return self._tomega_tensor(R, w)

    def U(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        return 0.5 * self._tomega_tensor(R, w) * (R - self._p("a")) ** 2

    def U_R(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        return self._tomega_tensor(R, w) * (R - self._p("a"))

    def U_RR(self, R: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        return self._tomega_tensor(R, w)


_S_SIGMA_REGISTRY: dict[str, Callable[[SSigmaSpec], SSigmaProvider]] = {
    "smooth_positive_placeholder_v1": _SmoothPositivePlaceholderProvider,
    "patha_static_mms_v1": _PathAStaticMMSProvider,
    "homogeneous_isotropic_hooke_v1": _HomogeneousIsotropicHookeProvider,
}


def registered_s_sigma_families() -> tuple[str, ...]:
    return tuple(sorted(_S_SIGMA_REGISTRY))


def resolve_s_sigma(spec: SSigmaSpec | Mapping[str, Any]) -> SSigmaProvider:
    resolved_spec = spec if isinstance(spec, SSigmaSpec) else SSigmaSpec.from_dict(spec)
    try:
        provider_type = _S_SIGMA_REGISTRY[resolved_spec.family]
    except KeyError as exc:
        raise ValueError(f"Unregistered S_Sigma family {resolved_spec.family!r}") from exc
    return provider_type(resolved_spec)


@dataclass(frozen=True)
class StaticBalanceTerms:
    flux_divergence: torch.Tensor
    gradient_square: torch.Tensor
    potential_gradient: torch.Tensor
    source: torch.Tensor
    lhs: torch.Tensor
    residual: torch.Tensor
    face_fluxes: torch.Tensor
    center_gradient: torch.Tensor


def wall_center_gradient(values: torch.Tensor, grid: WallGrid) -> torch.Tensor:
    """Second-order gradient at wall cell centers."""

    if values.ndim != 1 or values.shape[0] != grid.spec.nw:
        raise ValueError("wall_center_gradient expects a one-dimensional wall field")
    if grid.spec.nw < 3:
        raise ValueError("Need at least three wall cells")
    grad = torch.empty_like(values)
    grad[1:-1] = (values[2:] - values[:-2]) / (2.0 * grid.dw)
    grad[0] = (-3.0 * values[0] + 4.0 * values[1] - values[2]) / (2.0 * grid.dw)
    grad[-1] = (3.0 * values[-1] - 4.0 * values[-2] + values[-3]) / (2.0 * grid.dw)
    return grad


def static_balance_terms(
    values: torch.Tensor,
    grid: WallGrid,
    *,
    s_sigma: SSigmaSpec | SSigmaProvider | Mapping[str, Any],
    source: torch.Tensor,
    lower_bc: BoundaryCondition,
    upper_bc: BoundaryCondition,
) -> StaticBalanceTerms:
    """Evaluate the nonlinear static throat-balance residual on a wall grid."""

    if values.ndim != 1 or values.shape[0] != grid.spec.nw:
        raise ValueError("static_balance_terms expects a one-dimensional wall field")
    if source.shape != values.shape:
        raise ValueError("static_balance_terms source must match values")
    if lower_bc.kind != "dirichlet":
        raise ValueError("static balance requires a Dirichlet mouth boundary")
    upper_zero_traction = (
        upper_bc.kind == "robin"
        and upper_bc.alpha == 0.0
        and upper_bc.beta != 0.0
        and upper_bc.gamma == 0.0
    )
    if upper_bc.kind != "dirichlet" and not upper_zero_traction:
        raise ValueError(
            "static balance supports a Dirichlet exit or natural zero-traction exit"
        )
    provider = s_sigma if isinstance(s_sigma, SSigmaProvider) else resolve_s_sigma(s_sigma)

    lower_value = torch.as_tensor(lower_bc.gamma, dtype=values.dtype, device=values.device)
    upper_value = torch.as_tensor(upper_bc.gamma, dtype=values.dtype, device=values.device)

    face_values = torch.empty(grid.spec.nw + 1, dtype=values.dtype, device=values.device)
    face_values[0] = lower_value
    face_values[1:-1] = 0.5 * (values[1:] + values[:-1])
    face_values[-1] = values[-1] if upper_zero_traction else upper_value

    t_w_faces = provider.T_w(face_values, grid.w_faces)
    fluxes = torch.zeros(grid.spec.nw + 1, dtype=values.dtype, device=values.device)
    fluxes[0] = t_w_faces[0] * (
        (-46.0 / 15.0) * lower_value
        + (15.0 / 4.0) * values[0]
        - (5.0 / 6.0) * values[1]
        + (3.0 / 20.0) * values[2]
    ) / grid.dw
    fluxes[1:-1] = t_w_faces[1:-1] * (values[1:] - values[:-1]) / grid.dw
    if upper_zero_traction:
        fluxes[-1] = torch.zeros((), dtype=values.dtype, device=values.device)
    else:
        fluxes[-1] = t_w_faces[-1] * (
            (46.0 / 15.0) * upper_value
            - (15.0 / 4.0) * values[-1]
            + (5.0 / 6.0) * values[-2]
            - (3.0 / 20.0) * values[-3]
        ) / grid.dw

    center_gradient = wall_center_gradient(values, grid)
    flux_divergence = -(fluxes[1:] - fluxes[:-1]) / grid.cell_widths
    gradient_square = (
        0.5 * provider.T_w_R(values, grid.w_centers) * center_gradient * center_gradient
    )
    potential_gradient = provider.U_R(values, grid.w_centers)
    lhs = flux_divergence + gradient_square + potential_gradient
    residual = lhs - source
    return StaticBalanceTerms(
        flux_divergence=flux_divergence,
        gradient_square=gradient_square,
        potential_gradient=potential_gradient,
        source=source,
        lhs=lhs,
        residual=residual,
        face_fluxes=fluxes,
        center_gradient=center_gradient,
    )


def static_balance_operator(
    values: torch.Tensor,
    grid: WallGrid,
    *,
    s_sigma: SSigmaSpec | SSigmaProvider | Mapping[str, Any],
    source: torch.Tensor,
    lower_bc: BoundaryCondition,
    upper_bc: BoundaryCondition,
) -> torch.Tensor:
    """Return ``-div(T_w grad R0) + 0.5 T_w_R |grad R0|^2 + U_R - S``."""

    return static_balance_terms(
        values,
        grid,
        s_sigma=s_sigma,
        source=source,
        lower_bc=lower_bc,
        upper_bc=upper_bc,
    ).residual
