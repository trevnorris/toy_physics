from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path

import torch

from stage1_solver.backend import configure_backend
from stage1_solver.config import BackendConfig
from stage1_solver.patha_static_balance import (
    SSigmaSpec,
    registered_s_sigma_families,
    resolve_s_sigma,
)
from stage1_solver.patha_static_mms import (
    PathAStaticBalanceMMSConfig,
    run_patha_static_balance_mms,
)


def test_s_sigma_spec_registry_is_hashable_and_serializable() -> None:
    dtype = configure_backend(BackendConfig())
    spec = SSigmaSpec.smooth_positive_placeholder(w_min=0.0, w_max=1.0)
    round_tripped = SSigmaSpec.from_dict(json.loads(spec.to_json()))
    assert round_tripped == spec
    assert hash(round_tripped) == hash(spec)
    assert spec.digest() == round_tripped.digest()
    assert spec.family in registered_s_sigma_families()

    provider = resolve_s_sigma(spec)
    radius = torch.linspace(0.85, 1.25, 7, dtype=dtype)
    w = torch.linspace(0.0, 1.0, 7, dtype=dtype)
    for values in (
        provider.mu(radius, w),
        provider.T_w(radius, w),
        provider.T_Omega(radius, w),
        provider.U(radius, w),
    ):
        assert torch.all(torch.isfinite(values))
        assert torch.min(values).item() > 0.0
    for values in (
        provider.T_w_R(radius, w),
        provider.T_w_RR(radius, w),
        provider.U_R(radius, w),
        provider.U_RR(radius, w),
    ):
        assert torch.all(torch.isfinite(values))


def test_homogeneous_isotropic_hooke_family_values_and_derivatives() -> None:
    dtype = configure_backend(BackendConfig())
    tau = 2.5
    a = 1.25
    spec = {
        "family": "homogeneous_isotropic_hooke_v1",
        "parameters": {"tau": tau, "a": a, "w_min": 0.0, "w_max": 1.85},
    }
    assert "homogeneous_isotropic_hooke_v1" in registered_s_sigma_families()
    provider = resolve_s_sigma(spec)

    radius = torch.tensor([0.95, 1.25, 1.40], dtype=dtype)
    w = torch.tensor([0.20, 0.80, 1.30], dtype=dtype)
    t_omega = tau / (a * a)

    assert torch.allclose(provider.mu(radius, w), torch.full_like(radius, tau))
    assert torch.allclose(provider.T_w(radius, w), torch.full_like(radius, tau))
    assert torch.allclose(provider.T_Omega(radius, w), torch.full_like(radius, t_omega))
    assert torch.allclose(provider.U(radius, w), 0.5 * t_omega * (radius - a) ** 2)
    assert torch.allclose(provider.U_R(radius, w), t_omega * (radius - a))
    assert torch.allclose(provider.U_RR(radius, w), torch.full_like(radius, t_omega))

    eps = torch.as_tensor(1.0e-5, dtype=dtype)
    u_plus = provider.U(radius + eps, w)
    u = provider.U(radius, w)
    u_minus = provider.U(radius - eps, w)
    tw_plus = provider.T_w(radius + eps, w)
    tw = provider.T_w(radius, w)
    tw_minus = provider.T_w(radius - eps, w)

    assert torch.allclose(provider.U_R(radius, w), (u_plus - u_minus) / (2.0 * eps), atol=1.0e-10)
    assert torch.allclose(
        provider.U_RR(radius, w),
        (u_plus - 2.0 * u + u_minus) / (eps * eps),
        atol=1.0e-6,
    )
    assert torch.allclose(provider.T_w_R(radius, w), (tw_plus - tw_minus) / (2.0 * eps))
    assert torch.allclose(
        provider.T_w_RR(radius, w),
        (tw_plus - 2.0 * tw + tw_minus) / (eps * eps),
    )


def test_patha_static_balance_mms_is_dual_engine_second_order(tmp_path: Path) -> None:
    cfg = replace(
        PathAStaticBalanceMMSConfig(),
        run_root=str(tmp_path / "runs"),
        report_path=str(tmp_path / "patha_static_balance_report.md"),
    )
    result = run_patha_static_balance_mms(cfg)
    assert result["passed"] is True
    assert result["dual_engine"]["passed"] is True
    assert result["dual_engine"]["max_abs_diff"] <= cfg.dual_engine_abs_tol
    assert result["convergence"]["rows"][-1]["observed_order"] > cfg.min_observed_order
    assert result["target_token_scan"]["passed"] is True
    assert all(row["pass"] for row in result["genuine_checks"])
    assert all(
        not row["check"].endswith("not_a_physics_gate")
        for row in result["genuine_checks"]
    )
    assert {
        "dual_engine_forcing_agreement",
        "second_order_static_balance_mms",
        "flux_divergence_term_nonzero",
        "gradient_square_term_nonzero",
        "potential_gradient_term_nonzero",
        "target_token_scan_clean",
    } == {row["check"] for row in result["genuine_checks"]}
