from __future__ import annotations

import copy
import importlib.util
import math
import numbers
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import sympy as sp

from stage1_solver import patha_extraction as pe
from stage1_solver.patha_static_balance import SSigmaSpec, resolve_s_sigma


def _load_script_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _numeric_diff(left: Any, right: Any) -> tuple[float, float]:
    max_abs = 0.0
    max_rel = 0.0

    def visit(a: Any, b: Any, path: str = "$") -> None:
        nonlocal max_abs, max_rel
        if isinstance(a, dict) or isinstance(b, dict):
            assert isinstance(a, dict) and isinstance(b, dict), f"type mismatch at {path}: {type(a)} vs {type(b)}"
            left_keys = set(a)
            right_keys = set(b)
            assert left_keys == right_keys, (
                f"dict key mismatch at {path}: "
                f"left_only={sorted(left_keys - right_keys)!r}, "
                f"right_only={sorted(right_keys - left_keys)!r}"
            )
            for key in sorted(left_keys):
                visit(a[key], b[key], f"{path}.{key}")
            return
        if isinstance(a, list) or isinstance(b, list):
            assert isinstance(a, list) and isinstance(b, list), f"type mismatch at {path}: {type(a)} vs {type(b)}"
            assert len(a) == len(b), f"list length mismatch at {path}: {len(a)} != {len(b)}"
            for idx, (item_a, item_b) in enumerate(zip(a, b)):
                visit(item_a, item_b, f"{path}[{idx}]")
            return
        if isinstance(a, str) or isinstance(b, str):
            assert isinstance(a, str) and isinstance(b, str), f"type mismatch at {path}: {type(a)} vs {type(b)}"
            assert a == b, f"string mismatch at {path}: {a!r} != {b!r}"
            return
        if isinstance(a, (bool, np.bool_)) or isinstance(b, (bool, np.bool_)):
            assert bool(a) == bool(b), f"bool mismatch at {path}: {a!r} != {b!r}"
            return
        if isinstance(a, numbers.Real) and isinstance(b, numbers.Real):
            fa = float(a)
            fb = float(b)
            if math.isfinite(fa) and math.isfinite(fb):
                diff = abs(fa - fb)
                max_abs = max(max_abs, diff)
                max_rel = max(max_rel, diff / max(abs(fb), 1.0e-300))
            elif math.isnan(fa) and math.isnan(fb):
                return
            else:
                raise AssertionError(f"nonmatching numeric values at {path}: {a!r} vs {b!r}")
            return
        assert a == b, f"value mismatch at {path}: {a!r} != {b!r}"

    visit(left, right)
    return max_abs, max_rel


V21_PASS_FLAG_KEYS = (
    "open_gate_pass",
    "stability_gate_pass",
    "isotropic_D0_pass",
    "isotropic_N0_pass",
    "one_pole_pass",
    "normalization_pass",
    "constant_prefactor_P2_pass",
    "constant_prefactor_P4_pass",
    "tail_transport_pass",
    "target_packet_pass",
)
V22_PASS_FLAG_KEYS = (
    "open_gate_pass",
    "stability_gate_pass",
    "isotropic_D0_pass",
    "isotropic_N0_pass",
    "axisymmetric_P0_pattern_pass",
    "one_pole_pass",
    "normalization_pass",
    "constant_prefactor_P2_pass",
    "constant_prefactor_P4_pass",
    "tail_transport_pass",
    "target_packet_pass",
)
PACKET_KEYS = (
    "name",
    "input_hash",
    "open_gate",
    "lane_packets",
    "grouped",
    "target_values",
    "residuals",
    "stability",
    "pass_flags",
)
V22_LANE_OMIT_KEYS = {
    "constant_prefactor_N2_residual",
    "constant_prefactor_N4_residual",
}


def _project_lane_packets(
    lane_packets: dict[str, dict[str, Any]],
    *,
    omit_keys: set[str] | frozenset[str] = frozenset(),
) -> dict[str, dict[str, Any]]:
    return {
        lane: {key: value for key, value in packet.items() if key not in omit_keys}
        for lane, packet in lane_packets.items()
    }


def _project_packet(
    packet: dict[str, Any],
    *,
    pass_flag_keys: tuple[str, ...],
    allowed_extra_pass_flags: set[str] | frozenset[str] = frozenset(),
    omit_lane_keys: set[str] | frozenset[str] = frozenset(),
) -> dict[str, Any]:
    assert set(packet) == set(PACKET_KEYS)
    missing_pass_flags = set(pass_flag_keys) - set(packet["pass_flags"])
    assert not missing_pass_flags, f"missing pass flags: {sorted(missing_pass_flags)!r}"
    unexpected_pass_flags = set(packet["pass_flags"]) - set(pass_flag_keys) - set(allowed_extra_pass_flags)
    assert not unexpected_pass_flags, f"unexpected pass flags: {sorted(unexpected_pass_flags)!r}"
    return {
        "name": packet["name"],
        "input_hash": packet["input_hash"],
        "open_gate": packet["open_gate"],
        "lane_packets": _project_lane_packets(packet["lane_packets"], omit_keys=omit_lane_keys),
        "grouped": packet["grouped"],
        "target_values": packet["target_values"],
        "residuals": packet["residuals"],
        "stability": packet["stability"],
        "pass_flags": {key: packet["pass_flags"][key] for key in pass_flag_keys},
    }


def _v21_packet_view(packet: dict[str, Any]) -> dict[str, Any]:
    return _project_packet(
        packet,
        pass_flag_keys=V21_PASS_FLAG_KEYS,
        allowed_extra_pass_flags={"axisymmetric_P0_pattern_pass"},
    )


def _v22_packet_view(packet: dict[str, Any]) -> dict[str, Any]:
    return _project_packet(
        packet,
        pass_flag_keys=V22_PASS_FLAG_KEYS,
        omit_lane_keys=V22_LANE_OMIT_KEYS,
    )


def test_numeric_diff_rejects_mutated_module_output_structure_and_strings() -> None:
    v21 = _load_script_module(
        "stage_v2_21_branch_extraction_fixture_mutation_oracle",
        _repo_root() / "research/pde_audit/scripts/stage_v2_21_branch_extraction_fixture.py",
    )
    packet = _v21_packet_view(pe.extract_branch(v21.default_manifest()["branches"][1], tol=1.0e-9))

    extra_key = copy.deepcopy(packet)
    extra_key["open_gate"]["unexpected_gate"] = True
    with pytest.raises(AssertionError, match="dict key mismatch"):
        _numeric_diff(extra_key, packet)

    missing_key = copy.deepcopy(packet)
    del missing_key["open_gate"]["hard_cap_forbidden"]
    with pytest.raises(AssertionError, match="dict key mismatch"):
        _numeric_diff(missing_key, packet)

    changed_string = copy.deepcopy(packet)
    changed_string["name"] = f"{changed_string['name']}_mutated"
    with pytest.raises(AssertionError, match="string mismatch"):
        _numeric_diff(changed_string, packet)

    changed_list_length = copy.deepcopy(packet)
    changed_list_length["lane_packets"]["20"]["port_diagnostics"].append(
        copy.deepcopy(changed_list_length["lane_packets"]["20"]["port_diagnostics"][0])
    )
    with pytest.raises(AssertionError, match="list length mismatch"):
        _numeric_diff(changed_list_length, packet)


def _observed_orders(errors: list[float], spacings: list[float]) -> list[float]:
    orders = []
    for left, right, h_left, h_right in zip(errors, errors[1:], spacings, spacings[1:]):
        orders.append(math.log(left / right) / math.log(h_left / h_right))
    return orders


def _piecewise_linear_mode_error(w: np.ndarray, chi: np.ndarray, *, L: float) -> float:
    points, weights = np.polynomial.legendre.leggauss(6)
    total = 0.0
    for left in range(w.size - 1):
        x0 = w[left]
        x1 = w[left + 1]
        half = 0.5 * (x1 - x0)
        center = 0.5 * (x0 + x1)
        samples = center + half * points
        t = (samples - x0) / (x1 - x0)
        interp = (1.0 - t) * chi[left] + t * chi[left + 1]
        exact = pe.analytic_hooke_l2_mode(samples, L=L)
        total += half * float(np.sum(weights * (interp - exact) ** 2))
    return math.sqrt(total)


def test_l2_hooke_eigenproblem_matches_analytic_oracle_and_ignores_r0() -> None:
    tau = 1.7
    a = 1.0
    L = 37.0 / 20.0
    kappa_hat = pe.hooke_kappa_hat(a=a, L=L)
    spec = SSigmaSpec.homogeneous_isotropic_hooke(tau=tau, a=a, w_min=0.0, w_max=L)
    provider = resolve_s_sigma(spec)

    kappa_errors: list[float] = []
    chi_errors: list[float] = []
    spacings: list[float] = []
    for elements in (32, 64, 128, 256):
        w = np.linspace(0.0, L, elements + 1)
        R0 = a + 0.04 * (w / L) * (1.0 - w / L) + 0.01 * np.sin(2.0 * np.pi * w / L)
        result = pe.solve_l2_wall_eigenproblem(w, R0, provider)
        unknown = result.chi[1:]
        field_error = _piecewise_linear_mode_error(w, result.chi, L=L)
        kappa_errors.append(abs(result.K / tau - kappa_hat))
        chi_errors.append(field_error)
        spacings.append(L / elements)
        assert abs(float(unknown @ result.matrices.W @ unknown) - 1.0) <= 5.0e-13
        assert result.orientation_integral > 0.0
        assert abs(result.M - tau) <= 2.0e-12

    kappa_orders = _observed_orders(kappa_errors, spacings)
    chi_orders = _observed_orders(chi_errors, spacings)
    assert kappa_orders[-1] >= 1.95
    assert chi_orders[-1] >= 1.95
    assert kappa_errors[-1] <= 4.0e-5
    assert chi_errors[-1] <= 8.0e-6

    w = np.linspace(0.0, L, 129)
    flat_R0 = np.full_like(w, a)
    perturbed_R0 = a + 0.15 * np.sin(np.pi * w / L) + 0.02 * (w / L) ** 2
    flat = pe.solve_l2_wall_eigenproblem(w, flat_R0, spec)
    perturbed = pe.solve_l2_wall_eigenproblem(w, perturbed_R0, spec)
    assert abs(flat.K - perturbed.K) <= 1.0e-12


def test_v2_21_observable_fixture_reproduction() -> None:
    v21 = _load_script_module(
        "stage_v2_21_branch_extraction_fixture_oracle",
        _repo_root() / "research/pde_audit/scripts/stage_v2_21_branch_extraction_fixture.py",
    )
    max_abs = 0.0
    max_rel = 0.0
    for branch in v21.default_manifest()["branches"]:
        oracle = v21.extract_branch(branch, tol=1.0e-9)
        ours = pe.extract_branch(branch, tol=1.0e-9)
        abs_diff, rel_diff = _numeric_diff(_v21_packet_view(ours), _v21_packet_view(oracle))
        max_abs = max(max_abs, abs_diff)
        max_rel = max(max_rel, rel_diff)
    assert max_abs <= 1.0e-12
    assert max_rel <= 1.0e-12


def test_v2_22a_profile_adapter_fixture_reproduction() -> None:
    v22 = _load_script_module(
        "stage_v2_22a_profile_to_coefficient_adapter_oracle",
        _repo_root() / "research/pde_audit/scripts/stage_v2_22a_profile_to_coefficient_adapter.py",
    )
    profile_manifest = v22.default_profile_manifest()
    oracle_manifest, oracle_diag = v22.adapt_manifest(profile_manifest)
    ours_manifest, ours_diag = pe.adapt_profile_manifest(profile_manifest)

    abs_diff, rel_diff = _numeric_diff(ours_manifest, oracle_manifest)
    diag_abs, diag_rel = _numeric_diff(ours_diag, oracle_diag)
    max_abs = max(abs_diff, diag_abs)
    max_rel = max(rel_diff, diag_rel)

    for ours_branch, oracle_branch in zip(ours_manifest["branches"], oracle_manifest["branches"]):
        ours_packet = pe.extract_branch(ours_branch, tol=1.0e-9)
        oracle_packet = v22.extract_branch(oracle_branch, tol=1.0e-9)
        abs_packet, rel_packet = _numeric_diff(_v22_packet_view(ours_packet), _v22_packet_view(oracle_packet))
        max_abs = max(max_abs, abs_packet)
        max_rel = max(max_rel, rel_packet)

    assert max_abs <= 1.0e-12
    assert max_rel <= 1.0e-12


_SYMBOLIC_PREFACTOR_FUNCS: tuple[Any, Any, Any] | None = None


def _symbolic_prefactor_functions() -> tuple[Any, Any, Any]:
    global _SYMBOLIC_PREFACTOR_FUNCS
    if _SYMBOLIC_PREFACTOR_FUNCS is None:
        D0, D2, D4, N0, N2, N4, omega = sp.symbols("D0 D2 D4 N0 N2 N4 omega")
        D = D0 + D2 * omega**2 + D4 * omega**4
        N = N0 + N2 * omega**2 + N4 * omega**4
        prefactor = sp.series(D0 * N / D**2, omega, 0, 6).removeO().expand()
        args = (D0, D2, D4, N0, N2, N4)
        _SYMBOLIC_PREFACTOR_FUNCS = tuple(
            sp.lambdify(args, sp.simplify(prefactor.coeff(omega, power)), "math")
            for power in (0, 2, 4)
        )
    return _SYMBOLIC_PREFACTOR_FUNCS


def _symbolic_prefactor_values(
    D0: float,
    D2: float,
    D4: float,
    N0: float,
    N2: float,
    N4: float,
) -> tuple[float, float, float]:
    return tuple(float(fn(D0, D2, D4, N0, N2, N4)) for fn in _symbolic_prefactor_functions())


def _weighted_bar(values: dict[str, float]) -> float:
    return (values["20"] + 2.0 * values["21"] + 2.0 * values["22"]) / 5.0


def _assert_close(left: float, right: float, *, label: str) -> None:
    assert math.isclose(float(left), float(right), rel_tol=1.0e-12, abs_tol=1.0e-12), (
        f"{label}: {left!r} != {right!r}"
    )


def _assert_symbolic_observable_audit(branch: dict[str, Any]) -> None:
    packet = pe.extract_branch(branch, tol=1.0e-9)
    independent: dict[str, dict[str, float]] = {
        key: {} for key in ("D0", "D2", "D4", "P0", "P2", "P4", "R_pole", "A", "C")
    }

    for lane, lane_packet in packet["lane_packets"].items():
        D0 = lane_packet["K"] - lane_packet["B0"] - lane_packet["Z0"]
        D2 = -(lane_packet["M"] + lane_packet["B2"] + lane_packet["Z2"])
        D4 = -(lane_packet["B4"] + lane_packet["Z4"])
        P0, P2, P4 = _symbolic_prefactor_values(
            D0,
            D2,
            D4,
            lane_packet["N0"],
            lane_packet["N2"],
            lane_packet["N4"],
        )
        A = lane_packet["M"] + lane_packet["B2"] + lane_packet["Z2"]
        C = lane_packet["B4"] + lane_packet["Z4"]
        R_pole = D0 * C - 3.0 * A**2

        for key, expected in (
            ("D0", D0),
            ("D2", D2),
            ("D4", D4),
            ("P0", P0),
            ("P2", P2),
            ("P4", P4),
            ("R_pole", R_pole),
        ):
            _assert_close(lane_packet[key], expected, label=f"{packet['name']} lane {lane} {key}")
            independent[key][lane] = expected
        independent["A"][lane] = A
        independent["C"][lane] = C

    for key in ("D0", "D2", "D4", "P0", "P2", "P4", "R_pole"):
        _assert_close(packet["grouped"][key]["bar"], _weighted_bar(independent[key]), label=f"{packet['name']} {key}_bar")

    constants = branch.get("target", {}).get("constants", {})
    G = float(constants.get("G", 1.0))
    c_s = float(constants.get("c_s", 1.0))
    c = float(constants.get("c", 1.0))
    a = float(constants.get("a", 1.0))
    mhat0 = float(constants.get("mhat0", 1.0))
    S_port = float(constants.get("S_port", 1.0))
    P0_target = 54.0 * G * c_s**5 / (5.0 * a**5 * c**5)
    D0_bar = _weighted_bar(independent["D0"])
    A_bar = _weighted_bar(independent["A"])
    C_bar = _weighted_bar(independent["C"])
    _assert_close(packet["residuals"]["R_norm"], mhat0**2 * S_port * _weighted_bar(independent["P0"]) - P0_target, label=f"{packet['name']} R_norm")
    _assert_close(packet["residuals"]["R_pole"], D0_bar * C_bar - 3.0 * A_bar**2, label=f"{packet['name']} R_pole")
    _assert_close(packet["residuals"]["R_P2"], _weighted_bar(independent["P2"]), label=f"{packet['name']} R_P2")
    _assert_close(packet["residuals"]["R_P4"], _weighted_bar(independent["P4"]), label=f"{packet['name']} R_P4")


def test_independent_symbolic_observable_audit_matches_fixture_outputs() -> None:
    v21 = _load_script_module(
        "stage_v2_21_branch_extraction_fixture_symbolic_oracle",
        _repo_root() / "research/pde_audit/scripts/stage_v2_21_branch_extraction_fixture.py",
    )
    v22 = _load_script_module(
        "stage_v2_22a_profile_to_coefficient_adapter_symbolic_oracle",
        _repo_root() / "research/pde_audit/scripts/stage_v2_22a_profile_to_coefficient_adapter.py",
    )
    assert v21.run_symbolic_audit()["checks_failed"] == 0
    assert v22.run_symbolic_audit()["checks_failed"] == 0

    v22_manifest, _ = v22.adapt_manifest(v22.default_profile_manifest())
    for branch in [*v21.default_manifest()["branches"], *v22_manifest["branches"]]:
        _assert_symbolic_observable_audit(branch)


def _mms_exact_coefficients() -> dict[str, float]:
    sqrt3 = math.sqrt(3.0)
    I_eta_phi = sqrt3 * (1.0 / 2.0 + 1.0 / 4.0)
    I_eta_u = sqrt3 * (1.0 / 2.0 + 1.0 / 3.0)
    I_eta_w = sqrt3 * (1.0 / 2.0 + 1.0 / 5.0)
    I_u_w = 1.0 + 1.0 / 2.0 + 1.0 / 4.0 + 1.0 / 5.0

    lam_B = 0.37
    varpi = 2.3
    lam_U = 0.23
    lam_W = 0.31
    lam_R = 0.17
    OU = 2.7
    OW = 3.4
    c = lam_B * I_eta_phi
    gU = lam_U * I_eta_u
    gW = lam_W * I_eta_w
    R = lam_R * I_u_w
    Delta = OU**2 * OW**2 - R**2
    S = OU**2 + OW**2
    Q = gU**2 * OW**2 + 2.0 * gU * gW * R + gW**2 * OU**2
    H = gU**2 + gW**2
    P = OU**2 * gW + R * gU
    return {
        "B0": c**2 / varpi**2,
        "B2": c**2 / varpi**4,
        "B4": c**2 / varpi**6,
        "Z0": Q / Delta,
        "Z2": (Q * S - H * Delta) / Delta**2,
        "Z4": (Q * (S**2 - Delta) - S * H * Delta) / Delta**3,
        "N0": P**2 / Delta**2,
        "N2": 2.0 * P * (P * S - Delta * gW) / Delta**3,
        "N4": (
            Delta**2 * gW**2
            - 2.0 * Delta * P**2
            - 4.0 * Delta * P * S * gW
            + 3.0 * P**2 * S**2
        )
        / Delta**4,
    }


def test_manufactured_sampled_field_to_coefficient_round_trip_is_second_order() -> None:
    exact = _mms_exact_coefficients()
    errors: list[float] = []
    spacings: list[float] = []
    for nodes in (17, 33, 65, 129):
        w = np.linspace(0.0, 1.0, nodes)
        chi = math.sqrt(3.0) * w
        phi = 1.0 + w**2
        u_profile = 1.0 + w
        w_profile = 1.0 + w**3
        packet = pe.coefficients_from_sampled_fields(
            w,
            chi,
            K=1.8,
            M=0.6,
            bdg_modes=[
                {
                    "name": "poly_bdg",
                    "lambda_B": 0.37,
                    "varpi": 2.3,
                    "profile": phi,
                }
            ],
            mixed_ports=[
                {
                    "name": "poly_port",
                    "lambda_U": 0.23,
                    "lambda_W": 0.31,
                    "lambda_R": 0.17,
                    "Omega_U": 2.7,
                    "Omega_W": 3.4,
                    "u_profile": u_profile,
                    "w_profile": w_profile,
                }
            ],
        )
        errors.append(max(abs(packet[key] - exact[key]) for key in ("B0", "Z0", "N0")))
        spacings.append(1.0 / (nodes - 1))

    orders = _observed_orders(errors, spacings)
    assert orders[-1] >= 1.95
    assert errors[-1] <= 2.0e-5
