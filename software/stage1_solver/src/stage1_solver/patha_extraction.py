"""Path-A field-to-D0 extraction machinery.

This module implements the post-freeze decision-11 section-5 calculator without
running the closed physical solve or any calibration root find.
"""

from __future__ import annotations

from dataclasses import dataclass
import copy
import hashlib
import json
import math
from typing import Any, Callable, Iterable, Mapping, Optional, Sequence

import numpy as np
from scipy.linalg import eigh
import sympy as sp
import torch

from .patha_static_balance import SSigmaProvider, SSigmaSpec, resolve_s_sigma


LANES = ("20", "21", "22")
AXISYM_LAMBDA = {"20": 1.0, "21": 0.5, "22": -1.0}


def stable_hash(obj: Any) -> str:
    """Return the oracle-compatible SHA256 hash of a JSON-serializable object."""

    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class WallEigenMatrices:
    """Generalized l=2 wall eigenproblem matrices with the mouth node removed."""

    A2: np.ndarray
    W: np.ndarray
    M_mu: np.ndarray
    w_nodes: np.ndarray
    R0_nodes: np.ndarray
    potential_midpoints: np.ndarray
    tension_midpoints: np.ndarray


@dataclass(frozen=True)
class WallEigenResult:
    """Lowest positive l=2 wall mode, normalized by chi.T W chi = 1."""

    K: float
    M: float
    eigenvalue: float
    chi: np.ndarray
    w_nodes: np.ndarray
    matrices: WallEigenMatrices
    normalization: float
    orientation_integral: float
    mode_index: int


def hooke_kappa_hat(*, a: float = 1.0, L: float = 37.0 / 20.0) -> float:
    """Analytic lowest l=2 Hooke eigenvalue per unit tau.

    The BVP is -eta'' + 7 eta/a^2 = kappa eta, eta(0)=0, eta'(L)=0.
    Hence eta ~ sin(k w), cos(k L)=0, and k=pi/(2L) for the lowest mode.
    """

    if a <= 0.0 or L <= 0.0:
        raise ValueError("a and L must be positive")
    return 7.0 / (a * a) + (math.pi / (2.0 * L)) ** 2


def analytic_hooke_l2_mode(w: np.ndarray, *, L: float = 37.0 / 20.0) -> np.ndarray:
    """Continuum-normalized analytic Hooke l=2 lowest mode on [0, L]."""

    w = np.asarray(w, dtype=np.float64)
    if L <= 0.0:
        raise ValueError("L must be positive")
    return math.sqrt(2.0 / L) * np.sin(math.pi * w / (2.0 * L))


def _as_provider(s_sigma: SSigmaSpec | SSigmaProvider | Mapping[str, Any]) -> SSigmaProvider:
    return s_sigma if isinstance(s_sigma, SSigmaProvider) else resolve_s_sigma(s_sigma)


def _provider_eval(
    fn: Callable[[torch.Tensor, torch.Tensor], torch.Tensor],
    R: np.ndarray | float,
    w: np.ndarray | float,
) -> np.ndarray:
    R_t = torch.as_tensor(R, dtype=torch.float64)
    w_t = torch.as_tensor(w, dtype=torch.float64)
    values = fn(R_t, w_t).detach().cpu().numpy()
    return np.asarray(values, dtype=np.float64)


def _validate_nodes(w_nodes: np.ndarray, R0_nodes: np.ndarray) -> None:
    if w_nodes.ndim != 1 or R0_nodes.ndim != 1:
        raise ValueError("w_nodes and R0_nodes must be one-dimensional")
    if w_nodes.shape != R0_nodes.shape:
        raise ValueError("R0_nodes must match w_nodes")
    if w_nodes.size < 4:
        raise ValueError("Need at least four wall nodes")
    if not np.all(np.isfinite(w_nodes)) or not np.all(np.isfinite(R0_nodes)):
        raise ValueError("w_nodes and R0_nodes must be finite")
    if np.any(np.diff(w_nodes) <= 0.0):
        raise ValueError("w_nodes must be strictly increasing")


def assemble_l2_wall_matrices(
    w_nodes: Sequence[float],
    R0_nodes: Sequence[float],
    s_sigma: SSigmaSpec | SSigmaProvider | Mapping[str, Any],
    *,
    r0_prime_nodes: Sequence[float] | None = None,
) -> WallEigenMatrices:
    """Assemble the symmetric l=2 stiffness, norm, and inertial mass matrices.

    The discretization is linear finite elements on wall nodes. The mouth node is
    eliminated for eta(0)=0; the exit condition T_w eta'(L)=0 is the natural FEM
    boundary condition.
    """

    w = np.asarray(w_nodes, dtype=np.float64)
    R0 = np.asarray(R0_nodes, dtype=np.float64)
    _validate_nodes(w, R0)
    provider = _as_provider(s_sigma)
    n_nodes = w.size
    n_unknown = n_nodes - 1
    A2 = np.zeros((n_unknown, n_unknown), dtype=np.float64)
    W = np.zeros_like(A2)
    M_mu = np.zeros_like(A2)

    if r0_prime_nodes is None:
        R0_prime = np.gradient(R0, w, edge_order=2)
    else:
        R0_prime = np.asarray(r0_prime_nodes, dtype=np.float64)
        if R0_prime.shape != R0.shape:
            raise ValueError("r0_prime_nodes must match w_nodes")

    q_nodes = _provider_eval(provider.T_w_R, R0, w) * R0_prime
    potential_midpoints = np.empty(n_nodes - 1, dtype=np.float64)
    tension_midpoints = np.empty(n_nodes - 1, dtype=np.float64)

    for element in range(n_nodes - 1):
        h = float(w[element + 1] - w[element])
        w_mid = 0.5 * (w[element] + w[element + 1])
        R_mid = 0.5 * (R0[element] + R0[element + 1])
        R_slope = (R0[element + 1] - R0[element]) / h
        dq_dw = (q_nodes[element + 1] - q_nodes[element]) / h

        T_w = float(_provider_eval(provider.T_w, R_mid, w_mid))
        T_omega = float(_provider_eval(provider.T_Omega, R_mid, w_mid))
        mu = float(_provider_eval(provider.mu, R_mid, w_mid))
        K_eta = (
            float(_provider_eval(provider.U_RR, R_mid, w_mid))
            - dq_dw
            + 0.5 * float(_provider_eval(provider.T_w_RR, R_mid, w_mid)) * R_slope**2
        )
        potential = K_eta + 6.0 * T_omega
        potential_midpoints[element] = potential
        tension_midpoints[element] = T_w

        local_stiffness = (T_w / h) * np.array([[1.0, -1.0], [-1.0, 1.0]])
        local_mass = (h / 6.0) * np.array([[2.0, 1.0], [1.0, 2.0]])
        local_A = local_stiffness + potential * local_mass
        local_M_mu = mu * local_mass
        for a_local, node_a in enumerate((element, element + 1)):
            if node_a == 0:
                continue
            ia = node_a - 1
            for b_local, node_b in enumerate((element, element + 1)):
                if node_b == 0:
                    continue
                ib = node_b - 1
                A2[ia, ib] += local_A[a_local, b_local]
                W[ia, ib] += local_mass[a_local, b_local]
                M_mu[ia, ib] += local_M_mu[a_local, b_local]

    return WallEigenMatrices(
        A2=A2,
        W=W,
        M_mu=M_mu,
        w_nodes=w,
        R0_nodes=R0,
        potential_midpoints=potential_midpoints,
        tension_midpoints=tension_midpoints,
    )


def solve_l2_wall_eigenproblem(
    w_nodes: Sequence[float],
    R0_nodes: Sequence[float],
    s_sigma: SSigmaSpec | SSigmaProvider | Mapping[str, Any],
    *,
    positive_tol: float = 1.0e-12,
) -> WallEigenResult:
    """Solve A2 chi = K W chi for the lowest positive l=2 wall mode."""

    matrices = assemble_l2_wall_matrices(w_nodes, R0_nodes, s_sigma)
    evals, evecs = eigh(matrices.A2, matrices.W, check_finite=True)
    positive = np.flatnonzero(evals > positive_tol)
    if positive.size == 0:
        raise ValueError("No positive l=2 wall eigenvalue found")
    mode_index = int(positive[0])
    chi_unknown = np.asarray(evecs[:, mode_index], dtype=np.float64)
    normalization = math.sqrt(float(chi_unknown @ matrices.W @ chi_unknown))
    if normalization <= 0.0 or not math.isfinite(normalization):
        raise ValueError("Degenerate eigenvector normalization")
    chi_unknown = chi_unknown / normalization
    chi = np.zeros(matrices.w_nodes.size, dtype=np.float64)
    chi[1:] = chi_unknown
    orientation_integral = float(np.trapz(chi, matrices.w_nodes))
    if orientation_integral < 0.0:
        chi = -chi
        chi_unknown = -chi_unknown
        orientation_integral = -orientation_integral
    K = float(chi_unknown @ matrices.A2 @ chi_unknown)
    M = float(chi_unknown @ matrices.M_mu @ chi_unknown)
    return WallEigenResult(
        K=K,
        M=M,
        eigenvalue=float(evals[mode_index]),
        chi=chi,
        w_nodes=matrices.w_nodes,
        matrices=matrices,
        normalization=float(chi_unknown @ matrices.W @ chi_unknown),
        orientation_integral=orientation_integral,
        mode_index=mode_index,
    )


def _as_array(values: Sequence[float], shape: tuple[int, ...] | None = None) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if shape is not None and arr.shape != shape:
        raise ValueError(f"array shape {arr.shape} does not match expected {shape}")
    if not np.all(np.isfinite(arr)):
        raise ValueError("profile arrays must be finite")
    return arr


def trapezoid_product_integral(
    w: Sequence[float],
    profiles: Sequence[Sequence[float]],
    *,
    weight: Sequence[float] | None = None,
) -> float:
    """Composite trapezoid integral of weight times a product of sampled fields."""

    w_arr = _as_array(w)
    if w_arr.ndim != 1 or w_arr.size < 2:
        raise ValueError("w must be a one-dimensional grid with at least two nodes")
    product = np.ones_like(w_arr)
    if weight is not None:
        product *= _as_array(weight, w_arr.shape)
    for profile in profiles:
        product *= _as_array(profile, w_arr.shape)
    return float(np.trapz(product, w_arr))


def field_to_reduced_couplings(
    w: Sequence[float],
    chi: Sequence[float],
    *,
    bdg_modes: Sequence[Mapping[str, Any]] = (),
    mixed_ports: Sequence[Mapping[str, Any]] = (),
    weight: Sequence[float] | None = None,
) -> dict[str, Any]:
    """Compute v2_22a-style overlaps and reduced couplings from sampled fields."""

    w_arr = _as_array(w)
    chi_arr = _as_array(chi, w_arr.shape)
    reduced_bdg: list[dict[str, Any]] = []
    reduced_ports: list[dict[str, Any]] = []
    overlap_summary: dict[str, float] = {}

    for mode in bdg_modes:
        name = str(mode.get("name", f"bdg_{len(reduced_bdg)}"))
        phi = _as_array(mode["profile"], w_arr.shape)
        I_eta_phi = trapezoid_product_integral(w_arr, [chi_arr, phi], weight=weight)
        coupling = float(mode["lambda_B"]) * I_eta_phi
        overlap_summary[f"I_eta_phi::{name}"] = I_eta_phi
        reduced_bdg.append(
            {
                "name": name,
                "coupling": coupling,
                "varpi": float(mode["varpi"]),
                "overlap_I_eta_phi": I_eta_phi,
                "lambda_B_eff": float(mode["lambda_B"]),
            }
        )

    for port in mixed_ports:
        name = str(port.get("name", f"port_{len(reduced_ports)}"))
        u = _as_array(port["u_profile"], w_arr.shape)
        w_profile = _as_array(port["w_profile"], w_arr.shape)
        I_eta_u = trapezoid_product_integral(w_arr, [chi_arr, u], weight=weight)
        I_eta_w = trapezoid_product_integral(w_arr, [chi_arr, w_profile], weight=weight)
        I_u_w = trapezoid_product_integral(w_arr, [u, w_profile], weight=weight)
        overlap_summary[f"I_eta_u::{name}"] = I_eta_u
        overlap_summary[f"I_eta_w::{name}"] = I_eta_w
        overlap_summary[f"I_u_w::{name}"] = I_u_w
        reduced_ports.append(
            {
                "name": name,
                "Omega_U": float(port["Omega_U"]),
                "Omega_W": float(port["Omega_W"]),
                "R": float(port["lambda_R"]) * I_u_w,
                "g_U": float(port["lambda_U"]) * I_eta_u,
                "g_W": float(port["lambda_W"]) * I_eta_w,
                "overlap_I_eta_u": I_eta_u,
                "overlap_I_eta_w": I_eta_w,
                "overlap_I_u_w": I_u_w,
                "lambda_U_eff": float(port["lambda_U"]),
                "lambda_W_eff": float(port["lambda_W"]),
                "lambda_R_eff": float(port["lambda_R"]),
            }
        )

    return {
        "bdg_modes": reduced_bdg,
        "mixed_ports": reduced_ports,
        "overlaps": overlap_summary,
    }


def coefficients_from_sampled_fields(
    w: Sequence[float],
    chi: Sequence[float],
    *,
    K: float,
    M: float,
    bdg_modes: Sequence[Mapping[str, Any]] = (),
    mixed_ports: Sequence[Mapping[str, Any]] = (),
    weight: Sequence[float] | None = None,
) -> dict[str, Any]:
    """Compute the lane-level coefficient packet from sampled profiles."""

    reduced = field_to_reduced_couplings(
        w,
        chi,
        bdg_modes=bdg_modes,
        mixed_ports=mixed_ports,
        weight=weight,
    )
    lane = {
        "K": float(K),
        "M": float(M),
        "bdg_modes": reduced["bdg_modes"],
        "mixed_ports": reduced["mixed_ports"],
    }
    packet = lane_extract(lane)
    packet["overlaps"] = reduced["overlaps"]
    return packet


def grouped_decomposition(values: Mapping[str, float]) -> dict[str, float]:
    """Weighted grouped-P2 decomposition for values keyed by 20/21/22."""

    x20 = float(values["20"])
    x21 = float(values["21"])
    x22 = float(values["22"])
    xbar = (x20 + 2.0 * x21 + 2.0 * x22) / 5.0
    a_x = (2.0 * x20 - x21 - x22) / 10.0
    b_x = (x21 - x22) / 2.0
    return {
        "bar": xbar,
        "a": a_x,
        "b": b_x,
        "anisotropy_norm_sq": 4.0 * a_x * a_x + (4.0 / 5.0) * b_x * b_x,
        "axisymmetric_b_minus_3a": b_x - 3.0 * a_x,
    }


def lane_values(packet_by_lane: Mapping[str, Mapping[str, float]], key: str) -> dict[str, float]:
    return {lane: float(packet_by_lane[lane][key]) for lane in LANES}


def bdg_moments(modes: Iterable[Mapping[str, float]]) -> dict[str, float]:
    B0 = B2 = B4 = 0.0
    for mode in modes:
        coupling = float(mode["coupling"])
        varpi = float(mode["varpi"])
        if varpi <= 0.0:
            raise ValueError(f"BdG mode has nonpositive varpi={varpi}")
        B0 += coupling**2 / varpi**2
        B2 += coupling**2 / varpi**4
        B4 += coupling**2 / varpi**6
    return {"B0": B0, "B2": B2, "B4": B4}


def mixed_port_moments(
    ports: Iterable[Mapping[str, float]],
) -> tuple[dict[str, float], list[dict[str, float]]]:
    Z0 = Z2 = Z4 = N0 = N2 = N4 = 0.0
    diagnostics: list[dict[str, float]] = []
    for idx, port in enumerate(ports):
        OU = float(port["Omega_U"])
        OW = float(port["Omega_W"])
        R = float(port.get("R", 0.0))
        gU = float(port.get("g_U", 0.0))
        gW = float(port.get("g_W", 0.0))
        Delta = OU**2 * OW**2 - R**2
        if Delta == 0.0:
            raise ValueError(f"mixed port {idx} has Delta=0")
        S = OU**2 + OW**2
        Q = gU**2 * OW**2 + 2.0 * gU * gW * R + gW**2 * OU**2
        H = gU**2 + gW**2
        P = OU**2 * gW + R * gU
        Z0_r = Q / Delta
        Z2_r = (Q * S - H * Delta) / Delta**2
        Z4_r = (Q * (S**2 - Delta) - S * H * Delta) / Delta**3
        N0_r = P**2 / Delta**2
        N2_r = 2.0 * P * (P * S - Delta * gW) / Delta**3
        N4_r = (
            Delta**2 * gW**2
            - 2.0 * Delta * P**2
            - 4.0 * Delta * P * S * gW
            + 3.0 * P**2 * S**2
        ) / Delta**4
        Z0 += Z0_r
        Z2 += Z2_r
        Z4 += Z4_r
        N0 += N0_r
        N2 += N2_r
        N4 += N4_r
        diagnostics.append(
            {
                "port_index": float(idx),
                "Delta": Delta,
                "S": S,
                "Q": Q,
                "H": H,
                "P": P,
                "Z0": Z0_r,
                "Z2": Z2_r,
                "Z4": Z4_r,
                "N0": N0_r,
                "N2": N2_r,
                "N4": N4_r,
            }
        )
    return {"Z0": Z0, "Z2": Z2, "Z4": Z4, "N0": N0, "N2": N2, "N4": N4}, diagnostics


def lane_extract(lane: Mapping[str, Any]) -> dict[str, Any]:
    """Extract the v2_21 lane-level coefficient and observable algebra."""

    if "direct_coefficients" in lane:
        coeff = {k: float(v) for k, v in lane["direct_coefficients"].items()}
        required = {"K", "M", "B0", "B2", "B4", "Z0", "Z2", "Z4", "N0", "N2", "N4"}
        missing = required - set(coeff)
        if missing:
            raise ValueError(f"direct_coefficients missing {sorted(missing)}")
        port_diagnostics: list[dict[str, float]] = []
    else:
        K = float(lane["K"])
        M = float(lane["M"])
        b = bdg_moments(lane.get("bdg_modes", []))
        z_n, port_diagnostics = mixed_port_moments(lane.get("mixed_ports", []))
        coeff = {"K": K, "M": M, **b, **z_n}

    K = coeff["K"]
    M = coeff["M"]
    B0 = coeff["B0"]
    B2 = coeff["B2"]
    B4 = coeff["B4"]
    Z0 = coeff["Z0"]
    Z2 = coeff["Z2"]
    Z4 = coeff["Z4"]
    N0 = coeff["N0"]
    N2 = coeff["N2"]
    N4 = coeff["N4"]

    D0 = K - B0 - Z0
    D2 = -(M + B2 + Z2)
    D4 = -(B4 + Z4)
    u2 = -D2 / D0 if D0 != 0.0 else math.nan
    u4 = (D2**2 - D0 * D4) / D0**2 if D0 != 0.0 else math.nan
    P0 = N0 / D0 if D0 != 0.0 else math.nan
    P2 = (D0 * N2 - 2.0 * D2 * N0) / D0**2 if D0 != 0.0 else math.nan
    P4 = (
        (D0**2 * N4 - 2.0 * D0 * (D2 * N2 + D4 * N0) + 3.0 * D2**2 * N0)
        / D0**3
        if D0 != 0.0
        else math.nan
    )

    A = M + B2 + Z2
    C = B4 + Z4
    R_pole = D0 * C - 3.0 * A**2
    constant_prefactor_N2_residual = N2 - 2.0 * D2 * N0 / D0 if D0 != 0.0 else math.nan
    constant_prefactor_N4_residual = (
        N4 - N0 * (D2**2 + 2.0 * D0 * D4) / D0**2 if D0 != 0.0 else math.nan
    )

    stability = {
        "D0_positive": D0 > 0.0,
        "C_positive": C > 0.0,
        "N0_nonzero": N0 != 0.0,
        "wall_M_positive": M > 0.0,
        "all_Delta_positive": all(p["Delta"] > 0.0 for p in port_diagnostics),
    }
    return {
        **coeff,
        "D0": D0,
        "D2": D2,
        "D4": D4,
        "u2": u2,
        "u4": u4,
        "P0": P0,
        "P2": P2,
        "P4": P4,
        "A_MplusB2plusZ2": A,
        "C_B4plusZ4": C,
        "R_pole": R_pole,
        "constant_prefactor_N2_residual": constant_prefactor_N2_residual,
        "constant_prefactor_N4_residual": constant_prefactor_N4_residual,
        "port_diagnostics": port_diagnostics,
        "stability": stability,
    }


def observable_residuals(
    *,
    P0: float,
    D0: float,
    B4: float,
    Z4: float,
    M: float,
    B2: float,
    Z2: float,
    G: float = 1.0,
    c_s: float = 1.0,
    c: float = 1.0,
    a: float = 1.0,
    mhat0: float = 1.0,
    S_port: float = 1.0,
) -> dict[str, float]:
    """Decision-11 residual formulas from an already extracted lane/group."""

    R_norm = mhat0**2 * S_port * P0 - 54.0 * G * c_s**5 / (5.0 * a**5 * c**5)
    R_pole = D0 * (B4 + Z4) - 3.0 * (M + B2 + Z2) ** 2
    return {"R_norm": R_norm, "R_pole": R_pole}


def _input_hash_for_branch(branch: Mapping[str, Any], constants: Mapping[str, Any]) -> str:
    """Match the historical oracle input-hash object for each fixture schema."""

    if "profile_adapter" in branch:
        return stable_hash(
            {
                "name": branch["name"],
                "geometry": branch.get("geometry"),
                "target": branch.get("target"),
                "lanes": branch.get("lanes"),
            }
        )
    return stable_hash(
        {
            "branch_name": branch.get("name", "unnamed"),
            "geometry": branch.get("geometry", {}),
            "target_constants": constants,
            "lanes_input": branch["lanes"],
        }
    )


def extract_branch(branch: Mapping[str, Any], tol: float = 1.0e-9) -> dict[str, Any]:
    """Extract the full grouped v2_21 observable packet for one branch."""

    lanes_in = branch["lanes"]
    if not all(lane in lanes_in for lane in LANES):
        raise ValueError(f"branch lanes must include {LANES}")
    lane_packets = {lane: lane_extract(lanes_in[lane]) for lane in LANES}
    grouped: dict[str, dict[str, float]] = {}
    for key in [
        "K",
        "M",
        "B0",
        "B2",
        "B4",
        "Z0",
        "Z2",
        "Z4",
        "N0",
        "N2",
        "N4",
        "D0",
        "D2",
        "D4",
        "u2",
        "u4",
        "P0",
        "P2",
        "P4",
        "R_pole",
    ]:
        grouped[key] = grouped_decomposition(lane_values(lane_packets, key))

    constants = branch.get("target", {}).get("constants", {})
    G = float(constants.get("G", 1.0))
    c_s = float(constants.get("c_s", 1.0))
    c = float(constants.get("c", 1.0))
    a = float(constants.get("a", 1.0))
    mhat0 = float(constants.get("mhat0", 1.0))
    S_port = float(constants.get("S_port", 1.0))
    theta_tail = float(constants.get("theta_tail", 1.0))

    D0_bar = grouped["D0"]["bar"]
    B4_bar = grouped["B4"]["bar"]
    Z4_bar = grouped["Z4"]["bar"]
    A_bar = grouped["M"]["bar"] + grouped["B2"]["bar"] + grouped["Z2"]["bar"]
    P0_bar = grouped["P0"]["bar"]
    P2_bar = grouped["P2"]["bar"]
    P4_bar = grouped["P4"]["bar"]

    P0_target = 54.0 * G * c_s**5 / (5.0 * a**5 * c**5)
    gamma_eff = mhat0**2 * S_port * P0_bar * a**5 / (27.0 * c_s**5)
    gamma_GR = 2.0 * G / (5.0 * c**5)
    R_pole_iso = D0_bar * (B4_bar + Z4_bar) - 3.0 * A_bar**2
    R_norm = mhat0**2 * S_port * P0_bar - P0_target
    R_tail = theta_tail * (c / c_s) ** 3 - 1.0

    geometry = branch.get("geometry", {})
    open_gate = {
        "R_exit_positive": float(geometry.get("R_exit", 0.0)) > 0.0,
        "boundary_class_open_impedance": geometry.get("boundary_class") == "open_impedance",
        "hard_cap_forbidden": geometry.get("boundary_class") != "hard_cap"
        and float(geometry.get("R_exit", 0.0)) > 0.0,
    }
    lane_stability_all = {
        name: all(packet["stability"].get(name, False) for packet in lane_packets.values())
        for name in sorted({k for packet in lane_packets.values() for k in packet["stability"]})
    }
    residuals = {
        "R_pole": R_pole_iso,
        "R_norm": R_norm,
        "R_P2": P2_bar,
        "R_P4": P4_bar,
        "R_tail": R_tail,
        "gamma_eff_minus_gamma_GR": gamma_eff - gamma_GR,
    }
    pass_flags = {
        "open_gate_pass": all(open_gate.values()),
        "stability_gate_pass": all(lane_stability_all.values()),
        "isotropic_D0_pass": abs(grouped["D0"]["anisotropy_norm_sq"]) <= tol,
        "isotropic_N0_pass": abs(grouped["N0"]["anisotropy_norm_sq"]) <= tol,
        "axisymmetric_P0_pattern_pass": abs(grouped["P0"]["axisymmetric_b_minus_3a"]) <= tol,
        "one_pole_pass": abs(R_pole_iso) <= tol,
        "normalization_pass": abs(R_norm) <= tol,
        "constant_prefactor_P2_pass": abs(P2_bar) <= tol,
        "constant_prefactor_P4_pass": abs(P4_bar) <= tol,
        "tail_transport_pass": abs(R_tail) <= tol,
    }
    pass_flags["target_packet_pass"] = all(
        pass_flags[k]
        for k in [
            "open_gate_pass",
            "stability_gate_pass",
            "isotropic_D0_pass",
            "isotropic_N0_pass",
            "one_pole_pass",
            "normalization_pass",
            "constant_prefactor_P2_pass",
            "constant_prefactor_P4_pass",
            "tail_transport_pass",
        ]
    )
    return {
        "name": branch.get("name", "unnamed"),
        "input_hash": _input_hash_for_branch(branch, constants),
        "open_gate": open_gate,
        "lane_packets": lane_packets,
        "grouped": grouped,
        "target_values": {
            "P0_target": P0_target,
            "gamma_eff": gamma_eff,
            "gamma_GR": gamma_GR,
        },
        "residuals": residuals,
        "stability": lane_stability_all,
        "pass_flags": pass_flags,
    }


def extract_manifest(manifest: Mapping[str, Any], tol: float = 1.0e-9) -> list[dict[str, Any]]:
    return [extract_branch(branch, tol=tol) for branch in manifest.get("branches", [])]


def builtin_expr(name: str, s: sp.Symbol, L: sp.Symbol) -> sp.Expr:
    """Return a SymPy expression for a named v2_22a finite-throat profile."""

    if name in {"chi_eta_sin", "wall_sin", "u_wall"}:
        return sp.sqrt(2 / L) * sp.sin(sp.pi * s / L)
    if name in {"phi_DN", "w_DN", "half_wave_DN"}:
        return sp.sqrt(2 / L) * sp.sin(sp.pi * s / (2 * L))
    if name in {"phi_ND", "cos_half"}:
        return sp.sqrt(2 / L) * sp.cos(sp.pi * s / (2 * L))
    if name in {"constant_unit", "one"}:
        return sp.Integer(1)
    raise ValueError(f"unknown built-in profile name: {name}")


def analytic_expr(profile: Mapping[str, Any], s: sp.Symbol, L: sp.Symbol) -> Optional[sp.Expr]:
    kind = profile.get("kind", "builtin")
    if kind == "builtin":
        return builtin_expr(str(profile["name"]), s, L)
    if kind == "expr":
        namespace = {
            "s": s,
            "L": L,
            "pi": sp.pi,
            "sqrt": sp.sqrt,
            "sin": sp.sin,
            "cos": sp.cos,
            "exp": sp.exp,
        }
        return sp.sympify(profile["expr"], locals=namespace)
    if kind == "sampled":
        return None
    raise ValueError(f"unknown profile kind: {kind}")


def make_numeric_evaluator(profile: Mapping[str, Any], L_value: float) -> Callable[[float], float]:
    kind = profile.get("kind", "builtin")
    if kind in {"builtin", "expr"}:
        s_sym, L_sym = sp.symbols("s L", positive=True)
        expr = analytic_expr(profile, s_sym, L_sym)
        assert expr is not None
        f = sp.lambdify(s_sym, expr.subs(L_sym, L_value), "math")
        return lambda x: float(f(float(x)))
    if kind == "sampled":
        samples = sorted((float(p[0]), float(p[1])) for p in profile["samples"])
        if len(samples) < 2:
            raise ValueError("sampled profile requires at least two points")

        def interp(x: float) -> float:
            x = float(x)
            if x <= samples[0][0]:
                return samples[0][1]
            if x >= samples[-1][0]:
                return samples[-1][1]
            lo, hi = 0, len(samples) - 1
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if samples[mid][0] <= x:
                    lo = mid
                else:
                    hi = mid
            x0, y0 = samples[lo]
            x1, y1 = samples[hi]
            t = (x - x0) / (x1 - x0)
            return y0 * (1.0 - t) + y1 * t

        return interp
    raise ValueError(f"unknown profile kind: {kind}")


def integrate_product(
    profiles: Sequence[Mapping[str, Any]],
    weight: Mapping[str, Any],
    L_value: float,
    grid_points: int = 4001,
) -> dict[str, Any]:
    """v2_22a profile integral: exact SymPy when analytic, trapezoid otherwise."""

    s_sym, L_sym = sp.symbols("s L", positive=True)
    exprs: list[sp.Expr] = []
    analytic = True
    for spec in [weight, *profiles]:
        expr = analytic_expr(spec, s_sym, L_sym)
        if expr is None:
            analytic = False
            break
        exprs.append(expr)

    if analytic:
        integrand = sp.prod(exprs)
        exact = sp.integrate(integrand, (s_sym, 0, L_sym))
        exact_at_L = sp.simplify(exact.subs(L_sym, sp.nsimplify(L_value)))
        return {
            "method": "analytic_sympy",
            "exact": str(sp.simplify(exact)),
            "value": float(sp.N(exact_at_L, 30)),
        }

    evaluators = [make_numeric_evaluator(spec, L_value) for spec in [weight, *profiles]]
    n = int(grid_points)
    if n < 2:
        raise ValueError("grid_points must be >= 2")
    h = L_value / (n - 1)
    total = 0.0
    for i in range(n):
        x = i * h
        y = 1.0
        for ev in evaluators:
            y *= ev(x)
        total += (0.5 if i == 0 or i == n - 1 else 1.0) * y
    return {"method": "trapezoid", "grid_points": n, "value": total * h}


def profile_lookup(branch: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    profiles = branch["profiles"]
    if name not in profiles:
        raise KeyError(f"profile {name!r} is not defined in branch {branch.get('name', 'unnamed')}")
    return profiles[name]


def lane_scaled(base: float, branch: Mapping[str, Any], lane: str, key: str) -> float:
    weak = branch.get("weak_axisymmetric")
    if not weak:
        return float(base)
    eps = float(weak.get("epsilon", 0.0))
    slope = float(weak.get("slopes", {}).get(key, 0.0))
    return float(base) * (1.0 + eps * AXISYM_LAMBDA[lane] * slope)


def compute_branch_overlaps(branch: Mapping[str, Any]) -> dict[str, Any]:
    L_value = float(branch.get("geometry", {}).get("L", 1.0))
    grid_points = int(branch.get("integration", {}).get("grid_points", 4001))
    weight = branch["profiles"].get("weight", {"kind": "builtin", "name": "one"})
    wall_profile = profile_lookup(branch, "wall")

    overlaps: dict[str, Any] = {}
    for mode in branch["reduction"].get("bdg_modes", []):
        phi = profile_lookup(branch, mode["profile"])
        overlaps[f"I_eta_phi::{mode['name']}"] = integrate_product(
            [wall_profile, phi], weight, L_value, grid_points
        )
    for port in branch["reduction"].get("mixed_ports", []):
        u = profile_lookup(branch, port["u_profile"])
        w_profile = profile_lookup(branch, port["w_profile"])
        overlaps[f"I_eta_u::{port['name']}"] = integrate_product(
            [wall_profile, u], weight, L_value, grid_points
        )
        overlaps[f"I_eta_w::{port['name']}"] = integrate_product(
            [wall_profile, w_profile], weight, L_value, grid_points
        )
        uw_names = port.get("u_w_profiles", [port["u_profile"], port["w_profile"]])
        uw_0 = profile_lookup(branch, uw_names[0])
        uw_1 = profile_lookup(branch, uw_names[1])
        overlaps[f"I_u_w::{port['name']}"] = integrate_product(
            [uw_0, uw_1], weight, L_value, grid_points
        )
    return overlaps


def adapt_profile_branch_to_v21(branch: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Convert one v2_22a profile branch into a v2_21 branch manifest."""

    reduction = branch["reduction"]
    source_branch_hash = stable_hash(branch)
    if "direct_coefficients" in reduction:
        direct_coefficients = {k: float(v) for k, v in reduction["direct_coefficients"].items()}
        v21_branch = {
            "name": branch["name"],
            "geometry": {k: v for k, v in branch.get("geometry", {}).items() if k != "L"},
            "target": branch.get("target", {}),
            "lanes": {lane: {"direct_coefficients": copy.deepcopy(direct_coefficients)} for lane in LANES},
            "profile_adapter": {
                "source_branch_hash": source_branch_hash,
                "coefficient_path": "direct_derived_coefficients",
                "derived_coefficients_provenance": reduction.get(
                    "derived_coefficients_provenance", {}
                ),
                "weak_axisymmetric": branch.get("weak_axisymmetric", None),
            },
        }
        diagnostics = {
            "source_branch": branch["name"],
            "source_branch_hash": source_branch_hash,
            "coefficient_path": "direct_derived_coefficients",
            "direct_coefficients": direct_coefficients,
            "v21_branch_hash": stable_hash(v21_branch),
        }
        return v21_branch, diagnostics

    overlaps = compute_branch_overlaps(branch)
    lanes: dict[str, Any] = {}
    for lane in LANES:
        wall = reduction["wall"]
        lane_packet = {
            "K": lane_scaled(float(wall["K"]), branch, lane, "K"),
            "M": lane_scaled(float(wall["M"]), branch, lane, "M"),
            "bdg_modes": [],
            "mixed_ports": [],
        }
        for mode in reduction.get("bdg_modes", []):
            I = float(overlaps[f"I_eta_phi::{mode['name']}"]["value"])
            lam_B = lane_scaled(float(mode["lambda_B"]), branch, lane, "lambda_B")
            varpi = lane_scaled(float(mode["varpi"]), branch, lane, "varpi")
            lane_packet["bdg_modes"].append(
                {
                    "name": mode["name"],
                    "coupling": lam_B * I,
                    "varpi": varpi,
                    "overlap_I_eta_phi": I,
                    "lambda_B_eff": lam_B,
                }
            )
        for port in reduction.get("mixed_ports", []):
            I_eta_u = float(overlaps[f"I_eta_u::{port['name']}"]["value"])
            I_eta_w = float(overlaps[f"I_eta_w::{port['name']}"]["value"])
            I_u_w = float(overlaps[f"I_u_w::{port['name']}"]["value"])
            lam_U = lane_scaled(float(port["lambda_U"]), branch, lane, "lambda_U")
            lam_W = lane_scaled(float(port["lambda_W"]), branch, lane, "lambda_W")
            lam_R = lane_scaled(float(port["lambda_R"]), branch, lane, "lambda_R")
            OU = lane_scaled(float(port["Omega_U"]), branch, lane, "Omega_U")
            OW = lane_scaled(float(port["Omega_W"]), branch, lane, "Omega_W")
            lane_packet["mixed_ports"].append(
                {
                    "name": port["name"],
                    "Omega_U": OU,
                    "Omega_W": OW,
                    "R": lam_R * I_u_w,
                    "g_U": lam_U * I_eta_u,
                    "g_W": lam_W * I_eta_w,
                    "overlap_I_eta_u": I_eta_u,
                    "overlap_I_eta_w": I_eta_w,
                    "overlap_I_u_w": I_u_w,
                    "lambda_U_eff": lam_U,
                    "lambda_W_eff": lam_W,
                    "lambda_R_eff": lam_R,
                }
            )
        lanes[lane] = lane_packet

    v21_branch = {
        "name": branch["name"],
        "geometry": {k: v for k, v in branch.get("geometry", {}).items() if k != "L"},
        "target": branch.get("target", {}),
        "lanes": lanes,
        "profile_adapter": {
            "source_branch_hash": source_branch_hash,
            "overlap_summary": overlaps,
            "weak_axisymmetric": branch.get("weak_axisymmetric", None),
        },
    }
    diagnostics = {
        "source_branch": branch["name"],
        "source_branch_hash": source_branch_hash,
        "overlaps": overlaps,
        "v21_branch_hash": stable_hash(v21_branch),
    }
    return v21_branch, diagnostics


def adapt_profile_manifest(profile_manifest: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    v21_branches: list[dict[str, Any]] = []
    diagnostics: list[dict[str, Any]] = []
    for branch in profile_manifest.get("branches", []):
        v21_branch, diag = adapt_profile_branch_to_v21(branch)
        v21_branches.append(v21_branch)
        diagnostics.append(diag)
    v21_manifest = {
        "schema": "stage_v2_21_branch_extraction_fixture/v1",
        "description": "Generated by stage_v2_22a_profile_to_coefficient_adapter.py from profile overlap data.",
        "source_profile_manifest_hash": stable_hash(profile_manifest),
        "branches": v21_branches,
    }
    adapter_diag = {
        "profile_manifest_hash": stable_hash(profile_manifest),
        "v21_manifest_hash": stable_hash(v21_manifest),
        "branch_diagnostics": diagnostics,
    }
    return v21_manifest, adapter_diag
