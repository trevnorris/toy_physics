"""Path-A B2b Maxwell transfer derivation and validation.

This module derives the low-frequency Maxwell transfer moments
``{Z0,Z2,Z4,N0,N2,N4}`` on the self-consistent Path-A closed background.
It deliberately stops at the B1 ``direct_coefficients`` shape; downstream
``D0``/``R_norm`` assembly is B2c scope.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.linalg import solve as dense_solve

from . import patha_b2a_bdg as b2a


DEFAULT_RUN_ROOT = Path("software/stage1_solver/runs/patha_b2b_maxwell_transfer")
DEFAULT_REPORT_PATH = Path("software/stage1_solver/reports/patha_b2b_maxwell_transfer_derivation_report.md")
DEFAULT_B2A_RUN_ROOT = b2a.DEFAULT_RUN_ROOT
FROZEN_A = b2a.FROZEN_A
FROZEN_L = b2a.FROZEN_L
SMOKE_RESIDUAL_REFERENCE = b2a.SMOKE_RESIDUAL_REFERENCE
L_HARM = 6.0
XI_GAUGE = 1.0
CS_CONST = 1.0
COEFF_KEYS = ("Z0", "Z2", "Z4", "N0", "N2", "N4")
Z_KEYS = ("Z0", "Z2", "Z4")
N_KEYS = ("N0", "N2", "N4")
DIRECT_KEYS = ("K", "M", "B0", "B2", "B4", *COEFF_KEYS)
DEFAULT_TRANSFER_GRIDS = ((31, 13), (43, 17), (47, 17))
DEFAULT_WINDOWS = (0.058, 0.046, 0.036, 0.028)
DEFAULT_TRUNCATION_SCALES = (3.00, 4.00, 5.00)
DEFAULT_FINAL_GRID = (47, 17)
DEFAULT_FINAL_WINDOW = 0.028
DEFAULT_FINAL_TRUNCATION = 5.00
DEFAULT_MODE_COUNT = 30
DUAL_ENGINE_TOLERANCES = {"abs": 2.0e-9, "rel": 5.0e-2}
CONSUMER_TOLERANCES = {"abs": 2.0e-9, "rel": 2.0e-7}
CONVERGENCE_REL_TOL = 5.0e-2
CONVERGENCE_SHRINK_FACTOR = 0.98
EXTERIOR_SPONGE_WIDTH_FRACTION = 0.40
HIGH_ORDER_BOUNDARY_CLOSURE_MARGIN = 8


def _json_default(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, complex):
        return {"re": float(value.real), "im": float(value.imag)}
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _full_stable_hash(obj: Any) -> str:
    text = json.dumps(obj, sort_keys=True, separators=(",", ":"), default=_json_default)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )


def _format_tau(tau: float) -> str:
    return b2a._format_tau(tau)


def _parse_grid(text: str) -> tuple[int, int]:
    return b2a._parse_grid(text)


def _parse_grid_list(text: str) -> list[tuple[int, int]]:
    return [_parse_grid(piece) for piece in text.split(";") if piece]


def _parse_float_list(text: str) -> list[float]:
    return [float(piece) for piece in text.split(";") if piece]


def _to_float(value: Any) -> float:
    if isinstance(value, Mapping):
        return float(value.get("re", value.get("real", 0.0)))
    return float(value)


def _complex_parts(value: Any) -> tuple[float, float]:
    if isinstance(value, Mapping):
        return (
            float(value.get("re", value.get("real", 0.0))),
            float(value.get("im", value.get("imag", 0.0))),
        )
    return float(value), 0.0


def _interp1_clamped(x_old: np.ndarray, y_old: np.ndarray, x_new: np.ndarray | float) -> np.ndarray:
    return np.interp(x_new, x_old, y_old, left=y_old[0], right=y_old[-1])


def _interp2_clamped(
    r_old: np.ndarray,
    w_old: np.ndarray,
    values: np.ndarray,
    r_new: np.ndarray,
    w_new: np.ndarray,
) -> np.ndarray:
    rows = np.empty((r_old.size, w_new.size), dtype=np.result_type(values, np.float64))
    for i in range(r_old.size):
        rows[i, :] = np.interp(w_new, w_old, values[i, :], left=values[i, 0], right=values[i, -1])
    out = np.empty((r_new.size, w_new.size), dtype=np.result_type(values, np.float64))
    for j in range(w_new.size):
        out[:, j] = np.interp(r_new, r_old, rows[:, j], left=rows[0, j], right=rows[-1, j])
    return out


def _radial_exterior_taper(r_old: np.ndarray, r_new: np.ndarray) -> np.ndarray:
    """Smoothly continue compact background/source fields beyond the B2a export."""
    r_edge = float(np.max(r_old))
    span = max(float(np.max(r_old) - np.min(r_old)), 1.0e-12)
    width = EXTERIOR_SPONGE_WIDTH_FRACTION * span
    distance = np.maximum(np.asarray(r_new, dtype=np.float64) - r_edge, 0.0)
    return np.exp(-((distance / width) ** 2))


def _flatten_w_major(values: np.ndarray) -> np.ndarray:
    return np.asarray(values).T.reshape(-1)


def _cell_index(i: int, j: int, nr: int) -> int:
    return j * nr + i


def _diag(values: Sequence[complex] | np.ndarray) -> np.ndarray:
    return np.diag(np.asarray(values, dtype=np.complex128))


def _first_derivative_matrix(n: int, h: float) -> np.ndarray:
    mat = np.zeros((n, n), dtype=np.complex128)
    for i in range(n):
        if i < n - 1:
            mat[i, i + 1] += 1.0 / (2.0 * h)
        if i > 0:
            mat[i, i - 1] -= 1.0 / (2.0 * h)
    return mat


def _second_derivative_matrix(n: int, h: float) -> np.ndarray:
    mat = np.zeros((n, n), dtype=np.complex128)
    for i in range(n):
        mat[i, i] -= 2.0 / h**2
        if i < n - 1:
            mat[i, i + 1] += 1.0 / h**2
        if i > 0:
            mat[i, i - 1] += 1.0 / h**2
    return mat


def _finite_difference_weights(x0: float, xs: np.ndarray, order: int) -> np.ndarray:
    """Return finite-difference weights for derivative ``order`` at ``x0``."""
    m = int(order)
    n = int(xs.size)
    vand = np.vstack([(np.asarray(xs, dtype=np.float64) - float(x0)) ** k for k in range(n)])
    rhs = np.zeros(n, dtype=np.float64)
    rhs[m] = math.factorial(m)
    return np.linalg.solve(vand, rhs)


def _high_order_derivative_matrix(n: int, h: float, order: int) -> np.ndarray:
    """Fourth-order interior derivative matrix with stable second-order boundary closure."""
    if n < 5:
        raise ValueError("high-order derivative matrix requires at least five points")
    margin = min(HIGH_ORDER_BOUNDARY_CLOSURE_MARGIN, max(2, (n - 1) // 2))
    if order == 1:
        mat = _first_derivative_matrix(n, h)
        weights = np.asarray([1.0, -8.0, 0.0, 8.0, -1.0], dtype=np.float64) / (12.0 * float(h))
    elif order == 2:
        mat = _second_derivative_matrix(n, h)
        weights = np.asarray([-1.0, 16.0, -30.0, 16.0, -1.0], dtype=np.float64) / (12.0 * float(h) ** 2)
    else:
        raise ValueError(f"unsupported derivative order: {order}")
    for i in range(margin, n - margin):
        mat[i, :] = 0.0
        mat[i, i - 2 : i + 3] = weights
    return mat


def _block_row(blocks: Sequence[np.ndarray]) -> np.ndarray:
    return np.concatenate(blocks, axis=1)


def _weighted_gram(matrix: np.ndarray, weights: np.ndarray, sign: float) -> np.ndarray:
    return sign * matrix.conj().T @ (weights[:, None] * matrix)


def _gamma_port(a: float = FROZEN_A, c_s: float = CS_CONST) -> float:
    return 4.0 * a**5 / (27.0 * c_s**5)


def _background_arrays_on_transfer_grid(
    background: Mapping[str, Any],
    r_centers: np.ndarray,
    w_centers: np.ndarray,
) -> dict[str, np.ndarray]:
    r_old = np.asarray(background["grid"]["r_centers"], dtype=np.float64)
    w_old = np.asarray(background["grid"]["w_centers"], dtype=np.float64)

    def field2(name: str) -> np.ndarray:
        return _interp2_clamped(
            r_old,
            w_old,
            np.asarray(background["fields"][name], dtype=np.float64),
            r_centers,
            w_centers,
        )

    z_old = np.asarray(background.get("derived", {}).get("Z_w", np.ones_like(w_old)), dtype=np.float64)
    r0_old = np.asarray(background["fields"]["R0_w"], dtype=np.float64)
    taper = _radial_exterior_taper(r_old, r_centers)

    def compact_field(name: str) -> np.ndarray:
        return taper[:, None] * field2(name)

    return {
        "psi_R0": compact_field("psi_R0"),
        "psi_I0": compact_field("psi_I0"),
        "A_00": compact_field("A_00"),
        "A_r0": compact_field("A_r0"),
        "A_w0": compact_field("A_w0"),
        "Z_w": _interp1_clamped(w_old, z_old, w_centers),
        "R0_w": _interp1_clamped(w_old, r0_old, w_centers),
        "radial_exterior_taper": taper,
    }


def make_transfer_grid(
    background: Mapping[str, Any],
    *,
    nr: int,
    nw: int,
    radial_scale: float = DEFAULT_FINAL_TRUNCATION,
) -> dict[str, Any]:
    if nr < 3 or nw < 3:
        raise ValueError("Maxwell transfer grid must be at least 3x3")
    geometry = background["geometry"]
    a = float(geometry["a"])
    outer = float(background["grid"]["r_max"]) * float(radial_scale)
    r_min = 0.22 * a
    if outer <= r_min:
        raise ValueError(f"radial truncation outer={outer} is not beyond r_min={r_min}")
    L = float(geometry["L"])
    dr = (outer - r_min) / (nr + 1)
    dw = L / (nw + 1)
    r1 = np.asarray([r_min + (i + 1) * dr for i in range(nr)], dtype=np.float64)
    w1 = np.asarray([(j + 1) * dw for j in range(nw)], dtype=np.float64)
    r = np.asarray([r1[i] for j in range(nw) for i in range(nr)], dtype=np.float64)
    w = np.asarray([w1[j] for j in range(nw) for i in range(nr)], dtype=np.float64)
    arrays = _background_arrays_on_transfer_grid(background, r1, w1)
    z2d = np.tile(np.asarray(arrays["Z_w"], dtype=np.float64)[None, :], (nr, 1))
    z = _flatten_w_major(z2d)
    quad = dr * dw * r**2 * z
    return {
        "nr": nr,
        "nw": nw,
        "n": nr * nw,
        "dr": float(dr),
        "dw": float(dw),
        "r1": r1,
        "w1": w1,
        "r": r,
        "w": w,
        "z": z,
        "w0": quad,
        "r_min": float(r_min),
        "r_outer": float(outer),
        "radial_scale": float(radial_scale),
        "h": float(max(dr, dw)),
        "arrays": arrays,
    }


def make_staggered_transfer_grid(
    background: Mapping[str, Any],
    *,
    nr: int,
    nw: int,
    radial_scale: float = DEFAULT_FINAL_TRUNCATION,
) -> dict[str, Any]:
    if nr < 5 or nw < 5:
        raise ValueError("staggered Maxwell transfer grid must be at least 5x5")
    geometry = background["geometry"]
    a = float(geometry["a"])
    outer = float(background["grid"]["r_max"]) * float(radial_scale)
    r_min = 0.22 * a
    if outer <= r_min:
        raise ValueError(f"radial truncation outer={outer} is not beyond r_min={r_min}")
    L = float(geometry["L"])
    dr = (outer - r_min) / nr
    dw = L / nw
    r1 = np.asarray([r_min + (i + 0.5) * dr for i in range(nr)], dtype=np.float64)
    w1 = np.asarray([(j + 0.5) * dw for j in range(nw)], dtype=np.float64)
    r = np.asarray([r1[i] for j in range(nw) for i in range(nr)], dtype=np.float64)
    w = np.asarray([w1[j] for j in range(nw) for i in range(nr)], dtype=np.float64)
    arrays = _background_arrays_on_transfer_grid(background, r1, w1)
    z2d = np.tile(np.asarray(arrays["Z_w"], dtype=np.float64)[None, :], (nr, 1))
    z = _flatten_w_major(z2d)
    quad = dr * dw * r**2 * z
    return {
        "nr": nr,
        "nw": nw,
        "n": nr * nw,
        "dr": float(dr),
        "dw": float(dw),
        "r1": r1,
        "w1": w1,
        "r": r,
        "w": w,
        "z": z,
        "w0": quad,
        "r_min": float(r_min),
        "r_outer": float(outer),
        "radial_scale": float(radial_scale),
        "h": float(max(dr, dw)),
        "arrays": arrays,
        "mesh": "cell-centered staggered",
    }


def _boundary_matrix(grid: Mapping[str, Any], fields: Mapping[str, np.ndarray]) -> np.ndarray:
    nr = int(grid["nr"])
    nw = int(grid["nw"])
    n = int(grid["n"])
    w_exit_scalar = np.zeros(n, dtype=np.float64)
    w_exit_angular = np.zeros(n, dtype=np.float64)
    r_outer_scalar = np.zeros(n, dtype=np.float64)
    r_outer_angular = np.zeros(n, dtype=np.float64)
    r = np.asarray(grid["r"], dtype=np.float64)
    z = np.asarray(grid["z"], dtype=np.float64)
    dr = float(grid["dr"])
    dw = float(grid["dw"])
    for j in range(nw):
        for i in range(nr):
            cell = _cell_index(i, j, nr)
            if j == nw - 1:
                weight = dr * r[cell] ** 2 * z[cell]
                w_exit_scalar[cell] += weight
                w_exit_angular[cell] += L_HARM * weight
            if i == nr - 1:
                weight = dw * r[cell] ** 2 * z[cell]
                r_outer_scalar[cell] += weight
                r_outer_angular[cell] += L_HARM * weight
    return (
        fields["E_r"].conj().T @ (_diag(w_exit_scalar) @ fields["E_r"])
        + fields["E_E"].conj().T @ (_diag(w_exit_angular + r_outer_angular) @ fields["E_E"])
        + fields["E_B"].conj().T @ (_diag(w_exit_angular + r_outer_angular) @ fields["E_B"])
        + fields["E_w"].conj().T @ (_diag(r_outer_scalar) @ fields["E_w"])
    )


def assemble_vsh_maxwell_system(
    background: Mapping[str, Any],
    *,
    nr: int,
    nw: int,
    omega: float,
    radial_scale: float = DEFAULT_FINAL_TRUNCATION,
    discretization: str = "primary_second_order",
) -> dict[str, Any]:
    if discretization == "primary_second_order":
        grid = make_transfer_grid(background, nr=nr, nw=nw, radial_scale=radial_scale)
        dr1 = _first_derivative_matrix(nr, float(grid["dr"]))
        dw1 = _first_derivative_matrix(nw, float(grid["dw"]))
        drr1 = _second_derivative_matrix(nr, float(grid["dr"]))
        dww1 = _second_derivative_matrix(nw, float(grid["dw"]))
    elif discretization == "staggered_second_order":
        grid = make_staggered_transfer_grid(background, nr=nr, nw=nw, radial_scale=radial_scale)
        dr1 = _first_derivative_matrix(nr, float(grid["dr"]))
        dw1 = _first_derivative_matrix(nw, float(grid["dw"]))
        drr1 = _second_derivative_matrix(nr, float(grid["dr"]))
        dww1 = _second_derivative_matrix(nw, float(grid["dw"]))
    elif discretization == "staggered_high_order":
        grid = make_staggered_transfer_grid(background, nr=nr, nw=nw, radial_scale=radial_scale)
        dr1 = _high_order_derivative_matrix(nr, float(grid["dr"]), 1)
        dw1 = _high_order_derivative_matrix(nw, float(grid["dw"]), 1)
        drr1 = _high_order_derivative_matrix(nr, float(grid["dr"]), 2)
        dww1 = _high_order_derivative_matrix(nw, float(grid["dw"]), 2)
    else:
        raise ValueError(f"unknown Maxwell transfer discretization: {discretization}")
    n = int(grid["n"])
    total = 5 * n
    id_r = np.eye(nr, dtype=np.complex128)
    id_w = np.eye(nw, dtype=np.complex128)
    id_n = np.eye(n, dtype=np.complex128)
    z_n = np.zeros((n, n), dtype=np.complex128)
    dr_mat = np.kron(id_w, dr1)
    dw_mat = np.kron(dw1, id_r)
    drr_mat = np.kron(id_w, drr1)
    dww_mat = np.kron(dww1, id_r)

    r = np.asarray(grid["r"], dtype=np.float64)
    inv_r = 1.0 / r
    inv_r2 = inv_r**2
    r_diag = _diag(r)
    r2_diag = _diag(r**2)
    inv_r_diag = _diag(inv_r)
    inv_r2_diag = _diag(inv_r2)
    w0 = np.asarray(grid["w0"], dtype=np.float64)
    w0_diag = _diag(w0)
    lw0_diag = _diag(L_HARM * w0)
    w_lane_values = np.concatenate([w0, w0, L_HARM * w0, L_HARM * w0, w0])
    w_lane = _diag(w_lane_values)
    w_lane_inv = _diag(1.0 / w_lane_values)
    div_r = inv_r2_diag @ dr_mat @ r2_diag

    iomega = 1j * float(omega)
    m_g = _block_row([iomega * id_n, div_r, -L_HARM * inv_r_diag, z_n, dw_mat])
    m_er = _block_row([-dr_mat, iomega * id_n, z_n, z_n, z_n])
    m_ee = _block_row([-inv_r_diag, z_n, iomega * id_n, z_n, z_n])
    m_eb = _block_row([z_n, z_n, z_n, iomega * id_n, z_n])
    m_ew = _block_row([-dw_mat, z_n, z_n, z_n, iomega * id_n])
    m_cr = _block_row([z_n, -dw_mat, z_n, z_n, dr_mat])
    m_ce = _block_row([z_n, z_n, -dw_mat, z_n, inv_r_diag])
    m_cb = _block_row([z_n, z_n, z_n, -dw_mat, z_n])
    m_bb = _block_row([z_n, -inv_r_diag, inv_r_diag @ dr_mat @ r_diag, z_n, z_n])
    m_br = _block_row([z_n, z_n, z_n, -L_HARM * inv_r_diag, z_n])
    m_be = _block_row([z_n, z_n, z_n, -inv_r_diag @ dr_mat @ r_diag, z_n])

    k_action = (
        _weighted_gram(m_er, w0, 1.0)
        + _weighted_gram(m_ee, L_HARM * w0, 1.0)
        + _weighted_gram(m_eb, L_HARM * w0, 1.0)
        + _weighted_gram(m_ew, w0, 1.0)
        - _weighted_gram(m_cr, w0, 1.0)
        - _weighted_gram(m_ce, L_HARM * w0, 1.0)
        - _weighted_gram(m_cb, L_HARM * w0, 1.0)
        - _weighted_gram(m_bb, L_HARM * w0, 1.0)
        - _weighted_gram(m_br, w0, 1.0)
        - _weighted_gram(m_be, L_HARM * w0, 1.0)
        - (1.0 / XI_GAUGE) * _weighted_gram(m_g, w0, 1.0)
    )
    field_matrices = {
        "G": m_g,
        "E_r": m_er,
        "E_E": m_ee,
        "E_B": m_eb,
        "E_w": m_ew,
        "C_r": m_cr,
        "C_E": m_ce,
        "C_B": m_cb,
        "B_B": m_bb,
        "B_r": m_br,
        "B_E": m_be,
    }
    boundary_b = _boundary_matrix(grid, field_matrices)
    anchor = np.zeros(total, dtype=np.float64)
    anchor_cell = int(math.ceil(n / 2.0)) - 1
    anchor[anchor_cell] = 1.0e-8 * w0[anchor_cell]
    k_cons = k_action + _diag(anchor)
    gamma = _gamma_port(float(background["geometry"]["a"]), CS_CONST)
    k_ret = k_cons + 1j * gamma * float(omega) ** 5 * boundary_b
    return {
        **grid,
        "omega": float(omega),
        "total": total,
        "Dr": dr_mat,
        "Dw": dw_mat,
        "Drr": drr_mat,
        "Dww": dww_mat,
        "W0": w0_diag,
        "WLane": w_lane,
        "WLaneValues": w_lane_values,
        "WLaneInv": w_lane_inv,
        "KCons": k_cons,
        "KRet": k_ret,
        "ACons": w_lane_inv @ k_cons,
        "ARet": w_lane_inv @ k_ret,
        "BoundaryB": boundary_b,
        "AnchorMatrix": _diag(anchor),
        "ScalarGaugeAnchor": {"cell": int(anchor_cell), "strength": 1.0e-8},
        "FieldMatrices": field_matrices,
        "discretization": discretization,
    }


def _bdg_response_on_w(
    bdg_packet: Mapping[str, Any],
    w_values: np.ndarray,
    *,
    omega: float,
    mode_count: int = DEFAULT_MODE_COUNT,
) -> np.ndarray:
    profile_points = np.asarray(bdg_packet["grid"]["profile_points"], dtype=np.float64)
    response = np.zeros_like(w_values, dtype=np.float64)
    modes = list(bdg_packet["bdg_modes"])[:mode_count]
    for mode in modes:
        varpi = float(mode["varpi"])
        denom = varpi**2 - float(omega) ** 2
        if abs(denom) <= 1.0e-14:
            denom = math.copysign(1.0e-14, denom)
        profile = np.asarray(mode.get("profile_values", mode["profile"]), dtype=np.float64)
        response += float(mode["coupling"]) * _interp1_clamped(profile_points, profile, w_values) / denom
    return response


def build_matter_source(
    sys: Mapping[str, Any],
    background: Mapping[str, Any],
    bdg_packet: Mapping[str, Any],
    *,
    mode_count: int = DEFAULT_MODE_COUNT,
) -> dict[str, Any]:
    nr = int(sys["nr"])
    nw = int(sys["nw"])
    n = int(sys["n"])
    arrays = sys["arrays"]
    constants = background["constants"]
    gauge_charge = float(constants["gauge_charge"])
    particle_mass = float(constants["particle_mass"])
    omega = float(sys["omega"])
    r = np.asarray(sys["r"], dtype=np.float64)
    w = np.asarray(sys["w"], dtype=np.float64)
    psi_r = _flatten_w_major(np.asarray(arrays["psi_R0"], dtype=np.float64))
    psi_i = _flatten_w_major(np.asarray(arrays["psi_I0"], dtype=np.float64))
    a0 = _flatten_w_major(np.asarray(arrays["A_00"], dtype=np.float64))
    ar0 = _flatten_w_major(np.asarray(arrays["A_r0"], dtype=np.float64))
    aw0 = _flatten_w_major(np.asarray(arrays["A_w0"], dtype=np.float64))
    rho = psi_r**2 + psi_i**2
    axial_response = _bdg_response_on_w(
        bdg_packet,
        np.asarray(sys["w1"], dtype=np.float64),
        omega=omega,
        mode_count=mode_count,
    )
    z = np.asarray(sys["z"], dtype=np.float64)
    radial_shape = np.zeros(n, dtype=np.float64)
    for j in range(nw):
        cells = np.asarray([_cell_index(i, j, nr) for i in range(nr)])
        raw = rho[cells].copy()
        weight = float(sys["dr"]) * r[cells] ** 2 * z[cells]
        norm = float(np.sum(weight * raw))
        if not math.isfinite(norm) or norm <= 1.0e-300:
            raw = np.ones_like(raw)
            norm = float(np.sum(weight * raw))
        radial_shape[cells] = raw / norm
    mean_a0 = float(np.average(a0, weights=np.maximum(np.asarray(sys["w0"], dtype=np.float64), 1.0e-300)))
    scalar_mod = 1.0 + gauge_charge * (a0 - mean_a0)
    delta_rho = np.zeros(n, dtype=np.float64)
    for j in range(nw):
        cells = np.asarray([_cell_index(i, j, nr) for i in range(nr)])
        delta_rho[cells] = axial_response[j] * radial_shape[cells] * scalar_mod[cells]
    j_r = -(gauge_charge / particle_mass) * ar0 * delta_rho
    j_e = -1j * omega * r * delta_rho / L_HARM
    j_w = -(gauge_charge / particle_mass) * aw0 * delta_rho
    zero_b = np.zeros(n, dtype=np.complex128)
    j_vec = np.concatenate(
        [
            gauge_charge * delta_rho,
            gauge_charge * j_r,
            gauge_charge * j_e,
            zero_b,
            gauge_charge * j_w,
        ]
    ).astype(np.complex128)
    source_norm = math.sqrt(
        max(0.0, float(np.real(j_vec.conj() @ (np.asarray(sys["WLaneValues"]) * j_vec))))
    )
    continuity = (
        -1j * omega * delta_rho
        + np.asarray(sys["Dr"] @ j_r, dtype=np.complex128)
        + np.asarray(sys["Dw"] @ j_w, dtype=np.complex128)
        - L_HARM * j_e / r
    )
    continuity_residual = float(
        np.linalg.norm(continuity)
        / max(np.linalg.norm(delta_rho) + np.linalg.norm(j_r) + np.linalg.norm(j_w) + np.linalg.norm(j_e), 1.0e-30)
    )
    return {
        "jVec": j_vec,
        "deltaRho": delta_rho,
        "jr": j_r,
        "jE": j_e,
        "jB": zero_b,
        "jw": j_w,
        "sourceNorm": source_norm,
        "continuityResidual": continuity_residual,
        "usedB2aBdgModes": int(min(mode_count, len(bdg_packet["bdg_modes"]))),
        "usedReturnSourceSEtaA": False,
    }


def self_energy_solve(sys: Mapping[str, Any], j_vec: np.ndarray, mode: str) -> dict[str, Any]:
    op = np.asarray(sys["ARet" if mode == "retarded" else "ACons"], dtype=np.complex128)
    rhs = float(sys.get("mu0", 1.0)) * np.asarray(j_vec, dtype=np.complex128)
    solution = dense_solve(op, rhs, assume_a="gen", check_finite=False)
    residual = float(np.linalg.norm(op @ solution - rhs) / max(np.linalg.norm(rhs), 1.0e-30))
    w_lane_values = np.asarray(sys["WLaneValues"], dtype=np.float64)
    sigma = complex(np.vdot(j_vec, w_lane_values * solution))
    boundary_quadratic = float(np.real(np.vdot(solution, np.asarray(sys["BoundaryB"]) @ solution)))
    gamma = _gamma_port(FROZEN_A, CS_CONST)
    flux = gamma * float(sys["omega"]) ** 5 * boundary_quadratic / 2.0
    return {
        "solution": solution,
        "sigma": sigma,
        "residual": residual,
        "boundaryQuadratic": boundary_quadratic,
        "outgoingFlux": float(flux),
    }


def _fit_even_coefficients(omegas: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, float]:
    x = np.asarray(omegas, dtype=np.float64) ** 2
    mat = np.column_stack([np.ones_like(x), x, x**2])
    coeff, *_ = np.linalg.lstsq(mat, np.asarray(values, dtype=np.float64), rcond=None)
    cond = float(np.linalg.cond(mat))
    return coeff.astype(np.float64), cond


def _fit_even_coefficients_scaled(
    omegas: np.ndarray,
    values: np.ndarray,
    *,
    window: float,
) -> tuple[np.ndarray, float, float]:
    scale = max(float(window), float(np.max(np.abs(omegas))), 1.0e-300)
    x = (np.asarray(omegas, dtype=np.float64) / scale) ** 2
    mat = np.column_stack([np.ones_like(x), x, x**2])
    scaled_coeff, *_ = np.linalg.lstsq(mat, np.asarray(values, dtype=np.float64), rcond=None)
    residual = np.asarray(values, dtype=np.float64) - mat @ scaled_coeff
    coeff = np.asarray(
        [
            scaled_coeff[0],
            scaled_coeff[1] / scale**2,
            scaled_coeff[2] / scale**4,
        ],
        dtype=np.float64,
    )
    cond = float(np.linalg.cond(mat))
    rms = float(np.linalg.norm(residual) / max(math.sqrt(float(residual.size)), 1.0))
    return coeff, cond, rms


def transfer_for_grid(
    background: Mapping[str, Any],
    bdg_packet: Mapping[str, Any],
    *,
    nr: int,
    nw: int,
    window: float,
    radial_scale: float = DEFAULT_FINAL_TRUNCATION,
    mode_count: int = DEFAULT_MODE_COUNT,
    engine: str = "python_primary_second_order",
    discretization: str = "primary_second_order",
) -> dict[str, Any]:
    fractions = np.asarray([0.55, 0.70, 0.85, 1.00, 1.15], dtype=np.float64)
    omegas = float(window) * fractions
    rows: list[dict[str, Any]] = []
    for omega in omegas:
        sys = assemble_vsh_maxwell_system(
            background,
            nr=nr,
            nw=nw,
            omega=float(omega),
            radial_scale=radial_scale,
            discretization=discretization,
        )
        matter = build_matter_source(sys, background, bdg_packet, mode_count=mode_count)
        cons = self_energy_solve(sys, matter["jVec"], "conservative")
        ret = self_energy_solve(sys, matter["jVec"], "retarded")
        gamma = _gamma_port(float(background["geometry"]["a"]), CS_CONST)
        rows.append(
            {
                "omega": float(omega),
                "sigmaCons": {"re": cons["sigma"].real, "im": cons["sigma"].imag},
                "sigmaRet": {"re": ret["sigma"].real, "im": ret["sigma"].imag},
                "minusImSigmaRet": float(-ret["sigma"].imag),
                "NofOmega": float(-ret["sigma"].imag / (gamma * float(omega) ** 5)),
                "consResidual": cons["residual"],
                "retResidual": ret["residual"],
                "continuityResidual": matter["continuityResidual"],
                "outgoingFlux": ret["outgoingFlux"],
                "boundaryQuadratic": ret["boundaryQuadratic"],
                "sourceNorm": matter["sourceNorm"],
            }
        )
    z_coeff, z_cond, z_fit_rms = _fit_even_coefficients_scaled(
        omegas,
        np.asarray([float(row["sigmaCons"]["re"]) for row in rows], dtype=np.float64),
        window=window,
    )
    n_coeff, n_cond, n_fit_rms = _fit_even_coefficients_scaled(
        omegas,
        np.asarray([float(row["NofOmega"]) for row in rows], dtype=np.float64),
        window=window,
    )
    packet: dict[str, Any] = {
        "schema": "stage1_patha_b2b_maxwell_transfer/v1",
        "engine": engine,
        "background_content_hash": background["content_hash"],
        "bdg_packet_content_hash": bdg_packet["content_hash"],
        "tau": float(background["constants"]["tau"]),
        "geometry": {
            "a": float(background["geometry"]["a"]),
            "L": float(background["geometry"]["L"]),
            "radial_scale": float(radial_scale),
        },
        "constants": {
            "gauge_charge": float(background["constants"]["gauge_charge"]),
            "particle_mass": float(background["constants"]["particle_mass"]),
            "mu0": float(background["constants"]["mu0"]),
            "tau": float(background["constants"]["tau"]),
            "Gamma_port": _gamma_port(float(background["geometry"]["a"]), CS_CONST),
            "Gamma_port_formula": "4*a^5/(27*c_s^5) with c_s=1",
        },
        "grid": {
            "nr": int(nr),
            "nw": int(nw),
            "window": float(window),
            "radial_scale": float(radial_scale),
            "mode_count": int(min(mode_count, len(bdg_packet["bdg_modes"]))),
            "discretization": discretization,
        },
        "coefficients": {
            "Z0": float(z_coeff[0]),
            "Z2": float(z_coeff[1]),
            "Z4": float(z_coeff[2]),
            "N0": float(n_coeff[0]),
            "N2": float(n_coeff[1]),
            "N4": float(n_coeff[2]),
        },
        "rows": rows,
        "diagnostics": {
            "operator_matrix_dimension": int(5 * nr * nw),
            "max_green_residual": float(max(max(row["consResidual"], row["retResidual"]) for row in rows)),
            "max_continuity_residual": float(max(row["continuityResidual"] for row in rows)),
            "min_minus_im_sigma_ret": float(min(row["minusImSigmaRet"] for row in rows)),
            "min_outgoing_flux": float(min(row["outgoingFlux"] for row in rows)),
            "fit_condition_number_Z": z_cond,
            "fit_condition_number_N": n_cond,
            "fit_basis": "scaled even powers: 1,(omega/window)^2,(omega/window)^4; coefficients converted back",
            "fit_rms_Z": z_fit_rms,
            "fit_rms_N": n_fit_rms,
            "source_chain": (
                "decision-05 D3 Frechet source over B2a BdG-mode response -> axial density response -> Path-A A0-modulated "
                "charge/current source -> five-lane VSH Maxwell Green solve"
            ),
            "operator_term_map": _operator_term_map(),
        },
    }
    packet["content_hash"] = _full_stable_hash(packet)
    return packet


def _operator_term_map() -> list[dict[str, str]]:
    return [
        {
            "term": "five-lane VSH Maxwell weak form",
            "implementation": "lanes phi,a_r,a_E,a_B,a_w with L=l(l+1)=6 and H=Z gauge weight",
            "source": "software/stage1_solver/decisions/05_mixed_maxwell_spike_design.md D1-D2",
        },
        {
            "term": "open DtN loss",
            "implementation": "K_ret=K_cons+i Gamma_port omega^5 B_tan",
            "source": "notes/moving_throat_pde_program_compact.md:2559-2568",
        },
        {
            "term": "basis-invariant self energy",
            "implementation": "Sigma=<j,G j>; Z_n from Re Sigma_cons; N_n from -Im Sigma_ret/(Gamma_port omega^5)",
            "source": "software/stage1_solver/decisions/05_mixed_maxwell_spike_design.md D4",
        },
        {
            "term": "wall/BdG source",
            "implementation": (
                "canonical decision-05 D3 Frechet source driven by the B2a BdG-mode response; "
                "B2a overlaps c_j carry the chi coupling, no live U/W ports"
            ),
            "source": "software/stage1_solver/directives/pathA_11_chunk_b2b_maxwell_transfer.md requirements 3-4",
        },
    ]


def solve_python_command(
    *,
    background_path: Path,
    bdg_path: Path,
    nr: int,
    nw: int,
    window: float,
    radial_scale: float,
    mode_count: int,
    run_root: Path,
    out_path: Path | None = None,
    engine: str = "python_primary_second_order",
    discretization: str = "primary_second_order",
) -> tuple[Path, dict[str, Any]]:
    background = _load_json(background_path)
    bdg_packet = _load_json(bdg_path)
    packet = transfer_for_grid(
        background,
        bdg_packet,
        nr=nr,
        nw=nw,
        window=window,
        radial_scale=radial_scale,
        mode_count=mode_count,
        engine=engine,
        discretization=discretization,
    )
    if out_path is None:
        out_dir, prefix = _engine_folder_and_prefix(engine)
        out_dir = run_root / out_dir
        tau_label = _format_tau(float(background["constants"]["tau"]))
        out_path = out_dir / (
            f"patha_b2b_{prefix}"
            f"_tau_{tau_label}_nr_{nr}_nw_{nw}_w_{window:.6g}_rs_{radial_scale:.6g}.json"
        )
    _write_json(out_path, packet)
    return out_path, packet


def _default_background_path(b2a_run_root: Path, tau: float) -> Path:
    return b2a_run_root / "backgrounds" / f"patha_b2a_closed_background_tau_{_format_tau(tau)}_latest.json"


def _default_b2a_bdg_path(b2a_run_root: Path, tau: float, grid: tuple[int, int] = (10, 10)) -> Path:
    if abs(float(tau) - 1.0) <= 1.0e-15:
        bundle = b2a_run_root / "bundles" / "patha_b2a_validated_bdg_bundle_tau_1.json"
        if bundle.exists():
            return bundle
    return b2a_run_root / "python" / f"patha_b2a_python_tau_{_format_tau(tau)}_nr_{grid[0]}_nw_{grid[1]}.json"


def _engine_folder_and_prefix(engine: str) -> tuple[str, str]:
    if engine in {"python", "primary", "python_primary_second_order"} or engine.startswith("python_primary"):
        return "python", "python"
    if engine in {"independent", "staggered", "python_staggered_high_order"} or engine.startswith("python_staggered"):
        return "independent", "independent"
    if engine in {"mathematica", "mma", "mathematica_transcription"} or engine.startswith("mathematica"):
        return "mathematica", "mma"
    if engine.startswith("python"):
        return "python", "python"
    return "independent", "independent"


def _default_packet_path(
    run_root: Path,
    engine: str,
    tau: float,
    nr: int,
    nw: int,
    window: float,
    radial_scale: float,
) -> Path:
    folder, prefix = _engine_folder_and_prefix(engine)
    return run_root / folder / (
        f"patha_b2b_{prefix}_tau_{_format_tau(tau)}_nr_{nr}_nw_{nw}"
        f"_w_{window:.6g}_rs_{radial_scale:.6g}.json"
    )


def _radial_grid_for_scale(
    background: Mapping[str, Any],
    *,
    final_grid: tuple[int, int],
    final_truncation: float,
    radial_scale: float,
) -> tuple[int, int]:
    geometry = background["geometry"]
    r_min = 0.22 * float(geometry["a"])
    r_background_max = float(background["grid"]["r_max"])
    final_outer = r_background_max * float(final_truncation)
    final_dr = (final_outer - r_min) / (int(final_grid[0]) + 1)
    outer = r_background_max * float(radial_scale)
    nr = int(round((outer - r_min) / final_dr - 1.0))
    return max(5, nr), int(final_grid[1])


def _coeff_values(packet: Mapping[str, Any]) -> list[float]:
    return [float(packet["coefficients"][key]) for key in COEFF_KEYS]


def _max_diffs(left: Sequence[float], right: Sequence[float]) -> dict[str, float]:
    left_arr = np.asarray(left, dtype=np.float64)
    right_arr = np.asarray(right, dtype=np.float64)
    diff = np.abs(left_arr - right_arr)
    denom = np.maximum(np.abs(right_arr), 1.0e-300)
    return {
        "max_abs": float(np.max(diff)) if diff.size else 0.0,
        "max_rel": float(np.max(diff / denom)) if diff.size else 0.0,
    }


def compare_engine_packets(python_packet: Mapping[str, Any], mma_packet: Mapping[str, Any]) -> dict[str, Any]:
    per_coeff: dict[str, dict[str, float]] = {}
    for key in COEFF_KEYS:
        per_coeff[key] = _max_diffs(
            [float(python_packet["coefficients"][key])],
            [float(mma_packet["coefficients"][key])],
        )
    all_diff = _max_diffs(_coeff_values(python_packet), _coeff_values(mma_packet))
    return {"coefficients": all_diff, "per_coefficient": per_coeff}


def _strict_diff_pass(diff: Mapping[str, float], *, abs_tol: float, rel_tol: float) -> bool:
    return bool(float(diff["max_abs"]) <= abs_tol and float(diff["max_rel"]) <= rel_tol)


def _dual_engine_gate(dual_rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    max_abs = max(float(row["coefficients"]["max_abs"]) for row in dual_rows)
    max_rel = max(float(row["coefficients"]["max_rel"]) for row in dual_rows)
    per_coeff: dict[str, dict[str, float]] = {}
    for key in COEFF_KEYS:
        per_coeff[key] = {
            "max_abs": max(float(row["per_coefficient"][key]["max_abs"]) for row in dual_rows),
            "max_rel": max(float(row["per_coefficient"][key]["max_rel"]) for row in dual_rows),
            "abs": DUAL_ENGINE_TOLERANCES["abs"],
            "rel": DUAL_ENGINE_TOLERANCES["rel"],
        }
        per_coeff[key]["passed"] = _strict_diff_pass(
            per_coeff[key],
            abs_tol=DUAL_ENGINE_TOLERANCES["abs"],
            rel_tol=DUAL_ENGINE_TOLERANCES["rel"],
        )
    return {
        "passed": bool(
            max_abs <= DUAL_ENGINE_TOLERANCES["abs"]
            and max_rel <= DUAL_ENGINE_TOLERANCES["rel"]
            and all(row["passed"] for row in per_coeff.values())
        ),
        "max_abs": max_abs,
        "max_rel": max_rel,
        "per_coefficient": per_coeff,
        "criteria": (
            "AND gate: max coefficient abs <= "
            f"{DUAL_ENGINE_TOLERANCES['abs']:.1e} and rel <= {DUAL_ENGINE_TOLERANCES['rel']:.1e}"
        ),
    }


def basis_invariance_check(
    background: Mapping[str, Any],
    bdg_packet: Mapping[str, Any],
    *,
    broken: bool = False,
    inject_port_leak: bool = False,
) -> dict[str, Any]:
    omega = 0.039
    sys = assemble_vsh_maxwell_system(background, nr=4, nw=4, omega=omega, radial_scale=0.9)
    matter = build_matter_source(sys, background, bdg_packet, mode_count=min(8, len(bdg_packet["bdg_modes"])))
    j_vec = matter["jVec"]
    cons = self_energy_solve(sys, j_vec, "conservative")
    sigma_base = cons["sigma"]
    w_vals = np.asarray(sys["WLaneValues"], dtype=np.float64)
    w_sqrt = np.sqrt(w_vals)
    w_inv_sqrt = 1.0 / w_sqrt
    a_tilde = (w_sqrt[:, None] * np.asarray(sys["ACons"])) * w_inv_sqrt[None, :]
    b_tilde = (w_inv_sqrt[:, None] * np.asarray(sys["BoundaryB"])) * w_inv_sqrt[None, :]
    s_tilde = w_sqrt * j_vec
    u_base = dense_solve(a_tilde, s_tilde, assume_a="gen", check_finite=False)
    n_base = float(np.real(np.vdot(u_base, b_tilde @ u_base)))
    deltas = []
    tested_deltas = []
    dim = s_tilde.size
    n_cell = int(sys["n"])
    port = np.zeros(dim, dtype=np.complex128)
    port[:n_cell] = 1.0 / math.sqrt(float(n_cell))
    port[2 * n_cell : 3 * n_cell] = 0.35 / math.sqrt(float(n_cell))
    leaked_base = float(abs(np.vdot(port, u_base)) ** 2)
    tested_base = leaked_base if inject_port_leak else float(sigma_base.real)
    leaked_rebased_values = []
    lane_rotations = []
    for theta_scalar, theta_transverse in ((0.31, -0.22), (-0.17, 0.39), (0.24, 0.28)):
        lane = np.eye(5, dtype=np.complex128)
        cs = math.cos(theta_scalar)
        ss = math.sin(theta_scalar)
        lane[np.ix_([0, 4], [0, 4])] = np.asarray([[cs, -ss], [ss, cs]], dtype=np.complex128)
        ct = math.cos(theta_transverse)
        st = math.sin(theta_transverse)
        lane[np.ix_([2, 3], [2, 3])] = np.asarray([[ct, -st], [st, ct]], dtype=np.complex128)
        # Mild scalar/radial gauge-lane rebasis after the orthogonal rotations.
        lane[1, 0] += 0.08
        lane[0, 1] -= 0.08
        q, _ = np.linalg.qr(lane)
        lane_rotations.append(q)
    for idx, lane in enumerate(lane_rotations):
        transform = np.kron(lane, np.eye(n_cell, dtype=np.complex128))
        rotated_a = transform.conj().T @ a_tilde @ transform
        rotated_s = transform.conj().T @ s_tilde
        rotated_b = transform.conj().T @ b_tilde @ transform
        if broken:
            rotated_s = s_tilde
        sol = dense_solve(rotated_a, rotated_s, assume_a="gen", check_finite=False)
        sig = complex(np.vdot(rotated_s, sol))
        n_rot = float(np.real(np.vdot(sol, rotated_b @ sol)))
        leaked_wrong = float(abs(np.vdot(port, sol)) ** 2)
        tested_value = leaked_wrong if inject_port_leak else float(sig.real)
        leaked_rebased_values.append(leaked_wrong)
        tested_deltas.append(float(abs(tested_value - tested_base) / max(abs(tested_base), 1.0e-300)))
        deltas.append(
            {
                "rotation": idx + 1,
                "relative_Z_delta": float(abs(sig - sigma_base) / max(abs(sigma_base), 1.0e-300)),
                "relative_N_delta": float(abs(n_rot - n_base) / max(abs(n_base), 1.0e-300)),
                "wrong_fixed_port_relative_delta": float(
                    abs(leaked_wrong - leaked_base) / max(abs(leaked_base), 1.0e-300)
                ),
            }
        )
    max_z = max(row["relative_Z_delta"] for row in deltas)
    max_n = max(row["relative_N_delta"] for row in deltas)
    max_port_leak = max(row["wrong_fixed_port_relative_delta"] for row in deltas)
    max_tested_delta = max(tested_deltas)
    invariant_passed = bool(max_z <= 1.0e-10 and max_n <= 1.0e-10)
    leak_failed_as_expected = bool(max_port_leak >= 1.0e-3)
    tested_extraction_passed = bool(max_tested_delta <= 1.0e-10)
    branch_difference = float(abs(tested_base - float(sigma_base.real)) / max(abs(float(sigma_base.real)), 1.0e-300))
    return {
        "grid": "4x4",
        "omega": omega,
        "broken_control": bool(broken or inject_port_leak),
        "extraction_mode": "fixed_port_leak" if inject_port_leak else "basis_invariant_quadratic",
        "base_Z": float(sigma_base.real),
        "base_N_boundary_quadratic": n_base,
        "tested_base_value": tested_base,
        "rotations": deltas,
        "max_relative_Z_delta": max_z,
        "max_relative_N_delta": max_n,
        "max_wrong_fixed_port_relative_delta": max_port_leak,
        "max_tested_extraction_relative_delta": max_tested_delta,
        "tested_extraction_passed": tested_extraction_passed,
        "branch_difference": {
            "relative_to_invariant_Z": branch_difference,
            "different_from_invariant": bool(branch_difference >= 1.0e-3),
        },
        "port_leak_negative_control": {
            "failed_as_expected": leak_failed_as_expected,
            "relative_delta": max_port_leak,
            "wrong_answer": "a posited U/W-style fixed-port extraction that is not rebased with the physical lanes",
        },
        "passed": bool(invariant_passed and leak_failed_as_expected and tested_extraction_passed),
        "criteria": (
            "lane-local gauge/transverse rebasis in W-normalized coordinates leaves Sigma and boundary N invariant; "
            "a fixed-port U/W leak must move under the same rebasis"
        ),
    }


def _basis_invariance_gate(check: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "passed": bool(check["passed"]),
        "check": check,
        "criteria": check["criteria"],
    }


def _outgoing_physicality_gate(packet: Mapping[str, Any]) -> dict[str, Any]:
    coeff = packet["coefficients"]
    finite = all(math.isfinite(float(coeff[key])) for key in COEFF_KEYS)
    n0_positive = float(coeff["N0"]) > 0.0
    rows = packet["rows"]
    cons_real = all(
        abs(_complex_parts(row["sigmaCons"])[1])
        <= 1.0e-8 * max(abs(_complex_parts(row["sigmaCons"])[0]), 1.0)
        for row in rows
    )
    radiating = all(float(row["minusImSigmaRet"]) > 0.0 and float(row["outgoingFlux"]) > 0.0 for row in rows)
    return {
        "passed": bool(finite and n0_positive and cons_real and radiating),
        "finite": finite,
        "N0_positive": n0_positive,
        "Sigma_cons_real": cons_real,
        "minus_Im_Sigma_ret_positive": radiating,
        "min_minus_im_sigma_ret": min(float(row["minusImSigmaRet"]) for row in rows),
        "min_outgoing_flux": min(float(row["outgoingFlux"]) for row in rows),
        "criteria": "finite coefficients, N0>0, conservative self-energy real, and positive outgoing radiated power",
    }


def _n_channel_robustness(packet: Mapping[str, Any]) -> dict[str, Any]:
    rows = packet["rows"]
    minus_im = np.asarray([float(row["minusImSigmaRet"]) for row in rows], dtype=np.float64)
    sigma_scale = np.asarray(
        [abs(complex(*_complex_parts(row["sigmaRet"]))) for row in rows],
        dtype=np.float64,
    )
    residuals = np.asarray(
        [max(float(row["consResidual"]), float(row["retResidual"])) for row in rows],
        dtype=np.float64,
    )
    residual_floor = float(np.max(residuals * np.maximum(sigma_scale, 1.0e-300)))
    conservative_imag_floor = float(
        max(abs(_complex_parts(row["sigmaCons"])[1]) for row in rows)
    )
    fit_rms_n = float(packet["diagnostics"].get("fit_rms_N", math.inf))
    nof_omega = np.asarray([float(row["NofOmega"]) for row in rows], dtype=np.float64)
    n_fit_floor = fit_rms_n
    omegas = np.asarray([float(row["omega"]) for row in rows], dtype=np.float64)
    signal_floor = max(residual_floor, conservative_imag_floor, 1.0e-300)
    min_signal = float(np.min(minus_im))
    signal_to_floor = float(min_signal / signal_floor)
    coeff = packet["coefficients"]
    n_term_contributions = {
        "N0": abs(float(coeff["N0"])),
        "N2": float(np.max(np.abs(float(coeff["N2"]) * omegas**2))),
        "N4": float(np.max(np.abs(float(coeff["N4"]) * omegas**4))),
    }
    n_term_to_fit = {
        key: float(value / max(n_fit_floor, 1.0e-300))
        for key, value in n_term_contributions.items()
    }
    n_terms_above_fit_floor = {
        key: bool(value >= 1.0e3)
        for key, value in n_term_to_fit.items()
    }
    n_signal_to_fit = n_term_to_fit["N0"]
    n_fit_condition_ok = bool(float(packet["diagnostics"].get("fit_condition_number_N", math.inf)) <= 1.0e3)
    return {
        "passed": bool(
            min_signal > 0.0
            and signal_to_floor >= 1.0e3
            and all(n_terms_above_fit_floor.values())
            and n_fit_condition_ok
        ),
        "min_minus_im_sigma_ret": min_signal,
        "max_linear_solve_imag_floor": signal_floor,
        "signal_to_solve_floor": signal_to_floor,
        "fit_condition_number_N": float(packet["diagnostics"].get("fit_condition_number_N", math.inf)),
        "fit_condition_ok": n_fit_condition_ok,
        "fit_rms_NofOmega": n_fit_floor,
        "N0_to_fit_rms": n_signal_to_fit,
        "N_term_contributions_over_window": n_term_contributions,
        "N_term_to_fit_rms": n_term_to_fit,
        "N_terms_above_fit_floor": n_terms_above_fit_floor,
        "NofOmega_range": {
            "min": float(np.min(nof_omega)),
            "max": float(np.max(nof_omega)),
            "span": float(np.max(nof_omega) - np.min(nof_omega)),
        },
        "criteria": (
            "-Im Sigma_ret positive and at least 1e3 above the conservative/solve residual imaginary floor; "
            "scaled N fit condition <= 1e3 and each N0/N2*omega^2/N4*omega^4 contribution "
            "at least 1e3 above fit RMS"
        ),
    }


def _self_consistency_gate(backgrounds: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    rows = []
    for background in backgrounds:
        rows.append(
            {
                "tau": float(background["constants"]["tau"]),
                "residual_linf": float(background["residuals"]["closed_stationary_linf"]),
                "self_consistent": bool(background["residuals"]["self_consistent"]),
                "a": float(background["geometry"]["a"]),
                "L": float(background["geometry"]["L"]),
                "background_hash": background["content_hash"],
            }
        )
    passed = all(
        row["self_consistent"]
        and row["residual_linf"] < 1.0e-6
        and row["residual_linf"] < SMOKE_RESIDUAL_REFERENCE / 1.0e6
        and abs(row["a"] - FROZEN_A) <= 1.0e-15
        and abs(row["L"] - FROZEN_L) <= 1.0e-15
        for row in rows
    )
    return {
        "passed": bool(passed),
        "rows": rows,
        "smoke_residual_reference": SMOKE_RESIDUAL_REFERENCE,
        "criteria": "closed residual << smoke 243.39 and frozen geometry a=1,L=1.85",
    }


def _relative_change(new: float, old: float) -> float:
    return abs(float(new) - float(old)) / max(abs(float(new)), 1.0e-300)


def _coeff_rel_change(a: Mapping[str, Any], b: Mapping[str, Any]) -> dict[str, float]:
    return {
        key: _relative_change(float(a["coefficients"][key]), float(b["coefficients"][key]))
        for key in COEFF_KEYS
    }


def _sweep_table(
    packets: Sequence[Mapping[str, Any]],
    *,
    label_key: str,
) -> list[dict[str, Any]]:
    table: list[dict[str, Any]] = []
    previous: Mapping[str, Any] | None = None
    for packet in packets:
        row = {
            label_key: packet["grid"].get(label_key, f"{packet['grid']['nr']}x{packet['grid']['nw']}"),
            "grid": f"{packet['grid']['nr']}x{packet['grid']['nw']}",
            "window": float(packet["grid"]["window"]),
            "radial_scale": float(packet["grid"]["radial_scale"]),
            **{key: float(packet["coefficients"][key]) for key in COEFF_KEYS},
        }
        if previous is None:
            row["rel_change_max"] = None
            row["rel_change"] = None
        else:
            rel = _coeff_rel_change(packet, previous)
            row["rel_change"] = rel
            row["rel_change_max"] = max(rel.values())
        table.append(row)
        previous = packet
    return table


def _sweep_convergence_summary(
    table: Sequence[Mapping[str, Any]],
    *,
    tolerance: float = CONVERGENCE_REL_TOL,
    shrink_factor: float = CONVERGENCE_SHRINK_FACTOR,
) -> dict[str, Any]:
    if len(table) < 3:
        return {
            "passed": False,
            "reason": "requires at least three levels",
            "step_max": [],
            "final_max": math.inf,
            "previous_max": math.inf,
            "per_coefficient": {},
        }
    per_coeff: dict[str, dict[str, Any]] = {}
    for key in COEFF_KEYS:
        increments = [float(row["rel_change"][key]) for row in table[1:]]
        finite = bool(all(math.isfinite(value) for value in increments))
        under_tol = bool(increments[-1] <= tolerance) if increments else False
        shrinking_steps = [
            bool(increments[i] <= shrink_factor * max(increments[i - 1], 1.0e-300))
            for i in range(1, len(increments))
        ]
        shrinking = bool(shrinking_steps and all(shrinking_steps))
        per_coeff[key] = {
            "increments": increments,
            "final_increment": increments[-1] if increments else math.inf,
            "previous_increment": increments[-2] if len(increments) >= 2 else math.inf,
            "finite": finite,
            "shrinking": shrinking,
            "under_tolerance": under_tol,
            "passed": bool(finite and shrinking and under_tol),
            "tolerance": tolerance,
            "shrink_factor": shrink_factor,
        }
    step_max = [float(row["rel_change_max"]) for row in table[1:]]
    final_max = step_max[-1]
    previous_max = step_max[-2]
    finite = bool(all(math.isfinite(value) for value in step_max))
    return {
        "passed": bool(finite and all(row["passed"] for row in per_coeff.values())),
        "finite": finite,
        "shrinking": bool(all(row["shrinking"] for row in per_coeff.values())),
        "under_tolerance": bool(all(row["under_tolerance"] for row in per_coeff.values())),
        "step_max": step_max,
        "final_max": final_max,
        "previous_max": previous_max,
        "per_coefficient": per_coeff,
        "tolerance": tolerance,
        "shrink_factor": shrink_factor,
        "criteria": (
            f">=3 levels, every coefficient's final relative increment <= {tolerance:.2%}, "
            f"and every coefficient's increments shrink by factor <= {shrink_factor:.2f}; "
            "no max-over-coefficients aggregation is used for pass/fail"
        ),
    }


def _convergence_gate(
    *,
    grid_packets: Sequence[Mapping[str, Any]],
    window_packets: Sequence[Mapping[str, Any]],
    truncation_packets: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    grid_table = _sweep_table(grid_packets, label_key="grid_label")
    for row in grid_table:
        row["grid_label"] = row["grid"]
    window_table = _sweep_table(window_packets, label_key="window")
    trunc_table = _sweep_table(truncation_packets, label_key="radial_scale")
    grid_summary = _sweep_convergence_summary(grid_table)
    window_summary = _sweep_convergence_summary(window_table)
    trunc_summary = _sweep_convergence_summary(trunc_table)
    passed = bool(grid_summary["passed"] and window_summary["passed"] and trunc_summary["passed"])
    error_rel = {
        key: max(
            float(grid_table[-1]["rel_change"][key]),
            float(window_table[-1]["rel_change"][key]),
            float(trunc_table[-1]["rel_change"][key]),
        )
        for key in COEFF_KEYS
    }
    return {
        "passed": passed,
        "grid_table": grid_table,
        "omega_window_table": window_table,
        "truncation_table": trunc_table,
        "grid_summary": grid_summary,
        "omega_window_summary": window_summary,
        "truncation_summary": trunc_summary,
        "error_bars_rel": error_rel,
        "error_bars_max_rel": max(error_rel.values()),
        "chosen_resolution": {
            "grid": grid_table[-1]["grid"],
            "window": window_table[-1]["window"],
            "radial_scale": trunc_table[-1]["radial_scale"],
        },
        "criteria": grid_summary["criteria"],
    }


def _grid_convergence_gate(grid_packets: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    grid_table = _sweep_table(grid_packets, label_key="grid_label")
    for row in grid_table:
        row["grid_label"] = row["grid"]
    summary = _sweep_convergence_summary(grid_table)
    error_rel = {key: float(grid_table[-1]["rel_change"][key]) for key in COEFF_KEYS}
    return {
        "passed": bool(summary["passed"]),
        "grid_table": grid_table,
        "grid_summary": summary,
        "error_bars_rel": error_rel,
        "error_bars_max_rel": max(error_rel.values()),
        "criteria": summary["criteria"],
    }


def _tau_sensitivity_gate(tau1_packet: Mapping[str, Any], tau2_packet: Mapping[str, Any]) -> dict[str, Any]:
    diffs = _max_diffs(_coeff_values(tau1_packet), _coeff_values(tau2_packet))
    per = {
        key: _max_diffs([float(tau1_packet["coefficients"][key])], [float(tau2_packet["coefficients"][key])])
        for key in COEFF_KEYS
    }
    return {
        "passed": bool(diffs["max_rel"] >= 1.0e-6),
        "tau1": float(tau1_packet["tau"]),
        "tau2": float(tau2_packet["tau"]),
        "max_abs": diffs["max_abs"],
        "max_rel": diffs["max_rel"],
        "per_coefficient": per,
        "criteria": "tau=1 vs tau=2 max relative movement in Z/N coefficients >= 1e-6",
    }


def _synthetic_packet(value: float, *, grid: tuple[int, int], window: float, radial_scale: float) -> dict[str, Any]:
    return {
        "grid": {"nr": grid[0], "nw": grid[1], "window": window, "radial_scale": radial_scale},
        "coefficients": {key: float(value) for key in COEFF_KEYS},
    }


def _synthetic_coeff_packet(
    coeff: Mapping[str, float],
    *,
    grid: tuple[int, int],
    window: float,
    radial_scale: float,
) -> dict[str, Any]:
    return {
        "grid": {"nr": grid[0], "nw": grid[1], "window": window, "radial_scale": radial_scale},
        "coefficients": {key: float(coeff[key]) for key in COEFF_KEYS},
    }


def convergence_negative_control() -> dict[str, Any]:
    base = {key: 1.0 for key in COEFF_KEYS}
    shrink_one_grow_another = [
        {**base, "Z0": 1.0, "N4": 1.0},
        {**base, "Z0": 1.20, "N4": 1.05},
        {**base, "Z0": 1.26, "N4": 1.16},
    ]
    grid_packets = [
        _synthetic_coeff_packet(coeff, grid=grid, window=DEFAULT_FINAL_WINDOW, radial_scale=DEFAULT_FINAL_TRUNCATION)
        for coeff, grid in zip(shrink_one_grow_another, ((31, 13), (39, 15), (47, 17)))
    ]
    window_packets = [
        _synthetic_coeff_packet(coeff, grid=DEFAULT_FINAL_GRID, window=window, radial_scale=DEFAULT_FINAL_TRUNCATION)
        for coeff, window in zip(shrink_one_grow_another, (0.058, 0.046, 0.036))
    ]
    trunc_packets = [
        _synthetic_coeff_packet(coeff, grid=DEFAULT_FINAL_GRID, window=DEFAULT_FINAL_WINDOW, radial_scale=scale)
        for coeff, scale in zip(shrink_one_grow_another, DEFAULT_TRUNCATION_SCALES)
    ]
    gate = _convergence_gate(
        grid_packets=grid_packets,
        window_packets=window_packets,
        truncation_packets=trunc_packets,
    )
    return {
        "failed_as_expected": not gate["passed"],
        "wrong_answer": "a grow-one/shrink-another sequence hidden by a max-over-coefficients convergence check",
        "gate_passed": gate["passed"],
        "grow_one_shrink_another": True,
        "final_grid_increment": gate["grid_summary"]["final_max"],
        "final_truncation_increment": gate["truncation_summary"]["final_max"],
        "N4_truncation_increments": gate["truncation_summary"]["per_coefficient"]["N4"]["increments"],
        "Z0_truncation_increments": gate["truncation_summary"]["per_coefficient"]["Z0"]["increments"],
    }


def self_consistency_negative_control(background: Mapping[str, Any]) -> dict[str, Any]:
    stale = copy.deepcopy(dict(background))
    stale["content_hash"] = "negative-control-stale-smoke-background"
    stale["residuals"] = {
        "closed_stationary_linf": SMOKE_RESIDUAL_REFERENCE,
        "self_consistent": False,
    }
    gate = _self_consistency_gate([stale])
    return {
        "failed_as_expected": not gate["passed"],
        "wrong_answer": "a stale smoke-level background or wrong-geometry fallback",
        "gate_passed": gate["passed"],
    }


def tau_sensitivity_negative_control(packet: Mapping[str, Any]) -> dict[str, Any]:
    stale = copy.deepcopy(dict(packet))
    stale["tau"] = 2.0
    if "constants" in stale:
        stale["constants"] = copy.deepcopy(stale["constants"])
        stale["constants"]["tau"] = 2.0
    gate = _tau_sensitivity_gate(packet, stale)
    return {
        "failed_as_expected": not gate["passed"],
        "wrong_answer": "a frozen tau=1 transfer relabelled as tau=2",
        "gate_passed": gate["passed"],
        "max_rel": gate["max_rel"],
    }


def dual_engine_negative_control(packet: Mapping[str, Any]) -> dict[str, Any]:
    perturbed = copy.deepcopy(dict(packet))
    perturbed["coefficients"] = copy.deepcopy(packet["coefficients"])
    perturbed["coefficients"]["N2"] = float(perturbed["coefficients"]["N2"]) * 1.25
    diff = compare_engine_packets(packet, perturbed)
    gate = _dual_engine_gate([{**diff, "tau": packet.get("tau", 1.0), "grid": "negative"}])
    return {
        "failed_as_expected": not gate["passed"],
        "wrong_answer": "a second discretization that lands on a materially different coefficient",
        "gate_passed": gate["passed"],
        "max_rel": gate["max_rel"],
    }


def basis_invariance_negative_control(background: Mapping[str, Any], bdg_packet: Mapping[str, Any]) -> dict[str, Any]:
    leaked = basis_invariance_check(background, bdg_packet, inject_port_leak=True)
    return {
        "failed_as_expected": not leaked["passed"],
        "wrong_answer": leaked["port_leak_negative_control"]["wrong_answer"],
        "gate_passed": leaked["passed"],
        "fixed_port_relative_delta": leaked["max_wrong_fixed_port_relative_delta"],
        "tested_extraction_relative_delta": leaked["max_tested_extraction_relative_delta"],
        "branch_difference": leaked["branch_difference"],
    }


def outgoing_physicality_negative_control(packet: Mapping[str, Any]) -> dict[str, Any]:
    bad = copy.deepcopy(dict(packet))
    bad["coefficients"] = copy.deepcopy(packet["coefficients"])
    bad["coefficients"]["N0"] = -abs(float(bad["coefficients"]["N0"]))
    bad["rows"] = copy.deepcopy(packet["rows"])
    for row in bad["rows"]:
        row["minusImSigmaRet"] = -abs(float(row["minusImSigmaRet"]))
    gate = _outgoing_physicality_gate(bad)
    return {
        "failed_as_expected": not gate["passed"],
        "wrong_answer": "an anti-radiating or non-positive-N0 transfer",
        "gate_passed": gate["passed"],
    }


def consumer_negative_control(packet: Mapping[str, Any], bdg_packet: Mapping[str, Any]) -> dict[str, Any]:
    bad = copy.deepcopy(dict(packet))
    bad["coefficients"] = copy.deepcopy(packet["coefficients"])
    bad["coefficients"]["N4"] = float(bad["coefficients"]["N4"]) * 1.25
    gate = _consumer_compatibility_gate(
        primary_packet=packet,
        independent_packet=bad,
        bdg_packet=bdg_packet,
    )
    return {
        "failed_as_expected": not gate["passed"],
        "wrong_answer": "a B1 direct-coefficient lane with a cross-engine mismatch",
        "gate_passed": gate["passed"],
        "max_rel": gate["cross_engine_maxwell_diffs"]["max_rel"],
    }


def _direct_coefficients(
    *,
    maxwell_packet: Mapping[str, Any],
    bdg_packet: Mapping[str, Any],
) -> dict[str, float]:
    wall = bdg_packet["wall"]
    bdg = bdg_packet["bdg_moments"]
    coeff = {
        "K": float(wall["K"]),
        "M": float(wall["M"]),
        "B0": float(bdg["B0"]),
        "B2": float(bdg["B2"]),
        "B4": float(bdg["B4"]),
    }
    coeff.update({key: float(maxwell_packet["coefficients"][key]) for key in COEFF_KEYS})
    return coeff


def _consumer_compatibility_gate(
    *,
    primary_packet: Mapping[str, Any],
    independent_packet: Mapping[str, Any],
    bdg_packet: Mapping[str, Any],
) -> dict[str, Any]:
    primary_coeff = _direct_coefficients(maxwell_packet=primary_packet, bdg_packet=bdg_packet)
    independent_coeff = _direct_coefficients(maxwell_packet=independent_packet, bdg_packet=bdg_packet)
    missing = sorted(set(DIRECT_KEYS) - set(primary_coeff))
    finite = all(math.isfinite(float(primary_coeff[key])) for key in DIRECT_KEYS)
    n0_positive = float(primary_coeff["N0"]) > 0.0
    diffs = _max_diffs(
        [primary_coeff[key] for key in COEFF_KEYS],
        [independent_coeff[key] for key in COEFF_KEYS],
    )
    cross_pass = _strict_diff_pass(
        diffs,
        abs_tol=DUAL_ENGINE_TOLERANCES["abs"],
        rel_tol=DUAL_ENGINE_TOLERANCES["rel"],
    )
    return {
        "passed": bool(not missing and finite and n0_positive and cross_pass),
        "direct_coefficients": primary_coeff,
        "independent_direct_coefficients": independent_coeff,
        "missing": missing,
        "finite": finite,
        "N0_positive": n0_positive,
        "cross_engine_maxwell_diffs": diffs,
        "lane_extract_scope_note": (
            "Validated the direct_coefficients input contract only. The actual lane_extract call is not "
            "made here because that function immediately assembles D0/P0/R_norm-scope derived quantities."
        ),
        "criteria": "required direct keys present, finite, N0>0, and primary/independent Z/N match through the coefficient dict",
    }


def _load_packets_for_validation(
    *,
    run_root: Path,
    background: Mapping[str, Any],
    tau: float,
    grids: Sequence[tuple[int, int]],
    windows: Sequence[float],
    truncation_scales: Sequence[float],
    final_grid: tuple[int, int],
    final_window: float,
    final_truncation: float,
) -> dict[str, Any]:
    grid_packets = [
        _load_json(_default_packet_path(run_root, "python", tau, nr, nw, final_window, final_truncation))
        for nr, nw in grids
    ]
    window_packets = [
        _load_json(_default_packet_path(run_root, "python", tau, *final_grid, window, final_truncation))
        for window in windows
    ]
    truncation_packets = [
        _load_json(
            _default_packet_path(
                run_root,
                "python",
                tau,
                *_radial_grid_for_scale(
                    background,
                    final_grid=final_grid,
                    final_truncation=final_truncation,
                    radial_scale=scale,
                ),
                final_window,
                scale,
            )
        )
        for scale in truncation_scales
    ]
    independent_grid_packets = [
        _load_json(_default_packet_path(run_root, "independent", tau, nr, nw, final_window, final_truncation))
        for nr, nw in grids
    ]
    independent_window_packets = [
        _load_json(_default_packet_path(run_root, "independent", tau, *final_grid, window, final_truncation))
        for window in windows
    ]
    independent_truncation_packets = [
        _load_json(
            _default_packet_path(
                run_root,
                "independent",
                tau,
                *_radial_grid_for_scale(
                    background,
                    final_grid=final_grid,
                    final_truncation=final_truncation,
                    radial_scale=scale,
                ),
                final_window,
                scale,
            )
        )
        for scale in truncation_scales
    ]
    final_primary = _load_json(_default_packet_path(run_root, "python", tau, *final_grid, final_window, final_truncation))
    final_independent = _load_json(_default_packet_path(run_root, "independent", tau, *final_grid, final_window, final_truncation))
    return {
        "grid_packets": grid_packets,
        "window_packets": window_packets,
        "truncation_packets": truncation_packets,
        "independent_grid_packets": independent_grid_packets,
        "independent_window_packets": independent_window_packets,
        "independent_truncation_packets": independent_truncation_packets,
        "final_primary": final_primary,
        "final_independent": final_independent,
    }


def validate_and_report(
    *,
    run_root: Path,
    report_path: Path,
    b2a_run_root: Path,
    tau1: float,
    tau2: float,
    grids: Sequence[tuple[int, int]],
    windows: Sequence[float],
    truncation_scales: Sequence[float],
    final_grid: tuple[int, int],
    final_window: float,
    final_truncation: float,
) -> dict[str, Any]:
    background1 = _load_json(_default_background_path(b2a_run_root, tau1))
    background2 = _load_json(_default_background_path(b2a_run_root, tau2))
    bdg1 = _load_json(_default_b2a_bdg_path(b2a_run_root, tau1))
    bdg2 = _load_json(_default_b2a_bdg_path(b2a_run_root, tau2))
    packets = _load_packets_for_validation(
        run_root=run_root,
        background=background1,
        tau=tau1,
        grids=grids,
        windows=windows,
        truncation_scales=truncation_scales,
        final_grid=final_grid,
        final_window=final_window,
        final_truncation=final_truncation,
    )
    final_primary = packets["final_primary"]
    final_independent = packets["final_independent"]
    tau2_primary = _load_json(_default_packet_path(run_root, "python", tau2, *final_grid, final_window, final_truncation))
    tau2_independent = _load_json(_default_packet_path(run_root, "independent", tau2, *final_grid, final_window, final_truncation))

    dual_rows = []
    dual_specs = [(tau1, *final_grid, final_window, final_truncation), (tau2, *final_grid, final_window, final_truncation)]
    for tau, nr, nw, window, scale in dual_specs:
        primary = _load_json(_default_packet_path(run_root, "python", tau, nr, nw, window, scale))
        independent = _load_json(_default_packet_path(run_root, "independent", tau, nr, nw, window, scale))
        dual_rows.append(
            {
                "tau": float(tau),
                "grid": f"{nr}x{nw}",
                "window": float(window),
                "radial_scale": float(scale),
                **compare_engine_packets(primary, independent),
            }
        )
    dual_gate = _dual_engine_gate(dual_rows)
    self_gate = _self_consistency_gate([background1, background2])
    basis_check = basis_invariance_check(background1, bdg1)
    basis_gate = _basis_invariance_gate(basis_check)
    physical_gate = _outgoing_physicality_gate(final_primary)
    n_channel_gate = _n_channel_robustness(final_primary)
    convergence_gate = _convergence_gate(
        grid_packets=packets["grid_packets"],
        window_packets=packets["window_packets"],
        truncation_packets=packets["truncation_packets"],
    )
    independent_convergence_gate = _convergence_gate(
        grid_packets=packets["independent_grid_packets"],
        window_packets=packets["independent_window_packets"],
        truncation_packets=packets["independent_truncation_packets"],
    )
    tau_gate = _tau_sensitivity_gate(final_primary, tau2_primary)
    consumer_gate = _consumer_compatibility_gate(
        primary_packet=final_primary,
        independent_packet=final_independent,
        bdg_packet=bdg1,
    )
    independent_tau_gate = _tau_sensitivity_gate(final_independent, tau2_independent)
    negative_controls = {
        "self_consistency": self_consistency_negative_control(background1),
        "dual_engine_agreement": dual_engine_negative_control(final_primary),
        "basis_invariance": basis_invariance_negative_control(background1, bdg1),
        "outgoing_physicality": outgoing_physicality_negative_control(final_primary),
        "convergence": convergence_negative_control(),
        "tau_sensitivity": tau_sensitivity_negative_control(final_primary),
        "b1_consumer_compatibility": consumer_negative_control(final_primary, bdg1),
    }
    no_cant_fail_gate = {
        "passed": bool(all(control.get("failed_as_expected", False) for control in negative_controls.values())),
        "negative_controls": negative_controls,
        "criteria": "each validation gate has a concrete negative control that fails as expected",
    }
    legacy_7x7_path = _default_packet_path(run_root, "python", tau1, 7, 7, 0.036, 1.0)
    legacy_7x7 = _load_json(legacy_7x7_path) if legacy_7x7_path.exists() else packets["grid_packets"][0]
    moved_vs_legacy = {
        key: _relative_change(float(final_primary["coefficients"][key]), float(legacy_7x7["coefficients"][key]))
        for key in COEFF_KEYS
    }
    error_budget = {
        "Maxwell_ZN_rel": convergence_gate["error_bars_rel"],
        "Maxwell_ZN_max_rel": convergence_gate["error_bars_max_rel"],
        "independent_ZN_rel": independent_convergence_gate["error_bars_rel"],
        "dual_engine_max_rel": dual_gate["max_rel"],
        "components": {
            "grid_finest_vs_previous": convergence_gate["grid_table"][-1]["rel_change"],
            "omega_window_finest_vs_previous": convergence_gate["omega_window_table"][-1]["rel_change"],
            "radial_truncation_finest_vs_previous": convergence_gate["truncation_table"][-1]["rel_change"],
            "independent_grid_finest_vs_previous": independent_convergence_gate["grid_table"][-1]["rel_change"],
            "independent_omega_window_finest_vs_previous": independent_convergence_gate["omega_window_table"][-1]["rel_change"],
            "independent_radial_truncation_finest_vs_previous": independent_convergence_gate["truncation_table"][-1]["rel_change"],
            "dual_engine_converged_value": {
                key: dual_gate["per_coefficient"][key]["max_rel"] for key in COEFF_KEYS
            },
        },
    }
    direct_coefficients = consumer_gate["direct_coefficients"]
    final_bundle: dict[str, Any] = {
        "schema": "stage1_patha_b2b_validated_maxwell_transfer/v1",
        "source": "primary Python Path-A B2b output after genuinely independent staggered/high-order validation",
        "scope": "Maxwell transfer only; target-blind; no D0/P0/R_norm/root-find assembly",
        "background": {
            "path": str(_default_background_path(b2a_run_root, tau1)),
            "content_hash": background1["content_hash"],
            "residual_linf": float(background1["residuals"]["closed_stationary_linf"]),
            "smoke_residual_reference": SMOKE_RESIDUAL_REFERENCE,
        },
        "bdg_bundle": {
            "path": str(_default_b2a_bdg_path(b2a_run_root, tau1)),
            "content_hash": bdg1["content_hash"],
        },
        "tau": float(tau1),
        "geometry": final_primary["geometry"],
        "constants": final_primary["constants"],
        "grid": final_primary["grid"],
        "coefficients": final_primary["coefficients"],
        "direct_coefficients": direct_coefficients,
        "diagnostics": final_primary["diagnostics"],
        "error_budget": error_budget,
        "validation": {
            "passed": False,
            "gates": [],
            "dual_engine_rows": dual_rows,
            "dual_engine_gate": dual_gate,
            "self_consistency": self_gate,
            "basis_invariance": basis_gate,
            "outgoing_physicality": physical_gate,
            "N_channel_robustness": n_channel_gate,
            "convergence": convergence_gate,
            "independent_convergence": independent_convergence_gate,
            "tau_sensitivity": tau_gate,
            "independent_tau_sensitivity": independent_tau_gate,
            "b1_consumer_compatibility": consumer_gate,
            "no_cant_fail_gates": no_cant_fail_gate,
            "legacy_unconverged_7x7": {
                "path": str(legacy_7x7_path) if legacy_7x7_path.exists() else "grid ladder first level",
                "coefficients": legacy_7x7["coefficients"],
                "relative_move_to_converged": moved_vs_legacy,
            },
        },
    }
    gates = [
        {
            "gate": "self_consistency",
            "passed": self_gate["passed"],
            "catches": "running the old smoke R0/Z(w) background, stale tau bundle, or wrong M1a geometry",
        },
        {
            "gate": "dual_engine_agreement",
            "passed": dual_gate["passed"],
            "catches": "a discretization-specific VSH assembly, Green solve, source normalization, or abs-only tiny-number false pass",
        },
        {
            "gate": "basis_invariance",
            "passed": basis_gate["passed"],
            "catches": "basis-dependent U/W-port leakage or rotating the operator without the source/boundary quadratic",
        },
        {
            "gate": "outgoing_physicality",
            "passed": physical_gate["passed"],
            "catches": "wrong DtN sign, bound non-radiating solution, anti-radiating solution, or N0<=0",
        },
        {
            "gate": "N_channel_robustness",
            "passed": n_channel_gate["passed"],
            "catches": "shipping an N-channel below the solve/fit noise floor or from an ill-conditioned omega fit",
        },
        {
            "gate": "convergence_truncation",
            "passed": convergence_gate["passed"] and independent_convergence_gate["passed"],
            "catches": "under-resolved mesh, unstable omega-window fit, or radial truncation dependence",
        },
        {
            "gate": "tau_sensitivity",
            "passed": tau_gate["passed"] and independent_tau_gate["passed"],
            "catches": "reusing a frozen tau=1 background/BdG packet instead of re-solving at tau=2",
        },
        {
            "gate": "b1_consumer_compatibility",
            "passed": consumer_gate["passed"],
            "catches": "missing direct-coefficient keys, non-finite scalar, N0<=0, or cross-engine mismatch through B1 input shape",
        },
        {
            "gate": "no_cant_fail_gates",
            "passed": no_cant_fail_gate["passed"],
            "catches": "report-only validations without a concrete negative-control failure mode",
        },
    ]
    passed = bool(all(gate["passed"] for gate in gates))
    final_bundle["validation"]["passed"] = passed
    final_bundle["validation"]["gates"] = gates
    final_bundle["content_hash"] = _full_stable_hash(final_bundle)
    bundle_dir = run_root / "bundles"
    bundle_path = bundle_dir / f"patha_b2b_validated_maxwell_transfer_tau_{_format_tau(tau1)}.json"
    _write_json(bundle_path, final_bundle)
    write_report(
        report_path=report_path,
        bundle_path=bundle_path,
        final_bundle=final_bundle,
        background1=background1,
        background2=background2,
        bdg1=bdg1,
        bdg2=bdg2,
        gates=gates,
    )
    return {
        "passed": passed,
        "report_path": str(report_path),
        "bundle_path": str(bundle_path),
        "gates": gates,
        "final_bundle_hash": final_bundle["content_hash"],
    }


def generate_validation_packets(
    *,
    run_root: Path,
    b2a_run_root: Path,
    tau1: float,
    tau2: float,
    grids: Sequence[tuple[int, int]],
    windows: Sequence[float],
    truncation_scales: Sequence[float],
    final_grid: tuple[int, int],
    final_window: float,
    final_truncation: float,
    mode_count: int = DEFAULT_MODE_COUNT,
) -> dict[str, Any]:
    specs: dict[tuple[str, float, int, int, float, float], tuple[str, str]] = {}
    background1 = _load_json(_default_background_path(b2a_run_root, tau1))

    def add(engine: str, tau: float, grid: tuple[int, int], window: float, scale: float, discretization: str) -> None:
        specs[(engine, float(tau), int(grid[0]), int(grid[1]), float(window), float(scale))] = (
            engine,
            discretization,
        )

    for grid in grids:
        add("python_primary_second_order", tau1, grid, final_window, final_truncation, "primary_second_order")
        add("python_staggered_high_order", tau1, grid, final_window, final_truncation, "staggered_high_order")
    for window in windows:
        add("python_primary_second_order", tau1, final_grid, window, final_truncation, "primary_second_order")
        add("python_staggered_high_order", tau1, final_grid, window, final_truncation, "staggered_high_order")
    for scale in truncation_scales:
        radial_grid = _radial_grid_for_scale(
            background1,
            final_grid=final_grid,
            final_truncation=final_truncation,
            radial_scale=scale,
        )
        add("python_primary_second_order", tau1, radial_grid, final_window, scale, "primary_second_order")
        add("python_staggered_high_order", tau1, radial_grid, final_window, scale, "staggered_high_order")
    add("python_primary_second_order", tau2, final_grid, final_window, final_truncation, "primary_second_order")
    add("python_staggered_high_order", tau2, final_grid, final_window, final_truncation, "staggered_high_order")

    written = []
    for (engine, tau, nr, nw, window, scale), (_, discretization) in sorted(specs.items()):
        background = _default_background_path(b2a_run_root, tau)
        bdg = _default_b2a_bdg_path(b2a_run_root, tau)
        path, packet = solve_python_command(
            background_path=background,
            bdg_path=bdg,
            nr=nr,
            nw=nw,
            window=window,
            radial_scale=scale,
            mode_count=mode_count,
            run_root=run_root,
            engine=engine,
            discretization=discretization,
        )
        written.append(
            {
                "path": str(path),
                "engine": engine,
                "discretization": discretization,
                "tau": tau,
                "grid": f"{nr}x{nw}",
                "window": window,
                "radial_scale": scale,
                "content_hash": packet["content_hash"],
            }
        )
    return {"written": written, "count": len(written)}


def _fmt(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.12e}"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, Mapping):
        return ", ".join(f"{k}={_fmt(v)}" for k, v in value.items())
    return str(value)


def _markdown_table(headers: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(_fmt(row.get(header)) for header in headers) + " |")
    return "\n".join(lines)


def _increment_rows(table: Sequence[Mapping[str, Any]], *, label_key: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    previous_label: Any | None = None
    for row in table:
        label = row.get(label_key, row.get("grid"))
        if previous_label is not None and row.get("rel_change") is not None:
            rows.append(
                {
                    "from": previous_label,
                    "to": label,
                    **{key: float(row["rel_change"][key]) for key in COEFF_KEYS},
                }
            )
        previous_label = label
    return rows


def write_report(
    *,
    report_path: Path,
    bundle_path: Path,
    final_bundle: Mapping[str, Any],
    background1: Mapping[str, Any],
    background2: Mapping[str, Any],
    bdg1: Mapping[str, Any],
    bdg2: Mapping[str, Any],
    gates: Sequence[Mapping[str, Any]],
) -> None:
    validation = final_bundle["validation"]
    coeff = final_bundle["coefficients"]
    direct = final_bundle["direct_coefficients"]
    dual_gate = validation["dual_engine_gate"]
    convergence = validation["convergence"]
    independent_convergence = validation["independent_convergence"]
    tau_gate = validation["tau_sensitivity"]
    independent_tau_gate = validation["independent_tau_sensitivity"]
    basis = validation["basis_invariance"]["check"]
    physical = validation["outgoing_physicality"]
    n_channel = validation["N_channel_robustness"]
    consumer = validation["b1_consumer_compatibility"]
    negative_controls = validation["no_cant_fail_gates"]["negative_controls"]
    legacy = validation["legacy_unconverged_7x7"]
    passed = bool(final_bundle["validation"]["passed"])
    a00_linf = float(np.max(np.abs(np.asarray(background1["fields"]["A_00"], dtype=np.float64))))
    ar0_linf = float(np.max(np.abs(np.asarray(background1["fields"]["A_r0"], dtype=np.float64))))
    aw0_linf = float(np.max(np.abs(np.asarray(background1["fields"]["A_w0"], dtype=np.float64))))
    lines: list[str] = []
    lines.append("# Path-A B2b Maxwell Transfer Derivation")
    lines.append("")
    lines.append(f"Overall B2b gate: {'PASS' if passed else 'FAIL'}")
    lines.append(f"Validated bundle: `{bundle_path}`")
    lines.append(f"Bundle content hash: `{final_bundle['content_hash']}`")
    lines.append("")
    lines.append("## Scope")
    lines.append("")
    lines.append("This run derives only `{Z0,Z2,Z4,N0,N2,N4}` on the Path-A closed background and emits the B1 `direct_coefficients` shape. No `D0`, `P0`, `R_norm`, `R_pole`, `P2`, `P4`, root-find, or `mt15_05` preview was run.")
    lines.append("")
    lines.append("## Canonical Sources")
    lines.append("")
    lines.append("- Five-lane VSH Maxwell operator and open DtN design: `software/stage1_solver/decisions/05_mixed_maxwell_spike_design.md` D1-D4.")
    lines.append("- Localized Maxwell PDE and mixed invariants: `notes/moving_throat_pde_program_compact.md:592`, `:674`, `:769`.")
    lines.append("- Compact outgoing normalization: `Gamma_port=4a^5/(27 c_s^5)` from `notes/moving_throat_pde_program_compact.md:2559`.")
    lines.append("- The old U/W formulas were used only as a regression reference, not as the live extraction path: `research/pde_audit/notes/stage_v2_09_maxwell_mixed_kernel_derivation.md:117` and `:299`.")
    lines.append("")
    lines.append("## Background")
    lines.append("")
    lines.append(
        f"`tau=1` closed residual Linf `{float(background1['residuals']['closed_stationary_linf']):.12e}` "
        f"versus old smoke residual `{SMOKE_RESIDUAL_REFERENCE:.12e}`."
    )
    lines.append(
        f"`tau=2` closed residual Linf `{float(background2['residuals']['closed_stationary_linf']):.12e}`."
    )
    lines.append(
        f"Frozen geometry: `a={background1['geometry']['a']}`, `L={background1['geometry']['L']}`. "
        f"Background hash: `{background1['content_hash']}`."
    )
    lines.append(
        f"Rest-background gauge-field linf: `A_00={a00_linf:.12e}`, `A_r0={ar0_linf:.12e}`, "
        f"`A_w0={aw0_linf:.12e}`. The B2a rest defect has `A_r0 == A_w0 == 0`, as expected for "
        "no rest current; this calibration transfer therefore runs through scalar `A_00`, kinematic `j_E`, "
        "and `q*delta_rho` channels. Vector-current channels are held-out non-rest/excited-defect coverage."
    )
    lines.append("")
    lines.append("## Derived Tau-1 Transfer")
    lines.append("")
    lines.append(
        f"`Z0={float(coeff['Z0']):.12e}`, `Z2={float(coeff['Z2']):.12e}`, `Z4={float(coeff['Z4']):.12e}`."
    )
    lines.append(
        f"`N0={float(coeff['N0']):.12e}`, `N2={float(coeff['N2']):.12e}`, `N4={float(coeff['N4']):.12e}`."
    )
    lines.append(
        f"`Gamma_port={float(final_bundle['constants']['Gamma_port']):.12e}` using `4a^5/(27c_s^5)` at `a=1`, `c_s=1`."
    )
    lines.append("Relative movement from the old unconverged 7x7/window-0.036/radial-scale-1.0 packet:")
    lines.append("")
    lines.append(
        _markdown_table(
            ["Z0", "Z2", "Z4", "N0", "N2", "N4"],
            [legacy["relative_move_to_converged"]],
        )
    )
    lines.append("")
    lines.append("Direct coefficients emitted for B1 shape:")
    lines.append("")
    lines.append(
        _markdown_table(
            ["K", "M", "B0", "B2", "B4", "Z0", "Z2", "Z4", "N0", "N2", "N4"],
            [direct],
        )
    )
    lines.append("")
    lines.append("## Dual Engine")
    lines.append("")
    lines.append(
        "Correctness dual-engine comparison is primary Python second-order interior-grid FD versus a genuinely "
        "different Python staggered cell-centered high-order discretization. The second engine uses different point "
        "locations, cell quadrature, and fourth-order interior derivative rows with a second-order boundary closure; "
        "it does not read or reuse primary packet outputs. The old Python-Wolfram files remain useful only as a "
        "transcription check of the primary discretization."
    )
    lines.append("")
    lines.append(
        _markdown_table(
            ["tau", "grid", "window", "radial_scale", "abs", "rel"],
            [
                {
                    "tau": row["tau"],
                    "grid": row["grid"],
                    "window": row["window"],
                    "radial_scale": row["radial_scale"],
                    "abs": row["coefficients"]["max_abs"],
                    "rel": row["coefficients"]["max_rel"],
                }
                for row in validation["dual_engine_rows"]
            ],
        )
    )
    lines.append(
        f"Max dual-engine coefficient abs/rel diff: `{dual_gate['max_abs']:.12e}`/`{dual_gate['max_rel']:.12e}`. "
        f"Criterion: {dual_gate['criteria']}."
    )
    lines.append("Per-coefficient dual diffs:")
    lines.append("")
    lines.append(
        _markdown_table(
            ["coefficient", "max_abs", "max_rel", "passed"],
            [
                {
                    "coefficient": key,
                    "max_abs": row["max_abs"],
                    "max_rel": row["max_rel"],
                    "passed": row["passed"],
                }
                for key, row in dual_gate["per_coefficient"].items()
            ],
        )
    )
    lines.append("")
    lines.append("Independent staggered-high-order engine grid convergence:")
    lines.append("")
    lines.append(
        _markdown_table(
            ["grid", "window", "radial_scale", "Z0", "Z2", "Z4", "N0", "N2", "N4", "rel_change_max"],
            independent_convergence["grid_table"],
        )
    )
    lines.append(
        "Independent omega-window convergence:"
    )
    lines.append("")
    lines.append(
        _markdown_table(
            ["grid", "window", "radial_scale", "Z0", "Z2", "Z4", "N0", "N2", "N4", "rel_change_max"],
            independent_convergence["omega_window_table"],
        )
    )
    lines.append("")
    lines.append("Independent outward-radial convergence:")
    lines.append("")
    lines.append(
        _markdown_table(
            ["grid", "window", "radial_scale", "Z0", "Z2", "Z4", "N0", "N2", "N4", "rel_change_max"],
            independent_convergence["truncation_table"],
        )
    )
    lines.append(
        f"Independent full-axis max relative error bar `{independent_convergence['error_bars_max_rel']:.12e}` "
        f"under criterion: {independent_convergence['criteria']}."
    )
    lines.append("")
    lines.append("## Convergence, Window, Truncation")
    lines.append("")
    lines.append(f"Chosen resolution: `{convergence['chosen_resolution']}`.")
    lines.append(f"Criterion: {convergence['criteria']}.")
    lines.append("")
    lines.append("Mesh sweep:")
    lines.append("")
    lines.append(
        _markdown_table(
            ["grid", "window", "radial_scale", "Z0", "Z2", "Z4", "N0", "N2", "N4", "rel_change_max"],
            convergence["grid_table"],
        )
    )
    lines.append("")
    lines.append("Omega-window sweep:")
    lines.append("")
    lines.append(
        _markdown_table(
            ["grid", "window", "radial_scale", "Z0", "Z2", "Z4", "N0", "N2", "N4", "rel_change_max"],
            convergence["omega_window_table"],
        )
    )
    lines.append("")
    lines.append("Radial truncation sweep:")
    lines.append("")
    lines.append(
        _markdown_table(
            ["grid", "window", "radial_scale", "Z0", "Z2", "Z4", "N0", "N2", "N4", "rel_change_max"],
            convergence["truncation_table"],
        )
    )
    lines.append("")
    lines.append("Primary per-coefficient outward-radial increments:")
    lines.append("")
    lines.append(
        _markdown_table(
            ["from", "to", "Z0", "Z2", "Z4", "N0", "N2", "N4"],
            _increment_rows(convergence["truncation_table"], label_key="radial_scale"),
        )
    )
    lines.append("")
    lines.append("Independent per-coefficient outward-radial increments:")
    lines.append("")
    lines.append(
        _markdown_table(
            ["from", "to", "Z0", "Z2", "Z4", "N0", "N2", "N4"],
            _increment_rows(independent_convergence["truncation_table"], label_key="radial_scale"),
        )
    )
    lines.append("")
    lines.append("Recorded per-coefficient relative error bars:")
    lines.append("")
    lines.append(
        _markdown_table(
            ["Z0", "Z2", "Z4", "N0", "N2", "N4"],
            [convergence["error_bars_rel"]],
        )
    )
    lines.append("")
    lines.append(
        f"Primary max relative error bar `{convergence['error_bars_max_rel']:.12e}`; "
        f"independent full-axis max relative error bar `{independent_convergence['error_bars_max_rel']:.12e}`. "
        f"Cross-engine converged-value max relative difference is `{dual_gate['max_rel']:.12e}`; interpret that "
        "against these convergence bars, not as precision beyond the demonstrated discretization spread."
    )
    lines.append("")
    lines.append("## Tau Sensitivity")
    lines.append("")
    lines.append(
        f"`tau=1` vs `tau=2` primary transfer max abs/rel movement "
        f"`{tau_gate['max_abs']:.12e}`/`{tau_gate['max_rel']:.12e}` ({tau_gate['criteria']})."
    )
    lines.append(
        f"`tau=1` vs `tau=2` independent transfer max abs/rel movement "
        f"`{independent_tau_gate['max_abs']:.12e}`/`{independent_tau_gate['max_rel']:.12e}`."
    )
    lines.append(
        f"B2a input hashes moved from `{bdg1['content_hash']}` to `{bdg2['content_hash']}`; "
        f"background hashes moved from `{background1['content_hash']}` to `{background2['content_hash']}`."
    )
    lines.append("")
    lines.append("## Basis Invariance And Physicality")
    lines.append("")
    lines.append(
        f"Basis rotation check: max relative `Z` delta `{basis['max_relative_Z_delta']:.12e}`, "
        f"max relative `N` delta `{basis['max_relative_N_delta']:.12e}`."
    )
    lines.append(
        f"Port-leak negative control: fixed-port relative movement "
        f"`{basis['port_leak_negative_control']['relative_delta']:.12e}`; "
        f"failed as expected: `{basis['port_leak_negative_control']['failed_as_expected']}`. "
        f"Toggling to the fixed-port branch changes the tested base value by relative "
        f"`{negative_controls['basis_invariance']['branch_difference']['relative_to_invariant_Z']:.12e}` "
        f"and gives tested extraction movement "
        f"`{negative_controls['basis_invariance']['tested_extraction_relative_delta']:.12e}`."
    )
    lines.append(
        f"Outgoing physicality: `N0_positive={physical['N0_positive']}`, "
        f"`Sigma_cons_real={physical['Sigma_cons_real']}`, "
        f"`minus_Im_Sigma_ret_positive={physical['minus_Im_Sigma_ret_positive']}`, "
        f"min `-Im Sigma_ret={physical['min_minus_im_sigma_ret']:.12e}`, "
        f"min flux `{physical['min_outgoing_flux']:.12e}`."
    )
    lines.append("")
    lines.append("## N-Channel Robustness")
    lines.append("")
    lines.append(
        f"Minimum `-Im Sigma_ret={n_channel['min_minus_im_sigma_ret']:.12e}`; solve/imaginary floor "
        f"`{n_channel['max_linear_solve_imag_floor']:.12e}`; signal/floor "
        f"`{n_channel['signal_to_solve_floor']:.12e}`."
    )
    lines.append(
        f"Scaled N fit condition `{n_channel['fit_condition_number_N']:.12e}` with fit RMS "
        f"`{n_channel['fit_rms_NofOmega']:.12e}` and `N0/fit_RMS={n_channel['N0_to_fit_rms']:.12e}`."
    )
    lines.append("Individual N-term fit-floor checks over the fitted omega window:")
    lines.append("")
    lines.append(
        _markdown_table(
            ["coefficient", "window_contribution", "to_fit_rms", "above_floor"],
            [
                {
                    "coefficient": key,
                    "window_contribution": n_channel["N_term_contributions_over_window"][key],
                    "to_fit_rms": n_channel["N_term_to_fit_rms"][key],
                    "above_floor": n_channel["N_terms_above_fit_floor"][key],
                }
                for key in N_KEYS
            ],
        )
    )
    lines.append("")
    lines.append("## B1 Consumer Compatibility")
    lines.append("")
    lines.append(
        f"Required direct keys present: `{not consumer['missing']}`; finite: `{consumer['finite']}`; "
        f"`N0>0`: `{consumer['N0_positive']}`. Cross-engine Z/N abs/rel through the coefficient dict: "
        f"`{consumer['cross_engine_maxwell_diffs']['max_abs']:.12e}`/"
        f"`{consumer['cross_engine_maxwell_diffs']['max_rel']:.12e}`."
    )
    lines.append(consumer["lane_extract_scope_note"])
    lines.append("")
    lines.append("## Gates And Wrong Answers Caught")
    lines.append("")
    lines.append(
        _markdown_table(
            ["gate", "status", "catches"],
            [
                {
                    "gate": gate["gate"],
                    "status": "PASS" if gate["passed"] else "FAIL",
                    "catches": gate["catches"],
                }
                for gate in gates
            ],
        )
    )
    lines.append("")
    lines.append("Negative controls proving gates can fail:")
    lines.append("")
    lines.append(
        _markdown_table(
            ["gate", "failed_as_expected", "wrong_answer"],
            [
                {
                    "gate": key,
                    "failed_as_expected": value.get("failed_as_expected"),
                    "wrong_answer": value.get("wrong_answer"),
                }
                for key, value in negative_controls.items()
            ],
        )
    )
    lines.append("")
    lines.append("## Dual-Engine Feasibility Note")
    lines.append("")
    lines.append(
        "A genuine independent l=2 transfer discretization was built: the primary engine uses the original "
        "interior-node second-order FD grid; the independent engine uses staggered cell centers, cell-volume "
        "quadrature, and fourth-order interior derivative stencils with a separate boundary closure. Both read the "
        "same Path-A background and B2a BdG JSON inputs and use the same canonical self-energy definitions. The old "
        "Mathematica worker is now labelled transcription-only because it uses the same stencil/grid/quadrature as "
        "the primary path."
    )
    lines.append("")
    lines.append("## Adaptations")
    lines.append("")
    lines.append("- The old hardcoded smoke `R0(w)` and `Z(w)` were replaced by the exported Path-A `A0={A_00,A_r0,A_w0}`, `Z_w`, and frozen geometry.")
    lines.append("- The VSH operator and DtN/self-energy forms are retained from Spike-1/Spike-2; only the background interpolation and forward source changed.")
    lines.append("- The live transfer uses `Sigma=<j,Gj>` and never constructs U/W mixed ports. The U/W formulas remain only as the cited regression identity.")
    lines.append("- The forward source is the canonical decision-05 D3 Frechet source over the B2a BdG-mode response; the B2a overlaps `c_j` carry the chi coupling. It is not a separate live B1-chi adapter-overlap path.")
    lines.append("")
    lines.append("## Files Created Or Modified")
    lines.append("")
    files = [
        "software/stage1_solver/src/stage1_solver/patha_b2b_maxwell.py",
        "software/stage1_solver/mathematica/mt15_03_spike1_vsh_maxwell_operator.wls",
        "software/stage1_solver/mathematica/mt15_04_spike2_transfer_n0.wls",
        "software/stage1_solver/mathematica/mt15_03_patha_b2b_vsh_maxwell_operator.wls",
        "software/stage1_solver/mathematica/mt15_04_patha_b2b_maxwell_transfer.wls",
        "software/stage1_solver/directives/pathA_11_chunk_b2b_maxwell_transfer.md",
        "software/stage1_solver/tests/test_patha_b2b_maxwell.py",
        str(report_path),
        str(bundle_path),
    ]
    for path in files:
        lines.append(f"- `{path}`")
    run_root = bundle_path.parents[1]
    generated_files = sorted(
        path
        for path in run_root.rglob("*")
        if path.is_file()
        and path.parent.name in {"operator", "python", "independent", "mathematica", "bundles"}
    )
    for path in generated_files:
        if str(path) not in files:
            lines.append(f"- `{path}`")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Path-A B2b Maxwell transfer helpers")
    sub = parser.add_subparsers(dest="command", required=True)

    solve = sub.add_parser("solve-python")
    solve.add_argument("--background", type=Path, required=True)
    solve.add_argument("--bdg", type=Path, required=True)
    solve.add_argument("--nr", type=int, required=True)
    solve.add_argument("--nw", type=int, required=True)
    solve.add_argument("--window", type=float, default=DEFAULT_FINAL_WINDOW)
    solve.add_argument("--radial-scale", type=float, default=DEFAULT_FINAL_TRUNCATION)
    solve.add_argument("--mode-count", type=int, default=DEFAULT_MODE_COUNT)
    solve.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    solve.add_argument("--out", type=Path, default=None)
    solve.add_argument("--engine", default="python_primary_second_order")
    solve.add_argument("--discretization", default=None)

    validate = sub.add_parser("validate-report")
    validate.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    validate.add_argument("--report-path", type=Path, default=DEFAULT_REPORT_PATH)
    validate.add_argument("--b2a-run-root", type=Path, default=DEFAULT_B2A_RUN_ROOT)
    validate.add_argument("--tau1", type=float, default=1.0)
    validate.add_argument("--tau2", type=float, default=2.0)
    validate.add_argument(
        "--grids",
        default=";".join(f"{nr},{nw}" for nr, nw in DEFAULT_TRANSFER_GRIDS),
    )
    validate.add_argument(
        "--windows",
        default=";".join(f"{window:.6g}" for window in DEFAULT_WINDOWS),
    )
    validate.add_argument(
        "--truncation-scales",
        default=";".join(f"{scale:.6g}" for scale in DEFAULT_TRUNCATION_SCALES),
    )
    validate.add_argument("--final-grid", type=_parse_grid, default=DEFAULT_FINAL_GRID)
    validate.add_argument("--final-window", type=float, default=DEFAULT_FINAL_WINDOW)
    validate.add_argument("--final-truncation", type=float, default=DEFAULT_FINAL_TRUNCATION)

    generate = sub.add_parser("generate-validation")
    generate.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    generate.add_argument("--b2a-run-root", type=Path, default=DEFAULT_B2A_RUN_ROOT)
    generate.add_argument("--tau1", type=float, default=1.0)
    generate.add_argument("--tau2", type=float, default=2.0)
    generate.add_argument(
        "--grids",
        default=";".join(f"{nr},{nw}" for nr, nw in DEFAULT_TRANSFER_GRIDS),
    )
    generate.add_argument(
        "--windows",
        default=";".join(f"{window:.6g}" for window in DEFAULT_WINDOWS),
    )
    generate.add_argument(
        "--truncation-scales",
        default=";".join(f"{scale:.6g}" for scale in DEFAULT_TRUNCATION_SCALES),
    )
    generate.add_argument("--final-grid", type=_parse_grid, default=DEFAULT_FINAL_GRID)
    generate.add_argument("--final-window", type=float, default=DEFAULT_FINAL_WINDOW)
    generate.add_argument("--final-truncation", type=float, default=DEFAULT_FINAL_TRUNCATION)
    generate.add_argument("--mode-count", type=int, default=DEFAULT_MODE_COUNT)

    args = parser.parse_args()
    if args.command == "solve-python":
        discretization = args.discretization
        if discretization is None:
            discretization = (
                "staggered_high_order"
                if "staggered" in args.engine or "independent" in args.engine
                else "primary_second_order"
            )
        path, packet = solve_python_command(
            background_path=args.background,
            bdg_path=args.bdg,
            nr=args.nr,
            nw=args.nw,
            window=args.window,
            radial_scale=args.radial_scale,
            mode_count=args.mode_count,
            run_root=args.run_root,
            out_path=args.out,
            engine=args.engine,
            discretization=discretization,
        )
        print(f"packet_path: {path}")
        print(f"content_hash: {packet['content_hash']}")
        print(f"tau: {packet['tau']}")
        print(f"grid: {args.nr}x{args.nw}")
        print(f"coefficients: {packet['coefficients']}")
        print(f"max_green_residual: {packet['diagnostics']['max_green_residual']:.12e}")
        return 0
    if args.command == "validate-report":
        result = validate_and_report(
            run_root=args.run_root,
            report_path=args.report_path,
            b2a_run_root=args.b2a_run_root,
            tau1=args.tau1,
            tau2=args.tau2,
            grids=_parse_grid_list(args.grids),
            windows=_parse_float_list(args.windows),
            truncation_scales=_parse_float_list(args.truncation_scales),
            final_grid=args.final_grid,
            final_window=args.final_window,
            final_truncation=args.final_truncation,
        )
        print(f"report_path: {result['report_path']}")
        print(f"bundle_path: {result['bundle_path']}")
        print(f"final_bundle_hash: {result['final_bundle_hash']}")
        print(f"passed: {result['passed']}")
        for gate in result["gates"]:
            print(f"{gate['gate']}: {'PASS' if gate['passed'] else 'FAIL'}")
        return 0 if result["passed"] else 1
    if args.command == "generate-validation":
        result = generate_validation_packets(
            run_root=args.run_root,
            b2a_run_root=args.b2a_run_root,
            tau1=args.tau1,
            tau2=args.tau2,
            grids=_parse_grid_list(args.grids),
            windows=_parse_float_list(args.windows),
            truncation_scales=_parse_float_list(args.truncation_scales),
            final_grid=args.final_grid,
            final_window=args.final_window,
            final_truncation=args.final_truncation,
            mode_count=args.mode_count,
        )
        print(f"generated_packets: {result['count']}")
        for row in result["written"]:
            print(
                f"{row['engine']} {row['tau']} {row['grid']} "
                f"w={row['window']:.6g} rs={row['radial_scale']:.6g}: {row['path']}"
            )
        return 0
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
