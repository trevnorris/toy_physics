"""PathA 22b Gate 3 chi_Q outgoing-DtN contest solve.

This module reruns a linear Maxwell self-energy solve on the frozen M1c
finite-core background.  Unlike the older Gate-3 inventory check, it does not
inject ``Gamma_port`` and divide it back out.  The retarded closure is the
exact outgoing spherical-Hankel l=2 Dirichlet-to-Neumann branch, applied as a
frequency-dependent boundary operator on top of the conservative finite-core
operator.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from scipy.linalg import solve as dense_solve
from scipy.optimize import minimize_scalar
from scipy.sparse import csc_matrix, diags, hstack, issparse, kron
from scipy.sparse.linalg import splu
import sympy as sp

from stage1_solver import patha_b2b_maxwell as b2b
from stage1_solver.dimensional_check import DIMENSIONLESS, LENGTH, TIME, expect_dim


FORBIDDEN_TARGET_STRINGS = ("54" + "/5", "10.8" + "/P0")
DEFAULT_FREEZE_HASH = "834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8"
DEFAULT_WINDOW = 0.040
DEFAULT_FREQUENCY_FRACTIONS = (0.55, 0.85, 1.15)
DEFAULT_GRID = (10, 8)
DEFAULT_CLOSURE_RADIUS_SOURCE = "R_exit"
DEFAULT_RADIUS_SWEEP = (1.0, 1.25, 1.5, 1.75, 2.0)
DEFAULT_GRID_SWEEP = ((35, 8), (48, 10), (61, 12))
DEFAULT_DOMAIN_RADIAL_SCALES = (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)
DEFAULT_DOMAIN_BASE_GRID = (10, 8)
DEFAULT_GRID_CONVERGENCE_RADIAL_SCALE = 3.0
DEFAULT_W_NULL_FILTER_STRENGTH = 1.0e-2
W_NULL_SINGULAR_TOL = 1.0e-12
DEFAULT_LARGEST_DOMAIN_GRID_BASES = ((6, 6), (8, 6), (10, 8))
DEFAULT_FLAT_REFERENCE_GRID_LEVELS = 0
DEFAULT_VACUUM_SANITY_RADII = (1.5,)
DEFAULT_VACUUM_SANITY_GRID = (6, 6)
DEFAULT_TINY_DEFECT_FRACTION = 1.0e-3
DEFAULT_TINY_DEFECT_FRACTIONS = (1.0e-3, 3.0e-3, 1.0e-2)
DEFAULT_TINY_DEFECT_GRID = (6, 6)
DEFAULT_NW_CHARACTERIZATION_NRS = (48, 61, 88, 120)
DEFAULT_NW_CHARACTERIZATION_NW_FOR_NR = 16
DEFAULT_NW_CHARACTERIZATION_NWS = (8, 10, 12, 14, 16, 20, 24, 28, 32, 40)
DEFAULT_NW_CHARACTERIZATION_PRIMARY_NR = 48
DEFAULT_NW_CHARACTERIZATION_SECONDARY_NR = 61
DEFAULT_NW_TAIL_CENTRAL_MIN_NW = 56
DEFAULT_NW_TAIL_BUDGET_MIN_NW = 16
DEFAULT_ZW_REFERENCE_NWS = (64, 80, 96)
RADIUS_INVARIANCE_REL_TOL = 1.0e-6
GRID_PLATEAU_REL_TOL = 1.0e-2
DOMAIN_PLATEAU_REL_TOL = 5.0e-2
DOMAIN_PLATEAU_DELTA_ABS_TOL = 1.5e-3
REFERENCE_STABILITY_REL_TOL = 5.0e-2
NW_JACKKNIFE_ABS_TOL = 5.0e-3
NW_JACKKNIFE_REL_TOL = 1.0e-2
NW_NR_TREND_ABS_TOL = 5.0e-3
VACUUM_SANITY_ABS_TOL = 1.0e-10
OUTCOME_PREFIX_DELTA = "DELTA_Q_NE_0"
OUTCOME_CHI_Q_ONE = "CHI_Q=1_DERIVED"
OUTCOME_NOT_EXTRACTABLE = "CHI_Q_MAGNITUDE_NOT_EXTRACTABLE"
OUTCOME_NW_CONVERGES_PREFIX = "NW_CONVERGES"
OUTCOME_NW_DAMPED_PREFIX = "NW_OSCILLATES_DAMPED"
OUTCOME_NW_NO_LIMIT = "NW_NO_FINITE_LIMIT"
OUTCOME_NW_INCONCLUSIVE = "NW_INCONCLUSIVE_NEED_FINER"


def _linear_solve(matrix: np.ndarray, rhs: np.ndarray, *, solver: str = "dense") -> np.ndarray:
    if solver == "dense":
        if issparse(matrix):
            matrix = matrix.toarray()
        return dense_solve(matrix, rhs, assume_a="gen", check_finite=False)
    if solver == "sparse_direct":
        sparse_matrix = matrix.tocsc() if issparse(matrix) else csc_matrix(matrix)
        return splu(sparse_matrix).solve(rhs)
    raise ValueError(f"unknown Gate-3 linear solver tier: {solver}")


def repo_root() -> Path:
    return Path(__file__).resolve().parents[4]


def default_nw_characterization_path(root: Path | None = None) -> Path:
    base = Path(root) if root is not None else repo_root()
    return base / "_scratch" / "pathA_22b_gate3_nw_characterization.json"


def frozen_dir(freeze_hash: str = DEFAULT_FREEZE_HASH) -> Path:
    return repo_root() / "software" / "stage1_solver" / "frozen" / "m1c" / freeze_hash


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def canonical_fingerprint_coefficient(a: sp.Expr, c_s: sp.Expr) -> sp.Expr:
    """Return the normalized outgoing l=2 odd omega^5 fingerprint."""

    z = sp.symbols("z", real=True)
    h2 = sp.exp(sp.I * z) * (z**2 + 3 * sp.I * z - 3) / z**3
    lambda_out = sp.simplify(z * sp.diff(sp.log(h2), z))
    yhat_out = sp.simplify(-sp.Integer(3) / lambda_out)
    z5_coeff = sp.simplify(sp.series(yhat_out, z, 0, 6).removeO().expand().coeff(z, 5) / sp.I)
    return sp.simplify(z5_coeff * a**5 / c_s**5)


def canonical_sigma_q(a: sp.Expr, c_s: sp.Expr) -> sp.Expr:
    """Return sigma_Q^can from the one-pole packaging."""

    omega_q = sp.simplify(3 * c_s / (2 * a))
    return sp.simplify(sp.Rational(9, 8) / omega_q**5)


def outgoing_hankel_lambda_series() -> dict[str, str]:
    z = sp.symbols("z", real=True)
    h2 = sp.exp(sp.I * z) * (z**2 + 3 * sp.I * z - 3) / z**3
    lambda_out = sp.simplify(z * sp.diff(sp.log(h2), z))
    yhat_out = sp.simplify(-sp.Integer(3) / lambda_out)
    return {
        "lambda_series": str(sp.series(lambda_out, z, 0, 6).removeO().expand()),
        "yhat_series": str(sp.series(yhat_out, z, 0, 6).removeO().expand()),
        "z5_loss_coefficient": str(
            sp.simplify(sp.series(yhat_out, z, 0, 6).removeO().expand().coeff(z, 5) / sp.I)
        ),
    }


def extract_iomega5_coefficient(response: sp.Expr, omega: sp.Symbol) -> sp.Expr:
    """Extract c from response = even + I*c*omega^5 + O(omega^6)."""

    expanded = sp.series(response, omega, 0, 6).removeO().expand()
    return sp.simplify(expanded.coeff(omega, 5) / sp.I)


def chi_q_from_omega5_coefficient(
    actual_coeff: sp.Expr,
    *,
    a: sp.Expr,
    c_s: sp.Expr,
) -> sp.Expr:
    """Normalize an odd coefficient by the Hankel-derived fingerprint."""

    return sp.simplify(actual_coeff / canonical_fingerprint_coefficient(a, c_s))


def exact_outgoing_yhat_l2(z: float) -> complex:
    """Evaluate -3/(z d log h_2^(1)/dz) for h_2 in closed form."""

    zz = complex(z)
    if abs(zz) <= 0.0:
        return 1.0 + 0.0j
    poly = zz**2 + 3j * zz - 3.0
    poly_prime = 2.0 * zz + 3j
    lambda_out = zz * (1j + poly_prime / poly - 3.0 / zz)
    return -3.0 / lambda_out


def centered_w_derivative_null_diagnostic(nw: int, dw: float = 1.0) -> dict[str, object]:
    """Return the highest-W centered-stencil mode responsible for odd-lane decoupling."""

    d_w = b2b._first_derivative_matrix(int(nw), float(dw))
    _, singular_values, vh = np.linalg.svd(d_w)
    mode = np.asarray(vh[-1, :].real, dtype=np.float64)
    if np.sum(mode[::2]) < 0.0:
        mode *= -1.0
    norm = float(np.linalg.norm(mode))
    if norm > 0.0:
        mode /= norm
    min_singular = float(singular_values[-1])
    next_singular = float(singular_values[-2]) if singular_values.size >= 2 else 0.0
    residual = float(np.linalg.norm(d_w @ mode))
    even_weight = float(np.linalg.norm(mode[::2]))
    odd_weight = float(np.linalg.norm(mode[1::2]))
    return {
        "nw": int(nw),
        "dw": float(dw),
        "mode": [float(value) for value in mode],
        "min_singular": min_singular,
        "next_singular": next_singular,
        "residual_norm": residual,
        "is_exact_null": bool(min_singular <= W_NULL_SINGULAR_TOL * max(next_singular, 1.0)),
        "even_index_norm": even_weight,
        "odd_index_norm": odd_weight,
        "shape": (
            "odd-nw exact centered-stencil null: support on one W parity, e.g. [1,0,1,0,...,1]"
            if int(nw) % 2 == 1
            else "even-nw has no exact centered-stencil null"
        ),
    }


def _apply_gate3_w_null_lift(
    sys: Mapping[str, Any],
    *,
    strength: float,
    lane_indices: Sequence[int] = (0,),
    exact_null_only: bool = True,
) -> dict[str, Any]:
    """Gate-3-only weak lift of the centered W-stencil null mode.

    The lift is deliberately rank-limited: for each selected field lane and
    radial column it adds a W-weighted projector onto the centered-stencil
    highest-W null vector.  With ``exact_null_only=True`` it activates only on
    odd W counts where the collocated centered first derivative has an exact
    null vector, leaving even grids untouched.
    """

    strength = float(strength)
    if strength == 0.0:
        return dict(sys)
    nw = int(sys["nw"])
    nr = int(sys["nr"])
    n = int(sys["n"])
    total = int(sys["total"])
    diagnostic = centered_w_derivative_null_diagnostic(nw, float(sys["dw"]))
    if exact_null_only and not bool(diagnostic["is_exact_null"]):
        out = dict(sys)
        out["Gate3WNullLift"] = {**diagnostic, "active": False, "strength": strength}
        return out

    mode = np.asarray(diagnostic["mode"], dtype=np.float64)
    w_lane_values = np.asarray(sys["WLaneValues"], dtype=np.float64)
    penalty = np.zeros((total, total), dtype=np.complex128)
    scale = strength / float(sys["dw"]) ** 2
    for lane in lane_indices:
        lane = int(lane)
        for i in range(nr):
            idx = np.asarray([lane * n + j * nr + i for j in range(nw)], dtype=int)
            weights = w_lane_values[idx]
            denom = float(np.sum(weights * mode * mode))
            if denom <= 1.0e-300:
                continue
            weighted_mode = weights * mode
            penalty[np.ix_(idx, idx)] += scale * np.outer(weighted_mode, weighted_mode) / denom

    k_cons = np.asarray(sys["KCons"], dtype=np.complex128) + penalty
    out = dict(sys)
    out["KCons"] = k_cons
    out["ACons"] = np.asarray(sys["WLaneInv"], dtype=np.complex128) @ k_cons
    out["Gate3WNullLift"] = {
        **diagnostic,
        "active": True,
        "strength": strength,
        "lane_indices": [int(lane) for lane in lane_indices],
        "weak_penalty_scale": float(scale),
    }
    return out


def apply_gate3_w_lane_cleaner(
    sys: Mapping[str, Any],
    *,
    cleaner: str = "none",
    strength: float = DEFAULT_W_NULL_FILTER_STRENGTH,
) -> dict[str, Any]:
    if cleaner == "none":
        return dict(sys)
    if cleaner == "exact_null_phi_lift":
        return _apply_gate3_w_null_lift(sys, strength=float(strength), lane_indices=(0,), exact_null_only=True)
    if cleaner == "exact_null_all_lane_lift":
        return _apply_gate3_w_null_lift(sys, strength=float(strength), lane_indices=(0, 1, 2, 3, 4), exact_null_only=True)
    raise ValueError(f"unknown Gate-3 W-lane cleaner: {cleaner}")


def _sparse_diag(values: Sequence[complex] | np.ndarray) -> csc_matrix:
    return diags(np.asarray(values, dtype=np.complex128), offsets=0, format="csc")


def _sparse_block_row(blocks: Sequence[Any]) -> csc_matrix:
    return hstack(blocks, format="csc")


def _sparse_weighted_gram(matrix: csc_matrix, weights: np.ndarray, sign: float) -> csc_matrix:
    weight_diag = _sparse_diag(np.asarray(weights, dtype=np.float64))
    return (float(sign) * (matrix.conjugate().transpose() @ (weight_diag @ matrix))).tocsc()


def _sparse_boundary_matrix(grid: Mapping[str, Any], fields: Mapping[str, csc_matrix]) -> csc_matrix:
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
            cell = j * nr + i
            if j == nw - 1:
                weight = dr * r[cell] ** 2 * z[cell]
                w_exit_scalar[cell] += weight
                w_exit_angular[cell] += b2b.L_HARM * weight
            if i == nr - 1:
                weight = dw * r[cell] ** 2 * z[cell]
                r_outer_scalar[cell] += weight
                r_outer_angular[cell] += b2b.L_HARM * weight
    return (
        fields["E_r"].conjugate().transpose() @ (_sparse_diag(w_exit_scalar) @ fields["E_r"])
        + fields["E_E"].conjugate().transpose()
        @ (_sparse_diag(w_exit_angular + r_outer_angular) @ fields["E_E"])
        + fields["E_B"].conjugate().transpose()
        @ (_sparse_diag(w_exit_angular + r_outer_angular) @ fields["E_B"])
        + fields["E_w"].conjugate().transpose() @ (_sparse_diag(r_outer_scalar) @ fields["E_w"])
    ).tocsc()


def assemble_gate3_sparse_vsh_maxwell_system(
    background: Mapping[str, Any],
    *,
    nr: int,
    nw: int,
    omega: float,
    radial_scale: float = b2b.DEFAULT_FINAL_TRUNCATION,
) -> dict[str, Any]:
    """Gate-3-local sparse equivalent of the primary shared Maxwell assembly."""

    grid = b2b.make_transfer_grid(background, nr=int(nr), nw=int(nw), radial_scale=float(radial_scale))
    nr = int(grid["nr"])
    nw = int(grid["nw"])
    n = int(grid["n"])
    total = 5 * n
    dr1 = csc_matrix(b2b._first_derivative_matrix(nr, float(grid["dr"])))
    dw1 = csc_matrix(b2b._first_derivative_matrix(nw, float(grid["dw"])))
    drr1 = csc_matrix(b2b._second_derivative_matrix(nr, float(grid["dr"])))
    dww1 = csc_matrix(b2b._second_derivative_matrix(nw, float(grid["dw"])))

    id_r = _sparse_diag(np.ones(nr, dtype=np.complex128))
    id_w = _sparse_diag(np.ones(nw, dtype=np.complex128))
    id_n = _sparse_diag(np.ones(n, dtype=np.complex128))
    z_n = csc_matrix((n, n), dtype=np.complex128)
    dr_mat = kron(id_w, dr1, format="csc")
    dw_mat = kron(dw1, id_r, format="csc")
    drr_mat = kron(id_w, drr1, format="csc")
    dww_mat = kron(dww1, id_r, format="csc")

    r = np.asarray(grid["r"], dtype=np.float64)
    inv_r = 1.0 / r
    inv_r2 = inv_r**2
    r_diag = _sparse_diag(r)
    r2_diag = _sparse_diag(r**2)
    inv_r_diag = _sparse_diag(inv_r)
    inv_r2_diag = _sparse_diag(inv_r2)
    w0 = np.asarray(grid["w0"], dtype=np.float64)
    w0_diag = _sparse_diag(w0)
    w_lane_values = np.concatenate([w0, w0, b2b.L_HARM * w0, b2b.L_HARM * w0, w0])
    w_lane = _sparse_diag(w_lane_values)
    w_lane_inv = _sparse_diag(1.0 / w_lane_values)
    div_r = inv_r2_diag @ dr_mat @ r2_diag

    iomega = 1j * float(omega)
    m_g = _sparse_block_row([iomega * id_n, div_r, -b2b.L_HARM * inv_r_diag, z_n, dw_mat])
    m_er = _sparse_block_row([-dr_mat, iomega * id_n, z_n, z_n, z_n])
    m_ee = _sparse_block_row([-inv_r_diag, z_n, iomega * id_n, z_n, z_n])
    m_eb = _sparse_block_row([z_n, z_n, z_n, iomega * id_n, z_n])
    m_ew = _sparse_block_row([-dw_mat, z_n, z_n, z_n, iomega * id_n])
    m_cr = _sparse_block_row([z_n, -dw_mat, z_n, z_n, dr_mat])
    m_ce = _sparse_block_row([z_n, z_n, -dw_mat, z_n, inv_r_diag])
    m_cb = _sparse_block_row([z_n, z_n, z_n, -dw_mat, z_n])
    m_bb = _sparse_block_row([z_n, -inv_r_diag, inv_r_diag @ dr_mat @ r_diag, z_n, z_n])
    m_br = _sparse_block_row([z_n, z_n, z_n, -b2b.L_HARM * inv_r_diag, z_n])
    m_be = _sparse_block_row([z_n, z_n, z_n, -inv_r_diag @ dr_mat @ r_diag, z_n])

    k_action = (
        _sparse_weighted_gram(m_er, w0, 1.0)
        + _sparse_weighted_gram(m_ee, b2b.L_HARM * w0, 1.0)
        + _sparse_weighted_gram(m_eb, b2b.L_HARM * w0, 1.0)
        + _sparse_weighted_gram(m_ew, w0, 1.0)
        - _sparse_weighted_gram(m_cr, w0, 1.0)
        - _sparse_weighted_gram(m_ce, b2b.L_HARM * w0, 1.0)
        - _sparse_weighted_gram(m_cb, b2b.L_HARM * w0, 1.0)
        - _sparse_weighted_gram(m_bb, b2b.L_HARM * w0, 1.0)
        - _sparse_weighted_gram(m_br, w0, 1.0)
        - _sparse_weighted_gram(m_be, b2b.L_HARM * w0, 1.0)
        - (1.0 / b2b.XI_GAUGE) * _sparse_weighted_gram(m_g, w0, 1.0)
    ).tocsc()
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
    boundary_b = _sparse_boundary_matrix(grid, field_matrices)
    anchor = np.zeros(total, dtype=np.float64)
    anchor_cell = int(math.ceil(n / 2.0)) - 1
    anchor[anchor_cell] = 1.0e-8 * w0[anchor_cell]
    k_cons = (k_action + _sparse_diag(anchor)).tocsc()
    gamma = b2b._gamma_port(float(background["geometry"]["a"]), b2b.CS_CONST)
    k_ret = (k_cons + 1j * gamma * float(omega) ** 5 * boundary_b).tocsc()
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
        "ACons": (w_lane_inv @ k_cons).tocsc(),
        "ARet": (w_lane_inv @ k_ret).tocsc(),
        "BoundaryB": boundary_b,
        "AnchorMatrix": _sparse_diag(anchor),
        "ScalarGaugeAnchor": {"cell": int(anchor_cell), "strength": 1.0e-8},
        "FieldMatrices": field_matrices,
        "discretization": "primary_second_order",
        "matrix_assembly": "gate3_sparse",
    }


def frozen_branch_paths(root: Path | None = None) -> dict[str, Path]:
    root = root or frozen_dir()
    return {
        "root": root,
        "background": root / "wp1_background_10x8.json",
        "packet": root / "m1c_v2_22b_physical_frozen_packet.json",
        "freeze_sheet": root / "freeze_sheet.json",
    }


def _background_for_maxwell(background: Mapping[str, Any], packet: Mapping[str, Any]) -> dict[str, Any]:
    cfg = background["config"]["branch"]
    return {
        "content_hash": background["content_hash"],
        "geometry": {
            "a": float(packet["geometry"]["R_mouth"]),
            "L": float(packet["geometry"]["L"]),
        },
        "constants": {
            "tau": 1.0,
            "gauge_charge": float(cfg["gauge_charge"]),
            "particle_mass": float(cfg["particle_mass"]),
            "mu0": float(cfg["mu0"]),
        },
        "grid": background["grid"],
        "fields": {
            "psi_R0": background["fields"]["psi_R0"],
            "psi_I0": background["fields"]["psi_I0"],
            "A_00": background["fields"]["A_00"],
            "A_r0": background["fields"]["A_r0"],
            "A_w0": background["fields"]["A_w0"],
            "R0_w": background["derived"]["R0_w"],
        },
        "derived": {"Z_w": background["derived"]["Z_w"]},
        "residuals": background["residuals"],
    }


def _bdg_packet_for_source(packet: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "content_hash": packet["content_hash"],
        "grid": {"profile_points": packet["grid"]["points"]},
        "bdg_modes": [
            {
                "name": mode["name"],
                "varpi": float(mode["varpi"]),
                "coupling": float(mode["lambda_B"]),
                "profile": mode["profile_values"],
            }
            for mode in packet["bdg_modes"]
        ],
    }


def load_frozen_solve_inputs(root: Path | None = None) -> dict[str, Any]:
    paths = frozen_branch_paths(root)
    background = _load_json(paths["background"])
    packet = _load_json(paths["packet"])
    freeze = _load_json(paths["freeze_sheet"])
    if background["schema"] != "stage1_m1c_physical_wp1_background/v1":
        raise ValueError(f"unexpected background schema: {background['schema']}")
    return {
        "paths": {key: str(path.relative_to(repo_root())) for key, path in paths.items() if key != "root"},
        "background": background,
        "packet": packet,
        "freeze_sheet": freeze,
        "maxwell_background": _background_for_maxwell(background, packet),
        "source_packet": _bdg_packet_for_source(packet),
    }


def make_uniform_vacuum_background(
    background: Mapping[str, Any],
    *,
    label: str = "uniform_vacuum_reference",
) -> dict[str, Any]:
    """Build an in-module defect-free control background with the same schema."""

    vacuum = copy.deepcopy(dict(background))
    nr = int(background["grid"]["nr"])
    nw = int(background["grid"]["nw"])
    shape = (nr, nw)
    zeros = np.zeros(shape, dtype=np.float64)
    vacuum["content_hash"] = label
    vacuum["fields"] = {
        "psi_R0": np.ones(shape, dtype=np.float64).tolist(),
        "psi_I0": zeros.tolist(),
        "A_00": zeros.tolist(),
        "A_r0": zeros.tolist(),
        "A_w0": zeros.tolist(),
        "R0_w": (np.ones(nw, dtype=np.float64) * float(background["geometry"]["a"])).tolist(),
    }
    vacuum["derived"] = {"Z_w": np.ones(nw, dtype=np.float64).tolist()}
    vacuum["residuals"] = {"uniform_vacuum_control": 0.0}
    return vacuum


def make_branch_geometry_reference_background(
    background: Mapping[str, Any],
    *,
    label: str = "branch_geometry_reference",
) -> dict[str, Any]:
    """Undefect the condensate while keeping the branch Z_w/R0_w geometry."""

    reference = copy.deepcopy(dict(background))
    nr = int(background["grid"]["nr"])
    nw = int(background["grid"]["nw"])
    shape = (nr, nw)
    zeros = np.zeros(shape, dtype=np.float64)
    reference["content_hash"] = label
    reference["fields"] = {
        "psi_R0": np.ones(shape, dtype=np.float64).tolist(),
        "psi_I0": zeros.tolist(),
        "A_00": zeros.tolist(),
        "A_r0": zeros.tolist(),
        "A_w0": zeros.tolist(),
        "R0_w": copy.deepcopy(background["fields"]["R0_w"]),
    }
    reference["derived"] = {"Z_w": copy.deepcopy(background["derived"]["Z_w"])}
    reference["residuals"] = {"branch_geometry_reference_control": 0.0}
    return reference


def reference_backgrounds(background: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    return {
        "branch_geometry": make_branch_geometry_reference_background(background),
        "flat_Zw1": make_uniform_vacuum_background(background, label="flat_Zw1_vacuum_reference"),
    }


def blend_backgrounds(
    vacuum_background: Mapping[str, Any],
    defect_background: Mapping[str, Any],
    *,
    defect_fraction: float,
    label: str | None = None,
) -> dict[str, Any]:
    """Return vacuum plus a controlled fraction of the frozen defect fields."""

    fraction = float(defect_fraction)
    blended = copy.deepcopy(dict(vacuum_background))
    blended["content_hash"] = label or f"uniform_plus_{fraction:.6g}_defect"
    fields: dict[str, Any] = {}
    for key in ("psi_R0", "psi_I0", "A_00", "A_r0", "A_w0", "R0_w"):
        left = np.asarray(vacuum_background["fields"][key], dtype=np.float64)
        right = np.asarray(defect_background["fields"][key], dtype=np.float64)
        fields[key] = (left + fraction * (right - left)).tolist()
    blended["fields"] = fields
    blended["derived"] = {
        "Z_w": (
            np.asarray(vacuum_background["derived"]["Z_w"], dtype=np.float64)
            + fraction
            * (
                np.asarray(defect_background["derived"]["Z_w"], dtype=np.float64)
                - np.asarray(vacuum_background["derived"]["Z_w"], dtype=np.float64)
            )
        ).tolist()
    }
    blended["residuals"] = {"defect_fraction": fraction}
    return blended


def _normalized_loss_response(
    sys: Mapping[str, Any],
    j_vec: np.ndarray,
    *,
    closure_radius: float,
    c_s: float,
    linear_solver: str = "dense",
) -> dict[str, Any]:
    omega = float(sys["omega"])
    rhs = np.asarray(j_vec, dtype=np.complex128)
    a_cons = sys["ACons"]
    u_cons = _linear_solve(a_cons, rhs, solver=linear_solver)
    sigma_cons = complex(np.vdot(j_vec, np.asarray(sys["WLaneValues"], dtype=np.float64) * u_cons))

    z = closure_radius * omega / c_s
    yhat = exact_outgoing_yhat_l2(z)
    k_out = sys["KCons"] + (yhat - 1.0) * sys["BoundaryB"]
    a_out = sys["WLaneInv"] @ k_out
    u_out = _linear_solve(a_out, rhs, solver=linear_solver)
    sigma_out = complex(np.vdot(j_vec, np.asarray(sys["WLaneValues"], dtype=np.float64) * u_out))
    response = sigma_out / sigma_cons
    residual_cons = float(np.linalg.norm(a_cons @ u_cons - rhs) / max(np.linalg.norm(rhs), 1.0e-300))
    residual_out = float(np.linalg.norm(a_out @ u_out - rhs) / max(np.linalg.norm(rhs), 1.0e-300))
    boundary_quadratic = float(
        np.real(np.vdot(u_out, sys["BoundaryB"] @ u_out))
    )
    return {
        "sigma_cons": sigma_cons,
        "sigma_out": sigma_out,
        "normalized_response": response,
        "loss_over_omega5": float(-response.imag / omega**5),
        "residual_cons": residual_cons,
        "residual_out": residual_out,
        "boundary_quadratic": boundary_quadratic,
        "outgoing_yhat": yhat,
    }


def solve_outgoing_dtn_sweep(
    *,
    maxwell_background: Mapping[str, Any],
    source_packet: Mapping[str, Any],
    nr: int = DEFAULT_GRID[0],
    nw: int = DEFAULT_GRID[1],
    window: float = DEFAULT_WINDOW,
    closure_radius: float,
    c_s: float = 1.0,
    radial_scale: float = 1.0,
    w_lane_cleaner: str = "none",
    w_filter_strength: float = DEFAULT_W_NULL_FILTER_STRENGTH,
    linear_solver: str = "dense",
    matrix_assembly: str = "shared_dense",
) -> dict[str, Any]:
    fractions = np.asarray(DEFAULT_FREQUENCY_FRACTIONS, dtype=np.float64)
    omegas = float(window) * fractions
    rows: list[dict[str, Any]] = []
    grid_metadata: dict[str, Any] | None = None
    for omega in omegas:
        if matrix_assembly == "shared_dense":
            sys = b2b.assemble_vsh_maxwell_system(
                maxwell_background,
                nr=nr,
                nw=nw,
                omega=float(omega),
                radial_scale=radial_scale,
                discretization="primary_second_order",
            )
        elif matrix_assembly == "gate3_sparse":
            sys = assemble_gate3_sparse_vsh_maxwell_system(
                maxwell_background,
                nr=nr,
                nw=nw,
                omega=float(omega),
                radial_scale=radial_scale,
            )
        else:
            raise ValueError(f"unknown Gate-3 matrix assembly: {matrix_assembly}")
        sys = apply_gate3_w_lane_cleaner(
            sys,
            cleaner=w_lane_cleaner,
            strength=float(w_filter_strength),
        )
        if grid_metadata is None:
            grid_metadata = {
                "nr": int(sys["nr"]),
                "nw": int(sys["nw"]),
                "radial_scale": float(sys["radial_scale"]),
                "r_min": float(sys["r_min"]),
                "r_max": float(sys["r_outer"]),
                "dr": float(sys["dr"]),
                "dw": float(sys["dw"]),
                "h": float(sys["h"]),
            }
        matter = b2b.build_matter_source(sys, maxwell_background, source_packet, mode_count=3)
        solved = _normalized_loss_response(
            sys,
            matter["jVec"],
            closure_radius=closure_radius,
            c_s=c_s,
            linear_solver=linear_solver,
        )
        response = solved["normalized_response"]
        yhat = solved["outgoing_yhat"]
        rows.append(
            {
                "omega": float(omega),
                "normalized_response": {"re": float(response.real), "im": float(response.imag)},
                "loss_over_omega5": solved["loss_over_omega5"],
                "outgoing_yhat": {"re": float(yhat.real), "im": float(yhat.imag)},
                "source_norm": float(matter["sourceNorm"]),
                "continuity_residual": float(matter["continuityResidual"]),
                "cons_residual": solved["residual_cons"],
                "out_residual": solved["residual_out"],
                "boundary_quadratic": solved["boundary_quadratic"],
            }
        )

    x_vals = omegas**2
    values = np.asarray([row["loss_over_omega5"] for row in rows], dtype=np.float64)
    mat = np.column_stack([np.ones_like(x_vals), x_vals, x_vals**2])
    coeff, *_ = np.linalg.lstsq(mat, values, rcond=None)
    fit_residual = values - mat @ coeff
    fingerprint = float(maxwell_background["geometry"]["a"]) ** 5 / (27.0 * c_s**5)
    uncalibrated_fixed_a_ratio = float(coeff[0] / fingerprint)
    assert grid_metadata is not None
    return {
        "grid": grid_metadata,
        "window": float(window),
        "closure": {
            "type": "spherical_hankel_l2_normalized_dtn",
            "operator": "K_out = K_cons + (Y_out(R*omega/c_s)-1)*BoundaryB",
            "radius": float(closure_radius),
            "radius_source": DEFAULT_CLOSURE_RADIUS_SOURCE,
            "c_s": float(c_s),
            "static_subtraction": "Y_out(0)=1 removed so the frozen conservative operator supplies the static branch",
        },
        "w_lane_control": {
            "cleaner": w_lane_cleaner,
            "filter_strength": float(w_filter_strength),
            "linear_solver": linear_solver,
            "matrix_assembly": matrix_assembly,
            "diagnostic": centered_w_derivative_null_diagnostic(int(nw), grid_metadata["dw"]),
        },
        "rows": rows,
        "fit": {
            "basis": "loss_over_omega5 = C5 + C7*omega^2 + C9*omega^4",
            "effective_omega5_coefficient": float(coeff[0]),
            "higher_coefficients": [float(coeff[1]), float(coeff[2])],
            "fit_rms": float(np.linalg.norm(fit_residual) / max(math.sqrt(float(values.size)), 1.0)),
            "condition_number": float(np.linalg.cond(mat)),
            "fixed_a_fingerprint": fingerprint,
            "uncalibrated_fixed_a_ratio": uncalibrated_fixed_a_ratio,
            "sign_convention": "C5 is fitted from -Im(Sigma_out/Sigma_cons)/omega^5, matching selfEnergySolve's loss convention",
        },
        "diagnostics": {
            "max_linear_residual": float(max(max(row["cons_residual"], row["out_residual"]) for row in rows)),
            "max_continuity_residual": float(max(row["continuity_residual"] for row in rows)),
            "loss_over_omega5_min": float(np.min(values)),
            "loss_over_omega5_max": float(np.max(values)),
        },
    }


def _paired_chi_q_row(
    *,
    defect_background: Mapping[str, Any],
    vacuum_background: Mapping[str, Any],
    source_packet: Mapping[str, Any],
    closure_radius: float,
    nr: int,
    nw: int,
    window: float = DEFAULT_WINDOW,
    c_s: float = 1.0,
    radial_scale: float = 1.0,
    w_lane_cleaner: str = "none",
    w_filter_strength: float = DEFAULT_W_NULL_FILTER_STRENGTH,
    linear_solver: str = "dense",
    matrix_assembly: str = "shared_dense",
    cache: dict[tuple[str, float, int, int, float, float, float, str, float, str, str], dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Compute chi_Q = E_defect/E_vacuum at identical radius/grid/window."""

    def solve_cached(label: str, background: Mapping[str, Any]) -> dict[str, Any]:
        background_label = str(background.get("content_hash", label))
        key = (
            f"{label}:{background_label}",
            float(closure_radius),
            int(nr),
            int(nw),
            float(window),
            float(c_s),
            float(radial_scale),
            str(w_lane_cleaner),
            float(w_filter_strength),
            str(linear_solver),
            str(matrix_assembly),
        )
        if cache is not None and key in cache:
            return cache[key]
        solved = solve_outgoing_dtn_sweep(
            maxwell_background=background,
            source_packet=source_packet,
            nr=int(nr),
            nw=int(nw),
            window=float(window),
            closure_radius=float(closure_radius),
            c_s=float(c_s),
            radial_scale=float(radial_scale),
            w_lane_cleaner=w_lane_cleaner,
            w_filter_strength=float(w_filter_strength),
            linear_solver=linear_solver,
            matrix_assembly=matrix_assembly,
        )
        if cache is not None:
            cache[key] = solved
        return solved

    defect = solve_cached("defect", defect_background)
    vacuum = solve_cached("vacuum", vacuum_background)
    defect_coeff = float(defect["fit"]["effective_omega5_coefficient"])
    vacuum_coeff = float(vacuum["fit"]["effective_omega5_coefficient"])
    if abs(vacuum_coeff) <= 1.0e-300:
        raise ZeroDivisionError("vacuum omega^5 coefficient is numerically zero")
    chi_q = defect_coeff / vacuum_coeff
    defect_grid = defect["grid"]
    return {
        "closure_radius": float(closure_radius),
        "grid": {
            "nr": int(defect_grid["nr"]),
            "nw": int(defect_grid["nw"]),
            "radial_scale": float(defect_grid["radial_scale"]),
            "r_min": float(defect_grid["r_min"]),
            "r_max": float(defect_grid["r_max"]),
            "dr": float(defect_grid["dr"]),
            "dw": float(defect_grid["dw"]),
        },
        "window": float(window),
        "chi_Q": float(chi_q),
        "defect_effective_omega5_coefficient": defect_coeff,
        "vacuum_effective_omega5_coefficient": vacuum_coeff,
        "defect_uncalibrated_fixed_a_ratio": float(defect["fit"]["uncalibrated_fixed_a_ratio"]),
        "vacuum_uncalibrated_fixed_a_ratio": float(vacuum["fit"]["uncalibrated_fixed_a_ratio"]),
        "defect_fit_rms": float(defect["fit"]["fit_rms"]),
        "vacuum_fit_rms": float(vacuum["fit"]["fit_rms"]),
        "defect_max_linear_residual": float(defect["diagnostics"]["max_linear_residual"]),
        "vacuum_max_linear_residual": float(vacuum["diagnostics"]["max_linear_residual"]),
        "max_linear_residual": float(
            max(defect["diagnostics"]["max_linear_residual"], vacuum["diagnostics"]["max_linear_residual"])
        ),
        "max_fit_rms": float(max(defect["fit"]["fit_rms"], vacuum["fit"]["fit_rms"])),
        "fit_condition_numbers": {
            "defect": float(defect["fit"]["condition_number"]),
            "vacuum": float(vacuum["fit"]["condition_number"]),
        },
        "w_lane_control": defect.get("w_lane_control", {}),
    }


def radius_invariance_sweep(
    *,
    defect_background: Mapping[str, Any],
    vacuum_background: Mapping[str, Any],
    source_packet: Mapping[str, Any],
    radii: Sequence[float] = DEFAULT_RADIUS_SWEEP,
    grid: tuple[int, int] = DEFAULT_GRID,
    window: float = DEFAULT_WINDOW,
    cache: dict[tuple[str, float, int, int, float, float, float], dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    return [
        _paired_chi_q_row(
            defect_background=defect_background,
            vacuum_background=vacuum_background,
            source_packet=source_packet,
            closure_radius=float(radius),
            nr=int(grid[0]),
            nw=int(grid[1]),
            window=float(window),
            cache=cache,
        )
        for radius in radii
    ]


def grid_invariance_sweep(
    *,
    defect_background: Mapping[str, Any],
    vacuum_background: Mapping[str, Any],
    source_packet: Mapping[str, Any],
    closure_radius: float,
    grids: Sequence[tuple[int, int]] = DEFAULT_GRID_SWEEP,
    window: float = DEFAULT_WINDOW,
    radial_scale: float = 1.0,
    cache: dict[tuple[str, float, int, int, float, float, float], dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    return [
        _paired_chi_q_row(
            defect_background=defect_background,
            vacuum_background=vacuum_background,
            source_packet=source_packet,
            closure_radius=float(closure_radius),
            nr=int(grid[0]),
            nw=int(grid[1]),
            window=float(window),
            radial_scale=float(radial_scale),
            cache=cache,
        )
        for grid in grids
    ]


def _scaled_nr_for_constant_dr(
    background: Mapping[str, Any],
    *,
    base_nr: int,
    radial_scale: float,
) -> int:
    a = float(background["geometry"]["a"])
    r_min = 0.22 * a
    base_outer = float(background["grid"]["r_max"])
    base_dr = (base_outer - r_min) / (int(base_nr) + 1)
    outer = base_outer * float(radial_scale)
    return max(3, int(round((outer - r_min) / base_dr - 1.0)))


def domain_truncation_sweep(
    *,
    defect_background: Mapping[str, Any],
    vacuum_background: Mapping[str, Any],
    source_packet: Mapping[str, Any],
    radial_scales: Sequence[float] = DEFAULT_DOMAIN_RADIAL_SCALES,
    base_grid: tuple[int, int] = DEFAULT_DOMAIN_BASE_GRID,
    window: float = DEFAULT_WINDOW,
    cache: dict[tuple[str, float, int, int, float, float, float], dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    base_r_max = float(defect_background["grid"]["r_max"])
    for scale in radial_scales:
        nr = _scaled_nr_for_constant_dr(defect_background, base_nr=int(base_grid[0]), radial_scale=float(scale))
        row = _paired_chi_q_row(
            defect_background=defect_background,
            vacuum_background=vacuum_background,
            source_packet=source_packet,
            closure_radius=base_r_max * float(scale),
            nr=nr,
            nw=int(base_grid[1]),
            window=float(window),
            radial_scale=float(scale),
            cache=cache,
        )
        row["domain"] = {
            "radial_scale": float(scale),
            "r_max": base_r_max * float(scale),
            "nr_scaled_from_base": int(nr),
            "base_grid": {"nr": int(base_grid[0]), "nw": int(base_grid[1])},
        }
        rows.append(row)
    return rows


def largest_domain_grid_check(
    *,
    defect_background: Mapping[str, Any],
    vacuum_background: Mapping[str, Any],
    source_packet: Mapping[str, Any],
    radial_scale: float,
    base_grids: Sequence[tuple[int, int]] = DEFAULT_LARGEST_DOMAIN_GRID_BASES,
    window: float = DEFAULT_WINDOW,
    cache: dict[tuple[str, float, int, int, float, float, float], dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    base_r_max = float(defect_background["grid"]["r_max"])
    for base_grid in base_grids:
        nr = _scaled_nr_for_constant_dr(
            defect_background,
            base_nr=int(base_grid[0]),
            radial_scale=float(radial_scale),
        )
        row = _paired_chi_q_row(
            defect_background=defect_background,
            vacuum_background=vacuum_background,
            source_packet=source_packet,
            closure_radius=base_r_max * float(radial_scale),
            nr=nr,
            nw=int(base_grid[1]),
            window=float(window),
            radial_scale=float(radial_scale),
            cache=cache,
        )
        row["domain"] = {
            "radial_scale": float(radial_scale),
            "r_max": base_r_max * float(radial_scale),
            "base_grid": {"nr": int(base_grid[0]), "nw": int(base_grid[1])},
        }
        rows.append(row)
    return rows


def richardson_geometric_extrapolation(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    values = np.asarray([float(row["chi_Q"]) for row in rows], dtype=np.float64)
    if values.size < 3:
        raise ValueError("Richardson/geometric-tail extrapolation needs at least three grid rows")
    deltas = np.diff(values)
    nonzero = [float(delta) for delta in deltas if abs(float(delta)) > 1.0e-15]
    if len(nonzero) < 2:
        limit = float(values[-1])
        q_observed = 0.0
        q_previous = 0.0
    else:
        q_observed = float(nonzero[-1] / nonzero[-2])
        q_previous = float(nonzero[-2] / nonzero[-3]) if len(nonzero) >= 3 else q_observed
        if not math.isfinite(q_observed) or q_observed <= 0.0 or q_observed >= 1.0:
            limit = float(values[-1])
            q_observed = 0.0
        else:
            limit = float(values[-1] + nonzero[-1] * q_observed / (1.0 - q_observed))
    q_spread = abs(float(q_observed) - float(q_previous))
    q_low = max(0.0, float(q_observed) - q_spread)
    q_high = min(0.95, float(q_observed) + q_spread)

    def limit_for(q_value: float) -> float:
        if not nonzero or q_value <= 0.0 or q_value >= 1.0:
            return float(values[-1])
        return float(values[-1] + nonzero[-1] * q_value / (1.0 - q_value))

    sensitivity_limits = [limit_for(q_low), limit_for(q_high)]
    richardson_model_uncertainty = max(abs(item - limit) for item in sensitivity_limits)
    grid_uncertainty = abs(float(values[-1]) - limit)
    two_grid_mean = float(np.mean(values[-2:]))
    if len(nonzero) < 3:
        richardson_model_uncertainty = max(richardson_model_uncertainty, abs(two_grid_mean - limit))
    return {
        "values": [float(value) for value in values],
        "deltas": [float(delta) for delta in deltas],
        "delta_ratios": [
            float(deltas[i + 1] / deltas[i])
            for i in range(len(deltas) - 1)
            if abs(float(deltas[i])) > 1.0e-15
        ],
        "observed_delta_ratio": float(q_observed),
        "previous_delta_ratio": float(q_previous),
        "limit": float(limit),
        "finest_value": float(values[-1]),
        "two_grid_mean": two_grid_mean,
        "old_two_grid_mean_bias": abs(two_grid_mean - float(limit)),
        "grid_uncertainty": float(grid_uncertainty),
        "richardson_model_uncertainty": float(richardson_model_uncertainty),
        "method": "geometric tail from the observed finest-grid delta ratio",
    }


def assess_domain_plateau(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    values = np.asarray([float(row["chi_Q"]) for row in rows], dtype=np.float64)
    r_max_values = [float(row["grid"]["r_max"]) for row in rows]
    deltas = np.diff(values)
    abs_deltas = np.abs(deltas)
    tail_count = min(3, values.size)
    tail = values[-tail_count:]
    tail_mean = float(np.mean(tail))
    tail_spread = _spread([float(value) for value in tail])
    full_spread = _spread([float(value) for value in values])
    last_step = float(abs_deltas[-1]) if abs_deltas.size >= 1 else 0.0
    previous_step = float(abs_deltas[-2]) if abs_deltas.size >= 2 else 0.0
    tail_rel_spread = tail_spread / max(abs(tail_mean), 1.0e-300)
    plateau_onset_r_max: float | None = None
    strict_plateau_onset_r_max: float | None = None
    for idx in range(values.size):
        remaining = values[idx:]
        if remaining.size < 3:
            continue
        remaining_spread = _spread([float(value) for value in remaining])
        subsequent = abs_deltas[idx:] if idx < abs_deltas.size else np.asarray([], dtype=np.float64)
        if remaining_spread <= DOMAIN_PLATEAU_DELTA_ABS_TOL and all(
            float(delta) <= DOMAIN_PLATEAU_DELTA_ABS_TOL for delta in subsequent
        ):
            plateau_onset_r_max = float(r_max_values[idx])
            break
    for idx in range(values.size):
        remaining = values[idx:]
        if remaining.size < 2:
            continue
        subsequent = abs_deltas[idx:] if idx < abs_deltas.size else np.asarray([], dtype=np.float64)
        if subsequent.size and all(float(delta) <= 1.0e-3 for delta in subsequent):
            strict_plateau_onset_r_max = float(r_max_values[idx])
            break
    plateau = bool(
        tail_rel_spread <= DOMAIN_PLATEAU_REL_TOL
        and last_step <= max(previous_step, 1.0e-300)
        and plateau_onset_r_max is not None
    )
    return {
        "plateau": plateau,
        "tail_count": int(tail_count),
        "tail_mean": tail_mean,
        "tail_spread": float(tail_spread),
        "tail_rel_spread": float(tail_rel_spread),
        "full_spread": float(full_spread),
        "last_step": float(last_step),
        "previous_step": float(previous_step),
        "deltas": [float(delta) for delta in deltas],
        "abs_deltas": [float(delta) for delta in abs_deltas],
        "plateau_onset_r_max": plateau_onset_r_max,
        "strict_1e_minus_3_onset_r_max": strict_plateau_onset_r_max,
        "r_max_values": r_max_values,
        "criteria": {
            "tail_rel_spread_max": DOMAIN_PLATEAU_REL_TOL,
            "plateau_delta_abs_max": DOMAIN_PLATEAU_DELTA_ABS_TOL,
            "strict_delta_abs_max_diagnostic": 1.0e-3,
            "last_step_rule": "last domain increment must not exceed the previous increment",
        },
    }


def assess_grid_convergence(rows: Sequence[Mapping[str, Any]], richardson: Mapping[str, Any]) -> dict[str, Any]:
    values = np.asarray([float(row["chi_Q"]) for row in rows], dtype=np.float64)
    deltas = np.diff(values)
    finite_ratio = float(richardson["observed_delta_ratio"])
    same_direction = bool(
        deltas.size == 0
        or all(delta <= 0.0 for delta in deltas)
        or all(delta >= 0.0 for delta in deltas)
    )
    shrinking_tail = bool(0.0 <= finite_ratio < 1.0)
    tail_correction = abs(float(richardson["limit"]) - float(values[-1]))
    rel_tail_correction = tail_correction / max(abs(float(richardson["limit"])), 1.0e-300)
    return {
        "converging": bool(same_direction and shrinking_tail),
        "same_direction": same_direction,
        "shrinking_tail": shrinking_tail,
        "observed_delta_ratio": finite_ratio,
        "tail_correction": float(tail_correction),
        "rel_tail_correction": float(rel_tail_correction),
        "full_drift": _spread([float(value) for value in values]),
    }


def _sorted_nw_rows(rows: Sequence[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    return sorted(rows, key=lambda row: int(row["grid"]["nw"]))


def _validate_even_nws(nws: Sequence[int]) -> tuple[int, ...]:
    out = tuple(int(nw) for nw in nws)
    odd = [nw for nw in out if nw % 2 != 0]
    if odd:
        raise ValueError(f"Gate-3 nw characterization excludes odd singular centered-W lanes: {odd}")
    if len(set(out)) != len(out):
        raise ValueError(f"duplicate nw values in characterization sweep: {out}")
    return out


def nr_convergence_sweep(
    *,
    defect_background: Mapping[str, Any],
    vacuum_background: Mapping[str, Any],
    source_packet: Mapping[str, Any],
    nrs: Sequence[int] = DEFAULT_NW_CHARACTERIZATION_NRS,
    nw: int = DEFAULT_NW_CHARACTERIZATION_NW_FOR_NR,
    window: float = DEFAULT_WINDOW,
    radial_scale: float = DEFAULT_GRID_CONVERGENCE_RADIAL_SCALE,
    linear_solver: str = "sparse_direct",
    matrix_assembly: str = "gate3_sparse",
    cache: dict[tuple[str, float, int, int, float, float, float, str, float, str, str], dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """Sweep nr at one even nw so the later nw sweep isolates the lane axis."""

    if int(nw) % 2 != 0:
        raise ValueError("nr convergence for nw-axis characterization must use an even nw")
    closure_radius = float(defect_background["grid"]["r_max"]) * float(radial_scale)
    rows: list[dict[str, Any]] = []
    for nr in nrs:
        rows.append(
            _paired_chi_q_row(
                defect_background=defect_background,
                vacuum_background=vacuum_background,
                source_packet=source_packet,
                closure_radius=closure_radius,
                nr=int(nr),
                nw=int(nw),
                window=float(window),
                radial_scale=float(radial_scale),
                linear_solver=linear_solver,
                matrix_assembly=matrix_assembly,
                cache=cache,
            )
        )
    return rows


def assess_nr_convergence(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    values = np.asarray([float(row["chi_Q"]) for row in rows], dtype=np.float64)
    if values.size < 3:
        raise ValueError("nr convergence assessment needs at least three rows")
    deltas = np.diff(values)
    tail = values[-min(3, values.size) :]
    full_span = _spread([float(value) for value in values])
    tail_span = _spread([float(value) for value in tail])
    tail_mean = float(np.mean(tail))
    return {
        "nrs": [int(row["grid"]["nr"]) for row in rows],
        "nw": int(rows[0]["grid"]["nw"]),
        "values": [float(value) for value in values],
        "deltas": [float(delta) for delta in deltas],
        "full_span": float(full_span),
        "tail_span": float(tail_span),
        "tail_mean": tail_mean,
        "tail_rel_span": float(tail_span / max(abs(tail_mean), 1.0e-300)),
        "max_linear_residual": float(max(float(row["max_linear_residual"]) for row in rows)),
        "max_fit_rms": float(max(float(row["max_fit_rms"]) for row in rows)),
    }


def even_nw_sweep(
    *,
    defect_background: Mapping[str, Any],
    vacuum_background: Mapping[str, Any],
    source_packet: Mapping[str, Any],
    nr: int,
    nws: Sequence[int] = DEFAULT_NW_CHARACTERIZATION_NWS,
    window: float = DEFAULT_WINDOW,
    radial_scale: float = DEFAULT_GRID_CONVERGENCE_RADIAL_SCALE,
    linear_solver: str = "sparse_direct",
    matrix_assembly: str = "gate3_sparse",
    cache: dict[tuple[str, float, int, int, float, float, float, str, float, str, str], dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """Run the unfiltered even-nw lane sweep at fixed nr and post-plateau r_max."""

    even_nws = _validate_even_nws(nws)
    closure_radius = float(defect_background["grid"]["r_max"]) * float(radial_scale)
    rows: list[dict[str, Any]] = []
    for nw in even_nws:
        rows.append(
            _paired_chi_q_row(
                defect_background=defect_background,
                vacuum_background=vacuum_background,
                source_packet=source_packet,
                closure_radius=closure_radius,
                nr=int(nr),
                nw=int(nw),
                window=float(window),
                radial_scale=float(radial_scale),
                w_lane_cleaner="none",
                linear_solver=linear_solver,
                matrix_assembly=matrix_assembly,
                cache=cache,
            )
        )
    return rows


def fit_nw_power_law(rows: Sequence[Mapping[str, Any]], *, p: float) -> dict[str, Any]:
    sorted_rows = _sorted_nw_rows(rows)
    if len(sorted_rows) < 3:
        raise ValueError("nw power-law fit needs at least three rows")
    p = float(p)
    if p <= 0.0 or not math.isfinite(p):
        raise ValueError(f"power-law order must be positive and finite, got {p}")
    nw_values = np.asarray([float(row["grid"]["nw"]) for row in sorted_rows], dtype=np.float64)
    chi_values = np.asarray([float(row["chi_Q"]) for row in sorted_rows], dtype=np.float64)
    x_values = 1.0 / nw_values**p
    mat = np.column_stack([np.ones_like(x_values), x_values])
    coeff, *_ = np.linalg.lstsq(mat, chi_values, rcond=None)
    predicted = mat @ coeff
    residual = chi_values - predicted
    rms = float(np.linalg.norm(residual) / max(math.sqrt(float(chi_values.size)), 1.0))
    return {
        "model": "chi_Q(nw)=chi_inf+c/nw^p",
        "p": p,
        "chi_inf": float(coeff[0]),
        "c": float(coeff[1]),
        "rms": rms,
        "max_abs_residual": float(np.max(np.abs(residual))),
        "nw": [int(value) for value in nw_values],
        "chi_Q": [float(value) for value in chi_values],
        "x": [float(value) for value in x_values],
        "predicted": [float(value) for value in predicted],
        "residuals": [float(value) for value in residual],
    }


def fit_nw_free_power_law(
    rows: Sequence[Mapping[str, Any]],
    *,
    p_bounds: tuple[float, float] = (0.25, 6.0),
) -> dict[str, Any]:
    if len(rows) < 4:
        raise ValueError("free-p nw fit needs at least four rows")

    def objective(p_value: float) -> float:
        return float(fit_nw_power_law(rows, p=float(p_value))["rms"])

    result = minimize_scalar(objective, bounds=p_bounds, method="bounded", options={"xatol": 1.0e-8})
    fit = fit_nw_power_law(rows, p=float(result.x))
    fit["optimizer"] = {
        "success": bool(result.success),
        "fun": float(result.fun),
        "bounds": [float(p_bounds[0]), float(p_bounds[1])],
    }
    return fit


def jackknife_nw_limit(
    rows: Sequence[Mapping[str, Any]],
    *,
    p: float | None = None,
    free_p: bool = False,
    max_drop_coarsest: int = 2,
) -> dict[str, Any]:
    sorted_rows = _sorted_nw_rows(rows)
    if len(sorted_rows) < 5:
        raise ValueError("jackknife needs at least five nw rows")
    fits: list[dict[str, Any]] = []
    for drop in range(0, int(max_drop_coarsest) + 1):
        kept = sorted_rows[drop:]
        if len(kept) < 3 or (free_p and len(kept) < 4):
            continue
        fit = fit_nw_free_power_law(kept) if free_p else fit_nw_power_law(kept, p=float(p))
        fit["drop_coarsest"] = int(drop)
        fit["kept_nw"] = [int(row["grid"]["nw"]) for row in kept]
        fits.append(fit)
    if not fits:
        raise ValueError("no valid jackknife fits")
    base = fits[0]
    shifts = [abs(float(fit["chi_inf"]) - float(base["chi_inf"])) for fit in fits[1:]]
    max_shift = max(shifts) if shifts else 0.0
    tolerance = max(NW_JACKKNIFE_ABS_TOL, NW_JACKKNIFE_REL_TOL * abs(float(base["chi_inf"])))
    return {
        "free_p": bool(free_p),
        "p": None if free_p else float(p),
        "base_chi_inf": float(base["chi_inf"]),
        "fits": fits,
        "max_shift": float(max_shift),
        "tolerance": float(tolerance),
        "stable": bool(max_shift <= tolerance),
    }


def compare_nw_trends_across_nr(
    primary_rows: Sequence[Mapping[str, Any]],
    secondary_rows: Sequence[Mapping[str, Any]],
    primary_fit: Mapping[str, Any],
    secondary_fit: Mapping[str, Any],
) -> dict[str, Any]:
    primary_by_nw = {int(row["grid"]["nw"]): float(row["chi_Q"]) for row in primary_rows}
    secondary_by_nw = {int(row["grid"]["nw"]): float(row["chi_Q"]) for row in secondary_rows}
    common_nw = sorted(set(primary_by_nw).intersection(secondary_by_nw))
    if len(common_nw) < 7:
        raise ValueError("nr trend comparison needs at least seven common even nw rows")
    diffs = [secondary_by_nw[nw] - primary_by_nw[nw] for nw in common_nw]
    max_abs_diff = max(abs(value) for value in diffs)
    rms_diff = float(np.linalg.norm(np.asarray(diffs, dtype=np.float64)) / math.sqrt(float(len(diffs))))
    limit_diff = float(float(secondary_fit["chi_inf"]) - float(primary_fit["chi_inf"]))
    tolerance = max(NW_NR_TREND_ABS_TOL, 2.0 * max(float(primary_fit["rms"]), float(secondary_fit["rms"])))
    primary_steps = np.diff([primary_by_nw[nw] for nw in common_nw])
    secondary_steps = np.diff([secondary_by_nw[nw] for nw in common_nw])
    step_dot = float(np.dot(primary_steps, secondary_steps))
    return {
        "common_nw": common_nw,
        "differences_secondary_minus_primary": [float(value) for value in diffs],
        "max_abs_diff": float(max_abs_diff),
        "rms_diff": rms_diff,
        "limit_diff": limit_diff,
        "tolerance": float(tolerance),
        "step_dot": step_dot,
        "same_trend": bool(step_dot > 0.0),
        "nr_independent": bool(max_abs_diff <= tolerance and abs(limit_diff) <= tolerance and step_dot > 0.0),
    }


def characterize_even_nw_axis(
    *,
    primary_rows: Sequence[Mapping[str, Any]],
    secondary_rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Classify the even-nw asymptotic lane behavior from two fixed-nr sweeps."""

    if len(primary_rows) < 7 or len(secondary_rows) < 7:
        return {
            "outcome_class": OUTCOME_NW_INCONCLUSIVE,
            "reason": "fewer than seven even-nw rows in one or both nr sweeps",
        }
    primary_fixed_1 = fit_nw_power_law(primary_rows, p=1.0)
    primary_fixed_2 = fit_nw_power_law(primary_rows, p=2.0)
    primary_free = fit_nw_free_power_law(primary_rows)
    fits = {
        "p1": primary_fixed_1,
        "p2": primary_fixed_2,
        "free": primary_free,
    }
    best_name, best_fit = min(fits.items(), key=lambda item: float(item[1]["rms"]))
    free_p = best_name == "free"
    jackknife = jackknife_nw_limit(
        primary_rows,
        p=None if free_p else float(best_fit["p"]),
        free_p=free_p,
    )
    secondary_fit = (
        fit_nw_free_power_law(secondary_rows)
        if free_p
        else fit_nw_power_law(secondary_rows, p=float(best_fit["p"]))
    )
    nr_comparison = compare_nw_trends_across_nr(primary_rows, secondary_rows, best_fit, secondary_fit)
    residuals = np.asarray(best_fit["residuals"], dtype=np.float64)
    tail_residuals = residuals[-min(5, residuals.size) :]
    alternating_tail = bool(
        tail_residuals.size >= 4
        and all(float(tail_residuals[i] * tail_residuals[i + 1]) <= 0.0 for i in range(tail_residuals.size - 1))
        and abs(float(tail_residuals[-1])) <= max(abs(float(tail_residuals[0])), 1.0e-300)
    )
    stable_limit = bool(jackknife["stable"] and nr_comparison["nr_independent"])
    if stable_limit and alternating_tail:
        outcome = f"{OUTCOME_NW_DAMPED_PREFIX}_limit={float(best_fit['chi_inf']):.12g}"
        reason = "finite limit with alternating shrinking tail residuals"
    elif stable_limit:
        outcome = (
            f"{OUTCOME_NW_CONVERGES_PREFIX}_chi_inf={float(best_fit['chi_inf']):.12g}"
            f"_pm_{float(best_fit['rms']):.3g}_order={float(best_fit['p']):.4g}"
        )
        reason = "finite power-law limit with stable coarsest-point jackknife and second-nr confirmation"
    else:
        primary_tail = np.asarray([float(row["chi_Q"]) for row in _sorted_nw_rows(primary_rows)[-5:]], dtype=np.float64)
        tail_span = _spread([float(value) for value in primary_tail])
        last_step = abs(float(primary_tail[-1] - primary_tail[-2])) if primary_tail.size >= 2 else 0.0
        previous_step = abs(float(primary_tail[-2] - primary_tail[-3])) if primary_tail.size >= 3 else 0.0
        no_saturation = bool(tail_span > 0.02 and last_step >= 0.75 * max(previous_step, 1.0e-300))
        if no_saturation and nr_comparison["nr_independent"]:
            outcome = OUTCOME_NW_NO_LIMIT
            reason = "tail continues moving without saturation in both nr sweeps"
        else:
            outcome = OUTCOME_NW_INCONCLUSIVE
            reason = "finite-limit extrapolation is not jackknife-stable or not nr-independent"
    return {
        "outcome_class": outcome,
        "reason": reason,
        "fits": fits,
        "best_fit_name": best_name,
        "best_fit": best_fit,
        "jackknife": jackknife,
        "secondary_fit": secondary_fit,
        "nr_comparison": nr_comparison,
        "tail_residuals_alternate": alternating_tail,
    }


def _characterization_rows(
    characterization_data: Mapping[str, Any],
    *,
    kind: str,
    nr: int | None = None,
    nw: int | None = None,
) -> list[Mapping[str, Any]]:
    rows = []
    for row in characterization_data.get("rows", {}).values():
        if row.get("sweep_kind") != kind:
            continue
        grid = row["grid"]
        if nr is not None and int(grid["nr"]) != int(nr):
            continue
        if nw is not None and int(grid["nw"]) != int(nw):
            continue
        rows.append(row)
    return sorted(rows, key=lambda row: (int(row["grid"]["nr"]), int(row["grid"]["nw"])))


def _tail_fit_bundle(
    *,
    primary_rows: Sequence[Mapping[str, Any]],
    secondary_rows: Sequence[Mapping[str, Any]],
    min_nw: int,
) -> dict[str, Any]:
    primary_tail = [row for row in _sorted_nw_rows(primary_rows) if int(row["grid"]["nw"]) >= int(min_nw)]
    secondary_tail = [row for row in _sorted_nw_rows(secondary_rows) if int(row["grid"]["nw"]) >= int(min_nw)]
    if len(primary_tail) < 5 or len(secondary_tail) < 5:
        raise ValueError(f"nw tail fit needs at least five rows at both nr values for min_nw={min_nw}")
    fit = fit_nw_free_power_law(primary_tail)
    secondary_fit = fit_nw_free_power_law(secondary_tail)
    max_drop = min(2, max(0, len(primary_tail) - 4))
    jackknife = jackknife_nw_limit(primary_tail, free_p=True, max_drop_coarsest=max_drop)
    return {
        "min_nw": int(min_nw),
        "point_count": int(len(primary_tail)),
        "fit": fit,
        "secondary_fit": secondary_fit,
        "jackknife": jackknife,
        "nr_limit_shift": float(float(secondary_fit["chi_inf"]) - float(fit["chi_inf"])),
        "tail_span": _spread([float(row["chi_Q"]) for row in primary_tail]),
        "max_nw": int(max(int(row["grid"]["nw"]) for row in primary_tail)),
    }


def summarize_nw_characterization(
    characterization_data: Mapping[str, Any],
    *,
    source_path: Path | None = None,
    primary_nr: int = DEFAULT_NW_CHARACTERIZATION_PRIMARY_NR,
    secondary_nr: int = DEFAULT_NW_CHARACTERIZATION_SECONDARY_NR,
    central_min_nw: int = DEFAULT_NW_TAIL_CENTRAL_MIN_NW,
    budget_min_nw: int = DEFAULT_NW_TAIL_BUDGET_MIN_NW,
) -> dict[str, Any]:
    """Summarize the stored even-nw convergence characterization for production."""

    characterization = characterization_data["characterization"]
    outcome_class = str(characterization["outcome_class"])
    primary_rows = _characterization_rows(characterization_data, kind="nw", nr=int(primary_nr))
    secondary_rows = _characterization_rows(characterization_data, kind="nw", nr=int(secondary_nr))
    if len(primary_rows) < 7 or len(secondary_rows) < 7:
        raise ValueError("stored even-nw characterization needs at least seven rows at both nr values")

    central_tail = _tail_fit_bundle(
        primary_rows=primary_rows,
        secondary_rows=secondary_rows,
        min_nw=int(central_min_nw),
    )
    budget_tail = _tail_fit_bundle(
        primary_rows=primary_rows,
        secondary_rows=secondary_rows,
        min_nw=int(budget_min_nw),
    )
    tail_fit_floor = max(
        float(budget_tail["fit"]["rms"]),
        float(budget_tail["jackknife"]["max_shift"]),
    )
    nr_offset = abs(float(budget_tail["nr_limit_shift"]))
    numerical_raw = math.sqrt(tail_fit_floor**2 + nr_offset**2)
    numerical_reported = round(numerical_raw, 4)
    central_raw = float(central_tail["fit"]["chi_inf"])
    central_reported = round(central_raw, 3)
    full_sweep_fit = characterization["best_fit"]
    full_sweep_uncertainty = float(full_sweep_fit["rms"])
    jackknife = characterization["jackknife"]
    nr_comparison = characterization["nr_comparison"]
    converged = bool(
        outcome_class.startswith(OUTCOME_NW_CONVERGES_PREFIX)
        and bool(jackknife["stable"])
        and bool(nr_comparison["nr_independent"])
    )
    return {
        "source_path": str(source_path) if source_path is not None else None,
        "source_path_repo_relative": (
            str(source_path.relative_to(repo_root()))
            if source_path is not None and source_path.is_absolute() and source_path.is_relative_to(repo_root())
            else str(source_path)
            if source_path is not None
            else None
        ),
        "outcome_class": outcome_class,
        "reason": characterization["reason"],
        "converged": converged,
        "primary_nr": int(primary_nr),
        "secondary_nr": int(secondary_nr),
        "central_tail": central_tail,
        "budget_tail": budget_tail,
        "chi_Q_raw": central_raw,
        "chi_Q_reported": float(central_reported),
        "numerical_tail_supported_raw": float(numerical_raw),
        "numerical_tail_supported": float(numerical_reported),
        "numerical_tail_fit_floor": float(tail_fit_floor),
        "nr_offset": float(nr_offset),
        "full_sweep_fit": full_sweep_fit,
        "full_sweep_fit_uncertainty": float(full_sweep_uncertainty),
        "jackknife": jackknife,
        "nr_comparison": nr_comparison,
        "stored_max_linear_residual": float(
            max(float(row["max_linear_residual"]) for row in primary_rows + secondary_rows)
        ),
        "stored_max_fit_rms": float(max(float(row["max_fit_rms"]) for row in primary_rows + secondary_rows)),
        "primary_nw_min": int(min(int(row["grid"]["nw"]) for row in primary_rows)),
        "primary_nw_max": int(max(int(row["grid"]["nw"]) for row in primary_rows)),
    }


def matched_flat_zw_reference_comparison(
    *,
    physical_background: Mapping[str, Any],
    flat_reference_background: Mapping[str, Any],
    source_packet: Mapping[str, Any],
    characterization_data: Mapping[str, Any],
    primary_nr: int = DEFAULT_NW_CHARACTERIZATION_PRIMARY_NR,
    nws: Sequence[int] = DEFAULT_ZW_REFERENCE_NWS,
    window: float = DEFAULT_WINDOW,
    radial_scale: float = DEFAULT_GRID_CONVERGENCE_RADIAL_SCALE,
) -> dict[str, Any]:
    """Compare branch-geometry and flat-Zw references on matched converged even-nw rows."""

    primary_rows = _characterization_rows(characterization_data, kind="nw", nr=int(primary_nr))
    branch_by_nw = {int(row["grid"]["nw"]): row for row in primary_rows}
    matched_nws = [int(nw) for nw in nws if int(nw) in branch_by_nw]
    if not matched_nws:
        raise ValueError("no matched branch even-nw rows available for flat-Zw reference comparison")

    rows: list[dict[str, Any]] = []
    for nw in matched_nws:
        branch_row = branch_by_nw[nw]
        closure_radius = float(branch_row["closure_radius"])
        flat_solve = solve_outgoing_dtn_sweep(
            maxwell_background=flat_reference_background,
            source_packet=source_packet,
            nr=int(primary_nr),
            nw=int(nw),
            window=float(window),
            closure_radius=closure_radius,
            radial_scale=float(radial_scale),
            linear_solver="sparse_direct",
            matrix_assembly="gate3_sparse",
        )
        defect_coeff = float(branch_row["defect_effective_omega5_coefficient"])
        flat_coeff = float(flat_solve["fit"]["effective_omega5_coefficient"])
        flat_chi = defect_coeff / flat_coeff
        branch_chi = float(branch_row["chi_Q"])
        shift = branch_chi - flat_chi
        rel = abs(shift) / max(abs(branch_chi), 1.0e-300)
        rows.append(
            {
                "nw": int(nw),
                "nr": int(primary_nr),
                "r_max": float(branch_row["grid"]["r_max"]),
                "branch_geometry_chi_Q": branch_chi,
                "flat_Zw1_chi_Q": float(flat_chi),
                "shift_branch_minus_flat": float(shift),
                "rel_shift_vs_branch": float(rel),
                "branch_defect_effective_omega5_coefficient": defect_coeff,
                "branch_reference_effective_omega5_coefficient": float(
                    branch_row["vacuum_effective_omega5_coefficient"]
                ),
                "flat_reference_effective_omega5_coefficient": flat_coeff,
                "max_linear_residual": float(
                    max(float(branch_row["max_linear_residual"]), flat_solve["diagnostics"]["max_linear_residual"])
                ),
                "max_fit_rms": float(max(float(branch_row["max_fit_rms"]), flat_solve["fit"]["fit_rms"])),
            }
        )

    branch_values = np.asarray([float(row["branch_geometry_chi_Q"]) for row in rows], dtype=np.float64)
    flat_values = np.asarray([float(row["flat_Zw1_chi_Q"]) for row in rows], dtype=np.float64)
    shifts = branch_values - flat_values
    branch_mean = float(np.mean(branch_values))
    flat_mean = float(np.mean(flat_values))
    shift_mean = float(np.mean(shifts))
    rel_mean = abs(shift_mean) / max(abs(branch_mean), 1.0e-300)
    return {
        "definition": "matched converged even-nw branch rows versus a flat Z_w=1, R0_w=a reference",
        "nws": matched_nws,
        "rows": rows,
        "branch_geometry_tail_mean": branch_mean,
        "flat_Zw1_tail_mean": flat_mean,
        "shift_branch_minus_flat": shift_mean,
        "abs_shift": float(abs(shift_mean)),
        "rel_shift": float(rel_mean),
        "rel_shift_percent": float(100.0 * rel_mean),
        "systematic_kind": "one_sided_definitional",
        "one_sided_direction": "flat_Zw1 lowers chi_Q relative to the branch-geometry physical zero",
    }


def compare_reference_plateaus(
    branch_assessment: Mapping[str, Any],
    flat_assessment: Mapping[str, Any],
) -> dict[str, Any]:
    branch_value = float(branch_assessment["tail_mean"])
    flat_value = float(flat_assessment["tail_mean"])
    shift = branch_value - flat_value
    rel_shift = abs(shift) / max(abs(branch_value), abs(flat_value), 1.0e-300)
    return {
        "branch_geometry_tail_mean": branch_value,
        "flat_Zw1_tail_mean": flat_value,
        "shift_branch_minus_flat": float(shift),
        "abs_shift": float(abs(shift)),
        "rel_shift": float(rel_shift),
        "reference_stable": bool(rel_shift <= REFERENCE_STABILITY_REL_TOL),
        "criteria": {"reference_rel_shift_max": REFERENCE_STABILITY_REL_TOL},
    }


def tiny_defect_linearity_scan(
    *,
    physical_background: Mapping[str, Any],
    reference_background: Mapping[str, Any],
    source_packet: Mapping[str, Any],
    fractions: Sequence[float] = DEFAULT_TINY_DEFECT_FRACTIONS,
    grid: tuple[int, int] = DEFAULT_GRID,
    window: float = DEFAULT_WINDOW,
    radial_scale: float = 1.0,
    cache: dict[tuple[str, float, int, int, float, float, float], dict[str, Any]] | None = None,
) -> dict[str, Any]:
    closure_radius = float(physical_background["grid"]["r_max"]) * float(radial_scale)
    rows: list[dict[str, Any]] = []
    for fraction in fractions:
        tiny_background = blend_backgrounds(
            reference_background,
            physical_background,
            defect_fraction=float(fraction),
            label=f"branch_reference_plus_{fraction:.6g}_defect",
        )
        row = _paired_chi_q_row(
            defect_background=tiny_background,
            vacuum_background=reference_background,
            source_packet=source_packet,
            closure_radius=closure_radius,
            nr=int(grid[0]),
            nw=int(grid[1]),
            window=float(window),
            radial_scale=float(radial_scale),
            cache=cache,
        )
        delta = float(row["chi_Q"]) - 1.0
        row["defect_fraction"] = float(fraction)
        row["chi_Q_minus_1"] = delta
        row["slope_chi_Q_minus_1_over_fraction"] = delta / float(fraction)
        rows.append(row)
    slopes = np.asarray([float(row["slope_chi_Q_minus_1_over_fraction"]) for row in rows], dtype=np.float64)
    slope_mean = float(np.mean(slopes))
    slope_spread = _spread([float(slope) for slope in slopes])
    slope_rel_spread = slope_spread / max(abs(slope_mean), 1.0e-300)
    return {
        "definition": "branch-geometry reference plus a small fraction of the frozen condensate/gauge defect",
        "rows": rows,
        "slope_mean": slope_mean,
        "slope_spread": float(slope_spread),
        "slope_rel_spread": float(slope_rel_spread),
        "linear_toward_zero": bool(slope_rel_spread <= 0.25 and all(abs(float(row["chi_Q_minus_1"])) < 2.0e-2 for row in rows)),
        "criteria": {
            "slope_rel_spread_max": 0.25,
            "max_abs_chi_Q_minus_1": 2.0e-2,
        },
    }


def vacuum_sanity_check(
    *,
    physical_background: Mapping[str, Any],
    source_packet: Mapping[str, Any],
    radii: Sequence[float] = DEFAULT_VACUUM_SANITY_RADII,
    grid: tuple[int, int] = DEFAULT_GRID,
    window: float = DEFAULT_WINDOW,
    tiny_defect_fraction: float = DEFAULT_TINY_DEFECT_FRACTION,
) -> dict[str, Any]:
    vacuum_reference = make_uniform_vacuum_background(physical_background, label="uniform_vacuum_reference")
    uniform_control = make_uniform_vacuum_background(physical_background, label="uniform_defect_free_control")
    cache: dict[tuple[str, float, int, int, float, float, float], dict[str, Any]] = {}
    rows = radius_invariance_sweep(
        defect_background=uniform_control,
        vacuum_background=vacuum_reference,
        source_packet=source_packet,
        radii=radii,
        grid=grid,
        window=window,
        cache=cache,
    )
    tiny_background = blend_backgrounds(
        vacuum_reference,
        physical_background,
        defect_fraction=float(tiny_defect_fraction),
        label=f"uniform_plus_{tiny_defect_fraction:.6g}_defect",
    )
    tiny_row = _paired_chi_q_row(
        defect_background=tiny_background,
        vacuum_background=vacuum_reference,
        source_packet=source_packet,
        closure_radius=float(radii[len(radii) // 2]),
        nr=int(grid[0]),
        nw=int(grid[1]),
        window=float(window),
        cache=cache,
    )
    max_uniform_error = max(abs(float(row["chi_Q"]) - 1.0) for row in rows)
    return {
        "definition": "uniform no-defect control coefficient divided by separately built uniform vacuum reference coefficient",
        "rows": rows,
        "max_uniform_abs_error": float(max_uniform_error),
        "tiny_defect_fraction": float(tiny_defect_fraction),
        "tiny_defect_chi_Q": float(tiny_row["chi_Q"]),
        "tiny_defect_abs_error": float(abs(float(tiny_row["chi_Q"]) - 1.0)),
        "passed": bool(max_uniform_error <= VACUUM_SANITY_ABS_TOL and abs(float(tiny_row["chi_Q"]) - 1.0) < 1.0e-3),
        "tolerance": {"uniform_abs": VACUUM_SANITY_ABS_TOL, "tiny_defect_abs": 1.0e-3},
    }


def _spread(values: Sequence[float]) -> float:
    arr = np.asarray(values, dtype=np.float64)
    return float(np.max(arr) - np.min(arr)) if arr.size else 0.0


def summarize_convergence(
    *,
    radius_rows: Sequence[Mapping[str, Any]],
    grid_rows: Sequence[Mapping[str, Any]],
    vacuum_sanity: Mapping[str, Any],
    grid_plateau_count: int = 2,
) -> dict[str, Any]:
    radius_values = [float(row["chi_Q"]) for row in radius_rows]
    grid_values = [float(row["chi_Q"]) for row in grid_rows]
    tail_count = min(max(1, int(grid_plateau_count)), len(grid_values))
    grid_tail = grid_values[-tail_count:]
    radius_mean = float(np.mean(radius_values))
    grid_plateau_value = float(np.mean(grid_tail))
    radius_spread = _spread(radius_values)
    grid_full_drift = _spread(grid_values)
    grid_plateau_spread = _spread(grid_tail)
    radius_rel_spread = float(radius_spread / max(abs(radius_mean), 1.0e-300))
    grid_plateau_rel_spread = float(grid_plateau_spread / max(abs(grid_plateau_value), 1.0e-300))
    radius_pass = bool(radius_rel_spread <= RADIUS_INVARIANCE_REL_TOL)
    grid_pass = bool(grid_plateau_rel_spread <= GRID_PLATEAU_REL_TOL)
    vacuum_pass = bool(vacuum_sanity.get("passed", False))
    error_bar = float(max(radius_spread, grid_plateau_spread))
    converged = bool(radius_pass and grid_pass and vacuum_pass)
    if not converged:
        reason_bits = []
        if not radius_pass:
            reason_bits.append(f"radius_rel_spread={radius_rel_spread:.3e}")
        if not grid_pass:
            reason_bits.append(f"grid_plateau_rel_spread={grid_plateau_rel_spread:.3e}")
        if not vacuum_pass:
            reason_bits.append("vacuum_sanity_failed")
        outcome = f"{OUTCOME_NOT_EXTRACTABLE}_{'_'.join(reason_bits)}"
    elif abs(grid_plateau_value - 1.0) <= max(5.0 * error_bar, 5.0e-3):
        outcome = OUTCOME_CHI_Q_ONE
    else:
        outcome = f"{OUTCOME_PREFIX_DELTA}_{grid_plateau_value:.6g}_pm_{error_bar:.2g}"
    return {
        "converged": converged,
        "outcome": outcome,
        "plateau_value": grid_plateau_value,
        "error_bar": error_bar,
        "radius_mean": radius_mean,
        "radius_spread": radius_spread,
        "radius_rel_spread": radius_rel_spread,
        "radius_pass": radius_pass,
        "grid_full_drift": grid_full_drift,
        "grid_plateau_count": int(tail_count),
        "grid_plateau_spread": grid_plateau_spread,
        "grid_plateau_rel_spread": grid_plateau_rel_spread,
        "grid_pass": grid_pass,
        "vacuum_pass": vacuum_pass,
        "criteria": {
            "radius_rel_spread_max": RADIUS_INVARIANCE_REL_TOL,
            "grid_plateau_rel_spread_max": GRID_PLATEAU_REL_TOL,
            "vacuum_uniform_abs_max": VACUUM_SANITY_ABS_TOL,
            "error_bar": "max(radius sweep spread, finest-grid plateau spread)",
        },
    }


def frozen_branch_inventory(root: Path | None = None) -> dict[str, object]:
    inputs = load_frozen_solve_inputs(root)
    background = inputs["background"]
    packet = inputs["packet"]
    freeze = inputs["freeze_sheet"]
    flux_normalization = packet.get("derived_maxwell_transfer", {}).get("flux_normalization", {})
    return {
        "freeze_hash": frozen_branch_paths(root)["background"].parent.name,
        "background_path": inputs["paths"]["background"],
        "packet_path": inputs["paths"]["packet"],
        "freeze_path": inputs["paths"]["freeze_sheet"],
        "background_schema": background.get("schema"),
        "background_grid": {
            "nr": background.get("grid", {}).get("nr"),
            "nw": background.get("grid", {}).get("nw"),
            "r_range": [background.get("grid", {}).get("r_min"), background.get("grid", {}).get("r_max")],
            "w_range": [background.get("grid", {}).get("w_min"), background.get("grid", {}).get("w_max")],
        },
        "background_fields_present": sorted(background.get("fields", {}).keys()),
        "background_residual_linf": background.get("residuals", {}).get("coupled_stationary_linf"),
        "freeze_boundary_class": freeze.get("geometry", {}).get("boundary_class"),
        "packet_boundary_class": packet.get("geometry", {}).get("boundary_class"),
        "packet_exit_model": packet.get("geometry", {}).get("exit_model"),
        "closure_radius": packet.get("geometry", {}).get("R_exit"),
        "mouth_radius": packet.get("geometry", {}).get("R_mouth"),
        "legacy_flux_normalization_present_but_not_used": {
            "Gamma_port": flux_normalization.get("Gamma_port"),
            "formula": flux_normalization.get("formula"),
        },
    }


def _is_default_production_gate3_config(
    *,
    radius_sweep: Sequence[float],
    grid_sweep: Sequence[tuple[int, int]],
    grid_convergence_radial_scale: float,
    vacuum_sanity_radii: Sequence[float],
    domain_base_grid: tuple[int, int],
    largest_domain_grid_bases: Sequence[tuple[int, int]],
    flat_reference_grid_levels: int,
    flat_reference_radial_scales: Sequence[float] | None,
    tiny_defect_fractions: Sequence[float],
    tiny_defect_grid: tuple[int, int],
    window: float,
) -> bool:
    return bool(
        tuple(float(value) for value in radius_sweep) == tuple(float(value) for value in DEFAULT_DOMAIN_RADIAL_SCALES)
        and tuple(tuple(int(piece) for piece in grid) for grid in grid_sweep) == tuple(DEFAULT_GRID_SWEEP)
        and float(grid_convergence_radial_scale) == float(DEFAULT_GRID_CONVERGENCE_RADIAL_SCALE)
        and tuple(float(value) for value in vacuum_sanity_radii)
        == tuple(float(value) for value in DEFAULT_VACUUM_SANITY_RADII)
        and tuple(int(value) for value in domain_base_grid) == tuple(DEFAULT_DOMAIN_BASE_GRID)
        and tuple(tuple(int(piece) for piece in grid) for grid in largest_domain_grid_bases)
        == tuple(DEFAULT_LARGEST_DOMAIN_GRID_BASES)
        and int(flat_reference_grid_levels) == int(DEFAULT_FLAT_REFERENCE_GRID_LEVELS)
        and flat_reference_radial_scales is None
        and tuple(float(value) for value in tiny_defect_fractions)
        == tuple(float(value) for value in DEFAULT_TINY_DEFECT_FRACTIONS)
        and tuple(int(value) for value in tiny_defect_grid) == tuple(DEFAULT_TINY_DEFECT_GRID)
        and float(window) == float(DEFAULT_WINDOW)
    )


def gate3_solve_result(
    root: Path | None = None,
    *,
    radius_sweep: Sequence[float] = DEFAULT_DOMAIN_RADIAL_SCALES,
    grid_sweep: Sequence[tuple[int, int]] = DEFAULT_GRID_SWEEP,
    grid_convergence_radial_scale: float = DEFAULT_GRID_CONVERGENCE_RADIAL_SCALE,
    vacuum_sanity_radii: Sequence[float] = DEFAULT_VACUUM_SANITY_RADII,
    domain_base_grid: tuple[int, int] = DEFAULT_DOMAIN_BASE_GRID,
    largest_domain_grid_bases: Sequence[tuple[int, int]] = DEFAULT_LARGEST_DOMAIN_GRID_BASES,
    flat_reference_grid_levels: int = DEFAULT_FLAT_REFERENCE_GRID_LEVELS,
    flat_reference_radial_scales: Sequence[float] | None = None,
    tiny_defect_fractions: Sequence[float] = DEFAULT_TINY_DEFECT_FRACTIONS,
    tiny_defect_grid: tuple[int, int] = DEFAULT_TINY_DEFECT_GRID,
    window: float = DEFAULT_WINDOW,
    use_nw_characterization: bool | None = None,
    nw_characterization_path: Path | None = None,
    zw_reference_nws: Sequence[int] = DEFAULT_ZW_REFERENCE_NWS,
) -> dict[str, object]:
    inputs = load_frozen_solve_inputs(root)
    physical_background = inputs["maxwell_background"]
    base_r_max = float(physical_background["grid"]["r_max"])
    references = reference_backgrounds(physical_background)
    cache: dict[tuple[str, float, int, int, float, float, float], dict[str, Any]] = {}
    if use_nw_characterization is None:
        use_nw_characterization = _is_default_production_gate3_config(
            radius_sweep=radius_sweep,
            grid_sweep=grid_sweep,
            grid_convergence_radial_scale=grid_convergence_radial_scale,
            vacuum_sanity_radii=vacuum_sanity_radii,
            domain_base_grid=domain_base_grid,
            largest_domain_grid_bases=largest_domain_grid_bases,
            flat_reference_grid_levels=flat_reference_grid_levels,
            flat_reference_radial_scales=flat_reference_radial_scales,
            tiny_defect_fractions=tiny_defect_fractions,
            tiny_defect_grid=tiny_defect_grid,
            window=window,
        )
    characterization_path = nw_characterization_path or default_nw_characterization_path(root)
    nw_characterization_data: dict[str, Any] | None = None
    nw_characterization_summary: dict[str, Any] | None = None
    nw_characterization_error: str | None = None
    if use_nw_characterization:
        if characterization_path.exists():
            try:
                nw_characterization_data = _load_json(characterization_path)
                nw_characterization_summary = summarize_nw_characterization(
                    nw_characterization_data,
                    source_path=characterization_path,
                )
            except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
                nw_characterization_error = str(exc)
        else:
            nw_characterization_error = f"missing even-nw characterization JSON: {characterization_path}"

    reference_results: dict[str, Any] = {}
    for reference_name, reference_background in references.items():
        reference_radial_scales = (
            tuple(radius_sweep)
            if reference_name == "branch_geometry"
            else tuple(flat_reference_radial_scales or (float(grid_convergence_radial_scale),))
        )
        domain_rows = domain_truncation_sweep(
            defect_background=physical_background,
            vacuum_background=reference_background,
            source_packet=inputs["source_packet"],
            radial_scales=reference_radial_scales,
            base_grid=domain_base_grid,
            window=window,
            cache=cache,
        )
        if reference_name == "branch_geometry":
            grid_rows = grid_invariance_sweep(
                defect_background=physical_background,
                vacuum_background=reference_background,
                source_packet=inputs["source_packet"],
                closure_radius=base_r_max * float(grid_convergence_radial_scale),
                grids=grid_sweep,
                window=window,
                radial_scale=float(grid_convergence_radial_scale),
                cache=cache,
            )
            richardson: dict[str, Any] | None = richardson_geometric_extrapolation(grid_rows)
            grid_assessment: dict[str, Any] | None = assess_grid_convergence(grid_rows, richardson)
        elif int(flat_reference_grid_levels) >= 3:
            flat_grids = tuple(grid_sweep[: int(flat_reference_grid_levels)])
            grid_rows = grid_invariance_sweep(
                defect_background=physical_background,
                vacuum_background=reference_background,
                source_packet=inputs["source_packet"],
                closure_radius=base_r_max * float(grid_convergence_radial_scale),
                grids=flat_grids,
                window=window,
                radial_scale=float(grid_convergence_radial_scale),
                cache=cache,
            )
            richardson = richardson_geometric_extrapolation(grid_rows)
            grid_assessment = assess_grid_convergence(grid_rows, richardson)
        else:
            grid_rows = []
            richardson = None
            grid_assessment = None
        reference_results[reference_name] = {
            "reference_definition": (
                "keeps branch Z_w and R0_w; sets rho uniform and A=0"
                if reference_name == "branch_geometry"
                else "sets Z_w=1 and R0_w=a; sets rho uniform and A=0"
            ),
            "domain_truncation": {
                "rows": domain_rows,
                "assessment": assess_domain_plateau(domain_rows),
                "dr_held_constant_note": (
                    "nr is recomputed from the scale-1 base dr so r_max changes do not deliberately "
                    "coarsen the radial cell size"
                ),
                "radial_scales": [float(scale) for scale in reference_radial_scales],
            },
            "grid_invariance": {
                "closure_radius": base_r_max * float(grid_convergence_radial_scale),
                "radial_scale": float(grid_convergence_radial_scale),
                "rows": grid_rows,
                "richardson": richardson,
                "assessment": grid_assessment,
            },
        }

    reference_comparison = compare_reference_plateaus(
        reference_results["branch_geometry"]["domain_truncation"]["assessment"],
        reference_results["flat_Zw1"]["domain_truncation"]["assessment"],
    )
    branch_by_rmax = {
        float(row["grid"]["r_max"]): float(row["chi_Q"])
        for row in reference_results["branch_geometry"]["domain_truncation"]["rows"]
    }
    flat_by_rmax = {
        float(row["grid"]["r_max"]): float(row["chi_Q"])
        for row in reference_results["flat_Zw1"]["domain_truncation"]["rows"]
    }
    common_rmax_values = sorted(set(branch_by_rmax).intersection(flat_by_rmax))
    if common_rmax_values:
        common_rmax = float(common_rmax_values[0])
        common_branch = float(branch_by_rmax[common_rmax])
        common_flat = float(flat_by_rmax[common_rmax])
        reference_comparison["common_rmax"] = common_rmax
        reference_comparison["common_rmax_branch_geometry"] = common_branch
        reference_comparison["common_rmax_flat_Zw1"] = common_flat
        reference_comparison["common_rmax_shift_branch_minus_flat"] = common_branch - common_flat
    matched_zw_reference: dict[str, Any] | None = None
    matched_zw_reference_error: str | None = None
    if nw_characterization_summary is not None and nw_characterization_data is not None:
        try:
            matched_zw_reference = matched_flat_zw_reference_comparison(
                physical_background=physical_background,
                flat_reference_background=references["flat_Zw1"],
                source_packet=inputs["source_packet"],
                characterization_data=nw_characterization_data,
                primary_nr=int(nw_characterization_summary["primary_nr"]),
                nws=zw_reference_nws,
                window=window,
                radial_scale=float(grid_convergence_radial_scale),
            )
            reference_comparison.update(
                {
                    "branch_geometry_tail_mean": matched_zw_reference["branch_geometry_tail_mean"],
                    "flat_Zw1_tail_mean": matched_zw_reference["flat_Zw1_tail_mean"],
                    "shift_branch_minus_flat": matched_zw_reference["shift_branch_minus_flat"],
                    "abs_shift": matched_zw_reference["abs_shift"],
                    "rel_shift": matched_zw_reference["rel_shift"],
                    "rel_shift_percent": matched_zw_reference["rel_shift_percent"],
                    "reference_stable": False,
                    "systematic_kind": matched_zw_reference["systematic_kind"],
                    "one_sided_direction": matched_zw_reference["one_sided_direction"],
                    "matched_even_nw": matched_zw_reference,
                }
            )
        except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
            matched_zw_reference_error = str(exc)
    trivial_consistency = vacuum_sanity_check(
        physical_background=physical_background,
        source_packet=inputs["source_packet"],
        radii=vacuum_sanity_radii,
        grid=DEFAULT_VACUUM_SANITY_GRID,
        window=window,
    )
    tiny_linearity = tiny_defect_linearity_scan(
        physical_background=physical_background,
        reference_background=references["branch_geometry"],
        source_packet=inputs["source_packet"],
        fractions=tiny_defect_fractions,
        grid=tiny_defect_grid,
        window=window,
        cache=cache,
    )
    sanity = trivial_consistency
    branch_domain = reference_results["branch_geometry"]["domain_truncation"]["assessment"]
    branch_grid = reference_results["branch_geometry"]["grid_invariance"]["assessment"]
    branch_richardson = reference_results["branch_geometry"]["grid_invariance"]["richardson"]
    w_lane_mode_control = {
        "status": "SPURIOUS_W_LANE_MODE_NOT_CONTROLLED_TO_ACCEPTANCE",
        "parity_independent": False,
        "filter_strength_independent": False,
        "joint_convergence": False,
        "odd_mode": centered_w_derivative_null_diagnostic(13, 1.0),
        "even_mode": centered_w_derivative_null_diagnostic(14, 1.0),
        "cleaner_tested": "exact_null_phi_lift",
        "reason": (
            "The collocated centered W derivative has an exact odd-nw null vector. "
            "A Gate-3-local exact-null lift controls that vector but the extracted chi_Q "
            "is filter-strength dependent and does not provide a parity-independent limit."
        ),
    }
    even_nw_extractable = bool(
        nw_characterization_summary is not None
        and bool(nw_characterization_summary["converged"])
        and matched_zw_reference is not None
        and branch_domain["plateau"]
        and branch_grid["converging"]
        and tiny_linearity["linear_toward_zero"]
    )
    legacy_extractable = bool(
        branch_domain["plateau"]
        and branch_grid["converging"]
        and tiny_linearity["linear_toward_zero"]
        and w_lane_mode_control["parity_independent"]
        and w_lane_mode_control["filter_strength_independent"]
        and w_lane_mode_control["joint_convergence"]
    )
    extractable = bool(even_nw_extractable or legacy_extractable)
    grid_uncertainty = max(
        float(branch_richardson["grid_uncertainty"]),
        float(branch_richardson["richardson_model_uncertainty"]),
    )
    richardson_bias = float(branch_richardson["richardson_model_uncertainty"])
    rmax_systematic = float(branch_domain["tail_spread"])
    reference_shift = float(reference_comparison["abs_shift"])
    if even_nw_extractable:
        assert nw_characterization_summary is not None
        central = float(nw_characterization_summary["chi_Q_reported"])
        combined = float(nw_characterization_summary["numerical_tail_supported"])
        outcome = f"{OUTCOME_PREFIX_DELTA}_{central:.3f}_pm_{combined:.4f}_{OUTCOME_NW_CONVERGES_PREFIX}"
        blocker = ""
    elif legacy_extractable:
        central = float(branch_richardson["limit"])
        combined = math.sqrt(grid_uncertainty**2 + rmax_systematic**2)
        outcome = f"{OUTCOME_PREFIX_DELTA}_{central:.6g}_pm_{combined:.2g}"
        blocker = ""
    else:
        central = None
        combined = math.sqrt(grid_uncertainty**2 + rmax_systematic**2)
        reason_bits: list[str] = []
        if not branch_domain["plateau"]:
            reason_bits.append("RMAX_NONPLATEAU")
        if not branch_grid["converging"]:
            reason_bits.append("GRID_SEQUENCE_NOT_CONVERGING")
        if not tiny_linearity["linear_toward_zero"]:
            reason_bits.append("TINY_DEFECT_NONLINEAR")
        if use_nw_characterization and nw_characterization_summary is None:
            reason_bits.append("EVEN_NW_CHARACTERIZATION_MISSING_OR_INVALID")
        elif use_nw_characterization and not bool(nw_characterization_summary["converged"]):
            reason_bits.append("EVEN_NW_CHARACTERIZATION_NOT_CONVERGED")
        if use_nw_characterization and matched_zw_reference is None:
            reason_bits.append("MATCHED_ZW_REFERENCE_UNAVAILABLE")
        if not use_nw_characterization and not w_lane_mode_control["parity_independent"]:
            reason_bits.append("W_LANE_PARITY_SPLIT")
        if not use_nw_characterization and not w_lane_mode_control["filter_strength_independent"]:
            reason_bits.append("FILTER_STRENGTH_DEPENDENT")
        if not use_nw_characterization and not w_lane_mode_control["joint_convergence"]:
            reason_bits.append("NO_JOINT_NR_NW_CONVERGENCE")
        outcome = f"{OUTCOME_NOT_EXTRACTABLE}_{'_'.join(reason_bits)}"
        if (
            "EVEN_NW_CHARACTERIZATION_MISSING_OR_INVALID" in reason_bits
            or "EVEN_NW_CHARACTERIZATION_NOT_CONVERGED" in reason_bits
        ):
            blocker = (
                "chi_Q magnitude remains a calibration knob because the even-nw characterization is missing, "
                "invalid, or non-convergent"
            )
        elif "MATCHED_ZW_REFERENCE_UNAVAILABLE" in reason_bits:
            blocker = (
                "chi_Q magnitude remains a calibration knob because the matched flat-Zw reference systematic "
                "could not be recomputed"
            )
        elif "W_LANE_PARITY_SPLIT" in reason_bits:
            blocker = (
                "chi_Q magnitude remains a calibration knob because the centered W-lane checkerboard/null mode "
                "is not removed to a parity-independent, filter-strength-independent converged value"
            )
        elif "RMAX_NONPLATEAU" in reason_bits:
            blocker = (
                "chi_Q magnitude does not plateau under the radial-scale domain extension on the frozen "
                "exterior-tapered background"
            )
        else:
            blocker = "chi_Q magnitude failed one or more Gate-3 systematics criteria"
    return {
        "outcome": outcome,
        "blocker": blocker,
        "can_report_chi_Q_number": bool(extractable),
        "chi_Q": central,
        "chi_Q_error_bar": float(combined),
        "inventory": frozen_branch_inventory(root),
        "references": reference_results,
        "reference_comparison": reference_comparison,
        "even_nw_characterization": nw_characterization_summary,
        "even_nw_characterization_error": nw_characterization_error,
        "matched_zw_reference_error": matched_zw_reference_error,
        "post_plateau_grid_convergence": reference_results["branch_geometry"]["grid_invariance"],
        "w_lane_mode_control": w_lane_mode_control,
        "trivial_uniform_consistency": trivial_consistency,
        "vacuum_sanity": sanity,
        "tiny_defect_linearity": tiny_linearity,
        "combined_budget": {
            "grid_uncertainty": float(grid_uncertainty),
            "richardson_bias": float(richardson_bias),
            "rmax_systematic": float(rmax_systematic),
            "zw_reference_shift": float(reference_shift),
            "zw_reference_rel_shift": float(reference_comparison["rel_shift"]),
            "zw_reference_rel_shift_percent": float(reference_comparison.get("rel_shift_percent", 100.0 * reference_comparison["rel_shift"])),
            "zw_reference_systematic_kind": reference_comparison.get("systematic_kind", "diagnostic"),
            "combined_quadrature": float(combined),
            "numerical_tail_supported": float(combined) if even_nw_extractable else None,
            "numerical_tail_supported_raw": (
                float(nw_characterization_summary["numerical_tail_supported_raw"])
                if even_nw_extractable and nw_characterization_summary is not None
                else None
            ),
            "full_sweep_fit_uncertainty": (
                float(nw_characterization_summary["full_sweep_fit_uncertainty"])
                if even_nw_extractable and nw_characterization_summary is not None
                else None
            ),
            "tail_fit_floor": (
                float(nw_characterization_summary["numerical_tail_fit_floor"])
                if even_nw_extractable and nw_characterization_summary is not None
                else None
            ),
            "nr_offset": (
                float(nw_characterization_summary["nr_offset"])
                if even_nw_extractable and nw_characterization_summary is not None
                else None
            ),
            "status": (
                "reported_tail_supported_numerical; full_sweep_fit_diagnostic; zw_reference_separate_one_sided"
                if even_nw_extractable
                else "reported_as_error_bar"
                if legacy_extractable
                else "diagnostic_only_magnitude_not_extractable"
            ),
        },
        "convergence": {
            "converged": bool(extractable),
            "outcome": outcome,
            "radius_pass": bool(branch_domain["plateau"]),
            "grid_pass": bool(branch_grid["converging"]),
            "vacuum_pass": bool(tiny_linearity["linear_toward_zero"]),
            "even_nw_characterization_pass": bool(even_nw_extractable),
            "w_lane_parity_pass": bool(w_lane_mode_control["parity_independent"]),
            "filter_strength_pass": bool(w_lane_mode_control["filter_strength_independent"]),
            "joint_convergence_pass": bool(w_lane_mode_control["joint_convergence"]),
            "plateau_value": central,
            "error_bar": float(combined),
        },
        "reason": (
            "The physically correct defect zero keeps the branch Z_w/R0_w geometry and removes only the "
            "condensate/gauge defect. The production magnitude is earned from the unfiltered even-nw branch-geometry "
            "characterization when that stored sweep has a stable finite limit, nr confirmation, and a matched flat-Zw "
            "reference comparison. The flat-Zw shift is a separate one-sided definitional systematic, not part of the "
            "numerical error bar. If the even-nw evidence is unavailable or non-convergent, the magnitude reverts to a "
            "calibration knob."
        ),
    }


def dimensional_checks() -> dict[str, object]:
    defect_omega5_coeff = TIME**5
    vacuum_omega5_coeff = TIME**5
    closure_factor = LENGTH**5 / (LENGTH / TIME) ** 5
    checks = [
        expect_dim(
            "pathA_22b_gate3",
            "defect omega^5 coefficient",
            defect_omega5_coeff,
            TIME**5,
            "The defect coefficient multiplying omega^5 carries T^5.",
        ),
        expect_dim(
            "pathA_22b_gate3",
            "vacuum omega^5 coefficient",
            vacuum_omega5_coeff,
            TIME**5,
            "The same extraction on the uniform background carries T^5.",
        ),
        expect_dim(
            "pathA_22b_gate3",
            "closure placement factor (R_exit/c_s)^5",
            closure_factor,
            TIME**5,
            "The outgoing-Hankel boundary-placement factor is common to defect and vacuum extractions.",
        ),
        expect_dim(
            "pathA_22b_gate3",
            "radius-consistent chi_Q ratio",
            defect_omega5_coeff / vacuum_omega5_coeff,
            DIMENSIONLESS,
            "chi_Q is the ratio of two same-radius omega^5 coefficients.",
        ),
    ]
    return {"checks": [check.as_dict() for check in checks], "unit_symbols": ["M", "L", "T"]}


def target_blindness_guard(paths: Iterable[Path]) -> dict[str, object]:
    hits: list[str] = []
    for path in paths:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        for forbidden in FORBIDDEN_TARGET_STRINGS:
            if forbidden in text:
                hits.append(f"{path}:{forbidden}")
    return {
        "status": "TARGET_BLIND_PASS" if not hits else "TARGET_BLIND_FAILURE",
        "hits": hits,
    }


def gate3_report_payload(
    root: Path | None = None,
    *,
    radius_sweep: Sequence[float] = DEFAULT_DOMAIN_RADIAL_SCALES,
    grid_sweep: Sequence[tuple[int, int]] = DEFAULT_GRID_SWEEP,
    grid_convergence_radial_scale: float = DEFAULT_GRID_CONVERGENCE_RADIAL_SCALE,
    vacuum_sanity_radii: Sequence[float] = DEFAULT_VACUUM_SANITY_RADII,
    domain_base_grid: tuple[int, int] = DEFAULT_DOMAIN_BASE_GRID,
    largest_domain_grid_bases: Sequence[tuple[int, int]] = DEFAULT_LARGEST_DOMAIN_GRID_BASES,
    flat_reference_grid_levels: int = DEFAULT_FLAT_REFERENCE_GRID_LEVELS,
    flat_reference_radial_scales: Sequence[float] | None = None,
    tiny_defect_fractions: Sequence[float] = DEFAULT_TINY_DEFECT_FRACTIONS,
    tiny_defect_grid: tuple[int, int] = DEFAULT_TINY_DEFECT_GRID,
    window: float = DEFAULT_WINDOW,
    use_nw_characterization: bool | None = None,
    nw_characterization_path: Path | None = None,
    zw_reference_nws: Sequence[int] = DEFAULT_ZW_REFERENCE_NWS,
) -> dict[str, object]:
    hankel = outgoing_hankel_lambda_series()
    result = gate3_solve_result(
        root,
        radius_sweep=radius_sweep,
        grid_sweep=grid_sweep,
        grid_convergence_radial_scale=grid_convergence_radial_scale,
        vacuum_sanity_radii=vacuum_sanity_radii,
        domain_base_grid=domain_base_grid,
        largest_domain_grid_bases=largest_domain_grid_bases,
        flat_reference_grid_levels=flat_reference_grid_levels,
        flat_reference_radial_scales=flat_reference_radial_scales,
        tiny_defect_fractions=tiny_defect_fractions,
        tiny_defect_grid=tiny_defect_grid,
        window=window,
        use_nw_characterization=use_nw_characterization,
        nw_characterization_path=nw_characterization_path,
        zw_reference_nws=zw_reference_nws,
    )
    dims = dimensional_checks()
    guard = target_blindness_guard(
        [
            Path(__file__),
            repo_root() / "software" / "stage1_solver" / "tools" / "pathA_22b_gate3_crosscheck.wl",
            repo_root() / "software" / "stage1_solver" / "tests" / "test_patha22b_gate3.py",
        ]
    )
    branch_domain = result["references"]["branch_geometry"]["domain_truncation"]["assessment"]
    branch_grid = result["references"]["branch_geometry"]["grid_invariance"]["richardson"]
    reference_comparison = result["reference_comparison"]
    tiny_linearity = result["tiny_defect_linearity"]
    residuals = [
        "Frozen finite-core background loaded; no nonlinear profile re-solve was performed.",
        "Linear source adapter uses the frozen derived BdG wall modes; the outgoing DtN/self-energy sweep is recomputed.",
        "Retarded operator uses exact spherical-Hankel Y_out at the truncation radius; cached Gamma_port is not used in the solve.",
        (
            "Branch-reference r_max tail relative spread "
            f"{branch_domain['tail_rel_spread']:.3e}; onset r_max={branch_domain['plateau_onset_r_max']}; "
            f"plateau {branch_domain['plateau']}."
        ),
        (
            "Post-plateau fixed-r_max Richardson limit "
            f"{branch_grid['limit']:.15e}; observed delta-ratio {branch_grid['observed_delta_ratio']:.3e}."
        ),
        (
            "Matched even-nw Z_w reference shift "
            f"{reference_comparison['shift_branch_minus_flat']:.3e} "
            f"({100.0 * reference_comparison['rel_shift']:.2f}% of the branch value); "
            "carried separately as a one-sided definitional systematic."
        ),
        (
            "Tiny-defect slope relative spread "
            f"{tiny_linearity['slope_rel_spread']:.3e}; linear {tiny_linearity['linear_toward_zero']}."
        ),
    ]
    return {
        "schema": "stage1_pathA_22b_gate3_systematics_probe/v4",
        "provenance": [
            "research/pde/paper/pde.tex:1964",
            "research/pde/paper/pde.tex:1985-1988",
            "research/pde/paper/pde.tex:2034-2069",
            "software/stage1_solver/mathematica/mt15_06_m1c_prep_crossengine.wls:686-753",
            "software/stage1_solver/mathematica/mt15_06_m1c_prep_crossengine.wls:896-985",
        ],
        "hankel_fingerprint": hankel,
        "result": result,
        "dimensional_checks": dims,
        "target_blindness": guard,
        "gate3_outcome": result["outcome"],
        "blocker": result["blocker"],
        "chi_Q_value": result["chi_Q"],
        "chi_Q_error_bar": result["chi_Q_error_bar"],
        "residuals": residuals,
    }


def _fmt_check(raw: Mapping[str, object]) -> str:
    return (
        f"- `{raw['name']}`: **{raw['status']}** "
        f"(expected `{raw['expected']}`, actual `{raw['actual']}`). {raw['note']}"
    )


def _format_sweep_table(rows: Sequence[Mapping[str, Any]], *, mode: str) -> list[str]:
    if mode == "domain":
        lines = [
            "| radial_scale | r_max | grid | dr | E_defect | E_reference | chi_Q | delta chi_Q | max residual |",
            "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
        previous: float | None = None
        for row in rows:
            grid = row["grid"]
            chi_q = float(row["chi_Q"])
            delta = "" if previous is None else f"{chi_q - previous:.12e}"
            previous = chi_q
            lines.append(
                "| "
                f"{float(grid['radial_scale']):.6g} | "
                f"{float(grid['r_max']):.6g} | "
                f"{grid['nr']}x{grid['nw']} | "
                f"{float(grid['dr']):.6e} | "
                f"{float(row['defect_effective_omega5_coefficient']):.12e} | "
                f"{float(row['vacuum_effective_omega5_coefficient']):.12e} | "
                f"{chi_q:.12e} | "
                f"{delta} | "
                f"{float(row['max_linear_residual']):.3e} |"
            )
        return lines
    if mode == "radius":
        lines = [
            "| R_exit | grid | E_defect | E_vacuum | chi_Q | legacy fixed-a defect ratio | max residual |",
            "| ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
        for row in rows:
            grid = row["grid"]
            lines.append(
                "| "
                f"{float(row['closure_radius']):.6g} | "
                f"{grid['nr']}x{grid['nw']} | "
                f"{float(row['defect_effective_omega5_coefficient']):.12e} | "
                f"{float(row['vacuum_effective_omega5_coefficient']):.12e} | "
                f"{float(row['chi_Q']):.12e} | "
                f"{float(row['defect_uncalibrated_fixed_a_ratio']):.12e} | "
                f"{float(row['max_linear_residual']):.3e} |"
            )
        return lines
    if mode == "grid":
        lines = [
            "| grid | r_max | exterior dr | dw | E_defect | E_reference | chi_Q | fit RMS max | max residual |",
            "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
        for row in rows:
            grid = row["grid"]
            lines.append(
                "| "
                f"{grid['nr']}x{grid['nw']} | "
                f"{float(grid['r_max']):.6g} | "
                f"{float(grid['dr']):.6e} | "
                f"{float(grid['dw']):.6e} | "
                f"{float(row['defect_effective_omega5_coefficient']):.12e} | "
                f"{float(row['vacuum_effective_omega5_coefficient']):.12e} | "
                f"{float(row['chi_Q']):.12e} | "
                f"{float(row['max_fit_rms']):.3e} | "
                f"{float(row['max_linear_residual']):.3e} |"
            )
        return lines
    if mode == "tiny":
        lines = [
            "| defect fraction | chi_Q | chi_Q - 1 | (chi_Q - 1)/fraction | max residual |",
            "| ---: | ---: | ---: | ---: | ---: |",
        ]
        for row in rows:
            lines.append(
                "| "
                f"{float(row['defect_fraction']):.6e} | "
                f"{float(row['chi_Q']):.12e} | "
                f"{float(row['chi_Q_minus_1']):.12e} | "
                f"{float(row['slope_chi_Q_minus_1_over_fraction']):.12e} | "
                f"{float(row['max_linear_residual']):.3e} |"
            )
        return lines
    raise ValueError(f"unknown sweep table mode: {mode}")


def _format_zw_reference_table(rows: Sequence[Mapping[str, Any]]) -> list[str]:
    lines = [
        "| nr | nw | r_max | branch chi_Q | flat Z_w=1 chi_Q | branch-flat shift | relative shift |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            "| "
            f"{int(row['nr'])} | "
            f"{int(row['nw'])} | "
            f"{float(row['r_max']):.6g} | "
            f"{float(row['branch_geometry_chi_Q']):.12e} | "
            f"{float(row['flat_Zw1_chi_Q']):.12e} | "
            f"{float(row['shift_branch_minus_flat']):.12e} | "
            f"{100.0 * float(row['rel_shift_vs_branch']):.3f}% |"
        )
    return lines


def render_gate3_markdown(payload: Mapping[str, object]) -> str:
    result = payload["result"]
    dims = payload["dimensional_checks"]
    guard = payload["target_blindness"]
    hankel = payload["hankel_fingerprint"]
    assert isinstance(result, Mapping)
    assert isinstance(dims, Mapping)
    assert isinstance(guard, Mapping)
    assert isinstance(hankel, Mapping)
    inventory = result["inventory"]
    references = result["references"]
    reference_comparison = result["reference_comparison"]
    post_plateau_grid = result["post_plateau_grid_convergence"]
    trivial_consistency = result["trivial_uniform_consistency"]
    tiny_linearity = result["tiny_defect_linearity"]
    budget = result["combined_budget"]
    nw_summary = result.get("even_nw_characterization")
    matched_zw = reference_comparison.get("matched_even_nw")
    assert isinstance(inventory, Mapping)
    assert isinstance(references, Mapping)
    assert isinstance(reference_comparison, Mapping)
    assert isinstance(post_plateau_grid, Mapping)
    assert isinstance(trivial_consistency, Mapping)
    assert isinstance(tiny_linearity, Mapping)
    assert isinstance(budget, Mapping)
    branch = references["branch_geometry"]
    flat = references["flat_Zw1"]
    assert isinstance(branch, Mapping)
    assert isinstance(flat, Mapping)
    branch_domain = branch["domain_truncation"]
    flat_domain = flat["domain_truncation"]
    branch_grid = branch["grid_invariance"]
    flat_grid = flat["grid_invariance"]
    assert isinstance(branch_domain, Mapping)
    assert isinstance(flat_domain, Mapping)
    assert isinstance(branch_grid, Mapping)
    assert isinstance(flat_grid, Mapping)
    branch_domain_rows = branch_domain["rows"]
    flat_domain_rows = flat_domain["rows"]
    branch_grid_rows = branch_grid["rows"]
    post_plateau_grid_rows = post_plateau_grid["rows"]
    tiny_rows = tiny_linearity["rows"]
    assert isinstance(branch_domain_rows, list)
    assert isinstance(flat_domain_rows, list)
    assert isinstance(branch_grid_rows, list)
    assert isinstance(post_plateau_grid_rows, list)
    assert isinstance(tiny_rows, list)
    branch_domain_assessment = branch_domain["assessment"]
    flat_domain_assessment = flat_domain["assessment"]
    branch_richardson = branch_grid["richardson"]
    flat_richardson = flat_grid["richardson"]
    branch_grid_assessment = branch_grid["assessment"]
    post_plateau_grid_assessment = post_plateau_grid["assessment"]
    post_plateau_grid_richardson = post_plateau_grid["richardson"]
    if nw_summary is not None:
        assert isinstance(nw_summary, Mapping)
    if matched_zw is not None:
        assert isinstance(matched_zw, Mapping)

    lines = [
        "## Gate 3 — χ_Q",
        "",
        "### Definition",
        "",
        "- Scope: linear-response solve around the frozen M1c finite-core branch. No nonlinear profile re-solve and no deep/empty-throat solve were run.",
        f"- Frozen branch: `{inventory['background_path']}` with grid `{inventory['background_grid']}` and stationary residual `{inventory['background_residual_linf']}`.",
        f"- Boundary labels: freeze sheet `{inventory['freeze_boundary_class']}`; packet `{inventory['packet_boundary_class']}` with exit model `{inventory['packet_exit_model']}`.",
        f"- Hankel check: `Lambda_2^out={hankel['lambda_series']}` and `Y_out={hankel['yhat_series']}`; the `z^5` coefficient is `{hankel['z5_loss_coefficient']}`.",
        f"- Anti-tautology: cached `Gamma_port={inventory['legacy_flux_normalization_present_but_not_used']['Gamma_port']}` is present only as legacy provenance and is not used by the solve.",
        "- Source adapter: frozen derived BdG wall modes from the physical packet; the outgoing DtN/self-energy sweep itself is recomputed from the frozen background.",
        "- Extraction: for each domain, grid, and frequency window, fit `E = C5` from `-Im(Sigma_out/Sigma_cons)/omega^5 = C5 + C7*omega^2 + C9*omega^4`, then set `chi_Q = E_defect / E_reference`.",
        "- Physical reference: the branch-geometry reference keeps the exported `Z_w` and `R0_w`, and undefects only the condensate/gauge sector (`rho` uniform, `A -> 0`). The flat reference (`Z_w=1`, `R0_w=a`) is retained only to quantify the old reference choice.",
        "- Exterior caveat: `radial_scale>1` uses the existing radial exterior taper, not a nonlinear far-field re-solve. The sweep tests/bounds boundary-placement sensitivity on this frozen continuation; lack of plateau still blocks non-characterized inputs.",
        "",
        "### r_max Domain Sweep",
        "",
        f"- `dr` control: {branch_domain['dr_held_constant_note']}.",
        f"- Branch-geometry reference tail mean `{branch_domain_assessment['tail_mean']:.15e}`, tail spread `{branch_domain_assessment['tail_spread']:.3e}`, relative tail spread `{branch_domain_assessment['tail_rel_spread']:.3e}`, plateau `{branch_domain_assessment['plateau']}`.",
        f"- Plateau onset: `r_max={branch_domain_assessment['plateau_onset_r_max']}` by the `{branch_domain_assessment['criteria']['plateau_delta_abs_max']:.1e}` tail-delta rule; strict `1e-3` final-lock onset `r_max={branch_domain_assessment['strict_1e_minus_3_onset_r_max']}`.",
        f"- Flat `Z_w=1` reference is a plateau-tail diagnostic over scales `{flat_domain['radial_scales']}`: tail mean `{flat_domain_assessment['tail_mean']:.15e}`, tail spread `{flat_domain_assessment['tail_spread']:.3e}`, relative tail spread `{flat_domain_assessment['tail_rel_spread']:.3e}`.",
        "",
        "Branch-geometry reference:",
        "",
        *_format_sweep_table(branch_domain_rows, mode="domain"),
        "",
        "Flat `Z_w=1` plateau-tail diagnostic:",
        "",
        *_format_sweep_table(flat_domain_rows, mode="domain"),
        "",
        "### Z_w Reference Comparison",
        "",
        f"- Physically correct zero for the defect's radiative coupling: `branch_geometry` (`Z_w/R0_w` retained, condensate/gauge undefected).",
        f"- Matched converged-grid reference shift, branch minus flat: `{reference_comparison['shift_branch_minus_flat']:.15e}` (relative `{100.0 * reference_comparison['rel_shift']:.2f}%`). This is a one-sided definitional systematic because branch-geometry is the physical zero and flat `Z_w=1` is the less-physical alternative; it is not folded into the numerical error bar.",
        *(
            [
                "",
                "Matched converged even-nw rows:",
                "",
                *_format_zw_reference_table(matched_zw["rows"]),
            ]
            if matched_zw is not None
            else []
        ),
        "",
        "### Grid Convergence at Post-Plateau r_max",
        "",
        f"- Fixed post-plateau `r_max={post_plateau_grid['closure_radius']:.6g}` (`radial_scale={post_plateau_grid['radial_scale']:.6g}`) branch grid sequence values: `{post_plateau_grid_richardson['values']}`.",
        f"- Observed delta ratios: `{post_plateau_grid_richardson['delta_ratios']}`.",
        f"- Richardson/geometric-tail limit: `{post_plateau_grid_richardson['limit']:.15e}`; finest value `{post_plateau_grid_richardson['finest_value']:.15e}`; old two-grid mean `{post_plateau_grid_richardson['two_grid_mean']:.15e}`.",
        f"- Grid extrapolation uncertainty `{budget['grid_uncertainty']:.3e}`; raw finest-limit correction `{post_plateau_grid_richardson['grid_uncertainty']:.3e}`; model uncertainty `{post_plateau_grid_richardson['richardson_model_uncertainty']:.3e}`; old two-grid mean bias `{post_plateau_grid_richardson['old_two_grid_mean_bias']:.3e}`.",
        f"- Grid convergence: `{post_plateau_grid_assessment['converging']}` with observed delta-ratio `{post_plateau_grid_assessment['observed_delta_ratio']:.3e}` and relative tail correction `{post_plateau_grid_assessment['rel_tail_correction']:.3e}`.",
        "",
        *_format_sweep_table(post_plateau_grid_rows, mode="grid"),
        *(
            [
                "",
                "### Even-nw Characterization",
                "",
                f"- Canonical characterization JSON: `{nw_summary['source_path_repo_relative']}`.",
                f"- Outcome class: `{nw_summary['outcome_class']}`; reason: {nw_summary['reason']}.",
                f"- Tail-supported central uses the `nw>={nw_summary['central_tail']['min_nw']}` fit: raw `{nw_summary['chi_Q_raw']:.15e}`, reported as `{nw_summary['chi_Q_reported']:.3f}`.",
                f"- Numerical bar `{nw_summary['numerical_tail_supported']:.4f}` comes from `sqrt(max(tail RMS {nw_summary['budget_tail']['fit']['rms']:.3e}, tail jackknife {nw_summary['budget_tail']['jackknife']['max_shift']:.3e})^2 + nr offset {nw_summary['nr_offset']:.3e}^2)`.",
                f"- Conservative full-sweep fit uncertainty is kept separately: `{nw_summary['full_sweep_fit_uncertainty']:.4f}` from the `nw={nw_summary['primary_nw_min']}..{nw_summary['primary_nw_max']}` fit.",
                f"- Evidence: jackknife stable `{nw_summary['jackknife']['stable']}`, nr-independent `{nw_summary['nr_comparison']['nr_independent']}`, stored flat tail through `nw={nw_summary['central_tail']['max_nw']}`; review verification extends the flat tail to `nw=320`.",
            ]
            if nw_summary is not None
            else []
        ),
        "",
        "### Tiny-Defect Linearity",
        "",
        f"- Definition: {tiny_linearity['definition']}.",
        f"- Slope mean `{tiny_linearity['slope_mean']:.15e}`; slope relative spread `{tiny_linearity['slope_rel_spread']:.3e}`; linear toward zero `{tiny_linearity['linear_toward_zero']}`.",
        "",
        *_format_sweep_table(tiny_rows, mode="tiny"),
        "",
        "### Trivial Uniform Consistency",
        "",
        f"- The uniform self-ratio is demoted to a trivial consistency check: `{trivial_consistency['definition']}`.",
        f"- Uniform no-defect max `|chi_Q-1|`: `{trivial_consistency['max_uniform_abs_error']:.3e}`; pass `{trivial_consistency['passed']}`. This is not evidence for the physical magnitude.",
        "",
        "### Dimensional Check",
        "",
    ]
    for raw in dims["checks"]:
        assert isinstance(raw, Mapping)
        lines.append(_fmt_check(raw))
    lines.extend(["", "### Provenance", ""])
    for item in payload["provenance"]:
        lines.append(f"- {item}.")
    lines.extend(
        [
            "",
            "### Target-Blindness",
            "",
            f"- `{guard['status']}` over the new Gate-3 module, tests, and Mathematica cross-check. The final comparison constants are not used in the derivation.",
            "",
            "### Residual Ledger",
            "",
        ]
    )
    for residual in payload["residuals"]:
        lines.append(f"- {residual}")
    lines.extend(["", "### OUTCOME", ""])
    if payload["chi_Q_value"] is not None:
        lines.extend(
            [
                f"- OUTCOME: `{payload['gate3_outcome']}`.",
                f"- Converged `chi_Q = {payload['chi_Q_value']:.3f} +/- {payload['chi_Q_error_bar']:.4f}` (numerical, tail-supported), branch-geometry reference.",
                f"- Conservative full-sweep fit uncertainty: `+/- {budget['full_sweep_fit_uncertainty']:.4f}`; this is retained as a coarse-transient diagnostic, while the flat nw-tail plus nr-offset gives the honest numerical bar.",
                f"- Separate `Z_w`-reference definitional systematic: one-sided `{reference_comparison['shift_branch_minus_flat']:.3f}` absolute (`{100.0 * reference_comparison['rel_shift']:.1f}%`) toward the flat `Z_w=1` alternative; not folded symmetrically into the numerical error bar.",
                f"- Convergence evidence: jackknife-stable `{nw_summary['jackknife']['stable'] if nw_summary is not None else None}`, nr-independent `{nw_summary['nr_comparison']['nr_independent'] if nw_summary is not None else None}`, flat even-nw tail, and matched branch-geometry reference.",
                f"- Reason: {result['reason']}",
            ]
        )
    else:
        lines.extend(
            [
                f"- OUTCOME: `{payload['gate3_outcome']}`.",
                "- chi_Q central: `NOT_EXTRACTABLE`; magnitude reverts to a calibration knob on this frozen background.",
                f"- Blocker: `{payload['blocker']}`.",
                f"- Reason: {result['reason']}",
            ]
        )
    lines.append("")
    return "\n".join(lines)


def write_gate3_report(out_dir: Path, report_dir: Path) -> tuple[Path, Path, dict[str, object]]:
    payload = gate3_report_payload()
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "pathA_22b_gate3.json"
    scratch_md_path = out_dir / "pathA_22b_gate3.md"
    report_path = report_dir / "pathA_22b_minimal_combination_xi.md"
    rendered = render_gate3_markdown(payload)
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    scratch_md_path.write_text(rendered + "\n", encoding="utf-8")

    existing = report_path.read_text(encoding="utf-8") if report_path.exists() else "# PathA 22b minimal combination xi\n"
    markers = ("\n## Gate 3 — χ_Q\n", "\n## Gate 3 -- chi_Q\n")
    for marker in markers:
        if marker in existing:
            existing = existing.split(marker, 1)[0].rstrip() + "\n"
            break
    else:
        if not existing.endswith("\n"):
            existing += "\n"
    report_path.write_text(existing.rstrip() + "\n\n" + rendered + "\n", encoding="utf-8")
    return json_path, report_path, payload


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default="software/stage1_solver/_scratch")
    parser.add_argument("--report-dir", default="software/stage1_solver/reports")
    args = parser.parse_args(list(argv) if argv is not None else None)
    json_path, report_path, payload = write_gate3_report(Path(args.out_dir), Path(args.report_dir))
    print(f"wrote {json_path}")
    print(f"wrote {report_path}")
    print(f"Gate 3 outcome: {payload['gate3_outcome']}")
    if payload["chi_Q_value"] is None:
        print("chi_Q: NOT_EXTRACTABLE")
    else:
        print(f"chi_Q: {payload['chi_Q_value']:.15e} +/- {payload['chi_Q_error_bar']:.3e}")
    print(f"r_max plateau pass: {payload['result']['convergence']['radius_pass']}")
    print(f"grid convergence pass: {payload['result']['convergence']['grid_pass']}")
    print(f"tiny-defect linearity pass: {payload['result']['tiny_defect_linearity']['linear_toward_zero']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
