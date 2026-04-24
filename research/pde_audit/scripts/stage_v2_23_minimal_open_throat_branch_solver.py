#!/usr/bin/env python3
"""
Stage V2-23: Minimal open-throat branch solver and first real residual extraction.

This is a reduced branch-realization prototype, not a full nonlinear GNLS/Maxwell
moving-throat PDE solver. It solves a frozen one-dimensional open-throat geometry
and Sturm-Liouville support/gauge profiles, computes actual overlap integrals,
builds the V2-21/V2-22 observable packet, and reports residuals against the
isotropic full-bundle target surface.

The branch definition is target-blind: the geometry, potentials, and coupling
scales are declared before residuals are evaluated.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import math
from typing import Any, Dict, List, Tuple

import numpy as np

try:
    from scipy.linalg import eigh
except Exception as exc:  # pragma: no cover
    eigh = None
    _SCIPY_IMPORT_ERROR = repr(exc)
else:
    _SCIPY_IMPORT_ERROR = None


SCRIPT_DIR = Path(__file__).resolve().parent
ARTIFACT_DIR = SCRIPT_DIR / "output" / "artifacts"


def trapz(y: np.ndarray, x: np.ndarray) -> float:
    trapezoid = getattr(np, "trapezoid", np.trapz)
    return float(trapezoid(y, x))


def central_gradient(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    return np.gradient(y, x, edge_order=2)


def sha256_json(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), default=float).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def assert_finite_array(name: str, arr: np.ndarray) -> None:
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")


def assemble_fem_matrices(
    s: np.ndarray,
    mu: np.ndarray,
    T: np.ndarray,
    V: np.ndarray,
    robin_Y_right: float = 0.0,
    dirichlet_left: bool = True,
) -> Tuple[np.ndarray, np.ndarray, List[int]]:
    """Assemble 1D linear-FEM stiffness and mass matrices.

    Weak operator:
        ∫ T q' v' ds + ∫ V q v ds + Y_R q(L)v(L) = λ ∫ μ q v ds

    The left boundary is Dirichlet by default. The right boundary is open
    impedance/Robin; Y=0 gives the Neumann/free-end limit.
    """
    n = len(s)
    K = np.zeros((n, n), dtype=float)
    M = np.zeros((n, n), dtype=float)

    for e in range(n - 1):
        h = s[e + 1] - s[e]
        Tm = 0.5 * (T[e] + T[e + 1])
        Vm = 0.5 * (V[e] + V[e + 1])
        mum = 0.5 * (mu[e] + mu[e + 1])

        ke_grad = (Tm / h) * np.array([[1.0, -1.0], [-1.0, 1.0]])
        ke_pot = (Vm * h / 6.0) * np.array([[2.0, 1.0], [1.0, 2.0]])
        me = (mum * h / 6.0) * np.array([[2.0, 1.0], [1.0, 2.0]])

        idx = [e, e + 1]
        for a in range(2):
            for b in range(2):
                K[idx[a], idx[b]] += ke_grad[a, b] + ke_pot[a, b]
                M[idx[a], idx[b]] += me[a, b]

    if robin_Y_right != 0.0:
        K[-1, -1] += float(robin_Y_right)

    if dirichlet_left:
        keep = list(range(1, n))
    else:
        keep = list(range(n))

    return K[np.ix_(keep, keep)], M[np.ix_(keep, keep)], keep


def solve_sturm_liouville(
    s: np.ndarray,
    mu: np.ndarray,
    T: np.ndarray,
    V: np.ndarray,
    robin_Y_right: float,
    mode_index: int = 0,
    label: str = "mode",
) -> Dict[str, Any]:
    """Solve a generalized Sturm-Liouville eigenproblem and return one mode."""
    if eigh is None:
        raise RuntimeError(f"scipy.linalg.eigh unavailable: {_SCIPY_IMPORT_ERROR}")

    Kmat, Mmat, keep = assemble_fem_matrices(
        s=s, mu=mu, T=T, V=V, robin_Y_right=robin_Y_right, dirichlet_left=True
    )
    vals, vecs = eigh(Kmat, Mmat)
    vals = np.real(vals)
    positive = np.where(vals > 1e-12)[0]
    if len(positive) <= mode_index:
        raise RuntimeError(f"{label}: not enough positive eigenvalues")

    j = int(positive[mode_index])
    lam = float(vals[j])
    v_red = np.real(vecs[:, j])
    full = np.zeros(len(s), dtype=float)
    full[keep] = v_red

    # Fix an orientation convention: positive bulk average.
    if trapz(mu * full, s) < 0:
        full *= -1.0

    norm = math.sqrt(trapz(mu * full * full, s))
    if norm <= 0:
        raise RuntimeError(f"{label}: zero norm eigenvector")
    full /= norm

    # Residual using FEM matrices and the normalized reduced vector.
    red = full[keep]
    residual = Kmat @ red - lam * (Mmat @ red)
    residual_rel = float(np.linalg.norm(residual) / max(np.linalg.norm(Kmat @ red), 1e-30))

    return {
        "label": label,
        "eigenvalue": lam,
        "frequency": math.sqrt(lam),
        "profile": full,
        "norm_mu": trapz(mu * full * full, s),
        "left_value": float(full[0]),
        "right_derivative": float((full[-1] - full[-2]) / (s[-1] - s[-2])),
        "right_value": float(full[-1]),
        "fem_residual_relative": residual_rel,
    }


def solve_open_shape_profile(
    s: np.ndarray,
    a_mouth: float,
    R_exit_pref: float,
    T_R: float,
    K_R: float,
    Y_exit: float,
) -> Dict[str, Any]:
    """Solve a target-blind stationary open-throat shape profile.

    Energy:
        E[R] = 1/2 ∫ T_R R'^2 ds
             + 1/2 ∫ K_R (R - R_pref)^2 ds
             + 1/2 Y_exit (R(L) - R_exit_pref)^2

    with fixed mouth R(0)=a_mouth and open finite-radius exit.
    """
    n = len(s)
    L = float(s[-1] - s[0])
    x = s / L
    # Smooth preferred open conduit: finite exit radius and a small bulge.
    R_pref = a_mouth - (a_mouth - R_exit_pref) * x + 0.04 * np.sin(np.pi * x)

    K = np.zeros((n, n), dtype=float)
    b = np.zeros(n, dtype=float)

    for e in range(n - 1):
        h = s[e + 1] - s[e]
        ke_grad = (T_R / h) * np.array([[1.0, -1.0], [-1.0, 1.0]])
        ke_mass = (K_R * h / 6.0) * np.array([[2.0, 1.0], [1.0, 2.0]])
        # Consistent RHS for K_R R_pref. Use nodal interpolation.
        rp = np.array([R_pref[e], R_pref[e + 1]])
        be = ke_mass @ rp
        idx = [e, e + 1]
        for aa in range(2):
            b[idx[aa]] += be[aa]
            for bb in range(2):
                K[idx[aa], idx[bb]] += ke_grad[aa, bb] + ke_mass[aa, bb]

    # Exit penalty / open finite-radius preference.
    K[-1, -1] += Y_exit
    b[-1] += Y_exit * R_exit_pref

    # Apply Dirichlet mouth R(0)=a_mouth.
    keep = list(range(1, n))
    fixed = [0]
    Kuu = K[np.ix_(keep, keep)]
    Kuf = K[np.ix_(keep, fixed)]
    bu = b[keep] - Kuf[:, 0] * a_mouth
    R_unknown = np.linalg.solve(Kuu, bu)
    R = np.zeros(n, dtype=float)
    R[0] = a_mouth
    R[keep] = R_unknown

    Rp = central_gradient(R, s)
    residual = K @ R - b
    residual[0] = 0.0  # Dirichlet row is not part of the free stationarity.
    residual_rel = float(np.linalg.norm(residual[1:]) / max(np.linalg.norm(b[1:]), 1e-30))

    return {
        "R": R,
        "R_pref": R_pref,
        "R_prime": Rp,
        "R_exit": float(R[-1]),
        "R_min": float(np.min(R)),
        "R_mouth": float(R[0]),
        "stationary_residual_relative": residual_rel,
    }


def compute_coefficients(branch: Dict[str, Any]) -> Dict[str, Any]:
    s = np.array(branch["grid"]["s"], dtype=float)
    mu = np.array(branch["measure"]["mu_s"], dtype=float)

    chi = np.array(branch["profiles"]["chi_eta"], dtype=float)
    phi = np.array(branch["profiles"]["phi_bdg"], dtype=float)
    u = np.array(branch["profiles"]["u_gauge"], dtype=float)
    w = np.array(branch["profiles"]["w_mixed"], dtype=float)

    kappa = branch["couplings"]
    eig = branch["eigenvalues"]

    I_eta_phi = trapz(mu * chi * phi, s)
    I_eta_u = trapz(mu * chi * u, s)
    I_eta_w = trapz(mu * chi * w, s)
    I_uw = trapz(mu * u * w, s)

    cB = kappa["lambda_B"] * I_eta_phi
    gU = kappa["lambda_U"] * I_eta_u
    gW = kappa["lambda_W"] * I_eta_w
    Rmix = kappa["lambda_R"] * I_uw

    K = float(eig["lambda_wall_l2"])
    M = 1.0
    varpi2 = float(eig["lambda_bdg"])
    OU2 = float(eig["lambda_U"])
    OW2 = float(eig["lambda_W"])

    B0 = cB**2 / varpi2
    B2 = cB**2 / varpi2**2
    B4 = cB**2 / varpi2**3

    Delta = OU2 * OW2 - Rmix**2
    S = OU2 + OW2
    Q = gU**2 * OW2 + 2 * gU * gW * Rmix + gW**2 * OU2
    H = gU**2 + gW**2

    Z0 = Q / Delta
    Z2 = (Q * S - H * Delta) / Delta**2
    Z4 = (Q * (S**2 - Delta) - S * H * Delta) / Delta**3

    Pproto = OU2 * gW + Rmix * gU
    N0 = Pproto**2 / Delta**2
    N2 = 2 * Pproto * (Pproto * S - Delta * gW) / Delta**3
    N4 = (
        Delta**2 * gW**2
        - 2 * Delta * Pproto**2
        - 4 * Delta * Pproto * S * gW
        + 3 * Pproto**2 * S**2
    ) / Delta**4

    D0 = K - B0 - Z0
    D2 = -(M + B2 + Z2)
    D4 = -(B4 + Z4)

    u2 = -D2 / D0
    u4 = (D2**2 - D0 * D4) / D0**2

    P0 = N0 / D0
    P2 = (D0 * N2 - 2 * D2 * N0) / D0**2
    P4 = (D0**2 * N4 - 2 * D0 * (D2 * N2 + D4 * N0) + 3 * D2**2 * N0) / D0**3

    constants = branch["constants"]
    P0_target = 54 * constants["G"] * constants["c_s"]**5 / (
        5 * constants["a_throat"]**5 * constants["c_light"]**5
    )
    R_pole = D0 * (B4 + Z4) - 3 * (M + B2 + Z2) ** 2
    R_norm = constants["mhat0"] ** 2 * constants["S_port"] * P0 - P0_target
    R_tail = constants["Theta_tail"] * (constants["c_light"] / constants["c_s"]) ** 3 - 1.0

    gamma_eff = constants["mhat0"] ** 2 * constants["S_port"] * P0 * constants["a_throat"]**5 / (
        27 * constants["c_s"]**5
    )
    gamma_GR = 2 * constants["G"] / (5 * constants["c_light"]**5)

    coeffs = {
        "overlaps": {
            "I_eta_phi": I_eta_phi,
            "I_eta_u": I_eta_u,
            "I_eta_w": I_eta_w,
            "I_uw": I_uw,
        },
        "reduced_couplings": {
            "c_B": cB,
            "g_U": gU,
            "g_W": gW,
            "R_mix": Rmix,
        },
        "wall": {"K": K, "M": M},
        "bdg": {"varpi2": varpi2, "B0": B0, "B2": B2, "B4": B4},
        "mixed": {
            "Omega_U2": OU2,
            "Omega_W2": OW2,
            "Delta": Delta,
            "S": S,
            "Q": Q,
            "H": H,
            "Pproto": Pproto,
            "Z0": Z0,
            "Z2": Z2,
            "Z4": Z4,
            "N0": N0,
            "N2": N2,
            "N4": N4,
        },
        "operator": {"D0": D0, "D2": D2, "D4": D4},
        "response": {"u2": u2, "u4": u4},
        "prefactor": {"P0": P0, "P2": P2, "P4": P4},
        "targets": {
            "P0_target": P0_target,
            "gamma_eff": gamma_eff,
            "gamma_GR": gamma_GR,
            "R_pole": R_pole,
            "R_norm": R_norm,
            "R_P2": P2,
            "R_P4": P4,
            "R_tail": R_tail,
        },
    }
    return coeffs


def build_branch(grid_points: int = 181) -> Dict[str, Any]:
    # Frozen, target-blind branch packet.
    N = int(grid_points)
    if N < 25:
        raise ValueError("grid_points must be at least 25 for the reduced FEM branch")

    constants = {
        "G": 1.0,
        "c_s": 1.0,
        "c_light": 1.0,
        "a_throat": 1.0,
        "mhat0": 1.0,
        "S_port": 1.0,
        "Theta_tail": 1.0,
    }
    geometry_params = {
        "L": 1.85,
        "a_mouth": 1.0,
        "R_exit_pref": 0.42,
        "T_R": 0.35,
        "K_R": 6.0,
        "Y_exit": 2.5,
        "boundary_class": "open_impedance",
        "boundary_protocol": "open_impedance_AC_reflecting_DC_leaking",
        "AC_robin_Y": 0.0,
    }
    couplings = {
        "lambda_B": 0.60,
        "lambda_U": 0.55,
        "lambda_W": 0.75,
        "lambda_R": 0.25,
    }
    L = geometry_params["L"]
    s = np.linspace(0.0, L, N)

    shape = solve_open_shape_profile(
        s=s,
        a_mouth=geometry_params["a_mouth"],
        R_exit_pref=geometry_params["R_exit_pref"],
        T_R=geometry_params["T_R"],
        K_R=geometry_params["K_R"],
        Y_exit=geometry_params["Y_exit"],
    )
    R0 = shape["R"]
    Rp = shape["R_prime"]

    # Effective axial measure for an S^2 throat surface, normalized to mean 1.
    mu_raw = R0**2 * np.sqrt(1.0 + Rp**2)
    mu_s = mu_raw * (L / trapz(mu_raw, s))
    assert_finite_array("mu_s", mu_s)

    x = s / L

    # Densitized wall coefficients for the l=2 support lane.
    T_wall = 1.0 + 0.12 * x
    K_eta = 0.32 + 0.08 * x**2
    T_Omega = 0.045 / (R0**2 + 0.08)
    V_wall_l2 = K_eta + 6.0 * T_Omega

    # Stable support and internal Maxwell/mixed operators.
    T_bdg = 0.92 + 0.06 * (1.0 - x)
    V_bdg = 1.65 + 0.18 * np.cos(0.5 * np.pi * x) ** 2

    T_U = 1.08 + 0.03 * x
    V_U = 2.05 + 0.10 * x

    T_W = 0.86 + 0.05 * (1.0 - x)
    V_W = 2.45 + 0.12 * np.sin(0.5 * np.pi * x) ** 2

    robin = geometry_params["AC_robin_Y"]

    wall_mode = solve_sturm_liouville(s, mu_s, T_wall, V_wall_l2, robin, 0, "wall_l2_open")
    bdg_mode = solve_sturm_liouville(s, mu_s, T_bdg, V_bdg, robin, 0, "bdg_support_open")
    U_mode = solve_sturm_liouville(s, mu_s, T_U, V_U, robin, 0, "brane_like_U_open")
    W_mode = solve_sturm_liouville(s, mu_s, T_W, V_W, robin, 0, "mixed_W_open")

    branch: Dict[str, Any] = {
        "schema": "stage_v2_23_minimal_open_throat_branch_solver/v1",
        "branch_id": "v2_23_minimal_open_throat_frozen_demo",
        "status": "reduced_branch_solver_prototype_not_full_nonlinear_PDE",
        "freeze": {
            "pre_target_freeze": True,
            "target_blind": True,
            "no_post_residual_refit": True,
            "parent_action_reading": "GNLS + localized Maxwell + promoted effective wall action",
            "gauge_convention": "localized/controlled gauge bookkeeping from V2-02; mixed fields retained",
            "geometry_protocol": "finite open throat with impedance AC reflection and DC leakage",
        },
        "constants": constants,
        "geometry": {**geometry_params, **{k: shape[k] for k in ["R_exit", "R_min", "R_mouth", "stationary_residual_relative"]}},
        "grid": {"N": N, "L": L, "s": s.tolist()},
        "measure": {"mu_s": mu_s.tolist(), "mu_integral": trapz(mu_s, s)},
        "coefficients_1d": {
            "T_wall": T_wall.tolist(),
            "K_eta": K_eta.tolist(),
            "T_Omega": T_Omega.tolist(),
            "V_wall_l2": V_wall_l2.tolist(),
            "T_bdg": T_bdg.tolist(),
            "V_bdg": V_bdg.tolist(),
            "T_U": T_U.tolist(),
            "V_U": V_U.tolist(),
            "T_W": T_W.tolist(),
            "V_W": V_W.tolist(),
        },
        "profiles": {
            "R0": R0.tolist(),
            "R_prime": Rp.tolist(),
            "chi_eta": wall_mode["profile"].tolist(),
            "phi_bdg": bdg_mode["profile"].tolist(),
            "u_gauge": U_mode["profile"].tolist(),
            "w_mixed": W_mode["profile"].tolist(),
        },
        "eigenvalues": {
            "lambda_wall_l2": wall_mode["eigenvalue"],
            "lambda_bdg": bdg_mode["eigenvalue"],
            "lambda_U": U_mode["eigenvalue"],
            "lambda_W": W_mode["eigenvalue"],
        },
        "mode_diagnostics": {
            "wall": {k: v for k, v in wall_mode.items() if k != "profile"},
            "bdg": {k: v for k, v in bdg_mode.items() if k != "profile"},
            "U": {k: v for k, v in U_mode.items() if k != "profile"},
            "W": {k: v for k, v in W_mode.items() if k != "profile"},
        },
        "couplings": couplings,
    }
    branch["branch_freeze_hash"] = sha256_json({
        "schema": branch["schema"],
        "branch_id": branch["branch_id"],
        "freeze": branch["freeze"],
        "constants": constants,
        "geometry_params": geometry_params,
        "couplings": couplings,
        "coefficient_family": "linear_FEM_open_throat_reduced_branch_v1",
    })
    return branch


def evaluate_gates(branch: Dict[str, Any], coeffs: Dict[str, Any]) -> Dict[str, Any]:
    D0 = coeffs["operator"]["D0"]
    C = coeffs["bdg"]["B4"] + coeffs["mixed"]["Z4"]
    N0 = coeffs["mixed"]["N0"]
    Delta = coeffs["mixed"]["Delta"]
    R_exit = branch["geometry"]["R_exit"]

    tolerances = {
        "geometry_exit_min": 1e-8,
        "stability_abs_min": 1e-10,
        "target_rel_tol": 1e-6,
        "formula_abs_tol": 1e-8,
    }

    residuals = coeffs["targets"]
    P0_target = residuals["P0_target"]
    target_rel = {
        "R_pole_rel": abs(residuals["R_pole"]) / max(abs(D0 * C), 1.0),
        "R_norm_rel": abs(residuals["R_norm"]) / max(abs(P0_target), 1.0),
        "R_P2_abs": abs(residuals["R_P2"]),
        "R_P4_abs": abs(residuals["R_P4"]),
        "R_tail_abs": abs(residuals["R_tail"]),
    }

    open_gate = (branch["geometry"]["boundary_class"] == "open_impedance") and (R_exit > tolerances["geometry_exit_min"])
    stability_gate = (D0 > tolerances["stability_abs_min"]) and (C > 0.0) and (Delta > 0.0)
    transfer_gate = abs(N0) > tolerances["stability_abs_min"]

    target_gate = (
        target_rel["R_pole_rel"] < tolerances["target_rel_tol"]
        and target_rel["R_norm_rel"] < tolerances["target_rel_tol"]
        and target_rel["R_P2_abs"] < tolerances["target_rel_tol"]
        and target_rel["R_P4_abs"] < tolerances["target_rel_tol"]
        and target_rel["R_tail_abs"] < tolerances["target_rel_tol"]
    )

    return {
        "tolerances": tolerances,
        "open_gate_pass": bool(open_gate),
        "stability_gate_pass": bool(stability_gate),
        "outgoing_transfer_gate_pass": bool(transfer_gate),
        "target_packet_pass": bool(target_gate),
        "relative_residuals": target_rel,
        "failure_diagnosis": {
            "dominant_residual": max(target_rel.items(), key=lambda kv: kv[1])[0],
            "normalization_ratio_P0_over_target": coeffs["prefactor"]["P0"] / P0_target,
            "one_pole_ratio_D0C_over_3A2": (
                (coeffs["operator"]["D0"] * (coeffs["bdg"]["B4"] + coeffs["mixed"]["Z4"]))
                / max(3.0 * (coeffs["wall"]["M"] + coeffs["bdg"]["B2"] + coeffs["mixed"]["Z2"]) ** 2, 1e-30)
            ),
        },
    }


def sympy_formula_checks() -> Dict[str, Any]:
    import sympy as sp

    w = sp.symbols("w")
    D0, D2, D4, N0, N2, N4 = sp.symbols("D0 D2 D4 N0 N2 N4", nonzero=True)
    a, cs, G, c, mhat, Sport = sp.symbols("a cs G c mhat Sport", positive=True)

    D = D0 + D2*w**2 + D4*w**4
    N = N0 + N2*w**2 + N4*w**4

    Y = sp.series(D0/D, w, 0, 6).removeO().expand()
    u2 = -D2/D0
    u4 = (D2**2 - D0*D4)/D0**2

    Pref = sp.series(D0*N/D**2, w, 0, 6).removeO().expand()
    P0 = N0/D0
    P2 = (D0*N2 - 2*D2*N0)/D0**2
    P4 = (D0**2*N4 - 2*D0*(D2*N2 + D4*N0) + 3*D2**2*N0)/D0**3

    one_pole = sp.simplify((u4 - 4*u2**2) * D0**2 - (-3*D2**2 - D0*D4))
    N2_const = 2*D2*N0/D0
    N4_const = N0*(D2**2 + 2*D0*D4)/D0**2
    P2_const = sp.simplify(P2.subs(N2, N2_const))
    P4_const = sp.simplify(P4.subs({N2: N2_const, N4: N4_const}))

    P0target = 54*G*cs**5/(5*a**5*c**5)
    gamma_eff = Sport*mhat**2*P0target*a**5/(27*cs**5)
    gamma_GR = 2*G/(5*c**5)

    checks = {
        "response_u2": sp.simplify(Y.coeff(w, 2) - u2) == 0,
        "response_u4": sp.simplify(Y.coeff(w, 4) - u4) == 0,
        "prefactor_P0": sp.simplify(Pref.coeff(w, 0) - P0) == 0,
        "prefactor_P2": sp.simplify(Pref.coeff(w, 2) - P2) == 0,
        "prefactor_P4": sp.simplify(Pref.coeff(w, 4) - P4) == 0,
        "one_pole_residual_identity": one_pole == 0,
        "constant_prefactor_P2": P2_const == 0,
        "constant_prefactor_P4": P4_const == 0,
        "normalization_equivalence": sp.simplify(gamma_eff.subs({Sport:1, mhat:1}) - gamma_GR) == 0,
    }
    return {
        "checks": checks,
        "passed": int(sum(bool(v) for v in checks.values())),
        "total": len(checks),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage V2-23 minimal open-throat branch solver")
    parser.add_argument("--grid-points", type=int, default=181, help="Number of axial grid points for the reduced FEM solve")
    args = parser.parse_args()

    ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)

    branch = build_branch(grid_points=args.grid_points)
    coeffs = compute_coefficients(branch)
    gates = evaluate_gates(branch, coeffs)
    formula = sympy_formula_checks()

    observable_packet = {
        "schema": "stage_v2_23_observable_packet/v1",
        "source_branch_id": branch["branch_id"],
        "branch_freeze_hash": branch["branch_freeze_hash"],
        "coefficients": coeffs,
        "gates": gates,
        "formula_audit": formula,
    }
    observable_packet["packet_hash"] = sha256_json(observable_packet)

    # V2-21-compatible reduced manifest: three grouped lanes share isotropic data.
    lane_data = {
        "K": coeffs["wall"]["K"],
        "M": coeffs["wall"]["M"],
        "B0": coeffs["bdg"]["B0"],
        "B2": coeffs["bdg"]["B2"],
        "B4": coeffs["bdg"]["B4"],
        "Z0": coeffs["mixed"]["Z0"],
        "Z2": coeffs["mixed"]["Z2"],
        "Z4": coeffs["mixed"]["Z4"],
        "N0": coeffs["mixed"]["N0"],
        "N2": coeffs["mixed"]["N2"],
        "N4": coeffs["mixed"]["N4"],
    }
    v21_manifest = {
        "schema": "stage_v2_21_branch_manifest_generated_by_v2_23/v1",
        "source": "stage_v2_23_minimal_open_throat_branch_solver",
        "branch_id": branch["branch_id"],
        "branch_freeze_hash": branch["branch_freeze_hash"],
        "geometry": {
            "boundary_class": branch["geometry"]["boundary_class"],
            "R_exit": branch["geometry"]["R_exit"],
            "R_mouth": branch["geometry"]["R_mouth"],
        },
        "constants": branch["constants"],
        "grouped_lanes": {"20": lane_data, "21": lane_data, "22": lane_data},
        "isotropic_lane": lane_data,
    }
    v21_manifest["manifest_hash"] = sha256_json(v21_manifest)

    tolerance_report = {
        "schema": "stage_v2_23_tolerance_report/v1",
        "branch_id": branch["branch_id"],
        "branch_freeze_hash": branch["branch_freeze_hash"],
        "mode_residuals": {
            "wall_l2": branch["mode_diagnostics"]["wall"]["fem_residual_relative"],
            "bdg": branch["mode_diagnostics"]["bdg"]["fem_residual_relative"],
            "U": branch["mode_diagnostics"]["U"]["fem_residual_relative"],
            "W": branch["mode_diagnostics"]["W"]["fem_residual_relative"],
            "shape_stationarity": branch["geometry"]["stationary_residual_relative"],
        },
        "gates": gates,
    }
    tolerance_report["report_hash"] = sha256_json(tolerance_report)

    # Write JSON artifacts.
    (ARTIFACT_DIR / "stage_v2_23_reduced_solver_branch_manifest.json").write_text(
        json.dumps(branch, indent=2), encoding="utf-8"
    )
    (ARTIFACT_DIR / "stage_v2_23_generated_v21_manifest.json").write_text(
        json.dumps(v21_manifest, indent=2), encoding="utf-8"
    )
    (ARTIFACT_DIR / "stage_v2_23_observable_packet.json").write_text(
        json.dumps(observable_packet, indent=2), encoding="utf-8"
    )
    (ARTIFACT_DIR / "stage_v2_23_tolerance_report.json").write_text(
        json.dumps(tolerance_report, indent=2), encoding="utf-8"
    )

    # Human-readable output.
    lines: List[str] = []
    lines.append("Stage V2-23: minimal open-throat branch solver")
    lines.append("=" * 64)
    lines.append(f"formula_sympy_audit: {formula['passed']}/{formula['total']} checks passed")
    lines.append(f"branch_freeze_hash: {branch['branch_freeze_hash']}")
    lines.append(f"observable_packet_hash: {observable_packet['packet_hash']}")
    lines.append("")
    lines.append("Geometry / open-throat gate")
    lines.append(f"  boundary_class: {branch['geometry']['boundary_class']}")
    lines.append(f"  R_mouth: {branch['geometry']['R_mouth']:.12g}")
    lines.append(f"  R_exit: {branch['geometry']['R_exit']:.12g}")
    lines.append(f"  R_min: {branch['geometry']['R_min']:.12g}")
    lines.append(f"  shape_stationary_residual_relative: {branch['geometry']['stationary_residual_relative']:.3e}")
    lines.append(f"  open_gate_pass: {gates['open_gate_pass']}")
    lines.append("")
    lines.append("Mode eigenvalues")
    lines.append(f"  lambda_wall_l2 K: {coeffs['wall']['K']:.12g}")
    lines.append(f"  lambda_bdg varpi^2: {coeffs['bdg']['varpi2']:.12g}")
    lines.append(f"  lambda_U Omega_U^2: {coeffs['mixed']['Omega_U2']:.12g}")
    lines.append(f"  lambda_W Omega_W^2: {coeffs['mixed']['Omega_W2']:.12g}")
    lines.append("")
    lines.append("Overlap integrals")
    for k, v in coeffs["overlaps"].items():
        lines.append(f"  {k}: {v:.12g}")
    lines.append("")
    lines.append("Reduced coefficients")
    lines.append(f"  B0,B2,B4: {coeffs['bdg']['B0']:.12g}, {coeffs['bdg']['B2']:.12g}, {coeffs['bdg']['B4']:.12g}")
    lines.append(f"  Z0,Z2,Z4: {coeffs['mixed']['Z0']:.12g}, {coeffs['mixed']['Z2']:.12g}, {coeffs['mixed']['Z4']:.12g}")
    lines.append(f"  N0,N2,N4: {coeffs['mixed']['N0']:.12g}, {coeffs['mixed']['N2']:.12g}, {coeffs['mixed']['N4']:.12g}")
    lines.append(f"  D0,D2,D4: {coeffs['operator']['D0']:.12g}, {coeffs['operator']['D2']:.12g}, {coeffs['operator']['D4']:.12g}")
    lines.append("")
    lines.append("Response / prefactor")
    lines.append(f"  u2: {coeffs['response']['u2']:.12g}")
    lines.append(f"  u4: {coeffs['response']['u4']:.12g}")
    lines.append(f"  P0: {coeffs['prefactor']['P0']:.12g}")
    lines.append(f"  P2: {coeffs['prefactor']['P2']:.12g}")
    lines.append(f"  P4: {coeffs['prefactor']['P4']:.12g}")
    lines.append("")
    lines.append("Target residuals")
    lines.append(f"  P0_target: {coeffs['targets']['P0_target']:.12g}")
    lines.append(f"  R_pole: {coeffs['targets']['R_pole']:.12g}")
    lines.append(f"  R_norm: {coeffs['targets']['R_norm']:.12g}")
    lines.append(f"  R_P2: {coeffs['targets']['R_P2']:.12g}")
    lines.append(f"  R_P4: {coeffs['targets']['R_P4']:.12g}")
    lines.append(f"  R_tail: {coeffs['targets']['R_tail']:.12g}")
    lines.append(f"  gamma_eff: {coeffs['targets']['gamma_eff']:.12g}")
    lines.append(f"  gamma_GR: {coeffs['targets']['gamma_GR']:.12g}")
    lines.append("")
    lines.append("Gate verdict")
    lines.append(f"  stability_gate_pass: {gates['stability_gate_pass']}")
    lines.append(f"  outgoing_transfer_gate_pass: {gates['outgoing_transfer_gate_pass']}")
    lines.append(f"  target_packet_pass: {gates['target_packet_pass']}")
    lines.append(f"  dominant_residual: {gates['failure_diagnosis']['dominant_residual']}")
    lines.append(f"  P0/P0_target: {gates['failure_diagnosis']['normalization_ratio_P0_over_target']:.12g}")
    lines.append(f"  D0C/(3A^2): {gates['failure_diagnosis']['one_pole_ratio_D0C_over_3A2']:.12g}")

    output_text = "\n".join(lines) + "\n"
    print(output_text)


if __name__ == "__main__":
    main()
