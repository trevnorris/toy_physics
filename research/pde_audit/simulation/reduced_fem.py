"""Reduced finite-element primitives for target-blind simulation generation."""

from __future__ import annotations

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


def trapz(y: np.ndarray, x: np.ndarray) -> float:
    trapezoid = getattr(np, "trapezoid", np.trapz)
    return float(trapezoid(y, x))


def central_gradient(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    return np.gradient(y, x, edge_order=2)


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
    """Assemble 1D linear-FEM stiffness and mass matrices."""
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

    keep = list(range(1, n)) if dirichlet_left else list(range(n))
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

    if trapz(mu * full, s) < 0:
        full *= -1.0

    norm = math.sqrt(trapz(mu * full * full, s))
    if norm <= 0:
        raise RuntimeError(f"{label}: zero norm eigenvector")
    full /= norm

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
    """Solve a target-blind stationary open-throat shape profile."""
    n = len(s)
    L = float(s[-1] - s[0])
    x = s / L
    R_pref = a_mouth - (a_mouth - R_exit_pref) * x + 0.04 * np.sin(np.pi * x)

    K = np.zeros((n, n), dtype=float)
    b = np.zeros(n, dtype=float)

    for e in range(n - 1):
        h = s[e + 1] - s[e]
        ke_grad = (T_R / h) * np.array([[1.0, -1.0], [-1.0, 1.0]])
        ke_mass = (K_R * h / 6.0) * np.array([[2.0, 1.0], [1.0, 2.0]])
        rp = np.array([R_pref[e], R_pref[e + 1]])
        be = ke_mass @ rp
        idx = [e, e + 1]
        for aa in range(2):
            b[idx[aa]] += be[aa]
            for bb in range(2):
                K[idx[aa], idx[bb]] += ke_grad[aa, bb] + ke_mass[aa, bb]

    K[-1, -1] += Y_exit
    b[-1] += Y_exit * R_exit_pref

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
    residual[0] = 0.0
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
