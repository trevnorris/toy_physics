#!/usr/bin/env python3
"""pathA_26 SymPy engine: throat-soliton Derrick and open-drain stability.

This is a deliberately cheap existence test.  It derives the sharp-wall
collective energy, runs a coarse fixed-regulator envelope grid, checks the
pre-registered conservative rescue terms, and linearizes the slaved open
drain closure.  The script writes the required report and YAML results.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp
import yaml
from scipy.optimize import root


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
PROJECT_ROOT = STAGE1_ROOT.parents[1].parents[0]
REPORTS = STAGE1_ROOT / "reports"
SCRATCH = STAGE1_ROOT / "_scratch"
REPORT_OUT = REPORTS / "pathA_26_derrick.md"
YAML_OUT = REPORTS / "pathA_26_derrick_results.yaml"
JSON_OUT = SCRATCH / "pathA_26_derrick_sympy.json"
MMA_JSON = SCRATCH / "pathA_26_derrick_mathematica.json"


@dataclass(frozen=True)
class Params:
    # Natural-unit nondimensional execution point.  The dimensional audit in
    # the report restores units; these values only select the able-to-fail
    # numerical slice.
    rho0: float = 1.0
    K_eos: float = 0.002
    P_vac: float = 0.001
    sigma: float = 0.003
    I_wave: float = 1.0
    c_w: float = 1.0
    mu_in: float = 0.2
    mu_out: float = 2.5
    chi_perp: float = math.pi
    chi_DN: float = math.pi / 2.0
    hbar: float = 0.03
    m: float = 1.0
    delta_r: float = 0.05
    delta_parallel: float = 0.05
    eps_w: float = 0.01
    depletion_eta_min: float = 0.95
    numeric_grid_n: int = 48
    numeric_hessian_step: float = 0.05
    phase_b_max_coeff: float = 1.0e-4
    d_slab: float = 2.5
    c_damp_sample: float = 1.0
    c_damp_max: float = 10.0
    gain_min: float = 0.1
    gain_max: float = 10.0
    robust_axis_fraction: float = 0.30

    @property
    def P0(self) -> float:
        return self.K_eos * self.rho0**5

    @property
    def alpha_volume(self) -> float:
        return self.P0 + self.P_vac

    @property
    def mu_fluid(self) -> float:
        return 5.0 * self.K_eos * self.rho0**4 / 4.0


P = Params()


def clean_float(x: Any, digits: int = 12) -> float:
    return float(round(float(x), digits))


def clean_complex(z: complex, digits: int = 12) -> dict[str, float]:
    return {"re": clean_float(z.real, digits), "im": clean_float(z.imag, digits)}


def volume(a: Any, L: Any) -> Any:
    return sp.Rational(4, 3) * sp.pi * a**3 * L


def area_side(a: Any, L: Any) -> Any:
    return 4 * sp.pi * a**2 * L


def area_cap(a: Any) -> Any:
    return sp.Rational(8, 3) * sp.pi * a**3


def wave_omega(a: Any, L: Any) -> Any:
    return sp.sqrt(
        P.c_w**2 * (P.chi_perp**2 / a**2 + P.chi_DN**2 / L**2) + P.mu_in**2
    )


def wave_omega_num(a: float, L: float) -> float:
    return math.sqrt(
        P.c_w**2 * (P.chi_perp**2 / a**2 + P.chi_DN**2 / L**2) + P.mu_in**2
    )


def sharp_energy_expr(a: Any, L: Any) -> Any:
    return (
        P.alpha_volume * volume(a, L)
        + P.sigma * (area_side(a, L) + area_cap(a))
        + P.I_wave * wave_omega(a, L)
    )


def sharp_energy_num(q: np.ndarray) -> float:
    a, L = [float(v) for v in q]
    if a <= 0.0 or L <= 0.0:
        return float("inf")
    V = 4.0 * math.pi * a**3 * L / 3.0
    A = 4.0 * math.pi * a**2 * L + 8.0 * math.pi * a**3 / 3.0
    omega = wave_omega_num(a, L)
    return P.alpha_volume * V + P.sigma * A + P.I_wave * omega


def build_symbolic_baseline() -> dict[str, Any]:
    a, L = sp.symbols("a L", positive=True)
    E = sharp_energy_expr(a, L)
    grad = [sp.diff(E, a), sp.diff(E, L)]
    sol = sp.nsolve(grad, (a, L), (1.9, 1.8), tol=1e-28, maxsteps=100)
    q = np.array([float(sol[0]), float(sol[1])], dtype=float)
    H = sp.hessian(E, (a, L))
    H_num = np.array(H.subs({a: sol[0], L: sol[1]}), dtype=float)
    eig = np.linalg.eigvalsh(H_num)
    grad_num = [float(g.subs({a: sol[0], L: sol[1]})) for g in grad]
    omega_num = float(wave_omega(sol[0], sol[1]))

    Vn = float(volume(sol[0], sol[1]))
    As = float(area_side(sol[0], sol[1]))
    Ac = float(area_cap(sol[0]))
    E_vol = P.alpha_volume * Vn
    E_side = P.sigma * As
    E_cap = P.sigma * Ac
    wave_a = P.I_wave * P.c_w**2 * P.chi_perp**2 / (q[0] ** 2 * omega_num)
    wave_L = P.I_wave * P.c_w**2 * P.chi_DN**2 / (q[1] ** 2 * omega_num)
    virial_a = 3.0 * E_vol + 2.0 * E_side + 3.0 * E_cap - wave_a
    virial_L = E_vol + E_side - wave_L

    if np.linalg.norm(grad_num) > 1e-8 or np.any(q <= 0.0):
        label = "A_FORBIDDEN"
    elif float(np.min(eig)) <= 1e-10:
        label = "A_DEGENERATE"
    else:
        label = "A_CANDIDATE_MIN"

    return {
        "q_star": [clean_float(q[0]), clean_float(q[1])],
        "L_over_a": clean_float(q[1] / q[0]),
        "energy": clean_float(sharp_energy_num(q)),
        "omega": clean_float(omega_num),
        "gradient_residual": [clean_float(x, 14) for x in grad_num],
        "hessian": [[clean_float(x) for x in row] for row in H_num.tolist()],
        "hessian_eigenvalues": [clean_float(x) for x in eig.tolist()],
        "label": label,
        "virial": {
            "a_identity": "3*E_volume + 2*E_side + 3*E_cap = I*c_w^2*chi_perp^2/(a^2*omega)",
            "L_identity": "E_volume + E_side = I*c_w^2*chi_DN^2/(L^2*omega)",
            "a_residual": clean_float(virial_a, 13),
            "L_residual": clean_float(virial_L, 13),
        },
    }


def fluid_only_contrast() -> dict[str, Any]:
    # With no wave term and positive couplings, both first derivatives are
    # strictly positive for a,L>0.  The throat closes toward the boundary.
    a, L = sp.symbols("a L", positive=True)
    E = P.alpha_volume * volume(a, L) + P.sigma * (area_side(a, L) + area_cap(a))
    da = sp.factor(sp.diff(E, a))
    dL = sp.factor(sp.diff(E, L))
    positive_sample = [
        float(da.subs({a: 1.0, L: 1.0})),
        float(dL.subs({a: 1.0, L: 1.0})),
    ]
    no_stationary = all(x > 0.0 for x in positive_sample)
    return {
        "label": "FLUID_ONLY_COLLAPSE_NO_INTERIOR_STATIONARY" if no_stationary else "FLUID_ONLY_CHECK_FAILED",
        "dE_da": str(da),
        "dE_dL": str(dL),
        "positive_derivative_sample": [clean_float(x) for x in positive_sample],
    }


class SmoothEnvelopeGrid:
    def __init__(self, q0: list[float]):
        a0, L0 = q0
        self.nr = P.numeric_grid_n
        self.nw = P.numeric_grid_n
        self.rmax = a0 + 0.60
        self.wmax = L0 / 2.0 + 0.60
        r = np.linspace(0.0, self.rmax, self.nr)
        w = np.linspace(-self.wmax, self.wmax, self.nw)
        self.R, self.W = np.meshgrid(r, w, indexing="ij")
        wr = np.ones(self.nr)
        ww = np.ones(self.nw)
        wr[[0, -1]] = 0.5
        ww[[0, -1]] = 0.5
        dr = self.rmax / (self.nr - 1)
        dw = 2.0 * self.wmax / (self.nw - 1)
        self.weights = wr[:, None] * ww[None, :] * dr * dw * 4.0 * math.pi * self.R**2

    @staticmethod
    def sech2(x: np.ndarray) -> np.ndarray:
        return 1.0 / np.cosh(np.clip(x, -40.0, 40.0)) ** 2

    def fluid_excess(self, a: float, L: float) -> float:
        eta = P.depletion_eta_min
        u = (self.R - a) / P.delta_r
        Gr = 0.5 * (1.0 - np.tanh(u))
        dGr = -0.5 * self.sech2(u) / P.delta_r
        smooth_abs = np.sqrt(self.W**2 + P.eps_w**2)
        v = (smooth_abs - L / 2.0) / P.delta_parallel
        Gw = 0.5 * (1.0 - np.tanh(v))
        dGw = -0.5 * self.sech2(v) * (self.W / smooth_abs) / P.delta_parallel
        G = Gr * Gw
        dG2 = (dGr * Gw) ** 2 + (Gr * dGw) ** 2
        rho = P.rho0 * (1.0 - eta * G)
        rho_clip = np.maximum(rho, 1.0e-5)
        grand = P.K_eos * rho**5 / 4.0 - P.mu_fluid * rho
        grand -= P.K_eos * P.rho0**5 / 4.0 - P.mu_fluid * P.rho0
        kinetic = P.hbar**2 * (P.rho0 * eta) ** 2 * dG2 / (8.0 * P.m * rho_clip)
        return float(np.sum((grand + kinetic) * self.weights))

    def total_energy(self, q: np.ndarray) -> float:
        a, L = [float(x) for x in q]
        if a <= 0.0 or L <= 0.0:
            return float("inf")
        # The smooth grid explicitly computes the EOS pressure depletion part,
        # so only the vacuum part of the PV geometry is added here.
        V = 4.0 * math.pi * a**3 * L / 3.0
        A = 4.0 * math.pi * a**2 * L + 8.0 * math.pi * a**3 / 3.0
        omega = wave_omega_num(a, L)
        return self.fluid_excess(a, L) + P.P_vac * V + P.sigma * A + P.I_wave * omega

    def hessian_at(self, q: list[float]) -> np.ndarray:
        qv = np.array(q, dtype=float)
        h = P.numeric_hessian_step
        H = np.zeros((2, 2), dtype=float)
        for i in range(2):
            for j in range(2):
                ei = np.zeros(2)
                ej = np.zeros(2)
                ei[i] = h
                ej[j] = h
                H[i, j] = (
                    self.total_energy(qv + ei + ej)
                    - self.total_energy(qv + ei - ej)
                    - self.total_energy(qv - ei + ej)
                    + self.total_energy(qv - ei - ej)
                ) / (4.0 * h * h)
        return H


def numeric_envelope_check(q_star: list[float], analytic_eigs: list[float]) -> dict[str, Any]:
    grid = SmoothEnvelopeGrid(q_star)
    H = grid.hessian_at(q_star)
    eig = np.linalg.eigvalsh(H)
    sign_agree = bool(np.all(eig > 0.0) == np.all(np.array(analytic_eigs) > 0.0))
    a0, L0 = q_star
    a_values = np.linspace(a0 - 2.0 * P.numeric_hessian_step, a0 + 2.0 * P.numeric_hessian_step, 5)
    L_values = np.linspace(L0 - 2.0 * P.numeric_hessian_step, L0 + 2.0 * P.numeric_hessian_step, 5)
    mesh = []
    for av in a_values:
        row = []
        for Lv in L_values:
            row.append(clean_float(grid.total_energy(np.array([av, Lv])), 10))
        mesh.append(row)
    return {
        "method": "coarse smooth-gate fixed-mu depletion envelope; eta relaxed to constrained lower bound",
        "grid_shape": [P.numeric_grid_n, P.numeric_grid_n],
        "eta_relaxed": P.depletion_eta_min,
        "domain": {
            "rmax": clean_float(grid.rmax),
            "wmax": clean_float(grid.wmax),
        },
        "hessian": [[clean_float(x) for x in row] for row in H.tolist()],
        "hessian_eigenvalues": [clean_float(x) for x in eig.tolist()],
        "analytic_numeric_sign_agree": sign_agree,
        "candidate_interior_to_grid": bool(
            q_star[0] + 2.0 * P.numeric_hessian_step < grid.rmax
            and q_star[1] / 2.0 + 2.0 * P.numeric_hessian_step < grid.wmax
        ),
        "local_5x5_energy_mesh": mesh,
    }


def phase_b_energy_num(q: np.ndarray, coeffs: tuple[float, float, float]) -> float:
    a, L = [float(x) for x in q]
    b_flow, s_slab, lam_w = coeffs
    if a <= 0.0 or L <= 0.0:
        return float("inf")
    bend = (16.0 * math.pi / 9.0) * L
    return (
        sharp_energy_num(q)
        + b_flow / a**2
        + s_slab * (L / P.d_slab) ** 2
        + lam_w * bend
    )


def finite_hessian(func, q: np.ndarray, h: float = 1.0e-4) -> np.ndarray:
    H = np.zeros((2, 2), dtype=float)
    for i in range(2):
        for j in range(2):
            ei = np.zeros(2)
            ej = np.zeros(2)
            ei[i] = h
            ej[j] = h
            H[i, j] = (
                func(q + ei + ej)
                - func(q + ei - ej)
                - func(q - ei + ej)
                + func(q - ei - ej)
            ) / (4.0 * h * h)
    return H


def finite_grad(func, q: np.ndarray, h: float = 1.0e-6) -> np.ndarray:
    g = np.zeros(2, dtype=float)
    for i in range(2):
        ei = np.zeros(2)
        ei[i] = h
        g[i] = (func(q + ei) - func(q - ei)) / (2.0 * h)
    return g


def solve_phase_b_corner(coeffs: tuple[float, float, float], q_guess: list[float]) -> dict[str, Any]:
    func = lambda x: phase_b_energy_num(x, coeffs)
    sol = root(lambda x: finite_grad(func, x), np.array(q_guess, dtype=float), tol=1e-11)
    if not sol.success or np.any(sol.x <= 0.0):
        return {"success": False, "coeffs": list(coeffs), "message": sol.message}
    H = finite_hessian(func, sol.x)
    eig = np.linalg.eigvalsh(H)
    return {
        "success": bool(np.min(eig) > 0.0),
        "coeffs": [clean_float(x, 10) for x in coeffs],
        "q_star": [clean_float(x) for x in sol.x.tolist()],
        "hessian_eigenvalues": [clean_float(x) for x in eig.tolist()],
    }


def phase_b_analysis(q_star: list[float]) -> dict[str, Any]:
    maxc = P.phase_b_max_coeff
    corners = []
    for bf in (0.0, maxc):
        for ss in (0.0, maxc):
            for lw in (0.0, maxc):
                # The open region is the strict positive interior; zero corners
                # are included only as closure checks.
                corners.append(solve_phase_b_corner((bf, ss, lw), q_star))
    all_ok = all(corner.get("success", False) for corner in corners)
    hard_status = "interior_not_boundary" if q_star[1] < P.d_slab else "pushed_to_L_equals_d_slab"
    label = "B_RESCUABLE" if all_ok and hard_status == "interior_not_boundary" else "B_NOT_RESCUABLE"
    return {
        "label": label,
        "open_positive_region": {
            "B_flow": f"(0, {maxc:g})",
            "sigma_ret": f"(0, {maxc:g})",
            "lambda_W": f"(0, {maxc:g})",
            "d_slab": P.d_slab,
        },
        "corner_checks": corners,
        "hard_slab": {
            "constraint": f"L <= {P.d_slab:g}",
            "status": hard_status,
        },
        "soft_slab_form": "sigma_ret*(L/d_slab)^2",
        "flow_form": "Phi^2/(8*pi^2*rho*a^2)",
        "bending_form": "(16*pi/9)*lambda_W*L plus fixed-regulator edge corrections",
    }


def conservative_grad(q: np.ndarray) -> np.ndarray:
    a, L = [float(x) for x in q]
    omega = wave_omega_num(a, L)
    dVa = 4.0 * math.pi * a**2 * L
    dVL = 4.0 * math.pi * a**3 / 3.0
    dAa = 8.0 * math.pi * a * L + 8.0 * math.pi * a**2
    dAL = 4.0 * math.pi * a**2
    dwa = -P.c_w**2 * P.chi_perp**2 / (a**3 * omega)
    dwL = -P.c_w**2 * P.chi_DN**2 / (L**3 * omega)
    return np.array(
        [
            P.alpha_volume * dVa + P.sigma * dAa + P.I_wave * dwa,
            P.alpha_volume * dVL + P.sigma * dAL + P.I_wave * dwL,
        ],
        dtype=float,
    )


def conservative_hessian(q: np.ndarray) -> np.ndarray:
    return finite_hessian(sharp_energy_num, q, h=1.0e-5)


def phase_c_force(q: np.ndarray, g_open: float, sign: int) -> np.ndarray:
    a, L = [float(x) for x in q]
    return sign * g_open * np.array([2.0 * a**7 / L, a**8 / L**2], dtype=float)


def phase_c_force_jacobian(q: np.ndarray, g_open: float, sign: int) -> np.ndarray:
    a, L = [float(x) for x in q]
    return sign * g_open * np.array(
        [
            [14.0 * a**6 / L, -2.0 * a**7 / L**2],
            [8.0 * a**7 / L**2, -2.0 * a**8 / L**3],
        ],
        dtype=float,
    )


def collective_mass(q: np.ndarray) -> np.ndarray:
    a, L = [float(x) for x in q]
    V = 4.0 * math.pi * a**3 * L / 3.0
    return np.diag([P.rho0 * V, P.rho0 * V / 4.0])


def phase_c_fixed_point(g_open: float, sign: int, q_guess: list[float]) -> np.ndarray | None:
    guesses = [
        np.array(q_guess, dtype=float),
        np.array([1.0, 1.0]),
        np.array([2.0, 2.0]),
        np.array([2.0, 4.0]),
        np.array([1.5, 1.2]),
        np.array([3.0, 3.0]),
    ]
    for guess in guesses:
        sol = root(lambda x: conservative_grad(x) - phase_c_force(x, g_open, sign), guess)
        if sol.success and np.all(sol.x > 0.0):
            residual = np.linalg.norm(conservative_grad(sol.x) - phase_c_force(sol.x, g_open, sign))
            if residual < 1.0e-7 and np.all(sol.x < 50.0):
                return sol.x
    return None


def phase_c_jacobian(q: np.ndarray, g_open: float, sign: int, c_value: float) -> tuple[np.ndarray, np.ndarray]:
    H = conservative_hessian(q)
    K_open = -phase_c_force_jacobian(q, g_open, sign)
    stiffness = H + K_open
    C = np.diag([c_value, c_value])
    Minv = np.linalg.inv(collective_mass(q))
    J = np.block(
        [
            [np.zeros((2, 2)), np.eye(2)],
            [-Minv @ stiffness, -Minv @ C],
        ]
    )
    return J, stiffness


def phase_c_stable_at(g_open: float, sign: int, q_guess: list[float]) -> tuple[bool, dict[str, Any] | None]:
    q = phase_c_fixed_point(g_open, sign, q_guess)
    if q is None:
        return False, None
    J, stiffness = phase_c_jacobian(q, g_open, sign, P.c_damp_sample)
    eig = np.linalg.eigvals(J)
    sym_stiffness_eig = np.linalg.eigvalsh((stiffness + stiffness.T) / 2.0)
    stable = bool(np.max(eig.real) < -1.0e-9)
    return stable, {
        "q_star": q,
        "jacobian_eigenvalues": eig,
        "det_total_stiffness": float(np.linalg.det(stiffness)),
        "stiffness_symmetric_eigenvalues": sym_stiffness_eig,
    }


def phase_c_static_divergence(q_star: list[float], g_open: float) -> dict[str, Any]:
    q = np.array(q_star, dtype=float)
    signs: dict[str, Any] = {}
    for sign in (+1, -1):
        J, stiffness = phase_c_jacobian(q, g_open, sign, P.c_damp_sample)
        sym_stiffness_eig = np.linalg.eigvalsh((stiffness + stiffness.T) / 2.0)
        total_stiffness_eig = np.linalg.eigvals(stiffness)
        signs["plus" if sign == 1 else "minus"] = {
            "s": sign,
            "det_total_stiffness": clean_float(float(np.linalg.det(stiffness)), 6),
            "min_symmetric_stiffness_eigenvalue": clean_float(sym_stiffness_eig[0], 6),
            "symmetric_stiffness_eigenvalues": [
                clean_float(x, 6) for x in sym_stiffness_eig.tolist()
            ],
            "total_stiffness_eigenvalues": [
                clean_complex(complex(z), 6) for z in total_stiffness_eig.tolist()
            ],
            "max_jacobian_real_part": clean_float(float(np.max(np.linalg.eigvals(J).real)), 6),
        }
    return {
        "reference_point": "conservative_phase_A_q_star",
        "q_reference": [clean_float(x) for x in q.tolist()],
        "g_open": clean_float(g_open, 10),
        "damping_note": (
            "Static divergence: H+K_open has negative determinant and the symmetric "
            "stiffness part has a negative eigenvalue, so passive C>=0 can only set "
            "timescales; it cannot remove the negative-stiffness direction."
        ),
        "signs": signs,
    }


def find_gcrit(sign: int, q_star: list[float]) -> float:
    # Find the edge of the tiny stable corner in the scalar product
    # g_open=(kappa_c*mu_drive)^2/rho0.  It is enough for the robustness
    # falsifier to bracket the first loss of LHP stability or steady existence.
    lo = 1.0e-4
    hi = 1.0
    stable_lo, _ = phase_c_stable_at(lo, sign, q_star)
    if not stable_lo:
        return lo
    while True:
        stable_hi, _ = phase_c_stable_at(hi, sign, q_star)
        if not stable_hi:
            break
        hi *= 2.0
        if hi > 100.0:
            return hi
    for _ in range(50):
        mid = math.sqrt(lo * hi)
        stable_mid, _ = phase_c_stable_at(mid, sign, q_star)
        if stable_mid:
            lo = mid
        else:
            hi = mid
    return lo


def phase_c_analysis(q_star: list[float]) -> dict[str, Any]:
    signs: dict[str, Any] = {}
    g_sample = (P.gain_min * P.gain_min) ** 2 / P.rho0
    for sign in (+1, -1):
        stable_sample, sample = phase_c_stable_at(g_sample, sign, q_star)
        if sample is None:
            sample_payload: dict[str, Any] = {"steady_point": None, "stable": False}
        else:
            sample_payload = {
                "steady_point": [clean_float(x) for x in sample["q_star"].tolist()],
                "jacobian_spectrum": [
                    clean_complex(complex(z), 10) for z in sample["jacobian_eigenvalues"].tolist()
                ],
                "det_total_stiffness": clean_float(sample["det_total_stiffness"], 10),
                "stiffness_symmetric_eigenvalues": [
                    clean_float(x, 10) for x in sample["stiffness_symmetric_eigenvalues"].tolist()
                ],
                "stable": stable_sample,
            }
        gcrit = find_gcrit(sign, q_star)
        signs["plus" if sign == 1 else "minus"] = {
            "s": sign,
            "sample_gain": {
                "kappa_c": P.gain_min,
                "mu_drive": P.gain_min,
                "zeta": 1.0,
                "c_a": P.c_damp_sample,
                "c_L": P.c_damp_sample,
                "g_open": clean_float(g_sample, 10),
            },
            "sample": sample_payload,
            "stable_corner_gcrit": clean_float(gcrit, 10),
            "stable_corner_condition": f"(kappa_c*mu_drive)^2 < {gcrit:.10g}",
        }

    required_radius = P.robust_axis_fraction * (P.gain_max - P.gain_min)
    min_center_on_any_robust_ball = P.gain_min + required_radius
    robust_center_g = (min_center_on_any_robust_ball**2) ** 2 / P.rho0
    min_upper_on_any_robust_interval = P.gain_min + 2.0 * required_radius
    robust_forced_g = (min_upper_on_any_robust_interval**2) ** 2 / P.rho0
    max_gcrit = max(signs["plus"]["stable_corner_gcrit"], signs["minus"]["stable_corner_gcrit"])
    robust_exists = bool(max_gcrit >= robust_forced_g)
    label = "C_STABILIZABLE" if robust_exists else "C_GENERICALLY_UNSTABLE"
    static_divergence = phase_c_static_divergence(q_star, robust_center_g)
    return {
        "label": label,
        "state_variables": ["Phi_w", "Delta_h"],
        "slaving": {
            "G_c": "a^3/L",
            "Phi_w": "kappa_c*(a^3/L)*mu_drive",
            "Delta_h": "zeta*Phi_w",
        },
        "work_kernel": "G_work=a^2*L/ell0^5 (nondimensional execution: a^2*L)",
        "F_nc_after_slaving": "s*(kappa_c*mu_drive)^2/rho0 * (2*a^7/L, a^8/L^2)",
        "K_open": "-dF_nc/dq",
        "mass_matrix": "diag(rho0*V, rho0*V/4)",
        "damping": "diag(c_a,c_L), c_a,c_L>=0",
        "gain_box": {
            "kappa_c": [P.gain_min, P.gain_max],
            "zeta": [P.gain_min, P.gain_max],
            "mu_drive": [P.gain_min, P.gain_max],
            "c_a_c_L": [0.0, P.c_damp_max],
        },
        "robust_required_axis_radius": clean_float(required_radius),
        "robust_min_center_gain": clean_float(min_center_on_any_robust_ball),
        "robust_center_g_open_lower_bound": clean_float(robust_center_g),
        "robust_forced_g_lower_bound": clean_float(robust_forced_g),
        "robust_region_exists": robust_exists,
        "static_divergence_at_box_forced_gain": static_divergence,
        "signs": signs,
        "instability_proof": (
            "Stable fixed points exist only in a tiny near-zero-drain corner. "
            "The most favorable in-box 30%-axis-radius ball center already has "
            f"kappa_c=mu_drive={min_center_on_any_robust_ball:.6g}, hence "
            f"g_open={robust_center_g:.6g}, far above max(gcrit)={max_gcrit:.6g}; "
            "therefore no qualifying ball can be stable."
        ),
    }


def forms_and_dimensions() -> dict[str, Any]:
    return {
        "sharp_wall_forms": [
            {
                "name": "fluid_grand_potential_depletion",
                "form": "P(rho0)*(4*pi/3)*a^3*L",
                "exponents_a_L": [3, 1],
                "sign": "+",
                "citation": "notes/inner_throat/inner_throat_hard_mode.md:1544; notes/inner_throat/inner_throat_freeze_sheet.md:202",
            },
            {
                "name": "vacuum_pressure_geometry",
                "form": "P_vac*(4*pi/3)*a^3*L",
                "exponents_a_L": [3, 1],
                "sign": "+",
                "citation": "notes/inner_throat/inner_throat_freeze_sheet.md:191; notes/inner_throat/inner_throat_freeze_sheet.md:202",
            },
            {
                "name": "brane_tension_side",
                "form": "sigma*4*pi*a^2*L",
                "exponents_a_L": [2, 1],
                "sign": "+",
                "citation": "notes/inner_throat/inner_throat_freeze_sheet.md:191; notes/inner_throat/inner_throat_freeze_sheet.md:206",
            },
            {
                "name": "brane_tension_caps",
                "form": "sigma*(8*pi/3)*a^3",
                "exponents_a_L": [3, 0],
                "sign": "+",
                "citation": "notes/inner_throat/inner_throat_freeze_sheet.md:191; notes/inner_throat/inner_throat_freeze_sheet.md:209",
            },
            {
                "name": "fixed_action_wave",
                "form": "I*sqrt(c_w^2*(pi^2/a^2+(pi/2)^2/L^2)+mu_in^2)",
                "exponents_a_L": "algebraic; asymptotic pieces (-1,0) and (0,-1), not a monomial",
                "sign": "+",
                "citation": "notes/inner_throat/inner_throat_hard_mode.md:479; notes/lepton_mass_notes.md:303",
            },
            {
                "name": "phaseB_internal_4D_self_flow",
                "form": "Phi^2/(8*pi^2*rho*a^2)",
                "exponents_a_L": [-2, 0],
                "sign": "+",
                "citation": "notes/lepton_mass_notes.md:188; notes/lepton_mass_notes.md:201; notes/lepton_mass_notes.md:1047",
            },
            {
                "name": "phaseB_soft_slab",
                "form": "sigma_ret*(L/d_slab)^2",
                "exponents_a_L": [0, 2],
                "sign": "+",
                "citation": "docs/conceptual_foundation.md:321; docs/conceptual_foundation.md:325",
            },
            {
                "name": "phaseB_hard_slab",
                "form": "constraint L<=d_slab",
                "exponents_a_L": "inequality constraint",
                "sign": "constraint",
                "citation": "docs/conceptual_foundation.md:325; docs/conceptual_foundation.md:326",
            },
            {
                "name": "optional_Willmore_side",
                "form": "lambda_W*(16*pi/9)*L",
                "exponents_a_L": [0, 1],
                "sign": "+",
                "citation": "notes/inner_throat/inner_throat_hard_mode.md:586; notes/inner_throat/4d_next_steps.md:353",
            },
            {
                "name": "phaseC_conductance",
                "form": "G_c=a^3/L",
                "exponents_a_L": [3, -1],
                "sign": "+",
                "citation": "notes/inner_throat/4d_next_steps.md:1210; notes/inner_throat/inner_throat_freeze_sheet.md:200",
            },
            {
                "name": "phaseC_work_kernel",
                "form": "G_work=a^2*L/ell0^5",
                "exponents_a_L": [2, 1],
                "sign": "s in {+,-} enters F_nc",
                "citation": "notes/inner_throat/4d_next_steps.md:1224; notes/inner_throat/inner_throat_freeze_sheet.md:206",
            },
        ],
        "dimensional_check": {
            "fluid_terms": "P and P_vac have units E/L^4; V has L^4.",
            "surface_term": "sigma has E/L^3; area has L^3.",
            "wave_term": "I has action E*T; omega has T^-1.",
            "phaseB_flow": "rho*u^2 integrated over d^4X gives energy; cited closed form Phi^2/(8*pi^2*rho*a^2).",
            "phaseC": (
                "kappa_c carries the units needed for Phi_w=kappa_c*G_c*mu_drive; "
                "G_work is normalized by ell0^5 so dG_work/dq has L^-3 and "
                "(Phi_w^2/rho)*dG_work/dq is a generalized force E/L."
            ),
            "jacobian": "M^{-1}*(H+K_open) has T^-2 and M^{-1}*C has T^-1.",
        },
    }


def lambda_ray_regression() -> dict[str, Any]:
    a, A, B, C = sp.symbols("a A B C", positive=True)
    F = A / a + B / a**2 + C * a**3
    dF = sp.diff(F, a)
    identity = sp.factor(a * dF)
    expected = -A / a - 2 * B / a**2 + 3 * C * a**3
    if sp.simplify(identity - expected) != 0:
        raise AssertionError("Lambda-ray reduced-ledger identity failed")
    return {
        "object": "reduced ledger only, not the full 4D continuum",
        "F": "A/a + B/a^2 + C*a^3",
        "exponents": [-1, -2, 3],
        "stationarity_identity": "-E_w - 2*E_f + 3*E_PV = 0",
        "virial": "E_w + 2*E_f = 3*E_PV",
    }


def approximation_gate(baseline: dict[str, Any], numeric: dict[str, Any]) -> dict[str, Any]:
    a_star, L_star = baseline["q_star"]
    scale_ratio = min(a_star / P.delta_r, L_star / P.delta_parallel)
    binding_margin = (P.mu_out - baseline["omega"]) / P.mu_out
    checks = {
        "scale_separation": scale_ratio > 10.0,
        "binding": baseline["omega"] < P.mu_out and binding_margin > 0.10,
        "candidate_interior_to_numeric_grid": bool(numeric["candidate_interior_to_grid"]),
        "hessian_sign_agreement": bool(numeric["analytic_numeric_sign_agree"]),
    }
    return {
        "checks": checks,
        "passes": all(checks.values()),
        "scale_ratio_min": clean_float(scale_ratio),
        "binding_margin_fraction": clean_float(binding_margin),
        "remedy_if_failed": "shrink fixed regulators, raise outside threshold, enlarge grid, or run a higher-resolution envelope solve",
    }


def top_line_verdict(phase_a: str, phase_b: str, phase_c: str, gate: dict[str, Any]) -> str:
    cons = phase_a == "A_CANDIDATE_MIN" or phase_b == "B_RESCUABLE"
    drain_ok = phase_c in {"C_STABILIZABLE", "C_UNNEEDED"}
    if cons and not drain_ok and not gate["passes"]:
        return "INCONCLUSIVE"
    if cons and drain_ok and not gate["passes"]:
        return "INCONCLUSIVE"
    if (not cons) and (not drain_ok):
        return "THROAT_FORBIDDEN"
    if (not cons) and drain_ok:
        return "THROAT_NEEDS_DYNAMICS"
    if cons and (not drain_ok):
        return "THROAT_DRAIN_DESTABILIZED"
    return "THROAT_CANDIDATE"


def citations() -> dict[str, list[str]]:
    return {
        "fluid_GNLS_and_H": [
            "notes/inner_throat/inner_throat_hard_mode.md:331",
            "notes/inner_throat/inner_throat_hard_mode.md:365",
            "notes/inner_throat/inner_throat_hard_mode.md:1544",
        ],
        "wave_surrogate_and_H": [
            "notes/inner_throat/inner_throat_hard_mode.md:459",
            "notes/inner_throat/inner_throat_hard_mode.md:479",
            "notes/inner_throat/inner_throat_hard_mode.md:1540",
        ],
        "fixed_regulator_Vconf": [
            "notes/inner_throat/inner_throat_freeze_sheet.md:54",
            "notes/inner_throat/inner_throat_freeze_sheet.md:62",
            "notes/inner_throat/inner_throat_freeze_sheet.md:70",
            "notes/inner_throat/inner_throat_freeze_sheet.md:89",
            "notes/inner_throat/inner_throat_freeze_sheet.md:101",
            "notes/inner_throat/inner_throat_freeze_sheet.md:113",
        ],
        "geometry_measures": [
            "notes/inner_throat/inner_throat_freeze_sheet.md:196",
            "notes/inner_throat/inner_throat_freeze_sheet.md:200",
            "notes/inner_throat/inner_throat_freeze_sheet.md:202",
            "notes/inner_throat/inner_throat_freeze_sheet.md:206",
        ],
        "failure_and_hessian": [
            "notes/inner_throat/inner_throat_hard_mode.md:715",
            "notes/inner_throat/inner_throat_hard_mode.md:1714",
            "notes/inner_throat/inner_throat_hard_mode.md:1720",
            "notes/inner_throat/inner_throat_hard_mode.md:1727",
        ],
        "open_system": [
            "notes/inner_throat/4d_next_steps.md:1184",
            "notes/inner_throat/4d_next_steps.md:1210",
            "notes/inner_throat/4d_next_steps.md:1217",
            "notes/inner_throat/4d_next_steps.md:2045",
            "notes/inner_throat/4d_next_steps.md:3145",
            "notes/inner_throat/4d_next_steps.md:3246",
            "notes/inner_throat/4d_next_steps.md:3262",
        ],
        "flow_and_DN": [
            "notes/lepton_mass_notes.md:50",
            "notes/lepton_mass_notes.md:60",
            "notes/lepton_mass_notes.md:155",
            "notes/lepton_mass_notes.md:175",
            "notes/lepton_mass_notes.md:188",
            "notes/lepton_mass_notes.md:199",
            "notes/lepton_mass_notes.md:266",
            "notes/lepton_mass_notes.md:303",
            "notes/lepton_mass_notes.md:1047",
        ],
        "multi_brane_slab": [
            "docs/conceptual_foundation.md:321",
            "docs/conceptual_foundation.md:325",
            "docs/conceptual_foundation.md:326",
            "docs/conceptual_foundation.md:328",
        ],
        "scalar_surrogate_scope": [
            "docs/conceptual_foundation.md:274",
            "docs/conceptual_foundation.md:275",
            "docs/conceptual_foundation.md:280",
            "docs/conceptual_foundation.md:283",
        ],
    }


def compare_with_mathematica(current: dict[str, Any]) -> dict[str, Any]:
    if not MMA_JSON.exists():
        return {
            "status": "pending",
            "mathematica_json": str(MMA_JSON),
            "message": "Mathematica payload not present yet; rerun this script after pathA_26_derrick.wl.",
        }
    other = json.loads(MMA_JSON.read_text(encoding="utf-8"))
    mismatches: list[str] = []

    def cmp_exact(path: str, a: Any, b: Any) -> None:
        if a != b:
            mismatches.append(f"{path}: {a!r} != {b!r}")

    def cmp_float_list(path: str, a: list[float], b: list[float], tol: float = 5.0e-6) -> None:
        if len(a) != len(b):
            mismatches.append(f"{path}: length {len(a)} != {len(b)}")
            return
        for idx, (x, y) in enumerate(zip(a, b)):
            if abs(float(x) - float(y)) > tol:
                mismatches.append(f"{path}[{idx}]: {x} != {y}")

    cmp_exact("top_line", current["top_line_verdict"], other["top_line_verdict"])
    cmp_exact("phase_A", current["phase_A"]["label"], other["phase_A"]["label"])
    cmp_exact("phase_B", current["phase_B"]["label"], other["phase_B"]["label"])
    cmp_exact("phase_C", current["phase_C"]["label"], other["phase_C"]["label"])
    cmp_exact("lambda_ray", current["lambda_ray"]["virial"], other["lambda_ray"]["virial"])
    cmp_float_list(
        "baseline_hessian_eigs",
        current["phase_A"]["hessian_eigenvalues"],
        other["phase_A"]["hessian_eigenvalues"],
    )
    cmp_float_list(
        "numeric_hessian_eigs",
        current["numeric_envelope"]["hessian_eigenvalues"],
        other["numeric_envelope"]["hessian_eigenvalues"],
        tol=2.0e-5,
    )
    for sign_name in ("plus", "minus"):
        cmp_float_list(
            f"phase_C_{sign_name}_sample_q",
            current["phase_C"]["signs"][sign_name]["sample"]["steady_point"],
            other["phase_C"]["signs"][sign_name]["sample"]["steady_point"],
            tol=2.0e-5,
        )
        if abs(
            current["phase_C"]["signs"][sign_name]["stable_corner_gcrit"]
            - other["phase_C"]["signs"][sign_name]["stable_corner_gcrit"]
        ) > 2.0e-4:
            mismatches.append(f"phase_C_{sign_name}_gcrit mismatch")
        current_div = current["phase_C"]["static_divergence_at_box_forced_gain"]["signs"][sign_name]
        other_div = other["phase_C"]["static_divergence_at_box_forced_gain"]["signs"][sign_name]
        if abs(
            current_div["det_total_stiffness"] - other_div["det_total_stiffness"]
        ) > 1.0:
            mismatches.append(f"phase_C_{sign_name}_static_det mismatch")
        if abs(
            current_div["min_symmetric_stiffness_eigenvalue"]
            - other_div["min_symmetric_stiffness_eigenvalue"]
        ) > 1.0e-2:
            mismatches.append(f"phase_C_{sign_name}_static_min_eig mismatch")
    if abs(
        current["phase_C"]["robust_center_g_open_lower_bound"]
        - other["phase_C"]["robust_center_g_open_lower_bound"]
    ) > 1.0e-8:
        mismatches.append("phase_C_robust_center_g mismatch")

    return {
        "status": "AGREE" if not mismatches else "DISAGREE",
        "mathematica_json": str(MMA_JSON),
        "mismatches": mismatches,
    }


def build_report(results: dict[str, Any]) -> str:
    forms = results["forms"]["sharp_wall_forms"]
    form_lines = "\n".join(
        f"| {row['name']} | `{row['form']}` | `{row['exponents_a_L']}` | {row['sign']} | {row.get('citation', '')} |"
        for row in forms
    )
    phase_c = results["phase_C"]
    plus = phase_c["signs"]["plus"]
    minus = phase_c["signs"]["minus"]
    div = phase_c["static_divergence_at_box_forced_gain"]
    div_plus = div["signs"]["plus"]
    div_minus = div["signs"]["minus"]
    report = f"""{results['top_line_verdict']}

**pathA_26 Derrick + Open-System Stability Report**

Computed phase labels: Phase A = `{results['phase_A']['label']}`, Phase B = `{results['phase_B']['label']}`, Phase C = `{results['phase_C']['label']}`.
The top-line follows the directive decision table from those computed labels and the approximation gate.

**Functional And Sources**

The fluid GNLS/Hamiltonian and n=5 EOS are anchored at `notes/inner_throat/inner_throat_hard_mode.md:331`, `notes/inner_throat/inner_throat_hard_mode.md:365`, and `notes/inner_throat/inner_throat_hard_mode.md:1544`.
The scalar wave surrogate and wave Hamiltonian are anchored at `notes/inner_throat/inner_throat_hard_mode.md:459`, `notes/inner_throat/inner_throat_hard_mode.md:479`, and `notes/inner_throat/inner_throat_hard_mode.md:1540`.
The fixed smooth confinement gates/regulators are from `notes/inner_throat/inner_throat_freeze_sheet.md:54`, `notes/inner_throat/inner_throat_freeze_sheet.md:62`, `notes/inner_throat/inner_throat_freeze_sheet.md:70`, `notes/inner_throat/inner_throat_freeze_sheet.md:89`, `notes/inner_throat/inner_throat_freeze_sheet.md:101`, and `notes/inner_throat/inner_throat_freeze_sheet.md:113`.
The 4D tube measures are from `notes/inner_throat/inner_throat_freeze_sheet.md:196`, `notes/inner_throat/inner_throat_freeze_sheet.md:202`, `notes/inner_throat/inner_throat_freeze_sheet.md:206`, and `notes/inner_throat/inner_throat_freeze_sheet.md:213`.
Open-system flux/back-pressure diagnostics are from `notes/inner_throat/4d_next_steps.md:1184`, `notes/inner_throat/4d_next_steps.md:1210`, `notes/inner_throat/4d_next_steps.md:1217`, `notes/inner_throat/4d_next_steps.md:1224`, the fixed-density reservoir template from `notes/inner_throat/4d_next_steps.md:3145`, and the back-pressure observable from `notes/inner_throat/4d_next_steps.md:3246` and `notes/inner_throat/4d_next_steps.md:3262`.
The internal 4D self-flow and D/N ladder are from `notes/lepton_mass_notes.md:188`, `notes/lepton_mass_notes.md:201`, `notes/lepton_mass_notes.md:266`, and `notes/lepton_mass_notes.md:303`; the three distinct flux objects are separated at `notes/lepton_mass_notes.md:1043` and `notes/lepton_mass_notes.md:1047`.

**Object And Closures**

The tested object is `E*(a,L)=Omega_fluid,excess + I*omega(a,L) + P_vac*V + sigma*A`.
The fixed-mu fluid closure uses `mu=U'(rho0)` and the sharp-wall depletion contributes `+P(rho0)*V`; in the execution slice `P(rho0)={P.P0:g}`, `P_vac={P.P_vac:g}`, and `sigma={P.sigma:g}`.
The fixed-action wave branch is `omega=sqrt(c_w^2*(pi^2/a^2+(pi/2)^2/L^2)+mu_in^2)`, with `chi_perp=pi`, `chi_DN=pi/2`, `mu_in={P.mu_in:g}`, and the binding threshold `mu_out={P.mu_out:g}`.

**Forms**

| term | form | a,L exponents/form | sign | source |
| --- | --- | --- | --- | --- |
{form_lines}

The wave term is intentionally not monomialized.  The scalar surrogate, Maxwell `F^2`, and transverse brane-shear all share the same two-derivative dispersion structure for this scaling test; polarization count changes degeneracy, not the `a,L` powers.  The exception is a future vector/shear calculation with extra boundary-localized gauge constraints, which this necessary-condition test does not include.

**Phase A**

Sharp-wall stationary point: `a*={results['phase_A']['q_star'][0]}`, `L*={results['phase_A']['q_star'][1]}`, `L/a={results['phase_A']['L_over_a']}`.
Analytic envelope Hessian eigenvalues: `{results['phase_A']['hessian_eigenvalues']}`.
Virial residuals: `a={results['phase_A']['virial']['a_residual']}`, `L={results['phase_A']['virial']['L_residual']}`.

Coarse fixed-regulator envelope Hessian eigenvalues: `{results['numeric_envelope']['hessian_eigenvalues']}`.
The numeric grid is a cheap constrained smooth-gate depletion check, not a high-resolution PDE solve; it uses grid `{results['numeric_envelope']['grid_shape']}` and relaxed depletion `eta={results['numeric_envelope']['eta_relaxed']}`.
Fluid-only contrast: `{results['fluid_only']['label']}` because both positive-coupling derivatives are positive for `a,L>0`.

**Phase B**

Phase B label: `{results['phase_B']['label']}`.
The certified open positive region is `{results['phase_B']['open_positive_region']}`.
The hard slab check gives `{results['phase_B']['hard_slab']['status']}` for `L <= {P.d_slab:g}`.
The optional Willmore side term uses the averaged side mean curvature `H=2/(3a)`, giving `(16*pi/9)*lambda_W*L` plus fixed-regulator edge corrections, so it is not scale-invariant under independent `(a,L)` variation.

**Phase C**

The slaved drain closure is `Phi_w=kappa_c*(a^3/L)*mu_drive`, `Delta_h=zeta*Phi_w`.
The single work-kernel is `G_work=a^2*L/ell0^5`, so after slaving `F_nc=s*(kappa_c*mu_drive)^2/rho0*(2*a^7/L, a^8/L^2)`.
`K_open=-dF_nc/dq`, `M=diag(rho0*V,rho0*V/4)`, and passive damping is `diag(c_a,c_L)`.

For `s=+1`, the nonzero-drain sample fixed point is `{plus['sample']['steady_point']}` and the Jacobian spectrum is `{plus['sample']['jacobian_spectrum']}`.
For `s=-1`, the nonzero-drain sample fixed point is `{minus['sample']['steady_point']}` and the Jacobian spectrum is `{minus['sample']['jacobian_spectrum']}`.
The stable near-zero-drain thresholds are `gcrit_plus={plus['stable_corner_gcrit']}` and `gcrit_minus={minus['stable_corner_gcrit']}`, where `g=(kappa_c*mu_drive)^2/rho0`.

At the box-required center gain `g_open={div['g_open']}` evaluated at `q*={div['q_reference']}`, the Phase-C failure is static DIVERGENCE, not flutter: for `s=+1`, `det(H+K_open)={div_plus['det_total_stiffness']}` and `lambda_min(sym(H+K_open))={div_plus['min_symmetric_stiffness_eigenvalue']}`; for `s=-1`, `det(H+K_open)={div_minus['det_total_stiffness']}` and `lambda_min(sym(H+K_open))={div_minus['min_symmetric_stiffness_eigenvalue']}`.
Because the stiffness has a negative determinant and a negative symmetric-part eigenvalue at the required box gain, passive damping `C>=0` cannot remove the negative-stiffness direction; the damping value only sets timescales, so the damping axis need not be swept to `c_max`.
Robust stable region exists: `{phase_c['robust_region_exists']}`.
Instability proof: {phase_c['instability_proof']}
The conclusion is robust to the magnitude of `M`: `M` changes timescales in the Jacobian, while the no-robust-ball proof uses the steady point and open stiffness thresholds.

**Approximation Gate**

Gate checks: `{results['approximation_gate']['checks']}`.
Scale separation min ratio: `{results['approximation_gate']['scale_ratio_min']}`.
Binding margin fraction: `{results['approximation_gate']['binding_margin_fraction']}`.

**Lambda-Ray Regression**

The reduced-ledger regression reproduces `F=A/a+B/a^2+C*a^3` with exponents `{results['lambda_ray']['exponents']}` and virial `{results['lambda_ray']['virial']}`.
This is only the reduced ledger identity from `notes/lepton_mass_notes.md:50`, `notes/lepton_mass_notes.md:60`, and `notes/lepton_mass_notes.md:72`; it is not a claim that the full 4D continuum has the same Lambda-ray volume scaling.

**Engine Agreement**

Engine agreement status: `{results['engine_agreement']['status']}`.
Agreement details: `{results['engine_agreement']}`.

**Falsifier**

Named falsifier: `{results['falsifier']['name']}`.
Input that would flip the top-line: {results['falsifier']['description']}

**Scope**

This is a necessary-condition test over the two collective coordinates `(a,L)` and the named `j=0` wave branch.  A positive conservative result is not an existence proof.  The scalar wave stands in only for the two-derivative scaling of the trapped support sector.  The grid is coarse and regulator-fixed.  The Phase-C drain law is parameterized, not derived from a PDE.  Stability of the `(a,L)` family does not exclude off-family shape or field negative modes.
"""
    return report


def emit_outputs(results: dict[str, Any]) -> None:
    REPORTS.mkdir(parents=True, exist_ok=True)
    SCRATCH.mkdir(parents=True, exist_ok=True)
    JSON_OUT.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    YAML_OUT.write_text(yaml.safe_dump(results, sort_keys=False, width=120), encoding="utf-8")
    REPORT_OUT.write_text(build_report(results), encoding="utf-8")


def summarize(results: dict[str, Any]) -> str:
    plus = results["phase_C"]["signs"]["plus"]
    minus = results["phase_C"]["signs"]["minus"]
    div = results["phase_C"]["static_divergence_at_box_forced_gain"]
    div_plus = div["signs"]["plus"]
    div_minus = div["signs"]["minus"]
    return "\n".join(
        [
            f"top_line_verdict: {results['top_line_verdict']}",
            f"phase_labels: A={results['phase_A']['label']} B={results['phase_B']['label']} C={results['phase_C']['label']}",
            f"exponent_form_table: {[(row['name'], row['exponents_a_L']) for row in results['forms']['sharp_wall_forms']]}",
            f"envelope_hessian_eigenvalues: analytic={results['phase_A']['hessian_eigenvalues']} numeric={results['numeric_envelope']['hessian_eigenvalues']}",
            f"phase_C_jacobian_spectrum_plus: {plus['sample']['jacobian_spectrum']}",
            f"phase_C_jacobian_spectrum_minus: {minus['sample']['jacobian_spectrum']}",
            f"phase_C_stable_region: robust_exists={results['phase_C']['robust_region_exists']} gcrit_plus={plus['stable_corner_gcrit']} gcrit_minus={minus['stable_corner_gcrit']}",
            f"phase_C_static_divergence: g_open={div['g_open']} plus_det={div_plus['det_total_stiffness']} plus_min_sym_eig={div_plus['min_symmetric_stiffness_eigenvalue']} minus_det={div_minus['det_total_stiffness']} minus_min_sym_eig={div_minus['min_symmetric_stiffness_eigenvalue']}",
            f"engine_agreement: {results['engine_agreement']['status']}",
            f"named_falsifier: {results['falsifier']['name']}",
        ]
    )


def main() -> None:
    baseline = build_symbolic_baseline()
    fluid_only = fluid_only_contrast()
    numeric = numeric_envelope_check(baseline["q_star"], baseline["hessian_eigenvalues"])
    phase_b = phase_b_analysis(baseline["q_star"])
    phase_c = phase_c_analysis(baseline["q_star"])
    gate = approximation_gate(baseline, numeric)
    top = top_line_verdict(baseline["label"], phase_b["label"], phase_c["label"], gate)
    forms = forms_and_dimensions()
    lambda_ray = lambda_ray_regression()

    max_gcrit = max(
        phase_c["signs"]["plus"]["stable_corner_gcrit"],
        phase_c["signs"]["minus"]["stable_corner_gcrit"],
    )
    required_radius = P.robust_axis_fraction * (P.gain_max - P.gain_min)
    min_upper = P.gain_min + 2.0 * required_radius
    robust_g = (min_upper**2) ** 2 / P.rho0
    eta_required = max_gcrit / robust_g

    results: dict[str, Any] = {
        "engine": "sympy",
        "parameters": {
            key: clean_float(value) if isinstance(value, float) else value
            for key, value in P.__dict__.items()
        }
        | {
            "P0": clean_float(P.P0),
            "alpha_volume": clean_float(P.alpha_volume),
            "mu_fluid": clean_float(P.mu_fluid),
        },
        "citations": citations(),
        "forms": forms,
        "phase_A": baseline,
        "fluid_only": fluid_only,
        "numeric_envelope": numeric,
        "phase_B": phase_b,
        "phase_C": phase_c,
        "approximation_gate": gate,
        "lambda_ray": lambda_ray,
        "top_line_verdict": top,
        "falsifier": {
            "name": "FALSIFIER_C_ROBUST_DRAIN_BALL",
            "description": (
                "A measured non-conservative drain closure with a 30%-axis-radius LHP-stable "
                "gain ball for nonzero drain would change Phase C to C_STABILIZABLE.  "
                f"For this kernel that is equivalent to reducing the effective open-force "
                f"normalization below eta_open={eta_required:.6g}, or otherwise changing the "
                "measured work-kernel/Jacobian so the robust ball fits inside the O(1) box."
            ),
        },
    }
    results["engine_agreement"] = compare_with_mathematica(results)
    emit_outputs(results)
    print(summarize(results))


if __name__ == "__main__":
    main()
