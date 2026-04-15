#!/usr/bin/env python3
"""
Damped dressing-leg continuation of the reduced same-charge barrier model.

This script corrects the previous "perfect energy trap" treatment of the
dressing leg V by adding a Langevin-style shedding term. It keeps the same
reduced stationary barrier branch used in the earlier relaxed-constraint and
dynamic/helicity runs, and then does three things:

1) Cold crossing setup at v0 = 2.6 on the aligned-spin lowered-barrier branch.
2) Envelope-based damping correction for the collapse time:
       t_collapse,0 = sqrt(mu_eta / (g_UV * chi_lambda_peak))
       gamma_crit   = 1 / t_collapse,0
       t_collapse^damped(gamma_tot) = 1 / (gamma_crit - gamma_tot)
   for gamma_tot < gamma_crit, with unconditional stability for
   gamma_tot >= gamma_crit.
3) Time-domain dissipation audit of the V-leg along the actual v0 = 2.6
   trajectory using
       mu_eta Vdd + mu_eta gamma_tot Vdot + a_V V = - g_UV chi(r(t)) U(r(t)),
   so the event-dissipated energy is
       E_diss = mu_eta gamma_tot int Vdot(t)^2 dt,
   partitioned as vacuum vs lattice by the chosen rate ratio.

Important scope note
--------------------
The "critical shedding rate" requested by the assignment is an energy-envelope
criterion that matches the earlier undamped collapse-time estimate.
A purely mechanical damped inverted oscillator would not produce a finite
gamma_crit, so the script reports the envelope criterion explicitly and uses it
for the stability sweep.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Tuple

import math

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp


# ---------------------------------------------------------------------------
# Symbolic formulas
# ---------------------------------------------------------------------------

def derive_symbolic_formulas() -> Dict[str, sp.Expr]:
    gamma_tot, gamma_crit = sp.symbols("gamma_tot gamma_crit", positive=True, real=True)
    gamma_vac, gamma_lattice = sp.symbols("gamma_vac gamma_lattice", nonnegative=True, real=True)
    mu_eta, g_uv, chi_peak = sp.symbols("mu_eta g_UV chi_peak", positive=True, real=True)
    lam_eff, m_s, E, V_peak = sp.symbols("lambda_eff m_s E V_peak", positive=True, real=True)
    V0, v0 = sp.symbols("V0 v0", real=True)
    t = sp.symbols("t", real=True)
    Vdot = sp.Function("Vdot")(t)

    t_cross_E = sp.simplify(lam_eff * sp.sqrt(m_s / (2 * (E - V_peak))))
    t_cross_v = sp.simplify(lam_eff * sp.sqrt(m_s / (m_s * v0**2 + 2 * V0 - 2 * V_peak)))
    gamma_crit_expr = sp.simplify(sp.sqrt(g_uv * chi_peak / mu_eta))
    t_collapse_damped = sp.Piecewise(
        (1 / (gamma_crit - gamma_tot), gamma_tot < gamma_crit),
        (sp.oo, True),
    )
    gamma_safe = sp.simplify(gamma_crit - 1 / t_cross_v)
    E_diss = sp.simplify(mu_eta * gamma_tot * sp.Integral(Vdot**2, (t, 0, sp.Symbol("t_cross", positive=True))))
    E_vac = sp.simplify((gamma_vac / gamma_tot) * E_diss)
    E_lattice = sp.simplify((gamma_lattice / gamma_tot) * E_diss)

    return {
        "t_cross(E)": t_cross_E,
        "t_cross(v0)": t_cross_v,
        "gamma_crit": gamma_crit_expr,
        "t_collapse_damped": t_collapse_damped,
        "gamma_safe(v0)": gamma_safe,
        "E_diss": E_diss,
        "E_vac": E_vac,
        "E_lattice": E_lattice,
    }


# ---------------------------------------------------------------------------
# Parameter blocks
# ---------------------------------------------------------------------------

@dataclass
class StationaryParams:
    # geometric / leakage
    lam: float = 1.0
    r_reg: float = 0.25
    E0_amp: float = 0.18
    rho0: float = 1.0
    mu_w: float = 0.8
    q: float = 1.0

    # rigid-mouth U/V branch
    a_U: float = 2.5
    a_V: float = 3.0
    g_UV: float = 0.95
    s0: float = 0.9
    epsilon_ref: float = 0.30

    # one-port bundle
    K_star: float = 4.0
    Omega_U_sq: float = 9.0
    Omega_W_sq: float = 16.0
    G_U: float = 1.0
    G_W: float = 1.25
    R_mix: float = 1.35

    # static primitive source amplitudes
    beta_Q: float = 0.03
    beta_U0: float = 0.15
    beta_W0: float = 0.20
    kappa: float = 1.0

    # sign-changing multimode amplitudes
    a0: float = 2.2
    b0: float = -0.6
    r_sigma: float = 0.8
    xi_R: float = 0.9

    # open-system reductions into effective brane potential
    eta_leak: float = 0.03
    eta_uv: float = 0.22
    kk_amp: float = 0.0

    # reference Family-1 value
    r_F1: float = 1.77799353547498


@dataclass
class DynamicParams:
    # reduced-unit scattering setup carried from the earlier dynamic run
    m_s: float = 1.0
    hbar_eff: float = 1.0
    r0: float = 5.0
    r_contact: float = 0.18
    E_sub: float = 2.5

    # aligned-spin mixed-force closure
    xi_spin: float = 0.40
    C0_spin: float = 0.05
    r_spin: float = 0.80
    r_core: float = 0.15

    # integrator
    dt: float = 1.0e-4
    tmax: float = 10.0


@dataclass
class SheddingParams:
    # corrected cold run
    v0_cold: float = 2.6

    # dressing leg
    mu_eta: float = 1.0
    vacuum_fraction: float = 0.25  # default intrinsic-vacuum share of gamma_tot
    gamma_scan_factor_max: float = 1.40
    n_gamma_points: int = 321

    # width used in the characteristic crossing-time formula
    # "model" uses the raw reduced width lambda=1.0
    # "trigger" uses lambda_th extracted from the earlier subbarrier turning point
    lambda_choice: str = "model"


# ---------------------------------------------------------------------------
# Stationary reduced closure carried forward from the earlier runs
# ---------------------------------------------------------------------------

def gradient_numeric(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.gradient(y, x)


def local_slope(x: float, xs: np.ndarray, ys: np.ndarray) -> float:
    if x <= xs[0]:
        return float((ys[1] - ys[0]) / (xs[1] - xs[0]))
    if x >= xs[-1]:
        return float((ys[-1] - ys[-2]) / (xs[-1] - xs[-2]))
    i = int(np.searchsorted(xs, x))
    x0, x1 = xs[i - 1], xs[i]
    y0, y1 = ys[i - 1], ys[i]
    return float((y1 - y0) / (x1 - x0))


def find_turning_points(xs: np.ndarray, V: np.ndarray, E: float) -> List[float]:
    y = V - E
    roots: List[float] = []
    for i in range(len(xs) - 1):
        y0, y1 = y[i], y[i + 1]
        if y0 == 0.0:
            roots.append(float(xs[i]))
            continue
        if y0 * y1 < 0.0:
            x0, x1 = xs[i], xs[i + 1]
            roots.append(float(x0 + (0.0 - y0) * (x1 - x0) / (y1 - y0)))
    return roots


def compute_uv_branch(r: np.ndarray, p: StationaryParams) -> Dict[str, np.ndarray]:
    chi = p.lam / (r + p.r_reg)
    s = p.s0 * np.exp(-r / p.lam)
    den = p.a_U * p.a_V - p.g_UV**2 * chi**2
    if np.any(den <= 0):
        raise ValueError("U/V branch denominator crossed zero.")
    U = p.a_V * s / den
    V = -p.g_UV * chi * s / den
    epsilon_eta = p.epsilon_ref * np.exp(V)
    Xi_1 = U
    energy_drop = (chi**2 * p.g_UV**2 * s**2) / (2.0 * p.a_U * den)
    return {
        "chi_lambda": chi,
        "U": U,
        "V": V,
        "Xi_1": Xi_1,
        "epsilon_eta": epsilon_eta,
        "UV_energy_drop": energy_drop,
    }


def compute_leakage_diagnostics(r: np.ndarray, p: StationaryParams) -> Dict[str, np.ndarray]:
    E0 = p.E0_amp / (r + p.r_reg) ** 2
    S_leak = p.rho0 * p.mu_w * E0 / (np.sqrt(2.0 * np.pi) * p.lam**3)
    work_bulk = 2.0 * p.q * p.rho0 * p.mu_w * E0**2 / p.lam**2
    return {"E0": E0, "S_leak": S_leak, "Work_bulk": work_bulk}


def compute_sigma_branch(r: np.ndarray, p: StationaryParams) -> Dict[str, np.ndarray]:
    a_amp = p.a0 / (1.0 + (r / p.r_sigma) ** 2)
    b_amp = p.b0 / (1.0 + (r / p.r_sigma) ** 2)
    g_sigma = 2.0 / np.pi + 2.0 * a_amp / (3.0 * np.pi) - 2.0 * b_amp / (15.0 * np.pi)
    R_sigma = ((g_sigma - p.r_F1) ** 2) / (1.0 + p.r_F1**2)
    return {"a_amp": a_amp, "b_amp": b_amp, "g_sigma": g_sigma, "R_sigma": R_sigma}


def one_port_susceptibilities(p: StationaryParams) -> Dict[str, float]:
    Delta = p.Omega_U_sq * p.Omega_W_sq - p.R_mix**2
    Q = p.G_U**2 * p.Omega_W_sq + 2.0 * p.G_U * p.G_W * p.R_mix + p.G_W**2 * p.Omega_U_sq
    P = p.Omega_U_sq * p.G_W + p.R_mix * p.G_U
    P_U = p.G_U * p.Omega_W_sq + p.R_mix * p.G_W
    D0 = p.K_star - Q / Delta
    if Delta <= 0.0 or D0 <= 0.0:
        raise ValueError("One-port bundle left the admissible branch: require Delta>0 and D0>0.")
    chi_qq = 1.0 / D0
    chi_qU = P_U / (Delta * D0)
    chi_qW = P / (Delta * D0)
    chi_UU = (p.K_star * p.Omega_W_sq - p.G_W**2) / (Delta * D0)
    chi_UW = (p.K_star * p.R_mix + p.G_U * p.G_W) / (Delta * D0)
    chi_WW = (p.K_star * p.Omega_U_sq - p.G_U**2) / (Delta * D0)
    return {
        "Delta": Delta,
        "Q": Q,
        "P": P,
        "P_U": P_U,
        "D0": D0,
        "chi_qq": chi_qq,
        "chi_qU": chi_qU,
        "chi_qW": chi_qW,
        "chi_UU": chi_UU,
        "chi_UW": chi_UW,
        "chi_WW": chi_WW,
    }


def compute_static_mixed_potential(
    r: np.ndarray,
    p: StationaryParams,
    sigma_branch: Dict[str, np.ndarray],
    sus: Dict[str, float],
) -> Dict[str, np.ndarray]:
    beta_U = p.beta_U0 * (1.0 + p.xi_R * (sigma_branch["R_sigma"] - 0.25))
    beta_W = p.beta_W0 * (1.0 + p.xi_R * (sigma_branch["R_sigma"] - 0.25))
    beta_Q = p.beta_Q * np.ones_like(r)

    C6 = sus["chi_qq"] * beta_Q**2
    C4 = sus["chi_qU"] * beta_Q * beta_U + sus["chi_qW"] * beta_Q * beta_W
    C2 = sus["chi_UU"] * beta_U**2 + 2.0 * sus["chi_UW"] * beta_U * beta_W + sus["chi_WW"] * beta_W**2

    V_static = -0.5 * (
        C6 / r**6
        + 2.0 * C4 * np.exp(-2.0 * p.kappa * r) / r**4
        + C2 * np.exp(-4.0 * p.kappa * r) / r**2
    )
    return {
        "beta_Q": beta_Q,
        "beta_U": beta_U,
        "beta_W": beta_W,
        "V_static_mix": V_static,
    }


def compute_effective_potential(r: np.ndarray, p: StationaryParams) -> Dict[str, np.ndarray]:
    uv = compute_uv_branch(r, p)
    leak = compute_leakage_diagnostics(r, p)
    sigma_branch = compute_sigma_branch(r, p)
    sus = one_port_susceptibilities(p)
    mixed = compute_static_mixed_potential(r, p, sigma_branch, sus)

    V_coulomb = 1.0 / r
    V_kk = p.kk_amp * np.exp(-2.0 * r / p.lam) / r
    V_eff = (
        V_coulomb
        + V_kk
        + mixed["V_static_mix"]
        - p.eta_leak * leak["Work_bulk"]
        - p.eta_uv * uv["UV_energy_drop"]
    )
    return {
        "r": r,
        "V_coulomb": V_coulomb,
        "V_eff": V_eff,
        **uv,
        **leak,
        **sigma_branch,
        **sus,
        **mixed,
    }


# ---------------------------------------------------------------------------
# Dynamic radial trajectory on the aligned-spin branch
# ---------------------------------------------------------------------------

def reduced_spin_force(r: float, stat: StationaryParams, dyn: DynamicParams, s: int = +1) -> float:
    E0 = stat.E0_amp / (r + stat.r_reg) ** 2
    v_w = stat.mu_w * E0
    C_spin = -dyn.C0_spin * (1.0 + s * dyn.xi_spin * np.exp(-(r / dyn.r_spin) ** 2)) / (r + dyn.r_core) ** 2
    return stat.q * v_w * C_spin


def radial_acceleration(
    r: float,
    xs: np.ndarray,
    V: np.ndarray,
    add_spin: bool,
    s: int,
    stat: StationaryParams,
    dyn: DynamicParams,
) -> float:
    dVdr = local_slope(r, xs, V)
    accel = -(1.0 / dyn.m_s) * dVdr
    if add_spin:
        accel += reduced_spin_force(r, stat, dyn, s=s) / dyn.m_s
    return accel


def simulate_trajectory(
    r0: float,
    inward_speed: float,
    xs: np.ndarray,
    V: np.ndarray,
    stat: StationaryParams,
    dyn: DynamicParams,
    add_spin: bool = False,
    s: int = 0,
) -> Dict[str, np.ndarray | float | str]:
    r = float(r0)
    v = -float(inward_speed)
    t = 0.0

    ts = [t]
    rs = [r]
    vs = [v]
    status = "max_time"

    while t < dyn.tmax:
        a1 = radial_acceleration(r, xs, V, add_spin, s, stat, dyn)
        r_new = r + v * dyn.dt + 0.5 * a1 * dyn.dt * dyn.dt

        if r_new <= dyn.r_contact:
            t += dyn.dt
            v = v + a1 * dyn.dt
            r = r_new
            ts.append(t)
            rs.append(r)
            vs.append(v)
            status = "contact"
            break

        a2 = radial_acceleration(r_new, xs, V, add_spin, s, stat, dyn)
        v_new = v + 0.5 * (a1 + a2) * dyn.dt

        if v < 0.0 <= v_new:
            frac = 0.0 if v_new == v else (-v / (v_new - v))
            t = t + frac * dyn.dt
            r = r + frac * (r_new - r)
            ts.append(t)
            rs.append(r)
            vs.append(0.0)
            status = "turn"
            break

        t += dyn.dt
        r = r_new
        v = v_new
        ts.append(t)
        rs.append(r)
        vs.append(v)

    return {
        "status": status,
        "t": np.asarray(ts, dtype=float),
        "r": np.asarray(rs, dtype=float),
        "v": np.asarray(vs, dtype=float),
    }


# ---------------------------------------------------------------------------
# Damped V-leg model and cold sweep
# ---------------------------------------------------------------------------

def choose_lambda_eff(choice: str, station: StationaryParams, dyn: DynamicParams, r: np.ndarray, V_eff: np.ndarray) -> Tuple[float, str, float]:
    if choice == "model":
        roots_sub = find_turning_points(r, V_eff, dyn.E_sub)
        if len(roots_sub) < 2:
            raise RuntimeError("Could not identify the subbarrier turning point for the trigger-width cross-check.")
        lambda_trigger = float(abs(np.interp(roots_sub[-1], r, V_eff) / np.interp(roots_sub[-1], r, np.gradient(V_eff, r))))
        return station.lam, "raw model width", lambda_trigger

    if choice == "trigger":
        roots_sub = find_turning_points(r, V_eff, dyn.E_sub)
        if len(roots_sub) < 2:
            raise RuntimeError("Could not identify the subbarrier turning point used to define lambda_trigger.")
        r_turn = float(roots_sub[-1])
        dV = np.gradient(V_eff, r)
        lambda_trigger = float(abs(np.interp(r_turn, r, V_eff) / np.interp(r_turn, r, dV)))
        return lambda_trigger, "trigger width from the earlier subbarrier turning point", lambda_trigger

    raise ValueError("lambda_choice must be 'model' or 'trigger'.")


def solve_dressing_leg(
    traj: Dict[str, np.ndarray | float | str],
    gamma_tot: float,
    station: StationaryParams,
    U_traj: np.ndarray,
    chi_traj: np.ndarray,
    shed: SheddingParams,
) -> Dict[str, np.ndarray | float]:
    t = np.asarray(traj["t"], dtype=float)

    V_leg = np.zeros_like(t)
    Vdot_leg = np.zeros_like(t)
    E_diss = 0.0
    E_in = 0.0

    for i in range(len(t) - 1):
        dt = float(t[i + 1] - t[i])

        U_i = float(U_traj[i])
        chi_i = float(chi_traj[i])
        F_i = -station.g_UV * chi_i * U_i

        acc_i = (F_i - shed.mu_eta * gamma_tot * Vdot_leg[i] - station.a_V * V_leg[i]) / shed.mu_eta

        V_pred = V_leg[i] + Vdot_leg[i] * dt
        Vdot_pred = Vdot_leg[i] + acc_i * dt

        U_j = float(U_traj[i + 1])
        chi_j = float(chi_traj[i + 1])
        F_j = -station.g_UV * chi_j * U_j
        acc_j = (F_j - shed.mu_eta * gamma_tot * Vdot_pred - station.a_V * V_pred) / shed.mu_eta

        V_leg[i + 1] = V_leg[i] + 0.5 * (Vdot_leg[i] + Vdot_pred) * dt
        Vdot_leg[i + 1] = Vdot_leg[i] + 0.5 * (acc_i + acc_j) * dt

        E_diss += shed.mu_eta * gamma_tot * (Vdot_leg[i] ** 2) * dt
        E_in += F_i * Vdot_leg[i] * dt

    E_store = 0.5 * shed.mu_eta * Vdot_leg**2 + 0.5 * station.a_V * V_leg**2
    return {
        "t": t,
        "r": np.asarray(traj["r"], dtype=float),
        "V_leg": V_leg,
        "Vdot_leg": Vdot_leg,
        "E_store": E_store,
        "E_diss": float(E_diss),
        "E_in": float(E_in),
    }


def run_damped_shedding_cold_sweep(
    station: StationaryParams,
    dyn: DynamicParams,
    shed: SheddingParams,
) -> Dict[str, object]:
    r = np.linspace(dyn.r_contact, 8.0, 5000)
    stat_res = compute_effective_potential(r, station)
    V_eff = np.asarray(stat_res["V_eff"], dtype=float)
    dV = gradient_numeric(r, V_eff)
    dlnV = dV / V_eff

    i_peak = int(np.argmax(V_eff))
    r_peak = float(r[i_peak])
    V_peak = float(V_eff[i_peak])
    V0 = float(np.interp(dyn.r0, r, V_eff))
    vcrit_new = float(math.sqrt(2.0 * (V_peak - V0) / dyn.m_s))

    traj_new = simulate_trajectory(
        r0=dyn.r0,
        inward_speed=shed.v0_cold,
        xs=r,
        V=V_eff,
        stat=station,
        dyn=dyn,
        add_spin=True,
        s=+1,
    )
    traj_coul = simulate_trajectory(
        r0=dyn.r0,
        inward_speed=shed.v0_cold,
        xs=r,
        V=1.0 / r,
        stat=station,
        dyn=dyn,
        add_spin=False,
        s=0,
    )

    if str(traj_new["status"]) != "contact":
        raise RuntimeError("The chosen cold run did not cross on the lowered aligned-spin branch.")

    lambda_eff, lambda_label, lambda_trigger = choose_lambda_eff(shed.lambda_choice, station, dyn, r, V_eff)
    strip_mask = (r >= dyn.r_contact) & (r <= r_peak)
    chi_peak = float(lambda_eff * np.max(np.abs(dlnV[strip_mask])))

    E_cold = V0 + 0.5 * dyn.m_s * shed.v0_cold**2
    t_cross_char = float(lambda_eff * math.sqrt(dyn.m_s / (2.0 * (E_cold - V_peak))))
    gamma_crit = float(math.sqrt(station.g_UV * chi_peak / shed.mu_eta))
    t_collapse_0 = 1.0 / gamma_crit
    gamma_safe = max(0.0, gamma_crit - 1.0 / t_cross_char)

    U_traj = np.interp(np.asarray(traj_new["r"], dtype=float), r, np.asarray(stat_res["U"], dtype=float))
    chi_traj = np.interp(np.asarray(traj_new["r"], dtype=float), r, np.asarray(stat_res["chi_lambda"], dtype=float))

    safe_exact = solve_dressing_leg(traj_new, float(gamma_safe), station, U_traj, chi_traj, shed)
    crit_exact = solve_dressing_leg(traj_new, float(gamma_crit), station, U_traj, chi_traj, shed)

    gamma_grid = np.linspace(0.0, shed.gamma_scan_factor_max * gamma_crit, shed.n_gamma_points)
    t_collapse_damped = np.where(gamma_grid < gamma_crit, 1.0 / (gamma_crit - gamma_grid), np.inf)
    S_ratio = np.where(np.isfinite(t_collapse_damped), t_cross_char / t_collapse_damped, 0.0)

    E_diss = np.zeros_like(gamma_grid)
    E_store_final = np.zeros_like(gamma_grid)
    for i, g in enumerate(gamma_grid):
        sol = solve_dressing_leg(traj_new, float(g), station, U_traj, chi_traj, shed)
        E_diss[i] = float(sol["E_diss"])
        E_store_final[i] = float(sol["E_store"][-1])

    gamma_vac_grid = shed.vacuum_fraction * gamma_grid
    gamma_lat_grid = (1.0 - shed.vacuum_fraction) * gamma_grid
    E_vac = shed.vacuum_fraction * E_diss
    E_lat = (1.0 - shed.vacuum_fraction) * E_diss

    def idx_of_closest(x: np.ndarray, value: float) -> int:
        return int(np.argmin(np.abs(x - value)))

    i_safe = idx_of_closest(gamma_grid, gamma_safe)
    i_crit = idx_of_closest(gamma_grid, gamma_crit)

    t_total_cross = float(np.asarray(traj_new["t"], dtype=float)[-1])
    r_turn_coul = float(np.asarray(traj_coul["r"], dtype=float)[-1])

    return {
        "r": r,
        "stat_res": stat_res,
        "V_eff": V_eff,
        "dV": dV,
        "dlnV": dlnV,
        "r_peak": r_peak,
        "V_peak": V_peak,
        "V0": V0,
        "vcrit_new": vcrit_new,
        "traj_new": traj_new,
        "traj_coul": traj_coul,
        "r_turn_coul": r_turn_coul,
        "lambda_eff": lambda_eff,
        "lambda_label": lambda_label,
        "lambda_trigger": lambda_trigger,
        "chi_peak": chi_peak,
        "E_cold": E_cold,
        "t_cross_char": t_cross_char,
        "t_total_cross": t_total_cross,
        "gamma_crit": gamma_crit,
        "t_collapse_0": t_collapse_0,
        "gamma_safe": gamma_safe,
        "gamma_grid": gamma_grid,
        "t_collapse_damped": t_collapse_damped,
        "S_ratio": S_ratio,
        "E_diss": E_diss,
        "E_store_final": E_store_final,
        "gamma_vac_grid": gamma_vac_grid,
        "gamma_lat_grid": gamma_lat_grid,
        "E_vac": E_vac,
        "E_lat": E_lat,
        "i_safe": i_safe,
        "i_crit": i_crit,
        "safe_exact": safe_exact,
        "crit_exact": crit_exact,
    }


# ---------------------------------------------------------------------------
# Reporting and plotting
# ---------------------------------------------------------------------------

def build_report(
    station: StationaryParams,
    dyn: DynamicParams,
    shed: SheddingParams,
    formulas: Dict[str, sp.Expr],
    data: Dict[str, object],
) -> str:
    lines: List[str] = []
    lines.append("Damped dressing-leg correction for the reduced same-charge cold crossing")
    lines.append("=" * 92)
    lines.append("")
    lines.append("Scope")
    lines.append("-" * 92)
    lines.append("This run keeps the same reduced lowered-barrier branch used in the earlier")
    lines.append("dynamic/helicity continuation. The cold event is fixed at v0 = 2.6 on the")
    lines.append("aligned-spin branch, which still crosses the lowered barrier while pure Coulomb")
    lines.append("turns back. The only new ingredient is Langevin-style shedding on the V-leg.")
    lines.append("")
    lines.append("Symbolic formulas used")
    lines.append("-" * 92)
    lines.append(f"t_cross(E)         = {sp.simplify(formulas['t_cross(E)'])}")
    lines.append(f"t_cross(v0)        = {sp.simplify(formulas['t_cross(v0)'])}")
    lines.append(f"gamma_crit         = {sp.simplify(formulas['gamma_crit'])}")
    lines.append(f"t_collapse^damped  = {sp.simplify(formulas['t_collapse_damped'])}")
    lines.append(f"gamma_safe(v0)     = {sp.simplify(formulas['gamma_safe(v0)'])}")
    lines.append(f"E_diss             = {sp.simplify(formulas['E_diss'])}")
    lines.append(f"E_vac              = {sp.simplify(formulas['E_vac'])}")
    lines.append(f"E_lattice          = {sp.simplify(formulas['E_lattice'])}")
    lines.append("")
    lines.append("Parameter set")
    lines.append("-" * 92)
    for k, v in asdict(station).items():
        lines.append(f"{k:>20s} = {v}")
    for k, v in asdict(dyn).items():
        lines.append(f"{k:>20s} = {v}")
    for k, v in asdict(shed).items():
        lines.append(f"{k:>20s} = {v}")
    lines.append("")
    lines.append("Cold crossing setup")
    lines.append("-" * 92)
    lines.append(f"Barrier peak position r_peak                = {data['r_peak']:.8f}")
    lines.append(f"Barrier peak height V_peak                 = {data['V_peak']:.8f}")
    lines.append(f"Outer-point energy V_eff(r0)              = {data['V0']:.8f}")
    lines.append(f"Lowered-branch classical threshold v_crit = {data['vcrit_new']:.8f}")
    lines.append(f"Chosen cold speed v0_cold                 = {shed.v0_cold:.8f}")
    lines.append(f"Cold event energy E_cold                  = {data['E_cold']:.8f}")
    lines.append(f"Lowered-branch trajectory status          = {data['traj_new']['status']}")
    lines.append(f"Pure Coulomb trajectory status            = {data['traj_coul']['status']}")
    lines.append(f"Pure Coulomb turning radius at v0=2.6     = {data['r_turn_coul']:.8f}")
    lines.append(f"Actual lowered-branch time to contact     = {data['t_total_cross']:.8f}")
    lines.append("")
    lines.append("Characteristic crossing vs collapse")
    lines.append("-" * 92)
    lines.append(f"Width used in t_cross                      = {data['lambda_eff']:.8f} ({data['lambda_label']})")
    lines.append(f"Trigger-width cross-check                  = {data['lambda_trigger']:.8f}")
    lines.append(f"chi_lambda^peak                            = {data['chi_peak']:.8f}")
    lines.append(f"Undamped collapse time t_collapse,0       = {data['t_collapse_0']:.8f}")
    lines.append(f"Characteristic t_cross(v0=2.6)            = {data['t_cross_char']:.8f}")
    lines.append(f"Undamped stability ratio S(0)             = {data['t_cross_char'] / data['t_collapse_0']:.8f}")
    lines.append(f"Critical total shedding rate gamma_crit    = {data['gamma_crit']:.8f}")
    lines.append(f"Minimum shedding for this cold crossing    = {data['gamma_safe']:.8f}")
    lines.append("")
    lines.append("Interpretation:")
    lines.append("  - gamma_crit is the unconditional-stability threshold where the envelope")
    lines.append("    injection rate equals the shedding rate and t_collapse -> infinity.")
    lines.append("  - gamma_safe(v0=2.6) is the weaker condition needed only for this specific")
    lines.append("    cold event to satisfy S = t_cross / t_collapse^damped < 1.")
    lines.append("")
    lines.append("Vacuum vs lattice split (default fixed ratio)")
    lines.append("-" * 92)
    lines.append(f"gamma_vac / gamma_tot                      = {shed.vacuum_fraction:.8f}")
    lines.append(f"gamma_lattice / gamma_tot                  = {1.0 - shed.vacuum_fraction:.8f}")
    lines.append(f"At gamma_safe:   gamma_vac                 = {shed.vacuum_fraction * data['gamma_safe']:.8f}")
    lines.append(f"At gamma_safe:   gamma_lattice             = {(1.0 - shed.vacuum_fraction) * data['gamma_safe']:.8f}")
    lines.append(f"At gamma_crit:   gamma_vac                 = {shed.vacuum_fraction * data['gamma_crit']:.8f}")
    lines.append(f"At gamma_crit:   gamma_lattice             = {(1.0 - shed.vacuum_fraction) * data['gamma_crit']:.8f}")
    lines.append("")
    lines.append("Dissipated-energy audit along the actual v0 = 2.6 event")
    lines.append("-" * 92)
    lines.append(f"At gamma_safe:   E_diss,total              = {data['safe_exact']['E_diss']:.8f}")
    lines.append(f"At gamma_safe:   E_vac                     = {shed.vacuum_fraction * data['safe_exact']['E_diss']:.8f}")
    lines.append(f"At gamma_safe:   E_lattice                 = {(1.0 - shed.vacuum_fraction) * data['safe_exact']['E_diss']:.8f}")
    lines.append(f"At gamma_safe:   E_store,final             = {data['safe_exact']['E_store'][-1]:.8f}")
    lines.append("")
    lines.append(f"At gamma_crit:   E_diss,total              = {data['crit_exact']['E_diss']:.8f}")
    lines.append(f"At gamma_crit:   E_vac                     = {shed.vacuum_fraction * data['crit_exact']['E_diss']:.8f}")
    lines.append(f"At gamma_crit:   E_lattice                 = {(1.0 - shed.vacuum_fraction) * data['crit_exact']['E_diss']:.8f}")
    lines.append(f"At gamma_crit:   E_store,final             = {data['crit_exact']['E_store'][-1]:.8f}")
    lines.append("")
    lines.append("Bottom line")
    lines.append("-" * 92)
    lines.append("1) On the corrected damped-envelope branch, the cold v0 = 2.6 event is not")
    lines.append("   safe when the V-leg is treated as a perfect trap: the undamped ratio S is")
    lines.append("   far above one.")
    lines.append("2) A finite shedding window now exists. The cold event becomes stable once")
    lines.append("   gamma_tot exceeds gamma_safe, and it becomes unconditionally stable once")
    lines.append("   gamma_tot reaches gamma_crit.")
    lines.append("3) With the default 3:1 lattice-to-vacuum split, the lattice takes 75% of the")
    lines.append("   dissipated dressing-leg energy. The threshold-event lattice heat per crossing")
    lines.append("   is therefore the E_lattice value quoted above.")
    return "\n".join(lines)


def make_timescale_plot(outdir: Path, data: Dict[str, object]) -> Path:
    path = outdir / "damped_shedding_timescales.png"
    gamma = np.asarray(data["gamma_grid"], dtype=float)
    tc = np.asarray(data["t_collapse_damped"], dtype=float)
    mask = np.isfinite(tc)

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.plot(gamma[mask], tc[mask], label=r"$t_{\rm collapse}^{\rm damped}$")
    ax.axhline(float(data["t_cross_char"]), linestyle="--", linewidth=1.2, label=r"$t_{\rm cross}(v_0=2.6)$")
    ax.axvline(float(data["gamma_safe"]), linestyle=":", linewidth=1.0, label=r"$\gamma_{\rm safe}$")
    ax.axvline(float(data["gamma_crit"]), linestyle="-.", linewidth=1.0, label=r"$\gamma_{\rm crit}$")
    ax.fill_between(
        gamma[mask],
        1e-3,
        tc[mask],
        where=(tc[mask] > float(data["t_cross_char"])),
        alpha=0.15,
        label="stable cold-crossing region",
    )
    ax.set_xlabel(r"Total shedding rate $\gamma_{\rm tot}$")
    ax.set_ylabel("Timescale")
    ax.set_yscale("log")
    ax.set_title("Cold-crossing stability window after V-leg shedding is added")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def make_heat_partition_plot(outdir: Path, data: Dict[str, object]) -> Path:
    path = outdir / "damped_shedding_heat_partition.png"
    gamma = np.asarray(data["gamma_grid"], dtype=float)
    E_diss = np.asarray(data["E_diss"], dtype=float)
    E_vac = np.asarray(data["E_vac"], dtype=float)
    E_lat = np.asarray(data["E_lat"], dtype=float)

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.plot(gamma, E_diss, label=r"$E_{\rm diss,total}$")
    ax.plot(gamma, E_vac, label=r"$E_{\rm vac}$")
    ax.plot(gamma, E_lat, label=r"$E_{\rm lattice}$")
    ax.axvline(float(data["gamma_safe"]), linestyle=":", linewidth=1.0, label=r"$\gamma_{\rm safe}$")
    ax.axvline(float(data["gamma_crit"]), linestyle="-.", linewidth=1.0, label=r"$\gamma_{\rm crit}$")
    ax.set_xlabel(r"Total shedding rate $\gamma_{\rm tot}$")
    ax.set_ylabel("Dissipated energy per crossing event")
    ax.set_title("Vacuum vs lattice partition of V-leg shedding")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def main() -> None:
    outdir = Path(__file__).resolve().parent
    station = StationaryParams()
    dyn = DynamicParams()
    shed = SheddingParams()
    formulas = derive_symbolic_formulas()
    data = run_damped_shedding_cold_sweep(station, dyn, shed)

    report_text = build_report(station, dyn, shed, formulas, data)
    report_path = outdir / "damped_shedding_report.txt"
    report_path.write_text(report_text, encoding="utf-8")

    times_path = make_timescale_plot(outdir, data)
    heat_path = make_heat_partition_plot(outdir, data)

    print(report_text)
    print("")
    print("Artifacts")
    print("-" * 92)
    print(f"Report : {report_path}")
    print(f"Plot   : {times_path}")
    print(f"Plot   : {heat_path}")


if __name__ == "__main__":
    main()
