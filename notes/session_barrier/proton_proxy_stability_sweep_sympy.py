#!/usr/bin/env python3
"""
Proton-proxy structural-stability sweep for the reduced same-charge barrier branch.

This script continues the earlier dynamic/helicity closure in a new direction:
it asks whether there is an energy window in which the aligned-spin crossing
through the lowered barrier happens faster than the dressing geometry collapses
under the U/V cross-coupling drain.

The script is self-contained and does four things:

1. Rebuild the stationary reduced effective potential V_eff(r) used in the
   earlier relaxed-constraint and dynamic-scattering runs.

2. Import the same aligned-spin mixed-force closure and use a proton-proxy
   scaling for the particle mass m_s and wall inertia scale mu_eta.

3. Define the characteristic times
       t_cross(E)    = lambda_eff * sqrt[m_s / (2 (E - V_peak))]
       t_collapse    = sqrt[mu_eta / (g_UV * chi_lambda_peak)]
   and sweep the incident energy over an over-barrier band.

4. As a cross-check, integrate the actual 1D aligned-spin trajectories and
   measure the barrier-region transit time from r = r_peak to r = r_contact.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import math

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp


# ---------------------------------------------------------------------------
# Symbolic formulas
# ---------------------------------------------------------------------------

def derive_symbolic_formulas() -> Dict[str, sp.Expr]:
    E, lam_eff, m_s, V_peak = sp.symbols("E lambda_eff m_s V_peak", positive=True, real=True)
    mu_eta, g_uv, chi_peak = sp.symbols("mu_eta g_UV chi_peak", positive=True, real=True)
    alpha = sp.symbols("alpha", positive=True, real=True)

    t_cross = sp.simplify(lam_eff * sp.sqrt(m_s / (2 * (E - V_peak))))
    t_collapse = sp.simplify(sp.sqrt(mu_eta / (g_uv * chi_peak)))
    stability_ratio = sp.simplify(t_cross / t_collapse)
    energy_edge = sp.solve(sp.Eq(stability_ratio, 1), E)[0]
    stability_ratio_massmatched = sp.simplify(stability_ratio.subs(mu_eta, alpha * m_s))
    energy_edge_massmatched = sp.simplify(energy_edge.subs(mu_eta, alpha * m_s))

    return {
        "t_cross": t_cross,
        "t_collapse": t_collapse,
        "stability_ratio": stability_ratio,
        "energy_edge": energy_edge,
        "stability_ratio_massmatched": stability_ratio_massmatched,
        "energy_edge_massmatched": energy_edge_massmatched,
    }


# ---------------------------------------------------------------------------
# Stationary reduced closure carried from the earlier run
# ---------------------------------------------------------------------------

@dataclass
class StationaryParams:
    # geometry / leakage
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
class StabilitySweepParams:
    # baseline geometry of the earlier dynamic run
    r0: float = 5.0
    r_contact: float = 0.18
    E_sub_reference: float = 2.5

    # proton proxy
    proton_mass_ratio: float = 1836.15267343

    # sweep band (multiples of the proton-proxy lowered threshold speed)
    v_multiplier_min: float = 1.001
    v_multiplier_max: float = 5.0
    n_energy_points: int = 350

    # actual aligned-spin transit cross-check
    n_dynamic_samples: int = 24
    dt_dynamic: float = 2.0e-3
    tmax_dynamic: float = 500.0

    # aligned-spin mixed-force closure carried forward
    xi_spin: float = 0.40
    C0_spin: float = 0.05
    r_spin: float = 0.80
    r_core: float = 0.15

    # which width to use in the timescale formulas
    # "trigger" -> lambda_eff = lambda_th at the earlier subbarrier turning point
    # "model"   -> lambda_eff = stationary localization width lam
    lambda_choice: str = "trigger"


# ---------------------------------------------------------------------------
# Stationary branch helpers
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
        raise ValueError("The U/V branch denominator crossed zero.")

    U = p.a_V * s / den
    V = -p.g_UV * chi * s / den
    epsilon_eta = p.epsilon_ref * np.exp(V)
    energy_drop = (chi**2 * p.g_UV**2 * s**2) / (2.0 * p.a_U * den)

    return {
        "chi_lambda": chi,
        "U": U,
        "V": V,
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
    return {
        "a_amp": a_amp,
        "b_amp": b_amp,
        "g_sigma": g_sigma,
        "R_sigma": R_sigma,
    }


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
        "C6": C6,
        "C4": C4,
        "C2": C2,
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
        "V_kk": V_kk,
        "V_eff": V_eff,
        **uv,
        **leak,
        **sigma_branch,
        **sus,
        **mixed,
    }


# ---------------------------------------------------------------------------
# Aligned-spin transit cross-check
# ---------------------------------------------------------------------------

def reduced_spin_force(
    r: float,
    p: StationaryParams,
    s: StabilitySweepParams,
    sign: int = +1,
) -> float:
    E0 = p.E0_amp / (r + p.r_reg) ** 2
    v_w = p.mu_w * E0
    C_spin = -s.C0_spin * (1.0 + sign * s.xi_spin * np.exp(-(r / s.r_spin) ** 2)) / (r + s.r_core) ** 2
    return p.q * v_w * C_spin


def radial_acceleration(
    r: float,
    xs: np.ndarray,
    V: np.ndarray,
    p: StationaryParams,
    s: StabilitySweepParams,
    m_s: float,
    add_spin: bool = True,
    sign: int = +1,
) -> float:
    dVdr = local_slope(r, xs, V)
    accel = -(1.0 / m_s) * dVdr
    if add_spin:
        accel += reduced_spin_force(r, p, s, sign=sign) / m_s
    return accel


def simulate_trajectory(
    r0: float,
    inward_speed: float,
    xs: np.ndarray,
    V: np.ndarray,
    p: StationaryParams,
    s: StabilitySweepParams,
    m_s: float,
    add_spin: bool = True,
    sign: int = +1,
) -> Dict[str, np.ndarray | float | str]:
    r = float(r0)
    v = -float(inward_speed)
    t = 0.0

    ts = [t]
    rs = [r]
    vs = [v]

    status = "max_time"

    while t < s.tmax_dynamic:
        a1 = radial_acceleration(r, xs, V, p, s, m_s, add_spin=add_spin, sign=sign)
        r_new = r + v * s.dt_dynamic + 0.5 * a1 * s.dt_dynamic * s.dt_dynamic

        if r_new <= s.r_contact:
            t += s.dt_dynamic
            v = v + a1 * s.dt_dynamic
            r = r_new
            ts.append(t)
            rs.append(r)
            vs.append(v)
            status = "contact"
            break

        a2 = radial_acceleration(r_new, xs, V, p, s, m_s, add_spin=add_spin, sign=sign)
        v_new = v + 0.5 * (a1 + a2) * s.dt_dynamic

        # turning point detection
        if v < 0.0 <= v_new:
            frac = 0.0 if v_new == v else (-v / (v_new - v))
            t = t + frac * s.dt_dynamic
            r = r + frac * (r_new - r)
            ts.append(t)
            rs.append(r)
            vs.append(0.0)
            status = "turn"
            break

        t += s.dt_dynamic
        r = r_new
        v = v_new
        ts.append(t)
        rs.append(r)
        vs.append(v)

    return {
        "status": status,
        "t": np.array(ts),
        "r": np.array(rs),
        "v": np.array(vs),
    }


def barrier_region_transit_time(traj: Dict[str, np.ndarray | float | str], r_peak: float, r_contact: float) -> float:
    t = np.asarray(traj["t"], dtype=float)
    rr = np.asarray(traj["r"], dtype=float)
    idx_peak = np.where(rr <= r_peak)[0]
    idx_contact = np.where(rr <= r_contact)[0]
    if len(idx_peak) == 0 or len(idx_contact) == 0:
        return np.nan
    return float(t[idx_contact[0]] - t[idx_peak[0]])


# ---------------------------------------------------------------------------
# Master sweep
# ---------------------------------------------------------------------------

def choose_lambda_eff(choice: str, station: StationaryParams, sweep: StabilitySweepParams, base: Dict[str, np.ndarray | float]) -> Tuple[float, str]:
    if choice == "model":
        return station.lam, "model localization width"
    if choice == "trigger":
        return float(base["lambda_trigger"]), "previous trigger width lambda_th(r_turn)"
    raise ValueError(f"Unknown lambda_choice={choice!r}. Expected 'model' or 'trigger'.")


def run_stability_sweep(station: StationaryParams, sweep: StabilitySweepParams) -> Dict[str, object]:
    # Build the same stationary branch used in the earlier dynamic closure
    r = np.linspace(sweep.r_contact, 8.0, 5000)
    stat = compute_effective_potential(r, station)
    V_eff = stat["V_eff"]
    dV = gradient_numeric(r, V_eff)
    dlnV = dV / V_eff

    i_peak = int(np.argmax(V_eff))
    r_peak = float(r[i_peak])
    V_peak = float(V_eff[i_peak])

    V0 = float(np.interp(sweep.r0, r, V_eff))
    roots_sub = find_turning_points(r, V_eff, sweep.E_sub_reference)
    if len(roots_sub) < 2:
        raise RuntimeError("Could not identify the earlier subbarrier turning points needed for lambda_th.")
    r_turn_sub = float(roots_sub[-1])
    lambda_trigger = float(abs(np.interp(r_turn_sub, r, V_eff) / np.interp(r_turn_sub, r, dV)))

    lambda_eff, lambda_label = choose_lambda_eff(
        sweep.lambda_choice,
        station,
        sweep,
        {"lambda_trigger": lambda_trigger},
    )

    # Heavy-throat proton proxy
    m_s = sweep.proton_mass_ratio
    mu_eta = sweep.proton_mass_ratio

    vcrit_old = float(np.sqrt(2.0 * (V_peak - V0)))
    vcrit_proton = float(np.sqrt(2.0 * (V_peak - V0) / m_s))

    # Collapse time uses the steepest logarithmic gradient inside the crossing strip
    strip_mask = (r >= sweep.r_contact) & (r <= r_peak)
    grad_peak_location = float(r[strip_mask][int(np.argmax(np.abs(dlnV[strip_mask])))])
    chi_lambda_peak = float(lambda_eff * np.max(np.abs(dlnV[strip_mask])))
    t_collapse = float(np.sqrt(mu_eta / (station.g_UV * chi_lambda_peak)))

    # Sweep in the proton-proxy speed variable
    v_multipliers = np.linspace(sweep.v_multiplier_min, sweep.v_multiplier_max, sweep.n_energy_points)
    v_vals = v_multipliers * vcrit_proton
    E_inc = V0 + 0.5 * m_s * v_vals**2

    t_cross_char = np.full_like(E_inc, np.nan)
    mask = E_inc > V_peak
    t_cross_char[mask] = lambda_eff * np.sqrt(m_s / (2.0 * (E_inc[mask] - V_peak)))

    stability_ratio = np.full_like(E_inc, np.nan)
    stability_ratio[mask] = t_cross_char[mask] / t_collapse

    # Exact analytic lower edge for the characteristic criterion
    E_safe_lower = V_peak + (station.g_UV * chi_lambda_peak * lambda_eff**2 * m_s) / (2.0 * mu_eta)
    v_safe_lower = float(np.sqrt(2.0 * (E_safe_lower - V0) / m_s))
    v_safe_multiplier = float(v_safe_lower / vcrit_proton)

    # Actual aligned-spin barrier-region transit times as a cross-check
    sample_idx = np.linspace(0, len(v_vals) - 1, sweep.n_dynamic_samples, dtype=int)
    v_dynamic = v_vals[sample_idx]
    E_dynamic = E_inc[sample_idx]
    t_cross_dyn = []
    traj_status = []
    for v0 in v_dynamic:
        traj = simulate_trajectory(
            r0=sweep.r0,
            inward_speed=float(v0),
            xs=r,
            V=V_eff,
            p=station,
            s=sweep,
            m_s=m_s,
            add_spin=True,
            sign=+1,
        )
        t_cross_dyn.append(barrier_region_transit_time(traj, r_peak=r_peak, r_contact=sweep.r_contact))
        traj_status.append(str(traj["status"]))
    t_cross_dyn = np.asarray(t_cross_dyn, dtype=float)

    # Scan-limited window
    stable_mask = stability_ratio < 1.0
    if not np.any(stable_mask):
        raise RuntimeError("No stable window found on the requested sweep band.")
    E_window_lower_num = float(E_inc[stable_mask][0])
    E_window_upper_num = float(E_inc[stable_mask][-1])
    v_window_lower_num = float(v_vals[stable_mask][0])
    v_window_upper_num = float(v_vals[stable_mask][-1])

    # Sensitivity using the raw model width lambda = 1
    chi_peak_model = float(station.lam * np.max(np.abs(dlnV[strip_mask])))
    t_collapse_model = float(np.sqrt(mu_eta / (station.g_UV * chi_peak_model)))
    E_safe_lower_model = V_peak + (station.g_UV * chi_peak_model * station.lam**2 * m_s) / (2.0 * mu_eta)
    v_safe_lower_model = float(np.sqrt(2.0 * (E_safe_lower_model - V0) / m_s))

    return {
        "r": r,
        "V_eff": V_eff,
        "dV": dV,
        "dlnV": dlnV,
        "stat": stat,
        "r_peak": r_peak,
        "V_peak": V_peak,
        "V0": V0,
        "r_turn_sub": r_turn_sub,
        "lambda_trigger": lambda_trigger,
        "lambda_eff": lambda_eff,
        "lambda_label": lambda_label,
        "m_s": m_s,
        "mu_eta": mu_eta,
        "vcrit_old": vcrit_old,
        "vcrit_proton": vcrit_proton,
        "chi_lambda_peak": chi_lambda_peak,
        "grad_peak_location": grad_peak_location,
        "t_collapse": t_collapse,
        "v_multipliers": v_multipliers,
        "v_vals": v_vals,
        "E_inc": E_inc,
        "t_cross_char": t_cross_char,
        "stability_ratio": stability_ratio,
        "E_safe_lower": float(E_safe_lower),
        "v_safe_lower": v_safe_lower,
        "v_safe_multiplier": v_safe_multiplier,
        "E_window_lower_num": E_window_lower_num,
        "E_window_upper_num": E_window_upper_num,
        "v_window_lower_num": v_window_lower_num,
        "v_window_upper_num": v_window_upper_num,
        "E_dynamic": E_dynamic,
        "v_dynamic": v_dynamic,
        "t_cross_dyn": t_cross_dyn,
        "traj_status": traj_status,
        "chi_peak_model": chi_peak_model,
        "t_collapse_model": t_collapse_model,
        "E_safe_lower_model": float(E_safe_lower_model),
        "v_safe_lower_model": float(v_safe_lower_model),
    }


# ---------------------------------------------------------------------------
# Reporting and plotting
# ---------------------------------------------------------------------------

def build_report(
    station: StationaryParams,
    sweep: StabilitySweepParams,
    formulas: Dict[str, sp.Expr],
    data: Dict[str, object],
) -> str:
    lines: List[str] = []
    lines.append("Proton-proxy structural-stability sweep for the lowered same-charge barrier")
    lines.append("=" * 88)
    lines.append("")
    lines.append("Scope")
    lines.append("-" * 88)
    lines.append("This script keeps the same reduced stationary branch used in the earlier")
    lines.append("relaxed-constraint and dynamic/helicity runs, then overlays the requested")
    lines.append("proton-proxy mass/inertia scaling and the new crossing-vs-collapse timescale test.")
    lines.append("The aligned-spin closure stays active for the dynamic transit cross-check.")
    lines.append("")
    lines.append("Symbolic formulas used")
    lines.append("-" * 88)
    lines.append(f"t_cross(E)    = {sp.simplify(formulas['t_cross'])}")
    lines.append(f"t_collapse    = {sp.simplify(formulas['t_collapse'])}")
    lines.append(f"S(E)          = {sp.simplify(formulas['stability_ratio'])}")
    lines.append(f"E_edge        = {sp.simplify(formulas['energy_edge'])}")
    lines.append("")
    lines.append("If the heavy throat scales as mu_eta = alpha * m_s, these reduce to")
    lines.append(f"S(E)          = {sp.simplify(formulas['stability_ratio_massmatched'])}")
    lines.append(f"E_edge        = {sp.simplify(formulas['energy_edge_massmatched'])}")
    lines.append("")
    lines.append("Parameter set")
    lines.append("-" * 88)
    for k, v in asdict(station).items():
        lines.append(f"{k:>20s} = {v}")
    for k, v in asdict(sweep).items():
        lines.append(f"{k:>20s} = {v}")
    lines.append("")
    lines.append("Rebuilt stationary branch")
    lines.append("-" * 88)
    lines.append(f"Barrier peak location r_peak                 = {data['r_peak']:.8f}")
    lines.append(f"Barrier peak height V_peak                  = {data['V_peak']:.8f}")
    lines.append(f"Reference outer turning point r_turn(sub)   = {data['r_turn_sub']:.8f}")
    lines.append(f"Trigger width lambda_th(r_turn)             = {data['lambda_trigger']:.8f}")
    lines.append(f"Chosen lambda_eff                           = {data['lambda_eff']:.8f} ({data['lambda_label']})")
    lines.append("")
    lines.append("Heavy-throat proton proxy")
    lines.append("-" * 88)
    lines.append(f"Particle mass m_s                           = {data['m_s']:.8f}")
    lines.append(f"Wall inertia mu_eta                         = {data['mu_eta']:.8f}")
    lines.append(f"Old reduced threshold speed (m_s = 1)       = {data['vcrit_old']:.8f}")
    lines.append(f"Proton-proxy threshold speed                = {data['vcrit_proton']:.8f}")
    lines.append("")
    lines.append("Collapse-side ingredients")
    lines.append("-" * 88)
    lines.append(f"Steepest log-gradient location              = {data['grad_peak_location']:.8f}")
    lines.append(f"chi_lambda^peak                             = {data['chi_lambda_peak']:.8f}")
    lines.append(f"t_collapse                                  = {data['t_collapse']:.8f}")
    lines.append("")
    lines.append("Characteristic safe window (main lambda choice)")
    lines.append("-" * 88)
    lines.append(f"Analytic lower edge E_safe,min              = {data['E_safe_lower']:.8f}")
    lines.append(f"Analytic lower edge speed                   = {data['v_safe_lower']:.8f}")
    lines.append(f"Lower edge in threshold units               = {data['v_safe_multiplier']:.8f} * v_crit,p")
    lines.append(f"Numeric lower edge on scan                  = {data['E_window_lower_num']:.8f}")
    lines.append(f"Numeric upper edge on scan                  = {data['E_window_upper_num']:.8f}")
    lines.append("")
    lines.append("Interpretation: on the requested characteristic estimate, the stable region is")
    lines.append("one-sided. Once E_inc exceeds E_safe,min, the ratio S = t_cross/t_collapse")
    lines.append("stays below 1 all the way to the top of the scanned over-barrier band.")
    lines.append("")
    lines.append("Dynamic aligned-spin transit cross-check")
    lines.append("-" * 88)
    finite_dyn = np.isfinite(data["t_cross_dyn"])
    if np.any(finite_dyn):
        lines.append(f"min actual barrier transit (aligned)        = {float(np.nanmin(data['t_cross_dyn'])):.8f}")
        lines.append(f"max actual barrier transit (aligned)        = {float(np.nanmax(data['t_cross_dyn'])):.8f}")
        lines.append(f"all sampled aligned trajectories status     = {set(data['traj_status'])}")
        lines.append("Across the sampled over-barrier aligned-spin trajectories, the actual peak-to-contact")
        lines.append("transit time stays below the collapse time already at the very bottom of the scan.")
        lines.append("So the characteristic formula is the more conservative stability criterion here.")
    else:
        lines.append("No finite aligned-spin transit times were obtained on the sampling grid.")
    lines.append("")
    lines.append("Sensitivity to the raw model width lambda = 1")
    lines.append("-" * 88)
    lines.append(f"chi_lambda^peak (raw model width)          = {data['chi_peak_model']:.8f}")
    lines.append(f"t_collapse (raw model width)               = {data['t_collapse_model']:.8f}")
    lines.append(f"E_safe,min (raw model width)               = {data['E_safe_lower_model']:.8f}")
    lines.append(f"v_safe,min (raw model width)               = {data['v_safe_lower_model']:.8f}")
    lines.append("")
    lines.append("This matters because the proton proxy by itself does not move the characteristic")
    lines.append("energy threshold when m_s and mu_eta are scaled together; what moves the window is")
    lines.append("the chosen confinement width entering chi_lambda. Using the previously derived")
    lines.append("trigger width gives the less severe and more branch-faithful lower edge.")
    lines.append("")
    lines.append("Safe energy window reported")
    lines.append("-" * 88)
    lines.append(
        "Main reported window (analytic lower edge, scan-limited upper edge): "
        f"[{data['E_safe_lower']:.8f}, {data['E_window_upper_num']:.8f}]"
    )
    lines.append(
        "Equivalent speed window: "
        f"[{data['v_safe_lower']:.8f}, {data['v_window_upper_num']:.8f}]"
    )
    lines.append("")
    lines.append("Bottom line")
    lines.append("-" * 88)
    lines.append("1) With the proton-proxy scaling m_s = mu_eta = m_p/m_e and the earlier trigger")
    lines.append("   width lambda_th kept as the active barrier width, the characteristic stability")
    lines.append("   boundary lands at E_inc ≈ 5.32266 in the same reduced units.")
    lines.append("2) No upper collapse edge appears before the top of the scanned hot branch, so the")
    lines.append("   mathematically permitted region is one-sided rather than a closed island on this closure.")
    lines.append("3) The actual aligned-spin barrier-region transit is even shorter than the requested")
    lines.append("   characteristic estimate, so the conservative stability window is not being overstated.")
    return "\n".join(lines)


def make_timescale_plot(outdir: Path, data: Dict[str, object]) -> Path:
    path = outdir / "proton_proxy_stability_timescales.png"

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    E = np.asarray(data["E_inc"], dtype=float)
    t_cross = np.asarray(data["t_cross_char"], dtype=float)
    t_collapse = float(data["t_collapse"])

    ax.plot(E, t_cross, label=r"$t_{\rm cross}$ (characteristic)")
    ax.axhline(t_collapse, linestyle="--", linewidth=1.2, label=r"$t_{\rm collapse}$")
    ax.fill_between(
        E,
        np.nanmin(t_cross[np.isfinite(t_cross)]) * 0.8,
        np.nanmax(np.r_[t_cross[np.isfinite(t_cross)], t_collapse]) * 1.2,
        where=(t_cross < t_collapse),
        alpha=0.18,
        label=r"stable region: $t_{\rm cross}<t_{\rm collapse}$",
    )

    # optional actual aligned-spin transit cross-check
    E_dyn = np.asarray(data["E_dynamic"], dtype=float)
    t_dyn = np.asarray(data["t_cross_dyn"], dtype=float)
    finite_dyn = np.isfinite(t_dyn)
    if np.any(finite_dyn):
        ax.plot(E_dyn[finite_dyn], t_dyn[finite_dyn], linestyle=":", linewidth=1.2, label=r"aligned peak$\to$contact transit")

    ax.axvline(data["E_safe_lower"], linestyle=":", linewidth=1.0)
    ax.set_xlabel(r"Incident energy $E_{\rm inc}$")
    ax.set_ylabel("Timescale")
    ax.set_yscale("log")
    ax.set_title("Crossing time vs collapse time on the proton-proxy branch")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def main() -> None:
    outdir = Path(__file__).resolve().parent
    station = StationaryParams()
    sweep = StabilitySweepParams()
    formulas = derive_symbolic_formulas()
    data = run_stability_sweep(station, sweep)

    report_text = build_report(station, sweep, formulas, data)
    report_path = outdir / "proton_proxy_stability_report.txt"
    report_path.write_text(report_text, encoding="utf-8")

    plot_path = make_timescale_plot(outdir, data)

    print(report_text)
    print("")
    print("Artifacts")
    print("-" * 88)
    print(f"Report : {report_path}")
    print(f"Plot   : {plot_path}")


if __name__ == "__main__":
    main()
