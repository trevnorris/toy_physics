#!/usr/bin/env python3
"""
Reduced SymPy stress-test for the 4+1D moving-throat same-charge barrier model.

What this script does
---------------------
It does NOT claim to solve the full moving-throat PDE branch. Instead, it builds a
doc-grounded reduced stress-test closure that relaxes three constraints:

1) Transverse leakage: J^w != 0, with explicit S_leak^(s) and J^w E_w.
2) Rigid-mouth factorization: U/V transfer-shape and dressing legs are cross-coupled
   by a chi_lambda-scaled quadratic free-energy term.
3) Sign-changing multimode mouth sources: sigma(z) is allowed to oscillate and change sign,
   so the positive-source Family-1 closure is deliberately relaxed.

The code keeps the exact reduced bundle formulas wherever the file stack already fixes them,
and only uses simple phenomenological choices where the full moving-throat branch is still open.

Outputs
-------
- A text report with exact symbolic formulas and numerical diagnostics.
- A plot comparing Coulomb 1/r with the relaxed reduced V_eff(r).
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Symbolic derivations
# -----------------------------------------------------------------------------

def derive_symbolic_formulas() -> Dict[str, sp.Expr]:
    """Return the exact symbolic formulas used by the reduced stress-test."""
    formulas: Dict[str, sp.Expr] = {}

    # -------------------------------------------------------------------------
    # Task 1: leakage source and scalar-photon work channel
    # Choose a normalized Gaussian projection kernel and a minimal E_w-driven
    # transverse-current closure:
    #
    #   W(w)  = exp(-w^2/lambda^2)/(sqrt(pi) lambda)
    #   E_w   = E0(r) * W'(w)/W(w)
    #   j_w   = rho0 * mu_w * E0(r) * W'(w)
    #
    # Then:
    #   S_leak = \int W'(w) j_w dw
    #   Work   = q \int j_w E_w dw
    #
    # Both are exact.
    # -------------------------------------------------------------------------
    w, lam, E0, rho0, muw, q = sp.symbols("w lambda E0 rho0 mu_w q", positive=True, real=True)
    W = sp.exp(-w**2 / lam**2) / (sp.sqrt(sp.pi) * lam)
    Wp = sp.diff(W, w)
    Ew = sp.simplify(E0 * Wp / W)
    jw = sp.simplify(rho0 * muw * E0 * Wp)
    S_leak = sp.simplify(sp.integrate(Wp * jw, (w, -sp.oo, sp.oo)))
    Work_bulk = sp.simplify(q * sp.integrate(jw * Ew, (w, -sp.oo, sp.oo)))

    formulas["W(w)"] = W
    formulas["W'(w)"] = Wp
    formulas["E_w(w)"] = Ew
    formulas["j_w(w)"] = jw
    formulas["S_leak_exact"] = S_leak
    formulas["Bulk_work_exact"] = Work_bulk

    # -------------------------------------------------------------------------
    # Task 2: rigid-mouth U/V chart with chi_lambda cross-coupling
    #
    #   U = ln(T^2 / T_ref^2)
    #   V = ln(epsilon_eta / epsilon_eta_ref)
    #
    # Use a reduced quadratic free energy
    #
    #   F(U,V;r) = aU/2 U^2 + aV/2 V^2 + gUV * chi(r) * U V - s(r) U
    #
    # with chi(r) = lambda / (r + r_reg) and s(r) = s0 exp(-r/lambda).
    #
    # Stationarity gives exact U(r), V(r). The first-order packet scalar Xi_1 is
    # identified with U on the rigid-mouth chart. A finite nonlinear proxy is also
    # reported:
    #
    #   Xi_tilde = - ln(R_target/R_ref) - c_eta * ln(epsilon_eta/epsilon_eta_ref)
    #
    # with c_eta = epsilon_ref / (1 - epsilon_ref).
    # -------------------------------------------------------------------------
    r, rreg, aU, aV, gUV, s0 = sp.symbols("r r_reg a_U a_V g_UV s0", positive=True, real=True)
    eps_ref = sp.symbols("epsilon_ref", positive=True, real=True)

    chi = lam / (r + rreg)
    s = s0 * sp.exp(-r / lam)

    M = sp.Matrix([[aU, gUV * chi], [gUV * chi, aV]])
    rhs = sp.Matrix([s, 0])
    U_expr, V_expr = [sp.simplify(x) for x in M.LUsolve(rhs)]

    F = aU * sp.Symbol("U")**2 / 2 + aV * sp.Symbol("V")**2 / 2 + gUV * chi * sp.Symbol("U") * sp.Symbol("V") - s * sp.Symbol("U")
    # exact energy reduction relative to uncoupled branch
    coupling_energy_drop = sp.simplify((chi**2 * gUV**2 * s**2) / (2 * aU * (aU * aV - chi**2 * gUV**2)))

    c_eta = sp.simplify(eps_ref / (1 - eps_ref))
    R_ratio = sp.simplify(((1 - eps_ref * sp.exp(V_expr)) / (1 - eps_ref)) * sp.exp(-U_expr))
    Xi_tilde = sp.simplify(-sp.log(R_ratio) - c_eta * V_expr)

    formulas["chi_lambda(r)"] = chi
    formulas["s(r)"] = s
    formulas["U_exact(r)"] = U_expr
    formulas["V_exact(r)"] = V_expr
    formulas["epsilon_eta(r)"] = sp.simplify(eps_ref * sp.exp(V_expr))
    formulas["Xi_1_first_order(r)"] = U_expr
    formulas["Xi_tilde_finite(r)"] = Xi_tilde
    formulas["UV_coupling_energy_drop"] = coupling_energy_drop
    formulas["R_target_ratio(U,V)"] = sp.simplify(((1 - eps_ref * sp.exp(sp.Symbol("V"))) / (1 - eps_ref)) * sp.exp(-sp.Symbol("U")))

    # -------------------------------------------------------------------------
    # Task 3: sign-changing multimode source branch
    #
    # Put x = z/L in [0,1] and use
    #
    #   sigma(x) = 1 + a cos(pi x) + b cos(2 pi x)
    #
    # which is automatically normalized because the cosine modes have zero mean.
    #
    # The exact mouth-bias factor is
    #   g[sigma] = \int_0^1 sigma(x) cos(pi x/2) dx
    #
    # and the formal Family-1 mixed-to-shell ratio is
    #   R[sigma] = (g[sigma] - r_F1)^2 / (1 + r_F1^2)
    #
    # even though sign-changing sigma goes beyond the positive-source theorem.
    # -------------------------------------------------------------------------
    x, a, b, rF1 = sp.symbols("x a b r_F1", real=True)
    sigma = 1 + a * sp.cos(sp.pi * x) + b * sp.cos(2 * sp.pi * x)
    g_sigma = sp.simplify(sp.integrate(sigma * sp.cos(sp.pi * x / 2), (x, 0, 1)))
    R_sigma = sp.simplify((g_sigma - rF1) ** 2 / (1 + rF1 ** 2))

    n = sp.symbols("n", integer=True, positive=True)
    cosine_mode_projection = sp.simplify(sp.integrate(sp.cos(sp.pi * x / 2) * sp.cos(sp.pi * n * x), (x, 0, 1)))

    formulas["sigma(x)"] = sigma
    formulas["g_sigma_exact"] = g_sigma
    formulas["R_sigma_exact"] = R_sigma
    formulas["cos_mode_projection_exact"] = cosine_mode_projection

    # -------------------------------------------------------------------------
    # Optional dynamic theorem: linear outgoing phase-lag no-go
    #
    # If the outgoing port is passive:
    #   Pi(omega) = i Gamma_out(omega)
    #
    # and T_J(omega)^2 is real on the conservative branch, then
    #   delta V_mix = -(i/2) Gamma_out T_J^2
    #
    # so:
    #   Re delta V_mix = 0
    #   P_abs = (omega/2) Gamma_out T_J^2 >= 0
    # -------------------------------------------------------------------------
    omega, Gamma_out, T_J = sp.symbols("omega Gamma_out T_J", positive=True, real=True)
    deltaV_dyn = -sp.I * sp.Rational(1, 2) * Gamma_out * T_J**2
    P_abs = sp.simplify(-omega * sp.im(deltaV_dyn))

    formulas["deltaV_dynamic_first_order"] = deltaV_dyn
    formulas["P_abs_first_order"] = P_abs

    return formulas


# -----------------------------------------------------------------------------
# Numeric reduced closure
# -----------------------------------------------------------------------------

@dataclass
class ModelParams:
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
    kk_amp: float = 0.0  # set >0 to restore the Gaussian KK repulsive correction

    # reference Family-1 value
    r_F1: float = 1.77799353547498


def compute_uv_branch(r: np.ndarray, p: ModelParams) -> Dict[str, np.ndarray]:
    chi = p.lam / (r + p.r_reg)
    s = p.s0 * np.exp(-r / p.lam)
    den = p.a_U * p.a_V - p.g_UV**2 * chi**2

    # branch stability check
    if np.any(den <= 0):
        raise ValueError("U/V branch denominator crossed zero; choose a weaker coupling or larger r_reg.")

    U = p.a_V * s / den
    V = -p.g_UV * chi * s / den

    c_eta = p.epsilon_ref / (1.0 - p.epsilon_ref)
    Xi_1 = U
    Xi_tilde = U - np.log((1.0 - p.epsilon_ref * np.exp(V)) / (1.0 - p.epsilon_ref)) - c_eta * V
    epsilon_eta = p.epsilon_ref * np.exp(V)

    energy_drop = (chi**2 * p.g_UV**2 * s**2) / (2.0 * p.a_U * (p.a_U * p.a_V - p.g_UV**2 * chi**2))

    return {
        "chi_lambda": chi,
        "s": s,
        "U": U,
        "V": V,
        "Xi_1": Xi_1,
        "Xi_tilde": Xi_tilde,
        "epsilon_eta": epsilon_eta,
        "UV_energy_drop": energy_drop,
    }


def compute_leakage_diagnostics(r: np.ndarray, p: ModelParams) -> Dict[str, np.ndarray]:
    # minimal radial envelope for the transverse electric driver
    E0 = p.E0_amp / (r + p.r_reg) ** 2
    S_leak = p.rho0 * p.mu_w * E0 / (np.sqrt(2.0 * np.pi) * p.lam**3)
    work_bulk = 2.0 * p.q * p.rho0 * p.mu_w * E0**2 / p.lam**2

    return {
        "E0": E0,
        "S_leak": S_leak,
        "Work_bulk": work_bulk,
    }


def compute_sigma_branch(r: np.ndarray, p: ModelParams) -> Dict[str, np.ndarray]:
    a_amp = p.a0 / (1.0 + (r / p.r_sigma) ** 2)
    b_amp = p.b0 / (1.0 + (r / p.r_sigma) ** 2)

    g_sigma = 2.0 / np.pi + 2.0 * a_amp / (3.0 * np.pi) - 2.0 * b_amp / (15.0 * np.pi)
    R_sigma = ((g_sigma - p.r_F1) ** 2) / (1.0 + p.r_F1**2)

    # detect sign change numerically on x in [0,1]
    xgrid = np.linspace(0.0, 1.0, 400)
    sigma_min = []
    sigma_max = []
    sign_change = []
    for ai, bi in zip(a_amp, b_amp):
        sig = 1.0 + ai * np.cos(np.pi * xgrid) + bi * np.cos(2.0 * np.pi * xgrid)
        sigma_min.append(sig.min())
        sigma_max.append(sig.max())
        sign_change.append(sig.min() < 0.0)

    return {
        "a_amp": a_amp,
        "b_amp": b_amp,
        "g_sigma": g_sigma,
        "R_sigma": R_sigma,
        "sigma_min": np.array(sigma_min),
        "sigma_max": np.array(sigma_max),
        "sign_change": np.array(sign_change, dtype=bool),
    }


def one_port_susceptibilities(p: ModelParams) -> Dict[str, float]:
    Delta = p.Omega_U_sq * p.Omega_W_sq - p.R_mix**2
    Q = p.G_U**2 * p.Omega_W_sq + 2.0 * p.G_U * p.G_W * p.R_mix + p.G_W**2 * p.Omega_U_sq
    P = p.Omega_U_sq * p.G_W + p.R_mix * p.G_U
    P_U = p.G_U * p.Omega_W_sq + p.R_mix * p.G_W
    D0 = p.K_star - Q / Delta

    if Delta <= 0 or D0 <= 0:
        raise ValueError("One-port bundle left the admissible branch: need Delta>0 and D0>0.")

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
    p: ModelParams,
    sigma_branch: Dict[str, np.ndarray],
    sus: Dict[str, float],
) -> Dict[str, np.ndarray]:
    # minimal stress-test feedback from the sign-changing source into the primitive
    # Yukawa-sector amplitudes: this is the single phenomenological choice not fixed
    # by the file stack, and it is kept explicit on purpose.
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


def compute_effective_potential(r: np.ndarray, p: ModelParams) -> Dict[str, np.ndarray]:
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


# -----------------------------------------------------------------------------
# Reporting / plotting
# -----------------------------------------------------------------------------

def build_text_report(results: Dict[str, np.ndarray], p: ModelParams, formulas: Dict[str, sp.Expr]) -> str:
    r = results["r"]
    i_eval = int(np.argmin(np.abs(r - p.lam)))  # evaluate at r = lambda by default
    i_soft = int(np.argmin(results["V_eff"] / results["V_coulomb"]))

    lines = []
    lines.append("Reduced moving-throat same-charge stress-test")
    lines.append("=" * 72)
    lines.append("")
    lines.append("This report is a reduced closure, not a solved full moving-throat PDE theorem.")
    lines.append("It keeps the exact reduced formulas fixed by the file stack and makes only one")
    lines.append("explicit phenomenological choice: the sign-changing source branch feeds back into")
    lines.append("the primitive Yukawa amplitudes beta_U, beta_W through a linear function of R[sigma].")
    lines.append("")
    lines.append("Chosen parameter set")
    lines.append("-" * 72)
    for k, v in asdict(p).items():
        lines.append(f"{k:>14s} = {v}")
    lines.append("")
    lines.append("Exact symbolic formulas")
    lines.append("-" * 72)
    lines.append(f"S_leak^(s) = {sp.simplify(formulas['S_leak_exact'])}")
    lines.append(f"J^w E_w work = {sp.simplify(formulas['Bulk_work_exact'])}")
    lines.append("")
    lines.append(f"U(r) = {sp.simplify(formulas['U_exact(r)'])}")
    lines.append(f"V(r) = {sp.simplify(formulas['V_exact(r)'])}")
    lines.append(f"Xi_1(r) [first-order rigid-mouth packet] = {sp.simplify(formulas['Xi_1_first_order(r)'])}")
    lines.append(f"Xi_tilde(r) [finite nonlinear rigid-mouth proxy] = {sp.simplify(formulas['Xi_tilde_finite(r)'])}")
    lines.append(f"UV energy-drop = {sp.simplify(formulas['UV_coupling_energy_drop'])}")
    lines.append("")
    lines.append(f"g[sigma] for sigma=1+a cos(pi x)+b cos(2 pi x) = {sp.simplify(formulas['g_sigma_exact'])}")
    lines.append(f"R[sigma] = {sp.simplify(formulas['R_sigma_exact'])}")
    lines.append("")
    lines.append(f"Dynamic first-order outgoing correction = {sp.simplify(formulas['deltaV_dynamic_first_order'])}")
    lines.append(f"Average absorbed power = {sp.simplify(formulas['P_abs_first_order'])}")
    lines.append("")
    lines.append("One-port admissibility")
    lines.append("-" * 72)
    lines.append(f"Delta = {results['Delta']:.8f}  (must be >0)")
    lines.append(f"D0    = {results['D0']:.8f}  (must be >0)")
    lines.append("")
    lines.append("Evaluation at r = lambda")
    lines.append("-" * 72)
    lines.append(f"r_eval           = {r[i_eval]:.8f}")
    lines.append(f"Xi_1(r_eval)     = {results['Xi_1'][i_eval]:.8f}")
    lines.append(f"Xi_tilde(r_eval) = {results['Xi_tilde'][i_eval]:.8f}")
    lines.append(f"U(r_eval)        = {results['U'][i_eval]:.8f}")
    lines.append(f"V(r_eval)        = {results['V'][i_eval]:.8f}")
    lines.append(f"epsilon_eta      = {results['epsilon_eta'][i_eval]:.8f}")
    lines.append(f"S_leak(r_eval)   = {results['S_leak'][i_eval]:.8f}")
    lines.append(f"Jw*Ew(r_eval)    = {results['Work_bulk'][i_eval]:.8f}")
    lines.append(f"g_sigma(r_eval)  = {results['g_sigma'][i_eval]:.8f}")
    lines.append(f"R_sigma(r_eval)  = {results['R_sigma'][i_eval]:.8f}")
    lines.append(f"sigma_min        = {results['sigma_min'][i_eval]:.8f}")
    lines.append(f"sign_change?     = {bool(results['sign_change'][i_eval])}")
    lines.append("")
    lines.append("Strongest barrier softening in plotted range")
    lines.append("-" * 72)
    lines.append(f"r_soft          = {r[i_soft]:.8f}")
    lines.append(f"V_eff           = {results['V_eff'][i_soft]:.8f}")
    lines.append(f"V_coulomb       = {results['V_coulomb'][i_soft]:.8f}")
    lines.append(f"V_eff / V_coul  = {results['V_eff'][i_soft] / results['V_coulomb'][i_soft]:.8f}")
    lines.append(f"UV drop         = {results['UV_energy_drop'][i_soft]:.8f}")
    lines.append(f"bulk work       = {results['Work_bulk'][i_soft]:.8f}")
    lines.append(f"g_sigma         = {results['g_sigma'][i_soft]:.8f}")
    lines.append(f"R_sigma         = {results['R_sigma'][i_soft]:.8f}")
    lines.append("")
    lines.append("Interpretation")
    lines.append("-" * 72)
    lines.append("1) Leakage and J^w E_w both rise as r decreases because E0(r) was chosen to sharpen")
    lines.append("   at small separation. That is the open-system export channel.")
    lines.append("2) The UV cross-coupling makes V nonzero, so epsilon_eta is not invariant once")
    lines.append("   chi_lambda grows. The exact energy-drop formula quantifies the amount of")
    lines.append("   transfer-shape energy drained into the dressing leg.")
    lines.append("3) The sign-changing sigma branch can push g[sigma] outside the positive-source range.")
    lines.append("   That is the explicit compensated-branch stress test requested in Task 3.")
    lines.append("4) The plotted V_eff(r) is lower than Coulomb over part of the near-contact range")
    lines.append("   for this parameter set, but the reduction is still carried only by short-range")
    lines.append("   families plus open-system export; no new long-range attractive law is inserted.")
    return "\n".join(lines)


def make_plots(results: Dict[str, np.ndarray], outdir: Path) -> Tuple[Path, Path]:
    r = results["r"]

    potential_path = outdir / "relaxed_constraints_veff.png"
    diagnostics_path = outdir / "relaxed_constraints_diagnostics.png"

    # Potential comparison
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(r, results["V_coulomb"], label="Coulomb 1/r")
    ax.plot(r, results["V_eff"], label="Relaxed reduced V_eff(r)")
    ax.set_xlabel("Separation r")
    ax.set_ylabel("Effective potential")
    ax.set_title("Same-charge barrier stress test")
    ax.set_xlim(r.min(), r.max())
    ax.set_ylim(min(results["V_eff"].min(), results["V_coulomb"].min()) - 0.1,
                max(results["V_eff"][0], results["V_coulomb"][0]) + 0.5)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(potential_path, dpi=160)
    plt.close(fig)

    # Diagnostics
    fig, ax = plt.subplots(figsize=(8, 5))
    # normalize diagnostics so they share one axis
    work_n = results["Work_bulk"] / np.max(results["Work_bulk"])
    uv_n = results["UV_energy_drop"] / np.max(results["UV_energy_drop"])
    Xi_n = results["Xi_1"] / np.max(results["Xi_1"])
    R_n = results["R_sigma"] / np.max(results["R_sigma"])

    ax.plot(r, Xi_n, label="Xi_1 / max")
    ax.plot(r, uv_n, label="UV drain / max")
    ax.plot(r, work_n, label="J^w E_w / max")
    ax.plot(r, R_n, label="R[sigma] / max")
    ax.set_xlabel("Separation r")
    ax.set_ylabel("Normalized diagnostic")
    ax.set_title("Open-system and mouth-response diagnostics")
    ax.set_xlim(r.min(), r.max())
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(diagnostics_path, dpi=160)
    plt.close(fig)

    return potential_path, diagnostics_path


def main() -> None:
    outdir = Path("/mnt/data")
    formulas = derive_symbolic_formulas()

    p = ModelParams()
    r = np.linspace(0.18, 4.5, 600)
    results = compute_effective_potential(r, p)

    report_text = build_text_report(results, p, formulas)
    report_path = outdir / "stress_test_relaxed_constraints_report.txt"
    report_path.write_text(report_text)

    potential_path, diagnostics_path = make_plots(results, outdir)

    print(f"Wrote report: {report_path}")
    print(f"Wrote potential plot: {potential_path}")
    print(f"Wrote diagnostics plot: {diagnostics_path}")
    print("")
    print("Exact first-order Xi_1(r) formula:")
    print(sp.simplify(formulas["Xi_1_first_order(r)"]))
    print("")
    print("At r = lambda:")
    i_eval = int(np.argmin(np.abs(results['r'] - p.lam)))
    print(f"Xi_1(lambda) = {results['Xi_1'][i_eval]:.8f}")
    print(f"Xi_tilde(lambda) = {results['Xi_tilde'][i_eval]:.8f}")
    print(f"S_leak(lambda) = {results['S_leak'][i_eval]:.8f}")
    print(f"J^w E_w(lambda) = {results['Work_bulk'][i_eval]:.8f}")
    print(f"g_sigma(lambda) = {results['g_sigma'][i_eval]:.8f}")
    print(f"R_sigma(lambda) = {results['R_sigma'][i_eval]:.8f}")


if __name__ == "__main__":
    main()
