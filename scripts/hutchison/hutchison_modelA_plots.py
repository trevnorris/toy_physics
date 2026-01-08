#!/usr/bin/env python3
"""
Model A: Hutchison-module plots (analytic phase maps + dynamic cross-check + coupled capture proxy)

This script reproduces the figures from the latest modeling step:
  1) Analytic phase maps in (frequency, pA) for L = {5,10,20} mm
  2) Cross-check: dynamic envelope ODE vs analytic rolloff (L=10 mm)
  3) Coupled COM + internal-mode capture proxy vs amplitude (with damping feedback)

Notes:
- Thought-experiment model only.
- Units:
    pA in Pa, L in m, frequency in Hz (plots show kHz), g in m/s^2
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------------
# Baseline knobs (edit these)
# ----------------------------
g = 9.80665
rho_m = 2700.0     # aluminum density (kg/m^3)
rho_star = 1.2     # effective coupled density knob (kg/m^3)
c_star = 578.0     # effective coupled wave speed (m/s)
chi = 1.0          # shape/mode factor (O(1))
Q = 300.0          # internal-mode Q
epscrit = 1e-4     # critical strain proxy
alpha_over_beta = 1.0  # (α/β), geometry/overlap ratio

# Amplitude grid for maps (Pa)
pA_grid = np.linspace(0.0, 20_000.0, 201)  # 0..20 kPa


def f0(L: float) -> float:
    """Fundamental collective mode frequency (Hz)."""
    return chi * c_star / (2.0 * L)


def k_lock(L: float) -> float:
    """Standing-wave lock to size: k ≈ pi/L (1/m)."""
    return math.pi / L


def lift_capability(pA: float, L: float) -> float:
    """
    Ponderomotive capability metric: a0/g where
    a0 = pA^2 k / (4 rho_m rho_star c_star^2).
    """
    k = k_lock(L)
    a0 = (pA * pA * k) / (4.0 * rho_m * rho_star * c_star * c_star)
    return a0 / g


def eps_res_onres(pA: float) -> float:
    """
    On-resonance strain proxy:
      eps_res ≈ pA*Q / ((α/β) π^2 ρ_m c_*^2)
    """
    return (pA * Q) / (alpha_over_beta * math.pi**2 * rho_m * c_star**2)


def eps_rms_analytic(pA: float, f: float, L: float) -> float:
    """
    Analytic Lorentzian rolloff:
      eps(f) = eps_res_onres(pA) / sqrt(1 + (2Q delta)^2)
      delta = (f - f0)/f0
    """
    d = (f - f0(L)) / f0(L)
    roll = 1.0 / math.sqrt(1.0 + (2.0 * Q * d) ** 2)
    return eps_res_onres(pA) * roll


def p_jelly_onres() -> float:
    """
    Threshold pA at exact resonance for epscrit:
      pA = (α/β) π^2 ρ_m c_*^2 epscrit / Q
    """
    return alpha_over_beta * (math.pi**2) * rho_m * (c_star**2) * (epscrit / Q)


def p_jelly_of_f(f: float, L: float) -> float:
    """
    Threshold curve pA(f) such that eps_rms_analytic(pA,f,L)=epscrit.
    """
    d = (f - f0(L)) / f0(L)
    roll = 1.0 / math.sqrt(1.0 + (2.0 * Q * d) ** 2)
    return p_jelly_onres() / roll


def p_lev(L: float) -> float:
    """
    Solve lift_capability(pA,L)=1 for pA, under k=pi/L:
      pA_lev = 2 c_* sqrt(rho_m rho_* g / k)
    """
    k = k_lock(L)
    return 2.0 * c_star * math.sqrt(rho_m * rho_star * g / k)


def plot_phase_maps(Ls_mm=(5, 10, 20)) -> pd.DataFrame:
    """Make analytic phase maps for each L and return a summary table."""
    phase_rows = []
    for Lmm in Ls_mm:
        L = Lmm / 1000.0
        f0_Hz = f0(L)
        df_Hz = f0_Hz / Q
        f_grid = np.linspace(f0_Hz - 6 * df_Hz, f0_Hz + 6 * df_Hz, 401)

        pJ = np.array([p_jelly_of_f(fi, L) for fi in f_grid])  # Pa
        pL = p_lev(L)  # Pa

        # Region classification: 0 none, 1 jelly, 2 lift, 3 both
        PJ = pJ[None, :]
        PA = pA_grid[:, None]
        jelly = (PA >= PJ)
        lift = (PA >= pL)
        phase = jelly.astype(int) + 2 * lift.astype(int)

        plt.figure(figsize=(7.4, 4.8))
        plt.imshow(
            phase,
            aspect="auto",
            origin="lower",
            extent=[f_grid[0] / 1000, f_grid[-1] / 1000, pA_grid[0] / 1000, pA_grid[-1] / 1000],
            interpolation="nearest",
        )
        plt.colorbar(label="Phase: 0 none, 1 jelly, 2 lift, 3 both")
        plt.plot(f_grid / 1000, pJ / 1000, linewidth=1.5, label="jelly threshold")
        plt.axhline(pL / 1000, linestyle="--", linewidth=1.5, label="lift threshold")
        plt.axvline(f0_Hz / 1000, linestyle=":", linewidth=1.5, label="f0")
        plt.xlabel("Drive frequency (kHz)")
        plt.ylabel("Pressure amplitude pA (kPa)")
        plt.title(f"Model A phase map (analytic boundaries) — L={Lmm} mm, Q={int(Q)}")
        plt.legend(loc="upper right")
        plt.tight_layout()
        plt.show()

        phase_rows.append(
            {
                "L (mm)": Lmm,
                "f0 (kHz)": f0_Hz / 1000,
                "Δf (Hz)": df_Hz,
                "p_jelly,on-res (kPa)": p_jelly_onres() / 1000,
                "p_lev (kPa)": pL / 1000,
                "At f0: jelly easier than lift?": (p_jelly_onres() < pL),
                "Scaling check 2Lf0 (m/s)": 2 * L * f0_Hz,
            }
        )
    return pd.DataFrame(phase_rows)


def internal_envelope_sim(pA: float, f_drive: float, L: float, t_end=0.15, dt=2e-5) -> float:
    """
    Simulate only the internal-mode envelope quadratures (u,v) with constant local drive
    at an antinode (cos(kz)=1). Returns mean eps_rms over final half of the run.
    """
    w0 = 2 * math.pi * f0(L)
    w = 2 * math.pi * f_drive
    gamma_int = w0 / (2 * Q)
    Delta = (w0 * w0 - w * w) / (2 * w)

    # drive per mass; assumes (β/α)=1 and S/V~1/L
    f_drv = pA / (rho_m * L)

    u = 0.0
    vq = 0.0
    n = int(t_end / dt)
    burn = int(0.5 * n)
    eps_sum = 0.0
    count = 0

    for i in range(n):
        # midpoint RK2 is fine for slow envelope
        du1 = -gamma_int * u + Delta * vq + f_drv / (2 * w)
        dv1 = -gamma_int * vq - Delta * u
        u_mid = u + 0.5 * dt * du1
        v_mid = vq + 0.5 * dt * dv1
        du2 = -gamma_int * u_mid + Delta * v_mid + f_drv / (2 * w)
        dv2 = -gamma_int * v_mid - Delta * u_mid
        u += dt * du2
        vq += dt * dv2

        if i >= burn:
            amp = math.sqrt(u * u + vq * vq)
            eps = (amp / 1.41421356237) / L
            eps_sum += eps
            count += 1

    return eps_sum / count


def plot_crosscheck(Lmm=10, pAs=(2000.0, 5000.0, 12000.0)) -> None:
    """Compare envelope ODE simulation to analytic rolloff for several amplitudes."""
    L = Lmm / 1000.0
    f0_Hz = f0(L)
    df_Hz = f0_Hz / Q
    f_sweep = np.linspace(f0_Hz - 6 * df_Hz, f0_Hz + 6 * df_Hz, 81)

    plt.figure(figsize=(7.4, 4.8))
    for pA in pAs:
        eps_sim = np.array([internal_envelope_sim(pA, f, L) for f in f_sweep])
        eps_an = np.array([eps_rms_analytic(pA, f, L) for f in f_sweep])
        plt.plot(f_sweep / 1000, eps_sim, label=f"sim pA={int(pA)} Pa")
        plt.plot(f_sweep / 1000, eps_an, linestyle="--", label=f"analytic pA={int(pA)} Pa")

    plt.axhline(epscrit, linestyle=":", label="epscrit")
    plt.xlabel("Drive frequency (kHz)")
    plt.ylabel("ε_rms (strain proxy)")
    plt.title("Cross-check: dynamic envelope ODE vs analytic rolloff (L=10 mm)")
    plt.legend(ncol=2, fontsize=9)
    plt.tight_layout()
    plt.show()


def coupled_com_internal_sim(
    pA: float,
    L: float,
    f_drive: float | None = None,
    gamma_z: float = 8.0,
    feedback_kappa: float = 80.0,
    t_end: float = 0.25,
    dt: float = 2e-5,
) -> tuple[float, float, float]:
    """
    Minimal coupled model: COM (z,v) + internal envelope (u,v).

    Returns:
      z_std (m) over final half of run,
      eps_mean over final half,
      liftcap (a0/g) (capability, not realized lift)
    """
    if f_drive is None:
        f_drive = f0(L)

    k = k_lock(L)
    w0 = 2 * math.pi * f0(L)
    w = 2 * math.pi * f_drive
    gamma_int = w0 / (2 * Q)
    Delta = (w0 * w0 - w * w) / (2 * w)

    a0_base = (pA * pA * k) / (4.0 * rho_m * rho_star * c_star * c_star)

    z = 0.0
    vz = 0.0
    u = 0.0
    vq = 0.0

    n = int(t_end / dt)
    burn = int(0.5 * n)
    z_sum = 0.0
    z2_sum = 0.0
    eps_sum = 0.0
    count = 0

    for i in range(n):
        amp = math.sqrt(u * u + vq * vq)
        eps = (amp / 1.41421356237) / L
        gz = gamma_z * (1.0 + feedback_kappa * (eps / epscrit) ** 2) if epscrit > 0 else gamma_z

        a_rad = a0_base * math.sin(2 * k * z)

        # COM (Euler)
        z += dt * vz
        vz += dt * (-g + a_rad - gz * vz)

        # internal drive depends on position via cos(k z)
        p0 = pA * math.cos(k * z)
        f_drv = p0 / (rho_m * L)
        du = -gamma_int * u + Delta * vq + f_drv / (2 * w)
        dv = -gamma_int * vq - Delta * u
        u += dt * du
        vq += dt * dv

        if i >= burn:
            z_sum += z
            z2_sum += z * z
            eps_sum += eps
            count += 1

    z_mean = z_sum / count
    z_std = math.sqrt(max(0.0, z2_sum / count - z_mean * z_mean))
    return z_std, eps_sum / count, a0_base / g


def plot_capture_maps(Ls_mm=(5, 10, 20), pA_max=20_000.0) -> None:
    """Plot capture proxy z_std vs amplitude and internal excitation vs amplitude (L=10mm)."""
    pA_cap = np.linspace(0.0, pA_max, 81)

    # compute
    cap_rows = []
    for Lmm in Ls_mm:
        L = Lmm / 1000.0
        for pA in pA_cap:
            z_std, eps_mean, liftcap = coupled_com_internal_sim(pA, L)
            cap_rows.append(
                {"L (mm)": Lmm, "pA (kPa)": pA / 1000, "z_std (mm)": z_std * 1000, "eps_rms_mean": eps_mean, "lift_cap (a0/g)": liftcap}
            )

    cap_df = pd.DataFrame(cap_rows)

    plt.figure(figsize=(7.4, 4.8))
    for Lmm in Ls_mm:
        sub = cap_df[cap_df["L (mm)"] == Lmm]
        plt.plot(sub["pA (kPa)"], sub["z_std (mm)"], label=f"L={Lmm} mm")
    plt.axhline(0.5, linestyle="--", label="tight trap (z_std=0.5 mm)")
    plt.xlabel("pA (kPa)")
    plt.ylabel("COM position spread z_std (mm)")
    plt.title("Coupled model: trapping/capture proxy vs amplitude (with damping feedback)")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Internal excitation vs amplitude for L=10mm
    sub10 = cap_df[cap_df["L (mm)"] == 10]
    plt.figure(figsize=(7.4, 4.8))
    plt.plot(sub10["pA (kPa)"], sub10["eps_rms_mean"], label="ε_rms_mean (L=10 mm)")
    plt.axhline(epscrit, linestyle="--", label="εcrit")
    plt.xlabel("pA (kPa)")
    plt.ylabel("ε_rms_mean")
    plt.title("Coupled model: internal excitation vs amplitude (L=10 mm, on resonance)")
    plt.legend()
    plt.tight_layout()
    plt.show()


def main() -> None:
    # Phase maps + summary table (printed)
    df = plot_phase_maps()
    print("\n=== Phase-map summary (analytic) ===")
    print(df.to_string(index=False))

    # Cross-check plot
    plot_crosscheck()

    # Coupled capture proxy plots
    plot_capture_maps()


if __name__ == "__main__":
    main()
