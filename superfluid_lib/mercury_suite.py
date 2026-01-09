#!/usr/bin/env python3
"""
mercury_suite.py

Publish-grade perihelion-precession benchmark harness for the toy 1PN / superfluid-analog model.

Outputs (out_dir):
- summary.json      (all run metadata + results)
- sweep.csv         (one row per configuration)
- omega_angles.csv  (for the single-run case; omitted during sweeps)
- trace.csv         (optional, --save_trace)

Modes:
- newton            : pure Newtonian gravity
- toy_r3            : a = -[mu/r^2 + 3 mu^2/(cs^2 r^3)] r_hat
- beta_analytic     : β-inertia correction, analytic Φ=-mu/r
- beta_grid_poisson : β-inertia correction, Φ and ∇Φ from FFT Poisson on a periodic box

Notes for publishable results:
- Always run Newtonian baseline to establish a numerical-precession floor.
- Run steps_per_orbit sweeps to show convergence.
- Run cs sweeps to verify ~1/cs^2 scaling (or your model’s predicted scaling).
"""

from __future__ import annotations

import argparse, json, math, os, time
from dataclasses import dataclass
from typing import Optional, Tuple, Any, Dict, List

import numpy as np


# -------------------------
# Backend selection (CPU/GPU)
# -------------------------
def get_xp(device: str):
    if device == "cpu":
        return np, False
    if device == "gpu":
        import cupy as cp  # type: ignore
        return cp, True
    # auto
    try:
        import cupy as cp  # type: ignore
        return cp, True
    except Exception:
        return np, False


# -------------------------
# Constants
# -------------------------
ARCSEC_PER_RAD = 206264.80624709636
AU_M = 149_597_870_700.0
C_SI = 299_792_458.0
MU_SUN_SI = 1.32712440018e20  # m^3 / s^2

# Mercury orbital period in days (benchmark constant for arcsec/century mapping)
MERCURY_PERIOD_DAYS = 87.9691
ORBITS_PER_CENTURY_MERCURY = (100.0 * 365.25) / MERCURY_PERIOD_DAYS


@dataclass
class OrbitParams:
    mu: float
    a: float
    e: float
    cs: float
    beta: float


def kepler_period(mu: float, a: float) -> float:
    return 2.0 * math.pi * math.sqrt(a**3 / mu)


def initial_state_at_perihelion(mu: float, a: float, e: float):
    rp = a * (1.0 - e)
    vp = math.sqrt(mu * (1.0 + e) / (a * (1.0 - e)))
    # Start at perihelion on +x axis, velocity +y (prograde)
    r0 = np.array([rp, 0.0, 0.0], dtype=np.float64)
    v0 = np.array([0.0, vp, 0.0], dtype=np.float64)
    return r0, v0


# -------------------------
# Kepler universal-variable drift
# -------------------------
def stumpff_C(xp, z):
    # stable C(z)
    eps = 1e-8
    z = xp.asarray(z)
    C = 0.5 - z / 24.0 + z * z / 720.0 - z * z * z / 40320.0
    C = xp.where(z > eps, (1.0 - xp.cos(xp.sqrt(z))) / z, C)
    C = xp.where(z < -eps, (xp.cosh(xp.sqrt(-z)) - 1.0) / (-z), C)
    return C


def stumpff_S(xp, z):
    # stable S(z)
    eps = 1e-8
    z = xp.asarray(z)
    S = 1.0 / 6.0 - z / 120.0 + z * z / 5040.0 - z * z * z / 362880.0
    S = xp.where(
        z > eps,
        (xp.sqrt(z) - xp.sin(xp.sqrt(z))) / (xp.sqrt(z) ** 3),
        S,
    )
    S = xp.where(
        z < -eps,
        (xp.sinh(xp.sqrt(-z)) - xp.sqrt(-z)) / (xp.sqrt(-z) ** 3),
        S,
    )
    return S


def kepler_drift_universal(xp, r0, v0, dt, mu, iters=8):
    """
    Propagate (r0,v0) forward by dt under pure Kepler (two-body) dynamics.
    r0, v0: shape (M,3)
    Returns r1, v1 with same shape.
    """
    r0n = xp.sqrt(xp.sum(r0 * r0, axis=1))
    v02 = xp.sum(v0 * v0, axis=1)
    vr0 = xp.sum(r0 * v0, axis=1) / r0n

    alpha = 2.0 / r0n - v02 / mu  # 1/a

    # initial guess for chi
    sqrt_mu = xp.sqrt(mu)
    chi = sqrt_mu * xp.abs(alpha) * dt
    # handle near-parabolic
    chi = xp.where(xp.abs(alpha) < 1e-12, sqrt_mu * dt / r0n, chi)

    for _ in range(iters):
        z = alpha * chi * chi
        C = stumpff_C(xp, z)
        S = stumpff_S(xp, z)

        F = (
            (r0n * vr0 / sqrt_mu) * chi * chi * C
            + (1.0 - alpha * r0n) * chi**3 * S
            + r0n * chi
            - sqrt_mu * dt
        )
        dF = (
            (r0n * vr0 / sqrt_mu) * chi * (1.0 - z * S)
            + (1.0 - alpha * r0n) * chi * chi * C
            + r0n
        )

        chi = chi - F / dF

    z = alpha * chi * chi
    C = stumpff_C(xp, z)
    S = stumpff_S(xp, z)

    f = 1.0 - (chi * chi / r0n) * C
    g = dt - (chi**3 / sqrt_mu) * S

    r1 = f[:, None] * r0 + g[:, None] * v0
    r1n = xp.sqrt(xp.sum(r1 * r1, axis=1))

    fdot = (sqrt_mu / (r1n * r0n)) * (z * S - 1.0) * chi
    gdot = 1.0 - (chi * chi / r1n) * C

    v1 = fdot[:, None] * r0 + gdot[:, None] * v0
    return r1, v1


# -------------------------
# Acceleration models
# -------------------------
def accel_newton(xp, pos, vel, mu: float):
    r2 = xp.sum(pos * pos, axis=1)
    r = xp.sqrt(r2)
    eps = 1e-30
    r_safe = xp.maximum(r, eps)
    r_hat = pos / r_safe[:, None]
    return (-mu / (r_safe**2))[:, None] * r_hat


def accel_toy_r3(xp, pos, vel, mu: float, cs: float):
    r2 = xp.sum(pos * pos, axis=1)
    r = xp.sqrt(r2)
    eps = 1e-30
    r_safe = xp.maximum(r, eps)
    r_hat = pos / r_safe[:, None]
    a_mag = -(mu / (r_safe**2) + 3.0 * (mu**2) / (cs**2 * r_safe**3))
    return a_mag[:, None] * r_hat


def accel_beta_analytic(xp, pos, vel, mu: float, cs: float, beta: float):
    """
    β-inertia correction with analytic Φ = -mu/r
      σ = -β Φ / cs^2 = β mu/(cs^2 r)
      ∇σ = -β mu/(cs^2 r^2) r_hat
      g = -mu/r^2 r_hat
      a = 1/(1+σ) * [ g - (v·∇σ) v + 1/2 v^2 ∇σ ]
    """
    r2 = xp.sum(pos * pos, axis=1)
    r = xp.sqrt(r2)
    eps = 1e-30
    r_safe = xp.maximum(r, eps)
    r_hat = pos / r_safe[:, None]

    g = (-mu / (r_safe**2))[:, None] * r_hat

    sigma = -(beta * mu)/(cs**2*r)
    grad_sigma = (+(beta*mu)/(cs**2*r**2)) * r_hat

    v2 = xp.sum(vel * vel, axis=1)
    vdot_gs = xp.sum(vel * grad_sigma, axis=1)

    inv = 1.0 / (1.0 + sigma)
    return inv[:, None] * (g - vdot_gs[:, None] * vel + 0.5 * v2[:, None] * grad_sigma)


# -------------------------
# Grid/PDE helpers (Poisson FFT)
# -------------------------
def build_grid(xp, N: int, L: float):
    dx = L / N
    lin = (xp.arange(N, dtype=xp.float64) - N / 2.0) * dx
    X, Y, Z = xp.meshgrid(lin, lin, lin, indexing="ij")
    return dx, X, Y, Z


def poisson_solve_fft(xp, rho, G: float, L: float):
    """
    Solve ∇^2 Φ = 4π G ρ on a periodic box using FFT.
    """
    N = rho.shape[0]
    dx = L / N

    kfreq = xp.fft.fftfreq(N, d=dx) * 2.0 * xp.pi
    KX, KY, KZ = xp.meshgrid(kfreq, kfreq, kfreq, indexing="ij")
    k2 = KX * KX + KY * KY + KZ * KZ

    inv_k2 = xp.zeros_like(k2)
    mask = k2 > 0
    inv_k2[mask] = -1.0 / k2[mask]

    rho_k = xp.fft.fftn(rho)
    phi_k = (4.0 * xp.pi * G) * rho_k * inv_k2
    return xp.fft.ifftn(phi_k).real


def trilinear_sample_periodic(xp, field, pos, L: float):
    """
    field: (N,N,N), pos: (M,3) in centered box coordinates (periodic)
    """
    N = field.shape[0]
    dx = L / N

    coords = (pos + (L / 2.0)) / dx  # [0,N)
    base = xp.floor(coords).astype(xp.int64)
    frac = coords - base

    x0 = base[:, 0] % N
    y0 = base[:, 1] % N
    z0 = base[:, 2] % N
    x1 = (x0 + 1) % N
    y1 = (y0 + 1) % N
    z1 = (z0 + 1) % N

    fx = frac[:, 0]
    fy = frac[:, 1]
    fz = frac[:, 2]

    c000 = field[x0, y0, z0]
    c001 = field[x0, y0, z1]
    c010 = field[x0, y1, z0]
    c011 = field[x0, y1, z1]
    c100 = field[x1, y0, z0]
    c101 = field[x1, y0, z1]
    c110 = field[x1, y1, z0]
    c111 = field[x1, y1, z1]

    c00 = c000 * (1 - fz) + c001 * fz
    c01 = c010 * (1 - fz) + c011 * fz
    c10 = c100 * (1 - fz) + c101 * fz
    c11 = c110 * (1 - fz) + c111 * fz

    c0 = c00 * (1 - fy) + c01 * fy
    c1 = c10 * (1 - fy) + c11 * fy

    return c0 * (1 - fx) + c1 * fx


def accel_beta_from_grid(xp, pos, vel, cs: float, beta: float,
                         phi_grid, grad_phi_grids: Tuple[Any, Any, Any], L: float):
    """
    Use Φ and ∇Φ sampled from a grid.
    σ = -β Φ/cs^2, ∇σ = -β ∇Φ/cs^2
    g = -∇Φ
    """
    gphix, gphiy, gphiz = grad_phi_grids

    phi_p = trilinear_sample_periodic(xp, phi_grid, pos, L)
    dphix = trilinear_sample_periodic(xp, gphix, pos, L)
    dphiy = trilinear_sample_periodic(xp, gphiy, pos, L)
    dphiz = trilinear_sample_periodic(xp, gphiz, pos, L)

    g = xp.stack((-dphix, -dphiy, -dphiz), axis=1)

    sigma = -beta * phi_p / (cs**2)
    grad_sigma = xp.stack((-beta * dphix / (cs**2),
                           -beta * dphiy / (cs**2),
                           -beta * dphiz / (cs**2)), axis=1)

    v2 = xp.sum(vel * vel, axis=1)
    vdot_gs = xp.sum(vel * grad_sigma, axis=1)

    inv = 1.0 / (1.0 + sigma)
    return inv[:, None] * (g - vdot_gs[:, None] * vel + 0.5 * v2[:, None] * grad_sigma)


# -------------------------
# Integration + perihelion measurement
# -------------------------
def integrate_and_measure(
    xp,
    mode: str,
    params: OrbitParams,
    n_orbits: int,
    steps_per_orbit: int,
    device_is_gpu: bool,
    grid_cfg: Optional[Dict[str, Any]] = None,
    save_trace: bool = False,
):
    mu, a, e, cs, beta = params.mu, params.a, params.e, params.cs, params.beta
    T = kepler_period(mu, a)
    dt = T / steps_per_orbit
    n_steps = int(n_orbits * steps_per_orbit)

    # initial state
    r0, v0 = initial_state_at_perihelion(mu, a, e)
    pos = xp.asarray(r0[None, :], dtype=xp.float64)
    vel = xp.asarray(v0[None, :], dtype=xp.float64)

    # optional PDE/FFT field precompute
    phi_grid = None
    grad_phi_grids = None
    boxL = None

    if mode == "beta_grid_poisson":
        assert grid_cfg is not None
        N = int(grid_cfg["gridN"])
        boxL = float(grid_cfg["boxL"])
        G = float(grid_cfg.get("G", 1.0))
        sun_sigma_cells = float(grid_cfg.get("sun_sigma_cells", 1.5))

        dx, X, Y, Z = build_grid(xp, N, boxL)
        R2 = X * X + Y * Y + Z * Z
        sigma = sun_sigma_cells * dx

        # total mass M = mu / G
        M = mu / G
        rho = xp.exp(-0.5 * R2 / (sigma * sigma))
        rho = rho * (M / (xp.sum(rho) * dx**3))

        phi_grid = poisson_solve_fft(xp, rho, G, boxL)
        gphix, gphiy, gphiz = xp.gradient(phi_grid, dx)
        grad_phi_grids = (gphix, gphiy, gphiz)

    # accel dispatch
    def accel(pos_, vel_):
        if mode == "newton":
            return accel_newton(xp, pos_, vel_, mu)
        if mode == "toy_r3":
            return accel_toy_r3(xp, pos_, vel_, mu, cs)
        if mode == "beta_analytic":
            return accel_beta_analytic(xp, pos_, vel_, mu, cs, beta)
        if mode == "beta_grid_poisson":
            return accel_beta_from_grid(xp, pos_, vel_, cs, beta, phi_grid, grad_phi_grids, boxL)
        raise ValueError(f"Unknown mode: {mode}")

    trace = []

    def maybe_trace(step, t, pos_, vel_):
        if not save_trace:
            return
        # sparse trace (~200 samples per orbit)
        if step % max(1, steps_per_orbit // 200) != 0:
            return
        p = np.asarray(pos_.get() if device_is_gpu else pos_)[0]
        v = np.asarray(vel_.get() if device_is_gpu else vel_)[0]
        trace.append({
            "t": float(t),
            "x": float(p[0]), "y": float(p[1]),
            "vx": float(v[0]), "vy": float(v[1]),
            "r": float(np.linalg.norm(p)),
        })

    # integrate: kick-drift-kick with exact Kepler drift
    t = 0.0

    # Measure periapsis direction via eccentricity vector angle ω = atan2(e_y, e_x)
    # Sample ω at a sparse cadence, unwrap + regress on CPU once at the end.
    sample_stride = max(1, steps_per_orbit // 200)
    n_samples = n_steps // sample_stride + 1
    omega_samples = xp.empty((n_samples,), dtype=xp.float64)
    t_samples = xp.empty((n_samples,), dtype=xp.float64)
    sample_count = 0

    start_time = time.time()
    log_interval = max(1, n_steps // 100)  # log ~20 times during run
    for step in range(n_steps):
        # total accel (depends on mode)
        a_tot = accel(pos, vel)

        # Newtonian accel
        a_newt = accel_newton(xp, pos, vel, mu)

        # perturbation accel (what we want to attribute to the toy model)
        a_pert = a_tot - a_newt

        # KICK (half step)
        vel = vel + 0.5 * dt * a_pert

        # DRIFT (exact Kepler step)
        pos, vel = kepler_drift_universal(xp, pos, vel, dt, mu)

        # KICK (half step) using updated state
        a_tot = accel(pos, vel)
        a_newt = accel_newton(xp, pos, vel, mu)
        a_pert = a_tot - a_newt
        vel = vel + 0.5 * dt * a_pert

        t += dt

        # sparse series capture for diagnostics (kept small)
        if step % sample_stride == 0:
            # Compute eccentricity vector e = (v x h)/mu - r_hat
            r_vec = pos[0]
            v_vec = vel[0]
            rnorm = xp.sqrt(xp.sum(r_vec * r_vec))
            rhat = r_vec / rnorm

            h = xp.cross(r_vec, v_vec)
            e_vec = xp.cross(v_vec, h) / mu - rhat

            omega = xp.arctan2(e_vec[1], e_vec[0])  # raw in (-pi, pi]
            t_samples[sample_count] = t
            omega_samples[sample_count] = omega
            sample_count += 1

        maybe_trace(step, t, pos, vel)

        # progress logging with ETA
        if step > 0 and step % log_interval == 0:
            elapsed = time.time() - start_time
            frac_done = step / n_steps
            eta_sec = elapsed / frac_done * (1 - frac_done)
            if elapsed >= 60:
                elapsed_str = f"{elapsed / 60:.1f} min"
            else:
                elapsed_str = f"{elapsed:.0f} sec"
            if eta_sec >= 60:
                eta_str = f"{eta_sec / 60:.1f} min"
            else:
                eta_str = f"{eta_sec:.0f} sec"
            print(f"  step {step:>8}/{n_steps} ({100*frac_done:5.1f}%) | elapsed {elapsed_str} | ETA {eta_str}")

    omega_samples_np = np.asarray(omega_samples.get() if device_is_gpu else omega_samples)[:sample_count]
    t_samples_np = np.asarray(t_samples.get() if device_is_gpu else t_samples)[:sample_count]

    if sample_count >= 2:
        omega_u = np.unwrap(omega_samples_np)
        t0 = t_samples_np - t_samples_np.mean()
        w0 = omega_u - omega_u.mean()
        den = float(np.sum(t0 * t0))
        slope = float(np.sum(t0 * w0) / den) if den > 0 else float("nan")
    else:
        slope = float("nan")

    mean_adv = slope * T  # rad/orbit

    arcsec_per_orbit = mean_adv * ARCSEC_PER_RAD
    arcsec_per_century_mercury = arcsec_per_orbit * ORBITS_PER_CENTURY_MERCURY

    results = {
        "mode": mode,
        "n_orbits_target": n_orbits,
        "steps_per_orbit": steps_per_orbit,
        "dt": dt,
        "sim_period": T,
        "perihelia_found": None,
        "mean_advance_rad_per_orbit": mean_adv,
        "mean_advance_arcsec_per_orbit": arcsec_per_orbit,
        "mean_advance_arcsec_per_century_mercury": arcsec_per_century_mercury,
    }
    return results, t_samples_np, omega_samples_np, trace


def write_csv(path: str, rows: List[Dict[str, Any]]):
    import csv
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out_dir", default="out_mercury", help="Output directory")
    ap.add_argument("--device", choices=["auto", "cpu", "gpu"], default="auto")
    ap.add_argument("--units", choices=["si", "dimensionless"], default="dimensionless")

    ap.add_argument("--mode", choices=["newton", "toy_r3", "beta_analytic", "beta_grid_poisson"],
                    default="beta_analytic")

    ap.add_argument("--n_orbits", type=int, default=600)
    ap.add_argument("--steps_per_orbit", type=int, default=10000)

    # SI inputs (used in si mode)
    ap.add_argument("--mu", type=float, default=MU_SUN_SI)
    ap.add_argument("--a", type=float, default=0.38709927 * AU_M)
    ap.add_argument("--e", type=float, default=0.205630)

    # toy parameters
    ap.add_argument("--cs", type=float, default=C_SI)
    ap.add_argument("--beta", type=float, default=1.5)

    # grid params for PDE mode
    ap.add_argument("--gridN", type=int, default=256)
    ap.add_argument("--boxL", type=float, default=128.0)
    ap.add_argument("--G", type=float, default=1.0)
    ap.add_argument("--sun_sigma_cells", type=float, default=1.5)

    # sweeps (comma-separated)
    ap.add_argument("--steps_sweep", default="", help="Comma-separated steps_per_orbit values")
    ap.add_argument("--cs_sweep", default="", help="Comma-separated cs values")
    ap.add_argument("--beta_sweep", default="", help="Comma-separated beta values")

    ap.add_argument("--save_trace", action="store_true")
    args = ap.parse_args()

    xp, is_gpu = get_xp(args.device)

    # Choose units
    if args.units == "si":
        mu, a, e = args.mu, args.a, args.e
        cs_default = args.cs
    else:
        # dimensionless demo: mu=a=1, same e; choose cs so epsilon matches Mercury-with-c
        # epsilon = mu/(cs^2 a (1-e^2))  -> match Mercury epsilon computed from SI constants
        e = args.e
        mu = 1.0
        a = 1.0
        eps_mercury = MU_SUN_SI / (C_SI**2 * (0.38709927 * AU_M) * (1.0 - args.e**2))
        cs_default = 1.0 / math.sqrt(eps_mercury * (1.0 - e**2))

    steps_list = [args.steps_per_orbit]
    cs_list = [cs_default]
    beta_list = [args.beta]

    if args.steps_sweep.strip():
        steps_list = [int(x) for x in args.steps_sweep.split(",") if x.strip()]
    if args.cs_sweep.strip():
        cs_list = [float(x) for x in args.cs_sweep.split(",") if x.strip()]
    if args.beta_sweep.strip():
        beta_list = [float(x) for x in args.beta_sweep.split(",") if x.strip()]

    all_runs = []
    artifacts_single = None

    for steps_per_orbit in steps_list:
        for cs in cs_list:
            for beta in beta_list:
                params = OrbitParams(mu=mu, a=a, e=e, cs=cs, beta=beta)
                grid_cfg = None
                if args.mode == "beta_grid_poisson":
                    grid_cfg = {
                        "gridN": args.gridN,
                        "boxL": args.boxL,
                        "G": args.G,
                        "sun_sigma_cells": args.sun_sigma_cells,
                    }

                t0 = time.time()
                results, omega_t, omega_ang, trace = integrate_and_measure(
                    xp=xp,
                    mode=args.mode,
                    params=params,
                    n_orbits=args.n_orbits,
                    steps_per_orbit=steps_per_orbit,
                    device_is_gpu=is_gpu,
                    grid_cfg=grid_cfg,
                    save_trace=args.save_trace,
                )

                results.update({
                    "units": args.units,
                    "mu": mu, "a": a, "e": e,
                    "cs": cs, "beta": beta,
                    "runtime_s": time.time() - t0,
                    "device": "gpu" if is_gpu else "cpu",
                })

                print(json.dumps(results, indent=2))
                all_runs.append(results)

                # keep artifacts only if single configuration
                if len(steps_list) == len(cs_list) == len(beta_list) == 1:
                    artifacts_single = (omega_t, omega_ang, trace)

    os.makedirs(args.out_dir, exist_ok=True)

    with open(os.path.join(args.out_dir, "summary.json"), "w") as f:
        json.dump({"runs": all_runs}, f, indent=2)

    write_csv(os.path.join(args.out_dir, "sweep.csv"), all_runs)

    if artifacts_single is not None:
        omega_t, omega_ang, trace = artifacts_single
        rows = []
        for k in range(len(omega_t)):
            rows.append({
                "k": k,
                "t": float(omega_t[k]),
                "angle_rad": float(omega_ang[k]),
                "angle_deg": float(omega_ang[k] * 180.0 / math.pi),
            })
        write_csv(os.path.join(args.out_dir, "omega_angles.csv"), rows)

        if args.save_trace and trace:
            write_csv(os.path.join(args.out_dir, "trace.csv"), trace)


if __name__ == "__main__":
    main()
