#!/usr/bin/env python3
"""
mercury_suite.py

Publish-grade perihelion-precession benchmark harness for the toy 1PN / superfluid-analog model.

Outputs (out_dir):
- summary.json      (all run metadata + results)
- sweep.csv         (one row per configuration)
- perihelion_angles.csv  (for the single-run case; omitted during sweeps)
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

    sigma = (beta * mu) / (cs**2 * r_safe)
    grad_sigma = (-(beta * mu) / (cs**2 * r_safe**2))[:, None] * r_hat

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

    # perihelion detection buffers
    max_peri = n_orbits + 10
    peri_angles = xp.full((max_peri,), xp.nan, dtype=xp.float64)
    peri_times = xp.full((max_peri,), xp.nan, dtype=xp.float64)
    peri_count = 0

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

    # init radii buffers
    r_curr = float(np.linalg.norm(np.asarray(pos.get() if device_is_gpu else pos)[0]))
    r_prev = r_curr
    r_prevprev = r_curr

    # integrate: kick-drift-kick using v_half for v-dependent forces
    t = 0.0
    start_time = time.time()
    log_interval = max(1, n_steps // 100)  # log ~20 times during run
    for step in range(n_steps):
        a0 = accel(pos, vel)
        vel_half = vel + 0.5 * dt * a0
        pos_new = pos + dt * vel_half
        a1 = accel(pos_new, vel_half)
        vel_new = vel_half + 0.5 * dt * a1

        # perihelion detect at previous point (pos)
        r_new = xp.sqrt(xp.sum(pos_new * pos_new, axis=1))[0]
        r_new_f = float(r_new.get() if device_is_gpu else r_new)

        if (r_prev < r_prevprev) and (r_prev < r_new_f) and (step > 2):
            ang = xp.arctan2(pos[0, 1], pos[0, 0])
            peri_angles[peri_count] = ang
            peri_times[peri_count] = t
            peri_count += 1
            if peri_count >= max_peri:
                break

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

        r_prevprev, r_prev = r_prev, r_new_f
        pos, vel = pos_new, vel_new
        t += dt

    # fetch perihelia arrays
    peri_angles_np = np.asarray(peri_angles.get() if device_is_gpu else peri_angles)[:peri_count]
    peri_times_np = np.asarray(peri_times.get() if device_is_gpu else peri_times)[:peri_count]

    if peri_count >= 2:
        unwrapped = np.unwrap(peri_angles_np)
        diffs = np.diff(unwrapped)
        mean_adv = float(np.mean(diffs))
        std_adv = float(np.std(diffs))
    else:
        mean_adv = float("nan")
        std_adv = float("nan")

    arcsec_per_orbit = mean_adv * ARCSEC_PER_RAD
    arcsec_per_century_mercury = arcsec_per_orbit * ORBITS_PER_CENTURY_MERCURY

    results = {
        "mode": mode,
        "n_orbits_target": n_orbits,
        "steps_per_orbit": steps_per_orbit,
        "dt": dt,
        "sim_period": T,
        "perihelia_found": int(peri_count),
        "mean_advance_rad_per_orbit": mean_adv,
        "std_advance_rad_per_orbit": std_adv,
        "mean_advance_arcsec_per_orbit": arcsec_per_orbit,
        "mean_advance_arcsec_per_century_mercury": arcsec_per_century_mercury,
    }
    return results, peri_times_np, peri_angles_np, trace


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
                results, peri_t, peri_ang, trace = integrate_and_measure(
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
                    artifacts_single = (peri_t, peri_ang, trace)

    os.makedirs(args.out_dir, exist_ok=True)

    with open(os.path.join(args.out_dir, "summary.json"), "w") as f:
        json.dump({"runs": all_runs}, f, indent=2)

    write_csv(os.path.join(args.out_dir, "sweep.csv"), all_runs)

    if artifacts_single is not None:
        peri_t, peri_ang, trace = artifacts_single
        rows = []
        for k in range(len(peri_t)):
            rows.append({
                "k": k,
                "t": float(peri_t[k]),
                "angle_rad": float(peri_ang[k]),
                "angle_deg": float(peri_ang[k] * 180.0 / math.pi),
            })
        write_csv(os.path.join(args.out_dir, "perihelion_angles.csv"), rows)

        if args.save_trace and trace:
            write_csv(os.path.join(args.out_dir, "trace.csv"), trace)


if __name__ == "__main__":
    main()

