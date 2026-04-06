#!/usr/bin/env python3
"""
stepD_superfluid_acoustic_throat.py

Goal
----
Bridge stepC (static Poisson potential-flow from a throat sink) to *dynamics* derived from
(superfluid) hydrodynamics, without solving Poisson at each time step.

We evolve the linearized, irrotational, barotropic (acoustic) limit around a uniform background:

  v = ∇ψ
  ∂t ρ' + ρ0 ∇·v = Sρ
  ∂t v + (c_s^2/ρ0) ∇ρ' = 0

Eliminating ρ' gives a forced wave equation for ψ:

  ∂t^2 ψ - c_s^2 ∇^2 ψ = -(c_s^2/ρ0) Sρ

We optionally add a simple damping term γ ∂t ψ to absorb acoustic transients:

  ∂t^2 ψ + γ ∂t ψ - c_s^2 ∇^2 ψ = -(c_s^2/ρ0) Sρ

Static limit (∂t→0) gives Poisson:

  ∇^2 ψ = (Sρ/ρ0)  (with mean-subtraction for periodic solvability)

We model the "throat" as a localized sink with throughput Mdot and a uniform compensating
refill (a proxy for "dark-energy" leak-in) so that the net source integrates to zero in a periodic box.

Key outputs
-----------
- dpsi_clean_slope: fit of |ψ̄(r) - C_fit| ~ r^p. Target p ≈ -1
- g_slope: fit of |dψ̄/dr| ~ r^p. Target p ≈ -2
- lag_rms_over_psiP_rms: ||ψ - ψ_P|| / ||ψ_P|| (computed in k-space). Target → 0 as waves damp.

Backend
-------
Supports numpy (CPU) or cupy (GPU). Diagnostics are printed as JSON lines.

Notes on CuPy spherical averages
--------------------------------
We use cupyx.scatter_add in-place (it may return None in some CuPy versions), following
the robust pattern you posted.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from typing import Tuple, Optional

import numpy as np


def pick_backend(name: str):
    name = name.lower()
    if name == "numpy":
        return np
    if name == "cupy":
        import cupy as cp  # type: ignore
        return cp
    if name == "auto":
        try:
            import cupy as cp  # type: ignore
            _ = cp.zeros((1,))
            return cp
        except Exception:
            return np
    raise ValueError(f"unknown backend {name!r}")


def kgrid(xp, N: int, L: float, dtype):
    """Spectral wave numbers for periodic N^3 grid."""
    k1 = xp.fft.fftfreq(N, d=L / N).astype(dtype) * (2.0 * math.pi)
    kx = k1[:, None, None]
    ky = k1[None, :, None]
    kz = k1[None, None, :]
    k2 = kx * kx + ky * ky + kz * kz
    return kx, ky, kz, k2


def radial_bins_numpy(N: int, L: float, nbins: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Precompute radial bin index per cell for an N^3 periodic grid (centered)."""
    dx = L / N
    idx = np.arange(N)
    ii = idx - (N // 2)
    x = ii * dx
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rmax = (math.sqrt(3.0) * (L / 2.0))
    edges = np.linspace(0.0, rmax, nbins + 1)
    b = np.digitize(r.ravel(), edges) - 1
    b = np.clip(b, 0, nbins - 1)
    counts = np.bincount(b, minlength=nbins).astype(np.int64)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, b.astype(np.int32), counts


def spherical_average(xp, field, bin_idx_flat_np: np.ndarray, counts_np: np.ndarray, nbins: int, bin_idx_flat_gpu=None) -> np.ndarray:
    """Spherical average using precomputed bin indices; robust for cupy/numpy."""
    flat = field.ravel()
    if xp.__name__ == "cupy":
        import cupyx  # type: ignore

        if bin_idx_flat_gpu is None:
            bin_idx_flat_gpu = xp.asarray(bin_idx_flat_np, dtype=xp.int32)

        sums = xp.zeros(nbins, dtype=flat.dtype)
        cupyx.scatter_add(sums, bin_idx_flat_gpu, flat)

        denom = xp.asarray(counts_np)
        denom = xp.maximum(denom, 1).astype(flat.dtype)
        avg = sums / denom
        return avg.get()
    else:
        sums = np.bincount(bin_idx_flat_np, weights=np.asarray(flat), minlength=nbins)
        denom = np.maximum(counts_np, 1)
        return sums / denom


def fit_powerlaw(r: np.ndarray, y: np.ndarray, rmin: float, rmax: float, y_min_abs: float = 0.0) -> Tuple[float, float, int]:
    """Fit |y| ~ A r^p over [rmin,rmax] in log-log."""
    m = (r >= rmin) & (r <= rmax) & np.isfinite(y) & np.isfinite(r)
    rr = r[m]
    yabs = np.abs(y[m])
    if y_min_abs > 0:
        keep = yabs > y_min_abs
        rr = rr[keep]
        yabs = yabs[keep]
    if rr.size < 2:
        return float("nan"), float("nan"), 0
    logx = np.log(rr)
    logy = np.log(yabs)
    p, b = np.polyfit(logx, logy, 1)
    return float(p), float(b), int(rr.size)


def make_gaussian_kernel(xp, N: int, L: float, sigma: float, dtype):
    """Centered Gaussian kernel W normalized so that sum(W)*dx^3 = 1."""
    dx = L / N
    idx = xp.arange(N, dtype=dtype)
    ii = idx - (N // 2)
    x = ii * dx
    X, Y, Z = xp.meshgrid(x, x, x, indexing="ij")
    r2 = X * X + Y * Y + Z * Z
    W = xp.exp(-0.5 * r2 / (sigma * sigma))
    norm = xp.sum(W) * (dx**3)
    W = W / norm
    return W


def poisson_kernel_k(xp, qk_unit, k2):
    """Return psiP_k_unit solving -k2*psi = q (with k0 set to 0)."""
    psiP_k = xp.zeros_like(qk_unit)
    mask = k2 > 0
    psiP_k[mask] = -qk_unit[mask] / k2[mask]
    psiP_k[~mask] = 0
    return psiP_k


def l2_ratio_kspace(xp, a_k, b_k) -> float:
    """Return ||a||/||b|| in real-space RMS sense using Parseval (ratio cancels FFT scaling)."""
    # exclude k=0 automatically if b_k[0]=0; but keep anyway
    num = xp.sum(xp.abs(a_k) ** 2)
    den = xp.sum(xp.abs(b_k) ** 2)
    if den <= 0:
        return float("nan")
    return float(xp.sqrt(num / den).get() if xp.__name__ == "cupy" else math.sqrt(float(num / den)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", default="auto", choices=["auto", "numpy", "cupy"])
    ap.add_argument("--N", type=int, default=256)
    ap.add_argument("--L", type=float, default=200.0)
    ap.add_argument("--dtype", default="float32", choices=["float32", "float64"])

    # Throat/source
    ap.add_argument("--Mdot", type=float, default=20.0)
    ap.add_argument("--sigma", type=float, default=2.0)
    ap.add_argument("--rho0", type=float, default=1.0)
    ap.add_argument("--ramp_time", type=float, default=2.0)

    # EOS / sound speed
    ap.add_argument("--n", type=float, default=5.0)
    ap.add_argument("--P0", type=float, default=10.0)
    ap.add_argument("--cs", type=float, default=None)

    # Dynamics / damping
    ap.add_argument("--gamma", type=float, default=0.3, help="damping rate for ψ_t (absorber)")
    ap.add_argument("--steps", type=int, default=5000)
    ap.add_argument("--diag_every", type=int, default=200)
    ap.add_argument("--cfl", type=float, default=0.25)
    ap.add_argument("--dt", type=float, default=None)

    # Diagnostics
    ap.add_argument("--nbins", type=int, default=180)
    ap.add_argument("--fit_rmin", type=float, default=8.0)
    ap.add_argument("--fit_rmax", type=float, default=50.0)
    ap.add_argument("--sample_n", type=int, default=12)

    # Init
    ap.add_argument("--init", default="zero", choices=["zero", "poisson"], help="initial ψ; pi always 0")

    args = ap.parse_args()

    xp = pick_backend(args.backend)
    dtype = getattr(xp, args.dtype)
    cdtype = xp.complex64 if args.dtype == "float32" else xp.complex128

    N = int(args.N)
    L = float(args.L)
    dx = L / N
    V = L**3
    rho0 = float(args.rho0)

    cs = float(args.cs) if args.cs is not None else math.sqrt(float(args.n * args.P0 / rho0))

    # Conservative dt estimate: highest |k| ~ sqrt(3)*pi/dx -> omega_max = cs*|k|max
    omega_max = cs * (math.sqrt(3.0) * math.pi / dx)
    dt_max_est = 2.0 / omega_max  # loose RK stability heuristic
    dt = float(args.dt) if args.dt is not None else float(args.cfl * dt_max_est)

    # Precompute k-space operators
    kx, ky, kz, k2 = kgrid(xp, N, L, dtype)
    # Build throat kernel and mean-subtracted unit source q_unit = (Sρ/ρ0)/Mdot
    W = make_gaussian_kernel(xp, N, L, args.sigma, dtype=dtype)
    q_unit = (-W + (1.0 / V)) / rho0
    # Ensure exact mean-zero (periodic Poisson solvability)
    q_unit = q_unit - xp.mean(q_unit)

    qk_unit = xp.fft.fftn(q_unit.astype(dtype)).astype(cdtype)
    psiP_k_unit = poisson_kernel_k(xp, qk_unit, k2.astype(dtype)).astype(cdtype)

    # Initial conditions in k-space
    if args.init == "poisson":
        psi_k = (args.Mdot * psiP_k_unit).copy()
    else:
        psi_k = xp.zeros((N, N, N), dtype=cdtype)
    pi_k = xp.zeros_like(psi_k)

    # Enforce k=0 exactly zero
    psi_k[(0, 0, 0)] = 0
    pi_k[(0, 0, 0)] = 0

    # Radial bins (CPU) and optional GPU bin index for cupy
    r_centers, bin_idx_flat_np, counts_np = radial_bins_numpy(N, L, args.nbins)
    bin_idx_flat_gpu = None
    if xp.__name__ == "cupy":
        bin_idx_flat_gpu = xp.asarray(bin_idx_flat_np, dtype=xp.int32)

    init_msg = {
        "event": "init",
        "backend": xp.__name__,
        "N": N,
        "L": L,
        "dx": dx,
        "dtype": args.dtype,
        "params": {"rho0": rho0, "Mdot": args.Mdot, "sigma": args.sigma, "cs": cs, "gamma": args.gamma, "dt": dt, "dt_max_est": dt_max_est},
        "poisson": {"equation": "∇²ψ = (Sρ/ρ0) - <Sρ/ρ0>", "k0_mode": "set to 0"},
        "diag": {"nbins": args.nbins, "fit_rmin": args.fit_rmin, "fit_rmax": args.fit_rmax},
        "targets": {"dpsi_clean_slope": -1.0, "g_slope": -2.0, "lag_rms_over_psiP_rms": 0.0},
        "notes": "Evolve damped acoustic potential equation in k-space; compare to static Poisson solution.",
    }
    print(json.dumps(init_msg), flush=True)

    t = 0.0
    t0_wall = time.time()

    # Precompute constant operator for k-space RHS: -cs^2*k2
    k2_real = k2.astype(dtype)
    lap_op = (-cs * cs) * k2_real  # multiply psi_k by lap_op to get cs^2 ∇²ψ in k-space (since ∇² -> -k2)

    def mdot_of_time(tt: float) -> float:
        if args.ramp_time <= 0:
            return float(args.Mdot)
        x = max(0.0, min(1.0, tt / float(args.ramp_time)))
        # smoothstep
        s = x * x * (3.0 - 2.0 * x)
        return float(args.Mdot * s)

    def rhs_k(psi_k_local, pi_k_local, mdot: float):
        # psi_t = pi
        dpsi = pi_k_local
        # pi_t = -cs^2*k2*psi - gamma*pi - cs^2*mdot*qk_unit
        dpi = lap_op * psi_k_local - (args.gamma * pi_k_local) - (cs * cs) * (mdot * qk_unit)
        # Enforce k=0 mode stays zero
        dpsi = dpsi.copy()
        dpi = dpi.copy()
        dpsi[(0, 0, 0)] = 0
        dpi[(0, 0, 0)] = 0
        return dpsi, dpi

    for step in range(1, args.steps + 1):
        mdot0 = mdot_of_time(t)
        mdot1 = mdot_of_time(t + dt)
        mdot_half = 0.5 * (mdot0 + mdot1)

        # RK2 / Heun in k-space
        k1_psi, k1_pi = rhs_k(psi_k, pi_k, mdot_half)
        psi_k1 = psi_k + dt * k1_psi
        pi_k1 = pi_k + dt * k1_pi
        k2_psi, k2_pi = rhs_k(psi_k1, pi_k1, mdot_half)

        psi_k = psi_k + 0.5 * dt * (k1_psi + k2_psi)
        pi_k = pi_k + 0.5 * dt * (k1_pi + k2_pi)

        # pin k=0 again
        psi_k[(0, 0, 0)] = 0
        pi_k[(0, 0, 0)] = 0

        t += dt

        if (step % args.diag_every) == 0 or step == 1:
            # Convergence to Poisson reference in k-space (Parseval ratio)
            psiP_k_now = mdot_of_time(t) * psiP_k_unit
            lag_k = psi_k - psiP_k_now
            lag_ratio = l2_ratio_kspace(xp, lag_k, psiP_k_now)

            # Real-space ψ for radial fits
            psi = xp.fft.ifftn(psi_k).real.astype(dtype)

            psi_bar = spherical_average(
                xp, psi, bin_idx_flat_np, counts_np, args.nbins, bin_idx_flat_gpu=bin_idx_flat_gpu
            )

            outer_n = max(5, int(0.10 * args.nbins))
            psi_ref = float(np.mean(psi_bar[-outer_n:]))
            dpsi = psi_bar - psi_ref

            dpsidr = np.gradient(psi_bar, r_centers)
            g_r = np.abs(dpsidr)

            mask_fit = (r_centers >= args.fit_rmin) & (r_centers <= args.fit_rmax)
            dpsi_fit_abs = np.abs(dpsi[mask_fit])
            g_fit_abs = np.abs(g_r[mask_fit])
            dpsi_peak = float(np.nanmax(dpsi_fit_abs)) if dpsi_fit_abs.size else 0.0
            g_peak = float(np.nanmax(g_fit_abs)) if g_fit_abs.size else 0.0
            ymin_dpsi = 0.0 if dpsi_peak < 1e-30 else max(1e-30, dpsi_peak * 1e-3)
            ymin_g = 0.0 if g_peak < 1e-30 else max(1e-30, g_peak * 1e-3)

            p_dpsi, b_dpsi, n_dpsi = fit_powerlaw(r_centers, dpsi, args.fit_rmin, args.fit_rmax, y_min_abs=ymin_dpsi)
            p_g, b_g, n_g = fit_powerlaw(r_centers, g_r, args.fit_rmin, args.fit_rmax, y_min_abs=ymin_g)

            # Robust monopole fit: ψ̄ ≈ A/r + C
            A_fit = float("nan")
            C_fit = float("nan")
            rel_err_A = float("nan")
            A_expected = float(mdot_of_time(t) / (4.0 * math.pi * rho0))
            try:
                r_fit = r_centers[mask_fit]
                y_fit = psi_bar[mask_fit]
                good = np.isfinite(r_fit) & np.isfinite(y_fit) & (r_fit > 0.0)
                if np.count_nonzero(good) >= 2:
                    X = np.vstack([1.0 / r_fit[good], np.ones_like(r_fit[good])]).T
                    coeff, *_ = np.linalg.lstsq(X, y_fit[good], rcond=None)
                    A_fit = float(coeff[0])
                    C_fit = float(coeff[1])
                    if abs(A_expected) > 0:
                        rel_err_A = float((A_fit - A_expected) / A_expected)
            except Exception:
                pass

            dpsi_clean = psi_bar - C_fit
            dpsi_clean_fit_abs = np.abs(dpsi_clean[mask_fit])
            dpsi_clean_peak = float(np.nanmax(dpsi_clean_fit_abs)) if dpsi_clean_fit_abs.size else 0.0
            ymin_dpsi_clean = 0.0 if dpsi_clean_peak < 1e-30 else max(1e-30, dpsi_clean_peak * 1e-3)
            p_dpsi_clean, b_dpsi_clean, n_dpsi_clean = fit_powerlaw(
                r_centers, dpsi_clean, args.fit_rmin, args.fit_rmax, y_min_abs=ymin_dpsi_clean
            )

            # Optional: estimate Mdot from radial derivative (flux check)
            mdot_est = float("nan")
            try:
                vr = dpsidr  # for spherically symmetric ψ, v_r ≈ dψ̄/dr
                mdot_prof = -4.0 * math.pi * rho0 * (r_centers**2) * vr
                mdot_est = float(np.nanmean(mdot_prof[mask_fit]))
            except Exception:
                pass

            # Sample points
            idx = np.where(mask_fit)[0]
            stride = max(1, len(idx) // max(1, int(args.sample_n)))
            idx_s = idx[::stride]
            sample = [
                {
                    "r": float(r_centers[i]),
                    "dpsi_clean": float(dpsi_clean[i]),
                    "g": float(g_r[i]),
                    "psi_bar": float(psi_bar[i]),
                }
                for i in idx_s
            ]

            out = {
                "event": "diag",
                "step": step,
                "t": t,
                "dt": dt,
                "Mdot": mdot_of_time(t),
                "lag_rms_over_psiP_rms": lag_ratio,
                "refs": {"psi_ref": psi_ref},
                "fits": {
                    "dpsi_slope": p_dpsi,
                    "dpsi_npts": n_dpsi,
                    "dpsi_clean_slope": p_dpsi_clean,
                    "dpsi_clean_npts": n_dpsi_clean,
                    "g_slope": p_g,
                    "g_npts": n_g,
                },
                "monopole_fit": {
                    "A_fit": A_fit,
                    "C_fit": C_fit,
                    "A_expected": A_expected,
                    "rel_err_A": rel_err_A,
                    "Mdot_from_flux": mdot_est,
                },
                "fit_cutoffs": {
                    "ymin_dpsi": ymin_dpsi,
                    "ymin_dpsi_clean": ymin_dpsi_clean,
                    "ymin_g": ymin_g,
                    "dpsi_peak_fit": dpsi_peak,
                    "dpsi_clean_peak_fit": dpsi_clean_peak,
                    "g_peak_fit": g_peak,
                },
                "sample_points": sample,
                "wall_s": (time.time() - t0_wall),
            }
            print(json.dumps(out), flush=True)

    done = {"event": "done", "t": t, "steps": args.steps, "wall_s": (time.time() - t0_wall)}
    print(json.dumps(done), flush=True)


if __name__ == "__main__":
    main()
