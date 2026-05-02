#!/usr/bin/env python3
"""
single_throat_monopole.py

Minimal "Newtonian sector" test for the toy model:
- 3D conservative compressible Euler (barotropic) on a GPU using CuPy
- EOS: stiff polytrope P = P0 * (rho/rho0)^n, default n=5
- One fixed "throat" modeled as a normalized kernel W(x) removing superfluid density
  with a prescribed throughput Mdot(t) (ramped up smoothly).
- Bulk refill term (dark-energy leakage) drives rho back toward rho0 when rho < rho0.
- Reservoir boundary: ghost cells fixed to rho=rho0, v=0.

Diagnostics (prints copy/paste-friendly JSON blocks):
- Spherical averages of rho(r), P(r), ΔP(r)
- g_eff(r) = (1/rho(r)) dP/dr
- Power-law slope fits: ΔP ~ r^pP (target -1), g_eff ~ r^pg (target -2)

Notes:
- This is a robust first-order Rusanov (LLF) flux scheme with RK2 time stepping.
- Intended for subsonic/moderate flows (your "weak-field / acoustic-ish" regime).
- For N=512 this will be heavy but should fit comfortably on 80–100 GB GPUs.

Run:
  python single_throat_monopole.py --N 256 --steps 20000 --diag_every 200

"""

import argparse
import json
import math
import time

import cupy as cp
import numpy as np


# -----------------------------
# Utilities
# -----------------------------

def enforce_reservoir_bc(rho, mx, my, mz, rho0, g):
    """Set ghost cells to a reservoir: rho=rho0, momentum=0."""
    # x faces
    rho[:g, :, :] = rho0
    rho[-g:, :, :] = rho0
    mx[:g, :, :] = 0; mx[-g:, :, :] = 0
    my[:g, :, :] = 0; my[-g:, :, :] = 0
    mz[:g, :, :] = 0; mz[-g:, :, :] = 0

    # y faces
    rho[:, :g, :] = rho0
    rho[:, -g:, :] = rho0
    mx[:, :g, :] = 0; mx[:, -g:, :] = 0
    my[:, :g, :] = 0; my[:, -g:, :] = 0
    mz[:, :g, :] = 0; mz[:, -g:, :] = 0

    # z faces
    rho[:, :, :g] = rho0
    rho[:, :, -g:] = rho0
    mx[:, :, :g] = 0; mx[:, :, -g:] = 0
    my[:, :, :g] = 0; my[:, :, -g:] = 0
    mz[:, :, :g] = 0; mz[:, :, -g:] = 0


def eos_polytrope(rho, rho0, P0, n, rho_floor=1e-12):
    """Return (P, cs) for P = P0 * (rho/rho0)^n, cs^2 = dP/drho = n P / rho."""
    rho_safe = cp.maximum(rho, rho_floor)
    P = P0 * cp.power(rho_safe / rho0, n)
    cs2 = n * P / rho_safe
    cs = cp.sqrt(cp.maximum(cs2, 0.0))
    return P, cs


def compute_dt_cfl(vx, vy, vz, cs, dx, cfl, g):
    """Compute global stable dt = CFL * dx / max(|v|+cs) across interior."""
    # Interior slices
    sl = slice(g, -g)
    vabs_x = cp.abs(vx[sl, sl, sl]) + cs[sl, sl, sl]
    vabs_y = cp.abs(vy[sl, sl, sl]) + cs[sl, sl, sl]
    vabs_z = cp.abs(vz[sl, sl, sl]) + cs[sl, sl, sl]
    max_wave = cp.maximum(cp.maximum(cp.max(vabs_x), cp.max(vabs_y)), cp.max(vabs_z))
    max_wave_cpu = float(max_wave.get())
    if max_wave_cpu <= 1e-12:
        return 1e-6, 0.0
    return cfl * dx / max_wave_cpu, max_wave_cpu


def make_throat_kernel(N, L, dx, a_phys, g):
    """
    Create a normalized Gaussian kernel W over the *interior* (N^3),
    normalized such that sum(W)*dx^3 = 1.

    a_phys: Gaussian width parameter (same physical units as L).
    """
    # Interior coordinates (cell centers)
    x = (cp.arange(N, dtype=cp.float32) + 0.5) * dx
    X, Y, Z = cp.meshgrid(x, x, x, indexing="ij")
    c = 0.5 * L
    r2 = (X - c) ** 2 + (Y - c) ** 2 + (Z - c) ** 2
    W = cp.exp(-r2 / (a_phys ** 2)).astype(cp.float32)

    # Normalize: integral W dV = 1
    norm = cp.sum(W) * (dx ** 3)
    W /= norm
    return W


def precompute_radial_bins(N, L, dx, nbins, g):
    """
    Precompute r (bin centers), bin_idx for each interior cell, and counts per bin.
    Uses integer binning based on dr = r_max/nbins.

    Returns:
      r_centers (numpy),
      bin_idx_flat (cupy int32),
      counts (cupy float32)
    """
    x = (cp.arange(N, dtype=cp.float32) + 0.5) * dx
    X, Y, Z = cp.meshgrid(x, x, x, indexing="ij")
    c = 0.5 * L
    r = cp.sqrt((X - c) ** 2 + (Y - c) ** 2 + (Z - c) ** 2)

    r_max = 0.5 * L
    dr = r_max / nbins
    # bin index in [0, nbins-1], clamp
    b = cp.floor(r / dr).astype(cp.int32)

    # Put anything outside the inscribed sphere (r >= r_max) into an overflow bin = nbins
    b = cp.where(r < r_max, b, nbins)

    b_flat = b.ravel()
    ones = cp.ones_like(b_flat, dtype=cp.float64)

    # Allocate nbins+1 so overflow is counted separately, then discard overflow
    counts_all = cp.bincount(b_flat, weights=ones, minlength=nbins + 1).astype(cp.float64)
    counts = counts_all[:nbins]

    r_centers = (np.arange(nbins, dtype=np.float64) + 0.5) * dr
    return r_centers, b_flat, counts, dr


def spherical_average(field_interior, bin_idx_flat, counts, nbins):
    """Return spherical average over bins using bincount on GPU."""
    f_flat = field_interior.ravel().astype(cp.float64)
    sums_all = cp.bincount(bin_idx_flat, weights=f_flat, minlength=nbins + 1).astype(cp.float64)
    sums = sums_all[:nbins]
    avg = sums / cp.maximum(counts, 1e-20)
    return avg


def fit_powerlaw(x, y, x_min, x_max, y_min_abs=0.0):
    """Fit log-log slope y ~ x^p over window, returns (p, intercept, npts)."""
    mask = (x >= x_min) & (x <= x_max) & (np.isfinite(y)) & (np.abs(y) > y_min_abs)
    xx = x[mask]
    yy = y[mask]
    if len(xx) < 5:
        return float("nan"), float("nan"), int(len(xx))
    lx = np.log(xx)
    ly = np.log(np.abs(yy))
    p, b = np.polyfit(lx, ly, 1)
    return float(p), float(b), int(len(xx))


def _stop_or_none(stop: int):
    # In Python slicing, stop=0 means "stop before index 0" (empty if start>0).
    # We want "to the end" instead.
    return None if stop == 0 else stop


def save_projected_snapshot(
    path,
    *,
    rho,
    mx,
    my,
    mz,
    W,
    dt_rho,
    N,
    L,
    dx,
    g,
    Mdot,
    lam,
    rho0,
    P0,
    n,
    a,
    t_sim,
    step,
):
    """Write a projected 3D monitor snapshot for cfd_snapshot_adapters.py."""
    sl = slice(g, -g)
    x = (np.arange(N, dtype=np.float64) + 0.5) * dx
    W_cpu = cp.asnumpy(W).astype(np.float64, copy=False)
    V_domain = float(L) ** 3
    S_rho = (-float(Mdot)) * W_cpu + (float(lam) * float(Mdot) / V_domain)

    payload = {
        "rho": cp.asnumpy(rho[sl, sl, sl]).astype(np.float32, copy=False),
        "mx": cp.asnumpy(mx[sl, sl, sl]).astype(np.float32, copy=False),
        "my": cp.asnumpy(my[sl, sl, sl]).astype(np.float32, copy=False),
        "mz": cp.asnumpy(mz[sl, sl, sl]).astype(np.float32, copy=False),
        "S_rho": S_rho.astype(np.float32, copy=False),
        "W": W_cpu.astype(np.float32, copy=False),
        "x": x,
        "y": x.copy(),
        "z": x.copy(),
        "Mdot": np.array(float(Mdot), dtype=np.float64),
        "lambda_bulk": np.array(float(lam), dtype=np.float64),
        "V_domain": np.array(V_domain, dtype=np.float64),
        "rho0": np.array(float(rho0), dtype=np.float64),
        "P0": np.array(float(P0), dtype=np.float64),
        "n": np.array(float(n), dtype=np.float64),
        "a": np.array(float(a), dtype=np.float64),
        "dx": np.array(float(dx), dtype=np.float64),
        "t": np.array(float(t_sim), dtype=np.float64),
        "step": np.array(int(step), dtype=np.int64),
    }
    if dt_rho is not None:
        payload["dt_rho"] = cp.asnumpy(dt_rho).astype(np.float32, copy=False)
    np.savez_compressed(path, **payload)


# -----------------------------
# Rusanov flux divergence
# -----------------------------

def rusanov_divergence(
    rho, mx, my, mz, P, cs,
    dx, g
):
    """
    Compute -∇·F(U) on the interior using a memory-friendly Rusanov (LLF) scheme.
    Returns drho, dmx, dmy, dmz for interior cells shape (N,N,N).
    """
    slc = slice(g, -g)
    rho_c = rho[slc, slc, slc]
    mx_c = mx[slc, slc, slc]
    my_c = my[slc, slc, slc]
    mz_c = mz[slc, slc, slc]
    P_c = P[slc, slc, slc]
    cs_c = cs[slc, slc, slc]

    # Velocities (interior)
    rho_safe = cp.maximum(rho_c, 1e-12)
    vx_c = mx_c / rho_safe
    vy_c = my_c / rho_safe
    vz_c = mz_c / rho_safe

    # Allocate RHS (interior)
    drho = cp.zeros_like(rho_c)
    dmx = cp.zeros_like(mx_c)
    dmy = cp.zeros_like(my_c)
    dmz = cp.zeros_like(mz_c)

    stop_p = _stop_or_none(-g + 1)

    # ----- X direction -----
    rho_p = rho[g + 1: stop_p, slc, slc]
    mx_p  = mx [g + 1: stop_p, slc, slc]
    my_p  = my [g + 1: stop_p, slc, slc]
    mz_p  = mz [g + 1: stop_p, slc, slc]
    P_p   = P  [g + 1: stop_p, slc, slc]
    cs_p  = cs [g + 1: stop_p, slc, slc]

    rho_m = rho[g - 1: -g - 1, slc, slc]
    mx_m  = mx [g - 1: -g - 1, slc, slc]
    my_m  = my [g - 1: -g - 1, slc, slc]
    mz_m  = mz [g - 1: -g - 1, slc, slc]
    P_m   = P  [g - 1: -g - 1, slc, slc]
    cs_m  = cs [g - 1: -g - 1, slc, slc]

    vx_p = mx_p / cp.maximum(rho_p, 1e-12)
    vx_m = mx_m / cp.maximum(rho_m, 1e-12)

    s_p = cp.maximum(cp.abs(vx_c) + cs_c, cp.abs(vx_p) + cs_p)
    s_m = cp.maximum(cp.abs(vx_m) + cs_m, cp.abs(vx_c) + cs_c)

    # Fluxes at cell centers for x-direction
    # F_rho = mx
    F_rho_c = mx_c
    F_rho_p = mx_p
    F_rho_m = mx_m

    # F_mx = mx*vx + P
    F_mx_c = mx_c * vx_c + P_c
    F_mx_p = mx_p * vx_p + P_p
    F_mx_m = mx_m * vx_m + P_m

    # F_my = my*vx
    F_my_c = my_c * vx_c
    F_my_p = my_p * vx_p
    F_my_m = my_m * vx_m

    # F_mz = mz*vx
    F_mz_c = mz_c * vx_c
    F_mz_p = mz_p * vx_p
    F_mz_m = mz_m * vx_m

    # Rusanov interface fluxes (computed on interior grid)
    # plus face between c and p
    fluxp_rho = 0.5 * (F_rho_c + F_rho_p) - 0.5 * s_p * (rho_p - rho_c)
    fluxp_mx  = 0.5 * (F_mx_c  + F_mx_p ) - 0.5 * s_p * (mx_p  - mx_c)
    fluxp_my  = 0.5 * (F_my_c  + F_my_p ) - 0.5 * s_p * (my_p  - my_c)
    fluxp_mz  = 0.5 * (F_mz_c  + F_mz_p ) - 0.5 * s_p * (mz_p  - mz_c)

    # minus face between m and c
    fluxm_rho = 0.5 * (F_rho_m + F_rho_c) - 0.5 * s_m * (rho_c - rho_m)
    fluxm_mx  = 0.5 * (F_mx_m  + F_mx_c ) - 0.5 * s_m * (mx_c  - mx_m)
    fluxm_my  = 0.5 * (F_my_m  + F_my_c ) - 0.5 * s_m * (my_c  - my_m)
    fluxm_mz  = 0.5 * (F_mz_m  + F_mz_c ) - 0.5 * s_m * (mz_c  - mz_m)

    drho += -(fluxp_rho - fluxm_rho) / dx
    dmx  += -(fluxp_mx  - fluxm_mx ) / dx
    dmy  += -(fluxp_my  - fluxm_my ) / dx
    dmz  += -(fluxp_mz  - fluxm_mz ) / dx

    # ----- Y direction -----
    rho_p = rho[slc, g + 1: stop_p, slc]
    mx_p  = mx [slc, g + 1: stop_p, slc]
    my_p  = my [slc, g + 1: stop_p, slc]
    mz_p  = mz [slc, g + 1: stop_p, slc]
    P_p   = P  [slc, g + 1: stop_p, slc]
    cs_p  = cs [slc, g + 1: stop_p, slc]

    rho_m = rho[slc, g - 1: -g - 1, slc]
    mx_m  = mx [slc, g - 1: -g - 1, slc]
    my_m  = my [slc, g - 1: -g - 1, slc]
    mz_m  = mz [slc, g - 1: -g - 1, slc]
    P_m   = P  [slc, g - 1: -g - 1, slc]
    cs_m  = cs [slc, g - 1: -g - 1, slc]

    vy_p = my_p / cp.maximum(rho_p, 1e-12)
    vy_m = my_m / cp.maximum(rho_m, 1e-12)

    s_p = cp.maximum(cp.abs(vy_c) + cs_c, cp.abs(vy_p) + cs_p)
    s_m = cp.maximum(cp.abs(vy_m) + cs_m, cp.abs(vy_c) + cs_c)

    # Fluxes at cell centers for y-direction
    # G_rho = my
    G_rho_c = my_c
    G_rho_p = my_p
    G_rho_m = my_m

    # G_mx = mx*vy
    G_mx_c = mx_c * vy_c
    G_mx_p = mx_p * vy_p
    G_mx_m = mx_m * vy_m

    # G_my = my*vy + P
    G_my_c = my_c * vy_c + P_c
    G_my_p = my_p * vy_p + P_p
    G_my_m = my_m * vy_m + P_m

    # G_mz = mz*vy
    G_mz_c = mz_c * vy_c
    G_mz_p = mz_p * vy_p
    G_mz_m = mz_m * vy_m

    fluxp_rho = 0.5 * (G_rho_c + G_rho_p) - 0.5 * s_p * (rho_p - rho_c)
    fluxp_mx  = 0.5 * (G_mx_c  + G_mx_p ) - 0.5 * s_p * (mx_p  - mx_c)
    fluxp_my  = 0.5 * (G_my_c  + G_my_p ) - 0.5 * s_p * (my_p  - my_c)
    fluxp_mz  = 0.5 * (G_mz_c  + G_mz_p ) - 0.5 * s_p * (mz_p  - mz_c)

    fluxm_rho = 0.5 * (G_rho_m + G_rho_c) - 0.5 * s_m * (rho_c - rho_m)
    fluxm_mx  = 0.5 * (G_mx_m  + G_mx_c ) - 0.5 * s_m * (mx_c  - mx_m)
    fluxm_my  = 0.5 * (G_my_m  + G_my_c ) - 0.5 * s_m * (my_c  - my_m)
    fluxm_mz  = 0.5 * (G_mz_m  + G_mz_c ) - 0.5 * s_m * (mz_c  - mz_m)

    drho += -(fluxp_rho - fluxm_rho) / dx
    dmx  += -(fluxp_mx  - fluxm_mx ) / dx
    dmy  += -(fluxp_my  - fluxm_my ) / dx
    dmz  += -(fluxp_mz  - fluxm_mz ) / dx

    # ----- Z direction -----
    rho_p = rho[slc, slc, g + 1: stop_p]
    mx_p  = mx [slc, slc, g + 1: stop_p]
    my_p  = my [slc, slc, g + 1: stop_p]
    mz_p  = mz [slc, slc, g + 1: stop_p]
    P_p   = P  [slc, slc, g + 1: stop_p]
    cs_p  = cs [slc, slc, g + 1: stop_p]

    rho_m = rho[slc, slc, g - 1: -g - 1]
    mx_m  = mx [slc, slc, g - 1: -g - 1]
    my_m  = my [slc, slc, g - 1: -g - 1]
    mz_m  = mz [slc, slc, g - 1: -g - 1]
    P_m   = P  [slc, slc, g - 1: -g - 1]
    cs_m  = cs [slc, slc, g - 1: -g - 1]

    vz_p = mz_p / cp.maximum(rho_p, 1e-12)
    vz_m = mz_m / cp.maximum(rho_m, 1e-12)

    s_p = cp.maximum(cp.abs(vz_c) + cs_c, cp.abs(vz_p) + cs_p)
    s_m = cp.maximum(cp.abs(vz_m) + cs_m, cp.abs(vz_c) + cs_c)

    # Fluxes at cell centers for z-direction
    # H_rho = mz
    H_rho_c = mz_c
    H_rho_p = mz_p
    H_rho_m = mz_m

    # H_mx = mx*vz
    H_mx_c = mx_c * vz_c
    H_mx_p = mx_p * vz_p
    H_mx_m = mx_m * vz_m

    # H_my = my*vz
    H_my_c = my_c * vz_c
    H_my_p = my_p * vz_p
    H_my_m = my_m * vz_m

    # H_mz = mz*vz + P
    H_mz_c = mz_c * vz_c + P_c
    H_mz_p = mz_p * vz_p + P_p
    H_mz_m = mz_m * vz_m + P_m

    fluxp_rho = 0.5 * (H_rho_c + H_rho_p) - 0.5 * s_p * (rho_p - rho_c)
    fluxp_mx  = 0.5 * (H_mx_c  + H_mx_p ) - 0.5 * s_p * (mx_p  - mx_c)
    fluxp_my  = 0.5 * (H_my_c  + H_my_p ) - 0.5 * s_p * (my_p  - my_c)
    fluxp_mz  = 0.5 * (H_mz_c  + H_mz_p ) - 0.5 * s_p * (mz_p  - mz_c)

    fluxm_rho = 0.5 * (H_rho_m + H_rho_c) - 0.5 * s_m * (rho_c - rho_m)
    fluxm_mx  = 0.5 * (H_mx_m  + H_mx_c ) - 0.5 * s_m * (mx_c  - mx_m)
    fluxm_my  = 0.5 * (H_my_m  + H_my_c ) - 0.5 * s_m * (my_c  - my_m)
    fluxm_mz  = 0.5 * (H_mz_m  + H_mz_c ) - 0.5 * s_m * (mz_c  - mz_m)

    drho += -(fluxp_rho - fluxm_rho) / dx
    dmx  += -(fluxp_mx  - fluxm_mx ) / dx
    dmy  += -(fluxp_my  - fluxm_my ) / dx
    dmz  += -(fluxp_mz  - fluxm_mz ) / dx

    return drho, dmx, dmy, dmz


# -----------------------------
# Main
# -----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=256, help="Interior grid resolution per axis.")
    ap.add_argument("--L", type=float, default=100.0, help="Domain size.")
    ap.add_argument("--n", type=float, default=5.0, help="Polytropic exponent n.")
    ap.add_argument("--rho0", type=float, default=1.0, help="Background density.")
    ap.add_argument("--P0", type=float, default=10.0, help="Background pressure scale.")
    ap.add_argument("--steps", type=int, default=20000, help="Number of steps.")
    ap.add_argument("--cfl", type=float, default=0.25, help="CFL number.")
    ap.add_argument("--diag_every", type=int, default=200, help="Diagnostics interval (steps).")
    ap.add_argument("--nbins", type=int, default=160, help="Radial bins for spherical averages.")
    ap.add_argument("--a", type=float, default=1.0, help="Throat kernel width a (physical units).")
    ap.add_argument("--Mdot", type=float, default=2.0, help="Target suction throughput (density*vol/time).")
    ap.add_argument("--ramp_time", type=float, default=2.0, help="Ramp time for Mdot (physical time units).")
    ap.add_argument("--lambda_refill", type=float, default=0.5, help="Refill strength lambda.")
    ap.add_argument("--refill_momentum", choices=["comoving", "at_rest"], default="comoving",
                    help="Momentum added with refill: comoving (adds with local v) or at_rest (adds zero momentum).")
    ap.add_argument("--rho_floor", type=float, default=1e-4, help="Density floor for numerical safety.")
    ap.add_argument("--fit_rmin_factor", type=float, default=6.0, help="Fit window r_min = factor*a.")
    ap.add_argument("--fit_rmax_frac", type=float, default=0.35, help="Fit window r_max = frac*L.")
    ap.add_argument("--nu_drag", type=float, default=0.0, help="Linear momentum drag rate (1/time).")
    ap.add_argument("--snapshot_out", default=None, help="Write final projected 3D state as an NPZ for the runtime monitor.")
    args = ap.parse_args()

    N = args.N
    L = float(args.L)
    dx = L / N
    n = float(args.n)
    rho0 = float(args.rho0)
    P0 = float(args.P0)
    cfl = float(args.cfl)
    steps = int(args.steps)
    diag_every = int(args.diag_every)
    nbins = int(args.nbins)
    a = float(args.a)
    Mdot_target = float(args.Mdot)
    ramp_time = float(args.ramp_time)
    lam = float(args.lambda_refill)
    refill_momentum = args.refill_momentum
    nu_drag = float(args.nu_drag)
    rho_floor = float(args.rho_floor)

    # Ghost thickness (1 cell for 1st order)
    g = 1
    Nx = N + 2 * g

    # Allocate full arrays incl. ghost cells
    rho = cp.full((Nx, Nx, Nx), rho0, dtype=cp.float32)
    mx  = cp.zeros((Nx, Nx, Nx), dtype=cp.float32)
    my  = cp.zeros((Nx, Nx, Nx), dtype=cp.float32)
    mz  = cp.zeros((Nx, Nx, Nx), dtype=cp.float32)

    enforce_reservoir_bc(rho, mx, my, mz, rho0=rho0, g=g)

    # Precompute throat kernel W on interior
    W = make_throat_kernel(N=N, L=L, dx=dx, a_phys=a, g=g)

    # Precompute radial bins (interior)
    r_centers, bin_idx_flat, counts, dr = precompute_radial_bins(N=N, L=L, dx=dx, nbins=nbins, g=g)

    # Precompute volume element and some constants
    dV = dx ** 3

    # Tracking
    t_sim = 0.0
    wall0 = time.time()

    print(json.dumps({
        "event": "init",
        "N": N, "L": L, "dx": dx, "ghost": g, "nbins": nbins,
        "EOS": {"type": "polytrope", "n": n, "rho0": rho0, "P0": P0},
        "throat": {"a": a, "Mdot_target": Mdot_target, "ramp_time": ramp_time},
        "refill": {"lambda": lam, "momentum_mode": refill_momentum},
        "notes": "Single fixed throat monopole test. Targets: slope(ΔP)≈-1, slope(g_eff)≈-2."
    }, indent=None))

    # Main loop
    for step in range(steps):
        # EOS and primitive variables
        rho = cp.maximum(rho, rho_floor)
        enforce_reservoir_bc(rho, mx, my, mz, rho0=rho0, g=g)

        P, cs = eos_polytrope(rho, rho0=rho0, P0=P0, n=n, rho_floor=rho_floor)

        # Primitive velocities for dt and sources
        rho_safe_full = cp.maximum(rho, rho_floor)
        vx = mx / rho_safe_full
        vy = my / rho_safe_full
        vz = mz / rho_safe_full

        # CFL dt
        dt, max_wave = compute_dt_cfl(vx, vy, vz, cs, dx=dx, cfl=cfl, g=g)
        t_next = t_sim + dt

        # Mdot ramp
        if ramp_time <= 0:
            Mdot = Mdot_target
        else:
            # smooth ramp 0->1 using sin^2
            s = min(max(t_sim / ramp_time, 0.0), 1.0)
            ramp = math.sin(0.5 * math.pi * s) ** 2
            Mdot = Mdot_target * ramp

        # RHS function (interior), includes sources
        def rhs(rho_in, mx_in, my_in, mz_in, t_local):
            enforce_reservoir_bc(rho_in, mx_in, my_in, mz_in, rho0=rho0, g=g)
            rho_in = cp.maximum(rho_in, rho_floor)

            P_in, cs_in = eos_polytrope(rho_in, rho0=rho0, P0=P0, n=n, rho_floor=rho_floor)

            # Flux divergence term
            drho_f, dmx_f, dmy_f, dmz_f = rusanov_divergence(
                rho_in, mx_in, my_in, mz_in,
                P_in, cs_in,
                dx=dx, g=g
            )

            # Sources on interior
            sl = slice(g, -g)
            rho_c = rho_in[sl, sl, sl]
            mx_c = mx_in[sl, sl, sl]
            my_c = my_in[sl, sl, sl]
            mz_c = mz_in[sl, sl, sl]
            rho_c_safe = cp.maximum(rho_c, rho_floor)
            vx_c = mx_c / rho_c_safe
            vy_c = my_c / rho_c_safe
            vz_c = mz_c / rho_c_safe

            # sink: total integral removal is Mdot (per time)
            # S_rho_sink integrates to -Mdot
            S_rho_sink = (-Mdot) * W  # units: rho/time (because W integrates to 1)

            # --- Uniform bulk refill: total injection = lam*Mdot (per time) ---
            V_domain = L**3  # physical volume of interior domain (same L used to define dx)
            S_rho_ref = (lam * Mdot / V_domain) * cp.ones_like(rho_c, dtype=cp.float64)
            S_rho = S_rho_sink + S_rho_ref

            # momentum sources
            # sink removes comoving momentum:
            S_mx_sink = S_rho_sink * vx_c
            S_my_sink = S_rho_sink * vy_c
            S_mz_sink = S_rho_sink * vz_c

            if refill_momentum == "comoving":
                S_mx_ref = S_rho_ref * vx_c
                S_my_ref = S_rho_ref * vy_c
                S_mz_ref = S_rho_ref * vz_c
            else:
                S_mx_ref = 0.0
                S_my_ref = 0.0
                S_mz_ref = 0.0

            dmx_s = S_mx_sink + S_mx_ref
            dmy_s = S_my_sink + S_my_ref
            dmz_s = S_mz_sink + S_mz_ref

            # Optional bulk drag to suppress sink-driven inflow
            if nu_drag > 0.0:
                dmx_s = dmx_s - nu_drag * mx_c
                dmy_s = dmy_s - nu_drag * my_c
                dmz_s = dmz_s - nu_drag * mz_c

            return (drho_f + S_rho), (dmx_f + dmx_s), (dmy_f + dmy_s), (dmz_f + dmz_s)

        # --- RK2 (Heun) ---
        sl = slice(g, -g)

        k1_rho, k1_mx, k1_my, k1_mz = rhs(rho, mx, my, mz, t_sim)

        rho1 = rho.copy()
        mx1 = mx.copy()
        my1 = my.copy()
        mz1 = mz.copy()
        rho1[sl, sl, sl] = rho[sl, sl, sl] + dt * k1_rho
        mx1[sl, sl, sl]  = mx [sl, sl, sl] + dt * k1_mx
        my1[sl, sl, sl]  = my [sl, sl, sl] + dt * k1_my
        mz1[sl, sl, sl]  = mz [sl, sl, sl] + dt * k1_mz

        k2_rho, k2_mx, k2_my, k2_mz = rhs(rho1, mx1, my1, mz1, t_sim + dt)

        rho[sl, sl, sl] = rho[sl, sl, sl] + 0.5 * dt * (k1_rho + k2_rho)
        mx [sl, sl, sl] = mx [sl, sl, sl] + 0.5 * dt * (k1_mx  + k2_mx)
        my [sl, sl, sl] = my [sl, sl, sl] + 0.5 * dt * (k1_my  + k2_my)
        mz [sl, sl, sl] = mz [sl, sl, sl] + 0.5 * dt * (k1_mz  + k2_mz)

        # advance time
        t_sim = t_next

        # Diagnostics
        if step % diag_every == 0 or step == steps - 1:
            enforce_reservoir_bc(rho, mx, my, mz, rho0=rho0, g=g)
            rho = cp.maximum(rho, rho_floor)
            P, cs = eos_polytrope(rho, rho0=rho0, P0=P0, n=n, rho_floor=rho_floor)
            vx = mx / cp.maximum(rho, rho_floor)
            vy = my / cp.maximum(rho, rho_floor)
            vz = mz / cp.maximum(rho, rho_floor)

            # interior
            rho_c = rho[sl, sl, sl]
            P_c = P[sl, sl, sl]
            cs_c = cs[sl, sl, sl]
            vx_c = vx[sl, sl, sl]
            vy_c = vy[sl, sl, sl]
            vz_c = vz[sl, sl, sl]

            vmag = cp.sqrt(vx_c * vx_c + vy_c * vy_c + vz_c * vz_c)
            mach_max = float(cp.max(vmag / cp.maximum(cs_c, 1e-12)).get())

            # Spherical averages
            rho_bar = spherical_average(rho_c, bin_idx_flat, counts, nbins)
            rho_bar_np = rho_bar.get().astype(np.float64)

            # Compute pressure from the shell-mean density in float64 (avoids float32 pressure quantization)
            rho_bar_np_safe = np.maximum(rho_bar_np, rho_floor)
            P_bar_np = P0 * np.power(rho_bar_np_safe / rho0, n)

            # --- Reference-subtracted pressure perturbation ---
            # Use the outer shell as the reference "far field" pressure.
            # Pick bins near r_max (e.g., last 10% of bins).
            outer_frac = 0.10
            outer_n = max(5, int(len(r_centers) * outer_frac))

            P_ref = float(np.mean(P_bar_np[-outer_n:]))
            dP = P_bar_np - P_ref

            # Effective acceleration from radial pressure gradient
            dPdr = np.gradient(P_bar_np, r_centers)  # derivative of absolute P, not dP (same derivative anyway)
            g_eff = dPdr / np.maximum(rho_bar_np, 1e-12)

            # Fit windows
            fit_rmin = args.fit_rmin_factor * a
            fit_rmax = args.fit_rmax_frac * L

            # --- set fit cutoffs from the FIT WINDOW only (prevents npts=0 early-time) ---
            mask_fit = (r_centers >= fit_rmin) & (r_centers <= fit_rmax) & np.isfinite(dP) & np.isfinite(g_eff)

            dp_fit_abs   = np.abs(dP[mask_fit])
            geff_fit_abs = np.abs(g_eff[mask_fit])

            dp_peak_fit   = float(np.nanmax(dp_fit_abs))   if dp_fit_abs.size   else 0.0
            geff_peak_fit = float(np.nanmax(geff_fit_abs)) if geff_fit_abs.size else 0.0

            # If the signal is tiny, don't threshold it away.
            ymin_dp   = 0.0 if dp_peak_fit   < 1e-14 else max(1e-14, dp_peak_fit   * 1e-3)
            ymin_geff = 0.0 if geff_peak_fit < 1e-18 else max(1e-18, geff_peak_fit * 1e-3)

            pP, bP, nP = fit_powerlaw(r_centers, dP,    fit_rmin, fit_rmax, y_min_abs=ymin_dp)
            pg, bg, ng = fit_powerlaw(r_centers, g_eff, fit_rmin, fit_rmax, y_min_abs=ymin_geff)

            # Mass budgets (integrated rates)
            # In this model, W is normalized so the sink integral is -Mdot by construction.
            sink_rate = float(-Mdot)
            refill_rate = float(lam * Mdot)

            # Create a compact sample table for copy/paste (downsample bins)
            # Only bins in fit window, sample up to ~40 points
            mask = (r_centers >= fit_rmin) & (r_centers <= fit_rmax)
            idx = np.where(mask)[0]
            if len(idx) > 0:
                stride = max(1, len(idx) // 40)
                idx_s = idx[::stride]
                sample = [
                    {
                        "r": float(r_centers[i]),
                        "dP": float(dP[i]),
                        "geff": float(g_eff[i]),
                        "rho_bar": float(rho_bar_np[i]),
                        "P_bar": float(P_bar_np[i]),
                    }
                    for i in idx_s
                ]
            else:
                sample = []

            out = {
                "event": "diag",
                "step": int(step),
                "t": float(t_sim),
                "dt": float(dt),
                "max_wave": float(max_wave),
                "mach_max": float(mach_max),
                "Mdot": float(Mdot),
                "P_ref": float(P_ref),
                "mass_budget": {"sink_rate": float(sink_rate), "refill_rate": float(refill_rate)},
                "fit_window": {"rmin": float(fit_rmin), "rmax": float(fit_rmax)},
                "fits": {
                    "dP_slope": pP, "dP_npts": nP,
                    "geff_slope": pg, "geff_npts": ng
                },
                "fit_cutoffs": {"ymin_dp": float(ymin_dp), "ymin_geff": float(ymin_geff),
                "dp_peak_fit": float(dp_peak_fit), "geff_peak_fit": float(geff_peak_fit)},
                "sample_points": sample,
                "wall_s": float(time.time() - wall0),
            }
            print(json.dumps(out))

    if args.snapshot_out:
        enforce_reservoir_bc(rho, mx, my, mz, rho0=rho0, g=g)
        rho = cp.maximum(rho, rho_floor)
        dt_rho_snapshot = None
        if steps > 0:
            dt_rho_snapshot, _, _, _ = rhs(rho, mx, my, mz, t_sim)
        save_projected_snapshot(
            args.snapshot_out,
            rho=rho,
            mx=mx,
            my=my,
            mz=mz,
            W=W,
            dt_rho=dt_rho_snapshot,
            N=N,
            L=L,
            dx=dx,
            g=g,
            Mdot=Mdot if steps > 0 else 0.0,
            lam=lam,
            rho0=rho0,
            P0=P0,
            n=n,
            a=a,
            t_sim=t_sim,
            step=max(steps - 1, 0),
        )
        print(json.dumps({"event": "snapshot", "path": args.snapshot_out, "step": int(max(steps - 1, 0)), "t": float(t_sim)}))

    print(json.dumps({"event": "done", "wall_s": float(time.time() - wall0), "t_final": float(t_sim)}))


if __name__ == "__main__":
    main()
