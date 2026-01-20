#!/usr/bin/env python3
"""
stepE_euler_throat_longitudinal.py

Step E: Full (compressible) isothermal Euler + throat sink, analyzed in the
"potential-flow / Poisson" variables that actually matter for the toy model.

Why this exists
---------------
Earlier 3D Euler tests (single_throat_monopole.py) tried to read off an
effective "gravity" from pressure gradients. In a sustained sink problem,
the natural steady response is *through-flow* (potential flow), not a strictly
hydrostatic pressure profile. So the robust Newtonian-sector observable is the
longitudinal (irrotational) flow potential ψ defined by:

    v_L = ∇ψ

In the far-field, low-Mach, nearly-constant-density limit (ρ≈ρ0), sourced
continuity implies:

    ∂t ρ + ∇·(ρ v) = Sρ
    →  ρ0 ∇·v ≈ Sρ
    →  ∇²ψ ≈ Sρ/ρ0

So Poisson emerges as the static limit for ψ.

What this script does
---------------------
1) Evolves 3D *isothermal Euler* on a periodic box using a simple, robust
   Rusanov (LLF) finite-volume flux with RK2 time stepping.

   State variables:
     ρ, m = (mx,my,mz) = ρ v

   EOS (isothermal / constant sound speed):
     P = c_s^2 ρ

2) Adds a single "throat" modeled as a localized sink kernel W(x) with
   throughput Mdot(t):

     Sρ_sink = -Mdot(t) W(x),   with ∑ W dx^3 = 1

   and a uniform refill (bulk leak-in / "dark energy" bookkeeping) so that the
   net source integrates to ~0:

     Sρ_refill = +Mdot(t) / L^3

   Momentum source handling (selectable):
     - comoving: Sm = v Sρ   (preserve local velocity)
     - at_rest : Sm includes sink term v Sρ_sink, but refill adds zero momentum

   Optional linear drag on momentum: Sm_drag = -drag * m

3) Diagnostics focus on the *Poisson / potential-flow sector*:

   - Helmholtz / longitudinal potential ψ from velocity in k-space:
        ψ_k = -(i k·v_k)/k^2        (k=0 mode set to 0)
     (equivalently ψ solves ∇²ψ = ∇·v)

   - Poisson predictor from throat source (static limit):
        ψP_k = -(Sρ/ρ0)_k / k^2     (k=0 mode set to 0)
     (implemented by FFT of W at diag time; uniform refill only affects k=0)

   - "lag" ratio (Parseval-invariant):
        ||ψ - ψP|| / ||ψP||

   - Irrotational energy fraction using k-space projection identities:
        E_long(k) = |k·v_k|^2 / k^2
        E_tot (k) = |vx_k|^2 + |vy_k|^2 + |vz_k|^2
        frac_long = Σ E_long / Σ E_tot

   - Monopole scaling from spherical averages of ψ:
        ψ̄(r) ≈ A/r + C
        g(r) ≡ -dψ̄/dr  ~ 1/r^2

Outputs are JSON lines similar to stepB/C/D.

Notes on GPU memory
-------------------
A 512^3 float64 Euler run is very memory heavy (4+ large real arrays) and FFT
diagnostics add large complex arrays. For N=512, float32 is recommended unless
you have a very large GPU.

Examples
--------
GPU (recommended start):
  python stepE_euler_throat_longitudinal.py --backend cupy --dtype float32 --N 256 --L 200 \
    --Mdot 20 --sigma 2 --cs 7.0710678118654755 --steps 8000 --diag_every 200 --gamma_drag 0.0

GPU (heavier):
  python stepE_euler_throat_longitudinal.py --backend cupy --dtype float32 --N 512 --L 200 \
    --Mdot 20 --sigma 2 --cs 7.0710678118654755 --steps 6000 --diag_every 300

Run:
  python stepE_euler_throat_longitudinal.py \
    --backend cupy --dtype float64 \
    --N 512 --L 200 \
    --Mdot 20 --sigma 2 --ramp_time 2 \
    --refill_momentum at_rest \
    --gamma_drag 0.03 \
    --steps 12000 --diag_every 300


Output: {
    "event": "diag",
    "step": 10200,
    "t": 137.64123698929853,
    "dt": 0.013491859487148584,
    "Mdot": 20.0,
    "max_wave": 7.238160914217986,
    "mach_max": 0.02363053315259165,
    "mass": 7999999.9999999795,
    "mass_rel_drift": -2.561137080192566e-15,
    "mean_momentum": {
        "mx": -1.0587911840678754e-22,
        "my": -3.705769144237564e-22,
        "mz": -4.7356090068660834e-23},
        "fft_metrics": {
            "frac_long": 0.9999952877171262,
            "lag_rms_over_psiP_rms": 0.0316061719169052
        },
        "fits": {
            "dpsi_slope": -1.4854411345744267,
            "dpsi_npts": 76,
            "dpsi_clean_slope": -0.9953501929869432,
            "dpsi_clean_npts": 76,
            "g_slope": -2.044912024020909,
            "g_npts": 76
        },
        "monopole_fit": {
            "A_fit": 1.5767984094900636,
            "C_fit": -0.0222675921357526,
            "A_expected": 1.5915494309189535, 
            "Mdot_from_flux": 19.63040015313316,
            "n_fit": 76
        },
        "refs": {
            "psi_ref": -0.0020271460899637353
        },
        "fit_cutoffs": {
            "ymin_dpsi": 1.7596291693201083e-05,
            "ymin_g": 2.4839767730964108e-06,
            "dpsi_peak_fit": 0.17596291693201083,
            "g_peak_fit": 0.024839767730964107
        }, ...

"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass

import numpy as np


# ----------------------------
# Backend selection
# ----------------------------

@dataclass
class Backend:
    name: str
    xp: object
    fftn: callable
    ifftn: callable
    fftfreq: callable
    asnumpy: callable


def pick_backend(name: str) -> Backend:
    name = name.lower()
    if name not in {"auto", "numpy", "cupy"}:
        raise ValueError("--backend must be one of: auto, numpy, cupy")

    if name in {"auto", "cupy"}:
        try:
            import cupy as cp  # type: ignore

            return Backend(
                name="cupy",
                xp=cp,
                fftn=cp.fft.fftn,
                ifftn=cp.fft.ifftn,
                fftfreq=cp.fft.fftfreq,
                asnumpy=cp.asnumpy,
            )
        except Exception:
            if name == "cupy":
                raise

    return Backend(
        name="numpy",
        xp=np,
        fftn=np.fft.fftn,
        ifftn=np.fft.ifftn,
        fftfreq=np.fft.fftfreq,
        asnumpy=lambda x: x,
    )


def _dtype_from_str(s: str):
    s = s.lower()
    if s == "float32":
        return np.float32
    if s == "float64":
        return np.float64
    raise ValueError("--dtype must be float32 or float64")


# ----------------------------
# Geometry / binning
# ----------------------------

def precompute_radial_bins_cpu(N: int, L: float, nbins: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """CPU radial bin indices for spherical averages, with an overflow bin nbins.

    Returns:
      r_centers: (nbins,) float64
      bin_idx_flat: (N^3,) int32 with values in [0..nbins] (overflow bin=nbins)
      counts: (nbins,) float64 (excludes overflow bin)
      dr: bin width
    """
    dx = L / N
    x = (np.arange(N, dtype=np.float64) + 0.5) * dx - 0.5 * L
    x2 = x * x

    r_max = 0.5 * L
    dr = r_max / nbins

    bin_idx = np.empty((N, N, N), dtype=np.int32)
    for i in range(N):
        r2_plane = x2[i] + x2[:, None] + x2[None, :]
        r_plane = np.sqrt(r2_plane)
        b = np.floor(r_plane / dr).astype(np.int32)
        b = np.where(r_plane < r_max, b, nbins)
        bin_idx[i, :, :] = b

    bin_idx_flat = bin_idx.ravel()
    counts_all = np.bincount(bin_idx_flat, minlength=nbins + 1).astype(np.float64)
    counts = counts_all[:nbins]
    r_centers = (np.arange(nbins, dtype=np.float64) + 0.5) * dr
    return r_centers, bin_idx_flat, counts, dr


def spherical_average(xp, field, bin_idx_flat, counts, nbins: int):
    """Spherical average of a 3D field using precomputed bin indices.

    For CuPy: use cupyx.scatter_add (robust across versions).
    For NumPy: use np.bincount.

    bin_idx_flat must include the overflow bin index (=nbins), so we allocate
    sums_all of length nbins+1 and then discard the overflow slot.
    """
    flat = field.ravel()

    if xp.__name__ == "cupy":
        # cupyx.scatter_add is in-place and may return None depending on CuPy version.
        import cupyx  # type: ignore

        b = bin_idx_flat  # already cupy int32
        sums_all = xp.zeros(nbins + 1, dtype=flat.dtype)
        cupyx.scatter_add(sums_all, b, flat)
        sums = sums_all[:nbins]

        denom = counts
        denom = xp.maximum(denom, 1).astype(flat.dtype)  # avoid divide-by-zero
        avg = sums / denom
        return avg

    # NumPy
    b = np.asarray(bin_idx_flat, dtype=np.int32)
    sums_all = np.bincount(b, weights=np.asarray(flat), minlength=nbins + 1).astype(np.float64)
    sums = sums_all[:nbins]
    denom = np.maximum(np.asarray(counts), 1.0)
    return sums / denom


# ----------------------------
# Fitting helpers
# ----------------------------

def fit_powerlaw(x: np.ndarray, y: np.ndarray, x_min: float, x_max: float, y_min_abs: float = 0.0):
    """Fit log-log slope |y| ~ x^p over window."""
    mask = (x >= x_min) & (x <= x_max) & np.isfinite(y) & (np.abs(y) > y_min_abs)
    xx = x[mask]
    yy = y[mask]
    if xx.size < 6:
        return float("nan"), float("nan"), int(xx.size)
    lx = np.log(xx)
    ly = np.log(np.abs(yy))
    p, b = np.polyfit(lx, ly, 1)
    return float(p), float(b), int(xx.size)


def fit_A_over_r_plus_C(r: np.ndarray, psi: np.ndarray, rmin: float, rmax: float):
    """Least squares fit psi(r) ≈ A/r + C over [rmin,rmax]."""
    mask = (r >= rmin) & (r <= rmax) & np.isfinite(psi)
    rr = r[mask]
    yy = psi[mask]
    if rr.size < 6:
        return float("nan"), float("nan"), int(rr.size)
    X = np.vstack([1.0 / rr, np.ones_like(rr)]).T
    # Solve min ||X [A,C]^T - y||
    coeff, *_ = np.linalg.lstsq(X, yy, rcond=None)
    A, C = coeff[0], coeff[1]
    return float(A), float(C), int(rr.size)


# ----------------------------
# Throat kernel and source
# ----------------------------

def make_gaussian_kernel_W(xp, N: int, L: float, dx: float, sigma: float, dtype):
    """Return W(x) >=0 normalized so that sum(W)*dx^3 = 1."""
    x = (xp.arange(N, dtype=dtype) + 0.5) * dx - 0.5 * L
    x2 = x * x
    r2 = x2[:, None, None] + x2[None, :, None] + x2[None, None, :]
    W = xp.exp(-r2 / (sigma * sigma)).astype(dtype)
    norm = xp.sum(W) * (dx**3)
    W = W / norm
    return W


def smoothstep01(s: float) -> float:
    """C^1 smooth ramp from 0 to 1 on [0,1]."""
    if s <= 0.0:
        return 0.0
    if s >= 1.0:
        return 1.0
    return s * s * (3.0 - 2.0 * s)


# ----------------------------
# Euler RHS (isothermal) with Rusanov flux
# ----------------------------

def flux_divergence_isothermal(xp, rho, mx, my, mz, cs2: float, dx: float):
    """Compute -∇·F(U) for isothermal Euler using Rusanov flux in each direction.

    Returns:
      drho_dt, dmx_dt, dmy_dt, dmz_dt  (all same shape as rho)
    """
    # Protect against tiny/negative density
    rho_safe = rho
    vx = mx / rho_safe
    vy = my / rho_safe
    vz = mz / rho_safe
    P = cs2 * rho_safe

    drho = xp.zeros_like(rho)
    dmx_ = xp.zeros_like(rho)
    dmy_ = xp.zeros_like(rho)
    dmz_ = xp.zeros_like(rho)

    # Helper to add one-direction contribution
    def add_dir(axis: int):
        nonlocal drho, dmx_, dmy_, dmz_, rho, mx, my, mz, vx, vy, vz, P

        # Right states via periodic shift
        rhoR = xp.roll(rho, -1, axis=axis)
        mxR  = xp.roll(mx,  -1, axis=axis)
        myR  = xp.roll(my,  -1, axis=axis)
        mzR  = xp.roll(mz,  -1, axis=axis)

        vxR = xp.roll(vx, -1, axis=axis)
        vyR = xp.roll(vy, -1, axis=axis)
        vzR = xp.roll(vz, -1, axis=axis)
        PR  = xp.roll(P,  -1, axis=axis)

        # Normal velocities
        if axis == 0:
            vnL = vx
            vnR = vxR
        elif axis == 1:
            vnL = vy
            vnR = vyR
        else:
            vnL = vz
            vnR = vzR

        a = xp.maximum(xp.abs(vnL), xp.abs(vnR)) + math.sqrt(cs2)

        # Fluxes (L and R) for each conserved component
        if axis == 0:
            FrhoL = mx
            FrhoR = mxR
            FmxL = mx * vx + P
            FmxR = mxR * vxR + PR
            FmyL = my * vx
            FmyR = myR * vxR
            FmzL = mz * vx
            FmzR = mzR * vxR
        elif axis == 1:
            FrhoL = my
            FrhoR = myR
            FmxL = mx * vy
            FmxR = mxR * vyR
            FmyL = my * vy + P
            FmyR = myR * vyR + PR
            FmzL = mz * vy
            FmzR = mzR * vyR
        else:
            FrhoL = mz
            FrhoR = mzR
            FmxL = mx * vz
            FmxR = mxR * vzR
            FmyL = my * vz
            FmyR = myR * vzR
            FmzL = mz * vz + P
            FmzR = mzR * vzR + PR

        # Rusanov interface fluxes
        Frho = 0.5 * (FrhoL + FrhoR) - 0.5 * a * (rhoR - rho)
        Fmx  = 0.5 * (FmxL  + FmxR)  - 0.5 * a * (mxR  - mx)
        Fmy  = 0.5 * (FmyL  + FmyR)  - 0.5 * a * (myR  - my)
        Fmz  = 0.5 * (FmzL  + FmzR)  - 0.5 * a * (mzR  - mz)

        # Divergence contribution: (F_{i+1/2} - F_{i-1/2})/dx
        FrhoLft = xp.roll(Frho, 1, axis=axis)
        FmxLft  = xp.roll(Fmx,  1, axis=axis)
        FmyLft  = xp.roll(Fmy,  1, axis=axis)
        FmzLft  = xp.roll(Fmz,  1, axis=axis)

        drho -= (Frho - FrhoLft) / dx
        dmx_ -= (Fmx  - FmxLft)  / dx
        dmy_ -= (Fmy  - FmyLft)  / dx
        dmz_ -= (Fmz  - FmzLft)  / dx

    add_dir(0)
    add_dir(1)
    add_dir(2)

    return drho, dmx_, dmy_, dmz_


# ----------------------------
# k-space tools for ψ and irrotational fraction
# ----------------------------

def build_k_arrays(xp, fftfreq, N: int, dx: float, dtype):
    """Return 1D k array (2π * fftfreq) and inv_k2 3D array with k0 handled."""
    k1 = (2.0 * math.pi) * fftfreq(N, d=dx)
    k1 = k1.astype(dtype)
    kx = k1[:, None, None]
    ky = k1[None, :, None]
    kz = k1[None, None, :]
    k2 = (kx * kx + ky * ky + kz * kz).astype(dtype)
    # inv_k2 with k=0 set to 0 to avoid NaNs
    inv_k2 = xp.where(k2 > 0, 1.0 / k2, xp.array(0.0, dtype=dtype))
    return k1, inv_k2


# ----------------------------
# Main
# ----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", choices=["auto", "numpy", "cupy"], default="auto")
    ap.add_argument("--dtype", choices=["float32", "float64"], default="float32")

    ap.add_argument("--N", type=int, default=256)
    ap.add_argument("--L", type=float, default=200.0)

    ap.add_argument("--rho0", type=float, default=1.0)
    ap.add_argument("--cs", type=float, default=7.0710678118654755, help="Sound speed for isothermal EOS.")

    ap.add_argument("--sigma", type=float, default=2.0)
    ap.add_argument("--Mdot", type=float, default=20.0)
    ap.add_argument("--ramp_time", type=float, default=2.0)

    ap.add_argument("--refill_momentum", choices=["comoving", "at_rest"], default="comoving")
    ap.add_argument("--gamma_drag", type=float, default=0.0, help="Linear drag coefficient on momentum (0 disables).")
    ap.add_argument("--rho_floor", type=float, default=1e-6)

    ap.add_argument("--steps", type=int, default=8000)
    ap.add_argument("--diag_every", type=int, default=200)

    ap.add_argument("--cfl", type=float, default=0.25)
    ap.add_argument("--dt", type=float, default=None, help="Fixed dt (overrides CFL if set).")
    ap.add_argument("--dt_max", type=float, default=None, help="Optional cap on dt if CFL is used.")

    ap.add_argument("--remove_mean_momentum", action="store_true", help="Remove global mean momentum each step.")

    ap.add_argument("--nbins", type=int, default=180)
    ap.add_argument("--fit_rmin", type=float, default=8.0)
    ap.add_argument("--fit_rmax", type=float, default=50.0)

    ap.add_argument("--sample_points", type=int, default=12, help="How many radial samples to print in diagnostics.")
    args = ap.parse_args()

    be = pick_backend(args.backend)
    xp = be.xp

    dtype_np = _dtype_from_str(args.dtype)
    dtype = xp.float32 if dtype_np == np.float32 else xp.float64

    N = int(args.N)
    L = float(args.L)
    dx = L / N

    rho0 = float(args.rho0)
    cs = float(args.cs)
    cs2 = cs * cs

    sigma = float(args.sigma)
    Mdot_target = float(args.Mdot)
    ramp_time = float(args.ramp_time)

    steps = int(args.steps)
    diag_every = int(args.diag_every)

    nbins = int(args.nbins)
    fit_rmin = float(args.fit_rmin)
    fit_rmax = float(args.fit_rmax)

    # Precompute bins on CPU (large for N=512, but avoids runtime coordinate grids)
    r_centers_cpu, bin_idx_flat_cpu, counts_cpu, dr = precompute_radial_bins_cpu(N=N, L=L, nbins=nbins)

    # Move bin structures to backend
    bin_idx_flat = xp.asarray(bin_idx_flat_cpu.astype(np.int32))
    counts = xp.asarray(counts_cpu.astype(np.float64 if dtype_np == np.float64 else np.float32))

    # Precompute throat kernel W
    W = make_gaussian_kernel_W(xp=xp, N=N, L=L, dx=dx, sigma=sigma, dtype=dtype)

    # k arrays for FFT diagnostics (ψ and lag)
    k1, inv_k2 = build_k_arrays(xp, be.fftfreq, N=N, dx=dx, dtype=dtype)

    # Initial state (uniform density, zero momentum)
    rho = (rho0 * xp.ones((N, N, N), dtype=dtype))
    mx = xp.zeros_like(rho)
    my = xp.zeros_like(rho)
    mz = xp.zeros_like(rho)

    # Mass reference
    vol = L**3
    mass0 = float(be.asnumpy(xp.sum(rho) * (dx**3)))

    # Fixed dt if provided, else CFL
    dt_fixed = None if args.dt is None else float(args.dt)
    dt_max = None if args.dt_max is None else float(args.dt_max)

    wall0 = time.time()

    init = {
        "event": "init",
        "backend": be.name,
        "N": N,
        "L": L,
        "dx": dx,
        "dtype": str(np.dtype(dtype_np)),
        "params": {
            "rho0": rho0,
            "cs": cs,
            "gamma_drag": float(args.gamma_drag),
            "cfl": float(args.cfl),
            "dt": dt_fixed,
            "dt_max": dt_max,
            "rho_floor": float(args.rho_floor),
        },
        "throat": {"Mdot_target": Mdot_target, "sigma": sigma, "ramp_time": ramp_time},
        "refill": {"momentum_mode": args.refill_momentum, "method": "uniform + mean-zero periodic solvability"},
        "domain": {"bc": "periodic"},
        "diag": {"nbins": nbins, "fit_rmin": fit_rmin, "fit_rmax": fit_rmax},
        "targets": {"dpsi_clean_slope": -1.0, "g_slope": -2.0, "lag_rms_over_psiP_rms": 0.0, "frac_long": 1.0},
        "notes": "StepE: evolve isothermal Euler; diagnose longitudinal ψ and compare to Poisson predictor from Sρ.",
    }
    print(json.dumps(init))

    def compute_dt(mx_, my_, mz_, rho_):
        if dt_fixed is not None:
            return dt_fixed, float("nan"), float("nan")
        # velocity magnitude max
        rho_safe = xp.maximum(rho_, xp.array(float(args.rho_floor), dtype=dtype))
        vx = mx_ / rho_safe
        vy = my_ / rho_safe
        vz = mz_ / rho_safe
        vmag = xp.sqrt(vx * vx + vy * vy + vz * vz)
        vmax = float(be.asnumpy(xp.max(vmag)))
        max_wave = cs + vmax
        dt = float(args.cfl) * dx / max(max_wave, 1e-30)
        if dt_max is not None:
            dt = min(dt, dt_max)
        mach_max = vmax / max(cs, 1e-30)
        return dt, max_wave, mach_max

    def source_terms(t: float, rho_, mx_, my_, mz_):
        # Throughput ramp
        s = t / max(ramp_time, 1e-30)
        Mdot = Mdot_target * smoothstep01(s)

        # Density source: sink + uniform refill (net zero)
        S_sink = (-Mdot) * W
        S_ref = (Mdot / vol) * xp.ones_like(rho_)
        S_rho = S_sink + S_ref

        # Momentum source
        rho_safe = xp.maximum(rho_, xp.array(float(args.rho_floor), dtype=dtype))
        vx = mx_ / rho_safe
        vy = my_ / rho_safe
        vz = mz_ / rho_safe

        if args.refill_momentum == "comoving":
            # preserve v under both sink and refill
            S_mx = vx * S_rho
            S_my = vy * S_rho
            S_mz = vz * S_rho
        else:
            # preserve v under sink only; refill adds zero momentum
            S_mx = vx * S_sink
            S_my = vy * S_sink
            S_mz = vz * S_sink

        # Optional linear drag on momentum
        if args.gamma_drag != 0.0:
            S_mx = S_mx - float(args.gamma_drag) * mx_
            S_my = S_my - float(args.gamma_drag) * my_
            S_mz = S_mz - float(args.gamma_drag) * mz_

        return Mdot, S_rho, S_mx, S_my, S_mz

    def rhs(t: float, rho_, mx_, my_, mz_):
        # density floor in state for stability
        rho_safe = xp.maximum(rho_, xp.array(float(args.rho_floor), dtype=dtype))
        drho, dmx_, dmy_, dmz_ = flux_divergence_isothermal(
            xp=xp, rho=rho_safe, mx=mx_, my=my_, mz=mz_, cs2=cs2, dx=dx
        )
        Mdot, S_rho, S_mx, S_my, S_mz = source_terms(t, rho_safe, mx_, my_, mz_)
        drho = drho + S_rho
        dmx_ = dmx_ + S_mx
        dmy_ = dmy_ + S_my
        dmz_ = dmz_ + S_mz
        return Mdot, drho, dmx_, dmy_, dmz_

    def diag(step: int, t: float, dt: float, Mdot: float, max_wave: float, mach_max: float):
        # Mass / momentum bookkeeping
        mass = float(be.asnumpy(xp.sum(rho) * (dx**3)))
        mass_drift = (mass - mass0) / max(mass0, 1e-30)

        pmean = float(be.asnumpy(xp.mean(mx))) , float(be.asnumpy(xp.mean(my))) , float(be.asnumpy(xp.mean(mz)))

        # Build velocities
        rho_safe = xp.maximum(rho, xp.array(float(args.rho_floor), dtype=dtype))
        vx = mx / rho_safe
        vy = my / rho_safe
        vz = mz / rho_safe

        # --- FFT diagnostics: ψ from v, ψP from source ---
        # We compute dot = k·v_k by accumulating component FFTs.
        # Also accumulate total energy in k-space.
        dot_k = None
        Etot = 0.0

        # Pre-broadcast k components
        kx = k1[:, None, None]
        ky = k1[None, :, None]
        kz = k1[None, None, :]

        # vx
        vx_k = be.fftn(vx)
        if dot_k is None:
            dot_k = kx * vx_k
        else:
            dot_k = dot_k + kx * vx_k
        Etot = Etot + float(be.asnumpy(xp.sum(xp.abs(vx_k) ** 2)))

        # vy
        vy_k = be.fftn(vy)
        dot_k = dot_k + ky * vy_k
        Etot = Etot + float(be.asnumpy(xp.sum(xp.abs(vy_k) ** 2)))

        # vz
        vz_k = be.fftn(vz)
        dot_k = dot_k + kz * vz_k
        Etot = Etot + float(be.asnumpy(xp.sum(xp.abs(vz_k) ** 2)))

        # Longitudinal energy: |k·v_k|^2 / k^2  (skip k=0 where inv_k2=0)
        Elong = float(be.asnumpy(xp.sum((xp.abs(dot_k) ** 2) * inv_k2)))
        frac_long = Elong / max(Etot, 1e-30)

        # ψ_k from velocity: ψ_k = -(i k·v_k)/k^2 = -(1j * dot_k) * inv_k2
        psi_k = (-(1j) * dot_k) * inv_k2
        # k=0 handled by inv_k2=0

        # ψP_k from source: for k!=0, only sink contributes (uniform refill is k=0)
        # Compute W_k on the fly (one FFT) and scale.
        W_k = be.fftn(W)
        psiP_k = (Mdot / rho0) * (W_k * inv_k2)

        # lag ratio in k-space (Parseval invariant up to a constant)
        # Exclude k=0 automatically because inv_k2=0 ⇒ psiP_k[0]=0.
        num = float(be.asnumpy(xp.sum(xp.abs(psi_k - psiP_k) ** 2)))
        den = float(be.asnumpy(xp.sum(xp.abs(psiP_k) ** 2)))
        lag_ratio = math.sqrt(num / max(den, 1e-30))

        # Real-space ψ for radial profiles
        psi = be.ifftn(psi_k).real.astype(dtype, copy=False)

        # Spherical averages
        psi_bar = spherical_average(xp, psi, bin_idx_flat, counts, nbins)
        psi_bar_np = be.asnumpy(psi_bar).astype(np.float64)

        # Reference subtract (outer shells)
        outer_n = max(8, int(0.1 * nbins))
        psi_ref = float(np.mean(psi_bar_np[-outer_n:]))
        dpsi = psi_bar_np - psi_ref

        # g(r) = -dψ/dr (positive for sink monopole)
        g = -np.gradient(psi_bar_np, r_centers_cpu)

        # Fit cutoffs (avoid fitting numerical noise)
        mask_fit = (r_centers_cpu >= fit_rmin) & (r_centers_cpu <= fit_rmax)
        dpsi_peak = float(np.nanmax(np.abs(dpsi[mask_fit]))) if np.any(mask_fit) else 0.0
        g_peak = float(np.nanmax(np.abs(g[mask_fit]))) if np.any(mask_fit) else 0.0
        ymin_dpsi = 0.0 if dpsi_peak < 1e-20 else max(1e-20, dpsi_peak * 1e-4)
        ymin_g = 0.0 if g_peak < 1e-24 else max(1e-24, g_peak * 1e-4)

        dpsi_slope, _, dpsi_n = fit_powerlaw(r_centers_cpu, dpsi, fit_rmin, fit_rmax, y_min_abs=ymin_dpsi)
        g_slope, _, g_n = fit_powerlaw(r_centers_cpu, g, fit_rmin, fit_rmax, y_min_abs=ymin_g)

        # Monopole fit ψ̄ ≈ A/r + C
        A_fit, C_fit, n_fit = fit_A_over_r_plus_C(r_centers_cpu, psi_bar_np, fit_rmin, fit_rmax)

        dpsi_clean = psi_bar_np - C_fit
        dpsi_clean_slope, _, dpsi_clean_n = fit_powerlaw(
            r_centers_cpu, dpsi_clean, fit_rmin, fit_rmax, y_min_abs=ymin_dpsi
        )

        # Flux-based throughput estimate from g(r): Mdot ≈ 4π ρ0 r^2 g
        if np.any(mask_fit):
            mdot_prof = (4.0 * math.pi * rho0) * (r_centers_cpu[mask_fit] ** 2) * g[mask_fit]
            mdot_est = float(np.median(mdot_prof[np.isfinite(mdot_prof)])) if mdot_prof.size > 0 else float("nan")
        else:
            mdot_est = float("nan")

        # Optional sample points
        nsamp = max(0, int(args.sample_points))
        sample = []
        if nsamp > 0:
            # sample evenly in the fit window
            idx = np.where(mask_fit)[0]
            if idx.size > 0:
                pick = np.linspace(0, idx.size - 1, num=min(nsamp, idx.size)).astype(int)
                for j in idx[pick]:
                    sample.append({
                        "r": float(r_centers_cpu[j]),
                        "psi_bar": float(psi_bar_np[j]),
                        "dpsi": float(dpsi[j]),
                        "g": float(g[j]),
                    })

        out = {
            "event": "diag",
            "step": int(step),
            "t": float(t),
            "dt": float(dt),
            "Mdot": float(Mdot),
            "max_wave": float(max_wave) if math.isfinite(max_wave) else None,
            "mach_max": float(mach_max) if math.isfinite(mach_max) else None,
            "mass": mass,
            "mass_rel_drift": float(mass_drift),
            "mean_momentum": {"mx": float(pmean[0]), "my": float(pmean[1]), "mz": float(pmean[2])},
            "fft_metrics": {
                "frac_long": float(frac_long),
                "lag_rms_over_psiP_rms": float(lag_ratio),
            },
            "fits": {
                "dpsi_slope": float(dpsi_slope),
                "dpsi_npts": int(dpsi_n),
                "dpsi_clean_slope": float(dpsi_clean_slope),
                "dpsi_clean_npts": int(dpsi_clean_n),
                "g_slope": float(g_slope),
                "g_npts": int(g_n),
            },
            "monopole_fit": {
                "A_fit": float(A_fit),
                "C_fit": float(C_fit),
                "A_expected": float(Mdot / (4.0 * math.pi * rho0)),
                "Mdot_from_flux": float(mdot_est),
                "n_fit": int(n_fit),
            },
            "refs": {"psi_ref": float(psi_ref)},
            "fit_cutoffs": {"ymin_dpsi": float(ymin_dpsi), "ymin_g": float(ymin_g), "dpsi_peak_fit": dpsi_peak, "g_peak_fit": g_peak},
            "sample_points": sample,
            "wall_s": float(time.time() - wall0),
        }
        print(json.dumps(out))

    # Initial diagnostic
    t_sim = 0.0
    dt0, maxw0, mach0 = compute_dt(mx, my, mz, rho)
    Mdot0, *_ = source_terms(t_sim, rho, mx, my, mz)
    diag(step=0, t=t_sim, dt=dt0, Mdot=Mdot0, max_wave=maxw0, mach_max=mach0)

    # Time loop (RK2 / Heun)
    for step in range(1, steps + 1):
        dt, max_wave, mach_max = compute_dt(mx, my, mz, rho)

        Mdot, drho1, dmx1, dmy1, dmz1 = rhs(t_sim, rho, mx, my, mz)

        rho1 = rho + dt * drho1
        mx1 = mx + dt * dmx1
        my1 = my + dt * dmy1
        mz1 = mz + dt * dmz1

        # Density floor after stage
        rho1 = xp.maximum(rho1, xp.array(float(args.rho_floor), dtype=dtype))

        Mdot2, drho2, dmx2, dmy2, dmz2 = rhs(t_sim + dt, rho1, mx1, my1, mz1)

        rho = rho + 0.5 * dt * (drho1 + drho2)
        mx  = mx  + 0.5 * dt * (dmx1  + dmx2)
        my  = my  + 0.5 * dt * (dmy1  + dmy2)
        mz  = mz  + 0.5 * dt * (dmz1  + dmz2)

        # Density floor
        rho = xp.maximum(rho, xp.array(float(args.rho_floor), dtype=dtype))

        # Optional remove mean momentum (Galilean drift control)
        if args.remove_mean_momentum:
            mx = mx - xp.mean(mx)
            my = my - xp.mean(my)
            mz = mz - xp.mean(mz)

        t_sim += dt

        if (step % diag_every) == 0 or step == steps:
            diag(step=step, t=t_sim, dt=dt, Mdot=Mdot, max_wave=max_wave, mach_max=mach_max)

    print(json.dumps({"event": "done", "t_final": float(t_sim), "wall_s": float(time.time() - wall0), "backend": be.name}))


if __name__ == "__main__":
    main()
