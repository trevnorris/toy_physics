#!/usr/bin/env python3
"""stepB_scalar_poisson_lag.py

Step B: Scalar-field Newtonian sector regression test for the toy model.

This script implements exactly the scalar-field structure used in
`1pn_orbital_dynamics.tex`:

  (1) Instantaneous Poisson sector (0PN / Newtonian):
      ∇² Φ_P(x,t) = 4π G ρ(x,t)

  (2) Finite-speed sector via a sourced wave equation for the *total* potential:
      ∂²_t Φ(x,t) = S² [ ∇² Φ(x,t) - 4π G ρ(x,t) ]

and then defines the lag field diagnostically as:
      Φ_L ≡ Φ - Φ_P.

Why there is a damping knob here
-------------------------------
The analytic “Static Limit Theorem” assumes *outgoing / radiative* boundary
conditions so that any free wave content in Φ_L can leave the domain.

On a periodic FFT box, free waves do not leave — they become standing waves.
To emulate the “radiate away” condition while keeping the solver simple and
reviewable, we optionally add a *small linear absorber*:

      ∂²_t Φ + 2γ ∂_t Φ = S² [ ∇² Φ - 4π G ρ ]

Setting γ>0 causes Φ_L to decay, leaving the Poisson solution as the attractor.
If you set --gamma 0, Φ_L generally persists (as it should in a perfectly
lossless periodic box).

What this script measures
-------------------------
1) Static-limit decay:
     RMS(Φ - Φ_P) / RMS(Φ_P)  vs time

2) Monopole scaling in the weak/far field (inside an interior fit window):
     |ΔΦ(r)| ~ r^{-1}   and   |g(r)| = |∂_r Φ| ~ r^{-2}

Where “ΔΦ(r)” is reference-subtracted spherical average of Φ (constant offsets
are physically irrelevant in Newtonian gravity).

Backend
-------
- If CuPy is available, the script runs on GPU (FFT Poisson + FFT Laplacian).
- Otherwise it runs on CPU using NumPy.

Examples
--------
CPU:
  python stepB_scalar_poisson_lag.py --backend numpy --N 128 --steps 1500

GPU:
  python stepB_scalar_poisson_lag.py --backend cupy --N 256 --steps 2500 --gamma 0.3

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

    # Fallback
    return Backend(
        name="numpy",
        xp=np,
        fftn=np.fft.fftn,
        ifftn=np.fft.ifftn,
        fftfreq=np.fft.fftfreq,
        asnumpy=lambda x: x,
    )


# ----------------------------
# Helpers
# ----------------------------


def _dtype_from_str(s: str):
    s = s.lower()
    if s == "float32":
        return np.float32
    if s == "float64":
        return np.float64
    raise ValueError("--dtype must be float32 or float64")


def estimate_memory_bytes(N: int, dtype: np.dtype, include_poisson_temporaries: bool = True) -> dict:
    """Rough peak-memory estimate for this script.

    Notes:
    - FFT libraries allocate additional work buffers we cannot know exactly.
    - This estimate is meant to help you choose a GPU size conservatively.
    """
    N3 = int(N) ** 3
    item = np.dtype(dtype).itemsize

    bytes_real = N3 * item
    # complex has 2 reals
    bytes_cplx = N3 * 2 * item

    # Real arrays we keep resident during evolution:
    #   rho_src, phiP, phi, pi, k2
    n_real_resident = 5

    # Complex arrays:
    #   phi_k (during evolution)
    n_cplx_resident = 1

    # During the initial Poisson solve we briefly have rho_k and phiP_k.
    n_cplx_poisson = 2 if include_poisson_temporaries else 0

    est_peak = n_real_resident * bytes_real + (n_cplx_resident + n_cplx_poisson) * bytes_cplx

    return {
        "N": int(N),
        "dtype": str(np.dtype(dtype)),
        "bytes_real_array": int(bytes_real),
        "bytes_complex_array": int(bytes_cplx),
        "resident_real_arrays": int(n_real_resident),
        "resident_complex_arrays": int(n_cplx_resident),
        "poisson_complex_temporaries": int(n_cplx_poisson),
        "estimated_peak_bytes": int(est_peak),
        "estimated_peak_gib": float(est_peak) / (1024.0**3),
        "note": "FFT work buffers not included; for safety, budget ~1.5× to 2× this estimate on GPU.",
    }


def make_gaussian_defect_rho(xp, N: int, L: float, dx: float, sigma: float, mass: float, dtype):
    """Create a compact, smooth mass density rho(x) with total mass=mass.

    rho(x) = mass * W(x), where W is a normalized Gaussian kernel:
      sum(W) * dx^3 = 1.

    The defect is centered at the box center.
    """
    # 1D coordinates centered at 0
    x = (xp.arange(N, dtype=dtype) + 0.5) * dx - 0.5 * L
    x2 = x * x

    # r^2 via separability (broadcasting)
    r2 = x2[:, None, None] + x2[None, :, None] + x2[None, None, :]
    W = xp.exp(-r2 / (sigma * sigma)).astype(dtype)

    # Normalize: sum(W)*dx^3 = 1
    norm = xp.sum(W) * (dx**3)
    W = W / norm

    rho = (mass * W).astype(dtype)
    return rho


def build_k2(xp, fftfreq, N: int, L: float, dtype):
    """Return spectral k^2 grid for periodic box of side L."""
    dx = L / N
    k1 = (2.0 * math.pi) * fftfreq(N, d=dx)
    k1 = k1.astype(dtype)
    kx2 = (k1 * k1)[:, None, None]
    ky2 = (k1 * k1)[None, :, None]
    kz2 = (k1 * k1)[None, None, :]
    k2 = (kx2 + ky2 + kz2).astype(dtype)
    return k2


def poisson_solve_fft(xp, fftn, ifftn, rho_src, k2, G: float):
    """Solve ∇^2 phi = 4πG rho_src on a periodic domain, with k=0 set to 0."""
    four_pi_G = (4.0 * math.pi * G)
    rho_k = fftn(rho_src)

    # phi_k = -4πG rho_k / k^2   (since FFT(∇^2) = -k^2)
    # Handle k=0 mode explicitly.
    phi_k = (-four_pi_G) * rho_k / k2
    phi_k = phi_k.astype(rho_k.dtype, copy=False)
    phi_k[(0, 0, 0)] = 0.0

    phi = ifftn(phi_k).real
    return phi


def laplacian_fft(xp, fftn, ifftn, field, k2):
    """Compute periodic Laplacian via FFT: ∇^2 f = IFFT(-k^2 FFT(f))."""
    f_k = fftn(field)
    lap_k = (-k2) * f_k
    lap = ifftn(lap_k).real
    return lap


def precompute_radial_bins_cpu(N: int, L: float, nbins: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Precompute bin indices on CPU to keep code simple and memory predictable.

    Returns:
      r_centers: (nbins,) float64
      bin_idx_flat: (N^3,) int32  (bins 0..nbins-1, overflow bin=nbins)
      counts: (nbins,) float64
      dr: bin width

    Binning convention: inscribed sphere of radius r_max=L/2.
    Cells with r>=r_max go into the overflow bin (discarded in averages).
    """
    dx = L / N
    x = (np.arange(N, dtype=np.float64) + 0.5) * dx - 0.5 * L
    x2 = x * x

    r_max = 0.5 * L
    dr = r_max / nbins

    bin_idx = np.empty((N, N, N), dtype=np.int32)

    # Fill in x-slabs to avoid allocating X,Y,Z grids.
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
    """Spherical shell average via bincount."""
    f = field.ravel().astype(xp.float64)
    sums_all = xp.bincount(bin_idx_flat, weights=f, minlength=nbins + 1).astype(xp.float64)
    sums = sums_all[:nbins]
    avg = sums / xp.maximum(counts, 1e-30)
    return avg


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


# ----------------------------
# Main
# ----------------------------


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", choices=["auto", "numpy", "cupy"], default="auto")
    ap.add_argument("--dtype", choices=["float32", "float64"], default="float32")

    ap.add_argument("--N", type=int, default=192, help="Grid resolution per axis.")
    ap.add_argument("--L", type=float, default=200.0, help="Domain size (periodic box).")

    ap.add_argument("--G", type=float, default=1.0, help="Effective gravitational constant.")
    ap.add_argument("--S", type=float, default=10.0, help="Scalar propagation speed (sound speed analogue).")

    ap.add_argument("--mass", type=float, default=1.0, help="Total defect mass M = ∫rho dV.")
    ap.add_argument("--sigma", type=float, default=2.0, help="Gaussian width of defect density profile.")

    ap.add_argument("--gamma", type=float, default=0.3, help="Damping rate for wave absorber (0 disables).")
    ap.add_argument(
        "--cfl",
        type=float,
        default=0.25,
        help=(
            "CFL for wave equation time step (dt = cfl*dx/S). "
            "For a spectral Laplacian, a conservative stability guide is dt < 2*dx/(π S)."
        ),
    )
    ap.add_argument("--steps", type=int, default=2500)
    ap.add_argument("--diag_every", type=int, default=100)

    ap.add_argument("--nbins", type=int, default=180, help="Radial bins for spherical averages.")
    ap.add_argument("--fit_rmin", type=float, default=None, help="Fit window inner radius (default 4*sigma).")
    ap.add_argument("--fit_rmax", type=float, default=None, help="Fit window outer radius (default 0.25*L).")

    ap.add_argument("--seed", type=int, default=0, help="Random seed for initial lag noise.")
    ap.add_argument("--init", choices=["poisson_plus_noise", "zero_plus_noise"], default="poisson_plus_noise")
    args = ap.parse_args()

    be = pick_backend(args.backend)
    xp = be.xp

    dtype_np = _dtype_from_str(args.dtype)
    dtype = xp.float32 if dtype_np == np.float32 else xp.float64

    N = int(args.N)
    L = float(args.L)
    dx = L / N

    G = float(args.G)
    S = float(args.S)
    gamma = float(args.gamma)

    steps = int(args.steps)
    diag_every = int(args.diag_every)

    nbins = int(args.nbins)

    fit_rmin = float(args.fit_rmin) if args.fit_rmin is not None else 4.0 * float(args.sigma)
    fit_rmax = float(args.fit_rmax) if args.fit_rmax is not None else 0.25 * L

    # Memory estimate (informational)
    mem = estimate_memory_bytes(N=N, dtype=dtype_np)

    # Precompute bins on CPU for stable memory and reviewability
    r_centers_cpu, bin_idx_flat_cpu, counts_cpu, dr = precompute_radial_bins_cpu(N=N, L=L, nbins=nbins)

    # Move bin structures to backend
    bin_idx_flat = xp.asarray(bin_idx_flat_cpu)
    counts = xp.asarray(counts_cpu)

    # Build spectral grid k^2
    k2 = build_k2(xp, be.fftfreq, N=N, L=L, dtype=dtype)

    # Avoid division by zero at k=0 (handled explicitly in Poisson; here just prevent NaNs)
    k2 = xp.where(k2 == 0, xp.array(1.0, dtype=dtype), k2)

    # Build static defect density and mean-subtracted source (periodic solvability)
    rho = make_gaussian_defect_rho(
        xp=xp, N=N, L=L, dx=dx,
        sigma=float(args.sigma), mass=float(args.mass),
        dtype=dtype
    )

    rho_mean = xp.mean(rho)
    rho_src = rho - rho_mean

    # Solve Poisson once for reference
    phiP = poisson_solve_fft(xp, be.fftn, be.ifftn, rho_src=rho_src, k2=k2, G=G).astype(dtype)

    # Initialize total potential and its time derivative
    rng = np.random.default_rng(int(args.seed))
    noise = rng.standard_normal(size=(N, N, N)).astype(dtype_np)
    noise = noise / max(1.0, float(np.std(noise)))
    noise_amp = 0.05  # relative amplitude vs |phiP| (will be rescaled below)

    noise_xp = xp.asarray(noise)

    # Rescale noise to be ~noise_amp * RMS(phiP)
    phiP_rms = float(be.asnumpy(xp.sqrt(xp.mean(phiP.astype(xp.float64) ** 2))))
    if phiP_rms <= 0:
        phiP_rms = 1.0
    noise_xp = (noise_amp * phiP_rms) * noise_xp

    if args.init == "poisson_plus_noise":
        phi = phiP + noise_xp
    else:
        phi = noise_xp

    # Random initial phi_t as well (small)
    pi = 0.0 * phi + 0.0

    # Wave time step.
    # For a spectral Laplacian on a periodic grid, a conservative explicit stability
    # guide is: dt < 2/(S*k_max) where k_max ≈ π/dx  ⇒  dt < 2*dx/(π S).
    dt = float(args.cfl) * dx / max(S, 1e-30)
    dt_max_est = (2.0 * dx) / (math.pi * max(S, 1e-30))
    if dt > 0.95 * dt_max_est:
        print(json.dumps({
            "event": "warning",
            "kind": "dt_near_stability_limit",
            "dt": float(dt),
            "dt_max_est": float(dt_max_est),
            "hint": "Reduce --cfl or increase N if you see blow-ups / NaNs."
        }))

    # Precompute forcing term
    four_pi_G = (4.0 * math.pi * G)

    wall0 = time.time()

    print(json.dumps({
        "event": "init",
        "backend": be.name,
        "N": N,
        "L": L,
        "dx": dx,
        "dtype": str(np.dtype(dtype_np)),
        "params": {"G": G, "S": S, "gamma": gamma, "cfl": float(args.cfl), "dt": dt, "dt_max_est": float(dt_max_est)},
        "defect": {"mass": float(args.mass), "sigma": float(args.sigma)},
        "source": {"rho_mean": float(be.asnumpy(rho_mean)), "note": "rho_src=rho-<rho> for periodic Poisson solvability"},
        "diag": {"nbins": nbins, "fit_rmin": fit_rmin, "fit_rmax": fit_rmax},
        "memory_estimate": mem,
        "targets": {"dPhi_slope": -1.0, "g_slope": -2.0, "static_limit_ratio_to_zero": True},
        "notes": "Static-limit test of Φ=ΦP+ΦL with optional absorber. For γ=0 in a periodic box, ΦL generally persists as standing waves."
    }))

    def compute_phi_tt(phi_now, pi_now):
        lap = laplacian_fft(xp, be.fftn, be.ifftn, field=phi_now, k2=k2)
        # φ_tt + 2γ φ_t = S^2 (∇^2 φ - 4πG ρ_src)
        # => φ_tt = S^2 (...) - 2γ φ_t
        return (S * S) * (lap - four_pi_G * rho_src) - (2.0 * gamma) * pi_now

    # One initial diagnostics for slopes from Poisson
    def diag(tag: str, step: int, t: float, phi_now):
        # Static-limit lag measure
        lag = (phi_now - phiP).astype(xp.float64)
        lag_rms = xp.sqrt(xp.mean(lag * lag))
        phiP_rms_local = xp.sqrt(xp.mean(phiP.astype(xp.float64) ** 2))
        ratio = float(be.asnumpy(lag_rms / xp.maximum(phiP_rms_local, 1e-30)))

        # Spherical averages (do on backend, fit on CPU)
        phi_bar = spherical_average(xp, phi_now, bin_idx_flat, counts, nbins)
        phi_bar_np = be.asnumpy(phi_bar).astype(np.float64)

        # Reference subtract (constant offset irrelevant)
        outer_n = max(8, int(0.1 * nbins))
        phi_ref = float(np.mean(phi_bar_np[-outer_n:]))
        dphi = phi_bar_np - phi_ref

        g_r = -np.gradient(phi_bar_np, r_centers_cpu)

        # Fit thresholds
        mask_fit = (r_centers_cpu >= fit_rmin) & (r_centers_cpu <= fit_rmax)
        dp_peak = float(np.nanmax(np.abs(dphi[mask_fit]))) if np.any(mask_fit) else 0.0
        g_peak = float(np.nanmax(np.abs(g_r[mask_fit]))) if np.any(mask_fit) else 0.0
        ymin_dp = 0.0 if dp_peak < 1e-20 else max(1e-20, dp_peak * 1e-4)
        ymin_g = 0.0 if g_peak < 1e-24 else max(1e-24, g_peak * 1e-4)

        p_phi, _, n_phi = fit_powerlaw(r_centers_cpu, dphi, fit_rmin, fit_rmax, y_min_abs=ymin_dp)
        p_g, _, n_g = fit_powerlaw(r_centers_cpu, g_r, fit_rmin, fit_rmax, y_min_abs=ymin_g)

        out = {
            "event": "diag",
            "tag": tag,
            "step": int(step),
            "t": float(t),
            "lag_rms_over_phiP_rms": ratio,
            "fits": {
                "dPhi_slope": p_phi,
                "dPhi_npts": n_phi,
                "g_slope": p_g,
                "g_npts": n_g,
            },
            "refs": {"phi_ref": phi_ref},
            "fit_cutoffs": {"ymin_dPhi": ymin_dp, "ymin_g": ymin_g},
            "wall_s": float(time.time() - wall0),
        }
        print(json.dumps(out))

    diag(tag="init_phi", step=0, t=0.0, phi_now=phi)

    # Time loop (velocity-Verlet)
    t_sim = 0.0
    for step in range(1, steps + 1):
        phi_tt = compute_phi_tt(phi, pi)
        pi_half = pi + (0.5 * dt) * phi_tt
        phi_new = phi + dt * pi_half
        phi_tt_new = compute_phi_tt(phi_new, pi_half)
        pi_new = pi_half + (0.5 * dt) * phi_tt_new

        phi, pi = phi_new, pi_new
        t_sim += dt

        if (step % diag_every) == 0 or step == steps:
            diag(tag="evolve", step=step, t=t_sim, phi_now=phi)

    print(json.dumps({
        "event": "done",
        "t_final": float(t_sim),
        "wall_s": float(time.time() - wall0),
        "backend": be.name,
        "note": "If γ>0, lag_rms_over_phiP_rms should decay; slopes should approach (-1,-2) in the chosen fit window."
    }))


if __name__ == "__main__":
    main()
