#!/usr/bin/env python3
"""stepC_throat_potential_flow.py

Step C: Bridge test — does a throat/sink produce a Poisson monopole (1/r potential,
1/r^2 field) in the *potential-flow* (quasi-incompressible, irrotational) limit?

This script is intentionally *not* a hydrodynamics integrator.
It implements the analytic static-limit reduction:

    ∂tρ + ∇·(ρ v) = S_ρ
    (ρ≈ρ0 const, ∂tρ≈0)  ⇒  ρ0 ∇·v = S_ρ
    (irrotational, v=∇ψ) ⇒  ∇²ψ = S_ρ/ρ0

For a localized sink, S_ρ ≈ -Mdot δ(x), the solution is

    ψ(r) = - Mdot / (4π ρ0 r),
    v_r(r) = ∂r ψ = + Mdot / (4π ρ0 r²)

So the *potential* has slope -1 and the *field/velocity* has slope -2.

We use a finite-width Gaussian kernel W(x) normalized to ∫W dV = 1.
On a periodic grid, Poisson solvability requires mean(source)=0, so we
automatically subtract the mean (equivalently add a uniform compensator):

    q(x) = (S_ρ/ρ0) - <S_ρ/ρ0>

This is the same conceptual "subtract the cosmic mean" step used in
`stepB_scalar_poisson_lag.py` (rho_src = rho - <rho>). The uniform mode
represents a global background ("dark energy" / cosmological constant sector)
and is not part of the Newtonian 1/r monopole.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass
from typing import Tuple

import numpy as np


def get_xp(backend: str):
    backend = backend.lower().strip()
    if backend == "numpy":
        return np
    if backend == "cupy":
        import cupy as cp  # type: ignore
        return cp
    raise ValueError(f"Unknown backend: {backend!r} (expected 'numpy' or 'cupy').")


def asnumpy(x, xp):
    return x.get() if xp.__name__ == "cupy" else x


def make_gaussian_kernel(xp, N: int, L: float, sigma: float, dtype) -> "xp.ndarray":
    """Return W(x) on an N^3 periodic grid with ∫W dV = 1."""
    dx = L / N
    # Coordinates centered at 0 with periodic wrap.
    # Use integer grid indices to avoid floating drift.
    i = xp.arange(N, dtype=dtype)
    # map to [-N/2, N/2)
    ii = i - (N // 2)
    x = ii * dx
    X, Y, Z = xp.meshgrid(x, x, x, indexing="ij")
    r2 = X * X + Y * Y + Z * Z
    W = xp.exp(-0.5 * r2 / (sigma * sigma)).astype(dtype)
    # Normalize to integral 1
    dV = dx ** 3
    norm = xp.sum(W) * dV
    W = W / norm
    return W


def kgrid(xp, N: int, L: float, dtype) -> Tuple["xp.ndarray", "xp.ndarray", "xp.ndarray", "xp.ndarray"]:
    """Return Fourier wavenumber grids (kx,ky,kz,k2) for an N^3 periodic grid."""
    # xp.fft.fftfreq returns cycles per unit; multiply by 2π for radians.
    k1 = (2.0 * math.pi) * xp.fft.fftfreq(N, d=L / N).astype(dtype)
    kx, ky, kz = xp.meshgrid(k1, k1, k1, indexing="ij")
    k2 = kx * kx + ky * ky + kz * kz
    return kx, ky, kz, k2


def poisson_solve_periodic(xp, rhs: "xp.ndarray", k2: "xp.ndarray") -> "xp.ndarray":
    """Solve ∇²φ = rhs on a periodic grid with <rhs>=0, setting k=0 mode to 0."""
    rhs_k = xp.fft.fftn(rhs)
    phi_k = xp.zeros_like(rhs_k)
    mask = k2 > 0
    phi_k[mask] = -rhs_k[mask] / k2[mask]
    # k=0 mode left at 0
    phi = xp.fft.ifftn(phi_k).real
    return phi


def spectral_grad(xp, phi: "xp.ndarray", kx, ky, kz) -> Tuple["xp.ndarray", "xp.ndarray", "xp.ndarray"]:
    """Compute ∇phi using spectral differentiation on a periodic grid."""
    phi_k = xp.fft.fftn(phi)
    dphix = xp.fft.ifftn(1j * kx * phi_k).real
    dphiy = xp.fft.ifftn(1j * ky * phi_k).real
    dphiz = xp.fft.ifftn(1j * kz * phi_k).real
    return dphix, dphiy, dphiz


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
    # bin index 0..nbins-1
    b = np.digitize(r.ravel(), edges) - 1
    b = np.clip(b, 0, nbins - 1)
    counts = np.bincount(b, minlength=nbins).astype(np.int64)
    # centers for plotting/fits
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, b.astype(np.int32), counts


def spherical_average(xp, field: "xp.ndarray", bin_idx_flat_np: np.ndarray, counts_np: np.ndarray, nbins: int) -> np.ndarray:
    """Spherical average using precomputed numpy bin indices (works for cupy/numpy)."""
    # Flatten on device, then scatter-add into bins on device if cupy
    flat = field.ravel()
    if xp.__name__ == "cupy":
        import cupyx  # type: ignore
        b = xp.asarray(bin_idx_flat_np)
        sums = cupyx.scatter_add(xp.zeros(nbins, dtype=flat.dtype), b, flat)
        avg = sums / xp.asarray(counts_np).astype(flat.dtype)
        return avg.get()
    else:
        sums = np.bincount(bin_idx_flat_np, weights=flat, minlength=nbins)
        avg = sums / np.maximum(counts_np, 1)
        return avg


def fit_powerlaw(r: np.ndarray, y: np.ndarray, rmin: float, rmax: float, y_min_abs: float = 0.0) -> Tuple[float, float, int]:
    """Fit |y| ~ A r^p over [rmin,rmax] in log-log."""
    m = (r >= rmin) & (r <= rmax) & np.isfinite(y)
    yabs = np.abs(y[m])
    rr = r[m]
    if y_min_abs > 0:
        keep = yabs > y_min_abs
        yabs = yabs[keep]
        rr = rr[keep]
    if rr.size < 3:
        return float("nan"), float("nan"), int(rr.size)
    lx = np.log(rr)
    ly = np.log(yabs)
    p, b = np.polyfit(lx, ly, 1)
    return float(p), float(b), int(rr.size)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", choices=["numpy", "cupy"], default="cupy")
    ap.add_argument("--dtype", choices=["float32", "float64"], default="float64")
    ap.add_argument("--N", type=int, default=256)
    ap.add_argument("--L", type=float, default=200.0)
    ap.add_argument("--Mdot", type=float, default=20.0)
    ap.add_argument("--rho0", type=float, default=1.0)
    ap.add_argument("--sigma", type=float, default=2.0, help="Gaussian width of throat kernel.")
    ap.add_argument("--nbins", type=int, default=180)
    ap.add_argument("--fit_rmin", type=float, default=8.0)
    ap.add_argument("--fit_rmax", type=float, default=50.0)
    ap.add_argument("--sample_n", type=int, default=40)
    args = ap.parse_args()

    xp = get_xp(args.backend)
    dtype = np.float32 if args.dtype == "float32" else np.float64
    if xp.__name__ == "cupy":
        import cupy as cp  # type: ignore
        dtype = cp.float32 if args.dtype == "float32" else cp.float64

    N = int(args.N)
    L = float(args.L)
    dx = L / N
    Mdot = float(args.Mdot)
    rho0 = float(args.rho0)
    sigma = float(args.sigma)
    nbins = int(args.nbins)

    wall0 = time.time()

    # Precompute radial binning on CPU once.
    r_centers, bin_idx_flat_np, counts_np = radial_bins_numpy(N, L, nbins)

    # Build sink kernel and source term q = S_rho/rho0 with mean removed.
    W = make_gaussian_kernel(xp, N=N, L=L, sigma=sigma, dtype=dtype)
    S_rho = (-Mdot) * W  # ∫ S_rho dV = -Mdot
    q = S_rho / rho0
    q = q - xp.mean(q)  # enforce solvability (subtract mean)

    # Solve Poisson ∇²ψ = q
    kx, ky, kz, k2 = kgrid(xp, N, L, dtype)
    psi = poisson_solve_periodic(xp, q, k2)

    # Field/velocity from potential: v = ∇ψ
    dpsix, dpsiy, dpsiz = spectral_grad(xp, psi, kx, ky, kz)
    vmag = xp.sqrt(dpsix * dpsix + dpsiy * dpsiy + dpsiz * dpsiz)

    # Spherical averages (CPU arrays)
    psi_bar = spherical_average(xp, psi, bin_idx_flat_np, counts_np, nbins)
    v_bar = spherical_average(xp, vmag, bin_idx_flat_np, counts_np, nbins)

    # Reference-subtract psi using outer shells
    outer_n = max(5, int(0.10 * nbins))
    psi_ref = float(np.mean(psi_bar[-outer_n:]))
    dpsi = psi_bar - psi_ref

    # Effective field magnitude from radial derivative of spherical mean.
    dpsidr = np.gradient(psi_bar, r_centers)
    g_r = np.abs(dpsidr)

    # Fit power laws
    # Use cutoffs relative to peak in fit window.
    mask_fit = (r_centers >= args.fit_rmin) & (r_centers <= args.fit_rmax)
    dpsi_fit_abs = np.abs(dpsi[mask_fit])
    g_fit_abs = np.abs(g_r[mask_fit])
    dpsi_peak = float(np.nanmax(dpsi_fit_abs)) if dpsi_fit_abs.size else 0.0
    g_peak = float(np.nanmax(g_fit_abs)) if g_fit_abs.size else 0.0
    ymin_dpsi = 0.0 if dpsi_peak < 1e-30 else max(1e-30, dpsi_peak * 1e-3)
    ymin_g = 0.0 if g_peak < 1e-30 else max(1e-30, g_peak * 1e-3)
    ppsi, bpsi, npsi = fit_powerlaw(r_centers, dpsi, args.fit_rmin, args.fit_rmax, y_min_abs=ymin_dpsi)
    pg, bg, ng = fit_powerlaw(r_centers, g_r, args.fit_rmin, args.fit_rmax, y_min_abs=ymin_g)

    # Also fit ψ̄(r) ≈ A/r + C over the same fit window.
    # In a periodic, mean-subtracted Poisson solve, ψ is defined up to an additive constant,
    # and small residual offsets can bias a pure power-law fit of (ψ̄-ψ_ref). The A/r + C fit
    # is the more robust "monopole" diagnostic.
    A_fit = float("nan")
    C_fit = float("nan")
    A_expected = float(Mdot / (4.0 * math.pi * rho0))
    rel_err_A = float("nan")
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

    # "Cleaned" Δψ that removes the best-fit constant C before slope fitting
    dpsi_clean = psi_bar - C_fit
    dpsi_clean_fit_abs = np.abs(dpsi_clean[mask_fit])
    dpsi_clean_peak = float(np.nanmax(dpsi_clean_fit_abs)) if dpsi_clean_fit_abs.size else 0.0
    ymin_dpsi_clean = 0.0 if dpsi_clean_peak < 1e-30 else max(1e-30, dpsi_clean_peak * 1e-3)
    ppsi_clean, bpsi_clean, npsi_clean = fit_powerlaw(
        r_centers, dpsi_clean, args.fit_rmin, args.fit_rmax, y_min_abs=ymin_dpsi_clean
    )

    # Sample points
    idx = np.where(mask_fit)[0]
    stride = max(1, len(idx) // max(1, int(args.sample_n)))
    idx_s = idx[::stride]
    sample = [
        {
            "r": float(r_centers[i]),
            "dpsi": float(dpsi[i]),
            "g": float(g_r[i]),
            "psi_bar": float(psi_bar[i]),
            "v_bar": float(v_bar[i]),
        }
        for i in idx_s
    ]

    out_init = {
        "event": "init",
        "backend": args.backend,
        "N": N,
        "L": L,
        "dx": dx,
        "dtype": args.dtype,
        "throat": {"Mdot": Mdot, "rho0": rho0, "sigma": sigma},
        "poisson": {"equation": "∇²ψ = (S_ρ/ρ0) - <S_ρ/ρ0>", "k0_mode": "set to 0"},
        "diag": {"nbins": nbins, "fit_rmin": args.fit_rmin, "fit_rmax": args.fit_rmax},
        "targets": {"dpsi_slope": -1.0, "g_slope": -2.0},
    }
    print(json.dumps(out_init))

    out = {
        "event": "diag",
        "fits": {"dpsi_slope": ppsi, "dpsi_npts": npsi, "dpsi_clean_slope": ppsi_clean, "dpsi_clean_npts": npsi_clean, "g_slope": pg, "g_npts": ng},
        "refs": {"psi_ref": psi_ref, "C_fit": C_fit},
        "mono_fit": {"A_fit": A_fit, "A_expected": A_expected, "rel_err_A": rel_err_A, "Mdot_fit": float(4.0 * math.pi * rho0 * A_fit)},
        "fit_cutoffs": {"ymin_dpsi": ymin_dpsi, "ymin_dpsi_clean": ymin_dpsi_clean, "ymin_g": ymin_g, "dpsi_peak_fit": dpsi_peak, "dpsi_clean_peak_fit": dpsi_clean_peak, "g_peak_fit": g_peak},
        "sample_points": sample,
        "wall_s": float(time.time() - wall0),
    }
    print(json.dumps(out))


if __name__ == "__main__":
    main()
