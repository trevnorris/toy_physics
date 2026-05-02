#!/usr/bin/env python3
"""Convert solver-side state dumps into the runtime-monitor snapshot schemas."""
from __future__ import annotations

import argparse
import pathlib
from typing import Any

import numpy as np

from cfd_runtime_monitor_postprocess import first_derivative, uniform_spacing


FORWARD_OPTIONAL_KEYS = ("dt_rho", "phi3", "rho0", "N_probe", "Phi_eff")


def load_npz(path: pathlib.Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as raw:
        return {key: raw[key] for key in raw.files}


def save_npz(path: pathlib.Path, payload: dict[str, Any]) -> None:
    serializable = {key: np.asarray(value) for key, value in payload.items()}
    np.savez(path, **serializable)


def read_scalar(data: dict[str, np.ndarray], *keys: str, default: float | None = None) -> float:
    for key in keys:
        if key in data:
            return float(np.asarray(data[key]).reshape(()))
    if default is None:
        raise KeyError(f"missing scalar keys {keys}")
    return float(default)


def extract_complex_wavefunction(data: dict[str, np.ndarray]) -> np.ndarray:
    if "psi" in data:
        psi = np.asarray(data["psi"])
        if np.iscomplexobj(psi):
            return psi.astype(np.complex128, copy=False)
    if "psi_real" in data and "psi_imag" in data:
        return np.asarray(data["psi_real"], dtype=float) + 1j * np.asarray(data["psi_imag"], dtype=float)
    raise KeyError("wavefunction snapshot must contain either complex psi or psi_real/psi_imag")


def make_uniform_weight(w: np.ndarray) -> np.ndarray:
    span = float(w[-1] - w[0])
    if span <= 0.0:
        raise ValueError("w axis must span a positive interval to define a uniform projection weight")
    return np.full_like(w, 1.0 / span, dtype=float)


def adapt_wavefunction_4d(
    data: dict[str, np.ndarray],
    *,
    mass: float,
    hbar: float,
    periodic_xyz: bool = True,
    periodic_w: bool = False,
    allow_uniform_weight: bool = False,
) -> dict[str, Any]:
    for key in ("x", "y", "z", "w"):
        if key not in data:
            raise KeyError(f"wavefunction snapshot missing axis {key}")

    x = np.asarray(data["x"], dtype=float)
    y = np.asarray(data["y"], dtype=float)
    z = np.asarray(data["z"], dtype=float)
    w = np.asarray(data["w"], dtype=float)
    dx = uniform_spacing(x)
    dy = uniform_spacing(y)
    dz = uniform_spacing(z)
    dw = uniform_spacing(w)

    psi = extract_complex_wavefunction(data)
    if psi.ndim != 4:
        raise ValueError("wavefunction psi must have shape (nx, ny, nz, nw)")
    if psi.shape != (x.size, y.size, z.size, w.size):
        raise ValueError("wavefunction shape does not match x/y/z/w axis lengths")

    dpsi_dx = first_derivative(psi, dx, 0, periodic_xyz)
    dpsi_dy = first_derivative(psi, dy, 1, periodic_xyz)
    dpsi_dz = first_derivative(psi, dz, 2, periodic_xyz)
    dpsi_dw = first_derivative(psi, dw, 3, periodic_w)

    pref = float(hbar) / float(mass)
    rho = np.abs(psi) ** 2
    psi_conj = np.conjugate(psi)
    jx = pref * np.imag(psi_conj * dpsi_dx)
    jy = pref * np.imag(psi_conj * dpsi_dy)
    jz = pref * np.imag(psi_conj * dpsi_dz)
    jw = pref * np.imag(psi_conj * dpsi_dw)

    if "W" in data:
        W = np.asarray(data["W"], dtype=float)
    else:
        if not allow_uniform_weight:
            raise KeyError("wavefunction snapshot missing W; pass --allow-uniform-W to use a flat projection weight")
        W = make_uniform_weight(w)
    if W.shape != w.shape:
        raise ValueError("W must have the same shape as w")

    payload: dict[str, Any] = {
        "rho": rho,
        "jx": jx,
        "jy": jy,
        "jz": jz,
        "jw": jw,
        "W": W,
        "x": x,
        "y": y,
        "z": z,
        "w": w,
    }
    if "dWdw" in data:
        payload["dWdw"] = np.asarray(data["dWdw"], dtype=float)
    for key in FORWARD_OPTIONAL_KEYS:
        if key in data:
            payload[key] = data[key]
    return payload


def compute_monopole_source(data: dict[str, np.ndarray], x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    if "S_rho" in data:
        return np.asarray(data["S_rho"], dtype=float)
    if "W" not in data:
        raise KeyError("monopole snapshot needs either S_rho directly or W plus Mdot/lambda to reconstruct it")
    W = np.asarray(data["W"], dtype=float)
    Mdot = read_scalar(data, "Mdot")
    lam = read_scalar(data, "lambda_bulk", "lambda")
    dx = uniform_spacing(x)
    dy = uniform_spacing(y)
    dz = uniform_spacing(z)
    V_domain = read_scalar(data, "V_domain", default=float(x.size * y.size * z.size) * dx * dy * dz)
    return -Mdot * W + (lam * Mdot / V_domain) * np.ones_like(W, dtype=float)


def adapt_monopole_3d(data: dict[str, np.ndarray]) -> dict[str, Any]:
    for key in ("rho", "mx", "my", "mz", "x", "y", "z"):
        if key not in data:
            raise KeyError(f"monopole snapshot missing required key {key}")
    rho = np.asarray(data["rho"], dtype=float)
    mx = np.asarray(data["mx"], dtype=float)
    my = np.asarray(data["my"], dtype=float)
    mz = np.asarray(data["mz"], dtype=float)
    x = np.asarray(data["x"], dtype=float)
    y = np.asarray(data["y"], dtype=float)
    z = np.asarray(data["z"], dtype=float)
    if rho.shape != mx.shape or rho.shape != my.shape or rho.shape != mz.shape:
        raise ValueError("rho/mx/my/mz must share the same 3D shape")
    if rho.ndim != 3 or rho.shape != (x.size, y.size, z.size):
        raise ValueError("rho/mx/my/mz must have shape (nx, ny, nz) matching x/y/z")

    payload: dict[str, Any] = {
        "rho_brane": rho,
        "Jx_brane": mx,
        "Jy_brane": my,
        "Jz_brane": mz,
        "S_rho": compute_monopole_source(data, x, y, z),
        "x": x,
        "y": y,
        "z": z,
    }
    for key in FORWARD_OPTIONAL_KEYS:
        if key in data:
            payload[key] = data[key]
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    wave = sub.add_parser("wavefunction-4d", help="adapt a 4D complex wavefunction snapshot")
    wave.add_argument("input")
    wave.add_argument("output")
    wave.add_argument("--mass", type=float, default=None, help="particle mass for current extraction")
    wave.add_argument("--hbar", type=float, default=None, help="hbar for current extraction")
    wave.add_argument("--nonperiodic-xyz", action="store_true")
    wave.add_argument("--periodic-w", action="store_true")
    wave.add_argument("--allow-uniform-W", action="store_true", help="use a flat projection weight when W is absent")

    mono = sub.add_parser("monopole-3d", help="adapt a 3D single-throat monopole state dump")
    mono.add_argument("input")
    mono.add_argument("output")

    args = parser.parse_args()
    data = load_npz(pathlib.Path(args.input))

    if args.command == "wavefunction-4d":
        mass = args.mass if args.mass is not None else read_scalar(data, "mass", "m", default=1.0)
        hbar = args.hbar if args.hbar is not None else read_scalar(data, "hbar", default=1.0)
        payload = adapt_wavefunction_4d(
            data,
            mass=mass,
            hbar=hbar,
            periodic_xyz=not args.nonperiodic_xyz,
            periodic_w=args.periodic_w,
            allow_uniform_weight=args.allow_uniform_W,
        )
    else:
        payload = adapt_monopole_3d(data)

    save_npz(pathlib.Path(args.output), payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
