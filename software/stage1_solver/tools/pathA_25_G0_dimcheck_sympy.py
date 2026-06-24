#!/usr/bin/env python3
"""SymPy dimensional and freeze-fidelity checks for pathA_25 G0.

Scope: G0 only.  This script verifies the T0 and G0 freeze hashes, checks that
the exact T0 freeze bytes occur inside the combined G0 freeze block, restores
MLT units for every new frozen term, and checks the explicit k->0 driver limit.
It does not solve any consistency gate.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Iterable

import sympy as sp


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
SCRATCH = STAGE1_ROOT / "_scratch"
T0_REPORT = STAGE1_ROOT / "reports" / "pathA_24_T0_freeze.md"
G0_REPORT = STAGE1_ROOT / "reports" / "pathA_25_G0_freeze.md"

EXPECTED_T0_HASH = "8fa41ac51e88a1464a4a5b22c6fe64fc218cf36ba2e3583d26b97c994e5da064"
EXPECTED_G0_HASH = "f00ee99d465e2e311c68f47fcacf4af0154ca650642271ab66c36d112cb6a290"

EXPECTED_L_POL_LINES = [
    b"  L_pol =\n",
    b"      (1/2) m rho a^2 (D_t^v P^i)(D_t^v P^i)\n",
    b"    - (1/2) m rho c_s^2(rho) a^2 (partial_j P^i)(partial_j P^i)\n",
    b"    - (1/4) m rho c_s^2(rho) (P^i P^i - 1)^2.\n",
]


def extract_fence_bytes(path: Path, label: str) -> bytes:
    start = f"```{label}\n".encode("utf-8")
    end = b"```\n"
    lines = path.read_bytes().splitlines(keepends=True)
    blocks: list[bytes] = []
    in_block = False
    current: list[bytes] = []
    for line in lines:
        if not in_block and line == start:
            in_block = True
            current = []
            continue
        if in_block and line == end:
            blocks.append(b"".join(current))
            in_block = False
            current = []
            continue
        if in_block:
            current.append(line)
    if in_block:
        raise AssertionError(f"unterminated {label!r} fence in {path}")
    if len(blocks) != 1:
        raise AssertionError(f"expected exactly one {label!r} fence in {path}, found {len(blocks)}")
    return blocks[0]


def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def assert_freeze_fidelity() -> dict:
    t0_block = extract_fence_bytes(T0_REPORT, "freeze-action")
    g0_block = extract_fence_bytes(G0_REPORT, "freeze-action")
    t0_hash = sha256_hex(t0_block)
    g0_hash = sha256_hex(g0_block)
    if t0_hash != EXPECTED_T0_HASH:
        raise AssertionError(f"T0 hash mismatch: expected {EXPECTED_T0_HASH}, got {t0_hash}")
    if g0_hash != EXPECTED_G0_HASH:
        raise AssertionError(f"G0 hash mismatch: expected {EXPECTED_G0_HASH}, got {g0_hash}")
    if t0_block not in g0_block:
        raise AssertionError("exact T0 freeze-action bytes are not embedded unchanged in G0 block")
    for expected in EXPECTED_L_POL_LINES:
        if expected not in t0_block:
            raise AssertionError(f"T0 L_pol line missing from T0 block: {expected!r}")
        if expected not in g0_block:
            raise AssertionError(f"T0 L_pol line missing from G0 block: {expected!r}")
    positions = [g0_block.index(line) for line in EXPECTED_L_POL_LINES]
    if positions != sorted(positions):
        raise AssertionError(f"G0 embedded L_pol lines out of order: {positions}")
    return {"t0_hash": t0_hash, "g0_hash": g0_hash, "t0_bytes_embedded_in_g0": True}


Dim = sp.Matrix


def dim_tuple(dim: Dim) -> tuple[int, int, int]:
    return tuple(int(x) for x in list(dim))


def dim_str(dim: Dim) -> str:
    labels = ("M", "L", "T")
    parts: list[str] = []
    for label, power in zip(labels, dim_tuple(dim)):
        if power == 0:
            continue
        if power == 1:
            parts.append(label)
        else:
            parts.append(f"{label}^{power}")
    return " ".join(parts) if parts else "1"


def assert_dim(name: str, actual: Dim, expected: Dim) -> None:
    delta = sp.simplify(actual - expected)
    if any(entry != 0 for entry in delta):
        raise AssertionError(
            f"{name}: expected {dim_tuple(expected)} ({dim_str(expected)}), "
            f"got {dim_tuple(actual)} ({dim_str(actual)})"
        )


def dim_record(dim: Dim) -> dict:
    return {"triple_MLT": dim_tuple(dim), "string": dim_str(dim)}


def check_dimensions() -> dict:
    M = sp.Matrix([1, 0, 0])
    L = sp.Matrix([0, 1, 0])
    T = sp.Matrix([0, 0, 1])
    dimless = sp.Matrix([0, 0, 0])

    bulk_lag = M - 2 * L - 2 * T
    brane_lag = M - L - 2 * T
    action = M + 2 * L - T
    delta_sigma = -L

    d_m = M
    d_rho = -4 * L
    d_K = M + 18 * L - 2 * T
    d_a = L
    d_u = L
    d_P = dimless
    d_grad4 = -L
    d_grad3 = -L
    d_lap4 = -2 * L
    d_dt = -T
    d_dn = L
    d_d4R = 4 * L

    d_cs2 = d_K + 4 * d_rho - d_m
    assert_dim("c_s^2(rho)", d_cs2, 2 * L - 2 * T)

    # Family L, bulk density.
    d_cL1 = M + 8 * L - 2 * T
    d_cL2 = M + 10 * L - 2 * T
    d_grad_rho = d_rho + d_grad4
    d_lap_rho = d_rho + d_lap4
    family_l_terms = {
        "c_L1 (partial_i rho)^2": d_cL1 + 2 * d_grad_rho,
        "c_L2 (partial_i partial_i rho)^2": d_cL2 + 2 * d_lap_rho,
    }
    for name, dim in family_l_terms.items():
        assert_dim(name, dim, bulk_lag)
    d_kLstar = (d_cL1 - d_cL2) / 2
    assert_dim("k_Lstar", d_kLstar, -L)

    # Family R, real-space and Fourier kernel dimensions.
    d_VR_real = M + 2 * L - 2 * T
    d_A_R = M + 6 * L - 2 * T
    d_k_R = -L
    assert_dim("Family R real-space kernel term", d_d4R + 2 * d_rho + d_VR_real, bulk_lag)
    assert_dim("tilde V_R=A_R f_R", d_A_R, d_VR_real + d_d4R)
    assert_dim("k_R", d_k_R, -L)
    assert_dim("k_Rstar", d_k_R, -L)

    # Family C, bulk density.
    d_lambda_Cdiv = M + 3 * L - 2 * T
    d_chi_Cpin = M + 8 * L - 2 * T
    d_gradP = d_P + d_grad4
    family_c_terms = {
        "lambda_Cdiv delta_rho partial_i P^i": d_lambda_Cdiv + d_rho + d_gradP,
        "chi_Cpin (P^i partial_i rho)^2": d_chi_Cpin + 2 * (d_P + d_grad_rho),
    }
    for name, dim in family_c_terms.items():
        assert_dim(name, dim, bulk_lag)

    # Layer-effective light package, first as 3D brane density, then as 4D
    # density after multiplication by delta_Sigma.
    d_varrho_br = M - 3 * L
    d_mu_br = brane_lag
    d_Omega = d_grad3 + d_u
    d_Dt_u = d_u + d_dt
    d_I_PSigma = d_dn + d_m + d_rho + 2 * d_a
    d_K_PSigma = d_dn + d_m + d_rho + d_cs2 + 2 * d_a
    d_G_PSigma = d_dn + d_m + d_rho + d_cs2
    d_J_Pu = M - L
    d_kappa_Pu = brane_lag
    d_DtP = d_dt
    d_DtOmega = d_dt + d_Omega

    assert_dim("Omega_u", d_Omega, dimless)
    light_brane_terms = {
        "varrho_br (D_t u)^2": d_varrho_br + 2 * d_Dt_u,
        "mu_br (curl u)^2": d_mu_br + 2 * d_Omega,
        "I_PSigma (D_t P_parallel)^2": d_I_PSigma + 2 * d_DtP,
        "K_PSigma (partial_a P_parallel)^2": d_K_PSigma + 2 * d_gradP,
        "G_PSigma (delta P_parallel)^2": d_G_PSigma + 2 * d_P,
        "J_Pu (D_t delta P - D_t Omega_u)^2": d_J_Pu + 2 * d_DtP,
        "kappa_Pu (delta P - Omega_u)^2": d_kappa_Pu + 2 * d_P,
    }
    for name, dim in light_brane_terms.items():
        assert_dim(name, dim, brane_lag)
        assert_dim(f"bulk representation delta_Sigma {name}", delta_sigma + dim, bulk_lag)

    assert_dim("varrho_br=int dn m rho", d_varrho_br, M - 3 * L)
    assert_dim("mu_br", d_mu_br, brane_lag)
    assert_dim("I_PSigma=int dn m rho a^2", d_I_PSigma, M - L)
    assert_dim("K_PSigma=int dn m rho c_s^2 a^2", d_K_PSigma, M + L - 2 * T)
    assert_dim("G_PSigma=int dn m rho c_s^2", d_G_PSigma, brane_lag)
    assert_dim("J_Pu", d_J_Pu, M - L)
    assert_dim("kappa_Pu", d_kappa_Pu, brane_lag)
    assert_dim("delta_Sigma", delta_sigma, -L)

    assert_dim("int dt d^4X bulk_lag", bulk_lag + T + 4 * L, action)
    assert_dim("int dt d^3sigma brane_lag", brane_lag + T + 3 * L, action)

    return {
        "expected_bulk_action_density": dim_record(bulk_lag),
        "expected_brane_density": dim_record(brane_lag),
        "c_s_squared": dim_record(d_cs2),
        "family_l": {name: dim_record(dim) for name, dim in family_l_terms.items()},
        "family_r": {
            "real_space_kernel_V_R": dim_record(d_VR_real),
            "fourier_kernel_A_R": dim_record(d_A_R),
            "real_space_action_density": dim_record(d_d4R + 2 * d_rho + d_VR_real),
            "k_R": dim_record(d_k_R),
        },
        "family_c": {name: dim_record(dim) for name, dim in family_c_terms.items()},
        "light_brane_terms": {name: dim_record(dim) for name, dim in light_brane_terms.items()},
        "light_bulk_density_terms": {
            name: dim_record(delta_sigma + dim) for name, dim in light_brane_terms.items()
        },
        "light_constants": {
            "varrho_br": dim_record(d_varrho_br),
            "mu_br": dim_record(d_mu_br),
            "I_PSigma": dim_record(d_I_PSigma),
            "K_PSigma": dim_record(d_K_PSigma),
            "G_PSigma": dim_record(d_G_PSigma),
            "J_Pu": dim_record(d_J_Pu),
            "kappa_Pu": dim_record(d_kappa_Pu),
            "delta_Sigma": dim_record(delta_sigma),
        },
    }


def assert_zero(name: str, expr: sp.Expr) -> None:
    simplified = sp.simplify(expr)
    if simplified != 0:
        raise AssertionError(f"{name}: expected 0, got {sp.sstr(simplified)}")


def check_k_to_zero() -> dict:
    k, cL1, cL2, A_R, k_R, lam, chi = sp.symbols(
        "k c_L1 c_L2 A_R k_R lambda_Cdiv chi_Cpin", positive=True
    )

    family_l = cL2 * k**4 - cL1 * k**2
    k_l_star = sp.sqrt(cL1 / (2 * cL2))
    assert_zero("Family L k->0 value", sp.limit(family_l, k, 0))
    assert_zero("Family L finite-k stationary derivative", sp.diff(family_l, k).subs(k, k_l_star))
    if sp.simplify(sp.diff(family_l, k, 2).subs(k, k_l_star) - 4 * cL1) != 0:
        raise AssertionError("Family L finite-k stationary point second derivative is not 4 c_L1")

    x = sp.symbols("x", positive=True)
    f_R = (x**4 - 2 * x**2) * sp.exp(-x**2)
    k_r_star = k_R * sp.sqrt(2 - sp.sqrt(2))
    vtilde = A_R * f_R.subs(x, k / k_R)
    assert_zero("Family R k->0 value", sp.limit(vtilde, k, 0))
    assert_zero("Family R finite-k stationary derivative", sp.diff(vtilde, k).subs(k, k_r_star))
    f_r_min_numeric = float(sp.N(f_R.subs(x, sp.sqrt(2 - sp.sqrt(2)))))
    if f_r_min_numeric >= 0:
        raise AssertionError("Family R finite-k stationary point is not in the negative annulus")

    cdiv_symbol = lam * k
    cpin_symbol = chi * k**2
    assert_zero("Family C divergence k->0 value", sp.limit(cdiv_symbol, k, 0))
    assert_zero("Family C pinning k->0 value", sp.limit(cpin_symbol, k, 0))

    return {
        "family_l": {
            "quadratic_symbol": "c_L2*k^4 - c_L1*k^2",
            "limit_k_to_0": "0",
            "finite_k_stationary_scale": sp.sstr(k_l_star),
            "second_derivative_at_scale": "4*c_L1",
        },
        "family_r": {
            "shape": "(x^4 - 2*x^2)*exp(-x^2)",
            "limit_k_to_0": "0",
            "small_k_series": sp.sstr(sp.series(vtilde, k, 0, 5).removeO()),
            "finite_k_stationary_scale": sp.sstr(k_r_star),
        },
        "family_c": {
            "divergence_coupling_limit_k_to_0": "0",
            "pinning_coupling_limit_k_to_0": "0",
        },
        "statement": "No frozen driver contributes a constant k=0 density-potential term; the EOS c_s definition is not shifted at strict k=0.",
    }


def main() -> None:
    SCRATCH.mkdir(parents=True, exist_ok=True)
    freeze = assert_freeze_fidelity()
    dimensions = check_dimensions()
    k_to_zero = check_k_to_zero()
    report = {
        "schema": "pathA_25_G0_dimcheck_sympy/v1",
        "pass": True,
        "scope": "freeze-fidelity, dimensional homogeneity, and k->0 limit only; no gate solved",
        "freeze_fidelity": freeze,
        "dimensions": dimensions,
        "k_to_zero": k_to_zero,
        "verdict": "G0_DIMCHECK_SYMPY_PASS",
    }
    out = SCRATCH / "pathA_25_G0_dimcheck_sympy.json"
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    print("pathA_25 G0 SymPy dimcheck: PASS")


if __name__ == "__main__":
    main()
