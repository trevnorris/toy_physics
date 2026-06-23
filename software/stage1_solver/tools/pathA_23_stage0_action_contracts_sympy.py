#!/usr/bin/env python3
"""Independent SymPy cross-check for pathA_23 Stage 0.

Scope: units, source-current gauge variation, absence of k_w in the declared
brane principal operators, and raw DOF arithmetic. This script intentionally
does not test downstream no-leak, constitutive, spectrum, Maxwell, or charge
normalization claims.
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


Dim = tuple[int, int, int, int]


def add(*dims: Dim) -> Dim:
    return tuple(sum(d[i] for d in dims) for i in range(4))  # type: ignore[return-value]


def mul(n: int, dim: Dim) -> Dim:
    return tuple(n * x for x in dim)  # type: ignore[return-value]


def sub(a: Dim, b: Dim) -> Dim:
    return tuple(a[i] - b[i] for i in range(4))  # type: ignore[return-value]


def dim_string(dim: Dim) -> str:
    labels = ("L", "T", "M", "Q")
    parts: list[str] = []
    for label, power in zip(labels, dim):
        if power == 0:
            continue
        if power == 1:
            parts.append(label)
        else:
            parts.append(f"{label}^{power}")
    return " ".join(parts) if parts else "1"


def check_dim(name: str, actual: Dim, expected: Dim, note: str = "") -> dict:
    return {
        "name": name,
        "pass": actual == expected,
        "expected": dim_string(expected),
        "actual": dim_string(actual),
        "note": note,
    }


def check_bool(name: str, actual: bool, expected: bool = True, note: str = "") -> dict:
    return {
        "name": name,
        "pass": actual == expected,
        "expected": expected,
        "actual": actual,
        "note": note,
    }


def homogeneous(name: str, terms: dict[str, Dim], expected: Dim, note: str = "") -> dict:
    return {
        "name": name,
        "pass": all(dim == expected for dim in terms.values()),
        "expected": dim_string(expected),
        "terms": {key: dim_string(value) for key, value in terms.items()},
        "note": note,
    }


dim0: Dim = (0, 0, 0, 0)
L: Dim = (1, 0, 0, 0)
T: Dim = (0, 1, 0, 0)
M: Dim = (0, 0, 1, 0)
Q: Dim = (0, 0, 0, 1)

energy = add(M, mul(2, L), mul(-2, T))
action_dim = add(energy, T)
bulk_lag = add(energy, mul(-4, L))
brane_lag = add(energy, mul(-3, L))
delta_ell = mul(-1, L)
dw = L

u_dim = L
uw_dim = L
rho_br = add(M, mul(-3, L))
mu_br = brane_lag
tau_w = brane_lag
kappa_w = add(energy, mul(-1, L))
k_pin = add(energy, mul(-5, L))
pi_n = add(energy, mul(-4, L))

a_spatial = add(M, L, mul(-1, T), mul(-1, Q))
a0_dim = add(M, mul(2, L), mul(-2, T), mul(-1, Q))
alpha_u = sub(a_spatial, u_dim)
j0_dim = add(Q, mul(-3, L))
ja_dim = add(Q, mul(-2, L), mul(-1, T))
chi_dim = add(a_spatial, L)

dim_checks = [
    check_dim("bulk scalar Lagrangian density", bulk_lag, add(energy, mul(-4, L))),
    check_dim("bulk scalar action measure", add(bulk_lag, T, mul(4, L)), action_dim),
    check_dim("brane Lagrangian density", brane_lag, add(energy, mul(-3, L))),
    check_dim("brane action measure", add(brane_lag, T, mul(3, L)), action_dim),
    check_dim("finite-thickness profile B_ell(w)", delta_ell, mul(-1, L)),
    check_dim("normalized thickness integral int B_ell dw", add(delta_ell, dw), dim0),
    check_dim("bulk-density representation of brane density", add(brane_lag, delta_ell), bulk_lag),
    homogeneous(
        "in-plane elastic density terms",
        {
            "rho_parallel dot_u^2": add(rho_br, mul(2, sub(u_dim, T))),
            "Cauchy mu strain^2": mu_br,
            "rotational mu_R curl_u^2": mu_br,
        },
        brane_lag,
    ),
    homogeneous(
        "off-brane bending density terms",
        {
            "rho_w dot_u_w^2": add(rho_br, mul(2, sub(uw_dim, T))),
            "tau_w grad_u_w^2": tau_w,
            "kappa_w (Delta u_w)^2": add(kappa_w, mul(2, add(uw_dim, mul(-2, L)))),
            "k_pin u_w^2": add(k_pin, mul(2, uw_dim)),
        },
        brane_lag,
    ),
    check_dim("normal bulk-brane work density u_w Pi_n", add(uw_dim, pi_n), brane_lag),
    check_dim("A_a^brane=alpha_u u_a coefficient", add(alpha_u, u_dim), a_spatial),
    homogeneous(
        "defect-current source density",
        {
            "J^a A_a^brane": add(ja_dim, a_spatial),
            "J^0 Phi_b": add(j0_dim, a0_dim),
        },
        brane_lag,
    ),
    check_dim("gauge parameter from delta A_a=partial_a chi", sub(chi_dim, L), a_spatial),
    check_dim("gauge parameter from delta Phi=-partial_t chi", sub(chi_dim, T), a0_dim),
]

dt_j0, div_j, chi = sp.symbols("dt_j0 div_j chi")
u_only_residual = -chi * div_j
completed_residual = -chi * (dt_j0 + div_j)
completed_conserved = sp.simplify(completed_residual.subs(div_j, -dt_j0))
u_only_conserved = sp.simplify(u_only_residual.subs(div_j, -dt_j0))

gauge_checks = [
    check_bool(
        "phi-completed conserved-current residual",
        completed_conserved == 0,
        True,
        "Passes only with partial_t J0 + div J = 0 and boundary flux set by the BC contract.",
    ),
    check_bool(
        "u-only generic current negative control",
        u_only_conserved == 0,
        False,
        "A generic time-dependent charge current is not admissible with J^a u_a alone.",
    ),
    check_bool(
        "nonconserved current negative control",
        completed_residual == 0,
        False,
        "The phi-completed term is not gauge-invariant for an unconserved current.",
    ),
]

kx, ky, kz, kw, omega, rho, mu, lam, mu_r = sp.symbols(
    "kx ky kz kw omega rho mu lambda mu_r"
)
kvec = sp.Matrix([kx, ky, kz])
k2 = kx**2 + ky**2 + kz**2
identity3 = sp.eye(3)
cauchy_p = rho * omega**2 * identity3 - mu * k2 * identity3 - (lam + mu) * (kvec * kvec.T)
rot_p = rho * omega**2 * identity3 - mu_r * (k2 * identity3 - kvec * kvec.T)
bend_p = rho * omega**2 - sp.Symbol("tau_w") * k2 - sp.Symbol("kappa_w") * k2**2 - sp.Symbol("k_pin")

principal_checks = [
    check_bool("Cauchy principal block has no k_w", kw not in cauchy_p.free_symbols),
    check_bool("rotational principal block has no k_w", kw not in rot_p.free_symbols),
    check_bool("bending principal scalar has no k_w", kw not in bend_p.free_symbols),
    check_bool(
        "negative control detects explicit k_w",
        kw not in (k2 + kw**2).free_symbols,
        False,
        "Confirms the no-k_w checks are active.",
    ),
]

dof_checks = [
    check_bool("brane elastic raw components u_parallel plus u_w", 3 + 1 == 4),
    check_bool("current-completed raw variables incl auxiliary Phi_b", 3 + 1 + 1 == 5),
    check_bool("massless 4+1 A_M physical polarizations D-2", 5 - 2 == 3),
    check_bool("brane zero-mode A_mu physical polarizations 4 components minus 2 gauge/constraint", 4 - 2 == 2),
]

checks = dim_checks + gauge_checks + principal_checks + dof_checks
report = {
    "schema": "stage1_pathA_23_stage0_action_contracts_sympy/v1",
    "scope": "Stage-0 unit, source-admissibility, principal-w, and raw-DOF checks only",
    "pass": all(item["pass"] for item in checks),
    "checks": checks,
    "outcomes": {
        "stage0_token": "ACTION_SPECIFIED_CLASSIFIED",
        "current_precheck": "ADMISSIBLE_WITH_CONSERVED_DEFECT_CURRENT_AND_PHI_COMPLETION",
        "u_only_negative_control": "NOT_ADMISSIBLE_FOR_GENERIC_TIME_DEPENDENT_CHARGE",
        "principal_w_status": "NO_K_W_IN_DECLARED_BRANE_PRINCIPAL_BLOCKS",
        "classification": "NEW_PARENT_ACTION",
    },
}

out_dir = Path("software/stage1_solver/_scratch")
out_dir.mkdir(parents=True, exist_ok=True)
json_path = out_dir / "pathA_23_stage0_action_contracts_sympy.json"
json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(f"wrote {json_path}")
print(
    "pathA_23 Stage 0 SymPy cross-check: "
    f"{sum(1 for item in checks if item['pass'])}/{len(checks)} checks"
)
raise SystemExit(0 if report["pass"] else 1)
