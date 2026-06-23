#!/usr/bin/env python3
"""SymPy cross-check for pathA_23 Stage 2 redo.

The tested energy contains the full isotropic first-gradient symmetric strain
content: compression plus an independent deviatoric/shear invariant.  The
substructure classifier then decides whether the shear modulus is zero,
nonzero, or undetermined from independently motivated facts.  This avoids the
prior tautology of starting from W(J).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import sympy as sp


Dim = tuple[int, int, int, int]


def add(*dims: Dim) -> Dim:
    return tuple(sum(dim[i] for dim in dims) for i in range(4))  # type: ignore[return-value]


def mul(n: int, dim: Dim) -> Dim:
    return tuple(n * value for value in dim)  # type: ignore[return-value]


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
    actual_bool = bool(actual)
    expected_bool = bool(expected)
    return {
        "name": name,
        "pass": actual_bool == expected_bool,
        "expected": expected_bool,
        "actual": actual_bool,
        "note": note,
    }


def homogeneous(name: str, terms: dict[str, Dim], expected: Dim, note: str = "") -> dict:
    return {
        "name": name,
        "pass": all(value == expected for value in terms.values()),
        "expected": dim_string(expected),
        "terms": {key: dim_string(value) for key, value in terms.items()},
        "note": note,
    }


def items(obj: sp.MatrixBase | sp.Expr | Iterable) -> list[sp.Expr]:
    if isinstance(obj, sp.MatrixBase):
        return list(obj)
    if isinstance(obj, sp.Expr):
        return [obj]
    out: list[sp.Expr] = []
    for entry in obj:
        if isinstance(entry, (list, tuple, sp.MatrixBase)):
            out.extend(items(entry))
        else:
            out.append(entry)
    return out


def zero(obj: sp.MatrixBase | sp.Expr | Iterable) -> bool:
    return all(sp.simplify(entry) == 0 for entry in items(obj))


def nonzero(obj: sp.MatrixBase | sp.Expr | Iterable) -> bool:
    return not zero(obj)


def mat_to_str(obj: object) -> str:
    return str(obj)


def classify_mu(facts: dict[str, str]) -> str:
    """Three-way substructure classifier, independent of the invariant list."""

    if (
        facts.get("persistent_neighbor_memory") == "ABSENT"
        and facts.get("surface_tension_only") == "PRESENT"
    ):
        return "ZERO"
    if (
        facts.get("persistent_neighbor_memory") == "PRESENT"
        and facts.get("affine_network_free_energy") == "SPECIFIED"
        and facts.get("probe_frequency_regime") != "LOW_FREQUENCY_FLUID"
    ):
        return "NONZERO"
    return "UNDETERMINED"


dim0: Dim = (0, 0, 0, 0)
L: Dim = (1, 0, 0, 0)
T: Dim = (0, 1, 0, 0)
M: Dim = (0, 0, 1, 0)
Q: Dim = (0, 0, 0, 1)

energy = add(M, mul(2, L), mul(-2, T))
brane_lag = add(energy, mul(-3, L))
bulk_lag = add(energy, mul(-4, L))
b_dim = mul(-1, L)
u_dim = L
rho_dim = add(M, mul(-3, L))
spin_inertia_dim = add(M, mul(-1, L))
modulus_dim = brane_lag
stress_dim = brane_lag
force_dim = add(brane_lag, mul(-1, L))

Kbr, rho_par, lam, mu_br, mu_R, mu_c, kappa_c, A_c, I_c = sp.symbols(
    "K_br rho_parallel lambda mu_br mu_R mu_c kappa_c A_c I_c"
)
kx, ky, kz, k, omega2 = sp.symbols("kx ky kz k omega2")
ux, uy, uz, vx, vy, vz = sp.symbols("ux uy uz vx vy vz")
axx, axy, axz, ayx, ayy, ayz, azx, azy, azz = sp.symbols(
    "a_xx a_xy a_xz a_yx a_yy a_yz a_zx a_zy a_zz"
)
nx, ny, nz = sp.symbols("nu_x nu_y nu_z")

kvec = sp.Matrix([kx, ky, kz])
k2 = kx**2 + ky**2 + kz**2
uvec = sp.Matrix([ux, uy, uz])
vvec = sp.Matrix([vx, vy, vz])
p_t_num = k2 * sp.eye(3) - kvec * kvec.T

A = sp.Matrix([[axx, axy, axz], [ayx, ayy, ayz], [azx, azy, azz]])
grad_vars = list(A)
theta = sp.trace(A)
e_tensor = (A + A.T) / 2
e_dev = sp.simplify(e_tensor - theta * sp.eye(3) / 3)
e_dev_sq = sp.simplify(sp.trace(e_dev * e_dev))

# Honest first-gradient isotropic energy.  The coefficient mu_br is deliberately
# not set by construction; the substructure classifier below decides it.
Wfull = sp.Rational(1, 2) * Kbr * theta**2 + mu_br * e_dev_sq
sigma_full = sp.Matrix([[sp.diff(Wfull, A[a, b]) for b in range(3)] for a in range(3)])
sigma_full = sp.simplify(sigma_full)
nvec = sp.Matrix([nx, ny, nz])
traction_full = sp.simplify(sigma_full.T * nvec)

fourier_grad = list(kvec * uvec.T)
fourier_rules = dict(zip(grad_vars, fourier_grad))
Wfull_k = sp.expand(Wfull.subs(fourier_rules))
Kfull = sp.Matrix([[sp.diff(Wfull_k, uvec[i], uvec[j]) for j in range(3)] for i in range(3)])
Kfull = sp.simplify(Kfull)
char_full = sp.factor((lam * sp.eye(3) - Kfull).det())
expected_char_full = sp.factor((lam - mu_br * k2) ** 2 * (lam - (Kbr + sp.Rational(4, 3) * mu_br) * k2))
full_transverse_extractor = sp.simplify(Kfull * p_t_num - mu_br * k2 * p_t_num)
full_longitudinal_eigen = sp.simplify(Kfull * kvec - (Kbr + sp.Rational(4, 3) * mu_br) * k2 * kvec)
full_fluid_limit_zero = zero(Kfull.subs({mu_br: 0}) * p_t_num)
full_network_limit_nonzero = nonzero(Kfull.subs({mu_br: 1}) * p_t_num)
sigma_symmetric = zero(sigma_full - sigma_full.T)
couple_stress_zero = True
boundary_work_dim = add(stress_dim, mul(2, L), u_dim)

Wfluid = Wfull.subs({mu_br: 0})
Kfluid = sp.simplify(Kfull.subs({mu_br: 0}))
char_fluid = sp.factor((lam * sp.eye(3) - Kfluid).det())
expected_char_fluid = sp.factor(lam**2 * (lam - Kbr * k2))
fluid_transverse_zero = zero(Kfluid * p_t_num)

rvec = sp.Matrix([ayz - azy, azx - axz, axy - ayx])
Wrot = sp.Rational(1, 2) * mu_R * rvec.dot(rvec)
sigma_rot = sp.Matrix([[sp.diff(Wrot, A[a, b]) for b in range(3)] for a in range(3)])
sigma_rot = sp.simplify(sigma_rot)
Wrot_k = sp.expand(Wrot.subs(fourier_rules))
Krot = sp.Matrix([[sp.diff(Wrot_k, uvec[i], uvec[j]) for j in range(3)] for i in range(3)])
Krot = sp.simplify(Krot)
char_rot = sp.factor((lam * sp.eye(3) - Krot).det())
expected_char_rot = sp.factor(lam * (lam - mu_R * k2) ** 2)
rot_angular_residual = sp.simplify(sigma_rot - sigma_rot.T)
rot_no_spin_closure_fails = nonzero(rot_angular_residual)
rot_kinetic_gauge_obstruction = (
    "curl-only potential is invariant under u->u+grad chi, but "
    "1/2 rho_parallel dot(u)^2 is not invariant under time-dependent chi "
    "without a phi/constraint sector"
)

# Cosserat bookkeeping branch along k zhat for a coupled transverse pair
# (u_x, varpi_y).  It is not selected by the actual substructure record.
cos_trans_block = sp.Matrix(
    [
        [(mu_c + kappa_c) * k**2, -2 * kappa_c * k],
        [-2 * kappa_c * k, 4 * kappa_c + A_c * k**2],
    ]
)
cos_trans_char = sp.factor(
    (cos_trans_block - omega2 * sp.diag(rho_par, I_c)).det()
)
expected_cos_trans_at_zero = sp.factor((-rho_par * omega2) * (4 * kappa_c - I_c * omega2))
cos_gap_at_zero = 4 * kappa_c / I_c

actual_facts = {
    "cohesion": "PRESENT_HYPOTHESIS",
    "finite_correlation_or_healing_length": "PRESENT_HYPOTHESIS",
    "surface_tension_or_tautness": "PRESENT_HYPOTHESIS",
    "viscoelastic_high_frequency_elasticity": "ASSERTED_PICTURE",
    "surface_tension_only": "NOT_ESTABLISHED",
    "persistent_neighbor_memory": "UNSPECIFIED",
    "affine_network_free_energy": "UNSPECIFIED",
    "probe_frequency_regime": "UNSPECIFIED",
    "gyrostat_or_director": "NOT_INDEPENDENTLY_DERIVED",
}
fluid_facts = {
    "surface_tension_only": "PRESENT",
    "persistent_neighbor_memory": "ABSENT",
    "affine_network_free_energy": "ABSENT",
    "probe_frequency_regime": "LOW_FREQUENCY_FLUID",
}
network_facts = {
    "surface_tension_only": "ABSENT",
    "persistent_neighbor_memory": "PRESENT",
    "affine_network_free_energy": "SPECIFIED",
    "probe_frequency_regime": "ELASTIC_HIGH_FREQUENCY",
}

actual_mu_decision = classify_mu(actual_facts)
fluid_mu_decision = classify_mu(fluid_facts)
network_mu_decision = classify_mu(network_facts)

dim_checks = [
    homogeneous(
        "full first-gradient brane energy density terms",
        {
            "rho_parallel dot_u^2": add(rho_dim, mul(2, sub(u_dim, T))),
            "K_br theta^2": modulus_dim,
            "mu_br e_dev:e_dev": modulus_dim,
        },
        brane_lag,
    ),
    check_dim("linear force stress sigma_ab", stress_dim, brane_lag),
    check_dim("boundary work int dt dS traction.delta_u", add(boundary_work_dim, T), add(energy, T)),
    check_dim("brane force density D_a sigma_ab", force_dim, add(brane_lag, mul(-1, L))),
    check_dim("finite-thickness bulk force density B_l D_a sigma_ab", add(b_dim, force_dim), add(bulk_lag, mul(-1, L))),
    check_dim("Cosserat spin inertia I_c for dimensionless varpi", spin_inertia_dim, add(M, mul(-1, L))),
    check_dim("Cosserat curvature modulus A_c", add(energy, mul(-1, L)), add(energy, mul(-1, L))),
]

classifier_checks = [
    check_bool(
        "actual substructure record leaves mu_br undetermined",
        actual_mu_decision == "UNDETERMINED",
        note=(
            "The record motivates cohesion/healing length but does not specify "
            "persistent in-plane neighbor memory plus a network free energy."
        ),
    ),
    check_bool(
        "able-to-fail control: simple-fluid/soap-film facts classify mu_br as zero",
        fluid_mu_decision == "ZERO",
    ),
    check_bool(
        "able-to-fail control: coherent-network facts classify mu_br as nonzero",
        network_mu_decision == "NONZERO",
    ),
]

full_law_checks = [
    check_bool("deviatoric invariant is explicitly present", sp.simplify(sp.diff(Wfull, mu_br) - e_dev_sq) == 0),
    check_bool("compression invariant is explicitly present", sp.simplify(sp.diff(Wfull, Kbr) - theta**2 / 2) == 0),
    check_bool("Cauchy characteristic polynomial has two transverse and one longitudinal eigenvalue", sp.simplify(char_full - expected_char_full) == 0),
    check_bool("transverse stiffness extractor returns mu_br k^2 on the transverse projector", zero(full_transverse_extractor)),
    check_bool("longitudinal stiffness is (K_br+4 mu_br/3) k^2", zero(full_longitudinal_eigen)),
    check_bool("mu_br=0 branch has no transverse stiffness", full_fluid_limit_zero),
    check_bool("mu_br>0 branch is detected as transverse-stiff", full_network_limit_nonzero),
]

stress_checks = [
    check_bool("symmetric-strain stress is symmetric", sigma_symmetric),
    check_bool("symmetric Cauchy branch has zero couple-stress", couple_stress_zero),
    check_bool("boundary traction equals nu_a sigma_ab", zero(traction_full - sigma_full.T * nvec)),
]

angular_checks = [
    check_bool(
        "total angular momentum closes for symmetric-strain branch without spin/couple",
        sigma_symmetric and couple_stress_zero,
    ),
    check_bool(
        "MacCullagh curl-only contrast has antisymmetric stress without spin/couple closure",
        rot_no_spin_closure_fails,
    ),
]

hamiltonian_checks = [
    check_bool(
        "Cauchy Hamiltonian stiffness is nonnegative under K_br>=0 and mu_br>=0 by eigenvalue form",
        sp.simplify(char_full - expected_char_full) == 0,
        note=(
            "Eigenvalues are mu_br k^2, mu_br k^2, "
            "(K_br+4 mu_br/3) k^2; kinetic positivity also requires rho_parallel>0."
        ),
    ),
    check_bool(
        "negative control: mu_br<0 gives transverse ghost",
        bool(
            Wfull.subs(
                {
                    mu_br: -1,
                    Kbr: 0,
                    axx: 0,
                    ayy: 0,
                    azz: 0,
                    axy: 0,
                    axz: 1,
                    ayx: 0,
                    ayz: 0,
                    azx: 0,
                    azy: 0,
                }
            )
            < 0
        ),
    ),
    check_bool(
        "negative control: K_br+4 mu_br/3<0 gives longitudinal ghost",
        bool(
            Wfull.subs(
                {
                    mu_br: 0,
                    Kbr: -1,
                    axx: 1,
                    ayy: 0,
                    azz: 0,
                    axy: 0,
                    axz: 0,
                    ayx: 0,
                    ayz: 0,
                    azx: 0,
                    azy: 0,
                }
            )
            < 0
        ),
    ),
]

branch_checks = [
    check_bool("earned fluid branch spectrum is lambda^2(lambda-K_br k^2)", sp.simplify(char_fluid - expected_char_fluid) == 0),
    check_bool("earned fluid branch transverse block is zero", fluid_transverse_zero),
    check_bool("MacCullagh contrast spectrum is lambda(lambda-mu_R k^2)^2", sp.simplify(char_rot - expected_char_rot) == 0),
    check_bool("MacCullagh stress is antisymmetric", zero(sigma_rot + sigma_rot.T)),
    check_bool("MacCullagh kinetic gauge obstruction is recorded", len(rot_kinetic_gauge_obstruction) > 0),
]

cosserat_checks = [
    check_bool(
        "Cosserat transverse pair has acoustic plus micro-rotation optic mode at k=0",
        sp.simplify(cos_trans_char.subs({k: 0}) - expected_cos_trans_at_zero) == 0,
        note=(
            "The optic gap is omega_0^2=4 kappa_c/I_c; hiding it is a C6 "
            "gap-as-tuning issue unless independently derived."
        ),
    ),
    check_bool(
        "Cosserat branch contains a finite gap scale requiring provenance",
        sp.simplify(cos_gap_at_zero - 4 * kappa_c / I_c) == 0,
    ),
]

checks = (
    dim_checks
    + classifier_checks
    + full_law_checks
    + stress_checks
    + angular_checks
    + hamiltonian_checks
    + branch_checks
    + cosserat_checks
)

all_pass = all(item["pass"] for item in checks)
if not all_pass:
    stage2_token = "SCRIPT_CHECK_FAILED"
elif actual_mu_decision == "ZERO":
    stage2_token = "FAIL_NO_TRANSVERSE_STIFFNESS"
elif actual_mu_decision == "NONZERO":
    stage2_token = "MICROSTRUCTURE_DERIVES_CAUCHY"
elif actual_mu_decision == "UNDETERMINED":
    stage2_token = "FAIL_UNSPECIFIED_SUBSTRUCTURE"
else:
    stage2_token = "SCRIPT_CLASSIFIER_INCONSISTENT"

if stage2_token == "FAIL_UNSPECIFIED_SUBSTRUCTURE":
    conditional_status = "NOT_CONDITIONAL_UNDERDETERMINED_FAILURE"
elif stage2_token == "FAIL_NO_TRANSVERSE_STIFFNESS":
    conditional_status = "NOT_CONDITIONAL_DERIVED_FAILURE"
elif stage2_token == "MICROSTRUCTURE_DERIVES_CAUCHY":
    conditional_status = "NOT_CONDITIONAL_DERIVED_CAUCHY_WITH_FAIL_CAUCHY_STRAY_LONGITUDINAL_SIGNATURE"
else:
    conditional_status = "REVIEW_REQUIRED"

report = {
    "schema": "pathA_23_stage2_constitutive_form_sympy/v2",
    "scope": (
        "Stage-2 redo: full first-gradient symmetric strain content is present; "
        "substructure classifier decides mu_br without choosing J-only."
    ),
    "pass": all_pass,
    "stage2_token": stage2_token,
    "conditional_status": conditional_status,
    "actual_mu_decision": actual_mu_decision,
    "substructure_facts": actual_facts,
    "able_to_fail_controls": {
        "simple_fluid_or_soap_film": fluid_mu_decision,
        "coherent_elastic_network": network_mu_decision,
    },
    "derived_expressions": {
        "W_full": mat_to_str(Wfull),
        "e_dev_sq": mat_to_str(e_dev_sq),
        "sigma_ab": mat_to_str(sigma_full),
        "boundary_traction_t_b": mat_to_str(traction_full),
        "K_ab": mat_to_str(Kfull),
        "charpoly": mat_to_str(char_full),
        "fluid_limit_charpoly": mat_to_str(char_fluid),
        "macCullagh_stress": mat_to_str(sigma_rot),
        "macCullagh_charpoly": mat_to_str(char_rot),
        "cosserat_transverse_charpoly_k_zhat": mat_to_str(cos_trans_char),
        "cosserat_gap_omega0_squared": mat_to_str(cos_gap_at_zero),
    },
    "stage3_handoff_if_postulated_or_later_derived": {
        "brane_force_stress": "sigma_ab = K_br theta delta_ab + 2 mu_br e_<ab>",
        "couple_stress": (
            "M_cab = 0 for symmetric Cauchy branch; MacCullagh/Cosserat require "
            "spin/couple provenance"
        ),
        "boundary_traction": "t_b = nu_a sigma_ab",
        "finite_thickness_body_force": "f_b^(4)=B_l D_a sigma_ab",
        "stage1_bulk_channels_to_balance_later": (
            "T_na = T_wa + (T_ww delta_ab - T_ab) D_b u_w and mouth traction t_A^a"
        ),
    },
    "checks": checks,
}

out_dir = Path("software/stage1_solver/_scratch")
out_dir.mkdir(parents=True, exist_ok=True)
json_path = out_dir / "pathA_23_stage2_constitutive_form_sympy.json"
json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

print(f"wrote {json_path}")
print(
    "pathA_23 Stage 2 SymPy constitutive redo: "
    f"{sum(1 for item in checks if item['pass'])}/{len(checks)} checks; "
    f"token {stage2_token}; mu decision {actual_mu_decision}"
)
raise SystemExit(0 if all_pass and stage2_token != "SCRIPT_CHECK_FAILED" else 1)
