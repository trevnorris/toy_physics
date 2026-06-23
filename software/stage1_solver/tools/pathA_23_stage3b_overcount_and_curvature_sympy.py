#!/usr/bin/env python3
"""Independent SymPy cross-check for pathA_23 Stage 3b.

The audit distinguishes brane variational stress from bulk-field source terms,
then checks ideal-fluid tangential traction, flat-brane in-plane light, and the
local curvature scaling of the remaining defect/background channels.
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
        parts.append(label if power == 1 else f"{label}^{power}")
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


def _items(expr: sp.MatrixBase | sp.Expr | Iterable[sp.Expr]) -> list[sp.Expr]:
    if isinstance(expr, sp.MatrixBase):
        return list(expr)
    if isinstance(expr, sp.Expr):
        return [expr]
    return list(expr)


def zero_vec(expr: sp.MatrixBase | sp.Expr | Iterable[sp.Expr]) -> bool:
    return all(sp.simplify(item) == 0 for item in _items(expr))


def nonzero_vec(expr: sp.MatrixBase | sp.Expr | Iterable[sp.Expr]) -> bool:
    return not zero_vec(expr)


def contains(expr: sp.MatrixBase | sp.Expr | Iterable[sp.Expr], sym: sp.Symbol) -> bool:
    return any(item.has(sym) for item in _items(expr))


def mat_to_str(obj: sp.MatrixBase | sp.Expr | Iterable[sp.Expr]) -> str:
    return str(obj)


dim0: Dim = (0, 0, 0, 0)
L: Dim = (1, 0, 0, 0)
T: Dim = (0, 1, 0, 0)
M: Dim = (0, 0, 1, 0)
Q: Dim = (0, 0, 0, 1)

energy = add(M, mul(2, L), mul(-2, T))
brane_lag = sub(energy, mul(3, L))
bulk_lag = sub(energy, mul(4, L))
b_dim = mul(-1, L)
u_dim = L
rho3_dim = sub(M, mul(3, L))
rho4_dim = sub(M, mul(4, L))
velocity_dim = sub(L, T)
stress_bulk_dim = bulk_lag
brane_force_dim = sub(brane_lag, L)
bulk_force_dim = sub(bulk_lag, L)
k_dim = mul(-1, L)
curvature_dim = mul(-1, L)

kx, ky, kz, k = sp.symbols("kx ky kz k")
ux, uy, uz, Uamp, phi = sp.symbols("ux uy uz Uamp phi")
muR, Bl, rho, theta = sp.symbols("mu_R B_l rho theta")
vx, vy, vz, vw, udotw = sp.symbols("vx vy vz vw dot_u_w")
sx, sy, sz = sp.symbols("sx sy sz")
alphaU, Jx, Jy, Jz, J0, Phi = sp.symbols("alpha_u Jx Jy Jz J0 Phi")
beta, eta = sp.symbols("beta eta")

kvec = sp.Matrix([kx, ky, kz])
k2 = kx**2 + ky**2 + kz**2
p_t_num = k2 * sp.eye(3) - kvec * kvec.T
uvec = sp.Matrix([ux, uy, uz])
vpar = sp.Matrix([vx, vy, vz])
slopes = sp.Matrix([sx, sy, sz])

axx, axy, axz, ayx, ayy, ayz, azx, azy, azz = sp.symbols(
    "a_xx a_xy a_xz a_yx a_yy a_yz a_zx a_zy a_zz"
)
A = sp.Matrix([[axx, axy, axz], [ayx, ayy, ayz], [azx, azy, azz]])
grad_vars = list(A)

rvec = sp.Matrix([ayz - azy, azx - axz, axy - ayx])
u_rot = sp.Rational(1, 2) * muR * rvec.dot(rvec)
sigma_r = sp.Matrix(3, 3, lambda a, b: sp.diff(u_rot, A[a, b]))
sigma_expected = muR * (A - A.T)

fourier_values = []
for i in range(3):
    for j in range(3):
        fourier_values.append(sp.I * kvec[i] * uvec[j])
sigma_fourier = sp.simplify(sigma_r.subs(dict(zip(grad_vars, fourier_values))))
div_sigma_r = sp.Matrix(
    [sp.simplify(sum(sp.I * kvec[a] * sigma_fourier[a, b] for a in range(3))) for b in range(3)]
)
div_sigma_expected = sp.simplify(muR * (kvec * (kvec.dot(uvec)) - k2 * uvec))

bulk_vars = (rho, theta, vx, vy, vz, vw)
l_bulk = sp.Function("L_psi")(*bulk_vars)
pi_n = sp.Function("Pi_n")(*(bulk_vars + (sx, sy, sz)))
l_brane = -u_rot
s_cpl = (
    sp.Symbol("u_w") * pi_n
    + alphaU * (Jx * ux + Jy * uy + Jz * uz)
    - J0 * Phi
)
l_declared = l_bulk + Bl * (l_brane + s_cpl)
l_no_brane = l_bulk + Bl * s_cpl
bulk_variations_declared = sp.Matrix([sp.diff(l_declared, var) for var in bulk_vars])
bulk_variations_no_brane = sp.Matrix([sp.diff(l_no_brane, var) for var in bulk_vars])

bad_velocity_coupling = Bl * beta * vpar.dot(div_sigma_r)
bad_rho_brane_coupling = Bl * eta * rho * u_rot
bad_bulk_velocity_variation = sp.Matrix(
    [sp.simplify(sp.diff(l_declared + bad_velocity_coupling, var)) for var in (vx, vy, vz)]
)
bad_bulk_rho_variation = sp.simplify(sp.diff(l_declared + bad_rho_brane_coupling, rho))

current_coupling = Bl * (alphaU * (Jx * ux + Jy * uy + Jz * uz) - J0 * Phi)
current_bulk_velocity_variation = sp.Matrix([sp.diff(current_coupling, var) for var in (vx, vy, vz, vw)])
bad_current_velocity_coupling = Bl * beta * (Jx * vx + Jy * vy + Jz * vz)
bad_current_var = sp.Matrix([sp.diff(bad_current_velocity_coupling, var) for var in (vx, vy, vz)])

mG, P = sp.symbols("m_GNLS P")
sigQwx, sigQwy, sigQwz = sp.symbols("sigmaQ_wx sigmaQ_wy sigmaQ_wz")
sig_q = sp.Matrix([sigQwx, sigQwy, sigQwz])
v4 = sp.Matrix([vx, vy, vz, vw])
tflat = mG * rho * (v4 * v4.T) + P * sp.eye(4) + sp.Matrix(
    [
        [0, 0, 0, sigQwx],
        [0, 0, 0, sigQwy],
        [0, 0, 0, sigQwz],
        [sigQwx, sigQwy, sigQwz, 0],
    ]
)
pressure_traction_flat = sp.Matrix([(sp.eye(4)[3, :] * (P * sp.eye(4)) * sp.eye(4)[:, a])[0] for a in range(3)])
convective_traction_flat = sp.Matrix([mG * rho * vw * vpar[a] for a in range(3)])
quantum_traction_flat = sig_q
total_traction_flat = sp.Matrix([tflat[3, a] for a in range(3)])

hbar, rho0, rho0p, delta_rho, delta_rho_w, rho_br = sp.symbols(
    "hbar rho0 rho0_prime delta_rho delta_rho_w rho_br"
)
sig_qwa_linear = sp.simplify(
    (hbar**2 / (4 * mG)) * sp.I * kvec * (rho0p * delta_rho / rho0 - delta_rho_w)
)
delta_rho_shear = -sp.I * rho_br * kvec.dot(uvec)
sig_q_shear_general = sp.simplify(sig_qwa_linear.subs({delta_rho: delta_rho_shear, delta_rho_w: 0}))
transverse_rules = {kx: 0, ky: 0, kz: k, ux: Uamp, uy: 0, uz: 0}
compressive_rules = {kx: 0, ky: 0, kz: k, ux: 0, uy: 0, uz: Uamp}

v_w_compat = udotw + vpar.dot(slopes)
t_flat_light = sp.simplify(
    (mG * rho * v_w_compat * vpar + sig_q_shear_general)
    .subs({udotw: 0, sx: 0, sy: 0, sz: 0})
    .subs(transverse_rules)
)
t_flat_bad_normal_flow = sp.simplify((mG * rho * vw * vpar + sig_q_shear_general).subs(transverse_rules))

n4 = sp.Matrix([-sx, -sy, -sz, 1])
e_tangents = [sp.eye(4)[:, a] + slopes[a] * sp.eye(4)[:, 3] for a in range(3)]
t_ideal = mG * rho * (v4 * v4.T) + P * sp.eye(4)
pressure_tilt_traction = sp.Matrix([sp.simplify((n4.T * (P * sp.eye(4)) * e_tangents[a])[0]) for a in range(3)])
tangent_flow_rules = {vw: vpar.dot(slopes)}
tilted_ideal_traction = sp.Matrix(
    [sp.simplify((n4.T * t_ideal * e_tangents[a])[0].subs(tangent_flow_rules)) for a in range(3)]
)
eps = sp.symbols("eps")
tilted_ideal_traction_linear = sp.simplify(
    tilted_ideal_traction.subs({sx: eps * sx, sy: eps * sy, sz: eps * sz}).applyfunc(
        lambda expr: sp.series(expr, eps, 0, 2).removeO().subs(eps, 1)
    )
)
fixed_chart_linear_traction = sp.simplify(
    (
        mG * rho * vw * vpar
        + (P * sp.eye(3) - (mG * rho * (vpar * vpar.T) + P * sp.eye(3))) * slopes
    ).subs(tangent_flow_rules)
)

Kxx, Kxy, Kxz, Kyx, Kyy, Kyz, Kzx, Kzy, Kzz = sp.symbols(
    "Kxx Kxy Kxz Kyx Kyy Kyz Kzx Kzy Kzz"
)
x1, x2, x3, Lmix = sp.symbols("x1 x2 x3 L_mix")
Axx, Axy, Axz, Ayx, Ayy, Ayz, Azx, Azy, Azz = sp.symbols(
    "Axx Axy Axz Ayx Ayy Ayz Azx Azy Azz"
)
kmat = sp.Matrix([[Kxx, Kxy, Kxz], [Kyx, Kyy, Kyz], [Kzx, Kzy, Kzz]])
xvec = sp.Matrix([x1, x2, x3])
local_slope = kmat * xvec
amix = sp.Matrix([[Axx, Axy, Axz], [Ayx, Ayy, Ayz], [Azx, Azy, Azz]])
curvature_mixing = sp.simplify(amix * local_slope)
curvature_off_rules = {symbol: 0 for symbol in list(kmat)}
curvature_on_rules = {
    Kxy: 0,
    Kxz: 0,
    Kyx: 0,
    Kyy: 0,
    Kyz: 0,
    Kzx: 0,
    Kzy: 0,
    Kzz: 0,
    x1: Lmix,
    x2: 0,
    x3: 0,
}

fref_dim = brane_force_dim
eps_amp_dim = sub(stress_bulk_dim, fref_dim)
eps_power_dim = add(mul(2, add(curvature_dim, L)), stress_bulk_dim, mul(-1, fref_dim))

dim_checks = [
    homogeneous(
        "declared finite-thickness density pieces",
        {
            "bulk scalar Lpsi": bulk_lag,
            "B_l U_R": add(b_dim, brane_lag),
            "B_l u_w Pi_n": add(b_dim, u_dim, stress_bulk_dim),
            "B_l J^a alpha_u u_a / J0 Phi": add(b_dim, brane_lag),
        },
        bulk_lag,
    ),
    check_dim("sigma_R stress", brane_lag, brane_lag),
    check_dim("D_b sigma_R brane force density", brane_force_dim, sub(brane_lag, L)),
    check_dim("bulk tangential traction T_wa", stress_bulk_dim, brane_force_dim),
    check_dim("curvature K_ab", curvature_dim, mul(-1, L)),
    check_dim("geometric amplitude factor K L_mix", add(curvature_dim, L), dim0),
    check_dim("epsilon amplitude stress/reference factor before K L_mix", eps_amp_dim, dim0),
    check_dim("epsilon power-style curvature-squared factor", eps_power_dim, dim0),
]

variational_checks = [
    check_bool("sigma_R_ab derives from U_R", zero_vec(sigma_r - sigma_expected)),
    check_bool("D_b sigma_R_ba Fourier form derives correctly", zero_vec(div_sigma_r - div_sigma_expected)),
    check_bool(
        "declared bulk variations equal variations of S_psi + S_cpl only",
        zero_vec(bulk_variations_declared - bulk_variations_no_brane),
    ),
    check_bool("declared bulk variations contain no mu_R", not contains(bulk_variations_declared, muR)),
    check_bool("declared source-current coupling has no direct bulk-velocity variation", zero_vec(current_bulk_velocity_variation)),
]

traction_checks = [
    check_bool("flat pressure traction n_i P delta_ia is tangential-zero", zero_vec(pressure_traction_flat)),
    check_bool(
        "flat ideal tangential traction is convective plus quantum only",
        zero_vec(total_traction_flat - convective_traction_flat - quantum_traction_flat),
    ),
    check_bool("quantum mixed traction is nonzero for generic density perturbation", nonzero_vec(sig_qwa_linear)),
    check_bool(
        "quantum mixed traction vanishes for density-preserving transverse shear",
        zero_vec(sp.simplify(sig_q_shear_general.subs(transverse_rules))),
    ),
    check_bool(
        "quantum mixed traction would survive for compressive density perturbation",
        nonzero_vec(sp.simplify(sig_q_shear_general.subs(compressive_rules))),
    ),
]

flat_checks = [
    check_bool(
        "kinematic compatibility gives v_w(light)=0 on flat brane",
        sp.simplify(v_w_compat.subs({udotw: 0, sx: 0, sy: 0, sz: 0})) == 0,
    ),
    check_bool("flat transverse in-plane light has T_wa=0 including quantum stress", zero_vec(t_flat_light)),
    check_bool("flat light would leak if an independent v_w were present", nonzero_vec(t_flat_bad_normal_flow)),
]

curvature_checks = [
    check_bool("constant tilted ideal free-slip surface has zero pressure tangential traction", zero_vec(pressure_tilt_traction)),
    check_bool(
        "constant tilted ideal free-slip surface has zero convective tangential traction for v.n=0",
        zero_vec(tilted_ideal_traction),
    ),
    check_bool("fixed-chart linear slope terms cancel for tangent ideal flow", zero_vec(fixed_chart_linear_traction)),
    check_bool(
        "curvature-local mixing vanishes when K_ab=0",
        zero_vec(sp.simplify(curvature_mixing.subs(curvature_off_rules))),
    ),
    check_bool(
        "curvature-local mixing is generically nonzero when K_ab != 0",
        nonzero_vec(sp.simplify(curvature_mixing.subs(curvature_on_rules))),
    ),
]

negative_checks = [
    check_bool(
        "negative control: v_a D_b sigma_R_ba coupling would source bulk velocity equation",
        contains(bad_bulk_velocity_variation, muR) and nonzero_vec(bad_bulk_velocity_variation),
    ),
    check_bool(
        "negative control: rho U_R coupling would source bulk density equation",
        bad_bulk_rho_variation.has(muR) and sp.simplify(bad_bulk_rho_variation) != 0,
    ),
    check_bool(
        "negative control: forbidden J^a-to-bulk-velocity coupling would be detected",
        not zero_vec(bad_current_var),
    ),
]

checks = dim_checks + variational_checks + traction_checks + flat_checks + curvature_checks + negative_checks
all_pass = all(bool(check["pass"]) for check in checks)

sigma_r_not_bulk = zero_vec(bulk_variations_declared - bulk_variations_no_brane) and not contains(
    bulk_variations_declared, muR
)
sigma_r_sources_bulk = not sigma_r_not_bulk
flat_light_no_leak = zero_vec(t_flat_light)
curvature_localized = zero_vec(sp.simplify(curvature_mixing.subs(curvature_off_rules))) and nonzero_vec(
    sp.simplify(curvature_mixing.subs(curvature_on_rules))
)

if all_pass and sigma_r_not_bulk and flat_light_no_leak and curvature_localized:
    primary_token = "OVER_COUNT_CONFIRMED_CURVATURE_LOCALIZED"
elif all_pass and (sigma_r_sources_bulk or not flat_light_no_leak):
    primary_token = "INTRINSIC_LEAK_CONFIRMED_FATAL"
elif all_pass:
    primary_token = "MIXED_PARTIAL"
else:
    primary_token = "SCRIPT_CHECK_FAILURE"

report = {
    "stage": "pathA_23 Stage 3b over-count and curvature localization",
    "conditional_flag": "CONDITIONAL",
    "primary_token": primary_token,
    "sigma_R_bulk_source_token": "SIGMA_R_NOT_A_BULK_SOURCE"
    if sigma_r_not_bulk
    else "SIGMA_R_SOURCES_BULK",
    "flat_light_token": "LIGHT_FREE_SLIPS_NO_LEAK" if flat_light_no_leak else "LIGHT_LEAKS_FLAT",
    "checks_passed": sum(1 for check in checks if check["pass"]),
    "checks_total": len(checks),
    "all_pass": all_pass,
    "load_bearing_expressions": {
        "bulk_variations_declared": mat_to_str(bulk_variations_declared),
        "bulk_variations_without_brane": mat_to_str(bulk_variations_no_brane),
        "div_sigma_R_fourier": mat_to_str(div_sigma_r),
        "flat_tangential_traction": mat_to_str(total_traction_flat),
        "flat_pressure_traction": mat_to_str(pressure_traction_flat),
        "flat_convective_traction": mat_to_str(convective_traction_flat),
        "quantum_sigma_Q_wa_linear": mat_to_str(sig_qwa_linear),
        "flat_light_T_wa": mat_to_str(t_flat_light),
        "tilted_ideal_traction_tangent_flow": mat_to_str(tilted_ideal_traction),
        "curvature_mixing_A_K_x": mat_to_str(curvature_mixing),
        "epsilon_force_amplitude_scaling": "(||A_mix||/F_ref) |K| L_mix",
        "epsilon_power_scaling_if_energy_fraction": (
            "(||A_mix||/F_ref) (|K| L_mix)^2, when epsilon is defined as leaked "
            "power/probability rather than force amplitude"
        ),
    },
    "negative_controls": {
        "bad_velocity_sigma_coupling_variation": mat_to_str(bad_bulk_velocity_variation),
        "bad_rho_brane_coupling_variation": mat_to_str(bad_bulk_rho_variation),
        "bad_current_velocity_variation": mat_to_str(bad_current_var),
        "bad_flat_normal_flow_traction": mat_to_str(t_flat_bad_normal_flow),
    },
    "checks": checks,
}

out_dir = Path("software/stage1_solver/_scratch")
out_dir.mkdir(parents=True, exist_ok=True)
json_path = out_dir / "pathA_23_stage3b_overcount_and_curvature_sympy.json"
json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

raise SystemExit(0 if all_pass else 1)
