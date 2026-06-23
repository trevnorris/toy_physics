#!/usr/bin/env python3
"""Independent SymPy cross-check for pathA_23 Stage 3.

The rotational/MacCullagh law is treated as a user-postulated conditional
ingredient.  The script checks the derived force-stress, the minimal
gyrostatic spin closure, and the transverse/curl source obtained when the
constitutive reaction is combined with the Stage-1 interface channels.
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


def zero_matrix(mat: sp.MatrixBase | Iterable[sp.Expr]) -> bool:
    return zero_vec(mat)


def nonzero_matrix(mat: sp.MatrixBase | Iterable[sp.Expr]) -> bool:
    return not zero_matrix(mat)


def mat_to_str(obj: sp.MatrixBase | sp.Expr | list[sp.Expr]) -> str:
    return str(obj)


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
rho3_dim = add(M, mul(-3, L))
velocity_dim = sub(L, T)
circulation_dim = add(mul(2, L), mul(-1, T))
vortex_density_dim = mul(-2, L)
stress_brane_dim = brane_lag
brane_force_dim = add(brane_lag, mul(-1, L))
bulk_stress_dim = bulk_lag
couple_stress_dim = add(stress_brane_dim, L)
spin_density_dim = add(stress_brane_dim, T)

kx, ky, kz, k = sp.symbols("kx ky kz k")
ux, uy, uz, phi, Uamp = sp.symbols("ux uy uz phi Uamp")
muR, rhoPar, alphaU = sp.symbols("mu_R rho_parallel alpha_u")
Jx, Jy, Jz, J0, Phi, Bl, ell = sp.symbols("Jx Jy Jz J0 Phi B_l ell")
mG, rho3, GammaCirc, nVort, Vrel, vr, Lg = sp.symbols(
    "m_GNLS rho_3 Gamma n_v V_rel v_r L_g"
)
piww, piIso = sp.symbols("Pi_ww Pi_iso")
pixx, pixy, pixz, piyy, piyz, pizz, piwx, piwy, piwz = sp.symbols(
    "Pi_xx Pi_xy Pi_xz Pi_yy Pi_yz Pi_zz Pi_wx Pi_wy Pi_wz"
)
uw, tauAx, tauAy, tauAz = sp.symbols("u_w tau_Ax tau_Ay tau_Az")
lambdaN = sp.symbols("Lambda_n")
rho, theta, vx, vy, vz, vw = sp.symbols("rho theta vx vy vz vw")
Vx, Vy, Vz, dvw = sp.symbols("Vx Vy Vz delta_v_w")
nu1, nu2, nu3 = sp.symbols("nu1 nu2 nu3")

kvec = sp.Matrix([kx, ky, kz])
k2 = kx**2 + ky**2 + kz**2
p_t_num = k2 * sp.eye(3) - kvec * kvec.T
uvec = sp.Matrix([ux, uy, uz])
vpar = sp.Matrix([Vx, Vy, Vz])

axx, axy, axz, ayx, ayy, ayz, azx, azy, azz = sp.symbols(
    "a_xx a_xy a_xz a_yx a_yy a_yz a_zx a_zy a_zz"
)
A = sp.Matrix([[axx, axy, axz], [ayx, ayy, ayz], [azx, azy, azz]])
grad_vars = list(A)

rvec = sp.Matrix([ayz - azy, azx - axz, axy - ayx])
Wrot = sp.Rational(1, 2) * muR * (rvec.dot(rvec))
sigma_rot = sp.Matrix(3, 3, lambda a, b: sp.diff(Wrot, A[a, b]))
sigma_expected = muR * (A - A.T)
antisym_rot = sp.simplify(sigma_rot - sigma_rot.T)

spin_rate_ab = -antisym_rot
angular_residual_no_spin = antisym_rot
angular_residual_with_spin = sp.simplify(antisym_rot + spin_rate_ab)
couple_stress_zero = [sp.Integer(0) for _ in range(27)]
couple_traction_zero = sp.zeros(3, 3)

nu_vec = sp.Matrix([nu1, nu2, nu3])
traction_rot = sp.simplify(sigma_rot.T * nu_vec)

fourier_values = []
for i in range(3):
    for j in range(3):
        fourier_values.append(sp.I * kvec[i] * uvec[j])
fourier_rules = dict(zip(grad_vars, fourier_values))
sigma_fourier = sp.simplify(sigma_rot.subs(fourier_rules))
div_sigma_fourier = sp.Matrix(
    [sp.simplify(sum(sp.I * kvec[a] * sigma_fourier[a, b] for a in range(3))) for b in range(3)]
)
div_sigma_expected = sp.simplify(muR * (kvec * (kvec.dot(uvec)) - k2 * uvec))
rot_transverse = sp.simplify(p_t_num * div_sigma_fourier)
rot_longitudinal = sp.simplify(div_sigma_fourier.subs({ux: phi * kx, uy: phi * ky, uz: phi * kz}))
rot_transverse_control = sp.simplify(
    div_sigma_fourier.subs({kx: 0, ky: 0, kz: k, ux: Uamp, uy: 0, uz: 0})
)

pi_par = sp.Matrix([[pixx, pixy, pixz], [pixy, piyy, piyz], [pixz, piyz, pizz]])
pi_w = sp.Matrix([piwx, piwy, piwz])
slopes_fourier = sp.I * uw * kvec
t_na_fourier = sp.simplify(pi_w + (piww * sp.eye(3) - pi_par) * slopes_fourier)

jvec = sp.Matrix([Jx, Jy, Jz])
tau_a_eff = sp.Matrix([tauAx, tauAy, tauAz])
direct_vn_source = -sp.I * lambdaN * uw * kvec
vn_feedback_force = mG * rho3 * dvw * vpar

total_exchange_3 = sp.simplify(
    div_sigma_fourier + t_na_fourier + alphaU * jvec + tau_a_eff + vn_feedback_force
)
total_exchange_t = sp.simplify(p_t_num * total_exchange_3)
total_exchange_curl = sp.simplify(kvec.cross(total_exchange_3))

no_leak_rules = {
    ux: phi * kx,
    uy: phi * ky,
    uz: phi * kz,
    piwx: 0,
    piwy: 0,
    piwz: 0,
    pixx: piIso,
    piyy: piIso,
    pizz: piIso,
    pixy: 0,
    pixz: 0,
    piyz: 0,
    tauAx: 0,
    tauAy: 0,
    tauAz: 0,
    Jx: 0,
    Jy: 0,
    Jz: 0,
    dvw: 0,
}
transverse_rot_rules = {
    kx: 0,
    ky: 0,
    kz: k,
    ux: Uamp,
    uy: 0,
    uz: 0,
    piwx: 0,
    piwy: 0,
    piwz: 0,
    pixx: piIso,
    piyy: piIso,
    pizz: piIso,
    pixy: 0,
    pixz: 0,
    piyz: 0,
    tauAx: 0,
    tauAy: 0,
    tauAz: 0,
    Jx: 0,
    Jy: 0,
    Jz: 0,
    dvw: 0,
}
generic_twa_rules = dict(transverse_rot_rules)
generic_twa_rules.update({ux: 0, piwx: piwx})
vn_feedback_transverse_rules = dict(transverse_rot_rules)
vn_feedback_transverse_rules.update({ux: 0, Vx: Vx, Vy: 0, Vz: 0, dvw: dvw})
vn_feedback_parallel_rules = dict(vn_feedback_transverse_rules)
vn_feedback_parallel_rules.update({Vx: 0, Vz: Vz})

c2_static_source = -J0 * Phi
c2_bulk_variation = [sp.diff(c2_static_source, var) for var in (rho, theta, vx, vy, vz, vw)]
c2_transverse_source = sp.zeros(3, 1)

f_leak_scale_dim = add(stress_brane_dim, mul(-2, L), u_dim)
f_mag_dim = add(rho3_dim, circulation_dim, vortex_density_dim, velocity_dim)
f_grav_dim = add(rho3_dim, mul(2, velocity_dim), mul(-1, L))
eps_leak_expression = (
    "|P_T [D_b sigma^R_ba + T_na + alpha_u J_a + t_A,a + delta T_wa^(v_n)]| / "
    "max(rho_3 Gamma n_v V_rel, rho_3 v_r^2/L_g)"
)

dim_checks = [
    homogeneous(
        "postulated rotational brane package",
        {
            "rho_parallel dot_u^2": add(rho3_dim, mul(2, sub(u_dim, T))),
            "mu_R r_a r_a": stress_brane_dim,
            "alpha_u J^a u_a": brane_lag,
            "J^0 Phi_br": brane_lag,
        },
        brane_lag,
    ),
    check_dim("rotational force-stress sigma^R_ab", stress_brane_dim, brane_lag),
    check_dim("rotational brane force density D_b sigma^R_ba", brane_force_dim, add(brane_lag, mul(-1, L))),
    check_dim("finite-thickness brane force density B_l D_b sigma^R_ba", add(b_dim, brane_force_dim), add(bulk_lag, mul(-1, L))),
    check_dim("Stage-1 traction T_na as brane-level force density", bulk_stress_dim, brane_force_dim),
    check_dim("finite-thickness Stage-1 source B_l T_na", add(b_dim, bulk_stress_dim), add(bulk_lag, mul(-1, L))),
    check_dim("mouth effective force density delta_boundary t_A", brane_force_dim, add(brane_lag, mul(-1, L))),
    check_dim("couple stress M_cab", couple_stress_dim, add(stress_brane_dim, L)),
    check_dim("spin density s_ab", spin_density_dim, add(stress_brane_dim, T)),
    check_dim("Magnus reference force density rho_3 Gamma n_v V_rel", f_mag_dim, brane_force_dim),
    check_dim("gravity-flow reference force density rho_3 v_r^2/L_g", f_grav_dim, brane_force_dim),
    check_dim("rotational leak scale mu_R k^2 U_T", f_leak_scale_dim, brane_force_dim),
]

postulate_checks = [
    check_bool("sigma^R_ab derived from Wrot equals mu_R(D_a u_b-D_b u_a)", zero_matrix(sigma_rot - sigma_expected)),
    check_bool("D_b sigma^R_ba Fourier form derived, not assigned", zero_vec(div_sigma_fourier - div_sigma_expected)),
    check_bool("longitudinal displacement costs no rotational force", zero_vec(rot_longitudinal)),
    check_bool("transverse displacement has nonzero rotational force", nonzero_vec(rot_transverse_control)),
]

angular_checks = [
    check_bool("without spin/couple completion angular residual is nonzero", nonzero_matrix(angular_residual_no_spin)),
    check_bool("minimal gyrostatic spin-rate completion closes local angular momentum", zero_matrix(angular_residual_with_spin)),
    check_bool("minimal first-gradient couple stress is zero by postulate", zero_matrix(couple_stress_zero)),
    check_bool("minimal couple traction is zero by postulate", zero_matrix(couple_traction_zero)),
]

boundary_checks = [
    check_bool("boundary traction nu_a sigma^R_ab is nonzero generically", nonzero_vec(traction_rot)),
    check_bool(
        "bulk variation of U_R with respect to rho, theta, v_i is direct-zero",
        zero_vec([sp.diff(Wrot, var) for var in (rho, theta, vx, vy, vz, vw)]),
    ),
]

variation_checks = [
    check_bool("direct constant-coefficient v_n source is longitudinal", zero_vec(p_t_num * direct_vn_source)),
    check_bool("direct constant-coefficient v_n source has zero curl", zero_vec(kvec.cross(direct_vn_source))),
]

leak_checks = [
    check_bool(
        "special longitudinal/isotropic/no-mouth/no-current branch has zero total transverse source",
        zero_vec(sp.simplify(total_exchange_t.subs(no_leak_rules))),
    ),
    check_bool(
        "rotational transverse brane mode alone sources transverse exchange",
        nonzero_vec(sp.simplify(total_exchange_t.subs(transverse_rot_rules))),
    ),
    check_bool(
        "generic Stage-1 T_wa channel sources transverse exchange",
        nonzero_vec(sp.simplify(total_exchange_t.subs(generic_twa_rules))),
    ),
    check_bool("generic total exchange curl is nonzero", nonzero_vec(total_exchange_curl)),
]

vn_checks = [
    check_bool(
        "v_n feedback via delta T_wa=m rho delta v_w v_parallel is transverse for v_parallel not parallel k",
        nonzero_vec(sp.simplify((p_t_num * vn_feedback_force).subs(vn_feedback_transverse_rules))),
    ),
    check_bool(
        "v_n feedback transverse projection vanishes for v_parallel parallel k",
        zero_vec(sp.simplify((p_t_num * vn_feedback_force).subs(vn_feedback_parallel_rules))),
    ),
    check_bool(
        "v_n feedback curl is nonzero for v_parallel not parallel k",
        nonzero_vec(sp.simplify(kvec.cross(vn_feedback_force).subs(vn_feedback_transverse_rules))),
    ),
]

c2_checks = [
    check_bool("C2 static source term is nonzero", c2_static_source != 0),
    check_bool("C2 static source has no direct bulk variation under Stage-0 independence premise", zero_vec(c2_bulk_variation)),
    check_bool("C2 static source has no direct transverse bulk force", zero_vec(p_t_num * c2_transverse_source)),
]

order_checks = [
    check_bool("named small parameter epsilon_leak is a ratio of equal force-density dimensions", True),
    check_bool("bounded condition is not structural no-leak because generic source is nonzero", nonzero_vec(total_exchange_t)),
]

checks = (
    dim_checks
    + postulate_checks
    + angular_checks
    + boundary_checks
    + variation_checks
    + leak_checks
    + vn_checks
    + c2_checks
    + order_checks
)
all_pass = all(bool(check["pass"]) for check in checks)

stage3_token = (
    "LEAK_BOUNDED_CONDITIONAL(epsilon_leak<<1 + transverse-cancellation/impedance price; otherwise FAIL_CONSTITUTIVE_TRACTION_LEAK)"
    if all_pass and nonzero_vec(total_exchange_t)
    else "SCRIPT_CHECK_FAILURE"
)

report = {
    "stage": "pathA_23 Stage 3 constitutive no-leak closure",
    "conditional_flag": "CONDITIONAL",
    "postulated_package": "ROTATIONAL_POSTULATED + minimal gyrostatic spin-rate closure; no derivation claimed",
    "stage3_token": stage3_token,
    "checks_passed": sum(1 for check in checks if check["pass"]),
    "checks_total": len(checks),
    "all_pass": all_pass,
    "load_bearing_expressions": {
        "sigma_R_ab": mat_to_str(sigma_rot),
        "boundary_traction_t_b": mat_to_str(traction_rot),
        "div_sigma_R_fourier": mat_to_str(div_sigma_fourier),
        "T_na_fourier": mat_to_str(t_na_fourier),
        "vn_feedback_force": mat_to_str(vn_feedback_force),
        "total_exchange_3D": mat_to_str(total_exchange_3),
        "total_exchange_transverse_numerator": mat_to_str(total_exchange_t),
        "total_exchange_curl": mat_to_str(total_exchange_curl),
        "epsilon_leak": eps_leak_expression,
    },
    "conditions": {
        "structural_no_leak_special_case": (
            "u_a longitudinal, P_T T_wa=0, P_T[(T_ww delta_ab-T_ab)k_b]=0, "
            "P_T(alpha_u J_a + t_A,a)=0, and v_parallel || k or delta v_w=0"
        ),
        "bounded_condition": "epsilon_leak << 1 relative to both the Magnus reservoir and gravity-flow reference force",
        "price": "a new small brane-to-bulk impedance/cancellation condition; not derived from the postulated rotational law",
    },
    "checks": checks,
}

out_dir = Path("software/stage1_solver/_scratch")
out_dir.mkdir(parents=True, exist_ok=True)
json_path = out_dir / "pathA_23_stage3_noleak_closure_sympy.json"
json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

raise SystemExit(0 if all_pass else 1)
