#!/usr/bin/env python3
"""Independent SymPy cross-check for pathA_23 Stage 1.

This redo derives the source structures from the Stage-0 coupling and the
scalar-sector stress projection.  It deliberately does not assume Pi_n is a
free scalar.  Scope is pre-constitutive and linear in the brane slope unless a
term is explicitly labelled deferred/O(slope^2).
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


def _items(expr: sp.Matrix | sp.Expr) -> list[sp.Expr]:
    if isinstance(expr, sp.MatrixBase):
        return list(expr)
    return [expr]


def zero_vec(vec: sp.Matrix | sp.Expr) -> bool:
    return all(sp.simplify(component) == 0 for component in _items(vec))


def nonzero_vec(vec: sp.Matrix | sp.Expr) -> bool:
    return not zero_vec(vec)


def free_symbols(expr: sp.Matrix | sp.Expr) -> set[sp.Symbol]:
    out: set[sp.Symbol] = set()
    for item in _items(expr):
        out.update(item.free_symbols)
    return out


def free_of_all(expr: sp.Matrix | sp.Expr, symbols: Iterable[sp.Symbol]) -> bool:
    return free_symbols(expr).isdisjoint(set(symbols))


def lin_slope_expr(expr: sp.Expr) -> sp.Expr:
    scaled = expr.subs({sx: eps * sx, sy: eps * sy, sz: eps * sz})
    return sp.expand(sp.series(scaled, eps, 0, 2).removeO().subs(eps, 1))


def lin_slope(obj: sp.Matrix | sp.Expr) -> sp.Matrix | sp.Expr:
    if isinstance(obj, sp.MatrixBase):
        return obj.applyfunc(lin_slope_expr)
    return lin_slope_expr(obj)


def lin_brane_expr(expr: sp.Expr) -> sp.Expr:
    scaled = expr.subs({uw: eps * uw, sx: eps * sx, sy: eps * sy, sz: eps * sz})
    return sp.expand(sp.series(scaled, eps, 0, 2).removeO().subs(eps, 1))


def lin_brane(obj: sp.Matrix | sp.Expr) -> sp.Matrix | sp.Expr:
    if isinstance(obj, sp.MatrixBase):
        return obj.applyfunc(lin_brane_expr)
    return lin_brane_expr(obj)


def curl3(vec: sp.Matrix) -> sp.Matrix:
    return kvec.cross(vec)


def mat_to_str(obj: sp.Matrix | sp.Expr | list[sp.Expr]) -> str:
    return str(obj)


dim0: Dim = (0, 0, 0, 0)
L: Dim = (1, 0, 0, 0)
T: Dim = (0, 1, 0, 0)
M: Dim = (0, 0, 1, 0)
Q: Dim = (0, 0, 0, 1)

energy = add(M, mul(2, L), mul(-2, T))
bulk_lag = add(energy, mul(-4, L))
brane_lag = add(energy, mul(-3, L))
b_dim = mul(-1, L)

u_dim = L
uw_dim = L
stress_dim = add(energy, mul(-4, L))
velocity_dim = sub(L, T)
lambda_vn_dim = sub(brane_lag, velocity_dim)
a_spatial_dim = add(M, L, mul(-1, T), mul(-1, Q))
phi_dim = add(M, mul(2, L), mul(-2, T), mul(-1, Q))
alpha_u_dim = sub(a_spatial_dim, u_dim)
j0_dim = add(Q, mul(-3, L))
ja_dim = add(Q, mul(-2, L), mul(-1, T))

kx, ky, kz = sp.symbols("kx ky kz")
sx, sy, sz, eps = sp.symbols("sx sy sz eps")
uw, Bl, rho, theta, mG, Pbulk, dpinn_drho, lambda_n = sp.symbols(
    "uw B_l rho theta m_GNLS P_bulk dPi_nn_drho Lambda_n"
)
J0, Phi, alpha_u = sp.symbols("J0 Phi alpha_u")
ux, uy, uz = sp.symbols("ux uy uz")
Jx, Jy, Jz = sp.symbols("Jx Jy Jz")
tx, ty, tz = sp.symbols("tx ty tz")
lambda_br, mu_br, mu_r = sp.symbols("lambda_br mu_br mu_R")
lambda_c, mu_c, kappa_c, A_c, varpi_gap = sp.symbols(
    "lambda_c mu_c kappa_c A_c varpi_gap"
)
pixx, pixy, pixz, piyy, piyz, pizz, piwx, piwy, piwz, piww = sp.symbols(
    "Pi_xx Pi_xy Pi_xz Pi_yy Pi_yz Pi_zz Pi_wx Pi_wy Pi_wz Pi_ww"
)
pi_iso, pi_aniso = sp.symbols("Pi_iso Pi_aniso")
vx, vy, vz, vw, sig_nn = sp.symbols("vx vy vz vw sigma_nn")

kvec = sp.Matrix([kx, ky, kz])
k2 = kx**2 + ky**2 + kz**2
p_t_num = k2 * sp.eye(3) - kvec * kvec.T
slopes = sp.Matrix([sx, sy, sz])
fourier_slopes = sp.I * uw * kvec

# General scalar-sector matter stress representative:
# Pi_ij = m_GNLS*rho*v_i*v_j + delta_ij*P(rho) + sigma_Q,ij.
# The projection is first done for general components; the velocity variation
# then checks the convective part explicitly.
Pi4 = sp.Matrix(
    [
        [pixx, pixy, pixz, piwx],
        [pixy, piyy, piyz, piwy],
        [pixz, piyz, pizz, piwz],
        [piwx, piwy, piwz, piww],
    ]
)
pi_par = Pi4[:3, :3]
pi_w = sp.Matrix([piwx, piwy, piwz])
n4 = sp.Matrix([-sx, -sy, -sz, 1])
e_tangents = [
    sp.Matrix([1, 0, 0, sx]),
    sp.Matrix([0, 1, 0, sy]),
    sp.Matrix([0, 0, 1, sz]),
]

pi_nn_exact = (n4.T * Pi4 * n4)[0]
pi_na_exact = sp.Matrix([(n4.T * Pi4 * e_tangent)[0] for e_tangent in e_tangents])
pi_nn_lin = lin_slope(pi_nn_exact)
pi_na_lin = lin_slope(pi_na_exact)
expected_pi_nn_lin = piww - 2 * (slopes.dot(pi_w))
expected_pi_na_lin = pi_w + (piww * sp.eye(3) - pi_par) * slopes

t_na_fourier = sp.simplify(pi_na_lin.subs({sx: fourier_slopes[0], sy: fourier_slopes[1], sz: fourier_slopes[2]}))
t_na_transverse = sp.simplify(p_t_num * t_na_fourier)
t_na_curl = sp.simplify(curl3(t_na_fourier))

iso_rules = {
    piwx: 0,
    piwy: 0,
    piwz: 0,
    pixx: pi_iso,
    piyy: pi_iso,
    pizz: pi_iso,
    pixy: 0,
    pixz: 0,
    piyz: 0,
}
t_na_iso = sp.simplify(t_na_fourier.subs(iso_rules))

# Variation of B_l u_w Pi_nn. A scalar variation of Pi_nn gives a scalar
# source, so its Euler-force representation is a Fourier gradient.
density_scalar_variation = Bl * uw * dpinn_drho
density_euler_source = sp.I * density_scalar_variation * kvec
normal_work_scalar_longitudinal = sp.I * Bl * piww * uw * kvec

v4 = sp.Matrix([vx, vy, vz, vw])
v_n = (n4.T * v4)[0]
pi_nn_hydro_velocity = mG * rho * v_n**2 + Pbulk + sig_nn
d_pi_nn_dv = sp.Matrix([lin_slope(sp.diff(pi_nn_hydro_velocity, v)) for v in v4])
d_l_normal_dv_a = Bl * uw * d_pi_nn_dv[:3, 0]
d_l_normal_dv_a_linear_brane = lin_brane(d_l_normal_dv_a)

vn_scalar_coupling = Bl * lambda_n * v_n
d_vn_coupling_dv_a = sp.Matrix([lin_slope(sp.diff(vn_scalar_coupling, v)) for v in (vx, vy, vz)])
vn_velocity_in_plane_source = sp.simplify(
    d_vn_coupling_dv_a.subs({sx: fourier_slopes[0], sy: fourier_slopes[1], sz: fourier_slopes[2]})
)

current_coupling = Bl * (alpha_u * (Jx * ux + Jy * uy + Jz * uz) - J0 * Phi)
current_bulk_variations = [
    sp.diff(current_coupling, var) for var in (rho, theta, vx, vy, vz, vw)
]
current_brane_variation = sp.Matrix([sp.diff(current_coupling, var) for var in (ux, uy, uz)])
static_coulomb_coupling = -Bl * J0 * Phi
static_coulomb_bulk_source = sp.Matrix([0, 0, 0])

source_set = {
    "density scalar variation from B_l u_w Pi_nn": density_euler_source,
    "normal-work scalar force from B_l u_w Pi_nn": normal_work_scalar_longitudinal,
    "optional v_n scalar tangential velocity variation": vn_velocity_in_plane_source,
    "direct source-current bulk variation": static_coulomb_bulk_source,
}

candidate_coeffs = {
    "Cauchy/Navier": (lambda_br, mu_br),
    "rotational/MacCullagh": (mu_r,),
    "Cosserat/micropolar": (lambda_c, mu_c, kappa_c, A_c, varpi_gap),
}

dim_checks = [
    check_dim("stress projection Pi_nn or Pi_na", stress_dim, add(energy, mul(-4, L))),
    check_dim("normal work brane density u_w Pi_nn", add(uw_dim, stress_dim), brane_lag),
    check_dim(
        "finite-thickness normal work bulk density B_l u_w Pi_nn",
        add(b_dim, uw_dim, stress_dim),
        bulk_lag,
    ),
    check_dim("finite-thickness traction B_l Pi_na", add(b_dim, stress_dim), add(energy, mul(-5, L))),
    check_dim("velocity-intersection coefficient Lambda_n", add(lambda_vn_dim, velocity_dim), brane_lag),
    check_dim("bulk representation of Lambda_n v_n", add(b_dim, lambda_vn_dim, velocity_dim), bulk_lag),
    check_dim("A_a^br = alpha_u u_a", add(alpha_u_dim, u_dim), a_spatial_dim),
    homogeneous(
        "source-current brane density",
        {
            "J^a A_a^br": add(ja_dim, a_spatial_dim),
            "J^0 Phi_br": add(j0_dim, phi_dim),
        },
        brane_lag,
    ),
]

projection_derivation_checks = [
    check_bool(
        "derived Pi_nn projection matches T_ww - 2 s_b T_wb at O(slope)",
        sp.simplify(pi_nn_lin - expected_pi_nn_lin) == 0,
    ),
    check_bool(
        "derived Pi_na projection matches T_wa + (T_ww delta_ab - T_ab) s_b at O(slope)",
        zero_vec(pi_na_lin - expected_pi_na_lin),
    ),
    check_bool(
        "derived Pi_na transverse projector is generically nonzero",
        nonzero_vec(t_na_transverse),
        True,
        "This is why Stage 1 cannot claim no-leak without profile/traction assumptions.",
    ),
    check_bool(
        "derived Pi_na in-plane curl is generically nonzero",
        nonzero_vec(t_na_curl),
    ),
    check_bool(
        "isotropic in-plane stress with zero normal-tangent stress gives longitudinal Pi_na",
        zero_vec(p_t_num * t_na_iso),
    ),
    check_bool(
        "isotropic in-plane stress with zero normal-tangent stress gives zero curl",
        zero_vec(curl3(t_na_iso)),
    ),
]

variation_checks: list[dict] = []
for name, source in source_set.items():
    variation_checks.append(
        check_bool(
            f"{name} has zero transverse projector",
            zero_vec(p_t_num * source),
            True,
            "This applies only to scalar/normal pieces, not to the derived Pi_na traction channel.",
        )
    )
    variation_checks.append(
        check_bool(f"{name} has zero in-plane curl", zero_vec(curl3(source)))
    )

variation_checks.extend(
    [
        check_bool(
            "B_l u_w Pi_nn convective dL/dv_a is O(u_w slope), hence absent at linear brane order",
            zero_vec(d_l_normal_dv_a_linear_brane),
        ),
        check_bool(
            "source-current variation with respect to rho, theta, and v_i is zero under Stage-1 independence premise",
            all(sp.simplify(item) == 0 for item in current_bulk_variations),
        ),
        check_bool(
            "source-current brane variation is nonzero and proportional to J^a",
            zero_vec(current_brane_variation - Bl * alpha_u * sp.Matrix([Jx, Jy, Jz])),
        ),
    ]
)

c2_checks = [
    check_bool(
        "static charge coupling is not zeroed",
        static_coulomb_coupling != 0,
        True,
        "The C2 pass uses -B_l J^0 Phi_br as a nonzero static/longitudinal source.",
    ),
    check_bool(
        "static charge coupling contains J0",
        J0 in free_symbols(static_coulomb_coupling),
    ),
    check_bool(
        "static Coulomb direct bulk variation has zero transverse source",
        zero_vec(p_t_num * static_coulomb_bulk_source),
    ),
    check_bool(
        "static Coulomb direct bulk variation has zero in-plane curl",
        zero_vec(curl3(static_coulomb_bulk_source)),
    ),
]

candidate_checks = [
    check_bool(
        f"derived coupling/projection sources independent of {candidate} constitutive coefficients",
        all(free_of_all(source, coeffs) for source in list(source_set.values()) + [t_na_fourier]),
        True,
        "Stage 1 uses S_cpl plus scalar bulk stress only; constitutive brane traction/couple-stress remains later-stage.",
    )
    for candidate, coeffs in candidate_coeffs.items()
]

negative_checks = [
    check_bool(
        "negative control: arbitrary tangential source is detected by transverse projector",
        nonzero_vec(p_t_num * sp.Matrix([tx, ty, tz])),
    ),
    check_bool(
        "negative control: arbitrary tangential source is detected by curl",
        nonzero_vec(curl3(sp.Matrix([tx, ty, tz]))),
    ),
    check_bool(
        "negative control: forbidden direct J^a-to-bulk-velocity coupling is detected by transverse projector",
        nonzero_vec(p_t_num * (alpha_u * sp.Matrix([Jx, Jy, Jz]))),
    ),
    check_bool(
        "negative control: forbidden direct J^a-to-bulk-velocity coupling is detected by curl",
        nonzero_vec(curl3(alpha_u * sp.Matrix([Jx, Jy, Jz]))),
    ),
]

stress_channel_checks = [
    check_bool(
        "reachable FAIL condition: known generic T_wa would trigger transverse-source failure",
        nonzero_vec(
            t_na_transverse.subs(
                {
                    piwx: tx,
                    piwy: ty,
                    piwz: tz,
                    pixx: 0,
                    piyy: 0,
                    pizz: 0,
                    pixy: 0,
                    pixz: 0,
                    piyz: 0,
                }
            )
        ),
    ),
    check_bool(
        "reachable FAIL condition: anisotropic in-plane stress would trigger transverse-source failure",
        nonzero_vec(
            t_na_transverse.subs(
                {
                    piwx: 0,
                    piwy: 0,
                    piwz: 0,
                    pixx: pi_aniso,
                    piyy: 0,
                    pizz: 0,
                    pixy: 0,
                    pixz: 0,
                    piyz: 0,
                }
            )
        ),
    ),
]

bending_checks = [
    check_bool(
        "u_w slope is retained with Fourier factor I",
        any(component.has(uw) for component in list(fourier_slopes))
        and any(component.has(sp.I) for component in list(fourier_slopes)),
    ),
    check_bool(
        "normal scalar bending source is longitudinal when Pi_nn is scalar/isotropic",
        zero_vec(p_t_num * normal_work_scalar_longitudinal),
    ),
    check_bool(
        "optional v_n scalar source is longitudinal at constant Lambda_n",
        zero_vec(p_t_num * vn_velocity_in_plane_source),
    ),
    check_bool(
        "optional v_n scalar source has zero curl at constant Lambda_n",
        zero_vec(curl3(vn_velocity_in_plane_source)),
    ),
]

checks = (
    dim_checks
    + projection_derivation_checks
    + variation_checks
    + c2_checks
    + candidate_checks
    + negative_checks
    + stress_channel_checks
    + bending_checks
)

all_pass = all(item["pass"] for item in checks)
c2_pass = (
    static_coulomb_coupling != 0
    and zero_vec(p_t_num * static_coulomb_bulk_source)
    and zero_vec(curl3(static_coulomb_bulk_source))
)
scalar_sources_longitudinal = all(
    zero_vec(p_t_num * source) and zero_vec(curl3(source)) for source in source_set.values()
)
stress_projection_proven_no_leak = zero_vec(t_na_transverse) and zero_vec(t_na_curl)

mouth_traction_deferred = True
stress_projection_deferred = not stress_projection_proven_no_leak
vn_feedback_deferred = True
bending_o2_deferred = True
known_stress_leak = False
known_mouth_leak = False

if not all_pass:
    stage1_token = "SCRIPT_CHECK_FAILED"
elif known_stress_leak or known_mouth_leak:
    stage1_token = "FAIL_COUPLING_SOURCES_TRANSVERSE"
elif not c2_pass:
    stage1_token = "FAIL_COUPLING_SOURCES_TRANSVERSE"
elif not scalar_sources_longitudinal:
    stage1_token = "FAIL_COUPLING_LEAKS_INDEPENDENT_OF_CONSTITUTIVE"
elif (
    stress_projection_deferred
    or mouth_traction_deferred
    or vn_feedback_deferred
    or bending_o2_deferred
):
    stage1_token = "LEAK_CONDITIONS_DEFERRED"
else:
    stage1_token = "NO_KINEMATIC_LEAK_FOR_ALL_CANDIDATES"

report = {
    "schema": "pathA_23_stage1_kinematic_leak_audit_sympy/v2",
    "scope": "Stage-1 pre-constitutive coupling leak audit; stress projection derived; linear in brane slope",
    "pass": all_pass,
    "stage1_token": stage1_token,
    "derived_sources": {
        "Pi_nn_linear": mat_to_str(pi_nn_lin),
        "Pi_na_linear": mat_to_str(pi_na_lin),
        "Pi_na_fourier": mat_to_str(t_na_fourier),
        "Pi_na_transverse_projector": mat_to_str(t_na_transverse),
        "Pi_na_curl": mat_to_str(t_na_curl),
        "density_scalar_euler_source": mat_to_str(density_euler_source),
        "v_n_tangential_velocity_source": mat_to_str(vn_velocity_in_plane_source),
        "source_current_bulk_variations": mat_to_str(current_bulk_variations),
        "source_current_brane_variation": mat_to_str(current_brane_variation),
        "static_coulomb_coupling": mat_to_str(static_coulomb_coupling),
    },
    "no_leak_conditions": {
        "stress_projection": (
            "P_T^{ab} Pi_wb = 0 and P_T^{ab} Pi_bc k_c = 0, equivalently "
            "curl(k, Pi_w + i u_w (Pi_ww I - Pi_parallel) k)=0"
        ),
        "mouth": "P_T^{ab} t_A,b = 0 and curl(k,t_A)=0 at each mouth",
        "v_n_feedback": "normal-flow perturbation must not generate bulk vorticity through the Euler/vorticity evolution",
    },
    "deferred_quantities": [
        "actual scalar stress components Pi_wa and anisotropic Pi_ab on the selected throat/brane branch",
        "mouth tangential traction t_A^a from S_mouth",
        "O(slope^2) bending/curvature corrections",
        "v_n -> bulk-vorticity feedback through the scalar Euler/vorticity equations",
    ],
    "outcomes": {
        "scalar_sources_longitudinal": scalar_sources_longitudinal,
        "stress_projection_proven_no_leak": stress_projection_proven_no_leak,
        "stress_projection_deferred": stress_projection_deferred,
        "mouth_traction_deferred": mouth_traction_deferred,
        "c2_sector_decomposition": c2_pass,
        "reachable_fail_token_for_known_transverse_stress_or_mouth": "FAIL_COUPLING_SOURCES_TRANSVERSE",
    },
    "checks": checks,
}

out_dir = Path("software/stage1_solver/_scratch")
out_dir.mkdir(parents=True, exist_ok=True)
json_path = out_dir / "pathA_23_stage1_kinematic_leak_audit_sympy.json"
json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

print(f"wrote {json_path}")
print(
    "pathA_23 Stage 1 SymPy leak audit: "
    f"{sum(1 for item in checks if item['pass'])}/{len(checks)} checks; token {stage1_token}"
)
raise SystemExit(0 if all_pass and stage1_token != "SCRIPT_CHECK_FAILED" else 1)
