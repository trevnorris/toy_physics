#!/usr/bin/env python3
"""Machine checks for pathA_24 T0 frozen polar-OP arithmetic.

SymPy is the primary engine here.  This script checks:
  1. dimensions of the three frozen L_pol terms with units restored;
  2. the quadratic OP mode spectrum about a uniform |P0|=1 vacuum.

The frozen action block itself is not edited.  The literal L_pol lines below
are used as a fidelity guard before the algebra is run.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import sympy as sp


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
SCRATCH = STAGE1_ROOT / "_scratch"
FREEZE_REPORT = STAGE1_ROOT / "reports" / "pathA_24_T0_freeze.md"

EXPECTED_L_POL_LINES = [
    "  L_pol =",
    "      (1/2) m rho a^2 (D_t^v P^i)(D_t^v P^i)",
    "    - (1/2) m rho c_s^2(rho) a^2 (partial_j P^i)(partial_j P^i)",
    "    - (1/4) m rho c_s^2(rho) (P^i P^i - 1)^2.",
]


def extract_freeze_block(path: Path) -> list[str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    in_block = False
    block: list[str] = []
    for line in lines:
        if line == "```freeze-action":
            in_block = True
            continue
        if in_block and line == "```":
            return block
        if in_block:
            block.append(line)
    raise AssertionError(f"freeze-action block not found in {path}")


def assert_l_pol_fidelity() -> None:
    block = extract_freeze_block(FREEZE_REPORT)
    for expected in EXPECTED_L_POL_LINES:
        if expected not in block:
            raise AssertionError(f"frozen L_pol line missing or changed: {expected!r}")
    indices = [block.index(line) for line in EXPECTED_L_POL_LINES]
    if indices != sorted(indices):
        raise AssertionError(f"frozen L_pol lines are out of order: {indices}")


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


def assert_zero(name: str, expr: sp.Expr) -> None:
    simplified = sp.simplify(expr)
    if simplified != 0:
        raise AssertionError(f"{name}: expected 0, got {sp.sstr(simplified)}")


def check_dimensions() -> dict:
    # Dimension exponent triples are ordered (M, L, T).
    M = sp.Matrix([1, 0, 0])
    L = sp.Matrix([0, 1, 0])
    T = sp.Matrix([0, 0, 1])
    dimless = sp.Matrix([0, 0, 0])

    d_m = M
    d_rho = -4 * L
    d_K = M + 18 * L - 2 * T
    d_a = L
    d_DtP = -T
    d_gradP = -L
    d_P = dimless

    # c_s^2(rho) := 5 K rho^4 / m; the numeric 5 is dimensionless.
    d_cs2_from_eos = d_K + 4 * d_rho - d_m
    d_cs2_expected = 2 * L - 2 * T
    assert_dim("c_s^2 = 5 K rho^4 / m", d_cs2_from_eos, d_cs2_expected)

    # Faithful term structure from the freeze:
    #   (1/2) m rho a^2 (D_t^v P^i)(D_t^v P^i)
    # - (1/2) m rho c_s^2(rho) a^2 (partial_j P^i)(partial_j P^i)
    # - (1/4) m rho c_s^2(rho) (P^i P^i - 1)^2.
    kinetic = d_m + d_rho + 2 * d_a + 2 * d_DtP
    gradient = d_m + d_rho + d_cs2_from_eos + 2 * d_a + 2 * d_gradP
    potential = d_m + d_rho + d_cs2_from_eos + 4 * d_P

    lag_density = M - 2 * L - 2 * T
    action = M + 2 * L - T
    measure_dt_d4x = T + 4 * L

    term_dims = {
        "kinetic": kinetic,
        "gradient": gradient,
        "potential": potential,
    }
    for name, actual in term_dims.items():
        assert_dim(f"L_pol {name} term", actual, lag_density)

    unique_dims = {dim_tuple(dim) for dim in term_dims.values()}
    if len(unique_dims) != 1:
        raise AssertionError(f"L_pol terms are not homogeneous: {unique_dims}")

    assert_dim("int dt d^4X L_pol", lag_density + measure_dt_d4x, action)

    return {
        "cs2_from_eos": {
            "triple_MLT": dim_tuple(d_cs2_from_eos),
            "string": dim_str(d_cs2_from_eos),
        },
        "lagrangian_density_expected": {
            "triple_MLT": dim_tuple(lag_density),
            "string": dim_str(lag_density),
        },
        "term_dimensions": {
            name: {"triple_MLT": dim_tuple(dim), "string": dim_str(dim)}
            for name, dim in term_dims.items()
        },
        "action_dimension": {
            "triple_MLT": dim_tuple(lag_density + measure_dt_d4x),
            "string": dim_str(lag_density + measure_dt_d4x),
        },
    }


def coefficient(expr: sp.Expr, sym: sp.Symbol) -> sp.Expr:
    return sp.expand(expr).coeff(sym, 2)


def check_modes() -> dict:
    eps = sp.symbols("eps")
    m, rho0, K, a = sp.symbols("m rho0 K a", positive=True)
    omega, kx, ky, kz, kw = sp.symbols("omega kx ky kz kw", real=True)

    sigma, pi1, pi2, pi3 = sp.symbols("sigma pi1 pi2 pi3", real=True)
    sigma_t, pi1_t, pi2_t, pi3_t = sp.symbols(
        "sigma_t pi1_t pi2_t pi3_t", real=True
    )
    q = sp.Matrix([sigma, pi1, pi2, pi3])
    q_t = sp.Matrix([sigma_t, pi1_t, pi2_t, pi3_t])

    grad = sp.Matrix(
        4,
        4,
        lambda i, j: sp.symbols(f"g{i + 1}{j + 1}", real=True),
    )

    # Uniform background used for T0 arithmetic verification:
    # rho=rho0 const, v=0, P0=(1,0,0,0).  This is an O(4)-rotated
    # representative of any |P0|=1 vacuum.
    cs20 = sp.Rational(5) * K * rho0**4 / m
    I_P = m * rho0 * a**2
    K_P = m * rho0 * cs20 * a**2
    B = m * rho0 * cs20

    P = sp.Matrix([1 + eps * sigma, eps * pi1, eps * pi2, eps * pi3])
    DtvP = eps * q_t
    gradP = eps * grad

    # Transcription of the frozen density after the stated background
    # substitutions v=0, rho=rho0:
    L_pol = (
        sp.Rational(1, 2) * m * rho0 * a**2 * (DtvP.dot(DtvP))
        - sp.Rational(1, 2)
        * m
        * rho0
        * cs20
        * a**2
        * sum(gradP[i, j] ** 2 for i in range(4) for j in range(4))
        - sp.Rational(1, 4) * m * rho0 * cs20 * (P.dot(P) - 1) ** 2
    )
    L2 = sp.expand(sp.series(L_pol, eps, 0, 3).removeO().coeff(eps, 2))

    kinetic_matrix = sp.simplify(
        sp.Matrix([[sp.diff(L2, q_t[i], q_t[j]) for j in range(4)] for i in range(4)])
    )
    stiffness_matrix_by_dir = [
        sp.simplify(
            -sp.Matrix(
                [[sp.diff(L2, grad[i, direction], grad[j, direction]) for j in range(4)] for i in range(4)]
            )
        )
        for direction in range(4)
    ]
    mass_matrix = sp.simplify(
        -sp.Matrix([[sp.diff(L2, q[i], q[j]) for j in range(4)] for i in range(4)])
    )

    expected_kinetic = sp.eye(4) * I_P
    expected_stiffness = sp.eye(4) * K_P
    expected_mass = sp.diag(2 * B, 0, 0, 0)

    if sp.simplify(kinetic_matrix - expected_kinetic) != sp.zeros(4):
        raise AssertionError(f"unexpected kinetic matrix: {kinetic_matrix}")
    for direction, stiffness_matrix in enumerate(stiffness_matrix_by_dir, start=1):
        if sp.simplify(stiffness_matrix - expected_stiffness) != sp.zeros(4):
            raise AssertionError(
                f"unexpected stiffness matrix in direction {direction}: {stiffness_matrix}"
            )
    if sp.simplify(mass_matrix - expected_mass) != sp.zeros(4):
        raise AssertionError(f"unexpected potential Hessian/mass matrix: {mass_matrix}")

    k2 = kx**2 + ky**2 + kz**2 + kw**2
    stiffness = sp.Matrix([stiffness_matrix_by_dir[0][i, i] for i in range(4)])
    inertia = sp.Matrix([kinetic_matrix[i, i] for i in range(4)])
    mass_diag = sp.Matrix([mass_matrix[i, i] for i in range(4)])
    wave_operators = [
        sp.factor(-inertia[i] * omega**2 + stiffness[i] * k2 + mass_diag[i])
        for i in range(4)
    ]
    omega2 = [
        sp.factor(sp.simplify((stiffness[i] * k2 + mass_diag[i]) / inertia[i]))
        for i in range(4)
    ]

    transverse_indices = [i for i, value in enumerate(mass_diag) if sp.simplify(value) == 0]
    longitudinal_indices = [
        i for i, value in enumerate(mass_diag) if sp.simplify(value) != 0
    ]
    if longitudinal_indices != [0]:
        raise AssertionError(f"expected one longitudinal amplitude mode, got {longitudinal_indices}")
    if transverse_indices != [1, 2, 3]:
        raise AssertionError(f"expected three transverse modes, got {transverse_indices}")

    c_op_sq = sp.factor(sp.simplify(stiffness[1] / inertia[1]))
    c_s0_sq = cs20
    assert_zero("c_OP^2 - c_s0^2", c_op_sq - c_s0_sq)
    assert_zero("K_P/I_P - c_s0^2", K_P / I_P - c_s0_sq)

    gap_sq = sp.factor(sp.simplify(mass_diag[0] / inertia[0]))
    if not gap_sq.is_positive:
        raise AssertionError(f"longitudinal gap^2 is not positive under assumptions: {gap_sq}")
    if any(sp.simplify(omega2[i] - c_s0_sq * k2) != 0 for i in transverse_indices):
        raise AssertionError(f"unexpected transverse dispersion: {omega2}")

    return {
        "quadratic_lagrangian": sp.sstr(L2),
        "I_P": sp.sstr(I_P),
        "K_P": sp.sstr(sp.factor(K_P)),
        "c_s0_squared": sp.sstr(c_s0_sq),
        "c_OP_squared": sp.sstr(c_op_sq),
        "c_OP": f"sqrt({sp.sstr(c_op_sq)})",
        "transverse_count": len(transverse_indices),
        "longitudinal_count": len(longitudinal_indices),
        "longitudinal_gap_squared": sp.sstr(gap_sq),
        "longitudinal_gap": f"sqrt({sp.sstr(gap_sq)})",
        "wave_operators": [sp.sstr(op) for op in wave_operators],
        "dispersion_omega_squared": [sp.sstr(expr) for expr in omega2],
        "mass_matrix_diag": [sp.sstr(sp.factor(x)) for x in mass_diag],
    }


def print_term_dimensions(term_dimensions: dict[str, dict]) -> None:
    parts = [
        f"{name}={tuple(value['triple_MLT'])}"
        for name, value in term_dimensions.items()
    ]
    print("term_dimensions_MLT: " + " ".join(parts))


def main() -> int:
    assert_l_pol_fidelity()
    dimensions = check_dimensions()
    modes = check_modes()

    report = {
        "schema": "stage1_pathA_24_T0_dimcheck_sympy/v1",
        "engine": "sympy",
        "pass": True,
        "freeze_fidelity_guard": {
            "path": str(FREEZE_REPORT),
            "checked_lines": EXPECTED_L_POL_LINES,
        },
        "dimensions": dimensions,
        "modes": modes,
    }
    SCRATCH.mkdir(parents=True, exist_ok=True)
    out_path = SCRATCH / "pathA_24_T0_dimcheck_sympy.json"
    out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("T0_DIMCHECK_SYMPY: PASS")
    print_term_dimensions(dimensions["term_dimensions"])
    print(f"action_dimension_MLT: {tuple(dimensions['action_dimension']['triple_MLT'])}")
    print(f"c_OP_sympy: {modes['c_OP']}")
    print(f"longitudinal_gap_sympy: {modes['longitudinal_gap']}")
    print(f"transverse_count_sympy: {modes['transverse_count']}")
    print(f"wrote: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
