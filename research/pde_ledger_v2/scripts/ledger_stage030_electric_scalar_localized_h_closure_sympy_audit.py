#!/usr/bin/env python3
"""Ledger stage030 SymPy audit: electric scalar + localized-H closure.

Standalone, print-only, no arguments, no file output.  The optional internal
ablation switch is the environment variable ``LEDGER_STAGE030_MUTATION``;
its accepted values are the predicate names in ``ABLATIONS`` below.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
from typing import Any

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0
MUTATION_ENV = "LEDGER_STAGE030_MUTATION"


class AuditFailure(AssertionError):
    """A named ledger predicate failed."""

    def __init__(self, predicate: str, detail: str = "") -> None:
        super().__init__(predicate)
        self.predicate = predicate
        self.detail = detail


@dataclass(frozen=True)
class Ablation:
    primitive: str
    mutated_value: Any
    description: str


ABLATIONS: dict[str, Ablation] = {
    "PASS_TRANSVERSE_FACTORIZATION": Ablation(
        "factor_A_coefficient", sp.Integer(3),
        "A superpotential coefficient: 2 -> 3",
    ),
    "PASS_PARENT_TRANSVERSE_ZERO_EIGENVALUE": Ablation(
        "f0_power", sp.Integer(1),
        "f0 sech exponent: 2 -> 1",
    ),
    "PASS_UNIQUE_NORMALIZABLE_KERNEL": Ablation(
        "kernel_derivative_coefficient", sp.Integer(2),
        "coefficient of d/dw in the kernel equation A f=0: 1 -> 2",
    ),
    "PASS_POSITIVE_CONTINUUM_GAP": Ablation(
        "V_asymptote_coefficient", sp.Integer(0),
        "V_H asymptote coefficient: 4 -> 0",
    ),
    "PASS_ZERO_MODE_NORM": Ablation(
        "norm_inner_weight", sp.Integer(3),
        "parent-H norm weight: 2 -> 3",
    ),
    "PASS_NORMALIZED_ZERO_MODE_PROJECTION": Ablation(
        "projection_inner_weight", sp.Integer(3),
        "projection weight in h=P0 H: 2 -> 3",
    ),
    "PASS_PHYSICAL_H_NORM": Ablation(
        "physical_M4_star", -sp.Rational(3, 160),
        "M4*: 3/160 -> -3/160",
    ),
    "PASS_REDUCED_H_INERTIA": Ablation(
        "inertia_M4_star", sp.Rational(3, 320),
        "M4*: 3/160 -> 3/320 in M_h=N0 M4",
    ),
    "PASS_PARENT_STIFFNESS_DIMENSIONAL_CONSISTENCY": Ablation(
        "parent_K4_dimension", (1, -2, 1),
        "[K4]: M L^2 T^-2 -> M L T^-2",
    ),
    "PASS_REDUCED_H_STIFFNESS": Ablation(
        "declared_Kh_star", sp.Integer(2),
        "declared K_h*: 1 -> 2",
    ),
    "PASS_REDUCED_SPEED_PRESERVATION": Ablation(
        "speed_Kh_N0_power", sp.Integer(2),
        "K_h reduction: N0 K4 -> N0^2 K4",
    ),
    "PASS_STABILITY": Ablation(
        "stability_Chu_star", sp.Integer(2),
        "C_hu*: 1/2 -> 2",
    ),
    "PASS_POSITIVE_GENERALIZED_WAVE_SPEEDS": Ablation(
        "wave_Aeff_star", sp.Integer(2),
        "A_eff*: 1 -> 2 while the Sylvester margin stays positive",
    ),
    "PASS_REDUCED_H_MASSLESSNESS": Ablation(
        "reduced_h_mass_star", sp.Integer(1),
        "reduced k^0 h^2 coefficient*: 0 -> 1",
    ),
    "PASS_CONSERVATIVE_HESSIAN_SYMMETRY": Ablation(
        "u_equation_cross_star", sp.Rational(3, 4),
        "u_L Euler cross coefficient*: C_hu*=1/2 -> 3/4",
    ),
    "PASS_DIMENSIONAL_HOMOGENEITY": Ablation(
        "firewall_Chu_dimension", (0, -1, 1),
        "[C_hu] in the mix term: M T^-2 -> M T^-1",
    ),
}


ACTIVE_MUTATION = os.environ.get(MUTATION_ENV, "").strip()


def banner(title: str) -> None:
    print("")
    print("=" * len(title))
    print(title)
    print("=" * len(title))


def subbanner(title: str) -> None:
    print("")
    print(title)
    print("-" * len(title))


def compact(value: Any) -> str:
    if isinstance(value, sp.MatrixBase):
        return sp.sstr(value)
    try:
        return sp.sstr(sp.factor(sp.cancel(sp.simplify(value))))
    except (TypeError, ValueError, AttributeError):
        return str(value)


def assert_no_float(name: str, value: Any) -> None:
    if isinstance(value, (tuple, list, dict, str)):
        return
    clean = sp.sympify(value)
    floats = clean.atoms(sp.Float)
    if floats:
        raise AuditFailure(name, f"machine Float atom(s): {floats}")


def expect_zero(name: str, residual: sp.Expr | int, evidence: Any = None) -> None:
    global PASS_COUNT, FAIL_COUNT
    assert_no_float(name, residual)
    clean = sp.simplify(residual)
    assert_no_float(name, clean)
    if clean == 0:
        PASS_COUNT += 1
        print(f"PASS  {name}")
        return
    FAIL_COUNT += 1
    print(f"FIRST_FAILURE={name}")
    print(f"FAIL  {name}: residual = {compact(clean)}")
    if evidence is not None:
        print(f"      evidence = {evidence}")
    raise AuditFailure(name, compact(clean))


def expect_bool(name: str, condition: bool, evidence: Any = None) -> None:
    expect_zero(name, sp.Integer(0) if condition else sp.Integer(1), evidence)


def primitive(predicate: str, name: str, card_value: Any) -> Any:
    """Return a target-local copy mutation, or the clean card primitive."""
    if ACTIVE_MUTATION == predicate:
        ablation = ABLATIONS[predicate]
        if ablation.primitive == name:
            return ablation.mutated_value
    return card_value


@dataclass(frozen=True)
class Dim:
    """Exact exponent vector in [L,T,M] order."""

    l: sp.Rational | int = 0
    t: sp.Rational | int = 0
    m: sp.Rational | int = 0

    def __post_init__(self) -> None:
        object.__setattr__(self, "l", sp.Rational(self.l))
        object.__setattr__(self, "t", sp.Rational(self.t))
        object.__setattr__(self, "m", sp.Rational(self.m))

    def __mul__(self, other: "Dim") -> "Dim":
        return Dim(self.l + other.l, self.t + other.t, self.m + other.m)

    def __truediv__(self, other: "Dim") -> "Dim":
        return Dim(self.l - other.l, self.t - other.t, self.m - other.m)

    def __pow__(self, power: int | sp.Rational) -> "Dim":
        p = sp.Rational(power)
        return Dim(p * self.l, p * self.t, p * self.m)

    def vector(self) -> tuple[sp.Rational, sp.Rational, sp.Rational]:
        return (self.l, self.t, self.m)

    def __str__(self) -> str:
        return f"[{self.l},{self.t},{self.m}]_[L,T,M]"


ONE = Dim()
LENGTH = Dim(1, 0, 0)
TIME = Dim(0, 1, 0)
MASS = Dim(0, 0, 1)
ENERGY = Dim(2, -2, 1)
VELOCITY = LENGTH / TIME

DIM_H = LENGTH**-1
DIM_h = ONE
DIM_u = LENGTH
DIM_N0 = LENGTH**-1
DIM_M4 = MASS
DIM_K4 = ENERGY
DIM_AEFF = Dim(-3, 0, 1)
DIM_BEFF = Dim(-1, -2, 1)
DIM_MH = Dim(-1, 0, 1)
DIM_KH = Dim(1, -2, 1)
DIM_CHU = Dim(0, -2, 1)
DIM_CE = VELOCITY
BULK_DENSITY = ENERGY / LENGTH**4
BRANE_DENSITY = ENERGY / LENGTH**3
DIM_D = Dim(0, -4, 2)


def dim_residual(actual: Dim, expected: Dim) -> sp.Expr:
    return sp.simplify(sum((a - e) ** 2 for a, e in zip(actual.vector(), expected.vector())))


def matrix_zero(matrix: sp.MatrixBase) -> bool:
    return all(sp.simplify(entry) == 0 for entry in matrix)


def positive(expr: sp.Expr) -> bool:
    result = sp.simplify(expr > 0)
    return result is sp.true or result is True


def run_spectral_structure() -> None:
    subbanner("Localized parent H: factorization, zero mode, unique kernel, gap")
    w, ell = sp.symbols("w ell", positive=True, real=True)
    probe = sp.Function("probe")(w)
    potential = (4 - 6 / sp.cosh(w / ell) ** 2) / ell**2

    factor_coefficient = primitive(
        "PASS_TRANSVERSE_FACTORIZATION", "factor_A_coefficient", sp.Integer(2)
    )
    factor_superpotential = factor_coefficient * sp.tanh(w / ell) / ell
    a_probe = sp.diff(probe, w) + factor_superpotential * probe
    adag_a_probe = -sp.diff(a_probe, w) + factor_superpotential * a_probe
    o_probe = -sp.diff(probe, w, 2) + potential * probe
    factor_residual = sp.simplify(sp.expand_trig(adag_a_probe - o_probe))
    expect_zero(
        "PASS_TRANSVERSE_FACTORIZATION",
        factor_residual,
        {"O_perp-A^dagger A": factor_residual, "nonnegative_form": "<A psi,A psi> >= 0"},
    )
    print("      O_perp=A^dagger A>=0 with A=d/dw+2 tanh(w/ell)/ell")

    f0_power = primitive(
        "PASS_PARENT_TRANSVERSE_ZERO_EIGENVALUE", "f0_power", sp.Integer(2)
    )
    f0_zero_test = 1 / (ell * sp.cosh(w / ell) ** f0_power)
    zero_residual = sp.simplify(
        -sp.diff(f0_zero_test, w, 2) + potential * f0_zero_test
    )
    expect_zero(
        "PASS_PARENT_TRANSVERSE_ZERO_EIGENVALUE",
        zero_residual,
        {"f0_power": f0_power, "O_perp f0": zero_residual},
    )

    kernel_derivative_coefficient = primitive(
        "PASS_UNIQUE_NORMALIZABLE_KERNEL",
        "kernel_derivative_coefficient",
        sp.Integer(1),
    )
    kernel_shape = sp.cosh(w / ell) ** (
        -sp.Rational(2) / kernel_derivative_coefficient
    )
    canonical_f0 = 1 / (ell * sp.cosh(w / ell) ** 2)
    kernel_equation_residual = sp.simplify(
        kernel_derivative_coefficient * sp.diff(kernel_shape, w)
        + 2 * sp.tanh(w / ell) * kernel_shape / ell
    )
    ratio_to_f0 = sp.simplify(kernel_shape / canonical_f0)
    kernel_norm = sp.integrate(kernel_shape**2, (w, -sp.oo, sp.oo))
    unique_kernel = (
        kernel_derivative_coefficient != 0
        and sp.simplify(kernel_equation_residual) == 0
        and w not in ratio_to_f0.free_symbols
        and positive(kernel_norm)
    )
    expect_bool(
        "PASS_UNIQUE_NORMALIZABLE_KERNEL",
        bool(unique_kernel),
        {
            "first_order_kernel_dimension": 1 if kernel_derivative_coefficient != 0 else 0,
            "integrating_factor_basis": kernel_shape,
            "basis/f0": ratio_to_f0,
            "basis_norm": kernel_norm,
        },
    )
    print("      ker(A)=span{f0}; the first-order equation supplies one integration constant")

    asymptote_coefficient = primitive(
        "PASS_POSITIVE_CONTINUUM_GAP", "V_asymptote_coefficient", sp.Integer(4)
    )
    gap_potential = (
        asymptote_coefficient - 6 / sp.cosh(w / ell) ** 2
    ) / ell**2
    gap_plus = sp.limit(gap_potential, w, sp.oo)
    gap_minus = sp.limit(gap_potential, w, -sp.oo)
    gap = 4 / ell**2
    expect_bool(
        "PASS_POSITIVE_CONTINUUM_GAP",
        sp.simplify(gap_plus - gap) == 0
        and sp.simplify(gap_minus - gap) == 0
        and positive(gap),
        {"V_H(+infinity)": gap_plus, "V_H(-infinity)": gap_minus, "gap": gap},
    )
    print("      unique normalizable zero mode lies below the continuum threshold 4/ell^2>0")


def run_norm_projection_and_relations() -> None:
    subbanner("Zero-mode norm, normalized projection, and reduction relations")
    w, ell = sp.symbols("w ell", positive=True, real=True)
    f0 = 1 / (ell * sp.cosh(w / ell) ** 2)
    expected_n0 = 8 / (3 * ell)

    norm_weight = primitive(
        "PASS_ZERO_MODE_NORM", "norm_inner_weight", sp.Integer(2)
    )
    n0_norm_test = sp.integrate(norm_weight * f0**2, (w, -sp.oo, sp.oo))
    expect_bool(
        "PASS_ZERO_MODE_NORM",
        sp.simplify(n0_norm_test - expected_n0) == 0 and positive(n0_norm_test),
        {"N0": n0_norm_test, "target": expected_n0},
    )

    canonical_n0 = sp.integrate(2 * f0**2, (w, -sp.oo, sp.oo))
    projection_weight = primitive(
        "PASS_NORMALIZED_ZERO_MODE_PROJECTION",
        "projection_inner_weight",
        sp.Integer(2),
    )
    amplitude = sp.symbols("H0", real=True)
    p0_f0 = sp.simplify(
        sp.integrate(projection_weight * f0 * f0, (w, -sp.oo, sp.oo))
        / canonical_n0
    )
    projected_h = sp.simplify(
        sp.integrate(projection_weight * f0 * (amplitude * f0), (w, -sp.oo, sp.oo))
        / canonical_n0
    )
    idempotence_residual = sp.simplify(p0_f0**2 - p0_f0)
    expect_bool(
        "PASS_NORMALIZED_ZERO_MODE_PROJECTION",
        sp.simplify(p0_f0 - 1) == 0
        and sp.simplify(projected_h - amplitude) == 0
        and idempotence_residual == 0,
        {
            "h=P0 H": projected_h,
            "P0 f0": p0_f0,
            "P0^2 f0-P0 f0": idempotence_residual,
        },
    )

    ell_star = sp.Rational(1, 20)
    n0_star = sp.simplify(canonical_n0.subs(ell, ell_star))
    base_m4_star = sp.Rational(3, 160)

    physical_m4_star = primitive(
        "PASS_PHYSICAL_H_NORM", "physical_M4_star", base_m4_star
    )
    physical_norm_star = sp.simplify(physical_m4_star * n0_star)
    expect_bool(
        "PASS_PHYSICAL_H_NORM",
        positive(physical_norm_star)
        and dim_residual(DIM_M4 * DIM_N0, DIM_MH) == 0,
        {"M4* N0*": physical_norm_star, "[M4 N0]": str(DIM_M4 * DIM_N0)},
    )

    inertia_m4_star = primitive(
        "PASS_REDUCED_H_INERTIA", "inertia_M4_star", base_m4_star
    )
    mh_from_reduction_star = sp.simplify(n0_star * inertia_m4_star)
    mh_declared_star = sp.Integer(1)
    expect_bool(
        "PASS_REDUCED_H_INERTIA",
        sp.simplify(mh_from_reduction_star - mh_declared_star) == 0
        and positive(mh_declared_star)
        and dim_residual(DIM_N0 * DIM_M4, DIM_MH) == 0,
        {
            "N0*": n0_star,
            "M4*": inertia_m4_star,
            "N0* M4*": mh_from_reduction_star,
            "M_h*": mh_declared_star,
        },
    )

    c_e_star = sp.Integer(1)
    k4_derived_star = sp.simplify(base_m4_star * c_e_star**2)
    parent_k4_dimension_raw = primitive(
        "PASS_PARENT_STIFFNESS_DIMENSIONAL_CONSISTENCY",
        "parent_K4_dimension",
        DIM_K4.vector(),
    )
    parent_k4_dimension = Dim(*parent_k4_dimension_raw)
    expect_bool(
        "PASS_PARENT_STIFFNESS_DIMENSIONAL_CONSISTENCY",
        dim_residual(DIM_M4 * DIM_CE**2, parent_k4_dimension) == 0,
        {
            "[M4 c_E^2]": str(DIM_M4 * DIM_CE**2),
            "declared [K4]": str(parent_k4_dimension),
            "value dependency": "K4=M4 c_E^2 is definitional",
        },
    )

    speed_kh_n0_power = primitive(
        "PASS_REDUCED_SPEED_PRESERVATION",
        "speed_Kh_N0_power",
        sp.Integer(1),
    )
    kh_from_projection_star = sp.simplify(
        n0_star**speed_kh_n0_power * k4_derived_star
    )
    parent_speed_squared_star = sp.simplify(k4_derived_star / base_m4_star)
    reduced_speed_squared_star = sp.simplify(
        kh_from_projection_star / mh_from_reduction_star
    )
    expect_bool(
        "PASS_REDUCED_SPEED_PRESERVATION",
        sp.simplify(parent_speed_squared_star - c_e_star**2) == 0
        and sp.simplify(reduced_speed_squared_star - parent_speed_squared_star) == 0
        and dim_residual(DIM_MH * DIM_CE**2, DIM_KH) == 0,
        {
            "c_E*^2": c_e_star**2,
            "K4*/M4*": parent_speed_squared_star,
            "K_h*/M_h*": reduced_speed_squared_star,
            "M_h*=N0* M4*": mh_from_reduction_star,
            "K_h*=N0*^p K4*": kh_from_projection_star,
            "p": speed_kh_n0_power,
        },
    )

    kh_declared_star = primitive(
        "PASS_REDUCED_H_STIFFNESS", "declared_Kh_star", sp.Integer(1)
    )
    expect_bool(
        "PASS_REDUCED_H_STIFFNESS",
        sp.simplify(kh_from_projection_star - kh_declared_star) == 0
        and positive(kh_declared_star)
        and dim_residual(DIM_N0 * DIM_K4, DIM_KH) == 0,
        {
            "N0* K4*": kh_from_projection_star,
            "K_h*": kh_declared_star,
        },
    )

    print("      clean star chain: N0*=160/3, M4*=K4*=3/160, M_h*=K_h*=c_E*=1")


def run_coupled_kernel() -> None:
    subbanner("Coupled (u_L,h) scalar kernel: independent stability and speed teeth")
    b_eff_star = sp.Integer(2)
    k_h_star = sp.Integer(1)

    c_hu_stability_star = primitive(
        "PASS_STABILITY", "stability_Chu_star", sp.Rational(1, 2)
    )
    d_star = sp.simplify(b_eff_star * k_h_star - c_hu_stability_star**2)
    stiffness_star = sp.Matrix(
        [[b_eff_star, c_hu_stability_star], [c_hu_stability_star, k_h_star]]
    )
    sylvester_ok = (
        positive(b_eff_star)
        and positive(d_star)
        and sp.simplify(d_star - sp.Rational(7, 4)) == 0
        and stiffness_star.is_positive_definite is True
        and dim_residual(DIM_BEFF * DIM_KH, DIM_D) == 0
        and dim_residual(DIM_CHU**2, DIM_D) == 0
    )
    expect_bool(
        "PASS_STABILITY",
        bool(sylvester_ok),
        {
            "D*=B_eff* K_h*-C_hu*^2": d_star,
            "physical [D]": str(DIM_D),
            "K*": stiffness_star,
        },
    )

    a_eff_star = primitive(
        "PASS_POSITIVE_GENERALIZED_WAVE_SPEEDS",
        "wave_Aeff_star",
        sp.Integer(1),
    )
    m_h_star = sp.Integer(1)
    c_hu_star = sp.Rational(1, 2)
    inertia_star = sp.diag(a_eff_star, m_h_star)
    wave_stiffness_star = sp.Matrix(
        [[b_eff_star, c_hu_star], [c_hu_star, k_h_star]]
    )
    z_star = sp.symbols("z_star", real=True)
    determinant_star = sp.factor((wave_stiffness_star - z_star * inertia_star).det())
    expected_determinant_star = z_star**2 - 3 * z_star + sp.Rational(7, 4)
    roots_star = sp.solve(determinant_star, z_star)
    expected_roots_star = {
        (sp.Integer(3) - sp.sqrt(2)) / 2,
        (sp.Integer(3) + sp.sqrt(2)) / 2,
    }
    unchanged_margin_star = sp.simplify(b_eff_star * k_h_star - c_hu_star**2)
    wave_ok = (
        inertia_star == sp.eye(2)
        and sp.simplify(determinant_star - expected_determinant_star) == 0
        and set(roots_star) == expected_roots_star
        and all(positive(root) for root in roots_star)
        and unchanged_margin_star == sp.Rational(7, 4)
        and positive(unchanged_margin_star)
    )
    expect_bool(
        "PASS_POSITIVE_GENERALIZED_WAVE_SPEEDS",
        bool(wave_ok),
        {
            "M*": inertia_star,
            "det(K*-z* M*)": determinant_star,
            "c_pm*^2": roots_star,
            "unchanged D*": unchanged_margin_star,
        },
    )
    print("      det(K*-z* M*)=z*^2-3 z*+7/4; c_pm*^2=(3 +/- sqrt(2))/2>0")


def run_masslessness_and_hessian() -> None:
    subbanner("Reduced masslessness and conservative Hessian")
    omega, k = sp.symbols("omega k", real=True)
    reduced_h_mass_star = primitive(
        "PASS_REDUCED_H_MASSLESSNESS", "reduced_h_mass_star", sp.Integer(0)
    )
    q_scalar = sp.Matrix(
        [
            [omega**2 - 2 * k**2, -sp.Rational(1, 2) * k**2],
            [
                -sp.Rational(1, 2) * k**2,
                omega**2 - k**2 - reduced_h_mass_star,
            ],
        ]
    )
    q_at_origin = q_scalar.subs({omega: 0, k: 0})
    expect_bool(
        "PASS_REDUCED_H_MASSLESSNESS",
        reduced_h_mass_star == 0 and matrix_zero(q_at_origin),
        {"k^0 h^2 coefficient*": reduced_h_mass_star, "Q_s(0,0)": q_at_origin},
    )

    grad_u, grad_h = sp.symbols("grad_u grad_h", real=True)
    conservative_energy = (
        grad_u**2 + sp.Rational(1, 2) * grad_h**2
        + sp.Rational(1, 2) * grad_u * grad_h
    )
    action_hessian = sp.hessian(conservative_energy, (grad_u, grad_h))
    u_cross_star = primitive(
        "PASS_CONSERVATIVE_HESSIAN_SYMMETRY",
        "u_equation_cross_star",
        sp.Rational(1, 2),
    )
    displayed_euler_operator = sp.Matrix(
        [[2, u_cross_star], [sp.Rational(1, 2), 1]]
    )
    expect_bool(
        "PASS_CONSERVATIVE_HESSIAN_SYMMETRY",
        action_hessian == action_hessian.T
        and displayed_euler_operator == displayed_euler_operator.T
        and displayed_euler_operator == action_hessian,
        {
            "action Hessian": action_hessian,
            "displayed Euler operator": displayed_euler_operator,
        },
    )


def run_dimensional_firewall() -> None:
    subbanner("Units-restored electric-scalar dimensional firewall [L,T,M]")
    firewall_chu_raw = primitive(
        "PASS_DIMENSIONAL_HOMOGENEITY",
        "firewall_Chu_dimension",
        DIM_CHU.vector(),
    )
    firewall_chu = Dim(*firewall_chu_raw)

    terms = {
        "H_kin": DIM_M4 * (DIM_H / TIME) ** 2,
        "H_grad": DIM_K4 * (DIM_H / LENGTH) ** 2,
        "H_pot": DIM_K4 * DIM_H * (DIM_H / LENGTH**2),
        "u_kin": DIM_AEFF * (DIM_u / TIME) ** 2,
        "u_stiff": DIM_BEFF * (DIM_u / LENGTH) ** 2,
        "h_kin": DIM_MH * (DIM_h / TIME) ** 2,
        "h_stiff": DIM_KH * (DIM_h / LENGTH) ** 2,
        "mix": firewall_chu * (DIM_u / LENGTH) * (DIM_h / LENGTH),
    }
    targets = {
        **{name: BULK_DENSITY for name in ("H_kin", "H_grad", "H_pot")},
        **{
            name: BRANE_DENSITY
            for name in ("u_kin", "u_stiff", "h_kin", "h_stiff", "mix")
        },
    }
    residuals = {name: dim_residual(terms[name], targets[name]) for name in terms}
    expect_zero(
        "PASS_DIMENSIONAL_HOMOGENEITY",
        sp.simplify(sum(residuals.values())),
        {name: {"actual": str(terms[name]), "target": str(targets[name])} for name in terms},
    )
    print("      checked terms={H_kin,H_grad,H_pot,u_kin,u_stiff,h_kin,h_stiff,mix}")


def print_scope_and_verdict() -> None:
    subbanner("Scope and verdict")
    print("SCOPE: static algebraic prerequisites of the charge electric-scalar sector only.")
    print("  DEFERRED: mouth SOURCE mechanism -> stage031.")
    print("  EXCLUDED: assembled one-body/two-body/far-field claims -> Part VII.")
    print("  OUT-OF-SCOPE / UNDISCHARGED CROSS-SECTOR G0 PREDICATE (Part VII de-dup obligation):")
    print("    planar-wall factorization; bulk U/Fourier/phase predicates; drain kernel/controllers.")
    print("VERDICT_TOKEN: ELECTRIC_SCALAR_CLOSURE_STATIC")


def main() -> None:
    if ACTIVE_MUTATION and ACTIVE_MUTATION not in ABLATIONS:
        print(f"FIRST_FAILURE=UNKNOWN_MUTATION")
        print(f"FAIL  UNKNOWN_MUTATION: {ACTIVE_MUTATION!r}")
        raise AuditFailure("UNKNOWN_MUTATION", ACTIVE_MUTATION)

    banner("ledger_stage030_electric_scalar_localized_h_closure SymPy audit")
    if ACTIVE_MUTATION:
        ablation = ABLATIONS[ACTIVE_MUTATION]
        print(f"ACTIVE_MUTATION={ACTIVE_MUTATION}")
        print(f"MUTATED_CARD_PRIMITIVE={ablation.description}")

    run_spectral_structure()
    run_norm_projection_and_relations()
    run_coupled_kernel()
    run_masslessness_and_hessian()
    run_dimensional_firewall()
    print_scope_and_verdict()

    if ACTIVE_MUTATION:
        print(f"FIRST_FAILURE=MUTATION_DID_NOT_FIRE")
        print(f"FAIL  MUTATION_DID_NOT_FIRE: {ACTIVE_MUTATION}")
        raise AuditFailure("MUTATION_DID_NOT_FIRE", ACTIVE_MUTATION)

    print("")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print("OVERALL PASS: SymPy verified stage030 electric-scalar localized-H static closure")


if __name__ == "__main__":
    try:
        main()
    except AuditFailure as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
        print(f"OVERALL FAIL: SymPy stage030 audit did not close ({exc.predicate})")
        raise SystemExit(1)
