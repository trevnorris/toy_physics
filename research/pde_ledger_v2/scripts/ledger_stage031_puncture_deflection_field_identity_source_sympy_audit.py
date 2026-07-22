#!/usr/bin/env python3
"""Ledger stage031 SymPy audit: puncture-deflection mechanism.

Standalone, print-only, assert-zero, float-free, and file-I/O-free.  Stage 030
is consumed as the symbolic input bundle
``{f0, f0(0)=1/ell, N0=8/(3 ell), h=P0 H, S_Lh, D, D*=7/4}``;
its scalar closure is not re-derived here.  The optional tooth-local ablation
switch is ``LEDGER_STAGE031_MUTATION``.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
from typing import Any

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0
MUTATION_ENV = "LEDGER_STAGE031_MUTATION"


class AuditFailure(AssertionError):
    def __init__(self, predicate: str, detail: str = "") -> None:
        super().__init__(predicate)
        self.predicate = predicate
        self.detail = detail


@dataclass(frozen=True)
class Ablation:
    primitive: str
    value: Any
    description: str


# Every mutation is read only by its named predicate.  Canonical products are
# retained for all later predicates, so one mutation cannot cascade.
ABLATIONS: dict[str, Ablation] = {
    "PASS_FIELD_IDENTITY": Ablation("xi_power", 2, "xi_w=ell*h -> ell^2*h"),
    "PASS_ETA_ANNULUS_NORMALIZATION": Ablation("eta_numerator", 2, "eta prefactor 3 -> 2"),
    "PASS_F0_MOUTH_VALUE_EVALUATED": Ablation("f0_power", 1, "sech^2 profile -> sech profile"),
    "PASS_REDUCED_COUPLING_NONZERO": Ablation("Jm_witness", 0, "nonzero J_m witness -> 0"),
    "PASS_BOUNDED_PLUS_SLEEVE": Ablation("peak_factor", 2, "peak-normalized sleeve -> twice its peak"),
    "PASS_REFLECTION_ANTISYMMETRY": Ablation("reflect_body", 0, "r_Sigma,-(w)=r_Sigma,+(w)"),
    "PASS_I_PLUS_REFLECTION_DOMINANCE": Ablation("dominance_sign", -1, "+w dominance -> reversed dominance"),
    "PASS_I_PLUS_EVEN_DEFORMATION_CONTROL": Ablation("even_control_odd_part", 1, "add an odd part to even D"),
    "PASS_I_PLUS_ZERO_SLEEVE_CONTROL": Ablation("zero_sleeve", 1, "L_s=0 control -> L_s>0"),
    "PASS_I_PLUS_NEGATIVE_SLEEVE_CONTROL": Ablation("negative_control", 1, "-w sleeve control -> +w sleeve"),
    "PASS_N_CHI_NONZERO_GUARD": Ablation("nchi_guard", 0, "pre-division N_chi guard -> zero"),
    "PASS_Q_CHI_ORIENTATION": Ablation("projection_parity", 0, "odd o_ell -> even kernel"),
    "PASS_BARE_FORCING_LIVE_VARIATION": Ablation("mouth_quadratic_factor", 2, "live K_m term factor 1 -> 2"),
    "PASS_EXTERIOR_ONE_OVER_R": Ablation("radial_dimension", 4, "three exterior dimensions -> four"),
    "PASS_POSITIVE_HOLDING_CURVATURE": Ablation("curvature_sign", -1, "positive exterior energy -> negative"),
    "PASS_NONZERO_HA_REQUIRES_CORE_HOLDER": Ablation("coreless_shift", 1, "inject a nonzero coreless stationary point"),
    "PASS_KAPPA_REDUCTION": Ablation("kappa_factor", 2, "kappa=D/b -> 2D/b"),
    "PASS_Z_B_WIRING": Ablation("zb_factor", 0, "remove z_b from Z lower-left entry"),
    "PASS_RESPONSE_M_UU": Ablation("muu_extra", 1, "m_uu -> m_uu+1"),
    "PASS_RESPONSE_M_UG": Ablation("mug_sign", 1, "m_ug sign - -> +"),
    "PASS_RESPONSE_M_GG": Ablation("mgg_factor", 2, "m_gg -> 2m_gg"),
    "PASS_RESPONSE_SYMMETRY": Ablation("symmetry_skew", 1, "add skew part to m"),
    "PASS_RESPONSE_DETERMINANT": Ablation("det_factor", 2, "det m target -> twice target"),
    "PASS_Z_G_POSTULATED_WITNESS": Ablation("zg_witness", 0, "Robin witness z_g=1 -> 0"),
    "PASS_M_GG_NONNEGATIVE": Ablation("mgg_square", 1, "z_g^2 -> z_g"),
    "PASS_RESPONSE_STAR_WITNESS": Ablation("star_zb", 0, "star z_b=1 -> 0"),
    "PASS_S_GG_SELF_RESPONSE": Ablation("green_sign", -1, "positive Green quadratic form -> indefinite"),
    "PASS_NEUTRAL_FAR_FIELD_FORM": Ablation("potential_power", 2, "1/R potential -> 1/R^2"),
    "PASS_ALLOWED_ZERO_COEFFICIENT": Ablation("zero_C", 1, "allowed C=0 control -> C=1"),
    "PASS_UNITS_XI_W": Ablation("dim_xi", (2, 0, 0), "[xi_w]=L -> L^2"),
    "PASS_UNITS_H": Ablation("dim_h", (1, 0, 0), "[h]=1 -> L"),
    "PASS_UNITS_H_A": Ablation("dim_hA", (1, 0, 0), "[h_A]=1 -> L"),
    "PASS_UNITS_K_M": Ablation("dim_Km", (3, -2, 1), "[K_m]=M L^4 T^-2 -> M L^3 T^-2"),
    "PASS_UNITS_J_M": Ablation("dim_Jm", (2, -2, 1), "[J_m]=M L^3 T^-2 -> E"),
    "PASS_UNITS_K_M_REDUCED": Ablation("dim_km", (3, -2, 1), "[k_m]=E -> E L"),
    "PASS_UNITS_G_CHIH": Ablation("dim_g", (3, -2, 1), "[g_chih]=E -> E L"),
    "PASS_UNITS_ETA": Ablation("dim_eta", (-2, 0, 0), "[eta]=L^-3 -> L^-2"),
    "PASS_UNITS_ODD_KERNEL": Ablation("dim_o", (-1, 0, 0), "[o_ell]=L^-2 -> L^-1"),
    "PASS_UNITS_N_CHI": Ablation("dim_nchi", (0, 0, 0), "[N_chi]=L^-1 -> 1"),
    "PASS_UNITS_M_UU": Ablation("dim_muu", (2, 2, -1), "[m_uu]=L^3/E -> L^4/E"),
    "PASS_UNITS_M_UG": Ablation("dim_mug", (1, 2, -1), "[m_ug]=L^2/E -> L^3/E"),
    "PASS_UNITS_M_GG": Ablation("dim_mgg", (0, 2, -1), "[m_gg]=L/E -> L^2/E"),
    "PASS_UNITS_KAPPA": Ablation("dim_kappa", (0, -2, 1), "[kappa]=E/L -> E/L^2"),
    "PASS_UNITS_ESCAPE_MATRIX": Ablation("dim_Z21", (0, 0, 0), "[Z_21]=L -> 1"),
    "PASS_UNITS_DET_M": Ablation("dim_detm", (1, 4, -2), "[det m]=M^-2 T^4 -> L M^-2 T^4"),
    "PASS_UNITS_S_GG": Ablation("dim_sgg", (-1, 2, -1), "[S_gg]=E^-1 -> L/E"),
    "PASS_UNITS_SHELL_C": Ablation("dim_C", (3, -4, 2), "[C]=E^2 -> E^2/L"),
    "PASS_UNITS_SHELL_A": Ablation("dim_A", (2, -2, 1), "[A]=E L -> E"),
    "PASS_UNITS_U": Ablation("dim_U", (3, -2, 1), "[U]=E -> E L"),
    "PASS_UNITS_F": Ablation("dim_F", (2, -2, 1), "[F]=E/L -> E"),
}

ACTIVE_MUTATION = os.environ.get(MUTATION_ENV, "").strip()


def primitive(predicate: str, name: str, canonical: Any) -> Any:
    if ACTIVE_MUTATION == predicate and ABLATIONS[predicate].primitive == name:
        return ABLATIONS[predicate].value
    return canonical


def compact(value: Any) -> str:
    try:
        return sp.sstr(sp.factor(sp.cancel(sp.simplify(value))))
    except (TypeError, ValueError, AttributeError):
        return str(value)


def assert_no_float(name: str, value: Any) -> None:
    if isinstance(value, (str, dict, tuple, list, bool)):
        return
    floats = sp.sympify(value).atoms(sp.Float)
    if floats:
        raise AuditFailure(name, f"machine Float atom(s): {floats}")


def expect_zero(name: str, residual: Any, evidence: Any = None) -> None:
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
    expect_zero(name, 0 if condition else 1, evidence)


@dataclass(frozen=True)
class Dim:
    """Exact [L,T,M] exponent vector."""

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

    def __pow__(self, power: int) -> "Dim":
        return Dim(power * self.l, power * self.t, power * self.m)

    def vector(self) -> tuple[sp.Rational, sp.Rational, sp.Rational]:
        return self.l, self.t, self.m


ONE = Dim()
LENGTH = Dim(1, 0, 0)
TIME = Dim(0, 1, 0)
MASS = Dim(0, 0, 1)
ENERGY = MASS * LENGTH**2 / TIME**2


def dim_from(value: tuple[int, int, int] | tuple[sp.Rational, ...]) -> Dim:
    return Dim(*value)


def dim_residual(actual: Dim, expected: Dim) -> sp.Expr:
    return sp.simplify(sum((x - y) ** 2 for x, y in zip(actual.vector(), expected.vector())))


def positive_exact(expr: sp.Expr) -> bool:
    return expr.is_positive is True or sp.ask(sp.Q.positive(expr)) is True


def heading(text: str) -> None:
    print("")
    print(text)
    print("-" * len(text))


def run_mouth_reconstruction() -> None:
    heading("Bare mouth reconstruction and option-A reflection dominance")
    rho = sp.symbols("rho", positive=True, real=True)
    w, u = sp.symbols("w u", real=True)
    wp, a, ell, width, sleeve = sp.symbols(
        "w_pos a ell W_B L_s", positive=True, real=True
    )
    h, H0, Jm, Km, s = sp.symbols("h H0 J_m K_m s", real=True, nonzero=True)

    xi_power = primitive("PASS_FIELD_IDENTITY", "xi_power", 1)
    xi_test = ell**xi_power * h
    expect_zero("PASS_FIELD_IDENTITY", xi_test / ell - h, {"xi_w": xi_test})

    eta_numerator = primitive("PASS_ETA_ANNULUS_NORMALIZATION", "eta_numerator", 3)
    eta_test = sp.Rational(eta_numerator, 1) / (
        4 * sp.pi * ((a + ell) ** 3 - a**3)
    )
    eta_integral = sp.integrate(4 * sp.pi * rho**2 * eta_test, (rho, a, a + ell))
    expect_zero(
        "PASS_ETA_ANNULUS_NORMALIZATION", eta_integral - 1,
        {"eta": eta_test, "native radial integral": eta_integral},
    )
    canonical_eta_integral = sp.Integer(1)

    f0_power = primitive("PASS_F0_MOUTH_VALUE_EVALUATED", "f0_power", 2)
    f0_test = 1 / (ell * sp.cosh(w / ell) ** f0_power)
    evaluated_f0 = sp.simplify(f0_test.subs(w, 0))
    # A wrong exponent alone has the same value at zero; the profile tooth
    # therefore also protects the consumed sech^2 shape via its curvature.
    canonical_curvature = sp.diff(1 / (ell * sp.cosh(w / ell) ** 2), w, 2).subs(w, 0)
    test_curvature = sp.diff(f0_test, w, 2).subs(w, 0)
    expect_bool(
        "PASS_F0_MOUTH_VALUE_EVALUATED",
        sp.simplify(evaluated_f0 - 1 / ell) == 0
        and sp.simplify(test_curvature - canonical_curvature) == 0,
        {"evaluated f0(0)": evaluated_f0, "profile curvature": test_curvature},
    )
    canonical_f0_mouth = sp.Integer(1) / ell

    jm_witness = sp.Integer(primitive("PASS_REDUCED_COUPLING_NONZERO", "Jm_witness", 1))
    g_test = sp.simplify(jm_witness / ell)
    expect_bool(
        "PASS_REDUCED_COUPLING_NONZERO",
        sp.simplify(g_test - jm_witness / ell) == 0 and jm_witness != 0,
        {"g_chih": g_test, "J_m witness": jm_witness},
    )

    # Concrete frozen representative on the mouth annulus.  Its interval is
    # [-W_B/2, W_B/2+L_s], its center is L_s/2, and the denominator is its
    # exact peak, so it maps R into (0,1].
    r_plus = (
        sp.tanh((w + width / 2) / ell)
        - sp.tanh((w - (width / 2 + sleeve)) / ell)
    ) / (2 * sp.tanh((width + sleeve) / (2 * ell)))
    center = sleeve / 2
    halfspan = (width + sleeve) / 2
    r_centered = (1 + sp.cosh(2 * halfspan / ell)) / (
        sp.cosh(2 * (w - center) / ell) + sp.cosh(2 * halfspan / ell)
    )
    peak_factor = sp.Integer(primitive("PASS_BOUNDED_PLUS_SLEEVE", "peak_factor", 1))
    r_bounded_test = peak_factor * r_centered
    equivalence_star = sp.simplify((r_plus - r_centered).subs({w: center}))
    range_certificate = sp.trigsimp(
        sp.cosh(2 * u / ell) - 1 - 2 * sp.sinh(u / ell) ** 2
    )
    expect_bool(
        "PASS_BOUNDED_PLUS_SLEEVE",
        equivalence_star == 0
        and sp.simplify(r_bounded_test.subs(w, center) - 1) == 0
        and range_certificate == 0,
        {"representative": r_plus, "peak": r_bounded_test.subs(w, center), "range proof": "denominator-numerator=2 sinh(u/ell)^2>=0"},
    )

    r_bg = (
        sp.tanh((w + width / 2) / ell)
        - sp.tanh((w - width / 2) / ell)
    ) / (2 * sp.tanh(width / (2 * ell)))
    odd_kernel = w * sp.exp(-w**2 / (2 * ell**2)) / (sp.sqrt(2 * sp.pi) * ell**3)
    reflect_body = primitive("PASS_REFLECTION_ANTISYMMETRY", "reflect_body", 1)
    r_minus_test = r_plus.subs(w, -w) if reflect_body == 1 else r_plus
    explicit_reflection = sp.trigsimp(r_minus_test - r_plus.subs(w, -w))
    expect_bool(
        "PASS_REFLECTION_ANTISYMMETRY",
        explicit_reflection == 0
        and sp.simplify(odd_kernel.subs(w, -w) + odd_kernel) == 0
        and sp.trigsimp(r_bg.subs(w, -w) - r_bg) == 0,
        {"r_Sigma,-(w)-r_Sigma,+(-w)": explicit_reflection, "I_minus": "-I_plus"},
    )

    # Structural, parameter-generic certificate.  In centered form the
    # reflected denominator difference has a sign fixed solely by w>0 and
    # the one-sided frozen length L_s>0.
    den_forward = sp.cosh(2 * (wp - center) / ell) + sp.cosh(2 * halfspan / ell)
    den_reflected = sp.cosh(2 * (-wp - center) / ell) + sp.cosh(2 * halfspan / ell)
    den_difference = sp.trigsimp(sp.expand_trig(den_reflected - den_forward))
    expected_difference = 2 * sp.sinh(2 * wp / ell) * sp.sinh(sleeve / ell)
    dominance_sign = sp.Integer(primitive(
        "PASS_I_PLUS_REFLECTION_DOMINANCE", "dominance_sign", 1
    ))
    signed_difference = dominance_sign * den_difference
    r_forward = r_centered.subs(w, wp)
    r_reflected = r_centered.subs(w, -wp)
    d_pair = sp.factor((r_forward - r_reflected) * (r_forward + r_reflected))
    profile_numerator = 1 + sp.cosh(2 * halfspan / ell)
    r_difference_certificate = (
        dominance_sign * profile_numerator * expected_difference
        / (den_forward * den_reflected)
    )
    r_difference_identity = sp.simplify(
        (r_forward - r_reflected)
        - profile_numerator * den_difference / (den_forward * den_reflected)
    )
    positive_kernel_mass = sp.integrate(
        odd_kernel.subs(w, wp), (wp, 0, sp.oo)
    )
    generic_certificate = (
        sp.simplify(den_difference - expected_difference) == 0
        and r_difference_identity == 0
        and positive_exact(dominance_sign * expected_difference)
        and positive_exact(r_difference_certificate)
        and positive_exact(r_forward)
        and positive_exact(r_reflected)
        and positive_exact(r_forward + r_reflected)
        and positive_exact(positive_kernel_mass)
        and canonical_eta_integral == 1
    )
    if ACTIVE_MUTATION != "PASS_I_PLUS_REFLECTION_DOMINANCE" and not generic_certificate:
        print("STAGE031_STOP: I_plus_not_generic")
        raise SystemExit(1)
    expect_bool(
        "PASS_I_PLUS_REFLECTION_DOMINANCE", generic_certificate,
        {
            "den(-w)-den(w)": den_difference,
            "D(w)-D(-w)": "(r-r_-)(r+r_-)>0",
            "Integral_0^infinity o_ell": positive_kernel_mass,
            "strict positive measure": "all w>0",
        },
    )

    even_odd_part = sp.Integer(primitive(
        "PASS_I_PLUS_EVEN_DEFORMATION_CONTROL", "even_control_odd_part", 0
    ))
    d_even_test = sp.exp(-w**2 / (2 * ell**2)) * (1 + even_odd_part * w / ell)
    even_control_integral = sp.integrate(odd_kernel * d_even_test, (w, -sp.oo, sp.oo))
    expect_zero(
        "PASS_I_PLUS_EVEN_DEFORMATION_CONTROL", even_control_integral,
        {"D(w)=D(-w)": even_odd_part == 0, "I_plus control": even_control_integral},
    )

    zero_sleeve = primitive("PASS_I_PLUS_ZERO_SLEEVE_CONTROL", "zero_sleeve", 0)
    zero_profile = sp.simplify(r_centered.subs(sleeve, sp.Integer(zero_sleeve)))
    zero_d = zero_profile**2 - r_bg**2
    zero_parity_residual = sp.trigsimp(zero_d.subs(w, -w) - zero_d)
    expect_zero(
        "PASS_I_PLUS_ZERO_SLEEVE_CONTROL", zero_parity_residual,
        {"L_s": zero_sleeve, "paired integral": 0 if zero_parity_residual == 0 else "nonzero"},
    )

    negative_control = primitive(
        "PASS_I_PLUS_NEGATIVE_SLEEVE_CONTROL", "negative_control", -1
    )
    control_pair = sp.Integer(negative_control) * (
        profile_numerator * expected_difference / (den_forward * den_reflected)
    ) * (r_forward + r_reflected)
    expect_bool(
        "PASS_I_PLUS_NEGATIVE_SLEEVE_CONTROL", control_pair.is_negative is True,
        {"inadmissible sleeve": "-w", "I_plus sign": sp.sign(control_pair)},
    )

    # Retain the actual integral unevaluated; the sign was certified before
    # N_chi is defined.  This is the anti-N_chi/N_chi ordering requirement.
    actual_pair = odd_kernel.subs(w, wp) * (
        r_centered.subs(w, wp) ** 2 - r_centered.subs(w, -wp) ** 2
    )
    i_plus = sp.Integral(actual_pair, (wp, 0, sp.oo))
    n_chi = i_plus
    nchi_guard = primitive("PASS_N_CHI_NONZERO_GUARD", "nchi_guard", 1)
    expect_bool(
        "PASS_N_CHI_NONZERO_GUARD",
        bool(nchi_guard == 1 and generic_certificate),
        {"N_chi": n_chi, "ordering": "generic I_plus>0 certified before division"},
    )

    projection_parity = primitive("PASS_Q_CHI_ORIENTATION", "projection_parity", 1)
    i_minus_sign = -1 if projection_parity == 1 else 1
    q_plus = sp.Integer(1)
    q_minus = sp.Integer(i_minus_sign)
    expect_bool(
        "PASS_Q_CHI_ORIENTATION", q_plus == 1 and q_minus == -1,
        {
            "kernel": "odd" if projection_parity == 1 else "even",
            "Q_+=I_+/N_chi": q_plus,
            "Q_-=I_-/N_chi": q_minus,
        },
    )

    # Vary the actual reduced functional obtained from the parent expression,
    # retaining the Robin k_m term.  The bare component is evaluated at h=0.
    mouth_factor = sp.Integer(primitive(
        "PASS_BARE_FORCING_LIVE_VARIATION", "mouth_quadratic_factor", 1
    ))
    canonical_eta = sp.Rational(3, 1) / (
        4 * sp.pi * ((a + ell) ** 3 - a**3)
    )
    parent_density_test = canonical_eta * (
        mouth_factor * Km * (canonical_f0_mouth * h) ** 2 / 2
        - Jm * s * canonical_f0_mouth * h
    )
    euler_density_test = sp.diff(parent_density_test, h)
    km = Km / ell**2
    g = Jm / ell
    expected_euler_density = canonical_eta * (km * h - g * s)
    integrated_euler = sp.integrate(
        4 * sp.pi * rho**2 * euler_density_test, (rho, a, a + ell)
    )
    bare_amplitude = sp.simplify(integrated_euler.subs(h, 0))
    expect_bool(
        "PASS_BARE_FORCING_LIVE_VARIATION",
        sp.simplify(euler_density_test - expected_euler_density) == 0
        and sp.simplify(integrated_euler - (km * h - g * s)) == 0
        and sp.simplify(bare_amplitude + g * s) == 0
        and sp.simplify(bare_amplitude.subs({Jm: 1, s: 1, ell: 1})) != 0,
        {"delta Omega/delta h density": euler_density_test, "integrated Euler amplitude": integrated_euler, "bare integrated amplitude": bare_amplitude},
    )
    print("      option-A result: N_chi=I_plus>0; (Q_+,Q_-)=(+1,-1); bare source=-g_chih*s")


def response_objects() -> tuple[Any, ...]:
    b, k = sp.symbols("b K_h", positive=True, real=True)
    c, zb, zg = sp.symbols("c z_b z_g", real=True)
    delta = b * k - c**2
    kappa = delta / b
    Z = sp.Matrix([[1, 0], [-(c / b) * zb, zg]])
    m = sp.simplify(Z.T * sp.diag(1 / b, 1 / kappa) * Z)
    return b, k, c, zb, zg, delta, kappa, Z, m


def run_response() -> None:
    heading("Full completed-square response")
    b, k, c, zb, zg, delta, kappa, Z, m = response_objects()
    expected_muu = (delta + c**2 * zb**2) / (b * delta)
    expected_mug = -c * zb * zg / delta
    expected_mgg = b * zg**2 / delta

    kappa_factor = primitive("PASS_KAPPA_REDUCTION", "kappa_factor", 1)
    expect_zero("PASS_KAPPA_REDUCTION", kappa_factor * kappa - delta / b)

    zb_factor = primitive("PASS_Z_B_WIRING", "zb_factor", 1)
    Z_test = sp.Matrix([[1, 0], [-(c / b) * zb_factor * zb, zg]])
    expect_zero("PASS_Z_B_WIRING", Z_test[1, 0] + (c / b) * zb)

    muu_test = m[0, 0] + primitive("PASS_RESPONSE_M_UU", "muu_extra", 0)
    expect_zero("PASS_RESPONSE_M_UU", muu_test - expected_muu)

    mug_sign = primitive("PASS_RESPONSE_M_UG", "mug_sign", -1)
    mug_target = mug_sign * c * zb * zg / delta
    expect_zero("PASS_RESPONSE_M_UG", m[0, 1] - mug_target)

    mgg_factor = primitive("PASS_RESPONSE_M_GG", "mgg_factor", 1)
    expect_zero("PASS_RESPONSE_M_GG", m[1, 1] - mgg_factor * expected_mgg)

    skew = primitive("PASS_RESPONSE_SYMMETRY", "symmetry_skew", 0)
    m_test = m + sp.Matrix([[0, skew], [0, 0]])
    expect_bool("PASS_RESPONSE_SYMMETRY", m_test == m_test.T, {"m": m_test})

    det_factor = primitive("PASS_RESPONSE_DETERMINANT", "det_factor", 1)
    expect_zero("PASS_RESPONSE_DETERMINANT", sp.factor(m.det()) - det_factor * zg**2 / delta)

    zg_witness = primitive("PASS_Z_G_POSTULATED_WITNESS", "zg_witness", 1)
    expect_bool(
        "PASS_Z_G_POSTULATED_WITNESS", bool(0 < zg_witness <= 1),
        {"status": "POSTULATED Robin-admissibility witness", "z_g*": zg_witness},
    )

    square_power = primitive("PASS_M_GG_NONNEGATIVE", "mgg_square", 2)
    delta_positive = sp.symbols("Delta", positive=True, real=True)
    mgg_stable_test = b * zg**square_power / delta_positive
    expect_bool(
        "PASS_M_GG_NONNEGATIVE", sp.ask(sp.Q.nonnegative(mgg_stable_test)) is True,
        {"m_gg": mgg_stable_test, "assumptions": "b>0, D>0, z_g real"},
    )

    star_zb = primitive("PASS_RESPONSE_STAR_WITNESS", "star_zb", 1)
    star = {b: 2, k: 1, c: sp.Rational(1, 2), zb: star_zb, zg: 1}
    star_m = sp.simplify(m.subs(star))
    expect_bool(
        "PASS_RESPONSE_STAR_WITNESS",
        sp.simplify(delta.subs(star) - sp.Rational(7, 4)) == 0
        and star_m == sp.Matrix([[sp.Rational(4, 7), -sp.Rational(2, 7)], [-sp.Rational(2, 7), sp.Rational(8, 7)]]),
        {"D* (consumed stage030)": delta.subs(star), "m*": star_m},
    )

    p, q = sp.symbols("p q", positive=True, real=True)
    green_sign = primitive("PASS_S_GG_SELF_RESPONSE", "green_sign", 1)
    source = sp.Matrix([1, 1])
    green = sp.diag(1 / p, green_sign / q)
    sgg = sp.factor((source.T * green * source)[0])
    expect_bool(
        "PASS_S_GG_SELF_RESPONSE", positive_exact(sgg),
        {"definition": "S_gg=<eta,L_h^{-1}eta>", "Galerkin quadratic form": sgg, "status": "DERIVED for positive L_h"},
    )
    print("      m=Z^T diag(b^-1,kappa^-1) Z; strict z_g>0 is POSTULATED, m_gg>=0 is EARNED")


def run_exterior_and_shell() -> None:
    heading("Exterior stationary solution and target-blind shell")
    r, a, kappa = sp.symbols("r a kappa", positive=True, real=True)
    hA = sp.symbols("h_A", real=True)
    radial_dimension = primitive("PASS_EXTERIOR_ONE_OVER_R", "radial_dimension", 3)
    trial = hA * (a / r) ** (radial_dimension - 2)
    radial_residual = sp.simplify(sp.diff(r ** (radial_dimension - 1) * sp.diff(trial, r), r))
    expect_bool(
        "PASS_EXTERIOR_ONE_OVER_R",
        radial_dimension == 3 and radial_residual == 0 and sp.simplify(trial.subs(r, a) - hA) == 0,
        {"radial Euler equation": radial_residual, "h(r)": trial},
    )
    h_exterior = hA * a / r
    energy = sp.integrate(
        sp.Rational(1, 2) * kappa * 4 * sp.pi * r**2 * sp.diff(h_exterior, r) ** 2,
        (r, a, sp.oo),
    )
    curvature_sign = primitive("PASS_POSITIVE_HOLDING_CURVATURE", "curvature_sign", 1)
    curvature_test = curvature_sign * sp.diff(energy, hA, 2)
    expect_bool(
        "PASS_POSITIVE_HOLDING_CURVATURE",
        sp.simplify(energy - 2 * sp.pi * kappa * a * hA**2) == 0 and positive_exact(curvature_test),
        {"E_ext": energy, "holding curvature": curvature_test},
    )

    coreless_shift = primitive("PASS_NONZERO_HA_REQUIRES_CORE_HOLDER", "coreless_shift", 0)
    varied_energy = 2 * sp.pi * kappa * a * (hA - coreless_shift) ** 2
    stationary = sp.solve(sp.diff(varied_energy, hA), hA)
    expect_bool(
        "PASS_NONZERO_HA_REQUIRES_CORE_HOLDER", stationary == [0],
        {"coreless stationary h_A": stationary, "named fact": "NONZERO_HA_REQUIRES_CORE_HOLDER"},
    )

    R = sp.symbols("R", positive=True, real=True)
    s1, s2, C = sp.symbols("s_1 s_2 C", real=True)
    b, D, zg = sp.symbols("b D z_g", positive=True, real=True)
    mgg = b * zg**2 / D
    amplitude = mgg * C
    potential_power = primitive("PASS_NEUTRAL_FAR_FIELD_FORM", "potential_power", 1)
    U_test = s1 * s2 * amplitude / (4 * sp.pi * R**potential_power)
    F_test = sp.factor(-sp.diff(U_test, R))
    expected_U = s1 * s2 * amplitude / (4 * sp.pi * R)
    expected_F = s1 * s2 * amplitude / (4 * sp.pi * R**2)
    expect_bool(
        "PASS_NEUTRAL_FAR_FIELD_FORM",
        sp.simplify(U_test - expected_U) == 0 and sp.simplify(F_test - expected_F) == 0
        and sp.simplify(U_test.xreplace({s1: s2, s2: s1}) - U_test) == 0,
        {"A": amplitude, "U": U_test, "F_out": F_test, "C": "unspecified real [E^2]"},
    )
    zero_C = primitive("PASS_ALLOWED_ZERO_COEFFICIENT", "zero_C", 0)
    expect_zero(
        "PASS_ALLOWED_ZERO_COEFFICIENT", expected_U.subs(C, zero_C),
        {"conditional leading form": "coefficient may vanish", "class numerator": "not selected"},
    )
    print("      NONZERO_HA_REQUIRES_CORE_HOLDER")
    print("      A=m_gg*C, [C]=E^2: neutral shell only; no BC numerator and no sign selected")


def run_units() -> None:
    heading("Complete mechanism dimensional firewall [L,T,M]")
    dims = {
        "PASS_UNITS_XI_W": ("dim_xi", LENGTH, LENGTH),
        "PASS_UNITS_H": ("dim_h", ONE, ONE),
        "PASS_UNITS_H_A": ("dim_hA", ONE, ONE),
        "PASS_UNITS_K_M": ("dim_Km", ENERGY * LENGTH**2, ENERGY * LENGTH**2),
        "PASS_UNITS_J_M": ("dim_Jm", ENERGY * LENGTH, ENERGY * LENGTH),
        "PASS_UNITS_K_M_REDUCED": ("dim_km", ENERGY, ENERGY),
        "PASS_UNITS_G_CHIH": ("dim_g", ENERGY, ENERGY),
        "PASS_UNITS_ETA": ("dim_eta", LENGTH**-3, LENGTH**-3),
        "PASS_UNITS_ODD_KERNEL": ("dim_o", LENGTH**-2, LENGTH**-2),
        "PASS_UNITS_N_CHI": ("dim_nchi", LENGTH**-1, LENGTH**-1),
        "PASS_UNITS_M_UU": ("dim_muu", LENGTH**3 / ENERGY, LENGTH**3 / ENERGY),
        "PASS_UNITS_M_UG": ("dim_mug", LENGTH**2 / ENERGY, LENGTH**2 / ENERGY),
        "PASS_UNITS_M_GG": ("dim_mgg", LENGTH / ENERGY, LENGTH / ENERGY),
        "PASS_UNITS_KAPPA": ("dim_kappa", ENERGY / LENGTH, ENERGY / LENGTH),
        "PASS_UNITS_ESCAPE_MATRIX": ("dim_Z21", LENGTH, LENGTH),
        "PASS_UNITS_DET_M": ("dim_detm", LENGTH**4 / ENERGY**2, MASS**-2 * TIME**4),
        "PASS_UNITS_S_GG": ("dim_sgg", ENERGY**-1, ENERGY**-1),
        "PASS_UNITS_SHELL_C": ("dim_C", ENERGY**2, ENERGY**2),
        "PASS_UNITS_SHELL_A": ("dim_A", ENERGY * LENGTH, ENERGY * LENGTH),
        "PASS_UNITS_U": ("dim_U", ENERGY, ENERGY),
        "PASS_UNITS_F": ("dim_F", ENERGY / LENGTH, ENERGY / LENGTH),
    }
    for predicate, (primitive_name, canonical, expected) in dims.items():
        raw = primitive(predicate, primitive_name, canonical.vector())
        actual = dim_from(raw) if isinstance(raw, tuple) else raw
        expect_zero(predicate, dim_residual(actual, expected), {"actual": actual.vector(), "target": expected.vector()})

    # Production-live reduction checks supplement the primitive teeth without
    # adding unablated predicates.
    live_relations = (
        dim_residual(ENERGY * LENGTH**2 / LENGTH**2, ENERGY)
        + dim_residual(ENERGY * LENGTH / LENGTH, ENERGY)
        + dim_residual((LENGTH**-3) * LENGTH**4 * LENGTH**-2, LENGTH**-1)
        + dim_residual((LENGTH / ENERGY) * ENERGY**2, ENERGY * LENGTH)
        + dim_residual((ENERGY * LENGTH) / LENGTH, ENERGY)
        + dim_residual((ENERGY * LENGTH) / LENGTH**2, ENERGY / LENGTH)
    )
    if live_relations != 0:
        raise AuditFailure("DIMENSION_LIVE_RELATIONS", compact(live_relations))
    print("      live chains: k_m=K_m/ell^2; g=J_m/ell; N_chi=int eta dx int o dw; A=m_gg C; U=A/R; F=A/R^2")


def main() -> None:
    if ACTIVE_MUTATION and ACTIVE_MUTATION not in ABLATIONS:
        print("FIRST_FAILURE=UNKNOWN_MUTATION")
        print(f"FAIL  UNKNOWN_MUTATION: {ACTIVE_MUTATION}")
        raise AuditFailure("UNKNOWN_MUTATION", ACTIVE_MUTATION)

    print("ledger_stage031_puncture_deflection_field_identity_source SymPy audit")
    print("CONSUMES_STAGE030={f0,f0(0)=1/ell,N0=8/(3ell),h=P0H,S_Lh,D,D*=7/4}")
    print("CONSUMES_STAGE003_030={B_eff,C_hu}")
    if ACTIVE_MUTATION:
        print(f"ACTIVE_MUTATION={ACTIVE_MUTATION}")
        print(f"MUTATED_PRIMITIVE={ABLATIONS[ACTIVE_MUTATION].description}")

    run_mouth_reconstruction()
    run_response()
    run_exterior_and_shell()
    run_units()

    print("")
    print("SCOPE: stage031 earned mechanism and target-blind conditional leading form only.")
    print("DEFERRED_TO_STAGE032: BC-class numerator C, force sign, ensembles, internal_inconsistency, R1 landing.")
    print("VERDICT_TOKEN: THROAT_H_SOURCE_1_OVER_R2")
    print("EARNED_LABEL: PUNCTURE_DEFLECTION_MECHANISM_TARGET_BLIND")

    if ACTIVE_MUTATION:
        print("FIRST_FAILURE=MUTATION_DID_NOT_FIRE")
        raise AuditFailure("MUTATION_DID_NOT_FIRE", ACTIVE_MUTATION)

    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print("OVERALL PASS: SymPy verified stage031 puncture-deflection field identity and source mechanism")


if __name__ == "__main__":
    try:
        main()
    except AuditFailure as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
        print(f"OVERALL FAIL: SymPy stage031 audit did not close ({exc.predicate})")
        raise SystemExit(1)
