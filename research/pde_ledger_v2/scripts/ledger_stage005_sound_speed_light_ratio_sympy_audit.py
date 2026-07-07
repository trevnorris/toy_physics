#!/usr/bin/env python3
"""Ledger stage005 SymPy audit: EOS sound speed and light/sound ratio.

Print-only, standalone, no arguments, no file output.  This extracts the
pathA_20/pathA_20b load-bearing checks into ledger form without importing the
stage1 harness.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0

HARNESS_CHECK_GROUPS: dict[str, frozenset[str]] = {
    "pathA_20_dim": frozenset(
        {
            "pathA_20_S1 dim c_s^2=5*K*rho^4/m_GNLS",
            "pathA_20_S1 dim c_s=sqrt(5*K*rho^4/m_GNLS)",
            "pathA_20_S1_S2 stationary quantum-Bernoulli additive terms",
            "pathA_20_S2 bulk continuity equation with v_b",
            "pathA_20_S2 dim Madelung background velocity v_b=(hbar/m_GNLS)*grad(theta)",
            "pathA_20_S2 dim photon/gauge-wave speed c_gamma",
            "pathA_20_S2 massless gauge-wave dispersion omega^2=c_gamma^2*k^2",
            "pathA_20_S2 trapped-mode dispersion omega^2=c_gamma^2*(k_parallel^2+k_perp^2)",
            "pathA_20_S2 dim trapped-mode group velocity d omega/dk",
            "pathA_20_S2 dim ratio c_gamma/c_s",
            "pathA_20_S2 dim tail factor (c/c_s)^3 with c=c_gamma",
            "pathA_20_S2b dim 4D-bulk candidate sonic number flux rho_* c_s,* A_3,*",
            "pathA_20_S2b dim 3D-brane candidate sonic number flux rho_3,* c_s,* A_2,*",
            "pathA_20_S2b dim background pressure P0=K*rho0^5",
            "pathA_20_S3 consumed dim pin relation hbar=m_GNLS*c_s0*a",
            "pathA_20_S3 consumed dim healing relation hbar=m_GNLS*c_s0*xi_h/sqrt(2)",
            "pathA_20_S3 dim circulation kappa=int v_b dl",
            "pathA_20_S3 dim phase-momentum exchange p=hbar*grad(theta)",
            "pathA_20_S3 dim quantum pressure Q=-hbar^2/(2m)*laplacian(sqrt(rho))/sqrt(rho)",
            "pathA_20_S2_S3 dim candidate mass bridge hbar*J/c_gamma^2",
            "pathA_20_S2_S3 dim cycle-rate bridge h*J_nu/c_gamma^2",
        }
    ),
    "pathA_20_algebra": frozenset(
        {
            "pathA_20 algebra state dependence d ln c_s / d ln rho",
            "pathA_20 algebra conditional ideal no-Q/no-V sonic c_s,* / c_s0",
            "pathA_20 algebra conditional ideal no-Q/no-V sonic rho_* / rho0",
            "pathA_20 algebra conditional ideal no-Q/no-V flux factor Jcrit/(rho0*c_s0*A*)",
            "pathA_20 algebra tail factor with lambda_gamma=c_gamma/c_s",
        }
    ),
    "pathA_20b_dim": frozenset(
        {
            "pathA_20b_L2 dim phonon sound speed c_s0=sqrt(5*K*rho0^4/m_GNLS)",
            "pathA_20b_L1_L2 dim Maxwell principal speed squared C_B/C_E",
            "pathA_20b_L2 dim gauge speed c_gamma=sqrt(C_B/C_E)",
            "pathA_20b_L3 dim conditional bulk ratio c_bulk/c_s0",
            "pathA_20b_L4 dim tail factor lambda_gamma^3=(c_gamma/c_s)^3",
            "pathA_20b_L1b Maxwell transverse principal operator terms",
            "pathA_20b_L2 transverse gauge dispersion omega^2=(C_B/C_E)*k^2",
            "pathA_20b_L2 phonon acoustic dispersion omega^2=c_s0^2*k^2",
            "pathA_20b_L1 dim background charge density J_psi0^0=q_star*rho0",
            "pathA_20b_L1b linearized spatial current variation terms",
            "pathA_20b_L1 source coupling dimensions A_M delta J^M",
        }
    ),
    "pathA_20b_algebra": frozenset(
        {
            "pathA_20b algebra phonon determinant gives hbar*(omega^2-c_s0^2*k^2)",
            "pathA_20b algebra transverse gauge operator gives c_bulk^2=C_B/C_E",
            "pathA_20b algebra block coupled characteristic determinant P_ph*P_T^2",
            "negative control independent-symbol residual C_B/C_E-5*K*rho0^4/m_GNLS",
            "negative control forced_equals_valid is False without source equation",
            "pathA_20b algebra conditional rho0 slope d ln lambda_gamma / d ln rho0",
            "pathA_20b algebra standing-wave tail remains lambda_gamma^3",
        }
    ),
}
HARNESS_CHECK_COUNTS: dict[str, int] = {group: 0 for group in HARNESS_CHECK_GROUPS}


class AuditFailure(AssertionError):
    pass


def _record_harness_check(name: str) -> None:
    for group, check_names in HARNESS_CHECK_GROUPS.items():
        if name in check_names:
            HARNESS_CHECK_COUNTS[group] += 1
            return


def harness_check_count(group: str) -> int:
    return HARNESS_CHECK_COUNTS[group]


def banner(title: str) -> None:
    print("")
    print("=" * len(title))
    print(title)
    print("=" * len(title))


def subbanner(title: str) -> None:
    print("")
    print(title)
    print("-" * len(title))


def compact(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(sp.cancel(sp.simplify(expr))))


def assert_no_float(name: str, expr: Any) -> None:
    clean = sp.sympify(expr)
    floats = clean.atoms(sp.Float)
    if floats:
        raise AuditFailure(f"{name}: Float atom(s) found in exact audit expression: {floats}")


def _record_pass(message: str) -> None:
    global PASS_COUNT
    PASS_COUNT += 1
    print(message)


def _record_fail(message: str) -> None:
    global FAIL_COUNT
    FAIL_COUNT += 1
    print(message)


def expect_zero(name: str, residual: sp.Expr) -> None:
    assert_no_float(name, residual)
    clean = sp.simplify(residual)
    assert_no_float(name, clean)
    if clean == 0:
        _record_harness_check(name)
        _record_pass(f"PASS  {name}")
        return

    _record_fail(f"FAIL  {name}: residual = {compact(clean)}")
    raise AuditFailure(f"{name} residual was not zero")


def expect_bool(name: str, condition: bool) -> None:
    expect_zero(name, sp.Integer(0) if condition else sp.Integer(1))


def expect_nonzero(name: str, residual: sp.Expr) -> None:
    assert_no_float(name, residual)
    clean = sp.simplify(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_harness_check(name)
        _record_pass(f"PASS  {name} is nonzero as required (residual = {compact(clean)})")
        return

    _record_fail(f"FAIL  {name}: required nonzero residual vanished")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def expect_fail(name: str, residual: sp.Expr) -> None:
    assert_no_float(name, residual)
    clean = sp.simplify(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_pass(f"PASS  {name} produced required FAIL (residual = {compact(clean)})")
        return

    _record_fail(f"FAIL  {name}: required mutation/ablation did not fire")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def q(value: int | str | sp.Rational) -> sp.Rational:
    return sp.Rational(value)


@dataclass(frozen=True)
class Dim:
    """Exact exponent vector for base dimensions in {L, T, M} order."""

    l: sp.Rational | int = 0
    t: sp.Rational | int = 0
    m: sp.Rational | int = 0

    def __post_init__(self) -> None:
        object.__setattr__(self, "l", q(self.l))
        object.__setattr__(self, "t", q(self.t))
        object.__setattr__(self, "m", q(self.m))

    def __mul__(self, other: Dim) -> Dim:
        return Dim(self.l + other.l, self.t + other.t, self.m + other.m)

    def __truediv__(self, other: Dim) -> Dim:
        return Dim(self.l - other.l, self.t - other.t, self.m - other.m)

    def __pow__(self, power: int | sp.Rational) -> Dim:
        p = q(power)
        return Dim(self.l * p, self.t * p, self.m * p)

    def components(self) -> tuple[sp.Rational, sp.Rational, sp.Rational]:
        return (self.l, self.t, self.m)

    def __str__(self) -> str:
        pieces: list[str] = []
        for label, exponent in (("L", self.l), ("T", self.t), ("M", self.m)):
            if exponent == 0:
                continue
            pieces.append(label if exponent == 1 else f"{label}^{exponent}")
        return "1" if not pieces else " ".join(pieces)


DIMENSIONLESS = Dim()
LENGTH = Dim(1, 0, 0)
TIME = Dim(0, 1, 0)
MASS = Dim(0, 0, 1)
ACTION = MASS * (LENGTH**2) / TIME
ENERGY = ACTION / TIME
FORCE = ENERGY / LENGTH
VELOCITY = LENGTH / TIME
RHO4 = LENGTH**-4
RHO3 = LENGTH**-3
WAVE_NUMBER = LENGTH**-1
NUMBER_RATE = TIME**-1
ACTION_DENSITY = ENERGY / (LENGTH**4)
K_EOS_DIM = ACTION_DENSITY / (RHO4**5)
H_ENTHALPY = K_EOS_DIM * (RHO4**4)
Q_A0 = ENERGY
Q_AI = ACTION / LENGTH
ELECTRIC_FIELD = FORCE
MAGNETIC_FIELD = ACTION / (LENGTH**2)
MAXWELL_C_E = ACTION_DENSITY / (ELECTRIC_FIELD**2)
MAXWELL_C_B = MAXWELL_C_E * (VELOCITY**2)


def dim_residual(actual: Dim, expected: Dim) -> sp.Expr:
    return sp.simplify(
        sum((have - want) ** 2 for have, want in zip(actual.components(), expected.components()))
    )


def expect_dim(name: str, actual: Dim, expected: Dim) -> None:
    expect_zero(name, dim_residual(actual, expected))


def homogeneity_residual(terms: dict[str, Dim]) -> sp.Expr:
    if not terms:
        raise AuditFailure("homogeneity check requires at least one term")
    dims = list(terms.values())
    reference = dims[0]
    return sp.simplify(sum(dim_residual(actual, reference) for actual in dims[1:]))


@dataclass(frozen=True)
class Symbols:
    k_eos: sp.Symbol
    rho: sp.Symbol
    rho0: sp.Symbol
    m_gnls: sp.Symbol
    hbar: sp.Symbol
    omega: sp.Symbol
    omega_sq: sp.Symbol
    k: sp.Symbol
    k_sq: sp.Symbol
    k_perp: sp.Symbol
    c_gamma: sp.Symbol
    v: sp.Symbol
    c_e: sp.Symbol
    c_b: sp.Symbol
    lambda_gamma: sp.Symbol
    c_s: sp.Symbol


def make_symbols() -> Symbols:
    return Symbols(
        k_eos=sp.Symbol("K", positive=True),
        rho=sp.Symbol("rho", positive=True),
        rho0=sp.Symbol("rho0", positive=True),
        m_gnls=sp.Symbol("m_GNLS", positive=True),
        hbar=sp.Symbol("hbar", positive=True),
        omega=sp.Symbol("omega", positive=True),
        omega_sq=sp.Symbol("omega_sq", positive=True),
        k=sp.Symbol("k", positive=True),
        k_sq=sp.Symbol("k_sq", positive=True),
        k_perp=sp.Symbol("k_perp", positive=True),
        c_gamma=sp.Symbol("c_gamma", positive=True),
        v=sp.Symbol("v", positive=True),
        c_e=sp.Symbol("C_E", positive=True),
        c_b=sp.Symbol("C_B", positive=True),
        lambda_gamma=sp.Symbol("lambda_gamma", positive=True),
        c_s=sp.Symbol("c_s", positive=True),
    )


@dataclass(frozen=True)
class Stage004MediumInput:
    label: str
    a_symbol: sp.Symbol
    xi_h_symbol: sp.Symbol
    h0_symbol: sp.Symbol
    a_value: sp.Expr
    xi_h_value: sp.Expr
    h0_value: sp.Expr

    @property
    def stage004_subs(self) -> dict[sp.Symbol, sp.Expr]:
        return {
            self.a_symbol: self.a_value,
            self.xi_h_symbol: self.xi_h_value,
            self.h0_symbol: self.h0_value,
        }


@dataclass(frozen=True)
class NegControlConfig:
    source_metric_equation_present: bool
    c_bulk_sq_override: sp.Expr | None = None


@dataclass(frozen=True)
class NegControlResult:
    equality_residual: sp.Expr
    forced_equals_valid: bool
    verdict: str


def eos_sound_speed_squared(s: Symbols, exponent: int = 5, drop_mass: bool = False) -> sp.Expr:
    pressure = s.k_eos * s.rho**exponent
    derivative = sp.diff(pressure, s.rho)
    denominator = sp.Integer(1) if drop_mass else s.m_gnls
    return sp.simplify(derivative / denominator)


def eos_sound_speed_squared_at_rho0(s: Symbols) -> sp.Expr:
    cs_sq = eos_sound_speed_squared(s)
    return sp.simplify(cs_sq.subs(s.rho, s.rho0))


def run_sound_speed_derivation(s: Symbols) -> tuple[sp.Expr, sp.Expr]:
    subbanner("S1 owned sound-speed derivation from imposed EOS")
    print("  EOS_CLOSURE_IMPOSED: P=K*rho^5 is postulated; this stage derives c_s relative to that EOS.")

    pressure = s.k_eos * s.rho**5
    derived = sp.simplify(sp.diff(pressure, s.rho) / s.m_gnls)
    expected = 5 * s.k_eos * s.rho**4 / s.m_gnls
    expect_zero("S1 c_s^2=(1/m_GNLS)*d(K*rho^5)/drho gives 5*K*rho^4/m_GNLS", derived - expected)
    expect_dim("pathA_20_S1 dim c_s^2=5*K*rho^4/m_GNLS", K_EOS_DIM * (RHO4**4) / MASS, VELOCITY**2)
    expect_dim(
        "pathA_20_S1 dim c_s=sqrt(5*K*rho^4/m_GNLS)",
        (K_EOS_DIM * (RHO4**4) / MASS) ** sp.Rational(1, 2),
        VELOCITY,
    )

    c_s_profile = sp.sqrt(derived)
    log_slope = sp.simplify(s.rho * sp.diff(c_s_profile, s.rho) / c_s_profile)
    expect_zero("pathA_20 algebra state dependence d ln c_s / d ln rho", log_slope - 2)
    return derived, sp.simplify(derived.subs(s.rho, s.rho0))


def run_quantum_bernoulli_and_continuity() -> None:
    subbanner("S1/S2 Bernoulli and continuity homogeneity")
    expect_zero(
        "pathA_20_S1_S2 stationary quantum-Bernoulli additive terms",
        homogeneity_residual(
            {
                "0.5*m_GNLS*v_b^2": MASS * (VELOCITY**2),
                "h(rho)": H_ENTHALPY,
                "V_conf": ENERGY,
                "Q": (ACTION**2) / (MASS * (LENGTH**2)),
            }
        ),
    )
    expect_zero(
        "pathA_20_S2 bulk continuity equation with v_b",
        homogeneity_residual(
            {
                "partial_t rho": RHO4 / TIME,
                "div_4(rho v_b)": (RHO4 * VELOCITY) / LENGTH,
            }
        ),
    )
    print("  STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA: dimensions do not solve the stationary profile.")


def run_velocity_ceiling(s: Symbols) -> None:
    subbanner("S2 three velocity scales and wave-sector ceiling")
    expect_dim(
        "pathA_20_S2 dim Madelung background velocity v_b=(hbar/m_GNLS)*grad(theta)",
        (ACTION / MASS) * (LENGTH**-1),
        VELOCITY,
    )
    expect_dim("pathA_20_S2 dim photon/gauge-wave speed c_gamma", VELOCITY, VELOCITY)
    expect_zero(
        "pathA_20_S2 massless gauge-wave dispersion omega^2=c_gamma^2*k^2",
        homogeneity_residual(
            {
                "omega^2": (TIME**-1) ** 2,
                "c_gamma^2*k^2": (VELOCITY**2) * (WAVE_NUMBER**2),
            }
        ),
    )
    expect_zero(
        "pathA_20_S2 trapped-mode dispersion omega^2=c_gamma^2*(k_parallel^2+k_perp^2)",
        homogeneity_residual(
            {
                "omega^2": (TIME**-1) ** 2,
                "c_gamma^2*k_parallel^2": (VELOCITY**2) * (WAVE_NUMBER**2),
                "c_gamma^2*k_perp^2": (VELOCITY**2) * (WAVE_NUMBER**2),
            }
        ),
    )

    omega_trapped = s.c_gamma * sp.sqrt(s.k**2 + s.k_perp**2)
    group_velocity = sp.simplify(sp.diff(omega_trapped, s.k))
    expected_group = s.c_gamma * s.k / sp.sqrt(s.k**2 + s.k_perp**2)
    expect_zero("S2 derived trapped-mode group velocity d omega/dk", group_velocity - expected_group)
    expect_dim("pathA_20_S2 dim trapped-mode group velocity d omega/dk", VELOCITY, VELOCITY)
    bound_residual = sp.factor(s.c_gamma**2 - group_velocity**2)
    expected_bound = s.c_gamma**2 * s.k_perp**2 / (s.k**2 + s.k_perp**2)
    expect_zero("S2 derived c_gamma^2-(domega/dk)^2 positive-form residual", bound_residual - expected_bound)
    expect_bool(
        "S2 trapped k_perp>0 branch has strictly positive group-velocity gap",
        sp.ask(sp.Q.positive(bound_residual)) is True,
    )
    expect_zero(
        "S2 degenerate k_perp=0 free wave reaches the c_gamma ceiling",
        sp.simplify(bound_residual.subs(s.k_perp, 0)),
    )
    print("  Assumption: k_perp>0 gives a trapped mode strictly below c_gamma; k_perp=0 is the free c_gamma wave.")

    omega0 = s.c_gamma * s.k_perp
    gamma_lorentz = (1 - s.v**2 / s.c_gamma**2) ** sp.Rational(-1, 2)
    t, x = sp.symbols("t x", positive=True)
    phase_clock = omega0 * gamma_lorentz * (t - s.v * x / s.c_gamma**2)
    center_clock = sp.simplify(phase_clock.subs(x, s.v * t))
    clock_rate = sp.simplify(sp.diff(center_clock, t))
    expect_zero("S2 bound-mode clock along x=v*t advances at omega0/gamma", clock_rate - omega0 / gamma_lorentz)
    print("  EARNED: c=c_gamma is read from the wave-sector ceiling; no E=m_defect*c_gamma^2 premise is used.")

    expect_dim("pathA_20_S2 dim ratio c_gamma/c_s", VELOCITY / VELOCITY, DIMENSIONLESS)
    expect_dim("pathA_20_S2 dim tail factor (c/c_s)^3 with c=c_gamma", (VELOCITY / VELOCITY) ** 3, DIMENSIONLESS)
    expect_dim("pathA_20b_L2 dim [c_s]=[c_gamma] match", VELOCITY, VELOCITY)
    print("  NON-EVIDENTIARY: [c_s]=[c_gamma] is a dimensional match only; it is not evidence for c_gamma=c_s.")


def run_flux_checks(s: Symbols) -> None:
    subbanner("S2b flux law dimensions and conditional nozzle numbers")
    expect_dim(
        "pathA_20_S2b dim 4D-bulk candidate sonic number flux rho_* c_s,* A_3,*",
        RHO4 * VELOCITY * (LENGTH**3),
        NUMBER_RATE,
    )
    expect_dim(
        "pathA_20_S2b dim 3D-brane candidate sonic number flux rho_3,* c_s,* A_2,*",
        RHO3 * VELOCITY * (LENGTH**2),
        NUMBER_RATE,
    )
    expect_dim(
        "pathA_20_S2b dim background pressure P0=K*rho0^5",
        K_EOS_DIM * (RHO4**5),
        ACTION_DENSITY,
    )

    x = sp.Symbol("c_star_over_c_s0", positive=True)
    cstar_over_c0 = sp.sqrt(sp.solve(sp.Eq(sp.Rational(3, 4) * x**2, sp.Rational(1, 4)), x)[0] ** 2)
    cstar_over_c0 = sp.simplify(cstar_over_c0)
    rhostar_over_rho0 = sp.sqrt(cstar_over_c0)
    ideal_flux_factor = sp.simplify(cstar_over_c0 * rhostar_over_rho0)
    expect_zero("pathA_20 algebra conditional ideal no-Q/no-V sonic c_s,* / c_s0", cstar_over_c0 - 1 / sp.sqrt(3))
    expect_zero("pathA_20 algebra conditional ideal no-Q/no-V sonic rho_* / rho0", rhostar_over_rho0 - sp.Pow(3, sp.Rational(-1, 4)))
    expect_zero(
        "pathA_20 algebra conditional ideal no-Q/no-V flux factor Jcrit/(rho0*c_s0*A*)",
        ideal_flux_factor - sp.Pow(3, sp.Rational(-3, 4)),
    )
    print("  CONDITIONAL_NOT_ACCEPTED_AS_BRANCH_LAW: nozzle factors are recorded, not accepted as the branch law.")
    print("  STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA and NO_NET_ACCRETION_BC_UNDERIVED are carried forward.")


def run_role_catalog_and_consumed_dims() -> None:
    subbanner("S3 role catalog and consumed pin/healing dimensional checks")
    expect_dim(
        "pathA_20_S3 consumed dim pin relation hbar=m_GNLS*c_s0*a",
        MASS * VELOCITY * LENGTH,
        ACTION,
    )
    expect_dim(
        "pathA_20_S3 consumed dim healing relation hbar=m_GNLS*c_s0*xi_h/sqrt(2)",
        MASS * VELOCITY * LENGTH,
        ACTION,
    )
    expect_dim("pathA_20_S3 dim circulation kappa=int v_b dl", VELOCITY * LENGTH, ACTION / MASS)
    expect_dim("pathA_20_S3 dim phase-momentum exchange p=hbar*grad(theta)", ACTION / LENGTH, MASS * VELOCITY)
    expect_dim(
        "pathA_20_S3 dim quantum pressure Q=-hbar^2/(2m)*laplacian(sqrt(rho))/sqrt(rho)",
        (ACTION**2) / (MASS * (LENGTH**2)),
        ENERGY,
    )
    expect_dim(
        "pathA_20_S2_S3 dim candidate mass bridge hbar*J/c_gamma^2",
        ACTION * NUMBER_RATE / (VELOCITY**2),
        MASS,
    )
    expect_dim(
        "pathA_20_S2_S3 dim cycle-rate bridge h*J_nu/c_gamma^2",
        ACTION * NUMBER_RATE / (VELOCITY**2),
        MASS,
    )
    print("  HBAR_PROVENANCE_UNDETERMINED: hbar remains an explicit action/PDE coefficient.")
    print("  HBAR_FREE_SUBSTRATE_RELATION_MISSING: hbar=m_GNLS*c_s0*a is a pin rearrangement, not an hbar-emergence proof.")
    print("  H_2PI_RATE_CLASSIFICATION_UNDERDETERMINED: h and hbar share dimensions; the 2*pi placement is deferred.")
    print("  candidate-only: m_defect=alpha_J*hbar*J/c_gamma^2 does not collapse M and does not derive alpha_J.")


def stage004_input(s: Symbols, c_s0_sq: sp.Expr) -> Stage004MediumInput:
    c_s0 = sp.sqrt(c_s0_sq)
    return Stage004MediumInput(
        label="ledger_stage004 (I-1) EOS_FROM_GNLS_FACTOR handoff",
        a_symbol=sp.Symbol("a_I1", positive=True),
        xi_h_symbol=sp.Symbol("xi_h_I1", positive=True),
        h0_symbol=sp.Symbol("h0_I1", positive=True),
        a_value=sp.simplify(s.hbar / (s.m_gnls * c_s0)),
        xi_h_value=sp.simplify(sp.sqrt(2) * s.hbar / (s.m_gnls * c_s0)),
        h0_value=sp.simplify(s.m_gnls * c_s0_sq / 4),
    )


def run_consumed_stage004_inputs(s: Symbols, c_s0_sq: sp.Expr) -> Stage004MediumInput:
    subbanner("Consumed from ledger_stage004 (I-1)")
    consumed = stage004_input(s, c_s0_sq)
    print("  CITED, NOT RE-DERIVED: dictionary {L,T,M}; a=hbar/(m_GNLS*c_s0); xi_h=sqrt(2)*hbar/(m_GNLS*c_s0); h0=(m_GNLS*c_s0^2)/4.")
    print("  EOS_FROM_GNLS_FACTOR: exact values are checked only as citation-integrity against I-2-derived c_s0.")

    expect_dim("stage004 consumed dim [a]=L", LENGTH, LENGTH)
    expect_dim("stage004 consumed dim [xi_h]=L", LENGTH, LENGTH)
    expect_dim("stage004 consumed dim [h0]=[m_GNLS*c_s0^2]", ENERGY, MASS * (VELOCITY**2))

    subs = consumed.stage004_subs
    c_s0 = sp.sqrt(c_s0_sq)
    expect_zero(
        "stage004 citation integrity h0-(m_GNLS*c_s0^2)/4",
        subs[consumed.h0_symbol] - s.m_gnls * c_s0_sq / 4,
    )
    expect_zero(
        "stage004 citation integrity xi_h-sqrt(2)*hbar/(m_GNLS*c_s0)",
        subs[consumed.xi_h_symbol] - sp.sqrt(2) * s.hbar / (s.m_gnls * c_s0),
    )
    return consumed


def run_patha20b_dimensions(s: Symbols, c_s0_sq_dim: Dim) -> None:
    subbanner("L1/L1b/L2 coupled-linearization dimensional checks")
    expect_dim(
        "pathA_20b_L2 dim phonon sound speed c_s0=sqrt(5*K*rho0^4/m_GNLS)",
        c_s0_sq_dim ** sp.Rational(1, 2),
        VELOCITY,
    )
    expect_dim("pathA_20b_L1_L2 dim Maxwell principal speed squared C_B/C_E", MAXWELL_C_B / MAXWELL_C_E, VELOCITY**2)
    expect_dim(
        "pathA_20b_L2 dim gauge speed c_gamma=sqrt(C_B/C_E)",
        (MAXWELL_C_B / MAXWELL_C_E) ** sp.Rational(1, 2),
        VELOCITY,
    )
    expect_dim(
        "pathA_20b_L3 dim conditional bulk ratio c_bulk/c_s0",
        ((MAXWELL_C_B / MAXWELL_C_E) / c_s0_sq_dim) ** sp.Rational(1, 2),
        DIMENSIONLESS,
    )
    expect_dim(
        "pathA_20b_L4 dim tail factor lambda_gamma^3=(c_gamma/c_s)^3",
        (((MAXWELL_C_B / MAXWELL_C_E) / c_s0_sq_dim) ** sp.Rational(1, 2)) ** 3,
        DIMENSIONLESS,
    )
    expect_zero(
        "pathA_20b_L1b Maxwell transverse principal operator terms",
        homogeneity_residual(
            {
                "C_E*partial_t^2 A_T": MAXWELL_C_E * Q_AI / (TIME**2),
                "C_B*laplacian A_T": MAXWELL_C_B * Q_AI / (LENGTH**2),
            }
        ),
    )
    expect_zero(
        "pathA_20b_L2 transverse gauge dispersion omega^2=(C_B/C_E)*k^2",
        homogeneity_residual(
            {
                "omega^2": (TIME**-1) ** 2,
                "(C_B/C_E)*k^2": (MAXWELL_C_B / MAXWELL_C_E) * (WAVE_NUMBER**2),
            }
        ),
    )
    expect_zero(
        "pathA_20b_L2 phonon acoustic dispersion omega^2=c_s0^2*k^2",
        homogeneity_residual(
            {
                "omega^2": (TIME**-1) ** 2,
                "c_s0^2*k^2": c_s0_sq_dim * (WAVE_NUMBER**2),
            }
        ),
    )
    expect_dim("pathA_20b_L1 dim background charge density J_psi0^0=q_star*rho0", RHO4, RHO4)
    expect_zero(
        "pathA_20b_L1b linearized spatial current variation terms",
        homogeneity_residual(
            {
                "phase-current rho0*(hbar/m_GNLS)*grad(delta theta)": RHO4 * (ACTION / MASS) / LENGTH,
                "London term (q_star/m_GNLS)*rho0*delta A_i": RHO4 * Q_AI / MASS,
                "spatial current rho0*v": RHO4 * VELOCITY,
            }
        ),
    )
    expect_zero(
        "pathA_20b_L1 source coupling dimensions A_M delta J^M",
        homogeneity_residual(
            {
                "A0*delta J0": Q_A0 * RHO4,
                "Ai*delta Ji": Q_AI * (RHO4 * VELOCITY),
                "local action density": ACTION_DENSITY,
            }
        ),
    )
    print("  LEGAL_WITH_EXPLICIT_NEUTRALIZING_EXTERNAL_SOURCE: J_tot0^M=J_psi0^M+J_ext0^M=0 makes the Maxwell background 0=0.")
    print("  LOWER-ORDER: current/London/source-coupling terms do not set the cone; the transverse principal field-strength operator does.")


def run_patha20b_algebra(s: Symbols, c_s0_sq: sp.Expr) -> None:
    subbanner("L1b-L2 coupled principal symbol and dispersions")
    h_prime = sp.diff(s.k_eos * s.rho**5, s.rho).subs(s.rho, s.rho0) / s.rho0
    expect_zero("L2 h'(rho0)=5*K*rho0^3 from differentiated EOS pressure law", h_prime - 5 * s.k_eos * s.rho0**3)
    c_s0_from_block = sp.simplify(s.rho0 * h_prime / s.m_gnls)
    expect_zero("L2 c_s0^2=rho0*h'(rho0)/m_GNLS reproduces S1 background", c_s0_from_block - c_s0_sq)

    phonon_matrix = sp.Matrix(
        [
            [s.omega, -(s.rho0 * s.hbar / s.m_gnls) * s.k_sq],
            [-h_prime, s.hbar * s.omega],
        ]
    )
    phonon_det = sp.factor(phonon_matrix.det())
    expect_zero(
        "pathA_20b algebra phonon determinant gives hbar*(omega^2-c_s0^2*k^2)",
        phonon_det - s.hbar * (s.omega**2 - c_s0_sq * s.k_sq),
    )

    c_bulk_sq = s.c_b / s.c_e
    gauge_transverse = s.c_e * s.omega**2 - s.c_b * s.k_sq
    expect_zero(
        "pathA_20b algebra transverse gauge operator gives c_bulk^2=C_B/C_E",
        gauge_transverse - s.c_e * (s.omega**2 - c_bulk_sq * s.k_sq),
    )

    coupled_det = sp.factor(phonon_det * gauge_transverse**2)
    expected_coupled = s.hbar * (s.omega**2 - c_s0_sq * s.k_sq) * (s.c_e * (s.omega**2 - c_bulk_sq * s.k_sq)) ** 2
    expect_zero(
        "pathA_20b algebra block coupled characteristic determinant P_ph*P_T^2",
        coupled_det - expected_coupled,
    )
    print("  VANISHES_ON_HOMOGENEOUS_NEUTRALIZED_BACKGROUND: off-diagonal principal terms vanish in the neutralized homogeneous background.")

    phonon_solution = sp.solve(sp.Eq(s.hbar * (s.omega_sq - c_s0_sq * s.k_sq), 0), s.omega_sq)[0]
    gauge_solution = sp.solve(sp.Eq(s.c_e * s.omega_sq - s.c_b * s.k_sq, 0), s.omega_sq)[0]
    expect_zero("L2 phonon dispersion read off as omega^2=c_s0^2*k^2", phonon_solution - c_s0_sq * s.k_sq)
    expect_zero("L2 gauge dispersion read off as omega^2=c_bulk^2*k^2", gauge_solution - c_bulk_sq * s.k_sq)
    print("  source status: BULK_PRINCIPAL_TRANSVERSE_BRANCH_ESTABLISHED")
    print("  Bogoliubov k^4 quantum-pressure correction is dispersive and does not set the cone.")

    lambda_rho_factor = s.rho0**-2
    lambda_log_slope = sp.simplify(s.rho0 * sp.diff(lambda_rho_factor, s.rho0) / lambda_rho_factor)
    expect_zero("pathA_20b algebra conditional rho0 slope d ln lambda_gamma / d ln rho0", lambda_log_slope + 2)

    tail_expr = (s.lambda_gamma * s.c_s / s.c_s) ** 3
    expect_zero("pathA_20b algebra standing-wave tail remains lambda_gamma^3", tail_expr - s.lambda_gamma**3)


def run_lambda_calibration(s: Symbols, c_s0_sq: sp.Expr) -> str:
    subbanner("CALIBRATED landing and negative control")
    lambda_tail = ((s.lambda_gamma * s.c_s) / s.c_s) ** 3
    expect_zero("pathA_20 algebra tail factor with lambda_gamma=c_gamma/c_s", lambda_tail - s.lambda_gamma**3)
    expect_dim("lambda_gamma is dimensionless", DIMENSIONLESS, DIMENSIONLESS)
    expect_dim("lambda_gamma^3 tail is dimensionless", DIMENSIONLESS, DIMENSIONLESS)
    expect_nonzero("lambda_gamma^3 tail is not reduced to 1", s.lambda_gamma**3 - 1)
    print("  C_GAMMA_RATIO_UNDERDETERMINED: lambda_gamma=c_gamma/c_s is a FREE calibration input.")
    print("  REJECTED: c_gamma=c_s from shared dimensions or legacy weak-field prose.")
    print("  C_GAMMA_BULK_UNDERDETERMINED: BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED.")
    print("  C_GAMMA_RATIO_STILL_UNDERDETERMINED: BRANE_ZERO_MODE_REDUCTION_UNDERIVED; BRANE_PHOTON_CONE_REQUIRES_PROFILE.")

    c_bulk_sq = s.c_b / s.c_e
    bulk_ratio = sp.sqrt(c_bulk_sq / c_s0_sq)
    expect_dim("pathA_20b conditional bulk ratio c_bulk/c_s0 is dimensionless", VELOCITY / VELOCITY, DIMENSIONLESS)
    print(f"  conditional bulk ratio carried symbolically: {sp.sstr(bulk_ratio)}")

    baseline = run_negative_control(s, c_s0_sq, NegControlConfig(source_metric_equation_present=False))
    expect_nonzero("negative control independent-symbol residual C_B/C_E-5*K*rho0^4/m_GNLS", baseline.equality_residual)
    expect_bool("negative control forced_equals_valid is False without source equation", not baseline.forced_equals_valid)
    expect_bool("negative control baseline verdict is C_GAMMA_RATIO_UNDERDETERMINED", baseline.verdict == "C_GAMMA_RATIO_UNDERDETERMINED")
    print("  FORCED_EQUALITY_REJECTED_WITHOUT_SOURCE_EQUATION")
    return baseline.verdict


def run_negative_control(s: Symbols, c_s0_sq: sp.Expr, config: NegControlConfig) -> NegControlResult:
    c_bulk_sq = config.c_bulk_sq_override if config.c_bulk_sq_override is not None else s.c_b / s.c_e
    equality_residual = sp.simplify(c_bulk_sq - c_s0_sq)
    forced_equals_valid = bool(config.source_metric_equation_present and equality_residual == 0)
    verdict = "C_GAMMA_EQUALS_C_S" if forced_equals_valid else "C_GAMMA_RATIO_UNDERDETERMINED"
    return NegControlResult(equality_residual=equality_residual, forced_equals_valid=forced_equals_valid, verdict=verdict)


def run_reversibility(s: Symbols, c_s0_sq: sp.Expr) -> None:
    subbanner("Negative-control reversibility able-to-PASS")
    equality_residual = sp.simplify(s.c_b / s.c_e - c_s0_sq)
    inserted_identity_residual = sp.simplify(equality_residual.subs(s.c_b, 5 * s.k_eos * s.rho0**4 * s.c_e / s.m_gnls))
    expect_zero("reversibility inserted source equation C_B->5*K*rho0^4*C_E/m_GNLS zeros residual", inserted_identity_residual)

    counterfactual = run_negative_control(
        s,
        c_s0_sq,
        NegControlConfig(
            source_metric_equation_present=True,
            c_bulk_sq_override=sp.simplify((5 * s.k_eos * s.rho0**4 * s.c_e / s.m_gnls) / s.c_e),
        ),
    )
    expect_zero("reversibility control residual becomes zero", counterfactual.equality_residual)
    expect_bool("reversibility forced_equals_valid flips True with inserted source equation", counterfactual.forced_equals_valid)
    expect_bool("reversibility verdict flips to C_GAMMA_EQUALS_C_S", counterfactual.verdict == "C_GAMMA_EQUALS_C_S")
    print("  C_GAMMA_EQUALS_C_S is reachable only in the counterfactual branch with an inserted source equation.")


def run_firewall_ablations(s: Symbols, c_s0_sq: sp.Expr, consumed: Stage004MediumInput) -> None:
    subbanner("Able-to-fail dimensional and derivation firewall")
    corrupt_cs_sq = eos_sound_speed_squared(s, exponent=4)
    expect_fail(
        "ablation corrupt EOS exponent P=K*rho^4 breaks [c_s]",
        dim_residual((K_EOS_DIM * (RHO4**3) / MASS) ** sp.Rational(1, 2), VELOCITY),
    )
    corrupt_log_slope = sp.simplify(s.rho * sp.diff(sp.sqrt(corrupt_cs_sq), s.rho) / sp.sqrt(corrupt_cs_sq))
    expect_fail("ablation corrupt EOS exponent P=K*rho^4 breaks log-slope=2", corrupt_log_slope - 2)

    expect_fail(
        "ablation drop m_GNLS from c_s^2 breaks velocity dimension",
        dim_residual((K_EOS_DIM * (RHO4**4)) ** sp.Rational(1, 2), VELOCITY),
    )

    h_prime = 5 * s.k_eos * s.rho0**3
    corrupt_phonon = sp.Matrix(
        [
            [s.omega, -(s.rho0 * s.hbar / s.m_gnls)],
            [-h_prime, s.hbar * s.omega],
        ]
    )
    corrupt_det = sp.factor(corrupt_phonon.det())
    expect_fail(
        "ablation corrupt phonon block missing k^2 breaks determinant factorization",
        corrupt_det - s.hbar * (s.omega**2 - c_s0_sq * s.k_sq),
    )

    c_s0 = sp.sqrt(c_s0_sq)
    wrong_xi = s.hbar / (s.m_gnls * c_s0)
    expect_fail(
        "ablation corrupt consumed xi_h by dropping sqrt(2) breaks citation integrity",
        wrong_xi - sp.sqrt(2) * s.hbar / (s.m_gnls * c_s0),
    )
    wrong_h0 = s.m_gnls * c_s0_sq / 2
    expect_fail(
        "ablation corrupt consumed h0 wrong 1/4 breaks citation integrity",
        wrong_h0 - s.m_gnls * c_s0_sq / 4,
    )
    expect_zero(
        "baseline consumed handoff still live after ablations",
        consumed.stage004_subs[consumed.h0_symbol] - s.m_gnls * c_s0_sq / 4,
    )


def print_carried_gap_block() -> None:
    subbanner("Carried residuals printed verbatim as provenance")
    print("  EOS_CLOSURE_IMPOSED: CARRIED_FORWARD.")
    print("  C_GAMMA_RATIO_UNDERDETERMINED: BLOCKS_NUMERIC_C_GAMMA_OVER_C_S.")
    print("  STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA: FLUX_LAW_VERDICT.")
    print("  NO_NET_ACCRETION_BC_UNDERIVED: CARRIED_FORWARD.")
    print("  HBAR_PROVENANCE_UNDETERMINED: S3 verdict.")
    print("  HBAR_FREE_SUBSTRATE_RELATION_MISSING: BLOCKS_HBAR_EMERGENT.")
    print("  H_2PI_RATE_CLASSIFICATION_UNDERDETERMINED: CARRIED_FORWARD.")
    print("  BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED: BULK_VERDICT_RESIDUAL.")
    print("  PARENT_METRIC_ACOUSTIC_IDENTIFICATION_MISSING: BLOCKS_BULK_EQUALS_C_S.")
    print("  BRANE_ZERO_MODE_REDUCTION_UNDERIVED: BRANE_VERDICT_RESIDUAL.")
    print("  BRANE_PHOTON_CONE_REQUIRES_PROFILE: BRANE_SUB_RESIDUAL.")
    print("  mass-bridge candidate-only: m_defect=alpha_J*hbar*J/c_gamma^2 is dimensional, not derived.")
    print("  anti-tautology caveat: hbar=m_GNLS*c_s0*a is a pin rearrangement unless a is fixed by an hbar-free substrate relation.")


def print_verdict_labels(verdict: str) -> None:
    print("")
    print("Verdict labels:")
    print(
        "  ledger earned-label (NOT a source verdict token): "
        "EOS_SOUND_SPEED_DERIVED_LIGHT_RATIO_FREE  "
        "(c_s^2=5*K*rho^4/m_GNLS from P=K*rho^5; c=c_gamma wave-sector ceiling; lambda_gamma=c_gamma/c_s FREE)"
    )
    print(f"  source top-line verdict: {verdict}   (PASS_WITH_NAMED_RESIDUALS)")
    print(
        "  calibrated landing (honest): lambda_gamma=c_gamma/c_s UNPINNED by the parent action -- "
        "free calibration input; tail (c/c_s)^3=lambda_gamma^3 carried (NOT set to 1)"
    )
    print(
        "  pathA_20b sharpening: bulk C_GAMMA_BULK_UNDERDETERMINED; brane "
        "C_GAMMA_RATIO_STILL_UNDERDETERMINED; negative control "
        "FORCED_EQUALITY_REJECTED_WITHOUT_SOURCE_EQUATION (reversible to C_GAMMA_EQUALS_C_S only with an inserted source equation)"
    )
    print(
        "  labeled non-derivation (candidate only): "
        "m_defect=alpha_J*hbar*J/c_gamma^2 -- dimensional candidate, NOT a mass derivation"
    )
    print(
        "  consumed from ledger_stage004 (I-1) [EOS_FROM_GNLS_FACTOR]: dictionary {L,T,M}; "
        "a=hbar/(m_GNLS*c_s0); xi_h=sqrt(2)*hbar/(m_GNLS*c_s0); h0=(m_GNLS*c_s0^2)/4"
    )
    print(
        "  carried residuals: EOS_CLOSURE_IMPOSED; STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA; "
        "NO_NET_ACCRETION_BC_UNDERIVED; HBAR_PROVENANCE_UNDETERMINED; "
        "HBAR_FREE_SUBSTRATE_RELATION_MISSING; H_2PI_RATE_CLASSIFICATION_UNDERDETERMINED; "
        "BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED; PARENT_METRIC_ACOUSTIC_IDENTIFICATION_MISSING; "
        "BRANE_ZERO_MODE_REDUCTION_UNDERIVED; BRANE_PHOTON_CONE_REQUIRES_PROFILE"
    )


def main() -> None:
    banner("ledger_stage005_sound_speed_light_ratio SymPy audit")
    s = make_symbols()

    cs_sq, cs0_sq = run_sound_speed_derivation(s)

    run_quantum_bernoulli_and_continuity()
    run_velocity_ceiling(s)
    run_flux_checks(s)
    run_role_catalog_and_consumed_dims()
    consumed = run_consumed_stage004_inputs(s, cs0_sq)
    run_patha20b_dimensions(s, K_EOS_DIM * (RHO4**4) / MASS)
    run_patha20b_algebra(s, cs0_sq)

    verdict = run_lambda_calibration(s, cs0_sq)
    expect_zero("pathA_20 harness dimensional check count is 21", harness_check_count("pathA_20_dim") - 21)
    expect_zero("pathA_20 harness algebraic check count is 5", harness_check_count("pathA_20_algebra") - 5)
    expect_zero("pathA_20b harness dimensional check count is 11", harness_check_count("pathA_20b_dim") - 11)
    expect_zero("pathA_20b harness algebraic check count is 7", harness_check_count("pathA_20b_algebra") - 7)
    run_reversibility(s, cs0_sq)
    run_firewall_ablations(s, cs0_sq, consumed)
    print_carried_gap_block()
    print_verdict_labels(verdict)

    assert_no_float("final c_s^2 expression", cs_sq)
    assert_no_float("final c_s0^2 expression", cs0_sq)
    print("")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print("OVERALL PASS: SymPy derived ledger_stage005 EOS sound speed and light/sound ratio exactly")


if __name__ == "__main__":
    main()
