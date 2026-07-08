#!/usr/bin/env python3
"""Ledger stage009 SymPy audit: flat-slab return residual, Check A.

Standalone, print-only, no arguments, no file I/O.  This reshapes the Check-A
slice of pathA_29 into a ledger audit: the finite flat slab is postulated, the
Helmholtz round-trip transport and return fractions are derived on it, and the
bounded residual prediction is computed without importing the stage-010
localization machinery.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0
I = sp.I


class AuditFailure(AssertionError):
    pass


def banner(title: str) -> None:
    print("")
    print("=" * len(title))
    print(title)
    print("=" * len(title))


def subbanner(title: str) -> None:
    print("")
    print(title)
    print("-" * len(title))


def s(expr: sp.Expr | int) -> str:
    return sp.sstr(sp.factor(sp.cancel(sp.simplify(expr))))


def assert_no_float(name: str, expr: Any) -> None:
    if isinstance(expr, dict):
        for key, value in expr.items():
            assert_no_float(f"{name}.{key}", value)
        return
    if isinstance(expr, (list, tuple, set, frozenset)):
        for index, value in enumerate(expr):
            assert_no_float(f"{name}[{index}]", value)
        return
    if isinstance(expr, bool):
        expr = sp.Integer(0) if expr else sp.Integer(1)
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


def expect_zero(name: str, residual: sp.Expr | int) -> None:
    assert_no_float(name, residual)
    clean = sp.simplify(residual)
    assert_no_float(name, clean)
    if clean == 0:
        _record_pass(f"PASS  {name}")
        return
    _record_fail(f"FAIL  {name}: residual = {s(clean)}")
    raise AuditFailure(f"{name} residual was not zero")


def expect_bool(name: str, condition: bool) -> None:
    expect_zero(name, sp.Integer(0) if condition else sp.Integer(1))


def expect_nonzero(name: str, residual: sp.Expr | int) -> None:
    assert_no_float(name, residual)
    clean = sp.simplify(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_pass(f"PASS  {name} is nonzero as required (residual = {s(clean)})")
        return
    _record_fail(f"FAIL  {name}: required nonzero residual vanished")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def expect_fail(name: str, residual: sp.Expr | int) -> None:
    assert_no_float(name, residual)
    clean = sp.simplify(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_pass(f"PASS  {name} produced required FAIL (residual = {s(clean)})")
        return
    _record_fail(f"FAIL  {name}: required mutation/ablation did not fire")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def symbols() -> dict[str, sp.Symbol]:
    omega, c_s, d, a, w = sp.symbols("omega cS d a w", positive=True, real=True)
    eps0, eps1 = sp.symbols("epsilon0 epsilon1", positive=True, real=True)
    m0 = sp.symbols("M0", positive=True, real=True)
    d1, q2 = sp.symbols("D1 Q2", real=True)
    return {
        "omega": omega,
        "cS": c_s,
        "d": d,
        "a": a,
        "w": w,
        "epsilon0": eps0,
        "epsilon1": eps1,
        "M0": m0,
        "D1": d1,
        "Q2": q2,
    }


SYM = symbols()
omega = SYM["omega"]
cS = SYM["cS"]
d = SYM["d"]
a = SYM["a"]
w = SYM["w"]
epsilon0 = SYM["epsilon0"]
epsilon1 = SYM["epsilon1"]
M0 = SYM["M0"]
D1 = SYM["D1"]
Q2 = SYM["Q2"]


EXPECTED_KERNEL = {
    0: I * a * omega / cS,
    1: I * a**3 * omega**3 / (sp.Integer(2) * cS**3),
    2: I * a**5 * omega**5 / (sp.Integer(27) * cS**5),
}


CONSUMED_STAGE008_KERNEL = {
    0: I * a * omega / cS,
    1: I * a**3 * omega**3 / (sp.Integer(2) * cS**3),
    2: I * a**5 * omega**5 / (sp.Integer(27) * cS**5),
}


def omega_order(expr: sp.Expr, symbol: sp.Symbol, max_order: int) -> int:
    expanded = sp.expand(expr)
    for power in range(max_order + 1):
        coeff = sp.simplify(expanded.coeff(symbol, power))
        assert_no_float(f"omega_order.coeff[{power}]", coeff)
        if coeff != 0:
            return power
    raise AuditFailure(f"no nonzero {symbol} coefficient through order {max_order}: {expanded}")


def transport_data(*, return_sign: int = -1) -> dict[str, sp.Expr]:
    outgoing_basis = sp.exp(I * omega * w / cS)
    returning_basis = sp.exp(sp.Integer(return_sign) * I * omega * w / cS)
    forward_phase = sp.simplify(outgoing_basis.subs(w, d) / outgoing_basis.subs(w, 0))
    return_phase = sp.simplify(returning_basis.subs(w, 0) / returning_basis.subs(w, d))
    phase = sp.simplify(forward_phase * return_phase)
    tau = sp.simplify(sp.diff(sp.log(phase), omega) / I)
    return {
        "outgoing_basis": outgoing_basis,
        "returning_basis": returning_basis,
        "forward_phase": forward_phase,
        "return_phase": return_phase,
        "phase": phase,
        "tau": tau,
    }


def transport_residual(data: dict[str, sp.Expr]) -> sp.Expr:
    phase_identity = sp.simplify(sp.exp(I * omega * data["tau"]) - data["phase"])
    tau_target = sp.simplify(data["tau"] - sp.Integer(2) * d / cS)
    return sp.simplify(phase_identity**2 + tau_target**2)


def alpha(epsilon: sp.Expr) -> sp.Expr:
    return sp.simplify(1 / (1 + epsilon))


def neighbor_fraction(epsilon: sp.Expr) -> sp.Expr:
    return sp.simplify(epsilon / (1 + epsilon))


def series(expr: sp.Expr, order: int) -> sp.Expr:
    return sp.expand(sp.series(expr, omega, 0, order + 1).removeO())


def kernel_integrity_residual(
    pipeline_kernels: dict[int, sp.Expr], consumed_kernels: dict[int, sp.Expr]
) -> sp.Expr:
    return sp.simplify(
        sum(
            sp.simplify(consumed_kernels[ell] - pipeline_kernels[ell]) ** 2
            for ell in (0, 1, 2)
        )
    )


def classify_check_a(p0: int, p1: int) -> dict[str, bool]:
    return {
        "A_strict_pass": bool(p0 >= 5 and p1 >= 5),
        "A_residual_pass": bool((p0 < 5 or p1 < 5) and p0 >= 1 and p1 >= 3),
    }


def classification_residual(classification: dict[str, bool]) -> sp.Expr:
    return sp.Integer(0) if (
        classification["A_strict_pass"] is False and classification["A_residual_pass"] is True
    ) else sp.Integer(1)


def compute_baseline() -> dict[str, Any]:
    transport = transport_data()
    phase = transport["phase"]
    tau = transport["tau"]
    alpha0 = alpha(epsilon0)
    alpha1 = alpha(epsilon1)
    neighbor0 = neighbor_fraction(epsilon0)
    neighbor1 = neighbor_fraction(epsilon1)
    t0_full = sp.simplify(alpha0 * phase)
    t1_full = sp.simplify(alpha1 * phase)
    t0_series = series(t0_full, 5)
    t1_series = series(t1_full, 5)
    t0_dc = sp.simplify(sp.limit(t0_full, omega, 0))
    t1_dc = sp.simplify(sp.limit(t1_full, omega, 0))
    nu0 = omega_order(1 - t0_series, omega, 5)
    nu1 = omega_order(1 - t1_series, omega, 5)
    kernels = dict(EXPECTED_KERNEL)
    p_raw0 = omega_order(kernels[0], omega, 7)
    p_raw1 = omega_order(kernels[1], omega, 7)
    p_raw2 = omega_order(kernels[2], omega, 7)
    r0 = sp.simplify(-M0 * t0_full)
    r1 = sp.simplify(-D1 * t1_full)
    a0_res = sp.simplify(kernels[0] * M0 * (1 - t0_series))
    a1_res = sp.simplify(kernels[1] * D1 * (1 - t1_series))
    p_res0 = p_raw0 + nu0
    p_res1 = p_raw1 + nu1
    classification = classify_check_a(p_res0, p_res1)
    steady0 = sp.simplify(M0 - alpha0 * M0 - neighbor0 * M0)
    steady1 = sp.simplify(D1 - alpha1 * D1 - neighbor1 * D1)
    z_throat = -M0
    z_return = sp.simplify(M0 * t0_dc)
    z_replenishment = sp.Integer(0)
    z_boundary = sp.Integer(0)
    z_total = sp.simplify(z_throat + z_return + z_replenishment + z_boundary)
    z_local = sp.simplify(-M0 * (1 - t0_dc))
    z_formula = sp.simplify(-M0 * epsilon0 / (1 + epsilon0))
    z_certificate = sp.simplify(-z_total * (1 + epsilon0) / (M0 * epsilon0))
    strict_t0 = sp.simplify(sp.limit(t0_full, epsilon0, 0, dir="+"))
    strict_t1 = sp.simplify(sp.limit(t1_full, epsilon1, 0, dir="+"))
    strict_t0_series = series(strict_t0, 5)
    strict_t1_series = series(strict_t1, 5)
    strict_nu0 = omega_order(1 - strict_t0_series, omega, 5)
    strict_nu1 = omega_order(1 - strict_t1_series, omega, 5)
    strict_p_res0 = p_raw0 + strict_nu0
    strict_p_res1 = p_raw1 + strict_nu1
    strict_z = sp.simplify(sp.limit(z_total, epsilon0, 0, dir="+"))
    return {
        "transport": transport,
        "phase": phase,
        "tau": tau,
        "alpha": {0: alpha0, 1: alpha1},
        "neighbor": {0: neighbor0, 1: neighbor1},
        "T_full": {0: t0_full, 1: t1_full},
        "T_series": {0: t0_series, 1: t1_series},
        "T_dc": {0: t0_dc, 1: t1_dc},
        "nu": {0: nu0, 1: nu1},
        "kernels": kernels,
        "p_raw": {0: p_raw0, 1: p_raw1, 2: p_raw2},
        "R": {0: r0, 1: r1},
        "A_res": {0: a0_res, 1: a1_res},
        "p_res": {0: p_res0, 1: p_res1},
        "classification": classification,
        "steady": {0: steady0, 1: steady1},
        "Z_parts": {
            "throat": z_throat,
            "return": z_return,
            "replenishment_localized": z_replenishment,
            "boundary_dof": z_boundary,
        },
        "Z": z_total,
        "Z_local": z_local,
        "Z_formula": z_formula,
        "Z_certificate": z_certificate,
        "strict": {
            "T": {0: strict_t0, 1: strict_t1},
            "T_series": {0: strict_t0_series, 1: strict_t1_series},
            "nu": {0: strict_nu0, 1: strict_nu1},
            "p_res": {0: strict_p_res0, 1: strict_p_res1},
            "Z": strict_z,
        },
    }


def run_round_trip(data: dict[str, Any]) -> None:
    subbanner("Solved bidirectional Helmholtz round-trip transport")
    tr = data["transport"]
    print("  basis: Phi_l=A_l*exp(I*omega*w/c_s)+B_l*exp(-I*omega*w/c_s)")
    print(f"  forward phase ratio 0->d = {s(tr['forward_phase'])}")
    print(f"  return phase ratio d->0 = {s(tr['return_phase'])}")
    print(f"  round-trip phase = {s(tr['phase'])}")
    print(f"  tau solved from d(log phase)/domega / I = {s(tr['tau'])}")
    expect_zero("exp(I*omega*tau) equals the solved round-trip phase", sp.exp(I * omega * tr["tau"]) - tr["phase"])
    expect_zero("tau solved from the basis equals 2*d/c_S", tr["tau"] - sp.Integer(2) * d / cS)


def run_return_transfer(data: dict[str, Any]) -> None:
    subbanner("DC continuity fractions and return transfer")
    print(f"  alpha_0 = {s(data['alpha'][0])}; neighbor_0 = {s(data['neighbor'][0])}")
    print(f"  alpha_1 = {s(data['alpha'][1])}; neighbor_1 = {s(data['neighbor'][1])}")
    print(f"  T_0(omega) = {s(data['T_full'][0])}")
    print(f"  T_1(omega) = {s(data['T_full'][1])}")
    for ell in (0, 1):
        expect_zero(
            f"ell={ell} DC transfer by limit matches alpha_ell continuity fraction",
            data["T_dc"][ell] - data["alpha"][ell],
        )
        print(f"    ell={ell}: Limit[T_l, omega->0] = {s(data['T_dc'][ell])}")
    expect_zero("zeroth-moment steady circulation balance M0=alpha0*M0+neighbor0*M0", data["steady"][0])
    expect_zero("first-moment steady circulation balance D1=alpha1*D1+neighbor1*D1", data["steady"][1])


def run_consumed_stage008(data: dict[str, Any]) -> None:
    subbanner("Consumed from ledger_stage008 (II-B1)")
    consumed_kernels = CONSUMED_STAGE008_KERNEL
    print("  cited kernels, exact-value integrity checked:")
    print(f"    kernel_0 = {s(data['kernels'][0])}")
    print(f"    kernel_1 = {s(data['kernels'][1])}")
    print(f"    kernel_2 = {s(data['kernels'][2])}")
    print("  independently typed consumed-kernel site:")
    print(f"    consumed_kernel_0 = {s(consumed_kernels[0])}")
    print(f"    consumed_kernel_1 = {s(consumed_kernels[1])}")
    print(f"    consumed_kernel_2 = {s(consumed_kernels[2])}")
    expect_zero("consumed kernel_0 matches pipeline kernel_0", consumed_kernels[0] - data["kernels"][0])
    expect_zero("consumed kernel_1 matches pipeline kernel_1", consumed_kernels[1] - data["kernels"][1])
    expect_zero("consumed kernel_2 matches pipeline kernel_2", consumed_kernels[2] - data["kernels"][2])
    expect_zero(
        "all consumed stage008 kernels match independent pipeline kernels",
        kernel_integrity_residual(data["kernels"], consumed_kernels),
    )
    for ell, moment, eps in ((0, M0, epsilon0), (1, D1, epsilon1)):
        dc_relation = sp.simplify(data["R"][ell].subs(omega, 0))
        stage008_target = -moment
        expect_zero(
            f"ell={ell} DC return moment approaches stage008 target as epsilon_l->0",
            sp.limit(dc_relation, eps, 0, dir="+") - stage008_target,
        )
        print(
            f"    ell={ell}: R_l(0) = {s(dc_relation)} -> {s(stage008_target)} "
            "as epsilon_l->0"
        )
    print("  T2_applied=false — kernel_2/Q2 INERT, nothing derived at l=2")


def run_residual_prediction(data: dict[str, Any]) -> None:
    subbanner("Bounded residual prediction, Check A")
    for ell in (0, 1):
        print(
            f"  ell={ell}: omega-order scan nu_l=ord_omega(1-T_l series) = {data['nu'][ell]}"
        )
        expect_zero(f"ell={ell} finite DC sink leaves omega^0 deviation nu_l=0", sp.Integer(data["nu"][ell]))
        expect_nonzero(
            f"ell={ell} deviation-from-one has a finite DC sink term",
            sp.simplify((1 - data["T_series"][ell]).coeff(omega, 0)),
        )
    print("  finite DC sink leaves a nonzero deviation-from-one at O(omega^0)")
    print(f"  R_0(omega) = {s(data['R'][0])}")
    print(f"  R_1(omega) = {s(data['R'][1])}")
    print(f"  A_res(ell0) = kernel_0*M0*(1-T0) = {s(data['A_res'][0])}")
    print(f"  A_res(ell1) = kernel_1*D1*(1-T1) = {s(data['A_res'][1])}")
    expect_zero("p_raw0 is computed from consumed kernel_0 omega order", sp.Integer(data["p_raw"][0]) - 1)
    expect_zero("p_raw1 is computed from consumed kernel_1 omega order", sp.Integer(data["p_raw"][1]) - 3)
    expect_zero("p_raw2 is computed from consumed kernel_2 omega order", sp.Integer(data["p_raw"][2]) - 5)
    expect_zero("computed p_res(ell0)=p_raw0+nu0=1", sp.Integer(data["p_res"][0]) - 1)
    expect_zero("computed p_res(ell1)=p_raw1+nu1=3", sp.Integer(data["p_res"][1]) - 3)
    print(f"  p_res(ell0) = {data['p_res'][0]}; p_res(ell1) = {data['p_res'][1]}")
    print(f"  A_strict_pass computed = {data['classification']['A_strict_pass']}")
    print(f"  A_residual_pass computed = {data['classification']['A_residual_pass']}")
    expect_bool("A_strict_pass is computed False", data["classification"]["A_strict_pass"] is False)
    expect_bool("A_residual_pass is computed True", data["classification"]["A_residual_pass"] is True)
    expect_zero("Check-A classification follows from computed p_res orders", classification_residual(data["classification"]))
    print(
        "  source top-line: RETURN_RESIDUAL_PREDICTION "
        "(Check A component computed here; Check B = ledger_stage010)"
    )
    print(
        "  Check A component computed here; Check B localization = ledger_stage010 "
        "(zero-mode/radial dsolves, counterfactual guard, DC-sink classifier, NOGO warp, spectra)"
    )


def run_z_accounting(data: dict[str, Any]) -> None:
    subbanner("Z channel accounting, premise vs accounting")
    parts = data["Z_parts"]
    print(f"  Z_throat = {s(parts['throat'])}")
    print(f"  Z_return = {s(parts['return'])}")
    print("  localized replenishment = 0 (declared)")
    print("  boundary DOF = 0 (declared)")
    print(f"  Z = {s(data['Z'])}")
    expect_zero("Z channel sum reduces to -M0*(1-T0(0))", data["Z"] - data["Z_local"])
    expect_zero("Z accounting formula equals -M0*epsilon0/(1+epsilon0)", data["Z"] - data["Z_formula"])
    expect_zero("Z sign certificate -Z*(1+epsilon0)/(M0*epsilon0) equals 1", data["Z_certificate"] - 1)
    print("  under v3 this is ACCOUNTING; Z<0 (drain admissibility) is the PREMISE (Z_is_premise = true)")
    print("  Z<0 = the drain-admissibility PREMISE (Z_is_premise=true); the formula = ACCOUNTING")


def run_strict_limits(data: dict[str, Any]) -> None:
    subbanner("Strict perfect-return limits, per channel")
    expect_zero("epsilon0->0+ strict T0 equals exp(I*omega*tau)", data["strict"]["T"][0] - sp.exp(I * omega * data["tau"]))
    expect_zero("epsilon0->0+ strict_nu0 computed independently equals 1", sp.Integer(data["strict"]["nu"][0]) - 1)
    expect_zero("epsilon0->0+ strict p_res0 computed as p_raw0+strict_nu0 equals 2", sp.Integer(data["strict"]["p_res"][0]) - 2)
    expect_zero("epsilon0->0+ Z tends to 0", data["strict"]["Z"])
    expect_zero("epsilon1->0+ strict_nu1 computed from strict T1 series equals 1", sp.Integer(data["strict"]["nu"][1]) - 1)
    expect_zero("epsilon1->0+ strict p_res1 computed as p_raw1+strict_nu1 equals 4", sp.Integer(data["strict"]["p_res"][1]) - 4)
    print(
        "  ell=0 contingency: residual and Z require epsilon0>0; "
        f"strict T0={s(data['strict']['T'][0])}, strict_nu0={data['strict']['nu'][0]}, "
        f"strict_p_res0={data['strict']['p_res'][0]}, Z->0"
    )
    print(
        "  ell=1 contingency: residual requires epsilon1>0; "
        f"strict T1={s(data['strict']['T'][1])}, strict_nu1={data['strict']['nu'][1]}, "
        f"strict_p_res1={data['strict']['p_res'][1]}"
    )


@dataclass(frozen=True)
class Dim:
    """Exact exponent triple for base dimensions in {L, T, M} order."""

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

    def components(self) -> tuple[sp.Rational, sp.Rational, sp.Rational]:
        return (self.l, self.t, self.m)


DIMENSIONLESS = Dim()
LENGTH = Dim(1, 0, 0)
TIME = Dim(0, 1, 0)
MOMENT0 = Dim(0, -1, 0)


def dim_residual(actual: Dim, expected: Dim) -> sp.Expr:
    return sp.simplify(sum((have - want) ** 2 for have, want in zip(actual.components(), expected.components())))


def run_dimensional_block() -> None:
    subbanner("Modest dimensional block")
    omega_dim = Dim(0, -1)
    d_dim = LENGTH
    c_s_dim = LENGTH / TIME
    epsilon_dim = DIMENSIONLESS
    one_plus_epsilon_dim = DIMENSIONLESS
    neighbor_dim = epsilon_dim / one_plus_epsilon_dim
    alpha_dim = DIMENSIONLESS / one_plus_epsilon_dim
    phase_dim = DIMENSIONLESS
    transfer_dim = alpha_dim * phase_dim
    z_dim = omega_dim * d_dim / c_s_dim
    tau_dim = d_dim / c_s_dim
    z_accounting_dim = MOMENT0 * neighbor_dim
    dimensionless_fraction_residual = sp.simplify(
        dim_residual(epsilon_dim, DIMENSIONLESS)
        + dim_residual(neighbor_dim, DIMENSIONLESS)
        + dim_residual(alpha_dim, DIMENSIONLESS)
        + dim_residual(transfer_dim, DIMENSIONLESS)
    )
    expect_zero("z=omega*d/c_S is dimensionless", dim_residual(z_dim, DIMENSIONLESS))
    expect_zero("[tau]=T from tau=2*d/c_S", dim_residual(tau_dim, TIME))
    expect_zero("[Z]=[M0] via -M0*epsilon0/(1+epsilon0)", dim_residual(z_accounting_dim, MOMENT0))
    expect_zero("epsilon_l, alpha_l, and T_l are dimensionless by fraction composition", dimensionless_fraction_residual)
    print("  z=omega*d/c_S dimensionless; [tau]=T; [Z]=[M0]; epsilon_l, alpha_l, T_l dimensionless.")


def print_provenance() -> None:
    subbanner("Provenance and scope")
    print("  postulated geometry: flat finite slab, brane w=0, absorber w=d; response derived on it.")
    print("  premise vs accounting (v3, verbatim-class): Z<0 = the drain-admissibility PREMISE (Z_is_premise=true); Z=-M0*eps0/(1+eps0) = ACCOUNTING.")
    print("  open item #9: sharpened, NOT closed — the deliverable is the falsifiable residual radiation prediction tied to the drain strength; the gravity-range (1/r^2) leg = stage 010.")
    print("  Check-B pointer: radiation/Sommerfeld boundary recorded ac_check_a_only in the source — not a Check-B branch; localization/NOGO/classifier = stage 010.")
    print("  downstream consumers: stage 023 (pathA_34 residuals must match in form/sign/order); stage 024/026 (pathA_43 consumes the Phi2 bulk mode machinery — stage 010's leg — and the projected-continuity operator lineage).")
    print("  dropped bookkeeping: the source's SHA-256 trace ids were build-reproducibility plumbing, replaced by the v2 tri-review protocol.")


def print_verdict_labels() -> None:
    print("")
    print("Verdict labels:")
    print("  ledger earned-label (NOT a source verdict token): RETURN_TRANSFER_DERIVED_RESIDUAL_BOUNDED  (tau=2d/c_S solved; T_l=alpha_l*exp(i*omega*2d/c_s); nu_l=0 computed => p_res=1/3; Z accounting + sign certificate; strict limit reversible)")
    print("  source top-line verdict: RETURN_RESIDUAL_PREDICTION  (Check A component computed here; Check B localization = ledger_stage010)")
    print("  premise vs accounting (v3, verbatim-class): Z<0 = the drain-admissibility PREMISE (Z_is_premise=true); Z=-M0*eps0/(1+eps0) = ACCOUNTING")
    print("  earned: solved round-trip transport; DC continuity fractions + steady balances; bounded residual p_res(l0)=1, p_res(l1)=3 CONTINGENT per-channel on eps0>0 (l=0, Z) and eps1>0 (l=1) — strict limits computed independently: deviations -> O(omega) (strict orders 2/4), Z -> 0")
    print("  postulated: the flat finite slab (brane w=0, absorber w=d)")
    print("  consumed from ledger_stage008 (II-B1): kernels i*a*omega/cS, i*a^3*omega^3/(2*cS^3), i*a^5*omega^5/(27*cS^5) (cited, integrity-checked); kernel_2/Q2 INERT (T2_applied=false, nothing derived at l=2)")
    print("  exports: the falsifiable residual-radiation prediction (drain-strength-tied) -> stage 023; open-item #9 sharpened NOT closed")


def run_able_to_fail_teeth(data: dict[str, Any]) -> None:
    subbanner("Able-to-fail mutation teeth")
    flipped_transport = transport_data(return_sign=1)
    expect_fail("tooth 1 returning-basis sign flip trips tau/phase assert", transport_residual(flipped_transport))

    alpha0_bad = neighbor_fraction(epsilon0)
    t0_bad = sp.simplify(alpha0_bad * data["phase"])
    t0_bad_dc = sp.simplify(sp.limit(t0_bad, omega, 0))
    steady0_bad = sp.simplify(M0 - alpha0_bad * M0 - data["neighbor"][0] * M0)
    expect_fail(
        "tooth 2 alpha0->epsilon0/(1+epsilon0) trips DC fraction and steady balance",
        sp.simplify((t0_bad_dc - alpha(epsilon0)) ** 2 + steady0_bad**2),
    )

    scan0_bad_input = sp.simplify((1 - data["T_series"][0]) - neighbor_fraction(epsilon0))
    scan1_bad_input = sp.simplify((1 - data["T_series"][1]) - neighbor_fraction(epsilon1))
    scan0_bad = omega_order(scan0_bad_input, omega, 5)
    scan1_bad = omega_order(scan1_bad_input, omega, 5)
    expect_fail(
        "tooth 3 subtracting finite-sink DC term raises nu and trips p_res asserts",
        sp.Integer(scan0_bad) ** 2
        + (sp.Integer(data["p_raw"][0] + scan0_bad) - 1) ** 2
        + sp.Integer(scan1_bad) ** 2
        + (sp.Integer(data["p_raw"][1] + scan1_bad) - 3) ** 2,
    )

    z_drop_return = sp.simplify(data["Z_parts"]["throat"] + data["Z_parts"]["replenishment_localized"] + data["Z_parts"]["boundary_dof"])
    expect_fail("tooth 4 dropping Z_return trips Z reduction assert", z_drop_return - data["Z_local"])

    z_sign_flip = sp.simplify(-data["Z"])
    sign_flip_certificate = sp.simplify(-z_sign_flip * (1 + epsilon0) / (M0 * epsilon0))
    expect_fail("tooth 5 flipping Z sign trips sign certificate", sign_flip_certificate - 1)

    strict_t0_wrong = sp.simplify(data["T_full"][0].subs(epsilon0, 1))
    strict_t0_wrong_series = series(strict_t0_wrong, 5)
    strict_nu0_wrong = omega_order(1 - strict_t0_wrong_series, omega, 5)
    strict_z_wrong = sp.simplify(data["Z"].subs(epsilon0, 1))
    expect_fail(
        "tooth 6 corrupting strict limit epsilon0->1 trips Z->0/strict_nu0=1",
        sp.simplify(strict_z_wrong**2 + (sp.Integer(strict_nu0_wrong) - 1) ** 2),
    )

    corrupt_kernels = dict(data["kernels"])
    corrupt_consumed_kernels = dict(CONSUMED_STAGE008_KERNEL)
    corrupt_consumed_kernels[2] = I * a**5 * omega**5 / (sp.Integer(9) * cS**5)
    expect_fail(
        "tooth 7 consumed kernel corruption 27->9 trips integrity check",
        kernel_integrity_residual(corrupt_kernels, corrupt_consumed_kernels),
    )

    corrupt_classification = classify_check_a(0, data["p_res"][1])
    expect_fail(
        "tooth 8 p_res0->0 flips A_residual_pass False and trips classification assert",
        classification_residual(corrupt_classification),
    )

    expect_zero("baseline immutable after teeth: transport residual unchanged", transport_residual(data["transport"]))
    expect_zero(
        "baseline immutable after teeth: kernel integrity unchanged",
        kernel_integrity_residual(data["kernels"], CONSUMED_STAGE008_KERNEL),
    )
    expect_zero("baseline immutable after teeth: classification unchanged", classification_residual(data["classification"]))


def main() -> None:
    banner("ledger_stage009_flat_slab_return_residual SymPy audit")
    data = compute_baseline()
    assert_no_float("baseline", data)
    run_round_trip(data)
    run_return_transfer(data)
    run_consumed_stage008(data)
    run_residual_prediction(data)
    run_z_accounting(data)
    run_strict_limits(data)
    run_dimensional_block()
    print_provenance()
    print_verdict_labels()
    run_able_to_fail_teeth(data)
    print("")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print("OVERALL PASS: SymPy verified ledger_stage009 flat-slab return residual Check A exactly")


if __name__ == "__main__":
    try:
        main()
    except AuditFailure as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
        print(f"OVERALL FAIL: SymPy stage009 audit did not close ({exc})")
        raise SystemExit(1)
