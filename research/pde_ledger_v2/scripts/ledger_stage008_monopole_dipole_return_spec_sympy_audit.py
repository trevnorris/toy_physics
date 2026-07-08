#!/usr/bin/env python3
"""Ledger stage008 SymPy audit: monopole/dipole return constraint spec.

Standalone, print-only, no arguments, no file output.  This reshapes the
pathA_28 SymPy route into a ledger citizen with an earned Hankel DtN ladder,
honest return-target bookkeeping, live-source leakage anchors, and able-to-fail
in-memory mutation teeth.
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


def compact(expr: sp.Expr | int) -> str:
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
    _record_fail(f"FAIL  {name}: residual = {compact(clean)}")
    raise AuditFailure(f"{name} residual was not zero")


def expect_bool(name: str, condition: bool) -> None:
    expect_zero(name, sp.Integer(0) if condition else sp.Integer(1))


def expect_nonzero(name: str, residual: sp.Expr | int) -> None:
    assert_no_float(name, residual)
    clean = sp.simplify(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_pass(f"PASS  {name} is nonzero as required (residual = {compact(clean)})")
        return
    _record_fail(f"FAIL  {name}: required nonzero residual vanished")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def expect_fail(name: str, residual: sp.Expr | int) -> None:
    assert_no_float(name, residual)
    clean = sp.simplify(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_pass(f"PASS  {name} produced required FAIL (residual = {compact(clean)})")
        return
    _record_fail(f"FAIL  {name}: required mutation/ablation did not fire")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def nonzero(expr: sp.Expr) -> bool:
    return not bool(sp.simplify(expr) == 0)


def bool_residual(condition: bool) -> sp.Integer:
    return sp.Integer(0) if condition else sp.Integer(1)


def nonzero_assert_residual(expr: sp.Expr) -> sp.Integer:
    return sp.Integer(0) if nonzero(expr) else sp.Integer(1)


def same_string_residual(actual: str, expected: str) -> sp.Integer:
    return sp.Integer(0) if actual == expected else sp.Integer(1)


def s(expr: sp.Expr | int) -> str:
    return sp.sstr(sp.factor(sp.cancel(sp.simplify(expr))))


def symbols() -> dict[str, sp.Symbol]:
    a, k, omega, c_s, z, eta0, eta1 = sp.symbols(
        "a k omega cS z eta0 eta1",
        positive=True,
        real=True,
    )
    m0, d1, q2, r0, r1 = sp.symbols("M0 D1 Q2 R0 R1", real=True)
    return {
        "a": a,
        "k": k,
        "omega": omega,
        "cS": c_s,
        "z": z,
        "M0": m0,
        "D1": d1,
        "Q2": q2,
        "R0": r0,
        "R1": r1,
        "eta0": eta0,
        "eta1": eta1,
    }


SYM = symbols()
a = SYM["a"]
k = SYM["k"]
omega = SYM["omega"]
cS = SYM["cS"]
z = SYM["z"]
M0 = SYM["M0"]
D1 = SYM["D1"]
Q2 = SYM["Q2"]
R0 = SYM["R0"]
R1 = SYM["R1"]
eta0 = SYM["eta0"]
eta1 = SYM["eta1"]


def hankel_components(ell: int, z_symbol: sp.Symbol, *, corrupt_j1: bool = False) -> tuple[sp.Expr, sp.Expr]:
    one = sp.Integer(1)
    three = sp.Integer(3)
    if ell == 0:
        return sp.sin(z_symbol) / z_symbol, -sp.cos(z_symbol) / z_symbol
    if ell == 1:
        j1 = sp.sin(z_symbol) / z_symbol**2
        if not corrupt_j1:
            j1 -= sp.cos(z_symbol) / z_symbol
        y1 = -sp.cos(z_symbol) / z_symbol**2 - sp.sin(z_symbol) / z_symbol
        return j1, y1
    if ell == 2:
        j2 = ((three / z_symbol**3) - one / z_symbol) * sp.sin(z_symbol) - three * sp.cos(z_symbol) / z_symbol**2
        y2 = -((three / z_symbol**3) - one / z_symbol) * sp.cos(z_symbol) - three * sp.sin(z_symbol) / z_symbol**2
        return j2, y2
    raise AuditFailure(f"unsupported ell={ell}")


def first_radiating_power(y_norm: sp.Expr, max_order: int = 8) -> tuple[int, sp.Expr]:
    expanded = sp.expand(y_norm)
    for power in range(1, max_order + 1):
        coeff = sp.simplify(expanded.coeff(k, power))
        imag_coeff = sp.simplify(coeff / I)
        if coeff != 0 and imag_coeff != 0 and not imag_coeff.has(I):
            return power, imag_coeff
    raise AuditFailure(f"no imaginary radiating power through k^{max_order}: {expanded}")


def dtn_data(*, corrupt_j1: bool = False) -> dict[int, dict[str, sp.Expr | int]]:
    out: dict[int, dict[str, sp.Expr | int]] = {}
    for ell in (0, 1, 2):
        j_ell, y_ell = hankel_components(ell, z, corrupt_j1=(corrupt_j1 and ell == 1))
        h_ell = sp.simplify(j_ell + I * y_ell)
        lam = sp.simplify((k * sp.diff(h_ell, z) / h_ell).subs(z, k * a))
        lam_series = sp.expand(sp.series(lam, k, 0, 7).removeO())
        admittance_series = sp.expand(sp.series(1 / lam_series, k, 0, 7).removeO())
        static_value = sp.simplify(admittance_series.subs(k, 0))
        normalized = sp.expand(sp.series(admittance_series / static_value, k, 0, 7).removeO())
        p_raw, imag_coeff = first_radiating_power(normalized)
        kernel = sp.simplify(I * imag_coeff * (omega / cS) ** p_raw)
        out[ell] = {
            "h": h_ell,
            "Lambda_series": lam_series,
            "admittance_series": admittance_series,
            "static_value": static_value,
            "Y_norm_series": normalized,
            "p_raw": p_raw,
            "imag_coeff_k": imag_coeff,
            "radiation_kernel": kernel,
        }
    return out


EXPECTED_P = {0: 1, 1: 3, 2: 5}
EXPECTED_COEFF = {0: a, 1: a**3 / sp.Integer(2), 2: a**5 / sp.Integer(27)}
EXPECTED_KERNEL = {
    0: I * a * omega / cS,
    1: I * a**3 * omega**3 / (sp.Integer(2) * cS**3),
    2: I * a**5 * omega**5 / (sp.Integer(27) * cS**5),
}


def dtn_ladder_residual(dtn: dict[int, dict[str, sp.Expr | int]]) -> sp.Expr:
    residual = sp.Integer(0)
    for ell in (0, 1, 2):
        residual += (sp.Integer(dtn[ell]["p_raw"]) - EXPECTED_P[ell]) ** 2
        residual += sp.simplify(dtn[ell]["imag_coeff_k"] - EXPECTED_COEFF[ell]) ** 2
    return sp.simplify(residual)


def kernel_residual(dtn: dict[int, dict[str, sp.Expr | int]]) -> sp.Expr:
    return sp.simplify(
        sum(sp.simplify(dtn[ell]["radiation_kernel"] - EXPECTED_KERNEL[ell]) ** 2 for ell in (0, 1, 2))
    )


def copied_with_corrupt_ell2_kernel(dtn: dict[int, dict[str, sp.Expr | int]]) -> dict[int, dict[str, sp.Expr | int]]:
    mutated = {ell: dict(dtn[ell]) for ell in (0, 1, 2)}
    mutated[2]["radiation_kernel"] = I * a**5 * omega**5 / (sp.Integer(9) * cS**5)
    return mutated


def source_and_residual_data(dtn: dict[int, dict[str, sp.Expr | int]]) -> dict[str, dict[int, sp.Expr]]:
    sources = {0: M0, 1: D1, 2: Q2}
    raw = {ell: sp.simplify(dtn[ell]["radiation_kernel"] * sources[ell]) for ell in (0, 1, 2)}
    residual = {
        0: sp.simplify(dtn[0]["radiation_kernel"] * (M0 + R0)),
        1: sp.simplify(dtn[1]["radiation_kernel"] * (D1 + R1)),
    }
    without = {0: sp.simplify(residual[0].subs(R0, 0)), 1: sp.simplify(residual[1].subs(R1, 0))}
    with_return = {
        0: sp.simplify(residual[0].subs(R0, -M0)),
        1: sp.simplify(residual[1].subs(R1, -D1)),
    }
    derivative_vertex = {
        0: sp.simplify(eta0**2 * omega**2 * raw[0]),
        1: sp.simplify(eta1**2 * omega**2 * raw[1]),
    }
    return {
        "sources": sources,
        "raw": raw,
        "residual": residual,
        "without": without,
        "with": with_return,
        "derivative_vertex": derivative_vertex,
    }


def raw_amplitude_residual(src: dict[str, dict[int, sp.Expr]]) -> sp.Expr:
    expected_raw = {0: EXPECTED_KERNEL[0] * M0, 1: EXPECTED_KERNEL[1] * D1, 2: EXPECTED_KERNEL[2] * Q2}
    return sp.simplify(sum(sp.simplify(src["raw"][ell] - expected_raw[ell]) ** 2 for ell in (0, 1, 2)))


def condition_works(src: dict[str, dict[int, sp.Expr]]) -> bool:
    return all(nonzero(src["without"][ell]) for ell in (0, 1)) and all(src["with"][ell] == 0 for ell in (0, 1))


def residual_with_corrupt_return(src: dict[str, dict[int, sp.Expr]]) -> sp.Expr:
    corrupt0 = sp.simplify(src["residual"][0].subs(R0, -sp.Integer(2) * M0))
    corrupt1 = sp.simplify(src["residual"][1].subs(R1, -D1))
    return sp.simplify(corrupt0**2 + corrupt1**2)


def run_dtn_ladder() -> dict[int, dict[str, sp.Expr | int]]:
    subbanner("Earned outgoing Hankel DtN ladder")
    dtn = dtn_data()
    for ell in (0, 1, 2):
        print(f"  ell={ell}: h_ell = {s(dtn[ell]['h'])}")
        print(f"    Lambda_ell series = {s(dtn[ell]['Lambda_series'])}")
        print(f"    static-normalized admittance = {s(dtn[ell]['Y_norm_series'])}")
        print(
            "    coefficient-scan first radiation phase: "
            f"p_raw={dtn[ell]['p_raw']}, coeff={s(dtn[ell]['imag_coeff_k'])}"
        )
        expect_zero(f"ell={ell} scanned p_raw and radiation-phase coefficient match ladder", (
            sp.Integer(dtn[ell]["p_raw"]) - EXPECTED_P[ell]
        ) ** 2 + sp.simplify(dtn[ell]["imag_coeff_k"] - EXPECTED_COEFF[ell]) ** 2)
        expect_zero(f"ell={ell} kernel is i*(omega/c_S)^p times scanned coefficient", sp.simplify(dtn[ell]["radiation_kernel"] - EXPECTED_KERNEL[ell]))
    expect_zero("DtN ladder residual p=1/3/5 and coefficients a,a^3/2,a^5/27", dtn_ladder_residual(dtn))
    expect_zero("radiation kernels match i*a*w/cS, i*a^3*w^3/(2*cS^3), i*a^5*w^5/(27*cS^5)", kernel_residual(dtn))
    return dtn


def run_sources_and_bookkeeping(dtn: dict[int, dict[str, sp.Expr | int]]) -> dict[str, dict[int, sp.Expr]]:
    subbanner("Raw amplitudes and return-target bookkeeping")
    src = source_and_residual_data(dtn)
    print("  Moment definitions printed as target definitions:")
    print("    M0(omega)=int_brane S_leak(omega,x) d^3x")
    print("    D1_i(omega)=int_brane x_i S_leak(omega,x) d^3x + int_brane j_i(omega,x) d^3x, including the carried odd wake")
    print("    Q2 is a FREE ANCHOR symbol for the downstream quadrupole consumer, not derived ell=2 physics here")
    for ell, label in ((0, "M0"), (1, "D1"), (2, "Q2")):
        print(f"  ell={ell}: A_raw = kernel*{label} = {s(src['raw'][ell])}")
        expect_nonzero(f"ell={ell} raw amplitude is present", src["raw"][ell])
    expect_zero("raw amplitude closed forms match report lines 17-19", raw_amplitude_residual(src))

    print("  Residual label: x−x bookkeeping identity, NOT an earned cancellation result")
    for ell in (0, 1):
        print(f"    ell={ell} without return: {s(src['without'][ell])}")
        print(f"    ell={ell} with target return: {s(src['with'][ell])}")
        expect_nonzero(f"ell={ell} residual without return condition is nonzero", src["without"][ell])
        expect_zero(f"ell={ell} residual with bookkeeping target is exactly zero", src["with"][ell])
    print("  Required targets: R0(omega)=-M0(omega); R1_i(omega)=-D1_i(omega)")

    print("  Raw-vs-vertex note: derivative outlet vertex g_W0(omega)=eta*omega is branch_assumption; two vertices add eta^2*omega^2 and are NOT verdict-bearing.")
    for ell in (0, 1):
        print(f"    derivative-vertex branch ell={ell}: {s(src['derivative_vertex'][ell])}")
        expect_nonzero(f"branch_assumption derivative-vertex ell={ell} remains nonzero", src["derivative_vertex"][ell])

    steady_limit = sp.simplify(sp.limit(src["raw"][0], omega, 0))
    expect_zero("steady control lim_{omega->0} raw monopole amplitude is zero", steady_limit)
    expect_bool("dominance p(ell0)<p(ell2) from scanned powers", int(dtn[0]["p_raw"]) < int(dtn[2]["p_raw"]))
    return src


def stage_anchors() -> dict[str, sp.Expr]:
    subbanner("Stage-243/244 S_leak anchors by direct integration")
    w = sp.symbols("w", real=True)
    ellw, j0, E0 = sp.symbols("ellw j0 E0", real=True)
    lam, muw, rho0 = sp.symbols("lambda mu_w rho0", positive=True, real=True)

    w243 = sp.exp(-w**2) / sp.sqrt(sp.pi)
    jw243 = ellw * j0 * w * sp.exp(-w**2)
    sleak243 = sp.simplify(sp.integrate(sp.diff(w243, w) * jw243, (w, -sp.oo, sp.oo)))
    target243 = -sp.sqrt(2) * ellw * j0 / sp.Integer(4)

    w244 = sp.exp(-w**2 / lam**2) / (lam * sp.sqrt(sp.pi))
    phi = sp.Integer(2) * w * sp.exp(-w**2 / lam**2) / (sp.sqrt(sp.pi) * lam**3)
    ew = -E0 * phi
    jw244 = muw * rho0 * ew
    sleak244 = sp.simplify(sp.integrate(sp.diff(w244, w) * jw244, (w, -sp.oo, sp.oo)))
    target244 = sp.sqrt(2) * muw * rho0 * E0 / (sp.Integer(2) * sp.sqrt(sp.pi) * lam**3)

    print(f"  stage-243 direct integral int W' j^w dw = {s(sleak243)}")
    print(f"  stage-244 direct integral S_leak = {s(sleak244)}")
    print("  ledger_stage006 (I-3) OWNS the recovery-reduction derivation of the projected law; here the two closed forms serve only as live-source consistency anchors.")
    expect_zero("stage-243 anchor equals -sqrt(2)*ell_w*j0/4", sleak243 - target243)
    expect_zero("stage-244 anchor equals sqrt(2)*mu_w*rho0*E0/(2*sqrt(pi)*lambda^3)", sleak244 - target244)
    expect_nonzero("stage-243 leakage lane is symbolically nonzero before recovery shortcut", sleak243)
    expect_nonzero("stage-244 leakage lane is symbolically nonzero before recovery shortcut", sleak244)
    return {
        "stage243": sleak243,
        "target243": target243,
        "stage244": sleak244,
        "target244": target244,
        "ellw": ellw,
        "E0": E0,
        "muw": muw,
        "rho0": rho0,
        "lambda": lam,
    }


def compute_verdict(raw_present: bool, condition_ok: bool, cancellation_possible: bool) -> str:
    if raw_present and not cancellation_possible:
        return "MONOPOLE_RADIATION_UNAVOIDABLE"
    if raw_present and condition_ok and cancellation_possible:
        return "MONOPOLE_DIPOLE_RETURN_CONDITIONAL"
    return "INCONCLUSIVE"


def run_verdict(src: dict[str, dict[int, sp.Expr]]) -> tuple[str, bool, bool]:
    subbanner("Computed verdict with honest literal flag")
    raw_present = all(nonzero(src["raw"][ell]) for ell in (0, 1, 2))
    condition_ok = condition_works(src)
    cancellation_possible = True
    print(f"  raw_present computed from amplitudes = {raw_present}")
    print(f"  condition_works computed from residual pair = {condition_ok}")
    print("  SCOPE: parameter, not computed — track-3 decides: cancellation_possible=True")
    verdict = compute_verdict(raw_present, condition_ok, cancellation_possible)
    synthetic = compute_verdict(True, False, False)
    inconclusive = compute_verdict(False, False, True)
    print(f"  baseline verdict = {verdict}")
    print(f"  synthetic no-return control = {synthetic}")
    expect_bool("baseline computed verdict is MONOPOLE_DIPOLE_RETURN_CONDITIONAL", verdict == "MONOPOLE_DIPOLE_RETURN_CONDITIONAL")
    expect_bool("synthetic (True,False,False) reaches MONOPOLE_RADIATION_UNAVOIDABLE", synthetic == "MONOPOLE_RADIATION_UNAVOIDABLE")
    expect_bool("INCONCLUSIVE verdict rung remains reachable", inconclusive == "INCONCLUSIVE")
    return verdict, raw_present, condition_ok


def run_guards(dtn: dict[int, dict[str, sp.Expr | int]], src: dict[str, dict[int, sp.Expr]], anchors: dict[str, sp.Expr]) -> None:
    subbanner("Nine source-live guards")
    print("  These controls do NOT test whether suppression occurs. What they confirm is narrow: the source moments M0/D1 are kept live (no S_leak=0, no strict-recovery basis, no projection-locking that would zero them out by construction). Beyond keeping the source live, the controls pass by construction — they are not able-to-fail probes of the physical question. Treat them as guards against the obvious tautologies, not as evidence of suppression-vs-unavoidable.")
    raw_present = all(nonzero(src["raw"][ell]) for ell in (0, 1, 2))
    steady_limit = sp.simplify(sp.limit(src["raw"][0], omega, 0))
    m0_kill = sp.simplify(src["raw"][0].subs(M0, 0))
    recovery_243 = sp.simplify(anchors["stage243"].subs(anchors["ellw"], 0))
    recovery_244 = sp.simplify(anchors["stage244"].subs(anchors["E0"], 0))
    derivative_not_basis = nonzero(src["derivative_vertex"][0])
    quadrupole_survives = int(dtn[2]["p_raw"]) == 5 and nonzero(src["raw"][2])
    no_track3_bulk_kill = nonzero(src["without"][0]) and nonzero(src["without"][1])

    expect_bool("guard raw_monopole_present: raw monopole is live in same pipeline", nonzero(src["raw"][0]) and raw_present)
    expect_zero("guard steady_no_radiation: lim_{omega->0} raw0=0", steady_limit)
    expect_bool("guard quadrupole_survives: scanned p(ell2)=5 and raw2 nonzero", quadrupole_survives)
    expect_bool("guard return_necessity: without nonzero and with exactly zero", condition_works(src))
    expect_zero("guard anti_tautology_no_S_leak_zero: M0->0 kills raw0, so shortcut is visible", m0_kill)
    print("    declaration anti_tautology_no_S_leak_zero: no S_leak=0 shortcut is used as a verdict basis.")
    print("    declaration anti_tautology_no_strict_recovery_basis: ell_w->0 and E0->0 recovery slices are not used as verdict bases.")
    print("    observed_but_not_used — the strict-recovery limit exists but is NOT taken as a basis: sleak243|_{ell_w->0} = 0")
    expect_zero("observed_but_not_used strict-recovery slice sleak243|_{ell_w->0}=0", recovery_243)
    print("    observed_but_not_used — the strict-recovery limit exists but is NOT taken as a basis: sleak244|_{E0->0} = 0")
    expect_zero("observed_but_not_used strict-recovery slice sleak244|_{E0->0}=0", recovery_244)
    print("    declaration anti_tautology_no_projection_locking: M0 and D1 stay unconstrained source moments.")
    expect_bool("guard anti_tautology_derivative_vertex_not_basis: branch nonzero but not verdict-bearing", derivative_not_basis)
    expect_bool("guard anti_tautology_no_track3_bulk_kill: without-return residuals stay nonzero", no_track3_bulk_kill)
    print("    declaration anti_tautology_no_track3_bulk_kill: no track-3 bulk return kill is imported here.")


def run_frozen_slice_and_consumed_stage005() -> dict[str, sp.Expr]:
    subbanner("Frozen slice and consumed ledger_stage005 law")
    g_const = sp.Integer(1)
    c_const = sp.Integer(1)
    c_s_const = sp.Integer(1)
    k_eos = sp.Rational(1, 500)
    rho = sp.sqrt(10)
    m_gnls = sp.Integer(1)
    a_star = sp.Rational(4731, 2500)
    l_star = sp.Rational(18121, 10000)
    print("  Frozen slice (CALIBRATED, exact): G=c=c_s=1; K_eos=1/500; rho=sqrt(10); (a*,L*)=(4731/2500,18121/10000)")
    print("  Consumed from ledger_stage005 (I-2)")
    print("    cited symbolic law: c_s^2 = 5*K*rho^4/m_GNLS")
    print("    no EOS re-derivation is performed in stage008; this is citation integrity only.")
    consumed_residual = sp.simplify(5 * k_eos * rho**4 / m_gnls - 1)
    print(f"    frozen slice (CALIBRATED, cited): G={s(g_const)}")
    print(f"    frozen slice (CALIBRATED, cited): c={s(c_const)}")
    print(f"    frozen slice (CALIBRATED, cited): c_s={s(c_s_const)}")
    print(f"    frozen slice (CALIBRATED, cited): K_eos={s(k_eos)}")
    expect_zero("frozen rho=sqrt(10) symbolic exact", rho**2 - 10)
    print(f"    frozen slice (CALIBRATED, cited): a*={s(a_star)}")
    print(f"    frozen slice (CALIBRATED, cited): L*={s(l_star)}")
    expect_zero("consumed stage005 law exact-value citation-integrity 5*(1/500)*(sqrt(10))^4/1 - 1", consumed_residual)
    return {"K_eos": k_eos, "rho": rho, "m_GNLS": m_gnls}


@dataclass(frozen=True)
class Dim:
    l: sp.Rational | int = 0
    t: sp.Rational | int = 0

    def __post_init__(self) -> None:
        object.__setattr__(self, "l", sp.Rational(self.l))
        object.__setattr__(self, "t", sp.Rational(self.t))

    def __mul__(self, other: "Dim") -> "Dim":
        return Dim(self.l + other.l, self.t + other.t)

    def __truediv__(self, other: "Dim") -> "Dim":
        return Dim(self.l - other.l, self.t - other.t)

    def __pow__(self, power: int | sp.Rational) -> "Dim":
        p = sp.Rational(power)
        return Dim(self.l * p, self.t * p)

    def components(self) -> tuple[sp.Rational, sp.Rational]:
        return (self.l, self.t)


DIMENSIONLESS = Dim()
LENGTH = Dim(1, 0)
TIME = Dim(0, 1)


def dim_residual(actual: Dim, expected: Dim) -> sp.Expr:
    return sp.simplify(sum((have - want) ** 2 for have, want in zip(actual.components(), expected.components())))


def run_dimensional_block(dtn: dict[int, dict[str, sp.Expr | int]]) -> None:
    subbanner("Modest dimensional block")
    wave_number = LENGTH**-1
    sound_speed = LENGTH / TIME
    angular_frequency = TIME**-1
    z_dim = wave_number * LENGTH
    expect_zero("z=k*a is dimensionless", dim_residual(z_dim, DIMENSIONLESS))
    for ell in (0, 1, 2):
        p = sp.Integer(dtn[ell]["p_raw"])
        kernel_dim = ((LENGTH / sound_speed) ** p) * (angular_frequency**p)
        expect_zero(f"ell={ell} kernel structure (a/c_S)^p*omega^p is dimensionless", dim_residual(kernel_dim, DIMENSIONLESS))
    expect_bool(
        "scanned ladder powers are strictly increasing",
        int(dtn[0]["p_raw"]) < int(dtn[1]["p_raw"]) < int(dtn[2]["p_raw"]),
    )
    print("  No dimensions are invented here for M0/D1/Q2; they remain frozen-slice normalized source moments/anchor.")


def print_scope_and_provenance() -> None:
    subbanner("Scope caveat and provenance")
    print("  This is a VERIFIED CONSTRAINT-SPECIFICATION, not a falsifiable test of no-monopole-radiation. What is verified-solid: the raw ell=0/1/2 outgoing amplitudes and their radiation orders, and the exact moments/orders at which any brane<->bulk return must cancel. What this report does NOT establish: that monopole/dipole radiation is suppressed or unavoidable.")
    print("  The cancellation condition (R0=-M0, R1=-D1) is the algebraic bookkeeping of what the return must cancel. Its symbolic derivation is an identity (x - x), not a deep load-bearing result.")
    print("  The UNAVOIDABLE rung is NOT decidable in-scope: cancellation_possible is a parameter (a literal flag), not a computed quantity. Deciding suppression-vs-unavoidable requires the track-3 brane<->bulk return, which is out of scope here.")
    print("  The real falsification lives in track 3: whether an admissible return can actually deliver R0=-M0, R1=-D1. This audit only specifies the target it must hit.")
    print("")
    print("  Exports:")
    print("    R0=-M0 and R1=-D1 targets + raw amplitudes/kernels -> ledger stages 009/010 (pathA_29) and 022/023 (pathA_34).")
    print("    Q2 = FREE ANCHOR (not derived ell=2 physics) -> ledger stage 026.")
    print("  Provenance:")
    print("    DtN/outgoing expansions cite reuse of research/4d_2_5pn spherical-Hankel machinery; stage029 is the later formal DOI-cite home.")
    print("    M0/D1 are this gate's own constructions consistent with old-ledger Part-VIII projected continuity, NOT verbatim Part-VIII objects.")
    print("    ledger_stage006 (I-3) OWNS the recovery-reduction derivation of the projected law; stage008 uses stage-243/244 only as live-source consistency anchors.")
    print("    Raw-vs-vertex: derivative outlet vertex is branch_assumption, recorded but not used for the verdict.")
    print("    consumed from ledger_stage005 (I-2): c_s^2 = 5*K*rho^4/m_GNLS (cited, integrity-checked)")


def print_verdict_labels(verdict: str) -> None:
    print("")
    print("Verdict labels:")
    print("  ledger earned-label (NOT a source verdict token): DTN_LADDER_DERIVED_RETURN_TARGETS_SPECIFIED  (p_raw=1/3/5 scanned from the Hankel DtN series; raw amplitudes nonzero; cancellation targets bookkept)")
    print(f"  source top-line verdict: {verdict}")
    print("  honest scope (verbatim-class): VERIFIED CONSTRAINT-SPECIFICATION, not a falsifiable suppression test; R0=-M0 / R1=-D1 is x-x bookkeeping; cancellation_possible is a literal parameter flag (track-3 decides)")
    print("  earned: DtN ladder p=1/3/5 + kernels i*ak, i*a^3k^3/2, i*a^5k^5/27; steady limit; dominance ordering; stage-243/244 S_leak anchors (live-source consistency)")
    print("  exports: R0=-M0, R1=-D1 targets -> stages 009/010/022/023; Q2 = FREE ANCHOR (not derived l=2) -> stage 026")
    print("  guards: 9 carried (computed where they compute; declarations printed as declarations; pass-by-construction framing kept)")
    print("  frozen slice (CALIBRATED, cited): G=c=c_s=1; K_eos=1/500; (a*,L*)=(4731/2500,18121/10000)")
    print("  consumed from ledger_stage005 (I-2): c_s^2 = 5*K*rho^4/m_GNLS (cited, integrity-checked)")


def run_able_to_fail_teeth(
    dtn: dict[int, dict[str, sp.Expr | int]],
    src: dict[str, dict[int, sp.Expr]],
    anchors: dict[str, sp.Expr],
    consumed: dict[str, sp.Expr],
    raw_present: bool,
) -> None:
    subbanner("Able-to-fail mutation teeth")
    corrupt_dtn = dtn_data(corrupt_j1=True)
    expect_fail("tooth 1 hankel-form corruption changes scanned ladder", dtn_ladder_residual(corrupt_dtn))

    corrupt_kernel = copied_with_corrupt_ell2_kernel(dtn)
    expect_fail("tooth 2 kernel coefficient corruption a^5/27 -> a^5/9 trips kernel assert", kernel_residual(corrupt_kernel))

    expect_fail("tooth 3 M0->0 kills raw monopole and trips nonzero assert", nonzero_assert_residual(src["raw"][0].subs(M0, 0)))

    expect_fail("tooth 4 R0->-2*M0 breaks x-x bookkeeping identity", residual_with_corrupt_return(src))

    corrupt_verdict = compute_verdict(raw_present, False, True)
    expect_fail(
        "tooth 5 verdict machinery corruption falls to INCONCLUSIVE and trips baseline equality",
        same_string_residual(corrupt_verdict, "MONOPOLE_DIPOLE_RETURN_CONDITIONAL"),
    )

    wrong_244 = sp.sqrt(2) * anchors["muw"] * anchors["rho0"] * anchors["E0"] / (sp.sqrt(sp.pi) * anchors["lambda"] ** 3)
    expect_fail("tooth 6 stage-244 prefactor corruption trips anchor assert", anchors["stage244"] - wrong_244)

    steady_corrupt_dtn = {ell: dict(dtn[ell]) for ell in (0, 1, 2)}
    omega_power = omega ** sp.Integer(steady_corrupt_dtn[0]["p_raw"])
    kernel0 = sp.sympify(steady_corrupt_dtn[0]["radiation_kernel"])
    steady_corrupt_dtn[0]["radiation_kernel"] = sp.simplify(kernel0.subs(omega_power, 1))
    steady_corrupt_src = source_and_residual_data(steady_corrupt_dtn)
    expect_fail(
        "tooth 7 steady-limit corruption strips omega from scanned kernel0 and trips limit assert",
        sp.limit(steady_corrupt_src["raw"][0], omega, 0),
    )

    corrupt_consumed = sp.simplify(4 * consumed["K_eos"] * consumed["rho"] ** 3 / consumed["m_GNLS"] - 1)
    expect_fail("tooth 8 consumed-law corruption 5*K*rho^4 -> 4*K*rho^3 trips citation integrity", corrupt_consumed)

    expect_zero("baseline immutable after copy-mutation teeth: DtN ladder unchanged", dtn_ladder_residual(dtn))
    expect_zero("baseline immutable after copy-mutation teeth: kernels unchanged", kernel_residual(dtn))
    expect_zero("baseline immutable after copy-mutation teeth: bookkeeping with-target residuals unchanged", src["with"][0] ** 2 + src["with"][1] ** 2)


def main() -> None:
    banner("ledger_stage008_monopole_dipole_return_spec SymPy audit")
    dtn = run_dtn_ladder()
    src = run_sources_and_bookkeeping(dtn)
    anchors = stage_anchors()
    verdict, raw_present, _condition_ok = run_verdict(src)
    run_guards(dtn, src, anchors)
    consumed = run_frozen_slice_and_consumed_stage005()
    run_dimensional_block(dtn)
    print_scope_and_provenance()
    print_verdict_labels(verdict)
    run_able_to_fail_teeth(dtn, src, anchors, consumed, raw_present)
    print("")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print("OVERALL PASS: SymPy verified ledger_stage008 monopole/dipole return constraint spec exactly")


if __name__ == "__main__":
    try:
        main()
    except AuditFailure as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
        print(f"OVERALL FAIL: SymPy stage008 audit did not close ({exc})")
        raise SystemExit(1)
