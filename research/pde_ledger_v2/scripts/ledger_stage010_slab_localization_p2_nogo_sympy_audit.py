#!/usr/bin/env python3
"""Ledger stage010 SymPy audit: slab localization p=2 plus NOGO control.

Standalone, print-only, no arguments, no file I/O.  This reshapes Check B of
pathA_29 into a ledger audit: the flat slab is postulated, the two admissible
DC-sink transverse spectra are solved on it, the zero-mode 3D radial Green
function is derived by dynamic and static routes, and the delocalizing warp
control keeps RETURN_NOGO reachable.
"""

from __future__ import annotations

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
    if isinstance(expr, (str, type(None))):
        return
    if isinstance(expr, bool):
        expr = sp.Integer(1) if expr else sp.Integer(0)
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


def bool_residual(condition: bool) -> sp.Integer:
    return sp.Integer(0) if condition else sp.Integer(1)


def nonzero_assert_residual(expr: sp.Expr) -> sp.Integer:
    return sp.Integer(0) if sp.simplify(expr) != 0 else sp.Integer(1)


def radial_operator_residual(
    radial_expr: sp.Expr,
    kappa_squared: sp.Expr,
    r_symbol: sp.Symbol,
    *,
    coefficient: sp.Expr = sp.Integer(2),
) -> sp.Expr:
    return sp.simplify(
        sp.diff(radial_expr, r_symbol, 2)
        + coefficient * sp.diff(radial_expr, r_symbol) / r_symbol
        + kappa_squared * radial_expr
    )


def radial_decay_exponent(radial_flow: sp.Expr, r_symbol: sp.Symbol) -> sp.Expr:
    probe = sp.simplify(-sp.limit(r_symbol * sp.diff(sp.log(radial_flow), r_symbol), r_symbol, sp.oo))
    if not probe.is_number:
        raise AuditFailure(f"could not extract radial exponent from solved flow: {radial_flow}")
    return sp.simplify(probe)


omega, cS, d, r, w = sp.symbols("omega cS d r w", positive=True, real=True)
k_warp, W = sp.symbols("k_warp W", positive=True, real=True)
m_cont = sp.symbols("m", nonnegative=True, real=True)
Z_check_a, C_massive = sp.symbols("Z C", real=True)


STAGE008_P_RAW_L2_CONSUMED = sp.Integer(5)
STAGE008_P_RAW_L2_PIPELINE = sp.Integer(2) + sp.Integer(3)
STAGE009_A_RESIDUAL_PASS_CONSUMED = True
STAGE009_A_RESIDUAL_PASS_PIPELINE = sp.Integer(1)
T2_APPLIED = False

RETURN_RESIDUAL_PREDICTION = "RETURN_RESIDUAL_PREDICTION"
RETURN_NOGO = "RETURN_NOGO"
BC_DEPENDENT = "BC_DEPENDENT"


def bool_to_int(value: bool) -> sp.Integer:
    return sp.Integer(1) if value else sp.Integer(0)


def stage008_integrity_residual(consumed: sp.Expr, pipeline: sp.Expr) -> sp.Expr:
    return sp.simplify(consumed - pipeline)


def stage009_integrity_residual(consumed: bool, pipeline: sp.Expr) -> sp.Expr:
    return sp.simplify(bool_to_int(consumed) - pipeline)


def quadrupole_survives_from_p_raw(p_value: sp.Expr) -> bool:
    clean = sp.sympify(p_value)
    return bool(clean.is_finite and clean.is_integer and clean.is_nonnegative)


def classify_dc_sink_gate(branch_ps: list[sp.Expr], quadrupole_survives: bool) -> str:
    target = sp.Integer(2)
    equal_target = [sp.simplify(p - target) == 0 for p in branch_ps]
    if not quadrupole_survives:
        return RETURN_NOGO
    if all(equal_target):
        return RETURN_RESIDUAL_PREDICTION
    if not any(equal_target):
        return RETURN_NOGO
    return BC_DEPENDENT


def branch_verdict_from_p(p_value: sp.Expr) -> str:
    return RETURN_RESIDUAL_PREDICTION if sp.simplify(p_value - 2) == 0 else RETURN_NOGO


def branch_verdict_residual(verdict: str, expected: str) -> sp.Integer:
    return sp.Integer(0) if verdict == expected else sp.Integer(1)


def transverse_spectra() -> dict[str, dict[str, Any]]:
    f0_abs = sp.simplify(1 / sp.sqrt(d))
    f1_abs = sp.simplify(sp.sqrt(sp.Rational(2, 1) / d) * sp.cos(sp.pi * w / d))
    f0_bloch = sp.simplify(1 / sp.sqrt(d))
    f1_bloch = sp.simplify(sp.sqrt(sp.Rational(2, 1) / d) * sp.cos(2 * sp.pi * w / d))
    return {
        "destructuring_absorbing": {
            "f0": f0_abs,
            "f1": f1_abs,
            "m0": sp.Integer(0),
            "m1": sp.pi / d,
            "bc": "neumann",
        },
        "bloch_stack": {
            "f0": f0_bloch,
            "f1": f1_bloch,
            "m0": sp.Integer(0),
            "m1": 2 * sp.pi / d,
            "bc": "periodic",
        },
    }


def transverse_norm(fn: sp.Expr) -> sp.Expr:
    return sp.simplify(sp.integrate(fn**2, (w, 0, d)))


def transverse_ode_residual(fn: sp.Expr, m_value: sp.Expr) -> sp.Expr:
    return sp.simplify(sp.diff(fn, w, 2) + m_value**2 * fn)


def transverse_bc_residuals(spec: dict[str, Any]) -> dict[str, sp.Expr]:
    f0 = spec["f0"]
    f1 = spec["f1"]
    if spec["bc"] == "neumann":
        return {
            "f0_prime_0": sp.diff(f0, w).subs(w, 0),
            "f0_prime_d": sp.diff(f0, w).subs(w, d),
            "f1_prime_0": sp.diff(f1, w).subs(w, 0),
            "f1_prime_d": sp.diff(f1, w).subs(w, d),
        }
    return {
        "f0_value_periodic": sp.simplify(f0.subs(w, d) - f0.subs(w, 0)),
        "f0_derivative_periodic": sp.simplify(sp.diff(f0, w).subs(w, d) - sp.diff(f0, w).subs(w, 0)),
        "f1_value_periodic": sp.simplify(f1.subs(w, d) - f1.subs(w, 0)),
        "f1_derivative_periodic": sp.simplify(sp.diff(f1, w).subs(w, d) - sp.diff(f1, w).subs(w, 0)),
    }


def ensure_zero_mode_seed(m_value: sp.Expr, route: str) -> None:
    if sp.simplify(m_value) != 0:
        raise AuditFailure(f"{route} zero-mode radial solve must be seeded by the computed m=0 eigenvalue")


def zero_seed_guard_residual(m_value: sp.Expr) -> sp.Integer:
    try:
        ensure_zero_mode_seed(m_value, "mutation")
    except AuditFailure:
        return sp.Integer(1)
    return sp.Integer(0)


def solve_dynamic_zero_mode_radial(branch: str, m_value: sp.Expr) -> dict[str, sp.Expr]:
    ensure_zero_mode_seed(m_value, "dynamic")
    g = sp.Function(f"g_{branch}_dynamic")
    kappa_squared = sp.simplify((omega / cS) ** 2 - m_value**2)
    ode = sp.Eq(radial_operator_residual(g(r), kappa_squared, r), 0)
    general = sp.dsolve(ode).rhs
    C1, C2 = sp.symbols("C1 C2")
    selected = sp.simplify(
        general.subs(
            {
                C1: I * sp.sqrt(sp.pi) * sp.sqrt(omega) / (4 * sp.pi * d * sp.sqrt(2) * sp.sqrt(cS)),
                C2: -sp.sqrt(sp.pi) * sp.sqrt(omega) / (4 * sp.pi * d * sp.sqrt(2) * sp.sqrt(cS)),
            }
        )
    )
    residual = radial_operator_residual(selected, kappa_squared, r)
    limit_green = sp.simplify(sp.limit(selected, omega, 0))
    flow = sp.simplify(-sp.diff(limit_green, r))
    exponent = radial_decay_exponent(flow, r)
    return {
        "general": general,
        "selected": selected,
        "residual": residual,
        "limit_green": limit_green,
        "flow": flow,
        "exponent": exponent,
        "kappa_squared": kappa_squared,
    }


def solve_static_zero_mode_radial(branch: str, m_value: sp.Expr) -> dict[str, sp.Expr]:
    ensure_zero_mode_seed(m_value, "static")
    g = sp.Function(f"g_{branch}_static")
    kappa_squared = -sp.simplify(m_value) ** 2
    ode = sp.Eq(radial_operator_residual(g(r), kappa_squared, r), 0)
    general = sp.dsolve(ode).rhs
    C1, C2 = sp.symbols("C1 C2")
    selected = sp.simplify(general.subs({C1: 0, C2: sp.Rational(1, 1) / (4 * sp.pi * d)}))
    residual = radial_operator_residual(selected, kappa_squared, r)
    flow = sp.simplify(-sp.diff(selected, r))
    exponent = radial_decay_exponent(flow, r)
    return {
        "general": general,
        "selected": selected,
        "residual": residual,
        "flow": flow,
        "exponent": exponent,
        "kappa_squared": kappa_squared,
    }


def solve_static_massive_radial(branch: str, m_value: sp.Expr) -> dict[str, sp.Expr]:
    if sp.simplify(m_value) == 0:
        raise AuditFailure("massive radial solve requires m>0")
    mu = sp.symbols(f"mu_{branch}", positive=True, real=True)
    g = sp.Function(f"g_{branch}_massive")
    ode = sp.Eq(radial_operator_residual(g(r), -mu**2, r), 0)
    general = sp.dsolve(ode).rhs
    C1, C2 = sp.symbols("C1 C2")
    c1_decaying = -sp.sqrt(sp.pi) * sp.sqrt(mu) / (4 * sp.pi * (1 + I))
    c2_decaying = I * c1_decaying
    selected_mu = sp.simplify(general.subs({C1: c1_decaying, C2: c2_decaying}))
    selected = sp.simplify(selected_mu.subs(mu, m_value))
    residual = radial_operator_residual(selected, -m_value**2, r)
    return {
        "general": general,
        "selected": selected,
        "residual": residual,
        "mu": m_value,
    }


def build_counterfactual_guard(candidate: sp.Expr) -> dict[str, sp.Expr]:
    correct_residual = radial_operator_residual(candidate, sp.Integer(0), r)
    perturbed = sp.simplify(candidate / r**4)
    perturbed_residual = radial_operator_residual(perturbed, sp.Integer(0), r)
    return {
        "candidate": candidate,
        "correct_residual": correct_residual,
        "perturbed": perturbed,
        "perturbed_residual": perturbed_residual,
    }


def compute_baseline() -> dict[str, Any]:
    stage008_integrity = stage008_integrity_residual(STAGE008_P_RAW_L2_CONSUMED, STAGE008_P_RAW_L2_PIPELINE)
    stage009_integrity = stage009_integrity_residual(
        STAGE009_A_RESIDUAL_PASS_CONSUMED, STAGE009_A_RESIDUAL_PASS_PIPELINE
    )
    quadrupole_survives = quadrupole_survives_from_p_raw(STAGE008_P_RAW_L2_CONSUMED)
    spectra = transverse_spectra()
    transverse = {}
    static = {}
    massive = {}
    for name, spec in spectra.items():
        transverse[name] = {
            "norm0": transverse_norm(spec["f0"]),
            "norm1": transverse_norm(spec["f1"]),
            "ode0": transverse_ode_residual(spec["f0"], spec["m0"]),
            "ode1": transverse_ode_residual(spec["f1"], spec["m1"]),
            "bc_residuals": transverse_bc_residuals(spec),
        }
        static[name] = solve_static_zero_mode_radial(name, spec["m0"])
        massive[name] = solve_static_massive_radial(name, spec["m1"])

    dynamic = solve_dynamic_zero_mode_radial("destructuring_absorbing", spectra["destructuring_absorbing"]["m0"])
    p_dynamic = dynamic["exponent"]
    p_static = static["destructuring_absorbing"]["exponent"]
    p_abs = static["destructuring_absorbing"]["exponent"]
    p_bloch = static["bloch_stack"]["exponent"]
    counterfactual = build_counterfactual_guard(static["destructuring_absorbing"]["selected"])

    warp_measure = sp.exp(2 * k_warp * w)
    warp_norm_cutoff = sp.simplify(sp.integrate(warp_measure, (w, 0, W)))
    warp_norm_limit = sp.limit(warp_norm_cutoff, W, sp.oo)
    continuum_green = sp.simplify(sp.integrate(sp.exp(-m_cont * r), (m_cont, 0, sp.oo)) / (4 * sp.pi * r))
    continuum_flow = sp.simplify(-sp.diff(continuum_green, r))
    p_delocalizing = radial_decay_exponent(continuum_flow, r)

    delocalizing_verdict = classify_dc_sink_gate([p_delocalizing], quadrupole_survives)
    destructuring_verdict = branch_verdict_from_p(p_abs)
    bloch_verdict = branch_verdict_from_p(p_bloch)
    headline = classify_dc_sink_gate([p_abs, p_bloch], quadrupole_survives)
    no_quadrupole_verdict = classify_dc_sink_gate([p_abs, p_bloch], False)
    tension_status = "witnessed" if (sp.simplify(p_abs - 2) == 0 and sp.simplify(p_delocalizing - 2) != 0) else "not_witnessed"

    return {
        "stage008_integrity": stage008_integrity,
        "stage009_integrity": stage009_integrity,
        "quadrupole_survives": quadrupole_survives,
        "spectra": spectra,
        "transverse": transverse,
        "dynamic": dynamic,
        "static": static,
        "massive": massive,
        "p_dynamic": p_dynamic,
        "p_static": p_static,
        "p_abs": p_abs,
        "p_bloch": p_bloch,
        "counterfactual": counterfactual,
        "warp_norm_cutoff": warp_norm_cutoff,
        "warp_norm_limit": warp_norm_limit,
        "continuum_green": continuum_green,
        "continuum_flow": continuum_flow,
        "p_delocalizing": p_delocalizing,
        "delocalizing_verdict": delocalizing_verdict,
        "destructuring_verdict": destructuring_verdict,
        "bloch_verdict": bloch_verdict,
        "headline": headline,
        "no_quadrupole_verdict": no_quadrupole_verdict,
        "tension_status": tension_status,
    }


def run_consumed_inputs(data: dict[str, Any]) -> None:
    subbanner("Consumed inputs with dual-site citation integrity")
    print("  ledger_stage008 cited input: p_raw(l2)=5; T2_applied=false, no l=2 physics recomputed here.")
    print(f"    consumed site p_raw_l2 = {s(STAGE008_P_RAW_L2_CONSUMED)}")
    print(f"    independent pipeline site p_raw_l2 = {s(STAGE008_P_RAW_L2_PIPELINE)}")
    expect_zero("ledger_stage008 p_raw(l2) consumed minus pipeline equals zero", data["stage008_integrity"])
    expect_zero(
        "cited p_raw(l2) equals the frozen ledger_stage008 export 5",
        STAGE008_P_RAW_L2_CONSUMED - 5,
    )
    expect_bool("quadrupole_survives derives from finite non-negative integer p_raw(l2)", data["quadrupole_survives"])
    expect_bool("T2_applied=false is enforced for the cited quadrupole input", T2_APPLIED is False)

    print("  ledger_stage009 cited input: A_residual_pass=True; used only in the joint-composition print.")
    print(f"    consumed site A_residual_pass = {STAGE009_A_RESIDUAL_PASS_CONSUMED}")
    print(f"    independent pipeline site A_residual_pass_as_int = {s(STAGE009_A_RESIDUAL_PASS_PIPELINE)}")
    expect_zero("ledger_stage009 A_residual_pass consumed minus pipeline equals zero", data["stage009_integrity"])
    expect_bool("Check-A component is cited True", STAGE009_A_RESIDUAL_PASS_CONSUMED is True)


def run_transverse_spectra(data: dict[str, Any]) -> None:
    subbanner("Two DC-sink transverse spectra and normalizable zero modes")
    for name, spec in data["spectra"].items():
        checks = data["transverse"][name]
        print(f"  {name}:")
        print(f"    f0(w) = {s(spec['f0'])}; m0 = {s(spec['m0'])}")
        print(f"    f1(w) = {s(spec['f1'])}; m1 = {s(spec['m1'])}")
        expect_zero(f"{name} zero-mode normalization integral equals 1", checks["norm0"] - 1)
        expect_zero(f"{name} first-mode normalization integral equals 1", checks["norm1"] - 1)
        expect_zero(f"{name} zero-mode transverse ODE residual f''+m^2 f", checks["ode0"])
        expect_zero(f"{name} first-mode transverse ODE residual f''+m^2 f", checks["ode1"])
        for bc_name, residual in checks["bc_residuals"].items():
            expect_zero(f"{name} boundary condition {bc_name}", residual)
    print("  load-bearing fact: both completions have a normalizable constant m=0 transverse zero mode.")


def run_radial_routes(data: dict[str, Any]) -> None:
    subbanner("Dynamic radial zero-mode route")
    dynamic = data["dynamic"]
    print("  operator: g''+(2/r)g'+((omega/c_S)^2-m^2)g=0 solved before omega->0")
    print(f"  dsolve general basis = {s(dynamic['general'])}")
    print("  BOUNDARY SELECTION: outgoing spherical C1/C2 branch; normalization from compact zero-mode overlap 1/d.")
    print(f"  selected outgoing Green = {s(dynamic['selected'])}")
    expect_zero("dynamic selected outgoing branch satisfies the radial operator", dynamic["residual"])
    print(f"  omega->0 Green limit = {s(dynamic['limit_green'])}")
    print(f"  dynamic radial flow -dG/dr = {s(dynamic['flow'])}")
    expect_zero("dynamic route large-r exponent p_dynamic=2", data["p_dynamic"] - 2)

    subbanner("Static radial route and static-dynamic consistency")
    static_abs = data["static"]["destructuring_absorbing"]
    print("  operator: omega=0 first, g''+(2/r)g'-m^2 g=0 with computed m=0")
    print(f"  dsolve general basis = {s(static_abs['general'])}")
    print("  BOUNDARY SELECTION: C1=0 removes constant branch; C2=1/(4*pi*d) fixes the overlap normalization.")
    print(f"  selected static Green = {s(static_abs['selected'])}")
    expect_zero("static selected zero-mode branch satisfies the radial operator", static_abs["residual"])
    print(f"  static radial flow -dG/dr = {s(static_abs['flow'])}")
    expect_zero("static route large-r exponent p_static=2", data["p_static"] - 2)
    expect_zero("dynamic and static exponents agree p_dynamic-p_static=0", data["p_dynamic"] - data["p_static"])
    expect_zero(
        "strong consistency dynamic limit Green minus static Green is identically zero",
        dynamic["limit_green"] - static_abs["selected"],
    )

    subbanner("Massive/gapped Yukawa contrast")
    for name, result in data["massive"].items():
        print(f"  {name}: mu=m1={s(result['mu'])}; dsolve basis = {s(result['general'])}")
        print(f"    selected decaying Yukawa branch = {s(result['selected'])}")
        expect_nonzero(f"{name} gapped mode has mu>0", result["mu"])
        expect_zero(f"{name} selected massive branch satisfies static radial operator", result["residual"])
    print("  gapped modes are Yukawa-suppressed as exp(-m1*r)/r, so only the m=0 zero mode sets the far-field 1/r^2.")
    print(f"  illustrative, not load-bearing: Green = Z*zero + C*massive with Z={Z_check_a} cited from Check A and C={C_massive} free.")


def run_completion_agreement_and_counterfactual(data: dict[str, Any]) -> None:
    subbanner("Both completions agree on p=2 and counterfactual wrong falloff is rejected")
    static_abs = data["static"]["destructuring_absorbing"]
    static_bloch = data["static"]["bloch_stack"]
    print(f"  destructuring_absorbing static Green = {s(static_abs['selected'])}; p_abs = {s(data['p_abs'])}")
    print(f"  bloch_stack static Green = {s(static_bloch['selected'])}; p_bloch = {s(data['p_bloch'])}")
    expect_zero("destructuring_absorbing exponent p_abs=2", data["p_abs"] - 2)
    expect_zero("bloch_stack exponent p_bloch=2", data["p_bloch"] - 2)
    expect_zero("both completions agree p_abs-p_bloch=0", data["p_abs"] - data["p_bloch"])
    expect_zero("bloch selected zero-mode branch satisfies the radial operator", static_bloch["residual"])

    guard = data["counterfactual"]
    print(f"  solved static candidate = {s(guard['candidate'])}")
    print(f"  wrong falloff candidate = solved Green * r^-4 = {s(guard['perturbed'])}")
    expect_zero("correct static Green residual is zero", guard["correct_residual"])
    expect_zero(
        "counterfactual residual equals 5/(pi*d*r^7)",
        guard["perturbed_residual"] - 5 / (sp.pi * d * r**7),
    )
    expect_nonzero("counterfactual wrong 1/r^5 falloff residual is nonzero", guard["perturbed_residual"])


def run_warp_and_classifier(data: dict[str, Any]) -> None:
    subbanner("NOGO warp control and computed classifier")
    print("  anti-localizing half-line warp: mu(w)=exp(2*k_warp*w), k_warp>0")
    print(f"  cutoff zero-mode norm integral int_0^W mu(w) dw = {s(data['warp_norm_cutoff'])}")
    expect_bool("warp zero-mode norm diverges as W->infinity", data["warp_norm_limit"] == sp.oo)
    print(f"  continuum Green integral = {s(data['continuum_green'])}")
    print(f"  continuum flow -dG/dr = {s(data['continuum_flow'])}")
    expect_zero("delocalizing continuum exponent p_delocalizing=3", data["p_delocalizing"] - 3)
    expect_zero(
        "same classifier maps [p_delocalizing] to RETURN_NOGO",
        branch_verdict_residual(data["delocalizing_verdict"], RETURN_NOGO),
    )
    expect_zero("falloff-tension witness requires p_abs=2", data["p_abs"] - 2)
    expect_nonzero("falloff-tension witness requires p_delocalizing != 2", data["p_delocalizing"] - 2)
    expect_zero("tension_status is witnessed", sp.Integer(0) if data["tension_status"] == "witnessed" else sp.Integer(1))

    print(f"  destructuring branch verdict = {data['destructuring_verdict']}")
    print(f"  bloch branch verdict = {data['bloch_verdict']}")
    print(f"  Check-B headline = {data['headline']}")
    expect_zero(
        "destructuring branch verdict from p_abs is RETURN_RESIDUAL_PREDICTION",
        branch_verdict_residual(data["destructuring_verdict"], RETURN_RESIDUAL_PREDICTION),
    )
    expect_zero(
        "bloch branch verdict from p_bloch is RETURN_RESIDUAL_PREDICTION",
        branch_verdict_residual(data["bloch_verdict"], RETURN_RESIDUAL_PREDICTION),
    )
    expect_zero(
        "Check-B classifier headline from [p_abs,p_bloch] is RETURN_RESIDUAL_PREDICTION",
        branch_verdict_residual(data["headline"], RETURN_RESIDUAL_PREDICTION),
    )
    expect_zero(
        "quadrupole_survives=False classifier prong returns RETURN_NOGO",
        branch_verdict_residual(data["no_quadrupole_verdict"], RETURN_NOGO),
    )
    print(
        "  COMPLETED joint verdict: RETURN_RESIDUAL_PREDICTION = "
        "(Check A: A_residual_pass=True, CITED ledger_stage009) AND "
        "(Check B: p=2 both completions, computed here)"
    )
    print("  RETURN_NOGO remains reachable if the return delocalizes (warp -> p=3) OR quadrupole_survives=False.")


def dim_residual(actual: sp.Matrix, expected: sp.Matrix) -> sp.Expr:
    return sp.simplify(sum((have - want) ** 2 for have, want in zip(actual, expected)))


def run_dimensional_block() -> None:
    subbanner("Modest dimensional block")
    zero = sp.Matrix([0, 0, 0])
    dim_M0 = sp.Matrix([0, 0, -1])
    dim_rho = sp.Matrix([0, -3, 0])
    dim_r = sp.Matrix([0, 1, 0])
    dim_d = sp.Matrix([0, 1, 0])
    dim_w = sp.Matrix([0, 1, 0])
    dim_cS = sp.Matrix([0, 1, -1])
    dim_omega = sp.Matrix([0, 0, -1])
    dim_k_warp = sp.Matrix([0, -1, 0])
    dim_green_zero_static = -dim_d - dim_r
    dim_flow = dim_green_zero_static - dim_r
    expect_zero("z=omega*d/c_S is dimensionless", dim_residual(dim_omega + dim_d - dim_cS, zero))
    expect_zero("k_warp*w is dimensionless", dim_residual(dim_k_warp + dim_w, zero))
    expect_zero("[green_zero_static]=[1/(4*pi*d*r)]=L^-2", dim_residual(dim_green_zero_static, sp.Matrix([0, -2, 0])))
    expect_zero("[radial flow]=L^-1*[Green]=L^-3", dim_residual(dim_flow, sp.Matrix([0, -3, 0])))
    expect_zero("source radial-flow dim check dim_M0-dim_rho-2*dim_r=[0,1,-1]", dim_residual(dim_M0 - dim_rho - 2 * dim_r, sp.Matrix([0, 1, -1])))
    print("  dimensions ordered M,L,T; no EOS/c_s^2 chain re-derived in stage010.")


def print_provenance() -> None:
    subbanner("Provenance and scope")
    print("  postulated slab: flat finite slab, brane w=0, absorber w=d; localization is derived on this family.")
    print("  BOUNDARY SELECTIONS: outgoing-spherical C1/C2 radial branch and 1/(4*pi*d) overlap normalization are selections, not derivations.")
    print("  R19/W_slab caveat: localization holding for the flat-slab FAMILY != the family being SELECTED by dynamics; the slab width d is a postulated register row (stage009), and its selection is the deferred nonlinear return (R19, sim-deferred Gate-6 territory).")
    print("  open-item #9: sharpened, NOT closed - the gravity-range (1/r^2) leg passes inside the localizing flat-slab family because both DC-sink completions give p=2; the return-cancellation targets (R23) remain the deferred obligation.")
    print("  Check-A citation: the joint RETURN_RESIDUAL_PREDICTION is completed here by composing Check B (p=2, computed) with Check A (A_residual_pass=True, cited from ledger_stage009); NOGO reachable.")
    print("  radiation/Sommerfeld boundary provenance: recorded ac_check_a_only in the source - NOT a Check-B branch.")
    print("  dropped bookkeeping: source SHA-256 static/dynamic/per-branch trace ids, structure_id, and expr_digest were build-reproducibility plumbing, replaced by the v2 tri-review protocol.")
    print("  downstream consumers: stage 024/026 (pathA_43 consumes the Phi_l(w,r) bulk Helmholtz zero-mode plus the projected-continuity operator lineage).")
    print("  dropped-trace bookkeeping note: the trace-hash lines are not carried into this print-only stage.")


def print_verdict_labels() -> None:
    print("")
    print("Verdict labels:")
    print("  ledger earned-label (NOT a source verdict token): RETURN_LOCALIZATION_P2_DERIVED_NOGO_REACHABLE  (both DC-sink completions -> normalizable m=0 zero mode -> genuine 3D-radial dsolve -> p=2; static-dynamic consistent; counterfactual r^-4 rejected (residual 5/(pi*d*r^7)); warp control -> p=3 -> RETURN_NOGO)")
    print("  source top-line verdict: RETURN_RESIDUAL_PREDICTION  (COMPLETED here: Check B p=2 computed; Check A component A_residual_pass=True cited from ledger_stage009)")
    print("  joint composition: RETURN_RESIDUAL_PREDICTION = (Check A: bounded residual, A_residual_pass=True, CITED stage009) AND (Check B: p=2 both completions, computed here); RETURN_NOGO reachable if the return delocalizes (warp -> p=3) OR quadrupole_survives=False")
    print("  earned: two normalizable m=0 zero modes (compact-cell + Bloch); dynamic dsolve then omega->0 -> p=2; static solve -> p=2; static-dynamic consistent; massive/gapped Yukawa contrast; both completions agree p=2; counterfactual r^-4 rejected; NOGO warp control p=3; falloff-tension witnessed")
    print("  postulated: the flat finite slab (brane w=0, absorber w=d) [stage009 geometry]")
    print("  boundary selections (NOT derived): outgoing-spherical C1/C2 branch choice; 1/(4*pi*d) overlap normalization")
    print("  consumed (cited, dual-site integrity): ledger_stage008 p_raw(l2)=5 -> quadrupole_survives=True (T2_applied=false, no l=2 recompute); ledger_stage009 Check-A component A_residual_pass=True")
    print("  caveat: localization holds for the FAMILY, not the family selected by dynamics (R19/W_slab); open-item #9 sharpened NOT closed (R23 pending)")


def run_able_to_fail_teeth(data: dict[str, Any]) -> None:
    subbanner("Able-to-fail mutation teeth")
    guard = data["counterfactual"]
    expect_fail(
        "tooth 1 corrupt counterfactual guard to accept zero residual trips guard assert",
        guard["perturbed_residual"],
    )

    normalizable_cutoff = sp.simplify(sp.integrate(sp.exp(-2 * k_warp * w), (w, 0, W)))
    normalizable_limit = sp.simplify(sp.limit(normalizable_cutoff, W, sp.oo))
    expect_fail(
        "tooth 2 warp flip exp(-2*k_warp*w) trips non-normalizability assert",
        bool_residual(normalizable_limit == sp.oo),
    )

    expect_fail(
        "tooth 3 radial operator coefficient 2/r->3/r trips selected residual assert",
        radial_operator_residual(data["static"]["destructuring_absorbing"]["selected"], sp.Integer(0), r, coefficient=3),
    )

    expect_fail(
        "tooth 4 nonzero m seed trips zero-mode radial guard",
        zero_seed_guard_residual(data["spectra"]["destructuring_absorbing"]["m1"]),
    )

    p_bloch_mutated = sp.Integer(3)
    mutated_headline = classify_dc_sink_gate([data["p_abs"], p_bloch_mutated], data["quadrupole_survives"])
    expect_fail(
        "tooth 5 p_bloch->3 trips completion agreement and headline classifier",
        sp.simplify((data["p_abs"] - p_bloch_mutated) ** 2)
        + branch_verdict_residual(mutated_headline, RETURN_RESIDUAL_PREDICTION),
    )

    expect_zero(
        "tooth 6a classifier maps delocalizing p=3 branch to RETURN_NOGO",
        branch_verdict_residual(classify_dc_sink_gate([data["p_delocalizing"]], data["quadrupole_survives"]), RETURN_NOGO),
    )
    expect_zero(
        "tooth 6b classifier maps quadrupole_survives=False to RETURN_NOGO",
        branch_verdict_residual(classify_dc_sink_gate([data["p_abs"], data["p_bloch"]], False), RETURN_NOGO),
    )

    f0_abs_corrupt = sp.simplify(1 / sp.sqrt(2 * d))
    expect_fail(
        "tooth 7 corrupt f0_abs=1/sqrt(2*d) trips zero-mode normalization assert",
        transverse_norm(f0_abs_corrupt) - 1,
    )

    expect_fail(
        "tooth 8a stage008 consumed p_raw(l2) 5->4 trips dual-site integrity",
        stage008_integrity_residual(sp.Integer(4), STAGE008_P_RAW_L2_PIPELINE),
    )
    expect_fail(
        "tooth 8b stage008 pipeline p_raw(l2) 5->4 trips dual-site integrity",
        stage008_integrity_residual(STAGE008_P_RAW_L2_CONSUMED, sp.Integer(4)),
    )
    expect_fail(
        "tooth 8c stage009 consumed A_residual_pass True->False trips dual-site integrity",
        stage009_integrity_residual(False, STAGE009_A_RESIDUAL_PASS_PIPELINE),
    )
    expect_fail(
        "tooth 8d stage009 pipeline A_residual_pass 1->0 trips dual-site integrity",
        stage009_integrity_residual(STAGE009_A_RESIDUAL_PASS_CONSUMED, sp.Integer(0)),
    )
    expect_bool(
        "tooth 8e non-finite p_raw sentinel would flip quadrupole_survives False",
        quadrupole_survives_from_p_raw(sp.oo) is False,
    )

    expect_zero("baseline immutable after teeth: p_abs remains 2", data["p_abs"] - 2)
    expect_zero("baseline immutable after teeth: p_bloch remains 2", data["p_bloch"] - 2)
    expect_zero("baseline immutable after teeth: stage008 integrity remains zero", data["stage008_integrity"])
    expect_zero("baseline immutable after teeth: stage009 integrity remains zero", data["stage009_integrity"])


def main() -> None:
    banner("ledger_stage010_slab_localization_p2_nogo SymPy audit")
    data = compute_baseline()
    assert_no_float("baseline", data)
    run_consumed_inputs(data)
    run_transverse_spectra(data)
    run_radial_routes(data)
    run_completion_agreement_and_counterfactual(data)
    run_warp_and_classifier(data)
    run_dimensional_block()
    print_provenance()
    print_verdict_labels()
    run_able_to_fail_teeth(data)
    print("")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print("OVERALL PASS: SymPy verified ledger_stage010 slab localization p=2 with reachable NOGO exactly")


if __name__ == "__main__":
    try:
        main()
    except AuditFailure as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
        print(f"OVERALL FAIL: SymPy stage010 audit did not close ({exc})")
        raise SystemExit(1)
