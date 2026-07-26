#!/usr/bin/env python3
"""Ledger stage011 SymPy audit: frozen-reduction certificate.

Standalone, print-only, no arguments, no file I/O.  This builds the Part-II
pathA_30 II-G1a slice only: the dimensionally-safe frozen longitudinal
operator assembly, c_S^2 extraction, solved pinch-off domain, reduction
certificate, c_S^2 dimensional leg, and able-to-fail teeth.  Stage 012 owns the
D/N ladder, DtN, Robin counterfactual, pole ladder, and static series.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from ledger_dimensions import Dimension, DimensionBasis, dim_residual


PASS_COUNT = 0
FAIL_COUNT = 0

REDUCTION_CERTIFIED = "REDUCTION_CERTIFIED"
FAIL_DIMENSIONAL = "FAIL_DIMENSIONAL"
FAIL_OPERATOR_INTRUSION = "FAIL_OPERATOR_INTRUSION"
FAIL_WRONG_SPEED = "FAIL_WRONG_SPEED"
FAIL_WRONG_DOMAIN = "FAIL_WRONG_DOMAIN"
DN_UNITTEST_FAIL_DIMENSIONAL = "DN_UNITTEST_FAIL_DIMENSIONAL"


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


def fmt(expr: Any) -> str:
    if isinstance(expr, bool):
        return "True" if expr else "False"
    if isinstance(expr, str):
        return expr
    return sp.sstr(sp.factor(sp.cancel(sp.simplify(expr))))


def assert_no_float(name: str, expr: Any) -> None:
    if isinstance(expr, Dimension):
        for label, value in expr.exponents.items():
            assert_no_float(f"{name}.{label}", value)
        return
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
    _record_fail(f"FAIL  {name}: residual = {fmt(clean)}")
    raise AuditFailure(f"{name} residual was not zero")


def expect_bool(name: str, condition: bool) -> None:
    expect_zero(name, sp.Integer(0) if condition else sp.Integer(1))


def expect_nonzero(name: str, residual: sp.Expr | int) -> None:
    assert_no_float(name, residual)
    clean = sp.simplify(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_pass(f"PASS  {name} is nonzero as required (residual = {fmt(clean)})")
        return
    _record_fail(f"FAIL  {name}: required nonzero residual vanished")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def expect_fail(name: str, residual: sp.Expr | int) -> None:
    assert_no_float(name, residual)
    clean = sp.simplify(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_pass(f"PASS  {name} produced required FAIL (residual = {fmt(clean)})")
        return
    _record_fail(f"FAIL  {name}: required mutation/ablation did not fire")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def expr_equal(lhs: sp.Expr | int, rhs: sp.Expr | int = 0) -> bool:
    return sp.simplify(lhs - rhs) == 0


def bool_residual(condition: bool) -> sp.Integer:
    return sp.Integer(0) if condition else sp.Integer(1)


def nonzero_assert_residual(expr: sp.Expr | int) -> sp.Integer:
    return sp.Integer(0) if sp.simplify(expr) != 0 else sp.Integer(1)


def verdict_residual(actual: str, expected: str) -> sp.Integer:
    return sp.Integer(0) if actual == expected else sp.Integer(1)


DIMENSION_BASIS = DimensionBasis("L", "M", "T", render="tuple")
Dim = DIMENSION_BASIS
ZERO_DIM = Dim()


def dim_of(expr: sp.Expr, dims: dict[sp.Symbol, Dimension]) -> Dimension:
    clean = sp.sympify(expr)
    if clean == 0 or clean.is_Number:
        return ZERO_DIM
    if isinstance(clean, sp.Symbol):
        if clean not in dims:
            raise AuditFailure(f"missing dimension for symbol {clean}")
        return dims[clean]
    if isinstance(clean, sp.Mul):
        total = ZERO_DIM
        for arg in clean.args:
            total = total * dim_of(arg, dims)
        return total
    if isinstance(clean, sp.Pow):
        base, power = clean.args
        if not power.is_number:
            raise AuditFailure(f"non-numeric power in dimension expression {clean}")
        return dim_of(base, dims) ** sp.Rational(power)
    if isinstance(clean, sp.Add):
        arg_dims = [dim_of(arg, dims) for arg in clean.args if sp.simplify(arg) != 0]
        if not arg_dims:
            return ZERO_DIM
        first = arg_dims[0]
        if any(arg_dim != first for arg_dim in arg_dims[1:]):
            raise AuditFailure(f"dimension mismatch in sum {clean}")
        return first
    raise AuditFailure(f"unsupported dimension expression {clean}")


s, L0, omega = sp.symbols("s L0 omega", positive=True, real=True)
K, rho, rho_star, m, hbar = sp.symbols("K rho rho_star m hbar", positive=True, real=True)
A_perp0, R_mouth = sp.symbols("A_perp0 R_mouth", positive=True, real=True)
ell_c = sp.Symbol("ell_c", positive=True, real=True)
xi_symbol = sp.Symbol("xi", positive=True, real=True)
delta_wall = sp.Symbol("delta_wall", nonzero=True, real=True)
epsilon_rho = sp.Symbol("epsilon_rho", nonzero=True, real=True)
psi_hat = sp.Function("psi_hat")
A_perp_mut = sp.Function("A_perp0_mut")

y = psi_hat(s)
cS_squared_bulk = sp.simplify(5 * K * rho_star**4 / m)
bulk_N = sp.simplify(omega**2 / cS_squared_bulk)
linear_taper = sp.simplify(R_mouth * (1 - s / L0))


def r1_site_from_exponent(exponent: int) -> sp.Expr:
    e = sp.Integer(exponent)
    return sp.simplify(e * K * rho ** (e - 1) / m)


def r1_eos_site_from_exponent(exponent: int) -> sp.Expr:
    return sp.simplify(sp.diff(K * rho ** exponent, rho) / m)


def compute_verdict(
    *,
    dimensional_ok: bool,
    unsuppressed_operator_intrusion: bool,
    operator_is_helmholtz: bool,
    speed_is_cs: bool,
    domain_is_L0: bool,
) -> str:
    if not dimensional_ok:
        return FAIL_DIMENSIONAL
    if unsuppressed_operator_intrusion:
        return FAIL_OPERATOR_INTRUSION
    if not operator_is_helmholtz:
        return FAIL_OPERATOR_INTRUSION
    if not speed_is_cs:
        return FAIL_WRONG_SPEED
    if not domain_is_L0:
        return FAIL_WRONG_DOMAIN
    return REDUCTION_CERTIFIED


def xi_ell_c_firewall_ok(xi_expr: sp.Expr) -> bool:
    return (xi_symbol != ell_c) and (not xi_expr.has(ell_c))


def solve_cap_endpoint(taper_expr: sp.Expr) -> sp.Expr:
    roots = sp.solve(sp.Eq(taper_expr, 0), s)
    if not roots:
        raise AuditFailure(f"taper has no symbolic cap root: {taper_expr}")
    return sp.simplify(roots[0])


def build_reduction_case(
    *,
    sqrt_gamma0_expr: sp.Expr,
    rho0_expr: sp.Expr,
    taper_expr: sp.Expr,
    bdg_flag: int,
    bdg_deferred_by_smallness: bool,
    delta_v_conf_expr: sp.Expr = sp.Integer(0),
    speed_factor: sp.Expr = sp.Integer(1),
) -> dict[str, Any]:
    M = sp.simplify(sp.diff(sp.log(sqrt_gamma0_expr), s))
    rho0grad = sp.simplify(sp.diff(rho0_expr, s))
    cs_squared_local = sp.simplify(5 * K * rho0_expr**4 / m)
    csgrad = sp.simplify(sp.diff(sp.sqrt(cs_squared_local), s))
    N = sp.simplify(omega**2 / (cs_squared_local * speed_factor**2))

    q_background = sp.simplify(
        -hbar**2 * sp.diff(sp.sqrt(rho0_expr), s, 2) / (2 * m * sp.sqrt(rho0_expr))
    )
    qgrad = sp.simplify(sp.diff(q_background, s))

    B = sp.Integer(bdg_flag)
    bdg_operator_coeff = sp.simplify(hbar**2 / (4 * m**2 * cS_squared_bulk))
    L_s = sp.simplify(
        sp.diff(y, s, 2)
        + M * sp.diff(y, s)
        + N * y
        - B * bdg_operator_coeff * sp.diff(y, s, 4)
    )
    ideal = sp.simplify(sp.diff(y, s, 2) + bulk_N * y)
    residual = sp.simplify(L_s - ideal)
    operator_is_helmholtz = expr_equal(residual, 0)

    scalar_intrusion = (not expr_equal(delta_v_conf_expr, 0)) or (not expr_equal(qgrad, 0))
    bdg_intrusion = bool(B != 0 and not bdg_deferred_by_smallness)
    unsuppressed_operator_intrusion = (
        (not expr_equal(M, 0))
        or (not expr_equal(N - bulk_N, 0))
        or bdg_intrusion
        or scalar_intrusion
    )

    psi_coeff = sp.simplify(sp.expand(L_s).coeff(y))
    extracted_speed_squared = sp.simplify(omega**2 / psi_coeff)
    speed_is_cs = (
        expr_equal(psi_coeff - bulk_N, 0)
        and expr_equal(extracted_speed_squared - cS_squared_bulk, 0)
        and expr_equal(csgrad, 0)
    )

    cap_endpoint = solve_cap_endpoint(taper_expr)
    domain_is_L0 = expr_equal(cap_endpoint - L0, 0)
    verdict = compute_verdict(
        dimensional_ok=True,
        unsuppressed_operator_intrusion=unsuppressed_operator_intrusion,
        operator_is_helmholtz=operator_is_helmholtz,
        speed_is_cs=speed_is_cs,
        domain_is_L0=domain_is_L0,
    )

    return {
        "sqrt_gamma0": sqrt_gamma0_expr,
        "rho0": rho0_expr,
        "M": M,
        "rho0grad": rho0grad,
        "cs_squared_local": cs_squared_local,
        "csgrad": csgrad,
        "N": N,
        "delta_v_conf": sp.simplify(delta_v_conf_expr),
        "Q": q_background,
        "Qgrad": qgrad,
        "B": B,
        "bdg_operator_coeff": bdg_operator_coeff,
        "L_s": L_s,
        "ideal": ideal,
        "operator_residual": residual,
        "operator_is_helmholtz": operator_is_helmholtz,
        "unsuppressed_operator_intrusion": unsuppressed_operator_intrusion,
        "psi_coeff": psi_coeff,
        "extracted_speed_squared": extracted_speed_squared,
        "speed_is_cs": speed_is_cs,
        "taper": taper_expr,
        "cap_endpoint": cap_endpoint,
        "domain_is_L0": domain_is_L0,
        "verdict": verdict,
    }


def build_dimensional_block() -> dict[str, Any]:
    length_dim = Dim(1, 0, 0)
    energy_dim = Dim(2, 1, -2)
    four_volume_dim = length_dim**4
    pressure_dim = energy_dim / four_volume_dim
    rho_dim = length_dim**-4
    K_dim = pressure_dim / (rho_dim**5)
    dim_rules = {
        L0: length_dim,
        omega: Dim(0, 0, -1),
        K: K_dim,
        rho_star: rho_dim,
        m: Dim(0, 1, 0),
    }
    expected_cs_squared_dim = Dim(2, 0, -2)
    cs_squared_dim = dim_of(cS_squared_bulk, dim_rules)
    dimensional_ok = cs_squared_dim == expected_cs_squared_dim

    corrupt_rules = dict(dim_rules)
    corrupt_rules[K] = corrupt_rules[K] * Dim(1, 0, 0)
    corrupt_cs_squared_dim = dim_of(cS_squared_bulk, corrupt_rules)
    corrupt_ok = corrupt_cs_squared_dim == expected_cs_squared_dim
    probe_verdict = "NO_FAIL" if corrupt_ok else DN_UNITTEST_FAIL_DIMENSIONAL
    mutation_fires = probe_verdict == DN_UNITTEST_FAIL_DIMENSIONAL

    clean_verdict = compute_verdict(
        dimensional_ok=dimensional_ok,
        unsuppressed_operator_intrusion=False,
        operator_is_helmholtz=True,
        speed_is_cs=True,
        domain_is_L0=True,
    )
    mutated_verdict = compute_verdict(
        dimensional_ok=corrupt_ok,
        unsuppressed_operator_intrusion=False,
        operator_is_helmholtz=True,
        speed_is_cs=True,
        domain_is_L0=True,
    )
    fail_suppressed = (
        clean_verdict == REDUCTION_CERTIFIED
        and mutated_verdict == FAIL_DIMENSIONAL
        and mutation_fires
    )

    return {
        "length_dim": length_dim,
        "energy_dim": energy_dim,
        "four_volume_dim": four_volume_dim,
        "pressure_dim": pressure_dim,
        "rho_dim": rho_dim,
        "K_dim": K_dim,
        "cs_squared_dim": cs_squared_dim,
        "expected_cs_squared_dim": expected_cs_squared_dim,
        "dimensional_ok": dimensional_ok,
        "corrupt_K_dim": corrupt_rules[K],
        "corrupt_cs_squared_dim": corrupt_cs_squared_dim,
        "corrupt_ok": corrupt_ok,
        "probe_verdict": probe_verdict,
        "mutation_fires": mutation_fires,
        "clean_verdict": clean_verdict,
        "mutated_verdict": mutated_verdict,
        "fail_suppressed": fail_suppressed,
    }


def build_baseline() -> dict[str, Any]:
    reduction = build_reduction_case(
        sqrt_gamma0_expr=A_perp0,
        rho0_expr=rho_star,
        taper_expr=linear_taper,
        bdg_flag=0,
        bdg_deferred_by_smallness=True,
        delta_v_conf_expr=sp.diff(sp.Integer(0), s),
    )
    k = sp.simplify(omega / sp.sqrt(cS_squared_bulk))
    xi = sp.simplify(hbar / (m * sp.sqrt(cS_squared_bulk)))
    bdg_ratio_direct = sp.simplify(hbar**2 * k**2 / (4 * m**2 * cS_squared_bulk))
    bdg_ratio_window = sp.simplify((k * xi / 2) ** 2)
    site_a = r1_site_from_exponent(5)
    site_b = r1_eos_site_from_exponent(5)
    consumed = sp.simplify(site_a.subs(rho, rho_star))
    dim = build_dimensional_block()
    verdict = compute_verdict(
        dimensional_ok=dim["dimensional_ok"],
        unsuppressed_operator_intrusion=reduction["unsuppressed_operator_intrusion"],
        operator_is_helmholtz=reduction["operator_is_helmholtz"],
        speed_is_cs=reduction["speed_is_cs"],
        domain_is_L0=reduction["domain_is_L0"],
    )
    reduction["verdict"] = verdict

    firewall_ok = xi_ell_c_firewall_ok(xi)
    return {
        "reduction": reduction,
        "k": k,
        "xi": xi,
        "bdg_ratio_direct": bdg_ratio_direct,
        "bdg_ratio_window": bdg_ratio_window,
        "site_a": site_a,
        "site_b": site_b,
        "consumed": consumed,
        "dim": dim,
        "verdict": verdict,
        "firewall_ok": firewall_ok,
    }


def run_reduction_certificate(data: dict[str, Any]) -> None:
    red = data["reduction"]
    subbanner("Frozen background and reduction certificate")
    print("  POSTULATED geometry: straight finite throat, frozen eta=0, brane/mouth s=0, cap s=L0.")
    print("  POSTULATED fields: rho0(s)=rho_star, A_M0=0, sqrt(gamma0)=A_perp0, matter perturbation exp(-I*omega*t).")
    print("  CITED speed: c_S^2 = 5*K*rho_star^4/m from R1 at rho_star; EOS exponent-5 P=K*rho^5 IMPOSED.")
    print(f"  M=d_s log sqrt(gamma0) = {fmt(red['M'])}")
    expect_zero("projection measure coefficient M is computed zero", red["M"])
    print(f"  rho0grad=d_s rho_star = {fmt(red['rho0grad'])}")
    expect_zero("rho0grad is computed zero", red["rho0grad"])
    print(f"  c_s(s)^2 = {fmt(red['cs_squared_local'])}")
    print(f"  N=(omega/c_s(s))^2 = {fmt(red['N'])}")
    expect_zero("frozen N equals (omega/c_S)^2 with R1 c_S^2", red["N"] - bulk_N)
    print(f"  csgrad=d_s sqrt(c_s^2) = {fmt(red['csgrad'])}")
    expect_zero("csgrad is computed zero", red["csgrad"])
    print(f"  delta_V_conf witness = {fmt(red['delta_v_conf'])} (ell_c inert in frozen eta=0 test)")
    expect_zero("delta_V_conf witness is computed zero", red["delta_v_conf"])
    print(f"  Q(rho0) = {fmt(red['Q'])}; Qgrad=d_s Q = {fmt(red['Qgrad'])}")
    expect_zero("Qgrad is computed by differentiation as zero", red["Qgrad"])
    print(f"  BdG ratio direct = {fmt(data['bdg_ratio_direct'])}")
    print(f"  BdG ratio window form (k*xi/2)^2 = {fmt(data['bdg_ratio_window'])}")
    expect_zero("BdG smallness witness hbar^2*k^2/(4*m^2*c_S^2) equals (k*xi/2)^2", data["bdg_ratio_direct"] - data["bdg_ratio_window"])
    expect_zero("baseline BdG inclusion flag B=0 under k*xi<<1 deferral", red["B"])


def run_operator_and_speed(data: dict[str, Any]) -> None:
    red = data["reduction"]
    subbanner("Produced operator, speed extraction, and domain solve")
    print(f"  assembled L_s = {fmt(red['L_s'])}")
    print(f"  ideal Helmholtz target = {fmt(red['ideal'])}")
    print(f"  residual L_s - ideal = {fmt(red['operator_residual'])}")
    expect_zero("operator_is_helmholtz is produced by L_s assembly", bool_residual(red["operator_is_helmholtz"]))
    expect_zero("assembled L_s minus ideal Helmholtz operator is zero", red["operator_residual"])
    expect_bool("unsuppressed_operator_intrusion is computed false", red["unsuppressed_operator_intrusion"] is False)
    print("  ODE artifact: psi''(s) + (omega/c_S)^2 psi(s) = 0")
    print(f"  extracted psi_hat coefficient = {fmt(red['psi_coeff'])}")
    print(f"  extracted speed squared omega^2/coeff = {fmt(red['extracted_speed_squared'])}")
    print(f"  bulk c_S^2 trace = {fmt(cS_squared_bulk)}")
    expect_zero("extracted psi_hat coefficient equals (omega/c_S)^2", red["psi_coeff"] - bulk_N)
    expect_zero("extracted speed squared equals bulk R1 c_S^2", red["extracted_speed_squared"] - cS_squared_bulk)
    expect_bool("speed_is_cs is extracted from L_s and csgrad", red["speed_is_cs"])
    print("  POSTULATED reference taper: R0(s)=R_mouth*(1-s/L0), monotone pinch with R0(0)>0 and R0(L0)=0.")
    print(f"  solved cap endpoint solve(R0(s)=0,s) = {fmt(red['cap_endpoint'])}")
    print(f"  domain = [0, {fmt(red['cap_endpoint'])}]")
    expect_zero("domain_is_L0 is solved from the taper root", bool_residual(red["domain_is_L0"]))
    expect_zero("cap endpoint minus L0 is zero", red["cap_endpoint"] - L0)


def run_consumed_r1(data: dict[str, Any]) -> None:
    subbanner("Consumed R1 input with dual-site integrity")
    print("  CITED ledger_stage005 R1 at rho_star; pathA_30 bare m is the stage004 m_GNLS ACTION primitive.")
    print(f"  site A literal = {fmt(data['site_a'])}")
    print(f"  site B EOS route d(K*rho^5)/d rho / m = {fmt(data['site_b'])}")
    expect_zero("R1 site A minus site B equals zero", data["site_a"] - data["site_b"])
    expect_zero("R1 evaluated at rho_star equals c_S^2 bulk", data["consumed"] - cS_squared_bulk)
    expect_zero("explicit frozen-export anchor consumed - 5*K*rho_star^4/m equals zero", data["consumed"] - 5 * K * rho_star**4 / m)


def run_dimensional_block(data: dict[str, Any]) -> None:
    dim = data["dim"]
    subbanner("Stage011 c_S^2 dimensional leg and corrupt-[K] probe")
    print("  dimension order: (L,M,T)")
    print(f"  [energy] = {dim['energy_dim']}; [four-volume] = {dim['four_volume_dim']}; [P] = {dim['pressure_dim']}")
    print(f"  [rho] = {dim['rho_dim']}; [K]=[P]-5[rho] = {dim['K_dim']}")
    print(f"  [c_S^2=5*K*rho_star^4/m] = {dim['cs_squared_dim']}")
    expect_zero("c_S^2 dimensional leg equals (2,0,-2)", dim_residual(dim["cs_squared_dim"], dim["expected_cs_squared_dim"]))
    expect_bool("dimensional_ok for the 011 c_S^2 leg", dim["dimensional_ok"])
    print(f"  corrupt [K]+(1,0,0) gives [K] = {dim['corrupt_K_dim']}")
    print(f"  corrupt [c_S^2] = {dim['corrupt_cs_squared_dim']} -> {dim['probe_verdict']}")
    expect_zero("corrupt-[K] mutated c_S^2 dimension is exactly (3,0,-2)", dim_residual(dim["corrupt_cs_squared_dim"], Dim(3, 0, -2)))
    expect_bool("corrupt-[K] mutation_fires=True", dim["mutation_fires"])
    expect_zero("self-ablation with mutation gives FAIL_DIMENSIONAL", verdict_residual(dim["mutated_verdict"], FAIL_DIMENSIONAL))
    expect_zero("self-ablation without mutation gives REDUCTION_CERTIFIED", verdict_residual(dim["clean_verdict"], REDUCTION_CERTIFIED))
    expect_bool("self-ablation fail_suppressed=True", dim["fail_suppressed"])


def run_firewall_and_verdict(data: dict[str, Any]) -> None:
    red = data["reduction"]
    subbanner("xi != ell_c firewall and 011 verdict")
    print(f"  xi = hbar/(m*c_s) = {fmt(data['xi'])}")
    print("  ell_c is the confinement length in V_wall(Sigma/ell_c), inert here because delta_V_conf=0.")
    expect_bool("xi and ell_c are distinct symbols and never substituted", data["firewall_ok"])
    print(f"  011 scoped verdict = {data['verdict']}")
    expect_zero("011 verdict is REDUCTION_CERTIFIED", verdict_residual(data["verdict"], REDUCTION_CERTIFIED))
    print(
        "  DN_UNITTEST_BC_DEPENDENT (JOINT) = "
        "(011: REDUCTION_CERTIFIED, computed here) AND "
        "(012: D/N ladder + bc_derivation_emitted=False -> BC_DEPENDENT landing, stage 012)"
    )
    expect_bool("operator/speed/domain booleans are all computed true in baseline", red["operator_is_helmholtz"] and red["speed_is_cs"] and red["domain_is_L0"])


def print_provenance() -> None:
    subbanner("Provenance and scope")
    print("  POSTULATED geometry: straight finite throat + frozen eta=0; L0 = ACTION-geometry throat depth, NOT a medium constant; R0(s) taper POSTULATED.")
    print("  CITED-speed: c_S^2=5*K*rho_star^4/m is Part I edge R1 (stage005) evaluated at rho_star; EOS exponent-5 IMPOSED.")
    print("  de-rig: operator PRODUCED by assembly, speed EXTRACTED from L_s, domain SOLVED from the pinch-off; replaces X==X/L0==L0/literal-True checks; the pinch-off domain is a labeled SELECTION.")
    print("  validity-window: L_s is const-coeff Helmholtz conditional on {rho0'/rho0=0, sqrt(gamma0) const, delta_V_conf=0, grad_Q=0, k*xi<<1}; BdG k^4 is DEFERRED, not dropped unconditionally.")
    print("  firewall: xi=hbar/(m*c_s) healing length is distinct from ell_c confinement length; ell_c is INERT here because delta_V_conf=0.")
    print("  split: 011 carries the reduction-certificate component; D/N ladder + Robin + BC_DEPENDENT landing are stage 012; bc_derivation_emitted=False is 012's rung.")
    print("  dropped-bookkeeping: scratch-YAML/_sympy_exprs.wl export, MMA-YAML re-read, expression_digest, and engine_agreement plumbing are stripped.")
    print("  downstream consumers: stage 013 (harmonic beta lift) + stage 017 (calibration input) consume frozen L_s, domain, c_S, and validity window.")
    print("  register note: zero new counted knobs; L0 is POSTULATED ACTION-geometry; ell_c INERT; xi DERIVED; validity record + firewall are structural edge candidates.")


def print_verdict_labels() -> None:
    print("")
    print("Verdict labels:")
    print("  ledger earned-label (NOT a source verdict token): FROZEN_REDUCTION_HELMHOLTZ_CERTIFIED  (frozen wall eta=0 -> L_s assembled from the reduction (projection measure + every intruding coeff computed with its vanishing/deferral condition) -> const-coeff Helmholtz psi''+(omega/c_S)^2 psi=0; c_S^2=5*K*rho_star^4/m bulk; domain [0,L0] solved from the pinch-off R0(L0)=0; validity window {rho0'/rho0=0, sqrt(gamma0) const, delta_V_conf=0, grad_Q=0, k*xi<<1})")
    print("  source top-line verdict: DN_UNITTEST_BC_DEPENDENT  (JOINT; 011 carries the reduction-certificate component REDUCTION_CERTIFIED; the D/N ladder + bc_derivation_emitted=False -> BC_DEPENDENT landing = stage 012)")
    print("  joint composition: DN_UNITTEST_BC_DEPENDENT = (011: REDUCTION_CERTIFIED, computed here) AND (012: D/N pole ladder + Robin + BC_DEPENDENT landing, stage 012)")
    print("  earned (de-rig): operator_is_helmholtz PRODUCED by assembly (not X==X); speed_is_cs EXTRACTED from L_s (not literal True); domain_is_L0 SOLVED from R0(L0)=0 pinch-off (not L0==L0); unsuppressed_operator_intrusion COMPUTED; c_S^2 dim leg (2,0,-2) via [K]=[P]-5[rho] + corrupt-[K] probe fires")
    print("  postulated: straight finite throat, brane s=0, cap s=L0 (R0(L0)=0), rho0=rho_star, A_M0=0, frozen eta=0; L0 = ACTION-geometry throat depth (NOT a medium constant); R0(s) taper postulated")
    print("  cited (R1, stage005, dual-site integrity): c_S^2 = 5*K*rho^4/m evaluated at rho_star; EOS exponent-5 P=K*rho^5 IMPOSED")
    print("  deferred (validity window): BdG k^4 term hbar^2*k^4/(4*m^2), ratio hbar^2*k^2/(4*m^2*c_S^2)=(k*xi/2)^2, deferred only under k*xi<<1")
    print("  firewall: xi=hbar/(m*c_s) (healing length) != ell_c (confinement length) -- distinct symbols; ell_c INERT here (delta_V_conf=0)")


def run_able_to_fail_teeth(data: dict[str, Any]) -> None:
    subbanner("Able-to-fail mutation teeth")
    baseline = data["reduction"]

    measure_mut = build_reduction_case(
        sqrt_gamma0_expr=A_perp_mut(s),
        rho0_expr=rho_star,
        taper_expr=linear_taper,
        bdg_flag=0,
        bdg_deferred_by_smallness=True,
    )
    expect_fail("tooth 1 nonconstant sqrt(gamma0) makes M nonzero", measure_mut["M"])
    expect_fail("tooth 1 operator_is_helmholtz boolean flips false", bool_residual(measure_mut["operator_is_helmholtz"]))
    expect_zero("tooth 1 verdict is FAIL_OPERATOR_INTRUSION", verdict_residual(measure_mut["verdict"], FAIL_OPERATOR_INTRUSION))

    bdg_mut = build_reduction_case(
        sqrt_gamma0_expr=A_perp0,
        rho0_expr=rho_star,
        taper_expr=linear_taper,
        bdg_flag=1,
        bdg_deferred_by_smallness=False,
    )
    bdg_fourth_derivative = sp.diff(y, s, 4)
    baseline_bdg_fourth_derivative_coeff = sp.simplify(
        sp.expand(baseline["L_s"]).coeff(bdg_fourth_derivative)
    )
    retained_bdg_fourth_derivative_coeff = sp.simplify(
        sp.expand(bdg_mut["L_s"]).coeff(bdg_fourth_derivative)
    )
    expect_zero(
        "tooth 2 deferred BdG flag leaves fourth-derivative term absent in baseline",
        baseline_bdg_fourth_derivative_coeff,
    )
    expect_nonzero(
        "tooth 2 retained BdG flag injects fourth-derivative term into L_s",
        retained_bdg_fourth_derivative_coeff,
    )
    expect_fail("tooth 2 operator_is_helmholtz boolean flips false", bool_residual(bdg_mut["operator_is_helmholtz"]))
    expect_zero("tooth 2 verdict is FAIL_OPERATOR_INTRUSION", verdict_residual(bdg_mut["verdict"], FAIL_OPERATOR_INTRUSION))

    nonuniform_rho0 = sp.simplify(rho_star * (1 + epsilon_rho * s / L0))
    rho_mut = build_reduction_case(
        sqrt_gamma0_expr=A_perp0,
        rho0_expr=nonuniform_rho0,
        taper_expr=linear_taper,
        bdg_flag=0,
        bdg_deferred_by_smallness=True,
    )
    expect_fail("tooth 3 nonuniform rho0 makes N-(omega/c_S)^2 nonzero", rho_mut["N"] - bulk_N)
    expect_fail("tooth 3 operator_is_helmholtz boolean flips false", bool_residual(rho_mut["operator_is_helmholtz"]))
    expect_zero("tooth 3 verdict is FAIL_OPERATOR_INTRUSION", verdict_residual(rho_mut["verdict"], FAIL_OPERATOR_INTRUSION))

    delta_v_conf_mut = build_reduction_case(
        sqrt_gamma0_expr=A_perp0,
        rho0_expr=rho_star,
        taper_expr=linear_taper,
        bdg_flag=0,
        bdg_deferred_by_smallness=True,
        delta_v_conf_expr=delta_wall,
    )
    expect_zero(
        "tooth 3b verdict is FAIL_OPERATOR_INTRUSION",
        verdict_residual(delta_v_conf_mut["verdict"], FAIL_OPERATOR_INTRUSION),
    )
    expect_zero(
        "tooth 3b nonzero delta_V_conf witness makes unsuppressed_operator_intrusion true",
        bool_residual(delta_v_conf_mut["unsuppressed_operator_intrusion"]),
    )
    expect_zero(
        "tooth 3b operator_is_helmholtz stays true because the witness is not in L_s",
        bool_residual(delta_v_conf_mut["operator_is_helmholtz"]),
    )

    print("  note: FAIL_WRONG_SPEED is a defensive verdict branch not reachable by real operator corruption; intrusion dominates.")
    speed_mut_verdict = compute_verdict(
        dimensional_ok=True,
        unsuppressed_operator_intrusion=baseline["unsuppressed_operator_intrusion"],
        operator_is_helmholtz=baseline["operator_is_helmholtz"],
        speed_is_cs=False,
        domain_is_L0=baseline["domain_is_L0"],
    )
    expect_zero(
        "compute_verdict logic branch returns FAIL_WRONG_SPEED when only speed_is_cs is false",
        verdict_residual(speed_mut_verdict, FAIL_WRONG_SPEED),
    )

    taper_mut = sp.simplify(R_mouth * (1 - s / (2 * L0)))
    domain_mut = build_reduction_case(
        sqrt_gamma0_expr=A_perp0,
        rho0_expr=rho_star,
        taper_expr=taper_mut,
        bdg_flag=0,
        bdg_deferred_by_smallness=True,
    )
    expect_fail("tooth 5 corrupt taper root differs from L0", domain_mut["cap_endpoint"] - L0)
    expect_fail("tooth 5 domain_is_L0 boolean flips false", bool_residual(domain_mut["domain_is_L0"]))
    expect_zero("tooth 5 verdict is FAIL_WRONG_DOMAIN", verdict_residual(domain_mut["verdict"], FAIL_WRONG_DOMAIN))

    xi_conflated = sp.simplify(data["xi"].subs(hbar, ell_c * hbar))
    conflated_firewall_ok = xi_ell_c_firewall_ok(xi_conflated)
    expect_fail("tooth 6 ell_c -> xi conflation trips distinct-symbol firewall", bool_residual(conflated_firewall_ok))

    for exponent in (4, 6):
        expect_fail(
            f"tooth 7 site A exponent 5->{exponent} trips R1 dual-site integrity",
            r1_site_from_exponent(exponent) - data["site_b"],
        )
        expect_fail(
            f"tooth 7 site B exponent 5->{exponent} trips R1 dual-site integrity",
            data["site_a"] - r1_eos_site_from_exponent(exponent),
        )
    expect_fail(
        "tooth 7 coordinated both-site exponent drift still trips frozen-export anchor",
        r1_site_from_exponent(6).subs(rho, rho_star) - 5 * K * rho_star**4 / m,
    )

    dim = data["dim"]
    expect_fail(
        "tooth 8 corrupt-[K] probe trips dimensional gate",
        dim_residual(dim["corrupt_cs_squared_dim"], dim["expected_cs_squared_dim"]),
    )
    expect_zero("tooth 8 corrupt-[K] verdict is FAIL_DIMENSIONAL", verdict_residual(dim["mutated_verdict"], FAIL_DIMENSIONAL))
    expect_bool("tooth 8 self-ablation fail_suppressed remains true", dim["fail_suppressed"])

    expect_zero("baseline immutable after teeth: operator residual remains zero", baseline["operator_residual"])
    expect_zero("baseline immutable after teeth: speed extraction remains bulk", baseline["extracted_speed_squared"] - cS_squared_bulk)
    expect_zero("baseline immutable after teeth: cap endpoint remains L0", baseline["cap_endpoint"] - L0)
    expect_zero("baseline immutable after teeth: R1 site integrity remains zero", data["site_a"] - data["site_b"])
    expect_zero("baseline immutable after teeth: clean 011 verdict remains REDUCTION_CERTIFIED", verdict_residual(data["verdict"], REDUCTION_CERTIFIED))


def main() -> None:
    banner("ledger_stage011_frozen_reduction_certificate SymPy audit")
    data = build_baseline()
    assert_no_float("baseline", data)
    run_reduction_certificate(data)
    run_operator_and_speed(data)
    run_consumed_r1(data)
    run_dimensional_block(data)
    run_firewall_and_verdict(data)
    print_provenance()
    print_verdict_labels()
    run_able_to_fail_teeth(data)
    print("")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print("OVERALL PASS: SymPy verified ledger_stage011 frozen reduction certificate exactly")


if __name__ == "__main__":
    try:
        main()
    except AuditFailure as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
        print(f"OVERALL FAIL: SymPy stage011 audit did not close ({exc})")
        raise SystemExit(1)
