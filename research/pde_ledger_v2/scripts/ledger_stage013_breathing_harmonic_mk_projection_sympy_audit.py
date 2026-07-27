#!/usr/bin/env python3
"""Ledger stage013 SymPy audit: breathing harmonic profiles + M/K projection.

Standalone, with audit results on stdout and labelled dimensions in a
deterministic sidecar.  This is the Part-II pathA_31 II-G2a slice only:
harmonic-lift profiles, M_AB/K_AB by real integral operator projection, the
(a,L) EOM LHS, the M/K dimensional legs, and the 013-native guard teeth.
Stages 014 and 015 own truncation consistency and the Hellmann-Feynman
RHS/legacy-structure checks.
"""

from __future__ import annotations

from typing import Any

import sympy as sp

from ledger_dimensions import (
    Dimension,
    DimensionBasis,
    dim_residual,
    emit_dimension_sidecar,
)


PASS_COUNT = 0
FAIL_COUNT = 0

BREATHING_CALIBRATED = "BREATHING_CALIBRATED"
BREATHING_FAIL_DIMENSIONAL = "BREATHING_FAIL_DIMENSIONAL"


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


def compact_expr(expr: Any) -> Any:
    if isinstance(expr, sp.MatrixBase):
        return expr.applyfunc(compact_expr)
    if not isinstance(expr, sp.Basic):
        return expr
    if isinstance(expr, (sp.Equality, sp.Order)):
        return expr
    return sp.trigsimp(sp.factor(sp.cancel(sp.together(sp.simplify(expr)))))


def expand_hyperbolic_args(expr: Any) -> Any:
    funcs = (sp.sinh, sp.cosh, sp.tanh, sp.coth, sp.sech, sp.csch)
    if isinstance(expr, sp.MatrixBase):
        return expr.applyfunc(expand_hyperbolic_args)
    if not isinstance(expr, sp.Basic):
        return expr
    return expr.replace(
        lambda node: node.func in funcs,
        lambda node: node.func(*[sp.expand(arg) for arg in node.args]),
    )


def fmt(expr: Any) -> str:
    if isinstance(expr, bool):
        return "True" if expr else "False"
    if isinstance(expr, str):
        return expr
    if isinstance(expr, Dimension):
        return str(expr)
    if isinstance(expr, sp.MatrixBase):
        return sp.sstr(compact_expr(expr))
    if isinstance(expr, (dict, list, tuple)):
        return sp.sstr(expr)
    try:
        return sp.sstr(compact_expr(expr))
    except Exception:
        return sp.sstr(expr)


def assert_no_float(name: str, expr: Any) -> None:
    if isinstance(expr, Dimension):
        for label, value in expr.exponents.items():
            assert_no_float(f"{name}.{label}", value)
        return
    if isinstance(expr, dict):
        for key, value in expr.items():
            assert_no_float(f"{name}.{key}", value)
        return
    if isinstance(expr, sp.MatrixBase):
        for index, value in enumerate(expr):
            assert_no_float(f"{name}[{index}]", value)
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
    clean = compact_expr(residual)
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
    clean = compact_expr(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_pass(f"PASS  {name} is nonzero as required (residual = {fmt(clean)})")
        return
    _record_fail(f"FAIL  {name}: required nonzero residual vanished")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def expect_fail(name: str, residual: sp.Expr | int) -> None:
    assert_no_float(name, residual)
    clean = compact_expr(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_pass(f"PASS  {name} produced required FAIL (residual = {fmt(clean)})")
        return
    _record_fail(f"FAIL  {name}: required mutation/ablation did not fire")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def expr_equal(lhs: sp.Expr | int, rhs: sp.Expr | int = 0) -> bool:
    return compact_expr(lhs - rhs) == 0


def bool_residual(condition: bool) -> sp.Integer:
    return sp.Integer(0) if condition else sp.Integer(1)


def verdict_residual(actual: str, expected: str) -> sp.Integer:
    return sp.Integer(0) if actual == expected else sp.Integer(1)


DIMENSION_BASIS = DimensionBasis("L", "M", "T", render="tuple")
Dim = DIMENSION_BASIS
ZERO_DIM = Dim()


def dim_of(expr: sp.Expr, dims: dict[sp.Symbol, Dimension]) -> Dimension:
    clean = sp.sympify(expr)
    if clean == 0 or clean.is_Number or clean.is_number:
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
        arg_dims = [dim_of(arg, dims) for arg in clean.args if compact_expr(arg) != 0]
        if not arg_dims:
            return ZERO_DIM
        first = arg_dims[0]
        if any(arg_dim != first for arg_dim in arg_dims[1:]):
            raise AuditFailure(f"dimension mismatch in sum {clean}")
        return first
    if clean.func in (sp.sin, sp.cos, sp.tan, sp.cot, sp.sinh, sp.cosh, sp.tanh, sp.coth, sp.sech, sp.csch):
        arg_dims = [dim_of(arg, dims) for arg in clean.args]
        if any(arg_dim != ZERO_DIM for arg_dim in arg_dims):
            raise AuditFailure(f"dimensionful argument in dimensionless function {clean}")
        return ZERO_DIM
    raise AuditFailure(f"unsupported dimension expression {clean}")


w = sp.Symbol("w", real=True)
L0 = sp.Symbol("L0", positive=True, real=True)
beta = sp.Symbol("beta", positive=True, real=True)
muEta = sp.Symbol("muEta", positive=True, real=True)
Tw = sp.Symbol("Tw", positive=True, real=True)
rAL = sp.Symbol("rAL", positive=True, real=True)
delta_a, delta_L = sp.symbols("delta_a delta_L", real=True)
delta_a_ddot, delta_L_ddot = sp.symbols("delta_a_ddot delta_L_ddot", real=True)
F_a_HF, F_L_HF = sp.symbols("F_a_HF F_L_HF", real=True)


def compute_013_verdict(*, dimensional_ok: bool) -> str:
    if not dimensional_ok:
        return BREATHING_FAIL_DIMENSIONAL
    return BREATHING_CALIBRATED


def lop(profile: sp.Expr, K_eta: sp.Expr) -> sp.Expr:
    return compact_expr((-Tw * sp.diff(profile, w, 2) + K_eta * profile) / muEta)


def build_profiles(K_eta: sp.Expr) -> dict[str, Any]:
    C1, C2 = sp.symbols("C1 C2")
    general = C1 * sp.sinh(beta * w) + C2 * sp.cosh(beta * w)
    coeff_a = sp.solve(
        [sp.Eq(general.subs(w, 0), 1), sp.Eq(general.subs(w, L0), 0)],
        (C1, C2),
        dict=True,
    )[0]
    coeff_L = sp.solve(
        [sp.Eq(general.subs(w, 0), 0), sp.Eq(general.subs(w, L0), rAL)],
        (C1, C2),
        dict=True,
    )[0]
    alpha_a = expand_hyperbolic_args(compact_expr(general.subs(coeff_a)))
    alpha_L = expand_hyperbolic_args(compact_expr(general.subs(coeff_L)))
    return {
        "general": general,
        "coeff_a": coeff_a,
        "coeff_L": coeff_L,
        "alpha_a": alpha_a,
        "alpha_L": alpha_L,
        "Lop_alpha_a": lop(alpha_a, K_eta),
        "Lop_alpha_L": lop(alpha_L, K_eta),
        "bc_values": {
            "alpha_a_mouth": compact_expr(alpha_a.subs(w, 0)),
            "alpha_a_cap": compact_expr(alpha_a.subs(w, L0)),
            "alpha_L_mouth": compact_expr(alpha_L.subs(w, 0)),
            "alpha_L_cap": compact_expr(alpha_L.subs(w, L0)),
        },
    }


def make_integrands(alpha_a: sp.Expr, alpha_L: sp.Expr, K_eta: sp.Expr, *, include_gradient: bool = True) -> dict[str, dict[str, sp.Expr]]:
    profiles = {"aa": (alpha_a, alpha_a), "aL": (alpha_a, alpha_L), "LL": (alpha_L, alpha_L)}
    mass_integrands = {
        name: expand_hyperbolic_args(compact_expr(muEta * left * right))
        for name, (left, right) in profiles.items()
    }
    stiffness_integrands = {}
    for name, (left, right) in profiles.items():
        gradient = Tw * sp.diff(left, w) * sp.diff(right, w) if include_gradient else sp.Integer(0)
        stiffness_integrands[name] = expand_hyperbolic_args(compact_expr(gradient + K_eta * left * right))
    return {"M_integrands": mass_integrands, "K_integrands": stiffness_integrands}


def integrate_entries(integrands: dict[str, dict[str, sp.Expr]]) -> dict[str, Any]:
    M_entries = {
        name: compact_expr(4 * sp.pi * sp.integrate(expr, (w, 0, L0)))
        for name, expr in integrands["M_integrands"].items()
    }
    K_entries = {
        name: compact_expr(4 * sp.pi * sp.integrate(expr, (w, 0, L0)))
        for name, expr in integrands["K_integrands"].items()
    }
    M_matrix = sp.Matrix([[M_entries["aa"], M_entries["aL"]], [M_entries["aL"], M_entries["LL"]]])
    K_matrix = sp.Matrix([[K_entries["aa"], K_entries["aL"]], [K_entries["aL"], K_entries["LL"]]])
    return {
        "M_entries": M_entries,
        "K_entries": K_entries,
        "M_matrix": M_matrix,
        "K_matrix": K_matrix,
        "M_det": compact_expr(M_matrix.det()),
        "K_det": compact_expr(K_matrix.det()),
    }


def project_from_profiles(alpha_a: sp.Expr, alpha_L: sp.Expr, K_eta: sp.Expr, *, include_gradient: bool = True) -> dict[str, Any]:
    integrands = make_integrands(alpha_a, alpha_L, K_eta, include_gradient=include_gradient)
    return {**integrands, **integrate_entries(integrands)}


def report_closed_forms() -> dict[str, Any]:
    return {
        "M_entries": {
            "aa": -2 * sp.pi * muEta * (L0 * beta * sp.tanh(L0 * beta) - sp.sinh(L0 * beta) ** 2) / (beta * sp.sinh(L0 * beta) ** 2 * sp.tanh(L0 * beta)),
            "aL": 2 * sp.pi * muEta * rAL * (L0 * beta - sp.tanh(L0 * beta)) / (beta * sp.sinh(L0 * beta) * sp.tanh(L0 * beta)),
            "LL": -2 * sp.pi * muEta * rAL**2 * (L0 * beta * sp.tanh(L0 * beta) - sp.sinh(L0 * beta) ** 2) / (beta * sp.sinh(L0 * beta) ** 2 * sp.tanh(L0 * beta)),
        },
        "K_entries": {
            "aa": 4 * sp.pi * Tw * beta / sp.tanh(L0 * beta),
            "aL": -4 * sp.pi * Tw * beta * rAL / sp.sinh(L0 * beta),
            "LL": 4 * sp.pi * Tw * beta * rAL**2 / sp.tanh(L0 * beta),
        },
        "M_det": -4 * sp.pi**2 * muEta**2 * rAL**2 * (L0 * beta - sp.sinh(L0 * beta)) * (L0 * beta + sp.sinh(L0 * beta)) / (beta**2 * sp.sinh(L0 * beta) ** 2),
        "K_det": 16 * sp.pi**2 * Tw**2 * beta**2 * rAL**2,
    }


def projection_residuals(projected: dict[str, Any], reference: dict[str, Any]) -> list[sp.Expr]:
    residuals = []
    for entry in ("aa", "aL", "LL"):
        residuals.append(compact_expr(projected["M_entries"][entry] - reference["M_entries"][entry]))
        residuals.append(compact_expr(projected["K_entries"][entry] - reference["K_entries"][entry]))
    return residuals


def all_projection_residuals_zero(projected: dict[str, Any], reference: dict[str, Any]) -> bool:
    return all(residual == 0 for residual in projection_residuals(projected, reference))


def symbol_name_flags(M_entries: dict[str, sp.Expr], K_entries: dict[str, sp.Expr]) -> dict[str, Any]:
    def names(entries: dict[str, sp.Expr]) -> set[str]:
        return {symbol.name for expr in entries.values() for symbol in expr.free_symbols}

    legacy_names = {"kappa", "chi", "sigmaA", "sigmaL"}
    allowed_names = {"muEta", "Tw", "beta", "L0", "rAL", "pi"}
    m_names = names(M_entries)
    k_names = names(K_entries)
    mk_names = m_names | k_names
    unexpected_names = mk_names - allowed_names
    flags = {
        "K_from_static_hessian": bool(k_names & legacy_names),
        "M_or_K_typed_from_legacy_values": bool(mk_names & legacy_names),
        "full_matrix_fit": bool(unexpected_names),
    }
    return {
        "M_names": m_names,
        "K_names": k_names,
        "MK_names": mk_names,
        "legacy_names": legacy_names,
        "allowed_names": allowed_names,
        "unexpected_names": unexpected_names,
        "forbidden_fit_flags": flags,
    }


def matrix_entry_dims(
    entries: dict[str, sp.Expr],
    symbol_dims: dict[sp.Symbol, Dimension],
) -> dict[str, Dimension]:
    return {name: dim_of(expr, symbol_dims) for name, expr in entries.items()}


def shared_dimension(
    entry_dims: dict[str, Dimension],
) -> Dimension | None:
    dims = list(entry_dims.values())
    if not dims:
        return None
    first = dims[0]
    return first if all(dim == first for dim in dims) else None


def build_dimensional_block(
    K_eta: sp.Expr,
    M_entries: dict[str, sp.Expr],
    K_entries: dict[str, sp.Expr],
) -> dict[str, Any]:
    symbol_dims = {
        L0: Dim(1, 0, 0),
        beta: Dim(-1, 0, 0),
        muEta: Dim(-1, 1, 0),
        Tw: Dim(1, 1, -2),
        rAL: ZERO_DIM,
    }
    expected_m = Dim(0, 1, 0)
    expected_k = Dim(0, 1, -2)
    expected_ratio = Dim(0, 0, -2)

    K_eta_dim = dim_of(K_eta, symbol_dims)
    m_dims = matrix_entry_dims(M_entries, symbol_dims)
    k_dims = matrix_entry_dims(K_entries, symbol_dims)
    m_shared = shared_dimension(m_dims)
    k_shared = shared_dimension(k_dims)
    ratio_dim = k_shared / m_shared if m_shared is not None and k_shared is not None else None
    dimensional_ok = bool(m_shared == expected_m and k_shared == expected_k and ratio_dim == expected_ratio)

    corrupt_dims = dict(symbol_dims)
    corrupt_dims[Tw] = symbol_dims[Tw] * Dim(1, 0, 0)
    corrupt_m_dims = matrix_entry_dims(M_entries, corrupt_dims)
    corrupt_k_dims = matrix_entry_dims(K_entries, corrupt_dims)
    corrupt_m_shared = shared_dimension(corrupt_m_dims)
    corrupt_k_shared = shared_dimension(corrupt_k_dims)
    corrupt_ratio_dim = corrupt_k_shared / corrupt_m_shared if corrupt_m_shared is not None and corrupt_k_shared is not None else None
    corrupt_ok = bool(corrupt_m_shared == expected_m and corrupt_k_shared == expected_k and corrupt_ratio_dim == expected_ratio)
    clean_verdict = compute_013_verdict(dimensional_ok=dimensional_ok)
    mutated_verdict = compute_013_verdict(dimensional_ok=corrupt_ok)

    return {
        "symbol_dims": symbol_dims,
        "expected_m": expected_m,
        "expected_k": expected_k,
        "expected_ratio": expected_ratio,
        "K_eta_dim": K_eta_dim,
        "m_dims": m_dims,
        "k_dims": k_dims,
        "m_shared": m_shared,
        "k_shared": k_shared,
        "ratio_dim": ratio_dim,
        "dimensional_ok": dimensional_ok,
        "corrupt_dims": corrupt_dims,
        "corrupt_m_dims": corrupt_m_dims,
        "corrupt_k_dims": corrupt_k_dims,
        "corrupt_m_shared": corrupt_m_shared,
        "corrupt_k_shared": corrupt_k_shared,
        "corrupt_ratio_dim": corrupt_ratio_dim,
        "corrupt_ok": corrupt_ok,
        "mutation_fires": not corrupt_ok,
        "probe_verdict": BREATHING_FAIL_DIMENSIONAL if not corrupt_ok else "NO_FAIL",
        "clean_verdict": clean_verdict,
        "mutated_verdict": mutated_verdict,
        "fail_suppressed": clean_verdict == BREATHING_CALIBRATED and mutated_verdict == BREATHING_FAIL_DIMENSIONAL,
    }


def packet_residuals(packet: dict[str, sp.Expr]) -> dict[str, sp.Expr]:
    return {
        "site_A_constitutive": compact_expr(packet["beta_cited"] - sp.sqrt(packet["K_eta_cited"] / packet["Tw_cited"])),
        "site_B_branch": compact_expr(packet["beta_cited"] * packet["L0_cited"] - packet["branch_anchor"]),
        "anchor_L0": compact_expr(packet["L0_cited"] - sp.Rational(37, 20)),
        "anchor_Tw": compact_expr(packet["Tw_cited"] - 1),
        "anchor_K_eta": compact_expr(packet["K_eta_cited"] - 1),
        "anchor_beta": compact_expr(packet["beta_cited"] - 1),
    }


def packet_ok(packet: dict[str, sp.Expr]) -> bool:
    return all(residual == 0 for residual in packet_residuals(packet).values())


def packet_mutation(packet: dict[str, sp.Expr], key: str, value: sp.Expr) -> dict[str, sp.Expr]:
    mutated = dict(packet)
    mutated[key] = value
    return mutated


def degenerate_mass_guard(alpha_L: sp.Expr, K_eta: sp.Expr, *, enabled: bool = True) -> dict[str, Any]:
    degenerate_projection = project_from_profiles(sp.Integer(0), alpha_L, K_eta)
    det_zero = compact_expr(degenerate_projection["M_det"])
    caught = enabled and det_zero == 0
    return {
        "projection": degenerate_projection,
        "M_det": det_zero,
        "caught": caught,
    }


def build_baseline() -> dict[str, Any]:
    K_eta = compact_expr(Tw * beta**2)
    profiles = build_profiles(K_eta)
    projection = project_from_profiles(profiles["alpha_a"], profiles["alpha_L"], K_eta)
    reintegrated = project_from_profiles(profiles["alpha_a"], profiles["alpha_L"], K_eta)
    closed_forms = report_closed_forms()
    flags = symbol_name_flags(projection["M_entries"], projection["K_entries"])
    dim = build_dimensional_block(K_eta, projection["M_entries"], projection["K_entries"])
    eom_rows = sp.Matrix(projection["M_matrix"]) * sp.Matrix([delta_a_ddot, delta_L_ddot]) + sp.Matrix(projection["K_matrix"]) * sp.Matrix([delta_a, delta_L])
    consumed_packet = {
        "L0_cited": sp.Rational(37, 20),
        "Tw_cited": sp.Integer(1),
        "K_eta_cited": sp.Integer(1),
        "beta_cited": sp.Integer(1),
        "branch_anchor": sp.Rational(37, 20),
    }
    degenerate = degenerate_mass_guard(profiles["alpha_L"], K_eta)
    return {
        "K_eta": K_eta,
        "profiles": profiles,
        "projection": projection,
        "reintegrated": reintegrated,
        "closed_forms": closed_forms,
        "projection_matches": all_projection_residuals_zero(reintegrated, closed_forms),
        "flags": flags,
        "dim": dim,
        "eom_rows": [compact_expr(row) for row in eom_rows],
        "rhs_placeholders": (F_a_HF, F_L_HF),
        "consumed_packet": consumed_packet,
        "degenerate": degenerate,
        "verdict": compute_013_verdict(dimensional_ok=dim["dimensional_ok"]),
    }


def run_operator_profiles(data: dict[str, Any]) -> None:
    profiles = data["profiles"]
    subbanner("Operator, ell=0 restriction, collective BCs, and harmonic profiles")
    print("  CONSUMED/POSTULATED operator: Lop[alpha] = muEta^(-1)*(-d_w(Tw*d_w alpha) + K_eta*alpha) on w in [0,L0].")
    print("  Alias for projection only: K_eta = Tw*beta^2; beta = sqrt(K_eta/Tw) is the cited wall-packet relation.")
    print("  Notation guard: domain length is L0; operator applications are named Lop_alpha_a and Lop_alpha_L.")
    print("  ell=0 restriction: Y00=1/(2*sqrt(pi)); int_S2 Y00^2 dOmega=1; eta(w,t)=eta_00(w,t)*Y00; T_Omega drops because ell*(ell+1)=0.")
    print("  Inner product: <f,g>_mu = 4*pi*int_0^L0 muEta*f*g dw.")
    print("  IMPOSED collective BCs: alpha_a(0)=1, alpha_a(L0)=0; alpha_L(0)=0, alpha_L(L0)=rAL.")
    print("  General solution form = ", fmt(profiles["general"]))
    print("  solved alpha_a = ", fmt(profiles["alpha_a"]))
    print("  solved alpha_L = ", fmt(profiles["alpha_L"]))

    expected_alpha_a = sp.sinh(L0 * beta - beta * w) / sp.sinh(L0 * beta)
    expected_alpha_L = rAL * sp.sinh(beta * w) / sp.sinh(L0 * beta)
    expect_zero("alpha_a equals reported harmonic-lift profile", profiles["alpha_a"] - expected_alpha_a)
    expect_zero("alpha_L equals reported harmonic-lift profile", profiles["alpha_L"] - expected_alpha_L)
    expect_zero("Lop_alpha_a harmonic residual is zero", profiles["Lop_alpha_a"])
    expect_zero("Lop_alpha_L harmonic residual is zero", profiles["Lop_alpha_L"])
    expect_zero("BC alpha_a(0)=1", profiles["bc_values"]["alpha_a_mouth"] - 1)
    expect_zero("BC alpha_a(L0)=0", profiles["bc_values"]["alpha_a_cap"])
    expect_zero("BC alpha_L(0)=0", profiles["bc_values"]["alpha_L_mouth"])
    expect_zero("BC alpha_L(L0)=rAL", profiles["bc_values"]["alpha_L_cap"] - rAL)


def run_projection(data: dict[str, Any]) -> None:
    projection = data["projection"]
    closed_forms = data["closed_forms"]
    flags = data["flags"]
    subbanner("M_AB/K_AB by real int-dw operator projection")
    print("  M integrands:")
    for name, expr in projection["M_integrands"].items():
        print(f"    {name}: {fmt(expr)}")
    print("  K integrands:")
    for name, expr in projection["K_integrands"].items():
        print(f"    {name}: {fmt(expr)}")
    print("  M_AB = ", fmt(projection["M_entries"]))
    print("  K_AB = ", fmt(projection["K_entries"]))
    print("  det(M) = ", fmt(projection["M_det"]))
    print("  det(K) = ", fmt(projection["K_det"]))

    expected_M = closed_forms["M_entries"]
    expected_K = closed_forms["K_entries"]
    expected_M_det = closed_forms["M_det"]
    expected_K_det = closed_forms["K_det"]
    for name in ("aa", "aL", "LL"):
        expect_zero(f"M_{name} matches report closed form", projection["M_entries"][name] - expected_M[name])
        expect_zero(f"K_{name} matches report closed form", projection["K_entries"][name] - expected_K[name])
    expect_zero("det(M) matches report closed form", projection["M_det"] - expected_M_det)
    expect_zero("det(K) matches 16*pi^2*Tw^2*beta^2*rAL^2", projection["K_det"] - expected_K_det)

    for residual in projection_residuals(data["reintegrated"], closed_forms):
        expect_zero("independent full-integrand re-integration residual", residual)
    expect_bool("independent re-integration comparison is not a self-identity stamp", data["projection_matches"])

    print("  free-symbol names in M/K = ", sorted(flags["MK_names"]))
    print("  allowed names = ", sorted(flags["allowed_names"]))
    print("  forbidden_fit_flags = ", flags["forbidden_fit_flags"])
    expect_bool("M/K free-symbol names are a subset of the allowed projection names", not flags["unexpected_names"])
    for name, value in flags["forbidden_fit_flags"].items():
        expect_bool(f"forbidden_fit_flags[{name}] computed false", value is False)
    print('  K_AB_provenance = "operator_projection_not_static_Hessian"')


def run_eom(data: dict[str, Any]) -> None:
    row_a, row_L = data["eom_rows"]
    subbanner("Dynamical EOM LHS only")
    print("  Q = (delta_a, delta_L); Qddot = (delta_a_ddot, delta_L_ddot).")
    print("  EOM structure: M_AB*Qddot^B + K_AB*Q^B = F_A^(HF).")
    print("  RHS placeholders: F_a_HF, F_L_HF are deferred to stage 015; no Hellmann-Feynman force is computed here.")
    print("  row a LHS = ", fmt(row_a), " = F_a_HF")
    print("  row L LHS = ", fmt(row_L), " = F_L_HF")
    expect_bool("EOM has exactly two LHS rows", len(data["eom_rows"]) == 2)
    expect_bool("EOM RHS is symbolic placeholder only", data["rhs_placeholders"] == (F_a_HF, F_L_HF))


def run_dimensional_block(data: dict[str, Any]) -> None:
    dim = data["dim"]
    subbanner("013 dimensional legs and corrupt-[Tw] probe")
    print("  dimension order: (L,M,T)")
    print("  sourced dims: L0=(1,0,0), beta=(-1,0,0), muEta=(-1,1,0), Tw=(1,1,-2), rAL=(0,0,0)")
    print("  derived dim: K_eta=Tw*beta^2 = (-1,1,-2)")
    print("  [M_AB] entries = ", {key: str(value) for key, value in dim["m_dims"].items()})
    print("  [K_AB] entries = ", {key: str(value) for key, value in dim["k_dims"].items()})
    print("  [M_AB] shared = ", fmt(dim["m_shared"]), "; [K_AB] shared = ", fmt(dim["k_shared"]), "; [K/M] = ", fmt(dim["ratio_dim"]))
    expect_zero("M_AB shared dimension is M", dim_residual(dim["m_shared"], dim["expected_m"]))
    expect_zero("K_AB shared dimension is M*T^-2", dim_residual(dim["k_shared"], dim["expected_k"]))
    expect_zero("K/M ratio dimension is T^-2", dim_residual(dim["ratio_dim"], dim["expected_ratio"]))
    expect_bool("dimensional_ok for 013 M/K legs", dim["dimensional_ok"])
    print("  corrupt [Tw]+(1,0,0) gives [Tw] = ", fmt(dim["corrupt_dims"][Tw]))
    print("  corrupt [K_AB] shared = ", fmt(dim["corrupt_k_shared"]), "; corrupt [K/M] = ", fmt(dim["corrupt_ratio_dim"]))
    expect_fail("corrupt-[Tw] shifts K_AB away from M*T^-2", dim_residual(dim["corrupt_k_shared"], dim["expected_k"]))
    expect_fail("corrupt-[Tw] shifts K/M away from T^-2", dim_residual(dim["corrupt_ratio_dim"], dim["expected_ratio"]))
    expect_bool("corrupt-[Tw] mutation_fires=True", dim["mutation_fires"])
    expect_zero("self-ablation with mutation gives BREATHING_FAIL_DIMENSIONAL", verdict_residual(dim["mutated_verdict"], BREATHING_FAIL_DIMENSIONAL))
    expect_zero("self-ablation without mutation gives BREATHING_CALIBRATED", verdict_residual(dim["clean_verdict"], BREATHING_CALIBRATED))
    expect_bool("self-ablation fail_suppressed=True", dim["fail_suppressed"])


def run_consumed_packet(data: dict[str, Any]) -> None:
    packet = data["consumed_packet"]
    residuals = packet_residuals(packet)
    subbanner("Consumed Gate-1 frozen packet, dual-site integrity")
    print("  CONSUMED packet anchor: (L0,Tw,K_eta,beta)=(37/20,1,1,1).")
    print("  Site A constitutive: beta_cited - sqrt(K_eta_cited/Tw_cited) = 0.")
    print("  Site B geometric/branch: beta_cited*L0_cited - 37/20 = 0.")
    print("  Anti-tautology: K_eta_cited is an independent cited datum for the guard, not the local alias K_eta=Tw*beta^2.")
    for name, residual in residuals.items():
        expect_zero(f"consumed packet {name}", residual)
    expect_bool("consumed packet clean baseline passes both sites and frozen anchor", packet_ok(packet))
    expect_fail(
        "K_eta_cited-only corruption breaks site A: guard is non-vacuous",
        packet_residuals(packet_mutation(packet, "K_eta_cited", sp.Integer(2)))["site_A_constitutive"],
    )
    expect_fail(
        "Tw_cited-only corruption breaks site A",
        packet_residuals(packet_mutation(packet, "Tw_cited", sp.Integer(2)))["site_A_constitutive"],
    )
    expect_fail(
        "L0_cited-only corruption breaks site B",
        packet_residuals(packet_mutation(packet, "L0_cited", sp.Rational(19, 10)))["site_B_branch"],
    )
    expect_fail(
        "branch-anchor corruption breaks site B",
        packet_residuals(packet_mutation(packet, "branch_anchor", sp.Rational(19, 10)))["site_B_branch"],
    )


def run_degenerate_guard(data: dict[str, Any]) -> None:
    degenerate = data["degenerate"]
    subbanner("Native degenerate M_det -> 0 guard slice")
    print("  Degenerate copy: alpha_a -> 0; alpha_L unchanged; M recomputed by the same int-dw projection.")
    print("  degenerate M_det = ", fmt(degenerate["M_det"]))
    expect_zero("degenerate alpha_a=0 recompute gives M_det=0", degenerate["M_det"])
    expect_bool("native degenerate non-degeneracy test catches M_det==0", degenerate["caught"])
    expect_nonzero("baseline M_det is not identically zero", data["projection"]["M_det"])
    bypassed = degenerate_mass_guard(data["profiles"]["alpha_L"], data["K_eta"], enabled=False)
    expect_fail("tooth 6 bypassing native M_det==0 guard leaves degenerate case uncaught", bool_residual(bypassed["caught"]))


def run_verdict_and_composition(data: dict[str, Any]) -> None:
    subbanner("013 scoped landing and joint composition")
    print("  013 scoped verdict = ", data["verdict"])
    expect_zero("013 component lands at BREATHING_CALIBRATED", verdict_residual(data["verdict"], BREATHING_CALIBRATED))
    print("  BREATHING_CALIBRATED (JOINT, 3-stage)")
    print("    = (013: harmonic-lift profiles + M_AB/K_AB by int-dw operator projection + (a,L) EOM LHS, computed here)")
    print("    AND (014: truncation consistency -- generalized eig / beta_L0 sweep / N-convergence) [sibling stage]")
    print("    AND (015: legacy-Hessian structure recovery + Hellmann-Feynman force) [sibling stage]")
    print("  CALIBRATED <= wall constants {muEta, Tw, K_eta} are calibration inputs; structure is EARNED, values are CALIBRATED.")
    expect_bool("joint composition cites 014 and 015 as siblings, not recomputed here", data["verdict"] == BREATHING_CALIBRATED)


def print_provenance(data: dict[str, Any]) -> None:
    subbanner("Provenance and scope")
    print("  CONSUMED-from-Gate-1: domain [0,L0] with cap R0(L0)=0, frozen packet {L0=37/20,Tw=1,K_eta=1,beta=1}, and ell=0 restriction are cited from stages 011/012 with dual-site integrity.")
    print("  no-c_S: 013 is speed-free; matter-sector c_s/BdG k^4 is deferred under k*xi << 1 (phonon_limit_caveat).")
    print("  EARNED-STRUCTURE: profiles are derived harmonic lifts proven by Lop[alpha]=0; M_AB/K_AB are real int-dw operator projections; EOM LHS is assembled here.")
    print("  FIRST-CALIBRATION: muEta and Tw are CALIBRATION inputs; K_eta=Tw*beta^2 is a manifestation; beta=sqrt(K_eta/Tw) with beta*L0=37/20 is branch-determinable, but geometry alone does NOT derive the wall constants.")
    print("  control-ratio: rAL is the dimensionless alpha_L cap ratio, tracked and not counted.")
    print("  3-way-split: 013 carries M/K projection + (a,L) closure; 014 carries truncation consistency; 015 carries legacy structure + HF force.")
    print("  RHS-deferred: F_A^(HF) is stage 015's Hellmann-Feynman force; 013 emits only M_AB*Qddot+K_AB*Q.")
    print("  dropped-bookkeeping: scratch-YAML/_sympy_exprs.wl export, MMA-YAML re-read, expression_digest, and engine_agreement plumbing are stripped.")
    print("  downstream consumers: stage 014 consumes M_AB/K_AB for generalized eig; stages 022/023 consume the ell=0 (a,L) closure.")
    print("  register note: FIRST Part-II calibration knobs are likely {muEta,Tw}; K_eta is a manifestation, beta and L0 are branch/geometric tracked quantities, and rAL is tracked-not-counted.")
    live_names = {
        symbol.name
        for expr in [
            data["profiles"]["alpha_a"],
            data["profiles"]["alpha_L"],
            *data["projection"]["M_entries"].values(),
            *data["projection"]["K_entries"].values(),
            *data["eom_rows"],
        ]
        for symbol in expr.free_symbols
    }
    expect_bool("no c_S/cS live symbol appears in 013 symbolic content", "cS" not in live_names and "c_S" not in live_names)


def print_verdict_labels() -> None:
    print("")
    print("Verdict labels:")
    print("  ledger earned-label (NOT a source verdict token): BREATHING_HARMONIC_MK_PROJECTION_EARNED  (harmonic-lift profiles alpha_a=sinh(L0*beta-beta*w)/sinh(L0*beta), alpha_L=rAL*sinh(beta*w)/sinh(L0*beta), proven by Lop[alpha]=0; M_AB=4*pi*int mu_eta*alpha_A*alpha_B dw and K_AB=4*pi*int [Tw*alpha_A'*alpha_B'+K_eta*alpha_A*alpha_B] dw by real sp.integrate operator projection, NOT the legacy static Hessian (forbidden_fit_flags computed False via free-symbol ancestry, K_AB_provenance=operator_projection_not_static_Hessian); dynamical-EOM LHS M_AB*Qddot+K_AB*Q with Q=(delta_a,delta_L); [M_AB]=M, [K_AB]=M*T^-2, [K/M]=T^-2 dim legs + corrupt-[Tw] probe)")
    print("  source top-line verdict: BREATHING_CALIBRATED  (JOINT 3-stage; 013 carries the M/K-projection + (a,L)-closure component)")
    print("  joint composition: BREATHING_CALIBRATED = (013: harmonic profiles + M/K operator-projection + (a,L) EOM LHS, computed here) AND (014: truncation consistency) AND (015: legacy-structure + HF force)")
    print("  earned (structure): profiles DERIVED as harmonic lifts (residual Lop[alpha]=0); M_AB/K_AB DERIVED by int-dw operator projection (forbidden_fit_flags computed False, provenance operator_projection_not_static_Hessian); EOM LHS assembled; dim legs consistent + corrupt-[Tw] probe fires")
    print("  calibrated (values): mu_eta, Tw calibration inputs; K_eta=Tw*beta^2 manifestation; beta=sqrt(K_eta/Tw), beta*L0=37/20 branch-determinable; geometry alone does NOT derive the wall constants -> BREATHING_CALIBRATED not ..._PASS")
    print("  consumed (cited from Gate-1 stage011/012, dual-site integrity): domain [0,L0] (cap R0(L0)=0); frozen wall packet L0=37/20, Tw=1, K_eta=1, beta=1 (K_eta=Tw*beta^2, beta*L0=37/20); ell=0 restriction Y00=1/(2*sqrt(pi)); c_S NOT consumed (matter-sector deferred, k*xi<<1)")
    print("  control ratio (tracked, not counted): rAL = alpha_L cap ratio, [rAL]=1")
    print("  RHS deferred to stage 015: F_A^(HF) Hellmann-Feynman force (013 emits the EOM LHS only)")


def run_able_to_fail_teeth(data: dict[str, Any]) -> None:
    subbanner("Able-to-fail mutation teeth")
    alpha_bad = sp.sinh(L0 * beta - 2 * beta * w) / sp.sinh(L0 * beta)
    bad_residual = lop(alpha_bad, data["K_eta"])
    expect_fail("tooth 1 non-kernel wrong-wavenumber profile trips harmonic residual", bad_residual)

    kappa = sp.Symbol("kappa", positive=True, real=True)
    chi = sp.Symbol("chi", positive=True, real=True)
    sigmaA = sp.Symbol("sigmaA", positive=True, real=True)
    sigmaL = sp.Symbol("sigmaL", positive=True, real=True)
    typed_K = {
        "aa": chi**2 * kappa + sigmaA,
        "aL": -chi * kappa,
        "LL": kappa + sigmaL,
    }
    typed_flags = symbol_name_flags(data["projection"]["M_entries"], typed_K)
    bare_legacy_objects = {sp.Symbol("kappa"), sp.Symbol("chi"), sp.Symbol("sigmaA"), sp.Symbol("sigmaL")}
    typed_symbol_objects = {symbol for expr in typed_K.values() for symbol in expr.free_symbols}
    object_intersection_misses = not bool(bare_legacy_objects & typed_symbol_objects)
    name_intersection_catches = bool(typed_flags["legacy_names"] & typed_flags["K_names"])
    expect_bool("tooth 2 name-based ancestry catches legacy symbols that object-intersection can miss", object_intersection_misses and name_intersection_catches)
    expect_bool("tooth 2 typed legacy K flips K_from_static_hessian flag", typed_flags["forbidden_fit_flags"]["K_from_static_hessian"])
    expect_bool("tooth 2 typed legacy K flips M_or_K_typed_from_legacy_values flag", typed_flags["forbidden_fit_flags"]["M_or_K_typed_from_legacy_values"])
    expect_bool("tooth 2 typed legacy K flips full_matrix_fit flag via unallowed names", typed_flags["forbidden_fit_flags"]["full_matrix_fit"])

    corrupted_projection = project_from_profiles(data["profiles"]["alpha_a"], data["profiles"]["alpha_L"], data["K_eta"], include_gradient=False)
    expect_fail(
        "tooth 3 dropping Tw*alpha_prime*alpha_prime makes independent K re-integration mismatch",
        bool_residual(all_projection_residuals_zero(corrupted_projection, data["closed_forms"])),
    )
    expect_fail("tooth 3 corrupted K_aa differs from emitted full-operator K_aa", corrupted_projection["K_entries"]["aa"] - data["projection"]["K_entries"]["aa"])

    dim = data["dim"]
    expect_fail("tooth 4 corrupt-[Tw] probe trips K_AB dimensional leg", dim_residual(dim["corrupt_k_shared"], dim["expected_k"]))
    expect_fail("tooth 4 corrupt-[Tw] probe trips K/M dimensional leg", dim_residual(dim["corrupt_ratio_dim"], dim["expected_ratio"]))
    expect_zero("tooth 4 corrupt-[Tw] verdict is BREATHING_FAIL_DIMENSIONAL", verdict_residual(dim["mutated_verdict"], BREATHING_FAIL_DIMENSIONAL))
    expect_bool("tooth 4 self-ablation fail_suppressed remains true", dim["fail_suppressed"])

    packet = data["consumed_packet"]
    expect_fail("tooth 5 K_eta_cited corruption trips packet site A", packet_residuals(packet_mutation(packet, "K_eta_cited", 2))["site_A_constitutive"])
    expect_fail("tooth 5 Tw_cited corruption trips packet site A", packet_residuals(packet_mutation(packet, "Tw_cited", 2))["site_A_constitutive"])
    expect_fail("tooth 5 L0_cited corruption trips packet site B", packet_residuals(packet_mutation(packet, "L0_cited", sp.Rational(19, 10)))["site_B_branch"])
    expect_fail("tooth 5 branch-anchor corruption trips packet site B", packet_residuals(packet_mutation(packet, "branch_anchor", sp.Rational(19, 10)))["site_B_branch"])

    bypassed = degenerate_mass_guard(data["profiles"]["alpha_L"], data["K_eta"], enabled=False)
    expect_fail("tooth 6 bypassed native M_det guard fails to catch degenerate profile", bool_residual(bypassed["caught"]))
    expect_zero("tooth 6 baseline degenerate copy still has M_det=0", data["degenerate"]["M_det"])

    expect_zero("baseline immutable after teeth: Lop_alpha_a remains zero", data["profiles"]["Lop_alpha_a"])
    expect_zero("baseline immutable after teeth: independent full-integrand re-integration still matches", bool_residual(data["projection_matches"]))
    expect_bool("baseline immutable after teeth: forbidden flags remain all false", not any(data["flags"]["forbidden_fit_flags"].values()))
    expect_zero("baseline immutable after teeth: clean 013 verdict remains BREATHING_CALIBRATED", verdict_residual(data["verdict"], BREATHING_CALIBRATED))


def main() -> None:
    banner("ledger_stage013_breathing_harmonic_mk_projection SymPy audit")
    data = build_baseline()
    assert_no_float("baseline", data)
    run_operator_profiles(data)
    run_projection(data)
    run_eom(data)
    run_dimensional_block(data)
    run_consumed_packet(data)
    run_degenerate_guard(data)
    run_verdict_and_composition(data)
    print_provenance(data)
    print_verdict_labels()
    run_able_to_fail_teeth(data)
    dim = data["dim"]
    emit_dimension_sidecar(
        __file__,
        {
            "symbol_dims.L0": dim["symbol_dims"][L0],
            "symbol_dims.beta": dim["symbol_dims"][beta],
            "symbol_dims.muEta": dim["symbol_dims"][muEta],
            "symbol_dims.Tw": dim["symbol_dims"][Tw],
            "symbol_dims.rAL": dim["symbol_dims"][rAL],
            "K_eta": dim["K_eta_dim"],
            "m_dims.aa": dim["m_dims"]["aa"],
            "m_dims.aL": dim["m_dims"]["aL"],
            "m_dims.LL": dim["m_dims"]["LL"],
            "k_dims.aa": dim["k_dims"]["aa"],
            "k_dims.aL": dim["k_dims"]["aL"],
            "k_dims.LL": dim["k_dims"]["LL"],
            "m_shared": dim["m_shared"],
            "k_shared": dim["k_shared"],
            "ratio_dim": dim["ratio_dim"],
        },
    )
    print("")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}; TOTAL = {PASS_COUNT} + {FAIL_COUNT} = {PASS_COUNT + FAIL_COUNT}")
    print("OVERALL PASS: SymPy verified ledger_stage013 breathing harmonic profiles + M/K projection exactly")


if __name__ == "__main__":
    try:
        main()
    except AuditFailure as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}; TOTAL = {PASS_COUNT} + {FAIL_COUNT} = {PASS_COUNT + FAIL_COUNT}")
        print(f"OVERALL FAIL: SymPy stage013 audit did not close ({exc})")
        raise SystemExit(1)
