#!/usr/bin/env python3
"""Ledger stage024 SymPy audit: density-native l=2 two-port derivation.

Standalone, print-only, no arguments, and no file I/O.  This is only the
stage024 derivation slice: the full 2x2 inverse is load-bearing, while the
adjugate/determinant expression is an independent factorization oracle.
Sibling stages own vector-taint, continuity-lineage, and port closure checks.
"""

from __future__ import annotations

from typing import Any

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0

LOCAL_VERDICT = "DENSITY_TWO_PORT_DERIVED"
JOINT_PARTIAL = (
    "DENSITY_PORT_HOSTED (1/4, DERIVATION — the density-native two-port "
    "N0_den; 025 = vector-freedom taint, 026 = continuity lineage, "
    "027 = port checks + closure)"
)


class AuditFailure(AssertionError):
    """A named stage024 audit assertion failed."""


def banner(title: str) -> None:
    print("")
    print("=" * len(title))
    print(title)
    print("=" * len(title))


def subbanner(title: str) -> None:
    print("")
    print(title)
    print("-" * len(title))


def compact(expr: Any) -> Any:
    if isinstance(expr, sp.MatrixBase):
        return expr.applyfunc(compact)
    if not isinstance(expr, sp.Basic):
        return expr
    return sp.factor(sp.cancel(sp.together(sp.simplify(expr))))


def fmt(expr: Any) -> str:
    if isinstance(expr, str):
        return expr
    if isinstance(expr, sp.MatrixBase):
        return sp.sstr(compact(expr))
    return sp.sstr(compact(expr))


def _record_pass(name: str) -> None:
    global PASS_COUNT
    PASS_COUNT += 1
    print(f"PASS  {name}")


def _record_fail(name: str, detail: str) -> None:
    global FAIL_COUNT
    FAIL_COUNT += 1
    print(f"FAIL  {name}: {detail}")


def expect_zero(name: str, residual: sp.Expr | int) -> None:
    clean = compact(sp.sympify(residual))
    if clean == 0:
        _record_pass(name)
        return
    _record_fail(name, f"residual = {fmt(clean)}")
    raise AuditFailure(name)


def expect_bool(name: str, condition: bool) -> None:
    expect_zero(name, 0 if bool(condition) else 1)


# Density two-port symbols. I25 is a typed l=2 continuity-moment input;
# stage026, not this derivation slice, validates its lineage.
a, c_s, rho_eff = sp.symbols("a c_s rho_eff", positive=True, real=True)
I25, Xi_Q = sp.symbols("I25 Xi_Q", nonzero=True, real=True)
eta_q, eta_phi = sp.symbols("eta_q eta_phi", nonzero=True, real=True)
varpi_q2, varpi_Phi2 = sp.symbols(
    "varpi_q2 varpi_Phi2", positive=True, real=True
)
lambda_c = sp.Symbol("lambda_c", real=True)
q2, Phi2 = sp.symbols("q2 Phi2", real=True)

DENSITY_HOST_SYMBOLS = {
    q2,
    Phi2,
    c_s,
    a,
    rho_eff,
    I25,
    Xi_Q,
    eta_q,
    eta_phi,
    varpi_q2,
    varpi_Phi2,
    lambda_c,
}
DENSITY_HOST_ORDER = (
    "q2",
    "Phi2",
    "c_s",
    "a",
    "rho_eff",
    "I25",
    "Xi_Q",
    "eta_q",
    "eta_phi",
    "varpi_q2",
    "varpi_Phi2",
    "lambda_c",
)
FORBIDDEN_VECTOR_NAMES = {
    "A_w",
    "F_muw",
    "Jw",
    "U",
    "W",
    "Omega_U",
    "Omega_W",
    "R_mix",
    "g_U",
    "g_W",
}


def make_N0(response_argument: sp.Expr) -> sp.Expr:
    """Load-bearing response-to-numerator construction."""
    return compact(response_argument**2)


def guard_nonsingular(delta: sp.Expr, context: str) -> None:
    """Named domain guard; it must run before Matrix.inv()."""
    expect_bool(
        f"B Delta!=0 nonsingular guard before inversion ({context})",
        compact(delta) != 0,
    )


def derive_density_two_port(
    g_base_value: sp.Expr,
    *,
    context: str,
    operator_lambda: sp.Expr = lambda_c,
) -> dict[str, Any]:
    """Compute the Phi2 response using the full static-operator inverse."""
    g_q = compact(g_base_value * eta_q)
    g_phi = compact(g_base_value * eta_phi)
    static_operator = sp.Matrix(
        [[varpi_q2, -operator_lambda], [-operator_lambda, varpi_Phi2]]
    )
    source = sp.Matrix([g_q, g_phi])
    delta = compact(varpi_q2 * varpi_Phi2 - operator_lambda**2)

    # Tooth B is deliberately before the load-bearing inversion.
    guard_nonsingular(delta, context)
    response = compact((static_operator.inv() * source)[1])

    # Independent adjugate/determinant oracle; it does not construct N0_den.
    p_den = compact(varpi_q2 * g_phi + operator_lambda * g_q)
    n0_den = make_N0(response)
    return {
        "static_operator": static_operator,
        "source": source,
        "delta": delta,
        "p_den": p_den,
        "response": response,
        "N0_den": n0_den,
        "g_base": compact(g_base_value),
    }


def run_audit() -> dict[str, Any]:
    g_base = compact(
        sp.sqrt(rho_eff) * c_s**2 * I25 * Xi_Q / a ** sp.Rational(7, 2)
    )
    baseline = derive_density_two_port(g_base, context="baseline")

    subbanner("A. inverse factorization and response-to-N0 dataflow")
    expect_zero(
        "A inverse response equals independent P_den/Delta oracle",
        baseline["response"] - baseline["p_den"] / baseline["delta"],
    )
    delta_probe = sp.Symbol("delta_probe", nonzero=True, real=True)
    expect_bool(
        "A make_N0 runtime dataflow probe reads its response argument",
        compact(
            make_N0(baseline["response"] + delta_probe)
            - make_N0(baseline["response"])
        )
        != 0,
    )

    subbanner("C. coupling-vanishes boundary")
    expect_bool(
        "C baseline N0_den is nonzero for symbolic nonzero coupling",
        compact(baseline["N0_den"]) != 0,
    )
    zero_control = derive_density_two_port(sp.Integer(0), context="g_base=0 control")
    expect_zero(
        "C zero-control N0_den|g_base=0 is computed zero",
        zero_control["N0_den"],
    )

    subbanner("D. density-only host set")
    live_symbols = set(baseline["N0_den"].free_symbols)
    live_names = {symbol.name for symbol in live_symbols}
    expect_bool(
        "D N0_den host-set membership is density-only and vector-free",
        live_symbols <= DENSITY_HOST_SYMBOLS
        and live_names.isdisjoint(FORBIDDEN_VECTOR_NAMES),
    )

    baseline["live_symbols"] = live_symbols
    baseline["g_base"] = g_base
    return baseline


def emit(data: dict[str, Any]) -> None:
    subbanner("Derived density two-port transcript")
    print("port_picture: ii two-port(q2,Phi2)")
    print("method: SymPy full 2x2 Matrix.inv() response")
    print(f"static_operator: {fmt(data['static_operator'])}")
    print(f"g_base: {fmt(data['g_base'])}")
    print(f"Phi2_response: {fmt(data['response'])}")
    print(f"Delta_den: {fmt(data['delta'])}")
    print(f"P_den: {fmt(data['p_den'])}")
    print(f"N0_den (canonical factored): {fmt(data['N0_den'])}")
    print("Delta!=0 guard status: PASS (checked before each inversion)")
    print(
        "density-only host-set: {"
        + ", ".join(DENSITY_HOST_ORDER)
        + "}"
    )
    print(
        "N0_den live symbols: {"
        + ", ".join(
            name
            for name in DENSITY_HOST_ORDER
            if name in {symbol.name for symbol in data["live_symbols"]}
        )
        + "}"
    )

    subbanner("physical_relations provenance")
    print(
        "varpi_q2 <- pathA_32 wall: K2/M2 = (c_s/a)^2*kappa_q "
        "from the wall angular l=2 operator (stages016/017)"
    )
    print(
        "varpi_Phi2 <- pathA_29 bulk: (c_s/a)^2*(6+(m*a)^2) "
        "from the bulk Helmholtz l=2 mode (stages009/010/016)"
    )
    print(
        "lambda_c <- projected continuity/interface: "
        "(c_s/a)^2*lambda_hat_Q"
    )
    print(
        "I25: typed l=2 continuity-moment input; lineage validation is "
        "forward-referenced to stage026"
    )
    print(
        "rho_eff: effective reduced-3D mass density; its literal reduction "
        "is SIM_DEFERRED/GAP and it is not stage005 rho0"
    )

    subbanner("Verdict labels")
    print(f"LOCAL_AUDIT_VERDICT: {LOCAL_VERDICT}")
    print(f"JOINT_LANDING_LABEL (PARTIAL): {JOINT_PARTIAL}")


def main() -> int:
    banner("ledger_stage024_density_port_derivation_sympy_audit")
    print("Target stem confirmed: ledger_stage024_density_port_derivation")
    print("Engine: exact SymPy inverse derivation; no floats; zero file I/O.")
    data = run_audit()
    emit(data)
    return 0


if __name__ == "__main__":
    try:
        exit_code = main()
    except Exception as exc:
        if not isinstance(exc, AuditFailure):
            _record_fail("UNCAUGHT exception", repr(exc))
        banner("Tallies")
        total = PASS_COUNT + FAIL_COUNT
        print(f"TALLY sympy: {PASS_COUNT} pass + {FAIL_COUNT} fail = {total} checks")
        print("OVERALL FAIL")
        raise SystemExit(1) from exc

    banner("Tallies")
    total = PASS_COUNT + FAIL_COUNT
    print(f"TALLY sympy: {PASS_COUNT} pass + {FAIL_COUNT} fail = {total} checks")
    if exit_code == 0 and FAIL_COUNT == 0:
        print("OVERALL PASS")
        raise SystemExit(0)
    print("OVERALL FAIL")
    raise SystemExit(1)
