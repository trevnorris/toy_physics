#!/usr/bin/env python3
"""Ledger stage022 SymPy audit: cross-ell DtN fingerprints.

Standalone, print-only, exact, and zero-file-I/O.  This is the pathA_34
II-G5a EARNED-FIRST slice only: the dimensionless z-space outgoing
spherical-Hankel fingerprints for ell=0,1,2 and the stage018 quadrupole
non-regression.  The sibling stage023 owns residuals and the nullspace
departure that delivers the joint FAIL.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0

CROSS_L_FINGERPRINT_OK = "CROSS_L_FINGERPRINT_OK"
FAIL_FINGERPRINT = "FAIL_FINGERPRINT"
FAIL_QUAD_REGRESSION = "FAIL_QUAD_REGRESSION"
JOINT_PARTIAL = "FAIL_UNDERDETERMINED_NOT_PREDICTIVE (1/2)"


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


def compact(expr: Any) -> Any:
    if not isinstance(expr, sp.Basic):
        return expr
    if isinstance(expr, (sp.Equality, sp.Order)):
        return expr
    return sp.trigsimp(sp.factor(sp.cancel(sp.together(sp.simplify(expr)))))


def fmt(expr: Any) -> str:
    if isinstance(expr, bool):
        return "True" if expr else "False"
    if isinstance(expr, str):
        return expr
    try:
        return sp.sstr(compact(expr))
    except Exception:
        return sp.sstr(expr)


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
    clean = compact(residual)
    assert_no_float(name, clean)
    if clean == 0:
        _record_pass(f"PASS  {name}")
        return
    _record_fail(f"FAIL  {name}: residual = {fmt(clean)}")
    raise AuditFailure(f"{name} residual was not zero")


def expect_bool(name: str, condition: bool) -> None:
    expect_zero(name, sp.Integer(0) if bool(condition) else sp.Integer(1))


def expect_fail(name: str, residual: sp.Expr | int) -> None:
    assert_no_float(name, residual)
    clean = compact(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_pass(f"PASS  {name} produced required FAIL (residual = {fmt(clean)})")
        return
    _record_fail(f"FAIL  {name}: required mutation/ablation did not fire")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def bool_zero(residual: sp.Expr | int) -> bool:
    assert_no_float("bool_zero", residual)
    return compact(residual) == 0


def verdict_residual(actual: str, expected: str) -> sp.Integer:
    return sp.Integer(0) if actual == expected else sp.Integer(1)


z = sp.Symbol("z", real=True)
omega = sp.Symbol("omega", real=True)
a = sp.Symbol("a", positive=True, real=True)
c_s = sp.Symbol("c_s", positive=True, real=True)

ELLS = (0, 1, 2)
EXPECTED_RADIATIVE = {
    0: sp.Integer(1),
    1: sp.Rational(1, 2),
    2: sp.Rational(1, 27),
}
EXPECTED_ORDER = {ell: 2 * ell + 1 for ell in ELLS}
STAGE018_TYPED = {
    "u2_z": sp.Rational(1, 9),
    "u4_z": sp.Rational(4, 81),
    "v5_z": sp.Rational(1, 27),
}


def series_no_o(expr: sp.Expr, var: sp.Symbol, order: int) -> sp.Expr:
    return sp.expand(sp.series(expr, var, 0, order).removeO())


def spherical_j(lval: int) -> sp.Expr:
    if lval == 0:
        return sp.sin(z) / z
    if lval == 1:
        return sp.sin(z) / z**2 - sp.cos(z) / z
    if lval == 2:
        return (sp.Integer(3) / z**3 - sp.Integer(1) / z) * sp.sin(z) - sp.Integer(3) * sp.cos(z) / z**2
    raise ValueError(lval)


def spherical_y(lval: int) -> sp.Expr:
    if lval == 0:
        return -sp.cos(z) / z
    if lval == 1:
        return -sp.cos(z) / z**2 - sp.sin(z) / z
    if lval == 2:
        return (sp.Integer(1) / z - sp.Integer(3) / z**3) * sp.cos(z) - sp.Integer(3) * sp.sin(z) / z**2
    raise ValueError(lval)


def hankel_wave(lval: int, kind: str) -> sp.Expr:
    j_l = spherical_j(lval)
    y_l = spherical_y(lval)
    if kind == "outgoing_hankel1":
        return j_l + sp.I * y_l
    if kind == "incoming_hankel2":
        return j_l - sp.I * y_l
    raise ValueError(kind)


def derivation_copy_factor(lval: int, mutation: str) -> sp.Expr:
    if mutation == "none":
        return sp.Integer(1)
    if mutation == "radiative_slot_copy":
        return 1 + sp.I * z ** (2 * lval + 1)
    if mutation == "lower_imaginary_copy":
        return 1 + sp.I * z
    if lval != 2:
        raise ValueError(f"{mutation} is an ell=2-only derivation copy")
    if mutation == "u2_only_derivation_copy":
        return 1 + z**2 + z**4 / sp.Integer(18)
    if mutation == "u4_only_derivation_copy":
        return 1 + z**4
    if mutation == "v5_only_derivation_copy":
        return 1 + sp.I * z**5
    raise ValueError(mutation)


@lru_cache(maxsize=None)
def dtn_branch(lval: int, kind: str, mutation: str = "none") -> dict[str, Any]:
    h = compact(hankel_wave(lval, kind) * derivation_copy_factor(lval, mutation))
    lam = compact(z * sp.diff(h, z) / h)
    yhat = compact(-sp.Integer(lval + 1) / lam)
    series_order = max(8, 2 * lval + 4)
    lam_series = series_no_o(lam, z, series_order)
    y_series = series_no_o(yhat, z, series_order)
    radiative_power = 2 * lval + 1
    radiative_coeff = compact(y_series.coeff(z, radiative_power) / sp.I)
    imag_coefficients = {
        power: compact(sp.im(y_series.coeff(z, power)))
        for power in range(series_order)
    }
    nonzero_imaginary_powers = [
        power for power, coefficient in imag_coefficients.items() if not bool_zero(coefficient)
    ]
    first_nonzero_imag_order = nonzero_imaginary_powers[0] if nonzero_imaginary_powers else None
    lower_imaginary_zero = all(
        bool_zero(imag_coefficients[power]) for power in range(radiative_power)
    )
    return {
        "ell": lval,
        "kind": kind,
        "mutation": mutation,
        "h": h,
        "lambda": lam,
        "Y": yhat,
        "lambda_series": lam_series,
        "Y_series": y_series,
        "lambda_static": compact(lam_series.coeff(z, 0)),
        "static": compact(y_series.coeff(z, 0)),
        "radiative_power": radiative_power,
        "radiative_coeff_z": radiative_coeff,
        "imag_coefficients": imag_coefficients,
        "first_nonzero_imag_order": first_nonzero_imag_order,
        "all_lower_imag_zero": lower_imaginary_zero,
        "u2_z": compact(y_series.coeff(z, 2)),
        "u4_z": compact(y_series.coeff(z, 4)),
        "v5_z": compact(y_series.coeff(z, 5) / sp.I),
    }


@lru_cache(maxsize=None)
def pole_order_mutant_static(lval: int) -> sp.Expr:
    h_mut = compact(z * hankel_wave(lval, "outgoing_hankel1"))
    lam_mut = compact(z * sp.diff(h_mut, z) / h_mut)
    return compact(series_no_o(lam_mut, z, 3).coeff(z, 0))


@lru_cache(maxsize=None)
def build_fingerprints(radiative_mutation_ell: int | None = None) -> dict[str, Any]:
    outgoing = {
        ell: dtn_branch(
            ell,
            "outgoing_hankel1",
            "radiative_slot_copy" if ell == radiative_mutation_ell else "none",
        )
        for ell in ELLS
    }
    incoming = {ell: dtn_branch(ell, "incoming_hankel2") for ell in ELLS}
    matches: dict[str, dict[int, bool]] = {
        "radiative_coeff": {},
        "lambda_static": {},
        "scanned_order": {},
        "all_lower_imag_zero": {},
        "incoming_flips_sign": {},
        "incoming_lower_real_unchanged": {},
    }
    for ell in ELLS:
        out_l = outgoing[ell]
        in_l = incoming[ell]
        power = EXPECTED_ORDER[ell]
        matches["radiative_coeff"][ell] = bool_zero(
            out_l["radiative_coeff_z"] - EXPECTED_RADIATIVE[ell]
        )
        matches["lambda_static"][ell] = bool_zero(out_l["lambda_static"] + ell + 1)
        matches["scanned_order"][ell] = out_l["first_nonzero_imag_order"] == power
        matches["all_lower_imag_zero"][ell] = out_l["all_lower_imag_zero"]
        matches["incoming_flips_sign"][ell] = bool_zero(
            in_l["radiative_coeff_z"] + out_l["radiative_coeff_z"]
        )
        matches["incoming_lower_real_unchanged"][ell] = all(
            bool_zero(
                sp.re(in_l["Y_series"].coeff(z, k))
                - sp.re(out_l["Y_series"].coeff(z, k))
            )
            for k in range(power)
        )
    return {
        "outgoing": outgoing,
        "incoming": incoming,
        "matches": matches,
        "ok": all(all(family.values()) for family in matches.values()),
        "static_diagnostic": {ell: outgoing[ell]["static"] for ell in ELLS},
        "chi_Q_diagnostic": compact(outgoing[2]["v5_z"] / STAGE018_TYPED["v5_z"]),
    }


@lru_cache(maxsize=None)
def build_gate4_non_regression(
    break_gate4: bool = False,
    slot_mutation: str = "none",
) -> dict[str, Any]:
    kind = "incoming_hankel2" if break_gate4 else "outgoing_hankel1"
    out2 = dtn_branch(2, kind, slot_mutation)
    matches = {
        name: bool_zero(out2[name] - target)
        for name, target in STAGE018_TYPED.items()
    }
    return {
        "branch_used": out2["kind"],
        "mutation": slot_mutation,
        "derived": {name: out2[name] for name in STAGE018_TYPED},
        "typed_stage018": STAGE018_TYPED,
        "matches": matches,
        "chi_Q_diagnostic": compact(out2["v5_z"] / STAGE018_TYPED["v5_z"]),
        "ok": all(matches.values()),
    }


class TrackingGateMap(dict[str, bool]):
    def __init__(self, values: dict[str, bool]):
        super().__init__(values)
        self.reads: set[str] = set()

    def __getitem__(self, key: str) -> bool:
        self.reads.add(key)
        return super().__getitem__(key)


def local_audit_verdict(fingerprints: dict[str, Any], gate4: dict[str, Any]) -> dict[str, Any]:
    gates = TrackingGateMap(
        {
            "cross_l_fingerprints": bool(fingerprints["ok"]),
            "ell2_non_regression": bool(gate4["ok"]),
        }
    )
    if not gates["cross_l_fingerprints"]:
        verdict = FAIL_FINGERPRINT
    elif not gates["ell2_non_regression"]:
        verdict = FAIL_QUAD_REGRESSION
    else:
        verdict = CROSS_L_FINGERPRINT_OK
    return {
        "verdict": verdict,
        "verdict_read_set": frozenset(gates.reads),
        "gate_values": dict(gates),
    }


def run_local_gate(
    *,
    break_gate4: bool = False,
    fingerprint_mutation_ell: int | None = None,
) -> dict[str, Any]:
    fingerprints = build_fingerprints(fingerprint_mutation_ell)
    gate4 = build_gate4_non_regression(break_gate4)
    local = local_audit_verdict(fingerprints, gate4)
    return {"fingerprints": fingerprints, "gate4": gate4, "local": local}


def dynamic_gate4_ablation(break_gate4: bool) -> dict[str, Any]:
    with_context = run_local_gate(break_gate4=break_gate4)
    without_context = run_local_gate(break_gate4=False)
    with_mutation = with_context["local"]["verdict"]
    without_mutation = without_context["local"]["verdict"]
    verdict_trace = (
        (with_context["gate4"]["branch_used"], with_mutation),
        (without_context["gate4"]["branch_used"], without_mutation),
    )
    return {
        "rerun_gate_logic": verdict_trace[0][1] != verdict_trace[1][1],
        "verdict_trace": verdict_trace,
        "with_mutation": with_mutation,
        "without_mutation": without_mutation,
        "expected_fail": FAIL_QUAD_REGRESSION,
        "fail_suppressed": (
            with_mutation == FAIL_QUAD_REGRESSION
            and without_mutation == CROSS_L_FINGERPRINT_OK
            and with_mutation != without_mutation
        ),
    }


def run_cross_l_derivation(data: dict[str, Any]) -> None:
    fingerprints = data["fingerprints"]
    subbanner("Cross-ell outgoing DtN fingerprints (computed in z-space)")
    print("  SymPy route: explicit spherical j_ell+i*y_ell from exact sin/cos expressions.")
    print("  Lambda_ell=z*h_ell'/h_ell; Yhat_ell=-(ell+1)/Lambda_ell; coefficients come from computed series.")
    for ell in ELLS:
        out_l = fingerprints["outgoing"][ell]
        in_l = fingerprints["incoming"][ell]
        power = EXPECTED_ORDER[ell]
        print(
            f"  ell={ell}: radiative_coeff = {fmt(out_l['radiative_coeff_z'])} "
            f"at omega^{power} (z^{power}); Lambda_static = {fmt(out_l['lambda_static'])}; "
            f"static = {fmt(out_l['static'])} DIAGNOSTIC; "
            f"first_nonzero_imag_order = {out_l['first_nonzero_imag_order']} (SCANNED); "
            "incoming-flips-sign = "
            f"{fmt(fingerprints['matches']['incoming_flips_sign'][ell])}"
        )
        expect_zero(
            f"ell={ell} derived radiative_coeff matches typed cross-ell target",
            out_l["radiative_coeff_z"] - EXPECTED_RADIATIVE[ell],
        )
        expect_zero(
            f"ell={ell} Lambda_static derived from Hankel log-derivative is -(ell+1)",
            out_l["lambda_static"] + ell + 1,
        )
        expect_bool(
            f"ell={ell} scanned first nonzero imaginary power is 2*ell+1",
            out_l["first_nonzero_imag_order"] == power,
        )
        expect_bool(
            f"ell={ell} every lower imaginary coefficient vanishes",
            out_l["all_lower_imag_zero"],
        )
        expect_zero(
            f"ell={ell} incoming radiative coefficient flips only the sign",
            in_l["radiative_coeff_z"] + out_l["radiative_coeff_z"],
        )
        expect_bool(
            f"ell={ell} incoming lower real coefficients equal outgoing",
            fingerprints["matches"]["incoming_lower_real_unchanged"][ell],
        )
    print(
        "  static=1 is a de-counted diagnostic derived from -(ell+1)/Lambda_static; "
        "it is not a verdict tooth."
    )
    print(
        f"  chi_Q = 27*v5_z = {fmt(fingerprints['chi_Q_diagnostic'])} DIAGNOSTIC "
        "(de-counted because it is algebraically subsumed by v5_z=1/27)."
    )


def run_gate4_non_regression(data: dict[str, Any]) -> None:
    gate4 = data["gate4"]
    subbanner("Gate-4 ell=2 non-regression against stage018 typed fingerprint")
    print("  CHECKABLE consumption: stage018 independently earned the typed z-space tuple.")
    print(
        "  ell=2 derived tuple: "
        f"u2_z={fmt(gate4['derived']['u2_z'])}, "
        f"u4_z={fmt(gate4['derived']['u4_z'])}, "
        f"v5_z={fmt(gate4['derived']['v5_z'])}"
    )
    for name, target in STAGE018_TYPED.items():
        expect_zero(
            f"ell=2 derived {name} matches stage018 earned typed {fmt(target)}",
            gate4["derived"][name] - target,
        )
    print(
        f"  chi_Q = {fmt(gate4['chi_Q_diagnostic'])} DIAGNOSTIC only; "
        "the v5_z equality owns this content."
    )


def run_per_tooth_ablations(data: dict[str, Any]) -> None:
    fingerprints = data["fingerprints"]
    subbanner("Per-tooth derivation-copy ablations")
    for ell in ELLS:
        mutated = dtn_branch(ell, "outgoing_hankel1", "radiative_slot_copy")
        expect_fail(
            f"ell={ell} radiative derivation copy fires its own coefficient assert",
            mutated["radiative_coeff_z"] - EXPECTED_RADIATIVE[ell],
        )
        expect_fail(
            f"ell={ell} pole-order h_mut=z*h fires the derived Lambda_static assert",
            pole_order_mutant_static(ell) + ell + 1,
        )
        expect_zero(
            f"ell={ell} H1-to-H2 is inert for Lambda_static (not its mutant)",
            fingerprints["incoming"][ell]["lambda_static"]
            - fingerprints["outgoing"][ell]["lambda_static"],
        )
        claimed_outgoing = dtn_branch(ell, "incoming_hankel2")
        expect_fail(
            f"ell={ell} outgoing-to-incoming mutation fires sign-flip tooth",
            claimed_outgoing["radiative_coeff_z"]
            + fingerprints["incoming"][ell]["radiative_coeff_z"],
        )

    lower_order_mutant = dtn_branch(2, "outgoing_hankel1", "lower_imaginary_copy")
    print(
        "  +I*z lower-order Hankel-copy mutation: "
        f"first_nonzero_imag_order={lower_order_mutant['first_nonzero_imag_order']}"
    )
    expect_fail(
        "ell=2 +I*z corruption fires the SCANNED radiative-order assert",
        sp.Integer(lower_order_mutant["first_nonzero_imag_order"] - EXPECTED_ORDER[2]),
    )
    expect_bool(
        "ell=2 +I*z corruption makes a lower imaginary coefficient nonzero",
        not lower_order_mutant["all_lower_imag_zero"],
    )

    slot_mutations = {
        "u2_z": "u2_only_derivation_copy",
        "u4_z": "u4_only_derivation_copy",
        "v5_z": "v5_only_derivation_copy",
    }
    for slot, mutation in slot_mutations.items():
        mutant = build_gate4_non_regression(False, mutation)
        print(
            f"  {slot}-only copy tuple = "
            f"{{u2_z:{fmt(mutant['derived']['u2_z'])}, "
            f"u4_z:{fmt(mutant['derived']['u4_z'])}, "
            f"v5_z:{fmt(mutant['derived']['v5_z'])}}}"
        )
        expect_fail(
            f"ell=2 {slot}-only derivation copy fires its own non-regression assert",
            mutant["derived"][slot] - STAGE018_TYPED[slot],
        )
        for other_slot in STAGE018_TYPED:
            if other_slot != slot:
                expect_zero(
                    f"ell=2 {slot}-only copy leaves {other_slot} unchanged",
                    mutant["derived"][other_slot] - STAGE018_TYPED[other_slot],
                )


def run_local_verdict_and_3e(data: dict[str, Any]) -> None:
    subbanner("022-local verdict, read-set, and 3e dynamic self-ablation")
    baseline = data["local"]
    fingerprint_mutant = run_local_gate(fingerprint_mutation_ell=1)
    broken = run_local_gate(break_gate4=True)
    ablation = dynamic_gate4_ablation(True)
    neutered = dynamic_gate4_ablation(False)
    expected_reads = {"cross_l_fingerprints", "ell2_non_regression"}
    forbidden_reads = {
        "nullspace",
        "return_admittance",
        "base_verdict",
        "residuals",
        "selector",
    }
    print(f"  baseline gate values = {baseline['gate_values']}")
    print(f"  computed LOCAL_AUDIT_VERDICT read-set = {sorted(baseline['verdict_read_set'])}")
    expect_bool(
        "LOCAL_AUDIT_VERDICT reads exactly fingerprints plus ell=2 non-regression",
        baseline["verdict_read_set"] == expected_reads,
    )
    print(
        "  LOCAL_AUDIT_VERDICT read-set excludes nullspace/return/base_verdict objects = "
        f"{baseline['verdict_read_set'].isdisjoint(forbidden_reads)} DIAGNOSTIC "
        "(de-counted; the exact read-set equality tooth owns this content)."
    )
    expect_zero(
        "baseline 022-local verdict is CROSS_L_FINGERPRINT_OK",
        verdict_residual(baseline["verdict"], CROSS_L_FINGERPRINT_OK),
    )
    expect_zero(
        "corrupting an ell=1 fingerprint derivation reaches FAIL_FINGERPRINT",
        verdict_residual(fingerprint_mutant["local"]["verdict"], FAIL_FINGERPRINT),
    )
    print(
        "  3e_break_gate4: incoming ell=2 gives "
        f"v5_z={fmt(broken['gate4']['derived']['v5_z'])} -> "
        f"{broken['local']['verdict']}"
    )
    expect_zero(
        "3e incoming ell=2 reaches FAIL_QUAD_REGRESSION",
        verdict_residual(broken["local"]["verdict"], FAIL_QUAD_REGRESSION),
    )
    print(f"  3e self_ablation = {ablation}")
    expect_bool("3e self-ablation dynamically reruns 022-local gate logic", ablation["rerun_gate_logic"])
    expect_zero(
        "3e dynamic rerun with mutation is FAIL_QUAD_REGRESSION",
        verdict_residual(ablation["with_mutation"], FAIL_QUAD_REGRESSION),
    )
    expect_zero(
        "3e dynamic rerun without mutation is CROSS_L_FINGERPRINT_OK",
        verdict_residual(ablation["without_mutation"], CROSS_L_FINGERPRINT_OK),
    )
    expect_bool("3e dynamic self-ablation changes the local verdict", ablation["fail_suppressed"])
    expect_bool(
        "neutering 3e is detected as not able to fail",
        not neutered["fail_suppressed"]
        and neutered["with_mutation"] == neutered["without_mutation"],
    )


def run_scope_and_provenance(data: dict[str, Any]) -> None:
    subbanner("022 scope and PROVENANCE consumption")
    expressions = [
        branch[field]
        for branch in data["fingerprints"]["outgoing"].values()
        for field in ("h", "lambda", "Y", "lambda_series", "Y_series")
    ]
    live_names = {symbol.name for expr in expressions for symbol in expr.free_symbols}
    print(f"  live symbolic names in computed fingerprint expressions = {sorted(live_names)}")
    expect_bool("computed fingerprint algebra is dimensionless z-space only", live_names == {"z"})
    helper_names = set(globals())
    forbidden_helpers = {
        "base_" + "verdict",
        "build_" + "transfers",
        "build_" + "residuals",
        "build_" + "rank_audit",
        "build_" + "port_kernel",
    }
    expect_bool(
        "no 019 prefactor or 023 transfer/residual/nullspace helper is built",
        helper_names.isdisjoint(forbidden_helpers),
    )
    print("  z=a*omega/c_s is a provenance dictionary only; no units-restored dimensional leg is built in 022.")
    print("  CHECKABLE: stage018 typed {u2_z=1/9,u4_z=4/81,v5_z=1/27} is the non-regression reference.")
    print("  PROVENANCE: stage019/stage020/stage021 and the completed pathA_33 QUAD_CALIBRATED joint are cited DONE, not re-derived.")
    print("  PROVENANCE: stage008 raw ell=0/1/2 amplitudes; stage009/stage010 bulk Helmholtz mode.")
    print("  PROVENANCE: stage005 R1 supplies the c_s units symbol (distinct from frozen-wall c_S); no c_s value is consumed.")
    print("  EXCLUDED/DEFERRED TO 023: residual form/sign/order, return/nullspace/selector departure, and magnitude/nonzero prediction.")
    print("  dropped-bookkeeping: YAML/report writers, cross-engine scratch files, and all file I/O are absent.")


def print_verdict_labels(data: dict[str, Any]) -> None:
    subbanner("Verdict labels:")
    print("  ledger local exit-0 verdict: CROSS_L_FINGERPRINT_OK")
    print("  joint landing label (EARNED-first PARTIAL; FAIL not evaluated by 022): FAIL_UNDERDETERMINED_NOT_PREDICTIVE (1/2)")
    print("  earned here: computed cross-ell radiative coefficients/orders/Lambda_static plus the stage018 Gate-4 tuple non-regression.")
    print("  earned able-to-fail: per-ell coefficient, pole-order, scanned-order, sign-flip, per-slot Gate-4, and dynamic 3e teeth.")
    print("  de-counted diagnostics: static=1 and chi_Q=1; neither participates in the 022-local verdict.")
    print("  remaining earned + FAIL-delivering/deferred content: residual form/sign/order and native nullspace/magnitude work belong to 023.")
    print("")
    print(f"LOCAL_AUDIT_VERDICT: {data['local']['verdict']}")
    print("JOINT_LANDING_LABEL (PARTIAL; FAIL not evaluated by 022):")
    print(f"  {JOINT_PARTIAL}")


def main() -> int:
    banner("ledger_stage022_cross_l_fingerprints_sympy_audit")
    print("Target stem confirmed: ledger_stage022_cross_l_fingerprints")
    print("Engine: SymPy exact hand-built spherical j/y Hankel algebra; no floats/tolerances; zero file I/O.")
    data = run_local_gate()
    run_cross_l_derivation(data)
    run_gate4_non_regression(data)
    run_per_tooth_ablations(data)
    run_local_verdict_and_3e(data)
    run_scope_and_provenance(data)
    print_verdict_labels(data)
    return 0


if __name__ == "__main__":
    try:
        exit_code = main()
    except Exception as exc:
        if not isinstance(exc, AuditFailure):
            _record_fail(f"UNCAUGHT exception: {exc!r}")
        banner("Tallies")
        total = PASS_COUNT + FAIL_COUNT
        print(f"TALLY sympy: {PASS_COUNT} pass + {FAIL_COUNT} fail = {total} checks")
        print("OVERALL FAIL")
        raise SystemExit(1) from exc
    banner("Tallies")
    total = PASS_COUNT + FAIL_COUNT
    print(f"TALLY sympy: {PASS_COUNT} pass + {FAIL_COUNT} fail = {total} checks")
    if FAIL_COUNT == 0 and exit_code == 0:
        print("OVERALL PASS")
        raise SystemExit(0)
    print("OVERALL FAIL")
    raise SystemExit(1)
