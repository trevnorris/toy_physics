#!/usr/bin/env python3
"""S11c-b hand-coded cross-engine reconciliation for the parsed core families.

The imported T7 comparator remains the case extractor and accounting authority.
This layer adds the enumerated spelling map and the section 3a coefficient
bridge built from emitted invariant objects, then measures joined residuals.  It has no expected
zero/nonzero payload and never asserts on a measured result.

The section 3a invariants remain unfolded; only coefficient atoms acquire the
supplied source normalization. Other energy-basis representatives are compared
after spelling reconciliation.  The coupling adjointness cases use the
comparator's already compact-support-IBP-reduced context and receive no further
physics-bearing fold.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import sys
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Sequence

import sympy as sp
from sympy.core.function import AppliedUndef
from sympy.core.relational import Relational

import S11c_b_cross_engine_comparator as C
from S11b_cross_engine_comparator import ResidualAssociation


CORE_FAMILIES = (
    "COUPLING_KERNEL",
    "COUPLING_KERNEL_TERM_ORIGINS",
    "SLAB_OPERATOR",
    "SLAB_OPERATOR_TERM_ORIGINS",
    "MU_THETA_OPERATOR",
    "ADMISSIBILITY_OPERATOR_OPERAND",
    "ADMISSIBILITY_RESIDUAL",
    "ADMISSIBILITY_SUPPORT_OPERAND",
    "ENERGY_BASIS_VARIABLE",
    "ENERGY_BASIS_COUNT",
    "ENERGY_BASIS_NEW_INVARIANTS",
    "ENERGY_BASIS_OMISSIONS",
    "DIMENSIONS",
    "HOMOGENEITY_BASE_OPERAND",
    "HOMOGENEITY_CONTROL_OPERAND",
    "HOMOGENEITY_RESIDUAL",
)

NAMESPACE_INCOMPLETE_FAMILIES = (
    "CONTROL_FORM_ABLATED_OPERAND",
    "CONTROL_FORM_BASE_OPERAND",
    "CONTROL_FORM_RESIDUAL",
    "CONTROL_INDEPENDENCE_BASE_OPERAND",
    "CONTROL_INDEPENDENCE_CORRUPTED_OPERAND",
    "CONTROL_INDEPENDENCE_RESIDUAL",
    "REP_INVARIANCE_EULERIAN_OPERAND",
    "REP_INVARIANCE_MATERIAL_OPERAND",
    "REP_INVARIANCE_RESIDUAL",
    "UNIFORM_LIMIT_RESIDUAL",
    "UNIFORM_LIMIT_S11B_OPERAND",
    "UNIFORM_LIMIT_S11CB_OPERAND",
)


# HAND-VERIFIED spelling equivalences only.  Source citations are deliberately
# repeated per entry: each line is independently auditable.  ``muR`` maps to
# the uniform ``mu_R``; the distinct applied head ``modulusBackground`` maps to
# ``mu_R_bg`` (WL:204-214,294-295; PY:121-158; spec:176-181).
WL_TO_PY_RENAME = {
    "WZero": "W_0",  # WL:206-214,294; PY:121-126,150-158; spec:176-181 (uniform reference thickness).
    "muR": "mu_R",  # WL:206-214,295; PY:121-125; spec:176-181 (uniform curl modulus, not mu_R_bg).
    "rhoBr": "rho_br",  # WL:206-214,371-376; PY:121-127,179-190; spec:184-192 (integrated brane density).
    "rhoM": "rho_m",  # WL:206-214; PY:121-130; spec:135-143 (bulk mass density).
    "etaBg": "eta_bg",  # WL:294-295,323-327; PY:150-158; spec:176-203 (zero-jet contrast bookkeeper).
    "sigmaW": "sigma_W",  # WL:323-355; PY:150-177; spec:180-203 (first-background-jet bookkeeper).
    "LWidth": "L_W",  # WL:339-355; PY:150-158; spec:176-203 (independent profile length).
    "w1ProfileZero": "w1_profile",  # WL:292-295; PY:150-177; spec:176-181 (thickness-profile zero jet).
    "m1ProfileZero": "m1_profile",  # WL:292-295; PY:150-177; spec:176-181 (modulus-profile zero jet).
    "w1JetOne": "w1_profile_d1",  # WL:292-336; PY:160-177; spec:180-181 (thickness-profile first jet 1).
    "w1JetTwo": "w1_profile_d2",  # WL:292-336; PY:160-177; spec:180-181 (thickness-profile first jet 2).
    "w1JetThree": "w1_profile_d3",  # WL:292-336; PY:160-177; spec:180-181 (thickness-profile first jet 3).
    "m1JetOne": "m1_profile_d1",  # WL:292-337; PY:164-177; spec:180-181 (modulus-profile first jet 1).
    "m1JetTwo": "m1_profile_d2",  # WL:292-337; PY:164-177; spec:180-181 (modulus-profile first jet 2).
    "m1JetThree": "m1_profile_d3",  # WL:292-337; PY:164-177; spec:180-181 (modulus-profile first jet 3).
    "widthBackground": "W_bg",  # WL:283-290,366-369; PY:150-168; spec:176-192 (varying thickness field).
    "modulusBackground": "mu_R_bg",  # WL:283-295,366-369; PY:150-168; spec:176-192 (varying curl modulus).
    "eWBackground": "e_W_bg",  # WL:451-454; PY:150-158,649-651; spec:246-247 (local-background thickness fraction).
    "kappaW": "kappa_W",  # WL:471-478; PY:132-148,1138-1140; spec:188-197 (thickness-gradient coefficient).
    "gradientThetaEwCoefficient": "kappa_theta_W",  # WL:471-478; PY:141-148,1138-1140; spec:93-115 (mixed scalar-gradient coefficient).
    "gradientThetaCoefficient": "kappa_theta",  # WL:471-478; PY:141-148,1138-1140; spec:93-115 (densification-gradient coefficient).
    "cCross": "C",  # WL:471-479; PY:132-148,1133-1140; spec:93-115 (theta/e_W cross coefficient).
    "thetaDivUCoefficient": "G_theta_u",  # WL:471-479; PY:139-148,1130-1140; spec:93-115 (theta-divergence coupling).
    "ewDivUCoefficient": "G_W_u",  # WL:471-479; PY:139-148,1130-1140; spec:93-115 (thickness-divergence coupling).
    "kW": "k_W",  # WL:471-479; PY:132-148,1133-1140; spec:93-115 (thickness restoring coefficient).
    "muW": "mu_W",  # WL:838,923; PY:132-148,1490; spec:188-197 (thickness inertia).
    "modulusShear": "mu_S",  # WL:471-479; PY:132-148,1142-1145; spec:93-115 (symmetric-gradient modulus).
    "frequency": "omega",  # WL:216-219; PY:343-351 and inherited response laws; spec:137-143 (response frequency).
    "lambdaAZero": "Lambda_A_0",  # WL:216-219; PY:341-351; spec:137-143 (affinity face response).
    "lambdaVZero": "Lambda_V_0",  # WL:216-219; PY:341-351; spec:137-143 (velocity face response).
    "lambdaXZero": "Lambda_X_0",  # WL:216-219; PY:341-351; spec:137-143 (traction face response).
    "tauA": "tau_A",  # WL:216-219; PY:343-351; spec:137-143 (affinity relaxation time).
    "tauV": "tau_V",  # WL:216-219; PY:343-351; spec:137-143 (velocity relaxation time).
    "tauX": "tau_X",  # WL:216-219; PY:343-351; spec:137-143 (traction relaxation time).
    "theta_d1": "grad_theta_1",  # WL:239,246-252; PY:192-210; spec:244-250 (theta first jet 1).
    "theta_d2": "grad_theta_2",  # WL:239,246-252; PY:192-210; spec:244-250 (theta first jet 2).
    "theta_d3": "grad_theta_3",  # WL:239,246-252; PY:192-210; spec:244-250 (theta first jet 3).
    "thetaWave": "theta",  # WL:239,263; PY:192-218; spec:54-77,244-250 (Eulerian densification field).
    "eWave": "e_W",  # WL:240,263-267; PY:121-127,209-218; spec:54-77,244-250 (thickness fraction field).
    "eWField": "e_W",  # WL:206-214,240; PY:121-127; spec:54-77,244-250 (same e_W field in inherited table).
    "virtualEw": "delta_v_e_W",  # WL:237-244; PY:429-436; spec:54-73,135-143 (virtual e_W variation).
    "pressureUpper": "delta_p_plus",  # WL:266-267,702-703,767-769; PY:367-375; spec:139-143 (upper-face pressure; never differentiated in WL).
    "pressureLower": "delta_p_minus",  # WL:266-267,702-703,767-769; PY:367-375; spec:139-143 (lower-face pressure; never differentiated in WL).
    "longitudinalTrialPotential": "phi_L",  # WL:229-232; PY:254-272,1889-1912; spec:294-309 (longitudinal trial potential).
    "longitudinalTestPotential": "psi_L_s11cb",  # WL:229-232; PY:1889-1912; spec:294-309 (longitudinal test potential).
    "thetaTrial": "theta_probe",  # WL:233-236; PY:254-272; spec:294-309 (theta-sector trial probe).
    "ewTrial": "e_W_probe",  # WL:233-236; PY:254-272; spec:294-309 (thickness-sector trial probe).
    "thetaTest": "v_theta_s11cb",  # WL:233-236; PY:1889-1907; spec:294-309 (theta-sector test probe).
    "ewTest": "v_e_W_s11cb",  # WL:233-236; PY:1889-1907; spec:294-309 (thickness-sector test probe).
    "transverseTrialPotentialOne": "A_T_s11cb_1",  # WL:223-228; PY:1908-1912; spec:294-309 (transverse amplitude component 1).
    "transverseTrialPotentialTwo": "A_T_s11cb_2",  # WL:223-228; PY:1908-1912; spec:294-309 (transverse amplitude component 2).
    "transverseTrialPotentialThree": "A_T_s11cb_3",  # WL:223-228; PY:1908-1912; spec:294-309 (transverse amplitude component 3).
    "transverseTestPotentialOne": "A_T_s11cb_1",  # WL:226-228; PY:1908-1912; spec:294-309 (adjoint transverse amplitude component 1).
    "transverseTestPotentialTwo": "A_T_s11cb_2",  # WL:226-228; PY:1908-1912; spec:294-309 (adjoint transverse amplitude component 2).
    "transverseTestPotentialThree": "A_T_s11cb_3",  # WL:226-228; PY:1908-1912; spec:294-309 (adjoint transverse amplitude component 3).
    "forceHoldOne": "f_hold_u_1_0",  # WL:1361-1394 support emission; PY:274-289; spec:211-230 (held body force u1).
    "forceHoldTwo": "f_hold_u_2_0",  # WL:1361-1394 support emission; PY:274-289; spec:211-230 (held body force u2).
    "forceHoldThree": "f_hold_u_3_0",  # WL:1361-1394 support emission; PY:274-289; spec:211-230 (held body force u3).
    "forceHoldTheta": "f_hold_theta_0",  # WL:1361-1394 support emission; PY:274-289; spec:211-230 (held theta force).
    "forceHoldEw": "f_hold_e_W_0",  # WL:1361-1394 support emission; PY:274-289; spec:211-230 (held thickness force).
    "tractionHoldUpperOne": "t_hold_plus_0_1",  # WL:1361-1394; PY:278-289; spec:211-230 (upper traction 1).
    "tractionHoldUpperTwo": "t_hold_plus_0_2",  # WL:1361-1394; PY:278-289; spec:211-230 (upper traction 2).
    "tractionHoldUpperThree": "t_hold_plus_0_3",  # WL:1361-1394; PY:278-289; spec:211-230 (upper traction 3).
    "tractionHoldUpperW": "t_hold_plus_0_4",  # WL:1361-1394; PY:278-289; spec:211-230 (upper thickness traction).
    "tractionHoldLowerOne": "t_hold_minus_0_1",  # WL:1361-1394; PY:278-289; spec:211-230 (lower traction 1).
    "tractionHoldLowerTwo": "t_hold_minus_0_2",  # WL:1361-1394; PY:278-289; spec:211-230 (lower traction 2).
    "tractionHoldLowerThree": "t_hold_minus_0_3",  # WL:1361-1394; PY:278-289; spec:211-230 (lower traction 3).
    "tractionHoldLowerW": "t_hold_minus_0_4",  # WL:1361-1394; PY:278-289; spec:211-230 (lower thickness traction).
}

# These trial/test/support heads have bare PY images. Their derivatives are
# decoded first; their remaining base applications lose evaluation arguments.
# Live backgrounds are classified separately below and retain their arguments.  The cited engine definitions
# show the same independent trial/test/support object in both constructions.
TRIAL_TEST_SUPPORT_APPLIED = {
    "thetaWave",
    "eWave",
    "eWField",
    "virtualEw",
    "pressureUpper",
    "pressureLower",
    "longitudinalTrialPotential",
    "longitudinalTestPotential",
    "thetaTrial",
    "ewTrial",
    "thetaTest",
    "ewTest",
    "transverseTrialPotentialOne",
    "transverseTrialPotentialTwo",
    "transverseTrialPotentialThree",
    "transverseTestPotentialOne",
    "transverseTestPotentialTwo",
    "transverseTestPotentialThree",
    "forceHoldOne",
    "forceHoldTwo",
    "forceHoldThree",
    "forceHoldTheta",
    "forceHoldEw",
    "tractionHoldUpperOne",
    "tractionHoldUpperTwo",
    "tractionHoldUpperThree",
    "tractionHoldUpperW",
    "tractionHoldLowerOne",
    "tractionHoldLowerTwo",
    "tractionHoldLowerThree",
    "tractionHoldLowerW",
}

# WL:330-355 emits derivative-order counters; PY:353-360 uses explicit
# _dNdM symbols for the same profile jets; spec:180-181 defines those jets.
PROFILE_JET_HEADS = {
    "widthProfileJet": "w1_profile",
    "modulusProfileJet": "m1_profile",
}

# Section 3a coefficients are deliberately absent from the string table.
# Their images are built from emitted invariant objects by GammaBridge below.

# Inherited verbatim from the hand-verified S11c-a reconciliation: WL has no
# x-jets/Derivative of either branch operand, so its coordinate arguments are
# inert.  Current heads are different: w is physical and is always preserved;
# only the implicit x/time evaluation arguments are removed.
INERT_APPLIED = {"mu_theta_L", "mu_theta_M"}
LIVE_BACKGROUND_APPLIED = frozenset(C.LIVE_BACKGROUND_HEADS)
# Keep the complete decode domain separate from the removable applications.
BARE_APPLIED = frozenset(TRIAL_TEST_SUPPORT_APPLIED | INERT_APPLIED)
JET_APPLIED = BARE_APPLIED | LIVE_BACKGROUND_APPLIED
APPLIED_HEAD_CLASSES = {
    **{name: "INERT" for name in INERT_APPLIED},
    **{name: "LIVE_BACKGROUND" for name in LIVE_BACKGROUND_APPLIED},
    **{name: "TRIAL_TEST_SUPPORT" for name in TRIAL_TEST_SUPPORT_APPLIED},
}
CURRENT_KEEP_W = "delta_j_bulk"

W = C.BOUND_BINDER
COMBINE_BOUND_INTEGRALS = C.s11ca.combine_bound_integrals


def _association(entries: dict[str, object]) -> C.Association:
    return C.Association(tuple(entries.items()))


NS = {name: getattr(sp, name) for name in dir(sp) if not name.startswith("_")}
NS.update(
    {
        "Association": _association,
        "TextAtom": C.TextAtom,
        "Tuple": sp.Tuple,
        "nan": sp.nan,
        "oo": sp.oo,
    }
)


@dataclass
class CapturedCase:
    key: str
    operand_a: str | None = None
    operand_b: str | None = None
    comparator_residual: str | None = None
    contexts: dict[str, str] = field(default_factory=dict)


def parse(raw: str) -> object | None:
    """Parse one comparator operand without interpreting measured content."""
    raw = raw.strip()
    if raw == "<MISSING>" or raw.startswith("<PARSE_FAILED"):
        return None
    value = eval(raw, {"__builtins__": {}}, NS)
    return _normalise_container(value)


def _normalise_container(value: object) -> object:
    if isinstance(value, C.Association):
        return C.Association(
            tuple((key, _normalise_container(item)) for key, item in value.entries)
        )
    if isinstance(value, sp.Tuple):
        return tuple(_normalise_container(item) for item in value)
    if isinstance(value, tuple):
        return tuple(_normalise_container(item) for item in value)
    return value


def _jet_target(name: str, active_map: dict[str, str]) -> str | None:
    if "XJETX" not in name:
        return None
    base, suffix = name.split("XJETX", 1)
    target = active_map.get(base)
    if target is None or base not in JET_APPLIED:
        return None
    return C.s11ca.canon_jet_name(target + "_" + suffix.replace("X", "_"))


def _reconcile_basic(expr: sp.Basic, active_map: dict[str, str], *, coefficients: bool = True) -> sp.Basic:
    if any(name.startswith(("gammaWJET", "gammaMUJET", "gamma_s11cb_")) for name in active_map):
        raise C.InputError("section 3a coefficients cannot enter the string rename table")
    substitutions: dict[sp.Basic, sp.Basic] = {}
    for applied in expr.atoms(AppliedUndef):
        name = applied.func.__name__
        if name in INERT_APPLIED:
            substitutions[applied] = sp.Symbol(name)
            continue
        if name.startswith(CURRENT_KEEP_W) and W in applied.free_symbols and applied.args != (W,):
            substitutions[applied] = sp.Function(name)(W)
            continue
        if name in PROFILE_JET_HEADS and all(item.is_Integer for item in applied.args):
            base = PROFILE_JET_HEADS[name]
            suffix = "".join(
                f"d{index + 1}" * int(order)
                for index, order in enumerate(applied.args)
            )
            substitutions[applied] = sp.Symbol(base + ("_" + suffix if suffix else ""))
            continue
        jet_target = _jet_target(name, active_map)
        if jet_target is not None:
            substitutions[applied] = sp.Symbol(jet_target)
        elif name in active_map and name in BARE_APPLIED:
            substitutions[applied] = sp.Symbol(active_map[name])
        elif name in active_map:
            substitutions[applied] = sp.Function(active_map[name])(*applied.args)
    for symbol in expr.atoms(sp.Symbol):
        if isinstance(symbol, sp.Dummy):
            continue
        jet_target = _jet_target(symbol.name, active_map)
        if jet_target is not None:
            substitutions[symbol] = sp.Symbol(jet_target)
        elif symbol.name in active_map:
            substitutions[symbol] = sp.Symbol(active_map[symbol.name])
    output = expr.xreplace(substitutions) if substitutions else expr
    return gamma_bridge().apply(output) if coefficients and any(
        atom.name.startswith(("gammaWJET", "gammaMUJET")) for atom in output.atoms(sp.Symbol)
    ) else output


def reconcile(expr: object, active_map: dict[str, str] | None = None, *, coefficients: bool = True) -> object:
    """Apply only the enumerated spelling and verified jet/applied rules."""
    names = WL_TO_PY_RENAME if active_map is None else active_map
    if isinstance(expr, C.Association):
        return C.Association(
            tuple((key, reconcile(item, names, coefficients=coefficients)) for key, item in expr.entries)
        )
    if isinstance(expr, tuple):
        return tuple(reconcile(item, names, coefficients=coefficients) for item in expr)
    if isinstance(expr, sp.MatrixBase):
        return expr.applyfunc(lambda item: _reconcile_basic(item, names, coefficients=coefficients))
    if isinstance(expr, sp.Basic):
        return _reconcile_basic(expr, names, coefficients=coefficients)
    return expr


def _residual_value(left: object, right: object) -> object:
    if isinstance(left, C.Association) and isinstance(right, C.Association):
        left_items, right_items = left.as_dict(), right.as_dict()
        if set(left_items) != set(right_items):
            return C.residual(left, right, name="S11c_b_handcoded", leaf_budget_seconds=0.1)
        residuals = tuple(
            (key, _residual_value(left_items[key], right_items[key]))
            for key in sorted(left_items)
        )
        return ResidualAssociation(residuals)
    if isinstance(left, tuple) and isinstance(right, tuple):
        if len(left) != len(right):
            return C.residual(left, right, name="S11c_b_handcoded", leaf_budget_seconds=0.1)
        return tuple(_residual_value(a, b) for a, b in zip(left, right))
    if isinstance(left, Relational) and isinstance(right, Relational):
        if type(left) is not type(right):
            return C.residual(left, right, name="S11c_b_handcoded", leaf_budget_seconds=0.1)
        return (_residual_value(left.lhs, right.lhs), _residual_value(left.rhs, right.rhs))
    if isinstance(left, sp.Basic) and isinstance(right, sp.Basic) and left == right:
        return sp.S.Zero
    if isinstance(left, sp.Expr) and isinstance(right, sp.Expr):
        difference = COMBINE_BOUND_INTEGRALS(sp.expand(left - right))
        if difference == 0:
            return sp.S.Zero
        # A symbol present on only one side is an exact nonzero witness under
        # the engines' common free-symbol semantics.  This avoids spending
        # minutes simplifying a measured residual that is already proven
        # nonzero, while retaining the complete expanded candidate residual.
        if left.free_symbols != right.free_symbols:
            return difference
        try:
            return sp.simplify(difference)
        except Exception:
            return difference
    return C.residual(left, right, name="S11c_b_handcoded", leaf_budget_seconds=0.1)


def _structural_zero(value: object) -> bool:
    if isinstance(value, ResidualAssociation):
        return all(_structural_zero(item) for _, item in value.entries)
    if isinstance(value, tuple):
        return all(_structural_zero(item) for item in value)
    return isinstance(value, sp.Basic) and value == 0


def residual_zero(
    a_raw: str, b_raw: str, active_map: dict[str, str] | None = None
) -> bool | None:
    """Parse, reconcile, and exactly test the A-minus-B residual for zero."""
    parsed_a, parsed_b = parse(a_raw), parse(b_raw)
    if parsed_a is None or parsed_b is None:
        return None
    residual = _residual_value(
        reconcile(parsed_a, active_map), reconcile(parsed_b, active_map)
    )
    return _structural_zero(residual)


def _capture_cases(text: str) -> list[CapturedCase]:
    records: list[CapturedCase] = []
    current: CapturedCase | None = None
    for line in text.splitlines():
        if line.startswith("CASE family="):
            if current is not None:
                records.append(current)
            marker = " key="
            if marker not in line:
                raise C.InputError(f"CASE line has no key: {line[:200]}")
            current = CapturedCase(line.split(marker, 1)[1])
        elif current is not None and line.startswith("operand_A = "):
            current.operand_a = line.removeprefix("operand_A = ")
        elif current is not None and line.startswith("operand_B = "):
            current.operand_b = line.removeprefix("operand_B = ")
        elif current is not None and line.startswith("A_minus_B = "):
            current.comparator_residual = line.removeprefix("A_minus_B = ")
        elif current is not None and line.startswith("context."):
            label, value = line.split(" = ", 1)
            current.contexts[label.removeprefix("context.")] = value
    if current is not None:
        records.append(current)
    for record in records:
        if record.operand_a is None or record.operand_b is None or record.comparator_residual is None:
            raise C.InputError(f"incomplete captured case {record.key}")
    return records


def _comparison_raws(record: CapturedCase) -> tuple[str, str]:
    # Coupling adjointness cases are already reduced modulo the one supplied
    # compact-support IBP quotient by C.  Compare those emitted signatures,
    # never the unreduced raw density and never add another collapse.
    py_reduced = record.contexts.get("PY_DIVERGENCE_REDUCED")
    wl_reduced = record.contexts.get("WL_DIVERGENCE_REDUCED")
    if py_reduced is not None or wl_reduced is not None:
        if py_reduced is None or wl_reduced is None:
            raise C.InputError(f"one-sided divergence-reduced context in {record.key}")
        return py_reduced, wl_reduced
    if record.operand_a is None or record.operand_b is None:
        raise C.InputError(f"missing captured operands in {record.key}")
    return record.operand_a, record.operand_b


def _candidate_residual(a_raw: str, b_raw: str, active_map: dict[str, str]) -> str:
    parsed_a, parsed_b = parse(a_raw), parse(b_raw)
    if parsed_a is None or parsed_b is None:
        return "TextAtom('UNDEFINED_UNJOINED')"
    value = _residual_value(reconcile(parsed_a, active_map), reconcile(parsed_b, active_map))
    return C.serialise(value)


def _disable_imported_prepass_for_drop(name: str) -> None:
    """Make an ablation visible when C inherited/pre-applied the same rename."""
    for table in (C.EXTRA_HEAD, C.EXTRA_SYMBOL, C.s11ca.PARAM, C.s11ca.FIELD, C.s11ca.PROFILE):
        table.pop(name, None)


def _run_family(
    family: str,
    py_tags: dict[str, str],
    wl_tags: dict[str, str],
    active_map: dict[str, str],
) -> None:
    if family == "ENERGY_BASIS_NEW_INVARIANTS":
        emit_gamma_certificate(gamma_bridge())
        return
    cases = C.extract_family(family, py_tags.get(family), wl_tags.get(family))
    capture = io.StringIO()
    with contextlib.redirect_stdout(capture):
        accounting = C.compare_family(family, cases, budget=0.1)
    if accounting.parse_failed or accounting.duplicate_key:
        raise C.InputError(
            f"{family}: parse_failed={accounting.parse_failed} "
            f"duplicate_key={accounting.duplicate_key}"
        )
    records = _capture_cases(capture.getvalue())
    if not records:
        raise C.InputError(f"{family}: comparator emitted no cases")

    one_sided: list[CapturedCase] = []
    differing: list[tuple[CapturedCase, str]] = []
    for record in records:
        if record.operand_a is None or record.operand_b is None:
            raise C.InputError(f"missing captured operands in {record.key}")
        if record.operand_a == "<MISSING>" or record.operand_b == "<MISSING>":
            one_sided.append(record)
            differing.append((record, record.comparator_residual or "TextAtom('UNDEFINED_UNJOINED')"))
            continue
        a_raw, b_raw = _comparison_raws(record)
        if residual_zero(a_raw, b_raw, active_map) is not True:
            differing.append((record, _candidate_residual(a_raw, b_raw, active_map)))

    if len(one_sided) == len(records):
        keys = ", ".join(record.key for record in one_sided)
        print(f"COVERAGE {family} ({len(records)}/{len(records)} cases one-engine-only — {keys})")
        return
    if not differing:
        print(f"MATCH {family}")
        return

    keys = ", ".join(record.key for record, _ in differing)
    print(f"FLAG {family} ({len(differing)}/{len(records)} cases differ — {keys})")
    for record, residual in differing:
        print(f"A_minus_B {family} {record.key} = {residual}")


# ---- P2b: immutable coefficient-only bridge, built from invariant objects ----
# These factors are supplied premises. No ratio or representative fold can
# enter this closed source-to-symbol map.
P2B_FACTORS = MappingProxyType({"W_BG": sp.Symbol("W_0"), "MU_R_BG": sp.Symbol("mu_R")})


def p2b_emit(tag: str, value: sp.Basic) -> None:
    """Only CAS payloads; tags name quantities, independently of their values."""
    if not isinstance(value, sp.Basic):
        raise TypeError(type(value).__name__)
    print(f"{tag} = {sp.srepr(value)}", flush=True)


@dataclass(frozen=True)
class GammaRecord:
    branch: str
    source: str
    coefficient: sp.Symbol
    invariant: sp.Expr
    action: sp.Expr


@dataclass(frozen=True)
class GammaPair:
    wl: GammaRecord
    py: GammaRecord

    @property
    def factor(self) -> sp.Symbol:
        return P2B_FACTORS[self.wl.source]

    @property
    def image(self) -> sp.Expr:
        # Reuse the expression-valued field_symbol mechanism, never Symbol(str).
        return C.s11ca.field_symbol(self.py.coefficient.name, scale=self.factor.name)


@dataclass(frozen=True)
class GammaBridge:
    pairs: tuple[GammaPair, ...]

    def apply(self, value: sp.Basic) -> sp.Basic:
        # Only coefficient atoms are substitution keys. Invariant expressions,
        # factors and background jets have no representable fold in this API.
        substitutions = {pair.wl.coefficient: pair.image for pair in self.pairs}
        return value.xreplace(substitutions)


_GAMMA_BRIDGE: GammaBridge | None = None


def gamma_bridge() -> GammaBridge:
    global _GAMMA_BRIDGE
    if _GAMMA_BRIDGE is None:
        install_gamma_bridge(C.load_py(C.DEFAULT_PY), C.load_wl(C.DEFAULT_WL))
    return _GAMMA_BRIDGE


def gamma_records(py_tags: dict[str, str], wl_tags: dict[str, str]) -> tuple[tuple[GammaRecord, ...], tuple[GammaRecord, ...]]:
    family = "ENERGY_BASIS_NEW_INVARIANTS"
    py_records, wl_records = [], []
    for branch_raw, wrapped in C.py_top_pairs(py_tags[family]):
        key = C.decode_py_key(branch_raw, ("BRANCH",))
        branch = C.key_dict(key)["BRANCH"]
        rows = C.parse_py_value(C.py_field(wrapped), key, coefficient=False)
        for row in rows:
            source, label, invariant, coefficient, term, dimensions = row
            source = source.value
            invariant = reconcile(invariant, coefficients=False)
            term = reconcile(term, coefficients=False)
            # The PY stream carries the quadratic wave bookkeeper in ENERGY_TERM.
            action = sp.expand(term).coeff(C.s11ca.EPSILON, 2)
            py_records.append(GammaRecord(branch, source, coefficient, invariant, action))
    for branch, rows in C.wl_assoc_pairs(wl_tags[family]):
        key = C.make_key((("BRANCH", branch),))
        for label, raw in C.wl_assoc_pairs(rows):
            source = next((source for prefix, source in (("WJET_", "W_BG"), ("MUJET_", "MU_R_BG")) if label.startswith(prefix)), None)
            if source is None:
                raise C.InputError(f"unrecognised invariant source {label}")
            values = {
                name: reconcile(C.parse_wl_value(C.wl_field(raw, name), key), coefficients=False)
                for name in ("INVARIANT", "COEFFICIENT", "ENERGY_TERM")
            }
            wl_records.append(GammaRecord(branch, source, values["COEFFICIENT"], values["INVARIANT"], values["ENERGY_TERM"]))
    return tuple(py_records), tuple(wl_records)


def pair_gamma_records(py_records: tuple[GammaRecord, ...], wl_records: tuple[GammaRecord, ...], *, tag: str = "P2B_PAIRING") -> GammaBridge:
    for engine, records in (("PY", py_records), ("WL", wl_records)):
        keys = [(row.branch, row.coefficient) for row in records]
        duplicate_count = sp.Integer(len(keys) - len(set(keys)))
        p2b_emit(f"{tag}_{engine}_duplicate_count", duplicate_count)
        if duplicate_count != 0:
            raise C.InputError("duplicate emitted coefficient key")
    pairs = []
    degrees = {id(row): 0 for row in py_records}
    for wl in sorted(wl_records, key=lambda row: (row.branch, row.source, row.coefficient.name)):
        factor = P2B_FACTORS[wl.source]
        candidates = tuple(py for py in py_records if py.branch == wl.branch and py.source == wl.source)
        residuals = tuple(sp.expand(py.invariant - factor * wl.invariant) for py in candidates)
        matches = tuple(py for py, residual in zip(candidates, residuals) if residual == 0)
        p2b_emit(f"{tag}_forward_degree_{wl.branch}_{wl.coefficient}", sp.Integer(len(matches)))
        for py, residual in sorted(zip(candidates, residuals), key=lambda item: item[0].coefficient.name):
            p2b_emit(f"{tag}_candidate_{wl.branch}_{wl.coefficient}_{py.coefficient}", sp.Tuple(wl.invariant, py.invariant, factor, residual))
        if len(matches) != 1:
            raise C.InputError("invariant pairing degree")
        degrees[id(matches[0])] += 1
        pairs.append(GammaPair(wl, matches[0]))
    for py in sorted(py_records, key=lambda row: (row.branch, row.coefficient.name)):
        p2b_emit(f"{tag}_reverse_degree_{py.branch}_{py.coefficient}", sp.Integer(degrees[id(py)]))
        if degrees[id(py)] != 1:
            raise C.InputError("reverse invariant pairing degree")
    for branch in sorted(C.BRANCHES):
        for source in sorted(C.SOURCES):
            counts = sp.Tuple(*(sp.Integer(sum(row.branch == branch and row.source == source for row in records)) for records in (wl_records, py_records)))
            p2b_emit(f"{tag}_source_counts_{branch}_{source}", counts)
            if counts != sp.Tuple(15, 15):
                raise C.InputError("source cardinality")
        counts = sp.Tuple(*(sp.Integer(sum(row.branch == branch for row in records)) for records in (wl_records, py_records)))
        p2b_emit(f"{tag}_branch_counts_{branch}", counts)
        if counts != sp.Tuple(30, 30):
            raise C.InputError("branch cardinality")
    counts = sp.Tuple(sp.Integer(len(wl_records)), sp.Integer(len(py_records)))
    p2b_emit(f"{tag}_overall_counts", counts)
    if counts != sp.Tuple(60, 60):
        raise C.InputError("overall cardinality")
    bridge = GammaBridge(tuple(pairs))
    by_name = {}
    for pair in bridge.pairs:
        previous = by_name.setdefault(pair.wl.coefficient, pair.image)
        residual = sp.expand(previous - pair.image)
        p2b_emit(f"{tag}_branch_image_residual_{pair.wl.branch}_{pair.wl.coefficient}", residual)
        if residual != 0:
            raise C.InputError("branch-dependent coefficient image")
    return bridge


def install_gamma_bridge(py_tags: dict[str, str], wl_tags: dict[str, str]) -> GammaBridge:
    global _GAMMA_BRIDGE
    _GAMMA_BRIDGE = pair_gamma_records(*gamma_records(py_tags, wl_tags))
    return _GAMMA_BRIDGE


def gamma_certificate_objects(bridge: GammaBridge) -> tuple[tuple[str, sp.Tuple], ...]:
    rows = []
    for pair in bridge.pairs:
        wl, py = pair.wl, pair.py
        image = reconcile(wl.coefficient)
        invariant_residual = sp.expand(py.invariant - pair.factor * wl.invariant)
        # Actions are the engines' emitted energy terms, not reconstructed answers.
        left, right = py.action, reconcile(wl.action)
        term_residual = sp.expand(left - right)
        # Also bind each emitted action to its emitted invariant by differentiation.
        py_action_residual = sp.expand(sp.diff(py.action, py.coefficient) - py.invariant)
        wl_action_residual = sp.expand(sp.diff(wl.action, wl.coefficient) - wl.invariant)
        factor_lock = sp.sympify(image == pair.image and isinstance(image, sp.Mul))
        rows.append((f"{wl.branch}_{wl.coefficient}", sp.Tuple(wl.coefficient, py.coefficient, wl.invariant, py.invariant, pair.factor, image, invariant_residual, left, right, term_residual, py_action_residual, wl_action_residual, factor_lock)))
    return tuple(rows)


def emit_gamma_certificate(bridge: GammaBridge) -> None:
    for name, row in gamma_certificate_objects(bridge):
        for label, value in zip(("wl_coefficient", "py_coefficient", "wl_invariant", "py_invariant", "family_factor", "coefficient_image", "invariant_residual", "py_action", "wl_action_image", "term_residual", "py_action_invariant_residual", "wl_action_invariant_residual", "factor_lock"), row):
            p2b_emit(f"P2B_CERTIFICATE_{name}_{label}", value)
        if any(row[index] != 0 for index in (6, 9, 10, 11)) or row[12] != sp.true:
            raise C.InputError("coefficient bridge certificate")


def p2b_parse_controls() -> None:
    # Action/ansatz construction is the sole input to these full-parser probes.
    coordinates = sp.symbols("xOne xTwo xThree time")
    shifted = sp.symbols("chiOne chiTwo chiThree time")
    probe = sp.Symbol("p2b_variation")
    key = C.make_key((("BRANCH", "LAB_HELD"),))
    for head in sorted(LIVE_BACKGROUND_APPLIED | INERT_APPLIED | {"thetaWave"}):
        # Use engine spellings at the WL boundary, not internal underscore names.
        branch = next((branch for branch, image in C.s11ca.MU_HEAD_BY_BRANCH.items() if image == head), "LAB_HELD")
        key = C.make_key((("BRANCH", branch),))
        input_head = "muThetaOperand" if head in INERT_APPLIED else head
        left_action = probe * sp.Function(input_head)(*coordinates)
        right_action = probe * sp.Function(input_head)(*shifted)
        # Mathematica printer supplies concrete WL syntax; parsing includes
        # canonical_value -> _extra_basic -> _reconcile_basic without bypasses.
        from sympy.printing.mathematica import mathematica_code
        left = reconcile(C.parse_wl_value(mathematica_code(sp.diff(left_action, probe)), key))
        right = reconcile(C.parse_wl_value(mathematica_code(sp.diff(right_action, probe)), key))
        full_residual = reconcile(C.parse_wl_value(mathematica_code(sp.diff(left_action - right_action, probe)), key))
        residual = sp.expand(full_residual)
        parse_linearity_residual = sp.expand(residual - (left - right))
        p2b_emit(f"P2B_PARSE_{head}_action", sp.Tuple(left_action, right_action))
        p2b_emit(f"P2B_PARSE_{head}_operands", sp.Tuple(left, right))
        p2b_emit(f"P2B_PARSE_{head}_residual", residual)
        p2b_emit(f"P2B_PARSE_{head}_linearity_residual", parse_linearity_residual)
        if parse_linearity_residual != 0:
            raise C.InputError("full residual parse linearity")
        if (head in LIVE_BACKGROUND_APPLIED and residual == 0) or (head not in LIVE_BACKGROUND_APPLIED and residual != 0):
            raise C.InputError("applied argument fixture")
    for head, classification in sorted(APPLIED_HEAD_CLASSES.items()):
        p2b_emit(f"P2B_HEAD_{head}_class", sp.Symbol(classification))
    coverage = sp.sympify(set(APPLIED_HEAD_CLASSES) == set(JET_APPLIED))
    p2b_emit("P2B_HEAD_classification_coverage", coverage)
    if coverage != sp.true:
        raise C.InputError("applied-head classification coverage")
    # Print the engine's Derivative form from the differentiated AST. SymPy's
    # generic Mathematica printer uses Hold[D[...]], outside the stream grammar.
    from sympy.printing.mathematica import MCodePrinter

    class JetPrinter(MCodePrinter):
        def _print_Derivative(self, expression):
            orders = dict(expression.variable_count)
            application = expression.expr
            counters = ", ".join(str(orders.get(arg, 0)) for arg in application.args)
            arguments = ", ".join(self._print(arg) for arg in application.args)
            return f"Derivative[{counters}][{application.func.__name__}][{arguments}]"

    printer = JetPrinter()
    # Both derivative decoding entry points remain live, including backgrounds.
    for head in sorted(LIVE_BACKGROUND_APPLIED):
        action = probe * sp.Function(head)(*coordinates)
        derivative = sp.diff(action, probe, coordinates[0])
        from sympy.printing.mathematica import mathematica_code
        parsed = reconcile(C.parse_wl_value(printer.doprint(derivative), key))
        direct = reconcile(C.s11ca.canonical_value(C.parse_wolfram_payload(C.preprocess_wl(printer.doprint(derivative))), "WL", "LAB_HELD"))
        p2b_emit(f"P2B_JET_{head}_action", action)
        p2b_emit(f"P2B_JET_{head}_operands", sp.Tuple(parsed, direct))
        residual = sp.expand(parsed - direct)
        p2b_emit(f"P2B_JET_{head}_residual", residual)
        p2b_emit(f"P2B_JET_{head}_atom_test", sp.sympify(isinstance(parsed, sp.Symbol) and isinstance(direct, sp.Symbol)))
        if residual != 0 or not isinstance(parsed, sp.Symbol) or not isinstance(direct, sp.Symbol):
            raise C.InputError("background jet decode fixture")


def p2b_action_controls(bridge: GammaBridge) -> None:
    # Mutations are made while constructing the control action. Every output
    # invariant/term is subsequently extracted by symbolic differentiation.
    for branch in sorted(C.BRANCHES):
        for source in sorted(C.SOURCES):
            pairs = tuple(pair for pair in bridge.pairs if pair.wl.branch == branch and pair.wl.source == source)
            first, second = pairs[:2]
            for package in ("pairing", "coefficient_scale", "family_factor", "string_image"):
                pair = first
                if package == "pairing":
                    action = second.wl.action.xreplace({second.wl.coefficient: pair.wl.coefficient})
                    factor = pair.factor
                    image_action = reconcile(action)
                else:
                    # The fold control changes the action normalization before
                    # extracting its invariant. Its energy cancellation alone
                    # cannot certify the supplied source factor.
                    action = pair.factor * pair.wl.action if package == "family_factor" else pair.wl.action
                    factor = sp.S.One if package == "family_factor" else pair.factor
                    if package == "string_image":
                        image = sp.Symbol(str(pair.image))
                    else:
                        image = pair.py.coefficient
                    # This is a control action with a deliberately different
                    # coefficient ansatz, never a mutation of a residual.
                    image_action = reconcile(action.xreplace({pair.wl.coefficient: image}))
                invariant = sp.diff(action, pair.wl.coefficient)
                invariant_residual = sp.expand(pair.py.invariant - factor * invariant)
                term_residual = sp.expand(pair.py.action - image_action)
                factor_lock = sp.sympify(factor == P2B_FACTORS[source])
                name = f"P2B_CONTROL_{package}_{branch}_{source}"
                for label, value in (("actions", sp.Tuple(pair.py.action, action, image_action)), ("invariants", sp.Tuple(pair.py.invariant, invariant)), ("family_factor", factor), ("invariant_residual", invariant_residual), ("term_residual", term_residual), ("factor_lock", factor_lock)):
                    p2b_emit(f"{name}_{label}", value)
                if package == "pairing" and invariant_residual == 0:
                    raise C.InputError("pairing control")
                if package != "pairing" and term_residual == 0 and factor_lock == sp.true:
                    raise C.InputError("coefficient control")


def p2b_run(py_path: Path, wl_path: Path) -> int:
    py_records, wl_records = gamma_records(C.load_py(py_path), C.load_wl(wl_path))
    global _GAMMA_BRIDGE
    _GAMMA_BRIDGE = pair_gamma_records(py_records, wl_records)
    bridge = _GAMMA_BRIDGE
    emit_gamma_certificate(bridge)
    baseline = gamma_certificate_objects(bridge)
    for engine, py, wl in (("PY", tuple(reversed(py_records)), wl_records), ("WL", py_records, tuple(reversed(wl_records)))):
        permuted = pair_gamma_records(py, wl, tag=f"P2B_PERMUTATION_{engine}")
        # Compare complete certificates and mappings, not only their zero tests.
        pairing_equal = sp.sympify(permuted == bridge)
        certificate_equal = sp.sympify(gamma_certificate_objects(permuted) == baseline)
        p2b_emit(f"P2B_PERMUTATION_{engine}_pairing_equality", pairing_equal)
        p2b_emit(f"P2B_PERMUTATION_{engine}_certificate_equality", certificate_equal)
        if pairing_equal != sp.true or certificate_equal != sp.true:
            raise C.InputError("record-order dependence")
    for engine, py, wl in (("PY", py_records + py_records[:1], wl_records), ("WL", py_records, wl_records + wl_records[:1])):
        rejected = False
        try:
            pair_gamma_records(py, wl, tag=f"P2B_DUPLICATE_{engine}")
        except C.InputError as error:
            if str(error) != "duplicate emitted coefficient key":
                raise
            rejected = True
        p2b_emit(f"P2B_DUPLICATE_{engine}_rejection", sp.sympify(rejected))
        if not rejected:
            raise C.InputError("duplicate fixture")
    p2b_action_controls(bridge)
    p2b_parse_controls()
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--py", type=Path, default=C.DEFAULT_PY)
    parser.add_argument("--wl", type=Path, default=C.DEFAULT_WL)
    parser.add_argument(
        "--drop-rename",
        metavar="WLname",
        choices=tuple(WL_TO_PY_RENAME),
        help="remove one verified spelling equivalence for a one-sided ablation",
    )
    parser.add_argument("--p2b", action="store_true", help="emit the section 3a bridge certificate and action controls")
    arguments = parser.parse_args(argv)
    if arguments.p2b:
        return p2b_run(arguments.py, arguments.wl)
    try:
        active_map = dict(WL_TO_PY_RENAME)
        if arguments.drop_rename is not None:
            active_map.pop(arguments.drop_rename)
            _disable_imported_prepass_for_drop(arguments.drop_rename)
        py_tags = C.load_py(arguments.py)
        wl_tags = C.load_wl(arguments.wl)
        install_gamma_bridge(py_tags, wl_tags)
        for family in CORE_FAMILIES:
            _run_family(family, py_tags, wl_tags, active_map)
        for family in NAMESPACE_INCOMPLETE_FAMILIES:
            print(
                f"NAMESPACE_INCOMPLETE {family} "
                "(WL operand unparsed; cross-engine control comparison owed; "
                "each engine's internal control verified in the build legs)"
            )
        print(f"RENAME_MAP_SIZE={len(active_map)}")
        return 0
    except Exception as error:
        print(
            f"OPERATIONAL_ERROR {type(error).__name__}: {error}",
            file=sys.stderr,
            flush=True,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
