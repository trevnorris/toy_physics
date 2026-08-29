#!/usr/bin/env python3
"""S11c-b hand-coded cross-engine reconciliation for the parsed core families.

The imported T7 comparator remains the case extractor and accounting authority.
This layer adds only the explicitly enumerated, source-verified spelling map
below, then observes whether each joined residual is zero.  It has no expected
zero/nonzero payload and never asserts on a measured result.

The energy-basis quotient representatives are compared exactly as emitted
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
    "gammaWidthUTheta": "gamma_s11cb_w_bg_04",  # WL:390-435; PY:291-303,1073-1090,1270-1283; spec:242-269 (W-jet u*theta invariant).
    "gammaWidthUEw": "gamma_s11cb_w_bg_05",  # WL:390-435; PY:291-303,1073-1090,1270-1283; spec:242-269 (W-jet u*e_W invariant).
    "gammaWidthUDiv": "gamma_s11cb_w_bg_01",  # WL:390-435; PY:291-303,1073-1090,1270-1283; spec:242-269 (W-jet u*div(u) invariant).
    "gammaWidthShearGradTheta": "gamma_s11cb_w_bg_06",  # WL:390-435; PY:291-303,1073-1090,1270-1283; spec:242-269 (W-jet grad(u)*grad(theta)).
    "gammaWidthShearGradEw": "gamma_s11cb_w_bg_09",  # WL:390-435; PY:291-303,1073-1090,1270-1283; spec:242-269 (W-jet grad(u)*grad(e_W)).
    "gammaWidthThetaGradEw": "gamma_s11cb_w_bg_13",  # WL:390-435; PY:291-303,1073-1090,1270-1283; spec:242-269 (W-jet theta*grad(e_W)).
    "gammaModulusUTheta": "gamma_s11cb_mu_r_bg_04",  # WL:390-435; PY:291-303,1073-1090,1270-1283; spec:242-269 (mu-jet u*theta invariant).
    "gammaModulusUEw": "gamma_s11cb_mu_r_bg_05",  # WL:390-435; PY:291-303,1073-1090,1270-1283; spec:242-269 (mu-jet u*e_W invariant).
    "gammaModulusUDiv": "gamma_s11cb_mu_r_bg_01",  # WL:390-435; PY:291-303,1073-1090,1270-1283; spec:242-269 (mu-jet u*div(u) invariant).
    "gammaModulusShearGradTheta": "gamma_s11cb_mu_r_bg_06",  # WL:390-435; PY:291-303,1073-1090,1270-1283; spec:242-269 (mu-jet grad(u)*grad(theta)).
    "gammaModulusShearGradEw": "gamma_s11cb_mu_r_bg_09",  # WL:390-435; PY:291-303,1073-1090,1270-1283; spec:242-269 (mu-jet grad(u)*grad(e_W)).
    "gammaModulusThetaGradEw": "gamma_s11cb_mu_r_bg_13",  # WL:390-435; PY:291-303,1073-1090,1270-1283; spec:242-269 (mu-jet theta*grad(e_W)).
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

# The WL heads below are evaluation syntax for PY's algebraic field/jet
# symbols.  Their derivatives are decoded first; only the remaining base
# applications lose inert coordinate arguments.  The cited engine definitions
# show the same independent trial/test/support object in both constructions.
BARE_APPLIED = {
    "widthBackground",
    "modulusBackground",
    "eWBackground",
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

# Intentionally absent from WL_TO_PY_RENAME:
#   * WL bRho is B_rho, whereas PY B_rho_3 == B_rho*W_0
#     (WL:471-479; PY:132-138,1133-1140; spec:98-103).
#   * WL gamma{Width,Modulus}DivGrad{Theta,Ew} multiply omitted PY
#     candidates 08/11, not the selected PY quotient representatives 07/10
#     (WL:390-435; PY:972-1003,1073-1090,1270-1298; spec:105-114,242-269).
# A scale substitution or quotient identity would be physics-bearing, so all
# of these names remain visible in the measured residuals.

# Inherited verbatim from the hand-verified S11c-a reconciliation: WL has no
# x-jets/Derivative of either branch operand, so its coordinate arguments are
# inert.  Current heads are different: w is physical and is always preserved;
# only the implicit x/time evaluation arguments are removed.
INERT_APPLIED = {"mu_theta_L", "mu_theta_M"}
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
    if target is None or base not in BARE_APPLIED:
        return None
    return C.s11ca.canon_jet_name(target + "_" + suffix.replace("X", "_"))


def _reconcile_basic(expr: sp.Basic, active_map: dict[str, str]) -> sp.Basic:
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
    return expr.xreplace(substitutions) if substitutions else expr


def reconcile(expr: object, active_map: dict[str, str] | None = None) -> object:
    """Apply only the enumerated spelling and verified jet/applied rules."""
    names = WL_TO_PY_RENAME if active_map is None else active_map
    if isinstance(expr, C.Association):
        return C.Association(
            tuple((key, reconcile(item, names)) for key, item in expr.entries)
        )
    if isinstance(expr, tuple):
        return tuple(reconcile(item, names) for item in expr)
    if isinstance(expr, sp.MatrixBase):
        return expr.applyfunc(lambda item: _reconcile_basic(item, names))
    if isinstance(expr, sp.Basic):
        return _reconcile_basic(expr, names)
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
    arguments = parser.parse_args(argv)
    try:
        active_map = dict(WL_TO_PY_RENAME)
        if arguments.drop_rename is not None:
            active_map.pop(arguments.drop_rename)
            _disable_imported_prepass_for_drop(arguments.drop_rename)
        py_tags = C.load_py(arguments.py)
        wl_tags = C.load_wl(arguments.wl)
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
