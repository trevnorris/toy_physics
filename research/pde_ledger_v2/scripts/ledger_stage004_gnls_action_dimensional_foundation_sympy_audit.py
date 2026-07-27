#!/usr/bin/env python3
"""Ledger stage004 SymPy audit: GNLS action dimensional foundation.

Standalone, with audit results on stdout and labelled dimensions in a
deterministic sidecar.  This extracts the pathA_19 dimensional checks into
ledger form without importing the stage1 harness.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import sympy as sp

from ledger_dimensions import (
    Dimension,
    DimensionBasis,
    dim_residual,
    emit_dimension_sidecar,
)


PASS_COUNT = 0
FAIL_COUNT = 0


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


DIMENSION_BASIS = DimensionBasis("L", "T", "M", render="symbolic")
Dim = DIMENSION_BASIS
DIMENSIONLESS = Dim()
LENGTH = Dim(1, 0, 0)
TIME = Dim(0, 1, 0)
MASS = Dim(0, 0, 1)


def expect_dim(name: str, actual: Dimension, expected: Dimension) -> None:
    expect_zero(name, dim_residual(actual, expected))


def homogeneity_residual(terms: dict[str, Dimension]) -> sp.Expr:
    if not terms:
        raise AuditFailure("homogeneity check requires at least one term")
    dims = list(terms.values())
    reference = dims[0]
    return sp.simplify(sum(dim_residual(actual, reference) for actual in dims[1:]))


def factor_to_reach(expected: Dimension, actual: Dimension) -> Dimension:
    return expected / actual


def matrix_from_dims(dims: Iterable[Dimension]) -> sp.Matrix:
    return sp.Matrix([list(dim.components()) for dim in dims]).T


def signed_primitive_dictionary() -> dict[str, Dimension]:
    """Primitive input layer: posted, then checked through action usage."""

    return {
        "hbar": MASS * (LENGTH**2) / TIME,
        "m_GNLS": MASS,
        "rho0": LENGTH**-4,
    }


def derive_dictionary(
    primitives: dict[str, Dimension] | None = None,
) -> dict[str, Dimension]:
    """Derive the action dictionary from primitives and composition laws."""

    p = primitives or signed_primitive_dictionary()
    hbar = p["hbar"]
    m_gnls = p["m_GNLS"]
    rho0 = p["rho0"]

    action = hbar
    energy = action / TIME
    velocity = LENGTH / TIME
    action_density = energy / (LENGTH**4)

    # Kinetic density hbar*psi^2/T has the 4D action-density dimension.
    psi = (action_density / (hbar / TIME)) ** sp.Rational(1, 2)
    rho = psi**2
    k_eos = action_density / (rho**5)
    c_s0 = (k_eos * (rho0**4) / m_gnls) ** sp.Rational(1, 2)
    h0 = k_eos * (rho0**4)
    xi_h = ((hbar**2) / (m_gnls * h0)) ** sp.Rational(1, 2)

    force = energy / LENGTH
    q_a0 = energy
    q_ai = action / LENGTH
    electric_field = force
    magnetic_field = action / (LENGTH**2)
    lpsi_density = action_density
    maxwell_coeff = lpsi_density / (electric_field**2)

    return {
        **p,
        "1": DIMENSIONLESS,
        "L": LENGTH,
        "T": TIME,
        "M": MASS,
        "action": action,
        "energy": energy,
        "force": force,
        "velocity": velocity,
        "psi": psi,
        "rho": rho,
        "rho3_reduced": LENGTH**-3,
        "K": k_eos,
        "c_s0": c_s0,
        "h0": h0,
        "xi_h": xi_h,
        "lagrangian_density": lpsi_density,
        "q_A0": q_a0,
        "q_Ai": q_ai,
        "electric_field": electric_field,
        "magnetic_field": magnetic_field,
        "maxwell_coeff": maxwell_coeff,
        "mu_wall_restored_as_tau_over_c_s2": force / (velocity**2),
        "T_w": force,
        "U_Sigma_RR": force / (LENGTH**2),
        "G_3_spatial": Dim(3, -2, -1),
        "G_4_spatial": Dim(4, -2, -1),
    }


@dataclass(frozen=True)
class PinAnalysis:
    rank: int
    nullity: int
    relation_vector: tuple[sp.Integer, sp.Integer, sp.Integer, sp.Integer]
    relation_dim: Dimension


def normalize_null_vector(vector: sp.Matrix) -> tuple[sp.Integer, ...]:
    denominators = [sp.denom(component) for component in vector]
    scale = sp.ilcm(*[int(denom) for denom in denominators])
    entries = [sp.Integer(scale) * sp.Rational(component) for component in vector]
    gcd = abs(int(entries[0]))
    for entry in entries[1:]:
        gcd = sp.igcd(gcd, abs(int(entry)))
    entries = [sp.Integer(entry / gcd) for entry in entries]
    if entries[0] < 0:
        entries = [-entry for entry in entries]
    return tuple(entries)


def derive_pin_analysis(pin_dims: list[Dimension]) -> PinAnalysis:
    matrix = matrix_from_dims(pin_dims)
    nullspace = matrix.nullspace()
    if len(nullspace) != 1:
        return PinAnalysis(
            rank=int(matrix.rank()),
            nullity=len(pin_dims) - int(matrix.rank()),
            relation_vector=(sp.Integer(0), sp.Integer(0), sp.Integer(0), sp.Integer(0)),
            relation_dim=DIMENSIONLESS,
        )
    relation = normalize_null_vector(nullspace[0])
    relation_dim = DIMENSIONLESS
    for exponent, dim in zip(relation, pin_dims):
        relation_dim = relation_dim * (dim ** exponent)
    return PinAnalysis(
        rank=int(matrix.rank()),
        nullity=len(pin_dims) - int(matrix.rank()),
        relation_vector=relation,  # a, c_s0, hbar, m_GNLS
        relation_dim=relation_dim,
    )


def run_two_tier_dictionary_checks(d: dict[str, Dimension]) -> None:
    subbanner("Two-tier dictionary derivation")
    print("Primitive inputs (posted):")
    print(f"  [hbar] = {d['hbar']}  PRIMITIVE INPUT")
    print(f"  [m_GNLS] = {d['m_GNLS']}  PRIMITIVE INPUT")
    print(f"  [rho0] = {d['rho0']}  PRIMITIVE INPUT")
    print("Derived-by-composition dimensions:")
    print(f"  [psi] = {d['psi']}  from kinetic density")
    print(f"  [rho] = {d['rho']}  from psi^2")
    print(f"  [K] = {d['K']}  from EOS density")
    print(f"  [c_s0] = {d['c_s0']}  from 5*K*rho0^4/m_GNLS")
    print(f"  [h0] = {d['h0']}  from (5K/4)*rho0^4")
    print(f"  [xi_h] = {d['xi_h']}  from core balance")

    expect_dim("primitive hbar used as action in kinetic density", d["hbar"], d["action"])
    expect_dim("primitive m_GNLS used in Madelung velocity hbar/(m*L)", d["hbar"] / (d["m_GNLS"] * LENGTH), d["velocity"])
    expect_dim("primitive rho0 checked against derived rho=psi^2", d["rho"], d["rho0"])
    expect_dim("[psi] derived from action-density equation matches sqrt(rho0)", d["psi"], d["rho0"] ** sp.Rational(1, 2))
    expect_dim("[K] derived from EOS matches m*c_s0^2/rho0^4 target", d["K"], d["m_GNLS"] * (d["velocity"] ** 2) / (d["rho0"] ** 4))
    expect_dim("[c_s0] derived from sound-speed composition", d["c_s0"], d["velocity"])
    expect_dim("[h0] EOS scale matches m_GNLS*c_s0^2", d["h0"], d["m_GNLS"] * (d["c_s0"] ** 2))
    expect_dim("[xi_h] core-balance dimension", d["xi_h"], LENGTH)


def run_symbolic_core_relations() -> None:
    subbanner("Sound-speed and healing algebra")
    hbar, m, k_eos, rho0, c_s0, h0, xi_h = sp.symbols(
        "hbar m_GNLS K rho0 c_s0 h0 xi_h", positive=True
    )
    c_s0_sq = sp.solve(sp.Eq(c_s0**2, 5 * k_eos * rho0**4 / m), c_s0**2)[0]
    expect_zero("algebraic sound-speed law c_s0^2=5*K*rho0^4/m_GNLS", c_s0_sq - 5 * k_eos * rho0**4 / m)

    h0_from_eos = sp.Rational(5, 4) * k_eos * rho0**4
    h0_from_cs = sp.simplify(h0_from_eos.subs(5 * k_eos * rho0**4, m * c_s0**2))
    expect_zero("GNLS core balance h0=(5K/4)*rho0^4=(m_GNLS*c_s0^2)/4", h0_from_cs - m * c_s0**2 / 4)

    xi_solution = sp.solve(sp.Eq(xi_h**2, hbar**2 / (2 * m * h0)), xi_h)[0]
    xi_from_core = sp.simplify(xi_solution.subs(h0, m * c_s0**2 / 4))
    expect_zero("healing length xi_h=sqrt(2)*hbar/(m_GNLS*c_s0)", xi_from_core - sp.sqrt(2) * hbar / (m * c_s0))


def run_pin_analysis(d: dict[str, Dimension]) -> None:
    subbanner("Pin null-relation")
    pins = [LENGTH, d["c_s0"], d["hbar"], d["m_GNLS"]]
    analysis = derive_pin_analysis(pins)
    expect_zero("pin matrix rank is 3", sp.Integer(analysis.rank) - 3)
    expect_zero("four pins on three bases leave nullity 1", sp.Integer(analysis.nullity) - 1)
    expect_zero(
        "pin null vector is a*c_s0*hbar^-1*m_GNLS",
        sum((have - want) ** 2 for have, want in zip(analysis.relation_vector, (1, 1, -1, 1))),
    )
    expect_dim("derived pin relation a*c_s0*m_GNLS/hbar is dimensionless", analysis.relation_dim, DIMENSIONLESS)

    a, c_s0, hbar, m_gnls = sp.symbols("a c_s0 hbar m_GNLS", positive=True)
    relation_for_a = sp.solve(sp.Eq(a * c_s0 * m_gnls / hbar, 1), a)[0]
    expect_zero("derived pin relation a=hbar/(m_GNLS*c_s0)", relation_for_a - hbar / (m_gnls * c_s0))
    xi_h = sp.sqrt(2) * hbar / (m_gnls * c_s0)
    expect_zero("raw four pins give a/xi_h=1/sqrt(2)", relation_for_a / xi_h - 1 / sp.sqrt(2))
    print("  A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT: a is a branch collective moment, not a base invariant.")


# Harness mapping: the following 14 names preserve _patha19_foundation_checks
# load-bearing content, in order.  The 3 LT checks below preserve
# _patha19_lt_representation_checks.
def foundation_check_residuals(d: dict[str, Dimension]) -> dict[str, sp.Expr]:
    rho4 = d["rho0"]
    rho3 = d["rho3_reduced"]
    psi = d["psi"]
    velocity = d["velocity"]
    number_rate = TIME**-1
    energy = d["energy"]
    action = d["action"]
    lpsi_density = d["lagrangian_density"]
    electric_field = d["electric_field"]
    magnetic_field = d["magnetic_field"]
    maxwell_coeff = d["maxwell_coeff"]
    return {
        "pathA_19_F2: 4D-bulk closed-3-surface number flux J_bulk": dim_residual(
            rho4 * velocity * (LENGTH**3), number_rate
        ),
        "pathA_19_F2: 3D-brane reduced 2-sphere number flux J_brane": dim_residual(
            rho3 * velocity * (LENGTH**2), number_rate
        ),
        "pathA_19_F2: 4D-bulk volumetric flux Q_vol=rho^-1 J": dim_residual(
            (rho4**-1) * number_rate, (LENGTH**4) / TIME
        ),
        "pathA_19_F2: 3D-brane volumetric flux Q_vol=rho_3^-1 J": dim_residual(
            (rho3**-1) * number_rate, (LENGTH**3) / TIME
        ),
        "pathA_19_F2: constituent mass flux m_GNLS*J": dim_residual(
            d["m_GNLS"] * number_rate, MASS / TIME
        ),
        "pathA_19_F1: conditional defect rest-frequency conversion hbar*J/c_gamma^2": dim_residual(
            action * number_rate / (velocity**2), MASS
        ),
        "pathA_19_F2: bulk continuity equation": homogeneity_residual(
            {
                "partial_t rho": rho4 / TIME,
                "div_4(rho v)": (rho4 * velocity) / LENGTH,
            }
        ),
        "pathA_19_F3: sound-speed law 5*K*rho^4/m": dim_residual(
            d["K"] * (rho4**4) / d["m_GNLS"], velocity**2
        ),
        "pathA_19_F3: EOS enthalpy scale h0=(m_GNLS*c_s0^2)/4": dim_residual(
            d["m_GNLS"] * (velocity**2), energy
        ),
        "pathA_19_F3: GNLS healing length sqrt(hbar^2/(2*m_GNLS*h0))": dim_residual(
            (action**2 / (d["m_GNLS"] * energy)) ** sp.Rational(1, 2), LENGTH
        ),
        "pathA_19_F3: parent GNLS Lagrangian density terms": homogeneity_residual(
            {
                "i*hbar*psi*partial_t psi": action * (TIME**-1) * (psi**2),
                "hbar^2/(2m)*|D_i psi|^2": (action**2 / d["m_GNLS"]) * ((psi / LENGTH) ** 2),
                "V_conf*rho": energy * rho4,
                "U=K*rho^5/4": d["K"] * (rho4**5),
            }
        ),
        "pathA_19_F3: spatial gauge minimal-coupling dimension q*A_i/hbar": dim_residual(
            d["q_Ai"] / action, LENGTH**-1
        ),
        "pathA_19_F3: localized Maxwell sector with explicit c factors": homogeneity_residual(
            {
                "(Z/mu0)*E_i^2": maxwell_coeff * (electric_field**2),
                "(Z/mu0)*c^2*B_ij^2": maxwell_coeff * (velocity**2) * (magnetic_field**2),
                "A0*J0_ext": d["q_A0"] * rho4,
                "Ai*Ji_ext": d["q_Ai"] * (rho4 * velocity),
            }
        ),
        "pathA_19_F3: wall action density before dt*dw integration": homogeneity_residual(
            {
                "mu_eta*(partial_t eta)^2": d["mu_wall_restored_as_tau_over_c_s2"] * ((LENGTH / TIME) ** 2),
                "T_w*(partial_w eta)^2": d["T_w"],
                "K_eta*eta^2": d["U_Sigma_RR"] * (LENGTH**2),
            }
        ),
    }


def lt_representation_residuals() -> dict[str, sp.Expr]:
    rho = LENGTH**-4
    psi = LENGTH**-2
    velocity = LENGTH / TIME
    action_lt = Dim(2, -1, 0)
    energy_lt = Dim(2, -2, 0)
    mass_lt = DIMENSIONLESS
    force_lt = Dim(1, -2, 0)
    k_eos_lt = energy_lt / (rho**4)
    lpsi_density_lt = energy_lt / (LENGTH**4)
    electric_field_lt = force_lt
    magnetic_field_lt = action_lt / (LENGTH**2)
    maxwell_coeff_lt = lpsi_density_lt / (electric_field_lt**2)
    return {
        "pathA_19_LT_representation: local GNLS terms after projecting m_GNLS to dimensionless": homogeneity_residual(
            {
                "i*hbar*psi*partial_t psi": action_lt * (TIME**-1) * (psi**2),
                "hbar^2/(2m)*|D_i psi|^2": (action_lt**2 / mass_lt) * ((psi / LENGTH) ** 2),
                "V_conf*rho": energy_lt * rho,
                "U=K*rho^5/4": k_eos_lt * (rho**5),
            }
        ),
        "pathA_19_LT_representation: local Maxwell terms after M projection": homogeneity_residual(
            {
                "(Z/mu0)*E_i^2": maxwell_coeff_lt * (electric_field_lt**2),
                "(Z/mu0)*c^2*B_ij^2": maxwell_coeff_lt * (velocity**2) * (magnetic_field_lt**2),
            }
        ),
        "pathA_19_LT_representation: local wall terms after M projection": homogeneity_residual(
            {
                "mu_eta*(partial_t eta)^2": Dim(-1, 0, 0) * ((LENGTH / TIME) ** 2),
                "T_w*(partial_w eta)^2": force_lt,
                "K_eta*eta^2": Dim(-1, -2, 0) * (LENGTH**2),
            }
        ),
    }


@dataclass(frozen=True)
class FlaggedResidual:
    token: str
    actual: Dimension
    required_factor: Dimension


def flagged_residuals(d: dict[str, Dimension]) -> list[FlaggedResidual]:
    velocity = d["velocity"]
    formal_4d_target = d["G_4_spatial"] * (velocity**5) / ((LENGTH**5) * (velocity**5))
    observed_3d_target = d["G_3_spatial"] * (velocity**5) / ((LENGTH**5) * (velocity**5))
    lt_g3_target = Dim(3, -2, 0) * (velocity**5) / ((LENGTH**5) * (velocity**5))
    return [
        FlaggedResidual(
            "formal_4D_R_norm_target_not_dimensionless_without_conversion",
            formal_4d_target,
            factor_to_reach(DIMENSIONLESS, formal_4d_target),
        ),
        FlaggedResidual(
            "observed_3D_GR_target_not_dimensionless_without_conversion",
            observed_3d_target,
            factor_to_reach(DIMENSIONLESS, observed_3d_target),
        ),
        FlaggedResidual(
            "LT_R_norm_gate_fails_without_new_conversion_factor",
            lt_g3_target,
            factor_to_reach(DIMENSIONLESS, lt_g3_target),
        ),
    ]


def run_harness_mapped_checks(d: dict[str, Dimension]) -> None:
    subbanner("Harness-mapped 17 dimensional checks")
    for name, residual in foundation_check_residuals(d).items():
        expect_zero(name, residual)
    for name, residual in lt_representation_residuals().items():
        expect_zero(name, residual)
    expect_zero("harness check count is 14 foundation + 3 LT representation", len(foundation_check_residuals(d)) + len(lt_representation_residuals()) - 17)
    print("  hbar*J/c_gamma^2 = M check above is dimensional-only, not a mass derivation.")


def run_flagged_residual_gates(d: dict[str, Dimension]) -> None:
    subbanner("Carried flagged residual gates")
    for flag in flagged_residuals(d):
        print(f"  {flag.token}: actual {flag.actual}; factor needed {flag.required_factor}")
        expect_nonzero(
            f"{flag.token} remains non-dimensionless after M drop",
            dim_residual(flag.actual.without("M"), DIMENSIONLESS),
        )
        expect_dim(
            f"{flag.token} exact conversion factor",
            flag.required_factor * flag.actual,
            DIMENSIONLESS,
        )


def all_harness_checks_pass(d: dict[str, Dimension]) -> bool:
    return all(sp.simplify(residual) == 0 for residual in foundation_check_residuals(d).values()) and all(
        sp.simplify(residual) == 0 for residual in lt_representation_residuals().values()
    )


def lt_rejection_gate(
    d: dict[str, Dimension],
    *,
    repair_residuals: bool = False,
) -> bool:
    flags = flagged_residuals(d)
    if repair_residuals:
        residual_gate_fails = False
    else:
        residual_gate_fails = any(
            dim_residual(flag.actual.without("M"), DIMENSIONLESS) != 0
            for flag in flags
        )
    local_lt_passes = all(sp.simplify(residual) == 0 for residual in lt_representation_residuals().values())
    return all_harness_checks_pass(d) and local_lt_passes and residual_gate_fails


def verdict_from_predicate(
    d: dict[str, Dimension],
    *,
    m_defect_derived_here: bool,
    repair_lt_residuals: bool = False,
) -> str:
    retain_ltm = lt_rejection_gate(d, repair_residuals=repair_lt_residuals) and (not m_defect_derived_here)
    return "RETAIN_L_T_M" if retain_ltm else "NOT_RETAIN_L_T_M"


def run_verdict_and_provenance(d: dict[str, Dimension]) -> str:
    subbanner("Mass fork and computed verdict")
    m_defect_derived_here = False
    print(f"  m_defect_derived_here={m_defect_derived_here}")
    print("  INFLOW_MASS_SOURCE_MISSING: m_defect is NOT emergent at this foundation gate.")
    print("  hbar*J/c_gamma^2 = M is a dimensional conversion only, not a mass theorem.")
    verdict = verdict_from_predicate(d, m_defect_derived_here=m_defect_derived_here)
    expect_zero("computed verdict is RETAIN_L_T_M", 0 if verdict == "RETAIN_L_T_M" else 1)
    expect_bool(
        "verdict flips if LT residual gates are counterfactually repaired",
        verdict_from_predicate(d, m_defect_derived_here=False, repair_lt_residuals=True) != "RETAIN_L_T_M",
    )
    expect_bool(
        "verdict flips if m_defect_derived_here=True",
        verdict_from_predicate(d, m_defect_derived_here=True) != "RETAIN_L_T_M",
    )
    return verdict


def run_firewall_ablations(d: dict[str, Dimension]) -> None:
    subbanner("Able-to-fail dimensional firewall")
    psi = d["psi"]
    action = d["action"]
    energy = d["energy"]
    rho4 = d["rho0"]
    velocity = d["velocity"]

    expect_fail(
        "drop hbar from kinetic term breaks GNLS density homogeneity",
        homogeneity_residual(
            {
                "kinetic_without_hbar": (TIME**-1) * (psi**2),
                "gradient": (action**2 / d["m_GNLS"]) * ((psi / LENGTH) ** 2),
                "V_conf*rho": energy * rho4,
            }
        ),
    )

    wrong_k = Dim(17, -2, 1)
    expect_fail(
        "wrong K exponent M L^17 T^-2 breaks EOS density",
        homogeneity_residual(
            {
                "kinetic": action * (TIME**-1) * (psi**2),
                "U_wrong=K*rho^5/4": wrong_k * (rho4**5),
            }
        ),
    )

    expect_fail(
        "drop rho0^4 from c_s0^2=5*K*rho0^4/m_GNLS breaks velocity dimension",
        dim_residual((d["K"] / d["m_GNLS"]) ** sp.Rational(1, 2), velocity),
    )

    baseline_pin = derive_pin_analysis([LENGTH, d["c_s0"], d["hbar"], d["m_GNLS"]])
    dropped_pin_matrix = matrix_from_dims([LENGTH, d["hbar"], d["m_GNLS"]])
    expect_bool(
        "pin corruption by dropping c_s0 pin changes nullity",
        3 - int(dropped_pin_matrix.rank()) != baseline_pin.nullity,
    )
    corrupt_relation_dim = LENGTH * DIMENSIONLESS * (d["hbar"] ** -1) * d["m_GNLS"]
    expect_fail(
        "pin corruption c_s0 dimensionless breaks a=hbar/(m_GNLS*c_s0)",
        dim_residual(corrupt_relation_dim, DIMENSIONLESS),
    )


def print_verdict_labels(verdict: str) -> None:
    print("")
    print("Verdict labels:")
    print(
        "  ledger earned-label (NOT a source verdict token): "
        "DIMENSIONAL_FOUNDATION_LTM_RETAINED  "
        "(base {L,T,M}; a=hbar/(m*c_s0); xi_h=sqrt(2)*hbar/(m*c_s0))"
    )
    print(f"  source top-line verdict: {verdict}   (PASS_WITH_NAMED_RESIDUALS)")
    print(
        "  labeled non-derivation (carried gap): "
        "m_defect NOT emergent -- INFLOW_MASS_SOURCE_MISSING (m_defect_derived_here=False)"
    )
    print(
        "  carried residuals: LT_R_norm_gate_fails_without_new_conversion_factor "
        "(REJECTS_TRUE_LT_BASE); "
        "formal_4D_R_norm_target_not_dimensionless_without_conversion; "
        "observed_3D_GR_target_not_dimensionless_without_conversion"
    )
    print(
        "  carried forward (provenance, §3c): "
        "A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT; EOS_FROM_GNLS_FACTOR; "
        "NO_NET_ACCRETION_BC_UNDERIVED; M_TO_G_UNIFICATION; SCALE_MAP_INPUTS"
    )


def print_carried_gap_block() -> None:
    subbanner("Carried gaps printed as provenance")
    print("  A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT: a is branch geometry, not a base invariant.")
    print("  EOS_FROM_GNLS_FACTOR: h0=(m_GNLS*c_s0^2)/4 and xi_h=sqrt(2)*hbar/(m_GNLS*c_s0) carried to pathA_20.")
    print("  NO_NET_ACCRETION_BC_UNDERIVED: no-net-accretion is a boundary condition, not derived here.")
    print("  M_TO_G_UNIFICATION: defect-mass/back-reaction relation is deferred to pathA_21.")
    print("  SCALE_MAP_INPUTS: pathA_22 consumes J, a(branch), rho0, K, m_GNLS, hbar, 3D-reduction factors.")


def main() -> None:
    banner("ledger_stage004_gnls_action_dimensional_foundation SymPy audit")
    dictionary = derive_dictionary()

    run_two_tier_dictionary_checks(dictionary)
    run_symbolic_core_relations()
    run_pin_analysis(dictionary)
    run_harness_mapped_checks(dictionary)
    run_flagged_residual_gates(dictionary)
    verdict = run_verdict_and_provenance(dictionary)
    run_firewall_ablations(dictionary)
    print_carried_gap_block()
    print_verdict_labels(verdict)

    emit_dimension_sidecar(
        __file__,
        {
            "zero": dictionary["1"],
            "L": dictionary["L"],
            "T": dictionary["T"],
            "M": dictionary["M"],
            "action": dictionary["action"],
            "energy": dictionary["energy"],
            "force": dictionary["force"],
            "velocity": dictionary["velocity"],
            "rho3": dictionary["rho3_reduced"],
            "lagrangianDensity": dictionary["lagrangian_density"],
            "qA0": dictionary["q_A0"],
            "qAi": dictionary["q_Ai"],
            "electricField": dictionary["electric_field"],
            "magneticField": dictionary["magnetic_field"],
            "maxwellCoeff": dictionary["maxwell_coeff"],
            "muWall": dictionary["mu_wall_restored_as_tau_over_c_s2"],
            "Tw": dictionary["T_w"],
            "USigmaRR": dictionary["U_Sigma_RR"],
            "G3": dictionary["G_3_spatial"],
            "G4": dictionary["G_4_spatial"],
        },
    )
    print("")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print("OVERALL PASS: SymPy derived ledger_stage004 GNLS dimensional foundation exactly")


if __name__ == "__main__":
    main()
