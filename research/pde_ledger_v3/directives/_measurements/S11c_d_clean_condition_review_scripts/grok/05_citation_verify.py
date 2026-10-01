#!/usr/bin/env python3
"""Mechanical citation check: lookups excerpts vs source files vs A/B claims.
No physics predicates.
"""
from pathlib import Path

ROOT = Path("/var/projects/toy_physics")

def sed_n(path, start, end):
    lines = path.read_text().splitlines()
    # 1-indexed inclusive
    return "\n".join(lines[start - 1 : end])

checks = [
    (
        "V3_STEP_PLAN.md:1116-1126",
        ROOT / "research/pde_ledger_v3/V3_STEP_PLAN.md",
        1116,
        1126,
        "λγ = 1",
    ),
    (
        "assessment:11-15",
        ROOT / "research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md",
        11,
        15,
        "c_gamma/c_s=0.1",
    ),
    (
        "amendment:375-379",
        ROOT / "research/pde_ledger_v3/directives/S11c_d_SCATTERING_FORM_AMENDMENT.md",
        375,
        379,
        "v_bulk_normal_0=0",
    ),
    (
        "S11c_b spec:59-66",
        ROOT / "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md",
        59,
        66,
        "no w-component",
    ),
    (
        "S11c_b spec:87-92",
        ROOT / "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md",
        87,
        92,
        "appears in no derived operator",
    ),
    (
        "S11c_b spec:199-203",
        ROOT / "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md",
        199,
        203,
        "W_bg(y)",
    ),
    (
        "step record:32-36",
        ROOT / "research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md",
        32,
        36,
        "O(3)-Kronecker",
    ),
    (
        "ontology:26",
        ROOT / "docs/toy_model_ontology_summary.md",
        26,
        26,
        "trapped transverse brane-shear",
    ),
    (
        "ontology:362",
        ROOT / "docs/toy_model_ontology_summary.md",
        362,
        362,
        "helps hold the aperture open",
    ),
    (
        "ontology:100",
        ROOT / "docs/toy_model_ontology_summary.md",
        100,
        100,
        "localized throat drainage",
    ),
    (
        "ontology:1366",
        ROOT / "docs/toy_model_ontology_summary.md",
        1366,
        1366,
        "Gamma",
    ),
    (
        "native_light:113-116",
        ROOT / "docs/native_light_em_and_vortex_throat_interpretation.md",
        113,
        116,
        "trapped chiral shear",
    ),
    (
        "native_light:451",
        ROOT / "docs/native_light_em_and_vortex_throat_interpretation.md",
        451,
        451,
        "Spin-like behavior",
    ),
    (
        "native_light:1343",
        ROOT / "docs/native_light_em_and_vortex_throat_interpretation.md",
        1343,
        1343,
        "no comparable shear channel",
    ),
    (
        "assessment:5",
        ROOT / "research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md",
        5,
        5,
        "4097233",
    ),
    (
        "S11c_b spec:59-71 (B)",
        ROOT / "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md",
        59,
        71,
        "zeta_c is an independent face DOF",
    ),
    (
        "step record:15-31 (B)",
        ROOT / "research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md",
        15,
        31,
        "cross-engine-UNVALIDATED",
    ),
]

lookups = (ROOT / "research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_lookups.md").read_text()

print("=== FILE LINE CONTENT vs LOOKUPS CONTAINMENT ===")
for name, path, a, b, needle in checks:
    text = sed_n(path, a, b)
    has_needle = needle.lower() in text.lower() or needle in text
    # lookups may not contain B-only spans
    in_lookups = text.strip()[:80] in lookups or any(
        line in lookups for line in text.splitlines() if line.strip()
    )
    print(f"\n[{name}] path_ok={path.is_file()} needle={has_needle} some_line_in_lookups={in_lookups}")
    print("---FILE---")
    print(text)
    print("---END---")

# A's specific interpretive uses
print("\n=== A INTERPRETIVE USE vs CITED LINES ===")
amend = sed_n(
    ROOT / "research/pde_ledger_v3/directives/S11c_d_SCATTERING_FORM_AMENDMENT.md",
    375,
    379,
)
print("A claims: strict v_bulk_normal_0 = 0 (amendment :375-379)")
print("FILE 375-379:")
print(amend)

spec_p2 = sed_n(
    ROOT / "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md",
    59,
    66,
)
print("\nA P2 claims couplings through scalar face/bulk quantities, cites spec :59-66")
print("FILE 59-66:")
print(spec_p2)

record_p1 = sed_n(
    ROOT / "research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md",
    33,
    35,
)
print("\nA P1 cites record :33-35 for O(3)-Kronecker")
print("FILE 33-35:")
print(record_p1)

# engine 3-jet vs B coordinate claim
print("\n=== ENGINE BACKGROUND JETS (mechanical) ===")
eng = (ROOT / "research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py").read_text()
for s in ["w1_profile_d1", "w1_profile_d2", "w1_profile_d3", "DIRECTIONS = range(3)"]:
    print(f"engine_contains[{s}] = {s in eng}")
