#!/usr/bin/env python3
"""Check whether each exact line-range citation is reproduced in the lookup file."""

from pathlib import Path

root = Path("/var/projects/toy_physics")
lookup = (root / "research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_lookups.md").read_text()

checks = [
    ("A", "research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md", 11, 13),
    ("A", "research/pde_ledger_v3/directives/S11c_d_SCATTERING_FORM_AMENDMENT.md", 375, 379),
    ("A", "research/pde_ledger_v3/V3_STEP_PLAN.md", 1116, 1126),
    ("A", "research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md", 33, 35),
    ("A", "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md", 59, 66),
    ("A", "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md", 90, 91),
    ("A", "docs/toy_model_ontology_summary.md", 26, 26),
    ("A", "docs/toy_model_ontology_summary.md", 362, 362),
    ("A", "docs/toy_model_ontology_summary.md", 100, 100),
    ("A", "docs/toy_model_ontology_summary.md", 1366, 1366),
    ("A", "docs/native_light_em_and_vortex_throat_interpretation.md", 113, 116),
    ("A", "docs/native_light_em_and_vortex_throat_interpretation.md", 451, 451),
    ("A", "docs/native_light_em_and_vortex_throat_interpretation.md", 1343, 1343),
    ("A", "research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md", 5, 5),
    ("B", "research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md", 5, 5),
    ("B", "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md", 59, 71),
    ("B", "research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md", 15, 31),
    ("B", "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md", 90, 91),
]

for packet, relative, start, end in checks:
    lines = (root / relative).read_text().splitlines()
    excerpt = "\n".join(lines[start - 1 : end])
    print(f"{packet}|{relative}:{start}-{end}|verbatim_in_lookup={excerpt in lookup}")

print("B_SECTION_CITATIONS_WITHOUT_VERBATIM_LOOKUP_RANGES=spec section 2a; spec section 3c")
