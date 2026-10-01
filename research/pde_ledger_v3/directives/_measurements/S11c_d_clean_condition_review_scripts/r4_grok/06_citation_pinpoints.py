#!/usr/bin/env python3
"""Mechanical line-pinpoint checks for A/B citations against the source files."""
from pathlib import Path

ROOT = Path("/var/projects/toy_physics")


def lines(path, start, end):
    text = path.read_text().splitlines()
    return "\n".join(text[start - 1 : end])


checks = []

def check(name, path, start, end, must_contain):
    blob = lines(ROOT / path, start, end)
    missing = [s for s in must_contain if s not in blob]
    print(f"CHECK {name}")
    print(f"  RANGE {path}:{start}-{end}")
    for s in must_contain:
        print(f"  HAS {s!r}: {s in blob}")
    if missing:
        print(f"  MISSING {missing}")
    return blob

check("CHARTER_14_17", "research/pde_ledger_v3/CHARTER.md", 14, 17,
      ["REQUIREMENTS-FIRST", "falsification"])
check("SPEC_90_91", "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md", 90, 91,
      ["v_bulk_normal_0", "appears in no derived operator"])
check("SPEC_95_97", "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md", 95, 97,
      ["performs no curved-bulk response solve"])
check("SPEC_145_148", "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md", 145, 148,
      ["J_s = Λ_A", "t_s"])
check("SPEC_152_154", "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md", 152, 154,
      ["balance laws", "irreversible response kernel"])
check("RECORD_35_37", "research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md", 35, 37,
      ["O(3)-Kronecker", "40 = 10 uniform"])
check("U_214_219", "research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py", 214, 219,
      ['f"u_{a}"'])
check("ZETA_C_569_577", "research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py", 569, 577,
      ["zeta_c"])
check("WL_239_242", "research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl", 239, 242,
      ["spatialCoordinates = {xOne, xTwo, xThree}"])
check("ONTOLOGY_26", "docs/toy_model_ontology_summary.md", 26, 26,
      ["helps hold each throat open"])
check("ONTOLOGY_362", "docs/toy_model_ontology_summary.md", 362, 362,
      ["helps hold the aperture open"])
check("ONTOLOGY_315", "docs/toy_model_ontology_summary.md", 315, 315,
      ["not a complete nonlinear throat topology"])
check("ONTOLOGY_957", "docs/toy_model_ontology_summary.md", 957, 957,
      ["spectrally normalizable bound state"])
check("NATIVE_113_116", "docs/native_light_em_and_vortex_throat_interpretation.md", 113, 116,
      ["microrotation", "trapped chiral shear"])
check("NATIVE_1343", "docs/native_light_em_and_vortex_throat_interpretation.md", 1343, 1343,
      ["not sufficient by itself"])
check("NATIVE_2377_2381", "docs/native_light_em_and_vortex_throat_interpretation.md", 2377, 2381,
      ["zero time-averaged angular", "relative phase"])
check("ASSESS_3", "research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md", 3, 3,
      ["does not yet answer leakage in the calibrated, draining medium"])
check("ASSESS_5", "research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md", 5, 5,
      ["Do not launch another job alongside it"])
check("S11A_95_100", "research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py", 95, 100,
      ["delta_p_plus", "theta", "zeta_c"])
check("S11A_150_165", "research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py", 150, 165,
      ["delta_v_bulk", "d_w_delta_p_plus"])
check("AGENTS_8_12", "AGENTS.md", 8, 12,
      ["2 GiB memory cap", "s11c_guarded_run.py"])
check("AGENTS_24_33", "AGENTS.md", 24, 33,
      ["no wall-clock"])
check("V3_1107_1111", "research/pde_ledger_v3/V3_STEP_PLAN.md", 1107, 1111,
      ["λγ = 1", "calibrated / uncommitted"])
check("V3_1116_1126", "research/pde_ledger_v3/V3_STEP_PLAN.md", 1116, 1126,
      ["GW170817", "1 part in 10"])
check("ONTOLOGY_100", "docs/toy_model_ontology_summary.md", 100, 100,
      ["distributed return"])
check("ONTOLOGY_1366", "docs/toy_model_ontology_summary.md", 1366, 1366,
      ["converts ordered brane material"])
check("RECORD_112_114", "research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md", 112, 114,
      ["SURFACES them, does not normalize"])
check("DIRECTIONS_57", "research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py", 57, 57,
      ["DIRECTIONS = range(3)"])

# A's characterization: S11c-a 95-100 as trace definitions
blob = lines(ROOT / "research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py", 95, 100)
print("S11A_95_100_FULL")
print(blob)
print("S11A_95_100_TRACE_ONLY", blob.count("delta_p") == 2 and "delta_v_bulk" not in blob)
