#!/usr/bin/env python3
"""Mechanical citation/fold checks for the two v3 Markdown artifacts."""

from pathlib import Path
import subprocess


ROOT = Path("/var/projects/toy_physics")
A = ROOT / "research/pde_ledger_v3/directives/S11c_d_clean_condition.md"
B = ROOT / "research/pde_ledger_v3/directives/S11c_d_zinvariant_operator_blocks_directive.md"
GA = ROOT / "research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition.md"
GB = ROOT / "research/pde_ledger_v3/directives/_measurements/S11c_d_zinvariant_operator_blocks_directive.md"


def excerpt(relative, first, last=None):
    last = first if last is None else last
    rows = (ROOT / relative).read_text(encoding="utf-8").splitlines()
    return "\n".join(rows[first - 1:last])


checks = [
    ("A_CHARTER", "research/pde_ledger_v3/CHARTER.md", 14, 17, "REQUIREMENTS-FIRST"),
    ("A_ASSESSMENT_SCOPE", "research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md", 3, 3, "calibrated, draining medium"),
    ("A_ASSESSMENT_RATIOS", "research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md", 11, 13, "0.12247448713915873"),
    ("A_ASSESSMENT_PID", "research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md", 5, 5, "4097233"),
    ("A_PLAN_CLASS", "research/pde_ledger_v3/V3_STEP_PLAN.md", 1107, 1111, "calibrated / uncommitted"),
    ("A_PLAN_BOUND", "research/pde_ledger_v3/V3_STEP_PLAN.md", 1116, 1126, "1 part in 10¹⁵"),
    ("AB_RECORD_BASIS", "research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md", 35, 37, "40 = 10 uniform + 15"),
    ("AB_SPEC_DRAIN", "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md", 90, 91, "appears in no derived operator"),
    ("AB_SPEC_BULK", "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md", 95, 97, "no curved-bulk response solve"),
    ("AB_SPEC_FACE", "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md", 145, 148, "Λ_X"),
    ("B_SPEC_METHOD", "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md", 151, 154, "irreversible response kernel in an ordinary action"),
    ("B_SPEC_DOF", "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md", 59, 71, "ζ_c"),
    ("B_SPEC_BACKGROUND", "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md", 195, 196, "t_hold"),
    ("B_SPEC_WEAK_BLOCK", "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md", 312, 346, "S11CB_COUPLING_KERNEL"),
    ("B_RECORD_STATUS", "research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md", 15, 31, "cross-engine-UNVALIDATED"),
    ("B_RECORD_SIGNS", "research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md", 112, 114, "SURFACES them, does not normalize"),
    ("AB_FACE_LAWS", "research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md", 343, 354, "v_face"),
    ("AB_FACE_WORK", "research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md", 365, 366, "δ_vx"),
    ("A_ONTOLOGY_SUMMARY", "docs/toy_model_ontology_summary.md", 26, 26, "helps hold each throat open"),
    ("A_ONTOLOGY_DRAIN", "docs/toy_model_ontology_summary.md", 100, 100, "distributed return transfers material"),
    ("A_ONTOLOGY_GRAPHS", "docs/toy_model_ontology_summary.md", 315, 315, "not a complete nonlinear throat topology"),
    ("A_ONTOLOGY_SUPPORT", "docs/toy_model_ontology_summary.md", 362, 362, "helps hold the aperture open"),
    ("A_ONTOLOGY_MODE", "docs/toy_model_ontology_summary.md", 957, 957, "spectrally normalizable bound state"),
    ("A_ONTOLOGY_DRAIN_TERM", "docs/toy_model_ontology_summary.md", 1366, 1366, "de-structured bulk material"),
    ("A_NATIVE_FIELDS", "docs/native_light_em_and_vortex_throat_interpretation.md", 113, 116, "microrotation"),
    ("A_NATIVE_SHEAR", "docs/native_light_em_and_vortex_throat_interpretation.md", 1343, 1343, "not sufficient by itself"),
    ("A_NATIVE_SPIN", "docs/native_light_em_and_vortex_throat_interpretation.md", 2377, 2381, "relative phase"),
    ("B_ENGINE_DIRECTIONS", "research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py", 57, 57, "range(3)"),
    ("B_ENGINE_W1", "research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py", 190, 195, "w1_profile_d"),
    ("B_ENGINE_HIGHER_JET", "research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py", 420, 420, "w1_profile_d"),
    ("B_ENGINE_CENTER", "research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py", 2237, 2237, "CENTER_FACE_GENERALIZED_ROW"),
    ("B_WL_COORDS", "research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl", 241, 241, "xThree"),
]

ground_a = GA.read_text(encoding="utf-8")
ground_b = GB.read_text(encoding="utf-8")
ok = 0
for name, relative, first, last, token in checks:
    value = excerpt(relative, first, last)
    token_ok = token in value
    # Citations shared by A/B can occur in either grounding file.
    grounded = value in ground_a or value in ground_b
    if token_ok:
        ok += 1
    print("CITATION_CHECK", name, "TOKEN_OK", token_ok, "FULL_EXCERPT_GROUNDED", grounded)

print("CITATION_TOKEN_CHECK_COUNT", len(checks), "OK", ok)

engine = "research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py"
claimed_ranges = excerpt(engine, 190, 195) + "\n" + excerpt(engine, 420, 420)
actual_u_declaration = excerpt(engine, 214, 219)
print("B_INDEX_CLAIM_QUOTE", "symbol names such as u_1 and w1_profile_d1 (:190-195, :420)")
print("B_INDEX_CITED_TEXT", claimed_ranges.replace("\n", " | "))
print("B_INDEX_CITED_TEXT_CONTAINS_U_1", "u_1" in claimed_ranges)
print("B_ACTUAL_U_DECLARATION", actual_u_declaration.replace("\n", " | "))

charter_line14 = excerpt("research/pde_ledger_v3/CHARTER.md", 14, 14)
charter_context = excerpt("research/pde_ledger_v3/CHARTER.md", 14, 17)
print("A_CHARTER_CITED_LINE14", charter_line14)
print("A_CHARTER_REQUIRED_CONTEXT", charter_context.replace("\n", " | "))
print("A_GROUNDING_CONTAINS_REQUIRED_CONTEXT", charter_context in ground_a)

for commit in ("7b38e9dc", "e5f0dfda", "af560257"):
    result = subprocess.run(
        ["git", "log", "-1", "--format=%h %s", commit],
        cwd=ROOT,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    print("COMMIT_CHECK", commit, result.stdout.strip())

# Mechanical presence checks for the nine accepted round-2 folds.  These do
# not decide adequacy; they show what v3 actually contains.
a_text = A.read_text(encoding="utf-8")
b_text = B.read_text(encoding="utf-8")
fold_text = " ".join((a_text + "\n" + b_text).split())
fold_needles = {
    "R2-1": "symbol-name labels",
    "R2-2A": "full `O(2)`",
    "R2-2B": "tensor datum is invariant under rotations about direction 1",
    "R2-3": "emit **first**, for every entry, the raw `operand_PY`",
    "R2-4_CENTER": "CENTER_FACE_GENERALIZED_ROW",
    "R2-4_FACE": "components of `δ_v x_s`",
    "R2-5": "accepted 40-term stored-energy basis",
    "R2-6": "Two degenerate twist-type modes in quadrature",
    "R2-7": "**Nonlinear gates of R-LEAK-1**",
    "R2-8": "exponential for suitable analytic profiles",
    "R2-9": "record `:35–37`",
}
for name, needle in fold_needles.items():
    print("FOLD_TEXT_PRESENT", name, needle in fold_text, needle)
