#!/usr/bin/env python3
"""Mechanical citation check: A/B quoted phrases vs source files and grounding excerpts."""
from pathlib import Path

ROOT = Path("/var/projects/toy_physics")
A = (ROOT / "research/pde_ledger_v3/directives/S11c_d_clean_condition.md").read_text()
B = (ROOT / "research/pde_ledger_v3/directives/S11c_d_zinvariant_operator_blocks_directive.md").read_text()
GA = (ROOT / "research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition.md").read_text()
GB = (ROOT / "research/pde_ledger_v3/directives/_measurements/S11c_d_zinvariant_operator_blocks_directive.md").read_text()


def line(path, n):
    text = Path(path).read_text().splitlines()
    return text[n - 1] if 1 <= n <= len(text) else f"<missing line {n}>"


checks = []

def add(name, phrase, source_text, grounding_text=None):
    in_source = phrase in source_text
    in_ground = (phrase in grounding_text) if grounding_text is not None else None
    checks.append((name, phrase[:80], in_source, in_ground))


# A vs files
rec35 = line(ROOT / "research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md", 35)
add("A P1 Kronecker vs record:35 file", "O(3)-Kronecker field-bilinear invariant family", rec35)
add("A P1 Kronecker vs grounding excerpt", "O(3)-Kronecker field-bilinear invariant family", GA)

add("A assessment :3", "does not yet answer leakage in the calibrated, draining medium",
    line(ROOT / "research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md", 3),
    GA)
add("A spec :90-91 drain absent", "appears in no derived operator",
    line(ROOT / "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md", 91),
    GA)
add("A spec :97 no bulk solve", "S11c-b performs no curved-bulk response solve",
    line(ROOT / "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md", 97),
    GA)
add("A ontology :315 graphs", "not a complete nonlinear throat topology",
    line(ROOT / "docs/toy_model_ontology_summary.md", 315),
    GA)
add("A ontology :957 support mode", "spectrally normalizable bound state",
    line(ROOT / "docs/toy_model_ontology_summary.md", 957),
    GA)
add("A native :2377 heading", "A trapped standing wave is not automatically spinning",
    line(ROOT / "docs/native_light_em_and_vortex_throat_interpretation.md", 2377),
    GA)
add("A native :1343 shear channel", "de-structured bulk carries no comparable shear channel",
    line(ROOT / "docs/native_light_em_and_vortex_throat_interpretation.md", 1343),
    GA)
add("A native :113-116 microrotation", "independent\n  microrotation of the ordered substructure",
    "\n".join(
        Path(ROOT / "docs/native_light_em_and_vortex_throat_interpretation.md")
        .read_text().splitlines()[112:116]
    ),
    GA)
add("A ontology :26 trapped mode", "trapped transverse brane-shear standing mode helps hold each throat open",
    line(ROOT / "docs/toy_model_ontology_summary.md", 26),
    GA)
add("A ontology :362 structural support", "helps hold the aperture open",
    line(ROOT / "docs/toy_model_ontology_summary.md", 362),
    GA)
add("A ontology :1366 drain out", "converts ordered brane material into de-structured bulk material at throats",
    line(ROOT / "docs/toy_model_ontology_summary.md", 1366),
    GA)
add("A CHARTER.md in grounding", "CHARTER.md", GA)
add("A CHARTER requirements-first in CHARTER file", "REQUIREMENTS-FIRST",
    (ROOT / "research/pde_ledger_v3/CHARTER.md").read_text())

# B
add("B AGENTS 2 GiB", "whole-job 2 GiB memory cap",
    "\n".join(Path(ROOT / "AGENTS.md").read_text().splitlines()[7:12]),
    GB)
add("B spec :90-91", "appears in no derived operator",
    line(ROOT / "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md", 91),
    GB)
add("B record 40 terms vs file line 35", "40 = 10 uniform", rec35)
add("B record 40 terms vs grounding", "40 = 10 uniform", GB)
add("B kinetic sign convention", "kinetic −K PY vs +K WL",
    "\n".join(
        Path(ROOT / "research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md")
        .read_text().splitlines()[14:31]
    ),
    GB)
add("B af560257 exports", "af560257", GB)
add("B DIRECTIONS = range(3)", "DIRECTIONS = range(3)", GB)

print("=== citation checks (in_source, in_grounding) ===")
n_fail = 0
for name, phrase, in_source, in_ground in checks:
    status = f"source={in_source} ground={in_ground}"
    fail = (in_source is False) or (in_ground is False)
    mark = "FAIL" if fail else "ok  "
    if fail:
        n_fail += 1
    print(f"{mark} {name}: {status} :: {phrase!r}")
print(f"fail_count = {n_fail}")

print("\n=== record line 35 raw ===")
print(repr(rec35))
print("\n=== grounding A excerpt around record :35 ===")
idx = GA.find("sed -n 35p")
print(GA[idx:idx + 400] if idx >= 0 else "not found")
