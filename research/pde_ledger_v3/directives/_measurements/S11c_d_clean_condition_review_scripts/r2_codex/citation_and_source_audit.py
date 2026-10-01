#!/usr/bin/env python3
"""Literal line audit for the v2 packet and its two grounding files."""

from pathlib import Path

ROOT = Path("/var/projects/toy_physics")


def lines(rel, start, end=None):
    data = (ROOT / rel).read_text().splitlines()
    end = start if end is None else end
    return "\n".join(f"{i}: {data[i-1]}" for i in range(start, end + 1))


records = [
    (
        "A_P1_AND_B_40_PINPOINT",
        "research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md",
        35,
        37,
    ),
    (
        "A_SPIN_SOURCE",
        "docs/native_light_em_and_vortex_throat_interpretation.md",
        2377,
        2384,
    ),
    (
        "B_WEAK_FORM_SOURCE",
        "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md",
        312,
        337,
    ),
    (
        "B_COMPARATOR_CONVENTION_SOURCE",
        "research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md",
        106,
        115,
    ),
    (
        "ENGINE_OPERATOR_ROWS",
        "research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py",
        3095,
        3127,
    ),
]

for label, rel, start, end in records:
    print(f"=== {label} ===")
    print(lines(rel, start, end))

print("=== PACKET_CLAIMS ===")
print(lines("research/pde_ledger_v3/directives/S11c_d_clean_condition.md", 83, 88))
print(lines("research/pde_ledger_v3/directives/S11c_d_clean_condition.md", 110, 123))
print(lines("research/pde_ledger_v3/directives/S11c_d_zinvariant_operator_blocks_directive.md", 72, 87))
print(lines("research/pde_ledger_v3/directives/S11c_d_zinvariant_operator_blocks_directive.md", 127, 131))

print("=== GROUNDING_LITERAL_SNIPPETS ===")
print(lines("research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition.md", 45, 55))
print(lines("research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition.md", 132, 141))
print(lines("research/pde_ledger_v3/directives/_measurements/S11c_d_zinvariant_operator_blocks_directive.md", 89, 113))

