#!/usr/bin/env python3
"""Compare registry-derived medium dimensions with the existing fixture."""

from __future__ import annotations

import json

from registry_read import Registry, load_registry


# Copied verbatim only after the registry seed was loaded and its payload was
# computed.  This script never imports the old audit or its equation objects.
EXPECTED_MEDIUM_PAYLOAD = {
    "baseline": {"dim_before": 10, "dim_after": 5, "Delta": 5},
    "C-M1": {"dim_before": 10, "dim_after": 6, "Delta": 4},
    "C-M2": {"dim_before": 10, "dim_after": 5, "Delta": 5},
    "C-M3": {"dim_before": 10, "dim_after": 4, "Delta": 6},
}


def medium_cases(registry: Registry) -> dict[str, tuple]:
    """Recreate the old audit mutations on canonical registry expressions."""
    symbols = registry.symbols
    baseline = tuple(registry.admitted_constraint_set)
    r3 = registry.require_admitted("R3").residual
    assert r3 is not None
    without_r3 = tuple(expression for expression in baseline if expression != r3)
    lambda_gamma = symbols[registry.resolve_qid("lambda_gamma")]
    xi_h = symbols[registry.resolve_qid("xi_h")]
    a_pin = symbols[registry.resolve_qid("a")]
    big_k = symbols[registry.resolve_qid("K")]
    rho0 = symbols[registry.resolve_qid("rho0")]
    return {
        "baseline": baseline,
        "C-M1": without_r3 + (lambda_gamma - lambda_gamma,),
        "C-M2": baseline + (xi_h**2 - 2 * a_pin**2,),
        "C-M3": baseline + (big_k - rho0,),
    }


def compute_payload(registry: Registry) -> dict[str, dict[str, int]]:
    """Compute every number before consulting the expected payload."""
    ambient = len(registry.active_variables)
    result: dict[str, dict[str, int]] = {}
    for case_name, constraints in medium_cases(registry).items():
        after = registry.constraint_dimension(constraints)
        registry.certify_positive_real_dimension(constraints, dimension=after)
        result[case_name] = {
            "dim_before": ambient,
            "dim_after": after,
            "Delta": ambient - after,
        }
    return result


def main() -> int:
    registry = load_registry()
    observed = compute_payload(registry)

    print("OBSERVED_PAYLOAD:", json.dumps({"M": observed}, separators=(",", ":")))
    print(
        "EXPECTED_PAYLOAD:",
        json.dumps({"M": EXPECTED_MEDIUM_PAYLOAD}, separators=(",", ":")),
    )
    if observed != EXPECTED_MEDIUM_PAYLOAD:
        for case_name in EXPECTED_MEDIUM_PAYLOAD:
            if observed[case_name] != EXPECTED_MEDIUM_PAYLOAD[case_name]:
                print(
                    f"DISAGREEMENT {case_name}: "
                    f"observed={observed[case_name]} "
                    f"expected={EXPECTED_MEDIUM_PAYLOAD[case_name]}"
                )
        print("PHASE1_ACCEPTANCE: DISAGREEMENT_UNRECONCILED")
        return 1

    print("PHASE1_ACCEPTANCE: MATCH")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
