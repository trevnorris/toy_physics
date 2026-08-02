#!/usr/bin/env python3
"""Compare registry-derived medium dimensions with this module's literal fixture."""

from __future__ import annotations

import json

import sympy as sp

from registry_read import Registry, load_registry


# INDEPENDENT CONTROL FIXTURE.  Provenance of the POST-S11 values, stated at the
# coverage each party actually earned -- an earlier wording credited all three
# with the full triple and was corrected after a review leg checked it:
#   (1) the orchestrator, PRE-REGISTERED and committed (67d919bd) before this
#       script existed -- covers baseline and C-M2 only (predicted 12->7 for
#       both, and EXPLICITLY DECLINED to predict C-M3);
#   (2) a fresh review agent, full triple, from quantities.yaml and
#       relations.yaml alone;
#   (3) a second independent reviewer, full triple, same registry-only route.
# So {12,7,5}/{12,7,5} has three independent derivations and {12,6,6} has two.
# A control's own provenance claim is subject to overstatement like any other
# measurement; state coverage, not a party count.
# On any registry change, RECOMPUTE and independently re-derive; never copy forward.
EXPECTED_MEDIUM_PAYLOAD = {
    "baseline": {"dim_before": 12, "dim_after": 7, "Delta": 5},
    "C-M2": {"dim_before": 12, "dim_after": 7, "Delta": 5},
    "C-M3": {"dim_before": 12, "dim_after": 6, "Delta": 6},
}


def medium_cases(registry: Registry) -> dict[str, tuple]:
    """Build the acceptance mutations on canonical registry expressions."""
    symbols = registry.symbols
    baseline = tuple(registry.admitted_constraint_set)
    xi_h = symbols[registry.resolve_qid("xi_h")]
    h0 = symbols[registry.resolve_qid("h0")]
    hbar = symbols[registry.resolve_qid("hbar")]
    mass = symbols[registry.resolve_qid("mass")]
    c_s0 = symbols[registry.resolve_qid("c_s0")]
    big_k = symbols[registry.resolve_qid("K")]
    rho0 = symbols[registry.resolve_qid("rho0")]
    xi_residual = registry.require_admitted("R2.xi_h").residual
    h0_residual = registry.require_admitted("R2.h0").residual
    assert xi_residual is not None and h0_residual is not None
    entailed = 2 * mass * h0 * xi_h**2 - hbar**2
    ideal_combination = (
        2 * mass * xi_h**2 * h0_residual
        + mass
        * c_s0
        * (mass * c_s0 * xi_h + sp.sqrt(2) * hbar)
        * xi_residual
        / 2
    )
    assert sp.expand(entailed) != 0
    assert sp.simplify(entailed - ideal_combination) == 0
    return {
        "baseline": baseline,
        "C-M2": baseline + (entailed,),
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
