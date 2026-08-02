#!/usr/bin/env python3
"""Show admitted definitions reduced to the quantities that must be supplied."""

from __future__ import annotations

from typing import Mapping

import sympy as sp

from registry_read import EvaluationError, Registry, Relation, load_registry


def admitted_definitions(registry: Registry) -> dict[str, Relation]:
    """Return only admission-gated, explicit output definitions."""
    return {
        relation.designated_output: relation
        for relation in registry.admitted_relations
        if relation.designated_output is not None and relation.rhs is not None
    }


def fully_reduced_outputs(
    registry: Registry,
    definitions: Mapping[str, Relation] | None = None,
) -> dict[str, tuple[sp.Expr, frozenset[str]]]:
    """Reduce admitted outputs by dependency traversal, with cycle/stall checks."""
    selected = dict(admitted_definitions(registry) if definitions is None else definitions)
    memo: dict[str, tuple[sp.Expr, frozenset[str]]] = {}
    visiting: list[str] = []

    def reduce_qid(qid: str) -> tuple[sp.Expr, frozenset[str]]:
        if qid in memo:
            return memo[qid]
        relation = selected.get(qid)
        if relation is None:
            return registry.symbols[qid], frozenset({qid})
        if qid in visiting:
            raise EvaluationError(
                f"definition cycle prevents full reduction: {visiting + [qid]}"
            )
        if relation.rhs is None:
            raise EvaluationError(f"{relation.relation_id}: admitted definition has no RHS")

        visiting.append(qid)
        substitutions: dict[sp.Symbol, sp.Expr] = {}
        for input_qid in relation.input_qids:
            reduced_input, _ = reduce_qid(input_qid)
            substitutions[registry.symbols[input_qid]] = reduced_input
        visiting.pop()

        reduced = sp.simplify(relation.rhs.xreplace(substitutions))
        stalled = sorted(
            output_qid
            for output_qid in selected
            if registry.symbols[output_qid] in reduced.free_symbols
        )
        if stalled:
            raise EvaluationError(
                f"{relation.relation_id}: substitution stalled on admitted outputs {stalled!r}"
            )
        symbol_to_qid = {symbol: symbol_qid for symbol_qid, symbol in registry.symbols.items()}
        unknown_symbols = reduced.free_symbols - set(symbol_to_qid)
        if unknown_symbols:
            raise EvaluationError(
                f"{relation.relation_id}: reduced expression has unknown symbols "
                f"{sorted(map(str, unknown_symbols))!r}"
            )
        leaves = frozenset(symbol_to_qid[symbol] for symbol in reduced.free_symbols)
        memo[qid] = reduced, leaves
        return memo[qid]

    for output_qid in selected:
        reduce_qid(output_qid)
    return memo


def main() -> int:
    registry = load_registry()
    definitions = admitted_definitions(registry)
    reduced = fully_reduced_outputs(registry, definitions)
    candidate_outputs = {
        relation.designated_output: relation
        for relation in registry.relations.values()
        if relation.designated_output is not None
    }

    print("=" * 78)
    print("MUST BE SUPPLIED / SELECTED  (no admitted defining equation)")
    print("=" * 78)
    simulation_inputs = 0
    structural_choices = 0
    for quantity in registry.quantities.values():
        if quantity.qid in definitions:
            continue
        if quantity.counting_axis == "discrete-structural":
            role = "structural selection"
            structural_choices += 1
        else:
            role = "simulation input"
            simulation_inputs += 1
        detail = ""
        candidate = candidate_outputs.get(quantity.qid)
        if candidate is not None:
            decision = registry.admission_decisions[candidate.relation_id]
            detail = f"  [candidate refused: {decision.reason}]"
        print(
            f"  {quantity.symbol_name:<14} dim={list(quantity.dimension)}  "
            f"[{quantity.kind}; {quantity.counting_axis}; {role}]{detail}"
        )

    print()
    print("=" * 78)
    print("DERIVED  (admitted defining equation -> computed, not supplied)")
    print("=" * 78)
    for quantity in registry.quantities.values():
        relation = definitions.get(quantity.qid)
        if relation is None:
            continue
        reduced_expression, leaves = reduced[quantity.qid]
        print(f"\n  {quantity.symbol_name}   [{relation.provenance_status}]")
        print(f"      as written : {quantity.symbol_name} = {relation.rhs}")
        if reduced_expression != relation.rhs:
            print(f"      reduced    : {quantity.symbol_name} = {reduced_expression}")
        rests_on = ", ".join(
            registry.quantities[qid].symbol_name for qid in sorted(leaves)
        )
        print(f"      rests on   : {rests_on if rests_on else '(pure number)'}")

    print()
    print("=" * 78)
    print(
        f"SUMMARY: {simulation_inputs} simulation inputs, "
        f"{structural_choices} structural "
        f"{'selection' if structural_choices == 1 else 'selections'}, "
        f"{len(definitions)} derived, {len(registry.quantities)} tracked total"
    )
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
