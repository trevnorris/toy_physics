#!/usr/bin/env python3
"""Report bound-law transport coverage and the unresolved coefficient witness."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from engine_dimension_pin import (
    EngineDimensionPin,
    inspect_engine_dimension_pin,
    print_engine_dimension_pin,
)
from dimension_laws import BoundDimensionLaw, SymbolicDimension, as_symbolic_dimension
from registry_read import Registry, load_registry


# These sets state which registry entries must use each declaration form.  They
# intentionally contain no expected dimension component or coefficient.
REQUIRED_BOUND_LAW_QIDS = frozenset(
    {"Q.brane.rho_br", "Q.brane.mu_R", "Q.brane.B_comp"}
)
REQUIRED_CONSTANT_QIDS = frozenset({"Q.brane.c_gamma", "Q.brane.c_L"})


@dataclass(frozen=True)
class DeclarationFormResult:
    required_bound_qids: frozenset[str]
    declared_bound_qids: frozenset[str]
    missing_bound_qids: frozenset[str]
    required_constant_qids: frozenset[str]
    declared_constant_qids: frozenset[str]
    missing_constant_qids: frozenset[str]
    passed: bool


@dataclass(frozen=True)
class SymbolicWitnessCoverage:
    registry_symbolic_qids: frozenset[str]
    engine_symbolic_qids: frozenset[str]
    missing_engine_symbolic_qids: frozenset[str]
    pin: EngineDimensionPin
    passed: bool


def _render(values: SymbolicDimension) -> str:
    return "[" + ",".join(sp.sstr(sp.simplify(value)) for value in values) + "]"


def inspect_declaration_forms(registry: Registry) -> DeclarationFormResult:
    """Compute whether the required entries use the required declaration forms."""
    declared_bound = frozenset(
        qid
        for qid in REQUIRED_BOUND_LAW_QIDS
        if isinstance(registry.dimension_declaration(qid), BoundDimensionLaw)
    )
    declared_constant = frozenset(
        qid
        for qid in REQUIRED_CONSTANT_QIDS
        if not isinstance(registry.dimension_declaration(qid), BoundDimensionLaw)
    )
    missing_bound = REQUIRED_BOUND_LAW_QIDS - declared_bound
    missing_constant = REQUIRED_CONSTANT_QIDS - declared_constant
    passed = not missing_bound and not missing_constant
    return DeclarationFormResult(
        REQUIRED_BOUND_LAW_QIDS,
        declared_bound,
        missing_bound,
        REQUIRED_CONSTANT_QIDS,
        declared_constant,
        missing_constant,
        passed,
    )


def inspect_symbolic_witness_coverage(registry: Registry) -> SymbolicWitnessCoverage:
    """Compute coefficient coverage from committed engine output operands."""
    registry_symbolic = frozenset(
        qid
        for qid in REQUIRED_BOUND_LAW_QIDS
        if isinstance(registry.dimension_declaration(qid), BoundDimensionLaw)
    )
    pin = inspect_engine_dimension_pin(registry)
    engine_symbolic = pin.covered_qids
    missing = registry_symbolic - engine_symbolic
    passed = (
        not missing
        and registry_symbolic == REQUIRED_BOUND_LAW_QIDS
        and pin.passed
    )
    return SymbolicWitnessCoverage(
        registry_symbolic, engine_symbolic, missing, pin, passed
    )


def print_dimension_law_check(
    registry: Registry,
    forms: DeclarationFormResult,
    coverage: SymbolicWitnessCoverage,
) -> bool:
    """Print computed operands and residuals, then return the combined guard."""
    for qid in sorted(forms.declared_bound_qids):
        declaration = registry.dimension_declaration(qid)
        print(
            f"DECLARED_BOUND_LAW {qid}: "
            f"components={_render(as_symbolic_dimension(declaration))} "
            f"bindings={declaration.binding_map}"
        )
    print(f"DECLARATION_FORM_OPERAND required-bound: {sorted(forms.required_bound_qids)}")
    print(f"DECLARATION_FORM_OPERAND declared-bound: {sorted(forms.declared_bound_qids)}")
    print(f"DECLARATION_FORM_RESIDUAL missing-bound: {sorted(forms.missing_bound_qids)}")
    print(
        f"DECLARATION_FORM_OPERAND required-constant: "
        f"{sorted(forms.required_constant_qids)}"
    )
    print(
        f"DECLARATION_FORM_OPERAND declared-constant: "
        f"{sorted(forms.declared_constant_qids)}"
    )
    print(
        f"DECLARATION_FORM_RESIDUAL missing-constant: "
        f"{sorted(forms.missing_constant_qids)}"
    )
    print(
        f"D_COEFFICIENT_WITNESS_OPERAND registry-symbolic-qids: "
        f"{sorted(coverage.registry_symbolic_qids)}"
    )
    print(
        f"D_COEFFICIENT_WITNESS_OPERAND engine-symbolic-qids: "
        f"{sorted(coverage.engine_symbolic_qids)}"
    )
    print(
        f"D_COEFFICIENT_WITNESS_RESIDUAL missing-engine-qids: "
        f"{sorted(coverage.missing_engine_symbolic_qids)}"
    )
    print_engine_dimension_pin(coverage.pin, registry)
    print(f"DECLARATION_FORM_CHECK: {'PASS' if forms.passed else 'FAIL'}")
    print(f"D_COEFFICIENT_POLICED_IN_REDUCTION: {'YES' if coverage.passed else 'NO'}")
    passed = forms.passed and coverage.passed
    print(f"DIMENSION_LAW_CHECK: {'PASS' if passed else 'FAIL'}")
    return passed


def main() -> int:
    registry = load_registry()
    forms = inspect_declaration_forms(registry)
    coverage = inspect_symbolic_witness_coverage(registry)
    passed = print_dimension_law_check(registry, forms, coverage)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
