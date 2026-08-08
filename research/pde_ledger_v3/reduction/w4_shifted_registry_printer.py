#!/usr/bin/env python3
"""Drive the production pin printer with a valid shifted registry."""

from __future__ import annotations

from engine_dimension_pin import inspect_engine_dimension_pin, print_engine_dimension_pin
from registry_read import Registry, load_raw_documents


SHIFTED_QIDS = ("Q.brane.rho_br", "Q.brane.mu_R", "Q.brane.B_comp")


def _row(document: dict, qid: str) -> dict:
    return next(row for row in document["quantities"] if row["qid"] == qid)


def shifted_registry() -> Registry:
    quantities, relations, schema = load_raw_documents()
    for qid in SHIFTED_QIDS:
        dimension = _row(quantities, qid)["dimension"]
        dimension["exponents"][0] += 1
        component = dimension["law"]["components"][0]
        dimension["law"]["components"][0] = ["Add", component, 1]
    return Registry.from_documents(quantities, relations, schema)


def main() -> int:
    registry = shifted_registry()
    passed = print_engine_dimension_pin(inspect_engine_dimension_pin(registry), registry)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
