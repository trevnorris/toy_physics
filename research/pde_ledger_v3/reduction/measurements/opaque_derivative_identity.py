#!/usr/bin/env python3
"""Print identities, operands, residuals, equality, and hashes for five probes."""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from reduction.engine_output_checks import OpaqueDerivative, _map_tree


Identity = tuple[str, tuple[int, ...], str, tuple[str, ...]]
FIELDS = ("rendered", "orders", "function_name", "variables")


def identity(atom: OpaqueDerivative) -> Identity:
    return str(atom), atom.orders, atom.function_name, atom.variables


def residual(left: Identity, right: Identity) -> tuple[tuple[str, object, object], ...]:
    return tuple(
        (field, left_value, right_value)
        for field, left_value, right_value in zip(FIELDS, left, right)
        if left_value != right_value
    )


def print_probe(
    number: int,
    description: str,
    original: OpaqueDerivative,
    before: Identity,
    rebuilt: OpaqueDerivative,
) -> None:
    after = identity(original)
    rebuilt_identity = identity(rebuilt)
    print(
        f"PROBE {number}: {description} | identity_before={before!r} | "
        f"identity_after={after!r} | operand_original={after!r} | "
        f"operand_rebuilt={rebuilt_identity!r} | "
        f"residual_before_after={residual(before, after)!r} | "
        f"residual_original_rebuilt={residual(after, rebuilt_identity)!r} | "
        f"rebuilt_is_original={rebuilt is original} | "
        f"rebuilt_equals_original={rebuilt == original} | "
        f"hashes_equal={hash(rebuilt) == hash(original)}",
        flush=True,
    )


def main() -> int:
    original = OpaqueDerivative("probe1_rendered", (1, 0), "u", ("x", "t"))
    before = identity(original)
    mapped, _ = _map_tree(original, {"u": "v"})
    print_probe(1, "function-name rename u->v", original, before, mapped)

    original = OpaqueDerivative("probe2_rendered", (1, 0), "u", ("x", "t"))
    before = identity(original)
    mapped, _ = _map_tree(original, {"x": "y"})
    print_probe(2, "differentiation-variable rename x->y", original, before, mapped)

    original = OpaqueDerivative("probe3_rendered", (1, 0), "u", ("x", "t"))
    before = identity(original)
    mapped, _ = _map_tree(original, {"q": "r"})
    print_probe(3, "unrelated rename q->r", original, before, mapped)

    original = OpaqueDerivative("probe4_rendered", (1, 0), "u", ("x", "t"))
    _map_tree(original, {"u": "v"})
    before = identity(original)
    fresh = OpaqueDerivative("probe4_rendered", (1, 0), "u", ("x", "t"))
    print_probe(4, "fresh atom with the original identity", original, before, fresh)

    original = OpaqueDerivative("probe5_rendered", (1, 0), "u", ("x", "t"))
    before = identity(original)
    conflicting = OpaqueDerivative("probe5_rendered", (0, 1), "v", ("r", "s"))
    print_probe(5, "same rendered name and different attached identity", original, before, conflicting)
    return 0


if __name__ == "__main__":
    sys.exit(main())
