#!/usr/bin/env python3
"""List every emitted payload that is a bare word rather than a CAS object.

The defect the rebuild exists to remove is an engine printing a physics
conclusion as a typed sentence with no computed object behind it.  A typed
conclusion is a bare word; a CAS object almost always carries an operator, a
bracket, a digit, or a symbol name built from the step's own vocabulary.

This script partitions every emitted payload in both engines' committed outputs
into CAS-looking and bare-word, and prints the bare-word ones grouped by value
with their counts and one example tag each.  It classifies nothing as good or
bad; the reader decides which of the listed words are computed CAS booleans and
which are typed verdicts.
"""

from __future__ import annotations

import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUTPUTS = {
    "wl": ROOT / "mathematica" / "out" / "S10_brane_mode_spectrum_mathematica_audit.out",
    "py": ROOT / "scripts" / "out" / "S10_brane_mode_spectrum_sympy_audit.out",
}

TAG = re.compile(r"^[A-Z0-9_]+$")
# A bare word: letters and underscores only, no operator, bracket, digit or space.
BARE_WORD = re.compile(r"^[A-Za-z_]+$")


def main() -> int:
    for engine, path in OUTPUTS.items():
        values: Counter[str] = Counter()
        examples: dict[str, str] = {}
        by_tag_suffix: defaultdict[str, Counter[str]] = defaultdict(Counter)
        total = 0
        with path.open() as handle:
            for line in handle:
                head, separator, tail = line.partition(": ")
                if not separator or not TAG.fullmatch(head):
                    continue
                total += 1
                value = tail.strip()
                if BARE_WORD.fullmatch(value):
                    values[value] += 1
                    examples.setdefault(value, head)
                    by_tag_suffix[value][head.rsplit("_", 2)[-1]] += 1
        print(f"ENGINE={engine} SOURCE={path}")
        print(f"  total_payloads={total} bare_word_payloads={sum(values.values())}")
        print(f"  distinct_bare_words={len(values)}")
        for value, count in values.most_common():
            suffixes = ", ".join(
                f"{suffix}x{n}" for suffix, n in by_tag_suffix[value].most_common(4)
            )
            print(f"    {value!r}: count={count} example_tag={examples[value]} suffixes=[{suffixes}]")
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
