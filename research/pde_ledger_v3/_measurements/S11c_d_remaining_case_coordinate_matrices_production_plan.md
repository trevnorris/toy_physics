# Separate accepted and growing material matrix inventories

The completed input focus's `case-inventory.json` is an immutable copied input.
The prepared production coordinator would otherwise save its growing numerical
case summaries at that same path. This was found before any production launch;
no integration or matrix was repeated or lost.

The separate production wrapper changes exactly one filename constant inside
the whole original `construct` function, from `case-inventory.json` to
`production-case-inventory.json`. Reversing that one constant must reproduce
the entire original coordinator AST. The original helper file, main bytecode,
native material accumulator, basis preparation and independent assembly bodies
remain unchanged. The wrapper calls the original accepted-focus loader first,
then adds its own current/frozen source/plan and the exact inventory join to the
manifest. Every accepted focus artifact is still copied and checked byte-for-byte.

The independent saved-focus validator checks this routing join without running
production, integration, assembly or any solver. Production may launch only
after focus acceptance, with this wrapper, the same required whole-job guard
and silent completion/error hook. All original numerical settings and scope
from the material matrix plan remain. The accepted focus inventory is retained
under its original name; the growing numerical inventory has a separate name.
