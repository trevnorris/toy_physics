# S11 SymPy engine — full-rewrite build directive

## Authority and boundary

Rewrite `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py` in full. The binding physics
is `git show cf4a21a4:research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`; `CLAUDE.md` also binds. The
shared specification wins every conflict, which is reported. The existing Python file is not a build
premise.

The script imports only `LEDGER` from `research/pde_ledger_v3/scripts/S10_exports.py`, and only for Q6r's
read-only comparison. It imports nothing from `reduction/`, creates no sidecar, writes no ledger or export,
and modifies no file. Its product is the flushed stdout tag stream.

## Required end state

Every object that the binding specification directs this engine to emit, wherever that object is named, is
computed and emitted at its specified scope for every declared cell. Q9 precedes action assembly at each
dimension. Build `PD_TERM` from the live, computed `V6_BASIS`, emit it for every `(package, D)`, and feed
that same `P_D` object into each action that carries it. Controls enter at the action. Every operand and
residual required by the specification is emitted before any guard. There is no `VERDICT`, `PASS`, `FAIL`,
checks list, or summary judgement. A physics finding exits zero; nonzero exit is reserved for an operational
failure.

The route-selection object that supplies the matrix consumed downstream also supplies `M_ROUTE_USED`.
Read both from that object; do not reconstruct the route token from a second literal.

Each emission is one unique line `TAG: <payload>`, immediately flushed. Payloads follow the shared
specification's re-parseable SymPy form; ordered records are SymPy tuples of `(field-name, value)` pairs in
the field order fixed by the specification. Tokens remain symbolic tokens.

The shared prefix is `PY_S11_`; engine-local tags use `PY_S11_LOCAL_`. Use every tag and quantity spelling
stated by the shared specification. Coefficient ordering is lexicographic by exact symbol name. A differing
computed root or stratum enumeration remains a finding; it is not padded with invented rows to force an
instantiated tag-set match.

The local tag `PY_S11_LOCAL_TAG_NAMES` lists every emitted local tag, including itself. Solver-returned
conditions are retained and are also emitted, unconditionally for each solver object, under that object's
local `_CONDITIONS` companion; an empty condition collection is still emitted. `PREMISE_INVENTORY` is the
single per-cell premise declaration required by the specification.

## Q6r

A coefficient is resolved exactly when the specification's name map and the live two-step `LEDGER` lookup
reach its dimension-row value. Provenance does not control resolution. Catch lookup failures before direct
indexing and preserve their distinct stages, in `COEFFICIENT_ORDERING`, in the local
`Q6R_LOOKUP_DIAGNOSTICS`: coefficient row absent; `dimension_key` absent; or referenced dimension-row/value
lookup failed. None stops the sweep or receives a placeholder vector. Derive
`Q6R_RESOLVED_COEFFICIENTS` and `Q6R_UNRESOLVED_COEFFICIENTS` from those lookup outcomes and preserve
`COEFFICIENT_ORDERING`.

Every Q6r tag is engine-local under `PY_S11_LOCAL_`, including both coefficient inventories,
`Q6R_LOOKUP_DIAGNOSTICS`, every per-entry tag, and `Q6R_RESIDUAL_SCOPE`. For entry `<q>` in
`Q6R_RESOLVED_COEFFICIENTS`, use the local suffixes
`Q6R_RESOLVED<q>_DERIVED_VECTOR`, `Q6R_RESOLVED<q>_IMPORTED_VECTOR`,
`Q6R_RESOLVED<q>_RESIDUAL`, `Q6R_RESOLVED<q>_COEFFICIENT_ROW_PROVENANCE`, and
`Q6R_RESOLVED<q>_DIMENSION_ROW_PROVENANCE`. `<q>` is one-based. The provenance payloads carry the fields
required by Q6r in its stated order; an absent provenance field is identified in that provenance payload
and neither removes the resolved coefficient nor suppresses its comparison. `Q6R_RESIDUAL_SCOPE` is emitted
once per cell. A disagreement is emitted and the run continues.

## Completion and gaps

Tags are emitted and flushed as their objects become available; observable progress is the runtime
requirement. Maintain a per-cell completion registry keyed by every emission required at that cell's
instantiated package, dimension, root, stratum, point-evidence, and engine-local scopes. Register fixed-scope
objects before the cell runs and register all required objects of a dynamic scope when that scope is
instantiated. An emission completes its registry key only after its line has flushed. Accumulate a pair in
`RUN_PAIRS` only after normal cell return and an empty required-minus-completed registry difference; merely
reaching the end of cell control flow is not completion. After the sweep, emit `RUN_PAIRS` and
`SKIPPED_PAIRS` in declared sweep order from these observed registry results, not from the declaration.

An empty, underdetermined, or undecided object is content: emit the object as built in the specification's
form, whether or not that object has a pinned status token. In particular, absent interface content does not
turn Q11's closure inventory into an operational failure. Never make emission depend on a constructed
payload's content.

If the supplied material does not define an object at all, or an operational failure prevents its
construction or emission, do not invent one. Leave its completion-registry key open, include the cell in
`SKIPPED_PAIRS`, continue where safe, and report the exact object and absent definition, input, or failed
operation under §10. This omission responds only to the missing construction, never to the value of a
constructed payload. Do not modify a committed output while verifying the script.
