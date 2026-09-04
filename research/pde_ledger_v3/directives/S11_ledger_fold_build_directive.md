# Build directive — the ledger fold + under-export guard module

## Authority and boundary

Write `research/pde_ledger_v3/scripts/ledger_fold.py` in full, plus its self-test
`research/pde_ledger_v3/scripts/test_ledger_fold.py`. Those two files are the only writes.

`CLAUDE.md` binds. The physics-free **authority** is the committed design
`research/pde_ledger_v3/directives/export_ledger_bind_closure_design.md` (commit `c04e071f`), **§D2** (topology /
fold) and **§D3** (under-export guard), and the export-chain rules `S11_export_chain_decisions_v2.md` (**F9** —
write-time key routing; **F3** — a row carries its own evidence). This is **infrastructure, not physics**: it
manipulates LEDGER rows as opaque records; ⛔ it computes no physics, states no physics conclusion, and imports
no engine. Implement the contract below; where the design and this directive conflict, the **design wins**.

The row schema is the existing one (see any `scripts/*_exports.py`): a dict keyed by row-name, each value a dict
with `display`, `value` (a `_restore("<srepr>")` object), `value_kind`, `class`, `step`, optional
`dimension_key`, `f9_operands`, `corroborated_steps`. Treat rows as data.

## What to build — the contract (each function does exactly this; ⛔ add no physics, no expected value)

### 1 · `load_model(base_path, *delta_paths) -> Mapping[str, dict]`

Import the `LEDGER` of the base module and of each delta module (by file path; each is a `*_exports.py` exposing
`LEDGER`), and return the merged whole-model view:

- **Chronological last-wins on the exact key.** Apply the base first, then each delta in the given order; a
  later row with the same **exact key** replaces the earlier one. Preserve every non-colliding key.
- ⛔ **Do NOT re-apply F9.** The base is already F9-resolved internally and each delta already stores its rows
  under their final F9-resolved keys (bare on F9a/F9b, producer-prefixed on F9c — `F9`). The fold ⛔ never
  re-prefixes, re-compares, or re-routes; it merges by exact key only.
- **Route verification, not routing.** If a delta re-declares an exact key already present, that is an F9b
  overwrite the writer already decided — accept it but record the (key, prior `step`, new `step`) as a returned
  audit list. ⛔ Never silently pick between two *different* keys.
- Return the merged mapping plus an audit record (overwrites seen, per-source row counts). Deterministic and
  side-effect-free apart from the imports.

### 2 · The under-export guard (a consumer declares `IMPORT_KEYS`; ⛔ there is no `PASS_THROUGH`)

`check_consumer(fold, import_keys) -> report` must (raise on any failure, with the offending key(s) named):

- **(a) Manifest resolution.** Every key in `import_keys` is present in `fold` as an **exact key**. ⛔ Presence
  of a *bare* key does not satisfy a declared *prefixed* key or vice-versa — the declared key is the F9
  write-key and must match exactly (this is what closes the `face_response`/predecessor trap: a consumer that
  needs c1's curved response declares `s11c_c1_face_response`, and bare S11b `face_response` does not satisfy
  it).
- **(b) Recursive closure over every EDGE.** Compute the fixpoint closure of `import_keys` by, for each included
  row, adding (i) every **free-symbol name** of its `value` that names a row in `fold`, and (ii) its
  `dimension_key` target, recursively. Assert every closure member is present in `fold`.
- **(c) The F9c-ambiguity honesty check.** For any closure **symbol edge** whose symbol *name* matches **more
  than one** key in the fold (e.g. a bare `a1` row and a producer-prefixed `s11_a1` row both serialize
  `Symbol('a1')`), the srepr cannot disambiguate the intended producer — this is the unresolved routed-key
  contract (`S11_export_chain_decisions_v2.md:205`). ⛔ Do **not** silently pick one; **raise**, naming the
  ambiguous symbol and the competing keys, so the consumer must declare the intended write-key explicitly. (For
  a fold with no such collision the check is a no-op and must not fire.)
- **(d) Bidirectional smoke-test** `assert_lookups_equal_manifest(build_fn, fold, import_keys)`: wrap `fold` in
  an **access-recording proxy** (records every `__getitem__` key the consumer dereferences), run `build_fn`
  against the proxy, and assert the **set of recorded lookups equals `import_keys`** — an undeclared lookup
  (under-export risk) or a declared-but-unused key (stale manifest) each fails, naming the difference.

### 3 · Minimum-mode and the promotion-delta helper

- `assert_delta_is_minimal(delta_ledger, own_bind_closure, infra_keys=())`: the delta's exported key set equals
  `own_bind_closure` apart from explicitly named `infra_keys`; an accidental re-accumulation (extra rows) fails
  loudly, naming them.
- `promotion_delta(row_key, srepr, cls, evidence, ...) -> dict`: build a schema-complete row (the fields above,
  `value=_restore(srepr)`, `f9_operands`/route recorded) so an emit-only object can be **promoted** into the
  chain as a real delta row. ⛔ A manifest edit alone must not be treated as promotion; this helper (or a
  producer rebuild) is the only promotion path. Document that PY ⛔ never re-derives an inherited object at
  consume time (`N1`).

## `test_ledger_fold.py` — decisive tests that FAIL for the bad case (⛔ not tautologies)

Build small synthetic in-memory LEDGERs (a base + deltas) — ⛔ do not import a real engine export. Each test
must fail for the defect it guards and pass only when correct:

1. **Last-wins fold**: a delta re-declaring an exact key overwrites; a non-colliding delta row is added; the
   overwrite audit lists it. Corrupt to first-wins ⇒ test fails.
2. **No F9 re-apply**: a base row `x` and a delta row `s11c_c1_x` (a resolved F9c pair) both survive as distinct
   keys; the fold does not merge or re-prefix them.
3. **The predecessor trap**: a fold with base `face_response` (a `step:'S11b'` flat row) and delta
   `s11c_c1_face_response` (curved). `check_consumer(fold, {'face_response'})` must resolve to the S11b flat row
   (exact key) and `check_consumer(fold, {'s11c_c1_face_response'})` to the curved — declaring the wrong one
   silently binds the wrong object, so a test asserts the two are **distinct** objects and the manifest picks by
   exact key.
4. **Closure catches a missing symbol row**: a kept row whose `value` references `Symbol('foo')` where `foo` is
   **absent** from the fold ⇒ (b) raises naming `foo`. Add `foo` ⇒ passes. (This is the under-export D3-round-
   trip cannot catch.)
5. **Closure over `dimension_key`**: a kept row's `dimension_key` target absent ⇒ raises; present ⇒ passes.
6. **F9c ambiguity**: a fold with both `a1` and `s11_a1` (both serialize `Symbol('a1')`) and a kept row
   referencing `Symbol('a1')` ⇒ (c) raises naming both; remove the collision ⇒ no-op.
7. **Bidirectional smoke-test**: a `build_fn` that looks up an **undeclared** key ⇒ fails; one that omits a
   **declared** key ⇒ fails; exact match ⇒ passes.
8. **Minimum-mode**: a delta carrying an extra un-bound row ⇒ fails naming it; the exact bind-closure ⇒ passes.

Print each test's operands and the literal pass/fail per `CLAUDE.md` rule 2 (compute → emit → then assert);
exit nonzero only on a real failure. The build's §8-style report (≤ 25 lines) states files written, the
functions and tests, runtime, and any ambiguity or non-implementable clause.

## Digest / consumption note (for the consumers, not this module)

This module is a shared **executable input** to every consumer that imports it. Per the design build plan, a
consumer either inlines it or **adds `scripts/ledger_fold.py` to its `BUILD_INPUT_DIGESTS`** — ⭐ this directive
records that decision as **digest-pin the shared module** (do not duplicate it into each audit). `ledger_fold.py`
itself imports no export and pins nothing.

## Conflicts

Use the design (`c04e071f`) for any ambiguity; ⛔ do not fill a gap with physics, an expected value, or a
self-review mechanism. ⛔ Do not make a check audit the input it was built from. If a clause is not
implementable as stated (e.g. full F9c symbol routing without the routed-key contract), implement the honest
**raise-and-name** behaviour of (c) and report it in the build report rather than inventing a resolution.
