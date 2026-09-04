# Repair directive — ledger_fold guard (fold both review legs, one pass)

Repair `research/pde_ledger_v3/scripts/ledger_fold.py` and `research/pde_ledger_v3/scripts/test_ledger_fold.py`
(baseline committed `05f86947`). Two independent legs (fresh Claude + Grok, saved under
`directives/_legs/ledger_fold_review_{agent,grok}.out`) found the guard **not sound as-is**. Authority unchanged:
the design `export_ledger_bind_closure_design.md` (`c04e071f`) §D2/§D3. ⛔ Physics-free infra; ⛔ do not weaken any
passing behavior (last-wins fold, exact-key manifest, the four decisive tests). Fix exactly the following.

## Fix 1 — the closure must follow EVERY edge type the ledger uses, not only bare `Symbol`

`_free_symbol_names` (`ledger_fold.py:146-161`) collects only `sp.Symbol` atoms of `value.free_symbols`. But the
ledger represents **fields/profiles as undefined functions** (`Function('O_window')`, `u1..u5`,
`delta_j_bulk_4`, …) and may use `MatrixSymbol`. So a kept row that applies an **absent** profile-function
passes the guard — the guard cannot enforce its own §D1 clause "every symbol/knob/**profile** referenced in a
kept row's value". Extend the referenced-name extraction to the **union** of, for a sympy `Basic` value:
`value.atoms(Symbol)` names, `{f.func.__name__ for f in value.atoms(AppliedUndef)}`, and
`value.atoms(MatrixSymbol)` names — then (as now) intersect with fold keys to form closure edges. For a **non-`Basic`
payload** (a Python `str`, sympy `Str`, a bare container of strings) from which names cannot be extracted, ⛔ do
not silently treat it as edge-free if it could hide a dependency: document the limitation and, if the row is in
a consumer's closure and its payload is a non-walkable type, **raise a typed error naming the row** so the
consumer must resolve it explicitly (the honest raise-and-name, matching the F9c case). A row whose payload is a
plain declaration (`Symbol(...)`) remains a leaf.

## Fix 2 — the F9c-ambiguity check must detect EXPRESSION-valued collisions, not only `Symbol`-valued

`_symbol_rows` (`:169-174`) indexes a write-key only if its `value` is `isinstance(..., sp.Symbol)`, and
candidate assembly (`:210-213`) is "that index ∪ exact bare key". So an F9c pair of **computed** rows
(`face_response` vs `s11c_c1_face_response`, or any re-derived expression) never enters the index: a closure edge
naming `Symbol('a1')` resolves to the **bare predecessor alone**, silently, when the prefixed producer is also
present. Generalize the ambiguity detection so it does not depend on the colliding rows' value type: for a
closure edge name `n`, determine the set of fold keys that **represent object `n`** — the bare key `n` plus any
producer-routed key for `n` — using the row **route metadata** (`f9_operands` / the writer-stored route the fold
carries), ⛔ not a Symbol-value match. If that set has **more than one** member, **raise-and-name** (the routed-key
contract is unresolved from srepr — `S11_export_chain_decisions_v2.md:205`); a single member resolves; zero is
the missing-row error of Fix 1. Where the route metadata cannot decide, ⛔ raise rather than pick.

## Fix 3 — make the recursion DECISIVE in the tests (both legs: one-level ablation ships green)

Deleting the recursion (`pending.append` at `:228` and `:239`) currently leaves all 9 tests passing. Add tests
that FAIL for a one-level closure:
- a **depth-2 symbol chain** (`root` → `Symbol('mid')` row → `Symbol('leaf')` row) whose **leaf row is absent**
  must raise; present ⇒ passes;
- a **depth-2 `dimension_key` chain** (`root.dimension_key='mid'`, `mid.dimension_key='leaf'`) with the leaf
  absent must raise;
- a **function-edge** test (Fix 1): a kept row whose value is `Function('prof')(x)` with the `prof` row **absent**
  must raise; add `prof` ⇒ passes; and a `MatrixSymbol` edge likewise;
- an **expression-valued F9c collision** test (Fix 2): bare `a1` and `s11c_x_a1` both present as **expression**
  rows both representing object `a1`, a kept row naming `Symbol('a1')` ⇒ raise-and-name; remove the routed key ⇒
  resolves.
Each test prints operands + observed + literal pass/fail (rule 2), and must fail on the corresponding one-line
ablation. ⛔ Do not add a test that passes on the broken implementation.

## Fix 4 (nits — low priority, do if cheap)

Document in `_AccessRecordingLedger` that it only witnesses accesses **through the proxy** — a consumer capturing
the raw fold or reaching `_fold` evades it (out of contract, but state it). Add a test for
`assert_delta_is_minimal`'s missing-required-row branch and `promotion_delta`'s forbidden-field/missing-evidence
guards.

## Report

Re-run `test_ledger_fold.py`; all tests (old + new) pass. The ≤25-line report states the changed functions, the
new tests, that each new test fails on its one-line ablation, and any clause not fully implementable (e.g. the
residual routed-key-contract cases where it raises rather than resolves). ⛔ No physics, no expected value.
