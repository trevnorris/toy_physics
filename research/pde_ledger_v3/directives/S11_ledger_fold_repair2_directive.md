# Repair directive #2 — ledger_fold guard: resolve edges by IDENTITY, verified on the real ledger

Repair `research/pde_ledger_v3/scripts/ledger_fold.py` and `research/pde_ledger_v3/scripts/test_ledger_fold.py`.
The previous repair over-corrected: **verified against the real frozen base**, it false-positives — it demands
the structural window function `O_window` (0 rows, referenced 3,284× as `Function('O_window')`) be a row, and
it treats symbol **users** (`f_hold_e_W_0`, class `PREMISE`, which merely mentions `Symbol('W_0')`) as
**producers** of `W_0`. Authority: the **refined design §D3** (committed `fd8c89d0`) and the F9c-severity scan
`_measurements/f9c_pair_scan.py` (16 routed pairs; 11 identity-resolvable; 5 genuine same-srepr; **0 of the 5 on
the critical path** — only `omega`, resolvable, is). ⛔ Physics-free infra.

## Fix A — edge resolution by FULL SYMBOL IDENTITY, with structural-skip (§D3 step 2)

Resolve each referenced atom (from `atoms(Symbol)` ∪ `atoms(AppliedUndef)` function-heads ∪ `atoms(MatrixSymbol)`
of a `Basic` value) against the fold by the atom's **full identity** — its `srepr` (name **plus assumptions**),
⛔ **not** its bare name. Three outcomes:
- **exactly one** fold row has a `value` with that identity ⇒ **resolved**, recurse into it. (So
  `Symbol('omega', real=True)` resolves to the bare `omega` row, ⛔ not to `s11b_omega=Symbol('omega')`.)
- **zero** fold rows ⇒ the atom is a **free/structural** symbol or function (a coordinate; the window
  `O_window`) ⇒ it is **not a required dependency** — ⛔ **skip it, do not raise**. (This reverts the O_window
  false-positive.)
- **two or more** fold rows share that identity under **different keys** ⇒ **genuine ambiguity** ⇒
  raise-and-name (`AmbiguousSymbolError`). Measured: 0 such on the critical path.

⛔ **The identity→producer index keys on the row's OWN identity, not on symbols the row USES.** A row is a
candidate producer of atom `a` iff its `value` **is** `a` (a declaration) — ⛔ **not** iff `a` appears in its
value. (This reverts the `f_hold_e_W_0`-is-a-producer-of-`W_0` false-positive.) `dimension_key` closure is
unchanged (recurse into the exact target key).

## Fix B — the under-export guarantee is the SMOKE-TEST, not "assert every referenced name is a row"

Drop any requirement that a **bare/absent referenced name must be a row** — that is what false-positived, and it
is wrong: a symbol used purely symbolically (a knob as a free parameter, a coordinate, a structural function)
needs **no** row. The genuine under-export catch is the **bidirectional smoke-test** (`check_consumer` step 3,
unchanged): a consumer's build that **looks up** an absent required row fails there. The static closure's job is
only (i) expand `IMPORT_KEYS` to their identity-resolvable transitive deps for manifest completeness and (ii)
raise on genuine ambiguity — ⛔ it does not police free symbols.

## Fix C — the tests, reframed to the identity semantics, PLUS a mandatory real-ledger test

Keep the decisive tests that are still correct (last-wins fold; no-F9-reapply; exact-key `predecessor_trap`;
bidirectional smoke-test; minimum-mode; promotion-delta; the depth-2 recursion chain). **Replace** the prior
repair's now-wrong tests (a function/absent-name edge that "must raise"; the name-based expression-F9c raise)
with:
1. **identity-resolution**: a fold with bare `omega=Symbol('omega',real=True)` and `s11b_omega=Symbol('omega')`;
   a kept row referencing `Symbol('omega',real=True)` ⇒ **resolves** (no raise) to the bare row; referencing
   `Symbol('omega')` ⇒ resolves to `s11b_omega`.
2. **structural-skip**: a kept row referencing `Function('O_window')(x)` (and a `MatrixSymbol`) with **no** such
   row ⇒ **no raise** (skipped as structural); closure excludes it.
3. **genuine ambiguity**: two rows with the **same** identity `Symbol('z')` under keys `z` and `s11_z` ⇒ a kept
   row referencing `Symbol('z')` ⇒ **raise-and-name** both.
4. **user-is-not-a-producer**: a row `uses_w` whose value is `Symbol('W_0')+1` alongside a `W_0` declaration ⇒ a
   kept row referencing `Symbol('W_0')` ⇒ resolves to `W_0` only, **no** ambiguity from `uses_w`.
5. ⭐ **MANDATORY real-ledger test** (`test_ledger_fold.py`, guarded by presence of the file): import
   `scripts/S11c_b_exports.py`'s `LEDGER` and run `check_consumer` on each real critical-path root —
   `mu_theta_operator`, `slab_operator`, `coupling_kernel`, `face_normal`, `conormal_deriv`,
   `face_measure_shape_deriv`, `face_velocity`, `relative_flux`, `kinematic_balance`, `traction`, `face_shift`,
   `closure_shape_deriv` — and assert **each resolves without raising** (the `omega` edge resolves by identity;
   `O_window` is skipped). This is the test whose absence let the false-positives ship; it must pass.

Each test prints operands + observed + literal pass/fail (rule 2); the harness exits nonzero on a real failure.
Confirm each new test fails on the corresponding one-line ablation (identity→name; skip→raise; declaration→uses).

## Report

Re-run `test_ledger_fold.py`; all pass, including the real-ledger root test. The ≤25-line report states the
changed functions, the reverted false-positives, the new/removed tests, and that `check_consumer` resolves every
listed real critical-path root. ⛔ No physics, no expected value.
