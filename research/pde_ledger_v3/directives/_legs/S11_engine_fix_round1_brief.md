# S11 SymPy engine — fix round 1. ⭐ Four items. ⛔ ONE of them is physics.

## Authority and boundary

Edit `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py` **only**. ⛔ Change no other
file. `CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` is the sole physics
authority. `research/pde_ledger_v3/directives/S11_sympy_build_directive.md` is **closed** — ⛔ do not
re-open it, ⛔ do not change what is computed.

⭐ **The engine is not broken in its physics content.** It compiles, emits genuine `srepr` CAS objects, and
runs 21 of 22 declared cells. ⛔ Four runs have failed to publish `S11_exports.py`.

Evidence for every claim below: `directives/_measurements/S11_engine_fix_round1_brief.md`.

## ⛔⛔ 1 · THE ZERO TEST — ⭐ this is the blocker, and it is a CORRECTNESS defect

An isolated reproduction of `run_cell('MAIN', 5)` (`/home/trevnorris/.s11_build/repro_d5.py`) fails with
`MemoryError` at an exact frame:

```
run_cell:1380 → emit_q4:1239 → matrix_rank:870 → sympy rank() → _row_reduce → cross_cancel
```

```
$ sed -n '868,870p' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
def matrix_rank(matrix: sp.MatrixBase) -> int:
    simplified = sp.Matrix(matrix).applyfunc(lambda entry: sp.factor(sp.cancel(entry)))
    return int(simplified.rank(iszerofunc=lambda entry: entry == 0, simplify=False))
```

⚠ The **input** entries are pre-simplified. ⛔ The **intermediate** entries `cross_cancel` builds during
elimination are not (`simplify=False`), and those are what `iszerofunc` receives. `entry == 0` is SymPy
**structural** equality, so an intermediate that is zero only after reduction reads non-zero.

⛔⛔ **REVISED AFTER TWO LEGS — the paragraph above is WHERE IT SURFACES, ⛔ NOT WHY.** ⚠ Both legs measured
it independently: the matrix reaching `rank` is **small** (5×5, entries ≤22 ops), there were **zero** cases
of an algebraically-zero entry read as non-zero, and **a stronger zero predicate still `MemoryError`s**.
⇒ ⛔ **The structural zero test is NOT the cause.** ⭐ The cause is **expression-tree Gaussian elimination
swelling from small exact entries** during row reduction.
⚠ ⭐ The structural test **remains a real correctness hazard** — it never returns `None`, so SymPy never
falls back to simplification — ⛔ but it is a separate item, not the blocker.

⇒ ⭐⭐ **THE OBJECT IS Q4's `rank` AND `stacked_rank` OF THE LIVE `M_r`** (`:1239`, `:1242`) — ⭐ the exact
algebraic rank, over the coefficient field the cell already carries. ⛔ Not a zero-test property, ⛔ not a
generic rank at a numeric specialization, ⛔ not any cheaper object that happens to return an integer.
⛔ **Do not specify the construction** (rule 3). ⭐ **A terminating construction exists** — ⚠ both legs
found one — ⇒ ⛔ this is **not** a §10 unavailability, and ⛔ must not be reported as one.

⛔ **`:918-920` (`nullspace_basis`) is DEAD — it has no caller.** ⚠ Scoping it was the list's error. ⭐ The
live nullspace path is `generic_nullspace_vectors` (`:924`, called at `:1245`).

⛔⛔ **THE NEXT WALL, and the list missed it: `all_minors` (`:906-915`, called `:1289`).** ⚠ Measured at
D=5: ~26 s per 4×4 determinant × 25 minors, timing out at 120 s. ⇒ ⭐ **fixing rank alone does NOT complete
the cell.** ⭐ Same instruction: name the object Q8 requires; ⛔ do not invent the routine.

⛔⛔ **A CHANGED RESULT AT D2/D3/D4 IS A FINDING TO REPORT UNDER §10, ⛔ NOT A REGRESSION TO SUPPRESS.**
⚠ Those dimensions completed under the same zero test. ⛔ Do not tune anything to reproduce prior output,
and ⛔ do not treat prior output as a reference.

## ⛔ 2 · A CELL FAILURE MUST BE VISIBLE WHEN IT HAPPENS

`:1668` catches a failing `run_cell` into `ISSUES`, and `ISSUES` is only rendered in the §10 report **after
the whole package loop**. No run has ever reached it. ⇒ ⛔ four runs and ~14 hours produced **no
explanation of their own failure**, while the explanation existed in an unreachable place.

⇒ ⭐ A cell that fails must say so **at the time**, on the stream, with its exception type — without
altering §10's existing role.
⚠ There are bare `except Exception` handlers at `:756 :775 :955 :1317 :1450 :1668`. ⛔ **They are NOT
equal, and the list treated them as if they were.**
⛔⛔ **`:955` and `:1450` fall back with NO `ISSUES` entry at all** — ⚠ `:955` `inv` → `gauss_jordan_solve`;
`:1450` a failed comparison → `UNDECIDED`, which `F9` equality then consumes. ⇒ ⭐ **a silent physics
claim.** ⛔ `:1668` swallows the `MemoryError` and records it only where nothing ever reads it.
⚠⚠ **The two legs DISAGREED on `:756`** — one called it a wrong-payload risk (a forced Gröbner failure
returning `NOT_APPLICABLE`), the other called it soft because it does append to `ISSUES`. ⭐ **That
disagreement is a finding: settle it by computation and report which is right under §10.**
⛔ Do not silently delete any handler — ⚠ several encode genuine §10 unavailable constructions, and
deleting `:1668` without a stream emit reintroduces hard aborts that kill every later cell.

## ⛔ 3 · A COMPLETED `MAIN` MUST SURVIVE AN INTERRUPTION

```
$ sed -n '1680,1682p' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
    if main_completed == main_declared:
        ledger = merged_export(main_dim_data, run_pairs_payload, skipped_pairs_payload)
        write_exports(ledger)
```

⚠ This runs after the **whole** package loop, though its condition covers `MAIN` only. ⇒ a finished primary
package is held by later control cells, and any interruption discards it.

⛔⛔ **THIS IS NOT A FREE MOVE, and the brief said otherwise until this line was added.**
`run_pairs_payload` and `skipped_pairs_payload` are built **after** the loop from `completed_pairs`, and
`merged_export` consumes both. ⇒ ⭐ writing at MAIN-completion means those payloads are **not yet final**.

⇒ ⭐ **State what must hold: a completed `MAIN` survives an interruption, AND no row claims a run record it
does not have.** ⛔ Do not alter `F6`'s publish semantics, its MAIN-only condition, or the stale-export
deletion. ⚠ If those two cannot both hold, ⛔ do not guess — ⭐ report the conflict under §10 and leave the
placement as it is.

## ⭐ 4 · THE SCRIPT WRITES ITS OWN STDOUT

The script prints to stdout and has no `out/` writing; every run so far landed in `/tmp`, and the committed
`out/` file is a truncated casualty of an earlier failure. ⇒ ⭐ the run's tag stream must also reach
`research/pde_ledger_v3/scripts/out/S11_stray_longitudinal_sympy_audit.out`. ⛔ Do not change the tag
stream's content or ordering.

## Boundaries

- ⛔⛔ **No memory cap, no timeout, and no handler that swallows the swell.** Those force the failure rather
  than remove it, and were explicitly rejected by the user.
- ⛔ Do not change what is computed, the declared `PACKAGE_DIMS`, the `D` range, or any package.
- ⛔ No expected value and no acceptance criterion referencing one (`CLAUDE.md` rule 5).
- ⛔ Do not add a registry, completion map, or exit policy — the closed directive forbids parallel machinery.

## Acceptance

⛔⛔ **REVISED — the list's first acceptance was UNREACHABLE by its own item 1**, because `all_minors` walls
after rank. ⚠ That is the same defect the orchestrator was caught on earlier the same day: demanding an
end-state without measuring the artifact can deliver it.

⭐ **`/home/trevnorris/.s11_build/repro_d5.py` runs `run_cell('MAIN', 5)` to completion** — ⛔ not merely
past `:870` — ⭐ and the cell emits its terminal tags, so that in a full run it would enter
`completed_pairs`. ⭐ Report the literal stdout.
⚠ Use it before any full run: ~minutes against ~6 hours, and a full run has OOM-killed this machine.
⛔ No expected value is stated here and none may be introduced (rule 5).

## Deliverable

The edited script, plus a note saying for each item what changed, and — for item 1 — the property you gave
the zero test and whether any D2/D3/D4 value moved.
