# S11b comparator — fix round 2 directive (targeted; two verified defects + acceptance)

## Authority and boundary
Repair `research/pde_ledger_v3/scripts/S11b_cross_engine_comparator.py` IN PLACE (committed baseline
`bd97a571`). ⛔ Do not rewrite; fix round 1 correctly closed the three original holes and the precedence rule
`DISAGREE > UNCOMPARED > UNDECIDED > AGREE` is load-bearing and must stay. `CLAUDE.md` binds. Two defects
found by fix-round-1 review must be closed WITHOUT reopening any prior hole or opening a new one; the
comparator's one cardinal sin (a genuine difference reported as anything but DISAGREE — including buried as
UNCOMPARED) still governs.

## Defect 1 — the transliteration-collision guard is OVER-BROAD → false UNCOMPARED hides a real difference
`transliteration_collisions` / `_basic_source_labels` key injectivity by `KIND:name`, then group by
`mechanical_lower_camel(name)`. So a function head `fA` and a same-named symbol `fA` — which NEITHER rename
NOR ever collapse into each other (an applied function and a bare symbol are distinct objects after
transliteration) — are counted as two sources for target `fA` and flagged `TRANSLITERATION_COLLISION`. A
record `fA(x) + fA` vs `g(x) + h` is then reported UNCOMPARED, burying a genuine content difference.
Fix: a collision exists ONLY when transliteration merges two DISTINCT source objects OF THE SAME KIND into
one target — i.e. group by **(kind, target)** and flag a collision only when that group contains two or more
distinct SOURCE NAMES. `f_a` and `fA` in the same kind → collision (they collapse). A function `fA` and a
symbol `fA` → two different (kind,target) groups, one name each → ⛔ NOT a collision. A cross-kind pair of
DISTINCT names (`f_a` function + `fA` symbol) is likewise not a same-kind merge → not a collision; residual
it normally. ⛔ Do not flag a collision where no actual post-transliteration merge occurs.

## Defect 2 — classification is computed but not emitted until a huge operand renders → the run hangs
At the emit site, the whole `render_comparison` (which `sstr`-renders the full parsed operands
`PY_PARSED_OBJECT`/`WL_PARSED_OBJECT`) runs BEFORE the flush. For a large object (e.g. a ~120 KB payload) the
classification is computed in seconds but the operand `sstr` runs for minutes, so the already-decided
classification is never emitted and the stream appears to hang. The per-leaf residual budget does not
enclose rendering.
Fix: emit and FLUSH the object's NAME + CLASSIFICATION (+ reason) FIRST, before rendering the operands; then
render the operands under a size bound — if an operand's serialized length exceeds a fixed threshold, emit a
bounded/omitted form (e.g. a head + a `RENDER_TRUNCATED length=<n>` marker) instead of the full `sstr`.
⛔ The classification line must NEVER wait on operand rendering, and ⛔ the truncation may only affect the
DISPLAYED operand, never the classification.

## Acceptance strengthening (fix-round-1 review Finding A) — the committed battery must actually protect the fixes
The committed acceptance reads FROZEN `.stdout` snapshots and omits the directive-mandated rows, so reverting
a fix still passes. Fix: the acceptance must **re-run the comparator** on its fixtures (⛔ not read stale
snapshots) and assert the classifications, and its battery must include every extended row below plus the
fix-round-1 extended rows (function↔function and same-kind rename collisions; status precedence; symbol /
mixed / all-`Str` tuple keys; per-leaf budget). Regenerate any frozen snapshot from current code so no stale
`SOURCES=(...)` text remains.

## Acceptance — executable, value-free (rule 5); FROZEN before the real run
Keep the full existing battery passing, and ADD (each a fixture a defective fix FAILS, with an exact
per-name outcome), saving fixture + literal stdout to named build-scratch paths:
- **Same-name cross-kind, NOT a collision**: `fA(x) + fA` vs `g(x) + h` → DISAGREE (⛔ not UNCOMPARED); and
  `fA(x) + fA` vs `fA(x) + fA` → AGREE. Keep a genuine same-kind rename collision (`f_a` vs `fA`, same kind)
  → TRANSLITERATION_COLLISION/UNCOMPARED. Add a cross-kind DISTINCT-name pair (`f_a` function + `fA` symbol,
  genuinely different) → DISAGREE (not a false collision).
- **Classification before render**: a record whose operand serializes to a very large string AND whose
  classification is DISAGREE → the NAME+CLASSIFICATION line is emitted and flushed, and the operand appears
  only in bounded form; demonstrate (e.g. a timed run) that the classification is observable without waiting
  on the full operand render. A following object still classifies (the run continues).
- The strengthened acceptance, run against the CURRENT code, passes; and reverting either Defect-1 or
  Defect-2 fix (on a throwaway copy) makes it FAIL.
⛔ Do not run against the real `.out` pair as an acceptance gate and ⛔ do not tune anything to a real payload.

## Invariants
- ⛔ No behavioural change to the closed holes, the precedence rule, or the 7 core rules; the full existing
  synthetic acceptance and the repoint ablation must still pass.
- The three script clauses hold (PRINT; do not assert; interpretation → step record).
- ⛔ No fix may produce AGREE, UNCOMPARED, or UNDECIDED for a record that has a genuinely DISAGREEing leaf.

## Report (§13) — under 20 lines
Edits (file:line), each new/changed fixture + its stdout path, confirmation the full existing acceptance +
repoint still pass, confirmation the strengthened acceptance re-runs the tool and fails on a reverted fix,
and confirmation nothing was tuned to a real payload.
