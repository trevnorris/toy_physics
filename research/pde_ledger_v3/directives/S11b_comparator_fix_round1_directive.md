# S11b comparator — fix round 1 directive (folded once after two directive legs)

## Authority and boundary
Repair `research/pde_ledger_v3/scripts/S11b_cross_engine_comparator.py` IN PLACE (committed baseline). ⛔ Do
not rewrite; the 7 core delta rules are ablation-confirmed correct and stay behaviourally unchanged except as
sharpened below. `CLAUDE.md` binds. The comparator's ONE cardinal sin: ⛔ reporting a genuine
content/sign/dimension/structure difference as anything other than DISAGREE. All fixes below close a live
instance of that WITHOUT opening a new one.

## The unifying rule — classify by LEAF, then take the most severe (adopt this; it resolves most defects)
Classify every LEAF of a record, then the record's class is the **highest-severity leaf** under the strict
precedence **`DISAGREE > UNCOMPARED > UNDECIDED > AGREE`**. Consequences, all mandatory:
- Residual the FULL payload leaf-by-leaf — ⛔ never strip status/coverage keys before residualing, and ⛔
  never short-circuit the whole record on the presence of one status key.
- A differing sibling (DISAGREE) beside a status token, a boolean leaf, or a timed-out leaf ⇒ the record is
  DISAGREE.
- UNDECIDED is the record's class ONLY when its severest leaf is UNDECIDED — i.e. a status/coverage token
  with no disagreeing, uncomparable, or budget-exceeded sibling (a status-token-ONLY payload, or one whose
  every other leaf AGREEs).
- ⛔ No bucket (UNCOMPARED, UNDECIDED, BOOLEAN_NOT_RESIDUALABLE) may sit ABOVE DISAGREE in precedence.

## Defects to fix (each measured against the baseline)
1. **Function-head transliteration collision → false AGREE.** The injectivity guard
   (`transliteration_collisions`) inspects only `atoms(sp.Symbol)`, but `_transliterate_basic` also renames
   `AppliedUndef` function heads containing `_`. Fix: build ONE injective target map over BOTH
   `atoms(sp.Symbol)` AND `atoms(AppliedUndef)` heads, keyed by `mechanical_lower_camel(name)`; ANY target
   reached by two distinct sources — symbol↔symbol, function↔function, OR function↔symbol — flags the object
   `TRANSLITERATION_COLLISION` (an UNCOMPARED-class reason), ⛔ never a silent collapse.
2. **Status/UNDECIDED short-circuit buries a differing sibling.** `compare_records` currently returns
   UNDECIDED as soon as a `STATUS_TOKEN`/`COVERAGE` key appears, before residualing siblings. Fix: apply the
   unifying leaf-classification rule above — residual the full Association (status/coverage keys residual as
   ordinary leaves; two equal uncertain tokens are AGREE at that leaf, two DIFFERENT tokens are DISAGREE at
   that leaf, a lone uncertain token with no counterpart is UNDECIDED at that leaf), then take the most
   severe leaf. ⛔ A status key may never suppress a real difference in ANY sibling, at any nesting depth.
3. **Tuple-of-textual-pairs promotion must require EVERY key textual, never a bare Symbol.**
   `_convert_parsed_containers` promotes a 2-tuple when `_association_key(item[0])` is non-None, and
   `_text_value` returns a bare `Symbol`'s name — so `((p,1),(q,2))` and even a MIXED `((Str("p"),1),(q,2))`
   promote and match a WL Association, hiding the tuple↔Association STRUCTURE DISAGREE. Fix: promote ONLY when
   EVERY pair-key is an explicit textual atom of the SymPy Association encoding (`Str`/`_Str`/`TextAtom`/`str`)
   — ⛔ if ANY key is a bare `Symbol`, do not promote; the tuple stays a tuple and tuple↔Association is a
   STRUCTURE DISAGREE. A genuinely `Str`-keyed tuple still promotes and residuals normally.

## Runtime robustness — ⛔ neither item may ever produce AGREE or bury a DISAGREE
4. **Stream: compute → render → flush PER OBJECT.** ⛔ Do not compute the whole join then flush — that is 8
   minutes of silence. Each object's classification line is computed, printed, and flushed before the next
   object begins, so progress is observable while a later object is still computing.
5. **Per-LEAF residual budget, not per-object.** Bound the algebraic residual of a single LEAF
   (`factor`/`cancel`/`together`); on exceeding the budget, that LEAF is `ResidualFailure(RESIDUAL_BUDGET_EXCEEDED)`
   (an UNCOMPARED-class leaf). The record's class then follows the precedence rule — so a differing cheap
   sibling still makes the record DISAGREE. ⛔ Never bound the whole record's residual as one unit (that would
   bury a cheap differing sibling behind a heavy one).

## Acceptance — extend the battery so a defective fix FAILS (value-free, rule 5)
Keep the existing decisive fixtures (they must still pass) and ADD, each with an EXACT per-name outcome:
- **Collision, cross-kind, PARSED both sides**: PY `f_a(x) + fA(x)` vs a PARSED WL analogue with two heads
  colliding (`fA[x] + fA[x]` grammar) → the object's reason contains `TRANSLITERATION_COLLISION` and it is
  UNCOMPARED, ⛔ NOT AGREE. Include BOTH a function↔function and a function↔symbol collision row. Keep the
  existing symbol↔symbol collision fixture.
- **Status precedence**: `{STATUS_TOKEN:<uncertain>, RESULT:a}` vs `{…,RESULT:b}` → DISAGREE;
  `{STATUS_TOKEN:<uncertain-1>, RESULT:x}` vs `{STATUS_TOKEN:<uncertain-2>, RESULT:x}` (different tokens,
  equal RESULT) → DISAGREE; a NESTED status+differing-sibling → DISAGREE; a `COVERAGE_*` key case → same;
  and a status-token-ONLY payload vs a resolvable one → UNDECIDED.
- **Tuple keys**: PY `((p,1),(q,2))` (Symbol keys) vs WL Association → STRUCTURE DISAGREE; MIXED
  `((Str("p"),1),(q,2))` vs Association → STRUCTURE DISAGREE; ALL-`Str` `((Str("p"),1),(Str("q"),2))` vs the
  same Association → AGREE (legitimate promotion preserved).
- **Per-leaf budget**: a record with a budget-blowing heavy leaf AND a cheap DIFFERING sibling → DISAGREE
  (⛔ not UNCOMPARED, ⛔ not AGREE); and a following cheap DISAGREE row still classifies (the run continues).
- **Streaming**: demonstrate the first object's line is flushed while a later heavy object is still computing
  (⛔ not all-at-the-end).
Save every fixture + its literal stdout to named build-scratch paths, reported in §report. ⛔ Do not run
against the real `.out` pair as an acceptance gate and ⛔ do not tune anything to a real payload.

## Invariants
- ⛔ No behavioural change to the 7 core rules except as items 1-3 + the precedence rule sharpen them; the
  full existing synthetic acceptance must still pass.
- The three script clauses hold (PRINT residuals/classifications; do not assert; interpretation → step record).
- ⛔ No fix may create a path that reports AGREE, UNCOMPARED, or UNDECIDED for a record that has a genuinely
  DISAGREEing leaf.

## Report (§13) — under 20 lines
Edits (file:line), each new fixture + its stdout path, confirmation the full existing acceptance still passes,
and confirmation nothing was tuned to a real payload.
