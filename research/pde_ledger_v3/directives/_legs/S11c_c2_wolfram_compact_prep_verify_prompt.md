# Verify — the S11c-c2 compact-prep state (docs accurate? next-prompt correct? clear to compact?)

The orchestrator is about to compact context. Verify, against the REAL repo + git history, that (1) the state
docs are accurate and not overstated, (2) the carry-forward set is right, and (3) the hand-back prompt is correct
and complete. ⛔ Do not rubber-stamp; if something is wrong or overstated, say so plainly with file:line / commit
evidence. This is a document + light-git verification; ⛔ do not modify the working tree.

## What just happened (claimed)
Two things closed this session, both committed on `ledger-v3-rebuild`:
- **c2 fold physics** was independently 2-leg reviewed (fresh Claude agent + Grok) → **SOUND, 0 confirmed defects**,
  committed `8f3a017f`. Grok flagged three "defects" (F uniform-limit, G adjointness, E N6) that the orchestrator
  adjudicated as FALSE POSITIVES by its own computation (rule 13).
- **c2 export repair** (publication-only, astra build) shrank `scripts/S11c_c2_exports.py` **60,516,900 →
  22,441,522 bytes**, directive gated `a5f7a06c` (2 decision legs), committed `aa76105a`, re-reviewed clear.

## 1. STATUS + adjudication accuracy (verify against git + the files)
Read the new top clause of `STATUS.md` (2026-09-06, "c2 STEP A + STEP B/C DONE") and the two adjudication records
`_measurements/S11c_c2_physics_review_adjudication.md` and `_measurements/S11c_c2_export_repair_rereview_adjudication.md`.
Check against `git log --oneline -8`, `git show --stat 8f3a017f aa76105a`, and the artifacts:
- Are commits `8f3a017f` / `a5f7a06c` / `aa76105a` / `d2befb7c` real and as described?
- Is `scripts/S11c_c2_exports.py` actually ~22.4 MB, does it import, is `s11cc2SelfEnergyIncrement` absent, and are
  both `s11cc2ClosedSlabOperator` + `s11cc2ClosedCouplingKernel` present (4-case trees)?
- Is the diff `8f3a017f`→working of `scripts/S11c_c2_selfenergy_fold_sympy_audit.py` genuinely publication-only
  (only `EXPORT_ROOTS`, the `export_key` map in `run()`, `publish`, and the new `publication_compact` helper — NO
  construction function changed)?
- Is "physics SOUND, 0 defects" appropriately supported, or overstated? In particular: are the F/G/E dispositions
  (F genuine coupling decouples / G directional self-energy / E leading-order rep-invariance holds, σ_W deferred)
  defensible from the adjudication's own evidence, or does any read as the orchestrator rationalizing away a real
  leg finding?

## 2. The carry-forward set (right + complete?)
The step record is claimed to owe: (F) "genuine coupling decouples", not "increment vanishes" (+ a light §5e/§3c
spec-wording fix); (E) σ_W-sector rep-invariance remnant deferred; (G) directional self-energy; the two S11c-b
sign conventions surfaced by the WL comparator. Against the physics adjudication + `directives/S11c_c2_SHARED_PHYSICS.md`
(§5e, §3c, §7): is this set correct and complete, or is a real carry-forward missing / miscategorized (e.g. is any
of these actually a build defect rather than a record note)?

## 3. The hand-back prompt (correct + complete + clear to compact?)
Read `/tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/next_prompt_c2_wolfram.md`.
Verify: the STATE line matches reality; the NEXT sequence is right (STEP 1 light spec-wording fix = review-until-clear;
STEP 2 blind WL engine with 2 decision legs BEFORE + 2 build legs SERIALIZED; STEP 3 T7 comparator + reconcile →
step record); the model policy + ops notes are correct; nothing outstanding is being skipped. Is anything in the
prompt wrong, missing, or an overstatement? Is there any reason we are NOT clear to compact (an uncommitted
load-bearing artifact, an open review leg, a claimed-committed thing that isn't)?

## Output
For each of 1–3: your finding + evidence (commit / file:line). End with: **CLEAR TO COMPACT** or the exact list to
fix first, and (if any) the precise wording corrections for the STATUS clause / next-prompt.
