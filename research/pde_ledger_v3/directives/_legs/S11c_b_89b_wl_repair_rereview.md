# Independent re-review — S11c-b #89b WL operator REPAIR (Codex build)

You are one of two independent re-review legs (a fresh Claude agent + Grok) on a **Codex-written** repair of the
WL brane-operator engine. Whatever writes does not review. Derive everything yourself and demand computation —
a prose re-derivation is worth nothing; every claim carries the script and its literal stdout. Report a finding
only if it catches a way the physics could be wrong. ⛔ Do NOT trust the builder's own self-check numbers — rebuild
each check yourself.

## Artifact
`research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl` (working tree, uncommitted, just
repaired). Inspect the repair with `git diff HEAD -- <that file>` — but only ~28 net lines are the repair; the
bulk of the diff vs HEAD is the earlier #89b un-freeze + gate, already reviewed.

## Context — what the repair was supposed to fix (from the build review)
Three defects were found and repaired (record: `directives/_measurements/S11c_b_89b_wl_operator_build_review.md`;
directive: `directives/S11c_b_89b_wl_operator_repair.md`). Your job is to INDEPENDENTLY confirm each is genuinely
fixed AND that nothing correct regressed. ⚠ Full-depth is ~16 GB/case and OOMs the 30 GB box — work at REDUCED
scale (`basisRepresentativeIndices = {16}` or a small jet-bearing subset), extract definitions into a `/tmp`
harness, `timeout 600` per kernel, ONE kernel at a time, ablate `/tmp` copies only, kill leftover `WolframKernel`
by EXACT pid. Save every script + literal stdout to named absolute paths and report them.

## What you must independently establish
1. **Defect A — the emitted operator re-freeze is GONE (the decisive one).** The emitted `OPERATOR`/`ORIGINS` must
   now be the activate-THEN-reduce object (outer constrained `Inactive[Div]` activated BEFORE the final background
   reduction). Build your OWN reference by reducing `activateSpatialDivergences[operatorLive]` (from the
   un-reduced `KERNEL_SOURCE_OPERATOR`, ⛔ NOT by activating the emitted operator — reduce-then-activate is the
   frozen object and would falsely agree). On the U row: (a) the emitted operator equals your from-`operatorLive`
   reference; (b) the emitted U row now carries the mixed/higher background jets (e.g. the coefficient of
   `widthProfileJet[1,1,0]` is non-trivial, not a zero vector); (c) DECISIVENESS — a one-sided swap of one route
   to reduce-then-activate MOVES the residual off zero (⛔ if it does not, the check is not distinguishing the
   fix). Report the literal residuals and coefficients.
2. **Defect A ripples — comparisons are now like-with-like.** (a) The frozen Hessian-witness (`surfaceDifferenceJetAtoms`,
   deferred control) must compare derivative-normal-vs-derivative-normal for EVERY slot it differences
   (`OPERATOR`, `ORIGINS`, `MU_THETA`, `DIVERGENCE_FORM_SOURCE`, faces) — check both sides carry 0 leftover
   `Inactive[Div]` so the atom difference is not a Div-vs-expanded artifact. (b) The uniform-limit S11b residual
   (`SLAB_OPERATOR` AND `TRANSVERSE_DISPERSION`) must likewise compare like forms. Confirm the S11b reconstruction
   physics and the frozen replay/kernel were NOT changed to achieve this (only comparison-only surfaces added).
3. **Defect B — §5.E dimension walker now attaches.** Confirm `primitiveExpressionDimension` now walks composite
   `Times`/`Plus` invariants to their primitives. Decisive test: mutate ONE primitive atom's dimension and show
   the emitted residual MOVES for the invariants containing that atom and stays put for those that do not (⛔ not
   a blanket move — it must be selective). Confirm no invariant is special-cased or hard-coded.
4. **Defect C — independence control is no longer a `base − base = 0` tautology.** For a MATERIAL_ADVECTED branch,
   confirm the whole independence package (SLAB, ADMISSIBILITY, and COUPLING_KERNEL slots together) is either a
   genuine nonzero-capable corruption test or carries an explicit `VALIDATED_ON_REPRESENTATIVE_BRANCH` marker —
   ⛔ never a silent identical base/corrupted on some slots with a stale kernel on another. The `LAB_HELD` branch's
   live one-sided corruption must remain live (a real defect there still drives the residual nonzero).
5. **NO REGRESSION (ablate to confirm).** The repair must not have touched the physics that was correct:
   `KERNEL_SOURCE_OPERATOR`/`KERNEL_SOURCE_ORIGINS` remain the un-reduced LIVE operator/origins (still activatable
   to derivative-normal; kernel higher-jet depth unchanged); the `LIVE_DIVERGENCE_FORM_OPERAND` tower slot stays
   un-reduced-live; the `S11CB_SKIP_HEAVY_CONTROLS` deferral gate still toggles marker-vs-real correctly; the
   inner live-energy-EL un-freeze and the tractable activation (`operatorLive`-based) are unchanged. Name what you
   ablated to confirm each.

## Method + hygiene
- Derive independently; where you dispute a fix, write your own script, run it, save script + literal stdout to a
  named absolute path you report. A prose re-derivation is discarded.
- ⛔ Mathematica: wrap EVERY kernel run in `timeout 600`; ONE kernel at a time (2-seat licence); ablate `/tmp`
  copies only, never the working tree; kill any leftover `WolframKernel` by EXACT pid (never `pkill -f`); ⛔ do
  NOT run the full engine — extract a minimal reduced-scale harness.
- Report a numbered list: `severity — file:line — finding — the physics it gets wrong — the computation (script
  path + literal stdout) that shows it`. Severity ∈ {BLOCKER, SHOULD-FIX, NIT}. Then `VERDICT: N findings (B
  blockers)` or `VERDICT: CLEAR`. If you clear it, name the two or three ablations you actively ran to try to
  break each fix — a leg that finds nothing is weak evidence.
