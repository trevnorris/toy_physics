# Review record — "Committed ablation harnesses" standard (development_pipeline.md §4)

**Date:** 2026-09-07 · **Artifact reviewed:** the proposed §4 subsection (orchestrator-written process
proposal) · **Legs (G1, orchestrator-written → Codex + Grok):** codex-sol (gpt-5.6-sol, xhigh) + Grok
(grok-4.6, high), launched in parallel, identical prompt. Raw logs preserved in the session scratchpad
(`scratchpad/harness_review/{sol.log,grok.log}`); review prompt `scratchpad/harness_review/prompt.md`.

## Verdicts — both **SOUND-WITH-CHANGES** (convergent)
Both endorsed the core (S1-durable artifact; "both engines, never one" strengthens M1; the ephemeral
reviewer-chosen leg-ablation is retained, not replaced) and independently flagged the same defects in the
first draft. All findings verified by the orchestrator (G4) against CLAUDE.md + the pipeline doc + the N6
code before folding.

## Findings folded (each was a real defect in the draft)
1. **G2 waiver / target-list inversion.** Draft: "legs vet the manifest *in place of* a decision list" —
   contradicts §4 (target list is the orchestrator's, never the builder's) and skips G2. → Orchestrator writes
   the target/knife list, G2-gated before the harness builder launches; the harness only *mirrors* it.
2. **Self-confirming re-derivation.** "Reproduces the engine" invites a second derivation that agrees at
   baseline (copy hand-typed `H`, compute `H+δ` → byte-match + nonzero diff, proves nothing). → **Wrap the
   live production engine** (import/call/patch), never re-derive.
3. **Drift guard already false.** N6 SymPy emits RSS + timestamps, so byte-identical *stdout* fails every
   honest re-run. → Pin the guard to load-bearing **tagged payloads**, not full stdout; copy-identity, not the
   cut byte-identical-restore quarantine.
4. **Print, don't PASS; ablate the harness itself** (E1 + Phase 4). Emit `(baseline,corrupted,diff)` then
   guard; zero-diff under a claimed FORM knife = fail; a coefficient-rescale/dead-path mutation of the harness
   must fail to report a bite. Commit the run transcript to `_measurements/` (capability ≠ evidence).
5. **Scope + FORM/coefficient.** "Each load-bearing step" → the **claim-dependency frontier** (N6 was 5 knives,
   not 40 tags). FORM default (changes structure preserving type/dim/retained-grade); a sign-flip/rescale alone
   is coefficient; a **named coefficient exception** is allowed when a FORM knife can't see a channel (N6
   Φ-coefficient precedent).
6. **Backfill trigger = dependency, not disagreement.** Designed-to-agree / common-mode / shared freezes
   produce *agreement*; a disagreement-only trigger never fires on them. → Retrofit a directly-consumed prior
   premise that has no committed bite-triple, risk-ranked.

## Concrete downstream finding (both legs, independent)
The **S11c-b `slab_operator` pressure-slot carrier** — `∂(slab rows)/∂(δp±, ∂_w δp±)|_{P=0}` plus the #90
`closure_shape_deriv` fold — is the highest-value early retrofit: c2 binds it as a **supplied/unfalsifiable**
premise and it *becomes* c2's `C_E/C_M`, yet its #90 sign and cross-engine residual remain deferred. Decided
(user 2026-09-07): harden S11c-b before continuing. c1's response coefficients = second priority (already
cross-engine agreed).

## Disposition
Folded into `docs/development_pipeline.md` §4 (+ Roles table Codex row + §10 item 6). User signed off →
adopted. Backfill scope retargeted (user + orchestrator): not blanket S9→S11, but the **surviving
load-bearing chain**, risk-ranked, gap-audited at S11c-e closeout (dead/superseded branches excluded per
"ledger = surviving solution only").
