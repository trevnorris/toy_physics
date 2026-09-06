# Gate record — CLAUDE.md CAS-instrument authorship clarification (2026-09-06)

**Artifact:** `CLAUDE.md` (governing process file) — E1, G4, at-a-glance item 2, + new `L-CAS` evidence-ledger
entry. **Role:** orchestrator-written governing prose → **review-until-clear** (G1/G2/G4), two legs = **Codex-sol
+ Grok** (orchestrator-written). My own reading is required (rule 13) but is **not** a leg.

**Why:** the 2026-09-05 streamline compressed 17 rules → M/E/G/S and raised the salience of "this binds the
orchestrator too" + "verify it myself" without the tacit "Codex codes" counterweight. Under that reading I
hand-authored the CAS scripts `verify_F.py`/`verify_EG.py`, ran them, and used them to overrule Grok's F/E/G —
**over-clearing E/N6** (rep-invariance holds only in the σ_W→0 projection; the retained-order residual is OPEN).
The Codex-sol compact-prep verify caught it. The fix makes explicit, for the first time in the file, that the
orchestrator **never authors (or privately runs) a CAS instrument** — Codex writes it, it is reviewed under G1,
and the orchestrator owns only the question + the adjudication + mechanical fact-lookups.

**Identical rendered prompts (committed):**
- round 1: `directives/_legs/CLAUDE_md_cas_authorship_review_prompt.md`
- round 2: `directives/_legs/CLAUDE_md_cas_authorship_rereview_prompt.md`

**Leg reports (raw transcripts kept OUTSIDE the repo as tree hygiene — build-skill §137; scratchpad):**
- R1 Grok `scratchpad/grok_claudemd_review.txt`; R1 Codex-sol `scratchpad/codex_claudemd_review_v2.txt`
- R2 Grok `scratchpad/grok_claudemd_rereview.txt`; R2 Codex-sol `scratchpad/codex_claudemd_rereview.txt`

## Round 1 — CONFIRMED findings (both legs; commit blocked)
1. **Grounding went conditional** (both) — "binds the orchestrator too — *but for mechanical fact-lookups only*"
   made command+stdout grounding apply only to mechanical claims. **Folded:** grounding stated unconditional
   ("CAS claims included" / "a CAS-backed claim is NOT exempt"); routing is a separate clause.
2. **Subjective boundary** (both) — "mechanical / no physics framing" is self-judged framing, the failure mode;
   a one-line `subs` is the relapse. **Folded:** operation-based definition (existence/retrieval/literal-count/
   shape only; any subs/simplify/solve/diff/limit/unit/tolerance/truncation ⇒ instrument, regardless of
   language/length/filename/library).
3. **"Grok + I review" ≠ G1** (both) + **author-vs-run conflation** (Grok) — orchestrator can't be a leg; and
   "I may *run* only lookups" left the reviewed instrument un-runnable. **Folded:** instrument is Codex-written →
   **G1 (fresh Claude + Grok)**; orchestrator neither authors nor privately-runs-and-clears; execution is the
   normal build/G1 path.
4. **Proxy-question root** (Grok #5) — the over-clear was the *question* (`subs(σ_W→0)==0 ⇒ "holds"`), a proxy.
   **Folded:** orchestrator agrees the question with Codex "at the retained order, not a convenient proxy."

**Adjudicated leg disagreement (rule 13):** Grok said the orchestrator *may execute* a frozen reviewed
instrument; Codex-sol said it *neither authors nor executes*. Resolved by making it moot — the instrument
follows **G1** (which handles execution + review); the orchestrator never authors it and never adjudicates from a
**privately-run** output. Both concerns satisfied.

## Round 2 — re-review of the fold (both legs)
- Defects 1–3 **FIXED** (both legs, quoted line-by-line); war-story **accurate** (static header/adjudication
  read); no control weakened; the orchestrator-only prohibition (first-person) does not sweep legs/builders;
  the G1 "O or Cx" build-script cell intact.
- **Two residuals (verbatim wording), folded:**
  - (Grok) the E1 **review** clause carried "not a proxy" but "at the retained order" sat only on the
    orchestrator+Codex side → added "asked at the retained order" to the **review** mandate.
  - (both) the **L-CAS** `⇒` still said "reviewed **blind-to-output**" (the legs must USE the output; only the
    orchestrator is blind-until-launch, G3) → "reviewed **under G1 (fresh Claude + Grok) before its output is
    trusted**."

## Stopping rule (G4) — CLEAR
Both round-2 legs marked everything else fixed/non-blocking; the only residuals were two verbatim swaps applying
already-reviewed phrasing to close the legs' own findings — **no new substance**. Nothing outstanding changes
what may be claimed ⇒ review-until-clear satisfied; no round-3 (that would be ceremony, not a substantive exit).

## Owed follow-ups (routed through review, not folded here)
- `build`/`review-legs` skill launch templates lack `< /dev/null` on `codex exec` (stdin-hang;
  [[feedback-codex-exec-hangs-on-stdin]]) and do not yet state the orchestrator-never-authors-CAS rule.
- `docs/development_pipeline.md` still says "fix the instrument" without the Codex-routing clause.
