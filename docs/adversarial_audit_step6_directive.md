# Adversarial Audit — Step 6 Close-out Directive (Codex applies, Claude reviews)

**Status:** ACTIVE. Driver = `docs/adversarial_audit_execution_plan.md` Step 6. Phase C is COMPLETE
(no surviving fatal flaw across the entire layer-2 audit). This directive closes out the audit:
make the audit's coverage visible in the archive paper, mark layer-2 complete in the falsification
stack, and apply the non-fatal close-out fixes catalogued in the canonical punch list.

**Owner split (binding).** Codex APPLIES every content-bearing change in this directive (generator
scripts, generated `.tex` appendices, paper-build wiring, card/notes/script/provenance-YAML edits).
Claude reviews. Claude does not hand-apply. After each part: **leave the tree dirty — do NOT commit**
(Claude reviews, then Claude commits with explicit paths).

## Binding constraints (carried from the campaign contract)

- **Frozen:** `docs/adversarial_audit_directive.md` must never be edited.
- **Script timeout:** every script/build runs under `timeout 600`; exit 124 = failure → reformulate,
  never raise the cap.
- **Builds:** `latexmk` is unavailable. Build with `pdflatex <master>` run **twice** per the README
  fallback (`pdflatex pde_ledger.tex; pdflatex pde_ledger.tex` and likewise the reader build), from
  `research/pde_ledger/paper/`, under `timeout 600`. Iterate until exit 0.
- **No `git add -A`** and **do not commit** — Claude commits explicit paths after review.
- **Framing:** strict toy-analog framing in all paper-facing prose. NEVER write "validated / proven /
  evidence the model is true." The audit result is "no surviving fatal flaw," an operational-disclosure
  statement, not a truth claim.
- **No hand-authored content presented as generated.** The two new appendices MUST be emitted by a
  committed, deterministic, re-runnable generator that reads the audit artifacts. Re-running it must
  reproduce byte-identical output (modulo an embedded build timestamp, if any, which must be passed in,
  not read from the clock).
- **No silent caps.** If the generators exclude or group anything (aliases, scanner artifacts, empty
  bundles), the appendix must show the exclusion (footnote / "N aliases folded" line), never drop it.

## OFF LIMITS this step (user-gated — do NOT touch)

Two punch-list items are CONTINGENT on the user's `ANSATZ_LEDGER.md` §5 layer-3 call and are deferred:
- **A8** — `paper/stages/stage_162.tex` γ₀ `\StatusExactClosure` badge / γ₀=(1+𝔯²)/9 wording.
- **A9** — `paper/stages/stage_240.tex:15` "Exact selected loading ratio" c0 qualifier.

Do not edit stage_162's γ₀ badge or stage_240:15. They will be resolved (or confirmed no-edit) after
the user decides the §5 γ₀/c0 readings. Everything else in the punch list is in scope (Part B).

---

# PART A — Generated archive appendices + falsification-stack update

Goal: make the layer-2 audit's coverage and the parameter genealogy visible in the **archive** build,
as appendices GENERATED from the real `redteam_adversarial/` artifacts, and update the falsification
stack to show layer-2 complete. Design the format and the generator route yourself; the requirements
and acceptance criteria below are binding, the implementation is yours.

## A-1. Fit-insertion index appendix (generated)

- New file `paper/appendices/fit_insertion_index.tex`, emitted by a new committed generator.
- **Sources:** `redteam_adversarial/fit_insertion_points.yaml` (the 922 anchored candidates) +
  `redteam_adversarial/MANIFEST.yaml` (per-candidate `status` + any `verdict`) +
  `redteam_adversarial/provenance/*.yaml` (for each candidate's `constraint_kind`).
- **Content (a coverage disclosure — "what the adversarial audit examined and what it found"):**
  every Phase-A candidate, organized navigably (group by stage ascending). Each row carries at least:
  stage · parameter / candidate · source anchor (`file:line`) · `constraint_kind`
  (internal_consistency / free_choice / published_target) · audit disposition (the MANIFEST status:
  `audited` / `verdict_logged` + the verdict, `scanned` for the handful never advanced, etc.).
- **Exhaustive.** Every candidate appears, or is accounted for by a visible grouping line (aliases,
  scanner artifacts). A reader must be able to see that the audit's coverage is the full candidate set.
- Add a short generated header paragraph stating the totals (N candidates, the constraint_kind split,
  the audited/verdict_logged counts) and that the audit found no surviving fatal flaw — toy-analog
  framing, no truth claim.

## A-2. Parameter-provenance appendix (generated)

- New file `paper/appendices/parameter_provenance_audit.tex`, emitted by a (the same or a second)
  committed generator.
- **Sources:** `redteam_adversarial/provenance/*.yaml` (genealogy: origin stage, `constraint_kind`,
  derivation-vs-posit basis, `downstream_dependents`) + `redteam_adversarial/benchmarks.yaml` (the
  cited external matches for `published_target` entries).
- **Content:** per canonical parameter / family, its provenance slice: origin stage · `constraint_kind`
  · basis (the located derivation step for internal_consistency; "posited choice" for free_choice; the
  cited literature benchmark for published_target) · downstream dependents. Group by family or by stage
  — your call — but make it navigable and exhaustive over the canonicals (aliases shown as aliases, not
  dropped).
- For the `published_target` rows, surface the actual `benchmarks.yaml` citation (e.g. Peters 1964 for
  the GR quadrupole; CODATA 2018 for m_p/m_e) so the external-fit surface is explicitly disclosed.

## A-3. Wire into the ARCHIVE build only

- Add both appendices to `paper/pde_ledger.tex` after the existing audit-adjacent appendices (e.g.
  after `\input{appendices/source_file_index}`), in a consistent order, before `\input{appendices/fill_workflow}`
  or wherever reads best. **Archive build only** — do NOT add them to `pde_ledger_reader.tex`.

## A-4. Falsification-stack table → layer-2 complete

- `paper/appendices/reader_verification_summary.tex` (the "Falsification stack" `longtable`, ~lines
  60–84). Update the **"Adversarial fit-vs-derive audit"** row from "Designed in …; not yet run." to a
  COMPLETE status that accurately reflects the result:
  - layer-2 fit-vs-derive audit COMPLETE; no surviving fatal flaw;
  - free_choice + published_target + HIGH + internal-consistency completeness-critic all cleared;
  - 2 non-fatal `verdict_logged` PARTIALs (stage_192 card re-attribution; γ₀ badge scoping), both
    contingent on existing user-gated items.
  - Cross-reference the two new archive appendices and/or the `redteam_adversarial/` record so the
    disclosure is navigable.
- Update the table's "Status as of <date>" header to **2026-06-12**.
- Toy-analog framing only; no "validated."
- Check whether the **archive** build discloses the falsification stack anywhere (it currently uses the
  reader-only `reader_verification_summary`). If a parallel archive disclosure exists, update it
  consistently. If not, the new archive audit appendices ARE the archive's disclosure — leave it.

## A-5. Audit-status vocabulary (`fit-vs-derive-audited`) — REPORT, do not hand-tag

The `fit-vs-derive-audited` audit-status term is defined in
`paper/frontmatter/02_claim_status_firewall.tex` and is stated to be machine-generated, never
hand-tagged. **Do NOT hand-tag any stage.** Report only: is per-stage audit-status emission wired to a
generator? If yes and it merely needs re-running to reflect the completed audit, run it. If wiring it is
non-trivial, FLAG it as a Step-6 follow-up in your iteration report — do not build it under this part.

## Part A acceptance criteria

1. Both new appendices are emitted by committed, deterministic generators (re-run → identical output).
2. Both are exhaustive over their sources; every exclusion/grouping is visible (no silent drops).
3. A spot-sample of appendix rows matches the source YAML exactly (stage, anchor, constraint_kind,
   status).
4. Archive build (`pdflatex pde_ledger.tex` ×2) exits 0 with no new errors and warnings no worse than
   baseline. Reader build still exits 0.
5. Falsification table reflects layer-2 complete with accurate, non-overclaiming, toy-analog wording.
6. Tree left dirty; nothing committed. Iteration report lists files created/edited + the A-5 finding.

---

# PART B — Punch-list prose + metadata patches (apply after Part A reviewed)

**Canonical checklist:** `research/pde_ledger/redteam_adversarial/STEP6_CLOSEOUT_PUNCHLIST.md`. Apply
every item EXCEPT the OFF-LIMITS A8 / A9. The punch list carries the exact `file:line` anchor and the
fix for each item; follow it. Group the work as below. For each item, either APPLY it or, where the
punch list itself says "no fix required / informational only / log only," LEAVE it and record that in
your iteration report (do not invent edits for the audit-internal scanner-slug items).

**Non-negotiable for Part B:** these are all scoping / label / anchor / wording / metadata changes.
**None may alter what a claim asserts.** If applying any item would change the substance of a claim
(not just how it is scoped or labelled), STOP that item and flag it in your report as conceptual — do
not apply it.

- **B-i — Card-text overclaim scoping (A1–A7).** A1 (stage_192 re-attribute χ_Q/Δ_Q/N_Q to 194/195),
  A2 (stage_189 `\StatusExactClosure` + GR-target disclosure), A3 (stage_116:5 badge scoping), A4
  (stage_001:149 add the "ansatz / not yet frozen" hedge the notes carry), A5/A6 (stage_073 37/20 and
  ε_r=1/20 presentational tags), A7 (stage_115:1 "Theorem" wording, LOW). Apply A2 and A3 with one
  consistent badge-scoping rule (they are the same class). **A8/A9 are OFF LIMITS.**
- **B-ii — constraint_kind relabels (B1, B2)** in the provenance YAML. B1 (3 r_F1 bundles
  published_target→free_choice; USER-APPROVED 37/20=free_choice). B2 (fam_0212 V_known
  published_target→free_choice).
- **B-iii — Stale / misdirected stage-anchors (C1–C15).** Content-keyed numbering-drift fixes — for
  EACH item OPEN the cited notes/card/script and fix the specific stale token; **NEVER offset-sweep**
  (see `notes/LINEAR_STAGE_RENUMBERING_MANIFEST.json` for the arithmetic only; content decides). Items
  flagged "benign artifact / metadata-only" (e.g. C9/C11 introduced_at_stage encoding artifacts):
  correct if it is a real stale token in a file, else leave + log.
- **B-iv — Dedup / scanner metadata (D1, D2, D4, D7, D8).** Alias/dedup cross-reference notes in the
  provenance YAML. D3, D5, D6, D9 are explicitly **"no fix required"** audit-internal candidate_ids —
  leave, log only.
- **B-v — Disclosure notes (E1–E5).** Optional notes additions where a value lives only in a script
  (f'(0)=1, D_W_bare, the Ξ-prefactor/V_known/session-readback genealogies). Add the brief notes
  mention each item specifies.
- **B-vi — Script wording (F1, F2).** F1 (soften stage247 "independently derived / falsifiable"
  λ_L-closure wording). F2 (stage162 sympy script: replace the hardcoded r_F1 DECIMAL
  `1.77799353547498` with the closed form `sqrt(4107-100*pi**2)/(10*pi)`, matching the
  de-transcription already applied to 165/168/169). After F2, **re-run that script under `timeout 600`
  and confirm it still exits 0** (cosmetic form change must not break it).
- **B-vii — Completeness-critic items (CC-1, CC-3, CC-4).** CC-1 (downgrade the
  `fit_stage223_..._dual_engine_overclaim` bundle + its Phase-A provenance record from "dual_engine"
  to single-engine to match the honest card — provenance/record only, NO card/notes edit). CC-3
  (provenance-display completeness notes for 200 M_star / 228 K_compat / 094 c_geom). CC-4 (fill the
  EMPTY `provenance_findings` blocks for 126/127/128/130 — record the verified classification rationale;
  classifications already confirmed sound vs files). **CC-2 is handled by Claude in the catalog, not
  here.**
- **B-viii — Provenance citations (the "add citations" part of the close-out).** Where the audit found
  a card-level provenance gap, add the upstream anchor: notably C14 (stage_121 cites its Stage-073
  L/a=37/20 origin freeze), and any A-class item whose fix is "add the benchmark/provenance citation."

**G-class (G1–G3, Phase A/B tooling cosmetics):** LOWER priority, opportunistic, superseded by the
genealogy. Apply if quick; otherwise leave logged. Not required for close-out.

## Part B acceptance criteria

1. Every non-gated punch-list item is either applied or logged "left — <reason>" in the iteration
   report. Nothing silently skipped.
2. No substance change — scoping/label/anchor/wording/metadata only. Any item that would change what a
   claim asserts is stopped and flagged conceptual (not applied).
3. Numbering-drift fixes are content-keyed and verified against the cited notes before editing (no
   offset-sweep).
4. F2's script re-runs exit 0; the archive build still exits 0 (cards live in it).
5. A8 / A9 untouched.
6. Tree left dirty; nothing committed. Iteration report = per-item applied/left table.
