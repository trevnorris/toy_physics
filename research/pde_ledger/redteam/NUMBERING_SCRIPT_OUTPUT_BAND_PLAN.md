# Numbering reconciliation — SCRIPT / OUTPUT band (deferred dedicated pass)

**Status: ✅ COMPLETE 2026-06-10 — all three bands done, verified, and committed**
(`b25cb57` 001–090, `7a94cbc` 091–180, `857e255` 181–253). This was a dedicated pass,
SEPARATE from the red-team second pass (NOT inside the red-team fix loop), user-approved
2026-06-04 and unblocked once the red-team SECOND PASS completed (253/253, 2026-06-10).
Content-keyed, never offset-swept; every band's edits were digit-only (strip-digits proof —
HEAD vs working byte-identical after removing all digits, line counts unchanged), both engines
re-run exit 0, output diffs intended-only. The interim "defer to this plan" policy (§ below)
has ended. Discovered during pass-2 batches I.2 (stage 021) and II.1 (024–036); deferred items
accumulated across ALL pass-2 batches I.2→VIII.1 (per-batch inventories live in
`redteam/pass2/batches/batch_<ID>.md` + the handoff PASS-2 PROGRESS LOG). Companion to
`NUMBERING_BROAD_SWEEP_PLAN.md` (which covered the **notes** band, DONE). See memory
`[[numbering-drift-root-cause]]`.

> **Per-band records** (full evidence, applier, verification): band 001–090 and 091–180 in the
> commit bodies + `redteam/script_output_maps/` maps; band 181–253 in
> `redteam/script_output_maps/BAND_181-253_RESULT.md` (+ `apply_band_181_253.py`).
>
> **2 escalations LEFT for the user** (out of label-only scope — content/identifier judgment, NOT
> digit-only fixable):
> 1. **stage182.py:45** `# Stage-30 coherent branch definitions.` — ambiguous between an
>    old-epoch→047 (`coherent_kernel_map`, 30+17) reading and a loose concept-pointer; the vars
>    below it are defined locally at stage182. LEFT. Route to content review, or leave as provenance.
> 2. **stage203.wl** variable family `chiFromStage180` / `closureNumStage180` — the `.py` twin names
>    the same expression `chi_from_stage197` / `closure_num_stage197`; canonical owner is **197**
>    (`conditional_packetA_closure_theorem`; 197−17=180, a stale OLD-epoch number baked into a CODE
>    IDENTIFIER). A digit rename here changes a *variable* → OUT OF SCOPE for this label-only pass.
>    Recommend a Codex identifier-rename 180→197 in stage203.wl for cross-engine consistency.

> **One-line:** the "mechanically EXHAUSTED (001–253)" numbering reconciliation cleared the
> **notes** band but NOT the **script (`.py`/`.wl`) + committed-output (`.txt`)** band. This
> plan finishes that band — content-keyed, **never offset-sweep**, doc/label-only, math untouched.

---

## Why this exists (root cause)

The deterministic Phase-1 reconciliation (`e2a4780`, 2026-06-03) fixed `.py`/`.wl` **main
`banner(...)` headers + module docstrings** by keying on those specific forms, and the notes
broad/secondary sweeps fixed the notes band. But:

1. **Committed `.txt` OUTPUTS were never regenerated band-wide.** Many were captured
   *before* `e2a4780`, so they still print the pre-renumber self-banners even though the
   *source* banner is now canonical. Re-running the script refreshes them.
2. **The scan key missed several script self-label forms**, so stale tokens survive in source:
   - bare `print("All Stage NN checks passed.")` / `Print["All Stage NN checks passed."]`
     (prints INTO the transcript → a plain re-run alone leaves it stale);
   - module-docstring self-labels (`Moving-throat PDE Stage NN SymPy audit.`) and old
     **filename** self-refs in the docstring (`...stage8_...py`);
   - ledger / completion lines: `FINAL STAGE-007 LEDGER:`, `STAGE 13 AUDIT COMPLETE`;
   - **hyphenated** forms `STAGE-007` (the scan keyed on `STAGE␣NNN` with a space);
   - **sub-stage** comment refs `Stage 16.1 / 16.6`.
3. **Some refs are genuinely AMBIGUOUS** (self vs cross-ref) because a stage carries
   self-labels from MORE THAN ONE renumber epoch — e.g. stage024 shows both `Stage-6` and
   `STAGE-007` for itself, and `Stage-6 grouped bundle` might instead be a cross-ref to
   stage 023. These need per-reference content adjudication, exactly like the Class-D notes work.

This is the SAME multi-epoch, content-dependent drift as the notes band. **The canonical
numbering (filenames / MANIFEST / paper cards) is ground truth; NEVER offset-sweep.**

## Confirmed examples (from pass-2 I.2 / II.1 — non-exhaustive)

- **Stale committed OUTPUTS with pre-renumber self-banners (on stages the pass-2 auditor
  marked *clean*):** stage028 out `STAGE 011 — LOADED PROFILE SELECTION`; stage030 out
  `STAGE 13 AUDIT COMPLETE`; stage031 out `STAGE 14 AUDIT COMPLETE`; stage025 out
  `Stage 8 Mathematica audit passed.`; stage026/027 out `STAGE-8`/`STAGE-8/9` cross-refs +
  self `Stage 9`/`Stage 10` pass-lines.
- **Ambiguous self-vs-cross multi-epoch source refs (LEFT untouched in II.1):** stage024
  `.wl:293 FINAL STAGE-007 LEDGER` + `.py:10/472 "Stage-6 … grouped bundle / normalization
  stack"`; cross-refs stage033 `.py:50 Stage-15` (=stage032), stage036 `.py:106/127 Stage-18`
  (=stage035).
- **Already fixed in-loop (the unambiguous self-labels the red-team flagged)** — DO NOT redo:
  I.2 stage021 `.wl:195`; II.1 stage032/033/034/035/036 docstring + closing pass-line (and
  032 `.wl` Print, 033 `.wl` sub-stage comment), all `Stage {15..19}`→`{32..36}`; II.1 stage024
  notes typo `4π/122`→`4π/105` (a value typo, not a label — already done).
  **III.1 (037–048):** docstring/closing self-labels canonicalized (incl. 040:136/041:107
  inline-comment self-labels the arbiter grep caught); 047:121 `Stage-30` left AMBIGUOUS →
  this plan. **III.2 (049–060):** the 11 self-label stages + 049 `.wl:93` closing. **III.3
  (061–072):** 061–066 docstrings, 062/063 closing pass-lines, 068/070/071 filename-style
  docstrings (`stage51/53/54`→`stage068/070/071`), 070 `py:5` prose. The committed OUTPUTS
  for ALL of 061–072 were also refreshed in-loop (their self-banners are now canonical — the
  only remnants are deferred CROSS-refs, below).
- **Per-batch DEFERRED inventories (cross-refs / variable-names / ambiguous, LEFT untouched):**
  the exhaustive list for each pass-2 batch lives in its summary `redteam/pass2/batches/batch_<ID>.md`
  (§"Deferred to the dedicated SCRIPT/OUTPUT-band numbering pass") + the handoff PASS-2 PROGRESS
  LOG. III.3 examples: 063 `py:76`, 064 `py:25/122/180`+`wl:104`, 065 `py:22`, 066 `py:14/59`,
  069 `py:9/11/26/94/112/120/176/179`+`wl:99/116/164/167` (`py:26` "Stage 049" is the same ref
  in 3-digit form — ambiguous), 070 `py:57/59`+`wl:87/88` + the variable names `J1_stage48`/`J1Stage48`.
- **VI.1 (201–218), pass-2 batch close 2026-06-09 — 3 deferred items (content-keyed, NEVER offset-sweep):**
  (1) **206 `.py`** collapse-target label "Stage 239" should be the **Stage-205 quad-log predictor** collapse
  target (`py:131/133/142/145`) — the MATH is correct; only the cross-ref stage number is stale
  (cross-ref → content-map to the owning Stage 205). (2) **218 `.py`** dead `source_stage:198/200` dict
  fields + "Stage 249" comment (`py:347`) — NON-asserted (dead metadata + a comment), so label-only/cosmetic.
  (3) **203 `.wl`** variable names `chiFromStage180`/`closureNumStage180` carry a stale "Stage 180" attribution
  vs the correct **Stage-197** owner (variable-name self/cross-ref drift — adjudicate by the cited deliverable's
  owning stage). All three are doc/label-only with math untouched; LEFT for this dedicated pass.

## Scope

Project-wide over the **script + committed-output layer only**: `scripts/*.py`,
`mathematica/*.wl`, `scripts/output/*.txt`, `mathematica/output/*.txt`. The notes/`.tex`
band is already done (see `NUMBERING_BROAD_SWEEP_PLAN.md`); do not re-open it here.

## Method (per [[numbering-drift-root-cause]]: content-keyed, never offset-sweep)

**Phase A — INVENTORY (read-only, fan out to agents).** Build a worklist:
- A1. Every committed `.txt` whose printed banner/labels ≠ the stage's canonical number
  (compare output banner to `paper/stages/stage_NNN.tex` / filename canonical). These are
  "stale output, source may or may not be canonical."
- A2. Every `.py`/`.wl` SELF-label (docstring, filename self-ref, `print/Print` pass-line,
  ledger/AUDIT-COMPLETE line, banner, sub-stage `NN.k`) whose number ≠ canonical — using a
  BROAD regex that catches `Stage`/`STAGE` + any separator (space, `-`, `_`, none) + number,
  AND the old-filename `stageN_` stem form. (The narrow `STAGE␣NNN` key is what missed these.)
- A3. Every `.py`/`.wl` CROSS-ref to another stage whose number is stale (e.g. `Stage-15`
  inside stage033 = stage032). Map by described-deliverable→owning-stage, NOT by offset.
- A4. Flag the AMBIGUOUS self-vs-cross / multi-epoch cases for individual content review
  (e.g. stage024 "Stage-6").

**Phase B — REFRESH stale outputs (orchestrator).** For every A1 stage whose *source*
banner is already canonical: re-run both engines and refresh the committed `.txt`
(`$RT exec-sympy/exec-mathematica` then `sed '1,/^---$/d;/^# exit_code:/d' <log> > <out>`,
≤2 Mathematica seats, sequential). This alone fixes outputs that were stale only because
they predated `e2a4780`.

**Phase C — FIX unambiguous stale SOURCE self-labels (A2), then re-run.** Label-only:
change ONLY the number, matching each file's OWN canonical banner format (`.py` banners are
2-digit `STAGE NN.x`; `.wl` banners are 3-digit `STAGE 0NN` — match locally). Then re-run to
refresh the output (so pass-lines/ledger lines in the `.txt` become canonical too).

**Phase D — CONTENT-MAP A3 cross-refs + A4 ambiguous refs, per reference.** Exactly like the
Class-D notes pass: adjudicate each by content (which stage owns the cited deliverable?),
never by a uniform offset; leave genuine current-canonical refs alone; correct only the stale
ones; flag any truly unresolvable as "cannot-confirm → leave."

**WHO APPLIES.** This is a DEDICATED NUMBERING pass, not the red-team loop → the
numbering-reconciliation convention applies: **doc-only label/banner/string edits are
ORCHESTRATOR/agent-applied directly**, even inside `.py`/`.wl`, because no CODE/MATH changes
(per [[numbering-drift-root-cause]] "WHO APPLIES" + [[feedback-codex-is-fix-applier]] scoping).
Output refresh is the orchestrator exec re-run. (Inside the red-team loop, by contrast, Codex
applies in-loop self-label fixes — that distinction is intentional.)

## Verification (every phase)

- **Strip-all-digits / placeholder-collapse proof**: old vs new byte-identical except the
  Stage-citation number tokens (the c2/c3 appliers' guard — reuse `redteam/apply_*` style).
- Both engines **exit 0**; every `PASS:` / `= 0` residual unchanged; NO equation/value/variable
  /assertion/`\label`/logic byte changed. Spot-confirm a handful of `.py`/`.wl` run clean.
- Per-file residual-0 (no stale Stage token of the wrong number remains).
- Gate between batches (fan-out by band, e.g. 001–090 / 091–180 / 181–253), per
  `[[feedback-sequential-audit-chunks]]`.

## Cautions / out of scope

- **NEVER offset-sweep.** Offsets are multi-epoch and inconsistent (0/+17/+34/… and within a
  single file two epochs can coexist, as stage024 shows). Content decides every ref.
- Do NOT touch non-stage numbers (physics constants/literals like `Lambda_ell=37`, kernel
  symbols `\mathcal S_n`), `\label{}`, or any math.
- Do NOT re-open the notes/`.tex` band (done).
- A "cannot-confirm" ref → leave + log (e.g. the notes `stage229:5 "Stage 143/093"` already
  flagged in `NUMBERING_BROAD_SWEEP_PLAN.md`).

## Interim policy WHILE the red-team second pass continues (III.1 onward) — ⭐ ENDED 2026-06-10 (second pass COMPLETE; HISTORICAL, retained for context)

Until this dedicated pass runs, the red-team loop keeps the established policy: for a
findings stage, fix ONLY the audit-flagged UNAMBIGUOUS self-labels (matching the file's
canonical banner) + refresh its outputs; **defer** all other script/output numbering drift
(clean-stage stale outputs, ambiguous/cross refs) to THIS plan. The orchestrator's arbiter
grep on the source is the backstop, because the pass-2 audit agents apply an inconsistent
threshold to numbering labels (they flag some stages and miss the same class on clean ones).

## Pointers

- Companion (notes band, DONE): `redteam/NUMBERING_BROAD_SWEEP_PLAN.md`,
  `NUMBERING_RECON_CLASSD_MAP.md`, `NUMBERING_RECON_SCAN.md`.
- Applier style to reuse: `redteam/apply_secondary_maps.py`, `redteam/apply_cluster3.py`,
  `redteam/pad_neighbor_refs.py` (deterministic, per-edit-unique, strip-digits guard).
- Authoritative renumber arithmetic (for corroboration only, never to sweep):
  `notes/LINEAR_STAGE_RENUMBERING_MANIFEST.json` (linear +17).
- Memory: `[[numbering-drift-root-cause]]`, `[[project-full-second-pass]]`,
  `[[project-moving-throat-verification]]`.
