# Pass-2 augmentation — exhaustive script→doc value reconciliation

**User decision (2026-06-04):** the full SECOND PASS adds, on top of the standard
v2 audit, an explicit per-stage step that accounts for **every** result value the
SymPy/Mathematica scripts emit and confirms each is correctly reflected in the
stage paper card (`.tex`) and the per-stage notes (`.md`). The standard v2 audit
already flags load-bearing `paper_misalignment`s; this makes the script→doc
reflection check **exhaustive** rather than load-bearing-only.

Every pass-2 audit agent reads this file IN ADDITION to its rendered
`audit_prompt_<NNN>.md`. Do the standard audit exactly as that prompt specifies,
THEN do the reconciliation below.

## What to add to your report

Append a section titled `## Value Reconciliation (pass-2 augmentation)` to your
report (`redteam/pass2/reports/stage_<NNN>.md`), after the standard sections.

### Procedure

1. **Enumerate every RESULT/deliverable value the scripts emit.** Use the script
   source (`.py`, `.wl`) together with the committed saved outputs
   (`scripts/output/..._sympy_audit.txt`, `mathematica/output/..._mathematica_audit.txt`)
   as the authoritative record of what they produce. Do **not** run anything.
   Enumerate:
   - every named numeric constant the script computes or pins (e.g. `Pi_star =
     1.50882951…`, `gamma_0 = …`, benchmark/figure-of-merit numbers);
   - every boxed / closed-form **symbolic** result the script derives or asserts;
   - every value the script PRINTS as a *labeled result*.
   EXCLUDE verification scaffolding: pass/fail flags, residual-near-zero check
   values (`residual = 1e-58`), tolerances, iteration counts, timings, and
   intermediate quantities that exist only to drive an assertion.

2. **Locate each value in the docs** — the `.tex` card and the `.md` notes
   (whichever is the natural carrier). Classify each as:
   - **MATCH** — appears in `.tex` and/or `.md` and agrees. Give `file:line`.
   - **MISMATCH** — appears but the value differs. → raise a `paper_misalignment`
     finding, subtype `value_mismatch`. Quote BOTH sides verbatim with `file:line`.
   - **MISSING-DELIVERABLE** — the value is one of the stage's STATED deliverables
     (a boxed result, a `\stagefield{Output}` quantity, or a key constant the card
     / notes is meant to report) yet is absent from BOTH `.tex` and `.md`. → raise
     a `paper_misalignment` finding, subtype `script_missing_paper_claim`.
   - **INTERNAL** — genuine internal scaffolding not expected in prose. Account for
     it but raise NO finding.

3. **Present a table** in the new section:
   `| value | source (py/wl + output line) | .tex/.md location | status |`
   List the INTERNAL items separately as a short name-only list, so the main table
   is the deliverable-level reconciliation.

4. **Fold genuine problems into the standard findings.** Every MISMATCH and
   MISSING-DELIVERABLE is a `paper_misalignment` finding in the main `## Findings`
   section (ordered first, as the standard prompt requires), counts toward
   `findings_count`, sets `verdict: findings` and `needs_user_resolution: true` —
   even if the standard audit was otherwise clean. They route to the user gate
   exactly like standard `paper_misalignment`s; Codex applies nothing until the
   user resolves direction.

5. If every emitted deliverable value reconciles, say so explicitly and keep your
   standard verdict. Add a one-line `reconciliation: complete; N values checked, 0
   misaligned` note in the new section.

### Guards (avoid false-positive flooding)

- A terse `.tex` card legitimately omits intermediate quantities. Do NOT mark a
  value MISSING just because the card omits it — only **deliverables absent from
  BOTH** the card and the notes count. If the value lives correctly in the `.md`
  notes, that is a MATCH.
- Direction is never decided here. A MISMATCH only flags the discrepancy and
  quotes both sides; the user (with Claude+Codex) decides which side is right.
- Output `.txt` missing/stale: note it (it is the standard `stale_output` signal);
  base the reconciliation on the script source + your reasoning and say so.
- Status-only stages (no scripts, or Mathematica-only by design): reconcile
  whatever the present engine emits; note the single-engine status.
