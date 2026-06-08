---
unit_id: 162
batch: IV.6
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-08T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 162

## Per-finding outcomes

### F1 — stale numbering reference (script comment)

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage162_parent_compensation_rigidity_sympy_audit.py:39`
- Before: `# Exact parent family formulas from Stages 99 and 102.`
- After:  `# Exact parent family formulas from Stage 119 (with gamma0 from Stages 115-116).`

The captured diff (`exec_logs/stage_162_diff.patch`) shows exactly this one-line change and nothing else in the stage-162 `.py`. The surrounding code (`gamma0 = (1 + r**2)/9`, `Lratio`, `g_lower`, and all four checks) is byte-for-byte unchanged.

**Assessment:**
Correct and complete.

1. **Only change is the comment at line 39.** Reading the current `.py` in full: lines 1–38 and 40–83 are identical to the pre-fix structure; the only altered line is the comment at 39. No code, expression, assertion, or numeric literal changed. `git diff --stat` reports `2 insertions/deletions` on the `.py` consistent with a single line replacement.

2. **`.wl` untouched, no prose touched.** The stage-162 `.wl` does not appear in `git status` (only the stage-162 `.py` does, plus unrelated stage-158/161 files belonging to other units). `grep` for any provenance comment in the stage-162 `.wl` returns nothing — it never carried such a comment, as the report noted, so leaving it untouched is correct. No `paper/*.tex` or `notes/` file for stage 162 is modified.

3. **grep gate passes.** `grep -n "Stages 99 and 102"` on the `.py` returns nothing (exit 1). The new comment cites Stage 119, and explicitly attributes gamma0 to Stages 115–116.

4. **New attribution is correct.** Cross-checked against the stage-162 notes (read for provenance sanity only):
   - notes:40 — "Stage 119 rewrote the compensated throat-core condition in terms of the normalized parent ratios"; notes:46–51 give the exact balance law `1+𝔯²=4(𝔤−𝔯)²`, `𝔤±=𝔯±½√(1+𝔯²)`; notes:52–56 box the `L_W/a=(π/2)√((1+𝔯²)/3)` selection law — all owned by Stage 119.
   - notes:58–62 — "And Stages 115–116 fix the bare odd normalization by `γ₀=(1+𝔯²)/9`."
   The three formulas typed in the `.py` (`gamma0`, `Lratio`, `g_lower`) are exactly these Stage 119 / 115–116 closed forms, so the new comment's provenance is accurate. The old "Stages 99 and 102" was +17-era numbering-drift residue, now corrected.

Non-tautological status of the assertions is unchanged (the edit is comment-only): A1–A3 / B1–B3 still independently differentiate/decompose and assert the residual vanishes.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `similarity identity = 0`
- `lower-branch differential law = 0`
- `positive slope decomposition = 0`
- `(dg/dr)|_F1 = 0.564199521046342514…`, `dr/dg |_F1 = 1.77242263188284677…`

All three `expect_zero` checks report `= 0`; the run reached the carry-forward prints and exited 0.

**Mathematica:** exit=0. Notable lines:
- `PASS: similarity identity`
- `PASS: lower-branch differential law`
- `PASS: positive slope decomposition`
- `Stage 162 Mathematica audit passed.`

Both engines agree on the three zero-residual identities and agree to ~16 sig figs on the numeric `(dg/dr)|_F1` and `dr/dg|_F1`.

**Output freshness:** The committed sympy `.txt` (mtime 2026-06-08 11:56:32) and mathematica `.txt` (11:56:32) are both newer than the edited `.py` (11:48:37) and the `.wl` (2026-05-27), so the saved outputs were re-generated post-fix. As predicted, the committed `.txt` carries none of the comment text (`grep -c "from Stage"` → 0; comments are never printed), so the output is byte-identical to the pre-fix output — confirming the comment change could not and did not alter printed results. The committed `.txt` still shows `similarity identity = 0`, `lower-branch differential law = 0`, `positive slope decomposition = 0`.

## Material-change assessment

`material_change`: false.

The edit is a single non-load-bearing source comment. It feeds no assertion, changes no expression, and produces byte-identical output. No derived result changes, so no downstream unit (> 162) depends on anything that moved. This is pure provenance hygiene.

## Side observations (non-blocking)

- The working tree also carries modifications to stage-158 and stage-161 `.wl`/output files and the pass-2 `MANIFEST.yaml`; these belong to other audit units and are correctly outside the stage-162 diff patch (which contains only the one `.py` line). Not a stage-162 concern.

## Verdict justification

The sole finding (F1, low-severity stale-numbering comment) is resolved exactly as directed: the one comment line at `.py:39` was rewritten to cite Stage 119 plus Stages 115–116 for gamma0, with no collateral edits, no `.wl`/prose changes, and a confirmed-correct attribution against the notes. The grep gate passes, both engines still exit 0 with all three identities at zero, and the freshly re-generated outputs are byte-identical to the pre-fix transcripts (a comment cannot alter printed output). No regressions in the diff. Verdict: `verified`, `material_change: false`.
