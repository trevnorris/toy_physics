---
unit_id: 024
batch: II.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-04T23:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 024

## Per-finding outcomes

### F1 — paper_misalignment (value_mismatch)

**Classification:** resolved

**What changed:**
A single notes-only edit at `notes/stages/moving_throat_pde_stage024_overlap_isotropy.md:213`:
`(4 pi / 122)` → `(4 pi / 105)` in the sixth-moment-of-the-unit-sphere identity. The
git diff of that file shows exactly one hunk: one `-` line and one `+` line, differing only
in `122`→`105`; nothing else on the line or elsewhere in the notes changed. No script source
was touched — `git diff` of both the `.py` and `.wl` is empty, and the captured script diff
patch (`exec_logs/stage_024_diff.patch`) is empty, as expected.

**Assessment:**
Correct and complete. The fix matches the RESOLVED-F1 / Applied-F1 directive blocks exactly
(direction (a): notes-side typo, no script change). The corrected value `4π/105` is the
textbook 2-sphere 6th-moment prefactor `1/(3·5·7)=1/105`, which is what both engines already
use (SymPy py:128 `4*pi*s/105`; Mathematica direct surface integral wl:42-46) and which the
notes' own downstream `κ_* = √5/(7√π)` (notes:221) requires. The published `.tex` card is
silent on this prefactor, so the published card is unaffected. No collateral edit. No
assertion was changed, so no tautology concern.

### F2 — stale_output

**Classification:** resolved

**What changed:**
Both stage-024 transcripts were refreshed by the orchestrator re-run. The committed
`mathematica/output/...stage024...mathematica_audit.txt` diff is a single hunk: the top
banner line 2 `STAGE 007 — OVERLAP ISOTROPY` → `STAGE 024 — OVERLAP ISOTROPY`. The SymPy
transcript content is unchanged (all `= 0`). Output mtimes are now newer than the scripts.

**Assessment:**
Correct. Matches the Applied-F2 block (re-ran both engines under `timeout 600`, banner
corrected `STAGE 007`→`STAGE 024`, all checks still pass). No source edit was required and
none was made. Informational finding fully closed by the standard re-run.

## Exec log assessment

The orchestrator-captured `exec_logs/stage_024_sympy.log` and `stage_024_mathematica.log`
are empty/zero-length (the orchestrator routes the live run into the committed
`scripts/output/...` and `mathematica/output/...` transcripts, which both carry the
`# exit_code: 0` footer and an `argv` showing `timeout 600`). I read those committed
transcripts in lieu of the empty logs.

**SymPy:** exit=0. Notable lines:
- `Gram - I5 =` → zero 5×5 matrix (out 24-33).
- `M - M_target =` → zero 5×5 matrix, with `M^(20)` printed as
  `√5/(7√π) diag(1, 1/(2), 1/(2), -1, -1)` (out 131-160) — κ_*/M^(20) unchanged.
- All `B/Z/N/D` formula residuals `= 0`; all transport residuals `= 0`
  (`b_x - 3 a_x = 0`, `b_P - 3 a_P = 0`). No FAIL.

**Mathematica:** exit=0. Notable lines:
- `PASS: Gram - I5`, `PASS: M - M_target` (out 15, 91).
- `PASS: Z_full from matrix inverse matches paper rational`,
  `PASS: N_full from matrix inverse matches paper rational` (physics anchor, out 46-49).
- `PASS: Lane-breaking witness: collapse check is non-tautological` (out 85).
- Every line is `PASS`; no FAIL. Banner reads `STAGE 024 — OVERLAP ISOTROPY` (out 9).

**Output freshness:** confirmed. Script mtimes are `2026-06-03 15:59:11` (both `.py` and
`.wl`, unchanged); both output `.txt` mtimes are `2026-06-04 23:13:28`, newer than the
scripts. Outputs were re-generated post-fix.

## Material-change assessment

`material_change`: false.

The only edits are (i) a notes-side typo correction to a value the engines already computed
correctly (`4π/105`), and (ii) a transcript-banner refresh. No script source changed, no
assertion changed, and no derived result changed — κ_* = √5/(7√π) and
`M^(20)=√5/(7√π) diag(1,1/2,1/2,-1,-1)` are bit-identical to the pre-fix run, and all
residuals remain 0. Nothing downstream of unit 024 can have been invalidated by this fix.

## Side observations (non-blocking)

- The Mathematica `.wl` source line 293 (`Print["FINAL STAGE-007 LEDGER:"];`, surfacing at
  output line 121 `FINAL STAGE-007 LEDGER:`) still carries a stale `STAGE-007` label. This is
  NOT a regression — the committed output-txt diff shows this line unchanged (it was
  `STAGE-007` before and after the re-run; only the top banner at line 2 was updated). It is a
  pre-existing cosmetic label-only banner string with no effect on any assertion or result,
  and it is precisely the kind of stale "Stage NNN" label that is tracked separately under the
  known numbering-drift workstream. F2 as scoped (refresh stale transcripts so they match the
  current script + canonical top banner) is satisfied; flagging this purely so the orchestrator
  can fold the leftover `.wl:293` label into the numbering-drift cleanup if desired. Not a basis
  to fail verification, and out of the auditor-only scope to raise as a new finding.

## Verdict justification

Both findings are `resolved`. F1's authorized notes-only typo correction is in place exactly
as directed (`122`→`105`, single hunk, no collateral edit, no script change — empty diff
patch), and the corrected value matches both engines and the notes' own κ_*. F2's transcript
refresh is in place (banner `STAGE 007`→`STAGE 024`, fresh mtimes, content otherwise
unchanged). Both engines exit 0 with all checks `PASS`/`= 0`, and κ_*/M^(20) are unchanged,
so `material_change: false`. No regressions in the diff. Verdict: verified.
