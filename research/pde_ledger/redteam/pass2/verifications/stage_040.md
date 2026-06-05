---
unit_id: 040
batch: III.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 040

## Per-finding outcomes

### F1 — stale_output (self-label) [label-only]

**Classification:** resolved

**What changed:**
SymPy source `scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py`:
- line 3: docstring header `Moving-throat PDE — Stage 23 SymPy audit.` → `Stage 40 SymPy audit.`
- line 55: subbanner `23.1` → `40.1` (rest of string verbatim)
- line 81: subbanner `23.2` → `40.2`
- line 101: subbanner `23.3` → `40.3` (trailing `... continuum of Stage 22` cross-ref left intact)
- line 124: subbanner `23.4` → `40.4`

**Assessment:**
Correct and strictly label-only. Each edit touched only the leading stage-number token; the remainder of every line is byte-identical to the pre-fix text. The already-canonical `STAGE 40` banners (lines 53, 151) were left alone. No equation, value, or assertion changed — verified by reading the full current file: `lam_minus`, `alpha_req`, the eigenvector residual, `F_expected`/`G_expected`, the split-U map (`q_U=-sqrt(2/9)R_U`, `eta_U=(2/9)R_U`), `F_U`/`G_U`, the `F_stage18`/`G_stage19` recovery literals, and the `H_F`/`H_G` dual-route block are all unchanged. No collateral edits.

### F2 — stale_output (self-label, missed comment) [label-only]

**Classification:** resolved

**What changed:**
line 136 inline comment `... F_general, G_general (section 23.2),` → `(section 40.2),`. This is a self-reference to this stage's own subsection (now subbanner `40.2`). The captured `stage_040_diff.patch` was truncated before line 136, but the current file confirms the edit is in place.

**Assessment:**
Correct, label-only (only the numeric token changed). The comment is not printed, so it does not affect the transcript; the edit is internal-consistency only. `F_general`/`G_general` variable names in the same comment were correctly left untouched.

## Cross-ref / variable-name preservation check

All DO-NOT-TOUCH items confirmed intact in the current file:
- docstring cross-refs `Stage-18` (L9), `Stage 22` (L12), `Stage-18/19` (L14) — unchanged.
- ledger cross-refs `Stage-18` (L152), `Stage-18/19` (L160) — unchanged.
- `... of Stage 22` tail on subbanner L101 — unchanged.
- variable names `F_stage18` (L118,121), `G_stage19` (L119,122) — NOT renamed.
- provenance comments citing `stage035`/`stage036` paths (L114,116) — unchanged.

## Exec log assessment

**SymPy:** exit=0. Refreshed output (dated 2026-06-05T08:58, post-fix) now reads `STAGE 40 — ...` (L8), subbanners `40.1`/`40.2`/`40.3`/`40.4` (L12,21,31,39), and `STAGE 40 THEOREM LEDGER` (L49). All residuals `= 0`:
- `e1/e0 closed form = 0`, `eigenvector residual row 0/1 = 0`
- `F_general - expected = 0`, `G_general - expected = 0`
- `F_U(R_U=1) - Stage18 F = 0`, `G_U(R_U=1) - Stage19 G = 0`
- `H_F cross-check = 0`, `H_G cross-check = 0`

**Mathematica:** exit=0. The `.wl` was untouched (directive forbade editing it). Refreshed output (dated 2026-06-05T08:09) header `STAGE 040 — ...`, all eight checks `PASS`, closing `Stage 040 Mathematica audit passed.` The `.wl` uses numeric section headers (`1.`–`4.`) rather than `40.k` subbanner labels — consistent with it being already canonical and untouched. All result lines (`alpha_req`, `F_(q,eta)`, `G_q`, `F_U`, `G_U`, `H_F`, `H_G`) are algebraically identical to the SymPy transcript.

**Output freshness:** Both `.txt` were re-generated post-fix — SymPy stamp 2026-06-05T08:58, Mathematica 2026-06-05T08:09, both newer than the prior 2026-05-22 committed transcripts the auditor flagged. The stale `STAGE 23`/`STAGE 023` headers are gone; only banner labels changed, every result line matches the prior values, exactly as the finding's verification criterion required.

## Material-change assessment

`material_change`: false. The edits are stage-number label tokens in a docstring header, four subbanner strings, and one inline comment. No derived quantity, assertion, constant, or variable changed; both engines still pass all checks with identical numeric/symbolic results. Downstream units that depend on stage 040's results (`F_{q,eta}`, `G_q`, `F_U`, `G_U`, `H_F`, `H_G`) are unaffected.

## Side observations (non-blocking)

- The captured `stage_040_diff.patch` is truncated (ends mid-context at the L124 subbanner hunk, does not include the L136 F2 hunk). The file itself confirms F2 was applied, so this is a diff-capture artifact, not a missing edit — non-blocking.
- The `.wl`/SymPy transcripts use different section-header conventions (numeric `1.`–`4.` vs `40.k`). This is pre-existing and orthogonal to this fix; no action.

## Verdict justification

Both label-only findings are resolved. The diff (corroborated against the full current file) is strictly stage-number-token edits: docstring header, four subbanners, one inline comment — no equation, value, assertion, or variable rename, with `F_stage18`/`G_stage19` and all `Stage-18`/`Stage 22`/`Stage-18/19` cross-refs correctly left untouched. The refreshed transcripts now show `STAGE 40`/`STAGE 040` with `40.1–40.4` subbanners, every prior `= 0`/PASS residual remains, and both engines exit 0. No regressions; `material_change: false`.
