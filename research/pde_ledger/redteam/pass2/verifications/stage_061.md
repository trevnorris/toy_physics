---
unit_id: 061
batch: III.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T14:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 061

## Per-finding outcomes

### F1 — stale_output (stale self-label, numbering)

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage061_microscopic_gain_thresholds_sympy_audit.py:3` — the module
docstring self-label was changed from `Stage 44 SymPy audit — microscopic gain thresholds and
operator phase diagram.` to `Stage 61 SymPy audit — microscopic gain thresholds and operator
phase diagram.` (diff.patch lines 8-9). The captured diff shows exactly one changed line; the
removed and added lines are byte-identical once the digits are stripped:
`-Stage  SymPy audit — …` / `+Stage  SymPy audit — …`. NUMBER-only, not padded to 3 digits, as
the directive required.

**Assessment:**
Correct and complete. This is a strip-the-number-identical, label-only edit confined to the cited
line 3 with no collateral change. The directive specified `44`→`61` only, preserving "Stage" and
the rest of the docstring verbatim — exactly what was applied. No banner, cross-reference,
variable name, or `.wl` was touched (the banners were already canonical: py emits
`STAGE 61`/`STAGE 61 THEOREM LEDGER`, wl emits `STAGE 061`). The original finding was the
informational stale_output (committed `.txt` transcripts predated the banner fix and showed
`STAGE 44`/`STAGE 044`); it is closed both by this self-label fix and by the orchestrator's
re-run refreshing the two committed `.txt` outputs. A stale-label scan over the 061 `.py`, `.wl`,
and both `.txt` outputs finds no remaining `Stage 44`/`044` token. Not tautological — the edit is
documentation-only and does not alter any assertion; all algebraic/limit assertions in both
engines are unchanged and still pass.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 61 — MICROSCOPIC GAIN THRESHOLDS` (banner now canonical)
- `Xi_micro - kappa G_micro = 0`, `kappa*G_fail soft-support limit - Pe_req = 0`,
  `kappa*G_suff soft-support limit - 2 Pe_req = 0` (key identities residual-zero)
- `STAGE 61 THEOREM LEDGER` then `# exit_code: 0`

**Mathematica:** exit=0. Notable lines:
- `STAGE 061 — MICROSCOPIC GAIN THRESHOLDS` (banner canonical)
- `PASS: Xi_micro - kappa G_micro`, `PASS: G_fail^(inf) formula`,
  `PASS: stiff-support compliant-mouth limit: G_suff -> Pe_req`
- `Stage 061 Mathematica audit passed.` then `# exit_code: 0`
- The `Limit::alimv` warnings are pre-existing benign Mathematica notices about assumptions
  involving the limit variable; they do not affect any residual and are unrelated to this edit.

**Output freshness:** confirmed. Both committed `.txt` outputs (mtime 2026-06-05 13:58:32) are
newer than the edited `.py` (13:43:23) and the untouched `.wl` (2026-06-03 15:59:11). The sympy
`.txt` now reads `STAGE 61 — MICROSCOPIC GAIN THRESHOLDS` (line 3) and `STAGE 61 THEOREM LEDGER`
(line 38); the mathematica `.txt` reads `STAGE 061 — MICROSCOPIC GAIN THRESHOLDS` (line 3). Both
end `# exit_code: 0`.

## Material-change assessment

`material_change`: false.

The sole edit is a documentation-string stage number (`44`→`61`) plus the orchestrator's banner
refresh in the committed transcripts. No symbol, assertion, derived form, or numeric/limit result
changed. The refreshed transcripts carry the same algebra and the same residual-zero / PASS lines
as before; only the stage-label banners moved from stale `STAGE 44`/`STAGE 044` to canonical
`STAGE 61`/`STAGE 061`. Nothing downstream can depend on a docstring or a banner label, so no
unit > 061 is affected.

## Side observations (non-blocking)

- The git working tree shows pending changes for the whole III.3 band (062–072) as expected for
  an in-progress batch; for unit 061 specifically the only source change is the 061 `.py`
  docstring and the two 061 `.txt` outputs. The `stage070` `.wl` change in the working tree
  belongs to a different unit, not 061. No collateral 061 change.

## Verdict justification

The single low-severity finding (F1, stale_output / stale self-label) is fully resolved: the
docstring self-label was canonicalized `Stage 44`→`Stage 61` as a strip-the-number-identical,
label-only one-line edit with no collateral change and no `.wl` touch, and the orchestrator's
re-run refreshed both committed `.txt` outputs to the canonical `STAGE 61`/`STAGE 061` banners.
Both engines exit 0 with all algebraic and limit assertions still passing (the assertions are
unchanged, so non-tautological as established in the audit). No stale `44`/`044` label remains.
Outputs are fresh. No regression. material_change is false. Verdict: verified.
