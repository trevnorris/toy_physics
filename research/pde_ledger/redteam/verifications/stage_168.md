---
unit_id: 168
batch: V.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-28T16:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 168

## Per-finding outcomes

### F1 — stale_output (cosmetic banner mislabel)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage168_off_bundle_slippage_sympy_audit.py:30` — banner literal changed from `STAGE 151 — EXACT OFF-BUNDLE SLIPPAGE DECOMPOSITION` to `STAGE 168 — EXACT OFF-BUNDLE SLIPPAGE DECOMPOSITION`.
- `mathematica/moving_throat_pde_stage168_off_bundle_slippage_mathematica_audit.wl:26` — banner literal changed from `STAGE 151 — …` to `STAGE 168 — …`.

The `git diff` shows exactly these two one-line changes and nothing else; `151` → `168` only inside each banner string. No other line in either script was touched.

**Assessment:**
The edit matches the directive's "required change" precisely. The directive's `## Applied: F1` block lists both files and reports `deviation: none`, which the diff confirms — no collateral edits. Both regenerated transcripts now open with `STAGE 168 — EXACT OFF-BUNDLE SLIPPAGE DECOMPOSITION` (line 3 of each `.txt`; the report's "line 11" cite referred to the in-transcript header position before re-render, but the substance — the banner now reading "STAGE 168" — is what matters and is present). This is a label-only change: no assertion, symbol, RHS form, or numeric value was modified, and all five named checks still report zero residual / `PASS`. The fix is in place, correct, and confirmed by the post-fix exec logs. Non-tautological status is unaffected — the auditor already established all seven assertions are genuine cancellations, and the banner edit does not touch any of them.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 168 — EXACT OFF-BUNDLE SLIPPAGE DECOMPOSITION` (corrected banner)
- `bundle tangency (delta_perp on exact lower branch) = 0`
- `delta_perp + eps_perp = 0`
- `deltaPi transport identity = 0`; `dC + 16 sigma eps_perp / s = 0`; `dE2 … = 0`; `dE4 … = 0`; `DeltaQ … = 0`
No `AssertionError` / traceback appears, and `expect_zero` (L27-28) raises on any non-zero residual, so the clean run with all `= 0` lines confirms exit 0.

**Mathematica:** exit=0. Notable lines:
- `STAGE 168 — EXACT OFF-BUNDLE SLIPPAGE DECOMPOSITION` (corrected banner)
- `PASS: bundle tangency (delta_perp on exact lower branch)`
- `PASS: delta_perp + eps_perp`; `PASS: deltaPi transport identity`; `PASS: dC …`; `PASS: dE2 …`; `PASS: dE4 …`; `PASS: DeltaQ …`
- `Stage 168 Mathematica audit passed.` (final line)
All five `expectZero` checks pass; numeric coefficients (epsT `-0.7580…`, epsv `-1.0031…`, epsL `-1.8837…`, plus the three deltaPi coeffs) match the SymPy values to all shared digits.

**Output freshness:** confirmed. Script mtimes are 2026-05-28 15:54:51 (sympy `.py`) and 15:54:53 (mathematica `.wl`); the corresponding `.txt` outputs are 2026-05-28 16:10:14 (sympy) and 16:11:33 (mathematica) — both outputs post-date their scripts, so both transcripts were re-generated after Codex's fix. Neither log is truncated (both end cleanly with the carry-forward block / pass line).

## Material-change assessment

`material_change`: false.

The only edit is a printed-banner label string. No assertion, symbol definition, transport law, RHS form, numeric coefficient, or derived result changed. The five identity residuals and all numeric constants are byte-for-byte the same as pre-fix (apart from the banner line). No downstream unit can depend on a printed banner label, so units > 168 are not substantively affected by this change.

## Side observations (non-blocking)

- The SymPy transcript does not emit explicit `PASS:` lines (it relies on the `= 0` print plus the raise-on-failure contract in `expect_zero`), whereas the Mathematica transcript prints both `= 0` and `PASS:` lines. This asymmetry is purely cosmetic and pre-existing; it is not part of any finding and does not affect verification.
- The auditor noted (in the original report's independent-derivation section) that the `.wl` is a close structural port of the `.py` but diverges in its simplification path (free `s` with assumption `s == Sqrt[1+r^2]` vs. direct substitution), which the auditor judged corroborating rather than echoing. This was explicitly not raised as a finding and is unaffected by the banner fix. Noted for completeness only.

## Verdict justification

The single low-severity finding (F1, cosmetic banner mislabel) was applied exactly as the directive specified: a one-line `151` → `168` change in each of the two scripts, with no collateral edits per the git diff. Both engines were re-run; the regenerated transcripts now carry the correct `STAGE 168` banner, both exit 0, all five named identity checks pass with zero residual, and the numeric coefficients are unchanged and cross-engine consistent. Output mtimes post-date the script edits, confirming freshness. Every finding is `resolved`, exec logs pass, and there are no regressions. Verdict: `verified`.
