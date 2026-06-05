---
unit_id: 068
batch: III.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T14:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 068

The original report raised 2 findings, both low-severity and non-math:
- **F1** — stale committed `.txt` outputs (both engines: old `STAGE 51`/`STAGE 051` banners, obsolete `P_res*C_res^2 - 1 = 0` line, old `(A)/(B)` band labels).
- **F2** — stale `stage51` self-label in the SymPy docstring at `py:4` (plus a note that the `STAGE 68` banner at `py:37` is non-3-digit).

The directive folded these into a single applied finding (the `py:4` docstring self-label), with the `.txt` refresh handled by the orchestrator exec-refresh and the `STAGE 68` banner deliberately LEFT UNPADDED (banner padding deferred to the dedicated SCRIPT/OUTPUT-band plan). I verify both report findings below.

## Per-finding outcomes

### F1 — stale_output (committed `.txt` transcripts)

**Classification:** resolved

**What changed:**
No script edit; the two committed transcripts were re-generated post-fix. Confirmed via mtimes: both `.txt` files are 2026-06-05 13:58:32, newer than the `.py` (13:46:24) and `.wl` (06-03 15:59:11) sources.

**Assessment:**
The refreshed SymPy `.txt` (exec log) now reads `STAGE 68 — RESONANCE-CORRECTED THRESHOLDS` at the banner, contains the `P_res numeric residual = 5.69...E-16` line (line 13), and uses `(C-form)`/`(P-form)` band labels (lines 22-25) — the obsolete `P_res*C_res^2 - 1 = 0` line and `(A)/(B)` labels are gone. The Mathematica `.txt` reads `STAGE 068 — ...` (line 8), uses `(C-form)`/`(P-form)` labels (lines 33-36), and all `PASS:` lines are present. Both exit 0. The freshly-gained `P_res numeric residual` line and `(C-form)/(P-form)` labels are legitimate (prior June-3 script update), not collateral from this edit. This matches the report's F1 verification criteria exactly.

Note: the SymPy `.txt` banner reads `STAGE 68` (not `STAGE 088/068`) — this is the correct UNPADDED canonical number per the scope guard; only the engine cross-banner styling differs (Mathematica uses `068`), and that styling padding is explicitly deferred. The report's F1 verification text said the refreshed SymPy banner should read `STAGE 68`, which is what the log shows. Consistent.

### F2 — stale_output (self-label in SymPy docstring)

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py:4` docstring filename-label changed `moving_throat_pde_stage51_...` → `moving_throat_pde_stage068_...`. The git diff (`stage_068_diff.patch`) shows exactly one changed line; the two strings are byte-identical except `stage51`→`stage068`.

**Assessment:**
Correct strip-the-number filename-label edit. `git diff --stat` confirms 1 file changed, 1 insertion / 1 deletion — no collateral edits. I read the current `py:4` (now `stage068`, matching the on-disk filename, 3-digit as required for a filename-style self-label) and `py:37` (banner still `STAGE 68`, NOT padded) — exactly as the directive's scope guard mandated. The directive's required change listed padding the `py:37` banner to `STAGE 068`, but the orchestrator scope guard correctly overrode that (banner padding deferred to the dedicated plan); Codex applied only the `py:4` edit, which is the right call. SymPy exits 0 post-edit. Non-tautological: this is a pure label edit, touches no assertion.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 68 — RESONANCE-CORRECTED THRESHOLDS` (correct unpadded banner)
- `P_res numeric residual = 5.6958391724936524581E-16` (new line present)
- `success band C-form vs P-form (under Pres = 1/Cres^2) = 0` / `failure ... = 0` (band identities hold, new labels)

**Mathematica:** exit=0. Notable lines:
- `STAGE 068 — RESONANCE-CORRECTED THRESHOLDS`
- `PASS: W_res - C2 * W_wall (from gain decomposition)`
- `PASS: success band C-form vs P-form (under Pres = 1/Cres^2)` / `PASS: failure ...`
- All `PASS:` lines present; no failures or warnings.

**Output freshness:** confirmed. Both `.txt` mtimes (2026-06-05 13:58:32) are newer than the `.py` (13:46:24) and `.wl` (2026-06-03) source files.

## Material-change assessment

`material_change`: false.

No derived result changed. The `.py` edit is a docstring filename label; the `.wl` is untouched; the `.txt` refresh only updates banners/labels and surfaces the (pre-existing, June-3) `P_res numeric residual` line. All numeric results — `P_res = 1.005612487760576`, `C_res^2 = 0.994418836451529`, the threshold formulas, and the band widths — are unchanged. No downstream unit is affected.

## Side observations (non-blocking)

- The directive front-matter says `findings_count: 1` whereas the auditor report lists `findings_count: 2`; this is the expected F1/F2 fold (the .txt refresh is orchestrator-handled, not a Codex code edit), not a defect.
- The directive's `## F1` "required change" still includes the `py:37` `STAGE 68`→`STAGE 088`-style padding line, which the orchestrator scope guard correctly suppressed. The applied result follows the scope guard, not the stale required-change bullet. No action needed; just noting the directive text is internally slightly inconsistent with its own scope guard.

## Verdict justification

Both report findings are resolved. The source diff is exactly the single strip-the-number filename-label edit (`stage51`→`stage068` at `py:4`), the `STAGE 68` banner at `py:37` was correctly left unpadded per the scope guard, the `.wl` was not touched, and both committed `.txt` outputs were refreshed with canonical content (new banners, `P_res numeric residual` line, `(C-form)/(P-form)` labels) and are newer than their sources. Both engines exit 0 with all `PASS:` lines present. No regressions, no tautology introduced (label-only edit), `material_change` false.
