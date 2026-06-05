---
unit_id: 064
batch: III.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T20:05:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 064

## Per-finding outcomes

### F1 — stale_output

**Classification:** resolved

**What changed:**
Both committed transcripts were re-generated post-fix by the orchestrator's exec re-run, not by a script-logic edit (F1 explicitly required no logic change).

- SymPy `.txt` (`scripts/output/...sympy_audit.txt`): banner now `STAGE 64` (L3, was `STAGE 47`); contains the `GENERAL-H EQUILIBRIUM SOFTENING CHECK` section (L40) with `Theta (general, two-point)` (L42) and `Lambda (general, two-point)` lines; the obsolete `Lambda^2/Theta (closure form)` print is gone. Closing `STAGE 64 AUDIT PASSED` (L55).
- Mathematica `.txt` (`mathematica/output/...mathematica_audit.txt`): banner now `STAGE 064` (L3, was `STAGE 047`); contains the `CONTINUOUS CAUCHY BOUND CHECK` section (L37) with `continuous pair Cauchy gap` (L39) and `expected continuous gap` (L43); the obsolete discrete `((h1−h2)^2 w1 w2)/(h1^2 h2^2)` gap is gone.

**Assessment:**
Correct and complete. Each of F1's enumerated verification checks is satisfied (new sections present, obsolete prints absent, canonical banners). The new `GENERAL-H EQUILIBRIUM SOFTENING CHECK` and `CONTINUOUS CAUCHY BOUND CHECK` sections are not introduced by this fix — they reflect prior June-3 script revisions that the ~12-day-stale committed transcripts had never captured, exactly as the report and orchestrator note describe. The SymPy banner remains 2-digit `STAGE 64` because the L49/L176 banner padding was correctly deferred (per the scope guard); this is consistent between the script source and its refreshed `.txt`, so it is not a defect. mtimes confirm the refresh (below).

### F2 — paper_misalignment (stale self-label, numbering; routed under in-loop Reading-2 as a label-only self-label fix)

**Classification:** resolved

**What changed:**
The directive narrowed F2's full label list to the single unambiguous SELF-label and explicitly deferred the cross-refs and banner padding. Codex applied exactly that.

- `scripts/...sympy_audit.py:3` docstring: `Moving-Throat PDE — Stage 47 SymPy audit` → `Moving-Throat PDE — Stage 64 SymPy audit`. NUMBER-only (`47`→`64`), format preserved. This is the sole hunk in `stage_064_diff.patch` and the only changed line in `git diff --stat`.

Verified LEFT UNTOUCHED (deferred, per scope guard):
- `.py:25` `reproducing the Stage-45/46 best-alignment formulas.` — unchanged.
- `.py:122` assertion label `"matched-layer gain vs Stage-45 best-alignment formula"` — unchanged.
- `.py:180` `In the matched-layer limit this reproduces the Stage-45/46 best-alignment formulas.` — unchanged.
- `.py:49` / `.py:176` banners `STAGE 64` / `STAGE 64 AUDIT PASSED` — unchanged (no 3-digit padding).
- `.wl` — no edit at all (confirmed `git diff --stat`: only the `.py` changed, 1 insertion / 1 deletion; the `.wl:104` cross-ref is therefore untouched).

**Assessment:**
Correct and minimal. The change matches the directive's required edit verbatim, is a comment/string-only NUMBER edit (no expression, symbol, or assertion change), and the deferral scope was honored exactly — no collateral edits. `47 + 17 = 64` is the documented `+17` pre-renumber drift, so the canonicalization is right. Cross-refs to stages 062/063 (the `Stage-45/46` tokens, one inside an `expect_zero` label) are correctly left for the dedicated numbering plan. Label-only; no math touched.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 64 — PARENT EQUILIBRIUM SOURCE/SUPPORT ALIGNMENT` (banner)
- `two-point Cauchy gap identity = 0`
- `general equilibrium softening equals g_phi^2 I_1 = 0`
- `STAGE 64 AUDIT PASSED` — all `expect_zero` residuals print `= 0`.

**Mathematica:** exit=0. Notable lines:
- `PASS: closure law chi_sigma = g_phi chi_phi/H`
- `PASS: continuous Cauchy bound C^2 <= 1` (with `N_pp I2 - I1^2 = (L^2*Pi^2)/(64*hw^2)` matching `expected continuous gap`)
- `PASS: effective support softening`
- `Stage 064 Mathematica audit passed.` — every `expectZero` reports PASS.

No regressions: the diff is a single docstring line; both engines still pass all assertions. The substantive checks (closure, overlap identities, coherence form, matched-layer reductions, both Cauchy routes, general-H softening, eliminated-source softening) remain intact and non-tautological as the auditor established.

**Output freshness:** confirmed. mtimes:
- SymPy `.py` = 2026-06-05 13:43:33; SymPy `.txt` = 2026-06-05 13:58:32 (`.txt` newer ✓).
- Mathematica `.wl` = 2026-06-03 15:59:11; Mathematica `.txt` = 2026-06-05 13:58:32 (`.txt` newer ✓).
Both committed transcripts post-date their scripts; F1's mtime criterion is met.

## Material-change assessment

`material_change`: false.

Both edits are non-mathematical: F2 is a docstring NUMBER-only self-label correction, and F1 is an output transcript refresh (no script-logic change). No derived result, constant, closed form, or assertion was altered. All emitted deliverable values are unchanged from the auditor's reconciliation (0 misaligned). No downstream unit can depend on a docstring string or a transcript banner, so no unit > 064 is substantively affected.

## Side observations (non-blocking)

- The SymPy committed `.txt` still carries the 2-digit `STAGE 64` banner and the `Stage-45/46` cross-ref in its closing summary (line 65 of the log), and the `.wl` still emits `Stage-45 best-alignment formula` (mma log L38). These are the deferred cross-ref / banner-padding items, correctly out of scope for this unit and routed to the dedicated numbering plan. Flagged only for tracking; not a basis to fail verification.

## Verdict justification

Both findings are resolved. F2's single authorized SELF-label edit (`.py:3` `Stage 47`→`Stage 64`) is exactly the one hunk in the diff, NUMBER-only with format preserved, and every deferred cross-ref (`.py:25/122/180`, `.wl:104`) plus the `STAGE 64` banners (`.py:49/176`) were verifiably left untouched — no `.wl` change and no collateral edits. F1's transcripts were refreshed post-fix: canonical banners restored, the new `GENERAL-H` / `CONTINUOUS CAUCHY` sections present, obsolete prints gone, and both `.txt` mtimes newer than their scripts. Both engines exit 0 with all `expect_zero`/`expectZero` residuals at 0 / PASS. The fix is label-only; `material_change: false`.
