---
unit_id: 063
batch: III.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 063

The original report raised two low-severity findings (F1 stale `.txt` outputs, F2 stale SymPy self-labels). The directive folded the actionable self-label work into a single `F1` block (NUMBER-only `46`→`63` on `.py:3` and `.py:124`), with explicit SCOPE GUARD deferring `.py:31` banner padding and `.py:76` cross-reference to the dedicated numbering plan. Output refresh (report-F1) is handled by the orchestrator's exec re-run. I verify both report findings as closed under that decomposition.

## Per-finding outcomes

### F1 — stale numbering self-labels (SymPy)

**Classification:** resolved

**What changed:**
Two lines in `scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py`:
- `.py:3` docstring: `Stage 46 SymPy audit.` → `Stage 63 SymPy audit.`
- `.py:124` closing print: `print("\nAll Stage 46 symbolic checks passed.")` → `print("\nAll Stage 63 symbolic checks passed.")`

The git diff (`stage_063_diff.patch`) shows exactly these two hunks and nothing else — both strip-the-number identical (only `46`→`63`, all surrounding text preserved).

**Assessment:**
Correct and complete. I opened the post-edit file directly:
- `.py:3` reads `Stage 63 SymPy audit.` ✓
- `.py:124` reads `print("\nAll Stage 63 symbolic checks passed.")` ✓
- `.py:31` still reads `banner("STAGE 63 — PARENT-OVERLAP THRESHOLD THEOREM")` — unpadded `63`, UNCHANGED ✓ (deferred per SCOPE GUARD)
- `.py:76` still reads `# Insert Stage-44 threshold formulas` — UNCHANGED ✓ (cross-ref deferred per SCOPE GUARD)

Both edited lines are non-executable strings (docstring / print literal), so no assertion expression is touched and all residuals remain `= 0`. No collateral edit. The orchestrator's deferral of `.py:31`/`.py:76` is consistent with the standing in-loop Reading-2 policy (self-labels fixed in-loop; banner-padding and cross-references go to the dedicated numbering plan), so this is correctly scoped, not an incomplete fix.

### F2 — stale_output

**Classification:** resolved

**What changed:**
Both committed transcripts were re-generated post-fix. mtimes: `.../scripts/output/...sympy_audit.txt` = 2026-06-05 13:58 and `.../mathematica/output/...mathematica_audit.txt` = 2026-06-05 13:58, both newer than the SymPy script (13:43) and the `.wl` (2026-06-03 15:59, unedited). The SymPy transcript now opens `STAGE 63 — PARENT-OVERLAP THRESHOLD THEOREM` (line 3) and closes `All Stage 63 symbolic checks passed.` (line 21); the Mathematica transcript opens `STAGE 063 — PARENT THRESHOLDS` (line 3) and closes `Stage 063 Mathematica audit passed.` (line 32). No `Stage 46` / `STAGE 046` residual remains in either `.txt`.

**Assessment:**
Resolved. The stale `46`/`046` banners that triggered the freshness failure are gone; both files reflect the current scripts and every residual still prints `= 0` with passing footers. Note: the SymPy banner is `STAGE 63` (unpadded) rather than the `STAGE 063` the report's F2 "Verification" line anticipated — but that is the direct and correct consequence of the orchestrator deferring the `.py:31` padding to the dedicated plan, so the transcript faithfully mirrors the (intentionally still-unpadded) current script. Not a defect.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 63 — PARENT-OVERLAP THRESHOLD THEOREM` (banner now canonical-number)
- `coherence substitution in g_fail^2 = 0`, `C_suff^2/C_fail^2 - G_suff/G_fail = 0`, `G_micro - G_max*C^2 = 0`, `g_fail^2 from solve(G_micro=G_fail) matches hand-rearranged form = 0`, `G_max = G_micro at Cauchy saturation (O_sp^2 = N_ss N_pp) = 0`
- `All Stage 63 symbolic checks passed.`

**Mathematica:** exit=0. Notable lines:
- `STAGE 063 — PARENT THRESHOLDS`
- `PASS: coherence substitution in g_fail^2`, `PASS: C_suff^2/C_fail^2 - G_suff/G_fail`, `PASS: G_micro - G_max C^2`, `PASS: g_fail^2 from Reduce[gMicro==gFail, gphiSq>0] matches hand-rearranged form`, `PASS: G_max = G_micro at Cauchy saturation (oSP^2 = nSS nPP)`
- `Stage 063 Mathematica audit passed.`

The full identity battery (overlap↔coherence substitution, threshold ratio, factorization, four κ-insertions, both independent root re-derivations, Cauchy saturation) passes in both engines and the two engines agree symbol-for-symbol, exactly as the auditor inventoried (A1–A11 / B1–B11). The label edits left every assertion intact.

**Output freshness:** confirmed. Both `.txt` mtimes (13:58:32) are newer than the SymPy `.py` (13:43:21) and the `.wl` (2026-06-03; unedited, so its older mtime is expected and its transcript was still re-run as the 13:58 mtime shows). Banners/footers are canonical with no stale stage numbers.

## Material-change assessment

`material_change`: false. Both edits are NUMBER-only changes to a docstring string and a print literal; no symbolic expression, constant, or assertion was modified, and every residual remains `= 0`. No downstream unit can depend on a docstring/footer label. No upstream-stale cascade.

## Side observations (non-blocking)

- The report's F2 "Verification" line anticipated `STAGE 063` in the refreshed SymPy `.txt`, but the orchestrator deliberately deferred the `.py:31` padding, so the transcript correctly shows `STAGE 63`. The deferral is the standing policy; not blocking, and the Mathematica side already carries the padded `063` form, so the two engines differ only in zero-padding of the banner (cosmetic, tracked in the dedicated numbering plan).
- `.py:76` (`# Insert Stage-44 threshold formulas`) and the docstring `Stage-44`-implied cross-reference remain as the +17 cross-ref to upstream stage 061; correctly out of in-loop scope, deferred per SCOPE GUARD and the dedicated numbering plan.

## Verdict justification

The diff is exactly the two strip-the-number self-label edits the directive specified, applied correctly with no collateral change; lines 31 and 76 were verifiably left untouched per the SCOPE GUARD. Both transcripts were re-generated post-fix (mtimes 13:58 > script mtimes), now carry canonical banners/footers with no `46`/`046` residue, and both engines exit 0 with the entire non-tautological identity battery passing and agreeing symbol-for-symbol. Both original findings are resolved with no regression and no material change. Verdict: `verified`.
