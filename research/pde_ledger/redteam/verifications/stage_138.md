---
unit_id: 138
batch: IV.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 138

## Per-finding outcomes

### F1 — paper_misalignment (banner mismatch)

**Classification:** resolved

**What changed:**
This batch (IV.4) was orchestrator-direct: the user resolved the paper_misalignment via direction (a) from the directive's "Resolve before fix_loop" block (banner is a typo). The Mathematica banner at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage138_normalized_mouth_gain_family_mathematica_audit.wl:26` now reads `banner["STAGE 138 — NORMALIZED MOUTH-GAIN FAMILY AND COMPENSATION RATIO"];` (was `STAGE 121`). No other content changed in the `.wl` — surrounding context (L20-34: `expectZero` helper, `Clear`, `$Assumptions`, `mS`, `mQ`) is intact and untouched. The exec log at `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage138_normalized_mouth_gain_family_mathematica_audit.txt` was regenerated post-fix: L3 reads `STAGE 138 — NORMALIZED MOUTH-GAIN FAMILY AND COMPENSATION RATIO`, and all three `PASS:` lines (R_q exact formula, R_q on + branch, R_q on - branch) remain. The paper card (`paper/stages/stage_138.tex:1-2`) and notes header (`notes/stages/moving_throat_pde_stage138_normalized_mouth_gain_family.md:1`) confirm Stage 138 is the canonical owner — alignment is now consistent across paper, notes, filename, and script banner.

**Assessment:**
The edit is exactly what the directive's option (a) prescribed: single-token in-place banner fix, no math touched. The required "Verification" criterion from the original report (banner transcript says STAGE 138, all three PASS lines present, EXIT_CODE 0 preserved) is met. No collateral edits in the `.wl`. No new findings introduced.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- L11 `R_q = (g_c - r)**2/(r**2 + 1)`
- L12-13 `R_q on + branch = 1/4` / `R_q on - branch = 1/4`
- L17 `EXIT_CODE: 0`

**Mathematica:** exit=0. Notable lines:
- L3 `STAGE 138 — NORMALIZED MOUTH-GAIN FAMILY AND COMPENSATION RATIO` (banner now correct)
- L9 `PASS: R_q exact formula`
- L13 `PASS: R_q on + branch - 1/4`
- L15 `PASS: R_q on - branch - 1/4`
- L17 `Normalized mouth-gain family verified.`

**Output freshness:** Mathematica `.txt` mtime (1779925678) is newer than `.wl` mtime (1779924055) — log was regenerated after the banner fix. SymPy `.txt` mtime (1778525149) is newer than `.py` mtime (1775100991); SymPy was not touched by this batch but its log is still fresh relative to its script. Both engines agree on every quoted residual and final ratio (R_q = (gCore-rHat)^2/(1+rHat^2), both branches = 1/4), preserving the engine cross-check from the original audit.

## Material-change assessment

`material_change`: false.

The change is a banner string literal only. No derived quantities were modified; SymPy was not touched; both engines still produce the same residuals and the same `PASS:` outcomes as before the fix. Downstream units do not depend on this banner text.

## Side observations (non-blocking)

None. The math content of the audit (closed-form R_q, ±-branch checks, definitional Σ_0 substitution) is unchanged and remains aligned with the paper card and notes deliverables.

## Verdict justification

The sole low-severity finding F1 (banner mislabel `STAGE 121` → `STAGE 138`) is resolved: the `.wl` banner literal at L26 has been corrected, the Mathematica exec log was regenerated post-fix and now displays the correct stage banner with all three `PASS:` assertions and exit 0, and no math was altered. SymPy continues to exit 0 with matching residuals. No regressions, no new findings, material_change is false.
