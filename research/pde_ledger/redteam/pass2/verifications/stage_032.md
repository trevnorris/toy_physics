---
unit_id: 032
batch: II.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-04T23:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 032

## Per-finding outcomes

### F1 — stale_output (residual self-labels + stale transcripts)

**Classification:** resolved

**What changed:**
The diff (`stage_032_diff.patch`) is exactly three label-only line edits, nothing else:
- `scripts/...stage032..._sympy_audit.py:3` — docstring `Moving-throat PDE Stage 15 SymPy audit.` → `Moving-throat PDE Stage 32 SymPy audit.`
- `scripts/...stage032..._sympy_audit.py:217` — `print("All Stage 15 checks passed.")` → `print("All Stage 32 checks passed.")`
- `mathematica/...stage032..._mathematica_audit.wl:198` — `Print["All Stage 15 checks passed."];` → `Print["All Stage 32 checks passed."];`

The directive's `## Applied: F1` block reports these two files changed, summary "replaced stale Stage 15 self-labels … with Stage 32," deviation: none — consistent with the diff.

**Assessment:**
The edit matches the directive's "required change" precisely. Every changed line is a pure comment/print string; the only token altered is the stage number `15`→`32`. There is ZERO change to any equation, value, variable, assertion, or `expect_zero`/`expectZero` target — the diff hunks touch only the docstring line and the trailing summary `print`/`Print`, both outside any computation. No collateral edits. A grep over both live scripts for any `stage 15` / `stage-15` / `stage_15` variant returns empty, confirming no residual stale labels remain. The fix is label-only and complete; no assertion was weakened or made tautological (no assertion was touched at all).

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- Header banners now canonical: `STAGE 32.1 …`, `STAGE 32.2 …`, … `STAGE 32.5 …` (output lines 8, 28, 42, 108, 117).
- Every residual prints `= 0`: e.g. `sigma/kappa0^2 - 11/9 = 0`, `mhat_-^2(alpha=0) - 1 = 0`, `limit_{alpha->oo} mhat_-^2 - 11/9 = 0`, `Nprod(alpha=0) - beta0 * kappa0^2 / A = 0`, `limit_{alpha->oo} Nprod_nat = 0` (lines 25, 113, 114, 119, 120).
- Closing line: `All Stage 32 checks passed.` (line 121), `# exit_code: 0` (line 122). No FAIL.

**Mathematica:** exit=0. Notable lines:
- Header banners canonical `STAGE 32.1 … STAGE 32.5`; every check emits `PASS:` (e.g. `PASS: sigma/kappa0^2 - 11/9`, `PASS: Sigma - [Xi I + alpha vv^T]`, `PASS: s_check - s_minus_nat …`, `PASS: limit_{alpha->oo} Nprod_nat`).
- Closing line: `All Stage 32 checks passed.` (line 93), `# exit_code: 0` (line 94). No FAIL.
- Two `Limit::alimv` warnings (lines 78, 90) are the benign limit-variable-assumption warnings already noted by the auditor; both associated limits still evaluate and PASS (11/9 and 0).

The two engines' central `mhat_-^2` printouts remain algebraically identical (sympy out 110 vs wl out 70), unchanged by this label-only fix.

**Output freshness:** Confirmed regenerated post-fix. Both committed `.txt` outputs have mtime 2026-06-04 23:13:28, newer than both scripts (mtime 22:51:59). The transcripts now read the canonical `STAGE 32.x` banners and the corrected closing line `All Stage 32 checks passed.`, resolving the stale-transcript half of the finding.

## Material-change assessment

`material_change`: false.

The entire diff is cosmetic stage-label string edits (docstring + two closing summary prints). No derived result, residual, closed form, or assertion target changed. The `mhat_-^2` closed form, the `1 ≤ \widehat m_-^2 < 11/9` bound endpoints, the Schur decomposition, and the eliminated product `N_-` are all byte-for-byte unchanged in the refreshed transcripts. No downstream unit can be affected by a comment/banner relabel.

## Side observations (non-blocking)

None beyond what the auditor already noted (the unused `.wl` symbol declarations `gConst, cs, radius, cSpeed, mhat`, benign dead code, untouched by this fix). Not blocking.

## Verdict justification

The single finding was a label-only `stale_output`: stale `Stage 15` self-labels in the SymPy docstring/closing print and the `.wl` closing Print, plus pre-fix-dated transcripts. The diff is exactly the three named label edits with zero math/equation/assertion/value change; both engines re-ran to exit 0 with all residuals `0`/`PASS` and no FAIL; the refreshed transcripts read canonical `STAGE 32.x` banners and `All Stage 32 checks passed.`; outputs are confirmed newer than scripts; and no residual `Stage 15` token survives in either live script. Finding resolved, no regression, cosmetic only → `verified`, `material_change: false`.
