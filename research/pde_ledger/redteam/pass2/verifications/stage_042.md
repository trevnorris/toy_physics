---
unit_id: 042
batch: III.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 042

## Per-finding outcomes

### F1 — stale_output (self-label) [label-only]

**Classification:** resolved

**What changed:**
Codex edited only the SymPy source `scripts/moving_throat_pde_stage042_rank2_selected_mode_sympy_audit.py`. Per the captured diff (`exec_logs/stage_042_diff.patch`) and the current file state:
- line 3: docstring header `Moving-throat PDE — Stage 25 SymPy audit.` → `Stage 42` (file:3).
- line 53: subbanner `25.1` → `42.1` (file:53).
- line 75: subbanner `25.2` → `42.2` (file:75).
- line 103: subbanner `25.3` → `42.3` (file:103).
- line 114: subbanner `25.4` → `42.4` (file:114).
- line 132: subbanner `25.5` → `42.5` (file:132).
- line 140: subbanner `25.6` → `42.6` (file:140).

Exactly 7 changed lines; the diff shows only the leading numeric stage token replaced on each, with all surrounding text preserved verbatim.

**Assessment:**
Strictly label-only and matches the directive token-for-token. Confirmed:
- **No equation/value/assertion change.** Every changed line is either the docstring header or a `subbanner("...")` call (a print-string); no `expect_zero`, `sp.simplify`, symbol, or RHS `_expected` form is in any changed hunk. The diff's context lines (`n_req`, `D_qr`, `F_track`, etc.) are unchanged.
- **`F_stage23` NOT renamed.** Lines 106 and 112 still read `F_stage23` (CODE identifier left intact, as required).
- **Cross-ref tails preserved.** Line 53 retains `... after inserting the Stage-24 support loading`; line 103 retains `... collapse back to Stage 23`. The docstring cross-refs (`Stage-24` line 7, `Stage-23` line 12) and ledger cross-ref (`Stage-23` line 184) are untouched. These are deferred to the dedicated numbering pass per the directive.
- **Banners already canonical, left alone.** `STAGE 42 — ...` (line 51) and `STAGE 42 THEOREM LEDGER` (line 180) unchanged.
- **`.wl` untouched.** Mathematica source mtime is 2026-06-03 15:59 (unchanged from pre-fix); no `.wl` hunk in the diff. Correct — its labels are already canonical.

The finding is fully resolved at the level the directive scoped (UNAMBIGUOUS self-labels). Refreshed outputs confirm the new banner/subbanners are emitted with all residuals 0.

## Exec log assessment

**SymPy:** exit=0. Notable lines from `scripts/output/...sympy_audit.txt`:
- banner `STAGE 42 — SELECTED-MODE NORMALIZATION UNDER RANK-2 SUPPORT COMPLETION` and subbanners `42.1`–`42.6` all present.
- `row1 - expected = 0`, `row2 - expected = 0`, `Z_overlap - expected = 0`, `S_overlap - expected = 0`, `F_general - expected = 0`, `tracking collapse = 0`, `source-tied specialization = 0`, `F_src(R_U=1) - F_flat = 0`, `H_n^(src) - expected = 0`, `H_F^(src) - expected = 0`, `n_src - linear expansion = 0`, `F_src/F_flat - linear expansion = 0`.
- `STAGE 42 THEOREM LEDGER`; `# exit_code: 0`.

**Mathematica:** exit=0. Notable lines from `mathematica/output/...mathematica_audit.txt`:
- banner `STAGE 042 — RANK-2 SELECTED-MODE NORMALIZATION` (canonical `042`, unchanged).
- `PASS: row1 - expected`, `PASS: row2 - expected`, `PASS: Z_overlap - expected`, `PASS: S_overlap - expected`, `PASS: F_general - expected`, `PASS: tracking collapse`, `PASS: source-tied specialization`, `PASS: F_src(R_U=1) - F_flat`, `PASS: H_n^(src) - expected`, `PASS: H_F^(src) - expected`, `PASS: n_src - linear expansion`, `PASS: F_src/F_flat - linear expansion`.
- `Stage 042 Mathematica audit passed.`; `# exit_code: 0`.

Every prior PASS / `= 0` residual is retained; the SymPy banner now reads `STAGE 42` and subbanners `42.1`–`42.6` (the only label changes), and the `.wl`-side `STAGE 042` banner is correctly unchanged.

**Output freshness:** confirmed. Both `.txt` mtimes are 2026-06-05 08:11, newer than the `.py` (2026-06-05 07:44) and the `.wl` (2026-06-03 15:59). The `.wl` source mtime is unchanged from the original, corroborating that the Mathematica script was not edited and only re-run.

## Material-change assessment

`material_change`: false.

No derived result changed. The edits replace only stage-number tokens in a docstring and six print-only subbanner strings. All 12 residual-to-zero assertions (and their Mathematica mirrors) are byte-identical in math and still pass. The closed forms emitted (`e1/e0`, `D_{q,r}`, overlaps, `F_(q,r,t)`, `F_track`, `F_src`, `F_flat`, `H_n^(src)`, `H_F^(src)`, series) are unchanged. No downstream unit can depend on a label string, so nothing > 042 is invalidated.

## Side observations (non-blocking)

- The committed `.wl`/its output still carry the pre-renumber CROSS-refs `Stage-24` / `Stage 23` and the bare `1.`–`6.` Mathematica section indices (vs SymPy's `42.k`). These remain owned by the deferred SCRIPT/OUTPUT-band numbering pass (`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md`) per the directive's explicit DO-NOT-TOUCH list — out of scope here and correctly left.
- SymPy docstring/ledger cross-refs `Stage-24` (line 7) / `Stage-23` (lines 12, 184) likewise remain for the dedicated pass.

## Verdict justification

The single low-severity `stale_output` finding is resolved exactly as scoped: the SymPy self-labels (docstring header + six subbanner indices) were canonicalized `25/25.k → 42/42.k`, the diff is strictly label-only with no equation/value/assertion/variable-name change, `F_stage23` and all cross-ref tails were correctly left untouched, the `.wl` was not modified, and both engines re-ran to exit 0 with every prior PASS/`= 0` residual intact and refreshed outputs newer than their scripts. No regressions, no partials. Verdict: verified; material_change: false.
