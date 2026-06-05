---
unit_id: 041
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

# Verification — unit 041

## Per-finding outcomes

### F1 — stale_output (self-label) [label-only]

**Classification:** resolved

**What changed:**
In `scripts/moving_throat_pde_stage041_rank2_support_sympy_audit.py`, five stale
pre-renumber self-labels were canonicalized 24→41:
- py:3 docstring header `Moving-throat PDE — Stage 24 SymPy audit.` → `Stage 41 SymPy audit.`
- py:53 subbanner `24.1` → `41.1`
- py:84 subbanner `24.2` → `41.2`
- py:96 subbanner `24.3` → `41.3`
- py:104 subbanner `24.4` → `41.4`

**Assessment:**
Each diff hunk changes ONLY the leading numeric stage token; the trailing descriptive
string of every subbanner is byte-for-byte preserved (`— Exact rank-2 determinant…`,
`— Exact monotonicity…`, etc.). No equation, value, variable name, `expect_zero` call,
matrix construction, or assertion is touched. The already-canonical `STAGE 41` banners
(py:51, py:142) were left as-is, and the `Stage-23` CROSS-refs (py:13, py:151 — pointing
at the upstream one-direction geometry, now stage 040) were correctly LEFT untouched.
Strictly label-only and matches the directive exactly.

### F2 — stale_output (self-label, missed comment) [label-only]

**Classification:** resolved

**What changed:**
py:107 inline comment `# Derive n_src by substituting into the general n_expected from
section 24.1,` → `… from section 41.1,`.

**Assessment:**
Confirmed in the live file at line 107 (`from section 41.1`). This is a comment, so it
never reaches the transcript, but it is now internally consistent with subbanner `41.1`.
Single numeric-token change; the variable names `n_src`/`n_expected` referenced in the
same comment were left intact, as directed. The captured `stage_041_diff.patch` is
truncated mid-hunk and shows only the `-` (old) side of this line, but the current file
state confirms the `+` side landed correctly — a grep for any residual `24` token in the
whole `.py` returns ZERO hits, so no self-label was missed and no collateral `24` was
introduced.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 41 — RANK-2 SUPPORT COMPLETION` (now matches py:51) and subbanners `41.1`,
  `41.2`, `41.3`, `41.4` (txt:12,30,42,53), ledger `STAGE 41 THEOREM LEDGER` (txt:82).
- All eight residuals present and zero: `determinant decomposition = 0`,
  `n_req - expected = 0`, `dn/dm - expected = 0`, `tracking collapse = 0`,
  `source-tied formula = 0`, `source-tied dn/dm = 0`, plus the two printed thresholds.

**Mathematica:** exit=0. Notable lines:
- `STAGE 041 — RANK-2 SUPPORT COMPLETION` (txt:8) and `Stage 041 Mathematica audit
  passed.` (txt:47). The `.wl` was untouched (its labels were already canonical), so the
  Mathematica banner is unchanged and consistent.
- `PASS: determinant decomposition`, `PASS: n_req - expected`, `PASS: dn/dm - expected`,
  `PASS: tracking collapse`, `PASS: source-tied formula`, `PASS: source-tied dn/dm` — all
  six engine checks PASS; all `= 0` residuals preserved.

**Output freshness:** Both transcripts were regenerated post-fix. SymPy log dated
2026-06-05T08:58:41 and Mathematica 2026-06-05T08:10:07 (both well after the
2026-06-05 edit), each ending `# exit_code: 0`. The SymPy transcript now reads `STAGE 41`
/ `41.1–41.4` (resolving the prior `STAGE 24` staleness flagged in the report), and the
prior banner/subbanner internal disagreement is gone.

## Material-change assessment

`material_change`: false.

The diff touches only stage-number label tokens in one comment, one docstring line, and
four subbanner strings of the SymPy script. No derived symbolic result, threshold, matrix,
solve/diff operation, or assertion changed — every `expect_zero`/`expectZero` residual is
still literal `0`/`PASS` in both engines, and the printed closed forms (`D_sel`, `n_req`,
`dn/dm`, tracking collapse, `n_req^(src)`, thresholds, source-tied derivative) are
byte-identical to the pre-fix content. No downstream unit can depend on a banner string,
so nothing > 041 is affected by the math.

## Side observations (non-blocking)

- The captured `stage_041_diff.patch` is truncated (48 lines) and ends mid-hunk after the
  `-` side of the line-107 comment, so the F2 `+` line is not visible in the patch file
  itself. This is a capture artifact only; the live file and a zero-hit `24` grep confirm
  F2 landed correctly. No verification impact.

## Verdict justification

Both findings are `resolved`: the diff is strictly label-only (only stage-number tokens
on the docstring, four subbanners, and one comment), surrounding text is preserved
verbatim, the `Stage-23` cross-refs and the already-canonical `STAGE 41` banners are
untouched, and a whole-file grep finds no residual `24` token. Both engines re-ran to
exit 0 with every prior `= 0`/`PASS` residual intact, and the refreshed SymPy transcript
now reads `STAGE 41` / `41.1–41.4` while the untouched Mathematica transcript reads
`STAGE 041`. No math changed. Verdict: verified, material_change false.
