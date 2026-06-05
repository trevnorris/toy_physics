---
unit_id: 071
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

# Verification — unit 071

## Per-finding outcomes

### F1 — stale_output (committed `.txt` banners)

**Classification:** resolved

**What changed:**
Both saved transcripts were re-generated post-fix. SymPy output
`scripts/output/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.txt:3`
now reads `STAGE 71 — CANONICAL TANH-WALL BRANCH` and `:31` `STAGE 71 THEOREM LEDGER`.
Mathematica output
`mathematica/output/moving_throat_pde_stage071_tanh_wall_branch_mathematica_audit.txt:3`
reads `STAGE 071 — CANONICAL TANH-WALL BRANCH`, `:46` `STAGE 071 THEOREM LEDGER`.
The prior stale `STAGE 54`/`STAGE 054` banners are gone.

**Assessment:**
Correct. The refreshed banners read the canonical numbers in each engine's
established format (SymPy 2-digit `STAGE 71`, Mathematica zero-padded `STAGE 071`),
exactly as F1's required change and the directive's banner-policy specify. Output
mtimes (13:58:32) are newer than both scripts (py 13:46:26; .wl Jun 3 15:59),
confirming a genuine post-fix re-run. All result lines (I_f=1/3, I_g=4/15,
I_g/I_f=4/5, T_X, K_X, J_1, W_wall, K_m, eta, both reductions) are unchanged from
the audit's recorded values — only the stage-label banner moved.

### F2 — stale numbering self-labels (SymPy docstring)

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py`
line 3: `moving_throat_pde_stage54_...` → `moving_throat_pde_stage071_...` (filename-style label, 3-digit).
Line 5: `SymPy audit for Stage 54:` → `SymPy audit for Stage 71:` (prose label, 2-digit).
Git diff (`stage_071_diff.patch`) shows exactly these two changed lines and nothing else.

**Assessment:**
Correct and NUMBER-only. Stripping the digits from old vs new yields byte-identical
lines (the filename label keeps `stage`+3-digit to match the on-disk filename; the
prose label keeps 2-digit `Stage`+space). No collateral edit: the two banner strings
at py lines 24/91 remain `STAGE 71` — correctly LEFT UNPADDED (correct number,
reserved for the dedicated script/output-band plan). No `.wl` edit, no math, no
cross-reference, no variable name touched. Source now exits 0 (sympy log).

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- L8 banner `STAGE 71 — CANONICAL TANH-WALL BRANCH` (refreshed)
- `I_f - 1/3 = 0`, `I_g - 4/15 = 0`, `I_g/I_f - 4/5 = 0` (L16,18,20)
- `kappa_reduced - [4 chi_s^2 + 4/5 Lambda_ell^2] = 0` (L31); `W_wall_reduced - Upsilon_w Lambda_ell^2 = 0` (L33)
- `# exit_code: 0` (L44)

**Mathematica:** exit=0. Notable lines:
- L8 banner `STAGE 071 — CANONICAL TANH-WALL BRANCH`
- `PASS: I_f - 1/3`, `PASS: I_g - 4/15`, `PASS: I_g/I_f - 4/5`
- `PASS: T_X exact formula`, `PASS: K_X exact formula`, `PASS: J_1 exact formula`
- `PASS: kappa reduced law`, `PASS: W_wall reduced law`
- `# exit_code: 0` (L59)

**Output freshness:** confirmed. Both `.txt` outputs have mtime 2026-06-05 13:58:32,
newer than the SymPy script (13:46:26) and the `.wl` (Jun 3 15:59). Banners reflect
the post-fix canonical numbering in both engines.

## Material-change assessment

`material_change`: false.

All edits are docstring self-labels (py:3, py:5) and refreshed transcript banners.
No derived result, assertion, symbol, or numeric value changed; both engines emit
the identical math they did pre-fix. No downstream unit depends on a label string,
so no unit > 071 is affected on math grounds.

## Side observations (non-blocking)

- The two SymPy runtime banners (py:24, py:91) still print `STAGE 71` (2-digit)
  while the Mathematica banners are zero-padded `STAGE 071`. This banner-padding
  divergence is intentionally deferred to the dedicated script/output-band numbering
  plan, not in-scope here; flagging only for continuity.

## Verdict justification

Both findings are fully resolved. F2's git diff is precisely the two strip-the-number
self-label edits the directive specified — NUMBER-only, format-preserving, no collateral
change, and the `STAGE 71` ledger banner correctly left unpadded. F1's transcripts were
genuinely re-run post-fix (output mtimes newer than scripts) and now carry the canonical
`STAGE 71`/`STAGE 071` banners with all result lines unchanged. Both engines exit 0 with
every PASS/residual check holding. No regressions in the diff or logs; material_change
false. Verdict: verified.
