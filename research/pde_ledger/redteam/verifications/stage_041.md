---
unit_id: 041
batch: III.1
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 041

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage041_rank2_support_sympy_audit.py:104-116`:
  the `q2_src/r2_src/qr_src/qr_diff_sq` placeholder block is replaced with
  `t = sp.symbols("t", real=True)`, then
  `n_src = sp.simplify(n_expected.subs({q: t*R_U, r: t}))` and
  `n_src = sp.simplify(sp.expand(n_src).subs(t**2, lam0))`. The
  `n_src_expected` block (now lines 117–120) and the
  `expect_zero("source-tied formula", n_src - n_src_expected)` call
  (line 121) are unchanged, exactly as the directive specified.
- `mathematica/moving_throat_pde_stage041_rank2_support_mathematica_audit.wl:35-38`:
  `Clear[...]` and `$Assumptions` both extended to include `t`.
- `mathematica/.../stage041_..._mathematica_audit.wl:88-97`: the `q2Src/r2Src/qrSrc/qrDiffSq`
  block is replaced with `nSrc = FullSimplify[nExpected /. {q -> t rU, r -> t}, Assumptions -> $Assumptions]`
  followed by `nSrc = FullSimplify[Expand[nSrc] /. t^2 -> lambda0, Assumptions -> $Assumptions]`.
  `nSrcExpected`, `regThreshold`, `numZeroThreshold`, `dndmSrc`, and both `expectZero`
  calls are unchanged.

**Assessment:**
The edit matches the directive verbatim in both engines. The substantive
content change is exactly the one the auditor required: `n_src` is now
constructed by physical substitution into the general `n_expected` from
section 24.1, not as a fresh literal rational function. As a result,
`expect_zero("source-tied formula", n_src - n_src_expected)` and
`expectZero["source-tied formula", nSrc - nSrcExpected]` are no longer
character-for-character identities; they now verify that substituting
`q = t R_U, r = t, t^2 = lam0` into the general support-loading formula
reduces to the hand-written source-tied form. Both exec transcripts show
the residual as `0` (sympy) and `PASS` (Mathematica), so the substitution
does in fact close the identity. No collateral edits in either file; the
`dn_dm_src` derivative check (which the auditor explicitly said to leave
intact) is untouched.

The new sympy assertion is genuinely non-tautological: `n_expected` is
the `sp.simplify(...)` of the general formula derived independently of
the hand-written `n_src_expected`, so the equality after `subs` plus
`t**2 -> lam0` is a real algebraic identity. Same for the Mathematica
side: `nExpected` is the general formula and `nSrcExpected` is a separate
hand-written specialization. The displayed `n_req^(src)` in both outputs
shows the correctly simplified rational function (sympy: numerator
`2 R_U^2 m xi + 9 delta m - 9 delta xi + 9 m xi - 9 xi^2`, denominator
`2 R_U^2 m - 4 R_U m - 9 delta + 2 m - 11 xi`; Mathematica:
`(9 delta(-m+xi) + xi(-(m(9+2 rU^2))+9 xi))/(9 delta - 2 m (-1+rU)^2 + 11 xi)`),
which agrees up to overall sign and matches the hand-written form.

## Exec log assessment

**SymPy:** exit=0 (no FAIL/Traceback/AssertionError in the captured
output; the fix loop would have halted otherwise). Notable lines from
`scripts/output/moving_throat_pde_stage041_rank2_support_sympy_audit.txt`:
- line 15: `determinant decomposition = 0`
- line 22: `n_req - expected = 0`
- line 34: `dn/dm - expected = 0`
- line 45: `tracking collapse = 0`
- line 56: `source-tied formula = 0`  (post-fix, now substantive)
- line 74: `source-tied dn/dm = 0`

**Mathematica:** exit=0 (no `$Failed`/`fail:` lines; transcript ends with
`Stage 041 Mathematica audit passed.`). Notable lines from
`mathematica/output/moving_throat_pde_stage041_rank2_support_mathematica_audit.txt`:
- line 10–11: `determinant decomposition = 0` / `PASS: determinant decomposition`
- line 13–14: `n_req - expected = 0` / `PASS: n_req - expected`
- line 20–21: `dn/dm - expected = 0` / `PASS: dn/dm - expected`
- line 27–28: `tracking collapse = 0` / `PASS: tracking collapse`
- line 34–35: `source-tied formula = 0` / `PASS: source-tied formula`  (post-fix)
- line 39–40: `source-tied dn/dm = 0` / `PASS: source-tied dn/dm`

**Output freshness:** SymPy script mtime 2026-05-22 12:33:51; sympy
output mtime 2026-05-22 12:34:49 (newer, fresh). Mathematica script
mtime 2026-05-22 12:33:52; Mathematica output mtime 2026-05-22 12:34:58
(newer, fresh). Both outputs were re-generated after Codex's edits.

## Material-change assessment

`material_change`: false.

The edit changes only how the source-tied specialization is verified;
the symbolic formulas themselves (`n_expected`, `n_src_expected`,
`dn_dm_src_expected`, the regularity/positivity thresholds, the
tracking collapse `G_q - m`) are unchanged. No derived numerical value
or analytic expression that downstream units could consume has been
altered. Downstream units that may have referenced the unit's
conclusions (Theorem 5, the source-tied feasibility window) are
unaffected: the conclusions are the same; only the proof inside the
script is now algebraically anchored.

## Side observations (non-blocking)

- In the sympy script, the new symbol `t` is defined at line 109 with
  `real=True`. The rest of the symbols in this file declare `positive=True`
  where appropriate; `t` does not need positivity for the substitution
  to work (only `t**2 -> lam0` is used). Not a finding.
- The Mathematica `$Assumptions` now lists `t` but does not assert
  `t > 0`. Since `t` only appears as a transient substitution variable
  and is eliminated via `t^2 -> lambda0` before `nSrc` is compared, the
  missing positivity assumption does not affect the identity. Not a
  finding.
- The sympy displayed `n_req^(src)` (lines 50–55 of the .txt) is now
  written with an overall sign convention opposite to the Mathematica
  form, but the residual after subtracting `n_src_expected` is exactly
  zero in both engines, so the agreement is intact (as the original
  audit also noted for the `n_req` cross-engine comparison). Not a
  finding.

## Verdict justification

Codex applied the directive verbatim in both engines: the `n_src`
construction in section 24.4 now derives from the general `n_expected`
via the physical substitution `q -> t R_U, r -> t` followed by
`t^2 -> lam0`, replacing the previous character-for-character literal
that made the assertion tautological. Both exec transcripts show the
post-fix `source-tied formula = 0` (sympy) and `PASS: source-tied
formula` (Mathematica) without regressing any other assertion, and
both output files are newer than their source scripts. The single
finding is `resolved`; no `partial`, `regressed`, `blocked`, or
`not_attempted` classifications. Verdict: `verified`.
