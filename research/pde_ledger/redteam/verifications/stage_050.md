---
unit_id: 050
batch: III.2
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-26
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 050

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:63-67` — derivative check now asserts `dZetaNdx (2n+1)^2 (1 + n(n+1) x)^2 + n(n+1) = 0` instead of subtracting `dZetaNdxTarget`. `xEqClosedForm` removed from any `expectZero` argument; replaced at wl:70-73 with `(2n+1)^2 zetaReq (1 + n(n+1) xEq) - 1`. Ceiling check at wl:90-107 decomposes `ceilingDiff` via `Numerator[Together[...]]` and `Denominator[Together[...]]`, asserting each piece against the expected numerator `(1-eps)(2n+1)^2 n(n+1) x` and the expected denominator `((2n+1)^2 - eps)((2n+1)^2 (1 + n(n+1) x) - eps)` separately.

**Assessment:**
The three rewrites match the directive's "After" blocks (a), (b), (c) byte-for-byte. All three new residuals are non-tautological — they exercise the defining equation of `zetaN`/`xEq`/the ceiling difference directly rather than re-typing the SymPy targets. `dZetaNdxTarget` and `ceilingDiffTarget` are removed from the script entirely; `xEqClosedForm` is kept only as a local declaration referenced by the `Print` at wl:68 (which the directive explicitly permits). Output file (mtime 02:58 > script 02:56) shows all three new labels as `PASS`. No collateral edits beyond the three blocks.

### F2 — paper_misalignment

**Classification:** resolved

**What changed:**
Per the `## Applied: F2` block: user selected direction (a); paper-card edit landed at `paper/stages/stage_050.tex:43-51` (new boxed `S_n^{\rm twin}(x;\varepsilon) < S_n^{\rm max}(\varepsilon)` with label `eq:app-stage050-Sn-max`; Output line at line 50 updated). No script change made (correctly — F2 was a paper-card scope question, not a script defect).

**Assessment:**
Per the verifier prompt's special note, F2 is treated as `resolved` when the `## Applied: F2` block is present and consistent with direction (a). The Applied block records files_changed, direction, and `deviation: none`. The scripts retain the existing ceiling block (sympy:88-107, wl:85-107) as instructed by the "Codex must NOT re-apply F2" footer. Consistent.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:47-50` inserts `expect_zero("zeta_0^(twin) - 1 (anchors doubling to stage 049 import)", twin_support_ratio(sp.Integer(0), x) - 1)` immediately before the existing `S(1;eps) - 2` assertion. `mathematica/...wl:44-47` inserts the analogous `expectZero["zeta_0^(twin) - 1 (anchors doubling)", (1/((2 n + 1)^2 (1 + x n (n + 1))) /. n -> 0) - 1]` before the existing wl:48 `S(1;eps) - 2` assertion.

**Assessment:**
Both insertions match the directive's required code exactly. The SymPy anchor evaluates the imported `twin_support_ratio(0, x)` (not a hand-written `1`), so a future drift in the stage 049 module would surface here. The Mathematica analogue uses the literal `n -> 0` replacement pattern which fires before `FullSimplify` evaluates anything — the resulting residual `1/(1*(1+0)) - 1 = 0` reduces cleanly. New labels appear as PASS lines in both output transcripts (sympy output line 11, mathematica output lines 7-8). Non-tautological — the SymPy version anchors to the imported upstream symbol.

## Exec log assessment

**SymPy:** exit=n/a (exec log `redteam/exec_logs/stage_050_sympy.log` not present). Fallback: `scripts/output/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.txt` (mtime 02:58 > script mtime 02:56) shows all expected PASS lines including the new anchor at line 11 (`zeta_0^(twin) - 1 (anchors doubling to stage 049 import) = 0`) and the final tagline `All Stage-33 symbolic checks passed.`

**Mathematica:** exit=n/a (exec log `redteam/exec_logs/stage_050_mathematica.log` not present). Fallback: `mathematica/output/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.txt` (mtime 02:58 > script mtime 02:56) shows PASS for: `zeta_0^(twin) - 1 (anchors doubling)`, `d zeta_n / dx (denominator structure) ...`, `x_max from Solve satisfies (2n+1)^2 zeta_req (1 + n(n+1) x_max) - 1 = 0`, `Numerator of (S_n^(max) - S_n^(twin)) - (1-eps)(2n+1)^2 n(n+1) x`, `Denominator of (S_n^(max) - S_n^(twin)) - ((2n+1)^2 - eps) ((2n+1)^2 (1 + n(n+1) x) - eps)`, plus all retained PASS lines. Final tagline `Stage 050 Mathematica audit passed.`

**Output freshness:** confirmed. Sympy script mtime 02:56 < output 02:58. Mathematica script mtime 02:56 < output 02:58. Both outputs were regenerated post-fix.

## Material-change assessment

`material_change`: false.

F1 changes the algebraic pathway of three Mathematica checks but the mathematical content (the identities being verified) is identical to the pre-fix script — same `zetaN`, same `xEq` from `Solve`, same `ceilingDiff`. F3 adds new anchor assertions but does not alter any existing claim or derived expression; both new anchors evaluate to 0 trivially under their own algebra. F2 is paper-side only. No downstream unit's input/output card changes.

## Side observations (non-blocking)

- The `Print` line at wl:68 still references `xEqClosedForm` (declared at wl:60) for human readability of `x_max(n;zeta_req)`. This is explicitly permitted by the directive ("the line wl:56 declaration of `xEqClosedForm` can remain so the `Print` statement at wl:65 still references it"). Not a finding.
- Output line 15 of the Mathematica transcript prints the long label followed by `= 0 = 0` because the label string itself contains a literal `= 0` suffix. Cosmetic only; does not affect PASS logic.
- A previous verification file existed at this path from an earlier iteration with a different 4-finding structure (referring to `tautological_check` findings on cancellation, `xEq` self-definition, ceiling, and monotonicity). It has been overwritten with this 3-finding verification matching the current report and directive.

## Verdict justification

All three findings are resolved: F1's three load-bearing Mathematica checks now use independently-derived residuals matching the directive's "After" blocks; F2 was closed by a user-approved paper-card edit recorded in the `## Applied: F2` block per the verifier prompt's special instruction; F3's anchor lines are present in both scripts and PASS in the saved outputs. No exec logs were captured by the orchestrator, but the output `.txt` files are freshly regenerated (mtimes newer than script mtimes) and contain every expected PASS line plus the final pass tagline in each engine. No regressions in the diff. No material change to derived results downstream.
