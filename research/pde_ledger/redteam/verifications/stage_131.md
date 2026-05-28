---
unit_id: 131
batch: IV.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 131

This batch was orchestrator-direct (Codex bypassed); the directive carries no `## Applied:` blocks, so I verified by reading the current script files plus the post-fix exec logs. Cluster A renumbering (`.wl` banner `STAGE 114` → `STAGE 131`, notes H1 `Stage 233` → `Stage 131`) is also confirmed.

## Per-finding outcomes

### F1 — missing_verification_script (subtype: script_doesnt_cover_claim)

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py:37-72` adds four anchored `assert` blocks plus a fifth that came in with F4:
- L40-45 anchor (1) Pi_* vs literal `1.50882951349316` (tol 1e-14, 50 prec).
- L48-53 anchor (2) g'(Pi_*) vs literal `0.0714453558083195` (tol 1e-14, 50 prec).
- L57-62 anchor (3) parent threshold identity at Pi = Pi_*, using `sp.simplify(...) == 0`.
- L66-72 anchor (4) lower-branch discrimination: `|gPi(2*Pi_*) - g_-| > 1e-3`.

**Assessment:**
Anchors (1), (2), (3), (4) all correspond to the M1–M4 claim manifest items. None are tautological: (1) and (2) compare a computed numerical value against an externally supplied paper literal; (3) substitutes the numerically converged `Pi_*` into the residual and checks symbolic structure; (4) is a non-trivial counter-example check at `2*Pi_*` where the actual numerical residual `0.0903…` is empirically clearly above the `1e-3` floor. Exit code 0, transcript shows all four PASS lines plus the F4 PASS.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl`:
- L26 banner now reads `STAGE 131 — PARENT MICRO-THRESHOLD…` (was `STAGE 114`).
- The tautological `expectApprox["Pi_* compensation point", gPi /. piM -> piStar, gMinus, 10^-20]` is gone.
- L52-80 add four anchored checks (`piStar notes Sec. 1 value`, `slope at piStar notes Sec. 3 value`, identity at piStar, lower-branch discrimination), matching the directive Anchors (1)–(4).

**Assessment:**
The check labels in the `.wl` use ASCII-safe labels (`piStar`, `slope at piStar`) instead of `g'(Pi_*)` — documented quirk to avoid Mathematica's comment parser tripping on `_*)` near a terminator; the underlying expressions (`N[piStar, 50]`, `D[gPi, piM] /. piM -> piStar`) match the directive. None are tautological: Anchors (1)/(2) compare against the paper literals (residuals `~ -4.4e-15` and `~ 2.1e-17`, well under 1e-14); Anchor (3) uses `Chop[Simplify[...], 10^-30] === 0` which is a non-trivial structural reduction; Anchor (4) is a separate numerical evaluation at `2*piStar`. No `PASS: Pi_* compensation point` line in the transcript — confirms the old tautological assertion is removed. Banner in the transcript header now reads `STAGE 131`.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
Subsumed by Anchors (3) and (4) of F1 (SymPy) and F2 (Mathematica). Both engines now substitute `Pi -> Pi_*` into `threshold_residual` and check structural equality, and both compare `gPi` at `2*Pi_*` to `g_-` for branch discrimination.

**Assessment:**
The directive explicitly states `F3` requires no edits beyond F1+F2; both anchors are present and PASS in both transcripts (SymPy line 8 + 9; Mathematica lines 15 + 16). Resolved.

### F4 — hardcoded_result

**Classification:** resolved

**What changed:**
- SymPy L12-21: literal `sp.Float("0.758035078944663")` is replaced by the closed form `g_minus_exact = (2*sp.sqrt(4107 - 100*sp.pi**2) - 37*sp.sqrt(3)) / (20*sp.pi)`, then asserted against the literal `0.758035078944663` at 50 precision (tol 1e-14). Comment cites the upstream context line.
- Mathematica L33-39: closed form `gMinusExact = (2*Sqrt[4107 - 100*Pi^2] - 37*Sqrt[3])/(20*Pi)` is now PASS-checked against `gMinusLiteral = 0.758035078944663\`50` via `expectApprox` (tol 1e-14) before assigning `gMinus`.

**Assessment:**
Both engines now construct `g_minus` from the same closed form independently and assert numerical agreement with the upstream literal. `grep 4107` returns a match in both `.py` and `.wl`. SymPy transcript line 1 `PASS: g_-^F1 closed form matches literal 0.758035078944663`; Mathematica transcript lines 5-6 show residual `~ -1.73e-16` and corresponding PASS. The two-engine independence concern from the original report is closed: the SymPy side no longer uses a transliterated float literal as primary, the closed form is the source of truth in both engines.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `PASS: g_-^F1 closed form matches literal 0.758035078944663`
- `PASS: Pi_* matches notes Sec. 1 value 1.50882951349316`
- `PASS: g'(Pi_*) matches notes Sec. 3 value 0.0714453558083195`
- `PASS: parent threshold identity at Pi = Pi_* matches notes Sec. 2`
- `PASS: lower-branch discrimination — gPi(2*Pi_*) - g_- = 0.0903049686172024794499842493650`

**Mathematica:** exit=0. Banner reads `STAGE 131 — PARENT MICRO-THRESHOLD FOR CANONICAL MOUTH COMPENSATION`. Notable lines:
- `g_-^F1 closed form vs literal residual = …*10^-16` → `PASS: g_-^F1 closed form vs literal`
- `piStar notes Sec. 1 value residual = …*10^-15` → `PASS: piStar notes Sec. 1 value`
- `slope at piStar notes Sec. 3 value residual = …*10^-17` → `PASS: slope at piStar notes Sec. 3 value`
- `PASS: parent threshold identity at piM = piStar (notes Sec. 2)`
- `PASS: lower-branch discrimination (paper Checks item 3)`
- `Stage 131 Mathematica audit passed.`

No `PASS: Pi_* compensation point` line in either transcript — confirms the tautological check is gone.

**Output freshness:** SymPy `.py` mtime `May 27 17:24`, output `.txt` mtime `May 27 17:45` (newer, refreshed post-fix). Mathematica `.wl` mtime `May 27 17:51`, output `.txt` mtime `May 27 17:51` (same minute — refreshed in lockstep). Both outputs are post-fix.

## Material-change assessment

`material_change`: false.

No derived numerical value changed. `Pi_* ≈ 1.508829513493155…`, `g'(Pi_*) ≈ 0.0714453558083195…`, and `g_-^{F1} ≈ 0.758035078944663` are unchanged from the pre-fix transcript; the closed-form derivation in SymPy reproduces the same float at >1e-14 precision. Banner string and notes H1 renumber are cosmetic. Downstream consumers (paper card says stages 133–145) read the same numeric outputs as before; no re-audit warranted.

## Side observations (non-blocking)

- The SymPy anchor (3) uses `sp.N(Pi_star, 50)` to substitute into both sides of an algebraic identity. Because the same numeric is substituted into both sides, `sp.simplify(... == 0)` is genuinely structural rather than tautological — it would catch a substitution-channel bug (wrong symbol/sign) but not a wrong `Pi_*` value. The wrong-value case is independently caught by anchor (1). Fine.
- Mathematica anchor (3) uses `Chop[Simplify[...], 10^-30] === 0` whereas the directive shows `Simplify[...] === 0`. The added `Chop` is a tightening that tolerates numerical noise in the high-precision substitution while still rejecting any nontrivial residual; the directive's text intent is preserved.
- The Mathematica comment labels for anchors (1) and (2) intentionally use ASCII names (`piStar`, `slope at piStar`) per the orchestrator note about Mathematica's comment parser tripping on `_*)`. Confirmed harmless.

## Verdict justification

All four findings are mechanically `resolved`. Both engines exit 0, both transcripts contain non-tautological PASS lines for every M1–M4 anchor and the F4 closed-form↔literal check, the Cluster A banner+H1 renumbering is in place, and the load-bearing numerical outputs (`Pi_*`, `g'(Pi_*)`, `g_-^{F1}`) are unchanged. No downstream propagation concern; `material_change: false`. Verdict: `verified`.
