---
unit_id: 241
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-03T09:10:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 241

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
`scripts/...stage241..._sympy_audit.py:61-64`. The old line-61 assertion
`assert_zero(eps_sel - (1 - sp.Rational(3, 2) * varrho), ...)` (a verbatim
`X - X` against the line-58 `eps_sel` definition) is replaced by
`assert_zero((2 * (1 - eps_star) / 3 - varrho).subs(eps_star, eps_sel), "eps_* = 1 - 3 varrho / 2 solves Stage-240 support law")`.
The `eps_sel` definition at line 58 is retained (per directive); the companion
sigma check at line 65 is untouched.

**Assessment:**
Non-tautological and threads the Stage-240 support law (notes line 164,
`varrho = 2(1 - eps_*)/3`). The new residual contains NO second literal copy of
`1 - 3*varrho/2`; the only remaining `Rational(3, 2) * varrho` in the file is the
kept definition at line 58. The residual is built as `support_law - varrho` with
`eps_sel` substituted for `eps_star`: if `eps_sel` were the wrong branch law the
residual would not vanish, so the check is genuinely falsifiable. Hand-check:
`2(1 - (1 - 3varrho/2))/3 - varrho = varrho - varrho = 0`. Output line 1 now
reads `[ok] eps_* = 1 - 3 varrho / 2 solves Stage-240 support law`.

### F2 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
New file `mathematica/...stage241..._mathematica_audit.wl` (224 lines), exits 0,
38 `[ok]` lines covering M1-M8. New file (untracked), so it does not appear in
the captured `stage_241_diff.patch` (which only carries the `.py` hunk) —
verified directly on disk and via the mathematica exec log.

**Assessment:**
GENUINELY INDEPENDENT, not a transliteration. The independence-critical step M4
DERIVES the crossovers rather than importing the literals the `.py` uses:
- line 110 `epsWRoots = epsStar /. Solve[Together[wW - wLambda] == 0, epsStar]`,
  then `takeNonzeroRoot` drops the trivial `epsStar -> 0` root, then
  `checkZero[..., epsWCross - 1/(2 + beta^2)]` (line 115-118) — solve-for-the-zero.
- line 120 `epsURoots = ... Solve[Together[wUmag - wLambda] == 0, epsStar]`,
  nonzero root checked against `beta/(1 + beta + beta^2)` (line 125-128).
- `varrhoWLambda`/`varrhoULambda` are then obtained by `Solve[epsSelected == epsWCross, varrho]`
  (lines 130-137), i.e. by substituting the support law — not asserted from a
  pre-written form. This is a different decomposition than the SymPy script,
  which imports `1/(2+beta^2)` / `beta/(1+beta+beta^2)` as literals at .py
  lines 121-132 and checks consistency. Distinct var names (epsStar vs eps_star,
  wUmag vs w_Umag, den vs N, dPoly vs D), distinct primitives (Solve/Reduce-style
  derivation + FullSimplify vs simplify-assert choreography). Not a port.
- M7 independently computes `varrhoWLambda /. beta -> 2/11 == 125/369`
  (line 183, `[ok]` in log) — cross-engine corroboration of the F3 notes fix.
- FAIL path non-vacuous: `fail` calls `Exit[1]` (line 12); `checkZero` passes
  only if `res === 0` (line 24); `checkPositive` only if `FullSimplify[res > 0]`
  is `TrueQ` (line 31); `takeNonzeroRoot` calls `fail` unless exactly one nonzero
  root survives (line 40). ConditionalExpression stripped via `stripCE`; beta->0
  via direct substitution, not pole check. No silent parser-skip — all 38 named
  checks emit residual/value + `[ok]`.

### F3 — paper_misalignment (notes_contradicts_script)

**Classification:** resolved

**What changed:**
`notes/stages/...stage241..._sympy_audit.md:577`. The single authorized edit
`\frac{193}{369}` → `\frac{125}{369}`; the `\approx 0.338753` decimal is intact.
`git diff` on the notes file shows exactly one changed line (1 insertion,
1 deletion). The `.py` was NOT changed for F3 (the only `.py` hunk is the F1
edit). Nothing else in `notes/`/`paper/` for stage 241 was touched (the other
modified notes files in the working tree — 231, 232 — belong to other units in
this batch, not unit 241).

**Assessment:**
Correct. `125/369 = 0.338753`, matching the boxed decimal beside it; `193/369 =
0.523` was the typo. M7 in the new `.wl` independently confirms `125/369` exactly.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
`[ok] eps_* = 1 - 3 varrho / 2 solves Stage-240 support law` (F1 fix),
`[ok] varrho_WLambda(beta=2/11) = 125/369`, `All Stage 241 symbolic checks passed.`
32 `[ok]` lines (= A1-A32 in the report's assertion inventory).

**Mathematica:** exit=0. Notable lines:
`[ok] M4 Solve[wW - wLambda == 0] nonzero root`,
`[ok] M4 Solve[wUmag - wLambda == 0] nonzero root`,
`[ok] M7 varrhoWLambda(beta=2/11)`, `Stage 241 Mathematica audit passed.`
38 `[ok]` lines spanning M1-M8.

**Output freshness:** confirmed. sympy `.py` mtime 08:46:16, output `.txt` 08:56:23;
mathematica `.wl` mtime 08:46:17, output `.txt` 08:56:23. Both outputs newer than
their scripts.

## Material-change assessment

`material_change`: false. F1 strengthens a check without altering any derived
value (the same `eps_sel = 1 - 3varrho/2` is verified, now via the support law).
F2 is purely additive (a second engine confirming the same identities). F3 is a
notes-prose typo correction with no script change. No downstream-visible result
changed; downstream units depending on stage 241's thresholds are unaffected.

## Side observations (non-blocking)

- The captured `stage_241_diff.patch` contains only the `.py` hunk because the
  `.wl` (F2) and the notes file (F3) are/were untracked-or-prose; both edits were
  verified directly on disk and via the exec logs/git diff. No issue, just noting
  the patch is not the complete picture for this unit.

## Verdict justification

All three findings are resolved. F1 is now non-vacuous and threads the Stage-240
support law with no residual literal copy; F2 adds a genuinely independent
Mathematica route whose M4 derives the crossovers via `Solve` (not by importing
the `.py` literals) and whose M7 cross-corroborates `125/369`; F3 applies exactly
the one authorized notes edit with the decimal intact and no script change. Both
exec logs exit 0 with non-vacuous, fully-emitted check sets (32 sympy, 38
mathematica), outputs are fresh, and no regression appears in the diff.
material_change is false.
