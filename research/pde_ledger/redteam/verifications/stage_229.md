---
unit_id: 229
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T18:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 229

## Per-finding outcomes

### F1 — paper_misalignment (notes_contradicts_script)

**Classification:** resolved

**What changed:**
Notes-only edit (user-authorized direction (a), per the `## RESOLVED` block in the directive). In
`notes/stages/moving_throat_pde_stage229_..._sympy_audit.md:292` the crossover-cubic leading term
`189\xi^3` was changed to `121\xi^3`. No script change. Confirmed by grep: the file now contains
`121\xi^3+297\delta\xi^2+333\delta^2\xi+81\delta^3-72\delta^2` at :292 and the unchanged derivative
`363\xi^2+594\delta\xi+333\delta^2 > 0` at :302; no occurrence of `189` or `567` remains in the file.

**Assessment:**
Correct and minimal. Only the leading 189→121 changed; the rest of the cubic and the next-line
derivative are untouched (acceptance criteria met). The fix direction is corroborated three ways:
(i) the notes' own derivative `363\xi^2 = ∂(121\xi^3)` (not `567\xi^2 = ∂(189\xi^3)`), so the notes
are now internally consistent; (ii) the SymPy script re-derives the cubic in-script
(`numerator = expand(-together(R_ND-1)...)` at py:95) and asserts equality to `expected_P` with leading
`121` at py:96-97 — non-tautological, the literal is forced by the recomputed numerator; (iii) the NEW
Mathematica `.wl` independently derives the same cubic with leading coefficient `121` (see F2 below).
Cross-engine corroboration of `121` is decisive.

### F2 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
New independent `.wl` created at
`mathematica/moving_throat_pde_stage229_..._mathematica_audit.wl` covering M1–M10. It is NOT a
transliteration of the `.py`: it builds `sMinus`/`nMinus`/`F` from the raw constants
`kappa0sq=8/Pi^2`, `kappa1sq=16/(9 Pi^2)` (lines 61-78), uses native `D[Log[...],xi]` for the
log-slopes (lines 92-93), `Cancel[Together[...]]` for `RND` (line 103), one-sided
`Limit[..., Direction->"FromBelow"]` for the softening limits (lines 114-128), `Numerator[Together[RND-1]]`
for the cubic (line 137), `CoefficientList` for coefficient extraction, `Resolve[ForAll[...]]` for
monotonicity (lines 170-179), and `NSolve`/`Reals` root selection (lines 45-50). Different decomposition
order and idiom from the SymPy script.

**Assessment:**
Genuine independent verification. All 27 in-file checks PASS, exit 0.

- **M7 (crossover cubic) — DERIVED, not hardcoded.** `rawP = Expand[Numerator[Together[RND-1]]]`
  (line 137) is produced natively from `RND` (itself from `D[Log[FNum]]/D[Log[FDen]]`). The literal
  cubic appears only as `targetCoeffs` (lines 141-146) used for the coefficient comparison, never as
  the source. Sign-normalization (`If[Coefficient[rawP,xi,3]<0, -rawP, rawP]`, line 139) is applied
  because the raw numerator came out negative-leading (log: `72*delta^2 - ... - 121*xi^3`), yielding
  `+121*xi^3`. The log shows coefficient list `{-72δ²+81δ³, 333δ², 297δ, 121}` — leading `121`,
  independently corroborating the F1 typo direction.
- **M8 (∂_ξP>0 monotonicity) — genuine proof, not assumed.** `dP = Expand[D[P, xi]]` is the derivative
  of the *derived* P; its coefficients are checked, then positivity is proved by
  `Resolve[ForAll[{xiPos,deltaPos}, Implies[xiPos>=0 && deltaPos>0, dP>0]], Reals]` returning literal
  `True` (log line 51). Fresh quantifier symbols are used so the global `xi<1` assumption does not
  weaken the statement — positivity is proved for all `xi>=0`, stronger than required.
- **Token convention non-vacuous.** `fail` calls `Exit[1]` (line 13). `expectZero` `Exit[1]`s when the
  simplified residual `=!= 0` (line 30); `expectTrue` fails unless the result is literally `True`
  (`TrueQ`, line 36); `expectNumeric` fails unless `err < tol`; `stableRoot` fails unless exactly one
  root in (0,1); M10 side probes `fail[]` if the left/right classification is wrong. A bad check
  genuinely aborts with exit 1.
- **No fabricated literals.** classifier/log-slope/onset/threshold/cubic targets are all compared
  against natively-derived expressions; M10 sample roots are roots of the in-script derived polynomial
  `P`, matched to <5e-13 (log: errors ~1e-16).

## Exec log assessment

**SymPy:** exit=0. Notable lines: `P(xi,delta) = ... + 121*xi**3`, `dP/dxi = ... + 363*xi**2`,
`P(0,delta) = delta**2*(81*delta - 72)` (= `9δ²(9δ-8)`), sample roots `0.107223...`, `0.081847...`,
`0.032505...`, "All Stage 229 symbolic and numerical audits passed."

**Mathematica:** exit=0. 27 PASS, 0 FAIL. Notable: `M7 coefficient 3 residual = 0` / `PASS M7
coefficient 3` (leading `121`), `M8 derivative positivity = True` / `PASS`, `M9 P(0,delta) =
9*delta^2*(-8 + 9*delta)`, `M10 crossover root delta = 1/4 error = 5.55e-17`, `M10 delta=1 always
denominator slice = True`.

**Output freshness:** confirmed. Both output `.txt` mtimes (1780443634) are newer than the `.wl`
(1780443250) and the `.py` (1778522211). Note: the captured `stage_229_diff.patch` is 0 bytes — expected,
since the `.wl` and its output are new/untracked files (not in a tracked diff); the working tree shows
them present and the exec logs show the post-fix runs.

## Material-change assessment

`material_change`: false. F1 is a notes-only prose typo correction (189→121) where the script already
held the correct value `121` — no derived result changed. F2 adds a second engine that confirms the
existing SymPy results; no derived quantity changed. Downstream units depending on stage 229's cubic
already saw `121` (the script value). No re-audit of downstream units is warranted on account of 229.

## Side observations (non-blocking)

- The `.wl`'s M8 positivity proof is over `xi>=0` (all reals ≥0), broader than the stable-branch
  domain `0<=xi<1`; this is strictly stronger and fine.
- The `.wl` computes `M8 dP/dxi = 3*(111*delta^2 + 198*delta*xi + 121*xi^2)` (factored display) which
  expands to `333δ²+594δξ+363ξ²`, matching the SymPy `dP/dxi` and the notes' (unchanged) `363ξ²`
  derivative — consistent.
- M6 uses the reciprocal pole trick (`1/LDenSoft == 0`) per the Mathematica idiom guidance, rather than
  `=!= Infinity`. Good.

## Verdict justification

Both findings are `resolved`. F1's notes typo is corrected to `121\xi^3` with the rest of the cubic and
the derivative untouched, and the correction is corroborated by the notes' own derivative, the SymPy
in-script re-derivation, and now an independent Mathematica derivation. F2's new `.wl` is a genuine,
non-transliterated independent verification: the crossover cubic is derived (not hardcoded) and produces
leading `121`, the monotonicity is a real `Resolve[ForAll]` proof, the success helper genuinely `Exit[1]`s
on failure, and all 27 checks pass with exit 0. Both exec logs pass, outputs are fresh, no regressions.
`material_change` is false (notes-only typo + additive second engine; script value unchanged). Verdict:
**verified**.
