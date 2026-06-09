---
unit_id: 211
batch: VI.1
verifier_model: Opus 4.8 (1M context) [claude-opus-4-8[1m]]
verify_date: 2026-06-09T17:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 211

The directive carries one actionable finding (F1 = re-author the `.wl` to break the
`mathematica_transliteration`); F2 is the card-text-lag class, DEFERRED to PAPER_CLEANUP P4-51 with
no Codex action this batch. I confirmed the eliminants are now DERIVED via `Resultant`, that the
SymPy closed forms survive only as labeled comparison targets, that all original deliverables are
preserved, that the `.py` is untouched, and that both engines run clean — quoting `.wl`/`.py`/output
line numbers.

## Per-finding outcomes

### F1 — mathematica_transliteration (load-bearing; re-author the `.wl`)

**Classification:** resolved

**What changed (per `stage_211_diff.patch`, single `diff --git` block — the `.wl` only):**

- **§M1 (wl:85-89):** the posited `Nr = 2 Mr sqrtDelta + Lr` / `Ns` (old wl:87-92, mirroring
  `.py:53-58`) are DELETED. The stationary numerators are now obtained only from differentiation:
  `numR = Numerator[Together[D[Phi, r]]]`, `numS = Numerator[Together[D[Phi, s]]]` (wl:85-89). The
  old `expectZero["… numerator minus paper N_r", directNumerator[r] - Nr]` posit-comparison is
  replaced by a nonzero check (wl:93-94) + a reconstruction `D[Phi,r] - numR/stationaryDenominator[r]`
  (wl:95-96).
- **§M2 cross (wl:100-125):** the load-bearing quartic is now DERIVED — substitute `Sqrt[Delta] -> q`
  in the numerators (`numRq/numSq`, wl:100-101) and `crossDerived = primitivePartRS[Resultant[numRq,
  numSq, q]]` (wl:102). The SymPy form `Ms Lr - Mr Ls` appears ONLY at wl:113-117 under the comment
  `(* SymPy comparison target only. *)` as `crossTarget`, used in the labeled cross-check
  `crossDerived - crossRatio crossTarget == 0` (wl:125), NOT as the source of the tested polynomial.
- **§M3 squares (wl:129-159):** `srDerived/ssDerived = primitivePartRS[Resultant[numRq/numSq,
  q^2 - Delta, q]]` (wl:129-130) — the sextic square eliminants DERIVED by eliminating the radical
  symbol `q` against `q^2 - Delta`. The SymPy forms `Lr^2 - 4 Mr^2 Delta` / `Ls^2 - 4 Ms^2 Delta`
  appear ONLY as `squareTargetR/squareTargetS` (wl:143-144), again as labeled comparison targets in
  `srDerived - srRatio squareTargetR == 0` (wl:158-159).
- **§M4 Bezout (wl:163):** `bezoutBound = crossDegree srDegree` — the product of the DERIVED
  polynomials' degrees, not a hardcoded `4*6`.
- **§M5/§M6:** M5 (wl:177-188) keeps the iso substitution+`FullSimplify`/`PowerExpand` reduction
  (directive explicitly permits this — it is a substitution, not the load-bearing eliminant); M6
  (wl:203-207) now substitutes into the DERIVED `numR`/`numS` rather than a posited `Nr`/`Ns`.

**Assessment — INDEPENDENT? yes.** The single discriminating operation now present:
`Resultant[numRq, numSq, q]` (and `Resultant[numRq, q^2 - Delta, q]`), i.e. radical elimination of the
derived differentiated numerators. The `.py` takes the opposite route — it POSITS
`Nr = 2 Mr sqrtDelta + Lr` (`.py:57`), `Ccross = Ms Lr - Mr Ls` (`.py:83`), `Sr = Lr^2 - 4 Mr^2 Delta`
(`.py:84`) and checks them by algebraic identity (`.py:94-96`), with NO `Resultant`/`Eliminate`. So
the `.wl` now DERIVES the eliminants the `.py` POSITS — genuine derive-vs-posit independence, not a
relabel. The cross-check direction is the right one: the derived polynomial is the tested object, and
the SymPy closed form is only the comparison anchor (ratio `-1`, residual `0`).

**Non-tautology / genuine-reconstruction self-check.** The derived cross is a real degree-4 object
(output L29-30) not equal to the target by construction — the script must compute `crossRatio = -1`
and verify `crossDerived - crossRatio crossTarget == 0` (a wrong derivation or a sign error would
leave a nonzero residual). The sextic ratios are also `-1` (output L54-55, the conventional
Resultant-vs-`L^2-4M^2 Delta` sign), with residual `0`. The `M1` reconstruct asserts `D[Phi,r] -
numR/stationaryDenominator[r] == 0` against an INDEPENDENT `Denominator[Together[D[Phi,r]]]` rather
than a re-typed denominator, so it is non-tautological.

**Deliverables preserved — confirmed:**
- quartic cross degree 4 — output L30/L35 (`M2 total degree derived cross = 4`, PASS).
- sextic `S_r`/`S_s` degree 6 — output L48-49/L50-53 (`M3 total degree derived S_r/S_s = 6`, PASS) —
  now from the DERIVED polynomials, not posited forms.
- Bezout 24 — output L68/L73 (`M4 computed Bezout product = 24`, PASS), from derived degrees.
- M5 diag-iso reduction — output L75-78 (`Delta_iso = 0`, `tau_iso = 0`, PASS).
- M6 symmetric equal-mix `Nr(1,1)=Ns(1,1)=0` — output L83-86 (PASS), now on derived numerators.
No FAIL anywhere; both engines end "All Stage 211 identities verified."

## Exec log assessment

**SymPy:** exit=0 (`stage_211_sympy.log` L236 `# exit_code: 0`). The `.py` is the UNCHANGED reference
(the diff touches only the `.wl`); its refreshed output (banner `STAGE 211`, L42-43 `dPhi/dr/ds
compiler = 0`, L210-215 cross/square identities `= 0` + `deg 4/6/6`, L232 `Bezout … = 24`) is the
engine cross-check. The `.py` still POSITS the eliminants — the intended asymmetry vs the re-authored
`.wl`.

**Mathematica:** exit=0 (`stage_211_mathematica.log` L94 `# exit_code: 0`). Notable lines:
- L34-37: derived cross degree 4, PASS; L38-42: `derived cross / SymPy target = -1`, ratio nonzero
  constant PASS, `derived cross minus scaled SymPy target = 0` PASS.
- L51-68: derived `S_r`/`S_s` degree 6 PASS, ratios `-1`, residuals `0` PASS.
- L73-75: `M4 computed Bezout product = 24` PASS. L80-91: M5/M6 all `= 0`/PASS.
Every `expectZero`/`expectTrue` passes; no `FAIL`.

**Output freshness:** the committed `mathematica/output/…211…txt` mtime is 2026-06-09 16:51,
newer than the re-authored `.wl` (16:07); the exec-log re-run is dated 16:43. The `.py` mtime is
2026-06-02 (predates the re-author) — confirming it was NOT edited.

## Material-change assessment

`material_change`: false. Only the METHOD that produces the cross/square eliminants changed (posit →
`Resultant` derivation); every emitted deliverable value is identical — derived cross degree 4,
sextic degree 6, Bezout 24, the M5 reduction zeros, and the M6 equal-mix zeros all still hold, and the
derived polynomials match the SymPy targets up to the constant `-1`. No downstream-carried value
moves; no constant/identity relocated.

## Side observations (non-blocking)

- The derived `S_r`/`S_s` resultants carry the conventional overall sign `-1` relative to the SymPy
  `L^2 - 4 M^2 Delta` convention (output L54-55); the script correctly absorbs this via the
  `crossRatio`/`srRatio` scalar and asserts the scaled residual `= 0`, so it is a labeled convention
  difference, not a defect.
- F2 (card "Mathematica audit: none yet") is correctly left untouched per the directive's P4-51
  deferral — not a value/identity mismatch, out of scripts-only scope.

## Verdict justification

The re-author lands the one independence requirement with no relabeling: the quartic cross and both
sextic square eliminants are now DERIVED by `Resultant` over the radical symbol `q` applied to the
differentiated stationary numerators (`Numerator[Together[D[Phi,var]]]` with `Sqrt[Delta]->q`), the
single discriminating op the `.py` does not use — the `.py` instead POSITS `Nr`/`Ccross`/`Sr` and
checks algebraic identities, so the engines are now genuinely derive-vs-posit. The posited SymPy
closed forms survive only as clearly-labeled comparison targets (ratio `-1`, residual `0`), the
Bezout `24` falls out of the derived degrees rather than a hardcoded `4*6`, and every original
deliverable is preserved (deg 4, deg 6, Bezout 24, M5 diag-iso, M6 equal-mix). The `.py` is untouched
(single `.wl` diff block; `.py` mtime predates the re-author), both engines exit 0 with all PASS and
no FAIL, and the committed outputs are fresh. Verdict: verified; material_change: false.
