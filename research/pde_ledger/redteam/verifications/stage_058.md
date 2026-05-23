---
unit_id: 058
batch: III.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 058

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl:55-56` now derives `fc`, `fs` via `FullSimplify[Integrate[Exp[Pe*x]*Cosh[alpha*x], x], Assumptions -> $Assumptions && Pe != alpha]` (and the sinh mirror) rather than transcribing the SymPy closed-form ansatz. Lines 65-73 rebuild `delta = FullSimplify[Integrate[kernel*sigmaPe, {x, 0, 1}, ...]]` from the physical definition and add a `"delta independent integral matches combination form"` regression against the old SymPy-combination form. Antiderivative regressions are relabelled at lines 57-58 ("Ic/Is antiderivative regression (Mma re-derivation)"). Shared cross-engine labels are relabelled with `(Mma re-derivation)` at lines 47, 53, 82, 84, 125 (Kprime, Sigma normalization, Delta0 formula, Delta0 integral identity, weak-coupling constant term). Delta_inf label is updated in F3.

**Assessment:**
Edit matches the directive verbatim — independent `Integrate[]` calls now stand in for the transliterated ansatz, the new combination-form regression passes (output line 19 "PASS: delta independent integral matches combination form"), and all relabelled assertions still pass. No collateral edits beyond the directive's relabel list.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
SymPy `scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py:88-107` now contains the bracket-gap closed-form check, the positivity sweep over alpha in {1/10, 1, 3} × eta in {1/10, 1, 10}, and the `Delta_inf as Pe -> oo limit` assertion (the latter is computed in exp-form to handle SymPy's exponential representation; the deviation is documented in the directive's `## Applied: F2` block). Mathematica mirror at `.wl:101-121` adds `bracketGap`/`bracketGapExpected`, the numerical positivity sweep (with `If[AnyTrue[...]] → fail/pass`), and `deltaInfLimit = FullSimplify[Limit[delta, Pe -> Infinity, ...]]` with the `"Delta_inf as Pe -> oo limit"` assertion.

**Assessment:**
Both engines exit 0 with the three new assertions passing (SymPy output lines 22-25; Mathematica output lines 31-36). The Pe -> oo limit is the genuine physical claim ("sharp-bottom endpoint"), and the bracket-gap positivity sweep substantively tests that `Pe_lo <= Pe_hi`. New checks are non-tautological: `Delta_inf_limit` is derived by `sp.limit(Delta, Pe, oo)` and `Limit[delta, Pe -> Infinity, ...]` from the full `Delta` expression rather than from `K.subs(x,1)`, so a sign error in Ic/Is would produce a non-zero residual. SymPy's `.rewrite(sp.exp)` is a representational tweak (not a tautology hack) — the residual still must vanish symbolically.

### F3 — tautological_check

**Classification:** resolved

**What changed:**
SymPy line 80 label changed to `"Delta_inf direct substitution (sanity)"`; Mathematica line 94 label changed to `"Delta_inf direct substitution (sanity, Mma re-derivation)"`. The substantive replacement lives in the F2 Pe → ∞ limit check.

**Assessment:**
Labels match the directive verbatim; the underlying substitution check is unchanged (it still passes — see SymPy output line 19 and Mathematica output line 27-28). The docstring's "sharp-bottom endpoint" semantic claim is now carried by the F2 limit assertion, so the tautological substitution is correctly demoted to a sanity print.

### F4 — tautological_check

**Classification:** resolved

**What changed:**
SymPy `.py:113-120` adds `Pe1_coeff = sp.simplify(sp.expand(Delta_series).coeff(Pe, 1))`, prints it, and asserts non-vanishing at `alpha=eta=1` via `sp.N(...)` against a `1e-8` threshold. Mathematica `.wl:126-132` mirrors with `pe1Coeff = FullSimplify[SeriesCoefficient[deltaSeries, {Pe, 0, 1}], ...]`, prints it, and uses `If[Chop[pe1Val] === 0, fail[...], pass[...]]`.

**Assessment:**
The new first-order coefficient assertion is non-tautological: a sign error in Ic/Is or in the Delta combination would shift the Pe^1 coefficient and could cancel it at `alpha=eta=1`. Both engines print closed-form Pe^1 coefficients (SymPy line 28, Mathematica line 40) and both pass (`PASS: weak-coupling first-order coefficient nonvanishing at alpha=eta=1`). The original `"weak-coupling constant term"` check is retained (analyticity-guaranteed; directive only required *augmentation* via the Pe1 check, not replacement) — no regression. Note SymPy uses `Abs(Pe1_val) < 1/10^8` threshold rather than literal `!= 0`, which is the safer numeric form; this is a benign formatting choice consistent with the directive's code block.

## Exec log assessment

**SymPy:** exit=0 (the refresh step in `fix_batch_III.2.log` recorded successful regeneration at 17:48:39, and the saved output txt contains the trailing `"Stage 41 audit passed."` line). Notable lines from the saved output:

- L23 `bracket gap positivity sweep = PASS`
- L25 `Delta_inf as Pe -> oo limit = 0`
- L29 `weak-coupling first-order coefficient nonvanishing at alpha=eta=1: PASS`
- L31 `Stage 41 audit passed.`

**Mathematica:** exit=0 (refresh recorded at 17:49:08; saved output ends with `"Stage 058 Mathematica audit passed."`). Notable lines from the saved output:

- L19 `PASS: delta independent integral matches combination form`
- L32 `PASS: bracket gap closed form`
- L33 `PASS: bracket gap positivity sweep`
- L36 `PASS: Delta_inf as Pe -> oo limit`
- L41 `PASS: weak-coupling first-order coefficient nonvanishing at alpha=eta=1`
- L43 `Stage 058 Mathematica audit passed.`

**Output freshness:** Confirmed via `stat -c '%Y'`:
- SymPy `.py` mtime 1779493612 < SymPy txt mtime 1779493748 (newer by 136s)
- Mathematica `.wl` mtime 1779493485 < Mathematica txt mtime 1779493766 (newer by 281s)

Both outputs were regenerated after Codex's edits landed.

## Material-change assessment

`material_change`: false.

The substantive derived quantities (`Delta`, `Delta_0`, `Delta_inf`, `Ic`, `Is`, the bracket endpoints) are unchanged in closed form. The Mathematica re-derivation via `Integrate[]` produces the same antiderivative formulas and the same `delta` (validated by the `"delta independent integral matches combination form"` zero residual). The F2 additions are *new* assertions confirming previously-printed quantities; they do not alter any value used downstream. No downstream unit's inputs depend on the new labels or the augmented series check.

## Side observations (non-blocking)

- SymPy's `Delta(Pe -> oo)` is printed in exponential form (`exp(2*alpha)` factors) rather than the `cosh/sinh` form of `Delta_inf_expected`; the assertion is made zero via `.rewrite(sp.exp)`. This is the documented deviation in the `## Applied: F2` block and is mathematically equivalent. The Mathematica engine handles this without rewriting because `FullSimplify` produces the cosh/sinh form directly.
- The orchestrator's preemptive `ConditionalExpression -> e` strip in `expectZero` (lines 22-27 of the `.wl`) is a generic idiom fix and does not affect any specific assertion in this unit.
- The SymPy docstring still says "Stage 41" (line 3) and the banner prints "STAGE 41" (line 33). The Mathematica banner says "STAGE 041". This is a pre-existing labeling mismatch with the unit-058 filename, not introduced by this batch; flagged for future cleanup, not a verification blocker.

## Verdict justification

All four directive findings are `resolved`: Mathematica is now an independent re-derivation (F1), the bracket and sharp-bottom-endpoint claims have substantive assertions in both engines (F2), the tautological `Delta_inf formula` check is correctly demoted to a sanity label (F3), and the analyticity-guaranteed constant-term check is augmented with a genuine non-vanishing first-order coefficient check (F4). Both engines exit 0, all new assertion labels appear and pass, output txts are fresher than their sources, and no regressions are visible — `material_change: false` because all derived quantities preserve their closed forms.
