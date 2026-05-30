---
unit_id: 143
batch: IV.5
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T21:40:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 143

This is an integrity-remediation pass. The authoritative source of findings is the Codex
read-only review (`redteam/codex_reviews/stage_143.md`, verdict FINDINGS, 2 findings: R1
SymPy, R2 Mathematica). The directive `redteam/directives/stage_143.md` encodes those two
findings and carries Codex `## Applied:` blocks for both. I verified each is genuinely
CLOSED, that the reconciled PASS items were not disturbed, and that the new assertions are
non-tautological global positivity proofs (not the old cubic-coefficient corroboration).
This overwrites the prior 2026-05-27 verification (which graded the older 3-finding v2
directive).

## Per-finding outcomes

### R1 — insufficient_verification (SymPy)

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.py:52-61`.
The old single line `expect_equal("exp remainder leading term is Pi**3/6",
exp_rem_series.coeff(Pi, 3), sp.Rational(1, 6))` (plus its `sp.series(...)` setup) is gone.
It is replaced by a Taylor-remainder / monotonicity proof:
```
R_rem = sp.exp(Pi) - 1 - Pi - Pi**2/2
expect_equal("exp remainder R(0) == 0", R_rem.subs(Pi, 0), sp.Integer(0))
expect_equal("exp remainder R'(0) == 0", sp.diff(R_rem, Pi).subs(Pi, 0), sp.Integer(0))
expect_equal("exp remainder R''(0) == 0", sp.diff(R_rem, Pi, 2).subs(Pi, 0), sp.Integer(0))
expect_zero("exp remainder R'''(Pi) - exp(Pi) == 0", sp.diff(R_rem, Pi, 3) - sp.exp(Pi))
expect_positive("exp remainder R'''(Pi) = exp(Pi) > 0 for Pi>0", sp.exp(Pi))
```
The edit is confined to the named range; no collateral edits (helpers at 13-29, constant
pieces at 50-51, endpoint/constant gates at 63-86 all unchanged). I confirmed via grep that
no `coeff(Pi, 3)` / `leading term` text survives in the live `.py` or its committed `.txt`.

**Assessment:**
Correct and rigorous. The five assertions establish R(0)=R'(0)=R''(0)=0 and the closed-form
identity R'''(Pi) ≡ exp(Pi), with exp(Pi) > 0; by repeated FTC from 0 this forces R > 0 on
Pi>0 — the standard positive-third-derivative Taylor argument, a genuine GLOBAL proof rather
than a single-coefficient match. It is non-tautological: the comparison targets are
independent closed forms (0, 0, 0, exp(Pi)), not the decomposition primitive `num`/`decomp`
that the remainder feeds into. It would FAIL for the review's counterexample
`Pi**3/6 - Pi**4`: its third derivative is `1 - 24*Pi`, so `expect_zero("... R'''(Pi) -
exp(Pi) == 0", ...)` carries the nonzero residual `1 - 24*Pi - exp(Pi)` → AssertionError
(and `expect_positive(1 - 24*Pi)` is negative for Pi > 1/24, so it cannot conclude
`is_positive is True`). `expect_positive(sp.exp(Pi))` succeeds because `Pi` is declared
`positive=True` at line 31, making `sp.exp(Pi).is_positive is True`. Resolved.

### R2 — insufficient_verification (Mathematica)

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl:73-88`.
The old `expRemSeries = Normal[Series[...]]; expectEqual["exp remainder leading term is
piM^3/6", Coefficient[expRemSeries, piM, 3], 1/6]` is gone. It is replaced by BOTH a
real-domain `Reduce` proof and the same Taylor-remainder derivative backing:
```
expRemReduce = Reduce[Exp[piM] - 1 - piM - piM^2/2 > 0, piM, Reals] /.
  {(Element[piM, Reals] && piM > 0) -> True, (piM > 0) -> True};
... If[TrueQ[Simplify[expRemReduce === True || expRemReduce === (piM > 0)]], pass[...], fail[...]];
rRem = Exp[piM] - 1 - piM - piM^2/2;
expectEqual["exp remainder R(0) == 0", (rRem /. piM -> 0), 0];   (and R'(0), R''(0))
expectEqual["exp remainder R'''(piM) - Exp[piM] == 0", D[rRem, {piM, 3}], Exp[piM]];
expectPositive["exp remainder R'''(piM) = Exp[piM] > 0 for piM>0", Exp[piM]];
```

**Assessment:**
Correct and rigorous, and it is a SEPARATE Reduce from the pre-existing
`Reduce[num > 0, piM, Reals]`. I confirmed there are exactly two `Reduce[...]` calls in the
file: line 61 (the kept `num > 0` proof, review row 6) and line 75 (the new exp-remainder
proof). The wl:60-65 `num > 0` block is byte-for-byte the reconciled PASS item and was not
modified. The output shows `Reduce[exp remainder > 0, piM, Reals] = True` followed by
`PASS: exp remainder > 0 for piM > 0 via Reduce`, i.e. the Reduce actually reduced to True
(not an unresolved Inactive/interval form). The Taylor-remainder backing block is
independently rigorous and non-tautological for the same reasons as R1: derivatives of
`rRem` are compared against independent closed forms (0, 0, 0, Exp[piM]), not against
`num`/`decomp`. It would FAIL for `piM^3/6 - piM^4`: the Reduce would return a bounded
interval (≈ 0 < piM < 1/6) failing the `=== True || === (piM > 0)` test, and
`D[piM^3/6 - piM^4, {piM,3}] = 1 - 24*piM ≠ Exp[piM]` fails the `expectEqual` residual check.
Resolved.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `numerator - exact positive decomposition = 0` (decomposition lemma intact)
- `exp remainder R'''(Pi) - exp(Pi) == 0 = 0` and `exp remainder R'''(Pi) = exp(Pi) > 0 for Pi>0: exp(Pi)` (new global proof passes)
- `exp remainder R(0) == 0 / R'(0) == 0 / R''(0) == 0: lhs - rhs = 0` (derivative-at-0 chain)
No AssertionError; log footer `# exit_code: 0`.

**Mathematica:** exit=0. Notable lines:
- `Reduce[num > 0, piM, Reals] = True` → `PASS: num > 0 for piM > 0 via Reduce` (kept row-6 proof, untouched)
- `Reduce[exp remainder > 0, piM, Reals] = True` → `PASS: exp remainder > 0 for piM > 0 via Reduce` (NEW separate Reduce, genuinely reduced to True)
- `exp remainder R'''(piM) - Exp[piM] == 0: lhs - rhs = 0` → PASS (Taylor backing)
`Stage 143 Mathematica audit passed.` / log footer `# exit_code: 0`.

PASS-line counts: SymPy committed output shows 11 passing assertion lines (8 `lhs - rhs = 0`
equality lines including the new R(0)/R'(0)/R''(0) chain, 1 `expect_zero ... = 0`, 2
`expect_positive` value lines for the constant pieces + 1 for exp(Pi)), none raising
AssertionError, exit 0. Mathematica committed output shows 15 `PASS:` lines plus the two
`Reduce[...] = True` echoes, exit 0. Both match the directive's expected new PASS lines.

**Output freshness:** confirmed re-generated post-fix. Script mtimes are 2026-05-29
21:29:44 (both `.py` and `.wl`); committed outputs are newer — SymPy `.txt` 21:33:12,
Mathematica `.txt` 21:33:58. Exec logs are dated 21:32:56 (sympy) / 21:33:13 (mathematica),
consistent with the post-edit re-run.

## Material-change assessment

`material_change`: false. This is a verification-surface strengthening only. The carried
numeric results are unchanged: `R_infty = 0.145454452260420...` and `lim That/sqrt(Pi) =
0.725669130700713...` are identical in both engines' outputs to their pre-fix values, and
the decomposition identity / endpoint limits / limiting constants all still pass. No derived
constant moved; no downstream unit dependency is perturbed.

## Side observations (non-blocking)

- The SymPy `expect_positive` helper accepts a value via either `is_positive is True` or a
  numeric float check; for `pi**2 - 2*pi` and `pi**2/2 - 4` it resolves through the numeric
  branch (both are concrete `sp.pi`-based constants). This is pre-existing reconciled
  behavior (rows 1-2), not part of this fix.
- The ledger banner text prints `2/pi` in prose while the audit's stand-in symbol for the
  physical π is `Pi` and `piM` is the integration variable. Cosmetic, pre-existing, outside
  scope.

## Verdict justification

Both findings the Codex review raised — the exp-remainder positivity being gated only by the
cubic Taylor coefficient (which passes for the negative-going counterexample `Pi^3/6 -
Pi^4`) — are genuinely closed. In both engines the old single-coefficient line is removed and
replaced with a real global positivity proof: a Taylor-remainder/monotonicity argument
(R(0)=R'(0)=R''(0)=0, R'''≡exp>0) in SymPy, and that same backing plus a separate
`Reduce[Exp[piM]-1-piM-piM^2/2 > 0, piM, Reals] = True` in Mathematica. The new assertions
compare against independent closed forms (0/exp), are demonstrably wrong-input sensitive, and
the pre-existing `Reduce[num>0]` (row 6) and all reconciled PASS items (rows 1,2,4,5,6,7,9,10)
are untouched. Both scripts exit 0, outputs are fresh, and no derived numeric changed.
Verdict: verified.
