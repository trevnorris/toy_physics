---
unit_id: 143
batch: IV.5
created_at: 2026-05-29T00:00:00Z
supersedes: 2026-05-27T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-30T03:30:36Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 143 (rewrite encoding codex_review)

## What this rewrite does + why

The 2026-05-27 original directive (3 findings) predates the Codex read-only review
(`redteam/codex_reviews/stage_143.md`, verdict FINDINGS, 2 findings). That review found
that the exponential-remainder positivity is **not actually gated**: both engines only
check the cubic Taylor coefficient `coeff == 1/6`, which PASSES for a wrong remainder such
as `Pi^3/6 - Pi^4` that goes negative at large `Pi`. It does NOT prove
`exp(Pi) - 1 - Pi - Pi^2/2 > 0` for all finite `Pi > 0` — the inequality the paper claim
(`g_Pi < 1` decomposition lemma) actually rests on.

This rewrite encodes exactly the 2 review findings (R1 SymPy, R2 Mathematica). Every other
edit the original directive ordered (F1 endpoint/constant gates, F2 hardcoded-`sInf`
rework, F3 independent `Reduce[num>0]` proof) was already applied and the Codex review
graded all of those PASS — they are listed under "Reconcile (tainted-applied, do not
touch)" below and MUST NOT be re-edited. Only the leading-coefficient-only remainder check
is replaced, in each engine, with a genuine global positivity proof on `Pi > 0`.

Apply each finding below in order. After applying, append an `## Applied: R<n>` block under
that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append
`## Blocked: R<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line
ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>`
for Mathematica, each under `timeout 600`) and iterate until they exit 0 with all in-file
checks passing. Getting the scripts to run cleanly is your job; the orchestrator
independently re-runs afterward. If the Mathematica `Reduce` cannot resolve the
transcendental within the cap, do NOT raise the cap — reformulate the math using the
Taylor-remainder fallback documented in R2 below (the same argument R1 uses).

Do NOT touch paper.tex, notes/, or any prose documents.

---

## R1 — insufficient_verification (SymPy)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.py:52-54`

**Current code (confirmed live, lines 52-54):**
```python
# exp-remainder positivity via Taylor coefficient
exp_rem_series = sp.series(sp.exp(Pi) - 1 - Pi - Pi**2/2, Pi, 0, 5).removeO()
expect_equal("exp remainder leading term is Pi**3/6", exp_rem_series.coeff(Pi, 3), sp.Rational(1, 6))
```

**Issue:** `expect_equal(..., exp_rem_series.coeff(Pi, 3), sp.Rational(1, 6))` gates only the
cubic Taylor coefficient. As the review demonstrated, a wrong remainder like
`Pi**3/6 - Pi**4` has the same cubic coefficient and still PASSES, yet is negative for large
`Pi`. This does not prove `exp(Pi) - 1 - Pi - Pi^2/2 > 0` on `Pi > 0`, which is the positive
piece the `g_Pi < 1` decomposition lemma depends on.

**Required change — replace the 3 lines (52-54) with a genuine global positivity proof via
the Taylor-remainder / monotonicity argument:**

```python
# exp-remainder positivity: prove R(Pi) = exp(Pi) - 1 - Pi - Pi^2/2 > 0 for all Pi > 0.
# Rigorous argument (each line is an independent, can-fail assertion):
#   R(0) = 0, R'(0) = 0, R''(0) = 0, and R'''(Pi) = exp(Pi) > 0 for all Pi,
# so R is strictly increasing from 0 on Pi > 0, hence R(Pi) > 0 there.
R_rem = sp.exp(Pi) - 1 - Pi - Pi**2/2
expect_equal("exp remainder R(0) == 0", R_rem.subs(Pi, 0), sp.Integer(0))
expect_equal("exp remainder R'(0) == 0", sp.diff(R_rem, Pi).subs(Pi, 0), sp.Integer(0))
expect_equal("exp remainder R''(0) == 0", sp.diff(R_rem, Pi, 2).subs(Pi, 0), sp.Integer(0))
expect_zero("exp remainder R'''(Pi) - exp(Pi) == 0", sp.diff(R_rem, Pi, 3) - sp.exp(Pi))
expect_positive("exp remainder R'''(Pi) = exp(Pi) > 0 for Pi>0", sp.exp(Pi))
```

**Why this is a genuine global proof and CAN fail:**
- It establishes `R > 0` on `Pi > 0` rigorously: the three zero-derivative-at-0 facts plus
  `R''' = exp(Pi) > 0` everywhere force `R'' > 0`, then `R' > 0`, then `R > 0` for `Pi > 0`
  (repeated FTC from 0). This is the standard Taylor-with-positive-third-derivative argument,
  not a single-coefficient corroboration.
- It CAN FAIL for the review's counterexample: substitute the wrong remainder
  `Pi**3/6 - Pi**4` and `R'''` becomes `1 - 24*Pi`, so `expect_zero("... R'''(Pi) - exp(Pi) == 0", ...)`
  yields a nonzero residual and `expect_positive` on `1 - 24*Pi` (which is negative for
  `Pi > 1/24`) cannot conclude `is_positive is True` → AssertionError. The check is wrong-input
  sensitive precisely where the old coefficient check was blind.
- `R_rem` is built directly from `exp(Pi) - 1 - Pi - Pi^2/2`; the comparison targets
  (`0`, `0`, `0`, `exp(Pi)`) are independent closed forms, not reuses of the decomposition
  primitive `num`/`decomp`. (`expect_positive` here is the existing helper at lines 25-29;
  `expect_positive(sp.exp(Pi))` succeeds because `Pi` is declared `positive=True`, so
  `sp.exp(Pi).is_positive is True`.)

**Expected new PASS lines (replacing the single `exp remainder leading term is Pi**3/6` line):**
```
exp remainder R(0) == 0: lhs - rhs = 0
exp remainder R'(0) == 0: lhs - rhs = 0
exp remainder R''(0) == 0: lhs - rhs = 0
exp remainder R'''(Pi) - exp(Pi) == 0 = 0
exp remainder R'''(Pi) = exp(Pi) > 0 for Pi>0: exp(Pi)
```
(no AssertionError; script exits 0).

## Applied: R1

- files_changed:
  - `scripts/moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.py`
- summary: Replaced the cubic-coefficient exp-remainder check with derivative and positivity assertions proving the remainder positive on `Pi > 0`.
- deviation: none

---

## R2 — insufficient_verification (Mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl:73-74`

**Current code (confirmed live, lines 73-74):**
```mathematica
expRemSeries = Normal[Series[Exp[piM] - 1 - piM - piM^2/2, {piM, 0, 4}]];
expectEqual["exp remainder leading term is piM^3/6", Coefficient[expRemSeries, piM, 3], 1/6];
```

**Issue:** Identical defect to R1. `Coefficient[expRemSeries, piM, 3] == 1/6` gates only the
cubic coefficient and PASSES for a wrong finite-`piM` remainder; it does not prove
`Exp[piM] - 1 - piM - piM^2/2 > 0` on `piM > 0`.

**Required change — replace the 2 lines (73-74) with a direct global positivity proof.**
Primary route is Mathematica `Reduce` over the reals (this stage already runs a successful
`Reduce[num > 0, piM, Reals]` at lines 60-65, so a transcendental `Reduce` is feasible here);
back it with the same Taylor-remainder derivative argument used in R1 so the gate is
rigorous even if `Reduce` returns a form needing normalization:

```mathematica
(* exp-remainder positivity: prove R(piM) = Exp[piM] - 1 - piM - piM^2/2 > 0 for piM > 0. *)
(* Primary route: Reduce over the reals must reduce the inequality to piM > 0. *)
expRemReduce = Reduce[Exp[piM] - 1 - piM - piM^2/2 > 0, piM, Reals] /.
  {(Element[piM, Reals] && piM > 0) -> True, (piM > 0) -> True};
Print["Reduce[exp remainder > 0, piM, Reals] = ", fmt[expRemReduce]];
If[TrueQ[Simplify[expRemReduce === True || expRemReduce === (piM > 0)]],
  pass["exp remainder > 0 for piM > 0 via Reduce"],
  fail["exp remainder > 0 for piM > 0 via Reduce", expRemReduce]];
(* Independent backing route: Taylor-remainder monotonicity. *)
(* R(0)=0, R'(0)=0, R''(0)=0, R'''(piM)=Exp[piM]>0 => R strictly increasing from 0 => R>0. *)
rRem = Exp[piM] - 1 - piM - piM^2/2;
expectEqual["exp remainder R(0) == 0", (rRem /. piM -> 0), 0];
expectEqual["exp remainder R'(0) == 0", (D[rRem, piM] /. piM -> 0), 0];
expectEqual["exp remainder R''(0) == 0", (D[rRem, {piM, 2}] /. piM -> 0), 0];
expectEqual["exp remainder R'''(piM) - Exp[piM] == 0", D[rRem, {piM, 3}], Exp[piM]];
expectPositive["exp remainder R'''(piM) = Exp[piM] > 0 for piM>0", Exp[piM]];
```

**Why this is a genuine global proof and CAN fail:**
- `Reduce[Exp[piM] - 1 - piM - piM^2/2 > 0, piM, Reals]` resolves the inequality globally
  over the reals; under `$Assumptions = piM > 0` it should reduce to `piM > 0` (or `True`
  after the rule rewrite). For the wrong remainder `piM^3/6 - piM^4`, `Reduce` returns a
  bounded interval (roughly `0 < piM < 1/6`), so the `TrueQ[... === True || ... === (piM > 0)]`
  test is FALSE → `fail` → `Exit[1]`. The Taylor-remainder block then fails independently:
  `D[piM^3/6 - piM^4, {piM,3}] = 1 - 24*piM ≠ Exp[piM]`, so `expectEqual["... R'''(piM) - Exp[piM] == 0", ...]`
  has a nonzero residual and fails. Both routes are wrong-input sensitive.
- The Taylor-remainder block reuses the existing `expectEqual`/`expectPositive` helpers
  (lines 33-43) and compares `rRem` derivatives against independent closed forms
  (`0`, `0`, `0`, `Exp[piM]`), not against the decomposition primitive `num`/`decomp` —
  so it does not reuse the quantity it proves. `expectPositive[Exp[piM]]` succeeds because
  `$Assumptions` carries `piM > 0`, so `Simplify[Exp[piM] > 0]` is `True`.
- Mathematica idioms baked in: the Rule substitutions are parenthesized
  (`(rRem /. piM -> 0)`, `(D[rRem, piM] /. piM -> 0)`); the `Reduce` output is normalized
  by the same `/. {(...) -> True, (piM > 0) -> True}` rewrite the script already uses at
  line 61 (strips the `Element`/`ConditionalExpression`-style wrapper); no `*)` substrings
  appear inside any comment body (ASCII names `R`, `rRem` only).

**Expected new PASS lines (replacing the single `exp remainder leading term is piM^3/6` line):**
```
Reduce[exp remainder > 0, piM, Reals] = True            (* or piM > 0 *)
PASS: exp remainder > 0 for piM > 0 via Reduce
exp remainder R(0) == 0: lhs - rhs = 0
PASS: exp remainder R(0) == 0
exp remainder R'(0) == 0: lhs - rhs = 0
PASS: exp remainder R'(0) == 0
exp remainder R''(0) == 0: lhs - rhs = 0
PASS: exp remainder R''(0) == 0
exp remainder R'''(piM) - Exp[piM] == 0: lhs - rhs = 0
PASS: exp remainder R'''(piM) - Exp[piM] == 0
exp remainder R'''(piM) = Exp[piM] > 0 for piM>0: E^piM
PASS: exp remainder R'''(piM) = Exp[piM] > 0 for piM>0
```
(script must still exit 0).

## Applied: R2

- files_changed:
  - `mathematica/moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl`
- summary: Replaced the cubic-coefficient exp-remainder check with a real-domain `Reduce` proof and independent Taylor-remainder derivative checks.
- deviation: none

---

## Reconcile (tainted-applied, do NOT touch)

These items were applied orchestrator-direct under the 2026-05-27 directive and graded PASS
by the Codex review. They are tainted-applied (carried forward as accepted) and MUST NOT be
re-edited by this directive:

- SymPy `pi**2 - 2*pi > 0` (`...sympy_audit.py:50`) — review row 1, PASS.
- SymPy `pi**2/2 - 4 > 0` (`...sympy_audit.py:51`) — review row 2, PASS.
- SymPy endpoint limits `g_0`, `g_infty` (`...sympy_audit.py:61-62`) — review row 4, PASS.
- SymPy constants `R_infty`, `S_infty`, `That/sqrt(Pi)` (`...sympy_audit.py:77-79`) — review row 5, PASS.
- Mathematica independent `Reduce[num > 0, piM, Reals]` proof of `g_Pi < 1`
  (`...mathematica_audit.wl:60-65`) — review row 6, PASS. **Do NOT modify this `Reduce` block;
  R2 adds a SEPARATE `Reduce` for the exp-remainder.**
- Mathematica constant positivity pieces `Pi^2 - 2*Pi`, `Pi^2/2 - 4`
  (`...mathematica_audit.wl:71-72`) — review row 7, PASS.
- Mathematica endpoint limits `g_0`, `g_infty` (`...mathematica_audit.wl:82-83`) — review row 9, PASS.
- Mathematica dynamic limits / limiting constants `rInf`, `sInf`, `sigmaRatio`, `tHatRatio`
  (`...mathematica_audit.wl:91-102`, the hardcoded-`sInf` rework) — review row 10, PASS.

Edit ONLY the exp-remainder lines named in R1 (py:52-54) and R2 (wl:73-74). Leave the
helper definitions, banners, and every reconciled block byte-for-byte unchanged.

## Anti-tautology guard

Neither R1 nor R2 may compare the exp-remainder against any expression derived from the
same `num`/`decomp`/`gPi` primitive it is meant to certify. The positivity targets are the
independent closed forms `0, 0, 0, exp(Pi)` (derivative-at-0 and third-derivative identities)
and the standalone `Reduce` over the reals. A check is only valid if it FAILS for the
review's counterexample remainder `Pi^3/6 - Pi^4` (it does, via the `R'''` mismatch and the
bounded-interval `Reduce` result) and PASSES for the true `exp(Pi) - 1 - Pi - Pi^2/2`. Do
NOT fabricate any numeric literal; the only constants introduced (`0` and `1/6`-free) are the
exact derivative values of `R` at `0`, which the engine itself computes.

## For orchestrator/Codex consult

- **Feasibility of the Mathematica `Reduce` (R2 primary route):** the stage already runs
  `Reduce[num > 0, piM, Reals]` successfully at wl:60-65, and `num` contains `Exp[piM]`, so a
  transcendental `Reduce` is known-feasible on this machine within `timeout 600`. The new
  `Reduce[Exp[piM] - 1 - piM - piM^2/2 > 0, piM, Reals]` is a strictly simpler transcendental
  expression than `num`, so it should resolve at least as quickly. **Risk (low):** if `Reduce`
  returns an `Inactive`/interval form that the `/. {... -> True, (piM > 0) -> True}` rewrite
  does not normalize to `True`/`piM > 0`, the primary route's `If` will `fail`. **Mitigation
  already in the directive:** the Taylor-remainder backing block is fully self-contained and
  rigorous on its own; if `Reduce` proves troublesome, Codex may keep ONLY the Taylor-remainder
  block (drop the `Reduce` lines) — that still constitutes a genuine global positivity proof and
  satisfies R2. Codex must NOT raise the `timeout 600` cap; reformulate (i.e., drop to the
  Taylor route) instead.
- **SymPy `expect_positive(sp.exp(Pi))`:** relies on `Pi` being declared `positive=True`
  (live at `...sympy_audit.py:31`), which makes `sp.exp(Pi).is_positive is True`. Confirmed
  against the live declaration; no risk.
