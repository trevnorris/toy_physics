---
unit_id: 209
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 209

## Per-finding outcomes

### F1 — missing_verification_script (subtype: missing_mathematica)

**Classification:** resolved

**What changed:**
Codex authored a new Mathematica audit at the registered path
`mathematica/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_mathematica_audit.wl`
(226 lines, untracked/new — confirmed via `git status`). It follows the house
pattern: header comment (`:1`), `expectZero`/`expectTrue`/`fail` helpers that
`Exit[1]` on a nonzero/non-True residual (`:23-54`), a `stripConditional` guard
per house idiom (`:21`), and a final `STAGE 209 MATHEMATICA AUDIT PASSED` banner
with `Exit[0]` (`:224-225`). The exec log shows it appears, exits 0, with 18 PASS
lines and 0 FAIL.

**Assessment — manifest coverage.** All seven manifest items are present as
explicit, hard-failing assertions, none tautological:
- M1 `closureRoot[r] - coefficientRoot[r]` (`:102-105`)
- M2 (a) `abcByCoefficientExtraction - abcClosedForm` and (b) the polynomial
  residual `den[r](raySlope^2 - 2 H0 rayCurvature) - discFromCoefficients`
  (`:108-115`)
- M3 `closureRoot[r] - 2 H0/denominatorFunctional[r]` (`:118-121`)
- M4 non-triviality guard `D[denominatorFunctional[r],r] =!= 0` + scaled
  derivative numerator vs manifest N + derivative-law residual (`:140-154`)
- M5 degree-4 check on the resultant quartic + factorization identity +
  plus-radical-factor = derivative numerator (`:166-183`)
- M6 non-triviality guard `D[coefficientRoot[r],r] =!= 0`, ratio recovered by
  Solve, constant curvature, gradient-optimal stationarity (`:188-206`)
- M7 reciprocal invariance, differentiated reciprocal identity, equal-mix
  stationarity (`:211-222`)

**Assessment — genuine independence (load-bearing).** This is NOT a
transliteration of the SymPy `.py`; the `.wl` reaches the shared targets by
different Mathematica-native routes, matching the directive's anti-transliteration
guard:

1. **Quartic (M5), strongest evidence.** SymPy hand-squares:
   `Q = sp.expand(J**2 - 4*(kj-ki*r)**2*(A+B*r+C*r**2))` (`.py:97`).
   Mathematica eliminates the radical instead — exactly the directive's suggested
   alternative: `Resultant[polyPartFromDifferentiatedDiscriminant + 2(kj-ki r) z,
   z^2 - discFromCoefficients[r], z]` (`.wl:157-163`). Different algebraic engine
   (resultant elimination of an auxiliary `z` vs. literal product), same quartic.

2. **A/B/C discriminant reduction (M2).** SymPy subtracts the closed forms
   directly: `(1+r**2)*(kij**2 - 2*H0*kappa) - (A+B*r+C*r**2)` (`.py:69-72`).
   Mathematica instead *extracts* the coefficients with `CoefficientList[discPoly, r]`
   and compares the extracted `{c0,c1,c2}` against the closed forms
   `{ki^2-2H0u, 2ki kj-4H0v, kj^2-2H0w}` (`.wl:73-83,108-111`) — a coefficient-
   extraction decomposition, not the `.py`'s direct subtraction.

3. **Gradient ratio (M6).** SymPy substitutes the known answer `r_grad = kj/ki`
   and checks the derivative vanishes (`.py:115,122`). Mathematica additionally
   *solves* `Solve[D[raySlope[r],r]==0, r, Reals]` and checks `MemberQ[..., kj/ki]`
   (`.wl:187,192-198`) — recovering the optimizer rather than presupposing it.

The `.wl` also avoids the `.py`'s `kij -> kappa -> tau -> A,B,C -> S -> Phi -> N
-> J -> Q` variable choreography, defining instead `den/raySlope/rayCurvature/
discNumerator/closureRoot/coefficientRoot/denominatorFunctional` as functions of a
dummy `x`. Both M4-style non-triviality self-tests required by the directive are
present and pass (`M4 ... =!= 0 = True`, `M6 raw tau derivative ... =!= 0 = True`),
confirming the derivative checks are not trivially zero before reduction. No sign
assumption is imposed on the radicand: `$Assumptions` declares `u,v,w` only `Reals`
(`.wl:60-63`), keeping the identities branch-independent and matching the `.py`.

### F2 — symbol_assumption_error

**Classification:** resolved

**What changed:**
At the discriminant square root (`scripts/...sympy_audit.py:53-55`), Codex added a
three-line comment immediately above `S = sp.sqrt(A + B*r + C*r**2)`:
`# Radicand A + B r + C r^2 = Delta^sharp; admissible window requires` /
`# Delta^sharp >= 0 (notes sec. 3). Identities below are branch-independent,` /
`# so u, v, w are left unsigned.`
The diff (`stage_209_diff.patch`) confirms this is the ONLY change to the `.py` —
a pure 3-line comment insertion, no code touched.

**Assessment:**
Matches the directive's required change verbatim in substance (admissibility
premise `Delta^sharp >= 0`, notes §3 attribution, branch-independence note,
unsigned `u,v,w`). It is substantive, not vacuous: it records the domain premise,
the rationale for leaving `u,v,w` unsigned, and the branch-independence guarantee
that the Mathematica mirror reproduces. No `positive=True` was added to `u,v,w`
(`.py:43` unchanged), so the diagonal-neutral (`u=w=kappa_*`) and pair-symmetry
(`u=w`) reductions in §IV/§V are not over-constrained. No assertion was weakened or
changed — exec output is identical modulo the comment, and the script still exits 0
with every `expect_zero` reporting `= 0`.

## Exec log assessment

**SymPy:** exit=0. Notable lines: `explicit algebraic tau form = 0`,
`discriminant numerator reduction = 0`, `Phi derivative law = 0`,
`quartic degree minus 4 = 0`, `quartic factorization identity = 0`,
`gradient-optimal stationarity on diagonal-neutral branch = 0`,
`equal-mix stationarity on pair-symmetric branch = 0`,
`STAGE 192 SYMPY AUDIT COMPLETED SUCCESSFULLY`. (The "STAGE 192" banner is a
pre-existing internal-derivation-order label, already noted as benign by the
auditor; not part of either finding.)

**Mathematica:** exit=0, 18 PASS / 0 FAIL. Notable lines: `PASS: M1 closure root
minus coefficient-discriminant form`, `M2 discriminant coefficient extraction
residuals = {0, 0, 0}`, `PASS: M4 Phi derivative is not identically zero before
reductions` (= True), `M5 resultant quartic degree is 4 = True`, `PASS: M5
resultant factorization identity`, `PASS: M6 gradient ratio recovered from slope
stationarity` (= True), `PASS: M7 equal-mix stationarity`,
`STAGE 209 MATHEMATICA AUDIT PASSED`. Counted 18 PASS lines, matching the
orchestrator's independent re-run; M1-M7 all represented.

**Output freshness:** confirmed. Both saved `.txt` outputs have mtime 2026-06-02
10:54:33, newer than the `.wl` (10:47:39) and the `.py` (10:46:17) — re-generated
post-fix.

## Material-change assessment

`material_change`: false.

F1 is additive (a new second-engine verification of already-derived identities;
no result changed). F2 is a comment-only documentation edit. No derived result
that a downstream unit could depend on was altered.

## Side observations (non-blocking)

- The SymPy banner still self-labels "STAGE 192" (`.py:35,148`); the Mathematica
  `.wl` correctly says "STAGE 209". Cosmetic, pre-existing, outside both findings —
  flagged only, not blocking.

## Verdict justification

Both findings are resolved. F1's new `.wl` exists at the registered path, exits 0
with 18 PASS / 0 FAIL covering all seven manifest items as explicit hard-failing,
non-tautological checks, and is a genuinely independent derivation — it reaches the
quartic via `Resultant` (vs. the `.py`'s hand-squaring), the A/B/C reduction via
`CoefficientList` extraction (vs. direct subtraction), and the gradient ratio via
`Solve`+`MemberQ` (vs. substituting the known answer), with both required
non-triviality self-tests passing and no radicand sign assumption. F2 added the
exact required admissibility/branch-independence comment with no assertion change
and no over-constraining sign assumption. Diff is scoped to only the two intended
files, outputs are fresh, no regression. material_change is false.
