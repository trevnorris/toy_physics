---
unit_id: 107
batch: IV.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-06T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 107

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
The `.wl` route to the normalized deformed branch `Y_def` / branch coefficients
was re-authored away from the shared series-on-ratio primitive.

- REMOVED (old `.wl` L41): `yDirect = Expand[Normal[Series[l0/lambdaDef, {z, 0, 5}]]]`
  — the identical operator the `.py` uses (`sp.series(L0/Lambda_def, z, 0, 6)`, py L42).
- ADDED (new `.wl` L41–53): an **order-by-order polynomial inversion of the operator
  identity `Lambda_def · Y = L0`**. Concretely:
  - `branchJet = ell0 + ell2 z² + ell4 z⁴ + I ell5 z⁵` (symbolic deformed-branch jet),
  - `branchAnsatz = 1 + y2 z² + y4 z⁴ + I y5 z⁵` (undetermined normalized coefficients),
  - `branchResidual = Expand[branchJet*branchAnsatz - ell0]`,
  - `branchEquations = Thread[(CoefficientList[branchResidual, z][[#+1]] & /@ {2,4,5}) == 0]`,
  - `branchSol = Solve[branchEquations, {y2, y4, y5}]` with a unique-solution guard (L47).
  The branch coefficients `branchY2/branchY4/branchY5` are then read off the solve and
  substituted `ell→l` (`branchSubstitution`, L49).
- Downstream consumers rewired to the new route (diff L42–44):
  `m2 = branchY2 /. branchSubstitution`, `m4 = branchY4 /. branchSubstitution`
  (`branchY4 = branchY2² − ell4/ell0`), `chiQ = (branchY5 /. branchSubstitution)/(1/27)`.
  Previously these were direct re-types of `-l2/l0`, `l2²/l0²−l4/l0`, `(-l5/l0)/(1/27)`
  mirroring py L47–49.
- `Clear[...]` extended to the new helper symbols (diff L9–10); `yFormula` (the hand
  closed form) retained unchanged as the falsifiable comparison target at L64.

**Assessment:**
Genuinely independent. The `.py` obtains the normalized coefficients by expanding the
rational function `L0/Lambda_def` as a power series (a black-box series operator). The new
`.wl` instead posits an undetermined-coefficient ansatz and solves the linear system that
`Lambda_def · Y = L0` imposes order-by-order — a structurally different native primitive
(`CoefficientList`+`Solve` on an undetermined ansatz, NOT `Series` on the shared ratio),
with distinct intermediate structure (`branchJet/branchAnsatz/branchResidual/branchEquations`
have no `.py` analogue). This is the precise route the directive sanctioned ("polynomial
inversion of `Lambda_def*Y=L0` order-by-order", the stage-105/109 pattern). Grep confirms
`Series`/`series` appears NOWHERE in the `.wl`. The shared-black-box independence defect is
closed.

Non-tautology preserved/strengthened: the `expectZero["normalized expansion direct-formula",
yDirect − yFormula]` assertion (L64) now compares the **solve-derived** ansatz `yDirect`
against the **independently hand-typed** closed form `yFormula` — it can still fail (a wrong
sign in either route breaks it), and it now also cross-checks the new primitive against the
closed form. `chi_Q` flows from `branchY5` (the solved odd coefficient), not a re-typed
`-l5/l0`. No collateral edits beyond the route swap and the `Clear` extension; the `.py`
reference engine was untouched.

## Exec log assessment

**SymPy:** exit=0. Reference engine unchanged. Notable lines:
- `chi_Q = 3*(S*beta**5 + 9*Sigma5)/(3*S - Sigma0)`
- `Sigma2 exact formula = 0`, `Sigma4 exact formula = 0`,
  `chi_Q - 3(S beta^5 + 9 Sigma5)/(3S - Sigma0) = 0`

**Mathematica:** exit=0. All four checks PASS on the new route:
- `normalized expansion direct-formula = 0` → `PASS`
- `m2 = -((sigma2 + (beta^2*sNorm)/3)/(sigma0 - 3*sNorm))`
- `chi_Q = (-3*(9*sigma5 + beta^5*sNorm))/(sigma0 - 3*sNorm)` (≡ SymPy form, num/denom both negated)
- `Sigma2 exact formula = 0` PASS, `Sigma4 exact formula = 0` PASS,
  `chi_Q - 3(sNorm beta^5 + 9 sigma5)/(3 sNorm - sigma0) = 0` PASS
- `Stage 107 Mathematica audit passed.`

No-regression: `Y_def`, `m2`, `m4`, `chi_Q`, the `Sigma2/Sigma4` exact formulas, and the
general `chi_Q` factorization are all still asserted and all still pass with the SAME emitted
values as HEAD.

**Output freshness:** The committed
`mathematica/output/...stage107..._audit.txt` has NO diff vs HEAD (byte-identical, method-only
change — orchestrator-confirmed). Its mtime (1780732640) is newer than the `.wl` mtime
(1780732403), confirming regeneration post-fix. The committed `.txt` deliverable body
(`L0…L5`, `Y_def`, `m2`, `m4`, `chi_Q`, `Sigma2/4_evenmatch`, all PASS lines) matches the
fresh exec log line-for-line.

## Material-change assessment

`material_change`: false. This is a method-only re-author of the second engine. No derived
result changed — the emitted MMA output is byte-identical to HEAD, the SymPy reference engine
is untouched, and every constant/formula (`L0…L5`, `Σ2`, `Σ4`, `chi_Q`) is unchanged. No
downstream unit can be affected.

## Side observations (non-blocking)

- `branchY4` is reconstructed as `branchY2² − ell4/ell0` (L51) rather than read directly from
  the solve's `y4`; both equal `(l2² − l0 l4)/l0²` and the `yDirect − yFormula` assertion (L64)
  ties the full solved ansatz to the closed form, so the z⁴ coefficient is still independently
  cross-checked. Not a defect; noting for completeness.
- The unique-solution guard on the new branch solve (L47, `Length[branchSol] =!= 1 → fail`)
  is a good addition that makes the polynomial-inversion route self-validating.

## Verdict justification

The sole finding (F1, mathematica_transliteration) is fully resolved. The `.wl` no longer
performs `Series[l0/lambdaDef, …]` on the normalized ratio (grep-confirmed absent); it reaches
`Y_def`, `m2`, `m4`, and `chi_Q` via an order-by-order undetermined-coefficient solve of
`Lambda_def·Y = L0` — a structurally distinct native primitive, not a renamed transliteration of
the `.py`'s series expansion. All assertions still pass (mathematica exit 0, sympy exit 0), all
emitted values are unchanged, and the committed Mathematica output is byte-identical to HEAD.
Verdict: **verified**, `material_change: false`.
