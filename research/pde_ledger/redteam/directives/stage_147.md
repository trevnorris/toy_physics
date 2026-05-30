---
unit_id: 147
batch: IV.5
created_at: 2026-05-29T00:00:00Z
supersedes: 2026-05-27T00:00:00Z
findings_count: 5
stop_cold: null
applied: true
applied_at: 2026-05-29T21:40:52-06:00
findings_applied: 5
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 147 (rewrite encoding Codex review of 2026-05-29)

## What this rewrite does

The 2026-05-27 directive (now superseded) added assertions that PASS but are
**tautological / self-consistent** rather than independent. A Codex read-only review
(`redteam/codex_reviews/stage_147.md`) found 5 such problems. This directive replaces
each disguised tautology with a check that uses a genuinely **independent primitive**
and can FAIL for a wrong input.

The stage's actual mathematical content lives in
`paper/appendices/stage_appendix_part04.tex` (the owning prose for the kernel claim;
the notes only reference the moments implicitly). The operative definitions are:

- Canonical source profile (appendix line ~798):
  `Sigma_*(x) = Pi_* e^(-Pi_* x) / (1 - e^(-Pi_*))`.
- Source moments / inner products (appendix `eq:app-part04-gbar-Sbar`, lines ~812-816):
  `gbar[sigma] = integral_0^1 sigma(x) c(x) dx`,
  `Sbar[sigma] = integral_0^1 sigma(x) Kq(x) dx`.
- Kernels (appendix `eq:app-part04-c-Kq`, lines ~818-822):
  `c(x) = cos(pi x / 2)`, `Kq(x) = cosh((pi/2)(1-x)) / cosh(pi/2)`.
- Normalization (appendix `eq:app-part04-positive-deformation-family`):
  `integral_0^1  smallsigma(x) dx = 1` and `integral_0^1 Sigma_*(x) dx = 1`,
  hence `integral_0^1 (smallsigma - Sigma_*) dx = 0` (the "integrates to zero" claim
  in notes line 85 that licenses dropping the centering constants).
- Retuning identity (appendix `eq:app-part04-deltaPi-firstorder`):
  `deltaPi = -eps (gbar_smallsigma - g_*) / gp_*`, so the traction shift coefficient
  is `A_T = -(dT_m/dPi) / (dg/dPi)` evaluated at `Pi_*`. (This is the relation R1
  must test independently; confirmed below against the closed form in py:33-38.)

So `g_*` and `S_*` are NOT primitively "the closed-form `gPi`/`Sformula` evaluated at
`Pi_*`" — they are the **moment integrals of `Sigma_*` against `c` and `Kq`**. The
closed forms `gPi`, `Sformula` are the *analytic evaluations* of those integrals. That
distinction is exactly what makes a quadrature route independent: it can catch a
transcription error in the closed forms.

Apply each finding below in order. After applying, append an `## Applied: R<n>` block
under that finding with: `files_changed`, `summary` (one sentence), and `deviation`
(or "none"). If a finding's required change is ambiguous or unsafe to apply
mechanically, append `## Blocked: R<n>` with a question instead — skip that finding,
continue with the rest.

Do NOT introduce features/refactors/stylistic changes beyond the named edits. Do NOT
touch paper.tex, notes/, or any prose document. After editing, RUN the affected
scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate
until they exit 0 with all in-file checks passing, all within `timeout 600` (a timeout
is a failure — reformulate the math, never raise the cap). Getting the scripts to run
cleanly is your job; the orchestrator independently re-runs afterward.

---

## R1 — tautological_check (sympy: A_T "chain route" rebuilds the same algebra)

**Target:** `scripts/moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py:66-75`

**Current code (tautological):**
```python
# --- Audit assertion: chain-rule consistency for A_T (independent derivation route) ---
# d T_m / d Sigma_0 = (1/2) * sqrt(9/(20 Sigma_*)) = 9/(40 T_*); cross-check the
# closed-form A_T against this differential identity assembled from scratch.
dTm_dSigma = sp.Rational(9, 40) / T_star
dSigma_dPi_at_star = 1/(1 - S_star/4) + Pi_star * Sp_star / (4*(1-S_star/4)**2)
AT_chain = sp.N(-dTm_dSigma * dSigma_dPi_at_star / gp_star, 30)
AT_30 = sp.N(AT, 30)
assert abs(AT_chain - AT_30) < sp.Float("1e-20", 30), \
    f"A_T chain-rule route disagrees with closed form: {AT_chain} vs {AT_30}"
print("PASS: A_T closed form agrees with chain-rule decomposition (residual < 1e-20)")
```

**Why it is tautological:** `dSigma_dPi_at_star` re-types by hand the exact factor
`1/(1-S/4) + Pi S'/(4(1-S/4)^2)` that the closed-form `AT` (py:33-38) already
multiplies by `9/(40 T_*)/g'_*`. Both routes share the same hand-written quotient and
the same `9/(40 T_*)` constant. A copied sign/power error survives in both → still PASS.

**Replacement (insert in place of lines 66-75):**
```python
# --- Audit assertion: A_T from automatic differentiation of T_m(Pi) (independent) ---
# Independent primitive: build T_m as a SYMBOLIC function of the free symbol Pi from
# Sformula (NOT from the cached numeric S_star/Sp_star/T_star), and let SymPy derive
# dT_m/dPi automatically. The retuning identity (appendix eq:app-part04-deltaPi-firstorder)
# gives A_T = -(dT_m/dPi)/(dg/dPi) at Pi_*. sp.diff regenerates the (1-S/4)^2 factor and
# the S' factor on its own, so a hand-written sign/power error in the closed-form AT
# (py:33-38) cannot be reproduced here.
Tm_of_Pi = sp.sqrt(sp.Rational(9, 20) * (Pi / (1 - Sformula/4)))
dTm_dPi = sp.diff(Tm_of_Pi, Pi)
dg_dPi = sp.diff(gPi, Pi)
AT_autodiff = sp.N(-(dTm_dPi.subs(Pi, Pi_star)) / (dg_dPi.subs(Pi, Pi_star)), 30)
AT_30 = sp.N(AT, 30)
assert abs(AT_autodiff - AT_30) < sp.Float("1e-20", 30), \
    f"A_T autodiff route disagrees with closed form: {AT_autodiff} vs {AT_30}"
print("PASS: A_T closed form agrees with autodiff of T_m(Pi) (residual < 1e-20)")
```

**Independent primitive + why no shared mistake:** SymPy's symbolic
differentiation (`sp.diff`) of the composite `T_m(Pi) = sqrt((9/20) Pi/(1-Sformula/4))`.
The target closed form `AT` is a hand-expanded chain rule with literal factors
`1/(1-S/4)`, `Pi S'/(4(1-S/4)^2)`, `9/(40 T_*)`. The autodiff route never types any of
those factors — it builds them from `Sformula` and the square-root structure. The only
shared inputs are `gPi`, `Sformula`, `Pi_star` (the upstream definitions, which are
under audit elsewhere), so a mistake localized to the hand-written `AT` expansion is
caught. Confidence: HIGH.

**Expected new PASS line:**
`PASS: A_T closed form agrees with autodiff of T_m(Pi) (residual < 1e-20)`

## Applied: R1

- files_changed:
  - `scripts/moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py`
- summary: Replaced the hand-written A_T chain-rule consistency check with SymPy autodifferentiation of T_m(Pi).
- deviation: none

---

## R2 — tautological_check (sympy: centered kernel checked against its own typed form)

**Target:** `scripts/moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py:88-99`

**Current code (tautological):**
```python
# --- Audit assertion: centered kernel structure matches the notes' boxed form ---
# Notes (section 2): W_*(x) = A_T (c(x) - g_*) + B_T (K_q(x) - S_*).
# g_* is the value of gFormula at Pi_*; S_* is Sformula(Pi_*). Verify the
# constant offsets in Wcenter equal A_T*(-gminus) + B_T*(-S_*) by extracting
# the constant term.
Wcenter_const = sp.simplify(Wcenter.subs([(x, sp.Symbol("__dummy"))]) -
                            (AT*c.subs(x, sp.Symbol("__dummy")) +
                             BT*Kq.subs(x, sp.Symbol("__dummy"))))
Wcenter_const_expected = sp.simplify(-AT*gminus - BT*Sformula.subs(Pi, Pi_star))
assert sp.simplify(Wcenter_const - Wcenter_const_expected) == 0, \
    f"Centered kernel constant offset mismatch: {Wcenter_const} vs {Wcenter_const_expected}"
print("PASS: Centered kernel W_*(x) has form A_T(c - g_*) + B_T(K_q - S_*)")
```

**Why it is tautological:** `Wcenter` (py:84) is literally defined as
`AT*(c-gminus) + BT*(Kq-Sformula(Pi_*))`. The check then subtracts `AT*c + BT*Kq` and
compares the remainder to `-AT*gminus - BT*S_*` — i.e. it re-derives the constant it
just typed in. It tests nothing about the actual rigidity-kernel *projection* claim.

**Replacement (insert in place of lines 88-99):** Test the REAL claim from notes
section 2 / appendix `eq:app-part04` — that projecting `W_*` against a deformation
reproduces the two-moment traction formula:
`integral_0^1 W_*(x)[smallsigma(x) - Sigma_*(x)] dx == A_T(gbar_s - g_*) + B_T(Sbar_s - S_*)`
for an arbitrary normalized deformation `smallsigma`.
```python
# --- Audit assertion: rigidity-kernel projection identity (independent quadrature) ---
# Notes sec.2 / appendix eq:app-part04-deltaT-firstorder: the traction shift equals the
# projection of W_*(x) against (smallsigma - Sigma_*). Test it for a concrete, normalized,
# NON-canonical deformation smallsigma(x) = 2x (integral_0^1 2x dx = 1, positive, and a
# different shape from Sigma_*). LHS by numerical quadrature of the kernel projection;
# RHS by the two-moment coefficient formula. The two sides share NO derivation step:
# LHS integrates W_* numerically, RHS uses the algebraic A_T, B_T and the moment
# integrals. This also exercises the "integrates to zero" centering claim, since the
# -g_*, -S_* constants inside W_* drop out only because integral(smallsigma-Sigma_*)=0.
Sigma_star_x = Pi_star * sp.exp(-Pi_star * x) / (1 - sp.exp(-Pi_star))
smallsigma_x = 2*x
# normalization sanity (must integrate to 1 each):
norm_s = sp.N(sp.integrate(smallsigma_x, (x, 0, 1)), 40)
norm_Sigma = sp.N(sp.integrate(Sigma_star_x, (x, 0, 1)), 40)
assert abs(norm_s - 1) < sp.Float("1e-30", 40) and abs(norm_Sigma - 1) < sp.Float("1e-30", 40), \
    f"deformation/source not normalized: {norm_s}, {norm_Sigma}"
Wstar_x = AT*(c - g_star) + BT*(Kq - S_star)
lhs_proj = sp.N(sp.integrate(Wstar_x * (smallsigma_x - Sigma_star_x), (x, 0, 1)), 30)
gbar_s = sp.N(sp.integrate(smallsigma_x * c, (x, 0, 1)), 40)
Sbar_s = sp.N(sp.integrate(smallsigma_x * Kq, (x, 0, 1)), 40)
rhs_moment = sp.N(AT*(gbar_s - g_star) + BT*(Sbar_s - S_star), 30)
assert abs(lhs_proj - rhs_moment) < sp.Float("1e-22", 30), \
    f"kernel projection != two-moment formula: {lhs_proj} vs {rhs_moment}"
print("PASS: kernel projection of W_* reproduces two-moment traction shift (residual < 1e-22)")
# --- Audit assertion: source-centering of the rigidity kernel (independent) ---
# CONSULT Q5 (batch 6): the projection identity above is BLIND to the centering
# constants -A_T*g_*, -B_T*S_* -- they vanish against (smallsigma - Sigma_*) because
# both integrate to 1. The x-independence check (R3) only proves the offset is CONSTANT,
# not that it equals -A_T*g_*-B_T*S_*. The kernel's defining centering condition is
# orthogonality to the canonical source: integral_0^1 Sigma_*(x) W_*(x) dx == 0. This
# DOES test the constants: dropping them leaves integral Sigma_*(A_T c + B_T Kq) =
# A_T*g_* + B_T*S_* != 0, so a missing-centering bug now fails here.
center_resid = sp.N(sp.integrate(Sigma_star_x * Wstar_x, (x, 0, 1)), 30)
assert abs(center_resid) < sp.Float("1e-22", 30), \
    f"kernel not centered against Sigma_*: integral Sigma_* W_* = {center_resid}"
print("PASS: rigidity kernel W_* is source-centered (integral Sigma_* W_* = 0, residual < 1e-22)")
```

**Independent primitive + why no shared mistake:** LHS uses **numerical quadrature**
(`sp.integrate`) of the kernel-times-deformation integral; RHS uses the **algebraic
two-moment coefficient formula** with `A_T`, `B_T` and the moment integrals of the test
deformation. The previous tautology compared `W_*` to its own typed constant; here
`W_*` is exercised through an integral against an independent test profile `2x`, and the
result must match a separately-assembled moment expression. A wrong kernel (wrong sign,
missing centering term, swapped `c`/`Kq`) breaks the equality. Note the residual `1e-22`
absorbs quadrature precision at 30 digits. Confidence: HIGH that this is non-tautological.
Residual risk noted in consult block: this is an analytic-integral identity, so if both
sides are evaluated symbolically by the same SymPy integrator it could in principle share
an integrator bug — mitigated because RHS does not integrate `W_*` at all (only the
simple products `2x*c`, `2x*Kq`), so the kernel structure is genuinely cross-checked.

**Expected new PASS lines:**
`PASS: kernel projection of W_* reproduces two-moment traction shift (residual < 1e-22)`
`PASS: rigidity kernel W_* is source-centered (integral Sigma_* W_* = 0, residual < 1e-22)`

Leave the `Wcenter = sp.simplify(...)` print (py:84-86) intact — it is print-only
output, not an assertion, so it is not part of the tautology.

## Applied: R2

- files_changed:
  - `scripts/moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py`
- summary: Replaced the centered-kernel constant tautology with projection and source-centering quadrature checks.
- deviation: none

---

## R3 — insufficient_verification (mathematica: kernel checked at one sample point x=1/2)

**Target:** `mathematica/moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl:87-91`

**Current code (insufficient + tautological):**
```mathematica
(* --- Audit assertion: centered kernel structure --- *)
wCenterConst = FullSimplify[(wCenter - (aT*c + bT*kq)) /. x -> 1/2];
wCenterConstExpected = FullSimplify[-aT*gMinus - bT*sStar];
expectZero["W_*(x) centered form A_T(c-g_*) + B_T(K_q-S_*)",
  Chop[wCenterConst - wCenterConstExpected, 10^-25]];
```

**Why it is insufficient:** it samples only `x -> 1/2` after `wCenter` (wl:84) was
typed as the target form `aT*(c-gMinus)+bT*(kq-sStar)`. A nonconstant erroneous residue
that happens to vanish at `x=1/2` would pass; and like R2 it re-derives a typed-in
constant rather than testing the projection claim.

**Replacement (insert in place of lines 87-91):** mirror R2's projection identity using
Mathematica's independent numerical integrator, AND add a full-symbolic-in-x check that
`W_* - (aT*c + bT*kq)` is genuinely x-independent (no single-point sampling).
```mathematica
(* --- Audit assertion: full symbolic x-independence of the centering offset --- *)
(* W_*(x) - (aT c(x) + bT Kq(x)) must be CONSTANT in x (the centering shift). Check the *)
(* symbolic derivative in x is identically zero -- holds for all x, not one sample. *)
wStar = aT*(c - gMinus) + bT*(kq - sStar);
expectZero["W_* centering offset is x-independent",
  FullSimplify[D[wStar - (aT*c + bT*kq), x]]];

(* --- Audit assertion: rigidity-kernel projection identity (independent quadrature) --- *)
(* Mirror of the SymPy R2 check but via NIntegrate (independent engine + numerical    *)
(* primitive). smallsigma(x) = 2 x is normalized and non-canonical. LHS integrates the *)
(* kernel against (smallsigma - Sigma_*); RHS is the algebraic two-moment formula.     *)
sigmaStarX = pStar*Exp[-pStar*x]/(1 - Exp[-pStar]);
smallSigmaX = 2*x;
normS = NIntegrate[smallSigmaX, {x, 0, 1}, WorkingPrecision -> 40];
normSigma = NIntegrate[sigmaStarX, {x, 0, 1}, WorkingPrecision -> 40];
expectZero["deformation normalized", Chop[(normS - 1) + (normSigma - 1), 10^-30]];
lhsProj = NIntegrate[wStar*(smallSigmaX - sigmaStarX), {x, 0, 1}, WorkingPrecision -> 40];
gBarS = NIntegrate[smallSigmaX*c, {x, 0, 1}, WorkingPrecision -> 40];
sBarS = NIntegrate[smallSigmaX*kq, {x, 0, 1}, WorkingPrecision -> 40];
rhsMoment = aT*(gBarS - gStar) + bT*(sBarS - sStar);
expectZero["kernel projection reproduces two-moment traction shift",
  If[Abs[lhsProj - rhsMoment] < 10^-20, 0, lhsProj - rhsMoment]];

(* --- Audit assertion: source-centering of the rigidity kernel (independent) --- *)
(* CONSULT Q5 (batch 6): the projection identity is blind to the centering constants *)
(* (they vanish against (smallSigma - sigmaStar)); the D[...,x] check only proves the *)
(* offset is constant in x. The kernel's centering condition is orthogonality to the  *)
(* canonical source: integral Sigma_* W_* == 0. Dropping the constants leaves          *)
(* A_T g_* + B_T S_* != 0, so this assertion DOES test them.                           *)
centerResid = NIntegrate[sigmaStarX*wStar, {x, 0, 1}, WorkingPrecision -> 40];
expectZero["rigidity kernel W_* is source-centered (integral Sigma_* W_* = 0)",
  If[Abs[centerResid] < 10^-20, 0, centerResid]];
```

**Independent primitive + why no shared mistake:** two primitives. (1) Symbolic
`D[..., x]` proving x-independence of the centering offset for ALL x (replaces the
single-point `/. x -> 1/2`). (2) `NIntegrate` quadrature projection identity vs the
algebraic moment formula — the same content as R2 but computed by Mathematica's
numerical integrator, an implementation fully independent of SymPy's. A wrong kernel or
a residue nonconstant in x now fails. Confidence: HIGH.

**Expected new PASS lines:**
`PASS: W_* centering offset is x-independent`
`PASS: kernel projection reproduces two-moment traction shift`
`PASS: rigidity kernel W_* is source-centered (integral Sigma_* W_* = 0)`

Leave the `wCenter = FullSimplify[...]` print (wl:84-85) intact (print-only).

## Applied: R3

- files_changed:
  - `mathematica/moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl`
- summary: Replaced the single-point centered-kernel check with x-independence, projection, and source-centering checks.
- deviation: Used `Chop` on the zero derivative residual, higher-precision `NIntegrate` for the centering integral, and parser-safe comment wording without loosening thresholds.

---

## R4 — tautological_check (sympy: g_*, S_* resubstitution just repeats lines 25-26)

**Target:** `scripts/moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py:101-113`

**Current code (tautological):**
```python
# --- Audit assertion: source-moment definitions of g_*, S_* reproduce gFormula(Pi_*), Sformula(Pi_*) ---
# In the notes' inner product (lines 96-105), g_* is identified with the value
# gFormula takes at the canonical Pi_*, and S_* with Sformula at Pi_*.
# Verify the script's evaluation of g_*, S_* matches the symbolic substitution
# to high precision (this guards against an accidental redefinition of either
# moment between the family-1 anchor block and the kernel-assembly block).
g_star_resub = sp.N(gPi.subs(Pi, Pi_star), 40)
S_star_resub = sp.N(Sformula.subs(Pi, Pi_star), 40)
assert abs(g_star_resub - g_star) < sp.Float("1e-30", 40), \
    f"g_* resubstitution drift: {g_star_resub} vs {g_star}"
assert abs(S_star_resub - S_star) < sp.Float("1e-30", 40), \
    f"S_* resubstitution drift: {S_star_resub} vs {S_star}"
print("PASS: g_*, S_* moment values stable across audit (drift < 1e-30)")
```

**Why it is tautological:** `gPi.subs(Pi, Pi_star)` and `Sformula.subs(Pi, Pi_star)`
are byte-for-byte the definitions of `g_star`/`S_star` at py:25-26. This only catches
mutation, never tests whether the closed forms ARE the source moments.

**Replacement (insert in place of lines 101-113):** compute the moments from their
INTEGRAL definitions (appendix `eq:app-part04-gbar-Sbar`) and confirm the closed forms
`gPi`, `Sformula` reproduce them at `Pi_*`.
```python
# --- Audit assertion: g_*, S_* equal their source-moment integrals (independent quadrature) ---
# Appendix eq:app-part04-gbar-Sbar/c-Kq define g_* = integral_0^1 Sigma_*(x) c(x) dx and
# S_* = integral_0^1 Sigma_*(x) Kq(x) dx, with Sigma_*(x) = Pi_* e^(-Pi_* x)/(1-e^(-Pi_*)).
# The closed forms gPi, Sformula are the ANALYTIC evaluations of those integrals. Compute
# the integrals directly by quadrature and compare to gPi(Pi_*), Sformula(Pi_*): a
# transcription error in gPi/Sformula would be caught (the old resub check could not).
Sigma_star_x = Pi_star * sp.exp(-Pi_star * x) / (1 - sp.exp(-Pi_star))
g_star_moment = sp.N(sp.integrate(Sigma_star_x * c, (x, 0, 1)), 40)
S_star_moment = sp.N(sp.integrate(Sigma_star_x * Kq, (x, 0, 1)), 40)
assert abs(g_star_moment - g_star) < sp.Float("1e-25", 40), \
    f"g_* moment integral != gPi(Pi_*): {g_star_moment} vs {g_star}"
assert abs(S_star_moment - S_star) < sp.Float("1e-25", 40), \
    f"S_* moment integral != Sformula(Pi_*): {S_star_moment} vs {S_star}"
print("PASS: g_*, S_* equal their source-moment integrals (residual < 1e-25)")
```

**Independent primitive + why no shared mistake:** **numerical/analytic quadrature** of
`integral Sigma_*(x) c(x) dx` and `integral Sigma_*(x) Kq(x) dx`, where `Sigma_*(x)` is
built directly from `Pi_star` per the appendix. The target `g_star`/`S_star` are the
closed-form algebraic expressions `gPi`/`Sformula`. The integral route never references
`gPi`/`Sformula`, so a typo in either closed form (e.g. wrong `(4 Pi^2 + pi^2)`
denominator) is detected. Confidence: HIGH. (If `sp.integrate` of `Sigma_star_x * Kq`
is slow, Codex may instead evaluate via `Sigma_star_x * Kq` sampled with
`sp.Integral(...).evalf()` or `mpmath.quad` — both are still the quadrature primitive,
independent of the closed form. Stay within `timeout 600`.)

**Expected new PASS line:**
`PASS: g_*, S_* equal their source-moment integrals (residual < 1e-25)`

## Applied: R4

- files_changed:
  - `scripts/moving_throat_pde_stage147_first_order_rigidity_kernel_sympy_audit.py`
- summary: Replaced the g_*/S_* resubstitution drift check with source-moment integral checks.
- deviation: none

---

## R5 — transliteration (mathematica block is a direct port of the SymPy block)

**Target:** `mathematica/moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl:74-79`
(the A_T chain-rule block) and `:93-99` (the g_*/S_* resub block).

**Current code (chain-rule port, tautological):**
```mathematica
(* --- Audit assertion: chain-rule consistency for A_T (independent route) --- *)
dTmDSigma = 9/(40*tStar);
dSigmaDPi = 1/(1 - sStar/4) + pStar*sPrimeStar/(4*(1 - sStar/4)^2);
aTChain = N[-dTmDSigma*dSigmaDPi/gPrimeStar, 30];
expectZero["A_T closed form vs chain-rule route",
  If[Abs[aTChain - aT] < 10^-20, 0, aTChain - aT]];
```

**Current code (resub port, tautological):**
```mathematica
(* --- Audit assertion: source-moment values g_*, S_* stable under resubstitution --- *)
gStarResub = N[gFormula /. p -> pStar, 40];
sStarResub = N[sFormula /. p -> pStar, 40];
expectZero["g_* resubstitution drift",
  If[Abs[gStarResub - gStar] < 10^-30, 0, gStarResub - gStar]];
expectZero["S_* resubstitution drift",
  If[Abs[sStarResub - sStar] < 10^-30, 0, sStarResub - sStar]];
```

**Why it is a transliteration:** both blocks copy the SymPy hand-written chain rule and
the resubstitution drift check verbatim — same literals, same primitive — so a wrong
target copied into both engines passes in both.

**Replacement (a): chain-rule block at wl:74-79** — use Mathematica's own symbolic
differentiator of `T_m(p)`, mirroring R1.
```mathematica
(* --- Audit assertion: A_T from symbolic differentiation of T_m(p) (independent) --- *)
(* Build T_m as a symbolic function of p from sFormula and let D[] derive dT_m/dp. *)
(* Retuning identity: A_T = -(dT_m/dp)/(dg/dp) at pStar. Independent of the hand-  *)
(* written closed-form aT (wl:49-52). *)
tmOfP = Sqrt[(9/20)*(p/(1 - sFormula/4))];
dTmDp = D[tmOfP, p];
dgDp = D[gFormula, p];
aTAutodiff = N[-(dTmDp /. p -> pStar)/(dgDp /. p -> pStar), 30];
expectZero["A_T closed form vs autodiff of T_m(p)",
  If[Abs[aTAutodiff - aT] < 10^-20, 0, aTAutodiff - aT]];
```

**Replacement (b): resub block at wl:93-99** — replace with the moment-integral
quadrature route, mirroring R4 via `NIntegrate`.
```mathematica
(* --- Audit assertion: g_*, S_* equal their source-moment integrals (NIntegrate) --- *)
(* Appendix eq:app-part04-gbar-Sbar: g_* = integral Sigma_*(x) c(x) dx, S_* = integral *)
(* Sigma_*(x) Kq(x) dx, Sigma_*(x) = pStar e^(-pStar x)/(1-e^(-pStar)). Compare to the  *)
(* closed forms gFormula, sFormula at pStar via an independent numerical integrator.    *)
sigmaStarXm = pStar*Exp[-pStar*x]/(1 - Exp[-pStar]);
gStarMoment = NIntegrate[sigmaStarXm*c, {x, 0, 1}, WorkingPrecision -> 40];
sStarMoment = NIntegrate[sigmaStarXm*kq, {x, 0, 1}, WorkingPrecision -> 40];
expectZero["g_* equals source-moment integral",
  If[Abs[gStarMoment - gStar] < 10^-20, 0, gStarMoment - gStar]];
expectZero["S_* equals source-moment integral",
  If[Abs[sStarMoment - sStar] < 10^-20, 0, sStarMoment - sStar]];
```

**Independent primitive + why no shared mistake:** Mathematica now uses (a) its own CAS
`D[]` of the composite `T_m(p)` — independent of the hand-written `aT` and an
independent implementation from SymPy's `sp.diff`; (b) `NIntegrate` of the moment
integrals — independent of the closed forms `gFormula`/`sFormula` and an independent
implementation from SymPy's integrator. The two engines no longer share the
derivation primitive with the target. The *strategy* (CAS-diff + quadrature) is shared
across engines by design — see consult block; the anti-tautology requirement is
independence from the **target's** primitive (the hand-typed closed form), which is met.
Confidence: HIGH.

**Expected new PASS lines:**
`PASS: A_T closed form vs autodiff of T_m(p)`
`PASS: g_* equals source-moment integral`
`PASS: S_* equals source-moment integral`

Mathematica idioms reminder: no `*)` inside comment text (use `pStar`, `sStar`, not
`S_*` literally where it could close a comment — here comment text uses prose so it is
safe, but verify on apply); the `expectZero` wrapper already does `FullSimplify[Together[
Expand[...]]]` then `=== 0`; the `If[Abs[...] < tol, 0, residual]` pattern keeps
numerical residues from failing the `=== 0` test. Parenthesize any inline `(expr /. v ->
val)` you add. `Chop[..., 10^-30]` near-zero residues before comparing if FindRoot/
NIntegrate noise appears.

## Applied: R5

- files_changed:
  - `mathematica/moving_throat_pde_stage147_first_order_rigidity_kernel_mathematica_audit.wl`
- summary: Replaced the Mathematica chain-rule and resubstitution transliterations with D[] autodiff and NIntegrate moment checks.
- deviation: none

---

## Reconcile (tainted-applied, KEEP — genuine can-fail anchors)

These three numeric anchors are real, can FAIL against hard-coded literals, and were
applied by the prior orchestrator-direct edit. Mark them tainted-applied and KEEP as-is:

1. **A_T vs paper literal** — py:53,56 / wl:64,67-68. Literal `-4.27263956256927`.
   Sourced from `paper/appendices/stage_appendix_part04.tex:846` and
   `paper/parts/part04_geometry_retarded_mouth.tex:860` (both verbatim). Also quoted
   rounded as `-4.27264` in `paper/stages/stage_147.tex:16`. KEEP.
2. **B_T vs paper literal** — py:54,59 / wl:65,69-70. Literal `0.134875005736706`.
   Sourced from `stage_appendix_part04.tex:848` and `part04_geometry_retarded_mouth.tex:862`
   (both verbatim). Also rounded as `0.134875` in `stage_147.tex:16`. KEEP.
3. **|A_T|/B_T vs literal 31.6785** — py:55,62 / wl:66,71-72. WARNING: `31.6785` is NOT
   a verbatim paper literal — it does not appear anywhere in `paper/` (the paper quotes
   only `A_T`, `B_T` at `stage_appendix_part04.tex:843-848`). It is the script's own
   computed ratio (`31.678512...`) checked to 1e-3. **R6 (consult Q6, CONCUR): the live
   script currently MISLABELS this as a paper-quoted value** at
   `...stage147...py:62-64` (and the Mathematica analogue at `...wl:71-72`). KEEP the
   check (it can-fail: a corrupted A_T or B_T moves the ratio out of 1e-3), but CORRECT
   the comment/print label so it does NOT claim paper-sourcing — describe it as the
   script's own computed `|A_T|/B_T` ratio cross-check (e.g. label it
   `|A_T|/B_T computed ratio (cross-check, not a paper literal)`). This is a
   script-comment/label edit only (rides this fix loop; touches no paper/notes file).
   Do NOT fabricate or move the numeric value; only fix the misleading wording.

Optional strengthening (NOT required, do only if trivial): the appendix also quotes
`g'_* ≈ 0.0714453558083195` (`stage_appendix_part04.tex:831`, `part04...tex:845`),
`Sigma_0(Pi_*) ≈ 1.80594111095636` (line ~773) and `T_m(Pi_*) ≈ 0.901484054174205`
(line ~775). These are additional genuine can-fail anchors against the script's
`gp_star`, `Sigma_star`, `T_star`. If adding them is mechanical, add `gp_star` vs
`0.0714453558083195` (1e-12), `Sigma_star` vs `1.80594111095636` (1e-12), `T_star` vs
`0.901484054174205` (1e-12) as extra `assert`/`expectZero` anchors. If not trivial,
skip — they are not part of the 5 findings.

---

## Anti-X-X guard (read before applying)

Every replacement above is built so the compared quantities use DIFFERENT primitives:

- R1/R5(a): target = hand-expanded closed form `AT`; check = CAS `diff`/`D` of the
  composite `T_m(Pi)`. The check never types the chain-rule factors.
- R2/R3(quadrature)/R4/R5(b): target = algebraic closed forms (`gPi`, `Sformula`,
  typed kernel `W_*`); check = quadrature of the moment/projection integrals. The
  integral route never references the closed forms it is testing.
- R3(symbolic): the x-independence check uses `D[..., x]` over ALL x — no single sample.

For EACH new check, the script must be able to FAIL if the target is wrong: corrupt
`AT` (R1/R5a), corrupt the kernel sign or drop a centering term (R2/R3), corrupt the
denominator of `gPi`/`Sformula` (R4/R5b) — all break their respective assertions. Do
NOT "fix" a failing new check by editing the check to match the (possibly wrong) target;
if a new check fails, STOP and append a `## Blocked: R<n>` question — a failure here is
a genuine signal, not a tolerance to loosen.

Do NOT reuse `gPi`/`Sformula`/`AT`/`Wcenter` as the source of truth for any check whose
job is to test those very quantities. The moment integrals are sourced from `Pi_star`
and the kernel definitions `c`, `Kq`, `Sigma_*(x)` only.

---

## Consult resolution (batch-6 Claude+Codex, `redteam/codex_reviews/_consult_batch6.md`)

The read-only consult ran BEFORE this fix loop. Outcome on stage 147:
- **Q4 CONCUR** — `stage_appendix_part04.tex:798-824` is the operative source; the moment/
  kernel/normalization defs and `A_T=-(dT_m/dPi)/(dg/dPi)` are faithful; `T_m`/`T_star`
  match (py:30-31).
- **Q5 DISPUTE → resolved** — the projection identity is valid and catches sign/c↔Kq
  swaps, but is BLIND to the centering constants (they vanish against `(σ-Σ_*)`), and the
  `D[...,x]==0` check only proves the offset is constant, not its value. **Added the
  source-centering assertion `∫ Σ_* W_* dx == 0`** to R2 (SymPy) and R3 (Mathematica) —
  it fails if the `-A_T g_* - B_T S_*` constants are dropped. Non-conceptual strengthening.
- **Q6 CONCUR** — autodiff route genuinely independent; **but** the live script mislabels
  the `31.6785` ratio as paper-quoted → R6 label-correction added to Reconcile item 3.

Routes flagged at draft time (now settled by the consult above):

1. **Inner-product / source-moment definition location (HIGH importance).** The stage
   *notes* do NOT define the inner product — they only say "two source moments" and
   "Sigma_eps - Sigma_* integrates to zero" (notes lines 85, 123). The operative
   definitions (`gbar[sigma]=int sigma c`, `Sbar[sigma]=int sigma Kq`, `Sigma_*(x)=
   Pi_* e^(-Pi_* x)/(1-e^(-Pi_*))`, `c`, `Kq`, normalization to 1) live in
   `paper/appendices/stage_appendix_part04.tex` (lines ~798-822). I used the APPENDIX as
   the authoritative source. Confirm the appendix is the correct owning prose for stage
   147's kernel math before applying R2/R4 (the directive cites appendix line numbers
   that should be re-confirmed on apply, as they may drift).

2. **R2/R4 SymPy integrator sharing a bug (MEDIUM).** R2 (projection) and R4 (moments)
   both use `sp.integrate`. If `sp.integrate` had a systematic error it could affect
   both. Mitigation in the directive: R2's RHS does not integrate the kernel `W_*`
   (only the simple `2x*c`, `2x*Kq` products), and R4 cross-checks against the
   *closed-form* `gPi`/`Sformula` which use no integrator at all. So a shared integrator
   bug would have to coincidentally agree with an independent closed form — unlikely but
   not impossible. If Codex sees `sp.integrate` hang or return unevaluated for
   `Sigma_*(x)*Kq`, switch to `mpmath.quad` (pure numerical quadrature, a different
   primitive) — that strictly improves independence. Flag if you do.

3. **Cross-engine STRATEGY sharing (LOW, by design).** R5 has Mathematica mirror SymPy's
   strategy (CAS-diff for A_T, quadrature for moments/projection). The *implementations*
   are independent (Wolfram `D`/`NIntegrate` vs SymPy `diff`/`integrate`), and both are
   independent of the *target* (hand-written closed forms). I judged this satisfies the
   anti-tautology rule (independence from the target's primitive), but it does NOT give
   two structurally-different derivations of the kernel claim itself. If the orchestrator
   wants a stronger cross-engine guard, an alternative is to have one engine verify the
   projection identity symbolically (full `Integrate` in x, identity for symbolic
   `smallsigma`) while the other uses numerical quadrature — say so and I will revise.

4. **Ratio anchor 31.6785 is NOT a paper literal (MEDIUM, honesty flag).** See Reconcile
   item 3. It is the script's own computed ratio rounded; no `paper/` file contains
   `31.6785`. I kept it as a redundant cross-check but explicitly did NOT certify it as
   paper-sourced. Confirm the orchestrator is comfortable keeping a non-paper-sourced
   "anchor"; if not, it should be dropped rather than relabeled.

5. **Optional extra anchors (LOW).** `g'_* = 0.0714453558083195`, `Sigma_0(Pi_*) =
   1.80594111095636`, `T_m(Pi_*) = 0.901484054174205` are genuine appendix literals not
   currently asserted. Listed as optional strengthening; not required by the 5 findings.
