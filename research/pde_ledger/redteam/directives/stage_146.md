---
unit_id: 146
batch: IV.5
created_at: 2026-05-29T00:00:00Z
supersedes: 2026-05-27T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-29T21:41:23-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 146 (red-team v2 rewrite, paper-grounded)

## What this rewrite does

The 2026-05-27 directive PREDATES the Codex read-only review. That earlier
pass installed the integral-form affine checks correctly (F1/F2 of the old
directive) but the numeric *tolerances* were silently relaxed below the
required `1e-25`/`10^-25` fallback standard, and the Mathematica residuals
were `Chop`-ed to `10^-6` before being printed and asserted. The Codex
review (`redteam/codex_reviews/stage_146.md`) flagged exactly this: the
edited affine-law checks "can miss nontrivial residuals." This rewrite
encodes the review's three findings as concrete file:line edits.

The previously-installed direct-integral anchors for `g(Pi)` and `S_q(Pi)`
(SymPy lines 33-53, Mathematica lines 44-51) are GOOD independent checks and
the review rated them PASS — they are listed under "Reconcile
(tainted-applied)" below and are kept untouched.

Apply each finding in order. After applying, append an `## Applied: R<n>`
block under that finding with `files_changed`, `summary` (one sentence), and
`deviation` (or "none"). If a required change is ambiguous or unsafe to apply
mechanically, append `## Blocked: R<n>` with a question instead — skip that
finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond the
named edits. Do NOT touch paper.tex, notes/, or any prose documents.

After editing, RUN the affected scripts (`python3 <path>` for SymPy,
`math -script <path>` for Mathematica) under `timeout 600` and iterate until
they exit 0 with all in-file checks passing. A timeout (exit 124) is a
FAILURE — reformulate the math (integrate the closed-form integrand once,
then sample), never raise the cap.

---

## Root-cause analysis (read before applying R1/R2)

For both engines the affine-law residual algebraically collapses, by
linearity of the integral, to

```
g_eps residual  = (1 - eps) * ( <integral of (sigma|p=Pi_*) * cos(pi x/2)> - gminus )
               = (1 - eps) * ( gPi(Pi_*) - gminus )
```

i.e. it is controlled ENTIRELY by the root-finder residual
`gPi(Pi_*) - gminus`, NOT by machine epsilon. `Pi_*` is the transcendental
root of `gPi(Pi) - gminus = 0`, so:

- A **symbolic** zero is NOT achievable for the g_eps residual: `Pi_*` has no
  closed form, so `gPi(Pi_*) - gminus` cannot be `simplify`-ed to exact 0. It
  is genuinely a numeric quantity whose magnitude equals the root error.
- The `S_eps` residual collapses to `(1-eps)*( <integral of (sigma|p=Pi_*)*Kq>
  - Sformula(Pi_*) )`, where BOTH endpoints are the SAME quantity
  `Sformula(Pi_*)` (this is exactly what the F2 / direct-integral anchor
  proves). So the S_eps residual is pure quadrature noise — it reaches `~1e-31`
  in SymPy and literal `0`/`~1e-40` in Mathematica.

The decisive consequence for R1: in the current SymPy script the root is
`Pi_star = sp.N(sp.nsolve(gPi - gminus, 1.5), 30)` (line 73). **`sp.nsolve`
with no `prec` argument defaults to ~15 significant digits; the outer
`sp.N(..., 30)` only zero-pads — it does NOT improve root accuracy.** That is
why the saved g_eps residuals sit at `~3.6e-18` (transcript lines 156, 158) —
right at the 15-digit root floor, which the old `1e-15` tolerance was tuned to
just barely pass. To genuinely assert `< 1e-25` we must solve the root to
high precision so the residual drops well below `1e-25`, then evaluate the
residual at matching precision.

In Mathematica `pStar` is ALREADY solved at `WorkingPrecision -> 80,
AccuracyGoal -> 30, PrecisionGoal -> 30` (line 66), and the
`Pi_* compensation point diff` prints `0` to 39 digits (transcript line 25),
so `gFormula(pStar) - gMinus ~ 1e-40`. The Mathematica residuals are
genuinely below `10^-25` already; the only defect is that the script
`Chop`-s them at `10^-6` and asserts `< 10^-6`, hiding the true magnitude.
R2 is therefore a print/assert tightening, NOT a precision fix.

**Two-tier policy applied per residual (concrete below):**
- Tier 1 (symbolic exact zero): used ONLY where the engine can reduce the
  residual to literal `0`. For the affine-law residuals this is NOT possible
  (transcendental root), so they live in Tier 2.
- Tier 2 (high-precision numeric): print the RAW residual at >= 50 digits and
  assert `< 1e-25` / `< 10^-25`. No `Chop` before the assert. The precision
  budget must be high enough that the assertion has real headroom (root solved
  to >= 40 digits; residual evaluated at 50 digits).

---

## R1 — insufficient_verification (SymPy affine-law tolerance)

**Target file:**
`/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py`

**Issue (from review):** `if abs(float(g_res)) > 1e-15:` (current line 116)
relaxes the required `1e-25` fallback tolerance to `1e-15`. The g_eps residual
is `~3.6e-18` (transcript line 156), which passes `1e-15` but would FAIL
`1e-25` — so the check can miss a nontrivial first-order error. Per the
root-cause analysis the fix has two coupled parts: (a) solve `Pi_star` to high
precision so the residual genuinely drops below `1e-25`, and (b) evaluate the
residual at >= 50 digits and assert `< 1e-25` (raw residual printed, no
truncation).

### Edit 1a — raise the root-finder precision

**Current (line 73):**
```python
Pi_star = sp.N(sp.nsolve(gPi - gminus, 1.5), 30)
```

**Replace with:**
```python
Pi_star = sp.nsolve(gPi - gminus, 1.5, prec=50)
```

Justification: `nsolve(..., prec=50)` solves the root to ~50 significant
digits (the `prec` kwarg controls solver working precision; the bare call
defaults to ~15). `Pi_star` is then a 50-digit Float, so all downstream
`sp.N(..., 30)` consumers (lines 76-79) are unaffected in their printed
values, and the affine residual `(1-eps)*(gPi(Pi_*) - gminus)` drops from
`~1e-18` to `~1e-50`. This edit MUST land for Edit 1b to pass.

Anti-tautology: this strengthens, not weakens — `gPi`, `gminus`, and the
root condition are still independent constructions, and a bug in any of them
would now leave a residual `> 1e-25` that fails the assert. The check can
still fail.

### Edit 1b — tighten the residual print/assert to high-precision `< 1e-25`

**Current (lines 106-121):**
```python
# Symbolic simplify cannot reduce the residual (it contains both sqrt(3) and
# sqrt(4107 - 100*pi^2) factors that sp.simplify keeps). Fall back to numeric
# evaluation at two concrete eps samples, matching the Mathematica path.
g_eps_residual_expr = gbar_phys - (gminus + eps*(gbar_v - gminus))
S_eps_residual_expr = Sbar_phys - (Sformula.subs(Pi, Pi_star) + eps*(Sbar_v - Sformula.subs(Pi, Pi_star)))
for eps_val, label in [(sp.Rational(1, 10), "eps=1/10"), (sp.Rational(1, 2), "eps=1/2")]:
    g_res = sp.N(g_eps_residual_expr.subs(eps, eps_val), 30)
    S_res = sp.N(S_eps_residual_expr.subs(eps, eps_val), 30)
    print(f"g_eps affine law (integral form) at {label}: residual = {g_res}")
    print(f"S_eps affine law (integral form) at {label}: residual = {S_res}")
    if abs(float(g_res)) > 1e-15:
        raise AssertionError(f"g_eps affine law fails at {label}: {g_res}")
    if abs(float(S_res)) > 1e-15:
        raise AssertionError(f"S_eps affine law fails at {label}: {S_res}")
print("PASS: g_eps affine law (integral form) at eps=1/10 and eps=1/2")
print("PASS: S_eps affine law (integral form) at eps=1/10 and eps=1/2")
```

**Replace with:**
```python
# Two-tier check. Tier 1 (exact symbolic zero) is NOT reachable here: the
# residual collapses by linearity of the integral to (1-eps)*(gPi(Pi_*) -
# gminus), and Pi_* is a transcendental root with no closed form, so simplify
# cannot drive it to literal 0. Tier 2: evaluate the RAW residual at high
# precision (50 digits) and require it below the directive standard 1e-25.
# With Pi_star solved at prec=50 (line 73) the residual is ~1e-50 -- well
# inside the 1e-25 budget -- so the assertion has genuine headroom and is not
# tuned to the old 15-digit root floor.
g_eps_residual_expr = gbar_phys - (gminus + eps*(gbar_v - gminus))
S_eps_residual_expr = Sbar_phys - (Sformula.subs(Pi, Pi_star) + eps*(Sbar_v - Sformula.subs(Pi, Pi_star)))
TOL = sp.Float("1e-25", 50)
for eps_val, label in [(sp.Rational(1, 10), "eps=1/10"), (sp.Rational(1, 2), "eps=1/2")]:
    g_res = sp.N(g_eps_residual_expr.subs(eps, eps_val), 50)
    S_res = sp.N(S_eps_residual_expr.subs(eps, eps_val), 50)
    print(f"g_eps affine law (integral form) at {label}: residual = {g_res}")
    print(f"S_eps affine law (integral form) at {label}: residual = {S_res}")
    if sp.Abs(g_res) >= TOL:
        raise AssertionError(f"g_eps affine law fails at {label}: residual={g_res} >= 1e-25")
    if sp.Abs(S_res) >= TOL:
        raise AssertionError(f"S_eps affine law fails at {label}: residual={S_res} >= 1e-25")
print("PASS: g_eps affine law (integral form) at eps=1/10 and eps=1/2")
print("PASS: S_eps affine law (integral form) at eps=1/10 and eps=1/2")
```

Notes for Codex:
- Compare with `sp.Abs(...) >= TOL` on sympy `Float`s (not `abs(float(...))`)
  so the comparison itself is done at 50-digit precision, not collapsed to a
  ~15-16-digit Python `float` first.
- Keep the print labels `g_eps affine law (integral form) at <eps>: residual =`
  and `S_eps affine law (integral form) at <eps>: residual =`, and the two
  `PASS:` lines, BYTE-FOR-BYTE — the verify command greps them.
- If, after Edit 1a, the printed `g_res` is still NOT below `1e-25` (e.g. it
  stalls around `1e-18` again), do NOT widen the tolerance. That would mean
  `nsolve(prec=50)` did not take — stop and record a `## Blocked: R1` block
  noting the residual magnitude (see the consult block).

**Expected new PASS lines (SymPy transcript):**
```
g_eps affine law (integral form) at eps=1/10: residual = <something ~1e-50, < 1e-25>
S_eps affine law (integral form) at eps=1/10: residual = <something ~1e-31, < 1e-25>
g_eps affine law (integral form) at eps=1/2: residual = <something ~1e-50, < 1e-25>
S_eps affine law (integral form) at eps=1/2: residual = <something ~1e-31, < 1e-25>
PASS: g_eps affine law (integral form) at eps=1/10 and eps=1/2
PASS: S_eps affine law (integral form) at eps=1/10 and eps=1/2
```

**Verification command:** verifier runs `redteam exec-sympy 146` (exit 0). The
transcript must show all four residual lines with magnitudes `< 1e-25` and both
PASS lines.

## Applied: R1

- files_changed:
  - `scripts/moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py`
- summary: Raised the SymPy compensation root solve to 50-digit precision and tightened affine-law residual prints/assertions to raw 50-digit `< 1e-25` checks.
- deviation: none

---

## R2 — insufficient_verification (Mathematica affine-law Chop hides residual)

**Target file:**
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl`

**Issue (from review):** the affine-law block prints `Chop[..., 10^-6]` (hiding
the raw residual) and accepts anything `< 10^-6` (current lines 106-108,
116-118) — far weaker than `10^-25`, and could pass a nonzero first-order
error. Because `pStar` is already solved at `WorkingPrecision -> 80,
AccuracyGoal -> 30` (line 66) and the compensation residual is `0` to 39
digits, the residuals are genuinely tiny; this is a print/assert tightening,
not a precision fix. Tier-1 exact symbolic zero is NOT reachable (same
transcendental-root reason as SymPy), so use Tier 2: print the RAW residual at
50 digits and assert `< 10^-25`.

**Current (lines 95-121):**
```mathematica
(* Numeric-sample fallback: evaluate at two concrete eps values rather than
   simplifying an eps-polynomial residual (the integrate-with-numeric-pStar
   path produces complex near-zero coefficients that FullSimplify cannot
   reduce symbolically). *)
(* The numeric pStar substitution causes Integrate to produce complex near-zero
   residuals at low working-precision (~9-10 digits). Treat any value whose
   numerical magnitude is below 10^-6 (i.e., consistent with precision-9 zero)
   as satisfying the affine law. *)
gEpsRes = gBarPhys - (gMinus + eps*(gBarV - gMinus));
gEpsSample1 = N[gEpsRes /. eps -> 1/10, 40];
gEpsSample2 = N[gEpsRes /. eps -> 1/2, 40];
Print["g_eps affine law (integral form) at eps=1/10: ", fmt[Chop[gEpsSample1, 10^-6]]];
Print["g_eps affine law (integral form) at eps=1/2:  ", fmt[Chop[gEpsSample2, 10^-6]]];
If[NumericQ[gEpsSample1] && NumericQ[gEpsSample2] && Abs[gEpsSample1] < 10^-6 && Abs[gEpsSample2] < 10^-6,
  pass["g_eps affine law (integral form)"],
  fail["g_eps affine law (integral form)", {gEpsSample1, gEpsSample2}]
];

sEpsRes = sBarPhys - (sStar + eps*(sBarV - sStar));
sEpsSample1 = N[sEpsRes /. eps -> 1/10, 40];
sEpsSample2 = N[sEpsRes /. eps -> 1/2, 40];
Print["S_eps affine law (integral form) at eps=1/10: ", fmt[Chop[sEpsSample1, 10^-6]]];
Print["S_eps affine law (integral form) at eps=1/2:  ", fmt[Chop[sEpsSample2, 10^-6]]];
If[NumericQ[sEpsSample1] && NumericQ[sEpsSample2] && Abs[sEpsSample1] < 10^-6 && Abs[sEpsSample2] < 10^-6,
  pass["S_eps affine law (integral form)"],
  fail["S_eps affine law (integral form)", {sEpsSample1, sEpsSample2}]
];
```

**Replace with:**
```mathematica
(* Two-tier check. Tier 1 (exact symbolic zero) is NOT reachable: the residual
   collapses by linearity of the integral to (1-eps)*(gFormula[pStar]-gMinus),
   and pStar is a transcendental root with no closed form. Tier 2: evaluate the
   RAW residual at high precision (no Chop) and require Abs < 10^-25. pStar is
   solved at WorkingPrecision 80 / AccuracyGoal 30 (line 66), so the residual is
   ~10^-40 -- well inside the 10^-25 budget. *)
gEpsRes = gBarPhys - (gMinus + eps*(gBarV - gMinus));
gEpsSample1 = N[(gEpsRes /. eps -> 1/10), 50];
gEpsSample2 = N[(gEpsRes /. eps -> 1/2), 50];
Print["g_eps affine law (integral form) at eps=1/10: ", fmt[gEpsSample1]];
Print["g_eps affine law (integral form) at eps=1/2:  ", fmt[gEpsSample2]];
If[NumericQ[gEpsSample1] && NumericQ[gEpsSample2] && Abs[gEpsSample1] < 10^-25 && Abs[gEpsSample2] < 10^-25,
  pass["g_eps affine law (integral form)"],
  fail["g_eps affine law (integral form)", {gEpsSample1, gEpsSample2}]
];

sEpsRes = sBarPhys - (sStar + eps*(sBarV - sStar));
sEpsSample1 = N[(sEpsRes /. eps -> 1/10), 50];
sEpsSample2 = N[(sEpsRes /. eps -> 1/2), 50];
Print["S_eps affine law (integral form) at eps=1/10: ", fmt[sEpsSample1]];
Print["S_eps affine law (integral form) at eps=1/2:  ", fmt[sEpsSample2]];
If[NumericQ[sEpsSample1] && NumericQ[sEpsSample2] && Abs[sEpsSample1] < 10^-25 && Abs[sEpsSample2] < 10^-25,
  pass["S_eps affine law (integral form)"],
  fail["S_eps affine law (integral form)", {sEpsSample1, sEpsSample2}]
];
```

Notes for Codex:
- Removed BOTH `Chop[..., 10^-6]` wrappers and the two stale precision-9 /
  10^-6 justification comments. The RAW residual is now printed via `fmt[...]`.
- Parenthesize the replacement Rule before `N`: `N[(gEpsRes /. eps -> 1/10),
  50]` (the `/.` binds looser than function application; the parens make the
  argument unambiguous).
- `gBarPhys`/`sBarPhys`/`gBarV`/`sBarV` are computed by symbolic `Integrate`
  with the numeric `pStar` already substituted into `sigmaEps` (lines 90-94),
  so each `N[..., 50]` faithfully reflects `pStar`'s ~40-digit accuracy; no
  `Chop` is needed and none should be added.
- Do NOT change `fail[...]` to print a Chop-ed value — it must surface the raw
  residual if it ever fires.
- Keep the print labels `g_eps affine law (integral form) at eps=1/10: `,
  `... at eps=1/2:  ` (note: two trailing spaces, preserve them), and likewise
  for `S_eps`, plus the `pass[...]` argument strings BYTE-FOR-BYTE.
- If `Integrate[sigmaEps*Cos[Pi*x/2], {x,0,1}]` with numeric `pStar` is slow or
  returns unevaluated under the `timeout 600` cap: substitute `pStar` only into
  the `sigma` endpoint piece, integrate the closed-form pieces once
  (`Cos[Pi x/2]` against the exponential and against `kq` are the same
  integrands already proven in F2 / lines 44-51), and combine — do NOT raise
  the cap.

**Expected new PASS lines (Mathematica transcript):**
```
g_eps affine law (integral form) at eps=1/10: <raw residual, e.g. ...*10^-40, < 10^-25>
g_eps affine law (integral form) at eps=1/2:  <raw residual, < 10^-25>
PASS: g_eps affine law (integral form)
S_eps affine law (integral form) at eps=1/10: <raw residual, < 10^-25>
S_eps affine law (integral form) at eps=1/2:  <raw residual, < 10^-25>
PASS: S_eps affine law (integral form)
```
(The residuals will print as small numbers — possibly with a negligible
imaginary part from the numeric integral. If a tiny imaginary part appears,
compare `Abs[...]` as written — `Abs` of a complex near-zero is still the
magnitude, which the `< 10^-25` test handles correctly. Do NOT re-introduce
`Chop` to scrub it.)

**Verification command:** verifier runs `redteam exec-mathematica 146` (exit
0). The transcript must show all four raw residual lines with magnitude
`< 10^-25` and both PASS lines — and must NOT contain a literal Chop-to-`0`
for these four lines.

## Applied: R2

- files_changed:
  - `mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl`
- summary: Removed affine-law `Chop` masking and tightened Mathematica residual prints/assertions to raw 50-digit `< 10^-25` checks.
- deviation: Used the directive's symbolic endpoint-integral fallback for the affine residuals because numeric-`pStar` `Integrate` produced only low-precision zero residuals.

---

## R3 — paper_misalignment (Mathematica banner title carry-forward)

**Target file:**
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl:32`

**Issue (from review):** the banner literal says
`STAGE 146 — FINITE-CORRECTION EXPANSION FOR POSITIVE MOUTH-LAYER
DEFORMATIONS`, which does not match the paper-card title. This is a known
carry-forward (the prior pass corrected the stale stage number 129 -> 146 but
left "FINITE-CORRECTION" instead of "FIRST-ORDER"). This is a `.wl` SCRIPT
edit that rides this fix loop; it is flagged `paper_misalignment` only because
the script banner must match the paper-card title — the paper card needs NO
edit. `needs_user_resolution` stays false.

**Canonical title (CONFIRMED against the paper card):**
`paper/stages/stage_146.tex:1` reads
`\section[Stage 146]{Stage 146: First-Order Expansion for Positive
Mouth-Layer Deformations}`. So the canonical phrase is **"First-Order
Expansion for Positive Mouth-Layer Deformations"**, upper-cased for the banner
as **`FIRST-ORDER EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS`**.

**Current (line 32):**
```mathematica
banner["STAGE 146 — FINITE-CORRECTION EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS"];
```

**Replace with:**
```mathematica
banner["STAGE 146 — FIRST-ORDER EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS"];
```

Notes for Codex:
- Change ONLY `FINITE-CORRECTION` -> `FIRST-ORDER`. Keep the em-dash
  (`—`, U+2014) and `STAGE 146` exactly as-is.

**Verification command:** verifier confirms the new Mathematica transcript
contains the literal banner string `STAGE 146 — FIRST-ORDER EXPANSION FOR
POSITIVE MOUTH-LAYER DEFORMATIONS`.

## Applied: R3

- files_changed:
  - `mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl`
- summary: Updated the Mathematica banner title from `FINITE-CORRECTION` to `FIRST-ORDER`.
- deviation: none

---

## Reconcile (tainted-applied) — keep, do NOT re-edit

These checks were installed by the orchestrator-direct pass and rated PASS by
the Codex review as GOOD independent anchors. They are non-tautological
direct-integral-vs-closed-form comparisons that do NOT reuse the affine
primitive. Mark them tainted-applied (correct content, applied off-process)
and leave them UNCHANGED:

- SymPy `g(Pi) direct-formula` — `gPi_direct = sp.integrate(Sigma*cos(pi x/2),
  ...)` vs `gPi`, lines 33-42 (symbolic exact zero; transcript line 1: `= 0`).
- SymPy `S_q(Pi) direct-formula` — `Sq_direct = sp.integrate(Sigma*Kq, ...)`
  vs `Sformula`, lines 44-53 (symbolic, with documented 4-point numeric
  fallback at `Pi in {7/10,11/10,17/10,23/10}`, tol `1e-25`, transcript lines
  2-6).
- Mathematica `g(Pi) direct-formula` (line 48) and `S_q(Pi) direct-formula`
  (line 51) — `FullSimplify[Integrate[...]] - {g,s}Formula`, both symbolic
  exact zero (transcript lines 7-10).
- The pre-existing numeric kernel cross-check (`kernel check at Pi=...`, SymPy
  lines 62-68 / Mathematica lines 53-62, tol `1e-12`) and the
  `Pi_* compensation point` assertion (Mathematica line 79, tol `10^-20`) are
  independent and remain untouched.

These provide exactly the `gminus = int Sigma_*(x) cos(pi x/2)` /
`Sformula(Pi_*) = int Sigma_*(x) Kq(x)` anchoring that makes the R1/R2 affine
residual collapse to a pure root-error term — keeping them is what lets the
tightened `< 1e-25` assertions be meaningful rather than circular.

## Anti-tautology guard

Every check above can FAIL and exercises a distinct primitive. The affine-law
residuals are NOT compared against a re-statement of the affine form (the old
`expand(g_eps - (gminus + eps*(gbar-gminus)))` that was algebraically forced to
0); they integrate the physical convex profile `Sigma_eps = (1-eps)*Sigma_* +
eps*varsigma_test` against the physical kernels and compare to the affine RHS.
A bug in the kernels `cos(pi x/2)` / `Kq` or in the source `Sigma_*` perturbs
the residual above `1e-25` and trips the assert (those quantities do NOT feed
the root). Do NOT compare an integral against the same integral, and do NOT
reintroduce a tolerance loose enough (`1e-15` / `10^-6`) to absorb the very
error the check exists to catch.

**CONSULT Q3 caveat (batch 6, DISPUTE-resolved):** the affine residual collapses
to `(1-eps)*(gPi(Pi_*) - gminus)`, and `gminus` is ALSO used to compute
`Pi_star = nsolve(gPi - gminus, ...)`. So a WRONG `gminus` shifts `Pi_star` to
compensate and the residual stays small — this check is therefore **NOT** an
independent guard against `gminus`; it tests intercept-vs-direct-integral and the
kernels/source, given a correct `gminus`. `gminus` itself is the lower-branch
closed form `rF1 - sqrt(1+rF1^2)/2`, anchored at its owning branch stage; this
stage does not re-derive it. The verifier must NOT credit this stage with an
independent `gminus` check. (The fix here is exactly what `codex_reviews/stage_146.md`
R1/R2 asked for — tighten tolerance + root precision — which Codex confirmed
closes the finding.) No numeric literal here is
fabricated: `1e-25`/`10^-25` is the directive standard; `prec=50`/`N[...,50]`
is the evaluation precision; `varsigma_test = 6 x (1-x)` (normalized,
positive on (0,1)) is the existing test profile.

## For orchestrator/Codex consult

1. **The g_eps residual is NOT machine-zero — it equals the root-finder
   residual.** It collapses to `(1-eps)*(gPi(Pi_*) - gminus)`. The current
   SymPy run shows it at `~3.6e-18` only because `nsolve` ran at the 15-digit
   default. R1 Edit 1a (`nsolve(..., prec=50)`) should drop it to `~1e-50`.
   **Risk:** if `prec=50` does not actually tighten the SymPy root (e.g. a
   sympy version where `prec` interacts oddly with the `nsolve` start guess),
   the residual could stay near `1e-18` and the `< 1e-25` assert would FAIL.
   In that case Codex should NOT relax the tolerance — record `## Blocked: R1`
   with the printed residual and escalate. The honest reading would be that the
   "first-order affine law" holds only to root-finder precision, which the
   verifier must know. (Mathematica already solves the root at WP 80, so its
   side carries the high-precision confirmation independently.)

2. **The S_eps residual is genuine quadrature noise (~1e-31 SymPy / ~0
   Mathematica)** and comfortably clears `1e-25` — no risk expected there.

3. **Symbolic-zero is impossible for the affine residuals** (transcendental
   root `Pi_*`); Tier-1 cannot apply, so this stage legitimately rests on the
   Tier-2 high-precision numeric assertion. This is acceptable per the
   two-tier policy, but it means the stage's "exact" claim is "exact to the
   precision of the canonical-point root" — note this in the verification.

4. **High-precision-integrate cap risk (Mathematica):** `Integrate` of
   `sigmaEps*Cos[Pi x/2]` and `sigmaEps*kq` with a numeric 40-digit `pStar`
   substituted should be fast (the closed forms are the same ones F2 already
   evaluated symbolically at lines 44-51). If it nonetheless approaches
   `timeout 600`, reformulate by integrating the closed-form pieces once and
   substituting `pStar` after (do NOT raise the cap, do NOT switch to a
   low-WorkingPrecision `NIntegrate` that would re-introduce the very `~1e-9`
   noise the old `10^-6` Chop was hiding). SymPy `sp.N(..., 50)` of an already
   symbolic integral is cheap and carries no cap risk.
