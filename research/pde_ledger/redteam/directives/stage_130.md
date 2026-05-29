---
unit_id: 130
batch: IV.4
created_at: 2026-05-29T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-29T21:35:03Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 130 (red-team v2 remediation)

> **RESOLVED (Claude+Codex consult `019e7594`, 2026-05-29 — record `redteam/codex_reviews/_consult_batch4.md`).**
> Codex CONCUR: the FKG/Chebyshev symmetrized-covariance certificate (F1) is a SOUND global
> proof of `dg/dΠ>0` for all Π>0, and the `2/π < g_- < 1` bracket is a sound uniqueness
> argument — both strictly better than the rejected 6-point sweep. Double-integral
> feasibility under the 600s cap is "reasonable but not guaranteed"; the fallback ladder in
> `## NEEDS_APPROACH_REVIEW` is endorsed, **especially the no-sampling HALT**. Codex (the
> fix-applier) RUNS F1 and walks the ladder (single-integral deviation form → guarded
> `Reduce` → `## Blocked`/HALT) if the double integral stalls; under NO fallback revert to a
> finite sweep. F1 is now LIVE — apply it.

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward. Respect the 10-minute (`timeout 600`) cap on every script run; a timeout is a FAILURE, not a reason to raise the cap.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

> ORCHESTRATOR NOTE — read `## NEEDS_APPROACH_REVIEW` at the bottom BEFORE applying F1. F1 prescribes a rigorous symbolic certificate (FKG/Chebyshev symmetrized covariance integral) whose engine-decidability under the 600s cap is not guaranteed from reading alone. `needs_user_resolution: true` is set deliberately. The orchestrator should settle the F1 approach with Codex (run the candidate, confirm it terminates and produces `0`/PASS) BEFORE this directive is treated as final. The fallback options and their risks are spelled out there.

---

## F1 — insufficient_verification (STRENGTHENED — supersedes the prior six-point sweep)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py:38-44` (REPLACE the existing `# Strict monotonicity sweep ...` block — the six-point loop inserted by the prior pass)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl:61-70` (REPLACE the existing `(* Strict monotonicity sweep ... *)` `Module[...]` block inserted by the prior pass)

**Issue:**
The notes box the strict inequality `dg_Π/dΠ > 0` for all `Π > 0`
(`notes/stages/moving_throat_pde_stage130_mouth_bias_map.md:79-82`), and this
global monotonicity is exactly what makes the Family-1 compensation point
`Π_*` **unique** — the paper card claims "monotonicity, and **unique** Family-1
point `Π_*`" (`paper/stages/stage_130.tex:16`,
`notes/...:105 "the unique explicit compensation point"`).

The prior pass added only a finite six-point sweep
`{1/10, 1/2, 1, 15088/10000, 3, 10}` (SymPy ~line 40 / Mathematica ~line 63).
Codex review R1 (`redteam/codex_reviews/stage_130.md:22-25`) correctly flagged
this: finite sampling is **not** a proof of a global strict inequality. A wrong
implementation with a negative derivative *between* the samples, or for
`Π < 0.1` or `Π > 10`, would pass undetected — and a non-global monotonicity
result does not establish uniqueness of `Π_*`. Sampling-as-proof is the precise
defect.

This finding replaces the sweep with the notes' OWN proof structure, made
machine-checkable: the FKG / Chebyshev correlation-inequality (symmetrized
covariance) sign certificate. It is rigorous (covers all `Π > 0`) and avoids a
transcendental `Reduce` over `e^Π`, which is hang-prone under the 600s cap (see
`## NEEDS_APPROACH_REVIEW`).

**Why not a direct numerator certificate.** Writing `dg/dΠ = M(Π)/D(Π)^2` with
`D = (4Π²+π²)(e^Π−1)` (so `D² > 0` for `Π>0`, clear), the auditor expanded the
numerator `M(Π) = N'D − ND'` as a quadratic in `E = e^Π`:
`M = 8π²Π·E² + c₁(Π)·E + 2π(4Π²−π²)`, where the constant (E⁰) term
`2π(4Π²−π²)` is **negative** for `Π < π/2` and `c₁(Π)` is largely negative. So
`M` is NOT a sum of manifestly-positive monomials in `(Π, E)`; its positivity
genuinely relies on the magnitude of `e^Π`, not merely `E>1`. Hence a clean
"all coefficients positive" certificate (the stage-082 idiom) does NOT exist
here, and approach (a) in its naive form is rejected. The covariance route below
is the rigorous substitute.

**The certificate (notes §2, made checkable).** The script already verifies the
covariance identity
`dg/dΠ + (E[fz] − g·E[z])/L = 0` (`cov_id = 0`), i.e.
`dg/dΠ = −(1/L)·Cov_Π(f, z)` with `f = cos(πz/2L)` and probability density
`p(z) = σ_Π(z)` on `[0,L]` (note `∫₀ᴸ p dz = 1`). The remaining gap is
`Cov_Π(f, z) < 0` for all `Π > 0`. Use the symmetrized (FKG/Chebyshev) identity,
valid for ANY normalized density:
```
Cov_Π(f, z) = (1/2) ∫₀ᴸ ∫₀ᴸ ( f(z₁) − f(z₂) ) ( z₁ − z₂ ) p(z₁) p(z₂) dz₁ dz₂.
```
The integrand factor `( f(z₁) − f(z₂) )( z₁ − z₂ )` is `≤ 0` pointwise on
`[0,L]²` because `f = cos(πz/2L)` is strictly decreasing on `[0,L]` (so the two
factors always have opposite sign), strictly `< 0` off the diagonal, and
`p > 0`. Therefore `Cov_Π(f,z) < 0` strictly, hence `dg/dΠ > 0`, for ALL
`Π > 0` — not at sampled points. The script must verify the two checkable
pieces:
1. **Identity** (algebraic, decidable): the closed-form double integral equals
   the already-computed `Cov := E[fz] − g·E[z]`. Equivalently it equals
   `−L · dg/dΠ`.
2. **Integrand sign** (decidable on `[0,L]²`): `(f(z₁)−f(z₂))(z₁−z₂) ≤ 0` for
   `z₁, z₂ ∈ [0,L]`, with the only zeros on `z₁ = z₂`. This is the
   strict-decrease property of `f`; encode it as: for `0 ≤ z₂ < z₁ ≤ L`,
   `f(z₁) − f(z₂) < 0`, which `Reduce`/`solve` over the bounded box CAN decide
   (cosine on a bounded interval, no `e^Π`).

Combining (1)+(2): the covariance is a manifestly-non-positive integrand
integrated against a positive density, strictly negative because the integrand
is `<0` on a positive-measure set — so `dg/dΠ > 0` globally, and uniqueness of
`Π_*` follows (F1c below).

**Required change (SymPy):**
REPLACE lines 38-44 (the prior six-point sweep block, from the comment
`# Strict monotonicity sweep:` through the `raise AssertionError(... at Pi={val}.")`)
with:

```python
# Global strict monotonicity dg/dPi > 0 for all Pi > 0 (notes section 2, boxed).
# Proof structure: dg/dPi = -(1/L) Cov_Pi(f, z) (already checked: cov_id == 0).
# We certify Cov_Pi(f, z) < 0 via the FKG/Chebyshev symmetrized identity,
# valid for any normalized density p = sigma, and the pointwise sign of its
# integrand. This is a GLOBAL certificate on Pi > 0, not a finite sweep.
dgPi = sp.diff(gPi, Pi)

# p is a probability density on [0, L]: integral over [0, L] is 1.
norm_p = sp.simplify(sp.integrate(sigma, (z, 0, L)))
if sp.simplify(norm_p - 1) != 0:
    raise AssertionError("sigma_Pi is not a normalized density on [0, L].")

# Covariance as already used in cov_id: Cov = E[f z] - g E[z].
Cov = sp.simplify(Efz - gPi * Ez)

# (1) Symmetrized double-integral identity for the covariance (FKG/Chebyshev).
#     Cov = (1/2) ∫∫ (f(z1)-f(z2))(z1-z2) p(z1) p(z2) dz1 dz2.
z1, z2 = sp.symbols("z1 z2", positive=True, real=True)
f1 = f.subs(z, z1)
f2 = f.subs(z, z2)
p1 = sigma.subs(z, z1)
p2 = sigma.subs(z, z2)
integrand_sym = sp.Rational(1, 2) * (f1 - f2) * (z1 - z2) * p1 * p2
Cov_double = sp.simplify(
    sp.integrate(sp.integrate(integrand_sym, (z1, 0, L)), (z2, 0, L))
)
if sp.simplify(Cov_double - Cov) != 0:
    raise AssertionError("Symmetrized covariance identity failed.")

# (2) Pointwise sign of the symmetrizer factor on [0, L]^2:
#     for 0 <= z2 < z1 <= L, f(z1) - f(z2) < 0 (f strictly decreasing).
#     Hence (f(z1)-f(z2))(z1-z2) <= 0 everywhere on [0,L]^2, < 0 off-diagonal.
#     This is a bounded-domain cosine inequality with NO exp(Pi); it is decidable.
sign_ok = sp.reduce_inequalities(
    [f1 - f2 < 0],
    [z1],
)  # noqa: F841  (we instead assert via a direct positivity argument below)
# Direct, engine-light certificate: g(z) := -f(z) = -cos(pi z/2L) is strictly
# increasing on [0,L] because g'(z) = (pi/(2L)) sin(pi z/2L) > 0 for 0 < z < L.
gprime = sp.diff(-f, z)
# gprime = (pi/(2L)) sin(pi z/(2L)); positive on (0, L). Certify by showing it
# equals a manifestly-positive closed form on the open interval.
if sp.simplify(gprime - (sp.pi / (2 * L)) * sp.sin(sp.pi * z / (2 * L))) != 0:
    raise AssertionError("f'(z) closed form mismatch.")
# sin(pi z/(2L)) > 0 strictly for z in (0, L) since its argument lies in (0, pi/2).
# Therefore f is strictly decreasing, the symmetrizer is <= 0 (< 0 off-diagonal),
# p > 0, so Cov < 0 and dg/dPi = -Cov/L > 0 for ALL Pi > 0.
print("Cov_Pi(f,z) (symmetrized) =", Cov_double)
print("f'(z) =", sp.simplify(gprime))
# Sanity: the certified sign must be consistent with the verified identity,
# i.e. dg/dPi = -Cov/L. (Non-tautological: a wrong gPi or wrong Cov breaks it.)
if sp.simplify(dgPi + Cov / L) != 0:
    raise AssertionError("dg/dPi = -(1/L) Cov consistency failed.")
print("Global strict monotonicity certified: dg/dPi = -(1/L) Cov_Pi(f,z) > 0 for Pi>0.")

# (F1c) Uniqueness of the Family-1 compensation point Pi_* on (0, oo):
# g is strictly increasing from g(0+) = 2/pi (~0.6366) to g(oo) = 1, so the
# equation g(Pi) = g_minus has at most one root, and exactly one iff
# 2/pi < g_minus < 1. g_minus is the lower-branch target (notes section 3).
g_lo = sp.N(2 / sp.pi, 30)
g_hi = sp.Integer(1)
if not (g_lo < g_minus < g_hi):
    raise AssertionError(
        "g_minus not strictly inside (2/pi, 1); uniqueness of Pi_* not guaranteed."
    )
print(f"Bracket for unique Pi_*: 2/pi = {g_lo} < g_minus = {g_minus} < 1 = {g_hi}")
```

> Codex implementation note (SymPy): if `sp.reduce_inequalities([...], [z1])`
> raises or is slow, DELETE those three lines (`sign_ok = ...` through the
> `# noqa` comment) — they are not load-bearing; the load-bearing sign argument
> is the `gprime` closed-form + the `sin` argument range, which is the
> non-sampling certificate. Do NOT replace it with a numeric sweep. If the
> double integral `Cov_double` does not return in closed form / does not
> simplify to `Cov`, STOP and append `## Blocked: F1` describing where it stalled
> — do NOT fall back to sampling.

**Required change (Mathematica):**
REPLACE lines 61-70 (the prior six-point sweep `(* Strict monotonicity sweep ... *)`
comment through the closing `];` of the `Module[...]`) with:

```mathematica
(* Global strict monotonicity dg/dpiM > 0 for all piM > 0 (notes section 2, boxed).
   Proof structure: dg/dpiM = -(1/lM) Cov_piM(f, z), already checked covId == 0.
   Certify Cov < 0 via the FKG/Chebyshev symmetrized identity and the pointwise
   sign of its integrand. GLOBAL certificate on piM > 0, not a finite sweep. *)
dgPi = D[gPi, piM];

normP = FullSimplify[Integrate[sigma, {z, 0, lM}], Assumptions -> $Assumptions];
expectZero["sigma_piM normalized on [0,lM]", normP - 1];

cov = FullSimplify[eFZ - gPi*eZ, Assumptions -> $Assumptions];

(* (1) Symmetrized double-integral identity for the covariance. *)
Clear[z1, z2];
asm2 = Element[{z1, z2}, Reals] && lM > 0 && piM > 0 && 0 <= z1 <= lM && 0 <= z2 <= lM;
f1 = (f /. z -> z1);
f2 = (f /. z -> z2);
p1 = (sigma /. z -> z1);
p2 = (sigma /. z -> z2);
integrandSym = (1/2)*(f1 - f2)*(z1 - z2)*p1*p2;
covDouble = FullSimplify[
  Integrate[Integrate[integrandSym, {z1, 0, lM}], {z2, 0, lM}],
  Assumptions -> $Assumptions];
expectZero["symmetrized covariance identity", covDouble - cov];

(* (2) Pointwise sign of the symmetrizer factor: f strictly decreasing on [0,lM].
   f'(z) = -(Pi/(2 lM)) Sin[Pi z/(2 lM)] < 0 for 0 < z < lM (argument in (0,Pi/2)).
   Certify the closed form of f'(z); its sign then follows from Sin > 0 on (0,Pi/2).
   This is a bounded-domain trig statement with NO Exp[piM]; decidable. *)
fPrime = D[f, z];
expectZero["f'(z) closed form", fPrime + (Pi/(2*lM))*Sin[Pi*z/(2*lM)]];
sinPos = Reduce[Sin[Pi*z/(2*lM)] > 0 && 0 < z < lM && lM > 0, z, Reals];
Print["Sin[Pi z/(2 lM)] > 0 on (0,lM) decided as: ", fmt[sinPos]];
If[TrueQ[sinPos =!= False],
  pass["f strictly decreasing on (0,lM) -> symmetrizer <= 0"],
  fail["f strictly decreasing on (0,lM)", sinPos]];

Print["Cov_piM(f,z) (symmetrized) = ", fmt[covDouble]];
(* Consistency: dg/dpiM = -(1/lM) Cov. Non-tautological: a wrong gPi or Cov breaks it. *)
expectZero["dg/dpiM = -(1/lM) Cov consistency", dgPi + cov/lM];
Print["Global strict monotonicity certified: dg/dpiM > 0 for all piM > 0."];

(* (F1c) Uniqueness of Pi_* on (0,oo): g strictly increasing from g(0+)=2/Pi to
   g(oo)=1, so g(piM)=gMinus has exactly one root iff 2/Pi < gMinus < 1. *)
gLo = N[2/Pi, 40];
If[TrueQ[gLo < gMinus < 1],
  pass["g_minus strictly inside (2/Pi, 1): Pi_* unique"],
  fail["g_minus inside (2/Pi, 1)", gMinus]];
Print["Bracket for unique Pi_*: 2/Pi = ", fmt[gLo], " < g_minus = ", fmt[gMinus], " < 1"];
```

> Codex implementation note (Mathematica): NO `*)` substring may appear inside
> any comment. `(expr /. v -> x)` substitutions are parenthesized as required.
> If `Reduce[Sin[...] > 0 ...]` returns `ConditionalExpression[...]`, strip it
> (the `=!= False` test already tolerates a non-`False` symbolic result; if it
> returns a bare `ConditionalExpression`, wrap with
> `sinPos = sinPos /. ConditionalExpression[e_, _] :> e;` before the `If`).
> If `covDouble` (the double `Integrate`) does NOT return in closed form or
> `FullSimplify` of `covDouble - cov` does not reach `0` within the 600s cap,
> STOP and append `## Blocked: F1` — do NOT fall back to a finite `Do` sweep.

**Claim manifest:**
- **M1** — the symmetrized-covariance identity + integrand-sign + normalization
  checks jointly exercise the notes §2 boxed result `dg_Π/dΠ > 0` for ALL
  `Π > 0`. Notes anchors: `notes/stages/moving_throat_pde_stage130_mouth_bias_map.md:64-82`
  (boxed `dg_Π/dΠ = −(1/L)Cov_Π(cos, z)` and boxed `dg_Π/dΠ > 0`, with the
  strict-decrease-of-cosine justification at `:76`).
- **M2** — the `g_minus ∈ (2/π, 1)` bracket exercises the notes §3 / paper-card
  **uniqueness** of `Π_*`. The strict monotonicity (M1) gives at-most-one root;
  the bracket gives exactly-one. Notes anchors:
  `notes/stages/moving_throat_pde_stage130_mouth_bias_map.md:84-89`
  (the monotone range `2/π → 1`), `:98` (`g_-^{F1} ≈ 0.758...`), `:105-110`
  ("the unique explicit compensation point" `Π_* ≈ 1.50882951349316`); paper
  card `paper/stages/stage_130.tex:16` ("monotonicity, and unique Family-1
  point").

Why this is non-tautological (no X−X): the certificate does NOT differentiate
`gPi` and then re-test the same derivative against itself. It (a) re-derives the
covariance from an INDEPENDENT symmetrized double integral and asserts it equals
`E[fz] − g E[z]`; a wrong `gPi`, wrong `σ`, or wrong `f` makes the identity
FAIL; (b) certifies `f`'s monotonicity from `f'(z)`'s closed form and the bounded
argument range — a wrong `f` (e.g. an increasing profile) flips the symmetrizer
sign and the consistency check `dg/dΠ + Cov/L = 0` would expose any internal
inconsistency; (c) the uniqueness bracket FAILS if `g_minus` ever drifts outside
`(2/π, 1)`. None of these can pass by construction.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 130` and
`redteam exec-mathematica 130` and confirms:
- SymPy: prints `Cov_Pi(f,z) (symmetrized) = ...`, `f'(z) = ...`, the consistency
  and bracket lines; NO `AssertionError`; exits 0. The old six lines
  `dg/dPi at Pi=...` are GONE (replaced).
- Mathematica: `PASS` lines for `sigma_piM normalized on [0,lM]`,
  `symmetrized covariance identity`, `f'(z) closed form`,
  `f strictly decreasing on (0,lM) -> symmetrizer <= 0`,
  `dg/dpiM = -(1/lM) Cov consistency`, and `g_minus strictly inside (2/Pi, 1)`;
  the six `dg/dpiM > 0 at piM=...` PASS lines are GONE (replaced); exits 0.
- Neither transcript contains a finite-sample monotonicity sweep.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl`
- summary: Replaced the finite monotonicity sweeps with the global FKG/Chebyshev covariance certificate and uniqueness bracket checks.
- deviation: Removed the optional SymPy `reduce_inequalities` probe after it raised on the bounded cosine inequality; the directive identified it as non-load-bearing and allowed deletion.

---

## F2 — insufficient_verification (PASSED in prior pass — RECONCILE, do NOT weaken)

**Status:** This finding was applied in the prior pass and Codex review R1 found
it SUBSTANTIVE (verdict PASS for both the SymPy and Mathematica boxed-form
equality checks; see `redteam/codex_reviews/stage_130.md:16-17,27-28`). It is
ALREADY present in both scripts and must be PRESERVED unchanged. Codex should
NOT re-edit it; simply record it as reconciled.

**Target (already present — verify, do not modify):**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py:16-18`
  (the `gPi_boxed = ...` definition and the `if sp.simplify(gPi - gPi_boxed) != 0: raise AssertionError(...)` guard)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl:43-44`
  (the `gPiBoxed = ...` definition and `expectZero["g_Pi matches paper boxed closed form", gPi - gPiBoxed];`)

**Issue (closed):** The integral-evaluated `g_Π` is asserted exactly equal to
the boxed closed form
`g_Π = 2Π(2Π e^Π + π) / ((4Π² + π²)(e^Π − 1))`. This can fail (a misquoted
constant breaks the simplify-to-zero) and exercises the paper's exact-`g_Π`
claim.

**Claim manifest:**
- **M3** — `gPi == gPi_boxed` / `expectZero["g_Pi matches paper boxed closed form", ...]`
  exercises the notes §1 boxed closed form. Notes anchor:
  `notes/stages/moving_throat_pde_stage130_mouth_bias_map.md:33-40`.

**Required action for Codex:** Do NOT modify these lines. Confirm they are
present and that both scripts still exit 0 after the F1 replacement. Append
`## Applied: F2` recording `deviation: none (pre-existing, reconciled)`.

**Verification command:**
- SymPy transcript: the printed `g_Pi = ...` line (line 1 of the saved
  transcript) is unchanged; no `AssertionError`.
- Mathematica transcript: `PASS: g_Pi matches paper boxed closed form` present.
- Both scripts exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl`
- summary: Confirmed the pre-existing boxed-form equality checks remain present and pass after the F1 replacement.
- deviation: none (pre-existing, reconciled)

---

## NEEDS_APPROACH_REVIEW (orchestrator + Codex must settle BEFORE applying F1)

The R1 finding asks for a GENUINE global strict-monotonicity proof of a
transcendental `g_Π` (it contains `e^Π`). The auditor evaluated the three
candidate routes:

- **(a) Sign-definite numerator.** REJECTED as a standalone certificate. The
  cleared numerator `M(Π) = N'D − ND'` is a quadratic in `E = e^Π`:
  `M = 8π²Π·E² + c₁(Π)·E + 2π(4Π²−π²)`. Its constant term `2π(4Π²−π²)` is
  NEGATIVE for `Π < π/2`, so `M` is not a sum of manifestly-positive monomials
  in `(Π, E>1)`; positivity genuinely depends on the magnitude of `e^Π`. No
  clean "all-positive-coefficients" certificate exists (unlike stage 082, whose
  numerator was rational).

- **(b) `Reduce`/`solve` shows the numerator has no positive root.**
  HIGH RISK. `Reduce[M > 0 && Π > 0, Π, Reals]` (or `solve(M, Π)` with `M`
  containing `e^Π`, `e^{2Π}`) is a transcendental decision problem. Mathematica's
  `Reduce` over mixed polynomial × exponential is frequently UNDECIDED or
  hang-prone and may blow the 600s cap (exit 124 = FAILURE). Not prescribed as
  primary for that reason.

- **(c) Covariance (FKG/Chebyshev) symmetrized-integral certificate.** CHOSEN as
  primary (F1 above). It is the notes' OWN proof (§2), is rigorous for ALL
  `Π > 0`, and replaces the transcendental `Reduce` over `e^Π` with: a double
  integral that the CAS evaluates in closed form, an identity check
  (`covDouble − cov == 0`), and a BOUNDED-domain cosine-monotonicity statement
  (`Sin[πz/(2L)] > 0` on `(0,L)` — no `e^Π`, reliably decidable). Uniqueness of
  `Π_*` (F1c) then follows from monotonicity + the `(2/π, 1)` bracket.

**The residual risk in (c):** the symmetrized DOUBLE integral
`∫₀ᴸ∫₀ᴸ (f₁−f₂)(z₁−z₂) p₁ p₂` is an elementary (exp × linear × cosine) integral
and SHOULD evaluate and simplify to `cov`, but the auditor did NOT run it (read
+ reason only). If `Integrate`/`FullSimplify` of `covDouble − cov` does not reach
`0` within 600s, the certificate stalls.

**Resolution requested:** the orchestrator + Codex should, BEFORE finalizing,
run the F1 SymPy and Mathematica blocks and confirm (i) `covDouble` returns in
closed form, (ii) `covDouble − cov` simplifies to `0`, (iii) the
`Sin[...] > 0` `Reduce` decides, and (iv) both scripts exit 0 under the cap. If
(c) stalls, the agreed fallbacks, in order, are:
  1. Replace the explicit double integral with the equivalent single-integral
     deviation form `cov = ∫₀ᴸ (f(z) − E[f])(z − E[z]) p(z) dz` plus the same
     `f`-monotonicity sign argument (lighter integral; same FKG logic).
  2. Only if both integral forms stall: attempt (b) `Reduce` with an explicit
     `timeout 600` and accept it ONLY if it returns `True`/no-positive-root
     cleanly; otherwise HALT and escalate to the user.
Do NOT under any fallback revert to a finite-sample sweep and call it a proof.
