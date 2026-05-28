---
unit_id: 148
batch: IV.5
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage148_representative_positive_families.md]
  paper_appendix: present
---

# Audit unit 148 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_148.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage148_representative_positive_families.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only `\input{stages/stage_148}` row at line 1330; no additional appendix prose)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage148_representative_positive_families_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage148_representative_positive_families_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage148_representative_positive_families_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage148_representative_positive_families_mathematica_audit.txt`

## What the paper claims

The stage card states (free-quote): "Uniform and self-matched derivative families bracket the first-order correction; interpolation reproduces the earlier compensation fraction." The card lists three verification checks: (i) first-order profile deformations are centred before covariance formulas are used; (ii) the rigidity kernel, not branch ambiguity, controls non-exponential corrections; (iii) the one-step nonlinear correction remains within the reduced mouth-layer regime.

The notes are more explicit. They enumerate four boxed deliverables: (D1) uniform-family shifts `δΠ_u/ε ≈ +1.699414961314297`, `δT_{m,u}/ε ≈ +0.508756302215084`; (D2) self-matched-derivative-family shifts `δΠ_d/ε ≈ -0.382993186095928`, `δT_{m,d}/ε ≈ -0.116943802151811`; (D3) convex-interpolation linear formulas for `δΠ_λ/ε` and `δT_{m,λ}/ε`, with neutral points `λ_{Π,0} ≈ 0.816081594488460` and `λ_{T,0} ≈ 0.813099276577333`; and (D4) the non-trivial consistency claim that `1 − λ_{Π,0}` agrees numerically with the **earlier Stage 228 broadening fraction** whose **closed form** is given in the notes as `ξ_* = (−37√3 − 5π² + 2√(4107 − 168π²))/(5(8 − π²)) ≈ 0.183918405511538`. Deliverable D4 is the load-bearing claim — it bridges the linear-response analysis of this stage to the exact positive-family compensation result of Stage 228.

The notes' H1 title says "Stage 250" while the filename and stage card both label this as Stage 148; this is paper-side housekeeping noise, not a script-side defect.

## What the script claims to verify

Both scripts numerically evaluate `(ḡ, S̄)` for the uniform and self-matched derivative families, plug them into the linear-response formulas built from `(gπ(Π), S(Π))` and their derivatives at the canonical lower-compensated branch point `Π_*` (defined by `gπ(Π_*) = g_-^{F1}` where `g_-^{F1} = rF1 − √(1+rF1²)/2` and `rF1 = √(12·(37/20)²/π² − 1)`), report the four shift numbers for D1 and D2, expand the linear formula in `λ` for D3, and solve for the neutral points `λ_{Π,0}` and `λ_{T,0}`. The final consistency line compares `1 − λ_{Π,0}` to a quantity named `xi_star`. In sympy this `xi_star` is a hardcoded float `sp.Float("0.183918405511538")`; in Mathematica it is the algebraic combination `xiStar = ((Pi/4) − gMinus)/((Pi/4) − 2/Pi)` built from the same `gMinus` used to define `Π_*`.

The Mathematica banner line reads `STAGE 131 — REPRESENTATIVE NON-EXPONENTIAL POSITIVE FAMILIES`, a stale label from the script's progenitor.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1 uniform shifts (notes §1) | sympy lines 49-55, math lines 49-57 | match (numerics agree with notes to ~28-30 digits) |
| D2 derivative-family shifts (notes §2) | sympy lines 61-67, math lines 59-67 | match |
| D3 convex-interpolation linear formulas + neutral points (notes §3) | sympy lines 70-83, math lines 69-80 | match |
| D4 numerical agreement of `1 − λ_{Π,0}` with **Stage 228 closed-form `ξ_*`** | sympy lines 86-89, math lines 83-85 | mismatch / insufficient — see F1 |
| Card check (i) "first-order deformations centred" | — | missing (no explicit centring check) |
| Card check (ii) "rigidity kernel controls" | implicit via use of `gPi`, `Sformula` | partial (no explicit assertion) |
| Card check (iii) "one-step correction within reduced regime" | — | missing |

Set `paper_alignment: partial` — numerical bracketing deliverables (D1-D3) match, the load-bearing cross-stage consistency check (D4) is structurally tautological in Mathematica and hardcoded in sympy, and three card-level checks have no script-side counterpart but the card is so terse that they may not be intended as algebraic assertions.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 80-81 | `sp.solve(Eq(dPi_lam,0), lam)[0]`, ditto for dT | D3 (neutral points) | yes (numerics fall out and match notes) |
| A2 | sympy | 89 | print of `(1-lam_Pi_zero) - xi_star`, no `assert` | D4 | no — informational print only, no `assert`; and `xi_star` is hardcoded float |
| A3 | math | 77 | `lamPiZero = lam /. First[Solve[gLam == gMinus, lam, Reals]]` | D3 | yes |
| A4 | math | 83-85 | `xiStar = FullSimplify[((Pi/4) - gMinus)/((Pi/4) - 2/Pi)]; expectZero["(1-lamPiZero) - xiStar", (1-lamPiZero) - xiStar]` | D4 | no — tautological; see F1 |
| A5 | sympy/math (all of D1-D3) | n/a | no `assert` / `expectZero`; only `print` of computed numbers | D1, D2, D3 | partial — numbers are computed but not asserted against the notes' boxed values |

## Findings

### F1 — insufficient_verification

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage148_representative_positive_families_mathematica_audit.wl:77-85`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage148_representative_positive_families_sympy_audit.py:80-89`

**What's wrong:**
The load-bearing consistency check (deliverable D4) is structurally hollow in both engines.

In Mathematica, `lamPiZero` is *defined* as `lam /. First[Solve[gLam == gMinus, lam, Reals]]` where `gLam = (1-lam)*(2/Pi) + lam*(Pi/4)`. Solving yields `lam = (gMinus - 2/Pi)/((Pi/4) - 2/Pi)`, so `1 - lamPiZero = ((Pi/4) - gMinus)/((Pi/4) - 2/Pi)`. The script then defines `xiStar = ((Pi/4) - gMinus)/((Pi/4) - 2/Pi)` directly and asserts `(1 - lamPiZero) - xiStar == 0`. This is an algebraic tautology — the two sides are constructed from the same `gMinus` via the same affine inversion. `expectZero` cannot fail no matter what `gMinus` is.

The notes (lines 124-133) make the actual claim: `1 - lamPiZero` must agree with **Stage 228's broadening fraction** whose **closed form** is
```
ξ_* = (−37√3 − 5π² + 2√(4107 − 168π²)) / (5(8 − π²)).
```
This closed form does not appear in either script. Neither script imports or re-derives the Stage 228 result; both replace it with an in-stage algebraic combination that is identical to `1 - lamPiZero` by construction.

In sympy, the situation is similar in substance but expressed differently: `xi_star = sp.Float("0.183918405511538")` is a hardcoded numeric stand-in (see F2) and the comparison is only a `print`, not an `assert` (the residual `1.69e-15` appears as a printed difference, not a tested one). So the sympy script does not actually test D4 at all.

**Why this matters:**
The whole point of D4 is the cross-stage agreement — that an algebra derived from the linear-response machinery of Stages 248-249 reproduces the closed-form positive-family compensation result of Stage 228. The current scripts demonstrate only that the linear-response algebra is internally self-consistent (Mathematica) or that it numerically equals a typed-in literal (sympy). They do not exercise the bridge to Stage 228 that the paper card and notes are claiming.

**Required change:**
Replace the `xi_star` construction in both engines with the explicit closed-form from the notes:
- sympy line 87: replace `xi_star = sp.N(sp.Float("0.183918405511538"), 30)` with the symbolic expression
  ```
  xi_star_closed = (-37*sp.sqrt(3) - 5*sp.pi**2 + 2*sp.sqrt(4107 - 168*sp.pi**2)) / (5*(8 - sp.pi**2))
  ```
  evaluated numerically (`sp.N(xi_star_closed, 30)`), and add an `assert` that `abs(sp.N((1 - lam_Pi_zero) - xi_star_closed, 30)) < sp.Float("1e-25")` (or equivalent symbolic `simplify` check).
- math lines 83-85: replace `xiStar = FullSimplify[((Pi/4) - gMinus)/((Pi/4) - 2/Pi)]` with
  ```
  xiStarClosed = (-37*Sqrt[3] - 5*Pi^2 + 2*Sqrt[4107 - 168*Pi^2]) / (5*(8 - Pi^2));
  ```
  and `expectZero["(1 - lamPiZero) - xi_star_closed", FullSimplify[(1 - lamPiZero) - xiStarClosed]]`. This forces the consistency check to compare two *independently constructed* expressions — one from the linear-response algebra (`lamPiZero`), one from the Stage 228 positive-family compensation closed form — rather than two expressions sharing `gMinus`.

**Verification:**
After the change, the Mathematica `expectZero` will still need to succeed (the two expressions are equal in closed form, since the notes give the same numerical value to ~15 digits). The fact that `FullSimplify` reduces a difference of two genuinely distinct symbolic forms to zero is the actual paper claim. In sympy, the new `assert` exercises the same bridge numerically.

### F2 — hardcoded_result

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage148_representative_positive_families_sympy_audit.py:87`

**What's wrong:**
`xi_star = sp.N(sp.Float("0.183918405511538"), 30)` is a literal numeric constant. The notes provide the closed form (line 132 of the notes file: `ξ_* = (-37√3 - 5π² + 2√(4107-168π²))/(5(8-π²))`). The script ignores the closed form and types in a 15-digit number directly. Even at face value the check on line 89 (`print((1-lam_Pi_zero) - xi_star)`) is a `print`, not a tested assertion, so a regression that breaks D4 would not fail the script.

**Why this matters:**
A hardcoded numeric constant masks the algebraic structure and can drift silently. Combined with F1, this means the sympy script has *zero* substantive assertions tying its computation to Stage 228.

**Required change:**
Subsumed by F1 — the fix is the same edit (replace the hardcoded float with the symbolic closed form from the notes and add an `assert`). After applying F1's fix, F2 resolves.

**Verification:**
Same as F1.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage148_representative_positive_families_mathematica_audit.wl` (whole file)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage148_representative_positive_families_sympy_audit.py` (whole file)

**What's wrong:**
The Mathematica script mirrors the SymPy script line-for-line. Compare:

- sympy lines 16-20 (`gPi`, `Sformula`):
  ```
  gPi = 2*Pi*(2*Pi*sp.exp(Pi)+sp.pi)/((4*Pi**2+sp.pi**2)*(sp.exp(Pi)-1))
  Sformula = sp.simplify(Pi*(kap*sp.tanh(kap) + Pi*(sp.exp(-Pi)*sp.sech(kap) - 1))
      / ((1-sp.exp(-Pi))*(kap**2-Pi**2)))
  ```
  vs math lines 32-33:
  ```
  gFormula = 2*p*(2*p*Exp[p] + Pi)/((4*p^2 + Pi^2)*(Exp[p] - 1));
  sFormula = p*(kap*Tanh[kap] + p*(Exp[-p]*Sech[kap] - 1))/((1 - Exp[-p])*(kap^2 - p^2));
  ```
  identical algebraic forms.

- sympy lines 32-41 (`AT`, `BT`):
  ```
  AT = sp.N(-(sp.Rational(9,1)/(40*T_star)) * (
      1/(gp_star*(1-S_star/4)) + Pi_star*Sp_star/(4*gp_star*(1-S_star/4)**2)
      ), 30)
  BT = sp.N((sp.Rational(9,1)/(40*T_star)) * Pi_star/(4*(1-S_star/4)**2), 30)
  ```
  vs math lines 43-47:
  ```
  aT = N[-(9/(40*tStar))*(1/(gPrimeStar*(1 - sStar/4)) + pStar*sPrimeStar/(4*gPrimeStar*(1 - sStar/4)^2)), 30];
  bT = N[(9/(40*tStar))*pStar/(4*(1 - sStar/4)^2), 30];
  ```
  identical line-by-line algebra.

- Both scripts proceed in identical order: uniform block → derivative block → convex `gLam`, `sLam`, `dPiLam`, `dTLam` → solve for neutral points → final consistency.

There is no independent derivation: Mathematica's banner even reads `STAGE 131 — REPRESENTATIVE NON-EXPONENTIAL POSITIVE FAMILIES`, suggesting the .wl was templated from another stage's .wl and the assertions ported wholesale.

**Why this matters:**
The two-engine policy requires that the engines reach the same result by independent paths. As written, an algebraic error in the `gPi` / `sFormula` / `aT` / `bT` construction would propagate identically to both engines and would not be caught by engine cross-checking.

**Required change:**
Within the unit's existing scope (no new features), at minimum the Mathematica script should restructure the algebra so that the canonical-branch coefficients are derived from `(gFormula, sFormula)` by a different algebraic route — for example, by introducing the canonical retuning formula as `ΔΠ = -Δg / g'(Π_*)` and `ΔT = (9/(40·T_*)) · (Δσ - (1/4)·Π_*·ΔS/(1-S_*/4)) · ... ` worked through symbolically rather than copying the sympy's pre-baked `aT`, `bT` form. The simplest concrete change is to write the linear-response in the form
```
dPiLam = -(gLam - gStar)/gPrimeStar;
dTLam = (9/(40*tStar)) * ( (dPiLam) - (pStar/(4*(1 - sStar/4))) * (sLam - sStar)/(1 - sStar/4) ... )
```
i.e., expressing `dTLam` directly in terms of the upstream `dSigma = dPi/(1-S/4) + Pi · dS/(1-S/4)^2` rather than in terms of pre-collected `aT`/`bT`. The two engines then arrive at numerically agreeing `dT` values via algebraically distinct intermediate forms.

Also fix the banner string on line 26 to read `STAGE 148 — REPRESENTATIVE NON-EXPONENTIAL POSITIVE FAMILIES`.

**Verification:**
Auditor re-checks that the Mathematica script's path from `(gLam, sLam)` to `dTLam` no longer matches the sympy script's path symbol-for-symbol; numerical `dT_u`, `dT_d`, `dT_λ` outputs must still agree with sympy to within the printed precision (28-30 digits).

## Independent-derivation check (Mathematica)

The .wl is a transliteration — same coefficient algebra (`aT`, `bT`), same intermediate steps (`p*`, `gStar`, `sStar`, `gPrimeStar`, `sPrimeStar`, `sigmaStar`, `tStar`), same call structure. See F3 above for quoted side-by-side excerpts. Stale banner ("STAGE 131") corroborates the transliteration origin.

## Engine cross-check

Both engines produce numerics that agree to the precision they claim. Sample (uniform branch, `dPi/eps`):
- sympy: `1.69941496131429664915468699198`
- math:  `1.6994149613142967238719481043693877705099...`

`(1 - λ_{Π,0}) − ξ_*`:
- sympy: `1.6930901125533637241e-15` (printed, not asserted — this is a numeric-precision artifact since `xi_star` is a 15-digit float)
- math: `0` exactly (because both sides are symbolic and reduce to the same form via the construction in F1)

Engines agree at the level they claim — no `engine_disagreement` finding.

## Verdict justification

The four numerical deliverables D1-D3 (uniform, derivative, and interpolated linear-response shifts) compute and match the notes' boxed values, and the engines agree to the printed precision. The load-bearing claim D4 — that `1 − λ_{Π,0}` reproduces the Stage 228 closed-form broadening fraction `ξ_*` — is structurally hollow: in Mathematica the `expectZero` compares two expressions both constructed from `gMinus` via the same affine inversion (tautological by construction), and in sympy the comparison is a `print` against a hardcoded 15-digit literal with no `assert`. The Stage 228 closed-form `(−37√3 − 5π² + 2√(4107 − 168π²))/(5(8 − π²))` is given explicitly in the notes but appears nowhere in either script. Combined with the line-by-line transliteration between engines, the audit concludes `findings` — the load-bearing cross-stage bridge is not exercised, and the dual-engine independence policy is violated. Verdict is not `stop_cold`: the fix is local (replace the in-stage `ξ_*` construction with the closed form from the notes and add a real assertion; restructure the Mathematica algebra). Downstream Stages 154-163 consume the *bias-neutral interpolation* result, not the `ξ_*` consistency comparison itself, so a correction to F1 will tighten the verification without changing any numerical carry-forward.

## Self-test notes

- Verified the tautology in F1 by hand-tracing: `gLam = (1-λ)·(2/π) + λ·(π/4)`; setting `gLam = gMinus` gives `λ = (gMinus - 2/π)/((π/4) - 2/π)`, so `1 - λ = ((π/4) - gMinus)/((π/4) - 2/π)`, which is exactly `xiStar` as defined on line 83. So `(1 - lamPiZero) - xiStar = 0` identically — `FullSimplify` cannot return anything else.
- Variable-independence trap: no `sp.diff` over a constant; `sp.diff(gPi, Pi)` and `sp.diff(Sformula, Pi)` both depend on `Pi`, so non-trivial derivatives.
- Symmetry/parity: no integrals over unbounded domains here; all integrals are over `[0,1]` and evaluated symbolically by the notes; the scripts use the closed-form results of those integrals, not the integrals themselves.
- Trivial-case pre-check for F1's required change: substituting `gMinus → 2/π` makes `1 - lamPiZero → 0` and `xi_star_closed → ((π/4) - 2/π)/((π/4) - 2/π) = 1`. Wait — `xi_star_closed` is independent of `gMinus`, so substituting only affects the LHS. The check would then assert `0 - ((-37√3 -5π² + 2√(4107-168π²))/(5(8-π²))) == 0` which is **not** identically true; it equates to whether `(−37√3 − 5π² + 2√(4107 − 168π²))/(5(8 − π²)) = 0`. Since that closed-form numerically equals `~0.184`, not zero, the new assertion would now **fail** if `gMinus` were perturbed away from its true value. That confirms the new check is non-tautological and exercises the bridge. Original check passes regardless of `gMinus` — confirming F1.
- Path-specification: directive targets `.py` in `scripts/` and `.wl` in `mathematica/` — correct.
- Paper round-trip: the closed form `(-37√3 - 5π² + 2√(4107-168π²))/(5(8-π²))` is taken directly from notes line 132; no new constants introduced; no new paper_misalignment risk.
