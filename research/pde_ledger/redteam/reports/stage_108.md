---
unit_id: 108
batch: IV.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-28T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [notes/stages/moving_throat_pde_stage108_robustness_classes.md]
  paper_appendix: present
---

# Audit unit 108 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_108.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage108_robustness_classes.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage108_robustness_classes_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage108_robustness_classes_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage108_robustness_classes_mathematica_audit.txt`

## What the paper claims

Stage 108 classifies which explicit isotropic DtN outlet deformations preserve the canonical outgoing quadrupole normalization `chi_Q = 1` and which shift it. The boxed/quoted bottom line in the card is: "Pure scale is harmless; pure argument shift is harmless only if even moments remain canonical; additive throat-core data can move `chi_Q`." The notes enumerate Class A (pure scale `S`, gives `chi_Q=1`), Class B (scale+argument `beta`; even fingerprint forces `beta=1`, then `chi_Q=1`), Class C (additive core: `Sigma_2=-Sigma_0/9`, `Sigma_4=-Sigma_0/27`, `chi_Q = 3(S+9 Sigma_5)/(3S-Sigma_0)`), and an exact preservation submanifold `Sigma_5 = S(1-beta^5)/9 - Sigma_0/27`. The appendix (`stage_appendix_part04.tex` §"Outlet DtN and Robin outlet tests") additionally states the general formula `chi_Q = 3(S beta^5 + 9 Sigma_5)/(3S-Sigma_0)`. Critically, the stage card's `Checks` list ALSO requires: (2) "Check Robin and standalone mixed-pole limits before imposing compensation"; (3) "Check that the compensated branch preserves the even coefficients as well as the odd normalization." Those map to the appendix Robin result `chi_Q^R = 3/(3-rho_R)`, the standalone mixed-pole no-go (`kappa_W=-1/9` then `sigma_W=0`), and the compensated Robin–mixed branch (`chi_Q^hyb = (1-9 sigma_W gamma_W)/(1-sigma_W)`, preserved iff `gamma_W=1/9`).

## What the script claims to verify

Both scripts start from the literal outgoing DtN expansion `Lambda_out = -3 + z^2/3 + z^4/9 + i z^5/9` and verify: (A) pure-scale self-normalization leaves the normalized response unchanged; (B) the scale+argument coefficients are `m2=beta^2/9`, `m4=4 beta^4/81`, `chi=beta^5`, with even-match roots `beta∈{-1,1}` and `chi(beta=1)=1`; (C) the additive-core even match gives `Sigma_2=-Sigma_0/9`, `Sigma_4=-Sigma_0/27` and `chi_add = 3(S+9 Sigma_5)/(3S-Sigma_0)`, with preservation locus `Sigma_5=-Sigma_0/27`; and (D) a general beta-parameterized version reproducing `chi_gen = 3(S beta^5+9 Sigma_5)/(3S-Sigma_0)` and the notes submanifold `Sigma_5 = S(1-beta^5)/9 - Sigma_0/27`. The scripts do NOT contain any Robin, standalone mixed-pole, or compensated Robin–mixed checks.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Class A pure scale → `chi_Q=1` | sympy 30-32 / wl 35-37 (`Y_scale - Y_can`) | match (but tautological — see F3) |
| Class B scale+arg coeffs, `beta=1`, `chi_Q=1` | sympy 35-48 / wl 39-55 | match |
| Class C additive core `Sigma_2,Sigma_4`, `chi_add` formula | sympy 51-66 / wl 57-76 | match |
| Class C preservation locus `Sigma_5=-Sigma_0/27` | sympy 68-70 / wl 78-81 | partial (locus printed; assertion is a round-trip, see F4) |
| Appendix general `chi_Q=3(S beta^5+9 Sigma_5)/(3S-Sigma_0)` + submanifold | sympy 72-106 / wl 83-112 | match |
| Card Check #2: Robin limit `chi_Q^R=3/(3-rho_R)` | none | missing |
| Card Check #2: standalone mixed-pole no-go (`kappa_W=-1/9`, `sigma_W=0`) | none | missing |
| Card Check #3: compensated Robin–mixed even+odd preservation (`chi_Q^hyb`, `gamma_W=1/9`) | none | missing |

`paper_alignment: partial` — the deformation-class deliverables (the notes' content) are faithfully covered, but three explicit stage-card `Checks` items have no script-side counterpart.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 32 | `expect_zero(Y_scale - Y_can)` | Class A | no (S cancels by construction — tautology) |
| A2 | sympy | 46-47 | `set(roots)=={-1,1}` | Class B even-match roots | yes |
| A3 | sympy | 48 | `expect_zero(chi_arg(1)-1)` | Class B `chi_Q=1` | yes |
| A4 | sympy | 58-61 | unique solve for `Sigma_2,Sigma_4` | Class C even match | yes |
| A5 | sympy | 66 | `expect_zero(chi_add - 3(S+9 Sigma_5)/(3S-Sigma_0))` | Class C formula | yes |
| A6 | sympy | 70 | `expect_zero(chi_add.subs(Sigma_5,chi_pres)-1)` | Class C locus | no (solve→resubstitute round-trip) |
| A7 | sympy | 98-101 | `expect_zero(chi_pres_gen - (S(1-beta^5)/9 - Sigma_0/27))` | appendix/notes submanifold | yes |
| A8 | sympy | 102 | `expect_zero(chi_gen.subs(Sigma_5,chi_pres_gen)-1)` | submanifold | no (round-trip) |
| A9 | sympy | 103-106 | `expect_zero((chi_pres_gen+Sigma_0/27).subs(beta,1))` | D→C reduction | yes |
| B1 | mathematica | 37 | `expectZero[yScale - yCan]` | Class A | no (tautology, mirrors A1) |
| B2 | mathematica | 48-50 | `expectZero[m2Arg-..],[m4Arg-..],[chiArg-beta^5]` | Class B coeffs | yes |
| B3 | mathematica | 54-55 | roots `{-1,1}`, `chi(1)-1` | Class B | yes |
| B4 | mathematica | 71-76 | `Sigma_2+sigma0/9`, `Sigma_4+sigma0/27`, `chiAdd-...` | Class C | yes |
| B5 | mathematica | 80-81 | locus `+sigma0/27`; round-trip `(chiAdd/.sigma5->..)-1` | Class C locus | partial (80 anchored; 81 round-trip) |
| B6 | mathematica | 104-112 | submanifold, round-trip, D→C | appendix/notes submanifold | partial (104,109 anchored; 108 round-trip) |

## Findings

### F1 — paper_misalignment

**Severity:** high
**Subtype:** script_missing_paper_claim
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_108.tex:22-25`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage108_robustness_classes_sympy_audit.py` (entire file — no Robin/mixed/compensated checks)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl` (entire file)

**What's wrong:**
The stage-108 card lists three `Checks`:
> "Check pure scale and pure argument deformations separately. / Check Robin and standalone mixed-pole limits before imposing compensation. / Check that the compensated branch preserves the even coefficients as well as the odd normalization."

Only the first is exercised by the scripts. The Robin limit (`chi_Q^R = 3/(3-rho_R)`, appendix eq:app-part04-chi-robin), the standalone mixed-pole no-go (`kappa_W=-1/9` then `sigma_W=0`, eq:app-part04-mixed-pole-outlet), and the compensated Robin–mixed branch (`chi_Q^hyb = (1-9 sigma_W gamma_W)/(1-sigma_W)`, preserved iff `gamma_W=1/9`, eqs eq:app-part04-hybrid-chiQ / eq:app-part04-gammaW-canonical) are absent from both scripts.

However, the stage-108 NOTES file (`...stage108_robustness_classes.md`) only covers Classes A/B/C and the preservation submanifold — it does NOT mention Robin, mixed-pole, or compensated material. And the appendix block structure states (line 27) "Stages 107--113: low-frequency outlet deformations and the compensated Robin--mixed branch" and (line 86) "Stages~107--124 classify admissible low-frequency outlet/core deformations and identify the compensated Robin--mixed class." This strongly implies the Robin/mixed/compensated checks are owned by a SIBLING stage (≈109–113), and the stage-108 card's `Checks` list is a block-level checklist that over-scopes stage 108 specifically. The scripts and notes agree with each other (Classes only); the card's `Checks` items #2/#3 are the outlier.

Because the direction of resolution (re-scope the card vs. extend stage-108 scripts) changes the conceptual ownership of these checks across stages, this routes to the user, not Codex.

**Why this matters:**
A reader of the stage-108 card believes the unit's audit verifies the Robin/mixed/compensated no-go and survival results. If those are actually owned by a sibling stage, the card mis-advertises coverage; if they are genuinely stage-108's responsibility, the scripts have a real verification hole on load-bearing branch-selection results (`gamma_W=1/9` is the canonical-preservation condition for the compensated outlet).

**Required change:**
See `## Resolve before fix_loop` in the directive. Codex must not auto-resolve.

**Verification:**
After the user picks a direction: either the card `Checks` #2/#3 are re-scoped/annotated as block-level (paper-side edit, no script change), or new Robin/mixed/compensated checks appear in both scripts and exit 0.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:33-112`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage108_robustness_classes_sympy_audit.py:27-106`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`: same shared literal, same intermediate-quantity choreography, same solve targets, same truncation order, and even verbatim-copied inline comments. Three corresponding sections:

1. Class A. SymPy 30-32:
   `Y_scale = sp.series((-3*S)/(S*Lambda_out), z, 0, 6).removeO()` … `expect_zero('pure scale invariance', Y_scale - Y_can)`.
   Mathematica 35-37: `yScale = Expand[Normal[Series[(-3*sNorm)/(sNorm*lambdaOut), {z,0,5}]]]` … `expectZero["pure scale invariance", yScale - yCan]`. Same `-3*S`/`-3` numerators, same comparison.

2. Class C even-match. SymPy 56-58: `m2 = sp.simplify(-L2/L0); m4 = sp.simplify(L2**2/L0**2 - L4/L0); sol = sp.solve([sp.Eq(m2,1/9), sp.Eq(m4,4/81)], [Sigma2,Sigma4])`.
   Mathematica 63-65: `m2 = FullSimplify[-l2/l0,...]; m4 = FullSimplify[l2^2/l0^2 - l4/l0,...]; sol = Solve[{m2==1/9, m4==4/81}, {sigma2,sigma4}, Reals]`. Identical algebraic recipe.

3. Class D comment block. SymPy 72-73 comment: "# Class D: general scale + argument + additive (β-parameterized preservation submanifold). / # Notes box: Σ_5 = S(1 - β^5)/9 - Σ_0/27 (general locus); Class C is the β=1 reduction."
   Mathematica 83-84: the same two-line comment verbatim. The locus solve `chi_gen==1 for Sigma_5` and submanifold comparison are identical.

The Mathematica is NOT an independent re-derivation; it re-runs the same solve-and-assemble path. (It does add a few anchored sub-checks the SymPy only `print`s — `m2_arg-beta^2/9`, `Sigma2+sigma0/9`, etc. — but the derivation route is the same.)

**Why this matters:**
The second-engine policy exists so the two engines cross-check each other independently. A transliteration only confirms that the same algebra runs in two CAS dialects; it cannot catch an error baked into the shared choreography (e.g., a wrong truncation order or a wrong even-match target would be reproduced identically in both).

**Required change:**
Re-derive at least one class by a structurally different route in the `.wl` so the engines are not algorithmically identical. Minimal, in-scope option: for Class D, instead of porting the `Series`→`Coefficient`→`Solve` pipeline, build `chiGen` directly from the raw coefficients `L_0=-3 sNorm+sigma0`, `L_2=sNorm beta^2/3+sigma2`, `L_4=sNorm beta^4/9+sigma4`, `L_5=sNorm beta^5/9+sigma5`, impose the even-fingerprint conditions `-L_2/L_0 == 1/9` and `L_2^2/L_0^2 - L_4/L_0 == 4/81` as plain equations on those symbols (no `Series` of the rational function), assemble `chiGenAlt = 27*(-L_5/L_0)` under that solution, and assert `chiGenAlt - chiGen == 0`. See directive F2 for the concrete edit.

**Verification:**
The `.wl` Class D section no longer mirrors the `.py` `Series→Coefficient→Solve` order; a new independent-route quantity (`chiGenAlt`) is built from raw coefficients and compared, and the script exits 0 with the same `chi_gen` and submanifold output.

### F3 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage108_robustness_classes_sympy_audit.py:30-32`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:35-37`

**What's wrong:**
Class A computes `Y_scale = series((-3*S)/(S*Lambda_out))` and `Y_can = series((-3)/Lambda_out)`, then asserts `Y_scale - Y_can == 0`. The factor `S` cancels identically in `(-3*S)/(S*Lambda_out) = -3/Lambda_out` BEFORE any series is taken, so `Y_scale - Y_can` is structurally `0` for ANY `Lambda_out` and ANY nonzero `S`. The assertion has no failure mode and exercises no property of the canonical DtN expansion. The Mathematica mirrors this exactly (`yScale - yCan`).

**Why this matters:**
The check appears to verify Class A's claim (`chi_Q=1` under pure scale) but actually only verifies that a symbol cancels in a ratio. If `Lambda_out` were transcribed wrongly, this check would still pass.

**Required change:**
Anchor Class A to the paper's literal canonical fingerprint instead of to a self-construction. Replace the comparison target `Y_can` with the literal expansion from appendix eq:app-part04-Yout-dtn: `1 + z**2/9 + 4*z**4/81 + I*z**5/27` (sympy) and `1 + z^2/9 + (4 z^4)/81 + (I/27) z^5` (wl). This is falsifiable (it fails if `Lambda_out`'s coefficients are wrong) while still confirming the scale `S` drops out. See directive F3.

**Verification:**
The Class A assertion compares `Y_scale` against the literal canonical fingerprint and exits 0.

### F4 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage108_robustness_classes_sympy_audit.py:68-70,96,102`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:78,81,102,108`

**What's wrong:**
Two checks are solve-then-resubstitute round-trips that cannot fail:
- Class C (sympy 68-70 / wl 78,81): `chi_pres = solve(chi_add==1, Sigma5)` then `expect_zero(chi_add.subs(Sigma5, chi_pres) - 1)`. Substituting a value solved FROM `chi_add==1` back INTO `chi_add` and checking it equals 1 is guaranteed by `solve`.
- Class D (sympy 96,102 / wl 102,108): `chi_pres_gen = solve(chi_gen==1, Sigma5)` then `expect_zero(chi_gen.subs(Sigma5, chi_pres_gen) - 1)`. Same round-trip.

(The Class D submanifold assertion `chi_pres_gen - (S(1-beta^5)/9 - Sigma_0/27)` at sympy 98-101 / wl 104-107 IS a genuine, non-tautological anchor — keep it.)

**Why this matters:**
These two assertions look like independent confirmations of the preservation loci but verify nothing beyond `solve`'s own consistency. The Class C locus value `Sigma_5=-Sigma_0/27` is never anchored to the paper by an assertion (only `print`ed).

**Required change:**
- Class C: replace the round-trip `expect_zero(..., chi_add.subs(Sigma5, chi_pres) - 1)` with an anchor of the locus value to the paper: `expect_zero('Sigma5 locus (Class C) = -Sigma0/27', chi_pres - (-Sigma0/27))`. This matches the notes submanifold reduced to beta=1 and is falsifiable.
- Class D: the round-trip `expect_zero('general preservation locus check', chi_gen.subs(Sigma5, chi_pres_gen) - 1)` is redundant given the submanifold anchor already present (A7); demote it to a `print` (informational) rather than an assertion. See directive F4.
Mirror both in the `.wl`.

**Verification:**
Class C assertion now checks `chi_pres - (-Sigma0/27)` and exits 0; the Class D round-trip assertion no longer appears as a load-bearing `expect_zero`/`expectZero` (printed instead), and the submanifold anchor (A7) remains.

## Independent-derivation check (Mathematica)

Not independent — see F2. The `.wl` reproduces the `.py`'s exact `Series`→`Coefficient`→`Solve`→assemble pipeline with renamed symbols and verbatim-copied comments. It adds a handful of anchored sub-checks (`m2_arg - beta^2/9`, `Sigma2 + sigma0/9`, `Sigma4 + sigma0/27`, `Sigma5 locus + sigma0/27`) that the SymPy only prints, but the derivation route is identical, so it does not satisfy the independent-derivation requirement.

## Engine cross-check

Both engines agree at every shared quantity. SymPy output: `chi_add = 3*(S + 9*Sigma5)/(3*S - Sigma0)`, `chi_gen = 3*(S*beta**5 + 9*Sigma5)/(3*S - Sigma0)`, locus `-Sigma0/27`, general locus `S(1-beta^5)/9 - Sigma0/27`. Mathematica output: `chi_add = (-3*(9*sigma5 + sNorm))/(sigma0 - 3*sNorm)` (= same after sign-normalizing numerator/denominator), `chi_gen(beta) = (-3*(9*sigma5 + beta^5*sNorm))/(sigma0 - 3*sNorm)`, locus `-sigma0/27`, general locus `(-sigma0 - 3*(-1+beta^5)*sNorm)/27`. All final residual lines are `0` / `PASS`. No `engine_disagreement`.

Outputs are fresh: sympy output mtime (15:18) and mathematica output mtime (15:24) both postdate the script mtimes (15:09). No `stale_output`.

## Verdict justification

Verdict `findings` (4). The deformation-class core (Classes A/B/C and the general submanifold) is mathematically sound and well-anchored to the notes/appendix: I attacked the `chi_add`/`chi_gen` closed forms by re-deriving from `L_0,L_2,L_4,L_5` and they match the appendix formulas, and the `beta∈{-1,1}` root check and `chi(beta=1)=1` are genuine. What does NOT hold: (F1) the stage card's `Checks` #2/#3 demand Robin/mixed-pole/compensated coverage that neither script provides and that the notes do not mention — a paper-side scoping question routed to the user; (F2) the Mathematica is a transliteration, not an independent derivation, violating the second-engine policy; (F3) the Class A check is a pure symbol-cancellation tautology; (F4) two preservation-locus checks are solve-then-resubstitute round-trips. No `UNFIXABLE` (the math is consistent) and no `CRITICAL_DOWNSTREAM` (the verified class results are correct; F1 is a coverage/scoping question that must be user-resolved before any downstream impact is assessed). I confirm I read the stage card, the notes, and the appendix rows before opening the scripts. (Observation, not a script finding: the card's display title reads `\section[Stage~125]{Stage~125: ...}` while `\label{stage:108}` and the filename are 108 — a paper-prose label artifact outside script scope.)

## Self-test notes

I ran the required self-test on the prescribed F2/F3/F4 fixes. (1) Variable independence: no new `diff`/`D` introduced; the F2 `chiGenAlt` uses the raw `L`-coefficients (which genuinely depend on `sNorm,beta,sigma0,sigma2,sigma4,sigma5`) and the F3/F4 anchors depend only on already-present symbols. The Class A literal comparison reduces to 0 by hand-expansion of `-3/Lambda_out = 1 + z^2/9 + 4z^4/81 + i z^5/27`, and the Class C locus anchor reduces to `-Sigma_0/27 - (-Sigma_0/27) = 0`. (2) Parity: no unbounded integrals, N/A. (3) Trivial-case: each new `expect_zero` target reduces to literal 0 under the correct premises and becomes nonzero under a deliberately wrong `Lambda_out` coefficient or wrong even-match target, so all are falsifiable; `chiGenAlt` equals `chiGen` only because both encode the same physics by different routes, so the F2 assertion catches a route-specific error. (4)/(5) Paper round-trip: the F3 literal `1 + z^2/9 + 4z^4/81 + i z^5/27` is exactly eq:app-part04-Yout-dtn (appendix line 321), and the F4 Class C locus `-Sigma_0/27` is the beta=1 reduction of the notes submanifold box; F2 introduces no new constant; so none of the fixes create a new paper_misalignment.
