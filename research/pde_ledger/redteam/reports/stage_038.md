---
unit_id: 038
batch: III.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md
  paper_appendix: present
---

# Audit unit 038 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_038.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 54, summary-chain line 309, checklist row 345)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage038_dimensionless_continuum_placement_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage038_dimensionless_continuum_placement_mathematica_audit.txt`

## What the paper claims

Stage 038 compresses the Stage-037 continuum data into six dimensionless ratios `(eps_eta, eps_W, rho, Z_W, delta_0, Lambda)` and proves that the selected quadrupole placement splits into three structural lanes (geometry, product, redistribution). The card's `\stagefield{Output}` is verbatim: "The placement map (eq:app-stage038-placement-map)--(eq:app-stage038-Rtarget) and product law (eq:app-stage038-product-law)." The load-bearing identities boxed in the card are:

- `delta = delta_0 / (1 - eps_eta)`
- `M_mix = 8 Z_W (1+rho)^2 / [pi^2 (1-eps_eta)(1-eps_W)]`
- `R_target = Lambda (1-eps_eta)(1-eps_W)^2 / [Z_W (1+rho)^2]`
- `R_target M_mix = 8 Lambda (1-eps_W)/pi^2 = 54 G c_s^5 K_W^eff (1-eps_W) / (5 a^5 c^5 mu_W)`

The card's `\stagefield{Checks}` additionally requires that the product law be invariant under redistributions of `Z_W`, `rho`, and `eps_eta` that leave the product lane fixed. The notes file extends the .tex with (i) a corollary form for the outgoing transfer factor `beta_0 = (mu_W/mu_eta)(Keta_eff/KWeff) Z_W (1+rho)^2/(1-eps_W)^2` (Section 2) and (ii) nine "one-way parameter tendencies" giving the signs of `d{delta, M_mix, R_target}/d{eps_eta, eps_W, Z_W, rho}` on the natural branch `1+rho > 0` (Section 4). The part-3 appendix row 54 classifies Stage 038 as ExactClosure with deliverables "Ratios (eps_eta, eps_W, rho, Z_W, delta_0, Lambda), placement map, and product law"; the summary chain at line 309 names "Stage 038 placement product" as the link out of Stage 037.

## What the script claims to verify

The SymPy script declares the Stage-037 closed forms (`A`, `delta`, `M_mix`, `beta_0`, `R_target`) symbolically in microscopic variables, then defines the six dimensionless ratios via a substitution dictionary, and asserts via `expect_zero` that each microscopic expression collapses to the paper's compact form. Specifically it asserts the four placement identities (`delta`, `M_mix`, `R_target`, `beta_0`), the product law in both compact and `NQ K_W^eff(1-eps_W)/mu_W` forms, the algebraic factorization of each of the nine partial derivatives `d{M, R, delta}/d{eps_eta, eps_W, Z_W, rho}`, and the sign of each derivative (residual is multiplied by a manifestly-positive template before reducing to `±1`). The Mathematica script mirrors the same assertion set with `expectZero` and a `FullSimplify[Together[Expand[...]]]` simplifier, under the same positivity assumptions. Both scripts emit transcripts that show the simplified bottom line of every assertion as `0`.

## Paper / script cross-check

| Paper-side deliverable | Script-side coverage | Status |
|---|---|---|
| Ratios `(eps_eta, eps_W, rho, Z_W, delta_0, Lambda)` defined | substitution dictionary `subs_dimless` at sympy:84-91 / `applyDimless` at wl:77-95 | match |
| `delta = delta_0/(1-eps_eta)` (eq:app-stage038-placement-map, first) | sympy:109, wl:102 | match |
| `M_mix = 8 Z_W (1+rho)^2 / [pi^2 (1-eps_eta)(1-eps_W)]` (eq:app-stage038-placement-map, second) | sympy:110-113, wl:103-106 | match |
| `R_target = Lambda (1-eps_eta)(1-eps_W)^2 / [Z_W (1+rho)^2]` (eq:app-stage038-Rtarget) | sympy:114-117, wl:107-110 | match |
| Compact product law `R_target M_mix = 8 Lambda (1-eps_W)/pi^2` (eq:app-stage038-product-law, first equality) | sympy:134-135, wl:123-124 | match |
| Microscopic form `= 54 G c_s^5 K_W^eff (1-eps_W)/(5 a^5 c^5 mu_W)` (eq:app-stage038-product-law, second equality) | sympy:138-143, wl:125-133 | match |
| Invariance of product under `(Z_W, rho, eps_eta)` redistribution (Checks field) | implicit: the asserted product `8 Lambda (1-eps_W)/pi^2` contains none of `Z_W, rho, eps_eta`; the equality alone witnesses the invariance | match (implicit) |
| Notes Section 2: `beta_0 = (mu_W/mu_eta)(Keta_eff/KWeff) Z_W (1+rho)^2/(1-eps_W)^2` | sympy:118-121, wl:111-114 | extra (covered by notes) |
| Notes Section 4: nine derivative factorizations + nine signs | sympy:157-217, wl:143-181 | extra (covered by notes) |

`paper_alignment: aligned`. Every paper-card deliverable has a faithful script-side check, and every script-side check beyond the .tex (`beta_0` corollary and derivative-sign tendencies) is explicitly authorized by the notes file that informed the card. There is no `paper_missing_script_claim` and no `script_missing_paper_claim`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 109 | `expect_zero(delta_dimless - delta0/(1-eps_eta))` | placement map (delta) | yes |
| A2 | sympy | 110-113 | `expect_zero(M_dimless - 8 Z_W (1+rho)^2/(pi^2(1-eps_eta)(1-eps_W)))` | placement map (M_mix) | yes |
| A3 | sympy | 114-117 | `expect_zero(R_dimless - Lambda(1-eps_eta)(1-eps_W)^2/(Z_W(1+rho)^2))` | R_target identity | yes |
| A4 | sympy | 118-121 | `expect_zero(beta_dimless - (mu_W/mu_eta)(Keta_eff/KWeff) Z_W(1+rho)^2/(1-eps_W)^2)` | notes Section 2 corollary | yes |
| A5 | sympy | 134-135 | `expect_zero(R_dimless*M_dimless - 8 Lambda(1-eps_W)/pi^2)` | product law (compact) | yes |
| A6 | sympy | 141-143 | equivalence through `NQ K_W^eff(1-eps_W)/mu_W` | product law (microscopic) | yes |
| A7-A15 | sympy | 167-175 | nine derivative-factorization residuals reduced to `0` | notes Section 4 tendencies (factorizations) | yes |
| A16-A24 | sympy | 182-217 | nine sign-residual reductions to literal `±1 -> 0` | notes Section 4 tendencies (signs) | yes |
| B1-B24 | mathematica | 102-181 | structurally parallel assertion list under `FullSimplify` | same paper/notes claims as A1-A24 | yes |

Every assertion traces to either the paper card or the notes file. The substitution dictionary itself does the algebraic lifting; the assertions then test that the residual is zero — non-tautological because the residual is built from microscopic objects (`K_U Keta_eff - c_etaU^2`, `(K_U c_etaW + c_UW c_etaU)^2`, `K_U K_W^eff - c_UW^2 sigma`) that must rearrange across multiple identities to vanish. If the Stage-037 closed forms were wrong, or if the six ratio definitions were misnumbered, the residuals would not simplify to zero.

I verified the algebra by hand:

- `delta = pi^2 T_w K_U/(L^2 K_U Keta_eff (1-eps_eta))` with `T_w = delta_0 L^2 Keta_eff/pi^2` gives `delta_0/(1-eps_eta)`. Matches A1.
- `M_mix` numerator `(K_U c_etaW + c_UW c_etaU)^2` is left factored in SymPy's intermediate, so `c_UW*c_etaU -> rho K_U c_etaW` substitutes inside the square, yielding `K_U^2 c_etaW^2 (1+rho)^2 -> K_U^2 Z_W Keta_eff K_W^eff (1+rho)^2`. Denominator splits as `pi^2 K_U Keta_eff (1-eps_eta) * K_U K_W^eff (1-eps_W) / 1` (sigma cancels the `c_UW^2 sigma` term cleanly). Ratio = `8 Z_W (1+rho)^2 / [pi^2 (1-eps_eta)(1-eps_W)]`. Matches A2.
- `R_target = (54 G c_s^5/(5 a^5 c^5)) A (pi^2 / (8 beta_0))`. With `A = Keta_eff (1-eps_eta)/mu_eta`, `beta_0 = (mu_W/mu_eta) Z_W Keta_eff (1+rho)^2/(K_W^eff (1-eps_W)^2)`, and `G = 20 Lambda a^5 c^5 mu_W/(27 pi^2 c_s^5 K_W^eff)`, the prefactors collapse to `(54*20)/(8*5*27) = 1080/1080 = 1`. Matches A3.
- Product `R_target M_mix = [Lambda (1-eps_eta)(1-eps_W)^2/(Z_W(1+rho)^2)] * [8 Z_W (1+rho)^2/(pi^2 (1-eps_eta)(1-eps_W))] = 8 Lambda (1-eps_W)/pi^2`. The three "redistribution" variables `eps_eta, Z_W, rho` cancel — this is the paper's `Checks` invariance manifest in the symbolic form. Matches A5.
- The microscopic form `54 G c_s^5 K_W^eff (1-eps_W)/(5 a^5 c^5 mu_W)` matches `8 Lambda (1-eps_W)/pi^2` after `G -> 20 Lambda a^5 c^5 mu_W/(27 pi^2 c_s^5 K_W^eff)`: numerically `54/5 * 20/(27 pi^2) = 1080/(135 pi^2) = 8/pi^2`. Matches A6.

## Independent-derivation check (Mathematica)

The Mathematica script is **structurally parallel** to the SymPy script: same banner / subbanner / pass / fail scaffolding, same five-section layout (`1. Stage-20 continuum formulas`, `2. Dimensionless kernel substitutions`, `3. Exact product relation`, `4. Exact derivative factors`), same nine derivative names, same nine sign checks. The substitution dictionary is the same set of six identities (only renamed `eps_eta -> epsEta`, etc., to follow Mathematica naming conventions).

However, the two engines do not transliterate each other's algebra. They both **derive** the dimensionless map by applying the **same definitions** — but the substitution dictionary IS the definition of the six ratios. It is paper-side material, not engine-internal scaffolding. The simplification machinery is engine-specific:

- SymPy `apply_dimless` (sympy lines 93-102) uses `sp.expand(... .subs(...))` in two passes, first substituting the cross-product `c_UW * c_etaU -> rho K_U c_etaW`, then the squares.
- Mathematica `applyDimless` (wl lines 77-95) uses `Expand[PowerExpand[Factor[... /. ...]]]`, a different normal-form pathway, and additionally introduces explicit `cUW^4` and `cEtaW^(-2)` rewrite rules to handle higher powers that `Factor` can introduce.

The intermediate sentinel values printed in section 1 differ between engines (sympy: `M_mix = 288*L^2 * ...`; mathematica: `M_mix = -288*ell^2 * ...` due to ordering of `(cEtaU^2 - kU(kEta+6 tOmega))`), confirming the engines reach the answer through independent canonicalization passes rather than mirroring each other. Both arrive at zero residuals.

I checked specifically for the answer-baked rule that would betray a transliteration — e.g., a substitution like `(cEtaU*cUW + cEtaW*kU)^2 -> zW kEtaEff kWEff kU^2 (1 + rho)^2` — and **did not find one**. The Mathematica script substitutes only atomic monomials (`cUW*cEtaU`, `cUW^4`, `cEtaW^2`, `cEtaW^(-2)`, `cEtaU^2`, `cUW^2`, `tw`, `gNewton`); it does not bake in the factored numerator answer. Verdict: not a `mathematica_transliteration` finding. The parallel structure is in the paper-side definitions, which both engines must respect by construction.

## Engine cross-check

Both engines emit identical bottom-line assertions:

- SymPy output (lines 18-21): four placement assertions all `= 0`.
- Mathematica output (lines 18-25): same four, all `PASS:`.
- SymPy output (lines 30-32): product law `= 0`, equivalent form `= 0`.
- Mathematica output (lines 34-37): same two, `PASS:`.
- SymPy output (lines 37-54): all nine derivative + nine sign assertions `= 0`.
- Mathematica output (lines 43-78): same eighteen, `PASS:`.

The "Stage-20 formulas" pretty-printed in section 1 differ in cosmetic sign / ordering (e.g., SymPy presents `K_U*(K_eta+6*T_Omega) - c_etaU**2` while Mathematica presents `-cEtaU^2 + kU*(kEta + 6*tOmega)`); these are the same algebraic object up to canonicalization conventions. The residuals match. `engines_agree: true`.

## Freshness check

- SymPy script mtime: 2026-05-22 12:19; SymPy output mtime: 2026-05-22 12:23 → output newer than script.
- Mathematica script mtime: 2026-05-22 12:21; Mathematica output mtime: 2026-05-22 12:23 → output newer than script.

`outputs_fresh: true`. No `stale_output` finding.

## Adversarial probes that failed (and why)

I attempted the following attacks; none broke the audit:

1. **Algebraic independence of `eps_eta, eps_W, rho, Z_W`.** The four kernel ratios are not algebraically independent: from their microscopic definitions one derives `rho^2 sigma Z_W = eps_eta eps_W`. The script declares all four as independent `positive=True` symbols. If SymPy's substitution path applied to `c_UW^2 c_etaU^2` via the squares route (yielding `eps_eta eps_W K_U^2 Keta_eff K_W^eff/sigma`), the script would need this constraint to collapse to `rho^2 K_U^2 c_etaW^2`. **However**, the SymPy intermediates keep `(K_U c_etaW + c_UW c_etaU)^2` factored, so SymPy's `subs({c_UW*c_etaU: rho*K_U*c_etaW})` substitutes inside the square and never produces the bare `c_UW^2 c_etaU^2` factor. The constraint is therefore not required, and the identity holds for any positive values of the four ratios. Not a finding — the script's algebraic path is robust.

2. **`positive=True` on `rho`.** SymPy's `rho` is declared positive, stronger than the notes' "natural branch `1+rho > 0`". This is sufficient for the sign assertions (which only require `1+rho > 0` together with `Z_W, Lambda, delta_0 > 0`) and is consistent with the notes' branch choice. It does not introduce spurious cancellations. Not a finding.

3. **Missing `eps_eta < 1`, `eps_W < 1` constraints.** The placement formulas have `(1-eps_eta)` and `(1-eps_W)` in denominators, which would vanish at `eps_eta = 1` or `eps_W = 1`. The script's algebraic identities (substituted forms equal compact forms) hold as rational-function identities away from the divisor zero set; SymPy's `simplify` does not require sign information on `(1-eps_eta)` or `(1-eps_W)` to reach `0`. For the sign assertions, the residual is multiplied by `(1-eps_eta)^k (1-eps_W)^m` factors whose **squared** instances are manifestly positive; the residuals reduce to literal `±1` independent of where the ratios sit. The only sign claims that genuinely use `0 < eps < 1` are the prose lines at sympy:219-222 / wl:183-186, which are commentary, not assertions; the algebraic sign reductions hold without the bound. Not a finding.

4. **Docstring naming drift.** The SymPy docstring (line 3) says "Stage 21" and the Mathematica banner (line 33) says "STAGE 021"; the files live at `stage_038`. The notes file itself uses "Stage 21" interchangeably (e.g., "So Stage 21 turns the Stage-18/19 admissibility problem..."), and the paper card uses "Stage 038". This is stage-renumbering legacy across the project; the assertions and substitutions are intact. Not a `paper_misalignment` because no claim is affected.

5. **Variable name collision in Mathematica.** The `.wl` script reuses the name `a` for the Stage-037 amplitude `A` (line 47) while elsewhere using `aScale` for the gravitational length scale. Within the script the symbols are disjoint (the substitution `gNewton -> 20 lambda aScale^5 ...` does not touch `a`), but it is mildly confusing. No assertion misuses the name. Not a finding.

6. **Equivalence of the two product-law forms.** The script asserts both `R_target M_mix - 8 Lambda (1-eps_W)/pi^2 = 0` (compact form) and a chain through `NQ K_W^eff(1-eps_W)/mu_W` (microscopic form). The second uses `G -> 20 Lambda a^5 c^5 mu_W/(27 pi^2 c_s^5 K_W^eff)` and verifies `8 Lambda (1-eps_W)/pi^2 = (54 G c_s^5/(5 a^5 c^5)) K_W^eff (1-eps_W)/mu_W`. By hand: `(54/5)(G c_s^5/(a^5 c^5)) K_W^eff/mu_W = (54/5)(20 Lambda/(27 pi^2)) = 1080/(135 pi^2) Lambda = 8 Lambda/pi^2`. Holds. Not a finding.

## Verdict justification

Verdict `clean`. Stage 038's paper card and notes file claim (a) a six-ratio dimensionless reformulation of the Stage-037 continuum placement data, (b) three closed-form placement identities for `delta`, `M_mix`, `R_target` plus the `beta_0` corollary, (c) an exact product law with a microscopic-form equivalent, (d) invariance of the product under redistributions of `Z_W, rho, eps_eta`, and (e) nine signed parameter tendencies on the natural branch. Both engines verify all of (a)-(e) by symbolic substitution and simplification, with non-tautological residuals (each residual is built from microscopic objects that must rearrange across multiple identities to vanish). The Mathematica script reaches the same conclusions through engine-native simplification — `Expand[PowerExpand[Factor[...]]]` followed by `FullSimplify[Together[Expand[...]]]` — not by transliterating SymPy's intermediate algebra, and it does not bake any factored answer into a substitution rule. The transcripts are fresh.

I attempted six adversarial probes (algebraic independence of the four kernel ratios, positivity-of-rho strengthening, missing `eps < 1` constraints, docstring naming drift, Mathematica variable reuse, two-form equivalence in the product law) and none broke the audit. No `paper_misalignment` is present: every paper-card deliverable has a faithful script-side check, the script's "extras" (`beta_0`, derivative signs) are explicitly stated in the notes file that informed the card, and the only invariance check the card explicitly asks for is made manifest by the structural absence of `Z_W, rho, eps_eta` in the asserted product `8 Lambda (1-eps_W)/pi^2`.

## Self-test notes

I checked the relevant traps before finalizing:

1. **Variable independence** (trap 1): The four kernel ratios `eps_eta, eps_W, rho, Z_W` carry an algebraic redundancy `rho^2 sigma Z_W = eps_eta eps_W`. The SymPy substitution path avoids invoking the redundancy by keeping `(K_U c_etaW + c_UW c_etaU)^2` factored before substituting `c_UW*c_etaU` inside the square; treating the four ratios as independent positive reals is benign for the asserted identities.
2. **Symmetry/parity** (trap 2): No integrals over unbounded domains; the stage is pure rational-function algebra in six dimensionless ratios.
3. **Trivial-case substitution** (trap 3): At `eps_eta = 1/2, eps_W = 1/2, Z_W = 1, rho = 0, Lambda = 1`, the product law gives `R_target M_mix = 8(1/2)/pi^2 = 4/pi^2`, consistent with both sides. Plugging the same values into `dM/d eps_eta = 8 Z_W (1+rho)^2 / [pi^2 (1-eps_eta)^2 (1-eps_W)] = 8/(pi^2 (1/4)(1/2)) = 64/pi^2 > 0`, matching the asserted positive sign.
4. **Paper round-trip** (trap 5): Zero non-paper_misalignment fixes are prescribed because the verdict is clean; no risk of introducing new misalignment.
