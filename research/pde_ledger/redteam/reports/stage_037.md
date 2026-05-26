---
unit_id: 037
batch: III.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage037_continuum_kernel_extraction.md"]
  paper_appendix: present
---

# Audit unit 037 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_037.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage037_continuum_kernel_extraction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 52 only)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage037_continuum_kernel_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage037_continuum_kernel_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.txt`

## What the paper claims

The paper card's `\stagefield{Output}` reads: "Closed continuum formulas \eqref{eq:app-stage037-A-DK-continuum}--\eqref{eq:app-stage037-M-delta-continuum} and the stability gate \eqref{eq:app-stage037-stability}." Concretely the stage proves six explicit closed continuum-kernel formulas (boxed in the card and reproduced in section 5 of the notes): (i) `A = (K_U K_eta^eff - c_etaU^2)/(mu_eta K_U)`, (ii) `DeltaK_ax = pi^2 T_w/(mu_eta L^2)`, (iii) `alpha_mix = (K_U c_etaW + c_UW c_etaU)^2 / [mu_eta K_U (K_U K_W^eff - c_UW^2 sigma)]`, (iv) `beta_0 = (mu_W/mu_eta)(K_U c_etaW + c_UW c_etaU)^2 / (K_U K_W^eff - c_UW^2 sigma)^2`, (v) `M_mix = 8 (K_U c_etaW + c_UW c_etaU)^2 / [pi^2 (K_U K_eta^eff - c_etaU^2)(K_U K_W^eff - c_UW^2 sigma)]`, and (vi) `delta = pi^2 T_w K_U / [L^2 (K_U K_eta^eff - c_etaU^2)]`. The card also asserts the Schur-form rank-one wall self-energy `Sigma_wall(omega) = Xi(omega) I_2 + alpha(omega) v v^T` (eq. app-stage037-schur-form) and the stability inequalities `K_U K_eta^eff > c_etaU^2` and `K_U K_W^eff > c_UW^2 sigma`. Beyond the boxed `.tex`, the notes (section 4) further fix the explicit closed forms for `Xi(omega) = g_U^2/A_U(omega)` and `alpha(omega) = g_B^2/A_phi + (A_U g_W + g_R g_U)^2 / (A_U Delta_UW)`, with `Delta_UW = A_U A_W - g_R^2 sigma`. The Part III appendix row (line 52) summarises the same deliverables in a single line.

## What the script claims to verify

The scripts' docstring/banner ("STAGE 20 / 020 CONTINUUM-KERNEL EXTRACTION") states the intent to derive the six reduced quantities from the explicit continuum kernel and verify the Schur factorization. The SymPy script's load-bearing assertions are: (1) the N/N+D/N overlaps reproduce `kappa_0 = 2 sqrt(2)/pi`, `kappa_1 = -4/(3 pi)`, `sigma = 88/(9 pi^2)`; (2) `C B^{-1} C^T` equals the closed-form ansatz `Xi I_2 + alpha v v^T` with `Xi = g_U^2/A_U` and `alpha = g_B^2/A_phi + (A_U g_W + g_R g_U)^2/(A_U Delta_UW)`, where the full 2x2 matrix residual is asserted to vanish (line 144); (3) the abstract definitions `A = K_0 - g_U^2/Omega_U^2`, `Delta_0 = Omega_U^2 Omega_W^2 - g_R^2 sigma`, `Chi = Omega_U^2 g_W + g_R g_U`, `beta_0 = Chi^2/Delta_0^2`, `alpha_mix = Chi^2/(Omega_U^2 Delta_0)`, `M_mix = 8 alpha_mix/(pi^2 A)`, `delta = DeltaK_ax/A` each reduce to the boxed closed forms (i)-(vi) (lines 186-192). Section 5 re-verifies A and Delta_0 closed forms underwriting the printed stability inequalities. The Mathematica script asserts a parallel ledger; the one substantive deviation is that its Schur step verifies only the *rank-one+identity shape* of `Sigma_wall` via solve-then-consistency on a single entry, instead of asserting the closed forms for `Xi` and `alpha` as the SymPy script does.

## Paper - script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `A` formula (eq. app-stage037-A-DK-continuum) | sympy L186, L210; math L167, L168-174, L201 | match |
| `DeltaK_ax` formula (eq. app-stage037-A-DK-continuum) | sympy L102 (definition coincides with closed form); math L84 (same) | match |
| `alpha_mix` formula (eq. app-stage037-alpha-beta-continuum) | sympy L190; math L178 | match |
| `beta_0` formula (eq. app-stage037-beta-continuum) | sympy L189; math L177 | match |
| `M_mix` formula (eq. app-stage037-M-delta-continuum) | sympy L191; math L179 | match |
| `delta` formula (eq. app-stage037-M-delta-continuum) | sympy L192; math L180-187 | match |
| Schur factorization `Sigma_wall = Xi I_2 + alpha v v^T` with notes-stated closed forms `Xi = g_U^2/A_U`, `alpha = g_B^2/A_phi + (A_U g_W + g_R g_U)^2/(A_U Delta_UW)` | sympy L137-144 (full 2x2 matrix vs closed-form ansatz); math L114-125 (solves for `xi`, `alpha` from two entries of the computed `sigmaWall`; checks (2,2) consistency only; never asserts the recovered values match the notes' closed forms) | partial (sympy match; math partial) |
| Stability inequalities (eq. app-stage037-stability) | sympy L208-219 and math L199-211 print the inequalities; both re-verify A and Delta_0 closed forms but do not explicitly assert the equivalence `A > 0  <=>  K_U K_eta_eff > c_etaU^2` | partial (algebraically trivial corollary; not load-bearing) |

Dominant pattern: match. The single substantive partial is the Mathematica Schur check (see F1). The stability-gate "partial" is trivial: given the verified closed form `A = (K_U K_eta_eff - c_etaU^2)/(mu_eta K_U)` and the script's positivity declarations on `mu_eta` and `K_U`, the iff with `K_U K_eta_eff > c_etaU^2` is immediate; no new assertion required.

`paper_alignment` set to `aligned` - every paper-side deliverable is exercised in at least one engine. The notes' explicit closed forms `Xi = g_U^2/A_U` and the full `alpha(omega)` expression are directly cross-checked only in SymPy.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 65-73 | `expect_zero` mode normalization, orthogonality, BCs | basis prerequisites for eq. app-stage037-schur-form | yes |
| A2 | sympy | 83-85 | `expect_zero(kappa0 - 2 sqrt(2)/pi)`, etc. | Inputs line of stage card (kappa_0, kappa_1, sigma) | yes |
| A3 | sympy | 144 | `expect_zero(Sigma - Sigma_expected)` (full 2x2 matrix) | Schur form eq. app-stage037-schur-form + notes' closed Xi/alpha | yes |
| A4 | sympy | 186 | `expect_zero(A - A_expected)` | A formula | yes |
| A5 | sympy | 187 | `expect_zero(Delta0 - Delta0_expected)` | underwrites alpha_mix/beta_0 denominators | yes |
| A6 | sympy | 188 | `expect_zero(Chi - Chi_expected)` | Chi closed form (notes section 5) | yes |
| A7 | sympy | 189 | `expect_zero(beta0 - beta0_expected)` | beta_0 formula | yes |
| A8 | sympy | 190 | `expect_zero(alpha_mix - alpha_mix_expected)` | alpha_mix formula | yes |
| A9 | sympy | 191 | `expect_zero(M_mix - M_mix_expected)` | M_mix formula | yes |
| A10 | sympy | 192 | `expect_zero(delta - delta_expected)` | delta formula | yes |
| A11 | sympy | 208-214 | A and Delta_0 closed forms re-asserted | stability inputs (redundant with A4/A5) | yes (redundant) |
| B1 | math | 51-59 | mode normalization/orthogonality/BCs | basis prerequisites | yes |
| B2 | math | 68-70 | `expectZero` kappa_0/kappa_1/sigma | Inputs line | yes |
| B3 | math | 125 | `expectZero(sigmaWall[[2,2]] - xiSolved - alphaSolved*kappa1^2)` | Schur-form shape only (xi/alpha recovered from (1,1) and (1,2); (2,2) checked for consistency); does NOT verify `Xi = g_U^2/A_U` or the closed `alpha(omega)` | partial - see F1 |
| B4 | math | 167 | `expectZero(a - aDerived)` where `aDerived = FullSimplify[Together[k0 - gUCont^2/omegaU2]]` | tautological as a standalone check; covered substantively by B5 | partial (tautological in isolation, redundant given B5) |
| B5 | math | 168-174 | numerator/denominator of A matches `(kU*kEtaEff - cEtaU^2)/(muEta*kU)` | A formula | yes |
| B6 | math | 175 | Delta_0 closed form | underwrites alpha_mix/beta_0 | yes |
| B7 | math | 176 | Chi closed form | Chi (notes section 5) | yes |
| B8 | math | 177 | beta_0 closed form | beta_0 formula | yes |
| B9 | math | 178 | alpha_mix closed form | alpha_mix formula | yes |
| B10 | math | 179 | M_mix closed form | M_mix formula | yes |
| B11 | math | 180-187 | delta numerator/denominator match closed form | delta formula | yes |
| B12 | math | 199-205 | A and Delta_0 closed forms re-asserted | stability inputs (redundant with B5/B6) | yes (redundant) |

Note on B4: `aDerived` (math line 150) equals `FullSimplify[Together[k0 - gUCont^2/omegaU2]]`, which is identical to `a` (math line 142) up to `Together`. The assertion `a - aDerived == 0` is therefore algebraically tautological *as a standalone check*. It is not filed as a `tautological_check` finding because B5 (the next line) is the substantive A-closed-form assertion; B4 is dead scaffolding rather than a load-bearing check that hides a defect.

## Findings

### F1 - insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.wl:114-125`

**What's wrong:**

The Mathematica script's Schur-complement check does not verify the closed forms for `Xi` and `alpha` that the notes (and the SymPy script) treat as load-bearing. The notes (section 4) state:

> `Xi(omega) = g_U^2 / A_U(omega),`
>
> `alpha(omega) = g_B^2 / A_phi(omega) + ( A_U(omega) g_W + g_R g_U )^2 / [ A_U(omega) Delta_UW(omega) ],`
>
> `Delta_UW(omega) = A_U(omega) A_W(omega) - g_R^2 sigma.`

The SymPy script verifies this directly (lines 137-144):

> ```
> Xi = sp.simplify(gU**2 / A_U)
> alpha = sp.simplify(gB**2 / A_phi + (A_U * gW + gR * gU) ** 2 / (A_U * Delta_UW))
> Sigma_expected = sp.simplify(Xi * sp.eye(2) + alpha * sp.Matrix([
>     [kappa0**2, kappa0 * kappa1],
>     [kappa0 * kappa1, kappa1**2],
> ]))
> expect_zero("Sigma_wall - [Xi I + alpha v v^T]", Sigma - Sigma_expected)
> ```

The Mathematica script (lines 114-125) does the opposite: it computes `sigmaWall = cMat . LinearSolve[bMat, Transpose[cMat]]`, then *solves backwards* for `xi` and `alpha` from two entries of the computed Sigma, and only verifies a single consistency equation on the remaining (2,2) entry:

> ```
> alphaSolved = FullSimplify[
>   alpha /. First[Solve[sigmaWall[[1, 2]] == alpha*kappa0*kappa1, alpha]], ...];
> xiSolved = FullSimplify[
>   xi /. First[Solve[sigmaWall[[1, 1]] == xi + alphaSolved*kappa0^2, xi]], ...];
> expectZero["Sigma_wall (2,2) consistency with ansatz",
>            sigmaWall[[2, 2]] - xiSolved - alphaSolved*kappa1^2];
> ```

This check confirms only that `sigmaWall` has the rank-one-plus-identity *shape* — one consistency equation on three independent entries of a 2x2 symmetric matrix. It does NOT verify that the recovered `xiSolved` equals `g_U^2/A_U` nor that the recovered `alphaSolved` equals the closed form in the notes. The saved transcript (line 52 of `.../mathematica_audit.txt`) prints `alpha (recovered) = (-88*aU*gB^2*gR^2 + 9*(aPhi*gR^2*gU^2 + aU*(aU*aW*gB^2 + 2*aPhi*gR*gU*gW + aPhi*aU*gW^2))*Pi^2)/(aPhi*aU*(-88*gR^2 + 9*aU*aW*Pi^2))` — visually plausible as the closed form once `sigma = 88/(9*Pi^2)` is folded in, but the script never asserts that equality. A regression in which `sigmaWall` had the correct rank-one structure but with a wrong `alpha` numerator (e.g., a missing factor or sign-flipped cross term) would still pass the existing (2,2)-consistency check.

**Why this matters:**

The Schur factorization is one of the load-bearing structural claims of stage 037 (boxed equation app-stage037-schur-form), and the notes' explicit closed forms for `Xi` and `alpha` are the *content* of that claim. The two-engine policy requires both engines to verify the same content. Today only SymPy does; Mathematica verifies only the shape. The shape alone is a strictly weaker statement than the closed forms, and the existing Mathematica check would not catch a defect in either `Xi` or `alpha` that preserved the rank-one+identity structure.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.wl`, after computing `sigmaWall` (line 114), add closed-form `expectZero` checks that mirror the SymPy assertion at line 144 of the `.py`. Define independently:
- `xiExpected = gU^2 / aU`
- `deltaUW = aU*aW - gR^2 * sigma`
- `alphaExpected = gB^2/aPhi + (aU*gW + gR*gU)^2 / (aU*deltaUW)`
- `sigmaWallExpected = xiExpected * IdentityMatrix[2] + alphaExpected * {{kappa0^2, kappa0*kappa1}, {kappa0*kappa1, kappa1^2}}`

Then call `expectMatrixZero["Sigma_wall - [Xi I + alpha v v^T]", sigmaWall - sigmaWallExpected]`. The existing solve-then-consistency check on line 125 may remain (it provides a second route) but the new closed-form matrix-equality check must be present and must pass.

**Verification:**

After the patch, the Mathematica transcript at `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.txt` must contain a new `PASS: Sigma_wall - [Xi I + alpha v v^T]` line, demonstrating the full 2x2 matrix residual simplifies to zero against the closed-form ansatz. The script must still exit 0.

## Independent-derivation check (Mathematica)

The `.wl` script is structurally a close port of the `.py` script: same banner/subbanner functions, same five sections in the same order, same intermediate variable names (`k0`, `omegaU2`, `chi`, `delta0`, `mMix`, `beta0`, `alphaMix`, `delta`, etc.), and same input definitions. Sections 1, 2, 4, and 5 are essentially direct translations with no independent derivation step. Section 3 *does* differ in approach: SymPy declares Xi/alpha closed forms and compares the full Schur matrix against the ansatz, whereas Mathematica computes Sigma and solves backwards for Xi/alpha from two of its entries. The difference is the one piece of methodological independence in the file. It is not strong enough to constitute a fully independent derivation, but the broader two-engine cross-check is still meaningful because the final closed-form algebraic verifications in Section 4 use a different Mathematica idiom (`Numerator[Together[...]]*expected_denom - Denominator[Together[...]]*expected_numer`) than the SymPy `simplify(diff) == 0` style — that is, Mathematica is matching rational numerators/denominators rather than testing `FullSimplify[a - aExpected] == 0` on the unsimplified difference.

I do not file a `mathematica_transliteration` finding because (a) the Schur step uses a genuinely different algorithmic route (solve-then-consistency rather than declare-then-equate), (b) the rational num/den matching idiom in Section 4 is methodologically distinct, and (c) the substantive complaint about the Schur step is captured by F1 (`insufficient_verification`), which is the more precise category — the Mathematica Schur check is not a transliteration; it is a weaker check.

## Engine cross-check

Both engines verify the same set of closed-form identities for A, DeltaK_ax, alpha_mix, beta_0, M_mix, delta, and produce visually identical final expressions (modulo Mathematica's `Pi^2` placement and the order of summands in Delta_0; both forms are algebraically equal — confirmed by their respective `expectZero` checks against the same closed expressions). Sample comparisons from the saved transcripts:
- SymPy `A = (K_U*(K_eta + 6*T_Omega) - c_etaU**2)/(K_U*mu_eta)`
- Mathematica `A = (-cEtaU^2 + kU*(kEta + 6*tOmega))/(kU*muEta)` (identical up to ordering)
- SymPy `delta = pi**2*K_U*T_w/(L**2*(K_U*(K_eta + 6*T_Omega) - c_etaU**2))`
- Mathematica `delta = (kU*Pi^2*tw)/(ell^2*(-cEtaU^2 + kU*(kEta + 6*tOmega)))` (identical up to ordering)

Same agreement for Delta_0, Chi, beta_0, alpha_mix, M_mix. No engine disagreement.

The only mismatch in *scope of verification* is the Schur step (see F1): SymPy verifies the closed forms for `Xi` and `alpha`; Mathematica verifies only the rank-one+identity shape.

Both saved outputs are newer than their respective scripts (sympy script Apr 1 12:39, sympy output May 22 12:17; math script May 22 12:16, math output May 22 12:18), so no `stale_output` finding.

## Verdict justification

The paper card and notes prescribe six closed continuum-kernel formulas plus the Schur factorization plus two stability inequalities. Both engines verify the six formulas non-tautologically (the abstract definitions reduce to the boxed closed forms after substitution through mass-normalized variables; the assertions could fail under a sign flip, missing-coupling, or wrong-sigma error). The stability inequalities follow trivially from the verified closed forms together with the declared positivity of `mu_eta`, `K_U`; that does not need an extra assertion. The one substantive defect is that the Mathematica script verifies only the *structural shape* of the Schur factorization, not the explicit closed forms for `Xi` and `alpha` that the notes treat as load-bearing — SymPy covers this but Mathematica does not. Verdict: findings; one `insufficient_verification` for the Mathematica Schur step. Not `CRITICAL_DOWNSTREAM`: the closed forms for A, DeltaK_ax, alpha_mix, beta_0, M_mix, delta are correctly verified by both engines (these are the actual carry-forward data into stages 038 and beyond), so the fix tightens the second-engine policy without changing any quantity downstream stages consume.

The internal labels "STAGE 20" / "STAGE 020" in the script docstrings and banners are a holdover from pre-renumbering; the substance matches stage 037, so this is not a `paper_misalignment` finding (no math defect, only a cosmetic stale label that would be appropriate for a future polish pass).

I attacked: (a) tautological structure in the A-closed-form check (B4 *is* tautological standalone, but B5 is the load-bearing check that establishes the closed form; B4 is dead scaffolding, not a defect); (b) sign conventions in the B matrix off-diagonals (consistent with the Lagrangian's Hessian; the SymPy full-matrix `Sigma - Sigma_expected = 0` check would catch any sign flip); (c) symbol-domain errors (couplings are `real=True`, not `positive=True` — correct, since couplings may be negative and only squares enter the stability gates); (d) missing branches (the static conservative branch `Pi_out(0) = 0` is explicitly invoked in the notes; the script does not claim to cover the dynamic branch); (e) stale output (both outputs younger than their scripts); (f) hidden tautologies via `simplify` under aggressive assumptions (positivity is declared only on mass densities and stiffnesses, which is physically correct; FullSimplify with `Assumptions -> $Assumptions` does not have enough leverage to manufacture a false PASS on the closed-form identities). None of these produced a finding beyond F1.

## Self-test notes

For F1's required change I mentally substituted concrete profiles: `xiExpected = gU^2/aU` and `alphaExpected = gB^2/aPhi + (aU*gW + gR*gU)^2 / (aU*(aU*aW - gR^2*sigma))` are rational expressions in the same scalars that `sigmaWall` depends on; under `FullSimplify` the 2x2 difference `sigmaWall - sigmaWallExpected` reduces to a 2x2 zero matrix (SymPy already demonstrates this with identical algebra). Variable-independence and parity traps do not apply (no derivatives or symmetric-domain integrals in the proposed new assertion). The directive's target path correctly names the `.wl` file under `mathematica/`. Paper round-trip: the closed forms for `Xi` and `alpha` are stated verbatim in `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage037_continuum_kernel_extraction.md` section 4, so the new check has a clear paper-side anchor and does not introduce a new `paper_misalignment`.
