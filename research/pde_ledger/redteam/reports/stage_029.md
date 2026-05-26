---
unit_id: 029
batch: II.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage029_dynamic_loading.md"]
  paper_appendix: present
---

# Audit unit 029 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_029.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage029_dynamic_loading.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row at line 48; `\input` at line 96)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage029_dynamic_loading_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.txt`

## What the paper claims

Stage 029 gives the loading coefficient `alpha` a microscopic origin by computing the exact Schur complement of the coupled wall/BdG/Maxwell/mixed operator. The paper's `\stagefield{Output}` (line 106) is verbatim: *"Stage~029 outputs the Schur-complement split \eqref{eq:app-stage029-schur-split}, the static branch data \eqref{eq:app-stage029-static-data}, the outgoing transfer coefficient \eqref{eq:app-stage029-beta0}, and the selected odd coefficient \eqref{eq:app-stage029-selected-odd}."* The four deliverables are:

1. **Schur split:** `Sigma_wall(omega) = Xi(omega) I_2 + alpha(omega) v v^T` with `Xi = lambda_U^2 / A_U` and `alpha = lambda_B^2 / A_phi + (A_U lambda_W + lambda_R lambda_U)^2 / (A_U Delta_UW)`; `Delta_UW = A_U A_W - lambda_R^2 sigma`; `sigma = 88/(9 pi^2)` (eqs `Xi-alpha`, `DeltaUW`).
2. **Static branch data:** `Xi_0 = lambda_U^2/Omega_U^2`, `Delta_0 = Omega_U^2 Omega_W^2 - lambda_R^2 sigma`, `alpha_0 = lambda_B^2/varpi^2 + (Omega_U^2 lambda_W + lambda_R lambda_U)^2/(Omega_U^2 Delta_0)`, and the conservative static matrix `K_eff^(0)` (eq `static-data`).
3. **Outgoing transfer coefficient:** `beta(omega) = (A_U lambda_W + lambda_R lambda_U)^2 / Delta_cons^2` with `Delta_cons = A_U(Omega_W^2 - omega^2) - lambda_R^2 sigma`, and `beta_0 = (Omega_U^2 lambda_W + lambda_R lambda_U)^2 / Delta_0^2 >= 0` (eqs `beta`, `beta0`).
4. **Selected odd coefficient:** `delta D_-^odd(omega) = -i (a^5/(27 c_s^5)) beta_0 (v.e_-)^2 omega^5 + O(omega^7)` (eq `selected-odd`).

The notes provide the same content plus the additional named result `alpha_crit = (K_0 - Xi_0)(K_1 - Xi_0) / [(K_1 - Xi_0) kappa_0^2 + (K_0 - Xi_0) kappa_1^2]` (notes section 3.1). The `\stagefield{Checks}` (lines 108-113) emphasize: scalar U-sector -> I_2, mixed/support -> vv^T; positivity of alpha_0 under stable denominators; beta_0 as square-over-squared-determinant; projection onto e_- inserts (v.e_-)^2.

## What the script claims to verify

Both engines work in the same reduced operator (wall (q_0, q_1), internal block (u_0, u_1, W, phi)) with the same coupling pattern. Sympy uses the direct 4x4 Schur formula `C^T M_int^(-1) C`. Mathematica derives `Sigma_wall` by partial elimination - Sherman-Morrison-like inverse of the (U, W) block plus the phi contribution as a separate vv^T piece - and asserts equality with the closed-form `Xi I + alpha vv^T`. The static substitutions omega->0, Pi->0 give the printed expressions for Xi_0, Delta_0, alpha_0. The first-order Pi expansion verifies the closed-form `beta(omega)`. The `Pi -> i Gamma_port omega^5` substitution yields a series whose omega^5 coefficient is extracted and shown equal to `Gamma_port * beta(omega=0)`. `kappa_sel^2` is computed via Hellmann-Feynman `-d lambda_- / d alpha`; Mathematica additionally constructs `e_-` directly from `Eigensystem` and asserts that `(v.e_-)^2` equals the Hellmann-Feynman expression. Sympy only checks weak (alpha=0 -> kappa_0^2) and strong (alpha->infty -> sigma) limits.

## Paper <-> script cross-check

| # | Paper deliverable | Script-side check | Status |
|---|---|---|---|
| D1 | Schur split `Sigma = Xi I + alpha vv^T` | sympy line 127; mathematica line 116 | match |
| D2 | Static branch data Xi_0, Delta_0, alpha_0 | sympy lines 141-143 (printed only, no `expect_zero` against literal closed form); mathematica lines 121-123 (printed only) | partial - values are computed by substitution from the verified Sigma formula, so D2 is structurally implied by D1, but no explicit anchoring assertion exists |
| D3 | Outgoing transfer beta(omega) and beta_0 >= 0 | sympy line 216 asserts `beta(omega) = (A_U lambda_W + lambda_R lambda_U)^2/Delta_cons^2`; mathematica lines 175, 188 assert both `beta(omega)` and the explicit `beta_5 = Gamma_port (Omega_U^2 lambda_W + lambda_R lambda_U)^2/Delta_0^2` | match (sympy carries beta_0 implicitly through `beta_clean.subs(omega, 0)`; mathematica carries it explicitly) |
| D4 | Selected odd `delta D_-^odd = -i Gamma_port beta_0 (v.e_-)^2 omega^5` | sympy lines 263-265 print `odd_projection = -i beta_5 kappa_theta^2 omega^5` (with `kappa_theta^2 = (q.v)^2` as a function of theta, NOT specialized to theta_-) and prints `kappa_sel^2` separately; no `expect_zero` ties the two together or against the paper's combined formula. Mathematica line 228-229 forms an analogous template, also unasserted. | partial - building blocks present and verified, but the paper's combined claim is never an asserted identity in either engine |
| D5 (paper Checks) | beta_0 >= 0 inequality | not asserted; structurally true because formula is real^2 / real^2 | partial (structural) |
| extra | `alpha_crit = K0t*K1t/(K1t kappa_0^2 + K0t kappa_1^2)` | sympy line 194; mathematica line 166 | extra - verified by both engines, but not in paper card `\stagefield{Output}` (only in notes section 3.1) |

Dominant pattern: D1 and D3 are exactly anchored; D2 and D4 are implicit; D5 is structural; one extra (alpha_crit) lives only in the notes. Setting `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 127 | `expect_zero("Sigma - (Xi I + alpha vv^T)", Sigma - Sigma_expected)` | D1 (Schur split) | yes |
| A2 | sympy | 155 | `expect_zero("DeltaK_tilde - DeltaK_bare", DeltaK - DeltaK_ax)` | typo guard (script self-acknowledges via comment lines 152-154) | no (tautological by construction; explicitly flagged) |
| A3 | sympy | 184 | `expect_zero("dE/dtheta - stationarity/2", dE - stationarity/2)` | D2 stationarity / notes section 3 angle law | yes |
| A4 | sympy | 194 | `expect_zero("alpha_crit - expected", alpha_crit - alpha_crit_expected)` | extra (alpha_crit from notes section 3.1, not from paper card Output) | yes |
| A5 | sympy | 216 | `expect_zero("beta - clean transfer factor", beta - beta_clean)` | D3 (beta(omega) formula) | yes |
| A6 | sympy | 220 | `expect_zero("alpha - (alpha_cons + beta*Pi) at O(Pi)", ...)` | first-order outgoing expansion (paper eq `alpha-out`) | yes |
| A7 | sympy | 232 | `expect_zero("extracted beta_5 - expected beta_5", ...)` | D3 (Pi = i Gamma omega^5 projection); implicitly D4 (beta_0 via `beta_clean.subs(omega, 0)`) | yes |
| A8 | sympy | 278 | `expect_zero("weak-loading kappa_sel^2 -> kappa0^2", ...)` | D4 building block (e_- at alpha=0 = (1,0)) | partial - limit only |
| A9 | sympy | 279 | `expect_zero("strong-loading kappa_sel^2 -> sigma", ...)` | D4 building block (e_- at alpha->infty = v/||v||) | partial - limit only |
| M1 | math | 116 | `expectMatrixZero["Sigma_seq - (Xi I + alpha vv^T)", ...]` | D1 (Schur split, independent derivation) | yes |
| M2 | math | 117-119 | three `expectZero` for sigma, xi, eta literals | building blocks for v.v, etc. | yes |
| M3 | math | 132 | `expectZero["DeltaK_tilde - DeltaK_ax", ...]` | typo guard | no (tautological; explicitly flagged) |
| M4 | math | 148, 150 | `expectZero` for stationarity and stationarity-at-theta_- | D2 angle law | yes |
| M5 | math | 166 | `expectZero["det(alpha_crit)", ...]` | extra (notes section 3.1) | yes |
| M6 | math | 175, 179 | `expectZero` for beta(omega) and first-order Pi expansion | D3 | yes |
| M7 | math | 188, 189 | `expectZero` for beta_5 vs explicit `(Omega_U^2 lambda_W + lambda_R lambda_U)^2/Delta_0^2` AND series-coefficient extraction | D3 and direct D4 building block beta_0 closed form | yes |
| M8 | math | 223 | `expectZero["kappa_sel^2 closed-form vs eigenvector projection", kappaSelSq - kappaSelSqDirect]` | D4 (Hellmann-Feynman == direct projection onto e_-) | yes - independent cross-check |
| M9 | math | 225, 226 | weak/strong limit checks (mirror A8, A9) | D4 building block | partial - limit only, but redundant given M8 |

Sympy is missing the M8 cross-check (closed-form == eigenvector projection), and neither engine asserts D4 in fully-assembled form against the paper's combined formula.

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** notes_contradicts_script (docstring / banner refer to obsolete stage number)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:3` (first line of docstring "moving_throat_pde_stage12_dynamic_loading_sympy_audit.py")
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:5` ("SymPy audit for Stage 12 of the moving-throat PDE program.")
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:33` (`banner["STAGE 012 — DYNAMIC LOADING"]`)

**What's wrong:**
Paper card title (`paper/stages/stage_029.tex:1`): *"Stage 029: Coupled Response Operator"*; see also `\input{stages/stage_029}` at `paper/appendices/stage_appendix_part02.tex:96`. The sympy docstring states this is the audit for "Stage 12"; the Mathematica banner prints `STAGE 012 — DYNAMIC LOADING`. The closing line of the Mathematica script (line 232) correctly says "Stage 029 Mathematica audit passed." — so the file knows both numbers exist. The notes file similarly mixes legacy "Stage 11"/"Stage 12" references with the current "Stage 029" title; the notes are out of scope for this auditor (prose, read-only).

**Why this matters:**
Pure labeling drift from an earlier numbering scheme. The math is unaffected. The risk is operational: a future reader auditing by docstring or banner search may believe the script does not cover stage 029.

**Required change:**
Update sympy docstring lines 3 and 5 to reference "Stage 029" instead of "Stage 12". Update the Mathematica banner string at line 33 from `"STAGE 012 — DYNAMIC LOADING"` to `"STAGE 029 — DYNAMIC LOADING"`. Do not edit the notes file. Do not touch comments inside the script that reference physical results from "Stage 11"/"Stage 12" descriptively (those are content references, not file-identification labels).

**Verification:**
After Codex applies, the banner output (mathematica) should read `STAGE 029 — DYNAMIC LOADING`, and the sympy docstring header should reference Stage 029. Existing `Stage 029 Mathematica audit passed.` line is already correct.

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:246-283` (function `selected_mode_projection`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:191-229`

**What's wrong:**
Paper card eq `selected-odd` (lines 100-104) claims `delta D_-^odd(omega) = -i (a^5/(27 c_s^5)) beta_0 (v.e_-)^2 omega^5 + O(omega^7)`. Neither engine asserts this combined identity. Both compute the pieces but never bind them:

- sympy line 263: `odd_projection = -sp.I * beta5 * kappa_theta_sq * omega**5`, where `kappa_theta_sq = (q.T*(v*v.T)*q)[0]` is a generic theta-dependent quantity (line 258), *not* specialized to theta_-, and `beta5 = beta_clean(omega=0) * Gamma_port = Gamma_port * beta_0`. The script prints `odd_projection` and prints `kappa_sel_sq` (Hellmann-Feynman) separately, but never asserts e.g. `odd_projection_at_theta_minus - (-i * Gamma_port * beta_0 * kappa_sel^2(alpha_0) * omega^5) == 0`.
- mathematica line 228 forms an analogous `oddProjection = -I*beta5*kappaSelSq*omega^5` and only `Print`s it (line 229) without `expectZero`.

So both engines verify the building blocks ((v.e_-)^2 via two routes in Mathematica; one route in sympy; beta_0 closed form in Mathematica; beta_5 series extraction in both), but the combined paper claim is left implicit.

**Why this matters:**
Without an asserted equality binding the printed `odd_projection` to the paper's symbolic form, a future edit to one of the building blocks (e.g. a sign or factor in beta_5 or kappa_sel^2) could pass the existing assertions while quietly changing the combined paper formula. The paper claim is the load-bearing one; the building blocks are intermediate.

**Required change:**
Add one assertion per engine that ties the building blocks to the paper formula.

In sympy `selected_mode_projection`, after the current line 275:
```
delta_D_paper = (
    -sp.I * Gamma_port
    * (Omega_U**2 * lambda_W + lambda_R * lambda_U)**2
    / (Omega_U**2 * Omega_W**2 - lambda_R**2 * sigma)**2
    * kappa_sel_sq
    * omega**5
)
delta_D_script = sp.simplify(-sp.I * beta5 * kappa_sel_sq * omega**5)
expect_zero(
    "delta D_-^odd (script) - delta D_-^odd (paper formula)",
    delta_D_script - delta_D_paper,
)
```
(`beta5` and `kappa_sel_sq` are already defined; this is a trivial equality once `beta5` matches the paper closed form, which is the case by construction at line 224.)

In Mathematica, after line 228:
```
deltaDPaper = -I*GammaPort
              * (OmegaU^2*lambdaW + lambdaR*lambdaU)^2 / delta0^2
              * kappaSelSq
              * omega^5;
deltaDScript = -I*beta5*kappaSelSq*omega^5;
expectZero[
  "delta D_-^odd (script) - delta D_-^odd (paper formula)",
  deltaDScript - deltaDPaper
];
```
This binds the paper's eq `selected-odd` to the script-computed combination, anchoring deliverable D4 directly.

**Verification:**
New `expect_zero`/`expectZero` line appears in each script. After `redteam exec-sympy 029` and `redteam exec-mathematica 029`, the new line emits `= 0` / `PASS`. The saved outputs include the new line.

### F3 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:268-279`

**What's wrong:**
The Mathematica script at line 223 asserts the cross-check `kappa_sel^2 closed-form vs eigenvector projection`, comparing the Hellmann-Feynman `-d lambda_- / d alpha` expression against `(v.e_-_hat)^2` where `e_-_hat` is constructed directly from `Eigensystem[K_eff(al)]`. This is the genuine independent verification that the Hellmann-Feynman shortcut produces the same `(v.e_-)^2` the paper writes.

The sympy script computes only the Hellmann-Feynman form (`kappa_sel_template`, line 272) and then checks weak/strong limits (lines 278-279). It does NOT construct `e_-` explicitly via SymPy's eigenvector routines and project `v` onto it. So sympy alone could not catch a Hellmann-Feynman sign error or a wrong-eigenvalue branch identification - the limit checks pass for either branch up to sign in the alpha=0 limit (both kappa_0^2 and sigma are positive in the limits), but the formula being correct *between* the limits is unverified by sympy alone.

**Why this matters:**
The selected odd coefficient (D4) carries the sign that makes the eigenvalue's `-i omega^5`-shift damping rather than amplification. A wrong-branch `kappa_sel^2` would still satisfy the two limit checks at alpha->0 and alpha->infty but be wrong in the interior. Mathematica catches this via M8; sympy does not - so on sympy alone, the verification is asymmetric.

**Required change:**
Add a SymPy direct-eigenvector projection comparable to Mathematica's M8. In `selected_mode_projection` after line 275, add roughly:
```
K_eff_al = sp.Matrix(
    [
        [K0t - al * kappa0**2, -al * kappa0 * kappa1],
        [-al * kappa0 * kappa1, K1t - al * kappa1**2],
    ]
)
eigdata = K_eff_al.eigenvects()
vec_lo = None
for ev, _, vecs in eigdata:
    if sp.simplify(ev - lam_minus_template) == 0:
        vec_lo = vecs[0]
        break
assert vec_lo is not None, "lower eigenvector not found"
norm_sq = sp.simplify((vec_lo.T * vec_lo)[0])
kappa_sel_sq_direct = sp.simplify(((vec_lo.T * v)[0])**2 / norm_sq)
expect_zero(
    "kappa_sel^2 closed-form vs eigenvector projection",
    sp.simplify(kappa_sel_template - kappa_sel_sq_direct),
)
```
This mirrors mathematica lines 201-223 in SymPy.

**Note for Codex (self-test):** `sp.Matrix.eigenvects()` may return eigenvalues in either order. The loop's match against `lam_minus_template` (already defined at line 271) picks the correct branch. The squared projection `((vec.T * v)[0])**2 / (vec.T * vec)[0]` is sign-invariant in the eigenvector's overall scale. If `eigenvects()` returns a `CRootOf` form or stalls, fall back to constructing the eigenvector analytically from the nullspace of `K_eff_al - lam_minus_template * sp.eye(2)`:
```
N = (K_eff_al - lam_minus_template * sp.eye(2)).nullspace()
vec_lo = N[0]
```
Either path is acceptable; the equality is what matters.

**Verification:**
New `expect_zero` line emits `= 0` after `redteam exec-sympy 029`. The new sympy output file contains this new line under SECTION IV.

### F4 — paper_misalignment

**Severity:** low
**Subtype:** paper_missing_script_claim
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_029.tex:106` (`\stagefield{Output}` — does not mention `alpha_crit`)
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_029.tex:108-113` (`\stagefield{Checks}` — does not mention `alpha_crit`)
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage029_dynamic_loading.md:173-186` (section 3.1 - defines and discusses `alpha_crit`)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:189-194` (asserts `alpha_crit - expected == 0`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:159-166` (asserts `det(alpha_crit) == 0`)

**What's wrong:**
Both scripts verify a closed form for the refined softening threshold `alpha_crit = K0t*K1t / (K1t kappa_0^2 + K0t kappa_1^2)`. This appears in the notes section 3.1 but is **not** listed in the paper card's `\stagefield{Output}` (line 106) nor in `\stagefield{Checks}` (lines 108-113). The card's Output enumerates only Schur split, static branch data, beta_0, and selected odd coefficient. The scripts therefore verify an extra deliverable that the paper card does not state.

**Why this matters:**
Either the paper card should mention `alpha_crit` (since the notes treat it as a substantive result of the stage), or the script should drop the `alpha_crit` assertion as out-of-scope. Direction of resolution is the user's call. Without resolution, a future reader of the paper card would not know `alpha_crit` had been verified, and a future editor of the card could not gauge whether the script-side check needs to be retained or trimmed.

**Required change:**
See directive `## Resolve before fix_loop` block. Codex does not auto-resolve.

**Verification:**
After user resolution, either: (a) paper card line 106 expands `\stagefield{Output}` to mention `alpha_crit` (notes-source authority preserved), or (b) the script's `alpha_crit` assertion is removed (sympy lines 189-194; mathematica lines 159-166). No auto-fix is applied until the user picks a direction.

## Independent-derivation check (Mathematica)

The Mathematica script is **not** a transliteration of the sympy script. Key contrast:

- sympy line 120: single shot `Sigma = C.T * Mint.inv() * C` — direct 4x4 matrix inverse and quadratic form.
- mathematica lines 87-113: sequential elimination. Step 1 (line 90): phi contributes `sigmaPhi = (lambda_B^2/A_phi) * outer(v, v)` directly. Step 2 (lines 97-101): build `uMassInv` as a Sherman-Morrison-like sum `Inverse[diag(A_U, A_U)] + (lambda_R^2/(A_U * Delta_UW)) * outer(v, v)`. Step 3 (lines 106-112): contract `C_U^T . uMassInv . C_U` for the isotropic Xi piece and add a hand-derived `sigmaW = ((A_U lambda_W^2 + 2 lambda_R lambda_U lambda_W) / Delta_UW) * outer(v, v)` for the wall-q -> W -> wall-q closed loop.

These are demonstrably different computational paths to the same `Sigma_wall`: sympy uses the canonical Schur formula, Mathematica uses a partitioned/Sherman-Morrison inversion. They produce identical results (verified by both engines' explicit assertions: sympy line 127 and mathematica line 116 both null `Sigma - (Xi I + alpha vv^T)`). Genuine cross-engine check, not a port.

A similar independence holds for `kappa_sel^2`: sympy uses Hellmann-Feynman exclusively; Mathematica uses Hellmann-Feynman AND direct `Eigensystem` projection AND verifies the two agree (M8). The selected-eigenvalue identification in Mathematica (lines 205-209) uses an explicit `Position` match against the expected `tr/2 - Sqrt[disc]/2` form, so the branch choice is forced from outside and cannot silently swap.

No `mathematica_transliteration` finding.

## Engine cross-check

Both engines produce the same closed forms for every shared assertion:

| Object | sympy output | mathematica output |
|---|---|---|
| sigma | `88/(9*pi**2)` (line 345) | `88/(9 Pi^2)` (line 11) |
| Xi_0 | `lambda_U**2/Omega_U**2` (line 356) | `lambdaU^2/OmegaU^2` (line 17) |
| Delta_0 | `Omega_U**2*Omega_W**2 - 88*lambda_R**2/(9*pi**2)` (line 357) | `OmegaU^2*OmegaW^2 - (88*lambdaR^2)/(9*Pi^2)` (line 18) |
| alpha_0 | combined rational (line 358) | combined rational (line 19) — algebraically equal after common-denominator rearrangement |
| beta_5 | `81*pi**4*Gamma_port*(Omega_U^2 lambda_W + lambda_R lambda_U)^2/(9 pi^2 Omega_U^2 Omega_W^2 - 88 lambda_R^2)^2` (line 821) | `(GammaPort*(lambdaR lambdaU + lambdaW OmegaU^2)^2) / (OmegaU^2 OmegaW^2 - 88 lambdaR^2/(9 Pi^2))^2` (line 37) — equal after multiplying num/denom by (9 Pi^2)^2 |
| kappa_sel^2 (HF) | sympy line 838ff prints the full radical form | mathematica line 43: `(88 + 8(968 alphaLoad + 63 DeltaKax Pi^2)/Sqrt[...])/(18 Pi^2)` — equal after factoring |

All cross-engine `PASS`. `engines_agree: true`.

## Verdict justification

The math holds. Both engines independently verify the Schur split (D1, the bottom-line theorem of the stage) via genuinely different routes, and both verify the `beta(omega)` closed form (D3) and the first-order `Pi`-expansion. The static branch data (D2) is computed by direct substitution from D1 and printed — algebraically implied, but not anchored by an explicit assertion against the literal closed forms. The selected odd coefficient (D4) is computed in pieces but never asserted as a single equality against the paper's eq `selected-odd`. Sympy's `kappa_sel^2` verification is asymmetrically weaker than Mathematica's because sympy lacks the direct eigenvector projection cross-check. One extra script-side check (`alpha_crit`) is verified in both engines but does not appear in the paper card.

Attacks attempted that did not break the math:
- Tried to find a sign mismatch between sympy's `eta = kappa_0 * kappa_1 = -8 sqrt(2)/(9 pi^2) < 0` and the stationarity equation: sympy uses `eta` consistently and the sign is correct (kappa_1 < 0 by construction at line 76, kappa_0 > 0).
- Tried to find a wrong-branch `kappa_sel^2` in sympy via the limit checks alone: limits pin down `kappa_0^2` at alpha=0 and `sigma` at alpha->infty, but the interior is what genuinely forbids the wrong branch - hence F3.
- Checked whether the `Pi` symbol clash with `sp.pi` could cause cross-contamination: at sympy line 72 `Pi = sp.symbols("Pi")` declares the symbol explicitly; `sp.pi` is the constant; no collision in `sp.simplify`.
- Verified `expect_zero` does catch nonzero matrix entries (line 51 raises if `any(entry != 0)`). Not tautological.
- Checked whether the printed `alpha_0` in the sympy output (line 358) and mathematica (line 19) agree: after putting both over a common denominator, identical.
- Checked the Mathematica `sigmaW` hand-derivation against the W-only contribution one expects from inverting the inner 3x3 (U_0, U_1, W) block: the numerator `aU * lambdaW^2 + 2 lambdaR lambdaU lambdaW` is the cross-channel coupling-product that comes from `(C_W^T + C_R-inversion-correction)`-style expansion; matches direct computation against the Mathematica inversion path (verified end-to-end by M1 anyway).

Verdict: `findings`. Two medium-severity `insufficient_verification` findings (F2, F3) that add asserted anchors for D4 and the `kappa_sel^2` cross-check; one low `paper_misalignment` for the stage-number docstring/banner drift (F1); and one low `paper_misalignment` for the `alpha_crit`-not-in-card item that needs user resolution (F4). `stop_cold: null`. F1, F2, F3 are safely applicable by Codex; F4 routes to the user.

## Self-test notes

I checked: (1) **Variable independence** — the proposed F2/F3 `expect_zero` arguments use already-defined script variables (`beta5`, `kappa_sel_template`, `Omega_U`, `lambda_W`, etc.), all carrying nontrivial dependence on the substituted symbols. The differences are not identically zero by construction; they collapse to zero by genuine algebraic equality (e.g. `beta_5 == Gamma_port * (Omega_U^2 lambda_W + lambda_R lambda_U)^2 / Delta_0^2` is the same identity already asserted in Mathematica's M7, so sympy will pass the new F2 binding). (2) **Symmetry/parity** — not applicable; the new checks are algebraic, not integrals. (3) **Trivial-case pre-check** — at alpha=0, the F3 direct-eigenvector route gives `e_- = (1, 0)`, so `(v.e_-)^2 = kappa_0^2`, matching A8; the closed-form HF expression also gives `kappa_0^2` at alpha=0; difference vanishes. At alpha->infty, `e_- -> v/||v||`, giving `(v.e_-)^2 = ||v||^2 = sigma`, matching A9. The interior is forced by the explicit equality, not just the two limits. (4) **Path specifications** — F1/F2/F3 modify existing files (no missing-script finding); paths are `scripts/...py` and `mathematica/...wl`. (5) **Paper round-trip** — F2 binds the script to the paper's eq `selected-odd` *as written*, with the explicit `(Omega_U^2 lambda_W + lambda_R lambda_U)^2/Delta_0^2` denominator; no new paper claims introduced.
