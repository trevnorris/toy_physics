# Stage Verification Coverage Baseline

This document is the stage-coverage control sheet for the PDE ledger archive.

Use it together with `STAGE_PROVENANCE_INDEX.md`:

- `STAGE_PROVENANCE_INDEX.md` holds the exact raw note and audit artifact paths.
- `STAGE_VERIFICATION_COVERAGE.md` summarizes the current verification surface,
  exposes the main gaps, and gives us a stable baseline for audit planning.

Snapshot date: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 6** — re-verification of four IV.5-range stages whose first-pass fixes had been applied orchestrator-direct (Codex bypassed) and have now been Codex-reconciled + re-verified. Stages **143** (exp-remainder positivity), **144** (threshold-Pi_* independence), **146** (affine-law residual), **147** (rigidity-kernel projection) all REMEDIATION-`verified` (2026-05-29); both engines exit 0/0 on every stage; this status **supersedes the earlier orchestrator-direct "verified"** these four carried from their first-pass IV.5 batch. No new stages added to the verified count (all four were already counted from IV.5); this is a verification-hardening / integrity re-pass. **`material_change: false` on every stage**. **NO checkpoints in the batch-6 range** (IV.5 has no checkpoints; nearest checkpoints are 096 and 105). Findings closed per stage: **143** (2 findings, insufficient_verification) — replaced a cubic-Taylor-coefficient-only "positivity" gate (which passed for a wrong remainder like `Pi^3/6 − Pi^4`) with a GENUINE GLOBAL positivity proof of `exp(Pi)−1−Pi−Pi^2/2 > 0` on `Pi>0` (SymPy Taylor-remainder chain `R(0)=R'(0)=R''(0)=0, R'''(Pi)=exp(Pi)>0`; Mathematica primary `Reduce[...>0,piM,Reals]`→True with the same Taylor backing). **144** (1 finding, transliteration) — replaced a line-by-line SymPy-mirror `Pi_*` block with an INDEPENDENT cleared-denominator bracketed `FindRoot` on `gThresholdResidual[p]=2p(2p e^p+pi)−gMinus(4p^2+pi^2)(e^p−1)` (sign kept as `(e^p−1)`, NOT the §131-F3 `(1−e^p)` trap) + residual-near-zero witness, anchored to the OWNING stage-131 `Pi_*` value `1.50882951349315558300555075595`; MIRROR-policy `needs_user_resolution` cleared (consult Q2: a load-bearing transcendental root is NOT a sanctioned mirror → independent route, precedent §131-F3/§142). **146** (3 findings) — R1 SymPy `nsolve(prec=50)` + tightened affine-residual assert from `1e-15` to raw `<1e-25` (~1e-51); R2 Mathematica removed `Chop[...,10^-6]` and asserts raw `<10^-25` (~1e-58/exact 0) via the directive-authorized SYMBOLIC ENDPOINT-INTEGRAL fallback (sanctioned, genuinely hits `<1e-25`); R3 banner `FINITE-CORRECTION EXPANSION`→`FIRST-ORDER EXPANSION` (`.wl` script fix matching the paper-card title). **CONSULT caveat (146):** the affine residual collapses to `(1−eps)(gPi(Pi_*)−gMinus)` and `gMinus` also feeds `Pi_star`, so this check is NOT an independent guard against `gMinus` (documented scope limit — tests intercept-vs-direct-integral + kernels/source). **147** (5 findings, the heavy one) — R1/R5a CAS-autodiff of `T_m(Pi)=sqrt((9/20)Pi/(1−Sformula/4))` replacing a hand-typed `A_T` chain rule; R2/R3 projection identity `∫W_*(σ−Σ_*)=A_T(ḡ_σ−g_*)+B_T(S̄_σ−S_*)` [σ=2x] + full-symbolic x-independence check + (consult-Q5-ADDED) source-centering assertion `∫Σ_*W_*=0` (the projection identity alone is blind to the centering constants; this orthogonality check tests them); R4/R5b source-moment quadrature integrals vs the closed forms replacing `g_*`/`S_*` resubstitution; R5 independent Wolfram `D[]`/`NIntegrate` (not a SymPy port); R6 (consult Q6) the `|A_T|/B_T` ratio `31.6785` was MISLABELED paper-quoted → relabeled as the script's own computed cross-check (the genuine paper literals `A_T=−4.27263956256927`, `B_T=0.134875005736706` are sourced from `paper/appendices/stage_appendix_part04.tex:846,848`). Directions settled by ONE Claude+Codex read-only consult (4 CONCUR / 2 DISPUTE-resolved, none conceptual, none escalated), recorded at `redteam/codex_reviews/_consult_batch6.md`; NO new paper-cleanup items, logged in `PAPER_CLEANUP_TRACKER.md` P4-47. See `redteam/verifications/stage_{143,144,146,147}.md`. Previous batch-5 snapshot below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 5** — re-verification of four Family-1 mouth-gain/branch stages whose first-pass fixes had been applied orchestrator-direct (Codex bypassed) and have now been Codex-reconciled + re-verified. Stages **134, 137** (IV.4-range canonical_mouth_gain / schur_reduction) and **139, 142** (IV.5-range family1_actual_mouth_gains / mouth_gain_susceptibility) all REMEDIATION-`verified` (2026-05-29); both engines exit 0/0 on every stage; this status **supersedes the earlier orchestrator-direct "verified"** these four carried from their first-pass IV.4/IV.5 batches. No new stages added to the verified count (all four were already counted from IV.4/IV.5); this is a verification-hardening / integrity re-pass. **`material_change: false` on every stage**. **NO checkpoints in the batch-5 range** (`is_checkpoint: false` confirmed in MANIFEST for 134/137/139/142; IV.4/IV.5 have no checkpoints). New check types added per stage: **134** (canonical_mouth_gain) R1 tautological_check — REMOVED the X−X "canonical gain line" SymPy assert (`intercept = Pi_*` vs the literal that defined `Pi_*`), replaced with a provenance comment deferring outlet/gain-line consistency to Stage 135 (the paper-card downgrade itself was done in IV.4); kept the corroborated shell-limit + 3 S_q spot-checks + S_q(Π_*) checks; `.wl` gain-line block relabeled "printed only; not asserted". **137** (schur_reduction) R1/R2/R3 — added an INDEPENDENT **matrix-Schur reconstruction** of `rho_c, sigma_c` from the physical core matrix `M_core=[[K_s,λ],[λ,−K_q·D]]`, `v=(g_s,g_q)`, `δ=vᵀM_core⁻¹v` (SymPy `apart`+`limit`; Mathematica `Inverse[mCore]`+`Normal[Series]`, mirroring the verified OWNER stage 114); replaced the X−X static-limit check with a matrix-route-vs-reduced-envelope identity, and the `S_q=0` outlet check (which deleted M_q) with a nonzero-S_q mixed-channel check `M_q·S_q == −L·σ_c_schur·S_q/Θ` + sign non-vacuity guard — **9 PASS / 0 FAIL**. **139** (family1_actual_mouth_gains) R1/R2/R3/R4 — R1 replaced tautological outlet asserts with an independent `S_q(Π_*)` recompute from the Stage-134 closed-form kernel; R2 relabeled `R_q^comp=1/4` as definitional-consistency + added a **branch-discrimination value anchor** `g_-^F1=0.758035078944662826919680890414` (discriminates the lower from the upper branch), made the Mathematica `gMinus` the DIRECT closed form (sanctioned mirror); R3 (paper_misalignment) resolved no-script-fix (susceptibility closure is a Stage-140 deliverable; existing P4-42); R4 (provenance) resolved scripts-already-correct. **142** (mouth_gain_susceptibility) R1/R2 — R1 relabeled the self-solved `R_q(Π_*)=1/4` as solver-consistency + added a non-tautological **cross-stage independent-Π_* anchor** evaluating `R_q` at STAGE 131's independently-derived `Π_*=1.50882951349315558300555075595` (cleared-denominator FindRoot, batch-4-verified — NOT 142's own nsolve); R2 replaced a `Normal[Series]` self-comparison with an INDEPENDENT **projection integral** `∫₀¹ σ_Π(z)·cos(πz/2) dz` (σ_Π = Stage-129 source law) closing to the hardcoded `gPi` with symbolic-zero residual; 5 external-decimal anchors retained; tolerance 1e-15 kept — **11 PASS / 0 FAIL**. Directions settled by ONE Claude+Codex read-only consult (7/7 CONCUR, none conceptual, none escalated), recorded at `redteam/codex_reviews/_consult_batch5.md`; NO new paper-cleanup items, logged in `PAPER_CLEANUP_TRACKER.md` P4-46. See `redteam/verifications/stage_{134,137,139,142}.md`. Previous batch-4 snapshot below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 4** — re-verification of four stages whose first-pass fixes had been applied orchestrator-direct (Codex bypassed) and have now been Codex-reconciled + re-verified. Stages **125, 126** (IV.3-range) and **130, 131** (IV.4-range) all REMEDIATION-`verified` (2026-05-29); both engines exit 0/0 on every stage; this status **supersedes the earlier orchestrator-direct "verified"** these four carried from their first-pass IV.3/IV.4 batches. No new stages added to the verified count (all four were already counted from IV.3/IV.4); this is a verification-hardening / integrity re-pass. **`material_change: false` on every stage**. Fixed via a 2-parallel Codex run enabled by a user-approved `redteam.sh` MANIFEST flock. **NO checkpoints in the batch-4 range** (IV.3/IV.4 have no checkpoints). Findings closed + PASS-line counts per stage: **125** (positive_source_theorem) F1 insufficient_verification — strengthened the peaked-source check from a weak smallness test `abs(g_a(100))<0.05` (which accepted small NEGATIVE moments) to GENUINE paper-bound range checks `0 ≤ g ≤ 1` (the `≥0` line, NOT wrapped in abs, fails on a sign error in the closed-form moment), retaining `<1/20` as a subordinate trend witness; Mathematica mirrored with two range bounds at the a=100 proxy (deviation: extra precision + `Re[]`/`Chop` on the incomplete-gamma closed form, benign) — **SymPy 3 new check lines / Mathematica 15 PASS**. **126** (positive_source_families) F1 insufficient_verification — replaced endpoint/corner-only positivity SAMPLING with a GLOBAL box-positivity proof (exact affine-in-ξ check `d²σ_ξ/dξ²==0` PLUS both ξ-endpoint slices nonneg over `[0,L]` ⟹ nonneg on the (z,ξ) box; Mathematica uses `Resolve[ForAll[...]]` with the directive-authorized `lM→1` scale-covariant quantifier; an interior-curved perturbation fails the affineness check, non-tautological) — **SymPy +3 lines / Mathematica 7 PASS**. **130** (mouth_bias_map) F1 insufficient_verification — replaced the finite 6-point monotonicity sweep with the notes-§2 FKG/Chebyshev symmetrized-covariance GLOBAL certificate proving `dg_Π/dΠ>0 ∀Π>0` and Π_* uniqueness (density normalization, symmetrized double-integral identity, integrand sign `f'(z)=−(π/2L)sin(πz/2L)<0` on (0,L), consistency `dg/dΠ=−Cov/L`, uniqueness bracket `2/π<g_-<1`; the directive-sanctioned optional SymPy `reduce_inequalities` probe was removed as non-load-bearing, the load-bearing sign arg intact; the double integral simplified under the 600s cap, no fallback) + F2 (boxed `g_Π=2Π(2Πe^Π+π)/((4Π²+π²)(e^Π−1))` closed form) pre-existing reconciled-PASSED, unchanged — **Mathematica 11 PASS / SymPy print+assert guards**. **131** (parent_mouth_threshold) F1 tautological_check — the parent-threshold-identity-at-Π_* check was X−X (purely DEFINITIONAL); per consult option (ii) the Anchor-3 block was DROPPED entirely in both engines + F2 insufficient_verification — replaced weak `g(2Π_*)≠g_-` with 3 branch-discrimination checks (lower-branch membership, singular-branch exclusion `g_nat−g(Π_*)≈Δg_-=0.241964921055337`, upper-branch exclusion `|g(Π_*)−g_+^{F1}|>1` with `g_+^{F1}≈2.79795199200529`) + F3 transliteration — **NOT accepted as a policy-mirror**, added an INDEPENDENT Mathematica route for Π_* via the cleared-denominator residual solved by bracketing `FindRoot[{piM,1.4,1.6}]` (per the IV.5 139/143/144 transcendental-numerical-root precedent); reconciled-PASSED Anchor1 (Π_*=1.50882951349316), Anchor2 (g'(Π_*)=0.0714453558083195), F4 (g_-^{F1} closed form vs 0.758035078944663), banner STAGE 131 unchanged — **SymPy 6 / Mathematica 6 PASS**. The resolution directions came from a Claude+Codex consult (Codex base `019e7594`, all CONCUR, none conceptual), recorded at `redteam/codex_reviews/_consult_batch4.md`; paper-card cross-reference / informational-only status logged in `PAPER_CLEANUP_TRACKER.md` P4-45 (NO new paper-cleanup action items — note specifically that stage 130's paper-card claim "monotonicity, and unique Family-1 point Π_*" is now PROPERLY verified globally). See `redteam/verifications/stage_{125,126,130,131}.md`. Previous batch-3 snapshot below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 3** — re-verification of four IV.3-range stages whose first-pass IV.3 fixes had been **tainted-applied** (edits committed in the IV.3 pass `b4e02d8` while Codex was bypassed) and have now been Codex-reconciled + re-verified. Stages **117, 118, 119, 122** all REMEDIATION-`verified` (2026-05-29); both engines exit 0/0 on every scripted stage; this status **supersedes the earlier tainted "verified"** these four carried from their first-pass IV.3 batch. No new stages added to the verified count (all four were already counted from IV.3); this is a verification-hardening / integrity re-pass. **`material_change: false` on every stage**. **NO checkpoints in the batch-3 range** (IV.3 has no checkpoints). Findings closed + PASS-line counts per stage: **119** (parent_balance) F1 tautological_check de-tautologized (tube-length law `rc → rhat²` link ties L_W to the family `rhat`) + F2 insufficient_verification (`T_m (±)` branch match vs notes §5 closed form; `.wl` uses a `stripCE` ConditionalExpression-stripper) — **8 checks/PASS per engine**, genuine independence confirmed (not X−X). **118** (parent_core) F1 paper_misalignment (λ SIGN) → resolved **direction (a) = MINUS** (consistent with the script's own independently-derived section IV `sq_coeff = −q_* ϱ_s v0 A_q`, the notes' three boxed minus signs, and downstream stage 123's UN-squared λ; NOT conceptual — paper card states no λ sign, so no paper-cleanup item) + F2 insufficient_verification added three asserts (`K_q closed form`, `g_s closed form`, `lambda from bilinear`) — **12 sympy / 14 mathematica PASS**. **CAVEAT (118):** verify agent flagged `K_q closed form` as effectively X−X (target restates the line-82 definition) — **non-blocking** because K_q's value is independently anchored upstream (stage-118 sympy line 50); `g_s closed form` genuine; `lambda from bilinear` modest (sign/factorization consistency). **117** (outlet_core_status, a status-consolidation card) F1 mathematica_transliteration = **accepted policy-mirror** (pure rational series-coefficient classification of which outlet-deformation classes preserve Ŷ₂; no independent-route rewrite forced on a consolidation card) + F2 tautological_check de-tautologized by importing the stage-116 forward tube-length closed form `L_W = pi a sqrt((1+r_c)/3)/2` to build `kappa0 = 4 L_W²/(pi² a²)` (the target-inverting dead `Lw_required` was DELETED) + F3 tautological_check = comment-only provenance fix (γ₀=(1+r_c)/9 corrected from false "Stage 119"/"derived upstream" to a postulated pure-scale ANSATZ; κ₀'s "derived upstream → stage 116" retained; `gamma_c = 1/9` assertion byte-unchanged) + F4 tautological_check reconciles already-wired capstone classification flags (one nontrivial survivor `compensated Robin-mixed core realization`) — **12 PASS per engine**. **CAVEAT (117 F2 honest strength):** F2 is now a falsifiable consistency check against the in-script forward L_W law (a wrong coefficient breaks the O(z⁶) core-collapse residual), NOT an in-stage re-derivation of the D/N BVP eigenvalue (that lives at stage 116). **122** (mouth_source_compensation_test, **SymPy-ONLY, no `.wl`**) F1 tautological_check de-tautologized — bare reciprocal traction ratios (`T_ratio = 1/g±`, never-fail `g·(1/g)−1`) replaced by ratios DERIVED from the stage-119 proportionality law `g = C/T_m` (C a branch-independent constant that cancels symbolically) + the `g_nat = 1` natural-branch ansatz, compared against the notes §5 closed form `1/g±` — **7 `= 0` lines**, confirms the natural-defect uses `100π²` (NOT `168π²`). **CAVEAT (122 moderate strength):** a consistency check of `g_nat=1` + C-branch-independence; would fail if `g_nat≠1` or the normalization were branch-dependent. The 117/118 resolution directions came from a Claude+Codex consult (Codex base `019e74f7`, all CONCUR), recorded at `redteam/codex_reviews/_consult_batch3.md`; paper-card cross-reference / informational-only items logged in `PAPER_CLEANUP_TRACKER.md` P4-44 (NO new paper-cleanup action items). See `redteam/verifications/stage_{117,118,119,122}.md`. Previous batch-2 snapshot below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 2** — re-verification of four more stages whose first-pass fixes had been applied orchestrator-direct while Codex was bypassed. Stages **105, 106, 109, 112** all REMEDIATION-`verified` (Codex applied + clean Claude verify agents, 2026-05-29); both engines exit 0/0 on every stage; this status **supersedes the earlier tainted orchestrator-direct "verified"** these four carried from their first-pass IV.2 batch. No new stages added to the verified count (all four were already counted from IV.2); this is a verification-hardening / integrity re-pass. **`material_change: false` on every stage** (no derived constant moved). Findings closed per stage (from verification frontmatter): **105** (checkpoint) F1 mathematica_transliteration resolved — `.wl` chi_Q derivation rewritten along an INDEPENDENT residue/`Reduce`-witness path (single-ratio retarded module, `Im[]` projection, operator-product polynomial inversion; non-transliterated names) + F2 paper_misalignment (stale "Stage 88" labels → canonical "Stage 105"), 9 Mathematica PASS; **106** F1 paper_misalignment (doc-only citation fix: item (ii)→Stage 102, item (iii) DtN fingerprint→Stage 104, chi_Q=1→Stage 105) + F2 tautological_check resolved (by-construction K0/K2/K4 relation → INDEPENDENT hardcoded `*_target` literal consistency `K0_target·K4_target = 4·K2_target²`) + F3 mathematica_transliteration resolved (`.wl` re-authored on independent `Yret` one-pole omega-series path deriving chi_Q=1 from ω⁵ match) + F4 insufficient_verification (Delta_Q first-order sensitivity slope = −Gamma5_target on both engines), 11 Mathematica PASS; **109** F1/F2 tautological_check resolved ("preservation substitution" de-tautologized on BOTH engines: substitutes the INDEPENDENT closed form `a_5 = -5b/9 - a0/27`; `.wl` also re-derived via separate num/den linearization + denominator series-inversion, `coeff=(chiSeries-1)/eps` intermediate removed) + F3 paper_misalignment (NO script change), 4 Mathematica PASS; **112** F1 paper_misalignment (stale "Stage 95" → "Stage 112") + F2 mathematica_transliteration resolved (independent Stage-92 linearized (b,a0,a5)=(0,3σ,−σγ) cross-check deriving gamma_W=1/9, a route the `.py` lacks) + F3 symbol_assumption_error resolved (`sigma_W != 0` qualifier folded into BOTH engines: preservation `chi_B=1 ⟺ sigma_W(1−9 gamma_W)=0`, nontrivial branch sigma≠0 forces gamma_W=1/9 via `Reduce[...&&sigma!=0]`, degenerate sigma=0 gives chi_B=1 for any gamma), 14 Mathematica PASS. **Paper_misalignment resolutions** (Claude+Codex consult `019e748e`, recorded at `redteam/codex_reviews/_consult_batch2.md`) all non-conceptual: 105 F2 & 112 F1 = canonical-internal-stage-number labels (per IV.4/IV.5 convention); 109 F3 & 106 F1 = paper-card cross-references (no script change), logged in `PAPER_CLEANUP_TRACKER.md` P4-43. See `redteam/verifications/stage_{105,106,109,112}.md`. Previous batch-1 snapshot below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 1** — re-verification of the four stages whose first-pass fixes had been applied orchestrator-direct while Codex was bypassed. Stages **108, 116, 151, 170** all `verified` (116 & 151 in prior commit `e1cdfec`; 170 & 108 this session); both engines exit 0/0 on every stage. Findings closed per stage (from verification frontmatter): **116** (batch IV.3) verdict `verified`, sympy/math 0/0, F1 mathematica_transliteration + F2 insufficient_verification resolved (independent eigenvalue solve of `Cos[u]==0` on `(0,π)`; tautological renorm block → labeled prints); **151** (batch IV.6) verdict `verified`, sympy/math 0/0, F1 insufficient_verification resolved — **methodology: Mathematica = full all-`Pi_star` symbolic authority; SymPy = an EXACT 5-point cross-check at rational `Pi_star ∈ {1/2,1,3/2,2,5/3}` symbolic in `r1,r2,A_T,B_T,gprime`** (supersedes the IV.6-snapshot mpmath description; `material_change: false`); **170** (batch V.1) verdict `verified`, sympy/math 0/0, F1 tautological_check resolved (Section-5 lanes rerouted through the derived Section-2 maps `dkappa_from_du2`/`dgamma_from_dP0`); **108** (batch IV.2) verdict `verified`, sympy/math 0/0, F2 mathematica_transliteration + F3/F4 tautological_check resolved (independent non-`Series` `chiGenAlt` route; Class A anchored to fingerprint `1+z²/9+4z⁴/81+iz⁵/27`; Class C anchored to `-Σ0/27`). **108 F1 = paper-card cleanup, no script change** (direction (a), Claude+Codex concur — see `redteam/codex_reviews/_consult_108_f1.md`; deferred to `PAPER_CLEANUP_TRACKER.md` P4-42). 108 carries `material_change: true` (verification surface strengthened; no derived value changed; no downstream propagation); 116/151/170 `material_change: false`. Deferred paper/doc cleanups (106/139/148 + 112 precision caveat) logged in `PAPER_CLEANUP_TRACKER.md` P4-42. No new stages added to the verified count (all four were already counted from their first-pass batches); this is a verification-hardening / integrity re-pass. Previous V.1 snapshot below.)

Snapshot date prior: `2026-05-28` (batch V.1 close — first-pass paper-grounded audit under v2 prompt for stages 164-175, no checkpoints, no status-only units. 22 findings closed across all 12 (dirty) stages, 0 blocked. 1 material_change (stage 170 — additive new Sec. 5 weak-axisymmetric (1,1/2,-1) signature coverage; no carried/derived result changed, zero downstream propagation). **Cluster A only**: all 12 stages carried a stale `STAGE N-17` script banner (164->147 ... 175->158); 11 fixed inside their per-finding directives, 4 residuals (164.py, 165.py+.wl, 174.py) mass-fixed in place. **No Cluster B** (no body-text citation re-attribution) and **no Cluster C** (no paper-card Checks downgrade). One **non-cluster paper_misalignment** (170 F1, weak-axisymmetric (1,1/2,-1) signature unverified) resolved by user direction (a) = add the missing check to both engines — first non-cluster paper_misalignment in the v2 run. `mathematica_transliteration` on 9/12 (highest v2 share); 7 independent routes added (164/166/170/171/172/173/174), 2 accepted as policy mirrors (169 F2, 175 F3-step3). Three orchestrator catches: (1) 166 matrix round-trip vector residual -> `Total[(...)^2]` scalarization; (2) 175 F1 cross-check reduced to a simplify-commutes identity -> kept only the load-bearing `Sigma_N-dln(Lambda^2/K)`; (3) 171 F1 directive route was tautological (`zCombFormula` rebuilt from the same `D[...]` as `zCombExact`) -> reworked to collected-literal target + independent `Series` route, caught by the first verification (`needs_rework`) and re-verified clean. Thirteenth consecutive batch clear of stop-cold.)

Previous snapshot: `2026-05-27` (batch IV.5 close — first-pass paper-grounded audit under v2 prompt for stages 139-150, no checkpoints, zero material_change. Same renumbering pattern as IV.4 — 11 `.wl` banners + 6 `.py` banners off by -17 (Cluster A) and 4 notes H1 lines off by +102; plus 22 body-text forward-stage citations off by -51 or -102 re-attributed to current upstream (Cluster B). Stage 144 paper card Checks downgraded to carry-forward citations of 135/140 via Cluster C, mirroring IV.4's stage 134. Three status-only stages clean (141, 145, 149). Eleventh consecutive zero-redirection batch. Four orchestrator catches in the rework loop including: stage 148 directive's prescribed closed form was wrong (`168π²` typo from stage 148 notes; correct `100π²` confirmed at stage 126 upstream), stage 139 pitfall #13 recurrence on the `Pi_*)` comment substring, plus three SymPy tolerance/precision/Integrate-fallback adjustments.)

## Scope

- Canonical stage range: `001--253`
- Canonical stage source model:
  - Parts `I--VII`: compact stage cards
  - Part `VIII`: full derivation-stage files

## Coverage Totals

| Metric | Count |
|---|---:|
| Total stages in archive | 253 |
| SymPy audits present | 240 |
| Mathematica audits present | 165 |
| Numerical stress artifacts present | 15 |
| Stage-specific review notes present | 178 |
| Stages with no SymPy audit | 13 |
| Stages with no executable audit | 11 |
| Stages with Mathematica but no SymPy | 2 |

Mathematica counts above are presence counts only. They are **not**
independence counts. See `MATHEMATICA_MIRROR_POLICY.md` for the current rule
that distinguishes secondary replay coverage from genuinely independent
Mathematica mirrors.

The SymPy runner summary reports `241` passing files because it includes the
repo-level `moving_throat_pde_master_sympy_audit.py` in addition to the `240`
stage-level SymPy audits counted above.

### Red-team verification status

Distinct from "audit file present" (counted above) is "red-team verified" —
the stage has passed the `audit → directive → codex → verifier` pipeline run
out of `redteam/`, with both engines independently checking load-bearing
claims and a clean-context verifier agent confirming the directive's intent
was honored. See `redteam/BATCHES.md` for the live batch table.

As of `2026-05-28`: **187** of 253 stages red-team verified (73.9%). With V.1 closed (first audit pass for stages 164-175 under v2), the entire range 001-175 is now paper-aligned at v2 depth. No checkpoints, no status-only units in V.1. **One material_change** (stage 170 — additive Sec. 5 weak-axisymmetric (1,1/2,-1) signature coverage added per user direction (a); no carried/derived result changed, zero downstream propagation). **First non-cluster paper_misalignment in the v2 run** (170 F1), resolved by adding the missing script check rather than a carry-forward citation — content-disambiguation confirmed the deliverable is stage 170's own (notes Sec. 5), with 171/173 verifying the lane pattern only on other quantities. `mathematica_transliteration` hit 9/12 (highest v2 share); 7 independent routes added, 2 policy mirrors (169, 175). No Cluster B or Cluster C this batch.

As of `2026-05-28`: **175** of 253 stages red-team verified (69.2%). With IV.6 closed (first audit pass for these 13 stages under v2), the entire range 001–163 is now paper-aligned at v2 depth. No checkpoints in IV.6; three stages clean (153 status-only, 159 and 162 substantive-clean). **One material_change** in IV.6 (stage 151 SymPy rewrite from symbolic-integration to mpmath numerical because `sp.integrate` of `Pi_star * exp(-Pi_star*x) * cos(pi*x/2)` over [0,1] with free `Pi_star` hung >30 min CPU); downstream impact zero since the script verifies algebraic identities not numeric carry-forwards. **First forward (downstream) carry-forward in v2 batches** via Cluster C on stage 158 (items 2-3 cite `\ref{stage:159}` / `\ref{stage:162}` / `\ref{stage:163}`). The −51 (IV.4 refs) and −102 (IV.5 refs) body-citation offsets carried over from IV.4/IV.5; new in IV.6 is the **−85 offset** for IV.6 self-references (239-248 → 154-163), disambiguated from −102 by content cross-check (e.g., stage 158's "Stage 241, exact co-evolving compensated Family-1 point" matches stage 156 renormalized canonical branch, not stage 139 family1_actual_mouth_gains).

As of `2026-05-27`: **162** of 253 stages red-team verified. With IV.5 closed (first audit pass for these 12 stages under v2), the entire range 001–150 is now paper-aligned at v2 depth. No checkpoints in IV.5; three stages clean (141, 145, 149 — all status-only). **Zero material_change** in IV.5 — every fix was renumbering, notes re-attribution, paper-card downgrade, or script-substance addition; no derived numerical constants moved. Pitfall #13 (Mathematica `_*)` substring in comment closing prematurely) **re-confirmed** at stage 139 — workaround applied (ASCII-safe names `piStar`, `sQStar`). The two systematic stage-renumbering shifts and the body-text forward-citation pattern that surfaced in IV.4 repeated in IV.5 with the same offsets (banners −17, notes H1 +102, body citations −51 or −102).

| Batch | Range | Stages | Verified | Date |
|---|---|---:|---:|---|
| I.1 | `001--012` | 12 | 12 | 2026-05-21 (v1) / 2026-05-25 (v2 paper-grounded) |
| I.2 | `013--023` | 11 | 11 | 2026-05-21 (v1) / 2026-05-25 (v2 paper-grounded) |
| II.1 | `024--036` | 13 | 13 | 2026-05-22 (v1) / 2026-05-26 (v2 paper-grounded) |
| III.1 | `037--048` | 12 | 12 | 2026-05-22 (v1) / 2026-05-26 (v2 paper-grounded) |
| III.2 | `049--060` | 12 | 12 | 2026-05-22 (v1) / 2026-05-26 (v2 paper-grounded) |
| III.3 | `061--072` | 12 | 12 | 2026-05-22 (v1) / 2026-05-26 (v2 paper-grounded) |
| III.4 | `073--084` | 12 | 12 | 2026-05-25 (v1) / 2026-05-27 (v2 paper-grounded) |
| III.5 | `085--090` | 6 | 6 | 2026-05-27 (v2 paper-grounded, first-pass) |
| IV.1 | `091--102` | 12 | 12 | 2026-05-27 (v2 paper-grounded — first pass) |
| IV.2 | `103--114` | 12 | 12 | 2026-05-27 (v2 paper-grounded — first pass) |
| IV.3 | `115--126` | 12 | 12 | 2026-05-27 (v2 paper-grounded — first pass) |
| IV.4 | `127--138` | 12 | 12 | 2026-05-27 (v2 paper-grounded — first pass) |
| IV.5 | `139--150` | 12 | 12 | 2026-05-27 (v2 paper-grounded — first pass) |
| IV.6 | `151--163` | 13 | 13 | 2026-05-28 (v2 paper-grounded — first pass) |
| V.1 | `164--175` | 12 | 12 | 2026-05-28 (v2 paper-grounded — first pass) |
| V.2 onward | `176--253` | 78 | 0 | pending |

Cumulative findings closed: ~451 (~219 v1 + 10 v2 from I.1 + 10 v2 from I.2 + 18 v2 from II.1 + 13 v2 from III.1 + 16 v2 from III.2 + 13 v2 from III.3 + 14 v2 from III.4 + 15 from III.5 + 27 from IV.1 + 16 from IV.2 + 27 from IV.3 + 22 from IV.4 (21 resolved + 1 blocked_legitimate on 137 F3 matrix-Schur) + 31 from IV.5 (30 resolved + 1 blocked_legitimate on 144 F4 transliteration policy) + **19 from IV.6 (19 resolved + 0 blocked, 1 material_change on 151 mpmath rewrite)** + **22 from V.1 (22 resolved + 0 blocked, 1 material_change on 170 — additive Sec.5 (1,1/2,-1) signature coverage)**, plus 1 blocked-legitimate from IV.1). Of the 27 IV.3 closes, 7 were `paper_misalignment` (Cluster A 4 notes-side numerical typos at 121/122/123/126 — `168π²→100π²` for three stages and `228→160` for stage 123; Cluster B λ sign flip at 118 — internal script inconsistency; Cluster C integral inequality coverage gap at 125; **Cluster D 117 transliteration + 2 tautological resolved via cite-upstream-and-downgrade**), 3 `mathematica_transliteration` (115, 117 F1, 125 F2 — independent-derivation insertions), 6 `tautological_check` (116 F1, 116 F2, 117 F2, 117 F3, 117 F4, 119 F1), 5 `insufficient_verification` (116 F4, 119 F2, 121 F2, 122 F2, 126 F3), 2 `hardcoded_result` (116 F3, included in F3 count), 1 `script_missing_paper_claim` (121 F3 Ω_W), 1 `paper_missing_script_claim` (126 F2 positivity), 1 `stale_output` banner (126 F4). Plus a 10-site banner-relabel sweep across all IV.3 scripted stages (115, 116, 117, 121, 122, 123, 125, 126 — orchestrator-direct, matching IV.2's pattern).
`tautological_check` dominant overall, `mathematica_transliteration` second.
`hardcoded_result` rose sharply in III.4 to 12 because the Family-1 numerology
cluster 075-084 packs many literal constants; III.5 quieted again (1 hardcoded
in 089 F4). v2 added `paper_misalignment` as the 10th category — **31** items
total across the eight v2 batches (7 in I.1, 3 in I.2, 3 in II.1, 3 in III.1,
2 in III.2, 4 in III.3 — 2 substantive + 2 banner relabels, 7 in III.4 — 4
substantive + 3 audit-flagged banner relabels, plus an 8-stage
orchestrator-direct banner-relabel sweep when the global-renumber leftover
turned out to be pervasive across III.4; **2 in III.5 — both substantive (087 F1
status/checkpoint consolidation, 089 F1 Pe_req=0 chain closure), plus a 12-script
orchestrator-direct banner-relabel sweep**);
zero user redirections in II.1, III.1, III.2, III.3, III.4, III.5, IV.1, IV.2, IV.3, IV.4, IV.5, IV.6, **V.1** (13 consecutive
batches — Codex was bypassed from III.5 onward per the III.4 availability lesson; orchestrator-direct
math-authority worked cleanly because the audit + grep evidence was conclusive. V.1 surfaced the first non-cluster
paper_misalignment (170 F1); the user took the recommended direction (a) — add the missing check — so no conclusion
was redirected and the streak holds).
v2 surfaces `insufficient_verification` prominently — 8 in II.1, 5 in III.1, 8 in III.2,
4 in III.3, 1 in III.4, 1 in III.5, 7 in IV.1, 2 in IV.2, 5 in IV.3, 4 in IV.4 (130 F1, 130 F2, 131 F3, 137 F1), 11 in IV.5 (139 F1, 140 F2, 142 F2, 143 F1, 144 F2, 144 F3, 146 F2, 147 F3, 148 F1) = **56** cumulative.
Stage 060 (v1 `material_change: true`) returned **clean (0 findings)** under III.2 v2.
**Stage 068 (v1 `material_change: true`) returned clean at v2**.
Stages now carrying `material_change: true`: 001, 004
(I.1 v2); 013, 014, 015, 018 (I.2 v2); 045 (III.1 v2 — structural-only, F_tr
export value unchanged); 060 (v1, clean at v2); 068 (v1, clean at v2). II.1,
III.1, III.2, III.3, III.4, III.5 v2 each added **zero** value-changing material_change. **IV.1 added one structural material_change at stage 100** (closure derivation strengthened from tautological cross-check to substantive `mhat_0^2 Gamma_5 = Gamma_5_target` imposition; no derived value changed; downstream stages > 100 not marked `upstream_stale`). **IV.2 added one structural material_change at stage 108** (Cluster B β-parameterized preservation submanifold added; the β=1 reduction already verified previously is unchanged, only the verification surface widened; downstream stages > 108 not marked `upstream_stale`). **IV.3 added one substantive material_change at stage 118** (λ sign flip from + to − in section V's `lam_uniform = qstar*v0*I_sq` closure; internal script inconsistency where section IV's bilinear derivation and the notes' boxed form both had `−` but section V dropped the minus during integration. Resolved by aligning section V with section IV. Downstream Schur reductions use `K_s K_q + λ²` and `(K_s g_q − λ g_s)²` — both sign-invariant under squaring — so no numerical downstream propagation; upstream_stale NOT flagged. Future-batch caveat: revisit if a downstream stage uses λ in an unsquared cross-term).
III.3 v2 introduced one orchestrator hot-fix on stage 064 Mathematica
(`Integrate[]` with symbolic functions does not factor constants — verify
integrands first; pitfall #9 candidate). III.4 v2 introduced one orchestrator
fix on stage 082 SymPy (`sp.nsolve` is unstable for `y tan y = 37` near
`pi/2` and jumps to far-away roots — use `mpmath.findroot(..., solver="bisect")`
instead; pitfall #10 candidate). **III.5 introduced two orchestrator hot-fixes:
(a) on stage 088 SymPy, `Y_rho.subs(omega**2/Omega_Q**2, u)` failed silently
because `sp.simplify` reshapes the denominator into `(Omega_Q**2 - omega**2)`
form and the combined ratio is no longer a syntactic subexpression — fix:
substitute `omega**2 -> u * Omega_Q**2` then `sp.simplify`. (b) On stage 088
Mathematica, a comment containing the substring `stage085_*)` was prematurely
closed by the embedded `*)`, causing `Syntax::sntx` and silently skipping the
F1 assertion and regime trichotomy while still reaching `Exit[0]` (rc=0
masking a partial run) — fix: reword to avoid `*)` substrings in comment
text. New pitfall #11 candidate.**
Pitfall #8 was promoted from candidate to
documented in `codex.md` "Common cross-engine pitfalls" item #1 before
III.3 launched. Pitfalls #6, #7 remain candidates; #9 (Mathematica
`Integrate[]` constant factoring), #10 (SymPy `nsolve` near
singularities), and #11 (Mathematica `*)` substring inside comment body
closes prematurely; verifier must check that all expected PASS lines appear,
not just `rc=0`) added in III.5. **IV.1 added no new pitfall candidates** —
all orchestrator-direct edits applied first-attempt clean (Cluster A docstring
carry-forwards + Cluster B 100 closure derivation + Cluster C 23-site banner
sweep), with verifier PASS-line counting confirming all expected substantive
checks across the 12 stages. **IV.2 added no new pitfall candidates** — all
orchestrator-direct edits applied first-attempt clean (Cluster A 106/109
docstring carry-forwards + Cluster B 108 β-locus extension + Cluster C 24-site
banner sweep), plus a Mathematica parse-bug correction (`chiArg /. beta -> 1 - 1`
parsed as `beta -> 0` at 108 F2; fix: `(chiArg /. beta -> 1) - 1`) which
re-confirms pitfall #11 PASS-line discipline (the buggy line passed by accident).
**IV.3 added pitfall #12 candidate** (Mathematica `Solve[expr == 0, frakG]`
fails with "frakG is not a valid variable" when `frakG` is bound to its definition
`= gQ*Sqrt[kS]/(gS*Sqrt[kQ])` — always introduce a fresh symbol for Solve's
target variable, then substitute back). Also re-confirmed: Mathematica `Minimize[]`
on `cos` over an interval frequently returns the input unevaluated; use boundary-value
checks instead under monotone-decreasing assumptions. Directive-correction during
edit: (a) stage 116 directive's `gamma0_from_D = -I·coeff(z,5)` produced wrong sign;
orchestrator corrected to `+I·coeff(z,5)` (analytic check confirms `+I` is right);
(b) stage 115 directive's `parentFamilyResidual - balanceEq*(kS*kQ)/(gS^2)`
multiplicative factor was wrong because `balanceEq` carries denominator
`kS*(kQ*kS + lam^2)`; correct factor is `(kS*kQ + lam^2)/(gS²*kQ)`. Both
corrections caught by re-running the scripts; sub-agents would otherwise have
silently produced FAIL residuals. See per-batch summaries in
`redteam/batches/batch_<ID>.md`.

### Linear projected-EM update

Snapshot addendum: `2026-05-11`

Stages `004--021` are now canonical linear stages for the projection-first
Maxwell integration, parent-action packet, and retained reduced one-port normal
form.  They are no longer counted as Stage `004` substages.  The old compact
Stage `004` reduced Maxwell/mixed calculation is retained as Stage `021`.

Stages `004--020` have file-for-file SymPy migrations from the derivation-only
`notes/em_projected` scripts through
`step_18_parent_throat_action_weak_axisym_packet`.  The `step_19_*`
branch-export packet and `step_20+` computational/runtime diagnostics remain
excluded from paper coverage.  File-for-file Mathematica mirrors for
Stages `004--020` have not been independently derived yet; Stage `021` retains
the legacy reduced Mathematica audit.

## Coverage By Part

| Part | Stage Range | Total | SymPy | Mathematica | Numerical | Review |
|---|---|---:|---:|---:|---:|---:|
| I | `001--023` | 23 | 23 | 6 | 2 | 6 |
| II | `024--036` | 13 | 13 | 13 | 0 | 13 |
| III | `037--090` | 54 | 53 | 54 | 2 | 54 |
| IV | `091--163` | 73 | 61 | 59 | 8 | 69 |
| V | `164--200` | 37 | 37 | 25 | 1 | 25 |
| VI | `201--218` | 18 | 18 | 2 | 0 | 2 |
| VII | `219--242` | 24 | 24 | 3 | 0 | 3 |
| VIII | `243--253` | 11 | 11 | 3 | 2 | 3 |

## Coverage Classes

| Coverage class | Count | Stage ranges |
|---|---:|---|
| SymPy + Mathematica | 163 | See `STAGE_PROVENANCE_INDEX.md` for the file-by-file list. |
| SymPy only | 77 | Includes projected-EM Stages `004--020` and later stages whose Mathematica mirrors are not yet present. |
| Mathematica only | 2 | `084`, `093` |
| No executable audit | 11 | `103`, `113`, `120`, `124`, `128`, `132`, `136`, `141`, `145`, `149`, `153` |

## Constant Provenance Rule

Coverage counts are not trust grades.

Likewise, `Mathematica present` is not the same thing as `independent second
CAS derivation`. Repo-wide counts track coverage breadth; independence is now a
separate policy classification in `MATHEMATICA_MIRROR_POLICY.md`.

For this archive, an audit should be treated as insufficient until every
nontrivial constant used in it is classified as one of:

- `derived in audit`
- `carried forward with source anchor`
- `probe-only numeric value labeled`

Any unexplained literal should be treated as a verification defect, not a style
issue.

## Immediate Gaps

### 1. No executable audit

These stages currently have neither SymPy nor Mathematica nor numerical stress:

`103`, `113`, `120`, `124`, `128`, `132`, `136`, `141`, `145`, `149`, `153`

These are the first places where the archive has mathematical content without an
executable backstop.

### 2. Mathematica without SymPy

These stages have a Mathematica artifact but no SymPy mirror:

`084`, `093`

These need reconciliation before we can make strong claims about dual-CAS
coverage.

### 3. SymPy-only frontier

The current SymPy-only region is:

`121--123`, `188--199`, `201--202`, `204--217`, `219--220`, `222--238`, `240--241`, `244--247`, `249--252`

Operationally, the main contiguous Mathematica mirror still ends at Stage `187`,
with isolated later checkpoints now hardened at Stages `200`, `203`, `218`,
`221`, `239`, `242`, `243`, `248`, and `253`.

### 4. Review-coverage gap

Stage-specific review notes are missing for:

`121--124`, `188--199`, `201--202`, `204--217`, `219--220`, `222--238`, `240--241`, `244--247`, `249--252`

That means the late-stage archive is not only missing Mathematica coverage, but
also lacks the earlier review-note pattern that exists for most of Stages
`001--187`.

### 5. Numerical stress remains sparse

Numerical stress coverage exists only for the following stage families:

`003/021`, `045--048`, `058/170`, `106`, `125`, `142--144`,
`146--148`, `150--152`, `154`, `155--156`, `157`, `185--187`,
`248`, `253`

This is a narrow slice relative to the symbolic verification surface and should
be treated as targeted spot-checking, not broad numerical validation.

## How To Use This Baseline

1. Use `STAGE_PROVENANCE_INDEX.md` when you need the exact note or script path
   for a stage.
2. Use this baseline when deciding which verification gaps are structural and
   which are only missing a second engine or review note.
3. Use the part-level counts to prioritize the next audit wave.
4. Use the coverage classes to define the future-paper citation support set.
5. Use `CITATION_SUPPORT_SET.md` when deciding which gaps are most important for
   downstream citation hardening.
6. Use `CHECKPOINT_CONSTANT_PROVENANCE.md` for the growing no-magic-numbers log
   on the checkpoint subset.
7. Use `CHECKPOINT_TRUST_AUDIT.md` for the current stage-level trust baseline
   on that checkpoint subset.

## Recommended Next Verification Sequence

1. Use `CHECKPOINT_TRUST_AUDIT.md` as the current checkpoint trust baseline.
2. Reconcile the remaining repo-wide Mathematica-only outliers:
   `084`, `093`, `103`, `113` (last two are IV.2 status-only).
3. Use `MATHEMATICA_MIRROR_POLICY.md` when deciding whether a mirror gap is an
   execution-coverage gap or an independence gap.
4. Then widen the audit wave to the remaining repo-wide gaps.
5. Backfill executable audits for the remaining no-executable stages:
   `128`, `132`, `136`, `141`, `145`, `149`, `153`. (103, 113, 120, 124
   cleared as status-only consolidation cards under IV.2/IV.3.)
6. Use the new `248` / `253` harnesses as the template if numerical-stress
   coverage is widened beyond the current spot-check set.
