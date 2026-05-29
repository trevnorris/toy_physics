# EM projected integration tracker

Created: 2026-05-11
Updated: 2026-05-29 (IV.x/V.1 orchestrator-direct integrity remediation, batch 3 (stages 117/118/119/122) — **no change (batch 3 out of EM-projected range)**. These IV.3-range stages are core-balance / parent-core / parent-balance / mouth-source-compensation microphysics, NOT the projected-Maxwell core (004-021, closed at I.2); no projected-EM identities were touched, and `material_change: false` on every stage. Consistent with the batch-1 and batch-2 remediation snapshots, which were likewise out of the EM-projected range and produced no EM-tracker change. **Previous V.1 entry below.**)
Updated prior: 2026-05-28 (batch V.1 close — **V.1 closed the microscopic-logs / drifts / bundle-inversion range 164-175** with 12/12 verified, 22 findings closed across all 12 dirty stages, 0 blocked. No checkpoints, no status-only units. 1 material_change (stage 170 — additive Sec.5 weak-axisymmetric (1,1/2,-1) signature coverage; no carried/derived result changed, zero downstream propagation). **EM-projection relevance: low this batch** — V.1's content is grouped-P2 / bundle-inversion / outlet-map microphysics (log-imbalance channels, lower-branch drift laws, bundle inversion/transport, off-bundle slippage, grouped outlet map, slope collapse, weak-axisymmetric loading mismatch, static self-similarity, wall-normalized load shape), NOT the projected-Maxwell core (004-021, closed at I.2). No projected-EM identities were touched. **Cluster A only** (stale `STAGE N-17` script banner across all 12 stages; 11 fixed in per-finding directives, 4 residuals mass-fixed); no Cluster B (body-text citation re-attribution) and no Cluster C (paper-card Checks downgrade) this batch. One non-cluster paper_misalignment (170 F1) resolved by adding the missing weak-axisymmetric signature check to both engines (user direction a) — first non-cluster paper_misalignment in the v2 run. `mathematica_transliteration` 9/12 (7 independent routes added, 2 policy mirrors 169/175). Thirteenth consecutive batch clear of stop-cold. Three orchestrator catches: 166 round-trip vector-residual scalarization, 175 F1 minimal-resolution switch (cross-check was a simplify-commutes identity), 171 bundle-route rework (directive route was tautological). **Previous IV.6 entry below.**)

Updated prior: 2026-05-28 (batch IV.6 close — entire range 001-163 now paper-aligned at v2 depth. **IV.6 closed the correction / coevolution / traction / off-family range 151-163** with 13/13 verified, 19 findings closed across 10 dirty stages (151/152/154/155/156/157/158/160/161/163), 3 clean (153 status-only consolidation card, 159 and 162 substantive-clean). 19 resolved + 0 blocked. No checkpoints in IV.6. **One material_change** (stage 151 SymPy rewrite from symbolic-integration to mpmath numerical because `sp.integrate(Sigma_star·cos(πx/2)·R, (x,0,1))` with free `Pi_star` hung >30 min CPU; downstream impact zero since the script verifies algebraic identities not numeric carry-forwards). 3 user-gate clusters resolved (all "Recommended"): **Cluster A** systematic stage-number renumberings (18 banner edits across 9 stages × 2 engines at −17 offset; 151/152/157 had no offset; 2 notes H1 edits at −85 offset for stages 159/160; 20 mechanical edits via `/tmp/iv6_mass_renumber.py`); **Cluster B** body-text forward-stage citation re-attribution (53 citations across 13 notes; three offsets: −51 for IV.4 references 188-199, −85 for IV.6 internal cross-references 239-248 — **new offset specific to IV.6** disambiguated from −102 by content cross-check, −102 for 221 IV.3 reference and 249-250 IV.5 references; applied via `/tmp/iv6_reattribute.py`); **Cluster C** stage 158 paper-card Checks downgrade (items 2 even-preservation and 3 tangent motion δ⊥=0 rewritten as forward-carry citations of `\ref{stage:159}` and `\ref{stage:162}` / `\ref{stage:163}`; **first FORWARD (downstream) carry-forward downgrade in v2 batches** — IV.4 134 and IV.5 144 both carried upstream). Twelfth consecutive zero-redirection batch. Three orchestrator catches in the rework loop: (1) stage 151 SymPy hang → mpmath rewrite (first v2 batch where the auditor-prescribed engine approach was infeasible at the engine level and the orchestrator had to redesign the verification); (2) stage 154 multivariate `Series` retained cross-products → switched to single-epsilon parameterization `piExprEps /. epsLin -> 1`; (3) stage 161 directive variable-substitution typo (`depsg_direct` instead of `depsg_branch`) → fixed in both engines. **Previous IV.5 entry below.**

Updated prior: 2026-05-27 (batch IV.5 close — entire range 001-150 now paper-aligned at v2 depth. **IV.5 closed the susceptibility / branch / defect-transport range 139-150** with 12/12 verified, 31 findings closed across 9 dirty stages (139/140/142/143/144/146/147/148/150), 3 clean (141, 145, 149 — all status-only consolidation cards per MANIFEST `is_status_only_candidate`). 30 resolved + 1 `blocked_legitimate` (144 F4 transliteration policy, user-accepted). No checkpoints in IV.5. **Zero material_change** in IV.5 — every fix was renumbering, notes re-attribution, paper-card downgrade, or script-substance addition; no derived numerical constants moved. 3 user-gate clusters resolved (all "Recommended"): **Cluster A** systematic stage-number renumberings (11 `.wl` banners offset by −17 with LEDGER variants, 6 `.py` banner sites for stages 142/143/144 with LEDGER variants, 4 notes/ H1 lines offset by +102 for stages 146-149; 21 mechanical edits via `/tmp/iv5_mass_renumber.py`); **Cluster B** body-text forward-stage citation re-attribution (22 citations across 11 of 12 notes; two offsets: −51 for 188-199 range matching IV.4 body offset and −102 for 220-251 range matching IV.5 notes-H1 offset; applied via `/tmp/iv5_reattribute.py`); **Cluster C** stage 144 paper-card Checks downgrade (items (i) outlet-consistency and (ii) self-matched-susceptibility rewritten as carry-forward citations of `\ref{stage:135}` and `\ref{stage:140}`, mirroring IV.4's stage 134 pattern). Eleventh consecutive zero-redirection batch. **Pitfall #13 re-confirmed** at stage 139 (Mathematica comment `Pi_*)` substring parsed as comment-terminator; ASCII-safe rewrite required). Four orchestrator catches in the rework loop: (1) stage 148 directive-prescribed closed form `4107 - 168π²` was wrong (numeric ~1.547); auditor copied stage 148 notes' typo; orchestrator caught via Mathematica `FullSimplify` yielding a symbolic Sqrt residual, cross-referenced stage 126 upstream notes (correct `4107 - 100π²` ≈ 0.184), and corrected both engines + the stage 148 notes typo. (2) Stage 139 pitfall #13 recurrence (comment rewrite). (3) Stage 142 SymPy tolerance `1e-20` → `1e-15` for nsolve's actual `~1.95e-18` precision. (4) Stage 147 SymPy `sp.N(AT)` default-15-digit truncation → explicit `sp.N(AT, 30)`. Plus stage 146 needed two SymPy numeric-sample fallbacks (Pi-sample for F2 `Integrate` unevaluation, eps-sample for F1 `simplify` non-reduction). **Previous IV.4 entry below.**

Updated prior: 2026-05-27 (batch IV.4 close — entire range 001-138 paper-aligned at v2 depth. **IV.4 closed the penetration / mouth-boundary / fixedpoint range 127-138** with 12/12 verified, 22 findings closed across 10 dirty stages (127/130/131/132/133/134/135/136/137/138), 2 clean (128, 129 — 128 is status-only consolidation card; 129 had only a cosmetic banner note not raised as a formal finding). No checkpoints in IV.4. **Zero material_change** in IV.4 — every fix was renumbering, notes re-attribution, paper-card downgrade, or script-substance addition; no derived numerical constants moved. 4 paper_misalignment items + structural cluster resolutions handled via 3 user-gate clusters: **Cluster A** systematic stage-number renumberings (9 `.wl` banners offset by -17 and 10 notes/ H1 lines offset by +102, plus 5 `.py` banner sites for stages 127, 133×2, 134, 135; all fixed via `/tmp/iv4_mass_renumber.py` script in 19 mechanical edits); **Cluster B** status-only carry-forward re-attribution at stages 132 (`Stages 180-182` → `Stages 129-131`) and 136 (`Stages 184-186` → `Stages 133-135`), correcting the broken chain where notes attributed load-bearing constants to *downstream* stages; **Cluster C** stage 134 paper-card Checks downgrade (items 1 outlet-consistency and 2 susceptibility-closure rewritten as carry-forward citations of `\ref{stage:135}` and `\ref{stage:137}`). Tenth consecutive zero-redirection batch. **Pitfall #13 candidate added** (Mathematica parser fails on comment substring `g'(Pi_*)` adjacent to `*)` — workaround is ASCII labels). One directive-correction catch by orchestrator at edit time: stage 134 F1's literal target values for `S_q(1/2)`, `S_q(1)`, `S_q(2)` were FABRICATED by the auditor agent — orchestrator recomputed via mpmath at 50 digits and substituted the verified values (0.608336415687717…, 0.633127670034487…, 0.681366857005321…). One mechanical follow-up squashed in: stage 134 notes line 140 trailing-9 typo (`605429` → `605428` to match boxed forms at lines 86/92 and the scripts). **Previous IV.3 entry below.**

batch IV.3 close — entire range 001-126 now paper-aligned at v2 depth. **IV.3 closed the core-balance / DtN-mixed / outlet / positive-source range 115-126** with 12/12 verified, 27 findings closed across 10 dirty stages (115/116/117/118/119/121/122/123/125/126), 2 clean (120, 124 — status-only consolidation cards). No checkpoints in IV.3. 7 paper_misalignment items resolved via 4 user-gate clusters: **Cluster A** notes-side numerical typos (121/122/126 `168π² → 100π²`; 123 `228 → 160`) — fixed in place across 4 notes files; **Cluster B** substantive script-side λ sign flip at stage 118 (internal section IV vs section V disagreement caught by auditor — section V had dropped a minus during integration; flipped to align with section IV's bilinear and the notes); **Cluster C** parametric family extension at stage 125 (integral inequality coverage gap — added `σ_a(z) = (a+1)(z/L)^a/L` family in both engines with endpoint asserts); **Cluster D** stage 117 consolidation card (3 blocked items resolved via cite-upstream-and-downgrade — added comment block citing stages 115/116, replaced misleading `expect_zero` tautological wrappers with `print("carrying forward...")` lines, wired `classification_rows` booleans from sections 1-5 residuals). **One material_change at stage 118** (λ sign flip; downstream Schur reductions all use squared combinations, so no numerical propagation; upstream_stale NOT flagged). Ninth consecutive zero-redirection batch. Pitfall #12 candidate added (Mathematica `Solve[expr == 0, frakG]` fails when `frakG` is bound to its definition — use fresh symbol). Two directive-correction catches by orchestrator at edit time: 116's `-I·coeff(z,5)` → `+I·coeff(z,5)`, and 115's multiplicative factor `(kS·kQ)/gS²` → `(kS·kQ + lam²)/(gS²·kQ)`. **Previous IV.2 entry below.**

batch IV.2 close — entire range 001-114 was paper-aligned at v2 depth; **IV.2 closed the outgoing DtN / deformation / robustness / robin range 103-114 including checkpoint 105 (chi_Q fix from outgoing DtN, 10 PASS lines on Mathematica via fully independent path) and added a substantive β-parameterized preservation submanifold at stage 108 (Cluster B) extending the previously β=1-only verification to the full notes-boxed `Σ_5 = S(1-β^5)/9 - Σ_0/27` locus. 5 paper_misalignment items resolved via 3 user-gate clusters (Cluster A docstring carry-forwards at 106/109 naming upstream stages 102/104 + downstream stages 110/111/112; Cluster B substantive β-locus addition at 108; Cluster C 24-site banner sweep across all 10 IV.2 scripted stages with 12 paper card display titles deferred). One structural material_change at stage 108 (verification surface widened; no value change). Eighth consecutive zero-redirection batch. No new pitfall candidates; pitfall #11 PASS-line discipline re-confirmed by 108 F2 `chiArg /. beta -> 1 - 1` Mathematica parse bug catch.**; batch IV.1 close — entire range 001-102 was paper-aligned at v2 depth; **IV.1 closed the grouped-P_2 geometry / isotropic decoupling / contamination range 091-102 including checkpoint 096 (geometry-lane check verdict, 20 PASS lines = 15 orthogonality + 5 final ledger) and the headline `mhat_0^2 chi_Q N_Q = 1` outgoing-normalization factorization at stage 100 strengthened from tautological cross-check to substantive observable-condition closure derivation (impose `mhat_0^2 * Gamma_5 = Gamma_5_target` on script-derived Gamma_5(K_0, chi_Q, Omega), derive closure ratio = `mhat_0^2 chi_Q N_Q - 1`). 10 paper_misalignment items resolved via 3 user-gate questions (Cluster A docstring carry-forwards for 091/095/097/099 naming stage 094 orthogonality + stages 088-090 minimal-module + stages 091/092/096 static-limit anchors; Cluster B closure derivation + 100/101 docstrings; Cluster C 23-site script-side banner sweep with 12 paper section titles deferred to future paper-cleanup pass). One structural material_change at stage 100 (verification surface; no value change). Seventh consecutive zero-redirection batch. 100 F4 design-level Mathematica transliteration rewrite remains `blocked_legitimate`. No new pitfall candidates.**; the projected-EM core 004-021 was closed at I.2 v2; II.1 v2 closed the post-EM bookend 022-036; III.1 v2 closed the rank-2 / continuum-kernel range 037-048; III.2 v2 closed the tracking / zeta-threshold / asymmetry / boost range 049-060 including checkpoint stage 051 clean; III.3 v2 closed the microclosure / gain-thresholds / equilibrium / walls range 061-072 including checkpoint stage 069 clean and v1 material_change at stage 068 confirmed clean; III.4 v2 closed the Family-1 geometry / thresholds / quadrupole range 073-084 with zero material_change and 4 substantive paper-misalignment fixes — 074 alpha `128/sqrt(5)` → `111/sqrt(5)`, 075 Upsilon_w `117 Theta_w` → `100 Theta_w`, 082 closed-form `zeta_phys` pin + Family-1 instantiation, plus 11-stage banner-relabel sweep aligning all III.4 self-banners with the post-renumber numbering; **III.5 closed the quadrupole-cancellation / loading-ratio / verdict range 085-090 including both checkpoint stages 089 and 090 verified, two substantive paper-misalignment fixes — 087 F1 status/checkpoint consolidation (script docstrings now name upstream sources rather than re-deriving), 089 F1 chain closure adding Omega(Pe→0)=1 + zeta_F1(0)=A_F1 limits + explicit Pe_req=0 assertion locking the paper-side boxed Output — plus 12-script banner-relabel sweep ("STAGE 68/069/70/070/71/071/72/072/73/073" → "STAGE 085/086/087/088/089/090") and two orchestrator hot-fixes on stage 088 (SymPy `omega**2/Omega_Q**2 -> u` subs failed after sp.simplify reshape; Mathematica comment containing `stage085_*)` substring prematurely closed comment causing silent partial run — verifier caught it from missing PASS lines)**)

## Scope decision

Use `notes/em_projected` as derivation source only through the derivation notes
ending at step 18.  The directory intentionally has no step 06 file.  Do
not import the later computational/runtime material into the PDE ledger paper:

- step 19 branch-export / Galerkin negative-control packet;
- step 20 and later reduced-family scans/frontiers;
- electron/SI falsification, dimensional runtime maps, CFD postprocessors,
  fail-fast classifiers, snapshot adapters, monopole JSONL screens, and related
  runtime tools;
- current simulator/CFD verdicts.

Those later files may remain useful as private diagnostics, but they are not
paper derivation support in this integration pass.

## Current ledger structure

The projected Maxwell derivation is now represented as ordinary global stages:

| Stage | Source note | Ledger role |
|---:|---|---|
| `004` | `step_01_projected_maxwell_readme.md` | Bundle index and derivation boundary. |
| `005` | `step_02_projected_maxwell_covariant_notes.md` | Covariant projected Maxwell law. |
| `006` | `step_03_projected_maxwell_vector_notes.md` | Vector projected Maxwell split. |
| `007` | `step_04_projection_reduction_comparison_notes.md` | Projection/reduction comparison. |
| `008` | `step_05_projected_maxwell_extension_notes.md` | Extension terms. |
| `009` | `step_07_projected_maxwell_near_throat_notes.md` | Near-throat packet. |
| `010` | `step_08_projected_maxwell_push_bundle_master_notes.md` | Push-bundle master. |
| `011` | `step_09_projected_maxwell_p2_bridge_notes.md` | Grouped `P_2` bridge. |
| `012` | `step_10_projected_maxwell_stage4_primitive_bridge_notes.md` | Primitive bridge. |
| `013` | `step_11_projected_maxwell_mouth_taylor_master_notes.md` | Mouth-Taylor master. |
| `014` | `step_12_projected_maxwell_mouth_taylor_gate_bridge_notes.md` | Mouth-Taylor gate bridge. |
| `015` | `step_13_parent_throat_action_master_notes.md` | Parent throat action master. |
| `016` | `step_14_parent_throat_action_candidate_notes.md` | Parent throat action candidate. |
| `017` | `step_15_parent_throat_action_weak_axisym_notes.md` | Weak-axisymmetric parent limit. |
| `018` | `step_16_parent_throat_action_bundle_master_notes.md` | Parent action bundle master. |
| `019` | `step_17_parent_throat_action_isotropic_bundle_notes.md` | Isotropic parent bundle. |
| `020` | `step_18_parent_throat_action_weak_axisym_packet_notes.md` | Weak-axisymmetric packet. |
| `021` | `notes/stages/moving_throat_pde_stage021_reduced_one_port_normal_form.md` | Reduced one-port normal form retained from the prior EM reduction stage. |

Former post-EM ledger stages start at Stage `022`.

## Verification policy

- SymPy audits are top-level files named
  `scripts/moving_throat_pde_stageNNN_*_sympy_audit.py` and should match the
  paper stage order.
- Projected Maxwell Stages `004--012` have independent Mathematica mirrors
  landed during red-team audit batch I.1 (2026-05-21), and Stages `013--020`
  have independent Mathematica mirrors landed during red-team audit batch I.2
  (2026-05-21). Stage `021`'s reduced one-port Mathematica audit was
  substantially rewritten in batch I.2 (manual EL derivation replaced with
  `EulerEquations` from `VariationalMethods`; a recurrence of the I.1 stage
  003 multi-line `lRed = ...` continuation parsing defect was caught and
  patched via parenthesization; Sections II.2/III/V rewritten to use
  `LinearSolve`, an analytic-derivative route, and `SphericalHankelH1[2, z]`).
- Mathematica audits must be regenerated one script at a time.
- Paper `\StageFile{...}` entries should resolve to real source files and to
  regenerated transcript files after the audit runs.

## Completed verification checks

- Full SymPy audit regeneration completed:
  `TOTAL: 241  PASS: 241  FAIL: 0  SKIPPED: 0`, including the repo-level
  master audit and `240` stage-level audits.
- Current Mathematica transcript set regenerated one script at a time:
  `TOTAL: 165  PASS: 165  FAIL: 0`.
- Numerical stress transcript summaries are current for both Python and
  Mathematica harnesses.
- Missing `\StageFile{...}` and stale old-numbering scans were run against
  paper-facing files and active metadata after the linear renumbering.
- `paper/pde_ledger.tex` was rebuilt in the post-renumbering verification pass.
- Red-team audit batch I.1 (stages 001-012, "Part I.1 -- Geometry lift, BdG
  coupling, projected Maxwell setup") completed 2026-05-21 with all 12 stages
  reaching `verified` and both engines independently checking each load-bearing
  claim. `material_change: false` on every stage, so no upstream cascade. See
  `redteam/batches/batch_I1.md` plus per-stage reports at
  `redteam/reports/stage_NNN.md` and `redteam/verifications/stage_NNN.md`.
- Red-team audit batch I.2 (stages 013-023, "Part I.2 -- Maxwell bridge,
  parent throat action, reduced one-port") completed 2026-05-21 with all 11
  stages reaching `verified` and both engines independently checking each
  load-bearing claim. Stages 013-020 received new independent Mathematica
  mirrors; stage 021's mirror was substantially rewritten; stages 022-023
  verified without material changes. `material_change: false` on every stage,
  so no upstream cascade. See `redteam/batches/batch_I2.md` plus per-stage
  reports at `redteam/reports/stage_NNN.md` and
  `redteam/verifications/stage_NNN.md`.
- Red-team audit batch II.1 (stages 024-036, "Part II.1 -- Overlap isotropy
  through continuum kernel") completed 2026-05-22 with all 13 stages
  reaching `verified` and both engines independently checking each
  load-bearing claim. `material_change: false` on every stage, so no
  upstream cascade. See `redteam/batches/batch_II1.md` plus per-stage
  reports at `redteam/reports/stage_NNN.md` and
  `redteam/verifications/stage_NNN.md`.
- Red-team audit batch III.1 (stages 037-048, "Part III.1 -- Continuum
  kernel, generalized branch, rank-2") completed 2026-05-22 with all 12
  stages reaching `verified` and both engines independently checking
  each load-bearing claim. 10 of 12 stages required codex edits (27
  findings total: 10 `mathematica_transliteration`, ~11
  `tautological_check`, 4 `insufficient_verification`, 1
  `hardcoded_result`, 1 `symbol_assumption_error`); stages 042 and 048
  verified clean on first read. Zero codex iter-2 fixes needed.
  `material_change: false` on every stage, so no upstream cascade. See
  `redteam/batches/batch_III1.md` plus per-stage reports at
  `redteam/reports/stage_NNN.md` and `redteam/verifications/stage_NNN.md`.
- Red-team audit batch III.2 (stages 049-060, "Part III.2 -- Tracking,
  zeta thresholds, asymmetry, boost") completed 2026-05-22 with all 12
  stages reaching `verified` and both engines independently checking
  each load-bearing claim. 11 of 12 stages required codex edits (27
  findings total: 13 `tautological_check`, 6 `mathematica_transliteration`,
  5 `insufficient_verification`, 3 `hardcoded_result`); stage 056 verified
  clean on first read. Two batch-wide toolchain patches landed: an
  `expectZero` helper update to strip `ConditionalExpression[0, ...]`
  wrappers that `Solve`/`Reduce` introduce under aggressive `$Assumptions`,
  and a `1/pi1 == 0` infinity test replacing strict `pi1 =!= Infinity`
  in stage 051 because Mathematica's `Limit` returns inconsistent forms.
  Zero codex iter-2 fixes needed. One stage (060) flagged
  `material_change: true` due to a Csol restructuring; downstream Xi_micro
  consumers in batches III.3+ are still `pending` so no immediate cascade.
  See `redteam/batches/batch_III2.md` plus per-stage reports at
  `redteam/reports/stage_NNN.md` and `redteam/verifications/stage_NNN.md`.
- Red-team audit batch III.3 (stages 061-072, "Part III.3 -- Microclosure,
  gain thresholds, equilibrium, walls") completed 2026-05-22 with all 12
  stages reaching `verified` and both engines independently checking each
  load-bearing claim. 10 of 12 stages required codex edits (27 findings
  total: 14 `tautological_check`, 9 `mathematica_transliteration`, 3
  `insufficient_verification`, 1 `hardcoded_result`); stages 061 and 066
  verified clean on first read. Zero codex iter-2 fixes needed. One stage
  (068) flagged `material_change: true`: `Wfail_res`/`Wfail_match` now
  derived via `Solve` from explicit resonance-corrected premises rather
  than postulated, with derived expressions matching the prior postulated
  forms symbolically. III.3 is out of the linear projected-EM core range
  (004-021), so no downstream cascade into this tracker's scope. See
  `redteam/batches/batch_III3.md` plus per-stage reports at
  `redteam/reports/stage_NNN.md` and `redteam/verifications/stage_NNN.md`.
- Red-team audit batch III.4 (stages 073-084, "Part III.4 -- Family-1
  geometry, thresholds, quadrupole") completed 2026-05-25 with all 12
  stages reaching `verified` and both engines independently checking each
  load-bearing claim. All 12 of 12 stages required codex edits (first
  batch with no clean-on-first-read stages); 40 findings total (14
  `tautological_check`, 12 `hardcoded_result`, 7 `mathematica_transliteration`,
  7 `insufficient_verification`) — `hardcoded_result` rose sharply because
  the Family-1 numerology cluster 075-084 accumulates many literal
  constants. Zero codex iter-2 fixes needed; one orchestrator-applied
  mid-batch hot-fix on stage 081 (standard `ConditionalExpression[e_, _]
  :> e` strip retrofitted after `Solve` introduced the wrapper). Two
  acceptable codex deviations (stage 076: `P = K*rho^n_poly` over the
  directive's literal form; stage 078: removed a spurious `100` factor in
  the directive — both verified necessary). Zero `material_change: true`
  flags this batch — every derivation-route rewrite left printed symbolic
  and numeric content byte-identical. III.4 is out of the linear
  projected-EM core range (004-021), so no downstream cascade into this
  tracker's scope. See `redteam/batches/batch_III4.md` plus per-stage
  reports at `redteam/reports/stage_NNN.md` and
  `redteam/verifications/stage_NNN.md`.
- Red-team v2 paper-grounded re-audit of batch I.1 (stages 001-012, "Part
  I.1 -- Geometry lift, BdG coupling, projected Maxwell setup") completed
  2026-05-25. **First batch processed under the v2 auditor**, which reads
  `paper/stages/stage_NNN.tex`, per-stage notes under `notes/stages/`, and
  the part-level appendix BEFORE opening scripts. A 10th finding category
  `paper_misalignment` was added with subtypes; findings of that category
  route to user resolution rather than codex. All 12 stages reached
  `verified`. **5 substantive defects v1 missed**: 7 paper_misalignment
  items across stages 001/006/007/010/011, plus 3 new script-side
  tautologies (004 Faraday/Bianchi symbol-substitution, 008
  `Integrate[W*Z]` self-cancel, 012 M1 carried-forward primitive
  self-checks). All 7 paper_misalignment items resolved via a structured
  Codex-as-math-authority workflow (questions markdown → Codex fills
  recommendation blocks → user approves → Codex apply session edits paper
  + scripts). Two stages flagged `material_change: true`: 001 (source-
  coupling and gauge-fix sign flips per paper notation firewall
  `eta = diag(-1,+1,+1,+1,+1)`) and 004 (Faraday/Bianchi block now
  exercises a real cyclic identity via Schwarz commutativity instead of
  pure symbol substitution). Stage 010 paper card grew ~100 lines (δP_n
  display equations + 4 new cluster anchor paragraphs); stage 011 script
  trimmed (5 clusters moved attribution to stages 022-024 which already
  publish them); stage 007 scripts gained ~75 lines per engine (H(w)
  profile + `xi_eff^proj` checks for the gauge channel the paper
  promised). One new toolchain pitfall documented in `codex.md`
  (`Part[]`-on-pattern-parameter inside `Do[Module[...]]` silently drops
  half the body; discovered on stage 004; fix is precomputed
  immediate-valued `F_ij` expressions before any Do/Module scope opens).
  Since I.1 IS in the linear projected-EM core range (004-021), the v2
  sweep also serves as a paper-alignment certification for the first
  half of that range; the remaining 013-021 will be covered when batch
  I.2 is re-swept. See `redteam/batches/batch_I1_v2.md`,
  `redteam/resolutions/batch_I1_paper_alignment.md` (Codex's 7
  recommendations with rationales), and per-stage reports/verifications
  under `redteam/reports/stage_NNN.md` and `redteam/verifications/stage_NNN.md`.
- Red-team v2 paper-grounded re-audit of batch I.2 (stages 013-023, "Part
  I.2 -- Maxwell bridge, parent throat action, reduced one-port")
  completed 2026-05-25, **certifying paper-alignment for the back half
  of the linear projected-EM core range (013-021) plus the
  grouped-bundle bookends (022, 023)**. All 11 stages reached `verified`
  (3 clean: 016, 017, 022). 6 paper_misalignment items surfaced (Q1-Q6
  across stages 013, 014, 015, 018, 020, 021). Pattern notably different
  from I.1: dominated by **(b) trim scripts**, because cross-stage
  cross-check revealed duplication — the EM-projected scripts were
  file-for-file ports of `notes/em_projected/step_NN_*` master notes,
  and the later compact paper cards distributed content across multiple
  stages. Specific duplications resolved by trim: stage 013's
  `δP_2`/`δP_4`/sieve (destinations: stage 010 owns δP_n per I.1 v2
  paper add; stage 014 owns sieve); stage 014's
  Xi_load/`δP_n`/Compat (destinations: stage 013 owns Xi_load; stage
  010 owns δP_n + compatibility transport); stage 015's
  wall-only/Y20/grouped (~half of each engine — destination stage 017
  owns lane signature + b=3a + wall-only obstruction); stage 018's
  one-pole/gate/Xi_1 (destinations: stage 019 owns one-pole closure +
  normalization compatibility; stage 020 owns gate determinant +
  wall-slope solve + Xi_1 residual); stage 020's Y20 block
  (destinations: stages 010, 017 both own Y20 lane ratios). The sixth
  item (Q6 on stage 021) was the inverse direction —
  script_missing_paper_claim — paper Output enumerated three exports
  but scripts only asserted two-and-a-half; resolved by **adding** a
  composed `δD_2^(odd)(ω)` assertion in both engines. After verifier
  wave caught the initial composed assertion as tautological (used
  bare `N0` symbol on RHS instead of closed form), a Codex remediation
  pass replaced the RHS with the Section III closed form
  `(Ω_A² g_W + R g_A)² / (Ω_A² Ω_W² - R²)²` — assertion now
  non-tautologically requires both Section III's N(0) derivation and
  Section IV's `Γ_5^port = a^5/(27 c_s^5)` to be correct. Stage 023's
  Schur-complement F2 finding produced a real Codex iteration:
  directive's `+R_mix` prescription conflicted with existing
  `+2*g_U*g_W*R_mix` rational numerator; Codex deferred, then a
  remediation pass had Codex read paper's `eq:app-stage023-full-lagrangian`
  to derive the sign convention from physics — option α chosen
  (`+R U W` Lagrangian → `-R` off-diagonal in frequency-space spring
  matrix → existing numerator sign confirmed correct). 4 stages
  flagged `material_change: true` (013, 014, 015, 018 — all script-side
  scope reductions removing duplicates). User re-direction was load-
  bearing this batch: Codex initially recommended (c) acknowledgement
  for Q1/Q2; user pushed "each step builds on prior" principle; cross-
  check by orchestrator confirmed destination scripts already verified
  the content → Q1/Q2 flipped to (b) trim before apply. Apply prompt
  added a **destination-verification guardrail** ("grep destination
  script to confirm equivalent assertion exists before deleting from
  source") that worked as designed — no orphan trims occurred. With
  batch I.2 v2 closed, **the entire linear projected-EM core range
  (004-021) is now paper-aligned at v2 depth**, plus the grouped
  bookends 022-023. The 4 material_change cascades are contained within
  the same batch's verified stages (destinations 010/014/017/019/020
  are all in batches I.1 v2 or I.2 v2). See
  `redteam/batches/batch_I2_v2.md`,
  `redteam/resolutions/batch_I2_paper_alignment.md` (Codex's 6
  recommendations with rationales, plus user-revised Q1/Q2), and
  `redteam/resolutions/codex_remediation_batch_I2.md` (the math-decision
  remediation pass for 021 and 023).
