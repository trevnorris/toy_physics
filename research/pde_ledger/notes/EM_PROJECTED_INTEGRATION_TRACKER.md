# EM projected integration tracker

Created: 2026-05-11
Updated: 2026-05-25 (batch I.1 v2 close)

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
