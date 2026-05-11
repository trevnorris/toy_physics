# EM projected integration tracker

Created: 2026-05-11

## Scope decision

Use `notes/em_projected` as derivation source only through the parent-action
weak-axisymmetric packet, steps 1--18.

Do not import the later computational/runtime material into the PDE ledger paper:

- step 19 branch-export / Galerkin negative-control packet;
- step 20 and later reduced-family scans/frontiers;
- electron/SI falsification, dimensional runtime maps, CFD postprocessors,
  fail-fast classifiers, snapshot adapters, monopole JSONL screens, and related
  runtime tools;
- current simulator/CFD verdicts.

Those later files may remain useful as private diagnostics, but they should not
be cited as derivation support in the paper integration pass.

## Current PDE ledger structure

The paper is assembled from:

- `research/pde_ledger/paper/pde_ledger.tex`
- `research/pde_ledger/paper/main_parts.tex`
- `research/pde_ledger/paper/parts/part01_parent_geometry.tex`
- canonical stage files under `research/pde_ledger/paper/stages/`
- stage appendices under `research/pde_ledger/paper/appendices/`

The EM issue is concentrated in Part I:

- `part01_parent_geometry.tex` imports the localized Maxwell equation and has a
  short `Projection and reduction hooks` subsection.
- `stage_004.tex` is explicitly a one-lane reduced Maxwell/mixed model.
- `stage_005.tex` and `stage_006.tex` consume Stage 004 as the source of the
  outgoing transfer factor and grouped `D,N,P` bundle.
- the Stage Appendix Part I text describes Stage 004 as adding the localized
  Maxwell/mixed sector through a reduced outgoing bridge.

The current status is honest but incomplete: it marks the Maxwell/mixed bridge
as exact only inside the declared reduced model and reduced relative to the full
localized Maxwell PDE. The new projected work can strengthen the parent-to-bundle
chain by explaining how the projected mixed sector feeds the same `Z_n` and
`N_n` slots instead of treating the reduced one-lane block as the primary EM
derivation.

## Source material to import

Core projected Maxwell:

- `step_02_projected_maxwell_covariant_notes.md`: exact projected inhomogeneous
  law with boundary/leakage term and projected charge continuity.
- `step_03_projected_maxwell_vector_notes.md`: measured `(E,B)` versus
  source-coupled `(D,H)` fields; projected inhomogeneous laws include leakage
  and gauge-driver terms.
- `step_04_projection_reduction_comparison_notes.md`: projection-first
  effective coupling versus reduction-first coupling; equality requires matched
  source/localization structure.
- `step_05_projected_maxwell_extension_notes.md`: localized gauge fixing
  `H=Z`, effective `mu0` and `xi`, and exact matching channel
  `S=Z/Z_int`, `H=Z`.
- `step_07_projected_maxwell_near_throat_notes.md`: finite throat / mouth
  projection; the live term is `<partial_w(Z F^{w nu})>`; mouth kernels give
  first corrections at `O(ell)`.

Projection-to-PDE bundle bridge:

- `step_08_projected_maxwell_push_bundle_master_notes.md`
- `step_09_projected_maxwell_p2_bridge_notes.md`
- `step_10_projected_maxwell_stage4_primitive_bridge_notes.md`
- `step_11_projected_maxwell_mouth_taylor_master_notes.md`
- `step_12_projected_maxwell_mouth_taylor_gate_bridge_notes.md`

Parent throat action / parent packet:

- `step_13_parent_throat_action_master_notes.md`
- `step_14_parent_throat_action_candidate_notes.md`
- `step_15_parent_throat_action_weak_axisym_notes.md`
- `step_16_parent_throat_action_bundle_master_notes.md`
- `step_17_parent_throat_action_isotropic_bundle_notes.md`
- `step_18_parent_throat_action_weak_axisym_packet_notes.md`

## Main findings

1. Projection-first Maxwell is not ordinary reduced brane Maxwell. It is an
   open-system brane electrodynamics unless transverse mixed leakage and
   gauge-driver terms are suppressed.

2. The exact projected inhomogeneous law carries the mixed-sector term
   `<partial_w(Z F^{w nu})>`. This is precisely the parent-level channel that
   the existing reduced Stage 004 represents only indirectly by a mixed
   coordinate.

3. Reduction-first Maxwell is still useful as a matched limiting channel. With
   `H=Z` and a source profile matched to `Z/Z_int`, projection-first zero-mode
   Maxwell reproduces the reduction-first coupling and gauge parameter. Without
   those conditions, the projected coupling is observer/source-profile
   dependent.

4. Near the throat, the mouth-local projection is structurally stronger than an
   interior symmetric projection: the first correction is `O(ell)` rather than
   `O(sigma^2)`. This is the main mathematical reason the current ledger should
   not collapse the EM sector directly to the far-field reduction before the
   mouth/bundle bridge.

5. The projected mouth correction feeds exactly the existing grouped-bundle
   slots:
   `Z_0,Z_2,Z_4` and `N_0,N_2,N_4`.
   It therefore plugs into the ledger without changing the downstream
   `D,N,P`, `P_0`, `P_2`, `P_4`, and `Xi_1` bookkeeping.

6. The primitive one-port bridge gives a useful sorting map:
   - `Q', Delta', P'` move the weak-axisymmetric prefactor slope;
   - `Q', Delta', S_2', H_port'` move the conservative even gates;
   - `G_W'` first enters the constant-prefactor transport.

7. The mouth-Taylor gate bridge says a single EM scalar correction is not
   enough. A nontrivial conservative even-gate repair requires a mixed
   compensation surface; the projected EM contribution is multi-channel.

8. The parent throat action promotion is related but not the same issue. It
   turns wall coefficients into parent-action integrals. The projected EM work
   should be integrated alongside this because later branch tests need a clear
   wall, support, conservative mixed, and outgoing-transfer provenance, but the
   actual branch-export diagnostics are not part of this Stage 004 derivation.

## Recommendation

Do not delete the reduction-first Maxwell material outright. Replace its role.

The current Stage 004 should no longer be the foundational EM derivation. It
should become the matched/reduced one-port normal form used after the
projection-first parent law has been established. In other words:

- projection-first Maxwell should be the parent-level EM law;
- reduction-first Maxwell should be retained as a controlled matched limit and
  as the one-port algebraic normal form for `Z_n` and `N_n`;
- downstream Stages 005--006 can mostly stay because they consume only the
  abstract grouped `D,N,P` slots, but their inputs/provenance should be revised.

## Proposed edit plan

1. Add a new Part I subsection before the current localized Maxwell/mixed
   bridge:
   "Projection-first localized Maxwell law".
   This should contain the exact projected law, leakage term, field-layer split,
   projected continuity, and the `H=Z` gauge-localization result.

2. Rewrite the current "Projection and reduction hooks" subsection. It is too
   small for the role it now needs to play. It should explicitly distinguish:
   exact projection, matched zero-mode reduction, and mouth-local projection.

3. Replace the opening of the current Stage 004 / Part I Maxwell bridge. The
   one-lane reduced Lagrangian should be introduced as a normal form for the
   projected mixed-sector packet, not as the primary derivation of the EM
   sector.

4. Add Stage-004 substage files with zero-padded suffixes so `ls` preserves the
   intended local order without renumbering the global ledger.  Proposed file
   pattern:
   - `paper/stages/stage_004_001_projected_maxwell_readme.tex`
   - `paper/stages/stage_004_002_projected_maxwell_covariant.tex`
   - `paper/stages/stage_004_003_projected_maxwell_vector.tex`
   - `paper/stages/stage_004_004_projection_reduction_comparison.tex`
   - `paper/stages/stage_004_005_projected_maxwell_extension.tex`
   - `paper/stages/stage_004_006_projected_maxwell_near_throat.tex`
   - `paper/stages/stage_004_007_projected_maxwell_push_bundle_master.tex`
   - `paper/stages/stage_004_008_projected_maxwell_p2_bridge.tex`
   - `paper/stages/stage_004_009_projected_maxwell_stage4_primitive_bridge.tex`
   - `paper/stages/stage_004_010_projected_maxwell_mouth_taylor_master.tex`
   - `paper/stages/stage_004_011_projected_maxwell_mouth_taylor_gate_bridge.tex`
   - `paper/stages/stage_004_012_parent_throat_action_master.tex`
   - `paper/stages/stage_004_013_parent_throat_action_candidate.tex`
   - `paper/stages/stage_004_014_parent_throat_action_weak_axisym.tex`
   - `paper/stages/stage_004_015_parent_throat_action_bundle_master.tex`
   - `paper/stages/stage_004_016_parent_throat_action_isotropic_bundle.tex`
   - `paper/stages/stage_004_017_parent_throat_action_weak_axisym_packet.tex`
   - `paper/stages/stage_004_018_reduced_one_port_normal_form.tex`

   Matching audit files should use the same suffix:
   - `scripts/moving_throat_pde_stage004_001_projected_maxwell_readme_sympy_audit.py`
   - `scripts/moving_throat_pde_stage004_018_reduced_one_port_normal_form_sympy_audit.py`

   The canonical `paper/stages/stage_004.tex` should become the wrapper that
   inputs or summarizes these substage files and preserves the public Stage 004
   label.

5. Add a new derivation block after Stage 004 or inside Stage 004:
   "Projected Maxwell to grouped-bundle slots".
   It should state the shift
   `Z_n -> Z_n + eps z_n`, `N_n -> N_n + eps n_n`, the induced changes to
   `u_2,u_4,P_0`, and the key cancellation that `z_0` drops out of the
   eliminated isotropic compatibility surface.

6. Add the primitive/mouth-Taylor bridge as a theorem or proposition:
   map `Q',S_2',H_port',Delta',P',G_W'` to the live bottlenecks
   `Xi_load`, `K_1`, and `H_even`; state the mechanism sieve and mixed
   compensation surface at a summary level.

7. Update `stage_004.tex`, `stage_005.tex`, `stage_006.tex`, and the Part I
   appendix prose to reflect the new provenance. The downstream formulas can
   mostly survive, but the status language should shift from "reduced
   Maxwell/mixed model" to "projection-first parent law plus reduced/matched
   one-port normal form".

8. Update source/provenance appendices so `notes/em_projected` steps 1--18 are
   listed as derivation sources, while step 19 and steps 20+ are explicitly
   marked out of paper scope.

## Readiness checklist before implementation

The conceptual plan is ready, but the implementation pass should include these
mechanical checks so the new substages behave like first-class ledger entries.

### File ordering and LaTeX assembly

- Keep `paper/stages/stage_004.tex` as the public Stage 004 wrapper.
- Put substage files beside it using `stage_004_001_*`, `stage_004_002_*`, ...
  so `ls paper/stages/stage_004*` gives the local derivation order.
- Have the wrapper use explicit `\input{stages/stage_004_00N_*}` calls.  Do
  not rely on globbing or generated include order.
- Keep `\label{stage:004}` only in the wrapper.  Give substages local labels
  such as `subsec:stage004-projection-first-maxwell` to avoid duplicate stage
  anchors.

### Audit artifact naming

- Use the same suffix in SymPy and Mathematica filenames:
  `moving_throat_pde_stage004_001_*_sympy_audit.py` and
  `moving_throat_pde_stage004_001_*_mathematica_audit.wl`.
- Save transcripts under the existing output folders with matching basenames:
  `scripts/output/*.txt` and `mathematica/output/*.txt`.
- Add the new transcript paths to the Stage 004 verification field once the
  audits have been run.

### Audit harness check

- Before trusting the repo-wide SymPy runner, verify or repair
  `research/pde_ledger/scripts/run_all_audits.sh`.  Its current path logic still
  points at an older `scripts/moving_throat` layout, while the live audits are
  under `research/pde_ledger/scripts`.
- The Mathematica runner already scans `research/pde_ledger/mathematica` using
  the expected `moving_throat_pde_stage*_mathematica_audit.wl` pattern, so the
  proposed substage filenames should be discovered there.
- For the first implementation pass, run new audits individually by exact path
  before attempting any repo-wide audit sweep.

### Tables and maps to update

- `notes/STAGE_PROVENANCE_INDEX.md`: Stage 004 should list multiple substage
  note/audit artifacts rather than only the old reduced Maxwell note.
- `notes/STAGE_VERIFICATION_COVERAGE.md`: update counts if the new audits are
  accepted as Stage 004 coverage; otherwise add a note that they are Stage 004
  substage coverage and avoid double-counting the global stage count.
- `paper/appendices/stage_appendix_part01.tex`: revise the Part I narrative and
  Stage 004 row to mention projection-first Maxwell plus the reduced one-port
  normal form.
- `paper/frontmatter/04_result_anchor_map.tex`: split `MTDC-T4` subanchors so
  projection-first law, matching/reduction, and one-port normal form have
  distinct citation handles.
- `paper/frontmatter/05_dependency_map.tex`: update `MTDC-T4` dependencies so
  later blocks inherit projection-first EM unless they explicitly invoke the
  matched reduction.
- `paper/appendices/source_file_index.tex`: list `notes/em_projected` steps
  1--18 as derivation sources and mark step 19 plus steps 20+ as excluded
  computational or branch-diagnostic material.

### Status and claim firewall

- The new projection-first law can be `Exact` inside the localized parent
  Maxwell equation once the projection kernel is declared.
- The matched zero-mode reduction remains `Reduced / Controlled Reduction`.
- The mouth-local bundle bridge is `Exact Within Closure` for the declared
  mouth-kernel/primitive packet and remains `Open` for actual branch
  realization.
- Do not let the new projected derivation upgrade the downstream normalization
  target to realized.  It only improves provenance and identifies the correct
  exported branch packet.

### Exclusions

- Do not import `notes/em_projected` step 19 branch-export diagnostics or step
  20+ scans/frontiers/runtime tools into the paper derivation.
- If any later computational result is mentioned in a private planning note,
  keep it out of result anchors, stage verification fields, and source-file
  citation tables.

## My current opinion

The projected EM material is not a cosmetic add-on. It changes the correct
logical order of the EM sector.

The old reduction-first EM derivation is still useful and should remain in the
ledger, but it should be demoted from "this is the EM sector" to "this is the
matched zero-mode/one-port normal form used once the projected parent law is
collapsed to bundle coefficients." That preserves all the good Schur-complement
and transfer-factor algebra while fixing the conceptual problem: near the
throat, the mixed transverse projected term is live and cannot be erased before
the PDE has exported the branch packet.

The safest paper integration is therefore replacement-by-layering, not simple
deletion:

1. insert projection-first Maxwell as the controlling parent EM layer;
2. reinterpret the reduced Maxwell/mixed block as a derived normal form;
3. keep the downstream grouped-bundle algebra;
4. import the mouth-Taylor bridge as the new theorem gate explaining what the
   actual branch must export.

## Implementation status

Executed on 2026-05-11.  The canonical Stage 004 file is now a wrapper over
ordered substages `stage_004_001_*` through `stage_004_018_*`.  Substages
`001--017` are file-for-file migrations of the derivation-only
`notes/em_projected` SymPy scripts through step 18 into the root ledger order.
Substage `018` preserves the old reduced Maxwell/mixed calculation as the
matched one-port normal form.  Step 19 and steps `20+` from `notes/em_projected`
remain out of paper scope.

Verification results are recorded in
`research/pde_ledger/notes/EM_PROJECTED_INTEGRATION_PROGRESS.md`.
