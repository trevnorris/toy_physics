# Project status — S11c PARTIAL; S12 next

The light model has a useful uniform result, but no accepted number for light lost at a defect.
S11c is closed **PARTIAL**. No nonuniform build, numerical continuation or review resumes automatically.
Read the [short closeout](research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md) first; the
[canonical d record](research/pde_ledger_v3/steps/S11c_d_profile_conditioned_scattering.md) gives the evidence and limits.

## What later work may use

- **CONDITIONAL:** uniform light decouples from the bulk in the supplied linear model. SymPy and independent
  Wolfram work support it, with documented engine/repair reviews and comparison limits; stability needs
  nonnegative stiffness. [S11b, transverse mode and review provenance](research/pde_ledger_v3/steps/S11b_interface_coupling_law.md).
- **CONDITIONAL:** a later SymPy check found finite modes, nonzero energy current and zero face drives
  (light does not push on the bulk) at, above and below the tested light/sound speed matches, with the bulk
  at rest and other inputs fixed. Claude/Grok cleared the method, not a fresh review of the executed result.
  [Uniform completion, coverage/limits](research/pde_ledger_v3/_measurements/S11c_d_near_unity_uniform_continue_result.md).
- **UNRESOLVED:** the held-profile omega=3 benchmark at ratio approximately 0.12 found no numerically resolved
  loss. Its 1 ppm empirical floor is not a physical upper bound; it is not a matched-speed defect result.
  [Benchmark, signed-current table/controls](research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_balance_report.md).
- **CONDITIONAL:** four equation repairs have scoped support, not full composition clearance. Original b/c2
  checks were on pre-repair equations. [b changes](research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md#changed-after-close),
  [c2 changes](research/pde_ledger_v3/steps/S11c_c2_self_energy_fold.md#changed-after-close).
- **CONDITIONAL:** the bounded S10/S11 Lean contracts retain their reviewed proof status; they do not supply
  the missing physical scattering/work result. [S10 closure](research/pde_ledger_v3/lean/s10/FIDELITY_REVIEW.md),
  [S11 contracts and application limits](research/pde_ledger_v3/lean/s11/README.md).

## Open questions and owners

| Question | Owner / evidence |
| --- | --- |
| Replaying the uniform check needs its intermediate caches regenerated from the kept b/c2 engines and the S10/S11 export chain repaired (H4). The original caches are not in git. They stay on local disk under `_scratch/` (Phase 5). | S11c-d if reused; [H3](research/pde_ledger_v3/cleanup_2026_10/PHASE3_REVIEW.md). |
| S10/S11 exports pin their pre-2026-09-11 producers. Regenerate them from the current producers under review, or revert the producer edits. The S10 edit may be what the S10 Lean D3 link relies on, so a revert is not automatic. | S10/S11 maintenance; [H4](research/pde_ledger_v3/cleanup_2026_10/PHASE3_REVIEW.md). Not run in this cleanup. |
| Direct c2 `[0,2]` entry is set to zero without an order-counting justification; whether the full response needs a nonzero term is unresolved. | S11c-c2/d if reused; [joint disposition, mixed-term issue](research/pde_ledger_v3/_measurements/S11c_upstream_repair_review_joint_disposition.md). |
| N6 premises: physical correctness of Φ, the V transform and extracted-block leakage remain open. | S11c-c2/d if reused; [c2, review caveats](research/pde_ledger_v3/steps/S11c_c2_self_energy_fold.md), [N6 disposition, §7](research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_disposition.md). |
| Full cross-engine composition, b sign/control debts, c1 giants and accepted scattering/work/loss | S11c-b/c1/c2/d if reopened; [DEFERRED_HEAVY_RUNS](research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md), [d limits](research/pde_ledger_v3/steps/S11c_d_profile_conditioned_scattering.md). |
| Conversion observable and strong-edge interpretation | S11c-e; [split and d weak-limit obligation](research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md#actual-downstream-uses-and-next-step). |
| Dynamical drain/return source, boundary data and any change to the drain-frozen symmetry result | **S12**; [plan, S12](research/pde_ledger_v3/V3_STEP_PLAN.md). No automatic transfer of the one-engine no-leak rule. |
| Real defect profile, holder, response and charge-force sign | Q2/Q3/S22; [plan, charge phase/S22](research/pde_ledger_v3/V3_STEP_PLAN.md). |
| Physical light/sound speed ratio | S20a calibration; S22/R10 proposed derivation; [plan, S20a/S22](research/pde_ledger_v3/V3_STEP_PLAN.md). |
| Material admissibility and unbuilt substrate, including shear response and physical dimension | S1–S8; [SUBSTRATE_REQUIREMENTS, entries/passes](research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md). S9/S10 requirements pass 1 is recorded, but the substrate steps have not run. |
| Muonium/gravity dependency quarantine, PN geometry-independence and structural audit | S16/species matching; [correction note, §§2–5](notes/muonium_gravity_corrections.md). No withdrawn throat geometry may supply a prediction. |

The following are **OPEN** carry-overs, not instructions to resume. [Part F](research/pde_ledger_v3/cleanup_2026_10/FINDINGS.md#f-status-debt-carry-over) maps every grouped item to the old STATUS and its disposition.

| Other debt | Owner / register |
| --- | --- |
| μ_S shear normalization | S11b/S11c-b if reused; [inertia report, final paragraph](research/pde_ledger_v3/_measurements/S11c_inertia_repair_report.md). |
| S11c-a full-axis KEYING and control-form characterization | S11c-a/b; [a record, standing limits/owed items](research/pde_ledger_v3/steps/S11c_a_interface_shape_derivatives.md). |
| F/G interpretation re-grounding; unresolved historical publication-review closure; unreviewed sheet/N6 trace repairs | S11c-c2/d if reused; [c2 record](research/pde_ledger_v3/steps/S11c_c2_self_energy_fold.md), [Part F conflicts](research/pde_ledger_v3/cleanup_2026_10/FINDINGS.md#carry-over-conflicts-and-questions). |
| b control namespace/admissibility coverage, kernel optimization, five off-path export checks | S11c-b/export maintenance; [Part F, carry-over evidence](research/pde_ledger_v3/cleanup_2026_10/FINDINGS.md#f-status-debt-carry-over). |
| S10 CAS/export/naming debts and unreviewed record corrections; S11 deferred loci/T7 checks; S11b comparison/passivity limits | S10/S11/S11b maintenance; [DEFECT_REGISTER, C19/C20 and WL entries](research/pde_ledger_v3/DEFECT_REGISTER.md), [REBUILD_HANDOFF, owed-item sections](research/pde_ledger_v3/REBUILD_HANDOFF.md), [S11b scope](research/pde_ledger_v3/steps/S11b_interface_coupling_law.md). Lean closure does not discharge these. |
| Review owed for CLAUDE.md L control and lean/FORMALIZATION_POLICY.md (`c2bdb663`) | Policy maintainers; [Phase 2 review, F6](research/pde_ledger_v3/cleanup_2026_10/PHASE2_REVIEW.md). Non-author S10 fidelity review is confirmed separately. |
| ≥64 GB remote compute setup and outstanding licensing/acceptance questions | Infrastructure owner; [plan](notes/remote_compute_plan.md), [requirements, §§8–10](notes/remote_compute_requirements.md). Planning only; no run authorized here. |
| Earlier v2 dimension, provenance, integration and normalization debts | Dormant v2 maintenance; [DIMENSION_REWRITE, §§8–12](research/pde_ledger_v2/manifests/DIMENSION_REWRITE.md), [Part F](research/pde_ledger_v3/cleanup_2026_10/FINDINGS.md#f-status-debt-carry-over). Not a revived v2 queue. |
| Remaining model gaps, including species identity and old EM boundary-operator limitations | Named owners in [DEFECT_REGISTER](research/pde_ledger_v3/DEFECT_REGISTER.md); Q2/Q3/S22 and later matching, with [Part F](research/pde_ledger_v3/cleanup_2026_10/FINDINGS.md#f-status-debt-carry-over) retaining historical limits. |
| Memory-note compaction and temporary-path evidence portability | Documentation maintenance; [Part F](research/pde_ledger_v3/cleanup_2026_10/FINDINGS.md#f-status-debt-carry-over). Phase 3 keeps cited outputs; uniform replay caches remain on local disk; no live-memory size remeasurement. |

## Next and paused work

**S12 is next in the physics plan:** define bulk-to-brane order conversion and its separate boundary data.
It does not require a finished S11c leakage factor. [Plan, S12](research/pde_ledger_v3/V3_STEP_PLAN.md).
This cleanup does not start S12 or any scientific run.

The [throat/EM notes](docs/light_em_investigation_handoff.md) remain **exploratory and paused**, as input to
S22/Q2 (with S12/Q3 connections), not adopted constitutive laws. Mixed numerical recovery, centre-drive
implementation and new leakage workers remain parked. [Closeout, paused notes](research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md).

Cleanup: Phase 3 is prepared on local `cleanup/pruned` for its STOP review; **Phases 4–5 remain pending**. G1/G2 and the supplied AGENTS replacement are applied. Selection, checks and export-pin drift introduced during the inventoried period are in [FINDINGS, Phase 3](research/pde_ledger_v3/cleanup_2026_10/FINDINGS.md#g-phase-3-selection-and-checks). No history rewrite or push. The archive tag `archive/pre-cleanup-2026-10-04` is unchanged and already on both remotes ([directive](research/pde_ledger_v3/CLEANUP_2026-10_directive.md)).
