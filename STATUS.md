# Project status — light sector closed for now; S12 and S14a part 1 accepted; S14a part 2 requirement recorded for review

On 2026-10-19 the user closed the light sector **for now**. Its accepted results are conditional. What stays
open is carried as requirements in the [substrate register](research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md),
in its light-closing pass (accepted `e070c7f3`). S12's constructions (`a3da43d9`) and comparator
(`fe485b73`, with the engine-row note) are accepted. Its [step record](research/pde_ledger_v3/steps/S12_order_source.md)
and register pass are accepted at `ac909aaf` under "Accept, note both". Their known J1/J2 account/citation
items are repaired in the S14a part 1 record work.

S14a part 1's constructions (`94fddf2f`) and comparator (`de81b9fc`) are accepted. Its
[projection record](research/pde_ledger_v3/steps/S14a_1_projection.md) and register pass are accepted at
`4a4dfaf5` (round 1 clear), with [retrieval measurements](research/pde_ledger_v3/steps/_measurements/S14a_1_record_measurements.md).
S13 is folded into S14a (user, 2026-10-21). This part prints the projected balances and spherical
two-channel `J`, signed order loss, storage, finite-radius flux and conditioned velocity objects.
It preserves the comparator differences and coverage limits. It chooses no flow or conversion law;
its historical handoff left the controlled return law with later S14a and the profile-dependent
`J→Q` completion OPEN.

S14a part 2's [record](research/pde_ledger_v3/steps/S14a_2_slow_limit.md) and register pass are written
for source-first review, with [retrieval measurements](research/pde_ledger_v3/steps/_measurements/S14a_2_record_measurements.md).
The user chose **"Record as requirement"** (2026-10-21): no build ran. The required object is
the slow limit and far behaviour selected by the slab's own dynamics. The controlled
`ω→0` return law and profile-dependent `J→Q` map remain not derived. The filed witnesses retain
their separate pass/leg scopes as review instruments; they are not an accepted part 2 build.
The obligation merges into existing `R-S12-02`. The user's V2 choice **"All the slab's laws
(Recommended)"** makes selection owed jointly by S12's source/flux forms, S8/S22's material law
and the face/interior laws with no owner named ([decision](research/pde_ledger_v3/directives/_measurements/S14a_2_record_r0_review_disposition.md), `12aa65d5`). S14/S15/S20 remain
**CONDITIONAL, not carried results**, now on that named requirement. The two carried part-1
citation notes are repaired in its record; other part-1 claims retain their accepted standing.

**User decision, 2026-10-20: "Keep S12".** **"S12 stays open for the forms."** This record closes
only S12's requirements part and delivers none of the forms. The forms come once the brane is
defined, after the force sectors; every owner the sources and user gave is retained.

## Light sector at close

- **Accepted after S11c closed PARTIAL:**
  - O2, the steady brane balance: [record](research/pde_ledger_v3/steps/O2_steady_brane_balance.md), `72866fcf`;
  - S9b, conditional light bending and radar slope: [record](research/pde_ledger_v3/steps/S9b_light_bending.md), `ff3bac63`;
  - the S9b flow-only oracle check, `8d91f448`. Its [record](research/pde_ledger_v3/steps/S9b_flow_oracle_comparison.md) and the light-sector closeout are accepted at `e070c7f3`.
- **What the oracle check shows:** S9b's flow-only squared speed equals the river model's and Visser §7's.
  This holds only if S9b's orbital GM, clock and speed units are identified with each oracle's, and the check
  supplies none of these identifications. The signed velocity then matches for one radial orientation. Where
  the velocity is nonzero, no exchange makes the brane density fall as r^(−3/2). The check does not show
  that the brane is that fluid, and its optics are a supplied input.
- **New requirement R-S22-01:** how the energy a throat releases splits among light, the extra longitudinal
  wave and the bulk. It is needed to decide whether matter observes S11's extra wave, and whether the model conflicts
  with the blackbody, CMB and N_eff measurements (photon survey E32–E34, E38). S22 owns the split and
  Q2 the mouth/held datum; S12 retains the throat's order-conversion source and its boundary data,
  as the user chose on 2026-10-19 and retained under "Keep S12" on 2026-10-20. Those forms remain OPEN.
- **Parked:**
  - S11c (PARTIAL);
  - every heavy backlog run ([DEFERRED_HEAVY_RUNS](research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md));
  - the applicability checks of the other named oracles.

## S11c, closed PARTIAL

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
| Full cross-engine composition, b sign/control debts, c1 giants and accepted scattering/work/loss | S11c-b/c1/c2/d if reopened; the backlog runs are parked (2026-10-19), [DEFERRED_HEAVY_RUNS](research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md), [d limits](research/pde_ledger_v3/steps/S11c_d_profile_conditioned_scattering.md). |
| How a throat's released energy splits among light, the extra wave and the bulk | S22 with Q2 and S12, as chosen by the user on 2026-10-19 and retained on 2026-10-20; [R-S22-01](research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md). |
| Conversion observable and strong-edge interpretation | S11c-e; [split and d weak-limit obligation](research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md#actual-downstream-uses-and-next-step). |
| Dynamical drain/return source, boundary data and any change to the drain-frozen symmetry result | **S12, OPEN for the forms** under "Keep S12" (2026-10-20); [R-S12-02](research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md). This record delivers only the requirements part. No automatic transfer of the one-engine no-leak rule. |
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

**S12's record and register pass are accepted at `ac909aaf`:** its accepted constructions and comparator report the source
required by the committed order balance on live fields satisfying U1's support condition, including both material face
requirements and the labelled fast-ordering specialization. The six unrestricted face-coefficient
differences are labelled representation diagnostics confined off their respective delta supports,
not open differences in the S12 distribution; both raw expressions remain retrievable. The record
preserves all comparator coverage limits and the accepted engine-row
note. It supplies no drain law, throat solution or flow. Native source/controller forms, distinct
boundary data, native identifications and additional conversion partners retain their sourced S12
ownership and remain OPEN. None of their forms is delivered here; they come once the brane is
defined, after the force sectors.
The S12 pass left 18 OPEN entries,
with compatibility constraints merged into `R-S12-02` and the comparison of force sectors' flow
requirements owed by `R-S21-01`, consistently in its by-target row and requirement field.
[Record](research/pde_ledger_v3/steps/S12_order_source.md),
[measurements](research/pde_ledger_v3/steps/_measurements/S12_record_measurements.md).
The S14a part 1 pass also leaves 18 OPEN entries, adding scoped source/provenance to `R-S1-01`,
`R-S1-02`, `R-S12-02` and `R-S21-01`, with no new entry or status change. Its printed compatibility
constraints merge into `R-S12-02`; S21 retains the comparison of force sectors' flow requirements.
The S14a part 2 pass also leaves 18 OPEN entries, with no new entry or status change; its
selection obligation has the joint deliverers named by V2, while existing S12 obligations keep S12.
Only existing `R-S12-02` gains a source: the dynamically selected slow-limit/far-behaviour
requirement, its user-decision provenance and scoped filed-review evidence. The record/pass
await source-first review; no derivation or profile is supplied by recording this obligation.
The light closeout was accepted at `e070c7f3`; its review no longer gates S12. The oracle comparison's
questions about flow drive separate S12's source/boundary realization from O2's carried
force-balance handoffs. Material destination separates S12's drain-produced density/exchange
choice from S21's knit compatibility (`R-S21-01`). Each force states the flow it
needs; S21 compares them, and a mismatch reopens S9b (user, 2026-10-19).

The [throat/EM notes](docs/light_em_investigation_handoff.md) remain **exploratory and paused**, as input to
S22/Q2 (with S12/Q3 connections), not adopted constitutive laws. Mixed numerical recovery, centre-drive
implementation and new leakage workers remain parked. [Closeout, paused notes](research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md).

Cleanup is complete: history rebuilt and pushed at `bdc181c0`; archive tag `archive/pre-cleanup-2026-10-04` is on both remotes; nine uniform-cache runs are kept on local disk ([preservation list](research/pde_ledger_v3/cleanup_2026_10/FINDINGS.md#phase-5-preservation-list--keep-on-disk)).
