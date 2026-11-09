# Project status — S15 reading build accepted; item 4's Singular engine in its build (2026-11-09)

## 2026-11-09: S15, and where the work stands

**S16 is accepted.** Its record and register pass were accepted at `4f907e9d` (round 2 clear, fresh Claude + Grok).
The superseded clause below still calls the S16 record "written for review"; that line is out of date.

**What is committed for S15.** S15's results stay **conditional**, on `R-S14a-01`, as the clause below says.
- The S15 decision list was cleared at `91248b67`. The S15 build was accepted at `3192123b` (round 1 clear). The
  comparator's decision list was cleared at `e305f6d4`, and the comparator was accepted at `5a019394` (round 1).
- The [reading decision list](research/pde_ledger_v3/directives/S15_reading_decision_list.md) was cleared at
  `7f52b4d2` (scoped pass 3). Its items 1–6 ask for the separation-departure objects on the power-law profile family.
  Item 4 asks, for each departure object in closed form, for the exact and complete conditions under which it vanishes
  for every separation.
- **The reading build.**
  - Rounds 0–4 were preserved as reviewed, at `5f0fc03a`, `6c0fc66b`, `4c4e2558`, `90441b89` and `0d9675eb`.
  - Between rounds, repair briefs 1–3 and repair notes 4 and 4b–4p were gated, checked or accepted, many of them
    carrying user decisions.
  - It was accepted at `cc62dc9f` after the round-5 review, where both legs (fresh Claude and Grok) were CLEAR.
  - Its [run record](research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md) lists each
    object and how it ended. The SymPy and Wolfram run manifests are beside it.
- **Item 4 has one engine in the accepted build: Wolfram.** By a user decision, the SymPy engine attempts no item-4
  object (note 4o, `f7068356`). Six objects carry the label `user_decision_2026_11_06` (note 4n, `266f2acf`): the
  two-drain region momentum and rate, for each Cartesian component. Their complete conditions were not obtained, and
  they are queued for a larger-memory follow-up in each engine.
- A Singular side test was gated at `07b1e751`. It ran, but its output has not been reviewed, so it is not evidence.
- **Singular becomes item 4's second engine** (user, 2026-11-09: "Go with your recommendation. Sounds like Singular is
  a solid choice for this"). The engine directive
  ([directive](research/pde_ledger_v3/directives/S15_item4_singular_engine_directive.md)) was handled as follows:
  - its reviewed pass-1 version was preserved at `51fd98fa`;
  - it was accepted at `d4cd6d09` after two CLEAR pass-2 legs (Codex and Grok).
  Singular 4.4.1 is installed at `/var/projects/_tools/singular`.
- Since `de9efa9b` (user, 2026-11-02), each run commits a run manifest and no engine transcript.

**In flight (uncommitted).**
- A Codex builder (`gpt-6-astra`) is building the Singular engine. Its source and the README's Singular section exist.
  Its run record and manifest do not exist yet.
- Its first session ended during its first guarded run. That run kept going, but its guard reservation was left with
  no live owner. The guard's explicit reconciliation rule must clear that reservation before it admits another job,
  and the builder is then resumed with an operational note.
- Before acting on any of this, check the guard's unit, the guard registry and the builder's DONE markers. Do not
  assume any of it is still running.

**NEXT, in order.**
1. Finish the Singular build. Its stale guard reservation is reconciled by the guard's rule, and the builder is
   resumed with an operational note. The builder then finishes its runs, and writes its run record, its manifest and
   its final report.
2. Check that the deliverable exists: the run record, the manifest and the final transcript. Then launch the build
   review. Codex wrote the build, so the legs are a fresh Claude agent and Grok, with one identical prompt, reviewed
   until clear. The Grok launch must let Grok take the guard's lock. If it cannot, the launch is fixed before the Grok
   leg runs: `AGENTS.md` forbids any fallback to an unguarded launch.
3. Adjudicate the findings. Commit the reviewed baseline before any repair, and accept only on clearance.
4. Build the two comparators: item 4, Wolfram against Singular; and the reading build's items 1–3, 5 and 6. Each needs
   a decision list (two legs, folded once), a build and a review.
5. Write the S15 step record and its register pass. The record asks which returned conditions are compatible with a
   steady sink that conserves mass. Riemann's sink ether and Lagally are oracles for that question, not premises.
6. Close S15:
   - If the six labelled objects are still not obtained, they go to
     [DEFERRED_HEAVY_RUNS](research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md) through the record.
   - 76 committed S15 files cite paths under the disposable scratch directory: 72 measurement files and 4
     directives. Most were written before the user's 2026-11-02 rule. The evidence they rely on is moved out or quoted
     in a manifest. The scope of that move is for the user to decide.
   - The step's local transcripts and their preserved copies are deleted.

**Also owed.**
- The Grok launch for guarded review legs. In round 5, Grok's sandbox could not take the guard's lock (recorded in
  the round-5 disposition).
- The unsigned commits made while the user was remote, to be re-signed.
- The `gin` remote's git branch is behind `origin`. `gin` holds the legacy annexed transcripts.

The clause below is SUPERSEDED — kept for the artifact map.

## Light sector closed for now; S12 and S14a accepted; S16 record written for review [SUPERSEDED by the clause above]

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

S14a part 2's [record](research/pde_ledger_v3/steps/S14a_2_slow_limit.md) and register pass are accepted
at `79511272` (round 2 clear), with [retrieval measurements](research/pde_ledger_v3/steps/_measurements/S14a_2_record_measurements.md).
The user chose **"Record as requirement"** (2026-10-21): no build ran. The required object is
the slow limit and far behaviour selected by the slab's own dynamics. The controlled
`ω→0` return law and profile-dependent `J→Q` map remain not derived. The filed witnesses retain
their separate pass/leg scopes as review instruments; they are not an accepted part 2 build.
The obligation has its own entry, `R-S14a-01`; `R-S12-02` retains its existing S12 scope.
The user's V2 choice **"All the slab's laws
(Recommended)"** makes selection owed jointly by S12's source/flux forms, S8/S22's material law
and the face/interior laws with no owner named ([decision](research/pde_ledger_v3/directives/_measurements/S14a_2_record_r0_review_disposition.md), `12aa65d5`). S14/S15/S20 remain
**CONDITIONAL, not carried results**, now on that named requirement. The two carried part-1
citation notes are repaired in its record; other part-1 claims retain their accepted standing.

S16's constructions (`1d990ee6`) and comparator (`f08db34a`, scoped citation review with the
containment claim) are accepted. Its [moving-region record](research/pde_ledger_v3/steps/S16_moving_region.md)
and register pass are written for source-first review, with
[retrieval measurements](research/pde_ledger_v3/steps/_measurements/S16_record_measurements.md).
The supplied balances give live regional count/momentum rates, count-centre relations, received
momentum and count moments through second order. They supply no mass, inertial/passive factor,
surroundings field, potential, inter-body force law or throat solve. Exact/resolved residuals,
domain/text differences and unpaired coverage retain the record's separate qualifications.
The pass brings 19 → 25 OPEN entries: the selected user options file passive-factor needs and
the inertial relation; the list author separately files compactness, face momentum, density response
and the S22 reduction-debt account. Each entry names its own source and admission ground.
The user's **2026-10-24 "Require all four (Recommended)"** decision now requires admission of
those four list-author filings: **"You require all four on the register now, and each entry cites
today's decision as its reason."** This settles the record-list item-7/item-9 disagreement with
the register's earlier rest-on applications. The original filers and owners remain; no law, form
or material is selected ([decision](research/pde_ledger_v3/directives/_measurements/S16_record_r1_review_disposition.md),
`4e5a94bc`; S16 measurements M19).
S16 also adds scoped provenance to `R-S1-01` (coordinate dimension, S6), `R-S1-02` (inherited
support identification, S1), `R-S12-02` and `R-S21-01`. The Y1–Y4 repair is 25 → 25, with
no new entry; no existing requirement, target or status changes. The reviewed version is preserved
at `fa1fa413` (round 1, not accepted). This W1 repair is also 25 → 25: no new entry or S16
consumer source, with the four added admission-decision citations. It remains for review,
with no acceptance claimed.

The plan's **"S16 must be run BEFORE S19 is trusted"** gate and **"its result feeds back into
the charter"** handoff remain OPEN beside this regional delivery. The cited plan rows name
**no owner** for disposition of either; S19 retains the separate C10/tail question. This record
does not decide the gate or charter feedback, and supplies no charter amendment or scope-boundary
decision ([record §4](research/pde_ledger_v3/steps/S16_moving_region.md); measurements M04/M17).

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
| `R-S14a-01`: slow limit and far behaviour selected by the slab's own dynamics | Jointly S12 (source/flux forms), S8/S22 (material law), and **no owner named** (face/interior laws); [register entry](research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md#r-s14a-01--the-dynamically-selected-slow-limit-and-far-behaviour). Record/pass accepted `79511272`; requirement OPEN. |
| S16 passive-field/separation/normalization, compactness and face-momentum requirements | **No owner named** for delivery; `R-S16-01/-02/-03`, [S16 record §6](research/pde_ledger_v3/steps/S16_moving_region.md). Compactness is also a debt for S22's table. |
| S16 momentum–material relation, density response and reduction-debt inventory | S8 (`R-S16-04`); S8/S22 (`R-S16-05`); S22 for the debt-table inventory (`R-S22-02`, without inventing physical suppliers). [S16 register pass](research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md). |
| Plan's S19 trust gate: "S16 must be run BEFORE S19 is trusted" | **OPEN; no owner named** for disposition in the cited plan rows. S19's separate C10/tail ownership is retained; no gate decision here ([S16 record §4](research/pde_ledger_v3/steps/S16_moving_region.md), M04/M17). |
| Plan's charter feedback: "its result feeds back into the charter" | **OPEN; no owner named** for disposition in the cited plan rows. No charter or scope-boundary decision here ([S16 record §4](research/pde_ledger_v3/steps/S16_moving_region.md), M04/M17). |
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
The S14a part 2 pass brings 18 → 19 OPEN entries. Its new `R-S14a-01` selection obligation
has V2's joint deliverers; no existing entry gains a source or changes target/status.
No existing entry covers the same selection requirement with those deliverers, so it is
separate from `R-S12-02`'s S12 forms/data and compatibility objects. The record/pass
were accepted at `79511272` (round 2 clear); no derivation or profile is supplied by recording this obligation.
The light closeout was accepted at `e070c7f3`; its review no longer gates S12. The oracle comparison's
questions about flow drive separate S12's source/boundary realization from O2's carried
force-balance handoffs. Material destination separates S12's drain-produced density/exchange
choice from S21's knit compatibility (`R-S21-01`). Each force states the flow it
needs; S21 compares them, and a mismatch reopens S9b (user, 2026-10-19).

The [throat/EM notes](docs/light_em_investigation_handoff.md) remain **exploratory and paused**, as input to
S22/Q2 (with S12/Q3 connections), not adopted constitutive laws. Mixed numerical recovery, centre-drive
implementation and new leakage workers remain parked. [Closeout, paused notes](research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md).

Cleanup is complete: history rebuilt and pushed at `bdc181c0`; archive tag `archive/pre-cleanup-2026-10-04` is on both remotes; nine uniform-cache runs are kept on local disk ([preservation list](research/pde_ledger_v3/cleanup_2026_10/FINDINGS.md#phase-5-preservation-list--keep-on-disk)).
