# S11c-era inventory and findings

## Where the light sector stands

S11c asked how much light escapes into other motion when the brane's thickness or material properties change. It closes **PARTIAL**: there is no accepted loss number. **UNRESOLVED.** [Closeout, opening](../steps/S11c_PARTIAL_CLOSEOUT.md).

On a uniform brane, small light waves do not drive the surrounding bulk in the stated linear model; stable propagation requires nonnegative stiffness. SymPy and independently written Wolfram calculations support this, with engine/repair reviews and a documented comparison; one Wolfram review's full rerun was blocked. **CONDITIONAL.** [S11b, transverse mode, comparison and operational note](../steps/S11b_interface_coupling_law.md).

A later SymPy check found the selected light modes still finite and decoupled at the tested light/sound speed matches, with the bulk at rest and other inputs fixed. Claude and Grok cleared its method, not a fresh independent review of the executed result. **CONDITIONAL.** [Uniform check, coverage and limits](../_measurements/S11c_d_near_unity_uniform_continue_result.md).

Where the brane changes, the equations contain routes for conversion, but its size is unknown. The one completed scattering benchmark could not distinguish its final signal from numerical error; “within 1 ppm” is not a physical bound. **UNRESOLVED.** [Benchmark, signed-current table and controls](../_measurements/S11c_d_numerical_radiating_balance_report.md).

Four equation repairs have scoped review support, not approval of the whole nonuniform calculation. One direct mixed contribution is set to zero without an order-counting justification; whether the complete response needs a nonzero term remains unknown. **UNRESOLVED.** [Joint repair disposition, limits and mixed-term issue](../_measurements/S11c_upstream_repair_review_joint_disposition.md).

**OPEN:** S11c-b/c2/d owns equation/computation debts if reopened; S11c-e owns the conversion observable and strong-edge interpretation; Q2/Q3/S22 owns real defect support and charge response, S22/R10 the speed ratio, and a separate material audit the rotational-stiffness assumption. **S12 is next**, for drain/return dynamics; it does not need a finished leakage number. [Closeout, downstream uses and MacCullagh](../steps/S11c_PARTIAL_CLOSEOUT.md); [plan, S12/Q2/Q3/S22](../V3_STEP_PLAN.md).

## Inventory counts (metadata, not scientific findings)

**ESTABLISHED — inventory metadata:** [INVENTORY.tsv](INVENTORY.tsv) contains **4,434 A + 24 M = 4,458** paths in `dada3b7d..archive/pre-cleanup-2026-10-04`, whose endpoint is `2f56b303b144ffb3078d7b4f2b6f95b486bbc25b`. Reproduce with [inventory.py](inventory.py), functions `inventory`, `workstream`, `kind`; `--check` verifies the TSV without writing. It reads Git metadata, not scientific objects. Each row has one primary workstream and kind. Regular-file bytes are the archived blob size; **132** annex symlinks use payload size encoded in their keys, not pointer length. No annex content was fetched. This is the specified A/M interval, not a census of the whole repository or a keep list. Source: [directive, Phase 1a](../CLEANUP_2026-10_directive.md#phase-1-inventory-and-findings--stop), TSV and script.

Kind key: **K1** spec/amendment/contract; **K2** decision list or build directive; **K3** result report; **K4** review report or disposition; **K5** script that produces a reported result; **K6** worker/continuation/tooling; **K7** process record; **K8** output (`.out`/`.json` result); **K9** Lean proof/contract; **K10** Lean CAS bridge; **K11** doc/note; **K12** ledger front matter. These are the directive's roles, assigned by the script's explicit precedence rules; a role is not acceptance. JSON receipts remain K7 even when named `result_record`; scientific input/check payloads are grouped with K8. Bridge generators/bindings take K10. Source: [directive, Phase 1a](../CLEANUP_2026-10_directive.md), [script, `kind`](inventory.py).

| Workstream | K1 | K2 | K3 | K4 | K5 | K6 | K7 | K8 | K9 | K10 | K11 | K12 | Total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 01 Physics specification and contracts | 6 | 3 | 1 | 59 | 1 | 2 | 29 | 4 | 0 | 0 | 2 | 0 | 107 |
| 02 Symbolic build and omega1 premises | 1 | 224 | 148 | 70 | 95 | 199 | 449 | 347 | 0 | 0 | 36 | 0 | 1569 |
| 03 Upstream repairs and diagnostics | 0 | 7 | 15 | 19 | 13 | 56 | 98 | 118 | 0 | 0 | 7 | 0 | 333 |
| 04 Numerical scattering and omega3 benchmark | 2 | 3 | 6 | 5 | 4 | 44 | 67 | 23 | 0 | 0 | 6 | 0 | 160 |
| 05 Near-unity and defect prerequisites | 21 | 31 | 14 | 170 | 7 | 135 | 760 | 131 | 0 | 0 | 14 | 0 | 1283 |
| 06 Clean-condition packet | 2 | 2 | 0 | 11 | 19 | 59 | 14 | 77 | 0 | 0 | 3 | 0 | 187 |
| 07 Lean and S9-S11 formalization audits | 3 | 0 | 37 | 56 | 5 | 45 | 159 | 126 | 235 | 70 | 10 | 8 | 754 |
| 08 Exploratory throat and EM | 0 | 0 | 0 | 7 | 0 | 0 | 14 | 2 | 0 | 0 | 6 | 0 | 29 |
| 09 Muonium gravity corrections | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 2 |
| 10 Execution infrastructure | 0 | 0 | 2 | 0 | 0 | 9 | 9 | 1 | 0 | 0 | 1 | 3 | 25 |
| 11 Ledger front matter and other notes | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 6 | 9 |
| Total | 35 | 271 | 223 | 397 | 144 | 549 | 1599 | 829 | 235 | 70 | 89 | 17 | 4458 |

Workstream refinements are organizational: 02 includes the early finite omega=1 response that grew out of the symbolic build; 03 includes upstream diagnostic copies; 07 includes adjacent S9/S10 audits as well as Lean; 08 includes the exploratory documents' supporting reviews. The inventory does not infer a scientific status from a filename. Source: [script, `workstream`](inventory.py).

## B. Findings by workstream

### 01 — Physics specification, amendments and contracts

**Question:** what should the supplied-profile calculation mean? [Shared physics, §0](../directives/S11c_d_SHARED_PHYSICS.md).

- **ESTABLISHED — document review:** v10 received literal **SOUND** from fresh Claude/Opus and Grok; this is specification review. It supplies a localized, two-ended interface, independent retained grades and a full coupling vertex, while carrying c2 cross-engine operand debt. Sources: [Claude, opening](../directives/_legs/S11c_d_shared_physics_review_r10_opus.md), [Grok, verdict](../directives/_legs/S11c_d_shared_physics_review_r10_grok.md), [spec, §§0–1b](../directives/S11c_d_SHARED_PHYSICS.md).
- **WITHDRAWN/FAILED:** unrestricted pole/residue/projector identification was replaced by `nonlinearPoleV2`; the repair's controls are dimensionless synthetic examples, not physical S11c poles. Source: [pole repair report, paragraphs 1–5](../_measurements/S11c_d_nonlinear_pole_repair_report.md).
- **ESTABLISHED — scoped review:** Option B's amendment/inventory received Claude and Grok **CLEAR** in round 9; the total-loss amendment received **CLEAR FOR THIS PREPARATION AMENDMENT** from both in round 2. Neither clears execution/results. Sources: [Option B disposition, opening](../_measurements/S11c_d_scattering_form_review_round09_adjudication.md), [total-loss disposition, opening](../_measurements/S11c_d_total_transverse_loss_review_disposition.md).
- **OPEN — S11c-d/e if reopened; S22/R1:** actual response, supported power/normalization and finite-contrast interpretation remain owed; the spec explicitly separates a computable supplied-profile FORM from the real magnitude needing throat input. Sources: [spec, §0 out-of-scope and §3d](../directives/S11c_d_SHARED_PHYSICS.md), [Option B disposition, “Whole-plan assessment”](../_measurements/S11c_d_scattering_form_review_round09_adjudication.md).

### 02 — Symbolic build and omega=1 premises

**Question:** can the imported operator yield modes, currents and a two-ended response? [Program brief, §§1–2](../directives/S11c_d_sympy_build_PROGRAM_BRIEF.md).

- **CONDITIONAL:** the SymPy-led finite response reached four full-rank **645-unknown** cases and four incident directions, with case-specific maps and retained grades. Positive regulator, approximate boundaries and unresolved tiny signals remained. This is not a completed blind dual-engine d calculation. Sources: [builder report, “Four-case response completion” and “Retained user-approved solver/export contract”](../_measurements/S11c_d_sympy_builder_report.md).
- **CONDITIONAL:** at the approved omega=1 LAB_HELD/RHO4_CONSTANT input, open-thickness selector rank was **zero**, and the recorded real-momentum bulk-depth relation closed propagation. Zero open-thickness flux there is not a general no-leak theorem. Source: [continuum-current report, opening and final paragraph](../_measurements/S11c_d_continuum_currents_report.md).
- **WITHDRAWN/FAILED:** the fixed-point continuation supported selected REFERENCE/LEFT premises but stopped on RIGHT with `MemoryError`; exit zero was not an all-pass outcome. Its two build verdicts were **CLEAR FOR THIS PROGRESS-DEPENDENT FIXED-POINT CONTINUATION**, not result reviews. Source: [fixed-point continuation report, opening, failure and review paragraphs](../_measurements/S11c_d_transverse_face_fixed_point_continue_report.md).
- **OPEN — S11c-d if reopened:** final scattering/FORM integration, physical pole coverage, blind Wolfram/T7 and final export were not supplied by these milestones. Sources: [builder contract](../_measurements/S11c_d_sympy_builder_report.md#retained-user-approved-solverexport-contract), [Option B disposition, final section](../_measurements/S11c_d_scattering_form_review_round09_adjudication.md).

### 03 — Upstream repairs and their reviews

**Question:** did the inherited equations and their consumers represent the stated action/face laws? [Repair review disposition, repair table](../_measurements/S11c_upstream_repair_review_disposition.md).

- **CONDITIONAL:** saved repairs address inertia sign (`a74da30a`), mechanical-load orientation (`c643112a`), physical/reference pressure trace (`a05b05e3`) and thickness coordinate (b **`3b52afcb`**; c2 **`f618178a`**, export **`537d78fd`**). They have distinct domains and regeneration records; the coordinate diagnosis changed `W_bg e_W,t` to the physical `W0 e_W,t`. Sources: [inertia report](../_measurements/S11c_inertia_repair_report.md), [mechanical report, opening](../_measurements/S11c_mechanical_repair_report.md), [trace report, construction/checks](../_measurements/S11c_c2_trace_repair_report.md), [coordinate report, source repair](../_measurements/S11c_thickness_coordinate_repair_report.md).
- **UNRESOLVED:** overall nonuniform composition remains uncleared. Claude's scoped table says **“Clear on stated domain”**, **“Orientation clear; non-flat face-row content insufficient evidence”**, and **“Two-leg clear; three-leg mixed grade needs revision”** as applicable. Supplemental Grok's overall verdict is **“Needs revision”**, with scoped clear findings. Sources: [Claude disposition, table](../_measurements/S11c_upstream_repair_review_disposition.md), [joint disposition, opening/table](../_measurements/S11c_upstream_repair_review_joint_disposition.md).
- **UNRESOLVED — unreviewed repairs:** the separate Wolfram audit/repair corrected a reference-pressure double shift by extracting the actual face-law reference-pressure map and solving its ordered inverse on the tested N6 domain. Its report says no review or comparator ran; it is not a full self-energy/d engine. The d-internal sheet repair (`f9e28f5f`) resolves its bounded fixed-positive-frequency counterexample; its report explicitly calls the instruments **unreviewed**. Sources: [Wolfram repair, final paragraphs](../_measurements/S11c_wolfram_pressure_trace_repair_report.md), [sheet repair, opening and final section](../_measurements/S11c_d_sheet_repair_report.md).
- **OPEN — S11c-b/c2/d if reused:** direct mixed-term completeness and full cross-engine composition remain debts; see C2–C3 below. [Joint disposition, “Remaining mixed-term issue”](../_measurements/S11c_upstream_repair_review_joint_disposition.md).

### 04 — Numerical radiating track and omega=3 benchmark

**Question:** is the total transverse-current deficit numerically resolved for one held profile? [Benchmark, opening](../_measurements/S11c_d_numerical_radiating_balance_report.md).

- **UNRESOLVED:** strict rest bulk, LAB_HELD/RHO4_CONSTANT, omega=3, tangents **(1/5,1/10)**, selected phase-speed ratio approximately **0.12** (bare S9 ratio separately **0.1**), one fixed profile. Ten solves completed. Column 2/3 deficits in ppm were baseline **-1.773911 / +1.633609**, refined **-1.034869 / +0.942457**, larger interval **-0.608245 / +0.554788**, then regulator-halved **-0.608218 / +0.554792**. The **1 ppm** floor is empirical, not an upper bound. Source: [benchmark, scope and signed-current table](../_measurements/S11c_d_numerical_radiating_balance_report.md).
- **UNRESOLVED — review:** Claude and Grok both literally said **NEEDS REVISION** for the original-tolerance route. Actual matrix-route Gaussian, momentum-zero and saved trial joins remained unperformed. The coarse alternative was unused. Sources: [replacement disposition, verdict table and findings 1–6](../_measurements/S11c_d_numerical_radiating_replacement_review_disposition.md), [benchmark, controls](../_measurements/S11c_d_numerical_radiating_balance_report.md).
- **OPEN — S11c-d if reopened; S12 for flow:** no accepted gain/loss, physical bound, matched-speed defect result or draining-model conclusion transfers from these arrays. Tiny solve residuals do not establish physical accuracy. Evidence dependencies include the report's saved matrices/end maps and native source chain, not just its final table. [Benchmark, controls and ratio-1 handoff](../_measurements/S11c_d_numerical_radiating_balance_report.md).

### 05 — Near-unity, packet, receiving and centre prerequisites

**Question:** which concrete obstacles to a matched-speed calculation can be removed? [Uniform result, opening](../_measurements/S11c_d_near_unity_uniform_continue_result.md), [defect status, “Physics refocus”](../_measurements/S11c_d_defect_near_unity_applicability_status.md#physics-refocus--2026-10-04-user-directed).

- **CONDITIONAL:** omega=3 rest-bulk LAB_HELD/RHO4_CONSTANT, fixed tangents **(1/5,1/10)**: **12 exact speeds**, **48 sign checks**, two selected polarizations; LEFT match **sqrt(3/2)**, RIGHT **sqrt(150/101)**. Modes/current are finite and selected face drives vanish on those samples/limits. Both method verdicts: **CLEAR FOR THIS SELECTED UNIFORM METHOD**; no fresh worker/result review. [Uniform completion, coverage and cost/limits](../_measurements/S11c_d_near_unity_uniform_continue_result.md).
- **CONDITIONAL:** finite evidence includes a **38**-point inner bank and two local Gaussian actions, not a complete pressure/packet or light-loss action. The local build retained literal **Claude NEEDS / Grok CLEAR**, followed by a local tooling repair, not fresh clearance. [Defect status, “fixed inner-kernel bank accepted” and “two local Gaussian actions accepted”](../_measurements/S11c_d_defect_near_unity_applicability_status.md).
- **CONDITIONAL:** later first-order sources distinguish the two polarizations; finite matched amplitudes are not attenuation. The handoff records positive `G0/Gref`, `K1=0` and scoped C3 pole exclusion, not full work/field/leakage. [Source result](../_measurements/S11c_d_first_order_source_result.txt), [matching result](../_measurements/S11c_d_first_order_transverse_matching_result.txt), [handoff, “What is established”](../../../docs/light_em_investigation_handoff.md). Claude's literal source/build verdicts are **CLEAR FOR THIS NATIVE TRANSVERSE END-FLUX BUILD** and **CLEAR FOR THIS REAL-AXIS RECEIVING AND CROSS-FLUX BUILD**, with limited inspected evidence, not independent result clearance. [Flux review, opening/limits](../_measurements/S11c_d_first_order_transverse_flux_build_claude.md), [receiving review, opening](../_measurements/S11c_d_first_order_receiving_regular_build_r3_claude.md).
- **OPEN — S11c-d/e and Q2/S22:** physical work/support, missing orders and full loss remain; numerical recovery and centre implementation are parked. [Handoff, “Work state and next action”](../../../docs/light_em_investigation_handoff.md).

### 06 — Light no-leak clean-condition packet

**Question:** which symmetry/empty-channel assumptions could protect light on a nonuniform background? [Clean-condition packet, Hypothesis H and falsifiers](../directives/S11c_d_clean_condition.md).

- **CONDITIONAL:** round 5 reports a second single-engine measurement of H-planar on the SymPy operator: **`MIXED 0`** and **`WRONG_SIDE 0`** in all four cases, with responsive controls. This is review-leg evidence, not a blind Wolfram/Part 1 build. [Disposition, “Round 5”, measurement paragraph](../directives/_measurements/S11c_d_clean_condition_review_disposition.md#round-5--review-of-v5-codex-authored).
- **UNRESOLVED — packet review:** literal Claude **“not clear.”**, Grok **“not cleared.”**; selection-rule support did not clear Directive B or broader no-leak claims. The round-5 disposition says the empty odd-channel premise is rest-frame/drain-frozen, and discusses a minimal drain model and coherent trapped-mode deformation. Those scoped review arguments are not a computed native draining-throat result. [Disposition, round-5 verdicts and R5-1/R5-2](../directives/_measurements/S11c_d_clean_condition_review_disposition.md).
- **OPEN — S12 and Q2/S22:** the possible live conversion/drain channel and actual trapped support need their own specified dynamics. B was paused at `af9d4d55`; do not convert its conditional zero into general immunity. Sources: [disposition, final “Process decision”](../directives/_measurements/S11c_d_clean_condition_review_disposition.md), [plan, S12 and Q2](../V3_STEP_PLAN.md). The clean-condition A/B text, literal reports and cited probe outputs are the evidence set to consider, not A alone.

### 07 — Lean S9–S11, formalization policy and CAS bridge

**Question:** what is proved mathematically, rather than sampled or assumed? [Formalization policy, “Division of responsibility”](../lean/FORMALIZATION_POLICY.md).

- **CONDITIONAL:** S10 proves the supplied constant-coefficient curl action's full plane-wave classification for positive density/stiffness and nonzero wavevector, including **D−1** transverse dimensions and the D=1 exception. The six-family contract received Claude **CLEAR** and Grok **CLEAR**; this does not clear all CAS production/export claims. [S10 result, “Supplied action and checked conclusion”](../lean/s10/RESULT.md), [fidelity review, opening/review table](../lean/s10/FIDELITY_REVIEW.md).
- **CONDITIONAL:** the anisotropic audit found a missed perpendicular transverse stratum; focused SymPy/Wolfram checks cover D=3/4 witnesses, not a global sample-based classification. [Anisotropic report, “Finding and repair” and “Proof and integration boundaries”](../_measurements/S10_anisotropic_strata_report.md).
- **CONDITIONAL:** S11's homogeneous, invariant/bulk, variable-coefficient, pole and scattering-algebra contracts are bounded mathematical results, not the full S11c inverse or loss. Their review scopes are individually indexed. Examples: Claude/Grok **CLEAR** for NP1–NP4, T1–T4 and VC1–VC4. [S11 README, contract sections](../lean/s11/README.md), [pole closure](../lean/s11/POLE_FIDELITY_REVIEW.md), [tail closure](../lean/s11/ANALYTIC_ERROR_FIDELITY_REVIEW.md), [variable-coefficient closure](../lean/s11/VARIABLE_COEFFICIENT_FIDELITY_REVIEW.md).
- **OPEN — Phase 3 Lean pruning; S11c-d applications if reopened:** retain proof/coverage/mutations; evaluate removal of the byte bridge only with the later required build check. No proof is withdrawn here. [Policy, L1–L4](../lean/FORMALIZATION_POLICY.md), [CLAUDE.md, L-LEAN](../../../CLAUDE.md), [directive, Phase 3a](../CLEANUP_2026-10_directive.md).

### 08 — Exploratory throat/EM notes

**Question:** could a physical throat supply compatible support, conversion work and electric-force signs? [Core assessment, opening](../../../docs/oriented_throat_core_response_assessment.md).

- **UNRESOLVED:** Claude's core verdict is **REQUIRES A PHYSICAL CHOICE**. The assessment records a conditional sign obstruction in a restricted conservative comparator; it explicitly denies a universal no-go or solved native force. [Core assessment, opening and “Conditional paper diagnosis”](../../../docs/oriented_throat_core_response_assessment.md).
- **CONDITIONAL:** moving-interface and elastic-reference papers received Claude-only **COHERENT CONDITIONAL BRIDGE** and **COHERENT CONDITIONAL COMPARISON**, with partial coverage/qualifications. They are paper comparisons, not native runtime certificates or adopted constitutive laws. Return already means bulk re-ordering into the brane. [Conversion assessment, opening and qualifications](../../../docs/conversion_work_assessment_and_next_step.md), [reference assessment, opening and rejected review claims](../../../docs/elastic_reference_comparison_assessment.md).
- **OPEN — S12/Q2/Q3/S22; material reference requirements:** formation/reference transport, torque, core work and orientation-to-force mapping remain open. The notes and literal review disagreements remain exploratory and paused; no new throat or leakage claim follows. [Handoff, reading index/work state](../../../docs/light_em_investigation_handoff.md), [native interpretation, §§5.7, 14.3](../../../docs/native_light_em_and_vortex_throat_interpretation.md). Kept claims would need the frozen proposals and qualifications, not a review verdict alone.

### 09 — Muonium gravity corrections

**Question:** what may a later muonium prediction actually inherit? [Correction note, opening and §1](../../../notes/muonium_gravity_corrections.md).

- **WITHDRAWN/FAILED:** absolute/relative lepton throat geometries and geometry-derived gravity claims listed in §2 are do-not-import results. The note explicitly includes **`L/a≈1.85`** as real geometry and reverse-engineered mass/radius constructions. [Corrections, §2](../../../notes/muonium_gravity_corrections.md).
- **CONDITIONAL:** the chosen universal effective branch's far-field match survives only within its stated closure; **κ_ρ=1** is target-matched, not a first-principles species prediction. Passive, active and inertial mass must remain distinct. This is a Codex handoff distilled/reviewed by the orchestrator, described as **research-track + corrections**, not an independent two-engine prediction or literal CLEAR gate. [Corrections, provenance, §§1–2](../../../notes/muonium_gravity_corrections.md).
- **OPEN — S16 and later species/throat work:** geometry-dependency quarantine, the same-order momentum-exchange audit and species-to-worldline matching remain. The documented no-import search was a spot-check, not the full audit. No muonium prediction was supplied. [Corrections, §§3–6](../../../notes/muonium_gravity_corrections.md), [plan, S16](../V3_STEP_PLAN.md). Preserve the correction note and its [full handoff](../../../notes/muonium_gravity_research_track_handoff.md), not just old geometry claims.

### 10 — Execution infrastructure

**Question:** how was desktop/resource risk contained without losing saved work? [Host-freeze report, opening](../_measurements/S11c_d_host_freeze_report.md).

- **UNRESOLVED:** the September 20 whole-host freeze cause remains unconfirmed. The incident was not reproduced; containment cannot guarantee recovery from a hardware/kernel lockup. [Freeze report, opening and timeline](../_measurements/S11c_d_host_freeze_report.md).
- **ESTABLISHED — tooling record:** pooled admission tests and two systemd smoke workers are recorded, with aggregate **16 GiB** reservations, **4 GiB** host reserve, zero swap and disjoint CPU affinity; the report explicitly does not claim cpuset-controller enforcement. These are tooling tests, not scientific review or evidence. [Parallel preparation, tooling and pool paragraphs](../_measurements/S11c_parallel_preparation.md).
- **ESTABLISHED — instructions at the inventory endpoint:** archived AGENTS records no computational deadlines and temporary Claude-only reviews, superseding older paired/deadline instructions (`archive/pre-cleanup-2026-10-04:AGENTS.md`, named September 30/October 4 policies). Phase 3 replaces it with the supplied [current AGENTS](../../../AGENTS.md). Cleanup requires foreground metadata work and **no reviews**, so none is launched. [Cleanup directive, ground rules 1/9 and Phase 3b](../CLEANUP_2026-10_directive.md).
- **OPEN — Phase 3 infrastructure/dependency selection:** preservation dependencies must be traced before pruning. No scientific obligation is discharged by a guard, hash or successful process exit. [Directive, Phase 3a](../CLEANUP_2026-10_directive.md), [benchmark, “Completion and preservation”](../_measurements/S11c_d_numerical_radiating_balance_report.md).

### 11 — Ledger front matter, framing and other notes

**Question:** can a reader tell the current state without reconstructing launch history? [Directive, goal](../CLEANUP_2026-10_directive.md).

- **ESTABLISHED — disposition:** the closeout and STATUS's opening say **PARTIAL**; the plan still includes S11c **NEXT** and inherited historical framing. These are documentation conflicts for Phase 2, not authorization to restart. [Closeout, opening](../steps/S11c_PARTIAL_CLOSEOUT.md), [STATUS, lines 3–24](../../../STATUS.md), [plan, S11b/S11c table](../V3_STEP_PLAN.md).
- **UNRESOLVED:** the older MacCullagh essay describes confinement as derived and the longitudinal mode as controlled; the newer native interpretation keeps rotational reference/stress/material admissibility open. Record the disagreement rather than certify the historical or physics argument here. [Essay, “Four differentiators”](../../../docs/s11_maccullagh_differentiation.md), [native interpretation, §§5.7, 14.3](../../../docs/native_light_em_and_vortex_throat_interpretation.md), [closeout, MacCullagh section](../steps/S11c_PARTIAL_CLOSEOUT.md).
- **OPEN — deferred S11c-b/c1/c2 work and Phase 2 documentation:** remote-compute documents are planning requirements, not completed heavy comparisons; their hardware/cost estimates are historical. No provider research or deployment is performed. [Remote requirements, status and §§1–2](../../../notes/remote_compute_requirements.md), [deferred runs, S11c entries](../DEFERRED_HEAVY_RUNS.md). No independent physics verdict is claimed for these front-matter/planning notes.

## C. Conflicts and questions — left for the reviewer

These are the original Phase 1 **UNRESOLVED** questions, preserved as the review input. The [Phase 1 review, C1–C8](PHASE1_REVIEW.md) now supplies the editorial answers applied in Phase 2; it does not resolve the underlying open physics. This inventory did not adjudicate them independently ([ground rule 1](../CLEANUP_2026-10_directive.md)).

1. **Current status versus historical queue.** STATUS's PARTIAL opening coexists with “BUILD IN FLIGHT”, all-four-repairs “UNREVIEWED”, and instructions to finish the engine. The plan still says S11c NEXT. Which exact clauses should the Phase 2 rewrite replace? [STATUS, lines 3–24](../../../STATUS.md); [plan, S11b/S11c table](../V3_STEP_PLAN.md); [joint repair disposition, opening](../_measurements/S11c_upstream_repair_review_joint_disposition.md).
2. **Old b/c2 closes versus changed operators.** b's 2026-09-03 sign/stale-output narrative and c2's 2026-09-09 “VALUES are unaffected”/SOUND narrative predate the changed b/c2 exports. How should their historical verification be delimited against current repaired sources and uncleared composition? [b, status box](../steps/S11c_b_variable_coefficient_operator.md); [c2, status box/arc](../steps/S11c_c2_self_energy_fold.md); [trace repair, “Full c2 regeneration”](../_measurements/S11c_c2_trace_repair_report.md); [joint disposition](../_measurements/S11c_upstream_repair_review_joint_disposition.md).
3. **Mixed-term wording is stronger in the closeout.** Its prohibition says the old `[0,2]` slot “omits a direct mixed height–slope contribution.” The focused tanh report proves a bare upper-face source action but expressly does **not** by itself establish a missing term in the complete closed response. Later corrected constructions are recorded, without a corrected production scattering result. Which later evidence supports precisely which stronger statement? Do not resolve this by promoting the focused diagnostic or declaring it disproved. [Closeout, prohibited claims](../steps/S11c_PARTIAL_CLOSEOUT.md); [tanh result, opening/review limits](../_measurements/S11c_upstream_mixed_tanh_result.md); [defect applicability status, opening status table](../_measurements/S11c_d_defect_near_unity_applicability_status.md).
4. **No-leak scope versus live conversion.** Clean-condition round 5 limits its zero to drain-frozen empty-channel premises and raises an unanswered distributed-return question. S12 uses dynamical order conversion, not the old mass-sink picture. What can transfer from the minimal review model to that native law? No such transfer is accepted here. [Clean-condition disposition, R5-1 and final question](../directives/_measurements/S11c_d_clean_condition_review_disposition.md); [plan, S_leak and S12](../V3_STEP_PLAN.md).
5. **MacCullagh framing versus material debt.** The essay's title/“derived substrate” language and the native interpretation's open rotational-reference/torque questions need a common claim boundary. Neither literature history nor material admissibility is reassessed here. [Essay, differentiators 1–4](../../../docs/s11_maccullagh_differentiation.md); [native interpretation, §§5.7, 14.3](../../../docs/native_light_em_and_vortex_throat_interpretation.md).
6. **Literal reviewer claims versus their assessments.** The elastic-reference assessment rejects the report's conversion-funded energy, automatic load, transport and torque/sign conclusions; the core assessment similarly flags source/sign mismatches. These remain disagreements between preserved records, not newly established resolutions. [Elastic assessment, “Which review claims are not accepted”](../../../docs/elastic_reference_comparison_assessment.md) and its linked literal report; [core assessment, “Corrections and limits”](../../../docs/oriented_throat_core_response_assessment.md) and its linked literal report.
7. **Paper omission and hidden verification.** `paper/parts/part01_light.tex` inputs S9–S11b only; there is no S11c input there. `macros.tex` defaults `showstageverificationfalse`. Where should a visible PARTIAL account and its qualifications go? No PDF was built in Phase 1. [Part I, entire file](../paper/parts/part01_light.tex); [macros, lines 16–26](../paper/macros.tex); [directive, Phase 2 Paper](../CLEANUP_2026-10_directive.md).
8. **Completed formal contracts versus stale handoffs.** The pole-repair reception says Lean fidelity review was pending, whereas NP1–NP4's later closure reports both CLEAR. S11's README still heads its scattering observables “in progress” while specific bounded contracts say complete. Which application debts remain, distinct from completed theorem review? [Pole repair report, opening](../_measurements/S11c_d_nonlinear_pole_repair_report.md); [pole closure, opening](../lean/s11/POLE_FIDELITY_REVIEW.md); [S11 README, scattering observables](../lean/s11/README.md); [bookkeeping closure, opening/final paragraph](../lean/s11/SCATTERING_BOOKKEEPING_FIDELITY_REVIEW.md).

## D. Phase 1 proposals — approved with amendments in the review

The original proposals below were approved subject to [PHASE1_REVIEW.md, Part D](PHASE1_REVIEW.md), which governs the Phase 2 edits. They authorize no calculation or claim upgrade; the amended directions supersede the original wording where they differ.

1. **Create one canonical S11c-d step record** (proposed name `steps/S11c_d_profile_conditioned_scattering.md`). Fold B01–B05's spec/amendment hierarchy, omega=1 finite scope, exact omega=3 table, selected uniform match and first-order limits into it. Link evidence rather than copy the checkpoint chronology; label d PARTIAL. Evidence: the reports cited in B01–B05 and the [closeout](../steps/S11c_PARTIAL_CLOSEOUT.md).
2. **Reconcile `steps/S11c_PARTIAL_CLOSEOUT.md` with that record.** Keep the short programmer-facing summary and explicit owners. Resolve C3's wording only after reviewer direction; retain “no accepted matched-speed defect loss” and the empirical-floor qualification. [Closeout](../steps/S11c_PARTIAL_CLOSEOUT.md); [benchmark, controls](../_measurements/S11c_d_numerical_radiating_balance_report.md).
3. **Rewrite affected parts of the existing b and c2 step records.** Add the four repairs, separate the Wolfram N6 trace repair and d sheet repair, and attach their actual scoped review statuses. Preserve F/G as withdrawn interpretations and cross-engine composition as open; do not infer a corrected production loss result. [b status](../steps/S11c_b_variable_coefficient_operator.md); [c2 status](../steps/S11c_c2_self_energy_fold.md); [joint disposition](../_measurements/S11c_upstream_repair_review_joint_disposition.md). No new per-repair ledger cards.
4. **Rewrite `V3_STEP_PLAN.md`'s S11/S11c sections and `STATUS.md` as a short front door.** Replace superseded queue prose, link the canonical records and name S12 next. Preserve the plan's real dependencies: S11c-e's weak-limit connection, charge-phase motivation, and S22 items 1/4; do not make all S12 or Q-sector work depend on a loss number. [Plan, S12, charge-phase preamble and S22](../V3_STEP_PLAN.md); [d spec, §3d](../directives/S11c_d_SHARED_PHYSICS.md); [closeout, downstream uses](../steps/S11c_PARTIAL_CLOSEOUT.md).
5. **Bring the existing S10/S11 records into agreement with bounded formal closure.** Fold concise coverage/review links into the canonical records; distinguish exceptional-stratum repair, general proof and sampled CAS evidence. Propose a small replacement of the stale S11 Lean README scattering-status prose, not new theorem work or new cards. [S10 record, opening and Lean addenda](../steps/S10_two_transverse_photons.md); [S11 record](../steps/S11_stray_longitudinal.md); [B07 sources](#07--lean-s9s11-formalization-policy-and-cas-bridge).
6. **Add a visible S11c PARTIAL account to `paper/parts/part01_light.tex`**, consistent with the approved records; keep reader-critical limits outside the suppressed Verification field. Build/check the PDF in Phase 2 only. [Part I](../paper/parts/part01_light.tex); [macros, lines 16–26](../paper/macros.tex); [directive, Phase 2](../CLEANUP_2026-10_directive.md).
7. **Keep exploratory/throat and muonium ownership links in the existing closeout/plan.** Preserve frozen proposals and literal reviews. Propose a short qualification of the older MacCullagh framing essay after C5 is reviewed, without rewriting history as a solved material problem. No new S22 theory or gravity prediction. [B08–B09 sources](#08--exploratory-throatem-notes); [essay](../../../docs/s11_maccullagh_differentiation.md); [plan, S16/S22](../V3_STEP_PLAN.md).

## E. Keep rules — replaced by the approved tier policy

The [Phase 1 review, verdict](PHASE1_REVIEW.md) supersedes the original Part E with [directive Phase 3a](../CLEANUP_2026-10_directive.md). No files are selected or pruned in Phase 2.

- **Tier 1:** established/conditional claims that later work may use retain their report, producer, inputs, outputs and transitive dependencies, enough to rerun from the kept tree.
- **Tier 2:** unresolved, withdrawn, paused or exploratory work retains the final report/disposition and any literal reviews cited by the ledger. Its transitive chain stays only in `archive/pre-cleanup-2026-10-04`; its hash pins need not resolve in the pruned tree. **When unsure, use Tier 2.**
- **Process records:** prune unless a Tier 1 script reads them at runtime. An incidental hash link is not a reason to keep a process chain.
- **Directives:** keep the cleared governing spec/build directive, not drafts, superseded versions or per-round prompts.
- Scratch citations and Lean build dependencies follow the directive's specific Phase 3 rules. No blanket dependency closure applies to Tier 2.

The original Phase 1 proposal remains at `174f4fa9`. The [review, C1–C8 and Part D](PHASE1_REVIEW.md) supplies the editorial decisions for Phase 2; Part B above remains the detailed reading guide.

## F. STATUS debt carry-over

The short STATUS now keeps the unfinished work visible without restarting it. This sweep covers the entire
1,289-line `c8885c4f:STATUS.md`, including its superseded S11 and v2 sections. Repeated mentions of one
obligation are grouped into one row; counts below are **groups of obligations, not files or physics results**.
The old line numbers always refer to that revision (`git show c8885c4f:STATUS.md`), not today's STATUS.
This is an **ESTABLISHED document inventory**, implementing [Phase 2 review, F6](PHASE2_REVIEW.md#fixes).

Each row has exactly one editorial disposition. **carried** means OPEN under the owner or register linked
from today's [STATUS, open questions and carry-overs](../../../STATUS.md); it is not scheduled work.
**resolved** closes only the named historical task on the cited record's stated scope.
**superseded** retires a queue, procedure or framing in favor of a named replacement; it does not prove its
underlying physics. Conflicting closure claims are carried, not adjudicated in this cleanup.

| Item | Old STATUS line(s) | Disposition | Evidence / current owner |
| --- | --- | --- | --- |
| Held-defect omega=3 loss number | 5 | carried | S11c-d if reopened. [Benchmark, controls][benchmark]: unresolved within an empirical floor, not a bound. |
| Full d engine, response/FORM integration, blind WL/T7, final export and full physical work | 12, 17, 21–24, 31–32 | carried | S11c-d/e if reopened. [d record, intermediate results/ownership][d]; [builder, retained solver/export contract][builder]. Partial constructions do not close this list. |
| “BUILD IN FLIGHT,” resume-until-complete, stabilization-first and old d review/build sequences | 7, 11–12, 17, 19–24, 26, 32, 39 | superseded | [PARTIAL closeout, opening/next step][closeout] and [cleanup directive, Phase 2 STOP][directive]. S12 is next; no automatic S11c continuation. |
| Blanket “four upstream repairs UNREVIEWED” status | 13, 17 | superseded | [Joint disposition, scoped table and remaining issue][joint]; [b/c2 Changed after close][b]. Scoped reviews exist; full composition remains open, and the direct `[0,2]` entry is unjustified by order counting. |
| d sheet repair's missing independent review | 13, 17 | carried | S11c-d if reused. [Sheet report, scope][sheet]; [d record, unreviewed repairs][d]. The separate Wolfram N6 trace repair is likewise unreviewed (F1), not silently included in the four reviewed repairs. |
| μ_S shear-normalization debt | 17 | carried | S11b/S11c-b if reused. [Inertia report, final paragraph][inertia]; [end-resolvent report, lines 177–178][end]. It remains distinct from the inertia-sign repair. |
| ≥64 GB b residual/all-case WL output, c1 giants/full residual, full c2 self-energy comparison; b #88 hardenings/re-adjudication and #90 sign/survivor flags | 13, 17, 24, 32, 37, 39, 121, 124, 139, 142–160, 167, 173–175, 183–198, 205 | carried | S11c-b/c1/c2/e if reused. [DEFERRED_HEAVY_RUNS, cross-engine residual and per-engine entries][heavy]; [b, established versus owed][b]; [c1, status of close][c1]. Includes term origins, P2a/P2b, skipped in-band controls, traction/DtN/energy/flat-leg and density caveats; no all-family agreement is inferred. |
| b namespace/structure gaps, v2 N1/N2/N3 notes and absent emitted admissibility discriminator | 239 | carried | S11c-b/comparator maintenance. `c8885c4f:STATUS.md:239` lists 57 structure and 12 namespace gaps and §5 coverage; [b, owed cross-engine checks][b] is narrower than closure of them. No discharge found in this sweep. |
| b kernel-build optimization | 196–197 | carried | S11c-b tooling if reopened; explicitly deferred at `c8885c4f:STATUS.md:197`; [heavy, PY #89][heavy] retains the costly control path. No optimization is scheduled. |
| c2 carrier/source/Φ operand corroboration and inherited sign/density/traction/DtN/energy questions | 31, 37, 39, 42–43, 47, 106, 109 | carried | S11c-c2/d if reused. [c2, close status and cross-engine debt][c2]; [N6 disposition, Path B][n6]. Covariance or matched zeros do not settle operand agreement. |
| F/G numerical re-grounding and dependent §5e/§3c interpretation/wording | 31–32, 37–38, 47, 63, 90, 97, 102, 106, 109 | carried | S11c-c2 if needed; [F/G deferral, disposition][fg] and [c2, Changed after close][c2]. Indefinitely paused, nonblocking, interpretations withdrawn; the old “values unaffected” sentence is pre-repair only. |
| c2 publication repair's owed complete fresh-Claude leg | 101, 105–109 | carried | S11c-c2 record maintenance. **Conflicting closure statements**, retained below: [publication adjudication, verdict][export-review] still says owed; old STATUS line 101 says done. No new clearance asserted. |
| c2 N6 wrong-anchoring test/spec correction and later covariance construction | 84–103, 105–109 | resolved | [c2, review/repair arc][c2] and [N6 disposition, earned versus not earned][n6] record the corrected test and scoped covariance. This does not resolve the separate operand debt. |
| N6 premise caveats: physical correctness of Φ, V transform and extracted-block leakage | 89 | carried | S11c-c2/d if reused. [c2, review caveats (lines 180–184, 203 at Phase 2)][c2]; [N6 disposition, §7][n6]. Correcting the test and proving scoped covariance did not settle these premises. |
| b carrier-harness re-author/relaunch/commit queue | 59, 65–75 | superseded | Old STATUS line 59 records b closed after the harness work; [b record, original verification/close][b] is the replacement. The former one-full-leg/partial-leg attempt remains history, not relabelled paired clearance. |
| N6 harness output/budget repairs and stale production-output regeneration | 49–63, 69–82 | resolved | Old STATUS lines 50, 55–56 record the completed replacement; [N6 harness review, final verdict][n6-harness]; [c2, Wolfram engine and comparator provenance][c2]. Scoped historical task only. |
| c2 T7 comparator, reconcile disposition and step-record writeup | 34–57, 69, 75, 90, 97 | resolved | [N6 disposition][n6] and [c2 record][c2] complete the Path-B writeup; unresolved operand families were carried, not declared agreement. |
| c1 retro-review record corrections (density, grazing, independence and carry-forwards) | 112, 118 | resolved | [c1 corrections adjudication, corrected records and round 2][c1-corrections]. Record-only corrections; unchanged engines are not blanket closure of c1's remaining comparisons. |
| c1 split/build/export-migration/comparator queue | 120–139 | superseded | [c1, status of close and engines][c1] records the later build and scoped comparison. Remaining ≥64 GB work is carried above. |
| Five off-path export-chain pairs not checked by the new minimal consumer | 130 | carried | Export-chain maintenance; the deferred residual is recorded at `c8885c4f:STATUS.md:130`, citing `_measurements/f9c_pair_scan.py`. No new dependency closure or rerun is claimed. |
| S11c-a full-axis control KEYING and abstract-symbol control-form characterization | 239, 241 | carried | S11c-a/b; [a record, standing limits/owed items][a]. Missing density/face axes and other control bookkeeping remain separate from its physics-family comparison. |
| S11c-a FACE_SHIFT/current-freeze/projection comparator and representational adjudication queue | 241, 268–270, 293–305, 345–348, 377–396, 450 | resolved | [a record, fixes and cross-engine reconciliation][a], including `cccb4f9e` and close `3b552426`. The blanket-collapse “all agree” predecessor is withdrawn history, not the evidence. |
| b admissibility/kinetic/advective/coupling verdicts and basis/jet-depth repair queues (#84–90) | 162–239 | superseded | [b, walk and established versus owed][b] and its Changed after close section replace the interim diagnoses. Old “kinetic gap”/on-box feasibility claims are not current; composition/control debts remain above. |
| d profile-class/regime choice, governing spec and Fourier-reduction build plan | 27–32, 39 | resolved | [Shared physics, §§0–3][d-spec], [builder, retained contract][builder]. The design was supplied; completed scattering is not implied. |
| S11c-e conversion/strong-edge/bench-optics work and eventual family roll-up | 17, 24, 32, 39, 167, 398–400, 601–609 | carried | S11c-e; [closeout, actual downstream uses][closeout] and [a, standing limits][a]. A PARTIAL note is not completion of the full family. |
| Background-flow correction, frozen-wall and nonlinear-radiation limits | 485, 589–591, 601–608, 617 | carried | S12 for flow; [a, standing limits][a], [S11b, scope/limits][s11b], [plan, S12][plan]. The uniform result does not remove these qualifications. |
| S11b coverage/keying, exceptional denominator domains and non-passive power-source obligation | 497, 583–591 | carried | S11b/interface-model owner; [S11b, permeable response and cross-engine comparison][s11b]. Different solved/undecided sets and regular-domain agreement are not full coverage; a non-passive model needs a named power source. |
| S11b X-1 basis correction, unified step/card and export-integrity queue | 466–471, 537–538, 594–599 | resolved | [S11b, two engines/comparison][s11b] records `53fcd98d`, `8ddccb74`, `565b3fe8`. This is the uniform step's recorded closure. |
| Uniform grazing/interface question formerly assigned from S11 to S11b | 484–485, 617 | superseded | [S11, ownership boundary after grazing audit][s11] and [S11b, scope/comparison][s11b] now distinguish an allowed channel from its actual use; nonuniform loss is still carried. |
| S11 engine defects, undecided/omitted/spurious strata, witness gaps and obligation-4 coverage debt | 615–630, 726–730, 757–761, 783–785 | carried | S11 maintenance only if a downstream locus requires it; [DEFECT_REGISTER, WL entries][defects], [S11, strata-audit closure][s11]. Certified instruments and bounded Lean results do not erase CAS production gaps. |
| Frozen T7 comparator, complex-k rank-drop mismatch and native-boolean rejection | 627, 630, 653–655, 787–791 | carried | S11 comparator maintenance; [REBUILD_HANDOFF, open items/physics filter][rebuild], [DEFECT_REGISTER, F7 and WL entry 5][defects]. No comparison is run by cleanup. |
| S11 PY/WL sweep, D2/D4 memory-wall, spec/engine rebuild and census-instrument launch sequences | 621–625, 630–634, 640–724, 766–791, 807–845 | superseded | [S11, verification and strata-audit closure][s11] records the later completed campaigns and their limits. The historical compute-box fallback is not a new request; remaining defects are carried above. |
| S10 non-author fidelity-review confirmation | 14 | resolved | [S10 fidelity review, fixed revision and independent reviews][s10-fidelity]: Codex author, fresh non-author Claude and Grok sessions. Reviews inspected source/diagnostics, not a fresh Lean/CAS/generator rerun. |
| CLAUDE.md L control consistency review | 14, 17 | carried | Policy maintainers; [Phase 2 review, F6][p2-review]. Landed unreviewed in mixed-scope `c2bdb663`; cleanup edits neither CLAUDE.md nor its rules. |
| lean/FORMALIZATION_POLICY.md consistency review | 14, 17 | carried | Policy maintainers; [Phase 2 review, F6][p2-review]. Contract review is not review of this policy's adoption. |
| S10 production/export/ledger work beyond Lean, including three unlegged record corrections | 17, 653–655, 844 | carried | S10 maintenance; [REBUILD_HANDOFF, “What S10 still owes”][rebuild], [S10 fidelity review, stopping point][s10-fidelity]. No blanket discharge by Lean. |
| S10 D12/naming, F3/F4 row/regeneration, twelve dimension-key refs and split spectrum vocabulary | 653–655, 819–822, 844 | carried | S10/export maintenance; [DEFECT_REGISTER, C19/C20][defects], [REBUILD_HANDOFF, four-item pass and physics filter][rebuild]. Conflicting “owed”/“ceremony” descriptions are retained below; no export edited. |
| S9/S10 requirements-register pass 1 | 653–655, 844 | resolved | [SUBSTRATE_REQUIREMENTS, opening][substrate] explicitly records pass 1 complete. Later passes and the unbuilt substrate remain open; this closes the old capture task only. |
| Comprehensive substrate requirements, S1–S8 delivery and curl-only material admissibility | 609, 627, 855–859 | carried | S1–S8/material audit; [SUBSTRATE_REQUIREMENTS, entries and remaining passes][substrate]. No substrate execution is inferred from the requirements list. |
| Muonium dependency quarantine, PN geometry-independence and structural audit; species matching | 16–17, 883–885 | carried | S16 and later species/throat work; [muonium corrections, §§2–5][muonium], [DEFECT_REGISTER, D1][defects]. Numerical geometry is withdrawn; the audit is not done by this note. |
| Remaining force sectors and real throat/holder/charge mechanism | 609–610, 627, 879–885 | carried | Q2/Q3/S22 and the plan's named sectors; [plan][plan], [closeout, paused exploratory notes][closeout], [DEFECT_REGISTER, A/C/D][defects]. No throat mechanism is adopted. |
| Remote compute implementation, licensing/credit questions, smoke tests and acceptance | 13, 17, 717 | carried | Infrastructure owner; [remote requirements, §§8–10][remote-req] and [plan, §§7–9][remote-plan]. The former S11 D2 fallback is superseded above; the ≥64 GB S11c jobs remain deferred. Historical prices/configurations are not revalidated. |
| Commit-scope/annex review and local history rewrite (including c2bdb663 bundle) | 12, 15, 17 | carried | Cleanup Phases 3–4; [directive][directive]. Nothing rewritten or pruned in this fix pass. |
| Stale unpushed-commit counts and push/resume commands | 12, 17, 24, 27, 32, 39, 43, 47, 69, 80 | superseded | [Directive, archive and no-push rule][directive]; [Phase 2 review, F6][p2-review]. Current STATUS separates the reported work-branch push date from the already-pushed archive tag; no push performed here. |
| MEMORY.md compaction (retain live pointers) | 17, 24, 32, 63, 69, 75, 82, 97, 103, 109, 115 | carried | Documentation maintainer; old STATUS at these lines is the record. Historical size estimates are not a current measurement. This is maintenance, not a physics gate. |
| Temporary-only outputs/loaders and review/scratch evidence portability | 82, 90, 103, 109, 115, 627, 630, 795–801 | carried | Cleanup Phase 3a dependency selection; [directive, Tier 1/scratch rules][directive]. Keep-by-tier replaces blanket retention; no file moved or restored here. |
| Reverted watchdog/agent-roles/process detour and old resource/review schedules | 68, 139, 151, 762–764, 795–801, 1225–1243 | superseded | Old line 68 records `e406cde6`; [AGENTS, standing execution/review policies][agents] and [directive, ground rules][directive] govern now. No watcher, deadline or review launched. |
| v2 walkthrough/dimension-rewrite as the active front, and the old whole-ledger phase schedule | 834–846, 872, 895–917, 1001–1028, 1201–1212, 1250–1257 | superseded | [V3_STEP_PLAN, S11c close/S12][plan] and [PARTIAL closeout][closeout]. Distinct unresolved v2 contents are carried below; the old queue is not revived. |
| The a-pin/throat-radius identification and its repair-before-walkthrough prerequisite | 874–876, 930–938, 1010–1012 | resolved | [DEFECT_REGISTER, A1][defects] records removal `407eed94`; no other length identification is thereby settled. |
| v2 K_eta/T_Omega/mu_eta reduction-level adjudications and invisible T_w/Tw naming debt | 925–944 | carried | Dormant v2 dimension maintenance; [DIMENSION_REWRITE, §§7–8 and open routes][dim]. Same dimensions do not identify quantities. |
| v2 missing per-quantity route tables and broader inventory/independent-derivation limits | 948–980, 1050–1056 | carried | Dormant v2 physics-verification work; [DIMENSION_REWRITE, §4-c1 and §§8–12][dim]. “Converted” is not “verified”; provisional reachability counts do not prove completeness. |
| v2 stage027 shape/computed-vector reachability, then 027/021 conversion and harness hazards | 1001–1004, 1054–1067 | carried | Dormant v2 dimension maintenance; [DIMENSION_REWRITE, §§8–9][dim]. Includes the generator's axis-order failure before 035/036/037; no conversion scheduled. |
| v2 small ablation driver and shared Wolfram DIM emitter | 1013–1028 | carried | Dormant v2 tooling; [DIMENSION_REWRITE, §12b][dim] and its two requirements documents. The obsolete frozen-fixture machinery stays retired; old STATUS lines 1130–1131 also retain the `run_all_audits.sh` exit-zero limitation, not validation evidence. |
| Seven stage023 derivation routes | 1030–1045 | carried | Dormant v2: [DIMENSION_REWRITE, §12 WORK-023-*][dim] — MOMENT-CONVENTION, STAGE009-MOMENT0, D0-SEAM, STIFFNESS-REDUCTION, L1-L2-PROFILE-IDENTITY, CS-EVALUATION, SOURCED-PROVENANCE. W3/q_free are not added to this list. |
| v2 035/036 independent routes and old waiver/impossibility claims | 1088–1092, 1133–1135 | carried | Dormant v2 dimension maintenance; [DIMENSION_REWRITE, §3b][dim]. 037's prototype is not completion of the remaining routes. |
| v2 four surviving corpus findings (tier contradiction, false consumption, wrong loci, missing citations) | 1137–1143 | carried | Dormant v2 documentation; [v2 handoff, OPEN CORPUS FINDINGS][v2-handoff], [DEFECT_REGISTER, F][defects]. Retiring the census front did not discharge these findings. |
| v2 Part VII constant-provenance firewall and unbuilt stage046/input-source claim | 1144–1164 | carried | Later integration/v2 maintenance; [Part VII split, stage046 row][part7], [DIMENSION_REWRITE, charter/open work][dim]. Claimed consumption is not an implemented dependency. |
| v2 stage044-v2 dynamic sleeve and stage045 drain/return model choice | 1174–1183 | carried | Physical successors S7/S12/S22; [v3 plan][plan]. Original preparation remains [044][v2-044] / [045][v2-045]; no old implementation or pending material choice is silently approved. |
| v2 manifest semantic core, archive/per-file split and stage045–049 completion | 1184–1199, 1212 | carried | Dormant v2/integration maintenance; [v2 handoff][v2-handoff], [Part VII split][part7]. The old “archived”/“not yet executed” conflict remains below. |
| v2 stage043 exact-count/category discrepancy | 1214–1221 | carried | Dormant v2 record maintenance; [stage043 note, count discussion][v2-043] versus the script loci listed at `c8885c4f:STATUS.md:1215–1220`. Do not quote the note's range as reconciled. |
| Old EM boundary operator's 144/144 unresolved record | 1268–1269 | carried | Q2/Q3/S22 boundary-model history; [EM handoff, boundary-operator verdict][em-handoff]. Exploratory newer notes do not turn it into a solved native mechanism. |
| χ_Q / 54/5 reconciliation between pathA_22b and moving-throat ladder | 1278–1289 | carried | Later integration/Part VII; old STATUS lines 1286–1289 explicitly retain the mismatch. [Moving-throat ladder, Gate 4][ladder] and [Part VII split][part7] are context, not a new reconciliation. |
| Retired v2 fingerprint-join, sidecar-custody and module-digest-as-proof procedures | 863–870, 1005–1007, 1078–1086, 1094–1129, 1165–1172 | superseded | [DIMENSION_REWRITE, charter/§4/§9][dim] and [development pipeline, posture][pipeline]. Their actual drift/coverage limits remain; no new trust framework is proposed. |

**Disposition counts:** 42 carried · 10 resolved · 11 superseded = 63 grouped items.

### Carry-over conflicts and questions

These are **UNRESOLVED record conflicts**, not new physics findings or requests to restart work:

1. **c2 export-review closure:** `c8885c4f:STATUS.md:101` says the fresh-Claude leg was done, while the
   [publication adjudication, verdict][export-review] explicitly leaves it owed. Carry it until the reviewer
   identifies the closing report; this pass does not choose between the records.
2. **S10 export housekeeping:** [DEFECT_REGISTER, C20][defects] retains a four-item repair pass, while
   [REBUILD_HANDOFF, physics filter][rebuild] calls F3/F4 regeneration ceremony and separately keeps unlegged
   record claims. The short STATUS preserves the owner/pointers, not an instruction to execute either list.
3. **v2 manifest archive:** old STATUS lines 1193–1199 call it both archived and not yet executed; the old census
   paragraph at 1137–1143 has similar wording. Its surviving content is carried. A per-file decision belongs to
   the named dormant workstream, not this S11c cleanup or an edit to the off-limits v2 tree.

**Scope of this fix:** documentation and Git metadata only. No scientific functions, reviews, annex fetches,
pruning or history rewrite. Phases 3–5 remain pending under the [directive][directive].

[benchmark]: ../_measurements/S11c_d_numerical_radiating_balance_report.md
[d]: ../steps/S11c_d_profile_conditioned_scattering.md
[builder]: ../_measurements/S11c_d_sympy_builder_report.md
[closeout]: ../steps/S11c_PARTIAL_CLOSEOUT.md
[directive]: ../CLEANUP_2026-10_directive.md
[joint]: ../_measurements/S11c_upstream_repair_review_joint_disposition.md
[b]: ../steps/S11c_b_variable_coefficient_operator.md
[sheet]: ../_measurements/S11c_d_sheet_repair_report.md
[inertia]: ../_measurements/S11c_inertia_repair_report.md
[end]: ../_measurements/S11c_d_end_resolvent_report.md
[heavy]: ../DEFERRED_HEAVY_RUNS.md
[c1]: ../steps/S11c_c1_curved_bulk_closure.md
[c2]: ../steps/S11c_c2_self_energy_fold.md
[n6]: ../_measurements/S11c_c2_N6_reconcile_disposition.md
[fg]: ../_measurements/S11c_c2_FG_regrounding_deferred.md
[export-review]: ../_measurements/S11c_c2_export_repair_rereview_adjudication.md
[n6-harness]: ../_measurements/S11c_c2_N6_harness_build_review.md
[c1-corrections]: ../directives/_measurements/S11c_c1_record_corrections_review_adjudication.md
[a]: ../steps/S11c_a_interface_shape_derivatives.md
[d-spec]: ../directives/S11c_d_SHARED_PHYSICS.md
[s11b]: ../steps/S11b_interface_coupling_law.md
[s11]: ../steps/S11_stray_longitudinal.md
[defects]: ../DEFECT_REGISTER.md
[rebuild]: ../REBUILD_HANDOFF.md
[s10-fidelity]: ../lean/s10/FIDELITY_REVIEW.md
[p2-review]: PHASE2_REVIEW.md
[substrate]: ../SUBSTRATE_REQUIREMENTS.md
[muonium]: ../../../notes/muonium_gravity_corrections.md
[plan]: ../V3_STEP_PLAN.md
[remote-req]: ../../../notes/remote_compute_requirements.md
[remote-plan]: ../../../notes/remote_compute_plan.md
[agents]: ../../../AGENTS.md
[dim]: ../../pde_ledger_v2/manifests/DIMENSION_REWRITE.md
[v2-handoff]: ../../pde_ledger_v2/_scratch/NEXT_SESSION.md
[part7]: ../../pde_ledger_v2/notes/part7_integration_atomic_split.md
[v2-044]: ../../pde_ledger_v2/notes/stage044_v2_unfreeze_prep.md
[v2-045]: ../../pde_ledger_v2/notes/stage045_nonvariational_block_prep.md
[v2-043]: ../../pde_ledger_v2/notes/stages/ledger_stage043_irreducible_count_range.md
[em-handoff]: ../../../docs/em_analog_next_phase_handoff.md
[ladder]: ../../pde_ledger/notes/stages/moving_throat_pde_completion_ladder.md
[pipeline]: ../../../docs/development_pipeline.md

## G. Phase 3 selection and checks

**ESTABLISHED — repository cleanup, not new physics.** S11c remains PARTIAL. The uniform report and its cited
JSON results are kept, but the tree does **not** contain a ready-to-run copy of its intermediate caches.
Replaying it requires regeneration and export-chain repair. This section incorporates
[PHASE3_REVIEW, H1–H5](PHASE3_REVIEW.md), using H3's default, not option B.

G1/G2 were committed at `0dcc53dd`; the exact supplied AGENTS replacement is at `a811df3b`.
The earlier Phase 3 copies and reconstruction claim at `bd1119bb`/`9b98bd79` are superseded by these fixes.
The 33 ignored process files remain untracked and on disk (`bec92067`); the 20 older tracked ignored files
remain unchanged. No scientific producer, review, cache regeneration, export repair or push ran.

### Selection

[KEEP.tsv](KEEP.tsv) and [PRUNE.tsv](PRUNE.tsv) partition all 4,458 original inventory paths:
**507 kept / 3951 pruned**. PRUNE includes the 33 untrack-only paths.
KEEP has **571 unique paths**, including **64** baseline dependencies, cleanup files and newly
cited output files outside that inventory. The 106 prior copy/container paths were removed; they were
post-inventory additions and are not counted again in the original-inventory PRUNE table.
No relocation-map rows remain.

Same K1–K12 key as Part A. Each entry is **kept/pruned** for the original inventory.

| Workstream | K1 | K2 | K3 | K4 | K5 | K6 | K7 | K8 | K9 | K10 | K11 | K12 | Total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 01 Physics specification and contracts | 4/2 | 0/3 | 1/0 | 4/55 | 0/1 | 0/2 | 0/29 | 0/4 | 0/0 | 0/0 | 0/2 | 0/0 | 9/98 |
| 02 Symbolic build and omega1 premises | 0/1 | 2/222 | 4/144 | 0/70 | 17/78 | 26/173 | 2/447 | 2/345 | 0/0 | 0/0 | 0/36 | 0/0 | 53/1516 |
| 03 Upstream repairs and diagnostics | 0/0 | 0/7 | 9/6 | 3/16 | 4/9 | 3/53 | 2/96 | 7/111 | 0/0 | 0/0 | 0/7 | 0/0 | 28/305 |
| 04 Numerical scattering and omega3 benchmark | 0/2 | 0/3 | 1/5 | 1/4 | 1/3 | 3/41 | 0/67 | 1/22 | 0/0 | 0/0 | 0/6 | 0/0 | 7/153 |
| 05 Near-unity and defect prerequisites | 0/21 | 1/30 | 4/10 | 2/168 | 1/6 | 1/134 | 0/760 | 0/131 | 0/0 | 0/0 | 2/12 | 0/0 | 11/1272 |
| 06 Clean-condition packet | 1/1 | 0/2 | 0/0 | 1/10 | 0/19 | 0/59 | 0/14 | 0/77 | 0/0 | 0/0 | 0/3 | 0/0 | 2/185 |
| 07 Lean and S9-S11 formalization audits | 2/1 | 0/0 | 37/0 | 50/6 | 3/2 | 5/40 | 1/158 | 3/123 | 235/0 | 15/55 | 7/3 | 8/0 | 366/388 |
| 08 Exploratory throat and EM | 0/0 | 0/0 | 0/0 | 7/0 | 0/0 | 0/0 | 0/14 | 0/2 | 0/0 | 0/0 | 6/0 | 0/0 | 13/16 |
| 09 Muonium gravity corrections | 0/0 | 0/0 | 0/0 | 0/0 | 0/0 | 0/0 | 0/0 | 0/0 | 0/0 | 0/0 | 2/0 | 0/0 | 2/0 |
| 10 Execution infrastructure | 0/0 | 0/0 | 1/1 | 0/0 | 0/0 | 2/7 | 0/9 | 0/1 | 0/0 | 0/0 | 1/0 | 3/0 | 7/18 |
| 11 Ledger front matter and other notes | 0/0 | 1/0 | 0/0 | 0/0 | 0/0 | 0/0 | 0/0 | 0/0 | 0/0 | 0/0 | 2/0 | 6/0 | 9/0 |

### H1–H3: outputs kept, caches left on disk

- The eight `.blob` files and the six cache containers named `S11c_uniform_retained_*.out` are removed.
  Their original eight pickles and six cache/output files remain in ignored `_scratch/`; no original was deleted.
  Neither ignore nor annex rules were changed.
- The seventh `.out` was an existing transcript. It is restored at its original path,
  `research/pde_ledger_v3/scripts/out/S11c_d_mixing_scattering_sympy_audit.out`, because the kept c2 trace,
  sheet-repair and end-resolvent reports cite it. Its original annex key is preserved.
- The uniform worker and report remain. The report now directly cites **31 original JSON outputs** in
  [`S11c_d_near_unity_uniform_output`](../_measurements/S11c_d_near_unity_uniform_output/): 24 point files,
  four grazing-limit files, the aggregate checks, denominator certificates and certificate self-checks.
  Names and bytes are unchanged. These are plain JSON, not renamed caches or a replay reconstruction map.
- The rest of `_measurements/retained_uniform/` is removed, including the copied old export producers and
  frozen-method copy. Historical uniform launch/gate/census/input-copy records no longer supply a claimed
  replay bundle; they are in the archive. The report's historical preservation statements describe the original
  run, not the current tracked tree. Historical hashes inside the output JSON retain that same limitation.
- [STATUS, open questions](../../../STATUS.md#open-questions-and-owners) records the replay obligation:
  regenerate intermediate caches from the kept b/c2 engines and repair the S10/S11 export chain before reuse.
  The original caches were never in the archive tag under their scratch paths.

The earlier accepted Tier 2 decisions are unchanged: paused/intermediate and exploratory work retains its
final cited reports with their limits. The clean-condition result remains single-engine and drain-frozen.
The five compact S10 CAS modules (`PY`, `WL`, `Support`, `BasisCompletion`, `Bindings`) stay; removal of the
other 50 modules was accepted. This cleanup does not enlarge Lean's proof scope.

### H4: export drift introduced within the inventoried period

**OPEN — S10/S11 maintenance.** At `dada3b7d`, each export pin matched its producer. The subsequent producer
edits were `56595cf7` (2026-09-11, S10: `TRANSVERSE_RANK_DROP`, argparse and unlink removal) and `035bb654`
(2026-09-16, S11). The exports were not regenerated. This drift arose inside the inventoried period; calling
it “pre-existing” or saying the older producer snapshots repaired it was wrong. [Review, H4](PHASE3_REVIEW.md).

| Export / producer | Export pin, matching the producer at dada3b7d | Current producer SHA-256 |
| --- | --- | --- |
| S10 / S10_brane_mode_spectrum_sympy_audit.py | `7de1764c1c7ff7a87c70074db8d73b8a259aa53294e272e905531c901fd1bb33` | `f609d57afe9b4dedaaf93a5f49b08572779916cf3d14306ae7ea5da60be7f6fe` |
| S11 / S11_stray_longitudinal_sympy_audit.py | `352fc502f6e191afa512a2dec513685a48c9b32bc8b4e63ca84c5af9be629db4` | `dc77f30b3e0b2b72aa13c7531a586e6ef50167949e645d99378de8527e06bfe2` |

`git log dada3b7d..HEAD -- research/pde_ledger_v3/scripts/<producer>` supplies those edit commits. The checker
compares the literal export pins with the current files and the `dada3b7d` Git blobs; it does not import an
export or evaluate scientific expressions. The baseline blobs remain in Git history, not duplicate files.
Regeneration under review or reverting the edits is a future decision. The S10 Lean D3 link may rely on the
S10 edit, so no revert is automatic. STATUS carries this debt separately from the uniform replay debt.

### H5: frozen method

The only hash-pinned frozen-file exemption from the live-link check is
`research/pde_ledger_v3/directives/S11c_d_SCATTERING_FORM_AMENDMENT.md`.
It is restored byte-for-byte from `archive/pre-cleanup-2026-10-04` (SHA-256
`eb6f20a3fd0bce308e53901fa69fdf267b2917e5cc78b319a1e7700074a47560`).
Its old links are interpreted at that archive revision. The checker verifies its exact bytes and skips its
live links. Literal reviewer reports also remain unchanged historical text, as before.

### Phase 5 preservation list — keep on disk

The following scratch run directories are **KEEP ON DISK** in the future Phase 5 list. They contain the
original uniform inputs and intermediate caches; H3 does not authorize their deletion. This is a preservation
mark for Phase 5, not the full Phase 5 disk inventory or permission to clean scratch.

- `_scratch/s11c/s11c-frequency-source-20260919/finish-01` — **KEEP ON DISK**.
- `_scratch/s11c/s11c-thickness-coordinate-20260914/end_source_left` — **KEEP ON DISK**.
- `_scratch/s11c/s11c-thickness-coordinate-20260914/end_source_right` — **KEEP ON DISK**.
- `_scratch/s11c/s11c-thickness-coordinate-20260914/end_pairing_left` — **KEEP ON DISK**.
- `_scratch/s11c/s11c-thickness-coordinate-20260914/end_pairing_right` — **KEEP ON DISK**.
- `_scratch/s11c/s11c-uniform-source-20260919/production` — **KEEP ON DISK**.
- `_scratch/s11c/s11c-uniform-response-20260919/production` — **KEEP ON DISK**.
- `_scratch/s11c/s11c-parallel-near-unity-20261001/uniform-01` — **KEEP ON DISK**.
- `_scratch/s11c/s11c-parallel-near-unity-20261001/uniform-continuation-01` — **KEEP ON DISK**.

The required ignored-file check covers these eight pickles and six cache/output files at their original paths:

```text
_scratch/s11c/s11c-frequency-source-20260919/finish-01/complete/left-frequency-pencil.pickle
_scratch/s11c/s11c-frequency-source-20260919/finish-01/complete/right-frequency-pencil.pickle
_scratch/s11c/s11c-thickness-coordinate-20260914/end_pairing_left/complete.pickle
_scratch/s11c/s11c-thickness-coordinate-20260914/end_source_left/objects.pickle
_scratch/s11c/s11c-thickness-coordinate-20260914/end_source_right/objects.pickle
_scratch/s11c/s11c-uniform-response-20260919/production/complete/uniform-response.pickle
_scratch/s11c/s11c-uniform-source-20260919/production/complete/common.pickle
_scratch/s11c/s11c-uniform-source-20260919/production/complete/uniform-source.pickle
_scratch/s11c/s11c-parallel-near-unity-20261001/uniform-01/complete/operation-index.jsonl
_scratch/s11c/s11c-parallel-near-unity-20261001/uniform-continuation-01/complete/checks.json
_scratch/s11c/s11c-parallel-near-unity-20261001/uniform-continuation-01/complete/operations.sqlite
_scratch/s11c/s11c-thickness-coordinate-20260914/end_pairing_right/complete.pickle
_scratch/s11c/s11c-parallel-near-unity-20261001/uniform-01/complete/checks.json
_scratch/s11c/s11c-parallel-near-unity-20261001/uniform-01/complete/operations.sqlite
```

### Checks after H1–H5

All commands are metadata reads or compilation of existing sources. Logs remain under `/tmp/`; no process
records are added. Exit status is zero unless explicitly noted. Commands run from the repository root unless
another working directory is named.

```sh
python3 research/pde_ledger_v3/cleanup_2026_10/inventory.py --check-selection
```

```text
PASS: inventory partition 507 kept + 3951 pruned = 4458.
PASS: 571 kept paths; 31 byte-identical uniform JSON outputs; no relocation-map rows.
PASS: 0 other Tier1 JSON file/range pins; 26/28 current export pins; 2 compact Lean input hashes.
PASS: 1001 local record links and 279 archive references; 1 exact frozen method exempt from live-link checks.
KNOWN OPEN (H4): 2 export/producer mismatches introduced by 56595cf7 and 035bb654: [('S10_exports.py', 'S10_brane_mode_spectrum_sympy_audit.py'), ('S11_exports.py', 'S11_stray_longitudinal_sympy_audit.py')]
```

The selection check reports H4 as a known open debt, not a successful current-chain check.
`python3 research/pde_ledger_v3/cleanup_2026_10/inventory.py --check-live-exports` prints the same output and
exits **1** for those two mismatches. It does not accept the baseline producer as a substitute current file.
Output JSON is checked for its exact recorded byte digest and JSON syntax, not interpreted as scientific code.
Its historical cache/source receipts are not current-tree dependencies under H3.

```sh
python3 - <<'PY'
import csv,pathlib,py_compile,subprocess
with open('research/pde_ledger_v3/cleanup_2026_10/KEEP.tsv') as f:
    paths={r['path'] for r in csv.DictReader(f,delimiter='\t')}
paths |= set(subprocess.check_output(['git','ls-files','*.py'],text=True).splitlines())
paths=sorted(p for p in paths if p.endswith('.py') and pathlib.Path(p).exists())
for i,p in enumerate(paths):
    py_compile.compile(p,cfile='/tmp/phase3-h-pycache/'+str(i)+'.pyc',doraise=True)
print('PASS: python3 py_compile,',len(paths),'kept scripts; none imported or executed.')
PY
```

```text
PASS: python3 py_compile, 1427 kept scripts; none imported or executed.
```

In `research/pde_ledger_v3/lean/` (the existing S9/S10 default build; no proof edits):

```sh
LAKE_CACHE_DIR=.lake/cache lake build > /tmp/phase3-h-lake.log 2>&1 && tail -1 /tmp/phase3-h-lake.log
```

```text
Build completed successfully (3785 jobs).
```

In `research/pde_ledger_v3/paper/`:

```sh
pdflatex -interaction=nonstopmode -halt-on-error pde_ledger_v3.tex > /tmp/phase3-h-paper1.log 2>&1 &&
pdflatex -interaction=nonstopmode -halt-on-error pde_ledger_v3.tex > /tmp/phase3-h-paper2.log 2>&1 &&
tail -2 /tmp/phase3-h-paper2.log
```

```text
Output written on pde_ledger_v3.pdf (39 pages, 444438 bytes).
Transcript written on pde_ledger_v3.log.
```

The ignored-file command takes the exact 14 original paths listed above:

```sh
python3 - <<'PY'
from pathlib import Path
import subprocess
text=Path('research/pde_ledger_v3/cleanup_2026_10/FINDINGS.md').read_text()
section=text.split('The required ignored-file check covers',1)[1]
paths=section.split('```text\n',1)[1].split('```',1)[0].strip().splitlines()
assert len(paths)==14 and all(Path(p).is_file() for p in paths)
subprocess.run(['git','status','--ignored','--short','--untracked-files=all','--',*paths],check=True)
PY
```

```text
!! _scratch/s11c/s11c-frequency-source-20260919/finish-01/complete/left-frequency-pencil.pickle
!! _scratch/s11c/s11c-frequency-source-20260919/finish-01/complete/right-frequency-pencil.pickle
!! _scratch/s11c/s11c-parallel-near-unity-20261001/uniform-01/complete/checks.json
!! _scratch/s11c/s11c-parallel-near-unity-20261001/uniform-01/complete/operation-index.jsonl
!! _scratch/s11c/s11c-parallel-near-unity-20261001/uniform-01/complete/operations.sqlite
!! _scratch/s11c/s11c-parallel-near-unity-20261001/uniform-continuation-01/complete/checks.json
!! _scratch/s11c/s11c-parallel-near-unity-20261001/uniform-continuation-01/complete/operations.sqlite
!! _scratch/s11c/s11c-thickness-coordinate-20260914/end_pairing_left/complete.pickle
!! _scratch/s11c/s11c-thickness-coordinate-20260914/end_pairing_right/complete.pickle
!! _scratch/s11c/s11c-thickness-coordinate-20260914/end_source_left/objects.pickle
!! _scratch/s11c/s11c-thickness-coordinate-20260914/end_source_right/objects.pickle
!! _scratch/s11c/s11c-uniform-response-20260919/production/complete/uniform-response.pickle
!! _scratch/s11c/s11c-uniform-source-20260919/production/complete/common.pickle
!! _scratch/s11c/s11c-uniform-source-20260919/production/complete/uniform-source.pickle
```

`!!` means present, ignored and untracked. None of these original scratch files was deleted.
The final tracked-file census is reported at the STOP against the baseline command
`git ls-tree -r --name-only dada3b7d | wc -l` (**9,472**); the current-tree command is
`git ls-files | wc -l`. The earlier Phase 3 review counted 10,089; its newly committed review file brought
this fix pass's starting tree to 10,090.


### Older tracked-but-ignored files

`git ls-files -ci --exclude-standard` still lists exactly these 20; ignore files are unchanged:

```text
research/pde_ledger_v2/_scratch/DIRECTIVE_step1_retire_apin.md
research/pde_ledger_v2/_scratch/NEXT_SESSION.md
research/pde_ledger_v2/_scratch/PLAN_apin_repair_and_throat_restart.md
research/pde_ledger_v2/_scratch/REVIEW_PROMPT_apin_restart.md
research/pde_ledger_v2/_scratch/REVIEW_PROMPT_directive_step1.md
research/pde_ledger_v2/_scratch/prior_art/SURVEY_dimension_libs.md
research/pde_ledger_v2/_scratch/prior_art/SURVEY_mutation_testing.md
research/pde_ledger_v2/_scratch/prior_art/SURVEY_provenance.md
research/pde_ledger_v3/_scratch/S11_card_review_prompt.md
research/pde_ledger_v3/_scratch/S11_directive_review_prompt.md
research/pde_ledger_v3/_scratch/S11_sympy_directive.md
research/pde_ledger_v3/_scratch/S11_sympy_repair_directive.md
research/pde_ledger_v3/_scratch/S11_sympy_repair_review_prompt.md
research/pde_ledger_v3/_scratch/S11_sympy_review_prompt.md
research/pde_ledger_v3/_scratch/S11_tex_directive.md
research/pde_ledger_v3/_scratch/S11_wl_directive.md
research/pde_ledger_v3/_scratch/S11_wl_repair_directive.md
research/pde_ledger_v3/_scratch/S11_wl_repair_review_prompt.md
research/pde_ledger_v3/_scratch/S11_wl_script_review_prompt.md
research/pde_ledger_v3/_scratch/skills_directive.md
```

### Conflicts and questions at this STOP

The three record conflicts in Part F remain open. H4's two export/producer mismatches are explicitly carried;
the current export-chain check is not all-pass. Uniform replay needs cache regeneration and that chain repair.
No Phase 4 history rebuild, Phase 5 deletion, new science or new review has started.
