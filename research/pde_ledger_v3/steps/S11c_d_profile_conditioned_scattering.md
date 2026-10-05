# S11c-d — supplied-profile scattering: PARTIAL

## Where the light sector stands

S11c asked how much light escapes into other motion when the brane changes. There is **no accepted nonuniform light-loss number**. Reflected light counts as surviving light, along with transmitted light. **UNRESOLVED.** [Governing specification, §0](../directives/S11c_d_SHARED_PHYSICS.md); [benchmark, opening and observable](../_measurements/S11c_d_numerical_radiating_balance_report.md).

On a uniform brane, small light waves do not drive the bulk in the stated linear model; stable propagation requires nonnegative stiffness. Independent SymPy and Wolfram calculations support this, with engine/repair reviews and a documented comparison; one Wolfram review's full rerun was blocked. **CONDITIONAL.** [S11b, transverse mode, comparison and operational note](S11b_interface_coupling_law.md).

A later SymPy check found selected light modes finite and decoupled at the tested light/sound speed matches, with the bulk at rest and other inputs fixed. Claude and Grok cleared the method, not a fresh independent review of the executed result. **CONDITIONAL.** [Uniform result, coverage and limits](../_measurements/S11c_d_near_unity_uniform_continue_result.md).

Where the brane changes, conversion routes exist in the equations, but their effect on total light survival is unresolved. The completed benchmark's final signal is below its empirical numerical floor, not a physical bound. Equation repairs have only scoped review support; an unjustified zero in the direct mixed contribution is still open. **UNRESOLVED.** [Benchmark, controls](../_measurements/S11c_d_numerical_radiating_balance_report.md); [joint repair disposition, mixed-term issue](../_measurements/S11c_upstream_repair_review_joint_disposition.md).

**OPEN:** equation/response debts stay with S11c-b/c2/d if reopened, the observable and strong-edge interpretation with S11c-e, real throat/charge response with Q2/Q3/S22, the speed ratio with S22/R10, and rotational material admissibility with its separate audit. **S12 is next**, for drain/return dynamics; it needs no completed leakage number. [Plan, S12, charge preamble and S22](../V3_STEP_PLAN.md); [native interpretation, §§5.7, 14.3](../../../docs/native_light_em_and_vortex_throat_interpretation.md).

## Scope and evidence boundary

**ESTABLISHED — document disposition:** d closes PARTIAL by user direction. There is no automatic continuation, new review or defect sweep. The [umbrella closeout](S11c_PARTIAL_CLOSEOUT.md) links a/b/c1/c2 and assigns downstream debts. This record reports existing evidence, not a new calculation.

**ESTABLISHED — scoped specification review:** the v10 supplied-profile spec received Claude/Opus **SOUND** and Grok **SOUND**. It retains independent height/slope orders, a localized profile with constant ends and a full coupling vertex; c2 operand debt is explicit. The scattering-form amendment received Claude/Grok **CLEAR** for its scoped amendment/inventory; the total-transverse-loss amendment received **CLEAR FOR THIS PREPARATION AMENDMENT** from both. These do not clear a worker or result. [Spec, §§0–1b](../directives/S11c_d_SHARED_PHYSICS.md); [Claude verdict](../directives/_legs/S11c_d_shared_physics_review_r10_opus.md); [Grok verdict](../directives/_legs/S11c_d_shared_physics_review_r10_grok.md); [form disposition, opening](../_measurements/S11c_d_scattering_form_review_round09_adjudication.md); [total-loss disposition, opening](../_measurements/S11c_d_total_transverse_loss_review_disposition.md).

## What the uniform check establishes

**CONDITIONAL:** the fixed model is omega=3, strict rest bulk, LAB_HELD/RHO4_CONSTANT and tangential momenta (1/5,1/10), with a selected effective sound-speed family. Twelve exact speeds give 24 end evaluations and 48 normal-sign checks, each for two selected transverse polarizations. Modes and nonzero current stay finite; all 1,536 selected face-channel records are zero. The matching speeds are LEFT **sqrt(3/2)** and RIGHT **sqrt(150/101)**. The saved limits from the two sides supply the exact matches; singular raw substitutions are not ignored. [Uniform result, “Actual coverage and result”](../_measurements/S11c_d_near_unity_uniform_continue_result.md).

**CONDITIONAL — review scope:** both reviewers returned **CLEAR FOR THIS SELECTED UNIFORM METHOD**. The execution is SymPy, not a blind Wolfram replication; there was no fresh worker/result review. Finite samples and selected limits do not prove interval-wide smoothness, a complete mode census, flowing-background behavior or nonuniform confinement. [Uniform result, coverage and cost/limits](../_measurements/S11c_d_near_unity_uniform_continue_result.md).

## The omega=3 benchmark: numerical loss unresolved

**UNRESOLVED:** ten finite solves used one held profile, strict rest bulk, LAB_HELD/RHO4_CONSTANT, tangents (1/5,1/10), and selected transverse phase-speed ratio approximately **0.12**; the bare S9 coefficient ratio is separately **0.1**. No neighboring frequency or alternative anchoring was run. This is not the calibrated draining model. [Benchmark, opening](../_measurements/S11c_d_numerical_radiating_balance_report.md).

The measured quantity is one minus reflected-plus-transmitted transverse current divided by incident transverse current, including polarized cross terms. Positive is a deficit; negative is apparent gain. The original report gives:

| Full-contrast setting | Column 2 deficit (ppm) | Column 3 deficit (ppm) |
| --- | ---: | ---: |
| Baseline | -1.773911 | +1.633609 |
| Refined quadrature/basis/cutoff | -1.034869 | +0.942457 |
| Larger source/domain interval | -0.608245 | +0.554788 |
| Larger interval, regulator halved | -0.608218 | +0.554792 |

**UNRESOLVED:** the floor is empirical **1 ppm**, not a rigorous error bound. No column satisfies the report's positive-resolution criterion. Columns 1/4 stay within **0.000012 ppm** at these settings; that also is not a physical bound. Baseline apparent gain, substantial refinement/domain changes and unresolved contrast scaling prevent a loss interpretation. Small solve residuals and small regulator-halving changes do not certify continuum accuracy. [Benchmark, signed-current table and controls](../_measurements/S11c_d_numerical_radiating_balance_report.md).

**UNRESOLVED — review:** Claude and Grok both returned **NEEDS REVISION** for the original-tolerance route. Actual Gaussian matrix-route, momentum-zero and saved-operand row/trial checks remain unperformed; the coarse alternative was unused. No independent method/result clearance or corrected production scattering result is inherited. [Replacement disposition, verdict table](../_measurements/S11c_d_numerical_radiating_replacement_review_disposition.md); [benchmark, controls](../_measurements/S11c_d_numerical_radiating_balance_report.md).

## Equations: what was repaired and what remains open

**CONDITIONAL:** inertia sign, mechanical-load orientation, physical/reference pressure trace and thickness-coordinate repairs have distinct scoped reviews. The original b/c2 close applies to their pre-repair versions. Read [b, “Changed after close”](S11c_b_variable_coefficient_operator.md#changed-after-close) and [c2, “Changed after close”](S11c_c2_self_energy_fold.md#changed-after-close), rather than treating the old sign-offs as approval of current composition.

**UNRESOLVED:** `scripts/S11c_c2_selfenergy_fold_sympy_audit.py:400–418` sets the direct three-leg `[0,2]` entry to zero while retaining the iterated first-shape product. The retained grades include `eta*sigma_W`, so order counting does not justify that zero. **Whether the complete closed response has a nonzero direct term is unresolved.** Neither the focused source diagnostic nor later partial constructions quantify an effect on the benchmark. [Joint disposition, “The remaining mixed-term issue”](../_measurements/S11c_upstream_repair_review_joint_disposition.md); [focused tanh result, limits](../_measurements/S11c_upstream_mixed_tanh_result.md).

**UNRESOLVED — unreviewed repairs:** the separate Wolfram N6 pressure-trace repair corrected a reference-pressure double shift by extracting the reference-pressure map from the actual face law and solving its ordered inverse. Its checks are limited to the tested representation domain, not a full Wolfram self-energy engine; no review or comparator ran. The d-internal sheet repair resolves its fixed-positive-frequency counterexample but is also **unreviewed**. [Wolfram trace repair, opening and scope/limits](../_measurements/S11c_wolfram_pressure_trace_repair_report.md); [repair audit, line 16](../_measurements/S11c_wolfram_repair_audit_report.md); [sheet repair, scope](../_measurements/S11c_d_sheet_repair_report.md).

## Intermediate results are not a loss factor

**CONDITIONAL:** the SymPy symbolic/finite omega=1 work reached four full-rank 645-unknown cases and four incident directions with positive regulator and approximate boundaries. Its zero open-thickness selector rank holds at that input, not for all profiles/speeds. Full blind-engine scattering/export obligations remain. [Builder report, four-case completion and retained contract](../_measurements/S11c_d_sympy_builder_report.md); [continuum-current report, opening/limits](../_measurements/S11c_d_continuum_currents_report.md).

**WITHDRAWN/FAILED:** the omega=1 fixed-point continuation stopped at RIGHT with MemoryError; exit zero was not full success. Its two literal **CLEAR FOR THIS PROGRESS-DEPENDENT FIXED-POINT CONTINUATION** verdicts concerned the build, not results. [Fixed-point report, failure and review scope](../_measurements/S11c_d_transverse_face_fixed_point_continue_report.md).

**CONDITIONAL:** finite packet/source/matching evidence includes an inner-kernel bank, local Gaussian actions, polarization-dependent first-order forcing and matched transverse amplitudes. The handoff records positive end-current weights and first-order flux cancellation. The scoped receiving result excludes the recorded C3 real-axis poles; the zero-forcing U_B polarization has zero three-row response only in the stated weighted class, and the full five-field transverse poles remain. None supplies complete physical work, a full field or leakage. Reviews had different scopes: the local action retained Claude **NEEDS REVISION** / Grok **CLEAR** with a tooling repair; later Claude-only flux/receiving verdicts cleared source builds, not independent result execution. [Applicability status, accepted inner/local sections](../_measurements/S11c_d_defect_near_unity_applicability_status.md); [source result](../_measurements/S11c_d_first_order_source_result.txt); [matching result](../_measurements/S11c_d_first_order_transverse_matching_result.txt); [handoff, established results](../../../docs/light_em_investigation_handoff.md); [flux review, limits](../_measurements/S11c_d_first_order_transverse_flux_build_claude.md); [receiving review, limits](../_measurements/S11c_d_first_order_receiving_regular_build_r3_claude.md); [receiving method, §3 class/U_B restriction](../_measurements/S11c_d_first_order_receiving_regular_20261004_r2_plan.txt).

## Clean-condition result and its owner

**CONDITIONAL:** the clean-condition packet measured a symmetry selection rule on one engine, SymPy: `MIXED 0` and `WRONG_SIDE 0` in all four tested cases, with responsive controls. This was review-leg evidence, not a blind second-engine build. Claude literally said **“not clear.”**, Grok **“not cleared.”** for the broader packet. [Clean-condition disposition, round 5, verdicts and measurements](../directives/_measurements/S11c_d_clean_condition_review_disposition.md).

**OPEN — S12:** its no-leak conclusion needs the drain frozen and the relevant receiving channel empty. Nothing transfers automatically to live order conversion. The native drain/return functions and their boundary data belong to S12; trapped support belongs to Q2/S22. [Clean-condition disposition, R5-1/R5-2 and final question](../directives/_measurements/S11c_d_clean_condition_review_disposition.md); [plan, S12/Q2](../V3_STEP_PLAN.md).

## Handoff and prohibited claims

**OPEN:** S11c-e's observable/strong-edge calculation must recover d's weak limit if supplied; S22 inventory items 1/4 need nonuniform response and a supplied profile, with physical magnitude deferred to R1. Items 2/3 instead name S11b and S5–S7. The charge-phase preamble uses S11c as a question about linear support, not as proof of a holder. [Split, table/N5](../directives/S11c_decisions.md); [d spec, §3d](../directives/S11c_d_SHARED_PHYSICS.md); [plan, charge preamble/S22](../V3_STEP_PLAN.md).

**UNRESOLVED:** no claim of “below 1 ppm,” physical gain, matched-speed defect loss, harmless particles/strong edges, or unconditional confinement follows. Reflected light remains survival. Uniform decoupling, source zeros, first-order flux cancellation and conditional receiving regularity do not replace complete work and outgoing-power accounting. [Benchmark, controls](../_measurements/S11c_d_numerical_radiating_balance_report.md); [handoff, work state](../../../docs/light_em_investigation_handoff.md).

**OPEN — S22/Q2 inputs:** throat/EM notes remain exploratory and paused in [the investigation handoff](../../../docs/light_em_investigation_handoff.md). They do not select a support law, resolve material admissibility or authorize resuming this calculation. S12 can proceed as the next planned subject with these limits retained. [Phase 1 review, C4–C6/D4](../cleanup_2026_10/PHASE1_REVIEW.md).
