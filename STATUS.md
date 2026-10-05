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
| Direct c2 `[0,2]` entry is set to zero without an order-counting justification; whether the full response needs a nonzero term is unresolved. | S11c-c2/d if reused; [joint disposition, mixed-term issue](research/pde_ledger_v3/_measurements/S11c_upstream_repair_review_joint_disposition.md). |
| Full cross-engine composition and accepted scattering/work/loss | S11c-b/c1/c2/d if reopened; [canonical d limits](research/pde_ledger_v3/steps/S11c_d_profile_conditioned_scattering.md). |
| Conversion observable and strong-edge interpretation | S11c-e; [split and d weak-limit obligation](research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md#actual-downstream-uses-and-next-step). |
| Dynamical drain/return source, boundary data and any change to the drain-frozen symmetry result | **S12**; [plan, S12](research/pde_ledger_v3/V3_STEP_PLAN.md). No automatic transfer of the one-engine no-leak rule. |
| Real defect profile, holder, response and charge-force sign | Q2/Q3/S22; [plan, charge phase/S22](research/pde_ledger_v3/V3_STEP_PLAN.md). |
| Physical light/sound speed ratio | S20a calibration; S22/R10 proposed derivation; [plan, S20a/S22](research/pde_ledger_v3/V3_STEP_PLAN.md). |
| Physical admissibility of the rotational-stiffness material | Separate material audit; [native interpretation, §§5.7, 14.3](docs/native_light_em_and_vortex_throat_interpretation.md). |
| Muonium/gravity corrections, without withdrawn throat geometry | S16 and later species matching; [correction note, §§2–5](notes/muonium_gravity_corrections.md). |

## Next and paused work

**S12 is next in the physics plan:** define bulk-to-brane order conversion and its separate boundary data.
It does not require a finished S11c leakage factor. [Plan, S12](research/pde_ledger_v3/V3_STEP_PLAN.md).
This cleanup does not start S12 or any scientific run.

The [throat/EM notes](docs/light_em_investigation_handoff.md) remain **exploratory and paused**, as input to
S22/Q2 (with S12/Q3 connections), not adopted constitutive laws. Mixed numerical recovery, centre-drive
implementation and new leakage workers remain parked. [Closeout, paused notes](research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md).

History is preserved in `archive/pre-cleanup-2026-10-04`; the cleanup record is under [cleanup_2026_10/](research/pde_ledger_v3/cleanup_2026_10/FINDINGS.md).
