# PDE V2 Claim-Check Matrix

Snapshot date: 2026-04-24

This matrix is the referee-facing index for `research/pde_audit/`. It separates
what is checked exactly, what is checked by fixtures, and what remains a
branch-level assumption or non-claim.

Run everything with:

```bash
bash research/pde_audit/run_all.sh
```

The Python scripts are the primary executable audits. Mathematica mirrors are
secondary execution coverage for the load-bearing algebra and fixture contracts;
they are not claimed as independent derivations.

The `simulation/` bundle is the first target-blind branch-search layer. Its
default `operator_v1` protocol generates a fixed grid of frozen V2-22B packets
across geometry, coupling, and reduced-operator variants before target
evaluation, then classifies those packets through the existing V2-22B -> V2-22A
-> V2-21 audit chain.

## Summary

| Stage | Claim under audit | Primary Python check | Mathematica mirror | Expected result | Remaining assumption or non-claim |
|---|---|---|---|---|---|
| V2-01 | Parent wall action versus effective linear wall closure | `scripts/stage_v2_01_parent_wall_action_sympy_audit.py` | none | Split verdict: effective closure passes; strict parent dynamic wall not proved | Full parent-level moving throat action still needs promoted wall dynamics |
| V2-02 | Maxwell gauge localization and finite zero-mode gauge fixing | `scripts/stage_v2_02_maxwell_gauge_localization_sympy_audit.py` | none | Pass with wording patch | Either pre-reduction Lorenz constraint or localized gauge-fixing profile must be declared |
| V2-03 | Dimension ledger and normalization consistency | `scripts/stage_v2_03_dimension_ledger_sympy_audit.py` | none | Pass with port/source normalization warning | Raw mechanical transfer factor needs explicit gravitational port normalization |
| V2-04 | Open finite-radius exit replaces hard cap; D/N support field must be flow/displacement-like | `scripts/stage_v2_04_open_junction_impedance_sympy_audit.py` | `mathematica/stage_v2_04_open_junction_impedance_mathematica_audit.wl` | Pass, with negative controls for hard cap and generic scalar Neumann | Pressure/potential scalar would need an impedance-transforming mechanism |
| V2-06/07 | Projected continuity gives Poisson hook and inverse-square Green function | `scripts/stage_v2_06_07_poisson_newtonian_sympy_audit.py` | none | Conditional pass | Universal compact-source strength is not derived by this identity |
| V2-08 | Stable BdG support Schur complement and softening gate | `scripts/stage_v2_08_bdg_wall_schur_stability_sympy_audit.py` | none | Pass for positive-energy modes | Negative-Krein or ghost modes are outside this closure |
| V2-09 | Maxwell/mixed kernel, stability denominator, and nonnegative outgoing transfer | `scripts/stage_v2_09_maxwell_mixed_kernel_sympy_audit.py` | none | Pass in declared one-lane closure | Full nonlinear PDE branch realization is not proved |
| V2-10 | Hamiltonian stability gates for reduced conservative bundle | `scripts/stage_v2_10_hamiltonian_stability_sympy_audit.py` | none | Pass under positive-energy gates | Requires no ghost or Krein-sector modes |
| V2-11 | Grouped real `P2` projector calculus | `scripts/stage_v2_11_grouped_p2_projectors_sympy_audit.py` | none | Pass | Independent of radial/axial open-throat boundary data |
| V2-12 | STF angular source-map normalization | `scripts/stage_v2_12_stf_angular_source_map_sympy_audit.py` | none | Pass | Remaining normalization is radial/axial/port data |
| V2-13 | Grouped normalization ratio and constant-prefactor branch algebra | `scripts/stage_v2_13_grouped_normalization_ratio_sympy_audit.py` | `mathematica/stage_v2_13_grouped_normalization_ratio_mathematica_audit.wl` | Pass | Actual PDE branch must still hit the target surface |
| V2-14 | Outgoing `l=2` fingerprint and odd branch | `scripts/stage_v2_14_outgoing_l2_fingerprint_sympy_audit.py` | none | Pass with wording patch | Specific to the brane/worldtube 3-space STF outgoing port |
| V2-15 | 2.5PN to 4PN quadrupole/tail interface | `scripts/stage_v2_15_25pn_4pn_interface_sympy_audit.py` | none | Pass with tail transport gate | Requires `Theta_tail (c/c_s)^3 = 1` or a derived equivalent |
| V2-16 | Branch freeze and no-refit protocol | `scripts/stage_v2_16_branch_freeze_no_refit_sympy_audit.py` | `mathematica/stage_v2_16_branch_freeze_no_refit_mathematica_audit.wl` | Pass | Target comparison is not trustworthy without a pre-target freeze packet |
| V2-17 | Weak-axisymmetric grouped-P2 splitting and hidden-even condition | `scripts/stage_v2_17_weak_axisymmetric_splitting_sympy_audit.py` | `mathematica/stage_v2_17_weak_axisymmetric_splitting_mathematica_audit.wl` | Pass | Future actual branch must compute the remaining prefactor-loading scalar |
| V2-18 | Monomial quotient and similarity-orbit separation | `scripts/stage_v2_18_monomial_quotient_similarity_orbit_sympy_audit.py` | none | Pass | Does not prove an actual branch tangent stays in the zero-defect subspace |
| V2-19 | Isotropic full-bundle target surface and local codimension | `scripts/stage_v2_19_isotropic_full_bundle_target_surface_sympy_audit.py` | `mathematica/stage_v2_19_isotropic_full_bundle_target_surface_mathematica_audit.wl` | Pass | Frozen branch output must still satisfy the target packet |
| V2-20 | Weak-form and numerical branch-extraction preparation | `scripts/stage_v2_20_weak_form_branch_extraction_sympy_audit.py` | none | Pass | This is a scaffold, not a nonlinear PDE solve |
| V2-21 | Branch extraction fixture from frozen manifest | `scripts/stage_v2_21_branch_extraction_fixture.py` with `scripts/fixtures/stage_v2_21_sample_branch_manifest.json` | `mathematica/stage_v2_21_branch_extraction_fixture_mathematica_audit.wl` | Pass; calibrated control passes, primitive demo honestly misses target | Fixture does not rescue branches after residual evaluation |
| V2-22A | Profile-to-coefficient adapter | `scripts/stage_v2_22a_profile_to_coefficient_adapter.py` with `scripts/fixtures/stage_v2_22a_profile_input_manifest.json`; strict-profile rejection checked by V2-24 | `mathematica/stage_v2_22a_profile_to_coefficient_adapter_mathematica_audit.wl` | Pass | Adapter does not solve the PDE; it converts sampled or analytic profiles; trusted-local `expr` profiles require non-strict mode |
| V2-22B | Solver handoff schema and open-throat validation | `scripts/stage_v2_22b_solver_handoff_validator.py` with valid, hard-cap, and V2-24 negative fixtures | `mathematica/stage_v2_22b_solver_handoff_validator_mathematica_audit.wl` | Pass; hard-cap, grid, profile-length, freeze, target-blind, mixed-Delta, gauge, boundary-protocol, solver-metadata, and target-leakage controls are rejected | Additional malformed PDE-export cases can be added as the real solver schema grows |
| V2-22C | End-to-end smoke pipeline from solver packet to target residuals | `scripts/stage_v2_22c_end_to_end_smoke_pipeline.py` with V2-22C fixtures | `mathematica/stage_v2_22c_end_to_end_smoke_pipeline_mathematica_audit.wl` | Pass; stable smoke branch may target-fail, calibrated control can pass | A target-failing stable branch is a valid falsification, not a pipeline error |
| V2-23 | Minimal reduced open-throat solver and first real residual extraction | `scripts/stage_v2_23_minimal_open_throat_branch_solver.py`, `scripts/stage_v2_23_formula_sympy_audit.py`, and `scripts/stage_v2_23_mesh_convergence_audit.py` | `mathematica/stage_v2_23_formula_mathematica_audit.wl` | Pass; open/stable branch honestly misses target; reduced solver values converge under mesh refinement with stable freeze hash and decreasing pairwise deltas | Still a reduced linear-FEM branch prototype, not a full nonlinear PDE branch export |
| V2-24 | Negative-control suite for hardening failure behavior | `scripts/stage_v2_24_negative_controls.py` with `scripts/fixtures/negative_controls/*.json` | none | Pass when invalid packets and unstable branches are rejected | Extend with solver-specific controls once a full PDE exporter exists |

## Reproducibility Artifacts

| Artifact | Purpose |
|---|---|
| `scripts/fixtures/MANIFEST.json` | Canonical fixture filenames, schemas, roles, and SHA256 hashes |
| `scripts/output/_summary.txt` | Human-readable Python audit summary |
| `scripts/output/_summary.json` | Machine-readable Python audit summary with output and artifact hashes |
| `scripts/output/artifacts/*.json` | Generated downstream JSON artifacts from fixture-backed Python runs |
| `scripts/output/artifacts/stage_v2_23_mesh_convergence_report.json` | Mesh-refinement diagnostics for the reduced V2-23 solver |
| `scripts/output/artifacts/stage_v2_24_negative_controls_report.json` | Executable negative-control verdicts and expected rejection paths |
| `simulation/INTERPRETATION.md` | Human-readable interpretation of what the reduced simulation result does and does not establish |
| `simulation/NONLINEAR_PROTOCOL_V2.md` | Pre-target protocol for the next nonlinear PDE/continuation exporter |
| `simulation/output/reduced_fem_verification.json` | Manufactured-solution and matrix-structure checks for reduced FEM primitives |
| `simulation/output/nonlinear_protocol.json` | Machine-readable pre-target nonlinear readiness protocol and protocol hash |
| `simulation/output/nonlinear_readiness.json` | Nonlinear readiness report: import boundary, freeze hashes, Jacobian consistency, manufactured mesh convergence, and continuation sanity |
| `simulation/output/physical_model_status.json` | Evidence-backed executable guard that strict physical nonlinear export remains blocked until parent wall dynamics and coupled residual equations are frozen |
| `simulation/output/nonlinear_manifest.json` | Frozen target-blind nonlinear manufactured packet manifest |
| `simulation/output/nonlinear_target_blind_report.json` | Target-blind guard for the nonlinear packet exporter and frozen nonlinear packets |
| `simulation/output/nonlinear_evaluation_summary.json` | Post-hoc classification of frozen nonlinear packets through the V2-22B -> V2-22A -> V2-21 chain |
| `simulation/output/nonlinear_obstruction_report.json` | Post-hoc one-pole decomposition of frozen nonlinear packets |
| `simulation/output/nonlinear_required_deformation_report.json` | Post-hoc required-deformation map for frozen nonlinear manufactured packets |
| `simulation/output/nonlinear_packets/*.json` | V2-22B-compatible frozen packets emitted by the nonlinear manufactured exporter |
| `simulation/output/nonlinear_candidate_reports/*.json` | Per-candidate post-hoc evaluation reports for nonlinear packets |
| `simulation/output/manifest.json` | Frozen target-blind simulation packet manifest and protocol hash |
| `simulation/output/target_blind_report.json` | Executable guard that checks generator imports and confirms frozen packets do not contain target-output fields before evaluation |
| `simulation/output/evaluation_summary.json` | Post-hoc classification of frozen simulation candidates, score distribution, and dominant residual accounting |
| `simulation/output/obstruction_report.json` | Post-hoc one-pole decomposition of frozen candidates, including `D0*C/(3A^2)` and quartic-reservoir shortfall |
| `simulation/output/required_deformation_report.json` | Post-hoc map from frozen failures to required one-pole and normalization coefficient changes |
| `simulation/output/mechanism_gap_report.json` | Post-hoc summary connecting the simulated miss, required deformations, and blocked physical exporter requirements |
| `simulation/output/packets/*.json` | V2-22B-compatible frozen solver packets generated by the simulation phase |
| `simulation/output/candidate_reports/*.json` | Per-candidate V2-22B/V2-22A/V2-21 evaluation reports |
| `mathematica/output/_summary.txt` | Human-readable Mathematica mirror summary |
| `mathematica/output/_summary.json` | Machine-readable Mathematica mirror summary |
| `mathematica/output/_environment.txt` | Serial Mathematica kernel/environment capture |
| `output/_summary.txt` | Top-level referee harness summary |
| `output/_summary.json` | Combined referee harness summary |

## Current Referee Gaps

- Add adaptive-mesh and independent-quadrature repeats once the V2-23 branch
  solver becomes nonlinear.
- Promote and freeze `S_eta`/`S_Sigma` plus the coupled physical residual
  equations before replacing the physical-model guard with a real nonlinear
  PDE/continuation exporter.
- Add cross-CAS mirrors for the newer V2-23 mesh and V2-24 negative-control
  Python-only checks if they become central claims rather than harness guards.
- Add more malformed-export controls only as the real PDE exporter schema grows:
  continuation path provenance, multiple branch families, and solver-specific
  boundary-condition metadata.
