---
title: Charge and electromagnetism script catalog
type: script_catalog
status: current
sources:
  - software/em_charge_attribute/run_emergent_em_construction.sh
  - software/em_charge_attribute/puncture_deflection_electric_sign_check.py
  - software/em_charge_attribute/puncture_deflection_electric_sign_independent_recompute.py
  - software/em_charge_attribute/magnetism_moving_throat_check.py
  - software/em_charge_attribute/run_native_p_gate.sh
  - software/em_charge_attribute/run_u1_body_dynamics.sh
  - software/em_charge_attribute/run_u1_body_b2.sh
  - software/em_charge_attribute/run_u1_body_phaseC.sh
last_updated: 2026-09-03
---

# Charge and electromagnetism scripts

These are specialized scientific checks. Their success means the encoded
algebra and controls behaved as expected; it does not by itself establish a
physical electromagnetic sector.

## Emergent compact-U(1) construction

### `software/em_charge_attribute/run_emergent_em_construction.sh`

- **Role:** main runner for the spin-ice/flux-endpoint construction.
- **Runs:** `emergent_em_sympy.py`, `emergent_em_dual.wl`, then
  `compare_engines.py`, with a ten-minute timeout per process.
- **Checks:** incidence and charge algebra, endpoint dressing, hopping
  continuity, Maxwell/scalar static kernels, transverse mode count,
  dimensions, and declared negative controls.
- **Outputs:** JSON and logs under
  `software/em_charge_attribute/reports/artifacts/`, including
  `dimensional_firewall.log` and `engine_agreement.log`.
- **Interpretation:** an internally consistent possible realization. Phase
  existence is cited rather than computed, and the throat embedding is not
  dynamically derived.

Source: `software/em_charge_attribute/run_emergent_em_construction.sh` and
`software/em_charge_attribute/reports/emergent_em_construction.md` heading
`Dual-engine log`.

## Electric-sign calculation

### `software/em_charge_attribute/puncture_deflection_electric_sign_check.py`

- **Role:** primary SymPy calculation and able-to-fail campaign.
- **Inputs:** source assumptions embedded/read by the checker; no physical
  target value.
- **Checks:** field reduction, boundary classifier, coupled kernel, four
  ensemble functionals, dimensions, exhaustive decision table, and atomic
  mutations.
- **Outputs:** a structured human-readable calculation to stdout; `--json-only`
  emits the comparison payload.

### `software/em_charge_attribute/puncture_deflection_electric_sign_check.wl`

- **Role:** Wolfram implementation and cross-engine comparison.
- **Checks:** independently constructs the same kernel, ensembles, decision
  table, and mutations, and requests the Python JSON payload for comparison.
- **Outputs:** calculation, landing, mutation results, and `ENGINE_AGREE` to
  stdout.

### `software/em_charge_attribute/puncture_deflection_electric_sign_independent_recompute.py`

- **Role:** small independent SymPy check of the load-bearing result.
- **Checks:** (m_{gg}), determinant, V/M/J/MIXED pair coefficients and signs,
  mutation relations, and $1/R^2$ force falloff.
- **Outputs:** a short confirmation transcript; its stored summary is
  `puncture_deflection_electric_sign_independent_verdict.md`.
- **Interpretation:** the check confirms that the boundary ensembles disagree;
  it does not select one of them.

## Moving-throat magnetism

### `software/em_charge_attribute/magnetism_moving_throat_check.py`

- **Role:** primary symbolic calculation and mutation campaign.
- **Checks:** translated-source continuity, current identity, direct
  transverse kernel, Lorentz-anchor route, force derivatives, route
  independence, coefficient/cone ratios, dimensions, and decision table.
- **Outputs:** a human-readable calculation to stdout; `--json-only` emits the
  comparison payload.

### `software/em_charge_attribute/magnetism_moving_throat_check.wl`

- **Role:** Wolfram implementation and cross-engine comparison.
- **Checks/outputs:** reconstructs the two routes and controls, compares the
  Python payload, and prints the landing and `ENGINE_AGREE`.
- **Interpretation:** agreement establishes the conditional tensor and
  falloff, while the electric anchor and throat-current magnitude remain open.

## Native polar-sector constraint gate

### `software/em_charge_attribute/run_native_p_gate.sh`

- **Role:** main runner for the quadratic Hamiltonian/Dirac analysis.
- **Runs:** `native_p_gate_sympy.py`, six single-tooth Python ablations,
  `native_p_gate_dual.wl`, and `native_p_gate_compare.py`.
- **Checks:** whether native theories A/C contain a first-class Gauss
  descendant, plus six control theories that establish search sensitivity.
- **Outputs:** JSON/logs under `reports/native_p_gate_artifacts/` and the
  regenerated `reports/native_p_constraint_gate.md`.
- **Interpretation:** supports a generic quadratic no-Gauss result; the tuned
  rank-drop scan is not an exhaustive classification of every measure-zero
  sublocus.

## U1 throat-body dynamics

### `software/em_charge_attribute/run_u1_body_dynamics.sh`

- **Role:** staged Phase-A and Phase-B1 orchestrator.
- **Runs:** independent Phase-A and indexed-mechanics engines, their
  comparators, and batched mutation/liveness checks.
- **Outputs:** staged artifacts plus updates to
  `reports/u1_body_dynamics.md` and
  `reports/u1_body_dynamics_results.yaml`.
- **Operational note:** the first invocation deliberately halts after its
  stage-1 checkpoint; `--resume-stage2` performs the mutation/aggregation
  stage.

### `software/em_charge_attribute/run_u1_body_b2.sh`

- **Role:** heavy Phase-B2 intake/radiative workflow.
- **Behavior:** performs sandboxed, network-disabled, trace-authenticated
  stage/resume runs across the B2 engines, schemas, comparators, mutations, and
  merger.
- **Outputs:** B2 evidence directories and the aggregate report/results.
- **Operational note:** requires explicit startup/resume contract values and
  is not an ordinary quick validation command.

### `software/em_charge_attribute/run_u1_body_phaseC.sh`

- **Role:** thin dispatcher for the Phase-C stage runner and its shell probes.
- **Interpretation:** the presence of this runner is not evidence that Phase C
  completed; the aggregate report currently marks Phase-C gates not run.

Related source capsules:

- `memory/sources/software/em-emergent-construction.md`
- `memory/sources/software/em-electric-sign.md`
- `memory/sources/software/em-magnetism.md`
- `memory/sources/software/em-native-p.md`
- `memory/sources/software/em-u1-body-dynamics.md`
