# Nonlinear Protocol V2

This is the pre-target protocol for the next simulation layer.  It is not a
target-evaluation script and it should not import V2-21, V2-22A, or V2-22B
target extraction modules.

## Goal

Replace the current reduced and manufactured surrogate generators with a real
nonlinear PDE/continuation exporter that writes frozen V2-22B-compatible
packets.  The candidate set, continuation controls, solver tolerances, and
packet freeze hash must be fixed before any target residuals are computed.

The current ledger does not yet permit that strict physical export.  The
available materials support an effective wall/interface closure and reduced
BdG/Maxwell blocks, but the parent-level `S_eta`/`S_Sigma` promotion and full
coupled physical residual equations are not frozen.  Until that happens,
`verify_physical_model.py` is the executable guard: it records the blocker and
emits no physical branch packets.

## Unknowns

The nonlinear solve should at minimum evolve these coupled fields on an open
finite throat:

- throat radius/support state `R(s)` with `R(0)=R_mouth` and finite open exit
  impedance at `s=L`;
- wall/support profile `chi_eta(s)`;
- BdG support profile(s) `phi_B(s)`;
- mixed port profiles `u(s)` and `w(s)`;
- continuation parameters for geometry, couplings, and exit impedance.

## Boundary Conditions

- mouth: Dirichlet anchoring for support profiles and fixed positive
  `R_mouth`;
- exit: open impedance/Robin condition with finite `R_exit > 0`;
- forbidden: hard caps, closed-throat endpoint substitutions, and any
  post-target retuning.

## Required Pre-Evaluation Checks

Before target classification, the exporter must pass:

- manufactured-solution checks for each discretized operator;
- finite-difference or automatic-differentiation Jacobian consistency checks;
- Newton residual and continuation residual thresholds;
- mesh-refinement repeats on at least three grid sizes;
- branch-freeze hash verification;
- negative controls that reject hard-cap, nonpositive-frequency, missing-profile,
  and post-target-mutated packets.

The current executable readiness layer is:

- `nonlinear_protocol.py`: machine-readable pre-target protocol and candidate
  freeze-hash helpers;
- `verify_nonlinear_solver.py`: source import boundary, protocol/freeze hashes,
  manufactured nonlinear mesh convergence, analytic-vs-finite-difference
  Jacobian consistency, and continuation sanity.
- `generate_nonlinear_packets.py`: target-blind V2-22B packet exporter for the
  manufactured nonlinear readiness branches.
- `physical_nonlinear_model.py` and `verify_physical_model.py`: executable
  inventory of the strict physical-export blockers and source citations.

The readiness verifier intentionally emits no V2-22B candidate packets.  Packet
emission is isolated in `generate_nonlinear_packets.py`, and those packets are
then checked by the shared target-blind guard before post-hoc target evaluation.

The physical-model guard also emits no packets.  It should be replaced by a
physical exporter only after the parent action and coupled residual equations
are promoted and frozen in the ledger.

## Freeze Boundary

The exporter must write packets under `simulation/output/packets/` with:

- `freeze.pre_target_freeze = true`;
- `freeze.target_blind = true`;
- `freeze.no_post_residual_refit = true`;
- `freeze.candidate_freeze_hash` over the predeclared protocol, candidate
  parameters, solver tolerances, mesh, and source revision;
- no target residuals, target pass flags, target values, or post-hoc scores.

After the packet set is frozen, the existing chain remains:

```text
verify_target_blind.py
evaluate_frozen_sweep.py
diagnose_obstruction.py
diagnose_required_deformation.py
diagnose_mechanism_gap.py
```

The current reduced FEM results should remain in the bundle as a negative
control and baseline for the nonlinear exporter.
