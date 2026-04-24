# PDE Audit Simulation Bundle

This directory contains the executable simulation layer for the PDE V2 audit
release bundle.

The workflow is deliberately split into two phases:

1. `verify_reduced_fem.py` runs structural and manufactured-solution checks on
   the reduced FEM primitives.
2. `verify_nonlinear_solver.py` runs the pre-target nonlinear readiness checks:
   protocol/freeze hashing, source import boundary, Jacobian consistency,
   manufactured nonlinear mesh convergence, and continuation sanity.  It does
   not emit candidate packets.
3. `verify_physical_model.py` checks the ledger boundary for a strict physical
   nonlinear moving-throat exporter.  It passes only when the current bundle
   correctly reports that physical export is not yet permitted and that no
   physical packets were emitted.  The check also verifies source hashes and
   evidence phrases for every cited blocker and equation inventory item.
4. `generate_nonlinear_packets.py` exports a small target-blind nonlinear
   manufactured packet set under `output/nonlinear_packets/`, then the shared
   target-blind guard, evaluator, and obstruction diagnostic classify those
   frozen packets post-hoc.
5. `generate_reduced_sweep.py` builds frozen V2-22B solver packets from a
   predeclared reduced open-throat FEM sweep.  It uses only the local
   `reduced_fem.py` primitives and does not import the V2-21, V2-22A, V2-22B,
   or V2-23 target-evaluation modules.
6. `verify_target_blind.py` checks that the generator import boundary is still
   target-blind and that frozen packets do not contain target-output fields.
7. `evaluate_frozen_sweep.py` reads those frozen packets and runs the existing
   V2-22B -> V2-22A -> V2-21 audit chain.  Any target ranking in this phase is
   post-hoc classification only.
8. `diagnose_obstruction.py` performs post-hoc one-pole obstruction accounting
   on already-evaluated candidate reports.
9. `diagnose_required_deformation.py` maps the frozen failures to required
   coefficient changes, including one-pole `C`/`D0` multipliers, `A`
   reductions, and normalization multipliers.  This is post-hoc only and does
   not generate or retune candidates.
10. `diagnose_mechanism_gap.py` combines the deformation reports with the
    physical-model guard to state the current mechanism gap and the requirements
    for a future physical nonlinear exporter.

Run:

```bash
bash research/pde_audit/simulation/run_all.sh
```

Generated packets and summaries are written under `simulation/output/`.

The default protocol is `operator_v1`, a fixed grid of geometry, coupling, and
reduced-operator variants.  `smoke` remains available for a smaller quick check:

```bash
python research/pde_audit/simulation/generate_reduced_sweep.py --protocol smoke
```
