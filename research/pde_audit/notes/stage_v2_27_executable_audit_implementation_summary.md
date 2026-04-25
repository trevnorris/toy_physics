# Stage V2-27 - Executable Audit Implementation Summary

## Purpose

This note summarizes what the code in `research/pde_audit/` now accomplishes.
It is intended as a paper-writing handoff for the executable/reproducibility
side of the audit.

The math/status notes explain what the audit means.  This note explains what
was implemented, what artifacts are generated, and what the code does and does
not certify.

---

## 1. Top-level referee harness

The main entry point is:

```bash
bash research/pde_audit/run_all.sh
```

It runs the release bundle in this order:

1. fixture manifest verification;
2. Python audit scripts;
3. simulation bundle;
4. root JSON hygiene check;
5. serial Mathematica mirrors.

The combined summary is written to:

```text
research/pde_audit/output/_summary.txt
research/pde_audit/output/_summary.json
```

The last full run reported:

```text
PYTHON: 28/28 passed, 0 failed
SIMULATION: 16/16 passed, 0 failed
MATHEMATICA: 10/10 passed, 0 failed
ROOT_JSON_FILES: 0
REFEREE_PASS: True
```

This harness is a reproducibility check.  It does not prove a completed PDE
branch exists.  It proves that the stated audit scripts, fixtures, mirrors, and
simulation diagnostics execute consistently and produce the recorded verdicts.

---

## 2. Python audit layer

The Python audit layer lives in:

```text
research/pde_audit/scripts/
```

It is run by:

```bash
bash research/pde_audit/scripts/run_all_audits.sh
```

The current suite has 28 passing checks.  It covers:

- parent wall action versus effective linear wall closure;
- Maxwell gauge localization;
- dimension ledger and port/source normalization warnings;
- open finite-radius junction impedance;
- Poisson/Newtonian hook;
- stable BdG support Schur complement;
- Maxwell/mixed kernel;
- Hamiltonian positivity gates;
- grouped real `P2` projector calculus;
- STF angular source-map normalization;
- grouped normalization and constant-prefactor branch algebra;
- outgoing `l=2` fingerprint;
- 2.5PN/4PN interface;
- branch freeze and no-refit protocol;
- weak-axisymmetric splitting;
- monomial quotient and similarity-orbit separation;
- isotropic full-bundle target surface;
- weak-form branch-extraction scaffold;
- branch-extraction fixture;
- profile-to-coefficient adapter;
- solver handoff validator;
- end-to-end smoke pipeline;
- V2-23 formula, mesh convergence, and minimal reduced branch solver;
- V2-24 negative controls.

The layer writes text summaries to:

```text
research/pde_audit/scripts/output/*.txt
research/pde_audit/scripts/output/_summary.txt
research/pde_audit/scripts/output/_summary.json
```

and generated JSON artifacts to:

```text
research/pde_audit/scripts/output/artifacts/
```

Important generated artifacts include:

- V2-21 branch extraction packet;
- V2-22A generated V2-21 manifest and observable packet;
- V2-22B validation reports;
- V2-22C generated profile and V2-21 manifests;
- V2-23 reduced solver branch manifest, observable packet, tolerance report,
  and mesh convergence report;
- V2-24 negative-control report.

---

## 3. Fixture and negative-control layer

Fixtures live in:

```text
research/pde_audit/scripts/fixtures/
```

The fixture manifest is:

```text
research/pde_audit/scripts/fixtures/MANIFEST.json
```

The fixture verifier checks that fixture files are present and hash-stable.

The negative-control suite under

```text
research/pde_audit/scripts/fixtures/negative_controls/
```

tests malformed or invalid solver packets, including:

- bad boundary protocol;
- missing gauge convention;
- nonfinite solver residual;
- nonmonotone grid;
- nonpositive mixed denominator;
- missing pre-target freeze;
- profile length mismatch;
- target-blind false;
- target-output leakage.

The purpose is to show that the handoff pipeline rejects invalid branches
before target residuals are trusted.

---

## 4. Mathematica mirror layer

The Mathematica mirrors live in:

```text
research/pde_audit/mathematica/
```

They are run serially by:

```bash
bash research/pde_audit/mathematica/run_all_audits.sh
```

or through the top-level harness.  Mathematica is run one script at a time with
`math -script` because of licensing constraints.

The current mirror suite has 10 passing checks:

- V2-04 open junction impedance;
- V2-13 grouped normalization ratio;
- V2-16 branch freeze/no-refit;
- V2-17 weak-axisymmetric splitting;
- V2-19 isotropic full-bundle target surface;
- V2-21 branch extraction fixture;
- V2-22A profile-to-coefficient adapter;
- V2-22B solver handoff validator;
- V2-22C end-to-end smoke pipeline;
- V2-23 formula audit.

Outputs are written to:

```text
research/pde_audit/mathematica/output/
```

The mirrors are secondary execution coverage.  They are not claimed as
independent derivations of the entire program; they check load-bearing algebra
and fixture contracts through a separate CAS/runtime.

---

## 5. Simulation bundle

The simulation layer lives in:

```text
research/pde_audit/simulation/
```

It is run by:

```bash
bash research/pde_audit/simulation/run_all.sh
```

The simulation layer is target-blind during candidate generation.  Target
residuals are applied only after packets are frozen.

The current suite has 16 passing checks:

1. `verify_reduced_fem.py`
2. `verify_nonlinear_solver.py`
3. `verify_physical_model.py`
4. `diagnose_notes_intake.py`
5. `generate_nonlinear_packets.py`
6. nonlinear target-blind guard
7. nonlinear frozen sweep evaluator
8. nonlinear obstruction diagnostic
9. nonlinear required-deformation diagnostic
10. `generate_reduced_sweep.py`
11. reduced target-blind guard
12. reduced frozen sweep evaluator
13. reduced obstruction diagnostic
14. reduced required-deformation diagnostic
15. mechanism-gap diagnostic
16. projection-stress diagnostic

The simulation writes:

```text
research/pde_audit/simulation/output/_summary.txt
research/pde_audit/simulation/output/_summary.json
```

plus manifests, packets, candidate reports, and diagnostic reports.

---

## 6. Reduced FEM and frozen packet generation

The reduced FEM primitives are in:

```text
research/pde_audit/simulation/reduced_fem.py
```

`verify_reduced_fem.py` checks:

- matrix structure;
- manufactured D/N half-wave behavior;
- open-shape profile smoke behavior.

`generate_reduced_sweep.py` builds 192 frozen V2-22B-compatible reduced
open-throat packets from a predeclared `operator_v1` protocol.

Generated reduced packet artifacts are written under:

```text
research/pde_audit/simulation/output/packets/
research/pde_audit/simulation/output/manifest.json
```

The reduced packets are then classified post-hoc by `evaluate_frozen_sweep.py`.
The current result is:

```text
0/192 target-passing frozen candidates.
189/192 open and stable candidates.
```

This is an honest target-blind miss, not a pipeline error.

---

## 7. Nonlinear manufactured readiness lane

The nonlinear readiness protocol is implemented in:

```text
research/pde_audit/simulation/nonlinear_protocol.py
research/pde_audit/simulation/verify_nonlinear_solver.py
research/pde_audit/simulation/generate_nonlinear_packets.py
```

`verify_nonlinear_solver.py` checks:

- source import boundary;
- protocol and freeze hashes;
- Jacobian directional consistency;
- manufactured nonlinear mesh convergence;
- continuation sanity.

`generate_nonlinear_packets.py` emits a small target-blind manufactured
nonlinear packet set under:

```text
research/pde_audit/simulation/output/nonlinear_packets/
research/pde_audit/simulation/output/nonlinear_manifest.json
```

The current result is:

```text
0/3 target-passing nonlinear manufactured candidates.
3/3 open and stable candidates.
```

This lane verifies nonlinear mechanics and frozen export plumbing.  It is not
the final physical moving-throat exporter.

---

## 8. Physical-export boundary guard

The strict physical nonlinear exporter boundary is implemented in:

```text
research/pde_audit/simulation/physical_nonlinear_model.py
research/pde_audit/simulation/verify_physical_model.py
```

The current physical-model guard passes by confirming that physical export is
still blocked.

It verifies:

- no banned target-evaluation imports;
- status hash stability;
- strict parent dynamic wall status is not passed;
- effective wall closure is available;
- physical export is not permitted;
- no physical packets or physical manifest were emitted;
- cited ledger blockers and equation-inventory items have source hashes and
  evidence phrases.

The current summary is:

```text
physical_export_permitted: False
packets_emitted: False
blocker_count: 4
```

This guard prevents the manufactured nonlinear lane from being accidentally
described as a true physical branch realization.

---

## 9. Target-blind guards and frozen evaluation

`verify_target_blind.py` checks that packet generators do not import target
evaluation modules and that frozen packet manifests do not contain target
outputs.

The same guard is used for:

- reduced packet generation;
- nonlinear manufactured packet generation.

`evaluate_frozen_sweep.py` runs post-hoc classification through the existing
V2-22B -> V2-22A -> V2-21 chain after packets are frozen.

The target-blind separation is the key anti-overfitting feature of the
simulation code.

---

## 10. Post-hoc diagnostics

The simulation bundle includes post-hoc diagnostics that read already-frozen
candidate reports.  They do not generate, mutate, or refit candidates.

### Obstruction diagnostics

Implemented in:

```text
research/pde_audit/simulation/diagnose_obstruction.py
```

The diagnostic decomposes misses through the one-pole ratio:

```text
D0*C/(3*A^2).
```

Current reduced open-stable range:

```text
min 0.0033775383274364888
median 0.06134471930726503
max 0.1353664855760648
```

### Required-deformation diagnostics

Implemented in:

```text
research/pde_audit/simulation/diagnose_required_deformation.py
```

It reports the coefficient changes that would be required after the fact.
Current reduced open-stable values:

```text
C_multiplier_min 7.387352901601946
C_multiplier_median 16.30132163440465
P0_multiplier_median 171.65261223353198
local_continuation False
```

These are diagnostic numbers, not retuning instructions.

### Mechanism-gap diagnostic

Implemented in:

```text
research/pde_audit/simulation/diagnose_mechanism_gap.py
```

It combines deformation reports with the physical-model guard.  It classifies
the current miss as a large one-pole support deficit and records the physical
requirements for a future exporter.

### Projection-stress diagnostic

Implemented in:

```text
research/pde_audit/simulation/diagnose_projection_stress.py
```

It asks what would happen under post-hoc algebraic coefficient projections.
It records:

```text
one_pole_support_alone_is_insufficient: True
uniform_outgoing_amplitude_scale_is_insufficient: True
target_blind_hit_claimed: False
```

This is the code-level evidence that a future mechanism needs outgoing
moment-shape control, not only scalar support or scalar outgoing amplitude.

---

## 11. Notes-intake guard

The unincorporated-notes intake is implemented in:

```text
research/pde_audit/simulation/diagnose_notes_intake.py
```

It verifies 16 source anchors across:

- 5PN computational handoff;
- 5PN final packet shape;
- Family-1 support/source branch location;
- barrier audit summary;
- strict parent/effective wall status;
- no-refit status firewall;
- `U/V` dressing;
- microscopic export kernel;
- constant-prefactor `N2/N4` target equations;
- support-blind outgoing transfer split;
- leakage/work support-side lane;
- atomic finite-throat `P0/P2/P22` source;
- lepton scalar radiative `P0` flux hook;
- scalar hammer to `P22` veto;
- atom-work intrinsic `P22` bracing;
- simulation coefficient map.

Current output:

```text
pass: True
anchor_count: 16
failed_anchor_count: 0
actual_branch_packet_required: True
retune_current_candidates_allowed: False
outgoing_moment_shape_control_required: True
primary_next_artifact: actual_branch_protocol_v1
```

Generated artifacts:

```text
research/pde_audit/simulation/output/notes_intake_report.json
research/pde_audit/simulation/output/notes_intake_report.txt
```

This guard converts the unincorporated notes into evidence-backed source and
protocol constraints without claiming that they solve the branch.

---

## 12. Referee summary generator

The combined summary generator is:

```text
research/pde_audit/scripts/write_referee_summary.py
```

It reads:

```text
research/pde_audit/scripts/output/_summary.json
research/pde_audit/simulation/output/_summary.json
research/pde_audit/mathematica/output/_summary.json
```

and writes:

```text
research/pde_audit/output/_summary.json
research/pde_audit/output/_summary.txt
```

It also enforces the release hygiene condition that no stray root-level JSON
files remain under `research/pde_audit/`.

Recent fields added to the combined summary include:

- physical model export boundary status;
- notes-intake status;
- nonlinear export status;
- target-blind guard result;
- obstruction and required-deformation diagnostics;
- mechanism-gap classification;
- projection-stress result.

---

## 13. What the code establishes

The code establishes:

- the symbolic/reduced audit scripts execute successfully;
- fixtures and negative controls are stable and enforced;
- target-blind packet generation is separated from target evaluation;
- the current reduced and manufactured nonlinear branch families miss the
  target;
- the miss is quantitatively diagnosed;
- post-hoc coefficient projection is not claimed as evidence;
- strict physical nonlinear export remains blocked;
- unincorporated notes imply actual-branch extraction, not retuning.

The code does not establish:

- that a completed physical PDE branch exists;
- that known physics has been derived from one PDE;
- that the manufactured nonlinear lane is a physical branch;
- that the current miss can be repaired by tuning the existing candidates;
- that the unincorporated notes already provide calibrated `D0/C/P0/N2/N4`
  packets.

---

## 14. Paper-use summary

For paper writing, the implementation can be summarized as:

```text
The audit bundle contains a reproducible Python, Mathematica, fixture, and
simulation harness.  It verifies the stated reduced/effective algebra and
handoff contracts, rejects malformed branch packets, runs target-blind reduced
and manufactured nonlinear branch searches, and records the resulting miss
through obstruction, deformation, mechanism-gap, projection-stress, and
notes-intake diagnostics.  The final referee harness passes, but the current
simulation layer finds no target-passing frozen branch and explicitly blocks
reclassification of the manufactured nonlinear model as a completed physical
moving-throat exporter.
```

The short version:

```text
The code makes the audit reproducible and falsifiable.  It does not complete
the PDE.  It shows exactly where the current executable branch families fail
and what an actual physical exporter must supply next.
```
