# S11 Lean contracts

The homogeneous contract below is the first bounded S11 contract under
[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md). Read
[COVERAGE.md](COVERAGE.md) for H1–H4 and [FIDELITY.md](FIDELITY.md) for the precise
claims, native action/operator correspondence and limits. Local verification
passed; see [VERIFICATION.txt](VERIFICATION.txt). Both independent fidelity
reviews returned CLEAR, and **H1–H4 are complete**. The
[closure record](FIDELITY_REVIEW.md) retains the reviewed revision and findings.

The supplied homogeneous action adds compression to the existing curl stiffness.
The proofs reuse S10's calculus and geometry, classify full amplitude kernels
including frequency coincidence, recover B=0 and handle k=0 separately. The
bulk threshold theorem is kinematic; interfaces, bound states, leakage and
nonuniform scattering remain separate work.

| Source | Responsibility |
|---|---|
| [Action.lean](S11Homogeneous/Action.lean) | Actual density, integrated variation, local PDE, phase-average normalization and modal operator. |
| [Spectrum.lean](S11Homogeneous/Spectrum.lean) | Exhaustive full kernels, D3 counts, positive frequencies, coefficient coincidence and boundary limits. |
| [Threshold.lean](S11Homogeneous/Threshold.lean) | Phase-matching identity and below/on/above-grazing sign classification. |
| [S11Homogeneous.lean](S11Homogeneous.lean) | Selected load-bearing axiom audits. |

From the parent `lean/` directory, the resource-conscious verification command is

```sh
python3 ../_measurements/S11_lean_contract_check.py
```

It builds only the three S11 modules and the audit root, sequentially, then runs
isolated controls. Logs live in `s11/_scratch/verification/`; the durable result
is `../_measurements/S11_lean_contract_checks.json`. `--reuse-build` requires
identical local transitive sources, dependency pins, commands and output-module
hashes. The 4096 MiB setting limits Lean's allocator, not OS resident memory.

The compact source check is
`python3 ../_measurements/S11_lean_source_check.py`. It executes the selected
native D3 MAIN constructor and two small modal routes, not the production audit.
No new installation or dependency update is needed.


## D2 invariant-space contract I1–I4

The separately authorized [invariant contract](INVARIANT_NEXT.md) proves the
complete quadratic invariant spaces under SO(2) and O(2) conjugation, with
dimensions 4 and 3 and a one-dimensional reflection-odd complement. Local
verification passes and both independent fidelity reviews returned CLEAR.
**I1–I4 are complete**; see the
[review and closure record](INVARIANT_FIDELITY_REVIEW.md).
Read [INVARIANT_FIDELITY.md](INVARIANT_FIDELITY.md) for the object conventions,
minimal native Q9 transpose correction and its exact D2 span check, and
[INVARIANT_VERIFICATION.txt](INVARIANT_VERIFICATION.txt) for the evidence.

From `lean/`, its separate sequential verification command is

```sh
python3 ../_measurements/S11_lean_invariant_contract_check.py
```

The compact source command is
`python3 ../_measurements/S11_lean_invariant_source_check.py`. These use separate
reports from H1–H4. The original orientation-probe script is historical and
expects the pre-repair native source. The new contract does not certify D3–D5,
EL/total-divergence classes, production exports, or S11c calculations.

## D2 odd-invariant dynamics E1–E4

The next bounded increment is [DYNAMICS_COVERAGE.md](DYNAMICS_COVERAGE.md), with
its precise object and normalization map in
[DYNAMICS_FIDELITY.md](DYNAMICS_FIDELITY.md). **E1–E4 is complete: local verification passed and both independent reviews
returned CLEAR.** See [DYNAMICS_VERIFICATION.txt](DYNAMICS_VERIFICATION.txt) and
[DYNAMICS_FIDELITY_REVIEW.md](DYNAMICS_FIDELITY_REVIEW.md). It addresses the actual bulk first variation of the D2
odd term and its exhaustive longitudinal/transverse mixing criterion.

`S11OddDynamics/Action.lean` identifies `-beta P/2`, its derivative-defined
momenta, local PDE and modal operator. `Variation.lean` derives the finite
relative-action variation using existing S10 calculus. `Mixing.lean` identifies
the zero loci and supplies a smooth background with nonzero compact-test first
variation whenever beta is nonzero. `S11OddDynamics.lean` audits the new theorems.

Separate commands, run sequentially from `lean/`, are

```sh
python3 ../_measurements/S11_lean_dynamics_source_check.py
python3 ../_measurements/S11_lean_dynamics_contract_check.py
```

The latter binds the nine unchanged local dependency sources to their compiled
objects, then checks the new modules, axioms and mutations with one Lean worker.
It supports hash-checked `--reuse-build`. These checks do not regenerate pinned
S11c exports or certify the complete XFORM_EXTRA spectrum.

## D3 invariant completeness J1–J4

The bounded contract is [D3_COVERAGE.md](D3_COVERAGE.md), with the exact object
and native span map in [D3_FIDELITY.md](D3_FIDELITY.md). **J1–J4 is complete:
local verification passed and both independent fidelity reviews returned
CLEAR.** See [D3_VERIFICATION.txt](D3_VERIFICATION.txt) and
[D3_FIDELITY_REVIEW.md](D3_FIDELITY_REVIEW.md).

`S11D3Invariants` classifies quadratic forms on all real 3×3 gradients under
the full SO/O conjugation action. Its proved census is 3/3/0, with the unique
trace-square, trace-of-square and Frobenius-square parameterization. This is a
density classification before any EL or total-divergence quotient.

The compact native check and one-worker formal checks run sequentially:

```sh
python3 ../_measurements/S11_lean_d3_source_check.py
python3 ../_measurements/S11_lean_d3_contract_check.py
```

The formal check supports hash-checked `--reuse-build`. It also checks that the
finite coordinate certificate matches its generator. Neither command reruns a
production audit or modifies S11c exports.
