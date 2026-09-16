# Homogeneous S11 Lean work

This is the first bounded S11 contract under
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
