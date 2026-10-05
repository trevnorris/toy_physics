# S9 formalization pilot

Read [the Lean scope policy](../FORMALIZATION_POLICY.md) before extending this
pilot. Identify a specific theorem or fidelity/coverage-contract gap and reuse
the existing proofs; exhaustive CAS transcript translation is not the goal.

The original D=3 pilot proves the local action variation, integrated
stationarity/PDE equivalence, real plane-wave reduction, and mode census. See
[RESULT.md](RESULT.md) for its mathematical assumptions and scope, and
[VERIFICATION.txt](VERIFICATION.txt) for the audit and source hashes.

Sources are in `S9Pilot/`, with the audit entry point `S9Pilot.lean`. Smooth
backgrounds live on R^4, and test variations are smooth and compactly supported.
The curl-only action and its material coefficients remain supplied premises.

The bounded S9 completion work adds the scalar-phase velocity theorem in
[Madelung.lean](S9Pilot/Madelung.lean) and the compact original-engine source
connection in [FIDELITY.md](FIDELITY.md). The theorem rules out a nonzero
transverse velocity polarization for the stated smooth scalar-phase plane waves;
it does not prove a universal no-photon result. The bounded C1–C4 contract is
complete: verification and mutation controls pass, and both independent fidelity
reviews are CLEAR. See [CLOSURE_VERIFICATION.txt](CLOSURE_VERIFICATION.txt) for
the build record and [FIDELITY_REVIEW.md](FIDELITY_REVIEW.md) for closure.

The Lean environment is shared at [../](../README.md). Run from that parent:

```sh
lake build S9Pilot
lake env lean s9/S9Pilot.lean
```

When sharing resources with calculations, the closure instrument rebuilds just
the S9 modules sequentially, runs the axiom audit and isolated mutations, and
limits each Lean process to one thread and 4096 MiB of Lean-managed memory:

```sh
python3 ../_measurements/S9_lean_contract_check.py
```

Its durable report is `../_measurements/S9_lean_contract_checks.json`; scratch
logs are under `s9/_scratch/closure/`. This is a Lean allocator limit, not an OS
resident-memory cap. The compact native SymPy source check is separately
reproducible with `python3 ../_measurements/S9_lean_source_check.py`; it does not
run the full CAS audits or write their exports.

The S10 generalization is in [../s10/](../s10/README.md), including exact D=3
specialization theorems connecting its definitions to these S9 definitions.
See [COVERAGE.md](COVERAGE.md) for the completed pilot scope and the remaining
obligations in the full S9 ledger step.
