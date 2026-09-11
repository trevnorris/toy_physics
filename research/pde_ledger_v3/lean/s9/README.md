# S9 formalization pilot

The original D=3 pilot proves the local action variation, integrated
stationarity/PDE equivalence, real plane-wave reduction, and mode census. See
[RESULT.md](RESULT.md) for its mathematical assumptions and scope, and
[VERIFICATION.txt](VERIFICATION.txt) for the audit and source hashes.

Sources are in `S9Pilot/`, with the audit entry point `S9Pilot.lean`. Smooth
backgrounds live on R^4, and test variations are smooth and compactly supported.
The curl-only action and its material coefficients remain supplied premises.

The Lean environment is shared at [../](../README.md). Run from that parent:

```sh
lake build S9Pilot
lake env lean s9/S9Pilot.lean
```

The S10 generalization is in [../s10/](../s10/README.md), including exact D=3
specialization theorems connecting its definitions to these S9 definitions.
See [COVERAGE.md](COVERAGE.md) for the completed pilot scope and the remaining
obligations in the full S9 ledger step.
