# S15 separation departure — W1/W2 repair run record

Contract inputs: Parts A (`91248b67`) and B (`7f52b4d2`); reviewed repair baseline `6c0fc66b`. Each engine reconstructs its objects independently. No cross-engine output comparison or review was performed. No git, git-annex, or datalad command was run.

## Wolfram

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s15/guard/reading_w1_w2_wolfram_final_02 -- math -script /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S15_separation_departure.wl
```

Guard exit `0`; child exit `0`; guard reason `None`; `limitsVerified: true`; stderr `0` bytes. Guard wall time **1609.4741610650672 s**; peak recorded cgroup memory **3547086848 bytes**, including descendants.

- Source: `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S15_separation_departure.wl`; SHA-256 `ff784cacfae7369e9c9ad5e32cf073b677b9a730bfad7f52bb0c97c2250f32fa`.
- Transcript: `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S15_separation_departure.out`; SHA-256 `b133fdcba87d06aa314773f51a688f1009e7ee8bc21a704dd59d8256d923aaed`.
- Inventory: 336 Q sections; 336 item-3 sections; 336 item-4 statuses; 48 separation-dependent objects with finite-identity conditions. Transcript end marker present: `True`.
- Finite-identity objects missing a condition: `[]`; missing restricted Q: `[]`. Remaining separation occurrences in the printed condition definitions: `[]`.

Requested-object gaps indexed by this final engine: none. General-field family integrals remain unevaluated where item 3 expressly permits them.

## SymPy

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s15/guard/reading_w1_w2_sympy_final_01 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S15_separation_departure.py
```

Guard exit `0`; child exit `0`; guard reason `None`; `limitsVerified: true`; stderr `0` bytes. Guard wall time **12545.394144776976 s**; peak recorded cgroup memory **1705218048 bytes**, including descendants.

- Source: `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S15_separation_departure.py`; SHA-256 `ccd4324e6a227d6daaf97d55adf3927cedc44ff125b9e50c3eb7c42bdeb2d5c7`.
- Transcript: `/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S15_separation_departure.out`; SHA-256 `138391ad246470561aadc290b6fab74885b38a446f2d06b39502bf6802c83857`.
- Inventory: 336 Q sections; 336 item-3 sections; 336 item-4 statuses; 48 separation-dependent objects with finite-identity conditions. Transcript end marker present: `True`.
- Finite-identity objects missing a condition: `[]`; missing restricted Q: `[]`. Remaining separation occurrences in the printed condition definitions: `[]`.

Requested-object gaps indexed by this final engine: none. General-field family integrals remain unevaluated where item 3 expressly permits them.

## Construction and scope

W1: each engine independently eliminates the direction variables and constructs polynomial differential systems from its closed forms. A finite zero-multiplicity bound yields parameter equations as finite sums and ordered products of fully defined sparse matrices. The source gives the bound’s argument, including rank drops. Piecewise exponent supports remain in the definitions; no physical approximation or numerical sample is used.

Each engine solves its computed geometric scalar equation and inserts its closed Q at a regular symbolic reference, with the same instantaneous coefficients, to obtain Q under each condition. The family responses remain live.

W2: Wolfram differentiates a retained analytic Appell function using computed coefficient identities; `Inactive[AppellF1]` suppresses native parameter simplification. Actual Piecewise supports and denominator zero loci are printed for closed item-3 objects and Q operands; the final audit found no uncovered loci. SymPy’s item-3 expressions are unchanged; status text prints their actual branch inventory.

Items 1, 2, 5, and 6 retain the reviewed general-field construction. The final transcript links were removed before the corresponding runs; both outputs are new plain files. Runs were sequential, with the specified 50 GiB cap, zero swap, and no deadline. Wolfram retains both required thread options. No annex target or repository metadata was modified.
