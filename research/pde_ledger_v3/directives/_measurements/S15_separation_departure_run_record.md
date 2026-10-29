# S15 separation departure — V1/V2 repair run record

Contract inputs: Parts A (`91248b67`) and B (`7f52b4d2`); reviewed baseline `5f0fc03a`. The engines reconstruct their objects independently. No cross-engine output comparison, review or commit was performed.

## Wolfram

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s15/guard/reading_repair_wolfram_final_01 -- math -script /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S15_separation_departure.wl
```

Guard exit `0`; child outcome `{"exitCode": 0, "guardReason": null, "wallSeconds": 230.1380869220011, "childPid": 1073413}`; `limitsVerified: true`; stderr `0` bytes. Runtime **230.20257373800268 s**; peak recorded cgroup memory **1631453184 bytes**.

- [Script](../../mathematica/S15_separation_departure.wl); SHA-256 `407a40d480c0bda73ec125f6f1f0b6ee2a42654148689d13f0c7fc9e55a2d4e6`.
- [Transcript](../../mathematica/out/S15_separation_departure.out); SHA-256 `9495fdb12d6da41ff2b0ca5a907ac00fcaaf3ff8f44aa1ca0ca4366afbcf0360`.
- Inventory: 336 Q sections and 336 item-3 object sections. Missing item-3 section names: `[]`. End marker: `['END_OF_TRANSCRIPT: "S15 separation departure Wolfram objects"']`.

Item-3 emitted statuses: 144 closed forms; 192 permitted general-field integrals retained; 0 gaps.

Requested objects not obtained, with the exact method attempted (item-4 gaps include the complete conditions and Q restricted under every such condition):

**item4:** Indicator-polynomial direction elimination; SolveAlways on the separation-domain implication (Reduce for separation-independent conditions) returned unresolved.

- `s/carried/component_i/labels_g` (12 objects).
- `s/transport_centre/component_i/labels_g` (12 objects).
- `s/transport_radius/component_i/labels_g` (12 objects).
- `region/component_i/momentum_labels_g` (6 objects).
- `region/component_i/rate_labels_g` (6 objects).

`s` ranges over inner and outer; `i` over 1, 2 and 3; `g` over `{drain1}`, `{drain1, drain2}`. A pattern denotes every combination. Exact names also appear in this engine’s `EXACT_OBJECT_GAP` index.

## SymPy

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s15/guard/reading_repair_sympy_final_04 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S15_separation_departure.py
```

Guard exit `0`; child outcome `{"exitCode": 0, "guardReason": null, "wallSeconds": 2413.44689261599, "childPid": 1092872}`; `limitsVerified: true`; stderr `0` bytes. Runtime **2413.5108294290258 s**; peak recorded cgroup memory **682344448 bytes**.

- [Script](../../scripts/S15_separation_departure.py); SHA-256 `1a6a839d8425ebdd4a917ccdc77d7a4cbca1c1d2662a04dff1f1a5a7ee265bb4`.
- [Transcript](../../scripts/out/S15_separation_departure.out); SHA-256 `c88240f797d04fd1cda1eb56b75ae90f3632917fd9398bf55f49b8eeccd25013`.
- Inventory: 336 Q sections and 336 item-3 object sections. Missing item-3 section names: `[]`. End marker: `['END_OF_TRANSCRIPT: S15 separation departure SymPy objects']`.

Item-3 emitted statuses: 144 closed forms; 192 permitted general-field integrals retained; 0 gaps.

Requested objects not obtained, with the exact method attempted (item-4 gaps include the complete conditions and Q restricted under every such condition):

**item4:** NOT OBTAINED EXACTLY: Unit-sphere monomial reduction, exact kernel abbreviation retaining all arguments, and solve_undetermined_coeffs on each scalar coefficient were attempted. Kernel identities and exceptional exponent strata remain unclassified; complete conditions and all corresponding restricted Q objects were not obtained.

- `s/carried/component_i/labels_g` (12 objects).
- `s/transport_centre/component_i/labels_g` (12 objects).
- `s/transport_radius/component_i/labels_g` (12 objects).
- `region/component_i/momentum_labels_g` (6 objects).
- `region/component_i/rate_labels_g` (6 objects).

`s` ranges over inner and outer; `i` over 1, 2 and 3; `g` over `drain1`, `drain1+drain2`. A pattern denotes every combination. Exact names also appear in this engine’s `EXACT_OBJECT_GAP` index.

## Scope and interpretation

Item 3 uses computed angular moments, an exact squared-distance substitution, explicit logarithmic exponent branches, and finite-interval Appell integration. Region rates differentiate the computed closed momentum with live bounds and the computed separation/direction rates. General-field integrals permitted by item 3 remain unevaluated. Items 1, 2, 5 and 6 retain the original general-field construction.

Item 4 prints the actual solver operands and returns. SymPy uses exact kernel abbreviations local to each scalar matching attempt, with all argument dependencies retained; the definitions and restoration residuals are printed, and no analytic independence is inferred. Computed conditions and restrictions are printed where reached; unresolved attempts are indexed above. The reviewed record’s claim that an unevaluated zero equation was an obtained condition, with Q restricted to it, is withdrawn.

Each production transcript link was removed before its final run; each engine wrote a new plain file. Runs were sequential, with a 50 GiB cap, zero swap and no deadline. Runtime is the guard wallSeconds; peak includes descendants. Wolfram retains the specified thread options. Progress is emitted and flushed before integration/solve calls. No annex target or repository metadata was modified.
