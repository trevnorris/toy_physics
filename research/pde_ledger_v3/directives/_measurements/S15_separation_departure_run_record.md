# S15 separation departure — item-4 repair3 run record

Contract inputs: Part A (`91248b67`) and Part B (`7f52b4d2`); reviewed build supplied for this repair: `4c4e2558`. These are provenance supplied in the brief, not commit inspections. Artifact paths below are relative to `research/pde_ledger_v3/` unless absolute. Each engine reconstructs its objects independently. No cross-engine output comparison or review was performed. No git, git-annex, or datalad command was run.

This record supersedes the previous record's claims that finite-identity equations solved all requested conditions and supplied Q under every condition. Those claims were unsupported. The repaired engines print parameter branches obtained by exact attempts and their restricted Q, with explicit not-obtained labels for unresolved classifications. The repair does not establish all requested item-4 conditions.

## Wolfram

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s15/guard/repair3_wolfram_final_01 -- math -script /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S15_separation_departure.wl
```

Guard exit `0`; child exit `0`; guard reason `null`; `limitsVerified: true`; stderr `0` bytes. Guard wall time **653.0475128099788 s**; peak recorded cgroup memory **3055435776 bytes**, including descendants.

- Source: `mathematica/S15_separation_departure.wl`; SHA-256 `952923eda122dd92ecf5d9882f377696bac6fc5057df0de5c6881f53edd62fc8`.
- Transcript: `mathematica/out/S15_separation_departure.out`; SHA-256 `1dbe97d8e187fd550bbb79aa962dbfb4322dc045168a88a1d970eebccfe1ed6f`.
- The 336 item-4 statuses and transcript end marker are present. There are 210 numbered sufficient conditions and 210 corresponding Q restrictions across the 48 separation-dependent objects; none of these entries is `Missing`.

Missing requested objects: complete item-4 classifications, including all exceptional exponent strata, and Q under every condition, for all 48 separation-dependent objects. These are the three Cartesian components of each drain-1 and mixed drain-1/drain-2 group of carried current, centre transport and radius transport at each surface, and of region momentum and rate. Each object's `/item4/status` and final `EXACT_OBJECT_GAP` identify this omission.

Attempt: direction-polynomial reduction of the printed closed form, native `SolveAlways` on the separation-domain implication, and `Reduce` of computed physical-monomial equations. Native identity-solver returns remain candidates or unevaluated returns. The monomial construction yields sufficient branches only. Each printed branch is substituted into all scalar identities; their operands and residuals are printed. Q is separately substituted and simplified on that branch. No completeness follows from these sufficient branches.

## SymPy

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s15/guard/repair3_sympy_final_02 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S15_separation_departure.py
```

Guard exit `0`; child exit `0`; guard reason `null`; `limitsVerified: true`; stderr `0` bytes. Guard wall time **5262.722074120073 s**; peak recorded cgroup memory **808783872 bytes**, including descendants.

- Source: `scripts/S15_separation_departure.py`; SHA-256 `472cf2430877bd3f5d6216f7118c568a54d3f0d23532447aa769bc461bf907bb`. This matches the source pinned for the final run.
- Transcript: `scripts/out/S15_separation_departure.out`; SHA-256 `cf4267b547fd92ed46165b6ae150f3b8f5e87f375a07e43f300c74a01f7e829d`.
- The 336 item-4 statuses and transcript end marker are present. The 48 separation-dependent objects have 782 numbered parameter branches and 782 corresponding Q restrictions. These are printed branches, which can overlap; they are not a count of distinct solution components. Across all item-4 objects there are 890 numbered conditions and 890 corresponding Q restrictions.
- Nine separation-dependent objects have no unresolved branch entry: the three Cartesian components of the outer surface's drain-1 carried current, centre transport and radius transport. Their primitive-support partitions, necessary factor zero sets, full-identity residuals and Q substitutions are printed. This inventory observation is not an additional completeness theorem.

Missing requested objects: complete item-4 classifications and Q under every condition for the following **39** objects (each row includes all three Cartesian components):

| Objects | Groups | Count |
|---|---|---:|
| Inner carried current, centre transport, radius transport | drain 1 and mixed drain 1/drain 2 | 18 |
| Outer carried current, centre transport, radius transport | mixed drain 1/drain 2 | 9 |
| Region momentum and rate | drain 1 and mixed drain 1/drain 2 | 12 |

The exact unresolved parameter branches are indexed beside each object's `/item4/unresolved_parameter_branches`, with the attempted construction and its return. Final `EXACT_OBJECT_GAP` lines index all 39 objects. Conditions and Q are printed for the obtained branches; those do not fill the unclassified alternatives.

Attempt: a computed affine exponent coordinate, partition by the closed primitive's exceptional supports, differential determinant equations used only as necessary candidate equations, exact factor zero sets, and substitution into the full scalar identities. The fallback recursively expands in instantaneous polynomial parameters while retaining exact analytic coefficient definitions, then applies `solve_undetermined_coeffs` to a stronger formal-coefficient identity. Its returned branches are tested against every original scalar identity, and Q is computed on each verified branch. The fallback supplies sufficient conditions only; analytic coefficient cancellations and exceptional exponent alternatives remain unclassified where labelled. Amplitude values and time derivatives remain independent instantaneous slots.

The first full SymPy attempt (`repair3_sympy_final_01`) was stopped by the author through its live guard after a prolonged silent direct coefficient-matching call. Its unit became inactive and its reservation was released. The fallback was changed, its guarded probe completed, and the final run above was executed from the revised source. This was not an external termination.

## Scope and validation

The old finite-jet condition and reference-point Euler template were removed from item 4 in both engines. No finite multiplicity bound or equivalence is claimed. Intermediate necessary equations are labelled as intermediate candidates rather than parameter conditions.

The final outputs are checked separately against each engine's own pre-repair transcript. The audit excludes item-4 lines, item4-kernel lines, progress and final gap-index lines; it hashes every remaining line in order, including the outputs of items 1, 2, 3, 5 and 6. This is a preservation check, not a cross-engine comparison or a proof of the new item-4 classifications.

Wolfram: all **20438** retained lines match its own baseline, SHA-256 `430af6f5089e3c6b78ed13efb557637e39aad679e01882bc11250ed09e9e3571`.

SymPy: **18602** of **18608** retained lines match byte-for-byte. The six differing lines are `reading/outer/pressure/component_i/whole_minus_expansion/E_geo` and `/E_resp`, for `i = 1, 2, 3`. Exact parsing and bound-variable canonicalization alone retained differently factored integrands. Expanding the held integrands, without evaluating an integral, gave the literal residuals below. Thus these six differences are changes of written form, not value. The raw retained-output hashes are intentionally not reported as equal: baseline `c498b411fed7b180a345a3857e298ec59d9899156127fe54b2e1c191c705ec53`; final `835e957c8a8db08fe98b94368ea3b8e999b8faa438950073b6e25d78f00138a6`.

The additional same-engine preservation check used:

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s15/guard/repair3_sympy_preservation_02 -- python3 /var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/verify_sympy_retained.py
```

Guard and child exits `0`; `limitsVerified: true`; stderr `0` bytes; wall time `9.089649471919984 s`. The script SHA-256 is `673216b8d643d2b562178dde349ced20a95e6f2214fa7bb742701467d30fe58f`; guard `stdout` SHA-256 is `d7953eb31afdbec6437ea2a104f58d470c85356fc7c43f5fd4e9c71831d8cc19`. It prints both baseline and final operands before these residuals:

```text
reading/outer/pressure/component_1/whole_minus_expansion/E_geo/expanded_integrand_residual: 0
reading/outer/pressure/component_1/whole_minus_expansion/E_resp/expanded_integrand_residual: 0
reading/outer/pressure/component_2/whole_minus_expansion/E_geo/expanded_integrand_residual: 0
reading/outer/pressure/component_2/whole_minus_expansion/E_resp/expanded_integrand_residual: 0
reading/outer/pressure/component_3/whole_minus_expansion/E_geo/expanded_integrand_residual: 0
reading/outer/pressure/component_3/whole_minus_expansion/E_resp/expanded_integrand_residual: 0
```

The final engine sources were not changed after their reported final runs. The preservation check changed neither source nor transcript.

Both transcript links were removed before their respective runs, and the runs write new plain files. Runs are sequential with the specified 50 GiB cap, zero swap and no deadline. Wolfram retains both thread-limiting options. No annex target or repository metadata was modified. Source differences are confined to item 4. No commits or review verdicts are created by this record.
