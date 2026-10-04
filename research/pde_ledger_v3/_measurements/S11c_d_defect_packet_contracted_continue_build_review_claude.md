# NEEDS REVISION

I found no problem that stops the saved-prefix restoration. The blockers are in the new analytic layer, which is the only new science in this build. The AST test lists can pass while the coverage and provenance claims are still not enforced at runtime. I did not execute code, hash files, or calculate new sums.

## Blockers

**B1. The base window and carrier domains are never re-joined.** (`prepare.py:108, 111, 126`)
- `baseplan` is read from `preflight/tail-plan.json`, and nothing requires `(baseplan['K'], baseplan['T'])==(27,122)`. The original check (`original-prepare.py:182`) is only attested through control flow.
- The "both original window domains" check contains the literal `122>=27+4`, which is a tautology. Only the K=29/T=124 side is tested against saved data.
- Neither window's arguments are tied to the saved tail-plan K/T.
- The tail derivation assumes κ<3 and |p0|<3. The new code never requires these for the two `plannedCarriers`. It only emits `physical-plan.json`.
- Fix: require the base plan K/T from saved data. Require the carrier set (`sqrt(595)/10` and `0`) to be joined to the saved physical plan and to the κ<3 and |p0|<3 domain.

**B2. The coverage census cannot catch the errors it is meant to catch.** (`resume.py:19-40`, `prepare.py:143-163`)
- `expected` comes from `canonical_slots(entries,…)`, and `rows` is then built from `expected`. So `check_census(rows,expected)` compares a list with a copy of itself.
- The 20 entries are never joined to the saved preflight tail list. Nothing shows that the J/direct subset of `base['allAddresses']` has exactly these 20 ids. An extra mixed or direct id in the base list would be dropped silently.
- The face census is unenforced. `canonical_slots` accepts any mix of plus and minus, so "both faces present" and the per-face counts are not required.
- `total`, `outer` and `middle` are not in the census key. A mutant that puts K=29 values into a K=27 row, or swaps values between addresses, passes the census.
- The four runtime mutants are drop, duplicate, T+1 and face flip. Wrong K and wrong carrier are exercised only by the synthetic tests.
- Fix:
  - Require the entry ids to equal the saved preflight ids whose component is mixed-iteration or direct and non-zero.
  - Require 10 plus-face and 10 minus-face addresses, or the saved split.
  - Bind each row to its source record. For example, key each row on the saved `outer`/`middle` text for that K, or add a `source` field to the census key.
  - Add runtime mutants for wrong K, wrong carrier, and a swapped-window value.

**B3. The "same formula" premise for the regression and the H identity has no provenance join.** (`prepare.py:131-139`)
- The 20 saved enlarged values came from `original-prepare.py:200-211`. The base values came from the `preflight.py` fragment, lines 325-336.
- The new code compares the preflight fragment to a hand-typed reference expression. It never compares the executed enlarged block to the base fragment.
- The H identity `h·2^-124 == h·2^-122/4` is true for any symbol `h`. It is not bound to any saved value.
- The method accepts this as a dependency. Even so, a no-recomputation AST equality between the enlarged block and the base fragment is cheap and is the only available provenance gate.
- Fix: AST-join the outer/middle assignments in `original-prepare.py` (enlarged `K!=27` branch) to the `preflight.py` fragment, differing only in K and Tlim. Label the H identity as formal only, as the method already does.

**B4. The per-primitive D majorant dependency is not asserted.** (`prepare.py:70-71`)
- `new-per-primitive-absolute-tail-transport` is re-emitted but never checked.
- The saved record carries `eachPrimitiveUsesFullPositiveDEnvelope`. That flag is what justifies using the full D bound, including its outer `Cordinary` term, for each of Dr, Dh and Dq.
- Fix: require that flag, plus the record's presence and ancestry. Disclose that the claim that the individual triangle majorants are dominated by the D outer and middle bounds is an inherited lemma, not a new check.

**B5. The prior tree is not checked for completeness.** (`worker.py:70-78, 90-109`)
- `copy_inputs` copies and hashes only the listed files.
- Nothing asserts that the source tree contains exactly the 8271 files and 158,628,313 bytes. An unlisted extra file, or one deleted after the index was made, would not be noticed.
- Also missing are a count of 33 `new-tail-address-budget-*` records and the fact that the 33rd is the refusal.
- Fix: assert exact tree equality, the file count and the byte total, plus the 32 passing records and 1 refusal.

## Smaller changes before the gate

- Emit all four window and carrier sum records before any `require`. At present a failure at K=27 hides the K=29 records. The two carrier groups are identical copies by construction, so state that they are not independent evidence.
- The gate has a single boolean `independentBuildClearance`. Require both reviewers' literal verdicts in the review record and gate, as the earlier tail-method record did.
- `test_no_completed_function_calls` forbids named functions only. Also forbid any `N.`, `V.` or `G.` call inside `prepare.py`.
- Clarify "numerical identities unchanged". The request-identity construction code is unchanged (`numeric.py`, `request-index.py`, `evidence-store.py`). The identity values change anyway, because the namespace includes the new manifest hash and the new `tails` payload. The payload now has 160 carrier-labelled rows where the original had 80. This is harmless because no SQLite file exists.

## Checks that pass on source reading

- **Numerical run and geometry ASTs:** the `run` body in `worker.py:116-202` and the geometry statements in `geometry-tail.py` match the originals line for line. This is a visual reading of source plus the AST tests. Tests alone are not runtime proof.
- **No completed function is called:** `prepare.py` calls no producer, tail, rule or bank function. The restore step goes through `exact_constant` and `restore_scalar`, which are disclosed structural parser guards, not scientific predicates.
- **Chain and failure checks:** `restore_prefix` checks the chain sequence, hashes and sizes of all 5007 records. It checks the 8360/J/K27 refusal immediately before `failure`, and the failure text against stderr and stdout. It checks that no SQLite file, aggregate or geometry records exist, and that the posthash flags are intact. It requires 243 literal zeros and 20 enlarged observations.
- **Unit joins:** common unit [-2,-1,1] is required per address before any sum.
- **Sum rule:** each D bound is counted three times, J keeps its H overcount, and the comparison is strict `<1e-11` per group. Each group sum is persisted before its decision, and there is no cross-donation.
- **Readiness contract:** `verify_gate` pins the worker, manifest, 97 sources, guard, supervisor, authority, prior index and the earlier literal verdicts. The launcher arms the hook first, with no retry and no deadline.

## Coverage and runtime obligations

- **Gate and run:** a fresh READY gate must follow a clean revision and this review's CLEAR. The science must run once, hook first, under the unchanged guard. Resources are 4 GiB native/cgroup, 16 GiB pool, zero swap, one CPU/thread, 32 tasks, 4 GiB reserve and no deadlines.
- **Records to inspect after the run:**
  - the 8271 copy records and posthashes;
  - the full chain;
  - the 263 restored records;
  - the guard attestation;
  - the 20 unit joins;
  - the four census/sum record sets;
  - the 20 regression records;
  - the H record;
  - the first-failure state.
- **Aggregate:** the new aggregate is unknown. A failure is a preserved stop, not a retune.
- **Numerical stage:** only after the analytic guards pass does it run. Its runtime evidence must show the SQLite FULL chain, the independent A24, A48 and B50 comparisons, and the wrong-root and derivative controls. The 35 stdlib tests are tooling only.
- **Size, memory and time:** the restored context and the 544 address contexts are emitted whole. Their size, RSS and time are unknown and must be recorded.

## What I read

I read in full:
- `review-prompt.md`, `method.md`, `build.md`, `guide.md`;
- `worker.py`, `prepare.py`, `resume.py`, `geometry-tail.py`, `launcher.py`;
- `tooling-tests.py` and `.log`, `original-worker.py`, `original-prepare.py`, `restore-library.py`;
- `execution-authority.json` and `static-preparation.json`;
- `validation.py` and `constructor-reader.py` (function lists only).

Partial reads:
- `input-manifest.json` (scope, source pins, guard-source lists);
- `numeric.py` (grep for `tails`);
- the prior-run `preflight.py` tail section and `prior/complete/new-tail-address-budget-8360-J-27.json`;
- `checks.json` (key fields);
- `evidence-chain.jsonl` (line count 5007);
- the prior file listing.

## Gaps and scoped risks

- **Not independently read:**
  - `prior-result-record.json` (1.8 MB);
  - the rest of the `prior/complete` JSON files, so none of the 263 per-operation records;
  - the per-address and saved large operands;
  - the unchanged `numeric.py`, `geometry.py`, `request-index.py` and `evidence-store.py`;
  - the omitted peer-assessment files, which are opaque and not scientific operands.
- **Hashes:** I could not recompute any hash. I only compared the textual pins in the manifest and `static-preparation.json`.
- **Operands I did not see:** the preflight base `allAddresses` list, which is needed for B2. The saved `Fourier-envelope-constants` bounds, which are the CX/CY inputs, were also not inspected.
- **Inherited and unverified:**
  - the D majorant lemma (B4);
  - the Gaussian identity;
  - the original enlarged-formula arithmetic;
  - that unit [-2,-1,1] applies to the tail magnitudes (this is a declaration, not a derivation);
  - the guard claims for branches that did not run.
- **Pass is not guaranteed:** the preflight stopped on aggregate thirds. Counting each D bound three times can exceed 1e-11 even though the original aggregate passed. The one saved refused row alone is above the 1/80 allocation.

No pending-component budget, full-packet claim, leakage claim or science result is implied. J/direct remains a subset.