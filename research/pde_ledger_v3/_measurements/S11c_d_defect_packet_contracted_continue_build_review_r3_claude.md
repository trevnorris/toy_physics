I'm midway through the review. I've read the method, build notes, guide and all the new sources. I just checked the base preflight source to confirm where `cx` and `cy` come from. Now I'm writing up the verdict.

**CLEAR FOR THIS SAVED-PREFIX CONTRACTED NUMERICAL CONTINUATION BUILD**

This is a Claude-only, source-only verdict under the temporary user policy. It is not paired clearance, a runtime outcome, or a READY gate. The new aggregate is unknown, not a pass, and the original refusal at 8360/J/K27 stays a failure. J/direct is not full-packet or leakage. I executed nothing.

## What I read
- **Named files:** `method.md`, `build.md`, `guide.md`, `worker.py`, `prepare.py`, `resume.py`.
- **New and tooling sources:** `geometry-tail.py`, `launcher.py`, `tooling-tests.py`, `tooling-tests.log`.
- **Policy and authority:** `review-policy.json`, `policy-transition.json`, `execution-authority.json`, `static-preparation.json`.
- **Libraries:** `restore-library.py` and `constructor-reader.py`.
- **Original source:** `original-prepare.py`, lines 170–229 in full plus targeted greps.
- **Base preflight:** `prior/source/.../S11c_d_defect_packet_preflight.py`, by grep only.
- **Manifest and numeric:** `input-manifest.json` by grep only (source pins, attested guard sources, resources, `continuationSourceProof`), and `numeric.py` by grep only.

## Coverage gaps
I did not read these, so I make no claims about them:
- `numeric.py` beyond the grep, plus `geometry.py`, `request-index.py`, `evidence-store.py`, `validation.py` and `original-worker.py`.
- The 8,271 prior files and the per-address JSON operands.
- Most of `input-manifest.json` (about 49k lines), and `baseline-numerical-method.md` and `review-prompt.md`.

I could not check that the unchanged-engine claims, the run AST equality or the pinned hashes hold. I only confirmed that the code checks them. Hashes are taken from the packet's own labels.

## Assessment
- **Subset analytic rule:** It is separate from the numerical epsilon of 1e-11/80, which stays untouched. The ceiling is strict (`value<1e-11`) and per window and carrier (`prepare.py:240-242`). It covers 40 primitive rows, each D counted three times, J keeps its H overcount, and nothing is donated or cancelled.
- **Rows and units:** `canonical_slots` enforces 20 entries, 10 per face, 5 J and 5 D per face, and units `['-2','-1','1']` before any sum.
- **Census joins:**
  - `independent_slots` builds the expected rows from the full preflight list, the native rows and the published enlarged returns.
  - Those rows are joined by `check_census`, including the bound source and the encoded outer and middle operands.
  - Windows are cross-checked against the saved plan and the actual enlarged inputs.
- **Scalars:** `checked_totals` re-materializes the validated operands and returns the totals used by the sums, so the sums never consume mutable row totals.
- **Mutants:** All eleven mutants per group go through that same predicate: drop, duplicate, wrong T, K, carrier and face, swapped window values, swapped address values, and total-only, outer-only and middle-only. A refusal is required for each.
- **Ordering:** All four sum records are emitted before the final guard (`prepare.py:242-246`), and the 33 budget checkpoints are restored in chain order with the 33rd labelled refused.
- **Original code unchanged:** I found no call to a completed producer, prepare, identity, tail, rule or bank function in the new `prepare.py`, `resume.py` or `geometry-tail.py`. The runtime `continuation_source_joins` compares the original `run` and unfinished geometry tail by AST.
- **Formula joins:**
  - `tail_formula_join` compares the base and enlarged J, H and D ASTs, allowing only the two storage-variable renames.
  - `tail_binding_joins` joins `b`, `Cordinary`, `E30`, `F` and both moment helpers.
  - My own reading of the enlarged branch in `original-prepare.py:199-211` matches these pins.
- **Dependencies disclosed, not independent:**
  - The enlarged closed-form arithmetic, `weighted_tail`, `exponential_moment`, and the D triangle and shift-root lemma are inherited dependencies.
  - Restore-parser guards are structurally re-run, and the code labels that.
  - The enlarged-versus-base comparisons are argument-level regressions, not independent evidence.
  - The two carrier groups reuse the same uniform envelope.
- **Policy and gate:** The AGENTS exception is narrow: before and after hashes, an exact appendix, and old bytes joined from the frozen copy. The gate requires the literal Claude verdict string.

## Remaining scoped risks and obligations
None of these blocks the build, but each is owed at runtime or in the record:

1. **H relation is tautological.**
   - **Issue:** `prepare.py:206` is the identity h·2⁻¹²⁴ = h·2⁻¹²²/4 with hard-coded exponents and a free symbol.
   - **Why it's acceptable:** Its force comes only from the H-term AST equality at line 203, the T joins at line 195, and the base/enlarged join. It is disclosed as formal.
   - **Do not cite it as:** An independent bound on H.
2. **Base `cx` and `cy` are not AST-joined in code.**
   - **Issue:** `tail_binding_joins` does not pin `cx,cy=bounds[...]['X'],['Y']` on the base side.
   - **What I confirmed:** The base preflight at line 325 reads those same keys, and the saved enlarged `CX` and `CY` are compared to the envelope.
   - **Obligation:** Record this as inspected by reading only. Tighten it if the code is ever revised.
3. **Worker and module identity.**
   - **Issue:** The AST checks read pinned file paths, and `resume` and `tail` are imported by name through `sys.path`. Nothing asserts that the loaded modules' `__file__` equal the pinned paths.
   - **Obligation:** Treat the unchanged-AST claim as joined to pinned bytes. Add a `__file__` equality check if the code is revised.
4. **Spurious stops are possible, but a false pass is not.**
   - Operand re-encoding must round-trip against saved `srepr`.
   - The swapped-window and swapped-address mutants must actually change an operand.
   - A stop here is a stop, not a pass, and must be preserved as such.
5. **Tails list changed shape.** It now has carrier-tagged rows (160 versus the original 80). The unchanged `run` only stores it. I did not read `numeric.py` to rule out other consumers.
6. **Carrier duplication.** The two carrier groups are numerically identical by construction. They give no extra evidence and must not be reported as two independent passes.
7. **Tooling tests.** The 72 tests are stdlib, synthetic or AST only. Their pinned log is not scientific or runtime proof.

## Runtime obligations
- **Gate:** Create a fresh READY gate whose hashes match this reviewed worker, manifest and sources, and record this literal verdict. Launch hook-first under the unchanged guard: 4 GiB native/cgroup, 16 GiB pool, zero swap, one CPU/thread, 32 tasks, 4 GiB host reserve, no deadline, no retry.
- **Preserve and inspect, whatever the outcome:**
  - Full prior-tree and posthash receipts.
  - All 263 restored operation records.
  - The 33 budget checkpoints.
  - The four group sum records, the component regressions and the eleven refusals per group.
- **Interpretation of results:** If any analytic group fails, preserve it and stop, with no retuning. After a pass, SQLite FULL, all numerical identities, the controls and the A24/A48/B50 comparisons still require runtime evidence.