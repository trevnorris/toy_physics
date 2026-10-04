Source and JSON review only; nothing was run. The verdict is clear: I found no substantive mathematics or validation blocker, and only local tooling and wording notes.

CLEAR FOR THIS BOUNDED PACKET-GEOMETRY BUILD

## What I inspected

- **Documents:** `build.md`, `method.md`, `packet-action-method.md`.
- **Code:** all of `worker.py` and `library.py`, and all of `launcher.py`.
- **Runtime sources:** `runtime-source/supervisor.py` in full. For `shared-guard.py` I read only keyword-grep hits for the deadline, restart and memory properties. For `completion-hook.py` I read only grep hits for its arguments and `waiting` state. I did not read `AGENTS.md`.
- **Tests:** `tooling-tests.py`, `tooling-test-record.json` and `tooling-tests.log`.
- **Manifest and authority:** `execution-authority.json` and `inputs.json`. In the manifest I read the scope, the resources, the `preflight`/`rules`/`accepted` receipts, and `sourcePins` up to about line 1494.
- **Saved operands read in full:**
  - physical input and `local/context.json` (partly),
  - `preflight/checks.json` and `preflight/physical-plan.json`,
  - native scale join,
  - `B-G7-K15` rule,
  - one nonconstant and one constant field triple,
  - one minus-normal flat adapter input and arguments,
  - factor-12 normal and residual inputs,
  - the first rows of `fields.json` and `address-adapter-result.json`,
  - the first rows of `pressure-addresses.json`.
- **Grep counts over all 544 addresses:**
  - 102 live, 106 zero-source and 336 zero-consumer; no other status exists.
  - All five component names appear.
  - 288 flat addresses have `q(l)` input depth and `l-p` transfer; 256 non-flat have `q(k)` and `k-p`.
  - All 544 have `r-l` transfer and `q(l)` output depth.
  - The flat response-map sha is identical for the plus and minus faces (288 hits), so the worker's omission of `id` from the join is correct.
  - Proof 0 serves 144 addresses and proof 12 serves 72.
- **Coverage limits:**
  - I did not open every one of the 266 files.
  - The remaining field triples, adapter files, and factor files other than 0, 12 and 14 were not read.
  - I did not hash anything.
  - The gate and build-review record do not exist yet, so the exact gate key names are unverified.

## Evidence

- **Arithmetic and source bridge (`library.py:22-78`, `worker.py:264-278`).**
  - The Q(√(119/20)) arithmetic is exact, and the mixed-sign comparison by squares is correct.
  - The saved `kappa` string and the rational constant 595/100 both give κ² = 119/20; this also equals 6 − 1/20.
  - Float, bool and non-square-field inputs are refused.
- **Saved shapes.** Every key the worker reads matches the saved files I opened:
  - `coefficientId`, `transfer`, `fullFactorProof.proof`, `completeNormalMap`;
  - `actualResponseMap` and its eight joined keys, `requiredMap`, `addressNormalOriginal`, `savedNormal`, `mappedAddressFactor`;
  - `definitions` with `mapped`/`template`/`firstAddressId`/`original`, `addressJoins` with `adapter`, `arguments` with `address`/`actualCompleteFactor`;
  - the rule keys and the `frequencyOverride` shape.
- **Cancel proofs.** They are joined on `left`, `right` or `cancelled` only, with no left-equals-right requirement, as intended.
- **Normal signs.** The plus and minus normal signs enter through the complete numeric adapter and factor joins (`I*q` against `-I*q`).
- **Square arrangement (`library.py:108-192`, `195-213`).**
  - The line set matches the method: boundaries, ±κ, `l=-k+{0,±2}κ`, `l=k+d`, carrier horizontals.
  - Coalescing merges labels and requires unique labels.
  - All pairwise crossings are saved, in-box ones cut the k axis, and slab midpoint ordering is rechecked independently.
  - The audit demands nonnegative endpoint widths and positive midpoint widths and areas.
  - Box membership changes only at box-line crossings, which are cuts, so the argument is sound.
- **Height domain.**
  - It has `Q=0`, `U` and positive profile offsets.
  - `Q=±(v-k)` is built for ±κ and for every carrier-offset `v`, covering both consumer branches.
  - All crossings are included.
- **Controls.**
  - The missing `collision:sum:2` mutant still covers the box and is refused with exactly `required collision/resolution incidence`.
  - The reversed cell is refused with `oriented adjacent cell`.
  - Operands are saved before each refusal.
- **Counts (`library.py:216-224`, `worker.py:219-248`).**
  - (2n)² per cell and 2n per slab for A24/A48, and a 15² lower bound for B.
  - Flat and contact each get a single outer integral over the square slabs.
  - Flat, height-contact, height paired PV, slope, H, J and D are separate records, with per-grade address lists.
  - The 442 zero addresses stay in the dependency table and are not scheduled.
  - Families are interned on exact canonical JSON; the hash is only a label.
- **Claims.**
  - No false readiness, precision, normalization or unit claim appears; units are explicitly `PENDING`.
  - Old-bank reuse is a constant 0, never a computed match.
  - `numericalEvaluatorReady` is False.
- **Persistence and containment.**
  - Inputs are emitted before their guards in the field, numeric, address and factor loops, and the arrangement is emitted before the audit.
  - Each artifact is chained.
  - Failure records the active operation and completed list.
  - Posthashes run in `finally`.
  - Pins and argv are checked before the output directory is created.
  - The launcher arms the hook through a pipe handshake before the science starts.
  - The guard has `RuntimeMaxUSec=infinity`, `Restart=no` and no deadline, and `containment()` checks 4 GiB, zero swap, one CPU/thread, 32 tasks and no CPU rlimit.
  - The tests are synthetic and metadata-only, and the test record's worker and library hashes equal the manifest pins.

## Blockers

None.

## Local tooling fixes (not blocking)

1. `worker.py:277-278`: `geometry-basis` is emitted after the κ² guard. Swap the order.
2. `library.py:148`: `audit` takes its box from the plan under test. y is tied to K through the required box lines, but the x extent is not. The worker passes the spec box, so this is consistent. Add an explicit K-box check at `worker.py:288`.
3. `worker.py:240`: only the Y branch count (2) is flagged for height PV. The response is also evaluated at `k+Q` and `k-Q`, and `addressProductOccurrences` and the conditional storage sum are not scaled. Counts are per outer occurrence and understate PV records by up to 2× unless the multiplicity is applied. Add a response multiplicity field.
4. `worker.py:205`: `normalMultiplier` is not directly joined to `savedNormal`. The sign is enforced only through the complete factor and adapter equality.
5. `worker.py:150-152`: the 34 field joins check `left` and the zero return. `polynomial.json` is preserved but not cross-joined to the proof's `right` side, which would need CAS decoding, so this is a documented limit.
6. `worker.py:140-143`: there is no guard that `len(full)` is 2652, and none that the target grade equals the sum of the grade triple.
7. `worker.py:201-204`: the normalization strings are hard-coded. They match the saved field definition in the entries I read; add a join to `fields[fid]['definition']`.

## Optional text

`oldBankRequestMatches: 0` (`worker.py:215`) and the launcher message "Actual old-bank matches are zero" could read as though a match was attempted. None was.

## Runtime obligations

- Inspect stdout `checks.json`, strict stderr, guard and supervisor records, the evidence chain, the copy receipts and all posthashes. Exit code is not acceptance.
- Confirm each plan's box is [-K,K,-K,K] for the square and [-K,K,0,U] for the height domain, with K of 27 or 29 and U of 75.
- Confirm all 8 audits pass, with coverage, cell and slab counts and expected area equal to the box area.
- Confirm both controls refuse with the exact messages.
- Check that the 102/106/336 split is carried into the dependency table, and report per-grade and per-component emptiness explicitly.
- Apply the PV ×2 multiplicity by hand when reading storage counts.
- Treat all counts as static lower bounds or occurrences, not totals.
- Run the first actual read of the gate keys with a fresh pin check. A refusal there is a refusal, not clearance.

## Scope

This clears only the geometry, dependency and count build. It does not clear pressure unit certificates (still pending), any numerical route, adaptive cost, bytes or RSS, or runtime behavior.