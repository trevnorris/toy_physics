**NEEDS REVISION**

This is a source-only review. I read `build-guide.md`, `method.md`, `worker.py`, `launcher.py`, `inputs.json`, `tests.py` and the runtime helpers. I also read the native c2 and composition sources, the five partitions, binding context, physical input, authority, the pressure certificates, and the saved operands in the packet. I executed nothing.

## Substantive blockers

**B1. `worker.py:416` fails deterministically on the saved THETA_BALANCE operands, after all local work is done.**

The check is `require(oldrow['slotCoefficients'][slot]==full['left'] and full['left']==full['right'], ...)`. The second clause is a structural sympy `==`. The saved operands are only equal after `cancel`:

- `saved-inputs/inventory/THETA_BALANCE-delta_p_plus-full-coefficient-input.json` has left `-100*I*epsilon_shape*(3/100 - I/10)/109`. Its right is `-I*epsilon_shape*(3 - 10*I)/109`.
- The saved `...-return.json` records the raw residual and `cancelled` = 0. That is the cancel-based identity the old run used. It is not structural equality.
- The coefficient has more than two Mul factors, so sympy will not distribute it into the same form. `Basic.__eq__` is structural.

U0–U2 pass because every operand there is `Integer(0)`. The failure comes at THETA's first slot, in the loop at lines 410–419. That is after the worker has run:

- all 2,908 local children,
- all cells and the profile jets,
- the denominator joins.

The run then ends `FAILED_PRESERVED`. It never reaches the address joins, the pressure assembly, the controls, or `complete-weak-assembly`, and it uses up the single authorized execution (`scienceExecutionsAuthorized==1`).

Fix: drop the `full['left']==full['right']` clause. The `r['cancelled']` check on the next line already carries that identity. If an operand-level check is wanted, use `zero_record` (cancel) on left minus right.

The same operand-join pattern is untested elsewhere. `tests.py` checks JSON-level equalities only, and nothing in it covers line 416. Three joins fall in the same category and I could not check any of them:
- the structural comparison at `:402`,
- the `prior_manifest['scope']` comparison at `:390`,
- `old_bound['right']==old_affine['left']` at `:292`. The affine input files are in `inputs.json` but not in this packet.

**B2. The pressure ↔ native-partition join rests on labels, counts and past flags, not operand identity.**

What is joined:
- `:270` checks `cancelled==0` on the saved returns.
- `:292` joins raw to consumer raw and bound to consumer bound.
- `:416` joins the slot coefficient to the saved full-coefficient input.
- `:423-439` joins row, face, slot and grade labels, and checks that the 13,260 ids appear once.

What is not joined:
- **Address to partition child.** Every address carries `nativeRowChildHashes`. For THETA address 7956 the hash `be280bff…0607` is exactly the `sha256` of the pressure child at `THETA_BALANCE.json:1052-1053`. The worker never compares these hashes to `partition['children']`. Nothing checks that the 12 pressure children are covered by their row's addresses, or that no U-row address claims a child. "Pressure partition consumed exactly once" is therefore shown only by the id count.
- **Affine right-hand side.** The worker never compares `old_affine['right']` to the sum of `slotCoefficients[slot]` times the slot symbols. The bound-to-slot step depends only on the saved affine zero-return flag.
- **Whole-definition and normal joins.** `:395` compares a flag dictionary (`typedFactorsDistinct: true` and similar). The normal join is a string compare, `normalVariable=='response output l'`. Per-address `normalMultiplier`, `normalOriginal`, `waveMultiplier` and `fullFactorProof.completeNormalMap` are available but unused. The per-address `tagCount` check at `:430-434` is a real join. It supplies the multiplicity evidence for the whole-direct-once and native-iteration-once claims.

Minimum fix:
- Compare each address's `nativeRowChildHashes` to its row's pressure-child `sha256` values. Require that the union over the row's addresses equals the row's `pressureChildIndices`, and that U rows have no such hashes.
- Add one restored-operand identity per row: bound equals the sum of slot coefficients times slots. Use `zero_record`. This is a new join on saved operands, not a replay.

## Nonblocking suggestions

- Run all inherited-pressure operand joins before the 2,908-child derivation, or add stdlib tests that execute them against the saved operands. B1 shows why: today a cheap mismatch costs the whole run.
- `:274` reads the native line in text mode and re-encodes it. The test uses `rb`. Read it in binary to match the test.
- `executionStatus` stays `COMPLETED_…` when `integrityFailure` is set. The exit code and `checks.json` do flag it, so this is ambiguous but not false.
- `exact-evidence.jsonl` does one fsync per record. That is roughly 10⁵ records for the child loop, with no deadline, so the cost is time and disk. Disk size and the 4 GiB address-space margin are unmeasured. This is a limit, not a defect.
- `a['epsilonCount'] in (0,1)` is not tied to the address `status`. This is weak but not wrong.

## What I checked and found faithful

**Partition.**
- `validate_partition` joins each child's bytes to the AST of the full constructor. It checks disjoint, complete coverage and the pressure-name census.
- All 2,908 local children contain exactly one wave symbol, and none contains two.
- None contains `e_W_bg`, `W_bg`, `rho_*_bg`, or any speed or holding-force name. The `k_W` hits are a physical-input parameter.

**Bindings.**
- Extra gamma and material values come from the same physical input. The `omega=3` override is correct. Grades and `c_s0` stay unbound.
- The native `L`-power rule, the `wave_jet` rules and the pressure time-sign rule are joined by AST match against the c2 and composition sources.

**Time sign and denominators.**
- The time sign is checked: `den = scale·I·rho_br·(1 − I·ω·τ)`.
- The `rho_m` memory denominators occur only in the 8 pressure children, never in local ones. So the `rho_br`-only domain check does not trip on local children.
- The other denominator bases are `L_W` and `W_0`, which are constant.

**Epsilon, quotient, cells, controls.**
- Epsilon is degree one once per child.
- The quotient ring uses an exact recurrence with a constant zero-grade denominator, and records the remainders.
- Cells are enumerated over every row, field, order and grade, including zeros.
- The polynomial reconstruction, first-derivative and endpoint identities are exact.
- Control selection and Leibniz placement are sensitive as designed. The Leibniz case reports unavailable if no term applies.

**Launch and failure semantics.**
- The launcher is hook-first, and the coordinator refuses to start without the handshake byte.
- The `containment()` call comes before `import sympy`. The gate, authority and review pins are consistent across worker, launcher and `inputs.json`. Authority scope equals manifest scope.
- The supervisor's `--stage` is free-form, so `defect_full_weak` is accepted.
- Batch prefix and failed-child preservation, the hash-chained evidence log, and the final receipt and posthash paths match the stated semantics.

## Limits of this review

- Not verifiable from this packet:
  - the affine input and return files,
  - the 13,260-address arrays and their `tagCount`s,
  - the prior weak run's `input-manifest.json` scope strings,
  - the stdlib test results,
  - file hashes.
- The mathematical bounds (inequality (2), the derivative recurrence, Schwartz continuity) are method-level analysis, which the method itself says is not machine-checked.
- Local endpoints are not full nonlocal ends.
- This verdict is not a computed result.