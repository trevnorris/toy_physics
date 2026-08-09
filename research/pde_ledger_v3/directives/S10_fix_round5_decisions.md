# S9/S10 export chain — round 5. **A REMOVAL ROUND.** Decision list.

**You author the directive prose and apply the change.** Read `CLAUDE.md` first.

Scope: both SymPy engines · `scripts/extract_knob_inventory.py` · both generated modules · both
regenerated `scripts/out/*.out` · a new directive under `directives/`.
⛔ Out of scope: every `.wl`, every step record, every `.tex`. ⚠ S10 is slow: wrap runs in `timeout 1800`.

## ⛔⛔ THIS ROUND ADDS NOTHING. IT DELETES, IT FIXES ONE REGRESSION, AND IT RECORDS.

⛔⛔ **Do not add a single new check, residual, guard, count or emitted tag.** ⚠ Four rounds have now been
run and the orchestrator commissioned **five** in-run checks on the export writer against a standing rule
that says the export writer cannot verify itself. **Two are provably inert and one scores a perfect zero
when the thing it checks is deleted entirely.** ⭐ If you believe something here needs a new check,
**name it and stop** — ⛔ do not build it.

⭐ **The physics is untouched and stays untouched.** Every value in both transcripts is correct and four
legs have confirmed it. ⛔ Nothing below touches a derivation.

---

## R1 · Delete the two residuals that cannot be nonzero.

⭐ **`value_kind`** — `S9:904` writes the field with `export_value_kind(value)`; `S9:1116` checks it with
**the same function**. Measured: falsifying that function to return one constant leaves every residual 0
and publishes. ⇒ writer and checker are one source.

⭐ **`BUILD_INPUT_DIGEST_RESIDUAL`** — compares a `str→str` dict against itself across a `repr`→`eval`
round trip. Measured operand identity `True`.

**Decision: both residuals go.** ⛔ Do not replace them. ⭐ **Keep what they were attached to** — the
`value_kind` **field** and the `BUILD_INPUT_DIGESTS` **mapping** are both load-bearing and both work; a leg
used the digests to detect a stale module. ⚠ It is only the self-comparisons that are decoration.

## R2 · Importing an engine must not destroy the ledger.

⚠ **A regression this round's predecessor introduced.** The `unlink` moved to module scope (`S9:8`,
`S10:10`), so **merely importing either engine deletes the corresponding generated module**, even when the
import then fails. Measured: `S9_exports.py present before=True after=False`.

⛔ Nothing imports these engines today, so this is latent — ⭐ but a comparator, inventory scanner or doc
generator that imports one silently destroys the ledger, and a comparator is the next thing being built.

**Decision: importing an engine has no effect on any file.**

⚠ ⭐ **Keep the property the move was making true.** Measured by weakest-change control: putting the
`unlink` back inside `main()` makes the chain's own failure mode — upstream missing — leave a stale module
again. ⛔ **Do not simply revert it.** Both properties must hold at once; ⭐ if they cannot, **name that and
stop.**

## R3 · Record these as measured limits, in your directive. ⛔ Do NOT fix them.

⭐ Each is a real limit on an instrument, and ⛔ every repair anyone has proposed for them is either a
denylist or another self-comparison.

1. ⛔ **`value_kind` marks one carrier.** A sentence reaches the ledger tagged `COMPUTED_OBJECT` as
   `Symbol('…')`, `Function('…')(t)`, `Tuple(Str('…'))`, `Eq(Str,Str)`, or in the `display` field — which
   is a raw `str` on **every** record and is in no residual. ⚠ A `Symbol` whose name is a sentence gets an
   auto-created ledger record **whose key is the sentence.**
2. ⛔ **The authored-word count does not bound the population.** Emitted **24**; records carrying an
   authored word anywhere **36**; distinct authored strings **104**.
3. ⛔ **The assumption-channel residual is coverage-free.** Emptying the entire `Q` channel scores **0** and
   publishes. Strengthening a `Symbol` no `Q` mentions scores 0 and is inherited downstream. Coverage is
   **11 of 22** S9 symbols, **14 of 45** S10 symbols, `MAIN` packages only.
4. ⛔ **A parse failure leaves a stale module** — Python compiles before executing, so no module-scope
   statement can prevent it. ⭐ It is **detectable** via `BUILD_INPUT_DIGESTS`.
5. ⛔ **Every export guard is an `assert`** ⇒ under `PYTHONOPTIMIZE=1` a violated guard publishes.
6. ⛔ **`BUILD_INPUT_DIGESTS` omits the CAS** — no sympy or Python version, so two modules with identical
   digests can hold different values if the CAS changed underneath them.

⭐ **State plainly which of these a consumer must therefore check for itself.**

---

## Constraints

- ⛔ The derivation, action, ansatz and every computed physics value stay untouched.
- ⭐ Report the complete `.out` diff for both engines. ⭐ **Expect it to be deletions and nothing else** —
  ⚠ if anything else moves, that is a finding and you should report it rather than absorb it.
- ⭐ Every S9 record S10 does not overwrite stays identical between the two generated modules.
- ⛔ Do not state anywhere what a count, tally or partition size comes out as. Emit it and let it be read.

## Acceptance

⭐ Show, by running them, that the two deleted residuals could not have been nonzero — ⛔ not by argument.
⭐ For R2, demonstrate **both** properties at once: import destroys nothing, **and** the upstream-missing
failure leaves no stale module.

## Deliverables

The changed files, both literal stdouts, all diffs, and every ablation script with its stdout at named
absolute paths **outside the repository**.
