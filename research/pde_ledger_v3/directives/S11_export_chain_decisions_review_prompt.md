# Independent review — the export-chain decision list, before any builder launches

## Artifact under review

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_export_chain_decisions.md`

It is **orchestrator-written**. Nothing in it has been applied. You are one of two independent legs; the
other is not visible to you, and neither of you sees the other's report.

⚠ This list decides a **permanent convention** for an export chain that every later step consumes. A wrong
rule here corrupts artifacts that will be built for months and will not announce itself.

## ⭐⭐ REQUIRED READING ORDER — ⛔ do not read the decision list first

Reading it first anchors you to its framing, and its framing is the thing under test.

1. `research/pde_ledger_v3/S9_REWRITE_PLAN.md` — the chain design; `§3`, `§4`, and `D1`–`D12`.
2. `research/pde_ledger_v3/scripts/S10_exports.py` — ⭐ **import and inspect it programmatically.** It is
   4401 generated lines; ⛔ do not read it linearly.
3. `research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py` — the writer. `§write_exports`
   and `generated_ledger_key` are the relevant parts.
4. `directives/S10_SHARED_PHYSICS.md` `§8`, and `directives/S11_SHARED_PHYSICS.md` `§Q6r` and `§8`.
5. `research/pde_ledger_v3/DEFECT_REGISTER.md`, entries `C19` and `C20`.
6. **Write down, before step 7:** how *you* would key a chained export so that two steps computing
   different objects cannot silently overwrite one another, and what you would do about the three
   overwritten rows. ⭐ Keep that list; you will be asked to compare it at the end.
7. **Only now**, `directives/S11_naming_and_chain_plan.md`, then the decision list under review.

## ⛔⛔ FOUR UPSTREAM CLAIMS ARE IN SCOPE — ⛔ they are NOT settled background

The decision list rests on a plan and on four adjudications, **none of which has had independent
scrutiny**. Two prior legs reviewed an *earlier* version of that plan whose recommendation was the
**opposite** of what is now proposed. ⭐ Treat all four as open questions of this review.

### U1 · ⛔ Producer-scoping was proposed by a reviewer and reviewed by nobody

The `E2` shape — producer-scope new keys, carry imported rows unchanged — originated **inside a review**,
so exactly one party has ever considered it, and that party wrote it.
⭐ Is it the right shape? ⚠ Does it **destroy the cross-step contention the chain exists for**? Under `D5`,
a later step re-deriving an earlier object and writing the same key is what makes disagreement visible at
all; if every new key is producer-scoped, ⛔ can two steps still be made to contend over one object, and
by what mechanism? Answer with a concrete route through `S10_exports.LEDGER`.

### U2 · ⛔ The ordering argument is the orchestrator's, invented after the reviews

The plan says S10 must be re-keyed **before** S11 is built, because *"S11's PY engine binds to S10's export
keys."* The decision list now claims measurement refutes this (`M1`–`M4`) and `E3` leaves S10 alone.
⭐ **Settle it by running code, not by reasoning.** Establish S11's actual binding surface into
`S10_exports.LEDGER` from `S11_SHARED_PHYSICS.md` `§Q6r` — every name it spells, every pointer it
dereferences — and report what a producer-scoping of S10's 455 own keys would and would not move.
⛔ If anything S11 needs would break, `E3` is wrong and the retrofit is mandatory; ⭐ say so and name the row.

### U3 · ⛔ The "rule 2 violation in the chain" claim has been reversed, and the reversal is unreviewed

The plan asserted that `S9_REWRITE_PLAN`'s overwrite mechanism *destroys one operand and leaves cross-step
agreement merely asserted* — a violation of the project's rule that a script emits **both operands and the
residual, then guards.** `E5` now says the opposite, on the basis that the engine emits the operands and
the predecessor value is committed upstream.
⭐ **Adjudicate it against the artifacts.** Open what the engine emits for the overwritten rows and what
the committed output contains. ⚠ Then ask the harder question the reversal does **not** address: a
downstream consumer reads only `corroborated_steps`, a two-element tuple, from the merged file.
⛔ **Is that a typed claim of agreement with no operands behind it, inside the artifact that carries it?**
⭐ If a consumer cannot check it without opening two other files, say whether that is acceptable and what
the alternative costs.

### U4 · ⛔ A leg was overruled on a count, and the rebuttal was never shown to it

Two legs disputed how many keys would collide between S10 and S11 under a unified vocabulary. One said
**1**, one said **8**; the orchestrator adjudicated in favour of **8** and neither leg saw the argument.
The losing side's second point was never disposed of: the true post-retrofit count is **undecidable**
without the rename map, because `S11_SHARED_PHYSICS.md:940-942` permits an engine to name an object the
spec does not.
⭐ Re-run the count yourself and report the number with **the proxy you used**, explicitly.
⛔ Does the undecidability objection survive `E2` and `E4`, or do those make the count irrelevant? ⭐ If
irrelevant, say what is then load-bearing instead.

## What to check in the list itself

### 1. ⭐⭐ Are `M1`–`M13` reproducible, and does any of them fail to support the decision it is cited for?

⭐ Re-run them. ⛔ A measurement that is correct but does not entail its conclusion is the failure mode to
hunt: report the gap between the number and the claim it is used to justify.

### 2. ⭐⭐ Is `E2`'s boundary drawable by a mechanical rule?

`E2` splits keys into *created by this step* and *imported and untouched*. ⚠ The list does **not** choose
the scope syntax and does not state the rule.
⭐ **State one, and test it against all 617 entries of `S10_exports.LEDGER`.** Report how many fall each
side and whether any entry needs a human judgement. ⛔ If it does, `E2` defers its hardest problem into a
builder's hands, and you should say which entries.
⚠ Specifically: the three overwritten rows are *created by S10* **and** *carried from S9*. Which side?

### 3. ⭐ Does `E4` actually fire on the case it exists for?

Construct the collision it is meant to catch — a later step writing a key an earlier step already used for
a **different object** — and check that `E4` as worded catches it and that `D5`'s existing chain-integrity
check does not. ⛔ If `E4` also fires on a legitimate case, name it.

### 4. ⭐ Is `E1` decidable from the specs, or is it being invented here?

⚠ S11's shared spec appears not to require an export at all. ⭐ Does anything **oblige** S11 to write one,
and what breaks at S12 if it does not?

### 5. ⭐ `E6` and the naming question

The list claims S10's engines did not deviate from their spec, because that spec supplies a tag *grammar*
and no quantity vocabulary. ⭐ Verify that. ⚠ Then: does a spec that fixes a grammar but not the names
satisfy `D11`? ⛔ Can you construct a case where two engines both obey the grammar, agree on a name, and a
**wrong claim** reaches a step record or a **false agreement** reaches the comparator?

### 6. ⛔ Does the list leak, or commission damage?

⛔ Does anything in it state or imply what a computation will produce? ⚠ **The test is derivability, not
literal presence:** a stated justification from which a builder could *derive* a value is a leak even
though no value appears.
⚠ Does any decision touch a correct artifact more broadly than its defect requires? ⭐ Name what would be
damaged.

### 7. ⭐ What is missing?

⭐ Compare your step-6 list. What must be settled before a builder starts that this list leaves open?
⚠ Its *"does not decide"* section is a deliberate exclusion — ⛔ is anything in it actually load-bearing
for `E1`–`E7` and therefore wrongly deferred?

## Physics filter

Report a finding only if it catches a way the **physics could be wrong**, a way a **claim would be
unsupported**, or a way this list would cause **wrong or corrupted artifacts** to be built. ⛔ Not style,
⛔ not naming taste.

## Method

- ⭐⭐ **Quote both sides of every finding**: the list's text, and the source it fails against.
- ⭐ For any claim about a file, a line, or a count: **state what you opened or ran, and what you found.**
  ⛔ A path that resolves is not a source.
- ⭐⭐ `U2`, `U4`, and checks 1–3 are settled by **running code**. Write the script, run it, and paste the
  **literal stdout**. ⛔ Prose derivations are discarded for those — a hand-argued count is exactly the
  defect class this project exists to remove, relocated into the review.
- ⛔ Do **not** edit anything, and ⛔ do not execute any part of the list. Read-only.
- ⭐ End with: which of your step-6 items the list handles and which it misses; **whether S10 should be
  re-keyed before S11 is built, and the measurement that decides it**; and what you checked that could
  have failed and did not.
