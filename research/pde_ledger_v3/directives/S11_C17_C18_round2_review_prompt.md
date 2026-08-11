# Independent review — the orchestrator-authored `C17`/`C18` revision

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` (working tree, 1194
lines). ⭐ Diff it against its base:
```
git -C /var/projects/toy_physics diff 099be693 -- research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
```

⚠ **Shared spec both engines read.** An error makes both engines agree on the same wrong thing, or
**disagree for a reason that is not physics**.

⚠⚠ **A previous author's revision of this same file was REVERTED.** It closed eight items and introduced
six new defects into the material it had just written. ⭐ This revision is **orchestrator-authored**, does
**four** items instead of eight, and is the third attempt. ⇒ ⛔ **Assume nothing about it.**

You are one of two independent legs; the other is not visible to you.

## ⛔ DO NOT READ

- `directives/S11_C17_C18_fix_round1_directive.md` and `…_fix_round1_review_prompt.md`.
- Any file under `_scratch/`.

## ⭐⭐ READING ORDER — ⛔ artifact last

1. `DEFECT_REGISTER.md`, `C17` and `C18` in full.
2. `directives/S11_C17_C18_spec_repair_decisions_v2.md` — the decisions `T1`–`T7`.
3. The committed outputs `mathematica/out/S11_stray_longitudinal_mathematica_audit.out` and
   `scripts/out/S11_stray_longitudinal_sympy_audit.out`. ⭐ **Open them.** ⛔ A register line number is not
   a source; it has been stale twice.
4. **Write down** how *you* would make two blind engines emit comparable stratum information. Keep it.
5. **Only now** the artifact and the diff.

## ⭐⭐⭐ WHAT THIS REVISION CLAIMS TO DO — ⛔ judge each

| | claim |
|---|---|
| **A** | `§Q8c` is a two-pass protocol whose witness input is a **point plus a locus selector**, supplied by the **orchestrator**, that fires off the **source locus** so an empty `STRATUM_ORDERING` cannot suppress it, and whose `WITNESS<w>_DONOR_SCOPE` makes the receiver's row and the donor's row **one comparison** |
| **B** | `§4`'s quoted block is **untouched** and the permission lives outside it |
| **C** | `STRATUM<s>_<COUNT>` has **one payload form in all three statuses** |
| **D** | the component's **elimination is pinned**, so both engines describe one component in the same variables |
| **E** | `§Q10`'s two stale object references are repointed |
| **F** | corollary 4 admits the witness axis **explicitly**, and only because `WITNESS_ORDERING` is unconditional |

## ⭐⭐⭐ THE SIX DEFECTS THE PREVIOUS ATTEMPT BRED — ⛔ CHECK THIS ONE DID NOT REPEAT THEM

⭐ All six were **measured**, not argued. ⛔ Report any that recur here:

1. It **deleted** the engine-side rule *"never serialised as a host-language native boolean"* while `§6`
   still mandates a bare boolean per expression ⇒ a working class of rows became rejectable.
2. It left a cross-reference to a `§8` "prohibition" that bound nobody.
3. It made `_CANONICAL_LOCUS` fire by clearing premise-positive denominators — and on the step's central
   locus the canonical generator then **divides out the very factor** cutting the critical component,
   measured identically in both CASes.
4. It left the emitted basis **list order** unpinned; the two CASes were measured to disagree.
5. Its new component-wide tokens were spelled outside the only rule telling a comparator how to treat
   *undecided*.
6. The witness pass was a fourth emission axis corollary 4 did not permit.

⭐ **Diff against `099be693` and confirm 1–4 are absent** (the base has the engine-side boolean rule intact
and `_CANONICAL_LOCUS` inert). ⛔ If this revision reintroduced any, say so with the line.

## What to check

### 1. ⛔⛔ Does the file state or imply an answer?
⭐ Clause by clause over the added lines. **The test is DERIVABILITY, ⛔ not literal presence.** ⛔ Could a
builder reading only this file work out a count, value, sign, or named locus?

### 2. ⭐⭐ Does `Q8c` actually EXECUTE, and does it produce a cross-engine object?
⛔ **Trace it on `XFORM_EXTRA, D = 2` from the committed payloads.** ⭐ Does an input get supplied under the
source-locus trigger? What does the receiver emit? ⛔ **Can a compliant pair still reach silence?**
⭐ Then: does `DONOR_SCOPE` actually let a comparator pair the receiver's row with the donor's, or is it
still two per-engine computations? ⚠ Name the two tags that would be joined.

### 3. ⭐⭐ Is the witness input still exactly one record, and is the blind build intact?
⚠ It now carries a **selector** as well as a point. ⛔ Is the selector a route to receiving a value?
⛔ Could a builder read it as licence to receive anything else? ⭐ Is the Wolfram engine still importing
nothing of the other engine's **computation**?

### 4. ⭐⭐ Does `D`'s pinned elimination work in both CASes?
⭐ **Run it.** Take a real component from the committed outputs and apply the rule in each CAS. ⛔ Do both
engines retain the same variables? ⚠ What happens where the equation determines no variable, or several?

### 5. ⭐ Does `C` give `<COUNT>` one type — and is the `VARIES` pair actually the answer?
⛔ Under `VARIES`, is *"the count off the change sub-locus"* well defined and computable?

### 6. ⚠ Stale prose above the fix
⭐ For every edited section read its introduction. ⛔ Is any sentence now false? ⚠ **This file has produced
one such instance per round, in all three rounds.**

### 7. ⭐ What of `T1`–`T7` is still not discharged?

## Physics filter
Report a finding only if it catches a way the **physics could be wrong**, a **claim unsupported**, or
**wrong or incomparable artifacts**. ⛔ Not style.

## Method
- ⭐⭐ Quote both sides of every finding.
- ⭐ State what you opened or ran for every claim about a file, line or count.
- ⭐⭐ Checks 2 and 4 are settled by **running it** on a real case. ⛔ Prose is discarded there.
- ⭐ Any CAS run: `timeout 600`, ⛔ one kernel at a time (two seats), ⛔ copy to `/tmp` and run the copy.
- ⛔ Read-only. Do not edit the working tree.
- ⭐ End with: A–F held or not; whether any of the six bred defects recurs; whether the file leaks; what
  `T1`–`T7` still lacks; and what you checked that could have failed and did not.
