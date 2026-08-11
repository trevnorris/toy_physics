# Independent review — fix round 1 on the repaired `S11_SHARED_PHYSICS.md`

## Artifact under review

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` (working tree).

⭐ **The version this fix was applied to** is commit `8454b72f`. See exactly what moved with:
```
git -C /var/projects/toy_physics diff 8454b72f -- research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
```

⚠ This is the **shared spec both engines read**. An error makes both engines agree on the same wrong
thing, or **disagree for a reason that is not physics**.

⚠⚠ **This is the SECOND revision of this file today.** ⭐ Round 1 closed one defect and **introduced an
internal contradiction**. ⇒ ⛔ **The question is not only "are the eight items fixed" but "did fixing them
break something that was working."**

You are one of two independent legs; the other is not visible to you.

## ⛔⛔ DO NOT READ

- `directives/S11_C17_C18_fix_round1_directive.md` — ⚠ **the directive this fix was built from.**
  ⭐ An artifact can satisfy its directive completely and still fail its source.
- Any file under `_scratch/`.

## ⭐⭐ REQUIRED READING ORDER

1. `DEFECT_REGISTER.md` entries `C17` and `C18` **in full**.
2. `directives/S11_C17_C18_spec_repair_decisions_v2.md` — the decisions `T1`–`T7` that are owed.
3. The committed engine outputs `mathematica/out/S11_stray_longitudinal_mathematica_audit.out` and
   `scripts/out/S11_stray_longitudinal_sympy_audit.out`. ⭐ **Open them.** ⛔ A line number in a register
   is not a source — it has been stale twice today.
4. **Only now** the working-tree file and the diff against `8454b72f`.

## ⭐⭐⭐ THE ITEMS THIS ROUND OWED

Round 1's two legs and the orchestrator raised eight. ⭐ Judge each **discharged or not**, quoting the
sentence that settles it:

| | what was wrong |
|---|---|
| `R1` | `Q8c`'s received point is an expression over the physical symbols the receiving engine did **not** compute ⇒ ⛔ forbidden by `§4` ("verbatim, non-negotiable") and outside `§5`'s **closed** live-read exemption list |
| `R2` | `Q8c` cannot execute on a first pass, and the text **sanctioned** an empty witness ordering ⇒ ⛔ both engines compliantly emit nothing and the measured case reverts to silence |
| `R3` | the component count tag's payload was pinned for *cannot-compute* but **open** for `VARIES` |
| `R4` | the component's **parameterisation** was unpinned ⇒ ⛔ symbolic payloads differ for a non-physics reason |
| `R5` | `_CANONICAL_LOCUS` was inert wherever a premise-positive denominator appeared |
| `R6` | *"whose payload is an integer count"* was a **value-inspected** predicate ⇒ ⛔ non-parallel tag sets |
| `R7` | `§Q10` still named a pre-repair object, now ambiguous between two |
| `R8` | *"never a native boolean"* could not discriminate the two serialisations `§8` permits |

## ⭐⭐⭐ WHAT TO ATTACK HARDEST — ⛔ `R1` IS THE DANGEROUS FIX

`R1` was fixed by making a received point a **permitted input**. ⚠ That is one sentence away from
authorising an engine to read the other engine's output.

⛔⛔ **The Wolfram engine importing NOTHING is the only blind build in this project.** It is what makes the
cross-engine residual a measurement rather than a restatement.

⭐ **Therefore test, and report with quotes:**
- ⛔ Does the new exemption authorise **anything beyond** that one point? Could a builder read it as
  licence to receive equations, a matrix, a result, a status, or "whatever the orchestrator supplies"?
- ⛔ Does it weaken `§4`'s rule, or `§5`'s **closed** exemption list, for any **other** object?
- ⛔ Is the point still **evaluated against the receiver's own** objects, with no adjustment toward the
  counterpart?
- ⚠ Does `§4`'s corollary-1 test (*"if you deleted the computation, would this tag change?"*) now have an
  exception a builder could apply **elsewhere**?

## What else to check

### 1. ⛔⛔ Does the file state, or imply, an answer?

⚠ The fix directive quoted **measured payloads** from the committed outputs, and **none of it may appear
in this file.** ⭐ **The test is DERIVABILITY, ⛔ not literal presence.** Go clause by clause over the added
lines: ⛔ could a builder reading only this file work out a count, a value, a sign, or a named locus?

### 2. ⭐⭐ Did the fix break something that was working?

⭐ Round 1 established several things two legs confirmed: `§5`'s five original locus objects intact and
strengthened; `C17` closed on a real component; four **symmetric, exhaustive** status-token sets; no leak.
⛔ **Diff against `8454b72f` and check each still holds.** ⚠ Name anything weakened.

### 3. ⭐⭐ Is `C18`'s measured case now a COMPUTATION?

⭐ Trace `XFORM_EXTRA, D = 2` through the **current** text using the committed payloads. ⛔ Under `R2`'s new
run shape, what does each engine emit, and what can a comparator conclude? ⚠ *"Not compared"* is **not** a
pass where the mode count may move. ⭐ Say whether the outcome is a **computation**, a **finding**, or
**silence** — and ⛔ whether a compliant engine can still reach silence.

### 4. ⭐ Buildable by two engines that never speak?

⭐ For every obligation, can **each** CAS discharge it independently? ⚠ `§5` records a **measured**
asymmetry between the two CASes — check the new text against it. ⛔ Does `R5`'s denominator-clearing rule
produce the **same** object in both CASes? ⭐ Test it on equations copied from the committed outputs.

### 5. ⭐ Did the fix do more than `R1`–`R8` required?

⭐ Diff against `8454b72f`. ⛔ Name anything changed that no item required, and what it damaged.

### 6. ⚠ Stale prose above the fix

⭐ For every edited section read its **introduction** and any summary above it. ⛔ Is there a sentence the
edit made **false**? ⚠ This file has produced one such instance per round so far.

### 7. ⭐ What is still missing against `T1`–`T7`?

## Physics filter

Report a finding only if it catches a way the **physics could be wrong**, a way a **claim would be
unsupported**, or a way two engines built from this file would produce **wrong or incomparable** results.
⛔ Not style, ⛔ not wording taste.

## Method

- ⭐⭐ **Quote both sides of every finding**: the file's text, and the source it fails against.
- ⭐ For any claim about a file, a line, or a count: **state what you opened or ran, and what you found.**
- ⭐⭐ Checks 3 and 4 are settled by **running the instructions on a real case from the committed outputs**
  and reporting what an engine would emit. ⛔ A prose argument that it "should work" is discarded.
- ⭐ If you run a CAS: wrap **every** kernel run in `timeout 600`; a 600 s hit is a **failed** check.
  ⛔ Never raise it, ⛔ never more than one kernel at a time (two seats). ⛔ Copy to `/tmp` and run the copy.
- ⛔ Do **not** edit anything in the working tree. Read-only.
- ⭐ End with: `R1`–`R8` discharged or not; `T1`–`T7` discharged or not; whether the file leaks; whether
  anything that worked in `8454b72f` was broken; and what you checked that could have failed and did not.
