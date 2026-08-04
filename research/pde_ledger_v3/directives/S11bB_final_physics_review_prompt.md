# FINAL PHYSICS-ONLY REVIEW — S11b-B build directives, before the build runs

## Artifact

- `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11bB_wl_directive.md`
- `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11bB_py_directive.md`

Engine header + a shared physics block, byte-identical in both
(`sha256 227bbffa1c638b8b07ed24327e7a7f3dfb8f86163f27a066167cbf3c962eacc1`).

## ⛔⛔ THE ONLY QUESTION

> **Does this specification contain a defect that would make BOTH engines wrong about the PHYSICS?**

That is the whole review. ⚠ This directive has been through eleven revisions and many review rounds; the
remaining work is deliberately being closed out. ⇒ ⛔ **DO NOT report:**
- check coverage, check strength, or whether a check could be stronger;
- wording, emphasis, ordering, tag names, or formatting;
- read-bar / blindness gaps — **these are being handled operationally, outside the directive**;
- anything whose worst consequence is that the two engines **disagree** (that costs a round; it does not
  certify wrong physics);
- anything already recorded in the directive as a stated coverage limit or a known blind spot — ⭐ those are
  **deliberate and disclosed**, and a disclosed limit is not a defect.

⭐ **Report only what would land in BOTH engines and make their agreement certify something false.**

## What that looks like concretely

A finding qualifies if a competent engine following the text would derive **wrong physics**, and the other
engine following the same text would derive the **same wrong physics**. For example: a supplied equation
that is wrong; a sign or normalization error in something supplied rather than computed; a definition that
is internally inconsistent; a missing term in a relation the engines are told not to re-derive; a
prescription that is mathematically invalid.

⚠ **Weight your attention on everything the directive SUPPLIES rather than asks for**, because a supplied
object cannot be falsified by the engines' disagreement — it lands in both identically. Specifically:
the complex-frequency continuation (§1b), the affinity and the interfacial mass balance (§2b), the supplied
face-response acceptance value, the virtual-displacement rule and the balance-law route (§3b), the exact
mass balance (B1), the stored-energy list and symmetry group (§3), and the Onsager conversion procedure.

⭐ **Derive each of those independently and compare.** ⛔ **Do NOT verify any of them against the
directive's own power identity or its own internal relations — that is circular.** Use first principles.

## Method

⚠ **Computation required, not reading.** Use SymPy/numerics freely. Scratch scripts under
`/tmp/claude-1000/-var-projects-toy-physics/9640c755-fe17-4ce5-8856-7593e98346ca/scratchpad/`, prefixed
`final_`.

⚠ Where you believe a task will produce a specific answer, ⛔ do not report that answer — the engines must
compute it.

## Do not read

Anything named `PREREGISTERED`/`PREREG` · `research/pde_ledger_v3/_scratch/` (prior reviews and authoring
transcripts) · `steps/S11b_HANDOFF.md` · the other `S11bB_rev*` files · the git history of the directives ·
`research/pde_audit/`.

⭐ You may read `steps/S11bA_interface_response.md`, `V3_STEP_PLAN.md` and `reduction/`.

## Output

Under 30 lines.

1. ⭐ **A one-line verdict: BUILD or DO NOT BUILD.**
2. For each supplied object listed above: **VERIFIED / DISAGREES / NOT CHECKED**, with the disagreement
   stated if any.
3. Any qualifying finding — what is wrong, and the concrete physics consequence.

⭐ If nothing qualifies, **say so plainly and say BUILD.** ⚠ A clean result here is an expected and
acceptable outcome; ⛔ do not manufacture a finding to justify the review.
