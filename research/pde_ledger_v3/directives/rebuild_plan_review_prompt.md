# Independent review — a diagnosis, a new standing rule, and a stop-work plan

## Artifact

Commit `c48c64c0` in `/var/projects/toy_physics` (branch `ledger-v2-rebuild`). Read it with
`git show c48c64c0`. The central document is `research/pde_ledger_v3/REBUILD_HANDOFF.md`.

Also changed: `.claude/skills/{build,review-legs,step-run}/SKILL.md`, `STATUS.md`,
`research/pde_ledger_v3/LAUNCH_PROMPT.md`, and a new gate
`research/pde_ledger_v3/reduction/derived_or_declared.py`.

## The claim being made

Engine scripts `emit("TAG", "prose stating a conclusion")` where no computation produced that conclusion.
Measured across every SymPy engine in v3, only ~10–20% of tags depend on any computation. A new standing
rule follows (print objects, never state conclusions; print the residual, never only assert it), and **all
new ledger physics is stopped** until every engine is rebuilt.

## What to check — in priority order

**1. ⭐⭐ IS THE STOP-WORK CALL PROPORTIONATE? Attack this hardest, and it is genuinely contestable.**
The evidence so far is that the *physics* keeps surviving scrutiny (`K₀`, the passive region, the transverse
parity argument all check out by hand) while the *provenance claim* is false. If this is a *provenance*
crisis rather than a *physics* crisis, a full engine rebuild may be large over-investment compared with
correcting the claims and rebuilding only what is load-bearing. ⭐ Argue the strongest case **against** the
plan as written, then state which you actually believe and why. ⛔ Do not simply endorse it.

**2. ⭐⭐ IS THE MEASUREMENT SOUND, or is the gate misleading?** Run it yourself:
`cd research/pde_ledger_v3 && python3 reduction/derived_or_declared.py scripts/S11bB_interface_assembly_sympy_audit.py`
Its stated limits are in the handoff. Determine independently whether the headline "10 of 133" is defensible
or an artefact of a weak perturbation set. ⛔ The handoff's own numbers are what you are testing; do not
take them as given.

**3. ⭐ IS CLAUSE 2 CORRECT?** *"Print the residual; do not only assert it — an assert IS the builder writing
down the expected output."* Find the strongest counter-case: a situation where asserting is right and
printing is insufficient or harmful. Is `compute → emit → assert` actually adequate, or does it leave a
hole?

**4. ⭐ DOES THE INPUT/OUTPUT SPLIT LEAK?** The plan supplies every equation and withholds only "an
acceptance criterion that references an expected value." Is that a real boundary, or does supplying a
sufficiently complete specification already determine the answer such that the builder can pattern-match it
without computing? If the boundary is not real, say so plainly.

**5. ⛔ INTERNAL CONTRADICTIONS INTRODUCED BY THIS COMMIT.** New guidance was added to the skills **without
reconciling the old text around it**. Specifically check whether `step-run/SKILL.md` steps 2, 3 and 7
(pre-registration move-out, do-not-read list construction, quarantine of the whole answer-bearing set) now
contradict the new "do not blind the inputs; this runs in-repo, no quarantine" guidance in the same files
and in `build/SKILL.md`. Report every contradiction with both quotations. This is expected to find
something.

**6. Does the plan repeat traps it names?** It warns against building blindness apparatus, auditing 256
tags, and checks-on-checks. Does the plan itself do any of these?

**7. Dependency order.** `S9 → S10 → S11 → S11b-A → S11b-B` is asserted. Verify from the step records that
this is the actual dependency order and that nothing else feeds in.

## Required method

Derive and verify independently. Run the gate. Read the engine sources. ⛔ Do not accept the handoff's
framing as the frame of your review — its central claims are what is under test.

## Physics filter

Report a finding only if it bears on whether the diagnosis, the rule, or the plan is **wrong**. ⛔ Do not
report style, wording, or "a different reader might misread this."

## Do not read

`research/pde_ledger_v3/_scratch/`, any other reviewer's output, anything named `PREREGISTERED`/`PREREG`,
`research/pde_audit/`.

## Output — under 35 lines

1. ⭐ One-line verdict on the plan: **PROCEED**, **PROCEED WITH CHANGES**, or **WRONG CALL**.
2. Your answer to check 1, with the strongest counter-argument stated before your conclusion.
3. Per remaining check: pass/fail with specific evidence, file and line.
4. ⭐ Every internal contradiction found, quoted from both sides.

⭐ A clean result is acceptable; ⛔ do not manufacture a finding.
