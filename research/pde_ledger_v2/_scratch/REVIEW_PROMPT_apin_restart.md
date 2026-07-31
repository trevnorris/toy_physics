# REVIEW REQUEST — read-only, one pass

You are reviewing a research plan for a toy-physics project: a 4D compressible-superfluid analog model
intended as a mathematical bridge between the EM and gravity formalisms. It is an **analog**, not an
ontology claim.

**Read this file first, in full:**

    /var/projects/toy_physics/research/pde_ledger_v2/_scratch/PLAN_apin_repair_and_throat_restart.md

The repository root is `/var/projects/toy_physics`, branch `ledger-v2-rebuild`. The plan cites specific
files and line numbers throughout. **Open them.** Every claim in the plan is checkable against the
corpus, and the plan's author has already been wrong three times this session in ways that a reader who
opened the file would have caught (the plan's §9 lists them).

## Operating constraints

- **READ ONLY.** Do not modify, create, or delete any file in the repository. Do not commit.
- **One pass. Do NOT ask clarifying questions** — you are running non-interactively and a question will
  hang. If something is ambiguous, state the ambiguity in your review and proceed on your best reading.
- You may run SymPy / Python to check algebra. Do not write into the repo; use a temp location.
- Prefer opening the cited file over reasoning from the plan's summary of it.

## What to produce

Answer the plan's **§8 questions Q1–Q7** directly, then give an overall verdict. Structure:

```
# Review — <your name/model>

## VERDICT
<one of: PROCEED / PROCEED WITH CHANGES / DO NOT PROCEED>  + 2-3 sentences

## BLOCKING
<numbered. Each: the claim, the file:line evidence, why it is wrong or unsupported, what to do.
 Empty is a legitimate answer if you find nothing blocking.>

## Q1 unit-system rank claim
## Q2 the gravity-sector consumers of `a`
## Q3 the a ∝ m^(1/3) scaling
## Q4 the THROAT_DRAIN_DESTABILIZED result
## Q5 restart ordering
## Q6 apparatus creep
## Q7 what was missed
<each: a direct answer with file:line loci; say "cannot determine from the corpus" where true>

## MATH FLAGS
<table: claim | file:line | your independent result | agree/disagree>

## WHAT IS SOLID
<what you checked and found correct — this matters as much as the defects>
```

## Standards this project holds itself to — apply them

1. **A matching number is not evidence.** A check that agrees with its reference proves nothing until you
   know *why*. Two errors can cancel; this has been measured in this project three times. The delta
   column is usually the tell.
2. **Falsification is welcome.** A result that breaks the model is a first-class outcome, never something
   to rescue or soften. If the physics is wrong, say so.
3. **A check whose expected value or trigger lives inside the artifact it checks** is the artifact
   agreeing with itself.
4. **Apparatus growing above the physics has killed two previous efforts here.** Machinery built as a
   *precondition* for doing physics is the specific failure mode. Flag it.
5. **Dimensional agreement is not physical agreement.** Two different lengths agree on being lengths
   trivially. The plan exists because that error was made at the foundation of the model.
6. **Absence of a denial is not evidence.** A record that a route is *unexecuted* says nothing about
   whether it is *buildable*.

## The one thing most worth your effort

The plan's §6 argues that the model's own reduced energy ledger implies a **heavier lepton has a bigger
throat** (`a ∝ m^{1/3}`), contradicting an external estimate in the corpus that implies the opposite
(`a ∝ 1/m`). That derivation was done quickly and is load-bearing for the whole restart. Check it
properly — keep the sub-leading term the plan drops, and decide whether the two reduced models the corpus
contains are actually the same model.

Be specific, be adversarial, and cite loci. Do not be agreeable.
