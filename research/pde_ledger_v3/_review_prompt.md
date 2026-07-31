# REVIEW REQUEST — the v3 ledger plan set, before execution. Read-only, one pass.

You are reviewing the plan for a new ledger ("v3") in a toy-physics project, **before any of it is
executed**. Your job is to find what is wrong with it now, while it is cheap.

Repository root `/var/projects/toy_physics`, branch `ledger-v2-rebuild`, HEAD `2db79549`.

**Read these four, in order:**

1. `research/pde_ledger_v3/CHARTER.md` — scope and what it excludes
2. `research/pde_ledger_v3/DEFECT_REGISTER.md` — the known error surface
3. `research/pde_ledger_v3/V3_STEP_PLAN.md` — ⭐ **the main object of this review**
4. `research/pde_ledger_v3/SESSION_REASONING.md` — how the conclusions were reached

Method doc (referenced, not restated by the plan): `docs/derivation_walkthrough_plan.md`.

**Then open the loci.** Every claim in these documents is checkable against the corpus, and the author
has been wrong repeatedly this session in ways a reader who opened the file would have caught
(`SESSION_REASONING.md` §8 lists two; §4 lists a third).

## Operating constraints

- **READ ONLY.** Do not modify, create, or delete any file. Do not commit.
- **One pass. Do NOT ask clarifying questions** — you are non-interactive and a question will hang.
  State ambiguities in your review and proceed on your best reading.
- You may run Python/SymPy. Write scratch files outside the repository.
- Prefer opening the cited file over trusting the document's summary of it.

## Background

One 4D compressible superfluid (GNLS, `ρ=|ψ|²`, `P=Kρ⁵`); our 3D space is a domain-wall brane in it; a
particle is a "throat" puncturing that brane. Gravity is the drain flow between throats, light is brane
shear. The project is an **analog** — a mathematical bridge between EM and gravity formalisms — not an
ontology claim, and it is judged by calibrate-then-predict surplus on **dimensionless** held-out ratios.

v3 is scoped to **gravity + light + gravitomagnetism**, deliberately excluding the throat interior. Two
prior efforts in this project died from *apparatus growing above the physics*; that is the failure mode
most worth catching.

## What to produce

```
# v3 plan review — <your name/model>

## VERDICT
<SOUND / SOUND WITH CHANGES / DO NOT PROCEED> + 3-4 sentences

## BLOCKING
<numbered. Each: the claim or instruction, file:line evidence, why it is wrong, what to change.
 Empty is legitimate if you find nothing blocking.>

## Q1  Is the step decomposition right?
Are the 23 steps in V3_STEP_PLAN the right cuts? Is any step doing two things that should be
separated, or two steps doing one thing? Is anything in the gravity/light/gravitomagnetism scope
MISSING entirely? Check against the corpus, not against the plan's own narrative.

## Q2  Is the dependency order correct?
Does any step depend on something a later step introduces? Name the specific pair.

## Q3  Is the scope boundary honest?
The plan claims gravity is independent of the throat interior, resting on a worldtube reduction.
Verify that claim at its loci. Then: does any step in phases 2-4 actually need interior data that
the plan does not admit it needs?

## Q4  The defect register
Is any row wrong, mis-classified, or missing a locus? ⭐ More importantly: what is MISSING from it?
The register claims ten pin-shaped identifications (two same-dimension quantities silently equated) —
find one it does not list.

## Q5  Apparatus creep
Where does this plan rebuild census-grade machinery instead of doing physics? Be specific about which
step or instruction.

## Q6  The physics reasoning
SESSION_REASONING §5 claims there is NO surviving mass-radius slope, because the only route giving a
definite sign is the same construction that produces the falsified 1:9:25 tower. Verify or refute
that at `notes/lepton_mass_notes.md`. It is the session's headline conclusion.

## Q7  What will actually go wrong in the first three steps?
Concrete predictions.

## MATH FLAGS
<table: claim | file:line | your independent result | agree/disagree>

## WHAT IS SOLID
<what you checked and found correct — this matters as much as the defects>
```

## Standards to apply

1. **A matching number is not evidence** — a check agreeing with its reference proves nothing until you
   know why. Two errors can cancel; measured three times in this project.
2. **Dimensional agreement is not physical agreement.** Two different lengths agree on being lengths
   trivially. Ask whether both sides are indexed by the same thing — one number for the whole medium, or
   one per particle.
3. **Falsification is welcome and first-class.** If the plan's premise is wrong, or the model is broken,
   say so plainly. Never soften.
4. **Apparatus above physics has killed two efforts here.** Machinery built as a *precondition* for
   doing physics is the specific failure.
5. **A check whose expected value or trigger lives inside the artifact it checks** is the artifact
   agreeing with itself.
6. **Absence of a denial is not evidence** — that a route is unexecuted says nothing about whether it is
   buildable.

Be specific, adversarial, and cite loci. Do not be agreeable.
