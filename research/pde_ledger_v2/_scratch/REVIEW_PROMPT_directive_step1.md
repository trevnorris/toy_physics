# REVIEW REQUEST — a build directive, before it is executed. Read-only, one pass.

> ⛔⛔ **HISTORY — NOT INSTRUCTION.** This is a step-① review prompt, completed 2026-07-31 (`407eed94`). Its "OPEN
> pending a user decision" language about `R2.a_pin` describes a question that **no longer exists** —
> the relation was retired by removal. ⛔ Do not execute anything here.
> ⇒ **Live workstream: `research/pde_ledger_v3/` (start at `NEXT_SESSION.md`).**

You are reviewing a **directive** that is about to be handed to a builder agent. It has not been
executed. Your job is to find what is wrong with it **before** it runs.

Read, in this order:

1. `/var/projects/toy_physics/research/pde_ledger_v2/_scratch/DIRECTIVE_step1_retire_apin.md`
2. Then open every file the directive touches, and the corpus loci it cites.

Repository root `/var/projects/toy_physics`, branch `ledger-v2-rebuild`.

⛔ **Do not read** `_scratch/PLAN_apin_repair_and_throat_restart.md` or anything under
`_scratch/reviews/`. They contain a prior author's predicted numbers, and one of the things being asked
of you below is an **independent** derivation of those numbers.

## Operating constraints

- **READ ONLY.** Do not modify, create, or delete anything in the repository. Do not commit.
- **One pass. Do NOT ask clarifying questions** — you are non-interactive and a question will hang.
  State any ambiguity in your review and proceed on your best reading.
- You may run Python/SymPy to check things. Write scratch files outside the repo.

## Background you need

The project is a 4D compressible-superfluid analog model (a mathematical bridge between EM and gravity
formalisms — an analog, not an ontology claim). A "throat" is a defect that punctures a domain-wall
brane; its mouth radius is called `a`.

A stage note imposes four natural-unit pins `{a, c_s0, ħ, m_GNLS}` on three base dimensions and reads
the forced leftover dimensionless monomial as a relation `a = ħ/(m_GNLS c_s0)`. The directive retires
that relation on the grounds that it is a units artifact whose *form* is forced by dimensional analysis,
and whose application to a physical throat radius is a false claim in kind — `ħ/(m c_s0)` is built from
medium primitives alone and is one number for the whole medium, while throat radii form a
species-indexed family.

The v2 ledger is becoming reference material. The directive therefore repairs only what a future ledger
will **compute with** or **read as foundation** (seven files), and deliberately quarantines the rest.

## What to produce

```
# Directive review — <your name/model>

## VERDICT
<SAFE TO RUN / RUN WITH CHANGES / DO NOT RUN> + 2-3 sentences

## BLOCKING
<numbered; each: the instruction, the file:line evidence, why it is wrong, what to change.
 Empty is a legitimate answer.>

## Q-A  Is the scope right?
Is the seven-file set correct? Does anything OUTSIDE it still *compute* with the pin relation
(as opposed to merely mentioning it) in a way that would be broken or left inconsistent by these
edits? Grep for it; do not reason from the directive's own list. Conversely, is anything IN scope
unnecessary?

## Q-B  Remove the quantity, or reclassify it?
The directive removes `Q.medium.a_pin` from the registry entirely, arguing a throat radius is a
defect-sector quantity that does not belong in a medium block. The alternative is to keep it as a
free, undetermined quantity with no defining equation. Which is right, and what does each choice
imply for what the registry reports?

## Q-C  Independently derive the post-edit numbers
`reduction/acceptance_check.py` compares a computed payload against a literal expected dict for four
cases (`baseline`, `C-M1`, `C-M2`, `C-M3`). Work out, from the edited registry as the directive
specifies it, what `dim_before`, `dim_after` and `Delta` **should** become for each case, and show your
reasoning: how many quantities are in the ambient set, how many constraints are admitted, and what the
rank is. ⭐ Derive this yourself from the registry files. Do not read it off any existing expected value.
Note explicitly if a case can no longer be constructed.

## Q-D  The able-to-fail probe
The directive says one probe's premise is destroyed by the edit. Is that correct? Is the instruction
given for handling it adequate, or could it be satisfied in a way that silently weakens the harness?

## Q-E  Locus integrity
Two edited files are cited by line number from elsewhere. Are the directive's line-preservation
instructions sufficient? Enumerate every line-pinned citation into the two files that would be at risk.

## Q-F  Does the directive leak its own answer?
Does it anywhere tell the builder what result to produce, such that a builder could satisfy it by
writing down an expected value rather than deriving one? Quote anything that does.

## Q-G  What will go wrong when this runs?
Concrete failure predictions.

## MATH FLAGS
<table: claim | file:line | your independent result | agree/disagree>
```

## Standards to apply

1. **A matching number is not evidence.** A check agreeing with its reference proves nothing until you
   know why. Two errors can cancel — that has been measured in this project three times.
2. **A check whose expected value or trigger lives inside the artifact it checks** is the artifact
   agreeing with itself.
3. **A check whose only honest outcome is an invented value is worse than no check.**
4. **Apparatus growing above the physics has killed two previous efforts here.** Flag machinery built as
   a precondition for doing physics.
5. **Falsification is welcome.** If the directive's premise is wrong, say so plainly.
6. **Dimensional agreement is not physical agreement.**

Be specific, adversarial, and cite loci. Do not be agreeable.
