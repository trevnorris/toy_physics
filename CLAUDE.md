# How we work

Sixteen rules. Every one exists because ignoring it cost a session.
Full process: `docs/development_pipeline.md`. What we're building: `docs/development_plan.md`.
Where we are: `STATUS.md`.

## The method — the CAS answers, not me

1. **Two engines exist so they can disagree.** Independent construction, not hidden answers. The
   disagreement is the measurement.
2. **A script prints computed objects. It never states conclusions.** Emit both operands and the residual,
   then guard — a residual asserted zero always prints `0` and carries no information. Interpretation
   belongs to the step record.
3. **Name the object. Do not specify the recipe.** If a review is arguing about *how* to derive something —
   is this quotient well-defined, is this weight unique — the question was manufactured by specifying a
   derivation path. Ask for the object; let the engine hand over what it built.
4. **If I am deciding in prose what the engines should compute, I have inverted the method.** The tell:
   many turns reasoning toward an answer a script would settle in one. The cause, measured: the instrument
   was broken, so no measurement was available. Fix the instrument; don't reason around it.
5. **The spec says what to compute — never what anything equals**, is expected, or was measured. Withhold
   exactly one thing: an acceptance criterion referencing an expected value. A builder iterating to exit 0
   converges on any target it can see. Supply everything else, as equations.
6. **A disagreement is a finding.** Don't try to make divergence impossible with more careful prose; that
   defeats the reason there are two engines.

## The gates

7. **Whatever writes does not review.** Two legs on anything physics-bearing — and a spec both engines read
   is physics-bearing, because an error there makes both engines agree on the same wrong thing.
   Orchestrator-written → Codex + Grok. Codex-written → fresh Claude agent + Grok.
   **TRIGGER — no builder launches until its decision list has had two legs.** The list is
   orchestrator-written and is the one artifact the *builder* trusts: everything downstream is checked
   twice, the list is checked zero times. One pass, then fold and go — never iterated to green.
   Measured 2026-08-09: six spec defects, each costing a build round *plus* two legs, when two legs before
   the build would have caught them — three "level-above" misses, one exception named instead of the
   property (which bred a regression), measured counts stated four lines above "do not state the counts",
   and **an acceptance test that would have passed with the defect still in place.**
8. **Launch legs on sight, before I look at the result.** A self-check discharges the felt need for an
   independent one, and it is most convincing when it finds things.
9. **No commit before both legs report.** The commit is the last step. Reviewing the directive does not pay
   the tax for the build. My own verification is not a leg.
10. **Stop when nothing outstanding changes what is computed or what may be claimed** — not when both legs
    are green. "A leg that finds nothing is weak evidence" is my prior; put it in a leg's prompt and it
    becomes a quota.

## The traps

11. **Correctness is king; cost is never a reason** to drop a control, narrow a check, or skip a leg.
    Scaling work down is the user's call.
12. **Don't build blindness apparatus.** The measured failure is absence of computation, not anchoring.
    Quarantine is cut; rule 2 replaced it. A do-not-read list is a denylist, and a denylist means the
    architecture is wrong.
13. **A finding is not a mandate — verify it myself.** Legs have been wrong in both directions.
14. **Ablate to test; don't read.** A form control tests physics; a coefficient control tests arithmetic.
    Demand a script and its literal stdout — a prose re-derivation is the same defect relocated into the
    review.
15. **If successive revisions keep breeding defects in the material just changed, change the author.** Don't
    fold a fourth time.
16. **Prior art is an oracle, never a premise.** Check our computed result against it; never assume its
    result for our object. Its conditions may not be ours.
