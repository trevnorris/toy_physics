# Moving-Throat PDE — Physics Red-Team Status

## Current Status

The effective-physics red-team pass has now covered the main referee-risk slices:

- [PHYSICS_FOUNDATIONS_AUDIT.md](/var/projects/toy_physics/notes/moving_throat/review/PHYSICS_FOUNDATIONS_AUDIT.md)
- [PHYSICS_MOUTH_AUDIT.md](/var/projects/toy_physics/notes/moving_throat/review/PHYSICS_MOUTH_AUDIT.md)
- [PHYSICS_COEVOLUTION_AUDIT.md](/var/projects/toy_physics/notes/moving_throat/review/PHYSICS_COEVOLUTION_AUDIT.md)
- [PHYSICS_REALIZATION_INVARIANT_AUDIT.md](/var/projects/toy_physics/notes/moving_throat/review/PHYSICS_REALIZATION_INVARIANT_AUDIT.md)

## Bottom Line

On the standard you set:

- the derivation chain makes physical sense as an effective model,
- I have **not** found a place where the notes are obviously misapplying physical principles,
- I have **not** found evidence of hidden tuning or branch-by-branch rescue,
- and I have **not** found a critical physical contradiction that would make me say the derivation is broken.

That does **not** mean every closure is uniquely forced by the universe.
It means the chain currently survives an aggressive effective-model red-team pass.

## What Passed

### Foundations `001-004`

- The wall/interface, BdG reduction, and mixed outgoing bridge look physically legitimate as effective reduced sectors.
- The real limitations are scope limitations:
  effective closure, stable-sector assumptions, reduced outgoing-channel representation.

### Mouth law and susceptibility `112`, `120-124`

- The mouth boundary layer is standard positive drift-diffusion / Onsager physics.
- The self-matched susceptibility is a same-layer closure that removes freedom rather than adding a fit knob.
- No hidden tuning showed up there.

### Self-consistent and co-evolving mouth/core branch `125-127`, `137`, `139-140`

- This tranche gives some of the strongest anti-tuning evidence in the project.
- The old canonical point does **not** survive automatically once backreaction is turned on.
- Exact compensation is recovered only after a substantial renormalization cost, independently reproduced numerically.

### Realization / invariant bridge `108`, `153`, `168-170`

- The positivity selection theorem is physically clean inside its one-lane positive-source closure.
- The grouped outlet map is a sensible linear-response compression on the compensated isotropic branch.
- The monomial/orbit/quotient story is a legitimate reduced invariant classification, not fake physics, provided it is not oversold as a full PDE dynamical theorem.

## What Did Not Show Up

I did **not** find:

- a hidden parameter retuned later to rescue the branch,
- a sign convention inserted by hand to protect the endpoint,
- a place where a required value like `n=5` suddenly behaves like a context-specific fit knob,
- or a stage where the model only works because a clearly unphysical principle is being used incorrectly.

## What Still Needs Care

The remaining risk is now sharply localized.

It is **not**:

- more symbolic algebra,
- more SymPy vs Mathematica comparison,
- or more anti-tuning cleanup.

It **is**:

- the final realization question of whether the actual moving-throat PDE branch
  stays close enough to the reduced coherent closures and invariant structure for
  referee purposes.

In practice, that means the live objection is now:

- `these reduced closures may be coherent, but do you know the full PDE branch actually realizes them?`

That is the correct remaining red-team pressure point.

## Practical Verdict

If a referee attacks this successfully now, I do **not** think the winning attack is:

- “this is tuned,”
- “this is physically nonsensical,”
- or “you misused the basic transport / response / symmetry ideas.”

I think the strongest remaining attack is:

- “you have a strong reduced closure and invariant story, but the final full-PDE realization theorem is still not complete.”

That is a narrower and much more manageable position than where the project started.
