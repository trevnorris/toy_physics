# Moving-Throat PDE — Physics Mouth Audit (`112`, `120-124`)

## Scope

This pass targets the first explicit mouth-layer and mouth-gain closure chain:

- [stage112_mouth_boundary_layer.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage112_mouth_boundary_layer.md)
- [stage120_core_to_mouth_gain_map.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage120_core_to_mouth_gain_map.md)
- [stage121_normalized_mouth_gain_family.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage121_normalized_mouth_gain_family.md)
- [stage122_family1_actual_mouth_gains.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage122_family1_actual_mouth_gains.md)
- [stage123_selfmatched_mouth_susceptibility.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage123_selfmatched_mouth_susceptibility.md)
- [stage124_mouth_gain_status.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage124_mouth_gain_status.md)

with physical provenance from:

- [stage043_entropic_microclosure.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage043_entropic_microclosure.md)
- [stage045_parent_action_gain.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage045_parent_action_gain.md)
- [stage046_parent_thresholds.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage046_parent_thresholds.md)
- [stage047_equilibrium_alignment.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage047_equilibrium_alignment.md)

The standard here is not ontology. The standard is:

- does the mouth-layer derivation make physical sense as an effective model,
- are the transport and susceptibility choices coherent,
- is anything obviously being tuned by hand,
- and is there any point where a referee could fairly say the physics is being misapplied.

---

## Bottom Line

Current verdict:

- the `112`, `120-124` mouth chain looks **physically coherent as an effective
  mouth-layer closure**;
- I do **not** see evidence here of hidden tuning or branch-by-branch rescue;
- I do **not** see a point where a physicist should be able to dismiss the chain
  as “completely unphysical.”

The real referee-facing pressure points are narrower:

1. [stage112_mouth_boundary_layer.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage112_mouth_boundary_layer.md) uses a legitimate zero-flux Onsager mouth law, but the linear potential `V_m(z)≈V_1 z` is still a reduced mouth-layer approximation on the active interval, not a global theorem for the whole throat potential.
2. [stage123_selfmatched_mouth_susceptibility.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage123_selfmatched_mouth_susceptibility.md) is a **matched-layer closure**: it identifies the source susceptibility with the same active shell layer used for the shell/compliance mode. That is plausible and structurally motivated, but it is not yet the unique microscopic possibility.
3. [stage124_mouth_gain_status.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage124_mouth_gain_status.md) should always be read as a status statement **inside the explicit throat-core plus mouth-layer closure**, not as a full-PDE uniqueness theorem.

Those are important limitations, but they are not signs of hidden tuning or bad transport physics.

---

## Findings

## 1. No hidden tuning detected in the mouth-law chain

This is the main positive result of the pass.

Why:

- [stage112_mouth_boundary_layer.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage112_mouth_boundary_layer.md) collapses the source-profile freedom to the dimensionless bias `Pi_m = V_1 L / Theta_sigma`; it does not fit `Pi_m` to later output targets.
- [stage120_core_to_mouth_gain_map.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage120_core_to_mouth_gain_map.md) and [stage121_normalized_mouth_gain_family.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage121_normalized_mouth_gain_family.md) derive `(M_s,M_q)` and `R_q` from the already frozen core coefficients rather than inserting a gain pair by hand.
- [stage123_selfmatched_mouth_susceptibility.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage123_selfmatched_mouth_susceptibility.md) removes the last free susceptibility scale by identifying it with the already active shell layer. That move decreases freedom; it does not create a new fit knob.

The strongest anti-tuning evidence is downstream:

- when the later co-evolving solve is done in [stage139_renormalized_canonical_branch.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage139_renormalized_canonical_branch.md), the model is forced to pay a substantial traction increase to restore exact compensation.

If the mouth closure were secretly tuned to preserve the preferred branch cheaply, that is not what the later chain would look like. The later renormalization cost is good evidence that the mouth closure is constraining the model rather than rescuing it.

---

## 2. Stage `112` is physically standard drift-diffusion, not a misuse of transport theory

[stage112_mouth_boundary_layer.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage112_mouth_boundary_layer.md)
is physically conservative in the right way.

What is good:

- the free energy is the standard positive-density entropic form plus an external effective potential;
- the chemical potential is the corresponding variational derivative;
- the Onsager current `J_sigma = -M_sigma sigma partial_z mu_sigma` is the standard positivity-preserving drift-diffusion law;
- the stationary zero-flux branch gives the exact truncated exponential family.

That is all physically normal. A referee may dislike the modeling choice, but it is not a misuse of electrochemical or Onsager reasoning.

The real limitation is narrower:

- the stage linearizes the mouth potential to `V_1 z` and then uses that law over the active interval `[0,L]`.

That is acceptable as a reduced mouth-layer closure. It should not be sold as the exact global throat potential away from the mouth.

Also:

- normalizing `int_0^L sigma dz = 1` is physically fine here because the stage is clearly separating profile shape from total source amplitude. It is a reduced shape law, not a claim that the total source inventory is literally one particle.

---

## 3. The self-matched susceptibility closure is plausible and restrictive, not arbitrary

[stage123_selfmatched_mouth_susceptibility.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage123_selfmatched_mouth_susceptibility.md)
is the place where a referee is most likely to ask whether a hidden knob has been smuggled in.

My current answer is no.

Why:

- the closure `Theta_sigma = H_w J_s` is not chosen to hit a target value of `M_s`;
- it is inherited from the same-layer logic already visible in [stage047_equilibrium_alignment.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage047_equilibrium_alignment.md), where matched support/source layers drive the branch toward `C_(sigma phi)^2 = 1`;
- the source susceptibility is being tied to the already active shell layer instead of being left free.

So the closure is best described as:

- a **matched-layer approximation**,
- not a fitted parameterization,
- and not an arbitrary change of numerical constants to rescue the canonical branch.

The legitimate objection is different:

- if the actual mouth source lives on a materially different layer than the shell/compliance mode, the numerical prefactor `20/9` need not survive unchanged.

That would be a model-realization objection, not evidence of hidden tuning inside the present chain.

---

## 4. The derived gain formulas are structurally coherent

[stage120_core_to_mouth_gain_map.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage120_core_to_mouth_gain_map.md)
and [stage121_normalized_mouth_gain_family.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage121_normalized_mouth_gain_family.md)
use the same electrochemical combination that motivated Stage `112`:

- shell contribution enters with one sign,
- localized-Maxwell mixed contribution enters with the opposite sign,
- and the ratio law `R_q = (g_c-r)^2/(1+r^2)` is inherited from the exact core Schur complement.

That is physically coherent. The mouth layer is not inventing a new sign convention to get the desired branch.

This is one of the stronger anti-red-flag points in the mouth chain: once the core outlet is accepted, the mouth gains inherit its sign structure rather than redefining it locally.

---

## 5. The real unresolved issue is branch realization, not principle misuse

The strongest mouth-side objection is not “this uses bad physics.” It is:

- does the actual moving-throat mouth layer really stay on the zero-flux positive-density branch,
- and does the source susceptibility really track the same active shell layer used in the shell/compliance mode.

Those are legitimate realization questions.

But they are not the same as:

- using the wrong transport law,
- violating positivity,
- inserting branch-specific rescue parameters,
- or tuning coefficients after the fact.

So the right referee wording is:

- the mouth chain is a physically coherent explicit closure,
- but some same-layer and branch-realization assumptions remain closure assumptions rather than full-PDE theorems.

---

## Non-Findings

These are things I explicitly do **not** think are mouth-chain problems.

1. The truncated exponential source family is not an arbitrary profile guess anymore.
   Within the reduced mouth-layer model, it is the exact zero-flux equilibrium branch.

2. The use of an entropic logarithm for the positive source density is not suspicious.
   It is the standard convex free energy for positive drift-diffusion.

3. The susceptibility closure is not a hidden fitting parameter.
   It removes a free scale instead of introducing one.

4. The derived `~4%` traction shift in [stage123_selfmatched_mouth_susceptibility.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage123_selfmatched_mouth_susceptibility.md) is not evidence of tuning.
   It is just the output of the chosen closure before the later full co-evolving correction.

5. The fact that later stages demand much larger traction is not a contradiction.
   It is evidence that the chain is willing to pay a nontrivial self-consistency cost instead of protecting a preferred answer.

---

## Recommended Wording Discipline

To keep the mouth chain referee-hard:

- say `reduced mouth-layer closure` rather than implying the full throat potential is globally linear;
- say `self-matched same-layer closure` rather than implying `Theta_sigma = H_w J_s` is already unique microscopic truth;
- say `within the explicit throat-core plus mouth-layer closure` for status summaries;
- avoid language suggesting the mouth law has already been proved to be the only admissible full-PDE realization.

---

## Mouth Verdict

If I were trying to kill the project as a referee, I would **not** attack `112`,
`120-124` by saying the transport theory is nonsensical or the gains are being tuned by hand.

I would instead attack it by saying:

- the mouth law is still a reduced closure on the active interval,
- the self-matched susceptibility is a same-layer approximation rather than a unique microscopic derivation,
- and the full PDE realization of that branch is still not proved.

That is a much narrower objection than “bad physics” or “hidden tuning.”

So this slice currently passes the effective-model legitimacy test:

- **no hidden tuning found,**
- **no obvious misuse of physical principles found,**
- **but the same-layer and branch-realization assumptions must stay explicitly labeled as closure assumptions.**

---

## Next Slice

The next physics red-team target should be the co-evolving mouth/core branch:

- [stage125_selfconsistent_mouth_branch.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage125_selfconsistent_mouth_branch.md)
- [stage137_coevolving_core_mouth_map.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage137_coevolving_core_mouth_map.md)
- [stage139_renormalized_canonical_branch.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage139_renormalized_canonical_branch.md)

That is where a referee is most likely to ask whether the model is now paying a genuine self-consistency cost or still hiding a branch-selection bias.
