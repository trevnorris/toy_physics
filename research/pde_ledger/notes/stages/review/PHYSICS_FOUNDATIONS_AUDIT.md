# Moving-Throat PDE — Physics Foundations Audit (`001-004`)

## Scope

This is the first effective-physics red-team pass on the foundational derivation
chain:

- [stage001_geometry_lift.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage001_geometry_lift.md)
- [stage002_breathing_reduction.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage002_breathing_reduction.md)
- [stage003_bdg_coupling.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage003_bdg_coupling.md)
- [stage021_reduced_one_port_normal_form.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage021_reduced_one_port_normal_form.md)

The standard here is not ontology. The standard is:

- does the effective derivation make physical sense,
- are principles being used coherently,
- is anything obviously unphysical or misapplied,
- and is there any sign of hidden tuning.

---

## Bottom Line

Current verdict:

- the `001-004` foundation looks **physically legitimate as an effective-model
  setup**;
- I do **not** see evidence here of hidden tuning or “make it work by hand”
  parameter adjustment;
- I do **not** see a stage where a physicist should be able to say “this is just
  nonsense physics.”

But there are still three real referee-facing pressure points:

1. the wall/interface theory in [stage001_geometry_lift.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage001_geometry_lift.md) is a **controlled effective ansatz**, not a unique derivation from the parent theory;
2. the stable BdG elimination in [stage003_bdg_coupling.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage003_bdg_coupling.md) is only as good as the assumed branch stability/spectral separation;
3. the outgoing mixed-port attachment in [stage021_reduced_one_port_normal_form.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage021_reduced_one_port_normal_form.md) is a legitimate reduced representative, but not yet the unique microscopic radiation theorem.

Those are important limitations, but they are **not** the same as hidden tuning
or physical incoherence.

---

## Findings

## 1. No hidden tuning detected in `001-004`

This is the main positive result of the first red-team pass.

Why:

- [stage001_geometry_lift.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage001_geometry_lift.md) introduces constitutive functions `mu_eta, T_w, T_Omega, K_eta`, but keeps them symbolic rather than fitting them to later target values.
- [stage002_breathing_reduction.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage002_breathing_reduction.md) reduces those symbolic coefficients to overlap matrices; again, no later-output fitting enters.
- [stage003_bdg_coupling.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage003_bdg_coupling.md) derives Schur-complement self-energies structurally from couplings and stable support frequencies.
- [stage021_reduced_one_port_normal_form.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage021_reduced_one_port_normal_form.md) derives the transfer factor and odd `l=2` fingerprint from the reduced model algebra; the sign and `omega^5` structure are not inserted by hand.

So the early chain is not succeeding because of local retuning. It is
succeeding because once the reduced ingredients are chosen, the resulting
algebra is structurally forced.

---

## 2. The strongest legitimate objection is “effective closure,” not “misapplied physics”

The most serious foundational objection is not that the physics is wrong. It is
that parts of it are still effective-theory choices.

Most important case:

- [stage001_geometry_lift.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage001_geometry_lift.md) explicitly says the quadratic wall action is new and not yet frozen by the parent papers.

That is acceptable if kept honest. It becomes a problem only if later stages talk
as though this action was uniquely derived from first principles.

Current assessment:

- the note is mostly honest about this;
- the move is physically reasonable;
- it is the right kind of effective interface model for the program;
- but it must stay labeled as an effective interface theory, not a solved parent
  theorem.

This is a **scope objection**, not a coherence objection.

---

## 3. The Stage `003` reduction is physically clean but branch-conditional

[stage003_bdg_coupling.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage003_bdg_coupling.md)
uses stable BdG normal coordinates and integrates them out exactly at the reduced
level.

What is good:

- the sign structure is standard;
- the wall self-energy is the exact Schur complement of the reduced stable block;
- the pole repulsion formulas are physically standard;
- no fitted coefficient is inserted to force the result.

What is still conditional:

- the real microscopic branch must actually possess the stable discrete support
  sector being assumed here;
- dangerous continuum, unstable, or negative-signature sectors must not intrude
  into the regime where the reduction is being used.

So this stage is not unphysical. It is physically sensible **provided the branch
assumptions remain explicit**.

This is the sharpest current foundational assumption.

---

## 4. The Stage `004` outgoing bridge is legitimate, but not yet unique

[stage021_reduced_one_port_normal_form.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage021_reduced_one_port_normal_form.md)
is physically one of the strongest early stages, not one of the weakest.

Why:

- it keeps the mixed `A_w/F_{mu w}/J^w` sector alive instead of cheating with an
  early brane reduction;
- it identifies exact gauge-invariant mixed observables;
- it derives the outgoing wall coefficient from a passive port plus a positive
  transfer factor;
- it does not tune the sign of the odd part by hand.

The real limitation is narrower:

- the one-lane outgoing mixed-port attachment is a reduced representative of the
  microscopic radiation channel;
- it is not yet the full grouped, multimode, unique parent-theory radiation
  theorem.

That is a legitimate future objection, but it is not evidence of misuse of
Maxwell or passivity concepts.

---

## 5. The scalar-rescue mechanism is conditional, not automatic

One thing a referee could attack later if the language becomes careless:

- [stage021_reduced_one_port_normal_form.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage021_reduced_one_port_normal_form.md)
  shows that a derivative-coupled scalar outlet can delay the dangerous scalar
  odd term from `omega` to `omega^3`.

That is a useful mechanism, but it is **not yet a no-go theorem** against all
bad scalar contamination.

So later stages should not talk as if scalar contamination is settled globally on
the basis of Stage `004` alone.

This is not a current error. It is a future overclaim risk.

---

## 6. The `001-004` stack is conservative rather than suspicious

The foundational stages are, if anything, conservative in the right direction:

- they keep mixed channels instead of suppressing them;
- they flag new ansatzes explicitly;
- they defer full nonlinear claims;
- they do not smuggle in numerical target fits;
- and they derive reduced coefficients symbolically before any later branch data
  are inserted.

That is exactly what I would want to see if the goal is “effective model with no
hidden tuning.”

---

## Non-Findings

These are things I explicitly do **not** think are foundational problems.

1. The level-set/shape-field geometry lift is not “obviously unphysical.”
   It is a standard and coherent way to expose wall/interface modes.

2. The distributed wall action is not obviously a cheat.
   It is openly presented as the first minimal effective wall theory.

3. The BdG Schur-complement reduction is not a misuse of physics.
   It is the standard reduced response of a stable quadratic coupled sector.

4. The Maxwell/mixed outgoing transfer factor is not sign-tuned.
   The positive square structure is structural.

5. The grouped real `P2` isotropy logic is not cosmetic bookkeeping.
   On an isotropic branch it is the physically right degeneracy statement.

---

## Recommended Wording Discipline

To keep the foundations referee-hard:

- say `effective wall/interface theory` rather than implying unique derivation;
- say `stable reduced support sector` rather than implying the full spectrum has
  already been classified;
- say `reduced representative of the outgoing channel` rather than implying the
  full microscopic radiation theorem is done;
- keep scalar-rescue claims conditional unless a true no-go theorem is proved.

---

## Foundation Verdict

If I were trying to kill the project as a referee, I would **not** attack
`001-004` by saying the model is tuned or that the derivation is physically
nonsensical.

I would instead attack it by saying:

- parts of the effective closure are not yet uniquely derived,
- and some reduced-sector assumptions still need stronger branch justification.

That is a much narrower and more manageable objection.

So the foundations currently pass the effective-model legitimacy test:

- **no hidden tuning found,**
- **no obvious misuse of physical principles found,**
- **but several reduced-sector assumptions must stay explicitly labeled as such.**

---

## Next Slice

The next physics red-team target should be:

- [stage129_mouth_boundary_layer.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage129_mouth_boundary_layer.md)
- [stage140_selfmatched_mouth_susceptibility.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage140_selfmatched_mouth_susceptibility.md)

That is where a referee is most likely to ask whether an effective closure has
become a disguised tuning knob.
