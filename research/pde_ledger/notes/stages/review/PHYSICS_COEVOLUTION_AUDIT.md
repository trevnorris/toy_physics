# Moving-Throat PDE — Physics Co-Evolving Branch Audit (`142-144`, `154`, `156-157`)

## Scope

This pass targets the place where a referee is most likely to accuse the model of
quietly protecting a preferred branch:

- [stage142_selfconsistent_mouth_branch.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage142_selfconsistent_mouth_branch.md)
- [stage143_equal_normalized_singular_limit.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage143_equal_normalized_singular_limit.md)
- [stage144_unique_regular_canonical_branch.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage144_unique_regular_canonical_branch.md)
- [stage154_coevolving_core_mouth_map.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage154_coevolving_core_mouth_map.md)
- [stage156_renormalized_canonical_branch.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage156_renormalized_canonical_branch.md)
- [stage157_core_mouth_coevolution_status.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage157_core_mouth_coevolution_status.md)

with supporting numerical red-team artifacts:

- [stage142_144_mouth_branch_stress.txt](/var/projects/toy_physics/research/pde_ledger/scripts/numerical/output/stage142_144_mouth_branch_stress.txt)
- [stage154_coevolving_map_stress.txt](/var/projects/toy_physics/research/pde_ledger/scripts/numerical/output/stage154_coevolving_map_stress.txt)
- [stage155_156_fixedpoint_stress.txt](/var/projects/toy_physics/research/pde_ledger/scripts/numerical/output/stage155_156_fixedpoint_stress.txt)
- [stage155_156_fixedpoint_stress.txt](/var/projects/toy_physics/research/pde_ledger/mathematica/numerical/output/stage155_156_fixedpoint_stress.txt)

The standard here is still not ontology. The standard is:

- does the co-evolving branch logic make physical sense as a reduced closure,
- is the self-consistency cost real rather than manufactured,
- is any preferred branch being protected by hidden tuning,
- and is anything being overstated beyond what the analyzed branch window actually proves.

---

## Bottom Line

Current verdict:

- the `142-144`, `154`, `156-157` chain looks **physically coherent as a reduced
  self-consistent mean-field closure**;
- I do **not** see evidence that the canonical branch is being protected by hidden tuning;
- the later co-evolving solve is actually strong anti-tuning evidence, because it
  refuses to keep the old canonical point once backreaction is allowed.

The main referee-facing limitations are narrower:

1. [stage154_coevolving_core_mouth_map.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage154_coevolving_core_mouth_map.md) is still a quasi-static reduced fixed-point closure, not a dynamical existence/stability theorem for the full PDE.
2. [stage156_renormalized_canonical_branch.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage156_renormalized_canonical_branch.md) proves uniqueness only on the analyzed positive branch window, not across every conceivable non-positive or remote branch.
3. [stage157_core_mouth_coevolution_status.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage157_core_mouth_coevolution_status.md) should always be read as a reduced-closure status note, not a theorem that the full PDE has no other microscopic realizations.

Those are serious scope limits, but they are not evidence of hidden branch bias or misapplied physics.

---

## Findings

## 1. Stages `142-144` narrow the branch space instead of rescuing it

This is the first important anti-tuning result.

Inside the explicit positive mouth family:

- [stage142_selfconsistent_mouth_branch.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage142_selfconsistent_mouth_branch.md) identifies `g_c = g_Pi`, so the core-loading ratio is no longer a free label;
- [stage143_equal_normalized_singular_limit.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage143_equal_normalized_singular_limit.md) shows the naive equal-normalized branch is not a finite regular branch at all, but a singular `Pi -> infinity` limit with divergent traction;
- [stage144_unique_regular_canonical_branch.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage144_unique_regular_canonical_branch.md) shows that within the explicit positive exponential family the lower compensated branch is the unique regular finite-bias finite-traction branch.

That is the opposite of hidden tuning.

If the model were being massaged to keep many plausible rescue branches open, this is not the structure it would produce. Instead, the closure eliminates the upper branch, eliminates the finite equal-normalized branch, and leaves only one regular compensated branch in the explicit family.

The stress harness at [stage142_144_mouth_branch_stress.txt](/var/projects/toy_physics/research/pde_ledger/scripts/numerical/output/stage142_144_mouth_branch_stress.txt) reinforces this:

- finite positive bias always stays in `2/pi < g(Pi) < 1`,
- the upper branch never enters the positive-source range,
- and the equal-normalized branch is approached only as a singular large-`Pi` limit.

---

## 2. Stage `154` is an honest co-evolving closure, not a disguised branch prescription

[stage154_coevolving_core_mouth_map.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage154_coevolving_core_mouth_map.md)
is physically the right next move once the mouth profile and core ratio are no longer treated separately.

What is good:

- the same profile `Sigma` both determines the core moment `g[Sigma]` and feels the resulting potential `Phi_{Sigma_0}[Sigma]`;
- the fixed-point law `Sigma propto exp(-Phi)` is the natural quasi-static Boltzmann/Onsager closure for the reduced mouth free energy;
- the canonical relation `R=1/4 iff g=g_*` is not inserted locally but inherited from the exact Family-1 ratio law.

So the model is not being told what branch to land on. It is being asked to solve a self-consistent nonlinear closure where the preferred branch can fail if backreaction shifts the moment.

The strongest evidence for that is numerical:

- [stage154_coevolving_map_stress.txt](/var/projects/toy_physics/research/pde_ledger/scripts/numerical/output/stage154_coevolving_map_stress.txt) shows the slope law and first-order transport identities hold for canonical, broadened, and fixed-point profiles rather than only at one specially chosen profile.

The legitimate scope limit is different:

- this is still a reduced mean-field closure, not a proof that the time-dependent PDE dynamically converges to the same fixed point.

---

## 3. Stage `156` is the strongest anti-tuning evidence in the whole mouth/core chain

If the derivation were secretly protecting the preferred canonical branch, the old
canonical traction would continue to give `g = g_*` once co-evolution is switched on.

That is not what happens.

At the old canonical traction:

- the co-evolving fixed point moves to `g_fp ~= 0.693352`,
- the loading ratio shifts to `R_fp ~= 0.282714`,
- and compensation is lost.

Restoring exact compensation requires the much larger numerically located root

- `Sigma0_can ~= 4.65103355`,
- `T_hat_can ~= 1.44670837`,
- `Pi_can ~= 3.87156438`.

Relative to the original canonical point, that is roughly:

- `+174.54%` in `Sigma_0`,
- `+60.48%` in normalized traction,
- `+173.59%` in bias.

That is not a cheap rescue. It is a substantial self-consistency cost.

This is exactly the pattern I would want to see if I were testing for hidden tuning:

- the preferred branch is not preserved automatically,
- the model has to pay a real penalty to restore it,
- and the penalty is discovered by solving the co-evolving closure rather than chosen in advance.

The fixed-point stress outputs make this stronger:

- [stage155_156_fixedpoint_stress.txt](/var/projects/toy_physics/research/pde_ledger/scripts/numerical/output/stage155_156_fixedpoint_stress.txt)
- [stage155_156_fixedpoint_stress.txt](/var/projects/toy_physics/research/pde_ledger/mathematica/numerical/output/stage155_156_fixedpoint_stress.txt)

Both show:

- seed-independence from canonical-exponential, uniform, derivative-match, and broad-mix starts,
- monotone `g_fp(Sigma_0)` on the analyzed scan window,
- and the same renormalized root across independent brackets and both CAS stacks.

That makes it very hard to argue that the renormalized branch is just a numerical artifact or a handpicked seed effect.

---

## 4. The right objection is window/globality, not hidden bias

The fair referee objection here is not “you tuned the model to keep the canonical branch.”

The fair objection is:

- uniqueness is shown on the analyzed positive branch window,
- positivity and normalization are built into the reduced closure,
- and the full PDE could in principle have other branches or instabilities outside that reduced sector.

That is a real limitation and it should stay explicit.

But it is different from:

- a bad fixed-point prescription,
- using the wrong sign conventions,
- secretly choosing the branch after seeing the answer,
- or adjusting parameters by hand to preserve the desired endpoint.

---

## 5. The co-evolving branch is physically more credible than the frozen-core correction

There is another important positive point here.

The chain does not stop at the mouth-only correction where the fixed core is still helping the canonical branch survive. It goes one step further and allows the mouth profile to feed back into the core ratio itself.

That is physically more honest, not less.

And once that more honest closure is imposed, the model becomes more demanding:

- the old canonical traction is no longer enough,
- the fixed point broadens,
- and the exact compensated branch survives only at higher traction.

That is the behavior of a model under constraint, not a model being tuned to pass.

---

## Non-Findings

These are things I explicitly do **not** think are co-evolving-branch problems.

1. The self-consistent map in [stage142_selfconsistent_mouth_branch.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage142_selfconsistent_mouth_branch.md) is not a hidden fit.
   It removes a free branch label by tying it to the source profile.

2. The singularity of the equal-normalized branch is not a contradiction.
   It is a genuine narrowing result: that branch is not a regular finite solution of the explicit family.

3. The co-evolving fixed point is not secretly initialized at the desired answer.
   The stress harnesses converge to the same profile from multiple distinct seeds.

4. The renormalized canonical point is not suspiciously close to the old one.
   It is substantially farther away, which is exactly the opposite of a protected answer.

5. The large renormalization is not, by itself, evidence of inconsistency.
   It is evidence that backreaction matters.

---

## Recommended Wording Discipline

To keep this tranche referee-hard:

- say `co-evolving reduced closure` rather than implying a full PDE dynamical theorem;
- say `on the analyzed positive branch window` whenever uniqueness or monotonicity is invoked;
- say `numerically located renormalized canonical branch` rather than suggesting a closed-form exact theorem where there is only a bracketed numerical one;
- avoid language implying that no other microscopic realization can exist outside this closure.

---

## Co-Evolving Verdict

If I were trying to kill the project as a referee, I would **not** attack
`142-144`, `154`, `156-157` by saying the branch was protected through hidden tuning.

In fact, these stages are some of the strongest evidence against that accusation:

- branch options are removed rather than multiplied,
- the old canonical point fails once real backreaction is allowed,
- and the compensated branch survives only after paying a substantial, independently cross-checked renormalization cost.

The attack I would make instead is narrower:

- this is still a reduced co-evolving closure,
- uniqueness is only shown on the analyzed positive branch window,
- and full PDE realization/stability is not yet proved.

So this slice currently passes the effective-model legitimacy test:

- **no hidden tuning found,**
- **no obvious misuse of physical principles found,**
- **strong positive evidence that the model is paying a real self-consistency cost,**
- **but theorem scope must stay limited to the analyzed reduced closure window.**

---

## Next Slice

The next physics red-team target should be the remaining foundational realization bridges:

- [stage125_positive_source_theorem.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage125_positive_source_theorem.md)
- [stage176_outgoing_load_factorization.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage176_outgoing_load_factorization.md)
- [stage185_microscopic_monomials.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage185_microscopic_monomials.md)
- [stage187_orbit_quotient_closure.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage187_orbit_quotient_closure.md)

Those are where a referee is most likely to shift from “is this tuned?” to “does the full reduced invariant structure really correspond to the intended physical branch data?”
