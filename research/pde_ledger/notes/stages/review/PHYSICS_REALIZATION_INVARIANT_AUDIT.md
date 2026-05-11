# Moving-Throat PDE — Physics Realization and Invariant Audit (`125`, `170`, `185-187`)

## Scope

This pass targets the remaining referee-facing bridge between reduced symbolic
closure and physically interpretable branch data:

- [stage125_positive_source_theorem.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage125_positive_source_theorem.md)
- [stage170_linear_grouped_outlet_map.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage170_linear_grouped_outlet_map.md)
- [stage185_microscopic_monomials.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage185_microscopic_monomials.md)
- [stage186_similarity_orbit_closure.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage186_similarity_orbit_closure.md)
- [stage187_orbit_quotient_closure.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage187_orbit_quotient_closure.md)

with supporting red-team artifacts:

- [stage058_170_foundational_stress.txt](/var/projects/toy_physics/research/pde_ledger/scripts/numerical/output/stage058_170_foundational_stress.txt)
- [stage185_187_orbit_stress.txt](/var/projects/toy_physics/research/pde_ledger/scripts/numerical/output/stage185_187_orbit_stress.txt)

The standard here is:

- does the branch-selection/invariant logic make physical sense as a reduced model,
- are the remaining “orbits” and “invariants” being interpreted correctly,
- is there any sign that the model is still hiding tuning inside these reduced variables,
- and is the theorem scope kept separate from the still-open full PDE realization question.

---

## Bottom Line

Current verdict:

- the `125`, `170`, `185-187` chain looks **physically coherent as reduced
  branch-selection and invariant-structure analysis**;
- I do **not** see hidden tuning in this tranche;
- I do **not** see a physically nonsensical use of positivity, grouped anisotropy,
  or invariant/quotient language.

The main referee-facing limitations are narrower:

1. [stage125_positive_source_theorem.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage125_positive_source_theorem.md) is a theorem only inside the one-lane positive localized-source closure on the first D/N interval.
2. [stage170_linear_grouped_outlet_map.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage170_linear_grouped_outlet_map.md) is a linear-response theorem on the compensated isotropic branch, not a full anisotropic branch solution.
3. [stage185_microscopic_monomials.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage185_microscopic_monomials.md) through [stage187_orbit_quotient_closure.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage187_orbit_quotient_closure.md) are reduced invariant-structure theorems. They do not yet prove that the full PDE branch dynamically preserves those invariants.

Those are scope limits, not signs of bad physics or hidden rescue parameters.

---

## Findings

## 1. Stage `125` is a clean positivity-selection theorem, not an ad hoc branch choice

[stage125_positive_source_theorem.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage125_positive_source_theorem.md)
does something physically reasonable and quite restrictive:

- it takes a nonnegative normalized axial source density,
- projects it against the first cosine moment on the first D/N interval,
- and uses positivity of the weight and kernel to force `0 <= g[sigma] <= 1`.

That is not arbitrary. It is a direct consequence of positivity plus the mode shape on `[0,L]`.

The physical conclusion is therefore strong and clean:

- the upper compensated Family-1 branch is excluded structurally,
- the lower branch survives as the only admissible compensated branch inside this positive one-lane source model.

This is anti-tuning evidence. The branch is not chosen because it is desirable; it is chosen because positivity leaves no other compensated option in the stated closure.

The fair objection is simply:

- this theorem does not cover sign-changing, multimode, or nonlocalized mouth data.

That is real, but it is not a misuse of physical principles.

---

## 2. Stage `170` is a legitimate linear-response compression of grouped anisotropy

[stage170_linear_grouped_outlet_map.md](/var/projects/toy_physics/notes/moving_throat/moving_throat_pde_stage170_linear_grouped_outlet_map.md)
is physically a linearized symmetry/compression statement on the compensated isotropic branch.

What is good:

- Stage `169` has already removed scalar off-bundle slippages for pure grouped real `P2` anisotropy at linear order;
- the remaining linear problem is correctly reduced to the direct outlet coefficients `delta kappa_W` and `delta gamma_W`;
- the grouped microscopic bundle data then collapse to two exact combinations `(K_A, G_A)`, one even and one odd.

That is exactly what I would expect from a clean linear-response decomposition:

- even data feed the hidden even pole deformation,
- odd normalization data feed the hidden odd outgoing normalization,
- and consistency of the even channel imposes an explicit one-parameter compatibility law.

The stress harness at [stage058_170_foundational_stress.txt](/var/projects/toy_physics/research/pde_ledger/scripts/numerical/output/stage058_170_foundational_stress.txt) is useful here:

- two different grouped perturbation sets with the same `(K_A,G_A)` collapse to the same outlet coefficients,
- while an inconsistent even lane fails the hidden-even relation.

That is strong evidence that the map is genuinely compressing structured branch data, not hiding a fit parameter.

The fair limitation is narrower:

- this is a linearized map about the compensated isotropic branch, not yet a theorem about fully nonlinear anisotropic branch realization.

---

## 3. Stages `185-187` are mathematically strong reduced similarity theorems, not fake physics

This is the tranche where a referee could complain that the notes are sliding from
physics into arbitrary algebra.

I do not think that criticism is fair, provided the scope stays explicit.

What these stages actually do is:

- identify three direct microscopic monomials whose first log drifts are the branch-defect coordinates;
- show that the zero-defect ledger is equivalent, at first grouped weak-axisymmetric order, to preserving those three monomials;
- then upgrade that tangent-space statement to a finite similarity-orbit / quotient statement inside the positive coherent microscopic sector.

That is a clean reduced invariant-structure result.

The physically correct way to read it is:

- many microscopic co-scalings are redundant **for the reduced branch data being tracked**,
- the three monomials are the complete reduced invariants of that similarity action,
- and defect motion is exactly motion in those reduced invariant coordinates.

That is not “making up gauge symmetry.” It is an exact classification of the reduced coherent hierarchy being used.

The important wording discipline is:

- these are similarity directions in the reduced model,
- not a claim that the underlying PDE has a fundamental gauge redundancy in the same literal sense.

---

## 4. The orbit/invariant structure is being checked in the right way

The numerical stress harness at [stage185_187_orbit_stress.txt](/var/projects/toy_physics/research/pde_ledger/scripts/numerical/output/stage185_187_orbit_stress.txt) gives the kind of evidence I would want here.

Across interior and edge samples it shows:

- tangent directions lie in `ker(M_*)`,
- tangent directions preserve the invariant triple to first order,
- transverse directions produce genuinely nonzero invariant motion,
- and finite orbit pairs preserve `(C_tr, C_nt, epsilon_eta)` exactly while obeying the solved finite laws for `Delta_eta`, `Delta_T`, and `Delta_mu`.

That matters physically because it distinguishes:

- true reduced zero-cost similarity motion,
- from real defect-producing motion.

So the quotient language is not just aesthetic packaging. It is identifying which microscopic deformations are invisible to the reduced branch data and which ones are not.

---

## 5. The remaining risk is realization, not incoherence

The strongest remaining objection is not that these stages are unphysical.

It is that they stop at the reduced coherent hierarchy.

The still-open question is:

- does the actual moving-throat PDE branch preserve those three invariants closely enough that the reduced similarity-orbit classification is the right physical branch language?

That is a serious question, but it is exactly the question the notes now state openly.

It is not evidence that:

- positivity was misused in Stage `125`,
- grouped anisotropy was decomposed incorrectly in Stage `170`,
- or the invariant coordinates in `185-187` were introduced to hide tuning.

---

## Non-Findings

These are things I explicitly do **not** think are problems in this tranche.

1. Stage `125` is not secretly choosing the lower branch by hand.
   Positivity excludes the upper branch inside the stated closure.

2. Stage `170` is not compressing anisotropy arbitrarily.
   The collapse to `(K_A,G_A)` is exactly what the linear outlet map and even-consistency law enforce.

3. Stages `185-187` are not claiming the full PDE has been solved.
   When read correctly, they classify the reduced coherent invariant structure only.

4. The similarity orbit is not a hidden tuning family.
   It is the family of microscopic co-scalings that leave the tracked reduced invariants unchanged.

5. The quotient theorem is not physically meaningless.
   It sharply identifies which microscopic motions matter to the defect ledger and which do not, inside the reduced sector.

---

## Recommended Wording Discipline

To keep this tranche referee-hard:

- say `one-lane positive localized-source closure` for Stage `125`;
- say `linear grouped outlet map on the compensated isotropic branch` for Stage `170`;
- say `reduced coherent similarity orbit` rather than suggesting a fundamental full-PDE gauge symmetry in `186-187`;
- keep repeating that the open problem is full PDE branch realization, not reduced invariant classification.

---

## Realization / Invariant Verdict

If I were trying to kill the project as a referee, I would **not** attack
`125`, `170`, `185-187` by saying the branch logic is physically nonsensical or tuned.

I would attack them more narrowly by saying:

- the Stage `125` positivity theorem is only as broad as its one-lane positive-source closure,
- the Stage `170` anisotropy map is only linearized about the compensated isotropic branch,
- and the Stage `187` quotient theorem is still a reduced coherent invariant theorem rather than a full dynamical PDE branch theorem.

That is a much narrower objection than “bad physics.”

So this tranche currently passes the effective-model legitimacy test:

- **no hidden tuning found,**
- **no obvious misuse of physical principles found,**
- **the reduced branch-selection and invariant structure are physically interpretable,**
- **but the final realization question remains whether the actual PDE branch stays close to this reduced invariant hierarchy.**

---

## Remaining Gap

At this point the remaining red-team gap is sharply focused:

- not more algebra,
- not more CAS cross-checking,
- not branch tuning,
- but the final realization question of whether the full moving-throat PDE branch tracks the reduced coherent invariant structure closely enough for referee purposes.

That is the place where any further hardening should go.
