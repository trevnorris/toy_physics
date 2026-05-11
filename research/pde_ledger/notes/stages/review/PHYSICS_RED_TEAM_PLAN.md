# Moving-Throat PDE — Physics Red-Team Plan

## Purpose

The SymPy/Mathematica verification stack is now strong against algebraic error.
That does **not** settle the separate question a referee or skeptical physicist
will ask:

> even if the derivation is internally consistent, are the physical principles,
> closures, and reduced-sector moves actually justified?

This document turns that question into a concrete red-team audit.

The target is not only:

- symbolic correctness,
- branch bookkeeping,
- and theorem-scope wording.

The target is also:

- physical admissibility,
- variational consistency,
- thermodynamic consistency,
- gauge/boundary consistency,
- and whether the reduced model is being claimed more strongly than its physical
  derivation warrants.

This audit is **not** trying to prove ontology, meaning it is not trying to show
that the universe must literally realize this model. The practical standard is
weaker and more relevant:

- the derivation chain should make physical sense as an effective model,
- it should not rely on obviously misapplied principles,
- it should not hide ad hoc tuning to force success,
- and the parameter/closure choices should be structural rather than
  hand-adjusted to rescue isolated stages.

---

## Current Physical-Risk Assessment

Current state:

- no known algebraic blocker has been found,
- no known contradiction between SymPy and Mathematica has been found,
- no known stage-chain collapse has been found.

But there are still real referee-facing physics risks.

The relevant referee question for this project is therefore **not**

> “is this the one true ontology of nature?”

but rather

> “is this an internally coherent effective derivation, or is it using physical
> language incorrectly / tuning the model behind the scenes?”

The dominant remaining risks are:

1. **Closure ansatz risk**
   A stage introduces a physically plausible minimal closure, but it is not yet
   derived from the parent theory strongly enough to support theorem-level
   language about the actual PDE.

2. **Mode-elimination risk**
   A reduced model integrates out stable sectors or keeps only selected ports,
   but the spectral and stability assumptions that justify that reduction are
   not yet proved from the full linearized operator.

3. **Boundary/radiation risk**
   An outgoing/passive branch is represented correctly in the reduced model, but
   the map from the full parent PDE boundary conditions to that reduced outgoing
   port is not yet fully derived.

4. **Thermodynamic/kinetic risk**
   A mouth-layer equilibrium law is physically reasonable, but the exact free
   energy, mobility, conserved quantity, and zero-flux condition have not yet
   been shown to be the unique or natural reduction of the parent system.

5. **Realization risk**
   A reduced branch is characterized exactly, but it is not yet shown that the
   actual moving-throat PDE dynamically realizes that branch rather than merely
   permitting it.

6. **Hidden-tuning risk**
   A parameter value, closure relation, or branch restriction may look derived,
   but in practice may only have been selected because it makes the model work
   in one circumstance rather than as a structural consequence of the program.

---

## Effective-Model Acceptance Standard

For this project, the model should be treated as physically acceptable if all of
the following hold:

1. A physicist can see what effective principles are being used and they are not
   being misapplied.
2. The reductions and closures are clearly identified as effective-model moves,
   not disguised first-principles theorems.
3. The key parameter values and branch restrictions are fixed structurally
   rather than tuned case-by-case to salvage later stages.
4. When a stage says a certain value is required for the model to exist, that
   requirement propagates consistently through the later chain instead of being
   quietly relaxed.
5. The derivation does not succeed only because a hidden degree of freedom is
   repeatedly re-fit after the fact.

That is the red-team standard to apply below.

---

## Highest-Risk Stage Clusters

## 1. Foundations: Stages `001-004`

### Main concern

These stages contain the earliest physical modeling choices. A referee can accept
the algebra and still object that the model has not yet justified the reduction.

### Specific red-team issues

- **Stage `001`**
  The distributed wall action is explicitly introduced as a new minimal ansatz.
  That is honest, but it means the geometry PDE is not yet uniquely derived from
  the parent theory.

- **Stage `001`**
  The promotion `V_conf(X; a, L) -> V_conf(X; Sigma)` is physically natural, but
  it still needs a stronger derivation or matching argument to show it is the
  correct interface coupling rather than one reasonable choice.

- **Stage `002`**
  The two-mode breathing reduction is structurally clean, but the truncation is
  only as good as the separation between collective and neglected axisymmetric
  modes.

- **Stage `003`**
  Passing to a stable BdG normal-mode reduction is physically plausible, but a
  referee can ask whether unstable, continuum, or negative-Krein sectors have
  really been excluded on the claimed branch.

- **Stage `004`**
  The minimal localized-Maxwell + mixed-sector model is the right first reduced
  outlet model, but it is still a reduced model. The exact way the outgoing port
  attaches to the microscopic mixed sector remains a physical justification gap.

### What must be shown

- the wall action is either derived from, or explicitly declared as, the first
  controlled effective interface theory;
- the promoted confinement coupling is compatible with the parent stress/force
  law;
- the stable-mode elimination does not hide dangerous sectors on the claimed
  branch;
- the outgoing mixed port is the correct reduced representative of the parent
  radiative boundary condition.
- none of these ingredients is being quietly tuned to force the later
  normalization bridge.

---

## 2. Mouth Thermodynamics: Stages `129-140`

### Main concern

These stages are physically persuasive, but they are also exactly where a
physicist will ask whether the electrochemical mouth law is derived or posited.

### Specific red-team issues

- **Stage `129`**
  The positive source-density free energy and Onsager current produce the
  truncated exponential law cleanly, but the choice of free-energy functional,
  mobility law, and zero-flux stationary branch still needs to be defended as
  the right reduction of the parent GNLS + localized-Maxwell mouth physics.

- **Stage `129`**
  The variable normalized to one total mouth source is mathematically useful,
  but a referee may ask what conserved quantity is actually being normalized and
  whether the interval `[0,L]` carries the correct physical measure.

- **Stage `140`**
  The self-matched susceptibility closure `Theta_sigma = H_w J_s` is a sharp and
  useful simplification, but it is a closure choice. A referee can reasonably
  ask why the source susceptibility must equal the shell-layer susceptibility,
  rather than only being of the same order.

### What must be shown

- dimensional consistency of the mouth free energy, electrochemical slope, and
  normalized source law;
- that the stationary branch is an equilibrium/zero-flux branch of the intended
  parent transport law, not just a convenient profile family;
- that the self-matched susceptibility is either derivable, asymptotically
  controlled, or clearly labeled as the first effective closure.
- that the susceptibility closure is not functioning as a free tuning knob
  disguised as a derivation.

---

## 3. Core-Mouth Fixed Point: Stages `154-157`

### Main concern

These stages are mathematically coherent and numerically well-tested. The main
remaining question is physical realization.

### Specific red-team issues

- the co-evolving Boltzmann fixed-point law is exact inside the reduced closure,
  but not yet shown to be the unique physical completion of the parent PDE;
- the frozen-traction and renormalized compensated fixed points are numerically
  robust, but still reduced fixed points rather than a full dynamical theorem;
- monotonicity and uniqueness are established on analyzed windows, not globally.

### What must be shown

- what physical assumptions turn the reduced co-evolving fixed-point law into
  the actual mouth/core branch law;
- whether a Lyapunov/free-energy structure exists for the reduced fixed-point
  iteration;
- whether the parent PDE is expected to relax toward this branch, oscillate
  around it, or only pass near it transiently.
- whether the compensated branch survives because of structural identities, not
  because the closure has hidden enough flexibility to absorb any mismatch.

---

## 4. Microscopic Invariant Closure: Stages `185-187`

### Main concern

The invariant quotient structure is strong mathematics, but a physicist may still
say: this classifies a reduced coherent sector, not the realized PDE branch.

### Specific red-team issues

- Stage `185` is linearized/reference-branch rigidity.
- Stage `186` is tangent-space orbit structure.
- Stage `187` is a finite invariant-fibre theorem inside the positive coherent
  microscopic sector.

None of those by themselves prove that the actual moving-throat PDE preserves the
three invariants dynamically.

### What must be shown

- whether the preserved monomials correspond to actual conserved or adiabatically
  transported quantities of the parent system;
- whether the coherent positive sector is dynamically invariant under the parent
  evolution;
- whether the PDE branch really stays on a single similarity orbit rather than
  only the reduced finite fibre doing so.
- whether the invariant choice is structural, not selected post hoc just because
  it makes the defect collapse.

---

## Most Likely Referee Objections Right Now

These are the objections I would expect before simple algebra complaints.

1. **“You have an effective reduced model, not yet a derived full-PDE theorem.”**
   This is the strongest current objection.

2. **“The wall action and mouth susceptibility are plausible closures, but why
   are they the right physical ones?”**
   In this project that should be read more narrowly as:
   why are they legitimate effective closures rather than hidden tuning knobs?
   This is a real pressure point in Stages `001` and `140`.

3. **“Where do you prove the outgoing reduced port is the actual radiation
   channel of the microscopic model?”**
   This is the key Stage `004` objection.

4. **“Where do you prove the actual PDE realizes the selected reduced branch?”**
   This is the key Stage `154-157` and `185-187` objection.

5. **“Where are the conservation/passivity checks that tie the reduced equations
   back to the parent action or free energy?”**
   This is the global physical-consistency objection.

6. **“Are these required values structural, or were they effectively tuned to
   make the model survive?”**
   This is the project-specific anti-falsification objection.

---

## Physics Audit Workstreams

## 1. Variational Consistency Audit

Check:

- whether each reduced conservative equation can be traced back to a stated
  action or free-energy functional;
- whether wall, matter, and gauge couplings come from compatible variations;
- whether any reduced forcing term is inserted by hand rather than varied.

Deliverable:

- `notes/moving_throat/review/VARIATIONAL_CONSISTENCY_LEDGER.md`

## 2. Dimensional and Scaling Audit

Check:

- units/dimensions of every effective parameter introduced by closure;
- whether dimensionless reductions absorb measures consistently;
- whether asymptotic parameters are clearly small/large in physical terms.
- whether claimed required values such as branch or stiffness exponents are
  carried as structural consequences rather than rescued by local retuning.

Deliverable:

- `notes/moving_throat/review/DIMENSIONAL_AUDIT.md`

## 3. Passivity / Gauge / Boundary Audit

Check:

- gauge invariance of retained mixed observables;
- sign conventions for outgoing/passive odd terms;
- whether reduced port attachments correspond to proper parent boundary
  conditions;
- whether scalar monopole contamination is excluded for physical reasons, not
  only algebraic convenience.

Deliverable:

- `notes/moving_throat/review/PASSIVITY_GAUGE_BOUNDARY_AUDIT.md`

## 4. Thermodynamic Mouth-Law Audit

Check:

- whether the mouth free energy is the minimal physically admissible one or only
  one convenient choice;
- whether the normalized source law uses the right measure and conserved
  quantity;
- whether the self-matched susceptibility is controlled or only heuristic.

Deliverable:

- `notes/moving_throat/review/MOUTH_THERMODYNAMICS_AUDIT.md`

## 5. Realization Audit

Check:

- whether reduced branch existence is being conflated with actual PDE
  realization;
- whether there is a dynamical or energetic argument selecting the stated
  branch;
- where the current proof stops and where a physical completion theorem would
  begin.
- whether any apparent success actually depends on hidden local parameter
  readjustment.

Deliverable:

- `notes/moving_throat/review/PHYSICAL_REALIZATION_GAPS.md`

---

## Aggressive Red-Team Policy

For this audit, the default stance should be:

- treat every “minimal closure” as a potential referee target;
- treat every “natural” or “self-matched” identification as something that must
  either be derived, bounded, or downgraded in language;
- treat every reduced fixed point as nonphysical until its realization mechanism
  is stated clearly;
- treat every “required value” as suspect until it is shown to propagate
  structurally and not as a one-off tuning choice;
- prefer finding where the argument is weaker than advertised over defending the
  current presentation.

---

## Immediate Next Slice

If this physical-principle audit starts now, the best order is:

1. Stage `001-004`: identify exactly which ingredients are effective-theory
   ansatz versus parent-theory derivation.
2. Stage `129-140`: audit the thermodynamic mouth law and susceptibility
   closure.
3. Stage `154-157`: audit whether reduced fixed-point language outruns the
   realized physics.
4. Stage `185-187`: audit invariant-structure theorem versus actual PDE branch
   theorem.

That sequence attacks the most likely referee objections first.
