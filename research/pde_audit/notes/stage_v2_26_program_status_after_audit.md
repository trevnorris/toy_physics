# Stage V2-26 - Program Status After the PDE Audit

## Purpose

This note records the honest program status after the PDE audit and the first
target-blind simulation layer.

The audit was useful because it did not merely check algebra.  It separated:

- what has been derived cleanly;
- what is conditional or reduced;
- what has executable support;
- what the current simulations fail to realize;
- and what physical derivations are still missing before the PDE can be
  claimed as a completed one-equation account.

This note is intended to be paper-facing.  It should prevent the audit result
from being summarized too strongly or too weakly.

For the physical ontology that should accompany this status statement, see
`stage_v2_28_physical_picture_and_ontology.md`.  That note records the
finite brane-bulk throat/puncture picture, the open-conduit condition, the
mouth/interior distinction, flux distinctions, and common misreadings to avoid.

For the superfluid material-closure gap, see
`stage_v2_29_superfluid_material_closure_gap.md`.  That note records the
difference between existing reduced EOS/sound-speed formulas and a completed
branch-level derivation of density, sound speed, effective light-speed
behavior, and flux feedback.

---

## 1. What can be said solidly

The current project contains real mathematical work.

The audit verifies many nontrivial pieces of the reduced/effective framework:

- parent-action versus effective-wall status;
- Maxwell localization and gauge normalization;
- dimension and port/source normalization ledgers;
- open finite-radius throat boundary logic;
- Poisson/Newtonian hook;
- stable BdG support Schur complement;
- Maxwell/mixed conservative and outgoing kernels;
- Hamiltonian positivity gates;
- grouped real `P2` projector calculus;
- STF angular source-map normalization;
- grouped normalization and constant-prefactor branch algebra;
- outgoing `l=2` fingerprint;
- 2.5PN/4PN interface and tail transport gate;
- branch-freeze/no-refit discipline;
- weak-axisymmetric splitting;
- similarity-orbit separation;
- isotropic full-bundle target surface;
- weak-form and branch-extraction scaffolding;
- fixture-backed solver handoff and negative controls;
- reduced FEM and manufactured nonlinear simulation checks.

These are not just narrative assertions.  Many of them are backed by executable
Python audits, Mathematica mirrors, fixtures, generated artifacts, and
target-blind guards.

So the proper positive claim is:

```text
The program has derived a growing reduced/effective framework whose algebraic
targets, adapters, stability gates, and no-refit protocol are internally
consistent and executable under stated assumptions.
```

---

## 2. What cannot be claimed yet

The audit does not prove that known physics has been derived from one completed
PDE.

The PDE is still incomplete as a physical branch-realization system.  In
particular, the current audited exporter does not yet derive one frozen
moving-throat branch that simultaneously outputs the required conservative,
outgoing, and orbit-lock data.

The following claim would be too strong:

```text
Known physics has been derived from a single PDE.
```

The current honest replacement is:

```text
Several reduced/effective pieces of known-physics structure have been derived
or targeted under explicit assumptions, and the audit now identifies the exact
branch-output conditions that a completed PDE must satisfy.
```

The distinction matters.  The target algebra can be correct while the actual
PDE branch remains unsolved.

---

## 3. Why the clean algebra and failed simulations are not contradictory

The algebra answers a conditional question:

```text
If a moving-throat branch is going to reproduce the desired GR-like packet,
what coefficient relations must it satisfy?
```

The simulations answer a different question:

```text
Do the currently exportable reduced/manufactured branch families land on that
target surface without target feedback?
```

The current answer to the second question is no.

The referee run gives:

```text
0/192 target-passing reduced frozen candidates.
0/3 target-passing manufactured nonlinear candidates.
```

Among reduced open-stable candidates, the one-pole ratio

```text
D0*C/(3*A^2)
```

has maximum `0.1353664855760648`, far below the target value `1`.

The required-deformation diagnostics show that the miss is not a small local
continuation:

```text
minimum reduced C-or-D0 multiplier: 7.387352901601946
median reduced C-or-D0 multiplier: 16.30132163440465
median reduced P0 multiplier: 171.65261223353198
```

The projection-stress diagnostic then shows that even fixing one-pole support
and uniform outgoing amplitude after the fact is insufficient.  A successful
branch needs outgoing moment-shape control through `N2` and `N4`, not just a
scalar multiplier.

So the simulation miss means:

```text
The current exported branch families are not the desired physical branch.
```

It does not mean:

```text
The target algebra is inconsistent.
```

It also does not mean:

```text
The desired branch exists.
```

That existence question remains the central unsolved problem.

---

## 4. What the audit says is missing

The audit has turned a broad uncertainty into a concrete missing-physics list.

The completed branch must realize, on one frozen branch:

```text
D0 = K - B0 - Z0,
A  = M + B2 + Z2,
C  = B4 + Z4,
D0*C/(3*A^2) = 1,
P0 = N0/D0,
P2 = 0,
P4 = 0.
```

Equivalently, the outgoing transfer moments must satisfy

```text
N2 = -2*A*N0/D0,
N4 = N0*(A^2 - 2*D0*C)/D0^2.
```

On the one-pole surface this becomes

```text
N4 = -5*A^2*N0/D0^2.
```

The 5PN/orbit-lock side additionally requires the actual branch to satisfy:

```text
dln_R_tr = 0,
dln_R_target = 0,
dln_epsilon_eta = 0,
N_Q = 1.
```

The current materials do not yet provide the physical exporter that determines
all of these quantities from the throat equations.

---

## 5. Physical mechanisms not yet represented

The five audit readouts `D0`, `C`, `P0`, `N2`, and `N4` are not a full physical
description of the throat.  They are low-order projected response moments.

They do not directly encode:

- the actual throat geometry `R0(s)`;
- wall shape and axial boundary-layer structure;
- superfluid intake and mass balance;
- DC leakage and open/radiative mouth flux;
- reservoir depletion or backreaction;
- nonlinear saturation and branch selection;
- full BdG/Krein mode signatures;
- complete Maxwell/mixed port spectra;
- finite-throat `P0/P2/P22` mouth geometry;
- intrinsic `P22` bracing and orientation;
- half-flux/mixed-sector closure;
- orbit-lock tangent data;
- thermodynamic or heat/export partition;
- superfluid material closure: EOS, stationary density, `c_s(rho)`, effective
  light-speed behavior if density dependent, and flux feedback;
- higher dissipative PN information;
- the physical origin of source/port normalization.

Those mechanisms may affect the five readouts, but the readouts alone do not
derive the mechanisms.

This is why the final PDE is harder than a coefficient-fitting problem.  It
must physically select the branch whose projected response moments satisfy the
audit packet.

---

## 6. What the unincorporated notes contribute

The notes-intake audit found real candidate ingredients:

- finite-throat `P0/P2` forcing from Hessian/tidal mouth loading;
- scalar open/radiative `P0` flux hooks;
- intrinsic `P22` bracing from mixed-sector/half-flux structure;
- non-rigid `U/V` dressing;
- leakage/work lanes;
- microscopic odd passive/export kernels;
- support/source branch data that is no longer the active bottleneck in the
  reduced 5PN stack.

These are useful.  They should be carried into the next derivation phase.

But they are not yet a completed exporter.  They do not yet provide a calibrated
map:

```text
throat/source physics -> D0, C, P0, N2, N4
```

and they do not yet prove that the same branch satisfies the orbit-lock packet.

So the notes are best understood as:

```text
source physics and protocol constraints for the next branch derivation,
not a solved target-passing PDE branch.
```

---

## 7. Publication framing

The audit should be published as a boundary-setting document.

It should claim:

- the algebraic and executable audit layers are internally consistent;
- the reduced/effective claims have been separated from strict parent-PDE
  claims;
- the current target-blind branch families fail the target;
- the miss is quantitatively diagnosed;
- post-hoc rescue/tuning is explicitly disallowed;
- the remaining PDE completion problem is now sharply specified.

It should not claim:

- a completed derivation of known physics from one PDE;
- an actual physical branch that already passes the target packet;
- that support tuning alone can explain the simulation miss;
- that scalar outgoing normalization alone can repair the branch;
- that algebraic projection after seeing residuals is physical evidence.

Recommended paper language:

```text
The present audit does not complete the moving-throat PDE program.  It
establishes a reproducible reduced/effective ledger, identifies the exact
coefficient packet required of any successful branch, and shows that the
currently exportable branch families do not realize that packet.  The remaining
problem is the derivation and target-blind export of an actual moving-throat
branch containing the missing intake, finite-mouth, outgoing, and orbit-lock
physics.
```

---

## 8. Remaining derivation program

The next phase should be derivation-first and gap-driven.

Required derivations:

1. Decide the strict parent/effective status of the wall/throat dynamics.
   Either promote `S_Sigma` or keep an explicit effective closure label.

2. Derive the actual stationary open-throat branch equation, including the
   physical terms admitted from the finite-mouth, scalar-flux, `P22`, `U/V`,
   leakage/work, and export-kernel notes.

3. Derive the linearized wall/BdG/Maxwell/mixed/outgoing-port exporter that
   maps the branch to `K,M,B_n,Z_n,N_n`.

4. Derive the source/port normalization law that produces `P0` as branch data
   rather than as a fitted scale.

5. Derive independent outgoing moment-shape control for `N2` and `N4`.

6. Derive the weak-axisymmetric tangent and orbit-lock packet on the same
   branch.

7. Freeze the complete protocol before target evaluation.

8. Run the unchanged V2-22B -> V2-22A -> V2-21 residual chain and report the
   result without refitting.

Higher PN work may still be useful, especially dissipative half-integer orders
such as 3.5PN, 4.5PN, and 5.5PN, because they may constrain outgoing
moment-shape physics.  But more PN algebra is secondary until the actual branch
exporter exists.

---

## Current status statement

```text
Real mathematical audit work: yes.
Real reduced/effective physics bridges: yes, under stated assumptions.
Completed one-PDE derivation of known physics: no.
Current physical exporter: incomplete.
Current target-blind simulated branch families: fail.
Main remaining problem: derive and export the actual moving-throat branch.
```

The audit therefore moves the program forward by changing the question from

```text
Can the algebra be arranged to look right?
```

to

```text
Can the completed PDE physically select a stable open branch whose frozen
projected response packet satisfies the audited target conditions?
```
