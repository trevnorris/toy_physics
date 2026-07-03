# Codex handoff: brane/bulk material-state closure

## 0. Purpose of this handoff

We had a conceptual discussion about whether the brane/bulk model should be built from fine medium constituents rather than only a smooth projection field, and whether the throat drain requires an attractive force or can be pressure/chemical-potential driven.

The conclusion is:

> The brane/bulk model should be promoted to a **material-state closure**: the brane is an ordered, shear-supporting state of the one 4D medium; the bulk is the same medium in a de-structured / non-shear-supporting state. Throats are open defects where ordered brane material can de-structure into the bulk state. Return is bulk material re-ordering into the brane state. This removes the “bulk vacuum suction” paradox and gives a cleaner upstream source for the gravity/leakage ledger.

This is not meant to replace the existing PN work. The current `STATUS.md` says the conservative PN gravity ladder through 1PN→4PN plus 2.5PN radiation is already derived/audited/GR-matched inside the calibrated hierarchy, and explicitly says not to rederive that ladder. The active push is now the full nonlinear moving-throat PDE / brane↔bulk return closure and then the brane/light sector. ([GitHub][1])

So the goal for Codex is:

> Add a brane/bulk material-state layer that sits **upstream** of the current leakage, return, moving-throat, and brane-light derivations, while preserving the already-frozen projection/reduction firewall and PN ledgers.

---

## 1. Current source anchors to keep in mind

The conceptual foundation now treats the model as one compressible 4-spatial-dimensional superfluid medium, with no point particles. A particle is an extended throat/puncture, and the native mechanisms are flow, swirl, puncture direction, brane shear, surface tension, and cohesion. ([GitHub][2])

The current conceptual picture separates the four sectors: gravity is drain/inflow, magnetism is swirl/vorticity, electric charge is puncture direction `±w`, and light is brane rotational-elastic shear. ([GitHub][2])

The same doc says a particle is a finite throat: charge is the puncture direction, mass is trapped standing-wave energy, and throat size is a balance between brane tension/backpressure and the standing-wave/drain support. ([GitHub][2])

The current projection/reduction firewall must stay intact. Projection with (W(w)) defines what the brane observer measures; reduction is a controlled limit. The exact projected leakage law and exact longitudinal identity already exist, and the zero-mode Maxwell reduction suppresses (A_w,J^w,F_{\mu w}) only in the controlled far-field brane limit, not in the microscopic ontology. 

The moving-throat wall status must also stay honest: the strict parent GNLS/Maxwell action gives a wall force through (V_{\rm conf}(\Sigma)), but not an autonomous wall PDE unless a throat action (S_\Sigma) or effective (S_\eta^{(2)}) is included.  The quadratic wall action is a consistent effective closure and produces the scalar and grouped (P_2) modal PDEs, but its coefficients must be branch data rather than post-hoc rescue parameters. 

---

## 2. Core conceptual update

### 2.1 Old picture to retire

Avoid the wording:

> The bulk is a vacuum that pulls brane material in.

That picture creates a return paradox. If the bulk is an empty suction sink, it explains brane→bulk drainage but makes bulk→brane re-entry mysterious.

Also avoid:

> The escaped fine particle attracts the next fine particle behind it.

That adds an unnecessary new force law and risks double-counting the existing pressure, tension, and chemical-potential response.

### 2.2 Replacement picture

Use:

> The brane is an ordered, shear-supporting material state of the one 4D medium. The bulk is the same medium in a de-structured / shear-free state. A throat is an open puncture where brane-ordered material can transition into the bulk state. Return is the reverse transition: bulk material re-orders into the brane state.

That means the fundamental cycle is:

[
\text{ordered brane medium}
\rightarrow
\text{open throat drain}
\rightarrow
\text{bulk/de-structured medium}
\rightarrow
\text{return/re-ordering}
\rightarrow
\text{ordered brane medium}.
]

This matches the gravity-as-drain picture and the finite-slab/return picture: the current status says the brane↔bulk return must cancel the net monopole and dipole leakage rates to avoid forbidden GR monopole/dipole radiation, with constraints like (R_0=-M_0) and (R_1=-D_1). ([GitHub][1])

---

## 3. Do the fine particles need attraction?

Codex should encode this distinction carefully:

> The fine constituents need **cohesion/order interactions** so that the brane can have surface tension and shear structure.
> They do **not** need an extra throat-specific attraction law that pulls the next particle into the mouth.

A purely noninteracting pressure gas can flow through an opening, but it cannot by itself support a stable brane/shear surface, finite mouth, surface tension, or trapped-light pressure balance. So the model still needs constituent-level cohesion or an equivalent continuum free-energy term.

However, the throat drain itself can be pressure/chemical-potential driven:

[
\text{open mouth}+\nabla p,\ \nabla\mu,\ \text{surface tension/backpressure imbalance}
\Rightarrow
\text{flow into throat}.
]

The throat is therefore an **admittance / outlet**, not a suction force. The surrounding ordered medium is pushed/relaxed into it by the existing stress gradients.

In continuum language, the “push from behind” is:

[
\mathbf J_{\rm repair}
\sim
-M_n\nabla \mu_{\rm brane},
]

not a new pairwise pull.

---

## 4. Photon/throat rule to preserve

We also clarified the photon picture.

Keep this rule:

> Photons exist only as brane-shear excitations. They do not freely propagate in the bulk because the bulk is shear-free.

The trapped standing wave in a throat is still a brane-shear mode, but the brane is extrinsically extended into the throat. So the wave can occupy throat geometry while remaining intrinsically a brane mode. The conceptual doc states this directly: light is intrinsically a brane field, its “4D-ness” is the brane’s extrinsic curvature/embedding, and the trapped mode becomes larger than the neck, falls below waveguide cutoff, and becomes bound. ([GitHub][2])

Also preserve the current correction:

> The trapped photon standing wave is **mass/support**, not gravity itself. It holds the throat open. Gravity is the separate steady de-structuring drain through the open throat.

The current conceptual doc explicitly says the throat-soliton has no sloshing, (J_w=0) is expected on the exact trapped mode, the old AC→DC rectification idea is retired, and gravity is a separate background de-structuring drain. ([GitHub][2])

The D/N trapped-support notes reinforce this: the self-reproducing standing-wave branch has zero longitudinal standing-wave current, so it is a genuine trapped support channel, not a one-way transport state. 

So do **not** write:

[
\text{photon leaks into bulk and becomes gravity}.
]

Write:

[
\text{trapped brane-shear standing wave}
\rightarrow
\text{holds throat open}
\rightarrow
\text{open throat permits medium drain}
\rightarrow
\text{drain is gravity}.
]

---

## 5. Add a brane-order state variable

The current smooth GNLS/projection variables are not enough to express the new brane/bulk material picture. Add an effective brane-order field.

Use a temporary name such as:

[
\chi_B(\mathbf X,t)\in[0,1],
]

where:

[
\chi_B=1
\quad\text{means brane-ordered / shear-supporting medium},
]

[
\chi_B=0
\quad\text{means bulk-like / de-structured / shear-free medium}.
]

If `chi` conflicts with existing (\chi_Q,\chi_0) notation, rename this to something like:

[
\mathcal O_B,\quad b_{\rm ord},\quad \sigma_B,\quad \lambda_B.
]

The exact name is less important than the split:

[
n(\mathbf X,t)
==============

\text{total conserved constituent density},
]

[
n_B(\mathbf X,t)
================

# \chi_B n

\text{brane-ordered constituent density}.
]

The current (\rho=|\psi|^2) can remain the coarse-grained parent density. If Codex needs to avoid overload, introduce (n_c) as the constituent density and state that (\rho) is its GNLS coarse-grained representation.

---

## 6. New local conservation / conversion equations

Keep total constituent conservation exact:

[
\partial_t n+\nabla_4\cdot(n\mathbf u)=0.
]

Then add an order-state balance:

[
\partial_t(\chi_B n)
+
\nabla_4\cdot(\chi_B n\mathbf u+\mathbf J_\chi)
===============================================

n\Gamma_B.
]

Here:

* (\mathbf J_\chi) is optional relative order transport. Set it to zero in the simplest advective closure.
* (\Gamma_B<0) means ordered brane material de-structures into the bulk state.
* (\Gamma_B>0) means bulk material re-orders into brane material.

In advective form:

[
D_t\chi_B
=========

## \Gamma_B

\frac{1}{n}\nabla_4\cdot\mathbf J_\chi.
]

This is the clean mathematical form of the vacation idea:

> brane→bulk is (\chi_B:1\to0); bulk→brane is (\chi_B:0\to1).

Total material is conserved. Brane-order is not.

---

## 7. Replace the old (S_{\rm leak}) interpretation with a two-part source

The old projection identity remains valid. Projection with (W(w)) gives an exact open-system brane continuity law and leakage source from (j^w). 

But once we add brane-order, the brane-observed ordered density should be:

[
\rho_B(\mathbf x,t)
===================

\int W(w),\chi_B(\mathbf x,w,t)n(\mathbf x,w,t),dw.
]

The brane-ordered current should be:

[
\mathbf J_B(\mathbf x,t)
========================

\int W(w)
\left[
\chi_B n,\mathbf u_{xyz}
+
\mathbf J_{\chi,xyz}
\right]dw.
]

Then projection gives:

[
\partial_t\rho_B+\nabla_3\cdot\mathbf J_B
=========================================

S_{\rm flux}+S_{\rm convert},
]

where

[
S_{\rm flux}
============

-\left[
W\left(\chi_B n u^w+J_\chi^w\right)
\right]*{\partial}
+
\int W'(w)
\left(\chi_B n u^w+J*\chi^w\right)dw,
]

and

[
S_{\rm convert}
===============

\int W(w),n\Gamma_B,dw.
]

This is the key formula Codex should add.

The old (S_{\rm leak}) was only projected (w)-flux. The new brane/bulk material model says the brane observer can lose ordered brane material in two ways:

1. ordered material physically moves through (w);
2. material remains nearby but loses the brane-order state.

In the old limit:

[
\chi_B=1,\qquad
\Gamma_B=0,\qquad
\mathbf J_\chi=0,
]

this reduces back to the existing projected leakage law. That recovery test should be added.

---

## 8. Gravity source reinterpretation

The gravity source should not be phrased as “bulk vacuum suction.” It should be:

[
\text{steady brane-order loss through an open throat}
=====================================================

S_{\rm flux}+S_{\rm convert}<0
]

near the defect mouth.

The quasi-static gravity hook then uses the existing exact longitudinal identity, but with ordered brane density:

[
\rho_B\nabla_3^2\phi_B
======================

S_{\rm flux}
+
S_{\rm convert}
---------------

## \partial_t\rho_B

(\nabla_3\rho_B)\cdot(\nabla_3\phi_B+\mathbf v_T).
]

In the controlled Newtonian regime:

[
\nabla_3^2\phi_B
\simeq
\frac{S_{\rm eff}}{\rho_{B0}}.
]

This preserves the existing Poisson-hook logic: the old identity already says it is exact only before the controlled near-zone reduction, and only under quasi-static/slow-density assumptions does it become a Poisson equation. 

The new material-state layer simply says what (S_{\rm eff}) physically is:

[
S_{\rm eff}
===========

\text{projected ordered-material drain}
+
\text{order-conversion loss}
+
\text{return/re-ordering terms}.
]

---

## 9. Return/re-ordering law

Define separate local pieces:

[
\Gamma_B
========

## \Gamma_{\rm return}

\Gamma_{\rm drain}.
]

Near an open throat mouth:

[
\Gamma_{\rm drain}>0
\quad\Rightarrow\quad
\Gamma_B<0.
]

In a return layer / neighboring brane / finite slab:

[
\Gamma_{\rm return}>0
\quad\Rightarrow\quad
\Gamma_B>0.
]

The global return must satisfy the existing return constraints. At minimum:

[
R_0=-M_0,
]

[
R_1=-D_1,
]

with (M_0) the net mass-rate leakage and (D_1) the net dipole/momentum-rate channel. The current status doc says these are required to avoid GR-forbidden monopole/dipole gravitational radiation. ([GitHub][1])

So Codex should keep two ledgers:

1. **local gravity drain**: produces the brane-facing field;
2. **global return/re-ordering**: restores conservation and constrains forbidden radiation channels.

Do not collapse these into one local term.

---

## 10. Free-energy closure

Codex should add a minimal effective free-energy layer. Schematic:

[
F
=

\int d^4X
\left[
\frac12 n|\mathbf u|^2
+
U(n)
+
f_B(n,\chi_B)
+
\frac{\kappa_B}{2}|\nabla_4\chi_B|^2
+
\chi_B f_{\rm shear}
+
f_{\rm throat}(R,\chi_B,n)
+
f_{\rm mix}(A_M,\chi_B,n)
\right].
]

Interpretation:

* (U(n)): existing stiff GNLS/EOS bulk energy.
* (f_B(n,\chi_B)): brane-vs-bulk ordering energy.
* (\kappa_B|\nabla\chi_B|^2/2): interface/surface-tension cost.
* (\chi_B f_{\rm shear}): shear/light only exists in brane-ordered regions.
* (f_{\rm throat}): mouth/wall coupling and order loss near throat.
* (f_{\rm mix}): mixed Maxwell/gauge contributions.

Then use chemical potentials:

[
\mu_n=\frac{\delta F}{\delta n},
\qquad
\mu_\chi=\frac{\delta F}{\delta \chi_B}.
]

A simple dissipative order relaxation would be:

[
D_t\chi_B
=========

-M_\chi \mu_\chi
+
\Gamma_{\rm throat}
+
\Gamma_{\rm return}.
]

A conservative or Hamiltonian version may require adding a conjugate order momentum or explicit heat/latent-energy ledger. Do not silently dissipate energy without accounting for it.

---

## 11. Energy ledger requirement

Every brane↔bulk conversion term must have an energy counterpart.

If (n\Gamma_B) changes brane-order, then the free-energy change is approximately:

[
P_{\rm order}
=============

\int d^4X,
\mu_\chi,n\Gamma_B.
]

This has to be balanced by:

* wall/throat work,
* pressure/enthalpy work,
* mixed Maxwell work,
* outgoing flux,
* heat/internal reservoir,
* or another explicitly declared term.

The localized Maxwell/plasma energy ledger already has the exact mixed channels:

[
\partial_tu_{\rm EM}+\partial_A S_{\rm EM}^A
============================================

-(J^aE_a+J^wE_w),
]

with (J^wE_w) an explicit mixed work channel and (S_{\rm EM}^w) Poynting transport along (w). 

Therefore, if a throat drain uses mixed-sector energy, it should explicitly touch either:

[
J^wE_w,
\qquad
S_{\rm EM}^w,
]

or a wall/order energy term. Do not create or destroy energy by redefining (S_{\rm leak}).

---

## 12. What to do with pressure vs attraction

Codex should write this as a boxed conceptual clarification:

> Fine constituents need cohesion/order interactions to give the brane surface tension and shear support. But once an open throat exists, the drain does not require a special attractive law. Flow into the throat can be driven by pressure, chemical potential, surface tension, and backpressure gradients. The throat is an admittance/outlet, not a suction force.

In discrete-particle language:

* A constituent leaving/re-ordering changes local bonds/order.
* The nearby constituents respond to local stress and chemical-potential gradients.
* The surrounding brane pushes/relaxes into the available opening.
* No “escaped particle pulls the next particle” rule is needed.

In continuum language:

[
nD_t\mathbf u
=============

## -\nabla_4 p

## n\nabla_4 V_{\rm eff}

\chi_B n\nabla_4\mu_\chi
+
\text{gauge/wall/mixed terms}.
]

The existing GNLS Madelung/Euler identity already contains the relevant structure:

[
m(\partial_t+v_j\partial_j)v_i
==============================

## q_\star(E_i+v_jB_{ij})

\partial_i(V_{\rm conf}+h(\rho)+Q).
]

So Codex should extend the potential/free-energy terms, not add a new force law. 

---

## 13. Open throat, not hard cap

Use the V2 open-exit correction:

[
R(0)=a,
\qquad
R(L)>0.
]

The D/N half-ladder is an AC support-coordinate impedance/reflection result, not evidence for a hard geometric cap. DC/background throughput exits through the open finite-radius junction and is tracked by leakage/work terms. 

This matters for the material-state picture:

* trapped standing wave: AC support mode, zero net transport;
* steady drain: DC/background throughput through the open throat;
* return: separate re-ordering/completion law.

Do not make the standing wave and the drain the same current.

---

## 14. Where this fits in the “everything should boil down to brane/bulk” idea

The right hierarchy is:

[
\boxed{
\text{constituent brane/bulk material law}
}
]

[
\Downarrow
]

[
\boxed{
\text{stationary throat branch + drain/return law}
}
]

[
\Downarrow
]

[
\boxed{
S_{\rm flux},\ S_{\rm convert},\ P_0,\ N_0,\ D_n,\ \text{return moments}
}
]

[
\Downarrow
]

[
\boxed{
\text{Newtonian gravity, PN ledgers, radiation, EM/light, atom/lepton sectors}
}
]

The existing PN work is therefore a diagnostic suite for the brane/bulk material law, not the deepest foundation. The remaining 5PN/moving-throat notes say the remaining gap is not more reduced algebra but the actual moving-throat branch solve strongly enough to extract coherent placement and outgoing normalization data. 

So yes: the majority of the fundamental physics should collapse into the brane/bulk material closure. But the closure must be explicit enough that the existing exact/reduced ledgers can be recovered as consequences.

---

## 15. Proposed doc patch: brane/bulk material-state dictionary

Add a new section to the PDE ledger, probably near the projection/reduction firewall or early moving-throat ontology:

### “Brane/bulk material-state dictionary”

| Informal phrase                             | Ledger-safe phrase                                                                   |
| ------------------------------------------- | ------------------------------------------------------------------------------------ |
| fine particles of the brane enter the bulk  | medium constituents leave the brane-ordered state and/or flow through (w)            |
| reorganize into 4D                          | ordered brane medium de-structures into bulk-like non-shear-supporting state         |
| reorganize into 3D                          | bulk medium re-orders into brane shear-supporting state                              |
| the bulk is a vacuum pulling the brane      | open-throat admittance plus pressure/chemical-potential/tension gradients            |
| vacuum must be filled with new brane        | brane-order density is restored by return/re-ordering, not creation ex nihilo        |
| escaped particle attracts the next particle | local free-energy gradients, pressure, cohesion, and surface tension drive repair    |
| photons escape into the bulk                | no; photons are brane-shear modes; trapped standing waves support throats            |
| standing wave causes gravity                | no; standing wave holds throat open; steady ordered-medium drain is gravity          |
| (S_{\rm leak}) is conversion rate           | no; (S_{\rm leak}) is the brane-observer source. Split it into flux plus conversion. |

---

## 16. Suggested exact text for the paper

Codex can paste/adapt this:

> The brane/bulk exchange should not be pictured as a vacuum suction law. In the material-state reading, the brane is an ordered shear-supporting state of the same conserved 4D medium that fills the bulk. A throat is an open defect in this ordered state. Near the throat mouth, brane-ordered material may both advect through the transverse direction and de-structure into a bulk-like state. Conversely, return corresponds to bulk-like material re-ordering into the brane state. Thus the brane observer sees an open-system source, but the parent constituent density remains globally conserved.

Then the equations:

[
\partial_t n+\nabla_4\cdot(n\mathbf u)=0,
]

[
\partial_t(\chi_B n)
+
\nabla_4\cdot(\chi_B n\mathbf u+\mathbf J_\chi)
===============================================

n\Gamma_B.
]

Projection:

[
\rho_B=\int W\chi_B n,dw,
]

[
\mathbf J_B=\int W(\chi_B n\mathbf u_{xyz}+\mathbf J_{\chi,xyz}),dw,
]

[
\partial_t\rho_B+\nabla_3\cdot\mathbf J_B
=========================================

S_{\rm flux}+S_{\rm convert},
]

[
S_{\rm flux}
============

-\left[W(\chi_B n u^w+J_\chi^w)\right]*\partial
+
\int W'(\chi_B n u^w+J*\chi^w),dw,
]

[
S_{\rm convert}
===============

\int Wn\Gamma_B,dw.
]

Then the explanatory sentence:

> The old projected leakage law is recovered when (\chi_B=1), (\Gamma_B=0), and (\mathbf J_\chi=0). The new term (S_{\rm convert}) is not a violation of total conservation; it is the projected loss or gain of brane-order within the conserved medium.

---

## 17. Minimal Codex task list

### Task A — Add ontology section

Create or update a doc such as:

```text
research/pde_ledger/notes/brane_bulk_material_state.md
```

or the nearest equivalent.

It should state:

1. one conserved 4D constituent medium;
2. brane = ordered shear-supporting state;
3. bulk = same medium in non-shear-supporting/de-structured state;
4. throat = open defect where brane-order can de-structure;
5. return = re-ordering of bulk-like material into brane state;
6. photons = brane shear only;
7. trapped throat standing wave = support/mass, not drain current;
8. pressure/chemical-potential/tension drive the drain; no extra throat attraction law.

### Task B — Add (\chi_B) balance equations

Add the equations from Sections 6–7 above. Choose a symbol that does not collide with existing (\chi_Q,\chi_0). Add a notation warning if necessary.

### Task C — Re-derive projected ordered-brane continuity

Add a short derivation:

Start with

[
\partial_t(\chi_B n)+\nabla_4\cdot(\chi_B n\mathbf u+\mathbf J_\chi)=n\Gamma_B.
]

Project with (W(w)). Integrate the (w)-divergence by parts. Verify the boundary and (W') terms. Confirm recovery of the old formula.

### Task D — Add ordered-density longitudinal identity

Repeat the Helmholtz step with:

[
\rho_B,\quad \mathbf J_B,\quad \mathbf v_B=\mathbf J_B/\rho_B.
]

Derive:

[
\rho_B\nabla_3^2\phi_B
======================

## S_{\rm flux}+S_{\rm convert}

## \partial_t\rho_B

(\nabla_3\rho_B)\cdot(\nabla_3\phi_B+\mathbf v_T).
]

State clearly:

> This is exact after declaring (\chi_B) and (W). It is not yet a Poisson equation.

### Task E — Add energy bookkeeping

Add a placeholder order-energy ledger:

[
P_{\rm order}
=============

\int \mu_\chi n\Gamma_B,d^4X.
]

Connect it to wall work, mixed EM work, heat/internal energy, and boundary flux. Cite existing (J^wE_w) and (S_{\rm EM}^w) ledgers.

### Task F — Patch language around photons

Replace any lingering text that implies photon leakage causes gravity.

Use:

> The trapped brane-shear standing wave supplies mass/support pressure and holds the throat open. It has no net trans-brane current on the exact trapped branch. Gravity is the separate steady de-structuring drain of ordered brane medium through the open throat.

### Task G — Patch language around “vacuum”

Replace “bulk vacuum pulls” with:

> bulk-like state / de-structured reservoir / open-throat admittance.

### Task H — Keep wall-action status honest

Any moving-wall PDE written in this new brane/bulk section must say whether it is:

1. promoted parent action (S_\Sigma), or
2. effective closure (S_\eta^{(2)}).

The audit says confinement-only parent action does not supply (\eta_{tt}), (\partial_w(T_w\eta_w)), or angular wall stiffness; those come only after adding the wall action. 

---

## 18. Tests Codex should add

### Test 1 — Old leakage recovery

Set:

[
\chi_B=1,\quad \Gamma_B=0,\quad \mathbf J_\chi=0.
]

Verify:

[
S_{\rm flux}
============

-\left[Wnu^w\right]_\partial+\int W'nu^w,dw.
]

This must match the existing (S_{\rm leak}).

### Test 2 — Total constituent conservation

Integrate:

[
\partial_t n+\nabla_4\cdot(n\mathbf u)=0
]

over the full 4D domain. Verify total (N=\int n,d^4X) changes only by declared external boundaries, not by (\Gamma_B).

### Test 3 — Brane-order nonconservation

Integrate the (\chi_B n) equation. Verify ordered material can change by (\int n\Gamma_B) while total material stays conserved.

### Test 4 — Return constraints

Express:

[
M_0=\int d^3x,(S_{\rm flux}+S_{\rm convert}),
]

and the dipole/momentum channel consistently. Verify the return closure can impose:

[
R_0=-M_0,\qquad R_1=-D_1.
]

### Test 5 — Photon no-bulk rule

In the brane-shear sector, check that setting (\chi_B=0) removes shear/light propagation. Bulk-like regions should not carry free photons.

### Test 6 — Standing-wave support has no DC drain

For D/N trapped support, keep the zero-current check:

[
J_z\propto\Im(\psi_j^*\partial_z\psi_j)=0.
]

The drain/source term must come from background throughput/order conversion, not this standing-wave current.

### Test 7 — Energy closure

If (\Gamma_B\neq0), require a nonzero order-energy/work term. No unledgered energy creation.

### Test 8 — No special attraction law

Check that throat flux responds to pressure/chemical-potential/wall boundary terms. No new pairwise force should appear unless explicitly justified as the underlying cohesion that also gives surface tension.

---

## 19. Acceptance criteria

The new brane/bulk material-state section is acceptable if:

1. it recovers all existing projection/leakage formulas in the (\chi_B=1,\Gamma_B=0) limit;
2. it gives a mathematically explicit source for ordered-brane loss and return;
3. it preserves global constituent conservation;
4. it does not redefine photons as bulk modes;
5. it keeps trapped standing-wave support distinct from steady drain;
6. it avoids vacuum suction language;
7. it explains why cohesion is needed for the brane but no extra throat attraction is needed;
8. it keeps the moving-wall action status honest;
9. it identifies the downstream branch-realization data still required.

---

## 20. One-paragraph executive summary for future sessions

> The brane/bulk closure should be rewritten as a material-state model. The medium has conserved fine constituents. The brane is an ordered, shear-supporting state of those constituents; the bulk is the same medium in a de-structured/shear-free state. A throat is an open defect where brane-ordered material can advect through (w) and/or convert to the bulk state. Return is bulk→brane re-ordering, constrained by the monopole/dipole return ledger. This removes the false “bulk vacuum suction” picture. Constituent cohesion is required to make brane tension and shear possible, but the throat drain itself needs no extra attractive force: pressure, chemical potential, wall tension, and backpressure gradients drive flow into the open mouth. Photons remain brane-shear modes only; trapped standing waves hold throats open but do not become bulk photons or supply the drain current. The correct source term is (S_{\rm flux}+S_{\rm convert}), not an unledgered creation/destruction law.

[1]: https://raw.githubusercontent.com/trevnorris/toy_physics/master/STATUS.md "raw.githubusercontent.com"
[2]: https://raw.githubusercontent.com/trevnorris/toy_physics/master/docs/conceptual_foundation.md "raw.githubusercontent.com"

