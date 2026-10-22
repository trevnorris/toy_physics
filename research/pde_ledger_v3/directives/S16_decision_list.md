# S16: a region around one throat, moving through its surroundings (pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-21. **Author of the build:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`. The lookups behind this list are in
`directives/_measurements/S16_setup_lookups.md` ("L-n").

**Gate.** This is a pre-builder decision list. It gets one Codex + Grok pass and is folded once (G2). It also names the
objects the build computes, so that content is reviewed until clear (spec row). The builder starts only after that.

## Filing and framing
- **The plan's step.** S16, "The worldtube reduction" (`V3_STEP_PLAN.md:670–746`). The plan says: "Rewrite S16 as a
  conditional response-side monopole approximation — ⛔ not a source-generation or interior-independence theorem —
  and list `M`, `a`, the multipoles, compactness, source normalization and boundary-flux suppression as **supplied
  inputs** that belong in S22's debt table" (`:739–741`) (L1).
- **Where the user's decisions below differ from that sentence.** The region's mass and multipoles are printed as
  moments of live fields, and compactness and boundary-flux suppression are printed as conditions. None is
  supplied, and none gets a value here, so all stay debts for S22's table. No source normalization is used, because
  no potential is supplied.
- **The framing carries over** (user, 2026-10-19). This step states what the supplied balances give. It picks no
  material law, drain law or flow. Each force states the flow it needs, and S21 compares them
  (`SUBSTRATE_REQUIREMENTS.md:729`). Under "Keep S12" (user, 2026-10-20; `steps/S12_order_source.md:15`), the
  conversion law's forms stay S12's and stay OPEN (L2).
- **The register entry `R-S14a-01`** (accepted `79511272`) stays OPEN. This step takes no slow limit (L2).

## User decisions (2026-10-21)
The user selected these options. The wording is the orchestrator's option text.
- **Boundary — "Keep the flux live":** "Integrate over a region that moves with the body and keep every term:
  momentum carried across the boundary (including the material the drain removes), stresses on the boundary, and the
  change in the enclosed mass. Print the result. Any condition under which these terms drop out is printed as a
  condition, not imposed."
- **Momentum — "General balance":** "A general local momentum balance: the rate of change of momentum density plus
  the divergence of momentum flux equals a force density. Each of the three stays a general live field, with no
  material law chosen. Supplied as an equation and labelled supplied."
- **Outside field — "No potential":** "The outside field enters only as the live flow and stress around the region.
  The build prints the force on the region. It also prints the condition under which that force takes the form of a
  gradient times a mass, without assuming it. Inertial and passive mass are printed as separate objects."

## Supplied equations (unfalsifiable in this build)
On S14a part 1's projected domain: in-plane coordinates `x^i` (three) and time `t`, with the slab integrated in `w`
between its live faces.
- **The projected material balance,** as the accepted part 1 record prints it (`steps/S14a_1_projection.md:91`)
  (L3):

  ```
  ∂_t Σ + ∂_i F^i − n₁a₁,₋ + n₃a₃,₊ = 0
  ```

  - `Σ` and `F^i` are the projected number density and in-plane flux (`:84`).
  - The face terms are the material leaving through the faces (`:87`), with the material face factor `a_{r,s}` of
    `:62`.
  - It holds where the committed material law holds; the record's note at `:91` says so. Every field stays live.
- **The general projected momentum balance** (the user's "General balance"):

  ```
  ∂_t G^i + ∂_j Π^{ij} = f^i
  ```

  `G^i` is the projected momentum density, `Π^{ij}` its flux and `f^i` the force density, all general live fields of
  `x` and `t`. `f^i` includes every exchange of momentum with the faces and the bulk. `G^i` is not identified with
  any multiple of `F^i`, and `Π^{ij}` is not assumed symmetric.
- **`m`, a supplied constant mass per unit of `n`.** `n` is the material number density
  (`directives/S12_decision_list.md:45`), so `m Σ` is the projected mass density.

## The region
- **The throat's centre** is a live point `X_c(t)`. The excluded core lies inside an inner surface that moves with
  it. The core's content stays with S22.
- **The region** lies between that inner surface and an outer surface around it, and both surfaces move. The build
  may take them as spheres about `X_c(t)` with live radii `r_in(t)` and `R(t)`. That is a choice of control surface,
  not a symmetry of the fields, and it is labelled where printed. `r_in` is not identified with the mouth radius.
- **The fields are the total fields.** Items 3 and 4 need a field of the surroundings. For it, the build may write a
  field on the region as the sum of two parts:
  - a part that is regular across the region, the core included;
  - the rest.

  Both parts stay live, and neither gets a law. The rest's contributions stay printed. This is a labelled choice.

## What the build computes
Both engines construct these independently, with the `.wl` blind.

1. **The region's balances.** For the region's material mass `M(t)`, its centre of mass `X(t)` and its momentum
   `P(t)`, defined from `m Σ` and `G^i` over the region: their time rates, from the supplied balances, with the
   motion of both surfaces kept. The build prints separately:
   - at each surface, the material flux and the momentum flux relative to that surface's motion;
   - the face terms and the force density, integrated over the region.
2. **The centre-of-mass equation.** The relation the balances give for the acceleration of `X`, with every term
   printed. These include:
   - the change of `M`;
   - the motion of both surfaces;
   - the difference between `P` and `M dX/dt`;
   - the relation between `X` and `X_c`.

   No term is dropped. The factor that multiplies the acceleration is printed as the **inertial factor**. A factor
   fixed by the definitions of `X`, `M` and `P` alone is printed as defined, not as derived.
3. **Moments.**
   - The region's material moments about `X_c` through second order: the mass, the dipole and the second moment.
     `a_WT` is the width that the second moment defines. It is a different object from `r_in` and from the mouth
     radius, and none is substituted for another.
   - Each term of item 1 that involves the regular part is expanded about `X_c` through second order. The remainder
     is printed as a term.
4. **The force and its form.** The force on the region is every term of item 2's relation other than the inertial
   term.
   - Suppose its leading term is a scalar of the region times the gradient, at `X_c`, of a field of the regular part.
     The build prints that scalar as the **passive factor**.
   - It prints, as equations on the fields, the condition under which that form holds.
   - It prints the passive factor and the inertial factor separately, with their difference.
5. **Dimensional homogeneity** of items 1–4.

A name may name an object, never its value or a comparison outcome. The script prints objects and never states
conclusions. No acceptance criterion references an expected value.

## Model point
- **Regime:** none chosen. Time dependence is kept, and no steady, quasi-static or slow limit is taken.
- **Freezes:** none intended. Both surfaces move, and the region's mass may change. None of the following is imposed;
  each appears only as a printed condition:
  - compactness;
  - a vanishing dipole;
  - boundary-flux suppression;
  - a potential;
  - a symmetry of the fields.
- **Labelled choices:**
  - the control surfaces' shape;
  - the split into a regular part and the rest;
  - item 3's second-order truncation, with its remainder printed.
- **Transfer:** wherever part 1's projected balance holds, outside the excluded core.

## Inventories (for the record, under its own review; not computed, and not in the build directive)
The record points at each source row without restating it:
- **The plan's S16 section,** `V3_STEP_PLAN.md:670–746`. It covers the source statement and its transcription note
  (`:700–703`), the gravity-sector guardrail (`:705–730`) and the A11 gate (`:732–737`) (L1).
- **The source's assumptions:** `../4d_1pn_full/paper/4d_1pn_full.tex:695–720`.
  **Its reduced equation:** `:885–891`.
  **The 2.5PN statement as written:** `../4d_2_5pn/paper/4d_2_5pn.tex:605–614` (L4).
- **The charter's scope boundary:** `CHARTER.md:74–80` and `:109–119` (L5).
- **The S22 debt rows** for `a_WT`, the `a_WT`/`a_mouth` bridge, and `M` with its multipoles:
  `V3_STEP_PLAN.md:1303–1305` (L5).
- **The mass gap:** C3, `DEFECT_REGISTER.md:171` (L5).
- **O2's S16 interface row:** `directives/O2_SHARED_PHYSICS.md:607`. **Its scoping row:**
  `directives/O2_steady_brane_balance_scoping.md:164` (L6).

## Retained OPEN
- **No owner assigned here:** the forms of `G^i`, `Π^{ij}` and `f^i`, and the slab's material law. The accepted S12
  record gives slab material completion to S8/S22 (`steps/S12_order_source.md:292`) (L6).
- **No owner named:** the faces' law (same row).
- **S22:**
  - the core's interior, mass and momentum;
  - the `a_WT`/`a_mouth` bridge, under A11.
- **`R-S14a-01`:** the slow limit and far behaviour.
- **S15, conditional on S14a:** the source side, meaning what force a drain exerts.
- **C3, a gap:** the link between the drain rate and mass.
- **The C10 tail question:** the plan's S19 (`V3_STEP_PLAN.md:766–788`).
- **The plan's deferred branch** (`:705–730`), not run here: species-dependent response.

## Intended claim
What the supplied balances give for a region around one throat that moves through its surroundings, with every term
live:
- the region's balances;
- its centre-of-mass equation;
- its moments;
- the force on it;
- its inertial and passive factors, kept separate.

The conditions under which the point-particle form would follow are printed as conditions. The step gives no mass,
no potential, no force law between bodies, and no throat solve.

## Next
A build directive for Codex, derived from this list. Then:
- G1 review of the build (fresh Claude + Grok), until clear;
- a cross-engine comparison;
- the record and its register pass, each reviewed until clear.
