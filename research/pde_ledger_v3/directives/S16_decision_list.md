# S16: a region around one throat, moving through its surroundings (pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-21. **Author of the build:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`. The lookups behind this list are in
`directives/_measurements/S16_setup_lookups.md` ("L-n").

**Revision 2: the G2 fold.** Pass 1 (Codex + Grok) is preserved at `ffaa2824`, with its dispositions D1–D6 and the
user's two decisions on them in `directives/_measurements/S16_decision_list_review_p1_disposition.md`. The fold
changes what is computed, so its spec content is reviewed until clear before any builder.

## Filing and framing
- **The plan's step.** S16, "The worldtube reduction" (`V3_STEP_PLAN.md:670–746`). The plan says: "Rewrite S16 as a
  conditional response-side monopole approximation — ⛔ not a source-generation or interior-independence theorem —
  and list `M`, `a`, the multipoles, compactness, source normalization and boundary-flux suppression as **supplied
  inputs** that belong in S22's debt table" (`:739–741`) (L1).
- **Where the user's decisions differ from that sentence.**
  - The region's material count and its moments are printed from live fields.
  - The boundary fluxes are printed, each as its own object.
  - No mass is defined, so `M` stays an input the reduction would need.
  - Compactness and a field of the surroundings are filed as requirements, because the region's fields do not
    determine them (D1).
  - No source normalization is used, because no potential is supplied.

  All of these stay debts for S22's table.
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
- **After pass 1, on passive mass — "Record what it needs":** "The build prints the momentum the region receives:
  stresses and momentum carried across both surfaces, plus the force density. It prints this apart from the terms
  caused only by the surfaces moving. It computes no passive factor. What one needs (a named field of the
  surroundings, a rule separating it from the body's own field, and its normalization) goes to the register as a
  requirement. The missing owner is stated as missing."
- **After pass 1, on mass — "Count material; mass open":** "The region's content and centre come from the material
  count, which the supplied balance governs. No mass is defined. How the region's momentum relates to its motion
  (inertial mass) is printed as open. The law it needs goes to the register, as O2 keeps it (S8)."

The two later decisions narrow the last two sentences of "No potential". The gradient-form condition and the two
masses are filed as requirements, not computed. "The change in the enclosed mass" becomes the change in the
enclosed material count.

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
  - It holds where the committed material law holds, as the record's note at `:91` says. Every field stays live.
- **The general projected momentum balance** (the user's "General balance"):

  ```
  ∂_t G^i + ∂_j Π^{ij} = f^i
  ```

  - `G^i` is the projected momentum density, `Π^{ij}` its flux and `f^i` the force density. All three are general
    live fields of `x` and `t`.
  - `f^i` includes every exchange of momentum with the faces and the bulk.
  - `G^i` is not related to `F^i` or `Σ` by any supplied law, and `Π^{ij}` is not assumed symmetric.

## The region
- **The throat's centre** is a live point `X_c(t)`. The excluded core lies inside an inner surface that moves with
  it. The core's content stays with S22.
- **The region** lies between that inner surface and an outer surface around it, and both surfaces move. The build
  may take them as spheres about `X_c(t)` with live radii `r_in(t)` and `R(t)`. That is a choice of control surface,
  not a symmetry of the fields, and it is labelled where printed. `r_in` is not identified with the mouth radius.
- **The fields are the total fields on the region.** No object of the build depends on values inside the core or on
  a split of a field into parts (D1).

## What the build computes
Both engines construct these independently, with the `.wl` blind.

1. **The region's balances.** These are the rates of three objects, from the supplied balances, with the motion of
   both surfaces kept:
   - the material count `N(t)`, the integral of `Σ` over the region;
   - its first moment about `X_c`;
   - the momentum `P(t)`, the integral of `G^i`.

   For each rate, the build prints separately the face terms and the integrated force density. At each surface,
   it prints the material flux and the momentum flux. Each surface term is printed in two parts, kept apart:
   - the part present when that surface is at rest;
   - the part due to that surface's motion.
2. **The centre of the material count.** `X(t)` is the first moment over `N`, about `X_c`. The build prints the
   relations the material balance gives for the rate and the acceleration of `X`, with every term printed. These
   include:
   - the change of `N`;
   - the motion of both surfaces, with the terms that come only from that motion kept apart;
   - the relation between `X` and `X_c`.

   Any rate that no supplied balance governs stays a printed term.
3. **The momentum the region receives.** This is the rate of `P`, from item 1: the stresses and the momentum carried
   across both surfaces, plus the integrated force density. It is printed as its own object, apart from item 2's
   kinematic terms and from the parts due to the surfaces' motion.
4. **Moments.** The material count's moments about `X_c` through second order: the count, the dipole and the
   second-moment tensor. No scalar width is defined from the tensor (D5).
5. **Dimensional homogeneity** of items 1–4.

Each object normalized by `N` is stated on the domain where `N` is nonzero. The unnormalized balances and moments
keep the full domain (D6).

The build computes no inertial factor, no passive factor, no mass and no field of the surroundings. Items 1–3 print
`P` and the material-count objects side by side, each from its own balance. No relation between them is supplied or
derived.

A name may name an object, never its value or a comparison outcome. The script prints objects and never states
conclusions. No acceptance criterion references an expected value.

## Model point
- **Regime:** none chosen. Time dependence is kept, and no steady, quasi-static or slow limit is taken.
- **Freezes:** none intended. Both surfaces move, and the region's count may change. Nothing is imposed: no
  compactness, vanishing dipole, boundary-flux suppression, potential, mass-to-count ratio, or symmetry of the
  fields. Each boundary term is printed as its own object, so the condition under which it drops out is that the
  object vanishes.
- **Labelled choices:**
  - the control surfaces' shape;
  - moments printed through second order.
- **Transfer:** wherever part 1's projected balance holds, outside the excluded core.

## For the record and its register pass (not in the build directive)
- **Requirements the user's decisions file.** The record files each requirement once, against the owner the sources
  assign. Where no owner is assigned, it says so; no owner is invented.
  - **For a passive factor:** a named field of the surroundings, a rule separating it from the body's own field, and
    its normalization. Compactness belongs here too, as a statement about that field.
  - **For an inertial factor:** the law relating momentum density to the material fields, which the accepted O2
    spec keeps OPEN (`directives/O2_SHARED_PHYSICS.md:133`, S8). The slab's density response is also OPEN there
    (`:138`, S8 and S22).
- **Inventories.** The record points at each source row without restating it:
  - **The plan's S16 section,** `V3_STEP_PLAN.md:670–746`. It covers the source statement and its transcription note
    (`:700–703`), the gravity-sector guardrail (`:705–730`) and the A11 gate (`:732–737`) (L1).
  - **The source's assumptions:** `../4d_1pn_full/paper/4d_1pn_full.tex:695–720`.
  - **Its width, defined for an isotropic Gaussian:** `:761–765`.
  - **Its reduced equation:** `:885–891`.
  - **The 2.5PN statement as written:** `../4d_2_5pn/paper/4d_2_5pn.tex:605–614` (L4).
  - **The charter's scope boundary:** `CHARTER.md:74–80` and `:109–119` (L5).
  - **The S22 debt rows** for `a_WT`, the `a_WT`/`a_mouth` bridge, and `M` with its multipoles:
    `V3_STEP_PLAN.md:1303–1305` (L5).
  - **The mass gap:** C3, `DEFECT_REGISTER.md:171` (L5).
  - **O2's S16 interface row:** `directives/O2_SHARED_PHYSICS.md:607`.
  - **O2's scoping row:** `directives/O2_steady_brane_balance_scoping.md:164` (L6).
  - **The three-mass guardrail:** `../../notes/muonium_gravity_corrections.md:21–25`.

## Retained OPEN
- **No owner assigned here:**
  - the forms of `G^i`, `Π^{ij}` and `f^i`;
  - the slab's material law. The accepted S12 record gives slab material completion to S8/S22
    (`steps/S12_order_source.md:292`) (L6).
- **No owner named:** the faces' law (same row).
- **Filed by this step as requirements:**
  - the field of the surroundings, its separation rule and its normalization;
  - the momentum–material law and the density response (O2, S8/S22).
- **S22:**
  - the core's interior, count and momentum;
  - the `a_WT`/`a_mouth` bridge, under A11.
- **`R-S14a-01`:** the slow limit and far behaviour.
- **S15, conditional on S14a:** the source side, meaning what force a drain exerts.
- **C3, a gap:** the link between the drain rate and mass.
- **The C10 tail question:** the plan's S19 (`V3_STEP_PLAN.md:766–788`).
- **The plan's deferred branch** (`:705–730`), not run here: species-dependent response.

## Intended claim
What the supplied balances give for a region around one throat that moves through its surroundings, with every term
live:
- the region's material and momentum balances, with the parts due to the surfaces' motion kept apart;
- the relations for the centre of the material count;
- the momentum the region receives;
- the material moments.

What a passive factor and an inertial factor would need is filed as requirements. The step gives:
- no mass;
- no inertial or passive factor;
- no field of the surroundings;
- no potential;
- no force law between bodies;
- no throat solve.

## Next
Pass 2: the spec content of this revision, reviewed until clear. Then:
- a build directive for Codex, derived from this list;
- G1 review of the build (fresh Claude + Grok), until clear;
- a cross-engine comparison;
- the record and its register pass, each reviewed until clear.
