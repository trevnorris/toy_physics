# S14a part 2: steady throughput, the return and the far flux (pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-21. **Author of the build:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`. The setup lookups behind each pointer are in
`directives/_measurements/S14a_2_setup_lookups.md`.

**Gate.** This is a pre-builder decision list. It gets one Codex + Grok pass and is folded once (G2). It also names the
objects the build computes, so that content is reviewed until clear (spec row). The builder starts only after that.

**Revision 2.** Revision 2 is the one fold of the G2 pass (dispositions:
`directives/_measurements/S14a_2_decision_list_review_p1_disposition.md`; the pass-1 text is at `d72bb551`). It
defines the limit (Q1, Q2) and names the totals and budgets exactly (Q3).

## Filing and framing
- **The plan's step.** S14a, "The drain bridge" (`V3_STEP_PLAN.md:626`), asks for "the projected order-loss source ·
  the definition of `J` · the profile-dependent `J → Q` map · the controlled `ω → 0` return law" (`:636–637`).
  Until it succeeds, "S14, S15 and S20 are CONDITIONAL, not carried results" (`:639`).
- **Part 1 is accepted** (build `94fddf2f`, comparator `de81b9fc`, record `4a4dfaf5`). It prints the first two
  objects, with flux objects at a finite outer radius (`steps/S14a_1_projection.md:143–150`). It records "No
  steady/no-leakage, constant-density, radius-independent flux or far-radius limit is earned here" (`:140–141`).
  Its OPEN table leaves the return law to "A later part of S14a" and the `J → Q` map with no named later part
  (`:376–377`).
- **User decision, 2026-10-21: this part covers both remaining objects**, the `ω → 0` return law and the `J → Q` map.
- **The slow limit is derived, never posed.** The plan: "If a static limit is ever wanted it must be **derived** as
  the `ω→0` limit of the dynamical law, ⛔ never posed as an independent frozen solve" (`V3_STEP_PLAN.md:599–601`).
  The S12 decision list: "Any slow or static limit is derived from this law, never posed"
  (`directives/S12_decision_list.md:48–49`).
- **The framing carries over** (user, 2026-10-19). This part states what the committed law gives and picks no drain
  law and no return law's form. Each force states the flow it needs, and S21 compares them. Under "Keep S12" (user,
  2026-10-20), the conversion law's forms stay S12's and stay OPEN.

## Supplied equations (committed; unfalsifiable in this build)
The committed law, user decision of 2026-07-24 (`V3_STEP_PLAN.md:589–601`;
`../pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md:78–84`):

```
∂_t n + ∇₄·(n u) = 0                                   (total material)
∂_t(χ_B n) + ∇₄·(χ_B n u + J_χ) = n Γ_B ,   Γ_B = Γ_return − Γ_drain ,   [Γ_B] = T⁻¹
```
- **The fields** are part 1's, on the same 4D domain: in-plane `x^i` (three coordinates), bulk direction `w`, and
  time `t`. In this part `Γ_return` and `Γ_drain` are two separate live functions, as the committed law writes them;
  `Γ_B` is their difference. Stage 006's sign: "`Gamma_B < 0` = brane→bulk de-structuring (drain); `Gamma_B > 0` =
  bulk→brane re-ordering (return)" (`:86–87`). No sign is imposed on `Γ_B` or on either part.
- **U1 and its encoding** are as the S12 decision list states them (revision 4, cleared at `b08673a2`): "How U1
  enters", and item 1's "Fields at a face". The shear-holding region is bounded by two live faces
  `w₋(x,t) < w₊(x,t)`, and `χ_B` vanishes outside it. Both are adopted premises, not results of this build.
- **The material factor at a face**, as the accepted S12 record writes it (`steps/S12_order_source.md:79–80`), for
  region `r` at face `s`, using that region's one-sided trace at `w = w_s`:

  ```
  a_{r,s} = u_r^w − ∂_t w_s − u_r^i ∂_i w_s
  ```

  It is supplied and unfalsifiable in this build. The build's face terms are reached by its own integration.
- **Part 1's objects** are the accepted definitions this part takes the limit of: the projected balances (part 1
  item 1) and, in the spherical specialization, the face channel, the core channel, `J`, the order lost, storage and
  the outer flux (part 1 items 2–3; notation `steps/S14a_1_projection.md:80–127`). The SymPy engine may reuse the
  accepted part 1 SymPy construction (`scripts/S14a_1_projection.py`, accepted at `94fddf2f`). The `.wl` imports
  nothing and re-derives.

## User decisions on the limit and the far region (2026-10-21)
The user selected these two options. The wording is the orchestrator's option text.
- **The slow limit: "Steady throughput."** "The throat is at rest, and the drain and re-ordering run at constant
  rates. The pattern stops changing, except that the slab's steady spreading stays a live term. The flow runs
  throughout. No re-ordering form is picked; the forms stay S12's."
- **The far region: "Exchange live."** "Print how the far flux depends on radius, keeping the exchange and storage
  between the throat and the sphere. Take a far-radius limit only where the build shows one exists. Any condition it
  needs is printed, not imposed."

Not chosen: the slowly moving throat, with drive frequency `ω = k·V` (the G0 card's reading,
`../../software/em_charge_attribute/g0_closure_card_v0.md:396`); v2's exchange-free far region as a specialization;
and a split of the far flux into a throat share and a re-ordering share.

**The limit: the DC part (`ω → 0`).** This defines the limit for every item below. The plan calls the frozen
completions it replaces "DC sink/return completions" (`V3_STEP_PLAN.md:633`). No condition is placed on
`n Γ_return`, `n Γ_drain` or any other native field at a point.
- **What it acts on.** The projected objects of part 1 item 1, which are functions of in-plane position and time:
  the time rates of the projected material and order densities, the projected fluxes, the face terms, and the
  projected order source and its parts (below). In the specialization, also the shell objects of items 2 and 3, at
  fixed inner and outer radii. The integration in `w` between the live faces comes first, so no frame in `w` is
  chosen. The in-plane frame is the throat's: the throat is at rest, at the origin of the specialization.
- **The operation.** The DC part of such an object `X` is its zero-frequency part, the long-time average

  ```
  ⟨X⟩ = lim_{T→∞} (1/T) ∫_{t₀}^{t₀+T} X dt
  ```

  at fixed in-plane position, or at fixed radii. Whether it exists is printed as a condition on the live fields,
  never imposed. The DC part of a product or a quotient is the DC part of that whole object.
- **What the limit drops** is `X − ⟨X⟩`, for each object.
- **Steady throughput, as the user chose it.** The DC parts do not change in time, so the drain and re-ordering run
  at constant rates. The DC part of each projected density's time rate is the slab's steady spreading, and it stays
  live. No DC part is set to zero, and the slab may spread by thickening, by growing in area, or both.
- **Order of limits.** The average is taken at fixed position and fixed radii. Any far-radius limit (item 3) is
  taken after it, on the DC objects.

**The parts of the projected order source.** Part 1's projected order source is `Sχ = ∫n₂Γ_B,2 dw + C₋ + C₊`
(`steps/S14a_1_projection.md:89`), where `C₋` and `C₊` are its parts concentrated on the faces (`:87`). With
`Γ_B = Γ_return − Γ_drain`:
- the **return part** is `∫n₂Γ_return,2 dw` and the **drain part** is `∫n₂Γ_drain,2 dw`, each over the regular part
  inside the region only;
- the **face part** is `C₋ + C₊`, net. It is not split into return and drain: no such split is supplied.

Each concentrated part is counted once, in the face part.

## What the build computes
Both engines construct these independently, with the `.wl` blind.

1. **The DC part of the projected balances.** Apply the operation above to both time-dependent projected balances
   of part 1 item 1.
   - Print the DC part of each balance. The DC parts of the return, drain and face parts of the projected order
     source are printed separately.
   - **What the limit drops.** For each balance, print the time-dependent object, its DC part and what the limit
     drops, so nothing is dropped unseen.
   - Every field stays live.
2. **The return at DC.** In part 1's labelled specialization (spherical symmetry about the throat, the shell between
   an inner sphere around the excluded core and an outer sphere), print the DC parts of, separately:
   - **the re-ordering total, the drain total and the face total:** the return, drain and face parts integrated over
     the shell, `∫shell` as part 1 writes it (`steps/S14a_1_projection.md:112`);
   - **the order lost** over the shell, as part 1 defines it (`Lχ = −∫shell Sχ`, `:125`);
   - **the order budget over the shell:** part 1's projected order balance (`:91`) integrated over the shell. Its
     terms are the storage of projected order, `∂t∫shell Σχ`; the in-plane order flux `4πr²H_r` (`H` as at `:84`),
     outward through the outer sphere and inward through the inner sphere, as part 1 orients its outer flux and core
     channel (`:126`, `:123`); the face terms as that balance carries them, after the full concentrated part
     (`:88`); and `∫shell Sχ`, as the three totals;
   - **the material budget over the shell:** part 1's face channel, core channel, `J`, storage and outer flux, and
     the balance relating them.

   For each budget, print the time-dependent object, its DC part and what the limit drops.
3. **The far flux with exchange live** (user decision above). In the same specialization, print as functions of the
   outer radius:
   - the DC part of every object of item 2;
   - the DC part of the outer material flux, and, where the projected material density at the outer radius is
     nonzero over the average, the DC part of the volume flux of the projected velocity through the outer sphere.
     That condition is carried wherever the velocity or its volume flux is printed;
   - the dependence of each on the outer radius, with the exchange and storage between the inner and outer spheres
     kept.

   A far-radius limit is printed only where the build shows it exists, after the average. The condition on the live
   fields that it needs is printed as an object, never imposed.
4. **Dimensional homogeneity** of every new object in items 1–3.

A name may name an object, never its value or a comparison outcome. The script prints objects and never states
conclusions. No acceptance criterion references an expected value. Spherical symmetry is labelled as a
specialization wherever it is printed; it is not a limit earned from item 1. The DC part is labelled as the user's
steady-throughput limit wherever it is printed.

## Model point
- **Regime:** the DC part (`ω → 0`), the user's steady throughput (2026-10-21), reached from the time-dependent
  balances by the operation above. The time-dependent objects and what the limit drops are printed beside it.
- **Freezes:** none intended. The operation acts on whole projected objects and sets no DC part to zero. `J_χ`
  stays live, the faces keep their motion, and no exchange sign is imposed: U2 leaves the sign free. U3 and U4 enter no
  computed object, as in S12 (`steps/S12_order_source.md:154–155`). The user's premise that brane material leaves only
  in throats (2026-10-06) is not imposed.
- **Transfer:** wherever the committed law holds and the shear-holding region has two faces over the in-plane
  coordinates. The puncture's core is excluded and stays with S22.

## Inventories (for the record, under its own review; not computed, not used by the builder)
The record points at each source row, without restating it:
- **v2's `J → Q` map:**
  - `../../software/stage1_solver/reports/pathA_21_emergent_G_mass_bridge.md:19`, `:22` and `:62–66`: the rows for
    `S_leak`, `W_eff`, `N_∞,3`, `J` and `Θ_Q`;
  - `../pde_ledger_v2/notes/stages/ledger_stage004_gnls_action_dimensional_foundation.md:148–150`.
  The record compares what this build prints with these, by retrieval. Where the comparison needs computation, it is
  listed as open.
- **The picture the plan's item was written against:** the plan's text at `bcbc2716` (2026-07-31); the G0 card's
  return kernel, its stationary limit and its order of limits
  (`../../software/em_charge_attribute/g0_closure_card_v0.md:291–322` and `:387–398`). The G0 drain is ruled out
  (`V3_STEP_PLAN.md:597–601`); the card is a comparison only.
- **The steps that wait on S14a:** S14 (`V3_STEP_PLAN.md:643–646`), and S15's source-map row (`:663–665`). The record
  states which of S14a's four objects are now printed, and the standing of S14, S15 and S20 as the plan gives it.
- **Part 1's carried notes:** the stale status lines and the two notes in
  `directives/_measurements/S14a_1_record_r1_review_disposition.md:35–50`.

## Retained OPEN
- the link between `J` and mass (C3, a gap), and any link to orbital `GM`;
- the core's interior (S22);
- the flow profile (each force's requirement; S21);
- the slab's laws and S12's forms, with the owners the accepted S12 record gives (`steps/S12_order_source.md:292`).
  These include the forms of `Γ_return`, `Γ_drain` and `J_χ`, and the return controllers;
- the material branch, A13 (S1's gate, propagated through S5/S12);
- the readings the user did not choose (above). They are not computed, and no owner is named.

## Intended claim
What the committed law gives at DC, the user's steady-throughput limit, for the projected balances. For a throat at
rest with spherical symmetry: the DC re-ordering, drain and face totals over a shell, the order and material budgets,
`J`'s two channels and the outer flux. How each depends on the outer radius, with exchange live. What the limit
drops. A far-radius limit only where one is printed with its condition. It is not a mass, a drain law, a return form,
a throat solve or a flow.

## Next
A build directive for Codex derived from this list. Then:
- G1 review of the build (fresh Claude + Grok), until clear;
- a cross-engine comparison;
- the record and its register pass, each reviewed until clear.
