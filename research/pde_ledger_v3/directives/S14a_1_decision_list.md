# S14a part 1: the drain rate and its far-field flux (pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-21. **Author of the build:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`.

**Gate.** This is a pre-builder decision list. It gets one Codex + Grok pass and is folded once (G2). It also names the
object the build computes, so that content is reviewed until clear (spec row). The builder starts only after that.

## Filing and framing
- **The plan's steps.** This covers the plan's S13, "The drain rate `J`, and what it is not" (`V3_STEP_PLAN.md:619`),
  and the first part of S14a, "The drain bridge" (`V3_STEP_PLAN.md:626`). S14a asks for "the projected order-loss
  source · the definition of `J` · the profile-dependent `J → Q` map · the controlled `ω → 0` return law"
  (`:636–637`).
- **User decision, 2026-10-21: S13 is folded into S14a.** S13 has nothing of its own to compute: its definition of
  `J` is S14a's, and its mass link is the known gap C3 (`DEFECT_REGISTER.md:171`). This part derives `J` and its
  far-field flux. The `ω → 0` return law is a later part of S14a.
- **The framing carries over** (user, 2026-10-19). This part states what the committed law gives and picks no drain
  law. Each force states the flow it needs, and S21 compares them. Under "Keep S12" (user, 2026-10-20), the
  conversion law's forms stay S12's and stay OPEN.

## Supplied equations (committed; unfalsifiable in this build)
The committed law, user decision of 2026-07-24 (`V3_STEP_PLAN.md:589–601`;
`../pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md:78–84`):

```
∂_t n + ∇₄·(n u) = 0                                   (total material)
∂_t(χ_B n) + ∇₄·(χ_B n u + J_χ) = n Γ_B ,   Γ_B = Γ_return − Γ_drain ,   [Γ_B] = T⁻¹
```
- **The fields.** `n`, `u`, `χ_B`, `J_χ` and `Γ_B` are the objects S12 used, on the same 4D domain: in-plane `x^i`
  (three coordinates), bulk direction `w`, and time `t`.
- **U1 and its encoding** are as the S12 decision list states them (`directives/S12_decision_list.md`, revision 4,
  cleared at `b08673a2`): "How U1 enters", and item 1's "Fields at a face". The shear-holding region is bounded by
  two live faces `w₋(x,t) < w₊(x,t)`, and `χ_B` vanishes outside it. Both are adopted premises, not results of this
  build.
- **v2's labels**, as v2 defined them, for comparison in the record only. The build does not use them:
  - `J`: "Number-rate flux lim_{S_R} int n_3 v_brane.n dS in a no-leakage steady region"
    (`../../software/stage1_solver/reports/pathA_21_emergent_G_mass_bridge.md:65`);
  - `Q`, `Θ_Q`: the same report, `:19` and `:66`.

## User decision on what `J` counts (2026-10-21)
`J` counts **both channels** by which material leaves the slab around a throat:
- through its faces nearby;
- by flowing in-plane into the puncture's core.

The core is outside this domain, and its fate stays with S22. Each channel is printed separately, along with the
order lost over the same region, so nothing is hidden in the sum.

## What the build computes
Both engines construct these independently, with the `.wl` blind.

1. **The projected balances.** Integrate both committed balances in `w` over the closed shear-holding region, with
   each face's concentrated part counted in full.
   - The result is the projected total-material balance and the projected order balance, at each in-plane point and
     time.
   - Every field stays live, as S12 kept them.
   - The build prints, as objects:
     - the projected densities and their in-plane fluxes;
     - every term the integration leaves at each face;
     - the projected order source, including its face parts.
2. **The drain rate and its channels.** Work in a labelled specialization: spherical symmetry about the throat in
   the three in-plane coordinates, with every profile a live function of radius and time. The region is the shell
   between an inner sphere around the excluded core and an outer sphere. The build prints, separately:
   - **the face channel:** the net rate at which material leaves the shear-holding region through its faces over the
     shell;
   - **the core channel:** the net rate at which projected material crosses the inner sphere inward, in-plane;
   - **`J`**, the sum of the two channels (user decision above);
   - **the order lost** over the shell, including face parts;
   - **the outward in-plane flux of projected material** through the outer sphere;
   - **storage:** the rate of change of the projected material in the shell.

   It also prints the balance that relates them, and the difference between the face channel and the order lost.
   Each is a computed object, not a typed one.
3. **The far-field flux.** In the same specialization, the build prints how the outer-sphere flux depends on the
   outer radius, and the volume flux of the projected velocity through the outer sphere. The projected velocity is
   the projected material flux divided by the projected material density. Profiles stay live, and no limit in
   radius is taken.
4. **Dimensional homogeneity** of the projected balances and of objects 2 and 3.

A name may name an object, never its value or a comparison outcome. The script prints objects and never states
conclusions. No acceptance criterion references an expected value. Spherical symmetry is labelled as a
specialization wherever it is printed. It is not a limit earned from item 1.

## Model point
- **Regime:** none chosen. Time dependence is kept, with no steady assumption.
- **Freezes:** none intended. `J_χ` stays live, the faces move, and no exchange sign is imposed. U2 leaves the sign
  free. The user's premise that brane material leaves only in throats (2026-10-06) is not imposed here.
- **Transfer:** wherever the committed law holds and the shear-holding region has two faces over the in-plane
  coordinates. The puncture's core is excluded and stays with S22.

## Inventories (for the record, under its own review; not computed, not used by the builder)
The record points at each source row, without restating it:
- **the mass gap:** C3 (`DEFECT_REGISTER.md:171`), with its v2 sources:
  - `../../software/stage1_solver/reports/pathA_21_emergent_G_mass_bridge.md:97–104`;
  - `../pde_ledger_v2/notes/stages/ledger_stage005_sound_speed_light_ratio.md:216–220`;
  - `../pde_ledger_v2/notes/parameter_register.md:152–153`.
- **v2's `J`, `Q` and `Θ_Q`:**
  - the report rows above;
  - `../pde_ledger_v2/notes/stages/ledger_stage004_gnls_action_dimensional_foundation.md:142–150`.
- **the material-conservation conflict:** `../pde_ledger_v2/notes/stages/ledger_stage002_matter_stress_force_assembly.md:53`
  and `:111`, against `../pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md:81` and `:87`, and
  the plan, `V3_STEP_PLAN.md:630–634`.
- **other objects named `J`:**
  - `../pde_ledger_v2/notes/parameter_register.md:147`;
  - `directives/O2_SHARED_PHYSICS.md:483–488`;
  - `../pde_ledger_v2/notes/parameter_register.md:235`.
- **the plan's rows on S13, S14a and S15:** `V3_STEP_PLAN.md:619–641` and `:663–665`.
- **an earlier v3 projected material balance**, in the O2 spec (`directives/O2_SHARED_PHYSICS.md:485–488`). The record
  compares it with what this build prints. It is not supplied to the build.

## Retained OPEN
- the `ω → 0` return law (a later part of S14a);
- the link between `J` and mass (C3, a gap), and any link to orbital `GM`. S9b defines `GM` independently of any
  drain amplitude (`steps/S9b_light_bending.md:79–80`);
- the core's interior (S22);
- the flow profile (each force's requirement; S21);
- the conversion law's forms, the faces' law and the form of `J_χ` (S12, under "Keep S12");
- the material branch, A13 (S1's gate, propagated through S5/S12). The committed law is written for a real ordered
  fraction.

## Intended claim
What the committed law gives for the slab-projected balances. For a throat with spherical symmetry: its drain rate,
the two channels that make up that rate, the order lost over the same shell, and its flux at any outer radius, all
with every profile live. It is not a mass, a drain law, a throat solve or a flow.

## Next
A build directive for Codex derived from this list. Then:
- G1 review of the build (fresh Claude + Grok), until clear;
- a cross-engine comparison;
- the record and its register pass, each reviewed until clear.
