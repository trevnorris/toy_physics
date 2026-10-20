# S12: the throat's order-conversion source and its boundary data (pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-19. **Author of the build:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`.

**Gate.** This is a pre-builder decision list. It gets one Codex + Grok pass and is folded once (G2). It also names the
object the build computes, so that content is reviewed until clear (spec row). The builder starts only after that.

It was folded once on 2026-10-19, after a pass with six findings, all confirmed
(`directives/_measurements/S12_decision_list_review_p1_disposition.md`; reviewed version preserved at `cd505a87`).
The fold changes what the build computes, so the revised content gets its own review before the builder starts.

## Filing and framing
S12 is the plan's step "The throat's non-variational SOURCE, plus its boundary data" (`V3_STEP_PLAN.md:579`). The
user chose this framing on 2026-10-19: S12 states what any conversion law must satisfy and leaves the forms open. It
does not pick a drain law, and it does not define the brane: the user's requirements-first guide (2026-10-09)
collects what each sector needs and defines the brane last.

The flow is not fitted here. Each force states the flow it needs: light already has (`R-S21-01`), and gravity and EM
follow in their own steps. S21 compares them, and a mismatch reopens S9b (user, 2026-10-19). Galaxy-scale behaviour
is not a requirement at this step.

## Supplied equations (committed; unfalsifiable in this build)
The committed law, user decision of 2026-07-24 (`V3_STEP_PLAN.md:589–601`;
`../pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md:78–84`):

```
∂_t n + ∇₄·(n u) = 0                                   (total material)
∂_t(χ_B n) + ∇₄·(χ_B n u + J_χ) = n Γ_B ,   Γ_B = Γ_return − Γ_drain ,   [Γ_B] = T⁻¹
```

- `n` is the material number density on the 4D domain (in-plane `x^i`, bulk direction `w`), `u` its velocity,
  `χ_B` the ordered fraction, `J_χ` the order flux not carried by `u`. `Γ_B` is signed: positive is bulk material
  re-ordering onto the brane, negative is de-ordering (`../../docs/oriented_throat_return_source_correction.md:14`).
- The frozen-wall mass sink with a remote return stays ruled out. Any slow or static limit is derived from this law,
  never posed (`V3_STEP_PLAN.md:597–601`).

## User premises (2026-10-19)
Each is an **adopted substrate input to a conditional model**, not a derived result.

- **U1. Shear only inside the slab.** In the user's words: "the puncture must go past the slab thickness. This will
  be one key way of explaining why it loses the shear as it falls into the throat." Material holds shear only
  inside the slab's thickness. A throat, or a pinch hard enough, that carries material past that thickness
  de-orders it. Bulk material that settles into the slab re-orders. The user confirmed this one rule.
- **U2. Either direction, anywhere.** The exchange's sign is free at every point: "we shouldn't constrain that one
  to only one direction." De-ordering in throats and re-ordering elsewhere is the expected pattern, not a
  constraint.
- **U3. The bulk is an endless supply far away.** This is a statement about the bulk's IR condition. It does not hold
  the bulk still next to the slab.
- **U4. Re-ordering exceeds the drain overall.** The brane grows; the slab spreads instead of thickening. This is a
  sign premise on a global total, not a rate law.

**How U1 enters (an identification, flagged).** U1 is a condition on which order fields are allowed. It is not a
formula for the order.
- **The shear-holding region.** The slab has a region within which shear can be held at all. It is bounded by two
  live faces `w₋(x,t) < w₊(x,t)`, equivalently a live centre and thickness. The faces are properties of the slab, set
  by its own mechanics, which stay OPEN. They are not defined as wherever `χ_B` happens to be nonzero.
- **Outside it,** `χ_B` vanishes.
- **Inside it,** `χ_B` is a live field, between 0 and 1. It may depend on history or strain. No profile, rate or sign
  is chosen.
- **At a face,** `χ_B` may be continuous or may jump. Neither is chosen.

This reading is wrong if order can persist outside the shear-holding region, or if the region's extent is itself
defined by the order.

## What the build computes
Both engines, constructed independently, `.wl` blind:

1. **The source the committed law requires on order fields U1 allows.** The object is `n Γ_B`, with `n`, `u`, `J_χ`,
   the two faces and the interior order all live functions. It covers each part of the domain the condition
   distinguishes: outside the region, inside it, and at its faces, including any part concentrated on a face if the
   order jumps there. No flow profile is imposed, and the river profile does not enter.
2. **A labelled conditional specialization.** The same object where the interior order is a given function of
   position relative to the faces. This assumes ordering keeps pace with the slab. It is a specialization, not a limit
   earned from item 1, and it is labelled so wherever it is printed.
3. **Dimensional homogeneity** of both balances and of objects 1 and 2.

A name may name an object, never its value or a comparison outcome. The script prints objects and never states
conclusions. No acceptance criterion references an expected value.

## Inventories (record, sourced; not computed)
Kept apart, as the plan requires (`V3_STEP_PLAN.md:609–612`):
- **Source and controller functions.** `Γ_return`, `Γ_drain`, `J_χ`, and v2's return controllers. Whether the
  controllers have a role once the law conserves material locally is an open question for the record.
- **Boundary and domain data.** The mouth, the collar, the IR/return conditions and the `𝔅` branch selection
  (`V3_STEP_PLAN.md:612`). The IR/return scheme stays OPEN: one- or two-sided return, the return or absorber law, and
  far-field truncation (`../../docs/two_throat_simulation_handoff_spec.md:340`). U3 is a statement about the bulk's
  IR condition. No galactic or expansion setting is imposed.
- **The S9b/O2 domain, as recorded.** One isolated mass at rest with its drain flowing, Eulerian steady profiles,
  and lab time in the brane's far-field rest frame (`steps/O2_steady_brane_balance.md:74–75`). Its restrictions
  (isolation, steadiness, the far-field rest frame) are named as restrictions.
- **Momentum and energy partners of conversion.** O2 premise 3 closes the transported momentum only
  (`directives/O2_premise_decision_list.md:32–35`). OPEN with S12:
  - the additional momentum partners and the system that carries their reaction;
  - the bulk-direction carry, which is not `j_n` times one velocity (`directives/O2_SHARED_PHYSICS.md:243`);
  - the carried total energy, and the conversion-energy partners with their supply systems (`:280`).

  `R-S12-01` (the conditional reservoir) and `R-O2-01` (the net supplier) are separate obligations, and neither is
  assigned to S12 as a whole (`SUBSTRATE_REQUIREMENTS.md:521, 568`).
- **Native identifications (O6, `𝒥_map`).** The common exchanged material, and compatible mass, momentum, traction
  and work geometry and measures, between the 4D balance and the brane (`steps/O2_steady_brane_balance.md:399`;
  `directives/O2_SHARED_PHYSICS.md:605`). They are an OPEN S12 handoff, inventoried and not computed here.

## Model point
- **Regime:** none chosen. The objects are local, with every field live.
- **Freezes:** none intended. `J_χ` stays live (v2's `J_χ = 0` default is not adopted). Item 2 assumes fast ordering
  and is labelled so.
- **Transfer:** to wherever the committed law holds and the shear-holding region has two faces over the in-plane
  coordinates. The puncture's own interior is excluded and stays with S22.

## Retained OPEN
- the material branch, A13. The committed law is written for a real ordered fraction; a phase-carrying branch would
  add content not computed here;
- the interior order's law, including any dependence on history or strain; the faces' law; the form of `J_χ`;
- the projected order-loss source, `J`, the `J → Q` map and the return law (S14a; `V3_STEP_PLAN.md:636–637`);
- the link between drain strength and mass (S13);
- the flow profile (each force's requirement; S21);
- the slab's material law (S8/S22).

## Intended claim
The source the committed order law requires on order fields U1 allows, with every field live, and its labelled
fast-ordering specialization. Plus sourced inventories of the source functions, the boundary data, the conversion
partners and the native identifications. It is not a drain law, a throat solve or a flow.

## Next
A build directive for Codex derived from this list. Then G1 review of the build (fresh Claude + Grok, until clear),
the S12 record and its register pass, each reviewed until clear.
