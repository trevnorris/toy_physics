# S12: the throat's order-conversion source and its boundary data (pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-19. **Author of the build:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`.

**Gate.** This is a pre-builder decision list. It gets one Codex + Grok pass and is folded once (G2). It also names the
object the build computes, so that content is reviewed until clear (spec row). The builder starts only after that.

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
- **U3. The bulk is an endless supply far away.** This is a condition at the outer edge. It does not hold the bulk
  still next to the slab.
- **U4. Re-ordering exceeds the drain overall.** The brane grows; the slab spreads instead of thickening. This is a
  sign premise on a global total, not a rate law.

**Identification, flagged.** U1 is encoded as follows: the order a point *would* have is set by where it sits
relative to the slab, through a live function `F(w − ζ, W)` of its offset from the slab's centre `ζ(x,t)` and the
slab's thickness `W(x,t)`. The actual `χ_B` is `F` plus a live deviation from it. How fast `χ_B` approaches `F` is
**not** chosen. The case with no deviation, where order follows the slab instantly, is reached as a limit of the
general case. It is not posed. This encoding is wrong if order depends on something other than position relative to
the slab, such as the material's history or its strain.

## What the build computes
Both engines, constructed independently, `.wl` blind:

1. **The source the committed law requires under U1.** The object is `n Γ_B` on the U1 order field, with `n`, `u`,
   `J_χ`, `ζ`, `W`, `F` and the deviation all live functions. No flow profile is imposed, and the river profile does
   not enter. The same object with the deviation set to zero follows as the limit case.
2. **Dimensional homogeneity** of both balances and of object 1.

A name may name an object, never its value or a comparison outcome. The script prints objects and never states
conclusions. No acceptance criterion references an expected value.

## Inventories (record, sourced; not computed)
Kept apart, as the plan requires (`V3_STEP_PLAN.md:609–612`):
- **Source and controller functions.** `Γ_return`, `Γ_drain`, `J_χ`, and v2's return controllers. Whether the
  controllers have a role once the law conserves material locally is an open question for the record.
- **Boundary and domain data.** The mouth, the collar, the outer edge and the `𝔅` branch selection. The outer
  edge stays open: a body sits inside a galaxy's flow inside an expanding brane. "One body alone in a resting
  brane", the S9b/O2 domain, is a recorded freeze.
- **Momentum and energy partners of conversion.** Converted material carries the local brane velocity `V` (O2
  premise 3, `directives/O2_premise_decision_list.md`). The energy of ordering, and who supplies it, stay OPEN (`R-O2-01`, `R-S12-01`).

## Model point
- **Regime:** none chosen. The objects are local, with every field live, at any timescale.
- **Freezes:** none intended. `J_χ` stays live (v2's `J_χ = 0` default is not adopted).
- **Transfer:** to wherever the committed law holds and the slab has a centre and thickness over the in-plane
  coordinates. The puncture's own interior is excluded and stays with S22.

## Retained OPEN
- the material branch, A13. The committed law is written for a real ordered fraction; a phase-carrying branch would
  add content not computed here;
- the rate law for `χ_B` approaching `F`, the form of `F`, and the form of `J_χ`;
- the projection onto the slab and the far-field flux bridge, which tie `Γ_B` to the brane's exchange `j_n`
  (S14a);
- the link between drain strength and mass (S13);
- the flow profile (each force's requirement; S21);
- the slab's material law (S8/S22).

## Intended claim
The source the committed order law requires under U1, with every field live, plus sourced inventories of the source
functions, the boundary data and the conversion partners. It is not a drain law, a throat solve or a flow.

## Next
A build directive for Codex derived from this list. Then G1 review of the build (fresh Claude + Grok, until clear),
the S12 record and its register pass, each reviewed until clear.
