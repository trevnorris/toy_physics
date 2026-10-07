# S9b linked brane: spec-authoring directive (orchestrator)

**Author:** Claude (orchestrator), 2026-10-06. **Status:** folded once after one Codex + Grok pass
(`CLAUDE.md` G2; both NEEDS REVISION, on sources). The spec the author writes is then reviewed until clear by a
fresh Claude agent and Grok.

## Task
Revise `research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md` in place, as v9. The reviewed v8 is committed as
`c2f1cf2b`. Add one part and change nothing else in v8, except where the new part makes a v8 statement false.

## The object (user decision, 2026-10-06)
v8 treats the profiles `δ`, `V`, `ξ_w`, `ρ_br` and `j_n` as independent. In the model they are not necessarily
independent. `c_γ² = μ_⊥/ρ_br`, and near the mass the brane's own steady equations may tie its density, its
flow and its embedding together.

The new part asks for the Part B and Part C conditions, and the implied `j_n`, with every relation the model's
steady equations impose among those profiles included.

## What the author supplies
- **The relations, from the records.** Identify every steady-state relation the model's records write among
  `δ`, `V`, `ξ_w`, `ρ_br`, `μ_⊥` and `j_n` near the mass, **with `V` and `j_n` live**. Supply each one as an
  equation, labelled supplied, with its source lines. Sources to check:
  - `research/pde_ledger_v3/`: `steps/S11b_interface_coupling_law.md`, `directives/S11b_SHARED_PHYSICS.md`,
    `directives/S11c_a_SHARED_PHYSICS.md` §2d, `directives/S11c_b_SHARED_PHYSICS.md` §2b, the S11c step records,
    and `V3_STEP_PLAN.md` Q1–Q2;
  - the embedding sector, read-only: `research/pde_ledger_v2/notes/stages/ledger_stage030_electric_scalar_localized_h_closure.md`
    and `ledger_stage031_puncture_deflection_field_identity_source.md`. Verify field identity and applicability,
    and keep their postulated and static qualifications. Supply governing equations, not solved profiles.

  Three records are not the relation asked for here:
  - `S11b_SHARED_PHYSICS.md:353` is the linear momentum balance for the wave displacement `u`, not a balance for
    the background. S11b keeps `v_dr` out of every operator (`:103–111`).
  - S11c-a `:261` and S11c-b `:196` balance the background against an external support with
    `V_s⁰ = J_s⁰ = 𝒜_s⁰ = 0`, which freezes the flow.
  - The embedding holder is an unsettled debt (`V3_STEP_PLAN.md:896`).

  Do not adopt any of the three as a steady relation with live flow. List the support balance and the holder as
  OPEN.
- **OPEN premises.** A premise the records do not supply is listed as OPEN, with what is missing. Examples are
  how `μ_⊥` depends on density, the force that drives the flow, and the momentum carried by `j_n`. An OPEN
  premise enters the deliverables as a live symbol or response, like Part C's `s`. ⛔ Do not invent a law.
- **The order counting.** Keep v8's `ε` counting, and state where each new relation enters it.

## Constraints
- ⛔ No expected value, sign, cancellation or mechanism stated as an outcome (`CLAUDE.md` M2).
- Keep every varying quantity live (M3). A relation that holds only with a quantity frozen is reported, not
  adopted.
- Physics only. Implementation belongs to the build.

## Output
- **The v9 spec**, edited in place.
- **A note** at `research/pde_ledger_v3/directives/S9b_linked_brane_sources.md`, with each supplied relation, its
  source lines, and each OPEN premise.

Commit nothing. Stop after writing both.
