# S15 decision list, review pass 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S15_decision_list_review_p1_lookups.md`, beside this file
(generator `_scratch/s15/gen/s15_dl_p1_lookups.sh`). "P-n" cites a section of that file. Every computed result
quoted here is a leg's saved stdout; the legs' scripts were not re-run.

## The pass

- **Artifact.** `directives/S15_decision_list.md` (sha256 `6cc19b66…`) and its setup lookups (`0076cc5d…`),
  uncommitted on `4f907e9d`. Both matched the baseline after the pass. The working tree was unchanged apart from the
  lookups file (P0).
- **Role.** A pre-builder decision list: one Codex + Grok pass, folded once (G2). The list also names the objects
  the build computes, so that content is reviewed until clear (spec row).
- **Prompt and legs.** `_scratch/s15/dl/dl_review_prompt_p1.md`, sha256 `d348351d…` (P0). Both legs received it
  identically:
  - Codex (`gpt-6.1-sol`, xhigh): `/var/projects/_review_scratch/s15-dl-p1-codex/REVIEW.md`, exit 0;
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s15-dl-p1-grok/REVIEW.md`, exit 0.
- **Verdicts.** Codex: FINDINGS, 3. Grok: FINDINGS, 2. No two findings concern the same object.

## Findings

**D1. The list states a prior-art result and the options' reasons.** (Codex 1, HIGH.) **CONFIRMED.**
- **The list** (`:32–33`, `:41–42`): "The build directive carries the choices, not the options' reasons (M2)"; and
  "Keep each drain's velocity live. A moving drain that swallows material picks up a drag, the classic objection to
  'gravity as inflow' (prior-art row G13). Holding the drains still would hide it." (P1).
- **The sources.**
  - Row G13: "A moving drain feels a drag: Wang's formula (32) gives *"a universal damping force"*."
    (`docs/medium_requirements_and_prior_art.md:423`) (P1).
  - M2: "The spec says what to **compute** — ⛔ never what anything equals, is expected, or was measured"; "Prior
    art's result for our object is never a supplied premise" (`CLAUDE.md:81–82`, `:85–86`) (P1).
- **Effect.** The list is the spec the build directive is made from. A promise about the directive does not take a
  stated result out of the spec.
- **The other two option texts** carry reasons too: "S14a left the far flow open, so supplying it would freeze it",
  and "Keep it general" names the dependence the record looks for (P1). No leg reported these. They are the same kind
  of text, so the fold treats all three alike. This extension is the orchestrator's.
- **What must be true.** The list states each choice as what is computed or kept live, with no option's reason and
  no prior-art result. The option texts, verbatim, are kept below ("The user's decisions"), and the list points
  there. The record section may keep its reading instruction, marked as not for the build directive.

**D2. The pressure's expansion need not exist on the stated domain.** (Codex 2, HIGH.) **CONFIRMED.**
- **The list** (`:85–87`, `:104–105`, `:111`): "Nothing makes the remainders small or zero"; the pressure is "printed
  whole, and beside it, as a labelled choice, its expansion through first order in `ρ_1` about `ρ_2 + ϱ`"; the domain
  is "`|X_1 − X_2| > R > r_in > 0`, with both cores excluded" (P2).
- **The source.** "`p_br` is a general function of `ρ_br` only" (`directives/S9b_SHARED_PHYSICS.md:367–370`). The
  list states no domain for `p_br` (`EOS`: 0 matches) (P2).
- **The evidence.** Codex: with `p = log ρ`, a probe where `ρ_1 = 27`, `ρ_2 = 1`, `ϱ = −1` has total density `27`
  and `p(total) = log(27)`, while the expansion's base density is `0`; along the scaled path the pressure goes to
  `−oo` and its first derivative to `oo` (P2).
- **Effect.** For an allowed EOS and an allowed split, the whole pressure is defined and the expansion is not.
- **What must be true.** The expansion is printed with the condition under which it exists (`p_br` defined and
  differentiable at `ρ_2 + ϱ`), and with what it leaves out (the whole pressure minus the expansion), so the
  truncation stays visible. The whole pressure carries only its own condition (`p_br` defined at the total density).

**D3. The record's reading has no named object, and a control surface's motion alone can change it.** (Codex 3,
MEDIUM.) **CONFIRMED.**
- **The list** (`:99–100`, `:134–141`): transport "in two pieces: from the motion of `X_2`, and from the change of
  the surface's radius"; the record states "which conditions on the profiles, the remainders, `Π_o` and `f` give the
  crossing momentum a stated dependence on the separation, including which give an inverse-square one"; "The user's
  'Keep it general' files gravity's flow requirement" (P3).
- **The source.** The accepted S16 record: "`momentumReceived` is a three-vector made from the integrated `f_i` and
  the two **rest** momentum-flux surface pieces"; motion pieces print separately; "Momentum received is not a passive
  factor or the force a drain exerts on another body" (`steps/S16_moving_region.md:152–160`) (P3).
- **The evidence.** Codex: one set of fields at one instant, with the local balance satisfied
  (`local_balance_lhs_minus_f = 0`). The rest part plus the integrated force density is the same in both cases. The
  full crossing is `4π/d²` when `Rdot = 1/d²` and `4π` when `Rdot = 1` (P3). The legs' configurations are test
  cases, not S15 premises.
- **Effect.** "The crossing momentum" can gain or lose an inverse-square dependence through the outer surface's
  motion alone. A requirement filed on it would not be a requirement on the flow.
- **What must be true.** The record states each dependence on the separation for a named printed piece, or a named
  sum of printed pieces, with the conditions on the surfaces (both radii, their rates, and the motion of `X_2`) under
  which it holds. Gravity's flow requirement is filed as a condition on that named object, with those conditions. It
  is not called a pull or the force on a throat; that identification stays OPEN (`R-S16-03`, S22).

**D4. The surface normal has no orientation.** (Grok 1.) **CONFIRMED.**
- **The list** (`:97`): "`−Π^{ij}` contracted with the surface normal and integrated over the surface". The
  selected option says the flow "carries into a region" (P4).
- **The source.** P6's traction `t_br,inplane^i ≡ T^{ij} n_j` (`directives/S9b_SHARED_PHYSICS.md:366`) does not
  orient the region's surfaces. The list has no "outward" (0 matches) (P4).
- **The evidence.** Grok: for a radial flux `h(r) ê_r`, the outward-normal reading reproduces the volume integral of
  the divergence (`outward_minus_volume = 0`); taking `ê_r` on both spheres differs from it by `8π r_in² h(r_in)`. The
  inner surface's transport factors change sign between the two readings (P4).
- **What must be true.** One orientation is stated for every surface object: the region's outward unit normal, which
  points toward `X_2` on the inner surface. The transport pieces are stated as equations with that normal.

**D5. The grouping by drain is not a partition.** (Grok 2.) **CONFIRMED.**
- **The list** (`:102–103`): "the terms that contain fields of both drains are printed apart from the terms that
  contain one drain's fields only, and every term that contains a remainder is printed apart" (P5).
- **The evidence.** Grok: in the carried current, three terms carry both drains and a remainder
  (`2*V1*V2*rem_rho, 2*V1*W*rho2, 2*V2*W*rho1`). Two partitions that each put every term in one bin both add up to
  the product, and their both-drains bins differ by those three terms. Printing each term in every bin it fits
  over-counts by the same three terms (P5).
- **What must be true.** The grouping is a partition. Each term lies in exactly one printed group, the groups add up
  to the piece, and the rule for a term that carries both drains and a remainder is stated. A piece not built from
  the split fields is printed whole.

## What the legs checked and found sound
These do not change the objects. Grok's report has the scripts and stdout (`REVIEW.md`, "Checks that do not change
the objects"):
- P5 and P6 are quoted as adopted; the assembled balance puts P6's force on the force side;
- the three decisions are in the objects, and a uniform stream stays in the remainder;
- the split is not unique, but the whole pressure and the totals do not depend on it, and the list labels it;
- no labelled choice freezes a live quantity;
- the ownership statements match their sources, and both centres lie outside the region on the stated domain.

## The user's decisions (2026-10-24), verbatim
The user selected these options in conversation. The wording is the orchestrator's option text, as revision 1 of
the list carried it (`:34–42`). The record quotes them from here. The build directive carries the choices only.
- **Stress law, "Use P5 + P6":** "S15 adopts P5 and P6 as labelled premises. It computes the momentum the two drains'
  combined flow carries into a region around one drain, piece by piece (carried flow momentum, pressure, viscosity).
  The result is conditional on P5/P6 and on S14a's open far-flow requirement. It is not yet called the force on the
  throat."
- **Inflow, "Keep it general":** "Any radial inflow and density profile around each drain. S15 then says which
  profiles give an inverse-square pull, and that becomes gravity's flow requirement, compared with light's in S21.
  S14a left the far flow open, so supplying it would freeze it."
- **Motion, "Let them move":** "Keep each drain's velocity live. A moving drain that swallows material picks up a
  drag, the classic objection to 'gravity as inflow' (prior-art row G13). Holding the drains still would hide it."

## The fold
Revision 2 makes D1–D5 true, once (G2). It changes what is computed (D2, D4, D5) and what may be claimed (D3), so
its spec content gets another two-leg pass, reviewed until clear.
