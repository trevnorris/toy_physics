# S15: the momentum the flow of two moving drains carries across a region around one of them (pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-24. **Author of the build:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`. The lookups behind this list are in
`directives/_measurements/S15_setup_lookups.md` ("L-n").

**Gate.** This list gets one Codex + Grok pass and is folded once (G2). It also names the objects the build computes,
so that content is reviewed until clear before any builder (spec row).

**Revision 2** folds pass 1 once (`directives/_measurements/S15_decision_list_review_p1_disposition.md`, D1–D5). The
reviewed revision 1 is preserved at `d6e4f74a`.

## Filing and framing
- **The plan's step.** S15, "Two drains attract" (`V3_STEP_PLAN.md:652–669`) (L1). The plan:
  - takes the drain strengths as "given", with a "Noether stress at infinity" (`:660`);
  - makes the step "CONDITIONAL ON S14a" (`:653`; `:639`);
  - says S16 "does NOT license S15", which "carries **its own** Noether / control-surface assumptions"
    (`:655–658`);
  - requires the source map from `J` to `Q` to be "an explicit step, not an unremarked substitution" (`:663–665`),
    with `{J, Θ_Q, N_∞,3, Q}` and the normalization in S22's debt table (`:1306–1307`);
  - limits the sign claim to "the **leading matter-stress** contribution only" (`:666–668`).
- **Where the user's decisions differ from the plan's statement.**
  - No drain strength is given. Each drain's inflow and density are general radial profiles.
  - The stress is the one P5 and P6 supply, below.
  - The drains move.
  - No force law, coefficient or normalization is printed, and no printed object is called the force on a throat.
- **This step's premises are its own** (`:657–658`). It uses S14a part 1's projected domain and takes no printed object
  from S16.
- **The framing carries over** (user, 2026-10-19). Each force states the flow it needs, and S21 compares them
  (`SUBSTRATE_REQUIREMENTS.md:759`, `:765`). This step states what its supplied balance gives for general flows. It
  picks no flow, drain law or material law. `R-S14a-01` stays OPEN, and this step takes no slow limit or far behaviour
  (L2).

## User decisions (2026-10-24)
The user selected three options. Their verbatim text is kept in the pass 1 disposition ("The user's decisions"),
and the record quotes it from there. What each choice puts in this step:
- **Stress law, "Use P5 + P6":** P5 and P6 are adopted as labelled premises (below). The build computes the momentum
  the two drains' combined flow carries into a region around one drain, piece by piece. The result is conditional on
  P5 and P6 and on `R-S14a-01`. No printed object is called the force on the throat.
- **Inflow, "Keep it general":** each drain's inflow and density are general radial profiles (below). No far flow is
  supplied.
- **Motion, "Let them move":** both drains' centres move, and their velocities stay live in every object.

Which profiles give which dependence on the separation is the record's reading of the printed objects. The build
prints the objects.

## Supplied equations (unfalsifiable in this build)
On S14a part 1's projected domain: in-plane coordinates `x^i` (three) and time `t`, with the slab integrated in `w`
between its live faces (L4).
- **The projected momentum balance, with P5 and P6.** S16 supplied the general form
  (`directives/S16_decision_list.md:76–86`). It is supplied here for this step, with P5's momentum density and
  carried current and P6's stress (`directives/S9b_SHARED_PHYSICS.md:340–378`, adopted in
  `steps/S9b_light_bending.md:667–668`) (L3, L4):

  ```
  ∂_t G^i + ∂_j Π^{ij} = f^i
  G^i    = ρ_br V^i
  Π^{ij} = ρ_br V^i V^j − T^{ij} + Π_o^{ij}
  T^{ij} = −p_br(ρ_br) δ^{ij} + η (∂^iV^j + ∂^jV^i − (2/3) δ^{ij} ∂_kV^k) + ζ δ^{ij} ∂_kV^k
  ```

  - `ρ_br` is the brane density and `V^i` the velocity of the brane material, which carries the momentum (P5). Both
    are general live fields of `x` and `t`.
  - `p_br` is a general function of `ρ_br`. `η` and `ζ` are general live fields (P6).
  - `Π_o^{ij}` is any other in-plane momentum current, a general live field. P5 keeps it as "a general OPEN action";
    it contains neither the carried current nor P6's stress.
  - `f^i` is the force density, a general live field. It includes every exchange of momentum with the faces and the
    bulk, including the momentum carried by material that crosses the faces. It stays undivided, as in S16.
- **Two flagged identifications.**
  - P5 and P6 are applied to the projected fields, as S9b applied them to the brane's in-plane fields.
  - P6 was adopted as a steady form (`S9b_SHARED_PHYSICS.md:357`, `:372`). Here the drains move, so the flow is not
    steady. P6 is applied as written. Any departure from it is not modelled and lies in `Π_o`.

## The drains, their fields and the region
- **Two drains.** Each drain `k = 1, 2` has a live centre `X_k(t)` and an excluded core around it. The cores' content
  stays with S22.
- **Each drain's fields** (the user's "Keep it general") are a radial inflow and a density profile about its live
  centre:

  ```
  V_k^i = u_k(r_k, t) (x^i − X_k^i(t)) / r_k ,    ρ_k = ρ_k(r_k, t) ,    r_k = |x − X_k(t)| ,
  ```

  with `u_k` and `ρ_k` general live functions. The radial form about the moving centre is a labelled choice.
- **The total fields** are `V = V_1 + V_2 + W` and `ρ_br = ρ_1 + ρ_2 + ϱ`. The remainders `W` and `ϱ` are general
  live fields, and they hold whatever the two profiles do not, including any background flow. This is a labelled
  split of the total fields, not a superposition law. Nothing makes the remainders small or zero.
- **The region** lies around drain 2, between an inner surface around its core and an outer surface, both moving
  with it. The build may take them as spheres about `X_2(t)` with live radii `r_in(t)` and `R(t)`. That is a choice
  of control surface, labelled where printed. Drain 1 lies outside the outer surface.

## What the build computes
Both engines construct these independently, with the `.wl` blind.

1. **The momentum crossing each surface into the region.** For the inner and the outer surface, two parts, kept
   apart, in the coordinates `x` of the supplied equations. On each surface `n` is the region's outward unit normal,
   so on the inner surface it points toward `X_2`. `v_S` is the velocity of a point of the surface.
   - **The rest part:** `−∮ Π^{ij} n_j dS`, in four pieces kept apart: the carried current, the pressure, the viscous
     stress and `Π_o`.
   - **The transport of `G` by the surface's motion:** `∮ G^i (v_S · n) dS`, in two pieces: the part of `v_S` from
     the motion of `X_2`, and the part from the change of the surface's radius.

   **Grouping by drain.** Each field of the split carries one label: drain 1 (`ρ_1`, `V_1`), drain 2 (`ρ_2`, `V_2`)
   or remainder (`ϱ`, `W`). Within each piece built from the split fields, every term is put in the group named by
   exactly the set of labels its factors carry, with one group for each nonempty set. Each term then lies in one
   group, and the groups add up to the piece. `η`, `ζ`, `n` and `v_S` carry no label. `Π_o` is not built from the
   split fields and is printed whole, as is the whole pressure (below).

   **The pressure.** The pressure of the total density is printed whole, where `p_br` is defined at `ρ_br`. Beside
   it, as a labelled choice, its expansion through first order in `ρ_1` about `ρ_2 + ϱ`, grouped as above (a
   function of `ρ_2 + ϱ` carries the labels drain 2 and remainder), with its order printed, and the whole pressure
   minus the expansion. The expansion exists only where `p_br` is defined and differentiable at `ρ_2 + ϱ`; that
   condition is labelled where it is printed.
2. **The region's momentum balance.** From the supplied balance: the rate of the region's momentum, the integral of
   `G^i`, beside item 1's parts at both surfaces and the integrated force density, each printed apart. The integrated
   force density stays undivided.
3. **Dimensional homogeneity** of items 1 and 2.

Every object is stated on its domain: `|X_1 − X_2| > R > r_in > 0`, with both cores excluded, and each pressure
object where its condition above holds.

The build computes no force law, coefficient, normalization or mass. It gives no relation of any printed object to
S14a's `J` or to the plan's `Q`, `Θ_Q` or `N_∞,3`. It identifies no printed object with the force on a throat.

A name may name an object, never its value or a comparison outcome. The script prints objects and never states
conclusions. No acceptance criterion references an expected value.

## Model point
- **Regime:** none chosen. Time dependence is kept, and both drains move. No steady, slow, incompressible, inviscid,
  point-drain or far-separation limit is taken: `c_comp`, `η`, `ζ`, both radii and the separation stay live.
- **Freezes:** none intended. The remainders stay live, so no background is at rest and no superposition is imposed.
  The force density and the other current stay general.
- **Labelled choices:**
  - each drain's fields as radial profiles about its moving centre;
  - the split of the total fields into those profiles and the remainders, and the grouping of terms by it;
  - the control surfaces' shape;
  - the pressure's expansion through first order in `ρ_1`, where it exists.
- **Adopted premises:** P5 and P6, labelled, with the two flagged identifications above.
- **Transfer:** wherever S14a part 1's projected domain and P5 and P6 hold, outside both cores, with drain 1 outside
  the outer surface.

## For the record and its register pass (not in the build directive)
- **The user's decisions.** The record quotes them from the pass 1 disposition ("The user's decisions").
- **What the record reads from the printed objects.** From the printed objects alone, it states:
  - which conditions on the profiles, the remainders, `Π_o` and `f` give a named printed piece, or a named sum of
    printed pieces, a stated dependence on the separation, including which give an inverse-square one, as the user's
    "Keep it general" asks. Each such dependence is stated with the conditions on the surfaces under which it holds:
    both radii, their rates, and the motion of `X_2`;
  - the terms that depend on each drain's velocity;
  - the sign of a printed piece only on its stated conditions, and only for that piece (`:666–668`).
- **Requirements.** The record files each requirement as the register's rest-on test admits it, or on a user decision
  that it quotes. The user's "Keep it general" files gravity's flow requirement for S21's comparison, as a condition
  on the named object, with its conditions on the surfaces. It is not called a pull or the force on a throat (Retained
  OPEN). It merges into `R-S21-01` where the register's merge test says the two cover the same requirement.
- **Carried from S16's acceptance.** The status lines that still present S16 as written for review are updated to its
  acceptance at `4f907e9d` (`directives/_measurements/S16_record_r2_review_disposition.md`, "Carried to the next step
  record").
- **Prior art, as an oracle (M3).** The record checks the printed objects against the gravity section's rows `G7`,
  `G9` and `G13` (`../../docs/medium_requirements_and_prior_art.md:417`, `:419`, `:423`), quoting both sides. It does
  so only where their conditions match. Where they differ, for example steady, incompressible or inviscid flow, or a
  point sink, it says so and lists the comparison as open (L6).
- **Inventories.** The record points at each source row without restating it:
  - the plan's S15 section, `V3_STEP_PLAN.md:652–669`, and its S22 rows, `:1305–1307` (L1);
  - the earlier construction of the force, as earlier material and not a premise (L7):
    `../../software/stage1_solver/reports/pathA_21c_force_from_noether_stress_tensor.md:44–53`,
    `../../software/stage1_solver/reports/pathA_21b_force_closure_and_profile_bvp.md:34`, and
    `../pde_ledger_v2/notes/stages/ledger_stage002_matter_stress_force_assembly.md:141–144`, `:201–202`;
  - O2's drain premise and historical matter stress, `directives/O2_SHARED_PHYSICS.md:94`, `:250–257`, and its
    transfer limit, `:351` (L5);
  - the charter's conditional, `CHARTER.md:238–241` (L8).

## Retained OPEN
- **No owner assigned here:** the forms of `f^i` and `Π_o^{ij}`. The slab's material law goes to S8/S22
  (`steps/S12_order_source.md:292`).
- **Whether any printed momentum acts on the throat.** That needs the momentum exchanged at the faces (`R-S16-03`, no
  owner named) and the core's content (S22). This step does not make the identification.
- **S22:** the source map `{J, Θ_Q, N_∞,3, Q}`, the normalization and `I_F,12` (`V3_STEP_PLAN.md:1306–1307`), and
  the cores' content.
- **The pieces of the full sign that P5 and P6 do not supply:** the plan's quantum, confinement and Maxwell residuals
  (`:666–668`). No owner named.
- **P6 outside steady flow,** and P5's possible stress-carried momentum (S8; `S9b_SHARED_PHYSICS.md:351–355`).
- **`R-S14a-01`:** the slow limit and far behaviour. S15 stays conditional on it (`V3_STEP_PLAN.md:639`).
- **S19:** what the velocity-dependent terms mean for velocity-dependent gravity (`V3_STEP_PLAN.md:766–788`).
- **C3, a gap:** the link between the drain rate and mass (`DEFECT_REGISTER.md:171`).

## Intended claim
What P5, P6 and the supplied balance give for the momentum crossing a region around one of two moving drains, with
general radial profiles and every term live, piece by piece. It is conditional on P5 and P6 and on `R-S14a-01`. The
step gives no force law, no normalization, no mass, no map from `J` to a drain strength, and no force on a throat.

## Next
Pass 1 is folded once, here. The spec content changed, so it gets a second Codex + Grok pass, and is reviewed until
clear. Then the build directive, the build and its G1 review (fresh Claude + Grok) until clear, a cross-engine
comparison, and the record and its register pass, each reviewed until clear.
