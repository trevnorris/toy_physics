# O2 build r2: review dispositions (orchestrator)

**Artifacts:** the five O2 engine and harness files at `_scratch/s9b_build/o2_build_review_baseline_r2.sha256`,
the output of repair round 2. They are preserved here as r2, not accepted. Round 2 answered every r1 finding
routed to the engines (`O2_build_r1_review_disposition.md`).

**Legs (Codex-written → fresh Claude + Grok, identical prompt `o2_build_review_prompt_r2.md`).** This is build
review round 3. Both legs reported before any adjudication.
- **Fresh Claude (opus): NEEDS REVISION.** Two findings, both in the Wolfram engine
  (`_scratch/s9b_build/o2_build_review_r2_claude.md`). It found the SymPy engine clear.
- **Grok: NEEDS REVISION.** Three findings, all in the Wolfram engine
  (`_scratch/s9b_build/o2_build_review_r2_grok.txt`). It made no finding against the SymPy engine.

**Agreement.** Each leg constructed the objects independently and compared them with both engines. Every
explicit object matched with residual 0:
- the material velocity (`V^w = V·∇ξ_w`, tangent to the graph);
- `det g` and the graph normal;
- the coordinate-measure mass residual;
- the in-plane carried momentum `j_n V^i`;
- the graph-normal projection of the hold.

The carried bulk-direction component stays OPEN in both engines. Every K1–K13 site occurs exactly once, every
knife bites, and no knife is all-zero. Each leg's own FORM ablations bit in both engines, including:
- an untilted graph normal;
- the induced-measure mass law, where both engines give the same term `ρ V_r ξ′ξ″/(1+ξ′²)`;
- `∂_i(V^i ξ)` in place of `V·∇ξ`.

Neither leg found a defect in the SymPy engine in this round.

Each verification below is a mechanical lookup. The commands and their literal output are in
`O2_build_r2_review_disposition_lookups.md`. I ran no CAS (E1).

| # | Finding (leg) | Disposition | What must be true after repair |
|---|---|---|---|
| C1 | The WL bulk-normal load amplitude is `operand[TBulkNormalLive[s], …]`, a face-label constant multiplying a live normal. It has no native-point, time or history dependence, so no derivative of it can appear. SymPy's `action('T_bulk_n_s_live', native_point_context)` carries the face point. (Claude) | **ACCEPT, WL (M3).** Lookup: WL lines 57 and 104; PY line 229. Spec line 155 calls `𝒯_bulk,n,s^live` a "general signed OPEN normal-load amplitude", and lines 34–36 keep "loading and boundary data live, including their spatial derivatives and material-history dependence". The construction has been present since r0. | Every OPEN loading amplitude carries the live dependence the spec keeps: native point, time and history. No live operand is represented by a constant. |
| C2 = G2 | WL restricts native faces to single-valued height charts over the far-field `x` and prints the restriction. (Claude, Grok) | **ACCEPT, WL. Orchestrator error in the r1 disposition.** My r1 "what must be true" allowed a restriction to be "printed as a domain qualification". Spec lines 624–626 forbid that: "retain the named OPEN action where possible; otherwise report the missing restriction as a question for the user, without choosing it". The builder took the option I offered. | No restriction on the OPEN native face geometry beyond those the spec supplies. If an explicit action needs one, the OPEN action is kept. If it cannot be kept, the builder reports the restriction as a question and does not choose it. |
| G1 | WL emits `"Relation" -> Thread[holdBalance == ConstantArray[0,4]]` and `"Relation" -> (energyBalance == 0)`, asserting that the assembled momentum and energy balances vanish. SymPy emits the assembled sums with no relation. (Grok) | **ACCEPT, WL.** Lookup: WL lines 236 and 297; PY has no relation on the assembled balances. Spec lines 82–83: "supplies no assembled momentum or energy balance, expected residual, sign, cancellation"; line 273: "not a precomputed residual"; lines 621–622 defer "physical holder selection/solve". With `T_hold` OPEN, `== 0` presupposes a holder that closes the balance, which is a later step's question. The assertion has been present since r0, and both earlier rounds' legs missed it. | No emitted object asserts a value for the assembled momentum or energy balance. The constructed balance objects are emitted as constructed, oriented and unsolved. Relations the spec supplies, such as the steady mass balance, stay as supplied. |
| G3 | WL K10 replaces only `energyTransport` with 0. `JELive` (`energyCurrentOperand`) is still an operand of the energy density and the joint power occurrence. The directive's K10 removes `𝒥_E^live` "from the energy accounting". (Grok) | **ACCEPT, WL harness.** Lookup: directive line 62; WL harness lines 38–39; WL lines 61 and 265. The knife bites, but it is not the directive's knife, so the two engines' K10 triples do not mean the same thing. Present since r0. | K10 removes `𝒥_E^live` from every place it enters the energy accounting, and the harness output records what it removed. |

**SymPy engine.** No findings in round 3 from either leg. It is not accepted separately: acceptance is per
build, after the Wolfram repairs clear.

**Routing and the stop rule.** This was build review round 3. The r1 disposition committed the orchestrator to
stop and bring the build to the user if round 3 did not clear, rather than folding a fourth time.

None of the four findings was bred by a repair:
- C1, G1 and G3 have been present since r0 (provenance lookups); legs in two earlier rounds missed them, which is
  "legs wrong in both directions".
- C2 is the builder following my own r1 disposition.

So the author-change trigger (repairs breeding defects in the material just changed) is not met. All four
findings are in the Wolfram engine and its harness. **STOPPED: the decision goes to the user.**
