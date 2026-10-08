# O2 comparator build r3: review dispositions (orchestrator)

**Artifacts:** the comparator, tests and new ablation harness from the fresh author. The user chose an author change
after round 3 (2026-10-07). The new author worked from the cumulative ten-item brief
`_scratch/s9b_build/o2_comparator_fresh_author_prompt.md`. The artifacts also include the builder's full guarded run in
`_scratch/s11c/o2-comparator-build-r3/` and the round-3 builder report. They are pinned at
`_scratch/s9b_build/o2_comparator_build_review_baseline_r3.sha256`, and the baseline check prints `OK` for all nine files
(lookups).

**Legs (Codex-written → fresh Claude + Grok, identical prompt `o2_comparator_build_review_prompt_r3.md`).** This is
comparator build review round 4. Both legs reported before any adjudication.
- **Fresh Claude (opus): NEEDS REVISION**, 3 findings (report in the session transcript; leg scratch
  `/tmp/o2cmp_r3_freshleg_1791415867/`).
- **Grok: NEEDS REVISION**, 1 finding (`_scratch/s9b_build/o2_comparator_build_review_r3_grok.txt`).

**Agreed sound by both:**
- Coverage: every leaf is on exactly one accounted path, with 0 uncovered and 0 multiply covered.
- Tables injective.
- Every exact residual recomputes to 0 from the raw transcripts. Claude's check found a maximum relative difference of
  4.7e−40 over 276 keys, including derivative evaluation points, the transpose layout and every closed balance part.
- Relations: both sides compared.
- One-sided corruptions caught, including the carried-entry sign flip.
- No target; exit 0.
- Claude's independent extraction of OPEN roles, named operands and live objects for `InternalForce_0` and
  `CarriedMomentumW` matches the printed deltas.

Each verification is a mechanical lookup. Commands and literal output are in
`O2_comparator_build_r3_review_disposition_lookups.md`.

| # | Finding (leg) | Disposition | What must be true after repair |
|---|---|---|---|
| V1 | Join `coupled_embedding` (L407) pairs objects at the wrong grain. SymPy `COUPLED_INPUTS/embedding` (audit L395–396) is a 3-tuple: `E_h_live`, the graph identity `ξ=ℓh`, and the label `normal_equation_identity_and_count_unsettled`. Wolfram `B_HOLD_LIVE/O4Identity` (L239) is `UnresolvedIdentification[OPEN[EhLive], holdBalance]`, and `holdBalance` is the full hold (L235). The count token `O4EquationIdentityCount -> Unsettled` (L340) is unjoined, with a reason (L1265) that names the SymPy label as its counterpart inside the joined tuple. The row prints `not_formed`, but its structure delta lists the 40 hold-balance roles as one-engine-only. (Grok 1) | **ACCEPT.** Lookups V1. The tuple's components have different counterparts: `EhLive` inside `O4Identity`, the field identity already joined at `field_identity`, and the count token. U1 requires that an object is unjoined only when the other stream has no counterpart for it. | Each component is paired with its same-role counterpart, or unjoined with a reason that holds. No row pairs a container with an object of a different role. No structure delta reports, as one-engine-only, content that the other engine emits on another row. |
| V2 | OPEN-occurrence inventories are presence sets (`features()`/`merge_fields()`, L984–1010: `names[…]=1`, `live[…]=1`). In production both engines pass the same live object into an OPEN action through several argument slots (velocity, metric, profiles, history). Claude froze `V_r` (6 occurrences) and `ξ_w′` inside one slot of the WL carried action. Every residual, field delta, balance entry and closed residual stayed unchanged; only a whole-tree hash moved, and that hash is non-empty at baseline. The control fixtures have the frozen object in one slot only. (Claude 1) | **ACCEPT as a finding about what may be claimed.** Lookups V2. An empty live-argument delta cannot be read as "same live arguments in every slot", and a one-sided freeze inside OPEN content is undetected when the object recurs. **How far to fix it is a scope question for the user** (below). | One of these, per the user's decision: (a) a live object frozen or stripped in any single declared slot of a paired OPEN occurrence moves a cross-engine product; or (b) the comparator computes and prints, per paired OPEN occurrence, which live objects recur across slots, so the record states exactly which single-slot freezes this instrument cannot see, and nothing is claimed beyond that. |
| V3 | SymPy `Str` labels inside OPEN actions are parsed as `Text`, and the named inventory counts only `Name` heads (L989). `OPEN_CarriedMomentumW`'s labels `outward_native_relative_mass_current` and `premise_3_local_material_velocity_at_each_transfer` are therefore invisible. The printed delta shows the Wolfram `NativeRelativeMassCurrent` and `Premise3LocalMaterialVelocity` as Wolfram-only. (Claude 2) | **ACCEPT.** Lookups V2. The delta misstates what the SymPy action names. | Every operand an OPEN occurrence names, labels included, enters the named inventory, bound through the name table where the spec names one object. A control mutates a label. |
| V4 | Two FORM ablations of the comparator survive all 54 tests: `sp.sqrt`→`sp.cbrt` in the WL `Sqrt` translation (L647–648), and the named inventory restricted to declared operands (L989). The first moves production residuals off 0, so it would read as an engine difference. The second makes the deltas look closer, so it would read as agreement. Grok independently deleted the `Sqrt` case and found the 54 tests still green. (Claude 3; Grok, same Sqrt gap, judged visible in output) | **ACCEPT.** Lookups V3. | A test fails under each of these ablations. |

## The pattern across four rounds, and the question it raises

Four review rounds and two authors have produced 14 accepted findings: round 0: 4, round 1: 4, round 2: 2, round 3: 4.
The closed results have not moved in any round. Every leg in every round independently recomputes the exact
cross-engine residuals, and all are 0. The findings fall into three groups:
1. **Join grain:** R1 and V1. Both pair OPEN containers whose representations differ.
2. **OPEN-content comparison (directive item 4):** S1, S4, T1, V2 and V3. Each round finds a new way the cross-engine
   comparison of OPEN placeholder internals is unfaithful: positional pairing, spelling-dependent keys, presence-set
   masking, invisible labels. These engines represent OPEN content differently by accepted design (seven named
   representational differences, `O2_build_r3_review_disposition.md`).
3. **Control coverage:** R3, S2, S3, T2 and V4. Each round finds further comparator-logic ablations that pass the
   tests. The space of possible ablations of the comparator's own logic has no fixed end, and directive item 6 names a
   finite list that every round has exceeded.

Under the standing rule (more than two rounds ⇒ question the premise with the user), this goes to the user before any
repair. Narrowing what is compared or claimed would drop a control, which is the user's decision alone (R11).

The reviewed baseline is preserved in a commit before any repair (G4).

## Outcome

**NOT ACCEPTED.** All four findings are accepted. V2's repair form depends on the user's scope decision.
