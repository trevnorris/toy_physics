# O2 comparator build r1: review dispositions (orchestrator)

**Artifacts:** the comparator and tests after repair round 1 (`a702782d` findings R1–R4), the builder's full
guarded run in `_scratch/s11c/o2-comparator-build-r1/`, and the round-1 builder report. They are pinned at
`_scratch/s9b_build/o2_comparator_build_review_baseline_r1.sha256`, and the baseline check prints `OK` for all
eight files (lookups). Repair round 1 stopped once on the builder's own fixture bug (`StopIteration` in its new
test). With my go-ahead it fixed that fixture without weakening the assertion and completed the round:
- 35 tests pass;
- each of 7 comparator-logic ablation copies trips a test;
- the full run exits 0.

**Legs (Codex-written → fresh Claude + Grok, identical prompt `o2_comparator_build_review_prompt_r1.md`).** This
is comparator build review round 2. Both legs reported before any adjudication.
- **Fresh Claude (opus): NEEDS REVISION**, 4 findings (`_scratch/s9b_build/o2_comparator_build_review_r1_claude.md`,
  condensed by me; leg scratch `/tmp/o2cmp_leg_claude_20261007_152021_3084964/`).
- **Grok: CLEAR**, no findings (`_scratch/s9b_build/o2_comparator_build_review_r1_grok.txt`).

**Agreed sound by both:**
- Coverage is complete: 159 joined and 202 unjoined, with no unaccounted leaves and no catch-all reasons.
- Both tables are injective, and the name citations hit their engine lines.
- Every exact residual each leg recomputed from the raw transcripts matches the printed one. Claude's numeric
  check ran at 3 generic points to 40 digits on 22 operands.
- One-sided corruptions are caught on the unablated comparator.
- There is no target and no verdict token.

On the current streams, everything printed that either leg could recompute is faithful. The R1–R4 repairs are in
place: `mass_loss_orientation` is re-paired, and the held-aggregate sign propagates (lookups S1, comparator
L679–709).

**Where the legs differ.** Grok's own FORM ablation, a Wolfram `Sqrt[arg]` collapsed to a bare `Sqrt`, left the 35
tests green, but it moved 91 of 136 production residuals off 0. Grok therefore judged the gap visible in the
output. Claude's FORM-A is a symmetric collapse of a derivative's evaluation point. It leaves the production output
unchanged and makes a one-sided corruption print 0: invisible both to the tests and to the output. Legs have been
wrong in both directions (G4), so each finding is adjudicated by evidence below.

Each verification is a mechanical lookup. Commands and literal output are in
`O2_comparator_build_r1_review_disposition_lookups.md`.

| # | Finding (leg) | Disposition | What must be true after repair |
|---|---|---|---|
| S1 | Only OPEN heads get an orientation role (comparator L703–706), and the hold rows are `not_formed`. The OPEN-free premise-3 carried in-plane entry `j_n V^i` in `ℬ_hold` (PY audit L259 `(1, carry)`; WL L232 `CarriedExchange,1`) has no orientation entry and no residual. A one-sided sign flip of it changes only node-head counts. (Claude 1) | **ACCEPT.** Lookups S1. Without this, a record can claim neither that entry's orientation nor the agreement of the balance's closed remainder. Both are physics the spec orients (§9). | Every additive entry of each assembled balance (hold and energy) has its engine orientation printed and compared, OPEN-free entries included. The closed (OPEN-free) part of each balance component is subtracted across the engines as its own residual. A one-sided flip of such an entry changes a printed comparison. |
| S2 | Collapsing a derivative's evaluation point (WL `Derivative[n][F][arg]`, comparator L597–606; PY `Subs`, L617–637), symmetrically in both engines, leaves the 35 tests green and every production residual unchanged. It also makes a one-sided corruption of a `det_g` evaluation point print `exact 0`. The stripped-argument control covers only a plain application (`VR[x1]→VR`, test L97–98). The derivative tests change the order only (L124–125). (Claude 2) | **ACCEPT.** Lookups S2. This is the `X(args)→X` collapse class, in its derivative form, and it can manufacture agreement. R3 required controls for argument content on the exact path; derivative evaluation points were not covered. | A one-sided move or strip of a derivative's evaluation point, in either engine's serialization, changes the printed residual. A test fails if the comparator's own derivative handling stops comparing evaluation points. |
| S3 | Dropping orientation through a held SymPy `Derivative` body (comparator L720–721) leaves the 35 tests green. In production it manufactures a sign difference (PY `OPEN_MomentumDensity_0` −1→+1 against WL −1) and hides a one-sided storage-sign corruption. No test has a negated `Derivative` body (lookups S3). (Claude 3) | **ACCEPT.** Lookups S3. R2/R3 required held-aggregate orientation controls; the `Derivative` body form is not covered. | A test fails if orientation stops propagating through any held linear form the engines use, the SymPy `Derivative` body included. A one-sided sign flip there changes the printed orientation comparison. |
| S4 | Item 4's differences are not computed. The structure census is printed per engine only. `differences()` (L651–665) is a positional tree diff that pairs children by index, so OPEN subtrees that do not correspond get paired. `live_arguments_and_binders` is the constant string `complete mapped operand tree; no argument erasure` (L739). (Claude 4) | **ACCEPT.** Lookups S4. The constant is a typed conclusion inside the instrument (E1). A positional pairing of non-corresponding subtrees is not a faithful difference, so a record would have to compare the OPEN content in prose. | For each OPEN action and balance entry, the differences item 4 lists are computed and printed through declared, cited tables: role or head, named OPEN operands, live arguments, orientation. Nothing pairs objects by position alone. No field states a conclusion that the comparator did not compute. |

## Outcome

**NOT ACCEPTED.** Four findings are accepted. S1 and S4 change what a record may claim. S2 and S3 are control gaps
of the kind R3 required closed, and S2 is a manufactured-agreement class. Three of the four (S1, S2, S4) were
present at r0 and are not bred by repair 1. S3 is in a held form the repair did not cover.

The reviewed baseline is preserved in a commit before any repair overwrites it (G4).

**This was the second review round, and it is not clear.** Another repair means a third review round. Under the
standing >2-rounds rule, the next step goes to the user before any repair: continue (same or fresh author), or
change the scope. Narrowing the scope would drop a control, which is the user's decision alone (R11).
