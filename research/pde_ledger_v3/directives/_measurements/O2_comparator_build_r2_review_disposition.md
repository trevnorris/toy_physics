# O2 comparator build r2: review dispositions (orchestrator)

**Artifacts:** the comparator and tests after repair round 2 (`e4c3adfd` findings S1–S4), the builder's full guarded
run in `_scratch/s11c/o2-comparator-build-r2/` and the round-2 builder report. They are pinned at
`_scratch/s9b_build/o2_comparator_build_review_baseline_r2.sha256`; the baseline check prints `OK` for all ten files
(lookups). The user chose the round-2 repair by the same builder (2026-10-07).

**Legs (Codex-written → fresh Claude + Grok, identical prompt `o2_comparator_build_review_prompt_r2.md`).** This is
comparator build review round 3. Both legs reported before any adjudication.
- **Fresh Claude (opus): NEEDS REVISION**, 2 findings (`_scratch/s9b_build/o2_comparator_build_review_r2_claude.md`,
  condensed by me; leg scratch `/tmp/o2cmp_review_r2_claude_20261007163402/`).
- **Grok: CLEAR**, no findings (`_scratch/s9b_build/o2_comparator_build_review_r2_grok.txt`).

**Agreed sound by both:**
- Coverage: 159 joined and 202 unjoined, with 0 unaccounted.
- Both tables are injective.
- Every exact residual is 0, and each leg recomputes the ones it checked. Claude's 40-digit evaluation agrees on 14
  objects, including the closed parts of `hold_inplane` and `hold_normal`.
- The readers are lossless.
- A name match is never a zero, and there is no target.
- One-sided corruptions on the unablated comparator are caught, including the carried-entry flip K2.
- The output reproduces byte for byte, and 46/46 tests pass.

**Where the legs overlap.** Grok's own FORM ablation, in which the relation branch is replaced by an exact 0, leaves
all 46 tests green. Grok noted that "the current tests would stay green if a later edit zeroed every relation", but
judged the production rows correct. This is the same gap as Claude's F2/FA1.

Each verification is a mechanical lookup. Commands and literal output are in
`O2_comparator_build_r2_review_disposition_lookups.md`.

| # | Finding (leg) | Disposition | What must be true after repair |
|---|---|---|---|
| T1 | The OPEN-action `live_arguments`, `binders` and `named_OPEN_operands` deltas (`features()`, L863–886) are built from the name-mapped but non-canonical tree: the key is `json.dumps(data(v))`. Same-object occurrences therefore differ by engine spelling: `FunctionName` vs `Symbol` heads, `Pow(…,1/2)` vs `Sqrt`, and occurrence counts (`V_r` 9 vs 7, `ξ_w′` 105 vs 12). The recursion visits arguments only, never heads, so a Wolfram `Derivative[n][F][arg]` is never counted, and the `'wl::Derivative'` clause (L875) cannot fire. SymPy `Subs` dummies count as named operands. The tests (L430–460) assert only that a mutated delta differs from the base delta. The base delta for identical objects is already non-empty, and freezing a live derivative moves only the `argument_trees` hash. (Claude F1; Grok observed the non-empty deltas but did not test them) | **ACCEPT.** Lookups T1. Nothing can be claimed about whether the two engines' OPEN actions take the same live quantities, and an M3 freeze inside an OPEN action is invisible. **This defect is in the material repaired in round 2 (S4).** | For the same live object on both sides, the live-argument, binder and named-operand deltas are empty. Freezing, re-ordering or re-pointing a live profile or its derivative on either side makes them non-empty. Keys are canonical in the objects compared: the same derivative of a profile at the same argument is one key in both engines, whatever the spelling. Binder placeholders are not operands, and the comparison does not rest on occurrence counts that depend on spelling. Tests assert the empty delta for identical objects and include per-side freeze and order controls through `run()`. |
| T2 | Two result-bearing paths have no control. FA1: a relation compared on its left operand only (L1051–1061) leaves the tests green and drops a real right-side corruption (K3) from the output. No test mentions a relation (lookups T2: 0 hits). FA3: the closed part kept to its first OPEN-free term leaves the tests green. `hold_normal` has 3 OPEN-free terms per engine. (Claude F2; Grok's relation ablation, same gap) | **ACCEPT.** Lookups T2. The production output is correct today, but the instrument's controls do not cover these paths. FA3 is in the closed-part comparison added in round 2 (S1). | One-sided corruption of either side of a relation, and of any OPEN-free term of a balance with several such terms, changes the printed comparison. A test fails if the comparator stops comparing relation right sides or drops a later closed term. |

## Outcome

**NOT ACCEPTED.** Both findings are accepted. T1 changes what a record may claim about OPEN content, and T2 is a
control gap on result-bearing paths.

**The repairs are breeding defects.** T1 is in the S4 material repaired in round 2, and FA3 is in the S1 material
repaired in round 2. Round 2 had already found a gap in a held form that round 1's repair did not cover (S3). Under
G4, the response is to change the author, not to fold again with the same one. This is also the third review round
that is not clear, so under the >2-rounds rule the next step goes to the user.

The reviewed baseline is preserved in a commit before any repair (G4).

Across all three rounds, the closed results have not moved. The exact residuals are 0, and every leg independently
recomputes them. The defects found are all in what the instrument can detect or claim about OPEN content and in its
controls.
