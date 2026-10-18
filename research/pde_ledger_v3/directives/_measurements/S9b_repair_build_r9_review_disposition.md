# S9b repair build, review round 9 (scoped): dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S9b_repair_build_r9_review_lookups.md` beside this file
(generator `_scratch/s9b_build/gen/s9b_build_r9_lookups.sh`).

## The round
- **Artifacts.** Both ablation harnesses, the Wolfram checks script, the SymPy tag inventory and run measurements,
  the Wolfram round-9 measurement files and its rewritten round-8 measurements file, and both builders' reports.
  33 files are in the baseline `_scratch/s9b_build/s9b_repair_build_review_baseline_r9.sha256`, all OK at review.
- **Contract.** The build directive at `36ccec3e` (amendments 1–12), item 14 as amended by amendments 11 and 12. The
  spec is at `34a994d5`.
- **Builders.** The same `gpt-6.1-sol` session per engine as in round 8, resumed. **Item-15 deviation, recorded
  again:** item 15 names `gpt-6-astra` at high effort, in a fresh session per engine.
- **Run history.**
  - SymPy: run `full1` finished with exit 0 after 8.3 h. The builder then found a reader defect in its own harness:
    exported copies of shared objects were read as whole tuples. It fixed it, passed its software checks, and ran
    `full2` (exit 0, about 9.6 h). Only `full2` and `checks3` were handed to the legs.
  - Wolfram: it ran out of memory at 4 GiB in K9b (an item 11 stop), and the guard then refused launches at 9 and
    5 GiB while the SymPy job held 22 GiB. The final runs at 13 GiB, on a machine with no other job, all exited 0.
    Only those were handed.
- **Engines.** Both engines and `S9b_exports.py` are byte-identical to `86097a44`.
- **Prompt and legs.** The prompt is `_scratch/s9b_build/s9b_repair_build_review_prompt_r9.md` (sha256 `5d5d369b…`).
  The legs were a fresh Claude agent (opus) and Grok (`grok-4.7`, effort high). Their reports are
  `_scratch/s9b_build/s9b_repair_build_review_r9_claude.md` and `…_r9_grok.txt`, with evidence directories beside
  them. Both returned **NEEDS REVISION**.

## Findings

| # | Leg | Finding | Disposition |
|---|---|---|---|
| 1 | Claude 1 | **SymPy K10 prints "undecided" from attempts that cannot return a witness at a mathematical part, and the printed evaluations were not made at the printed inputs.** `finite_attempt` returns "witness" only for text parts (line 807). A mathematical part returns "no change" or "undecided", with a fixed sentence saying no witness is inferred (lines 883 and 893). K7 and K10 skip `direct_search` (line 904), and round 9 removed the `cautious_condition_search` K10 used at `86097a44`. So none of K10's 4,163 attempts could produce a mathematical witness. The leg's guarded script decided the gate of the harness's own preferred part (`PY_S9B_A_DEFLECTION_G000`) at one input. Every conjunct except the negated `real_domain` decides True in both copies, and the gate is True in the baseline and False in the corrupted copy. The attempt records also print a point assigning `b` and `c_0`, while the printed evaluations still hold both symbols inside `ConditionSet` conditions. `subs` does not enter a `ConditionSet` condition; the leg's `xreplace` does. Grok's "confirms" on this row checked the 1,141 no-change decisions, which both legs found sound; it did not test whether a witness could be found. | **Confirmed.** It is round-8 Finding 1 again, in the other engine: an outcome printed from a search that cannot succeed. The undecided outcome misses the user's bar ("R4 at every knife's target in both engines"). The printed evaluations state inputs that were not substituted (E1). |
| 2 | Claude 2 | **The Wolfram K10 mutation is not the negation of the computed branch-existence condition.** The directive defines K10 as "the computed branch-existence condition is replaced by its logical negation before the observables are gated on it" (line 393). The spec defines branch existence as the condition under which the branch is real and propagating everywhere along the ray (line 193). Wolfram negates the local predicate (engine line 166). That predicate enters `propagationGate` inside the along-ray `ForAll` (lines 174 and 261–265), so the mutated gate is ∀r ¬P, not ¬(∀r P). It also changes the local construction's `Domain` (line 203) and the emitted `BRANCH_EXISTENCE` (line 275). SymPy negates `real_domain`, the computed along-ray condition (line 696), at line 770, after emitting `BRANCH_EXISTENCE` at line 720. The Wolfram K10 witness (local `Domain`, False vs True) is correct as computed, but it is a witness for a different mutation. Grok checked the conjunct values, which are correct; it did not check the mutation against the knife's definition. | **Confirmed.** It changes what may be claimed: the two engines' K10 runs do not ablate the same check. It is a harness repair: item 14 runs a copy "with exactly that one mutation at the named construction site", and the knife names the computed condition. The engine stays unchanged. My round-4 note had already carried this as an "engine placement" reading. It is a logical-form difference, not only a placement. |
| 3 | Claude 3 | **The rewritten `S9b_repair_wolfram_round8_measurements.md` says the round-9 files "record the actual resource stop and incomplete coverage" (line 5).** The round-9 measurements say all four final runs exited 0, with no guard, memory, task-limit or item-12 event (line 116). | **Confirmed.** It is a false statement in a measurement file. Correct it in the repair. |
| 4 | Claude 4 | **Wolfram excludes `Profiles` as a name list extracted from an expression, but the engine writes it as a literal list in source.** The harness lists it in `extractedNameKeys` (line 397). The engine writes `"Profiles" -> {deltaProfile,velocityProfile,xiProfile,rhoBr,fProfile}` (line 401). Every `Profiles` entry is excluded with exclusion 1's reason. | **Confirmed, minor.** No knife outcome changes. The printed exclusion reason is untrue, and item 14 allows exactly two exclusions. Fix in the repair. |
| G1 | Grok 1 | **Wolfram `IntegrationConstants` should be excluded as a symbol list, like `FreeInputs`.** | **Rejected.** Item 14's first exclusion covers lists of names "extracted from expressions" (directive lines 303–305). The engine builds `IntegrationConstants` as `If[homogeneous === 0,{}, {velocityIntegrationConstant}]` (lines 537 and 542). That is a computed choice of whether an integration constant exists, not a list extracted from an expression. Keeping it as a value site part follows the text. Claude reached the same reading. |

## Checked and found sound by both legs
- The marks are literal in both engines and load-bearing: FORM ablations by both legs flip marks.
- Every printed witness for K1–K9b and K11 in both engines re-checked at its printed inputs (both legs' tables).
- The SymPy software checks print each fixture's mark or value before raising (round-8 Finding 2 is folded).
- The SymPy K10 no-change decisions: 1,141 of 1,141 have literally identical operands.
- K11 coverage is printed: SymPy reaches 3 tags and Wolfram 10, and the unreached tags are labelled as item 14 states.
- The SymPy inventory matches its run, with no payload copied.
- The Wolfram K6, K7 and K10 searches now exist and print what was evaluated (round-8 Finding 1).
- Engines unchanged; reuse is clean; nothing leaks.

## Outcome
**Not clear.** Findings 1 and 2 change what the harnesses print for K10 and what may be claimed for it in both
engines. Findings 3 and 4 change printed statements.

**What must be true after the repair.**
- **SymPy:** each K10 site part gets an attempt that can return a witness: an engine-CAS-decided difference at one
  printed input, with every attached condition decided. Every printed evaluation is made at the printed inputs,
  including conditions inside `ConditionSet`. K10's outcome is re-derived from those attempts, and nothing else
  changes beyond what that requires.
- **Wolfram:** the K10 mutation is the negation of the engine's computed branch-existence condition, the same object
  SymPy negates, before the observables are gated on it. `Profiles` is read as the engine writes it. The round-8
  measurements sentence is corrected.
- Both harnesses rerun from final code, and the measurements and reports are regenerated.

The round-9 harnesses and files, as reviewed, are preserved before the repair overwrites them.

## The user's decision (2026-10-17)

The user chose through the question tool; the option texts are mine. Question: "Both reviewers want revisions. The two
main problems are both in K10. SymPy's K10 search can't actually find a witness (a reviewer found one at the
harness's own preferred part), and Wolfram's K10 negates a different condition than SymPy's, so the two engines
don't test the same knife. How should we finish?" Chosen: "Repair both (Recommended)", whose description reads: "Fix
SymPy's K10 search so it can find a witness and records the real inputs. Fix Wolfram's K10 to negate the same
condition SymPy does, plus two small text fixes. SymPy and Wolfram rerun in parallel, then a fresh two-leg review.
About 12–14 h, mostly SymPy's 8–10 h run; another round is possible." The other option was to repair Wolfram only
and record SymPy's K10 as not shown.

So repair round 10 rebuilds both harnesses to the "What must be true after the repair" list above. The same builder
session continues per engine. For each engine this is its second repair after its round-8 build, so G4 does not
yet call for a new author. The two guarded jobs run in parallel, sized to fit together: SymPy at 18 GiB, from its measured
round-9 peak of 15.6 GiB, and Wolfram at 8 GiB.
