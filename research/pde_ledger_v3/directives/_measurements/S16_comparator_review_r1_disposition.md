# S16 comparator build review, round 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S16_comparator_review_r1_lookups.md`, beside this file
(generator `_scratch/s16/cmp_review/gen_r1_lookups.sh`). "Q-n" cites a section of that file. Round 0's dispositions
are in `S16_comparator_review_r0_disposition.md`.

## The round

- **Artifact.** The comparator, its tests, the run record and the builder's report, as repaired by Codex
  (`gpt-6-astra`, high, a fresh session) from the version reviewed in round 0 and preserved at `98936790`. The repair
  brief is `_scratch/s16/cmp_repair1/brief_s16_comparator_repair1.md`. It points at round 0's "What must be true" for
  K1–K4 and asks for nothing else.
  - The production output is `_scratch/s16/cmp_work/repair_r1_full_01.jsonl`, with its summary and the last test
    output `repair_r1_tests_final_01.json`.
  - The builder's report is the repair's final message, copied verbatim from its transcript.
- **Prompt and legs.** `_scratch/s16/cmp_review/s16_cmp_build_review_prompt_r1.md` (Q1). It is round 0's prompt
  with an item 0 added: K1–K4, and whether anything else changed against `98936790`. Codex wrote the artifact, so the
  legs are a fresh Claude agent (Opus) and Grok (`grok-4.7`, high) (G1). They had the identical prompt and ran in
  parallel. Grok exited 0.
- **Baseline.** All seven baseline files matched after the round. The working tree was unchanged (Q1).
- **Verdicts.** Claude: NEEDS REVISION, one finding. Grok: CLEAR, no findings (Q2).

## Round 0 findings

Both legs report each of K1–K4 repaired. Each leg restored the round 0 behaviour on a copy; the targeted test failed
and the production output moved. Neither leg reports a change beyond the repair (Q5).

- **K1.** On binding 182, every sample residual is now 0, in every stratum. The only comparisons with a sample
  residual other than 0 are `domains:0` and `domains:1`, the condition comparisons (Q4).
  - Claude: restoring the simultaneous substitution fails `volume_partial_time` ("partial time became chart time").
    It also makes bindings 182 and 356 nonzero again.
  - Grok: restoring it raises the same assertion.
- **K2.** In the zero-count stratum, both `N(t)` (`definitions:37`) and the printed count (`definitions:38`) are 0
  (Q4).
  - Grok: the density's derivatives stay live there (`Derivative(Sigma, x1)` stays `47281/42`).
  - Both legs: restoring a generic `Σ` in that stratum fails `count_zero_fields` ("printed count does not vanish").
- **K3.** Of 885 exact residuals, 655 are 0 as computed and 228 more resolve to 0. Two resolve to something else:
  `domains:0` and `domains:1`, the condition comparisons. Of the 114 as-computed residuals containing `Subs(0,`, none
  keeps it after resolution (Q4).
  - Both legs: replacing the resolver with the identity fails `wrapper_resolution`.
  - Grok: a slot change still resolves to a different derivative.
- **K4.** The name table cites SymPy 101, 102 and 103 for `N`, `D` and `P`, the lines that introduce them (Q4).

## Findings

**K5. Many binding rules cite SymPy lines that do not build the bound object.** (Claude 1.) **CONFIRMED; it changes
nothing computed.**
- The printed rule map (the configuration record, Q3) cites SymPy lines for each rule. The engine is 439 lines long
  and unchanged since `1d990ee6` (Q3).
- Cited lines that are outside the object (verbatim lines, Q3):
  - `definitions` cites 208–211. Those lines build the momentum rate parts. `count.definition` is emitted at 201, and
    the dipole and momentum definitions at 203–204.
  - `momentumReceived` and `momentumRate` cite 218–222, a comment and `def integral_time`. They are built at 211–214.
  - `momentum_surface_motion` cites 223–224, inside `integral_time`. It is emitted at 215–216.
  - `force_companion` cites 198, the moment index list. It is emitted at 193.
  - `partition_routes` cites 293, a comment. It is written at 289.
- Cited ranges past the end of the file: `coverage` cites 446–448, and `dimensions` cites 398–445 (Q3).
- Reason strings (Q3):
  - "py 205-207", for the moment index encoding, cites the count and dipole rate builds. The indices are at 198–200.
  - "py 418-445", for the dimension-sum audit, runs past line 439.
  - These strings appear in output records: 60279, 22 and 2 records for "py 418-445", "py 183-185" and "py 205-207".
    They do not appear literally in the run record.
  - Claude also reports "py 183-185" and several rules as partly wrong. The lookups above do not settle each of those
    one by one, and the repair covers them all the same way.
- The citations are only printed. They are stored in `RULES` and the surface table, and written to the configuration
  record. No comparison reads them (Q3).
- The same citations were in the round 0 version (Q3). Round 0 did not report them; the repair did not introduce them.
- Every Wolfram citation is correct (Claude).
- The list, item 2: "every join made through the table, generated from rules that each cite the object both ends
  denote and the line in each engine that builds it". Item 3: each entry cites "the line in each engine that
  introduces it".
- **What may not be claimed now.** That the printed binding table's SymPy citations ground the bindings. The bindings
  themselves rest on the legs' measurements, not on the citations: the independent comparisons, the form ablations
  under which the production output moves, and the one-sided corruption that moves only its own binding.
- **What must be true.** Each rule's and each reason string's SymPy citation names the lines that build or print the
  object, inside the engine.

## Not findings

- **The tests stay green when the production tables are ablated** (Claude): surfaces by index, shared names, `Π`
  transposed, shifted components, merged pieces, alias repointed, components by position. Under each of these, the
  production output moves. This is round 0's "not a finding": item 7's controls are fixture-level.
- **One production ablation exits 1.** Claude collapsed applied functions inside the reader; the production run raised
  `TypeError`. Collapsing at the comparison stage runs, and prints the affected samples as unresolved. The collapse is
  never silent.
- **One-sided corruption** is sound. Claude's change of a single SymPy piece (`count.surface1.radiusMotion`) moved
  only binding 61, that tag's leaf counts and the input sha. Grok's `+1` moved only its own binding.
- **Reproducibility** is sound. Each leg's guarded full rerun reproduced every record, apart from runtime and memory.
- **Compare first.** Claude's 306 evaluated residuals are all 0, and it reproduced the printed sample values for 12
  bindings. Grok's own comparison matched the production records for its objects.

## Outcome

One finding, K5, confirmed. It is the class of K4, in the binding table. It changes nothing computed and limits one
claim: what the SymPy citations ground. Both legs clear the comparison itself.

The physics is already verified, so how to finish was put to the user. The options were:
- a citation-only repair with a scoped check;
- accepting now with K5 carried;
- the repair with a full round 2.

The user chose the citation-only repair with a scoped check (2026-10-24).
- The reviewed version is preserved with these dispositions.
- Codex repairs K5 to what must be true and changes nothing else. It reruns the tests and the full comparison.
- A fresh Claude agent and Grok then check, scoped:
  - each SymPy citation against the engine;
  - that the code diff against this preservation commit touches only citations;
  - that the new output differs from `repair_r1_full_01.jsonl` only in citation strings and runtime fields.
