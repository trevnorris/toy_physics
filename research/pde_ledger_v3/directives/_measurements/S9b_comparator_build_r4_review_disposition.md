# S9b cross-engine comparator, build review round 4: dispositions (orchestrator)

The comparator is Codex-written, so the legs are a fresh Claude agent and Grok (`CLAUDE.md` G1). It is an instrument,
reviewed until clear (G4). This round reviews the fourth repair, by the third repair's author (round-3 dispositions:
`S9b_comparator_build_r3_review_disposition.md`). Lookups: `S9b_comparator_build_r4_review_disposition_lookups.md`
(mechanical; generator `_scratch/s9b_build/gen/s9b_comparator_r4_lookups.sh`). Operand values the legs computed stay
in their evidence and are not copied here.

**Artifacts.** The repaired comparator (sha256 `326d204b…`), its tests (`efe7249b…`), the builder's report with its
two appendices and helper scripts, and the guarded runs `repair-04/full-02` (exit 0, 11,116 s) and `repair-04/tests-02`
(65 tests, exit 0), staged under `_scratch/s9b_build/comparator_r4/` and `_scratch/s9b_build/comparator_output/repair-04/`
(generator `_scratch/s9b_build/gen/s9b_comparator_r4_stage.sh`). Baseline: `_scratch/s9b_build/s9b_comparator_build_review_baseline_r4.sha256`
(all 33 entries verify). Prompt: `_scratch/s9b_build/s9b_comparator_build_review_prompt_r4.md` (sha256 `5479dab5…`).

**Legs.**
- **Fresh Claude (opus): CLEAR, no findings.** Report: `_scratch/s9b_build/s9b_comparator_build_review_r4_claude.md`.
  Evidence: `_scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/`.
- **Grok: CLEAR, no findings.** Report: `_scratch/s9b_build/s9b_comparator_build_review_r4_grok.txt`. Evidence:
  `_scratch/s9b_build/s9b_comparator_build_review_r4_grok_evidence/`.

Both legs found the four round-3 requirements hold on the real streams and in the tests: test profiles of different
functional forms for every live profile, with the round-3 one-sided brane-density corruption now moving a formed
residual (C1); the six small-magnitude exponential samples forming (C2); the label map as name-table rows citing the
engine lines, with a repoint control (G1); a cache-key control (O1). Both found sound: their own enumerations (330
joined, 49 and 44 local, none unjoined) match the accounting, and the leaf identities hold on every joined row; the
name table is injective and its cited lines denote the same objects; each leg's recomputation of formed samples from
the raw transcripts, in both engines' languages for Claude, matches the printed lines; one-sided corruptions of
load-bearing values move the printed residuals; no target token appears, the run exits 0 and joins no ablation
transcript or export delta; the inputs named in the configuration record are the committed transcripts. Neither leg
found anything narrowed, reconciled or regressed by the fourth repair.

| # | Observation (leg) | Disposition |
|---|---|---|
| a | The tolerance is relative to the integral (line 1208). The fourth repair makes formation independent of the integrand's amplitude, but an integral that cancels to about zero over an O(1) integrand still never forms; a probe prints unresolved coverage on all 33 samples. (Claude 1) | **NO REPAIR; claim limit.** Lookups: line 1208; the round-3 requirement as dispositioned (magnitude); the probe line; 0 occurrences of the reason the probe prints in the full-02 output; 183 and 27 occurrences of the two pole reasons, the leg's 105 unresolved value samples, all in the vanishing-brane-density stratum. No full-02 sample is affected, and a cancelling sample prints unresolved, not a residual. The S9b record may state the C2 property as amplitude-independence only. |
| b | No test covers the Wolfram reader's operator precedence or chained comparisons: ablation F1 (unary minus binds looser than `+`) and F2 (a chained comparison drops its first link) pass all 65 tests; so does F6 (the two density-label rows swapped). (Claude 2) | **NO REPAIR.** Lookups: F1, F2 and F6 pass, F3–F5 fail; under F1 3 lines report a moved residual on the real streams, under F2 none (only the stratum inventory changes, over 3276 chained sites in domains already unresolved); the label rows' cited engine lines. The production reader is checked on the real streams by both legs' independent recomputations, so the full-02 output is unaffected. F6 passes by design: item 6's repoint and swap controls are relative, and the table's orientation rests on the cited lines, which both legs read. A reuse of the comparator on other transcripts does not inherit this check; it needs its own recomputation. |
| c | The order of the domain-stratum inventory differs between runs; the multisets are equal. (Claude 3) | **NO REPAIR.** Lookups: multiset equal on every listed path. Not a defect in what is computed. A rerun is therefore not byte-identical to full-02, so the reviewed full-02 run itself is filed as the production output. |
| — | Grok, checked: evaluating the right operand under the exponential profile when the inverse is requested passes the 65 tests and moves the printed residual. | **NO REPAIR,** as in round 3: not a required item-6 mutation, and the printed output catches it. |

**Record caveats, not repairs.** What the printout leaves unformed limits what the S9b record may claim: the record
claims a comparison only where full-02 prints one formed. Carried from earlier rounds: no along-ray domain predicate
forms, because the Wolfram domains carry free turning radii; the 18 `B_GAMMA` rows compare only their markers; the
vanishing-brane-density stratum leaves its rows unformed with an interior denominator zero. The accounting rows with
no formed value comparison, and each remaining reason with its count, are in the builder's coverage appendix.

**Route.** Clear (G4): nothing outstanding changes what the comparator computes or what may be claimed from full-02.
The fourth repair is **accepted**. The comparator, its tests and the builder's report are committed in the main tree,
and the reviewed full-02 output is filed as the production output (record: `S9b_comparator_production_run.md`).
