# S9b cross-engine comparator, build review round 1: dispositions (orchestrator)

The comparator is Codex-written, so the legs are a fresh Claude agent and Grok (`CLAUDE.md` G1). It is an instrument,
reviewed until clear (G4). This round reviews the first repair (round-0 dispositions:
`S9b_comparator_build_r0_review_disposition.md`). Lookups: `S9b_comparator_build_r1_review_disposition_lookups.md`
(mechanical; generator `_scratch/s9b_build/gen/s9b_comparator_r1_lookups.sh`). Operand values the legs computed stay
in their evidence and are not copied here.

**Artifacts.** The repaired comparator (sha256 `5f713350…`), its tests (`cdcbb953…`), the builder's report
(`3e9701c8…`) and the guarded full run `repair-01/full-03` (`comparison.jsonl` sha256 `24f5e65c…`; `accounting.jsonl`
`64e03060…`), staged under `_scratch/s9b_build/comparator_r1/` and `_scratch/s9b_build/comparator_output/repair-01/`.
Baseline: `_scratch/s9b_build/s9b_comparator_build_review_baseline_r1.sha256` (all 19 entries verify). Prompt:
`_scratch/s9b_build/s9b_comparator_build_review_prompt_r1.md`.

**Legs.**
- **Fresh Claude (opus): NEEDS REVISION, four findings.** Report: `_scratch/s9b_build/s9b_comparator_build_review_r1_claude.md`.
  Evidence: `_scratch/s9b_build/preserved/rev_r1_claude_evidence/`.
- **Grok: CLEAR.** Report: `_scratch/s9b_build/s9b_comparator_build_review_r1_grok.txt`. Evidence:
  `_scratch/s9b_build/preserved/rev_r1_grok_evidence/`.

Both legs found sound: their own enumerations of the two transcripts (330 joined, 49 and 44 local, none unjoined)
match the comparator's accounting; the name table is injective; value and domain are separated; the `ρ > 0` and
`ρ ≤ 0` parts and their labels are compared part by part; booleans are compared as values, and a one-sided flip moves
the output; `not_formed` outcomes are printed where representations differ; no target token appears, and the run
exits 0.

| # | Finding (leg) | Disposition |
|---|---|---|
| R1 | Whether a residual forms depends on its value. `_numeric_value` (lines 906–910) prints a quadrature value only when its 18- and 30-digit evaluations differ by at most `1e-14` of the 30-digit value, and the same rule is applied to the residual (line 1104). When the operands agree, the residual is near zero, the rule fails, and the sample prints `unresolved coverage` with the residual's computed values dropped; when they disagree, it forms. The run has 1,513 such samples on 61 rows (Part A's deflection, radar-slope, one-way and round-trip rows at several grades, and both `B_RESIDUAL_*` rows, across the base, forward and single-mechanism variants), and no formed nonzero residual. A one-sided coefficient corruption turned 24 of those samples into formed nonzero residuals. (Claude 1) | **ACCEPT.** Lookups: lines 906–910, 1038–1042 and 1103–1104; 61 rows of `comparison.jsonl` carry the reason; Claude's counts and its C1 ablation summary. Grok's own scan of `A_DEFLECTION_G001` records the same outcome (16 of 33 samples with matching operands and no formed residual); Grok did not classify it as a defect. The rule implements the round-0 D2 requirement, which was written for printed operand values, and the repair applied it to the residual as well. Required: whether a residual forms does not depend on its value. The residual's computed values are printed. Its reliability is judged at the scale of the operands it is the difference of, never at its own size. A control fails if two operands equal in value but different in form leave the residual unformed. |
| R2 | The sample points are not generic: at each of the three points `b` and `c_0` take the same value (lines 690–696). An engine difference that exchanges a factor of `b` for `c_0` is invisible: a one-sided FORM corruption of that kind left all 33 samples identical. (Claude 2) | **ACCEPT.** Lookups: lines 690–696; the coincidence census; the C4 ablation line; item 4 asks for generic rational points. Required: at each sample point, distinct symbols take distinct values, and a control fails unless exchanging one symbol for another in one operand moves the computed comparison. |
| R3 | `compared_leaves` counts every leaf the comparison step consumed, including leaves in comparisons where nothing formed (lines 1097, 1117, 1353–1355). About 99% of each engine's "compared" leaves are in comparisons where every sample is unresolved. (Claude 3) | **ACCEPT.** Lookups: the cited lines; Claude's totals, which reproduce the accounting. Item 6 asks for "the count actually compared". Required: the accounting separates leaves in a comparison that formed an outcome from leaves whose comparison was attempted but formed nothing, and a row whose comparisons formed nothing stays visible as such. |
| R4 | The maps from Wolfram function names to SymPy functions (`WL_BUILTINS`, lines 439–445, and the inline maps at lines 474 and 537–538) have no repoint control. Pointing `Sqrt` at a cube root left all 46 tests passing; the production output showed it as cross-engine differences that would be the comparator's own. (Claude 4) | **ACCEPT.** Lookups: the cited lines; item 2's scope ("Any map between the two engines' symbol and function names") and its repoint control; the FA1 test run and output counts. Required: each such map is an explicit, injective table under item 2, each entry names the function both names denote, and a control fails unless repointing any entry moves the computed comparison. |
| O1 | Grok: a Wolfram profile derivative is translated as `Subs(Derivative(…), _profile_argument, r)`. Where the SymPy side has the plain derivative, the symbolic residual stays unresolved on 12 rows, while every sample forms. Removing the `Subs` wrapper gave a symbolic zero, and the tests still passed. (Grok, under "checked and found sound") | **NO REPAIR.** Lookups: Grok's ablation lines. The unresolved symbolic residual is printed as unresolved and the samples form, so no claim rests on it. Removing the wrapper at a bare argument does not change the value, so the passing tests do not show a missing control. The record reports the 12 rows' symbolic outcome as printed. |
| O2 | Claude: the velocity integration-constant row binds two constants that multiply independently computed homogeneous solutions, so they need not denote the same object. `rho_bulk` and `fLocal` never meet a counterpart in a shared operand. Every row that carries them is not formed. (Claude, name table) | **ACCEPT IN PART.** Lookups: lines 282–283. Item 2 requires each table row to denote one spec object. Required: the repair report states, with its citation, whether the two integration constants denote the same object; if that is not shown, the row is removed and both names stay engine-qualified. `rho_bulk` and `fLocal` need no change. |

**Record caveats, not repairs.** Claude's "Coverage facts the record needs" stand as printed: no comparison formed for
the Part B and C conditions, the implied `j_n`, both effective `γ`s and their difference, the N-dependence,
`BRANCH_TYPE` and `F_MASS_SOLUTION`; every along-ray domain predicate is unresolved (the Wolfram domains carry free
turning radii); and 36 one-way and round-trip values at G021, G101, G120 and G121 are never evaluated on a live
profile. Whatever remains unformed after this repair limits what the S9b record may claim.

**Route.** The same builder repairs the comparator (second repair round). R1 is a defect in the material the first
repair added; the others predate it. The author-change condition (successive revisions breeding defects in the
material just changed) is not met by one instance; if the next review finds defects in what this repair changes, the
author changes. The repair is reviewed by fresh legs.
