# S9b cross-engine comparator, build review round 0: dispositions (orchestrator)

The comparator is Codex-written, so the legs are a fresh Claude agent and Grok (`CLAUDE.md` G1). It is an instrument,
reviewed until clear (G4). Lookups: `S9b_comparator_build_r0_review_disposition_lookups.md` (mechanical; generator
`_scratch/s9b_build/gen/s9b_comparator_r0_lookups.sh`). Operand values the legs computed stay in their evidence and are
not copied here.

**Artifacts.** `scripts/S9b_cross_engine_comparator.py`, `scripts/test_S9b_cross_engine_comparator.py`, the builder's
report `S9b_comparator_build_report.md`, and the guarded full run `full-04` (`comparison.jsonl` sha256 `bf312dc2…`,
1,154,909,769 bytes; `accounting.jsonl` sha256 `1bf34cbc…`), hard-linked under `_scratch/s9b_build/comparator_output/`.
Baseline: `_scratch/s9b_build/s9b_comparator_build_review_baseline_r0.sha256`. Prompt:
`_scratch/s9b_build/s9b_comparator_build_review_prompt_r0.md`.

**Legs.**
- **Fresh Claude (opus): NEEDS REVISION, four findings.** Report: `_scratch/s9b_build/s9b_comparator_build_review_r0_claude.md`.
  Evidence: `_scratch/s9b_build/s9b_comparator_build_review_r0_claude_evidence/`.
- **Grok: NEEDS REVISION, two findings.** Report: `_scratch/s9b_build/s9b_comparator_build_review_r0_grok.txt`. Evidence:
  `_scratch/s9b_build/s9b_comparator_build_review_r0_grok_evidence/`.

Both legs found sound: their own enumerations of the two transcripts (379 SymPy and 374 Wolfram tags; 330 joined, 49
and 44 local, none unjoined, no duplicates) match the comparator's accounting; the name table is injective, and each
pair denotes the same spec object; operands print before any guard; `not_formed` is distinct from a parse failure; no
target token appears, and the run exits 0. Their one-sided corruptions of a stream copy moved the printed output.

| # | Finding (leg) | Disposition |
|---|---|---|
| D1 | On nine G121 rows (`A_`, `F_FLOW_A_`, `F_LIVE_A_` of `ROUND_TRIP`, `ONE_WAY_ER`, `ONE_WAY_RE`), 312 samples print the SymPy operand as undefined and the Wolfram operand as a number. The SymPy object is a second derivative in `b` of an integral with a `b`-dependent endpoint. The comparator's `Subs(...).doit(deep=False)` differentiates under the integral sign, which leaves a `zoo` boundary term and an integrand that cannot be integrated; `_numeric_value` then prints undefined. Claude's own evaluation of the object is finite. Every such sample prints `numeric operands printed`, although no numeric comparison was formed, and the summary does not count them. (Claude 1, Grok 2) | **ACCEPT.** Lookups: lines 711, 741 and 829; the nine rows and their counts. A sample prints undefined only when the object is undefined there. When the comparator's evaluation cannot produce a value for an object that is defined, the sample prints as unresolved coverage. A sample's comparison outcome says what was formed. |
| D2 | Under `brane_zero_at_point`, the forward integrals of `F_LIVE` rows have a non-integrable singularity inside the range, yet print as finite values, with 18 nonzero residuals. The quadrature checks only that its result is finite. (Claude 2) | **ACCEPT.** Lookups: line 757; the census counts. A numeric value prints only where the integral exists and the quadrature result is shown reliable at the printed precision. Otherwise the sample prints as undefined or as unresolved coverage. |
| D3 | The mutation controls do not fail when the step they target is disabled. `mutation()` disables all of `Comparator.leaf`; `products()` keeps the symbolic `residual`, which for an unresolved functional echoes the operands. Three FORM ablations of the comparator pass the tests: numeric evaluation disabled (28 of 32 pass), profile derivatives dropped when sampling (32 of 32, and the printed output is blind to it as well), and the zero-density probes removed (32 of 32). The expected configuration set comes from the comparator's own `PROFILE_NAMES`. (Claude 3) | **ACCEPT.** Lookups: the test lines and the three ablated test runs. Each control fails when the comparison step it targets is disabled, through the sampled residuals, as comparator item 6 asks. The expected configuration set does not come from the comparator's own list. |
| D4 | On six `B_GAMMA_DIFFERENCE` rows, the SymPy on/outside gate is not recognized, because the reader turns the on-part into a record and the gate rule needs a sequence. The `NOT_ESTABLISHED` markers are never compared. (Claude 4) | **ACCEPT.** Lookups: reader lines 108–111, extraction lines 570–573, the six rows. The gate is recognized whatever form the reader gives it. |
| D5 | `A_NONRECIPROCAL_PATH_DEPENDENCE`: both engines emit `False`. Its condition piece prints as `unresolved coverage`, because a piece counts as sampled only when a probe finds it true. Grok's fix: mark a piece evaluated whenever its probes return booleans. (Grok 1) | **ACCEPT IN PART.** The piece is the literal condition `False`, which is empty, so "unresolved coverage" misstates it. Grok's fix is not adopted: a piece that no probe falls in and that is not shown empty is not covered, and item 4 requires the sample to meet every piece. A piece shown empty prints as empty, distinct from unresolved coverage. |
| O1 | No attached domain stated through a set or a quantifier is compared. The translator has no constructor for `ConditionSet`, `ForAll` or `Exists` and raises `unsupported constructor`, so 270 rows print a domain as `not_formed` for that reason, and one more for `Inactive/ForAll`. Comparator item 3 requires a domain comparison wherever a payload attaches a domain, and item 4 requires both conditions' truth values at each sampled point, with unresolved coverage only where the comparator cannot evaluate. A constructor the engines use for attached domains is not such a case. Neither leg raised it; Claude recorded it as scope. (orchestrator) | **ACCEPT** (orchestrator, from lookups). Lookups: line 543; no `ConditionSet`, `ForAll` or `Exists` in the source; the census of `not_formed` reasons. Each attached domain gets a domain comparison as items 3 and 4 say. Where a sampled truth value cannot be decided, that sample prints as unresolved coverage with its reason. The builder report lists each remaining `not_formed` reason with its count and the spec objects it leaves uncompared. |

**Record caveats, not repairs.** Claude found the printed `approx_30` operands accurate to about 23 significant
digits, and the residuals measuring the two integrands' agreement at the quadrature nodes rather than the operands'
accuracy. Whatever remains `not_formed` after the repair limits what the S9b record may claim.

**Route.** The same builder repairs the comparator (first repair round; the author-change condition is not met). The
repair is reviewed by fresh legs.
