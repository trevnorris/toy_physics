# S14a part 1 comparator build, review round 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S14a_1_comparator_build_review_r1_lookups.md`, beside this
file (generator `_scratch/s14a/gen/cmp_build_r1_lookups.sh`). "E-n" cites a section of that file. Every CAS result
quoted here is a leg's saved stdout.

## The round

- **Artifact.** The comparator after repair 1 of the round-0 findings H1–H3
  (`S14a_1_comparator_build_review_r0_disposition.md`; the reviewed round-0 build is preserved at `88bf7bf9`). The
  repair went to the same author in a resumed session, which exited 0 after 615,948 tokens (E1). Its guarded runs,
  `r1_tests_01` and `r1_full_01`, both exited 0 (E1). The comparator and test changed by 242 insertions and 17
  deletions against `88bf7bf9` (E1). Reviewed as the nine files of the round-1 baseline:
  - `scripts/S14a_1_cross_engine_comparator.py`;
  - `scripts/test_S14a_1_cross_engine_comparator.py`;
  - the run record and the builder report in `directives/_measurements/`;
  - `_scratch/s14a/cmp_work/r1_full_01.jsonl`, with its inventory and unpaired-components files, `r1_tests_01.out`
    and `r1_run_inventory.json`.

  All nine matched the baseline after the round. The working tree was unchanged apart from the lookups file itself
  (E1).
- **Prompt and legs.** `_scratch/s14a/cmp/review/s14a1_cmp_build_review_prompt_r1.md` (E1). It asks first whether H1–H3
  hold, by undoing each repair and by one-sided corruption, then for the whole comparator against items 1–8 of the
  decision list. The comparator is Codex-written, so the legs are a fresh Claude agent (Opus) and Grok (`grok-4.7`,
  effort high) (G1). The comparator runs no Wolfram kernel, so the legs ran in parallel, with the identical prompt.
  Reports: `/var/projects/_review_scratch/s14a1-cmp-r1-claude/REVIEW.md` and
  `/var/projects/_review_scratch/s14a1-cmp-r1-grok/REVIEW.md`.
- **Verdicts.** Claude: CLEAR, no findings. Grok: CLEAR, no findings (E1).

## The round-0 findings, checked

**H1. Nested integrals.** **Holds.**
- Binders are now numbered by nesting depth (l.719–735), so an inner integral cannot take an outer variable's name
  (E2).
- One-sided slot changes inside a nested integral each move exactly one binding: Claude C1 (Wolfram, order lost) moves
  b169, Claude C2 (SymPy, storage) moves b171 (E2). Grok's slot change on Wolfram storage moves b171 at 72 sample keys,
  in all three forms (E2, E7).
- Undoing the repair moves b169, b171, b173, b174 and b175 in operands and values, and fails
  `test_nested_integral_variable_identity` in both legs (E2, E7). With the repair, Grok counts `_int1` twice in each
  engine's storage operand; undone, it counts none (E2).

**H2. A derivative taken at another point.** **Holds.**
- The quadratic form now adds `argument_coupling` (l.870–880), a sum over every pair of arguments, so each first
  partial depends on every other argument. The affine and cubic forms are unchanged (E3).
- Moving storage's face derivative from `r` to `ri` (Claude C3, Grok's face-point case) moves only b171. Its samples
  leave zero in the quadratic form, at 36 sample keys (E3). Field and face derivatives moved to another point in
  three other bindings each move only that binding: C4 b128, C5 b50, C10 b168 (E3). b50, b128 and b171 are three of
  the seven parts whose exact residual is not resolved (E9).
- With the coupling removed, Grok's same face-point change moves no sample (`sample_equal undo_h2 face_h2 True`); the
  exact text still changes (E3, E7).
- `test_every_live_first_partial_at_another_argument_point` loops over every name-table row with arguments, crossing
  every derivative slot with every argument point (E3). Undoing the repair fails it in both legs (Claude: 919 failures)
  (E3, E7).

**H3. Unjoined reasons.** **Holds.**
- `is_label` now recognizes a comparison note through `is_comparison_note` (l.499–500). Wolfram line 9,
  `transverse integration / comparison`, is listed as an engine-local comparison note and read on the label path
  (`read_as_label: true`) (E4).
- Of the 983 tags, 327 are bound and 656 unjoined. Of the 656 unjoined reasons, 343 begin "Reference-table/source
  curation (not computed active-table or counterpart-completeness evidence)", and 313 begin "Engine-local comparison
  note" (E4).
- Undoing the repair restores the old reason at line 9 and fails `test_freed_component_bindings_have_current_accounting`
  in both legs (Claude: 126 failures) (E4, E7).

## The repair on the real streams

Claude compared the round-0 output (`full_02.jsonl`) with the round-1 output (E5):
- no binding was added or dropped;
- no residual moved;
- 190 bindings changed only in operand and sample values (the quadratic form changed);
- 656 tag reasons changed;
- no tag's accounting status or parsed or compared leaf count changed.

## The whole comparator, items 1–8

- **Ablations.** Both legs ran the nine FORM ablations the prompt asks for, and Claude added one more (A8b: the
  component selector ignored). Every ablation fails at least one test (E6, E7). For four of Claude's, no residual
  moves on the real streams: A1 raw coefficient, A2b non-separating forms, A3 shared test function and A5 CAS endpoint
  weight. The targeted test still fails in each. Grok reports that the streams contain no closed-endpoint integral
  (E6, E7, E9). Claude's A4 (range dropped) fails `test_integral_range_variable_integrand` among three tests (E6).
- **One-sided corruptions.** Each moves exactly the binding it targets: C6 b167, C8 b126 (upper delta only), C9 b177
  (condition only). C7 changes a delta coefficient away from its face only, and moves nothing (E8).
- **Reading and accounting.** SymPy's 619 lines = 600 tag starts + 13 progress markers + 6 continuation lines.
  Wolfram's 397 lines = 383 + 14 + 0. Compared leaves equal parsed leaves for each engine: 21714 SymPy, 24948 Wolfram
  (E10).
- **No target.** The output and the source contain no verdict token, and the configuration prints "No residual
  target was given." (E10).
- **Reproducibility.** Claude's guarded rerun of the final code reproduces `r1_full_01.jsonl` apart from the resources
  record (sha256 `640b95d1…` for both). Grok's rerun gives equal samples; its resources record differs (E10).

## Carried to the S14a part 1 record (not findings)

- **Exact residuals not literally 0 or False** (E9).
  - Seven integral parts: b50, b72, b128, b148, b171, b173 and b174. For each, the exact residual as computed is
    printed, and the samples beside it are the resolved check. A derivative moved to another point moves these
    samples through the quadratic form only (E3).
  - b165, the shell domain: dispositioned in round 0.
  - b188 and b189, the projected-velocity and volume-flux dimensions: their "external printed domain contexts" part
    is not literally zero. SymPy's operand there is `True`. Wolfram's is `Ne(0, Integral(...))`: the projected
    density at the outer radius is nonzero (E9).

  The record interprets these.
- **Completeness.** That every object both engines print and the S14a part 1 decision list names is compared rests on
  the binding table's curation. The 343 curated reasons say so (E4).
- **The endpoint-delta convention** (item 4) is checked by the tests only: the streams contain no closed-endpoint
  integral (Grok, E9).

## Operations

- Grok ran one guarded job in pool `s14a1-cmp-review-r1`: exit 0, wall 2043 s.
- Claude ran its jobs one at a time in the same pool. None was refused admission or killed.

## Outcome

Both legs are clear. H1–H3 hold, each checked by undoing its repair and by one-sided corruption. Nothing outstanding
changes what the comparator computes or what may be claimed from its output, beyond the limits carried to the record
above. The comparator is accepted as reviewed in round 1, at the baseline sha256s. Its production output is the
reviewed run `r1_full_01` itself, filed as `scripts/out/S14a_1_cross_engine_comparator.out`, with
`S14a_1_comparator_production_run.md`.
