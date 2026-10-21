# S14a part 1 comparator build, review round 0: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S14a_1_comparator_build_review_r0_lookups.md`, beside this
file (generator `_scratch/s14a/gen/cmp_build_r0_lookups.sh`). "E-n" cites a section of that file. Every CAS result
quoted here is a leg's saved stdout.

## The round

- **Artifact.** The comparator build to the decision list at `43174e8c`, as the builder left it. Codex
  (`gpt-6.1-sol`, effort xhigh, session `01a15f56…`) exited 0 after 302,055 tokens. Its final guarded runs,
  `tests_04` and `full_02`, both exited 0 (E1). Reviewed as the eight files of the round-0 baseline:
  - `scripts/S14a_1_cross_engine_comparator.py`;
  - `scripts/test_S14a_1_cross_engine_comparator.py`;
  - the run record and the builder report in `directives/_measurements/`;
  - `_scratch/s14a/cmp_work/full_02.jsonl`, with its inventory and unpaired-components files, and `tests_04.out`.

  All eight matched the baseline after the round. The working tree was unchanged apart from the lookups file itself
  (E1).
- **Prompt and legs.** `_scratch/s14a/cmp/review/s14a1_cmp_build_review_prompt_r0.md` (E1). The comparator is
  Codex-written, so the legs are a fresh Claude agent (Opus) and Grok (`grok-4.7`, effort high) (G1). The comparator
  runs no Wolfram kernel, so the legs ran in parallel, with the identical prompt. Reports:
  `/var/projects/_review_scratch/s14a1-cmp-r0-claude/REVIEW.md` and `/var/projects/_review_scratch/s14a1-cmp-r0-grok/REVIEW.md`.
- **Verdicts.** Claude: NEEDS REVISION (F1, F2, F3). Grok: NEEDS REVISION (1) (E1).

## Findings

**H1. Nested integrals: the reader builds an object neither engine printed.** (Claude F1.) **CONFIRMED.**
- **The list.** Item 1: "the reader builds the same object from both … The exact residual is then the residual of the
  printed objects", and "A different object is still read as different: another derivative slot, order or point;
  another integrand, variable or range."
- **The gap.**
  - Line 788 renames each integral's own variables `_int0`, `_int1`, …, restarting at every nesting level (E2). So an
    outer integral's variable gets the same name as the inner one, and the inner binder captures every place where the
    inner integrand uses the outer variable.
  - Claude's fixture goes through the production path. It reads two different nested integrals as the same tree, with
    residual 0; integrated directly, the two give different values (E2).
  - Five bindings have an outer and an inner bound variable that share one symbol after reading, in both engines (E2):
    - b169, order lost;
    - b171, storage;
    - b173 and b174, the shell balance;
    - b175, face channel minus order lost.
  - One-sided corruptions of the real Wolfram transcript swap a slot inside the nested integral of storage (c2a) and
    of order lost (c2b). Neither moves any binding (0 of 227). A copy that changes only line 788 moves exactly the
    corrupted binding in each case (E2).
  - The sample values printed for b169 and b175 differ from the values of the printed objects (E2).
  - Claude reports that the five bindings still agree once their bound variables stay distinct. Against production,
    that copy moves 9 bindings' printed output (E2).
- **What must be true.**
  - After renaming, an outer integral's variable stays distinct from every inner bound variable, so the reader builds
    the printed object.
  - A test fails when two nested integrals that differ only in which variable the inner integrand uses compare equal.

**H2. Evaluation cannot see a face derivative taken at another point, and for three parts evaluation is the only
resolved check.** (Claude F2.) **CONFIRMED.**
- **The list.**
  - Item 1: "another derivative slot, order or point" is read as different.
  - Item 5: evaluation is "added beside" an exact residual, "never in its place", and stands alone only where exact
    subtraction is not feasible.
- **The gap.**
  - The test functions are additively separable: each form is a constant plus a sum of one-argument terms. The faces
    are built from the same forms (E3). Claude's probe finds every mixed second partial of a face to be zero, in every
    form and stratum (E3).
  - So a face's time derivative takes the same value at every radius. A face derivative evaluated at another point
    changes no sample.
  - For storage (b171) and the shell balance (b173, b174), the exact residual prints as two integrals of opposite
    sign, over the same variable and range, left uncombined (E3). Their samples are all zero, so the samples are the
    only resolved check on these parts.
  - One-sided corruption c3 moves storage's face derivative from the integration variable to `inner_radius`. It
    changes the exact text and no sample (E3). The comparator's resolved output therefore cannot tell the corrupted
    storage from the accepted one.
  - Item 5's sentence "its derivatives in different arguments differ" is met as written. The gap is item 1's, wherever
    evaluation stands in for an exact residual that is not resolved.
- **What must be true.**
  - For every part whose exact residual is not resolved to zero, one of these holds:
    - a face or field derivative taken at another point, in any argument, moves the samples in at least one form and
      stratum; or
    - the exact residual resolves.
  - A test fails when a derivative taken at another point, in another argument, moves no sample.

**H3. An unjoined tag's reason is chosen by a pattern in its name, and can be false for that tag.** (Grok 1; Claude
F3.) **CONFIRMED, minor.**
- **The list.** Item 2:
  - "An engine's own checks (its comparison notes, operand pairs and residuals) … Otherwise each is listed as
    unjoined with that reason."
  - "It is left unjoined only where the cited lines build different objects, and its reason says so."
- **The gap.**
  - `unjoined_reason` (l.1305–1339) picks a reason by substrings of the tag name and falls back to an engine-wide
    default (E4). `is_label` (l.494–495) recognizes a comparison note by the suffix `/comparison` (E4).
  - Wolfram line 9, `transverse integration / comparison`, is a comparison note with a spaced separator. It is the only
    one in the stream (E4). It misses both suffix tests: it is not read as a label, and it is listed under the
    coordinate-construction reason (Grok, E4).
  - Claude's ablation v08b removes the 125 bindings that have a component end. The freed tags still get reasons that
    state their counterparts are bound or absent, and those statements are now false (E4).
- **Effect.** No compared value moves. What it affects is whether item 2's accounting supports what it states. A
  reason may assert a fact about the other stream that nothing computed. So "every object the S14a part 1 decision
  list names that both engines print is compared" rests on how the table was curated, not on the accounting.
- **What must be true.**
  - Each unjoined tag's printed reason is true of that tag. The note at Wolfram line 9 is listed as an engine-local
    comparison note, and read on the label path.
  - Where a reason states a fact about the other stream, that fact is either computed from the parse and the binding
    table, or the reason says it rests on the table's curation.
  - A test fails when a tag freed from a binding still receives a reason stating that its counterpart is bound or
    absent.

## Not findings

- **The shell domain (b165).** The engines print different domains: SymPy `inner_radius < outer_radius`, Wolfram
  `Inactive[Less][0, ri, ro]` (E5). The exact residual prints their exclusive-or; the samples read `False` because
  every point has a positive inner radius. Both legs saw this. The difference is printed, not hidden, as item 4 asks.
  The record interprets it.
- **Reproducibility.** Both legs' guarded reruns of the final code reproduce `full_02.jsonl`. The only difference is
  the resources record (E6).
- **Operations.** Two of the Claude leg's guarded probes were refused admission, because it launched one while its
  own queue was running. It reran the probe alone. No job was killed. No effect on the evidence.

## Outcome

Three findings, all confirmed: H1 and H2 change what the comparator computes or what may be claimed from its output,
and H3 limits what its accounting may be cited for. The reviewed build is preserved before the repair. The repair goes
to the same author, in a resumed session, with these dispositions. Then review round 1.

For the record: completeness of the comparison of named objects rests on the binding table's curation, unless the
repaired comparator computes it (H3).
