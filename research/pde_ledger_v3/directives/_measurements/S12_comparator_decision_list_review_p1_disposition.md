# S12 comparator decision list, review pass 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S12_comparator_decision_list_review_p1_lookups.md`, beside this
file (generator `_scratch/s12/gen/cmp_dl_p1_lookups.sh`). "E-n" cites a section of that file.

## The pass

- **Artifact.** `directives/S12_comparator_decision_list.md`, as first written, uncommitted at review on top of
  `a3da43d9`. It is preserved unchanged in the commit that adds this file. Its hash matched the baseline after the
  pass, and the working tree was unchanged apart from this file and the lookups (E1).
- **Prompt and legs.** `_scratch/s12/cmp/cmp_dl_review_prompt_p1.md`, sha256 `3588b049…` (E1). Both legs received it
  identically, in parallel.
  - **Codex:** `gpt-6.1-sol`, effort xhigh. Report `/var/projects/_review_scratch/s12-cmpdl-p1-codex/REVIEW.md`.
  - **Grok:** `grok-4.7`, effort high. Report `/var/projects/_review_scratch/s12-cmpdl-p1-grok/grok_stdout.txt`.
- **Verdicts.** Codex: FINDINGS, 3. Grok: FINDINGS, 2 (E2). Four are distinct.
- **Gate.** This is a decision list under G2: one pass, findings verified, folded once, then the builder starts.

## Findings

**K1. A face delta's coefficient is compared as a function, not as its value on the face.** (Codex 1 and Grok 1,
High.) **CONFIRMED.**
- **The list.** Item 3 compares an object's "coefficient of each face's delta" with its counterpart.
- **The gap.** A delta's coefficient is data about the distribution only on the delta's support. The engines print
  both the raw coefficient and its trace on the face (E3). A comparator following item 3 could compare the raw
  coefficients as ordinary functions. It would then report a difference between distributions that are equal on the
  supplied domain. Each leg measured this on the accepted transcripts (E5).
- **What must be true.** Each delta part is compared as the coefficient restricted to that delta's support, with its
  normalization in `w`, derived from the printed distribution itself. Raw coefficients may also be compared, but only
  as a labelled representation diagnostic. A control changes a coefficient away from the face, without changing the
  distribution, and the delta-part comparison does not move.

**K2. No control catches a comparator that only samples inside the supplied domain.** (Grok 2, High.) **CONFIRMED.**
- **The list.** Item 4 prints the residual "as computed" beside its value on the supplied domain. Item 6's controls
  each change a sample inside that domain.
- **The gap.** A factor that equals one on the supplied domain leaves every sample inside it unchanged. A comparator
  that only samples there shows a residual of zero, hides the factor, and still passes every required control (E5).
- **What must be true.** The as-computed residual is exact where feasible. A control changes only a factor that equals
  one on the supplied domain, and the as-computed residual must move.

**K3. The test functions need not separate a function's arguments.** (Codex 2, High.) **CONFIRMED.**
- **The list.** Item 4 asks for "test functions of several different forms" for each live function, `G` included.
- **The gap.** Several forms can all depend on one combination of the arguments. For `G(w − w₋, w₊ − w)`, forms that
  depend only on the sum of the offsets give equal partials. They therefore cannot see a term that depends on the
  partials' difference: Codex's diagnostic drops one such term and finds no residual under six such forms. One
  admissible form with unequal partials does see it (E5).
- **What must be true.** For every live function, the test functions separate each argument's dependence, so that
  derivatives in different arguments differ. That includes both offsets of `G`. Controls fail when the coverage is
  reduced to one point, to one form, or to forms that do not separate the arguments.

**K4. A shared tag name can bind payloads that are not the same object.** (Codex 3, Medium.) **CONFIRMED.**
- **The list.** Item 1 joins every shared tag name. Item 2's table maps symbols and functions only.
- **The gap.** Some shared tags hold differently constituted payloads. A current is a five-component spacetime tuple in
  SymPy and a four-component spatial tuple in Wolfram. `material_elimination_operand` is a bulk multiple in SymPy and
  time-derivative rules in Wolfram (E4, E5). Joining by name alone would let such pairs be compared position by
  position.
- **What must be true.** Before a joined binding is compared as values, the comparator establishes, with citations to
  both engines' sources, which components denote the same object. A binding without such a correspondence is listed as
  not compared, with its reason. A structural mismatch is printed as a comparison that could not be formed. Both
  operands stay printed under their original tags.

## Outcome

Four findings, confirmed. The list as reviewed is preserved, and the fold makes each "What must be true" hold. That is
the one fold G2 allows; the builder starts after it.
