# S9b repair build directive, amendment 10: review dispositions (orchestrator)

Amendment 10 answers the build review r6 (`S9b_repair_build_r6_review_disposition.md`, rows C1 and C2). It changes what
two ablation harnesses print, so it is reviewed as physics until clear. The orchestrator wrote it, so the legs are
Codex and Grok (`CLAUDE.md` G1). Lookups: `S9b_repair_build_directive_amend10_review_disposition_lookups.md`
(mechanical; generator `_scratch/s9b_build/gen/s9b_amend10_r0_lookups.sh`).

## Round 0

**Artifact.** The working-tree directive at sha256 `8dbf6025…`. Its exact bytes are preserved in
`_scratch/s9b_build/preserved/S9b_repair_build_directive_amend10_r0.md`. Prompt:
`_scratch/s9b_build/s9b_repair_build_directive_amend10_review_prompt_r0.md`.

**Legs.**
- **Codex (gpt-6.1-sol, xhigh): NEEDS REVISION, three findings.** Final message:
  `_scratch/s9b_build/s9b_amend10_review_r0_codex_final.txt`. Full report and scripts:
  `_scratch/s9b_build/s9b_amend10_review_r0_codex_evidence/`.
- **Grok: NEEDS REVISION, two findings.** Report: `_scratch/s9b_build/s9b_amend10_review_r0_grok.txt`. Scripts and
  stdout: `_scratch/s9b_build/s9b_amend10_review_r0_grok_evidence/`.

Both legs found:
- no leakage;
- the engines unchanged;
- that printing the witness data is consistent with "It prints no verdict".

| # | Finding (leg) | Disposition | Round-1 text |
|---|---|---|---|
| A1 | The zero, change and undecided cases are not exclusive, and nothing ranks them. An equal part could print as a change if "differ exactly" compares substituted forms: the inert integral of zero is not syntactically `0`, but is equal after `Activate` or `doit`. A part with a symbolic change certificate could still print undecided after a failed witness search. (Codex 3, Grok 2) | **ACCEPT.** Lookups: in round 0, "has looked for such a point and profile and found none" was the only condition for undecided, and the "so an unchanged part prints as zero" sentence did not rank the cases. Grok's guarded run prints `substituted_sameq_zero False` and `substituted_activate_sameq_zero True`. | A numbered rule, in order. Zero where the engine shows the part unchanged, including an inert integral of zero against `0`; no evaluation overrides it. Otherwise a change, where the engine shows one symbolically or by a witness. Otherwise undecided, once the harness has looked for a witness and found none. |
| A2 | "By more than their printed error bound" does not say whose bound. Two equal operands, each evaluated within a bound of 1, can differ by more than either single bound. (Codex 1; Grok 2, the machine value of `Sin²+Cos²−1`) | **ACCEPT.** | A value difference counts when the two evaluations differ exactly, or by more than a bound that covers the error of both. The bound prints when the arithmetic is not exact. |
| A3 | "Admissible profile" is not one set of functions. Requiring the profile to meet "the part's domain" can exclude the witness a domain change needs: no profile meets both a condition and its negation. Item 13's "general radial function" and "no profile family" sentences still read against a concrete witness profile. (Codex 2, Grok 1) | **ACCEPT.** Lookups: the round-0 definition and item 13's family sentence. Grok's guarded run prints `exists_profile_meeting_both False`. | An admissible profile assigns one specific radial function, with fixed coefficients, to each live profile the part carries. It is defined, with every derivative the part takes, at the witness point. It need not meet the condition being compared, nor the spec's untested order counting. A witness's point and profile are neither a profile choice nor a profile family under item 13. |
| A4 | A condition inside a quantifier can change while both closed statements evaluate the same. Grok's case: both are `False`, while their matrices at `r = 3/2` are `True` and `False`. Grok asks for the witness point to give values to the bound variables. (Grok 2, part) | **ACCEPT IN PART.** The ranking fix (A1) is accepted. The bound-variable instantiation is not. Item 14 calls a part unchanged when it is "equal to the baseline part as a mathematical object". Two quantified conditions that hold at the same points of their free variables define the same domain, so they are unchanged even when their matrices differ. Instantiating bound variables would print a change for an unchanged part, which is A1's defect in the other direction. A change in a quantified condition that shows at some free point, under some admissible profile, is found by the witness search, which A1 now requires before undecided. | The text states that a condition binding variables, in a set or a quantifier, is compared as the condition it states on its free variables, so its witness point gives values to those only. |

**Process.** Grok's first retry loop called the guard with its arguments out of order, so the guard admitted nothing.
Grok saw that and relaunched correctly. Its Wolfram run exited 0. Its SymPy run exited 1 because SymPy 1.14 has no
`ForAll`, which Grok reports. Codex's guarded checks completed.

## Round 1

**Artifact.** The working-tree directive at sha256 `5da3c7e2…`. Its exact bytes are preserved in
`_scratch/s9b_build/preserved/S9b_repair_build_directive_amend10_r1.md`. Prompt:
`_scratch/s9b_build/s9b_repair_build_directive_amend10_review_prompt_r1.md`.

**Legs.**
- **Codex (gpt-6.1-sol, xhigh): CLEAR, no findings.** Final message:
  `_scratch/s9b_build/s9b_amend10_review_r1_codex_final.txt`. Full report and scripts:
  `_scratch/s9b_build/s9b_amend10_review_r1_codex_evidence/` (`report.txt`).
- **Grok: NEEDS REVISION, one finding.** Report: `_scratch/s9b_build/s9b_amend10_review_r1_grok.txt`. Scripts and
  stdout: `_scratch/s9b_build/s9b_amend10_review_r1_grok_evidence/`.

| # | Finding (leg) | Disposition | Round-2 text |
|---|---|---|---|
| B1 | A changed set prints undecided. A set is equal to another when the two have the same members, so the point at which two sets differ is an element, which is a value of the coordinate the set binds. The round-1 sentence lets a point for "a condition that binds variables, in a set or a quantifier" give values to free variables only. Once the free variables are fixed, a builder has no coordinate left to assign, and a set change that membership shows prints undecided. Grok keeps the free-variable rule for quantifiers. The SymPy engine builds set-valued objects (`roots`, `traversed`). (Grok 1) | **ACCEPT.** The two legs agree on the mathematics and differ on the text. Codex's run of the round-6 harness prints `UNDECIDED` for the pair `ConditionSet(r, r > a, Reals)` and `ConditionSet(r, r > a + 1, Reals)`, and its own membership check at `3/2` after `a = 1` prints `(True, False)`. Codex counts that check as a showing. Grok reads the round-1 sentence as excluding it. A new author could read it Grok's way, so the text must say it. Lookups (`S9b_repair_build_directive_amend10_review_r1_lookups.md`, generator `_scratch/s9b_build/gen/s9b_amend10_r1_lookups.sh`): the round-1 sentence, Codex's `R6_SET_DIFFERENCE` and membership lines, Grok's `set_decision.py.stdout`, and the engine lines. | Two sets differ at an element that one contains and the other does not, so a point for a set gives values to its free variables and to that element. A statement that binds variables, such as a quantifier or a test of whether a set is empty, is compared as the statement it makes on its free variables, so its point gives values to those only. |

**Author.** B1 is in the sentence round 1 added, the first defect in material this amendment just changed. Under G4,
a further defect in material a fold just changed moves the text to a new author.

**Process.** Grok's guarded run directories date from 23:51 to 00:00, except `sat-bound-py-2`, from 01:45 (lookups).
Codex retried an admission refusal until admitted, and its guarded checks completed.

## Round 2

**Artifact.** The working-tree directive at sha256 `faf60e3c…`, the accepted text. Prompt:
`_scratch/s9b_build/s9b_repair_build_directive_amend10_review_prompt_r2.md`. Lookups:
`S9b_repair_build_directive_amend10_review_r2_lookups.md` (generator `_scratch/s9b_build/gen/s9b_amend10_r2_lookups.sh`).

**Legs.**
- **Codex (gpt-6.1-sol, xhigh): CLEAR, no findings.** Final message:
  `_scratch/s9b_build/s9b_amend10_review_r2_codex_final.txt`. Full report and scripts:
  `_scratch/s9b_build/s9b_amend10_review_r2_codex_evidence/` (`report.txt`).
- **Grok: CLEAR, no findings.** Report: `_scratch/s9b_build/s9b_amend10_review_r2_grok.txt`. Scripts and stdout:
  `_scratch/s9b_build/s9b_amend10_review_r2_grok_evidence/`. Grok checked B1's sentence on the set-valued and
  quantified objects both engines build, including sets nested in another set's condition, and the A4 rejection
  against the same text.

**Disposition.** Nothing outstanding changes what the harnesses compute or print. Amendment 10 is accepted.

**Process.** The round-0 and round-1 preservation commits (`331e748d`, `192adc31`, both 06:16:01) put those texts at
the directive's path for about a second, after which the round-2 text was restored. Codex was then in a guard
admission retry loop; its final baseline check (06:47) reads `OK`. Grok's directive quotes date from 06:29 and carry
the round-2 sentence.
