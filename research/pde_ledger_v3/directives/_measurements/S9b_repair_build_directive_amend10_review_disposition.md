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
