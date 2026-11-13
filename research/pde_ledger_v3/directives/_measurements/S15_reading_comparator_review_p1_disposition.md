# S15 reading comparator decision list, pass 1

The user deferred the independent item-4 check and requested a fresh Codex agent and Grok for the reading comparator. Both reviewed revision 1, preserved in `f683e570` (SHA-256 `96d877eae8c066d7f561fa13228cff87391ff32df7846dabdaf312011cc79ca1`). Codex authored the list. The Codex reviewer started with no inherited conversation and received the list, the source excerpts, and the input inventory. Grok 4.7 received the same list and a review request in an isolated directory, with tools and web search disabled. The user authorized disclosure of the decision list alone. Each reviewer worked independently.

- The fresh Codex reviewer's verdict was “Revision 1 needs a fold before builder handoff for the S15 reading comparison of items 1–3, 5 and 6.” [Codex report](S15_reading_comparator_review_p1_codex.md)
- Grok's verdict was “No blocking findings.” Grok limited that finding to the decision text and explicitly left the cited sources, implementation, and runs unverified. [Grok report](S15_reading_comparator_review_p1_grok.md)

## Findings and dispositions

| Finding | Local check and revision 2 disposition |
| --- | --- |
| Codex 1. Both engines could omit the same required quantity. | Verified that the lookup started from the two printed manifests. Derived a separate roster from the original S15 pieces and groups and the reading specification’s additional sums. All 336 required names appear in each current stream. Removing the same name from both streams’ sets is detectable. Decision 1 requires both directions of this check in the build. |
| Codex 2. Sampling could miss motion or indistinguishable derivative slots. | The earlier comparator specifies multiple function forms, admissible points, distinguishable derivative arguments, a generic moving case, and cases with selected velocities set to zero while higher derivatives stay live. Decision 6 requires that coverage when an exact residual is unresolved and identifies unexercised cases. |
| Codex clarification and Grok 4. Residual interpretation and sign. | Decision 6 specifies the native operands, arithmetic subtraction order, evaluator and translations, and the treatment of typed comparisons. The as-computed residual is printed first. Reductions using a supplied equation get a separate label. The step record interprets the output. |
| Grok 1. A binding must preserve branch and operand identity. | Decision 1 adds native stream, branch, drain, substitution role, derivative slot and evaluation point where relevant. Repeated records carry occurrence offsets. The accepted comparator pairs by source-defined meaning and forbids pairing by position alone. |
| Grok 2. Name the two velocity parts. | The accepted reading list (item 5) and both engine sources use the original quantity `Q` and compute `Q − Q at zero instantaneous velocity of drain k` once for each drain. Grok's suggested split of the two *departures* conflicts with these sources. Decision 5 explicitly names the two drain-indexed parts of `Q`. |
| Grok 3 and 5. Sum coverage, signs, and the open force normalization. | Decision 2 specifies inner, outer, and the signed sum over both surfaces. The measured 309 direct references and 27 new sums must be reconciled against the independently derived roster. Printed measures and signs carry through to the comparison. Force normalization remains open. |
| Grok 6. Excluded tags inside the hashed span. | The Wolfram source writes item-4 records as well as the reading outputs in its base pass. Inputs and scope requires identified exclusions to carry a role and reason. Unrecognized or potentially in-scope records remain unresolved. This preserves the distinction between a deliberate exclusion and a gap in required coverage. |

The [local checks and literal output](S15_reading_comparator_decision_list.md) establish the input hashes, roster, and velocity-part definitions. The [original comparator decisions](../S15_comparator_decision_list.md) give the earlier sampling and record-binding conventions. The [accepted reading list](../S15_reading_decision_list.md) defines the comparison objects.

Both Codex blocking findings are addressed in the single G2 fold. Revision 2 is ready for the comparator builder. The build and its results still require independent review. S15 remains open.
