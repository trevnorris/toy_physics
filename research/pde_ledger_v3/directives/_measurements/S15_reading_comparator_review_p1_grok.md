No blocking findings.

This verdict covers revision 1 of the S15 reading comparator pre-builder decision list in this packet, dated 2026-11-12: scope, evidence rules, and the builder decisions for items 1–3, 5 and 6, including item 1’s same-engine original reference. Item 4 stays deferred. The comparator implementation, engine runs, cited sources, and the other review are outside this verdict. I used only the packet.

The decisions are ready for the builder. The optional wording fixes below can be folded under G2 without reopening the comparison. Readiness still depends on the cited sources containing the definitions this list assigns to them. Those checks are unverified here.

## Blocking findings

None.

## Fidelity, evidence, domains, and branches

The list matches the authorized scope. It compares the SymPy and Wolfram readings accepted at `cc62dc9f` under P5, P6, and `R-S14a-01`, and it uses the originals at `3192123b` for the additional item-1 check (Inputs and scope, packet lines 7 and 9). Singular and the Wolfram–Singular item-4 comparison are deferred. Classification of profiles with inverse-square scaling at every separation is outside this pass. The handoff leaves throat acceleration, the full force sign, source normalization, and the item-4 completeness claim open (packet line 27).

Item coverage inside the list is coherent:

- Item 1 has two checks: cross-engine readings, and each reading against the original in that same engine. The 27 sums are rebuilt per engine from original surface terms, with the original measure and sign, whole pressure in the rest sum, centre and radius terms in transport, and both crossing contributions. Pieces stay in the record, and the compared object is the value of each sum, per Cartesian component (Decisions, item 2, packet line 15).
- Item 2 compares the geometric departure and the response departure separately, with the printed separation derivative, partition operands, residual, and existence conditions. Held quantities stay held. Live field variations and source-centre velocity stay live. Domain, moving surfaces, time, derivative slots, evaluation points, pressure domains, and measures stay with the comparison (item 3, packet line 17).
- Item 3 compares the family integral and the emitted family result on their common stated domain, including printed closed operands, exponent branches, endpoint and special-function conventions, existence and exceptional-status records, and per-engine gaps. Solving all-separation vanishing conditions is a separate problem and is preserved as separate (item 4, packet line 19).
- Items 5 and 6 stay on each object. Each drain is zeroed on its own. The other independent velocities and the higher derivatives remain. Dimensional operands, exponent sets, reference dimensions, differences, and argument checks use the respective dimension bases (item 5, packet line 21).

Reproducibility is specified at decision-list level. Transcripts are checked against accepted hashes and byte boundaries. The Wolfram reading is the prefix through the base-pass marker, because later passes append item-4 work. Each stream is enumerated independently. The input lookup is a cross-reference: 336 named quantities per engine is a located count, and every missing comparison needs a reason (item 1, packet line 13). The comparator at `5a019394` may supply notation and machinery within its documented scope, and its results stay with the original build (packet line 9). Decision 6 requires both native operands and their residual before a verdict, a report that can recover every component and outcome, and reproducible checks listing inputs, coverage, and unresolved objects (packet line 23).

Native-domain and branch handling is specified. Engines keep their native definitions and printed limitations, and accepted native sources are read with their own streams (packet lines 7 and 9). A bare zero and an unevaluated integral of zero keep their native forms, and the comparison is of what each denotes on the stated domain (packet line 17). Exact unevaluated integrals with general fields stay identifiable (packet line 19). One-sided domain restrictions remain in the record (packet line 23). An empty set of dimensional comparisons stays distinct from an evaluated zero vector (packet line 21).

Evidence that can support a verdict is also bounded. Exact comparison is the route where it is feasible. Evaluation, where used, records the common admissible domain, independent inputs for distinct live quantities, test functions, branches, strata, actual points, both values, and the residual. Numerical integration and special functions record precision, error assessment, branch convention, and failures. The outcome classes are exact identity, exact or numerical evaluation at specified points, difference, unresolved evaluation, and missing coverage (packet line 23). A closed-form label, or matching samples at generic exponents, cannot establish equivalence on an exceptional branch (packet line 19). Same-engine byte identity of a payload is evidence for that same-engine check. Decision 6 still applies to every comparison, including the cross-engine comparison of those quantities. Pointwise numerical agreement stays in the point-evaluation class and is not promoted to an exact identity on a domain. This review does not add a numerical tolerance and does not require a code ablation. The one-sided-change check, with other terms kept live, is already an acceptance condition for the builder’s checks and the later build review (packet line 23).

## Optional suggestions

1. Binding identity is shorter than the distinctions later items require. Decision 1, packet line 13, matches by region or surface, component, piece, drain group, and quantity. Branch and exceptional status, reading versus original stream, before versus after substitution, which drain was zeroed, and evaluation or derivative slot are specified later as real distinctions (items 2, 4, and 5). Impact: a matcher built only from the short key can attach a residual to the wrong occurrence. Smallest correction: add those qualifiers to the identity tuple and pair repeats by occurrence within that tuple.

2. “Both velocity parts” is unnamed. Decision 5, packet line 21, says each drain has “the resulting velocity-dependent part,” then asks for dimensional comparison of “both departures and both velocity parts.” Impact: the builder can invent a second velocity split. Smallest correction: state that the two parts are the velocity-dependent parts of the geometric departure and the response departure for that drain.

3. The item-1 partition and the surface union are left to arithmetic and to “over both.” Decision 2, packet line 15. The counts 309 and 27 add to the 336 named quantities, and “over” most naturally means both surfaces. Impact: a builder can treat the lookup counts as a verdict, or can omit the union-of-surfaces sums. Smallest correction: state that 336 = 309 + 27 is a lookup claim to reproduce, that the enumerated manifests govern if the count differs, and that each Cartesian component is checked on each surface and over both surfaces.

4. The residual’s operand order and host are unset. Decision 6, packet line 23, requires the residual before any verdict and cites original comparator practice for the method. Impact: signed differences in the later S15 step record can flip, and a residual taken after import can change a branch cut. Smallest correction: name the order for same-engine and cross-engine residuals, record the host system and any translation, and state that same-engine byte identity is still printed as operands plus residual.

5. Open physical questions can be misread as unfinished preconditions. The handoff, packet line 27, leaves the full force sign and source normalization open, while item 2 requires the original measure and sign. Impact: a builder can rescale or drop an item-1 operand while those questions are open. Smallest correction: state that those open questions stay open and that printed operands, measures, and signs are compared unchanged.

6. Out-of-scope disposition is “by role” without an action. Inputs and scope, packet line 9. Impact: an in-span record can be dropped silently, or an item-4 record inside the hashed prefix can be compared. Smallest correction: record and exclude Singular and item-4 roles, and record any other out-of-scope record inside a hashed span as unresolved coverage.

## Source-dependent checks not verified from this draft

These are limits of this review. They are not defects in the decision text. If a cited source does not contain the assigned rule, that failure belongs to the G2 source fold.

- Hashes, byte boundaries, the literal base-pass marker, and whether that marker byte is inside the Wolfram prefix: `_measurements/S15_reading_comparison_inputs.md`, cited at packet line 9.
- The unique sum formulas behind “whole pressure,” “centre and radius terms,” and “both contributions,” and that they agree with the item-1 prose constraints: same file, cited at packet line 15.
- That 309 and 27 are same-engine facts, partition the 336, and name the same objects in both engines: packet lines 13 and 15, defined only in the lookup.
- In-scope tags, premises P5, P6, and `R-S14a-01`, the separation family, held versus live quantities, and the original wording of items 1–3, 5 and 6: `S15_reading_decision_list.md`, cited at packet lines 7, 17, and 19. Fidelity of this list to those item statements is therefore unconfirmed.
- Native field names for before, after, the velocity-dependent part, exponent sets, reference dimensions, differences, argument checks, and the two dimension bases: the inputs file, cited at packet line 21.
- Roles of out-of-scope records and the residual convention in the original comparator practice: packet lines 9 and 23, via `S15_comparator_decision_list.md`.
- What the repository memory and zero-swap guard require of CAS and numerical work: packet line 27.
- The user-decision text that defers item 4: `../DEFERRED_HEAVY_RUNS.md#s15--independent-item-4-check`, cited at packet line 7. The deferral itself is stated in this list.

## Readiness

A builder who follows these six decisions, verifies the hashed inputs, treats the lookup as a cross-reference, and takes formulas and field names from the cited accepted sources can build the scoped comparator without inventing a scientific comparison. Deferred item 4, the all-separation classification, and the open physical questions stay out of the verdicts. No blocking change is required before that build.
