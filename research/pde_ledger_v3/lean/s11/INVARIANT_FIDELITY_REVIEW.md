# D2 S11 invariant fidelity review and closure

**I1–I4 are complete under the bounded contract in
[INVARIANT_NEXT.md](INVARIANT_NEXT.md).** Both independent non-author reviews
returned CLEAR with no blocking findings. Closure recorded 2026-09-16T12:27:34.421380+00:00.
This does not complete the whole S11 step or downstream CAS production work.

## Reviewed object and evidence

The fixed 28-file packet (plus its manifest) has aggregate SHA-256
`f95d1667e54b4c8139881feaa86267ebe60eb4e9d85f0f4f4bd2dff9dc2b2747`. The aggregate hashes the file-to-SHA256 map as UTF-8
`json.dumps(files, sort_keys=True)`. The
[manifest](../../_measurements/S11_lean_invariant_review_packet.json),
[archive](../../_measurements/S11_lean_invariant_review_packet_v1.tar.gz),
[local verification](INVARIANT_VERIFICATION.txt), and
[closure validation](../../_measurements/S11_lean_invariant_closure_validation.json)
retain the exact source and result identities. The user authorized transfer
of this fixed packet by saying “Approved”.

The result classifies all real quadratic forms on arbitrary 2×2 matrices under
all SO(2)/O(2) conjugations. The spaces have dimensions 4 and 3; the reflection-odd
subspace within SO invariants has dimension 1. The even/odd decomposition is
exhaustive and unique. The theorem does not take candidate-basis completeness
as an assumption. The compact native correspondence compares full D2 V1/V2/V6
spans after the minimal per-generator coefficient-action transpose correction.

Local verification passed five modules and the audit root, 46 standard-axiom
audits, nine mathematical mutations, and eight explicit positive controls.
No proof, native engine source, verification instrument, or terminal check
result changed after the reviewed packet was frozen. The only subsequent
changes to packet members are closure documentation.

## Independent reviewers

| Reviewer | Independent session | Verdict | Durable report |
|---|---|---|---|
| Claude, Opus 5 | `97d2cc9d-20df-4063-b08e-f90a22a9a340` | CLEAR; no blockers | [Claude report](../../_measurements/S11_lean_invariant_fidelity_claude_v1.md) |
| Grok, 4.6 | `b6a98ca4-8b3f-477c-8b89-b27f2a7617b6` | CLEAR; no blockers | [Grok report](../../_measurements/S11_lean_invariant_fidelity_grok_v1.md) |

The sessions ran sequentially against the same isolated snapshot, without
access to the other review in that packet. Both terminated with exit status 0
and `end_turn`, with substantive final reports identifying the correct packet.
The raw JSON, process records and stderr remain unchanged. The Grok markdown
report starts at its terminal verdict; its preceding progress prose is not
part of the extracted report. Model-usage metadata remains in the raw records.

These are source fidelity reviews. The reviewers did not independently rebuild
Lean or recompute cryptographic hashes; their hash discussion cross-checks the
provided records. The local validation separately recomputed the live,
snapshot, archive and transport hashes, checked all six compiled modules, and
rechecked all 15 dependency checkout pins with no tracked-source changes.
Grok's stderr contains nonfatal plugin/hook and `/tmp` repository-discovery
warnings; they did not prevent a completed review or supply its verdict.

## Finding dispositions

No substantive change or additional proof was required. Optional suggestions
were handled within the policy's stopping rule:

| Finding | Disposition |
|---|---|
| Claude: the source-check residual uses an explicit polynomial, rather than extracting it from the in-memory mutant basis. | Clarified in INVARIANT_FIDELITY.md: this residual is an independent polynomial check. Its provenance as a pre-repair emitted polynomial is established separately by the frozen native probe; the exact same polynomial and witness are proved in Controls.lean. The native full-span comparison is a separate check. No stronger linkage is claimed for that instrument line. |
| Both: the orientation mutant also emits an unused-simp-argument linter error. | Retained the inspected evidence. Acceptance specifically requires the false algebraic identity in `coefficient_action`; the linter is not the rejection criterion. Cosmetic isolation is optional and does not warrant a new proof/control cycle. |
| Claude: distinguish the missing-coefficient mutation from the actual completeness claim. | Clarified in INVARIANT_VERIFICATION.txt. `odd_pairing_coefficient` checks the definition/evaluation identity; `odd_pairing_is_required` rejects the claim that the odd pairing lies in the three-form even span. |
| Claude: the automated Wolfram guard is only a substring check. | Clarified that limit. The two reviewers inspected the actual per-generator `Table` context and found the transpose in the correct location. No fresh Wolfram execution or stronger automated source check is claimed. |
| Claude: transfer permission was still marked pending in the frozen documentation. | Updated current closure documentation and the authorization/review state. The reviewed snapshot remains unchanged as historical evidence. |
| Grok: SymPy does not transpose the reflection-difference block. | No defect for the specified diagonal reflection: that block is diagonal in the monomial basis. This was independently checked by both reviewers; D2 native span equality also passes. |
| Grok: original probe status still reports a source discrepancy. | Preserved intentionally against its pre-repair source hash. Current documentation identifies it as historical and directs current checks to the repaired-source instrument. |

Claude's report has two incidental prose count errors: the manifest contains
28 payload files plus MANIFEST, and the suite has seven concrete false-statement
mutants plus two source mutants, rather than six concrete controls. The exact
manifest and records, including all nine mathematical failures, were validated
locally. These prose errors change neither the mathematical review nor the
terminal verdict. Reviewer reports are preserved verbatim.

## Live correspondence, preservation and stopping rule

At closure validation, all 28 live packet files matched the reviewed snapshot;
there were **no shared-source changes to reconcile**. The current SymPy and
Wolfram hashes equal the reviewed hashes. All proof/dependency and compiled
module hashes match the validated local records. Documentation changes for
closure are enumerated in the final validation record.

H1–H4 remains its separately completed historical contract at checkpoint
`d80d02d2`; no homogeneous proof or historical check record was rewritten.
This closure made no writes to S11c calculations or to the pinned S11c-b/c1/c2
exports requested to remain unchanged for the other running job. Other-session
and unrelated changes were preserved. No commit is requested or made here.

D3–D5 completeness, V5 Euler–Lagrange/total-divergence classification, PD-package
production consequences, comparator/export reconciliation, and S11b/S11c work
remain outside I1–I4. The other calculation session's Q9 investigation remains
separate; it is not cleared by these reviews. The bounded Lean deliverable is
complete, so no additional coverage or CAS bridge expansion is added.
