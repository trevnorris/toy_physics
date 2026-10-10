# S9b repair build directive, amendment 6: review dispositions (orchestrator)

**Artifact.** Amendment 6 to `directives/S9b_repair_build_directive.md`, written by the orchestrator for the build's
repair round 3 (`S9b_repair_build_r2_review_disposition.md`, findings C1 and C2). It changes what the engines and
harnesses print, so it is reviewed as physics until clear, by Codex and Grok (`CLAUDE.md` G1, orchestrator-written).

## Round 0

**Reviewed version.** The working-tree directive, sha256 `46b942e4…`
(`_scratch/s9b_build/s9b_repair_build_directive_amend6_review_baseline_r0.sha256`), frozen as
`_scratch/s9b_build/S9b_repair_build_directive_amend6_reviewed_r0.md`. Identical prompt
`_scratch/s9b_build/s9b_repair_build_directive_amend6_review_prompt_r0.md` to both legs; both reported before I
adjudicated.
- **Codex (`gpt-6.1-sol`, xhigh): CLEAR.** Final message `_scratch/s9b_build/s9b_amend6_review_r0_codex_final.txt`;
  full report and scripts copied to `_scratch/s9b_build/s9b_amend6_review_r0_codex_evidence/`.
- **Grok (`grok-4.7`): NEEDS REVISION, two findings.** Report `_scratch/s9b_build/s9b_amend6_review_r0_grok.txt`;
  scripts and stdout copied to `_scratch/s9b_build/s9b_amend6_review_r0_grok_evidence/`.

**Verification.** Mechanical lookups in `S9b_repair_build_directive_amend6_review_disposition_lookups.md`, generated
by `_scratch/s9b_build/gen/s9b_amend6_r0_lookups.sh`.

| # | Finding (Grok) | Disposition | What must be true |
|---|---|---|---|
| 1 | Item 14's new bullet does not say which part of a payload is the value when the computed content is itself logical. Wolfram's `A_NONRECIPROCAL_PATH_DEPENDENCE` carries a `Closedness` predicate beside its `Domain`; the `γ` objects are a relation under a gate. A harness could put every logical part on the domain side, print zero for the value, and meet the bullet's words. The bullet also does not cover K11's compact differences and residual explicitly. | **ACCEPT.** Lookups: Wolfram audit L166–171 emits `"Closedness" -> …` and `"Domain" -> …` as sibling keys; the SymPy harness differences a bare boolean as one `Xor` (PA L213–214). Codex's probe tested a relation paired with a domain (`payload_review.wl` L37), not a bare logical value, so it does not bear on this case. Codex's CLEAR is weighed against a gap the text leaves open. | The bullet defines the domain as the condition a payload attaches to its computed object, and the value as the computed object itself, including a relation, a predicate or a condition. It covers K11's compact differences and method residual, and K11's sample domain is coverage. |
| 2 | The Part 2 line says "items 5 and 6 hold" for the `γ` objects "inside the branch domain". "The branch domain" is undefined, and item 6's stratum sentence covers only the `γ` difference. Dropping the `GM ≠ 0` gate while keeping a solved quotient would meet the line as written but not C1. | **ACCEPT. The cause is mine:** C1's required outcome says each `γ` object "prints … what its relation determines there"; my Part 2 line compressed that into a pointer to items 5 and 6, which lost it. Lookups: "branch domain" occurs once (L310); item 6 L135–136; the r2 disposition's C1 outcome; spec L193–198 states two conditions. Grok's computed values on the `GM = 0` stratum are reviewer measurements and stay out of the directive. | The line defines the branch domain as the set where both conditions of the spec's "Branch existence" paragraph hold, names the `γ` objects, and states C1's outcome: on every stratum inside the branch domain, including `GM = 0`, each object and its restriction and forward copies print what its relation determines there, computed by the engine; `NOT_ESTABLISHED` only outside. |

**Revision.** Both repairs are made in the working tree; the reviewed r0 bytes are kept in the frozen copy. Round 1
goes to fresh Codex and Grok legs on an identical prompt.

**Process.** Codex left an auxiliary Wolfram probe running when it exited (unit `s11c-guard-f4020c055d31`), with its
launcher dead. I stopped it and reconciled the stale reservation by the guard's procedure; the record is
`_scratch/s11c/s9b-amend6-review-r0/codex-payload-wl-01/RECONCILIATION.md`.
