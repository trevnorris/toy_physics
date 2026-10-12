# S9b repair build directive, amendment 8: review dispositions (orchestrator)

**Artifact.** Amendment 8 to `directives/S9b_repair_build_directive.md`, written by the orchestrator for the build's
repair round 5 (`S9b_repair_build_r4_review_disposition.md`, findings G1, C1 and C2). It changes what the engines and
their harnesses print, so it is reviewed as physics until clear, by Codex and Grok (`CLAUDE.md` G1,
orchestrator-written).

## Round 0

**Reviewed version.** The working-tree directive, sha256 `cf258688…`
(`_scratch/s9b_build/s9b_repair_build_directive_amend8_review_baseline_r0.sha256`), frozen as
`_scratch/s9b_build/S9b_repair_build_directive_amend8_reviewed_r0.md`. Identical prompt
`_scratch/s9b_build/s9b_repair_build_directive_amend8_review_prompt_r0.md` (sha256 `4fba4aca…`) to both legs; both
reported before I adjudicated. Item 9 in this version made the local bulk density's positivity a supplied-input domain
predicate of each row whose response reads it.
- **Codex (`gpt-6.1-sol`, xhigh): NEEDS REVISION, one finding.** Final message
  `_scratch/s9b_build/s9b_amend8_review_r0_codex_final.txt`; full report, scripts and stdout copied to
  `_scratch/s9b_build/s9b_amend8_review_r0_codex_evidence/`.
- **Grok (`grok-4.7`): NEEDS REVISION, two findings.** Report `_scratch/s9b_build/s9b_amend8_review_r0_grok.txt` (the
  report starts at "**Verdict: NEEDS REVISION**"; the text before it is progress narration); scripts and stdout
  copied to `_scratch/s9b_build/s9b_amend8_review_r0_grok_evidence/`.

**Verification.** Mechanical lookups in `S9b_repair_build_directive_amend8_review_disposition_lookups.md`, generated
by `_scratch/s9b_build/gen/s9b_amend8_r0_lookups.sh`.

**The user's direction during the round.** While the legs ran, the user read my report of the round-0 text and said:
"Hold up on deleting results. What if they are significant in ways we don't understand yet." The user then chose to
keep the results where the local bulk density is zero or negative, labelled, and to state the asymptotic bulk
density's positivity as a premise; both engines already assume it (`sp.symbols('rho_0 m', positive=True)`, and
Wolfram's `$Assumptions` with `rho0 > 0`), and the spec does not write it. The user also chose to let the round-0
legs finish. Both legs reached the same conclusion on the round-0 text independently.

| # | Finding (leg) | Disposition | Round-1 change |
|---|---|---|---|
| X1 | Strict positivity of `ρ₀(1+f)` is a premise the spec does not supply. A nonnegative density with a local zero leaves the fixed-ratio response defined and positive at `n = 1`, and the amendment excludes that point. (Codex 1) | **ACCEPT.** The spec writes no inequality on either density (Grok's quotes agree). | Item 9: the local density's sign is not a domain predicate, so no row excludes a point where it is zero or negative. Zero density falls in the labelled `ρ ≤ 0` part. |
| G1 | Item 9 states a sign the spec does not write, without item 10's explicit item-12 status, so a builder can stop or conjoin; at negative and zero density both responses are defined and positive, and `ρ₀(1+f) > 0` does not say both densities are positive. Fix: no engine conjoins a sign for `ρ₀(1+f)`; a sign belongs in the spec first, or with explicit item-12 status. (Grok 1) | **ACCEPT.** Lookups: Grok's script prints the responses defined at negative and zero density, and the sign predicates false there. | As X1. The one premise kept, `ρ₀ > 0`, is the user's, supplied in item 9 with explicit item-12 status, and printed by each engine as a supplied premise. |
| G2 | "Whose response reads `ρ`" does not fix the K9a–b copies: the knives rebuild the responses in forms that differ between engines, so whether the conjunct survives a knife is not one shared object. Stands only while the sign is kept. (Grok 2) | **ACCEPT the underlying point; the conjunct is gone.** Codex's report adds that a rule keyed to a row's name does not show that a corrupted response reads `ρ`; with no conjunct, nothing is attached on that basis. | The labelled split belongs to the row: the fixed-ratio and power rows, with their implied-`j_n` rows and item 10 copies, print the two parts, and each knife's corrupted copy prints the same two parts. |

**Item 14.** Neither leg has a finding on it. Codex: "No additional findings on item 14, coverage, or M2." Grok's
checks find that the new sentences match C1 and C2, decide the integral-of-zero, sign-flip, set and quantifier cases,
and keep value and domain separate; it notes that an undecided part records no reach. Item 14 is unchanged in round 1.

**Process.** My attempt to stop both legs after the user's direction was refused by the session's permission check,
and the user then chose to let them finish. Each leg left one supplementary Wolfram run stuck at its 32-task ceiling,
reported it, and did not use it as evidence. Both legs had exited and both owners were dead; I stopped the two units
and reconciled the registry (records in each run's `RECONCILIATION.md`; backup
`reservations.json.stale-amend8-r0-legs.bak`).

**Round 0 outcome.** Item 9 is revised for round 1 as above; the header names the user's choices. Not clear yet.
