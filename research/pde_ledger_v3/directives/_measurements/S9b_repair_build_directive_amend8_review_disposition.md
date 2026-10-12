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

## Round 1

**Reviewed version.** The working-tree directive, sha256 `8127b59a…`
(`_scratch/s9b_build/s9b_repair_build_directive_amend8_review_baseline_r1.sha256`), frozen as
`_scratch/s9b_build/S9b_repair_build_directive_amend8_reviewed_r1.md`. Item 9 carried the user's two choices: the
`ρ₀ > 0` premise, and no sign predicate on the local density, with the `ρ > 0` and `ρ ≤ 0` parts labelled. Identical
prompt `_scratch/s9b_build/s9b_repair_build_directive_amend8_review_prompt_r1.md` (sha256 `9c29891f…`) to both legs;
both reported before I adjudicated.
- **Codex: CLEAR, no findings.** Final message `_scratch/s9b_build/s9b_amend8_review_r1_codex_final.txt`; full
  report, scripts and stdout copied to `_scratch/s9b_build/s9b_amend8_review_r1_codex_evidence/`.
- **Grok: NEEDS REVISION, one finding.** Report `_scratch/s9b_build/s9b_amend8_review_r1_grok.txt` (from
  "**Verdict: NEEDS REVISION**"); scripts and stdout copied to `_scratch/s9b_build/s9b_amend8_review_r1_grok_evidence/`.

**Verification.** Mechanical lookups in the "Round 1" section of
`S9b_repair_build_directive_amend8_review_disposition_lookups.md`, generated by
`_scratch/s9b_build/gen/s9b_amend8_r1_lookups.sh`.

| # | Finding (leg) | Disposition | Round-2 change |
|---|---|---|---|
| G3 | On the `ρ ≤ 0` part, "the response being real and positive" does not pick out one set, because the text does not fix which expression is the response. The spec's fixed-ratio relation, the non-negative root of the ratio of sound-speed squares, and the identity `(ρ/ρ₀)^((n−1)/2)`, which holds only for positive density, differ on an open set of negative densities: one is positive where the other is negative. At zero density, the engines' evaluators read a zero exponent (`n = 1` in `c_s²`, `s = 0` in the power response) differently, SymPy as `0⁰ = 1` and Wolfram as indeterminate, so the two engines attach different domains. Fix: the fixed-ratio response is the non-negative root form, not the identity; a zero exponent reduces before the base is substituted. (Grok 1) | **ACCEPT.** Lookups: the spec writes the fixed-ratio relation as `c_γ/c_s = c₀/c_s0` (L245), `c_s² = nKρ^(n−1)/m` (L219–220), and `c_γ` as a non-negative root (L9–10). Both engines currently build the root form (PY L985–986, WL L74–75), so the open-set split is a risk the text leaves open, not a present difference. Grok's scripts print the two forms differing at negative density, and the two evaluators differing at zero density with a zero exponent (SymPy `1`; Wolfram `Indeterminate`). Keeping these results was the user's choice; a rewriting valid only for positive density would delete them again. Codex found no issue. | Item 9: on both parts, each row's response is the spec's relation as written, each sound speed the non-negative root of the spec's `c_s²`; a rewriting that holds only where `ρ > 0` is not used for the `ρ ≤ 0` part. A zero power of a density or density ratio is `1` at every value, zero included. |

**Process.** Grok left one Wolfram run stuck at its 32-task ceiling, reported it, and used only the lines it printed
before. The leg had exited and the owner was dead; I stopped the unit and reconciled the registry (record in the run's
`RECONCILIATION.md`; backup `reservations.json.stale-amend8-r1-grok.bak`).

**Round 1 outcome.** Item 9 is revised for round 2. Not clear yet.

## Round 2

**Reviewed version.** The working-tree directive, sha256 `8218600d…`
(`_scratch/s9b_build/s9b_repair_build_directive_amend8_review_baseline_r2.sha256`), frozen as
`_scratch/s9b_build/S9b_repair_build_directive_amend8_reviewed_r2.md`. Item 9 gained a bullet naming each response on
both parts and banning, on the `ρ ≤ 0` part, a rewriting that holds only where `ρ > 0`. Identical prompt
`_scratch/s9b_build/s9b_repair_build_directive_amend8_review_prompt_r2.md` (sha256 `c200e75c…`; it adds the
`--tasks-max 64` line for Wolfram runs) to both legs; both reported before I adjudicated.
- **Codex: CLEAR, no findings.** Final message `_scratch/s9b_build/s9b_amend8_review_r2_codex_final.txt`; full
  report, scripts and stdout copied to `_scratch/s9b_build/s9b_amend8_review_r2_codex_evidence/`.
- **Grok: NEEDS REVISION, two findings.** Report `_scratch/s9b_build/s9b_amend8_review_r2_grok.txt` (from
  "**Verdict: NEEDS REVISION**"); scripts and stdout copied to `_scratch/s9b_build/s9b_amend8_review_r2_grok_evidence/`.

**Verification.** Mechanical lookups in the "Round 2" section of
`S9b_repair_build_directive_amend8_review_disposition_lookups.md`, generated by
`_scratch/s9b_build/gen/s9b_amend8_r2_lookups.sh`.

| # | Finding (leg) | Disposition | Round-3 change |
|---|---|---|---|
| G4 | The ban covered only the `ρ ≤ 0` part, so on `ρ > 0` the fixed-ratio response could still be the identity `(ρ/ρ₀)^((n−1)/2)`. That is a different object: it carries no sound speed, so K9a's knife does not reach it; it is positive at `n = 0`, where the root ratio is `0/0`; and at `n = −1`, where the squares are negative, a principal root makes the ratio real and positive. (Grok 1) | **ACCEPT.** Lookups: Grok's stdout prints, at `n = 0`, `ρ = 4`, the root ratio `ZERO_OVER_ZERO` and the identity `1/2`; at `n = −1`, `ρ = 4`, no non-negative root, and the principal ratio and the identity both `1/4`. The scope "for the `ρ ≤ 0` part" was mine; the response is one object on both parts. | The fixed-ratio response is one object on both parts: the spec's relation solved for `c_γ`, with each sound speed the non-negative root of its square. It has a value only where both squares are real and non-negative and the asymptotic root is positive; a principal root of a negative square is not a value of it. No expression that differs from a response at some point of either part stands in for it. |
| G5 | Where the local root vanishes and the asymptotic one does not (`ρ = 0` at `n = 2` or `3`), the spec's quotient `c_γ/c_s = c₀/c_s0` has no value, while the bullet's form `c_γ/c₀ = c_s/c_s0` gives `0`. With item 1 making the spec win every conflict, the labelled `ρ ≤ 0` part could print `0` or "undefined". Positive-set membership agrees either way. (Grok 2) | **ACCEPT.** Lookups: Grok's stdout prints the root ratio `0` and `SPEC_QUOTIENT_UNDEFINED_CS_ZERO` at `n = 2`, `ρ = 0`; item 1 (L63–65) makes the spec win. Keeping a result there is the user's choice; the value is the relation solved for `c_γ`. | Item 9 supplies the value `0` where the local root vanishes and the asymptotic root is positive, so it is not a stop event under item 12 and item 1 does not replace it. |

**Premise check.** This is the amendment's third round. Each round has narrowed one object, which the user's choice
to keep non-positive densities brought into the build: the fixed-ratio response over every density and exponent. Round
3 names it whole instead of adding another scoped rule. Codex has found nothing in rounds 1 and 2; Grok's findings are
each backed by both CAS.

**Round 2 outcome.** Item 9's response bullet is rewritten for round 3. Not clear yet.
