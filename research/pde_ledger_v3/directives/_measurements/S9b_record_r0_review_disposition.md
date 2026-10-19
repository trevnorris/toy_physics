# S9b record, review round 0: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S9b_record_r0_review_lookups.md` beside this file (generator
`_scratch/s9b_build/gen/s9b_record_r0_lookups.sh`).

## The round
- **Artifacts.** The S9b step record `steps/S9b_light_bending.md`, its measurements file
  `steps/_measurements/S9b_record_measurements.md` and generator `scripts/S9b_record_measurements.py`, and the S9b
  pass in `SUBSTRATE_REQUIREMENTS.md` (against `caa0710d`). Four files in the baseline
  `_scratch/s9b_build/record/s9b_record_review_baseline_r0.sha256`, all OK at review.
- **Contract.** The record decision list `directives/S9b_record_directive.md` at `caa0710d`.
- **Author.** A fresh Codex session (`gpt-6-astra`, high effort) in the main checkout, brief
  `_scratch/s9b_build/record/brief_record_author.md`. Its STOP report is saved beside the brief.
- **Prompt and legs.** The prompt is `_scratch/s9b_build/record/s9b_record_review_prompt_r0.md` (sha256
  `3e851de0…`). It hands the sources the decision list names and no build directive. The legs were a fresh Claude
  agent (opus) and Grok (`grok-4.7`, effort high), launched identically, with scratch under
  `/var/projects/_review_scratch/s9b-record-r0-{claude,grok}/`. Reports:
  `_scratch/s9b_build/record/s9b_record_review_r0_claude.md` and `…_r0_grok.txt`. Both returned **NEEDS REVISION**.

## Findings

| # | Leg | Finding | Disposition |
|---|---|---|---|
| 1 | Claude 1 | **The within-engine deflection/radar statements claim more than the Wolfram payloads show, and no retrieval supports them.** The record says each engine prints "the same profile condition for deflection and radar slope" (lines 7–8), counts no independent profile test because "the within-engine conditions do not differ" (192–193), and says the copies "have the same native condition fields" (267); the register repeats it (824, 849). In SymPy the deflection and radar payloads are byte-identical. In Wolfram all 23 pairs differ in `InverseWeight`, `ReferenceInverse` and `VelocityReduction`; only the `ProfileEquation` and `GM` leaves are identical. The generator has no deflection/radar pairing retrieval. | **Confirmed.** The leg's pair census shows 23 differing Wolfram pairs. The generator's one `radar` match is an exclusion list, so no cited command supports the statements. Grok found the same Wolfram difference and read the record's claim as being about the profile equation only. But line 267's "same native condition fields" is false for Wolfram as written, and every one of these statements is a measured claim without its command (E1). Whether the differing Wolfram leaves state the same condition is a question beyond retrieval. The "no independent profile test" conclusion, which the register hands on, rests on it. |
| 2 | Claude 2 | **The build's K11 claim limit is not reproduced.** The round-7 build disposition (lines 72–74) states: "Neither reaches the shared `A_NONRECIPROCAL_*` tags. Claim limit: K11 supports claims about the local one-form and closedness construction only, not item 8's shared observables." Round 11 left K11 unchanged. The record reports the nonreciprocal path-dependence predicate (lines 96–98, 293) and carries only "K11 unreached parts remain 'not evaluated'" (370). | **Confirmed.** The decision list's item 3 requires the record to reproduce the limits named in the build's acceptance dispositions. The limit is absent (literal count 0), and the round-11 disposition (line 73) keeps K11 unchanged. The record must also stay accurate about what does reach the shared tag: the round-11 Wolfram K10 witness is at `WL_S9B_A_NONRECIPROCAL_PATH_DEPENDENCE`'s domain part (r11 disposition line 65). |
| 3 | Claude 3 | **The register's second route names an incomplete set of prior art.** The S9b pass says "the identified prior-art objects in this pass are the supplied GR/PPN comparison references" (register 845; record 452). The repository also identifies river/acoustic-flow prior art for the flow mechanism: `docs/medium_requirements_and_prior_art.md:246`, and the register's own O2 pass (788–790, "river/acoustic-flow candidate oracles"). The register's method says to run the second route "for every sector with identified prior art". | **Confirmed.** The lookups show both sources and no mention of the river/acoustic class in the record. M3 makes prior art an oracle. The record may not apply it (no new computation), but it must not hand on an incomplete list of the identified oracles, and it must say which were not applied. Whether to run that oracle check against the flow-only condition is new computation, and it goes to the user as a separate decision. |
| 4 | Claude 4 | **The census command reads a moving `HEAD`.** Generator line 179 runs `git show HEAD:…`, and the filed output prints `HEAD ENTRY_COUNT 16` and `WORKTREE ENTRY_COUNT 17` (lines 21288, 21305), which the register cites as "16 → 17 entries" (811). After the S9b commit the same command prints 17 → 17. Line 171 already pins explicit revisions; M1 reads its sources from the worktree. | **Confirmed.** `HEAD` is `caa0710d` now. Once the record is committed, the cited command no longer supports the census and the regenerated file no longer matches the filed one. |
| 5 | Grok 1 | **`R-S21-01`'s requirement applies the full mass-law exchange to every condition.** The requirement says each condition includes its printed implied `jn`, with the coordinate mass expression and both signed-velocity components (lines 558–561), while speed-, tilt- and bulk-only print `jn = 0`. | **Not confirmed.** The sentence follows the requirement's list of conditions: Part B's condition for each observable and Part C's constant-speed, fixed-ratio and power conditions (lines 543–556). In all of those, `V` is live. The single-mechanism restrictions are not in that list. The same entry's evidence clause states their exchange: "Speed-, tilt- and bulk-only exchange prints zero; flow-only retains the signed live-density exchange" (574–575). So nothing handed on changes. |
| 6 | Grok 2 | **`F_MASS_SOLUTION`'s printed residual is attributed to both engines.** The record says "`F_MASS_SOLUTION` prints `V=Φ/(4πr²ρbr)` and a mass residual 0" (248). Wolfram prints `"Residual" -> 0`; SymPy prints `Φ = 4π C1` and the flux velocity, and no residual field. | **Confirmed.** The SymPy tag's payload has no `Residual` (literal count 0); the Wolfram tag's does. |

**Sound, both legs, each with its own scripts and stdout.** Both legs reran the generator, and it reproduces the
filed measurements file byte for byte (sha256 `cb320c4b…`). Both legs' own enumerations of the comparator output
match the record's row classes, sample counts, `not_formed` reasons and coverage appendix, with nothing dropped or
merged. The single-mechanism restrictions are reported as printed, with their domains and implied `jn`, before the
combined conditions. The `ρ>0` and `ρ≤0` parts keep their labels, and `ρ₀>0` stays a premise. The forward case is
labelled conditional, and no matching outcome, expected value, sign or bending outcome is stated. Adopted premises
carry their source, date, status and dependants, and the density link has no stated effect. D5's `ξ_w''` obligation
is carried as undischarged. The dependences on S8, S11c and O1–O7 have owners, and none is closed. `R-S21-01`
follows the register's schema, passes its rest-on test, names no means and duplicates no entry.

**Not clear.** Findings 1–4 and 6 change what the record and the register claim or hand on. This is the first
repair, so the same author repairs it. The reviewed record is preserved before the repair overwrites it.

## What must be true after the repair
1. Every statement comparing an engine's deflection and radar conditions rests on a cited retrieval. It claims
   identity only for the fields and payloads the retrieval shows to be identical, and states which printed fields
   differ. Whether differing fields state the same condition is listed as open. What the record and the register
   conclude about counting the two observables as separate tests follows from that.
2. The record reproduces the build dispositions' K11 claim limit wherever it reports the shared nonreciprocal
   objects. It states which knives reach those objects, as the accepted harness outputs print.
3. The register's second route, and the record, name every prior-art object the repository identifies for this
   pass's mechanisms, with its source, and state whether each was applied. No oracle is applied in the record.
4. Every command that supports a claim about a file's contents at a revision reads that revision explicitly, so
   that regenerating the measurements file after the record is committed reproduces the filed output.
5. Each printed field the record attributes to an object is attributed to the engines that print it.
6. Nothing else changes. The measurements file is regenerated by the final generator.
