# Light-leakage scoping inventory: gate record

**Accepted artifact:** `steps/LIGHT_LEAKAGE_SCOPING.md`, a byte copy of v5 (sha256 `17293320…a028`,
`scoping_v5.sha256`). It is an inventory for a later worst-case/best-case leakage bracket, not a calculation or a
ledger result. **Role:** physics-bearing scoping input, reviewed until clear (G4). **Accepted:** 2026-10-07, under
the user's rule that a version is accepted when nothing outstanding changes what the bracket would compute or
claim (`scoping_r5_dispositions.md`).

**Directive:** `light_leakage_scoping_directive.md`, written by the orchestrator. It had one Codex + Grok pass
(`light_leakage_scoping_directive_review_prompt.md`; `directive_review_codex_final.txt`,
`directive_review_grok.txt`; both NEEDS CHANGE) and was folded once.

## Versions, authors and legs

Legs are chosen by authorship (G1): Codex-written versions go to fresh Claude + Grok; Claude-written versions go
to Codex + Grok. Each round used one identical prompt for both legs.

| Version (sha256 prefix) | Author | Review prompt | Legs and verdicts | Dispositions |
|---|---|---|---|---|
| v0 `0386753e` | Codex (gpt-6.1-sol, xhigh, web) | `light_leakage_scoping_review_prompt.md` | Claude NEEDS REVISION (`scoping_review_claude.md`); Grok NEEDS REVISION (`scoping_review_grok.txt`) | `scoping_r1_dispositions.md` |
| v1 `d6101683` | Codex, repair (`scoping_repair1_codex_final.txt`) | same | Claude NEEDS REVISION (`scoping_review_r1_claude.md`); Grok NEEDS REVISION (`scoping_review_r1_grok.txt`) | `scoping_r2_dispositions.md`; the Codex repair bred defects in changed material, so the **author changed** |
| v2 `7adcfdc1` | fresh Claude (`scoping_repair2_claude_author_report.md`) | `…_prompt_r2.md` | Codex NEEDS REVISION (`scoping_review_r2_codex_final.txt`); Grok NEEDS REVISION (`scoping_review_r2_grok.txt`) | `scoping_r3_dispositions.md`; stopped after more than 2 rounds (`scoping_r3_status_for_user.md`); the user chose one more fold |
| v3 `85f1d09c` | fresh Claude (`scoping_repair3_claude_author_report.md`) | `…_prompt_r3.md` | Codex NEEDS REVISION (`scoping_review_r3_codex_final.txt`); Grok CLEAR (`scoping_review_r3_grok.txt`) | `scoping_r4_dispositions.md`; two findings verified; the user chose a targeted fix |
| v4 `214d1fe6` | fresh Claude (`scoping_repair4_claude_author_report.md`) | `…_prompt_r4.md` (scoped to the diff) | Codex NEEDS REVISION (`scoping_review_r4_codex_final.txt`); Grok CLEAR (`scoping_review_r4_grok.txt`) | `scoping_r5_dispositions.md`; one finding verified, in a clause the v4 author added; the user chose a targeted fix by a different author |
| **v5 `17293320`** | a different fresh Claude (`scoping_repair5_claude_author_report.md`) | `…_prompt_r5.md` (scoped to the diff) | **Codex CLEAR** (`scoping_review_r5_codex_final.txt`); **Grok CLEAR** (`scoping_review_r5_grok.txt`) | `scoping_r5_dispositions.md`, section "v5 scoped check and acceptance": **ACCEPTED** |

The frozen reviewed baselines v0–v4 are here as `light_leakage_scoping_v0.md` … `_v4.md`. The hash files
`scoping_v0.sha256` … `scoping_v5.sha256` give the reviewed hash of each version.

## What this record does not include

- **Paths.** The dispositions and prompts cite paths under `_scratch/light_leakage/`, where the review ran. The
  files with those names are here.
- **Codex transcripts.** For each Codex run, `*_codex_final.txt` holds only the final report: the transcript text
  after its last `tokens used` count line (`tail -n +$((n+2))`). The full transcripts are raw logs and stay outside
  the tree, in `_scratch/light_leakage/` (gitignored).
- **Carried item.** The v3 author noted that the c1 record has no "Changed after close" note, though `b:39` says
  the inertia repair propagated. That is a records-housekeeping item for the S11c review bucket, not part of this
  inventory.
