# S9b record directive: G2 review dispositions (orchestrator)

The directive is the orchestrator-written pre-builder decision list for the S9b step record, so the legs are Codex and
Grok, one pass, folded once (`CLAUDE.md` G2). Lookups: `S9b_record_directive_review_disposition_lookups.md`
(mechanical; generator `_scratch/s9b_build/gen/s9b_record_dl_review_lookups.sh`). The directive's own claims about its
sources are grounded in `S9b_record_directive.md` in this directory (generator
`_scratch/s9b_build/gen/s9b_record_directive_lookups.sh`).

**Artifact.** The working-tree directive at sha256 `9ad9611a…`, preserved as
`_scratch/s9b_build/preserved/S9b_record_directive_reviewed_v0.md`. Prompt:
`_scratch/s9b_build/s9b_record_directive_review_prompt.md`. Baseline:
`_scratch/s9b_build/s9b_record_directive_review_baseline.sha256`.

**Legs.**
- **Codex (gpt-6.1-sol, xhigh): NEEDS REVISION, one finding.** Final message:
  `_scratch/s9b_build/s9b_record_directive_review_codex_final.txt`; log `s9b_record_directive_review_codex.txt`.
- **Grok: NEEDS REVISION, two findings.** Report: `_scratch/s9b_build/s9b_record_directive_review_grok.txt`. Grok found
  Parts A–C, the forward case, the single-mechanism restrictions, the cross-engine handling, the O2 obligation and the
  retrieval-only rule consistent with the sources.

| # | Finding (leg) | Disposition | Folded text |
|---|---|---|---|
| C1 | The sources list the build directive. The step-record row of `CLAUDE.md` puts no build directive in a record's packet. (Codex 1) | **ACCEPT.** Lookups: the reviewed lines 26–27 and the `CLAUDE.md` row. | The build directive is removed from the sources. A new section, "Build objects the record reports", states the forward case, its premise and `Φ`, the four restriction families and the density-sign premise, from build directive items 9 and 10 and its tag table (grounded in `S9b_record_directive.md`). The gate says no build directive is in the record author's or the record legs' packet. |
| G1 | Item 8 names the entry's object as what "the brane's steady in-plane force balance must deliver", which is a route, and "targets the knit". The register's `requirement` is the object, not a derivation path; `target` is a step; status is OPEN, DELIVERED or RETIRED; a forward consumer's need is not an entry; D5 admits an entry only for a sourced requirement. Accepted Part B is a condition on the profiles relative to `GM` with its implied `j_n`, and no momentum balance is supplied or established (spec, O2 §8). (Grok 1) | **ACCEPT.** Lookups: the reviewed item 8, the register's schema and calibration lines, spec Part B, D5, O2 §8. The "force balance must deliver" wording was the orchestrator's paraphrase of the user's 2026-10-09 decision, not a sourced object. | Item 8, "What Part D was to meet": if the rest-on test admits it, the entry's requirement is the printed Part B condition per observable, with its Part C rewrites, and its implied `j_n`, with `ρ_br` symbolic and every printed free quantity free. It names no means, states no effect of the density link, and takes a step id or "no owner named" as target and a schema status. If the test does not admit it, the record says why. |
| G2 | Item 8 carries P5, P6 and the density link onto the entry, and names power-law creep and a two-component supersolid as knit candidates. D2 scopes P5, P6 and the density link to Part D and keeps `μ_⊥(x)` and `ρ_br(x)` independent in Parts A–C. P6 excludes power-law creep. No listed source names the other candidates. (Grok 2) | **ACCEPT.** Lookups: D2 rows, spec lines 374 and 539, and a search of the sources for the candidate names (none). | Item 8 names no premise, material or candidate. Item 7 keeps P5, P6 and the density link as Part D-scoped premises that enter no Parts A–C object, and now states that Parts A–C keep the speed change and the brane density as independent live profiles. |

**Routing.** The density link and the material candidates the user has raised are not lost. They are user ideas,
not sourced requirements, and the light register pass after this record lists parked ideas as one-line candidates.
They go there, with their dates and no effect stated.

**Folded once.** The folded directive is at sha256 `cdead60d…`. No second pass (G2).
