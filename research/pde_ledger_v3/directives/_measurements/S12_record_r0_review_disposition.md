# S12 record, review round 0: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S12_record_r0_review_lookups.md`, beside this file (generator
`_scratch/s12/gen/rec_r0_lookups.sh`). "K-n" cites a section of that file.

## The round

- **Artifact.** Codex wrote it, on top of `560c6c51`, under `directives/S12_record_decision_list.md`. It is preserved
  unchanged in the commit that adds this file. It has six files:
  - `steps/S12_order_source.md`;
  - `steps/_measurements/S12_record_measurements.md`;
  - `scripts/S12_record_measurements.py`;
  - the S12 pass in `SUBSTRATE_REQUIREMENTS.md`;
  - the status-line edits in `../../STATUS.md` and `steps/S9b_flow_oracle_comparison.md`.

  All six matched the baseline after the round, and the working tree was unchanged (K1).
- **Prompt and legs.** `_scratch/s12/rec/s12_record_review_prompt_r0.md`, sha256 `a0f00d94…` (K1). Both legs
  received it identically, in parallel:
  - a fresh Claude agent (Opus), report `/var/projects/_review_scratch/s12-record-r0-claude/REVIEW.md`;
  - Grok (`grok-4.7`, effort high), report `/var/projects/_review_scratch/s12-record-r0-grok/REVIEW.md`.

  The record is Codex-written, so these are its two non-author legs.
- **Verdicts.** Claude: NEEDS REVISION, 2 findings. Grok: NEEDS REVISION, 1 finding (K2).
- **Both legs** reran the generator to their scratch directories and reproduced the filed measurements file byte for
  byte. Each enumerated the comparator output itself and found every row class, count and difference family in the
  record.

## Findings

**C1. What S12 does not deliver is still shown as S12's to deliver.** (Claude 1.) **CONFIRMED.**
- **The sources.**
  - The register's schema defines `target` as "the step that must deliver it". It calls a requirement "an obligation
    on a *future* step", and its method keys each entry "by the step that must deliver" (K3).
  - The record says "S12 record completion is not delivery of those functions" (K3).
- **What the pass hands on** (K3).
  - `R-S12-02` keeps `target` S12. It keeps the owner line "S12 owns the drain/return data … and is next for dynamical
    bulk-to-brane order conversion with separate boundary data".
  - The by-target table lists `R-S12-02` under **S12** with no qualifier. `R-S12-01` in the same table carries
    *(register inference; no live owner named)*.
  - The record gives "the S12 inventory" as owner of the source/flux forms, the controller-role assessment, the
    boundary laws and data, the native O6 identifications and the conversion partners (§5 l.261, §6 l.271 and l.276).
  - `STATUS.md:65` lists **S12** as owner of the dynamical drain/return source and its boundary data.
- **The gap.** Once this record is accepted, these items would show as owned by a step that has finished without
  delivering them.
- **Part of a handoff is lost.**
  - The S9b oracle record splits the material-destination question: what the drain produces "is again S12's object
    (R-S12-02), and whether the knit can meet it belongs to R-S21-01" (K3, oracle record l.166).
  - The record's §7 says only "material destination with S12/O2" (K3, l.315).
- **What must be true.**
  - In the register, the record and `STATUS.md`, no obligation S12 does not deliver shows S12 as the step that will
    deliver it.
  - Where no source names a later step, the obligation is OPEN with no live owner named, in the register's own
    convention. The record lists who should own it as a question for the user.
  - The disagreement between the earlier S12 assignment and the user's framing stays reported.
  - Each handoff the sources split keeps every part, with its owner as the source gives it.

**C2. Differences the comparator's contract settles are handed on as open.** (Claude 2.) **CONFIRMED.**
- **The record and its hand-on** (K4).
  - §4 says the six unrestricted coefficient operands differ "in `regular/above` for face 1 and `regular/below` for
    face 2". They are "explicitly **not delta parts**". It then says "No cited disposition settles the unrestricted
    expressions' difference: **OPEN cross-engine diagnostic differences; no owner named**".
  - §6 lists "adjudication of the six" as needing computation.
  - `SUBSTRATE_REQUIREMENTS.md:593` and `STATUS.md:90` carry the six forward.
- **The source.** The accepted comparator contract, item 3, says (K4):
  - "A delta's coefficient carries data about the distribution only on that delta's support."
  - An unrestricted raw coefficient "is a separately labelled representation diagnostic, never the delta's part."
  - Its equivalence control changes a coefficient "away from its face only, so the distribution is unchanged".
- **The filed output** (K4). In each of the six, the only stored parts whose as-computed and supplied-domain
  residuals are not the literal `"0"` are:
  - `regular/above` (outside above) for face 1, whose delta sits at `w₋`;
  - `regular/below` (outside below) for face 2, whose delta sits at `w₊`.

  Under the supplied strict order `w₋ < w₊`, each lies off its delta's support. Every other stored part is the
  literal `"0"` on both fields, the coefficient's own delta parts included.
- **The leg disagreement.** Grok listed the six as soundly kept open. The contract's wording and the stored parts
  settle it as above.
- **What must be true.**
  - The record classifies the six as differences in a labelled diagnostic, confined to parts off each delta's
    support. It cites the contract's statement that values there carry no data about the distribution.
  - Both raw expressions stay printed and retrievable.
  - No statement in the record, the register or `STATUS.md` hands them on as an open difference in the S12 object.
  - The record claims nothing about why the expressions differ beyond what a retrieval shows.

**G1. Component units are attributed to both transcripts, and one does not print them.** (Grok 1.) **CONFIRMED.**
- **The record.** The record's item 3 says "They assign density `(-4,0)`, order `(0,0)`, velocity `(1,−1)`, faces
  `(1,0)` and rates `(0,−1)`", citing the transcripts' retrievals M03 and M10 (K5).
- **The transcripts** (K5).
  - SymPy prints these in `primitive_function_dimensions`: `n` `[-4],[0]`, `u1` and `uw` `[1],[-1]`, `wminus` and
    `wplus` `[1],[0]`.
  - The comparator leaves that tag unjoined: "tag absent in other engine".
  - The Wolfram transcript has no such tag, and no `{-4, 0}` or `{1, -1}`.
  - The Wolfram source assigns the same pairs in its `unit` switch (l.29–32), which the record's M04 covers.
- **The hand-on.** `R-S1-01` gives "the density's assigned dimension is `L⁻⁴`" (K5).
- **What must be true.**
  - Each unit assignment is attributed to the construction or transcript that shows it: the sources' assignments,
    and SymPy's printed table, which is unjoined.
  - Only dimension objects that both transcripts print are presented as printed by both, or as compared.
  - `R-S1-01` rests the density's dimension on what shows it.

## Outcome

Three findings, all confirmed. None changes a computed object. They change what the record, the register and
`STATUS.md` hand on and attribute. The reviewed record is preserved as-is in this commit. A repair round by the same
author follows, and then round 1, with the same pairing (fresh Claude + Grok). Round 1 is review round 2 of the
record, so if anything remains after it, the user decides how to finish.
