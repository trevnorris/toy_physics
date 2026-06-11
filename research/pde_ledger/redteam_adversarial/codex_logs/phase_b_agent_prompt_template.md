# Phase B genealogy agent — prompt template (4C)

This is the reusable clean-agent prompt for the 28-band genealogy fan-out. The orchestrator fills `{BATCH}`, `{STAGE_LO}`, `{STAGE_HI}`, `{CANDIDATE_LIST}`, and any batch-specific chain notes, then passes it inline via the Agent tool. Stored here for resume/audit. Each agent is a fresh clean context (no campaign contamination).

---

You are a provenance-genealogy auditor for a 253-stage toy-physics derivation ledger (`research/pde_ledger/`). This is layer 2 of a falsification stack: your job is to classify, for each fitted/posited parameter value, **how it was constrained** — the fit-vs-derive discriminator. Internal-consistency red-teams cannot see claim-level overstatement; that is exactly what you are checking.

ALL paths below are relative to the project root `/var/projects/toy_physics/research/pde_ledger/`. Work only by READING files and REASONING. Do NOT run python/scripts to "compute" anything; do NOT edit any source file. You only WRITE the synthesis YAMLs described at the end.

## Your batch: {BATCH} (stages {STAGE_LO}–{STAGE_HI})

For each candidate id below you will produce one provenance synthesis YAML.

{CANDIDATE_LIST}

## Per-candidate procedure

For each candidate id `CID`:

1. **Read its pending bundle** at `redteam_adversarial/provenance/CID__<param>.yaml`. (A candidate may have >1 bundle if `multi_target` — handle each.) The bundle gives you:
   - `parameter_name` — the resolved primary target. **Use this string verbatim** as `parameter_name` in your output; do not invent your own.
   - `source_evidence` — pointers (path/line/excerpt) into the stage's paper card, per-stage notes, and audit scripts.
   - `anchor_stages`, `target_values` (the literal value(s) the scanner extracted), and `graph_context`.
2. **Open the per-stage `notes/` source.** This is a BINDING RULE: you may not assert an origin ("introduced at stage N, line L") from cards/scripts/outputs/graph alone — you must open and cite a file under `notes/` (per-stage notes live at `notes/stages/moving_throat_pde_stage<NNN>_*.md`). Also read the paper card `paper/stages/stage_<NNN>.tex` and, when the value's status is in question, the audit script `scripts/moving_throat_pde_stage<NNN>_*_sympy_audit.py`.
   - **Some early stages (notably 004–020) have NO per-stage `notes/stages/` file.** If the bundle's `source_evidence` lists no `notes_stage` entry for your stage, search the broader topical notes for where the value is documented before concluding a gap: `notes/moving_throat_notes_full.md`, `notes/moving_throat_pde_*.md`, `notes/pde/`, `notes/em_projected/`, `notes/barrier_audit/`, `notes/summaries/` (grep these for the parameter symbol / stage number). Cite whichever `notes/` file actually documents the value. Only if NO file under `notes/` documents it do you log a `provenance_gap`. A citation under `notes/` is required regardless of which notes file it comes from — the paper card and scripts alone never satisfy the binding rule.
3. **Determine the genealogy:**
   - **origin** — the stage+line where THIS parameter's value first ENTERED the ledger (may be an upstream stage that the current stage merely carries/restates). Cite the `notes/` line that establishes it.
   - **constraint_kind** — the central call (see taxonomy below). Cite the `notes/` line that justifies the classification.
   - **downstream_dependents** — stages that consume this value (from the card's "Downstream Use" field / graph_context / notes).
   - **findings** — any contradiction is itself a finding (see below); never silently resolve one.

## constraint_kind taxonomy (choose exactly one per constraint; this is the deliverable)

- `internal_consistency` — the value is FORCED by requiring consistency with already-derived ledger quantities (solved from a pole/normalization/determinant/matching condition among quantities the ledger itself produced). A genuine derivation. Evidence: the notes show the equation it solves and that all inputs are themselves ledger-internal.
- `published_target` — the value is chosen or back-solved to MATCH an EXTERNAL published number (a GR coefficient, a CODATA constant, a known 5PN/PN result, a textbook value). This is a FIT to an external target — the load-bearing case for this audit. Evidence: the notes name the external target and/or show a back-solve to hit it. (E.g. λ_L back-solved at 247, CODATA 1836.152… imported at 250, V_known barrier imported at 222–224, the 54/5 / 37/20 targets.)
- `free_choice` — the value is POSITED / assumed / a convention / a gauge or branch selection, with neither a derivation nor an external anchor. An ansatz. Evidence: the notes posit it without justification beyond convenience/convention.

If the evidence is genuinely ambiguous between two kinds, pick the best-supported and record the tension in `provenance_findings` (type `constraint_kind_ambiguous`) with the competing citations. Do NOT default everything to `internal_consistency`.

## When notes/ does not support an attribution

If no `notes/` text supports an origin or a constraint classification for a parameter, do NOT fabricate one. Leave that claim out and emit a `provenance_findings` entry of type `provenance_gap` (severity `needs_triage`) naming what is missing and where you looked. A candidate with only a logged gap is a valid, ingestable result — but use gaps only when the notes genuinely lack support, not to avoid work.

## Findings to raise (each is a `provenance_findings` entry)

- `self_contradictory_origin` — the sources disagree on where/what the value is (the χ_Q=1 / "097 vs 105" class of bug).
- `derive_vs_fit_mismatch` — the card/notes CLAIM a value is derived/forced, but the actual mechanism is a back-solve to an external number (claimed `internal_consistency` but really `published_target`). This is the highest-value find.
- `stale_provenance_anchor` — an "introduced at Stage N" / cross-reference that points at the wrong stage (numbering drift).
- `paper_card_overclaim` — the paper card states more certainty/derivation than the notes support.
- `provenance_gap` / `graph_gap` — missing notes source / missing atlas node.

Every finding needs `type`, `severity` (`needs_triage` | `low` | `high`), `summary`, and `citations` (a list of `{path, line, excerpt}`; for findings the path need not be under notes/ but should be a real file).

## Output — write one YAML file per (candidate, parameter)

Write to `redteam_adversarial/provenance/_synthesis/{BATCH}/CID__<param-slug>.yaml` (create the dir). Use the SAME `<param-slug>` as the bundle filename. Exact schema (the ingest tool validates it strictly):

```yaml
candidate_id: CID                      # verbatim
parameter_name: <param>                # verbatim from the bundle
origin_claims:                         # list; may be empty if logged as a gap
  - parameter: <param>
    introduced_at_stage: <int or NNN>
    introduced_at_line: <int>
    citation:                          # REQUIRED, path MUST be under notes/
      path: notes/stages/moving_throat_pde_stage<NNN>_<...>.md
      line: <int>
      excerpt: "<verbatim line text>"
constraints:                           # list; one per parameter you classify
  - parameter: <param>
    constraint_kind: internal_consistency | published_target | free_choice
    evidence_citation:                 # REQUIRED, path MUST be under notes/
      path: notes/stages/moving_throat_pde_stage<NNN>_<...>.md
      line: <int>
      excerpt: "<verbatim line text>"
downstream_dependents: [ <stage ints> ]
provenance_findings:                   # list; may be empty
  - type: derive_vs_fit_mismatch
    severity: high
    summary: "<one line>"
    citations:
      - { path: <real file>, line: <int>, excerpt: "<verbatim>" }
```

Write each file as you finish it (incrementally), so partial progress survives. Emit YAML only — no JSON.

## When done

Reply with a compact YAML report: per-candidate one line `CID: <constraint_kind or GAP> [+N findings]`; a list of any candidate you could not complete and why; and a tally of constraint_kind counts + findings-by-severity for the batch. Do not paste the full synthesis files back — they are on disk.
