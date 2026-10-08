# O2 sub-step 7: the O2 record (directive / pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-07. **Author of the record:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`.

**Gate:** this is a pre-builder decision list. It gets one Codex + Grok pass and is folded once (G2), and the record
author starts only after that. It was folded once on 2026-10-07: two findings, both accepted
(`directives/_measurements/O2_record_directive_review_disposition.md`). The record it produces is physics-bearing prose. Two non-author legs (fresh Claude +
Grok) review it, source-first, until clear on what may be claimed (O2 scoping §7 row 7).

## Deliverable

`steps/O2_steady_brane_balance.md`: the record O2 scoping §7 row 7 asks for
(`directives/O2_steady_brane_balance_scoping.md:253`). Its O2 pass entries go in `SUBSTRATE_REQUIREMENTS.md`.

## Sources (read-only)

- **Scope and premises:**
  - `directives/O2_steady_brane_balance_scoping.md`, §1 and §7;
  - `directives/O2_premise_decision_list.md` (`77d2c39a`);
  - `directives/O2_input_contract.md` (`217a92e9`).
- **The physics spec:** `directives/O2_SHARED_PHYSICS.md` (`4680e251`).
- **The accepted constructions** (`0f2e0af8`): the two engine sources, and their production transcripts and run
  record (`8cd40f59`; `directives/_measurements/O2_production_runs.md`).
  - The transcripts are git-annex content: run `datalad get` first.
  - The ablation harnesses' knife triples are evidence about their own engine only.
- **The engines' acceptance:** `directives/_measurements/O2_build_r3_review_disposition.md`. Its notes route
  named items to the comparator.
- **The accepted comparator** (`48c9bf33`): `scripts/O2_cross_engine_comparator.py`.
- **The comparator's production output** (`f435dd88`; run record `directives/_measurements/O2_comparator_production_run.md`):
  - `scripts/out/O2_cross_engine_comparator.out`, the comparison stream (git-annex content; run `datalad get`
    first);
  - `scripts/out/O2_cross_engine_comparator_accounting.jsonl`;
  - `scripts/out/O2_cross_engine_comparator_catalog.json`.
- **The comparator's acceptance:** `directives/_measurements/O2_comparator_build_r5_review_disposition.md`. Its
  notes name the limits carried to this record.
- **The register:** `SUBSTRATE_REQUIREMENTS.md`: its entry schema, its requirement/postulate/defect distinction, and
  its population method, including both routes.

## What must be true

1. **The object and its standing.** The record states:
   - what the two accepted constructions emit (spec §9);
   - its exact domain and model point (spec and input contract);
   - which inputs are user-adopted premises (premise decision list 1–4), which are supplied, and which stay OPEN.

   Each item is stated as the sources state it. An adopted premise is never presented as derived. O2 is a
   conditional balance, not a derived brane law.

2. **Cross-engine evidence, as printed.** Every statement about the two constructions' agreement or difference
   rests on the filed comparator output or the transcripts. Each statement cites the retrieval command, and the
   commands with their literal output go in `steps/_measurements/O2_record_measurements.md`, produced by a
   committed generator script.
   - Only retrieval is allowed: existence, verbatim retrieval, literal-match counts, and the shape of a named
     stored object.
   - A question that needs computation beyond retrieval is listed as open, not computed.
   - The record covers every row class the comparator prints, and every printed difference.
   - A difference is classified as representational only where a cited source settles it. Otherwise it is listed
     as an open cross-engine difference, with its owner or "no owner named". No difference is reconciled or
     explained away (`CLAUDE.md` M1).

3. **The instrument's limits travel with its output.** The record reproduces:
   - the comparator's printed statement of what it compares and does not compare;
   - the limits named in its acceptance disposition.

   For each item the engines' acceptance routed to the comparator, the record states what the comparator output
   shows about it, or that the comparator does not compare it.

4. **Remaining dependences.** O1 and O3–O7 are each stated with its relation to this result and its owner, per the
   scoping §1 table and the sources. No dependence is closed here.

5. **Register entries.** Populate an O2 pass in `SUBSTRATE_REQUIREMENTS.md` by the register's own schema, rest-on
   test and method, as written.
   - Keep the distinction between a derived requirement and an adopted substrate premise.
   - Merge into existing entries; don't duplicate them.
   - Where this directive and the register disagree, the register wins, and the disagreement is reported.

6. **What S9b Part D may use.** The record identifies which O2 objects and conditional inputs a later S9b Part D
   calculation may use, and under which conditions.

## Exclusions

- No new CAS, derivation, oracle check, closure, normalization or profile choice.
- No engine or comparator edit, and no rerun of either.
- No S9b calculation.
- Edit only the record, its measurements file and generator, and `SUBSTRATE_REQUIREMENTS.md`.
- No commit and no push.

## STOP report (at most one page)

- The record's headline claims, each with its source.
- Every printed cross-engine difference and its classification.
- Register entry counts before and after, new entries by target, and existing entries that gained a source.
- Disagreements between sources, with both sides quoted.
- Open questions that need computation beyond retrieval.
