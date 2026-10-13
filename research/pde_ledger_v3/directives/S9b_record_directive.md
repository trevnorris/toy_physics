# S9b record: the light-bending step record (directive / pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-13. **Author of the record:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`.

**Gate:** this is a pre-builder decision list. It gets one Codex + Grok pass and is folded once (G2). The record
author starts only after that, and only once the build and the comparator below are accepted. The record is
physics-bearing prose. Two non-author legs (fresh Claude + Grok) review it, source-first, until clear on what may be
claimed. No build directive is in the record author's or the record legs' packet (`CLAUDE.md`, step-record row). The
build objects the record reports are defined below instead.

It was folded once on 2026-10-13: three findings, all accepted
(`directives/_measurements/S9b_record_directive_review_disposition.md`).

**Commits marked ⟨at acceptance⟩** are the acceptance commits of work still under review. The orchestrator fills them
in mechanically when that work is accepted, before the record author starts. Filling them changes no decision here.

## Deliverable

- `steps/S9b_light_bending.md`: the S9b step record, Parts A–C.
- `steps/_measurements/S9b_record_measurements.md`, produced by a committed generator script.
- An S9b pass in `SUBSTRATE_REQUIREMENTS.md`.

## Sources (read-only)

- **The physics spec:** `directives/S9b_SHARED_PHYSICS.md` at `34a994d5`, Parts A–C. Part D is held and not accepted
  (`97580010`). It is not a source for any computed claim.
- **The governing decisions:** `directives/S9b_repair_decision_list.md` at `79055918`: D2 (adopted premises) and D5
  (the O2 record's returned obligation).
- **The accepted constructions** ⟨at acceptance⟩: the two engine sources, unchanged since `503736ff`, their
  production transcripts and the run record `directives/_measurements/S9b_production_runs.md`.
  - The transcripts are git-annex content: run `datalad get` first.
  - The two ablation harnesses and their outputs are evidence about their own engine only.
- **The build's acceptance:** its review disposition ⟨at acceptance⟩, and the earlier build-review dispositions it
  cites. Their notes name limits and items routed to the comparator or the record.
- **The accepted comparator** ⟨at acceptance⟩: `scripts/S9b_cross_engine_comparator.py`, its production output, its
  run record, and its acceptance disposition with the earlier comparator review dispositions.
- **O2:** `steps/O2_steady_brane_balance.md` at `72866fcf`, §8, for the obligation D5 routes here.
- **The register:** `SUBSTRATE_REQUIREMENTS.md`: its entry schema, its requirement/postulate/defect distinction, and
  its population method, including both routes.

## Build objects the record reports

The engines print these under the tags named. Their definitions and premises are stated here, as the build defined
them.
- **Forward case.** Its premise (user, 2026-10-06/07/08): brane material leaves the brane only in throats, so
  `j_n ≡ 0` in the far zone, and the same brane mass crosses every sphere around the body. `Φ ≡ ∮ ρ_br V·n_out dA`
  is the brane mass flux through a far-zone coordinate sphere with outward unit normal `n_out`. It has no supplied
  value or sign, and its relation to `GM` is open. Under the premise, `V` solves the supplied mass balance with
  `ρ_br` live, its integration constant expressed through `Φ` (`F_MASS_SOLUTION`). Part A's observables, Part B's
  objects and conditions, and the branch conditions follow with that `V` substituted, first with `δ ≡ 0`,
  `ξ_w ≡ 0` (`F_FLOW_`), then with `δ` and `ξ_w` live (`F_LIVE_`, which also carries Part C's rewritten
  conditions). `j_n ≡ 0` is the case's premise, so no implied `j_n` is printed for it.
- **Single-mechanism restrictions.** Each is the computed general object with the stated profiles set to zero, then
  reduced. It adds no premise.
  - `R_FLOW_ONLY_`: `δ ≡ 0`, `ξ_w ≡ 0`, `V` live.
  - `R_SPEED_ONLY_`: `V ≡ 0`, `ξ_w ≡ 0`, `δ` live.
  - `R_TILT_ONLY_`: `δ ≡ 0`, `V ≡ 0`, `ξ_w` live.

  These three restrict the deflection, the radar slope, both effective `γ`s and their difference, the residuals,
  and Part B's conditions with their implied `j_n`. `R_BULK_ONLY_` (bulk density alone: `V ≡ 0`, `ξ_w ≡ 0`)
  restricts Part C's fixed-ratio and power conditions and their implied `j_n`, with `ρ_br` symbolic.
- **Density sign.** `ρ₀ > 0`, the asymptotic bulk number density, is a supplied premise (user, 2026-10-11). The local
  bulk density's sign is not a domain predicate. The rows that read `ρ` print the part where `ρ > 0` and the part
  where `ρ ≤ 0` separately, each labelled.

## User decisions this directive records

- **Part D is not built in S9b** (user, 2026-10-08: build A–C, hold D; 2026-10-09: what Part D was to meet is carried
  by the register instead). Which brane material meets the needs light, the drain and gravity place on it is decided
  at the knit, after every force sector, not in S9b.
- **S9b's adopted brane premises are conditional** (user, 2026-10-09). Each is an adopted substrate input, re-tested at
  the knit. None is derived.

## What must be true

1. **The object and its standing.** The record states:
   - what the two accepted constructions emit for Parts A–C;
   - its exact domain and model point, as the spec and the section above state them;
   - which inputs are supplied, which are adopted premises, and which stay live or OPEN.

   An adopted premise or a supplied identification is never presented as derived. S9b states conditions on the local
   wave equation brane light sees near a mass; it does not derive that equation from the substrate.

2. **Cross-engine evidence, as printed.** Every statement about the two constructions' agreement or difference rests
   on the filed comparator output or the transcripts. Each statement cites its retrieval command, and the commands
   with their literal output go in the measurements file.
   - Only retrieval is allowed: existence, verbatim retrieval, literal-match counts, and the shape of a named stored
     object.
   - A question that needs computation beyond retrieval is listed as open, not computed.
   - The record covers every row class the comparator prints, and every printed difference and unresolved coverage.
   - A difference is classified as representational only where a cited source settles it. Otherwise it is listed as
     an open cross-engine difference, with its owner or "no owner named". No difference is reconciled or explained
     away (`CLAUDE.md` M1).
   - Results on the `ρ ≤ 0` parts are reported with their labels, as printed. `ρ₀ > 0` is a premise, not a finding.

3. **The instruments' limits travel with their output.** The record reproduces:
   - the comparator's printed statement of what it compares and does not compare, with every remaining `not_formed`
     reason and the objects it leaves uncompared;
   - the limits named in the comparator's and the build's acceptance dispositions.

   A claim reaches only as far as the compared objects reach.

4. **Each mechanism alone, first.** For the deflection and for the radar slope, the record reports the flow-only,
   speed-only and tilt-only restrictions, and Part C's bulk-density-alone restrictions, before any combined
   condition. Each is reported as printed, with its domain and the `j_n` it implies, and with its
   cross-engine status. Then the combined conditions.

5. **The two observables.** The record states, from printed objects, which combination of inputs each observable's
   Part B condition constrains, and the printed difference between the two effective `γ`s. It counts the deflection
   and the radar slope as separate tests only where the printed objects show their conditions differ.

6. **The forward case.** The record reports the forward case as conditional on its premise (above). Every forward
   object carries that label. Its conditions are computed relations between `Φ`, the live profiles and `GM`, not supplied
   ones, and no matching outcome is claimed.

7. **Conditional premises.** The record lists each adopted premise S9b uses or adopted: the forward case's
   no-far-zone-loss premise, w-parity, P5, P6 and the density link (decision list D2), and O2's P1–P4 as they enter.
   For each: its source and date, its status (adopted, conditional, re-tested at the knit), and which printed objects
   depend on it. A premise scoped to Part D enters no computed object, and the record says so. For the density link
   (`c_γ² = μ_⊥/ρ_br`, with `μ_⊥ ∝ ρ_br^α` and `α` live), the record also states that Parts A–C keep the speed
   change and the brane density as independent live profiles, and states no effect of the link.

8. **What Part D was to meet.** If the register's rest-on test admits it, the S9b pass adds an entry whose
   requirement is the printed Part B condition for each observable, with its Part C rewrites: the combination of
   live profiles it constrains relative to `GM`, and the `j_n` it implies, with `ρ_br` symbolic and every quantity
   the print leaves free kept free. The entry:
   - cites the printed conditions and restrictions it rests on;
   - names no means of meeting the requirement: no balance, material, premise or candidate;
   - states no effect of the density link;
   - has a step id or "no owner named" as its target, and a status from the register's schema.

   If the rest-on test does not admit it, the record says why.

9. **The O2 obligation.** The record carries D5's obligation as undischarged, owned by the orchestrator: the
   WL-only `ξ_w''` dependence with the seven first-derivative keys at the same locations. S9b, without Part D, does not
   discharge it. A register entry follows only if a sourced requirement is established.

10. **Remaining dependences.** Each dependence the spec lists as outside this step, and S8, S11c and O1–O7 where
    they enter, is stated with its relation to this result and its owner, as the sources state them. None is closed
    here.

11. **The S9b register pass.** Populate the S9b pass by the register's own schema, rest-on test and method, as
    written.
    - Keep the distinction between a derived requirement and an adopted substrate premise.
    - Merge into existing entries; don't duplicate them.
    - Where this directive and the register disagree, the register wins, and the disagreement is reported.

## Exclusions

- No new CAS, derivation, oracle check, closure, normalization or profile choice.
- No engine, harness or comparator edit, and no rerun of any.
- No Part D calculation.
- Edit only the record, its measurements file and generator, and `SUBSTRATE_REQUIREMENTS.md`.
- No commit and no push.

## STOP report (at most one page)

- The record's headline claims, each with its source.
- The single-mechanism results, each with its cross-engine status.
- Every printed cross-engine difference and unresolved coverage, and its classification.
- Register entry counts before and after, new entries by target, and existing entries that gained a source.
- Disagreements between sources, with both sides quoted.
- Open questions that need computation beyond retrieval.
