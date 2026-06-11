# Step 3 batching consult — agreed decisions (Claude+Codex, 2026-06-11)

Status: **AGREED — all five questions, no conceptual change, no user gate required.**
Consult shape: per the signed-off adjudication contract (`docs/adversarial_audit_execution_plan.md` §5, operating contract item 5).
Record: prompt `codex_logs/step3_batching_consult_prompt.md` (Claude positions) + full Codex session log `codex_logs/step3_batching_consult_codex.txt` (verdict verbatim at end; Codex independently verified manifest counts, the 408 same-stage/same-parameter duplicate groups, graph coverage stopping at stage 025, and the band prompts' full stage lists before agreeing).

Note: `BATCHES.md` is auto-generated (clobbered on every manifest write), so the decisions live here; extending the `BATCHES.md` renderer to carry batch assignments is part of the D1 tooling below.

## D1 — Dedup/alias pass gates Phase B entry (Q1: AGREE)

Before any Phase B genealogy is built, a flock-gated dedup/aliasing pass maps each of the 922 scanned candidates to a canonical insertion point (same anchor stage + same parameter family + same/nearby cited line ⇒ alias). NO deletions: aliased entries stay in the manifest with `duplicate_of: <canonical_id>`, modality attribution merges onto the canonical, aliases inherit the canonical verdict at Phase C close. Ambiguous cases (same stage+parameter, different value or insertion site) stay separate. Mechanism: Codex builds the tooling and the mechanical alias proposal (workspace-write session); a clean review agent audits the alias map BEFORE it is applied; applied under the manifest lock. Codex measured 408 same-stage/same-parameter duplicate groups, consistent with the ~400–550 post-dedup canonical estimate.

## D2 — Phase B unit = parameter-value family (Q2: AGREE)

Genealogies are built per parameter-value family (one genealogy for e.g. the 54/5 target covering every stage it touches), matching the deployment doc's per-parameter definition and the three-kind taxonomy. Family map derived mechanically (normalized parameter_names + shared literal values across stages), reviewed by a clean agent, recorded as `provenance/_family_map.yaml`. Every canonical candidate maps to ≥1 family or an explicit singleton. `benchmarks.yaml` entries are keyed per family.

## D3 — Phase C batch composition (Q3: AGREE)

Dependency-chain clusters first, band-topical second:
1. Dedicated chain batches: 54/5 + gamma_GR quadrupole-target chain; 37/20 → F1-branch cascade (incl. 4107 closed form); chi_Q/m0_hat/N_Q normalization; Sigma0_can/T_can transport chain (incl. the 867/876 digit transposition); barrier/kill-test chain 222–224 (incl. the external V_known import); calibration chain 245–253 (incl. the lambda_L back-solve + CODATA mass-ratio import); wall-action coefficients + m̂/P0 "fitted scales" gate.
2. Remaining canonical candidates: band-order topical batches of ~25–35.

Deviation from the plan's letter, agreed non-conceptual: chain membership comes from the Phase B family map, NOT `query_graph.py neighborhood` — the atlas graph indexes no stage above 025 (Codex-verified), so graph clustering cannot deliver; the family map preserves the dependency-cluster intent.

## D4 — Parallelism caps (Q4: AGREE)

- Phase C adversarial agents: up to **8 concurrent** clean agents per wave (read-only + report; zero Mathematica seats; orchestrator context is the binding constraint — agents return summaries only).
- Codex defense sessions: resumed per-parameter per the skill contract; at most **2 concurrent sessions that RUN `math -script`**; `.py`-only sessions unlimited; never overlapping an orchestrator Mathematica exec (user-clarified seat rule, 2026-06-05).

## D5 — Phase B entry gate (Q5: AGREE, operationalized)

No additional conceptual gate. Phase B entry = D1 alias map applied + clean-agent review of the map + manifest invariant check. Logged as non-blocking, surfaced to the user at step close:
- Atlas graph blindness above stage 025 (graph modality yielded 28 candidates / 471 gaps): fixing the atlas is out of campaign scope; Phase B builds genealogies from `notes/`, not the graph.
- The completeness critic's "8 silently dropped stages" claim was DISCONFIRMED by the orchestrator (stages present in all rendered band prompts; grep evidence in the consult log) — zero candidates there is honest vacuity, double-confirmed by the critic's own per-stage reads (007, 013, 136, 149, 153, 188, 205, 216).
