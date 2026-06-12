# Phase C (Step 5) execution decisions — Claude+Codex consult, 2026-06-11

Status: **AGREED — sound to launch with the refinements below. No conceptual change, no user gate** (batch-granularity class per execution plan §5 / operating contract item 4).
Consult record: prompt `codex_logs/step5_phase_c_granularity_consult_prompt.md` + full Codex session log `codex_logs/step5_phase_c_granularity_consult_codex.txt` (session `019eba04-8d63-7d61-9bc6-528f35b7b4a5`, read-only; Codex independently re-measured 1017 bundles / 60 published_target / 27 high / 65-candidate union before agreeing).

This refines, does not replace, `BATCHING_DECISIONS.md` (D1–D5). It settles HOW Phase C is executed.

## PC1 — Adversarial unit = disjoint audit-unit partition (parameter-value family), with per-member sub-verdicts

One clean adversarial agent per audit unit, handed the FROZEN directive + EVERY member candidate's provenance bundle + underlying sources + the unit's `benchmarks.yaml` entries + graph context. Deviation from the plan's literal "one per candidate" — non-conceptual (changes how the same claims are grouped/checked, not what is claimed; every candidate still gets a verdict).

**Codex refinement applied:** the units must be a *disjoint partition with per-member/per-parameter sub-verdicts*, NOT blind `candidate_family_map` inheritance, because the family map carries overlapping seeded overlays + primary families. Each agent MUST return a separate YES/NO line per member candidate (and per parameter for compound bundles), not a single family verdict. Specifically `chain_calibration_245_253` mixes distinct external-fit questions (λ_L back-solve at 247 vs CODATA mass proxy at 250) — kept in one agent ONLY because the output is per-member.

Disjointness verified: the 65 batch-1 candidates map to exactly 28 units, member counts sum to 65, no candidate double-counted. Roster: `_phase_c_batch1_roster.yaml`.

## PC2 — Batch 1 = the 28-unit `published_target ∪ HIGH` cluster (external-fit + central adjudication + named HIGH defects)

Trigger `published_target ∪ severity:high` confirmed correct for the first batch (Codex: no external-facing omission requires pulling `c0=3/4`, `n=5`, or `γ₀` forward — those carry no external-match claim for the fit-vs-derive test to bite on; they wait for later band-topical batches).

Explicit precedence partition (seeded chains first, then primary/small units): `chain_calibration_245_253` (6) · `chain_chi_Q_norm` (6) · `chain_quad_54_5` (20) · `chain_aspect_37_20` (4) · `chain_barrier_222_224` (1) · `chain_wall_action` (1) → then the 22 primary/small/compound units (incl. `fam_0290_xi_turn`, the compound `fit_stage247_session1_packet`, `fit_stage250_goldilocks_inputs`). NOTE Codex caught: `Xi_turn` is NOT in the calibration seeded chain — it is its own `fam_0290_xi_turn` (already separate in the roster).

## PC3 — Verdict routing

- YES on any member → Codex defense (per-parameter sessions, scoped to the implicated member/parameter, not the whole unit) → adjudication consult → FIND_STANDS / FIND_FAILS / PARTIAL. **Every FIND_STANDS → user; hard-stop that member's dependency cone until adjudicated.**
- All-NO unit → all members `audited`.
- `chain_aspect_37_20` (L/a = 37/20) goes to the adjudication consult + user gate **regardless** of the adversarial binary — it is a disclosure/classification call (free_choice vs published_target), FIND_STANDS-class by construction. Bundles classify it free_choice while flagging the unresolved published-target concern (073:32, 121:51).

## PC4 — Rendering mechanics (skill contract preserved)

`phase-c-render <representative_candidate_id>` produces the skill-standard prompt (frozen-directive embed + provenance + benchmark + graph) and moves the representative `provenance_built → audit_pending`. The agent is additionally handed ALL member bundle paths and must cover every one with a per-member sub-verdict. Status discipline (Codex Q4): advance member statuses only after the WHOLE unit's triggered bundles are covered — at dispatch move all members to `audit_pending`; at verdict move all together to `audited` (NO) or the implicated members onward to `defense_pending` (YES). State machine is single-step-forward, so this is scripted as serial flock-gated `set-status` loops.

## PC5 — Pilot first

Run `chain_quad_54_5` (the 20 batch-1 P0_target=54/5 members — the dominant external GR back-solve, duplicated across scanner angles, with sourced benchmarks) as a single pilot unit FIRST; verify verdict quality + that the rendered prompt carried everything; THEN fan out the remaining 27 units ≤8 concurrent.

## PC6 — Remaining Phase C scope after batch 1 (USER-DECIDED 2026-06-11): TARGETED sweep

Batch 1 covered the entire external-fit surface (all 60 `published_target`) + all 27 HIGH findings = everywhere the fit-vs-derive test bites, and found NO fatal flaw. The user chose the **targeted sweep** (not exhaustive, not stop):
- **DO:** Phase C adversarial agents over the remaining **~190 `free_choice` ansätze** (where a posit-dressed-as-derived overclaim could hide — exactly the stage_192 failure mode), batched band-topically ~25–35 → ~6–8 batches. PLUS a **completeness-critic spot-check** across the `internal_consistency` band.
- **SKIP:** exhaustive per-family adversarial agents on the ~685 honest `internal_consistency` derivations — already verified by the layer-1 red-team ×2; the fit-vs-derive test has near-zero surface there. The critic spot-check is the safety net; log what is intentionally not deep-audited (no silent caps).
- Same execution mechanics (PC1–PC5): disjoint family partition, per-member sub-verdicts, sourced benchmarks, ≤8 concurrent, /compact between batches.

**User adjudication sign-offs (2026-06-11): ACCEPT BOTH.** (1) L/a=37/20 = `free_choice` confirmed (apply the 3 r_F1 bundle label corrections `published_target→free_choice` in Step 6). (2) stage_192 card overclaim = PARTIAL, non-fatal card-text fix in Step 6. The full Step-6 close-out fix queue is in `verdicts/batch_c1_adjudication.md`.
