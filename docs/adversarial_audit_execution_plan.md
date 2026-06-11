# Adversarial Audit — Execution Plan

**Status:** ACTIVE. Authorized by the user 2026-06-10. This is the operational driver for the layer-2 adversarial fit-vs-derive campaign; design rationale and phase definitions live in `docs/adversarial_audit_deployment.md` (the deployment doc). Read that doc before executing any step here.

**How to use this file:** one step at a time, top to bottom. Each step has entry criteria, an exit gate, and a status line. Update the status line when a step closes. A fresh session resuming after `/compact` starts by reading this file, the deployment doc, and the memory anchor `[[project-adversarial-audit]]`, then continues at the first step not marked DONE.

---

## Operating contract (user-authorized 2026-06-10)

1. **Codex codes, Claude reviews.** Codex applies ALL file changes that carry content: paper TeX, notes, scripts, the skill's runtime code, provenance/benchmark YAML built by agents under directives. Claude orchestrates, writes directives/verifications/manifests/trackers (red-team-scaffolding class), reviews and verifies everything Codex produces. Claude never hand-applies a directive fix. This split is the calibrated process; deviating from it is never the orchestrator's call.
2. **Auto mode.** Claude completes each step autonomously, end to end, as it sees necessary — including Claude+Codex consults to settle math-level disputes.
3. **`/compact` pauses.** Pause for `/compact` between steps, and within Step 5 between Phase C batches. The pause is context hygiene, not a review gate — resume and continue per this file.
4. **User escalation policy.** The user is involved ONLY when (a) a fatal flaw survives the defense pass (any FIND_STANDS verdict), or (b) a change is conceptual — it alters what the program claims, not how a claim is checked. Everything else (tolerances, wrong constants, tautologies, which stage verifies what, batch granularity, parallelism caps) is Claude+Codex's to resolve, with written consult records.
5. **Adjudication shape — SIGNED OFF.** Phase C adjudication = structured Claude+Codex consult (`codex exec -s read-only`, written record under `redteam_adversarial/verdicts/`), user as final gate on every FIND_STANDS and anything conceptual. This resolves the deployment doc's §7 open question.
6. **Standing rules carried over from the red-team:** 10-minute script timeout (`timeout 600`; exit 124 = failure, reformulate the math, never raise the cap); ≤2 concurrent `math -script` (seat-gate only `.wl`-running sessions); YAML/markdown for anything an agent reads or writes, JSON only machine-to-machine; agent prompts rendered under the project root (sub-agents cannot read `/tmp/`); clean agent per review; offload broad reads to sub-agents to protect orchestrator context; squash small follow-up fixes into the prior commit.
7. **Frozen:** `docs/adversarial_audit_directive.md` is calibrated and must not be edited.

---

## Step 0 — Commit the design baseline

**Status: ✅ DONE (2026-06-10)**

Commit the revised deployment doc, the paper-integration directive, and this execution plan. Exit: clean tree, all three docs on `master`. Landed as the "design baseline" commit.

---

## Step 1 — Paper integration (Codex applies, Claude reviews)

**Status: ✅ DONE (2026-06-10)** — Codex applied all three items (session `019eb50d-7f6d-7d42-a964-4ec7bf1a55ed`, log `research/pde_ledger/redteam_adversarial/codex_logs/step1_iter1.txt`); both builds exit 0 (orchestrator independent re-run; `pdflatex` — `latexmk` not installed, README-documented fallback); clean review agent verified all 14 acceptance criteria PASS, zero edits under `stages/`, two LOW cosmetic notes only (unedited adjacent table vocabulary; pre-existing path-convention mix in provenance summary). The canonical-source fix that gates Phase B is landed.

- Entry: Step 0 committed.
- Hand `docs/adversarial_audit_paper_integration_directive.md` to Codex. Three items: (1) resolve the canonical-source contradiction (two-layer structure; `notes/stages/` named decisive audit layer) — **this gates Phase B**; (2) falsification-stack disclosure table in the reader verification summary; (3) audit-status vocabulary defined next to the claim-status firewall, values generated never hand-tagged. Out-of-scope list in the directive is binding (no appendix stubs, no prose patching, no stage-card edits).
- Codex iterates until both paper builds (reader + archive) compile clean.
- Claude review: acceptance criteria in the directive, checked item by item by a clean review agent + orchestrator spot-read.
- Exit: builds clean, review clean, committed.

## Step 2 — Build the skill

**Status: PENDING**

- Entry: Step 1 committed.
- New skill parallel to `.claude/skills/redteam-audit/`, mirroring its structure: config (`.adversarial-audit-config.yaml` or an `adversarial:` section of `.redteam-config.yaml`), `research/pde_ledger/redteam_adversarial/MANIFEST.yaml` with the per-candidate state machine from deployment doc §5, report tree (`fit_insertion_points.yaml`, `benchmarks.yaml`, `provenance/`, `reports/`, `defenses/`, `verdicts/`, `BATCHES.md`), `query_graph.py` wrappers, Codex invocation reusing `redteam-audit/lib/` patterns (flock-serialized MANIFEST writes included).
- Split per the contract: Claude writes the directive/spec for the skill (requirements + acceptance criteria); Codex designs and writes the runtime code; Claude reviews.
- Exit: skill scaffolding committed; a dry-run invocation against 2–3 known stages exercises the A→B→C plumbing end to end without producing binding verdicts.

## Step 3 — Phase A: fit-insertion-point scan

**Status: PENDING**

- Entry: Step 2 committed.
- Run the four blind modalities from deployment doc §3.1 (numeric-literal, claim-label, graph, existing-provenance cross-check vs `CHECKPOINT_CONSTANT_PROVENANCE.md` + pass-2 reconciliation outputs), union, then the completeness critic; anything the critic surfaces becomes a fifth scan.
- Exit: `fit_insertion_points.yaml` committed, every candidate anchored to stage + file:line; critic pass recorded. Then a Claude+Codex consult decides batch granularity (topical vs dependency-chain clusters via `query_graph.py neighborhood`) and the Phase C parallelism cap — recorded in `BATCHES.md`, escalated only if conceptual.
- Byproduct to surface to the user (non-blocking): the candidate list seeded with parameter values is the start of the long-requested **ansatz-value catalog** (`[[project-ansatz-value-catalog]]`; first entry γ₀ = (1+r_c)/9).

## Step 4 — Phase B: provenance ledger + benchmarks

**Status: PENDING**

- Entry: Step 3 committed AND Step 1's canonical-source fix landed (hard prerequisite — every genealogy entry cites a specific layer).
- Build per-parameter genealogies under `redteam_adversarial/provenance/` by traversing the atlas graph; ingest checkpoint provenance entries as high-confidence source-cited seeds (not exempt from the standard). Binding rules: no attribution without opening `notes/`; a self-contradictory origin claim is itself a finding; graph gaps logged as `graph_gap` and fixed or flagged per deployment doc §7.
- Build `benchmarks.yaml`: one cited literature value per claimed external match (web lookup / textbook / CODATA recorded per entry); never adjudicate a match from model memory.
- Exit: ledger + benchmarks committed; every Phase A candidate has a complete provenance slice or a logged gap.

## Step 5 — Phase C: adversarial deployment (batched)

**Status: PENDING**

- Entry: Step 4 committed.
- Per batch (membership from the Step 3 decision): clean-context adversarial Claude agents, one per candidate, each given the frozen directive + primary source + provenance slice + benchmark entries + graph context. Binary verdict per the directive's required-output format (NO verdicts must list what was checked).
- YES → Codex defense pass (sessions resumed **per-parameter**, not per-stage) → adjudication consult per the signed-off shape → FIND_STANDS / FIND_FAILS / PARTIAL. **Every FIND_STANDS goes to the user — hard stop on that candidate's dependency cone until adjudicated.** NO → checklist logged, archived.
- `/compact` between batches; MANIFEST + `BATCHES.md` updated at each batch close so resume is deterministic.
- Exit: every candidate carries a verdict.

## Step 6 — Close-out (branches on the verdicts)

**Status: PENDING**

- **No surviving fatal flaws:** Codex generates the two archive appendices FROM the real artifacts (fit-insertion index from `fit_insertion_points.yaml`; parameter-provenance appendix from `provenance/`) — generated, never hand-authored; post-audit prose patching per the Codex review's item 6 (demote overclaims, replace unearned "derived"/"forced" labels, add benchmark/provenance citations) — flag conceptual demotions to the user; falsification-stack table updated to show layer 2 complete; builds clean; Claude review; commit. Deliver the completed ansatz-value catalog to the user.
- **Surviving fatal flaws:** user adjudication → program revision → reconverge the algebra → re-run the affected slice of BOTH audit layers (red-team + adversarial) over the flaw's dependency cone → back to Step 5 for the affected candidates.
- Exit: paper and ledgers agree with the verdict record; trackers synced.

## Step 7 — Handoff to layer 3

**Status: PENDING**

Layer 2 passed = the program is cleared for Stage-1 branch realization (`docs/branch_realization_execution_plan.md`, its own plan — not driven by this file). The methodology paper remains gated on at least one Stage-1 result.

---

## Independent thread (interleave anywhere)

Stem-keyed numbering reconciliation (the one remaining post-253 follow-up, Phase 1 done) neither blocks nor is blocked by any step here.

## References

- `docs/adversarial_audit_deployment.md` — design doc (phases, rationale, lessons)
- `docs/adversarial_audit_directive.md` — per-pass adversarial directive (FROZEN)
- `docs/adversarial_audit_paper_integration_directive.md` — Step 1's Codex directive
- `.claude/skills/redteam-audit/` — structural template for Step 2
- Memory anchors: `[[project-adversarial-audit]]` (state), `[[feedback-codex-is-fix-applier]]`, `[[feedback-claude-reviews-codex-codes]]`, `[[feedback-never-alter-calibrated-process]]`, `[[feedback-claude-codex-resolve-math]]`
