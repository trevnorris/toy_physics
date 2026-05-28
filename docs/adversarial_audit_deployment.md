# Adversarial Audit Deployment Plan

**Status:** Design doc. Not yet a skill. Will be turned into one (parallel to `.claude/skills/redteam-audit/`) **after** the current self-consistency red-team finishes its first end-to-end pass.

**Purpose:** Capture the design decisions made during the 2026-05-27 smoke test so the eventual skill is built with empirical grounding, not from scratch.

**Read this if you are:** a future session that needs to know how the *second* audit layer (adversarial fit-vs-derive, on top of internal-consistency red-team) will be deployed across the 253-stage PDE ledger, and why the deployment is shaped the way it is.

---

## 1. Where this fits in the larger pipeline

The PDE ledger faces a three-layer falsification stack:

1. **Internal-consistency red-team** *(in progress)*. The existing `redteam-audit` skill at `.claude/skills/redteam-audit/`. Checks that scripts run, that SymPy and Mathematica independently agree, that no tautological asserts or hardcoded results sneak past. **Currently at 60/253 stages verified through batch III.2** (see `[[project-moving-throat-verification]]` in auto-memory).
2. **Adversarial fit-vs-derive audit** *(this doc)*. Checks that the analytical claims aren't internally consistent *fits dressed as derivations*. The kind of failure self-consistency cannot detect.
3. **Stage-1 branch realization** *(after both audits land)*. Numerical PDE solve testing whether the realized branch returns the targets target-blind. See `docs/branch_realization_execution_plan.md`.

Each layer catches a different failure mode. This doc is about layer 2.

The directive that operates the audit at the per-pass level already exists: `docs/adversarial_audit_directive.md`. What's missing — and what this doc plans — is the **deployment infrastructure** that decides which stages to audit, in what order, with what context, and how to process findings at scale across 253 stages.

---

## 2. Lessons from the 2026-05-27 smoke test

Smoke test: ran one Claude adversarial agent against the first 10 stages (paper TeX + notes + scripts), under the directive. Output: "NO fatal flaw, but the scope is structurally limited because stages 001-010 are pure scaffolding with no external numerical commitments." Smoke-test output was then deleted to avoid contaminating the active red-team pipeline.

What we learned:

### 2.1 The directive is well-calibrated against over-aggression

The "Caution against inventing flaws" line in the directive did real work. The agent in adversary mode returned a clean NO rather than fabricating findings to look productive. It explicitly cited each red-flag category from the directive and reported empty when the category honestly didn't apply.

**Design implication:** No edits to the directive itself. It's working. Calibration risk lives at the deployment-and-triage layer, not the directive layer.

### 2.2 Adversarial audit on scaffolding-only stages is structurally vacuous

Stages 001-010 introduce coordinates, fields, modal decompositions, Schur-complement reductions — pure mathematical setup. No `g=2`, no `α = 1/137`, no PN coefficient claimed against a literature value. The fit-vs-derive test cannot fire because there's nothing to fit. Every numerical constant in the slice (`6 = l(l+1)|_{l=2}`, `4π`, `1/(2√π)`, etc.) is a pure mathematical identity, not a phenomenological match.

**Design implication:** **Stage-order auditing is wrong.** Walking the 253 stages in order would produce a long string of cheap NOs on scaffolding stages, prove nothing, and burn agent-hours.

### 2.3 The audit has to target *fit-insertion points*

The real adversarial check fires only where:

- an external numerical value enters the program (a literature constant, a measured value, a PN coefficient claimed to match GR),
- *or* a previously-free parameter gets fixed to a numerical value,
- *or* a stage claims a derivation whose result is something the program is trying to reproduce.

These are the points where fit-dressed-as-derivation could plausibly hide. Everywhere else, the audit is empty.

**Design implication:** the deployment skill needs to identify fit-insertion points first, then audit those.

### 2.4 The defense step is conditional, not unconditional

If the adversarial verdict is NO, there's nothing to defend against — running the defense agent is wasted compute. The pipeline branches:

- **Verdict YES** → defense pass → adjudication step → either find-stands or find-fails.
- **Verdict NO** → log what was checked, archive, move on.

**Design implication:** the skill should fork after the adversarial verdict, not run a defense pass unconditionally.

### 2.5 Value provenance is the missing infrastructure

The fit-vs-derive test at stage N depends critically on knowing where every parameter in stage N's expressions was introduced and how it was constrained. Without that genealogy, each adversarial pass has to reconstruct the parameter history from scratch, which is expensive and error-prone — and will miss things.

`research/pde_ledger/notes/CHECKPOINT_CONSTANT_PROVENANCE.md` already does narrative provenance tracking for **selected checkpoint stages** in the existing red-team. It's not a queryable structure, doesn't cover all stages, and is keyed to the internal-consistency audit's findings rather than to adversarial deployment.

**Design implication:** the adversarial audit needs its own provenance ledger, queryable, covering every introduction of a numerical-constant or a previously-symbolic parameter, built incrementally during the audit run.

---

## 3. Deployment strategy

### 3.1 Three phases

**Phase A — Fit-insertion-point detection (one pass, non-adversarial).**

A single scanning pass across all 253 stages identifies:

1. Every stage that introduces a numerical literal claimed to be a derived value (not a pure mathematical identity).
2. Every stage that fixes a previously-free parameter to a specific value.
3. Every stage that claims a match against a known physics result (GR, QED, atomic spectrum, etc.).
4. Every stage labeled `derived`, `forced`, or `non-tunable` (the directive specifically warns about these labels — they are claims to check, not facts).

This phase is mechanical, not adversarial. Output: a structured list of *fit-insertion-point candidates*, each anchored to a stage and a file:line citation.

**Phase B — Provenance ledger build.**

A second scanning pass traverses the stage dependency chain (using `graph/fluid_universe_derivation_atlas_graph.yaml` as the dependency backbone) and builds, for each fit-insertion-point candidate, the genealogy:

- where each parameter in the candidate's expressions was first introduced,
- at what stage / line each was assigned a value,
- what constraint set the value (was it pinned by internal consistency, by a published target, or by a free choice?),
- which downstream stages depend on the candidate's claim.

Output: a queryable ledger (YAML) under `research/pde_ledger/redteam_adversarial/provenance/`. One entry per parameter or numerical commitment.

This phase is also non-adversarial — it's bookkeeping. But it's the bookkeeping that makes Phase C sharp.

**Phase C — Adversarial audit at fit-insertion points.**

Adversarial agents are deployed only at fit-insertion-point candidates from Phase A, each provided with:

- The directive (`docs/adversarial_audit_directive.md`).
- The candidate's primary source (TeX, notes, scripts).
- The relevant slice of the provenance ledger (Phase B) covering every parameter in the candidate's expressions.
- The dependency context from the atlas graph.

Each agent produces a binary verdict (YES fatal flaw / NO) per the directive. YES verdicts trigger a Codex-driven defense pass; NO verdicts are logged with their check list.

Defense pass: separate agent (Codex, per skill convention), reads the adversarial report, reads the same source plus provenance, produces DEFEND / CONCEDE / PARTIAL per finding. Human adjudicates.

### 3.2 Use the atlas graph as the navigation backbone

`graph/fluid_universe_derivation_atlas_graph.yaml` already exists with ~500+ nodes and edges encoding the program's derivation chain. The query tool `graph/query_graph.py` supports:

- `search` — find nodes by text
- `show` — node + direct edges
- `neighborhood` — bounded subgraph around a node
- `paths` — short graph paths between two nodes
- `source` — nodes tied to a file path
- `open-gates` — open gates / future-paper markers
- `context` — compact topic packet
- `stats` — graph summary counts

Phase B should consume this graph for dependency traversal rather than re-deriving the dependency structure from raw file scans. The graph is already maintained and validated by the existing tooling; using it directly avoids forking the dependency model.

Concrete use: an agent at stage 87 asks `query_graph.py paths stage_087 stage_023` to confirm whether stage 087's claim genuinely depends on stage 023's parameter assignment, or whether the dependency goes through a different intermediate.

### 3.3 Use the existing checkpoint provenance work as a seed, not a replacement

`CHECKPOINT_CONSTANT_PROVENANCE.md` already contains rigorous narrative provenance for the existing red-team's checkpoint stages. The adversarial provenance ledger should:

- ingest the existing checkpoint entries as authoritative (those stages have been hand-audited),
- treat non-checkpoint stages as still requiring a Phase B traversal,
- format the ledger as queryable YAML rather than prose so Phase C agents can ingest it programmatically.

The existing doc stays — it's the human-readable narrative. The new ledger is the machine-readable backbone.

---

## 4. Pipeline shape

```
Phase A: scan all stages → fit_insertion_points.yaml
                                    ↓
Phase B: traverse atlas graph → provenance/<param>.yaml (one per parameter)
                                    ↓
Phase C: for each fit-insertion-point candidate:
            ↓
       adversarial agent (Claude) + directive + source + provenance slice
            ↓
       verdict: YES or NO
            ↓
       ┌────────── YES ──────────┐         ┌──── NO ────┐
       ↓                          ↓         ↓             ↓
       defense agent (Codex)                log checklist, archive
            ↓
       findings adjudicated by human
            ↓
       FIND_STANDS  /  FIND_FAILS  /  PARTIAL
            ↓
       (FIND_STANDS → blocked; FIND_FAILS → logged with provenance for downstream)
```

Sequential batches, with human gate between batches — same shape as the existing red-team skill. Never roll forward autonomously.

---

## 5. Skill design sketch

When the time comes to turn this into a skill (after the current red-team finishes), the skill should mirror `redteam-audit`'s structure:

- **Config:** `.adversarial-audit-config.yaml` in the project root (or extend the existing `.redteam-config.yaml` with an `adversarial:` section).
- **Manifest:** `redteam_adversarial/MANIFEST.yaml` tracking per-fit-insertion-point state (`pending / scanned / provenance_built / audit_pending / audited / defense_pending / defended / adjudicating / verdict_logged / blocked`).
- **Batches:** group fit-insertion-points by topical or dependency clusters, not by stage number. Use `graph/query_graph.py neighborhood` to define batch membership.
- **Human gate:** between batches, same as red-team skill.
- **Agent invocation pattern:** Claude for adversarial (clean context per fit-insertion-point), Codex for defense (resumed session per finding so it accumulates defense context for the same parameter across multiple stages).
- **Reports under:** `research/pde_ledger/redteam_adversarial/`
  - `fit_insertion_points.yaml` — Phase A output
  - `provenance/` — Phase B per-parameter ledger
  - `reports/` — Phase C adversarial reports
  - `defenses/` — Phase C defense reports
  - `verdicts/` — adjudication outcomes
  - `BATCHES.md` — human-readable progress
- **Tooling:** wraps `graph/query_graph.py` for dependency lookups; reuses Codex invocation pattern from `redteam-audit/lib/`.

Total infrastructure: roughly half the complexity of `redteam-audit/`, because Phase A and B are one-shot bookkeeping rather than per-unit fix-loop work.

---

## 6. Sequencing

Strictly:

1. **Finish current red-team first.** All 253 stages through `verified` status under the existing `redteam-audit` skill. The adversarial audit shouldn't run on scripts that haven't been internally validated — finding a "fit dressed as derivation" in a script that's already algebraically broken is wasted effort.
2. **Run the user's planned full second pass** of the current red-team (see `[[project-full-second-pass]]` in auto-memory). The second pass catches anything the first pass missed; the adversarial audit then builds on stable internal-consistency ground.
3. **Phase A scan.** Build `fit_insertion_points.yaml`. One pass.
4. **Phase B provenance ledger build.** Build per-parameter genealogies. One pass.
5. **Phase C adversarial deployment.** Batched, human-gated, runs until every fit-insertion-point candidate has a verdict.
6. **Adjudication.** Human reviews FIND_STANDS verdicts; either accepts as fatal flaw (program revision required) or rejects (defense holds, finding fails).
7. **If no fatal flaws survive adjudication:** the program has passed both falsification layers and is ready for Stage 1 numerical realization (`docs/branch_realization_execution_plan.md`).
8. **If fatal flaws survive:** program revision; reconverge the algebra; re-run the affected portion of the red-team and adversarial layers; back to step 6.

---

## 7. Open questions / decisions deferred

Things this doc deliberately does **not** decide yet, because they're better resolved against the actual provenance ledger structure once Phase B runs:

- **Batch granularity.** Topical clusters vs dependency-chain clusters vs hybrid. Decide after seeing the shape of `fit_insertion_points.yaml`.
- **Parallel-agent count for Phase C.** The red-team skill caps parallelism via `parallel_audit_max`. The adversarial pass will need its own cap because adjudication is human-bottlenecked downstream; running 16 adversarial agents in parallel produces 16 findings the human now has to adjudicate. Tune empirically once Phase A is done.
- **Defense agent session management.** Codex sessions can be resumed across iterations (the red-team skill does this per-unit). For the adversarial audit, sessions might be resumed per-parameter instead of per-stage, since a defense's "is parameter X back-solved" reasoning is the same across every stage that uses X. Decide after Phase B reveals how parameter reuse actually looks across the program.
- **Whether the adversarial agent should also produce a sidecar provenance update.** Originally proposed in the smoke-test discussion; folded into Phase B as a separate dedicated step rather than a per-pass sidecar. Phase B handles it more cleanly because it runs once with full graph context rather than incrementally with partial context.
- **Atlas graph completeness.** The graph covers "most" of the program (per the user). Phase B may surface stages or parameters that aren't yet graph-indexed. Decision rule: if Phase A finds a fit-insertion-point candidate that has no corresponding graph node, log it as `graph_gap` and either (a) add it to the graph as part of Phase B's bookkeeping or (b) flag for human review. Treat graph gaps as data-quality issues to be fixed during Phase B rather than blocked work.

---

## 8. References

**Companion docs:**
- `docs/adversarial_audit_directive.md` — the per-pass directive itself (no edits needed; calibration confirmed by the smoke test)
- `docs/branch_realization_brief.md` — the Stage-1 numerical test gated on both audit layers passing
- `docs/branch_realization_execution_plan.md` — engineering plan for Stage 1
- `docs/methodology_paper_outline.md` — the methodology paper deferred until at least one Stage-1 result is in hand
- `docs/redteam_thoroughness.md` — first-layer red-team thoroughness write-up

**Infrastructure to be reused:**
- `.claude/skills/redteam-audit/` — first-layer red-team skill; provides the manifest / batches / agent-orchestration patterns the adversarial skill will mirror
- `graph/query_graph.py` — atlas-graph query tool; navigation backbone for Phase B and Phase C
- `graph/fluid_universe_derivation_atlas_graph.yaml` — the program's dependency graph
- `research/pde_ledger/notes/CHECKPOINT_CONSTANT_PROVENANCE.md` — existing narrative provenance to ingest in Phase B
- `research/pde_ledger/notes/CHECKPOINT_TRUST_AUDIT.md` — companion trust-audit doc

**Memory anchors (auto-memory):**
- `[[project-moving-throat-verification]]` — current red-team progress
- `[[project-full-second-pass]]` — planned full second red-team pass
- `[[project-analog-framework-goal]]` — project goal frame (operational analog, not ontology claim)
- `[[feedback-sequential-audit-chunks]]` — never roll forward across batches
- `[[feedback-review-agents]]` — clean agent per review
- `[[feedback-no-json-for-llm-io]]` — YAML/markdown only for LLM-readable files
