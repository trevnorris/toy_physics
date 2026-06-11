# Adversarial Audit Deployment Plan

**Status:** Design doc, revised 2026-06-10; campaign AUTHORIZED same day. The original gate — finish the internal-consistency red-team first — is **satisfied**: both end-to-end passes completed 253/253 on 2026-06-10. Execution is driven step-by-step by `docs/adversarial_audit_execution_plan.md` (operating contract, step gates, resume protocol); this doc remains the design reference.

**Purpose:** Capture the design decisions made during the 2026-05-27 smoke test so the eventual skill is built with empirical grounding, not from scratch. Revised after the red-team completed to fold in what the two finished passes taught us (§2.6) and to close several previously-deferred design gaps (§3.1 Phase A completeness, benchmark sourcing in Phase C, adjudication shape in §7).

**Read this if you are:** a future session that needs to know how the *second* audit layer (adversarial fit-vs-derive, on top of internal-consistency red-team) will be deployed across the 253-stage PDE ledger, and why the deployment is shaped the way it is.

---

## 1. Where this fits in the larger pipeline

The PDE ledger faces a three-layer falsification stack:

1. **Internal-consistency red-team** *(complete)*. The existing `redteam-audit` skill at `.claude/skills/redteam-audit/`. Checks that scripts run, that SymPy and Mathematica independently agree, that no tautological asserts or hardcoded results sneak past. **Both end-to-end passes finished 2026-06-10: first pass 253/253 plus the full second pass 253/253** (per-batch ledger in `research/pde_ledger/redteam/REMEDIATION_HANDOFF.md`; see `[[project-moving-throat-verification]]` in auto-memory).
2. **Adversarial fit-vs-derive audit** *(this doc)*. Checks that the analytical claims aren't internally consistent *fits dressed as derivations*. The kind of failure self-consistency cannot detect.
3. **Stage-1 branch realization** *(after both audits land)*. Numerical PDE solve testing whether the realized branch returns the targets target-blind. See `docs/branch_realization_execution_plan.md`.

Each layer catches a different failure mode. This doc is about layer 2.

**Direct evidence layer 2 is needed:** on 2026-06-10 a Claude+Codex provenance consult caught a paper-card overclaim in `research/pde_ledger/paper/stages/stage_100`/`101.tex` (fixed in commit `b052471`) that **both** complete red-team passes had walked over without firing. Internal consistency was never the problem with that claim — the algebra checked out — so layer 1 structurally could not see it. That is exactly the failure class this audit targets: claims and status badges, not algebra.

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

## 2.6 Lessons from the completed red-team (added 2026-06-10)

The smoke test predates the bulk of the red-team. Three things the two completed passes taught us that this design must absorb:

### 2.6.1 Pass 2 produced a ready-made Phase A seed

The second pass ran an exhaustive script→doc value-reconciliation augmentation (`research/pde_ledger/redteam/pass2/RECONCILIATION_AUGMENTATION.md`) alongside every per-stage audit. Its outputs are, in effect, a per-stage inventory of every numerical value and where it appears across script/notes/TeX. Phase A should ingest this rather than re-scanning cold — it both seeds the candidate list and serves as a cross-check on whatever the fresh scan finds.

### 2.6.2 Mis-attributed provenance is a real, observed failure — not a hypothetical

The χ_Q = 1 provenance episode (fixed in `b052471`): the value's stated origin ("upstream, stage 097") was self-contradictory, and the true derivation point was stage 105. The error survived both red-team passes because nothing in the internal-consistency layer queries provenance. Two consequences for this design:

- The per-parameter genealogy of Phase B is load-bearing, not nice-to-have. It is precisely the structure that would have caught this immediately.
- **Phase B agent rule: never commit an attribution without opening the per-stage `notes/` source.** In the numbering-escalation consults, the orchestrator made attribution calls from cards/scripts alone and was overturned every time by an agent that had read the notes. The decisive evidence for "where was this value really introduced" lives in `notes/stages/` more often than anywhere else.

### 2.6.3 Adversarial disagreement between engines works

The escalation-resolution consults (`e6c8d34`, `b052471`) ran Claude and Codex as independent analysts on the same disputed question and let the disagreement surface the truth. Codex disputed three orchestrator LEAVE calls and was right each time. This is evidence for the attacker/defender split in Phase C — and for using the same structured-consult shape at adjudication time (§7).

---

## 3. Deployment strategy

### 3.1 Three phases

**Phase A — Fit-insertion-point detection (multi-modal scan, non-adversarial).**

A scanning phase across all 253 stages identifies:

1. Every stage that introduces a numerical literal claimed to be a derived value (not a pure mathematical identity).
2. Every stage that fixes a previously-free parameter to a specific value.
3. Every stage that claims a match against a known physics result (GR, QED, atomic spectrum, etc.).
4. Every stage labeled `derived`, `forced`, or `non-tunable` (the directive specifically warns about these labels — they are claims to check, not facts).

This phase is mechanical, not adversarial. Output: a structured list of *fit-insertion-point candidates*, each anchored to a stage and a file:line citation.

**Phase A false negatives are the biggest leak in the whole design.** If the scan misses a fit-insertion point, no adversary ever looks there — a candidate dropped here is silently exempt from the entire layer. A single mechanical pass is not enough. Phase A must be a multi-modal sweep whose modalities are blind to each other, unioned at the end:

1. **Numeric-literal scan** — every non-identity numerical literal in TeX/notes/scripts.
2. **Claim-label scan** — every occurrence of `derived` / `forced` / `non-tunable` / `exact` / `matches` and cognates.
3. **Graph scan** — `query_graph.py` queries for nodes carrying numerical commitments or external-match claims.
4. **Existing-provenance cross-check** — every entry in `CHECKPOINT_CONSTANT_PROVENANCE.md` and every value in the pass-2 reconciliation outputs (§2.6.1) must map to a candidate or be explicitly classified as a pure mathematical identity.

After the union, a **completeness critic** pass asks the one question that matters: "which stage, value, or claim-class did no modality cover?" Anything it surfaces becomes a fifth scan. Only then is `fit_insertion_points.yaml` considered built.

**Phase B — Provenance ledger build.**

A second scanning pass traverses the stage dependency chain (using `graph/fluid_universe_derivation_atlas_graph.yaml` as the dependency backbone) and builds, for each fit-insertion-point candidate, the genealogy:

- where each parameter in the candidate's expressions was first introduced,
- at what stage / line each was assigned a value,
- what constraint set the value (was it pinned by internal consistency, by a published target, or by a free choice?),
- which downstream stages depend on the candidate's claim.

Output: a queryable ledger (YAML) under `research/pde_ledger/redteam_adversarial/provenance/`. One entry per parameter or numerical commitment.

This phase is also non-adversarial — it's bookkeeping. But it's the bookkeeping that makes Phase C sharp.

Two binding rules for Phase B agents, both bought with real failures (§2.6.2):

- **No attribution without the notes.** Every "introduced at stage N" claim must cite the per-stage `notes/` source, not just cards or scripts. Cards and scripts have repeatedly carried stale or mis-attributed origin claims that only the notes resolve.
- **Internal contradiction in a provenance claim is itself a finding.** If a stage's stated origin for a value is self-contradictory (the χ_Q pattern), log it — don't pick the more plausible reading and move on.

**Phase B prerequisite — resolved paper source hierarchy.** As of 2026-06-10 the paper contradicts itself about which TeX layer is canonical: `paper/README.md:43,46-51` calls `paper/stages/` "template inventory; not the active narrative source," while `paper/appendices/reader_provenance_summary.tex:38-42` calls it the canonical narrative source. The verified ground truth is a two-layer structure — the archive stage appendices carry inline theorem-block narrative *and* `\input` all 253 stage cards from `stages/` (which both red-team passes audited as live content, e.g. the `b052471` fix landed in `stages/stage_100`/`101.tex`). Every Phase B genealogy entry must cite a specific layer, so this contradiction must be fixed before Phase B starts — see `docs/adversarial_audit_paper_integration_directive.md`, item 1.

**Phase C — Adversarial audit at fit-insertion points.**

Adversarial agents are deployed only at fit-insertion-point candidates from Phase A, each provided with:

- The directive (`docs/adversarial_audit_directive.md`).
- The candidate's primary source (TeX, notes, scripts).
- The relevant slice of the provenance ledger (Phase B) covering every parameter in the candidate's expressions.
- The dependency context from the atlas graph.
- The sourced benchmark entry for any claimed external match (see below).

Each agent produces a binary verdict (YES fatal flaw / NO) per the directive. YES verdicts trigger a Codex-driven defense pass; NO verdicts are logged with their check list.

**Benchmark sourcing.** The directive tells the agent to "independently recall or compute" the literature value. Recall is reliable for g = 2 and α, but not for higher-order PN coefficients or obscure spectroscopic values — and a mis-recalled benchmark produces both false YESes and false NOs. The deployment layer therefore maintains a small curated benchmark file (`redteam_adversarial/benchmarks.yaml`): one entry per claimed external match, with the literature value, its source citation, and how it was obtained (web lookup at build time, textbook, CODATA). Phase C agents check against the sourced entry; "matches GR" is never adjudicated against model memory alone. The directive itself stays unedited — this is deployment-layer infrastructure, exactly where §2.1 said the calibration risk lives.

Defense pass: separate agent (Codex, per skill convention), reads the adversarial report, reads the same source plus provenance, produces DEFEND / CONCEDE / PARTIAL per finding. Adjudication shape is an open decision for the user (§7) — the original sketch said "human adjudicates," but the program's calibrated convention routes math-level disputes through a Claude+Codex structured consult first, with the human as final gate on anything conceptual.

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

- ingest the existing checkpoint entries as high-confidence *seed* entries, normalized and source-cited under the Phase B schema (those stages have been hand-audited, but no value is exempt from the parameter-value provenance standard — including the notes-citation rule),
- treat non-checkpoint stages as still requiring a Phase B traversal,
- format the ledger as queryable YAML rather than prose so Phase C agents can ingest it programmatically.

The existing doc stays — it's the human-readable narrative. The new ledger is the machine-readable backbone.

**Three kinds of provenance — keep them distinct** (taxonomy from the 2026-06-10 Codex review):

1. **File/source provenance** — which file a piece of text or a script lives in and descends from. Already tracked by `notes/STAGE_PROVENANCE_INDEX.md` (stage → TeX / note source / audit scripts) and the paper's reader provenance appendix.
2. **Result/stage provenance** — which stage a result was established at. Already tracked by the stage cards and `CHECKPOINT_CONSTANT_PROVENANCE.md`.
3. **Parameter-value provenance** — where a numerical value or parameter assignment was *introduced and what constrained it*. **This is the missing layer**, and the only one the fit-vs-derive test actually needs.

The Phase B ledger is layer 3 exclusively. "Provenance" in this doc always means parameter-value genealogy — never file lineage, never stage attribution alone. In particular, do not conflate the new ledger with `STAGE_PROVENANCE_INDEX.md`, which despite the name is a layer-1 file-mapping index.

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
       findings adjudicated (shape: §7)
            ↓
       FIND_STANDS  /  FIND_FAILS  /  PARTIAL
            ↓
       (FIND_STANDS → blocked; FIND_FAILS → logged with provenance for downstream)
```

Sequential batches with a `/compact` pause between them (context hygiene). Under the auto-mode contract the user authorized 2026-06-10 (`docs/adversarial_audit_execution_plan.md`), the between-batch pause is NOT a user review gate: Claude proceeds autonomously, and the user is pulled in only on a FIND_STANDS verdict or a conceptual change. This supersedes the original "human gate between batches" sketch for this campaign (red-team convention unchanged elsewhere).

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
  - `benchmarks.yaml` — sourced literature values for every claimed external match (§3.1 Phase C)
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

1. ~~**Finish current red-team first.**~~ ✅ **DONE 2026-06-10** — all 253 stages `verified` under the existing `redteam-audit` skill (first pass).
2. ~~**Run the user's planned full second pass.**~~ ✅ **DONE 2026-06-10** — second pass 253/253, isolated under `redteam/pass2/`, including the script→doc value-reconciliation augmentation (which now feeds Phase A, §2.6.1).
3. **Phase A scan.** Build `fit_insertion_points.yaml`. Multi-modal sweep + completeness critic (§3.1), seeded from the pass-2 reconciliation outputs.
4. **Phase B provenance ledger build.** Build per-parameter genealogies. One pass. **Prerequisite:** the paper's canonical-source contradiction must be resolved first (§3.1 Phase B; `docs/adversarial_audit_paper_integration_directive.md` item 1).
5. **Phase C adversarial deployment.** Batched, human-gated, runs until every fit-insertion-point candidate has a verdict. (Historical note: an earlier draft gated Phase C on a dual-engine backfill of stages 121/122/123 — that gap was already closed by the 2026-06-01 retro-sweep, 3/3 verified with independent-route `.wl` scripts; see `notes/STAGE_VERIFICATION_COVERAGE.md`. The only legitimately single-engine stages are the 11 status-only stages that compute nothing a `.wl` could check; Phase C agents should treat single-engine coverage there as expected, not as a red flag.)
6. **Adjudication.** Human reviews FIND_STANDS verdicts; either accepts as fatal flaw (program revision required) or rejects (defense holds, finding fails).
7. **If no fatal flaws survive adjudication:** the program has passed both falsification layers and is ready for Stage 1 numerical realization (`docs/branch_realization_execution_plan.md`).
8. **If fatal flaws survive:** program revision; reconverge the algebra; re-run the affected portion of the red-team and adversarial layers; back to step 6.

---

## 7. Open questions / decisions deferred

Things this doc deliberately does **not** decide yet, because they're better resolved against the actual provenance ledger structure once Phase B runs:

- **Batch granularity.** Topical clusters vs dependency-chain clusters vs hybrid. Decide after seeing the shape of `fit_insertion_points.yaml`.
- **Parallel-agent count for Phase C.** The red-team skill caps parallelism via `parallel_audit_max`. The adversarial pass will need its own cap because adjudication is human-bottlenecked downstream; running 16 adversarial agents in parallel produces 16 findings the human now has to adjudicate. Tune empirically once Phase A is done.
- **Defense agent session management.** *(Recommendation now firm: per-parameter.)* Codex sessions can be resumed across iterations (the red-team skill does this per-unit). For the adversarial audit, sessions should be resumed per-parameter rather than per-stage: values like χ_Q recur across many stages, and the "is parameter X back-solved" reasoning is identical at each occurrence. Phase B's ledger (one entry per parameter) makes this the natural unit. Confirm against the actual reuse pattern Phase B reveals, but per-parameter is the default.
- **Adjudication shape — ✅ RESOLVED (user sign-off 2026-06-10).** The original sketch had the human adjudicate every DEFEND/CONCEDE/PARTIAL outcome directly. The program's calibrated convention, however, routes math-level disputes through a Claude+Codex structured consult first, escalating to the user only when the outcome is conceptual — and that consult format has a strong empirical record (§2.6.3: three orchestrator calls disputed, three correctly overturned). **Signed-off shape:** adjudication = structured Claude+Codex consult (same shape as the numbering-escalation consults, `codex exec -s read-only`, written record under `redteam_adversarial/verdicts/`), with the user as final gate on every FIND_STANDS and on anything that changes the program's conceptual claims. Codex's 2026-06-10 review independently concurred; the user authorized it alongside the auto-mode operating contract in `docs/adversarial_audit_execution_plan.md`.
- **Whether the adversarial agent should also produce a sidecar provenance update.** Originally proposed in the smoke-test discussion; folded into Phase B as a separate dedicated step rather than a per-pass sidecar. Phase B handles it more cleanly because it runs once with full graph context rather than incrementally with partial context.
- **Atlas graph completeness.** The graph covers "most" of the program (per the user). Phase B may surface stages or parameters that aren't yet graph-indexed. Decision rule: if Phase A finds a fit-insertion-point candidate that has no corresponding graph node, log it as `graph_gap` and either (a) add it to the graph as part of Phase B's bookkeeping or (b) flag for human review. Treat graph gaps as data-quality issues to be fixed during Phase B rather than blocked work.

---

## 8. References

**Companion docs:**
- `docs/adversarial_audit_execution_plan.md` — the ACTIVE operational driver (operating contract, step gates, `/compact` resume protocol)
- `docs/adversarial_audit_directive.md` — the per-pass directive itself (no edits needed; calibration confirmed by the smoke test)
- `docs/adversarial_audit_paper_integration_directive.md` — Codex directive for the paper-side edits (source-hierarchy fix, falsification-stack disclosure) that gate Phase B
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
- `research/pde_ledger/redteam/pass2/RECONCILIATION_AUGMENTATION.md` — pass-2 value-reconciliation augmentation; Phase A seed (§2.6.1)
- `research/pde_ledger/redteam/REMEDIATION_HANDOFF.md` — durable per-batch record of both completed red-team passes
- `research/pde_ledger/notes/STAGE_PROVENANCE_INDEX.md` — *file-mapping* index only; see §3.3 naming-collision note

**Memory anchors (auto-memory):**
- `[[project-moving-throat-verification]]` — red-team record (both passes complete 2026-06-10)
- `[[project-full-second-pass]]` — planned full second red-team pass
- `[[project-analog-framework-goal]]` — project goal frame (operational analog, not ontology claim)
- `[[feedback-sequential-audit-chunks]]` — never roll forward across batches
- `[[feedback-review-agents]]` — clean agent per review
- `[[feedback-no-json-for-llm-io]]` — YAML/markdown only for LLM-readable files
