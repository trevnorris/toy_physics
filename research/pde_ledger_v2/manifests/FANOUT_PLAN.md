# Integration-test fanout plan (BANKED 2026-07-24, user-approved)

How the remaining ~38 stage manifests get extracted AFTER the pilot, keeping the
checker frozen. Runs only on explicit user opt-in (given for the Workflow approach).

## Delegation architecture (token-lean — the whole point)
The orchestrator (Claude) must NOT populate manifests inline — that absorbs every
11k-line Codex log into its context. Delegate so the logs stay out.

- **Driver = the Workflow tool** (user opted in 2026-07-24; first use). Per-stage
  pipeline, concurrency throttled to ~5 for monitorability, runs in the background so
  agent tool-output stays out of the orchestrator context, live progress via
  `/workflows`.
- **Per stage:** an agent extracts the v2.1 manifest against the FROZEN checker → runs
  `composite_build.py` → returns a COMPACT verdict. TWO coder options:
  - **(i) agent-launches-Codex** (Codex is the coder; agent busy-polls it). Preserves the
    Codex-codes contract. ⚠ Sub-agents are ONE-SHOT — they CANNOT await background
    notifications; the agent must busy-poll in a foreground loop with an emphatic "you
    will NEVER receive a notification; deciding to stop and wait is a BUG" rule
    ([[feedback-subagent-marathon-infra]]).
  - **(ii) agent-as-coder** (the agent does the extraction DIRECTLY — read sources, write
    the manifest, run the checker, iterate). Simpler process model, runs entirely
    in-harness (no external Codex process to stall).
  ⚠ RELIABILITY NOTE (2026-07-24): the Codex CLI stalled TWICE mid-`stage031`
  construction (~15 min silent, proc alive ~1.5% CPU = a hung model-API response), forcing
  a pivot to agent-as-coder for that stage. If Codex reliability stays flaky, option (ii)
  is the safer fanout default. Decide at the opt-in checkpoint based on how Codex behaves.
- **Fidelity review = a SEPARATE fresh agent per stage** (never the extractor —
  no self-review). Checks manifest faithfulness vs the sources — the thing the checker
  structurally cannot ([[feedback-review-agents]]). This is how the deferred
  "blind-rebuild" review (task #15) actually discharges. The stage030 pilot proved a
  fidelity agent is REQUIRED (it caught the `V_H`→`O_⊥` operator mislabel the checker
  passed).
- **Grok RESERVED** for the checker + genuine cross-stage integration surprises — NOT
  run per-manifest (would exhaust its usage limit for little gain; the per-manifest gate
  is fidelity-agent + mechanical checker).
- **Checker FROZEN during fanout:** agents REPORT checker/schema gaps, never patch.
  Serialized checker fixes (if any survive the pilot) are the orchestrator's, run
  BETWEEN waves — never by a parallel extraction agent (they would race on the shared
  `composite_build.py`; the pilot proved this).

## Mathematica / evidence — cap does NOT bind extraction
- Extraction does NOT execute Mathematica. A manifest cites the SYMPY audit
  (`engine: sympy`, actually re-run by the checker) + the `.wl` as a digest-pinned
  companion pointer (as stage030 does). So the ≤2 Mathematica-seat cap does NOT bind
  extraction → fanout runs at ~5.
- ⚠ FACT: there are **44 `.wl`, 0 `.out`** — the assumed saved-Mathematica-output does
  NOT exist on disk. DECISION (user, pending): (a) cite the `.wl` digest only [cap ~5,
  simplest, matches stage030]; or (b) a ONE-TIME ≤2-seat-capped batch
  (`math -script X.wl > X.out` ×44, ~1h Mathematica) to produce durable, digest-pinnable
  saved dual-engine evidence [stronger; realizes the original intent; NOT blocking].
  The ≤2 cap binds ONLY if we run (b).

## Sequencing (do NOT reorder)
1. **Pilot FIRST** (030/031/032/006/043 + one real 030→031→032 C7 mutation) — shakes out
   checker/export/citation gaps on ~5 stages so the 38-stage fanout hits far fewer
   stop-and-report snags. IN PROGRESS.
2. **User opt-in checkpoint** — confirm the fanout (+ the Mathematica (a)/(b) choice)
   before launching the Workflow.
3. **Workflow fanout** of ~38 stages (extract → checker; fidelity-agent second stage) →
   full composite build → ledger-wide integration report → triage real FAILs.
4. **C8 + prose/script consistency + wire** `composite_build.py --self-test` + a partial
   composite run into the per-commit flow so integration stays live.

Anchor: [[project-integration-test-manifest-system]]; companion `RE_PILOT_PLAN.md`.
