# Integration-test fanout plan — ⛔ SUPERSEDED / CANCELLED (2026-07-25)

> ⛔ **THE PARALLEL FANOUT IS CANCELLED.** Plan of record is now
> **`manifests/LEDGER_WIDE_PLAN.md`** — sequential, stage-by-stage, preceded by a ledger-wide
> dimension-unification pass.
> **Why:** a shared dimension-declaration module makes parallel extraction *actively wrong* (every
> stage writes one file — the `composite_build.py` race, moved from tooling into content), and
> sequential extraction in stage order builds the causal graph bottom-up so `ABSENT_PRODUCER`
> means "typo" immediately instead of everything sitting PARTIAL until the end. A Grok pre-freeze
> review also found the fanout would have failed outright: dimension recovery covers ~16 of 43
> scripts while the schema requires it for all, and `infer_closed_slice` can never close because
> stage029 has no audit script.
> Retained below for its still-valid parts: delegation/token-lean reasoning, the Mathematica
> evidence decision, and the export-completeness rationale.

How the remaining ~38 stage manifests get extracted AFTER the pilot, keeping the
checker frozen. Runs only on explicit user opt-in (given for the Workflow approach).

## Delegation architecture (token-lean — the whole point)
The orchestrator (Claude) must NOT populate manifests inline — that absorbs every
11k-line Codex log into its context. Delegate so the logs stay out.

- **Driver = the Workflow tool** (user opted in 2026-07-24; first use). Per-stage
  pipeline, concurrency throttled to ~5 for monitorability, runs in the background so
  agent tool-output stays out of the orchestrator context, live progress via
  `/workflows`.
- **Per stage: CODEX IS THE CODER — there is no second option.** An agent drives Codex
  (busy-polls it), Codex writes the manifest, then the checker runs → the agent returns a
  COMPACT verdict. This is the standing calibrated contract
  (`docs/development_pipeline.md` §Roles): agents are REVIEW instruments and never write
  code; delegating coding to an agent is the same violation as Claude doing it by hand.
  ⚠ Sub-agents are ONE-SHOT — they CANNOT await background notifications; the driving
  agent must busy-poll in a foreground loop with an emphatic "you will NEVER receive a
  notification; deciding to stop and wait is a BUG" rule
  ([[feedback-subagent-marathon-infra]]).

  > ⚠ CORRECTION (2026-07-24): this section previously offered "agent-as-coder" as a
  > co-equal option (ii) and recommended it, on a RELIABILITY NOTE claiming the Codex CLI
  > stalled twice mid-`stage031`. **That claim did not survive verification.** All 27 Codex
  > runs from that day terminate `___CODEX_BUILD_DONE___(exit=0)` — including five stage031
  > extractions (v1→v5) — with no marker-less log anywhere, which is what a real stall
  > leaves. codex-cli is 0.145.0 (current) and a fresh smoke run answered correctly in 20s.
  > The likely cause was a reaped background waiter or a marker match not anchored to
  > `tail -1`, misread as a stall. **If Codex genuinely breaks, we FIX Codex and HALT to the
  > user — we never promote an agent to coder.**
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
- ✅ RESOLVED (2026-07-24): the user chose **(b)** and it's DONE — a 2-seat-capped batch
  (`_scratch/run_math_batch.sh`) re-ran all 44 `.wl` → `mathematica/out/*.out`, committed
  `7f6d9481` (436K; all 43 stage audits `OVERALL PASS`, midway `CONSISTENT`, zero errors).
  These are the pinnable dual-engine evidence. **Follow-on (fold into #27 evidence pass):**
  update the manifests to CITE the `.out` (path + sha256 digest) alongside the sympy audit;
  a stage's Mathematica leg is then trusted from the saved `.out` while its `.wl` digest is
  unchanged. (`.out` are un-ignored via a `.gitignore` negation; the repo `*.out` rule is
  for LaTeX artifacts.)

## Export-completeness (✅ APPROVED by user 2026-07-24 — implement in #27 `PREFANOUT_ROUND.md`)
The pilot's systematic under-export (031 didn't export `S_gg` → 032 fell to opaque `C_V`;
030 lacked `O_perp` → had to be added) is a CURATION bug: extractors treated `exports` as
a highlight reel, but downstream cites intermediate results + quantity definitions, not just
headlines. **RECOMMENDED FIX — make it mechanical, not editorial: a stage exports EVERY
operative claim AND every ownership (`declare_*`) claim; only retired/superseded/departed
claims stay unexported.** Rationale: in a derivation ledger every result is legitimately
citable (an intermediate lemma is exactly what a downstream stage should cite rather than
re-derive); it eliminates the whole under-export failure class deterministically on the first
pass (vs a brittle demand-driven reconciliation pass); export-list size is free; and the
`NON_EXPORTED_CLAIM` guard survives where it matters (citing a RETIRED/departed claim is still
a real error → still caught). Lands as: (1) protocol rule-7 clarification ("export EVERY
operative + ownership claim; only retired/departed stay internal"); (2) retroactive complete
re-export of 030/031/032 in #27 → 031 exports `S_gg` → 032 drops the opaque `C_V`; (3) fanout
default so extractors export-complete from the start (zero under-export stops).

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
