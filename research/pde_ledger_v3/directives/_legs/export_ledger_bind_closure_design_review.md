# Decision-list review — the bind-closure LEDGER design

You are one of two independent legs reviewing an orchestrator-written **decision list** (rule 7: two legs before
any builder). It redesigns what the plain-git `*_exports.py` LEDGER carries and how it is stored. It synthesizes
a prior two-engine deliberation, but it makes **specific new decisions** (a topology, a migration, a guard) that
the deliberation did not — those are what you must pressure-test. ⛔ Do not defer to the design's own claims;
form your view from the sources first, then check each decision. Cite exact lines for every finding.

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/export_ledger_bind_closure_design.md`

## Read first, form your own view
- The two deliberation outputs this synthesizes: `directives/_legs/export_ledger_minimum_deliberation.md` (the
  question) and the saved leg outputs `scratchpad`-origin — if not present, re-derive from the sources below.
- The prior decision it reopens: `directives/chain_accumulate_or_generate_decision.md` (2026-08-11, generate vs
  accumulate — read the supersession section and the measurement).
- What it supersedes/keeps: `directives/S11_export_chain_decisions_v2.md` (F1–F10, esp. **F10** and the
  "What F9 does NOT decide" model/step note), `S9_REWRITE_PLAN.md` (**D1** membership, **D3** round-trip, **D5**,
  **D7**).
- The chain reality: `scripts/S11c_b_exports.py` (the accumulated 58.75 MiB file; the `_restore` reviver at
  the top; a `*_term_origins` row; `slab_operator`), `scripts/S11c_b_cross_engine_comparator.py` (line ~1206
  residual; header "not used as a tag stream"), the c1 spec `directives/S11c_c1_SHARED_PHYSICS.md` (§0 c2 scope,
  §3b/§7 c1 consume-set, §5c/§5d regression binds) and the c1 build directive.

## What to check — findings only where the design would MIS-BUILD the chain, UNDER-export (break a later build), fail to shrink (miss the goal), or contradict a live constraint

1. **D1 (bind-closure membership).** Is "export iff a later step binds it, + symbol/`dimension_key`/F9c-route
   closure, membership by bind not label" correct and complete? Does it wrongly **drop** any row a declared later
   step needs (probe c1→c2 specifically: μ_θ, T-substrate, the S11b regression refs, `slab_operator`/
   `coupling_kernel` for c2), or wrongly **keep** a dead one? Is the ~508-row / ~26 MiB estimate defensible?
2. **D2 (generate over a frozen base).** Is freezing the accumulated base through `S11c_b_exports.py` and
   writing own-rows deltas after it sound — does it truly avoid the `BUILD_INPUT_DIGESTS` sha cascade (check how
   digests pin upstream files), and does the `base ⊔ deltas` fold with `F9` collision semantics resolve
   correctly when a delta re-derives a base row? Is the hybrid (accumulated base + own-rows deltas) a real
   `F9`/fold ambiguity, or benign? Does it actually clear the 100 MiB cap for c1?
3. **D3 (manifest + closure guard).** The design claims **D3-round-trip does NOT catch under-export** because
   `_restore` builds `Symbol(...)` from srepr regardless of whether the symbol's row is present — verify that
   against the reviver. Is the 5-step guard (`IMPORT_KEYS`/`PASS_THROUGH`, F9-route assertion, closure
   assertion, bind smoke-test, minimum-mode) actually **closed** — does it miss a reachable under-export through
   F9c routing or `dimension_key` chains? Is **`PASS_THROUGH` declared at write time** the right mechanism, or
   should pass-through be **discovered by the fold's closure at consume time** (which would remove the
   forward-knowledge burden the generate topology exists to remove)? This is the design's own flagged question.
4. **D4 (representation out of scope).** Correct to separate it? The design says the comparator canonicalizes
   `A−B` so the LEDGER need not store expanded — verify against the comparator. Anything that makes representation
   *not* separable from the set decision?
5. **Reconciliation.** Does superseding F10 and reversing F10/N1's accumulation collide with any live consumer
   (the comparator, the WL engine, `extract_knob_inventory`, any importer)? Is anything in the
   supersede/keep list wrong (e.g. keeps something it should strike, or strikes something load-bearing)?
6. **Under- vs over-reach.** Is any decision here actually a **representation** or **naming** decision smuggled
   in? Does the design add physics or an expected value (it must not)? Is the no-manifest **D1-fallback** a
   loophole that lets the bloat return?

## Method + filter
- Derive your own answer to "what must c2 be able to bind, and can the design's fold+guard deliver it" **before**
  trusting the design. Cite lines.
- ⛔ Read-only; modify nothing.
- Report a finding only if it catches a concrete failure (mis-build, under-export, missed goal, contradiction) —
  not style. If a check passes, say so and name what you verified. Rank must-fix above nit.

## Output
Findings ranked most-severe first, each quoting (a) the source line your expectation comes from and (b) the
design line it violates or satisfies. End with a one-line verdict: is the design sound to build against as-is,
or are there must-fix items — and if the topology/pass-through question has a better answer, state it.
