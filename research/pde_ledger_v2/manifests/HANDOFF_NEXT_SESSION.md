# Integration-test thread — HANDOFF (written 2026-07-25, pre-compact)

⭐ **READ THIS FIRST**, then `docs/development_pipeline.md` (the process — substantially expanded
this session), then `manifests/FANOUT_PLAN.md`. Do NOT read vision docs; they re-confuse.

Branch `ledger-v2-rebuild`. All work below is COMMITTED unless stated.

---
## 1. WHERE WE ARE

**The physics ledger is essentially built** — 44 stage notes, 44 paper cards, 43 SymPy audits,
43 Mathematica audits, 44 saved `.out`. Each stage has UNIT tests (dual-engine, able-to-fail
teeth). Part VII is 2/7; paused threads at §7.

**The integration layer is 4 of 44 manifests** — 030, 031, 032, 043. Still the pilot. The fanout
has NOT started.

Checker `manifests/composite_build.py`, sha `85ad921b028212b0`, **104 self-test fixtures, all
green**. Current composite build (4 stages):

    FAIL: citations=18/21 claims=121/121 unresolved_producers=3 c7_edges=1/18
          closure=stage030,stage031,stage032,stage043 mathematica_outputs_checked=4

**⭐ The build is RED ON PURPOSE and that is a SUCCESS.** User's standing principle, set this
session: *findings are the product; green is not the goal.* The single new FAIL-class finding is a
REAL cross-stage inconsistency (§3). Baseline findings (3× `ABSENT_PRODUCER` into unextracted
stage003, `REGISTER_COVERAGE_PARTIAL`, 2× `LIFECYCLE_PRODUCER_UNEXTRACTED`, `C7_EDGE_UNCOVERED`)
are honest edges, not defects.

Commits this session: `4f60c9c4` (process contract restored) · `2030e344` (Phase 1: honest scope +
general basis) · `5ac5d8f1` (Phase 2: manifests reworked, drift fixtures) · `06b349fb` (pipeline
doc) · plus stage043 (commit pending fidelity leg).

---
## 2. ⛔ THE FANOUT IS BLOCKED — do not launch it before clearing these

A Grok pre-freeze review (the gate `FANOUT_PLAN.md` reserves for exactly this) found the fanout
would have failed. **Each needs a decision, most need work.**

1. **Dimension recovery covers ~16 of 43 scripts** (13 with a `Dim` class, 3 registered bare-tuple).
   ~18 have NO dimension machinery, yet the schema REQUIRES every symbol object to carry
   `dim`/`dim_source`/`dim_source_order`/`dim_source_tuple`, and `UNSUPPORTED` forces exit 1.
   Half the ledger cannot satisfy the schema honestly → extraction agents would invent loci.
   *Partial relief found at stage043:* `symbols` may be EMPTY (schema has no `minItems`), so a stage
   with no dimensional content can declare `symbols: []`. That does NOT help a stage that has real
   quantities but no dimension machinery in its script.
   **Proposed escape hatch (Codex): an explicit `not_recoverable_from_source` status that yields
   visible coverage and PARTIAL, instead of forcing invented loci. NEEDS A DECISION.**
2. **`infer_closed_slice` can never close.** It requires a continuous `1..max` stage range, but
   **stage029 has no audit script and will have no manifest**. So a full-ledger run stays
   permanently open and silently downgrades every closed-ledger FAIL to PARTIAL. Options: a
   symbol-free stage029 stub manifest; permit declared holes; or replace continuity inference with
   an explicit expected-stage set. **NEEDS A DECISION.**
3. **`BARE_TUPLE_DIM_ORDER_BY_SHA256` is a 3-digest allowlist frozen at freeze time.** Any other
   bare-tuple stage needs a checker edit, which the freeze forbids. Pre-register all such digests
   where order has independent evidence, or move the registry out of checker code.
4. **~42 of ~89 issue codes still have no fixture** (4 drift-critical ones were planted this
   session: `CITATION_DRIFT`, `EXPORT_DIGEST_MISMATCH`, `WRONG_FIRST_FAILURE`,
   `RANGE_ENDPOINT_DRIFT`). Remaining notable gaps: `TOKEN_DRIFT`, `SET_DRIFT`,
   `ADJUDICATION_DRIFT`, `WRONG_OWNER_REFERENCE`, `UNBACKED_SUBSTITUTION`,
   `MULTIPLY_CLASSIFIED_ROW`, `UNKNOWN_EXPORT`, `GENESIS_MISSING`, `MISSING_SOURCE`.
5. **Export-completeness is protocol-only, never enforced** — agents will under-export until a
   downstream consumer fails, causing rework cascades.
6. **`ABSENT_PRODUCER` is ambiguous** — "not extracted yet" vs "typo/wrong id" are
   indistinguishable in open mode. Expensive misdiagnosis at 41 stages.
7. **Vacuous `cas_equivalence`**: an ownership restatement (`S_gg = S_gg`) satisfies C2 while
   detecting nothing. Agents will think they linked content when they linked ids.
8. **Runtime**: ~23s for 3 stages; budget 1–5+ min at 44 (live dim recovery spawns one subprocess
   per unique script and EXECUTES the audit module; C7 mutators are sequential, 30s timeout each).

---
## 3. THE FIRST REAL CROSS-STAGE FINDING (do not "fix" casually)

`RANGE_ENDPOINT_DRIFT — census (39,49,10) != stage043/count_range (40,49,9)`

stage032 declares knob `Q_E` with `{low:0, high:1}`. stage043's audit counts `Q_E` at BOTH endpoints
as route-ful debt. Exact-once classification gives census `[39,49]` spread 10; stage043's range
claim says `[40,49]` spread 9. Codex's reading: **stage043 is likely right** ("pending debt stays
counted"), so stage032's lifecycle is the probable defect — but this was NOT adjudicated.

The schema has NO `reclassified`/`reconciled` lifecycle action (only introduced/inherited/retired/
discharged), so there is no honest way to express a supersession. Fabricating a `{low:+1,high:0}`
"introduced" event would go green and be a LIE — explicitly forbidden.

**Per the user: catalogue findings now, fix them in a dedicated later phase.** Do not resolve this
as a side-effect of making a build green.

---
## 4. THE PROCESS (this is the part to protect)

`docs/development_pipeline.md` grew 99 → 195 lines this session with four new sections. Read it.
The gauntlet per stage, proven on stage043 and now the standard:

> **draft directive → Codex DESIGN-REVIEW → fold → Codex CONFIRM-PASS → execute → checker + FRESH
> fidelity agent.** Grok is a GATE before any freeze, not a per-manifest reviewer.

**I skipped the design-review and confirm legs for most of this session and it cost real time.**
When finally run on stage043, the design-review caught that my directive would have encoded the
discrete count as a range component, producing `[51,60]` instead of `[40,49]` — an error that would
have passed the checker (internally consistent, just describing the wrong object) and probably a
fidelity audit too.

Hard-won specifics now in the pipeline doc: denylist⇒wrong-architecture; "who supplies the value vs
the alphabet"; specify the INVARIANT not the instance; measure empirical premises before building;
"independent of whom?" before calling a check redundant; write scope bounds INTO the directive;
green self-reports are not proof; verify what could BREAK not just what changed; fixtures must never
couple to production state; agents must NEVER delete directories; waiter-reaping ≠ a dead job.

Two defaults that emerged from stage043 and should apply to every fanout stage:
- **Expected-claim inventory** — enumerate from the sources what the stage must carry BEFORE
  writing, else omission passes silently (an executor could export one claim, attach all teeth to
  it, and stay green).
- **Negative control for path coverage** — a passing run emits no marker, so "the check ran" must be
  proven by deliberately tripping it (stage043 did this for the range path).

---
## 5. WHAT WAS PROVEN (don't re-litigate)

- **C7 anti-rederivation has teeth on real stages.** With `C7_FACET` injected the honest mutator
  emits `PASS_F0_MOUTH_VALUE_EVALUATED` (a real tooth in stage031's 50); `--decorative` emits `PASS`
  and is caught as `DECORATIVE_DEPENDENCY`. Note: running the mutator WITHOUT the env var emits
  `PASS` legitimately — that is not a failure, it is the unmutated baseline.
- **`DIMENSIONAL_CONSISTENCY` (formerly C4) is deliberately scoped** to manifest-internal algebra.
  Per-stage dimensional correctness belongs to the dual-engine unit audits. Locus shopping,
  composing dimensions from real bindings, the GL(3,Z) re-gauge, and manifest-supplied anchor
  digests are KNOWN AND ACCEPTED — they require a manifest misrepresenting its own source, and the
  operative risk is DRIFT. **Four rounds were spent here; do not reopen it.** Rationale is in the
  `2030e344` commit message.
- **The general dimensional basis works**: stage038's 4 axes `[M,L,T,E-charge]` and stage042's
  `charge=(1/2,3/2,-1)` over `[stiffness,length,time]` both express and check; cross-basis
  comparison is blocked on both code paths.
- **Codex is the coder** (`docs/development_pipeline.md` §Roles). A previous session's claim that
  "the Codex CLI stalled twice" was VERIFIED FALSE (all 27 runs that day ended exit=0) and had been
  used to justify promoting agents to coders. Agents review; they never write code.

---
## 6. IMMEDIATE NEXT STEPS

1. Commit stage043 once its fidelity leg reports (running at handoff time).
2. **Decide blockers §2.1 and §2.2** (recovery escape hatch; stage029/closed-slice). These gate the
   fanout and are user decisions.
3. Clear §2.3–§2.7 in ONE serialized pre-fanout checker round, then re-freeze and record the sha.
4. Optionally extract 006 (Codex's original slice-A recommendation) as a second no-C7 shakedown.
5. **Then the 44-stage fanout** — user opt-in, Workflow tool, per `FANOUT_PLAN.md`, Codex as coder,
   checker FROZEN, per-stage fidelity agent, export-all + `.out` citation + real loci + expected-
   claim inventory as defaults. Expect it to be RED and full of findings; that is the deliverable.
6. Then the ledger-wide integration report → catalogue every finding → a dedicated remediation phase.

---
## 7. ⏸ PAUSED (do not lose)

- **044-v2** — redo stage 044 with a dynamical-Σ sleeve (un-freeze `S_hold`);
  `notes/stage044_v2_unfreeze_prep.md`. ⚠ its v1 reference manifest was in `_scratch/` and was
  DESTROYED by an agent this session (gitignored, unrecoverable) — re-extract from sources.
- **045** — VII-2b non-variational drain/return block; `notes/stage045_nonvariational_block_prep.md`.
- **Long-term (user's goal):** generalize this whole apparatus into a portable pipeline and test it
  on published arXiv papers to find gaps current review misses. Decided: generalize AFTER the ledger
  is done — the ledger is the proving ground. When that starts, validate first against papers with
  KNOWN errata as positive controls; a pipeline that cannot rediscover a known error cannot be
  trusted to report a novel one.
