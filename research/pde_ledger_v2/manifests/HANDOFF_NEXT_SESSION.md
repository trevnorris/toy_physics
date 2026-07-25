# Integration-test thread — HANDOFF (written 2026-07-25, pre-compact)

⭐ **READ THIS FIRST**, then **`manifests/LEDGER_WIDE_PLAN.md`** (the PLAN OF RECORD — sequential
stage-by-stage work preceded by a ledger-wide dimension-unification pass), then
`docs/development_pipeline.md` (the process — substantially expanded this session).
⛔ `manifests/FANOUT_PLAN.md` is SUPERSEDED — the parallel fanout is CANCELLED.
Do NOT read vision docs; they re-confuse.

**Operating model for the whole run: Claude is a THIN CONDUCTOR.** Agents return ≤10-line VERDICTS
(full reports to disk, read only on a negative verdict); the orchestrator independently re-runs the
acceptance commands and greps for the specific line. Detail is caught by a correctly-scoped
reviewer, not by orchestrator attention — *you are your agents, as long as they are loaded up
correctly*. Holding one long session keeps process continuity STRUCTURAL rather than dependent on
docs surviving compaction. Details in `LEDGER_WIDE_PLAN.md` §1.

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

   ⭐ **USER DECISION (2026-07-25) — FIX THE CORPUS, NOT THE CHECK.** An earlier proposal here was an
   escape-hatch status (`not_recoverable_from_source`). **REJECTED, correctly:** it weakens the check
   because the corpus is inconvenient, when the corpus is what is actually wrong. Thirteen different
   dimension idioms across 43 scripts IS drift; the checker was reporting it accurately.
   **THE PLAN: a single shared dimension-declaration module, imported by every audit script.**
   Had this existed from the start, the whole class would have been caught immediately.
   What it buys beyond the bug fix:
   - **Deletes fragile checker machinery** — `BARE_TUPLE_DIM_ORDER_BY_SHA256`, the AST bare-tuple
     recovery, the live module-execution path, and the per-script order registry all exist ONLY to
     reverse-engineer dims out of inconsistent idioms. One import ⇒ one recovery path.
   - **May un-do the C4 compromise.** Dimensional certification was scoped to consistency-only
     BECAUSE identity→dimension resolved 0/92 symbols — there was no per-quantity binding to look
     up. A module keyed by quantity identity is precisely the missing structure. Re-evaluate the
     scope decision AFTER the module exists (do not assume it restores full certification; do not
     assume it cannot).
   - **stage038 (4 axes: M,L,T,E-charge) and stage042 (stiffness,length,time with fractional
     exponents) become DECLARED rather than discovered.**
   Honest costs / risks, all manageable but none skippable:
   - **Digest cascade:** every manifest pins its script's `source_digest`; touching 43 scripts
     invalidates those pins in all 4 existing manifests → re-pin. (The Mathematica `.out` evidence
     is unaffected — Python changes do not alter it.)
   - **The refactor MUST be behaviour-preserving, and that is testable:** every audit still exits 0
     with an IDENTICAL PASS count. A changed dimension VALUE is a regression, not a refactor.
   - **`ACTIVE_MUTATION` coupling:** every script reads it at module scope and the C7/ablation
     harnesses depend on that; the shared import must not disturb it.
   - **Not every script needs dims.** stage043 legitimately has none (integer-count stage). The
     survey must separate "no dimensional content" from "dimensional content, no machinery".
   - **This is ON THE CRITICAL PATH** — it must precede the fanout, or 41 manifests get extracted
     against the broken state and need redoing.
   Sequence: survey all 43 → design the module (multiple bases, exact rationals) → refactor with
   identical-PASS-count as the gate → re-pin digests → simplify the checker to the single recovery
   path → re-evaluate the C4 scope → fanout.
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
2. ✅ **BOTH §2 DECISIONS ARE MADE (user, 2026-07-25) — execute them, do not re-open:**
   - **§2.1 → build the shared dimension-declaration module and refactor all 43 audit scripts onto
     it.** Do NOT add an escape-hatch status. Start with the SURVEY (which stages have real
     dimensional content vs none), because it sizes everything downstream.
   - **§2.2 → replace continuity inference with an explicit canonical stage inventory** in
     `composite_config.json` (digest-pinned). This also resolves §2.6: a reference to a stage NOT in
     the inventory becomes a hard FAIL (typo), while "in the inventory, not yet extracted" stays
     PARTIAL. Give stage029 a real manifest too, but do not make closedness depend on it.
   - Other agreed opinions: §2.3 move the bare-tuple registry OUT of checker code into a
     digest-pinned data file (mostly moot once §2.1 lands); §2.4 fixture by USAGE — prioritise
     `TOKEN_DRIFT`, `SET_DRIFT`, `ADJUDICATION_DRIFT`, `MULTIPLY_CLASSIFIED_ROW`, not all 42;
     §2.5 ENFORCE export-completeness mechanically; §2.7 emit an ADVISORY (not FAIL) for tautological
     `cas_equivalence` so the report shows how many citations are actually protective.
3. Clear the remaining §2 items in ONE serialized pre-fanout checker round, then re-freeze + record sha.
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
