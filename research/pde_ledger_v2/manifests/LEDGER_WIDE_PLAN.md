# Ledger-wide plan — dimension unification + sequential manifest extraction

> **Supersession notice (2026-07-26): §2.1–§2.2 are historical and superseded by
> `manifests/DIMENSION_REWRITE.md`. Do not run the old survey-first/module-design
> sequence; use that document's measured inventory and per-stage rewrite loop.**

**Status: THE PLAN OF RECORD (user-decided 2026-07-25); the parallel fanout is
CANCELLED.** Read this after `manifests/DIMENSION_REWRITE.md`, alongside
`docs/development_pipeline.md`.

---
## 0. WHY NO FANOUT

1. **A shared dimension-declaration module makes parallel extraction actively wrong**, not merely
   risky — every stage writes to one file. That is the `composite_build.py` race the pilot already
   proved, moved from tooling into content.
2. **Sequential in stage order builds the causal graph bottom-up.** Every citation points at an
   already-extracted stage, so `ABSENT_PRODUCER` means "typo" *immediately*. The fanout's
   alternative was everything sitting PARTIAL until the end — the weakest possible error signal.
3. A frozen checker is required for a fanout; sequentially the checker can keep improving as gaps
   are found, which matters because ~42 issue codes are still unfixtured.

**Cost is not the constraint. Correctness is.** (User, 2026-07-25.) The corresponding rule:
**never adjust the process because the corpus is inconvenient** — an escape hatch for unrecoverable
dimensions was proposed and rejected for exactly this reason.

---
## 1. THE OPERATING MODEL — Claude is a THIN CONDUCTOR

The orchestrator holds the thread for the whole run so process continuity is *structural* rather
than dependent on documents surviving compaction. Docs become the backup, not the mechanism.

**"You are your agents, as long as they're loaded up correctly."** Detail is caught by a
correctly-scoped reviewer, NOT by the orchestrator reading transcripts. Every catch this session —
the `operator: psi` mislabel, the executed dimension forgery, the fanout blockers — came from
dispatch, not from orchestrator attention.

**Rules for staying thin:**
- **Agents return a VERDICT, not a report.** ≤10 lines: verdict, counts, severities, the one-line
  acceptance result. The full report is written to a file on disk. The orchestrator reads that file
  ONLY when the verdict is negative.
- **Codex likewise** — full report to disk; return acceptance lines + exit code.
- **The orchestrator independently re-runs the acceptance commands** and `grep`s for the specific
  line (never `head`/`cat` of a long output). This is what stops the thin model from decaying into
  rubber-stamping, and it is cheap.
- **Verdicts must be checkable, not adjectival.** "FAITHFUL" alone is useless. "FAITHFUL, 0
  findings, build line X, 20/20 teeth verified against TOOTH_ORDER" is checkable.
- **Per-stage detail will live in the Pass-2 progress ledger** (§6), created when Pass 2 begins and
  read as a `tail`. A stage whose
  row cannot show its fidelity leg is NOT done.
- Target ≈30–60 lines of orchestrator context per stage.

⛔ **Every agent prompt carries the file-safety clause:** never delete a directory or any file the
agent did not create; scratch goes in its own subdirectory. A review agent destroyed the shared
gitignored `_scratch/` tree this session — unrecoverable, and it took stage044's v1 reference
manifest with it.

---
## 2. PASS 1 — DIMENSION UNIFICATION (all 44 stages, before any manifest work)

### 2.1 Survey first (read-only, PARALLELIZABLE — safe, no shared writes)
For each of the 43 audit scripts record: does it have **real dimensional content** or none
(stage043 legitimately has none — it prints `DIMENSIONAL_HOMOGENEITY=N/A_INTEGER_COUNT_STAGE`);
which **basis** and axis order; which **quantities** carry dimensions; the current **idiom**; and
whether the `.wl` states dimensions **differently** from the `.py`.
This sizes everything downstream and is the module's requirements document.
Known going in: 13 scripts define a `Dim` class; 3 have registered bare-tuple digests (032/038/042);
~14 have no dimension machinery at all; several more have machinery recovery cannot see (type-alias
`Dim`, or bindings not named `dim_*`/`*_dim`, e.g. stage021's `SOURCED_N0_DIM`). **stage029 has no
script at all.**

### 2.2 Design the module against the EXTREMES first
**Do stage042 and stage038 BEFORE the easy LMT stages.** stage042's basis is
`(stiffness, length, time)` with fractional exponents (`charge=(1/2,3/2,-1)`); stage038 uses FOUR
axes `(M, L, T, E-charge)`. Designing for LMT, converting 37 stages, then discovering 038 breaks the
abstraction is the failure mode. **The extremes are the spec.**
Requirements: arbitrary named axis sets, arbitrary axis count, exact rational exponents (floats are
errors), and — critically — **keyed by quantity identity** so a symbol's dimension is locatable from
its own name/`quantity_id`. That last property is what may un-do the C4 compromise (§5).

### 2.3 Refactor, gated on behaviour preservation
**GATE: every audit script still exits 0 with an IDENTICAL PASS count.** A changed dimension VALUE
is a regression, not a refactor — investigate it, never absorb it.
- ⚠ **`ACTIVE_MUTATION` coupling:** every script reads it at module scope and the C7/ablation
  harnesses depend on that. The shared import must not disturb it.
- Separate "no dimensional content" (leave alone) from "dimensional content, no machinery" (add it).

### 2.4 Re-run Mathematica and regenerate `.out`
Cheap relative to everything else, and it catches dual-engine disagreement while we are already
inside each script.
- ⚠ **≤2 concurrent `math -script`** (2-seat licence). Never SIGKILL a healthy kernel.
- No live batch runner exists: the former scratch runner was destroyed. Follow
  `manifests/DIMENSION_REWRITE.md` §4 for the current per-stage rerun and
  reproducibility procedure; any future batch runner must still enforce the
  two-seat cap.
- New `.out` files ⇒ new digests, which every manifest's `mathematica_output` pins.

### 2.5 The digest cascade — settle it ONCE, before any manifest work
Three digest classes move in this pass: **script**, **shared module**, **`.out`**.
⭐ **The shared module becomes part of every stage's evidence and must itself be digest-pinned.** A
manifest pinning only its script's digest would no longer capture what the script computes. This is
the decisive reason for two passes: interleaving would have every module change invalidate every
already-written manifest. Re-pin the 4 existing manifests (030/031/032/043) at the end of Pass 1.

---
## 3. PASS 1b — CHECKER ROUND (serialized, then re-freeze)

With one uniform dimension source, **DELETE the machinery that existed only to reverse-engineer
dimensions out of inconsistent idioms**: `BARE_TUPLE_DIM_ORDER_BY_SHA256`, the AST bare-tuple
recovery, the live module-execution path, the per-script order registry. One import ⇒ one recovery
path. This is a strict simplification of the most fragile code in the checker.

Also in this round (all user-agreed):
- **Canonical stage inventory** in `composite_config.json`, digest-pinned, replacing
  `infer_closed_slice`'s continuity rule (which can never close, because stage029 has no script).
  This also fixes the `ABSENT_PRODUCER` ambiguity: a stage id NOT in the inventory is a typo →
  hard FAIL; "in the inventory, not yet extracted" → PARTIAL. Give stage029 a real manifest, but do
  not make closedness depend on it.
- **Enforce export-completeness mechanically** (every operative + ownership claim appears in
  `exports`) instead of leaving it a protocol hope.
- **Advisory (not FAIL) for tautological `cas_equivalence`** (`lhs - rhs ≡ 0`), so the report shows
  how many citations actually protect anything. Ownership restatements are legitimate; silent
  non-protection is not.
- **Fixture by USAGE, not exhaustively:** priority `TOKEN_DRIFT`, `SET_DRIFT`, `ADJUDICATION_DRIFT`,
  `MULTIPLY_CLASSIFIED_ROW`. (~42 codes are unfixtured; a cold-path gap is far cheaper than a
  hot-path one.)
Then **re-freeze and record the sha**, with a Grok gate before the freeze (§7).

---
## 4. PASS 2 — SEQUENTIAL MANIFEST EXTRACTION, 001 → 044

**Order matters: ascending stage number**, so every citation resolves against an already-extracted
producer.

**The per-stage directive is ONE REVIEWED TEMPLATE.** Gauntlet it ONCE, hard (draft → Codex
design-review → fold → Codex confirm), then execute it 44 times with a stage parameter. Running the
full gauntlet per stage would be ~176 review passes saying the same thing, and reviews that
repetitive degrade into rubber-stamps — worse than not running them.

**Per stage:** execute (Codex) → composite checker → **FRESH fidelity agent** → append the progress
ledger row → commit. Nothing rolls forward on a negative verdict.

**Template must carry (all now defaults):** export-completeness; **loci are REAL POINTERS** in
`evidence` as in `dim_source` (58 of 92 were prose restatements naming functions that did not
exist); an **expected-claim inventory enumerated BEFORE writing** (else omission passes silently —
an executor could export one claim, attach every tooth to it, and stay green); a **negative control**
proving each exercised checker path actually ran (a passing run emits no marker); mandatory `.out`
citation; real teeth from `TOOTH_ORDER`/`ABLATION_DESCRIPTIONS`; and **honest gap recording** —
never invent a plausible locus.

⭐ **FINDINGS ARE THE PRODUCT. Green is not the goal.** A faithful manifest that turns the build red
is a SUCCESS; a green build bought by contortion is a failed extraction. Catalogue findings as they
appear; fix them in a dedicated later phase (§6).

---
## 5. AFTER PASS 2 — RE-EVALUATE THE C4 SCOPE

`DIMENSIONAL_CONSISTENCY` was deliberately scoped to manifest-internal algebra **because**
identity→dimension resolved **0/92** symbols — there was no per-quantity binding to look up. If the
shared module is keyed by quantity identity, that blocker is gone and stronger certification may be
recoverable. **Re-evaluate then — assume nothing either way.** Do not reopen the scope decision
before the module exists; four rounds were spent there (rationale in commit `2030e344`).

---
## 6. FINDINGS LEDGER + REMEDIATION PHASE

Every finding gets catalogued as encountered, never silently fixed to make a build green.
**Already open:** `RANGE_ENDPOINT_DRIFT` — census `(39,49,10)` vs stage043's `count_range (40,49,9)`,
because stage032 declares knob `Q_E` as `{low:0,high:1}` while stage043 counts it at both endpoints
as route-ful debt. stage043 is probably right; NOT adjudicated. The schema has no
`reclassified`/`reconciled` lifecycle action, so supersession cannot yet be expressed honestly.

**Progress tracking:** Pass-1 progress lives in `manifests/DIMENSION_REWRITE.md`
§3. The planned Pass-2 ledger has not been created because Pass 2 has not
begun; create it when extraction starts, append one row per stage, and read it
as a `tail`. Each row records stage id, Pass-1 dimension status, Pass-2
manifest status, checker verdict and headline, fidelity verdict and finding
counts, commit sha, and open findings. **A row that cannot show its fidelity
leg means the stage is not done.**

After all 44: the ledger-wide integration report → catalogue every finding → a dedicated remediation
phase that adjudicates and fixes them (including `Q_E` and any schema gaps such as the missing
supersession action).

---
## 7. WHERE PARALLELISM IS STILL SAFE

Read-only work only, and it is where the Workflow tool still earns its keep:
- evidence gathering for the single stage currently in the canonical
  `DIMENSION_REWRITE.md` loop;
- post-hoc fidelity audits of already-completed stages;
- independent re-verification legs.
**Never** parallelize anything that writes shared state — the module, the checker, the config, or
the progress ledger.

**Grok is a GATE before any freeze**, not a per-stage reviewer. It found the fanout blockers in one
pass after hours of self-review had not. Run it before re-freezing in §3.

---
## 8. SEQUENCE SUMMARY

1. Use the already-measured inventory in `DIMENSION_REWRITE.md` §7; do not
   re-run the superseded survey.
2. Follow `DIMENSION_REWRITE.md` §4 one stage at a time; do not revive the
   superseded extremes-first module-design pass.
3. Complete the remaining per-stage rewrites, applying that document's three
   distinct gates to every stage.
4. Re-run all `.wl` (≤2 seats), regenerate `.out`.
5. Re-pin digests in the 4 existing manifests; settle script/module/`.out` digests.
6. Checker round: delete the reverse-engineering machinery, canonical stage inventory, enforce
   export-completeness, tautology advisory, priority fixtures → **Grok gate** → re-freeze, record sha.
7. Gauntlet the per-stage template ONCE.
8. Extract 001→044 sequentially: execute → checker → fidelity → ledger row → commit.
9. Ledger-wide integration report → findings catalogue.
10. Re-evaluate the C4 scope; then the remediation phase.
