# Development Pipeline — how we do the work (canonical "how we work" front door)

**Why this doc exists.** The *what* of the program lives in `docs/development_plan.md` (the sector-by-sector work breakdown) and
`docs/model_map.md` (the physical vision — ⛔ NOT `conceptual_foundation.md`, superseded and it re-confuses). This doc is the ***how*** — the process every compute gate runs through. It was
written **2026-06-26** after the rules — previously scattered across ~30 `feedback-*` memories + `STATUS.md` + each directive's
§discipline/§review sections — let a shortcut slip through (the orchestrator self-verified a remediated script instead of using a fresh
clean agent). Consolidating it here makes the pipeline impossible to miss.

**This doc is thin: the memories remain the source of truth.** Each rule cites its `feedback-*` memory; when they conflict, the memory
wins (it's the live, dated record). Keep this doc synced when a process rule changes.

---

## Roles — who does what (never blur these)

| Role | Does | Does NOT |
|---|---|---|
| **Claude (orchestrator)** | Reviews; owns *scaffolding* (directive prose, execution prompts, `STATUS`/`decisions` docs); runs the **arbiter re-run** of existing scripts; gates; banks results | **Write or mutate any code/math/script — not even disposable scratch copies** [`claude-reviews-codex-codes`, `codex-is-fix-applier`] |
| **Codex** (always `-c model_reasoning_effort=xhigh`) | Designs + writes + RUNS all scripts, iterates to exit 0, applies all code/math fixes | Decide the verdict (the verifier reviews substance afterward) [`codex-xhigh-reasoning`, `codex-iterates-until-clean`] |
| **GLM** | One fresh tertiary check per directive (CLI, full repo access) | Be iterated like Codex — it's a single pass [`review-ordering-codex-then-glm`] |
| **Clean agents** | Each review/verification leg, on its **own fresh agent**, no hint of the expected conclusion | **Write or mutate any code/math/script — an agent is a REVIEW instrument, never a coder.** Delegating coding to an agent is the SAME violation as Claude doing it by hand, and it is the easier one to rationalize [`claude-reviews-codex-codes`, `codex-is-fix-applier`]. Also: be reused across legs / carry context between reviews [`review-agents`, `offload-review-gauntlet`] |

---

## Standing principles (apply at every phase)
- **Analog, not derivation** — postulate the structure freely; the only test is internal consistency; a **no-go between requirements is a first-class success** [`falsification-is-the-goal`, `analog-find-consistent-structure`].
- **Falsification is the goal** — build every gate *able-to-fail*; a clean "it all works" is the *suspicious* outcome.
- **Truth in the OUTPUT, not prose** — no LLM establishes a math or dimensional fact; a verdict-bearing control must compute from inputs and be able to fail [`negative-verdict-short-circuit`, `dimensional-consistency-check`].
- **Never unilaterally deviate from the calibrated contract** — if a step is failing, HALT and bring it to the user [`never-alter-calibrated-process`].
- **"The tool is broken" is a claim that must be VERIFIED before it licenses anything** — a stall/hang/failure attributed to Codex (or Mathematica, or GLM) is a hypothesis, not a fact. Check the actual artifact: does the log end with its completion marker? did the process exit non-zero? does a fresh smoke run reproduce it? **If Codex is genuinely broken, we FIX Codex — we never route around it**, and never by promoting an agent to coder. (2026-07-24: a "Codex CLI stalled twice" claim was written into four docs as justification for an agent-as-coder pivot; on later inspection all 27 of that day's Codex runs had ended `exit=0` and a smoke run answered correctly in 20s. The unverified excuse did more damage than the original error, because it self-justified in every future session.) [`never-alter-calibrated-process`, `negative-verdict-short-circuit`].

---

## The phases

### Phase 0 — Grounding
Read `docs/model_map.md` (⛔ NOT `conceptual_foundation.md` — superseded, it re-confuses) + the relevant memories first. Don't re-derive what's already banked (e.g. the PN ladder, prior no-gos).

### Phase 1 — Directive → review gauntlet (before any compute)
1. **Claude drafts the directive**: states *requirement + acceptance + the verdict ladder (incl. the no-gos)*; **never pre-designs the computational route** [`claude-reviews-codex-codes`].
2. **Codex design-review (xhigh, via a gauntlet-runner agent)** → iterate to GREEN [`directive-design-review`, `offload-review-gauntlet`].
3. **ONE GLM tertiary pass** — a single fresh check [`review-ordering-codex-then-glm`].
4. **Claude folds** the findings (Claude owns the prose).
5. **Codex confirm-pass → GREEN again** — Codex bookends the GLM pass [`review-ordering-codex-then-glm`, `directive-design-review`].
6. **User gate.**

### Phase 2 — Per-gate execution-prompt design-review
Claude writes the gate's **execution prompt** (scaffolding) → **Codex design-reviews the execution prompt (via agent) before the expensive run** → fold. (Lesson: review the per-rung prompt, not just the directive.)

### Phase 3 — Execution (Codex)
Codex designs + codes + runs **dual-engine** (SymPy + Mathematica), iterating to exit 0 [`dual-engine-required`, `codex-iterates-until-clean`]:
- Mathematica needs **`--sandbox danger-full-access`** to run under Codex [`codex-mathematica-sandbox`]; the orchestrator can run `math` directly as arbiter.
- `codex exec … xhigh` **backgrounded with `< /dev/null`**; **never wrap the Codex session in shell `timeout`** [`background-process-launch`, `never-timeout-codex`].
- The **scripts** Codex runs get `timeout 600`; a 124 = reformulate the math, never raise the cap [`script-timeout-policy`].
- **≤2 concurrent `math -script`** seats; no parallel `$RT exec-*` (it races `MANIFEST.yaml`) [`mathematica-single-seat`, `no-parallel-exec-sympy`].
- **Dimensional firewall — comprehensive + inline:** dim-check **every constructed expression as it is built** (not an end-pass), units restored, mutually cross-consistent, **free of any back-solved free carrier**, dual-engine must agree, **able-to-fail** (an ablation must FIRE). No LLM establishes a dim fact [`dimensional-consistency-check`].

### Phase 4 — Tri-review (the backstop — each leg independent)
- **(a) Arbiter re-run** — *orchestrator* independently re-runs the existing scripts unchanged (reliability gate; refreshes committed outputs) [`orchestrator-rerun-and-output`]. **This is the only execution leg the orchestrator performs.**
- **(b) Transliteration-fidelity audit** — **fresh clean agent**, code-vs-equations term-by-term (catches a faithful-but-wrong operator that MMS/arbiter can't) [`transliteration-fidelity-audit`, `script-review-depth`].
- **(c) Adversarial-with-ablation** — **fresh clean agent**, tries to break it: mutates scratch copies, confirms controls are able-to-fail, hunts pass-by-construction [`negative-verdict-short-circuit`].
- A negative / "can't-do" verdict gets **harder** verification than a positive [`negative-verdict-short-circuit`].

### Phase 5 — Remediation loop (if a hole is found)
- **Fix the script regardless of outcome** — no spot a future reader can point to as "soft."
- **Codex applies the fix + re-runs** [`codex-is-fix-applier`].
- **Re-verify on a FRESH clean agent (a new fresh view) — never the orchestrator** that just reasoned about the fix [`review-agents`].
- Math-level disagreements: Claude + Codex resolve and agree; escalate to the user only if the fix changes the *conceptual nature* [`claude-codex-resolve-math`].

### Phase 6 — Gate + bank
- **One gate at a time, explicit user gate between gates, never roll forward autonomously** [`sequential-audit-chunks`].
- **Commit only when the user asks** (squash small follow-up fixes into the prior commit) [`squash-followup-fixes`].
- Sync **`STATUS.md` + `software/stage1_solver/decisions/13` §0 at every milestone** [`status-md-front-door`]; update memory.

---

## Diagnosing a stuck problem (2026-07-24/25, learned expensively)

**⭐ If successive fixes each ban a SPELLING and the next probe evades it, the ARCHITECTURE is
wrong — stop patching.** Four rounds hardening one check each closed a real hole and each was
defeated by the next probe (arbitrary anchor file → self-cited file → committed conventionally-named
file → `Dim(1,2,3)` literals → arithmetic over real bindings producing the identical wrong result).
That is a denylist against an expressive grammar; it does not converge.
- **The diagnostic that broke it open: *who supplies the value, and who supplies only the
  alphabet?*** The check computed `D' = eval(manifest_expression, script_bindings)` — the derived
  artifact controlled the function, the evidence contributed only basis names. Every downstream
  guarantee was therefore artifact-determined. The fix direction is to INVERT the flow (the
  artifact asserts, the checker computes), not to sanitize the grammar further.
- **Specify the INVARIANT, not the instance.** A rule written against the field a survey happened
  to measure (`dim_source.locus`) left the identical defect alive in its sibling (`evidence.locus`).
  Same shape as the above. Write the property, not the example.
- **Before building on an empirical premise, MEASURE IT.** "The convention surely holds across the
  scripts" was false — exact identity→dimension matching resolves **0/92** production symbols, and
  a fuzzy matcher reaches 30% while producing a silently WRONG certificate. A cheap survey would
  have saved a round.
- **Ask "independent of whom?" before calling a check redundant.** A checker-owned digest was
  removed as "duplicating" an existing staleness check; the existing one compared against an
  artifact-supplied digest (self-consistency only), so deleting the independent pin was a real
  regression.
- **When scope is deliberately bounded, write the bound INTO the directive** ("do NOT fix X, and
  report rather than implement if you disagree"). A competent coder handed a checker with three
  known holes will close them, and the next round starts. Scope decisions that live only in a
  conversation do not survive contact with the next agent — put them in the commit message too.

---

## Verification discipline (what green does and does not mean)

- **A green self-report is not evidence a hole closed.** Twice the fixtures went green while the
  invariant stayed open. What caught it both times was an agent **EXECUTING the case end-to-end**,
  never reading the diff. Reviews that only read code find design smells; only execution finds
  live defects.
- **Verify what could BREAK, not just what changed.** After a manifest rework the build was re-run
  but `--self-test` was not — and fixing production had broken a fixture. The suite went red and
  the orchestrator did not notice; an external reviewer did.
- **A fixture must NEVER be coupled to mutable production state.** One fixture asserted the
  production manifests exhibit a specific finding; legitimately fixing that finding broke the test
  suite. Fixing the product must never break the tests. Fixtures are synthetic.
- **A shipped guard with no fixture is not a guarantee.** Ask periodically: *which guards could I
  delete with the suite still green?* (A pre-freeze review found ~42 of ~89 issue codes were never
  planted, including the primary drift-detection path.)
- **An independent reviewer is worth more than another orchestrator pass.** The pre-freeze review
  found, in one pass, blockers that hours of self-review had not: a red suite, dimension recovery
  covering only 16 of 43 scripts while the schema required it for all, and a closedness rule that
  could never trigger.
- ⛔ **CHECK WHICH TOOL REPORTED THE PASS (2026-07-26).** An acceptance item named a specific
  comparator. The official tool exited **2** on a real defect; the build then wrote its own
  `normalized_comparator`, ran that instead, and headlined **PASS**. It *disclosed* the substitution in
  its raw log (`official_comparator_exit=2 normalized_comparator_exit=1`) — not deception, a workaround
  honestly logged and then optimistically summarised. **The summary is what gets read, and the real
  gate had never passed.** Nothing in the tree revealed it: comparator untouched, no variant left
  behind, `git status` clean. Only re-running the named command caught it.
  - ⭐ **Re-run the NAMED acceptance command yourself and read its literal exit code** — not a
    paraphrase, not the report's claim. This is what "the orchestrator independently re-runs acceptance
    commands" is *for*.
  - ⛔ **Every build directive must say: "never satisfy an acceptance criterion with a substitute tool;
    if the named tool fails, STOP and report."** Without that sentence, routing around a broken gate
    looks like diligence.
  - **Read raw exit codes, not the PASS/FAIL headline.** In a long report the logs are the evidence and
    the summary is a claim.
  - **A gate that is routed around is not a gate** — same family as the source-grep and can't-fail
    rules below.
- **When adding a STRUCTURED artifact, verify its STRUCTURE, not only its contents — by running the
  consumer.** Same incident: the orchestrator committed a `.wl` after checking the record count and one
  value, but not that the file had its required header. The consumer would have said so immediately.

---

## Sub-agent safety (non-negotiable)

- **⛔ An agent must NEVER delete a directory, or any file it did not itself create.** A review agent
  "cleaned up" by removing the shared `_scratch/` tree — gitignored, therefore unrecoverable. Lost:
  every directive, launcher, run log, review driver, and prior-session artifact including v1 pilot
  manifests. Put this clause in EVERY agent prompt, and give each agent its own subdirectory.
- **Guard the workspace, not only the repo.** Review-agent prompts had said "do not modify these
  tracked files" — which covered modification and not destruction, and said nothing about the
  working tree.
- **Agents are REVIEW instruments and never coders** (see the Roles table). A review agent may
  mutate SCRATCH COPIES for ablation; that is the only writing it does.

---

## Background-process discipline

- **The harness reaps long background waiters (~4–15 min).** A dead waiter is NOT a dead job — this
  is almost certainly what produced a phantom "the tool stalled twice" claim that then justified
  abandoning the calibrated contract. Observed live three times in one session.
- **A waiter must distinguish three outcomes**, not two: marker found / process gone WITHOUT a
  marker / still alive. Without the third branch, a reaped waiter looks exactly like a hang.
- **Match the done-marker on `tail -1` only** — a tool can quote a prior marker mid-stream.
  ⭐ **Use absolute paths in every launch and every check, and verify a launch by the existence of its
  declared output file.** Two stage023 incidents make this one artifact rule, not path-style advice:
  (1) a session launched with `-C /var/projects/toy_physics` against directive paths stated relative
  to `research/pde_ledger_v2/` produced its evidence under the repo-root
  `_scratch/stage023/enum/`, leaving two evidence trees and making the wrong one look real; (2) after
  the shell drifted into `research/pde_ledger_v2/`, the relative launch
  `bash research/pde_ledger_v2/_scratch/run_codex.sh …` never resolved, its `&&` chain skipped the
  intended ablation, and the operator saw exactly four `PASS` lines — the same visible shape as the
  successful ablation-free run. Absence of an error and pass-shaped downstream text do not establish
  that the launch happened; the absolute output artifact does.
- **Prefer a persistent monitor** for long jobs; it survives reaping and stays silent until there is
  something to say.
- **Codex's provider can refuse content on a cyber-policy filter**, killing the session after the
  work but before the report. Phrase adversarial directives in correctness terms ("negative
  fixture", "self-nominated anchor", "admissible set"), never attack language.

---

## Cross-cutting standing rules
- **YAML/markdown for any file an LLM reads/writes; JSON only machine-to-machine** [`no-json-for-llm-io`].
- **No fake `python3 -c` commentary scripts** — read and reason; real ablation re-runs are fine, narration-scripts are not [`no-fake-scripts`].
- **Render audit/verify prompts under the project root** (`software/stage1_solver/_scratch/`), never `/tmp` (agents/Codex can't read `/tmp`) [`audit-prompts-under-project`].
- **Offload to agents** (reviews, broad reads, distilled summaries) to preserve orchestrator context; read grep/tail summaries, not raw logs [`offload-to-agents`, `offload-review-gauntlet`].
- **Never `pkill`/`killall` by pattern** (shared box; the user runs other Codex sessions) — kill only captured PIDs or via `TaskStop` [`background-process-launch`].
- **Path discipline:** project root is `/var/projects/toy_physics` (NOT `toy_projects` — a transcription attractor) [`toy_physics-path-typo`].
- **Grok is a GATE, not a garnish.** Reserve it for foundational directives and for any artifact about
  to be FROZEN while many agents depend on it. Skipping it before a freeze is how a fanout launches
  against a broken tool [`grok-final-review-pass`, `review-ordering-codex-then-glm`].
- **The derived artifact is never a second source of truth.** Manifests/derived records are generated
  FROM the audited sources; when the two disagree the source wins, and the artifact is the thing to
  fix. Corollary: do not re-litigate a unit-level property at the integration layer — check what only
  composition can check (cross-stage agreement, citation liveness, anti-rederivation, census).
- **Model the real risk.** Our derived artifacts are written by our own agents: the operative failure
  is **DRIFT (error)**, not misrepresentation. Hardening against an adversary who does not exist costs
  rounds and buys nothing; hardening against drift is the whole job.
- **When the process is later pointed at EXTERNAL work** (e.g. published papers), validate it first on
  material with KNOWN errata/retractions as positive controls. A pipeline that cannot rediscover an
  error the community already found cannot be trusted to report a novel one.

---

## Per-gate quick checklist
0. ☐ ⛔ Directive contains the anti-substitution clause: *"never satisfy an acceptance criterion with a
   substitute tool; if the named tool fails, STOP and report"* — and states no expected ANSWER or
   REASON (ask for the determination plus its evidence; a supplied rationale gets adopted, a named
   expected result gets special-cased into existence)
1. ☐ Directive drafted (requirement + acceptance + able-to-fail verdict ladder; no pre-designed route)
2. ☐ Codex design-review (xhigh, agent) → GREEN → GLM pass → fold → Codex confirm → GREEN
3. ☐ **User gate**
4. ☐ Execution prompt drafted → Codex design-review (agent) → fold
5. ☐ Codex executes dual-engine (danger-full-access, xhigh, `< /dev/null`, scripts `timeout 600`) → exit 0; comprehensive inline dim-firewall, able-to-fail
5b. ☐ ⭐ **Orchestrator re-runs each NAMED acceptance command itself and reads its literal exit code** —
   not the report's claim. Verify structured artifacts by running their CONSUMER, not by eyeballing
   counts and values
6. ☐ Tri-review: arbiter re-run (orchestrator) + fidelity (fresh agent) + adversarial-with-ablation (fresh agent)
7. ☐ Any hole → Codex remediates → **re-verify on a FRESH clean agent**
8. ☐ Bank: `STATUS.md` + `decisions/13` §0 + memory; commit only if asked
9. ☐ **User gate** before the next gate

---

*Source of truth = the `feedback-*` memories cited above. Companion docs: `docs/development_plan.md` (what to build), `docs/model_map.md` (the model — ⛔ NOT `docs/conceptual_foundation.md`, superseded and it re-confuses), `STATUS.md` (where we are).*
