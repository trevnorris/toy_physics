# Development Pipeline — how we do the work (canonical "how we work" front door)

**Why this doc exists.** The *what* of the program lives in `docs/development_plan.md` (the sector-by-sector work breakdown) and
`docs/conceptual_foundation.md` (the physical vision). This doc is the ***how*** — the process every compute gate runs through. It was
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
Read `docs/conceptual_foundation.md` + the relevant memories first. Don't re-derive what's already banked (e.g. the PN ladder, prior no-gos).

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

## Cross-cutting standing rules
- **YAML/markdown for any file an LLM reads/writes; JSON only machine-to-machine** [`no-json-for-llm-io`].
- **No fake `python3 -c` commentary scripts** — read and reason; real ablation re-runs are fine, narration-scripts are not [`no-fake-scripts`].
- **Render audit/verify prompts under the project root** (`software/stage1_solver/_scratch/`), never `/tmp` (agents/Codex can't read `/tmp`) [`audit-prompts-under-project`].
- **Offload to agents** (reviews, broad reads, distilled summaries) to preserve orchestrator context; read grep/tail summaries, not raw logs [`offload-to-agents`, `offload-review-gauntlet`].
- **Never `pkill`/`killall` by pattern** (shared box; the user runs other Codex sessions) — kill only captured PIDs or via `TaskStop` [`background-process-launch`].
- **Path discipline:** project root is `/var/projects/toy_physics` (NOT `toy_projects` — a transcription attractor) [`toy_physics-path-typo`].

---

## Per-gate quick checklist
1. ☐ Directive drafted (requirement + acceptance + able-to-fail verdict ladder; no pre-designed route)
2. ☐ Codex design-review (xhigh, agent) → GREEN → GLM pass → fold → Codex confirm → GREEN
3. ☐ **User gate**
4. ☐ Execution prompt drafted → Codex design-review (agent) → fold
5. ☐ Codex executes dual-engine (danger-full-access, xhigh, `< /dev/null`, scripts `timeout 600`) → exit 0; comprehensive inline dim-firewall, able-to-fail
6. ☐ Tri-review: arbiter re-run (orchestrator) + fidelity (fresh agent) + adversarial-with-ablation (fresh agent)
7. ☐ Any hole → Codex remediates → **re-verify on a FRESH clean agent**
8. ☐ Bank: `STATUS.md` + `decisions/13` §0 + memory; commit only if asked
9. ☐ **User gate** before the next gate

---

*Source of truth = the `feedback-*` memories cited above. Companion docs: `docs/development_plan.md` (what to build), `docs/conceptual_foundation.md` (the vision), `STATUS.md` (where we are).*
