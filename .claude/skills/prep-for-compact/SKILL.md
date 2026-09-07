---
name: prep-for-compact
description: Prepare a long session for /compact so nothing is lost across the context summary. Brings the retention docs and memories current (STATUS.md, the memory files, MEMORY.md pointers), passes the changes to Codex (codex-sol) to verify they are accurate, leak-free, and free of over-claims — iterating review-until-clear — then hands the user a self-contained resume prompt to paste back after compacting.
allowed-tools: Bash, Read, Edit, Write, Monitor
user_invocable: true
---

# Prepare for Compact

Invoke as `/prep-for-compact`. One operation: bring the retention docs + memories current, have Codex verify
them **until clear**, and hand the user a resume prompt to paste back after `/compact`. ⭐ The test the result
must pass: a future you, after the context is summarized, resumes from the docs + the prompt with **nothing
lost** — every committed SHA, every in-flight operation, every owed item, and the exact NEXT ordering.

## Runbook

### 1. Update the plans, docs, and memories to the CURRENT state
Bring every retention artifact to what is true **right now**, including work still in flight:

- **`STATUS.md`** — this repo uses a **prepend-and-supersede** top clause: write a new dated top clause
  capturing the current state, and append `[SUPERSEDED by the clause above]` to the previous top header (end
  the new clause with "The clause below is SUPERSEDED — kept for the artifact map."). Capture: what is
  **committed** (with SHAs), what is **in flight**, what is **owed**, and the exact **NEXT ordering**.
- **Memories** (`~/.claude/projects/.../memory/`) — update the project/feedback files the work touched, and the
  matching **`MEMORY.md` pointers**. ⛔ Hunt down **stale pointers** — a pointer still saying "NEXT = X" after X
  is done is the exact defect the verify will catch. ⛔ Do NOT drop a live pointer to satisfy the size hook;
  the floor is one line per live entry.
- **Attribution precision** — say *which* legs did *which* review (standard-review vs a later gate); a wrong
  attribution is a finding.

### 2. Pass the changes to Codex to verify — DETACHED, review-until-clear
The point of this step is that **the verify catches real over-claims** — it has, more than once, caught a doc
claiming something was cleared/resolved/value-free that was not. Do not skip it, and do not accept the first
FIXES-NEEDED as "good enough" — iterate.

- Render a verify prompt for **codex-sol** (`-m gpt-5.6-sol -c model_reasoning_effort=xhigh -s danger-full-access`,
  so it can read `~/.claude`). Ask it to confirm, specifically: **(a)** every cited SHA resolves *and matches*
  its description (`git show --stat`); **(b)** no **over-claim** (nothing said built/cleared/resolved that is
  only in progress); **(c)** no leaked expected **value**; **(d)** no **stale** pointer / NEXT; **(e)** the
  **resume prompt** (embed it) is accurate and leak-free. It must end with **CLEAR TO COMPACT** or **FIXES-NEEDED**.
- **Launch detached** — `run_in_background` is reaped on this box, so use a `setsid`-detached launcher that
  writes a `DONE` marker on exit, log **outside the repo**, `< /dev/null` on the `codex exec`, and a **Monitor**
  on the DONE marker. (Same pattern the build skill uses for astra.)
- **Fold review-until-clear** — a finding is not a mandate: **verify each one yourself** (G4) against the real
  files/SHAs, fold the true ones, and **re-verify** until Codex returns CLEAR TO COMPACT. Trivial factual nits
  (a line anchor, a committed-vs-uncommitted word, a HEAD SHA) still get fixed — they are what a resumed session
  trips on.

### 3. Write the resume prompt
A single self-contained block the user pastes back **after** compacting. It is **orchestrator context**:

- **Header:** one line on where we are + `Trust git log --oneline -N` + the live-handoff memory + STATUS pointer.
- ⛔⛔ **If it carries physics VALUES or knife/answer details, mark it "ORCHESTRATOR CONTEXT — ⛔ NEVER give it
  to a blind builder"** (the builder gets only its directive). This is load-bearing — the whole point of
  blindness dies if the resume prompt leaks into a build.
- **STATE** — committed SHAs (make the **HEAD reference state-independent**: "HEAD = tip of `<branch>`, trust
  `git log`", ⛔ not a pinned SHA that goes stale), what each commit is, unpushed status.
- **CURRENT WORK + NEXT** — the exact ordered next steps. ⚠ **In-flight operations get a state-independent
  instruction** — "inspect the `DONE` markers + terminal `EXIT=` lines before folding", ⛔ never "it is still
  running" (it will have finished by the time the user resumes).
- **MODEL / OPS / GOVERNING / OWED** — the model+effort map, the ops gotchas that bit this session, the
  governing ⛔s that must survive the summary, and the standing owed items (unpushed commits, MEMORY.md size,
  anything paused).
- Write it to the **scratchpad** (it is ephemeral session context, ⛔ not committed).

### 4. Hand over + CLEAR TO COMPACT
- **Commit** the doc/memory changes (`STATUS.md` + any `_measurements/` records) on the working branch —
  explicit paths, never `-A`. Memories are auto-saved outside git. ⛔ Keep the resume prompt uncommitted.
- Paste the **finalized resume prompt** to the user in a fenced code block, say **CLEAR TO COMPACT**, and list
  the standing owed items. ⛔ Then stop — the user runs `/compact` and pastes the prompt back.

## Invariants (measured)
- ⭐ **The verify earns its keep** — it has caught genuine over-claims (a "value-free" that was not; a
  "resolved" that was only in progress) that would have propagated a false belief across the compact. Iterate
  it to CLEAR; ⛔ do not hand over on a FIXES-NEEDED.
- ⚠ **Detached launch or it dies** — `run_in_background` (and plain Monitors) are reaped mid-run; `setsid` +
  DONE-marker + re-armed Monitor is the pattern that survives.
- ⛔ **State-independent flags** — a resume prompt that says "the legs are still running" or pins `HEAD` to a
  SHA is stale the moment anything advances; write what to **inspect**, not what is currently true.
- ⛔ **Never leak the resume prompt into a build** — if it names values/answers, it is for the orchestrator only.
