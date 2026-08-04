---
name: build
description: Launch one Codex build and its two independent review legs as a single operation for the PDE ledger. Leak-gates the directive before launch, backgrounds an xhigh Codex run with absolute paths and no shell timeout, verifies the deliverable rather than the session, then automatically starts fresh-agent and Grok review before the caller reads the result.
allowed-tools: Bash, Read, Edit, Write, Agent
user_invocable: true
---

# Build With Automatic Review

Invoke as `/build <absolute-directive> <absolute-deliverable> --check "<physics to check>" --do-not-read <absolute-path> ...`.

This is one operation: build, verify the deliverable, and launch both review legs. Never return a
successful build without running `.claude/skills/review-legs/SKILL.md` yourself.

## Runbook

1. Resolve the directive, deliverable, log, review prompt, and review outputs to absolute paths before
   launch. Background shells retain their cwd between calls; a relative write followed by an absolute
   read has already produced a silent empty prompt.
2. Leak-gate the directive by grepping it with fixed-string `rg` searches for every pre-registered answer it must not
   contain. If anything matches, repair and re-run the gate **before launch**. Codex snapshots the prompt
   into argv, so editing the file after launch changes nothing.
3. Launch Codex at xhigh effort using the **Bash tool with `run_in_background: true`**. ⛔ Do NOT use a
   shell `&` — this harness detaches the job itself and notifies you when it exits; a `&` inside a
   foreground call leaves the build untracked and unreported.

   ```bash
   codex exec -c model_reasoning_effort=xhigh "$(</absolute/directive.md)" \
     > /absolute/path/OUTSIDE/the/repo/codex-build.log 2>&1
   ```

   ⛔⛔ **THE RAW TRANSCRIPT MUST BE WRITTEN OUTSIDE THE REPOSITORY.** ⚠ A transcript contains the
   engine's **complete tag values verbatim**, so one left in the tree is a blindness leak that no
   do-not-read list catches — it is not a `.wl`, not under `mathematica/`, not named `PREREGISTERED`, and
   not reachable by `git show`. ⚠ **Measured 2026-08-03: several such files were already sitting in
   `_scratch/`.** ⇒ `.claude/skills/review-legs/SKILL.md` § BLINDNESS IS ENFORCED BY ABSENCE.

   Add `--sandbox danger-full-access` when the build must run Mathematica. ⛔ Never wrap the command in a
   shell `timeout` — SIGKILL has cost 300k+ tokens.
4. Do not poll for completion; the harness re-invokes you when the job exits. On that notification,
   verify the **deliverable**, not exit status: require the
   requested artifact to exist and be non-empty, and require the transcript/token usage to be plausible
   for the requested build. Exit 0 plus `hook: Stop` has accompanied an empty prompt and no work; the
   measured tell was about 3k tokens instead of 37k+.
5. Do not open the deliverable or read its results. Immediately read and execute
   `.claude/skills/review-legs/SKILL.md` with the deliverable, `--check`, and `--do-not-read` arguments
   from this invocation. Launching the legs is the build skill's responsibility; do not tell the caller
   to invoke `/review-legs` later.
   ⭐ The deliverable here is **Codex-written**, so its two legs are **a fresh Claude agent + Grok**.
   ⛔ Codex does not review what Codex wrote — see that skill's authorship table.
6. Return only after both independent review legs have completed. Keep their findings separately
   attributed so the caller can filter them before acting.

## ⭐⭐⭐ EVERY SCRIPT DIRECTIVE CARRIES THESE THREE CLAUSES — non-negotiable

> **1. The script may PRINT computed objects. It may NOT state conclusions.** An `emit`/`Print` payload
> must be a CAS object — an expression, a solved root, a boolean from a symbolic test. ⛔ Never prose
> describing a result.
> **2. PRINT the residual; do NOT assert it.** `assert residual == 0` **is the builder writing down the
> expected output**, and it turns an informative value into a binary crash. Compute → emit → *then* assert.
> **3. Interpretation belongs to the STEP RECORD.** ⛔ The script does not editorialise.

⚠⚠ **Measured 2026-08-04, across three independently-built steps: only ~10–20% of emitted tags depended on
any computation.** The rest were typed sentences. ⛔ Cross-engine agreement on such a tag is **vacuous** —
both engines carry the same author's sentence — and **eight fidelity review legs did not catch it**, because
*"does it say the right thing"* and *"does it depend on anything"* are different questions.
⇒ `research/pde_ledger_v3/REBUILD_HANDOFF.md` · `[[feedback-scripts-print-never-assert]]`.

⭐ **Copy the one shape that worked** (`S11bB` lines 421–443): type a candidate, compute the object
**independently**, `emit` the symbolic difference, hard-stop if nonzero.

⭐ **After every script build, run `reduction/derived_or_declared.py` on the deliverable.** A script still
showing mostly CONSTANT has not been built. ⚠ It is a triage tool, ⛔ not a verdict — see its limits in the
rebuild handoff.

## ⭐⭐ SUPPLY WHAT IS ALREADY VERIFIED — withhold only the ACCEPTANCE CRITERION

⭐⭐ **Supply every already-verified object to the builders** — setup, field content, governing equations,
supplied premises. ⛔⛔ **DO NOT BLIND THE INPUTS.** ⚠ **Under-specification has cost this ledger far more
than contamination:** `∇·u` fell out of a spec **four times**, every time it was prose instead of an
equation; and a spec describing a system with **no linear coupling at all** would have had both engines
faithfully conclude *"it does not couple."* ⇒ `[[feedback-supply-verified-blind-only-open]]`.

⛔ **Withhold exactly one thing: an acceptance criterion that references an expected value.** ⭐ The risk is
**not** a faked computation — it is the **"fix until it matches" loop.** Codex iterates to exit 0, so if
*"matches the recorded value"* is the exit condition it will get there, and a genuine disagreement — the
most valuable output available — becomes silent confirmation.
⇒ **The builder's job ENDS at compute-and-print. The diff happens on our side**, where a mismatch is a
**finding**, ⛔ not a build failure.

⚠ **State IN the directive when an object is supplied** and therefore **unfalsifiable within the build**, so
a passing build does not read as if it verified it.

⚠ **Blindness is now a MINOR control, ⛔ not the architecture.** Clause 1 removes the slot a typed answer
goes in, which is structural; blindness is behavioural and has failed here repeatedly, including through
`git show`. ⭐ Keep only: **move** answer-bearing files out of the tree, and **audit the build log**
afterwards. ⛔ Stop there — enumerating git internals and backup archives is anti-adversary, and the engines
are not adversaries.

⛔ **Verify a builder's compliance declaration against the DIFF, never against its report.**

## Invariants

- Use absolute paths for every file a background job reads or writes.
- Never inspect results between deliverable verification and review-leg launch.
- Never treat process success, an arbiter re-run, or agreement with predictions as independent review.
