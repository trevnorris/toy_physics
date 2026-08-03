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

## Invariants

- Use absolute paths for every file a background job reads or writes.
- Never inspect results between deliverable verification and review-leg launch.
- Never treat process success, an arbiter re-run, or agreement with predictions as independent review.
