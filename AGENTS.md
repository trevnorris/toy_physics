# Instructions for Codex in this repository

Written by Claude (orchestrator) for the 2026-10 cleanup. These replace the earlier S11c-era instructions; those
are kept in the tag `archive/pre-cleanup-2026-10-04`. `CLAUDE.md` governs the research method. This file governs
how you work in the repository.

## Scope and stopping

- Do the task the user gives you, under the directive it names, and nothing else. Before starting, write down the
  deliverable in one sentence. When it is done, or a stop condition in the directive is hit, stop and report.
- If the work turns into prerequisites for prerequisites, stop and report:
  - a second method failure;
  - a new sub-problem the directive did not name;
  - a repair to a repair.

  Say what blocks the deliverable and what the options are. ⛔ Do not open a new line of work without the user's
  go-ahead. Four weeks of S11c work produced 1,200 commits and no answer this way.
- Questions that the plan assigns to a later step (`research/pde_ledger_v3/V3_STEP_PLAN.md`) go to that step. Do
  not answer them inside the current one.

## Git hygiene

- **Commit durable artifacts only:** specs, directives, scripts that produce reported results, result reports,
  review reports, step records, and outputs that a record relies on. Use one commit per meaningful unit of work,
  with a message that says what changed.
- ⛔ **Do not commit process records.** That means launch, admission, readiness, receipt, hook, guard, resume,
  journal or checkpoint files, resource samples, "record that X started" notes, and per-attempt review prompts.
  Run state lives under `_scratch/`, which is ignored. Move a run output into the tree only when a record cites it.
- ⛔ Never `git add -f`. A file that `.gitignore` excludes stays untracked.
- ⛔ Never edit `.gitattributes` or `.gitignore` unless the user asks. The annex policy there is the user's.
- Use `datalad save`. The large `out/*.out` files are annexed automatically. ⛔ Never annex a `*_exports.py`.
- ⛔ Never push, force-push, or move or delete tags. The user or Claude pushes after review.

## Reviews

- By default, Claude runs the reviews at the directive's STOP points. Do not launch, wait for or record review
  rounds yourself unless the user's current instruction asks you to.
- A review verdict is quoted literally, with its scope. "Scoped clear" is not "cleared", and "unresolved" is not a
  result.

## Running computations

- Run CAS and numerical jobs through `scripts/s11c_guarded_run.py`: per-job memory cap, zero swap.
  After the 2026-09-20 host freeze (cause unconfirmed), ⛔ never fall back to an unguarded launch.
- **No time limits** (user, 2026-09-30). Add no wall-clock, CPU-time or inactivity deadline unless the user asks
  for one for a specific run.
- The machine is dedicated to this work; use all of it (user, 2026-10-15). There is no host-memory reserve. Pooled
  jobs (`--pool`) share one machine-wide registry across every checkout. Each takes one logical CPU and keeps its own
  hard cap (`--memory-gib`, 1 GiB up to `MemTotal` rounded down to whole GiB). A job is admitted when the registry is
  clean, a CPU is free, the reserved caps plus its own fit within that budget, and available memory covers its cap
  plus the unused caps of jobs already admitted. Otherwise it is refused. A legacy (non-pooled) launch holds the
  machine-wide admission lock for its whole run: a pooled launch from another checkout waits until it ends, and one
  from the same checkout is refused.
- At most two Wolfram kernels at once (two licence seats). No automatic retries, and no replaying completed jobs.
- Don't poll the model while waiting.

## Reporting

At each stop, write a short report (at most one page) covering:
- what was done;
- the commit IDs;
- the files to read;
- open questions.

Quote numbers exactly and keep their qualifications. No process narrative.
