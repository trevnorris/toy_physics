# S15 item-4 Singular engine: resume note 1

Note from the orchestrator to the builder of `research/pde_ledger_v3/directives/S15_item4_singular_engine_directive.md`,
accepted at `d4cd6d09`. This note resumes your session. The directive stands, except where this note says otherwise.
The lookups behind this note are in
`research/pde_ledger_v3/directives/_measurements/resume_s15_singular_engine_1_lookups.md`.

## What happened

- Your session ended at about 03:45Z on 2026-11-09, while your guarded run `attempt01` was running.
  - The tool call that ran `attempt01` ended with your session, and the guard's launcher ended with it.
  - Your completion watcher waited on that launcher, so it fired early. The message it queued was never delivered.
- `attempt01` kept running in its unit, alone, at `--memory-gib 56`. It started at 03:42:21Z.
  - Its last output came at 04:14:41Z, on the first of the 144 objects in scope: `WORKING_ON: differential ideal
    round 3`. It printed nothing after that.
  - Its peak memory was 11.2 GiB, with no out-of-memory event.
- At 18:32:56Z the orchestrator stopped `attempt01` on the user's decision, below, after 14 hours 50 minutes.
  - Its guard record holds resource samples up to then.
  - Neither the child's outcome file nor the launcher's outcome file exists.
- The orchestrator reconciled the stale reservation by the guard's rule, keeping the old registry. The guard admits
  new jobs again.
- The orchestrator has kept a copy of your current source and README change. You may change them.

## User decision (2026-11-09)

The user was asked what should happen. They chose "Stop and rethink (Recommended)", whose text read:

> I stop the run and pause the resume script so the builder doesn't just restart it. I write a short note recording
> your decision and what was measured. The builder then picks a method that can finish, tries it on a few objects
> first, and runs objects in parallel across the cores. The note gets its two review legs before the builder resumes.

They then added: "let's let codex work through rethinking these."

## What must be true from here

- **The construction.** The construction is yours, as the directive says, and it may change. Whatever construction you
  use meets every line of the directive's section "What must be true".
- **A trial comes first.** Before any run over all objects in scope, run your construction on a few objects you
  choose. Between them, the trial objects include:
  - the first object in scope;
  - objects at a surface and over the region;
  - objects of each group kind in scope;
  - at least one object that carries both drains.

  Each trial attempt runs through the guard. For each, the report gives its wall time, its peak memory and how it
  ended.
- **A trial attempt is a measurement.** It may end in any of the directive's ways. You may also stop it when you judge
  that its method will not finish.
  - A stopped trial attempt labels nothing.
  - The run record reports a stopped trial attempt as stopped, with what it printed.
- **After the trial.**
  - If every trial object returns, run all objects in scope.
  - If one does not, you may change the construction and run the trial again.
  - If you cannot make every trial object return, do not start the full run. Write your report, saying what you tried
    and what you measured, and stop.
- **Run objects in parallel across the machine's cores.**
  - Each attempt runs through the guard, with `--pool s15-singular`.
  - Where attempts overlap, size each cap from the peaks the trial measured, so the jobs share the machine.
  - The directive's rule for out-of-memory endings stands. An attempt that runs out of memory at a shared cap is run
    again, alone, at the largest cap the guard admits.
- **Poll each job until it ends.**
  - For every guarded job you start, the guard's launcher stays alive until the job ends.
  - Your session waits for each job with a blocking poll.
  - Do not end your turn while a job of yours is running. Do not rely on a watcher or a queued message.
- **The records.** The run record and the manifest list `attempt01` and every trial attempt. For each they give its
  start, its end, its cap, its outcome, its wall time and its peak memory, as the guard records show them. For
  `attempt01` they also say what its record lacks, and why. They present no missing value as measured.
- **Stopping.** When the final run has ended, or when you stop without a full run as above, write the directive's
  final report, and stop. Commit nothing, and start no other agent or review.

Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
