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
  use meets every line of the directive's section "What must be true", except where this note changes how attempts
  run.
- **A trial comes first.** Before any run over all objects in scope, run your construction on a few objects you
  choose. Between them, the trial objects include:
  - the first object in scope;
  - objects at a surface and over the region;
  - objects of each group kind in scope;
  - at least one object that carries both drains.

  Each trial attempt runs through the guard. For each, the report gives its wall time, its peak memory and how it
  ended.
- **How a trial attempt ends.** A trial attempt ends in one of the directive's four ways, under the directive's rule
  for that way.
  - You may also stop a trial attempt when what it has printed or computed shows that its method will not finish.
  - You may not stop it for elapsed time or for a period without output (`AGENTS.md`, "No time limits").
  - A stopped trial attempt is not an ending for its object. It gives that object no label. The run record reports it
    as stopped, with the evidence you stopped it on and what it printed.
- **After the trial.**
  - If every trial object ends in one of the directive's four ways, run all objects in scope.
  - Otherwise, you may change the construction and run the trial again.
  - If you stop without a full run, write a report saying that the build is unfinished, what you tried and what you
    measured. That report is not the directive's final report. It closes no object and labels none, and nothing a
    trial printed is presented as the final transcript.
- **The full run.** Run objects in parallel across the machine's cores.
  - Each attempt runs through the guard, with `--pool s15-singular`.
  - Where attempts overlap, size each cap from the peaks the trial measured, so the jobs share the machine.
  - The directive's rules for how an attempt ends, and for the final transcript, stand unchanged for the full run. An
    attempt that runs out of memory at a shared cap is run again, alone, at the largest cap the guard admits.
- **Each job's launcher outlives your turns.**
  - For every guarded job you start, the guard's launcher stays alive until the job ends, whatever happens to your
    turn or your session.
  - To stop a job, stop its unit while its launcher is still alive. The launcher then records the outcome and
    releases the reservation.
  - Never edit the guard registry.
- **Waiting.**
  - Wait for your jobs without repeated model calls (`AGENTS.md`: "Don't poll the model while waiting"). You may
    check a trial attempt you are judging as often as the judgement needs.
  - Your session stays open until every job you started has ended and you have written your report.
  - In your first session, the message your watcher queued was never delivered, because the session had already
    ended. Whatever you use to wait must not depend on a message reaching a session that has ended.
- **The records.** The run record and the manifest list `attempt01` and every trial attempt. For each they give its
  start, its end, its cap, its outcome, its wall time and its peak memory, as the guard records show them. For a
  stopped trial attempt they also give the evidence it was stopped on. For `attempt01` they also say what its record
  lacks, and why. They present no missing value as measured.
- **Stopping.**
  - When the full run has ended, write the directive's final report, and stop.
  - If you stop without a full run, write the unfinished-build report described above, and stop.
  - Commit nothing, and start no other agent or review.

Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
