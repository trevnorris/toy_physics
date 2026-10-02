# Local execution safeguards

The user explicitly resumed S11c-d on September 21, 2026 after the desktop
freeze investigation and resource-guard tests. The freeze cause remains
unconfirmed. Preserve the incident record; do not replay completed jobs to
reproduce the freeze. All resumed scientific work must use the guard below.

For future S11c Python constructors, validators and export jobs, use
`scripts/s11c_guarded_run.py` around the existing supervisor. It requires a
host systemd user manager, verifies a whole-job 2 GiB memory cap, disables
job swap, limits process count and CPU affinity, lowers CPU/I/O priority,
records resource samples and fails closed. Keep native memory limits too. Never fall back to an unguarded launch if containment fails.
One job at a time; no overlapping CAS or automatic retries. Lightweight
read-only inspection and editing do not require the job runner.

Run silently with durable project logs and the existing local completion/error
hook. Waiting must not invoke the model; no recurring model polling. Preserve
every accepted operand, stream and result, and resume only unfinished work.

Large generated exports and machine manifests are excluded from automatic
text diffs by `.gitattributes`. Keep those exclusions: inspect summaries and
bounded file excerpts, not an entire branch's generated payload diff. Do not
rewrite accepted scientific files merely to compact their representation.
Add explicit diff exclusions for newly generated data files exceeding 1 MiB;
keep handwritten code, physical input files and concise reports reviewable.

# Standing execution-time policy (user instruction, September 30, 2026)

All runs have no wall-clock, native alarm, CPU-time or inactivity deadline unless
the user explicitly requests a limit for a particular run. Do not add automatic
time caps to save model tokens or split a computation into timed continuations.
Keep memory containment, zero swap, host-memory protection, process/thread/CPU
controls, overlap refusal, scientific failure checks and durable checkpoints.

Use `scripts/s11c_guarded_run.py` around the supervisor. It now runs without a
time cutoff; its legacy `--seconds` argument is recorded but imposes no deadline.
Before reusing a worker, remove obsolete computational alarms/time limits in a
new version; preserve historical source snapshots. This policy supersedes older
900/840-second caps and progress-stall timeouts in project instructions. No new review or permission is
needed solely to remove these obsolete timers. Pure tooling fixes remain local
test work; scientific method/equation/claim changes retain applicable review.

# Parallel resource authorization (user instruction, October 1, 2026)

The user explicitly authorized parallel scripts and use of available machine
resources for the near-unity uniform check and upstream repair review tracks.
This supersedes the one-job restriction above for explicitly pooled jobs; it
does not authorize unguarded execution or replay of completed calculations.
Use the tested opt-in pool in `scripts/s11c_guarded_run.py`, with disjoint CPU
assignments, an aggregate memory reservation, zero job swap, and a 4 GiB host
available-memory reserve. The initial pool budget is 16 GiB across all pooled
jobs, not 16 GiB for each job. Native scientific memory limits remain required.
Legacy exclusive jobs and pooled jobs must not overlap. Preserve prior pinned
guard/supervisor snapshots and all incident records. Pooled jobs use desktop-
managed priority; no scheduler exception or host configuration change is needed.

Reviewers may work independently in parallel on fixed source snapshots. Keep
reviewer outputs separate until both finish. Review submissions within the
authorized work use the standing authorization below. At most two Wolfram
kernels may hold the two license seats;
coordinate their use with the user's other Claude task. No elapsed-time limit
applies to review kernels either. Pure tooling fixes need local tests and a
rationale, not another physics-review round. Stop after the uniform evidence
and applicable upstream review findings before starting a defect sweep.

# Standing review-submission authorization (user instruction, October 1, 2026)

The user explicitly authorized removing repeated per-packet permission requests
and running the prepared review. For the ongoing authorized research, submit
necessary source/evidence packets and substantive revisions to the established
Claude and Grok reviewers without asking for permission again. Continuing work
and applicable scientific reviews within that scope do not need routine stage
approval. Ask only when a genuine user decision, a material expansion of scope,
or a new disclosure outside these established reviewers requires user input.

Freeze and hash each submitted packet, record its contents, recipients and the
standing user authorization, and preserve the exact reports and reviewed bytes.
Keep independent reviewers separate until both finish. Do not include credentials,
unrelated private material or peer reports. Scientific method/equation/claim
changes still receive applicable assessment; pure tooling fixes use local tests
and a recorded rationale. Review verdicts do not replace execution readiness or
result validation. Guarded execution, resource controls, saved-work preservation
and the no-deadline policy remain in force.

This instruction supersedes older exact-packet consent requirements in project
notes, preparation records and completion-hook messages. Preserve those records
as history; do not treat their superseded permission language as a new blocker.
