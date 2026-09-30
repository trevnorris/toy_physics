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
