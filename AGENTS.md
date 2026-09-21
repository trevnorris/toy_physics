# Local execution safeguards

The user explicitly resumed S11c-d on September 21, 2026 after the desktop
freeze investigation and resource-guard tests. The freeze cause remains
unconfirmed. Preserve the incident record; do not replay completed jobs to
reproduce the freeze. All resumed scientific work must use the guard below.

For future S11c Python constructors, validators and export jobs, use
`scripts/s11c_guarded_run.py` around the existing supervisor. It requires a
host systemd user manager, verifies a whole-job 2 GiB memory cap, disables
job swap, limits process count and CPU affinity, lowers CPU/I/O priority,
records resource samples and fails closed. Keep the native worker's existing
limits too. Never fall back to an unguarded launch if containment fails.
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
