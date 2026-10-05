# September 20 host freeze: work paused and containment added

The user reported a prolonged whole-computer freeze requiring a forced reboot.
The cause is **unconfirmed**. No scientific workload was restarted to reproduce
the incident. Accepted coordinate-source work remains at commit 3e167144.

The user subsequently confirmed that remote access also became unavailable,
preventing inspection while the laptop was frozen. Treat this as host-wide
unresponsiveness, not merely a frozen application window. A desktop workload
could still cause system-wide resource pressure, but an application-only
rendering stall is not an adequate diagnosis. The local cgroup limits and
runtime supervisor do not rely on a live remote/model connection; they cannot
guarantee recovery from a kernel or hardware lockup.

## Recorded timeline (America/Denver)

| Time | Evidence |
| --- | --- |
| 20:03:18–20:04:49 | Coordinate-source construction; exit zero, empty stderr. |
| 20:07:10–20:07:59 | Saved-source validation; 47.53 s, exit zero, empty stderr. |
| 20:13:59 | Source checkpoint written: 3,862,335 bytes. |
| 20:14:05 | Acceptance commit 3e167144. |
| 20:14:07–20:14:08 | Desktop app automatically refreshed branch diff/review statistics; logged cancellations. |
| 20:15:09 | Last available previous-boot user-journal entry; this task was planning the next binding stage. |
| 20:45:43 | Kernel text log records the new boot. |

There was no new binding or numerical construction launched after validation.
Host process inspection after reboot found no surviving S11c workload. The
machine has approximately 30 GiB RAM and 19 GiB swap; post-reboot free memory
does not establish the pre-freeze condition.

## Evidence and limits of the diagnosis

The user supplied previous-boot kernel journal and text kernel logs. The journal
returns `No entries`. The text log contains no incident-time OOM, GPU, disk or
kernel-stall report: its last pre-reboot messages are at 07:45 and it resumes
with boot messages at 20:45. A hard freeze may prevent logging. These absences
do not exclude memory pressure or a kernel/device failure.

The old Python workers used an address-space limit, but their supervisor did
not enforce a cgroup budget for the whole descendant tree. Validation imported
the scientific helpers before installing its in-process limit. No durable peak
RSS, swap, I/O or pressure measurements were recorded. The desktop application
and automatic Git diff work were outside those limits.

There are 63 tracked regular text files larger than one million bytes, totaling
358,569,770 bytes. Of these, 49 generated v3 exports and JSON manifests total
312,219,996 bytes. Automatic branch diffs coincided with the last commit, so
their extra load is a plausible contributor, **not an established cause**.
The investigation does not attribute unrelated application activity to S11c.

## Implemented safeguards

- `.gitattributes` excludes those 49 explicit generated paths from automatic
  text diffs. Their contents, scientific hashes and annex policy are unchanged.
  Handwritten code, physical input JSON and reports remain reviewable.
- `scripts/s11c_guarded_run.py` puts a job and all descendants in a separate
  systemd user service: 2 GiB memory, zero swap, 32 tasks, one CPU affinity,
  nice 15, idle I/O priority, bounded runtime and whole-group cleanup.
- The runner checks the effective cgroup limits before starting a workload,
  refuses overlapping guarded jobs, and does not fall back to an unguarded run.
  It retains the existing native limits and supervisor/checkpoint workflow.
- Resource samples are written locally every two seconds without model calls.
  Repeated host available-memory readings below 4 GiB stop the job. Source,
  stdout/stderr and completed packets remain available after a failure.
- Repository instructions and `_scratch/s11c/PAUSED_HOST_FREEZE.json` keep
  scientific launches paused pending explicit user continuation.

Small metadata-only controls verified the actual memory/swap/task/affinity/nice
limits, successful completion, preservation of child exit 7, and timeout exit
124 after about two seconds. A final smoke test passed after overlap/cleanup
changes. The pause control refused launch before creating its run directory.
No stress allocation, SymPy work, numerical solve or old validation was run.

These safeguards bound future scientific jobs and reduce automatic diff load;
they do not cap the desktop application itself or prove that an unrelated
hardware/kernel fault cannot recur. Future telemetry should make a resource
incident diagnosable without replaying accepted work.

Raw local evidence and test receipts are in
`_scratch/s11c/host-freeze-20260920/`; the concise machine record is
`S11c_d_host_freeze_checks.json`. Resume from the accepted coordinate sources
only after the user requests continuation, with the guarded runner around the
existing supervisor and silent completion/error hook.
