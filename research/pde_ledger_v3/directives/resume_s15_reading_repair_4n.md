Note from the orchestrator. The user has decided how item 4 finishes. I stopped both SymPy jobs that were still running.

## The stops

On 2026-11-06, on the user's decision:
- `repair4_series_unlimited_py_c1_momentum_01` (unit `s11c-guard-350cfe7325e2`). At about 05:40Z I ended your session. That also ended `launch_trial.py` and the guard launcher of this run, while its worker kept running in its unit. I then stopped the unit (`systemctl --user stop s11c-guard-350cfe7325e2.service`) at 05:41:09Z. The attempt had not returned its conditions and had not run out of memory. With its launcher gone, the guard wrote no `outcome.json`, and the run has no `_result.json`; its guard run directory keeps `stdout`, `stderr` and `resource-samples.jsonl`. Its reservation stayed in the guard registry. I removed it at 05:43:26Z and kept the old registry as `~/.local/state/s11c-guard/resource-pool/reservations.json.stale-20261106T054326Z`.
- `repair4_size_py_c1_rate_01` (unit `s11c-guard-fd3dd4540a3c`), with its launcher alive, at 05:42:51Z. The guard wrote its `outcome.json`.

I changed none of your files, in the repository or in your work directory.

## The user's decision

In this repair, none of the six region objects that carry both drains has returned its conditions, in either engine, by any construction tried. These are `region/component_1/momentum_labels_{drain1, drain2}`, `region/component_1/rate_labels_{drain1, drain2}`, `region/component_2/momentum_labels_{drain1, drain2}`, `region/component_2/rate_labels_{drain1, drain2}`, `region/component_3/momentum_labels_{drain1, drain2}` and `region/component_3/rate_labels_{drain1, drain2}`. The user decided:
- These six objects are not attempted again in this build. Each carries a not-obtained label, set by this decision, from the attempts already made on it.
- These six go on the list of runs for a machine with more memory. The orchestrator carries them there from your run record.
- Every other item-4 object gets its final runs.

## What must be true

Everything under "What must be true" in `resume_s15_reading_repair_4.md`, `resume_s15_reading_repair_4d.md`, `resume_s15_reading_repair_4f.md`, `resume_s15_reading_repair_4g.md`, `resume_s15_reading_repair_4h.md` and `resume_s15_reading_repair_4l.md` stands, and what `resume_s15_reading_repair_4l.md` supplies stays supplied. The trial rules end with this note: the 1-hour trial stop, the rule that a trial attempt runs with no time stop, the one-object limit, the rule against production source changes and the rule against starting the final runs (`resume_s15_reading_repair_4i.md`, `resume_s15_reading_repair_4l.md`, `resume_s15_reading_repair_4m.md`). Also:
- In each engine's final run, each of the six objects is printed with a not-obtained label that says it was set by this decision and that the object was not attempted in that run. This is a third kind of label. For these six objects it replaces the rules that a label is one of two kinds, the 6-hour stop or running out of memory (`resume_s15_reading_repair_4f.md`, `resume_s15_reading_repair_4g.md`, `resume_s15_reading_repair_4h.md`). For every other object those rules stand, and no other object carries this label.
- The run record lists, for each of the six objects in each engine, every item-4 attempt made on it under `resume_s15_reading_repair_4.md` and the notes after it: the run, the sha256 of the source it ran, its start and end, how it ended as its own files record it, whether another guarded job was running during it, the last step it reported and its peak memory. A trial attempt is listed as a trial. These attempts are listed as history: none of them gives any of the six a 6-hour or out-of-memory label.
- An item-4 attempt in a final run is labelled only when it runs out of memory, or when it is still running 6 hours after it started. An attempt that reaches its cap while another guarded job is running has not run out of memory; that object is attempted again, under the same rule, with the largest `--memory-gib` the guard admits when no other guarded job is running.
- The run record lists, for each engine, every object that carries a label, and the kind of each label.

## What to do now, and where to stop

- First write the report `resume_s15_reading_repair_4m.md` asks for, covering both of its attempts, from the files the runs left. Add the end of `repair4_size_py_c1_rate_01` from its guard files.
- Then make the final guarded run of each engine from the final code. Run the two at the same time, each with a `--memory-gib` near its own expected peak, the two together within what the guard admits; the user's rule for this machine is to run jobs in parallel when memory allows. Poll both until they end, and make every attempt again that the rule above requires.
- Then finish as the brief says: the updated run record and a short final report. Stop there; the orchestrator will reply.

Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
