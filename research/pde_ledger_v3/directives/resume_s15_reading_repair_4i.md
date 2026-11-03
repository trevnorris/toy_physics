Note from the orchestrator. I stopped `repair4_scalar_wl_c1_momentum_03`, and I have read your report of that stop (`method_trial_external_stop_03.md`). The user has made a further decision about the trial.

## The stop

On 2026-11-03 I stopped `repair4_scalar_wl_c1_momentum_03` through its unit (`systemctl --user stop s11c-guard-a21d67f54053`). The journal records the stop request at 21:05:48Z and the unit stopped at 21:05:54Z. The run had not run out of memory. Its attempt on `region/component_1/momentum_labels_{drain1, drain2}` started at 19:22:48Z. I changed none of your files, in the repository or in your work directory.

## The user's decision

The trial measures how fast your construction is. The user decided that in this trial, an attempt is stopped 1 hour after it started and is reported as not returned within 1 hour. The Wolfram trial of `region/component_1/momentum_labels_{drain1, drain2}` is finished: its three attempts are its result, so do not attempt it again in Wolfram in this trial.

A trial stop at 1 hour is not an item-4 label. The final runs keep the rule in `resume_s15_reading_repair_4g.md`: an exact item-4 attempt is stopped and labelled only when it runs out of memory, or when it is still running 6 hours after it started.

## What must be true

Everything under "What must be true" in `resume_s15_reading_repair_4.md`, `resume_s15_reading_repair_4d.md`, `resume_s15_reading_repair_4f.md`, `resume_s15_reading_repair_4g.md` and `resume_s15_reading_repair_4h.md` stands. Also:
- No trial outcome is printed as an item-4 label. A trial attempt ends when it returns its conditions, when it runs out of memory, or 1 hour after it started.

## What to do now, and where to stop

- Try the construction on the other eleven trial objects:
  - in Wolfram first, `region/component_1/rate_labels_{drain1, drain2}`, `region/component_2/momentum_labels_{drain1, drain2}`, `region/component_2/rate_labels_{drain1, drain2}`, `region/component_3/momentum_labels_{drain1, drain2}` and `region/component_3/rate_labels_{drain1, drain2}`;
  - then all six in SymPy, including `region/component_1/momentum_labels_{drain1, drain2}`.
- Run through the guard, one job at a time, with `--memory-gib 55`, and poll each job until it ends. Each attempt stops 1 hour after it started. If an attempt ends on an error in your code, you may fix it and run that object again; report every attempt.
- Then stop and write the short report that `resume_s15_reading_repair_4h.md` asks for. Cover all twelve objects, with every attempt. For each attempt give the run, its start and end, whether it returned its conditions, ran out of memory or reached 1 hour, its wall time and its peak memory. Also say what you changed for each of the two points in that note, and list the files you changed with their sha256. Do not start the final runs; the orchestrator will reply.

Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
