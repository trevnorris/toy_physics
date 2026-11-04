Note from the orchestrator. I have read your trial report (`method_trial_report.md`) and its attempt inventory (`method_trial_attempt_inventory.md`). The user has made a further decision about item 4.

## The user's decision

No trial attempt returned its conditions within 1 hour, in either engine. The user decided to measure how large this work is. Two objects, `region/component_1/momentum_labels_{drain1, drain2}` and `region/component_1/rate_labels_{drain1, drain2}`, are each attempted in both engines and allowed to run until the attempt ends on its own, however long that takes.

This is a measurement, not a final run. How an attempt ends here is not an item-4 label. The final runs keep the rule in `resume_s15_reading_repair_4g.md`.

## What must be true

Everything under "What must be true" in `resume_s15_reading_repair_4.md`, `resume_s15_reading_repair_4d.md`, `resume_s15_reading_repair_4f.md`, `resume_s15_reading_repair_4g.md`, `resume_s15_reading_repair_4h.md` and `resume_s15_reading_repair_4i.md` stands, except the 1-hour stop in `resume_s15_reading_repair_4i.md`, which this decision replaces for these four attempts. Also:
- No attempt is stopped for time. An attempt ends when it returns its conditions, when it runs out of memory at `--memory-gib 55` while running alone, or on an error in your code.
- No production source changes.

## What to do now, and where to stop

- Attempt the two objects in both engines, four attempts in all, with the construction of your last trial attempt in each engine. If an attempt ends on an error in your code, you may fix it and run that object again; report every attempt.
- Run through the guard. Run the two SymPy attempts alongside one Wolfram attempt at a time, in the order the objects are named above. While jobs share the machine, size each cap from that object's measured peaks so that the caps fit the guard pool together. An attempt stopped by its cap while sharing the machine has not run out of memory: run it again alone with `--memory-gib 55`.
- Poll each job until it ends.
- When all four objects have ended, stop and write a short report. For each attempt give the run, its start and end, whether it returned its conditions, ran out of memory or ended on an error in your code, its wall time, its peak memory, and its last progress lines. For an attempt that returned, also list every branch its construction set aside without determining its conditions, if there was any. List the files you changed with their sha256. Do not start the final runs; the orchestrator will reply.

Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
