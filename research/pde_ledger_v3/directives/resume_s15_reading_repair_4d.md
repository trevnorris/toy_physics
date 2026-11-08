Note from the orchestrator. I have read your report of the memory stop at 55 GiB (`repair4_resource_stop_55g.md`). 55 GiB is the most this machine can give one job, so the cap stays there. The user has made a decision about item 4 that changes how a run that runs out of memory is handled.

## The user's decision

For item 4, an exact attempt on an object that runs out of memory counts as an attempt that failed. That object then carries a not-obtained label, and the repair goes on to finish. The attempt must be of the kind `resume_s15_reading_repair_4.md` requires: one that can return conditions on the exponents. For that object, this replaces the rule that a construction which does not reach the conditions is a reason to change the construction. Label it and go on; this repair does not look for another construction for that object.

## What must be true

Everything under "What must be true" in `resume_s15_reading_repair_4.md` stands. For an object whose attempt ran out of memory, these are also true:
- its label says what was attempted, that the attempt ran out of memory, and where: the last step the attempt reported and the peak memory of the run;
- every condition printed now that is true stays printed for it, as before;
- the final run of each engine ends with exit 0 and prints every object, each item-4 object either with its conditions or with its label;
- the run record lists, for each engine, every object that carries such a label.

## Going on

- A run that ends because the kernel or the process ran out of memory is no longer a reason to stop and report; it is the case above. The other stopping reasons in `resume_s15_reading_repair_4.md` stand.
- The cap stays `--memory-gib 55`. Run the jobs one at a time.
- Finish as the brief says: the final guarded run of each engine from the final code, the updated run record, and a short final report.

Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
