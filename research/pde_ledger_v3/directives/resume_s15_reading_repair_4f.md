Note from the orchestrator. I stopped `repair3_90441_wolfram_final_18`, and I have read your report of that stop (`repair4_external_stop_final18.md`). The user has made a further decision about item 4.

## The stop

At 2026-11-03T03:26:17Z I stopped `repair3_90441_wolfram_final_18` through its unit (`systemctl --user stop s11c-guard-b084eb25cac6`), with its launcher alive. The run was printing, and it had not run out of memory. Its attempt on `region/component_2/momentum_labels_{drain1, drain2}`, the third of the 48 item-4 objects in its order, started at 2026-11-02T17:34:20Z; when I recorded the run's state at 03:26:08Z, its kernel had run 9 hours 51 minutes. Your report gives the run's recorded peak, 42911830016 bytes. Its last progress line is the last progress line in that run's stdout. Your runs `repair3_90441_wolfram_final_16`, `_17` and `_18` each reached this object, and none returned its conditions. I changed none of your files, in the repository or in your work directory.

## The user's decision

The user chose to label and finish. An exact item-4 attempt on an object that is still running 6 hours after it started counts as an attempt that failed, as one that runs out of memory does under `resume_s15_reading_repair_4d.md`. That object then carries a not-obtained label, and the repair goes on to finish; this repair does not look for another construction for that object. This holds in both engines.

The attempt on `region/component_2/momentum_labels_{drain1, drain2}` in `repair3_90441_wolfram_final_18` is such an attempt. Label that object from it, and do not attempt it again.

## What must be true

Everything under "What must be true" in `resume_s15_reading_repair_4.md` and in `resume_s15_reading_repair_4d.md` stands. For an object whose attempt was stopped at 6 hours, these are also true:
- the attempt was of the kind `resume_s15_reading_repair_4.md` requires: one that can return conditions on the exponents;
- its label says what was attempted, in which run, that the attempt was stopped after 6 hours or more without returning, the last step it reported, how long it ran, and the peak memory of the run;
- every condition printed now that is true stays printed for it, as before;
- the final run of each engine ends with exit 0 and prints every object, each item-4 object either with its conditions or with its label;
- the run record lists, for each engine, every object that carries such a label or the label for running out of memory, and which of the two it is.

## Going on

- A stop at 6 hours is not a reason to stop and report; it is the case above.
- For an item-4 object, an attempt that goes silent for a long stretch (`brief_s15_reading_repair_3.md` line 221) is handled the same way: label it, with the last step it reported, instead of changing the construction.
- The cap stays `--memory-gib 55`. Run the jobs one at a time.
- Finish as the brief says: the final guarded run of each engine from the final code, the updated run record, and a short final report.

Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
