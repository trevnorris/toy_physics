Note from the orchestrator. I stopped `repair3_90441_wolfram_final_19`, and I have read your report of that stop (`repair4_external_stop_final19.md`). The user has made a further decision about item 4.

## The stop

On 2026-11-03 I stopped `repair3_90441_wolfram_final_19` through its unit (`systemctl --user stop s11c-guard-4f8825ff30c7`), with its launcher alive; the journal records the stop at 04:40:21Z. The run was printing, and it had not run out of memory. Its attempt on `region/component_3/momentum_labels_{drain1, drain2}`, the fifth of the 48 item-4 objects in its order, started at 2026-11-03T03:56:20Z. Before that, the run had given `region/component_2/rate_labels_{drain1, drain2}` the stop kind `silent_stop`, after an attempt of 975.9 seconds, with a recorded run peak of 3583213568 bytes. I changed none of your files, in the repository or in your work directory.

## The user's decision

`resume_s15_reading_repair_4f.md` said that an item-4 attempt that goes silent for a long stretch is labelled, and it did not say how long. That sentence was mine; the user's decision was the 6-hour stop. The user has now decided: an exact item-4 attempt is stopped and labelled only when it runs out of memory, or when it is still running 6 hours after it started. Printing nothing for any stretch shorter than that is not, by itself, a reason to stop or label it. This holds in both engines.

So `region/component_2/rate_labels_{drain1, drain2}` carries no label from `repair3_90441_wolfram_final_19`; attempt it again. The two labels for running out of memory under `resume_s15_reading_repair_4d.md`, and the label for `region/component_2/momentum_labels_{drain1, drain2}` from `repair3_90441_wolfram_final_18` under `resume_s15_reading_repair_4f.md`, stand.

## What must be true

Everything under "What must be true" in `resume_s15_reading_repair_4.md`, `resume_s15_reading_repair_4d.md` and `resume_s15_reading_repair_4f.md` stands. This decision replaces the bullet in `resume_s15_reading_repair_4f.md` about an attempt that goes silent. Also:
- no item-4 object carries a label for going silent;
- every item-4 object that carries a label carries one of the two kinds, the 6-hour stop or running out of memory, and the run record says which.

## Going on

- For an item-4 attempt, this replaces the rule about silent runs in `brief_s15_reading_repair_3.md` line 221.
- The cap stays `--memory-gib 55`. Run the jobs one at a time.
- Finish as the brief says: the final guarded run of each engine from the final code, the updated run record, and a short final report.

Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
