Note from the orchestrator. I stopped `repair3_90441_wolfram_final_20`, and I have read your report of that stop (`repair4_external_stop_final20.md`). The user has made a further decision about item 4.

## The stop

On 2026-11-03 I stopped `repair3_90441_wolfram_final_20` through its unit (`systemctl --user stop s11c-guard-1af66ae632bc`), with its launcher alive; the journal records the stop at 17:57:00Z. The run had not run out of memory. Its attempt on `region/component_3/rate_labels_{drain1, drain2}` started at 16:51:08Z and had not reached 6 hours. Before that, the run had stopped `region/component_2/rate_labels_{drain1, drain2}` and `region/component_3/momentum_labels_{drain1, drain2}` at 6 hours. I changed none of your files, in the repository or in your work directory.

## The user's decision

None of the six region objects that carry both drains has returned its conditions in any run of this repair: two ran out of memory within minutes, three were stopped at 6 hours or more, and I stopped the sixth. The user asked whether there is a much faster exact way to compute item 4. The orchestrator obtained a read-only assessment, `research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md`. It reads the cost as coming from how item 4 is constructed rather than from what item 4 asks. It is advisory and unreviewed, and nothing in it was run.

The user chose: before any final run, look for an exact construction of item 4 that costs much less, try it on these six objects, and report how it did.

## What must be true

Everything under "What must be true" in `resume_s15_reading_repair_4.md`, `resume_s15_reading_repair_4d.md`, `resume_s15_reading_repair_4f.md` and `resume_s15_reading_repair_4g.md` stands, with these changes and additions:
- The object item 4 names, its contract and the stop rule are unchanged. An exact item-4 attempt is stopped and labelled only when it runs out of memory, or when it is still running 6 hours after it started.
- An object attempted with a changed construction carries the result of that attempt. A label from an earlier attempt stands only for an object whose construction has not changed since that attempt.
- The construction is yours, in each engine. Each engine constructs independently, and the `.wl` reads nothing from the SymPy side. The assessment is one reading of where the cost comes from, not a method to follow.
- An item-4 object either prints its complete conditions or carries one of the two labels; no label is printed for any other reason. The assessment reports that each engine can print a not-obtained label when its classification returns with branches unresolved (`S15_separation_departure.wl:1390–1393`, `S15_separation_departure.py:2143–2148`, and the same place in your SymPy working draft).
- Every printed condition is true of the whole object on its branch. The assessment reports, from reading the source only, that in the SymPy engine a substitution for the upper radius can reach the integrand without reaching the integration endpoints or a cached radial moment (`S15_separation_departure.py:444–462`, `:831–870`, `:1094–1114`, `:1235–1270`, `:1277–1287`, `:1307–1309`). It ran no case. Whatever a condition restricts is restricted alike in the integrand, the endpoints and any cached moment.

## What to do now, and where to stop

- Try the construction on the six objects `region/component_1/momentum_labels_{drain1, drain2}`, `region/component_1/rate_labels_{drain1, drain2}`, `region/component_2/momentum_labels_{drain1, drain2}`, `region/component_2/rate_labels_{drain1, drain2}`, `region/component_3/momentum_labels_{drain1, drain2}` and `region/component_3/rate_labels_{drain1, drain2}`, in each engine, Wolfram first. Run through the guard, one job at a time, with `--memory-gib 55`, and poll each job until it ends. Each attempt is under the stop rule above.
- Then stop and write a short report. For each object and engine, give: the run, its start and end, whether it returned its conditions or which of the two stops it reached, its wall time and its peak memory. Also say what you changed for each of the two points above. Do not start the final runs; the orchestrator will reply.

Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
