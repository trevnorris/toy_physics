# S15 item 4: a side test in Singular, on one object

Directive from the orchestrator. You are a new builder on this repository, and this is a side test. It is not part of the S15 reading build, and nothing it prints becomes an S15 result.

## Why

Item 4 of `research/pde_ledger_v3/directives/S15_reading_decision_list.md` asks, for each departure object of item 3 that depends on the separation, for the exact and complete conditions under which that object vanishes for every `D` on the domain, with every special exponent value treated. For the six region objects that carry both drains, no exact construction in Wolfram Language or SymPy has returned those conditions. The user wants to know whether a construction whose algebra runs in Singular returns them for one of these objects, and at what time and memory.

## The object

`region/component_1/momentum_labels_{drain1, drain2}`, as items 3 and 4 of the reading decision list define it, on the region and domain of `research/pde_ledger_v3/directives/S15_decision_list.md`. The production engines of the reading build, `research/pde_ledger_v3/scripts/S15_separation_departure.py` and `research/pde_ledger_v3/mathematica/S15_separation_departure.wl`, show how the object is formed from S15's fields. The Wolfram engine's output writes the object's name as above; the SymPy engine's output writes it `region/component_1/momentum_labels_drain1+drain2`.

## Supplied (unfalsifiable in this test)

The section "Supplied (unfalsifiable in this trial)" of `research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md`: the distance from drain 1 over the region, and its series in `z = r/d`. You may use it.

## What must be true

- The conditions are exact and complete, as item 4 asks: every special exponent value is treated, and each condition is printed as equations and inequalities in the named quantities. No limit or expansion is taken in place of an exact object. The supplied series is an identity on the whole region and may be used as one; a finite number of its terms is not the object.
- The algebra that decides the conditions runs in Singular: `/var/projects/_tools/singular/bin/Singular`, version 4.4.1, with its libraries in `/var/projects/_tools/singular/share/singular/LIB/`. Other code only prepares Singular's input and reads its output. Python 3 with SymPy is available for that.
- The construction is yours.
- The script prints computed objects and never states conclusions. Before each step it prints a line that starts with `WORKING_ON:` and names the step.
- A branch that the construction sets aside without determining its conditions is printed as not obtained.
- An attempt ends when it returns its conditions, when it runs out of memory at its cap, or on an error in your code. Nothing stops it for time or for a period without output. If an attempt ends on an error in your code, you may fix the code and run again; report every attempt.

## How to run it

- Write only in `/var/projects/toy_physics/_scratch/singular_side_test/` and in the guard log directories you name there. Change nothing else in the repository, in the guard registry or anywhere else, and run no `git` command that writes.
- Run each attempt through the guard, `scripts/s11c_guarded_run.py`, pooled as the S15 jobs are, one attempt at a time. Other guarded jobs may be running, among them the S15 final runs. Do not stop, signal or change them. Use the largest `--memory-gib` the guard admits beside them. If it admits less than 8, or refuses for lack of memory, try again every 15 minutes. Poll each job until it ends.

## The report, and where to stop

When the last attempt ends, write `/var/projects/toy_physics/_scratch/singular_side_test/report.md`, and stop. Give:
- your construction: which steps run in Singular and which run elsewhere, and the time each took;
- for each attempt: the run, its start and end, its `--memory-gib`, whether it returned its conditions, ran out of memory or ended on an error in your code, its wall time, its peak memory and its last `WORKING_ON:` lines;
- for an attempt that returned: every branch it set aside without determining its conditions, if any, and how the construction establishes that the printed conditions are complete;
- every file you wrote, with its sha256.

This task ends at your report; it includes no review.

Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
