Note from the orchestrator. The user has decided that SymPy attempts no item-4 object in this build. I stopped your SymPy final run and ended your session.

## The stops

On 2026-11-06, on the user's decision:
- `repair4_decision_final_sympy_02` (unit `s11c-guard-06d785790b54`). I stopped the unit (`systemctl --user stop s11c-guard-06d785790b54.service`) at 17:42:28Z, with its launcher alive. The guard wrote its `outcome.json`, and its reservation left the guard registry.
- At 17:42:36Z I ended your session. After that, no process of `run_parallel.py`, `run_sympy02.py` or `wait_event.py` was running. The guard launcher of `repair4_decision_final_wolfram_01` and its job kept running, and its reservation is in the guard registry.

I changed none of your files, in the repository or in your work directory.

## The user's decision

The user's words: "Feels like we need to stop using sympy for this work. It was great for symbolic work, but without being compiled it's just too slow". The user decided:
- SymPy attempts no item-4 object in this build.
- SymPy's final run still prints every object that items 1, 2, 3, 5 and 6 ask for.
- Wolfram's final run goes on as it is.

## What must be true

Everything under "What must be true" in `resume_s15_reading_repair_4n.md` stands, with these changes for SymPy:
- In SymPy's final run, no item-4 object is attempted. Each item-4 object is printed with a not-obtained label that says it was set by this decision and that the object was not attempted in SymPy. This is a fourth kind of label. In SymPy it replaces every other label an item-4 object could carry, the third kind of `resume_s15_reading_repair_4n.md` included. No Wolfram object carries it.
- The run record lists `repair4_decision_final_sympy_01` and `repair4_decision_final_sympy_02` as history, each with how it ended as its own files record it. Neither gives any object a label.
- `repair4_decision_final_wolfram_01` is Wolfram's final run. Do not stop, signal or change it or its job.

## What to do now, and where to stop

- Make SymPy's final run from its final code, guarded, beside Wolfram's, with a `--memory-gib` near its own expected peak and within what the guard admits. Poll it until it ends.
- Poll Wolfram's final run until it ends, and make every attempt again that the rule in `resume_s15_reading_repair_4n.md` requires.
- Then finish as `resume_s15_reading_repair_4n.md` says: the updated run record and a short final report. Stop there; the orchestrator will reply.

Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
