Note from the orchestrator. The user has made a further decision about the scope measurement in `resume_s15_reading_repair_4j.md`, and I stopped your session to carry it out.

## The stop

On 2026-11-04 at 21:24:59Z I sent SIGSTOP to your session's process (pid 236750). I meant to hold your session still while I stopped `run_measurements.py` (pid 237622) on its own. Instead, the shell that ran your session recorded exit status 147 and ended, and your session and `run_measurements.py` ended with it, before I sent `run_measurements.py` any signal. `run_measurements.py` last wrote `latest_health.json` at 21:24:11Z, and `controller.json` is as it last saved it, at 20:50:56Z. I changed none of your files, in the repository or in your work directory.

The two SymPy attempts were not stopped. `repair4_size_py_c1_momentum_01` (unit `s11c-guard-c4214c0da91a`) and `repair4_size_py_c1_rate_01` (unit `s11c-guard-fd3dd4540a3c`) are still running. Their guard launchers (pids 237624 and 237635) are alive. Their parent is now the user service manager instead of `run_measurements.py`, and their entries in the guard registry are unchanged.

## The user's decision

Both Wolfram attempts were stopped by their caps while sharing the machine: `repair4_size_wl_c1_momentum_01` at 32 GiB and `repair4_size_wl_c1_rate_01` at 20 GiB. Under `resume_s15_reading_repair_4j.md`, each runs again alone at `--memory-gib 55`, and your sequencer would have started them only after both SymPy attempts ended. The user decided not to wait. The Wolfram attempts run again now, one at a time, at `--memory-gib 55`, alongside the two SymPy attempts.

## What must be true

Everything under "What must be true" in `resume_s15_reading_repair_4j.md` stands, including what it carries over from the earlier notes, except that a Wolfram attempt at `--memory-gib 55` runs alongside the two SymPy attempts instead of alone. Also:
- The two SymPy attempts already running are not stopped or restarted.
- At most one Wolfram attempt runs at a time.
- An attempt ends when it returns its conditions, when it runs out of memory at `--memory-gib 55`, or on an error in your code.

## What to do now, and where to stop

- Attempt `region/component_1/momentum_labels_{drain1, drain2}` and then `region/component_1/rate_labels_{drain1, drain2}` again in Wolfram, through the guard with `--memory-gib 55`, with the construction of their first attempts. If an attempt ends on an error in your code, you may fix it and run that object again; report every attempt.
- Poll every job until it ends, including the two SymPy attempts already running.
- When all four objects have ended, stop and write the short report that `resume_s15_reading_repair_4j.md` asks for, covering every attempt in this measurement. Do not start the final runs; the orchestrator will reply.

Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
