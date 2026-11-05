Note from the orchestrator. I have read your report of the guard refusal (`measurement_guard_refusal_02.md`). The user has made a further decision about item 4.

## The refusal

The guard refused `repair4_size_wl_c1_momentum_02` at admission: `--memory-gib 55` beside the two running SymPy attempts is more than the host's available memory covers. That was my error; `resume_s15_reading_repair_4k.md` asked for a cap the guard could not admit.

## The user's decision

The user asked whether item 4 could be computed at much less cost, and chose this: supply the expansion below, try a construction of item 4 that uses it on one object in each engine, and report how it did, before any further measurement or final run.

- The Wolfram attempts that `resume_s15_reading_repair_4j.md` and `resume_s15_reading_repair_4k.md` ask for are not made.
- The SymPy measurement attempt `repair4_size_py_c1_rate_01` (unit `s11c-guard-fd3dd4540a3c`) keeps running under its guard launcher. Do not stop or restart it, and do not wait for it; the orchestrator watches it.
- The SymPy measurement attempt `repair4_size_py_c1_momentum_01` (unit `s11c-guard-c4214c0da91a`) ended on 2026-11-05 at 06:50:33Z, when the kernel's out-of-memory killer stopped it at its 2 GiB cap while it shared the machine. Under `resume_s15_reading_repair_4j.md` that is not running out of memory. It is not run again now, and its end is not a guard event for you to stop on.

## Supplied (unfalsifiable in this trial)

With the region's surfaces taken as spheres about `X_2` of radii `r_in` and `R` (`S15_decision_list.md:89–91`), write a point of the region as `x = X_2 + r e`, with `e` a unit vector and `r_in ≤ r ≤ R`, and write `z = r/d` and `μ = e·D/d`. On the domain, `d > R > r_in > 0` (`S15_decision_list.md:120`), so `0 < z ≤ R/d < 1` throughout the region, and `−1 ≤ μ ≤ 1`.

- `|x − X_1| = d (1 − 2μz + z²)^{1/2}`.
- For every real `s`, `(1 − 2μz + z²)^{s/2} = Σ_{n≥0} g_n(s, μ) z^n`, where

  `g_n(s, μ) = Σ_{k=0}^{⌊n/2⌋} (−1)^k (−s/2)_{n−k} (2μ)^{n−2k} / (k! (n−2k)!)`,

  with `(c)_0 = 1` and `(c)_m = c(c+1)⋯(c+m−1)`. The series converges absolutely and uniformly for `0 ≤ z ≤ z_0` and `−1 ≤ μ ≤ 1`, for every `z_0 < 1`.
- For `−s/2 > −1/2` and `s ≠ 0`, `g_n(s, μ)` is the Gegenbauer polynomial `C_n^{(−s/2)}(μ)` (DLMF 18.5.10 and 18.12.4). Here `g_n` is defined by the formula above for every real `s`.

## What must be true

Everything under "What must be true" in `resume_s15_reading_repair_4.md`, `resume_s15_reading_repair_4d.md`, `resume_s15_reading_repair_4f.md`, `resume_s15_reading_repair_4g.md`, `resume_s15_reading_repair_4h.md` and `resume_s15_reading_repair_4i.md` stands, the 1-hour trial stop of `resume_s15_reading_repair_4i.md` included. Also:
- Item 4 of `S15_reading_decision_list.md` and its sentence "No limit or expansion is taken in place of an exact object" bind this trial. An attempt that uses the supplied series still prints the exact and complete conditions that item 4 names.
- The construction is yours, in each engine, and it uses the supplied expansion. Each engine constructs independently, and the `.wl` reads nothing from the SymPy side.
- No production source changes.

## What to do now, and where to stop

- Write the construction in each engine, in your work directory, and attempt `region/component_1/momentum_labels_{drain1, drain2}` with it, Wolfram first. Run through the guard, one job at a time, with the largest `--memory-gib` the guard admits beside the running SymPy attempt, and poll each job until it ends. Each attempt stops 1 hour after it started. If an attempt ends on an error in your code, you may fix it and run that object again; report every attempt.
- Then stop and write a short report. For each attempt, give the run, its start and end, whether it returned its conditions, ran out of memory, reached 1 hour or ended on an error in your code, its wall time and its peak memory. For an attempt that returned, list every branch its construction set aside without determining its conditions, if there was any, and say how the construction establishes that the printed conditions are complete. List the files you changed with their sha256. Do not attempt any other object, and do not start the final runs; the orchestrator will reply.

Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
