# Diagnosis: a live Wolfram engine cell has been silent 15+ hours — name the computation and measure whether it is bounded

## The observables (all measured; do not take my word for anything else)

- Engine: `research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl`,
  invoked as `wolframscript -file <wl> XKIN_ANISO 4` (one cell, fresh kernel).
- The kernel (PID 2595102) has run **20h19m elapsed with 20h19m CPU** — a 100% duty cycle,
  never idle. It is STILL RUNNING. Do not touch it.
- RSS history: ~13.3 GB (yesterday afternoon) → ~9.6 GB (last night) → **3.7 GB now**.
  Memory is shrinking, not growing.
- Output file: `research/pde_ledger_v3/mathematica/out/S11_stray_longitudinal_mathematica_audit.out`
  (live, append-only; read it freely). Last three tags, emitted 16:58 yesterday, then silence:
  `WL_S11_XKIN_ANISO_D4_ROOT2_STACKED_DROP_K_EQUATIONS` / `..._K_SOLUTION` / `..._K_IDENTICALLY_SATISFIED`.
- Earlier same-cell timeline: `..._D4_ROOT2_RANK_DROP_JOINT_*` emitted ~14:19; the cell started 11:12.
- The SAME engine at the SAME package completed D3 — the entire cell, all stages — in **285 s**.
- The same engine at D2 was killed by a memory guard (>28 GB from a cold kernel).

## Your task

1. **Name the computation in progress.** From the engine source and the last-emitted tag, identify
   the exact expression being evaluated right now — file:line, the function, the equation set, and
   the solve variables. State what tags come after it in emission order for this cell, so the
   remaining queue is explicit.
2. **Measure bounded-vs-runaway.** Do not assert; compute. Reconstruct the operand objects for the
   in-progress stage at D4 in an independent tool (SymPy/Python) — the matrix, its stacked form,
   the maximal minors, their polynomial degrees and term counts in the relevant variables — and
   from those measured sizes, characterize the computation the engine is attempting. Compare the
   same measured quantities at D3 (where the cell took 285 s) so the D3→D4 growth is a computed
   ratio, not a guess. If you can bound the work, show the bound's operands. If you cannot, show
   what quantity explodes and its measured growth.
3. **Say what would change the outcome.** If the stage is a finite grind that will complete, give
   the evidence. If it will not complete in reasonable time, identify what property of the
   computation makes it so, and what a repair would have to compute instead — name the OBJECT a
   faster route would produce, not a code patch. Note whether the successor stages
   (the remaining queue from item 1) share the same property, measured the same way.

## Constraints — these bind you

- ⛔ **Do NOT launch any Mathematica kernel or `wolframscript` process.** The licence has two
  seats and one is held by the live run you must not disturb. All computation you do runs in
  Python/SymPy (available on this box).
- ⛔ Do not kill, signal, trace, or attach to PID 2595102 or any Wolfram process.
- ⛔ Never modify the working tree. Write scripts and outputs ONLY under a directory you create
  in /tmp.
- ⭐ A prose derivation is worth nothing: for every measured claim, save the script AND its
  literal stdout to named absolute paths and report those paths. Claims without them will be
  discarded.
- Wrap every script run in `timeout 600`. A timeout hit is itself a data point — report it, do
  not raise it.
- Long-running scripts must print observable progress; a silent long run is the failure mode
  under study, not a method.

## Report format

1. The computation in progress (file:line, equations, variables) and the remaining queue.
2. Measured operand table: D3 vs D4 for each quantity you computed, with script paths.
3. Verdict: bounded (evidence) / unbounded-in-practice (the exploding quantity, measured).
4. What a faster route would compute, as objects.
