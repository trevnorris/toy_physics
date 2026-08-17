# Diagnosis: the Wolfram engine's strata machinery exhausts memory at two deterministic sites — name the computations and measure whether a bounded route exists

## The observables (all measured; do not take my word for anything else)

- Engine: `research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl` at
  HEAD `dfc2568c`, invoked as `wolframscript -file <wl> XKIN_ANISO <D>` (one cell, fresh kernel).
  The QE time wall was fixed in round 1 (bounded QE route at three sites); these deaths are
  POST-fix and are a different mechanism: the kernel's working set grows to ~29 GB and a memory
  guard kills it when available memory crosses 1 GiB.
- **D4 death, reproduced 3×, fully deterministic**: every run ends at exactly **emission 999**,
  last tag `WL_S11_XKIN_ANISO_D4_STRATUM3_ROOT2_N2_RANK_CHANGE_LOCUS_EQUATIONS`, guard-killed
  ~21 min in (available memory fell from ~29 GB free to <1 GiB in ~12 min after a stable start).
- **D2 death, reproduced 2×, fully deterministic**: every run ends at exactly **emission 1,411**,
  last tag `WL_S11_XKIN_ANISO_D2_STRATUM6_ROOT1_N2_NULLITY_CHANGE_LOCUS_REAL_STATUS_OPERANDS` —
  note this is a COMPLETED record block; whatever the kernel was inside when killed emitted
  nothing. Max inter-emission gap ~330 s; guard-killed ~22 min in.
- **D3 is the control: the SAME engine at the SAME package completes the ENTIRE cell — all
  strata — in 255 s, 2,501 emissions.**
- Death records and per-emission timelines (read freely):
  - `/home/trevnorris/.s11_build/fix1_build/guarded_cells_record.log` (CELL_START/CELL_MEMGUARD/CELL_END lines)
  - `/home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D4.out` + `.emission_times.tsv` + `.progress.log`
  - `/home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D2.out` + `.emission_times.tsv` + `.progress.log`
  - `/home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D3.out` + `.emission_times.tsv` (the completed control)
  - `.prerepair2` / `.killed_run1` variants of the same names are earlier runs of the same engine
    revision family — byte-level differences vs the current runs are provenance payloads only.
- The SymPy sibling engine (`research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py`,
  record `research/pde_ledger_v3/scripts/out/S11_stray_longitudinal_sympy_audit.out`) COMPLETED
  all 21 cells including XKIN_ANISO D2 and D4. You may read both engines and both records freely.

## Your task

1. **Name both computations.** For each death site (D4 and D2), from the engine source and the
   emitted-tag order, identify the exact expression the kernel was evaluating when killed —
   file:line, the function, the operand objects, the variables. For D4 the last tag is an
   `_EQUATIONS` emission; state what the code does next before the following tag would print.
   For D2 the last tag closes a record block; walk the code order and name the NEXT expression —
   naming it is job #1 for that site, because nothing after it left any trace. State the remaining
   emission queue at each site so what's missing is explicit.
2. **Measure the operands.** Reconstruct, in Python/SymPy, the operand objects entering each
   named expression — the stratum-restricted systems: number of equations, number of unknowns,
   polynomial degrees, term counts, radical content. Compute the SAME quantities at D3 for the
   corresponding strata so D3→D4 and D3→D2 growth is a computed ratio, not a guess. D2 dying
   where D4 does not (and vice versa) is part of what your measurements must explain or explicitly
   fail to explain.
3. **Bounded-vs-runaway, per site.** From the emission timeline and the memory samples in
   `.progress.log`: is each death one expensive call or an accumulation across many? If the former,
   what quantity explodes inside it — measured, from your reconstruction.
4. **Name the object a bounded route would compute.** If the same mathematical object each dying
   expression is after (state exactly which object) is obtainable by a route whose intermediates
   stay bounded, demonstrate feasibility on the RECONSTRUCTED operands in SymPy with literal
   timing and size output. Name the OBJECT, not a code patch. If successor stages in the remaining
   queue share the same growth property, measure that too.

## Constraints — these bind you

- ⛔ **Do NOT launch any Mathematica kernel or `wolframscript` process.** All computation runs in
  Python/SymPy (available on this box). This is a 2-seat licence machine and memory is the
  resource under study.
- ⛔ Never modify the working tree. Create your own uniquely-named scratch directory under /tmp
  (`mktemp -d /tmp/s11_wall2_XXXX`), write scripts and outputs ONLY there, and report its path.
- ⭐ A prose derivation is worth nothing: for every measured claim, save the script AND its
  literal stdout to named absolute paths and report those paths. Claims without them will be
  discarded.
- Wrap every script run in `timeout 600`. A timeout hit is itself a data point — report it, do
  not raise it. This machine has 30 GB total; keep your own probes' memory modest and print
  observable progress from anything long — a silent long run is the failure mode under study,
  not a method.

## Report format

1. Per site: the named computation (file:line, operands, variables) and the remaining queue.
2. Measured operand table: D2 / D3 / D4 for each quantity, with script paths.
3. Per site: one-call-vs-accumulation verdict with the timeline/memory evidence.
4. The object(s) a bounded route would compute, with the feasibility measurement.
