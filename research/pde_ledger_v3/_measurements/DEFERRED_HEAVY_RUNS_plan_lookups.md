# Lookups — the 64 GB run plan in DEFERRED_HEAVY_RUNS.md (generated 2026-10-18 10:48)

Generator: `_scratch/backlog_plan/gen_lookups.sh` (sed/grep/git log/free only).

## The guard's memory ceiling
```
$ sed -n '10,12p;73,79p' scripts/s11c_guarded_run.py
checkouts. The registry uses that user's account-database home, independent
of $HOME, matching the scope of systemctl --user. The budget is
MemTotal rounded down to whole GiB, and each job takes one logical CPU
    for line in Path('/proc/meminfo').read_text().splitlines():
        if line.startswith('MemTotal:'):
            return int(line.split()[1]) * 1024 // GIB * GIB
    raise RuntimeError('MemTotal unavailable; refusing launch')


POOL_MEMORY_MAX = physical_memory()
```

```
$ free -g | head -2
               total        used        free      shared  buff/cache   available
Mem:              31           7           1           0          24          23
```

## Last commits of the scripts the backlog runs
```
$ git log --format='%h %ad %s' --date=short -1 -- research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
72d6fe0b 2026-10-05 Retain S11c equations, governing contracts and qualified scattering evidence
```

```
$ git log --format='%h %ad %s' --date=short -1 -- research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl
06048d15 2026-09-03 S11c-b P1-WL: WL S11CB_PRIMARIES_ONLY emit gate (fresh cross-engine .out) — leg-gated
```

```
$ git log --format='%h %ad %s' --date=short -1 -- research/pde_ledger_v3/scripts/S11c_c1_cross_engine_comparator.py
704308af 2026-09-05 S11c-c1 T7 comparator REPAIR applied + re-reviewed CLEAN (astra build)
```

## #89 and #90 name the same command
```
$ grep -n -A4 'Re-run command (full in-band' research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md | grep -n 'python3\|wolframscript\|S11CB_PRIMARIES_ONLY'
4:101-  python3 research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py \
10:126-  python3 research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py \
16:157-  wolframscript -file research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl \
```

## Memory figures quoted from the entries
```
$ grep -n 'OOMs ~3.5 h in\|~16 GB/case\|peak RSS ~301 MB' research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md
105:  ⚠ On the 30 GB box this OOMs ~3.5 h in, still inside the first control (the MATERIAL coupling-kernel
143:  un-frozen operator must hold its **full background jet tower un-reduced (~16 GB/case)** right up to the
161:  ⚠ The operator loop alone peaks ~16 GB/case on 30 GB (reaped there); needs ≥64 GB for the full 4-case run.
176:- **What ran on 30 GB (light+mid tier, 46 of 50 families):** peak RSS ~301 MB; the two engines AGREE — every
```

## P2a and P2b: decision lists committed, no build yet
```
$ git log --format='%h %ad %s' --date=short -- research/pde_ledger_v3/directives/S11c_b_p2b_gamma_bridge_directive.md research/pde_ledger_v3/directives/S11c_b_p2a_slab_row_join_directive.md
69bcc18d 2026-09-03 S11c-b P2b: §3a coefficient-bridge decision list (folded from 2 decision legs)
06048d15 2026-09-03 S11c-b P1-WL: WL S11CB_PRIMARIES_ONLY emit gate (fresh cross-engine .out) — leg-gated
```

```
$ git log --format='%h %ad %s' --date=short --all -i --grep='P2b\|P2a' | head -8
bcb9f7d7 2026-09-03 S11c-b step record (the close) — per-engine verified, cross-engine residual deferred (2 legs folded)
66e8d021 2026-09-03 S11c-b cross-engine residual DEFERRED to >=64 GB (STEP 0 overturned); S11c-b closes on per-engine verification
f04bdcf9 2026-09-03 STATUS: cross-engine integration pass in flight — P1-WL gate + both P2 decision lists committed
69bcc18d 2026-09-03 S11c-b P2b: §3a coefficient-bridge decision list (folded from 2 decision legs)
06048d15 2026-09-03 S11c-b P1-WL: WL S11CB_PRIMARIES_ONLY emit gate (fresh cross-engine .out) — leg-gated
```

## S11c-c2: no run scheduled
```
$ grep -n 'does not schedule a run' research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md
90:This entry carries the debt; it does not schedule a run or predict its result.
```

## Review findings: the sources behind each disposition
```
$ sed -n '251,253p;262,263p' scripts/s11c_guarded_run.py
    if (requested < 1024**3 or requested > POOL_MEMORY_MAX or requested % (1024**3)
            or reserved + requested > POOL_MEMORY_MAX):
        raise RuntimeError('aggregate physical-memory reservation budget exceeded')
    if available < AVAILABLE_MIN + headroom + requested:
        raise RuntimeError('available memory does not cover outstanding commitments plus this job')
```

```
$ sed -n '369,370p' scripts/s11c_guarded_run.py
    if spec.get('pool') is not None:
        resource.setrlimit(resource.RLIMIT_AS, (spec['memoryMax'], spec['memoryMax']))
```

```
$ grep -n 'outer_operator_algebra_signature\|exact_rational_canon_deferred_for_oversized_leaf' research/pde_ledger_v3/scripts/S11c_c1_cross_engine_comparator.py
607:            note="outer_operator_algebra_signature; raw_operand_emitted"
1948:            "exact_rational_canon_deferred_for_oversized_leaf",
```

```
$ grep -n 'G3 — Launch legs on sight' CLAUDE.md
161:**G3 — Launch legs on sight, before inspecting the result; no commit before both report.** *(was R8, R9)*
```

```
$ grep -n 'Do not rerun the full SymPy harness\|run it again with .--memory-gib 18' _scratch/s9b_build/s9b_repair_build_review_prompt_r10.md
122:- **Do not rerun the full SymPy harness** (about nine hours); read its output named above. You may rerun the
132:  If the guard kills a SymPy run at its memory limit, run it again with `--memory-gib 18`; the guard admits that
```

```
$ git diff --stat -- research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md
 research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md | 47 +++++++++++++++++++++++++++
 1 file changed, 47 insertions(+)
```

