# Premises behind `S11c_b_89_sympy_perf_directive.md` (measurements)

Rule 2: the perf directive's claims about the bottleneck carry the evidence that produced them. This file
records the diagnosis; it does not state the withheld corrected basis count.

## Context — why this directive exists

The #89-repaired engine build (Codex-authored, 1206 insertions, compiles) applied the jet-tower repair
correctly, but the background run of the engine was killed by a spurious background-task death **after** the
edits and **before** the run completed. Resuming the run mechanically (`python3 <engine> > <.out>`) showed:
the four primary tasks completed and emitted (including the corrected basis, which reproduced the #86
reference — value withheld here), then the run went **silent for >1 h at 99.9% CPU** inside the control
tasks (`.out` frozen at 12,743,108 bytes; `cputime ≈ elapsed`, ratio 1.000 — a continuous computation, not a
hang). The run was killed on the user's approval; the partial primaries are preserved at
`~/.s11_build/S11c_b_89_partial_primaries.out` and the repaired engine at
`~/.s11_build/S11c_b_89_engine_postbuild.snapshot.py`.

## P1 — builds are not the bottleneck; the retained-grade projection is

Profiler `~/.s11_build/S11c_b_89_profile.py` → `.out`:

    [init] import+module-init: 14.6s
    [time] build_operator EUL (cache miss): 13.0s
    [time] build_operator EUL (cache hit): 0.0s

Then the profiler died (no traceback) during `retained_grade(op[0])` — the next call — consistent with the
projection being the heavy step. `build_operator` at 13 s and a working cache rule out the builders.

## P2 — retained_grade = uncached map of a heavyweight per-scalar series

    $ sed -n '881,882p;713,725p' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py

`retained_grade(value) = map_object(value, first_shape_series)` (L881). `first_shape_series` (L713) does
`sp.subs(PROFILE_GRADE_SUBS)` then `sp.series(expression, eta_bg, 0, 2).removeO()` then `sp.expand` then a
per-term grade filter — on every scalar. No memoization on either function (grep for a cache near them
returns nothing; the builder caches at L70–74 do not cover the projection).

## P3 — the controls call retained_grade dozens of times; the primaries prove one call is tractable

    $ grep -nc 'retained_grade' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py

The primaries (`task_energy_basis` L2952, `task_slab_operator` L2977, `task_coupling_kernel` L2990,
`task_admissibility` L3005) each wrap their emitted objects in `retained_grade`/`retained_energy_basis` and
**completed** in the resumed run — so one projection per object is tractable. The controls repeat it at
volume: `task_rep_invariance` (L3033) projects 3 objects × 2 branches (`BRANCHES` L50) × 2 reps
(`DENSITY_REPS` L51) × 2 routes; `task_independence` (L3077) projects base + corrupt variants;
`task_form_control` (L3162) and the four new controls (`HESSIAN_FREEZE`/`COEFFICIENT`/`JET_ZERO`/
`TOWER_DEPTH`, L3229–3480) each rebuild-and-project. The coupling kernel is the largest object and is
projected many times. Death by volume on an uncached heavyweight projection, made worse because the jet-tower
repair enlarged every expression.

## Decision-leg fold (2 legs, one pass, rule 7) + scope decision

Reviewed by two concurrent legs (grok 1.0.13 concurrency confirmed): Codex (4 defects,
`~/.s11_build/S11c_b_89_perf_directive_codex.txt`) and Grok (3 defects,
`~/.s11_build/S11c_b_89_perf_directive_grok.txt`), converging. All verified against the bytes and folded once:
1. (Grok) stale line numbers — the repair moved `first_shape_series`; cite by function name, not L713.
2. (Codex C1 + profiler) builds are *also* a hotspot (kernel miss 62.7 s / MATERIAL 525 s, 18 distinct
   variants) — so memoization alone cannot carry it.
3. (Grok G2 + Codex C1) projection ~300 s/call, controls project *unique* variants ⇒ per-call speed matters.
4. (Codex C2) an object-level cache on a partial `(branch,rep,route)` key would alias distinct control
   variants (ablation/corrupt/depth/full-zero) → a control compares a baseline with itself, a silent false
   pass ⇒ scalar-level memo, or full builder-key only.
5. (Codex C4) `retained_grade` takes unhashable `dict`s ⇒ memoize at the scalar level.
6. (Grok G3 + Codex C3) the `sp.series` replacement is risky — `e_W_bg ∉ PROFILE_GRADE_SUBS`, many `/W_bg^n`
   forms across routes, plus a series-failure fallback ⇒ equivalence checked against the untouched original
   across all forms, fallback for unrecognized shapes.

**User scope decision (2026-08-31): "projection fix first."** The directive is narrowed to the projection
only (kernel-build cost deferred); the goal is a bounded, completing run, not a minute target. **User also
asked to use multiple cores.** The machine has 16 logical / 8 physical cores; the projection maps a pure
function over independent scalars, so parallelizing it (multiprocessing, fork, order-preserving reassembly)
is byte-identical and became lever **B (primary)** alongside scalar memo (A); the risky `sp.series`
replacement is demoted to optional lever **C**. Per rule 7 this is one fold; next gate is the build's two
legs, not a re-leg.

## Withholding note

The corrected live basis count is in the partial primaries and the #86 record; it is NOT restated here or in
the directive. Per rule 12, blindness is by absence of a matching acceptance criterion in the directive and by
the projection-equivalence self-check (§4), not by hiding files. The only target the directive exposes is the
public frozen 26. The orchestrator will separately verify the fixed run's primaries are byte-identical to the
saved partial (basis unchanged) — that anchor is not handed to the builder.
