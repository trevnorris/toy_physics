# Independent opinion — the S11 SymPy engine will not publish its export

Read-only. Change no file. This asks for your independent read of a problem, not confirmation of mine.

## What exists

Codex built `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py` from
`research/pde_ledger_v3/directives/S11_sympy_build_directive.md` (closed, 5 rounds / 10 legs). The physics
authority is `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`. `CLAUDE.md` binds.

The script is 1703 lines, compiles, and emits genuine SymPy `srepr` objects. It declares 22 cells:

```
$ sed -n '116,127p' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
PACKAGE_DIMS = {
    "MAIN": (2, 3, 4, 5),          "XFORM_CURLONLY": (2, 3, 4, 5),
    "XFORM_EXTRA": (2, 3, 4, 5),   "XFORM_DIVONLY": (3, 4),
    "XFORM_TRACELESS": (3, 4),     "XCOEF_BSCALE": (3,),
    "XCOEF_BSIGN": (3,),           "XKIN_ANISO": (2, 3, 4, 5),
}
PRIMARY_PACKAGE = "MAIN"
```

`MAIN` is the answer; the seven `X*` packages are form and coefficient controls.

## Four runs, none published `S11_exports.py`

| run | conditions | outcome |
|---|---|---|
| 1 | Codex's own build run, no cap | died ~36 min, mid-patch, no `tokens used` line |
| 2 | arbiter re-run, no cap | **22.3 GB RSS at 100% CPU**, OOM-killed at 856 tags |
| 3 | arbiter re-run, no cap | same, 856 tags |
| 4 | `ulimit -v 8000000` | **21 of 22 cells**, all of `MAIN` D2–D5, then **silent for ~3.5 h** at 99.7% CPU while RSS grew 208 → 459 MB. Stuck inside `XKIN_ANISO_D4`. Killed at 5h37m. |

Every OOM kill is a `SIGKILL`, so every one of these looked like a clean exit with empty stderr. That
misled me twice before I measured the python process rather than its shell wrapper.

Run 4's stdout is preserved at `/home/trevnorris/.s11_build/s11_stalled_FINAL.out` (14.9 MB, 5123 tags).

## Where the volume goes — measured on that file

```
$ awk -F: '/^PY_S11_/{...}' s11_stalled_FINAL.out    # bytes and tags per package
XKIN_ANISO      9.22 MB   737 tags        MAIN            0.97 MB   813 tags
XFORM_EXTRA     1.19 MB   813 tags        XFORM_TRACELESS 0.71 MB   503 tags
XFORM_CURLONLY  1.07 MB  1001 tags        XFORM_DIVONLY   0.50 MB   503 tags
XCOEF_BSCALE    0.28 MB   254 tags        XCOEF_BSIGN     0.24 MB   254 tags

$ # three largest SINGLE tags
0.82 MB  PY_S11_XKIN_ANISO_D4_ROOT_COINCIDENCE_R2_R3_K_REAL_ADMISSIBLE
0.55 MB  PY_S11_XKIN_ANISO_D3_ROOT_COINCIDENCE_R2_R3_K_REAL_ADMISSIBLE
0.41 MB  PY_S11_XKIN_ANISO_D4_ROOT_COINCIDENCE_R2_R3_K_INCONSISTENT

$ # MAIN by dimension
D2 249 tags   D3 254 tags   D4 249 tags   D5 61 tags
```

## Structural facts

```
$ sed -n '1678,1687p' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
    main_completed = {(package, n) for package, n in completed_pairs if package == PRIMARY_PACKAGE}
    main_declared = {(PRIMARY_PACKAGE, n) for n in PACKAGE_DIMS[PRIMARY_PACKAGE]}
    if main_completed == main_declared:
        ledger = merged_export(main_dim_data, run_pairs_payload, skipped_pairs_payload)
        write_exports(ledger)
    else:
        stale = SCRIPT_DIR / "S11_exports.py"
        if stale.exists(): stale.unlink()
        ISSUES.append("S11_exports.py not published because a declared MAIN cell did not complete")
```

- `write_exports` runs **after the whole package loop**, though its condition covers `MAIN` only.
- Bare `except Exception` at `:756 :775 :955 :1317 :1450 :1668`. `:1668` wraps `run_cell` in the main loop.
- The script has **no `out/` writing**; it prints to stdout and the caller pipes it.
- `S11_exports.py` has never existed.

## Three changes the user directed

1. Write the ledger as soon as `MAIN` completes, rather than after the whole loop.
2. Pipe stdout to `research/pde_ledger_v3/scripts/out/S11_stray_longitudinal_sympy_audit.out`.
3. Change the expensive calculation itself.

I had proposed instead lowering the memory cap to ~2 GB so the swelling cell would raise `MemoryError`,
be caught at `:1668`, and let the loop reach publication. The user rejected that as forcing a failure
rather than removing one. Say if you think that judgement was wrong.

## What I want your opinion on

1. **The stall.** What is actually happening in `XKIN_ANISO_D4` root coincidence, and is the fix a
   different formulation of the same object or a bound on it? `CLAUDE.md` rule 3 says name the object and
   never the recipe, so be careful about what you are proposing to change.
2. **`MAIN` D5 emitted 61 tags against ~250 at D2/D3/D4.** Is that genuine, or a swallowed `MemoryError`
   inside one of those bare handlers with the loop still counting the cell complete? What computation
   would distinguish them?
3. **Moving `write_exports` earlier.** Safe, or does it break something — `F6`'s publish semantics, the
   `run_pairs`/`skipped_pairs` payloads it consumes, or the deletion of a stale export?
4. **The bare `except Exception` handlers.** Which are legitimate and which mask a real failure?
5. **Anything the diagnosis gets wrong**, including the premise that this is fixable in-place rather than
   a spec-level `§10` unavailable-construction outcome.

## Method

Every claim carries the command that produced it and its literal output; a claim with no command behind it
will be discarded (`CLAUDE.md` rule 2). Do not propose new machinery, registries or checklists. Do not
edit any file. Report a finding only if it bears on whether the physics could come out wrong or whether a
wrong result could survive.
