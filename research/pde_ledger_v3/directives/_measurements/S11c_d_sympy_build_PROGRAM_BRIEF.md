# Measurements grounding `S11c_d_sympy_build_PROGRAM_BRIEF.md` (rule 2, orchestrator half)

Regenerated from the commands below (⛔ not transcribed). The brief is a POINTER/build-mandate directive; its
load-bearing artifact claims are grounded here + in the already-committed census
`directives/_measurements/S11c_d_sympy_build_directive_census.md` (the 3-parent fold: base 2441 + c1 44 + c2 70
→ 2555) and the builder report. ⚠ The prototype / `.out` are astra-live (user-driven build in flight) — the
state below is a point-in-time snapshot at prep.

## SHAs the brief relies on
$ git show -s --format="%h %s" 399a8516 64b85989
399a8516 S11c-d SHARED PHYSICS spec v10 — CLEARED (round-10 gate: BOTH legs SOUND) + post-clearance wording polish
64b85989 S11c-d SymPy build directive CLEARED — folded 4 round-1 findings + round-2 legs both SOUND

## F1 / the 5-slot payload (each = 8 = 4 (α,ρ) cases × 2 closed rows) — stable c2 exports
$ grep -oE "'(FOURIER_PROFILE_BINDINGS|COMPUTED_BRANCH_BINDINGS|VALUE|MULTIGRADE|DIMENSION_L_T_M)'" scripts/S11c_c2_exports.py | sort | uniq -c
      8 'COMPUTED_BRANCH_BINDINGS'
      8 'DIMENSION_L_T_M'
      8 'FOURIER_PROFILE_BINDINGS'
      8 'MULTIGRADE'
      8 'VALUE'

## Build state snapshot (astra-live, point-in-time): prototype exists, exports absent, work outstanding
$ ls -la scripts/S11c_d_mixing_scattering_sympy_audit.py scripts/S11c_d_exports.py 2>&1
ls: cannot access 'scripts/S11c_d_exports.py': No such file or directory
-rw-rw-r-- 1 trevnorris trevnorris 108581 Sep 10 09:45 scripts/S11c_d_mixing_scattering_sympy_audit.py

$ grep -ac "OUTSTANDING_CONSTRUCTIONS" scripts/out/S11c_d_mixing_scattering_sympy_audit.out   # tag present ⇒ work remains
2

## The cleared authorities the brief points at (exist)
$ ls directives/S11c_d_SHARED_PHYSICS.md directives/S11c_d_sympy_build_directive.md
directives/S11c_d_SHARED_PHYSICS.md
directives/S11c_d_sympy_build_directive.md
