# Measurements — S9b repair build directive amendment 3, review round 1 (generated 2026-10-09 09:11)

Generator: `_scratch/s9b_build/gen/s9b_amend3_r1_lookups.sh` (sha256sum/grep only).

## The reviewed version is the working-tree directive
```
$ sha256sum research/pde_ledger_v3/directives/S9b_repair_build_directive.md _scratch/s9b_build/S9b_repair_build_directive_amend3_reviewed_r1.md
0282e4e406d4ab685d077980630a045d660030b0b424f446ad8e5aea11032fe1  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
0282e4e406d4ab685d077980630a045d660030b0b424f446ad8e5aea11032fe1  _scratch/s9b_build/S9b_repair_build_directive_amend3_reviewed_r1.md
```

```
$ cat _scratch/s9b_build/s9b_repair_build_directive_amend3_review_baseline_r1.sha256
0282e4e406d4ab685d077980630a045d660030b0b424f446ad8e5aea11032fe1  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

## Both verdicts
```
$ grep -n -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend3_review_r1_codex_final.txt _scratch/s9b_build/s9b_amend3_review_r1_grok.txt
_scratch/s9b_build/s9b_amend3_review_r1_codex_final.txt:1:Verdict: CLEAR
_scratch/s9b_build/s9b_amend3_review_r1_grok.txt:1:Verdict: CLEAR
```

```
$ grep -n 'Findings:\|^No findings' _scratch/s9b_build/s9b_amend3_review_r1_codex_final.txt _scratch/s9b_build/s9b_amend3_review_r1_grok.txt
_scratch/s9b_build/s9b_amend3_review_r1_codex_final.txt:3:**Findings:** None within the requested scope.
_scratch/s9b_build/s9b_amend3_review_r1_grok.txt:3:No findings. The amendment gives each Parts A–C, branch-existence, and item 7, 8, and 10 object one shared name, keeps every profile those objects leave live as a general radial function, and its repair sentences match the round-0 disposition without stating a value.
```

## The round-0 repairs are in the reviewed text
```
$ grep -n 'emitter-to-reflector time minus\|deflection `γ` minus the radar' _scratch/s9b_build/S9b_repair_build_directive_amend3_reviewed_r1.md
83:   | `A_NONRECIPROCAL_G<abc>` | their nonreciprocal part, at that grade: half of (emitter-to-reflector time minus reflector-to-emitter time) |
86:   | `B_GAMMA_<OBS>`, `B_GAMMA_DIFFERENCE` | Part B's effective `γ`s, and their difference: the deflection `γ` minus the radar `γ` |
```

```
$ grep -n 'each profile the object leaves live\|apply their substitutions first\|it omits the object and reports it under item 17' _scratch/s9b_build/S9b_repair_build_directive_amend3_reviewed_r1.md
183:    - In every shared-vocabulary object (item 4), each profile the object leaves live is a general radial
185:      and forward stages, and Part C's responses, apply their substitutions first.
191:      an item 10 condition, that is not a stop event (item 10): it omits the object and reports it under item 17.
```

```
$ grep -n 'every graded `A_` name' _scratch/s9b_build/S9b_repair_build_directive_amend3_reviewed_r1.md
98:     `BRANCH_EXISTENCE`, `PATH_TRAVERSAL`, `BRANCH_TYPE`, every graded `A_` name (`A_…_G<abc>`),
```

```
$ grep -n 'solved on every stratum of the `ln(1/b²)` coefficient' -A 1 _scratch/s9b_build/S9b_repair_build_directive_amend3_reviewed_r1.md
293:  - The radar `γ` is solved on every stratum of the `ln(1/b²)` coefficient. The `γ` difference, and each
294-    restriction and forward copy of the radar `γ`, use that solution.
```

```
$ grep -c 'Each effective' _scratch/s9b_build/S9b_repair_build_directive_amend3_reviewed_r1.md
0
```

