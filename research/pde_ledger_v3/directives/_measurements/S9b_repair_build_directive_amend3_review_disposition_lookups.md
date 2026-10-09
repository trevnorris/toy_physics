# Measurements — S9b repair build directive amendment 3, review round 0 (generated 2026-10-09 08:50)

Generator: `_scratch/s9b_build/gen/s9b_amend3_r0_lookups.sh` (sha256sum/grep/sed only).

```
$ sha256sum _scratch/s9b_build/S9b_repair_build_directive_amend3_reviewed_r0.md
16f908a286353facf17d79090ad6994cefd7cfe1032ec451094dee7d91429886  _scratch/s9b_build/S9b_repair_build_directive_amend3_reviewed_r0.md
```

```
$ cat _scratch/s9b_build/s9b_repair_build_directive_amend3_review_baseline_r0.sha256
16f908a286353facf17d79090ad6994cefd7cfe1032ec451094dee7d91429886  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ grep -n -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend3_review_r0_codex_final.txt _scratch/s9b_build/s9b_amend3_review_r0_grok.txt
_scratch/s9b_build/s9b_amend3_review_r0_codex_final.txt:1:Verdict: NEEDS REVISION
_scratch/s9b_build/s9b_amend3_review_r0_grok.txt:3:Verdict: NEEDS REVISION
```

## C1: the spec's unordered differences and the vocabulary's names for them
```
$ grep -n 'half the' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
174:3. **The two one-way excess times** between the same endpoints, and their **nonreciprocal part** (half the
```

```
$ grep -n 'difference between the two' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
225:  - The difference between the two `γ`s, printed as an object.
```

```
$ grep -n 'A_NONRECIPROCAL_G\|B_GAMMA_DIFFERENCE' _scratch/s9b_build/S9b_repair_build_directive_amend3_reviewed_r0.md | head -4
83:   | `A_NONRECIPROCAL_G<abc>` | their nonreciprocal part, at that grade |
86:   | `B_GAMMA_<OBS>`, `B_GAMMA_DIFFERENCE` | Part B's effective `γ`s and their difference |
94:     `B_GAMMA_<OBS>`, `B_GAMMA_DIFFERENCE`, `B_RESIDUAL_<OBS>`, `B_CONDITION_<OBS>` and `B_IMPLIED_JN_<OBS>`.
98:     `BRANCH_EXISTENCE`, `PATH_TRAVERSAL`, `BRANCH_TYPE`, every `A_` name, `B_GAMMA_<OBS>`, `B_GAMMA_DIFFERENCE`,
```

## G1: item 13's general-profile sentence against item 10's substitutions
```
$ grep -n 'Every shared-vocabulary object' -A 1 _scratch/s9b_build/S9b_repair_build_directive_amend3_reviewed_r0.md
182:    - Every shared-vocabulary object (item 4) is computed for general radial profiles: `δ`, `V`, `ξ_w`, and,
183-      where they enter, `ρ_br` and `f`.
```

```
$ grep -n 'flow only:\*\*\|speed only:\*\*\|tilt only:\*\*\|first with' _scratch/s9b_build/S9b_repair_build_directive_amend3_reviewed_r0.md
128:      - Substitute that `V` into Part A's observables, first with `δ ≡ 0`, `ξ_w ≡ 0`, then with `δ` and `ξ_w`
150:      - **flow only:** `δ ≡ 0`, `ξ_w ≡ 0`, `V` live;
151:      - **speed only:** `V ≡ 0`, `ξ_w ≡ 0`, `δ` live;
152:      - **tilt only:** `δ ≡ 0`, `V ≡ 0`, `ξ_w` live.
```

```
$ grep -n 'not a stop event' _scratch/s9b_build/S9b_repair_build_directive_amend3_reviewed_r0.md
159:      that label only. Building this item is not a stop event under item 12, and it does not replace the live `j_n`
```

## G2: the forward prefixes reach item 8's object; item 10's forward print list
```
$ grep -n 'every `A_` name' _scratch/s9b_build/S9b_repair_build_directive_amend3_reviewed_r0.md
98:     `BRANCH_EXISTENCE`, `PATH_TRAVERSAL`, `BRANCH_TYPE`, every `A_` name, `B_GAMMA_<OBS>`, `B_GAMMA_DIFFERENCE`,
```

```
$ grep -n 'A_NONRECIPROCAL_PATH_DEPENDENCE' _scratch/s9b_build/S9b_repair_build_directive_amend3_reviewed_r0.md
84:   | `A_NONRECIPROCAL_PATH_DEPENDENCE` | item 8's object |
```

```
$ grep -n 'Print the deflection, the round trip' -A 2 _scratch/s9b_build/S9b_repair_build_directive_amend3_reviewed_r0.md
132:      - Print the deflection, the round trip's `ln(1/b²)` coefficient (item 7), both effective `γ`s, their
133-        difference and both Part B residuals against the references, each as a function of `b`. Also print the round-trip excess time, both
134-        one-way excess times and their nonreciprocal part, with `b`, `Z_E` and `Z_R` kept.
```

## G3: the Part 3 gamma sentence against the r0 disposition's G3 outcome
```
$ grep -n 'Each effective' -A 1 _scratch/s9b_build/S9b_repair_build_directive_amend3_reviewed_r0.md
290:  - Each effective `γ` is solved on every stratum of its observable's coefficient. The `γ` difference and every
291-    restriction and forward copy use that solution.
```

```
$ grep -o 'The radar `γ` is solved on every stratum[^|]*' research/pde_ledger_v3/directives/_measurements/S9b_repair_build_r0_review_disposition.md
The radar `γ` is solved on every stratum of the coefficient. The γ difference and every restriction and forward copy use it. 
```

```
$ sed -n '216,218p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
- **Part B.**
  - An effective `γ` from `Δθ`, and one from the coefficient of `ln(1/b²)` in the round-trip excess time for
    `Z_E, Z_R ≫ b`.
```

