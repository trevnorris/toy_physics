# Measurements — S9b spec amendment 1 and directive amendment 4, review round 0 (generated 2026-10-09 10:02)

Generator: `_scratch/s9b_build/gen/s9b_spec_amend1_r0_lookups.sh` (sha256sum/grep only).

## The builder's stop left its engine, harness and delta as reviewed in r0
```
$ sha256sum /var/projects/toy_physics-s9b-r2-py/research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py /var/projects/toy_physics-s9b-r2-py/research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py /var/projects/toy_physics-s9b-r2-py/research/pde_ledger_v3/scripts/S9b_exports.py
16b00f137623c897a5f890c37fa5ea0082e8c1156c9907409b2884603a6290d9  /var/projects/toy_physics-s9b-r2-py/research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
be3c83541572516228299cbdc58bc5cbd336d58ce2c2c9dcf680a1889da7c5fd  /var/projects/toy_physics-s9b-r2-py/research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
3f3f9ba07f52594aa803193f121db1ef9fc3e91e606123e0abe0971cd4c41b78  /var/projects/toy_physics-s9b-r2-py/research/pde_ledger_v3/scripts/S9b_exports.py
```

```
$ grep -E 'S9b_light_bending_sympy_audit.py|S9b_light_bending_sympy_ablation.py|S9b_exports.py' _scratch/s9b_build/s9b_repair_build_review_baseline_r0.sha256
16b00f137623c897a5f890c37fa5ea0082e8c1156c9907409b2884603a6290d9  research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
3f3f9ba07f52594aa803193f121db1ef9fc3e91e606123e0abe0971cd4c41b78  research/pde_ledger_v3/scripts/S9b_exports.py
be3c83541572516228299cbdc58bc5cbd336d58ce2c2c9dcf680a1889da7c5fd  research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
```

```
$ grep -n 'stopped under item 12\|engine, harness, and delta are unchanged' _scratch/s9b_build/s9b_py_amend3_stop_report.md
3:**Incomplete; stopped under item 12 after a second method failure.** I have not obtained the required general-profile Part B radar reduction. The guarded [probe](../../scripts/S9b_general_profile_log_probe.py) constructs the speed-only dispersion and first optical variation of the round-trip action, then tests endpoint-dilatio
18:**Deliverable status:** engine, harness, and delta are unchanged and retain their previous source pins. The vocabulary, general-profile engine, vector-dispersion repair, gate repair, and requested final demonstrations/full harness were **not completed** before the stop. No new export was published. The inventory enumerates al
```

## The object the builder could not extract, as the accepted texts named it
```
$ git show 97580010:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md | grep -n 'coefficient of `ln(1/b²)`\|logarithmic component'
172:   uses only the coefficient of `ln(1/b²)` in the regime `Z_E, Z_R ≫ b`, for every profile, including tails
173:   other than `1/r`. The radar claim is limited to this logarithmic component.
217:  - An effective `γ` from `Δθ`, and one from the coefficient of `ln(1/b²)` in the round-trip excess time for
```

```
$ git show 3699f4ef:research/pde_ledger_v3/directives/S9b_repair_build_directive.md | grep -n 'ln(1/b²)'
72:   - `<OBS>`: `DEFLECTION` (Δθ), or `RADAR` (the coefficient of `ln(1/b²)` in the round trip, item 7).
113:7. **The radar coefficient (D6 B3).** The coefficient of `ln(1/b²)` is obtained by computation from Part A's
133:      - Print the deflection, the round trip's `ln(1/b²)` coefficient (item 7), both effective `γ`s, their
145:    - For the deflection and for the radar `ln(1/b²)` coefficient, print each restriction below of:
207:    - **K6, radar source (a named coefficient knife):** the object from which the `ln(1/b²)` coefficient is
293:  - The radar `γ` is solved on every stratum of the `ln(1/b²)` coefficient. The `γ` difference, and each
```

## The reviewed versions are the working-tree files
```
$ sha256sum research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md research/pde_ledger_v3/directives/S9b_repair_build_directive.md _scratch/s9b_build/S9b_SHARED_PHYSICS_amend1_reviewed_r0.md _scratch/s9b_build/S9b_repair_build_directive_amend4_reviewed_r0.md
215ae8d7017fac46927db8064ce232639749b957097c68aebd49f76ca08deead  research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
77869b472754834ae0c609869a4d454d2f61def064c89ef05e2b838a5b1c41b4  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
215ae8d7017fac46927db8064ce232639749b957097c68aebd49f76ca08deead  _scratch/s9b_build/S9b_SHARED_PHYSICS_amend1_reviewed_r0.md
77869b472754834ae0c609869a4d454d2f61def064c89ef05e2b838a5b1c41b4  _scratch/s9b_build/S9b_repair_build_directive_amend4_reviewed_r0.md
```

```
$ cat _scratch/s9b_build/s9b_spec_amend1_review_baseline_r0.sha256
215ae8d7017fac46927db8064ce232639749b957097c68aebd49f76ca08deead  research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
77869b472754834ae0c609869a4d454d2f61def064c89ef05e2b838a5b1c41b4  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

## No text in either document still names the old object
```
$ grep -c 'ln(1/b²)' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md research/pde_ledger_v3/directives/S9b_repair_build_directive.md
research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md:2
research/pde_ledger_v3/directives/S9b_repair_build_directive.md:0
```

```
$ grep -n 'ln(1/b²)' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md research/pde_ledger_v3/directives/S9b_repair_build_directive.md
research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md:179:   𝒮_RT(b) ≡ lim_{Z_E, Z_R → ∞} ∂Δt_RT(b; Z_E, Z_R) / ∂ln(1/b²) ,      Z_E and Z_R held fixed in the derivative.
research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md:182:   This applies to every profile, including tails other than `1/r`. Where `Δt_RT` is `A·ln(1/b²)`, plus terms
```

## Both verdicts
```
$ grep -n -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_spec_amend1_review_r0_codex_final.txt _scratch/s9b_build/s9b_spec_amend1_review_r0_grok.txt
_scratch/s9b_build/s9b_spec_amend1_review_r0_codex_final.txt:1:Verdict: CLEAR
_scratch/s9b_build/s9b_spec_amend1_review_r0_grok.txt:1:Verdict: CLEAR
```

```
$ grep -n 'Findings:' _scratch/s9b_build/s9b_spec_amend1_review_r0_codex_final.txt
3:**Findings:** None within the amendments’ Part B comparison scope.
```

```
$ grep -n -o 'No finding changes what the engines compute[^.]*' _scratch/s9b_build/s9b_spec_amend1_review_r0_grok.txt
1:No finding changes what the engines compute, how the outputs are paired, or what may be claimed
```

