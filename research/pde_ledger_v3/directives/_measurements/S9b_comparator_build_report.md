- Files: `research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py` and `test_S9b_cross_engine_comparator.py`. All six repair findings addressed; 330 joined rows, every original stratum/profile/point retained. No residual target was given; no scientific interpretation or commits. Reused O2 lossless readers and the existing S9b bindings; complete output census rechecked against both production reader inventories.
- Numerics: held integral derivatives remain held and print unresolved when unevaluated; sample outcomes and completion counts distinguish formed, undefined and unresolved comparisons. Algebraic denominator-root and endpoint-convergence screening precedes quadrature; error estimates plus independent precision/panel checks and 18/30-digit result agreement gate 12 printed digits. This is numerical reliability evidence, not a certified bound. The first repair run exited 1 after 66 rows on an invalid NaN predicate, fixed with a regression. I intentionally interrupted repair full-02 after 91 rows to replace the observed general-coefficient root-expansion bottleneck with explicit unresolved coverage; its SIGUSR1 trace, builder-stop.txt and normal guard receipt are retained. Neither was a guard kill; no registry edits. Intentional-stop workload 2909.1531922659997s, observed peak 1257914368 bytes.
- Domains/controls: sequence and record gates expose marker/domain pairs; proved-empty pieces print empty; scoped ConditionSet/ForAll/Exists predicates are sampled (the exact univariate decision procedure supports rational coefficients; other fields remain unresolved), including one-sided attached domains, with explicit undecidable-truth reasons. Controls discard echoed functional residuals, target arithmetic/argument translation/sampling/derivatives/labels/markers/census/missing siblings, and use an independent configuration list. Numeric, derivative and zero/negative-density ablations are detected.
- Tests: `python scripts/s11c_guarded_run.py --pool s9b-comparator --memory-gib 16 --log-directory _scratch/s11c/s9b-comparator/repair-01/tests-07/guard -- python research/pde_ledger_v3/scripts/test_S9b_cross_engine_comparator.py` exited 0; 46 tests, unittest 61.947s; workload 62.48816336299933s; observed cgroup peak 99737600 bytes. Final source/test snapshots retained.
- Full: `python scripts/s11c_guarded_run.py --pool s9b-comparator --memory-gib 16 --log-directory _scratch/s11c/s9b-comparator/repair-01/full-03/guard -- python research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py --output _scratch/s11c/s9b-comparator/repair-01/full-03/comparison.jsonl --accounting _scratch/s11c/s9b-comparator/repair-01/full-03/accounting.jsonl` exited 0; workload 8249.764883671s; launcher 8250.492126230005s; observed cgroup peak 9331417088 bytes. Limits verified; no deadline, refusal or guard kill.
- Accounting/validation: `repair-01/full-03/accounting.jsonl` contains per-row parsed/compared/uncompared leaves and all reasons: 330 joined, 49 PY / 44 WL local, no unjoined; all 753 tag occurrences accounted for. `python _scratch/s11c/s9b-comparator/repair-01/validate_run.py` verified unchanged bindings/input hashes/profiles/points, source snapshots, 118454 samples and the repaired output classes. Sample outcomes: 2,959 residual formed; 1,485 boolean comparison formed; 114,010 unresolved coverage. See `repair-01/full-03/output_summary.json` for exact checks, sample outcomes and unresolved reasons; `not_formed_inventory.json` there lists every remaining gap by tag and path.
- Remaining `not_formed` reasons (counts are emitted occurrences, not distinct rows): `nested sibling has no bound counterpart`: 8722 occurrences; unpaired subfields of A deflection/radar/nonreciprocity/one-way/round-trip and branch/path objects; B gamma, reference matching and implied exchange; C constant/fixed-ratio/power responses, n-dependence and implied exchange; forward mass solution F_MASS_SOLUTION, including the emitted forward/restricted variants. Exact affected tags and field paths: [not_formed_inventory.json](repair-01/full-03/not_formed_inventory.json); `container/scalar representation differs`: 18 occurrences; B_GAMMA_DEFLECTION, B_GAMMA_DIFFERENCE, B_GAMMA_RADAR, F_FLOW_B_GAMMA_DEFLECTION, F_FLOW_B_GAMMA_DIFFERENCE, F_FLOW_B_GAMMA_RADAR, F_LIVE_B_GAMMA_DEFLECTION, F_LIVE_B_GAMMA_DIFFERENCE, F_LIVE_B_GAMMA_RADAR, R_FLOW_ONLY_B_GAMMA_DEFLECTION, R_FLOW_ONLY_B_GAMMA_DIFFERENCE, R_FLOW_ONLY_B_GAMMA_RADAR, R_SPEED_ONLY_B_GAMMA_DEFLECTION, R_SPEED_ONLY_B_GAMMA_DIFFERENCE, R_SPEED_ONLY_B_GAMMA_RADAR, R_TILT_ONLY_B_GAMMA_DEFLECTION, R_TILT_ONLY_B_GAMMA_DIFFERENCE, R_TILT_ONLY_B_GAMMA_RADAR. Spec objects: A branch/traversal/deflection/one-way/radar (shared spec 175–200,228–232), B gamma/reference matching (202–217,233–241), C bulk responses and implied exchange (219–224,242–250); exact affected quantity names/paths are in the linked inventory. Remaining gaps are coverage findings, not a verdict.
- Profiles/points: inverse, shifted, exponential, optical/delta/velocity/tilt-zero, bulk-zero/negative, brane-zero and stiffness-zero; (r,b)=(11/3,7/3),(17/5,11/5),(23/7,13/7), plus parameter probes 0,1,-1,1/2. Exact formulas/assignments/values: `repair-01/full-03/comparison.jsonl`. Admissibility: shared spec 52–57,80,124–149,175–176,219–224; repair item 9; boundary probes remain qualified outside propagation where false. Name-table PY/WL citations below refer to the two engine sources recorded in the comparator; spec citations are to S9b_SHARED_PHYSICS.md unless marked repair.

| SymPy → Wolfram | Spec object / lines | PY line | WL line |
|---|---|---:|---:|
| `s9b_r` → `r` | radial coordinate; 141,175 | 668 | 45 |
| `b` → `b` | flyby impact parameter; 177 | 668 | 45 |
| `Z_E` → `ZE` | emitter endpoint; 178-180 | 668 | 45 |
| `Z_R` → `ZR` | reflector endpoint; 178-180 | 668 | 45 |
| `r_E` → `rE` | emitter radius; 180 | 669 | 69 |
| `r_R` → `rR` | reflector radius; 180 | 669 | 69 |
| `s9b_b_min` → `bFar` | far-zone lower bound; 172-176 | 669 | 45 |
| `GM` → `GM` | independent orbital mass parameter; 202-204 | 670 | 46 |
| `c_0` → `c0` | asymptotic optical speed; 116 | 668 | 45 |
| `K` → `K` | bulk pressure coefficient; 219-220 | 670 | 46 |
| `n` → `n` | bulk pressure exponent; 219-220 | 670 | 46 |
| `m` → `m` | bulk particle mass; 219-220 | 671 | 46 |
| `s` → `s` | optical power exponent; 246 | 670 | 46 |
| `rho_0` → `rho0` | asymptotic bulk number density; 223 | 671 | 46 |
| `Phi` → `Phi` | outward sphere flux; repair directive item 10 | 1041 | 46 |
| `s9b_gamma` → `gamma` | reference-family gamma; 208-217 | 849 | 68 |
| `rho_bulk` → `rho` | local bulk number density; 219,248 | 672 | 92 |
| `s9b_fractional_bulk_change` → `fLocal` | local fractional bulk change; 222 | 1019 | 94 |
| `s9b_delta` → `deltaProfile` | live optical speed shift; 124,141 | 397 | 50 |
| `s9b_V` → `velocityProfile` | live radial material velocity; 60,141,175 | 398 | 50 |
| `s9b_xi` → `xiProfile` | live graph displacement; 83-87,141 | 399 | 51 |
| `rho_br_profile` → `rhoBr` | live brane mass density; 57,69,80 | 400 | 63 |
| `s9b_f` → `fProfile` | live fractional bulk density change; 222 | 401 | 53 |
| `mu_perp` → `muPerp` | live transverse stiffness; 52-57 | 402 | 63 |
| `j_n` → `jn` | live normal mass exchange; 69,80 | 673 | 67 |
| `s9b_mass_V` → `vForward` | forward no-loss velocity; repair item 10 | 1042 | 651 |
| `s9b_velocity_squared` → `yProfile` | velocity-square reduction; 236-240 | 447 | 504 |
| `s9b_velocity_integration_constant` → `velocityIntegrationConstant` | velocity reduction constant; repair item 6 | 470 | 524 |
