# S9b comparator builder report — full guarded comparison completed
Built `research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py` and `test_S9b_cross_engine_comparator.py`. No commits; no residual target was given. All 330 joined rows processed with the unchanged profile/point/stratum scope; findings and unresolved coverage remain for the record, without a verdict.
Measurements (`performance/`): A_DEFLECTION_G001 p0/inverse operands+residual took 72.74996133999957s before and 0.6720223069969506s after compiled quadrature under cProfile (108.25527751466393×). A_ONE_WAY_ER_G121 diagnostic was intentionally interrupted after 105.702 profiled seconds: 96.309s were integral evalf calls induced by printing; after the printer change both operands+residual completed in 26.337695738002368 profiled seconds. Commands/receipts: performance/{baseline-run,optimized-run,g121-run-02,g121-optimized-run}; profiler scripts and statistics retained.
Changes: compiled mpmath tanh-sinh quadrature at 18/30 digits with 15 guard digits; independent integral evaluation and residual arithmetic; bounded caches; per-expression substitution memoization; stored-order printing that does not evaluate integrals; zero-profile derivative/Subs reconstruction fix; object/sample/row timing, cache counters and on-demand stack dumps. No numerical error bound is certified. I intentionally stopped full-03 via systemctl for profiling (1106.293181s telemetry, peak 441126912 bytes); its stale reservation was reconciled by the orchestrator. The later diagnostic SIGINT saved profiling data and received a normal guard receipt; no registry edits by me.
Tests: `python scripts/s11c_guarded_run.py --pool s9b-comparator --memory-gib 16 --log-directory _scratch/s11c/s9b-comparator/tests-08 -- python research/pde_ledger_v3/scripts/test_S9b_cross_engine_comparator.py` exited 0: 32 tests; unittest 17.591s, workload 18.01690012099789s, observed cgroup peak 76664832 bytes. Mutation controls plus quadrature, cache, full-grid, derivative/Subs and non-evaluating-printer regressions. Final source matches tested snapshots.
Full: `python scripts/s11c_guarded_run.py --pool s9b-comparator --memory-gib 16 --log-directory _scratch/s11c/s9b-comparator/full-04/guard -- python research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py --output _scratch/s11c/s9b-comparator/full-04/comparison.jsonl --accounting _scratch/s11c/s9b-comparator/full-04/accounting.jsonl` exited 0; workload 2775.9440095520004s, launcher 2776.960058675999s; observed cgroup peak 8596099072 bytes. Limits verified; no deadline. All telemetry and progress retained.
Accounting: `full-04/accounting.jsonl` gives every row and both parsed/compared/uncompared leaf counts; 330 joined, 49 PY / 44 WL local (engine-local; not joined), no unjoined. `python _scratch/s11c/s9b-comparator/full-04/verify_outputs.py` verified all 753 production tag occurrences against the independent reader inventory, all emitted profile grids, framing, source snapshots and unchanged input hashes/name bindings/profiles/points. Summary: `full-04/output_summary.json`; 7376 sample records; 330 rows retain uncompared leaves. Mathematical outcomes remain uninterpreted.
Machinery: O2 lossless constructor AST/InputForm readers adapted for S9b framing, negative infinity and Alternatives; S9b extraction, separate values/domains/density labels, predicates/markers, symbolic/profile probes, independent census and unsupported-content accounting. Complete-stream validation: 379 PY / 374 WL tags, 9769542 / 1549036 parsed leaves in reader_inventory.json. Final output census matches that independent parse; final test/full-run code is identical. No open operational blocker; unresolved mathematical coverage is printed in the output.
Profiles/points unchanged: inverse-power, shifted rational/log, exponential, optical/delta/velocity/tilt-zero, bulk-zero/negative, brane-zero and stiffness-zero probes; (r,b)=(11/3,7/3),(17/5,11/5),(23/7,13/7), plus parameter probes 0,1,-1,1/2. Exact formulas, endpoint/parameter assignments and per-sample values are in full-04/comparison.jsonl. Admissibility: shared spec 52–57,80,124–149,175–176,219–224; repair item 9. Boundary probes are qualified outside propagation where false. Name-table citations below use the PY/WL engine sources and shared spec, except explicit repair references.

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
