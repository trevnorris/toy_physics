# Measurements — S11 WL engine fix round 2 brief (generated 2026-08-16 19:33)

Generator: ~/.s11_build/gen_twin_wlfix2.sh (this file is written by that script, never by hand).

## Guard record: every XKIN death and the D3 completion (floor kill, rc, wall, emissions, last tag)

```
$ grep -E '^CELL_(START|MEMGUARD|END)' /home/trevnorris/.s11_build/fix1_build/guarded_cells_record.log
CELL_START XKIN_ANISO_D3 command=/usr/bin/wolframscript -file /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl XKIN_ANISO 3 2026-08-16T10:03:04-06:00
CELL_END XKIN_ANISO_D3 rc=0 logger_rc=0 wall=255s guard=NONE completion=COMPLETE emissions=2501 max_inter_emission_gap=36.183270s last_tag=WL_S11_LOCAL_TAG_NAMES 2026-08-16T10:07:19-06:00
CELL_START XKIN_ANISO_D4 command=/usr/bin/wolframscript -file /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl XKIN_ANISO 4 2026-08-16T10:07:28-06:00
CELL_START XKIN_ANISO_D4 command=/usr/bin/wolframscript -file /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl XKIN_ANISO 4 2026-08-16T10:48:45-06:00
CELL_MEMGUARD XKIN_ANISO_D4 available_kb=905176 floor_kb=1048576 wrapper=162225 kernel=162586 2026-08-16T11:18:37-06:00
CELL_END XKIN_ANISO_D4 rc=137 logger_rc=0 wall=1798s guard=MEMORY completion=INCOMPLETE emissions=999 max_inter_emission_gap=137.157751s last_tag=WL_S11_XKIN_ANISO_D4_STRATUM3_ROOT2_N2_RANK_CHANGE_LOCUS_EQUATIONS 2026-08-16T11:18:43-06:00
CELL_START XKIN_ANISO_D2 command=/usr/bin/wolframscript -file /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl XKIN_ANISO 2 2026-08-16T11:20:31-06:00
CELL_MEMGUARD XKIN_ANISO_D2 available_kb=871348 floor_kb=1048576 wrapper=377009 kernel=377479 2026-08-16T11:42:58-06:00
CELL_END XKIN_ANISO_D2 rc=137 logger_rc=0 wall=1352s guard=MEMORY completion=INCOMPLETE emissions=1411 max_inter_emission_gap=331.620483s last_tag=WL_S11_XKIN_ANISO_D2_STRATUM6_ROOT1_N2_NULLITY_CHANGE_LOCUS_REAL_STATUS_OPERANDS 2026-08-16T11:43:03-06:00
CELL_START XKIN_ANISO_D4 command=/usr/bin/wolframscript -file /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl XKIN_ANISO 4 2026-08-16T16:27:05-06:00
CELL_START XKIN_ANISO_D4 command=/usr/bin/wolframscript -file /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl XKIN_ANISO 4 2026-08-16T16:34:03-06:00
CELL_MEMGUARD XKIN_ANISO_D4 available_kb=901140 floor_kb=1048576 wrapper=1002341 kernel=1002611 2026-08-16T16:55:25-06:00
CELL_END XKIN_ANISO_D4 rc=137 logger_rc=0 wall=1287s guard=MEMORY completion=INCOMPLETE emissions=999 max_inter_emission_gap=137.159875s last_tag=WL_S11_XKIN_ANISO_D4_STRATUM3_ROOT2_N2_RANK_CHANGE_LOCUS_EQUATIONS 2026-08-16T16:55:30-06:00
CELL_START XKIN_ANISO_D2 command=/usr/bin/wolframscript -file /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl XKIN_ANISO 2 2026-08-16T16:55:41-06:00
CELL_MEMGUARD XKIN_ANISO_D2 available_kb=876592 floor_kb=1048576 wrapper=1173416 kernel=1173693 2026-08-16T17:17:31-06:00
CELL_END XKIN_ANISO_D2 rc=137 logger_rc=0 wall=1315s guard=MEMORY completion=INCOMPLETE emissions=1411 max_inter_emission_gap=327.642844s last_tag=WL_S11_XKIN_ANISO_D2_STRATUM6_ROOT1_N2_NULLITY_CHANGE_LOCUS_REAL_STATUS_OPERANDS 2026-08-16T17:17:36-06:00
```

## D4 determinism: emission count and last tag across all three death artifacts

```
$ for f in /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D4.out /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D4.out.prerepair2 /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D4.out.killed_run1; do echo -n "$(basename $f): "; echo -n "$(grep -cE '^WL_S11_' $f) emissions, last="; grep -oE '^WL_S11_[A-Z0-9_]+' $f | tail -1; done
guarded_XKIN_ANISO_D4.out: 999 emissions, last=WL_S11_XKIN_ANISO_D4_STRATUM3_ROOT2_N2_RANK_CHANGE_LOCUS_EQUATIONS
guarded_XKIN_ANISO_D4.out.prerepair2: 999 emissions, last=WL_S11_XKIN_ANISO_D4_STRATUM3_ROOT2_N2_RANK_CHANGE_LOCUS_EQUATIONS
guarded_XKIN_ANISO_D4.out.killed_run1: 999 emissions, last=WL_S11_XKIN_ANISO_D4_STRATUM3_ROOT2_N2_RANK_CHANGE_LOCUS_EQUATIONS
```

## D2 determinism: emission count and last tag across both death artifacts

```
$ for f in /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D2.out /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D2.out.prerepair2; do echo -n "$(basename $f): "; echo -n "$(grep -cE '^WL_S11_' $f) emissions, last="; grep -oE '^WL_S11_[A-Z0-9_]+' $f | tail -1; done
guarded_XKIN_ANISO_D2.out: 1411 emissions, last=WL_S11_XKIN_ANISO_D2_STRATUM6_ROOT1_N2_NULLITY_CHANGE_LOCUS_REAL_STATUS_OPERANDS
guarded_XKIN_ANISO_D2.out.prerepair2: 1411 emissions, last=WL_S11_XKIN_ANISO_D2_STRATUM6_ROOT1_N2_NULLITY_CHANGE_LOCUS_REAL_STATUS_OPERANDS
```

## Guard-kill attribution: MEMGUARD record count per cell (D4's third death artifact was an external harness kill at the same frontier, not a guard kill)

```
$ for d in 4 2; do echo -n "D$d CELL_MEMGUARD records: "; grep -c "CELL_MEMGUARD XKIN_ANISO_D$d" /home/trevnorris/.s11_build/fix1_build/guarded_cells_record.log; done; echo -n 'D4 CELL_START records: '; grep -c 'CELL_START XKIN_ANISO_D4' /home/trevnorris/.s11_build/fix1_build/guarded_cells_record.log
D4 CELL_MEMGUARD records: 2
D2 CELL_MEMGUARD records: 2
D4 CELL_START records: 4
```

## Terminal silence: seconds from last emission to the guard's floor crossing (latest runs, sub-second from the tsv epoch)

```
$ for d in 4 2; do le=$(tail -1 /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D$d.emission_times.tsv | cut -f1); gt=$(grep "CELL_MEMGUARD XKIN_ANISO_D$d" /home/trevnorris/.s11_build/fix1_build/guarded_cells_record.log | tail -1 | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9:]{8}-[0-9:]{2}:[0-9]{2}'); ge=$(date -d "$gt" +%s); awk -v le="$le" -v ge="$ge" -v d="$d" 'BEGIN{printf "D%s: last_emit=%.3f guard=%d silence=%.1fs\n", d, le, ge, ge-le}'; done
D4: last_emit=1786920731.867 guard=1786920925 silence=193.1s
D2: last_emit=1786922228.687 guard=1786922251 silence=22.3s
```

## D2 accumulation span: first STRATUM5 emission to the guard kill

```
$ s5=$(awk -F'\t' '$3 ~ /STRATUM5/ {print $1; exit}' /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D2.emission_times.tsv); gt=$(grep 'CELL_MEMGUARD XKIN_ANISO_D2' /home/trevnorris/.s11_build/fix1_build/guarded_cells_record.log | tail -1 | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9:]{8}-[0-9:]{2}:[0-9]{2}'); ge=$(date -d "$gt" +%s); awk -v a="$s5" -v b="$ge" 'BEGIN{printf "stratum5_first_emit_to_guard=%.1fs\n", b-a}'
stratum5_first_emit_to_guard=1222.7s
```

## D4 death site source: the _EQUATIONS emit (wl:311) and the unbounded Solve (wl:312)

```
$ sed -n '305,316p' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl
   unrestrictedVariables, unrestrictedTestObject, unrestrictedExpression,
   unrestrictedDecision, unrestrictedReduction, inconsistentOperands,
   inconsistentTest, realAdmissible, admittedBranches, branch,
   branchDecision, branchOperands, branchTest,
   branchStatus, point, canonical, fullCondition, realTest, realStatus,
   realDecision, realWitness, realStatusOperands},
  emitCell[package, dimension, baseQuantity <> "_EQUATIONS", equations];
  solution = Quiet[Solve[And @@ equations, variables]];
  emitCell[package, dimension, baseQuantity <> "_SOLUTION", {
    "SOLVE_VARIABLES" -> variables,
    "SOLUTION_SET" -> solution
  }];
```

## D2 death site source: the stacked join and assumedRank call (wl:1076-1077)

```
$ sed -n '1070,1080p' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl
    If[TrueQ[componentCounts], AppendTo[rootCountStatuses, countStatus]];
    nullity = dimension - rank;
    countStatus = emitCountObject[package, dimension,
      rootQuantity[scope, pointEvidence, rootIndex, "N2_NULLITY"], nullity,
      componentCounts, freeParameters, rankDecision];
    If[TrueQ[componentCounts], AppendTo[rootCountStatuses, countStatus]];
    stacked = Join[matrixAtRoot, {wavevector}];
    stackedRank = assumedRank[stacked, assumptions];
    stackedDecision = If[TrueQ[componentCounts],
      matrixRankDecision[stacked, stackedRank, freeParameters, assumptions],
      tokenNotApplicable];
```

## The unbounded helpers: assumedRank / zeroTest / engineSimplify (wl:72-78)

```
$ sed -n '72,78p' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl
engineSimplify[expression_, assumptions_] :=
  FullSimplify[expression, Assumptions -> assumptions];
unrestrictedSimplify[expression_] := FullSimplify[expression];
zeroTest[assumptions_] := Function[value,
  TrueQ[engineSimplify[value == 0, assumptions]]];
assumedRank[matrix_, assumptions_] :=
  MatrixRank[matrix, ZeroTest -> zeroTest[assumptions]];
```

## D2 inter-emission gaps over 60 s (the unbounded FullSimplify-chain calls, incl. the 327.6 s completed twin of the killer)

```
$ awk -F'\t' 'NR>1{gap=$1-prev; if(gap>60) printf "%.1fs before %s\n", gap, $3} {prev=$1}' /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D2.emission_times.tsv
327.6s before WL_S11_XKIN_ANISO_D2_STRATUM5_ROOT1_N3_STACKED_RANK
202.8s before WL_S11_XKIN_ANISO_D2_STRATUM5_ROOT1_N5_M_DOT_K
109.6s before WL_S11_XKIN_ANISO_D2_STRATUM5_ROOT1_N6_DOT_K
244.2s before WL_S11_XKIN_ANISO_D2_STRATUM5_ROOT1_N6_RESIDUAL
```

## D2 stratum spans (accumulation: the last two reached strata dominate the run)

```
$ awk -F'\t' '{if (match($3,/STRATUM[0-9]+/)) {s=substr($3,RSTART,RLENGTH); if(!(s in first)) first[s]=$1; last[s]=$1}} END {for (s in first) printf "%s span %.0fs\n", s, last[s]-first[s]}' /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D2.emission_times.tsv | sort -V
STRATUM1 span 0s
STRATUM2 span 0s
STRATUM3 span 1s
STRATUM4 span 1s
STRATUM5 span 1034s
STRATUM6 span 143s
```

## D4 gaps over 60 s (the generic-pass unbounded Solves that completed; nothing after the death tag)

```
$ awk -F'\t' 'NR>1{gap=$1-prev; if(gap>60) printf "%.1fs before %s\n", gap, $3} {prev=$1}' /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D4.emission_times.tsv
136.6s before WL_S11_XKIN_ANISO_D4_ROOT2_RANK_DROP_JOINT_SOLUTION
137.2s before WL_S11_XKIN_ANISO_D4_ROOT3_RANK_DROP_JOINT_SOLUTION
```

## D4 death-site Solve operand size (from its own emitted _EQUATIONS payload: equation count and bytes) and the stratum's free-parameter count

```
$ line=$(grep 'STRATUM3_ROOT2_N2_RANK_CHANGE_LOCUS_EQUATIONS' /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D4.out | tail -1); echo "bytes=${#line}"; echo -n 'equation relations (==): '; echo "$line" | grep -oE '==' | wc -l; fp=$(grep 'WL_S11_XKIN_ANISO_D4_STRATUM3_FREE_PARAMETERS:' /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D4.out | tail -1); echo -n 'free parameters (unknowns): '; echo "$fp" | sed 's/^[^{]*{//; s/}.*//' | tr ',' '\n' | wc -l
bytes=1163
equation relations (==): 16
free parameters (unknowns): 7
```

## D2 sign-twin: STRATUM5 vs STRATUM6 defining equations (committed partial payloads; the two strata differ only in the radical's sign)

```
$ grep -E 'WL_S11_XKIN_ANISO_D2_STRATUM[56]_DEFINING_EQUATIONS:' /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D2.out
WL_S11_XKIN_ANISO_D2_STRATUM5_DEFINING_EQUATIONS: {k1 == -Sqrt[-((bComp*k2^2*muR)/(bComp^2 - 2*bComp*muR*sRho + muR^2*sRho^2)) - (bComp^2*k2^2*sRho)/(bComp^2 - 2*bComp*muR*sRho + muR^2*sRho^2) + (4*bComp*k2^2*muR*sRho)/(bComp^2 - 2*bComp*muR*sRho + muR^2*sRho^2) - (k2^2*muR^2*sRho)/(bComp^2 - 2*bComp*muR*sRho + muR^2*sRho^2) - (bComp*k2^2*muR*sRho^2)/(bComp^2 - 2*bComp*muR*sRho + muR^2*sRho^2) + Sqrt[(2*bComp*k2^2*muR + 2*bComp^2*k2^2*sRho - 8*bComp*k2^2*muR*sRho + 2*k2^2*muR^2*sRho + 2*bComp*k2^2*muR*sRho^2)^2 - 4*(k2^4*muR^2 - 2*bComp*k2^4*muR*sRho + bComp^2*k2^4*sRho^2)*(bComp^2 - 2*bComp*muR*sRho + muR^2*sRho^2)]/(2*(bComp^2 - 2*bComp*muR*sRho + muR^2*sRho^2))]}
WL_S11_XKIN_ANISO_D2_STRATUM6_DEFINING_EQUATIONS: {k1 == Sqrt[-((bComp*k2^2*muR)/(bComp^2 - 2*bComp*muR*sRho + muR^2*sRho^2)) - (bComp^2*k2^2*sRho)/(bComp^2 - 2*bComp*muR*sRho + muR^2*sRho^2) + (4*bComp*k2^2*muR*sRho)/(bComp^2 - 2*bComp*muR*sRho + muR^2*sRho^2) - (k2^2*muR^2*sRho)/(bComp^2 - 2*bComp*muR*sRho + muR^2*sRho^2) - (bComp*k2^2*muR*sRho^2)/(bComp^2 - 2*bComp*muR*sRho + muR^2*sRho^2) + Sqrt[(2*bComp*k2^2*muR + 2*bComp^2*k2^2*sRho - 8*bComp*k2^2*muR*sRho + 2*k2^2*muR^2*sRho + 2*bComp*k2^2*muR*sRho^2)^2 - 4*(k2^4*muR^2 - 2*bComp*k2^4*muR*sRho + bComp^2*k2^4*sRho^2)*(bComp^2 - 2*bComp*muR*sRho + muR^2*sRho^2)]/(2*(bComp^2 - 2*bComp*muR*sRho + muR^2*sRho^2))]}
```

## D2 death-site N1 matrix payload size and radical content (Sqrt/Abs/I node counts)

```
$ line=$(grep 'STRATUM6_ROOT1_N1:' /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D2.out | tail -1); echo "bytes=${#line}"; echo -n 'Sqrt nodes: '; echo "$line" | grep -oE 'Sqrt' | wc -l; echo -n 'Abs nodes: '; echo "$line" | grep -oE 'Abs' | wc -l; echo -n 'Complex I present: '; echo "$line" | grep -cE 'Complex|\bI\b'
bytes=871
Sqrt nodes: 9
Abs nodes: 12
Complex I present: 1
```

## D3 control N1 payloads for the same record family (bytes; radical content)

```
$ grep -E 'WL_S11_XKIN_ANISO_D3_STRATUM[0-9]+_ROOT[0-9]+_N1:' /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D3.out | awk '{print "bytes=" length($0)}' | sort -n | tail -3; echo -n 'Sqrt nodes across ALL D3 N1 payloads: '; grep -E 'WL_S11_XKIN_ANISO_D3_STRATUM[0-9]+_ROOT[0-9]+_N1:' /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D3.out | grep -oE 'Sqrt' | wc -l
bytes=199
bytes=227
bytes=235
Sqrt nodes across ALL D3 N1 payloads: 0
```

## Feasibility timings, analyst A (grok): staged bounded-route probes on reconstructed operands (names+durations only; result-stating name fragments redacted)

```
$ tar xzf /home/trevnorris/.s11_build/wall2_diag_scratches_grok.tar.gz -O s11_wall2_obBN/04_bounded_route.stdout | grep -E '^.progress. END' | sed 's/^.progress. //; s/confirm_radical_branch_implies_disc_zero/radical_branch_disc_membership_probe/; s/_solve_each_principal_factor/_per_factor_solve/'; tar xzf /home/trevnorris/.s11_build/wall2_diag_scratches_grok.tar.gz -O s11_wall2_obBN/05_d2_rank_mod_disc.stdout | grep -E '^dt_|^ELAPSED'
END D4_polynomialize_residuals status=OK dt_s=0.004024
END D4_factor_numerators status=OK dt_s=0.159452
END D4_per_factor_solve status=OK dt_s=0.127013
END D4_groebner_polynomialized status=OK dt_s=0.101449
END D4_unreduced_sympy_solve_manual status=OK dt_s=0.158657
END D2_discriminant_of_detM_in_omegaSquared status=OK dt_s=0.063487
END D2_radical_branch_disc_membership_probe status=OK dt_s=1.028076
END D2_generic_DomainMatrix_rank_of_M_and_stacked status=OK dt_s=0.008954
END D2_DomainMatrix_rank_on_polynomial_coincidence_slice status=OK dt_s=0.147040
END D2_abs_form_stacked_minors_size_only_no_simplify status=OK dt_s=0.112647
END D2_together_of_one_Abs_minor_no_FullSimplify status=OK dt_s=0.163269
END D3_discriminant_of_detM status=OK dt_s=2.944122
END D4_discriminant_of_detM status=OK dt_s=57.525401
dt_setup 0.07347346842288971
dt_sub 0.025721787475049496
dt_reduce 0.13890156894922256
dt_k1sq 0.5155844651162624
ELAPSED 0.759020384401083
```

## Feasibility timings, analyst B (fresh agent): exact-route durations (solution content stripped as answer-bearing)

```
$ tar xzf /home/trevnorris/.s11_build/wall2_diag_scratches_agent.tar.gz -O s11_wall2_AQ1c/08_fixed_parser_rerun.out | grep -oE 'factor-lattice time [0-9.]+s|completed in [0-9.]+s'
factor-lattice time 0.064s
completed in 0.207s
factor-lattice time 0.117s
completed in 0.162s
```

## No cross-engine oracle: the SymPy record's STRATUM_ORDERING payload on every cell

```
$ grep -E 'STRATUM_ORDERING' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11_stray_longitudinal_sympy_audit.out | grep -vE 'CANDIDATE|LOCAL_TAG' | cut -c1-60 | sort | uniq -c
      1 PY_S11_MAIN_D2_STRATUM_ORDERING: Tuple()
      1 PY_S11_MAIN_D3_STRATUM_ORDERING: Tuple()
      1 PY_S11_MAIN_D4_STRATUM_ORDERING: Tuple()
      1 PY_S11_MAIN_D5_STRATUM_ORDERING: Tuple()
      1 PY_S11_XCOEF_BSCALE_D3_STRATUM_ORDERING: Tuple()
      1 PY_S11_XCOEF_BSIGN_D3_STRATUM_ORDERING: Tuple()
      1 PY_S11_XFORM_CURLONLY_D2_STRATUM_ORDERING: Tuple()
      1 PY_S11_XFORM_CURLONLY_D3_STRATUM_ORDERING: Tuple()
      1 PY_S11_XFORM_CURLONLY_D4_STRATUM_ORDERING: Tuple()
      1 PY_S11_XFORM_CURLONLY_D5_STRATUM_ORDERING: Tuple()
      1 PY_S11_XFORM_DIVONLY_D3_STRATUM_ORDERING: Tuple()
      1 PY_S11_XFORM_DIVONLY_D4_STRATUM_ORDERING: Tuple()
      1 PY_S11_XFORM_EXTRA_D2_STRATUM_ORDERING: Tuple()
      1 PY_S11_XFORM_EXTRA_D3_STRATUM_ORDERING: Tuple()
      1 PY_S11_XFORM_EXTRA_D4_STRATUM_ORDERING: Tuple()
      1 PY_S11_XFORM_EXTRA_D5_STRATUM_ORDERING: Tuple()
      1 PY_S11_XFORM_TRACELESS_D3_STRATUM_ORDERING: Tuple()
      1 PY_S11_XFORM_TRACELESS_D4_STRATUM_ORDERING: Tuple()
      1 PY_S11_XKIN_ANISO_D2_STRATUM_ORDERING: Tuple()
      1 PY_S11_XKIN_ANISO_D3_STRATUM_ORDERING: Tuple()
      1 PY_S11_XKIN_ANISO_D4_STRATUM_ORDERING: Tuple()
```

## Guard contract the builder must use (floor, wall, per-PID kill) and its pinned hash

```
$ grep -nE 'memory_floor_kb|wall_limit_seconds|find_kernel_descendant' /home/trevnorris/.s11_build/fix1_build/run_guarded_cell.sh | head -8; sha256sum /home/trevnorris/.s11_build/fix1_build/run_guarded_cell.sh
14:memory_floor_kb=1048576
15:wall_limit_seconds=14400
17:find_kernel_descendant() {
76:  discovered_kernel=$(find_kernel_descendant || true)
80:  if [[ ${available_kb} -lt ${memory_floor_kb} ]]; then
84:      "${stem}" "${available_kb}" "${memory_floor_kb}" "${wrapper_pid}" "${kernel_pid}" "${guard_iso}" | tee -a "${record}"
96:  if [[ ${elapsed} -ge ${wall_limit_seconds} ]]; then
100:      "${stem}" "${elapsed}" "${wall_limit_seconds}" "${wrapper_pid}" "${kernel_pid}" "${guard_iso}" | tee -a "${record}"
ba17f9ab9016e40fa9932cfa875cbd2b40b5c7d731ae5257360cbdbfb1d0a664  /home/trevnorris/.s11_build/fix1_build/run_guarded_cell.sh
```

## The two registered out-of-scope defects exist in the register (entry heading lines)

```
$ grep -nE '^[0-9]+\. \*\*.*(IDENTICALLY_SATISFIED|ROOT_COINCIDENCE)' /var/projects/toy_physics/research/pde_ledger_v3/DEFECT_REGISTER.md
806:1. **`_IDENTICALLY_SATISFIED` is computed pointwise, not identically** (`wl:203`:
816:2. **`ROOT_COINCIDENCE_*` loci are one joint system over ALL root pairs, not per-pair loci**
```

## Diagnosis scratches archived out of /tmp before any builder launch

```
$ ls -la /home/trevnorris/.s11_build/wall2_diag_scratches_grok.tar.gz /home/trevnorris/.s11_build/wall2_diag_scratches_agent.tar.gz; tar tzf /home/trevnorris/.s11_build/wall2_diag_scratches_grok.tar.gz | head -6; tar tzf /home/trevnorris/.s11_build/wall2_diag_scratches_agent.tar.gz | head -6
-rw-rw-r-- 1 trevnorris trevnorris 22661 Aug 16 19:02 /home/trevnorris/.s11_build/wall2_diag_scratches_agent.tar.gz
-rw-rw-r-- 1 trevnorris trevnorris 44259 Aug 16 19:02 /home/trevnorris/.s11_build/wall2_diag_scratches_grok.tar.gz
s11_wall2_obBN/
s11_wall2_obBN/04_bounded_route.json
s11_wall2_obBN/02_timeline.py
s11_wall2_obBN/WL_S11_XKIN_ANISO_D4_STRATUM3_ROOT1_VALUE.wl
s11_wall2_obBN/WL_S11_XKIN_ANISO_D2_STRATUM4_ROOT2_N2_RANK_CHANGE_LOCUS_EQUATIONS.wl
s11_wall2_obBN/WL_S11_XKIN_ANISO_D4_STRATUM3_FREE_PARAMETERS.wl
s11_wall2_AQ1c/
s11_wall2_AQ1c/07_interval_certificates.out
s11_wall2_AQ1c/06_exact_certificates.py
s11_wall2_AQ1c/05_zero_discrimination.out
s11_wall2_AQ1c/02_operands.py
s11_wall2_AQ1c/04_bounded_routes.py
```

## HEAD at twin-generation time (the brief's commit is this commit's child)

```
$ cd /var/projects/toy_physics && git log --oneline -1
51f4018f Round-2 brief review prompt: defective-repair construction mandatory, both legs identical
```

