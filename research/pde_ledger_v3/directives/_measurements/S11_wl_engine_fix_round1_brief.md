# Measurements — S11 WL engine fix round 1 brief (generated 2026-08-16 09:17)

Generator: ~/.s11_build/gen_twin_wlfix1.sh (this file is written by that script, never by hand).

## Driver record: every cell's exit, the one memkill, the D4 kill

```
$ grep -E '^CELL_END|^CELL_MEMKILL|^DRIVER_END' /home/trevnorris/.s11_build/wl_percell_driver.log
CELL_END MAIN 2 rc=0 wall=15s out_bytes=67877 2026-08-15 10:40:45
CELL_END MAIN 3 rc=0 wall=15s out_bytes=161094 2026-08-15 10:41:00
CELL_END MAIN 4 rc=0 wall=30s out_bytes=308415 2026-08-15 10:41:30
CELL_END MAIN 5 rc=0 wall=30s out_bytes=636141 2026-08-15 10:42:00
CELL_END XFORM_CURLONLY 2 rc=0 wall=15s out_bytes=681176 2026-08-15 10:42:15
CELL_END XFORM_CURLONLY 3 rc=0 wall=15s out_bytes=744071 2026-08-15 10:42:30
CELL_END XFORM_CURLONLY 4 rc=0 wall=15s out_bytes=838555 2026-08-15 10:42:45
CELL_END XFORM_CURLONLY 5 rc=0 wall=30s out_bytes=999692 2026-08-15 10:43:15
CELL_END XFORM_EXTRA 2 rc=0 wall=31s out_bytes=1169320 2026-08-15 10:43:46
CELL_END XFORM_EXTRA 3 rc=0 wall=15s out_bytes=1267566 2026-08-15 10:44:01
CELL_END XFORM_EXTRA 4 rc=0 wall=15s out_bytes=1423053 2026-08-15 10:44:16
CELL_END XFORM_EXTRA 5 rc=0 wall=30s out_bytes=1755839 2026-08-15 10:44:46
CELL_END XFORM_DIVONLY 3 rc=0 wall=15s out_bytes=1819141 2026-08-15 10:45:01
CELL_END XFORM_DIVONLY 4 rc=0 wall=15s out_bytes=1912662 2026-08-15 10:45:16
CELL_END XFORM_TRACELESS 3 rc=0 wall=15s out_bytes=2016959 2026-08-15 10:45:31
CELL_END XFORM_TRACELESS 4 rc=0 wall=30s out_bytes=2184754 2026-08-15 10:46:01
CELL_END XCOEF_BSCALE 3 rc=0 wall=15s out_bytes=2284530 2026-08-15 10:46:16
CELL_END XCOEF_BSIGN 3 rc=0 wall=15s out_bytes=2358355 2026-08-15 10:46:31
CELL_MEMKILL XKIN_ANISO 2 available=891MB killing wrapper=2248514 kernel=2248549 2026-08-15 11:07:46
CELL_END XKIN_ANISO 2 rc=137 wall=1275s out_bytes=2774060 2026-08-15 11:07:46
CELL_END XKIN_ANISO 3 rc=0 wall=285s out_bytes=3594396 2026-08-15 11:12:31
CELL_END XKIN_ANISO 4 rc=137 wall=79127s out_bytes=3921488 2026-08-16 09:11:18
DRIVER_END 2026-08-16 09:11:18
```

## The hang site: lines 205, 211, 216 of the engine

```
$ sed -n '205p;211p;216p' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl
  emitCell[package, dimension, baseQuantity <> "_IDENTICALLY_SATISFIED", {
  unrestrictedReduction = Quiet[Reduce[And @@ equations, variables, Complexes]];
  emitCell[package, dimension, baseQuantity <> "_INCONSISTENT", {
```

## Last three tags emitted by XKIN D4 (committed .out at e2928b49)

```
$ grep -oE '^WL_S11_XKIN_ANISO_D4_[A-Z0-9_]+' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S11_stray_longitudinal_mathematica_audit.out | tail -3
WL_S11_XKIN_ANISO_D4_ROOT2_STACKED_DROP_K_EQUATIONS
WL_S11_XKIN_ANISO_D4_ROOT2_STACKED_DROP_K_SOLUTION
WL_S11_XKIN_ANISO_D4_ROOT2_STACKED_DROP_K_IDENTICALLY_SATISFIED
```

## The .out went silent at 16:58 on 2026-08-15; the cell was killed 09:11 on 2026-08-16

```
$ stat -c 'mtime %y' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S11_stray_longitudinal_mathematica_audit.out; grep 'CELL_END XKIN_ANISO 4' /home/trevnorris/.s11_build/wl_percell_driver.log
mtime 2026-08-15 16:58:40.424295301 -0600
CELL_END XKIN_ANISO 4 rc=137 wall=79127s out_bytes=3921488 2026-08-16 09:11:18
```

## Call-site recurrence of the quantifier-elimination class in the engine

```
$ grep -n 'Reduce\[' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl | head -20; echo ---; grep -n 'Resolve\[' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl | head -20
211:  unrestrictedReduction = Quiet[Reduce[And @@ equations, variables, Complexes]];
398:  v1Basis = If[v1Raw === {}, {}, RowReduce[v1Raw]];
414:  v2Basis = If[v2Raw === {}, {}, RowReduce[v2Raw]];
452:    v6Basis = If[v6Ambient === {}, {}, RowReduce[v6Ambient]]
---
137:      Resolve[Exists[quantifiedVariables, predicate], Reals]]]]
679:  attemptOutcome = Quiet[Resolve[testObject, Reals]];
```

## How many times the class ran in the completed D3 cell (record census)

```
$ grep -cE '^WL_S11_XKIN_ANISO_D3_[A-Z0-9_]*(INCONSISTENT|REAL_STATUS):' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S11_stray_longitudinal_mathematica_audit.out
328
```

## Spec: the object definitions and the pinned payload form (S11_SHARED_PHYSICS.md)

```
$ sed -n '247p;250,251p;253,260p;285,288p;290,291p' /var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
| `_INCONSISTENT` | the typed symbolic test of whether the system admits **no** solution at all |
⭐ `_IDENTICALLY_SATISFIED` and `_INCONSISTENT` are what separate the three degenerate cases, and they are
computed from `_EQUATIONS`, ⛔ never read off the solver's empty token.
⭐⭐ **THE BOOLEAN-TEST PAYLOAD FORM IS PINNED.** Each `_IDENTICALLY_SATISFIED` and `_INCONSISTENT`
payload has this ordered field sequence:

```
STATUS_TOKEN: <exactly one of PROVED_TRUE · PROVED_FALSE · UNDECIDED>
TEST_OBJECT:  <the live CAS boolean-valued object, in the §8 re-parseable CAS form>
OPERANDS:     <the live operands supplied to the test>
```
| `_CANONICAL_LOCUS` | when every residual `lhs − rhs` in `_EQUATIONS` is a polynomial expression in all symbols it contains, the reduced Gröbner basis of the ideal they generate in the named solve variables, for lexicographic order in the emitted solve-variable ordering and over the exact rational-function coefficient field generated by the remaining symbols; otherwise the single token `NOT_APPLICABLE` |
| `_REAL_STATUS` | exactly one of `PROVED_EMPTY` · `PROVED_NONEMPTY` · `UNDECIDED`, for the locus over the reals under the full joint premise set in force |
| `_REAL_WITNESS` | for `PROVED_NONEMPTY`, an exact point satisfying `_EQUATIONS` and the full joint premise set — ⭐ **an assignment to exactly the solve variables this locus names**, ⛔ and to nothing else; for either other status, the single token `NOT_APPLICABLE` |
| `_REAL_STATUS_OPERANDS` | the live `_EQUATIONS`, named solve variables and full joint premise set on which `_REAL_STATUS` was computed |
⚠ `_CANONICAL_LOCUS` names an object, ⛔ not an admissibility algorithm. No real-domain method is
prescribed. ⛔ Do not replace it with a complex-locus claim, and ⛔ do not turn a component returned by
```

## Diagnosis scratches archived out of /tmp before any builder launch

```
$ ls -la /home/trevnorris/.s11_build/d4_stall_diagnosis_scratches.tar.gz; tar tzf /home/trevnorris/.s11_build/d4_stall_diagnosis_scratches.tar.gz | head -8
-rw-rw-r-- 1 trevnorris trevnorris 76040 Aug 16 09:12 /home/trevnorris/.s11_build/d4_stall_diagnosis_scratches.tar.gz
s11_d4_reduce_diag/
s11_d4_reduce_diag/00_identify_stage.stderr
s11_d4_reduce_diag/03_measure_emitted_equations.stderr
s11_d4_reduce_diag/01_reconstruct_operands.py
s11_d4_reduce_diag/04_stacked_minors_laplace.stdout
s11_d4_reduce_diag/07_successor_queue.stderr
s11_d4_reduce_diag/07_successor_queue.py
s11_d4_reduce_diag/04_stacked_minors_laplace.stderr
```

## Diagnosis GB timings, both legs (content of bases withheld as answer-bearing)

```
$ tar xzf /home/trevnorris/.s11_build/d4_stall_diagnosis_scratches.tar.gz -O d4diag/s3_gb_probe.out 2>/dev/null | grep -E 'GROEBNER (start|done)|system:' | sed 's/: basis size=[0-9]*//' ; tar xzf /home/trevnorris/.s11_build/d4_stall_diagnosis_scratches.tar.gz -O d4diag/s5_successor_probe.out 2>/dev/null | grep -iE 'GROEBNER done|exit' | sed 's/: basis size=[0-9]*//' ; tar xzf /home/trevnorris/.s11_build/d4_stall_diagnosis_scratches.tar.gz -O s11_d4_reduce_diag/06_clean_elim.stdout 2>/dev/null | grep -E 'GROEBNER_SECONDS'
[     1.6s] D=3 system: 3 minor equations + radical relation; gens=[k1, k2, k3, r], coefficient field=QQ(bComp,muR,sRho)
[     1.6s] D=3 GROEBNER start (grevlex)
[     2.6s] D=3 GROEBNER done in 1.0s
[    24.2s] D=4 system: 4 minor equations + radical relation; gens=[k1, k2, k3, k4, r], coefficient field=QQ(bComp,muR,sRho)
[    24.2s] D=4 GROEBNER start (grevlex)
[    34.1s] D=4 GROEBNER done in 9.8s
exit=124
N3_GROEBNER_SECONDS 0.056144
N4_GROEBNER_SECONDS 0.692928
```

## Regression-bar operand: the 3 control-overlapping cells were byte-identical (per-cell rerun vs single-kernel control, physics tags)

```
$ CTL=/home/trevnorris/.s11_build/wl_sweep1_partial_singlekernel.out; OUTF=/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S11_stray_longitudinal_mathematica_audit.out; for cell in WL_S11_MAIN_D2 WL_S11_XFORM_CURLONLY_D2 WL_S11_XFORM_EXTRA_D2; do cmp -s <(grep "^${cell}_" $CTL | grep -vE '_TAG_NAMES' | sort) <(grep "^${cell}_" $OUTF | grep -vE '_TAG_NAMES' | sort) && echo "$cell IDENTICAL" || echo "$cell DIFFERS"; done
WL_S11_MAIN_D2 IDENTICAL
WL_S11_XFORM_CURLONLY_D2 IDENTICAL
WL_S11_XFORM_EXTRA_D2 IDENTICAL
```

## Head commit the brief targets

```
$ cd /var/projects/toy_physics && git log --oneline -1
e2928b49 WL per-cell sweep DRIVER_END: 19/21 cells complete, XKIN D2+D4 partial — the honest record
```

