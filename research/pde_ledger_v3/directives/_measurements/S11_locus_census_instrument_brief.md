# Measurements — S11 locus-census instrument brief (generated 2026-08-17T07:48:18-06:00 by gen_twin_census_instr.sh)

## Record identities
```text
$ git log -1 --oneline -- research/pde_ledger_v3/mathematica/out/S11_stray_longitudinal_mathematica_audit.out
a4cf6539 WL fix round 2 CLOSED: canonical record 21/21 complete, both legs reported, four defects registered
$ git log -1 --oneline -- research/pde_ledger_v3/scripts/out/S11_stray_longitudinal_sympy_audit.out
19591194 Sweep-5 complete: 21/21 cells, zero skips, 2h41m — the full S11 SymPy record
$ wc -l research/pde_ledger_v3/mathematica/out/S11_stray_longitudinal_mathematica_audit.out research/pde_ledger_v3/scripts/out/S11_stray_longitudinal_sympy_audit.out
    16587 research/pde_ledger_v3/mathematica/out/S11_stray_longitudinal_mathematica_audit.out
     5821 research/pde_ledger_v3/scripts/out/S11_stray_longitudinal_sympy_audit.out
    22408 total
```

## Old containment harness: parse failure on real records
```text
$ python3 ~/.s11_build/fix2_build/sympy_completeness.py /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D2.out /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D4.out 2>&1 | tail -3; echo rc=$?
SYMPY_LOCUS index=361 tag=WL_S11_XKIN_ANISO_D4_STRATUM9_ROOT_COUNT_DISTINCT_CHANGE_LOCUS_SOLUTION outcome=PARSE_ERROR detail=ValueError:not a braced object: "NOT_APPLICABLE"
SYMPY_COMPLETENESS_SUMMARY={'file':'/home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D4.out','loci':361,'completed':0,'independent_undecided':0,'missing_branches':0,'parse_errors':361}
SYMPY_COMPLETENESS_TOTALS={'loci':484,'completed':0,'missing_branches':0,'parse_errors':484}
rc=0
```

## Old probe census: death before first probe, exit 0
```text
$ cd ~/.s11_build/fix2_build && wolframscript -file probe_census.wl /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D2.out /home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D4.out 2>&1 | grep -v 'configuaration'; echo rc=$?

ToExpression::sntxi: Incomplete expression; more input is needed.

rc=0
```

## Sheet-consistent evaluation: the two contested branches of XKIN_ANISO D4 ROOT2_RANK_DROP_JOINT
```text
$ python3 (sheet-consistent minor evaluation; both global sheets; generic rationals)
{muR->bComp}: zero minors sheet(+)=15/16 sheet(-)=15/16
{sRho->1}: zero minors sheet(+)=16/16 sheet(-)=16/16
```

## NOT_APPLICABLE pairing (spec Q8b non-VARIES form)
```text
/home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D2.out: NA_SOLUTION=52 NA_EQUATIONS=52
/home/trevnorris/.s11_build/fix1_build/guarded_XKIN_ANISO_D4.out: NA_SOLUTION=178 NA_EQUATIONS=178
research/pde_ledger_v3/mathematica/out/S11_stray_longitudinal_mathematica_audit.out: NA_SOLUTION=392 NA_EQUATIONS=392
(no LIVE-EQ/NA-SOL lines = zero mismatches)
```

## Census surface: in-class record counts per committed record
```text
$ locus _SOLUTION / _EQUATIONS pairs and undecided-class records
--- research/pde_ledger_v3/mathematica/out/S11_stray_longitudinal_mathematica_audit.out
  LOCUS_SOLUTION lines: 1104
  LOCUS_EQUATIONS lines: 1169
  UNDECIDED status tokens: 814
  REAL_ADMISSIBLE lines: 1083
--- research/pde_ledger_v3/scripts/out/S11_stray_longitudinal_sympy_audit.out
  LOCUS_SOLUTION lines: 379
  LOCUS_EQUATIONS lines: 400
  UNDECIDED status tokens: 0
  REAL_ADMISSIBLE lines: 359
```

## Engine primary budgets (no-starvation floor for probe budgets)
```text
$ grep -nE 'BudgetSeconds = |BudgetBytes = ' research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl
29:qePrimaryBudgetSeconds = 30;
31:qeExactRationalBudgetSeconds = 30;
36:casAttemptMemoryBudgetBytes = 134217728;
39:simplifyPrimaryBudgetSeconds = 45;
43:algebraicPrimaryBudgetSeconds = 45;
45:simplifyFallbackBudgetSeconds = 15;
48:linearAlgebraPrimaryBudgetSeconds = 45;
49:linearAlgebraFallbackBudgetSeconds = 45;
52:linearAlgebraRadicalPrimaryBudgetSeconds = 45;
55:solvePrimaryBudgetSeconds = 30;
56:solveFallbackBudgetSeconds = 30;
59:minorPrimaryBudgetSeconds = 30;
60:minorFallbackBudgetSeconds = 30;
```

