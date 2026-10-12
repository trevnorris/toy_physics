# Measurements — S9b comparator decision list (generated 2026-10-12 08:22)

Generator: `_scratch/s9b_build/gen/s9b_comparator_directive_measurements.sh` (ls/grep/sed/sort/comm/wc/git only).
They ground what `directives/S9b_comparator_build_directive.md` says about its inputs.

## The inputs exist; the transcripts are annex content committed at 503736ff
```
$ git ls-files -s research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out | cut -c1-6; ls -l research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out | cut -c1-80
120000
120000
lrwxrwxrwx 1 trevnorris trevnorris 140 Oct 12 08:00 research/pde_ledger_v3/mathe
lrwxrwxrwx 1 trevnorris trevnorris 144 Oct 12 08:00 research/pde_ledger_v3/scrip
```

```
$ git log --oneline -1 -- research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out
503736ff Add the S9b production transcripts (reviewed r5 engines)
```

```
$ ls research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py research/pde_ledger_v3/scripts/test_O2_cross_engine_comparator.py
research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py
research/pde_ledger_v3/scripts/test_O2_cross_engine_comparator.py
```

```
$ ls research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py research/pde_ledger_v3/scripts/test_S9b_cross_engine_comparator.py
ls: cannot access 'research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py': No such file or directory
ls: cannot access 'research/pde_ledger_v3/scripts/test_S9b_cross_engine_comparator.py': No such file or directory
```

## The engines are unchanged since bfa034c4
```
$ git diff --quiet bfa034c4 -- research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl && echo ENGINES_UNCHANGED_SINCE_bfa034c4
ENGINES_UNCHANGED_SINCE_bfa034c4
```

## Both transcripts emit the same set of shared <QUANTITY> names (item 1)
```
$ grep -o '^PY_S9B_[A-Z0-9_]*' research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out | sed 's/^PY_S9B_//' | sort -u > /tmp/s9b_cmpdir_py.txt; grep -o '^WL_S9B_[A-Z0-9_]*' research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out | sed 's/^WL_S9B_//' | sort -u > /tmp/s9b_cmpdir_wl.txt; wc -l < /tmp/s9b_cmpdir_py.txt; wc -l < /tmp/s9b_cmpdir_wl.txt; comm -3 /tmp/s9b_cmpdir_py.txt /tmp/s9b_cmpdir_wl.txt | wc -l
330
330
0
```

```
$ grep -c '^PY_S9B_' research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out; grep -c '^WL_S9B_' research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out
330
330
```

```
$ grep -c '^PY_LOCAL_S9B_' research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out; grep -c '^WL_LOCAL_S9B_' research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out
49
44
```

## The build-directive items the list cites
```
$ grep -n '^4\. \*\*Tags\|^5\. \*\*No verdicts\|^9\. \|^10\. \|^14\. \*\*Ablation harness' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
97:4. **Tags (D6 B7).** Use the grammar `<ENGINE>_S9B_<QUANTITY>`, where `<ENGINE>` is `PY` or `WL`. Both engines emit
139:5. **No verdicts (D6 B5).** No `VERDICT`, `PASS` or `FAIL`. A boolean-valued test is emitted as the CAS object the
170:9. **Part C domains (D6 B6).** Each Part C row's domain predicates are built from that row's own response. They are
192:10. **Forward case and single-mechanism restrictions (user requests, 2026-10-08).**
264:14. **Ablation harness** (`docs/development_pipeline.md` §4).
```

```
$ git show 34a994d5 --stat | head -3
commit 34a994d5b3dbcdc294b688dd318a52371773ef49
Author: Trevor Norris <trev.norris@gmail.com>
Date:   Fri Oct 9 18:56:16 2026 -0600
```

