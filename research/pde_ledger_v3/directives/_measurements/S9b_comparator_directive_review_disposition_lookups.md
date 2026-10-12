# Measurements — S9b comparator decision-list review (generated 2026-10-12 08:21)

Generator: `_scratch/s9b_build/gen/s9b_comparator_directive_lookups.sh` (sha256sum/grep/sed/cat/wc only). The
reviewed text is kept verbatim at `_scratch/s9b_build/S9b_comparator_build_directive_reviewed.md`.

## The reviewed version
```
$ sha256sum _scratch/s9b_build/S9b_comparator_build_directive_reviewed.md; cat _scratch/s9b_build/s9b_comparator_directive_review_baseline.sha256
0f92d9dab4a59bd7c8602e5d7e5c6da0a26e49615df9cca7ebd95e3da4add137  _scratch/s9b_build/S9b_comparator_build_directive_reviewed.md
0f92d9dab4a59bd7c8602e5d7e5c6da0a26e49615df9cca7ebd95e3da4add137  research/pde_ledger_v3/directives/S9b_comparator_build_directive.md
```

```
$ sha256sum _scratch/s9b_build/s9b_comparator_directive_review_prompt.md
5909e39076c5360b44b49ed8be42146c8689aa71d980204521e267d7fd7e0d17  _scratch/s9b_build/s9b_comparator_directive_review_prompt.md
```

```
$ cat _scratch/s9b_build/s9b_comparator_directive_review_codex.DONE _scratch/s9b_build/s9b_comparator_directive_review_grok.DONE
EXIT=0
EXIT=0
```

```
$ grep -o 'Verdict:\*\* [A-Z ]*' _scratch/s9b_build/s9b_comparator_directive_review_codex_final.txt | head -1; grep -o 'Verdict:\*\* [A-Z ]*' _scratch/s9b_build/s9b_comparator_directive_review_grok.txt | head -1
Verdict:** NEEDS REVISION
Verdict:** NEEDS REVISION
```

## U1 (Codex 1, Grok 2): a control can pass on operand echo alone
```
$ grep -n 'print the SymPy operand, the' _scratch/s9b_build/S9b_comparator_build_directive_reviewed.md
52:4. **Three-valued residuals, printed before any guard.** For each compared value, print the SymPy operand, the
```

```
$ grep -n 'mutation changes the printed output' _scratch/s9b_build/S9b_comparator_build_directive_reviewed.md
74:   mutation changes the printed output; the assertion names no value. Required mutations:
```

```
$ grep -n 'A \*\*repoint control\*\* must show' _scratch/s9b_build/S9b_comparator_build_directive_reviewed.md
43:   itself a zero. A **repoint control** must show that pointing a mapped name at a different object moves the
```

## U2 (Codex 2): the sampling rule avoids every vanishing factor; the engines emit equality strata
```
$ grep -n 'no factor the objects contain vanishes there' _scratch/s9b_build/S9b_comparator_build_directive_reviewed.md
54:   generic rational points chosen so that no factor the objects contain vanishes there (E1; numeric-PIT practice).
```

```
$ grep -n 'including `GM = 0`' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
368:  existence" paragraph hold. On every stratum inside it, including `GM = 0`, each `γ` object (`B_GAMMA_DEFLECTION`,
```

```
$ grep -n 'A zero power of a density or a density ratio is' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
189:     - A zero power of a density or a density ratio is `1` at every value, zero included.
```

## U3 (Grok 1): boolean-valued objects are emitted as computed objects; the list rejects them as operands
```
$ grep -n 'A boolean-valued test is emitted as the CAS object' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
139:5. **No verdicts (D6 B5).** No `VERDICT`, `PASS` or `FAIL`. A boolean-valued test is emitted as the CAS object the
```

```
$ grep -n 'a native boolean, structurally mismatched' _scratch/s9b_build/S9b_comparator_build_directive_reviewed.md
60:   - A comparison that cannot be formed (a native boolean, structurally mismatched or unsupported content) is
```

```
$ grep -n 'a native boolean is rejected as an operand' _scratch/s9b_build/S9b_comparator_build_directive_reviewed.md
86:   Two further tests: a duplicate join or name row is rejected, and a native boolean is rejected as an operand while
```

```
$ grep -c '^WL_S9B_A_NONRECIPROCAL_PATH_DEPENDENCE' research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out
1
```

```
$ grep -A3 '^WL_S9B_A_NONRECIPROCAL_PATH_DEPENDENCE' research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out | grep -o 'ConditionalExpression\[False' | head -1
ConditionalExpression[False
```

```
$ grep '^PY_S9B_C_CONSTANT_N_DEPENDENCE' research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out
PY_S9B_C_CONSTANT_N_DEPENDENCE: Tuple(Tuple(Integer(0), Integer(0)), Tuple(Str('ATTACHED_DOMAIN'), true))
```

