# Measurements — O2 build r3 review dispositions (generated 2026-10-07 13:16)

Generator: `_scratch/s9b_build/gen/o2_build_r3_lookups.sh` (grep/sed lookups only; output lines cut at 400 characters; this file is written by it).

## Baseline
```
$ (cd /var/projects/toy_physics && sha256sum -c _scratch/s9b_build/o2_build_review_baseline_r3.sha256)
research/pde_ledger_v3/scripts/O2_exports.py: OK
research/pde_ledger_v3/scripts/O2_live_balance_sympy_audit.py: OK
research/pde_ledger_v3/scripts/O2_live_balance_sympy_ablation.py: OK
research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl: OK
research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_ablation.wl: OK
```

## Leg verdict lines
```
$ (cd /var/projects/toy_physics/_scratch/s9b_build && grep -o -m1 'Verdict: [A-Z ]*' o2_build_review_r3_claude.md o2_build_review_r3_grok.txt)
o2_build_review_r3_claude.md:Verdict: CLEAR
o2_build_review_r3_grok.txt:Verdict: CLEAR
```

## r2 findings: the WL constructs they named (counts in r3)
```
$ grep -c 'operand\[TBulkNormalLive\[s\]' mathematica/O2_live_balance_mathematica_audit.wl
0
```

```
$ grep -c '"Relation"' mathematica/O2_live_balance_mathematica_audit.wl
0
```

```
$ grep -c 'height-chart\|qFace\[s\]\[x1' mathematica/O2_live_balance_mathematica_audit.wl
0
```

```
$ grep -n 'TBulkNormalLive' mathematica/O2_live_balance_mathematica_audit.wl
68:bulkAmplitude = OpenNativeField[operand[TBulkNormalLive, "3.3,5", "7"]][
```

## Leg note — WL K7 replacement and the orientation symbol
```
$ grep -n '"K7"' -A1 mathematica/O2_live_balance_mathematica_ablation.wl
32: {"K7","faceNormal = Table[response[NativeUnitNormalComponent[a],{mapOperand},\n  Join[nativeState,<|\"TangentGeometry\" -> faceTangents|>]],{a,4}]; (* SITE K7 *)",
33-   "faceNormal = orientation[s] unitNormal[D[Append[x,0],#]& /@ x]; (* SITE K7 *)"},
```

```
$ grep -c 'orientation' mathematica/O2_live_balance_mathematica_audit.wl
0
```

```
$ grep -n 'K7' directives/O2_build_directive.md
57:   - **K7, tilted face normal:** the face normal used for bulk-normal loading is replaced by the untilted normal.
```

