# Measurements — amendment 11 review, round 7 (generated 2026-10-16 08:01)

Generator: `_scratch/s9b_build/gen/s9b_amend11_r7_lookups.sh` (sha256sum/grep/sed/cat/wc only). Line numbers in
`_scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r7.md` are the reviewed round-7 text's (saved before launch); the engines' are the committed files' at
`9af9f705`. The legs' stdout is quoted from their evidence directories.

## The reviewed round-7 text, its scope, and the baseline after both legs
```
$ grep 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend11_review_baseline_r7.sha256
263aa96240a96dc76294c36dc34e38f66559edfd5cc84b79f4825545212d1b30  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ sha256sum _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r7.md
263aa96240a96dc76294c36dc34e38f66559edfd5cc84b79f4825545212d1b30  _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r7.md
```

```
$ git diff -U0 4f56d054 -- research/pde_ledger_v3/directives/S9b_repair_build_directive.md | grep '^@@'
@@ -68,7 +68,3 @@ clear. Dispositions:
@@ -308,2 +304,6 @@ comparator, the production runs and the record.
@@ -429,3 +429 @@ comparator, the production runs and the record.
@@ -454,3 +452 @@ comparator, the production runs and the record.
```

```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_directive_amend11_review_baseline_r7.sha256 | grep -c ': OK$'
21
```

```
$ git diff --stat 9af9f705 -- research/pde_ledger_v3/scripts research/pde_ledger_v3/mathematica | wc -l
0
```

## Verdicts
```
$ grep -m1 'Verdict' _scratch/s9b_build/s9b_amend11_review_r7_claude.md
**Verdict:** NEEDS REVISION
```

```
$ grep -c '^\*\*[0-9]\. ' _scratch/s9b_build/s9b_amend11_review_r7_claude.md
1
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend11_review_r7_grok.txt
Verdict: CLEAR
```

```
$ grep -c 'No findings' _scratch/s9b_build/s9b_amend11_review_r7_grok.txt
1
```

## Claude 1: "shared" fixing equations
```
$ sed -n '306,309p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r7.md
      For value witnesses and compared domain conditions, pair solutions of a shared fixing equation at the same
      solution admitted by both copies, or print and compare the complete sets of (solution, value) pairs, using
      condition evaluations for domains. Never pair different solutions of that equation. Where the fixing
      constructions differ, keep each copy's own solutions.
```

```
$ sed -n '894p;903p;926p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    condition_mass={name:sp.solve(eq(qradial,target),GM)[0] for name,target in inverse_targets.items()}
    mass_density = rho  # K8
        solution=solve_velocity(mass,GM,r)
```

```
$ grep 'SQUARED_IDENTICAL\|SQUARED_FUNCS\|SQUARED_AT_INPUTS\|ZERO_SET\|BRANCH 3 SIGNS\|MATCHED\|CROSSED' _scratch/s9b_build/s9b_amend11_review_r7_claude_evidence/k4c_pairing.stdout
BASE SQUARED_FUNCS ['s9b_delta', 's9b_xi']
K4C SQUARED_FUNCS ['s9b_delta']
SQUARED_IDENTICAL_AS_EMITTED False
BASE SQUARED_AT_INPUTS (s9b_r - 4)**2*(s9b_r - 2)**2/s9b_r**5
BASE ZERO_SET z>b_min {2, 4}
K4C SQUARED_AT_INPUTS (s9b_r - 4)**2*(s9b_r - 2)**2/s9b_r**5
K4C ZERO_SET z>b_min {2, 4}
SQUARED_AT_INPUTS_DIFFERENCE 0
BASE BRANCH 3 SIGNS (-1, 1) z 2 left sqrt(2)/4 right sqrt(2)/4 join True | j_n -3*sqrt(2)/4
BASE BRANCH 3 SIGNS (-1, 1) z 4 left 1/16 right 1/16 join True | j_n -3/16
K4C BRANCH 3 SIGNS (-1, 1) z 2 left sqrt(2)/4 right sqrt(2)/4 join True | j_n -3*sqrt(2)/4
K4C BRANCH 3 SIGNS (-1, 1) z 4 left 1/16 right 1/16 join True | j_n -3/16
MATCHED branch 2 z 2 BASE 3*sqrt(2)/4 K4C 3*sqrt(2)/4 difference 0
MATCHED branch 2 z 4 BASE 3/16 K4C 3/16 difference 0
MATCHED branch 3 z 2 BASE -3*sqrt(2)/4 K4C -3*sqrt(2)/4 difference 0
MATCHED branch 3 z 4 BASE -3/16 K4C -3/16 difference 0
MATCHED branch 4 z 2 BASE 3*sqrt(2)/4 K4C 3*sqrt(2)/4 difference 0
MATCHED branch 4 z 4 BASE 3/16 K4C 3/16 difference 0
MATCHED branch 5 z 2 BASE -3*sqrt(2)/4 K4C -3*sqrt(2)/4 difference 0
MATCHED branch 5 z 4 BASE -3/16 K4C -3/16 difference 0
CROSSED branch 3 BASE z=2 vs K4C z=4: -3*sqrt(2)/4 -3/16 difference 3/16 - 3*sqrt(2)/4
```

```
$ grep 'WHOLE_DOMAIN_SAME_AFTER_INPUTS' _scratch/s9b_build/s9b_amend11_review_r7_claude_evidence/k4c_conjuncts.stdout
  WHOLE_DOMAIN_SAME_AFTER_INPUTS True
  WHOLE_DOMAIN_SAME_AFTER_INPUTS True
  WHOLE_DOMAIN_SAME_AFTER_INPUTS True
  WHOLE_DOMAIN_SAME_AFTER_INPUTS True
```

```
$ grep 'SQUARED_IDENTICAL\|DIFFERENCE_AT' _scratch/s9b_build/s9b_amend11_review_r7_claude_evidence/k9b_squared.stdout
SQUARED_IDENTICAL_AS_EMITTED False
DIFFERENCE_AT_f_EQUALS_0 0
```

```
$ sed -n '52p;572,573p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
xiSlope = Simplify[D[xiAnsatz, r], $Assumptions]; (* K4c *)
      "AtComponentBoundary" -> <|"VelocityLimit" -> Inactive[Limit][value,r -> rBoundary],
        "ExchangeLimit" -> Inactive[Limit][jValue,r -> rBoundary],
```

## Observation: review records in the builders' export
```
$ grep -c 'Under K3 the roots differ' research/pde_ledger_v3/directives/_measurements/S9b_repair_build_directive_amend11_review_disposition.md
1
```

```
$ git ls-tree -r --name-only HEAD research/pde_ledger_v3 | grep -cE '/_measurements/S9b_[^/]*review[^/]*$'
81
```

```
$ grep -n '_measurements/S9b_' _scratch/s9b_build/gen/s9b_r8_exports.sh
18:  grep -E '/(S9b_comparator_|S9b_cross_engine_comparator|test_S9b_cross_engine_comparator|S9b_production_runs|S9b_record_directive|S9b_repair_build_r[0-9]+_review_disposition)|/_measurements/S9b_[^/]*review[^/]*$')
35:  echo "remaining S9b cross-engine/review names: $(git ls-files | grep -cE 'S9b_comparator_|S9b_cross_engine|S9b_production_runs|S9b_record_directive|S9b_repair_build_r[0-9]+_review|_measurements/S9b_[^/]*review' || true)"
```

