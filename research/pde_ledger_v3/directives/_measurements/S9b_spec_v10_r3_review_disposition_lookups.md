# Measurements — S9b spec v10 review round 3 dispositions (generated 2026-10-08 13:52)

Generator: `_scratch/s9b_build/gen/s9b_spec_v10_r3_lookups.sh` (sed/grep/sha256sum/git only; lines cut at 330 characters).

```
$ sha256sum -c _scratch/s9b_build/s9b_spec_v10_review_baseline_r2.sha256
research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md: OK
research/pde_ledger_v3/directives/_measurements/S9b_v10_spec_lookups.md: OK
```

```
$ grep -n -o 'Verdict:[* ]*[A-Z][A-Z ]*' _scratch/s9b_build/s9b_spec_v10_review_r3_claude.md _scratch/s9b_build/s9b_spec_v10_review_r3_grok.txt
_scratch/s9b_build/s9b_spec_v10_review_r3_claude.md:3:Verdict: NEEDS REVISION
_scratch/s9b_build/s9b_spec_v10_review_r3_grok.txt:1:Verdict: CLEAR
```

## F1 — the density re-expression reading: undefined correction objects, exempted from the order restriction
```
$ sed -n '77,81p;90p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
    These are supplied identifications of the reading, using the supplied metric; they do not
    replace the mass law. Print the re-expressed densities and correction objects with their
    orders under the supplied metric and slope counting. No additional independent bound on the
    drain divergence is a premise of this density re-expression. The imposed-law restriction
    below does not apply to this reading.
    comparison (D4). The historical derivative condition is not imposed here.
```

```
$ grep -n 'unless it states a condition' research/pde_ledger_v3/directives/S9b_repair_decision_list.md
66:  unless it states a condition that bounds that correction relative to `j_n` itself.
```

```
$ git diff 478b5e7b -- research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md | grep -c '^+.*correction objects'
1
```

The leg's script output:
```
$ grep -n 'contains xi\|flat_div\|\[witness\] j_coord' _scratch/s9b_build/s9b_spec_v10_review_r3_claude_evidence/induced_measure_readings.stdout
10:[re-expression] relative correction contains xi''(R)? False
12:[re-expression] flat_div(rho_g V) + j_g = -V_r(R)*rho(R)*Derivative(xi(R), R)*Derivative(xi(R), (R, 2))/(Derivative(xi(R), R)**2 + 1)**(3/2)
13:[re-expression] flat_div(rho_g V) + j_g contains xi''(R)? True
16:[imposed] j_imp - j_coord contains xi''(R)? True
21:[witness] j_coord = 0
24:[witness] re-expression flat_div(rho_g V) + j_g on witness = -A*k/(2*R**(5/2)*(R + k)**(3/2))
26:[FORM control] modified re-expression correction contains xi''(R)? True
```

## Routing note — the spec no longer carries the build-deferral and stop sections
```
$ grep -c 'Deferred to the build\|Stop and report' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
0
```

```
$ git show c2f1cf2b:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md | grep -n 'Deferred to the build\|Stop and report'
185:- **Deferred to the build** (implementation, not new physics; the build directive owns it):
189:- **Stop and report**, without choosing, when any of these happens:
```

