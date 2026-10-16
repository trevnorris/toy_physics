# Measurements — amendment 11 review, round 0 (generated 2026-10-16 04:01)

Generator: `_scratch/s9b_build/gen/s9b_amend11_r0_lookups.sh` (sha256sum/grep/sed/cat/wc only). Line numbers in
`_scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r0.md` are the reviewed round-0 text's; the engines' are the committed engines'.

## The reviewed round-0 text
It was rebuilt from the Codex leg's captured `git diff 4b4db191` applied to the directive at `4b4db191`; its hash is
the round-0 baseline's.
```
$ grep 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend11_review_baseline_r0.sha256
0329a03d98f12563516c8c075d802ec6497fc52bb2c64574edc9efdb96498b82  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ sha256sum _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r0.md
0329a03d98f12563516c8c075d802ec6497fc52bb2c64574edc9efdb96498b82  _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r0.md
```

```
$ cat _scratch/s9b_build/s9b_amend11_review_r0_codex_evidence/amendment.command
cd /var/projects/toy_physics
git diff 4b4db191 -- research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

The other baseline entries, after both legs:
```
$ grep -v 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend11_review_baseline_r0.sha256 | sha256sum -c | grep -c ': OK$'
9
```

## Verdicts
```
$ cat _scratch/s9b_build/s9b_amend11_review_r0_codex_final.txt | head -3
**Verdict: CLEAR.**

**Findings:** None.
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend11_review_r0_grok.txt
Verdict: NEEDS REVISION
```

```
$ grep -c '^\*\*[0-9]\. ' _scratch/s9b_build/s9b_amend11_review_r0_grok.txt
3
```

## G1: the round-0 relaxation, the header, and rule 2
```
$ grep -n 'Any other part that differs may print undecided once the harness has looked for a witness' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r0.md
325:      Any other part that differs may print undecided once the harness has looked for a witness.
```

```
$ grep -n 'change. Any other part may print undecided' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r0.md
73:change. Any other part may print undecided. Parts 2 and 3 point at it. The engines do not change in this round.
```

```
$ grep -n 'Otherwise \*\*a change,\*\* where the engine shows the part changed, symbolically or by a witness' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r0.md
304:      2. Otherwise **a change,** where the engine shows the part changed, symbolically or by a witness. A witness
```

## G2: a value compared apart from its domain, and item 9's values
```
$ grep -n 'A value is compared apart from its domain' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r0.md
294:      function, such a point under some admissible profile. A value is compared apart from its domain: the point
```

```
$ grep -n 'has a value only where both' research/pde_ledger_v3/directives/S9b_repair_build_directive.md; grep -n 'A principal root of a negative square is not a value' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
197:       `c_s0` are the non-negative roots of the spec's `c_s²` at `ρ` and at `ρ₀`. It has a value only where both
198:       squares are real and non-negative and `c_s0 > 0`. A principal root of a negative square is not a value of
```

```
$ sed -n '97,102p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
responseInputDomains = {rho0 > 0,
  rho0 > 0 && Element[{localSoundSquared,asymptoticSoundSquared},Reals] &&
    localSoundSquared >= 0 && asymptoticSoundSquared > 0, rho0 > 0};
responseValues = {c0, c0 Sqrt[localSoundSquared]/Sqrt[asymptoticSoundSquared],
  c0 densityPower[localBulkDensity/rho0,s]};
responseInputs = MapThread[ConditionalExpression,{responseValues,responseInputDomains}];
```

```
$ grep -h -A2 'CASE sqrt(x\*\*2) vs x' _scratch/s9b_build/s9b_amend11_review_r0_grok_evidence/s9b-amend11-r0-defs/classify_defs.stdout
CASE sqrt(x**2) vs x, attached domain x>0 not required
 point 1 ('agree', ('real', 1), ('real', 1))
 point -1 ('differ', ('real', 1), ('real', -1))
```

```
$ grep -h -A1 'CASE item-9 principal root' _scratch/s9b_build/s9b_amend11_review_r0_grok_evidence/s9b-amend11-r0-defs/classify_defs.stdout
CASE item-9 principal root of a negative square
 sqrt(-1) = I is_real False
```

## G3: K10's gate, and K11's residual
```
$ grep -n 'K10, branch gate' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
350:    - **K10, branch gate:** the computed branch-existence condition is replaced by its logical negation before the
```

```
$ sed -n '719,720p;770p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    progress('BRANCH_EXISTENCE')
    emit('BRANCH_EXISTENCE',real_domain)
    observable_gate = real_domain  # K10: existence only
```

```
$ sed -n '166p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
branchGateLocal = propagatingLocal; (* K10 *)
```

```
$ grep -n 'the unmutated compact evaluation minus the full baseline payload' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
366:      - the method's residual: the unmutated compact evaluation minus the full baseline payload reduced by the same
```

```
$ grep -n 'A tag the method does not reach is printed as not evaluated' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
370:    - A tag the method does not reach is printed as not evaluated.
```

```
$ grep -n 'if any part the engine computes from the mutated' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r0.md
321:      - for each knife whose construction the engine contains: if any part the engine computes from the mutated
```

Lines of the directive containing "route" (two-route objects):
```
$ grep -c -i 'route' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
0
```

## Codex leg's evidence files
```
$ ls _scratch/s9b_build/s9b_amend11_review_r0_codex_evidence | grep -c ''
70
```

```
$ ls _scratch/s9b_build/s9b_amend11_review_r0_codex_evidence | grep -E '\.py$|\.wl$'
capture_files.py
k11_coverage.py
semantic_checks.py
semantic_checks_v2.wl
semantic_checks.wl
write_report.py
```

