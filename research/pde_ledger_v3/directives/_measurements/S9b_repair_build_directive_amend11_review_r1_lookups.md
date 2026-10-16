# Measurements — amendment 11 review, round 1 (generated 2026-10-16 04:26)

Generator: `_scratch/s9b_build/gen/s9b_amend11_r1_lookups.sh` (sha256sum/grep/sed/cat only). Line numbers in
`_scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r1.md` are the reviewed round-1 text's (saved before launch); the engines' are the committed engines'.

## The reviewed round-1 text, and the baseline after both legs
```
$ grep 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend11_review_baseline_r1.sha256
75f4d5625769db7651e1944532ffd47ff6585f233b2620a4ba64a21bbe8e75d0  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ sha256sum _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r1.md
75f4d5625769db7651e1944532ffd47ff6585f233b2620a4ba64a21bbe8e75d0  _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r1.md
```

```
$ grep -v 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend11_review_baseline_r1.sha256 | sha256sum -c | grep -c ': OK$'
12
```

## Verdicts
```
$ head -1 _scratch/s9b_build/s9b_amend11_review_r1_codex_final.txt
**Verdict: CLEAR.**
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend11_review_r1_grok.txt
Verdict: NEEDS REVISION
```

```
$ grep -c '^\*\*[0-9]\. ' _scratch/s9b_build/s9b_amend11_review_r1_grok.txt
2
```

## H1: "shows that it holds" and the domains the engines attach
```
$ grep -n 'A condition holds at a' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r1.md
292:      part that carries a profile function, such a point under some admissible profile. A condition holds at a
```

```
$ grep -n 'point only where the harness shows that it holds there' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r1.md
293:      point only where the harness shows that it holds there, whether it is stated as a formula or in words (such
```

```
$ sed -n '400,404p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
regularityDomain = <|"Radius" -> r > bFar, "Anchor" -> rAnchor > bFar,
  "Profiles" -> {deltaProfile,velocityProfile,xiProfile,rhoBr,fProfile},
  "FunctionalDomain" -> {"required derivatives exist", "displayed improper integrals converge",
    "Abel inverse exists on the exterior half-line"},
  "OpticalCounting" -> {"delta O(1/r)", "V/c0 O(1/Sqrt[r])", "xiPrimeSquared O(1/r)"}|>;
```

```
$ sed -n '262p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
propagationGate = rayConstruction && Inactive[ForAll][r, traversedRadii,
```

```
$ cat _scratch/s9b_build/s9b_amend11_review_r1_grok_evidence/resolve_wl.stdout
RESOLVE True
TRUEQ_RESOLVE True
TRUEQ_BARE False
```

```
$ grep -h 'summary' _scratch/s9b_build/s9b_amend11_review_r1_grok_evidence/classify_shows.stdout | sort | uniq -c
      3   summary CHANGED
      9   summary UNDECIDED
```

## H2: K10's two definitions of its site parts
```
$ grep -n 'A knife.s site parts are the parts of the emitted objects whose computation uses the mutated construction' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r1.md
324:      A knife's site parts are the parts of the emitted objects whose computation uses the mutated construction,
```

```
$ grep -n 'For K10 they are the' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r1.md
326:      which any edit of the source changes, or a part computed before the mutation enters. For K10 they are the
```

```
$ sed -n '166p;174p;182p;184p;275p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
branchGateLocal = propagatingLocal; (* K10 *)
physicalPropagation = branchGateLocal /. branchRules;
  {"unable to traverse in a required direction", branchGateLocal && Not[subcritical]}},
emit["LOCAL_BRANCH_TYPES", branchTypes];
  emit[prefix <> "BRANCH_EXISTENCE", wrapper[propagationGate /. rules]];
```

```
$ sed -n '720p;770p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    emit('BRANCH_EXISTENCE',real_domain)
    observable_gate = real_domain  # K10: existence only
```

```
$ grep -h '^  AGREE\|^  DIFFER' _scratch/s9b_build/s9b_amend11_review_r1_grok_evidence/k10_sites.stdout
  AGREE sympy BRANCH_EXISTENCE: uses=False gated_on=False
  AGREE sympy F_*_BRANCH_EXISTENCE: uses=False gated_on=False
  AGREE sympy A_NONRECIPROCAL_PATH_DEPENDENCE domain: uses=True gated_on=True
  AGREE sympy B_GAMMA_* domain: uses=True gated_on=True
  DIFFER wolfram LOCAL_BRANCH_TYPES: uses=True gated_on=False
  DIFFER wolfram BRANCH_EXISTENCE: uses=True gated_on=False
  DIFFER wolfram LOCAL_NONRECIPROCAL_CONSTRUCTION domain: uses=True gated_on=False
  AGREE wolfram observables passed through gated: uses=True gated_on=True
```

