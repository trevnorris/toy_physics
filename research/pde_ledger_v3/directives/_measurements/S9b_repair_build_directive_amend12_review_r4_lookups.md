# Measurements — amendment 12 review, round 4 (generated 2026-10-16 19:57)

Generator: `_scratch/s9b_build/gen/s9b_amend12_r4_lookups.sh` (sha256sum/grep/sed/cat/git show only). Line
numbers in `_scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r4.md` are the reviewed round-4 text's (saved before launch); the engines' and harnesses' are the
committed files' at `86097a44`.

## The reviewed text and the baseline
```
$ grep 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend12_review_baseline_r4.sha256 | cut -c1-64
498d85a0ff9ff165b627fc8abcfbecd2ca5351406beffa6a8ba200c7f77f1dde
```

```
$ sha256sum _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r4.md | cut -c1-64
498d85a0ff9ff165b627fc8abcfbecd2ca5351406beffa6a8ba200c7f77f1dde
```

```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_directive_amend12_review_baseline_r4.sha256 | grep -c ': OK$'
14
```

```
$ sha256sum _scratch/s9b_build/s9b_repair_build_directive_amend12_review_prompt_r4.md | cut -c1-64
358d6dbef7641bc3194a65723b4f5fdefcfb3dcc0b8a417b6f3c931793ea8e72
```

```
$ git diff --stat 6a78b956 -- research/pde_ledger_v3/directives/S9b_repair_build_directive.md
 research/pde_ledger_v3/directives/S9b_repair_build_directive.md | 9 ++++++---
 1 file changed, 6 insertions(+), 3 deletions(-)
```

## Verdicts
```
$ head -1 _scratch/s9b_build/s9b_amend12_review_r4_claude.md
**Verdict: NEEDS REVISION**: one finding about the rule, plus one smaller one about wording.
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend12_review_r4_grok.txt
Verdict: NEEDS REVISION
```

## The round-4 clauses
```
$ sed -n '300,325p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r4.md
      Include reached text only when the engine writes it as a literal status or label with no mathematical reading.
      Exclude exactly two additional kinds of part in both engines, printing each reason: lists of symbol, function
      or row names extracted from expressions, whether emitted as text or symbols; and text parts whose literal
      difference is the knife's replacement source. For text parts replace mathematical witness requirements with
      literal comparison: print both texts
      and location, count the part tried whether equal or different, and count a literal difference as a witness in
      both engines. One-copy text parts follow the one-copy rule.

      For mathematical parts count a part tried only when all witness requirements except difference hold. Decide
      values and conditions by engine CAS; substitution alone and unresolved evaluations do not count. Count an attempt
      only after running its witness test's CAS operations at printed inputs: solve fixing equations and evaluate values
      and attached conditions. Before listing an untried part, print the actual unresolved evaluation or requirement the
      CAS decided cannot be met; without those evaluations, continue the search. Print each copy's CAS-evaluated output
      of the construction defined by the knife list at trial inputs; count no attempt at inputs giving the same
      output. Limit each printed unmet requirement to the printed inputs.

      Try item 9's shared object on each sign even if stored once; label trials do not try it. Count a sign trial
      only at a printed concrete input or throughout a printed subrange of the part's stated interval, with printed
      density evaluations in that sign in each emitting copy. Respect the spec's profile definitions, including ρ₀
      as asymptotic bulk density. If no such input meets the sign on that interval, print the CAS decision and unmet
      requirement; retain the part as untried, with no witness.

      End each search after finitely many attempts with exactly one outcome: a witness; otherwise, after attempting
      every part, **none found** if all count as tried, or **undecided**, listing each untried part and unmet
      requirement. Print each attempt's tag, location, role, inputs, evaluations or **not emitted**, whether tried, and
      unmet requirements.
```

## My round-4 entry 2, which line 314 follows
```
$ sed -n '284,285p' research/pde_ledger_v3/directives/_measurements/S9b_repair_build_directive_amend12_review_disposition.md
2. **Inert inputs (R).** A trial at inputs where the knife-list construction gives the same output in both copies is
   not an attempt. A printed unmet requirement holds at the printed inputs; it does not stand for the part as a whole.
```

## Claude 1, Grok 1: the K10 gate is a conjunction with unmutated conjuncts
```
$ sed -n '166p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
branchGateLocal = propagatingLocal; (* K10 *)
```

```
$ sed -n '770,771p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    observable_gate = real_domain  # K10: existence only
    branch_domain=sp.And(observable_gate,travel_domain)
```

```
$ cat _scratch/s9b_build/s9b_amend12_review_r4_claude_evidence/wl_k10_masking.stdout
BASE_DOMAIN InputForm[r > 0 && muPerp[r]/rhoBr[r] > velocityProfile[r]^2*(1 + Derivative[1][xiProfile][r]^2) && Element[muPerp[r]/rhoBr[r], Reals] && muPerp[r]/rhoBr[r] > 0 && r^2*(1 + Derivative[1][xiProfile][r]^2) > 0 && c0*(1 + deltaProfile[r]) == Sqrt[muPerp[r]/rhoBr[r]]]
K10_DOMAIN InputForm[r > 0 && muPerp[r]/rhoBr[r] > velocityProfile[r]^2*(1 + Derivative[1][xiProfile][r]^2) &&  !(Element[muPerp[r]/rhoBr[r], Reals] && muPerp[r]/rhoBr[r] > 0 && r^2*(1 + Derivative[1][xiProfile][r]^2) > 0) && c0*(1 + deltaProfile[r]) == Sqrt[muPerp[r]/rhoBr[r]]]
INPUTS r=2 c0=1 delta=0 V=0 xi=0 rhoBr=1 muPerp=4
  CONSTRUCTION_OUTPUT base=True k10=False
  DOMAIN base=False k10=False
  SPEED_CHART_CONJUNCT False
INPUTS r=2 c0=1 delta=0 V=0 xi=0 rhoBr=1 muPerp=1
  CONSTRUCTION_OUTPUT base=True k10=False
  DOMAIN base=True k10=False
  SPEED_CHART_CONJUNCT True
DONE
```

```
$ grep 'K10_' _scratch/s9b_build/s9b_amend12_review_r4_grok_evidence/site_parts_and_attempts.stdout
K10_BASE_GATE real_domain & travel_domain
K10_COPY_GATE travel_domain & ~real_domain
K10_KIND_BASE And
K10_IS_NAME_LIST False
K10_NONINERT_FAIL {'construction_outputs': (True, False), 'gates': (False, False)}
K10_WITNESS_INPUT {'construction_outputs': (True, False), 'gates': (True, False)}
```

## Grok 1: an attempt at inputs where a requirement fails (K8, K9b)
```
$ sed -n '903,904p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    mass_density = rho  # K8
    exchange=sp.simplify(-sp.diff(area*mass_density*velocity,r)/area)
```

```
$ sed -n '44p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
    'K9b': ('power_density = density', 'power_density = rho0'),
```

```
$ grep 'K8_\|RHO_\|NONPOSITIVE\|FAR_RAY\|K9B_' _scratch/s9b_build/s9b_amend12_review_r4_grok_evidence/site_parts_and_attempts.stdout
K8_LIVE_BASELINE 1/(8*s9b_r**4)
K8_LIVE_COPY 0
K8_AT_R2_BASELINE 1/128
K8_AT_R2_COPY 0
K8_AT_R2_DIFFERENCE 1/128
K8_INERT_BASELINE 0
K8_INERT_COPY 0
K8_INERT_DIFFERENCE 0
RHO_LIMIT_INF rho0
RHO_AT_5 -2*rho0
NONPOSITIVE_SET Interval(5 - sqrt(2), sqrt(2) + 5)
FAR_RAY_SUBSET_OF_NONPOSITIVE False
NONPOSITIVE_NONEMPTY True
K9B_SAMPLE 5 RHO_OVER_RHO0 -2 SIGN_NONPOSITIVE True CONSTRUCTION_OUTPUTS (-2, 1) OUTPUTS_EQUAL False
K9B_SAMPLE 100 RHO_OVER_RHO0 9023/9026 SIGN_NONPOSITIVE False CONSTRUCTION_OUTPUTS (9023/9026, 1) OUTPUTS_EQUAL False
```

## Exclusions (round-4 entry 1)
```
$ grep -c 'lines containing full replacement source: 0' _scratch/s9b_build/s9b_amend12_review_r4_claude_evidence/replacement_source_grep.txt
25
```

```
$ grep 'CLASS' _scratch/s9b_build/s9b_amend12_review_r4_grok_evidence/site_parts_and_attempts.stdout
CLASS export_closure_row_names EXCLUDE_NAME_LIST SHOWN ['s9b_deflection', 's9b_constant_density']
CLASS wl_remaining_inputs EXCLUDE_NAME_LIST SHOWN [deltaProfile, xiProfile]
CLASS wl_free_inputs EXCLUDE_NAME_LIST SHOWN [velocityProfile, rhoBr, deltaProfile]
CLASS integration_constants_written KEEP_VALUE_OR_DOMAIN SHOWN [C]
CLASS k10_gate KEEP_VALUE_OR_DOMAIN SHOWN travel_domain & ~real_domain
CLASS k8_exchange KEEP_VALUE_OR_DOMAIN SHOWN 0
CLASS k5_factor KEEP_VALUE_OR_DOMAIN SHOWN 0
CLASS wl_density_value KEEP_VALUE_OR_DOMAIN SHOWN rhoBrConstant
CLASS k7_status KEEP_TEXT SHOWN ('DIFFERENTIAL', 'velocity undetermined by one fixed-b condition')
CLASS replacement_source_text EXCLUDE_REPLACEMENT_TEXT SHOWN ('observable_gate = sp.Not(real_domain)',)
```

