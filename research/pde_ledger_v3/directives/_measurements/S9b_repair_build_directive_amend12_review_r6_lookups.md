# Measurements — amendment 12 review, round 6 (generated 2026-10-16 21:41)

Generator: `_scratch/s9b_build/gen/s9b_amend12_r6_lookups.sh` (sha256sum/grep/sed/cat/git show only). Line
numbers in `_scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r6.md` are the reviewed round-6 text's (saved before launch); the engines' and harnesses' are the
committed files' at `86097a44`.

## The reviewed text and the baseline
```
$ grep 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend12_review_baseline_r6.sha256 | cut -c1-64
c050127bf547518275810c67cabd52122f0635834c1902a6a50eecd6d186ce8f
```

```
$ sha256sum _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r6.md research/pde_ledger_v3/directives/S9b_repair_build_directive.md | cut -c1-64
c050127bf547518275810c67cabd52122f0635834c1902a6a50eecd6d186ce8f
c050127bf547518275810c67cabd52122f0635834c1902a6a50eecd6d186ce8f
```

```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_directive_amend12_review_baseline_r6.sha256 | grep -c ': OK$'
16
```

```
$ sha256sum _scratch/s9b_build/s9b_repair_build_directive_amend12_review_prompt_r6.md | cut -c1-64
f367e62ac3a907b5d355fa055c2660073c07e99af2302fe4573fb474752a4bf7
```

```
$ git diff -U0 e506793f -- research/pde_ledger_v3/directives/S9b_repair_build_directive.md | grep '^@@'
@@ -284 +284,2 @@ comparator, the production runs and the record.
@@ -293,3 +294,4 @@ comparator, the production runs and the record.
@@ -370,2 +372,3 @@ comparator, the production runs and the record.
@@ -408,2 +411,3 @@ comparator, the production runs and the record.
@@ -412 +416 @@ comparator, the production runs and the record.
```

## Verdicts
```
$ head -1 _scratch/s9b_build/s9b_amend12_review_r6_claude.md
**Verdict: CLEAR.** I have no findings.
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend12_review_r6_grok.txt
Verdict: CLEAR
```

## Entry 1: K11's reach
```
$ sed -n '294,295p;372,374p;411,416p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r6.md
      Take the union across copies without filtering by structural marks. For K11 include every tag whose computation
      the mutation enters, whether or not the compact method reaches it. Print the set and all exclusions with reasons.
      object, and are neither a profile choice nor a profile family under item 13. Require K11 no-change decisions
      over every admissible input; use neither the method's sample domain nor truncation to restrict that coverage
      or supply an attached condition.
      Print both pairs as ordered operands with the tags' structural marks, without subtraction or mathematical
      equality or zero/change claims; print search decisions separately. Only a witness may claim change; the residual
      cannot furnish one.
    - Each pair carries its coverage: retained grades and truncation order or sample domain. Inexact evaluations
      carry their precision and proved error enclosures as above.
    - Print a tag the method does not reach as not evaluated and not compared, without a no-change decision.
```

```
$ sed -n '322,325p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r6.md
      Stop at the first witness and print it. Otherwise end after finitely many attempts at every site part with
      exactly one outcome: **none found** only with a no-change decision at every part; **undecided** otherwise,
      listing every part without one and what the CAS left unsettled or decided there. Print each attempt's tag,
      location, role, inputs, CAS operations and evaluations or **not emitted**, and unmet requirements.
```

## Entry 2: equality claims
```
$ sed -n '283,285p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r6.md
      equal; their printing preserves every distinction. Otherwise mark it **not compared**. Marks claim only
      structural identity; any claim of change requires a witness. Infer no mathematical equality from these marks.
      Print no verdict.
```

```
$ grep -c 'Claim no mathematical equality' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r6.md
0
```

```
$ grep -n 'zero/change' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r6.md
412:      equality or zero/change claims; print search decisions separately. Only a witness may claim change; the residual
```

## The legs' K11 reach evidence
```
$ grep 'aphi\|arr' _scratch/s9b_build/s9b_amend12_review_r6_grok_evidence/k10_k11_coverage.stdout
field arr mut_has_azimuthal False identical_to_baseline True
field aphi mut_has_azimuthal True identical_to_baseline False
SY 678: beta,arr,aphi = (geometry[k] for k in ('beta','arr','aphi'))
```

```
$ grep -c '' _scratch/s9b_build/s9b_amend12_review_r6_claude_evidence/py_k11_dataflow.stdout _scratch/s9b_build/s9b_amend12_review_r6_claude_evidence/wl_k11_dataflow.stdout
_scratch/s9b_build/s9b_amend12_review_r6_claude_evidence/py_k11_dataflow.stdout:47
_scratch/s9b_build/s9b_amend12_review_r6_claude_evidence/wl_k11_dataflow.stdout:20
```

