# Measurements — amendment 12 review, round 5 (generated 2026-10-16 21:19)

Generator: `_scratch/s9b_build/gen/s9b_amend12_r5_lookups.sh` (sha256sum/grep/sed/cat/git show only). Line
numbers in `_scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r5.md` are the reviewed round-5 text's (saved before launch); the engines' and harnesses' are the
committed files' at `86097a44`.

## The reviewed text and the baseline
```
$ grep 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend12_review_baseline_r5.sha256 | cut -c1-64
4583d7e32f05f9064b431503f2edba16fc4de4d5bfea74af85e52766f2534c70
```

```
$ sha256sum _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r5.md | cut -c1-64
4583d7e32f05f9064b431503f2edba16fc4de4d5bfea74af85e52766f2534c70
```

```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_directive_amend12_review_baseline_r5.sha256 | grep -c ': OK$'
15
```

```
$ sha256sum _scratch/s9b_build/s9b_repair_build_directive_amend12_review_prompt_r5.md | cut -c1-64
63825c09e6656fcb8654e95e7b0c11754b8a2be465aad585b0797d5aa25ef2d7
```

```
$ git diff --stat 10555f5d -- research/pde_ledger_v3/directives/S9b_repair_build_directive.md
 .../directives/S9b_repair_build_directive.md       | 34 ++++++++++------------
 1 file changed, 16 insertions(+), 18 deletions(-)
```

## Verdicts
```
$ head -1 _scratch/s9b_build/s9b_amend12_review_r5_claude.md
**Verdict: NEEDS REVISION.** There are two narrow findings. One is a contradiction that round 5 introduced. The other is a K11 rule gap that round 5's outcome rule makes load-bearing.
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend12_review_r5_grok.txt
Verdict: CLEAR
```

## The round-5 clauses
```
$ sed -n '303,325p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r5.md
      difference is the knife's replacement source. For text parts replace mathematical tests with engine CAS literal
      comparison: print both texts, location and decision; use a literal difference as a witness and literal equality
      as **no change** in both engines. One-copy text parts follow the one-copy rule.

      For mathematical parts report **no change** only on an engine CAS decision of agreement between the copies in
      the part and its attached conditions at every admissible input. Report **cannot be tried** only on an
      engine CAS decision that no admissible input meets the witness requirements apart from difference. For every
      no-change or cannot-be-tried decision, print its statement, CAS result and full input coverage: all inputs and
      live profiles admitted by the spec's definitions and supplied premises on the part's stated interval. Infer
      neither decision from samples or unresolved evaluations.

      Try item 9's shared object on each sign even if stored once; label trials do not try it. Count a sign trial
      only at a printed concrete input or throughout a printed subrange of the part's stated interval, with printed
      density evaluations in that sign in each emitting copy. Respect the spec's profile definitions, including ρ₀
      as asymptotic bulk density. If the CAS decides no such input meets the sign on that interval, retain the part
      as **cannot be tried**, printing the unmet requirement.

      Stop at the first witness and print it. Otherwise end after finitely many attempts at every site part with
      exactly one outcome: **none found** only with a no-change decision at every part; **undecided** otherwise,
      listing every part without one and what the CAS left unsettled or decided there. Print each attempt's tag,
      location, role, inputs, CAS operations and evaluations or **not emitted**, and unmet requirements.

      Attached conditions are the domains, gates and premises stating where the part holds, at every containing
```

## Claude 2: the older sentences against the new no-change report
```
$ sed -n '284p;408,409p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r5.md
      structural identity; any claim of change requires a witness. Claim no mathematical equality. Print no verdict.
      Both pairs print as ordered operands with the same structural marks as tags, without subtraction or any
      zero/change classification. Only a witness may claim change; the residual cannot furnish one.
```

```
$ git show 75c500db:research/pde_ledger_v3/directives/S9b_repair_build_directive.md | grep -n 'Claim no mathematical equality\|zero/change classification'
281:      structural identity; any claim of change requires a witness. Claim no mathematical equality. Print no verdict.
373:      zero/change classification. Only a witness may claim change; the residual cannot furnish one.
```

```
$ git show 10555f5d:research/pde_ledger_v3/directives/S9b_repair_build_directive.md | grep -c 'no change'
0
```

```
$ grep -c 'no change' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r5.md
2
```

## Claude 1: K11's site set and the compact method's reach
```
$ sed -n '289,293p;369,370p;410,412p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r5.md
      is the gate after negation. For K11, the witness compares compact corrupted/unmutated parts of a reached tag
      whose computation the mutation enters, never the method's residual.

      Enumerate all direct and transitive site parts, including engine-local tags and both item 9 signs.
      Take the union across copies (K11's compact pair), without filtering by structural marks. Print the set and
      all required derivatives at the witness point. Witness inputs are chosen for this test only, enter no engine
      object, and are neither a profile choice nor a profile family under item 13. K11's sample domain is coverage,
    - Each pair carries its coverage: retained grades and truncation order or sample domain. Inexact evaluations
      carry their precision and proved error enclosures as above.
    - A tag the method does not reach is printed as not evaluated and not compared.
```

```
$ git show 75c500db:research/pde_ledger_v3/directives/S9b_repair_build_directive.md | grep -n 'sample domain is coverage\|method does not reach'
334:      object, and are neither a profile choice nor a profile family under item 13. K11's sample domain is coverage,
376:    - A tag the method does not reach is printed as not evaluated and not compared.
```

```
$ sed -n '675p;678p;757p;759p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    geometry = vector_dispersion(c0,r)
    beta,arr,aphi = (geometry[k] for k in ('beta','arr','aphi'))
    aphij=evaluate_jet(aphi,mapping).tidy()
    optical_r=(aphij**sp.Rational(1,2)*c0).tidy()
```

```
$ sed -n '205p;232p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
(* COMPACT_K11_BOUNDARY *)
optical = Simplify[(c0^2 opticalMetric) /. localRules, $Assumptions];
```

```
$ grep -c '' research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round8_k11_evaluated_tags.txt research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round8_k11_not_evaluated_tags.txt
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round8_k11_evaluated_tags.txt:10
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round8_k11_not_evaluated_tags.txt:364
```

```
$ grep -n 'A_NONRECIPROCAL_PATH_DEPENDENCE$' research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round8_k11_not_evaluated_tags.txt
59:WL_S9B_A_NONRECIPROCAL_PATH_DEPENDENCE
```

```
$ grep 'evaluated 3 not_evaluated' _scratch/s9b_build/s9b_amend12_review_r5_claude_evidence/lookups.txt
evaluated 3 not_evaluated 376
```

```
$ grep -A2 'KEY aphi\|KEY beta' _scratch/s9b_build/s9b_amend12_review_r5_claude_evidence/k11_reach.stdout | grep 'KEY\|HAS_AZ'
KEY beta
  BASELINE_HAS_AZ False
  MUTANT_HAS_AZ True
KEY aphi
  BASELINE_HAS_AZ False
  MUTANT_HAS_AZ True
```

```
$ grep 'LITERAL_EQUAL' _scratch/s9b_build/s9b_amend12_review_r5_claude_evidence/wl_k11_fermat_pair.stdout
LITERAL_EQUAL False
```

## Grok: the K10 conjunction and the sign under the round-5 rule
```
$ cat _scratch/s9b_build/s9b_amend12_review_r5_grok_evidence/k10_sign_coverage.stdout
SAMPLE_TRAVEL_FALSE False False
POINT_TRAVEL_TRUE True False
AGREEMENT_ON_EVERY_BOOLEAN_INPUT False
EXISTS_INPUT_WHERE_THEY_DIFFER {travel: True, gate: False}
EXISTS_UNMASKED_DIFFERENCE {travel: True, gate: False}
AGREEMENT_IF_TRAVEL_FORCED_FALSE True
RHO_AT_5_OVER_RHO0 -1
RHO_AT_5_NEGATIVE True
RHO_AT_100_OVER_RHO0 1 - 2*exp(-9025)
RHO_AT_100_POSITIVE True
LIMIT_MINUS_RHO0 0
POINT_MEETS_SIGN_AND_FAR_POINT_DOES_NOT True
```

