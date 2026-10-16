# Measurements — amendment 12 review, round 0 (generated 2026-10-16 16:48)

Generator: `_scratch/s9b_build/gen/s9b_amend12_r0_lookups.sh` (sha256sum/grep/sed/cat/git show only). Line
numbers in `_scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r0.md` are the reviewed round-0 text's (saved before launch); the engines' and harness's are the
committed files' at `86097a44`.

## The reviewed text and the baseline
```
$ grep 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend12_review_baseline_r0.sha256 | cut -c1-64
f3cd38ed584857412dfbe3e71773e454da94bed45c2d49f4fd0ed07cfb0b7d46
```

```
$ sha256sum _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r0.md | cut -c1-64
f3cd38ed584857412dfbe3e71773e454da94bed45c2d49f4fd0ed07cfb0b7d46
```

```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_directive_amend12_review_baseline_r0.sha256 | grep -c ': OK$'
9
```

```
$ sha256sum _scratch/s9b_build/s9b_repair_build_directive_amend12_review_prompt_r0.md | cut -c1-64
5570a0b8a5cb2e24ac88da5fb3ad8f59719c0ee5f2ca0824df182a92407eae3d
```

## Verdicts
```
$ head -1 _scratch/s9b_build/s9b_amend12_review_r0_claude.md
**Verdict: NEEDS REVISION**
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend12_review_r0_grok.txt
Verdict: NEEDS REVISION
```

## The new search paragraph
```
$ grep -n 'look for one witness among its site parts' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r0.md
285:    - For each knife whose construction the engine contains, look for one witness among its site parts.
```

```
$ grep -n 'take the union of locations in the two copies' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r0.md
294:      parts. Trace direct and transitive computational uses, and take the union of locations in the two copies
```

```
$ grep -n 'after evaluating every site part in both copies at test inputs' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r0.md
296:      after evaluating every site part in both copies at test inputs, with its attached conditions. For **none found**,
```

```
$ grep -n 'treating each' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r0.md
292:      Enumerate all such values and domains recursively through payload containers and cases, treating each
```

```
$ grep -n 'both item 9 density-sign' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r0.md
293:      mathematical expression or predicate as one part. Include engine-local tags and both item 9 density-sign
```

## One-copy locations and positional keys (Grok 1, Claude 1)
```
$ sed -n '124p;127p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
        self.c = {g: v for g, v in self.c.items() if v != 0}
        return self.c.get(g, sp.S.Zero)
```

```
$ grep -n 'positions whose grade key differs' _scratch/s9b_build/s9b_amend12_review_r0_claude_evidence/out/py_jet_keys.stdout | head -3
2:  K1: 4 entries; positions whose grade key differs from the baseline: [(2, (0, 2, 0), (1, 0, 0)), (3, (0, 2, 1), (1, 0, 1))]
3:  K2: 4 entries; positions whose grade key differs from the baseline: [(1, (0, 0, 1), (0, 2, 0)), (2, (0, 2, 0), (1, 0, 0)), (3, (0, 2, 1), (1, 2, 0))]
4:  K3: 4 entries; positions whose grade key differs from the baseline: []
```

```
$ grep -n 'tags with location sets differing between copies' _scratch/s9b_build/s9b_amend12_review_r0_claude_evidence/out/py_location_skeleton.stdout | head -6
1:K1 R1: union expression leaves (excl. export roundtrip/digests) = 2544; tags with location sets differing between copies = 14
7:K1 R2: union expression leaves (excl. export roundtrip/digests) = 10380; tags with location sets differing between copies = 112
13:K2 R1: union expression leaves (excl. export roundtrip/digests) = 2504; tags with location sets differing between copies = 3
17:K2 R2: union expression leaves (excl. export roundtrip/digests) = 10214; tags with location sets differing between copies = 35
22:K3 R1: union expression leaves (excl. export roundtrip/digests) = 2504; tags with location sets differing between copies = 2
25:K3 R2: union expression leaves (excl. export roundtrip/digests) = 10216; tags with location sets differing between copies = 38
```

```
$ grep -n 'tags whose leaf-location sets differ between copies' _scratch/s9b_build/s9b_amend12_review_r0_claude_evidence/out/wl_location_skeleton.stdout | head -6
3:K1: tags 374; union leaf locations 6995; tags whose leaf-location sets differ between copies 20
7:K2: tags 374; union leaf locations 6853; tags whose leaf-location sets differ between copies 16
11:K3: tags 374; union leaf locations 7023; tags whose leaf-location sets differ between copies 30
15:K4a: tags 374; union leaf locations 6853; tags whose leaf-location sets differ between copies 12
19:K4b: tags 374; union leaf locations 6853; tags whose leaf-location sets differ between copies 20
23:K4c: tags 374; union leaf locations 6853; tags whose leaf-location sets differ between copies 16
```

```
$ grep -c '' _scratch/s9b_build/s9b_amend12_review_r0_grok_evidence/site_locations2.stdout
194
```

```
$ head -30 _scratch/s9b_build/s9b_amend12_review_r0_grok_evidence/site_locations2.stdout
DIRECTIVE_RULE
Enumerate all such values and domains recursively through payload containers and cases, treating each
      mathematical expression or predicate as one part. Include engine-local tags and both item 9 density-sign
      parts. Trace direct and transitive computational uses, and take the union of locations in the two copies
      (the compact pair for K11), without filtering by structural marks. Print one witness, or **none found** only
      after evaluating every site part in both copies at test inputs, with its attached conditions. For **none found**,
      name each part tried by tag, payload location and value or domain, and print its inputs and the evaluations
      actually performed.

SYMPY_JET_LINES
  124: self.c = {g: v for g, v in self.c.items() if v != 0}
  127: return self.c.get(g, sp.S.Zero)

OPTICAL_GRADE_KEYS
  baseline count=8 bytes=19896
    Tuple(Integer(0), Integer(0), Integer(0))
    Tuple(Integer(0), Integer(0), Integer(1))
    Tuple(Integer(0), Integer(2), Integer(0))
    Tuple(Integer(0), Integer(2), Integer(1))
    Tuple(Integer(1), Integer(0), Integer(0))
    Tuple(Integer(1), Integer(0), Integer(1))
    Tuple(Integer(1), Integer(2), Integer(0))
    Tuple(Integer(1), Integer(2), Integer(1))
  K1 count=4 bytes=2047
    Tuple(Integer(0), Integer(0), Integer(0))
    Tuple(Integer(0), Integer(0), Integer(1))
    Tuple(Integer(1), Integer(0), Integer(0))
    Tuple(Integer(1), Integer(0), Integer(1))
    only_baseline ['Tuple(Integer(0), Integer(2), Integer(0))', 'Tuple(Integer(0), Integer(2), Integer(1))', 'Tuple(Integer(1), Integer(2), Integer(0))', 'Tuple(Integer(1), Integer(2), Integer(1))']
    only_this []
```

## The unit of a part (Claude 1c)
```
$ sed -n '175p;178p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
    'K5':('PY_LOCAL_S9B_NONRECIPROCITY_ONEFORM','matrix00','path=nonreciprocity(r,c0,delta,velocity,xi,return_orientation)'),
    'K8':('PY_LOCAL_S9B_MASS_BALANCE','rhs','exchange=sp.simplify(-sp.diff(area*mass_density*velocity,r)/area)'),
```

```
$ sed -n '200,202p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
    if selector=='matrix00':return whole[0,0],(2,0)
    if selector=='whole':return whole,()
    if selector=='rhs':return whole.rhs,(1,)
```

## No floor on a try (Claude 2)
```
$ cat _scratch/s9b_build/s9b_amend12_review_r0_claude_evidence/out/py_degenerate_try.stdout | head -12
baseline rhs: -rho_br_profile(s9b_r)*Derivative(s9b_V(s9b_r), s9b_r) - s9b_V(s9b_r)*Derivative(rho_br_profile(s9b_r), s9b_r) - 2*rho_br_profile(s9b_r)*s9b_V(s9b_r)/s9b_r
K8 copy rhs : s9b_constant_density*(-s9b_r*Derivative(s9b_V(s9b_r), s9b_r) - 2*s9b_V(s9b_r))/s9b_r
r8 harness profiles: point {'s9b_r': 2, 's9b_constant_density': 1}; profiles {'s9b_V': 1/(8*q**2), 'rho_br_profile': 1 + 1/q}; baseline 1/128; K8 copy 0; difference 1/128
V = 0 (other profiles as r8): point {'s9b_r': 2, 's9b_constant_density': 1}; profiles {'s9b_V': 0, 'rho_br_profile': 1 + 1/q}; baseline 0; K8 copy 0; difference 0
rho_br = 1 (other profiles as r8): point {'s9b_r': 2, 's9b_constant_density': 1}; profiles {'s9b_V': 1/(8*q**2), 'rho_br_profile': 1}; baseline 0; K8 copy 0; difference 0
```

```
$ sed -n '461,469p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
        candidate=part.replace(lambda x:isinstance(x,AppliedUndef) and x.func.__name__ in PROFILES,
            lambda x:PROFILES[x.func.__name__].subs(RADIUS,x.args[0]))
        free=whole.free_symbols
        ranged={symbol for symbol in free if knife=='K7' and symbol.name=='s9b_r'}
        assignments={symbol:sp.sympify(POINT[symbol.name]) for symbol in free-ranged if symbol.name in POINT}
        prepared=candidate.subs(assignments,simultaneous=True)
        # A literal equality decision is sufficient only for this candidate's
        # local equation. It says nothing about its still quantified gates.
        decided=sp.true if prepared is sp.true else sp.false if prepared is sp.false else Str('UNRESOLVED')
```

```
$ grep -o "Str('partial_evaluation'), Str('UNRESOLVED'" _scratch/s9b_build/s9b_amend12_review_r0_claude_evidence/out/py_r8_K7_witness_line.txt _scratch/s9b_build/s9b_amend12_review_r0_claude_evidence/out/py_r8_K10_witness_line.txt | sort | uniq -c
      2 _scratch/s9b_build/s9b_amend12_review_r0_claude_evidence/out/py_r8_K10_witness_line.txt:Str('partial_evaluation'), Str('UNRESOLVED'
      2 _scratch/s9b_build/s9b_amend12_review_r0_claude_evidence/out/py_r8_K7_witness_line.txt:Str('partial_evaluation'), Str('UNRESOLVED'
```

## The density-sign parts (Claude 3)
```
$ sed -n '86,93p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
def density_parts(value, density):
    """One object with exhaustive sign labels, not a positivity gate."""
    if density is None:
        return cas(value)
    return cas(('LOCAL_BULK_DENSITY_PARTS',
                ('RHO_POSITIVE',sp.Gt(density,0)),
                ('RHO_NONPOSITIVE',sp.Le(density,0)),
                ('OBJECT_ON_BOTH_PARTS',value)))
```

```
$ sed -n '618,623p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
densityParts[i_, object_] := If[i == 1, object, <|
  "SuppliedBulkDensityPremise" -> (rho0 > 0),
  "DensityParts" -> {<|"Label" -> "rho > 0", "LabelCondition" -> (rho0 (1+fProfile[r]) > 0),
      "Object" -> object|>,
    <|"Label" -> "rho <= 0", "LabelCondition" -> (rho0 (1+fProfile[r]) <= 0),
      "Object" -> object|>}|>];
```

```
$ grep -n 'they add no attached condition' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r0.md
303:      shared object; they add no attached condition. A witness may lie on either part, including `ρ ≤ 0`; print
```

## "the guard" (Claude 4)
```
$ grep -n 'guard' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r0.md
73:part before **none found**. Part 2 also adds fixture marks or values before software-check guards.
248:    - Run every CAS demonstration through `scripts/s11c_guarded_run.py` with `--pool s9b --memory-gib 16`. Wolfram
251:    - If the guard refuses admission, kills a run, or a kernel hits the task limit, stop and report it. A higher
370:    **K11 evaluation (amendment 2; the user's choice, 2026-10-08).** The guard killed K11's full copy at 8 GiB and
373:    guard kills the compact run.
401:    - every guard refusal or kill;
443:  each software-check fixture's computed mark or value before the guard.
```

