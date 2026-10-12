# Measurements — S9b repair build, review round 4 (generated 2026-10-11 20:07)

Generator: `_scratch/s9b_build/gen/s9b_build_r4_lookups.sh` (sha256sum/grep/sed/cat/wc/ls only). `_scratch/s9b_build/s9b_repair_build_review_r4_claude_evidence/construct/` holds
the Claude leg's scripts and stdout, copied verbatim from `/tmp/s9b_r4_claudeleg/construct/`; `_scratch/s9b_build/s9b_repair_build_review_r4_grok_evidence/cas1/` holds
Grok's, copied verbatim from `/tmp/s9b_r4_grok_cas1/` (files under 2 MB).

## The reviewed artifacts are the working-tree files, unchanged after both legs
```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_review_baseline_r4.sha256
research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py: OK
research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py: OK
research/pde_ledger_v3/scripts/S9b_exports.py: OK
research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl: OK
research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_build_report.md: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_tag_inventory.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_build_report.md: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_tags.txt: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_k11_not_evaluated_tags.txt: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_k11_evaluated_tags.txt: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_not_established_tags.txt: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_executable_checks.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_payload_checks.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_validation.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round4_measurements.md: OK
research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md: OK
research/pde_ledger_v3/directives/S9b_repair_build_directive.md: OK
```

```
$ sha256sum _scratch/s9b_build/s9b_repair_build_review_prompt_r4.md
9c720533278d4ff004d571801c681fd3fa823da2b6363ea245eaab53dc41c846  _scratch/s9b_build/s9b_repair_build_review_prompt_r4.md
```

## Both verdicts
```
$ grep -n -o 'Verdict[*:]* *[A-Z][A-Z ]*' _scratch/s9b_build/s9b_repair_build_review_r4_claude.md _scratch/s9b_build/s9b_repair_build_review_r4_grok.txt
_scratch/s9b_build/s9b_repair_build_review_r4_claude.md:1:Verdict: NEEDS REVISION 
_scratch/s9b_build/s9b_repair_build_review_r4_grok.txt:1:Verdict: NEEDS REVISION
```

## The Grok leg ran once, to its report; it waited for the shared pool
```
$ cat _scratch/s9b_build/build_r4_grok.DONE
EXIT=0
```

```
$ ls _scratch/s11c/s9b-review-r4 | grep -c '^grok-probe-physics-'
830
```

```
$ ls _scratch/s11c/s9b-review-r4/grok-probe-physics-500 | wc -l
0
```

```
$ ls _scratch/s11c/s9b-review-r4/grok-probe-physics-830
child-outcome.json
duration-validation.json
effective-limits.json
invocation.json
limit-validation.json
outcome.json
resource-samples.jsonl
stderr
stdout
```

```
$ grep -n 'wallSeconds' _scratch/s11c/s9b-review-r4/claudeleg-py-harness/outcome.json
3:  "wallSeconds": 25434.276952803135,
10:    "wallSeconds": 25434.187951609027,
```

## Finding G1 (Grok 1): the engines attach different domains to the fixed-ratio and power rows of Part C
### What the spec supplies, and what directive item 9 allows
```
$ sed -n 218,224p research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
- **Part C bulk (supplied):**
  - `P = Kρ^n`, with `n` symbolic, `ρ` the bulk number density and `m` the particle mass (v3 convention). So
    `c_s² = nKρ^(n−1)/m`.
  - Define the fractional bulk-density change at the brane by
    `f(x) ≡ δρ(x)/ρ₀ = ρ(x)/ρ₀ − 1`. It is a symbolic radial profile owned by the gravity sector or S12.
    `ρ₀` is the asymptotic bulk number density, `c_s0` is the corresponding asymptotic sound speed, and Part C
    keeps first order in `f`.
```

```
$ sed -n 248,252p research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
  Here `ρ(x)` is the bulk number density at the brane and `ρ₀` is its asymptotic value. For each condition,
  print the `j_n` it implies through the supplied mass balance, with `ρ_br` symbolic. Print where `n` enters,
  if it enters at all.

- **Part D.** The Part B conditions re-expressed with the linked neutral-sector steady in-plane
```

```
$ grep -n 'Each predicate follows from the' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
156:   not built from another row, or from the amplitude before the response is applied. Each predicate follows from the
```

```
$ sed -n 155,158p research/pde_ledger_v3/directives/S9b_repair_build_directive.md
9. **Part C domains (D6 B6).** Each Part C row's domain predicates are built from that row's own response. They are
   not built from another row, or from the amplitude before the response is applied. Each predicate follows from the
   row's response being real and positive, or from the row's supplied inputs. The first-order counting in `f` is an
   order, not a numeric bound on any amplitude.
```

### SymPy adds the bulk-density conjunct to both non-constant rows; Wolfram's domain does not carry it
```
$ sed -n 994,997p research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
        local=response.subs(f_symbol,f)
        domain=sp.And(sp.Contains(local,sp.S.Reals),sp.Gt(local,0))
        if label!='CONSTANT': domain=sp.And(domain,sp.Gt(rho0*(1+f),0))
        if label=='FIXED_RATIO': domain=sp.And(domain,sp.Gt(cs0_squared,0),sp.Gt(cs_squared.subs(rho_bulk,rho0*(1+f)),0))
```

```
$ sed -n 579,583p research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
responseDomain[i_] := With[{resp = responseInputs[[i]]},
  Element[resp,Reals] && resp > 0 &&
    If[i == 2, Element[{Sqrt[localSoundSquared],
      Sqrt[bulkSoundSquared /. rho -> rho0]},Reals] &&
      (bulkSoundSquared /. rho -> rho0) > 0, True]] /. fLocal -> fProfile[r];
```

```
$ sed -n 46,48p research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
  GM, n, K, m, s, rho0, Phi, rAnchor}, Reals] && c0 > 0 && ell > 0 &&
  b > bFar > 0 && ZE > 0 && ZR > 0 && r > 0 && rr > 0 && impact > 0 &&
  rho0 > 0 && m > 0 && rAnchor > bFar;
```

```
$ sed -n 573,574p research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
responseDeltas = (Simplify[Normal[Series[#/c0 - 1,{fLocal,0,1}]],
  rho0 > 0 && c0 > 0 && Element[{n,s,fLocal},Reals]] &) /@ responseInputs;
```

### Grok's probes of the two engines (one guarded run)
```
$ sed -n 44,45p _scratch/s9b_build/s9b_repair_build_review_r4_grok_evidence/cas1/probe_jn_gamma.stdout
PROBE C_POWER_DOMAIN Contains((s9b_f(s9b_r) + 1)**s, Reals) & (rho_0*(s9b_f(s9b_r) + 1) > 0) & ((s9b_f(s9b_r) + 1)**s > 0)
PROBE POWER_AT_S2_F_MINUS2 1
```

```
$ sed -n 41,42p _scratch/s9b_build/s9b_repair_build_review_r4_grok_evidence/cas1/probe_follow.stdout
PROBE POWER_DOMAIN Element[c0*(1 + fProfile[r])^s, Reals] && c0*(1 + fProfile[r])^s > 0
PROBE POWER_DOMAIN_S2_FM2 True
```

```
$ sed -n 616p _scratch/s9b_build/s9b_repair_build_review_r4_grok_evidence/cas1/probe_follow.wl
probe["POWER_DOMAIN_S2_FM2", Simplify[responseDomain[3] /. {s -> 2, fProfile[r] -> -2}, c0 > 0 && rho0 > 0]];
```

### The conjunct predates the repairs (r1 baseline f4dc2d3d, r3 baseline c5c8fcc1)
```
$ git show f4dc2d3d:research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py | grep -n 'Gt(rho0\*(1+f),0)'
657:        if label!='CONSTANT': domain=sp.And(domain,sp.Gt(rho0*(1+f),0))
```

```
$ git show c5c8fcc1:research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py | grep -n 'Gt(rho0\*(1+f),0)'
932:        if label!='CONSTANT': domain=sp.And(domain,sp.Gt(rho0*(1+f),0))
```

## Findings C1 and C2 (Claude 1 and 2): unchanged parts print as nonzero differences
### Directive item 14
```
$ grep -n 'so an unchanged part prints as zero' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
240:      condition. Each part's difference is computed by the engine, so an unchanged part prints as zero. This holds
```

### The harnesses' difference functions
```
$ grep -n -A12 '^def _value_difference' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
242:def _value_difference(a,b):
243-    if a == b:
244-        return sp.S.Zero
245-    if isinstance(a,sp.Tuple) and isinstance(b,sp.Tuple) and len(a)==len(b):
246-        parts=sp.Tuple(*(_value_difference(x,y) for x,y in zip(a,b)))
247-        return sp.S.Zero if all(part==0 for part in parts) else parts
248-    if isinstance(a,sp.MatrixBase) and isinstance(b,sp.MatrixBase) and a.shape==b.shape:
249-        return b-a
250-    if isinstance(a,sp.Expr) and isinstance(b,sp.Expr):
251-        return b-a
252-    if isinstance(a,sp.logic.boolalg.Boolean) and isinstance(b,sp.logic.boolalg.Boolean):
253-        # Normalize arithmetic relations without evaluating functionals,
254-        # then compare their Boolean structure. This also removes spurious
```

```
$ grep -n -B4 -A2 'its membership problem while comparing' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
226-    replacements={}
227-    for relation in value.atoms(sp.core.relational.Relational):
228-        if not all(isinstance(side,sp.Expr) for side in (relation.lhs,relation.rhs)):
229-            # A quantified set relation is itself a predicate. Do not solve
230:            # its membership problem while comparing payload domains.
231-            continue
232-        sides=[side.replace(lambda node:isinstance(node,sp.Integral),_primitive_integral)
```

```
$ grep -n 'No assumption, sampling, global' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
66:   relations/predicates used as values. No assumption, sampling, global
```

```
$ grep -n 'True,a-b\|Boole\[a\]-Boole\[b\]' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
74:  logicalObjectQ[a] && logicalObjectQ[b], Boole[a]-Boole[b],
80:  True,a-b];
85:    (logicalObjectQ[b] || Head[b] === Symbol), Boole[a]-Boole[b],
```

### Which of these predate round 4
```
$ git show c5c8fcc1:research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py | grep -n -A10 '^def _value_difference' | grep -c 'return b-a'
2
```

```
$ git show c5c8fcc1:research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py | grep -c 'its membership problem while comparing'
1
```

```
$ git show c5c8fcc1:research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl | grep -n 'True,a-b\|Boole\[a\]-Boole\[b\]'
74:  logicalObjectQ[a] && logicalObjectQ[b], Boole[a]-Boole[b],
80:  True,a-b];
85:    (logicalObjectQ[b] || Head[b] === Symbol), Boole[a]-Boole[b],
```

### The legs' counts (no-op rows: a knife whose construction is absent from the row)
```
$ cat _scratch/s9b_build/s9b_repair_build_review_r4_claude_evidence/construct/c16_py_noop_rows.stdout
K1 value_only 0
K1 value_and_domain 2
K1 domain_only 56
K2 value_only 41
K2 value_and_domain 0
K2 domain_only 0
K3 value_only 1
K3 value_and_domain 37
K3 domain_only 84
K1 no-op rows 58 printed nonzero 58
K2 no-op rows 132 printed nonzero 41
K3 no-op rows 124 printed nonzero 122
```

```
$ head -3 _scratch/s9b_build/s9b_repair_build_review_r4_claude_evidence/construct/c15_wl_noop_rows.stdout
K1 {'VALUE_ZERO/DOMAIN_NONZERO': 28, 'VALUE_ZERO_INTEGRAL_ONLY/DOMAIN_NONZERO': 5, 'VALUE_ZERO/DOMAIN_ZERO': 20, 'VALUE_ZERO_INTEGRAL_ONLY/DOMAIN_ZERO': 5}
K2 {'VALUE_ZERO/DOMAIN_ZERO': 99, 'VALUE_ZERO_INTEGRAL_ONLY/DOMAIN_ZERO': 21, 'VALUE_NONZERO/DOMAIN_ZERO': 12}
K3 {'VALUE_ZERO/DOMAIN_NONZERO': 88, 'VALUE_ZERO_INTEGRAL_ONLY/DOMAIN_NONZERO': 21, 'VALUE_NONZERO/DOMAIN_NONZERO': 12, 'VALUE_NONZERO/DOMAIN_ZERO': 3}
```

```
$ grep -c 'differing_leaves_after_zero_functional_removal=0' _scratch/s9b_build/s9b_repair_build_review_r4_claude_evidence/construct/c10_K2_noop_rows.stdout
38
```

```
$ grep -c 'differing_leaves_after_zero_functional_removal=1' _scratch/s9b_build/s9b_repair_build_review_r4_claude_evidence/construct/c10_K2_noop_rows.stdout
3
```

```
$ grep -c 'simplified_residual=\[0\]' _scratch/s9b_build/s9b_repair_build_review_r4_claude_evidence/construct/c10_K2_noop_rows.stdout
3
```

```
$ cat _scratch/s9b_build/s9b_repair_build_review_r4_claude_evidence/construct/c14b_wl_conjunct_count.stdout
WL_LOCAL_S9B_ABLATION_K3_WL_S9B_F_FLOW_A_DEFLECTION_G000 Baseline count 2 Corrupted count 2
```

### Grok's harness run (the same two harnesses, the same baseline)
```
$ grep -n 'Every K1–K10 triple set has a nonzero value or a nonzero domain' _scratch/s9b_build/s9b_repair_build_review_r4_grok.txt | cut -c1-120
45:Harnesses, guard `s11c-guard-69b9ec333d86`, log `grok-harness-1`: Wolfram exit 0, SymPy exit 0. Baseline is the uncha
```

## Guard state after both legs
```
$ cat _scratch/s11c/resource-pool/reservations.json
{}
```

```
$ systemctl --user list-units 's11c-guard-*' --no-legend | wc -l
0
```

