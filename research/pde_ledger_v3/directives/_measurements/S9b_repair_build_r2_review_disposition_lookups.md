# Measurements — S9b repair build, review round 2 (generated 2026-10-10 10:14)

Generator: `_scratch/s9b_build/gen/s9b_build_r2_lookups.sh` (sha256sum/grep/sed/cat/wc only). Paths under
`_scratch/s9b_build/s9b_repair_build_review_r2_claude_evidence/` are the Claude leg's evidence, copied verbatim from `/tmp/s9b_r2_leg_claude/`.

## The reviewed artifacts are the working-tree files, unchanged after both legs
```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_review_baseline_r2.sha256
research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py: OK
research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py: OK
research/pde_ledger_v3/scripts/S9b_exports.py: OK
research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl: OK
research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_build_report.md: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_tag_inventory.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_staged_candidates.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_build_report.md: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_tags.txt: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_k11_not_evaluated_tags.txt: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_not_established_tags.txt: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_executable_checks.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_k11_checks.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_payload_checks.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round2_measurements.md: OK
research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md: OK
research/pde_ledger_v3/directives/S9b_repair_build_directive.md: OK
```

## Both verdicts
```
$ grep -n -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_repair_build_review_r2_claude.md _scratch/s9b_build/s9b_repair_build_review_r2_grok.txt
_scratch/s9b_build/s9b_repair_build_review_r2_claude.md:11:Verdict: NEEDS REVISION
_scratch/s9b_build/s9b_repair_build_review_r2_grok.txt:5:Verdict: CLEAR
```

### Grok's statements bearing on the two findings (report L28, L38, L42)
```
$ sed -n 28p _scratch/s9b_build/s9b_repair_build_review_r2_grok.txt | fold -s -w 300
SymPy emits the operator ratio `b/c₀`. Both engines emit the γ difference as 0 on GM ≠ 0, and `NOT_ESTABLISHED` off the ray gate. Wolfram’s endpoint limit is `<|"100" -> 0, "010" -> 0, "020" -> 0, "001" -> 0|>`. SymPy’s counting-certificate limits are `Integer(0)` on the same four grades, 
so the dropped endpoint remainder vanishes.
```

```
$ sed -n 38p _scratch/s9b_build/s9b_repair_build_review_r2_grok.txt | fold -s -w 300
5. **Harnesses.** Each knife is one replacement at its named site. The baseline is the unchanged engine. A reader can see a zero: SymPy prints `Integer(0)`; Wolfram’s dead-path self-diff is 0 on all 373 tags. Nonzero counts: SymPy K1 335, K2 247, K3 335, K4a 234, K4b 270, K4c 246, K5 35, K6 77, 
K7 51, K8 167, K9a 15, K9b 16, K10 320. Wolfram K1 316, K2 237, K3 337, K4a 104, K4b 104, K4c 234, K5 38, K6 54, K7 47, K8 184, K9a 13, K9b 12, K10 324. The bites sit on the construction each knife names (K5 on the odd nonreciprocal grades, K6 on the radar tags, K7 on the every-b conditions, 
K9a/K9b on the fixed-ratio and power rows). Stdout: `/var/projects/toy_physics/_scratch/s11c/s9b-review-r2/grok-py-harness-1/stdout` and `grok-wl-harness-1/stdout`.
```

```
$ sed -n 42p _scratch/s9b_build/s9b_repair_build_review_r2_grok.txt | fold -s -w 300
7. **Tags, export, repair text, leakage.** The two engines emit the same 330 shared names, and the objects checked above match under those names. SymPy’s fold binds `c_s0` to the upstream positive-symbol knob and anchors it with the bulk sound speed. No fold row stores a live profile. No 
`assert`, `VERDICT`, `PASS`, or `FAIL` appears in either engine. Branch type is computed from the signed c_γ² = μ_⊥/ρ radicand. Restrictions substitute into the ray gates. Wolfram keeps the GM = 0 stratum on the radar γ.
```

## The Claude leg's addendum: numeric ground truth for Part A (its guarded runs)
```
$ cat _scratch/s11c/s9b-review-r2/claude-num2-1/outcome.json | grep -E 'exitCode|unit'
  "exitCode": 0,
  "unit": "s11c-guard-5d7d9158ef84",
    "exitCode": 0,
```

```
$ grep -E '^CMP_(ROUND_TRIP|ONE_WAY_ER|DEFLECTION|NONRECIPROCAL)[^:]*_G(101|021|120|121):' _scratch/s11c/s9b-review-r2/claude-cmpnum-1/stdout | cut -c1-260
CMP_DEFLECTION_G021: NUM=0.017995734723637172341 | PY=0.017995734723264685991 | WL=0.01799573472326465034497802460035732491`20. | PY-NUM=-3.72486e-13 | WL-NUM=-3.72522e-13
CMP_DEFLECTION_G101: NUM=0.0054804344382852413078 | PY=0.0054804344382758510755 | WL=0.00548043443827584039518556884222409532`20. | PY-NUM=-9.39023e-15 | WL-NUM=-9.40091e-15
CMP_DEFLECTION_G120: NUM=0.027316769362330711839 | PY=0.027316769361656216291 | WL=0.0273167693616561638054258847402577794`20. | PY-NUM=-6.74496e-13 | WL-NUM=-6.74548e-13
CMP_DEFLECTION_G121: NUM=0.0047799256093445207135 | PY=0.004779925609146979522 | WL=0.0047799256091469684102388510032834018`20. | PY-NUM=-1.97541e-13 | WL-NUM=-1.97552e-13
CMP_NONRECIPROCAL_G021: NUM=0.0 | PY=0.0 | WL=0 | PY-NUM=0.0 | WL-NUM=0.0
CMP_NONRECIPROCAL_G101: NUM=0.0 | PY=0.0 | WL=0 | PY-NUM=0.0 | WL-NUM=0.0
CMP_NONRECIPROCAL_G120: NUM=0.0 | PY=0.0 | WL=0 | PY-NUM=0.0 | WL-NUM=0.0
CMP_NONRECIPROCAL_G121: NUM=0.0 | PY=0.0 | WL=0 | PY-NUM=0.0 | WL-NUM=0.0
CMP_ONE_WAY_ER_G021: NUM=-1.2016309412209786955 | PY=-1.2016309412816944316 | WL=-1.2016309412816950593911714064603702683`20. | PY-NUM=-6.07157e-11 | WL-NUM=-6.07164e-11
CMP_ONE_WAY_ER_G101: NUM=-0.9596998802256385997 | PY=-0.95969988009307371918 | WL=-0.95969988009307410862342600363530923521`20. | PY-NUM=1.32565e-10 | WL-NUM=1.32564e-10
CMP_ONE_WAY_ER_G120: NUM=-1.3520910833444936846 | PY=-1.3520910834077978869 | WL=-1.35209108340779862126867245519820399207`20. | PY-NUM=-6.33042e-11 | WL-NUM=-6.33049e-11
CMP_ONE_WAY_ER_G121: NUM=1.9727538941599717033 | PY=1.9727538937122121191 | WL=1.97275389371221284554815326729610205812`20. | PY-NUM=-4.4776e-10 | WL-NUM=-4.47759e-10
CMP_ROUND_TRIP_G021: NUM=-2.403261882441957391 | PY=-2.4032618825633888632 | WL=-2.40326188256339011878234281292074053659`20. | PY-NUM=-1.21431e-10 | WL-NUM=-1.21433e-10
CMP_ROUND_TRIP_G101: NUM=-1.9193997604512771994 | PY=-1.9193997601861474384 | WL=-1.91939976018614821724685200727061847041`20. | PY-NUM=2.6513e-10 | WL-NUM=2.65129e-10
CMP_ROUND_TRIP_G120: NUM=-2.7041821666889873691 | PY=-2.7041821668155957738 | WL=-2.70418216681559724253734491039640798415`20. | PY-NUM=-1.26608e-10 | WL-NUM=-1.2661e-10
CMP_ROUND_TRIP_G121: NUM=3.9455077883199434067 | PY=3.9455077874244242381 | WL=3.94550778742442569109630653459220411625`20. | PY-NUM=-8.95519e-10 | WL-NUM=-8.95518e-10
```

```
$ grep -E '^BORN_ROUND_TRIP[^:]*_G(101|021|120|121):' _scratch/s11c/s9b-review-r2/claude-born-1/stdout
BORN_ROUND_TRIP_b=40.0_ZE=300.0_ZR=500.0_G021: 0.75762045225313637586
BORN_ROUND_TRIP_b=40.0_ZE=300.0_ZR=500.0_G101: 0.10835461339668433076
BORN_ROUND_TRIP_b=40.0_ZE=300.0_ZR=500.0_G120: 1.0270963162294076345
BORN_ROUND_TRIP_b=40.0_ZE=300.0_ZR=500.0_G121: 0.056983214341411400456
```

## The findings' code is already in the review-r1 baseline (f4dc2d3d), so the r1→r2 repair did not introduce it
```
$ git show f4dc2d3d:research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py | grep -n 'sp.Ne(GM,0)'
622:            emit(prefix+'B_GAMMA_'+name,gated(substitute(gammas[name],replacement),sp.And(domain,sp.Ne(GM,0))))
624:        emit(prefix+'B_GAMMA_DIFFERENCE',gated(substitute(gammas['DEFLECTION']-gammas['RADAR'],replacement),sp.And(domain,sp.Ne(GM,0))))
```

```
$ git show f4dc2d3d:research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl | grep -n '_ConditionalExpression'
47:    !FreeQ[{a,b}, _String | _Association | _Rule | _Missing | _ConditionalExpression |
```

## Claude finding 1: SymPy labels the GM = 0 stratum of the γ objects NOT_ESTABLISHED
### What the directive asks for
```
$ sed -n 112,116p research/pde_ledger_v3/directives/S9b_repair_build_directive.md
5. **No verdicts (D6 B5).** No `VERDICT`, `PASS` or `FAIL`. A boolean-valued test is emitted as the CAS object the
   test returned. Emission never depends on a payload's value. Outside the branch-existence conditions the spec
   requires `NOT_ESTABLISHED`. That output, and the branch type, are produced from the computed conditions, never
   typed. The branch-existence condition is computed through `c_γ²` without presupposing its sign (spec, "Branch
   existence"), so each branch type the spec lists can be produced from it.
```

```
$ sed -n 129,130p research/pde_ledger_v3/directives/S9b_repair_build_directive.md
   - The `γ` difference is reduced by the engine's CAS, so its printed form shows its value on each stratum. It is
     not printed as an uncombined difference of two objects.
```

### SymPy: the γ tags are gated on (domain ∧ GM ≠ 0), and the complement is printed NOT_ESTABLISHED
```
$ sed -n 505,506p research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
def gated(value, domain):
    return cas(((domain,value),(sp.Not(domain),'NOT_ESTABLISHED')))
```

```
$ sed -n 820,829p research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
        for name in comparisons:
            progress(prefix+'B_GAMMA_'+name)
            emit(prefix+'B_GAMMA_'+name,gated(substitute(gammas[name],replacement),sp.And(domain,sp.Ne(GM,0))))
            progress(prefix+'B_RESIDUAL_'+name)
            emit(prefix+'B_RESIDUAL_'+name,gated(substitute(residuals[name],replacement),domain))
        progress(prefix+'B_GAMMA_DIFFERENCE: substitute and reduce integral kernels')
        difference=reduce_functional(substitute(gammas['DEFLECTION']-gammas['RADAR'],replacement))
        progress(prefix+'B_GAMMA_DIFFERENCE')
        emit(prefix+'B_GAMMA_DIFFERENCE',gated(difference,sp.And(domain,sp.Ne(GM,0))))
    def emit_a(prefix,replacement,domain,all_times=True):
```

### SymPy baseline print: the first slot's domain opens with GM ≠ 0; the second slot is NOT_ESTABLISHED
```
$ grep -m1 '^PY_S9B_B_GAMMA_DIFFERENCE:' _scratch/s9b_build/s9b_repair_build_review_r2_claude_evidence/repo/_scratch/py-harness-1/baseline/stdout | cut -c1-120
PY_S9B_B_GAMMA_DIFFERENCE: Tuple(Tuple(And(Unequality(Symbol('GM', real=True), Integer(0)), StrictLessThan(Symbol('s9b_r
```

```
$ grep -m1 '^PY_S9B_B_GAMMA_DIFFERENCE:' _scratch/s9b_build/s9b_repair_build_review_r2_claude_evidence/repo/_scratch/py-harness-1/baseline/stdout | rev | cut -c1-40 | rev
 Integer(0)))), Str('NOT_ESTABLISHED')))
```

### Wolfram baseline print: the γ difference carries a GM == 0 stratum; the γ tags are relations in γ with no GM gate
```
$ grep -m1 '^WL_S9B_B_GAMMA_DEFLECTION:' _scratch/s9b_build/s9b_repair_build_review_r2_claude_evidence/repo/_scratch/s9b_wl_repair2/ablation-final2/BASELINE.stdout | grep -o 'GM == 0' | wc -l
0
```

```
$ grep -m1 '^WL_S9B_B_GAMMA_DEFLECTION:' _scratch/s9b_build/s9b_repair_build_review_r2_claude_evidence/repo/_scratch/s9b_wl_repair2/ablation-final2/BASELINE.stdout | grep -o 'GM != 0' | wc -l
0
```

```
$ grep -m1 '^WL_S9B_B_GAMMA_DEFLECTION:' _scratch/s9b_build/s9b_repair_build_review_r2_claude_evidence/repo/_scratch/s9b_wl_repair2/ablation-final2/BASELINE.stdout | cut -c1-150
WL_S9B_B_GAMMA_DEFLECTION: <|"OnDomain" -> ConditionalExpression[((2*GM)/(b*c0^2) == Inactive[Integrate][(2*b*Derivative[1][deltaProfile][r])/Sqrt[-b^
```

```
$ grep -m1 '^WL_S9B_B_GAMMA_RADAR:' _scratch/s9b_build/s9b_repair_build_review_r2_claude_evidence/repo/_scratch/s9b_wl_repair2/ablation-final2/BASELINE.stdout | grep -o 'GM == 0' | wc -l
0
```

```
$ grep -m1 '^WL_S9B_B_GAMMA_RADAR:' _scratch/s9b_build/s9b_repair_build_review_r2_claude_evidence/repo/_scratch/s9b_wl_repair2/ablation-final2/BASELINE.stdout | grep -o 'GM != 0' | wc -l
0
```

```
$ grep -m1 '^WL_S9B_B_GAMMA_RADAR:' _scratch/s9b_build/s9b_repair_build_review_r2_claude_evidence/repo/_scratch/s9b_wl_repair2/ablation-final2/BASELINE.stdout | cut -c1-150
WL_S9B_B_GAMMA_RADAR: <|"OnDomain" -> ConditionalExpression[((2*GM)/c0^3 == Inactive[Integrate][(2*b^2*Derivative[1][deltaProfile][r])/(c0*Sqrt[-b^2 +
```

```
$ grep -m1 '^WL_S9B_B_GAMMA_DIFFERENCE:' _scratch/s9b_build/s9b_repair_build_review_r2_claude_evidence/repo/_scratch/s9b_wl_repair2/ablation-final2/BASELINE.stdout | grep -o 'GM == 0' | wc -l
1
```

```
$ grep -m1 '^WL_S9B_B_GAMMA_DIFFERENCE:' _scratch/s9b_build/s9b_repair_build_review_r2_claude_evidence/repo/_scratch/s9b_wl_repair2/ablation-final2/BASELINE.stdout | grep -o 'GM != 0' | wc -l
1
```

```
$ grep -m1 '^WL_S9B_B_GAMMA_DIFFERENCE:' _scratch/s9b_build/s9b_repair_build_review_r2_claude_evidence/repo/_scratch/s9b_wl_repair2/ablation-final2/BASELINE.stdout | cut -c1-150
WL_S9B_B_GAMMA_DIFFERENCE: <|"OnDomain" -> ConditionalExpression[(gammaDifference == 0 && GM != 0) || (GM == 0 && Inactive[Integrate][(b^2*(velocityPr
```

## Claude finding 2: the Wolfram harness prints a gated payload's difference unreduced
### What the directive asks for
```
$ sed -n 212,216p research/pde_ledger_v3/directives/S9b_repair_build_directive.md
14. **Ablation harness** (`docs/development_pipeline.md` §4).
    - Each engine gets a committed harness that runs the live engine unchanged as the baseline. For each knife
      below, it runs a copy with exactly that one mutation at the named construction site.
    - It prints the baseline payload, the corrupted payload and their difference for **every** tag, including tags
      whose difference is zero. It prints no verdict.
```

### Wolfram harness: a ConditionalExpression operand makes the whole difference Inactive[Subtract]
```
$ sed -n 47,53p research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
difference[a_, b_] := Which[SameQ[a,b], 0,
  StringQ[a] || StringQ[b] || MemberQ[{True, False}, a] || MemberQ[{True, False}, b] ||
    !FreeQ[{a,b}, _String | _Association | _Rule | _Missing | _ConditionalExpression |
      _Piecewise | _Equal | _Unequal | _Less | _LessEqual | _Greater | _GreaterEqual |
      _Inequality | _And | _Or | _Not | _Xor | _Element | _Exists | _ForAll],
    Inactive[Subtract][a,b],
  True, a-b];
```

### SymPy harness: boolean slots are differenced by Xor
```
$ sed -n 204,216p research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
def difference(a,b):
    if a == b:
        return sp.S.Zero
    if isinstance(a,sp.Tuple) and isinstance(b,sp.Tuple) and len(a)==len(b):
        return sp.Tuple(*(difference(x,y) for x,y in zip(a,b)))
    if isinstance(a,sp.MatrixBase) and isinstance(b,sp.MatrixBase) and a.shape==b.shape:
        return b-a
    if isinstance(a,sp.Expr) and isinstance(b,sp.Expr):
        return b-a
    if isinstance(a,sp.logic.boolalg.Boolean) and isinstance(b,sp.logic.boolalg.Boolean):
        return sp.Xor(a,b)
    return sp.sympify(a!=b)

```

### The leg's slot split of the Wolfram harness stream
```
$ grep -E '^WLSLOT_K[0-9a-z]+_(VALUE_CHANGED|DOMAIN_ONLY|UNCHANGED|NO_CE_SLOT):' _scratch/s9b_build/s9b_repair_build_review_r2_claude_evidence/wl_value_slots.stdout
WLSLOT_K1_VALUE_CHANGED: 122 (A-tags 29)
WLSLOT_K1_DOMAIN_ONLY: 132 (A-tags 35)
WLSLOT_K1_UNCHANGED: 16 (A-tags 0)
WLSLOT_K1_NO_CE_SLOT: 60 (A-tags 1)
WLSLOT_K2_VALUE_CHANGED: 99 (A-tags 23)
WLSLOT_K2_DOMAIN_ONLY: 93 (A-tags 41)
WLSLOT_K2_UNCHANGED: 78 (A-tags 0)
WLSLOT_K2_NO_CE_SLOT: 60 (A-tags 1)
WLSLOT_K3_VALUE_CHANGED: 99 (A-tags 23)
WLSLOT_K3_DOMAIN_ONLY: 171 (A-tags 41)
WLSLOT_K3_UNCHANGED: 0 (A-tags 0)
WLSLOT_K3_NO_CE_SLOT: 60 (A-tags 1)
WLSLOT_K4a_VALUE_CHANGED: 69 (A-tags 17)
WLSLOT_K4a_DOMAIN_ONLY: 0 (A-tags 0)
WLSLOT_K4a_UNCHANGED: 201 (A-tags 47)
WLSLOT_K4a_NO_CE_SLOT: 60 (A-tags 1)
WLSLOT_K4b_VALUE_CHANGED: 70 (A-tags 17)
WLSLOT_K4b_DOMAIN_ONLY: 0 (A-tags 0)
WLSLOT_K4b_UNCHANGED: 200 (A-tags 47)
WLSLOT_K4b_NO_CE_SLOT: 60 (A-tags 1)
WLSLOT_K4c_VALUE_CHANGED: 99 (A-tags 23)
WLSLOT_K4c_DOMAIN_ONLY: 93 (A-tags 41)
WLSLOT_K4c_UNCHANGED: 78 (A-tags 0)
WLSLOT_K4c_NO_CE_SLOT: 60 (A-tags 1)
WLSLOT_K5_VALUE_CHANGED: 36 (A-tags 12)
WLSLOT_K5_DOMAIN_ONLY: 0 (A-tags 0)
WLSLOT_K5_UNCHANGED: 234 (A-tags 52)
WLSLOT_K5_NO_CE_SLOT: 60 (A-tags 1)
WLSLOT_K6_VALUE_CHANGED: 28 (A-tags 3)
WLSLOT_K6_DOMAIN_ONLY: 0 (A-tags 0)
WLSLOT_K6_UNCHANGED: 242 (A-tags 61)
WLSLOT_K6_NO_CE_SLOT: 60 (A-tags 1)
WLSLOT_K7_VALUE_CHANGED: 0 (A-tags 0)
WLSLOT_K7_DOMAIN_ONLY: 0 (A-tags 0)
WLSLOT_K7_UNCHANGED: 270 (A-tags 64)
WLSLOT_K7_NO_CE_SLOT: 60 (A-tags 1)
WLSLOT_K8_VALUE_CHANGED: 56 (A-tags 0)
WLSLOT_K8_DOMAIN_ONLY: 82 (A-tags 0)
WLSLOT_K8_UNCHANGED: 132 (A-tags 64)
WLSLOT_K8_NO_CE_SLOT: 60 (A-tags 1)
WLSLOT_K9a_VALUE_CHANGED: 0 (A-tags 0)
WLSLOT_K9a_DOMAIN_ONLY: 0 (A-tags 0)
WLSLOT_K9a_UNCHANGED: 270 (A-tags 64)
WLSLOT_K9a_NO_CE_SLOT: 60 (A-tags 1)
WLSLOT_K9b_VALUE_CHANGED: 0 (A-tags 0)
WLSLOT_K9b_DOMAIN_ONLY: 0 (A-tags 0)
WLSLOT_K9b_UNCHANGED: 270 (A-tags 64)
WLSLOT_K9b_NO_CE_SLOT: 60 (A-tags 1)
WLSLOT_K10_VALUE_CHANGED: 0 (A-tags 0)
WLSLOT_K10_DOMAIN_ONLY: 270 (A-tags 64)
WLSLOT_K10_UNCHANGED: 0 (A-tags 0)
WLSLOT_K10_NO_CE_SLOT: 60 (A-tags 1)
WLSLOT_K11_VALUE_CHANGED: 0 (A-tags 0)
WLSLOT_K11_DOMAIN_ONLY: 0 (A-tags 0)
WLSLOT_K11_UNCHANGED: 0 (A-tags 0)
WLSLOT_K11_NO_CE_SLOT: 330 (A-tags 65)
```

## Process: the stuck review run
```
$ cat _scratch/s11c/s9b-review-r2/claude-wl-harness-split-1/outcome.json
{
  "exitCode": 1,
  "wallSeconds": 826.4182035229169,
  "unit": "s11c-guard-468f0d9e6910",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": null
}
```

```
$ grep -m2 -E 'Part::partw|OpenRead::noopen' _scratch/s11c/s9b-review-r2/claude-wl-harness-split-1/stdout
Part::partw: Part 2 of {} does not exist.
OpenRead::noopen: Cannot open {}[[2]].
```

```
$ cat _scratch/s11c/resource-pool/reservations.json
{}
```

