# Measurements — S9b repair build review, round 8 (generated 2026-10-16 15:56)

Generator: `_scratch/s9b_build/gen/s9b_build_r8_lookups.sh` (sha256sum/grep/sed/cat/wc/git show only). Line
numbers are the reviewed round-8 files' (the baseline below).

## Baseline, prompt and runs
```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_review_baseline_r8.sha256 | grep -c ': OK$'
28
```

```
$ wc -l < _scratch/s9b_build/s9b_repair_build_review_baseline_r8.sha256
28
```

```
$ sha256sum _scratch/s9b_build/s9b_repair_build_review_prompt_r8.md
c8b72a72668b09073c9418c09b6f78fb1a89d18b5277d5d5544be60b116bf625  _scratch/s9b_build/s9b_repair_build_review_prompt_r8.md
```

```
$ for f in research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py research/pde_ledger_v3/scripts/S9b_exports.py research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl; do echo "$(git show 9af9f705:$f | sha256sum | cut -c1-16) $(sha256sum $f | cut -c1-16) $f"; done
a230d052203aee81 a230d052203aee81 research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
0900676cd86f24ed 0900676cd86f24ed research/pde_ledger_v3/scripts/S9b_exports.py
2028badd3772d68d 2028badd3772d68d research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
```

```
$ grep -h '"exitCode"\|"wallSeconds"' _scratch/s9b_build/r8_harness_output/py/full1/outcome.json _scratch/s9b_build/r8_harness_output/wl/full-final/outcome.json | head -4
  "exitCode": 0,
  "wallSeconds": 18176.372185758024,
    "exitCode": 0,
    "wallSeconds": 18176.297543363005,
```

```
$ grep -o '"memory.peak": "[0-9]*"' _scratch/s9b_build/r8_harness_output/py/full1/resource-samples.jsonl | tr -dc '0-9\n' | sort -n | tail -1
20921196544
```

```
$ cat research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round8_source_hashes.json
{
  "research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl": "2028badd3772d68d43659e06d6ebc6644843dd0fb77563fc6fff37255f241561",
  "research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl": "740c3559a1c9ca12225cf85e64bff3dc5a2d54fa771d527422af03dd623f5a2f",
  "research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation_checks.wl": "5a84c10e5e016275aa7d98465dad1ae6124c12687e52219a8175c2222e9f3db8",
  "research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round8_validate.py": "a766848c6c07c31983802e9c69a2c3c5145324269a5379b519e3d7f2d6a8696a",
  "scripts/s11c_guarded_run.py": "c747c266e78ed519d3e31add8fffed85a0aafabc9e01528ffd3eb157f8f1dcdd"
}
```

```
$ grep 'ablation.py' _scratch/s9b_build/r8_harness_output/py/builder_scratch/full1/source-pins.json
  "research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py": "5209452b981deeb0011ba046cd64f585f06f2fe6a000711a71ccd33fe3cbab91",
```

## Verdicts
```
$ head -1 _scratch/s9b_build/s9b_repair_build_review_r8_claude.md
**Verdict: NEEDS REVISION**
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_repair_build_review_r8_grok.txt
Verdict: CLEAR
```

## Finding 1: the Wolfram K6, K7 and K10 searches
```
$ sed -n '201,204p;208p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
 "K4c"->{"slope",{"WL_LOCAL_S9B_PROFILE_ANSATZ","XiDerivative"},"\"XiDerivative\" -> xiSlope"},
 "K5"->{"odd",{"WL_LOCAL_S9B_NONRECIPROCAL_CONSTRUCTION","OneForm",1},"oddOneFormLocal = D["},
 "K6"->{"ray",{"WL_S9B_A_RADAR_LOG_G100","OnDomain"},"finiteSlope = Expand["},
 "K7"->{"ray",{"WL_S9B_B_CONDITION_DEFLECTION"},"If[quantifierMode === \"FIXED_B\""},
 "K10"->{"ray",{"WL_S9B_A_NONRECIPROCAL_PATH_DEPENDENCE","OnDomain"},"physicalPropagation = branchGateLocal"},
```

```
$ sed -n '213,225p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
  {kind,path}=Take[spec,2];note["witness begin",name];
  aa=part[a,path];bb=part[b,path];
  location=<|"Tag"->First[path],"PayloadLocation"->Rest[path],"SitePart"->"value",
    "Source"->engine,"SourceLine"->sourceLine[Last[spec]],"EntryConstruction"->Last[spec]|>;
  If[kind==="ray",
    (* Read the full candidate and fixing equations. The finite search stops
       at the unresolved full ray boundary problem; no guessed turning radius,
       impact, endpoint cutoff or inactive-domain decision is used. *)
    entry=<|"Result"->"none found","Site"->location,"Attempt"-><|
      "CandidateParts"->{aa,bb},"ProfilesTried"->Normal[profiles],
      "ConstructionEquationsAfterProfiles"->{a["State"]["RayConstruction"],b["State"]["RayConstruction"]}/.Normal[profiles],
      "UnassignedDependentQuantities"->{Global`rTurnFlyby,Global`rTurnRadar,Global`impactRadar},
      "Reason"->"no exact solution of the full endpoint/turning construction and all ranged conditions was obtained; no scalar point or bound-variable substitution was made"|>|>;
```

```
$ sed -n '112,113p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
finiteExactQ[z_] := (NumericQ[z] && Precision[z]===Infinity &&
  FreeQ[z,Indeterminate|Undefined|DirectedInfinity] && TrueQ[Im[z]==0]);
```

```
$ grep -c 'Solve\[\|FindInstance\[' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
0
```

```
$ sed -n '18p' research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round8_measurements.md
K6, K7 and K10 print **none found**. Their candidate parts and the endpoint/turning construction after the test-profile substitution are printed, together with the unresolved dependent radii and radar impact. The search obtained no exact solution of that whole ray construction and its ranged conditions. It assigns no guessed turning radius, substitutes no bound variable, treats no inactive conditi
```

```
$ sed -n '356p;374,375p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
radarSource = roundTrip; (* K6 *)
  Keys[radarKernels]];
emit["LOCAL_RADAR_SLOPE_DERIVATION", <|"Source" -> radarSource,
```

```
$ sed -n '174p;203p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
physicalPropagation = branchGateLocal /. branchRules;
  "Domain" -> (r > 0 && physicalTraversal && physicalPropagation && speedChart),
```

```
$ sed -n '154,166p' _scratch/s9b_build/s9b_repair_build_review_r8_claude_evidence/wl_out/wl_witness_recheck.stdout
K6cand emitted parts {corrupted, baseline}: {Inactive[Integrate][-(velocityProfile[integrationRadius]/c0^2), {integrationRadius, Sqrt[b^2 + ZE^2], Sqrt[b^2 + ZR^2]}], 0}
K6cand free symbols: {{b, c0, integrationRadius, ZE, ZR}, {}}
K6cand values {corrupted, baseline}: {-1/80*Log[29/20], 0}
K6cand finite exact: {True, True}
K6cand difference N[,30]: -0.0046445444554060379218506057027421353727238911366714481096`30.
K6cand attached PROFILE_ANSATZ domain {corrupted, baseline}: {True, True}
K6cand integrands after profiles and point (Inactive forms): {Inactive[Integrate][-1/40*1/integrationRadius, {integrationRadius, 2*Sqrt[5], Sqrt[29]}], 0}
K10cand emitted domains {corrupted, baseline}: {r > 0 && muPerp[r]/rhoBr[r] > velocityProfile[r]^2*(1 + Derivative[1][xiProfile][r]^2) &&  !(Element[muPerp[r]/rhoBr[r], Reals] && muPerp[r]/rhoBr[r] > 0 && r^2*(1 + Derivative[1][xiProfile][r]^2) > 0) && c0*(1 + deltaProfile[r]) == Sqrt[muPerp[r]/rhoBr[r]], r > 0 && muPerp[r]/rhoBr[r] > velocityProfile[r]^2*(1 + Derivative[1][xiProfile][r]^2) && Ele
K10cand free symbols: {{c0, r}, {c0, r}}
K10cand domain at inputs {corrupted, baseline}: {False, True}
K10cand conjuncts at inputs {corrupted, baseline}: {{True, True, False, True}, {True, True, True, True, True, True}}
K10cand PROFILE_ANSATZ domain {corrupted, baseline}: {True, True}
K10cand OneForm[[1]] at inputs {corrupted, baseline}: {-1207049/145800000, -1207049/145800000}
```

```
$ for d in _scratch/s11c/s9b-review-r8/claude-wl-witness-recheck-*; do echo "$d $(grep -o '"exitCode": [0-9-]*' $d/outcome.json | head -1)"; done
_scratch/s11c/s9b-review-r8/claude-wl-witness-recheck-1 "exitCode": 0
_scratch/s11c/s9b-review-r8/claude-wl-witness-recheck-2 "exitCode": 0
_scratch/s11c/s9b-review-r8/claude-wl-witness-recheck-3 "exitCode": 0
```

```
$ grep -n '| Wolfram | K6 \|| Wolfram | K7 \|| Wolfram | K10 ' _scratch/s9b_build/s9b_repair_build_review_r8_grok.txt
32:| Wolfram | K6 | none found; `WL_S9B_A_RADAR_LOG_G100` OnDomain; tried printed | confirms |
33:| Wolfram | K7 | none found; `WL_S9B_B_CONDITION_DEFLECTION`; tried printed | confirms |
37:| Wolfram | K10 | none found; path-dependence OnDomain, gate after negation inside the inactive quantifier; tried printed | confirms |
```

## Finding 2: the SymPy software-checks line
```
$ sed -n '606,611p;663,664p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
def self_checks():
    # Software checks, with no physics expected-value acceptance test.
    checks={
        'syntax_whitespace':structural_mark('Tuple(Integer(1), Integer(2))','Tuple( Integer(1),Integer(2) )')=='identical',
        'argument_order':structural_mark('Tuple(Integer(1), Integer(2))','Tuple(Integer(2), Integer(1))')=='not compared',
        'no_algebra':structural_mark('Add(Integer(1),Integer(1))','Integer(2)')=='not compared',
    emit('SOFTWARE_CHECKS',checks)
    if not all(checks.values()):raise ValueError(('software checks',checks))
```

```
$ git show 9af9f705:research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py | grep -c 'SOFTWARE_CHECKS'
0
```

```
$ grep 'SOFTWARE_CHECKS' _scratch/s9b_build/r8_harness_output/py/checks3/stdout | cut -c1-300
PY_LOCAL_S9B_ABLATION_SOFTWARE_CHECKS: Tuple(Tuple(Str('syntax_whitespace'), true), Tuple(Str('argument_order'), true), Tuple(Str('no_algebra'), true), Tuple(Str('float_decimal_tail'), true), Tuple(Str('float_precision'), true), Tuple(Str('python_float_token'), true), Tuple(Str('binder_names_literal
```

## SymPy K7 and K10: the printed none-found searches
```
$ sed -n '177p;181p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
    'K7':('PY_S9B_B_CONDITION_DEFLECTION','reduced_condition','local_mass=substitute(mass,replacement)'),
    'K10':('PY_S9B_A_DEFLECTION_G000','gate','branch_domain=sp.And(observable_gate,travel_domain)'),
```

```
$ grep -n 'K7\|K10' _scratch/s9b_build/s9b_repair_build_review_r8_claude.md | grep -i 'none found'
7:**1. Wolfram K6 and K10 print "none found" without any search, and a witness exists for each on a site part (K7 has the same unsearched branch).**
46:| K7 | none found; condition search of a site part, printed | record sound | none found, nothing tried | no witness shown either way |
50:| K10 | none found; domain search; the tried part holds the negated gate (`Eq`→`Ne`) | record sound | none found, nothing tried | domain witness exists (Finding 1) |
```

```
$ grep -n '| SymPy | K7 \|| SymPy | K10 ' _scratch/s9b_build/s9b_repair_build_review_r8_grok.txt
19:| SymPy | K7 | none found; value site `PY_S9B_B_CONDITION_DEFLECTION` payload `(4,)`; tried printed | confirms |
23:| SymPy | K10 | none found; domain site `PY_S9B_A_DEFLECTION_G000` payload `(0,0,1)`; tried printed | confirms |
```

