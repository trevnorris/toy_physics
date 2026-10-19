# Measurements: S9b record review, round 0 (generated 2026-10-19 02:30)

Generator: `_scratch/s9b_build/gen/s9b_record_r0_lookups.sh` (sha256sum/grep/sed/head/wc only). The legs'
evidence is hard-linked at `_scratch/s9b_build/record/s9b_record_review_r0_{claude,grok}_evidence/`.

## The reviewed files and the prompt
```
$ (cd /var/projects/toy_physics && sha256sum -c /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_baseline_r0.sha256 | grep -c ': OK$')
4
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_prompt_r0.md | cut -c1-64
3e851de02a342f3343e44f5615ffced8a171ab9eb4be7da8e3da80c762990a71
```

## Verdicts
```
$ head -1 /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r0_claude.md
**Verdict: NEEDS REVISION.** The record needs four repairs. The physics transcriptions and the comparison accounting are faithful, and the generator reproduces its filed file exactly.
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r0_grok.txt
Verdict: NEEDS REVISION
```

## Both legs: the generator reproduces the filed measurements
```
$ sha256sum /var/projects/toy_physics/research/pde_ledger_v3/steps/_measurements/S9b_record_measurements.md | cut -c1-64
cb320c4b23cc5c62f770e6c748385345e08454c5845f06efa55c9c269055eea9
```

```
$ grep -h 'REGENERATED_IDENTICAL_TO_FILED' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r0_claude_evidence/rerun_generator.stdout | head -2
REGENERATED_IDENTICAL_TO_FILED
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r0_grok_evidence/S9b_record_measurements.regen.md | cut -c1-64
cb320c4b23cc5c62f770e6c748385345e08454c5845f06efa55c9c269055eea9
```

## Claude 1: the within-engine deflection/radar statements
```
$ sed -n '7,8p;192,193p;267p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
substrate-derived light equation or a demonstrated match to GR. At the selected comparison grades,
each engine prints the same profile condition for deflection and radar slope within that engine.
M6 counts marker-only comparisons for all 18 γ-family rows. Since the within-engine conditions do
not differ, the record counts no independent profile test from the pair of observables at the
conditional deflection/radar copies have the same native condition fields; their cross-engine
```

```
$ sed -n '824p;849p' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
labels are retained without turning equal within-engine conditions into two independent tests.
remain OPEN; neither radial flow nor the printed within-engine equality of profile conditions
```

```
$ grep -c '^PAIR' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r0_claude_evidence/e7_wl_condition_fields.stdout
23
```

```
$ head -5 /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r0_claude_evidence/e7_wl_condition_fields.stdout
PAIR WL_S9B_B_CONDITION_DEFLECTION leaves 57 57 differing_leaf_paths 3
   DIFF Cases/0/InverseWeight | DEFL: (velocityProfile[r]^2 - 2*r*velocityProfile[r]*Derivative[1][velocityProfile][r] + c0^2*(2*r*Derivative[1][deltaProfile][r] + Derivative[1][xiProfile][r]^2))/(2*c0^2) | RADR: (velocityProfile[r]^2 - 2*r*velocityProfile[r]*Derivative[1][velocityProfile][r] + c0^2*(2*r*Derivative[1][deltaProfile][r] + Derivative[1][xiProfile][r]^2))/(2*c0^3)
   DIFF Cases/0/ReferenceInverse | DEFL: (2*GM)/(c0^2*r) | RADR: (2*GM)/(c0^3*r)
   DIFF Cases/0/VelocityReduction/SquaredVelocityEquation | DEFL: (-2*GM)/(c0^2*r) + yProfile[r]/(2*c0^2) + r*Derivative[1][deltaProfile][r] + Derivative[1][xiProfile][r]^2/2 - (r*Derivative[1][yProfile][r])/(2*c0^2) == 0 | RADR: (-2*GM)/(c0^3*r) + yProfile[r]/(2*c0^3) + (r*Derivative[1][deltaProfile][r])/c0 + Derivative[1][xiProfile][r]^2/(2*c0) - (r*Derivative[1][yProfile][r])/(2*c0^3) == 0
   FIELD Cases/0/ProfileEquation = velocityProfile[r]^2 + c0^2*(4*r*Derivative[1][deltaProfile][r] + Derivative[1][xiProfile][r]^2 + 2*r^2*Derivative[2][deltaProfile][r] + 2*r*Derivative[1][xiProfile][r]*Derivative[2][xiProfile][r]) == 2*r*(r*Derivative[1][velocityProfile][r]^2 + velocityProfile[r]*(Derivative[1][velocityProfile][r] + r*Derivative[2][velocityProfile][r]))
```

```
$ grep -h 'IDENTICAL_PAYLOAD' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r0_claude_evidence/e6_defl_vs_radar.stdout | head -3
py PY_S9B_B_CONDITION_DEFLECTION IDENTICAL_PAYLOAD len 231312
py PY_S9B_B_IMPLIED_JN_DEFLECTION IDENTICAL_PAYLOAD len 2500903
py PY_S9B_C_CONSTANT_CONDITION_DEFLECTION IDENTICAL_PAYLOAD len 96345
```

```
$ grep -h "count of 'RADAR'" -A1 /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r0_claude_evidence/e12_source_lines.stdout | head -3
== F1 generator: any deflection/radar within-engine pairing? count of 'RADAR' (case-insens) in generator:
1
```

```
$ grep -ic 'radar' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_record_measurements.py
1
```

## Claude 2: the K11 claim limit
```
$ sed -n '72,74p' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_build_r7_review_disposition.md
- **K11.** The SymPy harness evaluates 3 tags and prints 376 `NOT_EVALUATED`. The Wolfram harness reaches 10 tags and
  prints 364 `NOT_EVALUATED`. Neither reaches the shared `A_NONRECIPROCAL_*` tags. Claim limit: K11 supports claims
  about the local one-form and closedness construction only, not item 8's shared observables.
```

```
$ grep -n 'K1–K9b, K11' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_build_r11_review_disposition.md
73:- **Nothing else moved.** Every other knife (K1–K9b, K11), the identity marks and the baseline match round 10's
```

```
$ sed -n '96,98p;293p;370,371p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
They are comparison inputs, not a demanded construction outcome. The nonreciprocal path-dependence
predicate prints `False` on its domain in each engine: endpoint-only in this radial setting, not
zero one-way nonreciprocity in general; outside that domain it is `NOT_ESTABLISHED`.
| A_NONRECIPROCAL_PATH_DEPENDENCE | 33 value Boolean samples print `boolean_difference:false`; native on-domain predicate False | No along-ray domain equality or nonradial path statement |
  runs mentioned by the disposition source no result here. K11 unreached parts remain “not evaluated”
  and “not compared”; identity marks or witness absence cannot supply cross-engine equality.
```

```
$ grep -c 'local one-form and closedness construction only' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
0
```

```
$ head -6 /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r0_claude_evidence/e11_nonrecip_coverage.stdout
PY knives type list count 14
PY knife K1 mentions A_NONRECIPROCAL_PATH_DEPENDENCE 0
PY knife K2 mentions A_NONRECIPROCAL_PATH_DEPENDENCE 0
PY knife K3 mentions A_NONRECIPROCAL_PATH_DEPENDENCE 0
PY knife K4a mentions A_NONRECIPROCAL_PATH_DEPENDENCE 0
PY knife K4b mentions A_NONRECIPROCAL_PATH_DEPENDENCE 0
```

```
$ grep -n 'A_NONRECIPROCAL_PATH_DEPENDENCE' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_build_r11_review_disposition.md | head -3
65:- **K10's witness.** It is at `WL_S9B_A_NONRECIPROCAL_PATH_DEPENDENCE`, domain part, with evaluations
```

## Claude 3: the register's second route
```
$ sed -n '672p;879p' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
doubt"* — and says *"⇒ ⭐ **Run this route for every sector with identified prior art**, ⛔ not only the
⇒ ⭐ **Run this route for every sector with identified prior art**, ⛔ not only the light sector.
```

```
$ sed -n '788,790p' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
**Second route.** Read the O2 record against scoping §5 and its source qualifications. That section
identifies surface-mechanics, diffuse-interface/evolving-reference and river/acoustic-flow **candidate
oracles**, with explicit dimension, model and access limits. O2 reproduces no identified prior-art
```

```
$ sed -n '845,846p' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
**Second route.** The identified prior-art objects in this pass are the supplied GR/PPN
comparison references. S9b establishes conditional profile requirements, not a reproduced GR
```

```
$ sed -n '452,453p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
For the second route the identified prior-art objects are the supplied GR/PPN references. S9b
prints conditional matching requirements, not a reproduced GR solution or a reproduced mechanical
```

```
$ sed -n '246p' /var/projects/toy_physics/docs/medium_requirements_and_prior_art.md
| gravity = **medium flowing into** the defect | **KNOWN** separately | Cahill (2003); ⭐ Hamilton–Lisle *river model*, *Am. J. Phys.* **76** 519 (2008); Sbitnev; Unruh acoustic sinks — ⛔ all **3D**, no off-sheet bulk |
```

```
$ grep -ic 'river' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
0
```

## Claude 4: the census command's revision
```
$ sed -n '171p;179p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_record_measurements.py
            old=subprocess.run(['git','show',rev+':research/pde_ledger_v3/'+path],cwd=ROOT,check=True,capture_output=True).stdout
        old=subprocess.run(['git','show','HEAD:'+p],cwd=ROOT,check=True,capture_output=True,text=True).stdout
```

```
$ grep -n 'ENTRY_COUNT' /var/projects/toy_physics/research/pde_ledger_v3/steps/_measurements/S9b_record_measurements.md | head -4
21288:HEAD ENTRY_COUNT 16
21305:WORKTREE ENTRY_COUNT 17
```

```
$ sed -n '811p' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
`scripts/S9b_record_measurements.py`. **16 → 17 entries**, all OPEN. **New entry by target:
```

```
$ git -C /var/projects/toy_physics rev-parse --short HEAD
caa0710d
```

## Grok 1: the R-S21-01 requirement's implied exchange
```
$ sed -n '552,566p' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
  Part C carries the separately printed constant-speed condition
  `GM=r[c₀²ξw′²−2rVV′+V²]/4`, fixed-ratio condition
  `GM=r[c₀²(n−1)rf′+c₀²ξw′²−2rVV′+V²]/4`, and power condition
  `GM=r[2c₀²srf′+c₀²ξw′²−2rVV′+V²]/4`, with WL's corresponding anchor/profile equations.
  Fixed-ratio and power objects retain both `ρ>0` and `ρ≤0` labels and exact response domains;
  `ρ₀>0` remains supplied. The first-order-f expressions supply no finite-density or f–ε extension.
  Each condition includes its printed implied `jn`, with `ρbr` symbolic and differentiated:
  the coordinate mass expression is `jn=−ρbr V′−Vρbr′−2ρbr V/r`, with both signed-velocity
  components, zero-locus differentiability/flux joins, excluded nonreal branches and printed
  integration-constant domains retained. The complete implied-exchange values remain uncompared.
  The requirements keep δ, ξw, ρbr, f, n, s, GM and all remaining printed inputs free;
  they select no profile, exponent, flux value/sign or GM normalization. The deflection and
  radar conditions are not counted as independent constraints where each engine prints the
  same condition.
- **on failure** — The conditional reference-matching statement for the affected observable
```

```
$ sed -n '574,575p' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
  Speed-, tilt- and bulk-only exchange prints zero; flow-only retains the signed live-density
  exchange. Conditions/exchange/γ values have no formed cross-engine comparison; the domains
```

```
$ grep -A1 'R_SPEED_ONLY_B_IMPLIED_JN_DEFLECTION/Cases/0/ImpliedExchange' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r0_grok_evidence/extract_fields_stdout.txt | head -2
  R_SPEED_ONLY_B_IMPLIED_JN_DEFLECTION/Cases/0/ImpliedExchange
    jn[r] == 0
```

## Grok 2: F_MASS_SOLUTION's printed residual
```
$ sed -n '248,249p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
Under that premise, `F_MASS_SOLUTION` prints `V=Φ/(4πr²ρbr)` and a mass residual 0.
WL also prints the general `C[1]/(r²rhoBr[r])` solution; the flux fixes the displayed integration
```

```
$ grep 'PY_S9B_F_MASS_SOLUTION' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out | head -1 | grep -c 'Residual'
0
```

```
$ grep 'WL_S9B_F_MASS_SOLUTION' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out | head -1 | grep -o '"Residual" -> [^,|]*'
"Residual" -> 0
```

