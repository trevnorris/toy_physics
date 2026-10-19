# Measurements: round-0 review of the S9b flow-only oracle check (generated 2026-10-19 10:44)

Generator: `_scratch/s9b_oracle/gen/check_r0_lookups.sh` (sha256sum/grep/sed/cut/wc/cat only). The legs' scratch
directories are hard-linked at `_scratch/s9b_oracle/check_review/r0_{claude,grok}_evidence/`.

## Baseline and verdicts
```
$ (cd /var/projects/toy_physics && sha256sum -c /var/projects/toy_physics/_scratch/s9b_oracle/check_review/s9b_oracle_check_review_baseline_r0.sha256)
research/pde_ledger_v3/scripts/S9b_flow_oracle_check.py: OK
research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check.wl: OK
research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/identification_map.json: OK
research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/oracle_quotes.json: OK
research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out: OK
research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_run.md: OK
research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_identifications.md: OK
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s9b_oracle/check_review/s9b_oracle_check_review_prompt_r0.md | cut -c1-64
09d28745cbe089f510fd25dc66fa1d51e93b693a77adf392bd24e4a3af42da1e
```

```
$ grep -o 'Verdict: *[A-Z][A-Z ]*' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r0_claude.md | head -1
Verdict: NEEDS REVISION
```

```
$ grep -o 'Verdict: *[A-Z][A-Z ]*' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r0_grok.txt | head -1
Verdict: NEEDS REVISION
```

## Claude 1: the signed flow residuals and the domain
```
$ sed -n '17p;111p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check.wl
residual[a_, b_, assumptions_:True] := Simplify[a-b, assumptions];
 rad=eng["radius"]; rho=eng["density"]; ass=rad>0;
```

```
$ grep -n 'residual\[' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check.wl | head -20
17:residual[a_, b_, assumptions_:True] := Simplify[a-b, assumptions];
148:   Do[emit["FLOW_S9B_SQUARED_OPERAND",branch];emit["FLOW_SQUARED_RESIDUAL",residual[branch[[2]],Power[osq,2],ass]],{branch,{{"general",general},{"retained",squared}}}];
149:   Do[emit["FLOW_S9B_SIGNED_OPERAND",signed];emit["FLOW_ORACLE_AT_COMPONENT_POINT",osq/.rad->signed["point"]];emit["FLOW_SIGNED_RESIDUAL",residual[signed["value"],osq/.rad->signed["point"],ass]],{signed,allSigned}];
158:    Do[emit["FLOW_S9B_SQUARED_OPERAND",branch];emit["FLOW_SQUARED_RESIDUAL",residual[branch[[2]],Power[calibrated,2],ass]],{branch,{{"general",general},{"retained",squared}}}];
159:    Do[emit["FLOW_S9B_SIGNED_OPERAND",signed];emit["FLOW_ORACLE_AT_COMPONENT_POINT",calibrated/.rad->signed["point"]];emit["FLOW_SIGNED_RESIDUAL",residual[signed["value"],calibrated/.rad->signed["point"],ass]],{signed,allSigned}];
186:  rr=residual[ex,oracleExchange,ass];emit["EXCHANGE_RESIDUAL",rr];
236:    cmp=residual[sourceCone,coneOp,ass];emit["OPTICAL_CHARACTERISTIC_RESIDUAL",cmp];
244:     emit["OPTICAL_CHARACTERISTIC_RESIDUAL",residual[sourceCone,coneOp/.driftRules,ass]]];
```

```
$ sed -n '114p;123p;129p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
FLOW_SQUARED_RESIDUAL: 0
FLOW_SIGNED_RESIDUAL: (Sqrt[2]*(-Sqrt[py`GM] + py`c$0*Sqrt[py`GM/py`c$0^2]))/Sqrt[py`s9b$r]
FLOW_SIGNED_RESIDUAL: Sqrt[2]*py`c$0*Sqrt[py`GM/(py`c$0^2*py`s9b$velocity$zero$radius)]
```

```
$ grep -c '^FLOW_SIGNED_RESIDUAL' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
416
```

```
$ grep -c '^FLOW_SIGNED_RESIDUAL: 0$' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
0
```

```
$ grep -c "Symbol('c_0', positive=True)" /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out
345
```

```
$ grep '^WL_LOCAL_S9B_PROFILE_ANSATZ' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out | grep -o '"Domain" -> [^}]*}' | head -2
"Domain" -> Element[c0 | ell | b | bFar | ZE | ZR | r | rr | impact | GM | n | K | m | s | rho0 | Phi | rAnchor, Reals] && c0 > 0 && ell > 0 && b > bFar > 0 && ZE > 0 && ZR > 0 && r > 0 && rr > 0 && impact > 0 && rho0 > 0 && m > 0 && rAnchor > bFar, "Counting" -> {"delta O(epsilon)", "V/c0 O(Sqrt[epsilon])", "xiPrimeSquared O(epsilon)"}
```

```
$ grep -n 'domainRules' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check.wl | head -10
109:profileComparisons[eng_,route_] := Module[{tag,obj,meta,general,squared,cases,definition,vel,ex,ass,rad,roots,allSigned,signed,dom,oraw,ometric,osq,rr,rate,nativeSpeed,conditionalMetric,mappedMass,calibrated,zeroeq,solutions,rho,rawex,caseIndex,domainAss,domainRules,casePoint,oracleAtPoint},
179:  domainRules=If[!FreeQ[rawex,Inactive[Limit]],Solve[domainAss,eng["mass"]],{{}}];
180:  emit["EXCHANGE_DOMAIN_PARAMETER_RULES",domainRules];
183:   ex=Block[{$Assumptions=ass && domainAss},Simplify[Activate[ex/.First[domainRules],Limit],ass && domainAss]]];
```

```
$ cat /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r0_claude_evidence/derive/d6_domain_census.stdout
CENSUS {engine, printed_nonzero, reduces_to_0_on_domain, branch, zero_case} -> count
{PY, True, False, general, False} -> 64
{PY, True, False, retained, False} -> 56
{PY, True, False, retained, True} -> 96
{PY, True, True, retained, False} -> 8
{PY, True, True, retained, True} -> 32
{WL, True, False, general, False} -> 64
{WL, True, False, retained, False} -> 56
{WL, True, False, retained, True} -> 24
{WL, True, True, retained, False} -> 8
{WL, True, True, retained, True} -> 8
EXAMPLES (transcript line, oracle/mode, operand name, printed residual, reduced)
{123, {"HL2", "calibrated"}, {"retained", 1}, "(Sqrt[2]*(-Sqrt[pyQGM] + pyQcD0*Sqrt[pyQGM/pyQcD0^2]))/Sqrt[pyQs9bDr]", "on domain ->", 0}
{129, {"HL2", "calibrated"}, {"retained", 3}, "Sqrt[2]*pyQcD0*Sqrt[pyQGM/(pyQcD0^2*pyQs9bDvelocityDzeroDradius)]", "on domain ->", 0}
{132, {"HL2", "calibrated"}, {"retained", 4}, "Sqrt[2]*pyQcD0*Sqrt[pyQGM/(pyQcD0^2*pyQs9bDvelocityDzeroDradius)]", "on domain ->", 0}
{135, {"HL2", "calibrated"}, {"retained", 5}, "Sqrt[2]*pyQcD0*Sqrt[pyQGM/(pyQcD0^2*pyQs9bDvelocityDzeroDradius)]", "on domain ->", 0}
{138, {"HL2", "calibrated"}, {"retained", 6}, "Sqrt[2]*pyQcD0*Sqrt[pyQGM/(pyQcD0^2*pyQs9bDvelocityDzeroDradius)]", "on domain ->", 0}
{217, {"HL4", "calibrated"}, {"retained", 1}, "(Sqrt[2]*(-Sqrt[pyQGM] + pyQcD0*Sqrt[pyQGM/pyQcD0^2]))/Sqrt[pyQs9bDr]", "on domain ->", 0}
PY line-128-type example: 0
```

## Claude 2: comparison 2's oracle operand
```
$ sed -n '84,90p' /var/projects/toy_physics/research/pde_ledger_v3/directives/S9b_flow_oracle_check_decision_list.md
2. **Continuity and exchange.**
   - *S9b operand:* each engine's printed implied exchange, with `ρ_br` live.
   - *Oracle operands:* Visser's background continuity equation (5) and the density conditions of its §7 and §8
     backgrounds, quoted as stated, including §7's statement about which profiles continuity admits. Hamilton–Lisle
     states no continuity operand, so it is `UNMATCHED`.
   - *Residual:* S9b's exchange less the oracle's, as a function of live `ρ_br`, with the set of `ρ_br` on which it
     vanishes. No oracle density law is substituted into S9b's operand.
```

```
$ sed -n '37p;165p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check.wl
 Inactive[Div][densityField[ot,ox,oy,oz] {velocityX[ot,ox,oy,oz],velocityY[ot,ox,oy,oz],velocityZ[ot,ox,oy,oz]}, {ox,oy,oz}] == 0];
 oracleExchange=Last[First[visContinuity]];
```

```
$ grep -n 'oracleExchange\|visContinuity' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check.wl | head -12
36:visContinuity = HoldComplete[Inactive[D][densityField[ot,ox,oy,oz],ot] +
49: "V4"->visMetric,"V5"->visContinuity,"V6"->visEuler,"V7"->visForce,
165: oracleExchange=Last[First[visContinuity]];
169:  emit["EXCHANGE_ORACLE_OPERAND",oracleExchange];
186:  rr=residual[ex,oracleExchange,ass];emit["EXCHANGE_RESIDUAL",rr];
```

```
$ sed -n '559p;565p;566p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
EXCHANGE_ORACLE_OPERAND: 0
EXCHANGE_IDENTIFICATION_MAP: {"background_continuity_Visser", "section7_EOS_and_density", "section8_EOS_and_density", "no_off_sheet_exchange", "S9b density remains live; only zero source is compared"}
EXCHANGE_RESIDUAL: (Sqrt[py`GM]*(3*py`rho$br$profile[py`s9b$r] + 2*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r]))/(Sqrt[2]*py`s9b$r^(3/2))
```

```
$ grep -A1 F2_oracle_continuity_form /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r0_claude_evidence/ablate/summarize_diffs_1.stdout
==================== F2_oracle_continuity_form lines base/variant: 3033 3033
  block native Wolfram computation from both transcript representations -> {'PROVENANCE': 1, 'ORACLE_TRANSCRIBED_OBJECTS': 1}
```

```
$ grep -A1 F4_oracle_hl_spatial_form /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r0_claude_evidence/ablate/summarize_diffs_3.stdout
==================== F4_oracle_hl_spatial_form lines base/variant: 3033 108
LINE COUNT DIFFERS
```

## Claude 3: the PY s9b_A binding
```
$ sed -n '68p;99p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check.wl
baseHL=localMetric[lineMetric[hlLine],1];
 "openBindings"->{py`s9b$chi->(pySpeed[[1,1]]/.pyRestrict),py`s9b$A->baseHL[[2,2]],py`s9b$U->pyProfiles[[2]]},"fixedBindings"->{}|>,
```

```
$ grep -n 'localMetric\[' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check.wl | head -6
56:localMetric[metric_,factor_] := Module[{spatial,basis},
68:baseHL=localMetric[lineMetric[hlLine],1];
69:base53=localMetric[lineMetric[visLine53],visFactor53];
70:base55=localMetric[lineMetric[visLine55],1];
```

```
$ sed -n '1238p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
OPTICAL_IDENTIFICATION_MAP: <|"conditional_bindings" -> {}, "unbound_entries" -> {py`s9b$A -> 1, py`s9b$chi -> py`c$0^2, py`s9b$U -> py`s9b$V[py`s9b$r]}, "fixed_bindings" -> {}, "oracle" -> "HL2", "clock" -> {"open_clock", clockRate}, "frame" -> "UNMATCHED: conditional common coordinates only", "mass" -> "UNMATCHED: independent oracle parameters retained"|>
```

```
$ grep -n -A4 'PY_dispersion_chi_A_U' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/identification_map.json | head -12
19: {"entry":"PY_dispersion_chi_A_U", "status":"UNMATCHED", "reason":"The printed dispersion and advection do not explicitly bind chi, A and U to the speed-identification expression, spatial metric and material profile. Every bound/unbound combination is printed as a conditional substitution, with original symbols retained in unbound entries.", "source":"PY VECTOR_DISPERSION, VECTOR_ADVECTION, SPE
20- {"entry":"WL_dispersion_aMetric_vLocal", "status":"UNMATCHED", "reason":"WL Speed binds cSquared but the supplied Metric and profile fields do not explicitly bind aMetric or vLocal to them. The flow restriction rewrites only deltaProfile and xiProfile. All bound/unbound combinations are printed conditionally.", "source":"WL SUPPLIED_INPUTS, PROFILE_ANSATZ, RESTRICTION_FLOW_ONLY"}
21-]
```

```
$ cat /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r0_claude_evidence/derive/d8_localmetric_identity.stdout
localMetric[[2,2]] for general g, factor = factor
with factor = 1 (baseHL): 1
```

## Grok 1: Hamilton-Lisle's signed operand
```
$ sed -n '78,82p' /var/projects/toy_physics/research/pde_ledger_v3/directives/S9b_flow_oracle_check_decision_list.md
1. **Flow profile.**
   - *S9b operand:* the flow-only squared velocity and each signed component. It is taken with the general solution's
     integration constant live, and again on the retained branch.
   - *Oracle operands:* Hamilton–Lisle's river velocity profile, (2), and its mass-function form, (4). Visser's §7 and
     §8 background velocities (§7's background and (54)), each with the conditions the paper attaches to it.
```

```
$ sed -n '228,236p;243,248p;662,667p' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r0_grok_evidence/papers/hl/river.tex
the Schwarzschild
%\cite{Schwarzschild}
\cite{AL04,Rothman02}
metric can be expressed in the form
\begin{equation}
\label{gullstrandmetric}
  \dd s^2
  =
  - \, \dd t_\ff^2 + (\dd r + \vff \, \dd t_\ff)^2
in units of the speed of light, at radius $r$
from a spherical object of mass $M$
\begin{equation}
  \vff = \left( {2 G M \over r} \right)^{1/2}
\end{equation}
and $t_\ff$ is
can be expressed in the Gullstrand-Painlev\'e form~(\ref{gullstrandmetric})
with infall velocity
\begin{equation}
  \vff(r) = \left[ {2 M(r) \over r} \right]^{1/2}
\end{equation}
the free-fall time $t_\ff$ being
```

```
$ grep -n 'drift' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check.wl | head -12
64:drift[metric_] := Simplify[-metric[[1,2]]/metric[[2,2]]];
139:  oraw=drift[clockPullback[ometric,clockRate]];
145:   osq=If[mode=="raw_clock",oraw,drift[conditionalMetric]];
237:    (* A second explicitly conditional map identifies the drift components
239:    driftRules=If[model["name"]=="V4",Join[First[Solve[drift[transformed]==eng["velocity"],vxO]],{vyO->0,vzO->0}],{}];
240:    If[Length[driftRules]>0,
241:     emit["CONDITIONAL_RADIAL_VELOCITY_IDENTIFICATION",driftRules];
243:     emit["OPTICAL_ORACLE_OPERAND",coneOp/.driftRules];
244:     emit["OPTICAL_CHARACTERISTIC_RESIDUAL",residual[sourceCone,coneOp/.driftRules,ass]]];
```

```
$ sed -n '48p;78p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
FLOW_ORACLE_SIGNED_OPERAND: {"HL2", "raw_clock", -(Sqrt[2]*clockRate*Sqrt[gmH/py`s9b$r])}
FLOW_ORACLE_SIGNED_OPERAND: {"HL2", "speed_matched", -(Sqrt[2]*py`c$0*Sqrt[gmH/py`s9b$r])}
```

```
$ grep -n -B1 -A4 '"velocity_frame\|velocity' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/identification_map.json | head -30
1-[
2: {"entry":"velocity_and_frame_HL", "status":"UNMATCHED", "reason":"HL I-II describes an inward river and free-fall tetrads, whereas S9b V is the material velocity in the far-field lab frame. No physical frame-identification is supplied.", "source":"HL (1)-(2), II.A; S9b record section 2; native signed-component operands"},
3: {"entry":"velocity_and_frame_Visser", "status":"UNMATCHED", "reason":"Visser v is a Newtonian background fluid velocity; S9b V is the shear carrier velocity. Their physical identity is not supplied. Both radial orientations are retained conditionally.", "source":"Visser (4)-(6), (53)-(55); S9b record section 2"},
4- {"entry":"free_fall_versus_lab_time", "status":"UNMATCHED", "reason":"HL t_ff is infall proper time; S9b uses lab time. Visser uses Newtonian lab time, but a cross-model clock identification is not supplied. A positive clock-conversion parameter is retained in the oracle Jacobian.", "source":"HL p.1 following (2); Visser pp.4-5; WL SelectedAnchoring and BackgroundAnchoring; S9b record section 2
5- {"entry":"units_c0_cgamma_oracle_speed", "status":"UNMATCHED", "reason":"S9b c0 is asymptotic light speed and cgamma is the transverse speed. HL uses light-speed units (II.A: c=G=1); Visser uses sound speed (section 7: c=1). No physical equality to the oracle speed is supplied. Raw and explicitly conditional speed-matched cones are both printed.", "source":"PY SPEED_IDENTIFICATION; WL Speed/Non
6- {"entry":"spatial_dimension", "status":"SOURCE_FIXED", "reason":"Both comparisons use three spatial directions: PY VECTOR_ADVECTION has three entries; WL Metric has three rows; HL (1) has three spatial coordinates; Visser (4) uses a 3 by 3 identity. WL dispersion is printed on a two-component angular slice; no missing wave-vector component is invented.", "source":"Named native tags; HL (1); Vis
7- {"entry":"mass_calibration_HL", "status":"UNMATCHED", "reason":"S9b orbital GM has no supplied calibration to HL G M or M(r). Oracle parameters remain independent in the general residual; a common-parameter substitution is separately labelled conditional in common unit coordinates.", "source":"S9b record section 2; HL (2),(4); O2 scoping section 5"},
--
12- {"entry":"inviscid_Visser", "status":"UNMATCHED", "reason":"The S9b optical nondissipative limit does not supply the live Newtonian inviscid material law required by Visser.", "source":"Visser theorem and section 2.1; O2 scoping section 2"},
13: {"entry":"irrotational_Visser", "status":"UNMATCHED", "reason":"A radial exterior profile permits a conditional irrotational kinematic comparison, but the scalar velocity-potential fluid model is not supplied for the shear carrier.", "source":"Visser theorem and section 2.1; native radial profiles; O2 scoping section 5"},
14- {"entry":"section7_EOS_and_density", "status":"UNMATCHED", "reason":"Visser section 7 requires constant sound speed, its stated linear EOS and radially varying density. S9b retains live rho_br and supplies no corresponding brane EOS or density law.", "source":"Visser p.15 before (53); PY and WL density operands"},
15- {"entry":"section8_EOS_and_density", "status":"UNMATCHED", "reason":"Visser section 8 assumes incompressibility with position-independent density, pressure and sound speed. S9b does not supply that brane restriction.", "source":"Visser p.16 before (54); native live density and speed operands"},
16: {"entry":"scalar_acoustics_vs_transverse_optics", "status":"UNMATCHED", "reason":"Visser propagates a scalar velocity potential; S9b supplies a transverse optical dispersion with two polarizations. A common characteristic cone is not a wave-equation identification.", "source":"Visser theorem, (86); WL Polarizations; S9b record section 2"},
17- {"entry":"no_shear_bearing_sheet", "status":"UNMATCHED", "reason":"The oracle fluid examples contain no shear-bearing sheet; S9b's optical branch belongs to such a material sheet.", "source":"O2 scoping section 5; repository survey lines 246,263-266; Visser theorem"},
18- {"entry":"no_off_sheet_exchange", "status":"UNMATCHED", "reason":"The oracle backgrounds provide no off-sheet exchange operand or projection; S9b's j_n has that role.", "source":"Visser (5), sections 7-8; O2 scoping section 5; native flow-only exchanges"},
19- {"entry":"PY_dispersion_chi_A_U", "status":"UNMATCHED", "reason":"The printed dispersion and advection do not explicitly bind chi, A and U to the speed-identification expression, spatial metric and material profile. Every bound/unbound combination is printed as a conditional substitution, with original symbols retained in unbound entries.", "source":"PY VECTOR_DISPERSION, VECTOR_ADVECTION, SPE
20- {"entry":"WL_dispersion_aMetric_vLocal", "status":"UNMATCHED", "reason":"WL Speed binds cSquared but the supplied Metric and profile fields do not explicitly bind aMetric or vLocal to them. The flow restriction rewrites only deltaProfile and xiProfile. All bound/unbound combinations are printed conditionally.", "source":"WL SUPPLIED_INPUTS, PROFILE_ANSATZ, RESTRICTION_FLOW_ONLY"}
```

```
$ grep -E 'HL_QUOTED_PROFILE|HL_COORDINATE_DRIFT|SIGNED_RESIDUAL' /var/projects/toy_physics/_scratch/s9b_oracle/check_review/r0_grok_evidence/recompute/13_recompute.kernel.stdout
HL_QUOTED_PROFILE: Sqrt[2]*Sqrt[gmH/r]
HL_COORDINATE_DRIFT: -(Sqrt[2]*Sqrt[gmH/r])
SIGNED_RESIDUAL_NEG_MINUS_DRIFT: (Sqrt[2]*clockRate*Sqrt[gmH] - Sqrt[2*GM + C*r^2])/Sqrt[r]
SIGNED_RESIDUAL_NEG_MINUS_QUOTED: (-(Sqrt[2]*clockRate*Sqrt[gmH]) - Sqrt[2*GM + C*r^2])/Sqrt[r]
```

## The stalled F4 ablation (Claude leg)
```
$ sed -n '140,175p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check.py
  artifact=SCRIPT if f==SCRIPT.name else ASSETS/f
  data['provenance'][f]={'path':str(artifact),'sha256':sha(artifact)}
 (HERE/'inputs.json').write_text(json.dumps(data,indent=2)+'\n')
 print('WORKING_ON: native Wolfram computation from both transcript representations',flush=True)
 messages=[]
 completed=False
 with (HERE/'latest_kernel.out').open('w') as output:
  proc=subprocess.Popen(['math','-script',str(ASSETS/'S9b_flow_oracle_check.wl')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
  for line in proc.stdout:
   print(line,end='',flush=True)
   output.write(line)
   if re.match(r'^[A-Za-z$][A-Za-z0-9$]*::',line): messages.append(line.rstrip())
   if line.startswith('INSTRUMENT_COMPLETION: '): completed=True
  code=proc.wait()
 print('DRIVER_COMPLETION_GUARD:',json.dumps({'kernel_exit':code,'completion_record':completed,'messages':messages}),flush=True)
 return code or (0 if completed and not messages else 2)

if __name__=='__main__':
 try: sys.exit(main())
 except Exception as e:
  print('SOURCE_OR_DRIVER_ERROR:',repr(e),flush=True)
  print('OPERAND: UNMATCHED; RESIDUAL: UNMATCHED; GUARD: source/driver failure',flush=True)
  raise
```

