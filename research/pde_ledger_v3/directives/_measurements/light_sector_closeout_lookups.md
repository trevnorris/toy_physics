# Measurements: light-sector closeout (generated 2026-10-19 18:49)

Generator: `_scratch/s9b_oracle/gen/light_sector_closeout_lookups.sh`. Literal retrieval and counting only.
The oracle-check transcript is `scripts/out/S9b_flow_oracle_check.out` at `8d91f448`; its other outcome
lookups are in `S9b_flow_oracle_check_outcome_lookups.md` beside this file.

## C1. Transcript identity and the check's acceptance
```
$ sha256sum /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | cut -c1-64
e66ff63512e4f044bab2a60eb327665db865179531fe6077622de8fab553a25e
```

```
$ git -C /var/projects/toy_physics log --format='%h %ad %s' --date=short -1 8d91f448
8d91f448 2026-10-19 Accept the S9b flow-only oracle check (review clear)
```

```
$ sed -n '19,21p' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_record_r2_review_disposition.md
- **A user decision this round.** Asked whether a check of the flow-only condition against the river/acoustic-flow
  prior art should be queued, the user chose "Queue it after record". The check is new computation, done after the
  record is accepted and outside it: Codex writes the check and it is reviewed under G1.
```

## C2. Which oracle each exact squared-profile match (zero test True) comes from: nearest preceding oracle operand
```
$ grep -n -F '"residual_tag" -> "FLOW_SQUARED_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "value" -> True' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | cut -d: -f1 | while read n; do awk -v n=$n 'NR<n && /^FLOW_ORACLE_OPERAND/ {fo=$0} NR<n && /^RESIDUAL_OPERANDS/ {op=$0} NR==n {match(fo, /"equations" -> \{[^}]*\}/); print substr(fo, RSTART, RLENGTH) " | " op; exit}' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out; done | sort | uniq -c
      2 "equations" -> {"HL1", "HL2"} | RESIDUAL_OPERANDS: {(2*py`GM)/py`s9b$r, (2*py`GM)/py`s9b$r}
      2 "equations" -> {"HL1", "HL4"} | RESIDUAL_OPERANDS: {(2*py`GM)/py`s9b$r, (2*py`GM)/py`s9b$r}
      1 "equations" -> {"HL2"} | RESIDUAL_OPERANDS: {(2*py`GM)/py`s9b$r, (2*py`GM)/py`s9b$r}
      1 "equations" -> {"HL4"} | RESIDUAL_OPERANDS: {(2*py`GM)/py`s9b$r, (2*py`GM)/py`s9b$r}
      2 "equations" -> {"V7background"} | RESIDUAL_OPERANDS: {(2*py`GM)/py`s9b$r, (2*py`GM)/py`s9b$r}
      2 "equations" -> {"V7background", "V53", "V86", "V87"} | RESIDUAL_OPERANDS: {(2*py`GM)/py`s9b$r, (2*py`GM)/py`s9b$r}
```

```
$ grep -n -F '"residual_tag" -> "FLOW_SQUARED_RESIDUAL", "domain" -> "WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "value" -> True' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | cut -d: -f1 | while read n; do awk -v n=$n 'NR<n && /^FLOW_ORACLE_OPERAND/ {fo=$0} NR<n && /^RESIDUAL_OPERANDS/ {op=$0} NR==n {match(fo, /"equations" -> \{[^}]*\}/); print substr(fo, RSTART, RLENGTH) " | " op; exit}' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out; done | sort | uniq -c
      2 "equations" -> {"HL1", "HL2"} | RESIDUAL_OPERANDS: {(2*GM)/r, (2*GM)/r}
      2 "equations" -> {"HL1", "HL4"} | RESIDUAL_OPERANDS: {(2*GM)/r, (2*GM)/r}
      1 "equations" -> {"HL2"} | RESIDUAL_OPERANDS: {(2*GM)/r, (2*GM)/r}
      1 "equations" -> {"HL4"} | RESIDUAL_OPERANDS: {(2*GM)/r, (2*GM)/r}
      2 "equations" -> {"V7background"} | RESIDUAL_OPERANDS: {(2*GM)/r, (2*GM)/r}
      2 "equations" -> {"V7background", "V53", "V86", "V87"} | RESIDUAL_OPERANDS: {(2*GM)/r, (2*GM)/r}
```

### C2b. The signed comparison under the same substitution, and a clock-factor variant
```
$ sed -n '724,729p;741p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | grep -E '^(RESIDUAL_OPERANDS|FLOW_SIGNED_RESIDUAL|RESIDUAL_ZERO_TEST)'
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FLOW_SQUARED_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "value" -> True|>
RESIDUAL_OPERANDS: {-(Sqrt[2]*Sqrt[py`GM/py`s9b$r]), Sqrt[2]*py`c$0*Sqrt[py`GM/(py`c$0^2*py`s9b$r)]}
FLOW_SIGNED_RESIDUAL: -2*Sqrt[2]*Sqrt[py`GM/py`s9b$r]
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FLOW_SIGNED_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "value" -> False|>
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FLOW_SIGNED_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/2", "value" -> True|>
```

```
$ grep '^RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FLOW_SIGNED_RESIDUAL"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | grep -o '"value" -> [A-Za-z]*|>' | sort | uniq -c
     40 "value" -> False|>
    140 "value" -> True|>
```

```
$ grep -m1 -F 'py`GM == clockRate^2*gmH|>' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FLOW_SQUARED_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "value" -> py`GM == clockRate^2*gmH|>
```

```
$ grep '^RESIDUAL_ZERO_TEST: <|"residual_tag" -> "EXCHANGE_RESIDUAL"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | grep -o '"value" -> .*' | sort | uniq -c
      4 "value" -> 3*py`rho$br$profile[py`s9b$r] + 2*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r] == 0|>
      4 "value" -> 3*rhoBr[r] + 2*r*Derivative[1][rhoBr][r] == 0|>
     10 "value" -> True|>
```

## C3. The optical comparison: oracles in the identification map, and the oracle beside each printed zero
```
$ grep '^OPTICAL_IDENTIFICATION_MAP' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | grep -o '"oracle" -> "[^"]*"' | sort | uniq -c
     32 "oracle" -> "HL2"
     32 "oracle" -> "HL4"
     32 "oracle" -> "V4"
     32 "oracle" -> "V7_minus"
     32 "oracle" -> "V7_plus"
     32 "oracle" -> "V8_minus"
     32 "oracle" -> "V8_plus"
```

```
$ grep -n '^OPTICAL_CHARACTERISTIC_RESIDUAL: 0$' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | cut -d: -f1 | while read n; do awk -v n=$n 'NR>n-12 && NR<=n+3 && /^OPTICAL_IDENTIFICATION_MAP/' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | grep -o '"oracle" -> "[^"]*"' | sort -u | sed "s/^/line $n: /"; done
line 10866: "oracle" -> "V4"
line 10890: "oracle" -> "V4"
line 17075: "oracle" -> "V4"
```

```
$ grep -c '^OPTICAL_CHARACTERISTIC_RESIDUAL:' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
256
```

```
$ grep '^RESIDUAL_ZERO_TEST: <|"residual_tag" -> "OPTICAL_CHARACTERISTIC_RESIDUAL"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | sed -E 's/.*"value" -> (.*)\|>$/\1/' | awk '{ if ($0=="True"||$0=="False") print $0; else if ($0 ~ /s9b[$]V[[]|s9b[$]U|velocityProfile[[]|vLocal/) print "equation containing a live velocity symbol"; else print "equation with no velocity symbol" }' | sort | uniq -c
    240 equation containing a live velocity symbol
     13 equation with no velocity symbol
      3 True
```

```
$ grep '^RESIDUAL_ZERO_TEST: <|"residual_tag" -> "OPTICAL_CHARACTERISTIC_RESIDUAL"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | sed -E 's/.*"value" -> (.*)\|>$/\1/' | grep -v 's9b[$]V[[]\|s9b[$]U\|velocityProfile[[]\|vLocal\|^True$' | sort | uniq -c | sort -rn | cut -c1-200
      2 (py`s9b$chi - clockRate^2*soundO^2)*(otherQ^2 + radialQ^2 + tangentQ^2) == 0
      2 (py`c$0 - clockRate*soundO)*(py`c$0 + clockRate*soundO)*(otherQ^2 + radialQ^2 + tangentQ^2) == 0
      2 (py`c$0^2 - py`s9b$chi)*(otherQ^2 + radialQ^2 + tangentQ^2) == 0
      1 radialQ - radialQ/py`s9b$A == 0
      1 radialQ - radialQ/aMetric == 0
      1 (radialQ^2*(-py`s9b$chi + clockRate^2*py`s9b$A*soundO^2))/py`s9b$A == (py`s9b$chi - clockRate^2*soundO^2)*(otherQ^2 + tangentQ^2)
      1 py`c$0^2*(otherQ^2 + radialQ^2 + tangentQ^2) == py`s9b$chi*(otherQ^2 + radialQ^2/py`s9b$A + tangentQ^2)
      1 clockRate^2*soundO^2*(radialQ^2 + tangentQ^2) == (c0^2*(radialQ^2 + aMetric*tangentQ^2))/aMetric
      1 clockRate^2*soundO^2*(otherQ^2 + radialQ^2 + tangentQ^2) == (py`c$0^2*(radialQ^2 + py`s9b$A*(otherQ^2 + tangentQ^2)))/py`s9b$A
      1 (c0 - clockRate*soundO)*(c0 + clockRate*soundO)*(radialQ^2 + tangentQ^2) == 0
```

## C4. Quoted oracle equations and their stored conditions (the check's quote file)
```
$ python3 -I -c 'import json,sys; print(json.dumps({sys.argv[2]: json.load(open(sys.argv[1]))[sys.argv[2]]}))' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/oracle_quotes.json V4
{"V4": {"citation": "Visser gr-qc/9712010v2 (15 Dec 1997), local sources/gr-qc_9712010.pdf, pp.3-4, equation (4)", "quoted": "g_{mu nu}(t,x) = (rho/c) [[-(c^2-v^2), -v^T], [-v, I]]", "conditions": "I 
is the 3 by 3 identity; barotropic, inviscid fluid and irrotational flow; scalar acoustic velocity potential; Newtonian flat physical space and lab time. Overall factor rho/c is live and nonzero."}}
```

```
$ python3 -I -c 'import json,sys; print(json.dumps({sys.argv[2]: json.load(open(sys.argv[1]))[sys.argv[2]]}))' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/oracle_quotes.json V5
{"V5": {"citation": "Visser, same local copy, p.4, section 2.1, equation (5)", "quoted": "partial_t rho + nabla dot (rho v) = 0", "conditions": "Background continuity; steady exterior specialization o
f the stated equation. No off-sheet source appears."}}
```

```
$ python3 -I -c 'import json,sys; print(json.dumps({sys.argv[2]: json.load(open(sys.argv[1]))[sys.argv[2]]}))' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/oracle_quotes.json V6
{"V6": {"citation": "Visser, same local copy, p.4, section 2.1, equation (6)", "quoted": "rho d v/dt = rho [partial_t v + (v dot nabla) v] = F", "conditions": "Inviscid force accounting: pressure, New
tonian gravity and an external gradient force; their S9b counterparts are not supplied."}}
```

```
$ python3 -I -c 'import json,sys; print(json.dumps({sys.argv[2]: json.load(open(sys.argv[1]))[sys.argv[2]]}))' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/oracle_quotes.json V7
{"V7": {"citation": "Visser, same local copy, p.4, section 2.1, equation (7)", "quoted": "F = -nabla p - rho nabla phi - rho nabla Phi", "conditions": "phi is Newtonian potential; Phi is externally im
posed driving potential. No S9b force-law comparison is computed."}}
```

```
$ python3 -I -c 'import json,sys; print(json.dumps({sys.argv[2]: json.load(open(sys.argv[1]))[sys.argv[2]]}))' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/oracle_quotes.json V7background
{"V7background": {"citation": "Visser, same local copy, p.15, section 7, unnumbered background equations immediately preceding (53)", "quoted": "c = 1; v = sqrt(2GM/r); rho |v| proportional to 1/r^2; 
rho proportional to r^(-3/2); c^2 = dp/drho; p = p_infinity + c^2 rho; p-p_infinity proportional to c^2 r^(-3/2)", "conditions": "c is position independent. For constant rho and c with that velocity, 
the paper says: this assignment is incompatible with the continuity equation. Its displayed condition is nabla dot (rho v) != 0. Density-law normalizations remain unspecified. Neither density law is s
ubstituted into S9b's exchange."}}
```

```
$ python3 -I -c 'import json,sys; print(json.dumps({sys.argv[2]: json.load(open(sys.argv[1]))[sys.argv[2]]}))' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/oracle_quotes.json V54
{"V54": {"citation": "Visser, same local copy, p.16, section 8, equation (54)", "quoted": "v = c r_0^2/r^2", "conditions": "Incompressible, spherically symmetric fluid. rho is position independent; th
e barotropic assumption makes pressure and c position independent. Continuity gives v proportional to 1/r^2. For the bubble realization use the exterior of the bubble and the incompressible approximat
ion."}}
```

```
$ python3 -I -c 'import json,sys; print(json.dumps({sys.argv[2]: json.load(open(sys.argv[1]))[sys.argv[2]]}))' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/oracle_quotes.json HL1
{"HL1": {"citation": "Hamilton-Lisle gr-qc/0411060v2 (31 Aug 2006), local sources/gr-qc_0411060.pdf, p.1, section I, equation (1)", "quoted": "ds^2 = -dt_ff^2 + (dr + beta dt_ff)^2 + r^2 (dtheta^2 + s
in^2 theta dphi^2)", "conditions": "Schwarzschild river; beta is inward speed in light-speed units; t_ff is proper time of radial infall from rest at infinity. Overall factor 1."}}
```

```
$ python3 -I -c 'import json,sys; print(json.dumps({sys.argv[2]: json.load(open(sys.argv[1]))[sys.argv[2]]}))' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/oracle_quotes.json HL2
{"HL2": {"citation": "Hamilton-Lisle, same local copy, p.1, section I, equation (2)", "quoted": "beta = (2GM/r)^(1/2)", "conditions": "Spherical Schwarzschild mass M; beta in light-speed units. G M he
re is an oracle parameter, not S9b's orbital GM."}}
```

```
$ python3 -I -c 'import json,sys; print(json.dumps({sys.argv[2]: json.load(open(sys.argv[1]))[sys.argv[2]]}))' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/oracle_quotes.json HL4
{"HL4": {"citation": "Hamilton-Lisle, same local copy, p.3, section II.A, equation (4)", "quoted": "beta(r) = [2M(r)/r]^(1/2)", "conditions": "Spherical mass-function metric (3); units c = G = 1. M(r)
 is live, not assumed constant."}}
```

### C4b. Visser, local copy, page 15: the section-7 background and its conformal factor
```
$ pdftotext -f 15 -l 15 -layout /var/projects/toy_physics/_scratch/s9b_oracle_check/sources/gr-qc_9712010.pdf - | sed -n '15,20p;39,47p'
      Given the Painlevé–Gullstrand line element, it might seem trivial to force
the acoustic
       p       metric into this form: simply take ρ and c to be constants, and set
v = 2GM/r? While this certainly forces the acoustic metric into the Painlevé–
Gullstrand form the problem with this is that this assignment is incompatible
with the continuity equation ∇   ~ · (ρ~v ) 6= 0 that was used in deriving the acoustic
The net result is conformal to the Painlevé–Gullstrand form of the Schwarzschild
geometry but not identical to it. For many purposes this is quite good enough:
we have an event horizon, we can define surface gravity, we can analyze Hawking
radiation. Since surface gravity and Hawking temperature are conformal invari-
ants [34] this is sufficient for analyzing basic features of the Hawking radiation
process. The only way in which the conformal factor can influence the Hawk-
ing radiation is through backscattering off the acoustic metric. (The phonons
are minimally coupled scalars, not conformally coupled scalars so there will in
general be effects on the frequency-dependent greybody factors.)
```

## C5. The check's identification notes (every UNMATCHED row) and its controls record (opening paragraphs)
```
$ grep -c '^| `' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_identifications.md
18
```

```
$ grep '^| `' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_identifications.md | cut -d'|' -f2
 `velocity_and_frame_HL` 
 `velocity_and_frame_Visser` 
 `free_fall_versus_lab_time` 
 `units_c0_cgamma_oracle_speed` 
 `mass_calibration_HL` 
 `mass_calibration_Visser` 
 `background_continuity_Visser` 
 `force_accounting_Visser` 
 `barotropic_Visser` 
 `inviscid_Visser` 
 `irrotational_Visser` 
 `section7_EOS_and_density` 
 `section8_EOS_and_density` 
 `scalar_acoustics_vs_transverse_optics` 
 `no_shear_bearing_sheet` 
 `no_off_sheet_exchange` 
 `PY_dispersion_chi_A_U` 
 `WL_dispersion_aMetric_vLocal` 
```

```
$ sed -n '26p' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_identifications.md
Other unmatched comparisons: Hamilton–Lisle supplies no material continuity/exchange operand; S9b prints no line element. `spatial_dimension` remains `SOURCE_FIXED`.
```

```
$ sed -n '3,5p' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_controls.md
**Visser FORM sensitivity is reported as a change in computed values.** The value checker reads the baseline and control residuals, prints both operands and their domains, computes their difference, then prints its exact zero test under that same printed domain. It does not count textual inequality.

Of **18** balance comparisons, **10** computed deltas print exactly `0` and their zero tests print `True`. In the other **8** (cases 1 and 2 of each engine's deflection and radar families), the computed delta vanishes exactly when the live density's radial derivative at the comparison radius is zero. Thus the FORM change moves those eight values away from that derivative-zero locus. This is a conditional value-sensitivity result; no constant-density premis
```

```
$ grep -n 'The Hamilton–Lisle FORM control is dropped' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_controls.md
16:**The Hamilton–Lisle FORM control is dropped.** With exact zero tests it did not complete at 50 GiB in `run_023` and `run_025`. These are historical commands for the pre-removal code, not commands supported by the final driver:
```

```
$ grep -n 'optical comparison has no FORM control' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_review_r2_disposition.md
66:- **The optical comparison has no FORM control in the accepted version.** This follows from the user's choice, and
```

## C6. The S11 departure and its owner
```
$ sed -n '149,159p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S11_stray_longitudinal.md

| question | present owner |
|---|---|
| where the free homogeneous longitudinal root lies relative to simple bulk sound | S11 kinematic `KW_ZERO_LOCUS` |
| whether that mode is bound, a threshold state, a bound state in the continuum, or a resonance | S11b interface and full-slab spectrum |
| whether matter directly observes the second cone | derived matter/interface coupling |
| whether the mode participates in a finite-amplitude photon helper | nonlinear intensity coupling plus harmonic and sideband response |

At homogeneous linear order, the mode's bulk leakage and direct observability require the brane–bulk
interface law. Its possible finite-amplitude photon-helper role additionally requires nonlinear intensity
coupling, the complete slab spectrum, and a harmonic and sideband radiation audit.
```

```
$ sed -n '169,180p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S11_stray_longitudinal.md
- ⛔⛔ **`c_L` is NOT a gravitational wave.** Gravity in this model is *"the **FLOW** between draining
  defects — carried by the flow + Bernoulli pressure, **NOT** by ripples/radiation"*
  (`docs/conceptual_foundation.md:348`). `c_L` is an in-plane displacement of brane material.
  ⚠ Separately: the corpus contains **no mechanistic account of what a gravitational wave is**.
- ⭐⭐ **A second cone is only a departure if matter COUPLES to it.** The transverse wave equation is
  Lorentz-invariant with invariant speed `c_γ` automatically, so anything built from those modes
  inherits that invariance — and this model's matter **is** built from them (a throat held open by a
  trapped brane-shear standing wave). ⛔ The orchestrator's *"a second cone means observers see no single
  invariant speed"* was wrong as stated. ⭐ **Not circular:** the wave equation supplies the invariance,
  so a standing wave built from it must restructure under boost, with the cost diverging as `v → c_γ`.
- ⭐ **Uniform flow is unobservable; GRADIENTS are observable** — and the model needs exactly that split.
  If uniform flow were observable we would detect absolute motion; if gradients were not, there would be
```

## C7. The photon survey's thermal entries (classification rows, conditions cell)
```
$ grep -E '^\| E(32|33|34|38) \| \*\*In apparent conflict' /var/projects/toy_physics/research/pde_ledger_v3/PHOTON_POLARIZATION_SURVEY.md | cut -c1-700
| E32 | **In apparent conflict** | **Conditions:** the retained longitudinal branch couples to matter; for the thermal evidence it must also participate as an additional ordinary radiating mode under E33–E34's thermal occupation, dispersion and calibrated-coupling assumptions. **Experimental side:** “the two polarization degrees of freedom”, [James et al.](https://arxiv.org/pdf/quant-ph/0103121), and **E33's laboratory blackbody power** ([Quinn & Mar
| E33 | **In apparent conflict** | **Conditions:** the extra longitudinal branch couples to matter and participates as an **additional ordinary thermal radiation mode** at the laboratory temperature, with the usual quantum occupation/energy, dispersion and calibrated emission/absorption response used for the two-mode inference. **Experimental side:** [Quinn & Martin](https://doi.org/10.1098/rsta.1985.0058) report “\(\sigma=(5.66967\pm0.00076)\times10^{-8
| E34 | **In apparent conflict** | **Conditions:** the extra branch couples to matter and participates as an **additional ordinary thermal radiation mode** under the usual occupation, energy and dispersion assumptions; its power must enter the calibrated CMB radiance rather than cancel in the sky/reference comparison. **Experimental side:** [Fixsen et al.](https://arxiv.org/pdf/astro-ph/9605054) report RMS deviations “less than 50 parts per million of th
| E38 | **In apparent conflict** | **Conditions:** the retained extra longitudinal branch couples to matter and is cosmologically populated, remains relativistic at the relevant epochs, and contributes **additional gravitating radiation** beyond the standard photon/neutrino contents, with energy density and perturbation evolution mapped to E38's \(N_{\rm eff}\) parametrization. Its occupation, dispersion and temperature/decoupling history must give an extr
```

```
$ grep -n 'field quantization in a separate project' /var/projects/toy_physics/research/pde_ledger_v3/PHOTON_POLARIZATION_SURVEY.md | cut -c1-200
237:The current quantum boundary is explicit: **V3_STEP_PLAN.md:1198–1204** puts field quantization in a separate project [P8](#p8). Thus the quantum preparation/outcome entries below concern undefi
```

```
$ grep -n -o '.\{60\}earlier decoupling is not generically excluded' /var/projects/toy_physics/research/pde_ledger_v3/PHOTON_POLARIZATION_SURVEY.md
278:ond the fitted allowance; a sufficiently weak population or earlier decoupling is not generically excluded
```

```
$ grep -n -o '.\{60\}rather than cancel in the sky/reference comparison' /var/projects/toy_physics/research/pde_ledger_v3/PHOTON_POLARIZATION_SURVEY.md | head -2
274:sumptions; its power must enter the calibrated CMB radiance rather than cancel in the sky/reference comparison
```

```
$ sed -n '110,111p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S11_stray_longitudinal.md
| transverse | `ω² = (μ_R/ρ_br)k²`, nullity `D−1`, kernel ⊥ `k` ⭐ **unchanged from S9/S10** |
| longitudinal | `ω² = (B_comp/ρ_br)k²`, nullity `1`, kernel ∥ `k` |
```

```
$ sed -n '1198,1204p' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
| **a bound state** — a localized mode with a **discrete** frequency spectrum, trapped by a defect | ⭐ a linear wave equation in a **varying** medium (Sturm–Liouville) | ⭐⭐ **reachable, LINEAR.** No `ħ`, no nonlinearity |
| **a soliton** — one allowed lump, energy fixed by a size–amplitude relation | the nonlinear shear action | ⛔ **C6** — the real wall, and what the user actually wants |
| **`ħω` quanta** — energy in discrete lumps | field quantization | ⛔ **not a classical PDE sim at all**, at any effort |

⭐⭐ **⇒ Shelving `ħ_model` does NOT block the simulation.** A classical field has continuous energy no
matter how good the code is; row three was never reachable this way. ⇒ `ħ_model` connects the model to
quantum mechanics — a **separate project** from watching the brane move. ⛔ Do not re-open it here.
```

## C8. The owners named in R-S22-01: S12, Q2, S22 headings and S22's adjacent conversion tests
```
$ grep -n '^### S12 \|^### Q2 \|^### S22 ' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
579:### S12 · The throat's non-variational SOURCE, plus its boundary data
896:### Q2 · The puncture deflection and its source ⛔ the holder is a DEBT, not a result
1170:### S22 · ⭐⭐ HALF TWO — what only a simulation can settle
```

```
$ grep -n 'one nonlinear throat solve' /var/projects/toy_physics/research/pde_ledger_v3/CHARTER.md | head -1
147:interior debt list** (**S22**), because S22's organising spine is already *"one nonlinear throat solve"*
```

```
$ grep -n '^| 1 |\|^| 2 |' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md | cut -c1-220
1180:| 1 | ⭐⭐ **transverse → longitudinal conversion at a defect** — how much light converts to the stray mode passing a particle. ⚠ The families do **not** mix in a *homogeneous* brane (coupling `∝ k·a`); `
1181:| 2 | ⭐⭐ **does the longitudinal actually radiate into the bulk** | the interface law | **S11b** |
```

```
$ sed -n '508p' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md | cut -c1-200
S11c asked how much light converts into other motion where the brane changes. It has conditional operator
```

```
$ grep -n 'uniform light decouples from the bulk' /var/projects/toy_physics/STATUS.md | cut -c1-160
35:- **CONDITIONAL:** uniform light decouples from the bulk in the supplied linear model. SymPy and independent
```

## C9. Was the thermal conflict captured downstream? (before this packet: 8d91f448; after: working tree)
```
$ for f in SUBSTRATE_REQUIREMENTS.md V3_STEP_PLAN.md DEFECT_REGISTER.md; do printf '%s at 8d91f448: ' $f; git -C /var/projects/toy_physics show 8d91f448:research/pde_ledger_v3/$f | grep -c -i 'thermal\|blackbody\|Stefan\|N_eff'; done
SUBSTRATE_REQUIREMENTS.md at 8d91f448: 0
V3_STEP_PLAN.md at 8d91f448: 0
DEFECT_REGISTER.md at 8d91f448: 0
```

```
$ grep -c -i 'thermal\|blackbody\|Stefan\|N_eff' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
4
```

## C10. Register: acceptance commits of the earlier passes, and the census before and after
```
$ git -C /var/projects/toy_physics log --format='%h %ad %s' --date=short -1 22be0137
22be0137 2026-10-05 Accept S11 requirements pass 2 (repair 4b; both legs clear)
```

```
$ git -C /var/projects/toy_physics log --format='%h %ad %s' --date=short -1 72866fcf
72866fcf 2026-10-08 Accept the O2 record r4 (sub-step 7 closed; O2 closed)
```

```
$ git -C /var/projects/toy_physics log --format='%h %ad %s' --date=short -1 ff3bac63
ff3bac63 2026-10-19 Accept the S9b record and its register pass (review clear)
```

```
$ git -C /var/projects/toy_physics show 8d91f448:research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | grep -c '^### R-'
17
```

```
$ grep -c '^### R-' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
18
```

```
$ grep '^### R-' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | cut -c1-90
### R-S6-01 — `B_comp` must be retired or re-affirmed
### R-S7-01 — the slab-width flat direction, under a non-uniform width
### R-S1-01 — the brane's spatial dimension
### R-S1-02 — the substructure's shear response, as a function of the order parameter
### R-S1.5-01 — the GNLS's own mode content
### R-S6-02 — the reference state must be an equilibrium
### R-S8-01 — the form of the brane's quadratic stiffness functional
### R-S8-02 — the in-plane and out-of-plane sectors must decouple at quadratic order
### R-S8-03 — the SIGN of the physical transverse stiffness
### R-S8-04 — what carries the brane's internal angular momentum
### R-S8-05 — the frame the brane's rotational stiffness is measured against
### R-S1-03 — the substructure's microscopic time-reversibility
### R-S8-06 — material displacement and the slab's quadratic inertia
### R-S12-01 — the reservoir and its power budget
### R-S12-02 — background drain/return and separate boundary data
### R-O2-01 — the net-power supplier and its stated budget
### R-S21-01 — conditional deflection and radar-slope profile requirements
### R-S22-01 — how the energy a throat releases splits among light, the extra wave and t
```

```
$ grep -c '^| \*\*S[0-9.]*\*\*\|^| \*\*no owner' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
18
```

```
$ grep -o 'status\*\* [A-Z]*' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | sort | uniq -c
     18 status** OPEN
```

```
$ grep -n -A1 'A numerical observational bound' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
513:  that the model choose it. A numerical observational bound also needs matter-to-compression coupling;
514-  B says that is unbuilt and reports only a structural test, so no numerical-bound requirement is added.
```

## C11. Backlog runs: guard outcomes, caps and sampled peaks
```
$ grep MemTotal /proc/meminfo
MemTotal:       65779560 kB
```

```
$ sed -n '10,19p' /var/projects/toy_physics/scripts/s11c_guarded_run.py
checkouts. The registry uses that user's account-database home, independent
of $HOME, matching the scope of systemctl --user. The budget is
MemTotal rounded down to whole GiB, and each job takes one logical CPU
(distinct physical cores first, then their sibling threads). The machine is
dedicated to this work (user, 2026-10-15), so there is no host-memory reserve,
for pooled or legacy launches. A pooled job is admitted only when the registry
holds no stale, malformed or unaccounted record, a CPU and job slot are free,
the reserved caps plus its own fit the budget, and available memory covers its
cap plus every reserved job's unused cap. A live job's use is credited only
for memory the kernel cannot reclaim from it (its file cache and reclaimable
```

```
$ grep -o '"memoryMax": [0-9]*\|"startedUtc": "[^"]*"' /var/projects/toy_physics/_scratch/s11c/s11cb-89-90-full-1/invocation.json
"memoryMax": 32212254720
"startedUtc": "2026-10-18T21:52:58.175302+00:00"
```

```
$ grep -o '"exitCode": [0-9]*\|"wallSeconds": [0-9.]*\|"childOutcome": null' /var/projects/toy_physics/_scratch/s11c/s11cb-89-90-full-1/outcome.json
"exitCode": 1
"wallSeconds": 57485.648887393
"childOutcome": null
```

```
$ grep -o '"memory.peak": "[0-9]*"' /var/projects/toy_physics/_scratch/s11c/s11cb-89-90-full-1/resource-samples.jsonl | sort -t'"' -k4 -n | tail -1
"memory.peak": "32212254720"
```

```
$ tail -1 /var/projects/toy_physics/_scratch/s11c/s11cb-89-90-full-1/resource-samples.jsonl | grep -o '"utc": "[^"]*"\|"memory.current": "[0-9]*"\|"memory.events": "[^"]*"'
"utc": "2026-10-19T13:51:02.030873+00:00"
"memory.current": "32212062208"
"memory.events": "low 0\nhigh 0\nmax 43\noom 0\noom_kill 0\noom_group_kill 0\nsock_throttled 0"
```

```
$ grep -o '"memoryMax": [0-9]*\|"startedUtc": "[^"]*"' /var/projects/toy_physics/_scratch/s11c/s11cb-89b-wl-full-1/invocation.json
"memoryMax": 25769803776
"startedUtc": "2026-10-18T23:09:20.554028+00:00"
```

```
$ grep -o '"exitCode": [0-9]*\|"wallSeconds": [0-9.]*\|"childOutcome": null' /var/projects/toy_physics/_scratch/s11c/s11cb-89b-wl-full-1/outcome.json
"exitCode": 1
"wallSeconds": 4584.446006979999
"exitCode": 1
"wallSeconds": 4584.371123871999
```

```
$ grep -o '"memory.peak": "[0-9]*"' /var/projects/toy_physics/_scratch/s11c/s11cb-89b-wl-full-1/resource-samples.jsonl | sort -t'"' -k4 -n | tail -1
"memory.peak": "24155074560"
```

```
$ tail -1 /var/projects/toy_physics/_scratch/s11c/s11cb-89b-wl-full-1/resource-samples.jsonl | grep -o '"utc": "[^"]*"\|"memory.current": "[0-9]*"\|"memory.events": "[^"]*"'
"utc": "2026-10-19T00:25:45.009814+00:00"
"memory.current": "33406976"
"memory.events": "low 0\nhigh 0\nmax 0\noom 0\noom_kill 0\noom_group_kill 0\nsock_throttled 0"
```

```
$ tail -3 /var/projects/toy_physics/_scratch/s11c/s11cb-89b-wl-full-1/stdout
No more memory available.
Mathematica kernel has shut down.
Try quitting other applications and then retry.
```

```
$ grep -o '"label": "[^"]*"\|"memoryGiB": [0-9]*\|"wallSeconds": [0-9.]*\|"sampledPeakBytes": [0-9]*\|"systemdOomKill": [a-z]*' /var/projects/toy_physics-p2b/_scratch/s11c/p2b-r3-work/results.json
"label": "h-all"
"memoryGiB": 54
"wallSeconds": 588.0868241260032
"sampledPeakBytes": 55718088704
"systemdOomKill": true
"label": "a-all"
"memoryGiB": 54
"wallSeconds": 723.5162634169974
"sampledPeakBytes": 55622627328
"systemdOomKill": true
"label": "baseline"
"memoryGiB": 54
"wallSeconds": 1207.3433893829933
"sampledPeakBytes": 56083570688
"systemdOomKill": true
```

```
$ sed -n '13,23p' /var/projects/toy_physics-p2b/research/pde_ledger_v3/directives/_measurements/S11c_b_p2b_build_report.md
| Run | Entry; arguments | Cap (GiB) | LOG | Outcome | Runtime (s) | Peak memory |
|---|---|---:|---|---|---:|---|
| C | C; `--p2b` | 6 | `p2b-r1-c-controls-6-23` | exit 0 | 0.487 | 56877056 bytes |
| H | H; `--p2b` | 6 | `p2b-r1-h-controls-6-1` | exit 0 | 37.800 | 4551544832 bytes |
| A | A; `--p2b` | 6 | `p2b-r1-a-controls-6-1` | exit 0 | 62.461 | 4647460864 bytes |
| U48 | H; `--p2b-mutation unpaired` | 48 | `p2b-r2-h-all-48-1` | cgroup OOM-kill | 585.013 | ≥48 GiB |
| A48 | A; none | 48 | `p2b-r2-a-all-48-2` | cgroup OOM-kill | 723.713 | ≥48 GiB |
| B48 | B; none | 48 | `p2b-r2-baseline-48-2` | cgroup OOM-kill | 1194.346 | ≥48 GiB |
| U3 | H; `--p2b-mutation unpaired` | 54 | `p2b-r3-h-all-54-1` | cgroup OOM-kill | 588.087 | ≥54 GiB |
| A3 | A; none | 54 | `p2b-r3-a-all-54-1` | cgroup OOM-kill | 723.516 | ≥54 GiB |
| F3 | B; none | 54 | `p2b-r3-baseline-54-1` | cgroup OOM-kill | 1207.343 | ≥54 GiB |
```

```
$ for d in h-all a-all baseline; do printf '%s ' $d; grep -o '"memoryMax": [0-9]*\|"hostAvailableBytes": [0-9]*\|"poolMemoryMax": [0-9]*' /var/projects/toy_physics-p2b/_scratch/s11c/p2b-r3-$d-54-1/invocation.json | tr '\n' ' '; echo; done
h-all "memoryMax": 57982058496 "hostAvailableBytes": 60376612864 "poolMemoryMax": 66571993088 
a-all "memoryMax": 57982058496 "hostAvailableBytes": 61462528000 "poolMemoryMax": 66571993088 
baseline "memoryMax": 57982058496 "hostAvailableBytes": 61191802880 "poolMemoryMax": 66571993088 
```

```
$ git -C /var/projects/toy_physics-p2b log --format='%h %ad %s' --date=short -2
215726b 2026-10-19 Park the P2b build after repair round 3 (unreviewed, not accepted)
1131c39 2026-10-18 P2b build as reviewed in round 0 (main repo 77161062)
```

```
$ git -C /var/projects/toy_physics-p2b remote -v | wc -l
0
```

```
$ git -C /var/projects/toy_physics log --format='%h %ad %s' --date=short -1 77161062
77161062 2026-10-18 Preserve the S11c-b P2b coefficient-bridge build as reviewed (not accepted)
```

```
$ ls /var/projects/toy_physics/_scratch/s11c | grep -c -i 't7\|c1-t7\|s11cc1'
0
```

## C12. The light sector's goal and its prior art as the repository records them
```
$ sed -n '130,133p' /var/projects/toy_physics/research/pde_ledger_v3/CHARTER.md
**HALF ONE — the positive claim, and the user believes it is reachable:**

> ⭐⭐ **ONE medium supports the LINEAR part of every force — all the far-field effects of all the
> forces.**
```

```
$ sed -n '386,387p' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
   ⇒ `SUBSTRATE_REQUIREMENTS.md`. ⚠ And the algebra is **standard**: MacCullagh 1839
   ⇒ `docs/medium_requirements_and_prior_art.md`. ⛔ **Expected new: nothing** was right; ⛔ *"earned"* was
```

```
$ grep -n 'reproduces MacCullagh' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | cut -c1-200
1236:the records only capture what their authors thought to doubt — and the sector reproduces MacCullagh's
```

```
$ grep -n '^### R-S8-04\|^### R-S8-05' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
385:### R-S8-04 — what carries the brane's internal angular momentum
419:### R-S8-05 — the frame the brane's rotational stiffness is measured against
```

