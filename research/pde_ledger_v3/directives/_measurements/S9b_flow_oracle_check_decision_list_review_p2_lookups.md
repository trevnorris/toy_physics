# Measurements: review of the S9b flow-only oracle-check decision list, pass 2 (generated 2026-10-19 08:13)

Generator: `_scratch/s9b_oracle/gen/dl_review_p2_lookups.sh` (sha256sum/grep/sed/cut/wc only). The manuscripts' text
is the pass-1 Codex leg's download, hard-linked at `_scratch/s9b_oracle/dl_review_codex_evidence/`.

## The reviewed list and the prompt
```
$ (cd /var/projects/toy_physics && sha256sum -c /var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S9b_flow_oracle_check_decision_list_review_baseline_p2.sha256)
research/pde_ledger_v3/directives/S9b_flow_oracle_check_decision_list.md: OK
```

```
$ sha256sum /var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S9b_flow_oracle_check_decision_list_review_prompt_p2.md | cut -c1-64
92cc1a1646e43e3b3eca3e62ca48c9495789241f79515d981f7b1da58d85f936
```

## Verdicts
```
$ grep -o 'Verdict:[^A-Z]*[A-Z][A-Z ]*' /var/projects/toy_physics/_scratch/s9b_oracle/dl_review_p2_codex_evidence/review.md | head -1
Verdict: SOUND
```

```
$ grep -o 'Verdict:[^A-Z]*[A-Z][A-Z ]*' /var/projects/toy_physics/_scratch/s9b_oracle/dl_review_p2_grok.txt | head -1
Verdict: NEEDS REVISION
```

## Grok 1: the S9b optical operands as printed
```
$ sed -n '37,39p;79,85p' /var/projects/toy_physics/research/pde_ledger_v3/directives/S9b_flow_oracle_check_decision_list.md
  - **The supplied optical input.** The leading-eikonal dispersion `(ω − Vⁱkᵢ)² = cγ² gⁱʲkᵢkⱼ` (record lines
    68–73), as each engine prints it: `PY_LOCAL_S9B_VECTOR_DISPERSION`; `WL_LOCAL_S9B_SUPPLIED_INPUTS` (its metric)
    and `WL_LOCAL_S9B_FERMAT_OBJECT`. It is a supplied input of S9b, not an S9b result.
3. **Optical characteristics.**
   - *S9b operand:* the supplied leading-eikonal dispersion on the flow-only restriction. This comparison checks a
     supplied input against the oracle, not an S9b result, and its output says so.
   - *Oracle operands:* the characteristic cone of Hamilton–Lisle's river line element and of Visser's acoustic line
     element (§2, with §10's eikonal form), from the line elements as the papers state them.
   - *Residual:* the difference of the characteristic cones at leading eikonal order, with the overall conformal
     factor printed explicitly. S9b prints no line element, so a full line-element residual is `UNMATCHED`.
```

```
$ grep '^PY_LOCAL_S9B_VECTOR_DISPERSION' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out | cut -c1-300
PY_LOCAL_S9B_VECTOR_DISPERSION: Add(Mul(Integer(-1), Symbol('s9b_chi', real=True), Add(Pow(Symbol('s9b_k_t', real=True), Integer(2)), Pow(Symbol('s9b_k_z', real=True), Integer(2)), Mul(Pow(Symbol('s9b_A', real=True), Integer(-1)), Pow(Symbol('s9b_k_r', real=True), Integer(2))))), Pow(Add(Mul(Integer
```

```
$ grep '^PY_LOCAL_S9B_VECTOR_DISPERSION' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out | grep -o "Symbol('s9b_[a-zA-Z_]*'" | sort -u
Symbol('s9b_A'
Symbol('s9b_chi'
Symbol('s9b_k_r'
Symbol('s9b_k_t'
Symbol('s9b_k_z'
Symbol('s9b_omega'
Symbol('s9b_U'
```

```
$ grep '^WL_LOCAL_S9B_SUPPLIED_INPUTS' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out | grep -o '"Dispersion" -> [^,]*,[^,]*' | head -1
"Dispersion" -> -(cSquared*(kr^2/aMetric + kphi^2/r^2)) + (omega - kr*vLocal)^2 == 0, "Speed" -> cSquared == muPerp[r]/rhoBr[r]
```

```
$ grep '^WL_LOCAL_S9B_SUPPLIED_INPUTS' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out | grep -o '"[A-Za-z]*" -> ' | sort | uniq -c
      1 "BackgroundAnchoring" -> 
      1 "Balance" -> 
      1 "DensityPowerDefinition" -> 
      1 "Dispersion" -> 
      1 "Embedding" -> 
      1 "FarSpeed" -> 
      1 "Metric" -> 
      1 "NonnegativeRoot" -> 
      1 "Polarizations" -> 
      1 "Pressure" -> 
      1 "References" -> 
      1 "ReferenceSpeed" -> 
      1 "Responses" -> 
      1 "SelectedAnchoring" -> 
      1 "Speed" -> 
      1 "SuppliedBulkDensityPremise" -> 
      1 "UniformAnchor" -> 
```

```
$ grep -c 'aMetric ->\|vLocal ->\|aMetric ==\|vLocal ==' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out
0
```

```
$ grep '^WL_LOCAL_S9B_FERMAT_OBJECT' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out | grep -o '"[A-Za-z]*" -> ' | sort | uniq -c
      1 "EvenMetric" -> 
      1 "Odd" -> 
      1 "TimePolynomial" -> 
      1 "TimeRoot" -> 
```

```
$ grep '^WL_LOCAL_S9B_RESTRICTION_FLOW_ONLY' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out | cut -c1-200
WL_LOCAL_S9B_RESTRICTION_FLOW_ONLY: {deltaProfile -> Function[{r}, 0], xiProfile -> Function[{r}, 0]}
```

## Grok 2: the oracles' equation numbers
```
$ grep -n '(1)$\|(2)$\|(4)$' /var/projects/toy_physics/_scratch/s9b_oracle/dl_review_codex_evidence/hamilton_lisle.txt | head -8
37:                                      ds2 = − dt2ff + (dr + β dtff )2 + r2 (dθ2 + sin2 θ dφ2 ) (1)
44:                                                         β=                                 (2)
190:can be expressed in the Gullstrand-Painlevé form (1)
195:                    β(r) =                                 (4)
```

```
$ grep -n '(4)$\|(5)$\|(6)$\|(53)$\|(54)$\|(55)$\|(85)$\|(86)$\|(87)$' /var/projects/toy_physics/_scratch/s9b_oracle/dl_review_codex_evidence/visser.txt | head -14
155:                    gµν (t, ~x) ≡  · · · · · · · · · · · · · · · · · · ·  .     (4)
179:                           ∂t ρ + ∇ · (ρ ~v ) = 0,                        (5)
183:                             ρ       ≡ ρ [∂t~v + (~v · ∇)~v ] = F~ .                     (6)
757:                                                                                  (53)
785:                                          v=c          .                             (54)
793:                                                                                     (55)
1165:                        hµν ≡  · · · · · · · · · · · · · · · · · · ·          (85)
1181:                                 hµν ∂µ ϕ ∂ν ϕ = 0.                             (86)
1194:                       ⇐⇒          − ~v0 = c.                                  (87)
```

## Grok 1: where the PY dispersion's symbols occur
```
$ grep -n "Symbol('s9b_chi'" /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out | cut -d: -f1,2 | cut -c1-80
4:PY_LOCAL_S9B_VECTOR_DISPERSION
```

```
$ grep -n "Symbol('s9b_A'" /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out | cut -d: -f1,2 | cut -c1-80
4:PY_LOCAL_S9B_VECTOR_DISPERSION
```

```
$ grep -n "Symbol('s9b_U'" /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out | cut -d: -f1,2 | cut -c1-80
4:PY_LOCAL_S9B_VECTOR_DISPERSION
5:PY_LOCAL_S9B_VECTOR_ADVECTION
```

