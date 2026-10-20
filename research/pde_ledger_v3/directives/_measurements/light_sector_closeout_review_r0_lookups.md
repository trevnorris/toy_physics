# Measurements: light-closeout review round 0 (generated 2026-10-19 18:45)

Generator: `_scratch/light_closeout/gen/review_r0_lookups.sh`. Literal retrieval only.

## D1. The round: prompt, baseline and legs' verdict lines
```
$ sha256sum /var/projects/toy_physics/_scratch/light_closeout/review_prompt_r0.md | cut -c1-64
4138835633d9b3ab7296367e203acc658a9e380e9b4003dcac88824d0f5b67c4
```

```
$ cat /var/projects/toy_physics/_scratch/light_closeout/packet_baseline_r0.sha256 | sed 's#/var/projects/toy_physics/##'
d7430108409da53d74cd59ad71e3f1ee8c93b2ee73b3d983260f52444d54ee56  research/pde_ledger_v3/steps/S9b_flow_oracle_comparison.md
2db8c4b3ffddb4db87839cad7521e73a4abc6cdb2fec94a1c5f8ab95af8cd0c6  research/pde_ledger_v3/directives/_measurements/light_sector_closeout_lookups.md
1ae4b63efffd9f52c43a00a653a311059178e860ef3ed53b656d1950e320eec0  research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_outcome_lookups.md
b42dc759d6973ddfb562bc5d79a7796b06e8ff7815846c0ce9d436741a5dd968  research/pde_ledger_v3/directives/_measurements/light_sector_reassessment_lookups.md
581e28a7181170fdec3f09cd3c397ebb2df2e54b92cb397cbf8c939e4a6e5768  STATUS.md
5491a8cede990614d5ff3c4aa1a2f50259e57957ebdeca1e3f39c60f2fbec6de  research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
9b21c3f5efaa8487e1df32b2452a48e7e45b93cc8a08157c1b21efc485c13eed  research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md
```

```
$ grep -n '^\*\*FINDINGS\*\*\|^FINDINGS\|\*\*FINDINGS\*\*' /var/projects/_review_scratch/light-closeout-r0-codex/review_report.md /var/projects/_review_scratch/light-closeout-r0-grok/grok_review.txt
/var/projects/_review_scratch/light-closeout-r0-codex/review_report.md:150:**FINDINGS**
/var/projects/_review_scratch/light-closeout-r0-grok/grok_review.txt:15:**FINDINGS**
```

```
$ grep -n '^[0-9]\. \*\*' /var/projects/_review_scratch/light-closeout-r0-codex/review_report.md /var/projects/_review_scratch/light-closeout-r0-grok/grok_review.txt | cut -c1-200
/var/projects/_review_scratch/light-closeout-r0-codex/review_report.md:152:1. **High — equilibrium alone is not the survey's falsifier (R-S22-01:783–791).**
/var/projects/_review_scratch/light-closeout-r0-codex/review_report.md:180:2. **Medium — the mass-only correspondence summary loses the check's clock, orientation and component conditions.**
/var/projects/_review_scratch/light-closeout-r0-codex/review_report.md:207:3. **Medium — S22 test 2 is a longitudinal leakage test, not an incoming-light test (R-S22-01:795–798).**
/var/projects/_review_scratch/light-closeout-r0-codex/review_report.md:222:4. **Medium — a 54-GiB resource stop does not establish a machine-wide impossibility (DEFERRED_HEAVY_RUNS:36–37).**
/var/projects/_review_scratch/light-closeout-r0-grok/grok_review.txt:17:1. **High.** The conditional correspondence is stated from the mass identification alone. The signed velocity and the frame and 
/var/projects/_review_scratch/light-closeout-r0-grok/grok_review.txt:31:2. **Medium.** The register says exchange vanishes only on the `r^{-3/2}` density. The check also prints vanishing exchange with
/var/projects/_review_scratch/light-closeout-r0-grok/grok_review.txt:39:3. **High.** `R-S22-01`’s falsifier treats equilibrium with matter as a conflict with the blackbody, the CMB, and `N_eff`. The
/var/projects/_review_scratch/light-closeout-r0-grok/grok_review.txt:65:4. **Medium.** The comparison record names the constant-`c₀` question as the object of `R-S1-02`.
```

## D2. Finding on the correspondence: signed and clock-factor zero tests
```
$ sed -n '724,729p;741p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | grep -E '^(RESIDUAL_OPERANDS|FLOW_SIGNED_RESIDUAL|RESIDUAL_ZERO_TEST)'
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FLOW_SQUARED_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "value" -> True|>
RESIDUAL_OPERANDS: {-(Sqrt[2]*Sqrt[py`GM/py`s9b$r]), Sqrt[2]*py`c$0*Sqrt[py`GM/(py`c$0^2*py`s9b$r)]}
FLOW_SIGNED_RESIDUAL: -2*Sqrt[2]*Sqrt[py`GM/py`s9b$r]
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FLOW_SIGNED_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "value" -> False|>
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FLOW_SIGNED_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/2", "value" -> True|>
```

```
$ grep -c '^RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FLOW_SIGNED_RESIDUAL".*"value" -> False|>' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out
40
```

```
$ grep -m2 'clockRate^2\*gmH|>' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | cut -c1-220
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FLOW_SQUARED_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/general", "value" -> 2*py`GM + py`s9b$r^2*py`s9b$velocity$integration$constant == 2*clockRate^2*gmH
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FLOW_SQUARED_RESIDUAL", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "value" -> py`GM == clockRate^2*gmH|>
```

```
$ grep -E '^\| `(velocity_and_frame_HL|velocity_and_frame_Visser|free_fall_versus_lab_time|units_c0_cgamma_oracle_speed)`' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_identifications.md | cut -d'|' -f2,3
 `velocity_and_frame_HL` | HL I-II describes an inward river and free-fall tetrads, whereas S9b V is the material velocity in the far-field lab frame. No physical frame-identification is supplied. 
 `velocity_and_frame_Visser` | Visser v is a Newtonian background fluid velocity; S9b V is the shear carrier velocity. Their physical identity is not supplied. Both radial orientations are retained conditionally. 
 `free_fall_versus_lab_time` | HL t_ff is infall proper time; S9b uses lab time. Visser uses Newtonian lab time, but a cross-model clock identification is not supplied. An open clock-conversion parameter is retained in the oracle Jacobian. 
 `units_c0_cgamma_oracle_speed` | S9b c0 is asymptotic light speed and cgamma is the transverse speed. HL uses light-speed units (II.A: c=G=1); Visser uses sound speed (section 7: c=1). No physical equality to the oracle speed is supplied. Raw and explicitly conditional speed-matched cones are both printed. 
```

## D3. Finding on the exchange density law: the zero-velocity cases
```
$ grep '^RESIDUAL_ZERO_TEST: <|"residual_tag" -> "EXCHANGE_RESIDUAL"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_flow_oracle_check.out | grep -o '"value" -> .*' | sort | uniq -c
      4 "value" -> 3*py`rho$br$profile[py`s9b$r] + 2*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r] == 0|>
      4 "value" -> 3*rhoBr[r] + 2*r*Derivative[1][rhoBr][r] == 0|>
     10 "value" -> True|>
```

```
$ grep -n 'vanishes only on' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
633:  - In both engines the exchange `jn` vanishes only on `ρbr = C₁/r^(3/2)`, the density law Visser §7
```

## D4. Finding on the falsifier: the survey's conditions
```
$ grep -E '^\| E(33|34|38) \| \*\*In apparent conflict' /var/projects/toy_physics/research/pde_ledger_v3/PHOTON_POLARIZATION_SURVEY.md | cut -c1-900
| E33 | **In apparent conflict** | **Conditions:** the extra longitudinal branch couples to matter and participates as an **additional ordinary thermal radiation mode** at the laboratory temperature, with the usual quantum occupation/energy, dispersion and calibrated emission/absorption response used for the two-mode inference. **Experimental side:** [Quinn & Martin](https://doi.org/10.1098/rsta.1
| E34 | **In apparent conflict** | **Conditions:** the extra branch couples to matter and participates as an **additional ordinary thermal radiation mode** under the usual occupation, energy and dispersion assumptions; its power must enter the calibrated CMB radiance rather than cancel in the sky/reference comparison. **Experimental side:** [Fixsen et al.](https://arxiv.org/pdf/astro-ph/9605054) r
| E38 | **In apparent conflict** | **Conditions:** the retained extra longitudinal branch couples to matter and is cosmologically populated, remains relativistic at the relevant epochs, and contributes **additional gravitating radiation** beyond the standard photon/neutrino contents, with energy density and perturbation evolution mapped to E38's \(N_{\rm eff}\) parametrization. Its occupation, dis
```

```
$ sed -n '141p' /var/projects/toy_physics/research/pde_ledger_v3/PHOTON_POLARIZATION_SURVEY.md | cut -c1-400
**What “thermal count” means in E33–E34.** This is an inference from the Planck-law normalization, not a detector labelling each polarization. For \(g\) equally coupled radiating states with \(\omega=ck\), quantum energy \(h_{\rm P}\nu\), and thermal occupation at the same \(T\), the spectral radiance is \(B_{\nu,g}=g h_{\rm P}\nu^3/[c^2(e^{h_{\rm P}\nu/(k_{\rm B}T)}-1)]\). The usual photon 
```

```
$ grep -n 'ω² = \|omega^2\|B_comp/ρ_br\|mu_R/rho_br' /var/projects/toy_physics/research/pde_ledger_v3/steps/S11_stray_longitudinal.md | head -6 | cut -c1-200
9:S10 left one direction of `u` with `ω² = 0`. S11 asks what happens to that direction once the brane is
110:| transverse | `ω² = (μ_R/ρ_br)k²`, nullity `D−1`, kernel ⊥ `k` ⭐ **unchanged from S9/S10** |
111:| longitudinal | `ω² = (B_comp/ρ_br)k²`, nullity `1`, kernel ∥ `k` |
276:eigenvalues (transverse `μ_R/ρ_br` nullity `D−1`, longitudinal `B_comp/ρ_br` nullity 1, no cross-modulus,
```

## D5. Finding on R-S1-02's object
```
$ grep -n -A2 '^### R-S1-02' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | cut -c1-200
164:### R-S1-02 — the substructure's shear response, as a function of the order parameter
165-
166-- **source** S9b (brane optical shear only; `steps/S9b_light_bending.md`, §§2–4/8); O2 (`steps/O2_steady_brane_balance.md`, register handoff); S9 (`steps/S9_light_requires_shear.md`) · **targ
```

```
$ grep -n 'That is the object of R-S1-02' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_flow_oracle_comparison.md
79:  depends on how the brane's shear response varies with its density. That is the object of R-S1-02, and
```

## D6. Finding on S22's test 2
```
$ sed -n '1180,1181p' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md | cut -c1-200
| 1 | ⭐⭐ **transverse → longitudinal conversion at a defect** — how much light converts to the stray mode passing a particle. ⚠ The families do **not** mix in a *homogeneous* brane (coupling
| 2 | ⭐⭐ **does the longitudinal actually radiate into the bulk** | the interface law | **S11b** |
```

## D7. Finding on the P2b memory statement
```
$ for d in h-all a-all baseline; do printf '%s ' $d; grep -o '"memoryMax": [0-9]*\|"hostAvailableBytes": [0-9]*\|"poolMemoryMax": [0-9]*' /var/projects/toy_physics-p2b/_scratch/s11c/p2b-r3-$d-54-1/invocation.json | tr '\n' ' '; echo; done
h-all "memoryMax": 57982058496 "hostAvailableBytes": 60376612864 "poolMemoryMax": 66571993088 
a-all "memoryMax": 57982058496 "hostAvailableBytes": 61462528000 "poolMemoryMax": 66571993088 
baseline "memoryMax": 57982058496 "hostAvailableBytes": 61191802880 "poolMemoryMax": 66571993088 
```

```
$ grep -n 'need more memory than this box admits' /var/projects/toy_physics/research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md
37:runs need more memory than this box admits, or a script change, which would get its own review. Every
```

## D8. Codex's note: Visser on the section-7 conformal factor (local copy, page 15)
```
$ pdftotext -f 15 -l 15 -layout /var/projects/toy_physics/_scratch/s9b_oracle_check/sources/gr-qc_9712010.pdf - | sed -n '39,47p'
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

```
$ pdftotext -f 15 -l 15 -layout /var/projects/toy_physics/_scratch/s9b_oracle_check/sources/gr-qc_9712010.pdf - | sed -n '15,20p'
      Given the Painlevé–Gullstrand line element, it might seem trivial to force
the acoustic
       p       metric into this form: simply take ρ and c to be constants, and set
v = 2GM/r? While this certainly forces the acoustic metric into the Painlevé–
Gullstrand form the problem with this is that this assignment is incompatible
with the continuity equation ∇   ~ · (ρ~v ) 6= 0 that was used in deriving the acoustic
```

