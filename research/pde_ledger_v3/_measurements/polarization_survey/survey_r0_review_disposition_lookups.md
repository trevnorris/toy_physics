# Lookups — polarization survey review, round 0 (generated 2026-10-09 18:34)

Generator: `_scratch/polarization/gen/survey_r0_lookups.sh` (sha256sum/grep/sed only).

## The reviewed file is unchanged
```
$ sha256sum -c _scratch/polarization/survey_review_baseline_r0.sha256
sha256sum: POLARIZATION_SURVEY.md: No such file or directory
POLARIZATION_SURVEY.md: FAILED open or read
sha256sum: WARNING: 1 listed file could not be read
```

## Verdicts
```
$ grep -n -o 'Verdict:\*\* [A-Z][A-Z ]*\|Verdict: [A-Z][A-Z ]*' _scratch/polarization/survey_review_r0_claude.txt _scratch/polarization/survey_review_r0_grok.txt
_scratch/polarization/survey_review_r0_claude.txt:1:Verdict: NEEDS REVISION
_scratch/polarization/survey_review_r0_grok.txt:1:Verdict: NEEDS REVISION
```

## E01 and the survey's conflict rule (C1, G1)
```
$ grep -n '| E01' _scratch/polarization/POLARIZATION_SURVEY.md
71:| E01 | Polarization tomography using independent linear and circular analyzer settings | James, Kwiat, Munro & White, 2001: reconstruction in a two-dimensional single-photon polarization space; four Stokes quantities including intensity. Their two-photon example uses 16 projections. **No univers
183:| E01 | **Reproduced** | **S10**, conditionally: two exactly transverse in-plane directions on the supplied homogeneous, isotropic-inertia D=3 branch, **research/pde_ledger_v3/steps/S10_two_transverse_photons.md:73–76,182–197** [R10](#r10). This reproduces the restricted classical count, not
```

```
$ grep -n 'Apparent conflict requires' _scratch/polarization/POLARIZATION_SURVEY.md
179:“No mechanism” and “not addressed” refer to the audited current v3 results, with the cited exploratory leads where relevant. They do not claim that every historical file was silent. Apparent conflict requires an experimentally comparable model prediction. The unresolved physical coupling
```

```
$ grep -n 'apparent conflict' _scratch/polarization/survey_prompt_r0.md
53:   - in apparent conflict, quoting both sides;
```

```
$ sed -n 249,251p research/pde_ledger_v3/steps/S10_two_transverse_photons.md
directions remain in the spectrum. Within the light sector, Maxwell has no
counterpart to this surviving zero-frequency longitudinal direction. That is
a characterised departure, and its onward disposition belongs to S11
```

```
$ sed -n 172,176p research/pde_ledger_v3/steps/S10_two_transverse_photons.md
The physical selection D = 3 is not made in S10. The live S10 computation keeps
D symbolic for dimensions and evaluates an indexed sweep at D = 2, 3, 4, 5.
The new Lean baseline proof establishes the conditional map D ↦ D − 1 for
arbitrary finite D with a nonzero wavevector, extending the measured sweep.
Neither establishes which D nature selects.
```

```
$ ls research/pde_ledger_v3/steps | grep -i 's11'
S11bA_interface_response.md
S11bA_PREREG_ADDENDUM2_TAU.md
S11bA_PREREGISTERED_PREDICTION.md
S11bB_interface_assembly.md
S11b_HANDOFF.md
S11b_interface_coupling_law.md
S11b_PREREGISTERED_PREDICTION.md
S11b_RUN_CHECKLIST.md
S11b_wl_engine_review_disposition.md
S11c_a_interface_shape_derivatives.md
S11c_b_variable_coefficient_operator.md
S11c_c1_curved_bulk_closure.md
S11c_c2_self_energy_fold.md
S11c_d_profile_conditioned_scattering.md
S11c_PARTIAL_CLOSEOUT.md
S11c_SCOPE.md
S11_PREREGISTERED_PREDICTION.md
S11_stray_longitudinal.md
```

```
$ grep -n 'Lifts S10' research/pde_ledger_v3/steps/*.md
research/pde_ledger_v3/steps/S11_stray_longitudinal.md:15:Lifts S10's zero to a propagating longitudinal mode, enters the brane's compression modulus as a
```

```
$ grep -n 'only a departure if matter' research/pde_ledger_v3/steps/*.md
research/pde_ledger_v3/steps/S11_stray_longitudinal.md:173:- ⭐⭐ **A second cone is only a departure if matter COUPLES to it.** The transverse wave equation is
```

```
$ sed -n 128,136p research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
### R-S1-01 — the brane's spatial dimension

- **source** S10 (`steps/S10_two_transverse_photons.md`); S11
  (`steps/S11_stray_longitudinal.md`, moves 2–3 / finite census) · **target** S6 · **status** OPEN
- **requirement** — the brane's spatial dimension `D_brane`, as a derived quantity or an explicitly
  re-affirmed postulate.
- **on failure** — S10's headline reads *"light having exactly two polarisations is a statement that our
  space is three-dimensional."* Read backwards it says the opposite: `D_brane = 3` went in and `D−1 = 2`
  came out. Without a delivered `D_brane` the sentence is an assumption restated, ⛔ not a result.
```

## One standard across E10, E11, E13–E15 (C2)
```
$ grep -n '| E1[0-5]' _scratch/polarization/POLARIZATION_SURVEY.md
80:| E10 | Survival of gamma-ray linear polarization over a cosmological distance, testing energy-dependent helicity birefringence | Götz et al., 2014, GRB 140206A: redshift \(z=2.739\); polarization fraction above 28% at 90% confidence; dimension-five Lorentz-violating birefringence parameter boun
81:| E11 | Wavelength-dependent polarization changes from distant optical sources, testing vacuum Lorentz violation | Kostelecký & Mewes, 2002: birefringent dimension-four Standard-Model Extension coefficient combinations constrained at the \(2\times10^{-32}\) level. | A model-specific bound on bir
82:| E12 | Spatially varying rotation of CMB polarization | Bianchini et al., SPTpol, 2020: 500 square degrees at 150 GHz; scale-invariant rotation-power amplitude \(L(L+1)C_L^{\alpha\alpha}/(2\pi)<0.10\times10^{-4}\ \mathrm{rad}^2=0.033\ \mathrm{deg}^2\), 95% confidence. | Null result for **anisotr
83:| E13 | Uniform CMB polarization rotation, separating instrumental angle from foreground emission | Minami & Komatsu, Planck, 2020: \(\beta=0.35\pm0.14^\circ\), 68% confidence; 2.4σ preference. | **Reported hint, not established cosmic birefringence.** Its evidence is a nonzero fitted rotation u
84:| E14 | Uniform rotation with newer Planck PR4 maps and foreground/mask tests | Diego-Palazuelos et al., 2022: initial nearly full-sky \(\beta=0.30\pm0.11^\circ\), 68% confidence; the estimate decreases with more restrictive masks. | **Reported hint, not established.** The authors explicitly with
85:| E15 | Uniform cosmic rotation in ACT DR6, including calibration/systematic assessment | Diego-Palazuelos & Komatsu, 2025 preprint, **revised 2026-04-14, v2**: \(\beta=0.215\pm0.074^\circ\), 68% confidence, 2.9σ preference. | **Reported hint, not established.** Unknown instrumental systematics 
192:| E10 | **A place where the model could make a testable statement** | Homogeneous transverse degeneracy is present in the supplied S10 spectrum; nonuniform birefringence is a named **OPEN** future item, **S10:73–86** [R10](#r10), **V3_STEP_PLAN.md:1174–1184** [R19](#r19), **S11c_PARTIAL_CLOS
193:| E11 | **A place where the model could make a testable statement** | The same conditional common transverse root and open propagation-splitting item apply; **S10:73–86** [R10](#r10), **V3_STEP_PLAN.md:1183** [R19](#r19). The model has no supplied mapping to the measured SME coefficients; its 
194:| E12 | **Not addressed by the model** | No anisotropic cosmic-rotation power spectrum or CMB polarization transport result in the audited current records. S22 names defect birefringence, not a CMB \(C_L^{\alpha\alpha}\): **V3_STEP_PLAN.md:1174–1184** [R19](#r19); optical transport is not supp
195:| E13 | **Not addressed by the model** | No model prediction of a uniform cosmic rotation angle in the cited records; same current optical scope as E12 [R17](#r17), [R19](#r19). The tentative nonzero Planck fit is not a required established new effect. |
196:| E14 | **Not addressed by the model** | No foreground/calibration-dependent Planck rotation prediction; current birefringence item is defect-scoped and open, **V3_STEP_PLAN.md:1174–1184** [R19](#r19). A shared-data tentative signal does not create a reproduced or excluded brane prediction. |
197:| E15 | **Not addressed by the model** | No uniform ACT-band cosmic rotation law is delivered by **S11c_PARTIAL_CLOSEOUT.md:31** [R9](#r9) or **O2:125–127** [R17](#r17). The revised hint remains an observational possibility without a current model counterpart. |
```

```
$ sed -n 1183p research/pde_ledger_v3/V3_STEP_PLAN.md
| 4 | **birefringence near a defect** — the two polarisations are degenerate *by symmetry* in a homogeneous brane; a defect splits them | a defect profile | as (1) |
```

```
$ grep -n 'residual w^2_(1) - w^2_(2)\|Stokes residual' _scratch/polarization/survey_review_r0_claude_evidence/degenerate_transverse_polarization.stdout
8:residual w^2_(1) - w^2_(2) = 0
13:  Stokes residual out-in = [0, 0, 0, 0]
20:residual w^2_(1) - w^2_(2) = 45*k**2*mu*(s - 1)/(49*rho*s)
25:  Stokes residual out-in = [0, 0, -0.407159377378424, -0.514858388023709]
```

## The w-displacement in the plan (C3)
```
$ sed -n 453,455p research/pde_ledger_v3/V3_STEP_PLAN.md
⛔⛔ **DIFFERENT OBJECTS — ⛔ do not weld them, at any strength.** `u_L` is an **in-plane longitudinal
brane displacement**; **`±w`** is the **throat's normal/orientation direction**, and the charge the
rethink produced is carried by the **`h`-branon**. The committed sources keep them apart:
```

```
$ sed -n 460,461p research/pde_ledger_v3/V3_STEP_PLAN.md
  *"It is NOT the committed charge scalar (that's the `h`-branon, `h≠u_L`; `h` remains the committed
  mediator for the conditional `1/R²` falloff); `u_L` is a separate charge-odd density mode BC'd to
```

```
$ grep -n 'ξ_w = ℓh\|THROAT_H_SOURCE_1_OVER_R2' research/pde_ledger_v3/V3_STEP_PLAN.md | head -5
897:A **±w puncture geometrically bends the brane into ±w**: the field identity `ξ_w = ℓh` and the
898:orientation-odd mouth source. Token: **`THROAT_H_SOURCE_1_OVER_R2`**.
```

```
$ grep -n 'w-directed normal displacement' _scratch/polarization/POLARIZATION_SURVEY.md
115:For the requested geometry, the distinction is explicit: S10 supplies a D-component **in-plane** displacement and inherits its separation from out-of-plane fields; S11b supplies a finite thickness W along w; R-S8-02 asks for the full in-plane u plus normal h operator. With the stipulated three-d
```

```
$ grep -n '| E02' _scratch/polarization/POLARIZATION_SURVEY.md
72:| E02 | Departure from Coulomb's law, interpreted as photon mass with Proca equations | Williams, Faller & Hill, 1971: \(\mu^2=(1.04\pm1.2)\times10^{-19}\ \mathrm{cm}^{-2}\); equivalently force exponent \(q=(2.7\pm3.1)\times10^{-16}\) in \(r^{-(2+q)}\). Here \(\mu=m_\gamma c/\hbar\). | A laborato
184:| E02 | **A place where the model could make a testable statement** | S10 supplies a gapless transverse root, **S10_two_transverse_photons.md:237–252** [R10](#r10); that is not yet a Proca/Coulomb-law observable. The native Gauss-constraint status and extra electric-carrier issue remain unclea
```

## S9b's supplied polarization identification and S11c-d (C4)
```
$ git show ede8aa21:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md | sed -n 53,55p
  remain live. The same `ρ_br(x)` appears in the mass balance below. There is one isotropic speed, the same
  for both polarizations, measured relative to the shear-carrying material; the supplied polarization
  identification is `c_γ,1(x) ≡ c_γ,2(x) ≡ c_γ(x)`.
```

```
$ git show ede8aa21:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md | sed -n 106,109p
**Outside this step.** These were narrowed out with the user's approval, and are recorded as open:
- direction-dependent (radial versus tangential) stiffness;
- coupling to thickness or bulk fields;
- polarization-dependent propagation;
```

```
$ grep -c 'S9b_SHARED_PHYSICS' _scratch/polarization/POLARIZATION_SURVEY.md
0
```

```
$ ls research/pde_ledger_v3/steps | grep -ci s9b
0
```

```
$ grep -n 'polarization-dependent first-order forcing' research/pde_ledger_v3/steps/*.md
research/pde_ledger_v3/steps/S11c_d_profile_conditioned_scattering.md:58:**CONDITIONAL:** finite packet/source/matching evidence includes an inner-kernel bank, local Gaussian actions, polarization-dependent first-order forcing and matched transverse amplitudes. The handoff records positive end-curre
```

## Thermal mode count absent (C5)
```
$ grep -n -i 'stefan\|blackbody\|black-body\|N_eff\|Planck 2018' _scratch/polarization/POLARIZATION_SURVEY.md | head
1720:   149	| FD / CMB blackbody | Fixsen et al., [ApJ 473, 576 (1996), primary abstract](https://arxiv.org/abs/astro-ph/9605054): FIRAS RMS deviations **less than 50 ppm of CMB peak**; **`\|y\|<15×10⁻⁶`**, **`\|μ\|<9×10⁻⁵`**, **95% CL**. Bounds spectral distortion, not universal loss rat
1722:   151	| TL — Tolman brightness | Lubin & Sandage, [AJ 122, 1084 (2001), primary abstract](https://arxiv.org/abs/astro-ph/0106566): 34 early-type galaxies, clusters **`z=0.76,0.90,0.92`**; brightness exponent **`n=2.59±0.17` (R), `3.37±0.13` (I)** for **`q₀=1/2`**. With luminosity-evoluti
```

## Part 1 on why two states (C6)
```
$ grep -n 'does not mean three freely propagating' _scratch/polarization/POLARIZATION_SURVEY.md
48:**Spin angular momentum** is associated with polarization. An ideal circular plane/paraxial wave has angular momentum along its travel direction divided by energy equal to \(\pm1/\omega\). Quantum mechanically the corresponding photon carries spin projection \(\pm\hbar\); its **helicity** is that
```

```
$ grep -n -i 'massless' _scratch/polarization/POLARIZATION_SURVEY.md | head -5
158:| Historical normal-mode exclusion | **docs/conceptual_history.md:255–262** says the on-brane scalar uw “must stay gapped”; a phase-sector revival is “not a claim”. [R35](#r35) | **Historical hypothesis/exclusion statement**. It is not a current joint u/h spectrum. The old v2 masslessn
160:| v2 transverse count | **research/pde_ledger_v2/paper/stages/stage_003.tex:37–45**: “Transverse (earned)” and “two polarizations … physical_dof=2, massless”. [R24](#r24) | **“Earned” in the old supplied-action calculation**, an older lead. Current S10 supplies the stronger quali
161:| v2 normal/longitudinal calculation exists | **research/pde_ledger_v2/paper/stages/stage_030.tex:99–128** labels the coupled (uL,h) scalar block **“EARNED”**, with stiffness \([[B_{\rm eff},C_{hu}],[C_{hu},K_h]]\); at normalized sample values \(z_\pm=(3\pm\sqrt2)/2>0\); reduced-h massless
1558:    44	two polarizations \(\Rightarrow\) \(\mathrm{physical\_dof}=2\), massless. The
1620:   125	\paragraph{Reduced-\(h\) masslessness and conservative Hessian symmetry (EARNED).}
```

## E21 and E15 (G2, G3)
```
$ grep -n '| E21' _scratch/polarization/POLARIZATION_SURVEY.md
91:| E21 | Energy-dependent X-ray polarization of magnetar 4U 0142+61 | Taverna et al., 2022: **quoted from arXiv v1**, \(14\pm1\%\) at 2–4 keV and \(41\pm7\%\) at 5.5–8 keV, 1σ; angle changes by roughly \(90^\circ\) near 4–5 keV. | **Polarization structure observed; vacuum-birefringence inte
203:| E21 | **Not addressed by the model** | No magnetar energy-resolved mode-conversion or emission calculation in the cited current results; **S11c_PARTIAL_CLOSEOUT.md:31,37** [R9](#r9). The observed angle swing is not the model's already computed defect birefringence. |
```

```
$ grep -n '| E15' _scratch/polarization/POLARIZATION_SURVEY.md
85:| E15 | Uniform cosmic rotation in ACT DR6, including calibration/systematic assessment | Diego-Palazuelos & Komatsu, 2025 preprint, **revised 2026-04-14, v2**: \(\beta=0.215\pm0.074^\circ\), 68% confidence, 2.9σ preference. | **Reported hint, not established.** Unknown instrumental systematics 
197:| E15 | **Not addressed by the model** | No uniform ACT-band cosmic rotation law is delivered by **S11c_PARTIAL_CLOSEOUT.md:31** [R9](#r9) or **O2:125–127** [R17](#r17). The revised hint remains an observational possibility without a current model counterpart. |
```


## Correction (added 2026-10-09 19:00)
The first block above ran `sha256sum -c` from the wrong directory, so it could not open the file. The frozen copy taken before any repair carries the baseline hash, and both legs reported the hash matched at review time:
```
$ sha256sum POLARIZATION_SURVEY_reviewed_r0.md; cat survey_review_baseline_r0.sha256
d423244eac6f3c26981bbe3244c43c2d96205bfecda75b04bcead4eb1fa6913c  POLARIZATION_SURVEY_reviewed_r0.md
d423244eac6f3c26981bbe3244c43c2d96205bfecda75b04bcead4eb1fa6913c  POLARIZATION_SURVEY.md
$ grep -n -o "matches .survey_review_baseline_r0.sha256.\|matches the supplied baseline\|sha256 (.d423244e….) matches" survey_review_r0_grok.txt survey_review_r0_claude.txt
survey_review_r0_claude.txt:3:sha256 (`d423244e…`) matches
survey_review_r0_grok.txt:3:matches `survey_review_baseline_r0.sha256`
```
