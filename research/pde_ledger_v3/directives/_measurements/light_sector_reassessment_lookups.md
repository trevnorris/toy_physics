# Measurements: light-sector reassessment (generated 2026-10-19 17:51)

Generator: `_scratch/s9b_oracle/gen/light_sector_reassessment_lookups.sh` (grep/sed/git log/wc only).

## The goal as written
```
$ sed -n '130,135p' /var/projects/toy_physics/research/pde_ledger_v3/CHARTER.md
**HALF ONE — the positive claim, and the user believes it is reachable:**

> ⭐⭐ **ONE medium supports the LINEAR part of every force — all the far-field effects of all the
> forces.**

⭐ **This is the part that can be FULLY DERIVED**, which is why it is the ledger's body. Under
```

```
$ sed -n '354,356p;361,363p' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
## PHASE 2 — light (3 steps)

### S9 · The transverse shear speed
`{ρ_br, μ_R}` are **introduced and classified in S8**, and the **R10 debt starts in S8**.
⇒ S9 is **pure consequence**: `c_γ² = μ_R/ρ_br`. **Expected new: nothing.**
⭐⭐ **This closes O-01, the "universe hole."** v2's registry seed required `c_γ` as a supplied input
```

```
$ grep -n 'Expected new: nothing\|MacCullagh 1839' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md | sed -n '1,4p'
213:`c_s² = (1/m)dP/dρ = nKρ^(n-1)/m`. **Expected new: nothing** — pure consequence.
263:**Expected new: nothing.**
362:⇒ S9 is **pure consequence**: `c_γ² = μ_R/ρ_br`. **Expected new: nothing.**
386:   ⇒ `SUBSTRATE_REQUIREMENTS.md`. ⚠ And the algebra is **standard**: MacCullagh 1839
```

```
$ sed -n '392,394p' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
### S11 · The stray-longitudinal departure
Characterized, first-class. ⛔ Not a defect to fix — a recorded departure from the reference theory.

```

## When each light-sector artifact was first added
```
$ for f in steps/S9_light_requires_shear.md steps/S10_two_transverse_photons.md steps/S11_stray_longitudinal.md steps/S11b_interface_coupling_law.md steps/S11c_PARTIAL_CLOSEOUT.md steps/O2_steady_brane_balance.md steps/S9b_light_bending.md directives/S9b_SHARED_PHYSICS.md PHOTON_POLARIZATION_SURVEY.md SUBSTRATE_REQUIREMENTS.md; do printf '%s  ' $f; git -C /var/projects/toy_physics log --diff-filter=A --format='%ad %h' --date=short -- research/pde_ledger_v3/$f | tail -1; echo; done
steps/S9_light_requires_shear.md  2026-08-01 e33bd39f

steps/S10_two_transverse_photons.md  2026-08-02 0caab021

steps/S11_stray_longitudinal.md  2026-08-02 6a6e06ba

steps/S11b_interface_coupling_law.md  2026-08-23 3177832f

steps/S11c_PARTIAL_CLOSEOUT.md  2026-10-05 bdc181c0

steps/O2_steady_brane_balance.md  2026-10-07 243a530f

steps/S9b_light_bending.md  2026-10-19 36f59ff1

directives/S9b_SHARED_PHYSICS.md  2026-10-06 c2f1cf2b

PHOTON_POLARIZATION_SURVEY.md  2026-10-09 3cd7be20

SUBSTRATE_REQUIREMENTS.md  2026-08-06 1f6581a7

```

```
$ ls /var/projects/toy_physics/research/pde_ledger_v3/steps/ | grep -E '^S9|^S10|^S11|^O2' 
O2_steady_brane_balance.md
S10_PREREGISTERED_PREDICTION.md
S10_two_transverse_photons.md
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
S9b_light_bending.md
S9_light_requires_shear.md
S9_PILOT_ADJUDICATION.md
S9_REBUILD_QUESTION_DEFECTS.md
```

## Where it stands
```
$ sed -n '1,4p' /var/projects/toy_physics/STATUS.md
# Project status — S11c PARTIAL; S12 next

The light model has a useful uniform result, but no accepted number for light lost at a defect.
S11c is closed **PARTIAL**. No nonuniform build, numerical continuation or review resumes automatically.
```

```
$ sed -n '1,7p' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
# Substrate requirements — what the force sectors oblige the unbuilt steps to deliver

**Status: passes 1 (S9, S10) and 2 (S11, S11b-A/B, S11c PARTIAL) complete; pass-2 review pending.
O2 pass populated 2026-10-07; its record/pass review is pending.
S9b pass populated 2026-10-19; its record/pass review is pending.**
Seventeen entries, all OPEN. Pass 2 was populated 2026-10-05 from the kept results and their stated
conditions; it does not close the substrate or upgrade S11c's unresolved work. Later sectors remain to
```

```
$ sed -n '6,9p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
The accepted constructions give conditional optical functionals and profile requirements, not a
substrate-derived light equation or a demonstrated match to GR. At the selected comparison grades,
M13 retrieves byte-identical paired condition payloads in PY. In WL, `GM` and `ProfileEquation`
match literally but `InverseWeight`, `ReferenceInverse` and velocity-reduction fields differ (§3).
```

## The photon survey's classification of 38 experimental facts
```
$ grep -n '^Classification totals' /var/projects/toy_physics/research/pde_ledger_v3/PHOTON_POLARIZATION_SURVEY.md
280:Classification totals: **38 entries — reproduced 1; in apparent conflict 7; testable 10; required, no mechanism yet 5; not addressed 15.** E01 reproduces the conditional in-plane transverse count, without selecting the physical dimension or completing the radiation-state census. E12 is testable by the same supplied/computed propagation object used for E10–E11 and E13–E15; the CMB source/statistical mapping 
```

```
$ grep -E '^\| E(32|33|34|38) \|' /var/projects/toy_physics/research/pde_ledger_v3/PHOTON_POLARIZATION_SURVEY.md | cut -c1-330
| E32 | Constraint on an additional longitudinal radiating polarization, separated from the in-plane count | Experimental basis is E01's tomography, E02–E03's specified-theory null tests, and **E33–E34's conditional thermal count of radiating modes**. James et al. (2001) describe “the two polarization degrees of freedom”
| E33 | Laboratory blackbody total radiant exitance and its absolute Stefan–Boltzmann normalization | Quinn & Martin, 1985: cryogenic radiometry at 273.16 K and approximately 233–373 K; reported \(\sigma=(5.66967\pm0.00076)\times10^{-8}\ \mathrm{W\,m^{-2}\,K^{-4}}\), in the measurement's original unit realization. The quoted
| E34 | Cosmological blackbody radiance spectrum of the CMB | Fixsen et al., COBE/FIRAS, 1996: \(T=2.728\pm0.004\ \mathrm K\), 95% confidence; RMS spectral deviations below **50 parts per million of the peak CMB radiance**. | **Conditional thermal count:** the spectrum fits the ordinary Planck law with its two-polarization norma
| E38 | Cosmological effective count of relativistic species, \(N_{\rm eff}\), inferred from CMB anisotropies and BAO | Planck Collaboration, **Planck 2018 results. VI. Cosmological parameters** (2020 publication; opened preprint v4, 2021): **\(N_{\rm eff}=2.99\pm0.17\)**, 68% confidence, for **TT,TE,EE+lowE+lensing+BAO**; eq. (
| E32 | **In apparent conflict** | **Conditions:** the retained longitudinal branch couples to matter; for the thermal evidence it must also participate as an additional ordinary radiating mode under E33–E34's thermal occupation, dispersion and calibrated-coupling assumptions. **Experimental side:** “the two polarization deg
| E33 | **In apparent conflict** | **Conditions:** the extra longitudinal branch couples to matter and participates as an **additional ordinary thermal radiation mode** at the laboratory temperature, with the usual quantum occupation/energy, dispersion and calibrated emission/absorption response used for the two-mode inference. 
| E34 | **In apparent conflict** | **Conditions:** the extra branch couples to matter and participates as an **additional ordinary thermal radiation mode** under the usual occupation, energy and dispersion assumptions; its power must enter the calibrated CMB radiance rather than cancel in the sky/reference comparison. **Experime
| E38 | **In apparent conflict** | **Conditions:** the retained extra longitudinal branch couples to matter and is cosmologically populated, remains relativistic at the relevant epochs, and contributes **additional gravitating radiation** beyond the standard photon/neutrino contents, with energy density and perturbation evolutio
```

```
$ grep -E '^\| E(13|14|15) \|' /var/projects/toy_physics/research/pde_ledger_v3/PHOTON_POLARIZATION_SURVEY.md | cut -c1-200
| E13 | Uniform CMB polarization rotation, separating instrumental angle from foreground emission | Minami & Komatsu, Planck, 2020: \(\beta=0.35\pm0.14^\circ\), 68% confidence; 2.4σ preference. | **R
| E14 | Uniform rotation with newer Planck PR4 maps and foreground/mask tests | Diego-Palazuelos et al., 2022: initial nearly full-sky \(\beta=0.30\pm0.11^\circ\), 68% confidence; the estimate decreas
| E15 | Uniform cosmic rotation in ACT DR6, including calibration/systematic assessment | Diego-Palazuelos & Komatsu, 2025 preprint, **revised 2026-04-14, v2**: their Bayesian estimate is \(\beta=0.21
| E13 | **In apparent conflict** | **Conditions:** the selected homogeneous, isotropic-inertia in-plane \(D=3\) linear pair applies, and **parity/chiral constitutive content and polarization-transport
| E14 | **In apparent conflict** | **Conditions:** the selected homogeneous, isotropic-inertia in-plane \(D=3\) linear pair applies, and **parity/chiral constitutive content and polarization-transport
| E15 | **In apparent conflict** | **Conditions:** the selected homogeneous, isotropic-inertia in-plane \(D=3\) linear pair applies, and **parity/chiral constitutive content and polarization-transport
```

## Is the thermal conflict captured anywhere downstream?
```
$ grep -c -i 'thermal\|blackbody\|Stefan\|N_eff' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md /var/projects/toy_physics/research/pde_ledger_v3/DEFECT_REGISTER.md
/var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:0
/var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md:0
/var/projects/toy_physics/research/pde_ledger_v3/DEFECT_REGISTER.md:0
```

## Prior art already on record for the light results
```
$ grep -n 'MacCullagh' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md | head -3
386:   ⇒ `SUBSTRATE_REQUIREMENTS.md`. ⚠ And the algebra is **standard**: MacCullagh 1839
839:⛔⛔ **HIS STATED REASON FOR LEAVING LINEAR ELASTICITY:** *"MacCullagh's theory based on linear elasticity
```

```
$ sed -n '1,4p' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_outcome_lookups.md
# Measurements: outcome of the accepted S9b flow-only oracle check (generated 2026-10-19 16:16)

Generator: `_scratch/s9b_oracle/gen/check_outcome_lookups.sh`. Transcript `scripts/out/S9b_flow_oracle_check.out`
at commit `8d91f448`. Literal retrieval and counting only; the interpretation is in the report to the user, not here.
```

## Follow-up: can later sectors decide the extra-wave question?
```
$ grep -n 'whether matter directly observes the second cone' /var/projects/toy_physics/research/pde_ledger_v3/steps/S11_stray_longitudinal.md
154:| whether matter directly observes the second cone | derived matter/interface coupling |
```

```
$ sed -n '169,176p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S11_stray_longitudinal.md
- ⛔⛔ **`c_L` is NOT a gravitational wave.** Gravity in this model is *"the **FLOW** between draining
  defects — carried by the flow + Bernoulli pressure, **NOT** by ripples/radiation"*
  (`docs/conceptual_foundation.md:348`). `c_L` is an in-plane displacement of brane material.
  ⚠ Separately: the corpus contains **no mechanistic account of what a gravitational wave is**.
- ⭐⭐ **A second cone is only a departure if matter COUPLES to it.** The transverse wave equation is
  Lorentz-invariant with invariant speed `c_γ` automatically, so anything built from those modes
  inherits that invariance — and this model's matter **is** built from them (a throat held open by a
  trapped brane-shear standing wave). ⛔ The orchestrator's *"a second cone means observers see no single
```

```
$ grep -n 'BY_TUNING\|choice the model makes\|REAL and IRREMOVABLE' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md | head -5
413:⭐⭐ **The mode is REAL and IRREMOVABLE, and its role was CATALYTIC — ⛔ not constitutive** (user's own
484:⚠ **And the Maxwell locus IS reachable** — `pathA_36` emits `C5_RESOLVED_MAXWELL_BY_TUNING`
488:> `C5_RESOLVED_MAXWELL_BY_TUNING` — so the FAIL is not hardcoded — but `BY_TUNING`, not
491:⇒ The departure is a **choice the model makes**, ⛔ not an inability.
```

```
$ sed -n '508,510p' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
- **note** — this is the records' explicit condition for retaining that channel, not a requirement
  that the model choose it. A numerical observational bound also needs matter-to-compression coupling;
  B says that is unbuilt and reports only a structural test, so no numerical-bound requirement is added.
```

## Follow-up 2: what 'coupling' and 'excited' can mean in the model
```
$ grep -n 'S11c asked how much light converts' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
508:S11c asked how much light converts into other motion where the brane changes. It has conditional operator
```

```
$ grep -n 'symmetry selection rule' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md /var/projects/toy_physics/STATUS.md
/var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md:532:**OPEN — S12:** the clean-condition packet measured a symmetry selection rule on SymPy only; its no-leak
```

```
$ grep -n '### Q2' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
896:### Q2 · The puncture deflection and its source ⛔ the holder is a DEBT, not a result
```

```
$ grep -n 'one nonlinear throat solve' /var/projects/toy_physics/research/pde_ledger_v3/CHARTER.md | head -2
147:interior debt list** (**S22**), because S22's organising spine is already *"one nonlinear throat solve"*
```

```
$ grep -o 'field quantization in a separate project' /var/projects/toy_physics/research/pde_ledger_v3/PHOTON_POLARIZATION_SURVEY.md | head -1
field quantization in a separate project
```

