# Lookups — polarization survey review, round 2 (generated 2026-10-09 19:28)

Generator: `_scratch/polarization/gen/survey_r2_lookups.sh` (sha256sum/grep/sed only).

```
$ (cd _scratch/polarization && sha256sum -c survey_review_baseline_r2.sha256)
POLARIZATION_SURVEY.md: OK
```

```
$ grep -n -o 'Verdict:\*\* [A-Z][A-Z ]*\|Verdict: [A-Z][A-Z ]*' _scratch/polarization/survey_review_r2_claude.txt _scratch/polarization/survey_review_r2_grok.txt
_scratch/polarization/survey_review_r2_claude.txt:3:Verdict: NEEDS REVISION
_scratch/polarization/survey_review_r2_grok.txt:3:Verdict: CLEAR
```

## Finding 1: E12 and the rule-5 wording
```
$ grep -n 'no model source bears on it' _scratch/polarization/survey_repair2_prompt.md
22:   material interfaces, or emission processes, and no model source bears on it.
```

```
$ grep -n -o 'no source bears on that response' _scratch/polarization/POLARIZATION_SURVEY.md
203:no source bears on that response
```

```
$ grep -n '| E1[0-2] ' _scratch/polarization/POLARIZATION_SURVEY.md | sed -n '4,6p' | cut -c1-260
218:| E10 | **Testable** | **Conditions:** the selected homogeneous, isotropic-inertia in-plane \(D=3\) linear pair applies, and **parity/chiral constitutive content and polarization-transport terms outside that retained object make no additional splitting/rot
219:| E11 | **Testable** | **Conditions:** the selected homogeneous, isotropic-inertia in-plane \(D=3\) linear pair applies, and **parity/chiral constitutive content and polarization-transport terms outside that retained object make no additional splitting/rot
220:| E12 | **Not addressed** | **Conditions:** the tested anisotropic CMB rotation spectrum, not a uniform rotation or every possible sky spectrum. No anisotropic cosmic-rotation power spectrum or CMB polarization transport result in the audited current recor
```

```
$ sed -n 1183p research/pde_ledger_v3/V3_STEP_PLAN.md
| 4 | **birefringence near a defect** — the two polarisations are degenerate *by symmetry* in a homogeneous brane; a defect splits them | a defect profile | as (1) |
```

## Finding 2: how the survey reads S11:52–64
```
$ grep -n -o 'S11_stray_longitudinal.md:52–64[^.]*' _scratch/polarization/POLARIZATION_SURVEY.md | head -3
171:S11_stray_longitudinal.md:52–64** gives “`N_SO = {D=2 → 4, D=3 → 3, D=4 → 4, D=5 → 3}`” and “`N_O = 3 for every D`”; extras are “reflection-odd”
218:S11_stray_longitudinal.md:52–64** gives a dimension/object-specific proper-rotation/reflection count; **O2_steady_brane_balance
219:S11_stray_longitudinal.md:52–64** gives a dimension/object-specific proper-rotation/reflection count; **O2_steady_brane_balance
```

```
$ sed -n 52,55p research/pde_ledger_v3/steps/S11_stray_longitudinal.md
N_SO = {D=2 → 4, D=3 → 3, D=4 → 4, D=5 → 3}        N_O = 3 for every D
```

The extras are **reflection-odd**: `(tr G)·ε^{ij}G_{ij}` at `D=2`, `ε_{ijkl}G_{ij}G_{kl}` at `D=4`.
```

```
$ grep -c 'a defect splits them' _scratch/polarization/POLARIZATION_SURVEY.md
2
```

## Finding 3: N_eff
```
$ grep -c -i 'N_eff\|N_{\\rm eff}\|N_{eff}\|effective number of' _scratch/polarization/POLARIZATION_SURVEY.md
0
```

## Finding 4: rules 3 and 4, E01 and E22
```
$ grep -n 'A fact that follows only from a supplied input' _scratch/polarization/survey_repair2_prompt.md
20:   constrains. A fact that follows only from a supplied input that states it also belongs here; name that input.
```

```
$ grep -n -o 'not a computed tomography/physical-state-selection object for rules 1 or 3' _scratch/polarization/POLARIZATION_SURVEY.md
209:not a computed tomography/physical-state-selection object for rules 1 or 3
```

```
$ sed -n 73,76p research/pde_ledger_v3/steps/S10_two_transverse_photons.md
> For the supplied curl-only in-plane action, at nonzero wavevector and away
> from allowed exceptional strata, the nonzero root has D − 1 transverse null
> directions in every measured MAIN case D = 2, 3, 4, 5. Thus the D = 3 member
> has two transverse directions.
```

## Finding 5: missed sources
```
$ sed -n 74,76p research/pde_ledger_v3/steps/S11bB_interface_assembly.md
**The transverse mode is completely decoupled on a uniform background.** The coupling is **identically
zero**, the dispersion is `ρ_br⁰ω² = μ_R k²`, and the imaginary part is **zero**. Both engines, same
structural reason: in-plane parity admits no `e_W ↔ u_T` bilinear.
```

```
$ sed -n 425p research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
  or `μ_W`, and does not identify the thickness mode with S10's out-of-plane displacement. The cited
```

```
$ grep -c 'S11bB_interface_assembly.md:74\|S11bB:74' _scratch/polarization/POLARIZATION_SURVEY.md
0
```

```
$ grep -c 'SUBSTRATE_REQUIREMENTS.md:425\|R-S8-06' _scratch/polarization/POLARIZATION_SURVEY.md
0
```

