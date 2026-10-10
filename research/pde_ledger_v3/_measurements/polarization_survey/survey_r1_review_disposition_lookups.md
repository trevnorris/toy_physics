# Lookups — polarization survey review, round 1 (generated 2026-10-09 19:00)

Generator: `_scratch/polarization/gen/survey_r1_lookups.sh` (sha256sum/grep/sed only).

```
$ (cd _scratch/polarization && sha256sum -c survey_review_baseline_r1.sha256)
POLARIZATION_SURVEY.md: OK
```

```
$ grep -n -o 'Verdict:\*\* [A-Z][A-Z ]*\|Verdict: [A-Z][A-Z ]*' _scratch/polarization/survey_review_r1_claude.txt _scratch/polarization/survey_review_r1_grok.txt
_scratch/polarization/survey_review_r1_claude.txt:3:Verdict: NEEDS REVISION
_scratch/polarization/survey_review_r1_grok.txt:1:Verdict: CLEAR
```

## The survey's class rule and the rows it is applied to (findings 1, 2)
```
$ sed -n 192p _scratch/polarization/POLARIZATION_SURVEY.md
The task's five classes apply to the **specific fact in each E-row**. **Reproduced** requires a step to compute the fact from inputs that do not themselves supply it; agreement inherited from a supplied count or polarization identification is labelled with that input and uses another class. **Requir
```

```
$ grep -n '| E01 \|| E1[0-7] \|| E2[5-9] \|| E3[2-4] ' _scratch/polarization/POLARIZATION_SURVEY.md | awk -F'|' 'NF>3 {print $1"|"$2"|"$3}'
73:| E01 | Polarization tomography using independent linear and circular analyzer settings 
82:| E10 | Survival of gamma-ray linear polarization over a cosmological distance, testing energy-dependent helicity birefringence 
83:| E11 | Wavelength-dependent polarization changes from distant optical sources, testing vacuum Lorentz violation 
84:| E12 | Spatially varying rotation of CMB polarization 
85:| E13 | Uniform CMB polarization rotation, separating instrumental angle from foreground emission 
86:| E14 | Uniform rotation with newer Planck PR4 maps and foreground/mask tests 
87:| E15 | Uniform cosmic rotation in ACT DR6, including calibration/systematic assessment 
88:| E16 | Rotation of linear polarization transmitted through magnetized material 
89:| E17 | Electric-field-induced double refraction in dielectric material 
97:| E25 | Spin-dependent transverse beam displacement on refraction at an air/glass interface 
98:| E26 | Polarization of light reflected from transparent surfaces 
99:| E27 | Degree and angle of scattered skylight polarization 
100:| E28 | Polarization of the cosmic microwave background produced by scattering 
101:| E29 | Vector-field focusing of a radially polarized beam, including an axial electric component 
104:| E32 | Constraint on an additional longitudinal radiating polarization, separated from the in-plane count 
105:| E33 | Laboratory blackbody total radiant exitance and its absolute Stefan–Boltzmann normalization 
106:| E34 | Cosmological blackbody radiance spectrum of the CMB 
198:| E01 | **Required but with no mechanism in the model yet** 
207:| E10 | **A place where the model could make a testable statement** 
208:| E11 | **A place where the model could make a testable statement** 
209:| E12 | **Not addressed by the model** 
210:| E13 | **In apparent conflict** 
211:| E14 | **In apparent conflict** 
212:| E15 | **In apparent conflict** 
213:| E16 | **Required but with no mechanism in the model yet** 
214:| E17 | **Required but with no mechanism in the model yet** 
222:| E25 | **Not addressed by the model** 
223:| E26 | **Required but with no mechanism in the model yet** 
224:| E27 | **Required but with no mechanism in the model yet** 
225:| E28 | **Required but with no mechanism in the model yet** 
226:| E29 | **Not addressed by the model** 
229:| E32 | **In apparent conflict** 
230:| E33 | **Required but with no mechanism in the model yet** 
231:| E34 | **Required but with no mechanism in the model yet** 
```

```
$ grep -n 'required but with no mechanism\|not addressed by the model\|testable statement' _scratch/polarization/survey_prompt_r0.md
52:   - required but with no mechanism in the model yet;
54:   - not addressed by the model;
55:   - a place where the model could make a testable statement.
```

## Missed v3 sources (finding 3)
```
$ sed -n 52,64p research/pde_ledger_v3/steps/S11_stray_longitudinal.md
N_SO = {D=2 → 4, D=3 → 3, D=4 → 4, D=5 → 3}        N_O = 3 for every D
```

The extras are **reflection-odd**: `(tr G)·ε^{ij}G_{ij}` at `D=2`, `ε_{ijkl}G_{ij}G_{kl}` at `D=4`.
⚠ **And the walk's METHOD was the error, not its arithmetic** — *"decompose into pieces that do not mix
and count them"* misses **cross-pairings between isomorphic summands**, which is exactly what `D=2` and
`D=4` have (at `D=2` the trace and `Λ²` are both `SO(2)` scalars; at `D=4`, `Λ²` splits into self-dual ⊕
anti-self-dual). ⭐ Five independent derivations agree on the corrected counts.

⭐⭐ **The physical content of the extras, and it sharpens the step:** at `D=2` the extra invariant is
**NOT a total derivative** — its Euler–Lagrange operator is non-zero — so an `SO(2)`-invariant Lagrangian
**can mix the longitudinal and transverse sectors**. At `D=4` it **is** a total derivative and adds no
operator structure. ⇒ ⭐ **Sector separation is a `D=3` fact, ⛔ not a structural one.**
```

```
$ ls research/pde_ledger_v3/steps | grep -i '^O2'
O2_steady_brane_balance.md
```

```
$ sed -n 89,92p research/pde_ledger_v3/steps/O2_steady_brane_balance.md
exchange, loads, energy and power; optical `μ_⊥` uses the same measure as `ρ_br`. Induced metric and
native area factors remain explicit. The supplied geometric/optical/mass inputs are
`g_ij=δ_ij+∂_iξ_w∂_jξ_w`, its inverse, `ξ_w=ℓh`,
`c_γ²≡μ_⊥/ρ_br`, `c_γ=c₀(1+δ)` and `∇·(ρ_br V)=−j_n`.
```

```
$ sed -n 110,111p research/pde_ledger_v3/steps/O2_steady_brane_balance.md
loads and boundary data remain live with spatial derivatives. Radial profiles impose no constitutive
isotropy, parity, stress symmetry, absent couple/chiral content, derivative cutoff or finite history
```

```
$ sed -n 115,117p research/pde_ledger_v3/steps/O2_steady_brane_balance.md
Supplied counting is `ε=GM/(c₀²r)`, `δ=O(ε)`, `(∂ξ_w)²=O(ε)`, `V/c₀=O(ε^{1/2})`,
`(V/c₀)²=O(ε)` and the optical monomial box `0≤a≤1, 0≤b≤2, 0≤c≤1` in
`δ^a(V/c₀)^b((∂ξ_w)²)^c`, with nonnegative integer indices. Only the stiffness/density ratio
```

```
$ sed -n 53,54p research/pde_ledger_v3/steps/S11bB_interface_assembly.md
⚠ Non-reciprocal, "odd" constitutive couplings of exactly the excluded form are realised in driven
laboratory media (odd elasticity, odd viscosity, active and chiral fluids). ⛔ They are not unphysical.
```

```
$ git show ede8aa21:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md | sed -n 117,119p
- **Eikonal.** The retained object is the dispersion relation above, with its position-dependent
  coefficients, and its rays. Excluded from this step's claim and not computed: the explicit subprincipal
  terms of the underlying operator, which affect amplitude and polarization transport.
```

```
$ grep -c 'S11_stray_longitudinal.md:5[2-9]\|S11:5[2-9]\|S11:6[0-4]' _scratch/polarization/POLARIZATION_SURVEY.md
0
```

```
$ grep -c 'S11bB' _scratch/polarization/POLARIZATION_SURVEY.md
3
```

## Photon-mass bounds (finding 4)
```
$ grep -n -i 'ryutov\|bonetti\|fast radio burst\|FRB' _scratch/polarization/POLARIZATION_SURVEY.md | head
```

