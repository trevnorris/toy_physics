# Muonium Gravity as a New Falsification Track for the Fluid-Universe Toy Model

## Standalone handoff for independent review

**Status:** Research-track proposal, not a claimed prediction  
**Date:** 15 September 2026  
**Controlling correction:** All previous numerical throat-radius, throat-diameter, throat-depth, and aspect-ratio calculations for leptons are withdrawn and must be redone from first principles.

---

## 0. Purpose of this handoff

A newly demonstrated cold, directional muonium beam makes a direct measurement of the gravitational response of a second-generation lepton experimentally plausible. This matters for the Fluid-Universe toy model because the existing research stack has generally matched a universal point-particle gravitational response at the reduced far-field level, while the detailed map from an actual electron or muon throat solution to its inertial and gravitational coefficients has not yet been derived from the full moving-throat PDE.

The resulting opportunity is real but must be stated carefully:

1. The muon's rest energy and inertial properties are known extremely well from nongravitational experiments.
2. A direct measurement of the passive gravitational response of a muon has not yet been completed.
3. The proposed experiment uses neutral muonium, \(\mathrm{Mu}=\mu^+e^-\), because a bare charged muon would be overwhelmed by electromagnetic backgrounds.
4. The experiment will principally test how muonium responds to Earth's gravitational field. It will not directly measure the gravitational field sourced by a muon.
5. The toy model must not use assumed equivalence to infer throat geometry and then claim that geometry predicts equivalence. Electron and muon branches must instead be solved from nongravitational constraints, after which their passive gravitational responses must be calculated as outputs.

This document therefore adds a new branch to the research program: derive the electron and muon defect solutions without using gravitational universality as an input, compute their inertial and passive gravitational response coefficients separately, and compare the muonium prediction with the forthcoming measurement.

---

## 1. Experimental development

Researchers at ETH Zurich and the Paul Scherrer Institute have produced a high-intensity, superthermal muonium beam using superfluid helium. The beam is cold and directional enough to support a planned atom-interferometry measurement of muonium's response to Earth's gravity. The collaboration expects initial beam tests with the interferometer first, followed—if development proceeds successfully—by the gravity experiment in roughly two to three years.

Muonium contains a positive antimuon and an electron:

\[
\mathrm{Mu}=\mu^+e^-.
\]

It is electrically neutral and contains no proton or neutron, so nuclear finite-size and strong-interaction structure are absent at leading order. Its mass is dominated by the antimuon:

\[
m_{\mathrm{Mu}}\simeq m_\mu+m_e-E_{\mathrm{bind}}/c^2,
\qquad
\frac{m_\mu}{m_e}\simeq 206.77.
\]

The electron therefore contributes only about \(0.48\%\) of the constituent rest mass, apart from the tiny binding correction. A muonium free-fall measurement is consequently a particularly clean probe of the passive gravitational behavior of a second-generation charged lepton.

The experimental sources should be treated as follows:

- [ETH Zurich's announcement](https://ethz.ch/en/news-and-events/eth-news/news/2026/09/novel-particle-beam-could-challenge-einsteins-theory-of-gravity.html) states that this would be the first direct investigation of how a second-generation particle falls and describes the planned interferometer and two-to-three-year horizon.
- [The LEMING experiment page](https://lepp.ethz.ch/research/leming.html) describes the muonium composition, the approximately \(200{:}1\) muon-to-electron mass dominance, the lack of a completed direct gravity measurement, and the atom-interferometry concept.
- The beam result is reported as J. Zhang et al., *Generation of a high-intensity, superthermal muonium beam for gravity and laser spectroscopy experiments*, *Nature Physics* (2026), DOI: [10.1038/s41567-026-03433-x](https://doi.org/10.1038/s41567-026-03433-x).

The beam paper is an enabling result, not yet the gravitational measurement.

---

## 2. What is known and what is not

The phrase “the gravitational mass of the muon is unknown” is too loose. Several different quantities must be separated.

### 2.1 Inertial/rest energy

The muon's rest mass \(m_\mu\) is precisely known from particle physics. Its energy-momentum relation, response to electromagnetic forces, decay, magnetic moment, and its role in muonium spectroscopy provide strong nongravitational constraints.

In the toy model, the inertial mass of a stationary defect should ultimately be obtained from the fully renormalized defect energy and the translational collective-coordinate action. Schematically,

\[
E_\mu^{\rm rest}=m_{i,\mu}c^2,
\]

but the equality here defines the inertial coefficient of the effective branch; it does not by itself determine the throat radius.

### 2.2 Passive gravitational response

Passive gravitational response controls the force or acceleration of the particle in an externally produced weak gravitational field. Introduce

\[
L_{\rm test}^{(s)}
=
\frac12m_{i,s}v^2-m_{p,s}\Phi_N+\cdots,
\]

where \(s\) labels the defect species. Then

\[
\mathbf a_s=-\frac{m_{p,s}}{m_{i,s}}\nabla\Phi_N+\cdots.
\]

Define the species-dependent weak-equivalence coefficient

\[
\boxed{
\kappa_s^{\rm pass}
\equiv
\frac{m_{p,s}}{m_{i,s}}
}
\]

so that, in a sufficiently uniform terrestrial field,

\[
\frac{a_s}{g}=\kappa_s^{\rm pass}
\]

after conventional and composition-dependent corrections are included.

General relativity predicts \(\kappa_s^{\rm pass}=1\). The muonium experiment is designed to test this directly for a system dominated by a second-generation antilepton.

### 2.3 Active gravitational source strength

Active gravitational mass controls the field produced by the particle. Schematically,

\[
\nabla^2\Phi_s=4\pi G\,m_{a,s}\delta^3(\mathbf x)+\cdots.
\]

Define

\[
\kappa_s^{\rm act}\equiv\frac{m_{a,s}}{m_{i,s}}.
\]

The planned muonium interferometer does not directly measure \(\kappa_\mu^{\rm act}\). The project must not confuse a result about free fall with a result about how strongly the muon sources the far gravitational field.

### 2.4 Self-energy and geometry

The defect geometry—mouth radius, depth, wall profile, transverse localization, support-mode content, and radial throughput—is a fourth category. Denote a stationary branch by

\[
\mathcal B_s=
\{R_s(\Omega,w),\psi_s,A_{A,s},\chi_{B,s},\ldots\}.
\]

Then the relevant outputs are functionals of the complete solution:

\[
m_{i,s}=\mathcal M_i[\mathcal B_s],
\qquad
m_{p,s}=\mathcal M_p[\mathcal B_s],
\qquad
m_{a,s}=\mathcal M_a[\mathcal B_s].
\]

There is no justified one-variable rule of the form “larger mass means a throat of radius \(a=f(m)\)” until the actual branch equations establish it.

---

## 3. Mandatory reset of the throat-size work

The controlling instruction for all further work is:

> Every previous calculation assigning an absolute or relative throat radius, diameter, depth, volume, or aspect ratio to the electron, muon, or tau is invalid as a physical result and must not be used as an input.

This includes, at minimum:

- the preferred value \(L/a\approx1.85\) when treated as the actual lepton geometry rather than a provisional reduced branch;
- any numerical electron radius, muon radius, or tau radius inferred from rest mass;
- any claim that a heavier lepton must be narrower, wider, deeper, or shorter;
- the re-equilibrated law \(a_j\propto(2j+1)^{-1}\) when applied to real lepton families;
- reverse-engineered heavy-lepton “deep needle” geometries;
- geometric source dictionaries such as \(m_G\sim\kappa_m\rho_0\pi a^2L\) when used quantitatively before \(\kappa_m\), the profile, and the full branch have been derived;
- any gravitational conclusion whose only support is one of those geometric estimates.

### 3.1 What can survive provisionally

Withdrawing the throat-size calculations does not require discarding every symbolic result in the existing lepton program. Results may survive in one of three statuses.

**Geometry-independent identities.** Exact field identities, projection formulas, conservation laws, and algebraic statements that do not use the invalid geometry may remain intact.

**Conditional reduced theorems.** For example, stationarity of the assumed reduced functional

\[
F(a,\rho)=\frac{A(\rho)}a+\frac{B(\rho)}{a^2}+C(\rho)a^3
\]

implies

\[
E_w+2E_f=3E_{\rm PV}.
\]

That algebra remains correct for that functional. What is withdrawn is the claim that the functional and its coefficients are already the correct microscopic electron or muon branch.

**Branch-dependent results requiring revalidation.** The \(11{:}2{:}5\) energy partition, the breathing slope \(-57/64\), the absolute mass formula, the D/N family-radius scaling, and every numerical geometry conclusion depend on frozen response choices and/or the provisional geometry. They must be labeled conditional and rerun after the new branch is derived.

### 3.2 What this reset does not automatically invalidate

The separate far-field gravity program matched the Newtonian and post-Newtonian target through a declared reduced closure. Invalidating lepton throat-size calculations does not algebraically erase that far-field match. It changes the interpretation:

- the existing 0PN–4PN work remains evidence that a chosen universal effective point-particle branch can reproduce the desired gravitational dynamics within its closure;
- it does not prove that every microscopic defect species reduces to that same branch;
- the missing task is the defect-to-worldline matching calculation for each particle family.

The distinction prevents both overreaction and circular reasoning.

---

## 4. Correction to the status of \(\kappa_\rho=1\)

The current 1PN documents use the reduced scalar worldline ansatz

\[
L_{\rm sc}
=
-M_0\left(1+\kappa_\rho\frac{\Phi_N}{c^2}\right)c^2
\sqrt{1-
\frac{v^2/c^2}{1+(n-1)\Phi_N/c^2}}
.
\]

Expanding at Newtonian order gives

\[
L_{\rm sc}
=
-M_0c^2+\frac12M_0v^2-
\kappa_\rho M_0\Phi_N+\cdots.
\]

Comparison with the standard target

\[
L_N=\frac12M_0v^2-M_0\Phi_N
\]

fixes

\[
\kappa_\rho=1.
\]

This deserves a precise status statement:

- \(\kappa_\rho=1\) is not derived solely from projected continuity.
- It is not an arbitrary hidden number either; it is fixed by matching the reduced point-particle action to Newtonian gravity.
- That matching proves that the selected effective branch has ordinary Newtonian response by construction.
- It does not yet prove that a first-principles muon defect branch must produce the same coefficient.

For the new track, the species-specific coefficient should initially remain open:

\[
\kappa_{\rho,s}
=
\left.
\frac{1}{m_{i,s}c^2}
\frac{\partial E_s^{\rm eq}}{\partial(\Phi_N/c^2)}
\right|_{\Phi_N=0},
\]

with sign conventions checked against the effective action. Only after computing this derivative from the relaxed defect solution should the branch be compared with \(\kappa_{\rho,s}=1\).

The earlier global far-field work may retain \(\kappa_\rho=1\) as its GR-matching closure. The new lepton calculation must test whether the electron and muon branches independently reproduce it.

---

## 5. The clean research question

The central question is:

> Can the common brane–bulk parent theory produce electron and muon stationary defect branches that reproduce their known nongravitational observables while predicting their passive gravitational response without assuming universality?

More formally, solve for each species \(s\in\{e,\mu\}\) a branch \(\mathcal B_s\) constrained by nongravitational data \(\mathcal D_s^{\rm nongrav}\):

\[
\mathcal E[\mathcal B_s;\theta]=0,
\qquad
\mathcal O_k[\mathcal B_s;\theta]
=D_{k,s}^{\rm nongrav},
\]

where \(\theta\) contains only common parent-medium parameters and legitimate discrete branch labels. Then calculate

\[
m_{i,s},\qquad
m_{p,s},\qquad
m_{a,s},\qquad
\kappa_s^{\rm pass}=m_{p,s}/m_{i,s}
\]

without fitting to the muonium gravity result.

The experimentally relevant prediction is not a throat diameter by itself. It is

\[
\boxed{
\frac{a_{\rm Mu}}g-1
}
\]

or the equivalent interferometric phase shift after the apparatus transfer function and muonium composite corrections are included.

---

## 6. Nongravitational constraints for the branch solve

The throat geometry must be inferred through a simultaneous inverse problem, not from rest energy alone. At minimum, the candidate electron and muon branches should be confronted with the following classes of constraints.

### 6.1 Rest energy and translational inertia

The on-shell defect energy must reproduce the measured rest energy, while the coefficient of the translational collective-coordinate kinetic term must reproduce the same inertial mass to the required precision:

\[
E_s[\mathcal B_s]=m_s^{\rm obs}c^2,
\qquad
L_{\rm trans}^{(s)}=\frac12m_{i,s}\dot{\mathbf X}_s^2+\cdots.
\]

Equality of these two quantities should be derived from the model's symmetry and stress-energy structure, not imposed twice as two independent fits.

### 6.2 Electric charge

The electron and muon have the same charge magnitude. In the corrected ontology, charge sign and magnitude are not to be silently identified with circulation or radial throughput. The branch must reproduce

\[
|q_e|=|q_\mu|=e
\]

through the puncture orientation, microscopic coupling, and brane localization map while allowing different internal excitations.

This is a strong family constraint: a model that obtains the muon mass merely by changing a parameter that also changes observable charge has failed.

### 6.3 Magnetic moment and spin sector

The branch must eventually reproduce the appropriate spin label and magnetic moment, including the already-known limitations of the current reduced same-charge Berry/rotor construction. The existing work only provides a conditional route; it does not yet supply a completed free-particle spin theorem.

### 6.4 Muonium spectroscopy

Muonium energy levels and hyperfine structure strongly constrain any model-dependent change to electron–muon electromagnetic interaction, size, polarizability, or internal response. A large throat cannot be hidden if it introduces excluded finite-size or frequency-dependent effects.

### 6.5 Lifetime and decay channel

The muon lifetime must be treated as either an output of the open-system throat dynamics or a constraint on it. A candidate heavy-lepton branch that is stable in the model without a physically justified decay outlet is incomplete even if it matches the rest energy.

### 6.6 Scattering and effective size

The geometric mouth radius need not equal a conventional scattering radius, but the model must derive the map between them. Existing experimental bounds constrain form factors and contact-like behavior. Geometry cannot be declared observationally invisible without calculating its coupling to probes.

---

## 7. Deriving inertial mass from the defect branch

The most reliable definition of inertial mass is dynamical. Begin with a stationary solution \(\mathcal B_s^0\) and introduce a slowly varying collective position \(\mathbf X_s(t)\). Substitute

\[
\mathcal B_s(\mathbf x,w,t)
=
\mathcal B_s^0(\mathbf x-\mathbf X_s(t),w)
+\delta\mathcal B_s
\]

into the parent action, solve the constraint fields and induced wake to the required order in velocity, and integrate over the bulk and brane coordinates. The effective action should take the form

\[
S_{\rm eff}^{(s)}
=
\int dt\left[
-E_s^0
+\frac12M_{AB}^{(s)}\dot X_s^A\dot X_s^B
+O(v^4)
\right].
\]

For an isotropic brane particle,

\[
M_{AB}^{(s)}=m_{i,s}\delta_{AB}.
\]

This calculation automatically includes, if treated correctly:

- energy stored in the stationary support field;
- entrained or added inertia of the surrounding medium;
- the velocity-induced wake;
- wall deformation and retuning of trapped modes;
- gauge-field inertia;
- counterterms or background subtraction needed to isolate finite defect energy.

The old shortcut that assigned a fixed added-mass coefficient from an assumed geometry may remain a reduced benchmark, but it cannot close the species-specific inertial calculation.

---

## 8. Deriving passive gravitational response

Apply a weak, slowly varying external gravitational background generated by a distant ordinary source. In the model this must be represented by the appropriate perturbation of the common medium—not merely inserted as a Newtonian potential in the final worldline action.

Let \(\epsilon\) control the external field, with

\[
\Phi_N(\mathbf x)=\epsilon\,\widehat\Phi(\mathbf x).
\]

For each species, solve the constrained relaxation problem

\[
\mathcal B_s(\epsilon)
=
\arg\operatorname*{ext}_{\mathcal B}
E[\mathcal B;\epsilon]
\]

subject to fixed topological labels and the correct adiabatic or dynamical protocol. The on-shell energy is

\[
E_s^{\rm eq}(\epsilon)
=
E_s^0+m_{p,s}\Phi_N(\mathbf X_s)+O(\Phi_N^2,\nabla^2\Phi_N).
\]

The envelope theorem can simplify the first derivative because the background branch is stationary, but all explicit couplings of the external medium perturbation to matter, gauge, wall, support, and leakage sectors must be retained.

The resulting passive coefficient is

\[
\boxed{
m_{p,s}
=
\left.
\frac{\partial E_s^{\rm eq}}{\partial\Phi_N}
\right|_{\Phi_N=0}
}
\]

up to the chosen sign convention for \(\Phi_N\).

The important conceptual point is that internal geometry can change \(m_{p,s}\) and \(m_{i,s}\) differently. A different muon throat may therefore permit \(\kappa_\mu^{\rm pass}\neq1\), but only if the actual perturbation theory produces that result. Geometry alone is not evidence of nonuniversality.

---

## 9. Deriving active source strength separately

For completeness, the active gravitational coefficient should be computed from the far-field monopole sourced by the isolated branch. If the projected far field has

\[
\Phi_s(r)
=
-\frac{Gm_{a,s}}r+O(r^{-2}),
\]

then \(m_{a,s}\) is obtained from the coefficient of the \(1/r\) tail after all leakage and projection contributions have been included.

This calculation is logically independent of the weak-background response derivative. The model may ultimately enforce

\[
m_{a,s}=m_{p,s}=m_{i,s}
\]

through a common Ward identity, translational symmetry, energy conservation, or a virial theorem. If so, the equality should emerge as a theorem of the parent action and boundary conditions. Until then, the three coefficients remain distinct entries in the ledger.

---

## 10. From a muon prediction to a muonium prediction

The experiment measures a bound composite, not an isolated muon. The model therefore needs a composite calculation.

At leading order, one might write

\[
m_{i,\rm Mu}
=m_{i,\mu}+m_{i,e}+E_{\rm bind}/c^2,
\]

\[
m_{p,\rm Mu}
=m_{p,\mu}+m_{p,e}+m_{p,\rm bind},
\]

so that

\[
\kappa_{\rm Mu}^{\rm pass}
=
\frac{m_{p,\mu}+m_{p,e}+m_{p,\rm bind}}
{m_{i,\mu}+m_{i,e}+E_{\rm bind}/c^2}.
\]

Because \(m_\mu/m_e\simeq206.77\), the muon dominates. Nevertheless, a serious prediction must include:

- the electron's own passive coefficient;
- electromagnetic binding energy;
- any model-specific interaction energy between the two throats;
- polarization or tidal response in the external field;
- the fact that the constituent is an antimuon \(\mu^+\), which makes any matter–antimatter asymmetry relevant;
- apparatus-level electromagnetic and inertial systematics only when translating the theory result into the measured phase.

If the electron branch is already constrained to ordinary free fall, then approximately

\[
\kappa_{\rm Mu}^{\rm pass}-1
\simeq
\frac{m_\mu}{m_\mu+m_e}
\left(\kappa_\mu^{\rm pass}-1\right)
+\delta_{\rm bind},
\]

where \(\delta_{\rm bind}\) contains binding and composite corrections. The prefactor is about \(0.995\), so the muonium signal nearly preserves a muon anomaly.

---

## 11. Relationship to the current research track

This should be integrated as a bounded branch of the existing far-field-first program, not allowed to derail the main program.

### 11.1 Work that remains the main line

The existing track continues to prioritize:

1. completing and auditing the common brane–bulk equations;
2. deriving the remaining moving-throat PDE closures;
3. checking that the light, electromagnetic, and gravitational sectors arise from the same medium without contradictory parameter demands;
4. preserving the frozen far-field results where their derivations do not depend on invalid lepton geometry;
5. postponing strong-field particle claims until the underlying branch exists.

### 11.2 New parallel branch

The muonium branch should begin with definitions and sensitivity analysis now, but numerical prediction should wait until the stationary and moving defect machinery is sufficient. Its near-term deliverables are:

- a clean three-mass ledger \((m_i,m_p,m_a)\);
- a list of exactly where universality entered the existing 0PN–4PN reductions;
- a dependency audit of every lepton result that used \(L/a\), \(a\), \(L\), or \(\kappa_\rho=1\);
- a source-free definition of the electron and muon branch constraints;
- a perturbative weak-background calculation design;
- a blinded prediction protocol established before the muonium result is known.

### 11.3 Why this is not permission to tune

The experimental uncertainty does not create a free parameter that may be selected later. A scientifically useful result requires:

- common parent-medium constants for electron and muon;
- no muon-specific continuous knob introduced solely to alter gravity;
- branch labels justified by topology, boundary conditions, or a solved spectrum;
- all parameters fixed by nongravitational data before computing gravity;
- publication or timestamping of the predicted interval before the experimental result.

A deviation chosen after seeing the measurement would not count as a prediction.

---

## 12. Dependency audit of current results

The next session should build a machine-readable dependency graph. The initial classification is:

| Existing result | Current status after geometry reset | Required action |
|---|---|---|
| Exact 4D parent continuity and field identities | Retain unless separately contradicted | Verify no hidden geometry substitution |
| 4D-to-3D projection formulas | Retain conditionally on projection map | Audit species dependence and leakage |
| Newtonian point-particle match with \(\kappa_\rho=1\) | Retain as target-matched effective branch | Re-derive \(\kappa_{\rho,e}\) and \(\kappa_{\rho,\mu}\) microscopically |
| Conservative EIH/PN coefficient matching | Retain inside declared universal closure | Determine whether microscopic branches land on that closure |
| \(F=A/a+B/a^2+Ca^3\) | Conditional ansatz | Derive or replace from the full branch |
| Virial identity \(E_w+2E_f=3E_{\rm PV}\) | Algebraically valid for that ansatz | Do not call it a lepton theorem until ansatz is derived |
| \(11{:}2{:}5\) energy partition | Reopened | Recompute without importing species universality |
| \(d\ln a/d\ln\rho=-57/64\) | Reopened | Derive from the actual response Hessian |
| \(L/a\approx1.85\) | Withdrawn as physical geometry | Recompute from boundary-value problem |
| D/N spectrum on an assumed finite interval | Retain as exact spectrum of that boundary problem | Revalidate boundary conditions and geometry |
| \(a_j\propto(2j+1)^{-1}\) | Withdrawn as a lepton-radius prediction | Re-derive only after valid energy functional exists |
| Support-only mass ladder \(1{:}9{:}25\) fails | Retain as a no-go for that specific model | Do not generalize beyond the tested ansatz |
| Reverse-engineered deep heavy-lepton geometries | Withdrawn | Replace with forward branch solutions |
| Same-charge Berry/spin corridor | Conditional and incomplete | Keep separate from gravity until physical closure is derived |
| Active-source dictionary \(m_G\sim\rho_0\pi a^2L\) | Qualitative only | Derive normalization and profile dependence |

---

## 13. Proposed calculation sequence

### Phase A — Clean the ledger

1. Search every current paper, note, and CAS script for \(a\), \(L\), \(L/a\), diameter, radius, volume, \(m_G\), \(\kappa_\rho\), \(\kappa_{\rm PV}\), and the \(11{:}2{:}5\) partition.
2. Label each occurrence as exact, ansatz-dependent, fitted, target-matched, reverse-engineered, or obsolete.
3. Add explicit warnings to downstream documents so invalid geometry is not silently reused.
4. Create a “do not import” list for the fresh particle solve.

### Phase B — Define observables at the parent-action level

1. Define finite background-subtracted defect energy.
2. Define translational inertial mass from the moduli-space kinetic tensor.
3. Define passive mass from the on-shell weak-background energy derivative.
4. Define active mass from the projected far-field monopole.
5. Prove gauge and projection invariance of each definition.
6. Establish sign conventions with a simple ordinary-body benchmark.

### Phase C — Solve the electron branch first

1. Solve the stationary throat boundary-value problem without assigning an external radius.
2. Fix common medium parameters using the nonparticle sectors wherever possible.
3. Apply electron nongravitational constraints.
4. Test uniqueness and stability of the resulting branch.
5. Calculate \(m_{i,e}\), \(m_{p,e}\), and \(m_{a,e}\).
6. Check whether the branch independently gives the effective universal closure already used in the far-field papers.

### Phase D — Solve the muon branch

1. Use the same parent-medium constants.
2. Allow only justified discrete or dynamical branch differences.
3. Fit no gravitational datum.
4. Reproduce the muon rest energy, equal charge magnitude, magnetic constraints, scattering constraints, and decay behavior to the level supported by the model.
5. Determine the actual geometry as an output.
6. Calculate \(\kappa_\mu^{\rm pass}\) and \(\kappa_\mu^{\rm act}\).

### Phase E — Build muonium

1. Solve or consistently reduce the electron–antimuon bound state.
2. Include binding energy and interaction stress in both inertial and passive ledgers.
3. Derive \(a_{\rm Mu}/g\).
4. Translate the result into the interferometer's predicted phase shift only after the experimental transfer function is public.

### Phase F — Freeze the prediction

1. Produce a central value and uncertainty interval.
2. Separate numerical uncertainty, truncation error, parameter uncertainty, and branch ambiguity.
3. Timestamp the calculation and all scripts.
4. State in advance what experimental outcomes falsify the branch, the lepton construction, or the broader parent model.

---

## 14. Required symbolic and numerical checks

Each phase should produce a reproducible Mathematica and/or SymPy artifact with printed verification output. Minimum checks include:

### 14.1 Variational checks

- independently vary the stationary fields and the geometry field;
- verify boundary terms and junction conditions;
- confirm the stationary solution makes first-order internal variations vanish;
- verify that the passive-mass derivative includes all explicit external-field dependence.

### 14.2 Zero-mode and normalization checks

- identify translational zero modes;
- normalize collective coordinates consistently;
- distinguish gauge zero modes from physical modes;
- verify finite background subtraction;
- test coordinate and projection-map independence.

### 14.3 Consistency identities

- energy conservation for the closed conservative problem;
- correct accounting when leakage or decay ports are opened;
- equality or difference between rest energy and kinetic inertia;
- far-field Gauss-law normalization;
- recovery of the existing ordinary point-particle limit when the species-dependent structure is suppressed.

### 14.4 Numerical robustness

- convergence under grid refinement in \(r\) and \(w\);
- domain-size independence of renormalized quantities;
- continuation in external-field strength to verify the linear regime;
- branch tracking and detection of bifurcations;
- Hessian spectrum for stability;
- sensitivity to wall thickness and localization profiles;
- uncertainty propagation into \(\kappa_\mu^{\rm pass}\).

---

## 15. Possible outcomes and their meaning

### Outcome 1: Exact or effectively exact universality

\[
\kappa_e^{\rm pass}=\kappa_\mu^{\rm pass}=1
\]

emerges from the actual branches. This would strengthen the toy model because equivalence would no longer be merely inherited from target matching. The muonium experiment would become a consistency test rather than a discriminator from general relativity.

### Outcome 2: Small, parameter-free muon deviation

\[
\kappa_e^{\rm pass}\simeq1,
\qquad
\kappa_\mu^{\rm pass}=1+\Delta_\mu,
\qquad
\Delta_\mu\neq0.
\]

If \(\Delta_\mu\) is obtained after all nongravitational constraints are frozen, this would be a genuine novel prediction. It would also imply that the existing universal point-particle closure is not the exact microscopic reduction for every species and would require a controlled species-dependent extension of the PN ledger.

### Outcome 3: Large deviation already excluded indirectly

The branch may mathematically predict nonuniversality but conflict with muonium spectroscopy, particle kinematics, astrophysical energy-loss bounds, antimatter-gravity results, or other precision tests. Such a branch is falsified before LEMING runs.

### Outcome 4: No acceptable muon branch

If the parent equations cannot reproduce the known muon properties with the same common-medium constants used for the electron, the lepton-family construction fails. This is valuable falsification and should not be repaired by adding an unconstrained species knob.

### Outcome 5: Gravity cannot be calculated without importing \(\kappa=1\)

Then the current model has no muonium-gravity prediction. The honest conclusion would be that it reproduces the GR-like far field under a declared closure but has not derived equivalence from its particle ontology.

---

## 16. Falsification rules to freeze now

Before detailed solving begins, the following rules should be adopted.

1. No prior throat dimensions may be used as initial physical data. They may only be used as numerical seed guesses, and final solutions must be demonstrably seed-independent.
2. No gravitational measurement may be used to fit the electron or muon branch before the prediction is frozen.
3. The electron and muon must share all genuinely universal parent-medium parameters.
4. Any species-specific continuous parameter must be derived from a dynamical invariant or measured nongravitationally; otherwise it is a forbidden tuning knob.
5. Passive and active gravitational coefficients must be reported separately.
6. The muonium composite correction must be calculated rather than equated blindly to the bare-muon result.
7. If multiple stable branches survive all nongravitational constraints, the project must report a prediction set or interval rather than select the branch closest to the eventual measurement.
8. A post-experiment branch choice does not count as a successful prediction.

---

## 17. Questions for the reviewing session

The next review should answer these questions in order.

1. Where exactly does the parent action couple an external weak gravitational background to the throat, support field, gauge field, and wall geometry?
2. Is there a Noether or Ward identity that forces \(m_p=m_i\) for every localized defect, regardless of internal structure?
3. If such an identity exists, which assumptions—closed system, Lorentz symmetry, absence of leakage, common characteristic speed, boundary conditions—are required?
4. Can the electron and muon be distinct stationary branches of one parameter set without reverse-engineering their geometry from mass?
5. Which known muon observables most strongly constrain the allowed branch before gravity is considered?
6. Does the antimuon orientation introduce any passive-gravity sign or magnitude difference in this ontology, and is that compatible with existing antimatter constraints?
7. Which existing PN results depend only on the universal far-field sector, and which depend on the now-withdrawn particle geometry?
8. Can \(m_i\), \(m_p\), and \(m_a\) be extracted from one numerical solution with independent diagnostics?
9. What experimental precision and observable will LEMING actually publish, and how should the model prediction be mapped to it?
10. What result would falsify only the proposed muon branch, and what result would falsify the common parent model?

---

## 18. Recommended immediate deliverable

The next technical document should not attempt a new numerical throat diameter. It should be titled something like:

> **Defect Inertia and Passive Gravity: Observable Definitions and a Species-Dependence Audit**

Its job should be limited to:

1. deriving \(m_i\), \(m_p\), and \(m_a\) from the current parent action;
2. locating every current use of universal response;
3. proving which equalities follow from symmetry and which were target-matched;
4. stating the minimum boundary-value problem required for a forward electron and muon solve;
5. producing no particle-size number until that boundary-value problem is solved.

Only after that foundation is secure should the program calculate new electron or muon throat geometry.

---

## 19. Bottom-line assessment

The new muonium experiment does not show that the muon probably falls differently, nor does it by itself validate a variable-throat interpretation. It reveals that a foundational relationship—well motivated theoretically and supported broadly, but not yet directly tested for a second-generation lepton—will become experimentally accessible.

For this toy model, that creates a valuable falsification window precisely because the microscopic lepton geometry and the defect-to-worldline equivalence map remain unfinished. The correct response is not to adjust an old throat-radius estimate. The correct response is to discard those estimates, derive the electron and muon branches from nongravitational physics, calculate inertial and gravitational coefficients independently, and let the experiment decide.

The research opportunity can be summarized in one line:

\[
\boxed{
\text{nongravitationally fixed defect branch}
\;\Longrightarrow\;
\text{predicted }\kappa_{\rm Mu}^{\rm pass}
\;\Longrightarrow\;
\text{direct experimental test}
}
\]

That is a clean, bounded, and potentially decisive addition to the current research development track.
