# Polarization: experiments and the brane model

Survey date: 2026-10-09. Deliverable: explain polarization, assemble a cited experimental record, inventory the repository's claims, and classify the experimental facts against the current model.

The model under discussion is the user's finite-thickness ordered slab in four spatial bulk dimensions, centred at w = 0. Repository conclusions below use the v3 ledger as authority, as requested. “Two modes” always needs its field content and assumptions attached.

## Open items from the final review (read before using Part 4)

The body below is the round-3 text as reviewed (sha256 `f1baeff2…`). Its final review found two problems. By the
user's decision they are listed here rather than repaired. **Where an item applies, it overrides the body.** Record:
`_measurements/polarization_survey/survey_r3_review_disposition.md`, with its lookups beside it.

1. **E01 is not "Reproduced".** The ledger itself says the two-direction count restates a supplied input:
   "`D_brane = 3` went in and `D−1 = 2` came out. Without a delivered `D_brane` the sentence is an assumption
   restated, ⛔ not a result" (`research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:134–136`, R-S1-01, OPEN). S10
   computes the conditional map `D ↦ D − 1`, "conditional on the supplied action, the supplied `[u]`, and BOTH
   structural premises" (`research/pde_ledger_v3/V3_STEP_PLAN.md:374–376`). It also rests on the separation of `u`
   from out-of-plane fields being "inherited rather than tested" (`steps/S10_two_transverse_photons.md:182–183`).
   Under rule 4's precedence, E01 is **Required, no mechanism yet**, with those two supplied inputs named. The
   totals become: reproduced 0; in apparent conflict 7; testable 10; required, no mechanism yet 6; not addressed 15.
   The summary's sentence that E01 "reproduces the conditional in-plane transverse count" should read as "S10
   computes the conditional count from supplied inputs that state it".
2. **The parity/chirality condition is not located at `S11bB_interface_assembly.md:53–63`.** Those lines concern
   non-reciprocal, non-passive interface couplings (line 53 opens: Non-reciprocal, "odd" constitutive couplings),
   not parity, so splitting the two helicities does not, by that source, need a named reservoir. The interface specs
   *supply* in-plane parity: "In-plane isotropy _and_ parity — the full `O(3)` acting on the three in-plane
   directions" (`directives/S11b_SHARED_PHYSICS.md:284`; `directives/S11bB_SHARED_PHYSICS.md:358`), with
   "Time-reversal is NOT assumed" (`S11b_SHARED_PHYSICS.md:288`). So S11bB:76's "in-plane parity admits no
   `e_W ↔ u_T` bilinear" follows from that supplied symmetry. Read E10–E15, and Part 3's parity rows, as
   locating the open condition in what the ledger leaves open (`steps/O2_steady_brane_balance.md:110–111`), what it
   supplies (the specs' in-plane parity) and what it excludes (the S9b spec at `ede8aa21`, lines 117–119). No class
   changes, and the condition stays unresolved.

## 1. What polarization is

### Plain language

Imagine a wave travelling towards you. Its electric field can wiggle horizontally, vertically, or in a combination of those directions. **Polarization describes that wiggle**, rather than the wave's travel direction. A monochromatic wave's field tip traces a line, circle, or ellipse at a fixed position. These are classical electromagnetic-wave descriptions. The electric and magnetic fields are linked; they are not two independently selectable polarizations. [Maxwell (1865)](https://doi.org/10.1098/rstl.1865.0008).

- **Linear:** the tip moves back and forth along a fixed line.
- **Circular:** two perpendicular wiggles have equal amplitudes and a quarter-cycle phase difference; the tip rotates at constant magnitude.
- **Elliptical:** the general coherent combination; the tip traces an ellipse. Linear and circular are limiting cases. A phase difference means one component's cycle is offset from the other's. [Jones (1941)](https://doi.org/10.1364/JOSA.31.000488).
- **Unpolarized:** over the measurement's averaging interval, there is no preferred polarization. It is a statistical description, not a third vibration direction.
- **Partially polarized:** some preferred polarization survives averaging. Randomness, mixtures of independent sources, or averaging over different positions/frequencies can reduce the measured polarization. [Stokes (1852)](https://www.cambridge.org/core/books/abs/mathematical-and-physical-papers/on-the-composition-and-resolution-of-streams-of-polarized-light-from-different-sources/A021426E1F3E9BD366E8E1F967FBE3CE), [Wolf (1959)](https://doi.org/10.1007/BF02725127).

For a programmer, the coherent state resembles a two-component complex array. Horizontal and vertical are a basis; circular states are another basis of the same space. There are infinitely many possible polarization settings, but **two independent basis states** for a fixed free-photon momentum. Position, frequency and beam shape add other degrees of freedom; they do not increase that polarization dimension. [Jones (1941)](https://doi.org/10.1364/JOSA.31.000488), [Wigner (1939)](https://www.math.utoronto.ca/mgualt/courses/25-QM/docs/Wigner-1939.pdf).

### Precise description

For a source-free vacuum plane wave travelling along z,

\[
\mathbf E(z,t)=\mathrm{Re}\{(E_x\hat{\mathbf x}+E_y\hat{\mathbf y})e^{i(kz-\omega t)}\}.
\]

**Classical reason (Maxwell's theory):** Gauss's law, \(\nabla\cdot\mathbf E=\rho/\epsilon_0\), fixes the longitudinal electric field from charges; it is a constraint, not an independent propagating wave. For a source-free plane wave with nonzero \(\mathbf k\), \(\mathbf k\cdot\mathbf E=\mathbf k\cdot\mathbf B=0\) and \(\mathbf B=\mathbf k\times\mathbf E/\omega\). The plane perpendicular to \(\mathbf k\) in three spatial dimensions has two independent directions. This is a **plane-wave radiation statement**: near fields, media and tightly focused beams can have a component along a chosen optical axis without adding a third free-photon polarization. [Maxwell (1865)](https://doi.org/10.1098/rstl.1865.0008), [Dorn, Quabis & Leuchs (2003)](https://arxiv.org/abs/physics/0310007).

A **Jones vector** is \(J=(E_x,E_y)^T\). Examples, ignoring intensity and common phase, are \((1,0)^T\), \((0,1)^T\), \((1,1)^T/\sqrt2\), and \((1,\pm i)^T/\sqrt2\). Jones matrices describe coherent linear optical transformations. One Jones vector cannot represent a general incoherent mixture. [Jones (1941)](https://doi.org/10.1364/JOSA.31.000488).

**Stokes parameters** include mixtures. With an averaging operation appropriate to the detector, one convention is

\[
\begin{aligned}
S_0&=\langle |E_x|^2+|E_y|^2\rangle,&
S_1&=\langle |E_x|^2-|E_y|^2\rangle,\\
S_2&=2\,\mathrm{Re}\langle E_x^*E_y\rangle,&
S_3&=2\,\mathrm{Im}\langle E_x^*E_y\rangle .
\end{aligned}
\]

These are obtained from intensities measured in horizontal/vertical, diagonal/antidiagonal and opposite circular bases. The degree of polarization is
\(P=\sqrt{S_1^2+S_2^2+S_3^2}/S_0\).
On the **Poincaré sphere**, the normalized triple \((S_1,S_2,S_3)/S_0\) lies on the surface for fully polarized light, inside for partially polarized light, and at the centre for unpolarized light. Linear states lie on the equator, circular states at the poles, and other surface points are elliptical. Orthogonal pure states are antipodes. Signs and axis ordering vary between conventions. [Wolf (1959)](https://doi.org/10.1007/BF02725127).

**Spin angular momentum** is associated with polarization. An ideal circular plane/paraxial wave has angular momentum along its travel direction divided by energy equal to \(\pm1/\omega\). Quantum mechanically the corresponding photon carries spin projection \(\pm\hbar\); its **helicity** is that projection along its own momentum, conventionally labelled \(\pm1\). A linear photon state is a superposition of the two helicities, not a third helicity-zero photon. An elliptical pure state has unequal helicity probabilities. Photon “spin 1” therefore does not mean three freely propagating vacuum polarization states. [Beth (1936)](https://doi.org/10.1103/PhysRev.50.115), [Wigner (1939)](https://www.math.utoronto.ca/mgualt/courses/25-QM/docs/Wigner-1939.pdf).

**Handedness** describes the sense of rotation. “Right” and “left” need a viewing direction and time/phase convention: looking towards the source reverses the apparent rotation seen looking along propagation. Here the unambiguous states are \(J_\pm=(1,\pm i)^T/\sqrt2\), with the exponential convention above. The labels do not mean displacement into an extra spatial dimension. [Beth (1936)](https://doi.org/10.1103/PhysRev.50.115), [Wigner (1939)](https://www.math.utoronto.ca/mgualt/courses/25-QM/docs/Wigner-1939.pdf).

**Orbital angular momentum (OAM)** concerns spatial structure. A paraxial vortex mode with phase \(e^{i\ell\phi}\), for integer \(\ell\), has axial orbital angular momentum \(\ell\hbar\) per photon in the quantum description, or angular momentum/energy \(\ell/\omega\) classically. It can be linearly polarized and carry OAM while having zero mean axial spin. Spatial OAM states allow a larger state space than polarization's two dimensions. In structured or strongly focused fields, spin and orbital contributions can exchange; the simple paraxial formulas need their scope. [Allen et al. (1992)](https://doi.org/10.1103/PhysRevA.45.8185), [D'Ambrosio et al. (2013)](https://arxiv.org/abs/1306.1606), [Mair et al. (2001)](https://arxiv.org/abs/quant-ph/0104070).

### The quantum facts

**Quantum reason:** the photon is a **massless spin-1** particle. Its physical polarization states have helicity \(+1\) or \(-1\); there is no independent longitudinal helicity-zero photon. A massive spin-1 particle instead has three spin projections, \(-1,0,+1\), in its rest frame. Masslessness together with the photon's spin/helicity content explains the two-state quantum description; masslessness alone does not imply two states for every kind of particle. [Wigner (1939)](https://www.math.utoronto.ca/mgualt/courses/25-QM/docs/Wigner-1939.pdf).

For one photon in a fixed spatial/frequency mode, a pure polarization state is
\(|\psi\rangle=\alpha|H\rangle+\beta|V\rangle\), with \(|\alpha|^2+|\beta|^2=1\); a common phase is irrelevant. A mixed state uses a \(2\times2\) density matrix. Its normalized Stokes vector is the polarization qubit's Bloch vector. [James et al. (2001)](https://arxiv.org/pdf/quant-ph/0103121).

An ideal polarization measurement in an orthogonal basis gives **one outcome**, with probability equal to the squared overlap with that basis state. For linear preparation at angle \(\theta\) and analysis at angle \(a\), the probabilities are \(\cos^2(\theta-a)\) and \(\sin^2(\theta-a)\). Many photons give those frequencies. A photon is not detected as two fractional photons. An absorbing polarizer can remove it; a polarizing beam splitter routes its amplitudes into two paths, and detection yields a discrete event. Malus' intensity law also holds classically, so the cosine-squared curve alone does not prove photon quantization. Single-photon anticorrelation supplies different evidence. [D'Ambrosio et al. (2013)](https://arxiv.org/pdf/1306.1606), [Grangier, Roger & Aspect (1986)](https://doi.org/10.1209/0295-5075/1/4/004).

**Polarization entanglement** describes a joint state that cannot be factored into separate photon states, for example \((|HH\rangle+|VV\rangle)/\sqrt2\). Each photon alone is maximally mixed, yet their joint measurements have basis-dependent correlations. Classical correlated mixtures can imitate some same-basis correlations; polarization Bell tests distinguish the stronger quantum correlations under their stated locality and setting assumptions. Such correlations do not let one observer choose the other's local outcome or send a controllable message through it. Entanglement is a quantum-state fact, not a consequence of drawing two classical ellipses. [Giustina et al. (2015)](https://arxiv.org/abs/1511.03190), [Shalm et al. (2015)](https://arxiv.org/abs/1511.03189).

## 2. What experiments show

Each E-number is one experimental-fact entry and has one classification in part 4. E32 separates the longitudinal-mode issue using the evidence in E01–E03 and E33–E34; it is not another independent experiment. “Precision” below distinguishes an uncertainty, an exclusion bound, an instrumental sensitivity, and an experimental range. Where a historical demonstration or opened abstract supplies no scalar uncertainty, that limitation is explicit. The selection is representative, not a catalogue of every polarization experiment.

Citation access: **F** = primary full text opened, with relevant sections inspected; **A** = primary abstract/summary opened; **S** = primary publisher search extract inspected but original page not successfully opened; **R** = citation found in references/search metadata, original not opened. F does not mean every page was read. Preprint versions are identified where their numbers matter.

| ID | What was measured | Who, when; precision or bound | What the result establishes, and qualifications | Primary citation; access |
|---|---|---|---|---|
| E01 | Polarization tomography using independent linear and circular analyzer settings | James, Kwiat, Munro & White, 2001: reconstruction in a two-dimensional single-photon polarization space; four Stokes quantities including intensity. Their two-photon example uses 16 projections. **No universal bound on an arbitrary third mode is reported here.** | Operationally, ordinary polarization is described by two basis amplitudes, with superpositions and mixtures. Tomography assuming this space is not a model-independent census of every possible weakly coupled field. | [*Measurement of qubits*, DOI 10.1103/PhysRevA.64.052312](https://doi.org/10.1103/PhysRevA.64.052312); [paper, §II and experimental example](https://arxiv.org/pdf/quant-ph/0103121). F |
| E02 | Departure from Coulomb's law, interpreted as photon mass with Proca equations | Williams, Faller & Hill, 1971: \(\mu^2=(1.04\pm1.2)\times10^{-19}\ \mathrm{cm}^{-2}\); equivalently force exponent \(q=(2.7\pm3.1)\times10^{-16}\) in \(r^{-(2+q)}\). Here \(\mu=m_\gamma c/\hbar\). | A laboratory null test constrains a **specified massive-vector theory**. It does not measure the amplitude of every conceivable longitudinal material wave. The dimensionful result is an inverse Compton length squared, not kilograms squared. | [*New Experimental Test of Coulomb's Law*, DOI 10.1103/PhysRevLett.26.721](https://doi.org/10.1103/PhysRevLett.26.721). A |
| E03 | Solar-wind currents and magnetic-field curl as a test of Maxwell–Ampère/Proca dynamics | Retinò, Spallicci & Vaivads, 2016, Cluster data: estimated upper mass bounds range from \(1.4\times10^{-49}\) to \(3.4\times10^{-51}\ \mathrm{kg}\), depending on analysis assumptions. | These are analysis-dependent photon-mass limits, not one universal confidence interval. The authors question whether some stronger solar-wind limits are justified. A massive vector's additional longitudinal polarization and an independent scalar medium mode are different hypotheses. | [DOI 10.1016/j.astropartphys.2016.05.006](https://doi.org/10.1016/j.astropartphys.2016.05.006); [primary abstract](https://arxiv.org/abs/1302.6168). A |
| E04 | Mechanical torque from reversing circular polarization through a birefringent plate | Richard A. Beth, 1936: approximately 120 determinations by two independent observers, with magnitude and sign agreeing with the angular-momentum prediction. No single fractional error is assigned in this survey. | Circular light transfers angular momentum. Reversing helicity changes axial angular momentum by \(2\hbar\) per photon in the quantum interpretation, giving torque \(2P/\omega\) for ideal power \(P\). The mechanical experiment itself is a classical-wave torque test, not a measurement of photon-number discreteness. | [*Mechanical Detection and Measurement of the Angular Momentum of Light*, DOI 10.1103/PhysRev.50.115](https://doi.org/10.1103/PhysRev.50.115); [original scan](https://dielslab.unm.edu/sites/default/files/Beth36.pdf). F |
| E05 | Rotation of trapped absorbing particles driven by an optical vortex | He, Friese, Heckenberg & Rubinsztein-Dunlop, 1995: rotation with a linearly polarized doughnut beam carrying phase winding. Qualitative angular-momentum transfer demonstration; no calibrated per-photon uncertainty is quoted in the opened abstract. | Spatial OAM produces torque without circular polarization. It is not an additional polarization state. The \(\ell\hbar\) interpretation is developed in Allen et al.'s primary mode analysis. | [He et al., DOI 10.1103/PhysRevLett.75.826](https://doi.org/10.1103/PhysRevLett.75.826), A; [Allen et al. (1992), DOI 10.1103/PhysRevA.45.8185](https://doi.org/10.1103/PhysRevA.45.8185), A |
| E06 | Polarization analyzer count probabilities for heralded single photons, with polarization-only controls and spin/OAM “gears” | D'Ambrosio et al., 2013: ordinary \(\cos^2\theta\) control and amplified \(\cos^2(m\theta)\) fringes, with \(m\) up to 100; reported experimental fringe visibilities all exceed 0.73. | Repeated single-photon outcomes follow polarization-overlap statistics. The gear increases the rotation-dependent phase through spatial structure, not the number of intrinsic polarization basis states. Finite visibility is a measured imperfection, not a claimed exact probability curve. | [*Photonic polarization gears*, DOI 10.1038/ncomms3432](https://doi.org/10.1038/ncomms3432); [paper, Fig. 2 and pp. 5–6](https://arxiv.org/pdf/1306.1606). F |
| E07 | Coincidences at the two outputs of a beam splitter for a heralded photon; accompanying interference | Grangier, Roger & Aspect, 1986: anticorrelation parameter \(\alpha=0.18\pm0.06\), compared with classical-wave inequality \(\alpha\ge1\); herald-conditioned fringe visibility above 98%. | Evidence for indivisible detection events together with wave interference. This is load-bearing for the meaning of a “single photon”; the splitter experiment is not itself a polarization Bell test. | [DOI 10.1209/0295-5075/1/4/004](https://doi.org/10.1209/0295-5075/1/4/004); [primary scan](https://courses.physics.illinois.edu/phys513/sp2016/reading/week1/GrangierSinglePhoton1986.pdf). F |
| E08 | Bell correlations of polarization-entangled photon pairs, with high detection efficiency and spacelike separation | Giustina et al., 2015: local-realist null probability \(p\le3.74\times10^{-31}\), described as 11.5 standard deviations. | Simultaneous closure of the principal detection and locality loopholes. The conclusion retains the experiment's assumptions concerning setting independence and the tested local-causal class. Two classical wave components alone do not supply these joint detector statistics. | [DOI 10.1103/PhysRevLett.115.250401](https://doi.org/10.1103/PhysRevLett.115.250401); [primary abstract](https://arxiv.org/abs/1511.03190). A |
| E09 | Independent polarization Bell experiment with efficient photon detection and separated setting choices | Shalm et al., 2015: \(p=5.9\times10^{-9}\) before the setting-predictability correction; \(p=2.3\times10^{-7}\) with it. | Rejects the tested local-realist models without fair-sampling postselection. These p-values are not probabilities that quantum mechanics is true, and the experiment does not select a microscopic model of a brane. | [DOI 10.1103/PhysRevLett.115.250402](https://doi.org/10.1103/PhysRevLett.115.250402); [primary abstract](https://arxiv.org/abs/1511.03189). A |
| E10 | Survival of gamma-ray linear polarization over a cosmological distance, testing energy-dependent helicity birefringence | Götz et al., 2014, GRB 140206A: redshift \(z=2.739\); polarization fraction above 28% at 90% confidence; dimension-five Lorentz-violating birefringence parameter bounded at \(\xi<10^{-16}\) in their convention. | Different circular propagation phases would rotate and, across a bandwidth, wash out linear polarization. The bound depends on the dispersion parametrization and source/polarization assumptions. It is **not** a universal bound \(|c_R-c_L|/c<10^{-16}\). | [DOI 10.1093/mnras/stu1634](https://doi.org/10.1093/mnras/stu1634); [primary abstract](https://arxiv.org/abs/1408.4121). A |
| E11 | Wavelength-dependent polarization changes from distant optical sources, testing vacuum Lorentz violation | Kostelecký & Mewes, 2002: birefringent dimension-four Standard-Model Extension coefficient combinations constrained at the \(2\times10^{-32}\) level. | A model-specific bound on birefringent coefficients, including directional structure. It is not a single bound on all possible polarization-dependent group velocities or all terms in a different theory. | [DOI 10.1103/PhysRevD.66.056005](https://doi.org/10.1103/PhysRevD.66.056005); [primary abstract](https://arxiv.org/abs/hep-ph/0205211). A |
| E12 | Spatially varying rotation of CMB polarization | Bianchini et al., SPTpol, 2020: 500 square degrees at 150 GHz; scale-invariant rotation-power amplitude \(L(L+1)C_L^{\alpha\alpha}/(2\pi)<0.10\times10^{-4}\ \mathrm{rad}^2=0.033\ \mathrm{deg}^2\), 95% confidence. | Null result for **anisotropic** cosmic birefringence in the tested spectrum. It does not exclude a uniform rotation angle or every other angular spectrum. | [DOI 10.1103/PhysRevD.102.083504](https://doi.org/10.1103/PhysRevD.102.083504); [primary abstract](https://arxiv.org/abs/2006.08061). A |
| E13 | Uniform CMB polarization rotation, separating instrumental angle from foreground emission | Minami & Komatsu, Planck, 2020: \(\beta=0.35\pm0.14^\circ\), 68% confidence; 2.4σ preference. | **Reported hint, not established cosmic birefringence.** Its evidence is a nonzero fitted rotation under the foreground/calibration treatment. A rotation angle is not directly a photon arrival-speed measurement. | [DOI 10.1103/PhysRevLett.125.221301](https://doi.org/10.1103/PhysRevLett.125.221301); [primary abstract](https://arxiv.org/abs/2011.11254), A; [primary full text, abstract and introduction](https://arxiv.org/pdf/2011.11254), F. |
| E14 | Uniform rotation with newer Planck PR4 maps and foreground/mask tests | Diego-Palazuelos et al., 2022: initial nearly full-sky \(\beta=0.30\pm0.11^\circ\), 68% confidence; the estimate decreases with more restrictive masks. | **Reported hint, not established.** The authors explicitly withhold cosmological significance until polarized foreground emission is understood better. Shares Planck data with E13; not an independent replication. | [Authors' primary report](https://arxiv.org/abs/2203.04830), A; [primary full text, §1](https://arxiv.org/pdf/2203.04830), F. |
| E15 | Uniform cosmic rotation in ACT DR6, including calibration/systematic assessment | Diego-Palazuelos & Komatsu, 2025 preprint, **revised 2026-04-14, v2**: their Bayesian estimate is \(\beta=0.215\pm0.074^\circ\), 68% confidence, 2.9σ preference, with the baseline correlated calibration prior and residual leakage marginalized. The same paper describes the ACT collaboration's \(0.20\pm0.08^\circ\) as its **mean detector rotation** \(\langle\psi_i\rangle\), derived from EB with optics-calibration uncertainty in a frequentist analysis; ACT did not explicitly estimate \(\beta\) there. | **Reported hint, not established.** Unresolved instrumental systematics limit cosmological conclusions. The two numbers are different estimands and analyses, not earlier and revised values of these authors' \(\beta\). | [Primary v2 full text, §I, Table III and §VI](https://arxiv.org/html/2509.13654v2). F. The ACT-team comparison was read in this paper's account of its Ref. 27; that referenced original was not opened. |
| E16 | Rotation of linear polarization transmitted through magnetized material | Michael Faraday, discovery 1845, publication 1846: demonstrated magnetic rotation in transparent matter. Historical demonstration; this survey does not extract a modern Verdet-constant uncertainty or field threshold from the original. | **Faraday rotation in matter is established.** It concerns material response and opposite circular propagation phases; it should not be labelled rotation in empty vacuum. | [*Experimental Researches in Electricity.—Nineteenth Series*, DOI 10.1098/rstl.1846.0001](https://doi.org/10.1098/rstl.1846.0001). R: historical citation located; original fetch unsuccessful |
| E17 | Electric-field-induced double refraction in dielectric material | John Kerr, 1875: glass/dielectric medium becomes birefringent when electrified. Historical demonstration; no scalar precision or modern electric-field threshold is extracted here. | **Electro-optic birefringence in matter is established.** This row is the electric-field Kerr effect, rather than magnetic reflection rotation or intensity-driven optical Kerr measurements. It is not an empty-vacuum detection. | [*A new relation between electricity and light: Dielectrified media birefringent*, DOI 10.1080/14786447508641302](https://doi.org/10.1080/14786447508641302). S: publisher title/date extract; original full text not opened |
| E18 | Laboratory magnetic-vacuum birefringence and dichroism | Ejlli et al., PVLAS final results, 2020: at \(B=2.5\ \mathrm T\), \(\Delta n=(12\pm17)\times10^{-23}\); \(|\Delta\kappa|=(10\pm28)\times10^{-23}\). Predicted QED \(\Delta n\) at that field is \(2.5\times10^{-23}\). | **Null result at this sensitivity.** The uncertainty remains about seven times the QED prediction; this does not rule that prediction out. The laboratory observable is relative polarization phase/absorption, not automatically a loss of total light into another medium mode. | [*The PVLAS experiment*, DOI 10.1016/j.physrep.2020.06.001](https://doi.org/10.1016/j.physrep.2020.06.001); [primary abstract](https://arxiv.org/abs/2005.12913). A |
| E19 | Previously reported magnetic-vacuum optical rotation and its subsequent instrumental reanalysis | PVLAS, 2006: claimed \((3.9\pm0.5)\times10^{-12}\ \mathrm{rad/pass}\) at 5 T. Follow-up, 2008, excluded the earlier signal and identified instrumental artifacts. | **Withdrawn signal, not an established new vacuum effect.** It must not become a requirement that the toy model reproduce the original central value. | [Original DOI 10.1103/PhysRevLett.96.110406](https://doi.org/10.1103/PhysRevLett.96.110406), A; [follow-up DOI 10.1103/PhysRevD.77.032006](https://doi.org/10.1103/PhysRevD.77.032006), A |
| E20 | Optical polarization of isolated neutron star RX J1856.5−3754 | Mignani et al., 2017 publication, 2016 preprint: polarization degree \(16.43\pm5.26\%\); angle \(145.39\pm9.44^\circ\). | **Measured stellar polarization; model-dependent evidence for QED vacuum birefringence.** The inference uses surface emission, magnetic geometry and propagation models. It is not a direct refractive-index measurement or a standalone established vacuum-birefringence detection. | [DOI 10.1093/mnras/stw2798](https://doi.org/10.1093/mnras/stw2798); [primary abstract](https://arxiv.org/abs/1610.08323). A |
| E21 | Energy-dependent X-ray polarization of magnetar 4U 0142+61 | Taverna et al., **published version of record, 2022**, reports “a linear polarization degree of 13.5 ± 0.8%” averaged over 2–8 keV. The published band values are \(15.0\pm1.0\%\) at 2–4 keV and \(35.2\pm7.1\%\) at 5.5–8 keV; polarization falls below instrumental sensitivity around 4–5 keV, where the angle changes by roughly \(90^\circ\). | **Polarization structure observed; vacuum-birefringence interpretation remains indirect.** The published abstract reports consistency with surface thermal radiation reprocessed by magnetospheric charged-particle scattering. Emission and propagation modelling enter the interpretation. | [Version of record, DOI 10.1126/science.add0080](https://doi.org/10.1126/science.add0080); [published abstract in the authors' Caltech repository](https://authors.library.caltech.edu/records/q5fs0-fwh71). A: published abstract opened; publisher full text not opened. |
| E22 | Survival of GRB linear polarization as an inferred helicity-dependent gravitational delay test | Yang et al., 2017, GRB 110721A: claimed PPN difference \(\Delta\gamma_p<1.6\times10^{-27}\) between circular polarizations under their gravitational-potential/emission assumptions. | An **inferred and theory-dependent equivalence-principle bound**, not separately timed right/left pulses. Shapiro-delay comparisons need a specified observable and alternative theory; generic conversions can be problematic. This is not a universally established \(10^{-27}\) limit on every gravitational polarization effect. | [DOI 10.1093/mnrasl/slx045](https://doi.org/10.1093/mnrasl/slx045); [primary abstract](https://arxiv.org/abs/1706.00889), A; qualification: [Minazzoli et al. (2019)](https://arxiv.org/abs/1907.12453), A |
| E23 | Solar gravitational deflection of radio-source positions | Fomalont, Kopeikin, Lanyi & Benson, 2009 report of 2005 VLBA observations: PPN \(\gamma=0.9998\pm0.0003\), 68% confidence. | Precise **common deflection** test. The quoted parameter is not a right-minus-left or linear-polarization differential deflection measurement. Common lensing success does not establish polarization independence at this uncertainty. | [DOI 10.1088/0004-637X/699/2/1395](https://doi.org/10.1088/0004-637X/699/2/1395); [primary abstract](https://arxiv.org/abs/0904.3992). A |
| E24 | Solar gravitational time delay through spacecraft radio frequency shifts | Bertotti, Iess & Tortora, Cassini, 2003: \(\gamma-1=(2.1\pm2.3)\times10^{-5}\). | Precise **common Shapiro-delay** test, not a quoted helicity-differential delay constraint. The number must not be repurposed as \(|\gamma_R-\gamma_L|\). | [*A test of general relativity using radio links with the Cassini spacecraft*, DOI 10.1038/nature01997](https://doi.org/10.1038/nature01997). S: primary abstract search extract; publisher page fetch unsuccessful |
| E25 | Spin-dependent transverse beam displacement on refraction at an air/glass interface | Hosten & Kwiat, 2008: weak-measurement method with approximately ångström displacement sensitivity and nearly \(10^4\) amplification. | **Optical spin Hall effect observed in a material-interface experiment.** Sensitivity is not the value of every measured shift. This is not experimental detection of a gravitational spin Hall effect. | [DOI 10.1126/science.1152697](https://doi.org/10.1126/science.1152697); [primary paper](https://research.physics.illinois.edu/QI/Photonics/papers/My%20Collection.Data/PDF/Observation%20of%20the%20spin%20hall%20effect%20of%20light%20via%20weak%20measurements.pdf). F |
| E26 | Polarization of light reflected from transparent surfaces | David Brewster, 1815: polarizing incidence obeys \(\tan\theta_B=n_2/n_1\); historical measurements across materials. No modern angular uncertainty is extracted here. | Reflection preferentially selects polarization; at Brewster incidence an ideal dielectric interface's reflected monochromatic light is linearly polarized perpendicular to the incidence plane. The historical account also discusses dispersion and residual unpolarized white light. | [DOI 10.1098/rstl.1815.0010](https://doi.org/10.1098/rstl.1815.0010); [primary Royal Society proceedings summary, transcribed](https://en.wikisource.org/wiki/Proceedings_of_the_Royal_Society_of_London/Volume_2/On_the_laws_which_regulate_the_polarization_of_light_by_reflection_from_transparent_bodies). A, summary opened; full original not opened |
| E27 | Degree and angle of scattered skylight polarization | Pomozi, Horváth & Wehner, 2001: full-sky \(180^\circ\) imaging polarimetry at 450, 550 and 650 nm, under clear and partly cloudy conditions; compared with single-scattering Rayleigh patterns. No single scalar instrumental uncertainty is extracted from the available summary. | Scattering generates systematic linear polarization patterns, with cloud-dependent degree and angle structure. Actual sky polarization is not assumed to equal an ideal single-scattering formula everywhere. | [DOI 10.1242/jeb.204.17.2933](https://doi.org/10.1242/jeb.204.17.2933). S: primary publisher summary inspected in search output; original page not successfully opened |
| E28 | Polarization of the cosmic microwave background produced by scattering | Kovac et al., DASI, 2002: E-mode polarization detected at 4.9σ. | Observed CMB polarization agrees with the expected scattering origin and acoustic structure. Here “E-mode” is a sky-pattern decomposition, not a third electric-field polarization direction or the bulk coordinate w. | [DOI 10.1038/nature01269](https://doi.org/10.1038/nature01269); [primary abstract](https://arxiv.org/abs/astro-ph/0209478). A |
| E29 | Vector-field focusing of a radially polarized beam, including an axial electric component | Dorn, Quabis & Leuchs, 2003: measured focused spot area \(0.16(1)\lambda^2\), compared with the theoretical linear-polarization limit \(0.26\lambda^2\) under the stated conditions. | The focus can have electric field along its **mean optical axis**. This does not establish an independent longitudinal plane-wave photon: the beam combines nonparallel wavevectors. Its axial component must not be automatically identified with the toy model's independent compressional \(u_L\) branch. | [DOI 10.1103/PhysRevLett.91.233901](https://doi.org/10.1103/PhysRevLett.91.233901); [primary paper](https://arxiv.org/pdf/physics/0310007). F |
| E30 | Quantum correlations and entanglement in spatial OAM states | Mair, Vaziri, Weihs & Zeilinger, 2001: demonstrated entanglement in OAM mode space of dimension greater than two. No single global error or Bell significance is quoted from the opened abstract. | Photon state space extends beyond polarization when spatial modes are included. Higher-dimensional OAM correlations do not imply more than two intrinsic polarization basis states per momentum. | [DOI 10.1038/35085529](https://doi.org/10.1038/35085529); [primary abstract](https://arxiv.org/abs/quant-ph/0104070). A |
| E31 | Joint X-ray/radio polarization of magnetar 1E 1547.0−5408 | Stewart et al., Nature 2026; March–April 2025 observations, arXiv v5 revised 2026-08-17: X-ray 2–8 keV degree \(46\pm4\%\); fitted soft-band degree \(65\pm8\%\) at 2 keV, 1σ. | **Strong, model-dependent evidence, not an established direct vacuum-birefringence measurement.** Evidence includes high polarization, its energy/phase dependence, radio/X-ray geometry and worse Stokes fits in the authors' tested propagation models with vacuum birefringence off. The paper itself continues to describe the QED prediction as unconfirmed. | [DOI 10.1038/s41586-026-10859-z](https://doi.org/10.1038/s41586-026-10859-z); [primary v5, main text and methods](https://arxiv.org/html/2509.19446v5). F; Nature version-of-record fetch unsuccessful |
| E32 | Constraint on an additional longitudinal radiating polarization, separated from the in-plane count | Experimental basis is E01's tomography, E02–E03's specified-theory null tests, and **E33–E34's conditional thermal count of radiating modes**. James et al. (2001) describe “the two polarization degrees of freedom” in their experimental reconstruction. **No universal third-mode coupling/amplitude bound follows.** | Ordinary photon polarization is two-dimensional. The laboratory blackbody power and CMB spectrum in E33–E34 also constrain an extra branch **if it participates as an additional ordinary thermal radiation mode under their occupation, dispersion and coupling assumptions**. This row separates the longitudinal departure; it is not a new independent measurement. | [James et al.](https://arxiv.org/pdf/quant-ph/0103121), F; [Williams et al.](https://doi.org/10.1103/PhysRevLett.26.721), A; [Retinò et al.](https://arxiv.org/abs/1302.6168), A; **E33:** [Quinn & Martin](https://doi.org/10.1098/rsta.1985.0058), A; **E34:** [Fixsen et al.](https://arxiv.org/pdf/astro-ph/9605054), F. See those rows for assumptions and access details. |
| E33 | Laboratory blackbody total radiant exitance and its absolute Stefan–Boltzmann normalization | Quinn & Martin, 1985: cryogenic radiometry at 273.16 K and approximately 233–373 K; reported \(\sigma=(5.66967\pm0.00076)\times10^{-8}\ \mathrm{W\,m^{-2}\,K^{-4}}\), in the measurement's original unit realization. The quoted uncertainty is about \(1.34\times10^{-4}\) relative. | **Conditional thermal count:** absolute power agrees with the two-polarization Planck normalization within the reported measurement/calculated-constant uncertainties. Inferring a radiating count requires equilibrium occupation, the usual dispersion and energy per quantum, and known temperature/emissivity/coupling. It is not a count of decoupled material modes. | [Quinn & Martin, DOI 10.1098/rsta.1985.0058](https://doi.org/10.1098/rsta.1985.0058), A: publisher-deposited primary abstract [opened through Crossref](https://api.crossref.org/works/10.1098/rsta.1985.0058); full text not opened. The mode-count interpretation uses [Planck (1901), §8, eqs. (8), (11)–(12), primary-paper translation](https://jgluza.us.edu.pl/QM/planck1901.pdf), F. |
| E34 | Cosmological blackbody radiance spectrum of the CMB | Fixsen et al., COBE/FIRAS, 1996: \(T=2.728\pm0.004\ \mathrm K\), 95% confidence; RMS spectral deviations below **50 parts per million of the peak CMB radiance**. | **Conditional thermal count:** the spectrum fits the ordinary Planck law with its two-polarization normalization. FIRAS compares the sky with a calibrated blackbody; its residual is a spectral-distortion limit, **not** a universal 50-ppm limit on the population/coupling of any extra mode. Temperature, calibration, occupation, dispersion and how a mode couples to sky, reference and detector matter for that inference. | [*The Cosmic Microwave Background Spectrum from the Full COBE FIRAS Data Set*, DOI 10.1086/178173](https://doi.org/10.1086/178173); [primary paper, abstract and §6](https://arxiv.org/pdf/astro-ph/9605054). F. Mode-count interpretation: [Planck (1901)](https://jgluza.us.edu.pl/QM/planck1901.pdf), F. |
| E35 | Solar-wind dynamics at the scale of Pluto's orbit, interpreted as a photon-mass limit; PDG's adopted limit | Ryutov, 2007; **Particle Data Group, 2026 listing:** adopted \(m_\gamma<1\times10^{-18}\ \mathrm{eV}/c^2\) (approximately \(1.8\times10^{-54}\ \mathrm{kg}\)). The listing assigns **no confidence percentage** to this estimate. | **Assumptions:** Proca electrodynamics and its quasistatic, single-fluid MHD reduction; observed large-scale solar-wind magnetic field, density and largely radial flow around 40 AU. The opened author report bounds the additional magnetic force relative to the flow's inertia, with a safety factor \(q=3\). It is a dynamical/force-balance bound, not a model-independent radiation-state census. | [PDG 2026 photon listing, mass table and note 1](https://pdg.lbl.gov/2026/listings/rpp2026-list-photon.pdf), F; [Ryutov, DOI 10.1088/0741-3335/49/12B/S40](https://doi.org/10.1088/0741-3335/49/12B/S40); [primary author report, §§2–3](https://www.osti.gov/servlets/purl/940862), F: PDF opened and text extracted in memory. The numerical adopted limit here is PDG's, not an unlabelled substitution from the author-report draft. |
| E36 | Frequency-dependent arrival delays/dispersion measures of localized fast radio bursts | Lin, Tang & Zou, 2023, 17 FRBs: \(m_\gamma<4.8\times10^{-51}\ \mathrm{kg}\) at 1σ and \(<7.1\times10^{-51}\ \mathrm{kg}\) at 2σ, for **no redshift evolution of the host-DM distribution**. Allowing that evolution gives \(<6.7\times10^{-51}\) and \(<1.01\times10^{-50}\ \mathrm{kg}\), respectively. | **Assumptions:** massive-photon dispersion plus plasma dispersion, flat ΛCDM with Planck-2018 parameters, modelled IGM fluctuations and a lognormal host-DM distribution; Galactic DM uses NE2001 and halo DM is fixed at \(50\ \mathrm{pc\,cm^{-3}}\). Both mass and plasma delays scale as \(\nu^{-2}\) here; their redshift dependences and the adopted distributions underpin separation. These are Bayesian analysis-dependent bounds, not a detection or a handedness-speed limit. | [DOI 10.1093/mnras/stad228](https://doi.org/10.1093/mnras/stad228); [primary preprint v1, §§2–3 and Table 2](https://arxiv.org/html/2301.12103v1), F. Publisher full-text fetch unsuccessful; stated assumptions/numbers read in the opened primary preprint. |
| E37 | FRB dispersion-mass limit using an expansion history reconstructed from cosmic-chronometer data | Ran, Wang & Wei, 2024: 32 localized FRBs and 34 \(H(z)\) measurements; \(m_\gamma\le3.5\times10^{-51}\ \mathrm{kg}\) at 1σ and \(\le6.5\times10^{-51}\ \mathrm{kg}\) at 2σ, for the paper's Gaussian Galactic-halo prior. | **Assumptions:** plasma-plus-massive-photon delays, FLRW geometry, a smoothed ANN reconstruction of \(H(z)\), and modelled IGM/host DM distributions. The mass prior is uniform on \(0\)–\(10^{-49}\ \mathrm{kg}\); the halo prior is \(65\pm15\ \mathrm{pc\,cm^{-3}}\), restricted to 20–110. “Cosmology-independent” means no fixed ΛCDM expansion history in this fit, not assumption-free. With a flat halo prior the quoted limits become \(3.8\) and \(7.2\times10^{-51}\ \mathrm{kg}\). | [Primary paper v1, §§II–III and Table 3](https://arxiv.org/html/2404.17154v1), F; [primary abstract](https://arxiv.org/abs/2404.17154), A. |
| E38 | Cosmological effective count of relativistic species, \(N_{\rm eff}\), inferred from CMB anisotropies and BAO | Planck Collaboration, **Planck 2018 results. VI. Cosmological parameters** (2020 publication; opened preprint v4, 2021): **\(N_{\rm eff}=2.99\pm0.17\)**, 68% confidence, for **TT,TE,EE+lowE+lensing+BAO**; eq. (67b) gives **\(2.99^{+0.34}_{-0.33}\)** at 95%. | **Assumptions:** spatially flat ΛCDM with \(N_{\rm eff}\) as the single added parameter, adiabatic power-law primordial perturbations, baseline neutrino mass \(\sum m_\nu=0.06\ \mathrm{eV}\), standard recombination and BBN helium prediction, and noninteracting/nondecaying extra radiation. This counts **non-photon relativistic energy density in neutrino units**, not photon polarization states; it is consistent with 3.046. An extra branch constrains this count only through its cosmological population, dispersion, gravitating energy/perturbations and temperature/decoupling history; a mode census alone predicts no \(\Delta N_{\rm eff}\). | [Primary Planck paper, abstract, §§2.1, 7.5.2, 7.6.2 and §8](https://arxiv.org/pdf/1807.06209), F; [publication DOI 10.1051/0004-6361/201833910](https://doi.org/10.1051/0004-6361/201833910). The connection to an extra brane branch is a conditional survey inference, not a Planck test of this model. |

**PDG's nonzero FRB mass entry — not established.** The opened [PDG 2026 photon listing, mass table and note 3](https://pdg.lbl.gov/2026/listings/rpp2026-list-photon.pdf), F, also lists **LEMOS 25: \(m_\gamma=(16^{+3}_{-9})\times10^{-15}\ \mathrm{eV}/c^2\)** from fast radio bursts. PDG places it among results it excludes from averages, fits and limits; note 3 says the FRB/supernova inference uses several assumptions and gives no upper limit. The listing assigns no confidence percentage to this value. This is a reported nonzero value, not an established photon-mass detection; it is mentioned here as requested, without adding a separate experimental/classification row. LEMOS 25's original paper was not opened; the value and qualifications here were read in PDG.

**What “thermal count” means in E33–E34.** This is an inference from the Planck-law normalization, not a detector labelling each polarization. For \(g\) equally coupled radiating states with \(\omega=ck\), quantum energy \(h_{\rm P}\nu\), and thermal occupation at the same \(T\), the spectral radiance is \(B_{\nu,g}=g h_{\rm P}\nu^3/[c^2(e^{h_{\rm P}\nu/(k_{\rm B}T)}-1)]\). The usual photon law has \(g=2\); a third state under **all the same conditions** would multiply that normalization by \(3/2\). This mode-count extension is the survey's inference from [Planck's primary law, eqs. (8), (11)–(12)](https://jgluza.us.edu.pl/QM/planck1901.pdf), not a new bound reported by either measurement. The laboratory absolute-power test and FIRAS sky/reference comparison have the distinct qualifications stated in their rows.

<a id="beta-tests"></a>
**What the cited β papers test.** [Minami & Komatsu (2020), abstract/introduction](https://arxiv.org/pdf/2011.11254), F, search for parity-violating physics using CMB EB correlations and separate the cosmic rotation β from detector angle. [Diego-Palazuelos et al. (2022), §1](https://arxiv.org/pdf/2203.04830), F, describe β as the rotation from a difference between the phase velocities of opposite photon helicities, testing parity-violating physics. [Diego-Palazuelos & Komatsu (ACT v2), abstract/§I](https://arxiv.org/html/2509.13654v2), F, use parity-odd EB/TB correlations to estimate uniform cosmological rotation and test cosmological parity violation. These interpretations retain the calibration/foreground/systematic qualifications in E13–E15; β is not a separately measured photon arrival-speed difference.

**What “speed” means here.** Birefringence is a difference in propagation eigenmodes' phase accumulation. Frequency-dependent phase differences can imply different group delays and depolarize a broadband beam. A nearly frequency-independent fitted cosmic rotation is not, by itself, a measurement of an arrival-time difference. The SME/GRB parameters, anisotropic CMB power, and uniform CMB angle therefore remain separate observables. E10–E15 give the relevant primary citations; they do not supply one universal polarization-speed bound.

**Gravity and spin Hall effects.** Finite-wavelength gravitational spin Hall effects are a theoretical wave-propagation subject, distinct from leading geometric-optics lensing and from the measured air/glass effect in E25. [Oancea et al. (2020), DOI 10.1103/PhysRevD.102.024075](https://doi.org/10.1103/PhysRevD.102.024075), [primary abstract](https://arxiv.org/abs/2003.04553), A, derives polarization-dependent ray corrections in curved spacetime. The sources opened for this survey provide no established dedicated gravitational spin Hall measurement to enter as an experimental detection. That is a statement about this survey's retrieved evidence, not a universal experimental exclusion.

**Theory-citation access for part 1.** Jones: publisher abstract opened (A); Wolf: publisher summary opened (A); James: primary PDF opened (§II and relevant tomography discussion, F); Wigner: primary paper PDF opened (F); Allen: publisher abstract opened (A). Stokes: primary collected-paper chapter preview opened (A), not its complete text. Maxwell: original DOI identified but publisher fetch unsuccessful (R); the plane-wave relations above are the explicit vacuum-Maxwell specialization. Other part-1 sources have the access status recorded in E04–E09 and E29–E30. No unopened source is represented as a full-paper reading.

## 3. What the repository says

The following is an inventory of substantive records, statements and historical leads retrieved for this survey. It does not promote executable names, search hits or an old paper's “earned” label to a current v3 result. Every repository assertion is tied to a numbered evidence item in Appendix A, which contains the command and its literal output. The authority hierarchy is the one specified in the user's task.

The current substrate register says **“Sixteen entries, all OPEN”** and says the light sector rests on S1–S8, **“none of which have been run”**: **research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:3–19** [R28](#r28). Thus a computed spectrum of a supplied action and a physically supplied substrate are different statuses in the repository.

For the requested geometry, the distinction is explicit: S10 supplies a D-component **in-plane** displacement and inherits its separation from out-of-plane fields; S11b supplies a finite thickness W along w; R-S8-02 asks for the full in-plane u plus normal h operator. With the stipulated three-dimensional tangent space and propagation along x, the in-plane transverse displacement directions are y and z. A w-directed normal displacement is an additional field to audit, rather than one of those two directions. This is a coordinate restatement of the cited field scopes, not a completed guided-spectrum result: **research/pde_ledger_v3/steps/S10_two_transverse_photons.md:182–197** [R10](#r10); **research/pde_ledger_v3/steps/S11b_interface_coupling_law.md:14–16** [R30](#r30); **research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:263–305** [R28](#r28).

Two meanings of **“transverse”** must be kept separate. A displacement normal to the brane is transverse **to the brane**: the plan identifies its w component as \(\xi_w=\ell h\), with the h-branon as charge field/static mediator, **V3_STEP_PLAN.md:453–462,871–905** [P6](#p6). **Normal displacement also describes the brane's geometry and optical embedding:** O2 supplies \(g_{ij}=\delta_{ij}+\partial_i\xi_w\partial_j\xi_w\), retains \(\xi_w=\ell h\), and includes \((\partial\xi_w)^2\) in the optical counting, **O2_steady_brane_balance.md:89–92,115–117** [Q4](#q4). These are supplied geometric/optical inputs, not an independently derived geometry or a restriction of all w motion to charge. A light displacement transverse **to propagation** satisfies \(k\cdot u_T=0\); the guided-light proposal also imposes tangency \(N\cdot u_T=0\), **docs/light_guided_photon_soliton_research_plan.md:246–274** [R22](#r22). Perpendicular to propagation alone does not exclude w; tangency does. The source statuses are recorded below.

### Current authoritative ledger

| Subject | Source, literal quotation and evidence | Status stated by the source; scope |
|---|---|---|
| Light's empirical target | **research/pde_ledger_v3/steps/S9_light_requires_shear.md:21**: “Two transverse polarisations · one speed · no longitudinal mode.” [R6](#r6) | A target demanded of the light sector, not three experimental results derived from the medium. |
| Scalar-superfluid premise | **research/pde_ledger_v3/steps/S9_light_requires_shear.md:331–332**: “no script executes it” and “It remains a **supplied premise**.” [R6](#r6) | The claim that the scalar GNLS substrate has no transverse mode is supplied in this record; this survey does not upgrade its cited review into a new independent derivation. |
| Curl-only stiffness and bulk shear-freeness | **research/pde_ledger_v3/steps/S9_light_requires_shear.md:333–337**: “The absence of a propagating longitudinal wave is **ASSUMED, not derived**”; “It removes the restoring *force*, … **not the degree of freedom**”; “**Bulk shear-freeness** is postulated”. [R6](#r6) | **Postulated**, with an explicitly surviving longitudinal zero mode. The record's limits include a sharp zero-width sheet, rest background, no dissipation and frequency-independent moduli (lines 349–350). |
| Count actually measured | **research/pde_ledger_v3/steps/S10_two_transverse_photons.md:73–76**: “For the supplied curl-only in-plane action … the nonzero root has D − 1 transverse null directions … D = 2, 3, 4, 5. Thus the D = 3 member has two transverse directions.” [R10](#r10) | **Conditional computed count**, for nonzero wavevector and the stated generic strata. |
| Count is not unique to curl-only elasticity | **research/pde_ledger_v3/steps/S10_two_transverse_photons.md:80–84**: “FULLGRAD has the same nonzero root and the same D − 1 transverse nullity”; curl-only “leaves that direction at the zero root”. [R10](#r10) | **Computed control result**: the transverse count alone does not select the stiffness form; the longitudinal root distinguishes these controls. |
| Isotropy matters | **research/pde_ledger_v3/steps/S10_two_transverse_photons.md:98–101**: MAIN generic propagating count “2 / 2”; ANISO generic “2 / 1”; ANISO perpendicular stratum “2 / 2”. The column is total nullity / exactly-transverse nullity. [R10](#r10) | **Computed qualification**: anisotropic inertia can split propagating branches and reduce exact transversality on generic directions. Special axes differ. |
| D = 3 is not selected | **research/pde_ledger_v3/steps/S10_two_transverse_photons.md:172–176**: “The physical selection D = 3 is not made in S10”; the conditional Lean baseline extends the map to arbitrary finite D; “Neither establishes which D nature selects.” [R10](#r10) | **Conditional theorem/count**, not a dimensionality derivation. |
| Out-of-plane fields excluded from the count | **research/pde_ledger_v3/steps/S10_two_transverse_photons.md:182–197**: “a D-component in-plane displacement u, with its separation from every out-of-plane field inherited rather than tested”; “the real cosine plane-wave ansatz”; “The action is an input”. [R10](#r10) | **Supplied field content and premises**. This record alone does not establish normal-mode exclusion, a full finite-thickness spectrum, circular-wave angular momentum or quantum helicity. |
| Longitudinal zero is a remaining degree of freedom | **research/pde_ledger_v3/steps/S10_two_transverse_photons.md:245–252**: “The zero root retains a degree of freedom; it does not remove one”; “Maxwell has no counterpart to this surviving zero-frequency longitudinal direction”; “a characterised departure”. [R10](#r10) | **Derived within the supplied action**, kept as a departure, not erased by calling the pair “photons”. |
| Compression lifts that root | **research/pde_ledger_v3/steps/S11_stray_longitudinal.md:15–19**: “Lifts S10's zero to a propagating longitudinal mode”, with compression modulus “**postulated**”; “closes the homogeneous three-finite-mode census”. [R12](#r12) | **Computed conditional spectrum**; modulus postulated with a retirement condition. Actual leakage and helper-wave questions remain separately scoped. |
| uL, charge and w motion are not synonyms | **research/pde_ledger_v3/steps/S11_stray_longitudinal.md:23–25**: “does **not** identify \(u_L\) with charge”; “charge is the oriented throat/electric-odd boundary condition, while \(u_L\) is a separate in-plane longitudinal material mode.” **research/pde_ledger_v3/V3_STEP_PLAN.md:413–416** calls the earlier identification with “**±w** displacement” a “false identification”. [R12](#r12), [R19](#r19) | **Explicit corrected ontology**. The old FAIL token is described as misnamed, not a failure verdict on the retained mode. |
| Width change and compressional response | **research/pde_ledger_v3/steps/S11_stray_longitudinal.md:81–92**: “the brane **THICKENS**”; “the sheet bulges into **±w**”; “S11's cone below is computed with the wall width **FROZEN**”. [R12](#r12) | The record distinguishes in-plane strain from accommodating density/width channels. A long-wavelength flexural crossover is conditional future work, not an already obtained polarization spectrum. |
| Scope of homogeneous separation | **research/pde_ledger_v3/steps/S11_stray_longitudinal.md:74–76**: “homogeneous \(D=3\) quadratic transverse and longitudinal eigenbranches have zero linear cross-block under the selected action”; not established “at nonlinear order, on a nonuniform slab, at an interface or defect, or after additional allowed fields are introduced.” [R12](#r12) | **Derived, narrowly scoped decoupling**, not structural independence on every background. |
| Observability of the extra branch | **research/pde_ledger_v3/steps/S11_stray_longitudinal.md:154,161–163** assigns direct matter observation to “derived matter/interface coupling” and expressly does not deliver “observability or unobservability of the longitudinal branch”. [R12](#r12) | **Unresolved**. The ledger supplies neither a detector-level third-photon amplitude nor an observational exclusion of the branch. |
| Finite-thickness interface | **research/pde_ledger_v3/steps/S11b_interface_coupling_law.md:14–16**: “a slab of finite thickness W in the w direction with two faces meeting the bulk”. Lines 74–81: transverse-to-thickness coupling is “identically zero”, with stable real frequency “only where … \(\mu_\perp\ge0\)”; this “does NOT settle unconditional confinement”. [R30](#r30) | **Conditional uniform interface result**, including a sign/stability condition. Face motion and thickness response are explicit; general nonuniform confinement is not delivered. |
| Uniform transverse mode in the assembled interface | **research/pde_ledger_v3/steps/S11bB_interface_assembly.md:74–76**: “The transverse mode is completely decoupled on a uniform background”; coupling “identically zero”, dispersion “`ρ_br⁰ω² = μ_R k²`”, imaginary part “zero”; “Both engines” find that “in-plane parity admits no `e_W ↔ u_T` bilinear”. [Q5](#q5) | **Two-engine uniform-background result as stated by the source**, within the assembled interface's field content and assumptions. The stated parity exclusion concerns this thickness/transverse bilinear; it does not select all parity/chiral content outside that object. |
| Current confinement/conversion closeout | **research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md:3**: “It closes PARTIAL: no accepted nonuniform light-loss number was obtained.” Line 27: nonuniform confinement and physical admissibility of rotational stiffness “are not established”. [R9](#r9) | **PARTIAL**, with uniform decoupling conditional. This is the current closeout, not proof of universal photon stability. |
| Stiffness provenance | **research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:226–244**, R-S8-01: “status **OPEN**”; stiffness functional must be “delivered by the substructure rather than chosen”; curl form “postulated (structural)”. [R11](#r11) | **OPEN substrate requirement**. Matching two modes does not supply the material law. |
| Dimension provenance | **research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:128–145**, R-S1-01: Dbrane as a derived quantity or “explicitly re-affirmed postulate”, status **OPEN**; “Dbrane = 3 went in and D−1 = 2 came out”. [R28](#r28) | **OPEN**, target S6. The record does not prohibit deriving a codimension-one brane from a four-dimensional bulk elsewhere. |
| Normal h versus tangential u | **research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:263–278**, R-S8-02: “status **OPEN**”; “full displacement, in-plane u **and** out-of-plane h”; “h ≠ uL” is “confirmed as a picture, not computed”. [R28](#r28) | **OPEN joint-operator requirement**, target S8. The source distinguishes a scalar normal displacement from the in-plane longitudinal scalar. |
| The specific possible third transverse-sector mode | **research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:270–275**: h/uL scalar mixing can leave the transverse count D−1; “The condition that gives 3 is h being **DEGENERATE WITH THE TRANSVERSE PAIR**”, elasticity isotropic in D+1 rather than D. Lines 300–305: S11 excludes hBranon and does not state the census with it included. [R28](#r28) | **Open conditional failure condition**, not a proved third photon and not a proof that normal motion is forbidden. The inventory records the distinction without settling it. |
| Thickness inertia is not the out-of-plane field | **research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:412–425**, **R-S8-06**, asks for material-displacement identity and the quadratic inertia of material and thickness. Its note at **:424–425** says it asks for “field identity and inertia”, not numerical benchmark values, and “does not identify the thickness mode with S10's out-of-plane displacement”. [Q5](#q5) | **OPEN**, target S8; the source labels that owner a **register inference** for the original S11/S11b records (O2 names S8 for material identification). This is an unmet identity/inertia requirement, not an identification of breathing motion with the excluded normal field. |
| Positive stiffness | **research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:307–316**, R-S8-03: “status **OPEN**”; without the appropriate sign “two exponentially growing modes rather than two waves”. [R11](#r11) | **OPEN substrate sign requirement**. Counting roots does not establish stable light. |
| Angular-momentum carrier | **research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:334–345**, R-S8-04: “status **OPEN**”; asks what carries “internal angular momentum” or couple-stress; without it, curl-only stiffness is “not an admissible continuum mechanics”. [R18](#r18) | **OPEN material admissibility requirement**. This is not a derived spin-\(\hbar\) value for a photon. |
| Ruled-out polar-P route | **research/pde_ledger_v3/DEFECT_REGISTER.md:156**, B2: “route-failure for deriving the shear modulus \(\mu_R\) from a polar substructure P”; “not a finding that one medium cannot carry a longitudinal and a transverse mode”; status **FALSIFIED**, with stated imposed-axis/MacCullagh conditions. [R16](#r16) | **FALSIFIED narrow route**, not the entire shear-light ontology. P retirement and surviving bare postulated stiffness are the register's stated scope. |
| A different “spin problem” | **research/pde_ledger_v3/DEFECT_REGISTER.md:177**, C10: “a structural tension INSIDE the gravity sector”; “Inertia wants compact; spin wants extended”; status **OPEN**. [R16](#r16) | **OPEN throat/gravitomagnetic geometry problem**, not an established calculation of photon helicity. |
| Gravity record's optical limit | **research/pde_ledger_v3/steps/O2_steady_brane_balance.md:125–127**: “no … direction/polarization-dependent optical stiffness … or subprincipal/polarization transport result”. [R17](#r17) | **Explicitly outside this record's delivered results**. Its steady balance does not establish helicity-independent lensing, differential delay or gravitational spin Hall transport. |
| Birefringence named for later work | **research/pde_ledger_v3/V3_STEP_PLAN.md:1174–1184**: “proposed linear tests on supplied backgrounds”; item 4 is “birefringence near a defect”, with homogeneous degeneracy and defect splitting. **research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md:31**: S22 items 1/4 are **OPEN**, physical magnitude needs interior input. [R19](#r19), [R9](#r9) | **Proposed/OPEN**, not measured or a completed field-dependent optical law. The plan's future wording does not override the closeout. |
| Polarization phase is not conversion loss | **research/pde_ledger_v3/steps/LIGHT_LEAKAGE_SCOPING.md:36**: “PRESENT ANALOG, named future item”; “OPEN”; “no completed birefringence law”; magnetic-birefringence comparison has a “defect-observable map gap”. Line 37 quotes older uniform “unconditional” wording alongside current **CONDITIONAL** wording. [R27](#r27) | **OPEN comparison/observable map**, with an explicit historical wording tension. Polarization redistribution among surviving light modes is not automatically lost power. |
| Quantum photons are a separate scope | **research/pde_ledger_v3/V3_STEP_PLAN.md:1198–1204** distinguishes a discrete bound spectrum, a nonlinear soliton, and \(\hbar\omega\) quanta: the last needs “field quantization”; “not a classical PDE sim at all”; quantum mechanics is “a separate project”. [R19](#r19) | **Outside the classical simulation's present deliverable**, not a quantization, Born-rule or entanglement result. |
| Bounded search for polarization/measurement language | The exact keyword search in **R31** returned **empty stdout, exit 1**, over the plan, defect register, substrate register and top-level step Markdown. [R31](#r31) | **Search result only**: no matches for the listed Jones/Poincaré/helicity/handedness/entanglement/Malus/Faraday/Kerr/Bell/spin-Hall phrases in that scope. It does not prove absence in every repository file or rule out unnamed future constructions. |
| Plan's normal displacement and charge-field identity | **research/pde_ledger_v3/V3_STEP_PLAN.md:453–462**: “`u_L` is an **in-plane longitudinal brane displacement**”; “the charge the rethink produced is carried by the **`h`-branon**”; the held mouth datum is “`h_A = ξ_w\|_A/ℓ = P₀H\|_A` … **(distinct from `u_L`)**”; h “remains the committed mediator for the conditional `1/R²` falloff”. Lines **897–898** state “the field identity `ξ_w = ℓh`” and an orientation-odd mouth source. [P6](#p6) | **Plan-level committed identification**, explicitly h ≠ uL. The normal w displacement is the plan's charge field, not a forbidden displacement or an identification of uL with charge. This is not a completed current joint u/h spectrum; that requirement remains OPEN as quoted above. |
| Plan's static electric-scalar closure | **research/pde_ledger_v3/V3_STEP_PLAN.md:871–893**, Q1: the static scalar uses “localized-`H` / PT”; \(M_h=N_0M_4,\ K_h=N_0K_4=M_hc_E^2\) is “**EARNED — given the postulated action**”; “derived, conditional on the action”; regime “static”; “the action is postulated — a tier-1 item”. C6 is “Explicitly deferred, not closed”. [P6](#p6) | **Plan-level conditional derived reduction from a postulated action**, not from primitives. cE is an interior quantity; Q1 does not resolve the parent-action gap. |
| Plan's static falloff and source | **research/pde_ledger_v3/V3_STEP_PLAN.md:896–905**, Q2: “the holder is a DEBT, not a result”; a ±w puncture bends the brane into ±w; expected output is “source identity” and “far-field FORM”; “the **`1/R²` falloff** and the **`s₁s₂` product** are **target-blind EARNED**”, but “**within Q1's postulated G0 closure** … not from primitives”. [P6](#p6) | **Plan-level EARNED form within the postulated Q1 closure**; holder remains a DEBT. This is the conditional static charge-mediator source relevant to E02, not a Coulomb-law error or photon-mass bound computed here. |
| S9b's supplied polarization identification | **research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md at commit `ede8aa21`:53–55**: “There is one isotropic speed, the same for both polarizations”; “the supplied polarization identification is `c_γ,1(x) ≡ c_γ,2(x) ≡ c_γ(x)`”. [P4](#p4) | **Supplied input in the cited spec version**, not a derived polarization-independence result. The commit version is used because the working copy is under amendment, as specified in the repair request. |
| What that S9b spec excludes | **research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md at commit `ede8aa21`:106–109**: “Outside this step”; “recorded as open”: “direction-dependent (radial versus tangential) stiffness”; “coupling to thickness or bulk fields”; “polarization-dependent propagation”. [P4](#p4) | **Outside the step; OPEN in that spec version**. Thus its supplied common speed cannot establish those excluded effects or reproduce polarization independence from independent inputs. |
| S11c-d's polarization-dependent forcing | **research/pde_ledger_v3/steps/S11c_d_profile_conditioned_scattering.md:58** labels the finite packet/source/matching evidence **“CONDITIONAL”**, including “polarization-dependent first-order forcing” and “matched transverse amplitudes”; “None supplies complete physical work, a full field or leakage.” [P3](#p3) | **CONDITIONAL profile-conditioned evidence**, with weighted-class restrictions in the same line. It is not a completed polarization propagation, physical-work, full-field or leakage result. |
| Proper rotations, reflections and the stiffness count | **research/pde_ledger_v3/steps/S11_stray_longitudinal.md:40–41** defines a quadratic stiffness in \(\partial_i u_j\) [R12](#r12). **S11_stray_longitudinal.md:52–64** gives “`N_SO = {D=2 → 4, D=3 → 3, D=4 → 4, D=5 → 3}`” and “`N_O = 3 for every D`”. **At D = 3, proper rotations and full orthogonal symmetry therefore allow the same three invariants: there is no additional reflection-odd invariant in this quadratic first-gradient in-plane stiffness class.** The reflection-odd D=2 extra is not a total derivative and can mix longitudinal/transverse sectors; the D=4 extra is a total derivative and adds no operator structure. “Sector separation is a `D=3` fact, not a structural one.” [Q5](#q5) | **Corrected derived invariant count within that stiffness class**, agreed by five derivations as the source states. The open parity/chirality condition is **outside this class**: O2 does not select the wider constitutive content, S11bB sets a reservoir/power-budget condition for non-passive couplings, and S9b excludes subprincipal polarization transport, as quoted in the rows below [Q4](#q4). No general parity/chirality selection is resolved here. |
| O2's normal geometry and optical counting | **research/pde_ledger_v3/steps/O2_steady_brane_balance.md:89–92** calls \(g_{ij}=\delta_{ij}+\partial_i\xi_w\partial_j\xi_w\), \(\xi_w=\ell h\), \(c_\gamma^2=\mu_\perp/\rho_{\rm br}\) and the mass law “supplied geometric/optical/mass inputs”. **O2_steady_brane_balance.md:115–117** supplies \((\partial\xi_w)^2=O(\epsilon)\) and the optical monomial box; lines 118–120 leave other grades OPEN and say O2 is untruncated. [Q4](#q4) | **Supplied inputs/counting**, explicitly not earned results (lines 95–100). This uses normal displacement as geometry affecting the optical object as well as retaining the plan's h identity. |
| O2 does not select parity or absence of chiral content | **research/pde_ledger_v3/steps/O2_steady_brane_balance.md:110–111**: “Radial profiles impose no constitutive isotropy, parity, stress symmetry, absent couple/chiral content, derivative cutoff or finite history state.” [Q4](#q4) | **Explicitly unselected constitutive content** in a conditional model, not a derived parity-even/nonchiral substrate. Radiality is not declared to close those choices. |
| Odd/chiral interface couplings are not generally prohibited | **research/pde_ledger_v3/steps/S11bB_interface_assembly.md:53–54** says nonreciprocal “odd” couplings occur in driven laboratory media, including active/chiral fluids, and “They are not unphysical”. Lines 55–63 give the standing rule: “A non-passive coupling is admissible only with a NAMED reservoir and a STATED power budget.” [Q4](#q4) | **Correction and conditional admissibility rule**, not a blanket thermodynamic exclusion or a computed birefringence coefficient for this brane. These are the record's statements about admissibility, not adoption of an odd/chiral optical law. |
| S9b's eikonal versus polarization-transport scope | **research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md at commit `ede8aa21`:117–119**: “The retained object is the dispersion relation above … and its rays”; “Excluded … and not computed: the explicit subprincipal terms … which affect amplitude and polarization transport.” Lines **34–44** label the local optical dispersion/metric **supplied**; lines **165–186** name flyby deflection, round-trip excess time/logarithmic slope and one-way excess times as observables to compute. [Q4](#q4) | **Supplied eikonal object and specified computands in that spec version; subprincipal transport excluded/not computed.** This supplies an optical ray object for a comparison, not an earned helicity-transport law, parity selection or a numerical comparison with E22–E24. |

### Historical papers and exploratory documents: leads, not authority

The current S11c closeout calls the linked docs handoff material **“exploratory inputs for S22/Q2 … not adopted laws or solved throats”**: **research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md:35–39** [R9](#r9). The following sources therefore retain their own labels without gaining v3 authority.

| Subject | Source, literal quotation and evidence | Source's status and relation to the current ledger |
|---|---|---|
| Old exact-Maxwell and automatic-confinement prose | **docs/conceptual_foundation.md:372–394** says curl stiffness gives “exactly two transverse polarizations and no longitudinal mode” and “real electromagnetism”. Lines 474–478 label a “Working synthesis” in which light has “2 transverse polarizations automatically” and “never leaves the brane”. [R35](#r35) | **Working synthesis; stale relative to the narrower v3 statements** in S9/S10 and S11c. The opposing current quotations are retained above; no reconciliation is asserted. |
| Historical normal-mode exclusion | **docs/conceptual_history.md:255–262** says the on-brane scalar uw “must stay gapped”; a phase-sector revival is “not a claim”. [R35](#r35) | **Historical hypothesis/exclusion statement**. It is not a current joint u/h spectrum. The old v2 masslessness example and current R-S8-02 requirement below are distinct recorded statements. |
| Known-stale model map | **docs/model_map.md:16**: “This is a synthesis map, not the source of truth.”. Lines 111–113 say two polarizations establish three-dimensional space; lines 122–125 call cosmic photon stability a “CONSEQUENCE of bulk shear-freeness”. [R23](#r23), [R35](#r35) | **Navigational synthesis, known stale per the task**. S10 does not select D, and S11c has no accepted nonuniform loss number. These sentences are not used as classification authority. |
| v2 transverse count | **research/pde_ledger_v2/paper/stages/stage_003.tex:37–45**: “Transverse (earned)” and “two polarizations … physical_dof=2, massless”. [R24](#r24) | **“Earned” in the old supplied-action calculation**, an older lead. Current S10 supplies the stronger qualifications on field content, inertia and dimensionality. |
| v2 normal/longitudinal calculation exists | **research/pde_ledger_v2/paper/stages/stage_030.tex:99–128** labels the coupled (uL,h) scalar block **“EARNED”**, with stiffness \([[B_{\rm eff},C_{hu}],[C_{hu},K_h]]\); at normalized sample values \(z_\pm=(3\pm\sqrt2)/2>0\); reduced-h masslessness also **“EARNED”**. [R24](#r24) | **Old conditional calculation**, not no computation anywhere. Current **SUBSTRATE_REQUIREMENTS.md:281–305** explicitly cites this lead while saying v3 has not delivered its full u/h object [R28](#r28). A scalar-block example is not proof of an additional transverse photon. |
| Guided-light plan's geometry | **docs/light_guided_photon_soliton_research_plan.md:7** says it is a plan, not obtained results. Lines 246–274 impose \(N\cdot u_T=0\), \(k\cdot u_T=0\) and count \(4-1-1=2\), then say it “must be verified by the complete guided-mode spectrum”. [R22](#r22) | **Conditional target/proposed research**, not a demonstrated exclusion of normal modes. |
| Proposed material isotropy and nonlinear polarization criterion | **docs/light_guided_photon_soliton_research_plan.md:278–285** proposes no preferred x/y/z direction but a different w response. Lines 1788–1791 require only two independent light-like polarizations and no unacceptable linear/circular speed or energy split. [R22](#r22) | **Design/acceptance criteria**, not derived material laws or an experimental prediction already tested. |
| Circular shear as candidate helicity | **docs/native_light_em_and_vortex_throat_interpretation.md:2278–2282**: circular polarization “may … furnish two wave-helicity states”, but angular-momentum density and detector response “must be calculated”. [R21](#r21) | **Conditional candidate and explicit missing calculations**. Having two transverse displacement directions does not establish the Beth torque or a photon helicity representation. |
| Fluid helicity versus particle helicity | **docs/native_light_em_and_vortex_throat_interpretation.md:2284–2299** defines a projected \(\int v_{\rm br}\cdot(\nabla\times v_{\rm br})\,d^3x\) and says it “is not the same as the helicity of a particle representation” or a complete four-dimensional invariant. [R21](#r21) | **Terminological distinction**, not a derived photon-helicity identification. |
| Explicit photon spin and radiation status | **docs/native_light_em_and_vortex_throat_interpretation.md:4214**: acceleration radiation **“open”**, “two polarization amplitudes and extra emission channels are not derived”. Line 4218: “spin-1 helicity of the light quantum” is “outside the present classical derivation”; circular shear is a “candidate classical precursor only”. [R21](#r21) | **OPEN radiation source; quantum spin outside scope**, as the exploratory source states. |
| Handedness used for a different throat quantity | **docs/em_sector_reconsideration.md:3,9,40**: “WORKING DRAFT”, “not a settled result”; magnitude/helicity(±w) split is approximate with unknown mixing angle; in four dimensions “there is no unique 4D ‘right-hand rule’/scalar helicity”. [R35](#r35) | **Working hypothesis with a corrected earlier error**. Its ±w throat-orientation language is not a computation of the two photon helicities. |
| More careful ontology synthesis | **docs/toy_model_ontology_summary.md:947–961** calls the count conditional on the action, isotropic inertia and structural premises; confinement remains an interface question, with nonuniform coupling incomplete. [R36](#r36) | **Conditional synthesis with explicit OPEN items**, consistent in scope with the cited v3 count/closeout; still not an additional authority. |
| Extra electric radiation and polarization burden | **docs/cross_sector_research_burdens_and_compatibility_gates.md:211–230** distinguishes an unestablished first-class Gauss constraint from a possible extra dynamical carrier; lines 265–274 list a third freely radiating electromagnetic polarization among failure criteria; lines 973–979 include relative polarization speed/energy splitting in a later acceptance budget. [R36](#r36) | **Reviewer inference and acceptance criteria**, not a demonstrated third photon or a resolved coupling law. |
| Older EM paper's additional modes | **research/4d_em_fields/paper/4d_em_fields.tex:1415–1426**: above a massive-mode threshold, “additional propagating channels open”; “in principle, additional polarization content” controlled by the supplied tower and couplings. [R36](#r36) | **Prediction within that older controlled reduction**, not adopted as the current shear-sector spectrum. |
| Older gravity paper's polarization mention | **research/1pn_hybrid/paper/1pn_hybrid.tex:1798–1812** lists “potential signatures in high-energy emission or polarization patterns” among observational targets. [R36](#r36) | **Potential future observational handle**. These lines provide no computed helicity-differential deflection or delay. |
| Plasma circular wave and magnetic helicity | **research/4d_plasma/paper/4d_plasma.tex:2575–2590** lists circularly polarized Alfvén waves as a code benchmark. Lines 4291–4323 discuss \(A\cdot B\) helicity on a fixed-w 3+1 projection. [R26](#r26) | **Plasma benchmark proposal and projected-field identity**, not a v3 photon-polarization result. |
| Proposal chapter's unsigned normal axis and nonlinear targets | **research/research_proposal/paper/proposal_chapter_04_light_research_program.tex:40–49** says free light should use an unsigned normal axis at leading order and the complete guided spectrum must realize the two-direction count. Lines 329–337 call for deriving and comparing nonlinear linear/circular/elliptical responses. [R34](#r34) | **Research program/requirements**, not solved polarization energetics, handedness transport or a field-dependent speed law. |

## 4. Classification against the model

Apply these definitions **in this order; a row takes the first class whose definition it meets, except as stated in rule 4**. **(1) Reproduced.** A v3 step computes the measured fact from inputs that do not themselves state it. **(2) In apparent conflict.** A model source states something opposed to the measured fact, or the ledger itself records a departure from it. An open condition, such as an unresolved coupling, makes the conflict conditional; state the condition. Against a reported signal that is not established, the conflict is conditional on the signal too. **(3) Testable.** A model source supplies or computes the object the measurement constrains (for example, the spectrum whose splitting a birefringence bound limits), but no step has compared the two. **(4) Required, no mechanism yet.** The fact concerns light within the model's current scope (propagation in the brane, and interaction with the brane's own structures and defects), and no model source supplies the object it constrains. **A fact that follows only from a supplied input that states it belongs here, and this takes precedence over rule 3**; name that input. **(5) Not addressed.** The fact depends on physics the model has not defined, such as the response of ordinary matter, material interfaces, or emission processes, and **“no model source bears on it”**. State each row's conditions in its cell. A future item or missing calculation alone is not a supplied comparison object, and missing mechanisms do not demote an earlier ledger-recorded departure.

The current quantum boundary is explicit: **V3_STEP_PLAN.md:1198–1204** puts field quantization in a separate project [P8](#p8). Thus the quantum preparation/outcome entries below concern undefined emission/measurement physics, rather than a delivered classical propagation object. The classifications do not propose mechanisms or settle any open condition. E10–E15 retain both the homogeneous-degeneracy statement and the parity/chiral/transport qualifications in part 3; the [β note](#beta-tests) gives the papers' distinct physical target.

| ID | Classification | Conditions, object and source; reason the ordered rule selects this class |
|---|---|---|
| E01 | **Reproduced** | **Conditions:** classify the two-dimensional transverse branch content, under the selected homogeneous, isotropic-inertia in-plane \(D=3\) action, positive coefficients, nonzero wavevector and the stated generic strata; separation from out-of-plane fields is supplied. **Rule 1:** **S10_two_transverse_photons.md:73–76** computes “the nonzero root has D − 1 transverse null directions” and hence two at D=3 [Q5](#q5). Its action, field content and dimensional input do **not themselves state the two-direction result**; supplying D=3 is not supplying D−1=2. Thus rule 4's input exception does not apply to this computed count. This is a **conditional reproduction of that count**, not a derivation of D=3, the quantum tomography protocol or a complete radiation-state census. **S10:172–197** explicitly leaves physical dimension selection and out-of-plane separation unearned [Q5](#q5); R-S1-01 remains **OPEN**, **SUBSTRATE_REQUIREMENTS.md:128–145** [P5](#p5), and the joint u/h census remains **OPEN**, **:263–305** [R28](#r28). The separate longitudinal departure is classified in E32. |
| E02 | **Testable** | **Conditions:** comparison with this laboratory null result must specify the static mediator and charge/matter response; a Proca interpretation additionally needs an observable map to that alternative theory. **V3_STEP_PLAN.md:453–462,871–905**, Q1/Q2, identifies \(\xi_w=\ell h\), the h-branon charge scalar and mediator, with static \(1/R^2\) form “EARNED **within Q1's postulated G0 closure**”, not from primitives; its holder remains a DEBT and the parent action deferred. [P6](#p6) This conditional plan-level form is the relevant model statement, not a computed Coulomb-law error bar or photon-mass limit. S10's gapless transverse root, **S10:237–252** [R10](#r10), does not by itself supply that measurement. The native Gauss/extra-carrier status remains open in the exploratory sources, **docs/native_light_em_and_vortex_throat_interpretation.md:4203–4204** [R21](#r21), **docs/cross_sector_research_burdens_and_compatibility_gates.md:211–230** [R36](#r36). **Rule 3 object:** the conditional static charge-mediator falloff is supplied/computed at plan level, but no step compares its departures with this measured Coulomb-law bound. |
| E03 | **Testable** | **Conditions:** interpret the experiment's Proca photon-mass parameter as a gap of the native radiating branch, with the cited matter/electromagnetic identification; this identification is not established by counting modes. **Rule 3 object:** S10 computes a gapless transverse spectrum under supplied homogeneous inputs, **S10_two_transverse_photons.md:237–252** [R10](#r10), but no step compares it with the quoted solar-wind bound. This is testable only under that identification, not a computed Proca current law or MHD result; the matter/interface observable map remains unresolved, **S11_stray_longitudinal.md:154–163** [R12](#r12). |
| E04 | **Required, no mechanism yet** | **Conditions:** circular classical light transfers the measured torque; the quantum per-photon interpretation additionally requires quantization and receiver identification. Circular-wave angular-momentum density and detector/torque response remain to be calculated, **docs/native_light_em_and_vortex_throat_interpretation.md:2278–2282** [R21](#r21); the substrate's internal angular-momentum carrier is **OPEN**, **SUBSTRATE_REQUIREMENTS.md:334–345** [R18](#r18). Neither is a Beth torque result. **Rule 4 scope:** the intrinsic classical angular momentum carried by propagating light is within the light-wave scope; the source supplies no angular-momentum/flux object to compare. This does not claim that the laboratory plate's ordinary-matter response has been defined. |
| E05 | **Required, no mechanism yet** | **Conditions:** the structured beam and absorbing receiver must be represented together to compare the measured mechanical OAM transfer. Mechanical OAM transfer has no demonstrated native light angular-momentum/receiver law in these records. The open angular-momentum and radiation-response statements are **SUBSTRATE_REQUIREMENTS.md:334–345** [R18](#r18) and **docs/native_light_em_and_vortex_throat_interpretation.md:4214,4218** [R21](#r21). A spatially structured shear wave is not, by that description alone, a reproduced particle torque. **Rule 4 scope:** the beam's classical orbital angular momentum is a propagation property; no native orbital-angular-momentum/flux object is supplied. The experiment's absorbing-particle response also remains undefined. |
| E06 | **Not addressed** | **Conditions:** classify the measured Born probabilities and heralded preparation/detector statistics, not just classical wave amplitudes or a spatial mode count. **Rule 5:** that quantum emission/measurement physics is not defined in the present classical brane model; **V3_STEP_PLAN.md:1198–1204** explicitly assigns field quantization to a separate project [P8](#p8). No source here supplies the quantum state/outcome object this measurement constrains. A missing future quantum completion is not a requirement inside the current propagation/defect scope. |
| E07 | **Not addressed** | **Conditions:** classify the measured quantized photon production and indivisible detector outcomes, not just classical wave amplitudes or a spatial mode count. **Rule 5:** that quantum emission/measurement physics is not defined in the present classical brane model; **V3_STEP_PLAN.md:1198–1204** explicitly assigns field quantization to a separate project [P8](#p8). No source here supplies the quantum state/outcome object this measurement constrains. A missing future quantum completion is not a requirement inside the current propagation/defect scope. |
| E08 | **Not addressed** | **Conditions:** classify the measured quantum pair preparation and joint measurement statistics, not just classical wave amplitudes or a spatial mode count. **Rule 5:** that quantum emission/measurement physics is not defined in the present classical brane model; **V3_STEP_PLAN.md:1198–1204** explicitly assigns field quantization to a separate project [P8](#p8). No source here supplies the quantum state/outcome object this measurement constrains. A missing future quantum completion is not a requirement inside the current propagation/defect scope. |
| E09 | **Not addressed** | **Conditions:** classify the measured quantum pair preparation and the corrected joint detector statistics, not just classical wave amplitudes or a spatial mode count. **Rule 5:** that quantum emission/measurement physics is not defined in the present classical brane model; **V3_STEP_PLAN.md:1198–1204** explicitly assigns field quantization to a separate project [P8](#p8). No source here supplies the quantum state/outcome object this measurement constrains. A missing future quantum completion is not a requirement inside the current propagation/defect scope. |
| E10 | **Testable** | **Conditions:** the selected homogeneous, isotropic-inertia in-plane \(D=3\) linear pair applies, and **parity/chiral constitutive content and polarization-transport terms outside that retained object make no additional splitting/rotation**. **Rule 3 object:** the transverse spectrum/common root is computed conditionally in **S10:73–86** [R10](#r10); **V3_STEP_PLAN.md:1183** states “the two polarisations are degenerate *by symmetry* in a homogeneous brane” [P3](#p3). No step has compared that spectrum with this dimension-five GRB bound or computed its path-integrated rotation. **S11_stray_longitudinal.md:52–64** finds **N_SO=N_O=3 at D=3**, with **no additional reflection-odd invariant in its quadratic first-gradient in-plane stiffness class**; the D=2 extra can mix longitudinal/transverse sectors, while the D=4 extra is a total derivative and adds no operator structure [Q5](#q5). The unresolved parity/chirality condition lies **outside that class**: **O2_steady_brane_balance.md:110–111** leaves wider constitutive parity/chiral content unselected; **S11bB_interface_assembly.md:53–63** permits non-passive odd couplings only with a named reservoir/power budget; **S9b_SHARED_PHYSICS.md at `ede8aa21`:117–119** excludes subprincipal polarization transport [Q4](#q4). The [β papers' parity-violation/rotation target](#beta-tests) in E13–E15 is distinct from this specified dispersion bound. These open conditions are not resolved here. |
| E11 | **Testable** | **Conditions:** the selected homogeneous, isotropic-inertia in-plane \(D=3\) linear pair applies, and **parity/chiral constitutive content and polarization-transport terms outside that retained object make no additional splitting/rotation**. **Rule 3 object:** the transverse spectrum/common root is computed conditionally in **S10:73–86** [R10](#r10); **V3_STEP_PLAN.md:1183** states “the two polarisations are degenerate *by symmetry* in a homogeneous brane” [P3](#p3). No step has mapped that spectrum to the constrained dimension-four SME coefficients or compared the measured bound. **S11_stray_longitudinal.md:52–64** finds **N_SO=N_O=3 at D=3**, with **no additional reflection-odd invariant in its quadratic first-gradient in-plane stiffness class**; the D=2 extra can mix longitudinal/transverse sectors, while the D=4 extra is a total derivative and adds no operator structure [Q5](#q5). The unresolved parity/chirality condition lies **outside that class**: **O2_steady_brane_balance.md:110–111** leaves wider constitutive parity/chiral content unselected; **S11bB_interface_assembly.md:53–63** permits non-passive odd couplings only with a named reservoir/power budget; **S9b_SHARED_PHYSICS.md at `ede8aa21`:117–119** excludes subprincipal polarization transport [Q4](#q4). The [β papers' parity-violation/rotation target](#beta-tests) in E13–E15 is distinct from this specified dispersion bound. These open conditions are not resolved here. |
| E12 | **Testable** | **Conditions:** the selected homogeneous, isotropic-inertia in-plane \(D=3\) linear pair applies along the tested paths, and parity/chiral constitutive content and polarization-transport terms outside the retained stiffness class make no additional splitting/rotation. The comparison must use the tested **anisotropic** scale-invariant \(C_L^{\alpha\alpha}\), with a CMB polarization source and sky statistics to map propagation to that estimator; E13–E15 need the same source/statistical mapping for their uniform-rotation estimators. **Rule 3 object, the same as E10–E11 and the model side of E13–E15:** **S10:73–86** computes the common transverse root [Q5](#q5), and **V3_STEP_PLAN.md:1183** states “the two polarisations are degenerate *by symmetry* in a homogeneous brane” [P3](#p3). No step has compared that object with this anisotropic bound or computed its CMB rotation-power spectrum. **S11:52–64** finds **N_SO=N_O=3 at D=3**, with **no additional reflection-odd invariant in its quadratic first-gradient in-plane stiffness class**; the D=2 extra can mix sectors, while the D=4 extra is a total derivative and adds no operator structure [Q5](#q5). The unresolved condition lies **outside that class**: **O2:110–111** leaves wider parity/chiral content unselected, **S11bB:53–63** admits non-passive odd couplings only with a named reservoir/power budget, and **S9b at `ede8aa21`:117–119** excludes subprincipal polarization transport [Q4](#q4). Missing CMB source/statistical mapping does not select rule 5, because these model sources already bear on propagation rotation. No open condition is resolved. |
| E13 | **In apparent conflict** | **Conditions:** the selected homogeneous, isotropic-inertia in-plane \(D=3\) linear pair applies, and **parity/chiral constitutive content and polarization-transport terms outside that retained object make no additional splitting/rotation**. The reported signal must additionally be **real cosmic propagation rotation**, not calibration/foreground/systematic bias. Comparing its uniform CMB estimator also requires a cosmological polarization source and sky statistics, the same source/statistical mapping condition needed for E12's anisotropic estimator. **Experimental side:** “\(\beta=0.35\pm0.14^\circ\)”, [primary paper](https://arxiv.org/pdf/2011.11254), is a tentative nonzero β; these authors frame β as a test of **parity-violating cosmic rotation**, as detailed in the [source-specific β note](#beta-tests). **Model side:** “the two polarisations are degenerate *by symmetry* in a homogeneous brane”, **V3_STEP_PLAN.md:1183** [P3](#p3), with S10's common root conditional on supplied inputs, **S10:73–86** [R10](#r10). **Rule 2:** this is an apparent conflict only under the stated signal and parity/transport conditions. **S11_stray_longitudinal.md:52–64** finds **N_SO=N_O=3 at D=3**, with **no additional reflection-odd invariant in its quadratic first-gradient in-plane stiffness class**; the D=2 extra can mix longitudinal/transverse sectors, while the D=4 extra is a total derivative and adds no operator structure [Q5](#q5). The unresolved parity/chirality condition lies **outside that class**: **O2_steady_brane_balance.md:110–111** leaves wider constitutive parity/chiral content unselected; **S11bB_interface_assembly.md:53–63** permits non-passive odd couplings only with a named reservoir/power budget; **S9b_SHARED_PHYSICS.md at `ede8aa21`:117–119** excludes subprincipal polarization transport [Q4](#q4). No general parity/chirality selection or cosmic β transport law is settled by those records; no resolution is asserted. |
| E14 | **In apparent conflict** | **Conditions:** the selected homogeneous, isotropic-inertia in-plane \(D=3\) linear pair applies, and **parity/chiral constitutive content and polarization-transport terms outside that retained object make no additional splitting/rotation**. The reported signal must additionally be **real cosmic propagation rotation**, not calibration/foreground/systematic bias. Comparing its uniform CMB estimator also requires a cosmological polarization source and sky statistics, the same source/statistical mapping condition needed for E12's anisotropic estimator. **Experimental side:** “\(\beta=0.30\pm0.11^\circ\)”, [primary paper](https://arxiv.org/pdf/2203.04830), is a tentative nonzero β; these authors frame β as a test of **parity-violating cosmic rotation**, as detailed in the [source-specific β note](#beta-tests). **Model side:** “the two polarisations are degenerate *by symmetry* in a homogeneous brane”, **V3_STEP_PLAN.md:1183** [P3](#p3), with S10's common root conditional on supplied inputs, **S10:73–86** [R10](#r10). **Rule 2:** this is an apparent conflict only under the stated signal and parity/transport conditions. **S11_stray_longitudinal.md:52–64** finds **N_SO=N_O=3 at D=3**, with **no additional reflection-odd invariant in its quadratic first-gradient in-plane stiffness class**; the D=2 extra can mix longitudinal/transverse sectors, while the D=4 extra is a total derivative and adds no operator structure [Q5](#q5). The unresolved parity/chirality condition lies **outside that class**: **O2_steady_brane_balance.md:110–111** leaves wider constitutive parity/chiral content unselected; **S11bB_interface_assembly.md:53–63** permits non-passive odd couplings only with a named reservoir/power budget; **S9b_SHARED_PHYSICS.md at `ede8aa21`:117–119** excludes subprincipal polarization transport [Q4](#q4). No general parity/chirality selection or cosmic β transport law is settled by those records; no resolution is asserted. |
| E15 | **In apparent conflict** | **Conditions:** the selected homogeneous, isotropic-inertia in-plane \(D=3\) linear pair applies, and **parity/chiral constitutive content and polarization-transport terms outside that retained object make no additional splitting/rotation**. The reported signal must additionally be **real cosmic propagation rotation**, not calibration/foreground/systematic bias. Comparing its uniform CMB estimator also requires a cosmological polarization source and sky statistics, the same source/statistical mapping condition needed for E12's anisotropic estimator. **Experimental side:** “\(\beta=0.215\pm0.074^\circ\)”, [primary paper](https://arxiv.org/html/2509.13654v2), is a tentative nonzero β; these authors frame β as a test of **parity-violating cosmic rotation**, as detailed in the [source-specific β note](#beta-tests). **Model side:** “the two polarisations are degenerate *by symmetry* in a homogeneous brane”, **V3_STEP_PLAN.md:1183** [P3](#p3), with S10's common root conditional on supplied inputs, **S10:73–86** [R10](#r10). **Rule 2:** this is an apparent conflict only under the stated signal and parity/transport conditions. **S11_stray_longitudinal.md:52–64** finds **N_SO=N_O=3 at D=3**, with **no additional reflection-odd invariant in its quadratic first-gradient in-plane stiffness class**; the D=2 extra can mix longitudinal/transverse sectors, while the D=4 extra is a total derivative and adds no operator structure [Q5](#q5). The unresolved parity/chirality condition lies **outside that class**: **O2_steady_brane_balance.md:110–111** leaves wider constitutive parity/chiral content unselected; **S11bB_interface_assembly.md:53–63** permits non-passive odd couplings only with a named reservoir/power budget; **S9b_SHARED_PHYSICS.md at `ede8aa21`:117–119** excludes subprincipal polarization transport [Q4](#q4). No general parity/chirality selection or cosmic β transport law is settled by those records; no resolution is asserted. |
| E16 | **Not addressed** | **Conditions:** magnetized material transmission; a vacuum-degeneracy statement does not supply its material response. No magnetic-field-dependent material polarization-response law is supplied by the audited v3 optical results; bounded Faraday search [R31](#r31), and “no completed birefringence law” in **LIGHT_LEAKAGE_SCOPING.md:36** [R27](#r27). A two-mode homogeneous dispersion is not Faraday rotation in matter. **Rule 5:** Faraday rotation in ordinary magnetized matter depends on its dielectric/magnetic response, not the supplied brane shear spectrum. The cited open analogs or future labels supply no object for this specific measured response. |
| E17 | **Not addressed** | **Conditions:** electric-field-induced dielectric response in matter, under the measured material conditions. No derived electric-field-dependent dielectric polarization response in the current audited records; bounded Kerr search [R31](#r31), open optical law **LIGHT_LEAKAGE_SCOPING.md:36** [R27](#r27). Its relation to a completed native material/electric sector is unclear. **Rule 5:** electric-field Kerr birefringence in ordinary dielectric matter depends on its undefined material response. The cited open analogs or future labels supply no object for this specific measured response. |
| E18 | **Required, no mechanism yet** | **Conditions:** identify the applied magnetic background with the brane's own electromagnetic structures/defects so that this is a constraint on native light propagation, not a material dielectric response. **Rule 4:** no source supplies the needed \(B\)-dependent polarization spectrum/phase response; a named future item is not that object. **LIGHT_LEAKAGE_SCOPING.md:36,152** records a present analog and an OPEN defect-observable map, with “no completed birefringence law” [R27](#r27); physical magnitude remains OPEN, **S11c_PARTIAL_CLOSEOUT.md:31** [R9](#r9). Thus it is not rule-3 testable solely because a future comparison is named. Applicability of the magnetic-field identification is unclear; neither \(\Delta n(B)\) nor the QED coefficient is reproduced. |
| E19 | **Not addressed** | **Conditions:** the 2006 signal was withdrawn; its instrumental history is not an established vacuum effect. The withdrawn instrumental signal supplies no surviving new physical effect to classify as reproduced. Current model status remains the open optical-response map, **LIGHT_LEAKAGE_SCOPING.md:36** [R27](#r27); it has no demonstrated prediction of the historical apparatus artifact. **Rule 5:** the withdrawn signal is an apparatus artifact, not a surviving vacuum effect predicted by a native optical object. The cited open analogs or future labels supply no object for this specific measured response. |
| E20 | **Not addressed** | **Conditions:** stellar emission, magnetic geometry and propagation assumptions govern the indirect vacuum-birefringence inference. No native stellar surface-emission/magnetospheric polarization calculation is supplied in the audited current light records; **S11c_PARTIAL_CLOSEOUT.md:31,37** [R9](#r9). The model-dependent vacuum interpretation is not a derived brane result. **Rule 5:** the measured stellar polarization and its vacuum interpretation require undefined surface emission and magnetospheric response. The cited open analogs or future labels supply no object for this specific measured response. |
| E21 | **Not addressed** | **Conditions:** the published energy-resolved polarization is observed; its vacuum interpretation depends on emission and magnetospheric transport. No magnetar energy-resolved mode-conversion or emission calculation in the cited current results; **S11c_PARTIAL_CLOSEOUT.md:31,37** [R9](#r9). The observed angle swing is not the model's already computed defect birefringence. **Rule 5:** the observed magnetar spectrum/angle swing requires undefined surface emission and magnetospheric transport. The cited open analogs or future labels supply no object for this specific measured response. |
| E22 | **Required, no mechanism yet** | **Conditions:** E22 is an inferred helicity-delay bound under its stated gravitational-potential/emission assumptions, applied to light propagating near the brane's own mass structure; the observable map to \(\Delta\gamma_p\) must be specified. **Rule 4, taking precedence over rule 3:** the supplied polarization identification **`c_γ,1(x) ≡ c_γ,2(x) ≡ c_γ(x)`**, **S9b_SHARED_PHYSICS.md at commit `ede8aa21`:53–55**, states equal local speeds [P4](#p4). Polarization-independent eikonal delays that follow only from this identification belong here even though the spec supplies an optical ray object. S10 **does compute** the conditional two-direction transverse count and common homogeneous root, **S10_two_transverse_photons.md:73–86** [Q5](#q5); those results do not derive helicity-independent gravitational delay near a mass. The S9b common dispersion/rays are supplied, **:34–44**, while polarization-dependent propagation and subprincipal transport are outside/not computed, **:106–109,117–119** [P4](#p4), [Q4](#q4); **O2:125–127** likewise delivers no polarization transport [R17](#r17). No independently derived gravitational helicity-delay comparison is supplied by those records; neither the input's physical validity nor the measured PPN mapping is resolved here. |
| E23 | **Testable** | **Conditions:** the measured fact is common flyby deflection, not a helicity difference; use the spec's far-field, linear, stationary spherically symmetric setting and appropriate profiles/observable map. **Rule 3 object:** **S9b_SHARED_PHYSICS.md at commit `ede8aa21`:34–44,117–119,165–186** supplies the local optical dispersion/metric and its rays and names this observable [Q4](#q4). No step in the cited records compares that object with this experimental result. **Polarization independence is supplied**, `c_γ,1=c_γ,2=c_γ`, lines **53–55**, not reproduced; directional stiffness, thickness/bulk coupling and polarization-dependent propagation remain outside/OPEN, lines **106–109** [P4](#p4). **O2:104,125–127** does not supply the remaining optical transport [R17](#r17). The current physical profile/transport identification remains incomplete. |
| E24 | **Testable** | **Conditions:** the measured fact is common radar excess time/logarithmic slope, not a helicity difference; use the spec's far-field, linear, stationary spherically symmetric setting and appropriate profiles/observable map. **Rule 3 object:** **S9b_SHARED_PHYSICS.md at commit `ede8aa21`:34–44,117–119,165–186** supplies the local optical dispersion/metric and its rays and names this observable [Q4](#q4). No step in the cited records compares that object with this experimental result. **Polarization independence is supplied**, `c_γ,1=c_γ,2=c_γ`, lines **53–55**, not reproduced; directional stiffness, thickness/bulk coupling and polarization-dependent propagation remain outside/OPEN, lines **106–109** [P4](#p4). **O2:104,125–127** does not supply the remaining optical transport [R17](#r17). The current physical profile/transport identification remains incomplete. |
| E25 | **Not addressed** | **Conditions:** refraction at the measured material interface and the weak-measurement beam protocol, not gravitational spin Hall propagation. No derived polarization-dependent transverse beam position at a dielectric interface in the audited v3 results; bounded spin-Hall search [R31](#r31). The named S22 defect-birefringence item, **V3_STEP_PLAN.md:1183** [R19](#r19), does not establish this spatial shift. **Rule 5:** the measured spin Hall displacement depends on an ordinary air/glass optical interface and receiver protocol, not a supplied native defect/interface optical object. The cited open analogs or future labels supply no object for this specific measured response. |
| E26 | **Not addressed** | **Conditions:** optical reflection at a transparent dielectric interface; a brane/bulk material boundary is not already identified with that interface. No derived optical dielectric boundary law or Brewster-polarization condition is supplied by the current audited light records. S11b's faces are a brane/bulk material interface with uniform linear coupling, **S11b_interface_coupling_law.md:14–16,74–81** [R30](#r30); that record is not a computation of the optical experiment's transparent-material reflection. **Rule 5:** Brewster reflection depends on the undefined optical dielectric law of ordinary transparent materials; a brane/bulk material boundary is a different supplied object. The cited open analogs or future labels supply no object for this specific measured response. |
| E27 | **Not addressed** | **Conditions:** the measured wavelengths and clear/partly cloudy sky, with microscopic scattering/emission and receiver response needed for comparison. No polarized microscopic scattering/emission receiver law or sky-polarization pattern is obtained by the cited current records; acceleration radiation is explicitly open in the exploratory comparison, **docs/native_light_em_and_vortex_throat_interpretation.md:4214–4215** [R21](#r21), whose status is bounded by **S11c_PARTIAL_CLOSEOUT.md:37** [R9](#r9). **Rule 5:** the measured sky pattern depends on undefined atmospheric microscopic scattering/emission and cloud response. The cited open analogs or future labels supply no object for this specific measured response. |
| E28 | **Not addressed** | **Conditions:** the cosmological scattering/acoustic source and observed CMB E-mode sky pattern. No Thomson-scattering source, acoustic-to-CMB polarization prediction or E-mode spectrum is provided by the audited current results. The source/radiation and transport gaps are **docs/native_light_em_and_vortex_throat_interpretation.md:4214–4215** [R21](#r21) and **O2:125–127** [R17](#r17). **Rule 5:** the observed CMB E-mode pattern depends on undefined cosmological emission, Thomson scattering and acoustic-source physics. The cited open analogs or future labels supply no object for this specific measured response. |
| E29 | **Required, no mechanism yet** | **Conditions:** a tightly focused non-plane-wave light field propagating in the brane; the measured focal axial field is not an independent longitudinal plane-wave photon. **Rule 4:** the focal vector field/spot object has not been supplied, although light propagation lies within scope. **S10:182–197** supplies real cosine plane waves and in-plane field content, not that focused-beam field [R10](#r10). S11's separate compression branch/excluded normal field do not compute it, **S11:161–163** [R12](#r12), **SUBSTRATE_REQUIREMENTS.md:300–305** [R28](#r28). Hence this row does not claim a rule-3 comparison object merely from the plane-wave count. |
| E30 | **Not addressed** | **Conditions:** classify the measured quantum OAM preparation and entangled joint outcomes, not just classical wave amplitudes or a spatial mode count. **Rule 5:** that quantum emission/measurement physics is not defined in the present classical brane model; **V3_STEP_PLAN.md:1198–1204** explicitly assigns field quantization to a separate project [P8](#p8). No source here supplies the quantum state/outcome object this measurement constrains. A missing future quantum completion is not a requirement inside the current propagation/defect scope. |
| E31 | **Not addressed** | **Conditions:** the stellar emission, atmospheric scattering and magnetospheric assumptions in that tentative interpretation. No combined radio/X-ray stellar emission, polarized atmospheric scattering or strong-field magnetospheric transport calculation in the cited current results; **S11c_PARTIAL_CLOSEOUT.md:31,37** [R9](#r9), **O2:125–127** [R17](#r17). The tentative vacuum interpretation has no current numerical brane counterpart. **Rule 5:** the observed stellar X-ray/radio polarization depends on undefined atmospheric emission/scattering and magnetospheric transport. The cited open analogs or future labels supply no object for this specific measured response. |
| E32 | **In apparent conflict** | **Conditions:** the retained longitudinal branch couples to matter; for the thermal evidence it must also participate as an additional ordinary radiating mode under E33–E34's thermal occupation, dispersion and calibrated-coupling assumptions. **Experimental side:** “the two polarization degrees of freedom”, [James et al.](https://arxiv.org/pdf/quant-ph/0103121), and **E33's laboratory blackbody power** ([Quinn & Martin](https://doi.org/10.1098/rsta.1985.0058)) and **E34's CMB spectrum** ([Fixsen et al.](https://arxiv.org/pdf/astro-ph/9605054)) support the conditional ordinary two-mode description; they give no universal bound on uncoupled/unequilibrated fields. **Model side:** **S10:249–251** says “Maxwell has no counterpart to this surviving zero-frequency longitudinal direction” and calls it “a characterised departure” [P5](#p5). **S11:15** lifts it to a propagating mode; **161–163** leaves observability open; **173** says “A second cone is only a departure if matter COUPLES to it.” [P3](#p3) **Rule 2 takes priority:** the ledger-recorded departure gives conditional apparent conflict, rather than being demoted because thermal/measurement mechanisms are missing. No coupling or thermalization is established here. |
| E33 | **In apparent conflict** | **Conditions:** the extra longitudinal branch couples to matter and participates as an **additional ordinary thermal radiation mode** at the laboratory temperature, with the usual quantum occupation/energy, dispersion and calibrated emission/absorption response used for the two-mode inference. **Experimental side:** [Quinn & Martin](https://doi.org/10.1098/rsta.1985.0058) report “\(\sigma=(5.66967\pm0.00076)\times10^{-8}\)” in their original units (E33); the ordinary Planck normalization counts two radiating polarizations, under the assumptions stated in part 2. **Model side:** **S10:249–251** calls the surviving longitudinal direction “a characterised departure” from Maxwell [P5](#p5); **S11:15–17** gives a homogeneous three-finite-mode census, and **173** makes the departure conditional on matter coupling [P3](#p3). **Rule 2, the same as E32/E34:** a third ordinary thermal radiation channel is opposed to that conditional two-mode normalization. Missing field quantization/thermal response, **V3_STEP_PLAN.md:1198–1204** [P8](#p8), and the OPEN joint u/h object, **SUBSTRATE_REQUIREMENTS.md:263–305** [P8](#p8), do not override the earlier departure class. These sources do not establish the extra branch's thermal occupation/coupling; no unconditional blackbody exclusion or resolution is asserted. |
| E34 | **In apparent conflict** | **Conditions:** the extra branch couples to matter and participates as an **additional ordinary thermal radiation mode** under the usual occupation, energy and dispersion assumptions; its power must enter the calibrated CMB radiance rather than cancel in the sky/reference comparison. **Experimental side:** [Fixsen et al.](https://arxiv.org/pdf/astro-ph/9605054) report RMS deviations “less than 50 parts per million of the peak of the CMBR”; the fitted ordinary Planck law uses the two-polarization normalization (E34). **Model side:** **S10:249–251** calls the surviving longitudinal direction “a characterised departure” from Maxwell [P5](#p5); **S11:15–17** gives a homogeneous three-finite-mode census, and **173** makes the departure conditional on matter coupling [P3](#p3). **Rule 2, the same as E32/E33:** the ledger-recorded extra branch is in apparent conflict **conditional on the ordinary thermal-radiation identification**. FIRAS's differential sky/reference measurement does not set a universal 50-ppm extra-mode bound; source, reference and detection response matter. The quantization/thermal and full u/h status remains OPEN/outside the delivered classical calculation, **V3_STEP_PLAN.md:1198–1204; SUBSTRATE_REQUIREMENTS.md:263–305** [P8](#p8). Neither those conditions nor the conflict are resolved. |
| E35 | **Testable** | **Conditions:** the adopted Proca/MHD assumptions in E35 apply, and the constrained photon-mass parameter is identified with a gap of the native radiating branch; the observable identification is not earned by the supplied mode count. **Rule 3 object:** **S10_two_transverse_photons.md:237–252** supplies/computes the homogeneous gapless transverse spectrum [R10](#r10). No step compares it with this numerical mass/dispersion bound. This is not a reproduced photon-mass limit or a native plasma/emission calculation; matter/interface observability is unresolved, **S11_stray_longitudinal.md:154–163** [R12](#r12), and quantization is outside the delivered classical scope, **V3_STEP_PLAN.md:1198–1204** [P8](#p8). |
| E36 | **Testable** | **Conditions:** the paper's dispersion, plasma/host and cosmological assumptions in E36 apply, with its massive-photon dispersion interpreted as a gap/dispersion correction of the native radiating branch; the observable identification is not earned by the supplied mode count. **Rule 3 object:** **S10_two_transverse_photons.md:237–252** supplies/computes the homogeneous gapless transverse spectrum [R10](#r10). No step compares it with this numerical mass/dispersion bound. This is not a reproduced photon-mass limit or a native plasma/emission calculation; matter/interface observability is unresolved, **S11_stray_longitudinal.md:154–163** [R12](#r12), and quantization is outside the delivered classical scope, **V3_STEP_PLAN.md:1198–1204** [P8](#p8). |
| E37 | **Testable** | **Conditions:** the paper's dispersion, plasma/host and cosmological assumptions in E37 apply, with its massive-photon dispersion interpreted as a gap/dispersion correction of the native radiating branch; the observable identification is not earned by the supplied mode count. **Rule 3 object:** **S10_two_transverse_photons.md:237–252** supplies/computes the homogeneous gapless transverse spectrum [R10](#r10). No step compares it with this numerical mass/dispersion bound. This is not a reproduced photon-mass limit or a native plasma/emission calculation; matter/interface observability is unresolved, **S11_stray_longitudinal.md:154–163** [R12](#r12), and quantization is outside the delivered classical scope, **V3_STEP_PLAN.md:1198–1204** [P8](#p8). |
| E38 | **In apparent conflict** | **Conditions:** the retained extra longitudinal branch couples to matter and is cosmologically populated, remains relativistic at the relevant epochs, and contributes **additional gravitating radiation** beyond the standard photon/neutrino contents, with energy density and perturbation evolution mapped to E38's \(N_{\rm eff}\) parametrization. Its occupation, dispersion and temperature/decoupling history must give an extra density beyond the fitted allowance; a sufficiently weak population or earlier decoupling is not generically excluded. **Experimental side:** [Planck 2018 parameters, abstract and §7.5.2](https://arxiv.org/pdf/1807.06209) report “\(N_{\rm eff}=2.99\pm0.17\)” at 68% under E38's assumptions, with no established extra-radiation requirement. **Model side:** **S10:249–251** calls the surviving longitudinal direction “a characterised departure” [P5](#p5); **S11:15–17** supplies the homogeneous three-finite-mode census, **:161–163** leaves observability open, and **:173** says “A second cone is only a departure if matter COUPLES to it” [P3](#p3). **Rule 2, the same as E32–E34:** the ledger-recorded extra branch gives conditional apparent conflict if it participates as the additional relativistic radiation specified here. Missing quantization/occupation, **V3_STEP_PLAN.md:1198–1204** [P8](#p8), does not demote that earlier departure to rule 5. These sources supply no cosmological population, decoupling history or \(\Delta N_{\rm eff}\); no unconditional species exclusion or resolution is asserted. |

Classification totals: **38 entries — reproduced 1; in apparent conflict 7; testable 10; required, no mechanism yet 5; not addressed 15.** E01 reproduces the conditional in-plane transverse count, without selecting the physical dimension or completing the radiation-state census. E12 is testable by the same supplied/computed propagation object used for E10–E11 and E13–E15; the CMB source/statistical mapping remains a shared condition. The seven apparent conflicts are conditional: E13–E15 on a real cosmic signal and the selected parity/transport scope; E32–E34 on matter coupling and, for the blackbody comparison, ordinary thermal-radiation participation; E38 additionally on the extra branch's cosmological relativistic population, energy/perturbation mapping and decoupling history. E22's equal-speed input invokes rule 4's precedence over rule 3. No open condition or reported signal is resolved here.

## Appendix A. Repository commands and literal output

All commands below were run from /var/projects/toy_physics. These are read-only retrievals. Each R-, P- or Q-number linked in parts 3–4 resolves to its exact command and captured terminal output here. Line numbers are produced by nl -ba or rg -n. The evidence is shown as retrieved, including the source's own markup and historical qualifications; none of it is a new review verdict.

Only evidence used for repository claims is included. Empty search output is represented by an empty output block, with its exit code stated outside the block. No truncated retrieval is used as evidence.

<a id="r6"></a>
<details>
<summary>R6 — command and literal output</summary>

Command:

``````bash
nl -ba research/pde_ledger_v3/steps/S9_light_requires_shear.md | sed -n '1,120p;280,365p'
``````

Exit code: 0.

Literal output:

``````text
     1	# S9 · Light's requirement on the medium — shear on the brane, none in the bulk
     2	
     3	**Sector 1 (light), step 1.** Walked side by side with the user, 2026-08-01.
     4	⚠ First step banked under the **requirements-first** direction: light states what it needs; the medium's
     5	structure is not assumed in advance.
     6	
     7	---
     8	
     9	## What it is
    10	
    11	Light is a transverse wave. A GNLS superfluid cannot carry one. So light's first act is to state what
    12	the medium must have that the GNLS does not.
    13	
    14	## What it does
    15	
    16	Introduces the two constants the entire light sector rests on, and fixes the light cone as their ratio.
    17	Everything downstream in this sector is a consequence of these plus the curl-only form (S10, S11).
    18	
    19	## The argument, in three moves
    20	
    21	**1 — What Maxwell demands.** Two transverse polarisations · one speed · no longitudinal mode.
    22	
    23	**2 — What a GNLS has.** Linearising Madelung about uniform `ρ₀` gives exactly **one** mode:
    24	
    25	```
    26	ω² = c_s²k² + (ħ²/4m²)k⁴        c_s² = (1/m)·dP/dρ = 5Kρ₀⁴/m
    27	```
    28	
    29	and it is **longitudinal**. The obstruction is not incidental — it is representation-theoretic:
    30	
    31	> *"a single-component **SCALAR** superfluid cannot carry transverse light, period … `Ψ = √ρ e^{iθ}` is
    32	> one complex scalar. Its excitations are `δρ` and `δθ`, both spin-0/longitudinal; one broken `U(1)` →
    33	> one Goldstone phonon. **A scalar's spectrum cannot contain a spin-1 photon.** … This is
    34	> **dimension-independent** for a pure scalar"*
    35	> — `software/stage1_solver/decisions/15_em_medium_native_physical_picture.md:29-40`
    36	
    37	⚠ **Status of that argument:** cited from an external review; ⛔ **no script in this repo executes it.**
    38	A weaker in-repo form exists — `v = (ħ/m)∇θ ⇒ ∇×v ≡ 0`, so superfluid flow is curl-free — but it is
    39	never wired into a no-shear proof (`docs/model_map.md:57` states the flow form only).
    40	
    41	**3 — Therefore.** ⭐ **Light requires structure the GNLS does not contain.**
    42	
    43	> *"**Brane elasticity is NOT a consequence of the GP/NLS mean-field (C1).** A GP/NLS superfluid is a
    44	> *fluid* — zero shear modulus. The brane's shear rigidity therefore requires the **substructure**
    45	> (constituents/cohesion *beneath* the mean-field), which the GP/NLS equation does not contain. Honest
    46	> framing: **"GP/NLS as the effective/coarse-grained medium + a deeper substructure that supplies the
    47	> brane elasticity."***
    48	> — `software/stage1_solver/decisions/15_em_medium_native_physical_picture.md:230-233`
    49	
    50	⚠⚠ **ONE SOURCE, NOT TWO — state it or the argument reads as stronger than it is.** Moves 2 and 3 both
    51	come from `decisions/15`. ⛔ A draft of this record cited move 3 to `research/em_fields/paper/em_fields.tex:230`,
    52	which is `\label{eq:B-vorticity}` — the phrase occurs **zero** times in that file (F4 class: right line
    53	number, wrong file). ⇒ The whole *"a GNLS cannot carry light"* argument rests on **a single external-review
    54	document with no executing script**, ⛔ not on two independent sources and ⛔ not on the model's own paper.
    55	
    56	## ⭐⭐ The requirement, stated in full — and it is TWO-SIDED
    57	
    58	The medium has **substructure** whose coarse-grained description is the GNLS (charter-level standing
    59	postulate, ⛔ not light's to introduce). Light requires of that substructure:
    60	
    61	| | requirement | what it buys |
    62	|---|---|---|
    63	| **ordered phase** (`χ_B = 1`) | carries **MacCullagh curl-only** shear | photons **exist** |
    64	| **disordered phase** (`χ_B = 0`) | carries **no** shear | photons stay **confined** |
    65	
    66	⭐ **Both halves are light's own bill.** ⛔ The second is not magnetism's, as an earlier draft had it:
    67	without it a photon has somewhere else to go, light leaks off our 3D space, and that is an energy sink
    68	with no observational room. ⚠ **Stressed hardest at the throat**, where the brane is bent into `±w`.
    69	
    70	⭐ **A third consumer of the second half:** the recorded throat-support mechanism is an *outward trapped
    71	standing-wave pressure* whose trapped mode is a **brane-shear standing wave**. If the bulk carried
    72	shear that mode radiates away, the outward pressure vanishes, and the throat closes. ⇒ **Bulk
    73	shear-freeness is what makes the geon possible at all.**
    74	
    75	## The equations
    76	
    77	```
    78	L = ½ ρ_br (∂_t u)²  −  ½ μ_R (curl u)²                     (the brane's transverse sector)
    79	c_γ² = μ_R / ρ_br                                            (R4 — the light cone)
    80	```
    81	
    82	## Dimensions — derived here, then checked against the register
    83	
    84	`L` is an energy density on the 3D brane, `[L] = M L⁻¹T⁻²`. With `[u] = L`:
    85	
    86	```
    87	[∂_t u] = L T⁻¹                     [curl u] = L/L = 1  (dimensionless)
    88	
    89	ρ_br (∂_t u)²  :  [ρ_br]·L²T⁻² = M L⁻¹T⁻²   ⇒   [ρ_br] = M L⁻³      = [-3, 0, 1]
    90	μ_R  (curl u)² :  [μ_R]·1      = M L⁻¹T⁻²   ⇒   [μ_R] = M L⁻¹T⁻²    = [-1,-2, 1]
    91	
    92	R4 :  μ_R − ρ_br = [-1-(-3), -2-0, 1-1] = [2,-2,0] = 2×[1,-1,0]  ⇒  c_γ = [1,-1,0]  ✓
    93	```
    94	
    95	⭐ **Independent agreement:** these were derived from the Lagrangian *before* consulting the register,
    96	which records `μ_R` as `M L⁻¹ T⁻²` and `ρ_br` as `M L⁻³` (`research/pde_ledger_v2/notes/parameter_register.md:137`, `:138` — ⚠ **v2's** register, cited by locus; v3 does not have a `notes/` tree).
    97	
    98	### ⭐⭐ Rebuilt 2026-08-04: the dimensions are now solved SYMBOLICALLY IN `D`, from the action
    99	
   100	Both engines now extract the **derivative multi-order of every field factor from the Lagrangian's
   101	expression tree** and solve the resulting linear system. ⛔ Nothing is read off a table and ⛔ no operator
   102	is special-cased by name.
   103	
   104	```
   105	[ρ_br] = (−D,   0, 1)          [μ_R] = (2−D, −2, 1)
   106	[μ_R] − [ρ_br] = (2, −2, 0)    ⭐ INDEPENDENT OF D
   107	at D = 3:  (−3, 0, 1)  and  (−1, −2, 1)   ✓ registry
   108	```
   109	
   110	⭐⭐ **This retires a recorded S9 "blind spot" by explaining it.** The old audit's dimension check was
   111	insensitive to the assumed brane dimension, filed as a weakness. It is not a weakness but an **identity**:
   112	the speed's dimension **cannot** see `D`, because the difference is `(2,−2,0)` for every `D`.
   113	
   114	⭐ **And the block is now able-to-fail, which it was not.** A control with a **different derivative count**
   115	(flexural, `(∇²u)²`) moves the stiffness dimension to `(4−D,−2,1)`; a control with an **undifferentiated
   116	field** (`−½ μ_G u·u`) yields `(−D,−2,1)`. ⚠ Before the rebuild the block emitted **byte-identical output
   117	under a change of the action's form**, because the derivative counts were hand-encoded.
   118	
   119	### ⭐⭐ FIVE independently-built engines derive these two dimensions, and three of them read no register
   120	
   280	only that both engines typed `mu − rho`, rather than `rho − mu`; `SPEED_DIMENSION_DIFFERENCE` buys only
   281	that both posited the same velocity-squared reference.
   282	
   283	⚠ Under the gradient-elastic form control, **10 of the 12 standard-name rows are byte-identical**; only
   284	`FACTORED_DETERMINANT` and `FULL_ROOT_MULTISET` move. ⇒ The table's discriminating power against the
   285	ordinary-elastic alternative is **two rows**.
   286	
   287	⚠ **The arithmetic here was wrong in the source this record was written from** — a review leg reported
   288	*"9 of the 10 exported standard rows"*, which with two movers accounts for eleven rows out of ten. The
   289	builder flagged the inconsistency rather than silently picking a number, and it was **re-measured**: of the
   290	twelve standard names, **two move and ten do not**. (Separately, ten of the twelve reach the export; the
   291	flexural and bare-field rows belong to control packages `X7`/`X8`. Conflating those two tens produced the
   292	error.) ⭐ The two-row conclusion is unaffected.
   293	
   294	⭐⭐ **There is a measured common-mode blindness behind three dimension rows.** Both engines independently
   295	type `[u] = L`, and neither emits it into the standard-name comparison. Doubling it in a scratch copy of
   296	each engine moves `INERTIA_COEFFICIENT_DIMENSION`, `STIFFNESS_COEFFICIENT_DIMENSION`, and
   297	`BARE_FIELD_COEFFICIENT_DIMENSION` to the same wrong values; both engines exit 0 and still agree. Those
   298	three rows buy the derivative-multiorder extraction and the solve, ⛔ **not the dimension**.
   299	
   300	This blindness is an identity, not an empirical miss. If the field dimension is left symbolic,
   301	`[mu_R] − [rho_br]` is independent of it in all three axes. ⇒ ⛔ No computation inside S9, in either
   302	engine, can ever detect a wrong `[u]`. The error becomes detectable only when a consumer imports the entry
   303	instead of re-declaring it, and S10 currently re-declares it.
   304	
   305	⚠ **The rebuild also introduced three open defects:**
   306	
   307	- `wavevector_norm_dimension` names the wrong object. The name denotes `dim(|k|) = [-1,0,0]`, while the
   308	  exported value `[-2,0,0]` is `dim(k·k)`. The value is right and the name is wrong in both engines and in
   309	  the name-keyed file.
   310	- The placeholder-naming class contains **eight** entries, not one. In addition to that key, five keys
   311	  still carry the SymPy-only `q`, while `dim_energy_density` and `dim_squared_velocity` have exact
   312	  Wolfram counterparts under different names.
   313	- `q_dimension` is unpinned inside SymPy: typing it to anything lets the engine exit 0 with a wrong ledger
   314	  value. Two residuals detect the change but are unguarded; the guard needs scoping because the `X3`
   315	  residual is legitimately nonzero under a form control. ⇒ This remains a spec question.
   316	
   317	⚠ **The export has no provenance field.** `wavevector_norm_dimension` therefore reaches a consumer as a
   318	SymPy `PREMISE` with no indication that the Wolfram engine derives it. The flat record cannot express that
   319	asymmetry.
   320	
   321	⛔⛔ **Nothing in the repository performs the cross-engine standard-name lookup.** The naming convention
   322	makes comparison a join, but no artifact executes that join, so the export is ⛔ **not checked against the
   323	Wolfram engine**. Re-pointing one line in an engine's name→object table turns the light cone into an
   324	`omega^2`, wrong by `k²`, while the determinant, root multiset, and speed-dimension residual remain
   325	unmoved and the engine exits 0. The dimension check cannot catch the substitution because `q_dimension`
   326	is typed at exactly `−2L`; the Wolfram engine has no `q` and does not move. ⇒ The cross-engine comparison
   327	on that row is the only instrument that would see the error, and it does not exist.
   328	
   329	## ⛔ WHAT THIS STEP STILL DOES NOT ESTABLISH
   330	
   331	- ⛔ **P2 — that a scalar superfluid carries no transverse mode.** Cited from one external review; ⛔ **no
   332	  script executes it**, and the rebuild did not add one. It remains a **supplied premise**.
   333	- ⛔ **The absence of a propagating longitudinal wave is ASSUMED, not derived.** The curl-only action sets
   334	  the longitudinal restoring stiffness to **zero by construction**, so `ω² = 0` there is the postulate
   335	  restated. ⛔ It removes the restoring *force*, ⛔ **not the degree of freedom** — which is why the
   336	  longitudinal test is emitted at every root.
   337	- ⛔ **Bulk shear-freeness** is postulated; the bulk is absent from the action, so nothing here tests it.
   338	- ⚠ **The assumption set is computationally INERT.** Removing it entirely leaves every computed value
   339	  byte-identical. ⭐ The physics is nonetheless correctly bounded — by generic rank **plus explicitly
   340	  emitted exceptional loci**, and a reviewer verified that enumeration is **complete** for the anisotropic
   341	  control (roots collide exactly where `ρ_z = ρ_br` or `k_x = k_y = 0`, and both are emitted).
   342	  ⇒ ⭐ **Read the polarisation values as generic-locus values with the exceptions listed.**
   343	- ⚠ **The `.py` was written AFTER the `.wl`**, so this is ⛔ not the blind-first ordering. Construction was
   344	  independent (different naming, granularity, idioms; a hand-rolled Euler–Lagrange where SymPy has a
   345	  built-in) and the builder was barred from reading the `.wl` — ⭐ but agreement here is **weaker evidence**
   346	  than at a step where both engines predate any comparison.
   347	- ⚠ One tag remains **unparsed** by the consumer: the anisotropic third root's sign is a `Piecewise`,
   348	  because the sign is conditional on the parameter domain and the assumptions are inert.
   349	- ⚠ Limits taken: sharp zero-width sheet · `v₀ → 0` · no dissipation · frequency-independent moduli ·
   350	  continuum limit · amplitude → 0.
   351	
   352	## Registry additions (executed at this step)
   353	
   354	`Q.brane.rho_br` `[-3,0,1]` · `Q.brane.mu_R` `[-1,-2,1]` · `Q.brane.c_gamma` `[1,-1,0]` · relation
   355	**`R4`** `c_γ = √(μ_R/ρ_br)`. ⭐ `R4` is the corpus's existing name for this relation
   356	(`notes/parameter_register.md:271`). ⚠ `c_γ` re-enters here with provenance, having been deleted from
   357	the *medium* block at S0.5 — that round trip is the mechanism working.
``````

</details>

<a id="r9"></a>
<details>
<summary>R9 — command and literal output</summary>

Command:

``````bash
nl -ba research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md
``````

Exit code: 0.

Literal output:

``````text
     1	# S11c closure and handoff — PARTIAL
     2	
     3	S11c asked how much light converts into other motion where the brane's thickness or material properties change. **It closes PARTIAL: no accepted nonuniform light-loss number was obtained.** Both reflected and transmitted light count as survival. No S11c calculation or review resumes automatically. [Canonical d record, opening and scope](S11c_d_profile_conditioned_scattering.md).
     4	
     5	## Where we stand
     6	
     7	On a uniform brane, the stated linear model decouples light from the bulk. Independent SymPy/Wolfram work and scoped reviews support that result; stable propagation needs nonnegative stiffness. A later SymPy check found finite modes, nonzero energy current and zero face drives (light does not push on the bulk) at, above and below the tested light/sound speed matches, with the bulk at rest and other inputs fixed; Claude/Grok cleared its method, not a fresh review of its executed result. **CONDITIONAL.** [S11b, transverse mode/reviews](S11b_interface_coupling_law.md); [d, uniform check](S11c_d_profile_conditioned_scattering.md#what-the-uniform-check-establishes).
     8	
     9	Where the brane changes, the equations contain conversion routes, but their size remains unknown. The held-profile omega=3 benchmark at speed ratio approximately 0.12 is **UNRESOLVED**: final imbalances lie within its empirical 1 ppm floor, not a physical bound. This is not a matched-speed defect result. [d, benchmark and controls](S11c_d_profile_conditioned_scattering.md#the-omega3-benchmark-numerical-loss-unresolved).
    10	
    11	## Records and open ownership
    12	
    13	| Record / qualified result | Unresolved / owner |
    14	| --- | --- |
    15	| [a: interface shape derivatives](S11c_a_interface_shape_derivatives.md), per-family comparison and limits | Geometric/substrate dependencies remain conditional; affected S11c consumers if reused. |
    16	| [b: variable-coefficient operator](S11c_b_variable_coefficient_operator.md), original per-engine close plus later scoped repairs | Full cross-engine residual/composition: **S11c-b**. Original checks describe pre-repair equations. |
    17	| [c1: curved bulk response](S11c_c1_curved_bulk_closure.md), kernel/pressure/flux agreement on measured domains | Traction, energy, operand and deferred family comparisons: **S11c-c1/c2**. |
    18	| [c2: self-energy fold](S11c_c2_self_energy_fold.md), original wiring close plus later scoped repairs | Direct mixed entry and full composition: **S11c-c2/d**. Original “VALUES unaffected” describes the pre-repair operator. |
    19	| [d: supplied-profile scattering](S11c_d_profile_conditioned_scattering.md), selected uniform success and unresolved defect benchmark | Complete response/work/loss: **S11c-d if reopened**; observable/strong edges: **S11c-e**. |
    20	| [Clean-condition symmetry rule](S11c_d_profile_conditioned_scattering.md#clean-condition-result-and-its-owner), measured on SymPy only | No-leak interpretation needs drain frozen; live conversion belongs to **S12**, with no automatic transfer. |
    21	| Supplied-profile calculations do not determine a real particle. [Plan, Q2/Q3/S20a/S22](../V3_STEP_PLAN.md) | Holder/charge response: **Q2/Q3/S22**; speed-ratio calibration: **S20a**, derivation: **S22/R10**; drain/return data: **S12**. |
    22	
    23	## Claims future work must not inherit
    24	
    25	- **UNRESOLVED:** c2 sets the direct three-leg `[0,2]` entry to zero while keeping the iterated product. Retaining `eta*sigma_W` means order counting does not justify that zero. **Whether the complete closed response has a nonzero direct term is unresolved.** This is not proof of a missing physical term or its effect on the benchmark. [Joint disposition, mixed-term issue](../_measurements/S11c_upstream_repair_review_joint_disposition.md); [d, equations](S11c_d_profile_conditioned_scattering.md#equations-what-was-repaired-and-what-remains-open).
    26	- **UNRESOLVED:** the omega=3 benchmark establishes neither “leakage below 1 ppm” nor physical gain. Source identities, finite packet checks and first-order flux cancellation do not constitute a complete loss calculation. [d, benchmark and intermediate results](S11c_d_profile_conditioned_scattering.md).
    27	- **OPEN — material audit/Q2/S22:** uniform confinement is derived within the supplied linear model; nonuniform confinement and the physical admissibility of rotational stiffness are not established. That is the useful but limited MacCullagh distinction. [Native interpretation, §§5.7, 14.3](../../../docs/native_light_em_and_vortex_throat_interpretation.md); [framing essay, differentiators](../../../docs/s11_maccullagh_differentiation.md).
    28	
    29	## Actual downstream uses and next step
    30	
    31	**OPEN:** S11c-e owns the conversion observable and must recover d's weak limit in a future finite-contrast calculation. The plan's charge preamble asks whether linear nonuniform physics can support charge; PARTIAL closure answers neither that question nor all Q-sector work. **S22 items 1/4** concern conversion/birefringence and supplied defect profiles; physical magnitude needs interior input. Items 2/3 name S11b and S5–S7 instead. [Split, table/N5](../directives/S11c_decisions.md); [d spec, §3d](../directives/S11c_d_SHARED_PHYSICS.md); [plan, charge preamble/S22](../V3_STEP_PLAN.md).
    32	
    33	**S12 is next:** dynamical bulk-to-brane order conversion and separate boundary data. It does not require a numerical S11c loss factor and does not automatically settle light leakage. [Plan, S12](../V3_STEP_PLAN.md).
    34	
    35	## Exploratory notes, paused
    36	
    37	The [throat/EM handoff in docs/](../../../docs/light_em_investigation_handoff.md) indexes the preserved support, oriented-core, conversion/work and elastic-reference notes. These are **exploratory inputs for S22/Q2**, with S12/Q3 connections, not adopted laws or solved throats. Their review disagreements remain with their own assessments. [Phase 1 review, C6](../cleanup_2026_10/PHASE1_REVIEW.md).
    38	
    39	The full history lives in `archive/pre-cleanup-2026-10-04`; the [cleanup record](../cleanup_2026_10/FINDINGS.md) separates findings from unresolved work. No archived result is upgraded by this closure.
``````

</details>

<a id="r10"></a>
<details>
<summary>R10 — command and literal output</summary>

Command:

``````bash
nl -ba research/pde_ledger_v3/steps/S10_two_transverse_photons.md | sed -n '71,136p;172,197p;237,252p;645,715p'
``````

Exit code: 0.

Literal output:

``````text
    71	The measured generic result is:
    72	
    73	> For the supplied curl-only in-plane action, at nonzero wavevector and away
    74	> from allowed exceptional strata, the nonzero root has D − 1 transverse null
    75	> directions in every measured MAIN case D = 2, 3, 4, 5. Thus the D = 3 member
    76	> has two transverse directions.
    77	
    78	The form controls prohibit a stronger reading of that sentence. At both
    79	dimensions where the form-changing packages and MAIN were emitted, D = 3 and
    80	D = 4, FULLGRAD has the same nonzero root and the same D − 1 transverse
    81	nullity as MAIN. The transverse count therefore cannot be attributed to the
    82	curl-only form. What that form determines here is the disposition of the one
    83	remaining, longitudinal direction: unlike FULLGRAD, it leaves that direction
    84	at the zero root instead of putting it on the propagating root. DIVONLY, which
    85	reverses the two sectors, confirms that the stiffness form can change the
    86	count even though curl-only is not what uniquely produces D − 1.
    87	
    88	A separate condition comes from the inertia structure. `XFORM_ANISO` is not a
    89	stiffness-form control and does not multiply the inertia by one scalar: it
    90	keeps `S_curl` and changes the kinetic quadratic form on one distinguished
    91	axis to `s_ρ(∂_t u_1)² + Σ_{j=2..D}(∂_t u_j)²`, with `s_ρ > 0` and
    92	`s_ρ ≠ 1`. Thus the isotropic inertia value is not admitted by this control.
    93	Writing each root as `N2 / N3`, where N2 is total nullity and N3 is the
    94	basis-independent exactly-transverse nullity, the original and focused readings are:
    95	
    96	| D | case | roots in emitted order: N2 / N3 | total N2 / N3 over propagating roots |
    97	|---:|---|---|---:|
    98	| 3 | MAIN, generic | zero `1 / 0`; positive `2 / 2` | `2 / 2` |
    99	| 3 | ANISO, generic | zero `1 / 0`; positive `1 / 1`; positive `1 / 0` | `2 / 1` |
   100	| 3 | ANISO, parallel stratum | zero `1 / 0`; positive `2 / 2` | `2 / 2` |
   101	| 3 | ANISO, perpendicular stratum (new) | zero `1 / 0`; positive `1 / 1`; positive `1 / 1` | `2 / 2` |
   102	| 4 | MAIN, generic | zero `1 / 0`; positive `3 / 3` | `3 / 3` |
   103	| 4 | ANISO, generic | zero `1 / 0`; positive `2 / 2`; positive `1 / 0` | `3 / 2` |
   104	| 4 | ANISO, parallel stratum | zero `1 / 0`; positive `3 / 3` | `3 / 3` |
   105	| 4 | ANISO, perpendicular stratum (new) | zero `1 / 0`; positive `2 / 2`; positive `1 / 1` | `3 / 3` |
   106	
   107	Therefore the measured inertia-form change does not remove a propagating
   108	mode: at both shared dimensions its total N2 nullity over the positive roots
   109	remains D − 1, exactly as in MAIN. In the same breath, it does generically
   110	reduce the number of those modes that are exactly transverse from D − 1 to
   111	D − 2: ANISO puts N3 nullity D − 2 on one positive root and zero on its
   112	additional positive root, instead of MAIN's N3 nullity D − 1 on one positive
   113	root. MAIN itself is unchanged at every measured D; this is an added condition
   114	on the headline's generality, not a retraction of the MAIN result.
   115	
   116	The parallel ANISO stratum is `k_2 = ... = k_D = 0` with `k_1 ≠ 0`, so the
   117	wavevector has zero obliquity to the distinguished inertia axis. There the two
   118	generic positive branches coalesce, and the rerun restores N2 / N3 totals
   119	`2 / 2` at D = 3 and `3 / 3` at D = 4. This is distinct from restoring
   120	isotropic inertia: `s_ρ = 1` remains excluded. The scalar stiffness-coefficient
   121	control XCOEF_SCALE and the stiffness-sign control SIGNFLIP leave N3 unchanged
   122	in their measured cases, but ANISO shows that this control-specific result is
   123	not a general rule that only a stiffness-form change can move the
   124	exactly-transverse count. No control in this sweep varies the in-plane field
   125	content, the rest-background premise, dissipation, or the quadratic
   126	truncation; those remain supplied conditions.
   127	
   128	XCOEF_SCALE nevertheless moves the nonzero root's coefficient, and SIGNFLIP
   129	moves that root's sign and stability; their unchanged N3 counts do not mean
   130	that nothing physical moved.
   131	
   132	The perpendicular stratum is `k_1 = 0` with a nonzero remaining wavevector.
   133	Here the positive branches remain distinct, but the extra branch's N3 rises
   134	from zero to one. Its N2 is unchanged. Both exceptional directions have total
   135	propagating N2 / N3 equal to `(D − 1) / (D − 1)`; the Lean classification proves
   136	this throughout those directions under the supplied assumptions.
   172	The physical selection D = 3 is not made in S10. The live S10 computation keeps
   173	D symbolic for dimensions and evaluates an indexed sweep at D = 2, 3, 4, 5.
   174	The new Lean baseline proof establishes the conditional map D ↦ D − 1 for
   175	arbitrary finite D with a nonzero wavevector, extending the measured sweep.
   176	Neither establishes which D nature selects.
   177	
   178	## What was supplied and what was computed
   179	
   180	The shared specification supplies:
   181	
   182	- a D-component in-plane displacement u, with its separation from every
   183	  out-of-plane field inherited rather than tested;
   184	- the real cosine plane-wave ansatz;
   185	- positive inertia and stiffness coefficients, nonzero real wavevector,
   186	  unstrained rest background, no dissipation, and linear response;
   187	- the curl-only stiffness density
   188	
   189	      S_curl = (1/2) Σ_i Σ_j (∂_i u_j − ∂_j u_i)²
   190	
   191	  in
   192	
   193	      L = (ρ_br/2) Σ_j (∂_t u_j)² − (μ_R/2) S_curl.
   194	
   195	These premises and the action are at
   196	directives/S10_SHARED_PHYSICS.md:13-28, :30-47, and :82-107. The action is an
   197	input, not a result of the mode count.
   237	The raw root formula is 0 and μ_R |k|² / ρ_br. Under the supplied
   238	positive-coefficient and nonzero-wavevector assumptions, the latter is
   239	positive. The zero-root null direction is longitudinal and has no transverse
   240	nullity; the positive root has D − 1 total and D − 1 transverse null
   241	directions. This paired raw reading is not a blanket comparator pass for every
   242	Q3 expression; the joined integer families below are the cross-engine result
   243	claimed here.
   244	
   245	The zero root retains a degree of freedom; it does not remove one. What the
   246	curl-only stiffness removes is the restoring stiffness for the longitudinal
   247	direction. The per-root nullities make the distinction explicit:
   248	`1 + (D − 1) = D` in every MAIN dimension measured, so all D amplitude
   249	directions remain in the spectrum. Within the light sector, Maxwell has no
   250	counterpart to this surviving zero-frequency longitudinal direction. That is
   251	a characterised departure, and its onward disposition belongs to S11
   252	(`stray_longitudinal`). S10 assigns it no further interpretation.
   645	the inertia-coefficient dimension, stiffness-coefficient dimension, and their
   646	difference. They are not mutually independent. S10 freshly derives its two
   647	coefficient dimensions, but it imports from S9 the dimension symbol, the
   648	coefficient symbols, the unit objects, and the supplied field-dimension
   649	premise; its third record is then constructed from the first two. The guard is
   650	therefore blind to the whole class of mutations that leave the three dimension
   651	records equal on both sides. That class includes every dimension-preserving
   652	action or coefficient mutation, not only the coefficient ablation that was
   653	run, and every common-mode premise mutation accompanied by a consistent move
   654	of the upstream records. The form ablation and coefficient ablation are one
   655	caught member and one missed member of that class, not an exhaustive boundary.
   656	
   657	To reproduce either ablation, copy the five named files to a scratch directory,
   658	apply the single action mutation shown above to the scratch S9 engine, run the
   659	scratch S9 engine to regenerate its export, remove only the scratch S10 export
   660	target, and then run the scratch S10 engine. No ablated artifact is part of the
   661	committed ledger.
   662	
   663	## Live comparator: measured scope, not a blanket verdict
   664	
   665	The comparator joins only mechanically matched names, compares sequences in
   666	order, and canonicalizes null-space bases as spans
   667	(scripts/out/S10_cross_engine_comparator.out:1-5). Its measured summary at
   668	:9288-9313 is:
   669	
   670	| quantity | count |
   671	|---|---:|
   672	| Python names | 4233 |
   673	| Wolfram names | 2983 |
   674	| shared names | 562 |
   675	| compared shared names | 552 |
   676	| unparsed shared names | 10 |
   677	| agreements | 388 |
   678	| bare-integer agreements | 215 |
   679	| empty-container agreements | 14 |
   680	| symbolic or structured agreements | 159 |
   681	| disagreements | 164 |
   682	| naming-only disagreements | 23 |
   683	| representational disagreements | 13 |
   684	| content divergences | 128 |
   685	| shape or type mismatches | 109 |
   686	| route-token divergences | 13 |
   687	| numeric or algebraic residuals | 6 |
   688	| null-space bases compared | 26 |
   689	| duplicate rows | 0 |
   690	| format issues | 10 |
   691	
   692	For a value-comparable width, this record counts only 388 PASS rows plus the 6
   693	rows with a genuine numeric or algebraic residual: 394. The remaining 168
   694	shared names are 23 naming-only, 13 representational, 109 shape/type, 13
   695	route-token, and 10 unparsed rows. They do not establish value agreement.
   696	Likewise, the thousands of engine-only names are not cross-engine evidence.
   697	
   698	The six genuine residual rows are one MAIN sign row and five control rows:
   699	
   700	- MAIN D = 5 root 2 sign:
   701	  1 − sign(k1²+k2²+k3²+k4²+k5²), at comparator :938-943;
   702	- XCOEF_SCALE D = 3 Q6 unknown-coefficient count: −3, at :1065-1070;
   703	- ANISO D = 3 Q6 unknown-coefficient count: −3, at :1308-1313;
   704	- ANISO D = 3 root 3 sign:
   705	  undecidedUnderJointAssumptions − 1, at :1498-1503;
   706	- ANISO D = 4 Q6 unknown-coefficient count: −3, at :1592-1597;
   707	- ANISO D = 4 root 3 sign:
   708	  undecidedUnderJointAssumptions − sign(k1² sRho+k2²+k3²+k4²), at
   709	  :1782-1787.
   710	
   711	These are unresolved comparator rows. In particular, the MAIN D = 5 positive
   712	sign is not cross-engine established even though the raw premise
   713	Σ k_i² > 0 makes the mathematical sign positive; this record does not convert
   714	that inference into a comparator pass.
   715	
``````

</details>

<a id="r11"></a>
<details>
<summary>R11 — command and literal output</summary>

Command:

``````bash
nl -ba research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | sed -n '220,355p'
``````

Exit code: 0.

Literal output:

``````text
   220	  (`steps/S11c_PARTIAL_CLOSEOUT.md:7`; `steps/S11c_d_profile_conditioned_scattering.md:9,23`). This is an
   221	  unresolved tension between the driven reference state and the recorded rest-state calculations. The
   222	  rest-state obligation remains here; live drain/return data are `R-S12-02`, which also records the
   223	  uncarried correction and its record-named owner, and the power supply is `R-S12-01`. No equilibrium result
   224	  is promoted to a result about the driven state.
   225	
   226	### R-S8-01 — the form of the brane's quadratic stiffness functional
   227	
   228	- **source** S9, S10; S11 (`steps/S11_stray_longitudinal.md`, moves 2–3 / FORM control);
   229	  S11b-B (`steps/S11bB_interface_assembly.md`, energy basis / breathing stability); unified S11b
   230	  (`steps/S11b_interface_coupling_law.md`, energy quotient); S11c
   231	  (`steps/S11c_PARTIAL_CLOSEOUT.md`, conditional uniform result) · **target** S8 · **status** OPEN
   232	- **requirement** — the quadratic brane Lagrangian's **stiffness functional**, as delivered by the
   233	  substructure rather than chosen.
   234	- **on failure** — this is the one input the light sector's central claim is **measurably** sensitive to,
   235	  and S9's own controls prove it:
   236	  - the gradient-elastic form `−½ μ_R Σ(∂_i u_j)²` — an ordinary elastic solid — carries **the same two
   237	    transverse modes at the same `c² = μ_R/ρ_br`**, and *also* propagates the longitudinal;
   238	  - the divergence-only form makes the roles **swap**: the transverse pair drops to `ω² = 0` and the
   239	    longitudinal propagates.
   240	  ⇒ ⭐ **What curl-only buys is the ABSENCE of a propagating longitudinal mode, ⛔ not the presence of the
   241	  transverse ones.** If S8 delivers any other form, S9's and S10's transverse results may survive while
   242	  the no-longitudinal claim — Maxwell's third demand, and the whole reason the sector exists — does not.
   243	- **note** — S9 classifies the curl-only form as **postulated (structural)**, justified by *"forced by
   244	  Maxwell's no longitudinal mode"*. ⚠ That is a justification from the target, ⛔ not from the substrate.
   245	  Defect `B2` closes the route from a polar substructure `P` to `μ_R` — ⚠ **the whole quantity, not only
   246	  its magnitude**, as `DEFECT_REGISTER.md#B2` scopes it. ⛔ But it closes **one route**, and it says
   247	  nothing about the **form**, so this requirement is live.
   248	- **pass-2 consumers** — S11's unchanged transverse root rests on adding the trace invariant. As a
   249	  first-route source only, its record says: *"The form control's `4/3` independently reproduces the
   250	  corpus's rejected-Cauchy-branch coefficient"* (`steps/S11_stray_longitudinal.md:126–127`); the record
   251	  names no prior-art objection to that branch.
   252	  S11b-B needs the complete quadratic stored energy on `u, θ, e_W` under its whole stated symmetry group
   253	  (`directives/S11bB_SHARED_PHYSICS.md:353–362`). The unified record's ten-dimensional quotient is likewise
   254	  defined under its whole stated symmetry group (`directives/S11b_SHARED_PHYSICS.md:280–288`) and has different
   255	  valid representatives; no individual representative coefficient is the required object. Its
   256	  constrained breathing stiffness is `K₀ = B_ρ⁽³⁾ − 2CW₀ + k_W W₀²`. The stability claim is only the
   257	  `k = 0`, impermeable, zero-reciprocal-traction slice. On that slice B hands forward `K₀ > 0` as
   258	  *"an **empirical** constraint — the brane is here"* (`steps/S11bB_interface_assembly.md:192`). Failure
   259	  to supply this energy leaves the
   260	  decoupling and slice stability conditional; no sign or magnitude of an individual scalar coefficient
   261	  is inferred. S11c's uniform result inherits that qualification, not a nonuniform confinement law.
   262	
   263	### R-S8-02 — the in-plane and out-of-plane sectors must decouple at quadratic order
   264	
   265	- **source** S10; S11 (`steps/S11_stray_longitudinal.md:278–280`, in-plane, frozen-wall-width spectrum)
   266	  · **target** S8 · **status** OPEN
   267	- **requirement** — the quadratic operator on the brane's **full** displacement, in-plane `u` **and**
   268	  out-of-plane `h` together, and whether it is block-diagonal in that split.
   269	- **on failure** — ⛔ **S10's headline number changes** — but ⚠ **not by the mechanism an earlier draft of
   270	  this entry named, and a leg corrected it.** Mixing is **not** the failure mode: under the rotation group
   271	  that acts on the brane, `h` and `u_L` are both **scalars** and `u_T` is a **vector**, so `h` can mix with
   272	  `u_L` freely and the transverse count stays `D − 1` regardless. ⭐ **The condition that gives 3 is `h`
   273	  being DEGENERATE WITH THE TRANSVERSE PAIR** — the brane's elasticity being isotropic in `D+1` rather
   274	  than in `D`. That is what *"belonging to the same elastic sector"* has to mean, and it is what S8 must
   275	  settle.
   276	- **note** — ⭐ **This is the most load-bearing open item in the sector**, and it was flagged inside S10 at
   277	  the moment the identification was made rather than found later: `h ≠ u_L`, stated in three places in the
   278	  corpus, **user-confirmed as the picture** — but ⛔ confirmed as a picture, not computed. `V3_STEP_PLAN`
   279	  puts *"transverse and longitudinal sectors; the reduced `h`/`u_L` operator"* in S8, so the object is
   280	  already scheduled.
   281	  ⚠⚠ **And v3 is not starting from nothing — a leg found the computation already exists in v2 and this
   282	  entry did not cite it.** `research/pde_ledger_v2/paper/stages/stage_030.tex:99-125` builds the coupled
   283	  `(u_L, h)` scalar block with stiffness `K = [[B_eff, C_hu], [C_hu, K_h]]` — explicitly **not**
   284	  block-diagonal, with `C_hu` a registered free-unreduced parameter — and `R79` records that the mixed
   285	  poles are cone-coincident only when `C_hu = 0`. ⇒ ⛔ *"nothing computes that"* is true of **v3** and
   286	  misleading about the **corpus**; S8 should start from `stage_030`, ⛔ not from scratch.
   287	  ⭐ Note this is the `(u_L, h)` block — the scalar sector — so it bears on the **longitudinal** slot and
   288	  the charge anchor, ⛔ and not directly on the transverse count, which is what the corrected on-failure
   289	  above turns on.
   290	- **pass-2 source** — S11's result here is the one its record names as the in-plane, frozen-wall-width
   291	  spectrum with `hBranon` excluded: *"⚠ Scope caveats, none touching the decoupling: (i) this is the
   292	  **in-plane, frozen-wall-width** spectrum (`WALL_WIDTH_FIELDS={}`, `hBranon` excluded,
   293	  `INTERFACE_EQUATIONS_SUPPLIED={}`) — it does not decide inhomogeneous mode conversion,"*
   294	  (`steps/S11_stray_longitudinal.md:278–280`). S11's decoupling statement, with its own limit, reads:
   295	  *"Thus the homogeneous \(D=3\) quadratic transverse and longitudinal eigenbranches have zero linear
   296	  cross-block under the selected action. This does not establish decoupling at nonlinear order, on a
   297	  nonuniform slab, at an interface or defect, or after additional allowed fields are introduced."*
   298	  (`steps/S11_stray_longitudinal.md:74–76`). This entry's object is the joint `u`/`h` operator that the
   299	  spectrum leaves out.
   300	- **excluded field** — in that caveat S11's record freezes the wall width and excludes `hBranon`
   301	  (`steps/S11_stray_longitudinal.md:279–280`); its Mathematica audit binds the token as
   302	  `"OUT_OF_PLANE_FIELD_EXCLUDED" -> hBranon`
   303	  (`mathematica/S11_stray_longitudinal_mathematica_audit.wl:1082`). If S8 does not deliver the joint
   304	  operator, S11's spectrum is not extended beyond that stated scope. The record does not say what its
   305	  roots or census would be with `hBranon` included, and this register does not infer it.
   306	
   307	### R-S8-03 — the SIGN of the physical transverse stiffness
   308	
   309	- **source** S9, S10; S11 (`steps/S11_stray_longitudinal.md`, unchanged transverse branch);
   310	  unified S11b (`steps/S11b_interface_coupling_law.md`, transverse mode); S11c
   311	  (`steps/S11c_PARTIAL_CLOSEOUT.md`, conditional uniform result;
   312	  `steps/S11c_d_profile_conditioned_scattering.md`, uniform check) · **target** S8 · **status** OPEN
   313	- **requirement** — the **sign** of the physical transverse stiffness in the quadratic brane
   314	  Lagrangian, as delivered by the substructure: `μ_R` in S9/S10's selected action, `μ_⊥` in the enlarged
   315	  S11b action. `R-S8-01` asks for the stiffness **functional**; this asks for its transverse sign.
   316	- **on failure** — ⛔ **the transverse sector is two exponentially growing modes rather than two waves,
   317	  and every mode count in S9 and S10 is unchanged.** S10's `XFORM_SIGNFLIP` control measures exactly this:
   318	  with the sign flipped, `ω² = −μ_R k²/ρ_br` and **every nullity is identical to the baseline's** — the
   319	  count cannot tell a wave from an instability. What distinguishes them is a single emitted object,
   320	  `ROOT2_Q3_SIGN`.
   321	- **note** — ⭐ **Found by a review leg reading the register against S10's own controls**, ⛔ not by the
   322	  pass that wrote the register: `R-S8-01` is entirely about form and never mentions sign or positivity,
   323	  so the pass consolidated one and dropped the other.
   324	  ⚠ **It is not closed by defect `B2`**, which shuts the route for deriving `μ_R` — its sign as well as
   325	  its magnitude — from a polar substructure `P`. Other substrate routes remain open: a step that
   326	  delivers the stiffness functional delivers its sign with it, so the retirement condition is exactly
   327	  as live as `R-S8-01`'s.
   328	- **pass-2 object** — in the enlarged S11b action this is the physical transverse stiffness `μ_⊥`:
   329	  WL-`μ_R` and SymPy-`(μ_R + μ_S/2)` are its two representatives under the unified record's invertible
   330	  coefficient map. The obligation concerns the sign of `μ_⊥`, not a basis-dependent sign of `μ_S` or
   331	  either engine's bare `μ_R`. Stable real-frequency uniform modes require `μ_⊥ ≥ 0`; decoupling alone
   332	  does not exclude growth. S11c's selected finite-current check keeps its tested-input and review scope.
   333	
   334	### R-S8-04 — what carries the brane's internal angular momentum
   335	
   336	- **source** O2 (`steps/O2_steady_brane_balance.md`, register handoff); S9, S10, S11 · **target** S8 · **status** OPEN
   337	- **requirement** — the object in the substructure that carries **internal angular momentum** on the brane,
   338	  or the couple-stress it supports. ⛔ Not a mechanism, ⛔ not a model: the object, or a statement that there
   339	  is none.
   340	- **on failure** — ⛔⛔ **the curl-only stiffness functional is not an admissible continuum mechanics.** An
   341	  energy in `(∇×u)²` alone has an **antisymmetric** Cauchy stress, and balance of angular momentum forces
   342	  the Cauchy stress to be **symmetric** unless the medium carries distributed couples or internal spin. If
   343	  the substructure supplies neither, the light sector's central form is inadmissible **regardless of its
   344	  mode content** — S9's and S10's mode counts, dimensions and speeds would all be computed from a
   345	  functional no medium can have.
   346	- **note** — ⚠ **This is the objection that sank MacCullagh's aether**, and it is the one part of that
   347	  theory the 19th century never answered: Stokes pressed it, and Kelvin's gyrostatic models were attempts to
   348	  supply exactly this object. ⇒ ⭐ **prior art is the oracle here** — it tells us the obligation is real and
   349	  that answers exist, ⛔ it tells us nothing about whether **ours** delivers one ⇒ `CLAUDE.md` rule 16.
   350	  ⭐ Known families a delivered answer might fall into, ⛔ **none assumed and none prescribed**: continua
   351	  with an independent microrotation degree of freedom (Cosserat/micropolar), and media with stored internal
   352	  angular momentum.
   353	  ⛔ **One family is the wrong one and should not be reached for:** modern *odd elasticity* buys an
   354	  antisymmetric modulus tensor by making the solid **active and non-conservative**. MacCullagh's medium is
   355	  **conservative** — it has a genuine energy functional — so a non-conservative realisation would be
``````

</details>

<a id="r12"></a>
<details>
<summary>R12 — command and literal output</summary>

Command:

``````bash
nl -ba research/pde_ledger_v3/steps/S11_stray_longitudinal.md | sed -n '1,165p;215,250p'
``````

Exit code: 0.

Literal output:

``````text
     1	# S11 · The stray longitudinal — closed homogeneous census, open interface physics
     2	
     3	**Sector 1 (light), step 3.** Walked side by side with the user, 2026-08-02.
     4	
     5	---
     6	
     7	## What it is
     8	
     9	S10 left one direction of `u` with `ω² = 0`. S11 asks what happens to that direction once the brane is
    10	allowed to resist compression, what that costs the ledger, and what the resulting mode's physical status
    11	actually is.
    12	
    13	## What it does
    14	
    15	Lifts S10's zero to a propagating longitudinal mode, enters the brane's compression modulus as a
    16	**postulated** knob with a **named retirement condition**, closes the homogeneous three-finite-mode
    17	census, and derives the simple \(k_w=0\) grazing locus. The free mode's actual leakage remains an
    18	interface and spectral-boundary problem; its possible finite-amplitude photon-helper role additionally
    19	requires nonlinear coupling and a harmonic/sideband radiation audit.
    20	
    21	⛔ **This is not a defect being repaired.** The extra mode is a characterized departure from a
    22	transverse-only Maxwell field. Historical drum-head discussions used the compressional response as
    23	motivation for a material electric sector, but the current canonical ontology does **not** identify
    24	\(u_L\) with charge: charge is the oriented throat/electric-odd boundary condition, while \(u_L\) is a
    25	separate in-plane longitudinal material mode. ⛔ `FAIL_CAUCHY_STRAY_LONGITUDINAL` is a **misnamed** token;
    26	never read the prefix as a verdict.
    27	
    28	---
    29	
    30	## The walk, move by move
    31	
    32	**Move 1 — is there a compression channel at all?** ⭐ **The identification, flagged before use and
    33	user-confirmed:** `u` is the **material displacement of the stuff whose density is `ρ_br`**. That is
    34	forced by S9's own kinetic term — `½ρ_br(∂_t u)²` is only kinetic energy if `u̇` is the velocity of
    35	material carrying inertia `ρ_br`. ⇒ Linearised continuity applies, `δρ_br = −ρ_br(∇·u)`, and since the
    36	medium has an equation of state, **compression costs energy**. ⛔ Not a modelling choice.
    37	⚠ **What would have made it wrong:** if `u` were an orientation/director field, `∇·u` would be a splay,
    38	not a density change, and the channel would need separate argument.
    39	
    40	**Move 2 — what a quadratic stiffness on the brane can be.** Decompose `∂_i u_j` into trace,
    41	symmetric-traceless and antisymmetric parts. ⭐ **S10 did not forget a term — it kept ONE invariant of
    42	several.** Curl-only stiffness charges for twist and nothing else, which is why the longitudinal came out
    43	free rather than forbidden: a longitudinal plane wave has a symmetric gradient
    44	and hence no curl. That gradient generally contains both trace and symmetric-traceless
    45	parts; it is not pure trace. The selected curl-only action charges neither of those
    46	parts. (Clarified during the S11 homogeneous Lean fidelity pass, 2026-09-16 UTC.)
    47	
    48	⛔⛔ **CORRECTION — the orchestrator's count in this move was WRONG, and both engines caught it.** The
    49	walk asserted *"exactly three coefficients, for every `D ≥ 2`."* ⭐ **False under proper rotations.**
    50	
    51	```
    52	N_SO = {D=2 → 4, D=3 → 3, D=4 → 4, D=5 → 3}        N_O = 3 for every D
    53	```
    54	
    55	The extras are **reflection-odd**: `(tr G)·ε^{ij}G_{ij}` at `D=2`, `ε_{ijkl}G_{ij}G_{kl}` at `D=4`.
    56	⚠ **And the walk's METHOD was the error, not its arithmetic** — *"decompose into pieces that do not mix
    57	and count them"* misses **cross-pairings between isomorphic summands**, which is exactly what `D=2` and
    58	`D=4` have (at `D=2` the trace and `Λ²` are both `SO(2)` scalars; at `D=4`, `Λ²` splits into self-dual ⊕
    59	anti-self-dual). ⭐ Five independent derivations agree on the corrected counts.
    60	
    61	⭐⭐ **The physical content of the extras, and it sharpens the step:** at `D=2` the extra invariant is
    62	**NOT a total derivative** — its Euler–Lagrange operator is non-zero — so an `SO(2)`-invariant Lagrangian
    63	**can mix the longitudinal and transverse sectors**. At `D=4` it **is** a total derivative and adds no
    64	operator structure. ⇒ ⭐ **Sector separation is a `D=3` fact, ⛔ not a structural one.**
    65	
    66	**Move 3 — the spectrum.** The trace term enters strictly as `B_comp k(k·a)`:
    67	
    68	```
    69	ρ_br ω² a = μ_R[k²a − k(k·a)] + B_comp k(k·a)
    70	```
    71	
    72	⇒ transverse (`k·a = 0`): the new term vanishes **identically** in the selected homogeneous \(D=3\)
    73	quadratic action. Longitudinal (`a ∥ k`): the `μ_R` bracket vanishes and the zero is not perturbed but
    74	**replaced**. Thus the homogeneous \(D=3\) quadratic transverse and longitudinal eigenbranches have zero
    75	linear cross-block under the selected action. This does not establish decoupling at nonlinear order, on
    76	a nonuniform slab, at an interface or defect, or after additional allowed fields are introduced.
    77	
    78	**Move 4 — where `B_comp` comes from.** ⛔ **The orchestrator's first mechanism was wrong** — it had
    79	brane material "squeezing out into the bulk," which smuggles in a phase transition (bulk *is* the
    80	disordered phase) and the user rejected it: *"How can the brane become unordered just by compression?"*
    81	⭐ **The correct mechanism: the brane THICKENS.** In-plane compression raises the areal density, which
    82	can be absorbed as a higher 4D density (charged by the EOS) **or as a wider wall** (charged by the wall
    83	potential). ⛔ **Neither disorders anything** — the order parameter never changes; the sheet bulges into
    84	`±w`. ⇒ Two channels accommodating additive shares of one strain, i.e. **springs in series**, so
    85	compliances add and `B_comp` is **softer than either channel alone**.
    86	
    87	**Move 5 — the flat direction.** `S7` records that the double-well **selects no slab width**. ⇒ If that
    88	flatness survived, `B_wall = 0`, hence `B_comp = 0`, and S10's zero would never lift. ⭐ **It is lifted by
    89	GRADIENTS**: a wave modulates the width, tilting the interfaces and stretching them at cost
    90	`∝ σ_wall|∇W|²`. ⇒ flat at `k=0`, stiff as `k²`, predicting a **flexural crossover** — `ω ∝ k²` at long
    91	wavelength. ⚠ **S11's cone below is computed with the wall width FROZEN**; unfreezing it should soften
    92	the mode. ⛔ If it does not, **move 5 was wrong** — say so rather than reconciling.
    93	⭐ `σ_wall` is **S6-derived**, so this predicts `B_comp` may cost no knob at all → `DEFECT_REGISTER#c12`.
    94	
    95	**Move 6 — the bulk.** ⛔⛔ **THE ORCHESTRATOR OVERSTATED THIS MOVE TWICE, and both corrections stuck.**
    96	
    97	Phase matching gives `k_w² = k²(c_L²/c_s0² − 1)`, and the walk read that as *"`c_L > c_s0` ⇒ it
    98	radiates; `c_L < c_s0` ⇒ bound."* ⭐ **Neither direction is established.** Phase matching is **kinematic
    99	only** — it determines whether a propagating bulk channel **exists**, ⛔ not whether the mode uses it,
   100	nor whether an evanescent solution forms a **bound eigenmode**. Both require an interface coupling law
   101	that is absent. ⚠ A third case was also omitted: `k_w = 0`, grazing.
   102	
   103	---
   104	
   105	## ⭐⭐ The computed result — two engines, written independently
   106	
   107	| | value |
   108	|---|---|
   109	| dynamical matrix | `M_ij = μ_R(k·k)δ_ij + (B_comp − μ_R)k_i k_j` |
   110	| transverse | `ω² = (μ_R/ρ_br)k²`, nullity `D−1`, kernel ⊥ `k` ⭐ **unchanged from S9/S10** |
   111	| longitudinal | `ω² = (B_comp/ρ_br)k²`, nullity `1`, kernel ∥ `k` |
   112	| cross-sector | ⊥ root independent of `B_comp`; ∥ root independent of `μ_R` — **computed residuals, both zero** |
   113	| finite \(D=3\) census | characteristic polynomial degree 3 in `ω²`, leading coefficient `−ρ_br³ ≠ 0`; exactly three finite roots counted with multiplicity and no hidden fourth mode |
   114	| degeneracy | **exactly** `B_comp = μ_R` |
   115	| simple longitudinal/bulk threshold | `k_w² = k²(c_L²/c_s0² − 1)`; `KW_ZERO_LOCUS` is `B_comp = ρ_br c_s0²`, equivalently `c_L = c_s0`; kinematic only |
   116	| dimensions | `[ρ_br]=(−D,0,1)`, `[μ_R]=[B_comp]=(2−D,−2,1)`, both speeds `(1,−1,0)` |
   117	
   118	⭐ Every value matches predictions committed **before either script existed**
   119	(`steps/S11_PREREGISTERED_PREDICTION.md`, `67d919bd`) — except the move-2 invariant count, where the
   120	pre-registration was **wrong and the engines were right**.
   121	
   122	**Controls.** ⭐ **FORM** — replacing the trace invariant with the symmetric-traceless one moves **both**
   123	roots (`ω²_⊥ = (μ_R+μ_br)k²/ρ_br`, `ω²_∥ = 4μ_br k²/(3ρ_br)`) and the transverse root **acquires**
   124	`μ_br`: ⇒ the trace invariant is the **unique** one that lifts the longitudinal while leaving light
   125	untouched. ⛔ **COEFFICIENT** — rescaling `B_comp` moves only the parallel root and cannot test that
   126	shape claim. ⚠ The form control's `4/3` independently reproduces the corpus's rejected-Cauchy-branch
   127	coefficient, from an engine blind to that document.
   128	
   129	## What's new
   130	
   131	| item | class | why |
   132	|---|---|---|
   133	| `Q.brane.B_comp` | ⭐ **postulated**, with a **named retirement condition** | user's call: postulate now, retire visibly at S6. Knob count is an **upper bound that can only improve**. ⇒ `DEFECT_REGISTER#c12` |
   134	| `Q.brane.c_L` | **derived** | `R5`, from `B_comp` and `ρ_br` |
   135	| the mode itself | **derived** | nullity of the dynamical matrix |
   136	
   137	⚠ **Naming, and it is a live hazard:** `B_comp` is a **brane** compression modulus. ⛔ It is **not**
   138	`K_br` (a rejected elastic branch's bulk modulus), **not** `B_eff` (`= ρ_B0²/χ_c`), and not the bulk
   139	medium's modulus. Registered with `aliases: []` and no borrowed loci.
   140	
   141	**Registry:** ambient `10 → 12`, residue `6 → 7` (`{ħ, m, K, ρ0, ρ_br, μ_R, B_comp}`). Discrete payload
   142	`{n_eos=5, D_brane=3}` unchanged and off the continuous axis.
   143	
   144	---
   145	
   146	## ⭐⭐⭐ Departure — closed homogeneous census, open interface and nonlinear roles
   147	
   148	The audit now separates four questions that earlier prose compressed into one:
   149	
   150	| question | present owner |
   151	|---|---|
   152	| where the free homogeneous longitudinal root lies relative to simple bulk sound | S11 kinematic `KW_ZERO_LOCUS` |
   153	| whether that mode is bound, a threshold state, a bound state in the continuum, or a resonance | S11b interface and full-slab spectrum |
   154	| whether matter directly observes the second cone | derived matter/interface coupling |
   155	| whether the mode participates in a finite-amplitude photon helper | nonlinear intensity coupling plus harmonic and sideband response |
   156	
   157	At homogeneous linear order, the mode's bulk leakage and direct observability require the brane–bulk
   158	interface law. Its possible finite-amplitude photon-helper role additionally requires nonlinear intensity
   159	coupling, the complete slab spectrum, and a harmonic and sideband radiation audit.
   160	
   161	⛔ **What S11 does NOT deliver**, stated so it cannot be over-read: bound-versus-resonant classification ·
   162	nonzero interface overlap · an actual leakage rate · unconditional light confinement · observability or
   163	unobservability of the longitudinal branch · a localized or nonradiative nonlinear photon helper.
   164	
   165	## ⭐ What this settles about the light sector as a whole
   215	residual that records it, ablation-verified against four distinct corruptions; `R4` covered too.
   216	⇒ `DEFECT_REGISTER#f-r5`.
   217	
   218	**Gates:** acceptance `MATCH` (`12→7, 12→7, 12→6`) · dimensional homogeneity `PASS`, `R5` HOMOGENEOUS
   219	`[1,−1,0]`, 0 undetermined · able-to-fail `PASS` · 11 tests · both engine verdicts `PASS`.
   220	
   221	### ⚠ Known limits — recorded, ⛔ not fixed
   222	
   223	- **Two of six dimensional checks in the `.wl`** remain definitional identities `(X−2Y)+2Y == X` and
   224	  ⛔ **cannot fail**; only the `ρ_br` one got an independent route. ⭐ The guard as a whole **is** able to
   225	  fail — four independent ablations all caught, and a reviewer could not construct a single-point
   226	  corruption that changed a printed dimension and still passed — but those two entries ⛔ **must not be
   227	  counted as coverage**.
   228	- **`ASSERTION_12` in the `.py`** is `c ≔ √(X)` then asserting `c² − X = 0`, identically zero by
   229	  construction. Superseded in substance by the new registry control.
   230	- **The `.wl` disclaims scope on the transverse tag but not the parallel one**, though both come from the
   231	  identical substitution with the same absent operator. ⛔ Fixed **here**, in the record, rather than by a
   232	  third repair round on pinned output lines: **the scope limit applies to both channels.**
   233	- ⛔ **The SymPy script was repaired while UNTRACKED**, so no committed pre-repair baseline exists. A
   234	  reviewer's `/tmp` copy and a full hand re-derivation of every load-bearing value substituted for one.
   235	  ⚠ A straight violation of *commit before anything destructive*.
   236	- **Disclosure, volunteered by a review leg:** a `grep` for `K_br` incidentally returned two lines of
   237	  `V3_STEP_PLAN.md`, which was on its do-not-read list. The file was not opened and every reported result
   238	  was derived beforehand.
   239	
   240	---
   241	
   242	## ⭐⭐ Strata-audit closure (2026-08-19) — the conclusion is independent of the census; the census's one physical family belongs to S11b
   243	
   244	The exhaustive rank-drop / strata audit (spec §5 Q8a/Q8b — every `ρ×ρ` minor of each `M_r`, three
   245	solve-variable sets, the degenerate strata) was subsequently run through **certified** census instruments
   246	(`reduction/s11_*`; campaign closed `cbc49029`, four repair rounds, both engines, two legs + orchestrator
   247	each). Measured result (`~/.s11_build/census_build4/`): the two engines **under-decide 917 sub-cases** and
   248	carry finding-level gaps — 171 spurious branches, 72 omitted solve **records** (multiple missing
   249	memberships each), 104 membership witness failures — plus the 7 registered defects (`DEFECT_REGISTER`
   250	entries 4–7 + the obligation-4 instrument). ⚠ **Not all of these are "hard CAS questions":** the register
``````

</details>

<a id="r16"></a>
<details>
<summary>R16 — command and literal output</summary>

Command:

``````bash
rg -n -i 'polari[sz]|helicity|handedness|birefringence|two transverse|out.of.plane|spin problem' research/pde_ledger_v3/DEFECT_REGISTER.md
``````

Exit code: 0.

Literal output:

``````text
156:| **B2** | **Gate L returned a no-go — and the no-go is of a NARROWER thing than this row used to claim.** `FAIL_COUPLE_STRESS_NOGO`, the gate the medium survey had flagged as *"Highest risk; the most likely no-go"*, is a **route-failure for *deriving* the shear modulus `μ_R` from a polar substructure `P`** — ⛔ **not** a finding that one medium cannot carry a longitudinal and a transverse mode. The source says so outright: *"The no-go rules out only *deriving* `μ_R` from `P`; light stands on the bare postulated modulus (`pathA_36`/stage003 gets photons `P`-free)"*. ⚠ Its provenance is `CONDITIONAL_ON(both)` — *"conditional on the imposed axis and the postulated MacCullagh package"* — and the gauntlet was ⛔ **not** hardwired to fail: an able-to-**PASS** fixture `FREE_LIGHT_OK_CONDITIONAL` exists, *"the able-to-**PASS** tooth, proving the verdict machinery is NOT hardwired to fail"*. ⭐ **Countervailing computed fact:** the `μ_br > 0` branch carries *"two transverse modes `μ_br k²` **but also** a longitudinal mode `(K_br+4μ_br/3)k²`"* — simultaneously and consistently ⇒ the obstruction is to **suppressing** the longitudinal one (`FAIL_CAUCHY_STRAY_LONGITUDINAL`), ⛔ **never** to carrying both. ⛔ **The no-go itself stands undiluted**, and the field `P` is retired on a confirmed structural instability (`INSTABILITY_CONFIRMED_STRUCTURAL`, Decision 16) — but as *"retired-but-NOT-foreclosed (re-entry needs a NEW T0 freeze)"*. ⚠ The *"SUPERSEDED at the brane-existence level"* line is about the **GNLS-polar-smectic gate program**, ⛔ not about light | `research/pde_ledger_v2/notes/stage030_pathA35_gateL_source_map.md:240-244`, `:158`, `:16`; `software/stage1_solver/reports/pathA_35_gateL_light.md:11-13`; `software/stage1_solver/decisions/15_em_medium_native_physical_picture.md:266-272`; `software/stage1_solver/decisions/16_retire_brane_polar_field.md:12`; `research/pde_ledger_v2/notes/parameter_register.md:362`; `docs/conceptual_history.md:362` (quote), `:340`; `docs/medium_requirements_and_prior_art.md:172-177` (gate definition) | **FALSIFIED** | {#B2}
177:| **C10** | ⭐⭐ **The Spin Problem — a structural tension INSIDE the gravity sector, and it is in v3's scope.** The **inertial** sector (correct 1PN precession) constrains the throat to a **compact, "stubby"** geometry with a small aspect ratio. The **spin** sector needs the gravitomagnetic potential to fall off as a **dipole**, and a compact source gives the wrong scaling — a simple vortex yields a *"gravitomagnetic **monopole**"*, called *"physically inadmissible"*. Verdict in the note: *"You cannot get frame dragging from a compact 4D bubble; you **need** the tail."* ⇒ **Inertia wants compact; spin wants extended.** ⚠ The proposed fix — a composite *"Ion-Vortex Complex"* (stubby head + infinite vortex-filament tail) — is **uncited, unaudited, and referenced by nothing else in the corpus**; it also re-introduces **quantized circulation**, which the charge sector explicitly disclaims (*"not an additive winding"*). ⛔ Meanwhile `conceptual_foundation.md:589` still lists spin as *"**not yet placed in the picture**"* — so the corpus holds "unsolved" and "solved" simultaneously | `research/4d_1pn_bridge/notes/tadpole.md:1-19,:150-174`; `docs/conceptual_foundation.md:589` | **OPEN** |
288:*transverse-traceless* / *tensor mode* / *GW polarization* returns **nothing**.
298:the locus above. ⭐ Note also that GR's gravitational waves are **transverse with two polarizations** —
``````

</details>

<a id="r17"></a>
<details>
<summary>R17 — command and literal output</summary>

Command:

``````bash
nl -ba research/pde_ledger_v3/steps/O2_steady_brane_balance.md | sed -n '102,134p'
``````

Exit code: 0.

Literal output:

``````text
   102	| Premise | Adopted content and representation retained |
   103	|---|---|
   104	| 1 — material reference | Elastic in the optical shear regime, relaxing under steady load. Reference/strain evolution remains general and OPEN. Relaxation power is explicit, its sign OPEN, with any net power naming its supplier and budget. Optical consequences are deferred. |
   105	| 2 — drive | Dynamical order-conversion drain, with no separate external body force. The spec chooses the representation with `F_drive` absent as a separate entry; source/boundary provenance enters material, load and exchange accounting. No drain-force aggregate is added beside them. The link to orbital `GM` stays OPEN. |
   106	| 3 — exchange momentum | Converted material carries the local brane material velocity `V`; this closes transported momentum only. In the signed outward-loss convention, in-plane carry is `j_n V^i`. Bulk-direction carry remains an OPEN `𝒥_map`/normal-response reduction at native transfers. Additional non-variational partners/reactions stay OPEN with S12. This supplies no momentum-density identification or carried-total-energy formula. |
   107	| 4 — bulk loading | Retain the postulated shear-free scalar bulk, native face-normal mechanical loading and no independent tangential bulk stress. The signed amplitude, native geometry/map and full support partition stay OPEN. A tilted face-normal load can have in-plane coordinate projections. External support and carried exchange remain distinct. |
   108	
   109	`V`, `ρ_br`, `μ_⊥`, `ξ_w`, `h`, `δ`, `j_n`, bulk state, material responses, history/reference,
   110	loads and boundary data remain live with spatial derivatives. Radial profiles impose no constitutive
   111	isotropy, parity, stress symmetry, absent couple/chiral content, derivative cutoff or finite history
   112	state. No sharp-sheet or finite-slab material reduction, stress measure, relaxation family, passive
   113	sign or physical support is selected. Sources: spec §§1–6; contract §§1, 3–8 (M1).
   114	
   115	Supplied counting is `ε=GM/(c₀²r)`, `δ=O(ε)`, `(∂ξ_w)²=O(ε)`, `V/c₀=O(ε^{1/2})`,
   116	`(V/c₀)²=O(ε)` and the optical monomial box `0≤a≤1, 0≤b≤2, 0≤c≤1` in
   117	`δ^a(V/c₀)^b((∂ξ_w)²)^c`, with nonnegative integer indices. Only the stiffness/density ratio
   118	inherits the speed-change grade. Individual density/modulus, stress/inertia/normal response,
   119	exchange/source/load, relaxation/power, holder/embedding and derivative grades remain OPEN. O2 is
   120	untruncated; this optical box removes no mechanical or energy term. Bulk first order in `f` is
   121	separate, with no supplied `f`–`ε` relation. The mass law carries its recorded relative-`O(ε)`
   122	qualification for claims transferring `j_n` or `ρ_br` to induced measure; no induced-measure mass
   123	law is supplied or derived. Sources: spec §7; contract §9 (M1).
   124	
   125	Transfer limits are exact: no strong field, mouth/interior, moving/rotating mass, time-varying drain,
   126	swirl/angular profile, direction/polarization-dependent optical stiffness, coupled thickness/bulk
   127	branch, extra mixed `ωk`, or subprincipal/polarization transport result. Historical S11b rest-bulk,
   128	uniform quadratic/breathing-slice inputs and S11c supplied-profile, first-shape/background-jet,
   129	frozen-current/uniform-fold content retain those domains and unresolved composition. Static L3–L6
   130	postulated-parent, fixed-`ℓ`, constant-coefficient, frozen-sleeve/source-free/held-mouth relations
   131	are not live governing equations. S9's no-dissipation/frequency-independent limits are revisited by
   132	premise 1; its in-plane-flow and isotropic-strain freezes were lifted in v9. The retained/newly
   133	relaxed-reference comparison remains **EXPLORATORY / PAUSED**, with its earlier qualified verdict
   134	**COHERENT CONDITIONAL COMPARISON**; neither comparison law is adopted. Sources: spec §§7–8; contract
``````

</details>

<a id="r18"></a>
<details>
<summary>R18 — command and literal output</summary>

Command:

``````bash
nl -ba research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | sed -n '315,350p'
``````

Exit code: 0.

Literal output:

``````text
   315	  S11b action. `R-S8-01` asks for the stiffness **functional**; this asks for its transverse sign.
   316	- **on failure** — ⛔ **the transverse sector is two exponentially growing modes rather than two waves,
   317	  and every mode count in S9 and S10 is unchanged.** S10's `XFORM_SIGNFLIP` control measures exactly this:
   318	  with the sign flipped, `ω² = −μ_R k²/ρ_br` and **every nullity is identical to the baseline's** — the
   319	  count cannot tell a wave from an instability. What distinguishes them is a single emitted object,
   320	  `ROOT2_Q3_SIGN`.
   321	- **note** — ⭐ **Found by a review leg reading the register against S10's own controls**, ⛔ not by the
   322	  pass that wrote the register: `R-S8-01` is entirely about form and never mentions sign or positivity,
   323	  so the pass consolidated one and dropped the other.
   324	  ⚠ **It is not closed by defect `B2`**, which shuts the route for deriving `μ_R` — its sign as well as
   325	  its magnitude — from a polar substructure `P`. Other substrate routes remain open: a step that
   326	  delivers the stiffness functional delivers its sign with it, so the retirement condition is exactly
   327	  as live as `R-S8-01`'s.
   328	- **pass-2 object** — in the enlarged S11b action this is the physical transverse stiffness `μ_⊥`:
   329	  WL-`μ_R` and SymPy-`(μ_R + μ_S/2)` are its two representatives under the unified record's invertible
   330	  coefficient map. The obligation concerns the sign of `μ_⊥`, not a basis-dependent sign of `μ_S` or
   331	  either engine's bare `μ_R`. Stable real-frequency uniform modes require `μ_⊥ ≥ 0`; decoupling alone
   332	  does not exclude growth. S11c's selected finite-current check keeps its tested-input and review scope.
   333	
   334	### R-S8-04 — what carries the brane's internal angular momentum
   335	
   336	- **source** O2 (`steps/O2_steady_brane_balance.md`, register handoff); S9, S10, S11 · **target** S8 · **status** OPEN
   337	- **requirement** — the object in the substructure that carries **internal angular momentum** on the brane,
   338	  or the couple-stress it supports. ⛔ Not a mechanism, ⛔ not a model: the object, or a statement that there
   339	  is none.
   340	- **on failure** — ⛔⛔ **the curl-only stiffness functional is not an admissible continuum mechanics.** An
   341	  energy in `(∇×u)²` alone has an **antisymmetric** Cauchy stress, and balance of angular momentum forces
   342	  the Cauchy stress to be **symmetric** unless the medium carries distributed couples or internal spin. If
   343	  the substructure supplies neither, the light sector's central form is inadmissible **regardless of its
   344	  mode content** — S9's and S10's mode counts, dimensions and speeds would all be computed from a
   345	  functional no medium can have.
   346	- **note** — ⚠ **This is the objection that sank MacCullagh's aether**, and it is the one part of that
   347	  theory the 19th century never answered: Stokes pressed it, and Kelvin's gyrostatic models were attempts to
   348	  supply exactly this object. ⇒ ⭐ **prior art is the oracle here** — it tells us the obligation is real and
   349	  that answers exist, ⛔ it tells us nothing about whether **ours** delivers one ⇒ `CLAUDE.md` rule 16.
   350	  ⭐ Known families a delivered answer might fall into, ⛔ **none assumed and none prescribed**: continua
``````

</details>

<a id="r19"></a>
<details>
<summary>R19 — command and literal output</summary>

Command:

``````bash
nl -ba research/pde_ledger_v3/V3_STEP_PLAN.md | sed -n '369,416p;1170,1204p'
``````

Exit code: 0.

Literal output:

``````text
   369	### S10 · Two transverse photons
   370	⛔⛔ **REWRITTEN 2026-08-09 — ⚠ the previous two lines here asserted what the closed record REFUTES**, and
   371	they were still live while the paper card cited this very range. ⇒ ⭐ `steps/S10_two_transverse_photons.md`
   372	is the authority; ⛔ this entry may not outrun it.
   373	
   374	⭐ **What S10 measured:** the conditional map `D ↦ D − 1` transverse null directions **for the cases it
   375	swept** (`D = 2,3,4,5`), so the `D = 3` member has two — ⛔ **conditional on the supplied action, the
   376	supplied `[u]`, and BOTH structural premises.**
   377	
   378	⛔ **What this entry used to claim, and why each is wrong:**
   379	1. ⛔ *"brane shear gives exactly two transverse polarisations."* ⚠ **The form control refutes the
   380	   attribution**: `FULLGRAD` returns the **same root and the same `D − 1` count**. ⭐ What curl-only
   381	   determines is that the **longitudinal does not propagate** — a different claim about a different root.
   382	   ⚠ And the count needs a **second** premise the slogan never mentions: **isotropic inertia**. `ANISO`
   383	   keeps `S_curl`, satisfies every other condition, and drops the **exactly-transverse** count to `D − 2`.
   384	2. ⛔ *"A genuinely earned step — say so."* ⚠ The record says of itself that it does **NOT** establish a
   385	   viable physical light sector — only a conditional mode count — with **seven OPEN substrate obligations**
   386	   ⇒ `SUBSTRATE_REQUIREMENTS.md`. ⚠ And the algebra is **standard**: MacCullagh 1839
   387	   ⇒ `docs/medium_requirements_and_prior_art.md`. ⛔ **Expected new: nothing** was right; ⛔ *"earned"* was
   388	   the wrong word for why.
   389	
   390	⚠ **`D = 3` is NOT selected here** ⇒ `R-S1-01`, target **S6**, status OPEN.
   391	
   392	### S11 · The stray-longitudinal departure
   393	Characterized, first-class. ⛔ Not a defect to fix — a recorded departure from the reference theory.
   394	
   395	⛔⛔ **The token `FAIL_CAUCHY_STRAY_LONGITUDINAL` reads as a defect. It is not.** ⚠ Record the **naming
   396	defect** explicitly: the token is spelled `FAIL_…`, and that spelling already caused a reader this
   397	session to treat a **required feature** as a defect. ⛔ **Do not rename the token** — it is committed in
   398	engines and reports — ⭐ but flag at **every point of use** that the "failure" names a **characterized
   399	departure the ledger keeps as first-class** — ⛔ not a bug to erase. ⚠ ⛔ Do **not** upgrade that into
   400	*"the charge sector requires it"*: it does not (see the DIFFERENT OBJECTS block below).
   401	
   402	⭐ **The computed result: one medium carries both mode families consistently**
   403	(`software/stage1_solver/decisions/15_em_medium_native_physical_picture.md:266-272`):
   404	
   405	> *"Spectrum of the symmetric law is `(λ−μ_br k²)²(λ−(K_br+4μ_br/3)k²)` … **μ_br>0** (ordinary/Cauchy
   406	> symmetric shear …) → two transverse modes `μ_br k²` **but also** a longitudinal mode
   407	> `(K_br+4μ_br/3)k²` = a stray **"second photon"** → `FAIL_CAUCHY_STRAY_LONGITUDINAL` at Stage 4."*
   408	
   409	⇒ ⛔ **The obstruction is to SUPPRESSING the longitudinal mode** — i.e. to getting **transverse-only**
   410	Maxwell — ⛔ **never** to carrying both. Carrying two transverse modes **and** a longitudinal one
   411	simultaneously is exactly what the `μ_br > 0` branch *does*.
   412	
   413	⭐⭐ **The mode is REAL and IRREMOVABLE, and its role was CATALYTIC — ⛔ not constitutive** (user's own
   414	history, recalled 2026-07-31). ⛔⛔ **It is not charge; it does not become charge; it supplies no part of
   415	charge.** ⚠ An earlier version of this step said it "becomes the **±w** displacement" — ⛔ that was a
   416	false identification, corrected here. The history, in the order it happened:
  1170	### S22 · ⭐⭐ HALF TWO — what only a simulation can settle
  1171	
  1172	#### ⭐⭐ THE NEAR-TERM LINEAR SIM INVENTORY — banked 2026-08-02, ⛔ do not re-derive it
  1173	
  1174	These are proposed linear tests on supplied backgrounds, not a solved nonlinear throat. S11b supplies
  1175	conditional uniform interface data; S11c's nonuniform response remains partial. Do not infer that all
  1176	five tests are cleared or all depend on a completed S11c loss factor. [S11c closeout, downstream uses](steps/S11c_PARTIAL_CLOSEOUT.md).
  1177	
  1178	| # | what it measures | needs | blocked by |
  1179	|---|---|---|---|
  1180	| 1 | ⭐⭐ **transverse → longitudinal conversion at a defect** — how much light converts to the stray mode passing a particle. ⚠ The families do **not** mix in a *homogeneous* brane (coupling `∝ k·a`); `∇μ_R ≠ 0` mixes them | a defect profile | **response/form remains PARTIAL** (S11c-c2/d); physical **magnitude** needs the interior (`R1`) |
  1181	| 2 | ⭐⭐ **does the longitudinal actually radiate into the bulk** | the interface law | **S11b** |
  1182	| 3 | ⭐ **the flexural crossover** — `ω ∝ k²` at long wavelength once the wall width is dynamical. ⛔ **Able to fail:** a clean cone means move 5 was wrong | width mode inertia + stiffness | **S11b** + **S5–S7** |
  1183	| 4 | **birefringence near a defect** — the two polarisations are degenerate *by symmetry* in a homogeneous brane; a defect splits them | a defect profile | as (1) |
  1184	| 5 | **Test 5** — set `χ_B = 0`, confirm shear propagation dies (`notes/brane_bulk_handoff.md:880`) | nothing new | ⭐ **runnable now**, designed, never run |
  1185	
  1186	⛔⛔ **TAUTOLOGY GUARD:** a sim fed `μ_R` that reports `c_γ = √(μ_R/ρ_br)` has measured **nothing** — a
  1187	linear integrator already does exactly that. ⇒ Postulated parameters are fine to *run* with; they bound
  1188	what the run may **claim**. Conversion efficiency, mode structure, stability and confinement are genuine
  1189	outputs; ⛔ **the speeds are not.**
  1190	
  1191	#### ⭐⭐ "STAYS QUANTIZED IN A PACKET" IS THREE DIFFERENT QUESTIONS — ⛔ do not merge them
  1192	
  1193	⚠ The user's ambition is *"a photon travelling through space, self-sustaining."* Three readings, **three
  1194	different blockers**, and conflating them caused a false alarm that `ħ` was on the critical path:
  1195	
  1196	| reading | needs | status |
  1197	|---|---|---|
  1198	| **a bound state** — a localized mode with a **discrete** frequency spectrum, trapped by a defect | ⭐ a linear wave equation in a **varying** medium (Sturm–Liouville) | ⭐⭐ **reachable, LINEAR.** No `ħ`, no nonlinearity |
  1199	| **a soliton** — one allowed lump, energy fixed by a size–amplitude relation | the nonlinear shear action | ⛔ **C6** — the real wall, and what the user actually wants |
  1200	| **`ħω` quanta** — energy in discrete lumps | field quantization | ⛔ **not a classical PDE sim at all**, at any effort |
  1201	
  1202	⭐⭐ **⇒ Shelving `ħ_model` does NOT block the simulation.** A classical field has continuous energy no
  1203	matter how good the code is; row three was never reachable this way. ⇒ `ħ_model` connects the model to
  1204	quantum mechanics — a **separate project** from watching the brane move. ⛔ Do not re-open it here.
``````

</details>

<a id="r21"></a>
<details>
<summary>R21 — command and literal output</summary>

Command:

``````bash
nl -ba docs/native_light_em_and_vortex_throat_interpretation.md | sed -n '2267,2300p;4193,4220p'
``````

Exit code: 0.

Literal output:

``````text
  2267	### 9.4 Helicity has several non-equivalent native uses
  2268	
  2269	For a moving localized object on the brane, a particle-like helicity candidate
  2270	would be
  2271	
  2272	\[
  2273	\lambda_{\rm th}
  2274	=
  2275	\frac{S_{\rm br}\cdot P}{|P|}.
  2276	\]
  2277	
  2278	It is undefined for a throat at rest and changes meaning if mixed
  2279	\(S_{aw}\) components participate. For a freely propagating transverse shear
  2280	mode, right- and left-circular polarization may instead furnish two
  2281	wave-helicity states, but their angular momentum density and detector response
  2282	must be calculated.
  2283	
  2284	Fluid mechanics uses the word helicity for a different quantity, such as the
  2285	brane-projected integral
  2286	
  2287	\[
  2288	\mathcal H_v
  2289	=
  2290	\int d^3x\,
  2291	v_{\rm br}\cdot
  2292	\left(\nabla\times v_{\rm br}\right).
  2293	\]
  2294	
  2295	This can measure linking or twisting of brane-projected flow. It is not the
  2296	same as the helicity of a particle representation, and it is not a complete
  2297	four-dimensional invariant. A four-dimensional vorticity is a two-form, so
  2298	the model must state which contractions, projections, and boundary terms
  2299	define any proposed generalization.
  2300	
  4193	The comparison below asks whether a familiar implication is reproduced,
  4194	separately realized by material physics, open, or outside the present scope.
  4195	It does not require the native variables to share the standard ontology.
  4196	
  4197	| Familiar implication or observation | Native status | Present reading |
  4198	|---|---|---|
  4199	| two transverse light polarizations | conditionally reproduced | two exactly-transverse directions occur on the supplied homogeneous \(D=3\), isotropic-inertia branch; the count is not unique to curl-only stiffness and ANISO weakens exact transversality |
  4200	| approximately linear light dispersion | conditionally reproduced | \(\omega_T^2=(\mu_R/\rho_{\rm br})k^2\) on the selected background |
  4201	| absence of a physical longitudinal photon | separately realized, with a departure | the selected homogeneous \(D=3\) action separates two transverse branches from a physical longitudinal sound branch; that split is not structural in arbitrary dimension or background |
  4202	| microscopic origin of the curl stiffness | open after a route exclusion | the polar-\(P\) derivation failed; \(\mu_R\) remains postulated, and any different Cosserat/micropolar completion is new action content |
  4203	| microscopic electromagnetic gauge redundancy | not reproduced and not required as ontology | observable jobs must be supplied by material constraints, topology, source dynamics, or accepted departures |
  4204	| native Gauss constraint | not obtained on the tested polar-field route | exact emergent-\(U(1)\) route excluded within that quadratic scope |
  4205	| static \(1/R^2\) electric force range | conditionally reproduced | localized \(h\) carrier and nonzero core monopole required |
  4206	| like charges repel and opposites attract | open | fixed-value branch has the target sign; other admissible mouth ensembles differ |
  4207	| additive and conserved electric charge | open | \(s=\pm1\) supplies a binary orientation, not yet a conserved additive quantity |
  4208	| universal charge magnitude | open | core amplitude and normalization remain continuous or unresolved |
  4209	| current of a moving charge | partially reproduced | rigid translation gives \(I_s=s\eta V\); complete charge identity and field overlap remain open |
  4210	| magnetic Darwin tensor and velocity order | conditionally reproduced | follows from the supplied transverse moving-source coupling |
  4211	| common Coulomb/magnetic normalization | open | requires \(r_{BA}\) and the relevant cone relation from one throat branch |
  4212	| standard \(B\)-field time reversal | current characterized departure | audited \(b_T\) candidate is time-even in the active-drain construction |
  4213	| current-loop dipole forces | open held-out test | side-by-side and coaxial loops must use one derived coupling |
  4214	| transverse radiation from acceleration | open | two polarization amplitudes and extra emission channels are not derived |
  4215	| electric--magnetic induction structure | open | no complete dynamic source/carrier system currently establishes it |
  4216	| intrinsic magnetic moment | open | circulation does not guarantee a nonzero dipole |
  4217	| gyromagnetic relation | open | both \(S_{IJ}\) and \(\mu_{\rm th}\) must be derived |
  4218	| spin-\(1\) helicity of the light quantum | outside the present classical derivation | circular shear polarization is a candidate classical precursor only |
  4219	| spin-\(\tfrac12\) particle transformation | open stronger track | framed-vortex topology is a hypothesis, not a spinor derivation |
  4220	| fermionic exchange statistics | outside the current classical model | requires quantization and two-throat configuration-space analysis |
``````

</details>

<a id="r22"></a>
<details>
<summary>R22 — command and literal output</summary>

Command:

``````bash
nl -ba docs/light_guided_photon_soliton_research_plan.md | sed -n '1,16p;246,286p;1786,1800p'
``````

Exit code: 0.

Literal output:

``````text
     1	# Guided Light and the Photon-Soliton Research Program
     2	
     3	## Status and purpose
     4	
     5	This document records the current light ontology, the guided-mode hypothesis, the proposed route to a self-bound classical photon analog, and the calculations required to test it. It supplements, rather than replaces, the canonical [ontology and closure ledger](toy_model_ontology_summary.md) and the focused [opposite-orientation throat plan](opposite_orientation_throat_coincidence.md).
     6	
     7	It is a **research and derivation plan**, not a claim that the required guided modes, nonlinear couplings, or photon-like solutions have already been obtained.
     8	
     9	The model is a classical toy analog. The medium and its constitutive laws may be postulated so that they possess the properties required by the force and light sectors. The test is not whether those properties emerge accidentally from the few equations already written. The test is whether one frozen, mathematically coherent medium can support all required sectors at once without independent retuning.
    10	
    11	The following epistemic labels apply throughout:
    12	
    13	- **Postulated:** part of the candidate medium's definition.
    14	- **Operational definition:** a quantity defined by how a brane observer would measure it.
    15	- **Target relation:** a result the selected medium must produce; it may not be inserted as an answer.
    16	- **Derived under stated assumptions:** follows from a specified action, approximation, branch, or boundary condition.
   246	A guided light displacement should be tangent:
   247	
   248	\[
   249	N_A u_T^A=0.
   250	\]
   251	
   252	For a wave with brane-tangent propagation vector \(k^A\),
   253	
   254	\[
   255	N_Ak^A=0,
   256	\]
   257	
   258	transversality also requires
   259	
   260	\[
   261	k_Au_T^A=0.
   262	\]
   263	
   264	The parent medium has four spatial directions. Removing the one normal direction and the one propagation direction leaves two independent polarization directions:
   265	
   266	\[
   267	4-1-1=2.
   268	\]
   269	
   270	This gives a geometric interpretation of the two-polarization target:
   271	
   272	> A light displacement is forbidden from pointing through the normal axis and is transverse to its direction of travel inside the brane, leaving two independent tangential polarization directions.
   273	
   274	This count remains conditional on the selected action and constraints. It must be verified by the complete guided-mode spectrum rather than accepted from counting alone.
   275	
   276	### 2.3 Transverse isotropy
   277	
   278	The most natural material design is **transversely isotropic**:
   279	
   280	- no preferred direction among \(x,y,z\) in a homogeneous brane;
   281	- materially different response along the normal axis \(w\).
   282	
   283	This permits strong normal structure without selecting a fixed compass direction inside observed space.
   284	
   285	A permanently preferred direction within the brane could produce unacceptable anisotropy. By contrast, a unique normal axis is already part of the brane ontology.
   286	
  1786	### 14.3 Polarization and symmetry criteria
  1787	
  1788	11. There are only two freely selectable photon polarizations.
  1789	12. The helper field is slaved rather than independently selectable.
  1790	13. The packet has no leading electric-odd monopole on a reflection-symmetric background.
  1791	14. Linear and circular polarization branches do not acquire unacceptable speed or energy splitting.
  1792	
  1793	### 14.4 Propagation criteria
  1794	
  1795	15. The packet propagates at or extremely near \(c_\gamma\).
  1796	16. Speed dependence on amplitude, width, or frequency is acceptably small in the intended regime.
  1797	17. The energy-momentum relation is approximately light-like.
  1798	18. There is no persistent longitudinal, bulk, electric, thickness, or conversion wake from the DC helper, carrier harmonics, sidebands, or higher nonlinear source components.
  1799	19. Any metastable leakage lifetime is sufficiently long.
  1800	
``````

</details>

<a id="r23"></a>
<details>
<summary>R23 — command and literal output</summary>

Command:

``````bash
nl -ba docs/model_map.md | sed -n '1,20p'; rg -n -i 'polari[sz]|helicity|birefringence|two transverse|out.of.plane|handedness' docs/model_map.md
``````

Exit code: 0.

Literal output:

``````text
     1	# The Model Map — one-medium analog, derivation atlas + conceptual throughline
     2	
     3	> ⚠ **CURRENT FRONT (2026-07-31): `research/pde_ledger_v3/`** — start at `NEXT_SESSION.md`.
     4	> ⛔ The dimension rewrite named below is **not** the current front. {#current-front}
     5	
     6	
     7	**What this is.** The single high-level map of the whole toy model: the conceptual picture, every earned derivation with a one-line note + pointer to the full doc/scripts, the honest ledger of what is *predicted* vs *calibrated* vs *unresolved (R1)* vs *departure*, and a glossary. Read this to hold the model in your head; read the cited sources to act on any specific number.
     8	
     9	> ⛔ **SUPERSEDED 2026-07-31 — kept as history, ⛔ not as an instruction. The current front is the banner at the top of this file (`docs/model_map.md#current-front`): v3, `research/pde_ledger_v3/NEXT_SESSION.md`.** ⏸ **CURRENT FRONT (2026-07-27) — NOT stage 045.** The ledger build is PAUSED behind the **dimension
    10	> rewrite** (all 30 dimension-bearing SymPy audit scripts onto one shared module; **6 done** — 004,
    11	> 011, 012, 013, 016, 018, all waiver-free; ▶ NEXT = stage023 `.py` (step f), then 027/021).
    12	> Read `research/pde_ledger_v2/manifests/DIMENSION_REWRITE.md`, and `STATUS.md` for the front.
    13	> Every "▶ NEXT = stage 045" below is the *ledger-build* next, correct in its own sequence but not the
    14	> current action. Also paused: **stage 044-v2** (the dynamical-Σ un-freeze), which precedes 045.
    15	>
    16	> **This is a synthesis map, not the source of truth.** Assembled 2026-07-21 from a six-agent sector fan-out over the committed repo. For any discrepancy the **cited source files are authoritative** — the sector reports under `software/stage1_solver/reports/` and `software/em_charge_attribute/`, the v2-ledger stages under `research/pde_ledger_v2/notes/stages/`, the blueprint `notes/ledger_v2_blueprint.md`, the trackers `research/pde_ledger_v2/notes/{parameter_register,midway_knob_audit}.md`, and the resume doc `research/pde_ledger_v2/notes/RESUME_ROADMAP.md`. Re-read those before trusting any claim here.
    17	>
    18	> **Framing.** This is a *toy analog*: the goal is ONE self-consistent structure (a calibrated PDE) that reproduces GR-like and EM-like far-field behavior — a working math **bridge**, not an ontology claim. A result that **breaks** the concept is welcome and first-class. Magnitudes are calibrated; **among observational predictions**, only held-out **dimensionless** structure tests the model. ⚠ That is one falsification route, ⛔ **not the whole standard** — the governing standard (which also counts internal impossibility and treats surplus over an incomplete reference theory as expected) is `research/pde_ledger_v3/CHARTER.md#falsification-standard`; read it there, ⛔ not from this line. A clean "it all works" is suspicious.
    19	
    20	---
39:| **Light** | in-plane **MacCullagh shear** of the brane (energy in `curl u`, on the brane so it leaves the bulk/throat undisturbed) → two transverse photons | brane shear `u_T`, `c_γ` |
47:**Honest status in one breath:** all four far-field sectors now exist on one shared field set. Gravity is the most complete (form + fingerprints earned, magnitudes calibrated, one held-out falsifiable departure). Light earns two transverse photons but leaves a characterized stray-longitudinal departure. Charge and magnetism earn the *structure* (mechanism, falloff, tensor form) target-blind, but the electric **sign** — and with it the magnetic sign — is genuinely **unresolved (R1)**, waiting on a sim-deferred nonlinear throat solve. The model is a coherent, calibratable single-medium analog with a real, un-papered-over honesty ledger. ⛔ Gravity's "form … earned" and its "one held-out falsifiable departure" name the `1/r²` law, the attractive sign and the stage009 `RETURN_RESIDUAL_PREDICTION` residual — all three are **CONDITIONAL pending the S14a drain bridge, NOT earned** (the fingerprints are unaffected) — `research/pde_ledger_v3/CHARTER.md#conditional-s14a`.
107:*Light = in-plane MacCullagh shear of the brane; two transverse photons. The shared transverse sector — it establishes the `u_T`/`c_γ` foundation magnetism reuses. Under the surviving-solution rule Part III is stage003 alone; the `pathA_35` couple-stress no-go (retired-`P` post-mortem) → failures-paper backlog.*
111:- ⭐ **Two polarisations is a statement that our space is 3-dimensional.** The count is `D_brane − 1`,
``````

</details>

<a id="r24"></a>
<details>
<summary>R24 — command and literal output</summary>

Command:

``````bash
nl -ba research/pde_ledger_v2/paper/stages/stage_003.tex | sed -n '20,70p'; nl -ba research/pde_ledger_v2/paper/stages/stage_030.tex | sed -n '90,132p'
``````

Exit code: 0.

Literal output:

``````text
    20	\(c_\gamma^2=\mu_R/\rho_{br}\); verdict \texttt{PASS\_TRANSVERSE\_UNDISTURBED}.
    21	\emph{Top-line departure:} \texttt{FAIL\_CAUCHY\_STRAY\_LONGITUDINAL} -- on the
    22	provenance-fixed single-medium branch the \((\nabla\!\cdot\!u)+\theta\) sector
    23	carries one stray propagating longitudinal DOF, a Dirac--Bergmann
    24	\emph{second-class} pair rather than a Maxwell first-class Gauss gauge-removal.
    25	\emph{Reachability:} the tuned Maxwell locus (\(K_\theta=C_J^2/\rho_{br}\),
    26	\(B_{\mathrm{eff}}=0\), \(m_\theta^2=0\)) is reachable and emits
    27	\texttt{C5\_RESOLVED\_MAXWELL\_BY\_TUNING} -- so the FAIL is not hardcoded -- but
    28	\texttt{BY\_TUNING}, not \texttt{WITH\_PROVENANCE}.
    29	
    30	\paragraph{Derivation.}
    31	The primitive brane Lagrangian (no pre-completed square used as input) is
    32	\[
    33	  L=\tfrac12\rho_{br}(\partial_t u)^2-\tfrac12\mu_R(\nabla\times u)^2
    34	    -\tfrac12 B(\nabla\!\cdot\!u)^2+J(\partial_t\theta)\,\delta\rho_B
    35	    +\tfrac12 K_\theta(\nabla\theta)^2 .
    36	\]
    37	\resultanchor{Transverse (earned).} For a polarization \(\perp\mathbf k\),
    38	\(\nabla\!\cdot\!u=0\) and the \(\theta\) couplings vanish, so
    39	\(L_T=\tfrac12\rho_{br}(\partial_t u_T)^2-\tfrac12\mu_R k^2u_T^2\), giving
    40	\[
    41	  \omega^2=\frac{\mu_R}{\rho_{br}}\,k^2,\qquad
    42	  c_\gamma^2=\frac{\mu_R}{\rho_{br}},
    43	\]
    44	two polarizations \(\Rightarrow\) \(\mathrm{physical\_dof}=2\), massless. The
    45	speed is able-to-fail: \(\epsilon=2\rho_{br}\) shifts it to \(\mu_R/(2\rho_{br})\)
    46	(\texttt{FAIL\_TRANSVERSE\_DISTURBED}).
    47	\resultanchor{Longitudinal (departure).} Number conservation slaves
    48	\(\delta\rho_B=-\rho_{B0}\nabla\!\cdot\!u\); integration by parts gives the
    49	Maxwell-form cross-term with the \emph{derived} sign \(C_J=-J\rho_{B0}\), and the
    50	finite conjugate-density stiffness the Cauchy modulus
    51	\(B_{\mathrm{eff}}=\rho_{B0}^2/\chi_c\). Since \(\partial_t\theta\) enters only
    52	through the first-order cross-term, the momenta are
    53	\(p_u=\rho_{br}\partial_t u_L\), \(\pi_\theta=Jk\rho_{B0}u_L\), forcing the
    54	primary constraint \(\Phi_1=\pi_\theta-Jk\rho_{B0}u_L\); preservation gives
    55	\(\Phi_2=-k(Jp_u\rho_{B0}+k\kappa_{\mathrm{phase}}\rho_{br}\theta)/\rho_{br}\), and
    56	\[
    57	  \{\Phi_1,\Phi_2\}=\frac{k^2(J^2\rho_{B0}^2+\kappa_{\mathrm{phase}}\rho_{br})}{\rho_{br}}\neq0
    58	\]
    59	\(\Rightarrow\) second-class pair (\(0\) first-class, \(2\) second-class),
    60	\(N_{\mathrm{phys}}=(4-2)/2=1\). The finite-stiffness pole
    61	\(\omega^2=k^2\kappa_{\mathrm{phase}}\rho_{B0}^2/[\chi_c(J^2\rho_{B0}^2+\kappa_{\mathrm{phase}}\rho_{br})]\)
    62	has positive residue and a bounded reduced Hamiltonian: a real stray longitudinal
    63	mode, not a ghost and not a gauge-removal.
    64	
    65	\paragraph{Physical use, fidelity, and the departure.}
    66	The brane carries light: two transverse photons at \(c_\gamma^2=\mu_R/\rho_{br}\)
    67	from the medium's own moduli (light \(=\) in-plane brane shear). Honest scope: the
    68	order-parameter phase \(\theta\) does \emph{not}, on its frozen definitions,
    69	supply a gauge-removed Maxwell scalar -- it leaves one stray second-class
    70	longitudinal mode; the obstruction is pinned to the Cauchy modulus
    90	\texttt{PASS\_REDUCED\_SPEED\_PRESERVATION} checks that the \emph{separately}
    91	reduced ratio is speed-preserving, \(c_E^2=K_h/M_h=K_4/M_4\), so an inconsistent
    92	reduction (e.g.\ \(K_h=N_0^2K_4\)) makes \(K_h/M_h=160/3\ne1\) and the tooth
    93	fires on a \emph{real} error. \emph{Dependency direction:} \(M_4\) (a postulated
    94	\texttt{ACTION} primitive) and \(c_E\) (input) are the free data;
    95	\(K_4=M_4 c_E^2\), \(M_h=N_0 M_4\), \(K_h=N_0 K_4=M_h c_E^2\) are all
    96	\textbf{derived}; the reduced normalization \(M_h^*=1\) is a \texttt{CONV} choice
    97	imposed \emph{on} \(M_4\).
    98	
    99	\paragraph{The coupled \((u_L,h)\) scalar block --- two INDEPENDENT positivity facts (EARNED).}
   100	After eliminating the brane phase \(\theta_B\) (its Schur shift absorbed into
   101	\(A_{\rm eff}=\rho_{br}+C_J^2/\kappa_{\rm phase}\)), the reduced scalar action has
   102	inertia matrix \(M=\mathrm{diag}(A_{\rm eff},M_h)\) and stiffness
   103	\(K=\bigl[\begin{smallmatrix}B_{\rm eff}&C_{hu}\\C_{hu}&K_h\end{smallmatrix}\bigr]\).
   104	The stage splits kernel health into two teeth that are \emph{genuinely
   105	independent} in both engines. (a)~\emph{Sylvester positivity}
   106	(\texttt{PASS\_STABILITY}): \(B_{\rm eff}>0\) and \(D=B_{\rm eff}K_h-C_{hu}^2>0\)
   107	(physical \([D]=M^2T^{-4}\)), confirmed by SymPy \texttt{is\_positive\_definite}
   108	/ Wolfram \texttt{PositiveDefiniteMatrixQ}+\texttt{Eigenvalues}. (b)~\emph{Positive
   109	generalized wave speeds} (\texttt{PASS\_POSITIVE\_GENERALIZED\_WAVE\_SPEEDS}): the
   110	squared characteristic speeds \(z=c^2/c_E^2\) solve \(\det(K-zM)=0\); in the
   111	explicitly-normalized star units (\(M^*=I\), \(B_{\rm eff}^*=2\), \(K_h^*=1\),
   112	\(C_{hu}^*=1/2\)),
   113	\[
   114	  \det(K^*-z^*M^*)=z^{*2}-3z^*+\tfrac74,\qquad
   115	  z^*_\pm=c^{*2}_\pm=\frac{3\pm\sqrt2}{2}>0,
   116	\]
   117	so both coupled modes propagate. The two facts are independent: \(D\) depends on
   118	\(K\) only, while the speeds depend on both \(K\) and \(M\); the (b)-tooth
   119	mutates the inertia (\(A_{\rm eff}^*:1\to2\)), which \emph{keeps} the Sylvester
   120	margin \(D^*=7/4\) intact yet \emph{changes} the roots --- driving (b) to exit~1
   121	at its own assert while (a) stays green, so the roots check is not a corollary of
   122	the margin. (The physical \(\det(K-zM)\) is dimensionful and is never equated to
   123	the numeric polynomial --- the determinant is stated only through the star form.)
   124	
   125	\paragraph{Reduced-\(h\) masslessness and conservative Hessian symmetry (EARNED).}
   126	The time-dependent quadratic operator
   127	\(Q_s(\omega,k)=\bigl[\begin{smallmatrix}A_{\rm eff}\omega^2-B_{\rm eff}k^2 & -C_{hu}k^2\\ -C_{hu}k^2 & M_h\omega^2-K_h k^2\end{smallmatrix}\bigr]\)
   128	satisfies \(Q_s(0,0)=0\) at star values --- there is \emph{no} \(k^0h^2\) term
   129	(no bare \(h\) mass), the static prerequisite against
   130	\texttt{FAIL\_PINNED\_BRANON}/\texttt{FAIL\_YUKAWA}
   131	(\texttt{PASS\_REDUCED\_H\_MASSLESSNESS}; the tooth reinstates a \(k^0h^2\)
   132	coefficient and \(Q_s(0,0)\ne0\) fires). The gradient-energy Hessian of the
``````

</details>

<a id="r26"></a>
<details>
<summary>R26 — command and literal output</summary>

Command:

``````bash
nl -ba research/4d_plasma/paper/4d_plasma.tex | sed -n '2575,2590p;4291,4323p'
``````

Exit code: 0.

Literal output:

``````text
  2575	verify discrete \(\partial_t\rho_q+\gradthree\cdot\mathbf{J}=0\) and stability of
  2576	Gauss constraints.
  2577	\end{enumerate}
  2578	
  2579	\subsubsection*{(B) Two-fluid and ideal-MHD wave benchmarks}
  2580	\begin{enumerate}[leftmargin=1.6em]
  2581	\item \textbf{Circularly polarized Alfv\'en wave.}
  2582	Initialize a nonlinear Alfv\'en wave (a standard code verification target) and
  2583	verify amplitude preservation, phase speed, and convergence.
  2584	
  2585	\item \textbf{Fast/slow magnetosonic waves.}
  2586	Initialize linear perturbations about a uniform background and verify dispersion
  2587	relations in the two-fluid solver (and in the ideal-MHD closure when invoked).
  2588	
  2589	\item \textbf{Whistler/Hall branch (optional extended-MHD check).}
  2590	If the generalized Ohm law \eqref{eq:ohm_general} is used with the Hall term retained,
  4291	\subsection{Magnetic helicity identity at fixed \texorpdfstring{\(w\)}{w}}
  4292	\label{app:helicity_identity}
  4293	
  4294	Fix \(w\) and define the brane vector potential \(\mathbf{A}(t,\mathbf{x},w)\) and scalar
  4295	potential \(\Phi(t,\mathbf{x},w)\equiv A_0(t,\mathbf{x},w)\). Define the brane magnetic field
  4296	\(\mathbf{B}\equiv \nabla_3\times\mathbf{A}\) and brane electric field
  4297	\(\mathbf{E}\equiv -\nabla_3\Phi-\partial_t\mathbf{A}\). Then the standard local helicity identity
  4298	holds pointwise in \(w\):
  4299	\begin{equation}
  4300	\partial_t\!\big(\mathbf{A}\cdot\mathbf{B}\big)
  4301	+
  4302	\nabla_3\cdot\Big(\Phi\,\mathbf{B}+\mathbf{E}\times\mathbf{A}\Big)
  4303	=
  4304	-2\,\mathbf{E}\cdot\mathbf{B}.
  4305	\label{eq:helicity_local_app}
  4306	\end{equation}
  4307	This identity is purely \(3+1\) (with \(w\) acting as a parameter), and is independent of the
  4308	details of the \(4+1\) dynamics.
  4309	
  4310	\paragraph{Gauge note.}
  4311	Magnetic helicity is gauge invariant on a periodic domain, or for fields satisfying standard
  4312	boundary conditions (e.g.\ \(\mathbf{B}\cdot \mathbf{n}=0\) and fixed tangential \(\mathbf{A}\)
  4313	on \(\partial\Omega\)), or when formulated as \emph{relative helicity} with respect to a
  4314	reference field. In what follows, the boundary flux term is kept explicit.
  4315	
  4316	\subsection{Projected helicity budget and ``helicity leakage'' into transverse structure}
  4317	\label{app:proj_helicity_budget}
  4318	
  4319	Project \eqref{eq:helicity_local_app} with \(W(w)\) and integrate over \(w\). Define the
  4320	\emph{projected helicity density} and projected helicity flux:
  4321	\begin{equation}
  4322	\overline{h}(\mathbf{x})
  4323	\equiv
``````

</details>

<a id="r27"></a>
<details>
<summary>R27 — command and literal output</summary>

Command:

``````bash
nl -ba research/pde_ledger_v3/steps/LIGHT_LEAKAGE_SCOPING.md | sed -n '34,37p;148,153p'
``````

Exit code: 0.

Literal output:

``````text
    34	| 10 | Leaky modes/tunnelling — outgoing continuation or evanescent-continuum coupling (SL). | **PRESENT ANALOG:** acoustic legs/impedance/Fredholm loci, not direct bulk shear propagation (c1:35–44; S9:180–187). | Kernel/pressure/relative-flux **“ESTABLISHED”**; other c1 objects **“UNDECIDED”** (c1:100–110); physical loss unestablished (PC:3,9). | **YF/LP:** measured aggregate attenuation/cavity photon lifetime constrain optical-channel removal in their exposures; tunnelling-only component **gap**. | Propagating acoustic bulk; reactive near field not automatically loss. Direct bulk shear is O4 (A:41–43; c1:43–44; SR:144–153). |
    35	| 11 | Guided-mode coupling — perturbation connects supported modes (SL; Marcuse 1969 oracle). | **PRESENT ANALOG:** close-then-extract transverse↔`{θ,e_W,u_L}`, not two optical transverse modes (c2:59–70). | Kernel present (c2:63–70); repairs **“CONDITIONAL — repaired equations, scoped review only”** (c2:42); full composition/direct mixed entry **“UNRESOLVED”** (c2:51,55); “F/G remain withdrawn interpretations” (c2:53). | **YF:** aggregate optical removal, not redistribution among surviving transverse modes; scalar-conversion-specific measurement **gap**. Marcuse exchange example is theory. | **Another brane branch** `θ,e_W,u_L`, bulk through closure; purely transverse redistribution survives (c2:63–70; c1:43–44; PC:3). |
    36	| 12 | Polarization coupling/birefringence — symmetry breaking splits/mixes transverse polarizations (K/SL). | **PRESENT ANALOG, named future item:** S22 conversion/birefringence and defect inputs. d's polarization-dependent forcing is not identified with polarization coupling (PC:31; d:58). | **“OPEN”** (PC:31); no completed birefringence law. | **PB:** measured magnetic-birefringence upper limit on phase/polarization response; defect-observable map **gap**. | Same transverse branch **survives** for polarization redistribution; receiver/power law for the separately named conversion item **gap** (PC:31; PC:3). |
    37	| 13 | Lossy core/cladding overlap — absorbing/radiating second region (SL). | **ABSENT** as the records state the uniform transverse mode: “The transverse mode is completely decoupled on a uniform background. The coupling is **identically zero**” (B:74–75; B:191). The bulk's shear-freeness is postulated, not established; failure permits escape (S9:63–68,180–187; SR:151). No absorber follows from either. | B:191 **“unconditional on a uniform background”**; PC:7 and d:7 mark uniform decoupling **“CONDITIONAL”** (quoted side by side at E4, §2.1). R-S1-02 **“OPEN”** (SR:146); premise **“LIVE”**, **“POSTULATED”** (S9:183–185). | **YF:** measured aggregate attenuation per length; overlap-specific decomposition **gap**. | No second transverse region under the postulate; **bulk** if it fails (O4; SR:151). Absorptive receiver **gap**, not frequency exchange. |
   148	| V0 | B:61–63 says drain `v₀` **“bounded by cosmology”**. The plan names an observable for the DC leak: “it ties the DC leak rate to an **observable** (the expansion rate), which is a future calibration hook — ⛔ not claimed, not derived, and out of scope here” (PLAN:572–573). The plan calls the DC part of `S_leak` the dark-energy mechanism (PLAN:556–564); B:47 calls `v₀` “the dark-energy drain”. No record line read here states that the DC leak rate and `v₀` are one quantity, and no numerical bound is sourced: **gap**. SR:201–203 distinguishes normal drain from brane-rest `v₀=0`. | Drive/reservoir budget; numerical receiver map **gap**. | Cosmological time/column **gap**. |
   149	| FD / CMB blackbody | Fixsen et al., [ApJ 473, 576 (1996), primary abstract](https://arxiv.org/abs/astro-ph/9605054): FIRAS RMS deviations **less than 50 ppm of CMB peak**; **`\|y\|<15×10⁻⁶`**, **`\|μ\|<9×10⁻⁵`**, **95% CL**. Bounds spectral distortion, not universal loss rate. | Frequency/energy exchange producing those distortions (Part 1 rows 3–4; C3; O7); conditional CMB consistency test for tired light. Model→distortion map **gap**; not every spectrum-preserving frequency change is excluded by these numbers. | Accumulated sky-averaged CMB spectral history, over FIRAS observing band. |
   150	| TL — time dilation | DES, [MNRAS 533, 3365 (2024), primary abstract](https://arxiv.org/abs/2406.05050): 1504 SNe, **`0.1≲z≲1.2`**; `Δt_obs=Δt_em(1+z)^b`, **`b=1.003±0.005(stat)±0.010(sys)`**. Tests non-time-dilating cosmologies, not every leakage law. | Energy loss proposed **instead of expansion**, without corresponding time dilation; C3 and O7 still need an observable map. | Cosmological propagation measured through light-curve duration. |
   151	| TL — Tolman brightness | Lubin & Sandage, [AJ 122, 1084 (2001), primary abstract](https://arxiv.org/abs/astro-ph/0106566): 34 early-type galaxies, clusters **`z=0.76,0.90,0.92`**; brightness exponent **`n=2.59±0.17` (R), `3.37±0.13` (I)** for **`q₀=1/2`**. With luminosity-evolution assumptions, their tested nonexpanding tired-light model is excluded at **better than 10σ**. Published phenomenological finding, not a v3 verdict. | Energy loss replacing expansion, within surface-brightness/evolution assumptions. FD separately supplies CMB blackbody consistency. | Galaxy brightness versus redshift/zero-redshift fiducials. |
   152	| PB | Berceau et al., [PRA 85, 013837, primary paper](https://arxiv.org/pdf/1109.4792): vacuum magnetic birefringence **`Δn≤5.0×10⁻²⁰ T⁻²` per `4 ms` acquisition**. Phase/polarization response, not disappearance or generic defect birefringence. | Surviving transverse polarization response; corresponding magnetic/phase/defect map **gap** (PC:31). | Pulsed field with repeated cavity passes, stated acquisition. |
   153	| PS | Trapped shear can radiate if bulk shear exists (S9:70–73; SR:151–153). No record names a particle-stability or lifetime yardstick: the records assign owners only — holder/charge response to Q2/Q3/S22 (PC:21), and the throat/EM notes are “exploratory inputs for S22/Q2” (PC:37). Choosing a yardstick here is this inventory's own **gap**. No particle identity/lifetime chosen. | Trapped support decay into bulk/other receivers, distinct from free-photon propagation. | Trapped age/lifetime **gap** pending identity/support map. |
``````

</details>

<a id="r28"></a>
<details>
<summary>R28 — command and literal output</summary>

Command:

``````bash
nl -ba research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | sed -n '3,19p;128,148p;263,305p;351,369p'
``````

Exit code: 0.

Literal output:

``````text
     3	**Status: passes 1 (S9, S10) and 2 (S11, S11b-A/B, S11c PARTIAL) complete; pass-2 review pending.
     4	O2 pass populated 2026-10-07; its record/pass review is pending.**
     5	Sixteen entries, all OPEN. Pass 2 was populated 2026-10-05 from the kept results and their stated
     6	conditions; it does not close the substrate or upgrade S11c's unresolved work. Later sectors remain to
     7	be read as they close. The two 2026-08-07 prior-art entries retain their provenance; the second route is
     8	recorded under Population passes.
     9	O2 adds sources to six existing entries and one separate supplier/budget entry with no owner named;
    10	its four adopted premises retain their
    11	conditional-input labels in the O2 pass below. Entry counts are retrieved in
    12	`steps/_measurements/O2_record_measurements.md`, M11.
    13	
    14	## Why this exists
    15	
    16	v3 is requirements-first: each sector states what it needs, and the knit happens last. The light sector
    17	(S9, S10, S11, S11b) is built on **S1–S8, none of which have been run** — there are no step records for
    18	them. That is the design, not an oversight: the substrate gets built once every sector has said what it
    19	requires, so it can be checked against all of them at once rather than rebuilt per sector.
   128	### R-S1-01 — the brane's spatial dimension
   129	
   130	- **source** S10 (`steps/S10_two_transverse_photons.md`); S11
   131	  (`steps/S11_stray_longitudinal.md`, moves 2–3 / finite census) · **target** S6 · **status** OPEN
   132	- **requirement** — the brane's spatial dimension `D_brane`, as a derived quantity or an explicitly
   133	  re-affirmed postulate.
   134	- **on failure** — S10's headline reads *"light having exactly two polarisations is a statement that our
   135	  space is three-dimensional."* Read backwards it says the opposite: `D_brane = 3` went in and `D−1 = 2`
   136	  came out. Without a delivered `D_brane` the sentence is an assumption restated, ⛔ not a result.
   137	- **note** — ⚠ **The target is S6, not S1.** S1 owns `D = 4` and the two-phase split; `D_brane` is a
   138	  property of the **wall**, and S5–S6 are what construct it — S6 (the kink) is the first step at which a
   139	  codimension-1 surface exists to have a dimension.
   140	  ⛔⛔ **An earlier draft of this entry ruled out the obvious route on bad reasoning**, and a leg caught
   141	  it. It argued that because S10's *mode-count computation* contains no codimension, `D_bulk − 1` is not
   142	  the route. ⚠ **Non sequitur.** That the count does not *use* a codimension says nothing about whether
   143	  `D_brane` can be *derived* as one elsewhere — a domain wall in a `D = 4` bulk is exactly codimension 1,
   144	  and that is the natural derivation. ⇒ ⭐ **this entry does not prescribe a route**, per the schema; it
   145	  asks for the object.
   146	- **pass-2 consumer** — S11's selected homogeneous three-mode census and transverse/longitudinal
   147	  separation use `D_brane = 3`. Its proper-rotation invariant count explicitly warns that the sector
   148	  separation is not dimension-independent. This is the same dimension obligation, not a new knob.
   263	### R-S8-02 — the in-plane and out-of-plane sectors must decouple at quadratic order
   264	
   265	- **source** S10; S11 (`steps/S11_stray_longitudinal.md:278–280`, in-plane, frozen-wall-width spectrum)
   266	  · **target** S8 · **status** OPEN
   267	- **requirement** — the quadratic operator on the brane's **full** displacement, in-plane `u` **and**
   268	  out-of-plane `h` together, and whether it is block-diagonal in that split.
   269	- **on failure** — ⛔ **S10's headline number changes** — but ⚠ **not by the mechanism an earlier draft of
   270	  this entry named, and a leg corrected it.** Mixing is **not** the failure mode: under the rotation group
   271	  that acts on the brane, `h` and `u_L` are both **scalars** and `u_T` is a **vector**, so `h` can mix with
   272	  `u_L` freely and the transverse count stays `D − 1` regardless. ⭐ **The condition that gives 3 is `h`
   273	  being DEGENERATE WITH THE TRANSVERSE PAIR** — the brane's elasticity being isotropic in `D+1` rather
   274	  than in `D`. That is what *"belonging to the same elastic sector"* has to mean, and it is what S8 must
   275	  settle.
   276	- **note** — ⭐ **This is the most load-bearing open item in the sector**, and it was flagged inside S10 at
   277	  the moment the identification was made rather than found later: `h ≠ u_L`, stated in three places in the
   278	  corpus, **user-confirmed as the picture** — but ⛔ confirmed as a picture, not computed. `V3_STEP_PLAN`
   279	  puts *"transverse and longitudinal sectors; the reduced `h`/`u_L` operator"* in S8, so the object is
   280	  already scheduled.
   281	  ⚠⚠ **And v3 is not starting from nothing — a leg found the computation already exists in v2 and this
   282	  entry did not cite it.** `research/pde_ledger_v2/paper/stages/stage_030.tex:99-125` builds the coupled
   283	  `(u_L, h)` scalar block with stiffness `K = [[B_eff, C_hu], [C_hu, K_h]]` — explicitly **not**
   284	  block-diagonal, with `C_hu` a registered free-unreduced parameter — and `R79` records that the mixed
   285	  poles are cone-coincident only when `C_hu = 0`. ⇒ ⛔ *"nothing computes that"* is true of **v3** and
   286	  misleading about the **corpus**; S8 should start from `stage_030`, ⛔ not from scratch.
   287	  ⭐ Note this is the `(u_L, h)` block — the scalar sector — so it bears on the **longitudinal** slot and
   288	  the charge anchor, ⛔ and not directly on the transverse count, which is what the corrected on-failure
   289	  above turns on.
   290	- **pass-2 source** — S11's result here is the one its record names as the in-plane, frozen-wall-width
   291	  spectrum with `hBranon` excluded: *"⚠ Scope caveats, none touching the decoupling: (i) this is the
   292	  **in-plane, frozen-wall-width** spectrum (`WALL_WIDTH_FIELDS={}`, `hBranon` excluded,
   293	  `INTERFACE_EQUATIONS_SUPPLIED={}`) — it does not decide inhomogeneous mode conversion,"*
   294	  (`steps/S11_stray_longitudinal.md:278–280`). S11's decoupling statement, with its own limit, reads:
   295	  *"Thus the homogeneous \(D=3\) quadratic transverse and longitudinal eigenbranches have zero linear
   296	  cross-block under the selected action. This does not establish decoupling at nonlinear order, on a
   297	  nonuniform slab, at an interface or defect, or after additional allowed fields are introduced."*
   298	  (`steps/S11_stray_longitudinal.md:74–76`). This entry's object is the joint `u`/`h` operator that the
   299	  spectrum leaves out.
   300	- **excluded field** — in that caveat S11's record freezes the wall width and excludes `hBranon`
   301	  (`steps/S11_stray_longitudinal.md:279–280`); its Mathematica audit binds the token as
   302	  `"OUT_OF_PLANE_FIELD_EXCLUDED" -> hBranon`
   303	  (`mathematica/S11_stray_longitudinal_mathematica_audit.wl:1082`). If S8 does not deliver the joint
   304	  operator, S11's spectrum is not extended beyond that stated scope. The record does not say what its
   305	  roots or census would be with `hBranon` included, and this register does not infer it.
   351	  with an independent microrotation degree of freedom (Cosserat/micropolar), and media with stored internal
   352	  angular momentum.
   353	  ⛔ **One family is the wrong one and should not be reached for:** modern *odd elasticity* buys an
   354	  antisymmetric modulus tensor by making the solid **active and non-conservative**. MacCullagh's medium is
   355	  **conservative** — it has a genuine energy functional — so a non-conservative realisation would be
   356	  answering a different question.
   357	  ⚠ `R-S8-01` asks for the **form** and is silent on admissibility; a substructure could deliver the
   358	  curl-only form and still owe this.
   359	- **O2 scope** — OPEN `𝒜_rot^live` carries this existing S8 obligation, including whether internal
   360	  angular momentum/couple stress is present. O2 selects no stress symmetry or carrier and does not
   361	  establish a curl-only live stress. The original conservative-MacCullagh objection keeps its
   362	  original domain; it is not newly proved for an unspecified viscoelastic stress. The O2 addition
   363	  carries the explicitly required carrier/admissibility question within its named material and
   364	  corresponding power accounting; it adds no flowing constitutive response form. The accounting
   365	  remains conditional with that question OPEN
   366	  (`directives/O2_SHARED_PHYSICS.md`, §§3.2, 4, 6).
   367	
   368	### R-S8-05 — the frame the brane's rotational stiffness is measured against
   369	
``````

</details>

<a id="r30"></a>
<details>
<summary>R30 — command and literal output</summary>

Command:

``````bash
nl -ba research/pde_ledger_v3/steps/S11b_interface_coupling_law.md | sed -n '1,40p;70,112p'
``````

Exit code: 0.

Literal output:

``````text
     1	# S11b — the linear brane–bulk interface coupling law (step record)
     2	
     3	Slug `S11b_interface_coupling_law`. This is the unified step record for S11b, subsuming the historical A
     4	(bulk face response) + B (homogeneous assembly) execution stages into **one export-chain step**; **C — the
     5	non-uniform transverse coupling — is a separate later step** (deferred, scope ratified in the decision list
     6	`directives/S11b_unified_decisions.md` G14, `ddd0ae4c`). The physics authority is the unified spec
     7	`directives/S11b_SHARED_PHYSICS.md` (`1a2395a3`).
     8	
     9	⭐ This record is the **interpretation** layer; every result below names the computed object and the commit
    10	that produced it. The engines PRINT objects; they state no conclusions (the conclusions are here).
    11	
    12	## What the step computes
    13	
    14	The linear response of the brane–bulk interface: given a slab of finite thickness `W` in the `w` direction
    15	with two faces meeting the bulk, derive how a normal face motion couples to the bulk, what drives material
    16	across the interface, and the closed spectrum of the assembled system on a **uniform** background.
    17	
    18	## The two engines and their agreement
    19	
    20	- **SymPy engine** `scripts/S11b_interface_coupling_law_sympy_audit.py` — the packaging engine: imports the
    21	  S11 LEDGER (1663 rows), binds `c_s0`/`μ_R`/`ρ_br⁰=rho_br` to the imported objects, writes
    22	  `S11b_exports.py` (1958 rows). Committed `864d6f41`; X-1 repair `53fcd98d`.
    23	- **Wolfram engine** `mathematica/S11b_interface_coupling_law_mathematica_audit.wl` — **blind**: imports
    24	  nothing, re-derives every object from §1–§6 of the spec. Committed `ec89f9df`; repair `bd598ae7`. Its
    25	  blindness is the only cross-engine control (it cannot transcribe the `.py` because it never reads it).
    26	- **T7 comparator** `scripts/S11b_cross_engine_comparator.py` (reused from S10) — joins the two transcripts
    27	  by emitted-object name, residuals paired payloads, rejects a native boolean as a residual operand, is
    28	  three-valued, and was frozen before it saw either output. Re-run against both repaired engines: `fba6a34c`
    29	  (`scripts/out/S11b_cross_engine_comparison.out`).
    30	
    31	⭐ **No compared object is a physics contradiction.** After the WL repair (F-WL-1) and the SymPy X-1
    32	repair, the comparator's headline is AGREE 21 / DISAGREE 108 / UNCOMPARED 25 / UNPAIRED 71; the 108
    33	disagreements are **format, coefficient-basis, naming, convention, and a few coverage gaps — not conflicting
    34	values** (adjudicated below). ⚠
    35	**Finding — the engines are not emission-parallel:** a SymPy engine emits Python `Tuple`s and a Wolfram
    36	engine emits `Association`s, so the comparator marks 102 objects `STRUCTURE_DISAGREE(tuple vs Association)`
    37	even where the physics agrees (the parsed values coincide, or differ only by the invertible coefficient-basis
    38	remap of the representative split — see the adjudication below). This is a property of the two CAS output
    39	shapes and the basis choice, not a physics divergence; it is why `FINAL_OPERATIONAL_STATUS` reads `FAIL`
    40	while the physics agrees.
    70	**The velocity channel and its thermodynamic fate.** The velocity-driven conversion is an energy **source**
    71	unless it is instantaneous (`τ → 0`) — inside the passive region the second law forbids it, so the model
    72	refuses an unphysical process on its own terms. Outside the region it survives only against a named reservoir.
    73	
    74	**Transverse mode — decoupled on a uniform background.** The transverse-to-thickness coupling
    75	`∂²U/∂u_T ∂e_W` is **identically zero** and the transverse mode is **non-dissipative** (`TRANSVERSE_DISSIPATION`
    76	≡ 0, unconditionally) — a decoupled oscillator `ρ_br⁰ ω² = μ_⊥ k²` (`S11B_TRANSVERSE_DISPERSION`,
    77	`S11B_TRANSVERSE_COUPLING`). ⚠ It is a **stable real-frequency** mode (`Im ω = 0`) **only where the transverse
    78	stiffness `μ_⊥ ≥ 0`**: §5's moduli are free and ⛔ no positivity is assumed (§0), so `μ_⊥ = μ_R + μ_S/2 < 0`
    79	is admissible and gives a growing transverse root `ω = ±i k √(|μ_⊥|/ρ_br⁰)` — decay and growth are both
    80	admissible outcomes here, ⛔ not excluded. ⛔ This is the **uniform** limit only and does NOT settle
    81	unconditional confinement — that is the non-uniform question deferred to C. ⚠ The two
    82	engines express the transverse stiffness `μ_⊥` in **different coefficient bases** (a consequence of the
    83	energy-basis representative split, below): the blind WL engine emits `μ_⊥ = μ_R`, the SymPy engine
    84	`μ_⊥ = μ_R + μ_S/2`, with identical roots under the invertible coefficient identification WL-`μ_R` ≡
    85	SymPy-`(μ_R + μ_S/2)`. ⛔ They are NOT equal under the comparator's naive `muR↔mu_R` name-transliteration,
    86	which is why `TRANSVERSE_DISPERSION` shows a cross-engine difference (an extra `k²μ_S/2`, masked under its
    87	tuple-vs-Association STRUCTURE flag) — the same physical stiffness, not a physics disagreement.
    88	
    89	**Breathing-mode stability.** On the `k = 0` slice with impermeable faces (`Λ_A⁰ = Λ_V⁰ = 0`) **and no
    90	reciprocal traction (`Λ_X⁰ = 0`)** — all three cuts, since "impermeable" alone does not set `Λ_X⁰` — the
    91	bulk load is kept and `K₀ = B_ρ⁽³⁾ − 2C W₀ + k_W W₀² > 0`; growth iff `C > (B_ρ⁽³⁾ + k_W W₀²)/(2W₀)`. ⚠ The growing root is **not** an energy-conservation violation — the
    92	stored energy has no minimum in that direction and the accounting closes exactly. ⚠ This is a slice result,
    93	⛔ not the general breathing boundary.
    94	
    95	## The energy basis and the X-1 correction
    96	
    97	⭐ The spec (§5) carried **five** stored-energy terms but instructed each engine to CONSTRUCT the symmetry-
    98	allowed basis and check its closedness (`S11B_ENERGY_BASIS_COUNT`). The engines found more allowed
    99	invariants than the spec named — including `(∇·u)²`, the coupling that fell out of the specification four
   100	times while it lived in prose.
   101	
   102	⛔ **X-1 — a one-engine over-count, caught by the sibling.** The SymPy engine first emitted **11** basis
   103	invariants; the blind WL engine emitted **10**. §5's symmetry group states *"equivalence modulo total
   104	divergences — two densities differing by a total in-plane divergence are the same term; do not count both."*
   105	The SymPy independence test judged only pointwise polynomial independence over the field components and
   106	omitted that quotient; two of its `∇u`-sector invariants differ by a total in-plane divergence
   107	(`st² ≡ ½curl² + ⅔(∇·u)²`, since `(∇·u)² − tr((∇u)²)` is a total divergence). The spec-correct count is
   108	**10**; the SymPy engine over-counted. Verified independently three ways (the two build legs' own
   109	Euler–Lagrange enumerations, and the orchestrator's own `x1_independent_basis_count.py`).
   110	
   111	The X-1 repair (`53fcd98d`) made the SymPy independence judgment honor §5's quotient (Euler–Lagrange-signature
   112	rank) and eliminated the redundant invariant by **REWRITING** it into the retained curl²/strain invariants —
``````

</details>

<a id="r31"></a>
<details>
<summary>R31 — command and literal output</summary>

Command:

``````bash
rg -n -i '\b(jones|poincaré|poincare|helicity|handedness|entanglement|entangled|malus|faraday|kerr)\b|bell test|spin hall|circularly polarized|elliptical polarization' research/pde_ledger_v3/V3_STEP_PLAN.md research/pde_ledger_v3/DEFECT_REGISTER.md research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md research/pde_ledger_v3/steps/*.md
``````

Exit code: 1.

Literal output:

``````text
``````

</details>

<a id="r34"></a>
<details>
<summary>R34 — command and literal output</summary>

Command:

``````bash
nl -ba research/research_proposal/paper/proposal_chapter_04_light_research_program.tex | sed -n '35,68p;325,340p'
``````

Exit code: 0.

Literal output:

``````text
    35	transverse elastic branch. Light is identified with the lowest physical
    36	transverse branch guided by this finite material structure. It is not a second
    37	substance and no fundamental electromagnetic field is inserted solely to carry
    38	it.
    39	
    40	Let the brane possess a local normal axis represented by a unit vector
    41	\(N^A\). The sign of \(N^A\) distinguishes the two throat orientations used by
    42	the charge sector, but ordinary free light should depend at leading order only
    43	on the unsigned axis, through reflection-even structures such as \(N^A N^B\).
    44	For a wave vector tangent to the brane, a light displacement is required to be
    45	both tangent to the slab and transverse to the direction of propagation. In the
    46	four-dimensional parent space, excluding the normal direction and the
    47	propagation direction leaves two independent tangential polarization
    48	directions. The complete guided spectrum must realize that count dynamically;
    49	geometrical counting alone is not sufficient.
    50	
    51	Under the currently supplied homogeneous, isotropic quadratic brane action, the
    52	two transverse modes have
    53	
    54	\[
    55	\omega_T^2
    56	=
    57	\frac{\mu_R}{\rho_{\mathrm{br}}}k^2,
    58	\qquad
    59	c_\gamma^2
    60	=
    61	\frac{\mu_R}{\rho_{\mathrm{br}}},
    62	\]
    63	
    64	where \(\mu_R\) is the effective ordered-state shear modulus and
    65	\(\rho_{\mathrm{br}}\) is the effective inertia density of the brane mode. The
    66	same action contains one in-plane longitudinal branch,
    67	
    68	\[
   325	three-dimensional vector envelope and slaved helper. A one-dimensional success
   326	is therefore a precursor; it must survive transverse stability, full
   327	three-dimensional continuation, and reconstruction in the parent fields.
   328	
   329	The reduced vector system must preserve the symmetries of the parent slab. On a
   330	homogeneous isotropic background, no fixed in-brane direction may be introduced
   331	to collimate one chosen propagation axis. Any effective guide should align with
   332	the packet that creates it. Polarization-dependent nonlinear invariants must be
   333	derived explicitly so that linearly, circularly, and elliptically polarized
   334	families can be compared. A helper sourced only by total intensity may preserve
   335	leading polarization degeneracy, whereas other allowed invariants may split
   336	energies, speeds, widths, or stability. Such splitting is a prediction to be
   337	bounded, not a coefficient to be independently tuned for each polarization.
   338	
   339	Existence of one isolated numerical lump would also be insufficient. The desired
   340	result is a connected traveling family over a declared carrier-frequency range,
``````

</details>

<a id="r35"></a>
<details>
<summary>R35 — command and literal output</summary>

Command:

``````bash
nl -ba docs/conceptual_foundation.md | sed -n '366,394p;470,478p'; nl -ba docs/conceptual_history.md | sed -n '248,267p'; nl -ba docs/em_sector_reconsideration.md | sed -n '1,47p'; nl -ba docs/model_map.md | sed -n '105,127p'
``````

Exit code: 0.

Literal output:

``````text
   366	|---|---|---|
   367	| **Gravity** | The medium's **inflow / drain** toward defects. Test bodies are carried inward by the flow. Largely unobservable as a flow. | flow velocity `v_r` |
   368	| **Magnetism** | Two swirls, kept distinct (⭐ v8 block above): **EM magnetism** = the velocity-dependent (`O(V)`) part of the 4D throat-body interaction (the *moving* version of electric charge — charge-coupled, `±w` sign, strong) vs **gravitomagnetism** = the 3D flow-swirl of the mass-inflow (mass-coupled, weak frame-dragging; part of the *gravity* sector). The bulk is shear-free (it cannot "push"), which is *why* frame-dragging is tiny and light stays confined to the brane. | throat-body swirl (`±w`) / flow vorticity |
   369	| **Electric charge** | The **puncture direction** — which way a throat punctures the brane into the bulk (`+w` vs `−w`). A **binary orientation**. No rotation, no swirl, no winding. | puncture orientation `±w` |
   370	| **Light** | The brane's in-plane **rotational-elastic (MacCullagh) shear wave** — our 3D space resisting being *locally twisted*, the wave riding on the curl of an in-plane displacement (not the drumhead's up/down bowing, which is a separate scalar mode). Two transverse polarizations, no longitudinal mode; rigidity lives **on the brane, not the bulk**. | brane rotational-elastic shear |
   371	
   372	**Why two polarizations and no third (the MacCullagh point):** if the brane's stiffness is **rotational/curl-type** — it costs
   373	energy to *locally rotate* the medium, energy `∝ (∇×u)²` — then the wave equation has **exactly two transverse polarizations and
   374	no longitudinal mode**: precisely electromagnetism. (Ordinary "Cauchy" elasticity, which resists *stretching/sliding strain*
   375	`∝ (∂u)²`, would instead give a spurious *third*, longitudinal "photon" that real light does not have — so light needs the
   376	*rotational* kind of stiffness specifically, not generic shear.) This is MacCullagh's 1839 rotational ether, the historical
   377	ancestor of the EM field. It is the open "structured wall" question of §2 — see the deep-dive just below.
   378	
   379	**The summary image (keep this):** *Our 3D space is a taut elastic sheet in the 4D bulk; gravity is how that sheet drapes and how
   380	the bulk drains through it; light is the sheet shivering sideways within itself; and the bulk underneath stays a frictionless
   381	fluid so magnetism never feels the stiffness.*
   382	
   383	**What light actually is — and the two ways a wiggle can carry energy (read this; it's the part that's easy to get lost on).**
   384	We were never going to "derive photons" out of the bare scalar superfluid — a scalar has only compression waves (the `c_s` ripple,
   385	the gravity sector). We *postulate a mode* (the arrows) and ask whether it behaves like light. Picture the arrows lying flat in the
   386	brane and wiggle them. There are **two physically different ways the wiggle can store energy, and only one is light:**
   387	- **Frank / orientation wave (NOT light).** The arrows re-point while the medium stays put — like a stadium wave: the people
   388	  (material) don't move, only their "pointing" propagates. The energy cost is for *neighbours to disagree* — an orientation
   389	  gradient, `∝ (∇n)²` — and the restoring force is a **torque**. It can even have two transverse components and a fixed speed, so
   390	  it can *look* like light on a dispersion plot — **but it is not**, because nothing material is displaced and it exerts no
   391	  mechanical force across a surface.
   392	- **MacCullagh rotational-elastic wave (THIS is light).** Here the medium resists being *locally rotated at all*; energy
   393	  `∝ (∇×u)²`, the curl of a genuine displacement. This gives exactly two transverse polarizations, no longitudinal mode, Maxwell —
   394	  *real* electromagnetism. It is special because ordinary matter resists squeezing and sliding but **not** rotating; to resist
   470	  (one direction lower-energy or more stable), with the opposite orientation being the antiparticle. The observed cosmic
   471	  *abundance* asymmetry (why matter dominates) would then be a *separate* question — plausibly about a small asymmetry in the
   472	  bulk or in the two brane-states — layered on top of the structure↔direction preference. Both halves are open and interesting.
   473	
   474	**Sharpening (2026-06-24 — light's dimensionality, how the throat traps it, and the honest death). Working synthesis:**
   475	- **Light is INTRINSICALLY a (3+1)D brane field** — so it has **2 transverse polarizations automatically**, and it **never leaves the
   476	  brane**. *Why it can't escape:* the bulk is **shear-free** and light **is** shear → there is no medium for light off the brane. Its
   477	  "4D-ness" is the brane's **extrinsic curvature/embedding** (the surface bends into `w` at a throat), **not** a bulk excursion or a 3rd
   478	  polarization. (Like a 2D wave on a sheet of paper rolled into a 3D tube: intrinsically 2D, extrinsically 3D.) This resolves the old
   248	> slab at `w≈0`, de-structured outside) is a **genuinely different object** from both corpses:
   249	> - It **escapes T1** (the little-arrows wall that *spread and unwound* on a **connected** `S³` vacuum): a scalar `χ_B` with a real
   250	>   **double-well** free energy `f_B(n,χ_B)` has **disconnected** minima → a **φ⁴-kink-like, topologically stable wall** — exactly
   251	>   the structure T1's own negative control (the stable φ⁴ kink) showed it lacked.
   252	> - It **escapes the density no-go**: `χ_B` is a **single domain wall in an abstract order field, not a periodic density stack**, so
   253	>   it has **no layer normal to pin the arrows against** — the in-plane shear lives in a separate field that `χ_B` merely *gates on*.
   254	>
   255	> **The exciting sub-hypothesis — `χ_B`'s phase may supply the missing C5 `φ`.** The shear-surface death was precise: MacCullagh's
   256	> curl-only energy `½μ_R(∇×u)²` is gauge-invariant under `u→u+∇χ` but its kinetic term is not, forcing a **constrained physical
   257	> longitudinal zero mode** (`∂_t²(∇·u)=0`); Maxwell escapes via a scalar potential (`φ→φ−∂_tχ`) but MacCullagh has none, and the
   258	> only on-brane scalar (`u_w`) **must stay gapped** (a massless out-of-plane mode = an excluded fifth force), so it *cannot* be `φ`.
   259	> But if the brane-order parameter is **complex** (amplitude `χ_B` + **phase `θ`**), that phase is a **new on-brane scalar with
   260	> genuine mechanical provenance that is NOT `u_w`** — a natural candidate for the exact `φ`-analog the frozen light package lacked.
   261	> If `θ` couples to the displacement like Maxwell's scalar potential, it removes the longitudinal zero mode and **could revive
   262	> shear-surface light** as a *fresh G0*. Whether it actually does is the test — **not a claim**.
   263	>
   264	> **Honest costs (do not fool ourselves).** The double-well `f_B` is **still postulated** — the GNLS potential `U(ρ)∝ρ⁵` is
   265	> single-well, so "two coexisting phases" is the *same* degenerate-vacua obstacle the brane has always had (§ "The obstacle" below),
   266	> now conceded under the analog license (§0.6), **not derived**. It is **real drift**: `f_B` + the interface stiffness `κ_B` + the
   267	> `χ_B`-shear gating + (if complex) the phase sector — the gauntlet must count it and decide whether `χ_B` is
     1	# EM Sector Reconsideration (v2, panel-hardened) — magnetism from the trapped-light geon in the MacCullagh shear
     2	
     3	**Status: WORKING DRAFT v2 (2026-07-11).** For the EM-sector re-derivation this supersedes v1 and the `conceptual_foundation.md` §3 v8 magnetism block. v1 was reviewed by three independent AI panels (Codex, Grok, GLM); v2 folds in their objections, **corrects errors v1 made**, and installs a new leading mechanism (the circulating trapped-light **geon**) that the panel's own critiques pointed toward. Still **not** a settled result — a *sharpened hypothesis with a cheap decisive test*. **Falsification welcome; the EM sector may still break.**
     4	
     5	---
     6	
     7	## 0. What changed from v1 (honest ledger)
     8	
     9	- **Corrected error.** v1 called `pathA_39` stage-3 an "over-read." All three reviewers showed this is **wrong**: stage-3 *computed* that the parity split **leaks under motion** (27 combined-parity-mixing operators, 7 nonzero static witnesses). The magnitude/helicity split (v1 §2.3) is therefore **approximate, with an unknown mixing angle — not a clean legitimate split.**
    10	- **Demoted (panel-refuted as stated).** "Order-parameter elasticity sets the sign, free either way" (v1 §2.5) and "finite throat-body geometry projects a 4D interaction to `1/r²`" + "long-range Coulomb → `1/r`" (v1 §2.6). See §4.
    11	- **Installed as leading mechanism.** The **circulating trapped-light geon** carried in the **MacCullagh rotational-elastic (shear/light) sector** — which addresses the panel's *spin gap*, its *light↔magnetism unification* demand, and its *deepest objection* (the gauge / dual-sign problem). See §2–§3.
    12	- **Replaced the plan.** The decisive first test is now the **gauge-structure / Maxwell-dual-sign linear test** in the shear sector — not the falloff projection. See §5.
    13	
    14	---
    15	
    16	## 1. The problem (confirmed unanimously by the panel)
    17	
    18	`pathA_39` landed "like currents attract, correct EM sign." The panel (Codex, Grok, GLM) unanimously confirmed the headline is **not earned**:
    19	
    20	- **The sign is baked in.** "Attract" is fixed *before the integral* by three modeling choices: an EM-style current source `j = qV` (linear in velocity), positive-definite propagators, and the exchange rule `U = −jGj`. The ghost control (`μ_R → −μ_R` flips attract↔repel) confirms the sign tracks a *stability condition*, not a dynamical prediction. Verdict `CONDITIONAL_CURRENT_EXCHANGE_ATTRACTION; NATIVE_VORTEX_SIGN_UNDERIVED`.
    21	- **The definition welds two incompatible pictures** — a *swirl* (`v_φ ~ 1/r`, a vortex) given a *current's* force sign, called "FALLING OUT (not asserted)" when it was not.
    22	- **The `1/r²` was target-driven** — brane-localization "borrowed from gravity" and "must be shown, not assumed."
    23	- **The deeper structural finding (panel):** real EM's **dual sign** — like charges *repel*, like currents *attract* — is a **gauge-theoretic** fact. It is not a coincidence of stiffnesses; it follows from a constraint (Gauss) structure. Any single mediator without gauge structure generically gives the *same* sign for both channels. This is the organizing problem for v2 (§3, §5).
    24	
    25	---
    26	
    27	## 2. The reorganized picture (v2)
    28	
    29	### 2.1 Different-dimensional forces — sharpened, and softened where the panel corrected it
    30	- **Gravity = the drain = 3D-brane mass-flow.**
    31	- **Charge = the `±w` puncture** — genuinely bulk-*facing* (the throat body pokes into the bulk), **but its force is mediated by a brane mode** (the `h`-branon, a normalizable brane displacement into `w`; `pathA_38`). That brane-mode mediation is *why* Coulomb comes out `1/r²`.
    32	- **Magnetism = trapped circulating shear** — a **brane** phenomenon (§3).
    33	
    34	So the honest form of the different-dimensional hypothesis is: **EM's identity/orientation is bulk-facing (`±w`), but EM's forces are brane-mode-mediated** — which is exactly why they are `1/r²` *and* why they escape the "shear-free bulk can't carry a long-range interaction" trilemma. The clean "all of EM is a 4D-bulk force" of v1 is retired.
    35	
    36	### 2.2 EM is the ROTATIONAL-elastic (MacCullagh) sector — NOT a generic scalar/Frank texture
    37	Panel correction: a *generic* order-parameter texture (Frank/XY defect) gives like-defect **repulsion** and cannot, by itself, produce Maxwell's dual sign. So EM is **not** "any order parameter." It is specifically the **MacCullagh rotational-elastic shear** — the light sector — whose curl-only energy `∝ (∇×u)²` is the **historical ancestor of Maxwell** and carries a **gauge-like structure** a scalar `χ_B` lacks (§3). Gravity remains the mass-flow sector.
    38	
    39	### 2.3 The parity split is approximate, not clean (corrected)
    40	Stage-3 (`FAIL_UNPROTECTED_OPERATOR_PARITY_MIXING`) *computed* that a moving `P_w`-odd charge is symmetry-allowed to mix the odd (EM) and even (gravitomagnetic) sectors — the split **leaks** under motion. The magnitude(mass)/helicity(`±w`) decomposition therefore survives only as an **approximation with an unknown mixing angle**; angular momentum can flow between channels as ordered medium de-structures. (Also: in 4 spatial dimensions vorticity is a 2-form, not a 3D axial vector — there is no unique 4D "right-hand rule"/scalar helicity, so the "handedness" language must be made precise, not assumed.)
    41	
    42	---
    43	
    44	## 3. The leading mechanism (v2): the circulating trapped-light geon
    45	
    46	**The keystone: the standing wave holding the throat open is not static.** Mass = a trapped standing-wave geon (existing model). A standing wave is counter-propagating/circulating energy; a geon can carry genuine **angular momentum** — a **spin** that lives in the *circulating trapped light*, with **no rigid rotation of the medium at all.**
    47	
   105	### 3.3 Light — Part III ✅ DONE = stage003 (re-scoped 2026-07-22, surviving-solution rule)
   106	
   107	*Light = in-plane MacCullagh shear of the brane; two transverse photons. The shared transverse sector — it establishes the `u_T`/`c_γ` foundation magnetism reuses. Under the surviving-solution rule Part III is stage003 alone; the `pathA_35` couple-stress no-go (retired-`P` post-mortem) → failures-paper backlog.*
   108	
   109	⭐⭐ **v3 RE-DERIVED THIS SECTOR FORWARD (2026-08-02) — S9/S10/S11, and it says more than stage003 did.**
   110	⇒ `research/pde_ledger_v3/steps/` · ⛔ read `S11_stray_longitudinal.md` before citing the stray mode.
   111	- ⭐ **Two polarisations is a statement that our space is 3-dimensional.** The count is `D_brane − 1`,
   112	  **computed** as a nullity in both engines; ⛔⛔ the **bulk never enters**, so codimension is not the
   113	  wrong number, it is an **absent quantity**.
   114	- ⭐ **Sector separation is a `D=3` fact, ⛔ not a structural one** — at `D=2` a reflection-odd invariant
   115	  exists whose EL operator is non-zero, so compression and light *can* mix there.
   116	- ⭐⭐ **The stray mode's ENTIRE physical status — radiative, observable, Lorentz-breaking, or none —
   117	  reduces to ONE unbuilt object: the brane–bulk interface coupling law** (`V3_STEP_PLAN.md#s11b`). It is
   118	  **LINEAR**, so it is half-one work that was deferred by choice.
   119	- ⛔⛔ **`c_L` is NOT a gravitational wave** (`DEFECT_REGISTER.md#c13`) and ⛔ `u_L` is **not** `±w`.
   120	- ⭐ **A second cone is a departure only if matter COUPLES to it.** Matter here is built from the
   121	  transverse modes, whose wave equation is Lorentz-invariant with invariant speed `c_γ` automatically.
   122	- ⭐⭐ **Photon stability over cosmological distance is a CONSEQUENCE of bulk shear-freeness** — light can
   123	  only lose energy into a mode that exists to receive it, and the bulk has no transverse mode. ⇒ S9's
   124	  second requirement is an observational consequence, ⛔ not bookkeeping. ⚠ Light **gravitates** (a
   125	  co-moving disturbance) but does ⛔ **not dissipate**; the two were conflated once.
   126	- ⚠ **Lensing sign UNCHECKED** — the naive drain reading points the wrong way
   127	  (`DEFECT_REGISTER.md#c14`); ⛔ do not manufacture a reconciliation.
``````

</details>

<a id="r36"></a>
<details>
<summary>R36 — command and literal output</summary>

Command:

``````bash
nl -ba research/4d_em_fields/paper/4d_em_fields.tex | sed -n '1411,1434p'; nl -ba research/1pn_hybrid/paper/1pn_hybrid.tex | sed -n '1798,1815p'; nl -ba docs/toy_model_ontology_summary.md | sed -n '943,966p'; nl -ba docs/cross_sector_research_burdens_and_compatibility_gates.md | sed -n '207,230p;264,275p;969,979p'
``````

Exit code: 0.

Literal output:

``````text
  1411	These Yukawa corrections modify the field profile of a fixed topological charge
  1412	branch. They do not imply any dependence of electric charge on circulation,
  1413	throat radius, or geometry breathing.
  1414	
  1415	For time-dependent sources, each KK mode obeys a covariant massive wave equation
  1416	\cref{eq:mode_wave_eq} with retarded solution \cref{eq:mode_retarded_solution}.
  1417	Below the first brane-coupled threshold \(\omega<m_2\), massive modes contribute
  1418	only near-field/evanescent corrections and the response is extremely close to
  1419	Maxwell at distances \(r\gg\lam\). Above threshold, \(\omega>m_2\), additional
  1420	propagating channels open with dispersion
  1421	\begin{equation}
  1422	  \omega^2 = k^2 + m_n^2,
  1423	\end{equation}
  1424	and group velocity \(v_g = k/\omega < 1\). This predicts frequency-dependent
  1425	response (and, in principle, additional polarization content) controlled by the
  1426	discrete tower \(m_n^2=2n/\lam^2\) and couplings \(c_n\).
  1427	
  1428	\subsection{Falsifiable tests and observational handles}
  1429	The localization sector can be tested in several complementary ways. The
  1430	discussion below is organized by which assumption in the controlled reduction is
  1431	being probed.
  1432	
  1433	\paragraph{(i) Precision tests of Coulomb's law (static).}
  1434	The cleanest signature is a Yukawa departure with \emph{fixed coefficient pattern}
  1798	Finally, the model invites direct comparison with astrophysical data.  Key
  1799	targets include:
  1800	\begin{itemize}
  1801	  \item black--hole imaging (shadow size and shape, photon--ring
  1802	        structure) in systems where independent mass estimates are
  1803	        available,
  1804	  \item precision timing of compact binaries in regimes where 2PN and
  1805	        higher corrections are observable, and
  1806	  \item potential signatures in high--energy emission or polarization
  1807	        patterns arising from the unified superfluid description of
  1808	        gravity and electromagnetism.
  1809	\end{itemize}
  1810	In each case, the stiff \(n=5\) fluid and \(M \propto \rho\) relation
  1811	lead to concrete, quantitative predictions once the throat parameters
  1812	are fixed.
  1813	
  1814	\vspace{0.5em}
  1815	\noindent In summary, this paper has shown that the 1PN success of the
   943	\[
   944	\nabla\cdot\mathbf u_T=0.
   945	\]
   946	
   947	The supplied quadratic brane action supports two transverse modes in three brane dimensions. Their characteristic speed is
   948	
   949	\[
   950	c_\gamma^2=\frac{\mu_R}{\rho_{\rm br}},
   951	\]
   952	
   953	where \(\mu_R\) is the effective brane shear modulus and \(\rho_{\rm br}\) is its effective inertia density.
   954	
   955	The two-polarization count is conditional on the assumed action, isotropic inertia, and the relevant structural premises. The model has not yet derived the shear modulus from the underlying medium.
   956	
   957	Freely propagating light and the throat-support mode must be distinguished. Free light is a background-brane \(\mathbf u_T\) excitation. The support mode is a spectrally normalizable bound state or acceptably long-lived resonance of the complete variable-coefficient transverse operator; it is not a second photon substance and is not assumed to propagate as a shear wave through fully de-structured bulk material. Any \(\Omega_{\rm support}\) is a diagnostic energy/stress-localization region extracted from that solved mode, not a hard PDE domain imposed by \(\mu_R>0\).
   958	
   959	Whether transverse light is sufficiently confined to the brane is an interface question. Uniform transverse coupling to the bulk vanishes under the supplied assumptions, but the physically relevant nonuniform, gradient-driven coupling remains to be completed. The model must show that bulk leakage, longitudinal conversion, or defect-induced mixing does not spoil ordinary transverse propagation. It must also derive how an accelerating oriented throat excites the two far-zone transverse polarizations through the same coupling structure as static electricity and magnetism.
   960	
   961	**Status:** the conditional two-transverse-mode count follows from the supplied quadratic action; the origin of the shear modulus, nonuniform interface isolation, and propagation through defect backgrounds remain open.
   962	
   963	## 10. The stray longitudinal mode
   964	
   965	The same brane can also carry an in-plane longitudinal displacement \(u_L\). In an ordinary symmetric elastic response, the medium consistently supports both:
   966	
   207	\]
   208	
   209	and the static potential behaves as \(1/R\).
   210	
   211	The ontology also states that the present system does not possess an established first-class \(U(1)\) Gauss constraint. It therefore does not yet have the conventional electromagnetic mechanism by which a static Coulomb field is tied to a constrained gauge sector rather than an extra freely propagating scalar-like polarization. The moving-throat program correctly requires emitted powers into the odd sector, transverse light, even thickness, longitudinal, bulk, and conversion channels to be calculated separately.
   212	
   213	## Reviewer inference
   214	
   215	Once dynamics are supplied, a long-range static odd carrier must fall into some version of the following possibilities.
   216	
   217	### A. Ordinary local inertia
   218	
   219	With a conventional positive kinetic term, the same \(k^2\) static stiffness generally produces a gapless propagating branch,
   220	
   221	\[
   222	\omega_E(k)\simeq c_Ek.
   223	\]
   224	
   225	Then the model possesses an additional dynamical mode beyond the two transverse photon polarizations. Because throats must couple to it strongly enough to produce Coulomb forces, moving or accelerating throats may excite it. The resulting consequences can include:
   226	
   227	- extra radiation;
   228	- an additional polarization-like channel;
   229	- modified radiation reaction;
   230	- mode conversion near defects;
   264	
   265	The selected electric branch fails if every healthy dynamic completion that supplies the required static \(1/R\) response also necessarily produces one or more of the following at an unacceptable level:
   266	
   267	- unsuppressed extra radiation;
   268	- low-frequency drag;
   269	- an active or nonpassive reservoir response;
   270	- a ghost or negative-energy mode;
   271	- instability or strong coupling;
   272	- a third freely radiating electromagnetic polarization incompatible with the intended phenomenology;
   273	- independent electric, magnetic, and radiation normalizations;
   274	- violation of the conserved oriented-current target.
   275	
   969	- distinction between same-species annihilation and different-species neutral binding.
   970	
   971	---
   972	
   973	# 10. Predeclared acceptance budget
   974	
   975	The documents often and correctly use the phrase “acceptable.” Before broad parameter searches, each acceptable departure should be assigned a frozen numerical or scaling criterion. The exact numbers are a later program decision, but the quantities should include:
   976	
   977	- relative polarization speed or energy splitting;
   978	- transverse guided-mode leakage ratio \(\Gamma_0/\omega_0\);
   979	- brane-participation fraction and pole residue;
``````

</details>



<a id="p3"></a>
<details>
<summary>P3 — repair source, command and literal output</summary>

Command:

``````bash
nl -ba research/pde_ledger_v3/V3_STEP_PLAN.md | sed -n '453,461p;871,905p;1174,1184p'; nl -ba research/pde_ledger_v3/steps/S11_stray_longitudinal.md | sed -n '15,26p;154,176p'; nl -ba research/pde_ledger_v3/steps/S11c_d_profile_conditioned_scattering.md | sed -n '52,64p'
``````

Exit code: 0.

Literal output:

``````text
   453	⛔⛔ **DIFFERENT OBJECTS — ⛔ do not weld them, at any strength.** `u_L` is an **in-plane longitudinal
   454	brane displacement**; **`±w`** is the **throat's normal/orientation direction**, and the charge the
   455	rethink produced is carried by the **`h`-branon**. The committed sources keep them apart:
   456	
   457	- `research/pde_ledger_v2/notes/stages/ledger_stage031_puncture_deflection_field_identity_source.md#h-a-distinct-from-u-l`
   458	  — *"Define the held mouth datum `h_A = ξ_w|_A/ℓ = P₀H|_A`, `[h_A]=1` **(distinct from `u_L`)**"*.
   459	- `docs/em_analog_next_phase_handoff.md#u-l-clamp-not-the-charge-scalar` — of the `u_L`-clamp candidate:
   460	  *"It is NOT the committed charge scalar (that's the `h`-branon, `h≠u_L`; `h` remains the committed
   461	  mediator for the conditional `1/R²` falloff); `u_L` is a separate charge-odd density mode BC'd to
   871	### Q1 · The electric-scalar substrate
   872	The **static electric scalar**, closed by a **localized-`H` / PT** construction. ⚠ This is the phase's
   873	entry point onto PHASE 3's throat, ⛔ not onto PHASE 2's brane-shear apparatus.
   874	
   875	**Locus:** `research/pde_ledger_v2/notes/stages/ledger_stage030_electric_scalar_localized_h_closure.md:15`.
   876	
   877	The reduction
   878	
   879	```
   880	M_h = N₀M₄ ,    K_h = N₀K₄ = M_h c_E²
   881	```
   882	
   883	is **EARNED — *given the postulated action***.
   884	
   885	**Expected new:** the localized-`H`/PT closure, and the **throat Green speed `c_E`**.
   886	**Class:** `M_h`, `K_h` — **derived, conditional on the action**; `c_E` — an **interior** quantity
   887	(`parameter_register.md:135`).
   888	**Regime:** static.
   889	**Carry forward:** ⛔ **the action is postulated — a tier-1 item.** ⇒ "EARNED" here means earned *from a
   890	postulate*, not from primitives. ⛔ Do not promote it in transcription. → S22.
   891	**Defect register:** **C6** (*"No closed parent action."*) — ⛔ Q1 does **not** resolve it: `M_h`/`K_h`
   892	are earned *given* a postulated action, which is exactly the closure gap C6 names. ⇒ **Explicitly
   893	deferred**, not closed. → S22.
   894	**Parameter-register edges:** none.
   895	
   896	### Q2 · The puncture deflection and its source ⛔ the holder is a DEBT, not a result
   897	A **±w puncture geometrically bends the brane into ±w**: the field identity `ξ_w = ℓh` and the
   898	orientation-odd mouth source. Token: **`THROAT_H_SOURCE_1_OVER_R2`**.
   899	
   900	**Loci:** `ledger_stage031_puncture_deflection_field_identity_source.md:12`, `:18`, `:25`.
   901	
   902	**Expected new:** the **source identity** and the **far-field FORM** — ⛔ **not a holder.**
   903	**Class:** the **`1/R²` falloff** and the **`s₁s₂` product** are **target-blind EARNED** (`stage031:25`) —
   904	⭐ this is the **FORM** the phase preamble promises, and it is the part that does *not* wait on the
   905	interior. ⚠ EARNED ***within Q1's postulated G0 closure*** (`stage031:18`), ⛔ not from primitives. {#q2-earned-within-g0}
  1174	These are proposed linear tests on supplied backgrounds, not a solved nonlinear throat. S11b supplies
  1175	conditional uniform interface data; S11c's nonuniform response remains partial. Do not infer that all
  1176	five tests are cleared or all depend on a completed S11c loss factor. [S11c closeout, downstream uses](steps/S11c_PARTIAL_CLOSEOUT.md).
  1177	
  1178	| # | what it measures | needs | blocked by |
  1179	|---|---|---|---|
  1180	| 1 | ⭐⭐ **transverse → longitudinal conversion at a defect** — how much light converts to the stray mode passing a particle. ⚠ The families do **not** mix in a *homogeneous* brane (coupling `∝ k·a`); `∇μ_R ≠ 0` mixes them | a defect profile | **response/form remains PARTIAL** (S11c-c2/d); physical **magnitude** needs the interior (`R1`) |
  1181	| 2 | ⭐⭐ **does the longitudinal actually radiate into the bulk** | the interface law | **S11b** |
  1182	| 3 | ⭐ **the flexural crossover** — `ω ∝ k²` at long wavelength once the wall width is dynamical. ⛔ **Able to fail:** a clean cone means move 5 was wrong | width mode inertia + stiffness | **S11b** + **S5–S7** |
  1183	| 4 | **birefringence near a defect** — the two polarisations are degenerate *by symmetry* in a homogeneous brane; a defect splits them | a defect profile | as (1) |
  1184	| 5 | **Test 5** — set `χ_B = 0`, confirm shear propagation dies (`notes/brane_bulk_handoff.md:880`) | nothing new | ⭐ **runnable now**, designed, never run |
    15	Lifts S10's zero to a propagating longitudinal mode, enters the brane's compression modulus as a
    16	**postulated** knob with a **named retirement condition**, closes the homogeneous three-finite-mode
    17	census, and derives the simple \(k_w=0\) grazing locus. The free mode's actual leakage remains an
    18	interface and spectral-boundary problem; its possible finite-amplitude photon-helper role additionally
    19	requires nonlinear coupling and a harmonic/sideband radiation audit.
    20	
    21	⛔ **This is not a defect being repaired.** The extra mode is a characterized departure from a
    22	transverse-only Maxwell field. Historical drum-head discussions used the compressional response as
    23	motivation for a material electric sector, but the current canonical ontology does **not** identify
    24	\(u_L\) with charge: charge is the oriented throat/electric-odd boundary condition, while \(u_L\) is a
    25	separate in-plane longitudinal material mode. ⛔ `FAIL_CAUCHY_STRAY_LONGITUDINAL` is a **misnamed** token;
    26	never read the prefix as a verdict.
   154	| whether matter directly observes the second cone | derived matter/interface coupling |
   155	| whether the mode participates in a finite-amplitude photon helper | nonlinear intensity coupling plus harmonic and sideband response |
   156	
   157	At homogeneous linear order, the mode's bulk leakage and direct observability require the brane–bulk
   158	interface law. Its possible finite-amplitude photon-helper role additionally requires nonlinear intensity
   159	coupling, the complete slab spectrum, and a harmonic and sideband radiation audit.
   160	
   161	⛔ **What S11 does NOT deliver**, stated so it cannot be over-read: bound-versus-resonant classification ·
   162	nonzero interface overlap · an actual leakage rate · unconditional light confinement · observability or
   163	unobservability of the longitudinal branch · a localized or nonradiative nonlinear photon helper.
   164	
   165	## ⭐ What this settles about the light sector as a whole
   166	
   167	Established with the user, 2026-08-02, and load-bearing for how the departure above reads:
   168	
   169	- ⛔⛔ **`c_L` is NOT a gravitational wave.** Gravity in this model is *"the **FLOW** between draining
   170	  defects — carried by the flow + Bernoulli pressure, **NOT** by ripples/radiation"*
   171	  (`docs/conceptual_foundation.md:348`). `c_L` is an in-plane displacement of brane material.
   172	  ⚠ Separately: the corpus contains **no mechanistic account of what a gravitational wave is**.
   173	- ⭐⭐ **A second cone is only a departure if matter COUPLES to it.** The transverse wave equation is
   174	  Lorentz-invariant with invariant speed `c_γ` automatically, so anything built from those modes
   175	  inherits that invariance — and this model's matter **is** built from them (a throat held open by a
   176	  trapped brane-shear standing wave). ⛔ The orchestrator's *"a second cone means observers see no single
    52	## Intermediate results are not a loss factor
    53	
    54	**CONDITIONAL:** the SymPy symbolic/finite omega=1 work reached four full-rank 645-unknown cases and four incident directions with positive regulator and approximate boundaries. Its zero open-thickness selector rank holds at that input, not for all profiles/speeds. Full blind-engine scattering/export obligations remain. [Builder report, four-case completion and retained contract](../_measurements/S11c_d_sympy_builder_report.md); [continuum-current report, opening/limits](../_measurements/S11c_d_continuum_currents_report.md).
    55	
    56	**WITHDRAWN/FAILED:** the omega=1 fixed-point continuation stopped at RIGHT with MemoryError; exit zero was not full success. Its two literal **CLEAR FOR THIS PROGRESS-DEPENDENT FIXED-POINT CONTINUATION** verdicts concerned the build, not results. [Fixed-point report, failure and review scope](../_measurements/S11c_d_transverse_face_fixed_point_continue_report.md).
    57	
    58	**CONDITIONAL:** finite packet/source/matching evidence includes an inner-kernel bank, local Gaussian actions, polarization-dependent first-order forcing and matched transverse amplitudes. The handoff records positive end-current weights and first-order flux cancellation. The scoped receiving result excludes the recorded C3 real-axis poles; the zero-forcing U_B polarization has zero three-row response only in the stated weighted class, and the full five-field transverse poles remain. None supplies complete physical work, a full field or leakage. Reviews had different scopes: the local action retained Claude **NEEDS REVISION** / Grok **CLEAR** with a tooling repair; later Claude-only flux/receiving verdicts cleared source builds, not independent result execution. [Applicability status, accepted inner/local sections](../_measurements/S11c_d_defect_near_unity_applicability_status.md); [source result](../_measurements/S11c_d_first_order_source_result.txt); [matching result](../_measurements/S11c_d_first_order_transverse_matching_result.txt); [handoff, established results](../../../docs/light_em_investigation_handoff.md); [flux review, limits](../_measurements/S11c_d_first_order_transverse_flux_build_claude.md); [receiving review, limits](../_measurements/S11c_d_first_order_receiving_regular_build_r3_claude.md); [receiving method, §3 class/U_B restriction](../_measurements/S11c_d_first_order_receiving_regular_20261004_r2_plan.txt).
    59	
    60	## Clean-condition result and its owner
    61	
    62	**CONDITIONAL:** the clean-condition packet measured a symmetry selection rule on one engine, SymPy: `MIXED 0` and `WRONG_SIDE 0` in all four tested cases, with responsive controls. This was review-leg evidence, not a blind second-engine build. Claude literally said **“not clear.”**, Grok **“not cleared.”** for the broader packet. [Clean-condition disposition, round 5, verdicts and measurements](../directives/_measurements/S11c_d_clean_condition_review_disposition.md).
    63	
    64	**OPEN — S12:** its no-leak conclusion needs the drain frozen and the relevant receiving channel empty. Nothing transfers automatically to live order conversion. The native drain/return functions and their boundary data belong to S12; trapped support belongs to Q2/S22. [Clean-condition disposition, R5-1/R5-2 and final question](../directives/_measurements/S11c_d_clean_condition_review_disposition.md); [plan, S12/Q2](../V3_STEP_PLAN.md).
``````

</details>

<a id="p4"></a>
<details>
<summary>P4 — repair source, command and literal output</summary>

Command:

``````bash
git show ede8aa21:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md | nl -ba | sed -n '47,60p;101,114p'
``````

Exit code: 0.

Literal output:

``````text
    47	  ```
    48	  c_γ(x)² ≡ μ_⊥(x)/ρ_br(x) .
    49	  ```
    50	
    51	  Here `μ_⊥` is the basis-invariant transverse stiffness (`steps/S11b_interface_coupling_law.md:74–87`);
    52	  S9's basis writes it as `μ_R` (`steps/S9_light_requires_shear.md:78–79`). Both `μ_⊥(x)` and `ρ_br(x)`
    53	  remain live. The same `ρ_br(x)` appears in the mass balance below. There is one isotropic speed, the same
    54	  for both polarizations, measured relative to the shear-carrying material; the supplied polarization
    55	  identification is `c_γ,1(x) ≡ c_γ,2(x) ≡ c_γ(x)`.
    56	- **Advection (supplied).** `V(x)` is the in-plane velocity of the material that `u` displaces. Its
    57	  supplied kinetic symbol is
    58	
    59	  ```
    60	  K_adv(ω,k;x) ≡ (ω − V^i(x) k_i)² .
   101	  ```
   102	
   103	  `χ` is the inverse material map. MATERIAL_ADVECTED is not selected. LAB_HELD supplies neither a
   104	  material-reference evolution law nor a physical holder.
   105	
   106	**Outside this step.** These were narrowed out with the user's approval, and are recorded as open:
   107	- direction-dependent (radial versus tangential) stiffness;
   108	- coupling to thickness or bulk fields;
   109	- polarization-dependent propagation;
   110	- mixed `ωk` content other than advection.
   111	
   112	**Reference speed (supplied).** `c₀ ≡ lim_{r→∞} c_γ(r)`, identified with the measured light speed. The bulk
   113	sound speed `c_s` is separate and is not identified with `c₀`. Local light speed is not an observable here,
   114	because rulers and clocks are made of the same medium. Only the far-field observables below are compared.
``````

</details>

<a id="p5"></a>
<details>
<summary>P5 — repair source, command and literal output</summary>

Command:

``````bash
nl -ba research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | sed -n '128,148p'; nl -ba research/pde_ledger_v3/steps/S10_two_transverse_photons.md | sed -n '172,176p;245,252p'
``````

Exit code: 0.

Literal output:

``````text
   128	### R-S1-01 — the brane's spatial dimension
   129	
   130	- **source** S10 (`steps/S10_two_transverse_photons.md`); S11
   131	  (`steps/S11_stray_longitudinal.md`, moves 2–3 / finite census) · **target** S6 · **status** OPEN
   132	- **requirement** — the brane's spatial dimension `D_brane`, as a derived quantity or an explicitly
   133	  re-affirmed postulate.
   134	- **on failure** — S10's headline reads *"light having exactly two polarisations is a statement that our
   135	  space is three-dimensional."* Read backwards it says the opposite: `D_brane = 3` went in and `D−1 = 2`
   136	  came out. Without a delivered `D_brane` the sentence is an assumption restated, ⛔ not a result.
   137	- **note** — ⚠ **The target is S6, not S1.** S1 owns `D = 4` and the two-phase split; `D_brane` is a
   138	  property of the **wall**, and S5–S6 are what construct it — S6 (the kink) is the first step at which a
   139	  codimension-1 surface exists to have a dimension.
   140	  ⛔⛔ **An earlier draft of this entry ruled out the obvious route on bad reasoning**, and a leg caught
   141	  it. It argued that because S10's *mode-count computation* contains no codimension, `D_bulk − 1` is not
   142	  the route. ⚠ **Non sequitur.** That the count does not *use* a codimension says nothing about whether
   143	  `D_brane` can be *derived* as one elsewhere — a domain wall in a `D = 4` bulk is exactly codimension 1,
   144	  and that is the natural derivation. ⇒ ⭐ **this entry does not prescribe a route**, per the schema; it
   145	  asks for the object.
   146	- **pass-2 consumer** — S11's selected homogeneous three-mode census and transverse/longitudinal
   147	  separation use `D_brane = 3`. Its proper-rotation invariant count explicitly warns that the sector
   148	  separation is not dimension-independent. This is the same dimension obligation, not a new knob.
   172	The physical selection D = 3 is not made in S10. The live S10 computation keeps
   173	D symbolic for dimensions and evaluates an indexed sweep at D = 2, 3, 4, 5.
   174	The new Lean baseline proof establishes the conditional map D ↦ D − 1 for
   175	arbitrary finite D with a nonzero wavevector, extending the measured sweep.
   176	Neither establishes which D nature selects.
   245	The zero root retains a degree of freedom; it does not remove one. What the
   246	curl-only stiffness removes is the restoring stiffness for the longitudinal
   247	direction. The per-root nullities make the distinction explicit:
   248	`1 + (D − 1) = D` in every MAIN dimension measured, so all D amplitude
   249	directions remain in the spectrum. Within the light sector, Maxwell has no
   250	counterpart to this surviving zero-frequency longitudinal direction. That is
   251	a characterised departure, and its onward disposition belongs to S11
   252	(`stray_longitudinal`). S10 assigns it no further interpretation.
``````

</details>

<a id="p6"></a>
<details>
<summary>P6 — repair source, command and literal output</summary>

Command:

``````bash
nl -ba research/pde_ledger_v3/V3_STEP_PLAN.md | sed -n '453,464p;871,905p'
``````

Exit code: 0.

Literal output:

``````text
   453	⛔⛔ **DIFFERENT OBJECTS — ⛔ do not weld them, at any strength.** `u_L` is an **in-plane longitudinal
   454	brane displacement**; **`±w`** is the **throat's normal/orientation direction**, and the charge the
   455	rethink produced is carried by the **`h`-branon**. The committed sources keep them apart:
   456	
   457	- `research/pde_ledger_v2/notes/stages/ledger_stage031_puncture_deflection_field_identity_source.md#h-a-distinct-from-u-l`
   458	  — *"Define the held mouth datum `h_A = ξ_w|_A/ℓ = P₀H|_A`, `[h_A]=1` **(distinct from `u_L`)**"*.
   459	- `docs/em_analog_next_phase_handoff.md#u-l-clamp-not-the-charge-scalar` — of the `u_L`-clamp candidate:
   460	  *"It is NOT the committed charge scalar (that's the `h`-branon, `h≠u_L`; `h` remains the committed
   461	  mediator for the conditional `1/R²` falloff); `u_L` is a separate charge-odd density mode BC'd to
   462	  relax to 0."*
   463	- `docs/model_map.md#superseded-charge-routes` — the *"leftover-scalar `u_L`-clamp"* is listed among the
   464	  **superseded** charge routes.
   871	### Q1 · The electric-scalar substrate
   872	The **static electric scalar**, closed by a **localized-`H` / PT** construction. ⚠ This is the phase's
   873	entry point onto PHASE 3's throat, ⛔ not onto PHASE 2's brane-shear apparatus.
   874	
   875	**Locus:** `research/pde_ledger_v2/notes/stages/ledger_stage030_electric_scalar_localized_h_closure.md:15`.
   876	
   877	The reduction
   878	
   879	```
   880	M_h = N₀M₄ ,    K_h = N₀K₄ = M_h c_E²
   881	```
   882	
   883	is **EARNED — *given the postulated action***.
   884	
   885	**Expected new:** the localized-`H`/PT closure, and the **throat Green speed `c_E`**.
   886	**Class:** `M_h`, `K_h` — **derived, conditional on the action**; `c_E` — an **interior** quantity
   887	(`parameter_register.md:135`).
   888	**Regime:** static.
   889	**Carry forward:** ⛔ **the action is postulated — a tier-1 item.** ⇒ "EARNED" here means earned *from a
   890	postulate*, not from primitives. ⛔ Do not promote it in transcription. → S22.
   891	**Defect register:** **C6** (*"No closed parent action."*) — ⛔ Q1 does **not** resolve it: `M_h`/`K_h`
   892	are earned *given* a postulated action, which is exactly the closure gap C6 names. ⇒ **Explicitly
   893	deferred**, not closed. → S22.
   894	**Parameter-register edges:** none.
   895	
   896	### Q2 · The puncture deflection and its source ⛔ the holder is a DEBT, not a result
   897	A **±w puncture geometrically bends the brane into ±w**: the field identity `ξ_w = ℓh` and the
   898	orientation-odd mouth source. Token: **`THROAT_H_SOURCE_1_OVER_R2`**.
   899	
   900	**Loci:** `ledger_stage031_puncture_deflection_field_identity_source.md:12`, `:18`, `:25`.
   901	
   902	**Expected new:** the **source identity** and the **far-field FORM** — ⛔ **not a holder.**
   903	**Class:** the **`1/R²` falloff** and the **`s₁s₂` product** are **target-blind EARNED** (`stage031:25`) —
   904	⭐ this is the **FORM** the phase preamble promises, and it is the part that does *not* wait on the
   905	interior. ⚠ EARNED ***within Q1's postulated G0 closure*** (`stage031:18`), ⛔ not from primitives. {#q2-earned-within-g0}
``````

</details>

<a id="p8"></a>
<details>
<summary>P8 — repair source, command and literal output</summary>

Command:

``````bash
nl -ba research/pde_ledger_v3/V3_STEP_PLAN.md | sed -n '1191,1204p'; nl -ba research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | sed -n '263,305p'
``````

Exit code: 0.

Literal output:

``````text
  1191	#### ⭐⭐ "STAYS QUANTIZED IN A PACKET" IS THREE DIFFERENT QUESTIONS — ⛔ do not merge them
  1192	
  1193	⚠ The user's ambition is *"a photon travelling through space, self-sustaining."* Three readings, **three
  1194	different blockers**, and conflating them caused a false alarm that `ħ` was on the critical path:
  1195	
  1196	| reading | needs | status |
  1197	|---|---|---|
  1198	| **a bound state** — a localized mode with a **discrete** frequency spectrum, trapped by a defect | ⭐ a linear wave equation in a **varying** medium (Sturm–Liouville) | ⭐⭐ **reachable, LINEAR.** No `ħ`, no nonlinearity |
  1199	| **a soliton** — one allowed lump, energy fixed by a size–amplitude relation | the nonlinear shear action | ⛔ **C6** — the real wall, and what the user actually wants |
  1200	| **`ħω` quanta** — energy in discrete lumps | field quantization | ⛔ **not a classical PDE sim at all**, at any effort |
  1201	
  1202	⭐⭐ **⇒ Shelving `ħ_model` does NOT block the simulation.** A classical field has continuous energy no
  1203	matter how good the code is; row three was never reachable this way. ⇒ `ħ_model` connects the model to
  1204	quantum mechanics — a **separate project** from watching the brane move. ⛔ Do not re-open it here.
   263	### R-S8-02 — the in-plane and out-of-plane sectors must decouple at quadratic order
   264	
   265	- **source** S10; S11 (`steps/S11_stray_longitudinal.md:278–280`, in-plane, frozen-wall-width spectrum)
   266	  · **target** S8 · **status** OPEN
   267	- **requirement** — the quadratic operator on the brane's **full** displacement, in-plane `u` **and**
   268	  out-of-plane `h` together, and whether it is block-diagonal in that split.
   269	- **on failure** — ⛔ **S10's headline number changes** — but ⚠ **not by the mechanism an earlier draft of
   270	  this entry named, and a leg corrected it.** Mixing is **not** the failure mode: under the rotation group
   271	  that acts on the brane, `h` and `u_L` are both **scalars** and `u_T` is a **vector**, so `h` can mix with
   272	  `u_L` freely and the transverse count stays `D − 1` regardless. ⭐ **The condition that gives 3 is `h`
   273	  being DEGENERATE WITH THE TRANSVERSE PAIR** — the brane's elasticity being isotropic in `D+1` rather
   274	  than in `D`. That is what *"belonging to the same elastic sector"* has to mean, and it is what S8 must
   275	  settle.
   276	- **note** — ⭐ **This is the most load-bearing open item in the sector**, and it was flagged inside S10 at
   277	  the moment the identification was made rather than found later: `h ≠ u_L`, stated in three places in the
   278	  corpus, **user-confirmed as the picture** — but ⛔ confirmed as a picture, not computed. `V3_STEP_PLAN`
   279	  puts *"transverse and longitudinal sectors; the reduced `h`/`u_L` operator"* in S8, so the object is
   280	  already scheduled.
   281	  ⚠⚠ **And v3 is not starting from nothing — a leg found the computation already exists in v2 and this
   282	  entry did not cite it.** `research/pde_ledger_v2/paper/stages/stage_030.tex:99-125` builds the coupled
   283	  `(u_L, h)` scalar block with stiffness `K = [[B_eff, C_hu], [C_hu, K_h]]` — explicitly **not**
   284	  block-diagonal, with `C_hu` a registered free-unreduced parameter — and `R79` records that the mixed
   285	  poles are cone-coincident only when `C_hu = 0`. ⇒ ⛔ *"nothing computes that"* is true of **v3** and
   286	  misleading about the **corpus**; S8 should start from `stage_030`, ⛔ not from scratch.
   287	  ⭐ Note this is the `(u_L, h)` block — the scalar sector — so it bears on the **longitudinal** slot and
   288	  the charge anchor, ⛔ and not directly on the transverse count, which is what the corrected on-failure
   289	  above turns on.
   290	- **pass-2 source** — S11's result here is the one its record names as the in-plane, frozen-wall-width
   291	  spectrum with `hBranon` excluded: *"⚠ Scope caveats, none touching the decoupling: (i) this is the
   292	  **in-plane, frozen-wall-width** spectrum (`WALL_WIDTH_FIELDS={}`, `hBranon` excluded,
   293	  `INTERFACE_EQUATIONS_SUPPLIED={}`) — it does not decide inhomogeneous mode conversion,"*
   294	  (`steps/S11_stray_longitudinal.md:278–280`). S11's decoupling statement, with its own limit, reads:
   295	  *"Thus the homogeneous \(D=3\) quadratic transverse and longitudinal eigenbranches have zero linear
   296	  cross-block under the selected action. This does not establish decoupling at nonlinear order, on a
   297	  nonuniform slab, at an interface or defect, or after additional allowed fields are introduced."*
   298	  (`steps/S11_stray_longitudinal.md:74–76`). This entry's object is the joint `u`/`h` operator that the
   299	  spectrum leaves out.
   300	- **excluded field** — in that caveat S11's record freezes the wall width and excludes `hBranon`
   301	  (`steps/S11_stray_longitudinal.md:279–280`); its Mathematica audit binds the token as
   302	  `"OUT_OF_PLANE_FIELD_EXCLUDED" -> hBranon`
   303	  (`mathematica/S11_stray_longitudinal_mathematica_audit.wl:1082`). If S8 does not deliver the joint
   304	  operator, S11's spectrum is not extended beyond that stated scope. The record does not say what its
   305	  roots or census would be with `hBranon` included, and this register does not infer it.
``````

</details>


<a id="q4"></a>
<details>
<summary>Q4 — Repair 2 sources, command and literal output</summary>

Command:

``````bash
nl -ba research/pde_ledger_v3/steps/S11_stray_longitudinal.md | sed -n '52,64p'; nl -ba research/pde_ledger_v3/steps/O2_steady_brane_balance.md | sed -n '89,100p;109,123p'; nl -ba research/pde_ledger_v3/steps/S11bB_interface_assembly.md | sed -n '53,63p'; git show ede8aa21:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md | nl -ba | sed -n '117,119p;34,44p;165,186p'
``````

Exit code: 0.

Literal output:

``````text
    52	N_SO = {D=2 → 4, D=3 → 3, D=4 → 4, D=5 → 3}        N_O = 3 for every D
    53	```
    54	
    55	The extras are **reflection-odd**: `(tr G)·ε^{ij}G_{ij}` at `D=2`, `ε_{ijkl}G_{ij}G_{kl}` at `D=4`.
    56	⚠ **And the walk's METHOD was the error, not its arithmetic** — *"decompose into pieces that do not mix
    57	and count them"* misses **cross-pairings between isomorphic summands**, which is exactly what `D=2` and
    58	`D=4` have (at `D=2` the trace and `Λ²` are both `SO(2)` scalars; at `D=4`, `Λ²` splits into self-dual ⊕
    59	anti-self-dual). ⭐ Five independent derivations agree on the corrected counts.
    60	
    61	⭐⭐ **The physical content of the extras, and it sharpens the step:** at `D=2` the extra invariant is
    62	**NOT a total derivative** — its Euler–Lagrange operator is non-zero — so an `SO(2)`-invariant Lagrangian
    63	**can mix the longitudinal and transverse sectors**. At `D=4` it **is** a total derivative and adds no
    64	operator structure. ⇒ ⭐ **Sector separation is a `D=3` fact, ⛔ not a structural one.**
    89	exchange, loads, energy and power; optical `μ_⊥` uses the same measure as `ρ_br`. Induced metric and
    90	native area factors remain explicit. The supplied geometric/optical/mass inputs are
    91	`g_ij=δ_ij+∂_iξ_w∂_jξ_w`, its inverse, `ξ_w=ℓh`,
    92	`c_γ²≡μ_⊥/ρ_br`, `c_γ=c₀(1+δ)` and `∇·(ρ_br V)=−j_n`.
    93	LAB_HELD is supplied spatial speed anchoring, not a holder or reference-evolution law. The bulk
    94	inputs `P=Kρ^n`, `c_s²=nKρ^(n−1)/m`, `f=ρ/ρ₀−1` use bulk number density, particle mass and symbolic
    95	EOS exponent. They supply no bulk profile, brane-density response or live traction. These are inputs
    96	on their recorded domains, not results earned by engine agreement. Sources: spec §§1, 3.1, 7–8;
    97	contract §§5–7 (M1).
    98	
    99	Each premise below is exactly labelled **adopted substrate input to a conditional model (2026-10-06)**,
   100	from the user-selected premise decision list 1–4. None is derived (M1).
   109	`V`, `ρ_br`, `μ_⊥`, `ξ_w`, `h`, `δ`, `j_n`, bulk state, material responses, history/reference,
   110	loads and boundary data remain live with spatial derivatives. Radial profiles impose no constitutive
   111	isotropy, parity, stress symmetry, absent couple/chiral content, derivative cutoff or finite history
   112	state. No sharp-sheet or finite-slab material reduction, stress measure, relaxation family, passive
   113	sign or physical support is selected. Sources: spec §§1–6; contract §§1, 3–8 (M1).
   114	
   115	Supplied counting is `ε=GM/(c₀²r)`, `δ=O(ε)`, `(∂ξ_w)²=O(ε)`, `V/c₀=O(ε^{1/2})`,
   116	`(V/c₀)²=O(ε)` and the optical monomial box `0≤a≤1, 0≤b≤2, 0≤c≤1` in
   117	`δ^a(V/c₀)^b((∂ξ_w)²)^c`, with nonnegative integer indices. Only the stiffness/density ratio
   118	inherits the speed-change grade. Individual density/modulus, stress/inertia/normal response,
   119	exchange/source/load, relaxation/power, holder/embedding and derivative grades remain OPEN. O2 is
   120	untruncated; this optical box removes no mechanical or energy term. Bulk first order in `f` is
   121	separate, with no supplied `f`–`ε` relation. The mass law carries its recorded relative-`O(ε)`
   122	qualification for claims transferring `j_n` or `ρ_br` to induced measure; no induced-measure mass
   123	law is supplied or derived. Sources: spec §7; contract §9 (M1).
    53	⚠ Non-reciprocal, "odd" constitutive couplings of exactly the excluded form are realised in driven
    54	laboratory media (odd elasticity, odd viscosity, active and chiral fluids). ⛔ They are not unphysical.
    55	⭐ **They have a reservoir.**
    56	
    57	### ⭐⭐ THE STANDING RULE THIS REPLACES IT WITH
    58	
    59	⛔ *"Anything goes"* is not the correction — under it no mode is ever stable and the ledger predicts nothing.
    60	
    61	> ⭐⭐ **A non-passive coupling is admissible only with a NAMED reservoir and a STATED power budget.**
    62	> Unbounded growth fed by nothing is still a defect. Unbounded growth fed by `v₀` is **physics** — and it is
    63	> **quantitative**, because `v₀` is bounded by cosmology.
    34	## Governing object (supplied; the computation cannot test the input)
    35	
    36	Light is the transverse branch of the brane's in-plane displacement `u`, and `u` is the material displacement
    37	of the stuff whose density is `ρ_br` (`steps/S11_stray_longitudinal.md:32–35`; S10). Near the mass, at
    38	leading eikonal order, that branch obeys the supplied local dispersion relation
    39	
    40	```
    41	(ω − V^i k_i)² = c_γ(x)² · g^{ij}(x) k_i k_j ,      g_ij = δ_ij + ∂_iξ_w ∂_jξ_w .
    42	```
    43	
    44	Each piece is a supplied identification. Flag any result that depends on one.
   117	- **Eikonal.** The retained object is the dispersion relation above, with its position-dependent
   118	  coefficients, and its rays. Excluded from this step's claim and not computed: the explicit subprincipal
   119	  terms of the underlying operator, which affect amplitude and polarization transport.
   165	## The observables
   166	
   167	**Quantifier.** The Part B conditions to compute are the profile conditions for matching the reference
   168	for **every** far-zone `b`. Part D uses that same quantifier. No matching outcome is supplied.
   169	
   170	**Setting.** One isolated, spherically symmetric mass at rest, with its drain flowing: steady, not frozen,
   171	with `V` live. Far field; linear waves. Time is the lab time of the brane's far-field rest frame.
   172	1. **Δθ(b):** the total turning angle of a full flyby with impact parameter `b`.
   173	2. **Round-trip (radar) excess time.** An emitter at distance `Z_E` on one side of the mass and a reflector
   174	   at `Z_R` on the other, along the line. Subtract the flat round-trip time. Here
   175	   `r_E = √(b² + Z_E²)` and `r_R = √(b² + Z_R²)`. Part A prints the full excess time `Δt_RT(b; Z_E, Z_R)`.
   176	   The Part B comparison uses only its **logarithmic slope** in `b`, taken at fixed endpoints, for far endpoints:
   177	
   178	   ```
   179	   𝒮_RT(b) ≡ lim_{Z_E, Z_R → ∞} ∂Δt_RT(b; Z_E, Z_R) / ∂ln(1/b²) ,      Z_E and Z_R held fixed in the derivative.
   180	   ```
   181	
   182	   This applies to every profile, including tails other than `1/r`. Where `Δt_RT` is `A·ln(1/b²)`, plus terms
   183	   independent of `b`, plus terms that vanish with their `b`-derivatives as `Z_E, Z_R → ∞`, `𝒮_RT = A`. The
   184	   radar claim is limited to this slope. Terms independent of `b`, which carry the endpoints, are not compared.
   185	3. **The two one-way excess times** between the same endpoints, and their **nonreciprocal part** (half the
   186	   difference). Print whether that part depends on the path or only on the endpoints.
``````

</details>



<a id="q5"></a>
<details>
<summary>Q5 — Repair 3 sources, command and literal output</summary>

Command:

``````bash
nl -ba research/pde_ledger_v3/steps/S10_two_transverse_photons.md | sed -n '73,86p;172,197p'; nl -ba research/pde_ledger_v3/steps/S11_stray_longitudinal.md | sed -n '40,64p'; nl -ba research/pde_ledger_v3/steps/S11bB_interface_assembly.md | sed -n '74,76p'; nl -ba research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | sed -n '412,427p'
``````

Exit code: 0.

Literal output:

``````text
    73	> For the supplied curl-only in-plane action, at nonzero wavevector and away
    74	> from allowed exceptional strata, the nonzero root has D − 1 transverse null
    75	> directions in every measured MAIN case D = 2, 3, 4, 5. Thus the D = 3 member
    76	> has two transverse directions.
    77	
    78	The form controls prohibit a stronger reading of that sentence. At both
    79	dimensions where the form-changing packages and MAIN were emitted, D = 3 and
    80	D = 4, FULLGRAD has the same nonzero root and the same D − 1 transverse
    81	nullity as MAIN. The transverse count therefore cannot be attributed to the
    82	curl-only form. What that form determines here is the disposition of the one
    83	remaining, longitudinal direction: unlike FULLGRAD, it leaves that direction
    84	at the zero root instead of putting it on the propagating root. DIVONLY, which
    85	reverses the two sectors, confirms that the stiffness form can change the
    86	count even though curl-only is not what uniquely produces D − 1.
   172	The physical selection D = 3 is not made in S10. The live S10 computation keeps
   173	D symbolic for dimensions and evaluates an indexed sweep at D = 2, 3, 4, 5.
   174	The new Lean baseline proof establishes the conditional map D ↦ D − 1 for
   175	arbitrary finite D with a nonzero wavevector, extending the measured sweep.
   176	Neither establishes which D nature selects.
   177	
   178	## What was supplied and what was computed
   179	
   180	The shared specification supplies:
   181	
   182	- a D-component in-plane displacement u, with its separation from every
   183	  out-of-plane field inherited rather than tested;
   184	- the real cosine plane-wave ansatz;
   185	- positive inertia and stiffness coefficients, nonzero real wavevector,
   186	  unstrained rest background, no dissipation, and linear response;
   187	- the curl-only stiffness density
   188	
   189	      S_curl = (1/2) Σ_i Σ_j (∂_i u_j − ∂_j u_i)²
   190	
   191	  in
   192	
   193	      L = (ρ_br/2) Σ_j (∂_t u_j)² − (μ_R/2) S_curl.
   194	
   195	These premises and the action are at
   196	directives/S10_SHARED_PHYSICS.md:13-28, :30-47, and :82-107. The action is an
   197	input, not a result of the mode count.
    40	**Move 2 — what a quadratic stiffness on the brane can be.** Decompose `∂_i u_j` into trace,
    41	symmetric-traceless and antisymmetric parts. ⭐ **S10 did not forget a term — it kept ONE invariant of
    42	several.** Curl-only stiffness charges for twist and nothing else, which is why the longitudinal came out
    43	free rather than forbidden: a longitudinal plane wave has a symmetric gradient
    44	and hence no curl. That gradient generally contains both trace and symmetric-traceless
    45	parts; it is not pure trace. The selected curl-only action charges neither of those
    46	parts. (Clarified during the S11 homogeneous Lean fidelity pass, 2026-09-16 UTC.)
    47	
    48	⛔⛔ **CORRECTION — the orchestrator's count in this move was WRONG, and both engines caught it.** The
    49	walk asserted *"exactly three coefficients, for every `D ≥ 2`."* ⭐ **False under proper rotations.**
    50	
    51	```
    52	N_SO = {D=2 → 4, D=3 → 3, D=4 → 4, D=5 → 3}        N_O = 3 for every D
    53	```
    54	
    55	The extras are **reflection-odd**: `(tr G)·ε^{ij}G_{ij}` at `D=2`, `ε_{ijkl}G_{ij}G_{kl}` at `D=4`.
    56	⚠ **And the walk's METHOD was the error, not its arithmetic** — *"decompose into pieces that do not mix
    57	and count them"* misses **cross-pairings between isomorphic summands**, which is exactly what `D=2` and
    58	`D=4` have (at `D=2` the trace and `Λ²` are both `SO(2)` scalars; at `D=4`, `Λ²` splits into self-dual ⊕
    59	anti-self-dual). ⭐ Five independent derivations agree on the corrected counts.
    60	
    61	⭐⭐ **The physical content of the extras, and it sharpens the step:** at `D=2` the extra invariant is
    62	**NOT a total derivative** — its Euler–Lagrange operator is non-zero — so an `SO(2)`-invariant Lagrangian
    63	**can mix the longitudinal and transverse sectors**. At `D=4` it **is** a total derivative and adds no
    64	operator structure. ⇒ ⭐ **Sector separation is a `D=3` fact, ⛔ not a structural one.**
    74	**The transverse mode is completely decoupled on a uniform background.** The coupling is **identically
    75	zero**, the dispersion is `ρ_br⁰ω² = μ_R k²`, and the imaginary part is **zero**. Both engines, same
    76	structural reason: in-plane parity admits no `e_W ↔ u_T` bilinear.
   412	### R-S8-06 — material displacement and the slab's quadratic inertia
   413	
   414	- **source** O2 (material identification only; `steps/O2_steady_brane_balance.md`, register handoff); S11 (`steps/S11_stray_longitudinal.md`, move 1 / finite census); S11b-B
   415	  (`steps/S11bB_interface_assembly.md`, breathing quadratic) · **target** S8
   416	  (**register inference for the original S11/S11b records; O2 names S8 for the material identification**) · **status** OPEN
   417	- **requirement** — `u` as the material displacement of the stuff whose density is `ρ_br`, and the
   418	  quadratic kinetic form of that material and the thickness degree of freedom (B's `μ_W`) on the
   419	  original homogeneous/slab domain. O2 is a consumer of the material identification only.
   420	- **on failure** — if `u` is a director rather than that displacement, S11's continuity identification
   421	  `δρ_br = −ρ_br ∇·u` and its compression argument do not follow. Its finite census requires nonzero
   422	  `ρ_br`; B's breathing-root interpretation also rests on the stated inertial model. A stiffness
   423	  functional alone (`R-S8-01`) does not supply this identification or the kinetic form.
   424	- **note** — this asks for the field identity and inertia, not numerical benchmark values of `ρ_br`
   425	  or `μ_W`, and does not identify the thickness mode with S10's out-of-plane displacement. The cited
   426	  records name no future owner for this object; S8 is a register inference
   427	  (`steps/S11_stray_longitudinal.md:32–38`; `steps/S11bB_interface_assembly.md:80–85`).
``````

</details>
