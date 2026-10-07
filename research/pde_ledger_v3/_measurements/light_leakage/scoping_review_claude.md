# Light-leakage scoping inventory: fresh Claude leg (opus), 2026-10-07

Saved verbatim from the leg's final report. Prompt: `light_leakage_scoping_review_prompt.md` (identical to Grok's).
The leg's scratch files are in `/tmp/llrev/`.

**Verdict: NEEDS REVISION.** I read the directive and the cited records before the artifact. Every web claim below came from a fetch made in this review. Scratch files are in `/tmp/llrev/`, and nothing in the repo was modified. In the evidence lines, `sed -n Np <file>` means the record line was printed with that command. Record paths are relative to `research/pde_ledger_v3/` and the artifact is `A` (`_scratch/light_leakage/light_leakage_scoping.md`).

## Findings

**1. Part 1, column 4 does not name observations that bound the loss in question.** (A:11–23)
- **What's wrong:** Column 4 must bound the loss of the model's analog, which is light on the brane, i.e. vacuum light. That is why it says "Do not say whether the model must suppress." Instead it lists silica-fiber product specifications, material parameters and theoretical predictions.
- **Evidence from the artifact's own cells:**
  - A:11 "this bounds total line attenuation, not absorption separately"
  - A:18 "the bounded quantity is manufactured cladding diameter, not loss caused by it"
  - A:20 "does not isolate tunnelling"
  - Row 3 gives a Raman gain coefficient and row 4 an SBS threshold. Both are material parameters, not bounds.
- **Evidence that rows 7 and 11 cite theory, not observation:** both cite Marcuse's calculations. `curl` of the Bell Labs page for the 1969 paper returns "An rms deviation of one of the waveguide walls of 9A causes a radiation loss of 10 dB per kilometer (index difference 1 percent, guide width 2.5[μ])". That is a prediction. Putting it in a constraint column blurs the line between an oracle and a bound (M3).
- **Fix:** For each mechanism, name a free-space observation bounding that destination, or write "gap". Candidates:
  - cosmological transparency
  - CMB spectral distortion (FIRAS) for frequency exchange
  - astrophysical or cosmic birefringence limits for polarization coupling
  - vacuum cavity photon lifetime or finesse for per-length loss
  - the grating yardstick for edge conversion

  Move the Marcuse numbers to §2.2.

**2. The destination column is misclassified and contradicts Part 2.** (A:11–23)
- **Tally:** `awk -F'|' 'NR>=11&&NR<=23{print $2,"->",$7}' A` gives exchange of frequency ×4, bulk ×7, same transverse ×2, and **another brane branch ×0**.
- **Row 11:** The analog is the close-then-extract coupling, which c2:7 describes as re-extracting the "off-diagonal transverse↔`{θ,e_W,u_L}` coupling". Yet the row's destination is "same transverse branch (survives)". By the row's own mapping it should be another brane branch.
- **Rows 2, 6, 7, 8, 9:** These cite the S11c-b kernel, whose receivers are θ, e_W and u_L, but say only "bulk". Part 2 routes the same kernel "transverse → `θ,e_W,u_L` brane branches" (A:69).
- **Rows 1 and 13:** Absorption and lossy-cladding overlap deposit energy in the medium or the bulk. They are not a frequency exchange.
- **Row 2:** Elastic redirection within the brane counts as survival under the closeout rule (closeout:3).
- **Fix:** Reassign the destinations, and allow more than one where the analog has several receivers.

**3. Route C2 is misclassified and duplicates E1/E3.** (A:38)
- **Evidence:**
  - A:38 labels the gradient-driven thickness channel "Candidate beyond scope".
  - `sed -n 16p steps/S11c_SCOPE.md` gives "The gradient-driven channel is exactly what S11c computes."
  - The S11c-b and c2 receivers include e_W (c2:7).
- **Why it matters:** C2 is the e_W component of E1/E3, which are equation-level routes. Listing it separately invites double counting in the bracket.
- **Fix:** Merge C2 into E1/E3, or mark it as that component and give the scope line.

**4. Routes identified in the records are missing.**
- **(a) Direct leak into a bulk that carries shear.** S9:67 says "without it a photon has somewhere else to go, light leaks off our 3D space". SUBSTRATE_REQUIREMENTS:151 says "disordered phase carries shear ⇒ light is not confined to the brane, it leaks into the bulk". This is R-S1-02, and its status is OPEN, so it belongs under the third status. The Molz–Beamish oracle (fluid helium does not radiate SH(0); solid helium does) is the natural oracle for this route.
- **(b) The nonlinear DC/harmonic/sideband radiation audit.** S11c_SCOPE:50 reads "NOT S11c — the nonlinear light program: the DC/harmonic/sideband radiation audit and nonlinear…". Related lines: S11:157–159 and :297–304, and S11b:164. This is a candidate beyond scope whose destination is frequency exchange. O1 (A:39) covers mixing into brane fields only.
- **(c) Second-hop radiation from the scalar sector into the bulk.**
  - S11bA:72 reads "propagating `Re Z` survives as `ρ_m/a_q` — that is radiation resistance".
  - S11bB:93–94: for `K₀>0` "every nonzero imaginary part comes from propagating bulk radiation resistance".
  - S11c_SCOPE:17 says the longitudinal mode's "fate needs no gradients (S11b did it)".

  So C1's "kinematic only" framing leaves out the equation-level S11b result.
- **(d) Strong edges and finite contrast (S11c-e)** (closeout:19, :31; d:68), and **nonuniform confinement and rotational admissibility**, which are OPEN (closeout:27).

**5. The per-item v₀ quotes and named freezes are incomplete.** (A:33–41)
- **Missing v₀ quotes:** The directive asks, for each item, what the records say about the drain v₀. E1, C1, C2, O1 and O3 quote nothing about v₀, and E3 refers to E1, which has none.
- **Freezes the records name but the artifact omits:**
  - d:29 "one held profile, strict rest bulk, LAB_HELD/RHO4_CONSTANT"
  - c1:24 "seal 5 (background density) stays a **surfaced rule-17 freeze (UNDECIDED)**", which c2:192–195 re-adjudicates
  - c1:154 "strict-`v_dr=0` rest-frame qualification"
  - c1:194–196, convection dropped under the rest-frame limit
  - S9:349–350 "Limits taken: sharp zero-width sheet · `v₀ → 0` · no dissipation · frequency-independent moduli · continuum limit · amplitude → 0"
  - S11b:156–157, the frozen wall width that E4 inherits
- **Wrong citation:** E1's only listed freeze, "the flat uniform limit", cites S11bB:74–78, and those lines do not call it a freeze. S11c_SCOPE:21–22 does.

**6. E4 misstates the record.**
- **Artifact (A:36):** "the uniform transverse coupling is zero only in the passive uniform slice".
- **Records:** S11bB:191 says "coupling ≡ 0 … unconditional on a uniform background". S11b:75 says "`∂²U/∂u_T ∂e_W` is **identically zero**", and the record calls this unconditional. The zero comes from parity, not passivity.
- **Fix:** Say "zero in the uniform limit", citing :191.

**7. Model mapping and status are wrong in Part 1.**
- **Raman and Brillouin (rows 3–4) are labelled NOT ADDRESSED,** but the records do address them as an assigned open item. S11:158 reads "…additionally requires nonlinear intensity coupling…" and continues into the sideband audit at :159. S11c_SCOPE:50 says the same. These rows should be PRESENT ANALOG, with the brane's longitudinal and thickness branches as the acoustic partner, status OPEN.
- **Absorption (row 1)** should cite the "no dissipation" limit at S9:349.
- **Rows 6, 7 and 9 report CONDITIONAL** while their own justification cites lines marked UNRESOLVED (d:11 and :70, c2:51).
- **Rows 2, 5, 8 and 10 cite c1:99–110,** which marks the kernel ESTABLISHED and the other items UNDECIDED, not CONDITIONAL.
- **Row 12 maps the wrong record.** d:58's "polarization-dependent first-order forcing" is polarization-dependent conversion, not coupling between the two polarizations. The record line that names birefringence is closeout:31, which assigns it to S22 as OPEN.
- **Row 13's "OPEN"** cites S9, which marks the item "LIVE" and postulated. The OPEN mark is at SUBSTRATE_REQUIREMENTS:146.
- **Fix:** Quote each mark with its line.

**8. The Part 1 references were not opened and are not marked UNVERIFIED.** (A:7)
- **Artifact:** "the attenuation/loss chapter is visible in the Google Books record".
- **What I found:** Fetching the Keiser record returned "Only metadata and descriptions are visible—no chapter body text… specific subsection headings about loss… are not displayed." The Snyder–Love page returned "Only metadata and table of contents are visible." The *Optical Fiber Telecommunications* source is a shop page.
- **Why it matters:** The directive says to mark any source that could not be opened UNVERIFIED. The "union of loss headings" also cannot have come from Keiser's loss chapter.
- **Fix:** Mark all three UNVERIFIED, or open the chapters.

**9. The Gu & Fuller oracle quotes the wrong equations and maps them to the wrong route.** (A:59)
- **What was quoted:** Eqs. (31)–(33), which are the active-control cost function and the optimal control force.
- **Evidence:** The VTechWorks PDF abstract says "Feed-forward control is achieved by adding secondary line forces… cost function that integrates the far-field radiated acoustic intensity." The paper's passive result is the far-field pressure radiated when a subsonic wave scatters off a discontinuity (Eqs. 13–16; uncontrolled power `C = ∫s[b b*] ds`, Eq. 32). That passive result is the oracle for E2.
- **The mapping is wrong too:** It pairs "Active forces/controller ↔ E4's required named reservoir". A controller cancels radiation; E4's reservoir powers a non-passive interface. They are not analogous.
- **Fix:** Quote the passive radiation equations and map them to E2.

**10. The yardsticks are incomplete.** (A:89–94)
- **Friedland & Giannotti:**
  - The artifact gives only the supernova bounds, Eq. (21). Those numbers match the arXiv 0709.2164 text.
  - It omits the red-giant bounds, Eq. (19): "≃ 1.4 × 10²¹ TeV (n=1)… 5 × 10⁶ TeV (n=2)… 60 TeV (n=3)", which are the strongest. It also omits the horizontal-branch bounds, Eq. (22).
  - It omits the generic bounds on anomalous energy loss that the paper applies. Any invisible-energy route needs these, not a parameter of their model:
    - red giants: "not exceed about twice the standard neutrino luminosity"
    - supernova: "QMax = 3 × 10³³ erg cm⁻³ s⁻¹"
    - horizontal branch: "10 erg g⁻¹ s⁻¹"
  - It should also note that the paper's bounds apply to plasmons, i.e. photons that are massive inside a medium.
- **Tired light:** Only time dilation is given. The Tolman surface-brightness test and the CMB blackbody constraint are missing. The tired-light yardstick is assigned to no route; O1 gets "gap" (A:81).
- **Also missing:**
  - FIRAS spectral distortion, for the frequency-exchange destination
  - a lab photon-lifetime bound for vacuum loss per unit length
  - particle stability for the trapped brane-shear mode at the throat (S9:70–73; SUBSTRATE_REQUIREMENTS:151–153). Leave it out only if scope is explicitly deferred.
- **Undefined abbreviations:** The matrix uses BO, CT, SE and MD, but §2.4 never defines them.

**11. Minor issues**
- **Correlations column (A:69–83):** mostly uncited and generic, e.g. "jointly selected inputs". It omits the correlations the records actually note:
  - the passive-region and reciprocity relations (S11bB:23, :190)
  - the `σ_W` binding (SUBSTRATE_REQUIREMENTS:560; c1:165)
  - `K₀>0` (S11bB:192)
  - `μ_⊥≥0` (SUBSTRATE_REQUIREMENTS:316–319)
  - the density–shape term at order εη (c1:136–138)
- **E2 validity:** The records limit the first-shape coefficient to non-grazing conditions, and exact grazing is recorded as "NOT_ESTABLISHED_AT_FIRST_SHAPE_ORDER" (c1:146–155). The evanescent completion is a second-shape term (c1:191–193). The artifact states neither.
- **Wrong label (A:71):** "reactive/anti-Hermitian". c1:189 says "the spec's 'anti-Hermitian' label for `K_a` is wrong — fix" wherever it is reused.
- **Molz–Beamish link:** It points to a ResearchGate page for "Electrorheological Properties of Chiral Smectic-C Liquid Crystals…", a different paper. The content is nonetheless correct: the Crossref abstract says "For the SH(0) mode, only shear waves are generated in the helium but the L(1) mode can radiate both transverse and longitudinal waves." Replace the link.
- **Marcuse 1970:** The full text is open on archive.org (`bstj49-2-273`): "abrupt steps of a 2:1 ratio… can be kept below 1 percent". It can serve as an anchor value even if the equations are left UNVERIFIED.

## Checked and found sound
- The directive's seeds, the third-status placements (O1/O2/O3) and the counts (13 mechanisms; ABSENT 1, PRESENT ANALOG 9, NOT ADDRESSED 3; 9 routes) are all present.
- The rest-versus-driven tension is quoted, not resolved (A:43 matches S11bB:46–48 and SUBSTRATE_REQUIREMENTS:199–212).
- The survival rule is applied consistently in the matrix, and statuses are input assignments per SUBSTRATE_REQUIREMENTS:543–570.
- No verdicts or suppression requirements appear, and §2.2 disclaims any prior art as a v3 bound.
- Record lines confirmed: S11bA:79 (the Re Z formula), S11:95–101, d:60–64, SUBSTRATE_REQUIREMENTS:251–282, S11bB:61–68 and :150–154.
- Sources confirmed:
  - Marcuse 1969: abstract numbers, and the structure of Eqs. (72)–(76) in the archive.org OCR
  - Kubrusly Eqs. (1), (8), (14)–(16)
  - Peyton Eq. (2) with C = 1.38
  - Demma abstract
  - Molz–Beamish content
  - Friedland & Giannotti Eq. (21)
  - More, Bovy & Hogg "Δτ < 0.13 at 95% confidence"
  - DES "b=1.003±0.005 (stat) ±0.010 (sys)"
  - Optics Letters 34:1708 "(99.577±0.002)%", measured at a single point
  - PRX 13:011018: Q up to 10⁸, and 7 μHz. These are two different devices; the 7 μHz figure is for the torsion pendulum.
- I could not open the Corning specification sheet (HTTP 403), so the fiber spec numbers are unchecked.
