# Light-leakage scoping inventory r1: fresh Claude leg (opus), 2026-10-07

Saved verbatim from the leg's final report. Prompt: `light_leakage_scoping_review_prompt.md` (identical to Grok's).
Reviewed v1 sha256 `d6101683…aaa26` (`scoping_v1.sha256`). Leg scratch: `/tmp/llrev_r1leg.fgRA`.

**Verdict: NEEDS REVISION**

I reviewed `_scratch/light_leakage/light_leakage_scoping.md` at sha256 `d6101683…aaa26`, which matches `scoping_v1.sha256`. I read the directive and the records before the artifact, and I did not read any other review output. I modified no repository files. My scratch files are in `/tmp/llrev_r1leg.fgRA`.

One caveat: `/tmp/llrev` already held another session's downloads. My first curl wrote `fg.pdf`, `berceau.pdf` and `m70.txt` there before I noticed, and the overwrite of `m70.txt` replaced a file with HTML. I then moved to a private directory and read none of that session's files.

---

## Findings

**F1. The governing S11c drain decision and the definition of survival are missing.** This is substantive.
- **Where:** §2.1, the E1 and E3 drain cells, which say "b:54–57 says nothing about `v₀`" and "c2:63–70 says nothing separately about `v₀`". Also the §2.3 Validity column for E1–E3 and C4, and the definition of **S** ("energy in another brane branch or outgoing bulk does not [survive]"), which is attributed to PC:3.
- **Evidence:**
  - `awk 'NR>=179&&NR<=186' directives/S11c_decisions.md` returns: "S11c **inherits it as a standing rest-frame limit** … **every** S11c spectrum/leakage/confinement result is **conditional on the derived smallness domain** (`|q v_bulk_normal_0/ω| ≪ 1` …)".
  - At :130–133 (N13): "Conversion into a **bound** breathing/thickness mode kills the photon exactly as bulk radiation does".
  - PC:3 says only "Both reflected and transmitted light count as survival."
  - `grep -n "decisions\|N11\|N13\|smallness"` on the artifact returns nothing.
  - c1:195–196, which the artifact does cite, itself says "inheritance of N11a's standing rest-frame limit".
- **Fix:** Quote N11(a) as the drain freeze and as the validity domain for every item derived from S11c (E1–E3, C4, and O5's d evidence). Cite N13 for the half of **S** that says what does not survive.

**F2. Bends are mapped to an analog the records do not contain.** This is substantive.
- **Where:** Row 5 (macrobending) is mapped to "curved faces … (c1:35–44)". Row 6 (microbending, "small bends") and the face-odd part of row 7 (roughness) are mapped the same way.
- **Evidence:**
  - `awk 'NR>=235&&NR<=239' directives/S11c_a_SHARED_PHYSICS.md` gives the background as "`Q_bg ∈ {W_bg, μ_R,bg, ρ_4D,bg⁰, ρ_br,bg⁰}`". There is no background centre-line.
  - c1:181–182 says the kernel's shape displacement is "the face-EVEN outward displacement `a_s=(W_bg−W₀)/2`".
  - `ζ_c` appears only as a wave perturbation (S11c_a spec :41, :194).
  - The records do name brane bending: S9:68 says "**Stressed hardest at the throat**, where the brane is bent into `±w`" (also SR:157).
  - The Marcuse 1969 OCR distinguishes the two cases: "If a = 0, that is if the width of the guide changes sinusoidally, only even modes can be excited while sinusoidal deviations from straightness (a = π) couple the even fundamental mode only to odd spurious modes."
- **Fix:**
  - Map macro- and microbending, and face-odd roughness, to NOT ADDRESSED, with S9:68 and SR:157 as the record pointers.
  - Name the flat background centre-line as a freeze.
  - Add a §2.1 gap or open item for conversion at a bend.
  - List Marcuse's width-versus-straightness parity rule as an oracle.

**F3. The static background is never named as a freeze, so frequency exchange at linear order is missing.** This is substantive.
- **Where:** Rows 3 and 4 describe Raman and Brillouin only as "nonlinear … intensity dependence". "Frequency exchange" is routed only to C3, the nonlinear program.
- **Evidence:**
  - The S11c_a spec at :238 defines "`LAB_HELD: Q_bg^L(x,t) ≡ Q_bg(x)`", which is time-independent, so every equation-level route conserves `ω`.
  - `V3_STEP_PLAN.md:566–570` says "the wave is a perturbation on a background that is **not static** … hold the background fixed over a wave period … record `wave period ≪ leak timescale` as a stated validity condition". The artifact never quotes this.
  - Spontaneous Raman and Brillouin scattering are linear-in-field scattering off time-varying fluctuations. Only the stimulated forms are nonlinear. Wikipedia's Raman article (raw text, :73–77) reads "spontaneous Raman scattering … Stimulated Raman scattering is a nonlinear optical effect", and its Brillouin article (:24) defines SBS for "intense beams".
- **Fix:**
  - Quote the plan's hold as a freeze on E1–E5 and C4.
  - Add a NOT ADDRESSED or gap item for frequency exchange off a time-dependent background (moving or fluctuating profiles), with the FD and TL yardsticks.
  - Correct the physical-condition text in rows 3 and 4 to include the spontaneous processes.

**F4. Verified prior art that the records point to is missing from the oracle list for O4 and E4.**
- **Where:** §2.2 says Molz is the only O4 candidate (UNVERIFIED) and that E4 has a "gap".
- **Evidence:**
  - Friedland–Giannotti Eqs. (7)–(9), from `pdftotext` of arXiv:0709.2164: "Γ₀^vac = c_n m(m/k)^n", with c_n = πn/(2^{n+1}Γ[n/2+1]²). This is the tunnelling rate for a photon bound to a brane escaping into the bulk, which is the closest analog of S9:67, "light leaks off our 3D space". The artifact uses this paper only as a yardstick.
  - B:53–55 names prior art for E4's non-reciprocal driven coupling: "odd elasticity, odd viscosity, active and chiral fluids".
- **Fix:** Add both as oracle candidates, with their conditions (warped localization by gravity and a zero mode lifted by a mass; a driven medium). Neither is a bound.

**F5. The V0 yardstick says "no sourced observable", but a record names one.**
- **Evidence:** `V3_STEP_PLAN.md:572–573` says it "ties the DC leak rate to an **observable** (the expansion rate) … ⛔ not claimed, not derived".
- **Fix:** Quote this, name the expansion rate as the observable, and keep the record's caveat that the mapping is not claimed.

**F6. The PS row says the lifetime yardstick was "explicitly deferred", which the records do not say.**
- **Evidence:** PC:21 and PC:37 assign owners only. `grep -n -i "lifetime\|yardstick"` over PC, d, SR and S9 returns nothing.
- **Fix:** Label it a gap the inventory inferred.

**F7. Review material is cited as if it were a source.**
- **Where:**
  - The SE supernova row: "Claude review quotes `erg cm⁻³ s⁻¹`; unit discrepancy retained".
  - Line 37 (Claude finding 7 versus Grok finding 3), and the R7, R9 and R10 entries at line 143.
  - Leftover repair notes in the tables ("Gain coefficient removed", "Diameter tolerance removed", "Marcuse examples moved").
- **Evidence:** The preprint prints "3 × 10³³ erg−1 cm−3 s−1" but defines Q as "The energy loss rate per unit volume" (fg.txt :328 and :413). By the source's own definition, `erg⁻¹` is a typo.
- **Fix:** Quote the printed unit alongside the source's definition. Remove the review citations and the repair notes. The directive's disagreement rule covers records, not reviews.

**F8. A disagreement between records about status is shown from one side only.**
- **Where:** E4, which quotes "unconditional on a uniform background (B:191)".
- **Evidence:** PC:7 and d:7 mark the uniform decoupling "**CONDITIONAL**".
- **Fix:** Quote both statuses.

**Minor items**
- **Peyton:** the inventory says experiments used titanium. The paper says "Aluminium was chosen for the comparison tests"; titanium was the simulation material.
- **More et al.:** the abstract also publishes a per-length bound, "nσ < 2×10⁻⁴ h Mpc⁻¹". That is the form the matrix's per-length entries need, and the inventory omits it.
- **Marcuse 1969, Eqs. (75)–(76):** these assume "C₀ = 1, C₁ = 0 at z = 0" (OCR "Assuming C = \,C\ = at z ="). That condition is missing from the inventory's assumptions.
- **SR line numbers:** the inventory cites SR:207–212 for the naming distinction, which is at SR:201–203; 207–212 holds the tension.
- **S11 token:** the token "`INTERFACE_EQ={}`" is not verbatim; S11:280 has `INTERFACE_EQUATIONS_SUPPLIED={}`.
- **Corning SMF-28 ULL:** the page returned 403 to me, the row is not marked UNVERIFIED, and it is a product specification rather than an observation.
- **Row 13:** the ABSENT mapping is supported more directly by B:74–76 and B:191 (uniform coupling ≡ 0) than by the S9 shear postulate.

---

## Checked and found sound

- **Directive coverage:**
  - All 12 directive-listed mechanisms are present, plus one more. Each row has the five required columns.
  - The counts are correct: ABSENT 1, PRESENT ANALOG 9, NOT ADDRESSED 3; routes 5 / 3 / 5, by an awk count over the tables.
  - Seeded routes: b:54–57, c1:35–44 and c2:63–70 are present.
  - S11:74–76, d:64 and R-S8-02 are filed under the third status.
  - Folding the thickness channel into E1 is explicit and defensible, since b:50–57 is the gradient-driven channel.
- **Record fidelity:** I spot-checked about 60 citations across S9, S11, A, B, U, SC, b, c1, c2, d, PC and SR. Quotes and statuses match, including:
  - c1:100–110 (ESTABLISHED versus UNDECIDED);
  - c1:146–155 (`NOT_ESTABLISHED_AT_FIRST_SHAPE_ORDER`);
  - c1:188–190 (`K_a` is Hermitian);
  - A:79 (exact);
  - B:23–27, 46–48, 61–68, 80–94, 150–154, 196–202;
  - U:156–164;
  - d:29;
  - PC:3, 7, 9, 27 and 31.

  The c1 "ESTABLISHED" kernel survives the later inertia repair: `_measurements/S11c_inertia_repair_report.md:7` says "all 44 c1 values are identical".
- **Prior art:**
  - Crossref matches every DOI, title, volume and page range, including both Marcuse titles.
  - Matched to source text: Kubrusly Eqs. (1), (8), (15), (16); Peyton Eq. (2) with C = 1.38; Gu–Fuller Eqs. (13), (15), (16) (OCR); Marcuse 1969 (abstract and Eq. 72); and Marcuse 1970 ("2:1 … below 1 percent").
  - No prior-art result is used as a v3 bound or premise.
- **Yardstick numbers verified against primary text:**
  - Tamura: 0.1419 dB/km at 1560 nm.
  - Berceau: τ = 1.16 ms, L_c = 2.27 m, F = 481 000, Δn ≤ 5.0×10⁻²⁰ T⁻² per 4 ms.
  - Lu et al., Optics Letters 34:1708: (99.577 ± 0.002)%.
  - More et al.: Δτ < 0.13 at 95%.
  - Fixsen: 50 ppm, |y| < 15×10⁻⁶, |μ| < 9×10⁻⁵.
  - DES: b = 1.003 ± 0.005 ± 0.010.
  - Lubin–Sandage: n = 2.59 / 3.37, tired light excluded at >10σ.
  - Pratt: Q up to 10⁸, 7 μHz.
  - Friedland–Giannotti: "about twice the standard neutrino luminosity"; the red-giant, supernova and horizontal-branch bounds on k; k is the warped curvature and n the number of compact dimensions.
- **No verdicts:** I found no outcome, suppression requirement or new physics stated.
