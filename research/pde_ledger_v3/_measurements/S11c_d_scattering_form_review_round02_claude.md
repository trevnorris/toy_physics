# Independent review: S11c-d Option B scattering/FORM scope amendment (draft2) and scope inventory

**Reviewer:** fresh non-author Claude agent, DOCUMENT branch. I read the governing sources before the two artifacts. I make no blindness claim. This review ran no computation, CAS, script or rerun. Line numbers are packet lines. For the excerpted files, they are packet-excerpt lines, and the original line ranges are given in each excerpt header.

## Verdicts

- **AMENDMENT: NEEDS REVISION.** It has two blocking findings (A-1, A-2) and four optional improvements.
- **INVENTORY: NEEDS REVISION.** It has two blocking findings (I-1, I-2) and one optional improvement.

These verdicts cover document fidelity and whether the physics specification is adequate. I did not independently rerun or validate the reported computations.

---

## What the sources establish (my account, formed before reading the artifacts)

- **v10's scope.** v10 requires the complete two-ended S-matrix at fixed `(ω,k_∥)`. It is built from the §1c-reduced operator and kernel, with channels taken from the **full** end pencils (`S11c_d_SHARED_PHYSICS.md:372-377`, `:423-436`). It also requires:
  - currents derived from S11b, "possibly non-Hermitian" (`:440-457`);
  - the conversion `C_{T→H}=J_H,out/J_T,in` (`:459-466`);
  - the transverse survival ratio (`:527-537`);
  - independent `(ε,η,σ_W)` grades, with baseline, interference and quadratic slots, and a leading-order label only after the computed disposition (`:547-617`);
  - weak Taylor coefficients and the unsolved strong-edge handoff (`:633-676`);
  - the §5 controls.
- **The two photon-kill channels.** v10 names two distinct channels (`:481-486`). The first is "continuum conversion — transverse → the thickness **continuum / bulk escape** (the §3a amplitude projected on radiating/continuum channels…)". The second is the bound pole, whose unrestricted residue/projector prescription (`:495-525`) is superseded by nonlinearPoleV2 §§1–7.
- **N13.** N13 treats bulk radiation as a kill mechanism in its own right: "kills the photon exactly as bulk radiation does — the two are distinct emitted objects" (`S11c_decisions.md:131-133`).
- **Bulk radiation and damping.** S11b fixes that bulk radiation exists even with impermeable faces (`S11b_SHARED_PHYSICS.md:243-245`). The face closure is `Λ_I(ω)=Λ_I⁰/(1−iωτ_I)` with "DO NOT set any τ_I = 0" (`:159-160`, `:193-195`). So the thickness/face sector is generically non-conservative at real ω.
- **No second real frequency is mandated.** No existing clause requires a second real frequency. v10 §3a works "at fixed `(ω,k_∥)`" (`:423-424`). The retained contract asks only for "real continuum frequency and tangential momentum" (`S11c_d_sympy_builder_report.md:19-20`). The acceptance addendum asks for "the approved development example" (`S11c_d_EXPLORATORY_ACCEPTANCE.md:27-28`).
- **The FORM obligation.** v10 does require a *FORM* (`:815-818`) with the `s` and `k_aL_W` dependence kept live (`:314-316`). A single closed-channel point cannot exercise the conversion branch of that FORM.
- **Saved evidence at ω=1.** In all four cases there are 4 open transverse directions and 0 open thickness directions. Bulk depth is closed: `q_depth²=−k_n²−1/25` (`S11c_d_remaining_case_flux_report.md:9-15`; `S11c_d_continuum_currents_report.md:4-10`). The positive regulator is 0.1 (`S11c_d_remaining_case_response_report.md:54-56`). The baseline contour on `|ω−(1−0.01i)|=0.02` resolved no candidate and is not certified empty (`S11c_d_frequency_contour_refine_report.md:23-28`).

---

## Answers to the five questions

**1. Pole deferral through output, export and downstream: mostly coherent.**
- Scope, §2, §3b, §§4/7/8, nonlinearPoleV2, acceptance, the retained contract and N2 each get an explicit row (`AMENDMENT:37-47`).
- Export forbids zero or empty placeholders: "Do not export zero, an empty set, an empty result-bearing container…" (`:321`). Status metadata `DEFERRED_BY_SCOPE` is kept separate from values (`:322`).
- T7 keeps three values: "not residual zeros, evidence of equality or a fourth T7 truth value" (`:350-352`). Consumers receive an explicit unsupported-capability failure (`:340-341`).
- Local scattering needs stay in place: the resolvent, closed modes, radiation/current and branch checks are not deferred (`:259-262`), and neither are end-mode pairing safeguards (`:263-270`).
- Historical pole statements are limited correctly (see Q5).
- **Gap:** the deferral is framed as the *only* reason survival falls short of N13. The continuum half of N13, direct bulk-depth escape, is neither retained as an object nor deferred (**A-1**).

**2. A9 is adequately specified.**
- It keeps the full reduced local/nonlocal operator and vertex, the computed reference/left/right baselines, both-end channel spaces and derived signed currents (`:67-71`).
- It keeps independent grades and the mixed grade, with the homotopy as the only link (`:72-78`).
- It keeps the baseline, zero-jet, first-jet and induced amplitudes, the current form, the denominator, and the rule that the η² label must be computed rather than supplied (`:79-86`).
- It rejects "a frequency-one coefficient table, opaque response symbol, fingerprint, unevaluated missing construction or hand-typed expected scaling", but permits transparent operator/integral representations (`:88-95`).
- It keeps source/unit/domain joins (`:96-99`) and the distributional zero-jet step (`:101-104`).
- It explicitly bars a simple Born replacement without the §2 premises (`:41`). This matches v10 `:390-393` and retained contract items 4 and 6.
- No blocker. One optional point about like-for-like witness comparison is **A-5**.

**3. The open-thickness-channel treatment is honestly labelled but misses a premise.**
- The amendment correctly states that no existing clause mandates a second frequency (`:156-158`, `:181-182`, `:232-234`). I agree, and I quote the governing clauses under "What the sources establish" above.
- It correctly separates end channels of `𝓛_±^full` from bulk-depth radiation (`:171-175`, `:239-245`). It excludes leaky complex-wavenumber modes (`:209-211`) and keeps the c1 grazing limitation accurate (`:211-214`, matching `S11c_c1_SHARED_PHYSICS.md:56-61`).
- It does not require a global theorem (`:147-152`, `:216-218`).
- **Missing premise:** the check is only non-vacuous if a real-`k_n` thickness-like end channel can exist at all in the approved model. The inherited dissipative closure and bulk leakage can make that set empty. The regulator can also blur the open/closed decision. The amendment does not require this to be diagnosed first, and it gives no honest outcome for the case where it holds (**A-2**, which I label a **new clarification** of the amendment's own new requirement).
- **The bulk/end distinction has a claim-side gap.** C_{T→H}'s export meaning does not say that it excludes bulk-depth escape (**A-1**).

**4. Practical targets and stop rules are adequate for the limited claims, with one gap.**
- These are adequate: the 1% default for resolved nonzero observables; absolute goals confined to the channel frame (`:121-125`); selected independent comparisons; the rule that near-zero effects stay unresolved (`:140-143`); no assumed unitarity (`:143-145`); the stop rule of "do not scan … indefinitely" followed by a scope decision (`:221-230`); and "pole deferral alone cannot make the … handoff complete" (`:394-395`).
- **Gap:** the retained §5 controls produced only sub-resolution channel changes at ω=1. They therefore do not yet exercise the conversion FORM, and "relevant controls" at a witness is left undefined (**A-3**, optional).

**5. Saved results are represented accurately.**
- These amendment statements match the reports:
  - ω=1 channel counts and closed bulk depth (`:163-170`);
  - the limited status of the algebraic candidates (`:175-177`);
  - the historical search on the stated circle, "closed without a resolved candidate", not certified empty, finite positive-regulator pencil, no holomorphy proof (`:274-282`; compare `S11c_d_frequency_contour_refine_report.md:23-28` and `S11c_d_frequency_source_report.md:12-13`);
  - the Lean NP1–NP4 clearance scope (`:295-299`; compare `POLE_FIDELITY_REVIEW.md:7-9`).
- I found no obligation silently declared complete. The amendment states that the export is unfinished (`:312-313`).
- **Omitted from both lists:** bulk escape (A-1).
- **Contradiction:** in the inventory, B2/B4/B5 are not marked deferred (I-1).

---

## Amendment findings

### A-1 (BLOCKING): The bulk-escape part of v10 §3b(i)/N13 is neither retained nor deferred, and C_{T→H}'s export meaning is not restricted

**Sources:**
- v10 `:485-486`: "(i) continuum conversion — transverse → the thickness continuum / bulk escape".
- N13 `S11c_decisions.md:131-133`: "kills the photon exactly as bulk radiation does — the two are distinct emitted objects".
- S11b `:243-245`: "a bulk radiation loss (energy carried to infinity by §2's bulk, present even with impermeable faces)".
- v10 `:621-622`: "bare bulk Poynting are not substitutes for the channel-resolved currents".

**Artifact:**
- `:42`: "Retain continuum conversion and the transverse continuum-channel flux ratio … It is not the complete N13 two-mechanism confinement object."
- `:239-245`: "No standalone far-field bulk-escape FORM is introduced or claimed by this amendment … It cannot be substituted for the retained `J_H` conversion FORM."
- `:324-331`: the survival meaning is carried "in the object's export definition and consumer contract".

**Physics:** a localized interface scatters into a continuum of `k_n`. Where `ω²/c_s0² > k_∥²+k_n²`, part of that field radiates directly into bulk depth at first order in the vertex. v10's retained `C_{T→H}` is built only from asymptotic end channels (`:462-466`), so it cannot contain this loss.

At radiating parameters, `1 − P_T,surv − C_{T→H}` therefore mixes bulk escape, closure dissipation and truncation. The amendment calls the gap "two-mechanism", which attributes the whole shortfall to the bound-pole deferral. It also frames bulk escape as "not introduced", as if v10 never named it. No consumer-contract restriction is placed on `C_{T→H}`, so S11c-e could read it as the total continuum kill.

This is a scope change, possibly a legitimate resolution of an ambiguity in v10, that is not recorded as one.

**Minimal correction:**
1. Add a precedence row stating the change: v10 §3b(i)'s bulk-escape reading is not delivered as a separate object. Either defer it with a named owner or record it as unattributed current-balance content.
2. State in the `S11CD_CONTINUUM_CONVERSION` / `C_{T→H}` export definition and consumer contract that it means "end-channel thickness conversion only; excludes bulk-depth escape and closure dissipation".
3. Replace "two-mechanism" with wording that names both the deferred bound channel and the non-delivered bulk-escape object.

**Why the existing wording doesn't cover it:** `:239-245` only forbids *substituting* bulk escape for `J_H`. It does not stop a consumer from treating `J_H` as the whole continuum kill channel, and the only export-level meaning restriction applies to survival.

### A-2 (BLOCKING; a new clarification of the amendment's own new §3.3 requirement, not a v10 mandate): the openness premise for the A11 witness is not diagnosed, and the case where it structurally fails has no disposition

**Sources:**
- v10 `:433-434`: "Emit closed/evanescent modes as matching data, not flux channels."
- v10 `:455-457`: "multi-component, frequency-dependent, possibly non-Hermitian current".
- S11b `:159-160`: `Λ_I(ω) = Λ_I⁰/(1 − iωτ_I)`.
- S11b `:193`: "DO NOT set any τ_I = 0".
- c1 `:34-36`: the sound-cone branch points `ω=±c_s0|k|`.
- Saved state: "all six decaying thickness directions" (`continuum_currents_report:5`); "regulator 0.1" (`response_report:55`).

**Artifact:**
- `:209-211`: "Establish the real-normal-wavenumber, flux-carrying asymptotic channels … a leaky complex-wavenumber mode cannot be silently treated as one of those channels."
- `:205-206`: "Do not alter signs, constitutive parameters or boundary conditions simply to force leakage."
- `:221-230`: "record exactly what failed or remains untested … If neither non-vacuous route is established, A9 remains incomplete."

**Physics:** a thickness-like end root can fail to be a real-`k_n` channel for four different reasons:
- (a) it is below its cutoff and evanescent;
- (b) it lies above the bulk sound cone for its own in-plane wavevector and leaks through the DtN;
- (c) it is damped by the retained face-flux closure;
- (d) the positive regulator broadens it.

Only (a) can be escaped by moving to another admissible frequency. If (b) or (c) holds throughout the admissible real-ω domain of the approved parameter map, the amendment's criterion makes A11 **structurally unsatisfiable**. A11 would then stay "incomplete" indefinitely, and the parameter change that could fix it is forbidden. Meanwhile the physically real conversion into damped breathing excitation, which N13 counts as a photon kill, appears only as an unattributed survival deficit.

The reports do not say which mechanism closes the six thickness directions at ω=1. I could not determine this from the packet (see Coverage).

**Minimal correction:** in §3.3 bullet 2:
1. Require the bounded inspection first to classify the closure mechanism of each thickness-like end root on the full end pencils, as (a)–(d), with the regulator's effect on that classification removed or bounded.
2. State that roots whose `Im k_n` comes from retained physical leakage or dissipation are not open channels, consistent with v10 `:433-434`.
3. If (b) or (c) holds across the admissible domain, record the generic `J_H` branch as *structurally empty for the approved model*. That is a model statement with its mechanism, not an untested numerical gap, and conversion into damped thickness excitation is then reported only through survival and current balance. The scope/parameter decision of `:230` can then be taken on that basis.

**Why the existing wording doesn't cover it:** "record exactly what failed" does not require identifying the mechanism, and it does not distinguish a structural absence from an unfound witness.

### A-3 (optional): the retained §5 controls do not exercise the conversion FORM at ω=1

- **Sources:**
  - `S11c_d_remaining_case_profile_response_report.md:29-31`: "These scattering/current changes lie below the declared absolute reporting resolutions".
  - First-jet sensitivities are at most 8.243e-7 and 5.482e-7 (`…first_jet_response_report.md:23-25`).
  - CLAUDE.md `:42-44`: "a FORM control tests physics".
- **Artifact:** `:106-107`, "All four physical cases and the v10 §5 … controls remain"; `:190-191`, "an actual numerical witness with relevant controls".
- **Issue:** v10 prescribes no residual disposition (`:793-795`), so this is not a fidelity breach. However, the handoff must not imply that the §5 controls validated the profile or first-jet dependence of `J_H`.
- **Suggested fix:** either define "relevant controls" at a witness to include the §5c FORM ablation and one §5a first-jet probe evaluated on `J_H` (label this **new**), or state that those controls are untested on the nonempty conversion branch.

### A-4 (optional): the end-mode pairing safeguard leaves the spectral parameter ambiguous

- **Sources:** v10 `:446-454` normalizes with `⟨l_a,(∂_ω𝓛_e^full)r_b⟩` and requires the `∂_{k_n}`/`∂_ω` identity. The pole contract `:79-90` defines `D = W L'(ω_*) V` in the pencil's own parameter.
- **Artifact:** `:263-265`: "full left/right pairing safeguards for the actual modal pencil and spectral parameter".
- **Issue:** at fixed ω the end roots are roots in `k_n`. Whether a `k_n` root is semisimple (it fails at a thickness cutoff or zero-group-velocity point, exactly where an A11 witness would sit near onset) is governed by the `∂_{k_n}` pairing. The v10 normalization uses the `∂_ω` pairing.
- **Suggested fix:** name both, and require both to be nonsingular at any reported channel.

### A-5 (optional): witness comparisons should be like for like under the baseline disposition

- **Source:** v10 `:589-593`: "If the computed baseline `a₀` participates, omitted pure second-order amplitude/current terms can also enter the parent-theory `λ²` flux".
- **Artifact:** `:189-191`, together with `:77-78` ("A finite-contrast retained solve is not a higher-order continuum prediction").
- **Suggested fix:** state that a numerical witness checks the weak FORM against retained-grade extractions only. Agreement at order λ² with a finite-contrast solve may be claimed only under the computed `a₀`/interference disposition.

### A-6 (optional): record the N2/N10 family-level consequence

- **Sources:**
  - N2's S11c-e row: "confinement interpreted here (N13)" (`S11c_decisions.md:53`).
  - N10: "whether light's confinement is unconditional" is S11c scope, and "the family's card names each deferral's owner" (`:171-176`).
- **Artifact:** `:25-28` ("not … an assignment to S11c-e") and `:341-343`.
- **Suggested fix:** state that the S11c roll-up card and S11c-e's N13 interpretation must carry confinement as open, with the named package recorded as its owner.

### Verified as faithful (no finding)

- The precedence table's treatment of the incorrect v10 projector: "Do not reinstate v10's incorrect unrestricted projector/residue prescriptions" (`:44`), and `:305-306`.
- The resource guard is stricter than the acceptance addendum's four workers (`:45`, `:383-387`).
- Survival keeps v10's reflected-plus-transmitted definition (`:324-327`, compare v10 `:530-532`). Only its claim is narrowed.
- The c2 operand debt is still propagated (`:354-355`).
- Existing artifacts are not relabelled (`:33-35`).

---

## Inventory findings

### I-1 (BLOCKING): the pole-branch rows and the "What to defer" section do not reflect the adopted Option B

- **Evidence:**
  - B1 says "**0 new pole-support production under B**" (`inventory:125`).
  - B2 still carries cost "**M–L** for rows, assembly…" with no Option-B disposition (`:126`).
  - B4 is class "D, scoped to the chosen regions" and says "Other case searches remain … **L–XL**" (`:128`).
  - The "recommended scope" deferral list (`:144-157`) omits B2 and B4.
  - `:161-164` states: "Neither can the remaining casewise bounded pole-search obligation be silently relabeled as optional".
- **Artifact contrast:** the amendment's `:249-255` defers "new profile-frequency searches in any case; the eight unfinished material complex-frequency rows and two corresponding responses".
- **Why it matters:** the inventory's own stated purpose is to assign work to retained scattering versus deferred poles (`:18-20`). As written, the row-level retained/deferred assignment contradicts the adopted amendment and its own B1 row.
- **Minimal correction:** add an explicit Option-B disposition (deferred to the named package, 0 under B, with the v10-only cost kept for reference) to B2, B3 (packaging only), B4 and B5, and add B2 and B4 to the `:144` list. `:161-164` can stay as the unamended-v10 statement if it is labelled as such.

### I-2 (BLOCKING, depends on A-1): no row covers bulk-depth escape or the meaning of C_{T→H}

- The A table covers conversion and survival (A4) and A11. Nowhere does it say that the v10 §3b(i)/N13 bulk-escape object is neither delivered nor deferred.
- A4 (`:112`) says only "This establishes neither generic leakage absence nor complete N13 confinement."
- **Fix:** add a row, or extend A4/A11, to follow whatever disposition A-1 settles, including the export-meaning restriction.

### I-3 (optional): v10 §2/§3a baseline objects are not individually tracked

- v10 requires `K₀,K₋,K₊` to be emitted (`:372`, `:473`), along with the surfaced "reduced-operator-block-vs-reduced-kernel off-diagonal residual" (`:339-341`).
- A2 (`:110`) says only "Computed reduction, end modes, native current normalization…".
- Given A6's honest note that the coupling triplets are "not established by that report", map these objects explicitly rather than folding them into A2.

### Verified as faithful (no finding)

I checked each of the following against its saved report:

| Inventory item | Supporting report |
|---|---|
| A3: 645 unknowns, 4 incident columns, regulator 0.1, 687 s | `response_report:3-5,54-56` |
| A4: 4 open / 0 thickness, ≈278 s | `flux_report:4-10` |
| A5: 2.22e-13 / 2.97e-13 | `coordinate_report:90` |
| A6: triplet caveat | `uniform_report:18-24` |
| A7: sub-resolution changes, bump moment is not a scattering calculation | `profile_report:29-31,45` |
| A8: 8.243e-7 / 5.482e-7 | `first_jet_report:23-24` |
| A10: baseline-only domain set, 170 s | `domain_report:20,40-45` |
| B1: 24 of 32 = 9 LAB + 15 material | `…row_sensitivity_report:26`, `…rows_1d_report:26` |
| B3: condition ≈5655 and row-46 figures | Not re-verified line by line; consistent with the evidence index |
| B4: 16/32-point closure | `contour_refine_report:23-28` |
| Review states for v10 r10, nonlinearPoleV2, acceptance and NP1–NP4 | Evidence JSON `:736-743`, r10 legs, `POLE_FIDELITY_REVIEW.md:3-9`, pole contract `:9-10` |
| Absent `S11c_d_exports.py` | Evidence JSON `:729-731` |
| Paused draft | `paused-for-scope-assessment.json:2-5` |
| r10 post-clearance folds are present in v10 | v10 `:55`, `:345-349`, `:727`, `:831-833` |

Costs are clearly labelled as judgment (`:97-103`).

---

## Coverage limitations

- **Excerpts only.** S11b, S11c-a, S11c-c1, CLAUDE.md and the retained contract were supplied as excerpts. I could not check the retained-suffix SHA `f0102451…` against the excerpt; the whole-file SHA `4b083fe7…` matches the evidence index.
- **No parameter map.** The packet has no parameter map (`Λ_I⁰`, `τ_I`, `c_s0`, `k_∥`, contrasts, profiles) and no end-root data. A-2 is therefore stated conditionally. I could not determine the closure mechanism of the ω=1 thickness roots or where the regulator acts.
- **Missing documents.** The build directive, the round-1 review reports, the checkpoint JSONs (beyond evidence-index fields), and the c2 exports/step record were not supplied.
- **Not checked.** `POLE_HANDOFF.md` was only spot-checked. I verified no hashes or annex objects.
- **Status of reported numbers.** All reported numbers are taken from the saved reports, not recomputed.
