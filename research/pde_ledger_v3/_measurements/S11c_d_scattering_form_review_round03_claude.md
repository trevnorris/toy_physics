# Independent document review: S11c-d Option B scattering/FORM amendment (draft 3) and scope inventory

**Reviewer:** fresh non-author Claude agent. **Method:** document branch. I read the sources first, then the two artifacts. This is read-only: I ran no CAS, numerical work, reruns or builds. Where I did simple algebra on reported numbers, it is labelled as inference. **Scope of any CLEAR:** document fidelity and adequacy of the physics specification only. It does not validate the reported computations.

## Verdicts

- **AMENDMENT: NEEDS REVISION.** Two narrow blocking findings (F1, F2). Everything else is consistent and well bounded.
- **INVENTORY: NEEDS REVISION.** It carries the same two gaps (F6). Its source fidelity is otherwise good.

---

## Answers to the five questions

### Q1. Is the pole branch deferred coherently while scattering needs and historical limits are kept? Yes, for the pole branch. There is one gap in the loss-channel typing (F2).

**Deferral chain is complete.** The pole obligation is removed at every layer:
- scope: amendment table rows at lines 42 and 45;
- operator exemption: line 43;
- emit/export/comparator: line 47, §5 table lines 388–389, and the T7 paragraph at lines 431–436;
- nonlinearPoleV2: line 48;
- exploratory acceptance: line 49;
- retained contract: line 50;
- S11c-e boundary: lines 51 and 417–419.

Each is a labelled replacement ("Do not export zero, an empty set…", line 388; "not residual zeros … or a fourth T7 truth value", lines 433–434). It is not a silent deletion. This matches the no-placeholder rules in the sources:
- NONLINEAR_POLE_CONTRACT lines 201–202: "Unresolved constructions have explicit status and evidence, not zero or empty replacement payloads."
- EXPLORATORY_ACCEPTANCE lines 22–23: "No absent integral, mode, pole search, expansion or export may be replaced by … zero or an empty set."

**Local scattering needs are kept.** Lines 319–336 keep:
- the resolvent, closed modes and local branch checks;
- the nonlinearPoleV2 §3 full-pairing safeguards for end modes.

They also correctly separate the ∂_ω𝓛_e normalization pairing from the ∂_{k_n}𝓛_e semisimplicity pairing. This is physically right: end modes are k_n-roots at fixed ω. v10 §3a lines 447–448 fixes ⟨l,(∂_ω𝓛)r⟩=1, and a threshold coalescence shows up as a singular ∂_{k_n} pairing.

**Historical pole statement is accurate.** Amendment lines 340–348 say "circle centered at `1-0.01i` with radius `0.02` … closed without a resolved candidate … not a certified empty spectrum". This matches frequency_contour_refine_report lines 23–27: "no candidate is resolved by this bounded numerical search on |omega-(1-0.01i)|=0.02 … not a certified empty spectrum, a contour-interior holomorphy/exceptional-locus proof". The added caveat about the positive regulator is correct and prudent. Lines 361–365 describe the Lean NP1–NP4 status in line with POLE_FIDELITY_REVIEW lines 5–9.

**Confinement limits are accurate.** Lines 400–410 restrict survival to "continuum-channel transverse flux survival, not a complete N13 confinement answer". Line 408 notes that the near-unity ω=1 ratio "does not establish confinement". This is consistent with:
- flux report lines 9–15 (four open transverse directions, zero open thickness directions);
- response report lines 26–27 (ratios 0.9999995945–1.0000001045).

**Bulk escape vs end channels.** The separation of end-channel `J_H` from bulk-depth escape (lines 110–131, 391–398) is physically correct. v10 §3a lines 463–465 builds `J_H,out` only from end-channel currents `𝓙_n,H,e,out`, summed over e∈{−,+}. Radiation into depth is eliminated through the DtN closure in the reduced y_n problem, so it is not an asymptotic end channel.

**The gap (F2).** By narrowing `S11CD_CONTINUUM_CONVERSION`, the amendment explicitly names a third loss type, closure dissipation (line 394). It then gives that loss no retained or deferred deliverable.

### Q2. Does A9 keep the general computed FORM without a Born shortcut or placeholders? Yes. The bulk-flux obligation is correctly labelled as unconstructed, but it lacks a validity-domain rule (F1).

**Full distorted-wave construction is kept.** §2 items 1–3 (lines 71–90) keep:
- the full reduced local/nonlocal operator;
- computed reference/left/right baselines;
- both-end channel spaces and derived signed currents;
- evanescent matching;
- independent ε/η/σ_W grades including the mixed grade, with η–σ related only by the homotopy;
- baseline, zero-jet, first-jet and induced components;
- the outgoing current form and the incident denominator;
- the distinction between the total fraction, the induced-field quadratic diagnostic and the weak coefficients.

This reproduces v10 §3c lines 556–603 and §3d lines 637–639.

**No Born shortcut.** Line 44 ("does not license a simpler Born matrix element without §2's computed reduction premises") matches v10 §2 lines 390–393.

**No placeholders.** Item 4 (lines 91–99) excludes "a frequency-one coefficient table, opaque response symbol, fingerprint, unevaluated missing construction or hand-typed expected scaling". It also correctly allows a transparent operator/integral representation instead of demanding an elementary closed form. This matches retained contract items 6–7.

**Units and domains.** "Physical source/current units" (line 117) and item 5's operand/source/unit/domain joins are present.

**Bulk flux is not presented as computed.** Lines 117–118 say "No numerical value or closed expression for it is supplied here". Lines 135–136 say it is "not established by the existing real-frequency reports". This is honest. Two further points are correct:
- Line 128 forbids inferring bulk escape as `1−P_T,surv−C_{T→H}`.
- Lines 125–127 limit the S11b pressure-work discriminator to its stated subcase. This matches S11b lines 213–224: real ω, `q_out²>0`, impermeable faces, `Λ_X⁰=0`.

**Relation to the source scope.** v10 §3b(i) (lines 485–486) says: "continuum conversion — transverse → the thickness continuum / bulk escape (the §3a amplitude projected on radiating/continuum channels …)". So v10 did include bulk escape, but only through §3a channel projections. The amendment's separately constructed half-space flux FORM is a legitimate, correctly labelled clarification and addition ("explicit clarification … and an addition", lines 133–134).

**What is missing (F1).** The obligation is not tied to the inherited physical validity regime, and it has no non-vacuity or stop rule.

### Q3. Is the open-thickness-channel coverage non-vacuous and bounded? Yes for end channels. It is under-specified for bulk-depth radiation (F1).

**No existing clause requires a second real frequency.** I found none in the supplied sources:
- v10 §3a (lines 423–431) requires "every incoming open mode" at "fixed `(ω,k_∥)`".
- Retained contract item 3 requires only "Test actual reference/end channel availability … Do not manufacture an incident channel."
- Exploratory item 1 uses "the approved development example".

The thickness-coordinate repair's "retained two-frequency discrepancy" (repair report lines 327–328) is a pairing/current-identity diagnostic with independent frequency legs. It is not a second physical scattering frequency. The inventory says this correctly (A1, "not an automatic second open-thickness frequency mandate"). The amendment correctly labels §3.3 as new: "explicit proposed acceptance clarification, not a claim that v10 already specified an additional numerical frequency" (lines 224–225). The same point is repeated at lines 293–294.

**The end-channel check is well bounded.** §3.3 (lines 227–291):
- defines the checked object as v10's `C_{T→H}=J_H,out/J_T,in` on full-pencil channels;
- excludes zero-rank selectors, bulk-only radiation and current deficits as discharges;
- requires a diagnosis of the root-closure mechanism, including whether the regulator enters the end pencil at all;
- forbids changing signs, constitutive parameters or closure parameters to force leakage (consistent with S11b line 193: "DO NOT set any τ_I = 0");
- requires real-normal-wavenumber flux channels;
- treats local conditioning as practical evidence, not a proof that the domain is pole-free;
- separates "structural absence on a specified domain" from "no witness found", with a scope-decision stop.

This gives non-vacuous coverage without demanding a global theorem.

**Bulk-depth radiation is not equally specified.** Lines 300–305 say only that the bulk-escape FORM "must carry its own radiating-support/domain and current/measure evidence". There is no parallel non-vacuity or stop rule, and no reconciliation with the inherited rest-frame validity condition (F1).

### Q4. Are the tolerances, comparisons, current meanings and stop rules adequate for the retained claims? Mostly yes.

**Tolerances.** §3.1 (lines 164–195) faithfully merges exploratory items 3–5. The 1e-4 amplitude / 1e-6 current goals are applied only in their declared channel frame (lines 165–168). This matches the reports' "original amplitude 1e-4/current 1e-6 reporting resolution" (flux report lines 42–43).

**Near-zero handling and unitarity.** Near-zero effects and denominator floors (lines 184–186) and "do not impose lossless S-matrix unitarity by assumption" (lines 187–188) are correct.

**Current meanings.** These are carried in the export definitions, not only in prose (lines 391–404).

**Stop rules and completion.** Lines 471–477 correctly refuse to let the pole deferral stand in for A9 completion.

**Remaining gaps:**
- There is no bulk-escape stop rule (F1).
- The resonance-proximity status of the exported FORM domain should be an explicit field (F3, optional).

### Q5. Are the saved four-case results and review states represented accurately? Yes. One loss route is missing from both lists.

**Results are accurately represented.** I checked each claim against its report:
- Four-case ω=1 finite and independent-grade responses (645 unknowns, four incidents): response report lines 3–13.
- Current bookkeeping: flux report lines 9–24.
- q_depth² = −k_normal² − 1/25: flux report line 12; continuum report lines 8–9.
- Coordinate agreement ~2.22e-13 / 2.97e-13: coordinate report lines 90–91 and 279.
- Twelve uniform controls: uniform report lines 3–21.
- FORM control changes below resolution: profile report lines 23–31.
- First-jet changes ≤8.243e-7 / 5.482e-7: first-jet report lines 23–24.

None of these is described as reviewed or cleared.

**Review states are accurately represented.** Round-1 and round-2 outcomes, v10 round-10 clearance, nonlinearPoleV2 "not a claim of independent review" (contract lines 9–10), and the Lean CLEAR scope all match.

**Omission.** The closure/interface-absorption loss route is in neither the retained nor the deferred list (F2, F6b).

---

## Findings

### F1 — BLOCKING (amendment): the bulk-depth escape FORM has no validity-regime rule and no non-vacuity or stop rule

**Sources:**
- v10 §6 lines 839–841: "Every result inherits the c1/S11b smallness domain (`|q_out·v_bulk_normal_0/ω|≪1` + boundary-layer/subsonic; large `k c_s0/|ω|` is **necessary, ⛔ not sufficient**)".
- Decisions N11(a) line 183–185: "**every** S11c spectrum/leakage/confinement result is **conditional on the derived smallness domain** … large `k c_s0/|ω|` is necessary".
- c1 §2b (original lines 204–235): "the requested **grazing** behaviour (`q_out→0`) is the **strict `v_bulk_normal_0=0` result**".

**Kinematics.** Every Fourier component that radiates into depth has in-plane |k| < ω/c_s0. So on each such component k c_s0/|ω| < 1, and the radiating window always ends at grazing points where q_out→0.

At the approved inputs, the saved reports place bulk radiation only above a branch frequency:
- flux report line 12: q_depth² = −k_n² − 1/25 at ω=1;
- frequency_source report line 21: "bulk branch frequencies are plus/minus sqrt(5)".

My inference, assuming the standard q² = ω²/c_s0² − k_∥² − k_n² form: these two numbers imply k_∥c_s0/ω ≈ 2.24 at ω=1. So the radiating support opens only when k_∥c_s0/|ω| < 1. That is exactly where the inherited necessary condition, read literally, fails.

**Amendment.** Line 122 says only "Establish the domain where the flux integral and any baseline/excess separation are well defined". That is mathematical definedness, not physical validity. The c1 grazing clause appears only inside the `J_H` bullet, as "Any accompanying bulk-radiation calculation" (lines 263–266). v10 §6's k c_s0/|ω| condition is not cited anywhere in the amendment. I grepped both artifacts for `c_s0`, `N11`, `rest-frame` and `smallness`; only lines 257 and 266 matched.

§3.3 gives `J_H` a structural-absence and scope-decision stop (lines 283–291). Bulk escape has no analogue: lines 137–138 say "If a required object is missing, cost that work and complete it".

**Why this is a physics defect.** As written, A9/A12 could either:
- export a bulk-escape FORM whose entire nonempty support lies outside the model's stated validity domain, while the export presents it as a toy-model leakage prediction; or
- pursue an obligation that has no admissible non-vacuous domain, with no stop rule.

**Minimal correction.** Add to §2's bulk-escape paragraph:
1. The FORM's declared domain must be reconciled explicitly with v10 §6 / N11(a) and c1 §2b.
2. Components or regions with k c_s0/|ω| ≲ 1, and the grazing endpoints, must be labelled as strict `v_bulk_normal_0=0` / outside-inherited-domain results, not as validated leakage.
3. Apply §3.3's structure to bulk escape: a bounded admissibility inspection; if radiating support inside the inherited regime is structurally absent or unestablished, record the mechanism and domain and request a scope disposition. Do not export a zero. A12 remains unresolved.

This adds no new physics mandate. The meaning of "large k c_s0/|ω|" lives in `steps/S11b_interface_coupling_law.md:159-161`, which was not supplied (see coverage limitations). The correction requires that reading to be settled from source, not assumed.

### F2 — BLOCKING (amendment): the interface/closure-absorption loss route has no disposition

**Sources:**
- Decisions N13 lines 131–133: "Conversion into a **bound** breathing/thickness mode kills the photon exactly as bulk radiation does … a grating fails if the photon becomes a breathing excitation even with zero radiation."
- N2 table line 52: "the profile-conditioned **leakage rates** live here".
- S11b §4 lines 159–197: the closure uses `τ_I ≥ 0`, relaxation is never set to zero, and it is a source of real-axis complexity and dissipation.

**Amendment.**
- Lines 391–395: `S11CD_CONTINUUM_CONVERSION` … "exclude bulk-depth escape and closure dissipation as separately attributed losses … They are not the total continuum-loss functional."
- Bulk escape gets a retained root (§5 table, line 384). Bound capture goes to the pole package (line 388). Closure absorption gets neither.
- It appears only as a balance term: lines 186–188, and lines 303–305 ("a deficit cannot separate escape from dissipation … Report unresolved attribution explicitly").
- Lines 289–291 say "Damped thickness excitation is not thereby a separately measured escape or capture rate."

**Physics.** Near the interface, transverse → thickness conversion into closed or damped thickness excitation that the closure then dissipates is a weak-order loss at the same N12 order as `J_H` conversion. At ω=1 it is the only loss route not kinematically closed: end-thickness channels and bulk depth are both closed there. Line 241 of the amendment itself lists "retained closure dissipation" as a possible reason thickness roots are closed.

**Why this is blocking.** v10 already omitted this route. The amendment creates an explicitly named loss type with no owner. It would reach S11c-e as an unowned "unresolved attribution", even though N13 counts it as photon loss.

**Minimal correction.** Add a §5 table row and a precedence-table note that give closure/interface-exchange absorption an explicit disposition. The user can choose either:
- (a) retain it as a separately typed, incident-current-normalized, weak-graded balance observable built from S11b's signed exchange terms; or
- (b) name it as a deferred or out-of-scope item with an owner, and state in the consumer contract that `C_{T→H}` plus bulk escape plus survival do not bound total N13 loss.

Option (a) would be a new requirement. Option (b) costs nothing and resolves the omission.

### F3 — Optional (amendment): record resonance-proximity status in the FORM's exported domain

v10 §2 lines 397–399 says: "The continuum Born domain excludes thresholds, resonance enhancement, modal-gap closures".

Deferring the pole study removes the instrument that would detect a near-real-axis profile pole. The amendment handles this well for the §3.3 pilot (lines 267–272) and for S11c-e generally (lines 417–419). I recommend that §2 item 2 or the §5 FORM row also require a domain field on the exported FORM: "resonance proximity not assessed except at locally checked points". This keeps a bindable ω-domain from being read as resonance-free.

### F4 — Optional (amendment): clarify whether moving ω or k_∥ counts as a change to approved input

Line 260 says "A change to the approved physical input needs its applicable decision gate." Retained contract item 3 includes "real continuum frequency and tangential momentum" in the parameter map. It should be stated whether moving (ω, k_∥) for the A11 pilot counts as such a change.

### F5 — Optional (amendment): place the bulk-field reconstruction map under the reduced-representation rule and T7

The bulk-flux FORM needs a half-space field-reconstruction map. That map is not one of the two §1c reduced rows. v10 §2 lines 345–352 forbids any flux built from unreduced 3-D content.

I recommend naming this map as a reduced construction operand, with dimensions and grades, and listing it explicitly in the T7 retained-coverage sentence (lines 429–431). Inventory C3 already does this.

### F6 — BLOCKING (inventory): the two amendment gaps carry through

- **(a)** A12 (inventory line 121) says nothing about the inherited validity regime and has no stop or unresolved route. A11 has one. Its cost cell says "No extra frequency … mandated". It should add that, at the approved inputs, the saved reports place bulk radiating support only above the ±√5 branch (frequency_source report line 21). So non-vacuous bulk-escape evidence cannot come from ω=1, and its validity status follows F1.
- **(b)** A4 (line 113) repeats that `C_{T→H}` "excludes … closure dissipation", but no row owns that loss. Add a row or note matching whatever F2 disposition is chosen.

### F7 — Optional (inventory): cross-reference the open pairing diagnostic from A4

A1's open post-repair two-frequency endpoint-pairing diagnostic (repair report lines 327–328 and 388–390; no fresh pairing result recorded) is the check underlying v10 §3a's current identity (lines 454–455). A4's `J_T,in` denominators and survival ratios depend on that normalization. A1 currently refers it only to "A2/A11 input inspection". It should also be cross-referenced from A4.

### F8 — Optional (inventory): add the sub-resolution caveat to A8

A8 should state, as A7 already does, that the first-jet channel-level changes (≤8.243e-7 amplitude, 5.482e-7 current) are below the declared resolution. The amplitude change is also below the measured interval spread of 2.14e-6 (domain report line 29). So the one-sided probe's discriminating power at the reported observable is not established. The amendment already covers this at lines 147–149.

---

## Inventory source-fidelity checks that pass

- **B1 counts:** 32 − 24 = 8 missing owners; nine LAB (remainder-rows checkpoint) plus fifteen material 1D (row-2 pilot plus fourteen 1D owners).
- **B2 row numbers:** match material-row-summary.json (MATERIAL/RHO4 rows 60–63; MATERIAL/RHOBR rows 50–53).
- **A6:** correctly separates the response controls from the unestablished `K_uniform,end−K_uniform,reference` triplets (the uniform report contains no coupling triplets).
- **A1:** the Mathematica audit / pressure-trace repair dispositions match audit report lines 12–17 and repair report lines 43–49 and 75–81.
- **A10:** limited to the baseline domain study (domain report lines 26–30 and 40–45).
- **Other claims:** the frequency-source algebraic candidates keep their "not classified physical thresholds" status. Cost bands are labelled as planning judgments (lines 98–104). No row silently claims completion.

## Coverage limitations

- **Missing validity source:** `steps/S11b_interface_coupling_law.md:159-161`, the source of the "large k c_s0/|ω| necessary" condition, is not in the packet. Its exact reading, whether it applies per component or to the incident wave, could not be checked. F1 therefore requires the reconciliation rather than asserting an outcome.
- **√5 inference:** this is algebra on two reported numbers under the standard bulk-dispersion form. It is not a recomputation.
- **Other documents not supplied:** the build directive, program brief, focused completion plan, c1 §3a grazing request, c2 step records, and round-1/2 review reports. The round histories were checked only against the artifacts' own statements.
- **Reported results:** raw arrays and executables are absent. Every reported result was checked for document fidelity only; none was revalidated.
- **Unchecked claims:** A2's "computed reduction" coverage and v10's both-operand 3-D→1-D records cannot be checked from the supplied reports.
