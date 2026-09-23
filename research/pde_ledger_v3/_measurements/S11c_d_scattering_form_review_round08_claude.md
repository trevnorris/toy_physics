**AMENDMENT: NEEDS REVISION**
**INVENTORY: CLEAR** (scoped; it needs a mechanical sync if the amendment's A11/A12 text changes)

These verdicts cover document fidelity and whether the physics specification is adequate. I did not rerun, recompute or independently validate any saved result, and a CLEAR here does not clear any later build, export or physics.

---

## Method and source account

I read the governing sources first: v10 in full, `S11c_decisions.md`, nonlinearPoleV2, exploratoryAcceptanceV1, the full S11b step record, the S11b spec §§1–3, energy accounting and B9, the full c1 spec and record, the c1 and d source excerpts, the c1 export-key census, the retained contract excerpt and every saved report. Only then did I read the two artifacts. That order is a method, not a blindness control.

What the sources establish, briefly:
- **v10 §3a** defines `C_{T→H}` from full-end channels.
- **v10 §3b(i)** describes bulk escape as the §3a amplitude "projected on radiating/continuum channels and evaluated with each converted channel's own current" (`S11c_d_SHARED_PHYSICS.md:485-486`).
- **c1** supplies the operator Z, a closed face response and a first-shape power caveat. It exports no exterior solution map: the census lists only DtN/impedance, resolvent and face-input keys (`S11c_d_scattering_form_c1_export_key_index.json`).
- **The saved d current** is slab plus an infinite-depth integral times the bulk-normal current density. It is admitted only for pairs with Im q > 0 on both legs (`S11c_d_continuum_boundary.py`, excerpt line 40, original ≈240: `f.require(left['q'].imag>0 and right['q'].imag>0,'actual convergent bulk-current pair')`).
- **At ω=1:** "q_depth^2 = -k_normal^2 - 1/25" (`S11c_d_remaining_case_flux_report.md:11-12`). The bulk branch frequencies are ±√5 (`S11c_d_frequency_source_report.md:20-21`).

---

## Answers to the five questions

### Q1. Pole deferral, local needs, end-channel versus bulk typing, historical statements, closure attribution

**Mostly yes.**

The pole deferral is carried consistently through every layer:
- **Scope:** the table rows at `AMENDMENT:69-72`.
- **Emission and export:** `:668-669` says "Do not export zero, an empty set, an empty result-bearing container…"
- **Comparator:** `:756-758` says pole families "are not residual zeros, evidence of equality or a fourth T7 truth value".
- **Downstream:** `:716-720` requires "explicit unsupported/deferred capability failure".

Local scattering needs are kept (`:597-608`). The ∂_ω versus ∂_{k_n} pairing distinction (`:609-614`) is physically correct: a root solved in k_n is simple when the ∂_{k_n}𝓛 pairing is nonsingular.

The historical search statement is accurate. The contour report says "no candidate is resolved by this bounded numerical search on |omega-(1-0.01i)|=0.02 … not a certified empty spectrum" (`S11c_d_frequency_contour_refine_report.md:23-27`), and the amendment reproduces it at `:618-626`.

Survival is scoped to "continuum-channel transverse flux survival, not a complete N13 confinement answer" (`:701-703`). Its order status is incomplete (Finding 3).

Closure/interface attribution is handled honestly:
- Signed exchanges and balance checks stay current work (`:328-334`, `:664`, `:685-696`).
- A distinct absorption observable is neither promised nor zeroed. It goes to a named future decision with an escape clause: "If source inspection identifies an existing mandatory operand or check, it remains current work" (`:691-693`).
- No source I read requires d to produce a separately normalized closure-loss observable. N13 names bound capture and bulk radiation (`S11c_decisions.md:130-133`).

**Gap:** the end-channel/bulk split is typed carefully, but the amendment misses that on any radiating domain the end channel space itself contains the bulk continuum (Finding 1).

### Q2. A9, the bulk-flux FORM and the c1 linkage

**Largely yes.**

A9 requires the full reduced operator, both-end channel spaces, derived currents, independent grades including the mixed grade, and actual baselines and interference (`:101-122`). It rejects both failure modes:
- the simple-Born shortcut: "does not license a simpler Born matrix element without §2's computed reduction premises" (`:71`);
- placeholders: "A frequency-one coefficient table, opaque response symbol, fingerprint, unevaluated missing construction or hand-typed expected scaling is not the general computed FORM" (`:125-128`).

It also correctly allows a transparent operator/integral representation instead of an elementary closed form. That matches retained contract item 6.

The bulk FORM is correctly labelled as a new construction, consume-set and T7 duty, not an inherited object: "not an existing `IMPORT_KEYS` row or a map already exported by the slab-only operator" (`:164-166`, `:743-748`). This agrees with the census.

Several c1 details are source-correct:
- **The historical c1 far-field construction** uses equal legs. The excerpt at original lines 1380-1466 substitutes `{q_out_k: qprop, q_out_kp: qprop}` with an impermeable drive. Its three-term integrand omits the scattered–scattered term (original ≈807-812). The amendment describes both points accurately (`:167-170`, `:241-250`).
- **Depth handling:** "Only the in-plane content undergoes the §1c reduction; the exterior depth coordinate remains" (`:161`) is correct. So are disconnected half-spaces (c1 `:99-101`), the graph-versus-outward and Z-versus-N corrections (c1 record `:180-184`), and locate-first elimination of ζ_c (`:174-179`). The last is needed because c1 §1a forbids setting ζ_c=0, and the c2 closed operator is over {u,θ,e_W}.
- **Power orientations** match the record:
  - acoustic check: `P_into_bulk − P_outgoing − P_lateral` (`:214`);
  - traction check: `P_face + P_infinity`, with subtraction operand "minus outgoing power" (`:224-227`), matching c1 record `:174-179`;
  - `P_face` is the bulk-on-slab traction work, i.e. −P_into_bulk, so the two checks are consistent;
  - "Do not use acoustic pressure-work alone as that sign test" correctly reflects c1 spec `:322-324`.
- **The baseline/interference test is right:** "a computed baseline/interference disposition can establish a leading induced quadratic from first-order amplitudes without a second-order amplitude" (`:266-269`). Whether that disposition holds is left to computation rather than pre-classified (`:352-353`), which is correct.
- **Face-side versus far-field:** at the leading bulk grade, the face-side comparison needs the unavailable second-shape Hermitian part on the evanescent nullspace (c1 spec `:332-340`). The far-field |a₁|² does not. The amendment's "mark that comparison unavailable … A supported far-field claim still needs an applicable independent construction/comparison … through the retained blind-engine/T7 process" (`:232-239`) is the right structure.
- **Validity domains are kept separate.** Both-leg first-shape validity comes from c1 record `:146-155`. Omitted-flow validity is at `:355-366`. "Bare `N^-1`/impedance singularity does not prove the complete permeable closure resolvent singular" (`:340-342`) is faithful to c1 record `:150-153`.
- **c1 debt accounting** covers `:99-124` and the independence scoping at `:157-163` (`AMENDMENT:181-203`). **Gap:** it omits the density-as-multiplication-operator carry-forward (Finding 2).
- **Precision is declared per observable** before comparison (`:437-441`).

### Q3. Open-thickness coverage and validity reconciliation

**No existing clause requires a second numerical real frequency.**
- v10 §3a works "at fixed `(ω,k_∥)`" (`:423-424`).
- Acceptance item 1 names only "the approved development example" (`EXPLORATORY_ACCEPTANCE.md:27`).
- The only two-frequency object in the packet is the current-identity diagnostic: "Fresh endpoint sources, independent-frequency pairing and full current/adjoint normalization remain pending" (`S11c_thickness_coordinate_repair_report.md:291-293`). That is a pairing diagnostic, not a scattering frequency.

So the amendment is right to label §3.3 "an explicit proposed acceptance clarification, not a claim that v10 already specified an additional numerical frequency" (`:497-498`), and its stopping and structural-absence rules (`:549-567`) are bounded and demand no global theorem.

**The N11a reconciliation is source-correct.** The full S11b passage reads: "first order **fails** where `|q v₀ / ω| ≳ 1`; in the `k c_s0 ≫ ω` regime `|q| ≈ k`, so this needs `(|v₀|/c_s0)(k c_s0/|ω|) ≳ 1` — large `k c_s0/|ω|` is necessary, not sufficient" (`S11b_interface_coupling_law.md:158-161`). "Necessary" therefore refers to failure, not validity. The amendment's reading at `:368-376`, together with the explicit c1 §2b pair and the subsonic condition (`:357-360`, c1 `:216-219`), is faithful and supplies no new physical classification.

**Two gaps.** Non-vacuity does not exclude a structurally zero amplitude (Finding 4). The end-channel domain does not account for the bulk continuum (Finding 1).

### Q4. Practical targets, comparisons, current meanings and stop conditions

**Adequate in structure.** The amplitude 1e-4 and current 1e-6 goals are confined to the channel frame, each observable gets a pre-declared precision, near-zero results are handled honestly, no unitarity is imposed (`:461`), and there are bounded stops. One retained claim lacks a stated meaning: survival's order status (Finding 3).

### Q5. Representation of saved results

**Accurate.** The four open transverse and zero open thickness directions, the closed bulk depth at ω=1, the zero-rank selector versus retained closed amplitudes, the ≈1 ratio, and the V/P and pending review states all match the saved reports (see the Inventory section). Nothing is silently marked complete.

---

## Findings

### Blocking

**1. On radiating inputs the end channel space contains the bulk continuum. The amendment only guards against a failure of end-normalization convergence.**

*Source:*
- v10 `:485-486` makes bulk escape a projection onto "radiating/continuum channels" of the §3a S-matrix.
- c1 `:115-117` puts the q_out branch points at "the sound cone `ω = ±c_s0|k|`", with k the in-plane wavevector including k_n.
- The saved machinery is discrete and decaying-only: "36 candidate dispositions and 14 selected end directions" (flux report `:10`), "Finite modal boundaries" (domain report `:45`), and the Im q>0 requirement quoted above.
- At ω=1, q_depth² = −k_n² − 1/25, and the branch frequencies are ±√5. Radiating bulk support, which A12 needs, therefore exists only where the end symbol has **real** branch points in k_n. There the field along y_n includes a slowly (algebraically) decaying branch-cut radiation contribution that discrete-mode end closures do not represent.

*Artifact:*
- It addresses only "convergence of a saved infinite-depth end normalization is not established on a radiating bulk domain" (`:322-324`) and "a leaky complex-wavenumber mode cannot be silently treated as one of those channels" (`:537-538`).
- A9 lists only "open incident/output channels and evanescent matching contributions" (`:103-104`).

*Why it matters:* an A11 witness or A12 check run with the saved discrete-mode finite-boundary machinery on such an input would have an incomplete end closure and current accounting. The common-control-volume rule (`:315-326`) assumes end faces carry discrete channels only.

*Minimal correction:* in §2 ("Two continuum observables") and §3.3, state that:
- where ω/c_s0 > |k_∥|, the end asymptotics include the bulk radiation continuum;
- the saved end census, finite modal boundary and infinite-depth current are established only on decaying-bulk inputs;
- any witness or check there must either use an end closure and normalization that include, or bound at the declared precision, that continuum, or be restricted to inputs without real branch points, with the choice recorded.

**2. The density-as-multiplication-operator carry-forward is missing for direct reuse of the c1 closure in the bulk map.**

*Source:*
- c1 record `:185-187`: "Name the live `1/ρ_br,bg` as a multiplication operator … so the O(εη) channel … cannot be emitted as a bare constant."
- c1 record `:135-140`: "The channel is O(εη) … `d(μ_s)/dη|₀ = −μ_θ·w₁/ρ_br`."
- c1 spec `:150-154`: μ_s = μ_θ/ρ_br,bg⁰ enters J_s, and J_s enters the bulk drive `n̂_s·v_bulk,s = V_s + J_s/ρ_m`.

*Artifact:*
- It lists only "seal-5 density representation" as a cross-engine premise (`:181-183`).
- It allows "exact-input reuse" of per-engine-SOUND operands (`:193-195`).
- The carry-forward list at `:159-161` names only the graph/outward and Z/N corrections.

*Why it matters:* reusing the PY face response with a bare constant would drop a first-order term from the A12 face drive, whatever the cross-engine status.

*Minimal correction:* one sentence. Any d use of c1 closure or face-response operands must bind 1/ρ_br,bg as the live background multiplication operator, or consume c2's re-bound closure, and never the bare constant.

**3. Survival's weak-order meaning is unscoped.**

*Source:*
- v10 `:527-536` defines P_T,surv.
- The baseline/interference caveat at v10 `:589-592` is written for J_H.
- The saved report is explicit: "The retained rectangle does not provide omitted pure second-order terms that can interfere with a nonzero baseline; its quadratic transverse-current terms are not a physical tiny gain/loss claim" (`S11c_d_continuum_currents_report.md:32-34`).

*Artifact:*
- `:698-708` fixes survival's export meaning but not its order status.
- The general baseline sentence (`:263-266`) sits inside the bulk-coefficient section.

*Why it matters:* the transmitted transverse baseline is O(1). The leading departure of P_T,surv from unity therefore involves interference with omitted pure-second-order amplitudes. Survival is the one retained N13 object.

*Minimal correction:* in the survival consumer contract, state that its deficit and weak-order coefficients are retained-model data under the same three-status schema (`:257-261`), not a physical loss estimate.

**4. The non-vacuity rules do not exclude a structurally zero amplitude.**

*Artifact:*
- A11 requires an incoming transverse channel and a nonempty thickness channel space (`:500-508`), yet concedes "The presence of an open channel does not itself imply a nonzero conversion amplitude" (`:573-574`).
- A12 requires only "nonempty radiating support and an applicable nonzero incident current" (`:383-385`).

*Source:* v10 §5 says "only a **form** change tests physics" (`:708-709`), and v10 `:321-322` allows form-factor zeros.

*Why it matters:* a witness where the conversion or flux vanishes by symmetry or at a form-factor node compares zero with zero. It passes the stated gate without testing the coefficient, current normalization sign or measure.

*Minimal correction:* require at least one incident/outgoing pair, or bulk support region, where the computed quantity is resolved nonzero at the declared precision and is not forced to zero by an identified selection rule. A resolved zero is recorded as such and does not discharge the coverage. If no such point is admissible, apply the existing stop rule.

### Non-blocking

**5. Mixed-grade completeness in the bulk.** The amendment keeps the mixed amplitude grade (`:255-257`) while listing η·σ_W among the unsupported evanescent completions (`:271-273`, following c1 record `:191-193`). The build-plan dependency map should state explicitly whether c1's first-shape DtN/closure contains the height×slope terms that feed the ησ grade of c2's closed operator. The existing requirement to map dependencies "at each independent … grade" covers this generically.

**6. Process note (not a physics finding).** This is the seventh same-author fold, openly recorded as a departure from G4 (`:22-30`). It does not affect the physics assessment.

---

## Inventory

**Source fidelity: checked against the saved reports, no errors found.**
- **A3:** 645 unknowns, regulator 0.1, 687 s.
- **A4:** 8.82e-3 versus 3.22e-5, q_depth, and the survival/total-ratio coincidence.
- **A5:** 2.22e-13 and 2.97e-13.
- **A6:** coupling triplets correctly marked not established.
- **A7:** 2.23e-5, 7.73e-7 and 8.91e-7; 8.53e-3 versus 3.36e-5.
- **A8:** 8.243e-7 and 5.482e-7; ≈71 min.
- **A10:** 170 s, baseline case only.
- **B1–B4:**
  - 32 owners = 9 LAB + 23 material; the 23 material owners are 15 MAT/RHO4 1D plus 8 2D (rows 60–63 and 50–53), per the material-row summary;
  - condition 5655, 4.67e-14, 12.88 s, 31.29 s, 639.84 s.
- **A1:** 53 min and 10 min, and the 20,736-sample asymmetry.

**Other checks:**
- Done and remaining items are distinguished, and review state is honest: V/P, and the Lean NP1–NP4 clearance is CLEAR only within its scoped packet.
- Work is assigned correctly across retained scattering (A1–A13), deferred poles (B) and the A14 boundary.
- The owner gate for each deferral (N10) is kept.
- Costs are labelled as judgments.

**Optional improvement.** A1/A2 call post-repair normalization coverage "unestablished". The packet contains leads worth citing:
- the pressure-trace report says "accepted LEFT/RIGHT normalization … unchanged" and "REFERENCE normalization has resumed" (`S11c_wolfram_pressure_trace_repair_report.md:75-81`);
- the continuum consumer requires an "accepted complete end normalization" `…_thickness_repair_checkpoint.json` (`S11c_d_continuum_boundary.py`, excerpt lines 7-8).

This is locate-first, so the current wording is not wrong. If Findings 1 and 4 are adopted, A11 and A12 need matching one-line additions.

---

## Coverage limitations

- **Partial sources:**
  - The c1 and d scripts, `CLAUDE.md` and the builder report are partial excerpts.
  - Only S11b spec §§1–3, energy accounting and B9 were read.
  - c2 spec, step record and exports, the build directive, `DEFERRED_HEAVY_RUNS.md`, the round 1–7 reports and all checkpoint JSONs other than the two supplied were not in the packet.
  - I could not check whether c2 eliminated ζ_c, re-bound the density, or what its closed operator retains at ησ.
- **Unverifiable references:** commit hashes (`399a8516`, `7c98b8ee`, `0d77af53`, and others) and "later dispositions" cannot be checked from the packet.
- **Inferences:** the numeric link between the ±√5 branch frequencies and the radiating support in Finding 1 is my inference from the two saved formulas, not a recomputation.
- **What this review is not:** no numerical, CAS or script work was done, and no build or control was treated as executed.