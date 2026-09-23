# Independent review: S11c-d Option B scattering/FORM amendment (draft 7) and scope inventory (draft 7)

**AMENDMENT: NEEDS REVISION.** There are two small blocking corrections (findings 1–2). Both are one- or two-sentence fixes and neither changes the adopted Option B scope.

**INVENTORY: NEEDS REVISION.** Two small required corrections (findings 6–7) and one citation fix (finding 8). The inventory would also inherit the fix from finding 2.

Both verdicts cover document fidelity and whether the physics is specified adequately. I did not rerun, recompute or independently validate any reported calculation. Reading order was: governing sources first (v10, decisions, pole contract, exploratory acceptance, the S11b step record and the relevant part of the S11b spec, the c1 spec §§0–3, the c1 step record, the source excerpts and saved reports), then the two artifacts. This is a method statement, not a blindness claim.

---

## Answers to the five questions

### Q1. Pole deferral, typed obligations, closure attribution

**Pole deferral.** The deferral is threaded consistently through every layer:

- **Scope.** The precedence rows at amendment `:62`, `:65` and `:69` remove the pole deliverable.
- **Operator.** The pole-solve exemption is limited to historical or later work (`:63`).
- **Emission, export and comparator.** Pole families leave required membership, with no "zero, an empty set, an empty result-bearing container" (`:628`). A consumer gets "an explicit unsupported/deferred capability failure" (`:672–673`). In T7, pole families are "not residual zeros, evidence of equality or a fourth T7 truth value" (`:706–708`).
- **Downstream.** "S11c-e may use the scattering/weak FORM … while resonance/local-spectrum and complete confinement claims remain unavailable" (`:673–675`).

**Local scattering needs are kept.** The resolvent, closed modes, radiation/current and local branch checks stay (`:557–560`). The amendment also correctly separates the ∂_ω normalization pairing from the ∂_{k_n} pairing that governs semisimplicity of a root solved in k_n (`:569–575`). This matches v10 `:447–448` and pole contract §3 `:79–103`.

**Historical pole statements are accurate.** Amendment `:578–586` matches the contour report: "no candidate is resolved by this bounded numerical search on |omega-(1-0.01i)|=0.02 … not a certified empty spectrum" (contour report `:23–27`). The Lean NP1–NP4 scope (`:599–603`) matches POLE_FIDELITY_REVIEW `:7–9`: "not the analytic existence theory or a physical S11c pole/scattering calculation".

**Confinement limits are correct.**
- Survival is restricted to "continuum-channel transverse flux survival" in the export definition itself (`:653–658`).
- The near-unity ratio at ω=1 is explicitly not confinement (`:659–662`).
- `C_{T→H}` is "thickness end-channel conversion only" (`:631–638`).
- Bulk escape is not inferred as `1−P_T,surv−C_{T→H}` (`:265`).

This is the right reading of v10 `:620–622`: "A stationary bound state has no asymptotic outgoing current … total photon-loss probability is therefore not formed".

**Closure/interface attribution.** Signed exchange terms and their checks stay with d: "Per-face slab exchange and bulk acoustic power remain distinct signed operands … only separately normalized observable design is deferred, with no positive-loss or exclusive-mechanism inference" (`:290–294`, `:640–651`). This is consistent with S11b's discriminator 1: "Report every signed external exchange term — sink or source" (S11b spec `:458–459`). It is also consistent with the S11b standing rule that non-passive terms need a named reservoir (S11b step `:57–63`).

v10 never defined a separately normalized closure-absorption observable; §3b names only continuum conversion and bound capture. Treating that observable as a future decision is therefore honest: nothing is silently promised or eliminated.

**Gap:** the deferred packages are left without an owner, which conflicts with N10's card rule (finding 2).

### Q2. A9 and the bulk-flux FORM

**General FORM retained.**
- Amendment `:89–127` keeps the full reduced operator and vertex, the reference/left/right baselines, both-end channel spaces, derived signed currents, independent ε/η/σ grades including the mixed grade, and a transparent bindable representation.
- It forbids "A frequency-one coefficient table, opaque response symbol, fingerprint, unevaluated missing construction or hand-typed expected scaling" (`:117–120`).
- It allows transparent operator/integral forms rather than demanding elementary closed forms (`:120–123`).
- It rejects a simple Born substitute unless §2's premises are computed (`:64`). This matches v10 `:390–393`.

**Baseline/interference logic is correct.** The amendment states: "a baseline can make omitted second-order … terms enter a physical quadratic power through interference … Conversely, a computed baseline/interference disposition can establish a leading induced quadratic from first-order amplitudes without a second-order amplitude" (`:237–243`). This faithfully restates v10 `:589–593`.

It also fits the bulk physics. The far-field bulk current form is background-independent (uniform ρ_m, c_s0; c1 spec `:97–98`). Radiated components at k′≠k_in carry no zeroth-order baseline if the computed K₀ disposition gives none. The leading bulk power can then be B₀[a₁,a₁], so "quadratic power alone is not a deferral reason" is justified.

**Bulk map is labelled as new work, not inherited.** "This is an explicitly added d construction/consume-set and comparison duty, not an existing `IMPORT_KEYS` row" (`:156–157`). The sourcing is faithful:
- Disconnected half-spaces and the radiation condition (c1 spec `:98–113`).
- Face kinematics `n̂_s·v_bulk,s = V_s + J_s/ρ_m` (c1 `:154`, `:197–199`).
- Graph height versus outward displacement, and Z versus N (c1 record `:180–184`).
- In-plane reduction with depth kept (`:153–154`).

Requiring the actual c1/c2 face-drive map, "including any elimination of the independent centre displacement" (`:165–170`), is well-grounded. v10 types the closed operator over `{u,θ,e_W}` (v10 `:97–98`), while c1 insists "`ζ_c` is an independent face DOF" (c1 `:85–86`). Recording an upstream dependency finding instead of patching in d is the correct stop.

**The c1 far-field construction really is restricted.** The literal c1 excerpt evaluates `outgoing_farfield_poynting(inputs, source, qprop, qprop, velocity)` after `xreplace({q_out_k: qprop, q_out_kp: qprop})` (c1 script excerpt, original lines ~1380–1400). That is equal legs, a single propagating momentum, one impermeable velocity drive, and only first-shape cross terms. Amendment `:160–162` ("specified restricted subcase … does not establish the general solved-state map") is correct; see finding 3 (optional).

**Power orientations are correct.**
- Acoustic check: `P_into_bulk − P_outgoing − P_lateral` (`:199`).
- Traction check: `P_face + P_infinity`, with P_face as slab traction work (`:208–213`).

These match c1 record `:174–179` ("`P_face + P_∞ = 0` … define the bulk subtraction operand `B` as **minus** outgoing Poynting") and c1 spec `:321–330` (the face operand is not `½Re(δp V*)`). The source code agrees: `bulk_comparison = -outgoing_flux`. The rule that the traction identity must not be extended to permeable closure (`:214`) is right, because on permeable faces v_bulk,n ≠ V_s.

**Grade coverage is correct.** "Coverage must reach the grade carrying the claimed coefficient … If the face-side comparison at that grade needs deferred second-shape terms, mark that comparison unavailable" (`:217–224`).

This is the physically important case. At O(η²), face-side power on an incident mode that is evanescent in depth needs the second-shape Z₂. That is exactly c1's `NOT_ESTABLISHED_AT_FIRST_SHAPE_ORDER` nullspace (c1 spec `:332–338`; record `:191–193`). The far-field |φ₁|² side, by contrast, is available.

**Domain handling is correct.**
- Both-leg non-grazing validity is kept separate from omitted-flow validity (`:296–303`), matching c1 record `:146–155`, including "D⁻¹ 1/η-pole = False".
- The endpoint and partial-support rules (`:305–313`) are needed, because any radiating k′-integral ends at grazing.

**Current identity versus saved composition is described accurately.**
- The audit excerpt forms `total_current = slab_current+depth_integral*bulk_current`.
- The consumer keeps `slab`, `bulk` and `total` parts and joins `parts['slab']+parts['bulk']` against `parts['total']` (continuum_currents excerpt `:19–28`, `:55`).
- Amendment `:269–274` reports this without treating it as validation.

The saved infinite-depth current normalization exists only under `if decay_certified:` (audit excerpt original `:4353`). That is why "convergence domain" (`:278`) is a necessary item and not boilerplate: in a regime where the bulk radiates in depth, that normalization is undefined.

**Precision is declared in advance.** Amendment `:396–401` requires each observable's precision and absolute tolerance to be declared "before the acceptance comparisons", including weak coefficients and bulk contributions.

**c1 debt accounting.** It is complete for the c1 step `:99–124` statuses and the carry-forward list, except the independence-scoping status (finding 1).

### Q3. Conditional open-channel coverage and validity

**No existing clause requires a second numerical frequency.**
- v10 §3a requires "every incoming open mode" at fixed (ω, k_∥) (`:426–431`) and a derivative identity (`:454–455`), not a second frequency.
- Exploratory acceptance item 3 varies resolution, domain and regulator, not frequency.
- The thickness-repair "two-frequency endpoint pairing" is a pairing diagnostic that is still unfinished (repair report `:327–333`, `:388–390`), not an open-channel mandate. Inventory A1 says exactly this (`:119`).

The amendment correctly labels its requirement as new: "explicit proposed acceptance clarification, not a claim that v10 already specified an additional numerical frequency" (`:457–458`; also `:529–534`).

**Coverage is non-vacuous without a global theorem.**
- The object checked is v10's `C_{T→H}` on a nonempty thickness end-channel space with a nonzero incident denominator (`:460–468`).
- A symbolic route counts only on the same established physical branch; a continuation outside the physical domain does not count (`:510–516`).
- Structural absence is recorded separately from an unfound witness, and both trigger a scope decision (`:519–525`).
- Excluding leaky complex-k_n modes (`:497–499`) is correct.

**End-channel and bulk-depth availability are kept apart** (`:447–453`, `:536–543`). Physically, depth radiation from a localized interface requires |k_∥| < ω/c_s0. The saved point has none: `q_depth^2 = -k_normal^2 - 1/25` (flux report `:11–12`).

**The rest-frame reconciliation is source-correct.** The S11b step reads: "first order **fails** where `|q v₀ / ω| ≳ 1`; in the `k c_s0 ≫ ω` regime `|q| ≈ k`, so this needs `(|v₀|/c_s0)(k c_s0/|ω|) ≳ 1` — large `k c_s0/|ω|` is **necessary, not sufficient**" (S11b step `:158–161`). "Necessary" is necessary for *failure*. N11a (decisions `:183–185`) and v10 §6 (`:840–841`) compress this into a validity context, where it reads backwards.

Amendment `:328–336` restores the source meaning and uses the explicit c1 §2b pair plus the subsonic condition (c1 `:214–219`), including the strict v=0 grazing statement. It also correctly says that the clarification "establishes neither an admissible radiating domain nor its absence" and that algebraic candidates don't settle admissibility (`:334–336`; frequency-source report `:20–23`).

### Q4. Practical targets and stop conditions

These are adequate for the retained toy-model claims:
- 1% relative stability for resolved nonzero observables, plus absolute goals of 1e-4 (amplitude) and 1e-6 (current), each "in their declared channel normalization/unit frame" (`:393–401`).
- Near-zero effects stay unresolved without a denominator floor (`:417–419`).
- No unitarity is assumed (`:420–421`).
- Regulator effects must be separated before any claim about classification (`:474–479`).
- Bounded pilots require a cost/stop record (`:480–488`, `:509–518`).
- Certification theorems are not gates, but local validity is still required (`:423–428`).

These carry exploratory acceptance items 3–5 (`:37–55`) faithfully. One case-specific diagnostic should feed the per-observable precision declaration (finding 6).

### Q5. Saved four-case results and review state

These are represented accurately against the saved reports:
- ω=1: 645 unknowns and four incident columns; regulator 0.1.
- Four open transverse and zero open thickness channels; bulk-depth support closed.
- Anchoring changes below resolution.
- Coordinate differences 2.22e-13 / 2.97e-13.
- Twelve uniform controls; the reduced-coupling triplets are unestablished.
- First-jet sensitivities 8.243e-7 / 5.482e-7.
- B-branch counts: 24 of 32 owners (9 LAB + 15 material 1D), the Row46 spread of 4.67e-14, condition number 5655.

The review state is described as V/P or P, never cleared. I found no deliverable silently marked complete. Every v10 §§3–5 output family is either retained (A2–A13, C1–C4) or explicitly deferred (B2, B4–B6, A14). The §4 carrier-reconstruction item is only implicit in A2 and C3 (see coverage limitations). Omissions and fidelity issues are in findings 6–8.

---

## Findings

### Amendment

**1. Blocking (minor): c1 debt accounting omits the "independence is scoped" status, which bears directly on the new bulk map.**

*Source.* The c1 record lists the flat symbol as "ESTABLISHED — … cross-engine AGREE" (`:101–102`), and then limits that:

> "The c1 spec supplied the composition recipe and some expected structural values — rigid-shift cancellation, the flat `Z₀=ρ_m ω/q_out`, the zero-jet outcome … for THOSE objects part of the cross-engine 'agreement' is **fidelity to the supplied structure** … What IS independently confirmed is the **two-momentum DtN kernel**" (`:157–163`)

The carry-forward item at `:197–200` says the same.

*Artifact.* The amendment says: "Carry the full applicable c1 step 99–124 status" (`:175–176`). Its T7 paragraph lists "Direct c1 whole-DtN/flat-leg/ENERGY/traction/seal-5 debts and all other applicable c1 statuses listed in §2" (`:703–705`).

*Why this matters.* The new half-space solution map and the bulk-depth FORM rest on the flat outgoing symbol/propagator and the zero-jet limit. A build plan following the amendment would label the flat symbol "cross-engine AGREE" as an inherited independent confirmation. The existing "not exhaustive" / "no inherited agreement closes these debts" wording (`:175`, `:182–183`) does not catch this, because the flat symbol is in c1's *established* list, not its debt list.

*Minimal correction.* Add c1 record `:157–163` / `:197–200` to the carried status in §2 and §5. State that flat-symbol, rigid-shift and zero-jet agreement is fidelity-scoped and cannot be cited as independent validation of the half-space map or bulk FORM. Independence for those objects comes only from the blind d construction and T7.

**2. Blocking (minor): "unassigned" deferrals conflict with the family card's owner-naming rule, and the precedence table doesn't reconcile them.**

*Source.* N10: "⇒ the family's card names each deferral's owner" (decisions `:175`). N1: "A single roll-up ledger **card** closes the family" (decisions `:36–37`). N10's central in-scope question is "whether light's confinement is **unconditional**" (`:172–173`).

*Artifact.* "The family roll-up must record this work package as unassigned until an owner is actually appointed" (`:46–47`). The loss-attribution package has "no current production authorization or assignment to S11c-e" (`:624`). The precedence table (`:60–73`) has no N10 row.

*Why this matters.* The bound-capture half of N13 confinement, which is N10's headline question, could reach family closure with no owner. Under N10 the card cannot close in that state. The amendment says what the roll-up must *expose* (`:50`, `:676–688`), but not that an unowned deferral blocks card closure or requires an explicit user decision. As written, one could read "unassigned" as a valid final card state.

*Minimal correction.* Add an N10 row to the precedence table:
- N10's owner-naming duty is unchanged, and "unassigned" is an interim status only.
- The S11c roll-up card cannot close the family, or claim the confinement question is settled, while the pole package or loss-attribution question lacks a named owner, unless the user explicitly records that exception.

**3. Optional:** Name the concrete restriction of the c1 energy construction (`:160–162`, `:208–211`): equal input and output legs at one propagating `qprop`, with an impermeable single-velocity drive. That is the precise reason it cannot be the off-diagonal (k≠k′) radiation map that bulk escape from a localized interface needs. The present generic wording is adequate.

**4. Optional:** Specify how the bindable FORM (`:100–106`, `:115–123`) carries end-mode roots. The end relations are high-degree radical/polynomial eliminations (frequency-source report `:16–19`: "complex degree22 eliminations"). The export should include the defining end relations together with the outgoing/current-sign branch rule and the channel-availability and classification boundaries. Then binding a new (ω, k_∥) re-derives which channels are open rather than inheriting ω=1's census. Much of this is already implied by `:621` ("Carry valid domain and channel availability") and `:493–496`.

**5. Optional:** On radiating bulk support |q_out| ≤ ω/c_s0, so the c1 condition `|q_out·v/ω| ≪ 1` requires roughly v/c_s0 ≪ 1, which is stronger than subsonic. The amendment already requires both conditions (`:317–319`). A sentence saying so would stop "subsonic" from being read as sufficient.

### Inventory

**6. Required (minor): the MATERIAL_ADVECTED/RHO4_CONSTANT diagnostics that bear on A9/A10 are left out.**

*Source.*
- Response report `:29–33`: "At (eta,sigma)=(0.01,0.001), the new LAB/RHOBR and MAT/RHOBR field-coefficient differences are 3.22e-5; MAT/RHO4 is 8.82e-3."
- Profile report `:40–43`: "about 3.36e-5 for the RHOBR cases and 8.53e-3 for MATERIAL/RHO4".
- Profile report `:26`: the MAT/RHO4 finite amplitude change is 2.23e-5, against ≤8.9e-7 in the other cases.
- The pressure-trace repair report `:69–71` also singles out MAT/RHO4: zero raw R_N6 samples, against 20,736 in the other three cases.

*Artifact.* A3, A4 and A7 (`:121`, `:122`, `:125`) describe the four cases uniformly. A10 (`:128`) marks weak-coefficient assessment P but gives no case-specific signal.

*Why this matters.* A retained-polynomial remainder about 250× larger in one case is a direct indicator of Born-domain or resonance proximity for the first-priority weak coefficients. The amendment requires the export to state "where resonance proximity remains unassessed" (`:104–105`), and v10 excludes "resonance enhancement" from the continuum Born domain (`:398–399`).

*Minimal correction.* Record these values in A4 and A10 as open case-specific diagnostics that must be addressed when each weak-coefficient precision is declared. Do not interpret them further.

**7. Required (minor):** The same omission as finding 1 appears where the inventory lists debt: A12 "All applicable c1 step 99–124 statuses remain explicit" (`:130`) and C3 (`:151`). Apply the same correction.

**8. Minor citation fix:** In A12, the link labelled "c1 shape/power carry-forward" resolves to `[c1-retained]` → c1 record `:99` (`:282`), which is the established/owed list. The carry-forward corrections are at `:170`, as the amendment's `[c1-record]` correctly cites.

**9. Consequential:** "What to defer" (`:157–179`) and C4 (`:152`) need finding 2's owner-naming reconciliation for the pole package and the loss-attribution question.

---

## What is faithful and intentional

Most of the amendment reflects deliberate, clearly labelled scope choices that are consistent with the sources, not physical errors:
- **Pole deferral:** an explicit reduction, not a claim of completion (`:52–53`).
- **Weaker bulk-escape promise:** labelled as such (`:66`, `:680–683`).
- **Second-shape/threshold completion at the S11c-e boundary:** consistent with c1's "`O(η²)` leakage there belongs to S11c-e" (c1 spec `:338`), and labelled as extended to threshold completion (`:47–49`).
- **Non-vacuous A11 coverage:** labelled as a new requirement.

Honest labelling of the six same-author folds against G4 is a process matter outside this physics filter.

---

## Coverage limitations

- **Not supplied:** raw arrays, executable builds, the build directive and program brief, the c1 retro-review adjudication files, `S11c_c2_self_energy_fold.md`, and the full `S11c_d_sympy_builder_report.md` (only lines 488–532). Source code was supplied only as the listed excerpts. Claims about unsupplied source lines, including the meaning of `FREQUENCY_PENCIL_PLUS`, were not verified.
- **Read partially or by search only:** S11c-a spec, S11b spec outside §6 energy accounting, the full r10 review legs, POLE_HANDOFF, the inertia/mechanical/c2-trace/thickness repair reports (read only at the cited passages), the evidence JSON (status fields only), and the material-row summary. I confirmed MAT/RHO4 rows 60–62 in the summary. I did not verify the MAT/RHOBR rows 50–53 designation.
- **Carrier reconstruction:** v10 §4's "reconstruction of every c2 three-dimensional carrier" is only implicit in inventory A2 ("Computed reduction") and C3. I could not confirm its saved status.
- **What this review is not:** a document review only. It clears no build, export, cross-engine comparison or physics result, and no future control is treated as executed.