# Independent review: S11c-d Option B scattering/FORM amendment (draft 4) and scope inventory

**Reviewer:** fresh non-author Claude agent. **Method:** document review only. I read the governing sources before the artifacts. I did not recompute anything, run anything or change any file. A CLEAR verdict here would cover only document fidelity and physics-specification adequacy.

## Verdicts

- **AMENDMENT: NEEDS REVISION.** One blocking finding, F1.
- **INVENTORY: NEEDS REVISION.** It repeats F1 in A12 and C3 (finding I1). I2 is a smaller done/remaining inconsistency.

Everything else I checked in the amendment is sound, or needs only the optional clarifications listed below. That includes the pole deferral, the kept pole claims, the S11b/N11a validity reconciliation, the §3.3 open-channel treatment, the stop rules, and the historical search and survival statements.

---

## Source-grounded account (formed before reading the artifacts)

- **Scattering channels.** v10 §3a defines channels as left/right modes of the *full end pencils* `𝓛_±^full` at fixed `(ω,k_∥)`. It requires all "actually open" channels (v10:374–376, 426–434). The modal current is "the S11b quadratic energy current, evaluated on the §1c-reduced closed operator … including its closed/nonlocal bulk contribution" (v10:440–442).
- **Continuum conversion.** v10 §3b(i) names "transverse → the thickness **continuum** / bulk escape (the §3a amplitude projected on radiating/continuum channels…)" (v10:485–486). §3a/§3b define no separate out-of-plane bulk-radiation functional.
- **The reduced operator is slab-only.** The consumed closed operator is "over `{u,θ,e_W}`" (v10:97–98). The bulk has been eliminated through c1's per-face DtN. Unreduced rows are "operands of the §1c reduction **only**" (v10:350–352).
- **c1 bulk.** Bulk acoustics and the radiation condition are supplied. The two half-spaces are "disconnected … the impedance is a **per-face** operator" (c1 excerpt:17–20, 28–41).
- **Rest-frame validity.** c1 §2b sets out the away-from-grazing pair of conditions plus the subsonic condition, and says "the requested grazing behaviour … is the strict `v_bulk_normal_0=0` result". It adds that "this domain is recorded by the step record, not carried as a term" (c1 excerpt:56–63).
- **The S11b standing-limit passage** (S11b step:158–164) reads: "first order **fails** where `|q v₀/ω| ≳ 1`; in the `k c_s0 ≫ ω` regime `|q| ≈ k`, so this needs `(|v₀|/c_s0)(k c_s0/|ω|) ≳ 1` — large `k c_s0/|ω|` is **necessary, not sufficient**". So "necessary" is a necessary condition for *failure* in the evanescent high-k estimate. N11a (decisions:184–186) and v10 §6 (v10:839–841) compress this into a parenthetical attached to the *validity* domain, which invites a misreading.
- **Frequencies.** No existing clause I found requires a second real frequency:
  - v10 §3a is posed "at fixed `(ω,k_∥)`".
  - Exploratory acceptance item 1 names "the approved development example".
  - Retained contract item 3 says only "Test actual reference/end channel availability … Do not manufacture an incident channel" (builder excerpt:21–23).
  - The thickness-coordinate report's "retained two-frequency endpoint pairing" (report:327–333) is an independent-frequency *normalization diagnostic*, not a second physical scattering point.
- **Saved state at ω=1.** "four open transverse directions and zero open thickness directions … `q_depth^2 = -k_normal^2 - 1/25`, with empty real propagating domains" (flux report:9–12). The baseline contour "no candidate is resolved … on |omega-(1-0.01i)|=0.02 … not a certified empty spectrum" (contour report:23–26).

---

## Answers to the five questions

### 1. Pole deferral, the retained obligations and closure attribution

**The pole deferral is coherent through every layer:**
- scope (amendment:47, 50)
- output/export (table at 442–443: "Do not export zero, an empty set, an empty result-bearing container")
- comparator (503–505: "not residual zeros, evidence of equality or a fourth T7 truth value")
- downstream (482–486, 487–494)

**The local scattering needs are kept** (372–389). This includes the correct separation of the `∂_ω𝓛` normalization pairing from the `∂_{k_n}𝓛` semisimplicity pairing (384–389). That matches v10:447–455 and nonlinearPoleV2 §3.

**Historical statements are limited accurately:**
- 393–401 matches contour report:23–27.
- 414–418 matches POLE_FIDELITY_REVIEW:7–9 ("not the analytic existence theory or a physical S11c pole/scattering calculation").

**Survival is correctly restricted** (467–477). A stationary S-matrix cannot contain capture into a bound mode (v10:513–515, 620–622). So describing it as "continuum-channel transverse flux survival" is a faithful restriction, not a physics change.

**End-channel and bulk-depth obligations are separately typed** (445–452). This is correctly labelled as "clarification … and an addition to the detailed deliverable list" (168–169).

**Closure attribution** (438, 454–465) keeps signed terms. It explicitly neither promises nor zeroes a separate absorption observable, and it preserves any "existing mandatory operand or check".

Weaknesses: F1 (bulk reconstruction), and the optional points O2 and O3.

### 2. A9 as a general, bindable weak-order FORM

A9 does retain this. Item 1 keeps the full reduced operator, the vertex, the three baselines and both end channel spaces (78–82). Items 2–3 keep independent grades, the mixed grade, the actual baseline, interference and current forms, and the incident denominator (83–99). Item 5 keeps unit, source and domain joins (109–112).

It avoids both failure modes the prompt names:
- Born shortcut: "does not license a simpler Born matrix element without §2's computed reduction premises" (49). This matches v10:390–393.
- Placeholders: "A frequency-one coefficient table, opaque response symbol, fingerprint, unevaluated missing construction or hand-typed expected scaling is not the general computed FORM" (102–105). It is also reasonable in not demanding an elementary closed form (105–108).

**The bulk-flux obligation** is correctly presented as a requirement, not a result: "No numerical value or closed expression for the flux is supplied here" (130); "not established by the existing real-frequency reports" (170–171).

**However, the half-space reconstruction is not grounded in the supplied source scope. This is F1.**

### 3. Open-thickness coverage and the validity regime

**Frequency obligations.** The amendment correctly says no existing clause demands a second frequency (248–250, 345–346). It labels the new requirement as new: "This is an explicit proposed acceptance clarification, not a claim that v10 already specified an additional numerical frequency" (273–274). I found no clause that contradicts this.

**Non-vacuous coverage without a global theorem.** Coverage is required through either:
- an independent symbolic check on an established physical nonempty branch, or
- one bounded witness.

Both come with an explicit structural-absence-versus-no-witness distinction and a stop for a scope decision (325–343). The root-mechanism diagnosis (286–295) is physically apt. It covers sub-threshold evanescence, bulk leakage, *retained closure dissipation* (S11b §4 closures with finite `τ_I` are complex, S11b excerpt:159–161) and regulator effects. It forbids "setting a forbidden closure parameter to zero".

**End channels versus bulk depth.** These are kept separate (263–266, 352–358). Algebraic candidates are correctly denied physical status (164–166, 266–269; frequency-source report:20–23).

**S11b/N11a reconciliation.** Lines 158–166 read the S11b step correctly in context (quoted above). The operative conditions are the explicit c1 §2b pair plus the subsonic condition, with the grazing limit being strict rest frame (146–156). The rule that "an integral whose support reaches that excluded region needs an explicit strict-rest-frame interpretation" (151–154) is the right answer to c1 §2b's non-uniformity near grazing.

Optional: O1 (whether "admissible" means conditional or verified).

### 4. Practical targets, current meanings and stop conditions

These are adequate for the retained toy-model claims:
- The 1% relative target plus declared absolute targets are confined to the channel normalization frame (213–217). This matches exploratory acceptance items 4–5.
- Near-zero results stay unresolved; there is no denominator-floor inference (233–235).
- Unitarity must not be assumed (235–237).
- Weak comparisons must be like-for-like in retained grades (198–201).
- End-channel `C_{T→H}` and bulk-depth escape are not interchangeable (355–358).
- Bulk escape must not be inferred as `1−P_surv−C_{T→H}` (140).

One meaning needs tightening: O2 (double-counting risk in the bulk-current wording).

### 5. Saved results and completeness

The amendment's statements about saved results are accurate:
- 255–261 matches the flux and continuum reports.
- 473–476 is logically correct.
- The material-row accounting is consistent with the sources. The inventory's "24 of 32" is 9 LAB rows plus 15 material 1D rows; the remaining 8 material 2D rows are rows 60–63 and 50–53 in `material-row-summary.json`.

Cross-checking the v10 §§3–5 emit lists against amendment §5, I found nothing silently marked complete and no output list missing an item:
- The strong-edge obligation is kept at 203–207.
- The profile moments and c2 carrier reconstruction are kept at 440 and 498–499.

Contradictions: only F1's "reduced problem" wording, which is also repeated in the inventory (I1).

---

## Findings

### F1 — BLOCKING (amendment §2 and §5; inventory A12 and C3): the half-space field reconstruction has no source, and no check ties it to the operator

**What the artifact says:**
- "The half-space field-reconstruction map is a construction operand of this reduced problem … unreduced three-dimensional content cannot bypass v10 §2's representation rule" (amendment:127–129).
- T7 coverage lists "reconstruction from those actual reduced rows, the half-space field-reconstruction map and its bulk-depth flux" as items that "**remain** in T7 coverage" (498–502).

**Why this is a problem:**
- The reduced rows are slab-only. The closed operator is "over `{u,θ,e_W}`" (v10:97–98), with the bulk eliminated through the c1 per-face DtN (c1:22–26).
- No half-space solution operator is in the v10 §1a consume set (v10:97–125). So bulk fields in `(y_n,w)` cannot come from "this reduced problem". They require a new d construction, or a new import, from c1 §1b's supplied acoustics and radiation condition, evaluated at the first-shape-order *physical* (shifted and tilted) faces (N3; decisions:62–70).
- This is exactly where a known defect occurred. The native Mathematica c2 instrument "assigns its physical-face pressure response to a reference-pressure slot and shifts it again", which gave a "nonzero first-order height error on both faces" (Wolfram audit:16, 33–35).
- The checks the amendment lists would not catch an inconsistent reconstruction:
  - The S11b pressure-work discriminator is (correctly) limited to its impermeable, `Λ_X=0` subcase (137–139).
  - A slab-versus-bulk balance whose difference includes the unconstructed, attribution-deferred closure absorption (438) cannot separate reconstruction error from interface loss.
- "Remain in T7" also describes new coverage as if v10 §7 already contained it.

**Minimal correction.** In §2 (and mirrored in inventory A12/C3):
1. Name the source of the map. It is the c1 §1b supplied bulk acoustics and radiation condition (or an identified c1 export row), applied at the same first-shape-order physical-face geometry and reference/physical pressure convention as the repaired c1/c2 closure. Only in-plane content is §1c-reduced; the depth coordinate is not.
2. Record it as an added constructed or consumed operand for the build-directive consume-set/`IMPORT_KEYS` gate.
3. Require a literal face-trace consistency residual: the reconstruction's face traces versus the bulk-closure content of the reduced closed operator on the solved state.
4. Require the independent face-power ↔ far-field flux identity as the measure/normalization check.
5. Change "remain in T7 coverage" to "are added to".

**Why existing wording does not cover it.** "construction operand of this reduced problem" states a source that cannot supply the object. "sign, measure and closure checks" names no check that isolates reconstruction inconsistency while closure attribution is deferred.

### F2 — Not blocking. See the optional items O1–O4.

---

## Inventory findings

### I1 — BLOCKING, same issue as F1
- A12: "Construct the half-space fields/flux from the reduced problem" (inventory:127).
- C3: "reconstruction from those rows including the half-space field map" (inventory:148).

Both need F1's correction. A12's remaining items should add: source/consume-set reconciliation, the physical-face trace-consistency residual, and the face-power/far-field identity.

### I2 — Not blocking, but should be fixed: A2 overstates normalization status
- A2 says "native current normalization … exist" (inventory:117) without the caveat that A1 and A4 carry.
- The thickness-coordinate report says "Fresh endpoint sources, independent-frequency pairing and full current/adjoint normalization remain pending" (report:291–293). It ends with pairing still running (388–390).
- A4 does carry the caveat ("A1's unestablished endpoint-pairing coverage also applies to these current denominators and normalization claims", inventory:119).
- Correction: add the same cross-reference to A2.

### I3 — Optional
- A1's "found the pressure-trace defect" should say the defect was in the **native Mathematica c2 N6 instrument** (audit:8, 16). This affects which review (C2 versus C3) owns the repaired evidence.
- The blank line at inventory:128 detaches A13 from table A when rendered.

**Other inventory checks that passed:**
- Source fidelity of A3–A8, A10, B1–B4, the review-status rows, the cost bands and the retained/deferred assignment. Spot-checked figures all match:
  - 645 unknowns and regulator 0.1
  - 2.22e-13/2.97e-13 (coordinate)
  - 8.243e-7/5.482e-7 (first-jet)
  - rank 645, condition 5655, residual 9.38e-16 (LAB response)
  - 4.67e-14 (row46 sensitivity)
  - 639.84 s (contour midpoints)
  - 12.88 s and 31.29 s (LAB row/response timings)
  - about 53 and 10 minutes (Wolfram repair production and validation)
- The inventory is honest that A6's `K_uniform` triplets, A7's discriminator pairs and A8's §5a triplets are unestablished. Those match v10:782–783, 796–803 and 745–752.
- Costs are clearly labelled as planning judgments.

---

## Optional improvements (not blocking)

- **O1 (§2:172–176 and §3.3:305–308).** "Established admissible domain" could be read as requiring the c1 §2b inequalities to be *verified*. `v_bulk_normal_0` has no approved value, and c1 §2b says the domain is "recorded … not carried as a term". State that admissibility for the non-vacuity check means: the approved parameters, real propagating support, and retained-order/regularity, with the rest-frame conditions carried as explicit conditional labels. This keeps the stop rule from firing permanently.
- **O2 (140–143).** State that v10's end-channel modal current already includes its "closed/nonlocal bulk contribution" (v10:440–442). The depth-integrated interface-normal bulk tail current is therefore a *component* of `J_H`/`J_T`, not an additional loss. Only the outward half-space flux is an additional loss. As written, "do not count … as the same contribution" could invite double counting.
- **O3 (53, 438; inventory A13).** Pin the d-level source of the "retained" signed-balance duty:
  - acceptance item 2 (current and boundary checks);
  - the v10 §3a current identity;
  - the S11b discriminators, applied as method rather than as imported tasks.

  Also note that the signed per-face slab-exchange-minus-bulk-power operand can be computed from face data even while the design of a separate absorption observable is deferred.
- **O4.** N10 asks the family card to name "each deferral's owner" (decisions:175). The pole package deliberately has no owner (33–36). That is acceptable as long as the S11c family roll-up records it as unowned.

---

## Coverage limitations

- S11b spec, c1, S11c-a, the builder report and CLAUDE.md were supplied only as verbatim excerpts. In particular, c1's full perturbed two-face closure and any c1 export row that might supply a bulk solution map were not available. F1's claim that such a row is absent from the consume set rests on v10 §1a.
- The build directive and program brief were not supplied. The amendment itself defers line-by-line reconciliation of those (60–64).
- From the reports alone I could not confirm whether the four-case real-frequency normalization used post-thickness-coordinate-repair pairing (see I2).
- Partially read or grep-only: r10 Grok review, POLE_HANDOFF, the inertia/mechanical/c2-trace repair reports, and most of the evidence JSON.
- No numerical arrays or code were inspected. The reported numbers are taken as saved claims, not validated.
