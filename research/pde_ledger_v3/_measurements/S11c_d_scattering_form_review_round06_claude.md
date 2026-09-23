**AMENDMENT: NEEDS REVISION**
**INVENTORY: NEEDS REVISION**

# S11c-d Option B amendment (draft 6) and scope inventory (draft 6): independent document review

**Reviewer:** fresh non-author Claude agent (Opus 5.5). This is a document and specification review only. I executed nothing. I recomputed nothing and ran no CAS or numerical work. All reported numbers are taken from the saved reports as written.

Both verdicts are about what these documents specify. Neither clears the later build, the export or the physics. The amendment is largely coherent and correct about the sources. It needs two small fixes before clearance (Findings 1 and 2). The inventory repeats Finding 1.

---

## What the sources establish (formed before reading the artifacts)

- **v10 (`S11c_d_SHARED_PHYSICS.md`)**
  - The canonical object is the complete two-ended S-matrix on the §1c-reduced full operator. It uses S11b-derived modal currents and computed baselines `K₀, K₋, K₊` (§2, §3a :421–477).
  - The continuum response is re-expanded to the retained `(ε,η,σ_W)` rectangle.
  - §3c (:589–603) states the baseline/interference slots. It warns that "omitted pure second-order amplitude/current terms can also enter the parent-theory `λ²` flux through baseline interference". It also says that "under a computed disposition that removes the relevant baseline/interference slots, the leading induced-field `B₀[a₁,a₁]` term does not require an uncomputed second-order amplitude."
  - §3b(i) names "continuum conversion … into the thickness continuum / bulk escape" (:485–486). §3a defines `J_H` only through end-channel currents (:459–465).
  - §6 (:839–841) repeats N11a's parenthesis "large `k c_s0/|ω|` is **necessary, ⛔ not sufficient**".
  - No clause in v10 requires a second numerical frequency.
- **S11b step record :158–163**
  - Full context: "first order **fails** where `|q v₀/ω| ≳ 1`; in the `k c_s0 ≫ ω` regime `|q| ≈ k`, so this needs `(|v₀|/c_s0)(k c_s0/|ω|) ≳ 1` — large `k c_s0/|ω|` is necessary, not sufficient."
  - So "necessary" refers to **failure** of the expansion. It is not a condition for validity.
  - In the radiating regime `|q| ≤ ω/c_s0`, so `|q v/ω| ≤ v/c_s0`. That first condition cannot fail there at small Mach number. The live restriction is the c1 §2b boundary-layer term, which diverges at grazing.
- **c1**
  - Bulk acoustics, radiation condition and disconnected half-spaces: c1 SHARED §1b :88–122.
  - Kinematic drive `n̂·v_bulk = V_s + J_s/ρ_m`: §1d :154.
  - The restricted energy route (real ω, propagating, impermeable, `Λ_X⁰=0`) is §3b :321–330.
  - The first-shape power caveat, "the `O(η²)` leakage there belongs to S11c-e", is at :332–340.
  - The step record's carry-forward (:170–200) corrects the energy orientation to `P_face + P_∞ = 0` and gives four other clarifications:
    - graph height versus outward displacement;
    - `Z` versus `N`;
    - `K_a` is Hermitian;
    - the second-shape grades are η², ησ_W and σ_W².
  - The step record also scopes grazing to non-grazing asymptotics on both legs (:146–156).
  - It lists cross-engine debts at two levels (:23–31, :99–124):
    - UNDECIDED items;
    - **UNMEASURED/DEFERRED** giant families and full per-family residuals.
- **c1 script excerpt :765–831 and :1380–1466**
  - The historical far-field route sets `q_output = q_input = qprop`, which is the diagonal-leg subcase.
  - Its Poynting keeps only the flat×flat and flat×scattered terms. There is **no scattered×scattered term**.
  - This is exactly the restricted first-order construction. It cannot supply off-diagonal leading power, even as a template.
- **c1 export census** (44 keys)
  - It contains the DtN symbol, kernel, operators, `q_out` legs and response resolvents.
  - It contains no half-space field or far-field map.
- **Saved reports**
  - At ω=1: four open transverse directions and zero open thickness directions in all four cases.
  - `q_depth² = −k_normal² − 1/25`. This is closed for every real transfer momentum at this k_∥, so it is a structural closure at that input.
  - Current ratios are 0.9999995945–1.0000001045. Every profile, first-jet and anchoring difference is below the 1e-4/1e-6 resolution.
  - The baseline contour on |ω−(1−0.01i)|=0.02 closed "without a resolved candidate" and is "not a certified empty spectrum".
  - The frequency-source candidates (≈0.2739/0.2725, bulk branch ±√5) are "algebraic end candidates … not classified physical thresholds".

---

## Answers to the five questions

### Q1. Pole deferral, continuum obligations, confinement and the attribution disposition

**Mostly yes.**

**Pole deferral.** The deferral runs coherently through every layer:
- scope: §1 table rows for v10 §0 and §3b(ii);
- output and export: §5 table rows for the pole families and pole metadata, with "Do not export zero, an empty set…" (:585);
- comparator: "Pole families are outside the reviewed compared domain; they are not residual zeros … or a fourth T7 truth value" (:661–663);
- downstream consumers: an unsupported/deferred capability failure (:628–632).

**Local scattering needs.** These are explicitly kept: "This deferral does **not** remove the resolvent, closed modes, radiation/current construction or local branch/domain checks" (:514–517). So are the nonlinearPoleV2 §3 full-pairing safeguards for end modes. The draft correctly keeps the `∂_ω` pairing (normalization) separate from the `∂_{k_n}` pairing (semisimplicity in `k_n`) (:526–531).

**Historical pole statement.** The baseline contour statement (:535–543) matches the contour report word for word in substance.

**Confinement.** The limitation is correct. Survival is renamed "continuum-channel transverse flux survival, not a complete N13 confinement answer" (:613–614). This is an intentional scope change against v10 §3b's "This is the `N13` confinement object" (:536). The draft declares it openly.

**End channels versus bulk depth.** They are kept as separate types (:588–595, :404–407). The draft states that the depth-integrated interface-normal tail belongs to `J` and "is not another loss to add" (:247–251). This matches v10 §3a :440–442 ("including its closed/nonlocal bulk contribution") and the report text "Depth-integrated normal bulk current is not depth escape".

**Closure/interface attribution.**
- The draft keeps the signed exchanges and applicable checks (:252–255, :581, :597–608).
- It routes only a *separately normalized* absorption observable to an unowned future decision.
- It neither promises that observable nor sets it to zero.
- It keeps a catch clause: "If source inspection identifies an existing mandatory operand or check, it remains current work".

This is honest. Physically, per-face slab exchange minus acoustic power into the bulk is the permeable-port form of S11b check 3 (S11b SHARED :475–487). That means the required operands already determine the interface loss at the retained order. See Finding 5.

### Q2. A9 form, bulk-flux obligation and c1 linkage

**Yes on the core. Two gaps.**

**What A9 does correctly:**
- Keeps the full reduced two-ended distorted-wave problem (:95–99).
- Keeps independent grades and the mixed grade (:100–102, :211–213).
- Keeps baselines, interference and incident denominators (:109–116).
- Keeps transparent representations, and rejects "a frequency-one coefficient table, opaque response symbol, fingerprint, unevaluated missing construction or hand-typed expected scaling" (:117–125).
- Explicitly forbids an unjustified Born shortcut (§1 table row for §3a/§3c–§3d).

**Baseline/interference test.** The draft handles this correctly. It does not assume that every quadratic needs a second-order amplitude (:219–225), and this matches v10 §3c :589–593.

**Physics note (supports the draft).** At a radiating transfer momentum `k_out ≠ k_in`, the uniform baseline has no radiating component. The leading bulk-depth power is therefore `|first-order amplitude|²`. It is supported whenever the computed baseline face drive at the incident momentum vanishes. The c1 caveat (:332–338) concerns the face-side Hermitian form on an evanescent baseline. That is where Z₂ would enter.

**Bulk-flux obligation.** The draft labels this honestly as added work: "an explicitly added d construction/consume-set and comparison duty, not an existing `IMPORT_KEYS` row" (:158–160; T7 additions at :652–657). This matches the export census, which has no field map. It also correctly states that the c1 far-field route "has a specified restricted subcase" (:163–164). The script confirms this with its diagonal legs and missing scattered×scattered term.

**Source linkage.** The draft ties the construction to:
- the same physical faces, outward normals and true-area measures;
- the graph-versus-outward and `Z`-versus-`N` corrections (:151–154);
- in-plane-only reduction with depth kept (:155–156);
- the ζ_c elimination to be located rather than assumed (:166–169).

**Power orientations.** Both are consistent with the sources:
- The acoustic balance `P_into_bulk − P_outgoing − P_lateral` (:189) agrees with S11b check 2's `P̄_bulk = ½Re(δp V*)` as power *into* the bulk (S11b :466).
- The traction test `P_face + P_infinity` (:199–202) agrees with the c1 record :174–179.

**Grade statuses.** The three statuses are well posed (:214–217): retained-model coefficient, justified leading physical coefficient, and induced diagnostic.

**Grazing.** Both-leg first-shape validity is separated from omitted-flow validity (:257–264, :276–287). Endpoints of the flux integral are handled without deletion (:266–274). This fits a real feature of the problem: the radiating support of a transfer-continuum flux always ends at grazing points on the sound cone.

**Gaps:**
- The direct c1 debt list is incomplete (Finding 1).
- The face/far-field check can be satisfied only at vanishing grades (Finding 3, non-blocking).

### Q3. Open-thickness-channel coverage and the flow-validity reconciliation

**Yes.**
- No existing clause requires a second numerical real frequency, and I found none either. v10 §3a asks only for "every reflected/transmitted open channel" (:433). Acceptance item 1 asks for the approved development example. The thickness-coordinate report's "two-frequency endpoint pairing" is a normalization diagnostic, and the inventory A1 labels it correctly.
- The draft labels its stronger requirement as new: "explicit proposed acceptance clarification" (:414–415, :486–488).
- The coverage it requires is non-vacuous. A zero-rank selector, bulk-only radiation or a current deficit is excluded (:423–425). The bulk check likewise requires "nonempty radiating support and an applicable nonzero incident current" (:304–306).
- It demands no global theorem, and it has stop rules for both the end-channel and bulk sides (:316–331, :466–484).
- The N11a / v10 §6 reconciliation (:289–297) matches the full S11b passage quoted above. It correctly declines to turn that reading into a radiating domain ("establishes neither an admissible radiating domain nor its absence").
- Algebraic candidates are kept distinct from physical channels (:296–297, :407–410, :550–552).

### Q4. Tolerances, comparisons, current meanings and stop conditions

**Adequate apart from Finding 2.**
- The "roughly 1%" target, the 1e-4/1e-6 goals and the rule against applying them to "dimensionful source coefficients or arbitrary matrix entries" (:354–358) are consistent with the saved reports ("original amplitude 1e-4/current 1e-6 reporting resolution").
- The near-zero handling and the prohibition on imposing unitarity by assumption (:374–378) are correct.
- **Gap:** the consolidated text drops acceptance item 4's requirement to *pre-declare* precision and an absolute tolerance per observable before the acceptance comparisons. The new A9, A11 and A12 observables therefore have no declared frame.

### Q5. Saved-result representation and silent omissions

- The four-case real-frequency state is represented accurately in both artifacts. This covers the current record, the continuum report, the near-unity ratios, the closed bulk kinematics, the contour status, the Lean NP1–NP4 scope and the material-row counts (B1: 9+15 computed, 8 remaining; this matches the material-1D and remainder-row reports).
- Nothing is silently declared complete. Each item carries an explicit V/P or P status.
- **Omissions present in both the retained and the deferred lists:**
  - the c1 UNMEASURED debt families (Finding 1);
  - acceptance item 4's pre-declaration requirement (Finding 2).

---

## Findings

### 1. BLOCKING: the direct c1 debt list omits recorded UNMEASURED/DEFERRED families

**Amendment :171–173:** "Direct c1 use carries its actual unresolved cross-engine premises: the whole-form `dtn_operator`, off-diagonal flat-resolvent momentum leg labels, ENERGY audit, `t_s` traction leaf and seal-5 density representation." The same five-item list is repeated at :660–661, inventory A12 (:131) and inventory C3 (:152).

**The c1 record also lists:**
- "the **4 giant families** (PERMEABLE_PORT_HERMITIAN, PERMEABLE_DISSIPATION_VS_OMEGA_TAU, UNIFORM_LIMIT_S11CC1_OPERAND, UNIFORM_LIMIT_RESIDUAL) + the **full per-family symbolic residual** are **UNMEASURED — DEFERRED**" (c1 step :27–29, :118–124);
- "The Hermitian/reactive dissipation parts AGREE only at c1's first order … full symbolic DEFERRED" (:105–107);
- the HOMOGENEITY keying gap, "⛔ NOT AGREE-by-inheritance" (:116–117).

**Why it matters.** The amendment's retained signed-balance duty (:252–255) and the S11b check-3 method it invokes pair per-face slab exchange with the permeable-port form. That form is exactly PERMEABLE_PORT_HERMITIAN. The existing sentence "identify which claimed coefficient depends on each unresolved premise" (:177–178) only reaches the enumerated five. A balance or attribution claim could therefore rest on an UNMEASURED c1 family without that status being visible.

**Minimal fix.** Make the list non-exhaustive by pointing at c1 step :99–124. Add the four giant families, the deferred per-family residual (REP_INVARIANCE/CONTROL_*/DEGENERATE/DIMENSIONS) and the first-order-only Hermitian/reactive agreement. Optionally also carry the `K_a`-is-Hermitian label correction (:188–190) wherever dissipation objects are reused. Apply the same fix to inventory A12 and C3.

### 2. BLOCKING (low cost): consolidating the acceptance policy drops per-observable pre-declared precision

**Acceptance item 4 (:42–49):** "Set the intended reporting precision before the acceptance comparisons… Record absolute tolerances in the declared unit frame too."

**Amendment table (:73):** "Consolidate its practical numerical policy in §3 below." §3.1 (:354–369) keeps the legacy goals, restricts them to the channel frame, and forbids applying them to dimensionful coefficients. It never requires *new* observables to declare a precision or absolute tolerance before comparison. The new observables are the A9 weak coefficients, `J_H` at any open-channel pilot, and the bulk-depth flux FORM with its dimensionful source units. Every saved ω=1 effect is already below resolution, so these new claim-carrying observables are exactly where post-hoc tolerance choices would decide the outcome.

**Why the existing wording does not cover it.** "Small effects … need tighter targeted checks" says tighter checks are needed but does not fix *when* or *in which frame* the tolerance is set. It is also unclear whether "consolidate" supersedes item 4.

**Minimal fix.** Either state that acceptance items 1–5 remain in force except where explicitly changed, or add one sentence: each new reported observable declares its precision target and its absolute tolerance, in its own unit and normalization frame, before its acceptance comparisons.

### 3. Non-blocking: face/far-field check (ii) can be discharged at a vanishing grade

Check (ii) (:182–195) compares face-delivered acoustic power with far-field flux at "supported grades". When the incident baseline has no radiating component, both operands are zero at the η⁰ and η¹ grades. The claimed flux lives at the quadratic grade. If the computed baseline face drive at the incident momentum is nonzero (evanescent), the face-side operand at that grade needs second-shape face terms (c1 :332–338), so check (ii) is structurally unavailable there.

The general clause "Non-vacuous coverage must check each claimed supported bulk-flux FORM" (:303–306) arguably covers this.

**Suggested clarification:**
- State that check (ii) counts as coverage only at the grade carrying the claimed coefficient.
- When the face route is unsupported at that grade, name the independent route: the blind Wolfram far-field construction or T7. If no independent route exists, stop.

### 4. Non-blocking: face-drive map with no ζ_c elimination

The draft correctly requires locating "any elimination of the independent centre displacement" and forbids "an invented extra degree of freedom" (:166–169). The v10 consume-set operator is over `{u,θ,e_W}` (v10 :97–98). c1/S11c-a forbid setting `ζ_c = 0`, and c1 says a curvature-induced `δW↔ζ_c` coupling is a computed block (c1 :192–195).

If inspection finds that no upstream map determines ζ_c, the draft should say this is recorded as an **upstream c1/c2 finding** that also bears on the end-channel S-matrix. It is not a d construction task. The general "report as outstanding" clause (:707–708) covers the A12 status, but not the upstream escalation.

### 5. Non-blocking: the scope of the attribution deferral

The per-face slab exchange (S11b check 3's `−½ΣRe[(δp+Λ_X𝒜)V* + μ J*]`) and the acoustic power into the bulk are both required operands. Their difference is the interface port dissipation. The draft's catch clause (:603–605) keeps this difference as current work if it is found to be mandatory.

A one-line note would stop a reader from mistaking the future-scope item for deferral of that difference itself: what is deferred is only a separately *normalized* observable and its design.

### 6. Non-blocking: self-containment of the governing amendment

The amendment's obligations are indexed by inventory IDs (A9, A11, A12; e.g. :42, :497, :710). The inventory is a draft state document. The governing amendment should define these obligations in its own text, or pin the exact inventory version.

The amendment also does not carry inventory A1/A2's unestablished post-repair endpoint-pairing and current-adjoint normalization coverage (thickness-coordinate report :291–292, :388–390). That gap affects the saved current denominators that §3.2 relies on for its characterization.

### 7. Inventory: cosmetic

There is a duplicated heading: "## Evidence boundary and preservation## Evidence boundary and preservation" (inventory :214).

**Inventory fidelity otherwise.** I checked the following against their sources and found them accurate: A3/A4 against the response and flux reports, A5 against the coordinate report, A6 (coupling-triplet gap honestly marked P), A7/A8, A10 against the domain report, B1–B4, C2 against the Lean closure record, and the Wolfram approval chain. The inventory's only blocking issue is the shared debt-list omission in Finding 1.

---

## Coverage limitations

- **Not in the packet:**
  - the S11c-c2 spec and step record, so the ζ_c/face-drive map and c2 debt details cannot be checked;
  - the build directive and program brief (the amendment says it does not reconcile them);
  - `S11c_d_focused_completion_plan.md`;
  - round 1–5 drafts and reports;
  - the c1 reconcile and retro-review adjudication files;
  - `DEFERRED_HEAVY_RUNS.md`.
- **Partial excerpts only:** the builder report (lines 488–532), CLAUDE.md (47–65, 117–180) and the c1 script (765–831, 1380–1466).
- **Commit hashes** (399a8516, 96fa1f27, 0d77af53, etc.) and the amendment SHA in the evidence index cannot be verified from inside the packet.
- **Read only in part or by search:**
  - S11b SHARED: lines 420–500, plus grep;
  - S11c-a SHARED: grep only;
  - the thickness-coordinate report: sections only.
- **Not read:** POLE_HANDOFF, the nonlinear-pole-repair, inertia, mechanical, c2-trace and Wolfram-audit reports, the LAB operator/row-sensitivity reports, the material-input/1D-review reports, the frequency-matrix report, the r10 Grok report beyond its verdict, and the material-row-summary JSON. Findings that depend on them are not claimed.
- Nothing numerical was recomputed. Saved values are taken as reported.