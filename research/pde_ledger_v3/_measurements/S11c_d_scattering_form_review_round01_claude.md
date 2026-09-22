# Independent review: S11c-d Option B scattering/FORM amendment and scope inventory

**AMENDMENT: NEEDS REVISION**: one blocking finding (F1) and six non-blocking corrections (F2–F7).
**INVENTORY: NEEDS REVISION**: four fidelity findings (F8–F11). None is a physics blocker, but they misstate done/remaining status or assign work to the wrong branch.

These verdicts cover document fidelity and whether the physics specification is adequate. I did not rerun or validate any reported computation. This review does not clear the later build, export or physics. Historical CLEAR records are treated as historical review status only.

## Source-grounded account (formed before reading the artifacts)

- **v10 sets these requirements.** It requires a complete two-ended S-matrix over "every incoming open mode … every outgoing open mode" at "fixed `(ω,k_∥)`" (`S11c_d_SHARED_PHYSICS.md:423-431`). Its only defined conversion observable is the end-channel thickness current, `J_H,out[a] ≡ Σ_e s_e 𝓙_n,H,e,out[S_{H,e←T}a]` (`:464-465`).
- **Bulk escape is named but never defined as a functional.** §3b(i) names "the thickness **continuum** / bulk escape" (`:485-486`), and §0 names "(into the thickness continuum / bulk escape)" (`:48-49`). No bulk-escape flux object is specified anywhere in v10.
- **Validity limits carried forward.** The continuum Born domain "excludes thresholds, resonance enhancement, modal-gap closures" (`:397-399`). Every result inherits the c1/S11b smallness domain (`:839-841`). c1 §2b says the "grazing behaviour (`q_out→0`) is the **strict `v_bulk_normal_0=0` result**" (`S11c_c1_SHARED_PHYSICS.md` packet `:56-57`).
- **Pole contract.** nonlinearPoleV2 replaces v10's unrestricted projector prescriptions. Its §3 also governs modal normalization: "A singular or unresolved D cannot be made invertible by dropping basis directions or assigning a normalization" (`:93-95`). It says "Simple modes are the m=1 case with the derivative normalization retained from v10" (`:100-101`).
- **What the saved results establish at ω=1.** All four cases have zero open thickness directions, and bulk-depth radiation is closed ("`q_depth^2 = -k_normal^2 - 1/25`", `S11c_d_remaining_case_flux_report.md:11-14`). Because the −1/25 does not depend on `k_n`, bulk radiation is closed for the entire scattered spectrum at this `(ω,k_∥)`. The end-channel conversion is therefore kinematically empty at the development input.
- **Frequency dependence.** "All 80 nonlocal rows depend on frequency; a fixed-frequency matrix cannot stand in for the frequency pencil" (`S11c_d_frequency_source_report.md:7-8`).

## Answers to the five questions

**Q1: Pole deferral.** The deferral is mostly coherent across scope, emission, export and downstream use:
- The §1 table (`amendment:37-45`) defers the pole families.
- The §5 table forbids "zero, an empty set, an empty result-bearing container or a generic pole function" (`:263`).
- A consumer asking for pole results gets an "explicit unsupported/deferred capability failure" (`:271-272`).
- The comparator treats pole families as outside its domain, not as residual zeros (`:277-279`).
- Local scattering needs are explicitly kept: "does **not** remove the resolvent, closed modes, radiation/current construction or local branch/domain checks … An unresolved local singularity remains a limitation" (`:215-218`).
- The historical search statement (`:222-227`) matches the contour report (`S11c_d_frequency_contour_refine_report.md:23-27`).

Remaining gaps are non-blocking:
- The pole-contract row could be read as dropping nonlinearPoleV2's normalization safeguards for end modes (F4).
- The comparator sentence could be read as narrowing T7 coverage away from the reduced-row join (F3).
- The Born-domain premise loses one of its routes of support (F5).

**Q2: A9 FORM.** Adequate as written.
- §2 items 1–5 (`:62-94`) keep: the full reduced operator and vertex, reference/left/right baselines, both-end channel spaces and derived currents; independent ε/η/σ grades including the mixed grade; the amplitude components, incident denominator and baseline/interference/quadratic operands; and transparent, bindable representations.
- Item 4 excludes "A frequency-one coefficient table, opaque response symbol, fingerprint, unevaluated missing construction or hand-typed expected scaling" while allowing "a permitted transparent operator/integral representation."
- The table row at `:39` blocks an unjustified simple Born substitute ("does not license a simpler Born matrix element without §2's computed reduction premises"). This matches v10 `:390-393`.
- An optional clarification about the step's distributional transform versus the positive regulator is in F7.

**Q3: Open-thickness coverage.**
- I found **no existing clause requiring a second numerical real frequency.** The closest clauses are v10 §3a "at fixed `(ω,k_∥)`" (`:423-424`), acceptance item 1 "the approved development example" (`EXPLORATORY_ACCEPTANCE.md:27`), and retained item 3 "real continuum frequency … Test actual reference/end channel availability before attempting flux normalization" (`builder_report` packet `:19-23`).
- The amendment correctly labels its non-vacuous coverage requirement as proposed and new (`:167-168`, `:196-198`), and does not demand a global theorem (`:133-138`, `:188-194`).
- It correctly separates end-pencil thickness channels from bulk-depth radiation (`:157-160`).
- However, its coverage clause lets a *radiating* witness satisfy a check on the *thickness projection*, and leaves the admissibility of the symbolic route ambiguous. That is blocking (F1).

**Q4: Targets, current meanings, stop rules.** Mostly adequate.
- Absolute goals are tied to declared frames (`:111-115`).
- Unitarity is not imposed (`:129-131`).
- Near-zero results stay unresolved (`:127-129`).
- The stop rules are explicit (`:188-194`, `:320-322`).
- The consolidation silently drops two acceptance clauses (F2).

**Q5: Representation of saved results.** The amendment's use of the flux, continuum and frequency-source reports is accurate (`:149-163`). Nothing required in v10 §§3a, 3c, 3d or 5 is dropped from both the retained and deferred lists, with the caveat in F3. The inventory's issues are in F8–F11.

## Findings

### F1 (BLOCKING, amendment §3.3): a radiating witness can pay for thickness-projection coverage; the symbolic route's admissibility is ambiguous

- **Source.** v10's conversion functional counts only end-channel thickness current (`:464-465`). Bulk escape is named (`:485-486`) but has no defined functional or measure.
- **Artifact.** The requirement is to "Give the generic thickness projection, its current normalization and the weak-order FORM a **non-vacuous check** on an admissible nonempty channel domain" (`:170-172`). Yet the preferred witness needs only "outgoing thickness/radiating channel" (`:177-178`) with a "relevant radiation measure" (`:181`). The §5 table retains no bulk-escape family (`:259-264`).

**Why this is wrong.** A point where bulk radiation is open (|k_∥| < ω/c_s0) but no thickness end channel is open never exercises the thickness projector or its current normalization. Bulk loss there shows up only as a current-balance deficit, and that deficit is mixed with dissipation and truncation remainder (the amendment itself says so at `:129-131`). So "radiating" coverage could be recorded against the object it never tests.

**Two further physics premises are missing for any radiating witness:**
1. **Grazing directions.** The scattered bulk spectrum always contains grazing directions (q_out → 0). c1 §2b limits those to the strict rest-frame (v_bulk_normal_0 = 0) result (packet `:52-57`). The generic "away from a threshold" wording (`:184-185`) cannot be met pointwise across the whole scattered spectrum.
2. **Leaky end modes.** If the incident transverse or outgoing end modes fall inside the bulk cone, they become leaky (complex k_n). v10 §3a's flux-normalized two-ended S-matrix is then undefined for them. Retained item 3 already says "Do not manufacture an incident channel."

**The symbolic route is ambiguous.** Bullet 1 requires "an admissible nonempty channel domain" (`:171`). Bullet 4 then allows "A valid independent generic-domain check" after an admissible witness has failed (`:188-192`). That reads as permitting a check outside the approved or validity domain.

**Minimal correction:**
- (a) Require the non-vacuous check of the thickness projection and current to use a nonempty *thickness end-channel* domain. A bulk-only witness may support only a separately labelled statement.
- (b) Either state that bulk escape is visible only as the current-balance deficit and is not a retained, claimable FORM, or introduce a bulk-escape flux functional with its far-field measure. The second option is a NEW requirement and must be labelled as such.
- (c) At any radiating witness, require that the incident and outgoing end channels stay real-k_n (not leaky), and carry the c1 §2b grazing caveat.
- (d) Say that a symbolic check must also be on an admissible physical domain. A check on a mathematical continuation validates construction only, and the handoff must say so.

**Why existing wording does not cover it.** "Admissible" and "relevant radiation measure" never tie the witness to the object being checked. §3.2 separates the two channel types but §3.3 then merges them.

### F2 (non-blocking, consolidation fidelity): two acceptance clauses are silently dropped

- **Source.** "Small effects that carry the physical claim need tighter targeted checks" (`EXPLORATORY_ACCEPTANCE.md:48-49`). "Include a different numerical route where it is informative and affordable, targeting the dominant uncertainty" (`:40-41`).
- **Artifact.** §3.1 (`:117-131`) contains neither, although the table says it consolidates and supersedes that policy (`:43`).
- **Why it matters.** Weak-contrast leakage is exactly the kind of small, claim-carrying effect those clauses protect. As written, a sub-resolution witness would be recorded as "unresolved" rather than wrong, so it does not produce a false claim, but the requirement is quietly weakened.
- **Correction.** Restore both sentences.
- **Related fidelity nit.** `:112-113` calls amplitude 1e-4 and current 1e-6 "approved". The reports call them "declared" or "predeclared" (`first_jet:24`, `domain:40-41`), and the plan that set them is not in the packet. Say "declared", or cite the approval.

### F3 (non-blocking, clarification): T7 coverage of the reduced-row join

- **Source.** v10 §7's load-bearing comparator residuals include "the §1c **reduced-operator and reduced-kernel join**" and modal currents (`:879-885`). This is the control on Fourier and (2π) conventions.
- **Artifact.** Table row §§4/7/8: "Revise … comparator coverage … according to §5 below" (`:41`). §5 then keeps T7 only "for the retained scattering/FORM objects" (`:275`).
- **Why it matters.** The general clause "All other physics and review duties remain" (`:30`) probably preserves the join. But the §5 wording invites a narrowing that would let both engines agree on the same wrong convention.
- **Correction.** Add one sentence stating that the v10 §7 reduced-operator/kernel both-operand join, the reconstruction and the modal-current joins stay in T7 coverage.

### F4 (non-blocking, clarification): nonlinearPoleV2 normalization safeguards for scattering end modes

- **Artifact.** "For this step, replace the production mandate with the limited historical/claim safeguards in §4" (`:42`). §4.2 applies the distinctions "when interpreting existing records" (`:228-233`).
- **Source.** nonlinearPoleV2 §3 (`:93-95`, `:100-103`) and v10 §3a's degenerate-mode matrices (`:453`).
- **Why it matters.** New end-mode construction at a witness point needs the full-pairing and singular-pairing rejection for degenerate or defective k_n modes. The pole repair report says the existing code keeps these gates (`:55-57`), but the amendment text does not keep them for *new* scattering work.
- **Correction.** State that nonlinearPoleV2 §3's full-pairing and singular-pairing requirements continue to govern end-mode modal normalization in retained scattering. Also state that any later pole contract may not weaken nonlinearPoleV2.

### F5 (non-blocking, clarification): support for the Born-domain premise once poles are deferred

- **Source.** v10 `:397-399` (the Born domain excludes resonance enhancement).
- **Artifact.** The amendment says "away from a threshold/resonance/gap closure" (`:184-185`) but never says what evidence establishes that at the reported points once the pole tool is deferred.
- **Correction.** Name the local evidence to be used: for example, finite-system conditioning and a selected frequency-sensitivity comparison of the reported observable. The historical contour may be cited only as limited local context, and only within its §4.1 caveats.

### F6 (optional): completeness of the historical search statement

- The contour report also excludes "a contour-interior holomorphy/exceptional-locus proof" (`:25-27`).
- The search ran on the finite, positive-regulator, discretized pencil.
- The frequency-source report says "80 records contain nonanalytic constructors. They still need an explicitly chosen analytic continuation" (`:12-13`).
- Adding these to `:222-227` would make the statement complete. The existing text is accurate as far as it goes.

### F7 (optional): step transform in the symbolic FORM

- v10 §1c requires the exact delta/principal-value transform of the constant-plus-Heaviside part (`:225-227`). The saved numerics use a positive regulator.
- The amendment should state whether the A9 symbolic FORM carries the exact transform or a regulator-labelled surrogate.
- Conversion at nonzero transfer (s ≠ 0) is regulator-insensitive. A phase-matched limit (s → 0) is not.

### F8 (inventory): §5a advection probe and §5b residual triplets presented as complete without evidence

- **Source.** v10 §5a requires the one-sided omission/reversal of `u·∇ρ_4D,bg⁰/ρ_4D,bg⁰` and the RHO4 absence operands (`:747-752`, `:763`). §5b requires `K_uniform,end−K_uniform,reference` triplets (`:781-783`).
- **Evidence.** The first-jet report says "the other material routes, applicable advection controls … remain" (`:33-35`). The coordinate report records only "omitted material-covector forcing controls" (`:16`), which is an instrument control. The uniform report records S-matrix reflection and identity checks, not coupling triplets (`:18-24`).
- **Artifact.** A5 (`:109`) and A8 (`:112`) do not list the advection probe as remaining. A6 (`:110`) says "complete".
- **Correction.** Mark these items as unresolved or remaining until a record showing the exact operand is cited.

### F9 (inventory): Option B / A9 understates remaining work

- **Artifact.** A9 (`:113`) and Option B, "Saves B1–B5's unfinished production" (`:166`), omit the amendment's own §3.3 coverage obligation.
- **Why it matters.** Any new real-frequency witness has to recompute all 80 frequency-dependent nonlocal rows, plus end modes and currents, per case (frequency-source report `:7-8`). That is comparable in kind to the deferred complex-row work. The amendment's phrase "one small physical response pilot" (`:178-179`) carries the same understatement (an optional fix there).
- **Correction.** Add the obligation, labelled as proposed, together with this frequency-dependence caveat.

### F10 (inventory): threshold records assigned only to the deferred branch

- **Artifact.** B1 (`:120`) puts the end relations and algebraic threshold candidates only under the deferred pole branch.
- **Why it's wrong.** Amendment §3.3 names these same records as the first step of the retained witness inspection (`:175-177`, `:161-163`).
- **Correction.** Cross-list them as support for A9 and §3.3.

### F11 (inventory, minor)

- **A10.** A10 (`:114`) cites finite-response resolution evidence only (`domain:40-41`). The retained-grade weak coefficients, the priority deliverable, have no stated resolution or regulator comparison. Say so explicitly.
- **Survival.** A4 lists "survival" as computed. The reports show total-current ratios (`response:26-27`), which equal transverse survival only because no other channel is open. State that basis.
- **v10 clearance.** "v10 cleared at `399a8516`" (`:67`): the Opus leg records its nit as "[FOLDED post-clearance.]" (`r10_opus:13`), and Grok's three nit wordings appear in the pinned v10 text. Note that the pinned bytes include post-clearance wording folds.

**Accurately represented:** A3, A4 (except the survival point), A5's numbers, A6's publication, A7, B1's 32/24/8 count (it matches the checkpoint scopes), B3, B4, the nonlinearPoleV2 disclaimer, and the scoped Lean CLEAR status.

## Coverage limitations

- **Not supplied, so not checked:**
  - The build directive and program brief. The amendment's rows for build-directive §§4–5 and brief A–F are unverified, and neither can be searched for a second-frequency clause.
  - The focused completion plan (source of the 1e-4/1e-6 goals).
  - S11b step records (for the "large `k c_s0/|ω|`" clause).
  - The c2 N6 disposition.
  - The Lean fidelity reports themselves.
  - All checkpoint JSON contents beyond the evidence-index metadata.
- **Unverifiable from the packet:**
  - B2's row numbers (60–63, 50–53) and the 12.88 s LAB 2D benchmark.
  - Commits `399a8516` and `7c98b8ee`.
- **Only partially supplied (selected excerpts):**
  - S11b, S11c-a and S11c-c1, and CLAUDE.md lines 117–176. The inventory links review policy to CLAUDE.md:47, which is outside that excerpt.
