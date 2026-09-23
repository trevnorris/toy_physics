# S11c-d Option B amendment and scope inventory: independent document review

**Reviewer:** fresh non-author Claude agent (document branch). I did not run, rebuild or recompute anything. This review covers document fidelity and whether the physics specification is adequate. It does not clear any build, export or reported number.

## Verdicts

- **AMENDMENT: NEEDS REVISION.** Two blocking findings (A1, A2). Both concern the newly added bulk-depth half-space construction, which does not carry forward limits that c1 itself recorded. The pole deferral, the Born-shortcut prohibition, the reading of the N11a/S11b wording, the §3.3 open-channel coverage rule, and the attribution disposition are sound.
- **INVENTORY: NEEDS REVISION.** One blocking finding: its review-status chronology contradicts itself. One dependent fix: A12 must take on the A1/A2 limits. Its account of the saved results matched every report I checked.

## How I read the sources

I read the sources before the artifacts. From them I took the following.

**v10 (d spec)**
- It requires a complete two-ended S-matrix built on the §1c-reduced operator, with end channels taken from the full end pencils (§2, §3a).
- `J_H` is defined only from end-channel currents (lines 459–466).
- §3b(i) says continuum conversion covers "the thickness **continuum** / bulk escape", but defines it through "the §3a amplitude projected on radiating/continuum channels" (lines 485–486). That wording does conflate the two, as the amendment says.

**N11a and the S11b standing limit**
- N11a (`S11c_decisions.md:183-185`) compresses the S11b record into "large `k c_s0/|ω|` is necessary, ⛔ not sufficient".
- The full S11b passage (`steps/S11b_interface_coupling_law.md:158-161`) is about **failure**: "first order fails where `|q v₀ / ω| ≳ 1`; in the `k c_s0 ≫ ω` regime `|q| ≈ k`, so this needs `(|v₀|/c_s0)(k c_s0/|ω|) ≳ 1` — large `k c_s0/|ω|` is **necessary, not sufficient**".
- c1 §2b (`S11c_c1_SHARED_PHYSICS.md:209-221`) gives the operative conditions: two smallness conditions plus a subsonic condition, with grazing treated as the strict `v_bulk_normal_0=0` result.

**c1**
- c1 exports DtN/impedance operators, kernels, face response and resolvents. The mechanical census has 44 ledger keys, including `dtn_kernel`, `s11cc1_flat_normal_dtn_inverse` and `s11cc1_q_out_{input,output}`. There is **no** exported half-space field or far-field flux row. The energy objects are emit-only (c1 §7, lines 571–573).
- c1's energy construction is limited to one slice: real ω, propagating regime, impermeable faces, Λ_X⁰=0 (c1 §3b item 3, lines 321–330). The script excerpt also sets both legs equal: `xreplace({q_out_k: qprop, q_out_kp: qprop})` (excerpt lines 77–79).
- The c1 step record adds limits I rely on below: grazing non-analyticity, second-shape evanescent leakage, energy-residual orientation, and cross-engine UNDECIDED items.

**Saved results at ω=1**
- Four open transverse directions and zero open thickness directions.
- `q_depth^2 = -k_normal^2 - 1/25`, so bulk-depth radiation is closed for every real k_n.
- Bulk branch frequencies are `±sqrt(5)` (frequency-source report lines 20–23).
- The contour check on `|ω−(1−0.01i)|=0.02` resolved no candidate and is explicitly not a certified empty spectrum.

## Answers to the five questions

### 1. Pole deferral, retained scattering and attribution — mostly yes

The precedence table covers what the deferral needs to cover:
- v10 §§0, 2, 3b, 4, 7 and 8;
- nonlinearPoleV2 §§1–7;
- the acceptance addendum;
- the retained contract;
- the N2 handoff (amendment lines 54–67).

§5 is explicit on outputs:
- The two pole families are removed from output and export membership: "Do not export zero, an empty set, an empty result-bearing container or a generic pole function under these names" (line 492).
- They are recorded as outside T7 coverage, "not residual zeros, evidence of equality or a fourth T7 truth value" (lines 556–557).
- Consumers get "an explicit unsupported/deferred capability failure" (line 533).

The deferral does not drop anything local that scattering needs:
- "This deferral does **not** remove the resolvent, closed modes, radiation/current construction or local branch/domain checks" (lines 422–425).
- The distinction between the ∂_ω and ∂_{k_n} pairings is correct (lines 434–439).

The historical claims are limited correctly:
- The search summary (lines 443–451) matches the contour report exactly, including "not a certified empty spectrum" and "contour-interior holomorphy/exceptional-locus proof".
- Survival is restricted to "continuum-channel transverse flux survival, not a complete N13 confinement answer", and "a near-unity value there does not establish confinement" (lines 519–526).

End-channel conversion and bulk-depth escape are typed as two separate observables (lines 495–502). The closure/interface disposition keeps the signed exchange terms and the checks. It neither promises nor eliminates a separate absorption observable: "not a computed positive absorption rate and not a deferral of required signed exchange/balance checks … If source inspection identifies an existing mandatory operand or check, it remains current work" (lines 505–512). None of v10, N13 or S11b mandates a separately normalized d-level absorption observable, so nothing is silently eliminated.

Two points remain:
- **Blocking (A1):** c1 assigned the O(η²) evanescent-channel leakage to S11c-e. The new bulk-depth escape root claims it for d without reconciling that assignment.
- **Optional (A4):** the N10 owner requirement for the family card.

### 2. A9 and the bulk-depth map — A9 is sound; the bulk map misses c1 limits

A9 keeps the full retained problem:
- "Preserve the full reduced local/nonlocal operator and coupling vertex, computed reference/left/right baselines, both-end channel spaces and derived signed current forms" (lines 87–90).
- Independent grades: "Only the declared physical homotopy relates eta and sigma" (lines 92–94).
- Baselines, interference and current forms (lines 101–108).

It rejects both failure modes named in the prompt:
- No Born shortcut: "does not license a simpler Born matrix element without §2's computed reduction premises" (line 58).
- No tables or placeholders: "A frequency-one coefficient table, opaque response symbol, fingerprint, unevaluated missing construction or hand-typed expected scaling is not the general computed FORM. Do not demand elementary closed form where a permitted transparent operator/integral representation is appropriate" (lines 110–115).

The bulk map is correctly presented as a new duty rather than an inherited object:
- "This is an explicitly added d construction/consume-set and comparison duty, not an existing `IMPORT_KEYS` row" (lines 148–149).
- "The existing c1 far-field energy construction has a specified restricted subcase" (lines 152–154). The script excerpt confirms this: equal legs, impermeable, Λ_X⁰=0.

The in-plane and depth representations are handled correctly: "Only the in-plane content undergoes the §1c reduction; the exterior depth coordinate remains" (lines 145–146).

However, the bulk map's validity, truncation, orientation and inherited status do not carry c1's recorded corrections. The amendment cites the c1 carry-forward list only for "graph height" and "impedance from the normal DtN operator" (lines 142–144). See A1 and A2.

### 3. Conditional open-thickness-channel coverage — adequate

No existing clause requires a second numerical real frequency. I searched the whole packet:
- v10 §3a requires only "every incoming open mode" and "every outgoing open mode" at the chosen `(ω,k_∥)` (lines 427–431).
- The retained contract (item 3) says: "Test actual reference/end channel availability before attempting flux normalization. Do not manufacture an incident channel".
- The only "two-frequency" items in the packet are current-normalization pairing diagnostics. Examples: `S11c_c2_trace_repair_report.md:106` "physical current normalization still requires the subsequent two-frequency/end-mode construction", and `S11c_thickness_coordinate_repair_report.md:327-328`. These concern ∂_ω normalization, not open-channel coverage.

The amendment correctly labels its requirement as new: "This is an explicit proposed acceptance clarification, not a claim that v10 already specified an additional numerical frequency" (lines 323–324), and it repeats this at line 395.

The coverage rule is non-vacuous without demanding a global theorem:
- It requires a nonempty thickness-like end-channel space of `𝓛_±^full` with a nonzero incident denominator (lines 326–334).
- It requires diagnosing the closure mechanism (lines 336–345).
- It stops on structural absence or on no witness, and keeps those two cases distinct (lines 385–393).
- It excludes leaky complex-k_n modes (lines 363–365).

End channels and depth radiation are kept apart: "two distinct availability questions … Neither is inferred solely from the other" (lines 313–316); "It cannot discharge the `J_H` check, and an open `J_H` channel cannot discharge its bulk-flux check" (lines 404–406).

The N11a reconciliation (lines 203–211) matches the full S11b passage quoted above. The parenthetical is a condition for the flow correction to **fail**, not a lower-wavenumber requirement for valid radiation. The amendment then correctly defers to the explicit c1 §2b conditions (lines 191–201). This matters physically. At the approved k_∥, depth radiation needs ω > √5, which means `k c_s0/ω < 1`. Reading N11a as "large `k c_s0/ω` required for validity" would have made the bulk-escape obligation vacuous by misreading.

The amendment also does not let algebraic branch candidates stand in for a physical radiating domain: "Algebraic branch candidates, including those in the saved frequency-source report, do not settle physical admissibility" (lines 210–211).

### 4. Tolerances, current meanings and stop rules — adequate for the retained claims

- The 1% target and absolute 1e-4 amplitude / 1e-6 current goals are confined to their channel frame (lines 263–267).
- Near-zero and tiny effects are handled: "a denominator floor does not establish their relative accuracy, sign or absence" (lines 283–285).
- Unitarity is not imposed: "do not impose lossless S-matrix unitarity by assumption" (lines 286–287).
- Current meanings are typed separately (lines 495–502): `J_H` is end-channel only, depth-integrated bulk current is part of the end current (lines 180–184), and survival is restricted.
- The stop and scope-decision rules for A11 and A12 are bounded (lines 230–240, 375–393).

The one adequacy gap is the grade/truncation status of the bulk operands at the λ² slot (A1b).

### 5. Saved results and review state — the amendment is accurate; the inventory's data are accurate but its status text contradicts itself

I checked each saved number the inventory uses against its report:
- **Response report:** 645 unknowns, regulator 0.1, 687 s.
- **Flux report:** four open transverse / zero thickness directions, `q_depth²`, 277.79 s.
- **Uniform report:** twelve controls at `44de056f`, no coupling-triplet evidence.
- **Profile report:** sub-resolution changes, a bump moment only.
- **First-jet report:** 8.243e-7 / 5.482e-7, 4247 s.
- **Coordinate report:** 2.22e-13 / 2.97e-13.
- **Domain report:** 170 s.
- **B-branch reports:** 9 + 15 = 24 of 32 owners done, rows 60–63 and 50–53 remaining (matches `material-row-summary.json`), LAB condition 5655 and residual 9.38e-16, row46 spread 4.67e-14, contour 639.84 s.
- **Thickness-coordinate report:** the endpoint-pairing diagnostic is correctly shown as unestablished (the report ends at line 399 with "no fresh pairing result is claimed").

The inventory declares nothing complete without evidence. The contradictions are in its review-status narrative (I1).

## Findings

### A1. BLOCKING — the bulk-depth FORM and check (ii) omit c1's recorded grazing and truncation limits and its S11c-e assignment

**Sources**

- c1 step record, lines 146–155: "The **DtN inverse `N⁻¹`** (hence `Z=iρ_mω·N⁻¹`) is **nonanalytic as `q_out→0`** … at exact double grazing … carries a `~1/η` Laurent pole … so the first-shape coefficient `Z₁` is a valid **non-grazing asymptotic** coefficient ONLY — `‖N₀⁻¹N₁‖≪1` imposed on both legs … the exact-grazing *threshold response* is ⛔ **not claimed** at c1's first shape order".
- c1 spec §3b, lines 332–338: "A flat evanescent channel is a nullspace of the zeroth-order Hermitian form; the curvature-induced radiated power there is `O(η²)`, exactly the term the first-shape-order truncation omits … on its nullspace emit the typed token `NOT_ESTABLISHED_AT_FIRST_SHAPE_ORDER` (the `O(η²)` leakage there belongs to S11c-e)".
- c1 record, lines 191–193: the completion is "**second-shape order = η², η·σ_W, σ_W²**".
- c1 §2c, line 232: "`|∇_x h_s|²` is `O(σ_W²)`". This is the true-area correction.

**Artifact**

- Lines 158–164: "(ii) the signed true-area acoustic face power versus independently evaluated outgoing bulk control-surface/far-field flux. Match the tangential normalization, control-volume boundaries, retained grades and limit premises".
- The only grazing language (lines 194–197, 366–368) concerns the `v_bulk_normal_0` rest-frame limit of c1 §2b.
- Line 60 claims bulk escape for A9 with no mention of c1's S11c-e assignment.

**Why this is a physics defect**

- **Grazing.** The radiating support is `|k_n| < sqrt(ω²/c_s0² − k_∥²)`. Its endpoints are exactly where the output leg's `q_out → 0`, and there the first-shape Z₁ and N₀⁻¹ are invalid. This is separate from the flow-validity issue. Whether a λ² bulk-escape coefficient is uniform up to those endpoints is a claim that must be established; nothing currently requires it.
- **Truncation.** Suppose the computed baseline puts nonzero, bulk-evanescent face data on the faces (v10 insists `K₀` is computed, not assumed zero). Then:
  - the face-power operand at the λ² slot needs omitted second-shape terms: Z₂, and the O(σ_W²) true-area weight;
  - the far-field operand can be formed from first-order amplitudes;
  - so (ii) would give a residual equal to the omitted truncation, not a construction error.

  As written, that residual is "unexplained" (lines 168–169), and nothing types it as `NOT_ESTABLISHED_AT_FIRST_SHAPE_ORDER`. Likewise, when the baseline itself radiates, the bulk FORM's λ² coefficient contains a `2Re(a₀a₂*)` term needing omitted second-order amplitudes. v10 §3c states this truncation-status rule for `J_H` (lines 589–592). The amendment does not extend it to the bulk root.
- **Scope.** c1 assigned exactly this O(η²) curvature-induced leakage to S11c-e. The amendment moves it into d without saying so.

**Minimal correction.** Add a short paragraph to §2's bulk-escape block that:
1. Carries c1's non-grazing validity of Z₁ and N₀⁻¹ (both legs) into the bulk FORM's regularity domain, separately from c1 §2b. Grazing endpoints of the radiating support are excluded or labeled NOT_ESTABLISHED unless their contribution is shown to be of higher order.
2. Requires each λ-slot of check (ii) and of the bulk FORM to list its grade content. Where a slot needs omitted second-shape closure or measure terms, or second-order amplitudes via baseline interference, the slot emits the c1 token or v10 §3c truncation status instead of a residual verdict.
3. States explicitly, as an N2 boundary refinement, that the c1-assigned S11c-e O(η²) evanescent leakage is now owned by d/A12, or leaves it at S11c-e.

"Retained grades and limit premises" does not cover this. It never says the check can be structurally unable to close at the λ² slot.

### A2. BLOCKING — c1 inherited status and orientation corrections are not propagated to the new map

**Sources**

- c1 record, lines 108–111: "the **`dtn_operator` whole-form, the t_s traction, and the seal-5 density are UNDECIDED — c2 must NOT treat them as cross-engine-closed**".
- c1 record, lines 113–116, add: "the off-diagonal **flat-resolvent leg-labeling** (PY output-leg `q_out` vs WL input-leg …)", and the "**ENERGY** audit (PY closed-form … vs WL unevaluated far-field flux integral)" as UNDECIDED.
- v10 §1b, lines 166–167, carries these only as c2-level premises.
- c1 record, lines 174–179, energy-residual orientation: "The independent energy identity is `P_face + P_∞ = 0` … ⛔ NOT the literal `A−B` the spec wrote — a literal `A−B` … vanishes only for the WRONG `t` sign … define the bulk subtraction operand `B` as **minus** outgoing Poynting … or emit `A+B`. Carry to c2's energy control."

**Artifact**

- Line 559: "The existing c2 operand debt remains propagated". Nothing is said about c1 items that the new map consumes **directly**, bypassing c2.
- Line 161: "signed true-area acoustic face power versus … outgoing bulk … flux", with no orientation given.
- Lines 142–144 cite the c1 corrections only for graph height and Z/N terminology.

**Why this is a physics defect**

- The new half-space map is a direct consumer of c1 operands. The leg-labeling question decides which momentum's `q_out` carries the reconstructed field to depth. That is exactly the radiating-support and flux quantity A12 claims.
- Reusing c1's ENERGY operands, which the amendment explicitly allows, inherits an UNDECIDED cross-engine status.
- An unspecified "signed … versus" residual repeats the defect c1 recorded: a literal A−B that vanishes for the wrong sign.

**Minimal correction.**
1. List the c1 cross-engine UNDECIDED items (dtn_operator whole-form, flat-resolvent leg labeling, ENERGY audit, t_s leaf, seal-5 density) as named premises that the bulk-depth root and its checks propagate.
2. State the orientation of check (ii): power delivered outward into each half-space against positive outgoing flux, or equivalently the `P_face + P_∞` form, citing c1 record lines 174–179.

Line 559 covers c2 debt only. Lines 142–144 omit these two corrections.

### A3. Non-blocking (recommended) — per-face drive and independence of the trace check

- **Sources:** c1 §1a, lines 85–86: "⛔ `ζ_c` is an independent face DOF; ⛔ no S11c-c1 computation may … replace the two face variables by a thickness-only ansatz". c1 §2a, lines 191–195: parity and `{δW,ζ_c}` coupling are "a **computed** result … ⛔ not asserted here, in either direction". v10 §1a, line 117: the d field carriers are `θ, e_W, u{1,2,3}`, with no ζ_c.
- **Artifact:** lines 138–139 say "driven by the same per-face kinematic/closure data as the solved slab state".
- **Problem:** the solved d state has no centre-shift variable. Rebuilding two independent per-face drives, `V_s` and `J_s`, needs either the c2 elimination map or a cited computed parity-decoupling result. The w-reflection-symmetric background makes decoupling plausible, but c1 forbids asserting it.
- **Suggested fix:** require the build plan to state how per-face data are recovered.
- **Also:** say which operand in check (i) is independent, for example the reconstructed field evaluated through the shifted-trace map against the closure's physical-face pressure. This avoids an A−A velocity trace. v10 §6's "No tautological residual" already binds, so this is not blocking.

### A4. Non-blocking — the N10 owner requirement for the family card

- **Source:** N10 (`S11c_decisions.md:175`): "⇒ the family's card names each deferral's owner".
- **Artifact:** lines 43–44: "The family roll-up must record this work package as unassigned until an owner is actually appointed".
- **Suggested fix:** state that the S11c roll-up card cannot be finalized while the pole package, or the closure-attribution decision, has no recorded owner or no explicit user disposition.

### A5. Nits (optional)

- Table row at line 59 cites "v10 §3b(ii), including its survival definition". Survival is defined after (ii), at v10 lines 527–537.
- The control-volume paragraph (lines 162–164) could state the order of limits: lateral extent to infinity before depth, or the equivalent. Otherwise radiated energy can cross the lateral faces and be confused with the end-channel bulk tails that lines 180–184 correctly put in the end current.

### I1. INVENTORY BLOCKING — the review-status chronology contradicts itself

- **The header says four rounds and a fourth fold:** lines 14–19: "Draft 4 and round-4 reports are preserved at `116f6b36`: Grok cleared both documents … this fourth substantive fold".
- **Lines 200–203 say three rounds and a third fold:** "three document-review rounds have completed … This third substantive Codex fold uses fresh Claude/Grok non-author reviews".
- **The C2 row (line 150) is stale:** "Round 3 gave one amendment CLEAR and one NEEDS REVISION … Draft 4 remains P". Draft 4 has been reviewed, and the current draft is draft 5 (amendment line 3).
- **The preservation list omits draft 4:** lines 222–224 list `2e8ec5a3`, `75086f8a` and `134d033d`, but not `116f6b36`.

Review status is one of the inventory's stated functions, so these must be reconciled. **Minimal correction:** update lines 200–205, the C2 row and lines 222–224 to match the header.

### I2. INVENTORY, dependent on A1/A2

The A12 row (line 131) correctly marks the half-space map, trace and power checks, and support as unfinished. When the amendment is fixed, A12 should also carry:
- c1's non-grazing Z₁ validity;
- the second-shape `NOT_ESTABLISHED` / truncation status at the λ² slot;
- the S11c-e reassignment;
- the directly consumed c1 UNDECIDED items.

### I3. INVENTORY, optional

v10 §4 (lines 697–698) requires "reconstruction of every c2 three-dimensional carrier" and derived profile moments/form factors. A2 says "Computed reduction … exist", and C3 lists the T7 join. A2 could state the actual status of each engine's both-operand reduction and 3-D↔1-D round-trip record, or mark it "locate first", as A6–A8 already do for their own operands.

## Coverage limitations

- **Excerpts only:** the retained contract (lines 488–532), the c1 SymPy script (765–831 and 1380–1466) and CLAUDE.md (47–65 and 117–180) were supplied as excerpts. I did not see the rest of those sources.
- **Not supplied:** the c2 spec and step record. So whether and how ζ_c and the per-face data were eliminated in `s11cc2ClosedSlabOperator` is unverified (A3). I also cannot verify how the c1 flat-resolvent leg-labeling reaches c2 or d.
- **Partly read:**
  - S11b spec: lines 1–330, 440–500 and 680–940 were read; the §6 geometry section (330–440) and B0–B3 (500–680) were checked by grep only.
  - S11c-a: read to line 120.
  - Grok r10 review, inertia and mechanical repair reports, and POLE_HANDOFF beyond line 60: skimmed or grepped only.
- **Hashes:** I did not independently verify pinned hashes, including the v10 SHA and the suffix SHA.
- **Nothing recomputed:** all numbers are as reported. The export-key index is a declaration census and establishes nothing about nested payloads.
- **Process note, not a finding:** the amendment records a G4 author-change exception on the fourth fold under a user preference. I made no finding on it.