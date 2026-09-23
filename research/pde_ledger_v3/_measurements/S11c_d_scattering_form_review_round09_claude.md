# Independent non-author review: S11c-d Option B scattering/FORM amendment (draft 9) and scope inventory (draft 9)

**AMENDMENT: CLEAR**. I found no blocking defect. Seven optional improvements are listed below.
**INVENTORY: CLEAR**. I found no blocking defect. Three optional improvements are listed below.

Both verdicts cover only document fidelity and whether the physics specification is adequate. I did not rerun, recompute or validate any saved result. Neither verdict clears the later build directive, the build, the export, the blind Wolfram engine or T7 work, or any physics result. User adoption of Option B did not affect the verdicts; I reached each conclusion from the supplied sources.

## Method

I read the governing sources first:
- v10 in full;
- the decision list;
- nonlinearPoleV2;
- exploratoryAcceptanceV1;
- the S11b step record in full;
- c1 SHARED_PHYSICS in full;
- the c1 step record in full;
- the c2 excerpts;
- the retained contract excerpt;
- the S11b SHARED_PHYSICS passages at 95–164, 440–484 and 895–916;
- the script excerpts, read without executing them;
- the c1 export-key census.

I then read the saved reports and after that the two artifacts. This order was a method choice, not a blindness control. I make no blindness claim.

## Answers to the five questions

### Q1. Pole deferral, local scattering needs, typed channels, and closure/interface attribution

**The pole deferral is carried consistently through every layer:**
- scope and v10 §2: amendment table rows at 69–70;
- output, export and comparator: rows at 76–79 and the §5 table at 718–719;
- T7: "Pole families are outside the reviewed compared domain; they are not residual zeros…" (815–816);
- downstream consumers: "must receive an explicit unsupported/deferred capability failure, not an empty spectrum" (775–777).

This covers every place v10 carries pole objects:
- v10 §3b emits the pole set and overlap (540–542);
- §4 lists them (687–690);
- §7 has T7 join "bound Riesz data" (884, 893–894);
- §8 lists them as computed (923–924).

**The incorrect v10 formulas are not restored.** Row 77 says: "Do not reinstate v10's incorrect unrestricted projector/residue prescriptions." §4 item 2 keeps nonlinearPoleV2's distinctions (NP §2 lines 56–59 and §3). Examples: a nonlinear inverse residue is not generally a projector, and a zero residue does not exclude a higher-order pole.

**Local scattering needs are kept.** The amendment says the deferral "does **not** remove the resolvent, closed modes, radiation/current construction or local branch/domain checks" (647–650). It keeps the NP §3 full-pairing safeguards for end modes (651–658). It also correctly separates the ∂_ω normalization pairing from the ∂_{k_n} semisimplicity pairing (659–664). That matches v10 §3a (447–454) and NP §3 (79–93).

**The historical search is stated accurately.** The amendment (668–676) matches the contour report:
- "no candidate is resolved by this bounded numerical search on |omega-(1-0.01i)|=0.02";
- "not a certified empty spectrum, a contour-interior holomorphy/exceptional-locus proof, or a physical bound-pole set" (contour report 23–28).

**Confinement statements are correctly limited.** v10 calls P_T,surv "the `N13` confinement object" (536). N13 says bound capture "kills the photon exactly as bulk radiation does" (decisions 130–133). The amendment therefore narrows survival to "continuum-channel transverse flux survival, not a complete N13 confinement answer" (751–753) and keeps confinement open (780–792). That narrowing is faithful and made intentionally.

**End-channel and bulk-depth channels are typed separately.** The amendment's premise is accurate: v10 §3b(i) promises "the thickness **continuum** / bulk escape" (485–486), while §3a builds J_H only from full end channels along y_n (459–466). The amendment keeps them apart and forbids adding them without a derived common control volume (349–360, 721–733).

**Closure/interface attribution is handled honestly.** The amendment keeps:
- the signed exchange and balance operands: "Per-face slab exchange and bulk acoustic power remain distinct signed operands" (362–368);
- the S11b discriminators "as method", restricted to the subcases where they apply (322–324). This matches S11b 463–474, where the check is restricted "only this diagnostic" to the impermeable, propagating, Λ_X⁰=0 case.

A separately normalized absorption observable is openly named as a future scope decision with no owner (735–746). The amendment neither promises nor zeroes it. v10 never required such an observable, so this is a new question being named, not an existing obligation being dropped.

### Q2. A9, the bulk-flux map, grades and the c1 debts

**A9 keeps the full retained problem.** §2 items 1–5 (101–138) keep:
- the full reduced operator, the vertex, and the reference, left and right baselines;
- the derived currents and evanescent matching;
- independent ε, η and σ_W grades, including the mixed grade.

It rules out both failure modes:
- "does not license a simpler Born matrix element without §2's computed reduction premises" (71, per v10 390–393);
- "A frequency-one coefficient table, opaque response symbol, fingerprint, unevaluated missing construction or hand-typed expected scaling is not the general computed FORM" (128–131), consistent with retained-contract items 6–7.

It also keeps v10 §1c's zero-jet distributional transform (140–143).

**The bulk-depth flux is correctly labeled as a new duty, not an inherited object.** The amendment says: "an explicitly added d construction/consume-set and comparison duty, not an existing `IMPORT_KEYS` row" (189–190). That matches the source record:
- The c1 export census has no exterior-field solution-map row. It exports only the DtN, kernel, face-response and resolvent roots.
- c1 §7 makes the energy and far-field diagnostics emit-only.

**The historical c1 energy construction is described correctly.** Against the literal script:
- It uses equal legs, `{q_out_k: qprop, q_out_kp: qprop}` (script 1380–1381).
- It uses an impermeable drive: `closed_coefficients(... (0,0,0))`.
- The Poynting term keeps flat–flat plus the two flat–scattered terms and drops scattered–scattered (script 51–56).

The amendment's warning that this three-term integrand cannot replace a functional carrying B₀[a₁,a₁] (275–284) is correct.

**The power orientations match the carry-forward record.**
- For the acoustic identity, the amendment uses `P_into_bulk − P_outgoing − P_lateral` (248), which is correct for a stationary lossless control volume.
- For the traction test it uses `P_face + P_infinity`, with P_face as slab traction work (256–263). This matches c1 step 174–179, "traction work + positive outgoing far-field Poynting = 0", and the script's `bulk_comparison = -outgoing_flux`.
- It keeps the two checks distinct and restricted to their own subcases.

**The half-space reconstruction is tied to the physical faces and closure.** The amendment:
- keeps the half-spaces disconnected and the outgoing branch fixed (c1 §1b, 99–121);
- uses shifted/tilted faces, outward normals and true-area measures (c1 §2a);
- distinguishes graph height from outward displacement (c1 step 180–184);
- reduces only the in-plane content and keeps depth as a coordinate (164–165);
- requires the actual c1/c2 map from slab state to face drives, "including any elimination of the independent centre displacement… A parity assumption or an invented extra degree of freedom is not a substitute" (199–204).

That last point is physically essential. c1 §1a (85–86) forbids setting ζ_c=0, and the per-face drive is V_s = ∂_t(ζ_c + sδW/2). The trace check (i) correctly requires the pressure to be reconstructed from the bulk field rather than read back from the prescribed datum (242–246).

**Grades and promotion are handled correctly.** There are three statuses: retained-model coefficient, justified leading physical coefficient, and induced diagnostic (291–295). The baseline/interference test is stated correctly both ways:
- Omitted terms can enter through a baseline, per v10 589–593.
- "a computed baseline/interference disposition can establish a leading induced quadratic from first-order amplitudes without a second-order amplitude… quadratic power alone is not a deferral reason" (300–303).

This matters physically. In the uniform exterior fluid the current form has no dependence on the background, and a subsonic incident channel's baseline bulk field is evanescent. So the leading far-field coefficient may well need only first-order amplitudes. The amendment correctly leaves that to computation rather than supplying it (386–388).

**First-shape validity and omitted-flow validity are kept separate.**
- First-shape: "non-grazing asymptotic regime on both momentum legs, with `||N0^-1 N1|| << 1`. A strict rest frame removes no such shape-expansion restriction" (370–374, per c1 step 146–155).
- The bare singularity of N⁻¹ is not transferred to the permeable closure resolvent (374–376, per c1 step 150–154).
- For flux integrals whose support approaches an unsupported region: no deleted endpoints and no whole-total claims (379–388).

**The c1 debt list is complete against c1 step 99–124 and the carry-forward list.** It covers:
- whole-form DtN, flat-leg labeling, ENERGY, t_s and seal 5;
- the four UNMEASURED families;
- the deferred per-family residuals;
- Hermitian/reactive agreement at first order only;
- the HOMOGENEITY keying gap;
- K_a as Hermitian;
- the second-shape caveat on all three grades;
- the independence-scoping limit (206–228);
- the density multiplication operator (230–237, per c1 step 185–187 and c2 §3d.1).

c1 separately corrected the drain-projection wording. The amendment says nothing that conflicts with that correction; it keeps the rest-frame results, per N11a.

**Precision is declared per observable, before comparison** (474–477).

**Current composition is not taken on faith.** The saved construction really is `total_current = slab_current + depth_integral*bulk_current` (mixing excerpt 3378ff), and the consumer hard-requires `left['q'].imag>0 and right['q'].imag>0` (continuum_boundary 40). The amendment states that neither the wiring nor the phrase "closed/nonlocal bulk contribution" proves the general identity (327–334). It also correctly bars carrying decaying-bulk normalization over to a radiating domain (356–358).

### Q3. Open-thickness coverage and validity

**No existing clause requires a second physical real frequency.** v10 §3a works "at fixed `(ω,k_∥)`" (424). Retained-contract item 3 requires only "real continuum frequency and tangential momentum." The "two-frequency endpoint pairing" in the thickness-repair report (327–333) is a pair of frequency legs for the current/energy identity (`FREQUENCY_LEGS`/`BEAT_RATES` in the mixing excerpt). It is not a second scattering point. The amendment and inventory both say so and label A11 as new: "an explicit proposed acceptance clarification, not a claim that v10 already specified an additional numerical frequency" (533–534).

**A11 does real work without demanding a theorem.** The check must use a nonempty outgoing thickness-like space of the full end pencils, real flux-carrying channels, and "a leaky complex-wavenumber mode cannot be silently treated as one of those channels" (576–578). It has a bounded stop rule and keeps structural absence separate from "no witness found" (588–606). This is justified by acceptance item 2 and v10 §6's ban on tautological residuals.

**The source-context reconciliation is correct.** S11b step 158–161 reads: "first order **fails** where `|q v₀ / ω| ≳ 1`; in the `k c_s0 ≫ ω` regime `|q| ≈ k`, so this needs `(|v₀|/c_s0)(k c_s0/|ω|) ≳ 1` — large `k c_s0/|ω|` is **necessary, not sufficient**."

The parenthetical therefore concerns failure of the expansion. It says nothing about which radiating components are valid. Its compressed restatements in N11a (decisions 183–185) and v10 §6 (840–841) can be misread as a universal validity requirement, and the amendment's clarification (402–410) resolves that correctly.

The amendment then uses c1 §2b's explicit conditions: `|q_out·v/ω|≪1`, the boundary-layer bound `|ωv|/(c_s0²|q_out|)≪1` (c1 216–218), subsonic flow, and strict v=0 at grazing. It also correctly notes that no off-grazing estimate covers a radiating spectrum uniformly (393–398).

**Admissibility is kept separate from algebra.** "Algebraic branch candidates… do not settle physical admissibility" (409–410) matches the frequency-source report: "algebraic end candidates with denominator/opposite-sheet artifacts, not classified physical thresholds" (22–23). A12's nonempty-radiating-support check, its stop rule, and "do not use the deferral to waive all of A12" (429–446) hold together. Restricting to decaying bulk cannot discharge A12 (183–187).

### Q4. Resolution targets and stop conditions

These are adequate for the limited toy-model claims:
- 1% stability, amplitude 1e-4 and current 1e-6 in the declared frame only, with precision declared per observable before comparison (469–477);
- a small number of independent resolution, domain and regulator changes (479–488);
- near-zero effects stay unresolved, and there is no denominator floor (492–495);
- no assumed unitarity (496–497);
- certification is not a gate but is still required for theorem-strength claims (499–504).

This is consistent with acceptance 42–55. The current meanings are pinned in the export contract rather than in prose (721–733, 748–758).

### Q5. Saved results and completeness

The saved results are represented accurately. I spot-checked every inventory number against the supplied reports and found no discrepancy:
- 645 unknowns and four incident columns;
- 687 s, 278 s, 432 s, 53 s, about 71 min and 170 s;
- 8.82e-3 vs 3.22e-5; 8.53e-3 vs 3.36e-5; 2.23e-5;
- 8.243e-7 and 5.482e-7; 2.22e-13 and 2.97e-13;
- `q_depth^2=-k_normal^2-1/25`;
- LAB condition 5655 and residual 9.38e-16; row46 sensitivity 4.67e-14;
- 24 of 32 owners, with material 2D rows 60–63 and 50–53;
- Wolfram times 3,196.55 s and 601.61 s;
- `S11c_d_exports.py` absent.

Nothing is silently declared complete. A3/A4 are marked computed but only V/P. A6–A8 locate-first gaps are explicit. A9, A11 and A12 are unfinished. C1–C4 are pending.

Every v10 §§3–5 emit family is either retained, deferred with metadata, or named as outside scope. That includes the strong-edge handoff (459–463) and the §5d FORM (99). I found no contradiction between the two documents.

## Findings

None are blocking.

**Amendment**

1. **(Optional) Name the §1c reduction rules for new c1 operands.** The half-space map will consume c1 in-plane operands, such as the kernel and response resolvents, which carry 3-D momentum deltas. v10 §1a/§1c explicitly exclude those from the reduction-convention operand. The only rule the amendment states is "The new map cannot bypass the in-plane Fourier/measure rules" (165). Suggested fix: name v10 §1c's requirements for these new roots and their T7 join — a computed, both-operand reduction with no typed (2π) map, and stripping of δ²(Q_∥) before flux.

2. **(Optional) State the both-legs rule for the half-space map.** The ban on a single-k or left-quantized freeze (c1 §3a lines 247–259; §5e) now applies only through the citation of "c1 §§1d, 2a and 3" (159–160). A one-line explicit statement would help.

3. **(Optional) Consider a fresh key for the narrowed continuum-conversion root.** The amendment narrows the v10 tag `S11CD_CONTINUUM_CONVERSION`, which v10 §3b(i) defines to include bulk escape, to "thickness end-channel conversion only" (721–726). No d export exists yet, so no consumer is affected now. A fresh key or versioned schema would still guard against a consumer reading the v10 meaning.

4. **(Optional) Say whether d executes the traction-sign test.** The amendment "retain[s]" the corrected `P_face+P_infinity` test (258). It does not say whether d runs it at its own faces under a diagnostic-only impermeable specialization, as S11b 463 prescribes ("Restrict only this diagnostic"), or only carries c1's result with ENERGY still UNDECIDED.

5. **(Optional) Tighten citation targets.** The target of each change is still identifiable from its text, but several anchors are off:
   - [rest-source] points to `:155`; the passage is at 158–164.
   - [rest-decision] points to `:178`; N11a is at 179–186.
   - "§3b(ii), including its survival definition": the survival paragraph (v10 527–537) is not part of item (ii).
   - [spec-controls] points to §5a.
   - [spec-outputs] covers only §4.

**Inventory**

6. **(Optional) Cite an existing lead for the post-repair normalization.** A1/A2/A4 correctly call post-repair endpoint pairing and current-adjoint normalization "unestablished" in the repair report. However, a supplied consumer already requires `S11c_d_end_normalization_<end>_thickness_repair_checkpoint.json`, with `not cp['unaccountedResidualNormsAboveDiagnosticThreshold']` and an accepted pairing packet (continuum_boundary 65–73). Citing it would give a concrete starting point, without treating it as completion.

7. **(Optional) Mention the tracking of profile moments and reconstruction.** v10 §4 requires profile moments/form factors and the 3-D↔1-D reconstruction of every c2 carrier. The inventory tracks these only generically under A2 ("Computed reduction"), and the reconstruction appears only in C3.

8. **(Optional) Update a stale sentence.** "The first independent document-review round requested revision" (inventory 105–107) is historically true but stale after eight rounds.

9. **(Optional, both documents) Reconcile the S11c-e wording with c1.** c1 already assigns first-order-truncated O(η²) leakage to "S11c-e" (c1 §3b 338; c1 step 106–107). The amendment's refusal to claim an approved e owner or plan (48–51) is appropriately cautious, and it could cite that existing c1 assignment.

## Coverage limitations

- **Excerpts only:** c2 spec and record, the builder report, all scripts, and CLAUDE.md, which I inspected by grep.
- **Read in part:** S11b SHARED_PHYSICS (the ranges listed above); S11c-a (not read beyond references); the evidence JSON (grep only); the material-row summary (row indices only); the thickness-repair report (1–100 and 280–400).
- **Verdict line only:** the r10 Grok report.
- **Not read:** POLE_HANDOFF, and the inertia, mechanical and c2 trace repair reports.
- **Not supplied:** the round 1–8 drafts and reports. I cannot verify the recorded verdicts or the post-round-5 and post-round-8 user approvals. No arrays or checkpoints were available, so no numerical claim was independently validated.
- **Not assessed:** cost bands. They are planning judgments.