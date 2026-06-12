# Step 6 — Close-out punch list (consolidated, living)

Single canonical checklist for the layer-2 adversarial-audit close-out. Aggregates every non-fatal item
surfaced across Phase C batches + the Phase A/B tooling deferrals. **Append each new batch's items here**
(provenance stays in the per-batch verdicts files cited per item). Do NOT run Step 6 from the scattered
sources — run it from this list.

**Owner convention:** Codex APPLIES all content-bearing fixes (cards/notes/scripts/provenance YAML); Claude
REVIEWS + owns this tracker. Items tagged `[USER-GATED ✅]` were already approved by the user. There are now TWO
`verdict_logged` PARTIAL findings (both non-fatal, manifest verdict recorded): **A1** (stage_192 χ_Q/Δ_Q/N_Q
re-attribution, batch 1) and **A8** (stage_162 γ₀ exact-family-formula scoping, batch 5 — Claude+Codex adjudicated
non-fatal, contingent on the ANSATZ_LEDGER §5-1 user gate).

**Status legend:** `[ ]` open · `[x]` done. Nothing here affects verification substance — all non-fatal
(cosmetic / metadata / renumber / card-text scoping). No program revision required.

---

## A. Card-text overclaim & claim-status scoping (the stage_192 class — fix CONSISTENTLY)
These all concern a card/badge implying more than the notes support. The blanket-`StatusExactClosure` items
(A2/A3) are the SAME class — scope them with one consistent rule so the cleanup doesn't diverge.

- [ ] **A1.** `paper/stages/stage_192.tex:13/20/21` (+ the shared 192/193 boilerplate block) attributes
  χ_Q/Δ_Q/N_Q to stage 192, but those identities are genuinely downstream (χ_Q=1 at 194; Δ_Q:=χ_Q−1, N_Q=χ_Q⁻¹
  at 195). Re-attribute to 194/195. _[batch 1, Item 1 — PARTIAL, verdict_logged]_
- [ ] **A2.** `paper/stages/stage_189.tex` bare `\StatusExactClosure` without card-level GR-target disclosure. _[batch 1]_
- [ ] **A3.** `paper/stages/stage_116.tex:5` section-level blanket `\StatusExactClosure` sections over the posited
  pure-scale realization (γ₀ "A simple concrete realization is to take…"). Caveat-mitigated at `:27` but scope the
  badge. SAME class as A2 — apply one consistent scoping rule. _[batch 3, u05]_
- [ ] **A4.** `paper/stages/stage_001.tex:149` omits the "ansatz"/"not yet frozen" hedge the notes carry
  (under-disclosure asymmetry; bundles self-flag severity:low). _[batch 2, u01]_
- [ ] **A5.** `paper/stages/stage_073.tex` boxes L/a=37/20 under `\subsection*{Derivation}` (37/20 is the INPUT to
  the derivation, not derived) — mild presentational tension, self-flagged in `fit_stage073_paper_fixes_geometry_ratios`. _[batch 2, u11]_
- [ ] **A6.** `paper/stages/stage_073.tex` boxes ε_r = 1/20 without an inline "posited" tag (surrounding
  "selected/declared/reference" framing carries it). _[batch 2, u12]_
- [ ] **A7.** `paper/stages/stage_115.tex:1` title "Exact Core-Balance Compensation **Theorem**" is borderline-strong
  wording (it's about the balance surface/law, not a value; notes re-disclose free data). LOW priority. _[batch 3, u08]_
- [ ] **A8.** `paper/stages/stage_162.tex:5` `\StatusExactClosure` + `notes162:38/58/70` present γ₀=(1+𝔯²)/9 as one of
  "the two **exact** parent-family formulas" / "Stages 115–116 **fix** the bare odd normalization" with no ansatz tag —
  scope the badge + add a pure-scale-realization hedge distinguishing γ₀ (posited) from the genuinely-derived L_W/a
  relation it is grouped with. **SAME class as A3** (stage_116 γ₀ badge). **CONTINGENT** on the user's
  ANSATZ_LEDGER §5-item-1 call: if the user adopts the "derived (odd-preservation theorem)" reading, NO edit is owed
  and the current wording is correct. _[batch 5, u02 — YES→PARTIAL, verdict_logged; Claude+Codex adjudicated non-fatal]_

## B. constraint_kind / classification relabels (provenance YAML, label-only)
- [ ] **B1.** Relabel `constraint_kind` `published_target → free_choice` in the 3 stale r_F1 bundles:
  `fit_stage156_carried_constants_block__r_f1.yaml:191`, `fit_stage158_canonical_quartet_literals__r_f1.yaml:193`,
  `fit_stage161_rf1_family1_radius__r_f1.yaml:43`. (162/163 already free_choice.) **[USER-GATED ✅ — 37/20 = free_choice]** _[batch 1, Item 2]_
- [ ] **B2.** `fam_0212` V_known bundle `published_target → free_choice` (illustrative, not external). _[batch 1]_

## C. Stale / misdirected stage-anchors (the +17 EM-renumber drift family + misattributions)
Content-harmless label drift; see [[project-numbering-drift-root-cause]] — content-keyed, NEVER offset-sweep.
- [ ] **C1.** Pe_fail_chi (089) false "082" upstream anchor (script 61-62/92 + `CHECKPOINT_CONSTANT_PROVENANCE.md:111`);
  089-notes "086" anchor is the correct rho-window anchor. _[batch 1]_
- [ ] **C2.** `fam_0351` stale "Stage-18" anchor → Stage 035 (`CHECKPOINT_CONSTANT_PROVENANCE.md:376`). _[batch 1]_
- [ ] **C3.** g_UV "carried from 245" stale symbol-anchor. _[batch 1]_
- [ ] **C4.** `scripts/…stage090…sympy_audit.py:15` says c_contact "fixed upstream" but the 3/4 forcing is DOWNSTREAM
  at 091 (`notes091:34`) — misdirected forward-reference. _[batch 3, u02]_
- [ ] **C5.** Notes prose at 091/092/093 cites "Stage 71/74/75/76/77" (+17 drift → canonical 091-094). _[batch 3, u03]_
- [ ] **C6.** stage-109 notes cite "Stage-90" (:21) / "Stage 91" (:79); the χ_Q deformation formula originates at
  Stage 107/108 (+17/+18 drift). _[batch 3, u07]_
- [ ] **C7.** `paper/stages/stage_118.tex:7` reads "Stage~135 is a core outlet realization ledger step" inside a
  stage-118 card (cross-stage label drift). _[batch 3, u09]_
- [ ] **C8.** stage 070 establishes no own notes-level provenance for n=5 (inherits 062/Phase-0) — stale-anchor metadata. _[batch 2, u08]_
- [ ] **C9.** Ω_A_l bundle `introduced_at_stage:17` — benign EM-extension renumber artifact. _[batch 2, u05]_
- [ ] **C10.** Provenance `downstream_dependents` for 075→'062','063' and 076→'063','078' are stale/pre-renumber
  (cards route to 076/077/078). _[batch 2, u14]_
- [ ] **C11.** `n_eos introduced_at_stage:0` is a phase-0-spec encoding artifact (not a per-stage file). _[batch 3, u09]_
- [ ] **C12.** `paper/stages/stage_121.tex:7` Purpose field reads "Stage~138 is a core outlet realization ledger step"
  inside a Stage-121 card (wrong stage number in the boilerplate Purpose). _[batch 4, u01]_
- [ ] **C13.** `paper/stages/stage_139.tex:7` Purpose has the right number but borrows Stage-133's role description
  ("coupled mouth fixed point") — 139 is "Actual Family-1 Mouth Gains." Low priority. _[batch 4, u01]_
- [ ] **C14.** `paper/stages/stage_121.tex` boxes/uses L/a=37/20 without citing its Stage-073 origin freeze (card-level
  provenance gap; the value is honestly "carried" in notes but the card omits the upstream anchor). _[batch 4, u01]_
- [ ] **C15.** stage139 L_over_a bundle has no notes_stage line; scanner anchored only a Stage-140 script comment (origin
  recoverable from 073/121) — metadata-only. _[batch 4, u01]_

## D. Dedup / scanner-artifact metadata
- [ ] **D1.** `fit_stage_065_v_0f` is a tokenization artifact (scanner fused `V_0 f(…)` from `stage_065.tex:7`) —
  mark as alias/dedup of V_0, not a standalone free_choice candidate. _[batch 2, u09]_
- [ ] **D2.** The same 3/4 fraction is keyed c0 (088) / c_contact (090) / K_geom-vs-K_pole split (091) — add a
  cross-reference note (benign naming variation). _[batch 3, u02]_
- [ ] **D3.** Several scanner-generated candidate slugs overstate vs the actual card/notes content — `residual_freedom_fixed`
  (133, M_minus is OPEN), `bias_determined_by_fixedpoint` (134, M_q is FREE), `rf1_geometry_fixes`,
  `g_nat_equal_normalized_unity`, `parent_threshold_canonical`. These are audit-internal candidate_ids, NOT published
  paper/notes lines → **no paper fix required**; informational only (the slug-vs-content gap was batch 4's headline). _[batch 4]_
- [ ] **D4.** `fit_stage195_mhat0_natural_source_map` and `fit_stage195_natural_source_map_mhat0` are word-order-permuted
  slugs for the SAME parameter mhat_0 at stage 195 (the second bundle self-identifies "Duplicate candidate for m_hat_0
  at stage 195") — mark as alias/dedup, not two standalone free_choice candidates. _[batch 5, u07]_
- [ ] **D5.** Two more scanner-slug overstatements (audit-internal candidate_ids, NOT paper/notes lines → **no fix
  required**): `g2_common_tangent` (187 — no literal in notes; denotes the tangent space T_id G_* of which Lambda_1 is
  one free direction) and `doubled_slope_factors` (178 — the factor-2 lives in the derived identity ν_r=2(𝔭_r−𝔡_r),
  not in κ_1/d_r being "doubled"). Informational only. _[batch 5, u04/u05]_
- [ ] **D6.** Two scanner artifacts (audit-internal candidate_ids, NOT real parameters → **no fix required**, optional
  candidate-list cleanup): `fit_stage_224_hat` (the `hat` param is a LaTeX-command fragment scanner-extracted from
  `\hat m_0`; mhat_0 is an undetermined lower-bounded branch factor at 224, not a posited value) and
  `fit_stage_200_stale_output` (empty bundle, zero constraints). SAME class as D1 (V_0f token). _[batch 6, u05/u01]_
- [ ] **D7.** Stage-232 `L_over_a` scanner dedup + mis-name: `fit_stage232_family1_geometry_normalization` and
  `fit_stage232_lambda_ell_geometry_literals` are duplicate scanner candidates for the same literal `20`, and the param
  name `L_over_a` is a stale_provenance_anchor mislabel of the **inverse thin-wall ratio a/ℓ=20** (NOT the aspect ratio
  L/a=37/20). Audit-internal, **no content fix**; optional dedup/rename note. _[batch 6, u04]_

## E. Disclosure-completeness (optional notes additions)
- [ ] **E1.** f'(0)=1 is documented only in the sympy script (`…stage065…:118`), absent from notes/ — optional notes mention. _[batch 2, u09]_
- [ ] **E2.** stage 091 notes (:107/:115) say "K_geom=3K_pole … forced" while 092 frees it ("no longer forced … unless");
  tension is OPENLY logged at 092, not concealed — optional cross-note. _[batch 3, u06]_
- [ ] **E3.** `D_W_bare` appears only in the stage137 sympy script, absent from notes137 (bundle self-flags severity:low)
  — optional notes mention (parallels E1). _[batch 4, u06]_
- [ ] **E4.** Batch-6 provenance/disclosure gaps (non-overclaim, optional fills): `Xi_prefactor_100` (the literal 100 in
  Ξ_χ=100·Θ_w·Λ_ℓ² has no in-stage derivation; FoM structure traces to Stage 066); `V_known` float genealogy below 222
  is undocumented (disclosed illustrative, so not an overclaim); `beta_path` (base-2) / `KW_t` (exp-ramp) / `Rratio_base`
  (=7/6) live only in scripts, not notes (script-only illustrative path/sample coordinates). _[batch 6, u01/u02/u04]_

## F. Script wording
- [ ] **F1.** `scripts/…stage247…` (":239 independently derived / :248 falsifiable") overstates the tautological
  λ_L closure — soften wording. _[batch 1]_
- [ ] **F2.** `scripts/…stage162…sympy_audit.py` hardcodes the r_F1 DECIMAL (1.77799353547498) instead of the closed
  form `sqrt(4107-100π²)/(10π)`; the de-transcription applied to 165/168/169 was not applied to 162 (already noted in
  `CHECKPOINT_CONSTANT_PROVENANCE.md`). Cosmetic script-form consistency, not an overclaim. _[batch 5, u01]_

## G. Phase A/B tooling cosmetics (lower priority — superseded by genealogy; opportunistic)
- [ ] **G1.** Seeded-chain coverage check should resolve aliases→canonicals before flagging missing (6 false alarms). _[exec plan §line 65]_
- [ ] **G2.** `frac_54_5` / mention-fallback extractor tightening. _[exec plan §line 65]_
- [ ] **G3.** Tier-2/3 value-fill for the 813 low-confidence target records (opportunistic). _[exec plan §line 65]_

---

## Provenance index (where each batch's items were first recorded)
- Batch 1: `verdicts/batch_c1_adjudication.md` (Items 1-2 + "Other batch-1 close-out fixes")
- Batch 2: `reports/_batch_c2_verdicts.md` §OPTIONAL (8 items)
- Batch 3: `reports/_batch_c3_verdicts.md` §OPTIONAL (8 items)
- Batch 4: `reports/_batch_c4_verdicts.md` §OPTIONAL (6 items: C12-C15, E3, D3)
- Batch 5: `reports/_batch_c5_verdicts.md` §OPTIONAL (5 items: A8, F2, D4, D5, + the verdict_logged γ₀ PARTIAL = A8)
- Batch 6: `reports/_batch_c6_verdicts.md` §OPTIONAL (D6, D7, E4; NO new A-class — band 200–239 had zero card/notes
  overclaim. Also logged there: the in-band internal_consistency overclaim-flavored slugs = completeness-critic seed
  targets, NOT Step-6 fix items.)
- Tooling: `docs/adversarial_audit_execution_plan.md` line 65
- _Batch 7 (free_choice 240-259) + the internal_consistency completeness-critic: append here as they land._
