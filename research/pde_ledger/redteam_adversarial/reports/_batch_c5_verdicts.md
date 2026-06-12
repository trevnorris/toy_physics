# Phase C — Batch 5 verdicts (targeted free_choice sweep, stages 160-199)

Batch: `batch_c5_freechoice_160_199` (PC6 targeted sweep). Roster: `_phase_c_batch5_roster.yaml`.
18 free_choice candidates (provenance_built; 19 param-checks, stage191 carries chi_0_star + F_star) ->
7 disjoint stage-cluster units, per-member sub-verdicts. Adversarial framing = the **posit-dressed-as-derived**
test (the stage_192 failure mode): is a *chosen* value labeled derived/forced/exact/non-tunable/a-theorem/
fixed-by in any paper card or per-stage notes line? `Benchmarks: []` is correct/expected for free_choice.

## HEADLINE: 6/7 units NO FATAL FLAW (17/18 candidates NO, 18/19 param-checks NO). 1 candidate (u02 gamma_0) returned YES, DOWNGRADED to a NON-FATAL PARTIAL by the defense + Claude+Codex adjudication -> verdict_logged. No FIND_STANDS; no user gate triggered.

This is the **most sensitive band of the whole sweep**: it contains the stage 192/193/194/195 chain (the literal
neighborhood of the canonical stage_192 posit-dressed-as-derived failure + the open PARTIAL punch-list A1), AND
two of the four ANSATZ_LEDGER §5 borderline calls (gamma_0=(1+r_c)/9 at 162; n=5 wall EOS at 165), AND a
continuation of the batch-4 u03 r_F1 family (the closed-form-4107 radius at 162/163/165/168). The one YES landed
exactly where the design predicted it might — the gamma_0 borderline — and the adjudication confirmed it is the
KNOWN, already-user-owned ANSATZ_LEDGER item surfacing as a stage-162 presentation-scoping nit, NOT a new fatal
contradiction.

Key cross-batch reassurance: as in batch 4, the overclaim-loaded vocabulary lives in candidate-id SLUGS
(`slopes_fixed_by_d01`, `g2_common_tangent`, `doubled_slope_factors`, `closed_form_4107`, `rexact`,
`canonical_point_block`) and is ABSENT from the actual cards/notes, which disclose every value as carried /
posited / free / chosen / open. The single exception (gamma_0) is a card-BADGE / notes-WORDING over-scope, not a
slug artifact, and is mitigated three ways (see u02).

## Per-unit verdict table

| unit | stages | members | verdict | one-line basis |
|------|--------|---------|---------|----------------|
| u01_rF1_closed_form ⭐SENSITIVE | 162/163/165/168 | r_F1 ×5 | NO | r_F1 honestly "Family-1 value" / "fixed lower compensated Family-1 values" / "renormalized canonical Family-1 point", always with "≈"; closed form sqrt(4107-100pi^2)/(10pi) and "4107" label the FORMULA, whose numeric content is the carried L/a=37/20 (origin 121:18/58, "arithmetically forced" by the posit at CHECKPOINT_CONSTANT_PROVENANCE:1202). BIT-IDENTICAL 1.77799353547498 across all anchors (162:155, 163:112/341, 165:198, 168:192) -> NO re-fit drift |
| u02_gamma0_ninth ⭐SENSITIVE | 162 | gamma_0 ×1 | **YES -> PARTIAL (non-fatal, verdict_logged)** | gamma_0=(1+r^2)/9 is grouped as one of "the two **exact** parent-family formulas" (notes162:70) and "Stages 115-116 **fix** the bare odd normalization" (notes162:58, boxed) under `\StatusExactClosure` (card162:5) with no ansatz tag. ADJUDICATED non-fatal: "exact" attaches to the family LAW/relation form (not the value), the stage's result uses ONLY gamma_0's functional form (delta-ln identity, an exact consequence of the posited family), "fix" back-references the upstream POSIT (origin 116:62 "A simple concrete realization is to take..."), and the derived-vs-posited status is the already-user-owned ANSATZ_LEDGER §5-item-1 borderline. Disposition contingent on that gate |
| u03_n5_wall_eos ⭐SENSITIVE | 165 | n_wall_eos ×1 | NO | n=5 consistently "the **frozen** `n=5` wall-EOS" (notes165:27/328/419), applied "On the frozen GNLS branch" — a chosen freeze of the parent-action polytropic index, carried from 118/062 ("frozen \(n=5\) EOS"); "exact" attaches to the drift-law reduction/identity, never to n=5 being forced. Independently logged ANSATZ_LEDGER §5-3 borderline |
| u04_slopes_doubling | 173/178 | D_01 + kappa_1 + d_r (3) | NO | D_01 explicitly OPEN ("compute ... on the actual moving-throat branch", notes173:446); what is "fixed by D_01" is D_21=-D_01/9, D_41=-D_01/27 (derived relations among OTHER slopes expressed in terms of the free D_01). kappa_1 "wall-baseline slope" (notes178:37), d_r free port-detuning slope (notes178:79); the factor-2 "doubling" is a derived identity nu_r=2(p_r-d_r), not a value claim |
| u05_g2_tangent | 187 | Lambda_1 ×1 | NO | Lambda_1 = Delta_lambda = ln(lam~_W/lam_W), first of "the five **free** finite log-ratios" you "Choose" (notes187:217), "**pure gauge/similarity directions**" (notes187:358); card `\StatusExactClosure` + pervasive "exact" attach to the closure/orbit-fibre/quotient THEOREM, not the value. Slug token "g2_common_tangent" has no literal in notes (denotes the tangent space T_id G_*) |
| u06_P2_chain ⭐SENSITIVE | 191/192/193/194/195 | chi_0_star + F_star + P_2 ×4 (6 checks) | NO | P_2 = l=2 multipole-order CONVENTION (the '2'), never a fitted scalar; every status header reads "**Exact within the carried** ... grouped-`P2` / `l=2` ... hierarchy" (191:5, 192:5, 193:5, 194:5, 195:5) -> explicitly conditional/carried. chi_0_star/F_star "the **carried** coherent-branch coefficients" (notes191:253), opaque symbols only. The "forced into a one-parameter form" (193:217) attaches to the carrier SHAPE after imposing one-pole, not to P_2/chi_0_star/F_star. A1 chi_Q issue NOT re-litigated (separate, already logged) |
| u07_mhat0_source_map | 195 | mhat_0 ×2 (dedup) | NO | mhat_0 "The **carried** natural ... source-map branch ... m_{hat0}=1+O(a^2/r^2)" (notes195:202), "in the strict point-particle **limit** m_{hat0}->1" (notes195:210); "exact" attaches to the factorization THEOREM mhat_0^2 chi_Q N_Q=1, not the value. CONFIRMED dedup artifact: the two slugs (`mhat0_natural_source_map` / `natural_source_map_mhat0`) are word-order permutations of the same param; second bundle self-identifies "Duplicate candidate for m_hat_0 at stage 195" |

## Orchestrator spot-checks (distrust-all-clean backstop — mandatory on the 4 SENSITIVE units)
- **u02 gamma0_ninth** (the YES): VERIFIED directly. Read `paper/stages/stage_162.tex` full — `\StatusExactClosure{}`
  (:5), the boxed quote :16 "Staying on the exact parent compensation family gives ... Delta_Q=0 at first order"
  (a result on the family), :27 caveat "the card is a derivation ledger entry, not an unconditional actual-branch
  theorem." Read `notes162` §1-§2 (:38-104): line :58 boxes "Stages 115-116 fix the bare odd normalization by
  gamma_0=(1+r^2)/9"; :64 "both ... controlled by the **same single parent variable** r"; :70 "Differentiate the
  two exact parent-family formulas." Confirmed the result (:72-104) uses only ln gamma_0 = ln(1+r^2)-ln9 (the
  functional form). The overclaim is a real card-badge / notes-wording over-scope, NOT a fabricated value-derivation.
  Both the defense agent and Codex (read-only `codex exec -s read-only`) independently returned NON-FATAL and agreed
  on the contingency. Adjudication record below.
- **u01 rF1_closed_form**: VERIFIED. r_F1 = 1.77799353547498 is BIT-IDENTICAL across notes162:155, 163:112+341,
  165:198, 168:192 (grep-confirmed) -> a fit would drift; this does not. The closed form sqrt(4107-100pi^2)/(10pi)
  labels the algebraic surd. Agent NO CONFIRMED. (Non-fatal: stage_162 sympy script hardcodes the decimal instead
  of the closed form — de-transcription applied to 165/168/169 but not 162; already tracked in
  CHECKPOINT_CONSTANT_PROVENANCE; logged OPTIONAL below.)
- **u03 n5_wall_eos**: VERIFIED. grep confirms n=5 uniformly "frozen `n=5` wall-EOS" (notes165:27/328/419); "exact"
  scopes the reduction/identity. Agent NO CONFIRMED.
- **u06 P2_chain**: VERIFIED. grep confirms every status header 191-195 opens "Exact within the carried ...
  hierarchy" — the badges are explicitly conditional. P_2 honestly the carried l=2 sector convention. Agent NO
  CONFIRMED.
All four SENSITIVE spot-checks ratify the agents; the YES on u02 was EARNED (a genuine over-scope, correctly
surfaced) and the downgrade is EARNED (defense + Codex both grounded the non-fatal disposition in file:line).

## Claude+Codex adjudication record (u02 gamma_0 YES -> PARTIAL)
- Adversary (audit agent): YES — gamma_0=(1+r^2)/9 (free_choice) presented at stage 162 as an "exact parent-family
  formula" / "fixed", under `\StatusExactClosure` with no ansatz tag = the stage_192 posit-dressed-as-derived pattern.
- Defense agent: REFUTED as fatal. All 4 defense lines CONFIRMED — (1) "exact" attaches to the family LAW, (2)
  the stage result uses only the functional form, (3) "Stages 115-116 fix" back-references the upstream POSIT
  (116:62), (4) gamma_0 is logged ANSATZ_LEDGER §5-1 borderline awaiting a human gate.
- Codex (`codex exec -s read-only`, prompt `tmp_prompts/_codex_b5_u02_adjudication.md`): independently VERDICT =
  NON-FATAL; confirmed all 4 sub-questions with file:line; AGREES the disposition is contingent on the existing
  ANSATZ_LEDGER §5 gate. "Stage 162 over-scopes the badge/wording, but its actual theorem is conditional on staying
  inside the declared parent family and uses gamma_0's functional form."
- CONVERGENCE: NON-FATAL PARTIAL. Not FIND_STANDS (fatal did not survive defense); not a NEW conceptual change
  (the derived-vs-posited call is already user-owned). -> verdict_logged, punch-list item A8 (gamma_0 scoping,
  same class as A3), contingent on the user's ANSATZ_LEDGER §5-1 layer-3 call. NO user gate triggered.

## OPTIONAL close-out items (NON-fatal; cosmetic/metadata/scoping; -> Step-6 punch list)
1. **u02 gamma_0 scoping (-> punch-list A8)**: scope the `\StatusExactClosure` badge at `paper/stages/stage_162.tex:5`
   + add an ansatz hedge at `notes162:38/58/70` flagging gamma_0's pure-scale-realization status (distinct from the
   genuinely-derived L_W/a relation it is grouped with). SAME class as A3 (stage_116 gamma_0 badge). CONTINGENT on the
   user's ANSATZ_LEDGER §5-item-1 call: if the user adopts the "derived (odd-preservation theorem)" reading, no edit
   is owed and the current wording is correct.
2. **u01 stage_162 sympy script form (-> punch-list F2)**: `scripts/...stage162...sympy_audit.py` hardcodes the
   r_F1 DECIMAL (1.77799353547498) rather than the closed form sqrt(4107-100pi^2)/(10pi); the de-transcription that
   was applied to 165/168/169 was not applied to 162. Already noted in CHECKPOINT_CONSTANT_PROVENANCE; cosmetic
   script-form consistency, not an overclaim.
3. **u07 mhat0 dedup (-> punch-list D4)**: `fit_stage195_mhat0_natural_source_map` and
   `fit_stage195_natural_source_map_mhat0` are word-order-permuted slugs for the SAME param mhat_0 at stage 195
   (second bundle self-identifies as duplicate) — mark as alias/dedup, not two standalone candidates. Audit-internal.
4. **u05 slug-vs-content (-> punch-list D5, informational)**: `g2_common_tangent` (187) has no literal in notes
   (denotes the tangent space T_id G_*); audit-internal candidate_id, NOT a paper line, no fix required.
5. **u04 slug-vs-content (-> punch-list D5, informational)**: `doubled_slope_factors` (178) slightly mislabels — the
   factor-2 lives in the derived identity nu_r=2(p_r-d_r), not in kappa_1/d_r being "doubled factors." Audit-internal
   slug, NOT a paper line, no fix required.

## Coverage / no-silent-caps log
- DEEP-AUDITED this batch: all 18 free_choice candidates (19 param-checks) in band 160-199, stages 162/163/165/168/
  173/178/187/191/192/193/194/195. 0 free_choice candidates in band 160-199 left un-audited. (Bands 160/161, 164,
  166/167, 169-172, 174-177, 179-186, 188-190, 196-199 carry no free_choice candidates per the enumeration.)
- NOT in scope this batch (by PC6 design): internal_consistency and published_target candidates in 160-199 (covered
  by layer-1 red-team x2; the fit-vs-derive test has near-zero surface on internal_consistency). The
  internal_consistency completeness-critic spot-check is the planned safety net AFTER the free_choice sweep completes.

## Overall
Batch 5 extends the layer-2 result into the MOST sensitive band of the sweep and behaved exactly as the falsification
design intends: the one place a fatal overclaim could plausibly hide — the gamma_0 ANSATZ_LEDGER borderline — is the
one place an adversary fired YES, and the defense + Codex adjudication resolved it to a known, already-user-owned,
non-fatal presentation nit rather than a new fatal contradiction. The stage 192-195 canonical-failure neighborhood
came back clean on the P_2 convention (every badge is explicitly "Exact within the carried ... hierarchy"), the r_F1
closed-form family showed zero re-fit drift (bit-identical across four anchors), and n=5 is honestly "frozen." NO
program revision required; the 5 OPTIONAL items are cosmetic/metadata/scoping (A8 contingent on the user's existing
ANSATZ §5 gate). Manifest: 17 candidates -> audited, 1 (gamma_0) -> verdict_logged (walked audit_pending -> audited
-> defense_pending -> defended -> adjudicating -> verdict_logged per the protocol).
