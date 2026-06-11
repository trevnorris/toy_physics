# Step 4C genealogy fan-out — running findings ledger

Compact record of notable per-batch findings from the genealogy agents (full structured data lives in the per-candidate bundles + `_synthesis/<batch>/`). Newest waves appended. Severities: HIGH = real provenance defect; needs_triage = scanner artifact / ambiguity; low = genealogy completeness note.

## Cross-cutting (applies to multiple batches)
- **em_projected citation-root convention (stages 004–020):** these projected-Maxwell EM-extension stages have NO per-stage notes under `research/pde_ledger/notes/`; their derivation notes live at the **repo root** `/var/projects/toy_physics/notes/em_projected/step_*.md`. Agents correctly cited them as `notes/em_projected/...` (resolves from the repo root; transliteration-verified accurate on step_03:68). So Phase B citations have TWO roots: `notes/stages/...` = pde_ledger-relative, `notes/em_projected/...` = repo-root-relative. Both pass ingest (path starts with `notes/`). Document for any verifier. STAGE_PROVENANCE_INDEX points some 018/019/020 sources at `notes/em_projected/step_16/17/18_*.md` that do NOT exist in the checkout — those early stages fall back to downstream notes / CITATION_MAP.
- **Scanner false-positives** consistently logged as `provenance_gap` (not fabricated values): LaTeX macros (`\paragraph`, `StatusExact`), SymPy helper names (`assert_nonzero`, `boundary_value`, `tautological_check`, `det_eff`), prose words (`to`), meta/staleness flags (`stale_output`), and subscript/condition-RHS mis-extractions (a literal `0`/`1` that is a positivity/window condition, not a fitted value). These are correctly NOT classified.
- **Numbering-drift stale anchors** recur (the +17 EM-renumber pattern): notes cross-refs like "Stages 22–23"/"24–25" and "Stage-18"/"Stage 86" point at pre-renumber stage numbers. Logged `stale_provenance_anchor` (mostly low; see HIGH below for the load-bearing one).

## Wave 1 (batches 01–06, 08, 09) — 267 candidates → provenance_built

### HIGH findings (real provenance defects — Phase C / triage candidates)
- **batch_09 / Pe_fail_chi misattribution (stage 089) — self_contradictory_origin + stale_provenance_anchor (HIGH).** `Pe_fail_chi = 11220.5441626259` (and `Pe_suff_chi`) is anchored by stage-089 script + CHECKPOINT_CONSTANT_PROVENANCE to `scripts/output/...stage082...`, but stage 082 contains neither value; true origin is **stage 078** (`Theta_w^(chi)/3.62605617972939e-4`), carried via stage 080; meanwhile the 089 notes say "Stage 086." Three records disagree. This is the χ_Q-class bug the audit exists to catch.

### published_target chain (the load-bearing external fit) — batches 02, 03, 04
- The GR **quadrupole 54/5** coefficient (`= 2/5` of `γ_GR = 2G/(5c⁵)`, equivalently 27/20) is the single dominant external-fit anchor across stages **019–038**. Origin = **stage 022** (`NQ_target = solve(m̂²·P0·Γ5_port == γ_GR, P0)`); `P0_target`/`N_Q`/`K_req`/`λ_req`/`α_star`/`mhat_rad`/`Lambda` at 019/022/023/025/026/030/031/038 are all back-solved to it (8× `derive_vs_fit_mismatch` HIGH-severity-as-finding in batch_02, more in 03/04).
- **Honestly disclosed, NOT overclaimed:** the cards label these "a branch test, not a proof of branch realization" / `StatusOpen` / "the radiative demand scale" — so NO `paper_card_overclaim` was raised on the 54/5 chain. The fit is real and labeled as such. (This is the audit confirming the disclosure is adequate — the key layer-2 question.)
- batches 01, 05, 06, 08: **zero published_target** — foundational geometry-lift/projected-Maxwell (01), internal algebraic closure chains (05, 06, 08). All "derived" claims are genuine internal_consistency.

### free_choice ansätze surfaced (honest posits)
- Wall-action constitutive coeffs μ_η/T_w/T_Ω/K_η (stage 001, "new ansatz"); n=5 frozen EOS exponent (parent action, origin phase0/stage062); V_0 wall amplitude + f'(0)=1 normalization (065); ε_blk=0 unblocked limit; λ_μ=1; α_r=10 wall depth; the minimal-module coefficients c0=c_contact=3/4 (imported at 088/090 as "upstream-fixed" but actually DERIVED downstream at 091–094 — a forward-reference, needs_triage).

### batch tallies (internal_consistency / published_target / free_choice / gap-only)
- b01: 20 / 0 / 7 / 6   b02: 14 / 9 / 8 / 4   b03: 21 / 12 / 1 / 2   b04: 26 / 3 / 0 / 4
- b05: 31 / 0 / 1 / 3   b06: 28 / 0 / 4 / 1   b08: 26 / 0 / 4 / 3   b09: 19 / 0 / 15 / 0

## Wave 2 (batches 10–17) — 253 candidates → provenance_built (batch_11 needed a candidate_id case-fix; see ingest note below)

### LOAD-BEARING: the r_F1 / 37-20 dependency chain (batches 13, 15, 16, 17)
- The entire Family-1 canonical-constant block — **r_F1 ≈ 1.77799, g_star, Σ0_can, Π_can/Π_star, S_can, T_can, the 1/√(1+r_F1²) transport coeffs, Υ_Π, Δ_Q** — is `internal_consistency` at the SOLVE level but bottoms out on ONE posited input: **L/a = 37/20** (`r_F1 = √((12/π²)(37/20)² − 1)`, origin stage 121, freeze from stage 073). 37/20 is `free_choice` (notes: "reference-branch numerical freeze, not a new theorem"); NO external published number backs it. So these are **derived-from-a-posit, not first-principles** — agents logged ~12 `derive_vs_fit_mismatch` (needs_triage) across 121/148/152/154/155/156/158/161 cautioning against reading "forced by r_F1" as fundamental. Adjudicator question flagged: project memory groups 37/20 with the external 54/5 target class, but the in-band evidence says posit, not fit — decide free_choice vs published_target at curation.
- **r_F1 IS classified `published_target` by batches 16/17** (they treat 37/20's class as external-flagged) but `free_choice` by batch 13/15 — a cross-batch classification inconsistency on the SAME parameter to reconcile in 4C-curation (the chain view). The honest answer hinges on whether 37/20 itself is posited (free_choice) or an external fit (published_target).

### Σ0_can / T_can / Π_can precision drift (batches 16, 17) — the targeted transposition, surfaced organically
- `Sigma0_can` notes/fixture = ...168867 vs published appendix part04 = ...168**876**; `T_hat_can` ...7613 vs ...7624; `Pi_can` ...79002 vs ...790087. Notes and the numerical fixture AGREE; only the published appendix differs (sub-1e-13). `paper_card_overclaim`/`stale_provenance_anchor` (Σ0_can = needs_triage). This is the Σ0 …867/…876 drift the audit specifically targeted.
- Stage attribution: the `_can` quartet (Sigma0_can/Pi_can/S_can/T_can) belongs to **stage 156** (renormalized canonical point), NOT 155 — stage 155 explicitly DEFERS the retuned traction. 4 `stale_provenance_anchor`.

### chi_Q chain (batches 10, 11, 12) — confirms prior fix + ONE new live overclaim
- **chi_Q = 1 is `internal_consistency`, origin stage 105** (matched to the stage-104 outgoing-DtN h₂⁽¹⁾ fingerprint); the external GR Γ_5 = 2G/5c⁵ fit lives at **stage 106**, `published_target`, **disclosed honestly** (cards say "target values"/StatusOpen → no overclaim). The "097 vs 105" contradiction is confirmed RESOLVED in current sources (b052471).
- **NEW live finding (needs_triage): stage_104.tex:13 paper_card_overclaim** — carries the χ_Q=1 boilerplate that was corrected on the 100/101 cards but MISSED on 104 (104 keeps χ_Q symbolic, only derives the fingerprint).

### other
- **gamma_0 = (1+r_c)/9 framing tension (batch 12, needs_triage):** the ansatz-value catalog lists it as the first postulated value, but stage-115 notes present it as a DERIVED odd-preservation condition (free at 114, constrained at 115).
- numbering-drift stale anchors persist in prose (batches 12/13: "Stage-90/91/92" → 107/108/109; the +17 drift); all low.
- batches 10/12/13/14/15 mostly internal_consistency; **zero published_target** in the 108-155 internal-closure bands except the r_F1/37-20 question above.

### batch tallies (internal_consistency / published_target / free_choice / gap-only)
- b10: 27 / 3 / 3 / 3   b11: 27 / 3 / 1 / 2   b12: 28 / 0 / 6 / 0   b13: 18 / 0 / 14 / 1
- b14: 24 / 0 / 6 / 3   b15: 28 / 0 / 6 / 0   b16: 36 / 0 / 0 / 0   b17: 30 / 3 / 0 / 0

### ingest note
- batch_11: agent lowercased 3 candidate_ids (chiQ→chiq, copied from slugified bundle filenames). Pre-ingest validator mapped them back to the manifest's exact case (case-insensitive match). Added as a standing pre-ingest guard for remaining waves.

## Wave 3 (batches 18–25) — 263 candidates → provenance_built

### CORRECTION: the V_known barrier (222–224) is ILLUSTRATIVE, not a published_target (batches 23, 24, 25)
- Phase A flagged V_known≈1.18190922 as a likely external barrier import. The stage-222 notes **explicitly** call it illustrative ("the true local barrier requirement must be pulled back from the actual same-charge branch"; line 436) — no external source, no back-solve. Classified `free_choice`. The two stage-222 sample-slice blocks (20 params) are likewise illustrative placeholders (23 low `paper_card_overclaim`). **So the seeded `chain_barrier_222_224` is NOT an external-fit chain** — reclassify in curation. RQ_req (228) back-computes from V_known so it inherits the illustrative status (`published_target`+ambiguous, but the target is illustrative not published).

### The real external fit stays P0_target = 54G/5c⁵ (GR quadrupole)
- `published_target` consistently at stages **189/193/195/197/223/224** (P0_target / mhat0-normalization / T_quad), all tracing UPSTREAM to **stage 022** origin (`γ_GR = 2G/5c⁵`), all **honestly disclosed** ("universal target"/"target scale"/StatusOpen) → mismatch findings logged but no real overclaim. This is the legitimate layer-2 published-target.

### New HIGH / load-bearing findings
- **stage 192–195 paper-card boilerplate overclaim (batch_21, 4 HIGH `paper_card_overclaim`):** cards 192–195 share a block asserting "collapses Packet A to Δ_Q=χ_Q−1 / N_Q=χ_Q⁻¹", but at **stage 192** (an orbit/quotient projector stage) χ_Q/Δ_Q/N_Q appear NOWHERE in the notes (scanner conflated χ_{0,*} with χ_Q). Δ_Q:=χ_Q−1 is genuinely defined only at 195; χ_Q=1 fixed only at 194. Real card overclaim to fix.
- **stage163 R_q=1/4 mislabel (batch_18, HIGH self_contradictory_origin):** candidate `fit_stage163_Rstar_quarter` scanned value 1/4 labeled `r_F1` but is actually R_q=R_*=1/4 (derived from stage-119 balance). Scanner mislabel.

### Framing-risk (adjudication) — "5PN injection" is internal, not external
- stage 232 "Known 5PN Data Injection" (batch_25) + even-gate H_even/K_1 "imported strict 5PN" (batch_24): termed as external imports but the values (Theta_w → derived at 077; even-gate coeffs 1/9,2/3,1/27 → internal γ_0=(1+r_c)/9 algebra) are **ledger-internal**. The "5PN stack" = the EM-projected derivation chain, NOT a publication. Agents found NO external number matched. `constraint_kind_ambiguous`/needs_triage — adjudicate the framing, but no real external fit.

### corroborates curation inputs
- batch_19: `fit_stage_177_p_0` value '38' = scanner artifact (mtime "01:38" string) → confirms PRUNE from chain_quad_54_5 (curation input). batch_24: `fit_stage224_grouped_signature_exact` b_P0_slope = internal_consistency → confirms the genuine coverage-add. batch_18: r_F1/37-20 chain again (free_choice vs published_target cross-batch split persists — central curation call).

### batch tallies (internal_consistency / published_target / free_choice / gap-only)
- b18: 28 / 0 / 7 / 0   b19: 26 / 0 / 3 / 4   b20: 29 / 3 / 1 / 1   b21: 26 / 2 / 8 / 2
- b22: 24 / 1 / 4 / 4   b23: 25 / 0 / 23 / 0   b24: 28 / 5 / 5 / 0   b25: 20 / 1 / 9 / 1

### ingest note (wave 3)
- batch_21: agent rewrote `chiQ`→`chi_q` (case + added underscore) in 1 candidate_id. Extended the pre-ingest guard to alphanumeric-normalized matching (strip non-alnum + lowercase, unique-match only). Standing guard for wave 4.

## Wave 4 (batches 26–29) — 108 candidates → provenance_built. ALL 915 canonicals now provenance_built.

### published_target external imports confirmed (the last two of the three real external matches)
- **λ_L = 0.26971918 (stage 247) — back-solve to a recorded benchmark (HIGH).** Fixed by inverting the compiler so V_eff(r_soft) reproduces the recorded Session-I softening point **V_eff^sess = 1.74701126** (the deepest anchor, HIGH). Classified `published_target` across all 5 framings. The notes ARE honest ("audit point, not a general theorem"), but the SCRIPT wording ("independently derived... pinned to paper figures" / "falsifiable closure", lines 239/248) is misleading since the forward check is tautological by construction. Ξ_turn=0.34437471 / λ_th=0.42826825 (stage 248) are carried Session-II hardcodes (`published_target`; mild card overclaim — presented as diagnostics without disclosing they're imported). [Corroborates the curation coverage-add of stage248 Ξ_turn.]
- **CODATA 1836.15267343 proton/electron mass ratio (stage 250) — `published_target` (HIGH).** Imported as the "proton-proxy" benchmark (m_s = μ_η). Honestly disclosed as a proxy (no derive-vs-fit overclaim), but it IS an external published constant fit into the ledger. Surfaces under 3 (candidate,parameter) entries.
- **s_0 pinned to γ_crit (stage 251) — `derive_vs_fit_mismatch` (needs_triage).** s_0 is DEFINED as the internal growth rate √(κ_V/μ_η)=1/t_collapse,0, but its benchmark NUMBER is pinned to the Session-IV envelope-closure threshold γ_crit≈6.94311167 (which the notes themselves say is NOT a microscopic theorem); the script asserts equality only to 1e-6.

### ansatz roots (free_choice) at the calibration end
- c0_min = 3/4 (stage 240) — load-bearing posited "minimal isotropic conservative quadrupole module"; the whole ρ_α=4/3, ζ_req=1/3, Π_tr chain derives from it; card frames downstream ratios "Exact" without surfacing the 3/4 root (paper_card_overclaim).
- f_lat=3/4, μ_η=1, K_turn=2.73855812, the Session-I/II/III/IV/V benchmark readbacks (helicity exports, χ_peak, E_max, transit times, γ_crit, t_cross, etc.) — all un-derived benchmark conventions, generally honestly labeled "benchmark-only"/"reported".

### stage 252/253 (batch_29) — no external fit
- The final calibration slice has ZERO published_target: values are ledger-internal reduced-model outputs (s_c, s_0, r_turn, λ_ref carried from 248/251) or by-hand session conventions. The 4 long-id audit-file-name candidates + 3 LaTeX-token false-positives (mathcal, max, subsection) correctly logged as GAPs (bounded-filename ingest worked).

### batch tallies (internal_consistency / published_target / free_choice / gap-only)
- b26: 21 / 0 / 14 / 2   b27: 9 / 11 / 13 / 0   b28: 11 / 3 / 21 / 0   b29: 11 / 0 / 17 / 9

## ============ 4C COMPLETE: consolidated external-fit (published_target) inventory ============
The whole-ledger fit-vs-derive map. THREE genuine external published-target matches (all honestly disclosed in the notes; the audit's job was to confirm disclosure adequacy):
1. **GR quadrupole P0_target = 54·G·c_s⁵/(5·a⁵·c⁵)** (= 2/5 Einstein coeff × 27 branch factor), origin **stage 022** (`γ_GR = 2G/5c⁵`), carried 019–038 / 189–197 / 223–225. Disclosed (StatusOpen / "universal target"). The 54/5 ≡ 27/20 coefficient is THE dominant fit anchor.
2. **Recorded softening benchmark V_eff^sess = 1.74701126** via the λ_L=0.26971918 back-solve, **stage 247**. Disclosed in notes; script wording overclaims.
3. **CODATA proton/electron mass ratio 1836.15267343**, **stage 250** (proton-proxy μ_η/m_s). Disclosed as proxy.

NON-matches that LOOKED external but are NOT (audit corrections):
- **V_known≈1.18190922 (barrier 222–225)** = ILLUSTRATIVE placeholder, the notes say so explicitly. NOT a published_target. → reclassify seeded chain_barrier_222_224.
- **"Known 5PN data injection" (stage 232) / "imported strict 5PN" even-gate (224)** = the EM-projected internal derivation chain, NOT a publication. Values (Θ_w@077, even-gate 1/9,2/3,1/27) are ledger-internal.

Posited ansatz roots everything else rests on (free_choice): **L/a=37/20** aspect ratio (073/121 → r_F1 → the whole Family-1 canonical block); **c0=3/4** (240); **n=5** frozen EOS; **γ_0=(1+r_c)/9**; the wall-action constitutive coeffs (001); the Session benchmark readbacks (245–253).

CENTRAL ADJUDICATION FOR PHASE C / USER (the one cross-batch classification tension): **is L/a = 37/20 a `free_choice` (posited toy-geometry aspect ratio) or a `published_target` (external fit)?** The in-band genealogy evidence says free_choice (notes: "reference-branch numerical freeze, not a new theorem"; no external number named, no back-solve); project memory groups 37/20 with the external 54/5 class. Batches split (13/15/18/26 → free_choice; 16/17 → published_target). This is the highest-leverage FIND_STANDS-class question and is flagged for Phase C adversarial verification + user gate.

HIGH-severity provenance defects found (Phase C / triage candidates, NOT yet FIND_STANDS):
- Pe_fail_chi misattribution (089) — self_contradictory_origin (3 records disagree: 082 vs 086 vs true 078/080).
- stage 192–195 card boilerplate overclaim — χ_Q/Δ_Q/N_Q asserted where absent (esp. 192).
- stage163 R_q=1/4 mislabeled r_F1 — self_contradictory_origin (scanner mislabel).
- NEW: stage_104.tex:13 paper_card_overclaim — χ_Q=1 boilerplate not fixed (like the 100/101 fix).
- λ_L script "independently derived" wording (247) vs tautological forward-check.
