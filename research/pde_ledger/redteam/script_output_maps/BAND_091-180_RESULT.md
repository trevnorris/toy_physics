# Band 091-180 — RESULT (SCRIPT/OUTPUT numbering pass)

Status: COMPLETE, verified. Content-keyed, never offset-swept. All edits digit-only
(strip-digits proof). Both engines exit 0 on every re-run.

## Character of this band (vs band 001-090)
NOT the +17 renumber zone. Overwhelmingly genuine current-canonical refs: of the 78
actionable candidate rows, **53 LEAVE + 4 VAR-LEAVE** (genuine cross-/forward-refs and
code identifiers), leaving only **4 distinct content-stale cross-refs** (9 token edits).
Classified by 3 read-only content agents (sub_D/E/F → map_D/E/F), then every proposed fix
was re-verified by the orchestrator against the cited stage's filename-stem + paper card +
corroborating siblings before applying. Two agent claims were corrected on verification
(see below).

## Applied (9 source edits across 5 files; `apply_band_091_180.py`, orchestrator-vetted)
- **#2** `stage112.wl:55` `Stage-92` → `Stage-109` (×1). stage109 = `linearized_branch_selection`
  owns the linearized (b,a0,a5) cross-check with preservation `a0/3 + 9 a5 == 0`; stage092 =
  `dynamic_geometry_obstruction` does not. Old-epoch number (92+17=109).
- **#3** `stage116` `Stage 98` → `Stage 115` (×6: py:68/73/79, wl:67/70/76). stage115 =
  `core_balance_compensation`, line 47 `gamma0_can = (1+r_c)/9`, owns the gamma0
  "compensation requirement"; stage098 = `family1_support_is_automatic` does not. Old-epoch
  number (98+17=115). The 2 output rows (py-out:14, wl-out:21) self-corrected on re-run.
- **#4** `stage134` `Stage 137` → `Stage 140` (×2: py:90, wl:101). "susceptibility closure"
  is owned by stage140 = `selfmatched_mouth_susceptibility` (card "Self-Matched Mouth
  Susceptibility Closure"), corroborated by stage139:10 ("established at Stage 140, not here")
  and stage142:37 ("closure (Stage 140 / Sigma_0 = (20/9) That_m^2)"); stage137 =
  `core_to_mouth_gain_map` does not. Comment-only (outputs unchanged).

## Output refresh (`rerun_refresh.py`)
Re-ran py {116,134} + wl {112,116,134} — 5/5 exit 0. Only stage116's 2 outputs changed
(`Stage 98`→`Stage 115` in the printed gamma0 line); stage112/stage134 outputs byte-identical
(comment-only edits → idempotent re-run, confirming no spurious drift).

## Corrected agent claims (verification caught these)
- Agent D proposed `stage100.py:13 stage 097 → 105` (chi_Q). REJECTED — see LEAVE #1.
- Agent D cited "stage116 already says Stage 115 for gamma0" as corroboration; FALSE (grep
  found zero `Stage 115` in stage116 — it is uniformly `Stage 98`). The 98→115 fix still
  holds on the independent +17 + stem(`core_balance_compensation`) + deliverable(line 47
  gamma0=(1+r_c)/9) evidence.

## LEAVE + FLAG for the user (out of label-only scope / content-attribution judgment)
1. **chi_Q "Stage 097" family** — `stage100.py:13` "stage 097 (single_normalization_defect)",
   `stage100.py:26`, `stage100.wl:30/34`, `stage101.py:12` (+ several 3-digit CANON-BACK
   siblings). The deliverable "chi_Q = 1 identification by DtN comparison" matches stage105 =
   `chiQ_fix_from_outgoing_dtn` better than the cited stage097 = `single_normalization_defect`.
   BUT the number AND its parenthetical STEM are internally consistent (both name stage097), so
   any fix would have to change the stem too — NOT a digit-only edit — and whether 097 (first
   identification) or 105 (exact fix) "owns" it is a content/red-team call, not numbering drift.
   → out of scope for this label-only pass. **Recommend: route to a content review, or leave.**
2. **stage140 stress-output filename anomaly** — `scripts/output/stage140_core_mouth_coevolution_status_stress.txt`
   is a byte-identical stray duplicate of the canonical `scripts/numerical/output/stage157_core_mouth_coevolution_status_stress.txt`
   (numerical capstone of stage157), misfiled under a stale `stage140_` filename prefix. Its
   internal "Stage 155/156/157" strings are CORRECT for what the file is. The drift is in the
   FILENAME (a rename/delete), not in-file text → outside this pass's label-only scope.
   Referenced only in historical redteam codex_logs (not live). **Recommend: rename to
   `stage157_...` or delete as a duplicate (user file-op decision).**
3. **stage105.py:8/31 "Stage 088/074"** — 088 = `loading_ratio_from_minimal_module` is the
   clear owner of the "minimal isotropic quadrupole module"; 074 = `family1_healing_lock` does
   not obviously own the pole scale `Omega_Q = 3 c_s/(2a)`. But "074" is 3-digit canonical-form
   (not an obvious old-epoch artifact) and its correct replacement (if stale) is unconfirmed →
   LEAVE per the CANON-BACK default. **Flag only.**

## Genuine refs LEFT (not stale — verified, default LEAVE)
- All 157→158 forward refs (deferred-to-Stage-158 per the card; stage158 = linear_defect_transport).
- 164 py:50 "General Stage 118 formulas" (stage118 = parent_core owns the K_s/K_q/λ/g overlap formulas).
- 167→163 / 168→165 / 169→168 / 165→164 / 178→177 backward cross-refs (each cited stem owns the named deliverable).
- 179 "Stage 176/160/161 slippage formulas" compound — all THREE sub-numbers verified correct.
- 137 "Stage 97/114" compound (genuine 2-digit refs to stages 097+114; notes pointers confirm).
- 163 "Stage 162/146" compound (162 = g_lower match; 146 = positive-deformation benchmark).
- 109 "downstream stages 110/111/112", 134→135 forward refs, 121 `stage99TubeLength`/161 `stage160Prefactor` identifiers (VAR-LEAVE), 144 wl:78 file-path pointer (file exists).

## Verification
- 7 files changed (5 source: 2 `.py` + 3 `.wl`; 2 outputs `.txt`).
- **strip-digits label-only proof: 7/7 digit-only** (HEAD vs working byte-identical after
  removing all digits; line counts unchanged) → zero equation/value/variable/logic bytes changed.
- Residual scan: 0 remaining `Stage-92`/`Stage 98`/`Stage 137` (old tokens) in edited files;
  new tokens present (Stage-109 ×1, Stage 115 ×8, Stage 140 ×2).
- 5/5 script re-runs exit 0. No notes/`.tex` touched. `git grep` for the wrong-root path typo: empty.

## Housekeeping (this commit)
- Reworded 2 self-referential lines in `BAND_001-090_RESULT.md` that quoted the wrong-root
  path-typo string, so `git grep` for it stays empty (matches the band-1 reword precedent 14603bb).
