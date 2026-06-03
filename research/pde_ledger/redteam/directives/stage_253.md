---
unit_id: 253
batch: VIII.1
created_at: 2026-06-03T00:00:00-06:00
findings_count: 4
stop_cold: null
applied: true
verification_status: pending
needs_user_resolution: false
applied_at: 2026-06-03T13:32:44-06:00
findings_applied: 4
findings_blocked: 0
---

# Codex directive — unit 253

Apply F1's authorized notes edits (see `## RESOLVED — F1` below) AND each non-paper_misalignment finding (F2, F3, F4) in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

F1 (`paper_misalignment`) has been RESOLVED by the user (2026-06-03): correct the NOTES-file typos to the script values — see the `## RESOLVED — F1` block. Do NOT edit the published card `paper/stages/stage_253.tex` (it is clean) or the script/`.wl` literals to "fix" F1; the ONLY authorized prose edits are the two notes lines named in that block.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch `paper/stages/stage_253.tex` or any prose document, EXCEPT the two notes lines authorized in `## RESOLVED — F1` (user-approved 2026-06-03). Otherwise the red-team only modifies scripts.

---

## RESOLVED — F1 (paper_misalignment, USER-APPROVED 2026-06-03)

**Subtype:** value_mismatch — **NOTES-ONLY; the published paper card is UNAFFECTED.**

**Corrections (apply both verbatim):**
- `notes/stages/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_sympy_audit.md:274` — change `\frac{187.23361317}{\zeta_{\rm ep}}.` → `\frac{119.23361317}{\zeta_{\rm ep}}.`
- same notes file `:419` — change `10.95423247\,` → `10.95423248\,`

**Why (verified by orchestrator):** The published card `paper/stages/stage_253.tex` is only 140 lines and contains NEITHER `187.23361317` NOR lines 274/419 — the audit agent's `stage_253.tex:274`/`:419` citations were misattributions. Both disputed values live ONLY in the notes file. (i) The notes' own adjacent values `\gamma_{\rm lat,safe}^{\rm eq}≈65.45193926` (notes:266) and `s_c≈0.5489386551` (notes:259) give `65.45193926/0.5489386551 = 119.2336…`; both SymPy (py:192) and the existing Mathematica `.wl` (wl:254) independently assert `119.23361317476524` → cross-engine corroborated; `187` is a leading-digit notes typo (shared trailing `…23361317`). (ii) The a_int stiffness coefficient is `10.95423248 = 4·K_turn` in both engines (py:196 / wl:258); notes:419 `10.95423247` is a 1-ulp last-digit typo.

**Scope:** Edit ONLY notes:274 and notes:419. Do NOT touch `paper/stages/stage_253.tex` or the script/`.wl` literals. Codex applies the notes edits; Claude reviews.

## Applied: F1

- files_changed:
  - `notes/stages/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_sympy_audit.md`
- summary: Corrected the two authorized notes typos to `119.23361317` and `10.95423248`.
- deviation: none

---

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_mathematica_audit.wl` (whole file)

**Issue:** The `.wl` is a section-for-section, line-for-line transliteration of the `.py`: identical decomposition (I-VI mirroring py sections 1-6), identical intermediate-variable choreography (gammaLatSafeEq → tCrossPhys → gammaLatTurnPhys → thresholdLambda → thresholdProduct, exactly as py:49-60), identical substitution order (`/. VprimeTurnAbs -> KTurnSym rTurn/lambdaRef^2` then `/. lambdaPhys -> aInt/2`, mirroring py:104-105), and a 1:1 mapping of each `assert simplify(a-b)==0` to one `expectZero[name, a-b]`. This is an echo of the SymPy algebra, not an independent second-engine derivation, and cannot catch a shared algebra error.

**Required change:**
Re-author the `.wl` to reach the SAME final claims by a non-mirrored route. Do not change which claims are verified (engine agreement must stay clean). Specifically:
1. Build `r_turn_phys` and `chi_lambda` from the physical-unit-map premises directly. Define the harmonic potential `V[r_] = (1/2) kEff r^2`, take `D[Log[V[r]], r]` to obtain `2/r` natively (rather than transcribing the pre-cancelled `2/r_phys`), evaluate at `r = rTurnPhys`, and multiply by `lambdaPhys` to recover `chiLambdaLattice`; then assert it equals `2 lambdaRef/rTurn`. This makes the χ check a derivation, not a transcription.
2. Build `kEffReq` by setting up the force-match equation symbolically and `Solve`-ing for the stiffness:
   `Solve[ EStar (lambdaRef/lambdaPhys) VprimeTurnAbs == kEff rTurnPhys, kEff ]`, with `rTurnPhys = lambdaPhys rTurn/lambdaRef`, and confirm the solution equals `EStar lambdaRef^2 VprimeTurnAbs/(lambdaPhys^2 rTurn)`. (The `.py` writes the closed form directly at py:101; the `.wl` should instead derive it via `Solve`.)
3. For the threshold/Korringa/screening identities, reach them through a different grouping (e.g. construct `thresholdLambda` by dimensional substitution of `gammaLatTurnPhys = zetaEp lambdaEpOmegaD` solved for `lambdaEpOmegaD`, rather than copying py's `gammaLatTurnPhys/zetaEp`).
4. Keep all `expectZero`/`expectNear` final claims and the benchmark block intact (numbers unchanged — they already assert 119.23…/etc., consistent with both engines).

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 253` and confirms (i) the script still exits 0 with all checks passing, (ii) the construction now contains at least one native `Solve[...]` / `D[Log[...]]` step absent from the `.py`, so the route is no longer a line-by-line mirror.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_mathematica_audit.wl`
- summary: Reworked the Mathematica derivation to use native `Solve[...]` balances, `D[Log[V[r]], r]` for the harmonic geometry, and regrouped threshold/Korringa/screening constructions while preserving the final checks.
- deviation: none

---

## F3 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_sympy_audit.py:114` and `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_mathematica_audit.wl:176`

**Issue:** `r_turn_phys` is defined as `sp.simplify(lambda_phys * r_turn / lambda_ref)` (py:99) and then asserted equal to the byte-identical expression at py:114, i.e. `expr − expr ≡ 0`. The assertion cannot fail for any physics. The Mathematica twin (define wl:150, assert wl:176) is the same tautology. The dedicated radius-map check therefore proves nothing (the map is only indirectly exercised by the χ/k_eff cancellations downstream).

**Required change:**
Replace the tautological self-comparison with a round-trip consistency check against the *inverse* length map that §3.1 / py:101 / wl:153 actually use (`r = (lambda_ref/lambda_phys) r_phys`). Composing forward∘inverse must return `r_turn`:

- py:114 — change
  ```
  assert sp.simplify(r_turn_phys - lambda_phys * r_turn / lambda_ref) == 0
  ```
  to
  ```
  assert sp.simplify((lambda_ref / lambda_phys) * r_turn_phys - r_turn) == 0
  ```
- wl:176 — change
  ```
  expectZero["turning-point radius map", rTurnPhys - lambdaPhys rTurn/lambdaRef];
  ```
  to
  ```
  expectZero["turning-point radius round-trip", (lambdaRef/lambdaPhys) rTurnPhys - rTurn];
  ```

(Note for F2/F3 interaction: if F2 is applied, keep the F3 round-trip form; do not reintroduce the literal self-comparison.)

This residual is still 0 (forward∘inverse = identity) but the assertion is no longer `expr − expr`; it fails if the forward and inverse maps are inconsistent. Self-test confirmed: the expression depends on lambda_ref, lambda_phys, r_turn (no vacuous derivative), and reduces to 0 only by genuine cancellation through the inverse map.

**Verification command:**
The verifier runs `redteam exec-sympy 253` and `redteam exec-mathematica 253`; confirms the new round-trip assertion appears at py:114 / wl:176 (no longer restating `lambda_phys*r_turn/lambda_ref`) and both scripts exit 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_mathematica_audit.wl`
- summary: Replaced the radius-map self-comparison with the requested inverse-map round-trip check in both engines.
- deviation: none

---

## F4 — stale label

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_mathematica_audit.wl:65`

**Issue:** The Mathematica banner reads `banner["STAGE 236 — PHYSICAL CALIBRATION AND MATERIAL-THRESHOLD COMPANION"];` — a stale stage number; this unit is 253. Cosmetic (Print string only) but appears in the saved transcript (`output/...mathematica_audit.txt:11`) and mis-labels the artifact.

**Required change:**
wl:65 — change `banner["STAGE 236 — PHYSICAL CALIBRATION AND MATERIAL-THRESHOLD COMPANION"];` to `banner["STAGE 253 — PHYSICAL CALIBRATION AND MATERIAL-THRESHOLD COMPANION"];`

**Verification command:**
After `redteam exec-mathematica 253`, the transcript banner line reads "STAGE 253".

## Applied: F4

- files_changed:
  - `mathematica/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_mathematica_audit.wl`
- summary: Updated the Mathematica audit banner from stage 236 to stage 253.
- deviation: none
