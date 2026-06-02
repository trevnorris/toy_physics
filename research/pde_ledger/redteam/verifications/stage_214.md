---
unit_id: 214
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T12:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 214

## Per-finding outcomes

### F1 — paper_misalignment (notes_contradicts_script)

**Classification:** resolved

**What changed:**
F1 was resolved via the USER-authorized notes-only edit (`5·5·6 = 218` → `150`), applied by Codex and already orchestrator-reviewed as a single isolated boxed-expression edit. No `.py` change was authorized or needed for F1. Verified scripts-side: the diff (`stage_214_diff.patch`) touches ONLY the SymPy `.py` section-III block (lines 176–211 region); it contains no `bezout_projected`, no `218`, and no edit to the `!= 150` assertion. The current `.py` still reads `bezout_projected = deg_Crs * deg_Crt * deg_Sr` / `if bezout_projected != 150: raise AssertionError(...)` (lines 244–247), unchanged. SymPy log line 1165 prints `projected one-chart Bezout bound = 150`.

**Assessment:**
Correct and isolated. The `.py` F1 assertion was not modified; the script's internally-consistent `150` is preserved. `py_unchanged_for_f1: yes`.

### F2 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
New `.wl` created at `mathematica/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl` covering manifest M1–M7. Mathematica exec log exits 0 with all manifest items present and 32 `PASS:` lines; final line "All Stage 214 Mathematica audit checks passed."

**Assessment — INDEPENDENCE (load-bearing, not rubber-stamp):**
The `.wl` is a genuinely independent derivation, not a transliteration of the `.py`, on multiple axes:

1. **M4 eliminants via Resultant (strongest signal).** The `.wl` obtains the cross/square eliminants by `Resultant[liftedPolys[[i]], liftedPolys[[j]], y]` (eliminating `y` through the resultant) and then SEPARATELY confirms each equals the hand-built definition (`crossDefinitions = M_s·envNum_r − M_r·envNum_s`, etc.):
   `crossEliminants = {-Cancel[Resultant[liftedPolys[[1]], liftedPolys[[2]], y]/2], ...}`
   `expectZero["M4 C_rs resultant minus definition", crossEliminants[[1]] - crossDefinitions[[1]]]`
   The `.py` builds `Crs = Ms*Lr - Mr*Ls` directly. Re-deriving the same eliminant by resultant elimination and proving equivalence is exactly the anti-transliteration route the directive demanded.

2. **M1 derivative laws via a different decomposition.** The `.wl` builds the metric numerator as `den·D[linearK,var] − D[den,var]·linearK/2` and the envelope numerator as `den·D[quadEnvelope,var] − D[den,var]·quadEnvelope`, then forms `directNumerators = derivativeDen·D[phi,var]` and checks `directNumerators − stationaryNumerators == 0` AND `D[phi,var] − N/derivativeDen == 0`. This differentiates `Phi` directly rather than asserting against the `.py`'s pre-formed hand-expanded `M_r = (1+r²+s²+t²)k_j − r(...)`.

3. **Total degrees via Exponent/MonomialList.** `totalDegree[poly,vars]` uses `MonomialList`+`Exponent`+`Total`+`Max`, not `Poly(...).total_degree()`.

**M5 uses the COMPUTED degree product, not a literal.** `projectedProduct = crossDegrees[[1]] crossDegrees[[2]] squareDegrees[[1]]`, where `crossDegrees`/`squareDegrees` are measured via `totalDegree` on the resultant-derived eliminants. The `5*5*6` literal appears only in the check-name string and the comparison `projectedProduct - 5*5*6`; the value being tested is genuinely computed (= 150 from measured degrees 5,5,6) and would fail if any measured degree differed. No hardcoded `150` or `218`. The F1 dispute is not baked in. Reproduces card/appendix `54` (M3) and `(3,3,3,2)` (M2) independently.

### F3 — tautological_check

**Classification:** resolved

**What changed:**
The six definitional `expect_zero("cross-elimination identity ...")` / `expect_zero("square elimination identity ...")` asserts at lines 179–185 are GONE, replaced (diff lines 12–59) by:
1. **Non-vacuity:** `expect_true(f"{label} non-vacuous polynomial", not sp.Poly(eliminant, r, s, t).is_zero)` for all six eliminants (log lines 1143–1148, all `True`).
2. **Stationary-point vanishing:** a constructed numeric point evaluates all 10 polynomials — `F_r, F_s, F_t, F_Delta` AND `C_rs, C_rt, C_st, S_r, S_s, S_t` — to `0` (log lines 1149–1158, all `0`).
The degree checks (now lines 225–229; the directive's cited ~201–204) are retained.

**Assessment — FALSIFIABILITY (verified genuine, not a new tautology):**
The chosen point is the EXACT diagonal-isotropic envelope along the gradient-optimal ray, not an arbitrary point. Cross-checked the constants against the iso substitution with `ki=2,kj=3,kk=5,kl=7,H0=1`: off-diagonals match exactly (B=2·ki·kj=12, C=2·ki·kk=20, D=28, F=2·kj·kk=30, G=42, I=2·kk·kl=70), and the diagonal terms force a single consistent `u=116/3` (A=4−2u=−220/3, E=9−2u=−205/3, H=25−2u=−157/3, J=49−2u=−85/3 — all consistent). Because the four lifted polynomials `F_r=F_s=F_t=F_Delta=0` are independently asserted to vanish at this same point, the point is GENUINELY stationary (not arbitrary), so the eliminant-vanishing is meaningful. A sign error introduced upstream into `M`/`L`/`Delta` would change `F_*` and/or the eliminants and break the `=0` checks; the non-vacuity guard rules out a vacuous `0==0`. Genuinely falsifiable. `f3_now_falsifiable: yes`.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `1143–1148`: `C_rs/C_rt/C_st/S_r/S_s/S_t non-vacuous polynomial = True`
- `1149–1158`: all 10 lifted/eliminant polynomials `at diagonal-isotropic stationary point = 0`
- `1165`: `projected one-chart Bezout bound = 150`
- `1216`: `# exit_code: 0`

**Mathematica:** exit=0. Notable lines:
- `M4 C_rs resultant minus definition = 0` → `PASS` (resultant route matches hand definition)
- `M5 projected one-chart degree product = 150` → `PASS: M5 projected product minus 5*5*6`
- `M3 lifted degree product = 54` → `PASS`
- 32 `PASS:` lines; `# exit_code: 0`. Manifest M1–M7 all present.

**Output freshness:** confirmed. `.wl` mtime 2026-06-02 11:41:20; its `.txt` 11:59:40 (newer). `.py` mtime 11:38:56; its `.txt` 11:59:40 (newer). Both outputs regenerated post-fix.

## Material-change assessment

`material_change`: false. The F1 fix was notes-only prose (a typo correction, no carried value changed); the F3 replacement substitutes the verification mechanism for the same eliminants but introduces no new derived constant; the new `.wl` reproduces existing values (`54`, `(3,3,3,2)`, `150`) rather than introducing any. No downstream-relied result changed.

## Side observations (non-blocking)

- The Mathematica banner / SymPy "STAGE 197" labels are legacy stage-name strings inside otherwise-correct stage-214 audits; cosmetic only, not a finding.
- Section VI (winner/non-improvement integer ordering) remains SymPy-only by design per the directive (pure-logic transitivity model); the dual-engine requirement for this stage is satisfied by the M1–M7 re-derivation. Not blocking.

## Verdict justification

All three findings resolved. F1: the `.py` was confirmed untouched (notes-only, already reviewed). F2: a new Mathematica `.wl` independently derives M1–M7 via genuinely different decompositions — resultant-based eliminant elimination, direct differentiation of `Phi`, and Exponent/MonomialList degrees — and M5 carries the computed degree product (not a hardcoded literal), exiting 0 with 32 PASS lines. F3: the six definitional identities are replaced by non-vacuity plus vanishing at an exact, independently-stationary diagonal-isotropic point, which is genuinely falsifiable under an upstream sign error. Both engines exit 0, outputs are fresh, no regressions in the diff, no material change. Verdict: verified.
