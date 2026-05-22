---
unit_id: 032
batch: II.1
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 032

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
Stage 15.5 of both scripts was rewritten as directed.
- `scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py:193-210` now drops the fresh `kappa0_sq` and `P0_minus` symbols, applies `subs_nat` to `(beta0*s_minus/lambda_minus)` and `(s_minus/kappa0**2)`, then asserts `Nprod_nat(alpha0=0) - beta0*kappa0**2/A == 0` and `sp.limit(Nprod_nat, alpha0, sp.oo) == 0`.
- `mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl:179-189` mirrors with `p0MinusNat`, `nProdNat`, `nProdAt0 - beta0*kappa0^2/a`, and `Limit[nProdNat, alpha0 -> Infinity]`.

**Assessment:**
The new assertions are non-tautological. At `alpha0=0`, `subs_nat` makes `R = sqrt(DK^2) = DK` (DK>0), so `lambda_minus(0) = (A+B-DK)/2 = A` and `s_minus(0) = (sigma + DK*delta_kappa/DK)/2 = (sigma+delta_kappa)/2 = kappa0^2`, producing the target `beta0*kappa0^2/A`. The check would fail if the cross-term in `s_minus` or the `B = A + DK` definition were wrong. The alpha0→∞ limit check requires `lambda_minus ~ -alpha0*sigma` to actually diverge — also non-trivial. The unused `kappa0_sq` and `P0_minus` symbols were removed as required. Both engines saved-output PASS. The saved sympy output (lines 115-116) and mathematica output (lines 82-87) show the two new assertion names with residual = 0.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl:104-121, 135` adds an independent `LinearSolve`-based Schur derivation: builds `zVec = {z1,z2}`, solves `kInt . y = Transpose[bMat] . z`, forms `bMat . ySol`, reads off `sigmaMatViaSolve` via coefficient extraction in `z1, z2`, then asserts `Schur via Inverse vs LinearSolve == 0` and a separate `sigmaMatViaSolve - [Xi I + alpha vv^T] == 0`. The original `Inverse[kInt]` path is preserved alongside, as the directive's drop-in replacement intended.

**Assessment:**
This is a methodologically distinct algebraic path — `LinearSolve` solves a 4-component linear system rather than constructing a 4x4 inverse, and the Schur image is read off as coefficients of two formal external sources `z1, z2`. Both routes now agree with each other AND with the rank-1+identity target (three substantive checks where there was previously one). Minor deviation from the directive: Codex kept the existing `v`, `i2`, `kInt`, `bMat` definitions (lines 91-102) rather than re-declaring them inside the replaced block; substantively this is identical to the directive's intent because the directive's replacement code itself re-declares them with the same values. Saved Mathematica output (lines 52-60) shows all three Stage-15.3 checks PASS with zero matrices.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
Both scripts gained three new Stage 15.4 assertions between the alpha0=0 and alpha0→∞ endpoint checks.
- `scripts/...py:168-190` and `mathematica/...wl:153-176`:
  - `delta_kappa^2 + 4*Kprod - sigma^2 (natural)` — verifies the identity that makes Stage 15.4's algebra clean.
  - `s_minus_nat - s_minus_nat_simplified (interior identity)` — uses `R_nat = sqrt(DK^2 + 2*alpha0*DK*delta_kappa + alpha0^2*sigma_sym^2)` and `(sigma + (DK*delta_kappa + alpha0*sigma^2)/R_nat)/2` and asserts symbolic equality with `s_minus.subs(subs_nat)`.
  - `s_minus_nat at (alpha0=1, DK=1) interior point` — same identity evaluated at a concrete interior parameter point.

**Assessment:**
Non-tautological. The "interior identity" check expands `(DK + alpha0*delta_kappa)^2 + 4*alpha0^2*Kprod` to `DK^2 + 2*alpha0*DK*delta_kappa + alpha0^2*(delta_kappa^2 + 4*Kprod)` and similarly expands the numerator cross-term to `DK*delta_kappa + alpha0*(delta_kappa^2 + 4*Kprod)`. Equality with `s_minus_nat_simplified` REQUIRES that `delta_kappa^2 + 4*Kprod` collapses to `sigma^2` under `subs_nat` — i.e. the kappa-product identity is what makes the simplified closed form correct. If the original `s_minus` had a wrong coefficient on the cross term `((DK + alpha0*delta_kappa)*delta_kappa + 4*alpha0*Kprod)` (e.g. `4 -> 5` on the `alpha0*Kprod` term), then `s_minus_nat_simplified` (which presumes `delta_kappa^2 + 4*Kprod = sigma^2`) would no longer match `s_minus.subs(subs_nat)` and the interior identity would fail symbolically. The concrete `(alpha0=1, DK=1)` check is strictly redundant with the symbolic identity (if symbolic holds, concrete holds) but the directive asked for both and the assertion is per-spec. Saved sympy output (lines 107-109) and mathematica output (lines 68-73) confirm all three new assertions are zero residual.

The eigenvalue/eigenvector structural derivation discussed in the audit report was explicitly allowed by the directive as optional ("if Codex cannot reconstruct M from the script's existing content, it should ADD ... at minimum one interior-point numerical check"). The chosen identity+interior-point route is per the directive's accepted fallback.

## Exec log assessment

**SymPy:** exit=0 (inferred from saved output ending with "All Stage 15 checks passed." and `expect_zero` raising on non-zero residuals). The dedicated `redteam/exec_logs/stage_032_sympy.log` is absent, but the saved `scripts/output/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.txt` is the canonical evidence and is fresh post-fix.

Notable saved-output lines:
- `Nprod(alpha=0) - beta0 * kappa0^2 / A = 0` (sympy txt line 115)
- `limit_{alpha->oo} Nprod_nat = 0` (sympy txt line 116)
- `delta_kappa^2 + 4*Kprod - sigma^2 (natural) = 0` (line 107)
- `s_minus_nat - s_minus_nat_simplified (interior identity) = 0` (line 108)
- `s_minus_nat at (alpha0=1, DK=1) interior point = 0` (line 109)
- `All Stage 15 checks passed.` (line 117)

**Mathematica:** exit=0 (inferred from saved output ending with "All Stage 15 checks passed." and the script's `fail[]` calling `Exit[1]`). The dedicated `redteam/exec_logs/stage_032_mathematica.log` is absent, but `mathematica/output/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.txt` is fresh post-fix.

Notable saved-output lines:
- `PASS: Schur via Inverse vs LinearSolve` (math txt line 53)
- `PASS: sigmaMatViaSolve - [Xi I + alpha vv^T]` (line 60)
- `PASS: delta_kappa^2 + 4*Kprod - sigma^2 (natural)` (line 69)
- `PASS: s_minus_nat - s_minus_nat_simplified (interior identity)` (line 71)
- `PASS: s_minus_nat at (alpha0=1, dK=1) interior point` (line 73)
- `PASS: Nprod(alpha=0) - beta0 * kappa0^2 / a` (line 83)
- `PASS: limit_{alpha->oo} Nprod_nat` (line 87)
- Two `Limit::alimv` warnings appear (lines 75, 85); per the directive these are acceptable provided `expectZero` reports PASS, which it does.

**Output freshness:** confirmed. Script mtimes are `1779405849` (epoch, 2026-05-21 17:24); sympy output mtime is `1779405953` (~17:25:53, +104 s); mathematica output mtime is `1779405965` (~17:26:05, +116 s). Both saved outputs are newer than the corresponding script files, so they were regenerated after Codex's edits.

## Material-change assessment

`material_change`: false.

The edits add assertions and an alternative algebraic path; they do not change any derived closed form. `s_minus`, `lambda_minus`, `Sigma`, `Xi`, `alpha`, `mhat_-^2`, `kappa0`, `kappa1`, `sigma` retain their pre-fix values (verified by comparing the saved `Sigma`, `Xi`, `alpha`, `mhat_-^2`, `kappa0/kappa1/sigma` values in the post-fix output against the values quoted in the auditor's report Section "Engine cross-check"). No downstream unit needs to be re-audited on numerical grounds. The orchestrator will still mark units > 032 as `upstream_stale: true` as part of the standard policy.

## Side observations (non-blocking)

- The orchestrator's per-unit exec logs `redteam/exec_logs/stage_032_sympy.log` and `redteam/exec_logs/stage_032_mathematica.log` are missing; only the diff patch is present. Verification proceeded from the fresh saved `.txt` outputs (script mtime < output mtime). If a future audit expects dedicated exec-log files, the orchestrator may want to capture them alongside the diff.
- The `s_minus_nat at (alpha0=1, DK=1) interior point` check is strictly redundant with the symbolic `(interior identity)` check that precedes it (any symbolic identity holds at every concrete point). It does no harm and matches the directive verbatim, but it is not adding independent coverage. Flagging only for the auditor's awareness — not a basis for needs_rework.
- The Mathematica Stage 15.3 retains the direct `Inverse[kInt]` route alongside the new `LinearSolve` route. Per F2's directive this is intended (the directive's replacement code keeps both), but a reader skimming the file might initially miss that the `LinearSolve` path is the "independent" derivation; a brief header comment naming which route is canonical could help, though this is style not substance.

## Verdict justification

All three findings are resolved. F1's Stage 15.5 rewrite replaces the multiplicative-associativity tautology with two substantive checks (alpha0=0 endpoint with explicit `beta0*kappa0^2/A` target, and alpha0→∞ limit) under the natural D/N substitution. F2's Stage 15.3 adds a `LinearSolve`-based Schur derivation algebraically distinct from the 4x4 inverse, plus a cross-check assertion. F3's Stage 15.4 adds the kappa-product identity check and an interior algebraic-form consistency check (with a concrete-point spot-check on top). Saved outputs for both engines show every new assertion as PASS / residual=0, and both end with "All Stage 15 checks passed." The closed forms downstream units depend on (`Sigma`, `Xi`, `alpha`, `s_minus`, `lambda_minus`, `mhat_-^2`, `kappa0`, `kappa1`, `sigma`) are unchanged. Verdict: `verified`.
