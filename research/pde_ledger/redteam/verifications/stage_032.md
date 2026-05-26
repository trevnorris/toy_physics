---
unit_id: 032
batch: II.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26T00:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 032 (batch II.1 v2)

## Per-finding outcomes

### F1 — insufficient_verification (independent (v.e_-)^2 check)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py:167-194` — inserts the prescribed block: builds `v2 = sp.Matrix([kappa0, kappa1])`, constructs `M_loaded = sp.diag(A, A + DK) - alpha0 * (v2 * v2.T)`, calls `M_loaded.eigenvects()`, sorts by eigenvalue using the `(alpha0=1, A=1, DK=1)` probe rule, normalizes the lower eigenvector to unit length, computes `s_check = ((v2.T * e_minus)[0])**2`, and asserts `expect_zero("s_check - s_minus_nat (independent (v.e_-)^2 construction)", ...)` and `expect_zero("lam_minus_sym - lambda_minus (independent eigenvalue construction)", sp.simplify(lam_minus_sym - lambda_minus.subs(subs_nat)))`. The second assertion uses the `.subs(subs_nat)` "adjusted second assertion" form the directive offered.
- `mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl:152-174` — inserts the prescribed Mathematica analog: `v2 = {{kappa0},{kappa1}}`, `mLoaded = DiagonalMatrix[{a, a + dK}] - alpha0 * (v2 . Transpose[v2])`, `{eigVals, eigVecs} = Eigensystem[mLoaded]`, lower-index pick via `First[Ordering[N[eigValsProbe]]]` with `probeRule = {alpha0->1, a->1, dK->1}`, unit-normalized `eMinus`, then `expectZero["s_check - s_minus_nat (independent (v.e_-)^2 construction)", ...]` and an eigenvalue residual assertion.

**Assessment:**
The SymPy block is byte-identical to the directive's prescribed code. The Mathematica block is the directive's prescribed Mathematica idiom. Both new assertions print `0` / `PASS` (sympy log lines 111-112; mathematica log lines 71-74). The checks are non-tautological: `e_-` is constructed fresh from `eigenvects()`/`Eigensystem`, a route not used elsewhere in the script, and compared to the closed-form `s_minus_nat` / `lambda_minus` carried forward from stages 028-031. A bug in the upstream `s_minus` formula would surface here. The symbolic route worked cleanly in both engines, so the numeric-probe fallback was not needed. No collateral edits beyond the prescribed insertion.

### F2 — mathematica_transliteration (independent Eigensystem + Stage 15.5 cross-check)

**Classification:** resolved

**What changed:**
- `mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl:157` introduces an `Eigensystem[mLoaded]` call that does not appear in the SymPy script (SymPy uses `Matrix.eigenvects()`, a structurally distinct routine).
- `mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl:181-188` inserts the prescribed Stage 15.5 cross-check: `nProdIndep = FullSimplify[(sCheck/kappa0^2) * (beta0 * sCheck / lamMinusIndep), ...]` and `expectZero["nProdNat - nProdIndep (independent eigenvector path)", ...]`.

**Assessment:**
The directive said the F1 `Eigensystem` block doubles as F2's fix for Stage 15.4, plus the additional `nProdNat - nProdIndep` assertion for Stage 15.5. Both are present. The Mathematica `Eigensystem` route (with `Ordering`, `Flatten`, `FullSimplify`, native Mathematica idiom) is a genuinely different algebraic path from SymPy's `Matrix.eigenvects()` route. The Stage 15.5 cross-check uses the eigenvector-derived `sCheck` and `lamMinusIndep` rather than retracing the closed-form `sMinus`/`lamMinus`, so it adds independent coverage. Both new assertions pass (mathematica log lines 85-86). No collateral edits.

### F3 — insufficient_verification (drop algebraic-identity scaffolding)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py` — the three assertions `delta_kappa^2 + 4*Kprod - sigma^2 (natural)`, `s_minus_nat - s_minus_nat_simplified (interior identity)`, and `s_minus_nat at (alpha0=1, DK=1) interior point` are deleted (diff patch lines 80-93 show removal). Only `expect_zero("mhat_-^2(alpha=0) - 1", ...)` and `expect_zero("limit_{alpha->oo} mhat_-^2 - 11/9", ...)` remain in this region, alongside the new F1 block.
- `mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl` — the parallel `expectZero` calls for `delta_kappa^2 + 4*kappaProd - sigmaSym^2 (natural)`, `s_minus_nat - s_minus_nat_simplified (interior identity)`, and `s_minus_nat at (alpha0=1, dK=1) interior point` are deleted (diff patch lines 11-25 show removal).

**Assessment:**
The directive specified F3 deletion was conditional on F1's symbolic route passing cleanly (rather than falling back to numeric probes). The exec logs confirm F1's symbolic check passes (`s_check - s_minus_nat = 0` symbolically in both engines), so deletion is the correct branch. The four removed assertions no longer appear in the exec logs. The remaining `limit_{alpha->oo} mhat_-^2 - 11/9` assertion is preserved (sympy log line 114; mathematica log lines 79-80) per the directive's "Keep line 191 / Keep line 177" instruction. No collateral edits.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- Line 111: `s_check - s_minus_nat (independent (v.e_-)^2 construction) = 0` — new F1 assertion passes symbolically.
- Line 112: `lam_minus_sym - lambda_minus (independent eigenvalue construction) = 0` — second F1 assertion passes.
- Line 114: `limit_{alpha->oo} mhat_-^2 - 11/9 = 0` — preserved endpoint check.
- Line 121: `All Stage 15 checks passed.`
- F3-removed assertions are absent (no `delta_kappa^2 + 4*Kprod`, `s_minus_nat_simplified`, or `interior point` substrings appear).

**Mathematica:** exit=0. Notable lines:
- Lines 71-72: `s_check - s_minus_nat (independent (v.e_-)^2 construction) = 0` / `PASS`.
- Lines 73-74: combined-label eigenvalue check `= 0` / `PASS` (label is malformed but residual is zero — see side observation 1).
- Lines 85-86: `nProdNat - nProdIndep (independent eigenvector path) = 0` / `PASS` — F2 Stage 15.5 cross-check passes.
- Lines 78, 90: benign `Limit::alimv` warnings (Mathematica's `Limit` ignores assumptions on the limit variable; noted as benign in the v1 audit's "Engine cross-check" section, since residuals match SymPy).
- F3-removed assertions are absent.

**Output freshness:** The saved `.txt` outputs under `scripts/output/` (mtime 2026-05-21 17:25:53) and `mathematica/output/` (mtime 2026-05-21 17:26:05) are OLDER than the post-fix script mtimes (2026-05-25 23:46:25 sympy, 2026-05-25 23:48:32 mathematica). They were NOT regenerated post-fix. The authoritative post-fix transcripts are the exec logs in `redteam/exec_logs/stage_032_{sympy,mathematica}.log` (timestamps 2026-05-26 00:18-00:19), which are fresh and confirm both scripts exit 0 with all assertions passing. Per the user's policy memo against parallel `exec-*` calls, the `output/*.txt` artifacts should be regenerated by a separate direct `python3` / `math -script` invocation if downstream consumers need them. Flagging as a side observation, not a verification blocker.

## Material-change assessment

`material_change: false`.

All edits are additions of independent witnesses (F1, F2) or removals of redundant algebraic-identity scaffolding (F3). None changes a derived numerical or symbolic result that downstream units consume. The closed-form expressions for `lambda_minus`, `s_minus`, `s_minus_nat`, `mhat_sq`, `Sigma`, `Xi`, `alpha`, `Nprod_nat`, `kappa0`, `kappa1`, `sigma` are unchanged. Downstream stages that import or reference these quantities see the same values.

## Side observations (non-blocking)

1. The Mathematica `expectZero` label at `.wl:172` concatenates the SymPy label and the prescribed Mathematica label with a semicolon: `"lam_minus_sym - lambda_minus (independent eigenvalue construction); lam_minus_indep - lamMinus (independent eigenvalue construction)"`. The assertion semantics and residual (=0) are correct; only the printed label is malformed. Cosmetic, does not affect verification.
2. The saved `output/*.txt` artifacts are stale (older than the post-fix scripts). The exec logs under `redteam/exec_logs/` are fresh and authoritative, but if any downstream consumer reads `output/*.txt` directly it will see the pre-fix content.

## Verdict justification

All three findings are `resolved`. Both engines now exercise the headline claim `mhat_-^2 = s_-/kappa_0^2` via an independent `Eigensystem`/`eigenvects` construction of the lower eigenvector `e_-` and verify `(v.e_-)^2 = s_minus_nat` symbolically (F1). The Mathematica script now has a genuinely non-transliterated path through Stages 15.4-15.5 via `Eigensystem` plus the `nProdIndep` cross-check (F2). The redundant `(a-b)^2 + 4ab = (a+b)^2` algebraic scaffolding is removed cleanly (F3). Both exec logs exit 0 with every assertion printing `0` / `PASS`. No regressions in the diff. No new findings introduced.
