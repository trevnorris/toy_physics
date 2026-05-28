---
unit_id: 161
batch: IV.6
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-28T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 161

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:** sympy lines 70-85 and wl lines 65-76 now derive `depsg_*` from the exact `eps_gamma = 9*gamma0/(1+rc) - 1` via `sp.diff` / `D[]`. The assertion substitutes `dgamma0 -> (1+rc)*dln_gamma0/9` (the branch-point physical relation) into `depsg_branch` before subtracting `(dln_gamma0 - drc/(1+rc))`.

**Assessment:** Correct. Orchestrator deviated from the directive's typo (`depsg_direct`) and used `depsg_branch`, which is the right choice — `depsg_branch` is already evaluated on the branch `gamma0 -> (1+rc)/9` so the `dgamma0`-substitution test still exercises the branch-point relation. A sign flip in `epsg_exact` propagates through `sp.diff`. Non-tautological.

### F2 — tautological_check

**Classification:** resolved

**What changed:** sympy lines 46-52 and wl lines 37-43 now substitute `eps_kappa -> eps*deps_k, eps_gamma -> eps*deps_g, rc -> rc_star + eps*drc` into the exact `BW` (line 42 / 33) and take `sp.diff(..., eps).subs(eps, 0)` / `D[..., eps] /. eps -> 0`.

**Assessment:** Correct. `dBW` is now derived from the exact `BW`; sabotaging line 42 propagates to the assertion. Output shows `dB_W = (d_eps_g - d_eps_k)*(r_c_star + 1)/9` / `((depsG - depsK)*(1 + rcStar))/9`, matching the expected linearization.

### F3 — mathematica_transliteration

**Classification:** resolved

**What changed:** Banners on sympy line 37 and wl line 26 both read `"STAGE 161 — D/N SIMILARITY SLIPPAGE DECOMPOSITION"`. F1 and F2 fixes give the two engines distinct algebra: SymPy uses `sp.diff` + `sp.simplify` chains; Mathematica uses `D[]` + `FullSimplify` + `PolynomialRemainder` (already engine-distinct in the even-defect block).

**Assessment:** Correct. Banner labels propagated to both `.txt` outputs. Load-bearing blocks no longer share construction strategy.

## Exec log assessment

**SymPy:** exit=0. Notable lines: `exact similarity-defect decomposition = 0`, `linearized slippage law = 0`, `d eps_gamma = d ln gamma0 - d ln(1+r_c) = 0`, `difference identity = 0`, all collapse and preservation checks pass.

**Mathematica:** exit=0. `PASS:` printed for all 9 assertions including `linearized slippage law`, `d eps_gamma = d ln gamma0 - d ln(1+r_c)`, `difference identity`. Numeric `(1+r_F1^2)/9 = 0.4623623346878688...` matches SymPy.

**Output freshness:** confirmed. Output mtimes (1779989470, 1779989532) postdate script mtimes (1779989463, 1779989467).

## Material-change assessment

`material_change`: false. Edits only harden assertions; no derived numerical value or symbolic carry-forward formula changed. All printed quantities (`B_W`, `dB_W`, `eps_kappa`, `Upsilon_Pi`, `Delta_Q`, `(1+r_F1^2)/9`) match the original audit.

## Side observations (non-blocking)

None.

## Verdict justification

All three findings resolved with correct, non-tautological substitutions. SymPy and Mathematica both PASS with exit 0; outputs refreshed post-edit. Orchestrator's catch of the `depsg_direct`/`depsg_branch` directive typo is sound — using `depsg_branch` preserves the branch-point physics exercised by the `dgamma0 -> (1+rc)*dln_gamma0/9` substitution. F2 now linearizes the exact `BW`; F3 banners corrected and load-bearing blocks engine-distinct. Verdict: `verified`.
