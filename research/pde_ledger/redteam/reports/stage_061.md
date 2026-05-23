---
unit_id: 061
batch: III.3
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-22T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 061 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage061_microscopic_gain_thresholds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage061_microscopic_gain_thresholds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage061_microscopic_gain_thresholds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage061_microscopic_gain_thresholds_mathematica_audit.txt`

## What the script claims to verify

Starting from the Stage-41 endpoint functions `Delta_0(eta, kappa)` and `Delta_inf(eta, kappa)` with `alpha = sqrt(kappa)`, the scripts verify five claims of the operator phase diagram: (1) the microscopic coupling `Xi_micro = chi·Lambda^2·L^2/T_X` collapses to `kappa·G_micro` with `G_micro = chi·Lambda^2/K_X` under the kappa-defining substitution `T_X = K_X·L^2/kappa`; (2) the threshold parameter expressions `chi_fail`, `chi_suff`, `Lambda^2_fail`, `Lambda^2_suff` reduce consistently to `(K_X/Lambda^2)·G_fail`, etc., under the same substitution; (3) the soft-support limit `kappa -> 0+` gives `Delta_0 -> 1/2`, `Delta_inf -> 1`, hence `kappa·G_fail -> Pe_req` and `kappa·G_suff -> 2·Pe_req`; (4) the compliant-mouth limit `eta -> infinity` gives `Delta_0 -> (1-sech alpha)/kappa`, `Delta_inf -> tanh(alpha)/alpha`, hence `G_fail^inf = Pe_req/(alpha·tanh alpha)` and `G_suff^inf = Pe_req/(1-sech alpha)`; (5) the stiff-support limit of (4) at `z = alpha -> infinity` gives `z·G_fail^inf -> Pe_req` and `G_suff^inf -> Pe_req`.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 40-41 | `expect_zero(Xi_micro.subs(TX, KX*L^2/kappa) - kappa*G_def)` | yes (claim 1; ledger item 1) |
| A2 | sympy | 59-60 | `expect_zero(chi_fail.subs(...) - (KX/Lam^2)*G_fail)` | yes (claim 2) |
| A3 | sympy | 61-62 | `expect_zero(chi_suff.subs(...) - (KX/Lam^2)*G_suff)` | yes (claim 2) |
| A4 | sympy | 63-64 | `expect_zero(Lam2_fail.subs(...) - (KX/chi)*G_fail)` | yes (claim 2) |
| A5 | sympy | 65-66 | `expect_zero(Lam2_suff.subs(...) - (KX/chi)*G_suff)` | yes (claim 2) |
| A6 | sympy | 73 | `expect_zero(Delta0_k0 - 1/2)` | yes (claim 3) |
| A7 | sympy | 74 | `expect_zero(Delta_inf_k0 - 1)` | yes (claim 3) |
| A8 | sympy | 76-77 | `expect_zero(limit(kappa*G_fail, kappa, 0+) - Pe_req)` | yes (claim 3, ledger item 3) |
| A9 | sympy | 78-79 | `expect_zero(limit(kappa*G_suff, kappa, 0+) - 2*Pe_req)` | yes (claim 3, ledger item 3) |
| A10 | sympy | 86-87 | `expect_zero(Delta0_eta_inf - (1-sech alpha)/kappa)` | yes (claim 4) |
| A11 | sympy | 88-89 | `expect_zero(Delta_inf_eta_inf - tanh(alpha)/alpha)` | yes (claim 4) |
| A12 | sympy | 95-96 | `expect_zero(G_fail_inf - Pe_req/(alpha·tanh alpha))` | yes (claim 4, ledger item 4) |
| A13 | sympy | 97-98 | `expect_zero(G_suff_inf - Pe_req/(1-sech alpha))` | yes (claim 4, ledger item 4) |
| A14 | sympy | 106 | `expect_zero(limit(z·G_fail_inf_z, z->oo) - Pe_req)` | yes (claim 5, ledger item 5) |
| A15 | sympy | 107 | `expect_zero(limit(G_suff_inf_z, z->oo) - Pe_req)` | yes (claim 5, ledger item 5) |
| B1 | mathematica | 48 | `expectZero["Xi_micro - kappa G_micro", ...]` | yes (claim 1) |
| B2 | mathematica | 54 | `expectZero["chi_fail from G_fail", ...]` | yes (claim 2) |
| B3 | mathematica | 55 | `expectZero["chi_suff from G_suff", ...]` | yes (claim 2) |
| B4 | mathematica | 56 | `expectZero["Lambda^2_fail from G_fail", ...]` | yes (claim 2) |
| B5 | mathematica | 57 | `expectZero["Lambda^2_suff from G_suff", ...]` | yes (claim 2) |
| B6 | mathematica | 63 | `expectZero["Delta0 soft-support limit - 1/2", ...]` | yes (claim 3) |
| B7 | mathematica | 64 | `expectZero["Delta_inf soft-support limit - 1", ...]` | yes (claim 3) |
| B8 | mathematica | 65 | `expectZero["kappa G_fail soft-support - Pe_req", ...]` | yes (claim 3) |
| B9 | mathematica | 66 | `expectZero["kappa G_suff soft-support - 2 Pe_req", ...]` | yes (claim 3) |
| B10 | mathematica | 77 | `expectZero["Delta0 eta->Infinity formula", ...]` | yes (claim 4) |
| B11 | mathematica | 78 | `expectZero["Delta_inf eta->Infinity formula", ...]` | yes (claim 4) |
| B12 | mathematica | 79 | `expectZero["G_fail^(inf) formula", ...]` | yes (claim 4) |
| B13 | mathematica | 80 | `expectZero["G_suff^(inf) formula", ...]` | yes (claim 4) |
| B14 | mathematica | 85 | `expectZero["stiff-support compliant-mouth: sqrt(kappa) G_fail -> Pe_req", ...]` | yes (claim 5) |
| B15 | mathematica | 86 | `expectZero["stiff-support compliant-mouth: G_suff -> Pe_req", ...]` | yes (claim 5) |

## Findings

No findings.

Notes on attacks attempted (all rebuffed):

- **Tautology probe of A1-A5 / B1-B5.** These checks substitute `T_X -> K_X L^2/kappa` in the threshold parameters and assert equality with the macroscopic-gain form. The substitution is algebraic, so the residual identically reduces to 0 — but the substitution itself encodes a non-trivial physical claim (the Stage-1 definition of `kappa = K_X L^2/T_X`). Once that substitution is accepted as a definition (which the script asserts as ledger item 1), all five identities follow by algebra. This is a consistency-of-definitions check, not a tautology in the harmful sense (defining `x = expr` then asserting `x == expr`). The fact that A1 and B1 are stated as the ledger's first theorem makes the algebraic check substantive rather than vacuous.
- **Sign and factor probe on the soft-support limit.** Manually expanding `Delta_0` for small `alpha`: numerator `eta(cosh alpha - 1) ~ eta·alpha^2/2`, denominator `kappa(alpha sinh alpha + eta cosh alpha) ~ alpha^2·(alpha^2 + eta) ~ alpha^2·eta` to leading order. Ratio -> `1/2`. ✓ Matches A6/B6. Similar check for `Delta_inf` -> `1`. ✓ A7/B7.
- **Sign/factor probe on eta->infinity limits.** `Delta_0 -> eta(cosh alpha-1)/(eta·kappa·cosh alpha) = (cosh alpha-1)/(kappa·cosh alpha) = (1-sech alpha)/kappa`. ✓ `Delta_inf -> (eta/alpha)sinh alpha/(eta cosh alpha) = tanh(alpha)/alpha`. ✓ Match A10-A11/B10-B11.
- **G_fail^inf consistency.** `Pe_req/(kappa·Delta_inf^inf) = Pe_req/(kappa·tanh alpha/alpha) = Pe_req·alpha/(kappa·tanh alpha) = Pe_req/(alpha·tanh alpha)` since `kappa = alpha^2`. ✓ A12/B12.
- **G_suff^inf consistency.** `Pe_req/(kappa·Delta_0^inf) = Pe_req/(kappa·(1-sech alpha)/kappa) = Pe_req/(1-sech alpha)`. ✓ A13/B13.
- **Stiff-support attack.** `z·G_fail^inf|alpha=z = z·Pe_req/(z·tanh z) = Pe_req/tanh z -> Pe_req` as z->oo since tanh z -> 1. ✓ A14/B14. `G_suff^inf|alpha=z = Pe_req/(1-sech z) -> Pe_req` since sech z -> 0. ✓ A15/B15.
- **Branch probe.** With `kappa > 0`, `eta > 0`: `Delta_inf` numerator `cosh alpha - 1 + (eta/alpha)sinh alpha > 0`; denominator `alpha sinh alpha + eta cosh alpha > 0`. So `Delta_inf > 0` strictly. Similarly `Delta_0 > 0`. Thus `G_fail`, `G_suff` are well-defined and finite throughout the claimed parameter region. No division-by-zero or sign-flip risk.
- **Symbol assumption probe.** `kappa, eta, Pe_req` are declared `positive, real` in sympy and `Reals` plus positivity constraints in Mathematica. `alpha = sqrt(kappa)` derives positivity from kappa>0. The Mathematica warning `Limit::alimv` is benign (a documented behavior where assumptions on the limit variable are dropped) and does not affect correctness since direction is explicitly specified (`FromAbove`) or the limit is at infinity.
- **Engine cross-check.** Both engines produce identical residuals (`= 0`) on all corresponding checks. Final formulas (e.g., `G_fail^inf = Pe_req/(sqrt(kappa) tanh(sqrt(kappa)))` in sympy vs `peReq Coth[Sqrt[kappa]]/Sqrt[kappa]` in Mathematica) are equivalent under standard hyperbolic identities. The simplified `G_suff^inf` differs in displayed form (sympy: `Pe_req cosh(sqrt(kappa))/(cosh(sqrt(kappa)) - 1)`; Mathematica: `peReq (1 + 1/(Cosh[Sqrt[kappa]]-1))`) but both equal `peReq/(1-sech alpha)` — the `FullSimplify` residual against the target is 0 in both engines.
- **Freshness.** Sympy script mtime Apr 1 12:39, sympy output mtime May 11 12:44 (output newer). Mathematica script mtime May 11 11:56, output mtime May 11 12:54 (output newer). No stale_output.

## Independent-derivation check (Mathematica)

The Mathematica script does mirror the SymPy script's variable choreography closely: same closed-form definitions of `Delta_0` and `Delta_inf` (drawn from "Stage 41" per the docstring), same intermediate expressions for `G_fail`, `G_suff`, `chi_fail`, `chi_suff`, `Lambda^2_fail`, `Lambda^2_suff`, and the same substitution `T_X -> K_X L^2/kappa`. The order and content of `expectZero` checks is also identical.

However, the unit's claim is not "derive Delta_0 and Delta_inf from a PDE" (that is Stage 41's job) — it is "the limits and algebraic relations of these pre-derived endpoint functions match the operator phase diagram." The endpoint functions are *inputs* to this unit, not its outputs. Each engine independently evaluates the limits (`sp.limit` vs `Limit[...]`) and simplifications using its own algebra, and they reach the same conclusions. So while the surface choreography is parallel (which is unavoidable when both engines start from the same pre-derived inputs and compute the same limit relations), the limit-evaluation and simplification engines are independent.

The Mathematica `gFailInf` is simplified to `peReq Coth[Sqrt[kappa]]/Sqrt[kappa]` (output line 49) whereas SymPy emits `Pe_req/(sqrt(kappa)*tanh(sqrt(kappa)))` (output line 38) — both equivalent but obtained via different simplification paths. Similarly the stiff-support substitution differs structurally: SymPy uses `subs(alpha, z)`, Mathematica uses `kappa -> z^2`, demonstrating non-identical choices of representation. This is sufficient independence for the stage's claim. Not a `mathematica_transliteration` finding.

## Engine cross-check

Both engines pass all 15 corresponding assertions with residual 0. Side-by-side selected results:

| Quantity | SymPy | Mathematica |
|---|---|---|
| Delta_0 | `eta*(cosh(sqrt(kappa)) - 1)/(kappa*(eta*cosh(sqrt(kappa)) + sqrt(kappa)*sinh(sqrt(kappa))))` | `eta*(-1 + Cosh[Sqrt[kappa]])/(eta*kappa*Cosh[Sqrt[kappa]] + kappa^(3/2)*Sinh[Sqrt[kappa]])` |
| Delta_inf | `(eta*sinh(sqrt(kappa)) + sqrt(kappa)*(cosh(sqrt(kappa)) - 1))/(sqrt(kappa)*(eta*cosh(sqrt(kappa)) + sqrt(kappa)*sinh(sqrt(kappa))))` | `(-1 + Cosh[Sqrt[kappa]] + (eta*Sinh[Sqrt[kappa]])/Sqrt[kappa])/(eta*Cosh[Sqrt[kappa]] + Sqrt[kappa]*Sinh[Sqrt[kappa]])` |
| Delta_0(kappa->0+) | `1/2` | `1/2` |
| Delta_inf(kappa->0+) | `1` | `1` |
| Delta_0(eta->oo) | `(1 - 1/cosh(sqrt(kappa)))/kappa` | `(1 - Sech[Sqrt[kappa]])/kappa` |
| Delta_inf(eta->oo) | `tanh(sqrt(kappa))/sqrt(kappa)` | `Tanh[Sqrt[kappa]]/Sqrt[kappa]` |
| G_fail^inf | `Pe_req/(sqrt(kappa)*tanh(sqrt(kappa)))` | `peReq*Coth[Sqrt[kappa]]/Sqrt[kappa]` |
| G_suff^inf | `Pe_req*cosh(sqrt(kappa))/(cosh(sqrt(kappa)) - 1)` | `peReq*(1 + 1/(-1+Cosh[Sqrt[kappa]]))` |

All entries match algebraically and produce residual 0 against the target formulas in each engine's `expectZero` checks.

## Verdict justification

All fifteen sympy assertions and fifteen mathematica assertions pass with residual 0. The five "consistency of definitions" checks (A1-A5, B1-B5) are algebraic identities under the substitution `T_X = K_X L^2/kappa`, but that substitution is precisely the unit's claim 1, so the check is meaningful rather than vacuous. The eight limit-based checks (kappa->0+, eta->infinity, z->infinity) are non-trivial and exercise the exact regimes claimed in ledger items 3-5. Sign, factor-of-2, hyperbolic-identity, and branch attacks all reproduce the asserted target values. Both engines arrive at the same residuals via independent limit and simplification machinery, with displayed forms differing in non-trivial ways (e.g., `coth` vs `1/tanh`, `subs(alpha,z)` vs `kappa->z^2`). Outputs are fresh. No findings.

## Self-test notes

I checked: (1) factor-of-2 in `Delta_0(kappa->0+)` by expanding numerator and denominator to leading order in `alpha`, which gives 1/2 as asserted — no factor-of-2 trap; (2) parity/sign of `Delta_inf` and `Delta_0` for `kappa, eta > 0` — both strictly positive, no branch ambiguity; (3) the `eta->oo` limit by ratio of leading terms, which matches the target formulas; (4) the stiff-support `z->oo` limits using `tanh z -> 1` and `sech z -> 0` — both match `Pe_req`. No traps fired.
