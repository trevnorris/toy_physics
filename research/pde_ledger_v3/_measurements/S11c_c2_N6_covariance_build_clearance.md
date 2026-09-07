# S11c-c2 N6 covariance instrument (R_cov) build — CLEARED (both build legs) (2026-09-06)

Instrument `scripts/S11c_c2_N6_covariance_sympy.py` (astra-built per CLEAR-TO-BUILD directive `123d9a18`; baseline
`326e6123`). Two build legs (astra-written → **fresh Claude agent + Grok**), identical prompt
`directives/_legs/S11c_c2_N6_covariance_build_review_prompt.md`. **Both legs: BUILD CLEAR** (convergent, independent
derivation + FORM ablation). Reports: fresh-Claude agent transcript + `scratchpad/grok_N6cov_buildreview.log`.

## Both legs, convergent
- **Non-circularity does real work.** Each leg re-derived the Φ prolongation from scratch; `prolonged_phi` matches
  (images_equal); `ms_pred = source_terms(μ_E.subs(Φ), V_E)` calls NO `material_pullback` (params `('mu_e','phi')`
  only). Ablation making the prediction circular collapses SOURCE_PREDICTED onto SOURCE_ACTUAL (byte-identical) ⇒
  `R_cov ≡ number(0)` 312/312 by construction. Shipped: distinct constructions that agree in VALUE.
- **F1 rank-2 prolongation load-bearing.** `μ_E=EL(E)` carries the six second jets `theta_didi`/`e_W_didi`; naive 0+1
  substitution ≠ full prolongation, and the gap IS exactly the rank-2 map (`μ_full − μ_01.subs(Φ_rank2)) == 0`).
  Truncating Φ to 0+1 moves `R_COV` (SOURCE_ACTUAL unchanged). `PHI_DOMAIN_CENSUS`: `uncovered=[]`, covered 11/11,
  max rank 2, before-substitution — a real computed comprehension, raise AFTER emit.
- **`R_COV_INCREMENT` correctly pinned.** = `closed_response(m_coeff, R_cov)` sig 6/9/12; ablation to `build_increment`
  injects the affine sig-0 `−C_M·p` (216→288 cols) — the shipped filter is load-bearing.
- **Both knives bite (shipped uncorrupted `κ_a=1, κ_j=0`).** Φ-coefficient (`κ_a=2`, prediction held at 1): `R_COV`/
  `R_COV_CONTROL_DELTA` move, SOURCE_PREDICTED unchanged — bites though the `a_ρ+h_α` truth table is unchanged (what
  the truth table cannot see). Junk (`κ_j=1`, `J_μ·e_W`, RHO4): `R_COV`/`SOURCE_CONTROL_DELTA` move (4 cols) where
  baseline `R_cov=0`. Non-vacuous.
- **PIT sound + 0 asserts.** One joint `n.pit`; shared samples; 3 primes with regeneration; joint singular rejection;
  on-shell `q` from `k`; honest FN bound `δ=min(1, family·max(per_prime))=2.597e-22`; residual-zero never an exit.
  `V_E` SHA = `V_M` SHA (velocity channel null, μ isolated). No `assert` on a residual.
- Non-blocking (agent): the prediction keeps `V_E` (not `Φ(V_E)`) — a spec/premise question, null here since
  `V_E≡V_M` (see the RESOLVED caveats); cosmetic unconditional `cache_clear()` at `κ_j=0`; `R_COV_INCREMENT`
  double-protected (closed_response + filter). None change what is computed.

## The R_cov disposition (mechanical tally, `_measurements/S11c_c2_N6_reconcile_disposition_tally.py`-style over the 6 .out)
| object | LAB_HELD.RHO4 | LAB_HELD.RHOBR | MAT_ADV.RHO4 | MAT_ADV.RHOBR | a_knife | junk_knife |
|---|---|---|---|---|---|---|
| `SOURCE_ACTUAL` | 172/312 | 146/312 | 178/336 | 152/336 | (live) | (live) |
| `SOURCE_PREDICTED` | 172/312 | 146/312 | 178/336 | 152/336 | (held) | (held) |
| **`R_COV`** | **0** | **0** | **0** | **0** | **84** | **4** |
| `R_COV_INCREMENT` | 0 | 0 | 0 | 0 | 18 | 9 |
| `PHI census uncovered` | [] | [] | [] | [] | [] | [] |

`R_cov = 0` (no nonzero found at δ≈2.6e-22) in all 4 baseline cases with the knives biting (84, 4) ⇒ the material
construction **faithfully implements the declared field map Φ**; the source transforms covariantly. Disposition
adjudicated in `_measurements/S11c_c2_N6_RESOLVED.md`.
