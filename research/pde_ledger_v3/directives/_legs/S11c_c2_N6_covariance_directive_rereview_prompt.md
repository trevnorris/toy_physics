# Round-2 SCOPED re-review — the S11c-c2 N6 covariance (R_cov) directive, after folding 3 findings

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_covariance_directive.md` (v2). Round-1 (2 legs)
reported 3 folds; now applied. **Tight confirmation**, ⛔ not a fresh full review: verify each fold landed correctly +
faithfully against the sources, and scan ONLY for a regression. ⛔ No CAS, ⛔ no build, ⛔ no fictional-script ablation.
Working dir `/var/projects/toy_physics`. Adjudication of round-1:
`_measurements/S11c_c2_N6_covariance_directive_review_adjudication.md`.

## The 3 folds to confirm (re-verify against the sources, not just that text is present)
1. **F1 — Φ prolonged to the jet order in `μ_E`.** Confirm the directive now requires `μ_E.subs(Φ)` to prolong the
   declared field map through EVERY θ/`e_W` jet symbol occurring in imported `μ_E` (rank-2: `θ_I↦D_I(θ+a_ρ)`,
   `e_{W,I}↦D_I(e_W+h_α)`), generated via the LIVE `b.total_derivative(...,background_depth=3)`+`DERIVATIVE_MAP` chain,
   `simultaneous=True`, with a pre-substitution domain-coverage census — and explicitly forbids reusing the
   energy-level 0+1 `frozen_relations`/`material_pullback` dict. Grounds: `μ_E=EL(E)` writes second jets (diagnostic
   `constitutive` :345-347; `DERIVATIVE_MAP` `S11c_b_brane_operator_sympy_audit.py:733`; `theta_didi`/`e_W_didi`
   present at `S11c_b_exports.py:5629`). Is the fix complete (does it cover the ACTUAL jet ranks present in `μ_E`, and
   is the census the right guard against a silent uncovered atom)?
2. **F2a — `R_COV_INCREMENT` pinned to `closed_response`.** Confirm it is now the reconcile `closed_response` on
   `(m_coeff, R_cov)` (sig 6/9/12 only, `SOURCE_CHANNEL` keying), ⛔ not `build_increment`/`I(C_M,R_cov)` (which
   re-adds `−C_M·p` → certified-nonzero at `R_cov=0`).
3. **F2b — θ-independent-junk knife made concrete.** Confirm it is now an actual-only, wave-LINEAR perturbation
   (`κ_j·J_μ·e_W`, fresh nonzero `J_μ`, `DIMENSION_SCHEMA[J_μ]=(-1,-2,1)`, baseline `κ_j=0`, prediction unchanged),
   satisfying `wave_terms`' one-wave-per-term requirement (diagnostic :363), ⛔ not a literal θ-independent constant.

## Regression scan
Did the folds contradict any round-1 "held sound" content — non-circularity (imported `μ_E`+supplied `Φ`, forbidding
`material_pullback`/`ms`/`ΔS`), the `a_ρ`/`h_α` formulas, the object `R_cov=ms−ms_pred`, the μ-channel isolation, the
one-sided Φ-coefficient knife, the withheld strict-vs-covariance interpretation, no expected value / no residual-zero
exit, the three script clauses + corollaries?

## Physics filter
Report a finding only if it catches a way the sufficient test could be circular, vacuous, leak-prone, intractable, or
mislabel covariance — or a NEW defect/regression the folds introduced.

## Output
End with **DIRECTIVE SOUND — CLEAR TO BUILD** or the exact remaining fold (`file:line` + minimal fix). Brief.
