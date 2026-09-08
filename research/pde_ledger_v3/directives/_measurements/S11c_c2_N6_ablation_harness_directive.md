# Mechanical grounding — S11c-c2 N6 ablation-harness directive site claims

The directive `S11c_c2_N6_ablation_harness_directive.md` pins each knife to a unique construction fragment. Below is
the orchestrator's mechanical fact-lookup (fixed-string occurrence counts) confirming each pinned old-fragment occurs
exactly once in its engine — corroborating both round-2 review legs (Grok's `/tmp/n6_r2_review/knife_checks.py`
reported the same counts; the fresh Claude leg grep-confirmed each individually). Mechanical literal-match counting
only; no derived predicate.

Command (run 2026-09-08, working dir `research/pde_ledger_v3`):

```sh
grep -Fc "materialNormalKnife = 0;" mathematica/S11c_c2_N6_mathematica_audit.wl
grep -Fc "actualJunkCoefficient = 0;" mathematica/S11c_c2_N6_mathematica_audit.wl
grep -Fc "ACTUAL_JUNK = sp.Integer(0)" scripts/S11c_c2_N6_covariance_sympy.py
grep -Fc "rows=flatten({'U':folded['U'],'E_W':folded['E_W'],'THETA':mass-correction})" scripts/S11c_c2_N6_diagnostic_sympy.py
grep -Fc "slots=tuple(inputs.a(prefix+label) for label in ('plus','minus') for prefix in ('delta_p_','d_w_delta_p_'))" scripts/S11c_c2_N6_diagnostic_sympy.py
grep -Fc "rows, _, provenance = n.face_factory(a, b, inputs, alpha, rho, 'MATERIAL', mu_slot)" scripts/S11c_c2_N6_reconcile_sympy.py
grep -Fc "source = closed_response(comp, inputs, m_coeff, ds, kernels)" scripts/S11c_c2_N6_reconcile_sympy.py
```

Literal stdout:

```text
WL materialNormalKnife=0 : 1
WL actualJunkCoefficient=0 : 1
COV ACTUAL_JUNK=0 : 1
DIAG flatten E_W row : 1
DIAG slots tuple : 1
REC face_factory MATERIAL : 1
REC source=closed_response : 1
```

Full per-knife site + API + parse + cone verification (both legs): review record
`_measurements/S11c_c2_N6_ablation_harness_directive_review.md`; Grok leg + its check-script stdout
`_measurements/S11c_c2_N6_ablation_harness_directive_review_r2_grok.md`.
