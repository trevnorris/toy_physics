# Measurements behind S11b_sympy_x1_repair_directive.md (rule 2)

Every artifact claim in the directive is a `grep`/`sed` reading of the committed baseline (`864d6f41`) or the
spec (`1a2395a3`). Commands run from `research/pde_ledger_v3/`.

## Claim: the independence test is pointwise-polynomial and omits the divergence quotient
```
$ grep -n "def independent_columns" scripts/S11b_interface_coupling_law_sympy_audit.py
406:def independent_columns(expressions, variables):
# body (L406-416): builds sp.Poly(expr, *variables, domain="EX") per invariant, collects monomials,
# forms the coefficient matrix, returns tuple(matrix.rref()[1]).  => pointwise polynomial independence
# over the raw component symbols; NO total-divergence / integration-by-parts equivalence anywhere.
```

## Claim: §5 mandates equivalence modulo total divergences ("do not count both")
```
$ sed -n '286,287p' directives/S11b_SHARED_PHYSICS.md
- **Equivalence modulo total divergences** — two densities differing by a total in-plane divergence are the
  same term; ⛔ do not count both.
# (bullet sits inside "THE SYMMETRY GROUP, STATED IN FULL", L282+; the basis-count task is L279-281.)
```

## Claim: basis construction, emits, and coefficient assignment sites
```
$ grep -n "def construct_energy_basis\|ENERGY_BASIS\|full_energy\|independent_columns(" \
    scripts/S11b_interface_coupling_law_sympy_audit.py
463:def construct_energy_basis()
484:    pivots = independent_columns(basis, variables)
486:    carried_rank = len(independent_columns(carried, variables))
492:        trial_rank = len(independent_columns(trial, variables))
498:    full_energy = sp.expand( ... mu_R*basis[0]/2 + mu_S*basis[1]/2 + B_div*basis[2]/2 + ... )   # L498-510
618:    emit("ENERGY_BASIS", basis, ...)
619:    emit("ENERGY_BASIS_COUNT", sp.Integer(len(basis)), ...)
620:    emit("ENERGY_BASIS_OMISSIONS", ENERGY_DATA["omitted"], ...)
621:    emit("ENERGY_BASIS_INDEPENDENT_TERMS", ...)
624:    impermeable_reduced = ... ; 625: imp_pivots = independent_columns(impermeable_reduced, (eta, e_W))
632:    flux_reduced = ... ; 633: flux_pivots = independent_columns(flux_reduced, (eta, e_W))
```

## Claim: hand-enumerated coefficient rows that must be driven from the reduced basis
```
$ grep -n "mu_S\|MU_S_COEFFICIENT\|moduli = sp.Tuple" scripts/S11b_interface_coupling_law_sympy_audit.py
91:  mu_S = symbol("mu_S", "KNOB", "symmetric-traceless gradient coefficient", real=True)
500:      + mu_S * basis[1] / 2
1533: moduli = sp.Tuple(B_rho_3, C, k_W, kappa_W, B_div, mu_S, G_theta_u, G_W_u, kappa_theta, kappa_theta_W)
1808: ... B_div: "dim_b_div", mu_S: "dim_mu_s_coefficient",           # symbol_links
1827: mu_S: dimensions["MU_S_COEFFICIENT"][0], ...                    # symbol_dimensions
# plus the "MU_S_COEFFICIENT" dimension row in the `dimensions` dict (~L1795).
```

## Claim: the disagreement is already MEASURED and diagnosed (provenance only — values NOT restated)
The energy-basis disagreement between the two engines is recorded in the committed frozen comparator run
`scripts/out/S11b_cross_engine_comparison.out` (commit `17fe32c8`) and diagnosed against §5 in the committed
`steps/S11b_wl_engine_review_disposition.md`. ⛔ The withheld references — the resulting basis count, which
invariant is redundant, the identity among the invariants, and the fold coefficients — are deliberately NOT
restated in this twin or in the directive (rule 5, rule 12: blindness is enforced by ABSENCE, ⛔ not by a
prohibition-by-convention in a builder-readable file). The orchestrator holds those references in its own
independent derivation (scratchpad `x1_independent_basis_count.py`, run 2026-08-22: pointwise dimension vs
dimension modulo total in-plane divergence, difference one term) and diffs the engine's emitted objects
against them on our side; a mismatch is a finding, ⛔ not a build failure.
