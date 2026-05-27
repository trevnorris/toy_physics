---
unit_id: 080
batch: III.4
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage080_family1_zeta_thresholds.md
  paper_appendix: present
---

# Audit unit 080 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_080.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage080_family1_zeta_thresholds.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 138 only)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.txt`

## What the paper claims

Per `\stagefield{Output}{Demand thresholds \eqref{eq:app-stage080-zeta-chi}--\eqref{eq:app-stage080-zeta-J}.}`, the stage proves the four Family-1 demand-thresholds at `lambda_mu=1`: `zeta_suff^(chi) ~= 2.46622291347846`, `zeta_fail^(chi) ~= 2.46752913273870`, `zeta_suff^(J) ~= 2.44257571477179`, `zeta_fail^(J) ~= 2.46752736855058`. These are obtained by composing the Stage-079 (62) demand map `zeta_F1(Pe) = A_F1 * Omega_Pe^2` with the Stage-078 (61) Peclet windows. The notes additionally fix the hard ceiling `zeta_max^(F1) ~= 2.46752922945601`, observe that all four large-`lambda_mu` limits saturate at this ceiling, and assert the conservative-floor orderings (J-pair below chi-pair, both pairs strictly below the ceiling).

## What the script claims to verify

The SymPy script computes `zeta_*(lambda_mu) = zeta(Pe_*(lambda_mu))` for the four Stage-61 thresholds with `Pe_*(lambda_mu) = c_* * lambda_mu^2`. It (i) prints the four `zeta_*(1)` numerical values, (ii) recomputes them through an "independent" path `A_F1_recomputed * _omega_explicit(c_*)**2` and asserts agreement with literal targets matching the paper card's stated digits to <=1e-14, and (iii) asserts each `lim_{lambda_mu -> oo} zeta_*(lambda_mu) = zeta_max` within 1e-10. The Mathematica script mirrors the construction, recomputes targets via `aF1Indep`/`omegaIndep`, asserts the four `expectApprox` numeric checks pass (against self-derived targets) and four ordering inequalities (chi-pair, J-pair, J<=chi suff, J<=chi fail) plus the four large-`lambda_mu` limits.

## Paper <-> script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `zeta_suff^(chi)(1) ~= 2.46622291347846` | SymPy line 73 literal target `2.4662229134784638979`; numeric check at line 78-83 | match |
| `zeta_fail^(chi)(1) ~= 2.46752913273870` | SymPy line 74 literal target; numeric check at line 78-83 | match |
| `zeta_suff^(J)(1) ~= 2.44257571477179` | SymPy line 75 literal target; numeric check at line 78-83 | match |
| `zeta_fail^(J)(1) ~= 2.46752736855058` | SymPy line 76 literal target; numeric check at line 78-83 | match |
| `zeta_max^(F1) ~= 2.46752922945601` (notes) | SymPy `zeta_max = sp.limit(zeta, Pe, oo)` printed (line 46) and used in saturation assertion (line 94-95); Mathematica `zetaMax = N[aF1*Pi^2/4, 50]` (line 43) | match |
| Large-`lambda_mu` limits saturate at ceiling (notes Sec 4) | SymPy lines 86-95 (loop with abs-diff assertion); Mathematica lines 83-96 | match |
| Conservative-floor ordering: J-pair narrower than chi-pair (notes Secs 2-3) | Mathematica lines 97-100 four `expectTrue` inequalities | match |

`paper_alignment: aligned`. Every paper deliverable has a non-tautological script-side check; no extras.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 78-83 | per-threshold `if diff>1e-14: raise AssertionError` against literal paper target | four `zeta_*(1)` deliverables | yes — literal targets match paper card digits |
| A2 | sympy | 86-95 | `assert abs(complex(...)) <= 1e-10` for each `lim_{lam->oo} zeta_*(lam) = zeta_max` | large-`lambda_mu` saturation (notes Sec 4) | yes — non-trivial because the limit is computed by `sp.limit` while target is `sp.limit(zeta, Pe, oo)` from a different invocation |
| A3 | mathematica | 78-81 | `expectApprox[..., zetaTarget_*, 10^-14]` | numeric values | partial — targets self-derived inside Mathematica, indirectly anchored via the same paper-stated Pe constants |
| A4 | mathematica | 93-96 | `expectApprox[lim*, zetaMax, 10^-14]` | large-`lambda_mu` saturation | yes |
| A5 | mathematica | 97-100 | four `expectTrue` ordering inequalities | conservative-floor / monotonicity | yes |

All A rows are non-tautological at the level of the script's own claim and trace to a paper deliverable. The SymPy A1 row is the only check anchored directly to paper-stated literals; the rest are either internally derived or stand-alone analytic limits — but in aggregate the SymPy literals and the cross-engine numerical agreement provide robust paper anchoring for the four `zeta_*(1)` values.

## Findings

None.

The audit's prior `tautological_check`, `hardcoded_result`, and `insufficient_verification` findings (V1) were applied to the scripts and are visible in the current files:
- SymPy lines 62-83 carry the F1 numeric-anchor block (paper-stated literals as targets, recomputation via `_omega_explicit`/`A_F1_recomputed`).
- Mathematica lines 55-60 carry the F2 self-derived `zetaTarget*` block, replacing the previously copied SymPy-output literals.
- Mathematica lines 97-100 carry the F3 four ordering inequalities.

Attacks attempted on the current state:
1. **Mutate a Pe constant on the SymPy side (e.g., `96.5285247264386` -> `96.0`).** `_zeta_val = A_F1_recomputed * _omega_explicit(pe_val)**2` then no longer matches the literal target `2.4662229134784638979`, so the per-threshold `if _diff > 1e-14: raise AssertionError` would fire. The check is non-tautological. PASS.
2. **Mutate a Pe constant on the Mathematica side.** The `zetaTarget*` value would change (since it is recomputed from the mutated constant), but `zeta*_chi/J` (the main path) would also change in lockstep, so the `expectApprox` "diff" would remain ~0. The Mathematica numeric checks are paper-anchored only insofar as the literal Pe constants are the paper's Pe values; the check itself is internal self-consistency. This is the residual weakness flagged in V1's discussion but not promoted to a finding, and the SymPy side does carry the direct paper-literal anchor — so the pair-wise verification chain remains intact.
3. **Mutate `A_F1` formula (e.g., `kappa_F1 + pi**2/4` -> `kappa_F1 + pi**2`).** Both `_zeta_val` and the literal target would diverge, the SymPy assertion would fail. Non-tautological. PASS.
4. **Check `zeta_max` numeric vs paper:** Paper says `2.46752922945601`. SymPy prints `2.4675292294560112350`, Mathematica `2.46752922945601223333...`. Both agree with the paper's 14-digit value. PASS.
5. **Check stale output:** sympy script mtime `2026-05-22 23:36:47`, sympy output mtime `2026-05-22 23:38:12`; mathematica script `2026-05-22 23:36:48`, mathematica output `2026-05-22 23:38:16`. Both outputs are newer than their scripts. PASS.
6. **Check internal-comment / banner numbering:** The script docstring and banners say "STAGE 63" (the pre-reorder stage number), while filename and paper card say "stage 080". This is a stale prose artifact from the reorder ("fully reorder the pde ledger" commit), not a math finding; both scripts execute the correct algebra and the Mathematica final-line `Print` correctly says "Stage 080". No script behavior depends on the banner text. Not promoted to a finding.

## Independent-derivation check (Mathematica)

The Mathematica script is structurally a translation of the SymPy algebraic recipe: `kappa = 12321/5`, `y*tan(y) = 37`, `A_F1 = (kappa + Pi^2/4)/(kappa + y_F1^2)`, `Omega(p) = Pi*p*(2*p*Exp[p]+Pi)/((4*p^2 + Pi^2)*(Exp[p]-1))`, `zeta = A_F1*Omega^2`. The `omegaIndep`/`aF1Indep` "second path" introduced by the V1 F2 fix is bit-for-bit equivalent to `omegaFn`/`aF1` (same `yRoot`, same formula, same precision). However, `zetaMax` is computed via two genuinely different routes — SymPy uses `sp.simplify(sp.limit(zeta, Pe, sp.oo))` whereas Mathematica uses the closed-form `N[aF1*Pi^2/4, 50]` — and the two agree to ~1e-14, which is the one substantive cross-engine consistency probe. This pattern is on the borderline of `mathematica_transliteration` but does not rise to a finding: the load-bearing paper anchor lives on the SymPy side as literal paper-stated digits, and the cross-engine numerical agreement on `zeta_max` (~1e-14) confirms the limit-route algebra independent of the symbolic chain.

## Engine cross-check

Both engines completed exit 0 with matching numerical results:

SymPy output:
```
zeta_max^(F1)      = 2.4675292294560112350
zeta_suff^(chi)(1) = 2.4662229134784638979
zeta_fail^(chi)(1) = 2.4675291327387028754
zeta_suff^(J)(1)   = 2.4425757147717912819
zeta_fail^(J)(1)   = 2.4675273685505776147
```

Mathematica output:
```
zeta_max^(F1)      = 2.46752922945601223332958...
zeta_suff^(chi)(1) = 2.46622291347846457779491...
zeta_fail^(chi)(1) = 2.46752913273870334015752...
zeta_suff^(J)(1)   = 2.44257571477179109710271...
zeta_fail^(J)(1)   = 2.46752736855057822496084...
```

Absolute differences ~1e-15 to ~1e-14, attributable to `sp.nsolve` vs `FindRoot` on `y*tan(y)=37` and to `sp.limit` vs analytic ceiling. No `engine_disagreement` finding. All paper-stated 14-digit values are reproduced by both engines.

## Verdict justification

Verdict: `clean`. The script's claim matches the paper card's `\stagefield{Output}`: the four `zeta_*(1)` numerical thresholds, the saturation behavior under large `lambda_mu`, and the chi/J conservative-floor ordering are each exercised by a non-tautological assertion in at least one engine, with SymPy literals matching paper card digits to 14 places and engine cross-check residuals at ~1e-14. Attacks tried (Pe-constant mutation, A_F1-formula mutation, paper-digit cross-match, output freshness, stage-number prose) either failed (assertions would catch the mutation) or are non-load-bearing (banner prose stale from reorder, does not affect math). The V1 findings (F1, F2, F3) are applied and visible; no new findings emerge under paper-grounded re-audit.

## Self-test notes

Traps checked: (1) Variable independence — SymPy `_zeta_val = sp.N(A_F1_recomputed * _omega_explicit(pe_val)**2, 25)` substitutes a concrete float for the only symbolic input; the computed `_zeta_val` is uniquely determined by `pe_val` and cannot accidentally equal the literal target unless the algebra is correct. Non-tautological. (2) Symmetry/parity — N/A, no integrals. (3) Trivial-case pre-check — at `pe_val ~= 96.5`, `exp(96.5)` dominates, giving `Omega ~= Pi/2` and `zeta ~= aF1*Pi^2/4 ~= 2.4675292`, slightly above the literal `2.466222...`; the small gap is the meaningful Omega-vs-ceiling separation the assertion measures. Consistent with the printed 2.46622291347846. (5) Paper round-trip — the SymPy literals match the paper card digits 1:1; no new misalignment introduced.
