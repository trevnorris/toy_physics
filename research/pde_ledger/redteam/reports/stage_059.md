---
unit_id: 059
batch: III.2
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-26
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage059_operator_branch_residual_bounds.md
  paper_appendix: present
---

# Audit unit 059 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_059.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage059_operator_branch_residual_bounds.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 96 references this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.txt`

## What the paper claims

The paper card's `\stagefield{Output}` reads: "The exact gain thresholds \eqref{eq:app-stage059-Xi-thresholds}," referring to the boxed equations `Xi_fail = Pe_req / Delta_inf(kappa,eta)` and `Xi_suff = Pe_req / Delta_0(kappa,eta)`, together with the verdict rule `Xi < Xi_fail => fail`, `Xi > Xi_suff => succeed`. The notes elaborate seven deliverables: (1) the operator-selected physical ratio `zeta_phys(Xi,eta;kappa) = Omega_(Pe_*)^2 (kappa+pi^2/4)/(kappa+y^2)`, (2) the exact bracket `zeta_- <= zeta_phys <= zeta_+` driven by Stage-41's `Xi Delta_0 <= Pe_* <= Xi Delta_inf`, (3) the residual envelope `R_- <= R_phys <= R_+`, (4) the gate implications (success when `zeta_- >= zeta_req`, no-go when `zeta_+ < zeta_req`), (5) the `Pe_req` definition `Omega_(Pe_req)^2 = zeta_req (kappa+y^2)/(kappa+pi^2/4)`, (6) the threshold definitions `Xi_fail`/`Xi_suff` with the explicit ordering `Xi_fail <= Xi_suff` (because `Delta_inf >= Delta_0 > 0`), and (7) the weak-coupling expansion `zeta_phys = A_K [1 + ((4-pi)/pi) Xi Delta_0 + O(Xi^2)]`. The appendix row 96 simply restates the threshold deliverable.

## What the script claims to verify

The SymPy script (docstring lines 1-10) lists four checks: (a) symbolic forms of `zeta_-`, `zeta_+`, (b) symbolic forms of `R_-`, `R_+`, (c) symbolic forms of `Xi_fail`, `Xi_suff`, and (d) a weak-coupling expansion of `Omega_Pe^2` whose linear coefficient is `(4-pi)/pi`. The script then adds an additional numerical probe (lines 87-118) that solves the operator-equation `A_K Omega_Pe(Xi*Delta)^2 = zeta_req` at both the `Delta_inf` and `Delta_0` endpoints, asserting that the resulting `Xi*Delta = Pe_star` where `Pe_star` is defined by `A_K Omega_Pe(Pe_star)^2 = zeta_req` (i.e., that the threshold definitions saturate the bracket-boundary equations). The Mathematica script mirrors these checks. Neither script asserts the threshold ordering `Xi_fail <= Xi_suff`, even though both define `DeltaInf_ordered = Delta0 + delta_gap` and `Xi_fail_ordered`, `Xi_suff_ordered` precisely to set that comparison up.

## Paper ↔ script cross-check

| Paper deliverable | Script-side coverage | Status |
|---|---|---|
| `zeta_-`, `zeta_+` symbolic forms (notes 2) | SymPy L54-57; WL L55-56 (print only) | match |
| `R_-`, `R_+` definitions (notes 3) | SymPy L59-62; WL L57-58 (print only) | match |
| `Xi_fail`, `Xi_suff` definitions (paper Output) | SymPy L64-67 (symbolic) + L96-118 (numerical saturation); WL L59-67 + L87-106 | match |
| `Pe_req` defining equation `Omega_(Pe_req)^2 = zeta_req(kappa+y^2)/(kappa+pi^2/4)` (notes 5) | implicit in the FindRoot for `Pe_star` (SymPy L96-101; WL L87-91) | match |
| Threshold ordering `Xi_fail <= Xi_suff` (notes 4 closing) | dead scaffolding only (SymPy L68-72; WL L68-72) | missing |
| Bracket-implication chain (success/no-go gates, notes 3-4) | not directly asserted; relies on upstream monotonicity | partial (carried from Stage 39/40/41) |
| Weak-coupling slope `(4-pi)/pi` (notes 5) | SymPy L82-85; WL L78,111 | match |

The dominant pattern is `match` with one `missing` and one `partial` justified by upstream carry-forward. `paper_alignment: aligned` — the script's load-bearing checks faithfully exercise the paper's `\stagefield{Output}` and the notes' substantive deliverables.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 82-85 | `expect_zero(coeff(Omega_Pe^2, Pe, 1) - (4-pi)/pi)` | weak-coupling slope (notes 5) | yes |
| A2 | sympy | 117 | `expect_close(Xi_fail_solved * DeltaInf_probe, Pe_star)` | `Xi_fail = Pe_req/Delta_inf` saturates upper bracket | yes |
| A3 | sympy | 118 | `expect_close(Xi_suff_solved * Delta0_probe, Pe_star)` | `Xi_suff = Pe_req/Delta_0` saturates lower bracket | yes |
| A4 | mathematica | 105 | `expectApprox[xiFailSolved*deltaInfProbe, peStar, 1e-40]` | mirror of A2 | yes |
| A5 | mathematica | 106 | `expectApprox[xiSuffSolved*delta0Probe, peStar, 1e-40]` | mirror of A3 | yes |
| A6 | mathematica | 111 | `expectZero[omegaSqLinearCoeff - (4-Pi)/Pi]` | mirror of A1 | yes |

Note: the symbolic `zeta_lo`, `zeta_hi`, `R_lo`, `R_hi`, `Xi_fail`, `Xi_suff`, `Xi_fail_ordered`, `Xi_suff_ordered`, `zeta_req_branch` definitions (SymPy L54-72, WL L55-72) carry no `assert` — they are printed but not exercised. The `_ordered` scaffolding was clearly set up to support an ordering assertion that was never written. This is the F1 finding.

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py:68-72`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl:68-72`

**What's wrong:**
Both scripts define scaffolding to compare `Xi_fail` and `Xi_suff` under a positivity-ordered `Delta_inf = Delta_0 + delta_gap` (delta_gap > 0):

SymPy L68-72:
```
delta_gap = sp.symbols("delta_gap", positive=True, real=True)
DeltaInf_ordered = Delta0 + delta_gap
Xi_fail_ordered = sp.simplify(Pe_req / DeltaInf_ordered)
Xi_suff_ordered = sp.simplify(Pe_req / Delta0)
zeta_req_branch = sp.simplify(A_K * Omega(Pe_req) ** 2)
```

Mathematica L68-72:
```
Clear[deltaGap];
$Assumptions = $Assumptions && deltaGap > 0;
deltaInfOrdered = FullSimplify[delta0 + deltaGap, ...];
xiFailOrdered = FullSimplify[peReq/deltaInfOrdered, ...];
xiSuffOrdered = FullSimplify[peReq/delta0, ...];
```

Neither script actually asserts `Xi_fail_ordered <= Xi_suff_ordered` (or equivalently `Xi_suff_ordered - Xi_fail_ordered > 0`). The variables are computed and then abandoned. The paper notes section 4 ("Since `Delta_inf >= Delta_0 > 0`, one has `Xi_fail <= Xi_suff`") states this as one of the enumerated deliverables and as the structural condition that makes the three-zone `[Xi_fail, Xi_suff]` window well-defined. The SymPy `zeta_req_branch = A_K * Omega(Pe_req)^2` is also unused (no `expect_zero` references it).

**Why this matters:**
The threshold ordering is the structural invariant that makes the paper's verdict rule (`Xi < Xi_fail => fail`, `Xi > Xi_suff => succeed`, with a profile-sensitive band between them) coherent. If a future refactor swapped `Delta_0` and `Delta_inf` roles or introduced a sign error in `Pe_req`, nothing in this stage's verification would catch it because the relevant scaffolding sits idle. The current state is also misleading on inspection — a reader sees the `_ordered` variables and reasonably expects an ordering check.

**Required change:**
Add a positivity assertion `expect_positive("Xi_suff - Xi_fail (ordered)", Xi_suff_ordered - Xi_fail_ordered)` in the SymPy script immediately after line 72, and a corresponding `expectPositive["Xi_suff - Xi_fail (ordered)", xiSuffOrdered - xiFailOrdered]` in the Mathematica script immediately after line 72. Both `expect_positive`/`expectPositive` helpers are already defined in the respective scripts and operate under the positive-real `$Assumptions`/symbol declarations that cover `Pe_req`, `Delta0`, `delta_gap`. The unused `zeta_req_branch` line is collateral and can be left alone.

**Verification:**
After Codex applies, the SymPy transcript should show a new line of the form `Xi_suff - Xi_fail (ordered) = Pe_req*delta_gap/(Delta0*(Delta0 + delta_gap))` (or equivalent factored form) and the Mathematica transcript should show a `PASS: Xi_suff - Xi_fail (ordered)` line. The verifier will run `redteam exec-sympy 059` and `redteam exec-mathematica 059` and confirm both scripts still exit 0.

## Independent-derivation check (Mathematica)

Both scripts share the same algebra for the symbolic threshold definitions (`zeta_lo`, `zeta_hi`, `R_lo`, `R_hi`, `Xi_fail`, `Xi_suff`) and the same probe values (kappa=2, y=1, Delta0=3/5, delta_gap=2/5, zeta_req=2/5) in the numerical section. The symbolic definitions are not "derivations" in any meaningful sense — they are restatements of the paper's symbolic forms; both engines necessarily produce parallel printouts. The substantive checks use independent methods:

- SymPy line 80: `Omega_sq_series = sp.series(Omega_Pe**2, Pe, 0, 2).removeO()` then `.coeff(Pe, 1)`
- Mathematica line 78: `Limit[D[omegaPe^2, pe], pe -> 0]`

These are genuinely distinct algorithms (Taylor series coefficient extraction vs. analytic derivative-at-zero limit). The FindRoot/nsolve probes use the engines' native rootfinders. The Mathematica `Limit::alimv` warning indicates the positivity assumption was discarded by Limit; the numerical result `-1 + 4/Pi` is still correct because `Omega_Pe^2` is analytic at 0. Verdict: not a transliteration — parallel structure is forced by the stage's nature, and the load-bearing checks use independent engine-native methods.

## Engine cross-check

Both transcripts agree:

- SymPy: `Omega_Pe^2 small-Pe series = Pe*(-1 + 4/pi) + 1`, linear coefficient verified zero.
- Mathematica: `Omega_Pe^2 linear coefficient = -1 + 4/Pi`, verified equal to `(4-Pi)/Pi`.
- SymPy: `Xi_fail*DeltaInf saturates at Pe_star diff = 0`; `Xi_suff*Delta0 saturates at Pe_star diff = 7.24...e-71` (well below 1e-40 tolerance).
- Mathematica: both `Xi_fail*DeltaInf` and `Xi_suff*Delta0` saturate at `0``49.0959...` (a 50-digit zero representation), passing the 1e-40 tolerance.
- The symbolic `zeta_-`, `zeta_+`, `R_-`, `R_+`, `Xi_fail`, `Xi_suff` printouts agree up to obvious naming conventions (`Xi`/`capitalXi`, `Delta0`/`delta0`).

Engines agree.

## Verdict justification

The script's load-bearing assertions (weak-coupling slope `(4-pi)/pi`, threshold-saturation `Xi_fail*Delta_inf = Pe_star` and `Xi_suff*Delta_0 = Pe_star`) faithfully exercise the paper's `\stagefield{Output}` thresholds and the notes' Pe_req defining equation. Both engines pass non-trivially. Paper alignment is solid — the script even uses a `zeta_req` chosen independently of `Omega_Pe` to avoid construction-induced tautology (SymPy L95, WL L86). The one substantive gap is the dead `_ordered` scaffolding intended to verify `Xi_fail <= Xi_suff`, which the paper notes states explicitly. Attacks I tried: (a) could the FindRoot have multiple roots? — no, `Omega_Pe^2` is strictly monotone on `Pe > 0` by Stage 39, so the root used in the initial guess is the unique constructive one; (b) sign error in the bracket direction? — zeta_- uses Delta_0, zeta_+ uses Delta_inf; for `Delta_inf >= Delta_0`, monotonicity of `Omega_Pe` gives `Omega(Xi*Delta_inf) >= Omega(Xi*Delta_0)`, so `zeta_+ >= zeta_-` — consistent; (c) the docstring/banner say "Stage 42" while the file and paper say "Stage 059" — cosmetic-only, not load-bearing. The findings_count is 1 (insufficient_verification, medium).

## Self-test notes

Variable-independence: the `expect_positive` I prescribe operates on `Xi_suff_ordered - Xi_fail_ordered = Pe_req*delta_gap/(Delta0*(Delta0+delta_gap))` — depends on three declared positive symbols, manifestly positive under existing assumptions. Trivial-case check: substituting `delta_gap = 0` gives `Xi_suff - Xi_fail = 0` (boundary excluded by `delta_gap > 0`); any `delta_gap > 0` yields a strictly positive value. Paper round-trip: the prescribed assertion matches the inequality stated in notes section 4 ("Since `Delta_inf >= Delta_0 > 0`, one has `Xi_fail <= Xi_suff`") and is consistent with the paper card's three-zone verdict structure; the fix introduces no new paper_misalignment.
