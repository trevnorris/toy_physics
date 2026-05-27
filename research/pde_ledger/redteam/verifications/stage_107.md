---
unit_id: 107
batch: IV.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 107

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**

`/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage107_general_dtn_deformation_sympy_audit.py:60-67` now contains two new `expect_zero(...)` calls inserted immediately after the existing `print('Sigma4_evenmatch = ...')` line and before the `chi_even = sp.simplify(chiQ.subs(sol))` line:

```python
expect_zero(
    'Sigma2 exact formula',
    sol[Sigma2] - (-(3*S*beta**2 - 3*S + Sigma0)) / sp.Integer(9),
)
expect_zero(
    'Sigma4 exact formula',
    sol[Sigma4] - (-(3*S*beta**4 - 3*S + Sigma0)) / sp.Integer(27),
)
```

This matches the directive's "After" block character-for-character (same helper, same `sp.Integer(9)`/`sp.Integer(27)` usage, same paper-side RHS). No other line in the script was altered.

**Assessment:**

The new assertions are non-tautological: `sol[Sigma2]` and `sol[Sigma4]` come from `sp.solve` on the even-matching system at line 54, and the RHS is the paper's boxed closed form rewritten verbatim from `notes/stages/moving_throat_pde_stage107_general_dtn_deformation.md`. A wrong branch from `sp.solve` would produce a nonzero residual and trigger the `AssertionError` inside `expect_zero`. The SymPy exec log shows both new lines (`Sigma2 exact formula = 0`, `Sigma4 exact formula = 0`) appearing between `Sigma4_evenmatch = ...` and `chi_Q under canonical-even matching = ...`, exactly where the directive said they should land, and the pre-existing `normalized expansion direct-formula = 0` and `chi_Q - 3(S beta^5 + 9 Sigma5)/(3S - Sigma0) = 0` lines still appear. SymPy is now at engine parity with the Mathematica twin (which already had the corresponding `expectZero["Sigma2 exact formula", ...]` / `expectZero["Sigma4 exact formula", ...]` at `.wl:68-69`). No collateral edits.

## Exec log assessment

**SymPy:** exit=0. Notable lines:

```
Sigma2_evenmatch = -S*beta**2/3 + S/3 - Sigma0/9
Sigma4_evenmatch = -S*beta**4/9 + S/9 - Sigma0/27
Sigma2 exact formula = 0
Sigma4 exact formula = 0
chi_Q under canonical-even matching = 3*(S*beta**5 + 9*Sigma5)/(3*S - Sigma0)
chi_Q - 3(S beta^5 + 9 Sigma5)/(3S - Sigma0) = 0
```

**Mathematica:** exit=0. Notable lines:

```
PASS: normalized expansion direct-formula
Sigma2 exact formula = 0
PASS: Sigma2 exact formula
Sigma4 exact formula = 0
PASS: Sigma4 exact formula
chi_Q - 3(sNorm beta^5 + 9 sigma5)/(3 sNorm - sigma0) = 0
PASS: chi_Q - 3(sNorm beta^5 + 9 sigma5)/(3 sNorm - sigma0)
Stage 107 Mathematica audit passed.
```

**Output freshness:** both transcripts are newer than their scripts.
- SymPy script mtime 2026-05-27 15:08:45; SymPy output mtime 2026-05-27 15:18:03 (output newer).
- Mathematica script mtime 2026-05-27 15:08:50; Mathematica output mtime 2026-05-27 15:24:16 (output newer).

## Cluster C banner sweep

- SymPy `.py:25`: `banner('STAGE 107 — GENERAL ISOTROPIC DTN DEFORMATION ALGEBRA')` — STAGE 107 confirmed (no stale STAGE 90).
- Mathematica `.wl:26`: `banner["STAGE 107 — GENERAL ISOTROPIC DTN DEFORMATION ALGEBRA"];` — STAGE 107 confirmed (no stale STAGE 090).

Banner sweep clean in both engines.

## Material-change assessment

`material_change`: false.

The edit adds two new `expect_zero` assertions that anchor closed forms (`Sigma_2 = -(3 S beta^2 - 3 S + Sigma_0)/9`, `Sigma_4 = -(3 S beta^4 - 3 S + Sigma_0)/27`) that the Mathematica twin already asserted and that the existing SymPy `print` statements already exhibited. No derived value changes; only the verification surface area expands. Downstream Stage 108+ units rely on these closed forms, not on whether SymPy asserts them — engines were already in agreement on values prior to this fix.

## Side observations (non-blocking)

None.

## Verdict justification

Codex applied F1 exactly as directed: the two new `expect_zero` calls land at the right file:line range with the directive-prescribed `sp.Integer(9)` / `sp.Integer(27)` denominators and the paper-form RHS. Both engines exit 0 with all five assertions passing (SymPy: normalized expansion, Sigma2, Sigma4, chi_Q; Mathematica: same plus uniqueness). Outputs are fresher than scripts. Cluster C banner sweep already correct in both files. No collateral edits, no regressions, no tautologies. Verdict: `verified`.
