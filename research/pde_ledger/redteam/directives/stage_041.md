---
unit_id: 041
batch: III.1
created_at: 2026-05-22T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-22T18:33:33Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 041

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage041_rank2_support_sympy_audit.py:108-124`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage041_rank2_support_mathematica_audit.wl:88-112`

**Issue:**
Section 24.4 in both engines builds `n_src` from raw expressions instead of substituting into the general formula `n_expected` derived in section 24.1. The "expected" form `n_src_expected` is then constructed as the same expression with `lam0` left symbolic — character-for-character identical after evaluation of the `q2_src, r2_src, qr_diff_sq` placeholders. The assertion `expect_zero("source-tied formula", n_src - n_src_expected)` is therefore guaranteed to pass regardless of the algebra and does not verify that the source-tied formula follows from the general result under the substitution `q = t R_U, r = t, t^2 = lambda0`. The same flaw is repeated in the Mathematica script (`nSrc` vs `nSrcExpected`).

**Required change:**

### Edit 1 — SymPy script

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage041_rank2_support_sympy_audit.py`, replace lines 108–116 (the block defining `q2_src, r2_src, qr_src, qr_diff_sq` and constructing `n_src`) so that `n_src` is built from `n_expected` (the general formula derived on line 78) via the physical substitution `q -> t*R_U, r -> t, t^2 -> lam0`.

**Before** (current lines 108–116):

```python
# The physical source-tied specialization uses q = t R_U, r = t with t^2 = lam0.
# Because only q^2, r^2, q r, and (q-r)^2 appear, we can substitute directly.
q2_src = lam0 * R_U**2
r2_src = lam0
qr_src = lam0 * R_U
qr_diff_sq = lam0 * (R_U - 1)**2

n_src = sp.simplify(
    (xi * (delta + xi) - m * (delta + (1 + q2_src) * xi))
    / (delta + (1 + r2_src) * xi - m * qr_diff_sq)
)
```

**After**:

```python
# The physical source-tied specialization uses q = t R_U, r = t with t^2 = lam0.
# Derive n_src by substituting into the general n_expected from section 24.1,
# not by re-stating the substituted formula.
t = sp.symbols("t", real=True)
n_src = sp.simplify(
    n_expected.subs({q: t * R_U, r: t})
)
n_src = sp.simplify(sp.expand(n_src).subs(t**2, lam0))
```

Leave lines 118–124 unchanged (the `n_src_expected` block and the `expect_zero("source-tied formula", n_src - n_src_expected)` call). After this edit, that `expect_zero` call becomes a genuine cross-check: sympy must verify that the substituted general formula reduces to the hand-written source-tied form.

Leave lines 126–143 unchanged (the regularity/positivity thresholds and the `dn_dm_src` test). Those operate on `n_src_expected` and remain non-tautological derivative checks.

### Edit 2 — Mathematica script

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage041_rank2_support_mathematica_audit.wl`, replace lines 88–96 so that `nSrc` is derived from `nExpected` by substitution.

**Before** (current lines 88–96):

```
q2Src = lambda0 rU^2;
r2Src = lambda0;
qrSrc = lambda0 rU;
qrDiffSq = lambda0 (rU - 1)^2;

nSrc = FullSimplify[
  (xi (delta + xi) - m (delta + (1 + q2Src) xi))/(delta + (1 + r2Src) xi - m qrDiffSq),
  Assumptions -> $Assumptions
];
```

**After**:

```
(* Derive nSrc from nExpected by substituting the source-tied invariants.
   The physical specialization is q = t rU, r = t with t^2 = lambda0. *)
nSrc = FullSimplify[
  nExpected /. {q -> t rU, r -> t},
  Assumptions -> $Assumptions
];
nSrc = FullSimplify[
  Expand[nSrc] /. t^2 -> lambda0,
  Assumptions -> $Assumptions
];
```

Also update the `Clear[...]` statement on line 35 from

```
Clear[a0, delta, xi, m, n, q, r, rU];
```

to

```
Clear[a0, delta, xi, m, n, q, r, rU, t];
```

and update `$Assumptions` on lines 36–38 from

```
$Assumptions =
  Element[{a0, delta, xi, m, n, q, r, rU}, Reals] &&
  a0 > 0 && delta > 0 && xi > 0 && rU > 0;
```

to

```
$Assumptions =
  Element[{a0, delta, xi, m, n, q, r, rU, t}, Reals] &&
  a0 > 0 && delta > 0 && xi > 0 && rU > 0;
```

Leave lines 97–116 unchanged (the `nSrcExpected` definition, `regThreshold`, `numZeroThreshold`, `dndmSrc`, and the two `expectZero` calls). After this edit, `expectZero["source-tied formula", nSrc - nSrcExpected]` on line 112 becomes a genuine substitution-vs-hand-form cross-check.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 041` and `redteam exec-mathematica 041` and confirm:
- Both scripts exit 0.
- The line `source-tied formula = 0` still appears in both output files.
- The line `tracking collapse = 0` (substantive check at section 24.3) still appears in both.
- A diff of `n_src` construction shows it now begins from `n_expected` (sympy) / `nExpected` (Mathematica), not from a fresh literal expression.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage041_rank2_support_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage041_rank2_support_mathematica_audit.wl`
- summary: Re-derived the source-tied support formula by substituting into the general result in both engines.
- deviation: none
