---
unit_id: 050
batch: III.2
created_at: 2026-05-22T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-22T00:00:00Z
findings_applied: 4
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 050

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:71-74`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:57-59`

**Issue:**
The "suppression factor" assertion `(2n+1)^2 * zeta_n - 1/(1 + x n(n+1)) == 0` reduces to a `(2n+1)^2 / (2n+1)^2` cancellation because `zeta_n` is set to the closed form `1/((2n+1)^2 (1 + x n(n+1)))` (imported in SymPy, redeclared in Mathematica). The check cannot fail and does not exercise the impossibility-bound claim that the docstring attributes to it.

**Required change:**
Replace the trivial cancellation check with a substantive identity that pins down the impossibility bound `zeta_req <= 1/(2n+1)^2`. The substantive content is: the threshold formula `x_max = (1/((2n+1)^2 zeta_req) - 1)/(n(n+1))` is non-negative iff `(2n+1)^2 zeta_req <= 1`. Concretely:

SymPy (replace lines 69-74 of `moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py`):

Before (lines 69-74):
```
banner("3. Exact impossibility bound from higher-harmonic suppression")
# Suppression bound relative to overlap only.
supp = sp.simplify((2 * n + 1) ** 2 * zeta_n)
print("(2n+1)^2 zeta_n^(twin) =", supp)
expect_zero("suppression factor - 1/(1+x n(n+1))", supp - 1 / (1 + x * n * (n + 1)))
print("Necessary condition for success: zeta_req <= 1/(2n+1)^2")
```

After:
```
banner("3. Exact impossibility bound from higher-harmonic suppression")
# x_max is non-negative iff (2n+1)^2 zeta_req <= 1; the numerator (2n+1)^2 zeta_req - 1
# of -x_max * n(n+1) factors that sign condition out cleanly.
admissibility_num = sp.together(sp.simplify(-x_eq * n * (n + 1) - (1 - (2 * n + 1) ** 2 * zeta_req) / ((2 * n + 1) ** 2 * zeta_req)))
print("admissibility numerator residual =", admissibility_num)
expect_zero(
    "x_max non-negativity reduces to (2n+1)^2 zeta_req <= 1",
    admissibility_num,
)
print("Therefore the necessary condition is exactly zeta_req <= 1/(2n+1)^2.")
```

Mathematica (replace lines 57-59 of `moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl`):

Before:
```
supp = FullSimplify[(2 n + 1)^2 zetaN, Assumptions -> $Assumptions];
Print["(2n+1)^2 zeta_n^(twin) = ", fmt[supp]];
expectZero["suppression factor - 1/(1+x n(n+1))", supp - 1/(1 + x n (n + 1))];
```

After:
```
admissibilityNum = FullSimplify[
  -xEq n (n + 1) - (1 - (2 n + 1)^2 zetaReq)/((2 n + 1)^2 zetaReq),
  Assumptions -> $Assumptions
];
Print["admissibility numerator residual = ", fmt[admissibilityNum]];
expectZero[
  "x_max non-negativity reduces to (2n+1)^2 zeta_req <= 1",
  admissibilityNum
];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 050` and `redteam exec-mathematica 050` and confirm the new check appears AND the script exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl`
- summary: Replaced the suppression-factor tautology with an admissibility residual check for the threshold bound in both audit scripts.
- deviation: Corrected the residual sign to match the stated `-x_max * n(n+1)` numerator identity; the literal snippet would not simplify to zero.

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:47,52-55`

**Issue:**
The Mathematica script defines `xEq` as `(1/(((2n+1)^2) zetaReq) - 1)/(n (n + 1))` on line 47, and lines 52-55 assert that `xEq` equals exactly that expression. This is a definitional tautology. The SymPy counterpart at least uses `sp.solve(Eq(zeta_n, zeta_req), x)` and compares to the closed form; the Mathematica script should match that level of independence.

**Required change:**
Modify line 47 and lines 52-55 of `moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl` so that `xEq` is derived via `Solve` (mirroring SymPy's `sp.solve`) and the subsequent assertion compares the derived solution to the closed form. Concretely:

Before (line 47):
```
xEq = FullSimplify[(1/(((2 n + 1)^2) zetaReq) - 1)/(n (n + 1)), Assumptions -> $Assumptions];
```

After:
```
xEqSolution = x /. First[Solve[zetaN == zetaReq, x]];
xEq = FullSimplify[xEqSolution, Assumptions -> $Assumptions];
xEqClosedForm = (1/(((2 n + 1)^2) zetaReq) - 1)/(n (n + 1));
```

Before (lines 52-55):
```
expectZero[
  "x_max - [1/((2n+1)^2 zeta_req)-1]/[n(n+1)]",
  xEq - (1/(((2 n + 1)^2) zetaReq) - 1)/(n (n + 1))
];
```

After:
```
expectZero[
  "x_max (from Solve) - [1/((2n+1)^2 zeta_req)-1]/[n(n+1)]",
  xEq - xEqClosedForm
];
```

The line-51 assertion `(zetaN /. x -> xEq) - zetaReq == 0` remains as-is — it now substantively checks that the solved `xEq` satisfies the threshold equation.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 050` and confirm `Solve` is invoked, the new label "x_max (from Solve)" appears, and the script exits 0.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl`
- summary: Derived `xEq` via `Solve` and compared that solved value against the closed-form threshold expression.
- deviation: none

## F3 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:76-85`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:61-67`

**Issue:**
Claim 5 says `S_n^(max)` is a *ceiling* on `S_n^(twin)(x)` for `x >= 0`. The current assertions only check the value at `x = 0`, which is consistent with either monotonic decrease, monotonic increase, or non-monotonic behavior. The "ceiling" claim is not exercised.

**Required change:**
Add (after the existing `expect_zero("S_n^(twin)(x=0) - S_n^(max)", ...)` line) an assertion that `S_n^(max) - S_n^(twin)(x)` factors as a manifestly non-negative expression for `x >= 0`, `0 < eps < 1`, `n >= 1`.

SymPy — insert after line 82, before the `S_1^(max)`/`S_2^(max)` prints:
```
# Ceiling check: S_n^(max) - S_n^(twin)(x) >= 0 for x >= 0.
ceiling_diff = sp.simplify(S_n_max - S_n)
ceiling_diff_factored = sp.together(ceiling_diff)
print("S_n^(max) - S_n^(twin) =", ceiling_diff_factored)
# Expected closed form: (1-eps) * (2n+1)^2 * n(n+1) * x  /  [ ((2n+1)^2 - eps) * ((2n+1)^2 (1 + x n(n+1)) - eps) ]
ceiling_diff_target = (
    (1 - eps) * (2 * n + 1) ** 2 * n * (n + 1) * x
    / (((2 * n + 1) ** 2 - eps) * ((2 * n + 1) ** 2 * (1 + x * n * (n + 1)) - eps))
)
expect_zero(
    "S_n^(max) - S_n^(twin) factored form",
    sp.simplify(ceiling_diff - ceiling_diff_target),
)
# Under 0 < eps < 1, n >= 1, x >= 0 each factor of ceiling_diff_target is >= 0,
# so S_n^(max) is indeed an upper bound (ceiling).
```

Mathematica — insert after line 65 of the .wl file, before the `S_1^(max)`/`S_2^(max)` prints:
```
ceilingDiff = FullSimplify[sNMax - sN, Assumptions -> $Assumptions];
Print["S_n^(max) - S_n^(twin) = ", fmt[ceilingDiff]];
ceilingDiffTarget =
  ((1 - eps) (2 n + 1)^2 n (n + 1) x) /
  (((2 n + 1)^2 - eps) ((2 n + 1)^2 (1 + n (n + 1) x) - eps));
expectZero[
  "S_n^(max) - S_n^(twin) factored form",
  ceilingDiff - ceilingDiffTarget
];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 050` and `redteam exec-mathematica 050` and confirm the new ceiling assertion appears AND the script exits 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl`
- summary: Added factored ceiling-difference checks showing `S_n^(max) - S_n^(twin)(x)` has the required non-negative form.
- deviation: none

## F4 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:56-67`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:46-55`

**Issue:**
The threshold claim is directional ("`zeta_n^(twin) >= zeta_req` gives `x <= x_max`"). The current assertions only check equality at `x = x_max`. The monotonic-decrease of `zeta_n^(twin)` in `x` (which fixes the direction of the inequality) is never asserted.

**Required change:**
Add a single assertion that `d zeta_n^(twin) / dx` equals its closed-form negative-signed expression.

SymPy — insert after line 60 (`print("zeta_n^(twin) =", zeta_n)`):
```
# Monotonicity check: zeta_n^(twin) is strictly decreasing in x, so
# zeta_n^(twin)(x) >= zeta_req iff x <= x_max.
d_zeta_n_dx = sp.simplify(sp.diff(zeta_n, x))
d_zeta_n_dx_target = -n * (n + 1) / ((2 * n + 1) ** 2 * (1 + x * n * (n + 1)) ** 2)
expect_zero(
    "d zeta_n^(twin) / dx + n(n+1) / [(2n+1)^2 (1 + x n(n+1))^2]",
    sp.simplify(d_zeta_n_dx - d_zeta_n_dx_target),
)
```

Mathematica — insert after line 49 (`Print["zeta_n^(twin) = ", fmt[zetaN]];`):
```
dZetaNdx = FullSimplify[D[zetaN, x], Assumptions -> $Assumptions];
dZetaNdxTarget = -n (n + 1) / ((2 n + 1)^2 (1 + n (n + 1) x)^2);
expectZero[
  "d zeta_n^(twin) / dx + n(n+1)/[(2n+1)^2 (1 + x n(n+1))^2]",
  dZetaNdx - dZetaNdxTarget
];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 050` and `redteam exec-mathematica 050` and confirm the new derivative-sign assertion appears AND the script exits 0.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl`
- summary: Added derivative checks showing `zeta_n^(twin)` has the required negative derivative in `x`.
- deviation: none
