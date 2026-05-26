Stage 043 Mathematica audit FAILED on the F2-Insertion2 numeric-point sign anchor you added in iter1.

The failing line (from `redteam/exec_logs/stage_043_mathematica.log`):
```
mismatch sign at deltaU=1, sigma0=0, rho0=1 = (-1 - (2*gS*gU)/(gS*gU + gB*kU) + (2*gR*gU)/(gR*gU + gW*kU))/4
FAIL: mismatch sign at deltaU=1, sigma0=0, rho0=1
  residual -> (-1 - (2*gS*gU)/(gS*gU + gB*kU) + (2*gR*gU)/(gR*gU + gW*kU))/4
```

The issue: you substituted `sigma0=0, rho0=1, deltaU=1` directly into the Mathematica expression, but the expression's primary symbols are `gB, gU, gS, gR, gW, kU` — not `sigma0` and `rho0`. So `sigma0` and `rho0` are treated as fresh unrelated symbols, and the substitution doesn't actually reduce the mismatch formula to a numeric.

The closed forms (from notes and from sympy):
- `sigma_0 = gU * gS / (kU * gB)`
- `rho_0 = gU * gW / (kU * gR)`

To realize `sigma_0 = 0`: set `gS -> 0` (assuming `gU, kU, gB` nonzero). To realize `rho_0 = 1`: set `gW -> kU * gR / gU` (algebraically rearranged from `gU * gW / (kU * gR) = 1`).

# Required fix

Replace the failing assertion's substitution rules with explicit primitive substitutions. For the `(deltaU=1, sigma0=0, rho0=1)` test case:

```mathematica
mismatchAtTestPoint1 = FullSimplify[
  (rPhi - rU) /. {deltaU -> 1, gS -> 0, gW -> kU*gR/gU},
  Assumptions -> $Assumptions
];
expectZero[
  "mismatch sign at deltaU=1, sigma0=0, rho0=1",
  mismatchAtTestPoint1 - 0
];
```

(With `sigma_0 = 0`, the closed-form mismatch `R_phi - R_U = deltaU * (rho_0 - sigma_0)/[(1+deltaU)(1+rho_0)(1+sigma_0)] = 1 * 1 / (2 * 2 * 1) = 1/4`, NOT 0. So actually the test value should be `1/4`, not 0. Re-check what the auditor's directive said the expected numeric should be.)

Read the original directive's F2-Insertion2 block carefully — it specified THREE test triples each with a known numeric value of `R_phi - R_U`. Use those exact expected values; do NOT just assert "= 0".

If your iter1 implementation was sloppy on the expected RHS or on the substitution rules, fix both:
1. Substitution rules must use primitive symbols `gS, gW` (not derived `sigma0, rho0`).
2. Expected RHS should be the directive's prescribed numeric value for each triple — likely `+1/4`, `-1/4`, `0`, etc.

# Process

1. Read `redteam/directives/stage_043.md` § F2-Insertion2 to see the directive's exact prescribed test triples and their expected numeric values.
2. Re-derive the closed form: `R_phi - R_U = deltaU * (rho_0 - sigma_0) / [(1 + deltaU)(1 + rho_0)(1 + sigma_0)]`. At each triple, evaluate this number by hand and use it as the expected RHS.
3. Rewrite the Mathematica test points using `gS -> 0` (for sigma_0=0) and `gW -> kU*gR/gU` (for rho_0=1) substitution chains — or whatever combinations match the directive's three triples.
4. Run `math -script mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl`. Confirm exit 0. Single-seat — no other math is running.
5. Append `## Applied: F2-Insertion2-iter2` block to the directive (preserving the iter1 Applied block as history) with: files_changed, summary, deviation.

# Working directory

`/var/projects/toy_physics/research/pde_ledger`
