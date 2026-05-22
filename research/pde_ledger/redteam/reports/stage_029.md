---
unit_id: 029
batch: II.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-21T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 029 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage029_dynamic_loading_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.txt`

## What the script claims to verify

The script claims a first-principles dynamic-loading lift of the Stage-11 profile selection model. Concretely, it verifies (1) the exact Schur-complement elimination of the coupled u/W/phi internal-field block reproduces a wall self-energy of the form `Sigma(omega) = Xi(omega) I + alpha(omega) v v^T` with explicit closed-form `Xi` and `alpha`; (2) the conservative static data `Xi_0, Delta_0, alpha_0` are recovered at `omega = 0, Pi = 0`; (3) the refined effective stiffness `K_eff` produces a stationarity equation `dE/dtheta = (DeltaK + alpha_0 xi) sin(2theta) - 2 alpha_0 eta cos(2theta)` whose root gives `tan(2 theta_-) = 2 alpha_0 eta / (DeltaK + alpha_0 xi)`; (4) the softening determinant has a closed-form root `alpha_crit = K0t K1t / (K1t kappa0^2 + K0t kappa1^2)`; (5) the first-order expansion of `alpha(omega)` in the outgoing port variable `Pi` produces a clean transfer factor `beta(omega) = (A_U lambda_W + lambda_R lambda_U)^2 / Delta_cons^2`; (6) the `Pi = i Gamma_port omega^5` substitution yields the expected coefficient `beta_5` extractable from the `omega^5` term of the series; and (7) the selected lower-mode odd quadrupole coefficient obeys the Hellmann–Feynman identity `kappa_sel^2 = -d lambda_-/d alpha` with correct weak (`al -> 0 -> kappa0^2`) and strong (`al -> infty -> sigma`) limits.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1  | sympy | 127  | `expect_zero("Sigma - (Xi I + alpha vv^T)", ...)` | yes |
| A2  | sympy | 153  | `expect_zero("DeltaK_tilde - DeltaK_bare (Xi_0 cancellation)", DeltaK - DeltaK_ax)` | no (tautological) |
| A3  | sympy | 182  | `expect_zero("dE/dtheta - stationarity/2", dE - stationarity/2)` | yes |
| A4  | sympy | 192  | `expect_zero("alpha_crit - expected", alpha_crit - alpha_crit_expected)` | yes |
| A5  | sympy | 214  | `expect_zero("beta - clean transfer factor", beta - beta_clean)` | yes |
| A6  | sympy | 218  | `expect_zero("alpha - (alpha_cons + beta*Pi) at O(Pi)", ...)` | yes |
| A7  | sympy | 230  | `expect_zero("extracted beta_5 - expected beta_5", ...)` | yes |
| A8  | sympy | 276  | `expect_zero("weak-loading kappa_sel^2 -> kappa0^2", ...)` | yes |
| A9  | sympy | 277  | `expect_zero("strong-loading kappa_sel^2 -> sigma", ...)` | yes |
| A10 | mma   | 89   | `expectMatrixZero["Sigma - (Xi I + alpha vv^T)", ...]` | yes |
| A11 | mma   | 90-92| sigma/xi/eta literal-constant cross-checks | yes |
| A12 | mma   | 103  | `expectZero["DeltaK_tilde - DeltaK_ax", (k1t - k0t) - DeltaKax]` | no (tautological) |
| A13 | mma   | 119  | `expectZero["dE/dtheta - stationarity/2", ...]` | yes |
| A14 | mma   | 120-123 | `expectZero["-tan(2 theta_-) - manifestly positive form", -tan2Theta - 2*alpha0*(-eta)/(DKax+alpha0*xi)]` | no (tautological) |
| A15 | mma   | 133  | `expectZero["det(alpha_crit)", detTemplate /. al -> alphaCrit]` | yes |
| A16 | mma   | 142  | `expectZero["beta - clean transfer factor", beta - betaClean]` | yes |
| A17 | mma   | 146  | `expectZero["alpha - (alpha_cons + beta portPi) at O(portPi)", ...]` | yes |
| A18 | mma   | 155  | `expectZero["beta_5 - GammaPort (OmegaU^2 lambdaW + lambdaR lambdaU)^2/Delta0^2", ...]` | yes |
| A19 | mma   | 156  | `expectZero["extracted beta_5 - expected beta_5", ...]` | yes |
| A20 | mma   | 167  | `expectZero["weak-loading kappa_sel^2 - kappa0^2", ...]` | yes |
| A21 | mma   | 168  | `expectZero["strong-loading kappa_sel^2 - sigma", ...]` | yes |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:113-123`

**What's wrong:**
The Mathematica script defines

```
tan2Theta = FullSimplify[2*alpha0*eta/(DeltaKax + alpha0*xiConst), ...];
...
expectZero[
  "-tan(2 theta_-) - manifestly positive form",
  -tan2Theta - 2*alpha0*(-eta)/(DeltaKax + alpha0*xiConst)
];
```

Substituting the definition of `tan2Theta` into the residual gives

```
-(2*alpha0*eta/(DKax + alpha0*xiConst)) - 2*alpha0*(-eta)/(DKax + alpha0*xiConst)
  = -2 alpha0 eta / D + 2 alpha0 eta / D
  = 0
```

i.e. the residual is structurally `-x - (-x) = 0` for any value of `alpha0`, `eta`, `DeltaKax`, `xiConst`. The check is algebraically guaranteed by construction; no physics is being tested. The comment / label suggests it is verifying a "manifestly positive form" of `tan(2 theta_-)`, but the right-hand side of the residual is just a re-bracketing of the same expression with the minus sign moved from the leading factor into the numerator.

**Why this matters:**
A PASS here gives no evidence about the angle law, the sign of `eta`, or anything else physical. It dilutes the meaning of the script's overall PASS line. If the genuine claim is that the numerator `−eta = κ0|κ1|` is positive given the chosen `kappa0 > 0, kappa1 < 0`, then the check should verify positivity (e.g. `Sign[-eta] == 1` or `Refine[-eta > 0, Assumptions -> kappa0 > 0 && kappa1 < 0]`), not the identity `−x ≡ −x`.

**Required change:**
Replace the tautological residual on lines 120–123 with a check that actually exercises the angle-law content. Two viable substantive forms (pick one):

(a) Verify the stationarity-equation root by direct substitution:
```
expectZero[
  "stationarity at theta_-",
  ((DeltaKax + alpha0*xiConst)*Sin[2*theta] - 2*alpha0*eta*Cos[2*theta])
    /. theta -> ArcTan[2*alpha0*eta/(DeltaKax + alpha0*xiConst)]/2
];
```

(b) Or alternatively verify the sign claim — that the actual numerator `−eta = 8 Sqrt[2]/(3 Pi^2)` is positive in the chosen overlap convention:
```
expectZero["-eta - 8 Sqrt[2]/(3 Pi^2)", -eta - 8*Sqrt[2]/(3*Pi^2)];
```
which the script does NOT currently check (it only checks `eta + 8 Sqrt[2]/(3 Pi^2) = 0`, which is the same identity).

The SymPy script does not assert this identity at all, so removing the Mathematica check (without replacement) is also acceptable; if removed, also delete the `tan2Theta` symbol on line 113 if no longer used downstream.

**Verification:**
After patch, `redteam exec-mathematica 029` should still exit 0; the captured output should contain a new substantive PASS line (the replacement check) or — if the check was simply deleted — should NOT contain "manifestly positive form" anywhere.

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:146-153`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:97-103`

**What's wrong:**
Both engines define

```
K1 = K0 + DeltaK_ax        (sympy line 67 / mma line 47)
K0t = K0 - Xi0
K1t = K1 - Xi0
DeltaK = K1t - K0t         (sympy line 148 / implicit on mma line 103 via (k1t - k0t))
```

then assert `DeltaK - DeltaK_ax == 0`. Expanding the definitions in either engine:

```
K1t - K0t = (K0 + DeltaK_ax - Xi0) - (K0 - Xi0) = DeltaK_ax
```

so `DeltaK - DeltaK_ax = 0` is forced by the very definitions on the lines immediately above — no matter what `Xi_0` is, no matter what the physics is, the `Xi_0` shift cancels under subtraction. The accompanying comment "Verify the isotropic shift cancels from the stiffness splitting" overstates what the check accomplishes: the check is "subtraction cancels what is subtracted twice", which is trivially true.

**Why this matters:**
Lower severity than F1 because this is at least an honest sanity check on the construction (catches a typo where someone wrote `K1t = K1 + Xi0` instead of `K1t = K1 - Xi0`). But the script's comment claims it as verifying *physics* of the isotropic shift, when in fact it can only catch a typing error in the immediately preceding two lines.

**Required change:**
Either (a) downgrade the inline comments to make clear this is a typo-guard, not a physics check; or (b) replace the trivial check with a non-trivial one that actually exercises the isotropic-shift cancellation in a context where it's not algebraically forced.

The minimum-effort fix is (a). In SymPy at line 152, change the comment from:
```
# Verify the isotropic shift cancels from the stiffness splitting.
```
to:
```
# Sanity check: the K0t/K1t construction yields DeltaK_tilde = DeltaK_ax. This is
# enforced by the definitions on the two preceding lines (trivially true under
# the construction K1 = K0 + DeltaK_ax, K0t = K0 - Xi0, K1t = K1 - Xi0) and is
# kept here only as a typo guard.
```

In Mathematica at line 102 (immediately above the `expectZero[...]` on line 103), insert the same clarifying comment:
```
(* Sanity check: K0t = K0 - Xi0, K1t = K1 - Xi0, K1 = K0 + DeltaKax, so
   (K1t - K0t) - DeltaKax = 0 is algebraically forced. Kept as typo guard. *)
```

Do NOT delete the assertion itself; the typo guard has nonzero value.

**Verification:**
After patch, the assertion remains in both scripts (verifier still sees the PASS line in both output transcripts) but the docstring/comment honestly describes what it checks.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:74-89` (Schur block)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:161-168` (Hellmann–Feynman block)
- (versus the corresponding SymPy blocks at `scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:97-131` and `:265-277`)

**What's wrong:**
The Mathematica script is structurally a line-by-line port of the SymPy script. Side-by-side excerpts:

Schur block (SymPy `:101-127` vs Mathematica `:74-89`):
```
# SymPy
Mint = sp.Matrix([
    [AU, 0, -lambda_R * kappa0, 0],
    [0, AU, -lambda_R * kappa1, 0],
    [-lambda_R * kappa0, -lambda_R * kappa1, AW, 0],
    [0, 0, 0, Aphi],
])
C = sp.Matrix([
    [lambda_U, 0],
    [0, lambda_U],
    [lambda_W * kappa0, lambda_W * kappa1],
    [lambda_B * kappa0, lambda_B * kappa1],
])
Sigma = sp.simplify(C.T * Mint.inv() * C)
Sigma_expected = sp.simplify(Xi * I2 + alpha * (v * v.T))
expect_zero("Sigma - (Xi I + alpha vv^T)", Sigma - Sigma_expected)
```
```
(* Mathematica *)
mint = {
  {aU, 0, -lambdaR*kappa0, 0},
  {0, aU, -lambdaR*kappa1, 0},
  {-lambdaR*kappa0, -lambdaR*kappa1, aW, 0},
  {0, 0, 0, aphi}
};
cMat = {
  {lambdaU, 0},
  {0, lambdaU},
  {lambdaW*kappa0, lambdaW*kappa1},
  {lambdaB*kappa0, lambdaB*kappa1}
};
sigmaWall = FullSimplify[Transpose[cMat].LinearSolve[mint, cMat], ...];
sigmaExpected = FullSimplify[xiShift*i2 + alpha*Outer[Times, v, v], ...];
expectMatrixZero["Sigma - (Xi I + alpha vv^T)", sigmaWall - sigmaExpected];
```

The internal-field ordering `(u0, u1, W, phi)`, the coupling matrix `C`, the choice to compute `Sigma = C^T M^{-1} C` in one shot, and the comparison form `Xi*I + alpha*vv^T` are bit-for-bit identical. The only differences are surface syntax (`sp.Matrix` vs `{...}`, `Mint.inv()` vs `LinearSolve[mint, ...]`, `simplify` vs `FullSimplify`).

Hellmann–Feynman block (SymPy `:266-277` vs Mathematica `:161-168`):
```
# SymPy
disc = sp.simplify((DeltaK + al * xi) ** 2 + 4 * al**2 * eta**2)
tr_eff = sp.simplify((K0t + K1t) - al * sigma)
lam_minus_template = sp.simplify((tr_eff - sp.sqrt(disc)) / 2)
kappa_sel_template = sp.simplify(-sp.diff(lam_minus_template, al))
...
expect_zero("weak-loading kappa_sel^2 -> kappa0^2", ...kappa_sel_template.subs(al, 0) - kappa0**2)
expect_zero("strong-loading kappa_sel^2 -> sigma", ...sp.limit(kappa_sel_template, al, sp.oo) - sigma)
```
```
(* Mathematica *)
discTemplate = FullSimplify[(DeltaKax + al*xiConst)^2 + 4*al^2*eta^2, ...];
trTemplate = FullSimplify[k0t + k1t - al*sigma, ...];
lambdaMinusTemplate = FullSimplify[(trTemplate - Sqrt[discTemplate])/2, ...];
kappaSelSq = FullSimplify[-D[lambdaMinusTemplate, al], ...];
...
expectZero["weak-loading kappa_sel^2 - kappa0^2", (kappaSelSq /. al -> 0) - kappa0Sq];
expectZero["strong-loading kappa_sel^2 - sigma", FullSimplify[Limit[kappaSelSq, al -> Infinity], ...] - sigma];
```

Again: same `disc`/`tr_eff`/`lam_minus` form, same `−d/d(al)` Hellmann–Feynman formula, same two limit checks at exactly the same two limit points.

This violates the second-engine policy: the Mathematica script is echoing the SymPy script's algebra rather than re-deriving the result from independent intermediate forms. A genuine independent derivation would (a) eliminate the internal fields one at a time (eliminate `phi` first using the trivial `Aphi * phi = ...` block, then `W` using `aW*W = lambdaR*(kappa0 u0 + kappa1 u1)`, then the `U` doublet) to obtain `Sigma_wall` without ever inverting the 4×4, OR (b) derive `kappa_sel^2` by directly diagonalising `Keff0` at `alpha = alpha_0` and reading off the projection onto the lower eigenvector, OR (c) verify the angle law by an independent method such as completing the square on `(K0t - alpha0 kappa0^2) cos^2 theta + (K1t - alpha0 kappa1^2) sin^2 theta - 2 alpha0 kappa0 kappa1 sin theta cos theta`.

**Why this matters:**
If both engines simply echo the same algebra, a bug in either author's setup (e.g. a sign in `Mint`, a transposition of the coupling matrix, a wrong sign on `eta`) will be reproduced verbatim in both scripts and both will PASS. The second engine exists precisely to break this kind of correlated failure.

**Required change:**
Restructure the Mathematica script's Schur block (lines 74–89) and Hellmann–Feynman block (lines 161–168) to use independent derivation paths. Concretely:

(1) Replace lines 74–89 with a sequential-elimination derivation: eliminate `phi`, `W`, and the `U` doublet one at a time from the coupled 4×4 system, accumulating the effective wall self-energy as a sum of three rank-correction contributions (from `phi`, from `W`, from `U`). Compare the resulting `Sigma_wall_seq` to the same `Sigma_expected = Xi I + alpha vv^T` target. The check on line 89 then verifies the *sequentially-derived* `Sigma_wall_seq` agrees with the algebraic decomposition, not the one-shot 4×4 inverse.

(2) Replace lines 161–168 with a derivation that does NOT use the closed-form `lam_minus = (tr - Sqrt[disc])/2` and `-D[lam_minus, al]`. Instead, at `alpha = alpha_0`, diagonalise `Keff0` numerically over the symbolic ring via `Eigensystem[Keff0]` (Mathematica supports this for symbolic 2×2), take the eigenvalue with the smaller real part, and read `kappa_sel^2 = | <v_-, v> |^2` where `v_- = (cos theta_-, sin theta_-)`. Then the assertion is `(kappa_sel^2 - "closed-form formula" )`, with the closed form being the SymPy script's `−d lambda_- / d al` — and that becomes the genuine cross-engine check.

**Verification:**
After Codex applies, `redteam exec-mathematica 029` exits 0, the captured output (a) contains a new printed line `Sigma_wall_seq = ...` distinct from `Sigma_wall_oneshot`, and (b) the printed `kappa_sel^2` line is derived from `Eigensystem[Keff0]` rather than `-D[lambdaMinusTemplate, al]`. The verifier visually confirms the two derivation paths in the source.

## Independent-derivation check (Mathematica)

The Mathematica script is NOT an independent re-derivation. See F3 above for line-by-line correspondence with the SymPy script. The shared elements include: the identical Mint matrix layout with fields ordered `(u0, u1, W, phi)`, the identical coupling matrix C with the same rows in the same order, the one-shot Schur formula `C^T M^{-1} C` (vs sequential elimination), the same Hellmann–Feynman template `−D[lam_minus, al]`, and the same two limit checks (`al -> 0` and `al -> Infinity`). The only divergence is surface syntax.

## Engine cross-check

Both engines' final outputs agree at the level they claim to. Spot-check:

- `sigma`: sympy `88/(9*pi**2)` (line 353) = mma `88/(9 Pi^2)` (line 19 PASS). ✓
- `xi`: sympy `56/(9*pi**2)` (line 354) = mma `56/(9 Pi^2)` (line 21 PASS). ✓
- `eta`: sympy `-8*sqrt(2)/(3*pi**2)` (line 355) = mma `-8 Sqrt[2]/(3 Pi^2)` (line 23 PASS). ✓
- `Xi_0`: both `lambda_U^2/Omega_U^2`. ✓
- `alpha_0`: SymPy line 366 expansion equals the Mathematica line 27 form after clearing the `9 pi^2` factor between numerator and denominator. ✓
- `beta_5`: sympy `81*pi^4*Gamma_port*(...)^2/(9*pi^2*Omega_U^2*Omega_W^2 - 88*lambda_R^2)^2` (line 829) equals mma `Gamma_port*(...)^2/(Omega_U^2*Omega_W^2 - 88*lambda_R^2/(9*Pi^2))^2` (line 45) after factoring `9 pi^2` from the denominator and squaring. ✓
- `kappa(theta)^2`: sympy `(8/(9 pi^2))(11/2 + (7/2) cos(2 theta) - 3 sqrt(2) sin(2 theta))` = mma `(4/(9 Pi^2))(11 + 7 Cos[2 theta] - 6 Sqrt[2] Sin[2 theta])`. ✓
- All `expect_zero` / `expectZero` residuals print as `0` / `PASS` in both engines.

No engine disagreement.

## Verdict justification

Three findings: F1 (a tautological Mathematica check of the form `-x - (-x) = 0`, mislabeled as a positivity claim), F2 (a tautological check in both engines that subtraction of equal shifts cancels — kept honestly as a typo guard but mislabeled in comments as a physics check), and F3 (the Mathematica script is structurally a line-by-line port of the SymPy script, providing no genuine independent verification). The unit's substantive math (Schur decomposition, alpha_crit closed form, beta_5 outgoing dressing, Hellmann–Feynman selected-mode coefficient) holds up under attack — assertions are well-anchored, residuals genuinely vanish, and engine cross-check is clean. None of the findings is `UNFIXABLE` or `CRITICAL_DOWNSTREAM`: F1/F2 are local cosmetic/labeling fixes; F3 is a restructuring that does not change the proven result.

Attacks attempted that failed: (i) substituting the trivial profile `kappa0 = 1, kappa1 = 0` into the Schur decomposition — the residual still simplifies to zero, no hidden branch error in the `Mint^{-1}` computation; (ii) checking that the `assume(positive=True)` on `Omega_U, Omega_W, varpi, K0, M0, M1, DeltaK_ax` does not collapse a sign in `eta = kappa0 kappa1` — the assumptions affect only the physical-parameter symbols, while `eta`'s sign comes from the explicit numeric definition of `kappa1 = -4/(3*pi) < 0`; (iii) testing whether `sp.series(alpha.subs(Pi, eps), eps, 0, 2)` could agree with `alpha_cons + beta_clean * Pi` for trivial reasons even if `beta_clean` were wrong — but `beta_clean` is independently asserted to equal `D(alpha, Pi)|_{Pi=0}` on line 214 first, so the O(Pi) series check is consistent rather than vacuous; (iv) confirming the `omega -> 0` limit in `beta_5 = beta_clean.subs(omega, 0)*Gamma_port` does not silently drop a finite term — direct substitution shows `beta_clean(omega=0) = (Omega_U^2 lambda_W + lambda_R lambda_U)^2 / Delta_0^2`, which matches the `beta_5` reported in both outputs.

## Self-test notes

I checked each proposed directive edit before writing it. (a) For F1, I substituted the definition `tan2Theta = 2 a0 eta / D` into the residual `-tan2Theta - 2 a0 (-eta)/D` and obtained `-2 a0 eta/D + 2 a0 eta/D = 0` confirming the tautology; for the proposed replacement (a), `((DeltaKax + alpha0*xiConst) Sin[2t] - 2 alpha0 eta Cos[2t]) /. t -> ArcTan[2 alpha0 eta/(DKax+alpha0 xi)]/2` is a genuine equation `R sin(phi-phi0) = 0` whose root `2t = arctan(2 a0 eta/(DKax+a0 xi))` does make it vanish, so the replacement is non-trivial and PASS-able. (b) For F2, the assertion `(K1t - K0t) - DeltaKax = ((K0 + DeltaKax - Xi0) - (K0 - Xi0)) - DeltaKax = 0` is forced by definitions on the two preceding lines; my proposed fix only edits the comment, not the assertion, so no derivative/parity/trivial-case re-derivation is needed. (c) For F3, the directive proposes restructuring the Mathematica script; I verified that `Eigensystem` over a symbolic 2×2 in Mathematica does produce closed-form eigenvalues/eigenvectors, and that the sequential-elimination derivation of `Sigma_wall` (eliminate phi → W → U doublet) does land on the same `Xi I + alpha vv^T` form, so the proposed independent paths are not degenerate. (d) Path specifications: F1 and F3 target the `.wl` file under `mathematica/`; F2 targets both the `.py` under `scripts/` and the `.wl` under `mathematica/`. All paths are absolute and correctly placed.
