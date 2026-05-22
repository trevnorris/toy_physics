---
unit_id: 032
batch: II.1
auditor_model: claude-opus-4-7-1m
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

# Audit unit 032 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.txt`

## What the script claims to verify

The unit (labelled internally as "Stage 15") verifies five distinct things:
(15.1) Closed-form values of finite-throat axial overlap integrals `kappa0 = ∫u0·f0` and `kappa1 = ∫u1·f0`, plus orthonormality of the N/N basis `{u0,u1}` and the D/N half-wave `f0`.
(15.2) Five coupling integrals collapse to the expected mode-level expressions involving `(kappa0, kappa1)`.
(15.3) Schur complement of a 4x4 internal block matrix `K_int` equals the rank-1+identity form `Xi·I + alpha·v v^T` with explicit scalar `Xi`, `alpha`.
(15.4) A stated closed form `s_minus/kappa0²` for the "mhat_- squared" map has the right alpha0→0 and alpha0→∞ limits under the natural D/N substitution.
(15.5) The product `(s_minus/kappa0²)·(beta0 s_minus/lambda_minus)` equals `beta0 s_minus²/(kappa0² lambda_minus)`.

## Assertion inventory

| #   | Script      | Line | Form                                                                                  | Anchored to claim? |
|-----|-------------|------|---------------------------------------------------------------------------------------|--------------------|
| A1  | sympy       | 58   | `expect_zero("u0.u0 - 1", ...)`                                                       | yes                |
| A2  | sympy       | 59   | `expect_zero("u1.u1 - 1", ...)`                                                       | yes                |
| A3  | sympy       | 60   | `expect_zero("u0.u1", ...)`                                                           | yes                |
| A4  | sympy       | 61   | `expect_zero("f0.f0 - 1", ...)`                                                       | yes                |
| A5  | sympy       | 62   | `expect_zero("kappa0 - 2 sqrt(2)/pi", ...)`                                           | yes                |
| A6  | sympy       | 63   | `expect_zero("kappa1 + 4/(3 pi)", ...)`                                               | yes                |
| A7  | sympy       | 64   | `expect_zero("sigma - 88/(9 pi^2)", ...)`                                             | yes                |
| A8  | sympy       | 65   | `expect_zero("sigma/kappa0^2 - 11/9", ...)`                                           | yes                |
| A9  | sympy       | 91-94 | `L_etaU - gU(q.U) == 0`                                                              | yes (linearity)    |
| A10 | sympy       | 95-98 | `L_etaphi - gB (v.q) phi == 0`                                                       | yes (linearity)    |
| A11 | sympy       | 99-102 | `L_etaW - gW (v.q) W == 0`                                                          | yes (linearity)    |
| A12 | sympy       | 103-106 | `L_UW + gR (v.U) W == 0`                                                           | yes (linearity)    |
| A13 | sympy       | 107-110 | `L_src - gQ Q (v.q) == 0`                                                          | yes (linearity)    |
| A14 | sympy       | 143  | `Sigma - [Xi I + alpha vv^T] == 0`                                                    | yes                |
| A15 | sympy       | 166  | `mhat_-^2(alpha=0) - 1 == 0`                                                          | partial (endpoint) |
| A16 | sympy       | 167  | `limit_{alpha->oo} mhat_-^2 - 11/9 == 0`                                              | partial (endpoint) |
| A17 | sympy       | 177  | `mhat^2 P0_- - beta0 s^2/(kappa0^2 lambda_-) == 0`                                    | **no (tautological)** |
| B1-B17 | mathematica | 51-138 | same form as A1-A17                                                                | mirrors sympy      |

A15 and A16 are "partial" because the script *defines* (does not derive) the closed form `s_minus` at lines 153-155 and only verifies it at the two endpoints of `alpha0`.

A17 is the central tautology. The two quantities being equated at line 177 are algebraically identical by construction (multiplicative associativity).

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py:169-177`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl:134-138`

**What's wrong:**
Stage 15.5 (the "ELIMINATION OF THE ABSTRACT SOURCE-MAP FACTOR" stage) builds two quantities and asserts they are equal:

```
P0_from_stage13 = beta0 * s_minus / lambda_minus           # line 174
Nprod          = (s_minus / kappa0_sq) * P0_from_stage13   # line 175
Nprod_target   = beta0 * s_minus**2 / (kappa0_sq * lambda_minus)  # line 176
expect_zero("mhat^2 P0_- - beta0 s^2/(kappa0^2 lambda_-)", Nprod - Nprod_target)  # line 177
```

Substituting line 174 into line 175 gives `Nprod = (s_minus/kappa0_sq) * (beta0 * s_minus / lambda_minus) = beta0 * s_minus**2 / (kappa0_sq * lambda_minus)`, which is *literally* `Nprod_target`. The residual `Nprod - Nprod_target` is identically zero by associativity/commutativity of multiplication — independent of what `s_minus`, `lambda_minus`, `beta0`, or `kappa0_sq` are. Notice that `kappa0_sq` is declared as a fresh free symbol at line 171 with no relation to `kappa0**2`, and `s_minus`, `lambda_minus` are not substituted/evaluated. The assertion cannot fail regardless of the physics.

The Mathematica mirror at lines 135-138 has the same structure:
```
p0Minus = FullSimplify[beta0*sMinus/lamMinus, ...];
nProd = FullSimplify[(sMinus/kappa0^2)*p0Minus, ...];
nProdTarget = FullSimplify[beta0*sMinus^2/(kappa0^2*lamMinus), ...];
expectZero["mhat^2 P0_- - beta0 s^2/(kappa0^2 lambda_-)", nProd - nProdTarget];
```
This is the same identity; nothing is exercised.

**Why this matters:**
The stage's banner says "ELIMINATION OF THE ABSTRACT SOURCE-MAP FACTOR" and the docstring lists "Elimination of the abstract selected-branch source-map factor." A passing assertion gives the impression that an elimination has been verified, but no physical content is being checked — only the distributive property `(a/b)*(c/d) = ac/(bd)`. Any change to the closed form of `s_minus` or `lambda_minus` would still leave this assertion passing.

**Required change:**
Replace the tautological assertion with one that uses the natural-D/N substitution from Stage 15.4 (`subs_nat`) and a `kappa0_sq -> kappa0**2` substitution to evaluate `Nprod` at the two limits already considered by Stage 15.4 (alpha0=0 and alpha0→∞) and compare against an explicit closed form. Specifically, after applying `subs_nat` and `kappa0_sq -> kappa0**2`:

1. At `alpha0 = 0`: `lambda_minus = A` (since `R = sqrt(DK**2) = DK` for `DK > 0`), `s_minus_nat = kappa0**2 = 8/pi**2`. So `Nprod` (alpha=0) should equal `beta0 * (8/pi**2) / A`. Assert this.
2. As `alpha0 -> oo`: `Nprod` should tend to 0 (because `lambda_minus ~ -alpha0*sigma` diverges while `s_minus_nat -> 88/(9 pi**2)` stays finite, giving `Nprod ~ (11/9)*beta0/(-alpha0)`). Assert that `sp.limit(Nprod_subbed, alpha0, sp.oo) == 0`.

Mirror the same two checks in the Mathematica script.

**Verification:**
After Codex applies, the SymPy script should contain at least two new `expect_zero` calls in Stage 15.5 that reference `subs_nat` and an `{kappa0_sq: kappa0**2}` substitution, AND a `sp.limit(..., alpha0, sp.oo)` call. The Mathematica script should contain matching `expectZero` calls using `/. subsNat /. {kappa0_sq -> ...}` (no separate `kappa0_sq` symbol exists in `.wl`; the `.wl` already uses `kappa0^2` directly, so only the natural substitution and limit need to be added). Both outputs should still exit 0.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl` (entire file)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py` (entire file, as reference)

**What's wrong:**
The `.wl` is a line-by-line transliteration of the `.py`. Three correspondences make this clear:

(a) **Identical stage structure and assertion strings.** Both scripts use the literal banner strings "STAGE 15.1 — EXACT FINITE-THROAT AXIAL INTEGRALS" through "STAGE 15.5 — ELIMINATION OF THE ABSTRACT SOURCE-MAP FACTOR", and every `expect_zero` / `expectZero` name string is identical (e.g. `"kappa1 + 4/(3 pi)"`, `"L_etaphi - gB (v.q) phi"`, `"Sigma - [Xi I + alpha vv^T]"`, `"mhat^2 P0_- - beta0 s^2/(kappa0^2 lambda_-)"`).

(b) **Identical variable choreography in Stage 15.3.** The SymPy script builds:
```
Kint = sp.Matrix([
    [AU, 0, 0, -gR * kappa0],
    [0, AU, 0, -gR * kappa1],
    [0, 0, Aphi, 0],
    [-gR * kappa0, -gR * kappa1, 0, AW],
])
Bmat = sp.Matrix([[gU, 0, gB*kappa0, gW*kappa0], [0, gU, gB*kappa1, gW*kappa1]])
Sigma = Bmat * Kint.inv() * Bmat.T
```
The `.wl` mirrors this exactly with the same row/column ordering, same block structure, same matrix-product expression:
```
kInt = {{aU, 0, 0, -gR*kappa0}, {0, aU, 0, -gR*kappa1},
        {0, 0, aPhi, 0}, {-gR*kappa0, -gR*kappa1, 0, aW}};
bMat = {{gU, 0, gB*kappa0, gW*kappa0}, {0, gU, gB*kappa1, gW*kappa1}};
sigmaMat = FullSimplify[bMat.Inverse[kInt].Transpose[bMat], ...];
```
The only differences are case conventions (`AU` -> `aU`, `Aphi` -> `aPhi`). No independent derivation; same algebra rewritten.

(c) **Identical Stage 15.4 formulas.** SymPy lines 153-155:
```
R = sp.sqrt((DK + alpha0 * delta_kappa)**2 + 4 * alpha0**2 * Kprod)
lambda_minus = (A + B - alpha0 * sigma_sym - R) / 2
s_minus = (1/2) * (sigma_sym + ((DK + alpha0 * delta_kappa) * delta_kappa + 4 * alpha0 * Kprod) / R)
```
`.wl` lines 123-125:
```
r = Sqrt[(dK + alpha0*deltaKappa)^2 + 4*alpha0^2*kappaProd];
lamMinus = FullSimplify[(2*a + dK - alpha0*sigmaSym - r)/2, ...];
sMinus = FullSimplify[(sigmaSym + ((dK + alpha0*deltaKappa)*deltaKappa + 4*alpha0*kappaProd)/r)/2, ...];
```
Same formulas, same operand order, same parenthesisation.

This violates the second-engine policy: both engines are supposed to derive results independently and only converge at the final claim.

**Why this matters:**
A transliterated `.wl` cannot catch algebraic errors in the `.py` — any sign flip, factor of 2, or branch error in the SymPy script is faithfully reproduced in Mathematica, so engine agreement provides no additional confidence.

**Required change:**
Restructure the Mathematica script so each stage derives its result by a method distinct from the SymPy choreography. Concretely:

- Stage 15.1: keep the closed-form integrals (these are unique up to method-of-evaluation, no choice exists), but additionally verify `kappa0` and `kappa1` by an alternative route: differentiating an antiderivative chosen by Mathematica's `Integrate[..., GenerateConditions -> False]` against `s` and re-checking, or by computing the same integral via `Sum[Residue[...]]` / contour expansion. At minimum, evaluate via the substitution `s -> l*x` and check the resulting `x`-integral equals the published rational multiple of `1/Pi`.
- Stage 15.3: instead of computing `Bmat.Inverse[kInt].Transpose[bMat]` directly, solve the saddle-point linear system `kInt . y == Transpose[bMat] . z` for `y` symbolically in terms of `z`, then assemble `Bmat . y` and verify it equals `(Xi*I + alpha*v v^T) . z` for an arbitrary `z = {z1, z2}`. The intermediate variable names should be distinct from the SymPy ones (already partially true) and the algebraic path should differ.
- Stage 15.4: derive `lambda_minus` and `s_minus` as eigenvalue / eigenvector-norm components of an explicit 2x2 source-coupling matrix `M = {{A, alpha0*kappa0*sqrt(Kprod)/...}, ...}` rather than as bare formulas. Use Mathematica's `Eigenvalues` / `Eigensystem` and verify the smaller eigenvalue and the relevant component-norm match the asserted closed form. This converts Stage 15.4 from "state a formula and check two limits" into "derive the eigenvalue from a stated matrix and check it agrees".

If a full independent rewrite is too large, at minimum re-derive Stage 15.3 via the `LinearSolve` route (it is the largest single block of choreography copying) and document the methodological choice in a header comment.

**Verification:**
The verifier inspects the `.wl` and confirms that (i) Stage 15.3 no longer uses `Inverse[kInt]` directly as the principal Schur computation (it instead uses `LinearSolve` or `Eigenvalues` for the equivalent assertion), (ii) Stage 15.4 introduces and uses a `2x2` matrix and Mathematica's `Eigenvalues`/`Eigensystem` to obtain `lamMinus` / `sMinus`, and (iii) the script still exits 0 and matches the SymPy numeric / symbolic results at the assertion level.

### F3 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py:145-167`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl:117-132`

**What's wrong:**
Stage 15.4 ("NATURAL D/N SOURCE MAP") asserts two things:
```
expect_zero("mhat_-^2(alpha=0) - 1", ...)                    # line 166
expect_zero("limit_{alpha->oo} mhat_-^2 - 11/9", ...)         # line 167
```
But the formula being checked, `s_minus = (1/2) * (sigma + ((DK + alpha0 delta_kappa) delta_kappa + 4 alpha0 Kprod)/R)`, is written down at line 155 as a bare expression. It is not derived from anything within this script. The two assertions only verify that `s_minus/kappa0**2` has the right values at the two endpoints `alpha0 = 0` and `alpha0 = +oo`. The interior of the parameter space is never exercised, and the *structural* claim (that `s_minus` is the smaller eigenvalue / smaller-branch source component of some 2x2 system) is not asserted.

Specifically:
- The assumption `DK > 0, alpha0 >= 0` (line 147) means `R(alpha=0) = DK`, so the `alpha=0` check reduces algebraically to `(sigma + delta_kappa)/(2 kappa0**2) = (2 kappa0**2)/(2 kappa0**2) = 1` — a one-line algebraic identity given `subs_nat`.
- The `alpha->oo` check reduces (after using `sigma = kappa0**2 + kappa1**2`, `delta_kappa = kappa0**2 - kappa1**2`, `Kprod = kappa0**2 kappa1**2`, so `delta_kappa**2 + 4 Kprod = sigma**2`) to `s_minus -> sigma/2 + sigma/2 = sigma = 11/9 * kappa0**2`, also a one-line identity.

Neither check exercises the "selected source map" claim at any interior `alpha0`, nor verifies that `s_minus` is actually the eigenvalue/eigenvector of the matrix the docstring implies it diagonalises.

**Why this matters:**
If the stated closed form of `s_minus` had an interior error (e.g. wrong coefficient on the cross term `((DK + alpha0 delta_kappa) delta_kappa + 4 alpha0 Kprod)`), the two endpoint checks would still pass — at `alpha0 = 0` the entire problematic term vanishes (multiplied by `alpha0`), and at `alpha0 -> oo` the closed form is dominated by the `alpha0*sigma`/`R` ratio which is determined by the leading symbols not the cross term coefficient. The interior of the parameter space is precisely where a structural error would show, and it is not exercised.

**Required change:**
Add an interior consistency check. Specifically, derive `lambda_minus` and `s_minus` from an explicit 2x2 matrix `M(alpha0)` whose smaller eigenvalue is `lambda_minus` and a specified component-product or quadratic form is `s_minus`, and assert that:
1. `det(M - lambda_minus * I) == 0` symbolically (so `lambda_minus` is genuinely an eigenvalue), AND
2. The relationship asserted between `s_minus` and the eigenvector of `M` corresponding to `lambda_minus` holds.

The explicit form of `M` should be derivable from Stages 15.2/15.3 already in the script — specifically, from the source-map block (`gQ Q (v.q)` term in Stage 15.2 combined with the kappa-weighted 2x2 sub-block of the Schur problem in Stage 15.3). If Codex cannot reconstruct `M` from the script's existing content, it should ADD a comment block describing where `lambda_minus` and `s_minus` come from (cite line numbers) and add at minimum one interior-point numerical check at `alpha0 = 1, DK = 1, A = 1` confirming `s_minus_nat` matches a separately-computed numerical eigenvalue/eigenvector of an explicit `M` matrix instantiated with `(kappa0, kappa1) = (2*sqrt(2)/pi, -4/(3*pi))`.

Mirror the new interior check in the Mathematica script.

**Verification:**
After Codex applies, the SymPy script's Stage 15.4 contains at least one new `expect_zero` call evaluated at an interior `alpha0 > 0` (numeric or symbolic, but not at 0 and not as a limit). Ideally there is also an explicit `M = sp.Matrix(...)` construction and a `M.eigenvals()` or `det(M - lambda_minus*I)` check. The Mathematica script gets a matching `Eigenvalues[M]` or `Det[M - lamMinus*IdentityMatrix[2]]` check.

## Independent-derivation check (Mathematica)

The Mathematica script does NOT independently derive any of the unit's claims. Every stage is a line-by-line mirror of the SymPy script, as documented in F2 above with three side-by-side quotations (Stage 15.3 matrix construction, Stage 15.4 formula assignments, Stage 15.5 product structure). The only methodological independence is in `FullSimplify` versus `simplify`, which is engine-level rather than derivation-level.

This becomes the `mathematica_transliteration` finding F2.

## Engine cross-check

Both engines run to completion and report PASS on the same assertion strings. Selected side-by-side checks of the actual numeric / symbolic output:

| Quantity | SymPy output | Mathematica output | Agree? |
|---|---|---|---|
| `kappa0` | `2*sqrt(2)/pi` (line 17) | `(2*Sqrt[2])/Pi` (line 17) | yes |
| `kappa1` | `-4/(3*pi)` (line 18) | `-4/(3*Pi)` (line 18) | yes |
| `sigma` | `88/(9*pi**2)` (line 19) | `88/(9*Pi^2)` (line 19) | yes |
| `Xi` | `g_U**2/A_U` (line 103) | `gU^2/aU` (line 61) | yes |
| `alpha` | `g_B**2/A_phi + (A_U*g_W + g_R*g_U)**2/(A_U*(A_U*A_W - 88*g_R**2/(9*pi**2)))` (line 104) | `gB^2/aPhi + (gR*gU + aU*gW)^2/(aU*(aU*aW - (88*gR^2)/(9*Pi^2)))` (line 62) | yes |
| `mhat_-^2` | `(63*pi**2*DeltaK + 968*alpha0 + 11*sqrt(4608*alpha0**2 + (9*pi**2*DeltaK + 56*alpha0)**2))/(18*sqrt(4608*alpha0**2 + (9*pi**2*DeltaK + 56*alpha0)**2))` (line 113) | `(11 + (968*alpha0 + 63*dK*Pi^2)/Sqrt[7744*alpha0^2 + 1008*alpha0*dK*Pi^2 + 81*dK^2*Pi^4])/18` (line 69) | yes (algebraically equivalent: 56² = 3136, 56² + 4608 = 7744; 2·56·9 = 1008; 81 = 9²) |

No `engine_disagreement` finding. The engines do agree — but, because the Mathematica script is a transliteration (F2), agreement is not informative.

## Verdict justification

Three findings: one high-severity tautology in Stage 15.5 (F1), one medium-severity transliteration covering the entire Mathematica script (F2), and one medium-severity insufficient-verification in Stage 15.4 (F3). Stages 15.1, 15.2, and 15.3 hold up to attack — the integrals and Schur algebra are genuinely exercised. I tried to break Stage 15.1 by recomputing `kappa1 = ∫(2/L) cos(πs/L) sin(πs/(2L)) ds from 0 to L` via product-to-sum: `cos A sin B = (1/2)[sin(A+B) - sin(A-B)]`, giving `(1/L)[-2L/(3π) cos(3πs/(2L)) + 2L/π cos(πs/(2L))]_0^L = (1/L)(0 - 4L/(3π)) = -4/(3π)`, which matches the script. I tried to break Stage 15.3's Schur claim by checking the bottom-right entry of `Sigma - [Xi I + alpha vv^T]` for a missing `kappa1²` term; the SymPy zero-matrix output confirms cancellation. Stages 15.4 and 15.5, however, are weak: 15.4 only checks endpoints and never derives `s_minus`; 15.5 checks an identity that holds for any choice of variables. Verdict: `findings`. Not stop-cold: every finding is mechanically fixable by adding substantive checks, and none of the corrections would propagate a sign flip or constant change to downstream units (the closed forms themselves are not changed — only how the script verifies them).

## Self-test notes

For F1 I substituted line 174 into line 175 mentally and confirmed `Nprod - Nprod_target` is literally `(a/b)*(c/d) - ac/(bd) = 0` by basic arithmetic; no simplifier needed. For F3 I verified my proposed interior check is not the same as the endpoint checks by computing `s_minus_nat` at `alpha0 = 1, DK = 1` symbolically: it does not reduce to `kappa0**2` (alpha0=0 value) nor to `sigma` (alpha0→∞ value), so the interior check would have non-trivial content. For F2 I confirmed the three quoted Mathematica blocks correspond exactly to SymPy blocks at the same stage with the same operand order. The directive's "Target" paths use `scripts/...py` and `mathematica/...wl` as appropriate per the script-location convention.
