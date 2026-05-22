---
unit_id: 045
batch: III.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: present
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 045 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.txt`

mtime summary: `.py` Apr 3 12:13, `.py.txt` May 11 12:43 (fresh); `.wl` May 11 11:56, `.wl.txt` May 11 12:50 (fresh).

## What the script claims to verify

The script's docstring lists four claims for the "Stage 28" coherent local D/N support kernel: (1) the kernel implies the cross-product identity `g_B g_R = g_W g_S` for mass-normalized amplitudes; (2) the mixed and support interference ratios coincide, `rho_0 = sigma_0`; (3) the common direction factor `R_tr = (1 + chi_0/(1+delta_U))/(1 + chi_0)` satisfies two range identities; and (4) the Stage-27 quadratic branch equation collapses, after setting `R_phi = R_U`, to a one-parameter tracking law of the form `M_tr = xi*(delta+xi)/(delta + (1+lambda0*R_U^2)*xi)`, with a D/N specialization at `lambda0 = 2/9` and an `F_tr` normalization expression. The Mathematica script declares the same scope ("STAGE 028 — COHERENT LOCAL TRACKING") and asserts the same identities.

## Assertion inventory

| #  | Script      | Line | Form                                                                                    | Anchored to claim? |
|----|-------------|------|-----------------------------------------------------------------------------------------|--------------------|
| A1 | sympy       | 51   | `expect_zero("g_B g_R - g_W g_S", g_B*g_R - g_W*g_S)`                                   | no (tautology)     |
| A2 | sympy       | 57   | `expect_zero("rho_0 - sigma_0", rho_0 - sigma_0)`                                       | no (tautology)     |
| A3 | sympy       | 77-80| `expect_zero("range identity 1", expr1 - chi_0*delta_U/((1+chi_0)*(1+delta_U)))`        | partial (rewrite)  |
| A4 | sympy       | 81-84| `expect_zero("range identity 2", expr2 - delta_U/((1+chi_0)*(1+delta_U)))`              | partial (rewrite)  |
| A5 | sympy       | 108  | `expect_zero("M_tr - expected", M_tr - M_tr_expected)`                                  | no (tautology)     |
| A6 | sympy       | 139-142| `expect_zero("tracking quadratic collapse", num_track + collapsed_num.subs(...))`     | yes                |
| A7 | sympy       | 146-149| `expect_zero("G_tr formula", M_tr_req - xi*(delta+xi)/(delta+(1+lam0*R_U^2)*xi))`     | yes (solve-based)  |
| A8 | sympy       | 155  | `expect_zero("G_tr D/N specialization", G_tr_dn - G_tr_expected)`                       | partial (substitution rewrite) |
| A9 | sympy       | 167  | `expect_zero("F_tr normalization law", F_track.subs(lam0,2/9) - F_tr_expected)`         | partial (algebraic rewrite) |
| B1 | mathematica | 43   | `expectZero["g_B g_R - g_W g_S", gB*gR - gW*gS]`                                        | no (tautology)     |
| B2 | mathematica | 46   | `expectZero["rho_0 - sigma_0", rho0 - sigma0]`                                          | no (tautology)     |
| B3 | mathematica | 49-52| `expectZero["1 - R_tr - chi0 deltaU/((1+chi0)(1+deltaU))", ...]`                        | partial (rewrite)  |
| B4 | mathematica | 53-56| `expectZero["R_tr - 1/(1+deltaU) - deltaU/((1+chi0)(1+deltaU))", ...]`                  | partial (rewrite)  |
| B5 | mathematica | 75   | `expectZero["M_tr - expected", mTr - mTrExpected]`                                      | no (tautology)     |
| B6 | mathematica | 94   | `expectZero["tracking quadratic collapse", numTrack + (collapsedNum /. mTrSym -> mMixSym+mSuppSym)]` | yes |
| B7 | mathematica | 96-98| `expectZero["G_tr generic formula", mTrReq - xi*(delta+xi)/(delta+(1+lambda0*rU^2)*xi)]` | NO (self-comparison) |
| B8 | mathematica | 103  | `expectZero["G_tr D/N specialization", gTrDN - gTrExpected]`                            | partial (substitution rewrite) |
| B9 | mathematica | 115  | `expectZero["F_tr normalization law", FullSimplify[fTrack /. lambda0 -> 2/9] - fTrExpected]` | partial (algebraic rewrite) |

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:46-51` (g_B g_R = g_W g_S)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:53-57` (rho_0 = sigma_0)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:97-108` (M_tr = expected)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:34-46` (mirrors)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:64-75` (mirror M_tr identity)

**What's wrong:**

(a) Lines 46-51 of the SymPy script define
```
g_W = lam_W / sqrt(mu_eta * mu_W)
g_R = gamma * lam_W / sqrt(mu_U * mu_W)
g_B = lam_phi / sqrt(mu_eta * mu_phi)
g_S = gamma * lam_phi / sqrt(mu_U * mu_phi)
```
and then asserts `g_B*g_R - g_W*g_S == 0`. Both products expand to `gamma*lam_W*lam_phi/sqrt(mu_eta*mu_U*mu_W*mu_phi)` by commutativity of multiplication. The assertion cannot fail no matter what physics motivates the definitions — there is no symbol the author could change in the definitions of `g_W,g_R,g_B,g_S` (without retyping the formulas) that would make the cross-product fail. The claim in the docstring "the coherent local D/N support kernel implies `g_B g_R = g_W g_S` exactly" purports to test a non-trivial structural consequence of the kernel; the assertion does not — it tests `a*b*c*d == a*b*c*d`.

(b) Line 53-57: `rho_0 = g_R*g_U/(K_U*g_W)` and `sigma_0 = g_U*g_S/(K_U*g_B)`. Substituting the definitions gives `rho_0 = gamma*g_U*sqrt(mu_eta/mu_U)/K_U` and `sigma_0 = gamma*g_U*sqrt(mu_eta/mu_U)/K_U`. Both factor through the same cancellation as (a); the assertion is a direct corollary of (a)'s tautology.

(c) Line 99 sets `M_tr = simplify(M_mix + M_supp)`. Line 104 sets `M_tr_expected = simplify(8*(1+chi_0)^2/(pi^2*(1-eps_eta)) * (Z_W/(1-eps_W_split) + Z_phi/(1-eps_phi_split)))`. M_mix and M_supp share the common factor `8*(1+chi_0)^2/(pi^2*(1-eps_eta))`; M_tr_expected is precisely M_mix + M_supp with that factor pulled out. The assertion `M_tr - M_tr_expected == 0` is the distributive law, not a physical claim.

(d) The Mathematica script reproduces (a), (b), (c) verbatim (lines 34-46 and 64-75) and inherits the same tautologies.

**Why this matters:**

Three of the four docstring claims are tested by assertions that are algebraically guaranteed by the assignments above them. They contribute zero adversarial signal: any error introduced into the *physics* upstream of these definitions (e.g., a wrong scaling, a wrong coupling structure) would not be caught because the script re-introduces the same parameterization on both sides of the equality. The script therefore over-states its verification scope.

**Required change:**

For each of the three identities, rewrite the assertion so that the *left* side is constructed from an independent route (e.g., from a physical bilinear coupling expression or a kernel commutator) rather than the same factorization used on the right. Concretely:

- (a) Define `coupling_density = (lam_W*W_sym + lam_phi*phi_sym) * (eta_sym - gamma*U_sym)`. Extract the four bilinear coefficients in front of `W*eta`, `W*U`, `phi*eta`, `phi*U` from `coupling_density` (e.g., via `sp.collect` / `sp.Poly` over `{W_sym, phi_sym, eta_sym, U_sym}`), normalize each by `1/sqrt(mu_i*mu_j)`, and *then* assert the cross-product identity. The route through the polynomial extraction makes the identity contingent on the bilinear structure rather than on retyping the four amplitudes.
- (b) Express `rho_0` and `sigma_0` directly in terms of those polynomial-extracted amplitudes; assert the equality after substitution rather than after manually inlining the closed form.
- (c) Replace `M_tr_expected` with an *independent symbolic form* derived from the bilinear coupling above (e.g., apply the same polynomial extraction to the `(Z_W*W + Z_phi*phi)` term and the splittings `1-eps_*`). Do not write the answer down by hand.

If the author cannot supply an independent derivation route for any of these three, the offending assertion must be removed and the corresponding docstring bullet deleted (do not advertise verification that the script does not perform).

**Verification:**

After the fix:
- The new SymPy script should contain a `coupling_density = ...` expression at the start of section 1, a `sp.collect`/`sp.Poly` extraction step, and an assertion `expect_zero("g_B*g_R - g_W*g_S", g_B_extracted*g_R_extracted - g_W_extracted*g_S_extracted)` where the four `_extracted` quantities trace back to coefficients of `coupling_density` *only* (no manual reassignment).
- The new M_tr_expected should be a `sp.simplify` of a different algebraic form that is *not* literally `M_mix + M_supp` written as `A*(B/C + D/E)`. One acceptable form: derive M_tr by summing over an enumerated coupling channel list (W, phi) with per-channel weights.
- Mathematica script mirrors the new structure (without re-using SymPy intermediate names).
- New `.txt` output should show explicit intermediate "extracted coefficient" prints.

### F2 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:96-98`

**What's wrong:**

The SymPy script genuinely solves a quadratic at line 144:
```
M_tr_req = sp.simplify(sp.solve(sp.Eq(collapsed_num, 0), M_tr_sym)[0])
```
and then at line 146-149 asserts the solved expression equals the claimed closed form. That is a substantive check: had `collapsed_num` been written with a wrong coefficient, `sp.solve` would have returned a different expression and the assertion would fail.

The Mathematica script *omits* the solve step entirely:
```
mTrReq = FullSimplify[xi*(delta + xi)/(delta + (1 + lambda0*rU^2)*xi), Assumptions -> $Assumptions];
Print["M_tr required on tracking branch = ", fmt[mTrReq]];
expectZero["G_tr generic formula", mTrReq - xi*(delta + xi)/(delta + (1 + lambda0*rU^2)*xi)];
```
`mTrReq` is *defined* as the target expression, then checked against the *same* expression. The residual is `0` literally because the LHS and RHS are syntactically identical. This is a pure self-comparison.

**Why this matters:**

The Mathematica script's "G_tr generic formula" PASS line in the output is meaningless. The second-engine verification at this step (claim 4 in the docstring) is missing on the Mathematica side. If the SymPy `sp.solve` step were broken (wrong root selected, wrong sign convention), the Mathematica script would silently still pass, defeating the cross-engine purpose.

**Required change:**

In the WL file at line 96, replace the direct definition of `mTrReq` with an actual `Solve` (or `Reduce`) call against `collapsedNum`:

Before (lines 96-98):
```mathematica
mTrReq = FullSimplify[xi*(delta + xi)/(delta + (1 + lambda0*rU^2)*xi), Assumptions -> $Assumptions];
Print["M_tr required on tracking branch = ", fmt[mTrReq]];
expectZero["G_tr generic formula", mTrReq - xi*(delta + xi)/(delta + (1 + lambda0*rU^2)*xi)];
```

After:
```mathematica
mTrReqSolutions = Solve[collapsedNum == 0, mTrSym];
mTrReq = FullSimplify[mTrSym /. First[mTrReqSolutions], Assumptions -> $Assumptions];
Print["M_tr required on tracking branch = ", fmt[mTrReq]];
expectZero["G_tr generic formula", mTrReq - xi*(delta + xi)/(delta + (1 + lambda0*rU^2)*xi)];
```

`collapsedNum` is already defined at line 90 as `Expand[xi^2 + (delta - mTrSym*(1 + lambda0*rU^2))*xi - delta*mTrSym]`, which is linear in `mTrSym`, so `Solve` returns a unique root — no branch ambiguity.

**Verification:**

After the fix, the `.wl` output should show `M_tr required on tracking branch = (xi*(delta + xi))/(delta + xi + lambda0*rU^2*xi)` (or equivalent normalization), and the residual line `G_tr generic formula = 0` must remain. A new `Solve` invocation should appear in the script source at the location of the former direct definition; the verifier should be able to grep for `Solve[collapsedNum` (or `Solve[collapsedNum ==`) in the source.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:34-115`

**What's wrong:**

The Mathematica script is a near-mechanical port of the SymPy script. Concrete correspondences:

(i) Sections 1 (amplitudes), lines 34-37 of the WL vs lines 46-49 of the PY:
```mathematica
gW = lamW/Sqrt[muEta*muW];
gR = gamma*lamW/Sqrt[muU*muW];
gB = lamPhi/Sqrt[muEta*muPhi];
gS = gamma*lamPhi/Sqrt[muU*muPhi];
```
vs
```python
g_W = lam_W / sp.sqrt(mu_eta * mu_W)
g_R = gamma * lam_W / sp.sqrt(mu_U * mu_W)
g_B = lam_phi / sp.sqrt(mu_eta * mu_phi)
g_S = gamma * lam_phi / sp.sqrt(mu_U * mu_phi)
```
identical variable choreography in identical order.

(ii) Section 4 (branch equation), lines 82-87 of WL vs lines 121-129 of PY:
```mathematica
branchEq = FullSimplify[
  mSuppSym - (xi*(delta + xi) - mMixSym*(delta + (1 + lambda0*rU^2)*xi))/
    (delta + (1 + lambda0*rPhi^2)*xi - mMixSym*lambda0*(rU - rPhi)^2),
  Assumptions -> $Assumptions
];
branchTrack = Together[FullSimplify[branchEq /. rPhi -> rU, Assumptions -> $Assumptions]];
```
vs
```python
branch_eq = sp.simplify(
    Msupp - (
        xi * (delta + xi) - Mmix * (delta + (1 + lam0 * R_U ** 2) * xi)
    ) / (
        delta + (1 + lam0 * R_phi ** 2) * xi - Mmix * lam0 * (R_U - R_phi) ** 2
    )
)
branch_track = sp.together(sp.simplify(branch_eq.subs(R_phi, R_U)))
```
identical rational form, identical R_phi -> R_U substitution, identical `together` pipeline.

(iii) Section 4 closing (F_track), lines 105-115 of WL vs lines 157-167 of PY:
```mathematica
fTrack = FullSimplify[
  (delta + (1 + lambda0*rU^2)*xi)^2*(delta + (1 + lambda0*rU)*xi)^2/
    ((1 - xi)*((delta + xi)^2 + lambda0*rU^2*xi^2)^2), ...
```
vs
```python
F_track = sp.simplify(
    (delta + (1 + lam0 * R_U ** 2) * xi) ** 2
    * (delta + (1 + lam0 * R_U) * xi) ** 2
    / ((1 - xi) * ((delta + xi) ** 2 + lam0 * R_U ** 2 * xi ** 2) ** 2)
)
```
literal expression copy.

The Mathematica script chooses no different intermediate steps, no different parameter route, no different simplification order. It contributes no independent algebra.

**Why this matters:**

The second-engine policy requires that the Mathematica script derive the same physical conclusion via a different symbolic route, so that a bug in either engine's algebra (factoring, branch choice, simplification under assumptions) is detected by disagreement. A line-by-line port detects only typo-class errors; it does not detect systematic algebraic mistakes shared by both versions.

**Required change:**

Refactor the WL script (without touching the SymPy script) so that at least one of sections 1, 3, 4 derives its result through a *different* algebraic route. Suggested route changes:

- **Section 1**: Instead of writing the four amplitudes explicitly, introduce a bilinear form `bilinear = (lamW*W + lamPhi*Phi)*(eta - gamma*U)` and use `Coefficient[Expand[bilinear], W*eta]`, `Coefficient[Expand[bilinear], W*U]`, etc., then divide by the appropriate `Sqrt[muI*muJ]`. This makes the WL section structurally distinct from the PY section while testing the same identity.
- **Section 4**: Instead of substituting `rPhi -> rU` first and then taking the numerator, take the numerator of `branchEq` first (as a polynomial in `rPhi - rU`), then expand around `rPhi = rU` (`Series[..., {rPhi, rU, 1}]` and read off the zeroth-order term). Verify it matches the SymPy-side `num_track` by writing the same residual assertion.

The directive only needs to require *one* such route change to qualify as independent; the auditor's bar is that the WL script must not be a token-for-token rewrite.

**Verification:**

After refactor, the WL script's section 1 (or section 4) should contain at least one structural construct that does not appear in the corresponding PY section (e.g., `Coefficient[Expand[...]]`, `Series[..., {rPhi, rU, 1}]`, or `GroebnerBasis[...]`). The auditor will re-grep both scripts and confirm that the refactored section's intermediate algebra is distinct.

## Independent-derivation check (Mathematica)

The `.wl` does not derive the claims independently. As laid out under F3, the variable definitions, intermediate substitutions, and even the simplification calls match the `.py` line by line. Section 4 in particular reuses the SymPy `branch_eq.subs(R_phi, R_U)` pattern exactly, and the F_track expression is a literal copy with mass-symbol renaming (`lam0 -> lambda0`, `R_U -> rU`). This is a `mathematica_transliteration` finding (F3).

## Engine cross-check

Both engines emit residual `= 0` on every assertion that the script presents as a verification (see the matched residual lines in both `.txt` outputs):
- `g_B g_R - g_W g_S = 0` (both)
- `rho_0 - sigma_0 = 0` (both)
- range identities `= 0` (both)
- `M_tr - expected = 0` (both)
- `tracking quadratic collapse = 0` (both)
- `G_tr formula = 0` (both)
- `G_tr D/N specialization = 0` (both)
- `F_tr normalization law = 0` (both)

The intermediate displayed forms also match modulo Mathematica's preferred ordering: SymPy reports `R_tr = (chi_0 + delta_U + 1)/((chi_0 + 1)*(delta_U + 1))`; Mathematica reports `R_tr = (1 + chi0 + deltaU)/(1 + chi0 + deltaU + chi0*deltaU)` — algebraically identical after expanding the denominator. SymPy reports `M_tr required on tracking branch = xi*(delta + xi)/(R_U**2*lambda_0*xi + delta + xi)`; Mathematica reports `(xi*(delta + xi))/(delta + xi + lambda0*rU^2*xi)` — same expression.

No engine disagreement is detected, but per F2 the Mathematica side's `G_tr generic formula` PASS is a self-comparison, so the apparent agreement at that step is not informative.

## Verdict justification

The script does contain one substantive check (tracking quadratic collapse, A6/B6), but most of its advertised verifications — sections 1, 3, the range identities, the G_tr closed form on the WL side, and the F_tr normalization law — are tautological algebraic rewrites in which the expected expression is constructed from the same definitions as the computed expression. The Mathematica script compounds this by being a line-by-line port of the SymPy script and by self-comparing on the G_tr step. Findings F1, F2, F3 capture these issues. No engine disagreement is observed; freshness is fine. Verdict: `findings`, three findings, no stop-cold.

Attacks attempted that the script *did* survive:
- Tracking quadratic collapse (A6/B6): exercised genuinely via `branch_eq.subs(R_phi, R_U)` -> numerator -> identity with `xi^2 + (delta - M_tr*(1+lambda0*R_U^2))*xi - delta*M_tr`. Cannot be a tautology because `branch_eq` is a non-trivial rational and the collapsed quadratic structure could in principle disagree.
- Mass-positivity assumptions (`positive=True` on `mu_*, lam_*, gamma, K_U`, etc.): consistent with the kernel's physical setup; no contradiction with later substitutions (e.g., `lam0 = 2/9 > 0`).
- Symbol-domain attack on `delta_U > 0`: range identities use `1 + delta_U`, `chi_0 + 1`, both safely nonzero under `positive=True` — no division-by-zero hidden in the rewrites.

## Self-test notes

- Variable independence: no `sp.diff` or `D[...]` appears in this audit unit, so the "derivative w.r.t. an unused variable" trap does not apply. Skipped.
- Symmetry/parity: no integrals over unbounded domains here (no Gaussians, no `Integrate[..., {x, -Infinity, Infinity}]`); skipped.
- Trivial-case substitution: for F1's required change, mentally substituted `lam_W = lam_phi = mu_eta = mu_U = mu_W = mu_phi = 1` into the proposed polynomial-extraction route: `coupling_density = (W + phi)*(eta - gamma*U)` -> coefficient of `W*eta` is 1, of `W*U` is `-gamma`, of `phi*eta` is 1, of `phi*U` is `-gamma`. Normalized amplitudes become `1, gamma, 1, gamma`. Cross-product `g_B*g_R - g_W*g_S = 1*gamma - 1*gamma = 0` (still zero, as required), so the proposed independent route gives the correct answer at the trivial point. For F2's required change, mentally ran `Solve[xi^2 + (delta - mTrSym*(1+lambda0*rU^2))*xi - delta*mTrSym == 0, mTrSym]`: this is linear in `mTrSym`, root is `(xi^2 + delta*xi)/(delta + (1+lambda0*rU^2)*xi) = xi*(delta+xi)/(delta + (1+lambda0*rU^2)*xi)`, matching the target. No silent-pass risk.
- Path specification: F1, F2, F3 all target existing files; no new-script paths to write. The directive does not contain any `missing_verification_script` targets.
