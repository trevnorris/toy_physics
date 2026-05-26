---
unit_id: 028
batch: II.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
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
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage028_loaded_profile_selection.md
  paper_appendix: present
---

# Audit unit 028 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_028.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage028_loaded_profile_selection.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row at line 46, `\input` at line 94)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage028_loaded_profile_selection_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage028_loaded_profile_selection_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.txt`

Freshness: scripts mtime 2026-05-21 17:05; outputs mtime 2026-05-21 17:07. Both outputs are newer than their producing scripts. Fresh.

## What the paper claims

Quoting `\stagefield{Output}` verbatim: "Stage~028 outputs the rank-one wall matrix \eqref{eq:app-stage028-Keff}, exact eigenvalues \eqref{eq:app-stage028-eigenvalues}, the selected angle law \eqref{eq:app-stage028-angle-law}, and the softening threshold \eqref{eq:app-stage028-alpha-crit}." Concretely the four boxed deliverables are: (1) `K_eff = K_bare - alpha v v^T` with `v = (kappa0, kappa1)^T`, `kappa0 = 2 sqrt(2)/pi`, `kappa1 = -4/(3 pi)`, `K_bare = diag(K_0, K_1)`, `K_0 = K_eta + 6 T_Omega`, `K_1 = K_0 + DeltaK_ax`, `DeltaK_ax = pi^2 T_w/L^2`; (2) the exact eigenvalues `lambda_± = (1/2)[K_0 + K_1 - alpha(kappa0^2 + kappa1^2) ± sqrt((DeltaK_ax + alpha(kappa0^2 - kappa1^2))^2 + 4 alpha^2 kappa0^2 kappa1^2)]`; (3) the selected angle law `tan(2 theta_-) = 2 alpha kappa0 kappa1 / (DeltaK_ax + alpha(kappa0^2 - kappa1^2))`; (4) `alpha_crit = K_0 K_1 / (K_1 kappa0^2 + K_0 kappa1^2)` from `det K_eff = 0`. The `\stagefield{Checks}` block adds (a) sign of theta_- fixed by kappa0>0 and kappa1<0, and (b) the blind angle cannot be selected under positive loading. The notes (§§4.1, 4.2, 5) elaborate two consequence-level checks the script is expected to exercise: weak-loading expansion theta = alpha kappa0 kappa1/DeltaK_ax + O(alpha^2), and strong-loading limit tan(theta_max) = kappa1/kappa0 = -sqrt(2)/3.

## What the script claims to verify

Both engines construct `K_eff = K_bare - alpha v v^T` symbolically using the paper's literal overlap constants (`kappa0 = 2 sqrt(2)/pi`, `kappa1 = -4/(3 pi)`, `DeltaK_ax = pi^2 T_w/L^2`, `K_0 = K_eta + 6 T_Omega`, `K_1 = K_0 + DeltaK_ax`). They then assert: (i) trace and determinant of K_eff match the closed-form expressions; (ii) the characteristic polynomial factorizes through the constructed `lambda_±` (and, in Mathematica, that `Eigenvalues[kEff]` sum/product agree as an independent eigensolver path); (iii) `d/dtheta (q^T K_eff q / 2) = (DeltaK_ax + alpha(kappa0^2 - kappa1^2)) sin(2 theta)/2 - alpha kappa0 kappa1 cos(2 theta)`, recovering the tan(2 theta) relation; (iv) the kappa0 and kappa1 literals have the positive and negative signs the paper claims; (v) the weak-alpha leading coefficient of theta is `kappa0 kappa1/DeltaK_ax`; (vi) the strong-alpha limit of tan(2 theta) equals `tan(2 theta_max)` with `tan(theta_max) = -sqrt(2)/3`; (vii) `alpha_crit` obtained from `Solve[det K_eff = 0, alpha]` matches the paper's closed-form ratio and (Mathematica only) a re-simplified `9 pi^2 K_0 K_1 / (8(11 K_0 + 9 DeltaK_ax))`; and (viii) `det(K_eff)` evaluated at `alpha_crit*(1-eps)` factors as a positive multiple of `eps`, demonstrating that the determinant approaches zero linearly from above as `alpha` approaches `alpha_crit` from below.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) K_eff = K_bare - alpha vv^T with paper's kappa0, kappa1 | sympy L74-76 builds it with the literal kappa0=2 sqrt(2)/pi, kappa1=-4/(3 pi); math L33-41 mirrors | match |
| (2) Closed-form eigenvalues lambda_± | sympy L112-119 constructs and confirms via characteristic factorization; math L60-76 adds independent `Eigenvalues[kEff]` sum/product cross-check | match |
| (3) tan(2 theta_-) angle law | sympy L138-148 (derived from stationarity); math L81-94 (mirrors) | match |
| (4) alpha_crit = K0 K1 / (K1 kappa0^2 + K0 kappa1^2) | sympy L192-200 solves det=0 and matches expected; math L115-124 solves and matches both ratio and simplified closed form | match |
| Checks (a) sign of theta_- < 0 driven by kappa0 k1 < 0 | sympy L152-153 and math L97-98 assert kappa0 > 0 and kappa1 < 0 via `kappa0 - Abs[kappa0] == 0`, `kappa1 + Abs[kappa1] == 0`; sign of the rhs of tan(2 theta_-) thereby becomes manifest. Explicit `theta_- < 0` is not asserted, but the sign of the numerator/denominator is asserted via literal expansion: `2 kappa0 kappa1 = -16 sqrt(2)/(3 pi^2)` and `kappa0^2 - kappa1^2 = 56/(9 pi^2)` are printed (sympy L87-88, math L47-48). | partial-but-acceptable (paper claim is a sign consequence, not an independent identity) |
| Checks (b) blind angle dynamically disfavored | Indirectly: sympy L179-180 / math L112-113 prove strong-loading limit aligns the selected eigenvector with `v/|v|` (i.e., `tan(theta_max) = -sqrt(2)/3`), which is opposite-sign to `tan(theta_blind) = 3 sqrt(2)/2 > 0`. No explicit comparison line, but the strong-limit assertion plus the kappa0, kappa1 sign asserts establish the implication. | partial-but-acceptable |
| Notes §4.1 weak-loading coefficient | sympy L156-158, math L100-102: `weak_coeff - kappa0 kappa1/DeltaK == 0` | match |
| Notes §4.2 strong-loading limit | sympy L171-180, math L104-113: `strong_limit - tan(2 theta_max) == 0` AND `tan(theta_max) + sqrt(2)/3 == 0` | match |

`paper_alignment: aligned` — every load-bearing `\stagefield{Output}` deliverable has a hard-zero assertion against the paper's literal form; the two "partial-but-acceptable" rows verify consequences of asserted facts rather than re-asserting the consequence.

## Assertion inventory

| #   | Script      | Line  | Form                                                                                  | Exercises which paper claim?               | Anchored to claim? |
|-----|-------------|-------|---------------------------------------------------------------------------------------|--------------------------------------------|--------------------|
| A1  | sympy       | 109   | `expect_zero("trace - expected", tr_eff - tr_expected)`                                | (2) prerequisite (trace of K_eff)          | yes                |
| A2  | sympy       | 110   | `expect_zero("det - expected", det_eff - det_expected)`                                | (4) prerequisite (det of K_eff)            | yes                |
| A3  | sympy       | 119   | `expect_zero("characteristic-factorization check", char)`                              | (2) eigenvalues closed form                | yes                |
| A4  | sympy       | 146   | `expect_zero("dE/dtheta - stationarity_expected/2", ...)`                              | (3) angle law                              | yes                |
| A5  | sympy       | 152   | `expect_zero("kappa0 sign check", kappa0 - Abs(kappa0))`                               | (1) overlap literal sign / Checks (a)      | yes                |
| A6  | sympy       | 153   | `expect_zero("kappa1 sign check", kappa1 + Abs(kappa1))`                               | (1) overlap literal sign / Checks (a)      | yes                |
| A7  | sympy       | 158   | `expect_zero("weak-loading coefficient - k0*k1/DeltaK", ...)`                          | Notes §4.1                                 | yes                |
| A8  | sympy       | 179   | `expect_zero("strong-loading limit - tan(2 theta_max)", ...)`                          | Notes §4.2 / Checks (b)                    | yes                |
| A9  | sympy       | 180   | `expect_zero("tan(theta_max) + sqrt(2)/3", ...)`                                       | Notes §4.2 literal                         | yes                |
| A10 | sympy       | 200   | `expect_zero("alpha_crit - expected", ...)`                                            | (4) alpha_crit                             | yes                |
| A11 | sympy       | 207   | `expect_zero("det(alpha_crit_expected)", ...)`                                         | (4) det zero at threshold                  | yes                |
| M1  | mathematica | 57    | `expectZero["trace - expected", trEff - trExpected]`                                   | (2) prerequisite                           | yes                |
| M2  | mathematica | 58    | `expectZero["det - expected", detEff - detExpected]`                                   | (4) prerequisite                           | yes                |
| M3  | mathematica | 67-70 | `expectZero["Eigenvalues[kEff] sum vs trace", Total[eigvalsDirect] - (lambdaMinus + lambdaPlus)]` | (2) independent eigensolver path | yes                |
| M4  | mathematica | 71-74 | `expectZero["Eigenvalues[kEff] product vs determinant", Times @@ eigvalsDirect - lambdaMinus*lambdaPlus]` | (2) independent eigensolver path | yes      |
| M5  | mathematica | 76    | `expectZero["characteristic factorization", charResidual]`                             | (2) eigenvalues closed form                | yes                |
| M6  | mathematica | 93    | `expectZero["dE/dtheta - stationarity_expected/2", ...]`                               | (3) angle law                              | yes                |
| M7  | mathematica | 97    | `expectZero["kappa0 sign check (kappa0 > 0)", kappa0 - Abs[kappa0]]`                   | (1) sign / Checks (a)                      | yes                |
| M8  | mathematica | 98    | `expectZero["kappa1 sign check (kappa1 < 0)", kappa1 + Abs[kappa1]]`                   | (1) sign / Checks (a)                      | yes                |
| M9  | mathematica | 102   | `expectZero["weak-loading coefficient - kappa0 kappa1/deltaK", ...]`                   | Notes §4.1                                 | yes                |
| M10 | mathematica | 112   | `expectZero["strong-loading limit - tan(2 theta_max)", ...]`                           | Notes §4.2                                 | yes                |
| M11 | mathematica | 113   | `expectZero["tan(theta_max) + Sqrt[2]/3", ...]`                                        | Notes §4.2 literal                         | yes                |
| M12 | mathematica | 116-119 | `expectZero["alphaCrit solved vs ratio closed form", Solve[detEff==0, alpha] - K0 K1/(K1 kappa0^2 + K0 kappa1^2)]` | (4) alpha_crit | yes                |
| M13 | mathematica | 123   | `expectZero["alpha_crit - finite-throat closed form", ...]`                            | (4) alpha_crit alternate form              | yes                |
| M14 | mathematica | 124   | `expectZero["det(alpha_crit)", detEff /. alpha -> alphaCrit]`                          | (4) det zero at threshold                  | yes                |

Every assertion traces to a paper-side deliverable. No orphan checks. All "yes" anchored.

## Findings

(no findings — see Verdict justification)

## Independent-derivation check (Mathematica)

The .wl is structurally parallel to the .py: same constants in the same order (kappa0, kappa1, deltaK, K0, K1); same intermediate quantities (trEff, detEff, disc, lambdaMinus, lambdaPlus, q, energy, dEnergy, rhs, weakCoefficient, strongLimit, tMax, tan2TMax, alphaCrit); same final `det(alpha_crit*(1-eps))` factorization. Compare sympy lines 102-114:

```
tr_eff = sp.simplify(sp.trace(K_eff))
det_eff = sp.simplify(K_eff.det())
tr_expected = sp.simplify(K0 + K1 - alpha * (kappa0**2 + kappa1**2))
det_expected = sp.simplify(K0 * K1 - alpha * (K1 * kappa0**2 + K0 * kappa1**2))
...
disc = sp.simplify((DeltaK + alpha * (kappa0**2 - kappa1**2))**2 + 4 * alpha**2 * kappa0**2 * kappa1**2)
lam_minus = sp.simplify((tr_expected - sp.sqrt(disc)) / 2)
lam_plus  = sp.simplify((tr_expected + sp.sqrt(disc)) / 2)
```

with math lines 50-65:

```
trEff = FullSimplify[Tr[kEff], ...];
detEff = FullSimplify[Det[kEff], ...];
trExpected = FullSimplify[K0 + K1 - alpha*(kappa0^2 + kappa1^2), ...];
detExpected = FullSimplify[K0*K1 - alpha*(K1*kappa0^2 + K0*kappa1^2), ...];
disc = FullSimplify[(deltaK + alpha*(kappa0^2 - kappa1^2))^2 + 4*alpha^2*kappa0^2*kappa1^2, ...];
lambdaMinus = FullSimplify[(trExpected - Sqrt[disc])/2, ...];
lambdaPlus  = FullSimplify[(trExpected + Sqrt[disc])/2, ...];
```

The angle equation (sympy 138-144 vs math 81-91) and alpha_crit solve (sympy 196-200 vs math 115-119) are also parallel.

However, the Mathematica file does add two genuinely independent steps that the SymPy file lacks:

1. Lines 66-74: direct `Eigenvalues[kEff]` then `Total[eigvalsDirect] - (lambdaMinus + lambdaPlus)` and `Times @@ eigvalsDirect - lambdaMinus*lambdaPlus`. Mathematica's built-in `Eigenvalues` is its own closed-form path; the assertion that its sum and product match the hand-built `lambdaMinus`, `lambdaPlus` is an independent confirmation of the discriminant identity.
2. Lines 120-123: `alphaCritClosed = 9*Pi^2*K0*K1/(8*(11*K0 + 9*deltaK))` cross-check against the ratio form (which the sympy script does not perform).

Because both engines also independently `Solve` for `alpha_crit` (sympy `sp.solve(sp.Eq(det_eff, 0), alpha)` line 196; math `Solve[detEff == 0, alpha]` line 115), and the Mathematica file uses its built-in eigensolver as a parallel check, the .wl is no longer a pure transliteration — it adds substantive Mathematica-specific paths. The structural parallelism in the trace/det/discriminant ledger is essentially forced by the underlying 2x2-with-rank-1-perturbation problem. I do not raise a `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree at the level of fully simplified residuals on every load-bearing identity. Key shared closed forms (compare sympy output lines 76, 80, 82-84, 89 with math output lines 26, 32, 35-37, 43):

| Quantity | Sympy | Mathematica |
|---|---|---|
| `tan(2 theta)` | `-48 sqrt(2) L^2 alpha / (56 L^2 alpha + 9 pi^4 T_w)` | `(-48 Sqrt[2] alpha L^2)/(56 alpha L^2 + 9 Pi^4 Tw)` |
| weak coeff | `-8 sqrt(2) L^2 / (3 pi^4 T_w)` (from `weak/alpha`) | `(-8 Sqrt[2] L^2)/(3 Pi^4 Tw)` |
| `lim_{alpha→∞} tan(2 theta)` | `-6 sqrt(2)/7` | `(-6 Sqrt[2])/7` |
| `tan(theta_max)` | `-sqrt(2)/3` | `-1/3 Sqrt[2]` |
| `alpha_crit` | `9 pi^2 (K_eta^2 L^2 + 12 K_eta L^2 T_Omega + pi^2 K_eta T_w + 36 L^2 T_Omega^2 + 6 pi^2 T_Omega T_w) / (8 (11 K_eta L^2 + 66 L^2 T_Omega + 9 pi^2 T_w))` (expanded) | `9 Pi^2 (Keta + 6 TOmega)(L^2 (Keta + 6 TOmega) + Pi^2 Tw) / (88 L^2 (Keta + 6 TOmega) + 72 Pi^2 Tw)` (factored) |
| `det(alpha_crit (1-eps))` | `eps (K_eta + 6 T_Omega)(K_eta L^2 + 6 L^2 T_Omega + pi^2 T_w)/L^2` | `eps (Keta + 6 TOmega)(Keta + 6 TOmega + Pi^2 Tw/L^2)` |

The `alpha_crit` forms differ only by expansion vs. factorization (numerator factors as `9 pi^2 (K_eta + 6 T_Omega)((K_eta + 6 T_Omega) L^2 + pi^2 T_w)`; denominator `88 L^2 (K_eta + 6 T_Omega) + 72 pi^2 T_w` matches sympy's `8 (11 K_eta L^2 + 66 L^2 T_Omega + 9 pi^2 T_w)`). The `det(alpha_crit (1-eps))` forms differ by an overall `1/L^2` distribution. Engines agree.

## Verdict justification

The two scripts together cover every numbered deliverable in `\stagefield{Output}` (K_eff, lambda_±, tan(2 theta_-) angle law, alpha_crit) and the supporting weak/strong limits the notes call out. Constants used are exactly the paper's literal values. I attacked the following:

1. **Characteristic-factorization tautology**: sympy A3 / math M5 check `(x - lam_-)(x - lam_+) - (x^2 - tr_eff x + det_eff) == 0`. By construction `lam_± = (tr_expected ± sqrt(disc))/2`, so the residual reduces to `(tr_expected^2 - disc)/4 - det_eff = 0`. Re-deriving by hand: `tr_expected^2 = (K0+K1)^2 - 2 alpha (K0+K1)(k0^2+k1^2) + alpha^2 (k0^2+k1^2)^2`, and `disc = (K1-K0)^2 + 2 alpha (K1-K0)(k0^2-k1^2) + alpha^2 (k0^2-k1^2)^2 + 4 alpha^2 k0^2 k1^2 = (K1-K0)^2 + 2 alpha (K1-K0)(k0^2-k1^2) + alpha^2 (k0^2+k1^2)^2`. So `(tr_expected^2 - disc)/4 = [(K0+K1)^2 - (K1-K0)^2]/4 - alpha [(K0+K1)(k0^2+k1^2) + (K1-K0)(k0^2-k1^2)]/2 = K0 K1 - alpha [K1 k0^2 + K0 k1^2] = det_expected`. This is the substantive content the check enforces (the `disc` ansatz really is `tr^2 - 4 det`); not tautological.

2. **Angle equation**: re-derived `(q^T K_eff q)/2` by hand: `(1/2)[K_0 cos^2 theta + K_1 sin^2 theta - alpha(kappa0 cos theta + kappa1 sin theta)^2]`. `dE/dtheta = (K_1 - K_0) sin theta cos theta - alpha (kappa0 cos theta + kappa1 sin theta)(-kappa0 sin theta + kappa1 cos theta)`. Expanding and using double-angle identities: `(1/2)(K_1 - K_0) sin(2 theta) - (alpha/2)[2 kappa0 kappa1 cos(2 theta) + (kappa1^2 - kappa0^2) sin(2 theta)] = (1/2)[(K_1 - K_0) + alpha(kappa0^2 - kappa1^2)] sin(2 theta) - alpha kappa0 kappa1 cos(2 theta) = (DeltaK_ax + alpha(kappa0^2 - kappa1^2))/2 sin(2 theta) - alpha kappa0 kappa1 cos(2 theta)`. Matches the script's `stationarity_expected/2`. Both stationary points (lower and upper eigenvectors, differing by pi/2 in theta) satisfy the same tan(2 theta) equation; the paper's `tan(2 theta_-)` notation is the equation for the angle of the lower eigenvector, and the script's derivation is consistent.

3. **Strong-loading limit**: `lim_{alpha->oo} -48 sqrt(2) L^2 alpha / (56 L^2 alpha + 9 pi^4 T_w) = -48 sqrt(2)/56 = -6 sqrt(2)/7`. Double-angle of theta_max with tan(theta_max) = -sqrt(2)/3: `tan(2 theta_max) = 2(-sqrt(2)/3)/(1 - 2/9) = (-2 sqrt(2)/3)(9/7) = -6 sqrt(2)/7`. Match.

4. **alpha_crit**: `K_1 kappa0^2 + K_0 kappa1^2 = (K_0 + DeltaK_ax)(8/pi^2) + K_0 (16/(9 pi^2)) = (8/(9 pi^2))(9 K_0 + 9 DeltaK_ax + 2 K_0) = (8/(9 pi^2))(11 K_0 + 9 DeltaK_ax)`. So `K_0 K_1 / (K_1 kappa0^2 + K_0 kappa1^2) = 9 pi^2 K_0 K_1 / (8 (11 K_0 + 9 DeltaK_ax))`. Matches Mathematica's `alphaCritClosed`. Both engines also `Solve` independently and compare.

5. **Sign of theta_-**: rhs of tan(2 theta_-) has numerator `2 alpha kappa0 kappa1 < 0` (since kappa0 > 0, kappa1 < 0, alpha > 0, all asserted) and denominator `DeltaK_ax + alpha(kappa0^2 - kappa1^2) > 0` (DeltaK_ax > 0 declared, kappa0^2 > kappa1^2 shown literally). So `tan(2 theta_-) < 0`, hence `theta_- < 0` modulo branch. Paper's sign claim is verified by the asserted sign data. The script does not literally `assert theta_- < 0`, but the constituent sign assertions are present and the implication is immediate.

6. **Blind-angle disfavor**: Strong-loading limit aligns with `tan(theta_max) = -sqrt(2)/3 < 0`, opposite-sign to `tan(theta_blind) = 3 sqrt(2)/2 > 0` (mentioned in Notes §5). The script asserts the strong-limit alignment and the sign of kappa1/kappa0, so the disjointness from the blind branch follows.

All attacks failed. Outputs are fresh. Engines agree to identical algebraic residuals. Paper and script align on every load-bearing identity. Verdict: clean.

One non-blocking cosmetic note (not a finding): the SymPy docstring header (line 3) and the Mathematica banner (line 26) both label the work "Stage 11" / "STAGE 011" rather than "Stage 028"; the final print in the Mathematica file does say "Stage 028 Mathematica audit passed" (line 130). This is a stage-renumbering artifact, not a math/paper disagreement (the content fully matches Stage 028), so I do not raise it as `paper_misalignment`.

## Self-test notes

(1) Variable independence: the only `sp.diff` call is `sp.diff(E, theta)` at line 140 where `E = (q^T K_eff q)/2` with `q = (cos(theta), sin(theta))`; E genuinely depends on theta. (2) Symmetry/parity: the proposed stationarity identity uses `sin(2 theta)` (odd in theta) and `cos(2 theta)` (even in theta) — re-derived above by hand and matched the script. (3) Trivial-case checks: setting `alpha -> 0` reduces `tan(2 theta) -> 0` and `theta_- -> 0` (the bare flat branch, consistent with Notes §4.1); setting `kappa1 -> 0` reduces `tan(2 theta) -> 0` (no off-diagonal coupling, no rotation). (4) For the discriminant identity `tr^2 - 4 det = disc`, re-derived by hand: matches. (5) For `alpha_crit` denominator algebra, re-derived `K_1 kappa0^2 + K_0 kappa1^2 = (8/(9 pi^2))(11 K_0 + 9 DeltaK_ax)`; matches Mathematica's closed form. No traps tripped.
