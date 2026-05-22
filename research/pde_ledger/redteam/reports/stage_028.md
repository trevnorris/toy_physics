---
unit_id: 028
batch: II.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-21T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 028 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage028_loaded_profile_selection_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage028_loaded_profile_selection_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.txt`

Output mtimes (both newer than their script mtimes): sympy output 2026-05-11 12:40 vs script 2026-04-03 12:05; Mathematica output 2026-05-11 12:44 vs script 2026-05-11 11:56. Outputs are fresh.

## What the script claims to verify

In the {u0,u1} N/N basis, the scripts construct the loaded wall operator `K_eff = K_bare - alpha v v^T` with `v = (kappa0, kappa1) = (2 sqrt(2)/pi, -4/(3 pi))`, `K_bare = diag(K0, K1)`, `K1 = K0 + DeltaK`, and `DeltaK = pi^2 T_w / L^2`. They then assert (i) closed forms for `trace(K_eff)` and `det(K_eff)`, (ii) a characteristic-polynomial factorization with explicit eigenvalues `lambda_+/-` using a discriminant `disc = (DeltaK + alpha(kappa0^2-kappa1^2))^2 + 4 alpha^2 kappa0^2 kappa1^2`, (iii) the stationarity condition `dE/dtheta = stationarity_expected/2` and the closed-form `tan(2 theta) = 2 alpha kappa0 kappa1 / (DeltaK + alpha(kappa0^2-kappa1^2))`, (iv) the weak-loading leading coefficient `kappa0 kappa1/DeltaK`, (v) the strong-loading limit `tan(2 theta) -> tan(2 theta_max)` with `tan(theta_max) = kappa1/kappa0 = -sqrt(2)/3`, and (vi) the softening threshold `alpha_crit = K0 K1 / (K1 kappa0^2 + K0 kappa1^2)` with `det(K_eff)|_{alpha_crit} = 0` and a linear-in-eps determinant just below criticality.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | mathematica | 47 | `expectZero["kappa0^2 - kappa1^2 - 56/(9 Pi^2)", ...]` | no (tautology by literal substitution) |
| A2 | mathematica | 48 | `expectZero["2 kappa0 kappa1 + 16 Sqrt[2]/(3 Pi^2)", ...]` | no (tautology by literal substitution) |
| A3 | sympy | 109 | `expect_zero("trace - expected", tr_eff - tr_expected)` | yes |
| A4 | sympy | 110 | `expect_zero("det - expected", det_eff - det_expected)` | yes |
| A5 | mathematica | 57 | `expectZero["trace - expected", ...]` | yes |
| A6 | mathematica | 58 | `expectZero["det - expected", ...]` | yes |
| A7 | sympy | 119 | `expect_zero("characteristic-factorization check", char)` | yes |
| A8 | mathematica | 67 | `expectZero["characteristic factorization", ...]` | yes |
| A9 | sympy | 146 | `expect_zero("dE/dtheta - stationarity_expected/2", ...)` | yes |
| A10 | mathematica | 84 | `expectZero["dE/dtheta - stationarity_expected/2", ...]` | yes |
| A11 | sympy | 153 | `expect_zero("-tan(2 theta) - manifestly positive form", ...)` | no (tautology: `-rhs` vs `2 alpha (-k0 k1)/(...)` — identical by sign rewrite) |
| A12 | mathematica | 86-89 | `expectZero["-tan(2 theta) - manifestly positive form", ...]` | no (same tautology) |
| A13 | sympy | 159 | `expect_zero("weak-loading coefficient - k0*k1/DeltaK", ...)` | yes |
| A14 | mathematica | 93 | `expectZero["weak-loading coefficient - kappa0 kappa1/deltaK", ...]` | yes |
| A15 | sympy | 180 | `expect_zero("strong-loading limit - tan(2 theta_max)", ...)` | yes |
| A16 | sympy | 181 | `expect_zero("tan(theta_max) + sqrt(2)/3", ...)` | yes |
| A17 | mathematica | 103 | `expectZero["strong-loading limit - tan(2 theta_max)", ...]` | yes |
| A18 | mathematica | 104 | `expectZero["tan(theta_max) + Sqrt[2]/3", ...]` | yes |
| A19 | sympy | 201 | `expect_zero("alpha_crit - expected", ...)` | yes (verifies `sp.solve` returns the closed form) |
| A20 | mathematica | 109 | `expectZero["alpha_crit - finite-throat closed form", ...]` | yes (verifies two algebraic forms agree) |
| A21 | sympy | 208 | `expect_zero("det(alpha_crit)", det_at)` | no (`alpha_crit` was produced by `sp.solve(det_eff==0, alpha)`, so substitution back is identically zero) |
| A22 | mathematica | 110 | `expectZero["det(alpha_crit)", detEff /. alpha -> alphaCrit]` | yes (`alphaCrit` is the quoted closed form, not solved; substitution can fail if the closed form is wrong) |

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.wl:47-48`

**What's wrong:**
The two Mathematica assertions
```
expectZero["kappa0^2 - kappa1^2 - 56/(9 Pi^2)", kappa0^2 - kappa1^2 - 56/(9*Pi^2)];
expectZero["2 kappa0 kappa1 + 16 Sqrt[2]/(3 Pi^2)", 2*kappa0*kappa1 + 16*Sqrt[2]/(3*Pi^2)];
```
are algebraically guaranteed by the literal definitions on lines 33-34 (`kappa0 = 2 Sqrt[2]/Pi`, `kappa1 = -4/(3 Pi)`). Direct substitution gives `(2 Sqrt[2]/Pi)^2 - (-4/(3 Pi))^2 = 8/Pi^2 - 16/(9 Pi^2) = 56/(9 Pi^2)` and `2 * (2 Sqrt[2]/Pi)*(-4/(3 Pi)) = -16 Sqrt[2]/(3 Pi^2)`. The checks cannot fail unless `FullSimplify` itself malfunctions; they verify arithmetic, not physics. The sympy script (lines 87-88) only `print`s the same identities, which is the correct treatment for a sanity restatement.

**Why this matters:**
Tautological assertions inflate the apparent verification coverage without exercising the unit's claims. They also waste a slot in the Mathematica audit that could carry an independent check.

**Required change:**
Either (a) demote both lines to `Print[...]` (matching the sympy treatment), or (b) replace them with cross-engine consistency anchors that are not trivially true. The minimal fix is (a): remove the two `expectZero` wrappers and leave bare `Print` statements showing the values.

**Verification:**
After the fix, the Mathematica output should still show the numeric identities printed but no longer carry `PASS:` lines for them. The script must still exit 0.

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage028_loaded_profile_selection_sympy_audit.py:208`

**What's wrong:**
Line 197 sets `alpha_crit = sp.solve(sp.Eq(det_eff, 0), alpha)[0]`, i.e. by definition `det_eff` evaluates to zero at `alpha = alpha_crit`. Line 207-208 then substitutes `alpha_crit` back into `det_eff` and asserts the result is zero:
```python
det_at = sp.simplify(det_eff.subs(alpha, alpha_crit))
print("det(alpha_crit) =", det_at)
expect_zero("det(alpha_crit)", det_at)
```
This is algebraically guaranteed by `sp.solve`'s contract; it can only fail if sympy's simplification cannot collapse the resubstitution, in which case the failure is a simplification artefact rather than a physics violation. The Mathematica counterpart at line 110 is NOT tautological for the same reason because `alphaCrit` there is the quoted closed form `K0*K1/(K1*kappa0^2 + K0*kappa1^2)`, not a solve result.

**Why this matters:**
The substantive softening-threshold check is `expect_zero("alpha_crit - expected", alpha_crit - alpha_crit_expected)` on line 201, which can fail if `sp.solve` returns a different root form. The line-208 check adds nothing beyond confirming `sp.solve` did not lie about its return value.

**Required change:**
Replace the resubstitution check on line 208 with a substitution of the closed-form `alpha_crit_expected` (defined on line 198) into `det_eff`:
```python
det_at = sp.simplify(det_eff.subs(alpha, alpha_crit_expected))
print("det(alpha_crit_expected) =", det_at)
expect_zero("det(alpha_crit_expected)", det_at)
```
This mirrors the Mathematica path (substitute closed form, confirm zero) and removes the tautology.

**Verification:**
After the fix, the saved output should show the line label changed to `det(alpha_crit_expected)` and the value still zero. Script must exit 0.

### F3 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage028_loaded_profile_selection_sympy_audit.py:149-153`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.wl:86-89`

**What's wrong:**
Both scripts define
```
rhs = 2 alpha kappa0 kappa1 / (DeltaK + alpha (kappa0^2 - kappa1^2))
```
and then check that `-rhs` equals the same expression with the numerator written as `2 alpha (-kappa0 kappa1)` instead of `-(2 alpha kappa0 kappa1)`. Sympy line 150-152:
```python
negative_tan = sp.simplify(-rhs)
negative_tan_target = sp.simplify(
    2 * alpha * (-kappa0 * kappa1) / (DeltaK + alpha * (kappa0**2 - kappa1**2))
)
expect_zero("-tan(2 theta) - manifestly positive form", negative_tan - negative_tan_target)
```
Mathematica lines 86-89:
```
expectZero["-tan(2 theta) - manifestly positive form",
  -rhs - 2*alpha*(-kappa0*kappa1)/(deltaK + alpha*(kappa0^2 - kappa1^2))];
```
`(-1)*(kappa0*kappa1) === -(kappa0*kappa1)` algebraically — these two expressions are identical before any simplifier sees them. The check is guaranteed by construction. The comment "manifestly positive form" suggests the author intended to verify that `-kappa0 kappa1 > 0`, but the assertion does not actually test positivity.

**Why this matters:**
The "manifestly positive" framing is misleading. The check does not verify any sign or positivity property; it only restates `-rhs` with a sign moved from outside the fraction to inside the numerator factor.

**Required change:**
Replace the tautological assertion with one that actually tests the intended property — namely that `-rhs > 0` for `alpha > 0` (equivalently `-kappa0 kappa1 > 0`).

In sympy, line 150-153, replace with:
```python
print("-tan(2 theta) =", sp.simplify(-rhs))
# Manifest positivity: -kappa0*kappa1 > 0 since kappa0 > 0, kappa1 < 0.
expect_zero("kappa0 sign check (kappa0 > 0)", sp.simplify(kappa0 - sp.Abs(kappa0)))
expect_zero("kappa1 sign check (kappa1 < 0)", sp.simplify(kappa1 + sp.Abs(kappa1)))
```
In Mathematica, lines 86-89, replace with:
```
Print["-tan(2 theta) = ", fmt[FullSimplify[-rhs, Assumptions -> $Assumptions]]];
expectZero["kappa0 sign check (kappa0 > 0)", kappa0 - Abs[kappa0]];
expectZero["kappa1 sign check (kappa1 < 0)", kappa1 + Abs[kappa1]];
```
Both new checks can fail if `kappa0` or `kappa1` had been declared with the wrong sign (i.e., they test the substantive claim implied by the original comment).

**Verification:**
After the fix, the saved outputs should show both sign checks passing and no longer carry the tautological `-tan(2 theta) - manifestly positive form` line. Scripts must exit 0.

### F4 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.wl:33-110` (whole script)
- compare to `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage028_loaded_profile_selection_sympy_audit.py:60-201`

**What's wrong:**
The Mathematica script ports the sympy choreography line-for-line. Three matching segments:

(i) Symbol setup (sympy lines 55-64 vs Mathematica lines 28-37):
```python
kappa0 = sp.simplify(2 * sp.sqrt(2) / sp.pi)
kappa1 = sp.simplify(-4 / (3 * sp.pi))
DeltaK = sp.simplify(T_w * sp.pi**2 / L**2)
K0 = sp.simplify(K_eta + 6 * T_Omega); K1 = sp.simplify(K0 + DeltaK)
```
```mathematica
kappa0 = FullSimplify[2*Sqrt[2]/Pi, ...]; kappa1 = FullSimplify[-4/(3*Pi), ...];
deltaK = FullSimplify[Tw*Pi^2/L^2, ...];
K0 = Keta + 6*TOmega; K1 = K0 + deltaK;
```

(ii) Eigenvalue closed form (sympy lines 112-114 vs Mathematica lines 60-65):
```python
disc = sp.simplify((DeltaK + alpha*(kappa0**2 - kappa1**2))**2 + 4*alpha**2*kappa0**2*kappa1**2)
lam_minus = sp.simplify((tr_expected - sp.sqrt(disc)) / 2)
lam_plus  = sp.simplify((tr_expected + sp.sqrt(disc)) / 2)
```
```mathematica
disc = FullSimplify[(deltaK + alpha*(kappa0^2 - kappa1^2))^2 + 4*alpha^2*kappa0^2*kappa1^2, ...];
lambdaMinus = FullSimplify[(trExpected - Sqrt[disc])/2, ...];
lambdaPlus  = FullSimplify[(trExpected + Sqrt[disc])/2, ...];
```

(iii) Stationarity closed form (sympy line 144 vs Mathematica lines 79-82):
```python
stationarity_expected = sp.simplify((DeltaK + alpha*(kappa0**2 - kappa1**2)) * sp.sin(2*theta) - 2*alpha*kappa0*kappa1*sp.cos(2*theta))
```
```mathematica
stationarityExpected = FullSimplify[(deltaK + alpha*(kappa0^2 - kappa1^2))*Sin[2*theta] - 2*alpha*kappa0*kappa1*Cos[2*theta], ...];
```

Every closed-form ansatz the Mathematica script tests against was already hand-written in the same algebraic form in the sympy script. Mathematica does NOT call `Eigenvalues[kEff]` to derive `lambda_+/-` independently, and it does NOT call `Solve[detEff == 0, alpha]` to derive `alphaCrit` independently — both quantities are hand-coded as the same closed forms sympy uses. The only genuinely independent path is the implicit one inside each engine's `FullSimplify` / `sp.simplify`. The second-engine policy expects at least one substantive derivation by a path not used in the other engine.

**Why this matters:**
Both engines simplifying the same hand-written algebra is weaker than two independent derivations. If the hand-written closed forms (e.g. `disc`, `stationarityExpected`, `alpha_crit_expected`) are wrong in the same way in both engines, neither check catches the error.

**Required change:**
Add at least two genuinely independent derivations to the Mathematica script (do not modify the sympy script):

1. Compute the eigenvalues directly via Mathematica's `Eigenvalues[]` and assert they match `{lambdaPlus, lambdaMinus}` (up to ordering). Insert after line 65, before line 66:
```
eigvalsDirect = Sort[Eigenvalues[kEff], Less] /. {a_, b_} :> {a, b};
expectZero[
  "Eigenvalues[kEff] vs lambdaMinus",
  FullSimplify[First[Sort[{lambdaMinus, lambdaPlus}, Less]] - First[Sort[Eigenvalues[kEff], Less]], Assumptions -> $Assumptions]
];
expectZero[
  "Eigenvalues[kEff] vs lambdaPlus",
  FullSimplify[Last[Sort[{lambdaMinus, lambdaPlus}, Less]] - Last[Sort[Eigenvalues[kEff], Less]], Assumptions -> $Assumptions]
];
```
(If the `Sort[..., Less]` over symbolic eigenvalues does not reduce, fall back to comparing the unordered set: `expectZero["Eigenvalues[kEff] product vs lambdaMinus*lambdaPlus", Times @@ Eigenvalues[kEff] - lambdaMinus*lambdaPlus]` and `expectZero["Eigenvalues[kEff] sum vs lambdaMinus+lambdaPlus", Total[Eigenvalues[kEff]] - lambdaMinus - lambdaPlus]`. The sum/product variant is robust to ordering.)

2. Derive `alphaCrit` via `Solve[]` instead of quoting the closed form. Insert before line 106:
```
alphaCritSolved = alpha /. First[Solve[detEff == 0, alpha]];
expectZero["alphaCrit solved vs closed form", FullSimplify[alphaCritSolved - K0*K1/(K1*kappa0^2 + K0*kappa1^2), Assumptions -> $Assumptions]];
```
Then leave the existing `alphaCrit` line and its downstream uses unchanged.

**Verification:**
After the fix, the Mathematica output must show new `PASS:` lines for `Eigenvalues[kEff] sum vs lambdaMinus+lambdaPlus`, `Eigenvalues[kEff] product vs lambdaMinus*lambdaPlus` (or the ordered variants), and `alphaCrit solved vs closed form`. Script must exit 0.

## Independent-derivation check (Mathematica)

The Mathematica script is a near-line-by-line port of the sympy script for the symbol setup, eigenvalue closed forms, and stationarity ansatz (see F4 quotations above). The only divergent paths are: (a) Mathematica adds a closed-form `alphaCritClosed = 9 Pi^2 K0 K1 / (8 (11 K0 + 9 deltaK))` that sympy does not write, and confirms `alphaCrit - alphaCritClosed == 0`; (b) sympy uses `sp.solve(det_eff == 0, alpha)` whereas Mathematica quotes `alphaCrit` directly. Neither divergence is a substantive independent derivation of the central claim. Both engines re-simplify the same hand-written ansatze rather than re-deriving them.

## Engine cross-check

Both engines agree on every assertion they share (trace, determinant, characteristic factorization, stationarity, tan(2 theta), weak-loading coefficient, strong-loading limit, tan(theta_max), alpha_crit closed form, det(alpha_crit) vanishing). The outputs are visibly the same under different surface forms:

| Quantity | Sympy | Mathematica |
|---|---|---|
| tan(2 theta) | `-48 sqrt(2) L^2 alpha / (56 L^2 alpha + 9 pi^4 T_w)` | `(-48 Sqrt[2] alpha L^2)/(56 alpha L^2 + 9 Pi^4 Tw)` |
| weak coeff | `-8 sqrt(2) L^2 / (3 pi^4 T_w)` (in `weak/alpha` form) | `(-8 Sqrt[2] L^2)/(3 Pi^4 Tw)` |
| lim_{alpha→∞} tan(2 theta) | `-6 sqrt(2)/7` | `(-6 Sqrt[2])/7` |
| tan(theta_max) | `-sqrt(2)/3` | `-1/3 Sqrt[2]` |
| alpha_crit | `9 pi^2 (K_eta^2 L^2 + 12 K_eta L^2 T_Omega + pi^2 K_eta T_w + 36 L^2 T_Omega^2 + 6 pi^2 T_Omega T_w) / (8 (11 K_eta L^2 + 66 L^2 T_Omega + 9 pi^2 T_w))` | `(9 Pi^2 (Keta + 6 TOmega)(L^2 (Keta + 6 TOmega) + Pi^2 Tw))/(88 L^2 (Keta + 6 TOmega) + 72 Pi^2 Tw)` |
| det(alpha_crit (1-eps)) | `eps (K_eta + 6 T_Omega)(K_eta L^2 + 6 L^2 T_Omega + pi^2 T_w)/L^2` | `eps (Keta + 6 TOmega)(Keta + 6 TOmega + Pi^2 Tw/L^2)` |

The `alpha_crit` forms differ only by whether the numerator and denominator are expanded; `K_eta^2 + 12 K_eta T_Omega + 36 T_Omega^2 = (K_eta + 6 T_Omega)^2` matches directly, and the denominator factorizations also agree. No engine disagreement.

## Verdict justification

`findings: 4`. The substantive physics checks (trace/det closed forms, characteristic-polynomial factorization, stationarity ↔ closed-form tan(2 theta), weak-loading coefficient, strong-loading limit ↔ tan(2 theta_max), softening threshold) all hold under attack: I rederived the discriminant identity `tr^2 - 4 det = disc`, the closed-form `tan(2 theta)`, the strong-loading limit `-6 sqrt(2)/7`, and `alpha_crit = 9 pi^2 K0 K1 / (8 (11 K0 + 9 DeltaK))` by hand and all match the script's outputs.

The four findings are:
- F1, F2, F3 are low-severity tautologies that inflate the assertion count without exercising the claim.
- F4 is a medium-severity transliteration: the Mathematica script ports the sympy ansatze for eigenvalues, stationarity, and `alpha_crit` rather than deriving any of them by an independent Mathematica path. The fix is additive (add `Eigenvalues[kEff]` and `Solve[detEff==0]` cross-checks) and does not rewrite the existing structure.

No `stop_cold` warranted: the underlying math is consistent, and the fixes are surgical script-only edits.

## Self-test notes

I checked: (1) Variable-independence for the proposed `Solve[detEff == 0, alpha]` — `detEff` is built from `Det[kEff]` where `kEff` depends on `alpha`, so the solve is well-posed; (2) the proposed `Eigenvalues[kEff]` sum/product cross-check is robust to ordering and reduces to identities I confirmed by hand (sum = `trace(kEff) = trExpected`, product = `det(kEff) = detExpected`); (3) the sympy F2 fix substitutes `alpha_crit_expected` (already defined on line 198) into `det_eff`, and since `det_eff = K0 K1 - alpha (K1 kappa0^2 + K0 kappa1^2)` is linear in alpha, substituting `K0 K1/(K1 kappa0^2 + K0 kappa1^2)` gives `K0 K1 - K0 K1 = 0` symbolically — the check will pass after sympy's simplifier handles the rational; (4) the F3 sign checks `kappa0 - Abs[kappa0]` and `kappa1 + Abs[kappa1]` resolve to 0 only when `kappa0 >= 0` and `kappa1 <= 0` respectively — since `kappa0 = 2 Sqrt[2]/Pi > 0` and `kappa1 = -4/(3 Pi) < 0` are concrete numeric literals, both expressions simplify to 0 in both engines.
