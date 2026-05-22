---
unit_id: 028
batch: II.1
created_at: 2026-05-21T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-21T23:06:10Z
findings_applied: 4
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 028

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.wl:47-48`

**Issue:** Two `expectZero` assertions check numeric identities (`kappa0^2 - kappa1^2 = 56/(9 Pi^2)` and `2 kappa0 kappa1 = -16 Sqrt[2]/(3 Pi^2)`) that are guaranteed by the literal definitions of `kappa0` and `kappa1` on lines 33-34. The assertions cannot fail unless Mathematica's `FullSimplify` malfunctions. The sympy script handles the same identities correctly with bare `Print` statements on lines 87-88.

**Required change:**
Replace lines 47-48 with `Print` statements that display the values without wrapping them in `expectZero`. Before:
```
expectZero["kappa0^2 - kappa1^2 - 56/(9 Pi^2)", kappa0^2 - kappa1^2 - 56/(9*Pi^2)];
expectZero["2 kappa0 kappa1 + 16 Sqrt[2]/(3 Pi^2)", 2*kappa0*kappa1 + 16*Sqrt[2]/(3*Pi^2)];
```
After:
```
Print["kappa0^2 - kappa1^2 = ", fmt[FullSimplify[kappa0^2 - kappa1^2, Assumptions -> $Assumptions]]];
Print["2 kappa0 kappa1 = ", fmt[FullSimplify[2*kappa0*kappa1, Assumptions -> $Assumptions]]];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 028` and confirm the new `Print` lines appear in the saved output, the two `PASS: kappa0^2 - kappa1^2 ...` and `PASS: 2 kappa0 kappa1 ...` lines no longer appear, and the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.wl`
- summary: Replaced the tautological kappa identity assertions with plain printed simplified values.
- deviation: none

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage028_loaded_profile_selection_sympy_audit.py:204-208`

**Issue:** Line 197 sets `alpha_crit = sp.solve(sp.Eq(det_eff, 0), alpha)[0]`. Substituting `alpha_crit` back into `det_eff` (lines 204-208) is then identically zero by sympy's `solve` contract. The check verifies `solve`'s return-value, not the physics. The substantive check that `det_eff` vanishes at the closed-form `alpha_crit_expected` (defined on line 198) is not currently performed.

**Required change:**
Edit lines 204-208 of `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage028_loaded_profile_selection_sympy_audit.py`. Before:
```python
    # Check determinant sign switch.
    eps = sp.symbols("eps", positive=True, real=True)
    det_below = sp.simplify(det_eff.subs(alpha, alpha_crit * (1 - eps)))
    det_at = sp.simplify(det_eff.subs(alpha, alpha_crit))
    print("det(alpha_crit) =", det_at)
    expect_zero("det(alpha_crit)", det_at)
```
After:
```python
    # Check determinant sign switch.
    eps = sp.symbols("eps", positive=True, real=True)
    det_below = sp.simplify(det_eff.subs(alpha, alpha_crit_expected * (1 - eps)))
    det_at_expected = sp.simplify(det_eff.subs(alpha, alpha_crit_expected))
    print("det(alpha_crit_expected) =", det_at_expected)
    expect_zero("det(alpha_crit_expected)", det_at_expected)
```
(Both the substitution into `det_below` and the new substantive check now use `alpha_crit_expected` — the quoted closed form — instead of the `sp.solve` result.)

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 028` and confirm the output line label is now `det(alpha_crit_expected) = 0`, the script still exits 0, and the linear-in-eps factored form printed at line 209-210 is unchanged.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage028_loaded_profile_selection_sympy_audit.py`
- summary: Switched the determinant-at-threshold check and below-threshold substitution to use the quoted closed-form `alpha_crit_expected`.
- deviation: none

## F3 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage028_loaded_profile_selection_sympy_audit.py:149-154`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.wl:86-89`

**Issue:** Both scripts check that `-rhs` equals `2 alpha (-kappa0 kappa1)/(DeltaK + alpha (kappa0^2-kappa1^2))`. The two expressions are identical by elementary sign rewriting (`-2 alpha k0 k1` vs `2 alpha (-k0 k1)`), so the check is algebraically guaranteed. The comment "manifestly positive form" suggests the author intended to verify that `-kappa0 kappa1 > 0`, but no positivity test is performed.

**Required change:**

(A) In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage028_loaded_profile_selection_sympy_audit.py`, replace lines 149-154. Before:
```python
    print("tan(2 theta) =", rhs)
    negative_tan = sp.simplify(-rhs)
    negative_tan_target = sp.simplify(
        2 * alpha * (-kappa0 * kappa1) / (DeltaK + alpha * (kappa0**2 - kappa1**2))
    )
    expect_zero("-tan(2 theta) - manifestly positive form", negative_tan - negative_tan_target)
    print("-tan(2 theta) =", negative_tan)
```
After:
```python
    print("tan(2 theta) =", rhs)
    negative_tan = sp.simplify(-rhs)
    print("-tan(2 theta) =", negative_tan)
    # Manifest positivity of -kappa0*kappa1: kappa0 > 0, kappa1 < 0.
    expect_zero("kappa0 sign check (kappa0 > 0)", sp.simplify(kappa0 - sp.Abs(kappa0)))
    expect_zero("kappa1 sign check (kappa1 < 0)", sp.simplify(kappa1 + sp.Abs(kappa1)))
```

(B) In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.wl`, replace lines 85-89. Before:
```
Print["tan(2 theta) = ", fmt[rhs]];
expectZero[
  "-tan(2 theta) - manifestly positive form",
  -rhs - 2*alpha*(-kappa0*kappa1)/(deltaK + alpha*(kappa0^2 - kappa1^2))
];
```
After:
```
Print["tan(2 theta) = ", fmt[rhs]];
Print["-tan(2 theta) = ", fmt[FullSimplify[-rhs, Assumptions -> $Assumptions]]];
(* Manifest positivity of -kappa0*kappa1: kappa0 > 0, kappa1 < 0. *)
expectZero["kappa0 sign check (kappa0 > 0)", kappa0 - Abs[kappa0]];
expectZero["kappa1 sign check (kappa1 < 0)", kappa1 + Abs[kappa1]];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 028` and `redteam exec-mathematica 028`, and confirm both outputs show `PASS: kappa0 sign check (kappa0 > 0)` and `PASS: kappa1 sign check (kappa1 < 0)` (or sympy equivalent printed as zero residuals), no longer carry the `-tan(2 theta) - manifestly positive form` lines, and both scripts exit 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage028_loaded_profile_selection_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.wl`
- summary: Replaced the tautological positive-form tan check with explicit kappa sign checks in both scripts.
- deviation: none

## F4 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.wl:66, 106-109`

**Issue:** The Mathematica script ports the sympy script's hand-written ansatze for the eigenvalues (`lambdaMinus = (trExpected - Sqrt[disc])/2`) and for `alphaCrit = K0*K1/(K1*kappa0^2 + K0*kappa1^2)`. Neither quantity is derived by an independent Mathematica computation. The fix adds two Mathematica-native cross-checks (`Eigenvalues[]` and `Solve[]`) so at least these two central results have an independent path. The existing assertions remain.

**Required change:**

(C1) After line 65 of `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.wl`, insert the following block before the existing line 66 (`charResidual = ...`):
```
eigvalsDirect = Eigenvalues[kEff];
expectZero[
  "Eigenvalues[kEff] sum vs trace",
  FullSimplify[Total[eigvalsDirect] - (lambdaMinus + lambdaPlus), Assumptions -> $Assumptions]
];
expectZero[
  "Eigenvalues[kEff] product vs determinant",
  FullSimplify[Times @@ eigvalsDirect - lambdaMinus*lambdaPlus, Assumptions -> $Assumptions]
];
```
(The sum and product of the unordered eigenvalue set are robust to Mathematica's ordering convention and reduce to known identities that any wrong closed-form for `lambdaMinus`/`lambdaPlus` would break.)

(C2) Immediately before line 106 (which currently reads `alphaCrit = FullSimplify[K0*K1/(K1*kappa0^2 + K0*kappa1^2), ...]`), insert:
```
alphaCritSolved = alpha /. First[Solve[detEff == 0, alpha]];
expectZero[
  "alphaCrit solved vs ratio closed form",
  FullSimplify[alphaCritSolved - K0*K1/(K1*kappa0^2 + K0*kappa1^2), Assumptions -> $Assumptions]
];
```
Do not modify the existing `alphaCrit`, `alphaCritClosed`, or downstream `detEff /. alpha -> alphaCrit` lines.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 028` and confirm three new `PASS:` lines appear in the saved output (`Eigenvalues[kEff] sum vs trace`, `Eigenvalues[kEff] product vs determinant`, `alphaCrit solved vs ratio closed form`), all previously-passing checks still pass, and the script exits 0.

## Applied: F4

- files_changed:
  - `mathematica/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.wl`
- summary: Added Mathematica-native eigenvalue and solved-threshold cross-checks for the closed-form ledger quantities.
- deviation: Stripped Mathematica's `ConditionalExpression` wrapper from the `Solve` result so the inserted `expectZero` check reduces to literal zero.
