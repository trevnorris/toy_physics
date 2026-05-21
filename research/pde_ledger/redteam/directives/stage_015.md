---
unit_id: 015
batch: I.2
created_at: 2026-05-21T18:30:15Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-21T19:12:39Z
findings_applied: 3
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 015

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — missing_verification_script

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl` (new file)

**Issue:** Unit 015 is not status-only and not a checkpoint carve-out. Only the SymPy script exists; there is no Mathematica companion. A second-engine script is required so the engine cross-check is non-trivially satisfied. The new `.wl` must independently re-derive the unit's claims, not transliterate the `.py`.

**Required change:**
Create the new `.wl` file with the structure below. The script must terminate with `Print["STATUS: PASS"]` only after all claim checks pass; each individual check must call `Exit[1]` on failure (e.g. `If[FullSimplify[residual] =!= 0, Print["FAIL: ", label]; Exit[1]]`). Use Mathematica primitives `Series`, `Coefficient`, `D`, `Solve`, `Det`, `Simplify`, `FullSimplify`, `Integrate`, `ThreeJSymbol` (not `Gaunt`) rather than mirroring the SymPy choreography. Use distinct local variable names to avoid line-by-line correspondence — e.g., use `lagrangian`, `effectiveMass`, `gateMatrix`, `realY20Ratio[m_]` rather than `L`, `K_eta`, `wall_matrix`, `real_y20_square_ratio`.

In particular:
- For the K_eta derivation, expand the Lagrangian via `Series[lagrangian, {eps, 0, 2}] // Normal` and read the `eps^2` coefficient via `Coefficient[..., eps, 2]`, instead of `D[L, eps, 2]/2` style.
- For the wall-only gate Jacobian, build `K1full = D21full + D01full/9` and `Heven = D41full - (2/3) D21full - D01full/27` directly from the symbolic definitions of `D01, D21, D41` (do NOT skip the construction and just type the wall-only result; verify the specialization explicitly via `K1wall = K1full /. {b01 -> 0, b21 -> 0, b41 -> 0, z01 -> 0, z21 -> 0, z41 -> 0}`).
- For the Gaunt ratios, compute `gauntBase = ThreeJSymbol[{2, 0}, {2, 0}, {2, 0}]^2 * (some normalization)` or use the standard Wigner 3-j form directly; do NOT short-circuit the m=0 lane.

Claim manifest (each claim must have at least one explicit `If[..., Exit[1]]` check):

M1. **Generic IBP product-rule identity.** Define `aFun[w_] := A[w]`, `etaFun[w_] := eta[w]`. Verify symbolically that `-aFun[w] etaFun[w] etaFun'[w] - (D[-aFun[w] etaFun[w]^2/2, w] + aFun'[w] etaFun[w]^2/2) === 0` after `FullSimplify`. Also verify the sign-mutation `D[-(-aFun[w] etaFun[w]^2/2), w]` produces a nonzero residual.

M2. **Concrete Gaussian IBP boundary discharge.** With `aConcrete = Exp[-w^2]`, `etaConcrete = Exp[-w^2/2]`, evaluate `Limit[-aConcrete etaConcrete^2/2, w -> Infinity] - Limit[-aConcrete etaConcrete^2/2, w -> -Infinity]` and verify it is zero. Evaluate `Integrate[-aConcrete etaConcrete D[etaConcrete, w], {w, -Infinity, Infinity}] - Integrate[D[aConcrete, w] etaConcrete^2/2, {w, -Infinity, Infinity}]` and verify it is zero (the boundary term vanishes for Gaussians, so cross integral = bulk integral).

M3. **K_eta closed form from L expansion.** Set up symbols `mu0, Tw0, TO0, U0, TwR0, TwRR0, UR0, URR0, R0p, dTwRR0p` (use `dTwRR0p` for the script's `d_TwR_R0p`) and `eta, etat, etaw, grad2, eps`. Define:
```
Tw = Tw0 + eps TwR0 eta + eps^2 TwRR0 eta^2 / 2;
U  = U0  + eps UR0 eta  + eps^2 URR0 eta^2 / 2;
Rt = eps etat;
Rw = R0p + eps etaw;
lagrangian = mu0 Rt^2/2 - Tw Rw^2/2 - TO0 eps^2 grad2/2 - U;
L2raw = Coefficient[Series[lagrangian, {eps, 0, 2}] // Normal, eps, 2];
crossCoeff = D[D[L2raw, eta], etaw];
```
Verify `FullSimplify[crossCoeff + TwR0 R0p] === 0`. Verify that
```
canonicalL2 = mu0 etat^2/2 - Tw0 etaw^2/2 - TO0 grad2/2 - (URR0 - dTwRR0p + TwRR0 R0p^2/2) eta^2/2;
L2afterIBP  = Expand[L2raw - (-TwR0 R0p eta etaw) + dTwRR0p eta^2/2];
```
satisfies `FullSimplify[L2afterIBP - canonicalL2] === 0`. Mutate the sign on `dTwRR0p` (i.e., use `URR0 + dTwRR0p + TwRR0 R0p^2/2`) and verify the residual is nonzero.

M4. **Wall-only K1 and H_even specialization.** Define symbols `dK, dM, b01, b21, b41, z01, z21, z41`. Build:
```
D01full = dK - b01 - z01;
D21full = -(dM + b21 + z21);
D41full = -(b41 + z41);
K1full   = D21full + D01full/9;
Hevenfull = D41full - (2/3) D21full - D01full/27;
wallSpec  = {b01 -> 0, b21 -> 0, b41 -> 0, z01 -> 0, z21 -> 0, z41 -> 0};
K1wall    = Expand[K1full   /. wallSpec];
Hevenwall = Expand[Hevenfull /. wallSpec];
```
Verify `FullSimplify[K1wall - (-dM + dK/9)] === 0` and `FullSimplify[Hevenwall - (2/3 dM - dK/27)] === 0`.

M5. **Concrete overlap-integral wall-only check (independent of the SymPy script's symbolic-substitution tautology).** Use concrete Gaussian profiles to evaluate the overlap integrals defined symbolically in lines 110-114 of the SymPy script:
```
betaConcrete    = Exp[-w^2];
deltaMu         = Exp[-w^2];
deltaTw         = Exp[-w^2];
deltaTO         = Exp[-w^2];
deltaKeta       = Exp[-w^2];
dMoverlap = Integrate[deltaMu betaConcrete^2, {w, -Infinity, Infinity}];
dKoverlap = Integrate[deltaTw D[betaConcrete, w]^2 + (deltaKeta + 6 deltaTO) betaConcrete^2, {w, -Infinity, Infinity}];
```
These should evaluate to closed-form `Sqrt[Pi]`-multiples. Then verify:
```
K1wallNum   = K1wall   /. {dK -> dKoverlap, dM -> dMoverlap};
HevenwallNum = Hevenwall /. {dK -> dKoverlap, dM -> dMoverlap};
```
satisfy `FullSimplify[K1wallNum - (-dMoverlap + dKoverlap/9)] === 0` and the analogous H_even check. Additionally, mutate the `6 deltaTO` coefficient to `5 deltaTO` in `dKoverlap` and verify the resulting `K1wallNum` no longer equals `-dMoverlap + dKoverlap/9` of the unmutated form (this is the substantive coefficient guard).

M6. **Wall-only Jacobian determinant.** Build
```
wallMatrix = {{D[K1wall, dK], D[K1wall, dM]}, {D[Hevenwall, dK], D[Hevenwall, dM]}};
```
Verify `FullSimplify[Det[wallMatrix] - 1/27] === 0`. Then build
```
K1wallParam      = -dM + gateCoeff dK;
wallMatrixParam  = {{D[K1wallParam, dK], D[K1wallParam, dM]}, {D[Hevenwall, dK], D[Hevenwall, dM]}};
wallDetShift     = FullSimplify[
                     (Det[wallMatrixParam] /. gateCoeff -> 1/9 + eps)
                     - (Det[wallMatrixParam] /. gateCoeff -> 1/9)];
```
Verify `wallDetShift =!= 0` (it should equal `2 eps/3`).

M7. **Wall-only solve.** Use `Solve[{K1wall == 0, Hevenwall == 0}, {dK, dM}]` and verify the unique solution is `{dK -> 0, dM -> 0}`. Then `Solve[{K1wall + eps == 0, Hevenwall == 0}, {dK, dM}]` must give nonzero `dK` and nonzero `dM`.

M8. **Real-Y20 squared overlap ratios.** Define
```
gauntBase = ThreeJSymbol[{2, 0}, {2, 0}, {2, 0}];
realY20Ratio[m_] := (-1)^m * ThreeJSymbol[{2, 0}, {2, m}, {2, -m}] / gauntBase;
```
(Equivalent up to the standard `(2l+1)/(4 Pi)` prefactor cancellation in the ratio; the ratio of Gaunt coefficients to the m=m=0 baseline equals the ratio of Wigner 3-j coefficients squared divided by the base squared, but for the m=0 case used here the ratio simplifies to the 3-j ratio directly. If using `Gaunt`-equivalent formulas, the prefactor cancels in the ratio.) Verify `FullSimplify[realY20Ratio[0] - 1] === 0`, `FullSimplify[realY20Ratio[1] - 1/2] === 0`, `FullSimplify[realY20Ratio[2] + 1] === 0`. Also verify the same-sign cross terms vanish: `ThreeJSymbol[{2, 0}, {2, m}, {2, m}] === 0` for `m == 1, 2` (forced by the selection rule `m1 + m2 + m3 == 0`).

M9. **Grouped trace and `b = 3a` identity.** With `lam20 = 1, lam21 = 1/2, lam22 = -1` and free symbols `x0, eps1`:
```
x20 = x0 + eps1 lam20;
x21 = x0 + eps1 lam21;
x22 = x0 + eps1 lam22;
xbar = (x20 + 2 x21 + 2 x22)/5;
ax   = (2 x20 - x21 - x22)/10;
bx   = (x21 - x22)/2;
```
Verify `FullSimplify[xbar - x0] === 0` and `FullSimplify[bx - 3 ax] === 0`.

**Claim manifest** (for missing-script findings only):

| # | Claim (Mathematica form) |
|---|--------------------------|
| M1 | `FullSimplify[-A[w] eta[w] eta'[w] - (D[-A[w] eta[w]^2/2, w] + A'[w] eta[w]^2/2)] === 0` and mutated sign nonzero |
| M2 | Boundary value of `-Exp[-2 w^2]/2` at +/-Infinity is zero; concrete cross integral equals concrete bulk integral |
| M3 | `Coefficient[Series[L, {eps,0,2}]//Normal, eps, 2]` cross coefficient `+TwR0 R0p == 0`; `L2afterIBP - canonicalL2 == 0` with `K_eta = URR0 - dTwRR0p + TwRR0 R0p^2/2`; sign-mutation nonzero |
| M4 | `K1wall = -dM + dK/9`, `Hevenwall = (2/3) dM - dK/27` from the explicit `D01, D21, D41` construction |
| M5 | With Gaussian overlap profiles, `dKoverlap` and `dMoverlap` evaluate to closed forms; wall-only K1 and H_even on those concrete values satisfy the gate algebra; the `6 deltaTO` coefficient mutation breaks it |
| M6 | `Det[wallMatrix] == 1/27`; perturbing the `1/9` coefficient to `1/9 + eps` shifts the determinant by `2 eps/3` |
| M7 | Solving `K1wall == 0, Hevenwall == 0` gives `{dK -> 0, dM -> 0}`; perturbed K1+eps system gives nonzero solutions |
| M8 | `realY20Ratio[0] = 1`, `realY20Ratio[1] = 1/2`, `realY20Ratio[2] = -1` from Wigner 3-j; same-sign cross terms vanish by selection rule |
| M9 | Grouped trace `xbar - x0 == 0`; `b - 3 a == 0` for the (a, b) duo built from (lam20, lam21, lam22) |

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 015` and confirms the new `.wl` exists under `mathematica/` and exits 0, and that each of M1-M9 produces at least one PASS-style print line.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl`
- summary: Created the Mathematica companion audit script with independent M1-M9 checks and PASS/FAIL exits.
- deviation: none

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py:129-136`

**Issue:** Lines 127-128 already establish `K1_wall == -dM + dK/9` and `H_even_wall == (2/3) dM - dK/27`. Lines 129-136 then assert that substituting `{dK: dK_overlap, dM: dM_overlap}` into those expressions equals the literal substitution into the right-hand side — by sympy's substitution mechanics alone, with no integral ever evaluated. The labels `wall-only K1 from overlap-generated slots` and `wall-only H_even from overlap-generated slots` imply a coupling between the overlap integrals defined at lines 110-114 and the gate algebra, but the integrals are never concretely evaluated.

**Required change:**

Replace lines 129-136 entirely with a concrete-Gaussian overlap evaluation that produces a numerical (closed-form `sqrt(pi)`-multiple) test of the gate algebra, plus a mutation guard on the `6 delta_TO` coefficient. Concretely:

DELETE lines 129-136 (the two `assert_zero` blocks for `wall-only K1 from overlap-generated slots` and `wall-only H_even from overlap-generated slots`).

INSERT in their place:
```
    # Concretize the symbolic overlaps with Gaussian profiles so the wall-only
    # K1/H_even algebra is exercised against actual closed-form integrals, not
    # against a substitution rename of the same algebra.
    beta_concrete = sp.exp(-wall_w**2)
    delta_mu_concrete = sp.exp(-wall_w**2)
    delta_Tw_concrete = sp.exp(-wall_w**2)
    delta_TO_concrete = sp.exp(-wall_w**2)
    delta_Keta_concrete = sp.exp(-wall_w**2)
    dM_overlap_concrete = sp.integrate(
        delta_mu_concrete * beta_concrete**2,
        (wall_w, -sp.oo, sp.oo),
    )
    dK_overlap_concrete = sp.integrate(
        delta_Tw_concrete * sp.diff(beta_concrete, wall_w)**2
        + (delta_Keta_concrete + 6 * delta_TO_concrete) * beta_concrete**2,
        (wall_w, -sp.oo, sp.oo),
    )
    assert_zero(
        "wall-only K1 from concrete Gaussian overlap integrals",
        K1_wall.subs({dK: dK_overlap_concrete, dM: dM_overlap_concrete})
        - (-dM_overlap_concrete + dK_overlap_concrete / 9),
    )
    assert_zero(
        "wall-only H_even from concrete Gaussian overlap integrals",
        H_even_wall.subs({dK: dK_overlap_concrete, dM: dM_overlap_concrete})
        - (sp.Rational(2, 3) * dM_overlap_concrete - dK_overlap_concrete / 27),
    )
    # Coefficient guard: changing the 6*delta_TO coefficient to 5*delta_TO must
    # break the K1 identity against the unmutated dK_overlap.
    dK_overlap_mutated = sp.integrate(
        delta_Tw_concrete * sp.diff(beta_concrete, wall_w)**2
        + (delta_Keta_concrete + 5 * delta_TO_concrete) * beta_concrete**2,
        (wall_w, -sp.oo, sp.oo),
    )
    assert_nonzero(
        "wall-only K1 detects mutated 6*delta_TO coefficient",
        K1_wall.subs({dK: dK_overlap_mutated, dM: dM_overlap_concrete})
        - (-dM_overlap_concrete + dK_overlap_concrete / 9),
    )
```

**Verification:**
After patching, the output transcript prints `wall-only K1 from concrete Gaussian overlap integrals`, `wall-only H_even from concrete Gaussian overlap integrals`, and `wall-only K1 detects mutated 6*delta_TO coefficient`, all with PASS. The labels `wall-only K1 from overlap-generated slots` and `wall-only H_even from overlap-generated slots` no longer appear. Re-running the script still exits 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py`
- summary: Replaced the overlap-slot substitution tautology with concrete Gaussian overlap integrals and a 6-to-5 coefficient mutation guard.
- deviation: none

## F3 — hardcoded_result

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py:25-32`

**Issue:** The helper `real_y20_square_ratio(m)` short-circuits the `m == 0` case at line 27-28, returning `sp.Integer(1)` without invoking the Gaunt coefficient that the prose says defines the ratio. The downstream assertion at line 167 (`assert_zero("Y20 overlap lane 20", lam20 - 1)`) therefore checks `1 - 1 == 0` by construction, not by Wigner 3-j algebra. The two non-trivial lanes (`m = 1, 2`) are substantive; only the `m = 0` lane is hardcoded.

**Required change:**

Edit lines 25-32. Replace the function body so the m=0 branch uniformly computes the Gaunt ratio (it equals 1 only after the cancellation `gaunt(2,2,2,0,0,0)/gaunt(2,2,2,0,0,0)` actually evaluates), and the same-sign sanity check is skipped only when it is structurally degenerate (m=0 is its own negation, so `gaunt(2,2,2,0,m,m) = gaunt(2,2,2,0,m,-m)` and the check would compare with itself).

BEFORE (lines 25-32):
```
def real_y20_square_ratio(m: int) -> sp.Expr:
    base = sp.simplify(gaunt(2, 2, 2, 0, 0, 0))
    if m == 0:
        return sp.Integer(1)
    same_sign = sp.simplify(gaunt(2, 2, 2, 0, m, m))
    if same_sign != 0:
        raise AssertionError(f"Real-harmonic same-sign cross term should vanish for m={m}: {same_sign}")
    return sp.simplify((sp.Integer(-1) ** m) * gaunt(2, 2, 2, 0, m, -m) / base)
```

AFTER (lines 25-32):
```
def real_y20_square_ratio(m: int) -> sp.Expr:
    base = sp.simplify(gaunt(2, 2, 2, 0, 0, 0))
    if m != 0:
        same_sign = sp.simplify(gaunt(2, 2, 2, 0, m, m))
        if same_sign != 0:
            raise AssertionError(f"Real-harmonic same-sign cross term should vanish for m={m}: {same_sign}")
    return sp.simplify((sp.Integer(-1) ** m) * gaunt(2, 2, 2, 0, m, -m) / base)
```

The change removes the `if m == 0: return sp.Integer(1)` short-circuit and gates only the same-sign sanity check on `m != 0`. The final `return` line is unchanged in form but now also handles `m == 0`, returning `gaunt(2,2,2,0,0,0)/gaunt(2,2,2,0,0,0)` which `sp.simplify` reduces to 1 only after evaluating the nonzero Gaunt coefficient.

**Verification:**
After patching, line 167's assertion `lam20 - 1 == 0` still passes (because the Gaunt ratio of equal numerator and denominator is symbolically 1), but the residual is now obtained by sympy evaluating and dividing the actual Wigner 3-j machinery. A mutation to the `base` definition (e.g., multiplying it by 2) would flip `lam20` to `1/2` and surface the regression — which it would not under the current shortcut. The output transcript label `Y20 overlap lane 20` continues to print with PASS.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py`
- summary: Removed the m=0 shortcut so the real-Y20 ratio is always computed through the Gaunt coefficient ratio.
- deviation: none
