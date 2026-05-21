---
unit_id: 015
batch: I.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-21T18:30:15Z
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
---

# Audit unit 015 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.txt`
- mathematica output: (missing)

## What the script claims to verify

The script audits step 13 ("parent throat action master") by checking five algebraic structures. First, a generic and a concrete (Gaussian-weighted) integration-by-parts identity for the boundary discharge of an `-A η η_w` cross term against `-A η^2/2` and `(dA/dw) η^2/2`. Second, the quadratic limit of the promoted wall action `L = mu0 Rt^2/2 - Tw Rw^2/2 - TO0 eps^2 grad2/2 - U` reduces — after substituting the IBP image `-TwR0 R0p η η_w → +d_TwR_R0p η^2/2` — to the canonical Lagrangian with the closed-form effective mass `K_eta = URR0 - d_TwR_R0p + TwRR0 R0p^2/2`. Third, the wall-only specialization of the full even-gate combinations `K1 = D21 + D01/9` and `H_even = D41 - (2/3) D21 - D01/27` (with `B0_l = Z_l = 0`) collapses to `K1 = -dM + dK/9` and `H_even = (2/3) dM - dK/27`, and the corresponding 2x2 Jacobian in `(dK, dM)` has determinant `1/27` (and detects perturbations of the `1/9` coefficient). Fourth, the real-Y20 squared overlap ratios are `lam20 = 1, lam21 = 1/2, lam22 = -1` via Gaunt coefficients. Fifth, the grouped (weight 1,2,2) trace `xbar = (x20 + 2 x21 + 2 x22)/5` is invariant `x0`, and the `b = 3a` projection identity holds for the (a, b) duo computed from the same group weights. There is no Mathematica companion.

## Assertion inventory

| #   | Script | Line       | Form                                                                                  | Anchored to claim? |
|-----|--------|------------|---------------------------------------------------------------------------------------|--------------------|
| A1  | sympy  | 48-51      | `assert_zero -A eta eta_w - (d/dw(-A eta^2/2) + A' eta^2/2)`                          | partial (product rule + guarded by A2) |
| A2  | sympy  | 52-55      | `assert_nonzero` mutated-sign IBP residual                                            | yes (sign guard)   |
| A3  | sympy  | 58-59      | `assert_nonzero boundary_value(atan(w), w)` (= pi)                                    | yes (helper sanity) |
| A4  | sympy  | 69         | `assert_zero boundary_value(-exp(-2 w^2)/2, w)`                                       | yes (Gaussian decay) |
| A5  | sympy  | 70-73      | `assert_zero ∫(-A eta eta_w) - (boundary + ∫A' eta^2/2)`                              | yes (concrete IBP)  |
| A6  | sympy  | 94         | `assert_zero ∂^2 L2_raw/∂eta ∂eta_w + TwR0 R0p`                                       | yes                |
| A7  | sympy  | 98         | `assert_zero L2_after_ibp_derived - canonical_L2`                                     | yes (K_eta formula) |
| A8  | sympy  | 99-102     | `assert_nonzero` mutated K_eta sign                                                   | yes (sign guard)   |
| A9  | sympy  | 127        | `assert_zero K1_wall - (-dM + dK/9)`                                                  | partial (typed-form typo guard) |
| A10 | sympy  | 128        | `assert_zero H_even_wall - (2/3 dM - dK/27)`                                          | partial (typed-form typo guard) |
| A11 | sympy  | 129-132    | `assert_zero K1_wall.subs(dK→dK_overlap, dM→dM_overlap) - (-dM_overlap + dK_overlap/9)` | no (tautological)  |
| A12 | sympy  | 133-136    | `assert_zero H_even_wall.subs(...) - (2/3 dM_overlap - dK_overlap/27)`                | no (tautological)  |
| A13 | sympy  | 143        | `assert_zero det(wall_matrix) - 1/27`                                                 | yes                |
| A14 | sympy  | 145        | `assert_zero wall_even_solve[dK]`                                                     | yes (follows from det≠0) |
| A15 | sympy  | 146        | `assert_zero wall_even_solve[dM]`                                                     | yes (follows from det≠0) |
| A16 | sympy  | 148        | `assert_nonzero wall_even_solve_perturbed[dK]`                                        | yes (perturbation guard) |
| A17 | sympy  | 149        | `assert_nonzero wall_even_solve_perturbed[dM]`                                        | yes (perturbation guard) |
| A18 | sympy  | 162        | `assert_nonzero wall_det_shift` (= 2 eps/3)                                           | yes (coefficient guard) |
| A19 | sympy  | 167        | `assert_zero lam20 - 1`                                                               | no (hardcoded shortcut) |
| A20 | sympy  | 168        | `assert_zero lam21 - 1/2`                                                             | yes (Gaunt-derived) |
| A21 | sympy  | 169        | `assert_zero lam22 + 1`                                                               | yes (Gaunt-derived) |
| A22 | sympy  | 176        | `assert_zero xbar - x0`                                                               | yes                |
| A23 | sympy  | 177        | `assert_zero bx - 3 ax`                                                               | yes                |

## Findings

### F1 — missing_verification_script

**Severity:** high
**Files:**
- `(missing)` — no Mathematica script for unit 015

**What's wrong:**
The unit's manifest entry has `is_status_only_candidate: False`, so a second-engine companion is required. Only a SymPy script exists at `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py`; there is no matching `.wl` under `mathematica/`. Consequently `engines_agree = n/a`. Every claim of unit 015 — the IBP closure, the closed-form `K_eta = URR0 - d_TwR_R0p + TwRR0 R0p^2/2`, the wall-only even-gate Jacobian determinant `1/27`, the Gaunt-driven Y20 overlap weights `(1, 1/2, -1)`, and the grouped trace/`b = 3a` identity — is verified by a single algebra system.

**Why this matters:**
The closed-form effective mass `K_eta` and the even-gate Jacobian determinant `1/27` are quoted forward and consumed downstream. With only sympy, a wrong sign in `K_eta` or a coefficient slip in the gate Jacobian (e.g., `1/9` vs `1/3` vs `2/9` for the `D01_full` coefficient inside `K1_full`, or `1/27` vs `2/27` for `D01_full` inside `H_even_full`) would not be caught by an independent engine.

**Required change:**
Create `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl` that independently re-derives — not transliterates — the unit's claims using Mathematica primitives (`D`, `Series`, `Coefficient`, `Solve`, `Det`, `FullSimplify`, `Integrate`, `ThreeJSymbol` / `Gaunt`). It must terminate with an explicit pass/fail (`If[FullSimplify[...] =!= 0, Print["FAIL"]; Exit[1]]`) for each numbered claim. The Mathematica script must NOT mirror the sympy code line-for-line; it must derive `K_eta` by computing the second `eps`-derivative of `L` via Mathematica's `Series` (not `D[..., {eps, 2}]/2`), substitute the IBP image symbolically, and read the `eta^2` coefficient via `Coefficient`. Gaunt overlaps should use `ThreeJSymbol` directly (not `Gaunt`) and reduce via the standard `(2l+1)(2l'+1)(2l''+1)/(4 pi)` prefactor explicitly.

**Verification:**
After Codex creates the `.wl`, the verifier runs `redteam exec-mathematica 015` and expects (i) the file to exist under `mathematica/`, (ii) the script to exit 0, (iii) all eight numbered claims in the directive's claim manifest to print PASS lines.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py:129-132`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py:133-136`

**What's wrong:**
By the time the script reaches line 129, the previous assertion at line 127 has already established `K1_wall == -dM + dK/9` (and 128 established `H_even_wall == (2/3) dM - dK/27`). Line 129-132 then asserts:
```
K1_wall.subs({dK: dK_overlap, dM: dM_overlap}) - (-dM_overlap + dK_overlap / 9) == 0
```
Because `K1_wall` is literally the symbolic expression `-dM + dK/9`, applying the substitution `{dK: dK_overlap, dM: dM_overlap}` produces, by linear substitution, exactly `-dM_overlap + dK_overlap/9`. The right-hand side is the same expression typed out. The residual is `0` by sympy's substitution mechanics alone — no integral, no Gaussian moment, no physics. The same pattern repeats verbatim at line 133-136 for `H_even_wall`. The labels in the output (`wall-only K1 from overlap-generated slots`, `wall-only H_even from overlap-generated slots`) imply the script verifies a coupling between the overlap integrals (`dK_overlap`, `dM_overlap` defined at line 110-114) and the wall-only gate forms, but the integrals are never actually evaluated, simplified against any concrete `beta(w)` profile, or compared to anything beyond a rename.

**Why this matters:**
A reader of the output transcript sees two PASS lines that appear to anchor the overlap-integral generators against the gate-coefficient algebra. They anchor nothing — the overlap integrals enter only as variable names. If the overlap-generated story (`dK_overlap`, `dM_overlap`) had a sign flip, or a wrong coefficient on `delta_TO`/`delta_Keta`, or a `delta_Tw beta^2` instead of `delta_Tw (d beta/dw)^2`, this test would not detect it. Two assertions consume audit budget while exercising zero physics, and they sit in a stage that is the only verification of the wall-only gate collapse.

**Required change:**
Replace lines 129-136 with a check that exercises a concrete profile. Concretely:

(a) At line 104 the script declares `wall_w = sp.symbols("wall_w", real=True)` and `beta = sp.Function("beta")(wall_w)`. Introduce a concrete real-valued profile `beta_concrete = sp.exp(-wall_w**2)` and four concrete radial-shift profiles, e.g. `delta_mu_concrete = sp.exp(-wall_w**2)`, `delta_Tw_concrete = sp.exp(-wall_w**2)`, `delta_TO_concrete = sp.exp(-wall_w**2)`, `delta_Keta_concrete = sp.exp(-wall_w**2)`. (Choose Gaussians to keep the integrals closed-form.)

(b) Evaluate `dM_overlap_concrete = sp.integrate(delta_mu_concrete*beta_concrete**2, (wall_w,-sp.oo,sp.oo))` and `dK_overlap_concrete = sp.integrate(delta_Tw_concrete*sp.diff(beta_concrete, wall_w)**2 + (delta_Keta_concrete + 6*delta_TO_concrete)*beta_concrete**2, (wall_w,-sp.oo,sp.oo))`.

(c) Assert that the wall-only `K1_wall.subs({dK: dK_overlap_concrete, dM: dM_overlap_concrete})` equals the concretely-substituted closed form `-dM_overlap_concrete + dK_overlap_concrete/9`. This still has the substitution-tautology issue at the symbolic level, but now both sides are concrete numbers, so any wrong coefficient or integrand sign in lines 110-114 propagates into a concrete numerical mismatch.

(d) Additionally, mutate the `(delta_Keta + 6*delta_TO)*beta^2` term by changing the `6` to `5` and assert that the wall-only K1 (using the mutated `dK_overlap`) no longer matches the unmutated `dK_overlap_concrete/9` form. This guards the `6 delta_TO` coefficient.

After patching, lines 129-136 are removed and the new concrete-Gaussian block (with the mutation guard) sits in their place.

**Verification:**
The output transcript prints two new labels such as `wall-only K1 from concrete overlap integrals` and `wall-only K1 detects mutated 6 delta_TO coefficient`, both with PASS, and the labels `wall-only K1 from overlap-generated slots` and `wall-only H_even from overlap-generated slots` no longer appear (or have been re-purposed to fail under the same mutation).

### F3 — hardcoded_result

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py:25-32` (function `real_y20_square_ratio`)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py:167` (assertion `lam20 - 1 == 0`)

**What's wrong:**
The function `real_y20_square_ratio` defined at line 25-32 short-circuits the m=0 case:
```python
def real_y20_square_ratio(m: int) -> sp.Expr:
    base = sp.simplify(gaunt(2, 2, 2, 0, 0, 0))
    if m == 0:
        return sp.Integer(1)
    ...
```
When `m == 0`, the function returns the literal `Integer(1)` without computing any Gaunt coefficient — note the line above already evaluates `base = gaunt(2,2,2,0,0,0)`, but the returned value ignores it entirely. The subsequent assertion at line 167, `assert_zero("Y20 overlap lane 20", lam20 - 1)`, then reads `1 - 1 = 0` regardless of what the actual Gaunt structure is. For `m = 1, 2` the function does compute `(-1)^m gaunt(2,2,2,0,m,-m) / base`, so lines 168-169 are substantive. Only the `lam20 = 1` lane is hardcoded.

**Why this matters:**
The unit's claim is that the three real-Y20 overlap ratios `(lam20, lam21, lam22) = (1, 1/2, -1)` follow from Wigner 3-j algebra. The `lam20 = 1` lane is a coincidence of the normalization choice (dividing by `base = gaunt(2,2,2,0,0,0)`), not a free physics input — but the script never executes that division for `m=0`. If a future refactor swaps the normalization choice (e.g., divides by `(-1)^m gaunt(2,2,2,0,0,0)` or by some other lane), the `m=0` shortcut would silently continue to return 1 while the other lanes would shift relative to a different reference. The substantive content (`m=1: 1/2`, `m=2: -1`) is checked; the trivial m=0 lane is not.

**Required change:**
Remove the m=0 short-circuit. Edit lines 25-32 so the function uniformly computes the Gaunt ratio for all m, with the same-sign sanity check skipped only when it is structurally degenerate (m=0 is its own negation). Concretely:
```python
def real_y20_square_ratio(m: int) -> sp.Expr:
    base = sp.simplify(gaunt(2, 2, 2, 0, 0, 0))
    if m != 0:
        same_sign = sp.simplify(gaunt(2, 2, 2, 0, m, m))
        if same_sign != 0:
            raise AssertionError(
                f"Real-harmonic same-sign cross term should vanish for m={m}: {same_sign}"
            )
    return sp.simplify((sp.Integer(-1) ** m) * gaunt(2, 2, 2, 0, m, -m) / base)
```
For m=0 this returns `gaunt(2,2,2,0,0,0)/gaunt(2,2,2,0,0,0)` which simplifies to 1 only after the Gaunt coefficient is actually evaluated and divided. The assertion at line 167 then exercises the Gaunt machinery rather than the literal `Integer(1)` shortcut.

**Verification:**
After patching, line 167's `lam20 - 1` is residual-zero only because sympy's `gaunt(2,2,2,0,0,0)` actually computes a nonzero rational and divides cleanly. The output transcript still prints `Y20 overlap lane 20` with PASS, but a mutation to the `base` definition (e.g., replacing `gaunt(2,2,2,0,0,0)` with `2*gaunt(2,2,2,0,0,0)`) would now flip `lam20` to `1/2` and surface the regression.

## Independent-derivation check (Mathematica)

No `.wl` exists. See finding F1. There is no second engine to compare against.

## Engine cross-check

`engines_agree = n/a` — only sympy is present.

## Verdict justification

The SymPy script's substantive content is the closed-form `K_eta = URR0 - d_TwR_R0p + TwRR0 R0p^2/2` (anchored by A6-A8, including the mutated sign guard at A8), the wall-only Jacobian determinant `1/27` (A13) with two perturbation guards (A16, A17, A18), the Gaunt-derived `lam21 = 1/2`, `lam22 = -1` (A20, A21), and the grouped-trace / `b=3a` identities (A22, A23). These hold up under attack: I tried to break the K_eta sign by recomputing the `eps^2` coefficient of `L = mu0 Rt^2/2 - Tw Rw^2/2 - TO0 eps^2 grad2/2 - U` by hand and the script's `L2_raw - cross_term + cross_after_ibp` correctly yields `mu0 eta_t^2/2 - Tw0 eta_w^2/2 - TO0 grad2/2 - (URR0 - d_TwR_R0p + TwRR0 R0p^2/2) eta^2/2`. I tried to find a sign-error path in the Jacobian: with `K1_wall = -dM + dK/9` and `H_even_wall = (2/3) dM - dK/27`, the 2x2 matrix's columns are `[1/9, -1/27]^T` and `[-1, 2/3]^T`, giving det `(1/9)(2/3) - (-1)(-1/27) = 2/27 - 1/27 = 1/27` — the script's value. I checked the grouped-trace coefficients: `(2*1 - 1/2 - (-1))/10 = (2.5)/10 = 1/4 = a`, `(1/2 - (-1))/2 = 3/4 = b`, and `b - 3a = 0`. All substantive checks pass.

The three findings are: (1) the Mathematica engine is absent, which is required for `is_status_only_candidate: False`; (2) the wall-only-overlap-from-slots assertions at lines 129-136 reduce to symbolic substitution renames and exercise no physics or integral structure; (3) the `lam20 = 1` lane is a literal shortcut in the `real_y20_square_ratio` helper that never invokes the Gaunt coefficient it claims to. The unit is repairable; verdict is `findings`, not `stop_cold`. None of the corrections would propagate a sign or coefficient change to downstream units (the substantive results F1 must reproduce are unchanged).

## Self-test notes

I checked four traps. Variable independence: F1's required mathematica derivation uses `D[L, eps, 2]` on a polynomial in eps where every named symbol is independent of eps after the eps-expansion of Tw, U, Rw, Rt is performed, so each `D` exercises a real channel; no derivative is identically zero by construction. Symmetry/parity: F2's proposed Gaussian profiles `beta = exp(-w^2)` and `delta_X = exp(-w^2)` are all even in `wall_w`, so `delta_Tw (d beta/dw)^2` integrand is `exp(-w^2) * 4 w^2 exp(-2 w^2)` = even, and `(delta_Keta + 6 delta_TO) beta^2 = 7 exp(-3 w^2)` = even — both nonzero finite Gaussian moments; the mutation guard (changing 6 to 5) shifts the integral by `-exp(-3 w^2)` integrated, again nonzero. Trivial-case pre-check: for F3, with `gaunt(2,2,2,0,0,0)` being a known nonzero rational (Wigner 3-j with l1=l2=l3=2, m's all zero is sqrt(5/(14 pi)) * (-2/sqrt(35)) or similar nonzero closed form), `gaunt(2,2,2,0,0,0)/gaunt(2,2,2,0,0,0) = 1` symbolically; the assertion `lam20 - 1 == 0` still passes. Path specifications: F1's target path is `mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl` — under `mathematica/`, with `_mathematica_audit.wl` suffix, matching the existing `.py` filename pattern. No self-test step uncovered an error.
