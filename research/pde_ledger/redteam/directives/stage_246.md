---
unit_id: 246
batch: VIII.1
created_at: 2026-06-03T00:00:00-06:00
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-06-03T10:03:59-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 246

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named. Do NOT touch paper.tex, notes/, or any prose documents.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

There are NO paper_misalignment findings; all four are script-side fixes.

---

## F2 — tautological_check (SymPy: quadratic rewrite + σ_min self-test)

**Target:**
- `scripts/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.py:40,54` (quadratic rewrite)
- `scripts/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.py:62-83` (σ_min self-test)

**Issue:** Line 54 asserts `simplify(sigma_y - (1-b+a*y+2*b*y**2)) == 0`, but `sigma_y` was DEFINED on line 40 as exactly that expression — `expand(E)-E==0` is always true and never tests the substitution `y=cos(πx)`. Likewise lines 76-83 recompute `sigma_min_expected` from the SAME piecewise branch logic that defines `sigma_min_piece` (lines 62-65), so the loop only confirms SymPy's `Piecewise` dispatch matches Python's `if` dispatch on identical formulas — it never verifies the piecewise is the actual minimum of σ.

**Required change:**

(a) **Quadratic rewrite — connect σ(x) to the y-quadratic.** Immediately AFTER line 54, add a non-tautological assertion that the original `sigma` equals the quadratic with `y` replaced by `cos(πx)`, using `cos(2πx)=2cos²(πx)-1`:
```python
    # Non-tautological: the actual substitution y = cos(pi*x) must turn sigma into the quadratic.
    assert sp.simplify(sigma - (1 - b + a * sp.cos(sp.pi * x) + 2 * b * sp.cos(sp.pi * x) ** 2)) == 0
```
Do NOT delete line 54; leave it (it documents the y-form) and add this beside it.

(b) **σ_min — assert the piecewise equals the TRUE minimum of the quadratic.** Inside the existing loop (after line 81's `sigma_min_test`, before or replacing the existing assert at line 83), compute the true minimum directly from the quadratic candidates `sigma_y` at `y=±1` and (when admissible) the vertex, WITHOUT reusing the piecewise branch logic, then assert the piecewise matches it:
```python
        # True minimum of the quadratic sigma_y on y in [-1,1], computed independently
        # of the piecewise branch logic (boundary candidates always; vertex only when
        # the parabola opens upward, 2b>0, and the vertex lies in [-1,1]).
        cand = [sp.simplify(sigma_y.subs({a: aval, b: bval, y: 1})),
                sp.simplify(sigma_y.subs({a: aval, b: bval, y: -1}))]
        if bval > 0:
            ystar_val = sp.Rational(-aval, 4 * bval)
            if -1 <= ystar_val <= 1:
                cand.append(sp.simplify(sigma_y.subs({a: aval, b: bval, y: ystar_val})))
        sigma_min_true = sp.simplify(sp.Min(*cand))
        assert sp.simplify(sigma_min_test - sigma_min_true) == 0
```
Keep the existing `assert sp.simplify(sigma_min_test - sigma_min_expected) == 0` at line 83 as well (it is harmless documentation); the new `sigma_min_true` assert is the substantive one.

Note on the three existing test points (do not change them): `(1/2,-1/5)→3/10`, `(5/2,1/4)→-5/4`, `(1/2,1/2)→7/16`. The `sp.Min` construction must reproduce these. If `sp.Rational(-aval, 4*bval)` raises on a non-Rational `aval/bval`, use `sp.Rational(-aval, 4*bval)` only because all three test points are exact `Rational`s; if you prefer, write `ystar_val = -aval / (4 * bval)` and guard with `if ystar_val.is_real and -1 <= ystar_val <= 1`.

**Verification:** SymPy run exits 0; the two new asserts pass for all three test points. The verifier confirms a deliberately wrong vertex coefficient (e.g. `a**2/(7*b)`) would now fail.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.py`
- summary: Added the non-tautological trigonometric quadratic substitution check and an independent finite-candidate true-minimum assertion for the three sigma_min test points.
- deviation: none

---

## F3 — tautological_check (SymPy: transported σ_min branch selection)

**Target:** `scripts/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.py:189-196`

**Issue:** Line 190 sets `sigma_min_transport = simplify(1 + b_r - a_r)` and line 196 asserts it equals `1-(a0-b0)*s_r` — pure algebraic rearrangement of the definition, which cannot fail and never verifies that the Session-I orientation `a0>0,b0<0` actually selects the boundary-minimum branch of `sigma_min_piece` (Section 2).

**Required change:** Tie the transported minimum to `sigma_min_piece` under the Session-I orientation. In Section 8, introduce positively/negatively-signed orientation symbols and substitute them into `sigma_min_piece`, then assert the result equals `1-(a0-b0)*s_r`. Insert after line 196:
```python
    # Verify the Session-I orientation (a0>0, b0<0) selects the boundary-minimum branch
    # of the Section-2 piecewise, rather than asserting the boundary form by hand.
    a0p = sp.symbols("a0p", positive=True)   # stands in for a0 > 0
    b0n = sp.symbols("b0n", negative=True)   # stands in for b0 < 0
    a_r_or = a0p * s_r
    b_r_or = b0n * s_r
    sigma_min_branch = sp.piecewise_fold(
        sigma_min_piece.subs({a: a_r_or, b: b_r_or})
    )
    sigma_min_branch = sp.simplify(sigma_min_branch)
    assert sp.simplify(sigma_min_branch - (1 - (a0p - b0n) * s_r)) == 0
```
If `sigma_min_piece.subs(...)` does not auto-resolve the branch (SymPy may leave the `Piecewise` unevaluated even with signed assumptions), replace the assertion with an explicit branch-condition check instead: assert that the `b<=0` condition is satisfied (`b_r_or` is negative since `b0n<0`, `s_r>0`) and that the boundary piece `1 + b - sp.Abs(a)` evaluates to `1 + b_r_or - a_r_or` under `a_r_or>0`:
```python
    assert (b_r_or < 0) == True or b0n.is_negative  # b(r) < 0 under Session-I orientation
    boundary_piece = (1 + b - sp.Abs(a)).subs({a: a_r_or, b: b_r_or})
    assert sp.simplify(boundary_piece - (1 - (a0p - b0n) * s_r)) == 0
```
Use whichever of the two forms actually runs clean; prefer the first (it exercises the piecewise dispatch). Do NOT remove the existing line-196 assert.

**Verification:** SymPy run exits 0; the new assert references `sigma_min_piece` (or the explicit boundary piece) so a wrong branch boundary `4b` or a sign error in the boundary form surfaces. Output Section 8 printed values unchanged.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.py`
- summary: Added signed Session-I orientation substitutions through `sigma_min_piece` to verify the transported boundary-minimum branch.
- deviation: none

---

## F4 — insufficient_verification (SymPy: pin the S readback)

**Target:** `scripts/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.py:230-234`

**Issue:** `S_eval` (the Session-I shell-loading readback, ≈0.675847711465632 per the card) is computed and printed but never asserted, unlike `g`, `R`, `σ_min`.

**Required change:** After the existing line 232 (`assert abs(float(sigma_min_eval) - (-0.08979545)) < 5e-9`), add:
```python
    assert abs(float(S_eval) - 0.67584771) < 5e-9
```

**Verification:** SymPy run exits 0; new assert pins `S_eval` to the recorded value; output Section 9 `S[sigma](r_eval)` line unchanged.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_sympy_audit.py`
- summary: Added the missing numeric assertion pinning `S_eval` to the recorded Session-I readback value.
- deviation: none

---

## F1 — missing_verification_script (missing_mathematica)

**Target:** `mathematica/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_mathematica_audit.wl`
(NEW file — `.wl` lives in `mathematica/`, NOT `scripts/`.)

**Issue:** No Mathematica engine exists for this non-status, exact-closure unit. Add an INDEPENDENT-route audit using native Mathematica primitives — NOT a transliteration of the `.py`. In particular, verify the minimum claim with `Minimize`/`MinValue` (which the SymPy script does not do), and compute the moment integrals with `Integrate` directly (no by-hand product-to-sum).

**Claim manifest** (each must be verified by a distinct `expectZero`/`expectTrue`/`expectApprox`):

- **M1 — mean preservation:** `Integrate[1 + a Cos[Pi x] + b Cos[2 Pi x], {x, 0, 1}] == 1`.
- **M2 — mouth-bias functional:** `Integrate[(1 + a Cos[Pi x] + b Cos[2 Pi x]) Cos[Pi x/2], {x, 0, 1}]` equals `(2/Pi)(1 + a/3 - b/15)`. Use `FullSimplify[diff == 0]`.
- **M3 — shell-loading functional:** `Integrate[(1 + a Cos[Pi x] + b Cos[2 Pi x]) Cosh[Pi (1 - x)/2]/Cosh[Pi/2], {x, 0, 1}]` equals `(2 Tanh[Pi/2]/Pi)(1 + a/5 + b/17)`.
- **M4 — exact minimum / sign-change (INDEPENDENT of the piecewise):** for each of the three test points `{a->1/2,b->-1/5}`, `{a->5/2,b->1/4}`, `{a->1/2,b->1/2}`, compute `MinValue[{1 + a Cos[Pi x] + b Cos[2 Pi x], 0 <= x <= 1}, x]` (or `Minimize[...]`) and `expectApprox` it against the paper's piecewise values `3/10`, `-5/4`, `7/16` respectively. (This is the check the SymPy script omits — derive the minimum, do not re-encode the piecewise.)
- **M5 — two-moment map:** with `Msrc = {{1/3, -1/15}, {1/5, 1/17}}`, `Det[Msrc] == 14/425`; and `Inverse[Msrc]` reproduces the inverse coefficients — verify by `Inverse[Msrc] . {gt, St}` equals `{(85/42) St + (25/14) gt, (425/42) St - (85/14) gt}` via `FullSimplify`.
- **M6 — quarter-ratio theorem:** with `rF1 = Sqrt[(12/Pi^2)(37/20)^2 - 1]`, `gpm = rF1 ± (1/2) Sqrt[1 + rF1^2]`, verify `(gpm - rF1)^2/(1 + rF1^2) == 1/4` for both signs (`FullSimplify`).
- **M7 — compensation line:** `Solve[(2/Pi)(1 + a/3 - b/15) == gc, b]` gives `b == 5 a + 15 - (15 Pi/2) gc` (`FullSimplify` on the solved branch).
- **M8 — transported sign-change threshold (Session-I orientation):** with `s[r] = rsig^2/(r^2 + rsig^2)`, `a[r] = a0 s[r]`, `b[r] = b0 s[r]`, and assumptions `a0>0, b0<0, rsig>0, r>0`, verify via `Reduce`/`Simplify` that `Minimize[{1 + a[r] Cos[Pi x] + b[r] Cos[2 Pi x], 0<=x<=1}, x]` (or the boundary candidate `1 + b[r] - Abs[a[r]]`) equals `1 - (a0 - b0) s[r]`, and that the threshold solves to `r < rsig Sqrt[a0 - b0 - 1]`.
- **M9 — Session-I numeric readback:** with `a0=2.2, b0=-0.6, rsig=0.8, reval=1.00217028`, `expectApprox` (tol `5*10^-9`): `g==0.82823667`, `S==0.67584771`, `R==0.21677037`, `sigmaMin==-0.08979545`, and `expectTrue` `g[0] > 1` and `rThr > reval`.

Use the project's standard Mathematica audit harness helpers if present (`expectZero[expr_]`, `expectTrue[cond_]`, `expectApprox[lhs_, rhs_, tol_]`); otherwise define minimal versions that `Print` a labeled PASS and `Exit[1]` on failure. Wrap each `FullSimplify` over symbolic `a,b` with no false positivity assumptions (a,b are unrestricted reals — do NOT add `Assuming[a>0,...]` for M2/M3/M5/M6/M7). For M4/M8 use `Minimize`/`Reduce` so the minimum is DERIVED, not asserted.

**Self-test reminders honored in this manifest:**
- No `D[expr, var]` where `var` is absent from `expr` (the variable-independence trap) — M2/M3/M5/M7 use `FullSimplify[lhs-rhs==0]`, not differentiation.
- M2/M3 integrands over `[0,1]` are bounded, finite integrals — no symmetric-domain parity assumption is used; `Integrate` evaluates them directly.
- Trivial-case pre-check: M4 point `(1/2,-1/5)` → boundary min `min(1.3,0.3)=0.3=3/10`; `(5/2,1/4)` → `min(3.75,-0.75)=-1.25=-5/4`; `(1/2,1/2)` → interior vertex `7/16` — these match the `expectApprox` targets.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 246` and confirms the new `.wl` is present, exits 0, and M1-M9 checks appear and pass.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_mathematica_audit.wl`
- summary: Added an independent Mathematica audit covering M1-M9 with direct integrals, `MinValue`, matrix inversion, `Solve`, `Reduce`, and numeric readback checks.
- deviation: none
