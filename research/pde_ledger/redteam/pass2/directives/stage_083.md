---
unit_id: 083
batch: III.4
created_at: 2026-06-05T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-05T21:39:03Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 083

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py:57-93`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl:62-103, 137-154`

**Issue:**
The `Delta_0`/`Delta_inf` "defining-equation residual" checks are tautological. Both scripts define `Delta = numer/denom` (py:57-60, wl:55-60) and then assert `denom·Delta − numer == 0`, which is identically 0 by construction (`denom·(numer/denom) − numer ≡ 0`) for any closed form typed in — it cannot catch a Cosh↔Sinh paste error in the closed form. The SymPy comment (py:68-70) falsely calls this "non-tautological," and the Mathematica comment (wl:62-71) advertises an "Independent BVP derivation" that the code never performs. Two further `.wl` checks are tautologies of the same flavor: `A_F1 independent vs closed-form` (wl:102-103) compares `aF1` against `aF1Indep` that differs only by writing `(Pi/2)^2` for `Pi^2/4` (identical), and `Omega(Pe) identity residual` (wl:137-154) is `omega·denom − numer` with `omega := numer/denom`. The only real anchor on the closed forms is the Mathematica numeric battery (wl:169-170); SymPy has no numeric anchor at all.

**Required change:**

SymPy (`scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py`):
1. Add a numeric-tolerance helper near `expect_zero` (after py:30):
   ```python
   def expect_close(name: str, value, target, atol: float) -> None:
       diff = abs(sp.N(value, 40) - sp.N(target, 40))
       print(f"{name} diff = {diff}")
       if diff > sp.Float(atol):
           raise AssertionError(f"{name}: diff {diff} exceeds {atol}")
   ```
2. Replace the two tautological residual checks (py:66-81). Either delete the residual blocks and the false comment (py:68-70), OR keep them only as informational `print`s (no `expect_zero`). In their place, add NON-tautological numeric anchors against the external literals the notes report (notes:81,83):
   ```python
   expect_close("Delta_0(F1) numeric anchor", Delta0_F1, sp.Float("1.73302079021525e-4"), 1e-16)
   expect_close("Delta_inf(F1) numeric anchor", DeltaInf_F1, sp.Float("2.01447565540522e-2"), 1e-15)
   ```
   (These literals are exactly what the Mathematica `.wl` already pins at wl:169-170 and what notes:81/83 report — confirmed against the saved sympy output values 1.733020790…e-4 and 2.0144756554…e-2.) Fix or delete the comment at py:68-70 so it no longer claims non-tautology.

Mathematica (`mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl`):
3. The numeric anchors at wl:169-170 already provide real coverage and should STAY. The defect is the misleading comment block. EITHER (preferred, minimal) delete/correct the "Independent BVP derivation" comment (wl:62-71) so it accurately describes the residual as a redundant structural check rather than a BVP derivation; AND delete the tautological `A_F1 independent vs closed-form` check (wl:102-103) and the "independent" framing of its comment (wl:96-99), since `(Pi/2)^2 ≡ Pi^2/4` makes `aF1Indep` byte-identical to `aF1`. OR (if you prefer to keep an independence claim) replace wl:77-91 with a genuine `DSolve` of `u''[s] - alphaF1^2 u[s] == 0` under the Robin wall + unit-edge boundary conditions and `expectZero` the extracted boundary values against `delta0F1`/`deltaInfF1`. Do NOT silently leave the false "Independent BVP derivation" comment in place.
4. Correct the `Omega(Pe) identity residual` comment (wl:139-146): it already half-admits "this is the definition we typed in"; tighten it to state plainly that this is a structural sanity print, not an independent verification (the residual is tautological). The check itself may stay as a print/`expectZero` but the comment must not overstate.

**Self-test confirmation (auditor):**
- Variable independence: no `diff` in the new anchors; the perturbation test (swap `cosh`→`sinh` in `Delta0_F1`) would now flip the new numeric anchor from pass to fail — non-tautological. PASS.
- Trivial-case: `Delta0_F1` evaluated numerically gives 1.7330207902…e-4 (saved output py:15), which matches the target literal to <1e-16; `DeltaInf_F1` gives 2.0144756554…e-2, matches to <1e-15. The tolerances (1e-16, 1e-15) are exactly the `.wl`'s wl:169-170 and the saved `.wl` diffs were 1.6e-19 and 4.2e-17, both inside tolerance. PASS.
- Paper round-trip: literals come from notes:81/83; no new paper_misalignment. PASS.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 083` and `redteam exec-mathematica 083` and confirms: the new SymPy numeric anchors appear and the script exits 0; the SymPy script no longer contains the comment claiming the residual is "non-tautological"; the `.wl` no longer contains the unfulfilled "Independent BVP derivation" claim (or actually performs the BVP). Both scripts still exit 0 with all checks passing.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl`
- summary: Added SymPy numeric anchors for `Delta_0`/`Delta_inf`, removed misleading non-tautology/independence framing, and deleted the tautological Mathematica `A_F1` comparison.
- deviation: none

## F2 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py:111-162`

**Issue:**
The stage's `\stagefield{Output}` is the direct operator-selected support window (the four `Pe`, four `zeta`, four `Pi/C_mix` values and the ceiling). In the SymPy script every one of these is `print`ed but never asserted (py:111-122, 152-162). The SymPy engine therefore never verifies the stage's headline deliverable; only the Mathematica `.wl` pins them (wl:171-179). The SymPy "audit" listed in the card's `\stagefield{Verification}` confirms only the derivative identity, the `y_F1` root, and a 4-point monotonicity sample.

**Required change:**
Using the `expect_close` helper added in F1, add numeric assertions for the deliverable windows in the SymPy script, mirroring the Mathematica battery (wl:171-179) and the notes values (notes:111-137). Insert after the corresponding prints (after py:122 for the `Pe` block, after py:156 for the `zeta` block):
```python
expect_close("Pe_-^(chi) numeric check", Pe_minus_chi, sp.Float("96.5285247264385"), 1e-10)
expect_close("Pe_+^(chi) numeric check", Pe_plus_chi, sp.Float("11220.5441626259"), 1e-7)
expect_close("Pe_-^(J) numeric check",   Pe_minus_J,  sp.Float("22.0062226330754"), 1e-10)
expect_close("Pe_+^(J) numeric check",   Pe_plus_J,   sp.Float("2558.01892349205"), 1e-8)
...
expect_close("zeta_-^(chi) numeric check", zeta_minus_chi, sp.Float("2.46622291347846"), 1e-12)
expect_close("zeta_+^(chi) numeric check", zeta_plus_chi, sp.Float("2.46752913273870"), 1e-12)
expect_close("zeta_-^(J) numeric check",   zeta_minus_J,   sp.Float("2.44257571477179"), 1e-12)
expect_close("zeta_+^(J) numeric check",   zeta_plus_J,    sp.Float("2.46752736855058"), 1e-12)
expect_close("zeta_max^(F1) numeric check", zeta_max_F1,   sp.Float("2.46752922945601"), 1e-12)
```
Tolerances copy the `.wl` (wl:171-179). The `Pi/C_mix = 1 + zeta` values are then implied by the zeta checks; an explicit assert on them is optional and not required.

**Self-test confirmation (auditor):**
- Trivial-case: each target equals the saved SymPy output value (py:21-33) to full printed precision, so all new asserts pass as-is; perturbing `Theta_chi_coeff` (py:106) by 1e-3 would shift `Pe_-^(chi)` well outside 1e-10, flipping the assert — non-tautological. PASS.
- Paper round-trip: every literal is taken from notes:111-137; introduces no new paper_misalignment. PASS.
- Path: edit is entirely within `scripts/...sympy_audit.py`. PASS.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 083` and confirms the nine new window checks appear in the transcript and the script exits 0; `grep -c "expect_close" sympy_audit.py` is ≥ 11 (2 Delta anchors from F1 + 9 window checks).

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py`
- summary: Added SymPy numeric assertions for the four `Pe` window endpoints and five `zeta` window/ceiling values.
- deviation: none
