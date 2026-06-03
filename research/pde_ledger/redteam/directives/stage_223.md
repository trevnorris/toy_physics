---
unit_id: 223
batch: VII.1
created_at: 2026-06-02T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: true
verification_status: pending
needs_user_resolution: false
applied_at: 2026-06-02T16:26:36-06:00
findings_applied: 3
findings_blocked: 0
---

# Codex directive — unit 223

Apply each non-paper_misalignment finding below in order (F2, F3). After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

The `paper_misalignment` finding(s) (F1) have been RESOLVED by the user (2026-06-02) — see the `## RESOLVED` block at the END of this directive and apply the authorized notes-only edit there as part of this fix loop. (Codex applies notes/*.md; Claude reviews.)

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

Do NOT touch paper.tex or the appendix. The ONLY authorized prose/notes edit is the notes-only change specified in the `## RESOLVED` block at the END of this directive (user-authorized 2026-06-02); apply exactly that and make no other notes/prose edits.

## F1 — paper_misalignment

**Subtype:** value_mismatch

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md:351` quote: `| 0.2 | 0.000576970879843 | 29.3158464872314 | 206.814136942081 | 205.502546600713 |`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.py:371` quote: `(0.000576970879843045, 29.3158464872314, 138.814136942081, 137.502546600713),`
- corroborated by saved output line 63: `(0.2, 0.0005769708798430451, 29.31584648723137, 138.81413694208146, 137.50254660071312)`

## Resolve before fix_loop

At $\lambda_W=0.2$ the compatibility-family scan's two wall-like $\mathcal R_Q$ figures disagree: the notes table says 206.814136942081 (lower) / 205.502546600713 (upper); the SymPy script asserts and reproduces 138.814136942081 / 137.502546600713. The first two columns ($P_{0,\rm target,compat}$, $K_{\rm compat}$) agree exactly, and every other scan row (0.4–1.0) agrees on both sides, so this is an isolated discrepancy at one node. The fractional tails are byte-identical and the integer parts differ by exactly 68 on both entries (206→138, 205→137), which is the signature of a one-sided transcription typo rather than an algebraic divergence.

Which value set is correct, 206.81/205.50 (notes) or 138.81/137.50 (script)?

Possible directions (the user picks one):
- (a) Script is correct → the notes table line 351 should read `... 138.814136942081 | 137.502546600713 |` (paper-side edit; Claude applies, Codex does not). No script change.
- (b) Notes are correct → change the script `expected_scan[0]` at line 371 to `(0.000576970879843045, 29.3158464872314, 206.814136942081, 205.502546600713)` and re-run; the script's scan computation must then reproduce those figures (if it does not, the script's scan formula or `rq_func` is wrong and is the real defect).
- (c) Neither — both trace to a third source that contradicts both → flag for deeper review.

F1 is RESOLVED by the user (see the ## RESOLVED block at the end) — apply the authorized notes edit. F2/F3 are independent script-side fixes.

## Applied: F1

- files_changed:
  - `notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md`
- summary: Corrected the λ_W=0.2 lower and upper wall R_Q cells to the verified SymPy values.
- deviation: none

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.py:156-161`

**Issue:** The compatibility-surface asserts at lines 159-161 are circular. Line 156 defines `P0_target_compat = N0*(B4+Z4)/(3*(M+B2+Z2)**2)`, then line 159-160 asserts `cancel(N0/P0_target_compat - 3*(M+B2+Z2)**2/(B4+Z4)) == 0`, which is identically zero by construction (a/(a*b/c) = c/b). Line 161 asserts `cancel(D0_compat - N0/P0_target_compat) == 0` where `D0_compat` (line 157) equals `K_pole - B0 - Z0` and `K_pole` (line 154) is defined as `3*(M+B2+Z2)**2/(B4+Z4) + B0 + Z0`, so `D0_compat = 3*(M+B2+Z2)**2/(B4+Z4) = N0/P0_target_compat` — also guaranteed zero. These verify nothing about the physics; they verify fraction algebra. The card's headline Output (the compatibility surface) is therefore not independently exercised by these dedicated asserts. (The one-pole numerator identity at line 152 is a genuine, non-tautological check and must be left untouched.)

**Required change:**
Make the compatibility surface emerge from equating the two independently-built wall stiffnesses, so the assert can fail if they did not coincide.

Step 1. Keep lines 154-155 (`K_pole`, `K_norm`) as-is — they ARE the two independently-derived stiffnesses (one-pole side and normalization side).

Step 2. Replace the up-front definition at line 156 with a SOLVED quantity. Insert before the current line 156:

```python
    P0_target_compat_solved = sp.solve(sp.Eq(K_norm, K_pole), P0_target)[0]
```

This solves the actual compatibility condition $K_{\rm norm}=K_{\rm pole}$ for the target.

Step 3. Keep the boxed closed form as a separate symbol for comparison (this is the paper's claimed result):

```python
    P0_target_compat = sp.simplify(N0 * (B4 + Z4) / (3 * (M + B2 + Z2)**2))
```

Step 4. Replace the circular assert at line 159-160 with a NON-circular comparison of the solved target against the boxed closed form:

```python
    compatibility_identity = sp.cancel(P0_target_compat_solved - P0_target_compat)
    assert compatibility_identity == 0
```

This now fails if equating $K_{\rm norm}=K_{\rm pole}$ does NOT yield the boxed surface.

Step 5. Keep line 157 (`D0_compat = sp.simplify(K_pole - B0 - Z0)`) and line 161's intent, but it is acceptable to keep line 161 as a downstream consistency note since `D0_compat = N0/P0_target_compat` now follows from a solved (not defined) quantity. Leave line 161 as-is.

Step 6. Leave the primitive-specialization assert (lines 163-166) as-is — it is a definition-inline but harmless and matches paper deliverable (4); do not change it (changing it risks introducing churn without removing circularity that matters).

Do not alter any downstream numeric block; `P0_target_compat` keeps the same value, so all `assert_close` literals (lines 252, 414-417) remain valid.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 223` and confirms a new `sp.solve(sp.Eq(K_norm, K_pole), P0_target)` step appears, the `compatibility_identity` assert now compares solved-vs-closed-form, and the script exits 0 with all checks passing.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.py`
- summary: Replaced the circular compatibility assertion with a solved target from `K_norm == K_pole` compared against the boxed closed form.
- deviation: none

## F3 — missing_verification_script

**Subtype:** missing_mathematica

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_mathematica_audit.wl`

**Issue:** No Mathematica audit exists for unit 223. The unit is `is_checkpoint: False` and `is_status_only_candidate: False`. Under the dual-engine rule a `.wl` is REQUIRED wherever Mathematica CAN independently verify, and it can verify every claim here (integral, rational-function algebra, quartic, root census, numeric slice). Author an independent re-derivation — NOT a transliteration of the SymPy script.

**Anti-transliteration guard:**
- Use native Mathematica primitives: `Integrate` for $\kappa$; `Together`/`Cancel`/`Simplify` for the rational identities; `Solve` (symbolic) and `NSolve`/`Roots` for the quartic; `CoefficientList` for the quartic degree; `N[..., 30]` for slice values.
- Choose a DIFFERENT decomposition than the SymPy script: e.g. build the quartic $F(y)$ directly from the determinant/product form and read its degree via `Exponent[Fy, y]` (not `Poly(...).degree()`); obtain the compatibility surface by `Solve[Knorm == Kpole, P0target]` (the F2-corrected route), not by defining the ratio up front; select wall-like poles by physical magnitude ordering via `Sort` + explicit `TakeLargest[…, 2]` rather than the SymPy `roots[-2:]` slice convention, and cross-check that ordering against the internal/wall labels.
- Verify identities with `expectZero` that strips `ConditionalExpression[0, ...]`; for pole/non-pole checks use `1/expr == 0` form, not `=!= Infinity`.

**Claim manifest** (each must be independently verified):

- **M1** — Overlap constant: $\kappa = \int_0^L \frac{1}{\sqrt L}\sqrt{2/L}\,\sin\!\big(\frac{\pi s}{2L}\big)\,ds = \frac{2\sqrt2}{\pi}$.
- **M2** — One-pole numerator identity: with $D_0=K-B_0-Z_0$, $D_2=-(M+B_2+Z_2)$, $D_4=-(B_4+Z_4)$, $u_2=-D_2/D_0$, $u_4=(D_2^2-D_0D_4)/D_0^2$, show $u_4-4u_2^2 = \dfrac{D_0(B_4+Z_4)-3(M+B_2+Z_2)^2}{D_0^2}$ (exact, symbolic in all couplings).
- **M3** — Compatibility surface: solving $K_{\rm norm}=K_{\rm pole}$ for $P_{0,\rm target}$ (where $K_{\rm pole}=\frac{3(M+B_2+Z_2)^2}{B_4+Z_4}+B_0+Z_0$ and $K_{\rm norm}=\frac{N_0}{P_{0,\rm target}}+B_0+Z_0$) yields $P_{0,\rm target,compat}=\dfrac{N_0(B_4+Z_4)}{3(M+B_2+Z_2)^2}$, equivalently $\dfrac{N_0}{P_{0,\rm target,compat}}=\dfrac{3(M+B_2+Z_2)^2}{B_4+Z_4}$.
- **M4** — Primitive specialization: with $N_0=P^2/\Delta^2$, $B_4=C^2/\varpi^6$, $B_2=C^2/\varpi^4$, the surface equals $\dfrac{(P^2/\Delta^2)(C^2/\varpi^6+Z_4)}{3(M+C^2/\varpi^4+Z_2)^2}$.
- **M5** — Quartic: $F(y)=\big[(K-My)(\varpi^2-y)-C^2\big]\big[(\Omega_U^2-y)(\Omega_W^2-y)-R^2\big]-(\varpi^2-y)\big[G_U^2(\Omega_W^2-y)+2G_UG_WR+G_W^2(\Omega_U^2-y)\big]$ is degree 4 in $y$.
- **M6** — Sample slice (with $\kappa=2\sqrt2/\pi$, $\lambda_B=1/2,\lambda_U=3/10,\lambda_W=2/5,\lambda_R=1/4,\Omega_U=1,\Omega_W=7/5,\varpi=2,M=1,a=c_s=1$): $P_{0,\rm target,compat}\approx0.0020697923180628827$, $K_{\rm compat}\approx24.473754879290977$, $D_{0,\rm compat}\approx24.23730998862225$ (match SymPy `.txt` lines 47-49).
- **M7** — Pole census at $K=K_{\rm compat}$: positive-$\omega$ roots $\approx\{0.971575315129468,\,1.41651290122561,\,1.99753567893361,\,4.94905432364313\}$.
- **M8** — Residue/linewidth $\mathcal R_{Q,*}=\dfrac{27c_s^5}{a^5\omega_*^5 N(\omega_*)}$ with $N(\omega)=\dfrac{[(\Omega_U^2-\omega^2)G_W+RG_U]^2}{[(\Omega_U^2-\omega^2)(\Omega_W^2-\omega^2)-R^2]^2}$: four values $\approx\{0.159888393135835,\,0.000806281535937178,\,30.1999075602499,\,36.1711864832695\}$.
- **M9** — Survival windows (bisection along the $K=K_{\rm compat}(\lambda_W)$ family against thresholds $\mathcal R_Q^{\rm req}(\eta=0.1)\approx21.854566296358396$, $\eta=0.3\approx7.8618736841685335$): lower/upper wall $P_{0,\rm target,compat}$ ceilings $\approx\{0.00283133168555932,\,0.00359651058968466\}$ (10% loss) and $\{0.00817339430971383,\,0.0116633929790174\}$ (30% loss).
  - Note: for the $\lambda_W$ scan column values, hold pending the F1 user resolution at $\lambda_W=0.2$; the other scan rows (0.4–1.0) and the survival-window thresholds are unaffected by F1 and should be verified now.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 223` and confirms the new `.wl` exists at the exact target path, exits 0, exercises M1-M9, and its numeric values agree with the SymPy `.txt` (modulo the F1-pending $\lambda_W=0.2$ scan row).

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_mathematica_audit.wl`
- summary: Added an independent Mathematica audit covering M1-M9 with exact identities, the compatibility branch, pole census, residue values, scan rows, and survival-window bisections.
- deviation: none

## RESOLVED — F1 paper_misalignment (USER-AUTHORIZED 2026-06-02)

Direction: correct the notes to match the verified SymPy script (authoritative; output/*.txt). Notes-only; Codex applies, Claude reviews. Do NOT change the script, paper.tex, or appendix.
- In `notes/stages/moving_throat_pde_stage223_..._sympy_audit.md`, the λ_W=0.2 lower/upper wall R_Q cells read `206.814136942081` and `205.502546600713`; correct them to `138.814136942081` and `137.502546600713` (a +68 slip on each; fractional tails unchanged).
- Acceptance: `206.814136942081` and `205.502546600713` no longer appear in the notes; the two corrected values are present. Append `## Applied: F1`.
