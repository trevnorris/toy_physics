---
unit_id: 001
batch: I.1
created_at: 2026-05-25T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 001 (v2, paper-grounded re-audit)

Both findings below are `paper_misalignment`. **Do not edit any file** until the user has chosen a direction for each. The orchestrator is holding for user resolution. Codex must not silently flip signs in either `scripts/`, `mathematica/`, `paper/`, or `notes/` to "fix" the apparent mismatches.

The prior v1 directive (Mathematica transliteration finding) was already applied — that work stands and is not re-opened here. The two findings below are net-new from the paper-grounded v2 re-audit.

## F1 — paper_misalignment

**Subtype:** target_mismatch (sign convention on source coupling in modal-wall PDE)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_001.tex:156-163` quote:
  ```
  \boxed{
  \mu_\eta\,\partial_t^2 q_{lm}
   -\partial_w\!\bigl(T_w\partial_w q_{lm}\bigr)
   +\bigl[K_\eta+l(l+1)T_\Omega\bigr]q_{lm}
   = S_{lm}^{(\psi,A)}+f_{lm}^{\rm ext}}.
  ```
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage001_geometry_lift.md:358-361` quote:
  ```
  mu_eta q_{lm,tt} - partial_w( T_w partial_w q_{lm} )
   + [ K_eta + l(l+1) T_Omega ] q_{lm}
   = S_{lm}^{(psi,A)} + f_{lm}^{ext}.
  ```
- Both sources place the source on the RHS with **positive** sign.

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py:188-192` quote:
  ```python
  source_total = S_lm(t, w) + f_ext(t, w)
  ldens_forced = ldens - q(t, w) * source_total
  el_forced = euler_equations(ldens_forced, q(t, w), [t, w])[0]
  target_forced = target_dens - source_total
  expect_zero("sourced densitized Euler-Lagrange equation", el_forced.lhs - target_forced)
  ```
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl:184-192` quote:
  ```
  sourceTotal = Slm[t, w] + fext[t, w];
  ldensForced = ldens - qField*sourceTotal;
  elForcedEq = EulerEquations[ldensForced, q[t, w], {t, w}];
  elForced = FullSimplify[elForcedEq[[1]] - elForcedEq[[2]], Assumptions -> $Assumptions];
  targetForced = targetDens - sourceTotal;
  expectZero["sourced densitized Euler-Lagrange equation", elForced - targetForced];
  ```
- SymPy's `euler_equations` returns the EL identity in the form `∂L/∂q − Σ_i d/dx_i ∂L/∂(∂_{x_i} q) = 0`. With `ldens_forced = ldens − q · source_total`, the resulting equation is `μ q_tt − ∂_w(T_w q_w) + K_l q = −(S_lm + f_ext)` — the **negative** of the paper's RHS source.

## Resolve before fix_loop

The paper card and the notes both write the modal-wall PDE with `+(S_{lm}^{(psi,A)} + f_{lm}^{ext})` on the RHS. The SymPy and Mathematica scripts encode the source as `−q · source_total` in the Lagrangian, which produces `−(S_lm + f_ext)` on the equation RHS. The script's own summary line claims it verifies the paper's positive-RHS form, but the assertion exercises the opposite sign.

**Which sign is correct?**

Possible directions (the user picks one):

- **(a) Paper is correct (most likely)** → flip the script's source-coupling sign:
  - `scripts/.../sympy_audit.py:189` change `ldens - q(t, w) * source_total` → `ldens + q(t, w) * source_total`
  - `scripts/.../sympy_audit.py:191` change `target_dens - source_total` → `target_dens + source_total`
  - `mathematica/.../mathematica_audit.wl:188` change `ldens - qField*sourceTotal` → `ldens + qField*sourceTotal`
  - `mathematica/.../mathematica_audit.wl:191` change `targetDens - sourceTotal` → `targetDens + sourceTotal`

  Re-run sympy + mathematica; outputs should still show `= 0`.

- **(b) Script is correct** → update the paper card (line 161 of `stage_001.tex`) and the notes (line 361 of the markdown) to read `= -(S_{lm}^{(psi,A)} + f_{lm}^{ext})` on the RHS. (Unlikely — the standard physics convention places source on the RHS with positive sign.)

- **(c) Both are derived from a third source that contradicts both** → flag for deeper review (audit Stage 002+003 which reuse this modal PDE form to identify the upstream convention author).

The orchestrator will not invoke Codex on this unit until the user picks a direction.

## F2 — paper_misalignment

**Subtype:** target_mismatch (gauge-fix term sign in linearized localized Maxwell, contingent on unstated metric signature)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_001.tex:192-198` quote:
  ```
  \boxed{
  \partial_M\bigl(Z(w)\delta F^{MN}\bigr)
   +\frac1\xi\partial^N(\partial\!\cdot\!\delta A)
   =\mu_0\delta J^N}.
  ```
- The card does not state the metric signature. `partial^N` vs `partial_N` differs in sign for spatial indices under mostly-minus vs mostly-plus conventions.

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py:209-229` quote (relevant fragments):
  ```python
  lmax = (
      sp.Rational(1, 2) * Zloc(w) * Fwx**2
      - sp.Rational(1, 2) * divA**2 / gauge_xi
      + mu0 * (Jx(x, w) * Ax(x, w) + Jw(x, w) * Aw(x, w))
  )
  ...
  target_Ax = sp.diff(Zloc(w) * Fwx, w) - sp.diff(divA, x) / gauge_xi - mu0 * Jx(x, w)
  target_Aw = -sp.diff(Zloc(w) * Fwx, x) - sp.diff(divA, w) / gauge_xi - mu0 * Jw(x, w)
  ```
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl:204-211` mirrors with `lmax = (1/2) zloc fwx^2 - divA^2/(2 gaugeXi) + ...` and the same `targetAx`/`targetAw` constructions.

Setting `target_Ax = 0` gives `∂_w(Z F_wx) − (1/ξ) ∂_x(div A) = μ_0 J_x` — the gauge term enters the equation with a **negative** sign. Reading the paper's `+(1/ξ) ∂^N(∂·δA)` under a mostly-plus signature gives `+(1/ξ) ∂_x(div A)`, opposite to script. Under mostly-minus, `∂^x = −∂_x` and the paper's equation rewrites as `−(1/ξ) ∂_x(div A) = ...`, matching script.

## Resolve before fix_loop

The discrepancy is either a genuine script-side sign error or a metric-signature convention difference left unstated in the paper card.

**What metric signature does the (4+1)-dimensional parent theory use?**

Possible directions (the user picks one):

- **(a) Mostly-minus signature, e.g. `(+,−,−,−,−)`** → script is correct as-is; add a one-line metric-signature statement to either the paper card (preferred, near the localized-Maxwell box) or the script docstring. Suggested paper-side text: *"All raised-index spatial derivatives use mostly-minus signature; for spatial N, `partial^N = -partial_N`, so the gauge-fix term reads `-(1/xi) partial_N(partial . delta A)` in coordinate form."* No script change needed.

- **(b) Mostly-plus signature, e.g. `(−,+,+,+,+)`** → script's gauge-fix term sign is wrong. Flip the gauge-fix term in the Lagrangian and corresponding targets:
  - `sympy_audit.py:211` change `- sp.Rational(1, 2) * divA**2 / gauge_xi` → `+ sp.Rational(1, 2) * divA**2 / gauge_xi`
  - `sympy_audit.py:225` change `sp.diff(Zloc(w) * Fwx, w) - sp.diff(divA, x) / gauge_xi - mu0 * Jx(x, w)` → `... + sp.diff(divA, x) / gauge_xi - mu0 * Jx(x, w)`
  - `sympy_audit.py:226` change `-sp.diff(Zloc(w) * Fwx, x) - sp.diff(divA, w) / gauge_xi - mu0 * Jw(x, w)` → `... + sp.diff(divA, w) / gauge_xi - mu0 * Jw(x, w)`
  - `mathematica_audit.wl:205` change `divA^2/(2 gaugeXi)` (subtracted) → `+ divA^2/(2 gaugeXi)`
  - `mathematica_audit.wl:210` change `- D[divA, x]/gaugeXi` → `+ D[divA, x]/gaugeXi`
  - `mathematica_audit.wl:211` change `- D[divA, w]/gaugeXi` → `+ D[divA, w]/gaugeXi`

  Re-run sympy + mathematica; outputs should still show `= 0`.

- **(c) The parent theory's actual convention is neither pure mostly-plus nor pure mostly-minus, or `Z(w)` absorbs a non-trivial wall-direction metric factor** → flag for deeper review against the parent-theory definition document.

The orchestrator will not invoke Codex on this unit until the user picks a direction.
