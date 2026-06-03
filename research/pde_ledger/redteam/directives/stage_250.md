---
unit_id: 250
batch: VIII.1
created_at: 2026-06-03T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-03T12:58:45-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 250

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named. Do NOT touch paper.tex, notes/, or any prose document. There is no `paper_misalignment` here — every numeric literal already matches the notes/appendix; do not change any literal.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

## F1 — missing_verification_script (missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_mathematica_audit.wl` (NEW file — must live in `mathematica/`, filename exactly as given)

**Issue:** Stage 250 has only a SymPy audit; the unit is `is_checkpoint: false`, `is_status_only_candidate: false`, so the dual-engine rule requires a Mathematica audit wherever Mathematica CAN independently verify the stage. Every claim here is algebraic/symbolic and within native Mathematica primitives. Write a NEW, INDEPENDENT-route `.wl` — do NOT transliterate the `.py`. Specifically, derive the edge with `Solve`/`Reduce`, and prove the window and monotonicity GLOBALLY with `Reduce`/`Resolve[ForAll]` rather than evaluating at sample points (a deliberately different decomposition from the SymPy script).

**Required change:**
Create the `.wl` with checks covering the claim manifest below. Use the project's existing `.wl` helper idioms (e.g. `expectZero[expr_]`, `expectTrue[cond_]`, `expectApprox[a_,b_,tol_]`) consistent with neighbouring audits such as `mathematica/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl`; if you introduce helpers, strip `ConditionalExpression[0, ...]` from `Solve`/`Reduce` outputs before comparing, and prefer `Reduce`/`Resolve` over `Limit` for global/pole statements. Declare all physical parameters positive (`m_s, lambda_eff, g_UV, chi_peak, mu_eta, alpha > 0`) and the energy gap `x = E - Vpeak > 0`.

**Claim manifest** (the new `.wl` must independently verify each):

- **M1 — crossing-time form & monotonicity.** Define `tcross = lambda_eff*Sqrt[m_s/(2 x)]` with `x = E - Vpeak`. Verify `D[tcross, E] < 0` for all `x>0` via `Resolve[ForAll[{...}, params>0 && x>0, D[tcross,E] < 0]]` returning `True`. (LaTeX: `t_cross = lambda_eff sqrt(m_s/(2(E−Vpeak)))`, `dt_cross/dE < 0`.)
- **M2 — collapse time.** `tcollapse = Sqrt[mu_eta/(g_UV chi_peak)]`; confirm `Gamma_coll = 1/tcollapse = Sqrt[g_UV chi_peak/mu_eta]` via `expectZero[Gamma_coll - Sqrt[g_UV chi_peak/mu_eta]]`.
- **M3 — survival ratio & EDGE from Solve.** `S = tcross/tcollapse`. Solve `S^2 == 1` for `E` with `Solve`/`Reduce`; confirm the unique root equals `Vpeak + (m_s/mu_eta) g_UV chi_peak lambda_eff^2/2` via `expectZero`. (LaTeX: `E_edge = Vpeak + (m_s/mu_eta) g_UV chi_peak lambda_eff^2 / 2`.)
- **M4 — global one-sided window (the headline claim).** Prove `Resolve[ForAll[{E}, params>0 && E>Vpeak, (S<1) == (E > Eedge)]]` is `True` — i.e. the safe set `S(E)<1` is exactly the half-line `E>E_edge`. This is the independent-route counterpart to the SymPy single-point check; it must use `Reduce`/`Resolve`, NOT sampling. (LaTeX: `S(E)<1 ⟺ E>E_edge`.)
- **M5 — heavy-throat cancellation.** Substitute `mu_eta -> alpha m_s` into `Eedge`; confirm `D[Eedge /. mu_eta->alpha m_s, m_s] == 0` via `expectZero`, and that the result equals `Vpeak + g_UV chi_peak lambda_eff^2/(2 alpha)`. (LaTeX: `∂E_edge/∂m_s = 0` when `mu_eta = alpha m_s`.)
- **M6 — speed-space identity.** With `Elaunch = (1/2) m_s v0^2 + V0` and `vcrit^2 = 2(Vpeak−V0)/m_s`, confirm `2 (Eedge − V0)/m_s == vcrit^2 + lambda_eff^2 g_UV chi_peak/mu_eta` via `expectZero` (FullSimplify the difference to 0). (LaTeX: `v_safe,min^2 = v_crit,new^2 + lambda_eff^2 g_UV chi_peak / mu_eta`.)
- **M7 — Session-III benchmark numerics (paper-anchored regression).** Using the notes §7-§8 values `Vpeak=3.42933112, V0=0.19999794, lambda_eff=0.42826825, g_UV=0.95, chi_peak=21.73204372, m_s=mu_eta=1836.15267343`, and raw-width `lambda=1, chi_raw=50.74399964`, confirm with `expectApprox`: `t_collapse≈9.43066476`, `E_safe,min≈5.32265943`, `v_safe,min≈0.07469791`, `v_safe/v_crit≈1.25948037`, `t_collapse_raw≈6.17163516`, `E_safe,min_raw≈27.53273095`, and `S(E_edge)=1`, `E_safe,min_raw > E_safe,min`. Do NOT alter these literals — they match the published notes/appendix exactly.

**Self-test (variable independence / strictness, already checked by auditor):** M1 and M4 reference `E`/`x` explicitly, so the `Resolve[ForAll]` statements are NOT vacuous identically-zero-derivative tests. M5's `D[..., m_s]==0` is non-vacuous because `Eedge` depends on `m_s` until `mu_eta->alpha m_s` cancels it. M4 is a strict global window statement (`Reduce`/`Resolve`), not a sample-point check.

**Verification command:**
The verifier will run `redteam exec-mathematica 250` and confirm the new `.wl` exists, contains the M1-M7 checks, and exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_mathematica_audit.wl`
- summary: Added an independent Mathematica audit covering M1-M7 with global Resolve/ForAll monotonicity and window proofs, Solve/Reduce edge derivation, symbolic identities, and Session-III numeric regressions.
- deviation: none

## F2 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.py:53-93`

**Issue:** The paper's headline deliverable — the one-sided Goldilocks window with a unique lower edge (notes §3.1 "at most one lower survival edge", §3.2; appendix line 303; card line 67) — rests on the GLOBAL strict monotonicity `dS/dE < 0` for all `E > Vpeak` (and `dt_cross/dE < 0`, notes §1.2). The script only PRINTS `dS_dE` (line 83) and `dtcross_dE` (line 57) and checks the strict window inequality at a single sampled energy (`Smax_num < 1.0`, line 190) plus the limit `S_inf==0` (line 93). It never asserts the sign of either derivative globally, so a hypothetical non-monotone `S(E)` would still pass. Add global strict-monotonicity asserts. Note `E` is declared only `real=True` (line 37); the sign needs the gap `E - Vpeak > 0`.

**Required change:**
1. In section 2 (after line 57, where `dtcross_dE` is computed), add a strict-monotonicity assertion for the crossing time. Introduce a positive gap symbol and rewrite the derivative in terms of it, e.g.:

   ```python
   x = sp.symbols("x", positive=True)   # x = E - Vpeak > 0
   dtcross_dE_gap = sp.simplify(dtcross_dE.subs(E, Vpeak + x))
   # dt_cross/dE = -sqrt(2)*lambda_eff*sqrt(m_s)/(4*x**(3/2)) < 0 for x>0
   assert sp.simplify(dtcross_dE_gap) == sp.simplify(-sp.sqrt(2) * lam_eff * sp.sqrt(m_s) / (4 * x**sp.Rational(3, 2)))
   assert (dtcross_dE_gap < 0) == True   # strict, all symbols positive
   ```

   (If `(expr < 0) == True` does not reduce under SymPy's positivity for this form, instead assert the sign of the numerator coefficient: assert `sp.simplify(-dtcross_dE_gap * 4 * x**sp.Rational(3,2) / (sp.sqrt(2)*lam_eff*sp.sqrt(m_s))) == 1`, which forces the leading sign to be negative for all positive `x`. Choose whichever SymPy reduces cleanly; both establish global strict negativity.)

2. In section 4 (after line 89, where `dS_dE` is computed), add the analogous global strict-monotonicity assertion for the survival ratio:

   ```python
   dS_dE_gap = sp.simplify(dS_dE.subs(E, Vpeak + x))
   # dS/dE = -sqrt(2)*sqrt(chi_peak)*sqrt(g_UV)*lambda_eff*sqrt(m_s)/(4*sqrt(mu_eta)*x**(3/2)) < 0
   assert (dS_dE_gap < 0) == True
   ```

   (Same fallback as above: if `< 0` does not reduce, assert the leading coefficient sign is negative by dividing out the positive `x**(3/2)` and positive parameter factors and checking the residual equals `-1`.) Add a `print("dS/dE (x=E-Vpeak)     =", dS_dE_gap)` line so the output records the gap form.

Do NOT remove or weaken any existing assert. Do NOT touch the numeric benchmark literals (they are paper-anchored). Reuse the existing `lam_eff`, `m_s`, `gUV`, `chi_peak`, `mu_eta`, `Vpeak` symbols already declared.

**Self-test (variable independence / strictness, already checked by auditor):** Both new asserts reference `x` (= `E − Vpeak`), so they are genuine sign tests, not vacuous identically-zero derivatives. With every symbol positive and `x>0`, `dtcross_dE_gap = −(positive)/x^{3/2} < 0` and `dS_dE_gap = −(positive)/x^{3/2} < 0` — both strictly negative and satisfiable, establishing the global monotonicity that makes the window one-sided with a unique edge.

**Verification command:**
The verifier will run `redteam exec-sympy 250` and confirm: the script exits 0; new strict-monotonicity asserts for `dt_cross/dE` (section 2) and `dS/dE` (section 4) are present and reference `x`; no existing assert or numeric literal was changed.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_sympy_audit.py`
- summary: Added positive-gap derivative forms and strict global negativity assertions for `dt_cross/dE` and `dS/dE`.
- deviation: none
