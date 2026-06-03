---
unit_id: 237
batch: VII.2
created_at: 2026-06-02T00:00:00-06:00
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-06-02T22:24:56-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 237

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond what each finding names. Do NOT touch paper.tex, notes/, or any prose document.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

---

## F1 — missing_verification_script (missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_mathematica_audit.wl`

(`.wl` files live in `mathematica/`; saved output goes in `mathematica/output/` with the matching `.txt` name.)

**Issue:** Unit 237 is a non-checkpoint, non-status-only stage with only a SymPy script. The stage is entirely symbolic (log identities, two-term Taylor series in one parameter, a 2x2 determinant, and partial derivatives) and is therefore independently verifiable in Mathematica using native primitives. Sibling Part VII stages 239 and 242 already carry `_mathematica_audit.wl` files, confirming this class of algebra is Mathematica-verifiable. The dual-engine rule requires a second engine wherever it is POSSIBLE.

**Required change:** Write a NEW, independent-route Mathematica audit script verifying the claims M1–M7 below. Independent route means: derive each identity from the physical premises using a DIFFERENT decomposition than the `.py` (use `FullSimplify`, `Series`/`Normal`, `D`, `Det`, `LogExpand`/`PowerExpand` with explicit positivity assumptions via `Assuming[... > 0, ...]`), NOT a line-by-line transliteration of the SymPy variable choreography. Use a residual-zero pattern (e.g. a helper `expectZero[expr_]` that `FullSimplify`s under positivity assumptions and `Exit[1]`s on nonzero, after stripping any `ConditionalExpression[0, ...]` wrapper from Solve/Reduce output). For pole/invertibility checks prefer `1/x == 0` style guards over `Limit`.

**Claim manifest** (each must be an independent, non-tautological check; positivity assumptions: all `R*`, `eps_eta`, `eps_eta_ref` in `(0,1)` where they are dressings, all `K_*`, `c_etaU*`, `lambda_*` positive):

- **M1 — rigid-mouth reduction.** With `q_tr = -C* ln(Rtr/Rtr_ref)`, `q_nt = B* ln(Rtr/Rtr_ref) + ln((1-eps_eta)/(1-eps_eta_ref)) - ln(Rtarget/Rtarget_ref)`, `q_eta = ln(eps_eta/eps_eta_ref)`: on the rigid-mouth slice `Rtr = Rtr_ref`, show `q_tr = 0` and `q_nt = ln((1-eps_eta)/(1-eps_eta_ref)) - ln(Rtarget/Rtarget_ref)`.
- **M2 — finite static-blind curve and inverse.** Show that `q_nt = 0` (rigid-mouth form, with `Rratio := Rtarget/Rtarget_ref`) is equivalent to `Rratio = (1-eps_eta)/(1-eps_eta_ref)`; and that the parameterization `eps_eta = eps_eta_ref e^{q_eta}` gives `Rratio = (1 - eps_eta_ref e^{q_eta})/(1 - eps_eta_ref)`; and that the inverse `q_eta = ln((1 - (1-eps_eta_ref) Rratio)/eps_eta_ref)` round-trips both ways (composition reduces to the identity on each variable). Endpoint: `Rratio(q_eta=0) = 1`.
- **M3 — first-order packet compiler and tangent slope.** Linearizing `q_eta` and `q_nt` (rigid-mouth) to first order in a single parameter `t` along `eps_eta = eps_eta_ref e^{t E1}`, `Rratio = e^{t R1}`, show `q_eta = t E1` and `q_nt = -t (R1 + c_eta E1)` with `c_eta = eps_eta_ref/(1-eps_eta_ref)`. Independently show the tangent slope `d/dq_eta [ln Rratio(q_eta)] |_{q_eta=0} = -c_eta`. Show the packet matrix `M = {{-1, -c_eta}, {0, 1}}` has `Det[M] = -1`.
- **M4 — microscopic coherent compiler.** With `eps_eta = c_etaU^2/(K_U K_eta_eff)`, show `q_eta = ln((c_etaU^2/(K_U K_eta_eff)) (K_U_ref K_eta_eff_ref / c_etaU_ref^2))` equals `2 ln(c_etaU/c_etaU_ref) - ln(K_U/K_U_ref) - ln(K_eta_eff/K_eta_eff_ref)` (log-expansion under positivity).
- **M5 — microscopic first-order drift extractor.** Along `c_etaU = c_etaU_ref e^{t c1}`, `K_U = K_U_ref e^{t kappa_U}`, `K_eta_eff = K_eta_eff_ref e^{t kappa_eta}`, show the first-order term of `q_eta` in `t` equals `2 c1 - kappa_U - kappa_eta`.
- **M6 — support-blindness (NON-VACUOUS; see F2).** Do NOT differentiate a `q_eta` expression with respect to symbols it does not contain. Instead build `q_eta` in a frame where the support variables `{zeta, M_tr, lambda_phi, K_phi_eff}` ARE in scope, then show the derivative vanishes BECAUSE the dressing-sector quantities `{c_etaU, K_U, K_eta_eff}` are functionally independent of the support sector. Acceptance: the differentiated expression's variable list must (under the construction) plausibly include the support variable being differentiated, and the assert must be one that WOULD fail if a support variable leaked into `q_eta`. Show `D[q_eta, zeta] = 0`, `D[q_eta, M_tr] = 0`, `D[q_eta, lambda_phi] = 0`, `D[q_eta, K_phi_eff] = 0`, where `zeta = lambda_phi^2 K_W_eff/(lambda_W^2 K_phi_eff)` and `M_tr = M_mix (1 + zeta(1-eps)/(1-zeta eps))`.
- **M7 — dependent-plane ray and codimension-two orbit-lock.** Show `y_eta = -q_eta (0,1,1)^T` equals `-ln(eps_eta/eps_eta_ref)(0,1,1)^T`, and equals `-ln((1-(1-eps_eta_ref) Rratio)/eps_eta_ref)(0,1,1)^T` on the static-blind curve. Show the finite orbit-lock endpoint: at `Rratio = 1` (equivalently `eps_eta = eps_eta_ref`), `q_eta = 0`; and the first-order codim-two point `q_nt = 0 & q_eta = 0 <=> R1 = 0 & E1 = 0`.

**Verification command:** The verifier runs `redteam exec-mathematica 237` and confirms the new `.wl` contains the M1–M7 checks and exits 0, and that it is an independent decomposition (not a transliteration).

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_mathematica_audit.wl`
- summary: Added an independent Mathematica audit covering M1-M7 with residual-zero checks and support-sector independence rules.
- deviation: none

---

## F2 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_sympy_audit.py:227-241`

**Issue:** The support-blindness block (lines 227–241) differentiates `q_eta_micro` (lines 193–196, depending only on `{c_etaU, K_U, K_eta_eff}` and refs) with respect to `zeta`, `M_tr`, `lambda_phi`, `K_phi_eff` — symbols that are NEVER substituted into `q_eta_micro`. Each `sp.diff(q_eta_micro, VAR)` is identically zero simply because `VAR` is not a free symbol of `q_eta_micro`, for ANY expression. The headline support-blindness theorem is therefore verified vacuously: the assertion cannot fail no matter what `q_eta` is. `zeta_expr` and `M_tr_expr` (lines 235–236) are computed but only printed (lines 244–245), so they are dead with respect to the asserts.

**Required change:** Restructure the §5 block so the four `assert_zero(sp.diff(q_eta, VAR))` checks are non-vacuous — i.e. the support variables are in scope of the differentiated `q_eta` expression so that a counterfactual leak would make an assert fail. Codex designs the exact construction; the requirement is:

1. Keep the four deliverables: `∂_zeta q_eta = 0`, `∂_{M_tr} q_eta = 0`, `∂_{lambda_phi} q_eta = 0`, `∂_{K_phi_eff} q_eta = 0`, with `zeta = lambda_phi**2 K_W_eff/(lambda_W**2 K_phi_eff)` and `M_tr = M_mix (1 + zeta(1-eps)/(1-zeta eps))` (lines 235–236, already correct).
2. Build `q_eta` for this block in a form where the support variables COULD appear. The physically faithful construction: the dressing-sector quantities `c_etaU`, `K_U`, `K_eta_eff` are, by the stage's premise, functions of the dressing/microscopic sector ONLY and carry no dependence on the support variables `{lambda_phi, K_phi_eff, lambda_W, K_W_eff, eps, M_mix}`. Express `q_eta` as `sp.log(eps_eta/eps_eta_ref)` where `eps_eta` is a generic symbol, AND separately assert that the support-sector composites `zeta_expr`, `M_tr_expr` do not enter `eps_eta` — i.e. assert `eps_eta` (and hence `q_eta`) is free of `{zeta, M_tr, lambda_phi, K_phi_eff}` as an EXPLICIT precondition, then the vanishing derivative is the THEOREM, not a coincidence of symbol absence. The acceptable mechanical form: (a) add an explicit assertion that the set `{zeta, M_tr, lambda_phi, K_phi_eff, lambda_W, K_W_eff, eps, M_mix}` is disjoint from `q_eta_micro.free_symbols` AND from the support-sector definitions of the dressing quantities (documenting the sector split that makes the theorem true), and (b) keep the four `assert_zero(diff(...))`. Alternatively (preferred, more substantive): define `q_eta` via the chain through `c_etaU(c_etaU_indep)`, `K_U(...)`, `K_eta_eff(...)` where each dressing quantity is `Function`-wrapped or built so its `free_symbols` are dressing-only, then differentiate the resulting `q_eta` w.r.t. the support variables and the composite `zeta`/`M_tr` (using `sp.diff(q_eta.subs(M_tr, M_tr_expr).subs(zeta, zeta_expr), lambda_phi)` etc. so the chain rule actually has a path to traverse), and show it is zero. Choose whichever yields an assert that would FAIL if a support symbol were introduced into the dressing definitions.

3. Whichever construction you choose, ADD a deliberate negative-control comment (no new assert needed) noting that if `K_eta_eff` were redefined to depend on `K_phi_eff`, `∂_{K_phi_eff} q_eta` would be nonzero — this documents non-vacuity for the reviewer.

**Self-test the fix:** before finalizing, confirm with `q_eta.free_symbols` (mentally) that the differentiated expression in each of the four asserts has the support variable in a position where the chain rule could produce a nonzero term, so the zero result is the theorem. Do NOT replace the four physical-variable derivatives with anything weaker.

**Verification:** The four `assert_zero` calls for `zeta`, `M_tr`, `lambda_phi`, `K_phi_eff` remain; the differentiated `q_eta` is constructed so the support variables are reachable by the chain rule; script exits 0; output lines for `dq_eta/d{zeta,M_tr,lambda_phi,K_phi_eff}` still print `0`.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_sympy_audit.py`
- summary: Rebuilt the support-blindness block in a support-frame function chain and imposed dressing/support sector independence before asserting the four derivatives vanish.
- deviation: none

---

## F3 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_sympy_audit.py:170-174`

**Issue:** Line 173 defines `q_lin_vec = sp.simplify(M_packet * sp.Matrix([R1, E1]))` and line 174 asserts it equals `sp.Matrix([-R1 - c_eta * E1, E1])` — the hand-expansion of that same matrix product. The assert is guaranteed by construction and cannot fail; it does not anchor `M_packet` to any independently-derived quantity.

**Required change:** Make line 174 compare the packet action `M_packet * Matrix([R1, E1])` against the INDEPENDENTLY series-derived linearization `(q_nt_linear, q_eta_linear)` from lines 156–162, rather than against a literal restatement of the product. Both `q_nt_linear` and `q_eta_linear` carry a single factor of `t`; reconcile the `t` factor consistently (e.g. compare `t * (M_packet * Matrix([R1, E1]))` against `Matrix([q_nt_linear, q_eta_linear])`, or strip `t` from the series forms via `.subs(t, 1)` and compare against `M_packet * Matrix([R1, E1])`). Example acceptable form:
```python
assert_matrix_zero(
    M_packet * sp.Matrix([R1, E1])
    - sp.Matrix([q_nt_linear.subs(t, 1), q_eta_linear.subs(t, 1)])
)
```
(Verify the sign convention matches: `q_nt_linear` at line 165 satisfies `q_nt_linear = -t(R1 + c_eta E1)`, and `M_packet*[R1,E1]` first component is `-R1 - c_eta E1`, so with `t->1` these agree; the second component `E1` matches `q_eta_linear.subs(t,1) = E1`.) Remove or repurpose line 173's `q_lin_vec` if it is no longer needed.

**Verification:** Line 174 (or replacement) compares the packet action to the series-derived `(q_nt_linear, q_eta_linear)`, NOT to a restatement of the product; a wrong entry in `M_packet` would now make the assert fail. Script exits 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_sympy_audit.py`
- summary: Compared the packet matrix action against the independently series-derived first-order coordinates with the series parameter stripped.
- deviation: none

---

## F4 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_sympy_audit.py:265`

**Issue:** Line 265 `assert_zero(sp.simplify(qeta_param.subs(qeta_param, 0)))` substitutes `0` for the symbol `qeta_param` and asserts the result is `0` — i.e. `assert 0 == 0`. It exercises no physics and cannot fail. It is meant to be the `q_eta = 0` half of the orbit-lock endpoint but checks nothing about `q_eta`.

**Required change:** Replace line 265 with a substantive check that the actual dressing coordinate vanishes at the reference: assert that `q_eta` (line 95, `sp.log(eps_eta/eps_eta_ref)`) is `0` when `eps_eta = eps_eta_ref`:
```python
assert_zero(q_eta.subs(eps_eta, eps_eta_ref))
```
This makes the `q_eta = 0` endpoint a real check on the `q_eta` expression (it reduces to `log(1) = 0`, true and non-trivial — it would fail if `q_eta` were mis-defined).

**Verification:** Line 265 substitutes `eps_eta -> eps_eta_ref` into the `q_eta`-bearing expression and asserts zero; a mis-defined `q_eta` would fail. Script exits 0.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_sympy_audit.py`
- summary: Replaced the symbolic self-substitution with a check that the actual `q_eta` expression vanishes at `eps_eta = eps_eta_ref`.
- deviation: none
