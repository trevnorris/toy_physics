---
unit_id: 243
batch: VIII.1
created_at: 2026-06-03T09:30:47-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-03T09:48:25-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 243

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named. Do NOT touch paper.tex, notes/, or any prose document. The SymPy script is correct and fully paper-aligned — do NOT modify it.

After editing, RUN the affected script (`math -script <path>`) and iterate until it exits 0 with all in-file checks passing. Getting the script to run cleanly is your job; the orchestrator independently re-runs afterward.

## F1 — mathematica_transliteration (+ F2 banner, folded in)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_mathematica_audit.wl:59-214`

**Issue:**
The `.wl` is a line-by-line transliteration of the SymPy script (`scripts/moving_throat_pde_stage243_..._sympy_audit.py`): same block order, same variable choreography, same intermediate decomposition, and — critically — the SAME hardcoded `expected*` closed forms checked as residuals (`expectedSleak = -Sqrt[2] ellW j0/4` at wl:90, `Uexpected = kV fU/(kU kV - chiLam^2)` at wl:106, `drainExpected` at wl:115-118, etc.). Because both engines decompose the problem identically AND compare against the same hardcoded constants, the second engine provides no independent cross-check: a transcription error in a hardcoded form would not be caught. On a checkpoint stage this is a high-severity defect. Additionally the banner at wl:59 reads `STAGE 226` (a clone artifact) and must read `STAGE 243`.

Re-author the `.wl` so each block confirms its result by an INDEPENDENT mechanism — derive/verify from the defining premises by a different route than the SymPy script, and do NOT keep the same `expected*` literals as the sole reference forms. Keep all existing `expectZero`/`pass`/`fail` infrastructure and the `$Assumptions` block; replace only the verification *logic* of each block as specified by the claim manifest below. Keep using `math`-native primitives (Integrate, Solve/LinearSolve, Series/Limit, integration-by-parts via `D[W jw, w]`).

**Required change (step by step):**

1. wl:59 — change banner string `"STAGE 226 — RELAXED-CONSTRAINT BRANCH DECLARATION AND SHORT-RANGE OPEN-SYSTEM COMPILER"` to `"STAGE 243 — RELAXED-CONSTRAINT BRANCH DECLARATION AND SHORT-RANGE OPEN-SYSTEM COMPILER"`.

2. Block I (leakage/work, wl:73-97) — keep `W, jw, Ew` definitions and the direct `Integrate` computation of `Sleak`/`Wwork`. REMOVE reliance on the hardcoded `expectedSleak`/`expectedWwork` as the only reference: instead confirm `Sleak` by the integration-by-parts identity it is supposed to satisfy — i.e. verify `expectZero["IBP closure", Integrate[D[W jw, w], {w, -Infinity, Infinity}] - (boundary + Sleak)]` (the antiderivative-of-the-total-derivative equals boundary + interior, independent of any closed form). You MAY additionally keep the closed-form residual check, but the IBP closure must be present as the independent confirmation. For `Wwork`, confirm independently by a substitution route, e.g. evaluate via `Integrate[jw Ew /. w -> t, {t, -Infinity, Infinity}]` after rescaling, or confirm the value by `2 Integrate[jw Ew, {w, 0, Infinity}]` (exploiting even-integrand symmetry) and compare to the half-line result — a different decomposition than the SymPy full-line integral.

3. Block II (non-rigid U/V, wl:99-128) — REMOVE the hardcoded `Uexpected`/`Vexpected`/`drainExpected` as the reference forms. Confirm the solution two independent ways: (a) back-substitute into the stationarity system — `expectZero["U/V satisfy stationarity", stationarity /. uvSol]` (must give `{0,0}`); and (b) solve the SAME 2x2 linear system with a native matrix primitive `LinearSolve[{{kU, -chiLam}, {-chiLam, kV}}, {fU, 0}]` and check that its components equal `{Usol, Vsol}` via `expectZero`. The drain, det-Hessian, and `V/U` checks may stay as residuals against the symbolic `Usol/Vsol` (they are derived in-script, not hardcoded), but the drain's reference `drainExpected` literal should be replaced by `chiLam Usol Vsol` recomputed and checked to be `>= 0`-form `chiLam^2 kV fU^2 / detH^2` using the in-script `detH` rather than a re-typed denominator.

4. Block III (compensated source, wl:130-157) — keep `Integrate[varsigma, {z,0,1}]` for the mean (this is already a genuine native computation). For the quadratic rewrite, REMOVE the manual substitution rule `{Cos[Pi z]->y, Cos[2 Pi z]->2 y^2-1}` (which mirrors the SymPy `.subs`) and instead confirm the rewrite independently: verify that `varsigma` and the candidate quadratic `1 - b + a y + 2 b y^2` agree as functions of `z` under `y = Cos[Pi z]` by checking `expectZero["quadratic rewrite (functional)", (varsigma - ((1 - b + a y + 2 b y^2) /. y -> Cos[Pi z]))]` (uses `TrigExpand`/`FullSimplify` on the cosine double-angle, a different route than hand-substituting `2 y^2 - 1`). For the vertex, confirm `y_*` is the stationary point by `expectZero["vertex is stationary", D[1 - b + a y + 2 b y^2, y] /. y -> yStar]` (derivative vanishes) in ADDITION to the value check.

5. Block IV (recovery slice, wl:159-174) — keep as is; it is a genuine substitution on in-script-derived objects (`Usol`, `Vsol`, `drainUV`, `Sleak`, `Wwork`, `varsigma`), not a transliteration of a hardcoded form. No change required beyond the upstream `Usol/Sleak/...` now being independently confirmed.

6. Block V (short-range limits, wl:176-214) — keep `SQ, SY, QQ, QY, YY` and the `deltaVStat`/`VdynReal` definitions. For the limits, confirm the no-new-long-range result by an independent asymptotic route rather than only `Limit`: use `Series[x QQ, {x, Infinity, 1}]` (and the QY, YY analogues) and check the leading term -> 0, i.e. `expectZero["x QQ -> 0 (series)", Normal[Series[x QQ, {x, Infinity, 0}]]]` style, OR confirm `lim x*mode == 0` via the explicit power/exponential bound (e.g. `Limit` plus a `Series` leading-order check). The point is a second mechanism, not a re-run of the same `Limit` the SymPy script uses.

**Claim manifest** (the independent `.wl` route must verify each):
- M1 — leakage: `S_leak = integral_{-inf}^{inf} W'(w) j^w(w) dw = -(sqrt2/4) ell_w j0`, confirmed independently via IBP closure `int D[W jw, w] dw = [W jw] + S_leak`; boundary term `= 0`. Work: `W_w = int j^w E_w dw = (sqrt(2 pi)/8) ell_w j0 E0`, confirmed via even-integrand half-line decomposition.
- M2 — non-rigid: `(U, V)` satisfy `{k_U U - chi_lam V - f_U, k_V V - chi_lam U} = {0,0}` (back-substitution) AND equal `LinearSolve[{{k_U,-chi_lam},{-chi_lam,k_V}},{f_U,0}]`. `det H = k_U k_V - chi_lam^2`; `V/U = chi_lam/k_V`; `D_UV = chi_lam^2 k_V f_U^2 / (k_U k_V - chi_lam^2)^2`.
- M3 — compensated source: `int_0^1 varsigma dz = 1`; `varsigma == (1 - b + a y + 2 b y^2)|_{y=cos(pi z)}` as functions of z; `y_* = -a/(4b)` is the stationary point (derivative vanishes) with value `1 - b - a^2/(8b)`; boundary values `1 +/- a + b`.
- M4 — recovery slice: `(ell_w, f_U, a, b) = (0,0,0,0)` => `S_leak = W_w = U = V = D_UV = 0`, `varsigma == 1`.
- M5 — short-range span/limits: `QQ = x^-6`, `QY = e^{-2 kappa x}/x^4`, `YY = e^{-4 kappa x}/x^2`; and `lim_{x->inf} x*QQ = x*QY = x*YY = x*deltaV_stat = x*Re V_dyn = 0`, confirmed by a series/asymptotic route in addition to `Limit`.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 243` and confirms: (i) the literals `expectedSleak`, `Uexpected`, `Vexpected`, `drainExpected` no longer serve as the sole reference forms (independent IBP / LinearSolve / back-substitution / series checks are present); (ii) the banner reads `STAGE 243`; (iii) the script exits 0 with every `expectZero` reporting PASS.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_mathematica_audit.wl`
- summary: Replaced the hardcoded-reference Mathematica checks with IBP, half-line, stationarity back-substitution, LinearSolve, functional trig-rewrite, stationary-vertex, and asymptotic-series witnesses.
- deviation: The IBP residual uses the mathematically complete closure `Sleak + Integrate[W D[jw,w], ...] - boundary` because the literal sketch omitted the interior derivative term.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_mathematica_audit.wl`
- summary: Corrected the folded banner finding from `STAGE 226` to `STAGE 243`.
- deviation: none
