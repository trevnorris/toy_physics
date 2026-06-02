---
unit_id: 202
batch: VI.1
created_at: 2026-06-01T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-02T09:52:31-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 202

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes beyond what F1 requires.

After editing, RUN the new script (`math -script <path>`) and iterate until it exits 0 with all in-file checks passing. Getting the script to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose document, and do NOT modify the existing SymPy script. The red-team only adds the missing second-engine script here.

## F1 — missing_verification_script (subtype: missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage202_free_quintuple_target_graph_mathematica_audit.wl` (new file — `.wl` files live in `mathematica/`)

**Issue:** Stage 202 (non-status-only, non-checkpoint) is verified by SymPy only. It performs substantive symbolic algebra that Mathematica can independently verify, so the project dual-engine contract requires a second-engine `.wl`. You must author a NEW Mathematica script that verifies the stage's physical results independently. This directive states the requirement and acceptance criteria only — you design and write the implementation. The SymPy script is at `scripts/moving_throat_pde_stage202_free_quintuple_target_graph_sympy_audit.py` and the paper formulas are in `paper/appendices/stage_appendix_part06.tex` (§"Direct microscopic target graph", eqs `app-part06-Ctr-direct`, `Cnt-direct`, `epsilon-eta-direct`, `deltaU-graph`, `TU-graph`, `Keta-graph`, `muW-graph`, `main-graph-quotient`) and the notes file §§1–7. Use those as the source of truth for the symbolic forms; do not read the SymPy algebra as a template.

**Symbol conventions (match the paper):** positive reals `lamW, cetaU, gamma, KU, Keta, KW, muW, TU, L, sigma`; positive constants `chi0s ≡ 1+\chi_{0,*}`-related exponent `\chi_{0,*}`, `deltaUs ≡ \delta_{U,*}`, `Estar ≡ E_*`, `Fstar ≡ F_*`; positive target data `Ctrtgt, Cnttgt, epsEtatgt`. The monomials are
- `Ctr = (gamma*cetaU/KU)^(1+deltaUs) * (Pi^2 TU/(L^2 KU))^(1+chi0s)`
- `Cnt = (lamW^2 muW/(Keta KW^2)) * ((gamma^2 lamW^2 sigma)/(KU KW))^Estar * (Pi^2 TU/(L^2 KU))^(-Fstar)`
- `epsEta = cetaU^2/(KU Keta)`

**Claim manifest** (each must be independently verified by the new `.wl`):

- **M1** — Target-graph solve reconstructs the tracking monomial. Solve `Ctr = Ctrtgt` for the split-`U` ratio `delta_U := Pi^2 TU/(L^2 KU)` (treat the exponent `1+deltaUs` as an independent constant), define `deltaU_graph` from it, set `T_graph = L^2 KU deltaU_graph/Pi^2`, and verify `PowerExpand[Log[ (Ctr /. TU -> T_graph) / Ctrtgt ]] // FullSimplify == 0` under positivity.
- **M2** — Dressing solve reconstructs `epsEta`. Define `Keta_graph = cetaU^2/(KU epsEtatgt)` and verify `PowerExpand[Log[ (epsEta /. Keta -> Keta_graph) / epsEtatgt ]] // FullSimplify == 0`.
- **M3** — Nontracking solve reconstructs `Cnt`. Define `mu_graph = Cnttgt cetaU^2 KW^2/(epsEtatgt KU lamW^2) * ((gamma^2 lamW^2 sigma)/(KU KW))^(-Estar) * deltaU_graph^Fstar`, and verify `PowerExpand[Log[ (Cnt /. {TU -> T_graph, Keta -> Keta_graph, muW -> mu_graph}) / Cnttgt ]] // FullSimplify == 0`.
- **M4** — `mu_graph` is independent of `L` and `Pi`. Verify `FreeQ[FullSimplify[PowerExpand[mu_graph]], L] && FreeQ[FullSimplify[PowerExpand[mu_graph]], Pi]` is `True` (equivalently `D[Log[mu_graph], L] // PowerExpand // FullSimplify == 0`).
- **M5** — Graph-error identities equal the quotient packet. With `qtr = Log[Ctr/Ctrtgt]`, `qnt = Log[Cnt/Cnttgt]`, `qeta = Log[epsEta/epsEtatgt]`, and `E_T = Log[TU/T_graph]`, `E_K = Log[Keta/Keta_graph]`, `E_mu = Log[muW/mu_graph]`, verify under positivity (use `PowerExpand`+`FullSimplify`):
  - `E_T - qtr/(1+chi0s) == 0`
  - `E_K + qeta == 0`
  - `E_mu - (qnt - qeta + Fstar*qtr/(1+chi0s)) == 0`
- **M6** — Graph projection equals the Stage-201 canonical orbit projection. Build the canonical projection 8-vector independently from the monomial ratios `Rtr = Ctr/Ctrtgt`, `Rnt = Cnt/Cnttgt`, `Reta = epsEta/epsEtatgt`:
  `x_proj_can = {lamW, cetaU, gamma, KU, Keta*Reta, KW, muW*Rnt^(-1)*Reta*Rtr^(-Fstar/(1+chi0s)), TU*Rtr^(-1/(1+chi0s))}`,
  and the graph vector `x_graph = {lamW, cetaU, gamma, KU, Keta_graph, KW, mu_graph, T_graph}`, then verify componentwise `PowerExpand[Log[x_proj_can[[i]]/x_graph[[i]]]] // FullSimplify == 0` for all 8 entries.
- **M7** — Reduced-family packet and repair vector. Verify the repair-vector rewrite `{0,0,0,0,qeta,0, -qnt+qeta-Fstar*qtr/(1+chi0s), -qtr/(1+chi0s)} - {0,0,0,0,-E_K,0,-E_mu,-E_T} == 0` componentwise, and verify the reduced-family packet `{chiQ-1, E_T, E_K, E_mu}` evaluated at `{TU->T_graph, Keta->Keta_graph, muW->mu_graph, chiQ->1}` is the zero 4-vector (use `PowerExpand`+`FullSimplify`).

**Anti-transliteration guard:** The `.wl` must derive these results independently using native Mathematica primitives — `Solve`/`Reduce` to obtain the graph triple from the monomial equations (rather than pasting the closed-form `deltaU_graph = (Ctrtgt/(...))^(1/(1+chi0s))` straight across), `PowerExpand`/`FullSimplify`/`Assuming` under positivity for the log-ratio reductions, `D[]` or `FreeQ` for the `L`-independence check, and matrix/`Thread` ops for the vector comparisons. A line-by-line port of the `.py` variable choreography (same hand-substitution order, same intermediate names rewritten in Mathematica syntax) is rejected as `mathematica_transliteration`. Use a DIFFERENT decomposition than the `.py`: e.g. obtain `deltaU_graph` and `T_graph` from `Solve[Ctr == Ctrtgt, TU]` (then read off the split-`U` ratio) instead of writing the explicit fractional-power closed form by hand.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 202` and confirm the new `.wl` appears at the target path AND the script exits 0 with each of M1–M7 reporting a zero residual (or `True`).

## Iteration 2 (orchestrator) — TIMEOUT REWORK (your iter-1 .wl hangs)

The iter-1 `.wl` you wrote is CORRECT in intent but TOO SLOW: the orchestrator's independent re-run TIMED OUT at the 600s cap (exit 124). It hung in **Section I, at the first symbolic `Solve`** — `Solve` on the transcendental/fractional-power log-monomial equation (with symbolic exponents `1+deltaUs`, `1+chi0s`) does not terminate in reasonable time. Per project policy a timeout is a FAILURE → reformulate the math; NEVER rely on the run "just barely finishing," and never raise the cap.

**Requirement (iter 2):** reformulate so the ENTIRE script finishes well under 600s (target < ~120s) while still INDEPENDENTLY establishing M1–M7 (keep the claim manifest). Do NOT symbolic-`Solve` the transcendental log equations. You choose the route; acceptable strategies (any one, or a mix):

- **Linearize in log-space (recommended, fast + independent).** Introduce log-variables `LT=Log[TU]`, `LK=Log[Keta]`, `LMu=Log[muW]`, and the log-data `Log[Ctrtgt]`, etc. The monomial-match equations `Log[Ctr/Ctrtgt]=0`, `Log[epsEta/epsEtatgt]=0`, `Log[Cnt/Cnttgt]=0` become **LINEAR** in `{LT,LK,LMu}` (the exponents become coefficients). `Solve`/`LinearSolve` that linear system — instant — then exponentiate to recover `T_graph, Keta_graph, mu_graph`. This is genuinely different from the `.py`'s hand closed-form and from the slow nonlinear Solve.
- **OR high-precision numeric verification.** Pick several concrete positive-rational assignments of ALL free parameters and assert each M-item log-residual is `< 10^-30` via `N[..., 40]` at each sample (mpmath-style). Sidesteps symbolic Solve entirely; still an independent check.
- For the M-item reductions, replace any global `FullSimplify` on the full symbolic monomials that is slow with `PowerExpand` + targeted `Simplify`, or evaluate the residual at concrete positive rationals. Wrap any remaining symbolic `Simplify`/`FullSimplify` in `TimeConstrained[expr, 30, <numeric fallback>]` so a single slow reduction cannot stall the whole script.

**Acceptance:** `redteam exec-mathematica 202` exits 0 in well under 600s with all M1–M7 PASS lines present; the `.wl` remains an INDEPENDENT derivation (the linear-log-solve / numeric route is a different decomposition than the SymPy `.py`, not a transliteration). Update your `## Applied: F1` block to note the iter-2 reformulation (and add the block if it is missing).

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage202_free_quintuple_target_graph_mathematica_audit.wl`
- summary: Added a second-engine Mathematica audit for M1-M7, reformulated in iter 2 as a fast linear log-system solve for the dependent graph variables with exact residual checks.
- deviation: none
