# Directive pathA_14 — Chunk B2c: integrate the derived bundle → R_norm(τ)=0 calibration + §J error bars

**Date:** 2026-06-18
**Owner:** Codex (codes/designs). Claude reviews only.
**Builds on:** decision-12 §B2c (the architecture); B1 `patha_extraction.py` (committed `398ba27`); B2a
`patha_b2a_bdg.py` (committed `fa5712e`); B2b `patha_b2b_maxwell.py` (committed `c0c78e1`); GATE-A freeze
`ed358569…b1691c9`. Family = `homogeneous_isotropic_hooke_v1` (1 calibrated DOF `{τ}`).

This is the **decisive verdict chunk.** Calibrate the single wall DOF `τ` on the one anchor `R_norm(τ)=0`, on the
stable side `D0>0`, and REPORT the held-out predictions `R_pole/P2/P4` + naturalness. Build **target-blind in the
held-out sense:** `R_norm` is the calibration target (allowed); you may NOT tune anything to the held-out
`R_pole/P2/P4` targets.

---

## 0. Prep finding you must supersede, not assume (the whole reason B2c re-solves)
At `τ=1` the assembled bundle gives `R_norm = mhat0²·S_port·P0 − 54G c_s⁵/(5a⁵c⁵) = −10.800` exactly, because
`P0 = N0/D0 = 2.158e-8 / 7.7209 = 2.79e-9 ≪ 10.8`. Under a **frozen-background** sketch (hold `B0,Z0,N0` at their
τ=1 values, scale only `K=τκ̂`), `R_norm=0` would land at `τ* ≈ 5.06e-6` (`|ln τ*| ≈ 12.2`) sitting right on the
`D0→0` knife-edge (`D0 ≈ 2e-9`, cancellation ratio `K/(B0+Z0) → 1`). **This frozen-background number is a PRIOR, not
a result.** The entire verdict hinges on whether `{B0,Z0,N0}` (and `R0`) shift as the wall softens at small `τ`. B2c
MUST answer that by **real re-solves**, never by the frozen-background shortcut. A frozen-background root-find would
beg the question and is FORBIDDEN as the production path (it is allowed ONLY as an explicit negative-control contrast,
see §R6).

---

## 1. R-requirements (acceptance criteria in §A)

### R1 — Full per-τ re-solve of the derived bundle (NO frozen-background shortcut)
Every `R_norm(τ)` you evaluate must come from a genuine re-derivation at that `τ`:
1. Re-solve the **closed background** at `τ` on the frozen `homogeneous_isotropic_hooke_v1` family
   (`solve_closed_background(tau=…, grid_level=…)` / `make_background_bundle`), to the same convergence the B2a/B2b
   runs used (closed residual ~`2.7e-10` class at τ=1; report the achieved residual at each τ).
2. Re-derive the **BdG bundle** `{B0,B2,B4}` on that background (the B2a engine, **converged modal sweep K=30 of
   the positive spectrum** — not the truncated 3-mode sum that B2a's adversarial pass rejected).
3. Re-derive the **Maxwell bundle** `{Z0,Z2,Z4,N0,N2,N4}` on that background's `A0` (the B2b basis-invariant
   Green/self-energy transfer; production grid `47×17`, window `0.028`, radial_scale `5.0`, or finer).
4. Get `K=τκ̂` (`κ̂=7.7209…` analytic, B1), `M`, and the wall mode `χ` from B1 on the same background.
5. Assemble `direct_coefficients={K,M,B0,B2,B4,Z0,Z2,Z4,N0,N2,N4}` → `lane_extract` → `D0=K−B0−Z0`, `P0=N0/D0` →
   `observable_residuals` → `R_norm`.

You do NOT have to re-run the full convergence/dual-engine sweep at every τ in the root-find — that is for τ* only
(§J). Each root-find evaluation needs ONE converged bundle at the production grid. **Re-use the existing B2a/B2b
engines** (already dual-engine + fidelity validated); do not re-implement the BdG or Maxwell operators.

### R2 — Route is yours to design, under the ≤10-min/script cap
Each individual script must finish under `timeout 600` (exit 124 = reformulate the route, never raise the cap —
[[feedback-script-timeout-policy]]). A single monolithic 30-evaluation root-find loop will exceed that. You decide the
structure (examples, not a mandate): (a) full re-solve at a modest set of real τ anchors spanning the bracket, build a
**validated** smooth model `{B,Z,N,K}(τ)`, root-find on the model, then **CONFIRM with one more full re-solve at the
found τ\*** (residual gate, §R3); or (b) a direct bracketed root-find checkpointed across short runs. Whatever route:
the τ anchors must be REAL re-solves (R1), and the τ\* must be CONFIRMED by a real re-solve.

### R3 — Root-find: deterministic, stable side, confirmed
- Bracket on the **stable side `D0>0`** only. `D0` crosses 0 at `τ_crit=(B0+Z0)/κ̂`; the physical branch is
  `τ>τ_crit`. Reject any root with `D0≤0`.
- Deterministic method (bisection/Brent), fixed tolerance, no random seeds.
- **Confirmation gate:** at the reported `τ*`, a fresh full re-solve (R1) must give `|R_norm(τ*)| ≤ tol_Rnorm` AND
  the bundle `{B,Z,N,K}` used to locate `τ*` must agree with the re-solved bundle at `τ*` within the §J error bars.
  If the route used a fitted model (R2a), this gate is what proves the fit didn't lie.

### R4 — Soft-wall honesty / τ_floor (no silent truncation)
The wall is `τ×` softer than τ=1; near `τ_crit≈5e-6` the closed solve may fail to converge. If so:
- Report the lowest `τ` at which the closed solve **genuinely converges** (`τ_floor`), with the achieved residual.
- Do NOT extrapolate the background, BdG, or Maxwell bundle below `τ_floor`. If `τ*` would lie below `τ_floor`,
  REPORT THAT as the finding ("calibration root is below the resolvable wall-stiffness floor"), do not fabricate it.
- Log every dropped/failed τ ([[no silent caps]]); a non-converged solve is a failure, not a data point.

### R5 — §J error propagation into R_norm(τ*) and held-out
Propagate the **recorded** error budgets (do not invent new ones):
- B2a `error_budget`: modal-truncation `1.4e-5`, spatial-confirmation `9.3e-3`, spatial-ladder `2.0e-2`.
- B2b `error_budget`: per-coefficient `Maxwell_ZN_rel` (`Z0/Z2/N0/N2 ~0.4–0.6%`, `Z4 ~2.2%`, `N4 ~1.7%`),
  cross-engine `~4.4%`.
Method is yours (linearized sensitivity `∂R/∂coeff·δcoeff`, or re-evaluate at perturbed coefficients). Output error
bars on: `τ*`, `D0(τ*)`, `R_norm` margin, and the held-out `R_pole/P2/P4`. Note explicitly that the `R_norm`
calibration rides on the tight `N0`(0.4%)+`K`(exact), while the soft `Z4/N4` feed only held-out `P4/R_pole`.

### R6 — Naturalness + negative control (the verdict's substance)
- Report `τ*`, `|ln τ*|`, the cancellation ratio `K/(B0+Z0)` at τ* (+ how many leading digits cancel in
  `D0=K−B0−Z0`), and `D0(τ*)` with its error bar.
- **Negative control (proves the re-solve matters):** also compute the frozen-background root `τ*_frozen` (the §0
  prior) and SHOW how far the real per-τ re-solve moves it. If the background is τ-stable they nearly coincide (a
  reportable finding — then state it with the measured `{B,Z,N}(τ)` drift across the bracket as evidence); if it
  moves, `τ*` differs and you've shown the shortcut would have been wrong. Either way the contrast is REQUIRED
  evidence, not optional.

### R7 — Target-blind firewall (held-out)
- The ONLY knob is `τ`; it is calibrated SOLELY on `R_norm=0`.
- Compute and report `R_pole, P2, P4` at τ* as **predictions**. Do NOT compare them to their structural targets to
  adjust τ, the grid, the family, or any coefficient. No held-out target appears anywhere in the calibration path.
- No new GATE-A DOF, no rescue parameter. If `R_norm=0` cannot be hit on the stable side within `[τ_floor, τ_max]`,
  that is the finding (κ_PV discipline: a miss → DIAGNOSE missing physics in a LATER chunk, never a knob here).

### R8 — Dual-engine on the NEW surface
The new math in B2c is the **assembly + observable + root-find algebra** (the per-τ background/BdG/Maxwell re-solves
already carry their own dual-engine validation from B2a/B2b). Provide a **genuinely independent** recomputation of
`{D0,P0,R_norm,R_pole,P2,P4}` from a fixed `direct_coefficients` bundle — Mathematica (`mt15_05` clean-R_norm
assembly, re-pointed; it CAN independently verify this) is the preferred second engine; if you instead use a second
Python path it must be an independent derivation, NOT a transliteration (B2b's round-1 sin). Independent root-find on
the same `R_norm(τ)` samples too. Agreement to a stated abs+rel tol.

---

## A. Acceptance criteria (all must hold before you report DONE)
1. R1 satisfied: every reported `R_norm(τ)` traces to a real re-solve at that τ (background residual logged).
2. Root `τ*` on the stable side `D0>0`, confirmed by a fresh full re-solve to `|R_norm(τ*)|≤tol` (R3).
3. τ_floor reported with residual; nothing extrapolated below it; all failed τ logged (R4).
4. §J error bars on `τ*, D0(τ*), R_pole/P2/P4` from the RECORDED budgets (R5).
5. Naturalness reported (`|ln τ*|`, cancellation ratio, digit count) + the frozen-vs-resolved negative-control
   contrast with measured `{B,Z,N}(τ)` drift (R6).
6. Target-blind firewall intact: τ calibrated only on `R_norm`; held-out reported as predictions, never tuned to;
   no new DOF (R7).
7. Dual-engine agreement on the assembly/observable/root-find algebra (R8).
8. Tests pass (add B2c tests incl. real negative controls: a frozen-background root-find must give a DIFFERENT τ*
   from the full re-solve when `{B,Z,N}(τ)` drift is injected; the confirmation gate must FAIL on a deliberately
   mis-located τ*; the stable-side filter must reject a `D0≤0` root).
9. Report `reports/patha_b2c_calibration_report.md` written, honest about every soft spot.

## B. Discipline (unchanged)
- Codex codes/designs; Claude reviews. Do not have Claude pre-design the route.
- Iterate scripts until they exit 0 ([[feedback-codex-iterates-until-clean]]); the verifier reviews substance after.
- Mathematica ≤2 concurrent `math -script` seats; GPU off → CPU sparse-direct.
- YAML/markdown for any LLM-read file; JSON only for machine bundles.
- Write outputs under `runs/patha_b2c_calibration/` (gitignored bundles ok). Never touch
  `research/pde_audit/simulation/` or the `physical_export_permitted` guard.
- Commit only when the user asks (orchestrator handles it after review).

## C. Review stack (orchestrator runs after Codex reports clean)
1. Transliteration-fidelity audit (1 clean agent) on any NEW math module (assembly/observable/root-find/error-prop):
   code term-by-term vs decision-11 §5 / decision-12 formulas. (The re-used B2a/B2b operators are already audited.)
2. Adversarial audit (≥1 round, distrust-all-clean — B2b needed 3): hunt for frozen-background masquerading as
   re-solve, a can't-fail confirmation gate, a fitted-model that wasn't confirmed, silent τ truncation, a
   transliteration posing as the second engine, and any held-out leakage into the τ calibration.
