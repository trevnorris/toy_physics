# Directive pathA_13 — Chunk B2b remediation #2 (close the unmet full-rigor requirements)

**Owner:** Codex (codes + DESIGNS the route + runs + iterates). Claude reviews afterward. **Status:** B2b build +
remediation #1 (`pathA_12`) DONE. A fidelity audit + an adversarial pass (independently) confirm: the operator is
FAITHFUL and untouched, target-blind is CLEAN, the N-fit reconditioning is legitimate (value-preserving column
equilibration, cond 5.5e6→33), the forward-source provenance is corrected, and the grid + ω-window sweeps genuinely
converge — **keep all of that, do NOT touch it.** But **three of the user's "full rigor" requirements are NOT
actually met.** This directive closes exactly those. They are gate-quality / convergence fixes — they do NOT change
B2b's conceptual nature (still: derive `{Z_n,N_n}` on the Path-A background, target-blind).

Constraints unchanged: `timeout 600` per script (exit 124 = reformulate, never raise), ≤2 `math -script` seats, CPU,
firewall (no writes/imports under `research/pde_audit/simulation/`, don't touch `physical_export_permitted`;
`research/pde_audit/{scripts,notes}/` READ-ONLY), target-blind (NO `R_norm`/`R_pole`/`P2`/`P4`/`D0`/`P0`/root-find/
`mt15_05`), no `git add`/commit. Do NOT re-derive or alter the FAITHFUL operator/transfer/source math.

## F1 (PRIMARY — the radial truncation is NOT converged; the gate hides it)
Evidence (from `runs/.../python/patha_b2b_python_tau_1_nr_22_nw_17_w_0.028_rs_{2,2.25,2.5}.json`): the per-coefficient
radial-truncation increments **GROW** (`N0` 0.45%→0.49%, `N4` 2.12%→4.47%), all six coefficients march monotonically
in one direction, and **nothing beyond `radial_scale=2.5` was ever run** — so 2.5 is where the 3-level sweep was
truncated, not a demonstrated plateau. The convergence gate passes only because `_sweep_convergence_summary` tracks
the **max over coefficients** per step, and the dominant coefficient's identity flips (Z4→N4) between the two steps,
so `0.0447 ≤ 0.98·0.0479` passes while `N0` and `N4` are individually still growing. Fix:
- **Extend the radial domain OUTWARD over ≥3 NEW levels beyond 2.5** (e.g. `3.0, 4.0, 5.0, …` as needed) until EACH of
  `{Z0,Z2,Z4,N0,N2,N4}` has a per-step relative increment that **shrinks below the stated tolerance** AND the final
  step moves less than the tolerance (a genuine plateau — show the value stops moving).
- **Rewrite the convergence gate to require PER-COEFFICIENT shrinking increments** (each coefficient individually:
  monotone-decreasing increments converging below tolerance). A **max-over-coefficients aggregation is explicitly
  FORBIDDEN** — it is exactly what let a growing `N4` through. Add a negative-control TEST: a sequence where one
  coefficient shrinks while another GROWS must FAIL the gate.
- **If the coefficients do NOT plateau even at large domains**, the open-boundary / DtN is leaking (not a true
  non-reflecting condition — an exact DtN is domain-size-independent). Fix the absorbing/impedance treatment per
  `docs/boundary_and_noise_methods.md` (real sponge ≠ absorber; true outgoing-wave impedance on the radiative
  tangent). If that is a reformulation infeasible under the constraints, HALT and surface to the user — do NOT ship a
  marching tail as "converged."
- Carry the HONEST per-coefficient error bar = the residual at the demonstrated plateau (or the last increment if
  larger) into the bundle `error_budget`.

## F2 (the "independent engine" is NOT independent — wire in the genuinely-different one that already exists)
The live correctness-certifying "independent" engine is `staggered_second_order`, which shares the IDENTICAL derivative
stencils (`_first/_second_derivative_matrix`), quadrature form (`dr·dw·r²·z`), weak-form assembly, anchor, ω-window,
and fit with the primary — differing ONLY by node offset (`+0.5` cell-center vs `+1` interior) and spacing (`/nr` vs
`/(nr+1)`). That is a grid-shift sensitivity check, not a different discretization; the ~3.7% agreement carries
grid-placement information, not cross-scheme information. Fix:
- **Wire the genuinely-different `staggered_high_order` engine** (the true 4th-order stencils via
  `_high_order_derivative_matrix`, already present at ~py:389-394 but never invoked in the validation path — the
  default maps `staggered`→`staggered_second_order`) in as the correctness-certifying independent engine. Regenerate
  the independent packets at the converged resolution.
- **Subject the independent (4th-order) engine to the SAME full sweeps** (grid + outward-radial + ω-window) as the
  primary, so its error bar is comparable (currently it is grid-only — fix that). Report both engines' per-coefficient
  error bars on the same axes.
- Both engines must converge to the same `{Z_n,N_n}` to a tolerance reflecting combined discretization error (report
  abs AND rel per coefficient). Keep the dual-engine AND-gate, but its rel tolerance must be tight enough to still
  catch a real bug while honestly admitting the cross-scheme spread; do NOT loosen it to rubber-stamp a large
  disagreement. If `staggered_high_order` shares anything load-bearing with primary, build a genuinely different
  scheme (spectral/Chebyshev or FEM) or HALT to the user.
- The old Python↔MMA transliteration pair stays labelled transcription-only (already done — keep it).

## F3 (the basis-invariance gate is a dead conditional — make it a real test)
`basis_invariance_check` has byte-identical `if (broken or inject_port_leak): … else: …` branches (the flag injects
NOTHING), and the negative control "fails" only because `passed = … and not (broken or inject_port_leak)` — a
hardcoded boolean. The invariance result is still the algebraic identity: a unitary `Q` applied identically to
operator + source + boundary gives `⟨Qᴴs,(QᴴAQ)⁻¹Qᴴs⟩ = ⟨s,A⁻¹s⟩` for ANY operator, right or wrong. Fix:
- The `inject_port_leak` branch MUST construct an actual **basis-DEPENDENT** extraction — a posited-U/W-fixed-port
  functional that reads off a FIXED lane basis (the very thing decision-05 D4 rejects) — and the gate MUST FAIL on it.
  The two branches MUST compute genuinely different quantities (verify: toggling the flag changes the numbers).
- The invariance result MUST demonstrate that the basis-invariant `Z=⟨j,G j⟩` gives the SAME `{Z_n,N_n}` under a
  genuine lane re-basis (a change of lane basis / gauge frame), WHILE a fixed-port extraction MOVES under that same
  re-basis. The check must be capable of failing for a wrong (basis-dependent) extraction — not a global similarity
  transform that is invariant by algebra for any operator.
- The negative-control TEST must inject a real port-leak extraction and assert the gate FAILS (not a hardcoded
  boolean). State the specific wrong answer the gate now catches.

## F4 (report honesty + N-moment floor checks)
- State that the primary and independent error bars are measured on the SAME axis set (after F2); do NOT present a
  grid-only independent bar next to a full-sweep primary bar as if they were comparable.
- State the cross-engine agreement honestly relative to the convergence error (if both engines are ~X% converged and
  agree to ~X%, say so — do NOT frame it as a clean pass beyond what the numbers support).
- The N-channel robustness gate currently floor-checks only `N0`. Add individual above-floor / fit-stability checks
  for `N2` and `N4` (N4 is the least-converged, the quartic term of a near-flat ω-signal).

## Acceptance criteria (Codex iterates until ALL pass, exit 0)
1. Radial domain extended outward ≥3 new levels past 2.5; EACH coefficient's increment shrinks below tolerance with a
   demonstrated plateau (final step < tol). Convergence gate is PER-COEFFICIENT (max-over-coefficients forbidden) with
   a negative-control test that fails on grow-one/shrink-another. Honest per-coefficient error bars in the bundle.
   (Or a documented HALT-to-user if the open boundary cannot be made to plateau.)
2. The correctness-certifying independent engine is the genuinely-different 4th-order (`staggered_high_order`) scheme,
   run through the SAME full sweeps as primary; both converge to the same `{Z_n,N_n}` (abs AND rel per coefficient
   reported); dual-engine AND-gate tight enough to catch a real bug. (Or HALT if a genuinely different scheme is
   infeasible.)
3. Basis-invariance gate genuinely tests basis-dependence: the `inject_port_leak` branch computes a different
   (basis-dependent) quantity and the gate FAILS on it; toggling the flag changes the numbers; negative-control test
   is real (not a hardcoded boolean).
4. Report text honest on comparable error bars + cross-engine-vs-convergence; N2 and N4 individually floor/stability
   checked.
5. Full `pytest` for `test_patha*.py` green; target-blind intact (no `R_norm`/…/`mt15_05`); firewall untouched; no
   commit. Operator/transfer/source math unchanged from the FAITHFUL version.

## Report back
The converged-with-plateau `{Z0,Z2,Z4,N0,N2,N4}` at `τ=1` (and the outward-radial ladder showing each coefficient
plateau); the new per-coefficient convergence-gate criterion + its grow-one/shrink-another negative control; the
4th-order independent engine's full-sweep convergence + the two engines' per-coefficient converged-value agreement
(abs+rel); the genuine basis-invariance test + its real port-leak negative control (with proof the two branches
differ); the N2/N4 floor checks; the honest error-bar framing; and every file modified.
