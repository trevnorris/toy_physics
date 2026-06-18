# Directive — Path-A chunk 1c REMEDIATION (honesty relabel + one genuine-independence strengthening)

**Program:** Path-A closed solver, chunk 1c. The 1c build is COMPLETE and the operator is fidelity-correct; the
closed solve converges (3 levels, orders ~2.0/1.83), conservation is intact, Part B is clean, 62/62 tests pass.
Two independent clean-agent reviews (transliteration-fidelity + adversarial) found **NO fatal flaw** but flagged
**overclaim/labeling** in exactly the place the 1c directive made load-bearing. This directive fixes that. **You
(Codex) CODE; Claude re-verifies.** Sandbox workspace-write; every script you RUN gets `timeout 600`, iterate to
exit 0. CPU-only. Still target-blind: NO calibration / `R_norm` / `D0` / export / GATE-A forms; firewall
(`research/pde_audit/simulation/`, `physical_export_permitted`) untouched.

**This is honesty-relabeling + ONE new genuine check. Do NOT change the operator, the solver, the convergence
ladder, the conservation math, or the 1b harness.** Only: add the new flux check, re-route/relabel gates, fix
report wording, adjust two thresholds.

## The findings you are fixing (from the two reviews)
- The two "independent recomputes" (`independent_source_recompute` ≈L132, `independent_flux_divergence_recompute`
  ≈L144) are **verbatim copies** of the operator paths (`patha_radial_reduced_return_source` /
  `static_balance_terms` flux block) → diff ≡ 0.0 by construction. They catch edit-drift, NOT a shared-wrong
  operator — so calling them "independent reconstructions / not a tautology" (report + module `diagnostic_note`
  ≈L243, `method.non_tautology` ≈L940) OVERCLAIMS.
- `balance_residual_relative_below_solver_floor` (≈L228) and `closed_operator_gauss_residual_at_solver_floor`
  (≈L855) are true-by-construction once Newton converged (the wall residual is a slice of the Newton residual;
  the operator-Gauss residual is bounded by `final_residual_linf·volume`).
- Conservation envelope `conservation_baseline_factor=10.0` (≈L101) is arbitrary; 3 of 5 metrics are structurally
  0.0 so only the center-gradient Gauss metric is live.
- MODERATE: "R0-aware conservation" overstates — number/charge/Gauss ARE R0-aware (functions of the solved
  {ψ,A}), but the energy-density confinement term uses FROZEN geometry (`_energy_density` →
  `confinement_potential_torch` with no `radius=`). The closed-vs-frozen comparison is still apples-to-apples
  (both use frozen for that term).
- Nice-to-have: `nontrivial_term_floor=1e-10` (≈L97) is ~3e5× below the smallest real term; `null_floor_label`
  keys on transport-signal norm not residual (moot on this isotropic branch).

## What to change

### 1. Add a GENUINELY independent flux check (the real faithful-but-wrong-operator guard)
Add a second, DIFFERENT discretization of the continuum wall operator `−∂_w(T_{w,Σ} ∂_w R0)` and compare it to the
operator's conservative `StaticBalanceTerms.flux_divergence` — agreement that DECREASES under grid refinement is
the genuine check (two independent discretizations of the same continuum operator must converge to each other).
- Use the NON-conservative expanded form, built from the existing center-difference stencil `wall_center_gradient`
  (patha_static_balance.py): `g = wall_center_gradient(R0)`; `q = T_w(center, w_center)·g`;
  `nonconservative_flux_div = −wall_center_gradient(q)`. This shares NO code with the operator's conservative
  face-flux form (`F[j+½]=T_w(face)ΔR/dw`, one-sided mouth stencil, zero-traction exit).
- Compare to the operator's `flux_divergence` on INTERIOR wall cells only (the two forms legitimately differ at
  the mouth one-sided stencil and the zero-traction exit cell — exclude those and SAY SO in the report). The
  interior max-abs diff is O(h²) and must DECREASE under refinement.
- Make this the COUNTED gate `independent_flux_discretization_converges`: compute the interior diff at every ladder
  level, report the per-level diff + its observed order (reuse `observed_order_from_three`), and PASS iff the diff
  decreases (or is already at the solver floor) AND the finest-level interior diff ≤ a stated, refinement-aware
  tol. This is nonzero-but-converging, like the conservation center-gradient Gauss check.
- KEEP the existing verbatim `independent_source_recompute` and `independent_flux_divergence_recompute` but
  DEMOTE them to `_not_a_physics_gate` / identity diagnostics, honestly labeled "parallel reconstruction of the
  same kernel — edit-drift / wiring pin, diff==0 by construction; does NOT establish operator correctness." For
  the SOURCE (a pure exact reduction, not a differential operator — a different discretization is not meaningful),
  this verbatim pin + the physical reciprocity (same kernel as the forward wall→matter force) + the chunk-1a
  dual-engine MMS + the fidelity audit are what establish correctness; state that explicitly. Do NOT invent a fake
  "independent" source discretization.

### 2. Route the two solver-floor restatements out of the counted physics gates
Move `balance_residual_relative_below_solver_floor` and `closed_operator_gauss_residual_at_solver_floor` into
`identity_checks` (the `_not_a_physics_gate` bucket), labeled as solver-floor diagnostics (true-by-construction
once Newton converged). Keep reporting their numbers; just don't count them as independent physics gates. The
counted balance evidence becomes: `balance_terms_nontrivial` + the new `independent_flux_discretization_converges`
(+ the convergence-order and conservation gates).

### 3. Expose / justify the conservation envelope
For `_run_conservation_validation`: report, per metric, whether it is LIVE (nonzero on this branch) or STRUCTURALLY
ZERO, so the gate's actual discriminating power is visible. Tie the live-metric tolerance to something principled
or document the choice: replace the bare `conservation_baseline_factor=10.0` with a documented, smaller factor
(e.g. justify a modest factor as "closed and frozen equilibria differ only by the self-consistent R0 adjustment,
expected ≲ a few × the baseline") OR state in the report that the live Gauss metric's real evidence is its
decrease under refinement (already shown), with the baseline comparison a sanity bound. Feature the
center-gradient Gauss reconstruction as the genuinely-independent conservation signal.

### 4. Fix the "R0-aware conservation" wording (MODERATE)
In the report and the `method.conservation_framing` note: state that number/charge and the independent Gauss
closure ARE evaluated on the solved {ψ,A} (R0-aware), while the energy-density confinement term is evaluated on
frozen geometry (`confinement_potential_torch` is called without `radius=`), and the closed-vs-frozen comparison is
apples-to-apples because both branches use the same frozen confinement term for that channel. No code change to the
conservation math — wording only (unless you choose to thread `radius=R0` into the closed-state energy term, which
is OPTIONAL and must keep the comparison consistent; if unsure, just fix the wording).

### 5. Nice-to-haves
- Raise `nontrivial_term_floor` to a principled bar (e.g. a small fraction of the dominant term, or keep it but
  relabel it a "degeneracy tripwire" rather than a "non-degeneracy bar"), so the gate's role is honest.
- Add a one-line note that `null_floor_label` keys on the transport-signal norm (not residual magnitude), which is
  correct/moot on this isotropic stationary branch but would need revisiting for a current-carrying branch.

## Acceptance criteria
1. `independent_flux_discretization_converges` is a counted gate built from a genuinely different discretization
   (no shared code with the operator's conservative flux), interior-only, with per-level diff + observed order
   reported and decreasing/at-floor.
2. The verbatim source/flux recomputes are demoted to honestly-labeled `_not_a_physics_gate` edit-drift pins; the
   "independent reconstruction / not a tautology" overclaim is removed from the report and module notes.
3. The two solver-floor restatements are moved to `identity_checks`, labeled true-by-construction.
4. The conservation envelope is exposed (live vs structurally-zero metrics) and the factor justified or
   documented; the R0-aware wording is corrected.
5. `nontrivial_term_floor` raised/relabeled; the `null_floor_label` note added.
6. Re-run `python -m stage1_solver.patha_closed_validation` (regenerate the 1c report) and full
   `pytest software/stage1_solver/tests -q` — all pass; target-token scan clean over the 1c files; firewall +
   export guard untouched; 1b harness/report UNCHANGED by this directive.

## Stop conditions
If the new non-conservative flux discretization does NOT converge to the operator's conservative flux under
refinement (interior diff not decreasing), STOP and report it — that would be a real operator/discretization
discrepancy (a finding), not something to silence.

## End your final message
with: the new flux-convergence gate's per-level interior diffs + observed order; the final counted-gate list vs
the `_not_a_physics_gate`/identity list (showing what moved); the conservation live-vs-zero metric breakdown + the
factor decision; confirmation the report overclaim wording is gone and the R0-aware wording fixed; confirmation 1b
is untouched and full tests pass; any stop condition tripped.
