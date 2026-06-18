# Directive — Path-A build Chunk 1b: closed self-consistent Newton over {ψ, A, R0}

**Program:** Path-A self-consistent closed solver, per `decisions/10_pathA_solver_build_design.md` (design +
the RESOLVED pre-1b specs). **You (Codex) CODE this; Claude reviews + runs a transliteration-fidelity audit +
adversarial review.** Sandbox: workspace-write; scripts you RUN get `timeout 600`; iterate to exit 0. CPU-only.
**Target-blind ENGINEERING build:** placeholder `S_Σ` (the `smooth_positive_placeholder_v1` family), parameters
AWAY from the `D0→0` near-pole, **NO calibration, NO `R_norm`/target comparison, NO physical export, NO GATE-A
forms.** NEVER touch `research/pde_audit/simulation/` or the `physical_export_permitted` guard.

## Read first
- `decisions/10_pathA_solver_build_design.md` — esp. the **Pre-1b items RESOLVED** section (the source sign,
  radial reduction, exit BC) and the 1b spec.
- `derivations/pathA_01_return_source_and_balance.md` (D1 matter return, D2 gauge = matter-mediated/no separate
  kernel, D3 balance, D4 no-numerator-knob).
- The 1a operator you extend: `src/stage1_solver/patha_static_balance.py` (the static-balance operator + the
  `S_Σ` registry/providers).
- The existing solver to couple into: `src/stage1_solver/coupled_branch.py` (the `{ψ,A}` residual, the
  confinement coupling / `geometry_radius_torch:124`, `_matter_number_current:210`, Maxwell residual),
  `newton.py`, `preconditioners.py:195` (layouts `5*cells+1` / `5*cells+nw`), `p2_tangent.py:165` (the
  confinement R-derivative coefficient `4 V_radial r⁴/R0⁵`).

## What to build (1b only — the closed static solve; the fuller validation is 1c)

1. **Promote `R0(w)` into the coupled residual.** State vector `[ψ_R, ψ_I, A0, Ar, Aw, R0, μ]`. Add `R0` (1D
   wall-grid unknown on `grid.w_centers`) to pack/unpack, the initial guess, and resampling. Keep the
   effective-closure (frozen-geometry) path INTACT and selectable — this is an ADDITIVE closed-path mode.
2. **The wall residual = the 1a static-balance operator with the PHYSICAL return source on its RHS** (RESOLVED
   specs, decision-10):
   - source `S(r,w) = −V_conf,Σ(r−R0(w))·ρ0(r,w)` — the FULL `ρ0` (absolute, no `ρ_ref` subtraction);
   - radial reduction to the wall grid: `S_j = Σ_i ΔV_i^r·[−V_conf,Σ(r_i−R0(w_j)) ρ0_ij]`,
     `ΔV_i^r = (4π/3)(r_{i+½}³−r_{i−½}³)` (the solver's shell volumes, `grid.py`), in the flat-`dw` wall
     convention (`operators.py:485`);
   - **RECIPROCITY (load-bearing): use the SAME confinement derivative `V_conf,Σ` object that the matter
     equation already uses for the forward wall→matter force** (do NOT re-implement a second copy) — so forward
     and return share one kernel and the Jacobian's `R0↔ψ` block is the literal adjoint pair.
   - mouth BC `R0(0)=a` (FV Dirichlet); **exit BC = natural zero-traction `T_{w,Σ}(R0(L),L) R0'(L)=0`** (the
     `Y_L=0` open limit), `R0(L)>0`.
3. **Gauge return: NO separate kernel.** The matter-mediated `S_η^(A)` is captured by the monolithic coupled
   Newton (the existing matter/Maxwell/current equations) — confirm no double-counting of the matter current
   (`compact:627-630`). Do not add an explicit `η→A` source.
4. **Couple into the Newton/JFNK.** Extend the residual Jacobian with the new `R0` block + its `R0↔{ψ,A}`
   coupling; adapt the colored sparse preconditioner to the `5*cells+nw+1` layout (`preconditioners.py:195`).
   Keep the operator a true JVP (preconditioner is `M`, solution-preserving — the step-3b lesson).
5. **Run the closed solve** at a placeholder `S_Σ` away from the pole; report Newton/GMRES convergence + the
   final residual; confirm `R0(w)>0` throughout and a genuine converged equilibrium (not a partial iterate).

## Acceptance criteria (target-blind; genuine gates only, restatements split as `_not_a_physics_gate`)
1. **JVP consistency:** the analytic/AD residual JVP matches a finite-difference directional derivative
   (incl. the new `R0` rows/cols) to a stated tol — a GENUINE gate (would fail on a wrong Jacobian coupling).
2. **Closed-Newton convergence:** documented residual history dropping to a stated tol; GMRES iteration counts
   sane with the preconditioner; `R0>0`, `D0>0`/stability respected (you are away from the pole).
3. **Return-source fidelity (self-check):** the wall residual's source equals the radial-reduced
   `−V_conf,Σ·ρ0` built from the SAME confinement kernel as the matter force (reciprocity) — exposed as a
   diagnostic, not asserted.
4. **Placeholder-provider derivative validation (carries MINOR-2 from 1a):** directly validate the
   `smooth_positive_placeholder_v1` provider's `T_w_R/T_w_RR/U_R/U_RR` against finite-difference of its
   `T_w/U` (this family is what 1b actually solves with; 1a's MMS did NOT exercise it).
5. **Additive / no regression:** full `pytest software/stage1_solver/tests -q` still passes; effective-closure
   paths unchanged.
6. **Target-token scan clean**; no calibration; firewall + export guard untouched.

## Stop conditions (decision-10 — HALT and report, do NOT paper over)
- `R0` goes nonpositive anywhere; or the exit BC has to be changed from zero-traction to force convergence; or
  Newton converges only after hidden clamps/line-search hacks; or the source sign/radial reduction does not
  match the resolved spec. If any trips, STOP and report it (it is a finding, not a thing to silence).

## Recurring sins to avoid
No tautological/can't-fail gate in the genuine list; no circular self-check (the return-source diagnostic must
recompute from the kernel, not echo the residual); no hardcoded values; no silent fallback that degrades to
pass; honest labels (engineering smoke at placeholder `S_Σ`, NOT physics).

## End your final message
with: what you added (the `R0` coupling, the return source, the preconditioner-layout change), the closed-solve
convergence result (residual history, GMRES counts, `R0>0`, stability margin), the JVP-check + placeholder-
derivative-check results, the genuine-gate list, and any stop-condition that tripped or anything in the resolved
specs that was under-specified for coding. This is your report to the reviewer.
