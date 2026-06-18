# Directive — Path-A build Chunk 1a: S_Σ spec/registry + static-balance operator + MMS (dual-engine)

**Program:** Path-A self-consistent closed solver, per `software/stage1_solver/decisions/10_pathA_solver_build_
design.md`. **You (Codex) CODE this; Claude reviews + runs a transliteration-fidelity audit.** Target-blind
Stage-1 build: **NO calibration, NO `R_norm`/target comparison, NO physical export, NO GATE-A forms.** Do NOT
touch `research/pde_audit/simulation/` or the `physical_export_permitted` guard. Sandbox: workspace-write. Any
script you RUN gets `timeout 600`; iterate until exit 0. CPU-only (GPU off).

## Scope of 1a (ONLY this — do NOT start the closed Newton / 1b)

Build and MMS-validate the **static throat-balance OPERATOR** (the LHS) + the **parameterized `S_Σ` interface**.
1a uses a MANUFACTURED source `S_j` (MMS) and Dirichlet ends — it does NOT need the physical return source or the
exit-BC/`δρ` decisions (those are 1b). Read first: `decisions/10_*` (the design),
`derivations/pathA_01_return_source_and_balance.md` (D3 the balance + `K_η` reduction), and the existing code it
extends (`src/stage1_solver/operators.py:485` wall-FV pattern, `coupled_branch.py`, the step-2 MMS harness).

## Deliverables

### A. Serializable `S_Σ` constitutive spec + registry
A hashable, serializable spec (NOT raw callables in frozen config) that resolves to a provider of:
`mu(R,w), T_w(R,w), T_w_R(R,w), T_w_RR(R,w), T_Omega(R,w), U(R,w), U_R(R,w), U_RR(R,w)`
(`T_w_R = ∂T_w/∂R`, `U_R = ∂U/∂R`, `U_RR = ∂²U/∂R²`). Ship Chunk-1 defaults = **smooth, positive placeholder
families** + the MMS manufactured families. The registry must be reproducible/hashable (it will feed a freeze
hash later). NO GATE-A forms here.

### B. Discrete static-balance operator (conservative FV, faithful to decision-10 / derivation D3)
`balance_j(R0; S_Σ, S) = −(F[j+1]−F[j])/dw + ½ T_w_R(Rj,wj)(R0'_j)² + U_R(Rj,wj) − S_j`,
with face flux `F[j+½] = T_w(R_face, w_face)·(R[j+1]−R[j])/dw`. Match the existing wall-FV conventions
(`operators.py:485`) — same grid (`grid.w_centers`), same measure/`dw` convention — but this is a NEW nonlinear
`R0` operator (do not reuse the linear-`η` operator). Mouth BC: FV Dirichlet face value `R0(0)=a`. For 1a's MMS,
prescribe BOTH ends from the manufactured solution (the exit-BC choice is deferred to 1b). Keep it ADDITIVE — the
frozen-geometry helper + all effective-closure paths stay untouched.

### C. MMS validation (the load-bearing check — must be GENUINE, non-circular, dual-engine)
- Manufacture a **nonconstant** `R0(w)`, a **nonconstant** `T_w(R,w)` (so `T_w_R ≠ 0`), a **nonzero** `U(R,w)`
  (so `U_R ≠ 0`), and a regime with a **nonzero gradient-square** term `½T_w_R(R0')²`.
- Derive the continuum forcing `S(w) = −∂_w(T_w ∂_w R0) + ½ T_w_R (R0')² + U_R` **INDEPENDENTLY** — symbolically
  in **SymPy** AND in **Mathematica** (dual-engine; the two must agree). Do NOT copy the discrete operator into
  the MMS source (that is the circular-MMS sin — the source must come from the analytic continuum operator).
- Apply the discrete operator to the manufactured `R0` with that forcing; refine the grid (≥3 levels) and show
  **observed order ≈ 2** for the residual.
- **Term-by-term nonzero diagnostics:** each of the three terms (flux divergence, gradient-square, `U_R`) is
  individually exercised and nonzero in the test (so no term is silently untested).

### D. Honesty gates (this project's calibrated failure modes — avoid them)
- **No tautological / can't-fail checks presented as results.** Genuine gates (could fail) reported separately
  from any construction-restatements/guards (`_not_a_physics_gate`), exactly as the `pathA_01` verifier does.
- **Non-circular MMS** (source from the analytic operator, not the discrete one) — this is the #1 trap here.
- **No hardcoded values; no target leakage** — scan the new code + tests for `R_norm`, `54/5`, `10.8`,
  `P0_target`, GR-constant tokens (none).
- **Additive / no regressions:** the existing test suite (`pytest`) still passes; effective-closure paths
  unchanged.

## Acceptance criteria
1. `pytest` for the new module passes; the full existing suite still passes (additive).
2. MMS order ≈ 2 for the static-balance operator, with SymPy and Mathematica forcing agreeing.
3. Term-by-term nonzero diagnostics pass (every term exercised).
4. Genuine-vs-restatement gate split is explicit; no tautological gate in the genuine list.
5. Target-token scan clean.
6. A short report (committed-quality markdown at `software/stage1_solver/` top level or `mathematica/`) listing
   the genuine gates + results + the MMS orders; run artifacts under `runs/` (gitignored) are fine.

## End your final message
with a concise structured summary: what you built (the spec/registry + operator), the MMS result (orders, the
SymPy↔Mathematica agreement), the genuine-gate list, and anything in decision-10 / derivation D3 that was
under-specified for 1a (there should be none — the 3 open items are 1b). This message is your report to the
reviewer.
