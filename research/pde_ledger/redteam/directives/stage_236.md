---
unit_id: 236
batch: VII.2
created_at: 2026-06-02T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: true
verification_status: pending
needs_user_resolution: false
applied_at: 2026-06-02T22:14:34-06:00
findings_applied: 1
findings_blocked: 0
---

# Codex directive — unit 236

Apply the finding below. After applying, append an `## Applied: F1` block under it with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes beyond what the finding names. Do NOT touch paper.tex, notes/, or any prose documents.

After writing the script, RUN it (`math -script <path>`) and iterate until it exits 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

## F1 — missing_verification_script (subtype missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_mathematica_audit.wl`

**Issue:** Stage 236 has a SymPy audit but no Mathematica audit (the card's `\stagefield{Verification}` says "Mathematica audit: none yet"). The stage is non-checkpoint and non-status-only, and all of its claims are finite-dimensional linear algebra over the rational-function field in `eps_eta_star`, `chi0_star`, `F_star` — fully verifiable in native Mathematica. The dual-engine rule therefore requires a second-engine script. Write a NEW Mathematica script that verifies the claims by an INDEPENDENT route, not a line-by-line port of the `.py`.

**Independence requirement (mandatory):**
The existing SymPy script verifies the projectors by *constructing* `P_nt_dep = S_rm_dep . diag(1,0) . L_rm_dep` and `P_eta_dep = S_rm_dep . diag(0,1) . L_rm_dep` and comparing to explicit matrices, then checks idempotency by matrix multiplication. The Mathematica script must reach the same conclusions by a DIFFERENT decomposition. Acceptable independent routes (Codex chooses; do not transliterate the SymPy choreography):
- Verify the compiler/section claims by acting on a symbolic generic input vector and comparing component expressions (e.g., apply `C_rm_dep` to `{R1, E1}` and to a basis, rather than comparing the matrix product to a hardcoded matrix), and/or
- Establish the left-inverse and projector relations by solving the linear recovery `{q_nt, q_eta}` from `y_rm` via `LinearSolve`/`Solve` on the plane `Delta_T = 0` and showing the resulting recovery map equals `L_rm_dep`, then derive the projectors as `S_rm_dep . (selector) . (recovered map)` from that solved inverse rather than asserting `L_rm_dep` by hand, and/or
- Verify idempotency/complementarity spectrally (e.g., projector eigenvalues are `{0,1}` with the correct ranks, `P_nt_dep + P_eta_dep` restricted to the plane is the identity) instead of by repeating the same `P^2 - P = 0` matrix product.
Use native primitives: `Dot`, `Transpose`, `IdentityMatrix`, `DiagonalMatrix`, `LinearSolve`, `Eigenvalues`, `Together`/`FullSimplify`. Declare `eps_eta_star`, `chi0_star`, `F_star` real with `0 < eps_eta_star < 1` so `c_eta = eps_eta_star/(1 - eps_eta_star)` and `1/(1 - eps_eta_star)` are well-defined and `1 - eps_eta_star != 0`. Treat `q_nt, q_eta, q_tr, R1, E1` as free real symbols.

Use a zero-test helper that strips `ConditionalExpression[0, ...]` (per project Mathematica idioms) and exit nonzero on any nonzero residual.

**Claim manifest** (the new script must independently verify each):
- M1 — Rigid-mouth packet map: with `M_rm = {{-1, -c_eta}, {0, 1}}`, `M_rm . {R1, E1} = {-R1 - c_eta E1, E1}`, where `c_eta = eps_eta_star/(1 - eps_eta_star)`.
- M2 — Dependent section on the rigid-mouth slice (`q_tr = 0`): the general section `(Delta_T, Delta_Keta, Delta_mu) = (q_tr/(1+chi0_star), -q_eta, F_star q_tr/(1+chi0_star) + q_nt - q_eta)` reduces at `q_tr=0` to `y_rm = (0, -q_eta, q_nt - q_eta)`, and equals `S_rm_dep . {q_nt, q_eta}` with `S_rm_dep = {{0,0},{0,-1},{1,-1}}`.
- M3 — Direct-observable-to-dependent compiler: `C_rm_dep = S_rm_dep . M_rm = {{0,0},{0,-1},{-1, -1/(1 - eps_eta_star)}}`, and `C_rm_dep . {R1, E1} = {0, -E1, -R1 - E1/(1 - eps_eta_star)}`.
- M4 — Left inverse on the plane `Delta_T = 0`: `L_rm_dep = {{0,-1,1},{0,-1,0}}` satisfies `L_rm_dep . S_rm_dep = IdentityMatrix[2]`, and `L_rm_dep . y_rm = {q_nt, q_eta}`. (Per independence requirement, derive/confirm `L_rm_dep` by solving the recovery rather than asserting it by hand.)
- M5 — Dependent-plane projectors: `P_nt_dep = {{0,0,0},{0,0,0},{0,-1,1}}` and `P_eta_dep = {{0,0,0},{0,1,0},{0,1,0}}` are idempotent (`P^2 = P`), mutually annihilating (`P_nt_dep . P_eta_dep = 0` and `P_eta_dep . P_nt_dep = 0`).
- M6 — Plane completeness: `P_nt_dep + P_eta_dep = DiagonalMatrix[{0,1,1}]` (the identity on the plane `Delta_T = 0`).
- M7 — Decomposition: `y_nt = P_nt_dep . y_rm = {0,0,q_nt}`, `y_eta = P_eta_dep . y_rm = -q_eta {0,1,1}`, and `y_rm = y_nt + y_eta`.
- M8 — Static-strip equivalence and equal-drift ray: `Delta_mu - Delta_Keta = q_nt` (so `q_nt = 0 ⟺ Delta_mu = Delta_Keta`); on the strip `R1 = -c_eta E1`, `C_rm_dep . {R1, E1} = {0, -E1, -E1} = (R1/c_eta) {0,1,1}`; norms `y_eta . y_eta = 2 q_eta^2`, and on the strip the dependent vector's squared norm `= 2 E1^2 = 2 R1^2 / c_eta^2`.
- M9 — Microscopic correction compilers: `Delta_y_static = -y_nt = {0,0,-q_nt}` with `y_rm + Delta_y_static = y_eta`; `Delta_y_orbit = -y_rm = Delta_y_static + q_eta {0,1,1}` with `y_rm + Delta_y_orbit = 0`.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 236` and confirms the new `.wl` exists, asserts M1-M9 by an independent route, and exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_mathematica_audit.wl`
- summary: Added the missing independent Mathematica audit for M1-M9, deriving the plane recovery map with LinearSolve and checking the dependent projectors and correction identities.
- deviation: none
