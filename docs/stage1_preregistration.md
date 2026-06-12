# Stage-1 Branch-Realization — Pre-Registration Record

**STATUS: 🔒 FROZEN 2026-06-12** (user sign-off). This document is now **immutable**: no item below may be changed to rescue a missed target. A genuine spec error requires a new, openly-logged pre-registration, not an edit here. Freeze record in §L; non-rescue rules in §K.

**Purpose.** This record satisfies branch-realization prerequisite #3 (`docs/branch_realization_execution_plan.md` §8) and the brief's §3.4 pre-registration rule. It fixes — *in advance of any solve* — which branch is realized, the gauge/boundary/port conventions, the wall constitutive packet, the extraction formulas, the observables, the targets, and which observable is **primary**. It does **not** reconstruct any equation: every frozen formula is quoted from the source ledger with a citation (brief §2 mandate).

**Companion documents.** Scientific brief `docs/branch_realization_brief.md`; engineering plan `docs/branch_realization_execution_plan.md`; parent-status decision `docs/branch_realization_parent_status_decision.md`; freeze-boundary protocol `research/pde_audit/simulation/NONLINEAR_PROTOCOL_V2.md`; frozen-input template `research/pde_audit/simulation/ACTUAL_BRANCH_PROTOCOL_V1.md`. Canonical equation source: `notes/moving_throat_pde_program_compact.md` (cited as *compact* below).

---

## §A. Pre-registered parent-status decision

**`parent_action_status = effective_closure`** (Claude+Codex decision 2026-06-12, `docs/branch_realization_parent_status_decision.md`).

This Stage-1 test is the **effective-closure branch-realization test**: the quadratic wall closure `S_eta^(2)` is held as a fixed effective closure for this test; the strict parent action remains `S_current = S_psi[psi,A;Sigma] + S_EM[A]` with `Sigma/R` as a confinement-coupling argument (*compact* §13.1, L6783–6785). It is **NOT** a permanent program-wide abandonment of the autonomous-throat-field lane — promoting `S_Sigma[R]` (Path A) remains available later as a *separate, independently pre-registered* program with its own derivation + audit trail (non-rescue: a Path-A run may not be advertised as a continuation of this pre-registration).

---

## §B. Branch identity (what is solved)

- **WP1 — stationary isotropic moving-throat branch.** Stationary (drops `t`) + isotropic/spherically-symmetric (drops `Omega`), so fields are functions of `(r, w)` only on a 2D mesh parameterizing 3D-with-axial-symmetry physics (full `r^2` measure factor) — `branch_realization_execution_plan.md` §2.1.
- **WP3 — grouped real `P2` weak-axisymmetric tangent**, linearized around the converged WP1 branch; angular structure fixed by `Y_{2m}`, also on `(r, w)` (plan §2.1–2.2).
- Grouped real `P2` lanes: `A ∈ {20, 21, 22}`, full five-mode real bundle underneath (*compact* §12.3, L6182–6186).
- Geometry: open finite throat, level-set `Sigma(X,t) = r − R(Omega,w,t)`, surface `Sigma = 0`, `R(Omega,w,t) = R_0(w) + eta(Omega,w,t)` (*compact* §12.2, L6151–6163).

---

## §C. Fields / unknowns (frozen)

Parent fields (*compact* §12.1, L6118–6147): `psi(X,t)`, `rho = |psi|^2`, `A_M = (A_0, A_i)`, `F_{MN} = ∂_M A_N − ∂_N A_M`; brane velocity-potential `varphi`.

Coupled unknowns evolved on the open finite throat (`NONLINEAR_PROTOCOL_V2.md` "Unknowns", L21–31), on the `(r, w)` mesh:

- throat radius/support state `R(s)` with `R(0) = R_mouth` and finite open exit impedance at `s = L`;
- wall/support profile `chi_eta(s)`;
- BdG support profile(s) `phi_B(s)`;
- mixed port profiles `u(s)` and `w(s)`;
- continuation parameters for geometry, couplings, and exit impedance.

Madelung/axisymmetric working set (plan §2.3): `psi_real, psi_imag` (or amplitude/phase), `A_0, A_r, A_w`, `R_0(w)`, one gauge-fixing scalar / Lagrange multiplier — ≈5–7 coupled fields per cell.

---

## §D. Geometry, boundary, and gauge conventions (frozen)

| Item | Frozen choice | Source |
|---|---|---|
| Geometry | Open finite-throat branch; physical endpoint `R(L) > 0`; **no hard cap** | *compact* §13.1 L6788; `NONLINEAR_PROTOCOL_V2` Boundary Conditions; `ACTUAL_BRANCH_PROTOCOL_V1` Frozen-Inputs "geometry" |
| Mouth BC | Dirichlet anchoring for support profiles; fixed positive `R_mouth` | `NONLINEAR_PROTOCOL_V2.md` L35–37 |
| Exit BC | Open impedance / Robin condition with finite `R_exit > 0` | `NONLINEAR_PROTOCOL_V2.md` L38 |
| D/N ladder | AC impedance/reflection, **not** a hard cap | *compact* §13.1 L6789 |
| Mixed channels | `A_w, J^w, F_{mu w}, E_w, C_a` stay alive | *compact* §13.1 L6790 |
| Gauge | **H=Z** (finite localized gauge fixing, `xi_4 = xi`); clean-match `S = Z/Z_int → xi_eff = xi, mu_eff = mu0/Z_int` | *compact* §13.1 L6787, §13.1A L6803 |
| EM order | Projection first, reduction second (Stages 004–020 projected Maxwell + mouth transport) | *compact* §13.1A L6796–6804 |
| Outgoing-port convention | Natural point-particle source-map branch: `N_Q = chi_Q^{-1}`, `m̂0^2 · chi_Q · N_Q = 1` | *compact* §12.9 L6441–6445 |
| Forbidden | hard caps; closed-throat endpoint substitutions; any post-target retuning | `NONLINEAR_PROTOCOL_V2.md` L39–40 |

**Source-physics scope (frozen).** Active in this Stage-1 run: source/current terms generated by the frozen effective-closure WP1/WP3 PDE, projection-first Maxwell bookkeeping (`W, Z, H=Z, S=Z/Z_int`, boundary discharge, W' leakage, measured/flux split), the natural source-map export (`m̂0, S_port, chi_Q, N_Q`), and finite-throat `P0/P2` plus grouped `P2/P22` response only as extracted from the frozen WP3 tangent before residual evaluation. **Excluded** from this pre-registration: the post-242 relaxed/open-system companion lanes in *compact* §12.12 (`L_w, L_UV, L_varsigma`), including extra leakage/work lowering, non-rigid `U/V` dressing, compensated-source additions, microscopic export kernels, scalar-flux/open-radiative additions, and heat/material companions. They may not be added after residuals or used to reinterpret a miss.

---

## §E. Wall / support basis + constitutive packet (frozen)

Effective-closure wall action `S_eta^(2)` with constitutive functions held as **fixed posited inputs** for this test (not derived; red-team-classified `free_choice`): `mu_eta`, `T_w`, `T_Omega`, `K_eta` (`provenance/_synthesis/batch_01/fit_stage001_wall_action_constitutive_coefficients__mu_eta.yaml`, `__k_eta.yaml`). These are frozen **before** the solve and never refit after seeing residuals (§K).

Coherent reduced-kernel definitions consumed by the extraction (*compact* §12.5, L6245–6271):
```
epsilon_eta = c_etaU^2 / (K_U K_eta^eff)
epsilon_W   = gamma^2 lambda_W^2 sigma / (K_U K_W^eff)
Z_W         = lambda_W^2 / (K_eta^eff K_W^eff)
delta_U     = pi^2 T_U / (L^2 K_U)
chi_0       = gamma c_etaU / K_U
zeta        = lambda_phi^2 K_W^eff / (lambda_W^2 K_phi^eff)
R_tr        = (1 + chi_0/(1+delta_U)) / (1 + chi_0)
S(zeta;eps) = 1 + zeta(1-eps)/(1 - zeta eps),  eps = epsilon_W [1 - (2/11) delta_U/(1+delta_U)]
```
Free quintuple (*compact* §12.11 L6610): `y = (lambda_W, c_etaU, gamma, K_U, K_W^eff)`.

---

## §F. Extraction formulas (frozen, verbatim from ledger — applied AFTER the solve)

**Conservative coefficient bundle** (*compact* §13.2, L6809–6828):
```
D_0 = K − B_0 − Z_0
D_2 = −(M + B_2 + Z_2)
D_4 = −(B_4 + Z_4)
u_2 = −D_2/D_0
u_4 = (D_2^2 − D_0 D_4)/D_0^2
P_0 = N_0/D_0
P_2 = (D_0 N_2 − 2 D_2 N_0)/D_0^2
P_4 = (D_0^2 N_4 − 2 D_0 (D_2 N_2 + D_4 N_0) + 3 D_2^2 N_0)/D_0^3
```
Simulation aliases / one-pole surface (`ACTUAL_BRANCH_PROTOCOL_V1.md` L86–118):
```
D0 = K − B0 − Z0,   A = −D2 = M + B2 + Z2,   C = −D4 = B4 + Z4
one-pole surface:   D0*C/(3*A^2) = 1
constant-prefactor: N2 = −2 A N0/D0,   N4 = N0 (A^2 − 2 D0 C)/D0^2
on one-pole surface: N4 = −5 A^2 N0/D0^2
```
**Outgoing scalar / Packet-A finish line** (*compact* §12.9, L6436–6455):
```
chi_Q = 3(S beta^5 + 9 Sigma_5) / (3 S − Sigma_0)
N_Q = chi_Q^{-1},   m̂0^2 chi_Q N_Q = 1
Delta_Q := chi_Q − 1
Delta_norm = P0^target (chi_Q^{-1} − 1),   P0^target = 54 G c_s^5 / (5 a^5 c^5)
```
**Prefactor slope / weak-axisymmetric carrier** (*compact* §12.9 L6464–6476, §13.4 L6854–6856):
```
Xi_1 = P_1/P_0
a_{P0} = eps P̄0 Xi_1 / 4,   b_{P0} = 3 eps P̄0 Xi_1 / 4   (b = 3a, weak axisymmetric)
```
**Rigid-mouth physical chart `(U,V)`** (*compact* §12.8, L6404–6417):
```
U := ln(mathcal_T^2 / mathcal_T_ref^2),   V := ln(epsilon_eta / epsilon_eta,ref)
q_nt = U,   q_eta = V    (on the rigid-mouth slice q_tr = 0)
mathcal_T^2 = Z_W (1+chi_0)^2 / (Omega_W^2 (1-eps)^2)
R_target mathcal_T^2 = Lambda_0 (1 − epsilon_eta),   Lambda_0 = 27 pi^2 G c_s^5 / (20 a^5 c^5)
```
**Selected support-placement coordinate:** `varrho_phys`, computed from the stationary coherent support variable `epsilon = epsilon_W [1 − (2/11) delta_U/(1+delta_U)]` by
```
varrho_phys := pi^2 Pi_tr/(16 Lambda) = (2/3)(1 − epsilon)
```
(*compact* L4538–4546 for `Pi_tr = (16/pi^2) Lambda varrho`; L4334 for the closed form `varrho_phys = (2/3)(1−epsilon)`). The companion `sigma_phys = 2 epsilon/(1−epsilon) = 4/(3 varrho_phys) − 2` (*compact* L4336–4337) and the ranking region are derived readouts, not the registered coordinate.

---

## §G. Observables to extract (target-blind — brief §3.2)

At minimum: the coefficient bundle `(D0, D2, D4, N0, N2, N4, P0, P2, P4)`; source-map scalars `(m̂0, S_port, chi_Q, N_Q)`; the prefactor slope `Xi_1 = P_1/P_0`; the rigid-mouth orbit-lock chart `(U, V)`; the selected support-placement coordinate `varrho_phys`. Each carries a stated uncertainty (§J).

---

## §H. Targets (used ONLY for comparison, AFTER extraction)

**Do not insert into the solve, tune any parameter to, choose a branch by, or set a stopping criterion from these** (brief §3.3).

V2 target card (*compact* §13.3, L6832–6843):
```
R_pole = D_0(B_4+Z_4) − 3(M+B_2+Z_2)^2 = 0
R_norm = m̂0^2 S_port N_0/D_0 − 54 G c_s^5/(5 a^5 c^5) = 0
P_2 = P_4 = 0
R_tail = Theta_tail (c/c_s)^3 − 1 = 0   (if the tail sector is active)
```
**Tail sector:** `tail_sector_active = false` for this Stage-1 WP1/WP3 pre-registration. `R_tail` is quoted for ledger fidelity but is **not evaluated** in this run and may not be activated after residuals.

5PN actual-branch card (*compact* §13.4, L6848–6856) — WP3 tangent targets:
```
d ln R_tr = d ln R_target = d ln epsilon_eta = 0,   N_Q = 1
```

**PRIMARY observable:** `R_norm = 0` — the GR-quadrupole-normalized outgoing condition. Rationale: it is the only target carrying an *externally benchmarked* constant (the GR quadrupole `54 G c_s^5/(5 a^5 c^5)`, anchored to Peters 1964 / Maggiore in `benchmarks.yaml`), so a target-blind match reproduces *known physics* in the sense of brief §7, not merely internal self-consistency.

**Co-registered secondary targets** (reported with the same error budget, no re-ranking after data): `R_pole = 0`; `P_2 = P_4 = 0`; `chi_Q = 1` (`Delta_Q = 0`). **WP3 tangent targets:** `d ln R_tr = d ln R_target = d ln epsilon_eta = 0`, `N_Q = 1`.

---

## §I. Stability gates (frozen — `ACTUAL_BRANCH_PROTOCOL_V1.md` Frozen-Inputs + Required Packet)

Wall positivity; stable BdG/Krein certificate; mixed-sector positivity; `D0 > 0`; non-dark `N0`. A branch failing any gate is reported as a non-convergence/instability diagnostic, **not** retuned to pass. Every exported packet attaches `parent_action_status`, `boundary_protocol`, `stability_certificate`, `source_hashes`, `freeze_hash`.

---

## §J. Validation requirements (non-negotiable — brief §5; mandatory before any physics claim)

1. **Validate before trusting** — demonstrate the solver on ≥1 known-answer limit (free/linear limit, known GPE vortex/soliton, or standard GPE benchmark) before extracting any unknown observable (plan §7 step 1; `NONLINEAR_PROTOCOL_V2` manufactured-solution + Jacobian checks).
2. **Convergence study** — refine grid (≥3 levels, plan §2.3 / `NONLINEAR_PROTOCOL_V2`); report observed order of convergence per observable; an observable that drifts under refinement is not a measurement.
3. **Boundary control** — absorbing/PML/sponge layers; quantitatively demonstrate boundary reflections sit below the target signal (plan §7 step 5).
4. **Conservation diagnostics** — report mass/charge/energy drift over each run (plan §7 step 6).
5. **Error budget & noise floor** — state the numerical noise floor explicitly; every extracted observable carries a quantitative uncertainty; a match/mismatch is meaningful only relative to that floor.

A result without §J does not count as an answer (brief §6).

---

## §K. Freeze boundary + non-rescue rules (binding after freeze)

Frozen packet flags (`NONLINEAR_PROTOCOL_V2.md` Freeze Boundary, L74–83): `freeze.pre_target_freeze = true`; `freeze.target_blind = true`; `freeze.no_post_residual_refit = true`; `freeze.candidate_freeze_hash` over {predeclared protocol, candidate parameters, solver tolerances, mesh, source revision}; no target residuals / pass flags / target values / post-hoc scores in the frozen packet.

Non-rescue rules (`ACTUAL_BRANCH_PROTOCOL_V1.md` L131–139): do not project a realized `chi_Q ≠ 1` branch back to canonical; do not use `zeta`/support enhancement to explain an orbit-lock miss; do not use scalar `P0` mouth-hammering as direct `P22` moment control; do not mutate source support / boundary class / gauge convention / port normalization / extraction formulas after seeing residuals; do not report an algebraically projected zero-residual packet as a physical target-blind hit. Pre-registered branch only: a miss on this branch falsifies *this* branch (brief §3.4).

---

## §L. Freeze Block (to be completed at user sign-off)

```
frozen: YES
freeze_date: 2026-06-12
freeze_commit:          # SHA of the commit that sets frozen: YES — filled by the immediate follow-up commit
source_revision: 3c1365f6faae7b7c145554ea458339405727c916   # repo HEAD at freeze (ledger + protocol docs unchanged since)
candidate_freeze_hash: df13fb4747dbd4410661245ad3b5938f57a531a5c25d4af974f203dbebc8781f
candidate_freeze_hash_def: "SHA-256 over the text from the '## §C.' header up to (not including) the '## §J.' header (136 lines)"
signed_off_by: user
```

---

## §M. Convention-resolution log (Codex fidelity review — session `019ebd6d`, 2026-06-12)

Read-only Claude+Codex review (transcript `research/pde_ledger/redteam_adversarial/codex_logs/step7_prereg_review.txt`). Overall: **no sign/factor/power error in any coefficient/target formula** (`54/5`, `27/20`, `a^5 c^5`, `P_2/P_4` numerators all verified verbatim against the cited ledger lines). All five flagged choices resolved and folded into §C–§H above:

- **M-1 Gauge — RESOLVED:** H=Z confirmed correct/cleanest for the stationary isotropic branch (parent-clean localized gauge-fixing; clean projected match `S = Z/Z_int`).
- **M-2 Primary observable — RESOLVED:** `R_norm = 0` confirmed primary (only target carrying the externally-benchmarked GR quadrupole constant); §G extended to export `(m̂0, S_port, chi_Q, N_Q)` so `R_norm` is extractable at WP1.
- **M-3 Outgoing-port convention — RESOLVED:** natural point-particle source-map branch (`N_Q = chi_Q^{-1}`) confirmed; no post-solve canonical projection.
- **M-4 Support-placement coordinate — RESOLVED:** the single registered coordinate is `varrho_phys` (closed form folded into §F); `sigma_phys` / ranking-region are derived readouts.
- **M-5 Completeness — RESOLVED:** source-physics scope made explicit (relaxed/open-system §12.12 lanes excluded; block added after §D); `tail_sector_active = false` added to §H. One transcription disambiguation applied: `T^2` → `mathcal_T^2` in §F.

Codex verdict: after these corrections, the record is **fit to freeze**.
