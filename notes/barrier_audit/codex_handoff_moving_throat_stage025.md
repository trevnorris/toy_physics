# Codex Handoff — Moving-Throat PDE Completion and Same-Charge Stage 025

## 0. Mission

Build a **local-first reduced PDE / operator solver** that extracts the **actual stationary coherent moving-throat branch** and its **first weak-axisymmetric tangent**, then evaluate the final reduced finish-line conditions and the same-charge Stage-025 compiler packet.

This is **not** a request to reopen the reduced theorem side. The reduced algebra is already compressed. The missing object is the **actual PDE-selected orbit-lock / coherent placement point** on the coherent local D/N branch.

The job is therefore:

1. realize the stationary isotropic branch numerically,
2. realize the grouped weak-axisymmetric tangent on that same branch,
3. postprocess into the exact reduced packets already isolated in Stages 174–175 and the same-charge Stage-025 compiler,
4. return a final pass/fail verdict on the remaining conditions.

---

## 1. What is already closed and must **not** be re-derived

### 1.1 Parent theory and ontology

Keep fixed:

- the `4+1` parent field theory,
- the gauged GNLS matter sector,
- the stiff medium `P = K rho^5`,
- the corrected charge / ontology firewall,
- the moving-throat geometry lift,
- the exact geometry relation `Lambda_EM = sqrt(2) pi / x_01` where that relation is already carried forward.

Do **not** retune the parent medium or geometry just to make the final conditions pass.

### 1.2 Packetization theorem is already done

Stages 174–175 already reduced the endgame to two exact finite packets:

- **Packet A** = grouped low-frequency bundle
  - `((D_A0, D_A2, D_A4, N_A0, N_A2, N_A4) for A in {20,21,22}, mhat0)`
- **Packet B** = orbit-lock packet, in any equivalent form
  - `(q_tr, q_nt, q_eta)`
  - or `(m_T, m_K, m_mu)`
  - or `(R_tr_bar, R_nt_bar, R_eta_bar)` depending on notation layer.

The reduced closure criterion is already exact:

- `Delta_branch = 0`
- `Delta_orbit = 0`

Everything else is downstream algebra once those packets are available.

### 1.3 Same-charge front edge is already at Stage 025

Same-charge reduced algebra is closed through:

- **Stage 023**:
  - `rho_alpha = 4/3`
  - `zeta_req = 1/3`
  - `Pi_tr = (4/3) C_mix`
- **Stage 024**:
  - selected twin-support curve
    - `epsilon_* = 1 - 3 varrho / 2`
    - `sigma = 4/(3 varrho) - 2`
    - `0 < varrho < 2/3`
  - thresholds
    - `varrho_(WΛ) = 2(1+beta^2)/(3(2+beta^2))`
    - `varrho_(UΛ) = 2(1+beta^2)/(3(1+beta+beta^2))`
- **Stage 025**:
  - actual coherent branch placement on the selected twin-support curve
  - separate orbit-lock evaluation.

Do **not** spend time re-proving these algebraic facts.

---

## 2. Primary files Codex should load first

Load these first, in this order:

1. `moving_throat_pde_program_compact.md`
2. `5pn_notes_full.md`
3. `barrier_audit_full.md`
4. `moving_throat_pde_stage174_minimal_pde_data_packet.md`
5. `moving_throat_pde_stage175_orbit_quotient_projectors.md`
6. `same_charge_barrier_audit_stage025_actual_twin_support_and_orbit_lock_compiler.md`
7. `same_charge_barrier_audit_stage023_selected_branch_loading_ratio_from_minimal_isotropic_quadrupole_precursor_sympy_audit.py`
8. `same_charge_barrier_audit_stage024_exact_primitive_ranking_on_selected_twin_support_branch_sympy_audit.py`

Secondary / implementation-support references:

- `moving_throat_pde_stage174_minimal_pde_data_packet_sympy_audit.py`
- `moving_throat_pde_stage175_orbit_quotient_projectors_sympy_audit.py`
- `moving_throat_output_full.md`
- `moving_throat_pde_stage003_bdg_coupling.md`
- `moving_throat_pde_stage004_maxwell_mixed_response.md`
- `moving_throat_pde_stage005_grouped_p2_normalization_bridge.md`

---

## 3. Non-negotiable computational firewall

These rules are part of the handoff and should be enforced in code review.

### 3.1 Keep support/source and orbit-lock separate

The support ratio `zeta` and support enhancement factor `S(zeta; epsilon)` belong to the isotropic support / normalization lane.

They do **not** enter the coherent weak-axisymmetric orbit packet.

So the code path must keep two logically separate postprocessing packets:

- **support / normalization packet**
- **orbit / similarity packet**

Do not use support enhancement to explain an orbit miss.

### 3.2 Do not force canonical outgoing normalization

If the realized PDE branch gives `chi_Q != 1` or `N_Q != 1`, record that honestly.

Do **not** project the answer back onto the canonical outgoing branch just to make the theorem pass.

### 3.3 Save full branch state before reduction

Do not save only reduced scalars. Save:

- full stationary profiles,
- operator blocks,
- discrete matrices / vectors needed for reruns,
- the tangent eigenvector / drift data,
- and the reduced packet JSON.

### 3.4 Refinement stability is mandatory

Every reported scalar used in the final verdict must be checked at **at least two resolutions**.

No quantity is “finished” unless it is stable under refinement to an agreed tolerance.

---

## 4. Recommended numerical strategy

### 4.1 Overall approach

Do **not** start with a full nonlinear time-dependent cloud-scale run.

Start with a **harmonic-reduced, finite-throat, stationary + tangent** workflow:

1. stationary isotropic branch solve,
2. linear weak-axisymmetric tangent solve about that converged branch,
3. exact reduced postprocessing.

### 4.2 Local stack

Recommended local prototype stack:

- `scikit-fem` for finite element assembly on the axial finite-throat problem,
- `SciPy` sparse linear algebra / eigensolvers for low-mode extraction,
- plain Python / NumPy for postprocessing and regression tests.

Why:

- easiest local iteration loop,
- simple to instrument and debug,
- good fit for 1D axial finite-throat reduced problems,
- fast enough to validate the reduction chain before scale-up.

### 4.3 Scale-up stack

When the local prototype is stable, the preferred scale-up path is:

- `petsc4py` for scalable PDE linear / nonlinear solves,
- `slepc4py` for large sparse eigenproblems / tangent modes,
- optional `DOLFINx` if a higher-level UFL-based form language becomes useful.

### 4.4 Hardware recommendation

Start on **local CPU**, not GPU.

The first target is sparse stationary solve + sparse tangent/eigen solve on a reduced harmonic formulation. That favors iteration speed, reproducibility, and debugging over GPU throughput.

Escalate to cloud / cluster only after the local solver reproduces the exact reduction identities and the packet outputs are refinement-stable.

---

## 5. Exact quantities the solver must extract

There are two required extraction layers.

### 5.1 Stationary isotropic coherent placement state

Required stationary packet:

`(chi0, delta_U, Z_W, epsilon_W, epsilon_eta, Lambda, zeta, chi_Q or N_Q)`

Derived quantities to compute exactly:

```text
epsilon = epsilon_W * (1 - (2/11) * delta_U / (1 + delta_U))
R_tr = (1 + chi0 / (1 + delta_U)) / (1 + chi0)
R_target = Lambda * (1 - epsilon_eta) * (1 - epsilon)^2 / ( Z_W * (1 + chi0)^2 )
M_mix = 8 * Z_W * (1 + chi0)^2 / ( pi^2 * (1 - epsilon_eta) * (1 - epsilon) )
M_tr = M_mix * S(zeta; epsilon)
C_mix = 8 * Lambda * (1 - epsilon) / pi^2
```

Outgoing-side observable:

- either `chi_Q`
- or `N_Q`
- with natural-source-map relation `N_Q = 1 / chi_Q` when that branch is used.

### 5.2 Weak-axisymmetric tangent / drift packet

Required tangent packet:

```text
(dln_chi0, dln_delta_U, dln_Z_W, dln_epsilon_W, dln_epsilon_eta, dln_Lambda)
```

Optional but useful if easier to extract:

```text
(lambda_1, c_1, gamma_1, kappa_U, kappa_eta, kappa_W, mu_1, tau_1)
```

because Stage 174/175 machinery can convert it.

---

## 6. Exact same-charge Stage-025 compiler formulas

Once the stationary coherent packet is available, compute:

```text
varrho_phys = (2/3) * (1 - epsilon)
epsilon_*_phys = epsilon
sigma_phys = 2 * epsilon / (1 - epsilon)
```

Primitive threshold rewrite in the realized variable `epsilon`:

```text
epsilon_WLambda = 1 / (2 + beta^2)
epsilon_ULambda = beta / (1 + beta + beta^2)
```

Ranking decision:

- Region I if `epsilon > epsilon_WLambda`
- Region II if `epsilon_ULambda < epsilon < epsilon_WLambda`
- Region III if `0 < epsilon < epsilon_ULambda`

Support-lane classifier remains:

```text
Pi_tr <= C_mix                      -> mixed-only enough
C_mix < Pi_tr <= 2 * C_mix          -> lowest symmetric twin enough
Pi_tr > 2 * C_mix                   -> non-twin asymmetry required
```

But on the selected same-charge branch:

```text
Pi_tr = (4/3) * C_mix
```

so the realized support slice is automatically in the **lowest symmetric twin** window.

Important notation firewall:

- `zeta_req = 1/3` is the Stage-023 reduced demand-side loading ratio.
- physical coherent support ratio `zeta` is a different object.
- do **not** conflate them.

---

## 7. Exact orbit-lock and outgoing finish-line conditions

### 7.1 Direct orbit observables

Compute:

```text
dln_epsilon = dln_epsilon_W - [2 * delta_U / ((1 + delta_U) * (11 + 9 * delta_U))] * dln_delta_U
```

```text
dln_R_tr = - [chi0 * delta_U / ((1 + chi0) * (1 + delta_U) * (1 + chi0 + delta_U))]
           * ((1 + delta_U) * dln_chi0 + (1 + chi0) * dln_delta_U)
```

```text
dln_R_target = dln_Lambda - dln_Z_W
                - [2 * chi0 / (1 + chi0)] * dln_chi0
                - [epsilon_eta / (1 - epsilon_eta)] * dln_epsilon_eta
                - [2 * epsilon / (1 - epsilon)] * dln_epsilon
```

### 7.2 Finish-line conditions

The orbit-lock condition is exactly:

```text
dln_R_tr = 0
dln_R_target = 0
dln_epsilon_eta = 0
```

Outgoing finish-line is separate:

```text
N_Q = 1
```

or equivalently:

```text
chi_Q = 1
```

on the natural source-map branch.

### 7.3 Direct defect packet

For reporting, also compute:

```text
Theta_1 = dln_R_tr
Xi_1 = -dln_R_target - [epsilon_eta / (1 - epsilon_eta)] * dln_epsilon_eta
R_1 = dln_R_target
```

---

## 8. Packet A and Packet B postprocessing requirements

### 8.1 Packet A

Extract the grouped low-frequency coefficients for each grouped real `P2` lane `A in {20,21,22}`:

```text
D_A(omega) = D_A0 + D_A2 * omega^2 + D_A4 * omega^4 + O(omega^6)
N_A(omega) = N_A0 + N_A2 * omega^2 + N_A4 * omega^4 + O(omega^6)
```

and source-map factor `mhat0`.

Compile:

```text
u2^(A) = - D_A2 / D_A0
u4^(A) = (D_A2^2 - D_A0 * D_A4) / D_A0^2
P0^(A) = N_A0 / D_A0
P2^(A) = (D_A0 * N_A2 - 2 * D_A2 * N_A0) / D_A0^2
P4^(A) = (D_A0^2 * N_A4 - 2 * D_A0 * (D_A2 * N_A2 + D_A4 * N_A0) + 3 * D_A2^2 * N_A0) / D_A0^3
```

Grouped trace/anomaly decomposition for any grouped triple `(x20, x21, x22)`:

```text
xbar = (x20 + 2*x21 + 2*x22) / 5
a_x = (2*x20 - x21 - x22) / 10
b_x = (x21 - x22) / 2
```

Final branch residual packet:

```text
Delta_pole = u4_bar - 4 * u2_bar^2
Delta_norm = mhat0^2 * P0_bar - 54 * G * c_s^5 / (5 * a^5 * c^5)
Delta_branch = (a2, b2, a4, b4, aP0, bP0, Delta_pole, Delta_norm)
```

### 8.2 Packet B

If the microscopic eight-drift packet is available, use Stage 175 exactly.

Finite drift basis:

```text
Delta x = (Delta_lambda, Delta_c, Delta_gamma, Delta_U,
           Delta_Keta, Delta_W, Delta_mu, Delta_T)^T
```

Monomial-drift map:

```text
q = (q_tr, q_nt, q_eta)^T = M_* Delta x
```

Use Stage-175 quotient/orbit projector machinery:

- dependent triple `(Delta_T, Delta_Keta, Delta_mu)`
- exact quotient section `S_(T,Keta,mu)`
- projectors
  - `Q_quot = S M_*`
  - `O_orb = I - Q_quot`

Acceptance interpretation:

```text
Q_quot Delta x = 0  <=>  M_* Delta x = 0
```

That is the exact microscopic orbit-lock condition.

If you do **not** have the microscopic packet, the direct observable orbit packet
`(dln_R_tr, dln_R_target, dln_epsilon_eta)` is sufficient for the final finish-line verdict.

---

## 9. Implementation work packages

## WP0 — Freeze conventions and branch choice

Before any solve:

1. freeze `Lambda_EM` geometry relation,
2. freeze sign conventions for outgoing DtN branch,
3. freeze natural source-map convention,
4. freeze coherent local D/N branch definitions,
5. freeze whether `chi_Q` or `N_Q` is the primary outgoing observable in raw solver output,
6. define output JSON schemas.

Deliverable:

- `config/conventions.yaml`
- one short `README_conventions.md`

## WP1 — Build the local reduced solver skeleton

Repository shape recommendation:

```text
moving_throat_solver/
  README.md
  config/
    conventions.yaml
    solver_defaults.yaml
  src/
    mesh_axial.py
    geometry_wall.py
    stationary_branch.py
    tangent_branch.py
    bdg_block.py
    maxwell_mixed_block.py
    packetA.py
    packetB.py
    stage025_same_charge.py
    io_json.py
  tests/
    test_stage023.py
    test_stage024.py
    test_stage174.py
    test_stage175.py
    test_stage025.py
  outputs/
    stationary/
    tangent/
    packets/
```

Minimum objective:

- finite-throat axial mesh / discretization,
- geometry-only wall baseline,
- clean operator assembly plumbing,
- unit tests for exact algebraic compilers.

## WP2 — Conservative support stack

Add, in order:

1. geometry-only wall operator,
2. stable BdG support contribution via exact Schur-complement reduction,
3. localized Maxwell/mixed conservative block,
4. grouped real `P2` lane handling.

Minimum output of WP2:

- `D_A0, D_A2, D_A4`
- `N_A0, N_A2, N_A4`
- `mhat0`
- grouped isotropy diagnostics on the isotropic reference branch.

## WP3 — Stationary isotropic coherent branch solve

Goal:

- continue / solve the actual stationary isotropic coherent branch from a good reduced seed.

Recommended method:

- Newton or Newton-Krylov,
- pseudo-arclength continuation if branch folds,
- keep finite-throat D/N structure explicit,
- save full solution objects.

Required output:

```json
{
  "chi0": 0.0,
  "delta_U": 0.0,
  "Z_W": 0.0,
  "epsilon_W": 0.0,
  "epsilon_eta": 0.0,
  "Lambda": 0.0,
  "zeta": 0.0,
  "chi_Q": 0.0,
  "N_Q": 0.0,
  "Pi_tr": null,
  "mhat0": 1.0
}
```

## WP4 — Isotropic support / normalization postprocessing

From the stationary state, compute:

- `epsilon`
- `R_tr`
- `R_target`
- `M_mix`
- `M_tr`
- `C_mix`
- support regime classification
- `chi_Q - 1` and/or `N_Q - 1`

Also evaluate isotropic one-pole and outgoing-normalization identities whenever Packet A is available:

```text
D0 * (B4 + Z4) = 3 * (M + B2 + Z2)^2
mhat0^2 * N0 / D0 = 54 * G * c_s^5 / (5 * a^5 * c^5)
```

## WP5 — First weak-axisymmetric tangent solve

Goal:

- solve the grouped weak-axisymmetric tangent about the converged stationary branch.

Required output:

```json
{
  "dln_chi0": 0.0,
  "dln_delta_U": 0.0,
  "dln_Z_W": 0.0,
  "dln_epsilon_W": 0.0,
  "dln_epsilon_eta": 0.0,
  "dln_Lambda": 0.0
}
```

If microscopic drift is available, save it too.

## WP6 — Orbit packet and Stage-025 same-charge postprocessing

Compute:

- `dln_R_tr`
- `dln_R_target`
- `dln_epsilon_eta`
- `Theta_1, Xi_1, R_1`
- `varrho_phys`
- `sigma_phys`
- support ranking region
- `N_Q - 1`

Return the **smallest Stage-025 output packet**:

```text
(epsilon,
 varrho_phys,
 sigma_phys,
 ranking_region,
 R_tr,
 R_target,
 epsilon_eta,
 dln_R_tr,
 dln_R_target,
 dln_epsilon_eta,
 N_Q_minus_1)
```

## WP7 — Final verdict compiler

Definition of done:

1. one converged stationary isotropic branch solution,
2. one converged grouped weak-axisymmetric tangent on the same branch,
3. Packet A extracted,
4. Packet B or direct orbit packet extracted,
5. same-charge Stage-025 packet extracted,
6. final finish-line verdict written down honestly.

Final verdict fields:

```json
{
  "stationary_converged": true,
  "tangent_converged": true,
  "delta_branch": {
    "a2": 0.0,
    "b2": 0.0,
    "a4": 0.0,
    "b4": 0.0,
    "aP0": 0.0,
    "bP0": 0.0,
    "Delta_pole": 0.0,
    "Delta_norm": 0.0
  },
  "orbit_packet": {
    "dln_R_tr": 0.0,
    "dln_R_target": 0.0,
    "dln_epsilon_eta": 0.0,
    "N_Q_minus_1": 0.0
  },
  "same_charge_packet": {
    "epsilon": 0.0,
    "varrho_phys": 0.0,
    "sigma_phys": 0.0,
    "ranking_region": "",
    "R_tr": 0.0,
    "R_target": 0.0,
    "epsilon_eta": 0.0
  },
  "pass": false,
  "failure_mode": ""
}
```

---

## 10. Exact regression tests that must exist before any large run

### 10.1 Algebraic regression

Must pass locally:

- Stage-023 identities:
  - `rho_alpha = 4/3`
  - `zeta_req = 1/3`
  - `Pi_tr = (4/3) C_mix`
- Stage-024 identities:
  - twin-support curve
  - threshold formulas and ordering
- Stage-174 compiler formulas for Packet A and Packet B
- Stage-175 projector identities and orbit-lock equivalence
- Stage-025 support/orbit compiler formulas.

### 10.2 Numerical sanity checks

Before any cloud run, local solver must show:

1. isotropic grouped-lane degeneracy on the isotropic reference branch,
2. packet outputs stable across refinement,
3. stationary and tangent results saved reproducibly,
4. no manual enforcement of `chi_Q = 1` or `N_Q = 1`.

---

## 11. When to escalate from local machine to cluster

Stay local until all of the following are true:

1. exact regression suite passes,
2. stationary isotropic branch converges reliably,
3. tangent solve converges reliably,
4. packet extraction is refinement-stable,
5. final verdict script runs end-to-end.

Only then escalate if one of these is true:

- resolution required exceeds local RAM / time budget,
- continuation path becomes too stiff for local iterative solves,
- grouped tangent solve becomes large enough that SLEPc-class eigensolvers are justified.

---

## 12. Practical success criterion

The program is computationally complete when Codex can produce, from the actual branch rather than reduced notes alone:

1. one stationary coherent placement state,
2. one weak-axisymmetric tangent on that same branch,
3. the Stage-174/175 reduced packet verdict,
4. the Stage-025 same-charge packet,
5. and an honest final statement:
   - either the realized branch closes all four finish-line conditions,
   - or it fails at one explicit named condition.

Either result is useful.

---

## 13. Ready-to-paste Codex starter prompt

```text
You are continuing the moving-throat PDE completion and same-charge Stage-025 program.

Read these files first:
1. moving_throat_pde_program_compact.md
2. 5pn_notes_full.md
3. barrier_audit_full.md
4. moving_throat_pde_stage174_minimal_pde_data_packet.md
5. moving_throat_pde_stage175_orbit_quotient_projectors.md
6. same_charge_barrier_audit_stage025_actual_twin_support_and_orbit_lock_compiler.md
7. same_charge_barrier_audit_stage023_selected_branch_loading_ratio_from_minimal_isotropic_quadrupole_precursor_sympy_audit.py
8. same_charge_barrier_audit_stage024_exact_primitive_ranking_on_selected_twin_support_branch_sympy_audit.py

Non-negotiable rules:
- Do not rederive already-closed reduced algebra.
- Do not retune the parent PDE or EOS.
- Keep support/normalization and orbit packets separate.
- Do not force chi_Q = 1 or N_Q = 1.
- Save full branch state before reduction.
- Require refinement stability before calling anything done.

Primary goal:
Build a local-first reduced solver that extracts the actual stationary coherent branch and its first weak-axisymmetric tangent, then compute:
- dln R_tr
- dln R_target
- dln epsilon_eta
- N_Q - 1
and the Stage-025 same-charge packet:
- epsilon
- varrho_phys
- sigma_phys
- ranking region
- R_tr
- R_target
- epsilon_eta

Work packages:
- WP0 conventions and JSON schemas
- WP1 solver skeleton
- WP2 conservative support stack
- WP3 stationary isotropic branch solve
- WP4 isotropic packet extraction
- WP5 weak-axisymmetric tangent solve
- WP6 Stage-025 postprocessing
- WP7 final verdict compiler

Required exact formulas and JSON schemas are in codex_handoff_moving_throat_stage025.md.

Start by creating the repository skeleton, conventions file, output schemas, and the exact algebraic regression tests for Stage 023, Stage 024, Stage 174, Stage 175, and Stage 025 before attempting the stationary solve.
```

