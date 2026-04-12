# 5PN computational handoff — what still has to be done to finish the program

## 0. Executive status

The reduced theorem side is now essentially closed.

The current chain has already reduced the problem to four finish-line conditions on the actual PDE-selected branch:

1. `d ln R_tr = 0`
2. `d ln R_target = 0`
3. `d ln epsilon_eta = 0`
4. `N_Q = 1`  (equivalently `chi_Q = 1` on the natural source-map branch)

The explicit Family-1 support/source side has already been located on the refreshed exact `Lambda_EM` geometry and sits safely inside the lowest-symmetric-twin regime. So the support/source side is no longer the active bottleneck in the reduced hierarchy.

What is still missing is the **actual PDE-selected orbit-lock point** on the coherent local D/N branch, plus confirmation that the realized outgoing branch is the canonical passive/outgoing one rather than a nearby deformation.

So the remaining computational job is not “more reduced algebra.” It is:

- solve or continue the completed moving-throat branch strongly enough to extract the actual coherent placement state,
- compute the weak-axisymmetric tangent of that branch,
- evaluate the orbit packet and outgoing normalization on the same branch,
- and return a final four-condition verdict.

---

## 1. Non-negotiable computational firewall

These rules are part of the handoff. They should not be changed casually.

### 1.1 Keep the parent bulk PDE fixed

Do **not** retune the parent GNLS / bulk medium to force 5PN.

Keep:

- the `4+1` parent field theory,
- the `P = K rho^5` medium,
- the corrected charge / ontology firewall,
- and the moving-throat geometry lift.

The remaining gap is branch realization, not a bulk-theory refit.

### 1.2 Keep the exact geometry relation

Do **not** revert to the old `L/a = 37/20` shorthand.

Use the carried exact relation

`Lambda_EM = sqrt(2) pi / x_01`

and any downstream geometry derived from it.

### 1.3 Keep support/source and orbit-lock separate

The support ratio `zeta` and the support enhancement factor `S(zeta; epsilon)` affect the isotropic normalization / baseline lane.

They do **not** enter the coherent weak-axisymmetric orbit packet.

So the code path must keep two packets separate:

- **support / normalization packet**
- **orbit / similarity packet**

---

## 2. Exact remaining data packet to extract from the PDE branch

There are really two extraction layers.

### 2.1 Stationary isotropic coherent placement state

The first required packet is the actual stationary coherent-branch state

`(chi0, delta_U, Z_W, epsilon_W, epsilon_eta, Lambda, zeta)`

with derived quantities

`epsilon = epsilon_W * (1 - (2/11) * delta_U/(1 + delta_U))`

`R_tr = (1 + chi0/(1 + delta_U)) / (1 + chi0)`

`R_target = Lambda * (1 - epsilon_eta) * (1 - epsilon)^2 / ( Z_W * (1 + chi0)^2 )`

`M_mix = 8 * Z_W * (1 + chi0)^2 / ( pi^2 * (1 - epsilon_eta) * (1 - epsilon) )`

`M_tr = M_mix * S(zeta; epsilon)`

and, on the outgoing side,

- `chi_Q` and/or
- `N_Q`

with the natural-source-map relation `N_Q = 1/chi_Q` if that branch is used.

### 2.2 Weak-axisymmetric branch tangent / drift packet

The second required packet is the actual first grouped weak-axisymmetric drift of the same branch.

The clean reduced form is

`(d ln chi0, d ln delta_U, d ln Z_W, d ln epsilon_W, d ln epsilon_eta, d ln Lambda)`

Optionally, the microscopic eight-drift packet is also acceptable:

`(lambda_1, c_1, gamma_1, kappa_U, kappa_eta, kappa_W, mu_1, tau_1)`

because that can be pushed through the Stage-320 / Stage-321 / Stage-323 machinery.

### 2.3 Optional but valuable regime-classification packet

For the support regime classifier, it is also useful to extract one of the following:

- `Pi_tr`
- `rho_alpha = Pi_tr / C_mix`
- `xi_phys` together with the selected-branch parameters used in the coherent support notes

This is lower priority than the orbit packet, because the explicit Family-1 support/source side is already known to be non-bottlenecked on the canonical isotropic branch.

---

## 3. Work package sequence

## WP0 — Freeze conventions and branch choice

Before any solve:

1. fix the exact `Lambda_EM` geometry,
2. fix the sign conventions for the outgoing DtN branch,
3. fix the natural source-map convention,
4. fix the coherent local D/N branch definitions,
5. fix whether `chi_Q` or `N_Q` is the primary outgoing observable in the solver output.

Deliverable:

- one short conventions file or JSON manifest used by every downstream script.

## WP1 — Solve / continue the stationary isotropic moving-throat branch

Goal:

Compute the actual stationary isotropic branch state from the completed operator, not from the reduced prototype.

Minimum acceptable output:

`chi0, delta_U, Z_W, epsilon_W, epsilon_eta, Lambda, zeta, chi_Q or N_Q`

Recommended method:

- finite throat in axial coordinate,
- wall profile + support + localized Maxwell/mixed sectors retained,
- continuation from the explicit Family-1 / D/N reduced branch,
- Newton or Newton-Krylov solve on the stationary nonlinear system,
- pseudo-arclength continuation if the branch folds.

Numerical tasks:

1. choose the stationary throat profile family / discretization,
2. solve the isotropic stationary branch,
3. postprocess the solution into the coherent placement variables.

Stop condition:

- stationary residuals small enough that the placement variables are stable under refinement.

## WP2 — Extract isotropic support / normalization packet on the actual branch

Goal:

Evaluate the isotropic branch packet on the realized branch.

Compute:

- `epsilon`
- `R_tr`
- `R_target`
- `M_mix`
- `M_tr`
- `C_mix = 8 Lambda (1 - epsilon) / pi^2`
- `rho_alpha = Pi_tr / C_mix` if `Pi_tr` is available
- `chi_Q` and `N_Q`

Check:

1. support regime:
   - mixed-only,
   - lowest symmetric twin,
   - non-twin.
2. outgoing normalization:
   - `N_Q - 1`, or equivalently `chi_Q - 1`.

This is where the explicit Family-1 branch-location results become a benchmark, not the final answer.

## WP3 — Solve the first weak-axisymmetric tangent problem on that branch

Goal:

Compute the grouped weak-axisymmetric tangent of the **actual** realized stationary branch.

Minimum acceptable output:

`d ln chi0, d ln delta_U, d ln Z_W, d ln epsilon_W, d ln epsilon_eta, d ln Lambda`

Recommended method:

- linearize the completed stationary operator in the grouped real `P2` weak-axisymmetric sector,
- solve the tangent / response problem at zero or low frequency,
- extract the coherent log drifts by projection.

Optional richer output:

- microscopic drift vector `(lambda_1, c_1, gamma_1, kappa_U, kappa_eta, kappa_W, mu_1, tau_1)`.

## WP4 — Evaluate the orbit-lock packet on the realized branch

Goal:

Test whether the actual branch is tangent to the exact similarity orbit.

Compute either the finite quotient packet or its linearization:

### Finite packet

`q_tr, q_nt, q_eta`

or

`Q_tr = ln( C_tr / C_tr,ref )`

`Q_nt = ln( C_nt / C_nt,ref )`

`Q_eta = ln( epsilon_eta / epsilon_eta,ref )`

### Linear packet

`Theta_1, Xi_1, R_1`

with the direct observable formulas

`Theta_1 = d ln R_tr`

`Xi_1 = - d ln R_target - (epsilon_eta / (1 - epsilon_eta)) d ln epsilon_eta`

`R_1 = d ln R_target`

and inverse relations if needed.

Success condition:

`d ln R_tr = 0`

`d ln R_target = 0`

`d ln epsilon_eta = 0`

or equivalently

`q_tr = q_nt = q_eta = 0`.

## WP5 — Final four-condition verdict

Goal:

Issue the final reduced completion verdict on the actual realized branch.

The final conditions are

1. `d ln R_tr = 0`
2. `d ln R_target = 0`
3. `d ln epsilon_eta = 0`
4. `N_Q = 1`  (equivalently `chi_Q = 1` on the natural source-map branch)

Decision tree:

- if 1–3 fail and 4 passes: outgoing normalization is fine, orbit lock fails;
- if 1–3 pass and 4 fails: coherent branch is on-orbit, outgoing normalization fails;
- if support/source fails its regime test but 1–4 pass: check convention mismatch, because that would be inconsistent with the reduced theorem chain;
- only if all four pass do we have the reduced 5PN / 2.5PN / 4PN closure on the realized branch.

---

## 4. Recommended solver / post-processing architecture

### 4.1 Solver side

Use one code path for the stationary isotropic solve and one code path for the weak-axisymmetric tangent solve.

Suggested split:

- `solver_stationary.py` or equivalent:
  solves the isotropic stationary branch.

- `solver_tangent_p2.py` or equivalent:
  linearizes around the stored stationary branch and solves the grouped weak-axisymmetric tangent.

### 4.2 Post-processing side

Re-use the reduced compilers already built in the session chain instead of rewriting formulas by hand.

The most useful existing script anchors are:

- `5pn_stage201_home_stretch_theorem.py`
- `5pn_stage318_actual_branch_orbit_tester_api.py`
- `5pn_stage325_coherent_branch_two_packet_compiler.py`
- `5pn_stage332_actual_branch_theorem_gate.py`
- `5pn_stage339_nontwin_asymmetry_requirement.py`
- `5pn_stage349_actual_four_condition_extractor.py`
- `5pn_stage350_family1_exact_branch_locator.py`

These are best treated as validators / packet compilers once real branch data exist.

### 4.3 Suggested data format

Save one isotropic branch JSON and one tangent JSON.

Recommended isotropic JSON fields:

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

Recommended tangent JSON fields:

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

If the microscopic drift packet is easier to extract, also save it.

---

## 5. Practical numerical warnings

### 5.1 Do not force canonical outgoing normalization

If the realized PDE branch gives `chi_Q != 1`, record it.

Do **not** manually project it back onto the canonical outgoing branch just to make the theorem pass.

### 5.2 Do not mix support and orbit packets

Because `zeta` drops out of the orbit packet exactly, a change in support enhancement should never be used to “explain” an orbit-lock miss.

If the orbit packet misses, it is a real branch-realization miss.

### 5.3 Track refinement sensitivity

For every extracted quantity, compare at least two discretizations / resolutions.

The handoff is not complete unless the final packet is numerically stable under refinement.

### 5.4 Save the full branch state before reduction

Do not only save reduced numbers. Save the actual solved stationary profiles / operator data too, so later work can revisit the reduction if the verdict is surprising.

---

## 6. Definition of done

The computational side is complete when all of the following exist:

1. one converged stationary isotropic branch solution,
2. one converged grouped weak-axisymmetric tangent about that same branch,
3. extracted coherent placement state,
4. extracted orbit packet,
5. extracted outgoing normalization scalar,
6. one final four-condition verdict.

At that point the reduced program is either:

- **closed on the realized PDE branch**, or
- **falsified at one explicit condition**.

Either outcome is useful.
