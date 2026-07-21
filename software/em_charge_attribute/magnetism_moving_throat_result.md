# Magnetism as the moving throat — result

## Concise result (body; ≤2 pages)

**Scope and neutral object.** This is a target-blind Tier-A, medium-frame calculation for two slowly moving throats at (Rgg a). The electric input is kept branch-conditional:

\[
U_E=\frac{s_1s_2A_E}{4\pi R},\qquad
F_{E,r}=\frac{s_1s_2A_E}{4\pi R^2},
\]

where (A_E) denotes the selected electric ensemble's full far-field coefficient (its coupled kernel contains (m_{gg}=B_{\rm eff}z_g^2/D)). Its sign and normalization retain the committed `R1_REQUIRED(bc_selection)` status.

**Deliverable #1 — ledger-ready (G0+\delta) row.** G0 contains no (u_T). Importing only the pathA_36 transverse row and supplying the moving-dent coupling gives

\[
S_{T+\mathrm{move}}=\int dt\,d^3x\left[
 \frac{\rho_{\rm br}}2|\dot{\mathbf u}_T|^2
 -\frac{\mu_R}2|\nabla\!\times\!\mathbf u_T|^2
 +q_T\sum_i s_i\eta_a(\mathbf x-\mathbf X_i)
       \mathbf V_i\!\cdot\!\mathbf u_T\right],
\quad \nabla\!\cdot\!\mathbf u_T=0,
\]

\[
c_\gamma^2=\frac{\mu_R}{\rho_{\rm br}},\qquad
q_T=\lambda_T\tau_d.
\]

Here (eta_a) is a normalized finite throat profile and (	au_d) is the active throat's time-arrow (the full time reverse maps a drain to a source). The imported kinetic row and this one finite-profile coupling are the intended amendment, not G0 damage. With (ho_{\rm br},\mu_R>0), the transverse Hessian is stable and no pre-existing G0-v0 row changes: **`internal_inconsistency = none`**.

**Q-CURRENT and field identity.** Translating the actual signed dent,

\[
\sigma_i=s_i\eta_a(\mathbf x-\mathbf X_i(t)),\qquad
\partial_t\sigma_i+\nabla\cdot(\sigma_i\mathbf V_i)=0,
\]

fixes the unique isotropic local flux coefficient to one. Thus

\[
\mathbf I_i=s_i\eta_i\mathbf V_i,qquad
\mathbf J_{T,i}=q_T\mathbf I_i.
\]

No (N_u,a_T,a'_T,a_L,q_A^T,q_L), or pathA_39 source expression is used. The conditional compact limit is therefore **`CONVECTION_LIKE_CONDITIONAL`**: it reduces in tensor form to a profile-smeared (s\mathbf V), but (q_T) is not fixed. The preregistered field is the transverse shear displacement (mathbf u_T), ([mathbf u_T]=L); the curl candidate (mathbf b_T=\nabla\times\mathbf u_T) has unit dimension.

The derived parity census is:

| object | (R_w) | (P_w) | rotations | time reversal |
|---|---:|---:|---|---:|
| (s) | odd | odd | scalar | even |
| (mathbf V) | even | even | polar vector | odd |
| (	au_d), hence (q_T) | even | even | scalar | odd |
| (mathbf J_T=q_Ts\eta\mathbf V) | odd | odd | polar vector | even |
| (mathbf u_T) | odd | odd | polar vector | even |
| (mathbf b_T=\nabla\times\mathbf u_T) | odd | odd | axial vector | even |

Here (P_w) denotes a (w)-type reflection of the transverse coordinate, not ordinary three-dimensional spatial parity (\mathbf x\to-\mathbf x); under that reading, the (s), (\nabla), and (\mathbf b_T) assignments are self-consistent.

A passive time-reversal-even throat would not supply this (O(V)) action row; G0's active drain supplies the required time-arrow basis, while its realized overlap remains R1.

**Q-BOOST (Route A).** Put (mathbf n=(\mathbf X_2-\mathbf X_1)/R), (D_V=\mathbf V_1\cdot\mathbf V_2), and (A_V=(\mathbf V_1\cdot\mathbf n)(\mathbf V_2\cdot\mathbf n)). Independently differentiating the (k^{-4}) radial seed gives the transverse Darwin kernel

\[
I_{ij}(R)=\frac{\delta_{ij}+n_in_j}{8\pi R}.
\]

The Lorentz-completed conditional anchor, for independent ((\mathbf V_1,\mathbf V_2)), is

\[
U_A=\frac{s_1s_2A_E}{4\pi R}
\left[1-\frac{D_V+A_V}{2c_\gamma^2}\right]
+O(v^4/c_\gamma^4),
\]

\[
\mathbf F_{A,2}=\frac{s_1s_2A_E}{8\pi c_\gamma^2R^2}
\left[(\mathbf V_2\!\cdot\!\mathbf n)\mathbf V_1
 +(\mathbf V_1\!\cdot\!\mathbf n)\mathbf V_2
 -(D_V+3A_V)\mathbf n\right],
\quad
F_{A,2,r}=-\frac{s_1s_2A_E(D_V+A_V)}{8\pi c_\gamma^2R^2}.
\]

The (O(1)) and (O(v^2/c_\gamma^2)) orders are explicit; the next uncomputed term is (O(v^4/c_\gamma^4)). For side-by-side motion ((A_V=0)), the velocity term is opposite to the electric radial term for parallel velocities and reverses for antiparallel velocities. A single common side-by-side boost gives (F_r/F_{E,r}=1-v^2/(2c_\gamma^2)+O(v^4)), the required equal-velocity cross-check.

**Q-DIRECT (Route B, blind to Route A).** Route B is produced first, from Q-CURRENT and the transverse Euler equation. A separate tensor ansatz (G_{ij}=(a\delta_{ij}+bn_in_j)/R) gives (b=a) from transversality and (3a+b=2/(4\pi\mu_R)) from the two-polarization trace. Therefore

\[
U_B=-\frac{s_1s_2q_T^2}{8\pi\mu_RR}(D_V+A_V),
\]

\[
\mathbf F_B=\frac{s_1s_2q_T^2}{8\pi\mu_RR^2}
\left[(\mathbf V_2\!\cdot\!\mathbf n)\mathbf V_1
 +(\mathbf V_1\!\cdot\!\mathbf n)\mathbf V_2
 -(D_V+3A_V)\mathbf n\right],
\quad
F_{B,r}=-\frac{s_1s_2q_T^2(D_V+A_V)}{8\pi\mu_RR^2}.
\]

Thus Route B is (O(V_1V_2)), has (U\sim R^{-1}), (F\sim R^{-2}), and reverses between side-by-side parallel and antiparallel motion. Its production dependency record is exactly `{Q_CURRENT.source, pathA36.transverse_EOM, direct_transverse_tensor_ansatz}`. The illicit Route-A-copy mutation multiplies Route A by (q_T^2c_\gamma^2/A_E), leaving the copied value a factor (\mu_R) above Route B; `ROUTE_INDEPENDENCE` detects the illicit dependency tag, not equality of the copied value.

**Q-COMPARE / Q-MAG.** Tensor structure, falloff, and velocity order agree. The analytic coefficient comparison and preferred-frame diagnostics are

\[
r_{BA}=\frac{U_B}{U_{A,2}}
=\frac{q_T^2c_\gamma^2}{\mu_RA_E}
=\frac{q_T^2}{\rho_{\rm br}A_E},
\qquad
\delta_{BA}=r_{BA}-1,
\qquad
r_{\rm cone}=\frac{c_E^2}{c_\gamma^2}
=\frac{c_E^2\rho_{\rm br}}{\mu_R},
\]

\[
\Delta U=U_B-U_{A,2}
=-\frac{s_1s_2A_E}{8\pi c_\gamma^2R}(r_{BA}-1)(D_V+A_V).
\]

Because both the sign of (A_E) and (q_T)'s nonlinear throat normalization are open, the honest comparison fact is **`route_B_R1(nonlinear throat q_T normalization + electric boundary selection)`**, with **`relative_sign_anchor_conditional`**. The dimensionless census contains (q_T/\sqrt{\rho_{\rm br}|A_E|}); it is an unresolved output, not a declared second calibration and not a (q/g) knob. Fact: **`R1(magnitude)`**, with unbounded core-normalization uncertainty at Tier A, profile correction (O(a^2/R^2)), and kinematic remainder (O(v^4/c_\gamma^4)).

**SEALED §4 landing.** The first matching landing is **`R1_REQUIRED(electric_bc_selection)`**. All co-occurring honest blockers are **`R1_REQUIRED(direct_moving_throat)`**, **`R1_REQUIRED(magnitude)`**, and **`R1_REQUIRED(consistency)`**. No terminal consistency, departure, or exclusion label is emitted while either anchor is open. Only after electric branch selection does the direct sign become comparable; only after the throat solve can (r_{BA}=1) or a characterized nonunit delta be decided. On a selected like-repelling electric branch, the computed direct side-by-side parallel term is the downstream like-current attraction; that is not an upstream earn.

**DECIDED:** normalized dent continuity and conditional source identity; all field/source parities and units; the scoped (G0+\delta) row; transverse stability/no G0 damage; general independent-velocity Route A through (O(v^2/c_\gamma^2)); blind Route B's tensor, sign, falloff, and order; (r_{BA}), (Delta U), and (r_{\rm cone}); truth-table totality and action-level hooks. **R1:** electric boundary/variant selection; nonlinear moving-throat (q_T); whether it is tied to (A_E) or is a second factor; (c_E=c_\gamma); (O(v^4))/operator-contamination terms; and the active (F_{\rm flux}) (O(V_1V_2)) contribution plus full-force integrability.

---

## Appendix (evidence; exempt from body cap)

### A. Independent derivations and dimensions

Route A reconstructs

\[
\mathcal F^{-1}\!\left[\frac{P^T_{ij}}{k^2}\right]
=\frac{\delta_{ij}}{4\pi R}
-\mathcal F^{-1}\!\left[\frac{k_ik_j}{k^4}\right]
=\frac{\delta_{ij}+n_in_j}{8\pi R}
\]

by twice differentiating (mathcal F^{-1}[k^{-4}]=-R/(8\pi)). Route B does not call that result: it solves the most general isotropic (R^{-1}) tensor from divergence and trace constraints. Both engines then differentiate their own potentials component-by-component and compare the compact force identities.

Restored units are

| item | dimension |
|---|---:|
| (ho_{\rm br}) | (ML^{-3}) |
| (mu_R) | (EL^{-3}=ML^{-1}T^{-2}) |
| (mathbf u_T, \dot{\mathbf u}_T, \nabla\times\mathbf u_T) | (L, LT^{-1}, 1) |
| (q_T) | (MT^{-1}) |
| (eta_a, A_E) | (L^{-3}, EL) |
| each action density | (EL^{-3}) |
| (U_A,U_B; F_A,F_B) | (E; EL^{-1}) |
| (r_{BA},\delta_{BA},r_{\rm cone}) | (1) |

The amendment activates only the absent transverse DOF and the finite-profile moving `Q_chi V.u_T` component of G0's transverse zero row. Scalar, drain, return, wall, geon, and all other declared-zero entries remain unchanged. The G0 active (F_{\rm flux}) ledger is not silently absorbed into this conservative exchange.

### B. Parallel, antiparallel, and zero controls

For (A_V=0), write (D_V=\pm v_1v_2):

| geometry | Route-A (F_{2,r}/F_{E,r}) | Route-B (F_{B,r}) |
|---|---:|---:|
| parallel | (-v_1v_2/(2c_\gamma^2)) | (-s_1s_2q_T^2v_1v_2/(8\pi\mu_RR^2)) |
| antiparallel | (+v_1v_2/(2c_\gamma^2)) | (+s_1s_2q_T^2v_1v_2/(8\pi\mu_RR^2)) |
| either velocity zero | (0) | (0) |

This table is a neutral sign comparison. The electric branch maps it to charge-language only inside sealed §4.

### C. Exhaustive sealed truth table

The Cartesian domain is

\[
4\ \text{current identities}\times3\ \text{comparison facts}
\times4\ \text{relative-sign facts}\times3\ \text{magnitude facts}
\times2\ \text{tiers}\times2\ \text{electric-anchor states}
\times2\ \text{inconsistency states}=1152.
\]

Current identities are `CONVECTION_LIKE_CONDITIONAL`, `CHARACTERIZED_SOURCE_DEPARTURE`, `NULL_SOURCE`, and `R1_SOURCE_BASIS`. Comparison values are `routes_agree`, `routes_differ`, and `route_B_R1`. Relative-sign values are match, opposite, leading tensor conflict, and anchor-conditional. Magnitude, tier, and electric-anchor values are exactly the directive's enumerated alternatives.

The grouped table below is lossless: both engines enumerate every row in the order above, apply the same first-match precedence through independent implementations, reject every defensive fall-through, and hash the full row stream.

| first matching predicate / landing | cells |
|---|---:|
| `NO_GO(sector)` | 576 |
| `MAGNETISM_LORENTZ_CONSISTENT` | 3 |
| `MAGNETISM_CONSISTENT_FREE_MAGNITUDE(R1)` | 3 |
| `MAGNETISM_DEPARTURE_CHARACTERIZED` | 6 |
| `AMENDMENT_EXCLUDED(wrong_relative_sign)` | 12 |
| `AMENDMENT_EXCLUDED(routes_leading_conflict)` | 12 |
| `R1_REQUIRED(electric_bc_selection)` | 288 |
| `R1_REQUIRED(direct_moving_throat)` | 144 |
| `R1_REQUIRED(magnitude)` | 48 |
| `R1_REQUIRED(consistency)` | 60 |
| defensive `R1_REQUIRED(unclassified)` | 0 |

Truth-table digest:

`983556935e50f12670fef24f17a23c048e295ddf0b6952aa0bd1618e9f179619`

The production tuple is `CONVECTION_LIKE_CONDITIONAL × route_B_R1 × relative_sign_anchor_conditional × R1(magnitude) × tier_A_conditional × R1_REQUIRED(bc_selection) × no inconsistency`. Its first match is `R1_REQUIRED(electric_bc_selection)`; the complete blocker collector emits all four R1 rows stated in the body.

### D. Hooks

- **Characterized Maxwell-parity departure — `B_TIME_REVERSAL_EVEN`:** the candidate magnetic field (\mathbf b_T=\nabla\times\mathbf u_T) is time-reversal even, whereas a physical Maxwell (\mathbf B) is time-reversal odd. This follows because (\mathbf u_T) is a physical brane displacement (T-even) and the source construction includes the active-drain time-arrow (\tau_d); it is internally self-consistent with the derived parity census and consistent with the R1/`DEPARTURE_CHARACTERIZED` landing. It is a concrete way this model's magnetism departs from exact Maxwell, not a defect.
- **Emergent Lorentz invariance:** `UNDETERMINED`. It requires (delta_{BA}=0), (r_{\rm cone}=1), closed higher orders, and no unmatched active-flux term.
- **Preferred-frame signature:** the medium-frame (O(V_1V_2)) tensor is decided; coefficient, cone, (O(v^4)), and pathA_39 operator-parity leakage remain R1.
- **Conservative observable caveat:** G0's non-conservative (F_{\rm flux}) can contribute at (O(V_1V_2)); its value and full-force integrability are explicitly R1.
- **Hierarchy capstone:** flag only. Once both sector couplings come from one throat action, (F_e/F_g) is the held-out dimensionless test; it is not computed here.

### E. Atomic production-path mutation campaigns

Both engines fire exactly one failed assertion for each one-tooth mutation, restore the production value, and finish with every tooth passing:

`SOURCE_NOT_IMPORTED`, `ROUTE_INDEPENDENCE`, and `TARGET_BLINDNESS` demonstrate able-to-fail via flag-flip mutations; the underlying properties are enforced more strongly by live free-symbol disjointness, execution-order independence (Route B is built with `foreign_payload=None` before Route A), and a type-enforced neutral-enum seal, respectively.

`SOURCE_TRANSLATION_CONTINUITY`, `SOURCE_NOT_IMPORTED`, `SOURCE_BASIS`, `PARITY_RW`, `PARITY_PW`, `PARITY_ROTATION`, `PARITY_TIME_REVERSAL`, `FIELD_IDENTITY_UNITS`, `ACTION_KINETIC`, `ACTION_COUPLING`, `ACTION_STABILITY`, `G0_DAMAGE`, `ROUTE_INDEPENDENCE`, `BOOST_PROJECTOR`, `BOOST_GENERAL_VELOCITIES`, `BOOST_NEXT_ORDER`, `BOOST_COMMON_VELOCITY`, `DIRECT_SOURCE`, `DIRECT_PROJECTOR`, `DIRECT_EXCHANGE_SIGN`, `DIRECT_FALLOFF`, `DIRECT_VELOCITY_ORDER`, `COMPARE_COMPUTED`, `DELTA_RATIO`, `CONE_RATIO`, `QMAG_R1`, `UNITS_RESTORED`, `ACTIVE_FLUX_CAVEAT`, `HOOK_LORENTZ`, `LEDGER_READY_ROW`, `TRUTH_TOTALITY`, `TRUTH_PRECEDENCE`, `LANDING_OWNERSHIP`, `TARGET_BLINDNESS`, `DUAL_ENGINE_TERMS`.

SymPy and Mathematica independently construct the source continuity equation, both transverse kernels, both force derivatives, coefficient/delta/cone ratios, units, sealed truth table, and mutation results. Their analytic payloads compare symbolically rather than by text and both normal runs exit 0 with `ENGINE_AGREE=PASS`.
