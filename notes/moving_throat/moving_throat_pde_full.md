
=== moving_throat_pde_stage001_geometry_lift.md ===

# Moving-Throat PDE — Stage 1 Geometry Lift and Linearized PDE Skeleton

## 0. What this note does

This is the first actual derivation move after the roadmap and Phase-0 theorem-target sheet.
It does **not** try to solve the full nonlinear moving-throat problem.
It fixes the first missing ingredient needed before any honest PDE work can proceed:

1. an explicit **distributed throat geometry variable** replacing the current finite-dimensional `a(t), L(t)` closure,
2. a **minimal linearized coupled bulk/interface system** for matter, gauge, and geometry,
3. and a concrete **mode decomposition** in which the scalar lane, grouped real `P2` lane, and geometry lane become explicit.

The key decision is to use a **hybrid level-set / shape-field representation**:

- represent the finite throat surface by a level-set
  
  `Sigma(X,t)=0`,

- but parameterize that level-set by a brane-spherical mouth shape field
  
  `r = R(Omega,w,t)`,

  where `r = sqrt(x^2+y^2+z^2)` is the brane radius and `w` is the bulk axial direction.

This is the smallest upgrade that keeps the 4D ontology intact while making the grouped real `P2` modes explicit.

---

## 1. Why this is the right geometry lift

The existing parent theory already gives:

- exact bulk GNLS matter,
- exact localized `4+1` Maxwell,
- and a geometry sector carried only by `a(t), L(t)` through `V_conf(X; a, L)`.

That is enough for the reduced conservative hierarchy, but not enough for a genuine moving-throat PDE.

The missing PDE must do three jobs simultaneously:

1. make the current `a,L` variables emerge as collective moments,
2. expose the conservative grouped real `P2` support/geometry channels explicitly,
3. and provide the right mouth/worldtube response operator whose low-frequency data match the 2PN / 3PN / 2.5PN / 4PN ledgers.

A pure collective-coordinate closure cannot do that.
A pure abstract level-set `Sigma(X,t)` can do it, but hides the multipole structure.
A pure moving-wall field `R(theta,phi,w,t)` makes the multipoles visible, but is less clean as a bulk coupling variable.

So the right compromise is:

`Sigma(x,y,z,w,t) = r - R(Omega,w,t)`

with

- `r = sqrt(x^2+y^2+z^2)`,
- `Omega = x/r in S^2`,
- the throat surface given by `Sigma=0`,
- the exterior by `Sigma>0`,
- the interior/support region by `Sigma<0`.

This makes the mouth `S^2` harmonics explicit while preserving a level-set coupling to the bulk fields.

---

## 2. Reference stationary throat background

Take a stationary finite-throat reference surface

`Sigma_0(X) = r - R_0(w) = 0`,

with:

- mouth at `w=0`,
- finite throat depth `0 <= w <= L_0`,
- reference mouth radius
  
  `a_0 = R_0(0)`,

- and closure at the far end through either
  
  `R_0(L_0)=0`
  
  or an equivalent regular bottom condition.

This is the minimal finite-throat geometry compatible with the already preferred finite branch and the D/N-style internal support picture.

The current collective variables are recovered as monopole moments:

`a(t) = (1/4pi) int_{S^2} R(Omega,0,t) dOmega`,

and, if `W_b(Omega,t)` is the bottom location defined by `R(Omega,W_b(Omega,t),t)=0`,

`L(t) = (1/4pi) int_{S^2} W_b(Omega,t) dOmega`.

So the old `a(t),L(t)` survive, but only as the lowest geometry moments.

---

## 3. Promoted confinement coupling

The old confinement potential is `V_conf(X; a, L)`.
Promote it to a true moving-surface coupling by writing

`V_conf(X; Sigma) = V_wall( Sigma(X,t) / ell_c )`,

where `ell_c` is a wall-thickness / smoothing scale and `V_wall` is a steep profile.

Then the stationary background is

`V_0(X) = V_wall( Sigma_0(X) / ell_c )`.

Linearizing gives

`delta V_conf = (V_wall'(Sigma_0/ell_c)/ell_c) delta Sigma`.

Since

`Sigma = r - R(Omega,w,t)`,

one has simply

`delta Sigma = - delta R(Omega,w,t)`,

so

`delta V_conf = - (V_wall'(Sigma_0/ell_c)/ell_c) delta R`.

This is the cleanest way to couple the moving throat directly to the bulk fields.

---

## 4. Mode decomposition of the throat geometry

Write the geometry perturbation as

`R(Omega,w,t) = R_0(w) + eta(Omega,w,t)`.

Expand `eta` in real spherical harmonics on the brane mouth sphere:

`eta(Omega,w,t) = eta_0(w,t) Y_00(Omega)`
`                + sum_{m in P2(real)} q_{2m}(w,t) Y_{2m}^{real}(Omega)`
`                + eta_{>=3}(Omega,w,t)`.

The grouped real `P2` set is

- `20`,
- `21c`, `21s`,
- `22c`, `22s`.

This is the first place where the conservative grouped real `P2` lane becomes a literal geometry/support mode rather than a symbolic residual.

Normalization note: with the standard normalized real harmonic

`Y_00 = 1 / (2 sqrt(pi))`,

all grouped real `P2` harmonics have zero sphere average, while

`(1/4pi) int_{S^2} Y_00 dOmega = 1/(2 sqrt(pi))`.

So the physical mouth-average monopole shift `delta a` is related to the normalized monopole coefficient by

`q_00 = 2 sqrt(pi) delta a`.

A useful further split is:

`eta_0(w,t) = alpha_a(w) delta a(t) + alpha_L(w) delta L(t) + g(w,t)`,

where:

- `delta a` is the mouth-monopole radius shift,
- `delta L` is the axial-depth / bottom shift,
- `g(w,t)` is the remaining axisymmetric geometry lane orthogonal to the `a` and `L` moments.

So the geometry sector separates into:

1. `l=0` collective throat motion (`delta a`, `delta L`),
2. `l=0` residual geometry lane `g`,
3. grouped real `l=2` lanes `q_{2m}`,
4. higher `l>=3` lanes.

At the level of the current PN problem, the first serious nontrivial conservative payload is exactly the `l=2` bundle plus the axisymmetric geometry lane.

---

## 5. Minimal distributed geometry action (new ansatz)

This part is **new** and is not yet frozen by the current papers. It is the smallest passive local action I would try first.

Let `eta` be the radial wall displacement field above. Then take the quadratic geometry action

`S_eta^(2) = (1/2) int dt dw dOmega sqrt(gamma_0)`
`           [ mu_eta(w) (partial_t eta)^2`
`             - T_w(w) (partial_w eta)^2`
`             - T_Omega(w) eta (-Delta_{S^2}) eta`
`             - K_eta(w) eta^2 ]`.

Here:

- `mu_eta(w)` is the effective wall inertia density,
- `T_w(w)` is axial wall stiffness,
- `T_Omega(w)` is angular stiffness on the mouth sphere,
- `K_eta(w)` is a local restoring potential,
- `gamma_0` is the metric determinant on the background throat surface.

This is the minimal distributed upgrade that:

- reduces to a genuine PDE in `w`,
- makes the `S^2` multipole structure explicit,
- and gives separate scalar and quadrupole restoring lanes automatically.

Because

`-Delta_{S^2} Y_{lm} = l(l+1) Y_{lm}`,

this action yields the modal operator

`mu_eta q_{lm,tt} - partial_w( T_w partial_w q_{lm} ) + [ K_eta + l(l+1) T_Omega ] q_{lm}`
`= source_{lm}^{(psi,A)}`.

So the scalar lane `l=0` and grouped real `P2` lane `l=2` are already split before any further closure.

This is exactly the structural feature the current hierarchy has been asking for.

---

## 6. Stationary background fields

Choose a stationary finite-throat background

`psi(X,t) = exp(-i mu_0 t / hbar) psi_0(X)`,

`A_M(X,t) = A_{M0}(X)`,

`R(Omega,w,t) = R_0(w)`.

The background solves the exact parent bulk equations with the promoted confinement potential `V_0(X)`.

At this stage, mixed channels must stay alive:

- `A_w`,
- `J^w`,
- `F_{mu w}`,
- `E_w`,
- `C_a = F_{aw}`.

They may decouple later on a far-field brane branch, but they must **not** be removed from the PDE that is supposed to generate the response data.

---

## 7. Linearized matter sector

Write

`psi = exp(-i mu_0 t/hbar) [ psi_0 + delta psi ]`,

and define

`delta rho = psi_0^* delta psi + psi_0 delta psi^*`.

The linearized GNLS problem is naturally of Bogoliubov-de Gennes type:

`i hbar partial_t [delta psi, delta psi^*]^T`
`= L_BdG [delta psi, delta psi^*]^T`
`  + C_A[delta A_M]`
`  + C_eta[eta]`.

A convenient explicit form is

`H_0 = -(hbar^2/2m) D_{0i} D_0^i + V_0 + h(rho_0) - mu_0`,

and then schematically

`i hbar partial_t delta psi`
`= H_0 delta psi + h'(rho_0) delta rho psi_0 + delta V_conf psi_0 + C_A^(+)[delta A]`,

`-i hbar partial_t delta psi^*`
`= H_0^* delta psi^* + h'(rho_0) delta rho psi_0^* + delta V_conf psi_0^* + C_A^(-)[delta A]`.

Since

`delta V_conf = - (V_wall'(Sigma_0/ell_c)/ell_c) eta`,

the wall-displacement field enters the BdG system directly as the first moving-throat coupling.

---

## 8. Linearized gauge sector

Write

`A_M = A_{M0} + delta A_M`.

Linearizing the exact localized Maxwell equation gives

`partial_M [ Z(w) delta F^{MN} ] + (1/xi) partial^N (partial . delta A)`
`= mu_0 delta J^N`.

The linearized current is the exact first variation of the minimally coupled matter current. In schematic form,

`delta J^0 = q_* delta rho`,

`delta J^i = q_* delta j^i`,

with

`delta j^i = (hbar/m) Im( psi_0^* D_0^i delta psi + delta psi^* D_0^i psi_0 )`
`           - (q_*/m) rho_0 delta A^i`
`           + ...`.

The important point is not the final algebraic form but the channel content:

- `delta A_w` survives,
- `delta J^w` survives,
- `delta F_{mu w}` survives.

So the PDE can actually test whether the mixed channels that were suppressed in the brane limit are part of the microscopic response law.

---

## 9. Linearized geometry equation

Varying the new quadratic geometry action and coupling to the matter/gauge sectors gives the geometry PDE

`mu_eta partial_t^2 eta - partial_w( T_w partial_w eta ) - T_Omega Delta_{S^2} eta + K_eta eta`
`= S_eta^(psi)[delta psi,delta psi^*] + S_eta^(A)[delta A_M] + f_ext`.

At the linearized level, the first source term from the confinement coupling is generated by

`delta V_conf = - (V_wall'(Sigma_0/ell_c)/ell_c) eta`,

so the matter back-reaction contributes through the linearized generalized force

`S_eta^(psi) ~ - (V_wall'(Sigma_0/ell_c)/ell_c) delta rho + ...`.

Projecting onto spherical harmonics gives mode equations

`mu_eta q_{lm,tt} - partial_w( T_w partial_w q_{lm} )`
`+ [ K_eta + l(l+1) T_Omega ] q_{lm}`
`= S_{lm}^{(psi,A)} + f_{lm}^{ext}`.

This is the first usable linearized moving-throat PDE family.

---

## 10. The first mode sectors that matter

For the present program, the important sectors are:

### `l = 0` axisymmetric sector

This contains:

- `delta a(t)`,
- `delta L(t)`,
- the residual geometry lane `g(w,t)`.

This is the lane that should eventually reproduce the symbolic geometry completion / scalar-response structures left open in the lower-order hierarchy.

### grouped real `l = 2` sector

This contains

- `q_20(w,t)`,
- `q_21c(w,t)`, `q_21s(w,t)`,
- `q_22c(w,t)`, `q_22s(w,t)`.

On an exactly isotropic reference throat, these are degenerate at the geometry-only level.
Any splitting among the grouped channels then measures anisotropy induced by matter support, gauge support, or wall structure.

This is exactly the conservative information the 3PN grouped real `P2` closure needs.

### higher `l >= 3` sector

These modes should be kept in the full PDE, but they are not needed in the first reduced theorem target.

---

## 11. Mouth/worldtube response operator in the new variables

At the mouth `w = 0`, choose a real harmonic port basis

- `P_00 = Y_00`,
- `P_20`, `P_21c`, `P_21s`, `P_22c`, `P_22s`.

Then define effort variables by projecting the drive on the mouth/worldtube,
for example:

- scalar enthalpy drive `u_0(omega)`,
- grouped real `P2` drives `u_{2m}(omega)`,
- optional geometry drive `u_g(omega)`.

Define measured fluxes as the matching projected normal fluxes / generalized tractions,
for example:

- matter flux through the mouth,
- leakage flux in `w`,
- wall traction `T_w partial_w q_{lm}` at the mouth,
- grouped support amplitudes.

Then the response operator is

`j_A(omega) = sum_B Z_eff,AB(omega) u_B(omega)`.

For an isotropic background this should block-diagonalize into

- scalar/geometry block,
- grouped real `P2` block,
- higher-`l` blocks.

The grouped `P2` block is the one that eventually feeds

- `u2^(20), u2^(21), u2^(22)`,
- isotropy defects,
- and the canonical pair `(Kbar_0, Kbar_2)` on the passive/outgoing quadrupole branch.

---

## 12. How the open observables are supposed to emerge

This geometry lift gives a direct interpretation of the still-open response data.

### Dynamic pole scales

The mode-resolved DtN / response operator in each sector should have poles.
Those are the natural microscopic origin of the open dynamic scales

- `Omega_{1 perp}`,
- `Omega_10`,
- `Omega_0`,
- `Omega_20`,
- `Omega_21`,
- `Omega_22`,
- `Omega_g`.

### Conservative grouped real `P2` data

The low-frequency even expansion of the grouped `P2` response block yields the coefficients

- `u2^(20)`,
- `u2^(21)`,
- `u2^(22)`.

### Geometry completion

The axisymmetric residual geometry lane `g(w,t)` is the most natural microscopic carrier for the separate geometry completion that survives in the 3PN conservative split.

### Quadrupole normalization

On the passive/outgoing isotropic quadrupole branch, the same mode-resolved response should determine the canonical invariant pair

- `Kbar_0`,
- `Kbar_2`,

and hence the effective quadrupole normalization `gamma_quad^eff`.

---

## 13. Why I think this is the right first move

This lift has four advantages.

### Advantage 1 — it matches the actual ontology

The parent theory already has exact bulk matter and gauge fields, and it already treats the defect as a finite throat/cavity object rather than a point source.
This lift preserves that ontology instead of collapsing it too early.

### Advantage 2 — it makes the grouped `P2` lane literal

The grouped real quadrupole channels are no longer symbolic coefficients. They are actual `l=2` wall/support modes on the mouth sphere.

### Advantage 3 — it keeps `a` and `L` honest

The old `a(t),L(t)` remain present, but only as the lowest collective moments of a genuine distributed geometry field.
That is the right relationship between the old reduced hierarchy and the new PDE.

### Advantage 4 — it is the smallest passive distributed geometry action that can work

The quadratic wall action above is the smallest honest PDE upgrade I know that automatically separates:

- scalar/geometry,
- grouped real `P2`,
- and higher multipoles.

That is precisely the structure the current hierarchy has been trying to expose.

---

## 14. Immediate next derivation step

The next real calculation should be:

1. choose one explicit stationary reference profile `R_0(w)`,
2. project the linearized geometry PDE onto `l=0` and grouped real `l=2`,
3. compute the geometry-only mouth DtN operator as the baseline passive branch,
4. then add the BdG + Maxwell couplings and track how the pole locations and low-frequency coefficients shift.

That is the first point where the moving-throat PDE stops being a slogan and starts producing the observables the current hierarchy still leaves open.

=== moving_throat_pde_stage002_breathing_reduction.md ===

# Moving-Throat PDE — Phase 2A Breathing Reduction Back to the Old `(a,L)` Closure

## Purpose
This note is the first continuation after the Phase-1 scaffold and the finite-throat D/N benchmark.
Its job is to show that the new distributed wall field is not floating free from the old hierarchy.
At the conservative linear level, the axisymmetric part of the distributed wall theory reduces back to the same kind of matrix system that the parent 4D files already use for the effective throat variables `a(t)` and `L(t)`.

That reduction matters because it tells us we have not broken continuity with the existing program.
The distributed wall field is a **lift** of the old closure, not a replacement that forgets it.

---

## 1. Axisymmetric wall ansatz

Use the normalized real harmonic
\[
Y_{00}(\Omega)=\frac{1}{2\sqrt\pi},
\qquad
\int_{S^2}Y_{00}^2\,d\Omega=1,
\qquad
\frac{1}{4\pi}\int_{S^2}Y_{00}\,d\Omega=\frac{1}{2\sqrt\pi}.
\]
So if the physical mouth-average radius shift is `delta a(t)`, the corresponding monopole harmonic coefficient is
\[
q_{00}(0,t)=2\sqrt\pi\,\delta a(t).
\]
This is the normalization bridge between the distributed wall field and the old collective coordinate.

Take the axisymmetric wall perturbation to be the two-mode truncation
\[
\eta_{0}(w,t)=2\sqrt\pi\Big[\alpha_a(w)\,\delta a(t)+\alpha_L(w)\,\delta L(t)\Big],
\]
where:
- `alpha_a(w)` is the profile that changes the mouth radius,
- `alpha_L(w)` is the profile that changes the axial throat extent,
- higher axisymmetric wall modes are deferred.

This is the minimal truncation that can recover the old `(a,L)` closure.

---

## 2. Insert the two-mode ansatz into the quadratic wall action

On the axisymmetric branch, the quadratic wall action is
\[
S^{(2)}_{\eta,0}
=
\frac12\int dt\,dw\,d\Omega\,
\Big[
\mu_\eta(w)(\partial_t\eta_0)^2
-
T_w(w)(\partial_w\eta_0)^2
-
K_0(w)\eta_0^2
\Big],
\]
with
\[
K_0(w)=K_\eta(w).
\]
Using the two-mode ansatz and the normalized harmonic integral
\[
\int_{S^2}Y_{00}^2\,d\Omega=1,
\]
the reduced Lagrangian becomes
\[
L_{\rm red}^{(0)}
=
\frac12 M_{AB}\,\dot Q^A\dot Q^B
-
\frac12 K_{AB}\,Q^A Q^B,
\qquad
Q^A=(\delta a,\delta L),
\]
with the effective matrices
\[
M_{AB}=
4\pi\int dw\,\mu_\eta(w)\,\alpha_A(w)\alpha_B(w),
\]
\[
K_{AB}=
4\pi\int dw\,
\Big[
T_w(w)\,\alpha_A'(w)\alpha_B'(w)
+
K_0(w)\,\alpha_A(w)\alpha_B(w)
\Big].
\]
Indices `A,B` run over `a,L`.

So the distributed wall does exactly what it should do: it produces effective inertia and stiffness matrices by overlap integrals of the wall profiles.

---

## 3. Euler–Lagrange reduction

The reduced equations of motion are
\[
M_{AB}\,\ddot Q^B+K_{AB}\,Q^B=0.
\]
This is the conservative linearized version of the old geometry equations. In the parent 4D files the geometry sector is written schematically as
\[
M_{AB}\,\ddot Q^B+\Gamma_{AB}\,\dot Q^B=-\frac{\partial H_{\rm tot}}{\partial Q^A}.
\]
So the new distributed-wall reduction reproduces the old closure at the expected level:
- the conservative part comes directly from the wall action,
- any damping matrix `Gamma_AB` is an effective/open-system completion that would enter after coupling to matter/gauge/exterior channels.

This is exactly the relationship we wanted. The old `(a,L)` geometry equations are the lowest-mode truncation of the new distributed wall theory.

---

## 4. What happens to the grouped real `P2` sector at the same level

Now take one grouped real quadrupole component,
\[
\eta_{2m}(\Omega,w,t)=\beta_2(w)\,q_{2m}(t)\,Y_{2m}^{\rm real}(\Omega),
\]
with
\[
-\Delta_{S^2}Y_{2m}^{\rm real}=6Y_{2m}^{\rm real}.
\]
Using the same quadratic wall action, the reduced one-mode Lagrangian is
\[
L_{2m}
=
\frac12 M_2\,\dot q_{2m}^2
-
\frac12 K_2\,q_{2m}^2,
\]
with
\[
M_2=
\int dw\,\mu_\eta(w)\,\beta_2(w)^2,
\]
\[
K_2=
\int dw\,
\Big[
T_w(w)\,\beta_2'(w)^2
+
\bigl(K_\eta(w)+6T_\Omega(w)\bigr)\beta_2(w)^2
\Big].
\]
So before any symmetry breaking or matter/gauge coupling, every real `P2` component has the same uncoupled operator
\[
M_2\,\ddot q_{2m}+K_2\,q_{2m}=0.
\]
That is the microscopic reason the grouped real `P2` block is degenerate on the isotropic reference throat.

---

## 5. Why this reduction matters for the roadmap

This reduction gives three concrete answers.

First, it proves that the distributed wall lift is compatible with the old geometry sector. The old `a,L` equations are not being discarded; they are being reinterpreted as the lowest axisymmetric collective truncation.

Second, it shows that the grouped real `P2` bundle is on exactly the same footing as the old scalar geometry variables. It is not an artificial add-on. It is the next harmonic family of the same wall PDE.

Third, it tells us what the next coupled calculation has to do. Once the BdG matter sector, the localized Maxwell sector, and the mixed channels are turned back on, they must:
- renormalize the effective matrices,
- split or preserve the `P2` degeneracy,
- shift the pole data,
- and produce the passive/outgoing odd parts that the uncoupled wall theory cannot generate on its own.

---

## 6. Script-backed status

The accompanying SymPy audit verifies the concrete algebraic claims used here:
- the real-harmonic normalization and zero-average facts,
- the axisymmetric mouth-average extraction rule,
- the reduced two-mode Euler–Lagrange matrix form,
- and the one-mode grouped-`P2` reduction.

Supporting file:
- `moving_throat_pde_master_sympy_audit.py`

=== moving_throat_pde_stage003_bdg_coupling.md ===


# Moving-Throat PDE — Stage 3: Minimal BdG–Wall Coupling and First Pole Shifts

## Purpose

This note is the first genuinely coupled continuation of the moving-throat PDE program.

The earlier stages did three things already:

1. lifted the old finite-dimensional geometry variables `(a,L)` into a distributed wall field,
2. established the geometry-only modal structure and finite-throat D/N benchmark,
3. and showed that the axisymmetric lowest-mode truncation reduces back to an `(a,L)`-type matrix system.

The next missing step is to turn on the **matter support dynamics** in a form that is still honest enough to be useful and still simple enough to audit symbolically.

So the goal here is not yet “solve the full coupled GNLS + localized Maxwell + moving-wall PDE.”
The goal is narrower:

- reduce the linearized GNLS/BdG sector to stable normal coordinates in the `l=0` and grouped real `l=2` branches,
- couple those matter modes to the wall modes through the already-derived confinement variation,
- integrate the matter modes out exactly at the reduced level,
- and identify the first pole shifts and the first conservative grouped-`P2` isotropy/splitting diagnostics.

This is the first point where the moving-throat PDE begins to produce the sort of quantities the 2PN/3PN/2.5PN stack still leaves open.

---

## 1. Parent-theory anchor and why this is the right next step

The current stack already gives a fully action-based bulk parent with three sectors:

- gauged 4D GNLS matter,
- localized `4+1` Maxwell,
- and explicit geometry variables `(a,L)` / confinement structure.

That parent is strong enough to support reduced-sector work, but it is still **not** a fully solved moving-throat PDE theorem. The unresolved sectors still include the completed moving-throat dynamics and the dynamic derivation of the pole/compliance data. The 2PN and 3PN summaries make the bottleneck especially sharp: the remaining conservative-to-radiative bridge is no longer algebraic bookkeeping, but the microscopic derivation of the grouped real `P2` constitutive law and the geometry completion, followed by the narrow passive/outgoing quadrupole-normalization gap. See `atom_work.md`, `4d_2pn.tex`, `4d_3pn_summary.md`, and `4d_1pn_full_summary.md`. 

What this means operationally is:

- the geometry-only wall PDE is not enough,
- the next live conservative data must come from **coupling wall motion to actual matter support modes**,
- and the first useful reduced question is therefore:

  > how do stable BdG modes renormalize the wall response, shift the pole data, and preserve or split the grouped real `P2` degeneracy?

That is exactly what this stage computes.

---

## 2. Linearized coupling from the promoted confinement law

From the previous stage we already have the promoted confinement lift
\[
V_{\rm conf}(X;\Sigma)=V_{\rm wall}\!\left(\frac{\Sigma(X,t)}{\ell_c}\right),
\qquad
\Sigma=r-R(\Omega,w,t),
\]
and its linearization around the stationary reference throat
\[
\delta V_{\rm conf}
=
-\frac{V_{\rm wall}'(\Sigma_0/\ell_c)}{\ell_c}\,\eta(\Omega,w,t).
\]

So the wall displacement enters the linearized GNLS/BdG problem as a direct source.
Once the linearized matter sector is decomposed into normal modes, the wall–matter coupling is bilinear.

On an isotropic reference throat, the angular dependence is diagonal in the spherical harmonic labels because the coupling is scalar in the angles and the grouped real harmonics are orthonormal. That means:

- the axisymmetric wall sector couples only to axisymmetric matter modes,
- each grouped real `l=2` wall channel couples only inside the real `l=2` matter bundle,
- and no linear `l=0 <-> l=2` mixing appears on the isotropic background.

This selection rule is one of the important reasons the present reduction is clean.

---

## 3. Stable-mode reduction of the linearized matter sector

At the exact PDE level the linearized GNLS problem is BdG-like and first order in time.
For the present step, we pass to the **stable normal-mode reduction** of that BdG system.

That is a controlled reduced-sector move, not a full theorem of the unsolved PDE.
But it is the natural one if we want the first conservative pole shifts and low-frequency response data.

### 3.1 Axisymmetric scalar sector

Let
\[
Q^A(t)=\bigl(\delta a(t),\delta L(t)\bigr),
\qquad A\in\{a,L\},
\]
be the two lowest axisymmetric wall coordinates, and let
\[
X_\alpha(t)
\]
be stable scalar BdG normal coordinates with frequencies
\[
\varpi_{0\alpha}>0.
\]

The reduced quadratic Lagrangian is
\[
L^{(0)}_{\rm red}
=
\frac12 \dot Q^T M_0 \dot Q
-\frac12 Q^T K_0 Q
+\frac12\sum_\alpha \left(\dot X_\alpha^2-\varpi_{0\alpha}^2 X_\alpha^2\right)
+\sum_{A,\alpha} C_{A\alpha} Q^A X_\alpha .
\]

Here:

- `M_0` and `K_0` are the geometry-only axisymmetric inertia/stiffness matrices from the previous reduction,
- `C_{A\alpha}` are wall–matter overlap integrals generated by the confinement variation and the BdG mode profiles.

At this level, the scalar matter sector does not yet radiate or dissipate. It produces a conservative self-energy.

### 3.2 Grouped real `P2` sector

For each grouped real quadrupole channel
\[
A\in\{20,21,22\},
\]
let \(q_A(t)\) be the wall amplitude and let \(X_{A\alpha}(t)\) be stable BdG quadrupole coordinates with frequencies \(\varpi_{A\alpha}\).

Then the reduced quadratic Lagrangian is
\[
L^{(2)}_{\rm red}
=
\sum_A
\left[
\frac12 M_A \dot q_A^2
-\frac12 K_A q_A^2
+\frac12\sum_\alpha\left(\dot X_{A\alpha}^2-\varpi_{A\alpha}^2 X_{A\alpha}^2\right)
+\sum_\alpha g_{A\alpha} q_A X_{A\alpha}
\right].
\]

On the isotropic reference branch one expects
\[
M_{20}=M_{21}=M_{22},
\qquad
K_{20}=K_{21}=K_{22},
\]
and, if the matter support is also isotropic,
\[
g_{20,\alpha}=g_{21,\alpha}=g_{22,\alpha},
\qquad
\varpi_{20,\alpha}=\varpi_{21,\alpha}=\varpi_{22,\alpha}.
\]

That is the precise reduced meaning of “grouped real `P2` degeneracy survives matter coupling on the isotropic branch.”

---

## 4. Exact elimination of the stable BdG modes

The reduced matter modes can be integrated out exactly.

## 4.1 Axisymmetric matrix kernel

In frequency space, the scalar matter coordinates satisfy
\[
(\varpi_{0\alpha}^2-\omega^2)X_\alpha
=
C_{A\alpha} Q^A.
\]
So
\[
X_\alpha
=
\frac{C_{A\alpha}Q^A}{\varpi_{0\alpha}^2-\omega^2}.
\]

Substituting back into the wall equations gives the exact effective axisymmetric kernel
\[
\boxed{
D^{\rm eff}_0(\omega)
=
K_0-\omega^2 M_0
-
C(\Omega_0^2-\omega^2 I)^{-1}C^T,
}
\]
where \(C\) is the `2 x N` coupling matrix and
\[
\Omega_0^2=\mathrm{diag}(\varpi_{0\alpha}^2).
\]

This is the first explicit reduced formula in which matter support renormalizes the geometry lane.

Its low-frequency expansion is
\[
D^{\rm eff}_0(\omega)
=
K^{\rm eff}_0
-\omega^2 M^{\rm eff}_0
-\omega^4 N^{\rm eff}_0
+O(\omega^6),
\]
with
\[
\boxed{
K^{\rm eff}_0
=
K_0-C\Omega_0^{-2}C^T,
}
\]
\[
\boxed{
M^{\rm eff}_0
=
M_0+C\Omega_0^{-4}C^T,
}
\]
\[
\boxed{
N^{\rm eff}_0
=
C\Omega_0^{-6}C^T.
}
\]

So the conservative scalar BdG support does three things immediately:

1. it lowers the effective static stiffness,
2. it increases the effective inertial coefficient,
3. and it generates the next higher even coefficient automatically.

This is already enough to see how the axisymmetric geometry-completion lane can become a real microscopic dynamical object rather than a symbolic leftover.

## 4.2 Grouped real `P2` kernels

Channel by channel,
\[
(\varpi_{A\alpha}^2-\omega^2)X_{A\alpha}
=
g_{A\alpha} q_A,
\]
so
\[
X_{A\alpha}
=
\frac{g_{A\alpha}}{\varpi_{A\alpha}^2-\omega^2}q_A.
\]

Substituting back gives the exact channelwise grouped-`P2` kernels
\[
\boxed{
D_A^{\rm eff}(\omega)
=
K_A-M_A\omega^2
-\sum_\alpha \frac{g_{A\alpha}^2}{\varpi_{A\alpha}^2-\omega^2},
\qquad
A\in\{20,21,22\}.
}
\]

Expanding at low frequency,
\[
D_A^{\rm eff}(\omega)
=
d_{0A}^{\rm eff}
+d_{2A}^{\rm eff}\omega^2
+d_{4A}^{\rm eff}\omega^4
+O(\omega^6),
\]
with
\[
d_{0A}^{\rm eff}
=
K_A-\sum_\alpha \frac{g_{A\alpha}^2}{\varpi_{A\alpha}^2},
\]
\[
d_{2A}^{\rm eff}
=
-\left(
M_A+\sum_\alpha \frac{g_{A\alpha}^2}{\varpi_{A\alpha}^4}
\right),
\]
\[
d_{4A}^{\rm eff}
=
-\sum_\alpha \frac{g_{A\alpha}^2}{\varpi_{A\alpha}^6}.
\]

So stable BdG support automatically feeds the entire even low-frequency ladder of each grouped real `P2` channel.

That is exactly the kind of mechanism the conservative 3PN and 2.5PN bridge has been waiting for.

---

## 5. First exact pole shift formula

The simplest nontrivial example is one wall mode \(q\) coupled to one stable BdG mode \(X\):
\[
L
=
\frac12 M \dot q^2
-\frac12 K q^2
+\frac12 \dot X^2
-\frac12 \varpi^2 X^2
+g q X.
\]

The exact dispersion relation is
\[
(K-M\omega^2)(\varpi^2-\omega^2)-g^2=0.
\]

Writing the uncoupled wall frequency as
\[
\Omega_\eta^2=\frac{K}{M},
\]
the exact poles are
\[
\boxed{
\omega_\pm^2
=
\frac{\Omega_\eta^2+\varpi^2
\pm
\sqrt{(\Omega_\eta^2-\varpi^2)^2+4g^2/M}}{2}.
}
\]

If the matter mode is above the wall mode, \(\varpi^2>\Omega_\eta^2\), and the coupling is weak, then the wall-like pole shifts by
\[
\boxed{
\delta\Omega_\eta^2
=
-\frac{g^2/M}{\varpi^2-\Omega_\eta^2}
+O(g^4),
}
\]
while the matter-like pole shifts upward by
\[
\delta\varpi^2
=
+\frac{g^2/M}{\varpi^2-\Omega_\eta^2}
+O(g^4).
\]

This is the first explicit closed formula for the way matter support shifts a wall pole.

In the language of the earlier hierarchy, this is the reduced-sector prototype for how the open pole data are supposed to arise once the wall is no longer treated in isolation.

---

## 6. Grouped real `P2` isotropy diagnostics

The 3PN grouped-`P2` summary already uses the standard grouped trace/anomaly variables
\[
\bar u_2=\frac{u_2^{(20)}+2u_2^{(21)}+2u_2^{(22)}}{5},
\qquad
a_2=\frac{2u_2^{(20)}-u_2^{(21)}-u_2^{(22)}}{10},
\qquad
b_2=\frac{u_2^{(21)}-u_2^{(22)}}{2}.
\]

At the present reduced level, the same algebra applies to any fixed normalization of the channelwise low-frequency coefficients. Using the `d_{2A}^{eff}` defined above,
\[
\bar d_2=\frac{d_{220}^{\rm eff}+2d_{221}^{\rm eff}+2d_{222}^{\rm eff}}{5},
\]
\[
a_2=\frac{2d_{220}^{\rm eff}-d_{221}^{\rm eff}-d_{222}^{\rm eff}}{10},
\qquad
b_2=\frac{d_{221}^{\rm eff}-d_{222}^{\rm eff}}{2}.
\]

So the first conservative isotropy theorem on the matter-coupled branch is immediate:

- if the wall parameters, matter frequencies, and couplings are channel-independent,
  \[
  K_{20}=K_{21}=K_{22},\quad
  M_{20}=M_{21}=M_{22},\quad
  g_{20,\alpha}=g_{21,\alpha}=g_{22,\alpha},\quad
  \varpi_{20,\alpha}=\varpi_{21,\alpha}=\varpi_{22,\alpha},
  \]
  then
  \[
  \boxed{a_2=b_2=0.}
  \]

Conversely, any nonzero \(a_2\) or \(b_2\) must come from anisotropy in one or more of:

- the wall baseline,
- the matter support frequencies,
- or the wall–matter overlap integrals.

This is the first microscopic conservative criterion for the grouped real `P2` constitutive law.

---

## 7. What this stage achieves physically

This stage does **not** yet solve the full moving-throat PDE.
But it does something important and concrete.

### 7.1 It gives the first microscopic origin of the open pole data

The open 2PN dynamic pole scales were always telling us that the zero-frequency conservative ledger was not the full story.
The present reduced kernels show how those pole shifts arise structurally once the wall is coupled to stable matter support modes.

### 7.2 It makes the grouped real `P2` bundle dynamical rather than symbolic

In the 3PN summary the grouped real `P2` sector is already the live conservative payload that still needs a microscopic derivation.
The present coupled kernels show exactly how a grouped real `P2` bundle acquires low-frequency coefficients and isotropy/splitting data from actual support dynamics.

### 7.3 It explains why the axisymmetric geometry lane can carry the geometry completion

The axisymmetric scalar kernel is matrix-valued in `(a,L)` even before any further complication is added.
Once matter support renormalizes that matrix, a separate geometry lane is no longer mysterious.
It is the natural leftover after the low-dimensional wall sector is coupled back to scalar support modes.

---

## 8. What is still missing

This note is intentionally only the **first matter-coupled step**.

The following are still outside its scope:

1. the localized Maxwell sector on the same footing,
2. the mixed channels
   \[
   A_w,\qquad J^w,\qquad F_{\mu w},
   \]
   which remain part of the microscopic ontology outside the strict brane limit,
3. the passive/outgoing completion of the coupled wall–matter system,
4. and the final normalization step that converts the surviving quadrupole branch into the 2.5PN/4PN universal coefficient.

So the right reading is:

- geometry-only wall theory gave the baseline,
- this stage adds the first conservative matter self-energy,
- the next step is to add the localized Maxwell + mixed sector on the same reduced footing.

---

## 9. SymPy-backed status

The accompanying script verifies:

- the axisymmetric Euler–Lagrange reduction with wall coordinates coupled to stable scalar BdG modes,
- exact elimination of the matter modes and the resulting effective matrix kernel,
- the low-frequency renormalization formulas for the scalar/geometry lane,
- the exact two-pole spectrum and perturbative wall-pole shift in the one-mode case,
- the channelwise grouped-`P2` self-energies,
- the grouped trace/anomaly formulas,
- preservation of `a_2=b_2=0` on the isotropic matter-coupled branch,
- and the harmonic selection rule enforcing the `l=0` / grouped `l=2` block structure.

Supporting files:
- `moving_throat_pde_stage3_bdg_sympy_audit.py`
- `moving_throat_pde_stage3_bdg_sympy_output.txt`

---

## 10. Immediate next derivation step

The next clean move is now well defined:

1. keep the present BdG–wall reduction,
2. add the localized Maxwell and mixed-channel normal coordinates on the same footing,
3. derive the enlarged effective kernel,
4. and identify which part of that enlarged system can generate the passive/outgoing branch needed for the final quadrupole normalization.

That is the first point where the moving-throat PDE program can start touching the actual remaining 2.5PN/4PN bridge instead of only its conservative front end.

=== moving_throat_pde_stage004_maxwell_mixed_response.md ===

# Moving-Throat PDE — Stage 4: Localized Maxwell + Mixed-Sector Reduction and the First Honest Outgoing Quadrupole Bridge

## Purpose

This note is the next reduced-sector derivation after the first conservative BdG–wall coupling stage.

The previous step did something important but incomplete:

- it coupled the distributed wall to stable BdG support modes,
- it generated exact conservative self-energies and pole shifts,
- and it gave the first microscopic grouped real `P2` low-frequency formulas on the isotropic branch.

What it still did **not** do was introduce the sector that can honestly carry a passive/outgoing channel.
That sector is the localized Maxwell + mixed-channel block,

a) because the parent 4D theory already contains a genuine localized `4+1` Maxwell sector,

b) because the mixed variables
`A_w`, `F_{\mu w}`, and `J^w`
are suppressed only in the strict far-field brane limit rather than removed from the microscopic ontology,

and c) because the 2.5PN and 4PN packages already tell us that the surviving universal radiative route is the passive/outgoing quadrupole branch rather than another conservative closure.

So the right next question is:

> how do localized Maxwell modes and the mixed `A_w/F_{\mu w}/J^w` sector feed into the wall response, and how does an outgoing `l=2` branch transfer its `i \omega^5` signature to the wall/worldtube operator?

That is the whole purpose of this stage.

---

## 1. Parent-theory anchor

The exact parent theory already gives a localized Maxwell action and, from it, a bulk field equation of the form

`∂_M(Z(w) F^{MN}) + (1/ξ) ∂^N(∂·A) = μ0 J^N`,

with the mixed channels retained microscopically even when the controlled far-field brane Maxwell reduction later suppresses them. The 4D plasma extension then makes the mixed fields explicit:

`E_w = F_{w0} = -∂_t A_w - ∂_w A_0`,

`C_a = F_{aw} = ∂_a A_w - ∂_w A_a`,

and emphasizes that `A_w` is a genuine dynamical degree of freedom whose mixed forces can drive `j^w` and therefore brane leakage. This is exactly the sector we need if the moving-throat PDE is ever going to produce an honest passive/outgoing channel rather than only more conservative support self-energy.

The earlier 2.5PN audit then narrows the radiative goal further:

- the universal point-particle route is the orbital/worldtube STF quadrupole branch,
- the compact outgoing `l=2` fingerprint starts at `+ i ω^5`,
- and the remaining gap is the passive/outgoing normalization on the actual moving-throat branch.

So the localized Maxwell + mixed block is not an optional embellishment. It is the first microscopically available place where that branch can actually live.

---

## 2. Exact mixed-field gauge invariants

The first thing to check is that the mixed fields are real observables, not gauge artifacts.

Using the scalar-potential convention already implicit in the project summaries,

`A_0 -> A_0 - ∂_t χ`,

`A_a -> A_a + ∂_a χ`,

`A_w -> A_w + ∂_w χ`,

the mixed fields

`E_w = -∂_t A_w - ∂_w A_0`,

`C_a = ∂_a A_w - ∂_w A_a`,

are exactly gauge invariant.

That matters because it means the reduced mixed-sector coordinates are not bookkeeping artifacts; they are reduced representatives of honest parent-theory observables.

---

## 3. Minimal reduced Maxwell/mixed model in one harmonic lane

The smallest useful reduced model is one wall amplitude `Q_l`, one brane-like localized Maxwell coordinate `A_l`, and one mixed `A_w/F_{\mu w}/J^w`-active coordinate `W_l`.

The reduced quadratic Lagrangian is

`L_l`
`= (1/2) M_l Qdot_l^2 - (1/2) K_l Q_l^2`
`  + (1/2) Adot_l^2 - (1/2) Ω_{A,l}^2 A_l^2`
`  + (1/2) Wdot_l^2 - (1/2) Ω_{W,l}^2 W_l^2`
`  + R_l A_l W_l`
`  + g_{A,l} Q_l A_l + g_{W,l} Q_l W_l`.

Interpretation:

- `A_l` is a localized brane-like Maxwell support coordinate,
- `W_l` is the mixed-sector coordinate that is active in the `A_w/F_{\mu w}/J^w` block,
- `R_l` mixes the two internal gauge sectors,
- `g_{A,l}` and `g_{W,l}` are wall couplings generated by the moving-throat geometry.

This is already enough to show how the wall can inherit both conservative gauge self-energy and, after dressing the mixed block, a passive/outgoing odd term.

---

## 4. Exact conservative Maxwell/mixed self-energy

With the frequency convention `exp(-i ω t)`, define

`A_l(ω) = Ω_{A,l}^2 - ω^2`,

`W_l(ω) = Ω_{W,l}^2 - ω^2`,

`Δ_l(ω) = A_l(ω) W_l(ω) - R_l^2`.

Eliminating the two internal gauge coordinates gives the exact conservative self-energy

`Σ_l^(EM+mix,cons)(ω)`
`= [ g_{A,l}^2 W_l(ω) + 2 g_{A,l} g_{W,l} R_l + g_{W,l}^2 A_l(ω) ] / Δ_l(ω)`.

So the wall operator becomes

`D_l^(cons)(ω) = K_l - M_l ω^2 - Σ_l^(EM+mix,cons)(ω)`

before any BdG or outgoing dressing is added.

This is the Maxwell/mixed analog of the Stage-3 BdG self-energy.
It is purely even/conservative at this stage.

### 4.1 Low-frequency coefficients

Write

`D0_l = Ω_{A,l}^2 Ω_{W,l}^2 - R_l^2`,

`S2_l = Ω_{A,l}^2 + Ω_{W,l}^2`,

`N0_l = g_{A,l}^2 Ω_{W,l}^2 + 2 g_{A,l} g_{W,l} R_l + g_{W,l}^2 Ω_{A,l}^2`,

`G2_l = g_{A,l}^2 + g_{W,l}^2`.

Then the conservative self-energy expands as

`Σ_l^(EM+mix,cons)(ω)`
`= z_{l,0} + z_{l,2} ω^2 + z_{l,4} ω^4 + O(ω^6)`

with

`z_{l,0} = N0_l / D0_l`,

`z_{l,2} = (N0_l S2_l - G2_l D0_l) / D0_l^2`,

`z_{l,4} = [ N0_l (S2_l^2 - D0_l) - S2_l G2_l D0_l ] / D0_l^3`.

So the localized Maxwell + mixed block already contributes to the conservative low-frequency moments that the grouped real `P2` and geometry lanes are asking the moving-throat PDE to supply.

---

## 5. Dressing the mixed block by an outgoing port

Now comes the new ingredient.

The mixed coordinate is the natural place to attach the exterior passive/outgoing channel.
So replace

`W_l(ω) -> W_l(ω) - Π_l^(out)(ω)`.

The full self-energy becomes

`Σ_l^(full)(ω)`
`= [ g_{A,l}^2 ( W_l(ω) - Π_l^(out)(ω) ) + 2 g_{A,l} g_{W,l} R_l + g_{W,l}^2 A_l(ω) ]`
`  / [ A_l(ω) ( W_l(ω) - Π_l^(out)(ω) ) - R_l^2 ]`.

Expanding to first order in the outgoing port gives

`Σ_l^(full)(ω) = Σ_l^(cons)(ω) + Π_l^(out)(ω) N_l(ω) + O((Π_l^(out))^2)`

with the exact transfer factor

`N_l(ω) = [ A_l(ω) g_{W,l} + R_l g_{A,l} ]^2 / [ A_l(ω) W_l(ω) - R_l^2 ]^2`.

At zero frequency,

`N_l(0) = [ Ω_{A,l}^2 g_{W,l} + R_l g_{A,l} ]^2 / [ Ω_{A,l}^2 Ω_{W,l}^2 - R_l^2 ]^2 >= 0`.

This is one of the main results of the stage.

It says that if the outgoing branch is passive, then the wall inherits that odd part with a manifestly nonnegative transfer coefficient.
The moving-throat PDE does not have to guess the sign at this reduced level; the sign is carried by the outgoing port, and the internal Maxwell/mixed block only multiplies it by a positive conservative factor.

### 5.1 Wall-operator sign convention

Because the wall operator is written as

`D_l = K_l - M_l ω^2 - Σ_l`,

a positive-imaginary outgoing port

`Π_l^(out)(ω) = + i Γ_l^(port) ω^{2l+1} + ...`

appears in the wall operator as

`δ D_l^(odd)(ω) = - i N_l(0) Γ_l^(port) ω^{2l+1} + ...`.

That is the correct passive damping sign in the wall-operator convention used here.
If instead one rewrites the result in normalized response/admittance language, the same branch appears with the positive `+ i ω^{2l+1}` sign used in the 2.5PN theorem audit.

---

## 6. Compact outgoing `l=2` fingerprint

The 2.5PN audit already isolated the compact outgoing `l=2` branch as the universal quadrupole route.
The exact reduced fingerprint can be recovered from the outgoing spherical Hankel solution.

For the normalized outgoing `l=2` response one finds

`Yhat_2^(out)(ω)`
`= 1 + a^2 ω^2/(9 c_s^2) + 4 a^4 ω^4/(81 c_s^4) + i a^5 ω^5/(27 c_s^5) + O(ω^6)`.

So the outgoing port coefficient is

`Γ_2^(port) = a^5 / (27 c_s^5)`.

This is exactly the coefficient already isolated in the 2.5PN referee audit.

---

## 7. First honest wall-level odd quadrupole coefficient

Combining the last two sections gives the first wall-level odd quadrupole coefficient generated by the localized Maxwell + mixed sector:

`δ D_2^(odd)(ω)`
`= - i N_2(0) [ a^5 / (27 c_s^5) ] ω^5 + O(ω^7)`.

So the first honest mixed-sector quadrupole bridge is already visible.

The structure is exactly what we hoped for:

1. the compact outgoing branch supplies the universal `ω^5` fingerprint,
2. the internal mixed-sector dynamics supplies a positive conservative transfer factor `N_2(0)`,
3. and the remaining normalization problem is no longer diffuse.

It is now the concrete task of computing the static conservative transfer factor on the true moving-throat branch.

That is the Stage-4 version of the 2.5PN normalization gap.

---

## 8. Scalar compatibility check

One of the main 2.5PN worries was that a naive scalar outlet could generate a dangerous `i ω` monopole term before the quadrupole branch ever matters.

The present mixed-sector reduction gives a clean compatibility check.
If the scalar/breathing lane enters the port-active mixed sector only through a **derivative** coupling,
for example

`g_{W,0}(ω) = η ω`,

and there is no direct non-derivative wall coupling into the same port-active combination, then

`N_0(ω) ~ const * ω^2`.

So even if the scalar outgoing port law itself begins as

`Π_0^(out)(ω) = i γ_1 ω + ...`,

the wall-level correction becomes

`δ D_0^(odd)(ω) ~ i const * ω^3 + ...`.

So the mixed-sector reduction is **compatible** with the scalar rescue pattern isolated in the 2.5PN audit: the dangerous linear odd scalar term need not reappear at wall level if the scalar outlet is derivative-coupled.

This is not yet a full scalar no-go theorem, but it is the first explicit reduced model showing how the rescue can happen inside the same Maxwell/mixed-sector language rather than by hand.

---

## 9. What this stage accomplishes physically

This stage is the first place where the moving-throat PDE program honestly connects three things at once:

1. **conservative internal support data** from the wall + BdG + localized Maxwell/mixed block,
2. **mixed-sector transport structure** carried by the `A_w/F_{\mu w}/J^w`-active coordinate,
3. and the **passive/outgoing quadrupole fingerprint** needed by the 2.5PN and 4PN bridges.

Before this stage, the program had

- conservative wall poles,
- conservative BdG self-energy,
- and a symbolic statement that the outgoing branch had to appear somewhere.

After this stage, we now have

- an exact conservative Maxwell/mixed self-energy,
- an exact outgoing-transfer factor,
- an explicit wall-level odd quadrupole coefficient,
- and a concrete scalar derivative-coupling mechanism that delays the scalar odd term to `ω^3`.

That is a real narrowing of the theorem gap.

---

## 10. What is still missing

This is still **not** the full moving-throat PDE theorem.
Several things remain open.

### 10.1 The actual moving-throat branch data

The reduced coefficient

`N_2(0) = [ Ω_A^2 g_W + R g_A ]^2 / [ Ω_A^2 Ω_W^2 - R^2 ]^2`

is only the one-mode reduced prototype.
The completed PDE must determine its actual multi-mode generalization on the true moving-throat branch.

### 10.2 Grouped real `P2` channelization with the mixed sector included

The present stage was carried out lane by lane.
To feed the 3PN and 2.5PN grouped real `P2` programs all the way, the same Maxwell/mixed reduction has to be performed in the full grouped real `20/21/22` bundle and then checked for isotropy defects.

### 10.3 Full normalization against the universal coefficient

The 2.5PN theorem package says the final target is still the universal normalization condition.
The present stage isolates the product that must be computed, but does not yet prove that it equals the GR value on the actual moving-throat branch.

---

## 11. Immediate next step

The right next derivation is now much sharper.

Do **not** jump yet to the full nonlinear PDE.
Instead:

1. keep the wall + BdG + Maxwell/mixed conservative block derived here,
2. lift it from the one-lane prototype to the grouped real `P2` bundle,
3. compute the static and `ω^2` grouped coefficients including the mixed-sector contribution,
4. and then evaluate the outgoing quadrupole normalization product on the isotropic branch.

That is the smallest next theorem target that directly attacks the remaining 2.5PN/4PN bottleneck rather than only reorganizing the conservative input data.

=== moving_throat_pde_stage005_grouped_p2_normalization_bridge.md ===

# Moving-Throat PDE — Stage 5: Grouped Real `P2` Bundle and the Normalized Outgoing-Quadrupole Bridge

## Purpose

Stage 4 gave the first honest one-lane bridge from

- the conservative wall/BdG/Maxwell/mixed block,
- to a passive/outgoing `l=2` fingerprint,
- with a positive static transfer factor.

That was enough to show **where** the 2.5PN quadrupole route can live, but not yet enough to speak the same language as the grouped real `P2` conservative packages or the 2.5PN normalization stack.

The next missing bridge is therefore not “solve the whole PDE.”
It is much sharper:

> lift the one-lane Stage-4 result to the grouped real `20/21/22` bundle, convert the Stage-3/4 conservative operator moments into the normalized grouped-response moments used by the 2.5PN package, and isolate the exact normalization product that still has to hit the universal quadrupole target.

That is the whole point of this stage.

---

## 1. Grouped real `P2` lane setup

Work lane by lane with

`A in {20,21,22}`.

After the conservative Stage-3 + Stage-4 reductions, each grouped lane is described by a conservative wall/worldtube operator

`D_A^(cons)(omega) = D_{A0} + D_{A2} omega^2 + D_{A4} omega^4 + O(omega^6)`.

This notation is intentionally broad.
All of the following are already absorbed into the coefficients `D_{A0}, D_{A2}, D_{A4}`:

- the geometry-only finite-throat wall backbone,
- the stable BdG support self-energy,
- the conservative localized-Maxwell / mixed-sector self-energy.

So Stage 5 does **not** reopen the earlier reductions. It takes them as the conservative front end.

The outgoing mixed-sector dressing from Stage 4 is encoded in a transfer factor

`N_A(omega) = N_{A0} + N_{A2} omega^2 + N_{A4} omega^4 + O(omega^6)`.

The key new point is that the grouped real `P2` bundle is now described by **two** low-frequency objects per lane:

1. the conservative operator `D_A^(cons)`,
2. the outgoing transfer factor `N_A`.

---

## 2. Exact bridge from conservative operator moments to grouped-response moments

Stage 3 normalized the conservative operator itself.
The 2.5PN grouped-quadrupole package, however, is phrased most naturally in terms of the normalized **response**

`Y_A(omega) = D_{A0} / D_A^(cons)(omega)`.

Expanding,

`Y_A(omega) = 1 + u_2^(A) omega^2 + u_4^(A) omega^4 + O(omega^6)`.

The exact conversion formulas are

`u_2^(A) = - D_{A2} / D_{A0}`,

`u_4^(A) = ( D_{A2}^2 - D_{A0} D_{A4} ) / D_{A0}^2`.

This is the first important bookkeeping bridge of the stage.
It tells us exactly how the conservative wall/BdG/Maxwell data from Stages 3–4 enter the grouped real `P2` response moments used by the 2.5PN summary.

So from this point onward:

- `D_{A0}, D_{A2}, D_{A4}` are the conservative operator moments,
- `u_2^(A), u_4^(A)` are the normalized grouped-response moments.

---

## 3. Exact grouped inverse map and isotropy defects

For any grouped triple `(x_20, x_21, x_22)`, define

`xbar = (x_20 + 2 x_21 + 2 x_22) / 5`,

`a_x = (2 x_20 - x_21 - x_22) / 10`,

`b_x = (x_21 - x_22) / 2`.

Then the exact inverse map is

`x_20 = xbar + 4 a_x`,

`x_21 = xbar - a_x + b_x`,

`x_22 = xbar - a_x - b_x`.

Applied to the grouped response moments, this gives

`ubar_2 = (u_2^(20) + 2 u_2^(21) + 2 u_2^(22)) / 5`,

`a_2 = (2 u_2^(20) - u_2^(21) - u_2^(22)) / 10`,

`b_2 = (u_2^(21) - u_2^(22)) / 2`,

and similarly for `(ubar_4, a_4, b_4)`.

The corresponding anisotropy norm is

`A_2^2 = 4 a_2^2 + (4/5) b_2^2`,

and similarly for `A_4^2`.

So the grouped real `P2` isotropy gate is once again completely sharp:

`a_2 = b_2 = a_4 = b_4 = 0`.

If the three grouped lanes share the same conservative moments, then the anisotropy defects vanish automatically.

---

## 4. The outgoing bridge is controlled by an exact internal prefactor

To first order in the outgoing branch, the response correction is controlled by the internal prefactor

`Pref_A(omega) = D_{A0} N_A(omega) / D_A^(cons)(omega)^2`.

Expand

`Pref_A(omega) = P_{A0} + P_{A2} omega^2 + P_{A4} omega^4 + O(omega^6)`.

Then the exact coefficients are

`P_{A0} = N_{A0} / D_{A0}`,

`P_{A2} = ( D_{A0} N_{A2} - 2 D_{A2} N_{A0} ) / D_{A0}^2`,

`P_{A4}`
`= [ D_{A0}^2 N_{A4} - 2 D_{A0}( D_{A2} N_{A2} + D_{A4} N_{A0} ) + 3 D_{A2}^2 N_{A0} ] / D_{A0}^3`.

This is the central algebraic result of the stage.

It cleanly separates the job of the conservative moving-throat block into two pieces:

- `D_A^(cons)` tells us the grouped response moments `u_2^(A), u_4^(A)`,
- `Pref_A` tells us how strongly the outgoing `l=2` branch is transferred into each grouped lane.

---

## 5. Multiplying by the compact outgoing `l=2` fingerprint

The compact outgoing `l=2` branch already has the normalized fingerprint

`Yhat_2^(out)(omega)`
`= 1 + A omega^2 + B omega^4 + i G5 omega^5 + O(omega^6)`

with

`A = a^2 / (9 c_s^2)`,

`B = 4 a^4 / (81 c_s^4)`,

`G5 = a^5 / (27 c_s^5)`.

So the outgoing contribution in grouped lane `A` is

`delta Y_A^(out)(omega) = Pref_A(omega) * Yhat_2^(out)(omega)`.

Write the branch coefficients as

`delta Y_A^(out)(omega)`
`= K_{A0} + K_{A2} omega^2 + K_{A4} omega^4 + i Gamma_5^(A) omega^5 + ...`.

Then

`K_{A0} = P_{A0}`,

`K_{A2} = P_{A2} + A P_{A0}`,

`K_{A4} = P_{A4} + A P_{A2} + B P_{A0}`,

`Gamma_5^(A) = G5 P_{A0}`.

This formula is extremely informative.
It shows that only the **static prefactor** `P_{A0}` enters the first odd `i omega^5` coefficient.
The higher prefactor moments `P_{A2}, P_{A4}` matter for the even branch bookkeeping, but not for the leading 2.5PN odd coefficient.

So the 2.5PN normalization problem is narrower than the full grouped conservative problem.
It only needs the correct isotropic value of `P_{A0}` on the true moving-throat branch.

---

## 6. Stage-4 one-lane prototype gives explicit `N0/N2/N4`

The Stage-4 one-lane Maxwell/mixed model has

`N(omega) = (P0_proto - g_W omega^2)^2 / (Delta0 - S2 omega^2 + omega^4)^2`,

where

`Delta0 = Omega_A^2 Omega_W^2 - R^2`,

`S2 = Omega_A^2 + Omega_W^2`,

`P0_proto = Omega_A^2 g_W + R g_A`.

Its exact low-frequency coefficients are

`N0 = P0_proto^2 / Delta0^2`,

`N2 = 2 P0_proto ( P0_proto S2 - Delta0 g_W ) / Delta0^3`,

`N4`
`= [ Delta0^2 g_W^2 - 2 Delta0 P0_proto^2 - 4 Delta0 P0_proto S2 g_W + 3 P0_proto^2 S2^2 ] / Delta0^4`.

So the one-lane Maxwell/mixed prototype already gives explicit microscopic data for the first three outgoing-transfer moments.

These are not yet the full coupled moving-throat bundle values, but they are the exact reduced prototype formulas that the true grouped-lane computation must generalize.

---

## 7. Isotropic grouped-lane collapse

If the full grouped real `P2` bundle is isotropic at the conservative level, then

`D_{20,n} = D_{21,n} = D_{22,n}`,

and if the outgoing transfer is isotropic as well, then

`N_{20,n} = N_{21,n} = N_{22,n}`

for each low-frequency coefficient we keep.

In that case,

`u_2^(20) = u_2^(21) = u_2^(22) = ubar_2`,

`u_4^(20) = u_4^(21) = u_4^(22) = ubar_4`,

`P_{20,n} = P_{21,n} = P_{22,n} = P_n`,

`K_{20,n} = K_{21,n} = K_{22,n} = K_n`,

and therefore all grouped anisotropy defects vanish:

`a_2 = b_2 = a_4 = b_4 = 0`.

So Stage 5 makes the 3PN and 2.5PN grouped stories meet cleanly:

- the 3PN grouped real `P2` isotropy gate is a statement about the conservative operator moments,
- the 2.5PN grouped quadrupole normalization is a statement about the isotropic static outgoing prefactor `P_0`.

---

## 8. The exact normalization product still to be hit

The universal GR quadrupole target is

`gamma_GR = 2 G / (5 c^5)`.

From the Stage-5 branch formula,

`Gamma_5 = P_0 * a^5 / (27 c_s^5)`.

Including the source-map factor `mhat_0`, the invariant normalization condition is therefore

`mhat_0^2 * P_0 * a^5 / (27 c_s^5) = 2 G / (5 c^5)`.

Equivalently,

`mhat_0^2 * P_0 = 54 G c_s^5 / (5 a^5 c^5)`.

This is the sharpest expression of the remaining theorem gap I can give at this stage.

Everything that came before was needed to make this quantity well defined.
Now the problem is no longer diffuse:

> compute the isotropic moving-throat value of the invariant product `mhat_0^2 P_0`.

On the natural source-map branch `mhat_0 = 1 + O(a^2/r^2)`, this reduces to the familiar 2.5PN target

`P_0 = 54 G c_s^5 / (5 a^5 c^5)`

at leading point-particle order.

If, in addition, the outgoing branch is carried by a **constant isotropic prefactor** (`P_2 = P_4 = 0`), then the even branch coefficients automatically become

`K_2 = P_0 a^2 / (9 c_s^2)`,

`K_4 = 4 P_0 a^4 / (81 c_s^4)`.

That is exactly the constant-prefactor branch isolated in the 2.5PN normalization package.

---

## 9. What this stage accomplishes physically

This stage does three things that were missing before.

### 9.1 It unifies the Stage-3/4 conservative data with the 2.5PN grouped language

Before this note, the conservative moving-throat derivation was naturally phrased in terms of operator moments.
The 2.5PN grouped summary, however, is phrased in terms of normalized grouped-response coefficients and branch-normalization data.

Stage 5 now gives the exact conversion formulas between those languages.

### 9.2 It identifies the true internal quantity that controls 2.5PN

The leading odd coefficient does **not** depend on every conservative detail equally.
At the universal point-particle 2.5PN level, the only internal quantity that matters is the isotropic static outgoing prefactor `P_0 = N_0 / D_0` (or `mhat_0^2 P_0` in invariant form).

That is a much sharper theorem target than “solve the whole quadrupole branch.”

### 9.3 It explains exactly what the full PDE still has to provide

The full moving-throat PDE still has to do two hard things:

1. produce the true grouped-lane conservative coefficients `D_{A0}, D_{A2}, D_{A4}` and outgoing transfer coefficients `N_{A0}, N_{A2}, N_{A4}` on the coupled isotropic branch,
2. and land the invariant product `mhat_0^2 P_0` on the target value above.

That is now a concrete job description for the missing PDE.

---

## 10. Immediate next step

The right next derivation is now narrower than before.

Do **not** jump yet to the full nonlinear theorem.
Instead:

1. keep the grouped-lane formulas from this stage,
2. compute the actual coupled grouped real `20/21/22` conservative coefficients for the wall/BdG/Maxwell/mixed system,
3. evaluate the corresponding `P_0, P_2, P_4` prefactor data on the isotropic branch,
4. and test whether the invariant normalization product `mhat_0^2 P_0` lands on the universal target.

That is the smallest next theorem gate that directly attacks the remaining 2.5PN/4PN bottleneck.

=== moving_throat_pde_stage006_full_grouped_bundle.md ===

# Moving-Throat PDE — Stage 6: Full Grouped Real `P2` Bundle, Exact Projectors, and the Isotropic Normalization Test

## Purpose

Stage 5 isolated the exact algebraic object that still has to hit the universal quadrupole target:

`mhat_0^2 P_0 = 54 G c_s^5 / (5 a^5 c^5)`.

But Stage 5 still spoke mostly in one-lane language. It showed how a single grouped lane talks to the outgoing `l=2` port, and it translated the result into the grouped-`P2` normalization stack. What it did **not** yet do was write the full coupled `20/21/22` bundle in one place and show exactly how the microscopic wall/BdG/Maxwell/mixed data feed

- the conservative grouped coefficients `D_{A0}, D_{A2}, D_{A4}`,
- the outgoing-transfer coefficients `N_{A0}, N_{A2}, N_{A4}`,
- the grouped isotropy defects,
- and the final isotropic normalization product.

That is the whole job of this stage.

The main new outputs are:

1. an exact weighted-projector calculus for the grouped real `P2` bundle,
2. explicit full-bundle formulas for the coupled low-frequency coefficients,
3. a clean separation of the microscopic sources of anisotropy,
4. the exact isotropic-branch normalization test in terms of the full coupled bundle,
5. and first-order anisotropy transport laws showing how conservative or transfer anisotropy leaks into the grouped response and outgoing prefactor.

So this stage is not yet the full nonlinear PDE theorem. It is the first exact bookkeeping layer in which the **entire** reduced `20/21/22` bundle is carried at once.

---

## 1. Exact weighted projector calculus for the grouped real `P2` bundle

The grouped real `P2` bundle is not an ordinary Euclidean triple. The natural grouped metric is

`Ggrp = diag(1,2,2)`.

This is the same weighting that already appeared in the grouped trace/anomaly formulas:

`xbar = (x_20 + 2 x_21 + 2 x_22)/5`,

`a_x = (2 x_20 - x_21 - x_22)/10`,

`b_x = (x_21 - x_22)/2`.

A particularly useful `Ggrp`-orthogonal basis is

`e_bar = (1,1,1)`,

`e_a   = (4,-1,-1)`,

`e_b   = (0,1,-1)`.

These vectors satisfy

`e_bar^T Ggrp e_a = 0`,

`e_bar^T Ggrp e_b = 0`,

`e_a^T   Ggrp e_b = 0`,

with norms

`e_bar^T Ggrp e_bar = 5`,

`e_a^T   Ggrp e_a   = 20`,

`e_b^T   Ggrp e_b   = 4`.

So for any grouped coefficient vector

`x = (x_20, x_21, x_22)^T`,

we have the exact decomposition

`x = xbar e_bar + a_x e_a + b_x e_b`.

Equivalently, the three exact projectors are

`P_bar = (1/5) [[1,2,2],[1,2,2],[1,2,2]]`,

`P_a   = (1/20) [[16,-8,-8],[-4,2,2],[-4,2,2]]`,

`P_b   = (1/4) [[0,0,0],[0,2,-2],[0,-2,2]]`,

so that

`P_bar + P_a + P_b = I3`,

`P_i P_j = delta_ij P_i`.

This is the cleanest exact way to organize grouped-lane isotropy.

- `P_bar` extracts the isotropic grouped trace,
- `P_a` and `P_b` extract the two independent anisotropy defects.

So at any stage of the coupled derivation, the grouped-lane isotropy condition is simply

`P_a x = P_b x = 0`.

---

## 2. Full coupled lane-by-lane reduced bundle

At linear order on an isotropic reference throat, the grouped real `20/21/22` channels do not mix directly with one another. Each lane is coupled internally, but the three lanes remain separate copies of the same reduced mechanism unless anisotropy is introduced.

So for each grouped lane

`A in {20,21,22}`

we keep:

- a wall/worldtube amplitude `q_A`,
- stable BdG support modes `X_(A,alpha)` with frequencies `varpi_(A,alpha)`,
- one or more localized brane-like gauge coordinates `U_(A,r)` with frequencies `Omega_(U,A,r)`,
- one or more mixed `A_w/F_{mu w}/J^w` coordinates `W_(A,r)` with frequencies `Omega_(W,A,r)`,
- internal gauge-sector mixing `R_(A,r)`.

The reduced quadratic Lagrangian in lane `A` is therefore

`L_A`
`= (1/2) M_A qdot_A^2 - (1/2) K_A q_A^2`
`  + sum_alpha [ (1/2) Xdot_(A,alpha)^2 - (1/2) varpi_(A,alpha)^2 X_(A,alpha)^2 + c_(A,alpha) q_A X_(A,alpha) ]`
`  + sum_r [ (1/2) Udot_(A,r)^2 - (1/2) Omega_(U,A,r)^2 U_(A,r)^2`
`          + (1/2) Wdot_(A,r)^2 - (1/2) Omega_(W,A,r)^2 W_(A,r)^2`
`          + R_(A,r) U_(A,r) W_(A,r)`
`          + g_(U,A,r) q_A U_(A,r) + g_(W,A,r) q_A W_(A,r) ]`.

This is the full reduced grouped-lane bundle we need for the next theorem gate.

It combines, in one place,

- the geometry wall backbone,
- the BdG support self-energy,
- the conservative localized-Maxwell/mixed self-energy,
- and the outgoing transfer channel.

---

## 3. Exact full-bundle low-frequency coefficients

### 3.1 BdG support contributions

For each lane `A`, define the exact BdG moments

`B_(A,0) = sum_alpha c_(A,alpha)^2 / varpi_(A,alpha)^2`,

`B_(A,2) = sum_alpha c_(A,alpha)^2 / varpi_(A,alpha)^4`,

`B_(A,4) = sum_alpha c_(A,alpha)^2 / varpi_(A,alpha)^6`.

These are just the Stage-3 conservative support moments collected lane by lane.

### 3.2 Conservative Maxwell/mixed contributions

For each port-active internal pair `(U_(A,r),W_(A,r))`, define

`Delta_(A,r) = Omega_(U,A,r)^2 Omega_(W,A,r)^2 - R_(A,r)^2`,

`S_(A,r) = Omega_(U,A,r)^2 + Omega_(W,A,r)^2`,

`Q_(A,r) = g_(U,A,r)^2 Omega_(W,A,r)^2 + 2 g_(U,A,r) g_(W,A,r) R_(A,r) + g_(W,A,r)^2 Omega_(U,A,r)^2`,

`G_(A,r) = g_(U,A,r)^2 + g_(W,A,r)^2`.

Then the conservative low-frequency coefficients contributed by that port are

`Z_(A,0)^(r) = Q_(A,r) / Delta_(A,r)`,

`Z_(A,2)^(r) = [ Q_(A,r) S_(A,r) - G_(A,r) Delta_(A,r) ] / Delta_(A,r)^2`,

`Z_(A,4)^(r) = [ Q_(A,r)(S_(A,r)^2 - Delta_(A,r)) - S_(A,r) G_(A,r) Delta_(A,r) ] / Delta_(A,r)^3`.

Summing over all internal Maxwell/mixed branches gives

`Z_(A,n) = sum_r Z_(A,n)^(r)`

for `n=0,2,4`.

### 3.3 Outgoing-transfer coefficients

For each port-active internal pair, the exact outgoing transfer factor begins with

`N_(A,0)^(r) = [ Omega_(U,A,r)^2 g_(W,A,r) + R_(A,r) g_(U,A,r) ]^2 / Delta_(A,r)^2`.

At the next two even orders,

`N_(A,2)^(r)`
`= 2 P_(A,r) [ P_(A,r) S_(A,r) - Delta_(A,r) g_(W,A,r) ] / Delta_(A,r)^3`,

`N_(A,4)^(r)`
`= [ Delta_(A,r)^2 g_(W,A,r)^2 - 2 Delta_(A,r) P_(A,r)^2 - 4 Delta_(A,r) P_(A,r) S_(A,r) g_(W,A,r) + 3 P_(A,r)^2 S_(A,r)^2 ] / Delta_(A,r)^4`,

where

`P_(A,r) = Omega_(U,A,r)^2 g_(W,A,r) + R_(A,r) g_(U,A,r)`.

Summing over outgoing ports gives

`N_(A,n) = sum_r N_(A,n)^(r)`

for `n=0,2,4`.

### 3.4 Total conservative wall operator coefficients

Putting everything together, the full coupled conservative grouped-lane operator is

`D_A^(cons)(omega) = D_(A,0) + D_(A,2) omega^2 + D_(A,4) omega^4 + O(omega^6)`

with

`D_(A,0) = K_A - B_(A,0) - Z_(A,0)`,

`D_(A,2) = -[ M_A + B_(A,2) + Z_(A,2) ]`,

`D_(A,4) = -[ B_(A,4) + Z_(A,4) ]`.

These are the exact full-bundle low-frequency coefficients the previous stage was asking for.

Nothing symbolic is left hidden here:

- `B_(A,n)` is the BdG support part,
- `Z_(A,n)` is the conservative localized-Maxwell/mixed part,
- `N_(A,n)` is the outgoing-transfer part.

So each grouped lane now has fully explicit microscopic reduced coefficients.

---

## 4. Exact grouped trace/anomaly bookkeeping for the microscopic coefficients

For any lane family `X_(A,n)` — for example `K_A`, `M_A`, `B_(A,0)`, `Z_(A,2)`, or `N_(A,0)` — define

`Xbar_n = (X_(20,n) + 2 X_(21,n) + 2 X_(22,n))/5`,

`a_(X,n) = (2 X_(20,n) - X_(21,n) - X_(22,n))/10`,

`b_(X,n) = (X_(21,n) - X_(22,n))/2`.

Then the full coupled coefficients decompose exactly as

`Dbar_0 = Kbar - Bbar_0 - Zbar_0`,

`a_(D,0) = a_K - a_(B,0) - a_(Z,0)`,

`b_(D,0) = b_K - b_(B,0) - b_(Z,0)`,

`Dbar_2 = -[ Mbar + Bbar_2 + Zbar_2 ]`,

`a_(D,2) = -[ a_M + a_(B,2) + a_(Z,2) ]`,

`b_(D,2) = -[ b_M + b_(B,2) + b_(Z,2) ]`,

`Dbar_4 = -[ Bbar_4 + Zbar_4 ]`,

`a_(D,4) = -[ a_(B,4) + a_(Z,4) ]`,

`b_(D,4) = -[ b_(B,4) + b_(Z,4) ]`.

Similarly, for the outgoing-transfer bundle,

`Nbar_n = (N_(20,n) + 2 N_(21,n) + 2 N_(22,n))/5`,

`a_(N,n) = (2 N_(20,n) - N_(21,n) - N_(22,n))/10`,

`b_(N,n) = (N_(21,n) - N_(22,n))/2`.

This is the first exact microscopic answer to the question:

> where can grouped-lane anisotropy enter?

Answer: only through anisotropy in the wall baseline (`K_A,M_A`), the BdG support data (`B_(A,n)`), or the Maxwell/mixed conservative and outgoing-transfer bundles (`Z_(A,n), N_(A,n)`).

That is the complete reduced anisotropy ledger.

---

## 5. Isotropic branch and the exact full-bundle normalization test

The natural isotropic grouped-lane branch is the one on which all grouped anisotropy defects vanish:

`a_(D,0)=b_(D,0)=a_(D,2)=b_(D,2)=a_(D,4)=b_(D,4)=0`,

`a_(N,0)=b_(N,0)=a_(N,2)=b_(N,2)=a_(N,4)=b_(N,4)=0`.

Then the three grouped lanes collapse to common coefficients

`D_(20,n)=D_(21,n)=D_(22,n)=D_n`,

`N_(20,n)=N_(21,n)=N_(22,n)=N_n`.

On that branch, the exact Stage-5 formulas become

`u_2 = - D_2 / D_0`,

`u_4 = (D_2^2 - D_0 D_4) / D_0^2`,

`P_0 = N_0 / D_0`,

`P_2 = (D_0 N_2 - 2 D_2 N_0) / D_0^2`,

`P_4 = [ D_0^2 N_4 - 2 D_0 ( D_2 N_2 + D_4 N_0 ) + 3 D_2^2 N_0 ] / D_0^3`.

The universal normalization condition is therefore

`mhat_0^2 P_0 = 54 G c_s^5 / (5 a^5 c^5)`.

Using the explicit full-bundle coefficient definitions above, this becomes

`mhat_0^2 * N_0 / ( K - B_0 - Z_0 ) = 54 G c_s^5 / (5 a^5 c^5)`.

This is the exact isotropic normalization test for the full coupled reduced bundle.

So the remaining theorem gap is now completely sharp in reduced microscopic language:

> does the true moving-throat branch produce coupled bundle data `(K, M, B_n, Z_n, N_n)` whose isotropic static ratio `N_0/(K-B_0-Z_0)` lands on the universal target?

That is the whole question.

---

## 6. Constant-prefactor branch conditions

The 2.5PN normalization package repeatedly singles out the especially simple branch on which the outgoing prefactor is constant at the orders we retain.

In the present full-bundle notation, that means

`P_2 = 0`,

`P_4 = 0`.

The exact algebraic conditions are therefore

`N_2 = 2 D_2 N_0 / D_0`,

`N_4 = [ 2 D_0 ( D_2 N_2 + D_4 N_0 ) - 3 D_2^2 N_0 ] / D_0^2`.

So the constant-prefactor outgoing branch does **not** require the higher coefficients to vanish separately. It requires the outgoing-transfer moments `N_2,N_4` to be correlated with the conservative bundle moments `D_0,D_2,D_4` in exactly the way above.

That is another exact reduced theorem target for the completed PDE.

---

## 7. First-order anisotropy transport laws

One of the most useful things we can do at this stage is linearize around the isotropic branch.

Write

`D_(A,n) = D_n + delta D_(A,n)`,

`N_(A,n) = N_n + delta N_(A,n)`

with `delta D` and `delta N` small grouped-lane anisotropy corrections.

Then the first two objects that matter most are:

### 7.1 Response-moment anisotropy

Since

`u_2^(A) = - D_(A,2) / D_(A,0)`,

the first-order variation is

`delta u_2^(A) = -[ delta D_(A,2) + u_2 delta D_(A,0) ] / D_0`.

So the grouped anisotropy defects of the conservative response obey

`a_(u,2) = -[ a_(D,2) + u_2 a_(D,0) ] / D_0`,

`b_(u,2) = -[ b_(D,2) + u_2 b_(D,0) ] / D_0`.

This shows that even if the `omega^2` coefficient itself were isotropic, anisotropy in the static coefficient `D_0` would still leak into the normalized grouped response.

### 7.2 Outgoing-prefactor anisotropy

Since

`P_0^(A) = N_(A,0) / D_(A,0)`,

the first-order variation is

`delta P_0^(A) = [ delta N_(A,0) - P_0 delta D_(A,0) ] / D_0`.

Therefore

`a_(P,0) = [ a_(N,0) - P_0 a_(D,0) ] / D_0`,

`b_(P,0) = [ b_(N,0) - P_0 b_(D,0) ] / D_0`.

This is the most useful anisotropy diagnostic of the stage.

It says the outgoing-prefactor isotropy can fail in only two ways:

1. direct anisotropy in the outgoing transfer bundle `N_(A,0)`,
2. conservative anisotropy in `D_(A,0)` that is reweighted by the isotropic prefactor `P_0`.

So when the next PDE step is computed, we will know immediately whether any anisotropy is entering through the support side, the Maxwell/mixed transfer side, or both.

---

## 8. Expert reading of the normalization formula

The exact isotropic normalization product

`P_0 = N_0 / D_0`

has a simple physical reading once the branch is stable (`D_0 > 0`).

- Increasing the static outgoing-transfer weight `N_0` raises `P_0` directly.
- Increasing the conservative support renormalizations `B_0` or `Z_0` lowers the denominator `D_0 = K-B_0-Z_0` and therefore also raises `P_0`.

So if the natural branch undershoots the universal target, there are only two broad ways to move toward it:

1. increase the outgoing transfer strength,
2. or soften the static wall operator through conservative support dressing while keeping the branch stable.

This does **not** solve the theorem, but it tells us how the microscopic bundle parameters have to move if the completed PDE is to land on the right normalization.

---

## 9. What this stage achieves physically

This stage closes the first full-bundle bookkeeping gap.

### 9.1 It gives the exact coupled grouped-lane coefficients

The bundle coefficients `D_(A,0), D_(A,2), D_(A,4)` and `N_(A,0), N_(A,2), N_(A,4)` are now written explicitly in terms of

- wall data,
- BdG support data,
- conservative localized-Maxwell/mixed data,
- and outgoing-transfer data.

So the next theorem gate is not missing algebra anymore.

### 9.2 It localizes the sources of anisotropy exactly

The projector calculus makes it impossible to hide where grouped-lane anisotropy comes from.
Each microscopic sector carries its own exact `(bar,a,b)` defects, and the full bundle inherits them linearly.

### 9.3 It turns the final normalization problem into one ratio

On the natural isotropic branch, everything collapses to the exact ratio

`P_0 = N_0 / D_0`.

That is the reduced microscopic version of the remaining 2.5PN/4PN normalization bottleneck.

---

## 10. Immediate next step

The right next derivation is now even narrower than before.

Do **not** reopen the whole theory at once.
Instead:

1. compute the actual overlap integrals that determine the full-bundle data `B_(A,n), Z_(A,n), N_(A,n)` on the true moving-throat branch,
2. project them with `(P_bar, P_a, P_b)`,
3. check whether the grouped anisotropy defects vanish on the natural branch,
4. and then test the single ratio

`mhat_0^2 N_0 / (K-B_0-Z_0)`

against the universal target.

That is now the smallest honest next theorem gate.

=== moving_throat_pde_stage007_overlap_isotropy.md ===

# Moving-Throat PDE — Stage 7: Explicit Overlap Integrals, the `O(3)` Isotropy Theorem, and the First Axisymmetric Splitting Law

## Purpose

Stage 6 reduced the entire grouped real `20/21/22` bundle to one exact question:

`mhat_0^2 N_0 / (K - B_0 - Z_0) = 54 G c_s^5 / (5 a^5 c^5)`.

But Stage 6 still treated the microscopic inputs

- `B_(A,n)`,
- `Z_(A,n)`,
- `N_(A,n)`

as formal coefficients.

The next honest step is therefore **not** another abstract ratio manipulation. It is to write the first actual overlap-integral model for the grouped real `P2` bundle and ask what the angular geometry already forces before the unsolved radial/axial PDE is touched.

This stage does exactly that.

The main outputs are:

1. an exact normalized real-STF harmonic basis on the throat sphere,
2. an exact angular source-map identity for the grouped `20/21/22` quadrupole ports,
3. an `O(3)` isotropy theorem showing that every isotropic reduced kernel collapses the three grouped lanes to one common value,
4. explicit common radial/axial overlap formulas for `B_n`, `Z_n`, and `N_n`,
5. the first exact symmetry-breaking law for a weak **axisymmetric quadrupole** perturbation,
6. and the resulting first-order transport law for the normalization ratio `P_A = N_A / D_A`.

So Stage 7 is the first point where the words “actual overlap integrals on the moving-throat branch” become mathematically concrete.

---

## 1. Normalized real STF harmonics on the throat sphere

Let `n^i` be the unit direction on the throat sphere and let

`E_(20), E_(21c), E_(21s), E_(22c), E_(22s)`

be the canonical real STF basis already used in the 2.5PN audit. Define the real angular functions

`Y_A(n) = sqrt(15/(8 pi)) E_A^{ij} n_i n_j`,

with `A in {20,21c,21s,22c,22s}`.

The exact fourth-moment identity on the unit sphere gives

`int dOmega n_i n_j n_k n_l = (4 pi / 15) (delta_ij delta_kl + delta_ik delta_jl + delta_il delta_jk)`.

Since the STF basis satisfies `Tr(E_A E_B) = delta_AB`, the normalized real harmonics obey

`int dOmega Y_A Y_B = delta_AB`.

So this basis is orthonormal without any further angular normalization choices.

That one fact already matters for the theorem gate, because it means the grouped real quadrupole bundle has a **canonical angular basis** on the isotropic throat.

---

## 2. Exact angular source-map identity

Write the orbital/worldtube STF quadrupole source on the throat sphere as

`S(n) = sum_A S_A Y_A(n)`.

Projecting onto the same port basis gives

`S_A^(port) = int dOmega Y_A(n) S(n) = sum_B S_B int dOmega Y_A Y_B = S_A`.

So on the natural isotropic source map the **angular** matching matrix is exactly the identity.

This means that the still-open source normalization `mhat_0` factors as

`mhat_0 = mhat_rad * mhat_ang`,

with

`mhat_ang = 1`

exactly on the canonical isotropic branch.

So the remaining normalization gap is not an angular mismatch between the orbital STF source and the grouped real throat ports. It is a **radial/axial and port-amplitude** issue.

That is a real narrowing of the theorem problem.

---

## 3. Explicit overlap-integral factorization on the isotropic branch

Take the grouped real wall deformation to be expanded as

`eta_A(s,Omega,t) = q_A(t) chi_eta(s) Y_A(Omega)`.

Take the stable BdG support sector in the same `l=2` bundle as

`X_(alpha,A)(s,Omega,t) = X_(alpha,A)(t) phi_alpha(s) Y_A(Omega)`.

Take the conservative brane-like gauge and mixed `A_w/F_(mu w)/J^w` sectors as

`U_(r,A)(s,Omega,t) = U_(r,A)(t) u_r(s) Y_A(Omega)`,

`W_(r,A)(s,Omega,t) = W_(r,A)(t) w_r(s) Y_A(Omega)`.

Assume the reference throat and the reduced kernels are `O(3)` invariant, so that the angular dependence enters only through scalar contractions.

Then all angular overlaps collapse by orthonormality, and every microscopic coupling becomes lane-independent:

`c_(A,alpha) = C_alpha = lambda_(B,alpha) I_(eta,alpha)`,

`g_(U,A,r) = G_(U,r) = lambda_(U,r) I_(eta,u,r)`,

`g_(W,A,r) = G_(W,r) = lambda_(W,r) I_(eta,w,r)`,

`R_(A,r) = R_r = lambda_(R,r) I_(u,w,r)`.

The radial/axial overlap integrals are

`I_(eta,alpha) = int ds mu_s(s) chi_eta(s) phi_alpha(s)`,

`I_(eta,u,r)  = int ds mu_s(s) chi_eta(s) u_r(s)`,

`I_(eta,w,r)  = int ds mu_s(s) chi_eta(s) w_r(s)`,

`I_(u,w,r)    = int ds mu_s(s) u_r(s) w_r(s)`.

So on the isotropic branch the full Stage-6 coefficients become true scalar lane-independent objects:

`B_(A,0) = B_0 = sum_alpha C_alpha^2 / varpi_alpha^2`,

`B_(A,2) = B_2 = sum_alpha C_alpha^2 / varpi_alpha^4`,

`B_(A,4) = B_4 = sum_alpha C_alpha^2 / varpi_alpha^6`.

For each conservative Maxwell/mixed pair `r`, define

`Delta_r = Omega_(U,r)^2 Omega_(W,r)^2 - R_r^2`,

`S_r = Omega_(U,r)^2 + Omega_(W,r)^2`,

`Q_r = G_(U,r)^2 Omega_(W,r)^2 + 2 G_(U,r) G_(W,r) R_r + G_(W,r)^2 Omega_(U,r)^2`,

`H_r = G_(U,r)^2 + G_(W,r)^2`,

`P_r = Omega_(U,r)^2 G_(W,r) + R_r G_(U,r)`.

Then

`Z_0 = sum_r Q_r / Delta_r`,

`Z_2 = sum_r [ Q_r S_r - H_r Delta_r ] / Delta_r^2`,

`Z_4 = sum_r [ Q_r (S_r^2 - Delta_r) - S_r H_r Delta_r ] / Delta_r^3`,

and

`N_0 = sum_r P_r^2 / Delta_r^2`,

`N_2 = sum_r 2 P_r (P_r S_r - Delta_r G_(W,r)) / Delta_r^3`,

`N_4 = sum_r [ Delta_r^2 G_(W,r)^2 - 2 Delta_r P_r^2 - 4 Delta_r P_r S_r G_(W,r) + 3 P_r^2 S_r^2 ] / Delta_r^4`.

So the conservative wall operator is exactly

`D_A(omega) = D_0 + D_2 omega^2 + D_4 omega^4 + O(omega^6)`

with

`D_0 = K - B_0 - Z_0`,

`D_2 = -(M + B_2 + Z_2)`,

`D_4 = -(B_4 + Z_4)`.

The grouped-lane consequence is immediate:

`D_(20,n) = D_(21,n) = D_(22,n) = D_n`,

`N_(20,n) = N_(21,n) = N_(22,n) = N_n`,

and therefore

`a_(D,n) = b_(D,n) = 0`,

`a_(N,n) = b_(N,n) = 0`.

So inside any truly `O(3)`-invariant reduced kernel, the grouped real `20/21/22` bundle is forced to be isotropic.

This is the first honest PDE-side isotropy theorem reached in the program, even though it is still a reduced-sector theorem rather than a full nonlinear moving-throat theorem.

---

## 4. The first allowed symmetry-breaking channel: a weak axisymmetric quadrupole background

Once isotropy is understood, the next question is the first way it can fail.

The leading symmetry-breaking correction that talks directly to the grouped real `P2` bundle is a weak quadrupolar anisotropy. In the axisymmetric frame, write the perturbation as

`delta K = eps kappa(s) Y_20(Omega)`.

The angular matrix that perturbs any `l=2` overlap is then

`M_AB^(20) = int dOmega Y_A Y_20 Y_B`.

Using the exact sixth moment of the unit sphere,

`int dOmega n_i n_j n_k n_l n_m n_n = (4 pi / 105) sum_pairings delta delta delta`,

one finds the exact five-mode result

`M^(20) = kappa_* diag(1, 1/2, 1/2, -1, -1)`,

with

`kappa_* = sqrt(5) / (7 sqrt(pi))`.

So the weak axisymmetric quadrupole background does **not** produce an arbitrary lane splitting. It produces one universal angular fingerprint:

- `20` lane shift proportional to `+1`,
- `21c,21s` shifts proportional to `+1/2`,
- `22c,22s` shifts proportional to `-1`.

After regrouping the `c/s` pairs into the three grouped lanes, the pattern is

`lambda_(20) = 1`,

`lambda_(21) = 1/2`,

`lambda_(22) = -1`.

Any first-order microscopic coefficient on that branch therefore has the form

`x_(20) = x^(0) + eps x^(1)`,

`x_(21) = x^(0) + (eps/2) x^(1)`,

`x_(22) = x^(0) - eps x^(1)`.

The weighted grouped trace/anomaly variables are then

`xbar = x^(0)`,

`a_x = (eps/4) x^(1)`,

`b_x = (3 eps / 4) x^(1)`.

So the first axisymmetric symmetry-breaking law is

`b_x = 3 a_x`.

This is a strong and very usable result. It means that if a future PDE computation shows a weak grouped-lane anisotropy produced by an axisymmetric `l=2` distortion of the isotropic branch, the defects are **not** free. They must sit on the one-dimensional line

`b = 3 a`.

If that relation fails, the symmetry breaking is not a pure axisymmetric quadrupole perturbation. It must involve either

- non-axisymmetric `m != 0` structure,
- higher-rank angular content,
- or a more complicated non-separable reduced kernel.

---

## 5. First-order normalization transport on the weak axisymmetric branch

Now apply the same axisymmetric law to the Stage-6 normalization ratio.

Suppose

`D_A = D_0 + eps lambda_A D_1 + O(eps^2)`,

`N_A = N_0 + eps lambda_A N_1 + O(eps^2)`,

with

`lambda_(20)=1`,

`lambda_(21)=1/2`,

`lambda_(22)=-1`.

Then the grouped-lane prefactor becomes

`P_A = N_A / D_A = P_0 + eps lambda_A P_1 + O(eps^2)`

with

`P_0 = N_0 / D_0`,

`P_1 = (N_1 D_0 - N_0 D_1) / D_0^2`.

So the first normalization anisotropies obey the same exact defect law:

`abar_P = 0`,

`a_P = (eps/4) P_1`,

`b_P = (3 eps / 4) P_1`,

`b_P = 3 a_P`.

This is the first actual transport law for the grouped normalization test once symmetry is weakly broken.

It says that the isotropic normalization target from Stage 6 is stable in the expected way: weak axisymmetric anisotropy does not create an arbitrary three-parameter deformation of the normalization stack. It creates one universal first-order lane pattern.

---

## 6. What Stage 7 changes in the theorem problem

Stage 7 narrows the remaining gap in three important ways.

### 6.1 The angular part of the source map is no longer open

On the natural isotropic grouped real basis,

`mhat_ang = 1`

exactly.

So the remaining source normalization issue is radial/axial and dynamical, not angular.

### 6.2 The isotropy theorem is now explicit on the reduced PDE side

If the reference throat and reduced kernels are truly `O(3)` invariant, then the grouped real `20/21/22` bundle is forced to collapse exactly:

`u_2^(20) = u_2^(21) = u_2^(22)`,

`u_4^(20) = u_4^(21) = u_4^(22)`,

`P_0^(20) = P_0^(21) = P_0^(22)`.

So on that branch the remaining theorem gap is not isotropy itself. It is the actual scalar value of the common radial/axial overlap data.

### 6.3 The first symmetry-breaking pattern is now diagnostic

A weak axisymmetric `l=2` deformation predicts one exact grouped signature:

`(20,21,22) ~ (1, 1/2, -1)`,

or equivalently

`b = 3 a`.

So future PDE data can be classified immediately:

- if weak anisotropy obeys `b = 3 a`, it is consistent with a pure axisymmetric quadrupole perturbation of the isotropic branch;
- if not, the symmetry breaking must be more complicated.

---

## 7. Best current summary after Stage 7

The moving-throat PDE problem is now narrower than it was at the end of Stage 6.

The remaining higher-order bridge is no longer:

> somehow compute all grouped-lane coefficients.

It is now:

1. compute the **common radial/axial isotropic overlap integrals** on the natural branch,
2. insert them into
   `B_n`, `Z_n`, `N_n`,
3. test the single isotropic ratio
   `mhat_rad^2 N_0 / (K - B_0 - Z_0)`,
4. and then, only after that, study symmetry breaking around that branch using the exact Stage-7 angular fingerprints.

That is a real tightening of the theorem target.

=== moving_throat_pde_stage008_minimal_isotropic_normalization.md ===

# Moving-Throat PDE — Stage 8: Minimal Isotropic Single-Mode Closure and the Explicit Normalization Formula

## Purpose

Stage 7 removed the angular ambiguity.
On the natural isotropic grouped basis,

`mhat_ang = 1`

exactly, and the grouped `20/21/22` bundle collapses to a single common lane.

That means the next honest simplification is not angular. It is **radial/axial**.

This stage therefore freezes the minimal isotropic closure with

- one wall/worldtube quadrupole mode,
- one stable BdG support mode,
- one conservative localized-Maxwell/mixed internal pair,
- and one passive outgoing port.

The point is not to claim that the true moving-throat branch is literally one mode deep.
The point is to collapse the remaining theorem gap to the smallest exact formula that still carries all the physics.

The main output is an explicit isotropic normalization law:

`mhat_rad^2 P^2 / [ Delta (K Delta - Delta C^2 / varpi^2 - Q) ] = 54 G c_s^5 / (5 a^5 c^5)`.

So after Stage 8, the minimal isotropic branch is no longer “some unknown PDE normalization.”
It is one exact scalar equation in the radial/axial overlap amplitudes.

---

## 1. Minimal isotropic radial/axial data

Keep the Stage-7 isotropic grouped collapse and truncate the common lane to

- one BdG support mode with frequency `varpi`,
- one brane-like internal gauge coordinate `U`,
- one mixed `A_w/F_(mu w)/J^w` coordinate `W`.

Define the common radial/axial overlap amplitudes

`C   = lambda_B I_(eta,phi)`,

`G_U = lambda_U I_(eta,u)`,

`G_W = lambda_W I_(eta,w)`,

`R   = lambda_R I_(u,w)`.

The mode frequencies are

`varpi`,

`Omega_U`,

`Omega_W`.

As before, define

`Delta = Omega_U^2 Omega_W^2 - R^2`,

`Q = G_U^2 Omega_W^2 + 2 G_U G_W R + G_W^2 Omega_U^2`,

`P = Omega_U^2 G_W + R G_U`.

Then the common zero-frequency isotropic coefficients are

`B_0 = C^2 / varpi^2`,

`Z_0 = Q / Delta`,

`N_0 = P^2 / Delta^2`.

So the common conservative wall operator is

`D_0 = K - C^2 / varpi^2 - Q / Delta`.

---

## 2. Exact minimal isotropic normalization ratio

The grouped-lane normalization prefactor is still

`P_0 = N_0 / D_0`.

Substituting the minimal isotropic coefficients gives the exact closed expression

`P_0 = P^2 / [ Delta^2 ( K - C^2 / varpi^2 - Q / Delta ) ]`.

Equivalently,

`P_0 = P^2 / [ Delta ( K Delta - Delta C^2 / varpi^2 - Q ) ]`.

Because Stage 7 already proved `mhat_ang = 1`, the remaining source normalization is purely radial/axial:

`mhat_0 = mhat_rad`.

So the full isotropic target becomes

`mhat_rad^2 P^2 / [ Delta ( K Delta - Delta C^2 / varpi^2 - Q ) ] = 54 G c_s^5 / (5 a^5 c^5)`.

This is the sharpest explicit reduced normalization formula reached so far.

---

## 3. Stability and positivity conditions on the minimal branch

The minimal isotropic closure is physically admissible only if

`Delta > 0`

and

`D_0 > 0`.

These are just the statements that

- the conservative internal `(U,W)` block has not crossed its static instability, and
- the wall/worldtube quadrupole mode has not been driven through zero stiffness by the support and internal electromagnetic self-energies.

In compact form the stability condition is

`K Delta - Delta C^2 / varpi^2 - Q > 0`.

On that stable branch,

`N_0 = P^2 / Delta^2 >= 0`.

So the prefactor sign is controlled entirely by the sign of the stable denominator. On the admissible branch,

`P_0 > 0`

whenever

`P != 0`.

So the minimal isotropic closure gives a very clear physical criterion for the existence of a nontrivial quadrupole bridge:

- the port-transfer combination `P` must not vanish,
- and the stable denominator must remain positive.

---

## 4. Exact monotonicity on the minimal isotropic branch

The minimal formula already exposes how the remaining normalization can move.

Define the support-softening variable

`X = C^2 / varpi^2`.

Then

`P_0 = N_0 / (K - X - Q / Delta)`.

So the exact derivatives are

`partial_K P_0 = - N_0 / D_0^2`,

`partial_X P_0 = + N_0 / D_0^2`.

Thus, on the stable branch,

- increasing the bare wall stiffness `K` decreases the normalization prefactor,
- increasing the BdG support softening `X` increases it.

The second statement is useful because it makes the normalization target operational: stronger support dressing pushes the bridge **up**, while a stiffer bare wall pushes it **down**.

The internal Maxwell/mixed sector affects both numerator and denominator. Writing

`P_0 = P^2 / [ Delta (K Delta - Delta X - Q) ]`,

the numerator contribution is controlled by the transfer combination

`P = Omega_U^2 G_W + R G_U`,

while the denominator load is controlled by

`Q = G_U^2 Omega_W^2 + 2 G_U G_W R + G_W^2 Omega_U^2`.

So the normalization target is not asking for generic “more coupling.”
It is asking for the right balance between

- port transfer into the outgoing quadrupole branch,
- and conservative self-loading of the wall.

---

## 5. What this means for the actual PDE task

After Stage 8, the minimal isotropic theorem gate is no longer vague.
The completed moving-throat PDE does not need to output an unlimited family of unknown coefficients before the normalization question can even be asked.

On the minimal isotropic branch it needs only to determine the radial/axial amplitudes entering

`X = C^2 / varpi^2`,

`Delta`,

`Q`,

`P`,

and `mhat_rad`.

Then the normalization question is exactly

`mhat_rad^2 P^2 / [ Delta (K Delta - Delta X - Q) ] ?= 54 G c_s^5 / (5 a^5 c^5)`.

That is a real reduction of the theorem gap.

It means the next PDE computation can be judged immediately:

- if it lands on the target, the passive/outgoing grouped quadrupole bridge is closed on the minimal isotropic branch;
- if not, the failure is not mysterious — it is in one of the radial/axial quantities `X`, `Delta`, `Q`, `P`, or `mhat_rad`.

---

## 6. Best current summary after Stage 8

The road to the moving-throat PDE is now split cleanly into two layers.

### Angular layer
Already closed on the natural isotropic branch:

- canonical real STF harmonics,
- exact source-map identity,
- exact grouped isotropy theorem,
- exact weak axisymmetric splitting law.

### Radial/axial layer
Still open, but now reduced to an explicit formula.
On the minimal isotropic branch, the remaining higher-order bridge is the one scalar equation

`mhat_rad^2 P^2 / [ Delta (K Delta - Delta C^2 / varpi^2 - Q) ] = 54 G c_s^5 / (5 a^5 c^5)`.

That is the sharpest reduced theorem target reached so far.

=== moving_throat_pde_stage009_concrete_axial_overlaps.md ===

# Moving-Throat PDE — Stage 9: Concrete Finite-Throat Axial Modes, Exact Overlap Constants, and the First Branch-Level Normalization Test

## Purpose

Stage 8 reduced the remaining theorem gap to one scalar expression in the radial/axial overlap data,

`mhat_rad^2 P^2 / [ Delta ( K Delta - Delta C^2 / varpi^2 - Q ) ] = 54 G c_s^5 / (5 a^5 c^5)`.

But the quantities

- `I_(eta,phi)`,
- `I_(eta,u)`,
- `I_(eta,w)`,
- `I_(u,w)`

were still formal integrals.

The next honest step is therefore to pick a **concrete finite-throat axial mode problem**, compute those overlaps exactly, and substitute them into the Stage-8 formula.

The main result of this stage is that on the first natural mixed axial branch the whole overlap problem collapses to a single geometric constant

`kappa = 2 sqrt(2) / pi`.

That turns the Stage-8 normalization target into one explicit algebraic relation among

- the constant wall quadrupole stiffness,
- the BdG support coupling,
- the brane-like internal gauge frequency,
- the mixed-channel frequency,
- and the three reduced couplings `(lambda_B, lambda_U, lambda_W, lambda_R)`.

So after this stage the remaining gap is not “unknown overlaps in general.”
It is the value of a small number of branch parameters on one concrete finite-throat mode family.

---

## 1. Concrete finite-throat axial mode problem

Take the finite throat interval

`s in [0,L]`.

Use the flat axial measure

`mu_s(s) = 1`.

### 1.1 Constant N/N zero mode

For collective wall shape and the brane-like internal gauge coordinate, take the lowest Neumann/Neumann mode

`u_0(s) = 1 / sqrt(L)`.

It obeys

`-u_0'' = 0`,

`u_0'(0) = u_0'(L) = 0`,

and is normalized by

`int_0^L ds u_0^2 = 1`.

This is the simplest axial profile for a throat deformation or brane-like zero mode that remains nonzero at the mouth.

### 1.2 D/N half-wave ladder

For the trapped support channel and the mixed `A_w / F_(mu w) / J^w` channel, use the exact finite-throat D/N ladder already isolated in the frozen-wall benchmark:

`f_n'' + k_n^2 f_n = 0`,

`f_n(0) = 0`,

`f_n'(L) = 0`,

with normalized solutions

`f_n(s) = sqrt(2/L) sin( (n + 1/2) pi s / L )`,

`k_n = (n + 1/2) pi / L`.

So the lowest trapped axial support/mixed profile is

`f_0(s) = sqrt(2/L) sin( pi s / (2L) )`.

The associated support and mixed frequencies are taken to be

`varpi_n^2 = M_phi^2 + c_phi^2 k_n^2`,

`Omega_(W,n)^2 = M_W^2 + c_W^2 k_n^2`.

On the minimal branch we keep only `n=0`.

### 1.3 Wall and brane-like profile choice on the first concrete branch

The first concrete isotropic branch is therefore

- wall quadrupole axial profile: `chi_eta = u_0`,
- BdG support axial profile: `phi = f_0`,
- brane-like internal gauge profile: `u = u_0`,
- mixed-channel axial profile: `w = f_0`.

On that branch we write

`varpi = varpi_0`,

`Omega_W = Omega_(W,0)`,

and keep `Omega_U` as the reduced internal restoring frequency of the brane-like zero mode `u_0` (so its frequency need not come from an axial gradient term).

This is the simplest branch that

- keeps the wall deformation nonzero at the mouth,
- reuses the exact D/N finite-throat support ladder,
- and keeps the brane-like gauge coordinate on the natural zero mode.

---

## 2. Exact axial overlap constants

The exact overlap of the constant zero mode with the D/N ladder is

`kappa_n = int_0^L ds u_0(s) f_n(s) = sqrt(2) / ( (n + 1/2) pi )`.

So the lowest branch value is

`kappa = kappa_0 = 2 sqrt(2) / pi`.

Numerically,

`kappa ~= 0.900316316`.

The needed minimal-branch overlap integrals are therefore

`I_(eta,phi) = int_0^L ds u_0 f_0 = kappa`,

`I_(eta,u)   = int_0^L ds u_0 u_0 = 1`,

`I_(eta,w)   = int_0^L ds u_0 f_0 = kappa`,

`I_(u,w)     = int_0^L ds u_0 f_0 = kappa`.

So on this concrete branch all radial/axial overlap information collapses to the single geometric constant `kappa`.

There is also an exact axial selection rule lurking here:

- the constant wall/brane mode couples to the `n`th D/N throat wave with strength `kappa_n ~ 1/(n+1/2)`,
- so the lowest trapped half-wave `n=0` is automatically the dominant branch,
- and higher support/mixed D/N levels are axially suppressed.

That is the first real finite-throat hierarchy statement for the overlap problem.

---

## 3. Explicit minimal-branch coefficients

Define the reduced couplings exactly as in Stage 8,

`C   = lambda_B I_(eta,phi)`,

`G_U = lambda_U I_(eta,u)`,

`G_W = lambda_W I_(eta,w)`,

`R   = lambda_R I_(u,w)`.

On the concrete branch this becomes

`C   = kappa lambda_B`,

`G_U = lambda_U`,

`G_W = kappa lambda_W`,

`R   = kappa lambda_R`.

So the minimal-branch Stage-8 quantities are

`Delta = Omega_U^2 Omega_W^2 - kappa^2 lambda_R^2`,

`Q = lambda_U^2 Omega_W^2 + 2 kappa^2 lambda_U lambda_W lambda_R + kappa^2 lambda_W^2 Omega_U^2`,

`P = kappa ( Omega_U^2 lambda_W + lambda_R lambda_U )`.

The BdG softening becomes

`X = C^2 / varpi^2 = kappa^2 lambda_B^2 / varpi^2`.

Therefore

`B_0 = kappa^2 lambda_B^2 / varpi^2`,

`Z_0 = Q / Delta`,

`N_0 = kappa^2 ( Omega_U^2 lambda_W + lambda_R lambda_U )^2 / Delta^2`.

The conservative wall coefficient is

`D_0 = K - kappa^2 lambda_B^2 / varpi^2 - Q / Delta`.

---

## 4. First branch-level normalization test

Substituting the concrete overlaps into the Stage-8 target gives

`mhat_rad^2 kappa^2 ( Omega_U^2 lambda_W + lambda_R lambda_U )^2`
`/ [ Delta ( K Delta - Delta kappa^2 lambda_B^2 / varpi^2 - Q ) ]`
`= 54 G c_s^5 / (5 a^5 c^5)`.

This is the first fully explicit branch-level normalization equation derived from an actual finite-throat mode problem.

It can be solved exactly for the required wall stiffness:

`K_req = kappa^2 lambda_B^2 / varpi^2 + Q / Delta`
`      + mhat_rad^2 kappa^2 ( Omega_U^2 lambda_W + lambda_R lambda_U )^2`
`        / [ (54 G c_s^5 / (5 a^5 c^5)) Delta^2 ]`.

So on this branch the GR quadrupole target is equivalent to one concrete statement:

> the constant wall quadrupole stiffness must equal the support softening plus the conservative Maxwell/mixed self-load plus one explicit outgoing-normalization load.

That is a much sharper theorem target than we had before.

---

## 5. Interpreting the wall stiffness on this branch

On the constant wall profile `chi_eta = u_0`, the axial gradient term vanishes. In the same modal normalization used in the earlier wall reduction, the bare quadrupole wall stiffness is therefore

`K = K_eta + 6 T_Omega`.

So the normalization test can be rewritten directly as

`K_eta + 6 T_Omega = K_req`.

This is important for the roadmap because the conservative 3PN program already isolated a separate geometry lane, and the Stage-8 gap still involved the unsolved wall stiffness. On the present branch, the higher-order quadrupole normalization target directly constrains that same geometry-side combination.

So the moving-throat PDE job has tightened again.
It is not to produce an arbitrary number called `K`.
It is to determine whether the actual quadrupolar wall stiffness of the throat satisfies the explicit algebraic relation above after support and internal gauge loading are included.

---

## 6. What is exact and what is still conditional here

### Exact in this stage

The following are exact on the chosen finite-throat branch:

- the N/N zero mode `u_0`,
- the D/N ladder `f_n`,
- the overlap law `kappa_n = sqrt(2) / ( (n+1/2) pi )`,
- the minimal-branch value `kappa = 2 sqrt(2) / pi`,
- and the fully substituted Stage-8 normalization equation.

### Still conditional

What is **not** yet fixed from first principles are the branch parameters

- `lambda_B`,
- `lambda_U`,
- `lambda_W`,
- `lambda_R`,
- `Omega_U`,
- `Omega_W`,
- `varpi`,
- and `mhat_rad`.

So this is not yet the final normalization theorem.
But it is no longer an abstract overlap problem either.
It is one explicit finite-throat branch equation in a short list of reduced parameters.

---

## 7. Best current summary after Stage 9

The overlap problem is now materially tighter.

Before this stage, the remaining gap was

> compute the radial/axial overlap integrals on the true branch.

After this stage, the first concrete finite-throat branch shows that

- the overlap problem can collapse to one exact geometric constant,
- the lowest D/N half-wave is automatically the dominant support/mixed branch,
- and the quadrupole normalization target can be rewritten as an explicit condition on the wall quadrupole stiffness.

So the next honest step is no longer to “invent more overlap algebra.”
It is to decide whether the real moving-throat branch is better approximated by this N/N–D/N axial family or whether a different wall-profile family is forced by the full PDE.

Either way, the theorem target is now much sharper than it was at the end of Stage 8.

=== moving_throat_pde_stage010_nonconstant_axial_family.md ===

# Moving-Throat PDE — Stage 10: First Nonconstant Finite-Throat Wall/Brane Family, Exact Overlap Law, and the Profile-Selection Theorem Gate

## Purpose

Stage 9 showed that on the simplest finite-throat branch — constant N/N wall and brane-like zero modes coupled to the lowest D/N half-wave — the whole radial/axial overlap problem collapses to one exact geometric constant,

`kappa_0 = 2 sqrt(2) / pi`.

That was useful, but it still left one obvious concern:

> was the Stage-9 normalization equation an artifact of choosing the completely flat axial zero mode for the wall and brane-like internal coordinate?

So the next honest step is to replace the constant profile by the **first nontrivial finite-throat profile family** and see what survives.

The cleanest way to do that is to stay in the exact N/N and D/N finite-throat bases and turn on the first N/N axial excitation. The main result is that the Stage-9 branch does survive, but it deforms in a very specific way:

- the whole branch is controlled by a single profile-dependent overlap
  `kappa(theta)`,
- the wall stiffness picks up the expected axial-gradient penalty,
- there is an exact **blind angle** where the D/N support/mixed route is shut off,
- and there is also an exact **max-coupling angle** where the wall/support overlap is stronger than on the constant branch.

So after this stage the next theorem gate is no longer “does any nonconstant profile destroy the branch?”
It is:

> which axial profile is selected by the real moving-throat eigenproblem, and does it lie near the Stage-9 branch, the max-coupling branch, or the blind-angle no-go branch?

---

## 1. First nonconstant finite-throat family

Keep the finite throat interval

`s in [0,L]`.

Use the exact N/N basis for the wall and brane-like internal coordinate,

`u_0(s) = 1 / sqrt(L)`,

`u_1(s) = sqrt(2/L) cos(pi s / L)`.

These satisfy

`-u_0'' = 0`,

`-u_1'' = (pi^2 / L^2) u_1`,

with

`u_0'(0)=u_0'(L)=0`,

`u_1'(0)=u_1'(L)=0`,

and are orthonormal on `[0,L]`.

Keep the same exact D/N half-wave branch for the trapped support and mixed channels,

`f_0(s) = sqrt(2/L) sin(pi s / (2L))`.

The minimal coherent nonconstant family is then

`chi_theta(s) = cos(theta) u_0(s) + sin(theta) u_1(s)`.

In this stage we use that same profile family for both

- the wall quadrupole axial shape, and
- the brane-like internal gauge/support coordinate,

because it is the smallest nonconstant replacement of the old “constant wall/brane zero-mode” branch that still preserves one common coherent axial deformation parameter.

So the branch choice is

`chi_eta = chi_theta`,

`u = chi_theta`,

`phi = f_0`,

`w = f_0`.

This has one immediate exact benefit:

`int_0^L ds chi_eta u = 1`.

So the direct wall/brane overlap does not degrade just because we turned on the first N/N profile correction.

---

## 2. Exact overlap law with the D/N half-wave

Define the two basic constants

`kappa_0 = int_0^L ds u_0 f_0 = 2 sqrt(2) / pi`,

`kappa_1 = int_0^L ds u_1 f_0 = -4 / (3 pi)`.

Then the coherent profile family has the exact overlap

`kappa(theta) = int_0^L ds chi_theta f_0`
`             = kappa_0 cos(theta) + kappa_1 sin(theta)`
`             = 2 ( -2 sin(theta) + 3 sqrt(2) cos(theta) ) / (3 pi)`.

So the whole nonconstant branch still collapses to one profile function, but now it is angle dependent rather than constant.

A useful amplitude form is

`kappa(theta) = rho cos(theta - theta_max)`,

with

`rho = sqrt(kappa_0^2 + kappa_1^2) = 2 sqrt(22) / (3 pi)`.

That means the family can reach a maximal overlap magnitude

`|kappa|_max = 2 sqrt(22) / (3 pi) ~= 0.994031`.

This is larger than the constant-branch value

`kappa_0 = 2 sqrt(2) / pi ~= 0.900316`.

So the flat Stage-9 branch is not the strongest-coupled member of the first nonconstant family.

### 2.1 Exact blind angle

The support/mixed D/N route shuts off when `kappa(theta)=0`, i.e.

`tan(theta_blind) = 3 sqrt(2) / 2`.

So

`theta_blind = arctan( 3 sqrt(2) / 2 )`.

At this angle the coherent N/N profile is orthogonal to the lowest D/N half-wave.

That is the first exact profile-selection no-go theorem of the moving-throat PDE program:

> if the actual wall/brane profile is driven onto the blind angle, the Stage-9 support/mixed quadrupole branch cannot realize the outgoing quadrupole normalization target.

### 2.2 Exact max-coupling angle

The positive overlap is maximized at

`tan(theta_max) = kappa_1 / kappa_0 = - sqrt(2) / 3`,

so

`theta_max = - arctan( sqrt(2) / 3 )`.

At that angle,

`sin^2(theta_max) = 2 / 11`,

and

`kappa(theta_max) = 2 sqrt(22) / (3 pi)`.

So the strongest-coupled member of the first N/N family is obtained by a small **negative** admixture of the first excited N/N mode, and the price paid is a fixed axial-gradient occupancy of `2/11`.

---

## 3. Exact wall stiffness on the nonconstant family

The linear wall operator carried from the earlier stages is

`G_eta = -T_w d_s^2 + K_eta + 6 T_Omega`

in the `l=2` quadrupole sector.

Evaluating it on the coherent family gives

`K_geo(theta) = int_0^L ds chi_theta G_eta chi_theta`
`             = K_eta + 6 T_Omega + T_w (pi^2 / L^2) sin^2(theta)`.

So the first nonconstant family changes the Stage-9 wall branch in exactly the way one would expect physically:

- the constant branch is recovered at `theta=0`,
- the first excited N/N admixture adds an axial-gradient penalty,
- and the gradient cost is quadratic in the profile angle.

At the max-coupling point,

`K_geo(theta_max) = K_eta + 6 T_Omega + 2 T_w pi^2 / (11 L^2)`.

So the first nonconstant family makes the trade-off completely explicit:

- larger D/N overlap,
- but larger axial wall stiffness.

---

## 4. Exact substitution into the Stage-8/9 minimal isotropic module

Using the same reduced couplings as before,

`C   = lambda_B int chi_eta phi`,

`G_U = lambda_U int chi_eta u`,

`G_W = lambda_W int chi_eta w`,

`R   = lambda_R int u w`,

the coherent nonconstant family gives

`C   = lambda_B kappa(theta)`,

`G_U = lambda_U`,

`G_W = lambda_W kappa(theta)`,

`R   = lambda_R kappa(theta)`.

So the Stage-8/9 quantities become

`Delta(theta) = Omega_U^2 Omega_W^2 - lambda_R^2 kappa(theta)^2`,

`Q(theta) = lambda_U^2 Omega_W^2`
`         + 2 lambda_U lambda_W lambda_R kappa(theta)^2`
`         + lambda_W^2 Omega_U^2 kappa(theta)^2`,

`P(theta) = kappa(theta) ( Omega_U^2 lambda_W + lambda_R lambda_U )`.

Therefore

`B_0(theta) = lambda_B^2 kappa(theta)^2 / varpi^2`,

`Z_0(theta) = Q(theta) / Delta(theta)`,

`N_0(theta) = kappa(theta)^2 ( Omega_U^2 lambda_W + lambda_R lambda_U )^2 / Delta(theta)^2`,

`D_0(theta) = K_geo(theta) - B_0(theta) - Z_0(theta)`.

This is the cleanest possible outcome: the whole Stage-9 branch survives with the single replacement

`kappa_0  ->  kappa(theta)`

and the single geometry lift

`K_eta + 6 T_Omega  ->  K_eta + 6 T_Omega + T_w pi^2 sin^2(theta) / L^2`.

So the Stage-9 branch was not an algebraic accident of the constant profile. It is the `theta=0` member of a whole exact finite-throat family.

---

## 5. Exact branch-level normalization condition

Substituting the nonconstant family into the Stage-8 target gives

`mhat_rad^2 kappa(theta)^2 ( Omega_U^2 lambda_W + lambda_R lambda_U )^2`
`/ [ Delta(theta) ( K_geo(theta) Delta(theta)`
`                - Delta(theta) lambda_B^2 kappa(theta)^2 / varpi^2`
`                - Q(theta) ) ]`
`= 54 G c_s^5 / (5 a^5 c^5)`.

Solving for the required geometry-side stiffness gives

`K_req(theta)`
` = lambda_B^2 kappa(theta)^2 / varpi^2`
` + Q(theta) / Delta(theta)`
` + mhat_rad^2 kappa(theta)^2 ( Omega_U^2 lambda_W + lambda_R lambda_U )^2`
`   / [ (54 G c_s^5 / (5 a^5 c^5)) Delta(theta)^2 ]`.

So the full theorem gate on this family is now

`K_eta + 6 T_Omega + T_w pi^2 sin^2(theta) / L^2 = K_req(theta)`.

This is sharper than the Stage-9 gate, not weaker.

Stage 9 asked whether the constant wall stiffness matched one algebraic target.
Stage 10 asks whether the actual wall eigenprofile chosen by the PDE lands on an angle `theta` for which the geometry-side stiffness exactly equals the profile-dressed support/mixed normalization load.

---

## 6. Exact consequences of the first nonconstant family

### 6.1 Stage 9 is recovered exactly

At `theta=0`,

`chi_theta -> u_0`,

`kappa(theta) -> kappa_0 = 2 sqrt(2) / pi`,

`K_geo(theta) -> K_eta + 6 T_Omega`.

So the entire Stage-9 branch is recovered exactly.

### 6.2 Blind-angle no-go theorem

At the blind angle,

`kappa(theta_blind)=0`.

Then

`C = G_W = R = 0`,

`P(theta_blind)=0`,

`N_0(theta_blind)=0`.

Therefore the left-hand side of the normalization equation vanishes, while the GR target on the right-hand side remains strictly positive.

So the blind-angle branch is an exact no-go for the outgoing quadrupole normalization bridge.

This is important: it means not every nonconstant wall profile is acceptable. The PDE must choose a profile that keeps nonzero projection onto the D/N support/mixed half-wave.

### 6.3 Max-coupling branch is allowed but not free

At `theta=theta_max`, the D/N overlap is maximized,

`kappa = 2 sqrt(22) / (3 pi)`.

But the wall stiffness also increases to

`K_geo = K_eta + 6 T_Omega + 2 T_w pi^2 / (11 L^2)`.

So the first nonconstant family exposes an exact trade-off:

- stronger coupling to the D/N support/mixed branch,
- but a larger geometry-side axial-gradient cost.

That is the first place where profile selection becomes a real dynamical competition rather than just an overlap calculation.

---

## 7. Best current summary after Stage 10

The first nonconstant finite-throat family does **not** destroy the Stage-9 route.
Instead it does something better: it turns Stage 9 into one member of a whole exact finite-throat profile family and isolates the first genuine profile-selection theorem gate.

The remaining question is no longer

> do nonconstant wall profiles matter?

They do.

The sharper question is now

> what axial profile does the real moving-throat eigenproblem actually select, and does that profile keep enough overlap with the D/N support/mixed half-wave to satisfy the outgoing quadrupole normalization equation?

The next honest derivation step is therefore to stop choosing the profile family by hand and solve the first actual axial wall eigenproblem with the matter/gauge loading included, so that the profile angle `theta` is no longer free but is an output of the reduced moving-throat operator itself.

=== moving_throat_pde_stage011_loaded_profile_selection.md ===

# Moving-Throat PDE — Stage 11: Loaded Axial Profile Selection, Exact 2x2 Eigenproblem, and Why the Blind-Angle Branch Is Dynamically Disfavored

## Purpose

Stage 10 replaced the constant axial wall/brane branch by the first nonconstant N/N family and showed that the whole normalization problem depends on one profile angle `theta` through

`kappa(theta) = kappa_0 cos(theta) + kappa_1 sin(theta)`.

That clarified the overlap algebra, but it still left one unsatisfactory feature:

> the profile angle `theta` was being chosen by hand rather than produced by a reduced moving-throat operator.

The next honest step is therefore to make `theta` an **output** of the first loaded axial eigenproblem.

The main result of this stage is that the minimal loaded wall operator already selects the profile direction in a very rigid way.

Using the first two N/N wall modes,

- the bare wall geometry is diagonal,
- the support/mixed branch enters as a rank-1 attractive load proportional to the overlap vector `(kappa_0, kappa_1)`,
- the profile angle satisfies one exact formula,
- and for positive loading it is driven toward the **max-coupling angle**, not toward the Stage-10 blind-angle no-go branch.

So after this stage the theorem gap narrows again:

- the blind-angle branch is still an exact algebraic no-go,
- but in the minimal positive-loading model it is also dynamically disfavored,
- and the real PDE now needs to tell us only how large the effective loading is and whether it stays in the same sign/class as this minimal model.

---

## 1. Minimal loaded wall basis

Work in the first two exact N/N wall modes,

`u_0(s) = 1 / sqrt(L)`,

`u_1(s) = sqrt(2/L) cos(pi s / L)`.

Write the wall quadrupole profile as

`chi(s) = q_0 u_0(s) + q_1 u_1(s)`,

with normalized coefficient vector

`q = (q_0, q_1)^T`,

`q_0^2 + q_1^2 = 1`.

The natural angle parameterization is

`q_0 = cos(theta)`,

`q_1 = sin(theta)`.

The bare wall operator carried from the earlier stages is

`G_eta = -T_w d_s^2 + K_eta + 6 T_Omega`.

So its matrix in the `{u_0,u_1}` basis is exactly diagonal,

`K_bare = [[K_0, 0], [0, K_1]]`,

with

`K_0 = K_eta + 6 T_Omega`,

`K_1 = K_eta + 6 T_Omega + Delta K_ax`,

`Delta K_ax = T_w pi^2 / L^2`.

The D/N half-wave coupling vector from Stage 10 is

`v = (kappa_0, kappa_1)^T`,

with

`kappa_0 = 2 sqrt(2) / pi`,

`kappa_1 = -4 / (3 pi)`.

So the Stage-10 profile overlap is simply

`kappa(theta) = v . q`.

---

## 2. Minimal effective loading from the support/mixed branch

At low frequency, the support/mixed branch lowers the wall energy in proportion to the square of the overlap with the D/N half-wave. The minimal reduced model is therefore

`E_eff[q] = (1/2) q^T K_bare q - (alpha/2) (v . q)^2`,

with one positive loading susceptibility

`alpha > 0`.

In matrix form,

`K_eff = K_bare - alpha v v^T`.

So explicitly,

`K_eff = [[K_0 - alpha kappa_0^2,      -alpha kappa_0 kappa_1],`
`         [    -alpha kappa_0 kappa_1,  K_1 - alpha kappa_1^2 ]].`

This is the first actual loaded moving-throat axial operator in the program.

It is still reduced and minimal, but it already converts the free profile angle `theta` into a dynamical eigenvector problem.

---

## 3. Exact eigenvalues and exact profile-angle equation

The effective trace and determinant are

`tr(K_eff) = K_0 + K_1 - alpha (kappa_0^2 + kappa_1^2)`,

`det(K_eff) = K_0 K_1 - alpha ( K_1 kappa_0^2 + K_0 kappa_1^2 )`.

So the exact eigenvalues are

`lambda_(+-)`
` = (1/2) [ K_0 + K_1 - alpha (kappa_0^2 + kappa_1^2)`
`          +- sqrt( (Delta K_ax + alpha (kappa_0^2 - kappa_1^2))^2`
`                   + 4 alpha^2 kappa_0^2 kappa_1^2 ) ]`.

The lower branch `lambda_-` is the physically relevant loaded wall stiffness.

Parameterizing the eigenvector by `theta`, the exact stationarity condition is

`tan(2 theta)`
` = 2 alpha kappa_0 kappa_1 / ( Delta K_ax + alpha (kappa_0^2 - kappa_1^2) )`.

Since

`kappa_0^2 - kappa_1^2 = 56 / (9 pi^2) > 0`,

and

`2 kappa_0 kappa_1 = -16 sqrt(2) / (3 pi^2) < 0`,

the sign structure is rigid.

For any positive loading `alpha > 0`, the denominator is positive and the numerator is negative, so

`theta < 0`.

This already proves an important selection statement:

> the first loaded wall eigenprofile is driven away from the flat Stage-9 branch in the **negative** `u_1` direction.

That is exactly the direction of the Stage-10 max-coupling angle, not the positive blind-angle branch.

---

## 4. Small-loading and strong-loading limits

### 4.1 Weak loading

For small loading,

`theta = alpha kappa_0 kappa_1 / Delta K_ax + O(alpha^2)`.

Because `kappa_0 kappa_1 < 0`, the first correction to the constant branch is a small negative angle.

So the Stage-9 flat branch is perturbatively stable, but it is not the exact loaded eigenprofile once the support/mixed load is turned on.

### 4.2 Strong loading

For large positive loading,

`tan(2 theta) -> 2 kappa_0 kappa_1 / (kappa_0^2 - kappa_1^2)`.

That is precisely the two-angle relation satisfied by the normalized coupling vector `v / |v|`.

So in the strong-loading limit the selected profile approaches

`q -> v / |v|`,

which means

`theta -> theta_max`,

with

`tan(theta_max) = kappa_1 / kappa_0 = - sqrt(2) / 3`.

This is a major simplification.

The Stage-10 max-coupling angle is **not** an arbitrary hand-picked member of the family. It is the exact strong-loading eigenvector of the first reduced loaded wall operator.

So the minimal loaded model already explains why the physically preferred profile should move toward the max-coupling branch rather than stay flat forever.

---

## 5. Why the blind-angle no-go branch is dynamically disfavored

Stage 10 showed that the blind angle

`tan(theta_blind) = 3 sqrt(2) / 2`

forces

`kappa(theta_blind)=0`,

and therefore kills the outgoing quadrupole normalization bridge.

Stage 11 sharpens that conclusion.

Because the loaded wall angle satisfies

`tan(2 theta)`
` = 2 alpha kappa_0 kappa_1 / ( Delta K_ax + alpha (kappa_0^2 - kappa_1^2) )`,

and the right-hand side is always **negative** for positive loading,

the selected `theta` always lies in the negative-coupling direction.

But the blind angle is positive,

`theta_blind > 0`.

So in the minimal attractive-loading model,

> the blind-angle branch is not only an algebraic no-go for normalization; it is dynamically disfavored by the profile-selection eigenproblem itself.

That is the strongest theorem-like statement we have so far about profile selection.

---

## 6. Exact softening threshold

The lower eigenvalue crosses zero when

`det(K_eff)=0`.

So the exact loading threshold is

`alpha_crit = K_0 K_1 / ( K_1 kappa_0^2 + K_0 kappa_1^2 )`.

This is the first concrete softening threshold in the moving-throat PDE program.

Interpretation:

- for `alpha < alpha_crit`, the loaded quadrupole wall mode remains stable,
- for `alpha = alpha_crit`, the branch goes marginal,
- for `alpha > alpha_crit`, the minimal reduced wall operator predicts a quadrupole instability/condensation.

So the future PDE does not just need to say which angle is chosen.
It must also place the physical branch on the stable or near-softened side of this exact threshold.

---

## 7. What this means for the normalization bridge

The outgoing-normalization equation from Stage 10 depends on

- the geometry-side wall stiffness,
- and the profile-dressed overlap `kappa(theta)`.

Stage 11 shows that these are not independent in the minimal loaded model.
They are linked because the same support/mixed loading that wants a larger `|kappa|` also rotates the wall toward `theta_max` and softens the lower wall eigenvalue.

So the next theorem question is no longer

> can I pick a `theta` that makes the algebra work?

The sharper question is

> does the actual loaded wall eigenmode selected by the throat operator carry the right pair `(lambda_-, kappa(theta))` to satisfy the outgoing quadrupole normalization target while remaining on the stable branch?

That is a much more physical and much more falsifiable question.

---

## 8. Best current summary after Stage 11

The first reduced loaded wall eigenproblem already resolves a major ambiguity.

It shows that:

- the Stage-10 profile angle is not a free nuisance parameter,
- the support/mixed load selects a definite direction in the N/N wall basis,
- the selected profile rotates away from the constant branch toward the max-coupling branch,
- the blind-angle no-go branch is dynamically disfavored rather than merely algebraically bad,
- and there is an exact softening threshold `alpha_crit` for the loaded wall quadrupole mode.

So the next honest derivation step is now very sharply defined:

> derive the actual effective loading strength `alpha` and the associated loaded wall eigenmode from the coupled wall/BdG/Maxwell/mixed operator, then insert that eigenpair into the Stage-10 normalization equation.

At this point the roadmap is no longer “invent a PDE.”
It is starting to look like a specific spectral problem.

=== moving_throat_pde_stage012_dynamic_loading.md ===

# Moving-Throat PDE — Stage 12: Exact Dynamic Loading from the Coupled Wall/BdG/Maxwell/Mixed Operator, and the First Selected-Mode Quadrupole Damping Coefficient

## Purpose

Stage 11 deliberately introduced a reduced loading parameter `alpha` to turn the free profile angle `theta` into a real eigenproblem. That was the right move at the time, but it also left the central spectral question unfinished:

> what does `alpha` actually come from in the coupled wall/BdG/Maxwell/mixed operator, and how does the same operator feed the outgoing `l=2` odd term into the *selected* wall quadrupole mode?

This stage answers that question for the first nontrivial reduced operator that is still microscopically faithful to the project ontology.

The main result is that the coupled operator is much cleaner than it looked before.
Using the first two N/N wall modes, the lowest D/N support/mixed half-wave, and a brane-like internal zero-mode doublet, the full reduced wall self-energy splits **exactly** into two pieces:

- an isotropic shift `Xi(omega) I_2`,
- and a rank-1 directional load `alpha(omega) v v^T`.

So Stage 11 was not just a plausible toy model. Its rank-1 loading structure is the exact Schur complement of the first coupled wall/BdG/Maxwell/mixed block.

Even better, once the mixed channel is dressed by the passive/outgoing port, the same exact decomposition survives and the first odd coefficient of the selected wall mode becomes explicit.
So the bottleneck narrows again:

- the profile-selection problem is no longer carried by a free parameter,
- the conservative loading is now an exact rational function of the reduced couplings,
- and the outgoing `i omega^5` coefficient of the selected wall mode is now a one-line projection formula.

---

## 1. Minimal coupled axial operator in the `{u_0,u_1}` wall basis

Work in the first two exact N/N wall modes,

`u_0(s) = 1 / sqrt(L)`,

`u_1(s) = sqrt(2/L) cos(pi s / L)`,

and write the wall quadrupole shape as the coefficient vector

`q = (q_0, q_1)^T`.

The D/N overlap vector from Stages 10–11 is

`v = (kappa_0, kappa_1)^T`,

with

`kappa_0 = 2 sqrt(2) / pi`,

`kappa_1 = -4 / (3 pi)`,

`sigma = v.v = kappa_0^2 + kappa_1^2 = 88 / (9 pi^2)`.

So `sigma` is exactly the Stage-10 max-coupling value `|kappa|_max^2`.

Let the bare wall operator in this basis be

`D_eta^(bare)(omega) = diag( K_0 - M_0 omega^2, K_1 - M_1 omega^2 )`,

with

`K_0 = K_eta + 6 T_Omega`,

`K_1 = K_0 + Delta K_ax`,

`Delta K_ax = T_w pi^2 / L^2`.

Now keep the smallest internal block that still matches the earlier reduced hierarchy:

- one BdG support mode `phi` on the D/N half-wave,
- one brane-like internal doublet `u = (u_0^{int}, u_1^{int})^T` in the same N/N basis,
- one mixed `A_w/F_(mu w)/J^w` mode `W` on the same D/N half-wave.

Use the reduced frequency-domain kernels

`A_phi(omega) = varpi^2 - omega^2`,

`A_U(omega)   = Omega_U^2 - omega^2`,

`A_W(omega)   = Omega_W^2 - omega^2 - Pi_out(omega)`.

The exact reduced couplings are taken to be

- wall–BdG: `lambda_B (v.q) phi`,
- wall–U:   `lambda_U q.u`,
- wall–W:   `lambda_W (v.q) W`,
- U–W:      `- lambda_R (v.u) W`.

The last sign is the one compatible with the earlier conservative self-energy convention in which the mixed-sector cross term enters with the *positive* `+ 2 G_U G_W R` numerator after elimination.

This reduced operator is the minimal dynamic refinement of Stage 11:

- `q` is the wall profile we want to select,
- `phi` is the BdG support lane already turned on at Stage 3,
- `u` is the lowest brane-like internal gauge lane,
- `W` is the first mixed `A_w/F_(mu w)/J^w` lane,
- and `Pi_out` is the passive/outgoing dressing that can carry the `i omega^5` branch.

---

## 2. Exact Schur-complement decomposition of the wall self-energy

Eliminating the internal fields `(u, W, phi)` gives an exact effective wall operator of the form

`D_wall(omega) = D_eta^(bare)(omega) - Sigma_wall(omega)`.

The first key theorem of this stage is that the wall self-energy splits **exactly** as

`Sigma_wall(omega) = Xi(omega) I_2 + alpha(omega) v v^T`,

where

`Xi(omega) = lambda_U^2 / A_U(omega)`

and

`alpha(omega) = lambda_B^2 / A_phi(omega)`
`             + ( A_U(omega) lambda_W + lambda_R lambda_U )^2`
`               / [ A_U(omega) Delta_UW(omega) ]`.

Here

`Delta_UW(omega) = A_U(omega) A_W(omega) - lambda_R^2 sigma`.

So the full coupled wall/BdG/Maxwell/mixed operator really does reduce to exactly the Stage-11 geometry:

`D_wall(omega) = [ D_eta^(bare)(omega) - Xi(omega) I_2 ] - alpha(omega) v v^T`.

This is a stronger result than the old Stage-11 ansatz because it identifies the two physically distinct effects:

1. `Xi(omega)` is an **isotropic internal zero-mode shift** coming from the brane-like `U` doublet,
2. `alpha(omega)` is the **directional loading strength** along the D/N overlap vector `v`.

So the first loaded profile-selection problem is no longer phenomenological.
It is the exact first Schur complement of the reduced wall/BdG/Maxwell/mixed operator.

---

## 3. Conservative static loading and the refined Stage-11 profile-selection law

On the conservative static branch set

`Pi_out(0) = 0`.

Then

`Xi_0 = lambda_U^2 / Omega_U^2`,

`Delta_0 = Omega_U^2 Omega_W^2 - lambda_R^2 sigma`,

`alpha_0 = lambda_B^2 / varpi^2`
`         + ( Omega_U^2 lambda_W + lambda_R lambda_U )^2 / ( Omega_U^2 Delta_0 )`.

So the conservative loaded wall matrix is

`K_eff^(0) = [[K_0 - Xi_0, 0], [0, K_1 - Xi_0]] - alpha_0 v v^T`.

The isotropic shift moves both wall levels equally,
so the profile-angle equation keeps the same Stage-11 structure with

`K_0 -> K_0 - Xi_0`,

`K_1 -> K_1 - Xi_0`,

`alpha -> alpha_0`.

In particular, the exact conservative angle law remains

`tan(2 theta_-)`
`= 2 alpha_0 kappa_0 kappa_1 / ( Delta K_ax + alpha_0 (kappa_0^2 - kappa_1^2) )`.

So the same sign theorem survives:
for positive `alpha_0`, the selected lower eigenvector still rotates in the negative `u_1` direction, i.e. toward the Stage-10 max-coupling branch and away from the blind-angle branch.

### 3.1 Refined conservative softening threshold

Because the internal `U` doublet contributes the isotropic shift `Xi_0`, the exact softening threshold is refined to

`alpha_crit`
`= (K_0 - Xi_0)(K_1 - Xi_0)`
`  / [ (K_1 - Xi_0) kappa_0^2 + (K_0 - Xi_0) kappa_1^2 ]`.

So the internal zero-mode sector makes quadrupole softening easier in a very concrete sense:

- it does **not** change the direction-selection formula,
- but it lowers the stability margin by shifting both bare wall levels downward before the directional rank-1 loading is applied.

This is the first exact place where the BdG/Maxwell/mixed operator changes the Stage-11 theorem content rather than merely justifying it.

---

## 4. Exact outgoing dressing of the dynamic loading strength

Now restore the passive/outgoing dressing in the mixed channel,

`A_W(omega) = Omega_W^2 - omega^2 - Pi_out(omega)`.

The exact decomposition survives unchanged:

`Sigma_wall(omega) = Xi(omega) I_2 + alpha(omega) v v^T`,

with the same closed formula for `alpha(omega)` but with the dressed `Delta_UW(omega)`.

Expanding to first order in the outgoing port gives

`alpha(omega) = alpha_cons(omega) + beta(omega) Pi_out(omega) + O(Pi_out^2)`

with the exact transfer factor

`beta(omega) = [ A_U(omega) lambda_W + lambda_R lambda_U ]^2 / Delta_cons(omega)^2`,

where

`Delta_cons(omega) = A_U(omega) (Omega_W^2 - omega^2) - lambda_R^2 sigma`.

So the directional loading strength inherits the outgoing branch with the same kind of positive transfer factor already seen in the one-lane Stage-4 mixed-sector analysis.

At zero frequency,

`beta_0 = ( Omega_U^2 lambda_W + lambda_R lambda_U )^2 / Delta_0^2 >= 0`.

This is the second key theorem of the stage:

> once the passive/outgoing port is attached to the mixed channel, the *same* wall-profile loading parameter that selects the quadrupole eigenmode also inherits the odd `i omega^5` branch with a manifestly nonnegative conservative transfer factor.

---

## 5. First selected-mode odd quadrupole coefficient

On the natural compact outgoing `l=2` branch,

`Pi_out(omega) = + i Gamma_2^(port) omega^5 + O(omega^7)`

with

`Gamma_2^(port) = a^5 / (27 c_s^5)`

from the earlier outgoing Hankel/D/N audit.

Therefore the directional loading strength has the low-frequency form

`alpha(omega) = alpha_even(omega) + i beta_5 omega^5 + O(omega^7)`

with

`beta_5 = Gamma_2^(port) beta_0`
`      = Gamma_2^(port) ( Omega_U^2 lambda_W + lambda_R lambda_U )^2 / Delta_0^2`.

Equivalently,

`beta_5 = [ a^5 / (27 c_s^5) ]`
`        * ( Omega_U^2 lambda_W + lambda_R lambda_U )^2 / ( Omega_U^2 Omega_W^2 - lambda_R^2 sigma )^2`.

Now project this odd piece onto the *conservative* lower eigenvector `e_-` of `K_eff^(0)`.
Since the odd operator is rank-1,

`delta D_-^(odd)(omega)`
`= - i beta_5 ( v.e_- )^2 omega^5 + O(omega^7)`

in the wall-operator convention used throughout these notes.

So the selected-mode damping strength is controlled by exactly two ingredients:

1. the dynamic transfer coefficient `beta_5`,
2. the selected profile overlap factor `(v.e_-)^2`.

The second factor is just the Stage-10 overlap evaluated on the Stage-11 selected eigenmode. In angle language,

`(v.e_-)^2 = kappa(theta_-)^2`.

So the operator-selected odd coefficient is

`delta D_-^(odd)(omega)`
`= - i beta_5 kappa(theta_-)^2 omega^5 + O(omega^7)`.

This is the first exact formula in the moving-throat PDE program that feeds the passive/outgoing `l=2` odd term directly into the **selected** wall quadrupole mode rather than into a hand-chosen profile.

---

## 6. Relation to the earlier Stage-4/5 transfer factors

The conservative combination

`P_0 = Omega_U^2 lambda_W + lambda_R lambda_U`

is the same `P`-type numerator that already appeared in the earlier isotropic normalization work.
So the Stage-12 directional-loading transfer factor

`beta_0 = P_0^2 / Delta_0^2`

is not a new unrelated object.
It is the same mixed-sector transfer numerator that the Stage-4/5 bridge had already isolated, now appearing on the first real wall-profile selection operator.

That is an important unification point:

- Stage 4 first isolated the outgoing mixed-sector transfer factor on one lane,
- Stage 8 reduced the normalization problem to the radial/axial quantities `P`, `Delta`, `Q`, `X`,
- Stage 11 turned the wall profile into a real eigenproblem,
- and Stage 12 shows that the same mixed-sector numerator controls the directional loading and the selected-mode odd term.

So the different branches of the roadmap are starting to meet.

---

## 7. Best current summary after Stage 12

The bottleneck has narrowed again.

The coupled wall/BdG/Maxwell/mixed operator now tells us that the first loaded profile-selection problem is governed by two exact reduced objects:

- an isotropic shift
  `Xi_0 = lambda_U^2 / Omega_U^2`,
- and a directional load
  `alpha_0 = lambda_B^2 / varpi^2 + ( Omega_U^2 lambda_W + lambda_R lambda_U )^2 / ( Omega_U^2 Delta_0 )`.

The outgoing odd term of the selected lower quadrupole mode is then

`delta D_-^(odd)(omega)`
`= - i [ a^5 / (27 c_s^5) ]`
`    * [ ( Omega_U^2 lambda_W + lambda_R lambda_U )^2 / Delta_0^2 ]`
`    * kappa(theta_-)^2 omega^5 + O(omega^7)`.

So the Stage-11 “free parameter” bottleneck is gone.
What remains is now much sharper:

1. compute the conservative branch data `(Xi_0, alpha_0)` on the true moving-throat branch,
2. check that `alpha_0 < alpha_crit` so the selected quadrupole wall mode stays on the stable side,
3. compute the selected overlap `kappa(theta_-)^2`,
4. and then test whether the resulting operator-selected odd coefficient lands on the required normalization branch.

That is a much smaller theorem gap than we had before.
It is no longer “derive some PDE loading somehow.”
It is a specific spectral-transfer problem.

=== moving_throat_pde_stage013_selected_mode_normalization.md ===

# Moving-Throat PDE — Stage 13: Selected-Mode Normalized Response, Exact Static Prefactor, and the Selected-Branch Quadrupole Target

## Purpose

Stage 12 reached the first exact operator-selected odd quadrupole coefficient,

`delta D_-^(odd)(omega) = - i beta_5 (v.e_-)^2 omega^5 + O(omega^7)`.

That was the right operator-level result, but it still stopped one step short of the language used by the earlier grouped-real `P2` normalization bridge.

The next missing step is therefore very specific:

> translate the selected lower wall mode into the same normalized-response convention used in Stages 5 and 8, identify the exact selected-branch static prefactor, and write the surviving `2G/(5 c^5)` test as a single spectral condition on the selected conservative eigenvalue.

This stage does exactly that.

The main result is simple and important.
If the conservative selected lower wall eigenvalue is written as `lambda_-(omega)`, then the normalized selected-mode response has the low-frequency form

`Y_-(omega) = 1 + u_{2,-} omega^2 + u_{4,-} omega^4 + i Gamma_{5,-} omega^5 + ...`,

with

`Gamma_{5,-} = Gamma_2^(port) * P_{0,-}`,

where the exact selected static prefactor is

`P_{0,-} = beta_0 (v.e_-)^2 / lambda_-(0)`

and, equivalently,

`P_{0,-} = - beta_0 d(ln lambda_-)/d alpha_0`.

So the selected-branch normalization problem is no longer “some outgoing coefficient somewhere.”
It is one exact ratio of a conservative overlap factor to the conservative selected wall stiffness.

---

## 1. Selected conservative eigenvalue and notation

Keep the Stage-12 conservative wall matrix

`K_eff^(0) = [[A, 0], [0, B]] - alpha_0 v v^T`,

with

`A = K_0 - Xi_0`,

`B = K_1 - Xi_0 = A + Delta K_ax`,

`v = (kappa_0, kappa_1)^T`,

`sigma = kappa_0^2 + kappa_1^2`,

`delta_kappa = kappa_0^2 - kappa_1^2`,

`KappaProd = kappa_0^2 kappa_1^2`.

The exact selected lower conservative eigenvalue is

`lambda_-`
`= ( A + B - alpha_0 sigma - R ) / 2`,

where

`R = sqrt( (Delta K_ax + alpha_0 delta_kappa)^2 + 4 alpha_0^2 KappaProd )`.

The corresponding upper branch is

`lambda_+ = ( A + B - alpha_0 sigma + R ) / 2`.

Throughout this stage,

`D_{-0} = lambda_-`

is the selected conservative wall stiffness at zero frequency.

---

## 2. Exact selected overlap from Hellmann–Feynman

Stage 12 already established the exact Hellmann–Feynman identity

`(v.e_-)^2 = - d lambda_- / d alpha_0`.

Writing

`s_- = (v.e_-)^2`,

the explicit closed form is

`s_-`
`= (1/2) [ sigma`
`          + ( (Delta K_ax + alpha_0 delta_kappa) delta_kappa + 4 alpha_0 KappaProd ) / R ]`.

This formula immediately passes the expected checks:

- weak loading `alpha_0 -> 0`:
  `s_- -> kappa_0^2`,
- strong loading `alpha_0 -> +infinity`:
  `s_- -> sigma = |v|^2`.

So the selected overlap interpolates smoothly from the flat branch to the Stage-10 max-coupling branch.

---

## 3. From the selected operator to the normalized selected response

Write the selected wall operator in the same low-frequency form used in Stage 5,

`D_-(omega)`
`= D_{-0} + D_{-2} omega^2 + D_{-4} omega^4 - i C_{5,-} omega^5 + O(omega^6)`.

Here the Stage-12 odd coefficient is

`C_{5,-} = beta_5 s_-`.

The normalized selected-mode response is therefore

`Y_-(omega) = D_{-0} / D_-(omega)`
`          = 1 + u_{2,-} omega^2 + u_{4,-} omega^4 + i Gamma_{5,-} omega^5 + O(omega^6)`.

Exactly as in Stage 5, the even coefficients are

`u_{2,-} = - D_{-2} / D_{-0}`,

`u_{4,-} = ( D_{-2}^2 - D_{-0} D_{-4} ) / D_{-0}^2`.

The selected odd coefficient is

`Gamma_{5,-} = C_{5,-} / D_{-0}`
`           = beta_5 s_- / lambda_-`.

So the Stage-12 operator-level result is now fully translated into the Stage-5 normalized-response language.

---

## 4. Exact selected static prefactor `P_{0,-}`

Using the compact outgoing fingerprint

`beta_5 = Gamma_2^(port) beta_0`,

`Gamma_2^(port) = a^5 / (27 c_s^5)`,

we obtain

`Gamma_{5,-} = Gamma_2^(port) P_{0,-}`

with the exact selected static prefactor

`P_{0,-} = beta_0 s_- / lambda_-`.

This is the selected-mode analogue of the isotropic Stage-5 prefactor `P_0 = N_0 / D_0`.

Using the Hellmann–Feynman identity, the same prefactor can be rewritten as

`P_{0,-} = - beta_0 d(ln lambda_-)/d alpha_0`.

This is a very useful compression of the theorem gap.
It says the passive/outgoing quadrupole normalization is controlled by how sensitively the conservative selected wall stiffness responds to the directional loading parameter.

---

## 5. The selected-branch `2G/(5 c^5)` target

The invariant 2.5PN quadrupole condition is still

`mhat_-^2 Gamma_{5,-} = 2 G / (5 c^5)`.

Substituting the selected prefactor gives

`mhat_-^2 * Gamma_2^(port) * P_{0,-} = 2 G / (5 c^5)`.

Using

`Gamma_2^(port) = a^5 / (27 c_s^5)`

this is equivalent to

`mhat_-^2 P_{0,-} = 54 G c_s^5 / (5 a^5 c^5)`.

So the selected branch is required to hit exactly the same normalization stack as the isotropic Stage-5/8 branch, but now with the selected-mode prefactor in place of the old isotropic lane prefactor.

Equivalently, the target becomes a direct conservative spectral condition:

`lambda_- = mhat_-^2 beta_0 s_- / N_Q^(target)`,

where

`N_Q^(target) = 54 G c_s^5 / (5 a^5 c^5)`.

On the natural source-map branch,

`mhat_- = 1 + O(a^2/r^2)`,

so at leading point-particle order,

`lambda_- = beta_0 s_- * 5 a^5 c^5 / (54 G c_s^5)`.

This is the selected-mode generalization of the Stage-9 stiffness test.
The difference is that the “stiffness” is no longer a hand-picked wall constant.
It is the exact selected conservative eigenvalue of the loaded wall/BdG/Maxwell/mixed operator.

---

## 6. Exact determinant identity and a useful spectral rewrite

Because the selected mode comes from a rank-1 update of a diagonal `2 x 2` wall block, its determinant factor is exact:

`lambda_- lambda_+ = A B - alpha_0 ( B kappa_0^2 + A kappa_1^2 )`.

So the exact softening threshold is

`alpha_crit = A B / ( B kappa_0^2 + A kappa_1^2 )`,

which is just the Stage-12 refined threshold written in the compact `A,B` notation.

This means the selected-branch quadrupole target can also be read as a condition on **how close** the physical branch sits to the softening surface.
If `lambda_-` is too large, the selected outgoing coefficient is too small.
If `lambda_-` is too small, the branch is close to instability.

So the remaining PDE question is now sharply spectral:

> what conservative selected eigenvalue `lambda_-` and selected overlap `s_-` does the physical moving-throat branch produce, and does their ratio land on the universal target before the wall mode softens?

---

## 7. Best current summary after Stage 13

The theorem gap has narrowed again.

We no longer need to talk vaguely about an operator-level odd coefficient.
The selected-branch 2.5PN problem is now the single exact quantity

`P_{0,-} = beta_0 (v.e_-)^2 / lambda_-`

or, equivalently,

`P_{0,-} = - beta_0 d(ln lambda_-)/d alpha_0`.

And the full selected-branch normalization theorem is simply

`mhat_-^2 P_{0,-} = 54 G c_s^5 / (5 a^5 c^5)`.

That is the cleanest selected-mode formulation we have reached so far.
The next honest step is no longer to invent another reduced parameter.
It is to determine whether the physical stable branch can actually reach this target in a controlled way.

=== moving_throat_pde_stage014_selected_branch_reachability.md ===

# Moving-Throat PDE — Stage 14: Selected-Branch Reachability, Monotonicity, and the Stable-Side Normalization Window

## Purpose

Stage 13 translated the selected lower quadrupole mode into the same normalized-response language used by the earlier grouped-real `P2` work and reduced the remaining theorem gap to one exact spectral ratio,

`P_{0,-} = beta_0 (v.e_-)^2 / lambda_-`.

The next question is now unavoidable:

> on the stable branch, is this selected prefactor a wildly oscillatory function of the loading, or is the normalization target a clean one-parameter crossing problem?

The answer is much better than it had any right to be.
For the exact `2 x 2` loaded wall problem, the selected prefactor is **strictly monotone increasing** on the stable side.
It starts at the flat-branch value, diverges at the softening threshold, and therefore hits any target above its starting value exactly once before instability.

So the normalization question is no longer an arbitrary tuning problem.
It is a stable-side spectral crossing theorem.

---

## 1. Exact selected overlap derivative

Keep the Stage-13 notation

`A = K_0 - Xi_0`,

`B = K_1 - Xi_0 = A + Delta K_ax`,

`s_- = (v.e_-)^2`,

`lambda_- = ( A + B - alpha_0 sigma - R ) / 2`,

`R = sqrt( (Delta K_ax + alpha_0 delta_kappa)^2 + 4 alpha_0^2 KappaProd )`.

The first new exact identity is

`d s_- / d alpha_0 = 2 Delta K_ax^2 KappaProd / R^3`.

Because

`Delta K_ax > 0`,

`KappaProd = kappa_0^2 kappa_1^2 > 0`,

and `R > 0`,

this derivative is strictly positive on the whole stable branch.

So the selected overlap with the loading vector does not meander.
It grows monotonically as the directional loading is increased.

---

## 2. Exact monotonicity of the selected static prefactor

The selected prefactor is

`P_{0,-} = beta_0 s_- / lambda_-`.

Using the Hellmann–Feynman identity

`d lambda_- / d alpha_0 = - s_-`,

its derivative is exactly

`d P_{0,-} / d alpha_0`
`= beta_0 [ (d s_- / d alpha_0) lambda_- + s_-^2 ] / lambda_-^2`.

On the stable branch,

- `beta_0 > 0`,
- `lambda_- > 0`,
- `d s_- / d alpha_0 > 0`,
- `s_-^2 > 0`.

Therefore

`d P_{0,-} / d alpha_0 > 0`.

So the selected static prefactor is strictly monotone increasing all the way up to the softening threshold.

This is the core theorem of the stage.
It turns the selected-branch normalization problem from an opaque reduced-parameter search into a one-dimensional ordered crossing problem.

---

## 3. Starting value at zero loading

At `alpha_0 = 0`, the lower wall mode is just the flat `u_0` branch, so

`lambda_-(0) = A = K_0 - Xi_0`,

`s_-(0) = kappa_0^2`.

Therefore the selected prefactor starts at

`P_{0,-}(0) = beta_0 kappa_0^2 / (K_0 - Xi_0)`.

This is the exact stable-side starting value against which the universal target must be compared.

---

## 4. Exact softening threshold and divergence of `P_{0,-}`

The exact determinant identity is

`lambda_- lambda_+ = A B - alpha_0 ( B kappa_0^2 + A kappa_1^2 )`.

So the refined softening threshold is

`alpha_crit = A B / ( B kappa_0^2 + A kappa_1^2 )`.

At that point,

`lambda_-(alpha_crit) = 0`,

while the selected overlap `s_-` stays finite and positive.

Because

`P_{0,-} = beta_0 s_- / lambda_-`,

it follows that

`P_{0,-} -> +infinity`

as `alpha_0 -> alpha_crit^-` from the stable side.

So the stable branch spans the full interval

`P_{0,-}(0) <= P_{0,-} < +infinity`.

---

## 5. Unique stable-side crossing theorem

Combine the last two sections.

- `P_{0,-}(alpha_0)` is continuous on `0 <= alpha_0 < alpha_crit`,
- strictly increasing on that interval,
- starts at `P_{0,-}(0)`,
- and diverges at `alpha_crit^-`.

Therefore:

> for every target value `P_target > P_{0,-}(0)`, there exists a **unique** stable-side loading `alpha_*` with `0 < alpha_* < alpha_crit` such that `P_{0,-}(alpha_*) = P_target`.

Applied to the 2.5PN normalization target,

`P_target = N_Q^(target) / mhat_-^2`

with

`N_Q^(target) = 54 G c_s^5 / (5 a^5 c^5)`,

this means the selected-branch quadrupole condition is not a discrete impossibility.
It is a unique crossing problem on the stable branch.

The only remaining question is whether the physical moving-throat branch puts the system in the right region of parameter space so that the crossing occurs on the natural passive/outgoing branch before other approximations fail.

---

## 6. Best current theorem gate after Stage 14

After Stage 14, the normalization bottleneck can be stated in its sharpest current form.

The selected-branch 2.5PN problem is no longer:

- “derive some outgoing coefficient,”
- or “guess the right axial profile,”
- or “hope the branch lands near the right value.”

It is now the following exact test:

1. compute the physical stable-branch data `(Xi_0, beta_0, alpha_0)` from the coupled moving-throat operator,
2. evaluate the selected prefactor
   `P_{0,-}(alpha_0) = beta_0 (v.e_-)^2 / lambda_-`,
3. compare it to the target `N_Q^(target) / mhat_-^2`,
4. and check whether the resulting stable-side crossing sits on the natural passive/outgoing branch with `alpha_0 < alpha_crit`.

So the theorem gap is now a controlled spectral-placement problem, not an algebraic unknown.

=== moving_throat_pde_stage015_source_map_from_mode_integrals.md ===

# Moving-Throat PDE — Stage 15: Explicit Finite-Throat Mode Integrals, Kernel-Level Couplings, and Elimination of the Abstract Selected-Branch Source Map

## Purpose

After Stage 14 the normalization bottleneck had been reduced to the selected-branch quantity

`mhat_-^2 P_{0,-}`,

but the source-map factor `mhat_-` was still being carried as an abstract branch datum.
That was honest, but it also meant one important part of the selected-branch normalization problem had not yet been tied back to an explicit finite-throat mode model.

The next natural question is therefore:

> in the first explicit finite-throat axial source model, does the same D/N overlap structure that generates the mixed-sector loading also determine the selected source map, so that `mhat_-` stops being an independent parameter?

This stage answers yes.
In the first local isotropic kernel model built on the exact N/N wall basis and the exact D/N half-wave, every relevant wall–internal and source–wall coupling is determined by the same overlap vector

`v = (kappa_0, kappa_1)^T`.

As a result, on the natural D/N source branch the selected-mode source map is

`mhat_- = (v.e_-)/kappa_0`,

so

`mhat_-^2 = s_- / kappa_0^2`,

with `s_- = (v.e_-)^2`.

This removes the abstract source-map factor from the selected-branch quadrupole target.
The full selected-branch normalization product becomes

`mhat_-^2 P_{0,-} = beta_0 s_-^2 / (kappa_0^2 lambda_-)`,

which is now completely determined by the conservative selected-mode spectral data and the outgoing transfer factor.

There is also a useful side result:
under this natural axial source model the source-map amplification is modest and exact,

`1 <= mhat_-^2 <= 11/9`.

So the normalization burden is not being carried by a large hidden source renormalization.
It is being carried mainly by the selected stiffness `lambda_-` and the mixed-sector transfer factor `beta_0`.

---

## 1. Exact finite-throat axial basis and the overlap vector `v`

Keep the same finite-throat axial interval `s in [0,L]` and the same wall basis used from Stage 10 onward.
The first two exact N/N wall modes are

`u_0(s) = 1 / sqrt(L)`,

`u_1(s) = sqrt(2/L) cos(pi s / L)`.

The natural compact outgoing/internal half-wave is the exact D/N mode

`f_0(s) = sqrt(2/L) sin(pi s / (2L))`.

The wall-to-D/N overlap vector is therefore

`v_i = int_0^L u_i(s) f_0(s) ds`.

The exact values are

`kappa_0 = 2 sqrt(2) / pi`,

`kappa_1 = - 4 / (3 pi)`,

`v = (kappa_0, kappa_1)^T`.

So the squared norm is

`sigma = v.v = kappa_0^2 + kappa_1^2 = 88 / (9 pi^2)`.

It is also useful to record the flat-branch normalization constant

`kappa_0^2 = 8 / pi^2`,

and therefore the exact ratio

`sigma / kappa_0^2 = 11 / 9`.

That number will become the maximal source-map enhancement on the selected branch.

---

## 2. Local isotropic kernel model and exact reduced couplings

Now make the first explicit finite-throat bilinear-kernel choice that is still compatible with the project ontology:

- a wall field `eta(s)` expanded on `{u_0,u_1}`,
- a brane-like internal doublet `U(s)` expanded on the same N/N basis,
- one BdG support half-wave `phi(s) = phi f_0(s)`,
- one mixed `A_w/F_(mu w)/J^w` half-wave `W(s) = W f_0(s)`.

Take the simplest local isotropic couplings,

`L_(eta U) = g_U int_0^L eta(s) U(s) ds`,

`L_(eta phi) = g_B int_0^L eta(s) phi(s) ds`,

`L_(eta W) = g_W int_0^L eta(s) W(s) ds`,

`L_(U W) = - g_R int_0^L U(s) W(s) ds`.

Expanding in modes gives the exact reduced couplings

`L_(eta U) = g_U q.u`,

`L_(eta phi) = g_B (v.q) phi`,

`L_(eta W) = g_W (v.q) W`,

`L_(U W) = - g_R (v.u) W`.

So in the first explicit finite-throat kernel model,

`lambda_U = g_U`,

`lambda_B = g_B`,

`lambda_W = g_W`,

`lambda_R = g_R`,

and the entire directional structure is carried by the same overlap vector `v`.

This is exactly the pattern that had been assumed abstractly in Stage 12; now it is derived from the explicit finite-throat basis and local isotropic kernels.

---

## 3. Exact Schur-complement reduction with explicit kernel-level couplings

Let the reduced frequency-domain kernels be

`A_phi(omega) = varpi^2 - omega^2`,

`A_U(omega)   = Omega_U^2 - omega^2`,

`A_W(omega)   = Omega_W^2 - omega^2 - Pi_out(omega)`.

With the local couplings above, elimination of the internal block `(u, phi, W)` gives the same exact wall self-energy structure found in Stage 12,

`Sigma_wall(omega) = Xi(omega) I_2 + alpha(omega) v v^T`,

but now with explicit kernel-level meaning:

`Xi(omega) = g_U^2 / A_U(omega)`,

`alpha(omega) = g_B^2 / A_phi(omega)`
`             + ( A_U(omega) g_W + g_R g_U )^2`
`               / [ A_U(omega) Delta_UW(omega) ]`,

`Delta_UW(omega) = A_U(omega) A_W(omega) - g_R^2 sigma`.

So the static selected-branch data are no longer abstract at this level.
They are

`Xi_0 = g_U^2 / Omega_U^2`,

`Delta_0 = Omega_U^2 Omega_W^2 - g_R^2 sigma`,

`alpha_0 = g_B^2 / varpi^2 + ( Omega_U^2 g_W + g_R g_U )^2 / ( Omega_U^2 Delta_0 )`,

`beta_0 = ( Omega_U^2 g_W + g_R g_U )^2 / Delta_0^2`.

This is the first point where the Stage-12/14 selected-branch spectral data have been written directly in terms of explicit finite-throat mode integrals and local bilinear kernel strengths.

---

## 4. Natural D/N source branch and the exact selected source map

Now attach the orbital/worldtube STF quadrupole source through the same D/N mouth load,

`L_src = g_Q Q_STF int_0^L eta(s) f_0(s) ds`.

Expanding `eta(s) = q_0 u_0(s) + q_1 u_1(s)` gives

`L_src = g_Q Q_STF (v.q)`.

So the external source vector in wall-coordinate space is simply

`J_src = g_Q Q_STF v`.

Projecting onto the selected lower eigenvector `e_-` therefore gives the selected source amplitude

`J_- = g_Q Q_STF (v.e_-) = g_Q Q_STF sqrt(s_-)`,

where

`s_- = (v.e_-)^2`.

At zero loading the lower mode is just the flat branch, so

`J_-(0) = g_Q Q_STF kappa_0`.

This defines the natural selected-branch source map,

`mhat_- = J_-(alpha_0) / J_-(0) = (v.e_-) / kappa_0`,

and hence the exact squared source-map factor is

`mhat_-^2 = s_- / kappa_0^2`.

So on the natural D/N source branch the abstract factor `mhat_-` is gone.
It is completely fixed by the same selected overlap `s_-` that already controls the wall loading.

---

## 5. Exact bound on the source-map factor

Because the selected overlap grows monotonically from the flat branch to the max-coupling branch,

`kappa_0^2 <= s_- < sigma`,

and therefore

`1 <= mhat_-^2 < sigma / kappa_0^2 = 11/9`.

So the natural source-map factor stays in the exact window

`1 <= mhat_-^2 < 11/9`,

or

`1 <= mhat_- < sqrt(11/9)`.

This is a useful structural simplification.
The selected-branch normalization problem cannot be hidden inside a huge undetermined source-map amplification on the natural D/N source branch.
The source factor is real, positive, monotone, and modest.

---

## 6. Elimination of the abstract source-map factor from the quadrupole target

Stage 13 wrote the selected-branch normalization quantity as

`mhat_-^2 P_{0,-}`

with

`P_{0,-} = beta_0 s_- / lambda_-`.

Substituting the source-map result from Section 4 gives the exact product

`mhat_-^2 P_{0,-} = beta_0 s_-^2 / ( kappa_0^2 lambda_- )`.

So the invariant 2.5PN target becomes

`beta_0 s_-^2 / ( kappa_0^2 lambda_- ) = 54 G c_s^5 / (5 a^5 c^5)`.

This is sharper than the Stage-13 formulation because there is no longer an independent `mhat_-` datum on the natural D/N source branch.
Everything is now carried by:

- the explicit mixed-sector transfer factor `beta_0`,
- the selected overlap `s_-`,
- the selected conservative eigenvalue `lambda_-`,
- and the known flat-branch overlap `kappa_0`.

---

## 7. Best current summary after Stage 15

The selected-branch theorem gap has narrowed again.

The first explicit finite-throat mode-integral model now does three things at once:

1. it derives the Stage-12 loading structure from explicit local isotropic kernels,
2. it writes `Xi_0`, `alpha_0`, and `beta_0` directly in terms of those kernel strengths,
3. and it removes the abstract selected-branch source-map factor by showing that on the natural D/N source branch
   `mhat_-^2 = s_- / kappa_0^2`.

So the remaining normalization problem is no longer

- “selected stiffness plus an unknown source map,”

but rather

- “selected stiffness plus selected overlap plus the explicit mixed-sector transfer factor.”

The next honest step is therefore to write the full selected-branch normalization equation directly in microscopic coupling language and see what exact stability and reachability constraints it imposes.

=== moving_throat_pde_stage016_microscopic_normalization_equation.md ===

# Moving-Throat PDE — Stage 16: Microscopic Selected-Branch Normalization Equation, Exact Stability Gate, and the First Coupling-Level Onset Criterion

## Purpose

Stage 15 removed the abstract selected-branch source-map factor on the natural D/N source branch and rewrote the invariant quadrupole target as

`beta_0 s_-^2 / (kappa_0^2 lambda_-) = 54 G c_s^5 / (5 a^5 c^5)`.

That was a real narrowing, but it still left the theorem gate partly hidden inside the shorthand spectral data `(beta_0, alpha_0, Xi_0)`.

The next step is therefore to write the selected-branch normalization problem directly in the microscopic couplings of the first explicit finite-throat kernel model and to separate three logically different questions:

1. **existence of a stable selected branch**,
2. **entry into the normalization window**,
3. **exact hit of the universal target at the physical loading**.

This stage does that.

The main result is that the selected-branch 2.5PN test is now one exact microscopic equation together with two exact necessary inequalities.
The stable natural branch must satisfy

`Delta_0 > 0`,

`A > 0`,

`alpha_0 < alpha_crit`,

and its normalization must obey

`N_-(g_B,g_U,g_W,g_R,...) = 54 G c_s^5 / (5 a^5 c^5)`

with an explicit closed formula for `N_-`.

There is also a useful exact necessary onset condition:
if the branch starts above the universal target at zero directional loading, then the natural positive-loading branch can never come back down to hit it, because the selected normalization product is strictly monotone increasing.
So the first coupling-level onset test is already exact.

---

## 1. Microscopic coupling abbreviations

For the first explicit finite-throat kernel model, define

`A = K_0 - g_U^2 / Omega_U^2`,

`Delta_0 = Omega_U^2 Omega_W^2 - g_R^2 sigma`,

`Chi = Omega_U^2 g_W + g_R g_U`.

Then the Stage-15 data are

`Xi_0 = g_U^2 / Omega_U^2`,

`beta_0 = Chi^2 / Delta_0^2`,

`alpha_0 = g_B^2 / varpi^2 + Chi^2 / ( Omega_U^2 Delta_0 )`.

So the physical branch is controlled by:

- one isotropic wall shift `Xi_0`,
- one directional loading strength `alpha_0`,
- one outgoing transfer factor `beta_0`,
- and the finite-throat overlap constants `kappa_0`, `kappa_1`, `sigma` already fixed exactly.

---

## 2. Exact microscopic selected-branch normalization equation

On the natural D/N source branch,

`mhat_-^2 = s_- / kappa_0^2`,

`P_{0,-} = beta_0 s_- / lambda_-`,

so the full invariant selected-branch product is

`N_-(alpha_0) := mhat_-^2 P_{0,-}`
`              = beta_0 s_-(alpha_0)^2 / ( kappa_0^2 lambda_-(alpha_0) )`.

The exact 2.5PN target is therefore

`N_-(alpha_0) = N_Q^(target)`

with

`N_Q^(target) = 54 G c_s^5 / (5 a^5 c^5)`.

Substituting the microscopic couplings gives

`Chi^2 / Delta_0^2 * s_-(alpha_0)^2 / ( kappa_0^2 lambda_-(alpha_0) )`
`= 54 G c_s^5 / (5 a^5 c^5)`

where

`alpha_0 = g_B^2 / varpi^2 + Chi^2 / ( Omega_U^2 Delta_0 )`.

So the selected-branch normalization theorem is now a single explicit equation in the first finite-throat kernel-model couplings.

---

## 3. Exact stability gate in coupling language

Before the normalization equation can even be asked, the selected lower wall mode has to exist on a stable branch.
That imposes three exact requirements.

### 3.1 Internal passive/conservative positivity

The conservative internal block must satisfy

`Delta_0 = Omega_U^2 Omega_W^2 - g_R^2 sigma > 0`.

If `Delta_0 <= 0`, the reduced `U/W` block is already singular or overmixed before the wall problem is even formed.

### 3.2 Positive flat-branch stiffness

The flat branch must start stable,

`A = K_0 - g_U^2 / Omega_U^2 > 0`.

So the isotropic internal zero-mode shift cannot already overwhelm the bare wall stiffness.

### 3.3 Selected-branch softening bound

The directional loading must remain below the exact refined threshold

`alpha_0 < alpha_crit`,

with

`alpha_crit = A(A + Delta K_ax)`
`             / [ (A + Delta K_ax) kappa_0^2 + A kappa_1^2 ]`.

Using the exact finite-throat constants,

`kappa_0^2 = 8 / pi^2`,

`kappa_1^2 = 16 / (9 pi^2)`,

this becomes

`alpha_crit = 9 pi^2 A(A + Delta K_ax) / [ 8(11 A + 9 Delta K_ax) ]`.

So the exact coupling-level stability gate is

`g_B^2 / varpi^2 + Chi^2 / ( Omega_U^2 Delta_0 )`
`< 9 pi^2 A(A + Delta K_ax) / [ 8(11 A + 9 Delta K_ax) ]`.

This is the first fully explicit stability window for the selected moving-throat quadrupole branch in the current reduced program.

---

## 4. Exact monotonicity and the branch-onset criterion

From Stage 14, the selected normalization product is strictly monotone increasing on the stable branch.
With the Stage-15 source-map reduction, that exact statement becomes

`d N_- / d alpha_0`
`= beta_0 [ 2 s_- (d s_- / d alpha_0) lambda_- + s_-^3 ]`
`  / ( kappa_0^2 lambda_-^2 )`
`> 0`

for every stable branch point.

So the zero-loading value is now an exact **necessary onset condition** for the universal target.
At `alpha_0 = 0`,

`s_-(0) = kappa_0^2`,

`lambda_-(0) = A`,

and therefore

`N_-(0) = beta_0 kappa_0^2 / A = Chi^2 kappa_0^2 / ( A Delta_0^2 )`.

Because `N_-` only increases with positive loading, a necessary condition for the physical natural branch to hit the universal target is

`N_-(0) <= N_Q^(target)`.

Equivalently,

`Chi^2 <= N_Q^(target) A Delta_0^2 / kappa_0^2`.

This may be read as an onset stiffness bound,

`K_0 >= g_U^2 / Omega_U^2 + kappa_0^2 Chi^2 / ( N_Q^(target) Delta_0^2 )`.

If this inequality fails, then the natural positive-loading branch starts **above** the universal target and can never come back down to hit it.
That does not yet solve the normalization equation, but it is an exact necessary branch-admissibility test.

---

## 5. Weak-loading expansion of the microscopic normalization product

For a stable branch with small directional loading,

`alpha_0 << alpha_crit`,

the exact selected normalization product has the expansion

`N_-(alpha_0)`
`= beta_0 kappa_0^2 / A`
`  + alpha_0 beta_0 kappa_0^2 ( 4 A kappa_1^2 + Delta K_ax kappa_0^2 ) / ( A^2 Delta K_ax )`
`  + O(alpha_0^2)`.

Using the exact finite-throat constants gives

`N_-(alpha_0)`
`= beta_0 * 8 / (pi^2 A)`
`  + alpha_0 * 64 beta_0 (8A + 9 Delta K_ax)`
`    / ( 9 pi^4 A^2 Delta K_ax )`
`  + O(alpha_0^2)`.

Now substitute the microscopic loading strength

`alpha_0 = g_B^2 / varpi^2 + Chi^2 / ( Omega_U^2 Delta_0 )`.

Then the first explicit weak-loading approximation to the physical selected-branch product is

`N_-^(phys)`
`= Chi^2 * 8 / ( pi^2 A Delta_0^2 )`
`  + [ g_B^2 / varpi^2 + Chi^2 / ( Omega_U^2 Delta_0 ) ]`
`    * 64 Chi^2 (8A + 9 Delta K_ax)`
`      / ( 9 pi^4 A^2 Delta K_ax Delta_0^2 )`
`  + O(alpha_0^2)`.

This is not yet the final theorem, but it is the first concrete approximation that lets us diagnose which microscopic lane is doing the main work:

- `g_U` lowers the wall through `A`,
- `g_R` and `g_W` feed the transfer factor through `Chi` and `Delta_0`,
- `g_B` pushes the branch along the monotone loading direction,
- and `Delta K_ax` controls how fast the selected source map can grow before softening.

---

## 6. Best current theorem gate after Stage 16

The selected-branch normalization bottleneck is now sharply microscopic.
The first explicit finite-throat kernel model has reduced the theorem gap to:

1. compute the physical coupling set
   `(g_B, g_U, g_W, g_R, varpi, Omega_U, Omega_W, K_0, Delta K_ax)`
   from the completed moving-throat PDE,
2. verify the exact stability window
   `Delta_0 > 0`, `A > 0`, `alpha_0 < alpha_crit`,
3. check the exact necessary onset inequality
   `N_-(0) <= N_Q^(target)`,
4. and then test the full microscopic normalization equation
   `N_-(alpha_0) = N_Q^(target)`.

So the open problem is no longer hidden in a free source map or a generic tuning story.
It is an explicit coupling-level spectral-placement problem on the selected stable quadrupole branch.

=== moving_throat_pde_stage017_softening_depth_normal_form.md ===

# Moving-Throat PDE — Stage 17: Softening-Depth Normal Form, Exact Secular Reduction, and Elimination of the Selected-Mode Eigenvector Algebra

## Purpose

Stage 16 reduced the selected-branch 2.5PN normalization problem to one explicit microscopic equation,

`N_-(alpha_0) = N_Q^(target)`,

with `N_-` written in terms of the selected wall eigenvalue `lambda_-`, the selected overlap `s_-`, and the total directional loading `alpha_0`.

That was already sharp, but the formulation still carried unnecessary spectral baggage:
- the selected eigenvector `e_-`,
- the selected overlap `s_- = (v.e_-)^2`,
- and the branch stiffness `lambda_-`.

The next natural step is therefore to trade the selected eigenvalue for the more physical quantity

`x := A - lambda_-`,

namely the **softening depth** of the selected lower wall mode below the flat-branch stiffness `A`.

This removes the eigenvector algebra entirely.
In the first rank-1 selected-branch model, the full normalization problem becomes an explicit scalar problem in `x`.
The main outputs are:

1. an exact secular law `alpha_0(x)`,
2. an exact selected overlap `s_-(x)`,
3. an exact normalization product `N_-(x)`,
4. a manifestly monotone loading map `d alpha_0 / dx > 0`,
5. and an exact required support-loading formula `g_B^2/varpi^2` once `x` is known.

So Stage 17 replaces the spectral unknowns `(lambda_-, e_-)` by a single scalar deformation coordinate `x`.

---

## 1. Exact softening-depth variable

Recall the selected lower wall mode of the Stage-12/16 rank-1 reduced operator,

`D_wall = diag(A, A + DeltaK_ax) - alpha_0 v v^T`,

with

`v = (kappa_0, kappa_1)^T`.

On the stable branch,

`0 < lambda_- < A`,

so define the exact softening depth

`x := A - lambda_-`.

Then

`lambda_- = A - x`,

and stability is simply

`0 <= x < A`.

The selected branch starts at `x=0` when the directional loading vanishes and approaches softening as `x -> A^-`.

---

## 2. Exact secular equation in softening-depth form

The selected branch is the nontrivial rank-1-loaded solution of

`1 = alpha_0 [ kappa_0^2 / (A - lambda_-) + kappa_1^2 / (A + DeltaK_ax - lambda_-) ]`.

Substituting `lambda_- = A - x` gives the exact softening-depth secular equation

`1 = alpha_0 [ kappa_0^2 / x + kappa_1^2 / (x + DeltaK_ax) ]`.

Solving for the total directional loading gives

`alpha_0(x) = x (x + DeltaK_ax) / [ kappa_0^2 (x + DeltaK_ax) + kappa_1^2 x ]`.

This is the first key simplification:

the selected-branch loading is now an explicit rational function of the softening depth.

---

## 3. Exact selected overlap in softening-depth form

The selected source/wall overlap may be written through the exact secular derivative identity,

`s_- = - d lambda_- / d alpha_0`.

Because `lambda_- = A - x`, this becomes

`s_- = d x / d alpha_0`.

Using the secular law above, the exact overlap is

`s_-(x) = [ kappa_0^2 (x + DeltaK_ax) + kappa_1^2 x ]^2`
`         / [ kappa_0^2 (x + DeltaK_ax)^2 + kappa_1^2 x^2 ]`.

So the selected overlap is also a rational function of the same scalar variable `x`.
No explicit eigenvector is needed anymore.

---

## 4. Exact selected-branch normalization product in softening-depth form

Stage 15/16 gave the invariant selected-branch normalization quantity as

`N_-(alpha_0) = beta_0 s_-^2 / ( kappa_0^2 lambda_- )`.

Now substitute

`lambda_- = A - x`,

and the exact overlap formula from Section 3.
This gives

`N_-(x) = beta_0 [ kappa_0^2 (x + DeltaK_ax) + kappa_1^2 x ]^4`
`         / { kappa_0^2 (A - x) [ kappa_0^2 (x + DeltaK_ax)^2 + kappa_1^2 x^2 ]^2 }`.

So the full selected-branch quadrupole normalization theorem is no longer an eigenvector problem.
It is one scalar equation in the softening depth:

`N_-(x) = N_Q^(target)`.

---

## 5. Exact monotonicity of the loading map

Differentiate the secular law.
One obtains

`d alpha_0 / dx`
`= [ kappa_0^2 (x + DeltaK_ax)^2 + kappa_1^2 x^2 ]`
`  / [ kappa_0^2 (x + DeltaK_ax) + kappa_1^2 x ]^2`
`> 0`.

So the total directional loading grows strictly monotonically with softening depth.
This means:

- every stable selected-branch loading corresponds to exactly one softening depth,
- and vice versa.

This is useful because it turns the Stage-14 monotonicity statement into a direct scalar branch parameterization.

---

## 6. Exact required support loading once the softening depth is known

Stage 16 already split the total directional loading into

`alpha_0 = g_B^2 / varpi^2 + alpha_mix`,

with

`alpha_mix = Chi^2 / ( Omega_U^2 Delta_0 )`,

`Chi = Omega_U^2 g_W + g_R g_U`,

`Delta_0 = Omega_U^2 Omega_W^2 - g_R^2 sigma`.

So once the required softening depth `x` is known, the exact total loading is `alpha_0(x)`, and the exact required support contribution is

`g_B,req^2 / varpi^2 = alpha_0(x) - alpha_mix`

or explicitly

`g_B,req^2 / varpi^2`
`= x (x + DeltaK_ax) / [ kappa_0^2 (x + DeltaK_ax) + kappa_1^2 x ]`
`  - Chi^2 / ( Omega_U^2 Delta_0 )`.

So the support/BdG loading needed to hit the universal 2.5PN target is now an explicit function of one scalar deformation coordinate.

---

## 7. Best current theorem gate after Stage 17

The selected-branch theorem bottleneck is now even narrower than at Stage 16.
Instead of solving for a loaded eigenvector and its overlap, the reduced program only needs to determine the physical softening depth `x` of the selected lower wall mode.
Once that is known,

- the total directional loading is fixed by `alpha_0(x)`,
- the selected source map is fixed by `s_-(x)`,
- the normalization product is fixed by `N_-(x)`,
- and the required support loading follows directly.

So the remaining gap is no longer hidden in eigenvector algebra.
It is a scalar branch-placement problem.

=== moving_throat_pde_stage018_dimensionless_normalization_locus.md ===

# Moving-Throat PDE — Stage 18: Dimensionless D/N Shape Function, Unique Normalization Locus, and Exact Required Support Coupling

## Purpose

Stage 17 rewrote the selected-branch normalization problem as one scalar equation in the softening depth `x`.
That was already a real simplification, but the formulas still carried the dimensional parameters `(A, DeltaK_ax, beta_0)` in a way that partly obscured the geometry of the branch.

The next step is therefore to pass to the exact D/N finite-throat constants and reduce the entire problem to a **dimensionless shape function**.

This is worthwhile because it gives the strongest selected-branch result so far:

1. the whole normalization problem collapses to one equation
   `F(xi,delta) = R_target`,
2. the shape function `F` is strictly monotone on the stable branch,
3. it maps the stable branch exactly from `1` to `+infinity`,
4. so the normalization locus is unique whenever the onset condition is satisfied,
5. and the exact required total loading and required support loading then follow immediately.

So Stage 18 turns the microscopic normalization problem into a one-dimensional universal branch-geometry problem plus a final support-coupling feasibility check.

---

## 1. Exact D/N constants and dimensionless variables

Insert the exact finite-throat D/N overlap constants,

`kappa_0^2 = 8 / pi^2`,

`kappa_1^2 = 16 / (9 pi^2)`,

so the exact overlap ratio is

`eta := kappa_1^2 / kappa_0^2 = 2/9`.

Now define dimensionless selected-branch variables

`xi := x / A`,

`delta := DeltaK_ax / A`.

Then stable selected branches satisfy

`0 <= xi < 1`,

while `delta > 0` is the axial anisotropy of the bare wall block in units of the flat stiffness.

---

## 2. Exact dimensionless shape function

Using the Stage-17 softening-depth normal form, the exact normalization product may be written as

`N_-(x) = N_-(0) F(xi,delta)`,

with

`N_-(0) = beta_0 kappa_0^2 / A`.

After inserting the D/N constants, the exact dimensionless shape function is

`F(xi,delta)`
`= (9 delta + 11 xi)^4`
`  / [ 81 (1 - xi) (9 delta^2 + 18 delta xi + 11 xi^2)^2 ]`.

So the full selected-branch target becomes

`F(xi,delta) = R_target`,

where

`R_target := N_Q^(target) / N_-(0)`
`          = N_Q^(target) A / ( beta_0 kappa_0^2 )`.

This is the cleanest reduced theorem form yet.
All dependence on the microscopic transfer factor and the overall wall scale has collapsed into the single dimensionless target ratio `R_target`.

---

## 3. Exact monotonicity, onset, and softening limits

Differentiate the exact D/N shape function.
One obtains

`dF/dxi`
`= (9 delta + 11 xi)^3`
`  [ 81 delta^3 + 72 delta^2 + 189 delta^2 xi + 297 delta xi^2 + 121 xi^3 ]`
`  / [ 81 (1 - xi)^2 (9 delta^2 + 18 delta xi + 11 xi^2)^3 ]`
`> 0`.

So `F` is strictly monotone increasing on the entire stable branch `0 <= xi < 1`.

Its exact endpoint values are

`F(0,delta) = 1`,

`lim_{xi -> 1^-} F(xi,delta) = +infinity`.

More precisely,

`F(xi,delta) ~ C(delta) / (1 - xi)`

as `xi -> 1^-`, with

`C(delta) = (9 delta + 11)^4 / [ 81 (9 delta^2 + 18 delta + 11)^2 ]`.

This is a theorem-level improvement over Stage 16:

- if `R_target < 1`, the target is impossible on the stable branch,
- if `R_target = 1`, the only hit is the unloaded onset point `xi = 0`,
- if `R_target > 1`, there is one and only one stable selected-branch softening depth `xi_req` that hits the target.

So the Stage-16 onset inequality is now upgraded to a full uniqueness theorem for the selected normalization locus.

---

## 4. Exact required total loading

The Stage-17 total loading law becomes, in D/N dimensionless form,

`alpha_req(xi,delta)`
`= 9 pi^2 A xi (xi + delta) / [ 8 (9 delta + 11 xi) ]`.

As `xi -> 1^-`, this tends to the exact stable-branch softening threshold

`alpha_crit = 9 pi^2 A (1 + delta) / [ 8 (11 + 9 delta) ]`,

which is the same refined threshold already found in Stage 16.

So once the unique `xi_req` solving `F(xi,delta)=R_target` is known, the unique required total directional loading follows immediately from the formula above.

---

## 5. Exact required support coupling

The physical total loading still splits as

`alpha_0 = g_B^2 / varpi^2 + alpha_mix`,

with

`alpha_mix = Chi^2 / ( Omega_U^2 Delta_0 )`.

Therefore the exact support/BdG loading needed to realize the selected normalization locus is

`g_B,req^2 / varpi^2 = alpha_req(xi_req,delta) - alpha_mix`

or explicitly

`g_B,req^2 / varpi^2`
`= 9 pi^2 A xi_req (xi_req + delta) / [ 8 (9 delta + 11 xi_req) ]`
`  - Chi^2 / ( Omega_U^2 Delta_0 )`.

This gives one more exact feasibility gate:

`g_B,req^2 >= 0`.

So even after the unique normalization locus is found, the completed moving-throat PDE still has to place the physical support coupling on the nonnegative side of this equation.

---

## 6. Useful asymptotics

### 6.1 Near-onset regime

For `xi << 1`, the exact D/N shape function expands as

`F(xi,delta)`
`= 1 + ( 1 + 8/(9 delta) ) xi`
`  + ( 1 + 8/(9 delta) - 28/(27 delta^2) ) xi^2 + O(xi^3)`.

So if the target is only slightly above onset,

`R_target = 1 + eps_R`,

then the unique selected branch point is approximately

`xi_req ~= eps_R / ( 1 + 8/(9 delta) )`.

The required total loading then begins as

`alpha_req ~= pi^2 A xi_req / 8`.

### 6.2 Near-softening regime

For very large target ratio `R_target`, the unique solution lies close to softening and obeys

`1 - xi_req ~= C(delta) / R_target`,

with

`C(delta) = (9 delta + 11)^4 / [ 81 (9 delta^2 + 18 delta + 11)^2 ]`.

So large normalization demand pushes the selected branch into the thin neighborhood just below softening.

---

## 7. Best current theorem gate after Stage 18

The selected-branch quadrupole normalization bottleneck has now been reduced to the smallest exact reduced form reached so far.

The completed moving-throat PDE has to determine the microscopic quantities entering

`R_target = N_Q^(target) A / ( beta_0 kappa_0^2 )`,

and the bare anisotropy ratio

`delta = DeltaK_ax / A`.

Once those are known,

1. the unique stable normalization locus `xi_req` is determined by `F(xi,delta)=R_target`,
2. the unique required total loading is `alpha_req(xi_req,delta)`,
3. and the unique required support coupling is `g_B,req^2 / varpi^2 = alpha_req - alpha_mix`.

So the open problem is no longer a vague multi-parameter search.
It is now a one-dimensional universal D/N branch-shape problem plus one exact support-coupling feasibility check.

=== moving_throat_pde_stage019_support_feasibility_frontier.md ===

# Moving-Throat PDE — Stage 19: Dimensionless Support-Feasibility Frontier for the Selected Quadrupole Branch

## Purpose

Stage 18 reduced the selected-branch normalization problem to the unique D/N locus

`F(xi,delta) = R_target`,

and gave the exact total loading required once the unique solution `xi_req` is found.
But there was still one extra feasibility step left outside the locus itself:

the support/BdG channel has to supply a **nonnegative** additional loading,

`g_B,req^2 >= 0`.

The natural next move is therefore to isolate the exact dimensionless function that measures how much of the total directional loading can be carried by the selected branch itself.
That is what this stage does.

The result is an exact second branch function,

`G(xi,delta) = 9 xi (xi + delta) / (9 delta + 11 xi)`,

such that

`g_B,req^2 / varpi^2 = (pi^2 A / 8) [ G(xi_req,delta) - M_mix ]`,

with

`M_mix = 8 alpha_mix / (pi^2 A)`
`      = 8 Chi^2 / (pi^2 A Omega_U^2 Delta_0)`.

So the selected quadrupole branch is support-feasible iff

`M_mix <= G(xi_req,delta)`.

This turns the final support check into an exact geometric inequality in the same dimensionless branch parameter `xi` that already controlled the normalization locus.

---

## 1. Exact support-feasibility function

From Stage 18 the required total loading is

`alpha_req(xi,delta) = 9 pi^2 A xi (xi + delta) / [ 8 (9 delta + 11 xi) ]`.

Split the total loading into the mixed-sector baseline and the support contribution,

`alpha_req = alpha_mix + g_B,req^2 / varpi^2`,

with

`alpha_mix = Chi^2 / ( Omega_U^2 Delta_0 )`.

Now define the dimensionless mixed baseline

`M_mix := 8 alpha_mix / (pi^2 A)`
`      = 8 Chi^2 / (pi^2 A Omega_U^2 Delta_0 )`.

Then the exact branch-dependent support-feasibility function is

`G(xi,delta) := 8 alpha_req / (pi^2 A)`
`           = 9 xi (xi + delta) / (9 delta + 11 xi)`.

The required support loading becomes

`g_B,req^2 / varpi^2 = (pi^2 A / 8) [ G(xi_req,delta) - M_mix ]`.

So the support/BdG sector is feasible exactly when

`G(xi_req,delta) >= M_mix`.

---

## 2. Exact monotonicity and endpoint values

Differentiate `G`.
One finds

`dG/dxi = 9 [ 9 delta^2 + 18 delta xi + 11 xi^2 ] / (9 delta + 11 xi)^2 > 0`.

So `G` is strictly monotone increasing on the stable branch.
Its exact endpoint values are

`G(0,delta) = 0`,

`G_max(delta) := lim_{xi -> 1^-} G(xi,delta)`
`             = 9 (1 + delta) / (9 delta + 11)`.

This is useful because it turns the support-feasibility condition into a sharp branch window.
For fixed `delta`, the selected branch can support at most

`M_mix < G_max(delta)`.

That is of course equivalent to the refined stability bound `alpha_mix < alpha_crit`, but the present stage makes the same statement directly in the dimensionless selected-branch geometry.

---

## 3. The exact admissible region in `(R_target, M_mix)` space

Stages 18 and 19 together now give two exact branch functions driven by the same parameter `xi`:

`R_target = F(xi,delta)`,

`M_crit = G(xi,delta)`.

For fixed `delta`, the stable selected quadrupole branch therefore traces an exact parametric admissibility frontier in the `(R_target, M_mix)` plane:

`xi in [0,1)  ->  ( F(xi,delta), G(xi,delta) )`.

Because both `F` and `G` are strictly monotone increasing,

- the normalization target picks out a unique `xi_req`,
- and the support feasibility test is then simply whether the actual baseline `M_mix` lies below the critical value `G(xi_req,delta)`.

So the combined theorem gate is now:

1. `R_target >= 1` to enter the unique normalization locus,
2. `M_mix <= G(xi_req,delta)` to keep the required support loading nonnegative.

That is the first exact combined reachability-plus-feasibility statement for the selected moving-throat quadrupole branch.

---

## 4. Near-onset asymptotics

For `xi << 1`, the support-feasibility function expands as

`G(xi,delta) = xi - 2 xi^2 / (9 delta) + O(xi^3)`.

Combined with the Stage-18 onset relation

`xi_req ~= (R_target - 1) / (1 + 8/(9 delta))`,

this gives the first near-onset support-feasibility estimate,

`M_crit ~= (R_target - 1) / (1 + 8/(9 delta))`

up to the first nonlinear correction.

So just above onset, the admissible mixed baseline grows linearly with the excess normalization demand.

---

## 5. Best current theorem gate after Stage 19

The selected moving-throat quadrupole problem has now split cleanly into two exact scalar branch functions:

1. the **normalization function** `F(xi,delta)`,
2. the **support-feasibility function** `G(xi,delta)`.

For any fixed `delta`, the completed moving-throat PDE must determine a physical point `(R_target, M_mix)` such that

- `R_target` lands on the unique stable normalization locus,
- and `M_mix` lies below the exact support-feasibility frontier.

So the open theorem gap is now even sharper:

not “find the right eigenvector and source map,”
not even “find the right branch stiffness,”

but simply:

> does the completed moving-throat PDE place the physical defect on the exact admissible region of the universal `(F,G)` branch geometry?

=== moving_throat_pde_stage020_continuum_kernel_extraction.md ===

# Moving-Throat PDE — Stage 20: Exact Continuum-Kernel Extraction of `A`, `DeltaK_ax`, `beta_0`, and `M_mix`

## Purpose

Stage 19 reduced the selected-branch theorem gate to the exact admissibility geometry
driven by

- the normalization ratio `R_target`,
- the mixed baseline `M_mix`,
- and the bare anisotropy ratio `delta`.

But those quantities were still being treated as microscopic inputs.

The right next move is therefore to derive them from the first explicit **linearized finite-throat continuum operator** that is still consistent with the reduced hierarchy already used in Stages 10–19.

That is what this stage does.

The main result is that the first continuum operator already produces all of the key reduced quantities in closed form. After exact mode projection and mass normalization, one finds

`A = [ K_U (K_eta + 6 T_Omega) - c_(eta U)^2 ] / (mu_eta K_U),`

`DeltaK_ax = pi^2 T_w / (mu_eta L^2),`

`beta_0 = (mu_W / mu_eta) ( K_U c_(eta W) + c_(UW) c_(eta U) )^2`
`         / [ ( K_U K_W^(eff) - c_(UW)^2 sigma )^2 ],`

`alpha_mix = ( K_U c_(eta W) + c_(UW) c_(eta U) )^2`
`            / [ mu_eta K_U ( K_U K_W^(eff) - c_(UW)^2 sigma ) ],`

`M_mix = 8 alpha_mix / (pi^2 A),`

with

`K_W^(eff) := K_W + pi^2 T_W / (4 L^2),`
`sigma = 88 / (9 pi^2).`

So the Stage-17/19 branch variables are no longer abstract. They are exact low-mode functionals of one explicit continuum kernel.

---

## 1. Minimal linearized continuum operator

Work on the finite throat interval

`s in [0,L]`.

Keep the same wall/support ontology already used in the reduced stages:

- a quadrupole wall field `eta(s,t)` with N/N boundaries,
- a brane-like internal field `U(s,t)` on the same interval,
- a support/BdG field `phi(s,t)` on the D/N half-wave branch,
- a mixed `W(s,t)` field representing the `A_w/F_(mu w)/J^w` lane on the same D/N branch.

The minimal quadratic Lagrangian density is taken to be

`L_eta = (mu_eta/2) dot(eta)^2 - (T_w/2) (eta')^2 - (K_eta + 6 T_Omega) eta^2 / 2,`

`L_U   = (mu_U/2) dot(U)^2 - K_U U^2 / 2,`

`L_phi = (mu_phi/2) dot(phi)^2 - (T_phi/2) (phi')^2 - K_phi phi^2 / 2,`

`L_W   = (mu_W/2) dot(W)^2 - (T_W/2) (W')^2 - K_W W^2 / 2,`

with local bilinear couplings

`L_(eta U)   = - c_(eta U) int_0^L eta U ds,`

`L_(eta phi) = - c_(eta phi) int_0^L eta phi ds,`

`L_(eta W)   = - c_(eta W) int_0^L eta W ds,`

`L_(UW)      = + c_(UW) int_0^L U W ds.`

Two comments matter.

First, the wall operator is exactly the same quadrupole wall operator already carried from the earlier stages.

Second, the internal `U` field is kept in the same **brane-like flat-doublet limit** already used implicitly in Stages 12–16: at this minimal order there is no axial-gradient penalty in `U`, so the first two N/N channels stay degenerate.
That is the smallest continuum choice that reproduces the earlier reduced hierarchy.

---

## 2. Exact mode basis and overlaps

Use the exact first two N/N modes

`u_0(s) = 1 / sqrt(L),`

`u_1(s) = sqrt(2/L) cos(pi s / L),`

and the exact lowest D/N half-wave

`f_0(s) = sqrt(2/L) sin(pi s / (2L)).`

The exact D/N overlap vector is

`v_i = int_0^L u_i(s) f_0(s) ds,`

so

`kappa_0 = 2 sqrt(2) / pi,`

`kappa_1 = - 4 / (3 pi),`

`v = (kappa_0, kappa_1)^T,`

`sigma = v.v = 88 / (9 pi^2).`

This is the same overlap data already driving the Stage-10/15 profile and source maps, but now it is being used directly as the projection datum of a continuum operator.

---

## 3. Mass-normalized modal coordinates

Expand

`eta(s,t) = q_0(t) u_0(s) + q_1(t) u_1(s),`

`U(s,t)   = u_0^(int)(t) u_0(s) + u_1^(int)(t) u_1(s),`

`phi(s,t) = phi_0(t) f_0(s),`

`W(s,t)   = W_0(t) f_0(s).`

Now pass to mass-normalized modal coordinates

`Q_i   = sqrt(mu_eta) q_i,`

`U_i   = sqrt(mu_U) u_i^(int),`

`Phi   = sqrt(mu_phi) phi_0,`

`Wbar  = sqrt(mu_W) W_0.`

Then the reduced bare kernels are

`K_0 = (K_eta + 6 T_Omega) / mu_eta,`

`DeltaK_ax = pi^2 T_w / (mu_eta L^2),`

`varpi^2 = ( K_phi + pi^2 T_phi / (4 L^2) ) / mu_phi,`

`Omega_U^2 = K_U / mu_U,`

`Omega_W^2 = ( K_W + pi^2 T_W / (4 L^2) ) / mu_W`
`          = K_W^(eff) / mu_W.`

The mass-normalized couplings are

`g_U = c_(eta U) / sqrt(mu_eta mu_U),`

`g_B = c_(eta phi) / sqrt(mu_eta mu_phi),`

`g_W = c_(eta W) / sqrt(mu_eta mu_W),`

`g_R = c_(UW) / sqrt(mu_U mu_W).`

So the projected couplings become

`L_(eta U)   = - g_U Q.U,`

`L_(eta phi) = - g_B (v.Q) Phi,`

`L_(eta W)   = - g_W (v.Q) Wbar,`

`L_(UW)      = + g_R (v.U) Wbar.`

That is exactly the reduced pattern used abstractly in Stage 12. The difference is that every symbol in it has now been tied back to the continuum operator.

---

## 4. Exact Schur complement of the internal block

In frequency space the internal kernels are

`A_U(omega)   = Omega_U^2 - omega^2,`

`A_phi(omega) = varpi^2 - omega^2,`

`A_W(omega)   = Omega_W^2 - omega^2 - Pi_out(omega).`

Eliminating `(U, Phi, Wbar)` gives the exact wall self-energy

`Sigma_wall(omega) = Xi(omega) I_2 + alpha(omega) v v^T,`

with

`Xi(omega) = g_U^2 / A_U(omega),`

`alpha(omega) = g_B^2 / A_phi(omega)`
`             + ( A_U(omega) g_W + g_R g_U )^2`
`               / [ A_U(omega) Delta_UW(omega) ],`

`Delta_UW(omega) = A_U(omega) A_W(omega) - g_R^2 sigma.`

So the Stage-11 rank-1 loading law is not phenomenological.
It is the exact first Schur complement of the continuum operator.

On the conservative static branch `Pi_out(0)=0`, define

`Delta_0 = Omega_U^2 Omega_W^2 - g_R^2 sigma,`

`Chi = Omega_U^2 g_W + g_R g_U.`

Then the key reduced quantities are

`A = K_0 - g_U^2 / Omega_U^2,`

`alpha_mix = Chi^2 / ( Omega_U^2 Delta_0 ),`

`beta_0 = Chi^2 / Delta_0^2,`

`M_mix = 8 alpha_mix / (pi^2 A),`

`delta = DeltaK_ax / A.`

This is exactly the reduced data set that was still abstract at Stage 19.

---

## 5. Closed continuum formulas

Now substitute the mass-normalized projections back into the definitions above.

Introduce

`K_eta^(eff) := K_eta + 6 T_Omega,`

`K_W^(eff)   := K_W + pi^2 T_W / (4 L^2).`

Then

`A = [ K_U K_eta^(eff) - c_(eta U)^2 ] / ( mu_eta K_U ),`

`DeltaK_ax = pi^2 T_w / ( mu_eta L^2 ),`

`Delta_0 = [ K_U K_W^(eff) - c_(UW)^2 sigma ] / ( mu_U mu_W ),`

`Chi = [ K_U c_(eta W) + c_(UW) c_(eta U) ]`
`      / [ mu_U sqrt(mu_eta mu_W) ].`

So the static mixed loading is

`alpha_mix`
`= ( K_U c_(eta W) + c_(UW) c_(eta U) )^2`
`  / [ mu_eta K_U ( K_U K_W^(eff) - c_(UW)^2 sigma ) ],`

the outgoing transfer factor is

`beta_0`
`= (mu_W / mu_eta)`
`  ( K_U c_(eta W) + c_(UW) c_(eta U) )^2`
`  / [ ( K_U K_W^(eff) - c_(UW)^2 sigma )^2 ],`

and the dimensionless mixed baseline is

`M_mix`
`= 8 ( K_U c_(eta W) + c_(UW) c_(eta U) )^2`
`  / [ pi^2 ( K_U K_eta^(eff) - c_(eta U)^2 )`
`          ( K_U K_W^(eff) - c_(UW)^2 sigma ) ].`

Finally, the bare anisotropy ratio is

`delta`
`= pi^2 T_w K_U`
`  / [ L^2 ( K_U K_eta^(eff) - c_(eta U)^2 ) ].`

So the branch data entering Stage 17–19 are now explicit continuum-kernel functions.

---

## 6. Exact continuum stability inequalities

The selected-branch stability gates now become exact inequalities on the continuum operator.

### 6.1 Wall/internal flat-branch stability

`A > 0`
is equivalent to

`K_U K_eta^(eff) > c_(eta U)^2.`

So the wall/brane-like internal coupling cannot over-soften the flat wall branch.

### 6.2 Internal mixed-sector positivity

`Delta_0 > 0`
is equivalent to

`K_U K_W^(eff) > c_(UW)^2 sigma.`

So the internal `U/W` block must stay below its continuum overmixing threshold.

These are the first exact continuum-kernel positivity surfaces that underwrite the Stage-16/19 reduced branch geometry.

---

## 7. Best current theorem gate after Stage 20

The open theorem gap has now narrowed again.

The selected quadrupole branch is no longer parameterized by abstract reduced inputs.
The first explicit finite-throat continuum operator already determines

- `A`,
- `DeltaK_ax`,
- `alpha_mix`,
- `beta_0`,
- `M_mix`,
- and therefore `delta`.

So the remaining task is no longer “invent microscopic inputs.”
It is:

> determine whether the completed moving-throat PDE places the actual defect on the admissible Stage-18/19 branch geometry generated by these continuum-kernel quantities.

=== moving_throat_pde_stage021_dimensionless_continuum_placement.md ===

# Moving-Throat PDE — Stage 21: Dimensionless Continuum Placement Map, Exact Product Relation, and the Three-Lane Factorization of the Selected Quadrupole Branch

## Purpose

Stage 20 derived the reduced branch data directly from one explicit finite-throat continuum kernel.
That was already a meaningful advance, but it still left the placement problem written in a somewhat redundant microscopic language.

The next step is therefore to compress the Stage-20 continuum operator into the smallest useful **dimensionless kernel ledger**.

The main result is that the continuum placement problem collapses to five exact kernel ratios:

`eps_eta = c_(eta U)^2 / ( K_U K_eta^(eff) ),`

`eps_W   = c_(UW)^2 sigma / ( K_U K_W^(eff) ),`

`rho     = c_(UW) c_(eta U) / ( K_U c_(eta W) ),`

`Z_W     = c_(eta W)^2 / ( K_eta^(eff) K_W^(eff) ),`

`delta_0 = pi^2 T_w / ( L^2 K_eta^(eff) ).`

In terms of these, the actual selected-branch placement coordinates are

`delta = delta_0 / (1 - eps_eta),`

`M_mix = 8 Z_W (1 + rho)^2 / [ pi^2 (1 - eps_eta) (1 - eps_W) ],`

`R_target = Lambda (1 - eps_eta) (1 - eps_W)^2 / [ Z_W (1 + rho)^2 ],`

with

`Lambda := 27 pi^2 G c_s^5 K_W^(eff) / (20 a^5 c^5 mu_W).`

So the continuum PDE-side placement problem is now fully compressed to one geometric anisotropy ratio, one mixed-sector stability ratio, one interference ratio, one wall-to-mixed overlap ratio, and one radiative demand scale.

Even better, these variables satisfy an exact product law,

`R_target M_mix = 8 Lambda (1 - eps_W) / pi^2`
`               = 54 G c_s^5 K_W^(eff) (1 - eps_W) / (5 a^5 c^5 mu_W).`

That means three apparently independent microscopic lanes only redistribute the defect **along** a fixed product curve, while the mixed-sector stability ratio `eps_W` and the radiative demand scale `Lambda` set the product itself.

So Stage 21 turns the Stage-18/19 admissibility problem into a very small continuum-kernel map.

---

## 1. Exact dimensionless kernel ratios

Keep the Stage-20 continuum effective stiffnesses

`K_eta^(eff) = K_eta + 6 T_Omega,`

`K_W^(eff)   = K_W + pi^2 T_W / (4 L^2).`

Now define the exact dimensionless kernel ratios

`eps_eta := c_(eta U)^2 / ( K_U K_eta^(eff) ),`

`eps_W   := c_(UW)^2 sigma / ( K_U K_W^(eff) ),`

`rho     := c_(UW) c_(eta U) / ( K_U c_(eta W) ),`

`Z_W     := c_(eta W)^2 / ( K_eta^(eff) K_W^(eff) ),`

`delta_0 := pi^2 T_w / ( L^2 K_eta^(eff) ),`

and the radiative demand scale

`Lambda := 27 pi^2 G c_s^5 K_W^(eff) / (20 a^5 c^5 mu_W).`

The exact stability window is now simply

`0 < eps_eta < 1,`

`0 < eps_W   < 1,`

together with the natural nonvanishing transfer branch

`1 + rho != 0.`

---

## 2. Exact continuum placement formulas

Substituting the Stage-20 continuum formulas into the Stage-18/19 branch variables gives

`delta = delta_0 / (1 - eps_eta),`

`M_mix = 8 Z_W (1 + rho)^2 / [ pi^2 (1 - eps_eta) (1 - eps_W) ],`

`R_target = Lambda (1 - eps_eta) (1 - eps_W)^2 / [ Z_W (1 + rho)^2 ].`

So the selected branch is now placed by only five exact dimensionless kernel ratios.

There is also a useful corollary for the outgoing transfer factor:

`beta_0`
`= (mu_W / mu_eta) [ K_eta^(eff) / K_W^(eff) ]`
`  Z_W (1 + rho)^2 / (1 - eps_W)^2.`

So among the inertial masses only the ratio `mu_W / mu_eta` survives in the conservative-to-outgoing transfer factor, while the geometry coordinate `delta` and the mixed baseline `M_mix` are purely stiffness/coupling ratios.

---

## 3. Exact product relation

Multiplying the exact placement formulas yields the simplest structural identity of the stage:

`R_target M_mix`
`= 8 Lambda (1 - eps_W) / pi^2`
`= 54 G c_s^5 K_W^(eff) (1 - eps_W) / (5 a^5 c^5 mu_W).`

Equivalently,

`R_target M_mix = N_Q^(target) Delta_0 / Omega_U^2.`

This is important because the wall-U dressing `eps_eta`, the wall-to-mixed overlap `Z_W`, and the interference ratio `rho` cancel completely from the product.
Those quantities do still matter, but they only redistribute the physical defect along a fixed product curve.

So the continuum kernel separates into two distinct tasks:

1. the pair `(eps_W, Lambda)` sets the product scale,
2. the trio `(eps_eta, Z_W, rho)` decides where on that product curve the actual defect lands.

That is the cleanest factorization seen so far in the moving-throat branch-placement problem.

---

## 4. Exact one-way parameter tendencies

The dimensionless continuum map also makes the monotone tendencies completely explicit.

### 4.1 Wall-U dressing

`d delta / d eps_eta = delta_0 / (1 - eps_eta)^2 > 0,`

`d M_mix / d eps_eta = 8 Z_W (1 + rho)^2 / [ pi^2 (1 - eps_eta)^2 (1 - eps_W) ] > 0,`

`d R_target / d eps_eta = - Lambda (1 - eps_W)^2 / [ Z_W (1 + rho)^2 ] < 0.`

So stronger wall-U softening pushes the defect toward larger anisotropy and larger mixed baseline while lowering the normalization demand ratio.

### 4.2 Internal mixed blocking

`d M_mix / d eps_W = 8 Z_W (1 + rho)^2 / [ pi^2 (1 - eps_eta) (1 - eps_W)^2 ] > 0,`

`d R_target / d eps_W = - 2 Lambda (1 - eps_eta) (1 - eps_W) / [ Z_W (1 + rho)^2 ] < 0.`

So stronger `U/W` blocking simultaneously raises the mixed baseline and lowers the normalization demand ratio.
It is the only stability ratio that also enters the exact product relation.

### 4.3 Wall-to-mixed overlap strength

`d M_mix / d Z_W = 8 (1 + rho)^2 / [ pi^2 (1 - eps_eta) (1 - eps_W) ] > 0,`

`d R_target / d Z_W = - Lambda (1 - eps_eta) (1 - eps_W)^2 / [ Z_W^2 (1 + rho)^2 ] < 0.`

So increasing the wall-to-mixed overlap pushes the physical point upward in `M_mix` and downward in `R_target`.

### 4.4 Interference ratio

`d M_mix / d rho = 16 Z_W (1 + rho) / [ pi^2 (1 - eps_eta) (1 - eps_W) ],`

`d R_target / d rho = - 2 Lambda (1 - eps_eta) (1 - eps_W)^2 / [ Z_W (1 + rho)^3 ].`

So on the natural nonvanishing transfer branch `1 + rho > 0`, constructive interference increases `M_mix` and decreases `R_target`, while destructive interference does the opposite.

---

## 5. Three-lane factorization of the continuum placement problem

The continuum map now splits cleanly into three structural lanes.

### 5.1 Geometry lane

`delta = delta_0 / (1 - eps_eta)`

depends only on the bare wall anisotropy and the wall-U softening ratio.
So the bare geometry lane is controlled entirely by the wall sector plus the brane-like internal doublet.

### 5.2 Mixed-stability/product lane

`R_target M_mix = 8 Lambda (1 - eps_W) / pi^2`

depends only on the radiative demand scale and the internal mixed-sector stability ratio.
So the overall product scale is set in the mixed lane, not by the wall overlap or interference bookkeeping.

### 5.3 Redistribution lane

At fixed product, the trio `(eps_eta, Z_W, rho)` redistributes the point between `R_target` and `M_mix`.
So these parameters control branch placement **along** the product curve but not the curve itself.

This is the sharpest structural decomposition reached so far.

---

## 6. Best current theorem gate after Stage 21

The moving-throat selected quadrupole branch is now described by two exact layers.

### Layer 1 — universal branch geometry

From Stages 18–19, the stable selected branch is controlled by the universal functions

`R_target = F(xi,delta),`

`M_mix <= G(xi,delta).`

### Layer 2 — continuum placement map

From Stages 20–21, the continuum operator places the actual defect at

`delta = delta_0 / (1 - eps_eta),`

`M_mix = 8 Z_W (1 + rho)^2 / [ pi^2 (1 - eps_eta) (1 - eps_W) ],`

`R_target = Lambda (1 - eps_eta) (1 - eps_W)^2 / [ Z_W (1 + rho)^2 ].`

So the remaining theorem gap is now very narrow.

It is no longer:

- “derive more microscopic couplings somehow,”
- or “guess which branch the PDE lands on.”

It is:

> compute the dimensionless kernel ratios `(eps_eta, eps_W, rho, Z_W, delta_0, Lambda)` from the completed moving-throat PDE and check whether the resulting point lies inside the exact Stage-18/19 admissible region.

=== moving_throat_pde_stage022_split_u_sector.md ===

# Moving-Throat PDE — Stage 22: First Non-Flat `U` Doublet, Exact Split Continuum Map, and the Direction-Splitting Theorem

## Purpose

Stage 21 compressed the selected quadrupole placement problem to the dimensionless continuum ratios

`eps_eta, eps_W, rho, Z_W, delta_0, Lambda`,

but it still relied on one structural simplification:

> the internal `U` sector was taken in the flat-doublet limit, so the first two N/N `U` channels were degenerate.

That was the right minimal starting point, but it is not the first genuinely nontrivial continuum operator.

The next exact question is therefore:

> what survives, and what breaks, once the first axial structure of the internal `U` sector is turned on?

This stage answers that question exactly.

The main result is that the **scalar placement map survives**, but the **directional simplicity does not**.

More concretely:

1. the direct wall softening becomes mode-dependent,
2. the bare anisotropy ratio shifts to a new exact value `delta_split`,
3. the mixed-sector blocking ratio becomes an exact split quantity `eps_W_split`,
4. the mixed loading vector rotates away from the source/support direction,
5. and the old Stage-21 one-direction picture survives **iff** the new direction-splitting invariant vanishes.

So the flat-`U` assumption was not harmless bookkeeping. It was precisely what made source, support, and mixed loading all point in the same wall-basis direction.

---

## 1. Turning on the first axial `U` structure

Keep the wall, `phi`, and `W` sectors exactly as in Stage 20, but replace the flat internal `U` block by the first nontrivial N/N continuum operator

`L_U = (mu_U/2) dot(U)^2 - (T_U/2) (U')^2 - (K_U/2) U^2.`

Then the first two N/N internal modes have exact stiffnesses

`K_(U0) = K_U,`

`K_(U1) = K_U + pi^2 T_U / L^2 = K_U (1 + delta_U),`

where the new internal axial-splitting ratio is

`delta_U := pi^2 T_U / (L^2 K_U).`

So the flat-doublet Stage-20 limit is exactly the special case

`delta_U = 0.`

The wall basis itself is unchanged:

`u_0(s) = 1/sqrt(L),`

`u_1(s) = sqrt(2/L) cos(pi s/L),`

and the same D/N overlap data remain

`kappa_0 = 2 sqrt(2)/pi,`

`kappa_1 = -4/(3 pi),`

`sigma = kappa_0^2 + kappa_1^2 = 88/(9 pi^2),`

`lambda_0 := kappa_1^2 / kappa_0^2 = 2/9.`

---

## 2. Exact split direct softening and shifted anisotropy

The wall-U coupling is still diagonal in the N/N basis, but the denominators are no longer the same. The direct softening therefore becomes

`A_0 = [K_eta^(eff) - c_(etaU)^2 / K_U] / mu_eta,`

`A_1 = [K_eta^(eff)(1 + delta_0) - c_(etaU)^2 / K_(U1)] / mu_eta,`

with

`K_eta^(eff) := K_eta + 6 T_Omega,`

`delta_0 := pi^2 T_w / (L^2 K_eta^(eff)),`

`eps_eta := c_(etaU)^2 / (K_U K_eta^(eff)).`

After exact rearrangement,

`A_0 = K_eta^(eff) (1 - eps_eta) / mu_eta,`

`A_1 = A_0 (1 + delta_split),`

with the exact shifted anisotropy ratio

`delta_split = [ delta_0 + eps_eta delta_U / (1 + delta_U) ] / (1 - eps_eta).`

So internal axial structure raises the bare selected-branch anisotropy even if the wall sector itself is unchanged.

For small `delta_U`,

`delta_split = delta_0/(1 - eps_eta) + [eps_eta/(1 - eps_eta)] delta_U + O(delta_U^2).`

---

## 3. Exact split mixed blocking ratio

The `U/W` block also feels the split directly through the overlap-weighted inverse kernel

`S_U = kappa_0^2 / K_U + kappa_1^2 / K_(U1).`

Using the Stage-21 flat-doublet ratio

`eps_W := c_(UW)^2 sigma / (K_U K_W^(eff)),`

with

`K_W^(eff) := K_W + pi^2 T_W / (4 L^2),`

one finds the exact split blocking ratio

`eps_W_split = eps_W [ 1 - (2/11) delta_U / (1 + delta_U) ].`

So the first axial `U` structure lowers the mixed blocking ratio relative to the flat-doublet value.

For small `delta_U`,

`eps_W_split = eps_W [ 1 - (2/11) delta_U ] + O(delta_U^2).`

---

## 4. Exact mixed loading vector and the direction-splitting theorem

The decisive new object is the mixed loading vector.

Define the flat-mode interference ratio

`rho_0 := c_(UW) c_(etaU) / (K_U c_(etaW)).`

Then the two wall-basis components of the mixed loading vector are

`z_0 = kappa_0 g_W (1 + rho_0),`

`z_1 = kappa_1 g_W (1 + rho_0/(1 + delta_U)),`

where

`g_W = c_(etaW) / sqrt(mu_eta mu_W).`

Equivalently,

`z_1 / z_0 = (kappa_1 / kappa_0) R_U,`

with the exact direction-ratio factor

`R_U := [ 1 + rho_0/(1 + delta_U) ] / (1 + rho_0).`

This isolates the first real failure mode of the flat-`U` simplification.

In Stage 20–21 the mixed loading direction was proportional to the source/support vector `v = (kappa_0,kappa_1)^T`.
Once `delta_U != 0`, the mixed loading vector is instead proportional to

`z = ( z_0, z_1 )^T,`

which is generally **not collinear** with `v`.

The exact direction-splitting invariant is

`D_dir := kappa_0 z_1 - kappa_1 z_0`

`      = - kappa_0 kappa_1 g_W rho_0 delta_U / (1 + delta_U).`

So the exact collinearity theorem is

`D_dir = 0  <=>  delta_U = 0  or  rho_0 = 0,`

which means:

- flat internal `U` doublet, or
- zero `U/W` interference,

are the only ways to keep the old one-direction Stage-21 geometry exactly intact.

For small `delta_U`,

`R_U = 1 - [rho_0/(1 + rho_0)] delta_U + O(delta_U^2).`

So on the natural constructive branch `rho_0 > 0`, the first axial `U` structure pushes the mixed loading direction away from the flat-doublet direction in a controlled linear way.

---

## 5. Exact split continuum placement map

Even though the directions split, the scalar placement data themselves still factorize cleanly.

Using the Stage-21 dimensionless overlap ratio

`Z_W := c_(etaW)^2 / ( K_eta^(eff) K_W^(eff) ),`

and the same radiative demand scale

`Lambda := 27 pi^2 G c_s^5 K_W^(eff) / (20 a^5 c^5 mu_W),`

the exact split placement formulas are

`M_mix^(split U)`
`= 8 Z_W (1 + rho_0)^2`
`  / [ pi^2 (1 - eps_eta) (1 - eps_W_split) ],`

`R_target^(split U)`
`= Lambda (1 - eps_eta) (1 - eps_W_split)^2`
`  / [ Z_W (1 + rho_0)^2 ].`

And the exact product law survives:

`R_target^(split U) M_mix^(split U)`
`= 8 Lambda (1 - eps_W_split) / pi^2.`

So the Stage-21 factorization survives at the scalar-placement level.
What changes is the **directional geometry** seen by the selected branch.

For small `delta_U`,

`M_mix^(split U) = M_mix^(flat) [ 1 - 2 eps_W delta_U / (11(1 - eps_W)) ] + O(delta_U^2),`

`R_target^(split U) = R_target^(flat) [ 1 + 4 eps_W delta_U / (11(1 - eps_W)) ] + O(delta_U^2).`

So positive internal axial splitting lowers the mixed baseline and raises the normalization demand ratio before any selected-branch/source-map correction is even considered.

---

## 6. Best current theorem statement after Stage 22

The first non-flat `U` doublet does **not** destroy the continuum placement map.
It does something subtler and more important.

It separates two statements that were accidentally fused in Stages 20–21:

1. **scalar placement factorization**, which survives, and
2. **source/loading collinearity**, which does not survive generically.

That means the real structural role of the flat-`U` assumption is now clear.
It was not merely making the formulas shorter. It was enforcing the coincidence of

- the source map direction,
- the support direction,
- and the mixed loading direction.

Once the first axial `U` structure is turned on, the next theorem problem is no longer “what are the scalar placement ratios?”
Those are still exact.

The next theorem problem is:

> how does the selected-branch normalization law deform when the source vector and the loading vector are no longer the same?

That is the target of Stage 23.

=== moving_throat_pde_stage023_generalized_selected_branch.md ===

# Moving-Throat PDE — Stage 23: Generalized Selected-Branch Normalization with Source/Loading Mismatch

## Purpose

Stage 22 showed that the first axial structure of the internal `U` sector does **not** destroy the scalar continuum placement map, but it **does** generically rotate the mixed loading vector away from the source/support direction.

That means the old Stage-18/19 branch functions cannot be assumed blindly once `delta_U != 0`.

The right next question is therefore:

> how does the selected-branch normalization law deform when the source vector and the loading vector are no longer the same?

This stage answers that exactly.

The main result is that the selected-branch geometry survives, but it becomes a one-parameter deformation of the flat-`U` branch.

More concretely:

1. the selected lower wall mode can still be solved exactly,
2. the old normalization function `F(xi,delta)` is replaced by a two-vector function `F_(q,eta)(xi,delta)`,
3. the required baseline loading remains one-dimensional through a deformed function `G_q(xi,delta)`,
4. for the split-`U` continuum the whole deformation collapses to one exact ratio `R_U`,
5. and setting `R_U = 1` reproduces the Stage-18/19 flat-`U` branch exactly.

So the first non-flat `U` structure does not kill the theorem geometry — it deforms it in a controlled way.

---

## 1. Exact rank-1 loaded 2x2 branch solve

Consider the selected wall basis after the Stage-22 direct `U` softening has already been absorbed into the diagonal baseline matrix

`K_base = diag(A_0, A_0(1 + delta)).`

Now add a single directional loading

`K_loaded = K_base - alpha z z^T,`

with loading vector

`z = (z_0, z_1)^T.`

Write the lower selected eigenvalue as

`lambda_- = A_0 (1 - xi),`

with stable branch parameter

`0 <= xi < 1.`

Then the exact required loading is

`alpha_req = A_0 xi (xi + delta) / [ z_0^2 ( delta + (1 + q^2) xi ) ],`

where

`q := z_1 / z_0.`

So the basic one-parameter selected-branch geometry survives exactly even before any source map is specified.

The exact lower-mode eigenvector ratio is

`e_1 / e_0 = q xi / (delta + xi).`

That is the universal 2x2 rank-1-loaded backbone behind everything that follows.

---

## 2. Exact overlap formulas with distinct source and loading directions

Now separate two vectors that coincided in the flat-`U` stages:

- the **loading vector** `z`,
- the **source vector** `s`.

Define the signed mismatch parameter

`eta := (s_1 z_1) / (s_0 z_0).`

Then the exact normalized overlaps of the selected lower mode are

`(z.e_-)^2 / z_0^2`
`= [ delta + (1 + q^2) xi ]^2`
`  / [ (delta + xi)^2 + q^2 xi^2 ],`

`(s.e_-)^2 / s_0^2`
`= [ delta + (1 + eta) xi ]^2`
`  / [ (delta + xi)^2 + q^2 xi^2 ].`

So the selected-branch normalization product is no longer governed by the old single-vector function.
It is governed by the exact two-vector shape function

`F_(q,eta)(xi,delta)`
`= [ delta + (1 + q^2) xi ]^2 [ delta + (1 + eta) xi ]^2`
`  / [ (1 - xi) ( (delta + xi)^2 + q^2 xi^2 )^2 ].`

At the same time, the required baseline loading stays one-dimensional:

`G_q(xi,delta) = xi (xi + delta) / [ delta + (1 + q^2) xi ].`

This is the exact generalization of the Stage-18/19 pair `(F,G)`.

---

## 3. Specialization to the split-`U` continuum

For the split-`U` continuum of Stage 22, the source vector is still the original D/N overlap direction

`s = v = (kappa_0, kappa_1)^T,`

while the loading vector obeys

`z_1 / z_0 = (kappa_1 / kappa_0) R_U.`

Using

`lambda_0 := kappa_1^2 / kappa_0^2 = 2/9,`

this means

`q = - sqrt(lambda_0) R_U = - (sqrt(2)/3) R_U,`

and therefore

`eta = (s_1 z_1)/(s_0 z_0) = lambda_0 R_U = (2/9) R_U.`

So the full source/loading mismatch collapses to one exact parameter `R_U`.

Substituting that into the general formulas gives the exact split-`U` branch functions

`F_U(xi,delta;R_U)`
`= [ 9 delta + (9 + 2 R_U^2) xi ]^2 [ 9 delta + (9 + 2 R_U) xi ]^2`
`  / [ 81 (1 - xi) ( 9 delta^2 + 18 delta xi + (9 + 2 R_U^2) xi^2 )^2 ],`

`G_U(xi,delta;R_U)`
`= 9 xi (xi + delta) / [ 9 delta + (9 + 2 R_U^2) xi ].`

This is the exact one-parameter deformation of the old flat-`U` selected-branch geometry.

---

## 4. Exact recovery of the Stage-18/19 flat-`U` branch

If `R_U = 1`, then `q = kappa_1/kappa_0` and the source/loading mismatch disappears.
The exact branch functions collapse back to the Stage-18/19 formulas:

`F_U(xi,delta;1)`
`= (9 delta + 11 xi)^4 / [ 81 (1 - xi) (9 delta^2 + 18 delta xi + 11 xi^2)^2 ],`

`G_U(xi,delta;1)`
`= 9 xi (xi + delta) / (9 delta + 11 xi ).`

So the flat-`U` selected-branch theorem is not replaced. It is recovered as the exact collinear limit of the more general split-`U` theory.

---

## 5. Exact small-deformation expansion about the flat-`U` branch

Write

`R_U = 1 + eps,`

with `eps` small.

Then the exact normalized deformation of the selected-branch functions is

`F_U / F_flat = 1 + eps H_F + O(eps^2),`

`G_U / G_flat = 1 + eps H_G + O(eps^2),`

with

`H_F`
`= 4 xi ( 27 delta^2 + 36 delta xi + 11 xi^2 )`
`  / [ (9 delta + 11 xi) (9 delta^2 + 18 delta xi + 11 xi^2) ],`

`H_G = - 4 xi / (9 delta + 11 xi).`

So the deformation is smooth and exact.

Now combine this with the Stage-22 result

`R_U = 1 - [ rho_0 / (1 + rho_0) ] delta_U + O(delta_U^2).`

On the natural constructive branch `rho_0 > 0`, positive internal axial splitting implies `R_U < 1`, so:

- `F_U` is shifted **below** the flat-`U` normalization function at fixed `(xi,delta)`,
- `G_U` is shifted **above** the flat-`U` baseline-loading function at fixed `(xi,delta)`.

That is the first exact robustness statement for the selected branch once source/loading collinearity is relaxed.

---

## 6. Best current theorem statement after Stage 23

The exact theorem picture is now much sharper.

### What survives intact

- the Stage-22 scalar continuum placement map,
- the existence of a one-parameter selected branch `xi`,
- the exact lower-eigenvalue softening law,
- and the exact recovery of the Stage-18/19 functions when the split-`U` deformation is removed.

### What changes structurally

- the selected-branch normalization no longer depends on one vector only,
- the Stage-18 function `F(xi,delta)` is replaced by `F_(q,eta)(xi,delta)`,
- and in the physical split-`U` continuum this becomes a one-parameter deformation `F_U(xi,delta;R_U)`.

### The new bottleneck

The flat-`U` program reduced the selected branch to one collinear direction shared by

- source,
- support,
- and mixed loading.

Stages 22–23 show that the first non-flat `U` structure breaks that coincidence generically.

So the next honest derivation step is no longer “compute more scalar placement ratios.”
Those are already exact.

The next step is:

> determine how the additional support/BdG loading enters the now noncollinear selected-branch geometry, and whether the physical support direction tracks the deformed loading vector or remains tied to the original source vector.

That is the first place where the old Stage-19 support-feasibility theorem may need a true rank-2 extension.

=== moving_throat_pde_stage024_rank2_support_completion.md ===

# Moving-Throat PDE — Stage 24: Rank-2 Support Completion and the Exact Support-Loading Theorem

## Purpose

Stage 23 showed that once the first axial structure of the internal `U` sector is turned on, the mixed loading vector `z` is no longer generically collinear with the source/support direction `v`.

That leaves one honest reduced-theorem bottleneck:

> what happens when the selected wall branch is loaded by **two** directions rather than one?

This stage answers that exactly.

The main result is that the rank-2 completion is still analytically tractable because the determinant stays **linear** in the support loading. That gives an exact support-loading theorem.

More concretely:

1. the selected lower branch with mixed loading `z` and support loading `y` is still solvable in closed form,
2. the exact required support loading is a rational function `n_req(xi,delta;m,q,r)`,
3. increasing the mixed baseline always lowers the support needed to hit a given selected branch,
4. if the support tracks the deformed mixed direction, the whole rank-2 problem collapses back to Stage 23,
5. if the support stays tied to the original source direction, the branch geometry becomes genuinely new, with an exact support-feasibility window.

So the flat-`U` simplification did not just hide a numerical correction. It hid the first place where the support direction itself becomes a real theorem variable.

---

## 1. Exact rank-2 loaded selected branch

Work in the 2D wall basis of Stages 22–23 and write the diagonal baseline as

`K_base = A_0 diag(1, 1 + delta).`

Now allow two independent rank-1 loadings:

- a mixed baseline loading carried by `z`,
- a support/BdG loading carried by `y`.

Normalize the first wall-basis component of each direction and define

`z = z_0 (1, q)^T,`

`y = y_0 (1, r)^T,`

with dimensionless loadings

`m := alpha_mix z_0^2 / A_0,`

`n := alpha_supp y_0^2 / A_0.`

The full selected-branch matrix is therefore

`K_loaded / A_0 = diag(1,1+delta) - m (1,q)^T(1,q) - n (1,r)^T(1,r).`

As before, parameterize the lower selected eigenvalue by the softening depth

`lambda_- = A_0 (1 - xi),`

with stable branch parameter

`0 <= xi < 1.`

Then the exact determinant condition is

`D_sel = xi (delta + xi)`
`        - m [ delta + (1 + q^2) xi ]`
`        - n [ delta + (1 + r^2) xi ]`
`        + m n (q - r)^2`
`      = 0.`

This is the structural point of the stage:

> even with two loading directions, the determinant remains **linear** in the support loading `n`.

So the rank-2 completion is not a high-dimensional numerical problem yet. It is still an exact reduced theorem.

---

## 2. Exact support-loading theorem

Solving the determinant condition for the support loading gives

`n_req(xi,delta;m,q,r)`
`= [ xi(delta + xi) - m( delta + (1 + q^2) xi ) ]`
`  / [ delta + (1 + r^2) xi - m (q - r)^2 ].`

This is the exact rank-2 support-loading theorem.

The denominator shows the first genuinely new effect of noncollinearity:

`delta + (1 + r^2) xi - m (q - r)^2.`

If the two directions differ, the mixed baseline modifies the *support* denominator directly. That effect was absent in the collinear Stage-18/19/23 branch.

A second exact theorem follows immediately by differentiation:

`d n_req / d m`
`= - [ delta + (1 + q r) xi ]^2`
`  / [ delta + (1 + r^2) xi - m (q - r)^2 ]^2`
`< 0.`

So, whenever the branch is regular, **increasing the mixed baseline always lowers the support loading needed to reach the same softening depth**.

That monotonicity is exact and independent of the detailed direction mismatch.

---

## 3. Tracking theorem: support follows the deformed mixed direction

The first natural hypothesis is

`y || z`,

so that

`r = q.`

Then the exact support theorem collapses immediately to

`n_req = xi (delta + xi) / [ delta + (1 + q^2) xi ] - m`
`      = G_q(xi,delta) - m.`

This is exactly the old Stage-19/23 geometry with a split between

- mixed baseline loading `m`,
- and additional support loading `n`.

So the first sharp theorem is:

> **If the physical support/BdG loading follows the deformed mixed direction, the rank-2 completion introduces no new selected-branch geometry.**

The whole problem collapses back to the Stage-23 one-parameter deformation.

That makes `y || z` a very strong closure hypothesis.

---

## 4. Source-tied theorem: support remains aligned with the original source direction

The second natural hypothesis is the opposite one:

`y || v,`

where `v = (kappa_0, kappa_1)^T` is the original D/N source/support direction.

Write the source ratio as

`t := kappa_1 / kappa_0,`

so that

`t^2 = lambda_0 = 2/9.`

From Stage 22, the mixed vector obeys

`q = t R_U,`

where `R_U` is the exact split-`U` direction ratio. Under the source-tied hypothesis,

`r = t.`

The exact required support loading becomes

`n_req^(src)`
`= [ xi(delta + xi) - m( delta + (1 + lambda_0 R_U^2) xi ) ]`
`  / [ delta + (1 + lambda_0) xi - m lambda_0 (R_U - 1)^2 ].`

This is the first genuinely new branch formula beyond Stage 23.

Two exact support-feasibility conditions follow.

### Regularity condition

The branch stays regular only if

`delta + (1 + lambda_0) xi - m lambda_0 (R_U - 1)^2 > 0,`

or equivalently

`m < [ delta + (1 + lambda_0) xi ] / [ lambda_0 (R_U - 1)^2 ]`

for `R_U != 1`.

### Positive-support condition

Assuming the regularity denominator is positive, nonnegative support loading requires

`m <= [ xi(delta + xi) ] / [ delta + (1 + lambda_0 R_U^2) xi ].`

So the source-tied rank-2 branch has a sharp mixed-loading ceiling before the support channel can even remain physical.

That ceiling disappears in the collinear flat-`U` limit `R_U = 1`, where the denominator correction vanishes and the branch collapses back to Stage 19.

---

## 5. Exact comparison of the two hypotheses

The rank-2 support completion therefore splits into two sharply testable reduced hypotheses.

### Hypothesis A — support tracks the mixed vector

`y || z`

Then

`n_req = G_q - m,`

and no new branch geometry appears.

### Hypothesis B — support stays tied to the original source vector

`y || v`

Then

`n_req = n_req^(src)`

with the exact denominator correction

`- m lambda_0 (R_U - 1)^2`.

This is the first place where the noncollinearity of Stage 22 becomes a true structural effect instead of just a deformation of Stage-23 constants.

So the exact theorem picture is now:

- **tracking support** preserves the Stage-23 geometry,
- **source-tied support** creates a new rank-2 support-feasibility problem.

That is the sharpest reduced statement yet of the support-direction bottleneck.

---

## 6. Best current theorem statement after Stage 24

The rank-2 support completion is no longer vague.

### What is now exact

- the selected-branch determinant with two loading directions,
- the exact required support loading `n_req`,
- the exact monotonicity theorem `d n_req / d m < 0`,
- the exact collapse to Stage 23 when `y || z`,
- and the exact source-tied branch formula when `y || v`.

### What remains open

The unresolved PDE-side question is now very specific:

> does the physical support/BdG kernel align with the deformed mixed vector `z`, or does it remain tied to the original source direction `v`?

If it aligns with `z`, the Stage-23 one-parameter deformation is already the whole selected-branch geometry.
If it stays tied to `v`, the branch acquires a new exact rank-2 support-feasibility window.

That is the next theorem gate the completed moving-throat operator has to decide.

=== moving_throat_pde_stage025_rank2_selected_mode_normalization.md ===

# Moving-Throat PDE — Stage 25: Selected-Mode Normalization Under Rank-2 Support Completion

## Purpose

Stage 24 closed the rank-2 *existence* problem for the selected wall branch: with mixed loading `z` and support loading `y`, the exact support load needed to reach a chosen softening depth is known.

But the 2.5PN bridge does not depend only on branch existence. It depends on the **selected-mode normalization product** that combines

- the outgoing mixed overlap,
- the source-map overlap,
- and the selected stiffness.

So the right next question is:

> how does the selected-mode normalization law deform once the branch is supported by two different directions?

This stage answers that exactly.

The main result is that the selected-mode geometry remains closed-form, but the normalization now depends on **three** directions:

- the mixed/outgoing direction `z`,
- the support direction `y`,
- the source direction `s`.

That gives an exact generalized normalization function `F_(q,r,t)(xi,delta;m)`. The old Stage-23 law is recovered as the special case `r=q`.

So the rank-2 bottleneck is no longer conceptual. It is now a sharply defined source/support/outgoing alignment problem.

---

## 1. Exact selected-mode geometry with two loading directions

Carry forward the Stage-24 wall operator

`K_loaded / A_0 = diag(1,1+delta) - m (1,q)^T(1,q) - n (1,r)^T(1,r),`

with the lower selected eigenvalue parameterized by

`lambda_- = A_0 (1 - xi),`

and with the support loading fixed to the exact Stage-24 value

`n = n_req(xi,delta;m,q,r).`

Then the lower selected eigenvector ratio is still exact:

`e_1 / e_0`
`= [ m(q-r) + r xi ] / [ delta + xi - m q(q-r) ].`

So the two-direction selected branch is still described by one scalar deformation variable `xi`, but the eigenvector now depends explicitly on the baseline mixed loading `m` whenever the support direction differs from the mixed direction.

That dependence was absent in the Stage-23 one-direction problem.

---

## 2. Exact overlap formulas

Now introduce a general source vector

`s = s_0 (1,t)^T.`

Using the exact selected eigenvector above, the normalized overlaps are

`(z.e_-)^2 / z_0^2`
`= [ delta + (1 + q r) xi ]^2`
`  / D_(q,r)(xi,delta;m),`

`(s.e_-)^2 / s_0^2`
`= [ delta + (1 + r t) xi - m (q-r)(q-t) ]^2`
`  / D_(q,r)(xi,delta;m),`

with the exact denominator

`D_(q,r)(xi,delta;m)`
`= [ delta + xi - m q(q-r) ]^2 + [ m(q-r) + r xi ]^2.`

So the normalization data no longer factor through one direction only. The mixed/outgoing overlap and the source overlap are deformed differently once `q != r`.

---

## 3. Exact rank-2 normalization function

The selected-branch normalization product is therefore

`F_(q,r,t)(xi,delta;m)`
`= [ delta + (1 + q r) xi ]^2`
`  [ delta + (1 + r t) xi - m(q-r)(q-t) ]^2`
`  / [ (1 - xi) D_(q,r)(xi,delta;m)^2 ].`

This is the exact rank-2 generalization of the Stage-23 selected-mode normalization law.

It is the object that the outgoing quadrupole bridge now depends on whenever the support direction and the mixed direction differ.

---

## 4. Exact collapse to Stage 23 when support tracks the mixed vector

If the support tracks the mixed direction,

`r = q,`

then

`D_(q,q) = (delta + xi)^2 + q^2 xi^2,`

and the exact rank-2 normalization function collapses to

`F_(q,q,t)(xi,delta;m)`
`= [ delta + (1 + q^2) xi ]^2 [ delta + (1 + q t) xi ]^2`
`  / [ (1 - xi) ((delta + xi)^2 + q^2 xi^2)^2 ].`

So the second sharp theorem is:

> **If the support/BdG loading follows the deformed mixed direction, the entire rank-2 normalization law collapses exactly to the Stage-23 two-vector function.**

In that case the rank-2 completion adds no new normalization geometry.

---

## 5. Source-tied specialization for the physical split-`U` continuum

Now take the physically interesting opposite hypothesis:

- support direction tied to the original source vector,
- source vector itself also equal to that original D/N direction.

Write

`t := kappa_1 / kappa_0,`

so that

`t^2 = lambda_0 = 2/9.`

From Stage 22 the mixed direction obeys

`q = t R_U,`

while the source-tied support hypothesis gives

`r = t.`

Then the exact normalization function becomes

`F_src(xi,delta;m,R_U)`
`= [ delta + (1 + lambda_0 R_U) xi ]^2`
`  [ delta + (1 + lambda_0) xi - m lambda_0 (R_U - 1)^2 ]^2`
`  / [ (1 - xi) D_src(xi,delta;m,R_U)^2 ],`

with

`D_src(xi,delta;m,R_U)`
`= [ delta + xi - m lambda_0 R_U (R_U - 1) ]^2`
`  + lambda_0 [ xi + m (R_U - 1) ]^2.`

This is the first exact selected-mode normalization law beyond Stage 23 that depends explicitly on the mixed baseline loading `m`.

That dependence is the clean reduced signature of the source/support/outgoing mismatch.

---

## 6. Exact flat-`U` recovery and first-order deformation

When `R_U = 1`, all three directions become collinear and the source-tied rank-2 completion collapses exactly to the old flat-`U` law:

`F_src(xi,delta;m,1) = F_flat(xi,delta),`

where

`F_flat(xi,delta)`
`= [ delta + (1 + lambda_0) xi ]^4`
`  / [ (1 - xi) ( (delta + xi)^2 + lambda_0 xi^2 )^2 ].`

So the source-tied rank-2 branch is a genuine deformation of Stage 23, not a different disconnected object.

Write

`R_U = 1 + eps.`

Then the exact first-order support-loading deformation is

`n_src = G_flat(xi,delta) - m + eps H_n^(src) + O(eps^2),`

with

`H_n^(src)`
`= - 2 lambda_0 m xi / [ delta + (1 + lambda_0) xi ].`

Likewise the exact first-order normalization deformation is

`F_src / F_flat = 1 + eps H_F^(src) + O(eps^2),`

with

`H_F^(src)`
`= 2 lambda_0`
`  [ xi ( (delta + xi)^2 + lambda_0 xi^2 )`
`    + 2 m delta ( delta + (1 + lambda_0) xi ) ]`
`  / [ ( delta + (1 + lambda_0) xi ) ( (delta + xi)^2 + lambda_0 xi^2 ) ].`

On the constructive split-`U` branch of Stage 22, one has `R_U < 1`, i.e. `eps < 0`, so:

- the source-tied support requirement is **raised** above the simple flat subtraction `G_flat - m`,
- the selected-mode normalization function is **lowered** below the flat-`U` value.

So the source-tied hypothesis makes the selected-branch normalization test strictly harder on the natural split-`U` branch.

---

## 7. Best current theorem statement after Stage 25

The rank-2 support bottleneck is now completely explicit.

### If support tracks the mixed vector

- Stage 24 gives `n_req = G_q - m`,
- Stage 25 gives exact collapse back to the Stage-23 normalization law.

So the old one-parameter deformation survives intact.

### If support remains tied to the original source direction

- Stage 24 gives the new exact source-tied support-feasibility formula,
- Stage 25 gives the new exact source-tied normalization function `F_src`.

So the selected branch becomes a true two-layer deformation:

1. Stage 22 deforms the mixed direction through `R_U`,
2. Stages 24–25 decide whether the support kernel follows that deformation or resists it.

The open PDE-side theorem question is therefore no longer diffuse:

> compute the actual projected support/BdG loading direction from the moving-throat operator and determine whether the physical kernel lands on the tracking closure or the source-tied closure.

That is now the sharp reduced theorem gate for the next stage.

=== moving_throat_pde_stage026_support_direction_extraction.md ===

# Moving-Throat PDE — Stage 26: Continuum Extraction of the Actual Support/BdG Loading Direction

## Purpose

Stage 24 and Stage 25 reduced the next bottleneck to one sharp question:

> does the physical support/BdG kernel follow the deformed mixed direction `z`, or does it stay tied to the original source direction `v`?

The minimal Stage-20 continuum operator implicitly answered that only in the most degenerate case: because the support mode `phi` coupled directly to the wall through the same overlap vector `v`, the support direction was exactly source-tied.

But that was still not the first genuinely nontrivial support operator.

The first symmetry-allowed extension is to turn on the bilinear `U/phi` coupling. Once that is done, the support direction can be extracted exactly from the continuum kernel rather than guessed.

The main result is:

1. the effective support/BdG loading remains rank-1,
2. its actual wall-basis direction is an exact rotated vector `y`,
3. the rotation is controlled by one new interference ratio `sigma_0`,
4. the source-tied closure is the exact minimal-kernel limit `sigma_0 = 0`,
5. the tracking closure is an exact codimension-one interference-match condition,
6. and the generic first extended continuum kernel lands on an **intermediate support direction**, neither source-tied nor tracking.

So the Stage-24 binary fork is now resolved more sharply:

- the **minimal** continuum kernel is source-tied,
- the **first extended** continuum kernel is generically intermediate,
- and exact tracking requires one special interference match.

---

## 1. First symmetry-allowed support extension of the continuum kernel

Keep the Stage-22 split-`U` continuum operator, but add the first local bilinear support coupling

`L_(Uphi) = + c_(Uphi) int_0^L U phi ds.`

After the same N/N and D/N mode projection used in Stages 20–22, the reduced static couplings are

`L_(eta U)   = - g_U Q.U,`

`L_(eta phi) = - g_B (v.Q) Phi,`

`L_(Uphi)    = + g_S (v.U) Phi,`

with wall basis `Q = (Q_0,Q_1)^T`, D/N overlap vector

`v = (kappa_0, kappa_1)^T,`

and the same split internal stiffnesses

`A_(U0) = K_U,`

`A_(U1) = K_U (1 + delta_U).`

So the support sector now sees the wall both

- directly through `g_B v`, and
- indirectly through `U` via `g_U g_S`.

That is the first honest continuum mechanism that can rotate the support direction away from the original source vector.

---

## 2. Exact effective support vector

Eliminate the split `U` doublet first.

The exact effective wall-to-support coupling becomes

`y = g_B v + g_U g_S D_U v,`

with

`D_U = diag( 1/K_U, 1/[K_U(1+delta_U)] ).`

So the wall-basis components are

`y_0 = kappa_0 [ g_B + g_U g_S / K_U ],`

`y_1 = kappa_1 [ g_B + g_U g_S / ( K_U (1+delta_U) ) ].`

Define the exact support-interference ratio

`sigma_0 := g_U g_S / (K_U g_B)`

or equivalently, in continuum couplings,

`sigma_0 = c_(etaU) c_(Uphi) / ( K_U c_(etaphi) ).`

Then the support vector takes the exact factorized form

`y_0 = kappa_0 g_B (1 + sigma_0),`

`y_1 = kappa_1 g_B [ 1 + sigma_0/(1+delta_U) ].`

So the actual support/BdG direction ratio is

`y_1 / y_0 = (kappa_1 / kappa_0) R_phi,`

with the exact support-direction factor

`R_phi := [ 1 + sigma_0/(1+delta_U) ] / (1 + sigma_0).`

This is the exact support analogue of the mixed-direction factor from Stage 22,

`R_U = [ 1 + rho_0/(1+delta_U) ] / (1 + rho_0).`

---

## 3. Exact support direction-splitting theorem

The support direction is collinear with the original source vector `v` iff

`R_phi = 1.`

The exact support direction-splitting invariant is

`D_phi := kappa_0 y_1 - kappa_1 y_0`

`      = - kappa_0 kappa_1 g_B sigma_0 delta_U / (1 + delta_U).`

Therefore

`D_phi = 0  <=>  sigma_0 = 0  or  delta_U = 0.`

So the exact theorem is:

> **The support/BdG loading remains source-tied only in the minimal kernel limit `sigma_0 = 0` (or in the unsplit limit `delta_U = 0`).**

That identifies the precise structural role of the minimal Stage-20 operator: it was not just simple; it was exactly the source-tied special case.

---

## 4. Exact support pole shift and actual support baseline

Eliminating `U` also renormalizes the support pole.

The exact effective support denominator is

`A_phi^(eff) = K_phi^(eff) - g_S^2 v.D_U.v,`

with

`K_phi^(eff) := K_phi + pi^2 T_phi / (4 L^2).`

Using

`sigma = kappa_0^2 + kappa_1^2 = 88 / (9 pi^2),`

and `kappa_1^2/sigma = 2/11`, the overlap contraction is

`v.D_U.v = (sigma / K_U) [ 1 - (2/11) delta_U/(1+delta_U) ].`

So if we define the flat support-blocking ratio

`eps_phi := c_(Uphi)^2 sigma / ( K_U K_phi^(eff) ),`

then the exact split support-blocking ratio is

`eps_phi^(split) = eps_phi [ 1 - (2/11) delta_U/(1+delta_U) ].`

The physical support loading baseline on the selected wall problem is therefore

`M_supp`
`= 8 Z_phi (1 + sigma_0)^2`
`  / [ pi^2 (1 - eps_eta) (1 - eps_phi^(split)) ],`

with

`Z_phi := c_(etaphi)^2 / ( K_eta^(eff) K_phi^(eff) ).`

So the support sector is now no longer just a direction question. Its actual baseline strength is an exact continuum output, parallel to the mixed baseline `M_mix` of Stage 22.

---

## 5. Exact tracking theorem relative to the mixed vector

Stage 22 already gave the mixed loading vector

`z_0 = kappa_0 g_W (1 + rho_0),`

`z_1 = kappa_1 g_W [ 1 + rho_0/(1+delta_U) ].`

So the support vector `y` tracks the mixed vector `z` iff `R_phi = R_U`.

The exact mixed-support collinearity invariant is

`D_(phi z) := y_0 z_1 - y_1 z_0`

`          = - delta_U kappa_0 kappa_1 g_U ( g_B g_R - g_W g_S ) / [ K_U (1+delta_U) ].`

Therefore

`y || z  <=>  g_B g_R = g_W g_S`

or equivalently

`sigma_0 = rho_0.`

So exact tracking is not generic. It is a codimension-one interference match between

- the `eta/W`–`U/W` mixed lane, and
- the `eta/phi`–`U/phi` support lane.

---

## 6. Exact mismatch formula and branch interpretation

The support and mixed direction factors differ by

`R_phi - R_U`
`= delta_U (rho_0 - sigma_0)`
`  / [ (1+delta_U)(1+rho_0)(1+sigma_0) ].`

So the sign of the direction mismatch is set by the sign of `rho_0 - sigma_0`.

This yields the clean continuum interpretation:

- `sigma_0 = 0`  ->  exact source-tied closure,
- `sigma_0 = rho_0`  ->  exact tracking closure,
- generic `sigma_0`  ->  intermediate support direction.

So the physical kernel does **not** generically choose between the two Stage-24 extremes. It generically interpolates between them.

---

## 7. Best current theorem statement after Stage 26

The support-direction bottleneck is no longer open in principle.

### What is now exact

- the first symmetry-allowed continuum extension of the support sector,
- the actual effective support vector `y`,
- the exact support-direction factor `R_phi`,
- the exact support direction-splitting invariant `D_phi`,
- the exact support-blocking ratio `eps_phi^(split)`,
- the exact physical support baseline `M_supp`,
- and the exact tracking condition `g_B g_R = g_W g_S`.

### What this means physically

- The **minimal** continuum kernel lands on the **source-tied** closure.
- The **first extended** continuum kernel is generically **intermediate**.
- Exact **tracking** requires a special interference match.

So the next theorem step is no longer “guess the support direction.”
It is:

> insert the exact continuum-selected quantities `(M_mix, M_supp, R_U, R_phi)` into the Stage-24/25 selected-branch formulas and determine the physical selected branch.

=== moving_throat_pde_stage027_continuum_selected_rank2_closure.md ===

# Moving-Throat PDE — Stage 27: Continuum-Selected Rank-2 Closure and the Exact Quadratic Branch Equation

## Purpose

Stage 26 extracted the actual support/BdG loading data from the continuum operator:

- the support direction factor `R_phi`,
- and the physical support baseline `M_supp`.

That means the Stage-24/25 rank-2 problem is no longer a family of possible closures. The continuum kernel now supplies the data needed to write the **actual physical selected-branch equations**.

This stage does that.

The main result is that the continuum-selected wall branch is pinned by one exact quadratic equation for the softening depth `xi`, and the 2.5PN normalization test then becomes one exact residual equation evaluated on that root.

So the selected-branch problem is now no longer “support-tied or tracking?” in the abstract. It is an explicit continuum-selected closure with two exact special surfaces:

- the **minimal-kernel source-tied surface**, and
- the **interference-matched tracking surface**.

Generically, the physical kernel lands in between them.

---

## 1. Exact continuum-selected branch data

Collect the actual continuum outputs from Stages 22 and 26.

### Mixed lane

`M_mix`
`= 8 Z_W (1 + rho_0)^2`
`  / [ pi^2 (1 - eps_eta) (1 - eps_W^(split)) ],`

`R_U = [ 1 + rho_0/(1+delta_U) ] / (1 + rho_0).`

### Support/BdG lane

`M_supp`
`= 8 Z_phi (1 + sigma_0)^2`
`  / [ pi^2 (1 - eps_eta) (1 - eps_phi^(split)) ],`

`R_phi = [ 1 + sigma_0/(1+delta_U) ] / (1 + sigma_0).`

### Wall/source geometry

The source vector remains the original D/N direction

`s = s_0 (1,t)^T,`

with

`t = kappa_1 / kappa_0,`

`t^2 = lambda_0 = 2/9.`

So the actual mixed and support direction ratios entering the selected-branch problem are

`q = t R_U,`

`r = t R_phi.`

The anisotropy ratio is still the Stage-22 shifted wall value `delta = delta_split`.

---

## 2. Exact physical selected-branch equation

Stage 24 gave the exact support-loading theorem

`n_req(xi,delta;m,q,r)`
`= [ xi(delta + xi) - m( delta + (1 + q^2) xi ) ]`
`  / [ delta + (1 + r^2) xi - m (q-r)^2 ].`

The actual physical branch is obtained by setting

`m = M_mix,`

`n_req = M_supp,`

`q = t R_U,`

`r = t R_phi.`

So the exact continuum-selected branch equation is

`M_supp`
`= [ xi(delta + xi) - M_mix( delta + (1 + lambda_0 R_U^2) xi ) ]`
`  / [ delta + (1 + lambda_0 R_phi^2) xi - M_mix lambda_0 (R_U - R_phi)^2 ].`

This is already a complete exact statement of the selected wall branch.

---

## 3. Exact quadratic theorem for the softening depth

The key structural simplification is that the continuum-selected branch equation is exactly quadratic in `xi`.

Rearranging gives

`xi^2 + B_cont xi + C_cont = 0,`

with

`B_cont`
`= delta`
`  - M_mix (1 + lambda_0 R_U^2)`
`  - M_supp (1 + lambda_0 R_phi^2),`

`C_cont`
`= - delta (M_mix + M_supp)`
`  + lambda_0 M_mix M_supp (R_U - R_phi)^2.`

So the physical softening depth is

`xi_phys`
`= [ - B_cont + sqrt( B_cont^2 - 4 C_cont ) ] / 2,`

where the `+` branch is selected because it reduces continuously to

`xi_phys = 0`

when `M_mix = M_supp = 0`.

This is the exact quadratic selected-branch theorem.

So the physical branch is now no longer a search over a continuous family. It is an explicit algebraic root of the continuum kernel data.

---

## 4. Exact continuum-selected normalization test

Stage 25 gave the exact rank-2 normalization law

`F_(q,r,t)(xi,delta;m)`
`= [ delta + (1 + q r) xi ]^2`
`  [ delta + (1 + r t) xi - m(q-r)(q-t) ]^2`
`  / [ (1 - xi) D_(q,r)^2 ],`

with

`D_(q,r)`
`= [ delta + xi - m q(q-r) ]^2 + [ m(q-r) + r xi ]^2.`

Substituting the continuum-selected data gives

`F_cont(xi)`
`= [ delta + (1 + lambda_0 R_U R_phi) xi ]^2`
`  [ delta + (1 + lambda_0 R_phi) xi`
`    - lambda_0 M_mix (R_U - R_phi)(R_U - 1) ]^2`
`  / [ (1 - xi) D_cont(xi)^2 ],`

with

`D_cont(xi)`
`= [ delta + xi - lambda_0 M_mix R_U (R_U - R_phi) ]^2`
`  + lambda_0 [ M_mix (R_U - R_phi) + R_phi xi ]^2.`

The actual normalization theorem gate is now simply

`R_target = F_cont( xi_phys ).`

So the whole selected-branch problem has collapsed to two exact scalar equations:

1. the quadratic branch-selection equation for `xi_phys`,
2. the normalization residual equation `R_target - F_cont(xi_phys) = 0`.

---

## 5. Exact special surfaces

### 5.1 Minimal-kernel source-tied surface

If the support-interference ratio vanishes,

`sigma_0 = 0,`

then

`R_phi = 1,`

and Stage 26 shows that the support direction is exactly source-tied.

So the **minimal** continuum kernel lands on the source-tied Stage-24/25 closure.

### 5.2 Interference-matched tracking surface

If

`g_B g_R = g_W g_S,`

equivalently

`sigma_0 = rho_0,`

then

`R_phi = R_U,`

so the support and mixed directions coincide exactly.

In that case the continuum-selected branch equation collapses to

`M_mix + M_supp = G_q(xi,delta),`

with

`q = t R_U,`

`G_q(xi,delta) = xi(delta+xi) / [ delta + (1 + q^2) xi ].`

So on the tracking surface the two rank-1 loadings merge into a single one-direction branch with exact total baseline

`M_tot = M_mix + M_supp.`

---

## 6. Exact mismatch penalty and why the generic branch is intermediate

The direction mismatch is

`Delta_R := R_U - R_phi`
`= delta_U (sigma_0 - rho_0)`
`  / [ (1+delta_U)(1+rho_0)(1+sigma_0) ].`

The quadratic branch equation shows that the genuine rank-2 penalty enters through

`lambda_0 M_mix M_supp Delta_R^2`

in `C_cont`.

So the mismatch penalty is exact, positive, and quadratic in the support/mixed direction split.

This yields the structural conclusion:

- source-tied is one exact special surface,
- tracking is another exact special surface,
- and the generic extended continuum kernel sits between them with a genuine positive rank-2 mismatch penalty.

So the physical kernel does not generically choose one of the two Stage-24 extremes. It defines a continuum-selected intermediate closure.

---

## 7. Best current theorem statement after Stage 27

The support-direction bottleneck is now resolved at the reduced-theorem level.

### What is now exact

- the actual continuum-selected support direction `R_phi`,
- the actual continuum-selected support baseline `M_supp`,
- the exact quadratic physical branch equation for `xi_phys`,
- the exact continuum-selected normalization function `F_cont`,
- the minimal-kernel source-tied surface,
- the interference-matched tracking surface,
- and the exact mismatch penalty `lambda_0 M_mix M_supp (R_U - R_phi)^2`.

### The sharp current conclusion

> The minimal Stage-20 continuum kernel lands on the source-tied closure.
> The first symmetry-allowed extended continuum kernel lands generically on an exact intermediate closure.
> Exact tracking is a special interference-match surface, not the generic outcome.

So the next theorem step is no longer to decide which abstract closure to use.
It is to evaluate the continuum-selected residual

`R_target - F_cont(xi_phys)`

on concrete kernel data from the moving-throat operator.

=== moving_throat_pde_stage028_coherent_local_tracking.md ===

# Moving-Throat PDE — Stage 28: Coherent Local D/N Support Kernel and the Exact Tracking Reduction

## Purpose

Stage 27 showed that the first extended continuum kernel generically lands on an exact **intermediate** rank-2 closure, with the source-tied and tracking laws appearing only as special surfaces.

That was the correct generic reduced statement, but it still left one obvious question open:

> what happens if the support field `phi` and the mixed field `W` really are two D/N channels of the **same local throat-support density** rather than two unrelated reduced couplings?

That is the first honest “concrete kernel” question after Stage 27.

This stage answers it.

The main result is that the first coherent local D/N kernel does **not** land on the generic intermediate closure. It lands **exactly on the tracking surface**.

So the Stage-27 continuum-selected problem collapses back to a one-parameter selected branch, now with a physically identified total loading

`M_tr = M_mix + M_supp`

and a single common direction factor

`R_tr = R_U = R_phi`.

This is the first explicit evaluation of the Stage-27 residual on a concrete moving-throat kernel family.

---

## 1. Coherent local D/N support kernel

Keep the Stage-22/26 split-`U` continuum operator on the finite throat interval `s in [0,L]`, with

- wall basis `(u_0,u_1)`,
- split internal doublet `U`,
- mixed D/N half-wave `W`,
- support/BdG D/N half-wave `phi`.

Now impose the first concrete local-kernel hypothesis:

> the mixed lane `W` and the support lane `phi` couple to the **same local wall/U support density**.

The minimal local interaction density is

`L_int^(coh)`
`= - int_0^L ds [ lambda_W W(s,t) + lambda_phi phi(s,t) ] [ eta(s,t) - gamma U(s,t) ].`

So the two D/N lanes differ only by their amplitudes `lambda_W` and `lambda_phi`; their wall/U source combination is the same.

After the same modal projection used in Stages 20, 22, and 26, this gives

`L_(eta W)   = - g_W (v.Q) Wbar,`

`L_(UW)      = + g_R (v.U) Wbar,`

`L_(eta phi) = - g_B (v.Q) Phi,`

`L_(Uphi)    = + g_S (v.U) Phi,`

with the exact amplitude pattern

`g_W = lambda_W / sqrt(mu_eta mu_W),`

`g_R = gamma lambda_W / sqrt(mu_U mu_W),`

`g_B = lambda_phi / sqrt(mu_eta mu_phi),`

`g_S = gamma lambda_phi / sqrt(mu_U mu_phi).`

This is more specific than the Stage-26 generic extension: the same local support density forces the wall/U ratio to be identical in the mixed and support lanes.

---

## 2. Exact tracking theorem

Stage 26 proved that support and mixed directions coincide iff

`g_B g_R = g_W g_S`.

For the coherent local kernel above,

`g_B g_R`
`= gamma lambda_W lambda_phi / sqrt(mu_eta mu_phi mu_U mu_W),`

`g_W g_S`
`= gamma lambda_W lambda_phi / sqrt(mu_eta mu_W mu_U mu_phi),`

so

`g_B g_R - g_W g_S = 0`

exactly.

Therefore the Stage-26 codimension-one interference-match surface is now automatic.

### Exact conclusion

> **The first coherent local D/N support kernel lands exactly on the tracking surface.**

Equivalently,

`R_phi = R_U`.

So the generic Stage-27 intermediate closure is not the first coherent local-kernel outcome. The first coherent local-kernel outcome is the exact tracking reduction.

---

## 3. Exact common interference ratio and direction factor

Stage 26 defines

`rho_0 = g_R g_U / (K_U g_W),`

`sigma_0 = g_U g_S / (K_U g_B)`.

Under the coherent local kernel,

`rho_0 = sigma_0`.

Define the common value by

`chi_0 := rho_0 = sigma_0.`

Then the exact common direction factor is

`R_tr := R_U = R_phi`
`     = [ 1 + chi_0/(1+delta_U) ] / (1 + chi_0).`

Two exact identities are useful:

`1 - R_tr`
`= chi_0 delta_U / [ (1 + chi_0)(1 + delta_U) ],`

`R_tr - 1/(1+delta_U)`
`= delta_U / [ (1 + chi_0)(1 + delta_U) ].`

So for the constructive branch `chi_0 > 0`, `delta_U > 0`,

`1/(1+delta_U) < R_tr < 1.`

That gives the exact physical range of the tracking factor on this branch.

---

## 4. Exact total loading on the coherent branch

Stage 22 and Stage 26 give the mixed and support baselines

`M_mix`
`= 8 Z_W (1 + chi_0)^2`
`  / [ pi^2 (1 - eps_eta) (1 - eps_W^(split)) ],`

`M_supp`
`= 8 Z_phi (1 + chi_0)^2`
`  / [ pi^2 (1 - eps_eta) (1 - eps_phi^(split)) ].`

Because the branch is exactly tracking, the physical selected-mode problem depends only on the total baseline

`M_tr := M_mix + M_supp,`

so

`M_tr`
`= 8 (1 + chi_0)^2 / [ pi^2 (1 - eps_eta) ]`
`  * [ Z_W / (1 - eps_W^(split)) + Z_phi / (1 - eps_phi^(split)) ].`

So the first concrete local-kernel evaluation of the Stage-27 residual is already much simpler than the generic intermediate closure:

- one common direction factor `R_tr`,
- one total baseline `M_tr`.

---

## 5. Exact collapse of the Stage-27 quadratic branch equation

Stage 27 gave the continuum-selected quadratic branch equation

`xi^2 + B_cont xi + C_cont = 0,`

with mismatch penalty proportional to

`lambda_0 M_mix M_supp (R_U - R_phi)^2.`

On the coherent local-kernel branch,

`R_U - R_phi = 0,`

so the mismatch penalty vanishes identically. The physical branch equation reduces exactly to

`xi^2 + B_tr xi + C_tr = 0,`

with

`B_tr = delta - M_tr ( 1 + lambda_0 R_tr^2 ),`

`C_tr = - delta M_tr.`

Equivalently, the selected branch obeys the exact one-parameter law

`M_tr = G_tr(xi,delta;R_tr)`,

where

`G_tr(xi,delta;R)`
`= 9 xi (xi + delta) / [ 9 delta + (9 + 2 R^2) xi ].`

This is exactly the Stage-23 tracking/source-loading branch, now derived from a concrete local kernel rather than postulated as a special surface.

---

## 6. Exact normalization law on the coherent branch

Because the coherent local kernel lands on the tracking surface, the Stage-27 normalization residual also collapses to the Stage-23 tracking form:

`R_target = F_tr(xi,delta;R_tr)`,

with

`F_tr(xi,delta;R)`
`= [ 9 delta + (9 + 2 R^2) xi ]^2 [ 9 delta + (9 + 2 R) xi ]^2`
`  / [ 81 (1 - xi) ( 9 delta^2 + 18 delta xi + (9 + 2 R^2) xi^2 )^2 ].`

So the full Stage-27 continuum-selected residual now becomes

`R_target - F_tr( xi_phys, delta; R_tr ),`

where `xi_phys` is the physical root of the reduced quadratic above.

This is the first exact “concrete-kernel” form of the normalization test.

---

## 7. Best current theorem statement after Stage 28

### What is now exact

- the first coherent local D/N support kernel,
- the automatic tracking identity `g_B g_R = g_W g_S`,
- the common interference ratio `chi_0`,
- the exact tracking factor `R_tr`,
- the exact total baseline `M_tr`,
- the exact collapse of the Stage-27 quadratic branch equation to the one-parameter tracking law,
- and the exact collapse of the normalization residual to the tracking function `F_tr`.

### What this means physically

The generic Stage-27 intermediate closure is still the correct reduced theorem for the first unrestricted continuum extension. But the first **coherent local D/N kernel** is more special than that: it lands exactly on the tracking surface.

So the next theorem step is no longer to resolve the rank-2 closure ambiguity. That is now done for the first concrete local kernel. The next step is to evaluate how this physical tracking branch compares with the old flat branch and whether the split-`U` deformation helps or hurts the normalization test.

=== moving_throat_pde_stage029_tracking_branch_bounds.md ===

# Moving-Throat PDE — Stage 29: Exact Tracking-Branch Bounds, Monotonicity in `R_tr`, and the Residual Comparison Theorem

## Purpose

Stage 28 reduced the first coherent local D/N kernel to an exact one-parameter tracking branch:

`M_tr = G_tr(xi,delta;R_tr),`

`R_target = F_tr(xi,delta;R_tr).`

So the rank-2 closure ambiguity is gone on that branch.

The next honest question is not “what closure do we use?” It is:

> does the physical split-`U` deformation help or hurt the normalization test relative to the old flat branch?

This stage answers that exactly.

The main result is that the coherent local tracking branch is **ordered** by the tracking factor `R_tr`.

At fixed `(xi,delta)`:

- the required total loading `G_tr` decreases with `R_tr`,
- the normalized selected-branch response `F_tr` increases with `R_tr`.

Because the constructive split-`U` branch satisfies `R_tr < 1`, the physical local D/N branch therefore

- requires **more** total loading than the flat branch,
- but delivers **less** normalized response than the flat branch.

So the split-`U` deformation makes the Stage-18/19 normalization target harder, not easier.

---

## 1. Tracking-branch functions

On the coherent local D/N branch from Stage 28,

`G_tr(xi,delta;R)`
`= 9 xi (xi + delta) / [ 9 delta + (9 + 2 R^2) xi ],`

`F_tr(xi,delta;R)`
`= [ 9 delta + (9 + 2 R^2) xi ]^2 [ 9 delta + (9 + 2 R) xi ]^2`
`  / [ 81 (1 - xi) ( 9 delta^2 + 18 delta xi + (9 + 2 R^2) xi^2 )^2 ].`

The old flat branch is recovered at

`R = 1,`

so

`G_flat(xi,delta) = G_tr(xi,delta;1) = 9 xi (xi + delta) / (9 delta + 11 xi),`

`F_flat(xi,delta) = F_tr(xi,delta;1)`
`                 = (9 delta + 11 xi)^4`
`                   / [ 81 (1 - xi) (9 delta^2 + 18 delta xi + 11 xi^2)^2 ].`

---

## 2. Exact monotonicity in the tracking factor `R`

Differentiate the branch functions at fixed `(xi,delta)`.

### 2.1 Loading monotonicity

The exact derivative is

`dG_tr/dR`
`= - 36 R xi^2 (delta + xi) / [ 2 R^2 xi + 9 delta + 9 xi ]^2.`

So for the physical branch `R > 0`,

`dG_tr/dR < 0.`

Thus lowering `R` increases the total loading required to realize the same softening depth `xi`.

### 2.2 Normalization monotonicity

The exact derivative is

`dF_tr/dR`
`= 4 xi (2 R xi + 9 delta + 9 xi) (2 R^2 xi + 9 delta + 9 xi) P_R`
`  / [ 81 (1 - xi) ( 2 R^2 xi^2 + 9 delta^2 + 18 delta xi + 9 xi^2 )^3 ],`

with the positive polynomial

`P_R`
`= 4 R^4 xi^3`
`  + 54 R^2 delta^2 xi + 90 R^2 delta xi^2 + 36 R^2 xi^3`
`  + 162 R delta^3 + 324 R delta^2 xi + 162 R delta xi^2`
`  + 81 delta^3 + 243 delta^2 xi + 243 delta xi^2 + 81 xi^3.`

Every coefficient in `P_R` is positive, so for the physical branch `0 < xi < 1`, `delta > 0`, `R > 0`,

`dF_tr/dR > 0.`

Thus lowering `R` decreases the normalized selected-branch response at fixed `(xi,delta)`.

---

## 3. Exact comparison with the flat branch

Since the constructive split-`U` branch has `R_tr < 1`, the derivative theorems above already imply

`G_tr(xi,delta;R_tr) > G_flat(xi,delta),`

`F_tr(xi,delta;R_tr) < F_flat(xi,delta).`

But the comparison can be written in exact closed form.

### 3.1 Exact loading excess over the flat branch

`G_tr - G_flat`
`= 18 xi^2 (1 - R^2) (delta + xi)`
`  / [ (9 delta + 11 xi) (2 R^2 xi + 9 delta + 9 xi) ].`

So for `0 < R < 1`,

`G_tr - G_flat > 0.`

### 3.2 Exact normalization deficit relative to the flat branch

`F_flat - F_tr`
`= 4 xi (1 - R) P_1 P_2`
`  / [ 81 (1 - xi) (9 delta^2 + 18 delta xi + 11 xi^2)^2`
`      (2 R^2 xi^2 + 9 delta^2 + 18 delta xi + 9 xi^2)^2 ],`

with

`P_1`
`= 18 R^2 delta^2 xi + 36 R^2 delta xi^2 + 22 R^2 xi^3`
`  + 81 R delta^3 + 180 R delta^2 xi + 99 R delta xi^2`
`  + 162 delta^3 + 423 delta^2 xi + 360 delta xi^2 + 99 xi^3,`

`P_2`
`= 18 R^3 delta^2 xi^2 + 36 R^3 delta xi^3 + 22 R^3 xi^4`
`  + 81 R^2 delta^3 xi + 324 R^2 delta^2 xi^2 + 459 R^2 delta xi^3 + 220 R^2 xi^4`
`  + 81 R delta^3 xi + 243 R delta^2 xi^2 + 261 R delta xi^3 + 99 R xi^4`
`  + 729 delta^4 + 3078 delta^3 xi + 4959 delta^2 xi^2 + 3600 delta xi^3 + 990 xi^4.`

Every coefficient in `P_1` and `P_2` is positive, so for `0 < R < 1`,

`F_flat - F_tr > 0.`

So the split-`U` tracking branch is **strictly below** the flat-branch normalization function at fixed `(xi,delta)`.

---

## 4. Exact endpoint formulas and bounds

The formal strong-split endpoint is `R -> 0`.

At that endpoint the branch functions simplify exactly to

`G_tr(xi,delta;0) = xi,`

`F_tr(xi,delta;0) = 1 / (1 - xi).`

So for `0 <= R <= 1` the exact bounds are

`G_flat(xi,delta) <= G_tr(xi,delta;R) <= xi,`

`1/(1 - xi) <= F_tr(xi,delta;R) <= F_flat(xi,delta).`

These are absolute tracking-branch bounds. The actual constructive branch has the tighter Stage-28 range

`1/(1+delta_U) < R_tr < 1,`

so the physical kernel sits strictly between the strong-split endpoint and the flat branch.

---

## 5. Exact residual comparison theorem

Define the tracking-branch normalization residual

`E_tr(xi) := R_target - F_tr(xi,delta;R_tr),`

and the flat-branch residual

`E_flat(xi) := R_target - F_flat(xi,delta).`

Then at fixed `(xi,delta)`,

`E_tr - E_flat = F_flat - F_tr > 0`

on the constructive split-`U` branch.

So the local split-`U` deformation worsens the normalization residual relative to the old flat branch.

Likewise the exact loading excess theorem says that the local split-`U` branch requires a larger total baseline to reach the same `xi`.

This is the sharpest exact comparison yet between the abstract Stage-18/19 branch and the first coherent concrete local kernel.

---

## 6. Best current theorem statement after Stage 29

### What is now exact

- the tracking-branch functions `G_tr` and `F_tr`,
- their exact derivatives with respect to `R`,
- the exact loading excess `G_tr - G_flat`,
- the exact normalization deficit `F_flat - F_tr`,
- the strong-split endpoint formulas,
- and the exact residual ordering `E_tr > E_flat` at fixed `(xi,delta)` on the constructive split-`U` branch.

### What this means physically

The first coherent local D/N kernel does not rescue the Stage-18/19 normalization target by hidden split-`U` assistance. It does the opposite.

The constructive split-`U` deformation

- lowers the tracking factor from `1` to `R_tr < 1`,
- increases the total loading required to reach a given softening depth,
- and decreases the normalized quadrupole response at the same point.

So the remaining theorem gap is now even sharper:

> the completed moving-throat PDE must supply kernel data strong enough to overcome this exact tracking-branch deficit, not merely to match the old flat-branch heuristic.

=== moving_throat_pde_stage030_coherent_kernel_map.md ===

# Moving-Throat PDE — Stage 30: Exact Coherent-Kernel Dimensionless Map and the Support-Enhancement Factor

## Purpose

Stage 29 left the reduced theorem gap in its sharpest form so far:

> the first coherent local D/N kernel lands on the physical tracking branch, but the constructive split-`U` deformation lowers `R_tr` below `1` and therefore worsens the normalization target relative to the old flat branch.

That was already useful, but the actual kernel-strength question was still spread over too many dimensional couplings.

The next honest move is therefore to compress the coherent local kernel to the smallest exact **dimensionless placement map** and identify how the support lane enters the branch.

This stage does that.

The main result is that the first coherent local D/N kernel collapses to a very small parameter set:

- one wall/internal dressing ratio `eps_eta`,
- one split-`U` axial ratio `delta_U`,
- one common interference ratio `chi_0`,
- one mixed blocking ratio `eps_W`,
- one wall-to-mixed overlap ratio `Z_W`,
- one support-to-mixed amplitude ratio `zeta`,
- and one radiative demand scale `Lambda`.

In those variables the coherent branch is fully described by

`R_tr = [ 1 + chi_0/(1+delta_U) ] / (1 + chi_0),`

`eps = eps_W^(split) = eps_W [ 1 - (2/11) delta_U/(1+delta_U) ],`

`delta = [ delta_0 + eps_eta delta_U/(1+delta_U) ] / (1 - eps_eta),`

`M_mix = 8 Z_W (1 + chi_0)^2 / [ pi^2 (1 - eps_eta) (1 - eps) ],`

`M_supp = 8 zeta Z_W (1 + chi_0)^2 / [ pi^2 (1 - eps_eta) (1 - zeta eps) ],`

`R_target = Lambda (1 - eps_eta) (1 - eps)^2 / [ Z_W (1 + chi_0)^2 ].`

So the entire support lane enters through one exact factor

`S(zeta;eps) := M_tr / M_mix = 1 + zeta(1-eps)/(1-zeta eps),`

with

`M_tr = M_mix S(zeta;eps).`

That turns the support problem into a one-parameter enhancement problem rather than a diffuse multidimensional closure question.

---

## 1. Coherent local interaction data

On the coherent local D/N branch from Stage 28 the mixed and support lanes couple through the same local source density,

`L_int^(coh)`
`= - int_0^L ds [ lambda_W W + lambda_phi phi ] [ eta - gamma U ].`

So the continuum couplings obey

`c_(etaW) = lambda_W,`

`c_(UW)   = gamma lambda_W,`

`c_(etaphi) = lambda_phi,`

`c_(Uphi)   = gamma lambda_phi.`

The wall/internal coupling `c_(etaU)` remains independent.

The split-`U` stiffness is

`K_(U1) = K_U (1 + delta_U),`

with

`delta_U = pi^2 T_U / (L^2 K_U).`

The effective wall and D/N support stiffnesses are

`K_eta^(eff) = K_eta + 6 T_Omega,`

`K_W^(eff)   = K_W   + pi^2 T_W   / (4 L^2),`

`K_phi^(eff) = K_phi + pi^2 T_phi / (4 L^2).`

As before,

`sigma = kappa_0^2 + kappa_1^2 = 88 / (9 pi^2).`

---

## 2. Exact coherent dimensionless ratios

The natural coherent-kernel dimensionless ratios are

`eps_eta := c_(etaU)^2 / ( K_U K_eta^(eff) ),`

`eps_W   := gamma^2 lambda_W^2 sigma / ( K_U K_W^(eff) ),`

`Z_W     := lambda_W^2 / ( K_eta^(eff) K_W^(eff) ),`

`delta_0 := pi^2 T_w / ( L^2 K_eta^(eff) ),`

`chi_0   := gamma c_(etaU) / K_U,`

`zeta    := lambda_phi^2 K_W^(eff) / ( lambda_W^2 K_phi^(eff) ),`

`Lambda  := 27 pi^2 G c_s^5 K_W^(eff) / (20 a^5 c^5 mu_W).`

Two exact coherent identities then follow immediately:

`rho_0 = sigma_0 = chi_0,`

`eps_phi = zeta eps_W,`

`Z_phi   = zeta Z_W.`

So the support lane is not independent of the mixed lane at the level of dimensionless couplings. It differs only by the one ratio `zeta`.

---

## 3. Exact split data on the coherent branch

The coherent tracking factor is

`R_tr = [ 1 + chi_0/(1+delta_U) ] / (1 + chi_0),`

with the same exact range as Stage 28:

`1/(1+delta_U) < R_tr < 1`

on the constructive branch `chi_0 > 0`, `delta_U > 0`.

The split mixed and support blocking ratios are

`eps = eps_W^(split)`
`    = eps_W [ 1 - (2/11) delta_U/(1+delta_U) ],`

`eps_phi^(split) = zeta eps.`

The geometry anisotropy is still

`delta = [ delta_0 + eps_eta delta_U/(1+delta_U) ] / (1 - eps_eta).`

So once the coherent kernel is imposed, the branch direction, the blocking, and the anisotropy are all explicit functions of the same small dimensionless set.

---

## 4. Exact mixed/support baselines and support-enhancement factor

The mixed baseline becomes

`M_mix = 8 Z_W (1 + chi_0)^2 / [ pi^2 (1 - eps_eta) (1 - eps) ].`

Using `Z_phi = zeta Z_W` and `eps_phi^(split) = zeta eps`, the support baseline is

`M_supp = 8 zeta Z_W (1 + chi_0)^2 / [ pi^2 (1 - eps_eta) (1 - zeta eps) ].`

Therefore the exact total baseline is

`M_tr = M_mix + M_supp`
`     = M_mix S(zeta;eps),`

with the exact support-enhancement factor

`S(zeta;eps)`
`:= 1 + zeta(1-eps)/(1-zeta eps).`

This is the cleanest structural simplification of the stage.

The support lane does not change the radiative demand ratio directly. It only enhances the total baseline by the one factor `S`.

### Exact properties of `S`

`S(0;eps) = 1,`

`dS/dzeta = (1-eps)/(1-zeta eps)^2 > 0`

for `0 < eps < 1` and `0 <= zeta < 1/eps`.

So the coherent support lane is a strictly monotone baseline enhancer.

---

## 5. Exact target ratio and product law on the coherent branch

The normalization-demand ratio remains purely mixed/outgoing:

`R_target = Lambda (1 - eps_eta) (1 - eps)^2 / [ Z_W (1 + chi_0)^2 ].`

So `R_target` is exactly independent of `zeta`.

Multiplying by the total coherent baseline gives

`R_target M_tr`
`= 8 Lambda (1 - eps) / pi^2 * S(zeta;eps).`

Equivalently,

`R_target M_tr`
`= 8 Lambda (1 - eps) / pi^2`
`  * [ 1 + zeta(1-eps)/(1-zeta eps) ].`

So the mixed lane still sets the bare product scale, while the support lane simply multiplies that scale by the one enhancement factor `S`.

This is the exact coherent-kernel replacement of the Stage-21 product law.

---

## 6. Best current theorem statement after Stage 30

The first coherent local D/N kernel is now compressed to an exact dimensionless placement map.

### What is exact now

- the coherent branch depends on one common interference ratio `chi_0`,
- the support lane differs from the mixed lane by only one ratio `zeta`,
- the exact tracking factor is `R_tr(chi_0,delta_U)`,
- the exact mixed blocking is `eps = eps_W^(split)`,
- the exact total baseline is `M_tr = M_mix S(zeta;eps)`,
- the target ratio `R_target` is independent of `zeta`,
- and the support lane is a strictly monotone enhancer of the total baseline.

### What this means physically

The reduced theorem gap is no longer “how do we represent the support lane?”
It is now much sharper:

> for the first coherent local D/N kernel, the support lane enters only by increasing the baseline through one exact enhancement factor, while the retarded demand ratio is still fixed by the mixed/outgoing lane.

So the next honest question is no longer about closure choice.
It is:

> **how much coherent support enhancement is needed to overcome the exact tracking-branch deficit before the stable branch softens out?**

=== moving_throat_pde_stage031_support_compensation_theorem.md ===

# Moving-Throat PDE — Stage 31: Exact Support-Compensation Theorem on the Physical Tracking Branch

## Purpose

Stage 30 compressed the coherent local D/N kernel to one exact support-enhancement factor,

`M_tr = M_mix S(zeta;eps),`

with

`S(zeta;eps) = 1 + zeta(1-eps)/(1-zeta eps).`

So the reduced problem is no longer “what support closure do we use?”
It is now a direct kernel-strength question:

> can the coherent support lane overcome the exact tracking-branch normalization deficit before the selected branch reaches its softening limit?

This stage answers that exactly.

The main result is that there is **no reduced-level support no-go** on the physical tracking branch.

For any finite target ratio `R_target > 1`, any stable geometry `(delta,R_tr)`, and any mixed-only baseline below the tracking critical load, there is a unique coherent support ratio `zeta_req` that hits the target before softening.

So the remaining theorem gap is not a reduced closure obstruction. It is whether the completed moving-throat PDE supplies a physical `zeta` large enough to reach that exact required value.

---

## 1. Universal tracking branch and its critical load

On the coherent branch the selected mode obeys the exact tracking law

`M_tr = G_tr(xi,delta;R_tr),`

with

`G_tr(xi,delta;R)`
`= 9 xi (xi + delta) / [ 9 delta + (9 + 2 R^2) xi ].`

The stable selected branch lives on

`0 < xi < 1.`

The exact derivative is

`dG_tr/dxi`
`= 9 [ 2 R^2 xi^2 + 9 delta^2 + 18 delta xi + 9 xi^2 ]`
`  / [ 2 R^2 xi + 9 delta + 9 xi ]^2 > 0.`

So `G_tr` is strictly increasing along the stable branch.

Its exact softening-limit value is therefore

`M_crit(delta,R)`
`:= G_tr(1,delta;R)`
`= 9 (1 + delta) / [ 9 delta + 9 + 2 R^2 ].`

Equivalently,

`M_crit - G_tr(xi,delta;R)`
`= 9 (1-xi) [ 2 R^2 xi + 9 delta^2 + 9 delta xi + 9 delta + 9 xi ]`
`  / [ (2 R^2 + 9 delta + 9) (2 R^2 xi + 9 delta + 9 xi) ] > 0`

for `0 < xi < 1`.

So every stable selected-branch point lies strictly below the critical load.

---

## 2. Exact support-enhancement factor and its inverse

From Stage 30,

`M_tr = M_mix S(zeta;eps),`

`S(zeta;eps) = 1 + zeta(1-eps)/(1-zeta eps),`

with the stability window

`0 < eps < 1,`

`0 <= zeta < 1/eps.`

The exact derivative is

`dS/dzeta = (1-eps)/(1-zeta eps)^2 > 0,`

so `S` is strictly increasing.

The endpoint values are

`S(0;eps) = 1,`

`lim_(zeta -> 1/eps^-) S(zeta;eps) = +infinity.`

So the support lane can realize any finite enhancement factor above `1` while staying below the blocking pole.

The inverse map is exact. If `S_req > 1`, then the unique support ratio producing that enhancement is

`zeta_req = (S_req - 1) / [ 1 + eps (S_req - 2) ].`

It also obeys the exact stability margin

`1/eps - zeta_req`
`= (1-eps) / [ eps ( 1 + eps (S_req - 2) ) ] > 0.`

So every finite required enhancement lies strictly inside the stability window.

---

## 3. Existence of a stable-side target point

The exact tracking normalization function is

`F_tr(xi,delta;R)`
`= [ 9 delta + (9 + 2 R^2) xi ]^2 [ 9 delta + (9 + 2 R) xi ]^2`
`  / [ 81 (1 - xi) ( 9 delta^2 + 18 delta xi + (9 + 2 R^2) xi^2 )^2 ].`

Two endpoint facts are exact:

`F_tr(0,delta;R) = 1,`

`lim_(xi -> 1^-) F_tr(xi,delta;R) = +infinity.`

Therefore, by continuity, for every finite target ratio

`R_target > 1`

there exists at least one stable-side root

`xi_req in (0,1)`

such that

`F_tr(xi_req,delta;R_tr) = R_target.`

Define the corresponding required total load by

`M_req := G_tr(xi_req,delta;R_tr).`

Because `0 < xi_req < 1` and `G_tr` is strictly increasing,

`0 < M_req < M_crit(delta,R_tr).`

So the normalization target, if finite, is always associated with a finite stable-side required load.

---

## 4. Exact support-compensation theorem

Suppose the mixed-only coherent branch is stable,

`0 < M_mix < M_crit(delta,R_tr).`

There are then two cases.

### Case A — mixed-only branch already strong enough

If

`M_mix >= M_req,`

then the target is already reachable with

`zeta_req = 0.`

### Case B — support enhancement is needed

If

`M_mix < M_req,`

define the exact required enhancement factor

`S_req := M_req / M_mix > 1.`

Then the unique coherent support ratio that hits the target is

`zeta_req = (S_req - 1) / [ 1 + eps (S_req - 2) ].`

And because `M_req < M_crit`, one has

`S_req < S_crit := M_crit / M_mix,`

so the support ratio lies strictly below the critical support ratio

`zeta_crit = (S_crit - 1) / [ 1 + eps (S_crit - 2) ].`

Indeed,

`zeta_crit - zeta_req`
`= (S_crit - S_req)(1-eps)`
`  / [ (1 + eps(S_crit - 2))(1 + eps(S_req - 2)) ] > 0.`

So the target is reached **before** the selected branch softens out.

---

## 5. Exact monotone wall-softening response to support enhancement

Combining Stages 30 and 31 gives the exact implicit selected-branch law

`M_mix S(zeta;eps) = G_tr(xi_phys,delta;R_tr).`

Differentiating implicitly gives

`dxi_phys/dzeta`
`= M_mix (dS/dzeta) / (dG_tr/dxi)|_(xi_phys)`
`> 0.`

So coherent support enhancement always drives the selected branch to larger softening depth.

This is the exact reduced statement behind the compensation theorem:

- support does not alter `R_target`,
- it only increases the available baseline,
- and that increase moves the physical branch monotonically deeper into the tracking family.

---

## 6. Best current theorem statement after Stage 31

The coherent local D/N branch now has an exact support-feasibility theorem.

### What is exact now

- the tracking branch has a finite critical load `M_crit(delta,R_tr)`,
- the coherent support-enhancement factor `S(zeta;eps)` is strictly increasing and invertible,
- every finite target ratio `R_target > 1` corresponds to at least one stable-side branch point `xi_req in (0,1)`,
- the corresponding required load `M_req` is always strictly below `M_crit`,
- and whenever `M_mix < M_req`, there is a unique coherent support ratio `zeta_req < zeta_crit < 1/eps` that reaches the target before softening.

### What this means physically

The exact tracking-branch deficit from Stages 28–29 is real, but it is **not** a reduced-level no-go.

The first coherent local D/N support lane can compensate that deficit exactly.
The remaining question is no longer whether compensation is possible in principle.
It is narrower:

> **does the actual moving-throat PDE produce a physical support ratio `zeta` large enough to reach the exact required value `zeta_req` on the real branch?**

=== moving_throat_pde_stage032_dn_overlap_zeta.md ===

# Moving-Throat PDE — Stage 32: Explicit D/N Overlap Extraction of the Physical Coherent Support Ratio `zeta`

## Purpose

Stage 31 reduced the whole support-compensation problem to one exact inequality:

`zeta_phys >= zeta_req`.

But `zeta` was still a reduced parameter. The next honest move is therefore to derive it from a more explicit finite-throat support operator instead of leaving it phenomenological.

This stage does that for the first concrete local D/N family.

The main result is that, once the mixed lane `W` and the support lane `phi_n` are both generated by the **same local throat source density** on the finite throat interval, the physical coherent support ratio is no longer free. It is exactly

`zeta_n^(phys)`
`= (K_W^(eff) / K_(phi,n)^(eff)) * (I_n / I_0)^2,`

where

`I_n := int_0^L ds sigma(s) chi_n(s)`

is the source-to-mode overlap on the D/N ladder

`chi_n(s) = sqrt(2/L) sin[(n+1/2) pi s / L].`

For the first uniform local source density `sigma(s)=1`, the overlap hierarchy is exact:

`I_n = sqrt(2L) / [ (n+1/2) pi ],`

so

`I_n / I_0 = 1 / (2n+1).`

Therefore the physical coherent support ratio becomes

`zeta_n^(phys)`
`= (K_W^(eff) / K_(phi,n)^(eff)) / (2n+1)^2.`

This is the first direct microscopic comparison of the support lane against the exact Stage-31 support threshold.

---

## 1. Explicit finite-throat coherent source operator

Take the finite throat interval `s in [0,L]` and the exact D/N axial family

`chi_n(s) = sqrt(2/L) sin[(n+1/2) pi s / L],`

with

`k_n = (n+1/2) pi / L.`

Now let the mixed lane `W` and the support lane `phi_n` couple to the same local source density `sigma(s)` through

`L_int^(coh)`
`= - lambda_* int_0^L ds sigma(s) [ W(s,t) + phi(s,t) ] Sigma_src(s,t),`

where `Sigma_src` is the same wall/internal source combination already used in the coherent branch.

Expanding into the D/N basis,

`W(s,t)   = w_0(t) chi_0(s),`

`phi(s,t) = sum_n ph_n(t) chi_n(s),`

the reduced source amplitudes are

`lambda_W = lambda_* I_0,`

`lambda_(phi,n) = lambda_* I_n,`

with the exact overlap integrals

`I_n := int_0^L ds sigma(s) chi_n(s).`

So the Stage-30 definition of the coherent support ratio,

`zeta = lambda_phi^2 K_W^(eff) / [ lambda_W^2 K_phi^(eff) ],`

becomes the exact microscopic formula

`zeta_n^(phys)`
`= (K_W^(eff) / K_(phi,n)^(eff)) * (I_n / I_0)^2.`

This is the correct finite-throat replacement of the earlier phenomenological `zeta`.

---

## 2. Uniform local source density and the exact D/N overlap hierarchy

For the first explicit local kernel, choose the uniform source density

`sigma(s)=1.`

Then

`I_n = int_0^L ds sqrt(2/L) sin[(n+1/2) pi s / L]`
`    = sqrt(2L) / [ (n+1/2) pi ].`

Therefore

`I_n / I_0 = 1 / (2n+1).`

So the physical support ratio is now completely explicit:

`zeta_n^(phys)`
`= (K_W^(eff) / K_(phi,n)^(eff)) / (2n+1)^2.`

This is already enough to compare directly against the exact Stage-31 threshold `zeta_req`.

---

## 3. Same-operator twin-lane specialization

The cleanest concrete branch is when `W` and `phi_n` are not just sourced coherently, but belong to the **same D/N operator family** with common base stiffness `K_X` and tension `T_X`.

Then

`K_W^(eff) = K_X + pi^2 T_X / (4 L^2),`

`K_(phi,n)^(eff) = K_X + (n+1/2)^2 pi^2 T_X / L^2`
`                = K_W^(eff) + pi^2 T_X n(n+1) / L^2.`

So the exact physical support ratio becomes

`zeta_n^(twin)`
`= 1 / [ (2n+1)^2 ( 1 + x n(n+1) ) ],`

with the positive stiffness parameter

`x := pi^2 T_X / [ L^2 K_W^(eff) ].`

This yields three immediate facts.

### Lowest twin half-wave

For `n=0`,

`zeta_0^(twin) = 1.`

So the symmetric lowest D/N twin lane is an exact `zeta=1` branch.

### Higher D/N harmonics are suppressed

For every `n>=1`,

`zeta_n^(twin) < 1 / (2n+1)^2 <= 1/9.`

So higher support harmonics are automatically much weaker support enhancers than the lowest twin lane.

### Microscopic branch ordering

Because `zeta_n^(twin)` is strictly decreasing in both `n` and `x`, the physical local support hierarchy is:

- the lowest twin half-wave is the strongest support branch,
- every higher D/N harmonic is penalized first by overlap,
- and then penalized again by its larger axial stiffness.

This is the first explicit microscopic reason the lowest coherent twin lane should dominate the support-compensation problem.

---

## 4. Direct comparison with the Stage-31 support threshold

Stage 31 says the support lane succeeds exactly when

`zeta_phys >= zeta_req`.

Stage 32 now turns that into an explicit microscopic inequality.

### General coherent D/N branch

For a source density `sigma(s)`, the exact success condition is

`(K_W^(eff) / K_(phi,n)^(eff)) * (I_n / I_0)^2 >= zeta_req.`

### Uniform local source

On the first explicit local kernel `sigma=1`, this becomes

`K_(phi,n)^(eff) <= K_W^(eff) / [ (2n+1)^2 zeta_req ].`

### Same-operator twin family

For the same-operator branch,

`1 / [ (2n+1)^2 (1 + x n(n+1)) ] >= zeta_req.`

For `n>=1`, this is equivalent to

`x <= x_max(n; zeta_req)`
`   := [ 1 / ( (2n+1)^2 zeta_req ) - 1 ] / [ n(n+1) ].`

So the comparison with Stage 31 is now exact.

If `x_max < 0`, that harmonic is impossible even before stiffness details are examined. Equivalently,

`zeta_req > 1 / (2n+1)^2`

already rules out the `n`th D/N support harmonic on the first coherent local branch.

---

## 5. Best current theorem statement after Stage 32

The support-threshold problem is no longer phenomenological at the level of the first explicit finite-throat kernel.

### What is now exact

- the microscopic coherent support ratio is
  
  `zeta_n^(phys) = (K_W^(eff)/K_(phi,n)^(eff)) (I_n/I_0)^2`,
- for the first uniform local source, `I_n/I_0 = 1/(2n+1)`,
- for the same-operator twin family,
  
  `zeta_n^(twin) = 1 / [ (2n+1)^2 (1 + x n(n+1)) ]`,
- the lowest twin branch has `zeta_0^(twin)=1`,
- and higher D/N harmonics are strongly suppressed.

### What this means physically

The remaining theorem gap is now much narrower.

It is no longer “what is `zeta`?”
It is:

> **which explicit support harmonic is physically realized by the moving-throat PDE, and does its exact microscopic `zeta_n^(phys)` beat the exact Stage-31 threshold `zeta_req`?**

That question is now concrete enough to answer branch by branch.

=== moving_throat_pde_stage033_zeta_threshold_comparison.md ===

# Moving-Throat PDE — Stage 33: Exact Comparison of the Physical Coherent `zeta` Against `zeta_req`

## Purpose

Stage 31 proved that the support-compensation problem is exactly

`zeta_phys >= zeta_req`,

while Stage 32 derived the physical coherent D/N support ratio on the first explicit finite-throat kernel,

`zeta_n^(phys)`
`= (K_W^(eff) / K_(phi,n)^(eff)) / (2n+1)^2,`

and, on the same-operator twin family,

`zeta_n^(twin)`
`= 1 / [ (2n+1)^2 ( 1 + x n(n+1) ) ].`

So the comparison that was still abstract at Stage 31 can now be carried out exactly.

The main results are:

1. the symmetric lowest twin lane `n=0` has
   
   `zeta_0^(twin)=1`,
   
   so its support-enhancement factor is exactly
   
   `S_0 = 2`,
   
   independent of `eps`;
2. therefore the lowest symmetric twin branch succeeds **iff**
   
   `S_req <= 2`,
   
   equivalently
   
   `zeta_req <= 1`;
3. every higher D/N support harmonic is bounded by
   
   `zeta_n^(twin) < 1/(2n+1)^2`,
   
   so if
   
   `zeta_req > 1/(2n+1)^2`,
   
   that harmonic is ruled out immediately;
4. when it is not ruled out immediately, the exact stiffness threshold is
   
   `x <= x_max(n; zeta_req)`
   `   = [ 1 / ( (2n+1)^2 zeta_req ) - 1 ] / [ n(n+1) ].`

So the exact support-threshold comparison is now sharp enough to distinguish the physically viable twin lane from the strongly suppressed higher-harmonic branches.

---

## 1. Exact support-threshold inequality in Stage-31 variables

Stage 31 defines the required enhancement factor

`S_req := M_req / M_mix`

and the exact required support ratio

`zeta_req = (S_req - 1) / [ 1 + eps (S_req - 2) ],`

with `0 < eps < 1` and `S_req > 1` whenever support is really needed.

Because the coherent support-enhancement factor is

`S(zeta;eps) = 1 + zeta (1-eps) / (1 - eps zeta),`

the support branch succeeds exactly when

`zeta_phys >= zeta_req.`

This is the comparison formula every explicit kernel has to meet.

---

## 2. Lowest symmetric twin lane: exact doubling theorem

For the same-operator twin branch, Stage 32 gave

`zeta_0^(twin) = 1.`

Substituting into the exact support-enhancement factor gives

`S_0 = S(1;eps)`
`    = 1 + (1-eps)/(1-eps)`
`    = 2.`

So the symmetric lowest twin lane doubles the mixed baseline exactly, independent of `eps`.

This immediately yields the exact comparison theorem.

### Exact equivalence

`zeta_req <= 1`
`<=> S_req <= 2.`

So the lowest symmetric twin lane succeeds exactly when the Stage-31 required enhancement does not exceed a factor of two.

This is the cleanest physical interpretation yet of the coherent support problem.

The first explicit finite-throat twin lane is not just “supportive.” It is a universal **doubling** branch.

---

## 3. Higher D/N support harmonics: exact impossibility bound

For `n>=1`, Stage 32 gives

`zeta_n^(twin)`
`= 1 / [ (2n+1)^2 ( 1 + x n(n+1) ) ],`

with `x>0`.

Therefore

`zeta_n^(twin) < 1 / (2n+1)^2.`

So an exact necessary condition for the `n`th support harmonic to succeed is

`zeta_req <= 1 / (2n+1)^2.`

If this fails, the `n`th D/N support harmonic is ruled out before any finer stiffness modeling is done.

This gives the first explicit microscopic no-go tower:

- `n=1` impossible if `zeta_req > 1/9`,
- `n=2` impossible if `zeta_req > 1/25`,
- `n=3` impossible if `zeta_req > 1/49`,
- and so on.

So the physical support problem is already strongly biased toward the lowest coherent twin lane.

---

## 4. Exact stiffness threshold when a higher harmonic is not yet ruled out

When

`zeta_req <= 1 / (2n+1)^2`,

one can solve the exact inequality

`zeta_n^(twin) >= zeta_req`

for the stiffness parameter `x`.

For `n>=1`, the result is

`x <= x_max(n; zeta_req)`
`   = [ 1 / ( (2n+1)^2 zeta_req ) - 1 ] / [ n(n+1) ].`

So higher support harmonics are viable only when the corresponding support lane is sufficiently soft.

This is the first exact microscopic support-feasibility threshold written directly in terms of the Stage-31 support demand.

---

## 5. Exact enhancement bounds for the higher harmonics

The physical enhancement produced by the `n`th twin harmonic is

`S_n^(twin)`
`= 1 + zeta_n^(twin) (1-eps) / (1 - eps zeta_n^(twin)).`

Because `S` is strictly increasing in `zeta`, the upper bound on `zeta_n^(twin)` gives the exact enhancement ceiling

`S_n^(twin) < S_n^(max)`
`:= 1 + (1-eps) / [ (2n+1)^2 - eps ].`

So, for example,

`S_1^(twin) < 1 + (1-eps)/(9-eps) = (10 - 2 eps)/(9 - eps),`

which is only a modest enhancement above `1`.

By contrast,

`S_0^(twin) = 2`

exactly.

So the lowest twin lane is not just the strongest support branch. It is qualitatively different from the higher harmonics in the size of the enhancement it can deliver.

---

## 6. Best current theorem statement after Stage 33

The support-compensation problem is now resolved down to an explicit harmonic-selection test on the first finite-throat coherent kernel.

### What is now exact

- the lowest symmetric twin lane has `zeta=1` and therefore `S=2` exactly,
- `zeta_req <= 1` is equivalent to `S_req <= 2`,
- every higher harmonic obeys the exact impossibility bound `zeta_req <= (2n+1)^(-2)`,
- when that bound is satisfied, the exact softness threshold is `x <= x_max(n;zeta_req)`,
- and the higher-harmonic support enhancement is strictly bounded above by `S_n^(max)`.

### What this means physically

The first explicit moving-throat kernel does more than produce a physical `zeta`.
It also orders the entire support tower.

- If the exact Stage-31 demand satisfies `zeta_req <= 1`, the symmetric lowest twin D/N lane is already sufficient.
- If `zeta_req > 1`, the lowest symmetric twin lane fails and stronger-than-twin asymmetry or a different kernel family is required.
- If `zeta_req > 1/9`, the first excited D/N support harmonic is already impossible.
- And the higher harmonics become rapidly less viable.

So the next honest PDE question is now extremely narrow:

> **does the completed moving-throat operator place the physical coherent support lane on the lowest twin branch with `zeta_req <= 1`, or does the real branch need non-twin asymmetry to overcome the exact support threshold?**

=== moving_throat_pde_stage034_lowest_twin_criterion.md ===

# Moving-Throat PDE — Stage 34: Exact Lowest-Twin Sufficiency Criterion on the Physical Tracking Branch

## Purpose

Stage 33 reduced the explicit D/N support comparison to one sharp question:

> does the physical moving-throat branch lie in the regime `zeta_req <= 1`, so that the symmetric lowest twin support lane is already sufficient?

That statement is still useful, but it is written in terms of the reduced support ratio `zeta_req`.
The next honest step is to eliminate `zeta_req` and translate the test into the continuum variables already frozen in Stages 30–33.

This stage does that exactly.

The main result is that the lowest symmetric twin lane is sufficient **iff**
the actual selected-branch point satisfies

`Pi_tr(xi_req,delta;R_tr) <= 16 Lambda (1-eps) / pi^2,`

where

`Pi_tr := F_tr G_tr`

is the exact tracking-branch product function. Equivalently, the full twin-sufficiency question collapses to a single radiative-demand inequality.

So the support problem is now no longer “do we have enough support in some vague sense?”
It is a one-line theorem test on the physical branch product.

---

## 1. Exact tracking-branch product

From the earlier stages,

`R_target = F_tr(xi_req,delta;R_tr),`

`M_req    = G_tr(xi_req,delta;R_tr),`

with

`G_tr(xi,delta;R)`
`= 9 xi (xi + delta) / [ 9 delta + (9 + 2 R^2) xi ],`

`F_tr(xi,delta;R)`
`= [ 9 delta + (9 + 2 R^2) xi ]^2 [ 9 delta + (9 + 2 R) xi ]^2`
`  / [ 81 (1 - xi) ( 9 delta^2 + 18 delta xi + (9 + 2 R^2) xi^2 )^2 ].`

Multiplying them gives the exact branch product

`Pi_tr(xi,delta;R)`
`:= F_tr(xi,delta;R) G_tr(xi,delta;R)`

`= xi (xi + delta) [ 9 delta + (9 + 2 R) xi ]^2 [ 9 delta + (9 + 2 R^2) xi ]`
`  / [ 9 (1 - xi) ( 9 delta^2 + 18 delta xi + (9 + 2 R^2) xi^2 )^2 ].`

This is the unique product that matters for the lowest-twin question.

Two endpoint facts are exact:

`Pi_tr(0,delta;R) = 0,`

`lim_(xi -> 1^-) Pi_tr(xi,delta;R) = +infinity.`

So for every finite microscopic support scale there is at least one stable-side depth where the lowest-twin criterion is saturated.

---

## 2. Elimination of `zeta_req`

Stage 31 gave the exact support-demand map

`S_req := M_req / M_mix,`

`zeta_req = (S_req - 1) / [ 1 + eps (S_req - 2) ],`

with the mixed-only product law from Stage 30

`R_target M_mix = 8 Lambda (1 - eps) / pi^2.`

Since

`R_target = F_tr(xi_req,delta;R_tr),`

`M_req    = G_tr(xi_req,delta;R_tr),`

it follows that

`S_req = M_req / M_mix`
`     = [ R_target M_req ] / [ R_target M_mix ]`
`     = Pi_tr(xi_req,delta;R_tr) / [ 8 Lambda (1 - eps) / pi^2 ].`

Now impose the exact Stage-33 twin threshold

`zeta_req <= 1 <=> S_req <= 2.`

Then the lowest symmetric twin lane is sufficient **iff**

`Pi_tr(xi_req,delta;R_tr) <= 16 Lambda (1 - eps) / pi^2.`

That is the sharpest version yet of the support theorem.

---

## 3. Exact threshold scales

The same inequality can be written as exact thresholds for several microscopic variables.

### Radiative-demand threshold

Define the exact lowest-twin required radiative scale

`Lambda_twin,req`
`:= (pi^2 / [16 (1 - eps)]) Pi_tr(xi_req,delta;R_tr).`

Then

`Lambda >= Lambda_twin,req`
`<=>`
`zeta_req <= 1.`

### Mixed-baseline threshold

Since `S_req <= 2` is equivalent to `M_req <= 2 M_mix`, the exact mixed-only baseline threshold is simply

`M_mix >= M_mix^(twin,req)`
`:= (1/2) G_tr(xi_req,delta;R_tr).`

So the lowest twin succeeds exactly when the mixed branch already supplies at least half of the total load demanded by the physical selected point.

### Wall-to-mixed overlap threshold

Using the Stage-30 coherent map

`M_mix = 8 Z_W (1 + chi_0)^2 / [ pi^2 (1 - eps_eta) (1 - eps) ],`

the exact wall/mixed overlap threshold is

`Z_W >= Z_W^(twin,req)`

with

`Z_W^(twin,req)`
`:= pi^2 (1 - eps_eta) (1 - eps) G_tr(xi_req,delta;R_tr)`
`   / [ 16 (1 + chi_0)^2 ].`

So even before the full PDE is solved, the lowest-twin support test is already equivalent to an explicit lower bound on the wall-to-mixed overlap strength.

---

## 4. Exact twin-saturation depth at fixed mixed baseline

The twin-sufficiency boundary corresponds to

`M_req = 2 M_mix,`

i.e.

`G_tr(xi,delta;R_tr) = 2 M_mix.`

Because `G_tr` is strictly increasing on the stable branch, this equation has a unique stable-side root. Solving the resulting quadratic gives

`xi_(2x)(M_mix,delta;R)`
`= [ 2 M_mix (9 + 2 R^2) - 9 delta`
`    + sqrt( [2 M_mix (9 + 2 R^2) - 9 delta]^2 + 648 M_mix delta ) ] / 18.`

This is the exact depth at which the lowest twin lane is just saturated by the mixed baseline.

It is useful because it shows the twin threshold is not a loose scaling estimate.
It lives on an exact algebraic branch of the same selected-mode family.

---

## 5. Best current theorem statement after Stage 34

The lowest-twin question is no longer phrased in terms of an abstract support ratio.

### What is exact now

- the physical selected-branch product is

  `Pi_tr = F_tr G_tr`,
- the mixed-only product scale is

  `8 Lambda (1 - eps) / pi^2`,
- the lowest symmetric twin lane succeeds **iff**

  `Pi_tr(xi_req,delta;R_tr) <= 16 Lambda (1 - eps) / pi^2`,
- the exact required radiative scale is

  `Lambda_twin,req = pi^2 Pi_tr / [16(1-eps)]`,
- the exact required mixed baseline is

  `M_mix^(twin,req) = G_tr/2`,
- and the exact twin-saturation depth at fixed mixed baseline is the closed quadratic root `xi_(2x)` above.

### What this means physically

The support problem has collapsed again.

The first explicit coherent twin lane is sufficient or insufficient according to a **single branch product inequality**.
So the next honest question is no longer “how do we parameterize support?”
It is:

> **does the completed moving-throat operator produce enough mixed baseline / radiative scale to satisfy the exact twin-sufficiency product test, or must the physical lowest support lane become non-twin before the normalization target can be met?**

=== moving_throat_pde_stage035_nontwin_asymmetry_threshold.md ===

# Moving-Throat PDE — Stage 35: Exact Non-Twin Asymmetry Requirement Once the Symmetric Lowest Twin Fails

## Purpose

Stage 34 compressed the lowest-twin question to one exact product test:

`Pi_tr(xi_req,delta;R_tr) <= 16 Lambda (1-eps) / pi^2.`

That tells us **when** the symmetric lowest twin lane is enough.
But the next physical question is just as important:

> if the symmetric lowest twin fails, what exact non-twin deformation is required to beat the support threshold?

This stage answers that.

The main result is that the support problem now splits into three exact regimes controlled by the same branch product `Pi_tr` and mixed product scale

`C_mix := 8 Lambda (1-eps) / pi^2.`

They are:

1. **mixed-only enough**

   `Pi_tr <= C_mix;`
2. **symmetric lowest twin enough**

   `C_mix < Pi_tr <= 2 C_mix;`
3. **non-twin asymmetry required**

   `Pi_tr > 2 C_mix.`

And once the third regime is entered, the exact required lowest-lane support ratio is

`zeta_req = (Pi_tr - C_mix) / [ C_mix - eps (2 C_mix - Pi_tr) ].`

So the residual theorem gap is no longer just “do we need asymmetry?”
It is the much sharper question:

> **what exact overlap boost and/or support softening is required once the symmetric twin doubling branch is no longer enough?**

---

## 1. Exact regime classification in branch-product form

Define again

`C_mix := 8 Lambda (1-eps) / pi^2,`

`Pi_tr := F_tr G_tr.`

Then the exact support-demand factor is

`S_req = Pi_tr / C_mix,`

and the exact required coherent support ratio is

`zeta_req = (S_req - 1) / [ 1 + eps (S_req - 2) ]`
`         = (Pi_tr - C_mix) / [ C_mix - eps (2 C_mix - Pi_tr) ].`

This formula is relevant once support is actually needed, i.e. once `Pi_tr > C_mix`.

The exact derivative is

`dzeta_req / dPi_tr`
`= C_mix (1-eps) / [ C_mix - eps (2 C_mix - Pi_tr) ]^2 > 0,`

so `zeta_req` grows strictly with the branch product.

The key anchor values are

`Pi_tr = C_mix   => zeta_req = 0,`

`Pi_tr = 2 C_mix => zeta_req = 1.`

So the support regimes are exact:

- `Pi_tr <= C_mix`: mixed-only branch already enough;
- `C_mix < Pi_tr <= 2 C_mix`: support is needed, but the symmetric lowest twin lane still suffices;
- `Pi_tr > 2 C_mix`: the symmetric lowest twin fails, and non-twin asymmetry is required.

This is the cleanest classification yet of the support problem.

---

## 2. Exact excess beyond the symmetric twin

When `Pi_tr > 2 C_mix`, the extra support demand beyond the symmetric twin is

`Delta_zeta := zeta_req - 1.`

Substituting the exact formula above gives

`Delta_zeta`
`= (1-eps) (Pi_tr - 2 C_mix)`
`  / [ C_mix - eps (2 C_mix - Pi_tr) ].`

So `Delta_zeta` is strictly positive exactly when the symmetric twin fails.

This means the “distance beyond the twin branch” is not heuristic.
It is an explicit rational function of the same branch product and the same blocking parameter `eps`.

---

## 3. General lowest-lane support ratio and the exact asymmetry threshold

The lowest support lane need not remain perfectly twin-symmetric.
For a general lowest support lane, define the relative overlap boost

`Omega_0 := I_(phi,0) / I_W,`

and keep the effective stiffness ratio explicit. Then the exact physical lowest-lane support ratio is

`zeta_0^(phys) = (K_W^(eff) / K_(phi,0)^(eff)) * Omega_0^2.`

So the exact success condition is

`(K_W^(eff) / K_(phi,0)^(eff)) * Omega_0^2 >= zeta_req.`

This yields two equivalent threshold forms.

### Overlap-asymmetry threshold at fixed stiffness

If the stiffness ratio is held fixed, the exact required overlap boost is

`Omega_0^2 >= Omega_(0,req)^2`
`:= zeta_req * K_(phi,0)^(eff) / K_W^(eff).`

### Support-softening threshold at fixed overlap

If the overlap ratio is held fixed, the exact required lowest-lane softening is

`K_(phi,0)^(eff) <= K_(phi,0)^(req)`
`:= K_W^(eff) * Omega_0^2 / zeta_req.`

So once `zeta_req > 1`, the lowest support lane must become non-twin in one of two exact ways:

- enhanced coupling to the lowest support mode (`Omega_0 > 1`),
- softer effective lowest-lane support stiffness (`K_(phi,0)^(eff) < K_W^(eff)`),
- or a combination of both.

---

## 4. Exact twin-failure diagnostics

The perfectly symmetric same-operator twin lane has

`Omega_0 = 1,`

`K_(phi,0)^(eff) = K_W^(eff),`

so

`zeta_0^(twin) = 1.`

If `zeta_req > 1`, this branch fails identically.

Two immediate exact diagnostics follow.

### Pure-overlap rescue at equal stiffness

If `K_(phi,0)^(eff) = K_W^(eff)`, then the exact overlap threshold is

`Omega_0 >= sqrt(zeta_req).`

So any `zeta_req > 1` forces a true overlap asymmetry above the symmetric twin value.

### Pure-softening rescue at equal overlap

If `Omega_0 = 1`, then the exact stiffness threshold is

`K_(phi,0)^(eff) <= K_W^(eff) / zeta_req.`

So any `zeta_req > 1` forces the lowest support lane to be physically softer than the mixed lane.

The exact softening fraction in that equal-overlap case is

`1 - K_(phi,0)^(req) / K_W^(eff)`
`= 1 - 1 / zeta_req`
`= (1-eps) (Pi_tr - 2 C_mix) / (Pi_tr - C_mix).`

So the size of the required non-twin softening is also explicit.

---

## 5. Best current theorem statement after Stage 35

The support bottleneck is no longer vague and no longer merely qualitative.

### What is exact now

- the required coherent support ratio is

  `zeta_req = (Pi_tr - C_mix) / [ C_mix - eps (2 C_mix - Pi_tr) ],`
- `zeta_req` is strictly increasing in the branch product `Pi_tr`,
- the exact regime split is

  `Pi_tr <= C_mix`  : mixed-only enough,

  `C_mix < Pi_tr <= 2 C_mix` : symmetric lowest twin enough,

  `Pi_tr > 2 C_mix` : non-twin asymmetry required,
- the exact excess beyond the symmetric twin is `Delta_zeta = zeta_req - 1`,
- and the exact lowest-lane rescue thresholds are

  `Omega_0^2 >= zeta_req K_(phi,0)^(eff) / K_W^(eff),`

  `K_(phi,0)^(eff) <= K_W^(eff) Omega_0^2 / zeta_req.`

### What this means physically

The first explicit coherent kernel has now answered the question posed at the end of Stage 33.

If the physical branch lands in the window `Pi_tr <= 2 C_mix`, then the symmetric lowest twin lane is enough.

If it lands above that, then the missing theorem gap is no longer “find some better support somehow.”
It is:

> **derive from the completed moving-throat operator whether the lowest support lane acquires the exact overlap boost and/or stiffness softening required by the non-twin threshold formulas above.**

=== moving_throat_pde_stage036_overlap_boost_window.md ===

# Moving-Throat PDE — Stage 36: Exact Overlap-Boost Window for the Lowest Support Lane

## Purpose

Stage 35 reduced the non-twin rescue problem to two concrete lowest-lane resources:

- overlap boost,
- support softening.

The next honest step is to ask whether the overlap side can be made explicit from a finite-throat operator instead of being left as an abstract factor
`Omega_0 = I_(phi,0)/I_W`.

This stage does that.

The main results are:

1. for any nonnegative lowest-lane source density `sigma_phi(s)` on the finite throat with the same total source strength as the mixed lane,
   the overlap boost satisfies the exact bound

   `0 <= Omega_0 <= pi/2,`

   so

   `A_I := Omega_0^2 <= pi^2/4;`

2. the symmetric uniform source is the baseline point `Omega_0 = 1`;
3. an explicit exponentially bottom-biased source family continuously deforms the overlap from `1` up to the sharp upper limit `pi/2`;
4. therefore **pure overlap asymmetry alone** can beat the Stage-35 support threshold only if

   `zeta_req <= pi^2/4.`

So this stage turns the phrase “maybe the support lane overlaps better” into an exact finite-throat theorem window.

---

## 1. Exact lowest-mode overlap for a general coherent source density

Keep the same D/N lowest mode on the finite throat interval `s in [0,L]`:

`chi_0(s) = sqrt(2/L) sin(pi s / (2L)).`

The mixed lane in Stage 32 used the uniform source density

`sigma_W(s) = 1,`

so its lowest-mode source overlap is

`I_W = int_0^L ds chi_0(s) = 2 sqrt(2L) / pi.`

Now let the physical lowest support lane couple through a general nonnegative source density `sigma_phi(s)` with the same total source strength,

`sigma_phi(s) >= 0,`

`int_0^L ds sigma_phi(s) = L.`

Define the support overlap

`I_(phi,0) := int_0^L ds sigma_phi(s) chi_0(s),`

and the overlap boost

`Omega_0 := I_(phi,0) / I_W.`

Then the lowest-lane coherent asymmetry factor from Stage 35 is

`A_I := Omega_0^2.`

---

## 2. Exact finite-throat overlap bound

Because `chi_0(s)` is nonnegative on `[0,L]` and satisfies

`0 <= chi_0(s) <= max chi_0 = sqrt(2/L),`

one has the exact bound

`I_(phi,0) <= (max chi_0) int_0^L ds sigma_phi(s) = L sqrt(2/L) = sqrt(2L).`

Dividing by the uniform mixed overlap gives

`Omega_0 <= sqrt(2L) / (2 sqrt(2L)/pi) = pi/2.`

So the overlap boost window is exactly

`0 <= Omega_0 <= pi/2,`

and therefore

`0 <= A_I <= pi^2/4.`

The upper bound is sharp: it is approached by a source density that concentrates at the antinode of `chi_0`, i.e. near the D/N bottom end of the finite throat.

---

## 3. Uniform branch and explicit constructive family

The symmetric mixed branch uses the uniform source density `sigma_W=1`, so

`Omega_0 = 1.`

To exhibit a constructive non-twin family, choose the bottom-biased exponential density

`sigma_alpha(s) = alpha exp(alpha s/L) / (exp(alpha) - 1),`

with `alpha >= 0`.

It has the same total source strength,

`int_0^L ds sigma_alpha(s) = L,`

and its exact overlap is

`I_alpha = int_0^L ds sigma_alpha(s) chi_0(s)`
`       = 2 sqrt(2L) alpha (2 alpha exp(alpha) + pi)`
`         / [ (4 alpha^2 + pi^2) (exp(alpha)-1) ].`

Therefore the exact overlap boost is

`Omega_exp(alpha)`
`= I_alpha / I_W`
`= pi alpha (2 alpha exp(alpha) + pi)`
`  / [ (4 alpha^2 + pi^2) (exp(alpha)-1) ].`

Its exact endpoint values are

`Omega_exp(0) = 1,`

`lim_(alpha -> +infinity) Omega_exp(alpha) = pi/2.`

So this explicit family interpolates continuously from the symmetric twin value to the sharp finite-throat upper bound.

For small asymmetry,

`Omega_exp(alpha)`
`= 1 + (2/pi - 1/2) alpha + O(alpha^2),`

and since `2/pi - 1/2 = (4-pi)/(2pi) > 0`, the constructive branch immediately moves upward from the symmetric point.

---

## 4. Exact pure-overlap rescue criterion

Stage 35 showed that, at equal stiffness,

`A_I = Omega_0^2 >= zeta_req`

is the exact rescue condition.

Stage 36 now bounds the left-hand side by

`A_I <= pi^2/4.`

Therefore **pure overlap asymmetry alone** is possible only if

`zeta_req <= pi^2/4.`

Equivalently:

- if `zeta_req <= pi^2/4`, a non-twin lowest lane can in principle beat the threshold using overlap enhancement alone;
- if `zeta_req > pi^2/4`, then overlap enhancement alone is impossible, and support softening must contribute.

This is the first exact no-go/sufficiency split for the overlap resource by itself.

---

## 5. Best current theorem statement after Stage 36

### What is exact now

- the finite-throat lowest-mode overlap boost satisfies

  `0 <= Omega_0 <= pi/2`,

- the corresponding asymmetry factor satisfies

  `0 <= A_I <= pi^2/4`,

- the symmetric uniform source gives `Omega_0 = 1`,
- the explicit exponential bottom-bias family gives

  `Omega_exp(alpha)`
  `= pi alpha (2 alpha exp(alpha) + pi) / [ (4 alpha^2 + pi^2)(exp(alpha)-1) ]`,

  with endpoints `1` and `pi/2`,
- and pure overlap rescue is possible only if

  `zeta_req <= pi^2/4.`

### What this means physically

The overlap side of the Stage-35 non-twin budget is no longer vague.

The moving-throat operator may indeed produce `Omega_0 > 1`, but even the most favorable finite-throat concentration can supply only a finite factor `pi^2/4` in `A_I`.

So the remaining question is already sharper:

> if the physical branch demands `zeta_req > pi^2/4`, then the lowest support lane must become softer as well; overlap enhancement alone cannot rescue it.

=== moving_throat_pde_stage037_robin_softening_support_lane.md ===

# Moving-Throat PDE — Stage 37: Robin-Compliance Softening of the Lowest Support Lane

## Purpose

Stage 36 put an exact ceiling on the overlap resource:

`A_I <= pi^2/4.`

So the next honest step is to make the second Stage-35 resource explicit:

> can the lowest support lane become softer than the mixed D/N lane by a concrete finite-throat mechanism, and by how much?

This stage answers that with the simplest explicit support-lane deformation:

- keep the same bulk operator and the same Neumann bottom,
- but replace the Dirichlet mouth by a finite Robin compliance.

The main results are:

1. the lowest support wavenumber is no longer fixed at `pi/(2L)` but is the unique root `y/L` of

   `y tan y = eta,`

   with `eta := h L` the dimensionless mouth-compliance parameter;
2. the support-lane stiffness becomes

   `K_(phi,0)^(eff) = K_X + T_X y^2 / L^2,`

   which is strictly smaller than the mixed-lane D/N value whenever `0 < eta < +infinity`;
3. the exact softening factor is

   `A_K(eta) := K_W^(eff) / K_(phi,0)^(eff)`
   `= 1 / [ 1 - x/4 + x y^2 / pi^2 ],`

   with `x := pi^2 T_X / (L^2 K_W^(eff))` and `0 < x < 4`;
4. it ranges from

   `A_K = 1` at the symmetric D/N point to the exact maximum

   `A_K,max = 4/(4-x)`

   in the soft-mouth limit `eta -> 0^+`;
5. therefore **pure support softening alone** can rescue the Stage-35 threshold only if

   `zeta_req <= 4/(4-x)`.

So this stage makes the support-softening side of the non-twin budget just as explicit as the overlap side.

---

## 1. Explicit compliant lowest support lane

Take the same finite interval `s in [0,L]` and the same bulk support operator, but impose

- Robin condition at the mouth `s=0`,
- Neumann condition at the bottom `s=L`.

Write the axial support mode as `psi(s)` satisfying

`psi'' + k^2 psi = 0,`

`psi'(0) = h psi(0),`

`psi'(L) = 0,`

with `h > 0` the mouth-compliance coefficient.

Solving the boundary-value problem gives the exact characteristic equation

`k tan(kL) = h.`

Defining the dimensionless variables

`y := kL,`

`eta := hL,`

the lowest support branch is the unique root

`y tan y = eta,`

with

`0 < y < pi/2`.

The symmetric D/N limit is recovered as `eta -> +infinity`, where `y -> pi/2`.
The fully soft-mouth limit is `eta -> 0^+`, where `y -> 0`.

---

## 2. Exact lowest-lane support stiffness and softening factor

Let the undeformed same-operator family have base stiffness `K_X` and axial tension `T_X`.

Then the mixed D/N lane has the frozen stiffness

`K_W^(eff) = K_X + pi^2 T_X / (4L^2),`

while the Robin-deformed lowest support lane has

`K_(phi,0)^(eff) = K_X + T_X y^2 / L^2.`

So the exact support-softening factor is

`A_K(eta) := K_W^(eff) / K_(phi,0)^(eff).`

It is convenient to introduce the dimensionless stiffness/tension ratio

`x := pi^2 T_X / (L^2 K_W^(eff)).`

Since `K_X > 0`, one has

`K_X = K_W^(eff) (1 - x/4),`

hence

`0 < x < 4.`

Substituting gives the exact reduced form

`A_K(eta)`
`= 1 / [ 1 - x/4 + x y(eta)^2 / pi^2 ].`

---

## 3. Exact softening window

The endpoint values are immediate.

### Symmetric D/N point

At `eta -> +infinity`, `y -> pi/2`, so

`A_K = 1.`

### Soft-mouth limit

At `eta -> 0^+`, `y -> 0`, so

`A_K,max = 1 / (1 - x/4) = 4 / (4 - x).`

Thus the exact softening window is

`1 <= A_K < 4/(4-x)`

for finite `eta > 0`, with the upper endpoint reached only in the closure `eta -> 0^+`.

Because the map `y -> eta = y tan y` is strictly increasing on `(0,pi/2)`, and `A_K` is strictly decreasing in `y`, the softening factor is strictly decreasing in `eta`.
So weaker mouth pinning means stronger support softening.

---

## 4. Exact pure-softening rescue criterion

Stage 35 showed that, at equal overlap,

`A_K >= zeta_req`

is the exact rescue condition.

Stage 37 now bounds the left-hand side by

`A_K <= 4/(4-x)`

(with closure equality as `eta -> 0^+`).

Therefore pure support softening alone is possible only if

`zeta_req <= 4/(4-x).`

Equivalently, for `zeta_req > 1`, the support-tension parameter must satisfy

`x >= 4 - 4/zeta_req.`

So once the demanded support ratio exceeds the symmetric twin value, the first explicit softening branch supplies an exact compliance floor.

---

## 5. Exact Robin threshold at fixed `zeta_req`

From

`A_K = 1 / [ 1 - x/4 + x y^2 / pi^2 ] >= zeta_req,`

one obtains the exact eigenvalue bound

`y^2 <= y_req^2`

with

`y_req^2 := (pi^2/x) ( 1/zeta_req - 1 + x/4 ).`

Therefore the exact Robin-compliance threshold is

`eta <= eta_req := y_req tan(y_req)`

whenever the right-hand side is nonnegative.

If

`1/zeta_req - 1 + x/4 < 0,`

then `y_req^2 < 0` and pure softening is impossible on this branch, regardless of the Robin parameter.

So the support-softening question is no longer qualitative either. It is an exact root-placement problem on the finite throat.

---

## 6. Best current theorem statement after Stage 37

### What is exact now

- the Robin-deformed lowest support branch is determined by

  `y tan y = eta`,

- the exact support-softening factor is

  `A_K = 1 / [ 1 - x/4 + x y^2 / pi^2 ]`,

- it obeys the exact window

  `1 <= A_K <= 4/(4-x)`

  (closure at `eta -> 0^+`),
- and pure softening rescue is possible only if

  `zeta_req <= 4/(4-x)`.

### What this means physically

The second Stage-35 resource is now explicit.

The lowest support lane can indeed become softer than the mixed D/N lane, but only within a finite operator-controlled window set by `x`.

So the remaining question is already narrower:

> if the physical branch needs more than the Stage-36 overlap ceiling `pi^2/4`, does the moving-throat operator also supply enough Robin-type support softening to bridge the remaining gap?

=== moving_throat_pde_stage038_explicit_lowest_lane_reachability.md ===

# Moving-Throat PDE — Stage 38: Explicit Reachability of the Non-Twin Lowest Support Lane

## Purpose

Stages 36 and 37 isolated the two exact lowest-lane resources:

- overlap boost, with ceiling `A_I <= pi^2/4`,
- Robin-compliance softening, with ceiling `A_K <= 4/(4-x)`.

The next honest question is therefore the combined one:

> for the first explicit moving-throat operator family that includes both effects, what is the exact reachable window of the physical lowest support ratio?

This stage answers that.

Using the explicit exponential source asymmetry from Stage 36 and the explicit Robin-compliance softening from Stage 37, the lowest support ratio is

`zeta_0^(exp+R)(alpha,eta)`
`= Omega_exp(alpha)^2 / [ 1 - x/4 + x y(eta)^2 / pi^2 ],`

where `y(eta)` solves `y tan y = eta`.

The main result is that this family has the exact closure range

`1 <= zeta_0^(exp+R) <= pi^2 / (4 - x),`

so the explicit family reaches a Stage-35 demand `zeta_req` **iff**

`zeta_req <= pi^2 / (4 - x)`

(in the closure; strictly `<` for finite `alpha,eta`).

This gives the first exact operator-level reachability theorem for the non-twin lowest support lane.

---

## 1. Explicit combined lowest-lane support ratio

Carry forward the explicit constructive pieces.

### Overlap branch

From Stage 36,

`Omega_exp(alpha)`
`= pi alpha (2 alpha exp(alpha) + pi)`
`  / [ (4 alpha^2 + pi^2) (exp(alpha)-1) ],`

with

`Omega_exp(0)=1,`

`lim_(alpha -> +infinity) Omega_exp(alpha)=pi/2.`

### Softening branch

From Stage 37,

`A_K(eta)`
`= 1 / [ 1 - x/4 + x y(eta)^2 / pi^2 ],`

with `y tan y = eta`,

`A_K(+infinity)=1,`

`lim_(eta -> 0^+) A_K(eta)=4/(4-x).`

So the explicit combined family is

`zeta_0^(exp+R)(alpha,eta)`
`= Omega_exp(alpha)^2 A_K(eta)`
`= Omega_exp(alpha)^2 / [ 1 - x/4 + x y(eta)^2 / pi^2 ].`

---

## 2. Exact closure range of the explicit family

The lower endpoint is the symmetric twin point:

`alpha = 0,`

`eta = +infinity,`

so

`zeta_0^(exp+R)=1.`

The upper endpoint is the closure of the combined constructive family:

`alpha -> +infinity,`

`eta -> 0^+,`

so

`zeta_0,max^(exp+R)`
`= (pi/2)^2 * 4/(4-x)`
`= pi^2 / (4 - x).`

Therefore

`1 <= zeta_0^(exp+R) <= pi^2/(4-x)`

in the closure of the family.

For finite `alpha` and finite positive `eta`, the inequality is strict at the top end.

---

## 3. Exact reachability criterion for the Stage-35 threshold

Stage 35 demands

`zeta_0^(phys) >= zeta_req.`

The explicit family can satisfy that exactly when

`zeta_req <= pi^2/(4-x)`

(in closure), equivalently

`x >= 4 - pi^2/zeta_req.`

This is the exact support-compliance floor for the first explicit non-twin lowest-lane family.

It immediately gives three exact sub-regimes.

### Regime A — overlap alone is enough

If

`zeta_req <= pi^2/4,`

then Stage 36 already shows that overlap enhancement alone can meet the target, so no compliance floor is required.

### Regime B — overlap ceiling exceeded but explicit combined family still works

If

`pi^2/4 < zeta_req <= pi^2/(4-x),`

then overlap enhancement alone is not enough, but the explicit combined family can still reach the target using both overlap and softening.

### Regime C — even the explicit combined family fails

If

`zeta_req > pi^2/(4-x),`

then neither exponential overlap enhancement nor Robin-compliance softening, nor any combination of the two within this first explicit family, can beat the threshold. A stronger operator deformation would be required.

---

## 4. Equivalent stiffness-ratio form

Since

`K_X = K_W^(eff) (1 - x/4),`

the reachability floor

`x >= 4 - pi^2/zeta_req`

is equivalent to

`K_X / K_W^(eff) <= pi^2 / (4 zeta_req).`

So once the support demand exceeds the pure-overlap ceiling, the explicit combined family requires the background support stiffness to sit below a precise fraction of the mixed-lane stiffness.

That is the cleanest mechanical interpretation of the combined theorem.

---

## 5. Best current theorem statement after Stage 38

### What is exact now

- the first explicit non-twin lowest-lane family is

  `zeta_0^(exp+R)(alpha,eta)`
  `= Omega_exp(alpha)^2 / [ 1 - x/4 + x y(eta)^2 / pi^2 ],`

  with `y tan y = eta`,
- its closure range is

  `1 <= zeta_0^(exp+R) <= pi^2/(4-x)`,

- and it reaches the Stage-35 threshold exactly when

  `zeta_req <= pi^2/(4-x)`

  (strict `<` for finite parameters).

### What this means physically

The support question is no longer abstract even at the operator level.

For the first explicit constructive moving-throat family, the non-twin lowest support lane succeeds only inside one exact reachability window.

So the remaining gap is now extremely narrow:

> does the completed moving-throat PDE generate a physical lowest-lane deformation whose effective `x` and source-shape asymmetry place it inside the exact Stage-38 reachability window, or does the real branch require an even stronger non-twin mechanism than exponential overlap bias plus Robin-compliance softening?

=== moving_throat_pde_stage039_transport_source_asymmetry.md ===

# Moving-Throat PDE — Stage 39: Transport Origin of the Lowest-Lane Source Asymmetry

## Purpose

Stage 38 left one sharp operator-level question open:

> what is the **physical** origin of the source-shape asymmetry parameter that had previously been written as the abstract exponential-bias variable `alpha`?

This stage answers that with the simplest conservative axial source-transport law on the finite throat.

The main result is that the Stage-36 exponential family is not ad hoc. It is the exact stationary zero-flux branch of a drift-diffusion transport operator,

`partial_t sigma + partial_s J = 0,`

`J = -D_sigma partial_s sigma + v_sigma sigma,`

on the finite throat interval `s in [0,L]`.

On the stationary recirculating branch `J=0`, the normalized source density is

`sigma_Pe(s)`
`= Pe exp(Pe s/L) / (exp(Pe)-1),`

with the single physical asymmetry parameter

`Pe := v_sigma L / D_sigma.`

So the old abstract bias parameter is exactly the axial **Peclet number** of the lowest support-source transport problem.

This converts the overlap side of Stages 36–38 from a free deformation into a concrete operator output.

---

## 1. Explicit axial source-transport operator

Let `sigma(s,t)` be the coherent source density feeding the lowest support lane along the finite throat interval `s in [0,L]`.

Take the minimal conservative axial transport law

`partial_t sigma + partial_s J = 0,`

`J = -D_sigma partial_s sigma + v_sigma sigma,`

with

- `D_sigma > 0` the axial spreading coefficient,
- `v_sigma` the directed recirculation / drift speed along the throat.

This is the simplest operator that can bias the support source toward one end of the throat without introducing explicit nonconservative loss.

---

## 2. Stationary recirculating branch and exact exponential profile

On a stationary branch,

`partial_t sigma = 0,`

so `J` is constant in `s`.

For the lowest closed recirculating support branch, the natural condition is zero net axial source flux,

`J = 0.`

Then the transport equation collapses exactly to

`D_sigma sigma' = v_sigma sigma,`

hence

`sigma(s) = C exp(v_sigma s / D_sigma).`

Normalizing to the same total source strength used in Stage 36,

`int_0^L ds sigma(s) = L,`

gives

`sigma_Pe(s)`
`= Pe exp(Pe s/L) / (exp(Pe)-1),`

with

`Pe := v_sigma L / D_sigma.`

So the Stage-36 constructive family is exactly the stationary zero-flux branch of the minimal transport operator.

### Physical interpretation of the sign

- `Pe = 0`: symmetric uniform source, i.e. the twin baseline;
- `Pe > 0`: constructive branch, with source weight shifted toward the D/N bottom antinode;
- `Pe < 0`: destructive branch, with source weight shifted toward the mouth where the D/N mode is smallest.

So the same operator explains both the helpful and harmful non-twin deformations.

---

## 3. Exact overlap boost on the transport branch

Keep the D/N lowest mode

`chi_0(s) = sqrt(2/L) sin(pi s/(2L)).`

Its mixed-lane uniform overlap is

`I_W = int_0^L ds chi_0(s) = 2 sqrt(2L)/pi.`

For the transport branch,

`I_Pe = int_0^L ds sigma_Pe(s) chi_0(s)`

evaluates exactly to

`I_Pe = 2 sqrt(2L) Pe (2 Pe exp(Pe) + pi)`
`       / [ (4 Pe^2 + pi^2) (exp(Pe)-1) ].`

Therefore the physical overlap boost is

`Omega_Pe := I_Pe / I_W`
`= pi Pe (2 Pe exp(Pe) + pi)`
`  / [ (4 Pe^2 + pi^2) (exp(Pe)-1) ].`

So the abstract Stage-36 factor `Omega_exp(alpha)` is now identified exactly as `Omega_Pe` on a concrete operator branch.

---

## 4. Exact monotonicity identity on the constructive branch

Introduce the probability density

`p_Pe(s) := sigma_Pe(s)/L.`

Then

`I_Pe = L E_Pe[chi_0],`

and differentiation gives the exact covariance identity

`dI_Pe/dPe = Cov_Pe(chi_0,s),`

hence

`dOmega_Pe/dPe = Cov_Pe(chi_0,s) / I_W.`

Because `chi_0(s)` is strictly increasing on `[0,L]`, the covariance is positive on the constructive branch.

Therefore

`dOmega_Pe/dPe > 0`

for `Pe >= 0`.

So the physical transport branch is not merely continuous. It is an exact **monotone** route from the twin point toward maximal overlap.

This is important later because it means the required transport bias is unique once the target overlap is known.

---

## 5. Endpoint and asymptotic structure

The exact endpoint values are

`Omega_Pe(0) = 1,`

`lim_(Pe -> +infinity) Omega_Pe = pi/2.`

So the constructive transport branch reproduces the full Stage-36 overlap window

`1 <= Omega_Pe <= pi/2.`

Useful asymptotics are:

### Small transport bias

`Omega_Pe`
`= 1 + ((4-pi)/(2 pi)) Pe + O(Pe^2).`

So the constructive branch leaves the symmetric point with a strictly positive linear slope.

### Strong transport bias

`Omega_Pe`
`= pi/2 - pi^3/(8 Pe^2) + O(Pe^-3).`

So the approach to the sharp finite-throat ceiling is algebraic, not exponential.

---

## 6. Best current theorem statement after Stage 39

### What is exact now

- the Stage-36 exponential source family is the exact stationary zero-flux solution of

  `partial_t sigma + partial_s(-D_sigma sigma' + v_sigma sigma) = 0`,

- the asymmetry parameter is the physical Peclet number

  `Pe = v_sigma L / D_sigma`,

- the overlap boost on that branch is

  `Omega_Pe`
  `= pi Pe (2 Pe exp(Pe) + pi) / [ (4 Pe^2 + pi^2)(exp(Pe)-1) ]`,

- and on the constructive branch `Pe >= 0` it is strictly increasing from

  `Omega_Pe(0)=1`

  to

  `lim_(Pe->+infinity) Omega_Pe = pi/2`.

### What this means physically

The source-shape asymmetry is no longer an abstract deformation parameter.
It is the axial transport Peclet number of the lowest support-source channel.

So one of the two Stage-38 “unknown physical inputs” has now been converted into a concrete moving-throat operator datum.

The remaining operator-level task is to combine that transport bias with the physical support-compliance ratios, so that the whole lowest-lane reachability problem is written directly in terms of real throat operator parameters rather than the old abstract pair `(alpha,x)`.

=== moving_throat_pde_stage040_physical_parameter_map.md ===

# Moving-Throat PDE — Stage 40: Physical `(Pe, kappa, eta)` Placement Map for the Lowest Support Lane

## Purpose

Stage 39 converted the old abstract overlap-bias variable `alpha` into the physical transport Peclet number

`Pe = v_sigma L / D_sigma.`

The last remaining Stage-38 unknown is therefore the support-compliance variable previously written as `x`.

This stage converts the whole explicit non-twin family into a directly physical parameter map.

The main results are:

1. the old abstract support parameter is exactly

   `x = pi^2 / (kappa + pi^2/4),`

   where

   `kappa := K_X L^2 / T_X`

   is the baseline support stiffness ratio;
2. if the Robin mouth coefficient is written as

   `eta := hL = K_m L / T_X`,

   with `y(eta)` the unique root of `y tan y = eta`, then the support softening factor becomes

   `A_K(eta;kappa) = (kappa + pi^2/4)/(kappa + y(eta)^2);`
3. the physical explicit lowest-lane ratio is therefore

   `zeta_0^(Pe+R)(Pe,eta;kappa)`
   `= Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2);`
4. on the constructive branch `Pe >= 0`, this map is monotone

   - increasing in `Pe`,
   - decreasing in `eta`,
   - decreasing in `kappa`;
5. its exact closure ceiling is

   `zeta_max(kappa) = (pi^2/4)(kappa + pi^2/4)/kappa;`
6. therefore the Stage-35 demand is reachable on this first physical family **iff**

   `zeta_req <= zeta_max(kappa),`

   equivalently, when `zeta_req > pi^2/4`,

   `kappa <= pi^4 / [4(4 zeta_req - pi^2)].`

So the Stage-38 reachability problem is now fully rewritten in terms of three physical throat-operator ratios:

- `Pe` — axial support-source transport bias,
- `eta` — mouth Robin compliance,
- `kappa` — baseline support stiffness ratio.

---

## 1. Exact physical support parameters

Carry forward the finite-throat support operator with base stiffness `K_X`, axial tension `T_X`, and Robin mouth coefficient `h`.

Define the two physical dimensionless support ratios

`kappa := K_X L^2 / T_X,`

`eta := hL = K_m L / T_X`

(the last equality if the Robin boundary arises from a mouth spring `K_m`).

The mixed D/N lane has

`K_W^(eff) = K_X + pi^2 T_X/(4L^2)`
`          = (T_X/L^2) (kappa + pi^2/4),`

so the Stage-37 support parameter becomes

`x = pi^2 T_X / (L^2 K_W^(eff))`
`  = pi^2 / (kappa + pi^2/4).`

Thus the old abstract `x` is exactly the inverse stiffness ratio of the physical support operator.

---

## 2. Exact Robin softening in physical form

Let `y(eta)` be the unique root of

`y tan y = eta,`

with `0 < y < pi/2`.

Then the lowest Robin support lane has

`K_(phi,0)^(eff) = K_X + T_X y^2 / L^2`
`                = (T_X/L^2) (kappa + y^2).`

So the exact softening factor simplifies to

`A_K(eta;kappa)`
`= K_W^(eff)/K_(phi,0)^(eff)`
`= (kappa + pi^2/4)/(kappa + y(eta)^2).`

This is algebraically cleaner than the Stage-37 `x`-form and makes the physical support meaning explicit.

---

## 3. Exact physical explicit lowest-lane family

Combining Stage 39 with the Robin support branch gives the first fully physical explicit family:

`zeta_0^(Pe+R)(Pe,eta;kappa)`
`= Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2),`

where

`Omega_Pe`
`= pi Pe (2 Pe exp(Pe) + pi)`
`  / [ (4 Pe^2 + pi^2)(exp(Pe)-1) ]`

and `y tan y = eta`.

This is the direct physical replacement for the Stage-38 abstract family

`zeta_0^(exp+R)(alpha,eta)`
`= Omega_exp(alpha)^2 / [1 - x/4 + x y(eta)^2/pi^2].`

The old variables are now identified as

`alpha = Pe,`

`x = pi^2/(kappa + pi^2/4).`

---

## 4. Exact monotone placement map

The physical family has a clean monotone structure.

### Increasing in transport bias `Pe`

From Stage 39,

`dOmega_Pe/dPe > 0`

on the constructive branch `Pe >= 0`, so

`partial_Pe zeta_0^(Pe+R) > 0.`

### Decreasing in compliance parameter `eta`

Because `y tan y = eta` has

`dy/deta = 1 / (tan y + y sec^2 y) > 0,`

one has

`partial_eta zeta_0^(Pe+R) < 0.`

So weaker mouth pinning (smaller `eta`) always helps the support lane.

### Decreasing in stiffness ratio `kappa`

Differentiating directly gives

`partial_kappa zeta_0^(Pe+R)`
`= Omega_Pe^2 ( y^2 - pi^2/4 ) / (kappa + y^2)^2 < 0,`

since `0 < y < pi/2`.

So a softer baseline support branch (smaller `kappa`) always helps.

This means the physical branch placement is completely ordered:

- more axial source drift helps,
- more mouth compliance helps,
- less baseline support stiffness helps.

---

## 5. Exact closure window in physical form

On the constructive branch `Pe >= 0`, the exact baseline point is

`Pe = 0,`

`eta = +infinity,`

so

`zeta_0^(Pe+R) = 1.`

The closure ceiling is reached at

`Pe -> +infinity,`

`eta -> 0^+,`

which gives

`zeta_max(kappa)`
`= (pi^2/4) (kappa + pi^2/4) / kappa.`

Therefore the exact constructive-branch window is

`1 <= zeta_0^(Pe+R)(Pe,eta;kappa) <= zeta_max(kappa)`

in closure.

This is the first explicit reachability ceiling written entirely in physical operator variables.

---

## 6. Exact reachability criterion and stiffness ceiling

The Stage-35 support demand is reachable on this physical family exactly when

`zeta_req <= zeta_max(kappa)`

that is,

`zeta_req <= (pi^2/4)(kappa + pi^2/4)/kappa.`

When `zeta_req > pi^2/4`, this is equivalent to the exact stiffness ceiling

`kappa <= kappa_max(zeta_req)`
`:= pi^4 / [4(4 zeta_req - pi^2)].`

So above the pure-overlap ceiling, the question is no longer “maybe enough compliance exists.”
It is the exact physical inequality above.

---

## 7. Exact parameter thresholds inside the physical family

Because the map is monotone in each physical parameter, the threshold surfaces are exact.

### Required overlap/transport threshold at fixed `(kappa,eta)`

The support demand requires

`Omega_Pe^2 >= Omega_req^2(kappa,eta;zeta_req)`

with

`Omega_req^2 := zeta_req (kappa + y(eta)^2)/(kappa + pi^2/4).`

Since the constructive `Omega_Pe` branch runs continuously from `1` to `pi/2`, a physical transport solution exists whenever

`1 <= Omega_req^2 <= pi^2/4`,

and the least constructive transport bias is the least nonnegative root of

`Omega_Pe^2 = Omega_req^2.`

### Required compliance threshold at fixed `(Pe,kappa)`

Equivalently,

`y(eta)^2 <= y_req^2(Pe,kappa;zeta_req)`

with

`y_req^2 := (Omega_Pe^2/zeta_req)(kappa + pi^2/4) - kappa.`

So if `0 <= y_req < pi/2`, the exact Robin threshold is

`eta <= eta_req := y_req tan(y_req).`

### Required stiffness threshold at fixed `(Pe,eta)`

If overlap alone is not already enough, i.e. if `Omega_Pe^2 < zeta_req`, then the exact stiffness ceiling is

`kappa <= kappa_req(Pe,eta;zeta_req)`

with

`kappa_req`
`= [ Omega_Pe^2 pi^2/4 - zeta_req y(eta)^2 ] / [ zeta_req - Omega_Pe^2 ],`

provided the numerator is nonnegative.

So all three physical threshold surfaces are now explicit.

---

## 8. Best current theorem statement after Stage 40

### What is exact now

- the old abstract overlap parameter is the physical transport Peclet number
  `Pe = v_sigma L / D_sigma`,
- the old abstract support parameter is
  `x = pi^2/(kappa + pi^2/4)`
  with `kappa = K_X L^2/T_X`,
- the Robin mouth parameter is
  `eta = hL = K_m L / T_X`,
- the physical explicit lowest-lane family is

  `zeta_0^(Pe+R)(Pe,eta;kappa)`
  `= Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2)`,

- this map is monotone increasing in `Pe` and monotone decreasing in both `eta` and `kappa`,
- and the exact constructive-branch ceiling is

  `zeta_max(kappa) = (pi^2/4)(kappa + pi^2/4)/kappa`.

### What this means physically

The lowest-lane support problem is no longer formulated in abstract deformation parameters at all.

It has collapsed to three concrete moving-throat operator ratios:

- axial source Peclet number `Pe`,
- mouth compliance `eta`,
- baseline support stiffness ratio `kappa`.

So the remaining gap is now as narrow as it can be without the full PDE:

> compute the physical branch point `(Pe, eta, kappa)` from the completed moving-throat operator and check whether it satisfies the exact threshold `zeta_req <= zeta_max(kappa)` and the corresponding monotone placement inequalities above.

=== moving_throat_pde_stage041_coupled_support_source_operator.md ===

# Moving-Throat PDE — Stage 41: Coupled Support/Source Operator and the Exact `Pe` Branch Equation

## Purpose

Stage 40 rewrote the explicit lowest-lane reachability problem in terms of the physical variables

`Pe, eta, kappa.`

But those were still being treated as independent branch inputs.

The next honest derivation is therefore:

> build the first explicit **coupled** axial support/source operator that generates the transport bias `Pe` from the same support lane that already carries the Robin compliance `eta` and the baseline stiffness ratio `kappa`.

This stage does that.

The main results are:

1. the minimal coupled axial operator is

   `partial_t sigma + partial_s J = 0,`

   `J = -D_sigma partial_s sigma + v_sigma sigma,`

   `-T_X partial_s^2 phi + K_X phi = Lambda_phi sigma,`

   with the support boundary conditions

   `T_X phi_s(0) = K_m phi(0),`

   `phi_s(L) = 0;`
2. after nondimensionalization,

   `kappa := K_X L^2 / T_X,`

   `eta := K_m L / T_X,`

   `Pe := v_sigma L / D_sigma,`

   and the new coupled-strength parameter is

   `Xi := mu_sigma Lambda_phi^2 L^2 / (D_sigma T_X);`
3. on the stationary zero-flux branch, the source density is still the exact drift-diffusion family

   `Sigma_Pe(x) = Pe exp(Pe x)/(exp(Pe)-1),`

   now with `x=s/L in [0,1]`;
4. the support field generates an exact dimensionless end-to-end drop

   `Delta(Pe; kappa, eta) = Phi(1)-Phi(0)`

   with kernel

   `K_(kappa,eta)(x)`
   `= [ cosh(alpha x) + (eta/alpha) sinh(alpha x) - cosh(alpha(1-x)) ]`
   `  / [ alpha sinh(alpha) + eta cosh(alpha) ],`

   where `alpha := sqrt(kappa)`;
5. the drift closure produces the exact fixed-point equation

   `Pe = Xi Delta(Pe; kappa, eta);`
6. the kernel is strictly increasing,

   `dK_(kappa,eta)/dx`
   `= [ alpha sinh(alpha x) + eta cosh(alpha x) + alpha sinh(alpha(1-x)) ]`
   `  / [ alpha sinh(alpha) + eta cosh(alpha) ] > 0,`

   so `Delta(Pe;kappa,eta)` is strictly increasing on the constructive branch by the same covariance identity used in Stage 39;
7. therefore every constructive branch point `Pe_*` obeys the exact bracket

   `Xi Delta_0(kappa,eta) <= Pe_* <= Xi Delta_inf(kappa,eta),`

   where

   `Delta_0(kappa,eta)`
   `= eta (cosh(alpha)-1) / [ alpha^2 ( alpha sinh(alpha) + eta cosh(alpha) ) ],`

   `Delta_inf(kappa,eta)`
   `= [ cosh(alpha) + (eta/alpha) sinh(alpha) - 1 ]`
   `  / [ alpha sinh(alpha) + eta cosh(alpha) ].`

So the branch point is no longer an uncontrolled free variable. It is the root of one explicit support/source fixed-point equation, and even before solving it numerically it is trapped inside one exact operator-determined interval.

---

## 1. Minimal coupled axial support/source operator

Let `s in [0,L]` be the throat axis.

Take the coherent source density `sigma(s,t)` to obey the conservative transport law

`partial_t sigma + partial_s J = 0,`

`J = -D_sigma partial_s sigma + v_sigma sigma.`

This is exactly the Stage-39 drift-diffusion operator.

Now add the lowest support field `phi(s,t)` through the minimal static axial support law

`-T_X partial_s^2 phi + K_X phi = Lambda_phi sigma,`

with support boundary conditions

`T_X phi_s(0) = K_m phi(0),`

`phi_s(L) = 0.`

The physical interpretation is straightforward:

- `T_X` is the axial support tension,
- `K_X` is the baseline support stiffness,
- `K_m` is the mouth spring / Robin compliance,
- `Lambda_phi` is the source-to-support loading strength.

To close the drift, use the minimal averaged support-force law

`v_sigma = mu_sigma Lambda_phi [ phi(L)-phi(0) ] / L,`

with `mu_sigma` the transport mobility.

Then the axial drift is no longer independent of the support field. It is generated by the same lane that determines `eta` and `kappa`.

---

## 2. Exact nondimensionalization

Set `x=s/L in [0,1]` and define the dimensionless support field

`Phi(x) := T_X phi(s) / (Lambda_phi L^2).`

Then

`-Phi_xx + kappa Phi = Sigma,`

with

`Phi_x(0) = eta Phi(0),`

`Phi_x(1)=0,`

and

`kappa := K_X L^2 / T_X,`

`eta := K_m L / T_X.`

The transport side keeps the physical Peclet number

`Pe := v_sigma L / D_sigma.`

Using the support-driven drift law gives one additional exact dimensionless coupling

`Xi := mu_sigma Lambda_phi^2 L^2 / (D_sigma T_X).`

The closure is therefore

`Pe = Xi [ Phi(1)-Phi(0) ].`

So the three Stage-40 variables are no longer peers. `eta` and `kappa` remain support-shape data, but `Pe` is now an output of the same coupled operator through the single coupling strength `Xi`.

---

## 3. Exact support-drop kernel

On the stationary zero-flux branch,

`Sigma_Pe(x) = Pe exp(Pe x)/(exp(Pe)-1).`

The Green-kernel difference controlling the support drop is

`K_(kappa,eta)(x)`
`= G(1,x)-G(0,x)`
`= [ cosh(alpha x) + (eta/alpha) sinh(alpha x) - cosh(alpha(1-x)) ]`
`  / [ alpha sinh(alpha) + eta cosh(alpha) ],`

where `alpha := sqrt(kappa)`.

Therefore the exact dimensionless support drop is

`Delta(Pe;kappa,eta)`
`:= Phi(1)-Phi(0)`
`= int_0^1 dx K_(kappa,eta)(x) Sigma_Pe(x).`

For later explicit work it is convenient to write this in closed form using

`I_c(Pe,alpha)`
`= int_0^1 dx exp(Pe x) cosh(alpha x)`
`= [ exp(Pe)(Pe cosh(alpha)-alpha sinh(alpha)) - Pe ] / (Pe^2-alpha^2),`

`I_s(Pe,alpha)`
`= int_0^1 dx exp(Pe x) sinh(alpha x)`
`= [ exp(Pe)(Pe sinh(alpha)-alpha cosh(alpha)) + alpha ] / (Pe^2-alpha^2).`

Then

`Delta(Pe;kappa,eta)`
`= Pe/(exp(Pe)-1)`
`  * [ (1-cosh(alpha)) I_c + (eta/alpha + sinh(alpha)) I_s ]`
`    / [ alpha sinh(alpha) + eta cosh(alpha) ].`

This is the first exact coupled-operator formula for the transport-selected branch point.

---

## 4. Exact monotonicity and endpoint data

The kernel derivative is exact:

`dK_(kappa,eta)/dx`
`= [ alpha sinh(alpha x) + eta cosh(alpha x) + alpha sinh(alpha(1-x)) ]`
`  / [ alpha sinh(alpha) + eta cosh(alpha) ] > 0`

for `alpha>0`, `eta>0`, and `x in [0,1]`.

So the support-drop kernel is strictly increasing toward the D/N bottom end.

Because `Sigma_Pe` is the same exponential family as in Stage 39,

`dDelta/dPe = Cov_Pe(K_(kappa,eta), x) > 0`

on the constructive branch `Pe>=0`.

So the coupled support drop is strictly increasing with the source-transport bias.

Two exact endpoint values are especially useful.

### Uniform-source point

At `Pe=0`,

`Delta_0(kappa,eta)`
`:= Delta(0;kappa,eta)`
`= int_0^1 dx K_(kappa,eta)(x)`
`= eta (cosh(alpha)-1)`
`  / [ alpha^2 ( alpha sinh(alpha) + eta cosh(alpha) ) ].`

This is strictly positive for `eta>0`.

### Sharp-bottom-source point

As `Pe -> +infinity`, the source weight concentrates at `x=1`, so

`Delta_inf(kappa,eta)`
`:= lim_(Pe->+infinity) Delta(Pe;kappa,eta)`
`= K_(kappa,eta)(1)`
`= [ cosh(alpha) + (eta/alpha) sinh(alpha) - 1 ]`
`  / [ alpha sinh(alpha) + eta cosh(alpha) ].`

So the exact support-drop window is

`Delta_0(kappa,eta) <= Delta(Pe;kappa,eta) <= Delta_inf(kappa,eta)`

on the constructive branch.

---

## 5. Exact `Pe` branch equation and operator-selected interval

The coupled support/source closure is now one exact fixed-point equation:

`Pe = Xi Delta(Pe;kappa,eta).`

Because `Delta(Pe;kappa,eta)` is continuous, positive, and bounded by `Delta_inf`, every constructive branch point obeys

`Pe_* = Xi Delta(Pe_*;kappa,eta)`

and therefore

`Xi Delta_0(kappa,eta) <= Pe_* <= Xi Delta_inf(kappa,eta).`

The proof is immediate.

Define

`F(Pe) := Pe - Xi Delta(Pe;kappa,eta).`

Then

`F(Xi Delta_0) <= 0,`

because `Delta(Pe) >= Delta_0`, and

`F(Xi Delta_inf) >= 0,`

because `Delta(Pe) <= Delta_inf`.

So by continuity there is at least one constructive root in the interval above.

This is the first exact branch-point bracket derived from a coupled operator rather than imposed as a guessed placement.

---

## 6. Weak-coupling branch law

For small `Xi`, the fixed-point equation immediately gives

`Pe_* = Xi Delta_0(kappa,eta) + O(Xi^2).`

So the first physical branch leaves the symmetric point in a completely explicit way.

The lowest transport bias is not arbitrary. It is proportional to the uniform-source support drop of the same support operator.

---

## 7. Best current theorem statement after Stage 41

### What is exact now

- the first coupled axial support/source operator has been written explicitly;
- the source-transport bias `Pe` is no longer an independent knob;
- it is the root of the exact support/source fixed-point equation

  `Pe = Xi Delta(Pe;kappa,eta)`;
- the support-drop kernel is exact and strictly increasing;
- and every constructive physical branch point obeys the operator-selected interval

  `Xi Delta_0(kappa,eta) <= Pe_* <= Xi Delta_inf(kappa,eta)`.

### What this means physically

The Stage-40 map has now been upgraded from a placement map on three external variables to a partially self-consistent branch law.

The open problem is no longer “what abstract `Pe` should we try?”
It is:

> for the real moving-throat operator, what is the actual coupling strength `Xi`, and where inside the exact interval `[Xi Delta_0, Xi Delta_inf]` does the selected branch point `Pe_*` land?

=== moving_throat_pde_stage042_operator_branch_residual_bounds.md ===

# Moving-Throat PDE — Stage 42: Exact Residual Bounds on the Operator-Selected Branch

## Purpose

Stage 41 turned the transport bias `Pe` into the root of a coupled support/source fixed-point equation,

`Pe = Xi Delta(Pe;kappa,eta),`

and trapped every constructive branch point inside the exact interval

`Xi Delta_0(kappa,eta) <= Pe_* <= Xi Delta_inf(kappa,eta).`

That means the Stage-40 physical placement map can now be evaluated on a **true operator-selected branch** rather than on an arbitrary external `Pe`.

This stage uses that bracket to derive exact success/no-go tests for the normalization residual without solving the full root.

The main results are:

1. the operator-selected physical lowest-lane ratio is

   `zeta_phys(Xi,eta;kappa)`
   `= zeta_0^(Pe_*+R)`
   `= Omega_(Pe_*)^2 (kappa + pi^2/4)/(kappa + y(eta)^2),`

   with `y tan y = eta`;
2. because `zeta_0^(Pe+R)` is monotone increasing in `Pe`, the exact Stage-41 branch bracket immediately gives

   `zeta_-(Xi,eta;kappa) <= zeta_phys(Xi,eta;kappa) <= zeta_+(Xi,eta;kappa),`

   where

   `zeta_- = Omega_(Xi Delta_0)^2 (kappa + pi^2/4)/(kappa + y^2),`

   `zeta_+ = Omega_(Xi Delta_inf)^2 (kappa + pi^2/4)/(kappa + y^2);`
3. therefore the exact residual

   `R_phys := zeta_req - zeta_phys`

   obeys

   `R_- <= R_phys <= R_+,`

   with

   `R_- = zeta_req - zeta_+,`

   `R_+ = zeta_req - zeta_-;`
4. this yields two exact theorem gates **before** solving the full branch root:

   - guaranteed success if `zeta_- >= zeta_req`,
   - guaranteed failure inside this operator family if `zeta_+ < zeta_req`;
5. inside the reachable Stage-40 window, define the unique constructive branch point `Pe_req` by

   `Omega_(Pe_req)^2 = zeta_req (kappa + y^2)/(kappa + pi^2/4)`

   (`Pe_req` exists exactly when the Stage-40 ceiling holds);
6. then the exact operator-coupling thresholds are

   `Xi_fail = Pe_req / Delta_inf(kappa,eta),`

   `Xi_suff = Pe_req / Delta_0(kappa,eta),`

   and they satisfy

   `Xi_fail <= Xi_suff`;
7. so the coupled operator now has a precise three-zone structure:

   - `Xi <= Xi_fail` : impossible,
   - `Xi >= Xi_suff` : guaranteed,
   - `Xi_fail < Xi < Xi_suff` : full root solve needed, but only inside this narrow interval.

This is the first exact residual-bracketing theorem on the operator-selected branch.

---

## 1. Operator-selected physical lowest-lane ratio

Carry forward the Stage-40 physical family

`zeta_0^(Pe+R)(Pe,eta;kappa)`
`= Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2),`

with `y tan y = eta` and the constructive-branch monotonicity

`dOmega_Pe/dPe > 0.`

Now evaluate it on the Stage-41 branch point `Pe_*`, defined implicitly by

`Pe_* = Xi Delta(Pe_*;kappa,eta).`

The operator-selected physical support ratio is therefore

`zeta_phys(Xi,eta;kappa)`
`= Omega_(Pe_*)^2 (kappa + pi^2/4)/(kappa + y(eta)^2).`

This is the first branch-selected version of the Stage-40 placement map.

---

## 2. Exact branch brackets for the physical ratio

Because Stage 41 proved

`Xi Delta_0 <= Pe_* <= Xi Delta_inf`

and Stage 40 already proved that `zeta_0^(Pe+R)` is strictly increasing in `Pe`, one gets the exact operator-side support bracket

`zeta_-(Xi,eta;kappa) <= zeta_phys(Xi,eta;kappa) <= zeta_+(Xi,eta;kappa),`

where

`zeta_-(Xi,eta;kappa)`
`= Omega_(Xi Delta_0)^2 (kappa + pi^2/4)/(kappa + y^2),`

`zeta_+(Xi,eta;kappa)`
`= Omega_(Xi Delta_inf)^2 (kappa + pi^2/4)/(kappa + y^2).`

So the full root solve is not needed to get a rigorous placement interval for the physical branch.

---

## 3. Exact residual bracket

Define the physical support residual

`R_phys(Xi,eta,kappa;zeta_req)`
`:= zeta_req - zeta_phys(Xi,eta;kappa).`

Then the branch bracket gives

`R_- <= R_phys <= R_+,`

with

`R_-(Xi,eta,kappa;zeta_req)`
`:= zeta_req - zeta_+(Xi,eta;kappa),`

`R_+(Xi,eta,kappa;zeta_req)`
`:= zeta_req - zeta_-(Xi,eta;kappa).`

This is the exact operator-selected residual envelope.

It has two immediate consequences.

### Guaranteed success

If

`zeta_-(Xi,eta;kappa) >= zeta_req,`

then the lower branch bracket already clears the target, so

`R_phys <= 0`

for every constructive branch root.

### Guaranteed failure in this operator family

If

`zeta_+(Xi,eta;kappa) < zeta_req,`

then even the upper branch bracket cannot reach the target, so

`R_phys > 0`

for every constructive branch root in this coupled operator family.

So the full fixed-point solve is only needed in the intermediate window.

---

## 4. Exact coupling thresholds `Xi_fail` and `Xi_suff`

Inside the Stage-40 reachability window, define `Pe_req` as the unique constructive solution of

`Omega_(Pe_req)^2 = zeta_req (kappa + y^2)/(kappa + pi^2/4).`

Because `Omega_Pe` is strictly increasing from `1` to `pi/2`, this `Pe_req` exists iff the Stage-40 ceiling is satisfied.

Then the exact coupling thresholds are

`Xi_fail = Pe_req / Delta_inf(kappa,eta),`

`Xi_suff = Pe_req / Delta_0(kappa,eta).`

Their meaning is exact.

### No-go threshold

If `Xi <= Xi_fail`, then

`Xi Delta_inf <= Pe_req,`

hence

`zeta_+(Xi,eta;kappa) <= zeta_req`,

so the coupled operator family cannot reach the target.

### Guaranteed-success threshold

If `Xi >= Xi_suff`, then

`Xi Delta_0 >= Pe_req,`

hence

`zeta_-(Xi,eta;kappa) >= zeta_req`,

so every constructive branch root reaches the target.

Since `Delta_inf >= Delta_0 > 0`, one has

`Xi_fail <= Xi_suff.`

So the physical theorem gap is now reduced to a bounded coupling window, not a wide-open parameter hunt.

---

## 5. Weak-coupling branch law for the residual

From Stage 39,

`Omega_Pe = 1 + ((4-pi)/(2pi)) Pe + O(Pe^2),`

so

`Omega_Pe^2 = 1 + ((4-pi)/pi) Pe + O(Pe^2).`

Combining this with the Stage-41 weak-coupling branch law

`Pe_* = Xi Delta_0 + O(Xi^2)`

gives

`zeta_phys(Xi,eta;kappa)`
`= A_K(eta;kappa)`
`  [ 1 + ((4-pi)/pi) Xi Delta_0(kappa,eta) + O(Xi^2) ],`

where

`A_K(eta;kappa) = (kappa + pi^2/4)/(kappa + y(eta)^2).`

So the operator-selected physical branch departs from the symmetric support point with a completely explicit linear slope in the coupling strength `Xi`.

---

## 6. Best current theorem statement after Stage 42

### What is exact now

- the physical support branch is no longer just the Stage-40 family; it is that family evaluated on a genuine coupled-operator branch point `Pe_*`;
- every constructive branch root obeys the exact interval

  `Xi Delta_0 <= Pe_* <= Xi Delta_inf`;
- the physical support ratio and the exact residual therefore obey rigorous lower and upper bounds;
- and the full support question has collapsed to one narrow coupling window between the exact thresholds

  `Xi_fail = Pe_req / Delta_inf,`

  `Xi_suff = Pe_req / Delta_0.`

### What this means physically

The remaining theorem gap is no longer “derive everything about the lowest support lane.”
It is now:

> compute the real moving-throat coupling strength `Xi`, compare it to the exact interval `[Xi_fail, Xi_suff]`, and then solve the fixed-point root only if the real branch lands inside that already narrow intermediate window.

=== moving_throat_pde_stage043_entropic_microclosure.md ===

# Moving-Throat PDE — Stage 43: Entropic Source Microclosure and the Microscopic Support/Source Gain

## Purpose

Stage 41 introduced the coupled support/source branch equation

`Pe = Xi Delta(Pe;kappa,eta),`

but the coupling `Xi` was still phenomenological:

`Xi = mu_sigma Lambda_phi^2 L^2 / (D_sigma T_X).`

The next honest step is to derive that quantity from the first explicit microscopic closure that is still compatible with the parent 4D ontology:

- exact projected continuity for the source density,
- a positive scalar source density carried along the throat axis,
- an explicit support/source free-energy coupling,
- and an Onsager/Fokker–Planck transport law that preserves positivity.

This stage does that.

The main results are:

1. the first explicit microscopic source/support free energy can be written as

   `F[sigma,phi]`
   `= int_0^L ds [ Theta_sigma sigma (log(sigma/sigma_*) - 1) - Lambda_phi sigma phi`
   `               + (T_X/2) phi_s^2 + (K_X/2) phi^2 ] + (K_m/2) phi(0)^2;`
2. its Euler–Lagrange variations are

   `mu_sigma^(chem) = delta F / delta sigma = Theta_sigma log(sigma/sigma_*) - Lambda_phi phi,`

   `-T_X phi_ss + K_X phi = Lambda_phi sigma,`

   with Robin/Neumann support boundary conditions

   `T_X phi_s(0)=K_m phi(0),` `phi_s(L)=0;`
3. the minimal positive-density Onsager current is

   `J = -M_sigma sigma partial_s mu_sigma^(chem)`
   `  = -D_sigma partial_s sigma + M_sigma Lambda_phi sigma partial_s phi,`

   where the exact Einstein relation is

   `D_sigma = M_sigma Theta_sigma;`
4. under the same affine-drop reduction already implicit in the Stage-41 average-drift closure,

   `phi(s) ~= phi(0) + [Delta phi] s/L,`

   the stationary zero-flux branch is exactly the Stage-39 exponential family,

   `sigma(s) = C exp[(Lambda_phi Delta phi)/(Theta_sigma L) s],`

   or, in normalized coordinates `x=s/L`,

   `Sigma_Pe(x) = Pe exp(Pe x)/(exp(Pe)-1);`
5. the transport bias is therefore no longer phenomenological:

   `Pe = (Lambda_phi/Theta_sigma) Delta phi;`
6. using the Stage-41 support normalization

   `Phi = T_X phi/(Lambda_phi L^2),`
   `Delta phi = (Lambda_phi L^2 / T_X) Delta(Pe;kappa,eta),`

   the branch equation becomes

   `Pe = Xi_micro Delta(Pe;kappa,eta),`

   with the exact microscopic coupling

   `Xi_micro = Lambda_phi^2 L^2 / (Theta_sigma T_X)`
   `        = chi_sigma Lambda_phi^2 L^2 / T_X,`

   where `chi_sigma := 1/Theta_sigma;`
7. equivalently, the phenomenological Stage-41 coupling is now explained as

   `Xi = mu_sigma Lambda_phi^2 L^2/(D_sigma T_X)`
   `  = Lambda_phi^2 L^2/(Theta_sigma T_X)`

   by the Einstein relation `D_sigma = mu_sigma Theta_sigma;`
8. the source dynamics is passive: if `phi` is held fixed or slaved quasi-statically by the support Euler–Lagrange equation, then

   `dF/dt = - int_0^L ds J^2/(M_sigma sigma) <= 0`

   under no-flux boundaries.

So the operator strength is no longer an arbitrary ratio of mobility and diffusion. On the first explicit microscopic closure, it collapses to one entropic-support gain

`Xi_micro = chi_sigma Lambda_phi^2 L^2/T_X.`

---

## 1. Microscopic source/support free energy

Take the throat-axis source density `sigma(s,t) > 0` on `s in [0,L]` and the support field `phi(s,t)`.

The first explicit free-energy functional that preserves positivity of `sigma` and reproduces the Stage-41 support operator is

`F[sigma,phi]`
`= int_0^L ds [ Theta_sigma sigma (log(sigma/sigma_*) - 1) - Lambda_phi sigma phi`
`               + (T_X/2) phi_s^2 + (K_X/2) phi^2 ] + (K_m/2) phi(0)^2.`

Here:

- `Theta_sigma` is the source entropic/chemical scale,
- `sigma_*` is a reference density,
- `Lambda_phi` is the support/source coupling,
- `T_X` is the axial support tension,
- `K_X` is the baseline support stiffness,
- `K_m` is the Robin mouth spring.

Variation with respect to `sigma` gives the exact chemical potential

`mu_sigma^(chem) = Theta_sigma log(sigma/sigma_*) - Lambda_phi phi.`

Variation with respect to `phi` gives the support equation

`-T_X phi_ss + K_X phi = Lambda_phi sigma,`

with boundary conditions

`T_X phi_s(0)=K_m phi(0),`

`phi_s(L)=0.`

So the Stage-41 support law is now embedded directly into one explicit free-energy closure.

---

## 2. Positive-density Onsager transport law

Use the minimal gradient-flow current

`J = -M_sigma sigma partial_s mu_sigma^(chem),`

with mobility `M_sigma > 0`.

Expanding the chemical potential gradient gives

`J = -M_sigma sigma [ Theta_sigma partial_s log sigma - Lambda_phi partial_s phi ]`
`  = -M_sigma Theta_sigma partial_s sigma + M_sigma Lambda_phi sigma partial_s phi.`

So the transport law takes the drift-diffusion form

`J = -D_sigma partial_s sigma + M_sigma Lambda_phi sigma partial_s phi,`

with exact Einstein relation

`D_sigma = M_sigma Theta_sigma.`

This is the microscopic origin of the Stage-41 phenomenological pair `(D_sigma,mu_sigma)`.

---

## 3. Recovery of the Stage-39 exponential family

Stage 41 already used an average-drift closure in which the support field enters only through the end-to-end drop.

Within that same lowest-lane closure, replace the support profile by its affine-drop reduction

`phi(s) ~= phi(0) + [Delta phi] s/L,`

so that `partial_s phi ~= Delta phi / L` is constant.

On the stationary zero-flux branch `J=0`, the exact transport law gives

`partial_s sigma = (Lambda_phi Delta phi)/(Theta_sigma L) sigma.`

Therefore

`sigma(s) = C exp[(Lambda_phi Delta phi)/(Theta_sigma L) s].`

Normalize on `[0,L]` and pass to `x=s/L in [0,1]`. Then

`Sigma_Pe(x) = L sigma(Lx)`
`           = Pe exp(Pe x)/(exp(Pe)-1),`

with the exact Peclet number

`Pe = (Lambda_phi/Theta_sigma) Delta phi.`

So the Stage-39 exponential branch is not just a convenient guess. It is the exact stationary branch of the first positive-density Onsager closure after the same affine-drop reduction already used implicitly in Stage 41.

---

## 4. Exact microscopic coupling `Xi_micro`

Stage 41 defined the dimensionless support field by

`Phi = T_X phi/(Lambda_phi L^2),`

so the end-to-end physical support drop is

`Delta phi = (Lambda_phi L^2/T_X) Delta(Pe;kappa,eta).`

Insert this into the exact Peclet law above:

`Pe = (Lambda_phi/Theta_sigma) Delta phi`
`   = (Lambda_phi^2 L^2/(Theta_sigma T_X)) Delta(Pe;kappa,eta).`

Therefore the Stage-41 branch equation is recovered with the exact microscopic coupling

`Xi_micro = Lambda_phi^2 L^2/(Theta_sigma T_X).`

Equivalently, with the source susceptibility

`chi_sigma := 1/Theta_sigma,`

one has

`Xi_micro = chi_sigma Lambda_phi^2 L^2/T_X.`

And because `D_sigma = M_sigma Theta_sigma`, this also reproduces the Stage-41 phenomenological form:

`Xi_micro = mu_sigma Lambda_phi^2 L^2/(D_sigma T_X).`

So the apparent mobility/diffusion ambiguity is gone. The coupling is one entropic-support gain.

---

## 5. Exact passivity / Lyapunov identity

The same free-energy closure supplies an exact dissipation law.

At fixed `phi` one has

`dF/dt = int_0^L ds mu_sigma^(chem) partial_t sigma.`

Using continuity `partial_t sigma = -partial_s J` gives

`dF/dt = - int_0^L ds mu_sigma^(chem) partial_s J`
`     = - [ mu_sigma^(chem) J ]_0^L + int_0^L ds (partial_s mu_sigma^(chem)) J.`

Now insert the Onsager current

`J = -M_sigma sigma partial_s mu_sigma^(chem).`

Then

`(partial_s mu_sigma^(chem)) J = - J^2/(M_sigma sigma).`

So the exact free-energy identity is

`dF/dt = - [ mu_sigma^(chem) J ]_0^L - int_0^L ds J^2/(M_sigma sigma).`

Under no-flux boundaries `J(0)=J(L)=0`,

`dF/dt = - int_0^L ds J^2/(M_sigma sigma) <= 0.`

If `phi` is slaved quasi-statically by the support Euler–Lagrange equation, the same identity applies to the full coupled free energy because the `phi`-variation term vanishes on shell.

So this microscopic closure is automatically passive.

---

## 6. What Stage 43 changes

The operator problem has now advanced in a concrete way.

Before Stage 43, the support/source strength was the phenomenological ratio

`Xi = mu_sigma Lambda_phi^2 L^2/(D_sigma T_X).`

After Stage 43, the same quantity is an explicit microscopic gain:

`Xi_micro = chi_sigma Lambda_phi^2 L^2/T_X.`

That is a much tighter result because:

- it is independent of the separate values of mobility and diffusion,
- it is tied to one explicit free-energy closure,
- it preserves positivity of the source density,
- it reproduces the Stage-39 exponential family exactly in the same lowest-lane reduction already used earlier,
- and it is automatically passive.

So the remaining theorem gap is now sharper again:

> compute the actual moving-throat values of `chi_sigma`, `Lambda_phi`, `T_X`, and `L`, form `Xi_micro`, and compare it directly to the exact support thresholds from Stage 42.

=== moving_throat_pde_stage044_microscopic_gain_thresholds.md ===

# Moving-Throat PDE — Stage 44: Microscopic Gain Thresholds and the Exact Operator Phase Diagram

## Purpose

Stage 43 replaced the phenomenological support/source strength by the explicit microscopic gain

`Xi_micro = chi_sigma Lambda_phi^2 L^2 / T_X.`

The next honest step is to compare that gain directly with the exact Stage-42 thresholds

`Xi_fail = Pe_req / Delta_inf(kappa,eta),`

`Xi_suff = Pe_req / Delta_0(kappa,eta).`

This stage does that and pushes the result one step further by isolating a cleaner dimensionless microscopic control parameter.

The main results are:

1. using

   `kappa = K_X L^2 / T_X,`

   the microscopic coupling can be rewritten as

   `Xi_micro = kappa G_micro,`

   where the exact dimensionless support/source gain is

   `G_micro := chi_sigma Lambda_phi^2 / K_X;`
2. the exact no-go / sufficiency thresholds therefore become

   `G_fail = Pe_req / [ kappa Delta_inf(kappa,eta) ],`

   `G_suff = Pe_req / [ kappa Delta_0(kappa,eta) ],`

   with the exact phase diagram

   - fail if `G_micro <= G_fail`,
   - succeed if `G_micro >= G_suff`,
   - only the interval `G_fail < G_micro < G_suff` requires the full root solve;
3. equivalently, the original microscopic parameters satisfy the exact threshold surfaces

   `chi_sigma <= T_X Pe_req / [ Lambda_phi^2 L^2 Delta_inf ]`  -> fail,

   `chi_sigma >= T_X Pe_req / [ Lambda_phi^2 L^2 Delta_0 ]`    -> succeed,

   and, at fixed `chi_sigma`,

   `Lambda_phi^2 <= T_X Pe_req / [ chi_sigma L^2 Delta_inf ]`  -> fail,

   `Lambda_phi^2 >= T_X Pe_req / [ chi_sigma L^2 Delta_0 ]`    -> succeed;
4. in the soft-support limit `kappa -> 0^+`, the exact endpoint values are

   `Delta_0 -> 1/2,`

   `Delta_inf -> 1,`

   so the microscopic gain thresholds diverge as

   `G_fail ~ Pe_req / kappa,`

   `G_suff ~ 2 Pe_req / kappa;`
5. in the highly compliant-mouth limit `eta -> +infty`, the endpoint data simplify exactly to

   `Delta_0^(inf) = (1-sech(sqrt(kappa)))/kappa,`

   `Delta_inf^(inf) = tanh(sqrt(kappa))/sqrt(kappa),`

   and therefore

   `G_fail^(inf) = Pe_req / [ sqrt(kappa) tanh(sqrt(kappa)) ],`

   `G_suff^(inf) = Pe_req / [ 1 - sech(sqrt(kappa)) ];`
6. in the combined stiff-support / compliant-mouth regime `kappa >> 1`, `eta >> 1`, these reduce to

   `G_fail ~ Pe_req / sqrt(kappa),`

   `G_suff ~ Pe_req.`

So the explicit microscopic closure now has a very sharp interpretation:

- too-soft support (`kappa << 1`) requires unrealistically large gain,
- a sufficiently compliant mouth makes the sufficiency threshold saturate at an `O(Pe_req)` microscopic gain,
- and only the bounded intermediate band between `G_fail` and `G_suff` still requires the full fixed-point solve.

---

## 1. Exact microscopic gain `G_micro`

Stage 43 gave

`Xi_micro = chi_sigma Lambda_phi^2 L^2 / T_X.`

But Stage 40 already uses the support stiffness ratio

`kappa = K_X L^2 / T_X.`

Therefore

`Xi_micro = kappa ( chi_sigma Lambda_phi^2 / K_X ).`

Define the exact dimensionless support/source gain

`G_micro := chi_sigma Lambda_phi^2 / K_X.`

Then

`Xi_micro = kappa G_micro.`

This is the cleanest microscopic control parameter so far because it removes the explicit length and tension scales from the branch-strength comparison.

---

## 2. Exact operator phase diagram in microscopic variables

Stage 42 already proved that the operator-selected branch succeeds or fails according to the comparison of `Xi` with the exact thresholds

`Xi_fail = Pe_req / Delta_inf(kappa,eta),`

`Xi_suff = Pe_req / Delta_0(kappa,eta).`

Substitute `Xi_micro = kappa G_micro`. Then the exact microscopic thresholds are

`G_fail(kappa,eta) = Pe_req / [ kappa Delta_inf(kappa,eta) ],`

`G_suff(kappa,eta) = Pe_req / [ kappa Delta_0(kappa,eta) ].`

So the exact phase diagram is now

- `G_micro <= G_fail`  -> impossible inside this operator family,
- `G_micro >= G_suff`  -> guaranteed success,
- `G_fail < G_micro < G_suff` -> only then is the full root solve needed.

This is the first operator phase diagram written directly in microscopic support/source variables.

---

## 3. Threshold surfaces for `chi_sigma` and `Lambda_phi`

Undoing the definition of `G_micro` gives exact threshold surfaces in the original microscopic variables.

At fixed `Lambda_phi`, success and failure are controlled by

`chi_sigma <= T_X Pe_req / [ Lambda_phi^2 L^2 Delta_inf(kappa,eta) ]`  -> fail,

`chi_sigma >= T_X Pe_req / [ Lambda_phi^2 L^2 Delta_0(kappa,eta) ]`    -> succeed.

At fixed `chi_sigma`, they are controlled by

`Lambda_phi^2 <= T_X Pe_req / [ chi_sigma L^2 Delta_inf(kappa,eta) ]`  -> fail,

`Lambda_phi^2 >= T_X Pe_req / [ chi_sigma L^2 Delta_0(kappa,eta) ]`    -> succeed.

So the microscopic theorem gap is no longer “somehow get a large enough `Xi`.” It is now a concrete competition between source susceptibility and support/source coupling on one side, and the exact support geometry functions `Delta_0`, `Delta_inf` on the other.

---

## 4. Exact soft-support limit `kappa -> 0^+`

Using the exact endpoint formulas from Stage 41,

`Delta_0(kappa,eta)`
`= eta (cosh(sqrt(kappa))-1)`
`  / [ kappa ( sqrt(kappa) sinh(sqrt(kappa)) + eta cosh(sqrt(kappa)) ) ],`

`Delta_inf(kappa,eta)`
`= [ cosh(sqrt(kappa)) + (eta/sqrt(kappa)) sinh(sqrt(kappa)) - 1 ]`
`  / [ sqrt(kappa) sinh(sqrt(kappa)) + eta cosh(sqrt(kappa)) ],`

the exact soft-support limits are

`lim_(kappa->0+) Delta_0 = 1/2,`

`lim_(kappa->0+) Delta_inf = 1.`

Therefore

`G_fail ~ Pe_req / kappa,`

`G_suff ~ 2 Pe_req / kappa`

as `kappa -> 0^+`.

So a very soft baseline support channel is strongly disfavored: the microscopic gain required for success diverges like `1/kappa`.

---

## 5. Exact highly compliant-mouth limit `eta -> +infty`

At fixed `kappa`, the mouth-compliant limit is exact.

From the same endpoint formulas,

`lim_(eta->+infty) Delta_0(kappa,eta) = (1-sech(sqrt(kappa)))/kappa,`

`lim_(eta->+infty) Delta_inf(kappa,eta) = tanh(sqrt(kappa))/sqrt(kappa).`

So the exact microscopic gain thresholds become

`G_fail^(inf)(kappa) = Pe_req / [ sqrt(kappa) tanh(sqrt(kappa)) ],`

`G_suff^(inf)(kappa) = Pe_req / [ 1 - sech(sqrt(kappa)) ].`

This is useful because it isolates the best-case mouth-compliance regime of the same operator family.

Two consequences are immediate.

### Small `kappa` inside the compliant-mouth branch

If `kappa << 1`, then

`tanh(sqrt(kappa)) ~ sqrt(kappa),`

`1 - sech(sqrt(kappa)) ~ kappa/2,`

so the exact soft-support divergence is recovered:

`G_fail^(inf) ~ Pe_req / kappa,`

`G_suff^(inf) ~ 2 Pe_req / kappa.`

### Stiff-support inside the compliant-mouth branch

If `kappa >> 1`, then

`tanh(sqrt(kappa)) -> 1,`

`sech(sqrt(kappa)) -> 0,`

so

`G_fail^(inf) ~ Pe_req / sqrt(kappa),`

`G_suff^(inf) ~ Pe_req.`

That means sufficiently strong baseline support stiffness is not itself the main obstacle once the mouth is compliant enough; the sufficiency threshold saturates at an `O(Pe_req)` microscopic gain.

---

## 6. What Stage 44 changes

The theorem gap has narrowed again.

Before Stage 44, the support question was:

> what is the physical `Xi`, and is it above or below the exact interval `[Xi_fail, Xi_suff]`?

After Stage 44, it is:

> what is the physical dimensionless gain `G_micro = chi_sigma Lambda_phi^2/K_X`, and does it lie above or below the exact geometry-controlled thresholds `G_fail(kappa,eta)` and `G_suff(kappa,eta)`?

That is a stronger result because it separates the operator problem into two transparent pieces:

- the support geometry sector `(kappa,eta)`, which sets the exact threshold functions,
- and the microscopic source/support gain `G_micro`, which the completed moving-throat operator must supply.

The sharpest practical lesson is now clear:

- very soft support (`kappa << 1`) is strongly unfavorable,
- a highly compliant mouth (`eta >> 1`) helps substantially,
- and on the best-case compliant branch, success is controlled by whether `G_micro` can reach an order-`Pe_req` value.

So the next honest step is to compute `G_micro = chi_sigma Lambda_phi^2/K_X` from a more explicit moving-throat branch and compare it directly to the exact threshold surfaces derived here.

=== moving_throat_pde_stage045_parent_action_gain.md ===

# Moving-Throat PDE — Stage 45: Parent-Action Projection of the Microscopic Support/Source Gain

## Purpose

Stage 44 reduced the support/source theorem gap to one microscopic gain,

`G_micro = chi_sigma Lambda_phi^2 / K_X,`

but that quantity was still written in the language of the reduced axial source/support operator.

The next honest step is therefore to derive that gain from the **parent 4D action** rather than leave it as a phenomenological effective constant.

This stage does that by projecting the linearized GNLS/confinement sector of the parent theory onto one compressional source channel and one throat-support channel.

The main results are:

1. starting from the parent matter energy

   `H_psi = int d^4X [ (hbar^2/2m)|D_i psi|^2 + V_conf rho + U(rho) ],`

   with the frozen `n=5` EOS

   `U(rho) = K rho^5 / 4,`

   `h(rho) = dU/drho = (5K/4) rho^4,`

   the exact local compressional stiffness is

   `h'(rho_*) = 5 K rho_*^3 = m c_(s*)^2 / rho_*;`
2. after linearizing about a static throat branch `rho = rho_* + delta rho`, the compressional part of the matter energy is

   `(1/2) h'(rho_*) (delta rho)^2;`
3. if the support field enters the confinement as the first linear throat-shape perturbation

   `delta V_conf(s,y) = - g_phi chi_phi(y) phi(s),`

   and the source density is truncated to one transverse compression channel

   `delta rho(s,y) = sigma(s) chi_sigma(y),`

   then the parent 4D energy reduces exactly to the Stage-43 form

   `F_red[sigma,phi]`
   `= int_0^L ds [ (Theta_sigma/2) sigma^2 - Lambda_phi sigma phi`
   `               + (T_X/2) phi_s^2 + (K_X/2) phi^2 ] + (K_m/2) phi(0)^2,`

   with

   `Theta_sigma = h'(rho_*) N_(sigma sigma),`

   `Lambda_phi = g_phi O_(sigma phi);`
4. here the parent overlap invariants are

   `N_(sigma sigma) = int d^3y chi_sigma(y)^2,`

   `N_(phi phi)     = int d^3y chi_phi(y)^2,`

   `O_(sigma phi)   = int d^3y chi_sigma(y) chi_phi(y);`
5. the effective source susceptibility is therefore

   `chi_sigma^(eff) = 1 / Theta_sigma`
   `               = rho_* / [ m c_(s*)^2 N_(sigma sigma) ];`
6. and the microscopic gain becomes the explicit parent-action quantity

   `G_micro`
   `= chi_sigma^(eff) Lambda_phi^2 / K_X`
   `= rho_* g_phi^2 O_(sigma phi)^2 / [ m c_(s*)^2 K_X N_(sigma sigma) ];`
7. introducing the source/support coherence factor

   `C_(sigma phi)^2 := O_(sigma phi)^2 / [ N_(sigma sigma) N_(phi phi) ],`

   one obtains the exact factorization

   `G_micro = [ rho_* g_phi^2 N_(phi phi) / (m c_(s*)^2 K_X) ] C_(sigma phi)^2,`

   with `0 <= C_(sigma phi)^2 <= 1` by Cauchy–Schwarz;
8. finally, inserting `kappa = K_X L^2 / T_X` gives the microscopic fixed-point strength

   `Xi_micro = kappa G_micro`
   `        = rho_* g_phi^2 O_(sigma phi)^2 L^2 / [ m c_(s*)^2 T_X N_(sigma sigma) ].`

So the Stage-44 gain is no longer a free microscopic placeholder. It is the compressional susceptibility of the `n=5` GNLS medium times the square of one parent confinement/support overlap, divided by the baseline support stiffness.

---

## 1. Parent 4D matter energy and the `n=5` compressional stiffness

The parent 4D theory already fixes the matter sector to the gauged GNLS form with confinement and the frozen stiff EOS,

`H_psi = int d^4X [ (hbar^2/2m)|D_i psi|^2 + V_conf(X;a,L) rho + U(rho) ],`

`U(rho) = K rho^5 / 4,`

`h(rho) = dU/drho = (5K/4) rho^4,`

`c_s^2(rho) = (1/m) dP/drho = (5K/m) rho^4.`

Therefore

`h'(rho) = 5K rho^3 = m c_s^2(rho) / rho.`

Linearizing about a static throat branch `rho = rho_* + delta rho` and subtracting the equilibrium linear term leaves the local compressional quadratic energy

`delta H_comp = (1/2) int d^4X h'(rho_*) (delta rho)^2.`

So the parent matter sector already provides the exact local scalar stiffness needed by the Stage-43 source entropy term.

---

## 2. Parent-action reduction to one source channel and one support channel

Let the support field be the first axial throat-support amplitude `phi(s)` and let its leading linear effect on the confinement be

`delta V_conf(s,y) = - g_phi chi_phi(y) phi(s),`

where:

- `s in [0,L]` is the throat axis,
- `y` denotes the transverse/cross-sectional coordinates,
- `chi_phi(y)` is the transverse support profile,
- `g_phi` is the parent confinement-loading amplitude.

Now truncate the density perturbation to one transverse compression channel,

`delta rho(s,y) = sigma(s) chi_sigma(y),`

with source profile `chi_sigma(y)`.

Insert this into the parent linearized energy. The compressional term becomes

`(1/2) h'(rho_*) int_0^L ds sigma(s)^2 int d^3y chi_sigma(y)^2,`

while the linear support/source coupling becomes

`- int_0^L ds sigma(s) phi(s) g_phi int d^3y chi_sigma(y) chi_phi(y).`

Define the overlap invariants

`N_(sigma sigma) = int d^3y chi_sigma^2,`

`N_(phi phi)     = int d^3y chi_phi^2,`

`O_(sigma phi)   = int d^3y chi_sigma chi_phi.`

Then the parent energy reduces exactly to

`F_red[sigma,phi]`
`= int_0^L ds [ (Theta_sigma/2) sigma^2 - Lambda_phi sigma phi`
`               + (T_X/2) phi_s^2 + (K_X/2) phi^2 ] + (K_m/2) phi(0)^2,`

with

`Theta_sigma = h'(rho_*) N_(sigma sigma),`

`Lambda_phi = g_phi O_(sigma phi).`

This is precisely the Stage-43 reduced support/source free energy, now derived as a one-channel projection of the parent 4D action.

---

## 3. Exact parent formula for the effective source susceptibility

The reduced source susceptibility is

`chi_sigma^(eff) = 1 / Theta_sigma`
`               = 1 / [ h'(rho_*) N_(sigma sigma) ].`

Using the exact `n=5` identity above,

`h'(rho_*) = m c_(s*)^2 / rho_*,`

one gets

`chi_sigma^(eff) = rho_* / [ m c_(s*)^2 N_(sigma sigma) ].`

So the source compliance is no longer a free “entropy constant.” It is fixed by the local GNLS compressibility of the parent medium, dressed only by the chosen transverse source normalization.

---

## 4. Exact parent formula for the microscopic gain

Stage 44 defined the microscopic gain by

`G_micro = chi_sigma Lambda_phi^2 / K_X.`

Insert the projected parent coefficients:

`G_micro`
`= [ 1 / ( h'(rho_*) N_(sigma sigma) ) ] [ g_phi^2 O_(sigma phi)^2 ] / K_X`
`= rho_* g_phi^2 O_(sigma phi)^2 / [ m c_(s*)^2 K_X N_(sigma sigma) ].`

This is the first explicit parent-action formula for the gain that controls the operator phase diagram.

It says that the support/source drive is large when:

- the local medium is easily compressible (`rho_* / m c_(s*)^2` large),
- the confinement/support loading amplitude `g_phi` is large,
- the source and support transverse channels overlap strongly,
- and the baseline support stiffness `K_X` is small.

---

## 5. Coherence factor and the exact overlap decomposition

Define the source/support coherence factor

`C_(sigma phi)^2 := O_(sigma phi)^2 / [ N_(sigma sigma) N_(phi phi) ].`

Then

`0 <= C_(sigma phi)^2 <= 1`

by Cauchy–Schwarz, and the microscopic gain factorizes as

`G_micro`
`= [ rho_* g_phi^2 N_(phi phi) / (m c_(s*)^2 K_X) ] C_(sigma phi)^2.`

This is useful because it separates the microscopic problem into two independent pieces:

1. a **strength scale**

   `rho_* g_phi^2 N_(phi phi) / (m c_(s*)^2 K_X),`

2. a **profile-coherence factor**

   `C_(sigma phi)^2 in [0,1].`

So the unresolved PDE-side question is no longer a single opaque parameter. It is the combination of a parent confinement-loading strength and the coherence with which the source and support channels line up on the true branch.

---

## 6. Microscopic fixed-point strength `Xi_micro`

Stage 43 used

`Xi_micro = kappa G_micro,`

with

`kappa = K_X L^2 / T_X.`

Therefore the parent projected formula becomes

`Xi_micro`
`= rho_* g_phi^2 O_(sigma phi)^2 L^2 / [ m c_(s*)^2 T_X N_(sigma sigma) ].`

So the full operator-selected branch point

`Pe = Xi_micro Delta(Pe;kappa,eta)`

is now determined by three parent-action ingredients only:

- the local compressional stiffness of the `n=5` GNLS medium,
- the confinement/support loading overlap,
- and the support tension/Robin geometry sector entering `(T_X,kappa,eta)`.

That is the cleanest microscopic restatement yet of the support/source theorem gap.

=== moving_throat_pde_stage046_parent_thresholds.md ===

# Moving-Throat PDE — Stage 46: Parent-Overlap Threshold Theorem and Exact Microscopic Success/No-Go Tests

## Purpose

Stage 45 expressed the Stage-44 microscopic gain directly in parent 4D variables:

`G_micro = rho_* g_phi^2 O_(sigma phi)^2 / [ m c_(s*)^2 K_X N_(sigma sigma) ]`

or, equivalently,

`G_micro = [ rho_* g_phi^2 N_(phi phi) / (m c_(s*)^2 K_X) ] C_(sigma phi)^2.`

The next honest step is to combine that exact parent formula with the Stage-44 operator phase diagram and turn the support/source theorem gap into explicit thresholds on:

- the parent confinement-loading amplitude `g_phi`,
- and the source/support coherence factor `C_(sigma phi)^2`.

This stage does that.

The main results are:

1. inserting the parent gain into the Stage-44 thresholds gives the exact parent fail/succeed conditions

   `g_phi^2 <= g_(phi,fail)^2`  -> fail,

   `g_phi^2 >= g_(phi,suff)^2`  -> succeed,

   with

   `g_(phi,fail)^2`
   `= m c_(s*)^2 K_X N_(sigma sigma) G_fail / [ rho_* O_(sigma phi)^2 ],`

   `g_(phi,suff)^2`
   `= m c_(s*)^2 K_X N_(sigma sigma) G_suff / [ rho_* O_(sigma phi)^2 ];`
2. equivalently, in coherence-factor form,

   `g_(phi,fail)^2`
   `= m c_(s*)^2 K_X G_fail / [ rho_* N_(phi phi) C_(sigma phi)^2 ],`

   `g_(phi,suff)^2`
   `= m c_(s*)^2 K_X G_suff / [ rho_* N_(phi phi) C_(sigma phi)^2 ];`
3. solving instead for the required source/support alignment gives the exact coherence thresholds

   `C_fail^2 = m c_(s*)^2 K_X G_fail / [ rho_* g_phi^2 N_(phi phi) ],`

   `C_suff^2 = m c_(s*)^2 K_X G_suff / [ rho_* g_phi^2 N_(phi phi) ],`

   with

   `C_fail^2 <= C_suff^2;`
4. because `0 <= C_(sigma phi)^2 <= 1`, there is an immediate Cauchy-based no-go theorem:

   if

   `rho_* g_phi^2 N_(phi phi) / (m c_(s*)^2 K_X) < G_fail,`

   then the branch fails **for every possible** source profile in the chosen support channel;
5. the exact “best-case” reachable gain at fixed `g_phi` is therefore

   `G_max(g_phi) = rho_* g_phi^2 N_(phi phi) / (m c_(s*)^2 K_X),`

   attained only at perfect alignment `C_(sigma phi)^2 = 1`;
6. inserting `kappa = K_X L^2 / T_X` and the Stage-44 formulas

   `G_fail = Pe_req / [ kappa Delta_inf(kappa,eta) ],`

   `G_suff = Pe_req / [ kappa Delta_0(kappa,eta) ],`

   the parent amplitude thresholds become

   `g_(phi,fail)^2`
   `= m c_(s*)^2 T_X N_(sigma sigma) Pe_req / [ rho_* L^2 O_(sigma phi)^2 Delta_inf(kappa,eta) ],`

   `g_(phi,suff)^2`
   `= m c_(s*)^2 T_X N_(sigma sigma) Pe_req / [ rho_* L^2 O_(sigma phi)^2 Delta_0(kappa,eta) ];`
7. in coherence-factor form this is

   `g_(phi,fail)^2`
   `= m c_(s*)^2 T_X Pe_req / [ rho_* L^2 N_(phi phi) C_(sigma phi)^2 Delta_inf ],`

   `g_(phi,suff)^2`
   `= m c_(s*)^2 T_X Pe_req / [ rho_* L^2 N_(phi phi) C_(sigma phi)^2 Delta_0 ],`

   so the prefactor `K_X` cancels exactly, surviving only through the shape function `Delta_(0,inf)(kappa,eta)`.

So the remaining PDE-side theorem gap is no longer “find some support/source gain.” It is now the exact comparison of one parent confinement-loading amplitude and one source/support coherence factor against the explicit fail/succeed thresholds above.

---

## 1. Exact parent thresholds on the confinement-loading amplitude `g_phi`

Stage 44 proved the microscopic operator phase diagram

`G_micro <= G_fail(kappa,eta)`  -> fail,

`G_micro >= G_suff(kappa,eta)`  -> succeed.

Insert the parent formula from Stage 45,

`G_micro = rho_* g_phi^2 O_(sigma phi)^2 / [ m c_(s*)^2 K_X N_(sigma sigma) ].`

Then the fail/succeed thresholds on the parent loading amplitude are exact:

`g_(phi,fail)^2`
`= m c_(s*)^2 K_X N_(sigma sigma) G_fail / [ rho_* O_(sigma phi)^2 ],`

`g_(phi,suff)^2`
`= m c_(s*)^2 K_X N_(sigma sigma) G_suff / [ rho_* O_(sigma phi)^2 ].`

Therefore

- if `g_phi^2 <= g_(phi,fail)^2`, the chosen parent branch cannot reach the target inside this operator family,
- if `g_phi^2 >= g_(phi,suff)^2`, the branch is guaranteed to reach it,
- only the interval `g_(phi,fail)^2 < g_phi^2 < g_(phi,suff)^2` still requires the full fixed-point solve.

---

## 2. Exact threshold on the profile coherence factor

Using

`O_(sigma phi)^2 = N_(sigma sigma) N_(phi phi) C_(sigma phi)^2,`

the gain can be rewritten as

`G_micro = [ rho_* g_phi^2 N_(phi phi) / (m c_(s*)^2 K_X) ] C_(sigma phi)^2.`

So the exact required coherence floors are

`C_fail^2 = m c_(s*)^2 K_X G_fail / [ rho_* g_phi^2 N_(phi phi) ],`

`C_suff^2 = m c_(s*)^2 K_X G_suff / [ rho_* g_phi^2 N_(phi phi) ].`

These control the same three-zone structure:

- if `C_(sigma phi)^2 <= C_fail^2`, fail;
- if `C_(sigma phi)^2 >= C_suff^2`, succeed;
- only the narrow intermediate interval still needs the full root solve.

Because `G_fail <= G_suff`, one has

`C_fail^2 <= C_suff^2.`

So the unresolved PDE data have been split very cleanly: one parent strength scale and one profile-alignment factor.

---

## 3. Exact Cauchy no-go theorem

The coherence factor satisfies

`0 <= C_(sigma phi)^2 <= 1`.

So at fixed parent loading amplitude `g_phi`, the **largest possible** gain in the chosen support channel is

`G_max(g_phi) = rho_* g_phi^2 N_(phi phi) / (m c_(s*)^2 K_X),`

obtained only at perfect alignment `C_(sigma phi)^2 = 1`.

Therefore there is an exact no-go theorem:

if

`G_max(g_phi) < G_fail(kappa,eta),`

then the branch fails for **every possible** source profile in that support channel.

Equivalently,

if

`rho_* g_phi^2 N_(phi phi) / (m c_(s*)^2 K_X) < G_fail(kappa,eta),`

then no profile engineering of `chi_sigma` can rescue the branch.

This is the first exact parent-overlap no-go theorem in the program.

---

## 4. Exact amplitude thresholds in terms of `Pe_req`, `Delta_0`, and `Delta_inf`

Stage 44 already gave

`G_fail = Pe_req / [ kappa Delta_inf(kappa,eta) ],`

`G_suff = Pe_req / [ kappa Delta_0(kappa,eta) ],`

with

`kappa = K_X L^2 / T_X.`

Insert these into the parent amplitude thresholds.

For failure:

`g_(phi,fail)^2`
`= [ m c_(s*)^2 K_X N_(sigma sigma) / (rho_* O_(sigma phi)^2) ]`
`  [ Pe_req / (kappa Delta_inf) ]`
`= m c_(s*)^2 T_X N_(sigma sigma) Pe_req`
`  / [ rho_* L^2 O_(sigma phi)^2 Delta_inf(kappa,eta) ].`

For sufficiency:

`g_(phi,suff)^2`
`= m c_(s*)^2 T_X N_(sigma sigma) Pe_req`
`  / [ rho_* L^2 O_(sigma phi)^2 Delta_0(kappa,eta) ].`

So the explicit prefactor `K_X` cancels exactly, surviving only inside the geometry shape functions `Delta_inf` and `Delta_0` through `kappa`.

This is a very clean structural result: the amplitude threshold is tension/length/compressibility controlled, while the baseline support stiffness enters only through the support-shape response.

---

## 5. Matched-profile and normalized-lane special case

If the source and support channels coincide and are normalized,

`chi_sigma = chi_phi,`

`N_(sigma sigma)=N_(phi phi)=1,`

`O_(sigma phi)=1,`

`C_(sigma phi)^2 = 1,`

then the parent thresholds simplify to

`g_(phi,fail)^2 = m c_(s*)^2 K_X G_fail / rho_*,`

`g_(phi,suff)^2 = m c_(s*)^2 K_X G_suff / rho_*.`

Equivalently,

`g_(phi,fail)^2 = m c_(s*)^2 T_X Pe_req / [ rho_* L^2 Delta_inf(kappa,eta) ],`

`g_(phi,suff)^2 = m c_(s*)^2 T_X Pe_req / [ rho_* L^2 Delta_0(kappa,eta) ].`

So in the best possible aligned lowest-lane closure, the theorem gap reduces to a single confinement-loading amplitude compared against two explicit geometry-controlled thresholds.

---

## 6. What Stage 46 changes

After Stage 46, the remaining theorem gap is no longer phrased in terms of the abstract microscopic gain `G_micro`.

It is now:

1. compute the parent confinement-loading amplitude `g_phi`,
2. compute the profile coherence `C_(sigma phi)^2` on the true moving-throat branch,
3. compare them against the exact thresholds above.

That is stronger than the Stage-44 statement because it pushes the support/source theorem gap all the way back to parent-action overlap data.

In other words, the unresolved PDE question is now no longer “is the microscopic gain big enough?” It is

> does the completed moving-throat branch generate sufficient **parent confinement loading** and sufficient **source/support coherence** to cross the explicit fail/succeed surfaces above?

=== moving_throat_pde_stage047_equilibrium_alignment.md ===

# Moving-Throat PDE — Stage 47: Parent Equilibrium Source/Support Alignment and the Exact Matched-Layer Gain

## Purpose

Stages 45–46 reduced the support/source theorem gap to two parent quantities:

- the confinement-loading amplitude `g_phi`,
- and the source/support coherence `C_(sigma phi)^2`.

But both were still being treated as independent branch data.

The next honest step is to ask whether the **parent equilibrium equations themselves** tie the source channel to the support channel.

This stage shows that they do.

Within the same local compressional linearization already used in Stages 43–46, a quasi-static wall/support displacement does not excite an arbitrary source profile. It excites a very specific one:

`chi_sigma(y) = g_phi chi_phi(y) / H(y),`

where

`H(y) := h'(rho_*(y))`

is the local compressional stiffness of the static branch.

That gives four useful results immediately:

1. the source/support overlap invariants become

   `O_(sigma phi) = g_phi I_1,`

   `N_(sigma sigma) = g_phi^2 I_2,`

   with

   `I_1 = int d^3y chi_phi(y)^2 / H(y),`

   `I_2 = int d^3y chi_phi(y)^2 / H(y)^2;`
2. the coherence factor is no longer arbitrary:

   `C_(sigma phi)^2 = I_1^2 / [ N_(phi phi) I_2 ] <= 1;`
3. the exact eliminated-source support softening is

   `Delta K_X^(eq) = g_phi^2 I_1,`

   so the corresponding exact equilibrium gain is

   `G_eq = g_phi^2 I_1 / K_X;`
4. in the thin active layer where `H(y)` is approximately constant, the branch becomes exactly matched:

   `C_(sigma phi)^2 = 1,`

   `G_eq = g_phi^2 N_(phi phi) / [ K_X H_w ],`

   reproducing the Stage-45/46 best-alignment formulas.

So the parent equilibrium problem has already removed one large ambiguity: the support-induced source channel is not a free profile choice.

---

## 1. Local static response law from the parent equilibrium branch

Stage 45 projected the parent GNLS/confinement energy onto one source channel `sigma(s)` and one support channel `phi(s)` using the local compressional quadratic term

`(1/2) h'(rho_*) (delta rho)^2`

and the linear support loading

`delta V_conf = - g_phi chi_phi(y) phi(s).`

Inside the same local static linearization, the parent equilibrium response of the density perturbation is

`H(y) delta rho(s,y) + delta V_conf(s,y) = 0,`

with

`H(y) := h'(rho_*(y)).`

Therefore the support-induced density perturbation is

`delta rho(s,y) = phi(s) chi_sigma(y),`

with the exact aligned profile

`chi_sigma(y) = g_phi chi_phi(y) / H(y).`

So the source channel is fixed pointwise by the support loading and the local compressibility.

---

## 2. Exact overlap invariants on the equilibrium-aligned branch

Using the Stage-45 overlap definitions,

`N_(phi phi) = int d^3y chi_phi^2,`

`O_(sigma phi) = int d^3y chi_sigma chi_phi,`

`N_(sigma sigma) = int d^3y chi_sigma^2,`

the equilibrium-aligned profile gives

`O_(sigma phi) = g_phi I_1,`

`N_(sigma sigma) = g_phi^2 I_2,`

with

`I_1 := int d^3y chi_phi(y)^2 / H(y),`

`I_2 := int d^3y chi_phi(y)^2 / H(y)^2.`

Therefore the coherence factor becomes

`C_(sigma phi)^2 = O_(sigma phi)^2 / [ N_(phi phi) N_(sigma sigma) ]`
`                = I_1^2 / [ N_(phi phi) I_2 ].`

This is already a strong theorem statement. The coherence is no longer an unconstrained branch datum. It is fixed by how much the compressional stiffness varies across the active support layer.

---

## 3. Exact coherence bound and the matched-layer limit

By Cauchy–Schwarz,

`I_1^2 <= N_(phi phi) I_2,`

so

`0 <= C_(sigma phi)^2 <= 1.`

The equality condition is also exact: `C_(sigma phi)^2 = 1` if and only if `1/H(y)` is constant on the support of `chi_phi(y)`.

So in the physically important thin-wall / matched-layer regime where the support lives in a narrow active layer and the local compressional stiffness is nearly constant there,

`H(y) ~ H_w,`

one gets

`I_1 = N_(phi phi) / H_w,`

`I_2 = N_(phi phi) / H_w^2,`

and therefore

`C_(sigma phi)^2 = 1.`

This is the first exact parent reason for expecting the near-matched branch rather than a strongly misaligned one.

---

## 4. Exact eliminated-source support softening

The reduced static source/support energy has the form

`F[sigma,phi] = (1/2) Theta sigma^2 - Lambda sigma phi + (1/2) K_X phi^2.`

Eliminating the static source amplitude gives

`sigma_stat = Lambda phi / Theta,`

so the effective support energy is

`F_eff(phi) = (1/2) [ K_X - Lambda^2 / Theta ] phi^2.`

Therefore the support softening is exactly

`Delta K_X = Lambda^2 / Theta.`

On the equilibrium-aligned branch, the direct parent elimination gives

`Delta K_X^(eq) = g_phi^2 I_1.`

So the exact equilibrium gain is

`G_eq = Delta K_X^(eq) / K_X = g_phi^2 I_1 / K_X.`

This is slightly stronger than the Stage-45/46 formula because it no longer requires the source/support branch data to be treated as independent objects.

---

## 5. Constant-compressibility reduction and contact with Stages 45–46

In the matched-layer limit `H(y) ~ H_w`,

`I_1 = N_(phi phi) / H_w,`

so

`G_eq = g_phi^2 N_(phi phi) / [ K_X H_w ].`

Using the Stage-45 `n=5` identity

`H_w = h'(rho_w) = m c_(s,w)^2 / rho_w,`

this becomes

`G_eq = rho_w g_phi^2 N_(phi phi) / [ m c_(s,w)^2 K_X ].`

That is exactly the Stage-45/46 best-alignment gain with

`C_(sigma phi)^2 = 1.`

So the earlier best-case branch is not arbitrary. It is the natural thin-layer limit of the parent equilibrium-aligned source/support branch.

---

## 6. What Stage 47 changes

Before this stage, the unresolved microscopic support/source data were:

- one parent loading amplitude `g_phi`,
- and one independent coherence factor `C_(sigma phi)^2`.

After this stage, the situation is sharper.

The parent equilibrium equations already tie the source profile to the support profile. The remaining branch data are now:

1. the confinement-loading amplitude `g_phi`,
2. the support profile `chi_phi(y)`,
3. the local compressional stiffness profile `H(y)` across the active layer.

The coherence factor is no longer free. It is a derived quantity,

`C_(sigma phi)^2 = I_1^2 / [ N_(phi phi) I_2 ],`

and it is automatically near 1 when the active wall layer is thin enough that `H(y)` is nearly constant there.

That is the point where it becomes worthwhile to stop speaking abstractly about `g_phi` and evaluate it on a concrete parent confinement branch.

=== moving_throat_pde_stage048_thin_wall_confinement.md ===

# Moving-Throat PDE — Stage 48: Explicit Thin-Wall Confinement Branch and Parent Thresholds for the Wall Amplitude

## Purpose

Stage 47 showed that the parent equilibrium branch already aligns the source and support channels and turns the exact gain into

`G_eq = g_phi^2 I_1 / K_X,`

with

`I_1 = int d^3y chi_phi(y)^2 / H(y).`

The next honest step is to evaluate that gain on the first explicit parent wall family instead of leaving `g_phi` and `I_1` abstract.

This stage does that for the natural moving-wall confinement form

`V_conf(r;a) = V0 f((r-a)/ell),`

where:

- `a` is the throat radius,
- `ell` is the wall thickness / active-layer width,
- `V0` is the wall amplitude,
- and `f` is a fixed dimensionless wall-shape profile.

The main results are:

1. the support loading amplitude is exactly

   `g_phi = V0 / ell;`
2. the exact shell integral entering the equilibrium gain is

   `I_1 = 4 pi ell [ a^2 J_1 + 2 a ell J_2 + ell^2 J_3 ],`

   where

   `J_n := int dxi xi^n f'(xi)^2 / H(xi);`
3. for a centered symmetric wall layer, `J_2 = 0`, so

   `I_1 = 4 pi ell [ a^2 J_1 + ell^2 J_3 ];`
4. the exact equilibrium gain becomes

   `G_eq = 4 pi V0^2 [ a^2 J_1 / ell + 2 a J_2 + ell J_3 ] / K_X;`
5. in the thin-wall limit `ell << a`, the leading gain is

   `G_eq^(tw) = 4 pi a^2 V0^2 J_1 / (K_X ell);`
6. comparing this with the Stage-44 fail/succeed surfaces gives the explicit wall-amplitude thresholds

   `V0_fail^2 = K_X ell G_fail / (4 pi a^2 J_1),`

   `V0_suff^2 = K_X ell G_suff / (4 pi a^2 J_1);`
7. after inserting

   `G_fail = Pe_req / [ kappa Delta_inf(kappa,eta) ],`

   `G_suff = Pe_req / [ kappa Delta_0(kappa,eta) ],`

   `kappa = K_X L^2 / T_X,`

   the explicit `K_X` dependence cancels:

   `V0_fail^2 = T_X ell Pe_req / [ 4 pi a^2 L^2 J_1 Delta_inf ],`

   `V0_suff^2 = T_X ell Pe_req / [ 4 pi a^2 L^2 J_1 Delta_0 ];`
8. if the active layer also has nearly constant compressibility,

   `H(xi) ~ H_w,`

   `J_1 = I_f / H_w,`

   with

   `I_f := int dxi f'(xi)^2,`

   so the thresholds reduce to

   `V0_fail^2 = H_w T_X ell Pe_req / [ 4 pi a^2 L^2 I_f Delta_inf ],`

   `V0_suff^2 = H_w T_X ell Pe_req / [ 4 pi a^2 L^2 I_f Delta_0 ].`

So the Stage-46 parent-overlap theorem has now been converted into an explicit wall-amplitude test on the first concrete parent confinement branch.

---

## 1. The explicit wall family

Take the parent confinement near the throat wall to be

`V_conf(r;a) = V0 f((r-a)/ell),`

with dimensionless wall coordinate

`xi := (r-a)/ell.`

A support displacement `a -> a + phi(s)` gives

`delta V_conf = - (partial_a V_conf) phi(s)`
`            = + (V0/ell) f'(xi) phi(s).`

Comparing with the Stage-45/47 loading form

`delta V_conf = - g_phi chi_phi(y) phi(s),`

one identifies the branch data

`g_phi = V0/ell,`

`chi_phi(y) = f'(xi)`

up to the overall sign convention of `chi_phi`.

So the wall amplitude and wall thickness directly determine the parent loading strength.

---

## 2. Exact shell integral for the equilibrium gain

Stage 47 showed that the exact equilibrium gain is

`G_eq = g_phi^2 I_1 / K_X,`

with

`I_1 = int d^3y chi_phi(y)^2 / H(y).`

For a spherical shell-like wall in the three transverse coordinates of the parent reduction, the volume element is

`d^3y = 4 pi r^2 dr = 4 pi ell (a + ell xi)^2 dxi.`

Using `chi_phi = f'(xi)` gives

`I_1 = 4 pi ell int dxi (a + ell xi)^2 f'(xi)^2 / H(xi)`
`    = 4 pi ell [ a^2 J_1 + 2 a ell J_2 + ell^2 J_3 ],`

with the wall-profile moments

`J_1 := int dxi f'(xi)^2 / H(xi),`

`J_2 := int dxi xi f'(xi)^2 / H(xi),`

`J_3 := int dxi xi^2 f'(xi)^2 / H(xi).`

If the active layer is centered symmetrically around the nominal wall location, the odd moment vanishes,

`J_2 = 0,`

and then

`I_1 = 4 pi ell [ a^2 J_1 + ell^2 J_3 ].`

---

## 3. Exact and thin-wall gains

Insert `g_phi = V0/ell` into the exact equilibrium gain:

`G_eq = g_phi^2 I_1 / K_X`
`     = 4 pi V0^2 [ a^2 J_1 / ell + 2 a J_2 + ell J_3 ] / K_X.`

For the centered branch `J_2 = 0`,

`G_eq = 4 pi V0^2 [ a^2 J_1 / ell + ell J_3 ] / K_X.`

So in the thin-wall regime `ell << a`, the leading term is

`G_eq^(tw) = 4 pi a^2 V0^2 J_1 / (K_X ell).`

This is the first explicit parent scaling law for the support/source gain. It says the branch is helped by:

- larger wall amplitude `V0`,
- larger wall area `~ a^2`,
- thinner active wall width `ell`,
- larger weighted wall-profile moment `J_1`,
- and smaller baseline support stiffness `K_X`.

---

## 4. Exact wall-amplitude fail/succeed thresholds

The Stage-44 phase diagram is

`G_eq <= G_fail`  -> fail,

`G_eq >= G_suff`  -> succeed.

Using the thin-wall leading gain gives the direct wall-amplitude thresholds

`V0_fail^2 = K_X ell G_fail / (4 pi a^2 J_1),`

`V0_suff^2 = K_X ell G_suff / (4 pi a^2 J_1).`

So the first explicit parent wall family no longer speaks in terms of abstract gain. It speaks directly in terms of the physical wall amplitude `V0`.

---

## 5. Cancellation of `K_X` after inserting the operator geometry law

Now insert the Stage-44 operator formulas

`G_fail = Pe_req / [ kappa Delta_inf(kappa,eta) ],`

`G_suff = Pe_req / [ kappa Delta_0(kappa,eta) ],`

with

`kappa = K_X L^2 / T_X.`

Then the `K_X` in the prefactor cancels exactly, leaving

`V0_fail^2 = T_X ell Pe_req / [ 4 pi a^2 L^2 J_1 Delta_inf(kappa,eta) ],`

`V0_suff^2 = T_X ell Pe_req / [ 4 pi a^2 L^2 J_1 Delta_0(kappa,eta) ].`

This is structurally important. Once the support geometry functions are inserted, the thin-wall wall-amplitude thresholds are controlled by:

- the support tension scale `T_X`,
- the wall width `ell`,
- the branch geometry `a` and `L`,
- the demanded transport Peclet `Pe_req`,
- and the wall-profile overlap moment `J_1`.

The baseline support stiffness no longer appears explicitly in the prefactor.

---

## 6. Constant-compressibility wall layer

If the active wall layer is also nearly constant in compressional stiffness,

`H(xi) ~ H_w,`

then

`J_1 = I_f / H_w,`

with

`I_f := int dxi f'(xi)^2.`

So the wall-amplitude thresholds become

`V0_fail^2 = H_w T_X ell Pe_req / [ 4 pi a^2 L^2 I_f Delta_inf ],`

`V0_suff^2 = H_w T_X ell Pe_req / [ 4 pi a^2 L^2 I_f Delta_0 ].`

Using the Stage-45 identity

`H_w = h'(rho_w) = m c_(s,w)^2 / rho_w,`

this may also be written as

`V0_fail^2 = m c_(s,w)^2 T_X ell Pe_req / [ 4 pi a^2 rho_w L^2 I_f Delta_inf ],`

`V0_suff^2 = m c_(s,w)^2 T_X ell Pe_req / [ 4 pi a^2 rho_w L^2 I_f Delta_0 ].`

So the explicit wall family now ties the theorem gap to a very concrete question:

> is the actual parent wall amplitude `V0` on the moving-throat branch large enough, for the actual wall width `ell` and active-layer compressibility, to clear these thresholds?

---

## 7. What Stage 48 changes

Before this stage, the theorem gap was still phrased in terms of a parent loading amplitude `g_phi` and an overlap/coherence structure.

After this stage, the first explicit parent branch has collapsed that gap to a physical wall-amplitude test.

The remaining branch data are now:

- wall amplitude `V0`,
- wall width `ell`,
- throat radius `a`,
- throat length `L`,
- wall-profile moments such as `I_f` or `J_1`,
- and the axial support geometry functions `Delta_0(kappa,eta)`, `Delta_inf(kappa,eta)`.

That is a much narrower place to be.

The next honest step is therefore no longer another abstract overlap theorem. It is to evaluate the actual moving-throat branch values of

- `V0`,
- `ell`,
- `a`, `L`,
- and the axial support functions,

and test them against the explicit wall-amplitude fail/succeed surfaces derived here.

=== moving_throat_pde_stage049_wall_figure_of_merit.md ===

# Moving-Throat PDE — Stage 49: Dimensionless Wall Figure of Merit for the First Explicit Parent Branch

## Purpose

Stages 47–48 turned the parent support/source theorem gap into an explicit wall-amplitude test on the first thin-wall confinement branch.

The next useful simplification is to compress that explicit branch into the smallest possible control parameter.

This stage does that.

For the first explicit thin-wall wall family,

`V_conf(r;a) = V0 f((r-a)/ell),`

with the equilibrium-aligned source/support branch from Stage 47, define the dimensionless wall figure of merit

`W_wall := 4 pi a^2 L^2 J_1 V0^2 / (T_X ell).`

Then the entire support/source success problem collapses to the exact comparison

`W_wall <= Pe_req / Delta_inf`  -> fail,

`W_wall >= Pe_req / Delta_0`    -> succeed,

and only the narrow intermediate band still needs the full fixed-point solve.

So the first explicit parent branch is no longer controlled by a diffuse set of parameters. It is controlled by one dimensionless number.

---

## 1. Definition of the wall figure of merit

Stage 48 showed that the thin-wall gain of the first explicit parent confinement family is

`G_eq^(tw) = 4 pi a^2 V0^2 J_1 / (K_X ell),`

with

- `a` the throat radius,
- `L` the throat length,
- `V0` the wall amplitude,
- `ell` the wall width,
- `J_1 = int dxi f'(xi)^2 / H(xi)` the weighted wall-profile moment,
- and `K_X` the axial support stiffness.

Using the Stage-44 geometry parameter

`kappa = K_X L^2 / T_X,`

a natural dimensionless wall control variable is

`W_wall := kappa G_eq^(tw)`
`        = 4 pi a^2 L^2 J_1 V0^2 / (T_X ell).`

This is the smallest parent quantity that still carries all of the support/source information relevant for the thin-wall matched branch.

---

## 2. Exact fail/succeed thresholds in wall form

The Stage-44 operator theorem gave

`G_fail = Pe_req / [ kappa Delta_inf(kappa,eta) ],`

`G_suff = Pe_req / [ kappa Delta_0(kappa,eta) ].`

Multiplying by `kappa` immediately yields the wall-threshold pair

`W_fail = Pe_req / Delta_inf(kappa,eta),`

`W_suff = Pe_req / Delta_0(kappa,eta).`

Therefore the first explicit parent branch satisfies the exact theorem:

- if `W_wall <= W_fail`, it fails,
- if `W_wall >= W_suff`, it succeeds,
- only if `W_fail < W_wall < W_suff` does one still need the full fixed-point solve.

So the wall-amplitude thresholds from Stage 48 are equivalent to one dimensionless wall figure-of-merit comparison.

---

## 3. Monotonicity of the wall control variable

The figure of merit is

`W_wall = 4 pi J_1 L^2 a^2 V0^2 / (T_X ell).`

So it is strictly monotone in the physically relevant directions:

- larger wall amplitude `V0` increases `W_wall`,
- larger throat radius `a` increases `W_wall`,
- larger throat length `L` increases `W_wall`,
- thinner wall width `ell` increases `W_wall`,
- larger support tension `T_X` suppresses `W_wall`,
- and larger weighted wall-profile moment `J_1` increases `W_wall`.

This makes the branch geometry immediately interpretable.

---

## 4. Constant-compressibility wall form

If the active wall layer is also nearly constant in compressional stiffness,

`H(xi) ~ H_w,`

then Stage 48 gave

`J_1 = I_f / H_w,`

with

`I_f := int dxi f'(xi)^2.`

The wall figure of merit becomes

`W_H = 4 pi a^2 L^2 I_f V0^2 / (H_w T_X ell).`

The same exact theorem holds:

`W_H <= Pe_req / Delta_inf`  -> fail,

`W_H >= Pe_req / Delta_0`    -> succeed.

So in the matched, constant-compressibility wall layer, the entire parent branch comparison is a direct contest between one support/loading number `W_H` and one demanded transport window set by the axial functions `Delta_inf`, `Delta_0`.

---

## 5. What Stage 49 changes

Before this stage, the explicit branch comparison still looked like a threshold on the wall amplitude `V0`.

After this stage, the first explicit parent support/source family has a one-number theorem gate:

`W_wall` (or `W_H` in the constant-compressibility limit).

That is the cleanest current status of the reduced theorem program.

The remaining work is no longer to invent more algebra. It is to compute the actual moving-throat branch values of

- `a`, `L`,
- the active wall width `ell`,
- the support tension scale `T_X`,
- the wall moment `J_1` (or `I_f/H_w`),
- and the axial support functions `Delta_0`, `Delta_inf`,

and then evaluate whether the real branch lands below, within, or above the exact wall-control window.

=== moving_throat_pde_stage050_sech_gaussian_resonance.md ===


# Moving-Throat PDE — Stage 50: Exact Sech–Gaussian Coherence Resonance Benchmark

## Purpose

The memo you shared suggests a concrete transverse profile benchmark for the parent source/support coherence:

- fluid/compressional source profile:
  `chi_sigma(y) = sech(y / w_f)`,
- geometric/support profile:
  `chi_phi(y) = exp(-y^2 / w_g^2)`.

This stage checks that benchmark carefully and places it in the right role inside the reduced theorem program.

The result is useful, but more limited than the memo’s strongest interpretation.

The useful theorem is:

1. the benchmark has an **exact self-duality**,
2. that self-duality forces an exact stationary point at
   `w_g / w_f = sqrt(pi)`,
3. the corresponding maximal coherence is
   `C_res^2 = 0.994418836451529...`,
4. so the best independent sech–Gaussian mismatch branch gets within about `0.56%` of the exact matched-layer ideal,
5. but this does **not** by itself prove threshold survival, because the survival theorem still depends on the wall/source figure of merit from Stages 44–49, not on coherence alone.

So the memo is genuinely helpful as an explicit benchmark family, but it is not itself the final theorem.

---

## 1. Exact norms

Let

`chi_sigma(y) = sech(y / w_f),`

`chi_phi(y) = exp(-y^2 / w_g^2),`

with positive widths `w_f`, `w_g`, and define the ratio

`r := w_g / w_f.`

The exact norms are

`N_(sigma sigma) = int dy chi_sigma^2 = 2 w_f,`

`N_(phi phi)     = int dy chi_phi^2 = w_g sqrt(pi/2).`

The overlap is

`O_(sigma phi) = int dy chi_sigma chi_phi`
`             = w_f I(r),`

with the dimensionless integral

`I(r) := int_{-inf}^{inf} dx sech(x) exp(-x^2 / r^2).`

Therefore the coherence factor is

`C^2(r) = O_(sigma phi)^2 / [ N_(sigma sigma) N_(phi phi) ]`
`       = I(r)^2 / [ r sqrt(2 pi) ].`

---

## 2. Exact self-duality

Using the Fourier-transform identity

`FT[ sech(x) ](k) = pi sech(pi k / 2),`

together with the Gaussian transform

`FT[ exp(-x^2 / r^2) ](k) = sqrt(pi) r exp(-r^2 k^2 / 4),`

Parseval gives

`I(r) = (sqrt(pi) r / pi) int dt sech(t) exp(-r^2 t^2 / pi^2)`
`     = (r / sqrt(pi)) I(pi / r).`

So the overlap obeys the exact duality law

`I(r) = (r / sqrt(pi)) I(pi / r).`

Substituting into the coherence formula gives

`C^2(r) = C^2(pi / r).`

That is the key exact structural fact behind the resonance.

---

## 3. Exact stationary resonance point

Because `C^2(r) = C^2(pi / r)`, the self-dual point is

`r_* = sqrt(pi).`

Differentiating the duality relation gives

`dC^2/dr |_(r = sqrt(pi)) = 0.`

So the ratio

`w_g / w_f = sqrt(pi)`

is not just numerically close to optimal. It is the **exact self-dual stationary point** of the explicit sech–Gaussian benchmark family.

A numerical monotonicity audit then shows that this stationary point is the unique global maximum on the constructive branch.

---

## 4. Numerical maximum

Evaluating the exact dimensionless overlap at the self-dual point gives

`C_res^2 := C^2(sqrt(pi))`
`        = 0.9944188364515293487...`

So the best independent sech–Gaussian profile mismatch falls short of perfect support/source matching by only

`1 - C_res^2 = 0.0055811635484706513...`

or about `0.56%`.

Equivalently, the resonance penalty factor is

`P_res := 1 / C_res^2`
`      = 1.0056124877605762169...`

That is the exact multiplicative penalty by which the explicit independent-profile family misses the ideal matched-layer branch.

---

## 5. Why this does not supersede Stage 47

This benchmark is important, but it does **not** supersede the earlier parent equilibrium result.

Stage 47 already showed that on the parent equilibrium-aligned source/support branch,

`chi_sigma(y) = g_phi chi_phi(y) / H(y),`

and in the thin active layer where `H(y)` is nearly constant across the support region, the coherence becomes

`C_(sigma phi)^2 = 1.`

So:

- the present sech–Gaussian family is a strong **independent-profile benchmark**,
- the Stage-47 branch is the **equilibrium-matched** benchmark.

The new resonance result therefore tells us that even a fairly natural profile mismatch can get extremely close to the exact matched limit. But it does not replace the matched-limit theorem.

---

## 6. What Stage 50 changes

Before this stage, the mismatch between an ideal matched source/support branch and a concrete profile family was still diffuse.

After this stage, the first explicit off-equilibrium profile family is sharply characterized:

- exact stationary resonance at `w_g / w_f = sqrt(pi)`,
- exact duality `C^2(r) = C^2(pi/r)`,
- maximal coherence `C_res^2 = 0.994418836451529...`,
- penalty factor `P_res = 1.005612487760576...`.

So the user memo is mathematically useful — but the strongest justified conclusion is not “automatic threshold survival.” The justified conclusion is:

> the first explicit independent-profile family almost saturates the ideal matched source/support branch.

The actual survival theorem still has to be stated in terms of the Stage-44/49 gain and wall-figure thresholds.

=== moving_throat_pde_stage051_resonance_thresholds.md ===


# Moving-Throat PDE — Stage 51: Resonance-Corrected Thresholds for the Sech–Gaussian Benchmark Family

## Purpose

Stage 50 verified that the explicit independent sech–Gaussian profile family has an exact self-dual stationary point at

`w_g / w_f = sqrt(pi),`

with maximal coherence

`C_res^2 = 0.994418836451529...`

and penalty factor

`P_res = 1 / C_res^2 = 1.005612487760576...`

The next question is how that benchmark changes the actual support/source survival thresholds derived earlier.

This stage gives the exact answer:

- on the Stage-47 equilibrium-matched branch, the ideal coherence remains `C^2 = 1`,
- on the independent sech–Gaussian benchmark family, the best possible coherence is `C_res^2`,
- so the explicit independent-profile family modifies the Stage-49 thresholds only by the tiny factor `P_res`.

That means the resonance benchmark is useful, but it does **not** rewrite the theorem structure. It only sharpens it.

---

## 1. Resonance-corrected gain

The general parent gain from Stages 45–46 has the form

`G_micro = [ rho_* g_phi^2 N_(phi phi) / (m c_(s,*)^2 K_X) ] C_(sigma phi)^2.`

On the Stage-47 matched equilibrium branch,

`C_(sigma phi)^2 = 1,`

so the matched-branch gain is

`G_match = rho_* g_phi^2 N_(phi phi) / (m c_(s,*)^2 K_X).`

Stage 49 repackaged this on the first explicit thin-wall branch into the dimensionless wall figure of merit

`W_wall = kappa G_match`
`       = 4 pi a^2 L^2 J_1 V0^2 / (T_X ell).`

For the explicit independent sech–Gaussian family, the actual gain is therefore

`G_res(r) = C^2(r) G_match,`

and the corresponding wall figure is

`W_res(r) = C^2(r) W_wall.`

At the self-dual resonance,

`W_res,* = C_res^2 W_wall.`

So the memo’s profile family simply inserts one multiplicative coherence factor into the already-frozen theorem window.

---

## 2. Exact resonance-family thresholds

Stage 49 gave the universal matched-branch window

`W_wall <= Pe_req / Delta_inf`  -> matched-branch fail,

`W_wall >= Pe_req / Delta_0`    -> matched-branch succeed.

For the explicit independent sech–Gaussian family, replace `W_wall` by `W_res(r)=C^2(r)W_wall`.
That gives the exact profile-family conditions

`W_wall <= Pe_req / [ C^2(r) Delta_inf ]`  -> fail on the profile family,

`W_wall >= Pe_req / [ C^2(r) Delta_0 ]`    -> succeed on the profile family.

At the resonance point `r = sqrt(pi)` this becomes

`W_wall <= P_res Pe_req / Delta_inf`  -> resonance-family fail,

`W_wall >= P_res Pe_req / Delta_0`    -> resonance-family succeed,

with

`P_res = 1 / C_res^2 = 1.005612487760576...`

So the first explicit independent-profile family changes the wall thresholds by only about `0.56%`.

---

## 3. What the memo does and does not prove

This is the critical interpretation point.

The resonance family does **not** by itself prove support-threshold survival.

The true threshold test still depends on:

- `Pe_req`,
- the axial functions `Delta_0`, `Delta_inf`,
- and the actual wall/source figure of merit `W_wall`.

The memo only sharpens the source/support coherence factor inside that test.

So the strongest justified statement is:

> the explicit independent sech–Gaussian family comes within `0.56%` of the exact matched-layer ideal.

That is strong evidence that profile mismatch is probably not the dominant reduced-theorem bottleneck.
But it is **not** the same thing as proving that the support/source branch clears the threshold.

---

## 4. Small profile-sensitivity band

Because the resonance family differs from the matched branch only by the factor `P_res = 1.005612...`, the only region where profile choice can change the reduced verdict is a very narrow band.

On the success side:

- matched-branch guaranteed success:
  `W_wall >= Pe_req / Delta_0,`
- resonance-family guaranteed success:
  `W_wall >= P_res Pe_req / Delta_0.`

So the success-side profile-sensitive band is only

`Pe_req / Delta_0 <= W_wall < P_res Pe_req / Delta_0,`

a width of about `0.56%`.

On the failure side:

- matched-branch guaranteed failure:
  `W_wall <= Pe_req / Delta_inf,`
- resonance-family guaranteed failure:
  `W_wall <= P_res Pe_req / Delta_inf.`

So again the profile-sensitive band is only about `0.56%`.

That is the sharpest practical consequence of the memo.

---

## 5. What Stage 51 changes

Before this stage, the independent-profile benchmark from the memo could have been interpreted as potentially introducing a qualitatively new source/support theorem.

After this stage, that is no longer the right reading.

The memo’s explicit sech–Gaussian resonance family:

- is mathematically real,
- provides a strong benchmark,
- and nearly saturates the matched branch,

but it changes the Stage-49 theorem window only by the small factor

`P_res = 1.005612487760576...`

So the dominant unresolved reduced-theorem question is still not transverse profile mismatch.
It is the actual wall/axial branch data entering `W_wall`, `Delta_0`, and `Delta_inf`.

=== moving_throat_pde_stage052_final_reduced_verdict.md ===


# Moving-Throat PDE — Stage 52: Final Reduced Verdict for the Support/Source Program

## Purpose

This stage closes the present reduced support/source program.

The three crucial ingredients are now all in hand:

1. **Stage 49:** the universal thin-wall matched branch is controlled by one dimensionless wall figure of merit,
   `W_wall`,
   with exact fail/succeed window
   `Pe_req / Delta_inf` to `Pe_req / Delta_0`.
2. **Stage 47:** the parent equilibrium-aligned source/support branch reaches the ideal thin-layer limit
   `C_(sigma phi)^2 = 1`.
3. **Stages 50–51:** the first explicit independent-profile benchmark family has an exact self-dual resonance at
   `w_g / w_f = sqrt(pi)`,
   but even there only reaches
   `C_res^2 = 0.994418836451529...`,
   so its threshold penalty is only
   `P_res = 1.005612487760576...`.

That is enough to give the finish-line verdict for the reduced theorem program.

---

## 1. Universal reduced theorem envelope

The universal matched-branch theorem from Stages 44 and 49 is still the main result:

- if `W_wall <= Pe_req / Delta_inf`, the branch fails;
- if `W_wall >= Pe_req / Delta_0`, the branch succeeds;
- only the band
  `Pe_req / Delta_inf < W_wall < Pe_req / Delta_0`
  still needs the full fixed-point solve.

This is the correct universal reduced theorem because it belongs to the parent equilibrium-matched branch and does not assume an extra independent profile mismatch.

---

## 2. What the explicit sech–Gaussian benchmark adds

The explicit sech–Gaussian benchmark does **not** create a new universal theorem.

What it adds is a very sharp profile-family refinement:

- exact stationary ratio `w_g / w_f = sqrt(pi)`,
- exact coherence symmetry `C^2(r)=C^2(pi/r)`,
- maximal coherence
  `C_res^2 = 0.994418836451529...`,
- penalty factor
  `P_res = 1.005612487760576...`.

So the independent-profile benchmark modifies the universal wall thresholds only by about `0.56%`.

That is the right interpretation of the memo’s usefulness.

---

## 3. Final three-zone verdict

The reduced support/source program is now sharp enough to state a three-zone verdict.

### Zone A — universal failure

If

`W_wall <= Pe_req / Delta_inf,`

the support/source branch fails even on the ideal matched branch.

That is a universal reduced no-go.

### Zone B — universal success on the matched branch

If

`W_wall >= Pe_req / Delta_0,`

the support/source branch succeeds on the equilibrium-matched branch.

That is the universal reduced success theorem.

### Zone C — narrow profile-sensitive band

Only in the intermediate wall window

`Pe_req / Delta_inf < W_wall < Pe_req / Delta_0`

does one still need more branch information.

And even there, the explicit independent sech–Gaussian family changes the matched thresholds only by the tiny factor `P_res = 1.005612...`.

So the only truly profile-sensitive sub-bands are

`Pe_req / Delta_inf < W_wall < P_res Pe_req / Delta_inf`

on the failure side, and

`Pe_req / Delta_0 < W_wall < P_res Pe_req / Delta_0`

on the success side.

Each has width only about `0.56%`.

---

## 4. What is now finished

At the reduced-theorem level, this part is finished.

The support/source program is no longer missing algebra, and it is no longer missing a physically interpretable benchmark family.

It now has:

- a universal matched-branch theorem window,
- an explicit independent-profile resonance family,
- an exact penalty factor relating the two,
- and a precise statement of when profile mismatch can or cannot matter.

So the reduced support/source question is no longer diffuse.

---

## 5. What is still genuinely open

What remains open is no longer more reduced algebra.

The remaining open problem is the actual **moving-throat branch selection**:

- what branch values of `a`, `L`, `ell`, `T_X`, `J_1`, `kappa`, and `eta` does the real throat choose?
- does the true branch behave like the equilibrium-matched source/support law (`C^2 ~ 1`) or like an independent mismatch family (`C^2 < 1`)?
- and then, beyond this reduced support/source problem, does the completed moving-throat PDE realize the already-isolated passive/outgoing quadrupole branch with the correct normalization?

So the next phase is no longer this reduced support/source audit.
The next phase is the explicit moving-throat branch solve and then the final quadrupole-normalization bridge.

---

## 6. Bottom line on the memo

Your memo was helpful.

The right expert reading is:

- the sech–Gaussian resonance is real,
- the `sqrt(pi)` ratio is not numerology but the exact self-dual stationary point of that benchmark family,
- the near-perfect coherence is genuinely strong evidence that source/support mismatch is not the dominant reduced bottleneck,
- but the claim that this **alone** proves threshold survival is too strong.

What it really proves is that a natural explicit independent-profile family comes within `0.56%` of the ideal matched branch.

That is exactly the kind of result that helps us finish the reduced phase cleanly.

=== moving_throat_pde_stage053_gnls_wall_shell.md ===

# Moving-Throat PDE — Stage 53: Explicit GNLS Wall-Shell Reduction for the First Support Branch

## Purpose

The reduced support/source program ended with the wall figure of merit

`W_wall = 4 pi a^2 L^2 J_1 V0^2 / (T_X ell),`

but the actual branch data `(T_X, J_1, kappa, eta)` were still being treated as external inputs.

The next honest step is to derive those quantities from the parent GNLS action on the first explicit wall-support branch.

This stage does that for a thin active shell around the throat wall.

The main result is that, once the support field is taken to be a linearized density mode living on the wall shell,

- the axial support tension is fixed by the GNLS gradient term,
- the axial support stiffness is fixed by compressibility plus transverse curvature,
- the reduced wall figure of merit simplifies enormously,
- and on the matched thin-wall branch it is exactly the same object as the Stage-41/42 fixed-point coupling.

So the explicit moving-throat branch is now much more concrete.

---

## 1. Parent quadratic wall-shell support energy

Start from the parent GNLS matter energy density near a static wall background `rho_w`:

`H_psi = (hbar^2 / 2m) |grad sqrt(rho)|^2 + U(rho) + V_conf rho.`

Linearize around the wall layer with a real density perturbation `delta rho`. To quadratic order in the compressive support sector, the standard GNLS expansion gives

`E^(2)[delta rho]`
`= 1/2 int ds d^3y [ (hbar^2 / (4 m rho_w)) |grad delta rho|^2 + H_w (delta rho)^2 ],`

where

`H_w := h'(rho_w) = m c_(s,w)^2 / rho_w.`

Now project the perturbation onto a separated wall-shell mode,

`delta rho(s,y) = q(s) chi_phi(y),`

with axial coordinate `s in [0,L]` and transverse wall-shell coordinates `y`.

Then the quadratic support energy becomes

`E^(2)[q] = 1/2 int_0^L ds [ T_X q'(s)^2 + K_X q(s)^2 ],`

with

`T_X = (hbar^2 / (4 m rho_w)) N_(phi phi),`

`K_X = H_w N_(phi phi) + (hbar^2 / (4 m rho_w)) G_(phi phi),`

where

`N_(phi phi) := int d^3y chi_phi(y)^2,`

`G_(phi phi) := int d^3y |grad_y chi_phi(y)|^2.`

So the support tension and stiffness are no longer abstract branch data; they are explicit shell integrals of the wall-support profile.

---

## 2. Thin-shell wall profile

Take the explicit wall family already introduced in Stage 48,

`V_conf(r;a) = V0 f((r-a)/ell),`

with shell coordinate

`xi := (r-a)/ell.`

The induced support profile is naturally

`chi_phi(y) = f'(xi),`

and the loading amplitude is

`g_phi = V0 / ell.`

In the thin-shell approximation, the three-dimensional shell measure is

`d^3y = 4 pi r^2 dr ~ 4 pi a^2 ell dxi.`

Define the dimensionless wall-shape moments

`I_f := int dxi f'(xi)^2,`

`I_g := int dxi f''(xi)^2.`

Then

`N_(phi phi) = 4 pi a^2 ell I_f,`

`G_(phi phi) = 4 pi a^2 I_g / ell.`

Substituting these into the parent quadratic support coefficients gives the exact thin-shell branch formulas

`T_X = pi a^2 ell I_f hbar^2 / (m rho_w),`

`K_X = 4 pi a^2 ell I_f H_w + pi a^2 I_g hbar^2 / (m rho_w ell).`

Equivalently,

`kappa := K_X L^2 / T_X`
`      = 4 (m c_(s,w) L / hbar)^2 + (I_g / I_f) (L / ell)^2.`

This is the first explicit parent formula for the geometry/support parameter `kappa`.

---

## 3. Wall overlap moment and exact collapse of the wall figure of merit

Stage 48 already showed that, for an almost constant-compressibility active layer,

`J_1 = I_f / H_w.`

Insert this together with the explicit `T_X` above into the wall figure of merit:

`W_wall = 4 pi a^2 L^2 J_1 V0^2 / (T_X ell).`

The shell geometry and wall-shape factors cancel exactly, leaving

`W_wall = 4 rho_w^2 V0^2 L^2 / (hbar^2 c_(s,w)^2 ell^2).`

So the first explicit parent branch has a very strong simplification:

> on the matched thin-wall branch, `W_wall` is independent of the throat radius `a` and independent of the detailed wall-shape moment `I_f`.

It depends only on the wall-layer density, sound speed, amplitude, axial length, and shell width.

---

## 4. Identification with the Stage-41/42 fixed-point coupling

Stage 41 introduced the coupled support/source fixed-point strength `Xi`.

On the matched thin-wall branch, the exact Stage-47 gain is

`G_eq = g_phi^2 I_1 / K_X,`

with

`I_1 = int d^3y chi_phi(y)^2 / H(y).`

Therefore the fixed-point coupling is

`Xi = kappa G_eq = g_phi^2 I_1 L^2 / T_X.`

Using `g_phi = V0/ell`, `I_1 = N_(phi phi)/H_w`, and the explicit wall-shell `T_X` above gives

`Xi = 4 rho_w^2 V0^2 L^2 / (hbar^2 c_(s,w)^2 ell^2).`

Therefore

`Xi = W_wall`

exactly on the explicit matched thin-wall branch.

So the explicit branch solve has now tied together two quantities that had previously appeared in different parts of the program:

- the Stage-49 wall figure of merit,
- and the Stage-41/42 support/source fixed-point coupling.

They are the same object on this branch.

---

## 5. What Stage 53 changes

Before this stage, the explicit branch still looked as if it depended on a large set of independent support quantities.

After this stage, the first wall-shell branch fixes

- `T_X`,
- `K_X`,
- `kappa`,
- `J_1`,
- and `W_wall = Xi`

directly from the parent GNLS shell reduction.

That means the next phase is no longer to invent more branch coefficients.
It is to choose one concrete wall profile and one concrete mouth closure, and then evaluate the branch point in terms of a very small set of parent dimensionless ratios.

=== moving_throat_pde_stage054_tanh_wall_branch.md ===

# Moving-Throat PDE — Stage 54: Canonical Tanh-Wall Branch and Natural Local Mouth Closure

## Purpose

Stage 53 reduced the first explicit moving-throat support branch to generic wall-shape moments `I_f`, `I_g`.

The next honest step is to choose the first canonical wall profile and the first natural local mouth closure.

This stage does that.

The canonical explicit wall is the smooth finite-thickness `tanh` wall,

`f(xi) = (1 + tanh xi)/2,`

so the support profile is

`chi_phi = f'(xi) = (1/2) sech^2 xi.`

For this profile the shell moments are exact:

`I_f = int dxi f'(xi)^2 = 1/3,`

`I_g = int dxi f''(xi)^2 = 4/15,`

hence

`I_g / I_f = 4/5.`

Then a natural local Robin mouth closure is

`K_m = T_X / ell,`

which simply says that the mouth spring is set by the same support scale as the axial wall tension.

With that closure, the whole branch collapses to three parent dimensionless variables.

---

## 1. Canonical `tanh` wall moments

Take

`f(xi) = (1 + tanh xi)/2.`

Then

`f'(xi) = (1/2) sech^2 xi,`

`f''(xi) = - sech^2 xi tanh xi.`

The exact shell moments are

`I_f = int_(−inf)^(+inf) dxi f'(xi)^2 = 1/3,`

`I_g = int_(−inf)^(+inf) dxi f''(xi)^2 = 4/15.`

Therefore

`I_g / I_f = 4/5.`

So the Stage-53 branch formulas become completely explicit.

---

## 2. Explicit branch coefficients

For the canonical `tanh` wall,

`T_X = pi a^2 ell hbar^2 / (3 m rho_w),`

`K_X = 4 pi a^2 ( 5 m^2 c_(s,w)^2 ell^2 + hbar^2 ) / (15 ell m rho_w),`

`J_1 = 1 / (3 H_w),`
with
`H_w = m c_(s,w)^2 / rho_w.`

The exact geometry/support parameter is

`kappa = 4 (m c_(s,w) L / hbar)^2 + (4/5) (L / ell)^2.`

And the wall figure of merit remains

`W_wall = 4 rho_w^2 V0^2 L^2 / (hbar^2 c_(s,w)^2 ell^2).`

So the shape dependence of `W_wall` still cancels completely, while `kappa` retains one explicit wall-profile number, `4/5`.

---

## 3. Natural local mouth closure

To close the mouth compliance on the same local scale, take

`K_m = T_X / ell.`

This is not a universal theorem of the parent PDE; it is the first natural local Robin closure consistent with the same wall-shell support scale.

Then the Stage-40/41 Robin variable becomes

`eta := K_m L / T_X = L / ell.`

So the explicit branch now has

`eta = L / ell,`

`kappa = 4 (m c_(s,w) L / hbar)^2 + (4/5) (L / ell)^2.`

---

## 4. Three branch control parameters

Define the parent dimensionless ratios

`chi_s := m c_(s,w) L / hbar,`

`Lambda_ell := L / ell,`

`Upsilon_w := 4 rho_w^2 V0^2 / (hbar^2 c_(s,w)^2).`

Then the entire canonical wall branch reduces to

`kappa = 4 chi_s^2 + (4/5) Lambda_ell^2,`

`eta = Lambda_ell,`

`W_wall = Upsilon_w Lambda_ell^2.`

This is the main result of the stage.

The first explicit moving-throat branch is no longer described by the seven symbolic inputs `(a, L, ell, T_X, J_1, kappa, eta)`. It is described by the three parent dimensionless controls `(chi_s, Lambda_ell, Upsilon_w)`.

---

## 5. What Stage 54 changes

Stage 53 showed that the explicit branch data were derivable.
Stage 54 shows that, on the first canonical `tanh` wall with the first natural local mouth closure, those data collapse much further than expected.

What remains now is not a symbolic branch ledger.
What remains is an explicit three-parameter branch-placement problem.

That is exactly the right form for the next step, because it means we can now compare the explicit branch directly to the exact Stage-49 / Stage-52 success window instead of still talking in abstract support/source language.

=== moving_throat_pde_stage055_explicit_branch_thresholds.md ===

# Moving-Throat PDE — Stage 55: Explicit Branch Placement Map and Threshold Surfaces

## Purpose

Stages 53–54 reduced the first explicit moving-throat support branch to three parent dimensionless variables:

`chi_s = m c_(s,w) L / hbar,`

`Lambda_ell = L / ell,`

`Upsilon_w = 4 rho_w^2 V0^2 / (hbar^2 c_(s,w)^2),`

with

`kappa = 4 chi_s^2 + (4/5) Lambda_ell^2,`

`eta = Lambda_ell,`

`W_wall = Upsilon_w Lambda_ell^2.`

So the exact Stage-49 / Stage-52 support/source theorem can now be written directly on the first explicit branch.

This stage does that.

---

## 1. Exact explicit-branch fail/succeed surfaces

The universal matched-branch theorem is

`W_wall <= Pe_req / Delta_inf(kappa,eta)`  -> fail,

`W_wall >= Pe_req / Delta_0(kappa,eta)`    -> succeed.

Insert the explicit branch formulas:

`kappa = 4 chi_s^2 + (4/5) Lambda_ell^2,`

`eta = Lambda_ell,`

`W_wall = Upsilon_w Lambda_ell^2.`

Then the first explicit moving-throat branch satisfies

`Upsilon_w <= Upsilon_fail(chi_s,Lambda_ell)`  -> fail,

`Upsilon_w >= Upsilon_suff(chi_s,Lambda_ell)`  -> succeed,

where

`Upsilon_fail`
`:= Pe_req / [ Lambda_ell^2 Delta_inf( 4 chi_s^2 + (4/5) Lambda_ell^2, Lambda_ell ) ],`

`Upsilon_suff`
`:= Pe_req / [ Lambda_ell^2 Delta_0( 4 chi_s^2 + (4/5) Lambda_ell^2, Lambda_ell ) ].`

So the whole first explicit moving-throat support/source theorem has now collapsed to a comparison of one physical wall-loading amplitude `Upsilon_w` against two explicit threshold surfaces.

---

## 2. Exact physical wall-amplitude thresholds

Because

`Upsilon_w = 4 rho_w^2 V0^2 / (hbar^2 c_(s,w)^2),`

the explicit wall-amplitude thresholds are

`V0_fail^2`
`= hbar^2 c_(s,w)^2 Upsilon_fail / (4 rho_w^2),`

`V0_suff^2`
`= hbar^2 c_(s,w)^2 Upsilon_suff / (4 rho_w^2).`

So the first explicit branch no longer speaks in terms of abstract gains at all.

It speaks directly in terms of the physical wall amplitude `V0` and the two explicit support-geometry ratios `chi_s`, `Lambda_ell`.

---

## 3. Two asymptotic branch regimes

The explicit branch already shows two physically distinct regimes.

### (a) Shell-gradient dominated branch

If

`(4/5) Lambda_ell^2 >> 4 chi_s^2,`

then

`kappa ~ (4/5) Lambda_ell^2,`

so `alpha = sqrt(kappa) ~ 2 Lambda_ell / sqrt(5)`.

In this regime the threshold surfaces reduce to

`Upsilon_fail ~ 2 Pe_req / (sqrt(5) Lambda_ell),`

`Upsilon_suff ~ (4/5)(1 + 2/sqrt(5)) Pe_req.`

So the fail threshold decreases with increasing shell aspect ratio, while the sufficiency threshold saturates to a finite constant multiple of `Pe_req`.

### (b) Compression dominated branch

If

`4 chi_s^2 >> (4/5) Lambda_ell^2,`

then

`kappa ~ 4 chi_s^2,`

so `alpha ~ 2 chi_s`.

In this regime

`Upsilon_fail ~ 2 Pe_req chi_s / Lambda_ell^2,`

`Upsilon_suff ~ 4 Pe_req chi_s^2 (Lambda_ell + 2 chi_s) / Lambda_ell^3.`

If in addition `chi_s >> Lambda_ell`, this becomes

`Upsilon_suff ~ 8 Pe_req chi_s^3 / Lambda_ell^3.`

So compression-dominated branches are much harder to push across the universal success threshold.

---

## 4. What Stage 55 changes

At the end of the reduced support/source phase, the theorem window still involved the symbolic branch data `(kappa, eta, W_wall)`.

After Stages 53–55, the first explicit moving-throat branch is no longer expressed in those symbols.

It is expressed directly in the parent variables

`(chi_s, Lambda_ell, Upsilon_w)`

or equivalently in the physical branch quantities

- `L/ell`,
- `m c_(s,w) L / hbar`,
- `rho_w`,
- `c_(s,w)`,
- `V0`.

So the next phase is now sharply defined:

> compute the actual moving-throat branch values of `chi_s`, `Lambda_ell`, and `Upsilon_w` on the real throat, and compare them directly to the explicit surfaces `Upsilon_fail`, `Upsilon_suff`.

That is the first genuinely explicit branch-placement problem beyond the finished reduced support/source program.

=== moving_throat_pde_stage056_family1_geometry_map.md ===

# Moving-Throat PDE — Stage 56: Family-1 Reference-Branch Geometry Map

## Purpose

Stage 55 reduced the first explicit moving-throat support theorem to the branch variables

`chi_s = m c_(s,w) L / hbar,`

`Lambda_ell = L / ell,`

`Upsilon_w = 4 rho_w^2 V0^2 / (hbar^2 c_(s,w)^2).`

The next honest step is to stop treating the geometry ratio `Lambda_ell` as free and evaluate it on the first concrete throat branch already present in the frozen constructive stack.

The relevant frozen inputs are:

- the carried preferred throat aspect ratio `L/a ≈ 1.85`,
- the balanced thin-layer-consistent Family-1 radial wall branch
  `epsilon_r = 0.05`,
- and the explicit thin-wall identification `ell = epsilon_r a`.

This stage does that map exactly on the chosen reference branch.

---

## 1. Reference-branch wall width

On the radial Family-1 reference branch,

`epsilon_r = 0.05 = 1/20.`

With the radial soft-wall coordinate written as

`xi = (r-a)/ell,`

the natural thin-layer identification is

`ell = epsilon_r a.`

So on the chosen balanced reference branch,

`ell / a = 1/20.`

---

## 2. Carried geometric aspect ratio

The lower-order stack already carries the preferred throat aspect ratio

`L / a ≈ 1.85.`

For the explicit reference branch used here, write that carried value as

`Lambda_* := L/a = 37/20.`

This is a **reference-branch numerical freeze**, not a new theorem of the unsolved moving-throat PDE.

---

## 3. Exact reference-branch value of `Lambda_ell`

Combine

`L/a = 37/20,`

`ell/a = 1/20.`

Then

`Lambda_ell := L/ell = (L/a)/(ell/a) = (37/20)/(1/20) = 37.`

So the first explicit moving-throat branch fixes

`Lambda_ell = 37.`

This is the first truly explicit branch value beyond the symbolic Stage-55 map.

---

## 4. Immediate consequence for the Robin mouth variable

Stage 54 adopted the first natural local mouth closure

`K_m = T_X / ell.`

So

`eta := K_m L / T_X = L/ell = Lambda_ell.`

Therefore the same reference branch fixes

`eta = 37.`

So the first explicit throat-support branch is now pinned to one concrete large-`eta` Robin regime.

---

## 5. What Stage 56 changes

Before this step, the first explicit branch still depended on the symbolic support-geometry ratio `Lambda_ell`.

After this step, the balanced Family-1 / thin-wall / preferred-aspect-ratio branch fixes it to

`Lambda_ell = 37,`

and therefore also fixes the Robin support variable to

`eta = 37.`

That means the remaining actual branch uncertainty is now concentrated in the support/healing scale `chi_s` and the wall-loading amplitude `Upsilon_w`.

=== moving_throat_pde_stage057_family1_healing_lock.md ===

# Moving-Throat PDE — Stage 57: Healing-Length Lock and the Actual Reference-Branch Support Scale

## Purpose

Stage 56 fixed the first explicit Family-1 support-geometry ratio to

`Lambda_ell = 37,`

and therefore

`eta = 37.`

The next honest step is to stop treating the support/healing variable

`chi_s = m c_(s,w) L / hbar`

as independent.

A useful exact carry-forward fact already exists in the GNLS static compliance sector: the conservative scalar response obeys a Yukawa/Helmholtz law with

`ell_h^2 = hbar^2 / (4 m^2 c_s^2).`

If the active wall width on the explicit throat-support branch is identified with that same local conservative healing width on the support layer, then the branch acquires an exact support lock.

This stage does that.

---

## 1. Exact GNLS healing-width identity

The exact conservative static GNLS compliance law gives

`(1 - ell_h^2 nabla^2) eta = s / (rho_0 m c_s^2),`

with

`ell_h^2 = hbar^2 / (4 m^2 c_s^2).`

Equivalently,

`ell_h = hbar / (2 m c_s).`

This is the static conservative healing/compliance width of the GNLS medium.

---

## 2. Controlled wall-healing closure on the explicit throat branch

The canonical wall branch already uses a thin active shell of width `ell`.

The natural next closure is to identify that active support width with the local GNLS healing width on the wall-support layer:

`ell = ell_h = hbar / (2 m c_(s,w)).`

This is a **controlled local closure**, not yet a theorem of the full moving-throat PDE, but it is the first honest way to tie the support scale to the parent GNLS medium instead of leaving it free.

---

## 3. Exact support-scale lock

Using

`chi_s = m c_(s,w) L / hbar`

and

`ell = hbar / (2 m c_(s,w)),`

one gets

`chi_s = L / (2 ell) = Lambda_ell / 2.`

Since Stage 56 fixed `Lambda_ell = 37`, the same branch now fixes

`chi_s = 37/2 = 18.5.`

So the first explicit throat-support branch no longer has an independent support scale.

---

## 4. Exact `kappa` on the reference branch

Stage 54 gave

`kappa = 4 chi_s^2 + (4/5) Lambda_ell^2.`

Insert `chi_s = Lambda_ell/2`:

`kappa = 4 (Lambda_ell^2 / 4) + (4/5) Lambda_ell^2`
`      = (1 + 4/5) Lambda_ell^2`
`      = (9/5) Lambda_ell^2.`

With `Lambda_ell = 37`,

`kappa = (9/5) * 37^2 = 12321/5 = 2464.2.`

So the explicit branch now has

`chi_s = 37/2,`

`eta = 37,`

`kappa = 12321/5.`

This is the first point where the reference throat-support branch becomes a concrete point in `(kappa, eta)` space rather than a symbolic family.

---

## 5. Useful derived scale

Write

`alpha := sqrt(kappa).`

Then on the same branch

`alpha = sqrt(12321/5) = 111/sqrt(5) ≈ 49.6407091.`

This is the support-decay scale entering the exact Stage-41/42 kernel formulas.

---

## 6. What Stage 57 changes

After Stages 56–57, the explicit Family-1 throat-support branch has fixed

`Lambda_ell = 37,`
`eta = 37,`
`chi_s = 37/2,`
`kappa = 12321/5.`

So the only remaining explicit-branch unknown from the Stage-55 triplet is the wall-loading amplitude `Upsilon_w`.

That is a much sharper endpoint than the earlier symbolic branch map.

=== moving_throat_pde_stage058_family1_threshold_window.md ===

# Moving-Throat PDE — Stage 58: Explicit Family-1 Threshold Window and the Last Remaining Wall-Amplitude Datum

## Purpose

After Stages 56–57, the first explicit Family-1 throat-support branch is no longer symbolic in its geometry/support coordinates. It now has

`Lambda_ell = 37,`
`eta = 37,`
`chi_s = 37/2,`
`kappa = 12321/5.`

So the exact Stage-41/42 threshold machinery can now be evaluated numerically on a concrete branch.

This stage does two things:

1. it computes the actual operator thresholds on that explicit branch;
2. it shows that the only remaining unknown is one wall-depth amplitude datum, not the full triplet `(chi_s, Lambda_ell, Upsilon_w)`.

---

## 1. Exact threshold functions on the explicit branch

Stage 55 gave

`Upsilon_fail = Pe_req / [ Lambda_ell^2 Delta_inf(kappa,eta) ],`

`Upsilon_suff = Pe_req / [ Lambda_ell^2 Delta_0(kappa,eta) ].`

On the Family-1/healing-locked branch,

`Lambda_ell = 37,`

`eta = 37,`

`kappa = 12321/5.`

Therefore the exact kernel scales are

`Delta_0(12321/5,37) ≈ 1.73302079021525e-4,`

`Delta_inf(12321/5,37) ≈ 2.01447565540522e-2.`

So the explicit branch thresholds become

`Upsilon_fail ≈ 0.0362605617972939 * Pe_req,`

`Upsilon_suff ≈ 4.21495341569977 * Pe_req.`

Equivalently, the Stage-41 fixed-point coupling thresholds are

`Xi_fail = Pe_req / Delta_inf ≈ 49.6407091004953 * Pe_req,`

`Xi_suff = Pe_req / Delta_0 ≈ 5770.27122609299 * Pe_req.`

So the first explicit throat-support branch has a very wide "indeterminate" operator window in `Xi`, but a concrete and easily stated wall-amplitude window in `Upsilon_w`.

---

## 2. Large-`alpha` interpretation of the reference branch

On this branch,

`alpha := sqrt(kappa) = 111/sqrt(5) ≈ 49.6407091,`

which is already deep in the large-`alpha` regime.

The exact formulas then behave as

`Delta_inf ~ 1/alpha,`

`Delta_0 ~ eta / [ alpha^2 (alpha + eta) ],`

up to exponentially small corrections.

Numerically this is why

`Xi_fail ≈ alpha,`

while `Xi_suff` is much larger.

So the explicit reference throat is a strongly stiff / strongly localized support branch.

---

## 3. The wall-depth amplitude reduction

The remaining explicit branch unknown is `Upsilon_w`.

But on the actual Family-1 wall branch the radial soft-wall profile already carries a fixed dimensionless depth parameter

`alpha_r = 10.`

If the physical radial wall amplitude is written as

`V0 = alpha_r mu_*`

relative to the local Thomas–Fermi enthalpy/chemical-potential scale `mu_*`, then

`Upsilon_w = 4 rho_w^2 V0^2 / (hbar^2 c_(s,w)^2)`
`          = alpha_r^2 Theta_w,`

where the only remaining microscopic amplitude datum is

`Theta_w := 4 rho_w^2 mu_*^2 / (hbar^2 c_(s,w)^2).`

So on the balanced Family-1 reference branch,

`Upsilon_w = 100 Theta_w.`

This is the sharpest explicit reduction of the wall branch so far: after fixing the geometry/support coordinates, the branch no longer depends on a free wall amplitude and a free wall depth separately. It depends on one dimensionless microscopic wall-depth datum `Theta_w`.

---

## 4. Explicit threshold window for `Theta_w`

Since `Upsilon_w = 100 Theta_w`, the explicit branch theorem becomes

`Theta_w <= Theta_fail`  -> fail,

`Theta_w >= Theta_suff`  -> succeed,

with

`Theta_fail = Upsilon_fail / 100`
`           ≈ 3.62605617972939e-4 * Pe_req,`

`Theta_suff = Upsilon_suff / 100`
`           ≈ 4.21495341569977e-2 * Pe_req.`

So the moving-throat placement problem is now no longer a three-parameter branch hunt.

It is one explicit microscopic question:

> does the actual Family-1 wall depth datum `Theta_w` on the true branch lie below, within, or above this window?

---

## 5. What Stage 58 changes

Before this stage, the explicit branch still looked as if it required solving for all three Stage-55 controls.

After this stage:

- the geometry ratio is fixed,
- the support/healing ratio is fixed,
- the support operator scales are fixed,
- and the only remaining branch datum is the microscopic wall-depth amplitude `Theta_w`.

That means the explicit-branch phase is now essentially finished.

What remains is one last microscopic closure question on this branch:
derive or estimate `Theta_w` from the real wall/throat PDE and compare it directly to the explicit threshold window above.

=== moving_throat_pde_stage059_n5_wall_depth_lock.md ===

# Moving-Throat PDE — Stage 59: Exact `n=5` Wall-Depth Lock for the Family-1 Branch

## Purpose

Stage 58 reduced the explicit Family-1 reference branch to one remaining microscopic datum,

`Theta_w := 4 rho_w^2 mu_*^2 / (hbar^2 c_(s,w)^2),`

through

`Upsilon_w = alpha_r^2 Theta_w,`

with the frozen reference value `alpha_r = 10`.

The next honest step is to stop treating `Theta_w` as an opaque constant.

The present stage does that by combining two already-carried exact identities:

1. the frozen `n=5` GNLS thermodynamic relation,
2. the Stage-57 healing-width lock.

The result is that on the reference branch the wall-depth datum is not an arbitrary amplitude. It is directly proportional to the square of the effective wall density on the active shell.

---

## 1. Exact `n=5` enthalpy–sound-speed identity

For the frozen GNLS EOS

`P(rho) = K rho^5,`

we have

`c_s^2(rho) = (1/m) dP/drho = (5K/m) rho^4,`

and

`U(rho) = (K/4) rho^5,`

`h(rho) = dU/drho = (5K/4) rho^4.`

Therefore the local enthalpy is exactly

`h(rho) = m c_s(rho)^2 / 4.`

So on the `n=5` branch, enthalpy and sound speed are not independent.

---

## 2. Local enthalpy lock for the Family-1 wall amplitude

Stage 58 wrote the radial wall depth as

`V0 = alpha_r mu_*`

relative to a local Thomas–Fermi enthalpy / chemical-potential scale `mu_*`.

The natural exact closure on the `n=5` branch is to write

`mu_* = lambda_mu h_w,`

where

`h_w := h(rho_w) = m c_(s,w)^2 / 4,`

and `lambda_mu` keeps track of whether one chooses the wall enthalpy itself (`lambda_mu = 1`) or a nearby local chemical-potential normalization.

Then

`Theta_w = 4 rho_w^2 mu_*^2 / (hbar^2 c_(s,w)^2)`
`        = lambda_mu^2 m^2 rho_w^2 c_(s,w)^2 / (4 hbar^2).`

So the only remaining microscopic input is now the effective wall density on the active shell.

---

## 3. Exact healing-lock reduction

Stage 57 already fixed the local healing-width closure

`ell = hbar / (2 m c_(s,w)).`

Insert this into the previous formula:

`m c_(s,w) / hbar = 1 / (2 ell),`

hence

`Theta_w = lambda_mu^2 rho_w^2 / (16 ell^2).`

So on the explicit moving-throat branch,

> the wall-depth datum is exactly the square of the active-shell density, measured in healing-width units, up to the one local normalization factor `lambda_mu`.

---

## 4. Reference-branch form

On the explicit Family-1 reference branch,

`ell / a = 1/20.`

In the dimensionless Family-1 wall coordinates already used by the coupled profile,

`x = r/a,`

so `a = 1` in the normalized wall-shape description and therefore

`1 / ell^2 = 400.`

Hence the branch-level datum becomes

`Theta_w = 25 lambda_mu^2 rho_w^2.`

This is the cleanest algebraic reduction of `Theta_w` so far.

---

## 5. What Stage 59 changes

Before this step, the explicit Family-1 branch still had one unresolved microscopic wall datum with no clean parent formula.

After this step, that datum is no longer opaque:

`Theta_w = 25 lambda_mu^2 rho_w^2`

in the normalized Family-1 wall variables.

So the only real remaining task on this branch is now to choose the correct effective wall density `rho_w` on the active shell and compare the resulting `Theta_w` to the explicit Stage-58 threshold window.

=== moving_throat_pde_stage060_family1_theta_extraction.md ===

# Moving-Throat PDE — Stage 60: Shell-Weighted Extraction of `Theta_w` on the Explicit Family-1 Wall

## Purpose

Stage 59 reduced the explicit Family-1 wall-depth datum to

`Theta_w = 25 lambda_mu^2 rho_w^2`

in normalized Family-1 wall variables.

So the only remaining question on this branch is what effective wall density `rho_w` should actually be used on the active shell.

This stage answers that by using the concrete Family-1 radial Thomas–Fermi wall profile already frozen in the coupled-wall appendix together with the canonical wall-support weight used throughout Stages 53–58.

---

## 1. Explicit midplane radial wall profile

On the balanced Family-1 reference branch,

`alpha_r = 10,`

`epsilon_r = 0.05,`

`p_r = 2.`

At the symmetry midplane the endcap term is exponentially negligible, so the coupled full-profile wall reduces to the radial slice

`rho_r(x) = [ 1 - alpha_r S((x-1)/epsilon_r)^2 ]_+^(1/4),`

with

`S(xi) = (1 + tanh xi)/2.`

In the local shell coordinate `xi = (x-1)/epsilon_r`, this becomes

`rho_r(xi) = [ 1 - alpha_r S(xi)^2 ]_+^(1/4).`

The support layer only overlaps the interior flank up to the exact cut point

`xi_* = artanh( 2/sqrt(alpha_r) - 1 ).`

For `alpha_r = 10`,

`xi_* ≈ -0.3855810692.`

So the active support weight lies mostly on the inner edge of the wall, not at the formal center of the soft switch.

---

## 2. Canonical support weight carried from the explicit branch

To stay on the same branch as Stages 53–58, keep the canonical support profile

`chi_phi(xi) = S'(xi) = (1/2) sech^2 xi.`

Its exact normalization moment is

`I_f = int dxi chi_phi(xi)^2 = 1/3.`

So the natural shell-weighted average of any wall quantity `Q(xi)` is

`<Q>_chi := [ int dxi chi_phi(xi)^2 Q(xi) ] / [ int dxi chi_phi(xi)^2 ].`

---

## 3. Why the correct effective wall datum uses `<rho^2>_chi`

Stage 59 showed that the local wall-depth datum scales as

`Theta(xi) = 25 lambda_mu^2 rho_r(xi)^2.`

So the natural effective branch datum is the support-weighted average of `Theta`, not `Theta` evaluated on a point and not `Theta` built from the averaged density after the square is taken.

Therefore

`Theta_w^(chi) = 25 lambda_mu^2 < rho_r^2 >_chi.`

This is the quantity that preserves the actual quadratic wall-depth weighting of the explicit branch.

---

## 4. Numerical extraction on the `alpha_r = 10` Family-1 branch

Using the exact canonical support weight and the explicit Family-1 radial profile gives

`<rho_r>_chi ≈ 0.192619005556493,`

`<rho_r^2>_chi ≈ 0.162745294003265.`

So the effective wall-depth datum is

`Theta_w^(chi) = 25 lambda_mu^2 <rho_r^2>_chi`
`             ≈ 4.06863235008162 lambda_mu^2.`

This is the first concrete branch value of `Theta_w` derived directly from an explicit wall-support profile rather than left symbolic.

---

## 5. Conservative lower-envelope estimator

For bookkeeping it is also useful to record the stricter Jensen-style lower-envelope obtained by averaging the density first and squaring afterward:

`Theta_w^(J) := 25 lambda_mu^2 <rho_r>_chi^2`
`            ≈ 0.927552032539308 lambda_mu^2.`

Because

`<rho_r^2>_chi >= <rho_r>_chi^2,`

this satisfies

`Theta_w^(chi) >= Theta_w^(J).`

So `Theta_w^(J)` can be used as a very conservative branch floor, while `Theta_w^(chi)` is the natural quadratic-energy matching value.

---

## 6. What Stage 60 changes

Before this step, the explicit Family-1 branch still had one unresolved microscopic wall-depth amplitude.

After this step, that datum is no longer abstract. On the canonical explicit branch,

`Theta_w^(chi) ≈ 4.06863235008162 lambda_mu^2`

with the conservative lower envelope

`Theta_w^(J) ≈ 0.927552032539308 lambda_mu^2.`

So the final branch-level question is no longer “what is `Theta_w`?”
It is now only:

> where do these explicit branch values sit relative to the exact Stage-58 threshold window?

=== moving_throat_pde_stage061_family1_branch_verdict.md ===

# Moving-Throat PDE — Stage 61: Explicit Family-1 Branch Comparison and Closing Verdict for This Subprogram

## Purpose

Stage 58 reduced the explicit Family-1 reference branch to the threshold window

`Theta_fail ≈ 3.62605617972939e-4 * Pe_req,`

`Theta_suff ≈ 4.21495341569977e-2 * Pe_req.`

Stage 60 then extracted the actual explicit branch datum

`Theta_w^(chi) ≈ 4.06863235008162 lambda_mu^2,`

with conservative lower envelope

`Theta_w^(J) ≈ 0.927552032539308 lambda_mu^2.`

This stage compares them directly and states the branch-level verdict.

---

## 1. Natural quadratic branch comparison

Use the natural support-weighted datum `Theta_w^(chi)`.

The exact explicit branch inequalities are

`Theta_w >= Theta_suff`  -> guaranteed success,

`Theta_w <= Theta_fail`  -> guaranteed failure.

Insert `Theta_w = Theta_w^(chi)`.

Then the explicit Family-1 branch is guaranteed to succeed whenever

`Pe_req <= Pe_suff^(chi) := Theta_w^(chi) / 4.21495341569977e-2`
`                     ≈ 96.5285247264386 lambda_mu^2,`

and it is guaranteed to fail only if

`Pe_req >= Pe_fail^(chi) := Theta_w^(chi) / 3.62605617972939e-4`
`                     ≈ 11220.5441626259 lambda_mu^2.`

So on the natural explicit wall-depth extraction, the reference Family-1 branch already clears the branch theorem unless the required quadrupole demand is anomalously large.

---

## 2. Conservative lower-envelope comparison

If one insists on the stricter lower-envelope estimator `Theta_w^(J)`, then guaranteed success still holds whenever

`Pe_req <= Pe_suff^(J) := Theta_w^(J) / 4.21495341569977e-2`
`                   ≈ 22.0062226330754 lambda_mu^2,`

while guaranteed failure would require

`Pe_req >= Pe_fail^(J) := Theta_w^(J) / 3.62605617972939e-4`
`                   ≈ 2558.01892349205 lambda_mu^2.`

So even the conservative floor says the explicit Family-1 branch is not obviously wall-depth starved unless the required quadrupole demand is already very large.

---

## 3. Branch-level verdict

The explicit Family-1 subprogram can now be closed cleanly.

The dominant result is:

> on the first concrete moving-throat branch, the wall-depth datum is no longer the natural bottleneck.

Within the natural `n=5` enthalpy lock and the explicit shell-weighted extraction, the branch supplies an `O(1)` wall-depth strength,

`Theta_w^(chi) ≈ 4.069 lambda_mu^2,`

while the Stage-58 success window only demands

`Theta_w >= 4.2149534e-2 Pe_req.`

So the branch succeeds automatically for modest required demand and fails only for very large demand.

That means the explicit-branch phase has reached its natural finish line:

- the geometry is fixed,
- the healing/support scale is fixed,
- the operator thresholds are fixed,
- the wall-depth datum has been extracted from a concrete parent-wall closure,
- and the remaining unresolved quantity is now the quadrupole-side demand `Pe_req`, not the wall-depth supply.

---

## 4. What remains open after Stage 61

This closes the Family-1 explicit-branch placement subprogram.

The remaining serious problem is now the one the higher-order stack had already isolated:

> determine the actual quadrupole-side requirement `Pe_req` from the completed moving-throat / outgoing-normalization branch and compare it to the explicit success/failure bands above.

So from this point on, the real bottleneck is no longer the wall-depth amplitude. It is the final quadrupole normalization bridge.

=== moving_throat_pde_stage062_family1_quadrupole_pe_map.md ===

# Moving-Throat PDE — Stage 62: Exact Family-1 Map from Quadrupole Demand `zeta_req` to the Required Transport Bias `Pe_req`

## Purpose

Stage 61 closed the explicit Family-1 wall-depth subprogram and showed that the natural branch is not wall-depth starved. The remaining serious bottleneck was therefore no longer the supply side `Theta_w`, but the demand side `Pe_req` coming from the quadrupole-normalization branch.

The next honest step is to eliminate the abstract `Pe_req` in favor of the explicit support-ratio demand `zeta_req` carried from Stages 35 and 42.

This stage does that for the concrete Family-1 branch.

The main result is that, on the explicit Family-1 geometry/support branch,

`zeta_F1(Pe) = A_F1 Omega_Pe^2`

with a fixed branch factor

`A_F1 = (kappa_F1 + pi^2/4)/(kappa_F1 + y_F1^2),`

where

`kappa_F1 = 12321/5,`

`y_F1 tan(y_F1) = 37,  0 < y_F1 < pi/2.`

So the remaining demand problem is now the unique inverse relation

`zeta_F1(Pe_req) = zeta_req.`

That already yields a hard Family-1 quadrupole-demand ceiling:

`zeta_req <= zeta_max^(F1) := A_F1 pi^2 / 4`

with

`zeta_max^(F1) ≈ 2.46752922945601.`

So even before the final moving-throat quadrupole normalization is solved, the explicit Family-1 branch can only support the selected quadrupole branch if its required support ratio stays below this concrete `O(1)` ceiling.

---

## 1. Fixed Family-1 support-compliance factor

Stages 56–58 fixed the explicit Family-1 branch to

`kappa_F1 = 12321/5,`

`eta_F1 = 37.`

Let `y_F1` denote the unique Robin root on the constructive branch,

`y_F1 tan(y_F1) = 37,`

`0 < y_F1 < pi/2.`

Then the exact Robin support factor carried from Stage 40 is

`A_F1 := (kappa_F1 + pi^2/4) / (kappa_F1 + y_F1^2).`

Numerically,

`y_F1 ≈ 1.52948248371470,`

`A_F1 ≈ 1.00005192880220.`

So the explicit Family-1 branch is already slightly above the symmetric-twin baseline even at zero transport bias.

---

## 2. Exact Family-1 support-ratio map

Stage 39 gave the exact constructive transport overlap factor

`Omega_Pe`
`= pi Pe (2 Pe exp(Pe) + pi)`
`  / [ (4 Pe^2 + pi^2) (exp(Pe)-1) ],`

with continuous extension `Omega_0 = 1`.

Therefore the explicit Family-1 lowest-lane support ratio is

`zeta_F1(Pe) = A_F1 Omega_Pe^2.`

This is the exact demand-side bridge that Stage 61 was still missing.

Its endpoint values are

`zeta_F1(0) = A_F1 ≈ 1.00005192880220,`

`lim_(Pe -> +infinity) zeta_F1(Pe) = A_F1 pi^2/4`
`                                  ≈ 2.46752922945601.`

So the full Family-1 constructive branch spans only the compact interval

`1.0000519288... <= zeta_F1(Pe) <= 2.4675292294... .`

That is the first exact explicit ceiling on the quadrupole-demand side.

---

## 3. Exact inversion problem for `Pe_req`

Given a selected quadrupole-branch demand `zeta_req`, the required transport bias is the unique constructive root of

`A_F1 Omega_(Pe_req)^2 = zeta_req.`

Equivalently,

`Omega_(Pe_req)^2 = zeta_req / A_F1.`

Because the carried constructive branch has `Omega_Pe` strictly increasing from `1` to `pi/2`, this inverse exists iff

`A_F1 <= zeta_req <= A_F1 pi^2/4.`

So for the explicit Family-1 branch:

- if `zeta_req < A_F1`, the demand is already met at zero transport bias;
- if `A_F1 <= zeta_req <= zeta_max^(F1)`, there is a unique constructive `Pe_req`;
- if `zeta_req > zeta_max^(F1)`, the Family-1 branch fails irrespective of wall depth.

This is a theorem-level sharpening of Stage 61: the wall-depth supply may be generous, but the branch still has a hard support-ratio ceiling.

---

## 4. Small-demand expansion

Using the carried Stage-39 expansion

`Omega_Pe = 1 + ((4-pi)/(2pi)) Pe + O(Pe^2),`

the Family-1 support-ratio map begins as

`zeta_F1(Pe)`
`= A_F1 [ 1 + ((4-pi)/pi) Pe + O(Pe^2) ].`

So for demands only slightly above the zero-bias baseline,

`zeta_req = A_F1 (1 + eps_z),`

the required transport bias is approximately

`Pe_req ~= (pi/(4-pi)) eps_z.`

This shows that the initial transport cost is linear in the excess quadrupole demand above the zero-bias Family-1 baseline.

---

## 5. What Stage 62 changes

Before this step, the remaining bottleneck was written as the abstract Peclet requirement `Pe_req`.

After this step, the explicit Family-1 branch no longer needs `Pe_req` as an independent concept.
It has an exact demand map

`zeta_req <-> Pe_req`

with a hard ceiling

`zeta_max^(F1) ≈ 2.46752922945601.`

So the remaining serious question is now much narrower:

> does the final selected quadrupole branch demand a support ratio `zeta_req` that stays below this explicit Family-1 ceiling?

=== moving_throat_pde_stage063_family1_zeta_thresholds.md ===

# Moving-Throat PDE — Stage 63: Explicit Family-1 Conversion of the Stage-61 `Pe_req` Window into Quadrupole-Demand Thresholds `zeta_req`

## Purpose

Stage 61 reduced the explicit Family-1 branch to a wall-depth comparison window in the transport-bias variable `Pe_req`:

`Pe_req <= Pe_suff^(chi)`  -> guaranteed success on the natural shell-weighted branch,

`Pe_req >= Pe_fail^(chi)`  -> guaranteed failure,

with the conservative lower-envelope pair `Pe_suff^(J)`, `Pe_fail^(J)` defined analogously.

Stage 62 then replaced the abstract `Pe_req` by the exact Family-1 support-ratio map

`zeta_F1(Pe) = A_F1 Omega_Pe^2.`

So the next honest step is to push the explicit Family-1 verdict entirely into the quadrupole-demand variable `zeta_req`.

This stage does that.

The main result is that the Family-1 branch now carries exact quadrupole-demand windows

`zeta_req <= zeta_suff^(chi)(lambda_mu)`  -> guaranteed success,

`zeta_req >= zeta_fail^(chi)(lambda_mu)`  -> guaranteed failure,

with corresponding conservative lower-envelope functions `zeta_suff^(J)`, `zeta_fail^(J)`.

For `lambda_mu = 1`, the natural Family-1 branch already reaches

`zeta_suff^(chi)(1) ≈ 2.46622291347846,`

while the hard Family-1 ceiling is only slightly larger,

`zeta_max^(F1) ≈ 2.46752922945601.`

So on the explicit branch the wall-depth supply is indeed not the dominant unresolved issue. The remaining serious question is whether the final selected quadrupole branch demands a support ratio above or below this narrow `O(2.46)` window.

---

## 1. Exact Stage-61 transport thresholds

Stage 61 gave the explicit Family-1 transport-bias windows

`Pe_suff^(chi) = 96.5285247264386 lambda_mu^2,`

`Pe_fail^(chi) = 11220.5441626259 lambda_mu^2,`

for the natural shell-weighted datum, and

`Pe_suff^(J) = 22.0062226330754 lambda_mu^2,`

`Pe_fail^(J) = 2558.01892349205 lambda_mu^2,`

for the conservative lower envelope.

---

## 2. Exact conversion to `zeta_req` thresholds

Using the explicit Family-1 demand map from Stage 62,

`zeta_F1(Pe) = A_F1 Omega_Pe^2,`

define

`zeta_suff^(chi)(lambda_mu) := zeta_F1(Pe_suff^(chi)),`

`zeta_fail^(chi)(lambda_mu) := zeta_F1(Pe_fail^(chi)),`

and similarly

`zeta_suff^(J)(lambda_mu) := zeta_F1(Pe_suff^(J)),`

`zeta_fail^(J)(lambda_mu) := zeta_F1(Pe_fail^(J)).`

Then the explicit Family-1 branch theorem becomes

`zeta_req <= zeta_suff^(chi)(lambda_mu)`  -> guaranteed success,

`zeta_req >= zeta_fail^(chi)(lambda_mu)`  -> guaranteed failure,

with the conservative version obtained by replacing `(chi)` with `(J)`.

So the Stage-61 wall-depth result is now written entirely in the quadrupole-demand language.

---

## 3. Numerical values at `lambda_mu = 1`

For the natural shell-weighted datum,

`zeta_suff^(chi)(1) ≈ 2.46622291347846,`

`zeta_fail^(chi)(1) ≈ 2.46752913273870.`

For the conservative lower envelope,

`zeta_suff^(J)(1) ≈ 2.44257571477179,`

`zeta_fail^(J)(1) ≈ 2.46752736855058.`

Compare these to the exact Family-1 ceiling

`zeta_max^(F1) ≈ 2.46752922945601.`

So on the natural explicit branch with `lambda_mu = 1`, the guaranteed-success threshold already lies less than `0.00131` below the hard ceiling, and the guaranteed-failure threshold is essentially saturated.

This is the sharpest numerical statement yet of the Stage-61 verdict.

---

## 4. Large-`lambda_mu` limit

Because all four Stage-61 transport thresholds scale as `lambda_mu^2` and `Omega_Pe` approaches `pi/2`, the corresponding quadrupole-demand thresholds satisfy

`lim_(lambda_mu -> +infinity) zeta_suff^(chi) = zeta_max^(F1),`

`lim_(lambda_mu -> +infinity) zeta_fail^(chi) = zeta_max^(F1),`

and likewise for the `(J)` branch.

So increasing the wall-depth normalization beyond `O(1)` does not open an unlimited quadrupole-demand window.
It only drives the branch toward the same hard Family-1 ceiling found in Stage 62.

This is another precise sense in which wall-depth supply is no longer the dominant open issue.

---

## 5. What Stage 63 changes

Before this step, the explicit Family-1 result still spoke in the transport-bias variable `Pe_req`.

After this step, the explicit branch verdict is fully phrased in the same support-ratio variable `zeta_req` that the quadrupole normalization branch actually demands.

So the remaining serious theorem question is now extremely narrow:

> does the final selected quadrupole branch require `zeta_req` below the explicit Family-1 support window `zeta_suff^(chi)(lambda_mu)` (or at least below the hard ceiling `zeta_max^(F1)`), or does it demand more than the explicit branch can ever supply?

=== moving_throat_pde_stage064_family1_pi_thresholds.md ===

# Moving-Throat PDE — Stage 64: Final Explicit Family-1 Quadrupole-Demand Window in the Branch-Product Variable `Pi_tr`

## Purpose

Stage 63 translated the explicit Family-1 wall-depth verdict into the support-ratio demand variable `zeta_req`.

The last honest step in this explicit-branch phase is to put the result back into the exact quadrupole branch-product language used by Stages 34–35.

This stage does that.

The main result is that, for the explicit Family-1 branch, the final selected quadrupole branch is constrained by exact thresholds of the form

`Pi_tr <= Pi_suff^(chi)(lambda_mu;eps_blk)`  -> guaranteed success,

`Pi_tr >= Pi_fail^(chi)(lambda_mu;eps_blk)`  -> guaranteed failure,

with a hard explicit Family-1 ceiling

`Pi_tr < Pi_max^(F1)(eps_blk)`

required for the branch to be reachable at all.

At vanishing blocking (`eps_blk = 0`) and with the natural wall normalization `lambda_mu = 1`, the explicit branch already gives

`Pi_suff^(chi) / C_mix ≈ 3.46622291347846,`

while the hard ceiling is only

`Pi_max^(F1) / C_mix ≈ 3.46752922945601.`

So the entire remaining explicit Family-1 theorem gap is now compressed to a very narrow `Pi_tr / C_mix` window.

---

## 1. Exact inversion of the Stage-35 support-demand law

Stage 35 gave

`zeta_req = (Pi_tr - C_mix) / [ C_mix - eps_blk (2 C_mix - Pi_tr) ],`

with `eps_blk` the blocking ratio carried from the selected quadrupole branch.

Solving this exactly for `Pi_tr` gives

`Pi_tr = C_mix Q(zeta_req;eps_blk),`

where

`Q(zeta;eps_blk)`
`:= [ 1 + (1 - 2 eps_blk) zeta ] / [ 1 - eps_blk zeta ].`

This map is exact.

Its anchor values are

`Q(0;eps_blk) = 1,`

`Q(1;eps_blk) = 2,`

and its derivative is

`dQ/dzeta = (1 - eps_blk) / (1 - eps_blk zeta)^2 > 0`

whenever `eps_blk zeta < 1`.

So the product demand grows strictly with the support-ratio demand.

---

## 2. Explicit Family-1 success and failure thresholds in `Pi_tr`

Insert the Stage-63 Family-1 thresholds.

Define

`Pi_suff^(chi)(lambda_mu;eps_blk)`
`:= C_mix Q( zeta_suff^(chi)(lambda_mu); eps_blk ),`

`Pi_fail^(chi)(lambda_mu;eps_blk)`
`:= C_mix Q( zeta_fail^(chi)(lambda_mu); eps_blk ),`

and similarly for the conservative floor,

`Pi_suff^(J) := C_mix Q( zeta_suff^(J); eps_blk ),`

`Pi_fail^(J) := C_mix Q( zeta_fail^(J); eps_blk ).`

Then the explicit Family-1 branch theorem is

`Pi_tr <= Pi_suff^(chi)`  -> guaranteed success,

`Pi_tr >= Pi_fail^(chi)`  -> guaranteed failure,

with the conservative version obtained by replacing `(chi)` with `(J)`.

So the Family-1 explicit-branch result is now written directly in the same branch-product variable that the quadrupole normalization program actually uses.

---

## 3. Hard explicit Family-1 ceiling in product form

Stage 62 gave the hard support-ratio ceiling

`zeta_req <= zeta_max^(F1) ≈ 2.46752922945601.`

Therefore the explicit Family-1 branch can only be reached if

`Pi_tr < Pi_max^(F1)(eps_blk)`

with

`Pi_max^(F1)(eps_blk)`
`:= C_mix Q( zeta_max^(F1); eps_blk )`
` = C_mix [ 1 + (1 - 2 eps_blk) zeta_max^(F1) ] / [ 1 - eps_blk zeta_max^(F1) ].`

This requires the denominator to remain positive, i.e.

`eps_blk < 1 / zeta_max^(F1) ≈ 0.405263689711371.`

So even before the final PDE-side quadrupole normalization is solved, the explicit Family-1 branch has an exact finite product ceiling.

---

## 4. Numerical illustration at `eps_blk = 0`

If the final quadrupole branch stays close to the unblocked limit `eps_blk = 0`, then

`Q(zeta;0) = 1 + zeta.`

So at `lambda_mu = 1` the explicit thresholds become

`Pi_suff^(chi) / C_mix = 1 + zeta_suff^(chi)(1)`
`                      ≈ 3.46622291347846,`

`Pi_fail^(chi) / C_mix = 1 + zeta_fail^(chi)(1)`
`                      ≈ 3.46752913273870,`

while the hard Family-1 ceiling is

`Pi_max^(F1) / C_mix = 1 + zeta_max^(F1)`
`                    ≈ 3.46752922945601.`

For the conservative lower envelope,

`Pi_suff^(J) / C_mix ≈ 3.44257571477179,`

`Pi_fail^(J) / C_mix ≈ 3.46752736855058.`

So the natural shell-weighted Family-1 branch at `lambda_mu = 1` already pushes the explicit product threshold essentially all the way to the hard ceiling.

---

## 5. Final explicit verdict for this phase

The explicit Family-1 branch now has a completely closed reduced theorem statement.

- Stage 61 showed that wall-depth supply is not the dominant open issue.
- Stage 62 converted the remaining bottleneck into a hard support-ratio ceiling `zeta_max^(F1)`.
- Stage 63 translated the Stage-61 wall verdict into the demanded support-ratio variable `zeta_req`.
- This stage finally expresses the result directly in the selected-branch product variable `Pi_tr`.

So the remaining explicit theorem gap is now as narrow as it can be without solving the final outgoing quadrupole normalization branch itself:

> does the completed moving-throat quadrupole branch produce a selected product `Pi_tr` that stays below the explicit Family-1 ceiling `Pi_max^(F1)(eps_blk)` and, more sharply, below the natural success threshold `Pi_suff^(chi)(lambda_mu;eps_blk)`?

That is the clean finish line for the present explicit-branch phase.

=== moving_throat_pde_stage065_master_quadrupole_residual.md ===

# Moving-Throat PDE — Stage 65: Master Quadrupole Residual of the Full Reduced Moving-Throat PDE

## Purpose

The explicit Family-1 branch phase ended with a hard support ceiling in the selected-branch product variable `Pi_tr`, but that result was still being carried through several intermediate variables (`Pe_req`, `zeta_req`, `Pi_tr/C_mix`).

The next honest step toward a full moving-throat PDE write-up is to state the entire remaining quadrupole-normalization problem as **one reduced PDE residual** and then show how the earlier exact threshold theorems sit inside it.

That is what this stage does.

The main result is that the whole reduced moving-throat PDE, on the surviving passive/outgoing quadrupole branch, is now summarized by one scalar residual

`R_quad(Pi_tr,C_mix,eps_blk ; Xi,eta,kappa)`
`:= zeta_req(Pi_tr,C_mix,eps_blk) - zeta_phys(Pe_*(Xi,eta,kappa),eta;kappa),`

where

`zeta_req(Pi_tr,C_mix,eps_blk)`
`= (Pi_tr - C_mix) / [ C_mix - eps_blk (2 C_mix - Pi_tr) ]`

is the exact selected-branch demand from Stages 34–35,

and

`zeta_phys(Pe,eta;kappa)`
`= Omega_Pe^2 (kappa + pi^2/4) / (kappa + y(eta)^2),`

`y tan y = eta,`

is the exact explicit lowest-lane support ratio from Stages 39–40,
with `Pe_*` the operator-selected transport bias solving the Stage-41 fixed-point equation

`Pe_* = Xi Delta(Pe_*;kappa,eta).`

So the remaining reduced theorem question is no longer diffuse:

> **does the completed moving-throat quadrupole branch make `R_quad` nonpositive on the physical branch?**

---

## 1. Full reduced system entering the quadrupole residual

The reduced moving-throat PDE hierarchy now has five live layers.

### 1.1 Parent bulk matter/gauge sector

The exact parent theory still carries

`i hbar D_t psi`
`= [ -(hbar^2/2m) D_A D_A + V_conf(X;a,L) + h(|psi|^2) ] psi,`

and the localized `4+1` Maxwell sector

`partial_M ( Z(w) F^(MN) ) + (1/xi) partial^N(partial·A) = mu_0 J^N,`

with the geometry variables `(a,L)` or their distributed lift already embedded in the confinement sector.

### 1.2 Lowest support/source throat operator

On the reduced axial support/source branch, the same moving-throat program already isolated the coupled operator

`partial_t sigma + partial_s J = 0,`

`J = -D_sigma partial_s sigma + v_sigma sigma,`

`-T_X partial_s^2 phi + K_X phi = Lambda_phi sigma,`

with Robin/Neumann support boundaries

`T_X phi_s(0) = K_m phi(0),`

`phi_s(L) = 0.`

### 1.3 Stationary operator-selected transport bias

On the stationary zero-flux branch, the source density is

`Sigma_Pe(x) = Pe exp(Pe x)/(exp(Pe)-1),`

with `x=s/L in [0,1]`, and the support/source operator selects the transport bias through the exact fixed-point law

`Pe_* = Xi Delta(Pe_*;kappa,eta),`

where

`kappa = K_X L^2 / T_X,`

`eta = K_m L / T_X,`

`Xi = mu_sigma Lambda_phi^2 L^2 / (D_sigma T_X).`

### 1.4 Exact explicit lowest-lane support ratio

Once the support/source branch is selected, the explicit lowest support lane carries the exact ratio

`zeta_phys(Pe,eta;kappa)`
`= Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2),`

with the transport overlap factor

`Omega_Pe`
`= pi Pe (2 Pe exp(Pe) + pi)`
`  / [ (4 Pe^2 + pi^2) (exp(Pe)-1) ]`

and `y(eta)` the unique Robin root `y tan y = eta`, `0<y<pi/2`.

### 1.5 Exact selected-branch quadrupole demand

On the outgoing quadrupole side, the selected-branch support demand remains

`zeta_req(Pi_tr,C_mix,eps_blk)`
`= (Pi_tr - C_mix) / [ C_mix - eps_blk (2 C_mix - Pi_tr) ],`

which is just Stage 35 rewritten as the exact demand variable of the branch-product formulation.

So the reduced moving-throat PDE closes to one scalar comparison:

`zeta_req ? zeta_phys`.

---

## 2. The master quadrupole residual

Define the exact residual

`R_quad(Pi_tr,C_mix,eps_blk ; Xi,eta,kappa)`
`:= zeta_req(Pi_tr,C_mix,eps_blk)`
` - zeta_phys(Pe_*(Xi,eta,kappa),eta;kappa).`

Then the selected outgoing quadrupole branch satisfies:

- `R_quad < 0`  ->  support supply exceeds demand;
- `R_quad = 0`  ->  exact saturation of the surviving quadrupole branch;
- `R_quad > 0`  ->  the explicit lowest-lane support branch fails.

This is the cleanest full reduced-PDE statement reached so far.

It absorbs all of the earlier stage variables into one exact scalar diagnostic.

---

## 3. Exact bounded version using the support/source operator brackets

Stage 41 already proved the exact operator-selected bracket

`Xi Delta_0(kappa,eta) <= Pe_* <= Xi Delta_inf(kappa,eta),`

with

`Delta_0(kappa,eta)`
`= eta (cosh(alpha)-1) / [ alpha^2 ( alpha sinh(alpha) + eta cosh(alpha) ) ],`

`Delta_inf(kappa,eta)`
`= [ cosh(alpha) + (eta/alpha) sinh(alpha) - 1 ]`
`  / [ alpha sinh(alpha) + eta cosh(alpha) ],`

`alpha = sqrt(kappa).`

Since `zeta_phys(Pe,eta;kappa)` is strictly increasing in `Pe`, the physical support ratio obeys the exact bounds

`zeta_-(Xi,eta;kappa)`
`:= zeta_phys( Xi Delta_0(kappa,eta), eta; kappa )`
`<= zeta_phys(Pe_*,eta;kappa)`
`<= zeta_phys( Xi Delta_inf(kappa,eta), eta; kappa )`
`: = zeta_+(Xi,eta;kappa).`

Therefore the exact residual is trapped between

`R_-(...) := zeta_req - zeta_+ <= R_quad <= zeta_req - zeta_- =: R_+(...).`

So the reduced PDE already has theorem-level success/failure bounds:

- if `R_+ <= 0`, success is guaranteed;
- if `R_- > 0`, failure is guaranteed;
- only the strip `R_- <= 0 < R_+` requires the full fixed-point solve.

---

## 4. Direct thresholds in the branch-product variable

Because `zeta_req` is strictly increasing in `Pi_tr`, the bounded residual can be translated back into exact product thresholds.

Write

`Q(zeta;eps_blk)`
`:= [ 1 + (1 - 2 eps_blk) zeta ] / [ 1 - eps_blk zeta ].`

Then the inverse map is

`Pi_tr = C_mix Q(zeta_req;eps_blk).`

So the bounded reduced-PDE theorem is:

`Pi_tr <= Pi_suff(Xi,eta,kappa ; C_mix,eps_blk)`  -> guaranteed success,

`Pi_tr >= Pi_fail(Xi,eta,kappa ; C_mix,eps_blk)`  -> guaranteed failure,

with

`Pi_suff := C_mix Q( zeta_-(Xi,eta;kappa) ; eps_blk ),`

`Pi_fail := C_mix Q( zeta_+(Xi,eta;kappa) ; eps_blk ).`

This is the first exact product-window theorem derived directly from the coupled support/source operator rather than through the intermediate `Pe_req` bookkeeping.

---

## 5. Family-1 specialization of the master residual

The explicit Family-1 branch now drops out immediately by inserting the fixed geometry/support data

`kappa_F1 = 12321/5,`

`eta_F1 = 37,`

and the wall/source strength relation

`Xi_F1 = W_wall = Upsilon_w Lambda_ell^2 = 1369 Upsilon_w = 136900 Theta_w.`

So the Family-1 residual is simply

`R_quad^(F1)`
`= zeta_req(Pi_tr,C_mix,eps_blk)`
` - zeta_phys(Pe_*(Xi_F1,37,12321/5), 37; 12321/5).`

That is exactly the residual the later explicit branch stages had been approaching indirectly.

---

## 6. Why Stage 65 matters

This stage is the first place where the whole reduced moving-throat PDE can be written in one line without hiding the actual theorem gap.

The remaining unresolved problem is not another support/source reduction and not another wall-depth estimate. It is the sign of one scalar residual:

`R_quad(Pi_tr,C_mix,eps_blk ; Xi,eta,kappa).`

Everything else in the reduced program now feeds this one object.

=== moving_throat_pde_stage066_family1_direct_operator_window.md ===

# Moving-Throat PDE — Stage 66: Direct Operator-Selected Family-1 Window for the Surviving Quadrupole Branch

## Purpose

Stage 65 compressed the whole remaining reduced moving-throat PDE into one master residual

`R_quad = zeta_req - zeta_phys(Pe_*),`

with `Pe_*` selected by the exact support/source fixed-point equation.

The next honest step is to evaluate that master residual on the first explicit branch itself, rather than continuing to carry the Family-1 window through intermediate handoffs.

This stage does that.

The main results are:

1. the operator-selected transport bias is monotone in the wall/source strength `Xi`,
2. therefore the explicit Family-1 support ratio and the corresponding branch-product thresholds are also monotone in `Xi`,
3. inserting the natural Family-1 shell-weighted and Jensen wall data reproduces the Stage-61/63/64 windows directly from the coupled operator,
4. and the natural branch already lies extremely close to the hard Family-1 ceiling.

So the explicit Family-1 branch is now fully closed at the operator-selected level.

---

## 1. Exact monotonicity of the operator-selected branch

The stationary branch is selected by

`Pe_* = Xi Delta(Pe_*;kappa,eta).`

Differentiating implicitly gives

`dPe_*/dXi = Delta(Pe_*;kappa,eta) / [ 1 - Xi partial_Pe Delta(Pe_*;kappa,eta) ].`

On the constructive stable branch,

`Delta > 0,`

`partial_Pe Delta > 0,`

and the fixed-point branch remains stable only while

`1 - Xi partial_Pe Delta > 0.`

Therefore

`dPe_*/dXi > 0.`

Because the explicit lowest-lane support ratio

`zeta_phys(Pe,eta;kappa)`
`= Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2)`

is strictly increasing in `Pe`, the operator-selected support ratio is strictly increasing in `Xi`.

So stronger wall/source gain always pushes the explicit Family-1 branch upward toward its hard ceiling.

---

## 2. Exact Family-1 operator data

For the explicit Family-1 branch,

`kappa_F1 = 12321/5,`

`eta_F1 = 37,`

and the exact support/source strength reduces to

`Xi_F1 = W_wall = 1369 Upsilon_w = 136900 Theta_w.`

Using the exact Stage-41 support/source brackets,

`Pe_-^(F1) = Xi_F1 Delta_0(kappa_F1,eta_F1),`

`Pe_+^(F1) = Xi_F1 Delta_inf(kappa_F1,eta_F1),`

with

`Delta_0(kappa_F1,eta_F1) ≈ 1.73302079021525e-4,`

`Delta_inf(kappa_F1,eta_F1) ≈ 2.01447565540522e-2.`

So the explicit operator-selected Family-1 support ratio lies between

`zeta_-^(F1)(Xi_F1) = zeta_phys(Pe_-^(F1),37;12321/5),`

and

`zeta_+^(F1)(Xi_F1) = zeta_phys(Pe_+^(F1),37;12321/5).`

---

## 3. Natural shell-weighted and conservative floor windows

Stage 60 gave the exact wall-depth data

`Theta_w^(chi) ≈ 4.06863235008162 lambda_mu^2,`

`Theta_w^(J)   ≈ 0.927552032539308 lambda_mu^2.`

Therefore the operator-selected Family-1 strengths are

`Xi_F1^(chi) ≈ 556995.768726174 lambda_mu^2,`

`Xi_F1^(J)   ≈ 126981.873254631 lambda_mu^2.`

At `lambda_mu = 1`, the exact bracketed transport windows become

`Pe_-^(chi) ≈ 96.5285247264385,`

`Pe_+^(chi) ≈ 11220.5441626259,`

`Pe_-^(J)   ≈ 22.0062226330754,`

`Pe_+^(J)   ≈ 2558.01892349205.`

These are exactly the Stage-61 transport windows, now derived directly from the master support/source operator.

---

## 4. Direct operator-selected support-ratio windows

Inserting those operator-selected transport bounds into the exact Family-1 support map gives

`zeta_-^(chi)(1) ≈ 2.46622291347846,`

`zeta_+^(chi)(1) ≈ 2.46752913273870,`

`zeta_-^(J)(1)   ≈ 2.44257571477179,`

`zeta_+^(J)(1)   ≈ 2.46752736855058.`

Compare these with the hard Family-1 ceiling

`zeta_max^(F1) ≈ 2.46752922945601.`

So on the natural shell-weighted branch the guaranteed-success lower window already sits less than about `1.31e-3` below the hard ceiling, while the guaranteed-failure upper window differs from the ceiling only in the seventh decimal place.

This is the operator-selected version of the earlier explicit-branch verdict:

> the Family-1 branch is not support-starved; it already drives the explicit support ratio essentially to its maximal constructive value.

---

## 5. Direct operator-selected branch-product windows

Using the exact inverse map

`Pi = C_mix Q(zeta;eps_blk),`

`Q(zeta;eps_blk) = [1 + (1 - 2 eps_blk) zeta] / [1 - eps_blk zeta],`

these operator-selected support windows become direct selected-branch product windows

`Pi_suff^(F1) = C_mix Q( zeta_-^(F1); eps_blk ),`

`Pi_fail^(F1) = C_mix Q( zeta_+^(F1); eps_blk ).`

At vanishing blocking `eps_blk = 0`, these reduce simply to

`Pi_suff^(F1) / C_mix = 1 + zeta_-^(F1),`

`Pi_fail^(F1) / C_mix = 1 + zeta_+^(F1).`

So for `lambda_mu = 1`,

`Pi_suff^(chi)/C_mix ≈ 3.46622291347846,`

`Pi_fail^(chi)/C_mix ≈ 3.46752913273870,`

`Pi_suff^(J)/C_mix   ≈ 3.44257571477179,`

`Pi_fail^(J)/C_mix   ≈ 3.46752736855058,`

while the hard ceiling remains

`Pi_max^(F1)/C_mix ≈ 3.46752922945601.`

So the explicit operator-selected Family-1 window in the actual quadrupole branch variable is now fully closed.

---

## 6. What Stage 66 changes

This stage removes the last remaining ambiguity in the explicit Family-1 branch phase.

The earlier stages had already shown that support depth and source asymmetry could in principle reach the required window. What was still missing was a direct operator-selected statement.

Now we have it.

The coupled support/source operator itself selects a branch window that is already essentially saturated at the Family-1 ceiling once the natural shell-weighted wall datum is inserted.

So the remaining reduced theorem gap is no longer on the explicit support/source side at all.
It is entirely on the outgoing quadrupole-normalization side:

> what value of `Pi_tr` does the actual passive/outgoing moving-throat quadrupole branch produce, and does it stay below the explicit Family-1 ceiling? 

=== moving_throat_pde_stage067_full_reduced_pde_writeup.md ===

# Moving-Throat PDE — Stage 67: Full Reduced Moving-Throat PDE Write-Up Skeleton and the Final Remaining Theorem Gap

## Purpose

Stages 1–66 have now done two different kinds of work:

1. they built the reduced moving-throat PDE hierarchy layer by layer, and
2. they compressed the remaining support/source branch question to one explicit quadrupole residual.

So this stage is not another isolated calculation. It is the first honest **write-up skeleton** of the full reduced moving-throat PDE program as it currently stands.

The point is to say clearly what we now have, what is already exact or SymPy-backed, and what the actual remaining theorem gap is.

---

## 1. The reduced moving-throat PDE hierarchy now in hand

### 1.1 Exact parent bulk system

The parent theory already carries:

- the gauged 4D GNLS matter equation,
- the localized `4+1` Maxwell equation,
- the geometry variables `(a,L)` or their distributed wall lift,
- and the open-system projection/leakage machinery.

That is the exact parent backbone inherited from the action-based 4D program.

### 1.2 Linearized support/source throat operator

The explicit reduced throat-support sector is now an actual PDE/ODE operator, not just a symbolic closure:

`partial_t sigma + partial_s J = 0,`

`J = -D_sigma partial_s sigma + v_sigma sigma,`

`-T_X partial_s^2 phi + K_X phi = Lambda_phi sigma,`

with Robin/Neumann boundaries and the operator-selected fixed-point law

`Pe_* = Xi Delta(Pe_*;kappa,eta).`

### 1.3 Exact explicit lowest-lane support family

The explicit lowest support lane is now fully reduced to physical variables

`(Pe, eta, kappa)`

through

`zeta_phys(Pe,eta;kappa)`
`= Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2),`

`y tan y = eta.`

This family is no longer heuristic. Its monotonicity, closure window, overlap ceiling, Robin softening law, and exact operator-selected bounds are all SymPy-backed.

### 1.4 Exact selected-branch quadrupole demand map

On the outgoing quadrupole side, the selected-branch demand is now compressed to

`zeta_req(Pi_tr,C_mix,eps_blk)`
`= (Pi_tr - C_mix) / [ C_mix - eps_blk (2 C_mix - Pi_tr) ]`

with exact inversion

`Pi_tr = C_mix Q(zeta_req;eps_blk),`

`Q(zeta;eps_blk) = [1 + (1 - 2 eps_blk) zeta]/[1 - eps_blk zeta].`

So the outgoing branch has also been reduced to one exact scalar demand function.

---

## 2. The full reduced theorem gate

Putting those two sides together, the reduced moving-throat PDE is now governed by one scalar residual:

`R_quad(Pi_tr,C_mix,eps_blk ; Xi,eta,kappa)`
`= zeta_req(Pi_tr,C_mix,eps_blk)`
` - zeta_phys(Pe_*(Xi,eta,kappa),eta;kappa).`

This is the current write-up version of the full reduced theorem gate.

- `R_quad < 0`  ->  explicit support/source supply exceeds selected quadrupole demand;
- `R_quad = 0`  ->  exact saturation of the explicit branch;
- `R_quad > 0`  ->  the explicit branch fails.

And because the support/source fixed-point already obeys exact brackets,

`Xi Delta_0 <= Pe_* <= Xi Delta_inf,`

the reduced PDE also has exact success/failure bounds before any full root solve.

---

## 3. The explicit Family-1 specialization is now closed

The first real branch is now no longer abstract.
On Family-1 we have

`kappa_F1 = 12321/5,`

`eta_F1 = 37,`

`Xi_F1 = 1369 Upsilon_w = 136900 Theta_w.`

Using the exact wall-depth extraction from the frozen `n=5` branch and the carried Family-1 wall data, the natural shell-weighted and conservative-floor branches produce direct operator-selected windows in:

- transport bias `Pe`,
- support ratio `zeta`,
- and selected branch product `Pi_tr/C_mix`.

Numerically, the natural shell-weighted branch at `lambda_mu = 1` already pushes the explicit support window to

`zeta_-^(chi) ≈ 2.46622291347846,`

`zeta_+^(chi) ≈ 2.46752913273870,`

against the hard ceiling

`zeta_max^(F1) ≈ 2.46752922945601.`

So the explicit Family-1 support/source branch is effectively saturated.

That is why the support/source side is no longer the serious unresolved issue.

---

## 4. What is still not solved

Even after all this reduction, the project still does **not** yet have a first-principles theorem of the full moving-throat PDE.

What remains open is narrower:

1. the actual passive/outgoing quadrupole branch still has to produce its physical selected product `Pi_tr`,
2. the physical mixed baseline `C_mix` and blocking ratio `eps_blk` still have to be fixed on the true moving-throat branch,
3. and that final branch point must then be inserted into `R_quad`.

So the unresolved part is no longer the support/source operator, no longer the wall-depth supply, and no longer the explicit Family-1 reduction.

It is the **outgoing quadrupole-normalization branch itself**.

---

## 5. Honest current status

So we are indeed very close to the point where a full moving-throat PDE write-up becomes natural — but in the precise sense below.

### What is now ready to write up

A full **reduced moving-throat PDE write-up** is now ready:

- exact parent system,
- exact reduced support/source operator,
- exact explicit support family,
- exact selected-branch demand map,
- exact master residual `R_quad`,
- exact Family-1 specialization,
- and exact success/failure bounds.

### What is not yet honest to write up as finished

A full **first-principles theorem of the complete moving-throat quadrupole branch** is not yet ready, because the final outgoing-normalization branch is still the one unresolved ingredient.

So the correct statement is:

> the reduced moving-throat PDE program is now fully organized and almost fully written; the one remaining theorem gap is the actual passive/outgoing quadrupole-normalization branch that fixes `Pi_tr` on the true moving-throat solution.

---

## 6. The clean next move from here

The next derivation is now as narrow as it can reasonably get:

> derive the physical outgoing quadrupole product `Pi_tr` (and its accompanying `C_mix`, `eps_blk`) from the actual passive/outgoing moving-throat branch and evaluate the sign of `R_quad`.

That is the remaining finish line for the present reduced moving-throat PDE program.

=== moving_throat_pde_stage068_quadrupole_demand_cancellation.md ===

# Moving-Throat PDE — Stage 68: Exact Cancellation of the Outgoing-Normalization Factors in the Selected Quadrupole-Demand Product

## Purpose

Stage 67 left the reduced moving-throat PDE in the form of one master residual

`R_quad = zeta_req(Pi_tr,C_mix,eps_blk) - zeta_phys(Pe_*),`

but the selected branch still looked as if it depended on several separate outgoing-normalization objects at once:

- the target quadrupole normalization `N_Q^(target)`,
- the selected static prefactor `P_{0,-}`,
- the transfer factor `beta_0`,
- the selected overlap `s_-`,
- and the selected conservative eigenvalue `lambda_-`.

The next honest step is to ask whether those are really independent in the final support test.

This stage shows that they are not.

The main result is the exact pair of identities

`Pi_tr = (N_Q^(target)/beta_0) alpha_req,`

`C_mix = (N_Q^(target)/beta_0) alpha_mix,`

and therefore

`Pi_tr / C_mix = alpha_req / alpha_mix.`

Using the selected-mode normalization relation from Stage 13,

`N_Q^(target) = mhat_-^2 beta_0 s_- / lambda_-,`

the same identities become

`Pi_tr = mhat_-^2 (s_- / lambda_-) alpha_req,`

`C_mix = mhat_-^2 (s_- / lambda_-) alpha_mix.`

So once the outgoing quadrupole branch is required to hit the `2G/(5 c^5)` normalization target, the explicit support theorem no longer depends separately on `(mhat_-, beta_0, s_-, lambda_-)`.
It depends only on the **loading ratio** `alpha_req / alpha_mix`.

That is the cleanest bridge yet between the Stage-13 selected quadrupole normalization stack and the explicit Family-1 support/source branch.

---

## 1. Exact product identities

From Stages 18–19, the selected-branch target ratio is

`R_target = N_Q^(target) A / ( beta_0 kappa_0^2 ),`

with the exact D/N constant

`kappa_0^2 = 8 / pi^2.`

The exact dimensionless total and mixed loadings are

`G_tr = 8 alpha_req / (pi^2 A),`

`M_mix = 8 alpha_mix / (pi^2 A).`

By construction,

`Pi_tr = R_target G_tr,`

`C_mix = R_target M_mix.`

Inserting the formulas above and using `kappa_0^2 = 8/pi^2` gives the exact cancellations

`Pi_tr = (N_Q^(target)/beta_0) alpha_req,`

`C_mix = (N_Q^(target)/beta_0) alpha_mix.`

So the selected demand ratio is simply

`Pi_tr / C_mix = alpha_req / alpha_mix.`

This identity is exact.

---

## 2. Spectral form using the selected-mode normalization stack

Stage 13 gave the selected-mode outgoing normalization relation

`mhat_-^2 P_{0,-} = N_Q^(target),`

with

`P_{0,-} = beta_0 s_- / lambda_-`.

Therefore

`N_Q^(target) / beta_0 = mhat_-^2 s_- / lambda_-`.

Substituting this into the exact product identities yields the spectral forms

`Pi_tr = mhat_-^2 (s_- / lambda_-) alpha_req,`

`C_mix = mhat_-^2 (s_- / lambda_-) alpha_mix.`

So the same common spectral factor multiplies both the demanded product and the mixed baseline.

That is why the explicit support test loses all separate dependence on the outgoing normalization amplitudes once the branch is normalized.

---

## 3. Exact selected-demand law in pure loading-ratio form

Stage 35 gave the exact selected support demand

`zeta_req = (Pi_tr - C_mix) / [ C_mix - eps_blk (2 C_mix - Pi_tr) ].`

Using the cancellation identities above, the entire right-hand side reduces to

`zeta_req`
`= (alpha_req - alpha_mix)`
`  / [ alpha_mix - eps_blk (2 alpha_mix - alpha_req) ].`

Equivalently,

`zeta_req`
`= [ alpha_req/alpha_mix - 1 ]`
`  / [ 1 - eps_blk ( 2 - alpha_req/alpha_mix ) ].`

So the explicit support demand is now a function only of the loading ratio

`rho_alpha := alpha_req / alpha_mix`.

In the unblocked limit `eps_blk = 0`, this collapses to the especially simple law

`zeta_req = rho_alpha - 1.`

So on the unblocked branch the support-ratio demand is literally just the selected total loading divided by the mixed baseline, minus one.

---

## 4. What this changes conceptually

This stage removes one layer of apparent complexity from the remaining theorem gap.

Before this step, the outgoing quadrupole branch seemed to require us to know separately

- the selected conservative stiffness `lambda_-`,
- the overlap factor `s_-`,
- the transfer factor `beta_0`,
- the source map `mhat_-`,
- the demanded product `Pi_tr`,
- and the mixed baseline `C_mix`.

After the exact cancellations, the explicit support theorem depends on those objects only through the single ratio

`rho_alpha = alpha_req / alpha_mix`.

So the moving-throat support/source side is not being asked to solve the full outgoing normalization problem independently. Once the normalized quadrupole branch is imposed, it is being asked only one thing:

> **how much larger is the selected total directional loading than the mixed baseline?**

That is the right direct bridge between the selected quadrupole branch and the explicit Family-1 support/source program.

---

## 5. Best current theorem statement after Stage 68

Once the outgoing quadrupole branch satisfies the selected-mode normalization target,

`mhat_-^2 Gamma_{5,-} = 2 G / (5 c^5),`

the explicit support-demand problem is equivalent to the single ratio test

`rho_alpha = alpha_req / alpha_mix`.

More precisely,

`zeta_req = zeta_req(rho_alpha, eps_blk)`

with

`zeta_req(rho_alpha, eps_blk)`
`= (rho_alpha - 1) / [ 1 - eps_blk (2 - rho_alpha) ].`

So the remaining support theorem gap is no longer “determine several separate outgoing coefficients.”
It is:

> **determine the normalized selected-branch loading ratio `rho_alpha` of the actual passive/outgoing moving-throat quadrupole branch.**

=== moving_throat_pde_stage069_family1_loading_ratio_window.md ===

# Moving-Throat PDE — Stage 69: Exact Family-1 Success/Failure Window in the Pure Loading-Ratio Variable

## Purpose

Stage 68 showed that, once the selected outgoing quadrupole branch is normalized, the explicit support-demand problem depends only on the loading ratio

`rho_alpha := alpha_req / alpha_mix`.

The next honest step is to rewrite the entire explicit Family-1 theorem directly in that variable.

This stage does that.

The main result is that the explicit Family-1 branch no longer needs the product variables `(Pi_tr, C_mix)` at all.
It is governed exactly by the ratio window

`rho_alpha <= rho_suff^(chi)(lambda_mu;eps_blk)`  -> guaranteed success,

`rho_alpha >= rho_fail^(chi)(lambda_mu;eps_blk)`  -> guaranteed failure,

with a hard constructive ceiling

`rho_alpha < rho_max^(F1)(eps_blk).`

At `eps_blk = 0` and `lambda_mu = 1`, the natural shell-weighted Family-1 window is simply

`rho_alpha <= 3.46622291347846`  -> guaranteed success,

`rho_alpha >= 3.46752913273870`  -> guaranteed failure,

with the hard ceiling

`rho_alpha < 3.46752922945601.`

So the whole explicit Family-1 support theorem has collapsed to a very narrow pure loading-ratio window.

---

## 1. Exact ratio form of the Stage-64 product map

Stage 64 gave the exact product inversion

`Pi_tr = C_mix Q(zeta;eps_blk),`

with

`Q(zeta;eps_blk)`
`= [ 1 + (1 - 2 eps_blk) zeta ] / [ 1 - eps_blk zeta ].`

By Stage 68,

`Pi_tr / C_mix = rho_alpha = alpha_req / alpha_mix.`

Therefore the explicit Family-1 thresholds are immediately

`rho_suff^(chi)(lambda_mu;eps_blk)`
`:= Q( zeta_suff^(chi)(lambda_mu) ; eps_blk ),`

`rho_fail^(chi)(lambda_mu;eps_blk)`
`:= Q( zeta_fail^(chi)(lambda_mu) ; eps_blk ),`

and similarly for the conservative lower envelope,

`rho_suff^(J) := Q( zeta_suff^(J) ; eps_blk ),`

`rho_fail^(J) := Q( zeta_fail^(J) ; eps_blk ).`

The hard Family-1 constructive ceiling is

`rho_max^(F1)(eps_blk)`
`:= Q( zeta_max^(F1) ; eps_blk )`
` = [ 1 + (1 - 2 eps_blk) zeta_max^(F1) ] / [ 1 - eps_blk zeta_max^(F1) ].`

So the explicit support theorem is now entirely a statement about `rho_alpha`.

---

## 2. Exact unblocked Family-1 window

In the unblocked limit,

`eps_blk = 0,`

so

`Q(zeta;0) = 1 + zeta.`

Therefore the explicit Family-1 branch reduces to the exact unblocked loading-ratio thresholds

`rho_suff^(chi) = 1 + zeta_suff^(chi),`

`rho_fail^(chi) = 1 + zeta_fail^(chi),`

`rho_max^(F1) = 1 + zeta_max^(F1).`

Using the Stage-63/64 values at `lambda_mu = 1`,

`zeta_suff^(chi)(1) ≈ 2.46622291347846,`

`zeta_fail^(chi)(1) ≈ 2.46752913273870,`

`zeta_suff^(J)(1)   ≈ 2.44257571477179,`

`zeta_max^(F1)      ≈ 2.46752922945601,`

gives

`rho_suff^(chi)(1;0) ≈ 3.46622291347846,`

`rho_fail^(chi)(1;0) ≈ 3.46752913273870,`

`rho_suff^(J)(1;0)   ≈ 3.44257571477179,`

`rho_max^(F1)(0)     ≈ 3.46752922945601.`

So on the natural shell-weighted branch the guaranteed-success threshold lies only about

`1.30631597755e-3`

below the hard constructive ceiling, and the guaranteed-failure threshold differs from the ceiling only in the seventh decimal place.

---

## 3. Exact blocking condition in ratio form

The same denominator condition from Stage 64 remains necessary:

`eps_blk < 1 / zeta_max^(F1) ≈ 0.405263689711371.`

So the explicit Family-1 loading-ratio ceiling exists only while the blocked branch stays below that exact limit.

On any admissible branch,

`d rho_max^(F1) / d eps_blk`
`= zeta_max^(F1) [ zeta_max^(F1) - 1 ] / [ 1 - eps_blk zeta_max^(F1) ]^2 > 0,`

so blocking always raises the required loading ratio for a given support-ratio demand.

In other words, blocking hurts the support theorem exactly the way the product formulation suggested, but the statement is now written directly in the true reduced variable `rho_alpha`.

---

## 4. Final explicit Family-1 theorem in its cleanest form

The explicit Family-1 moving-throat support/source branch succeeds on the normalized outgoing quadrupole branch iff the selected loading ratio satisfies

`rho_alpha < rho_max^(F1)(eps_blk),`

and, more sharply, it is guaranteed to succeed once

`rho_alpha <= rho_suff^(chi)(lambda_mu;eps_blk)`.

At the natural shell-weighted normalization `lambda_mu = 1` and in the unblocked limit, that means

`rho_alpha <= 3.46622291347846`  -> guaranteed success,

`rho_alpha >= 3.46752913273870`  -> guaranteed failure,

with the absolute constructive ceiling

`rho_alpha < 3.46752922945601.`

So the explicit Family-1 theorem gap is no longer about the outgoing normalization factors separately and no longer about the product variables separately.
It is just this:

> **does the actual passive/outgoing moving-throat quadrupole branch produce a normalized loading ratio `rho_alpha = alpha_req/alpha_mix` below about `3.4675`?**

=== moving_throat_pde_stage070_outgoing_branch_loading_ratio_finish.md ===

# Moving-Throat PDE — Stage 70: Final Reduced Finish-Line for the Explicit Family-1 Outgoing Quadrupole Branch

## Purpose

Stages 65–69 compressed the full reduced moving-throat PDE to one surviving quadrupole residual, then to one explicit product window, and finally to one pure loading-ratio criterion.

This stage states the cleanest final verdict of that reduction chain.

The main conclusion is:

> **for the explicit Family-1 support/source branch, the remaining reduced theorem gap is exactly the normalized loading ratio**
>
> `rho_alpha = alpha_req / alpha_mix`
>
> **of the actual passive/outgoing moving-throat quadrupole branch.**

Everything else on the explicit support/source side is now fixed.

---

## 1. What has completely dropped out

Once the selected outgoing quadrupole branch is required to satisfy

`mhat_-^2 Gamma_{5,-} = 2 G / (5 c^5),`

Stages 68–69 show that the explicit Family-1 support theorem no longer depends separately on

- the selected conservative overlap `s_-`,
- the selected conservative stiffness `lambda_-`,
- the outgoing transfer factor `beta_0`,
- the source-map normalization `mhat_-`,
- the product variables `Pi_tr` and `C_mix`,
- or the intermediate Peclet demand `Pe_req`.

All of those collapse to the single ratio

`rho_alpha = alpha_req / alpha_mix`.

So the support/source side of the reduced moving-throat PDE is no longer a multivariable closure problem.
It is a one-number test.

---

## 2. Exact final Family-1 criterion

For the explicit Family-1 branch, the exact normalized success/failure theorem is

`rho_alpha <= rho_suff^(chi)(lambda_mu;eps_blk)`  -> guaranteed success,

`rho_alpha >= rho_fail^(chi)(lambda_mu;eps_blk)`  -> guaranteed failure,

with the hard constructive ceiling

`rho_alpha < rho_max^(F1)(eps_blk)`.

At the natural shell-weighted normalization `lambda_mu = 1` and in the unblocked limit,

`rho_alpha <= 3.46622291347846`  -> guaranteed success,

`rho_alpha >= 3.46752913273870`  -> guaranteed failure,

`rho_alpha < 3.46752922945601`   -> absolute constructive ceiling.

So the explicit support/source branch is already fully solved in reduced form.

---

## 3. What is still genuinely open

The only remaining reduced question is not on the support/source side anymore.
It is on the outgoing quadrupole side:

> **what value of `rho_alpha = alpha_req/alpha_mix` is actually selected by the passive/outgoing moving-throat quadrupole branch?**

That value still requires the actual branch solve.

So the program is now in the cleanest state it has reached so far:

- the explicit support/source side is finished,
- the explicit Family-1 ceiling is finished,
- the wall-depth side is finished,
- the outgoing normalization factors have been shown to cancel out of the explicit support theorem,
- and the surviving reduced theorem gap is exactly one normalized loading ratio on the physical quadrupole branch.

---

## 4. Best current expert verdict

At this point, the full reduced moving-throat PDE is no longer missing a large qualitative ingredient.
It is missing one sharp quantitative datum.

The strongest honest statement is therefore:

> **within the explicit Family-1 support/source branch, the reduced moving-throat PDE program is complete up to one final passive/outgoing quadrupole loading-ratio calculation.**

That is as close as the reduced program can get to a full explicit write-up without solving the last outgoing branch itself.

=== moving_throat_pde_stage071_loading_ratio_from_minimal_module.md ===

# Moving-Throat PDE — Stage 71: Loading-Ratio Extraction from the Minimal Isotropic Quadrupole Module

## Purpose

Stage 70 left the explicit Family-1 support/source theorem in its cleanest reduced form:

`rho_alpha = alpha_req / alpha_mix`

is the only quantity still needed from the actual passive/outgoing quadrupole branch.

The 2.5PN quadrupole audit had already isolated the smallest viable isotropic conservative precursor as

`Y_Q^cons(omega) = c0 + c1 / (1 - omega^2 / Omega_Q^2),`

with the unique positive-real coefficients

`c0 = 3/4,`
`c1 = 1/4,`
`Omega_Q = 3 c_s / (2 a).`

The next honest calculation is therefore to identify how that conservative precursor maps onto the support/source loading ratio.

This stage shows that under the natural **contact-plus-pole** reading of the explicit support/source branch,

- the mixed baseline contributes the static contact fraction,
- the extra support lane contributes the finite conservative pole,
- and the loading ratio is therefore fixed exactly by the precursor coefficients.

The main result is

`rho_alpha = 4/3,`
`zeta_req = 1/3,`
`Pi_tr = (4/3) C_mix.`

So the minimal isotropic quadrupole precursor lands in the exact **symmetric-lowest-twin** support regime and does not require non-twin asymmetry.

---

## 1. Natural contact-plus-pole identification

On the explicit support/source branch, the mixed baseline loading `alpha_mix` is the conservative directional loading already present before any extra support pole is added.
If the final passive/outgoing quadrupole branch is represented by one extra conservative support pole, then the natural normalized conservative precursor is

`Y_Q^cons(omega)`
`= alpha_mix / alpha_req`
`  + (alpha_req - alpha_mix)/alpha_req * 1/(1 - omega^2/Omega_Q^2).`

Equivalently, in the pure loading-ratio variable

`rho_alpha := alpha_req / alpha_mix,`

the same precursor is

`Y_Q^cons(omega)`
`= 1/rho_alpha`
`  + (rho_alpha - 1)/rho_alpha * 1/(1 - omega^2/Omega_Q^2).`

So the contact fraction and pole residue are

`c0 = 1/rho_alpha,`
`c1 = (rho_alpha - 1)/rho_alpha,`

with

`c0 + c1 = 1`

as required by the normalized static limit.

This gives the exact inverse formulas

`rho_alpha = 1/c0 = 1/(1-c1),`

`zeta_req = rho_alpha - 1 = c1/c0.`

So the support/source loading ratio is directly encoded in the static-contact / pole-residue split of the conservative quadrupole precursor.

---

## 2. Matching to the minimal isotropic quadrupole module

The 2.5PN quadrupole audit already fixed the smallest viable isotropic conservative precursor to

`c0 = 3/4,`
`c1 = 1/4.`

Inserting these into the exact inverse formulas above gives immediately

`rho_alpha = 1 / (3/4) = 4/3,`

`zeta_req = c1/c0 = (1/4)/(3/4) = 1/3.`

So the natural contact-plus-pole interpretation of the minimal isotropic quadrupole branch fixes the explicit support demand exactly:

`alpha_req / alpha_mix = 4/3,`
`(alpha_req - alpha_mix)/alpha_mix = 1/3.`

In product language, since Stage 68 proved

`Pi_tr / C_mix = alpha_req / alpha_mix,`

the same result is

`Pi_tr = (4/3) C_mix.`

---

## 3. Regime classification

Stage 35 split the support regimes into:

- `Pi_tr <= C_mix`              : mixed-only already enough,
- `C_mix < Pi_tr <= 2 C_mix`    : symmetric lowest twin enough,
- `Pi_tr > 2 C_mix`             : non-twin asymmetry required.

Because the minimal isotropic branch gives

`Pi_tr = (4/3) C_mix,`

it lands exactly in the middle regime:

`C_mix < Pi_tr < 2 C_mix.`

So the minimal isotropic passive/outgoing quadrupole branch:

- does require extra support beyond the mixed baseline,
- but requires only the symmetric lowest-twin lane,
- and does **not** require non-twin asymmetry.

Equivalently, in support-ratio language,

`zeta_req = 1/3 < 1,`

so the exact symmetric lowest twin already suffices.

---

## 4. What this changes

Before this step, the remaining reduced question was still “what value of `rho_alpha` does the actual passive/outgoing branch select?”

After this step, the minimal isotropic quadrupole precursor provides a very concrete answer:

> if the actual passive/outgoing branch is the natural contact-plus-pole realization of the minimal isotropic conservative module, then it selects
>
> `rho_alpha = 4/3`.

That is a much stronger statement than the earlier support-side ceiling alone.

It says the explicit Family-1 branch is not merely *capable* of surviving the outgoing quadrupole demand.
On the natural minimal isotropic branch, it only has to accommodate a support ratio of one third above the mixed baseline.

---

## 5. Best current theorem statement after Stage 71

Under the natural unblocked contact-plus-pole identification of the passive/outgoing isotropic quadrupole branch,

`Y_Q^cons(omega) = c0 + c1/(1 - omega^2/Omega_Q^2),`

with

`c0 = 3/4,`
`c1 = 1/4,`

the explicit Family-1 support theorem is driven by the exact loading ratio

`rho_alpha = 4/3,`

equivalently

`zeta_req = 1/3,`
`Pi_tr = (4/3) C_mix.`

So the reduced moving-throat PDE has advanced from a vague outgoing-branch loading question to a sharp statement:

> **if the actual passive/outgoing quadrupole branch realizes the minimal isotropic conservative precursor in the natural contact-plus-pole way, then the explicit support/source side is comfortably compatible with it.**

=== moving_throat_pde_stage072_family1_minimal_isotropic_verdict.md ===

# Moving-Throat PDE — Stage 72: Explicit Family-1 Verdict for the Minimal Isotropic Passive/Outgoing Quadrupole Branch

## Purpose

Stage 71 extracted the exact loading ratio selected by the natural contact-plus-pole realization of the minimal isotropic quadrupole precursor:

`rho_alpha = 4/3,`
`zeta_req = 1/3,`
`Pi_tr = (4/3) C_mix.`

The next honest step is to compare that branch directly against the explicit Family-1 support/source windows already frozen in Stages 62–70.

This stage shows that the match is not merely possible.
It is **strongly inside** the explicit Family-1 success region.

The main result is:

- the minimal isotropic passive/outgoing branch lies well below the explicit Family-1 loading-ratio ceiling,
- it lies in the exact symmetric-lowest-twin regime,
- and on the explicit Family-1 transport map it already succeeds at **zero** transport bias.

So under the natural contact-plus-pole identification, the explicit Family-1 support/source side is not the remaining reduced bottleneck anymore.

---

## 1. Exact comparison to the Family-1 loading-ratio window

From Stage 69, at the natural shell-weighted normalization `lambda_mu = 1` and in the unblocked limit,

`rho_suff^(chi) ≈ 3.46622291347846,`

`rho_fail^(chi) ≈ 3.46752913273870,`

`rho_max^(F1)   ≈ 3.46752922945601.`

From Stage 71, the minimal isotropic passive/outgoing quadrupole branch selects

`rho_alpha^(min) = 4/3 ≈ 1.33333333333333.`

So the exact margins are

`Delta_suff := rho_suff^(chi) - 4/3`
`           ≈ 2.13288958014513,`

`Delta_fail := rho_fail^(chi) - 4/3`
`           ≈ 2.13419579940537,`

`Delta_max  := rho_max^(F1) - 4/3`
`           ≈ 2.13419589612268.`

These are not small margins. The minimal isotropic branch sits far below even the guaranteed-success threshold.

---

## 2. Exact support-ratio and regime comparison

In the unblocked limit,

`zeta_req = rho_alpha - 1,`

so the same branch gives

`zeta_req^(min) = 1/3.`

The explicit Family-1 support ceiling from Stage 63/64 is

`zeta_max^(F1) ≈ 2.46752922945601.`

So the support-ratio margin is

`zeta_max^(F1) - zeta_req^(min)`
`≈ 2.13419589612268.`

And since

`0 < zeta_req^(min) = 1/3 < 1,`

the branch lies in the exact Stage-35 symmetric-lowest-twin window:

`C_mix < Pi_tr < 2 C_mix.`

So the minimal isotropic passive/outgoing branch does not need non-twin asymmetry at all.

---

## 3. Zero-transport-bias result on the explicit Family-1 branch

Stage 62 proved that the explicit Family-1 transport map obeys

`zeta_F1(Pe) = A_F1 Omega_Pe^2,`

with

`A_F1 ≈ 1.00005192880220,`

and that:

- if `zeta_req < A_F1`, the demand is already met at zero transport bias,
- if `A_F1 <= zeta_req <= zeta_max^(F1)`, a unique constructive `Pe_req` is needed.

For the minimal isotropic passive/outgoing branch,

`zeta_req^(min) = 1/3 < A_F1 ≈ 1.00005192880220.`

So the required transport bias is exactly

`Pe_req = 0`

on the explicit Family-1 branch.

Equivalently, the branch succeeds before any additional transport-driven overlap boost is turned on.

That is a very strong statement: the explicit Family-1 support/source branch already supports the minimal isotropic outgoing quadrupole demand in its zero-bias state.

---

## 4. What this changes

Before this step, the explicit Family-1 reduction had one unresolved outgoing-branch datum left: the loading ratio `rho_alpha`.

Stage 71 fixed that datum on the natural minimal isotropic passive/outgoing branch.
Stage 72 shows what that means numerically:

> **the explicit Family-1 support/source branch passes the minimal isotropic passive/outgoing quadrupole test with large margin.**

So for the explicit Family-1 branch, the remaining reduced theorem gap is no longer support sufficiency.
It is now the deeper question of whether the actual moving-throat grouped-`P2` / geometry branch really realizes the minimal isotropic contact-plus-pole conservative quadrupole module.

---

## 5. Best current theorem statement after Stage 72

Under the natural unblocked contact-plus-pole realization of the minimal isotropic passive/outgoing quadrupole branch,

`Y_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2),`

the explicit Family-1 support/source branch satisfies:

`rho_alpha = 4/3,`
`zeta_req = 1/3,`
`Pe_req = 0,`

and therefore lies safely inside the exact Family-1 success region.

So the explicit branch-level reduced theorem is now:

> **if the actual passive/outgoing moving-throat quadrupole branch realizes the minimal isotropic conservative precursor in the natural contact-plus-pole way, then the explicit Family-1 support/source side already succeeds without requiring transport bias or non-twin asymmetry.**

=== moving_throat_pde_stage073_updated_reduced_status.md ===

# Moving-Throat PDE — Stage 73: Updated Reduced Theorem Status After the Loading-Ratio Extraction

## Purpose

Stages 71–72 complete the explicit Family-1 support/source side under the natural minimal isotropic passive/outgoing quadrupole branch.

This note records the clean status update.

The key new fact is that the explicit Family-1 branch no longer merely has a large admissible ceiling.
Under the natural contact-plus-pole realization of the minimal isotropic conservative quadrupole precursor, it selects

`rho_alpha = 4/3,`
`zeta_req = 1/3,`
`Pe_req = 0,`

and therefore succeeds with large margin.

So the remaining reduced theorem gap is no longer on the explicit support/source side.

---

## 1. What is now settled inside the reduced Family-1 branch

The following are now fixed:

1. the explicit Family-1 support/source branch has a hard constructive ceiling;
2. the outgoing normalization factors cancel from the explicit support theorem;
3. the support theorem depends only on `rho_alpha = alpha_req/alpha_mix`;
4. the natural contact-plus-pole interpretation of the minimal isotropic conservative quadrupole precursor gives

   `rho_alpha = 4/3`;

5. that value lies strictly inside the exact Family-1 success region;
6. and on the explicit Family-1 transport map it already succeeds at zero transport bias.

So the explicit Family-1 support/source branch is no longer the place where the reduced program can fail first.

---

## 2. What is still genuinely open

The remaining reduced question is now narrower and more structural:

> does the actual moving-throat grouped-`P2` / geometry branch really realize the minimal isotropic conservative quadrupole precursor in the natural contact-plus-pole way?

Equivalently, the remaining reduced task is not to strengthen the explicit support/source side again.
It is to derive the conservative quadrupole module itself from the real grouped-`P2` / geometry branch of the moving throat.

That is exactly the point where the reduced program now reconnects to the already-frozen 3PN and 2.5PN theorem ledgers.

---

## 3. Best current expert verdict after Stage 73

For the explicit Family-1 branch, the reduced moving-throat PDE program has advanced one step further than Stage 70:

- the explicit support/source side is finished,
- the outgoing branch no longer leaves a free loading-ratio datum under the natural minimal isotropic contact-plus-pole identification,
- and the remaining reduced bottleneck is now the derivation of that minimal isotropic conservative quadrupole module from the actual grouped-`P2` / geometry throat branch.

That is the right place for the next derivation phase.

=== moving_throat_pde_stage074_grouped_p2_static_geometry_derivation.md ===

# Moving-Throat PDE — Stage 74: Deriving the `3/4 + 1/4` Conservative Quadrupole Module from the Grouped-`P2` + Geometry Split

## Purpose

Stage 71 extracted

`rho_alpha = 4/3`,

`zeta_req = 1/3`,

from the minimal isotropic conservative quadrupole module, but it did so by taking the contact-plus-pole representation as the natural reading of the explicit support/source branch.

The next honest step is sharper:

> derive that same `3/4 + 1/4` split directly from the conservative **grouped real `P2` + geometry** organization already frozen by the 3PN program.

This stage does that.

The main result is:

- if the isotropic grouped-`P2` conservative branch is carried by one effective pole,
- and if the 3PN geometry completion is genuinely static through `O(omega^4)`,

then the minimal isotropic 2.5PN branch identity forces

`K_pole = K0/4`,

`K_geom = 3 K0/4`,

and therefore

`Yhat_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)`.

So the `3/4 + 1/4` conservative module is no longer just a plausible parametrization. Under the minimal grouped-`P2` + static-geometry realization, it is forced.

---

## 1. Frozen input from the conservative hierarchy

The 3PN conservative split already says that the higher conservative payload is organized as

`free compiler image + grouped real P2 middle block + unique geometry completion`.

In the compact summary notation,

`Delta L3^GR = Delta l1 v^8 + L_(P2)^mid + Delta l15^(g) U^4`.

The important structural point is that the geometry lane appears there as a **static** completion, while the grouped-`P2` middle block carries the nontrivial higher-order conservative quadrupole structure.

The 2.5PN program, independently, already fixed the minimal isotropic quadrupole branch identity

`K0 K4 = 4 K2^2`,

equivalently in normalized language,

`u4 = 4 u2^2`.

So the only remaining move is to combine those two facts in the smallest consistent way.

---

## 2. Minimal grouped-`P2` + geometry realization

Take the isotropic conservative quadrupole module in the form

`K_Q^cons(omega) = K_geom + K_pole /(1 - omega^2/Omega_Q^2)`.

Here:

- `K_geom` is the static geometry completion,
- `K_pole` is the isotropic grouped-`P2` pole residue,
- `Omega_Q` is the effective isotropic grouped-`P2` pole.

This is the smallest realization compatible with the 3PN conservative split if the grouped-`P2` side is the only dynamic quadrupole lane and geometry contributes only the static completion.

Expanding at low frequency gives

`K_Q^cons(omega) = K0 + K2 omega^2 + K4 omega^4 + O(omega^6)`

with exact coefficients

`K0 = K_geom + K_pole`,

`K2 = K_pole / Omega_Q^2`,

`K4 = K_pole / Omega_Q^4`.

---

## 3. The branch identity forces the `3/4 + 1/4` split

Insert those coefficients into the minimal isotropic branch identity

`K0 K4 = 4 K2^2`.

This gives

`(K_geom + K_pole) * (K_pole / Omega_Q^4)`
`= 4 * (K_pole / Omega_Q^2)^2`.

Assuming the branch is nontrivial (`K_pole != 0`), the common factor `K_pole / Omega_Q^4` cancels and one finds

`K_geom + K_pole = 4 K_pole`.

So

`K_geom = 3 K_pole`.

Equivalently,

`K_pole = K0 / 4`,

`K_geom = 3 K0 / 4`.

Therefore the normalized conservative quadrupole response is forced to be

`Yhat_Q^cons(omega)`
`= K_Q^cons(omega) / K0`
`= 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)`.

This is exactly the minimal isotropic conservative quadrupole module previously isolated from the outgoing `l=2` moment matching.

So the earlier `3/4 + 1/4` structure is now recovered directly from the grouped-`P2` + static-geometry split.

---

## 4. Immediate corollary for the support/source loading ratio

Stage 71 already proved that on the explicit support/source branch the static contact fraction and the finite conservative pole map to the loading ratio as

`c0 = alpha_mix / alpha_req`,

`c1 = (alpha_req - alpha_mix)/alpha_req`,

with `c0 + c1 = 1`.

Using the newly derived grouped-`P2` + geometry split,

`c0 = 3/4`,

`c1 = 1/4`.

So

`alpha_mix / alpha_req = 3/4`,

which gives exactly

`rho_alpha = alpha_req / alpha_mix = 4/3`,

and

`zeta_req = (alpha_req - alpha_mix)/alpha_mix = 1/3`.

So the Stage-71 loading-ratio extraction is now a direct corollary of the grouped-`P2` + static-geometry realization.

---

## 5. What this actually proves

This stage does **not** prove that the full moving-throat PDE has already produced the minimal isotropic branch.

What it proves is narrower and more useful:

> if the actual conservative grouped-`P2` branch is one isotropic pole and the geometry lane contributes only the static completion through `O(omega^4)`, then the `3/4 + 1/4` conservative quadrupole module is forced algebraically.

So the remaining reduced theorem gap is now extremely sharp:

- either the real moving-throat branch obeys that minimal grouped-`P2` + static-geometry realization,
- or the missing PDE must generate extra dynamic geometry moments or a richer isotropic grouped-`P2` pole structure.

That is exactly the right question to carry into the next phase.

=== moving_throat_pde_stage075_dynamic_geometry_obstruction.md ===

# Moving-Throat PDE — Stage 75: Exact Obstruction Formula if the Geometry Lane Carries Dynamic Even Moments

## Purpose

Stage 74 showed that the grouped-`P2` + geometry split forces the `3/4 + 1/4` conservative quadrupole module **if**

- the grouped-`P2` side is one isotropic conservative pole,
- and the geometry lane is static through `O(omega^4)`.

The natural next question is:

> what exactly changes if the geometry lane is **not** purely static, but carries its own `omega^2` and `omega^4` moments?

This stage answers that exactly.

The main result is that the pole fraction is then no longer fixed to `1/4`.
Instead it becomes

`c_pole = (1 + eps_4) / [ 4 (1 + eps_2)^2 ]`,

where

`eps_2 = Omega_Q^2 K_(g,2) / K_pole`,

`eps_4 = Omega_Q^4 K_(g,4) / K_pole`.

So the `3/4 + 1/4` split is recovered **iff** the geometry lane is static at the relevant orders.

This is the cleanest reduced obstruction formula for the next phase.

---

## 1. Generalized isotropic grouped-`P2` + geometry ansatz

Allow the geometry lane to carry its own even moments:

`K_g(omega) = K_(g,0) + K_(g,2) omega^2 + K_(g,4) omega^4 + O(omega^6)`.

Keep the grouped-`P2` side as one isotropic conservative pole:

`K_P2(omega) = K_pole /(1 - omega^2/Omega_Q^2)`.

Then the total conservative isotropic quadrupole module is

`K_Q^cons(omega) = K_g(omega) + K_P2(omega)`.

Its low-frequency coefficients are

`K0 = K_(g,0) + K_pole`,

`K2 = K_(g,2) + K_pole/Omega_Q^2`,

`K4 = K_(g,4) + K_pole/Omega_Q^4`.

---

## 2. Exact branch identity with dynamic geometry

Imposing the same minimal isotropic branch identity

`K0 K4 = 4 K2^2`

now gives the exact relation

`(K_(g,0) + K_pole)( K_(g,4) + K_pole/Omega_Q^4 )`
`= 4 ( K_(g,2) + K_pole/Omega_Q^2 )^2.`

So the geometry-contact term is no longer forced to equal `3 K_pole` unless the dynamic geometry moments vanish.

Solving for `K_(g,0)` gives

`K_(g,0)`
`= 4 ( K_(g,2) + K_pole/Omega_Q^2 )^2 / ( K_(g,4) + K_pole/Omega_Q^4 )`
`  - K_pole.`

That is the exact obstruction formula.

---

## 3. Pole fraction in dimensionless contamination variables

Define the dimensionless geometry-contamination parameters

`eps_2 = Omega_Q^2 K_(g,2) / K_pole`,

`eps_4 = Omega_Q^4 K_(g,4) / K_pole`.

Then the total static normalization on the minimal branch is

`K0 = 4 K_pole (1 + eps_2)^2 / (1 + eps_4)`.

So the grouped-`P2` pole fraction becomes

`c_pole = K_pole / K0`
`       = (1 + eps_4) / [ 4 (1 + eps_2)^2 ]`.

And therefore the geometry contact fraction is

`c_geom = 1 - c_pole`.

In the strict static-geometry limit

`eps_2 = eps_4 = 0`,

this reduces exactly to

`c_pole = 1/4`,

`c_geom = 3/4`.

So the `3/4 + 1/4` split is not a generic identity of “grouped-`P2` plus geometry” in the abstract.
It is the exact consequence of the **static-geometry** realization.

---

## 4. Small-contamination expansion

For small geometry contamination,

`|eps_2| << 1`,

`|eps_4| << 1`,

the pole fraction expands as

`c_pole`
`= 1/4 [ 1 + eps_4 - 2 eps_2 + O(eps^2) ]`.

So:

- positive `eps_4` raises the pole fraction,
- positive `eps_2` lowers it twice as strongly at first order.

This is the cleanest reduced sensitivity formula for the next theorem gate.

---

## 5. What this changes

Stage 74 already showed that the minimal isotropic `3/4 + 1/4` module is forced if the geometry lane is static.

Stage 75 sharpens the remaining gap:

> the real moving-throat derivation now only needs to answer whether the geometry lane is static through `O(omega^4)` on the natural isotropic branch, or else compute the two contamination numbers `(eps_2, eps_4)`.

If both vanish, the Stage-74 result is exact.
If they do not, the contact/pole fractions are still fixed — but by the obstruction formula above rather than by the simple `3/4 + 1/4` split.

=== moving_throat_pde_stage076_grouped_p2_status_update.md ===

# Moving-Throat PDE — Stage 76: Updated Status After the Direct Grouped-`P2` + Geometry Derivation

## Purpose

Stages 74–75 close the exact next step that was left open after Stage 73.

They do two things:

1. derive the minimal isotropic `3/4 + 1/4` conservative quadrupole module directly from the grouped-`P2` + geometry split;
2. isolate the exact obstruction formula if the geometry lane carries dynamic `omega^2` or `omega^4` moments.

This note records the status update.

---

## 1. What is now settled

Inside the reduced hierarchy, the following implication is now exact:

- if the actual isotropic grouped-`P2` conservative branch is one pole,
- and if the geometry lane is static through `O(omega^4)`,

then the conservative quadrupole module is forced to be

`Yhat_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)`,

and therefore

`rho_alpha = 4/3`,

`zeta_req = 1/3`.

So the Stage-71/72 Family-1 verdict is now a direct corollary of the grouped-`P2` + static-geometry realization rather than an additional standalone assumption.

---

## 2. What remains genuinely open

The remaining reduced theorem gap is now extremely narrow.

It is no longer support-source sufficiency.
It is no longer the existence of the minimal isotropic module in the abstract.
It is now exactly this:

> does the real moving-throat branch realize
>
> - one isotropic grouped-`P2` pole,
> - and a geometry lane that is static through `O(omega^4)`?

If yes, the `3/4 + 1/4` split and the `rho_alpha = 4/3` verdict follow immediately.
If not, the exact deviation is still controlled by the Stage-75 obstruction formula

`c_pole = (1 + eps_4) / [ 4 (1 + eps_2)^2 ]`.

So the next real PDE-side task is no longer broad. It is to determine the two geometry-contamination numbers `eps_2` and `eps_4`, or prove that they vanish on the natural branch.

---

## 3. Best current expert verdict after Stage 76

At this point the reduced program is in its sharpest state yet.

- The explicit Family-1 support/source side is finished.
- The minimal isotropic quadrupole module is no longer just imported from the outgoing-moment analysis; it is directly recovered from the grouped-`P2` + static-geometry split.
- The only remaining reduced ambiguity is whether the real geometry lane is dynamically inert through `O(omega^4)` on the natural branch.

So the next derivation phase should be aimed squarely at the geometry lane, not at reopening the already-solved support/source side.

=== moving_throat_pde_stage077_isotropic_geometry_decoupling.md ===

# Moving-Throat PDE — Stage 77: Isotropic Geometry-Decoupling Theorem

## Purpose

Stage 75 showed that the only remaining reduced obstruction to the `3/4 + 1/4` conservative quadrupole module is dynamic contamination from the geometry lane through the two numbers

`eps_2 = Omega_Q^2 K_(g,2) / K_pole`,

`eps_4 = Omega_Q^4 K_(g,4) / K_pole`.

So the next honest check is not broad anymore:

> on the actual **isotropic** moving-throat branch, can the geometry lane contribute nonzero `omega^2` or `omega^4` moments to the grouped real `P2` conservative module at linear order?

This stage answers that directly.

The main result is:

- on an isotropic reference throat, the scalar/geometry `l=0` lane and the grouped real `l=2` lanes are exactly orthogonal in the quadratic wall theory,
- therefore the dynamic `l=0` geometry lane cannot contribute to the grouped real `P2` conservative module at linear order,
- and the only geometry contribution allowed on that branch is the already-identified **static** completion.

Equivalently,

`K_(g,2) = 0`,

`K_(g,4) = 0`,

so

`eps_2 = eps_4 = 0`

on the natural isotropic branch.

---

## 1. Quadratic wall action on the isotropic reference throat

From the distributed wall lift, the quadratic geometry action is

`S_eta^(2) = (1/2) int dt dw dOmega [`
`              mu_eta(w) (partial_t eta)^2`
`              - T_w(w) (partial_w eta)^2`
`              - T_Omega(w) eta (-Delta_(S^2)) eta`
`              - K_eta(w) eta^2 ]`.

All coefficients depend only on the axial coordinate `w`, not on the angles. So the angular operator is exactly `O(3)`-invariant.

Expand the wall field into scalar and quadrupole pieces,

`eta(Omega,w,t) = g(w,t) Y_00(Omega) + sum_A q_A(w,t) Y_(2A)(Omega) + ...`,

where `A in {20,21c,21s,22c,22s}`.

Because the operator is isotropic, every bilinear cross term between `Y_00` and any `Y_(2A)` is proportional to one of the angular integrals

`int Y_00 Y_(2A) dOmega`,

`int grad_(S^2) Y_00 . grad_(S^2) Y_(2A) dOmega`,

or equivalently

`int Y_00 (-Delta_(S^2)) Y_(2A) dOmega`.

All of them vanish exactly:

- `Y_00` is orthogonal to every `l=2` harmonic,
- `grad_(S^2) Y_00 = 0` because `Y_00` is constant,
- and `(-Delta_(S^2)) Y_(2A) = 6 Y_(2A)` reduces the third integral back to the first.

So the quadratic wall action is block diagonal in the `(l=0)` and `(l=2)` sectors.

---

## 2. Exact block structure of the isotropic reduced wall theory

After angular integration, the quadratic reduced action takes the schematic form

`S_red^(2) = (1/2) int dt dw [ g D_0 g + sum_A q_A D_2 q_A ]`,

with no bilinear mixing term of the form `g M_A q_A`.

So, on the isotropic branch,

`M_(0<->2) = 0`.

That means the scalar/geometry lane can renormalize only the scalar sector, while the grouped real `P2` channels evolve independently through their own isotropic quadrupole operator.

This is stronger than a small-coupling statement. It is an exact linear selection rule of the isotropic quadratic wall theory.

---

## 3. Consequence for the conservative grouped-`P2` quadrupole module

Now interpret the 3PN static geometry completion exactly the way the 3PN result says it should be interpreted:

- the grouped real `P2` lane carries the dynamic quadrupole pole structure,
- the leftover geometry completion is a scalar/pair-side static remainder.

Because the isotropic quadratic wall theory has no `l=0 <-> l=2` bilinear mixing, the scalar geometry lane cannot feed any dynamic even moments into the isotropic grouped-`P2` conservative quadrupole module.

So on the natural isotropic branch the effective conservative quadrupole module has the form

`K_Q^cons(omega) = K_(g,0) + K_pole /(1 - omega^2/Omega_Q^2)`

with

`K_(g,2) = 0`,

`K_(g,4) = 0`.

Therefore the contamination numbers of Stage 75 vanish:

`eps_2 = Omega_Q^2 K_(g,2) / K_pole = 0`,

`eps_4 = Omega_Q^4 K_(g,4) / K_pole = 0`.

And the exact obstruction formula collapses back to

`c_pole = 1/4`,

`c_geom = 3/4`.

So the `3/4 + 1/4` conservative quadrupole module is not obstructed on the natural isotropic branch.

---

## 4. What would be required to make `eps_2` or `eps_4` nonzero?

This theorem also makes the failure channels precise.

Nonzero geometry contamination at `O(omega^2)` or `O(omega^4)` requires at least one of the following:

1. **explicit angular anisotropy** in the quadratic wall operator, so that `l=0` and `l=2` cease to be orthogonal,
2. **a genuine second dynamic `l=2` geometry pole** independent of the grouped-`P2` pole already being used as the quadrupole carrier,
3. **nonlinear/higher-order backreaction** beyond the linear isotropic branch, which can induce contamination only beyond the exact linear theorem.

None of those are present in the natural minimal isotropic branch frozen by the present reduced hierarchy.

So the actual check requested after Stage 76 comes out cleanly:

> on the isotropic linear moving-throat branch, the geometry lane is dynamically inert through `O(omega^4)` with respect to the grouped real `P2` conservative quadrupole module.

---

## 5. Best current theorem statement after Stage 77

Inside the present reduced hierarchy,

- isotropy makes the quadratic wall operator block diagonal in `l`,
- the scalar/geometry lane is `l=0`,
- the conservative quadrupole carrier is the grouped real `l=2` bundle,
- so the geometry lane contributes only the already-known static completion and no dynamic `omega^2` or `omega^4` contamination.

Therefore

`eps_2 = eps_4 = 0`

on the natural isotropic branch,

and the Stage-74 `3/4 + 1/4` split is recovered exactly.

=== moving_throat_pde_stage078_second_order_geometry_contamination.md ===

# Moving-Throat PDE — Stage 78: First Nonzero Geometry Contamination Appears Only at Second Order in Anisotropy/Mixing

## Purpose

Stage 77 proved that on the natural isotropic branch the geometry lane is dynamically inert through `O(omega^4)` with respect to the grouped real `P2` conservative quadrupole module:

`eps_2 = eps_4 = 0`.

The next useful question is then:

> if isotropy is weakly broken and the scalar/geometry lane mixes with the grouped-`P2` quadrupole carrier, how do the contamination numbers turn on?

This stage answers that in the smallest exact reduced model.

The main result is:

- the first nonzero geometry contamination is **quadratic** in the mixing parameter,
- so the isotropic result is stable,
- and any deviation from the Stage-74 `3/4 + 1/4` split begins only at `O(chi^2)`.

---

## 1. Minimal mixed scalar-geometry / grouped-`P2` model

Take one grouped-`P2` quadrupole carrier `q` and one scalar/geometry mode `g`.

Write the conservative reduced kernel as

`D_q(omega) = K_stat + K_pole /(1 - omega^2/Omega_Q^2)`,

`D_g(omega) = G_0 + G_2 omega^2 + G_4 omega^4 + O(omega^6)`,

and let weak anisotropy or weak operator non-commutation generate a bilinear mixing term

`chi M_0 q g`.

Then the quadratic action is

`L = (1/2) q D_q q + (1/2) g D_g g + chi M_0 q g + J q`.

Integrating out the scalar/geometry mode gives the exact effective quadrupole kernel

`D_eff(omega) = D_q(omega) - chi^2 M_0^2 / D_g(omega)`.

So the whole contamination problem is encoded in the low-frequency expansion of the Schur-complement term.

---

## 2. Exact low-frequency contamination coefficients

Expand

`1 / D_g(omega)`

through `O(omega^4)`:

`1 / D_g(omega)`
`= 1/G_0`
`  - (G_2/G_0^2) omega^2`
`  + (G_2^2/G_0^3 - G_4/G_0^2) omega^4`
`  + O(omega^6)`.

So

`- chi^2 M_0^2 / D_g(omega)`
`= - chi^2 M_0^2 / G_0`
`  + chi^2 M_0^2 G_2 / G_0^2 * omega^2`
`  + chi^2 M_0^2 (G_0 G_4 - G_2^2) / G_0^3 * omega^4`
`  + O(omega^6)`.

The first term is only a static renormalization and can be absorbed into the geometry contact slot.
The dynamic contaminations are therefore

`K_(g,2)^eff = chi^2 M_0^2 G_2 / G_0^2`,

`K_(g,4)^eff = chi^2 M_0^2 (G_0 G_4 - G_2^2) / G_0^3`.

So the obstruction numbers of Stage 75 are

`eps_2 = Omega_Q^2 K_(g,2)^eff / K_pole`,

`eps_4 = Omega_Q^4 K_(g,4)^eff / K_pole`,

and both satisfy

`eps_2 = O(chi^2)`,

`eps_4 = O(chi^2)`.

That is the exact first stability theorem beyond Stage 77.

---

## 3. Consequence for the contact/pole fractions

Insert the small contamination into the Stage-75 obstruction formula

`c_pole = (1 + eps_4) / [ 4 (1 + eps_2)^2 ]`.

Expanding for small `chi` gives

`c_pole = 1/4 [ 1 + eps_4 - 2 eps_2 + O(chi^4) ]`.

Because both `eps_2` and `eps_4` are already `O(chi^2)`, the deviation from the `1/4` pole fraction is itself `O(chi^2)`.

So the Stage-74 `3/4 + 1/4` module is not just exact on the isotropic branch — it is also **perturbatively stable** against weak anisotropy/mixing.

---

## 4. Best current interpretation

Stages 77–78 together now give the exact reduced answer to the geometry-lane check:

1. on the natural isotropic branch,
   
   `eps_2 = eps_4 = 0`;

2. the first nonzero geometry contamination requires an explicit `l=0 <-> l=2` mixing source,
   and even then it appears only at
   
   `O(chi^2)`.

So the current reduced hierarchy supports the clean conservative reading:

- the grouped real `P2` branch is the dynamic quadrupole carrier,
- the geometry lane is static at the relevant order on the natural branch,
- and the Stage-74 `3/4 + 1/4` split is the correct actual branch value unless a genuine symmetry-breaking or extra `l=2` geometry pole is later found.

---

## 5. Best current theorem statement after Stage 78

Inside the present reduced hierarchy,

`eps_2 = eps_4 = 0`

exactly on the isotropic branch, and for weak symmetry breaking

`eps_2, eps_4 = O(chi^2)`.

So the dynamic geometry obstruction is absent on the actual isotropic branch and only enters at second order once an explicit mixing mechanism is turned on.

=== moving_throat_pde_stage079_geometry_lane_check_verdict.md ===

# Moving-Throat PDE — Stage 79: Actual Branch Verdict for the Geometry-Lane Check

## Purpose

Stages 77–78 were the direct answer to the open question left at Stage 76:

> on the actual moving-throat branch, does the geometry lane stay dynamically inert through `O(omega^4)`?

This note records the verdict.

---

## 1. Actual reduced-branch answer

Inside the present reduced hierarchy, the actual natural branch is the isotropic branch selected by:

- the isotropic quadratic wall operator,
- the grouped real `P2` quadrupole carrier,
- and the 3PN result that the leftover scalar contribution is a **static** geometry completion rather than a second dynamic quadrupole pole.

On that branch, Stage 77 proves exactly that the scalar/geometry `l=0` lane and the grouped real `l=2` bundle are block diagonal in the quadratic wall theory.
So the geometry lane cannot contribute dynamic `omega^2` or `omega^4` moments to the isotropic quadrupole module.

Therefore

`eps_2 = eps_4 = 0`

on the actual isotropic branch in the present hierarchy.

---

## 2. Consequence for the conservative quadrupole module

With the contamination numbers zero, the Stage-75 obstruction formula collapses to

`c_pole = 1/4`,

`c_geom = 3/4`.

So the actual branch realizes

`Yhat_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)`.

And therefore the already-carried support/source conclusions follow directly:

`rho_alpha = 4/3`,

`zeta_req = 1/3`.

So the geometry-lane check comes out **clean** on the actual isotropic branch.

---

## 3. What is still open after the check

This does **not** solve the full moving-throat PDE.
It removes one remaining reduced ambiguity.

What remains open is no longer whether the geometry lane contaminates the minimal isotropic quadrupole module. It does not, on the current actual branch.

The remaining serious open question is deeper and more physical:

> does the completed moving-throat PDE really realize the natural isotropic grouped-`P2` one-pole branch itself, with the passive/outgoing normalization required by the 2.5PN bridge?

So the geometry-lane check is now finished at reduced level.
The live gap is again the same narrow one isolated by the 2.5PN program: the final passive/outgoing quadrupole normalization on the true moving-throat branch.

=== moving_throat_pde_stage080_single_normalization_defect.md ===

# Moving-Throat PDE — Stage 80: The Actual Isotropic Passive/Outgoing Branch Collapses to a Single Normalization Defect

## Purpose

Stage 79 finished the last serious **conservative** ambiguity on the actual isotropic branch:

- the grouped real `P2` carrier is the only dynamic quadrupole lane,
- the geometry lane is static through `O(omega^4)`,
- and the conservative module is therefore

`Yhat_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)`.

The next honest step is to combine that actual conservative result with the already-frozen 2.5PN outgoing theorem target.

This stage shows that, once the actual branch is both

1. isotropic,
2. grouped-`P2` one-pole,
3. and passively/outgoingly completed in the natural `l=2` way,

then the entire reduced 2.5PN normalization problem collapses to **one scalar defect**.

---

## 1. Actual branch conservative module

Write the actual isotropic conservative quadrupole response in canonical invariant form as

`Kbar_Q^cons(omega)`
`= Kbar_0 [ 3/4 + (1/4)/(1 - omega^2/Omega_Q^2) ]`.

Then its low-frequency coefficients are exactly

`Kbar_2 = Kbar_0 / (4 Omega_Q^2),`

`Kbar_4 = Kbar_0 / (4 Omega_Q^4).`

So once `Kbar_0` and `Omega_Q` are known, the entire even conservative ledger is fixed.

---

## 2. Outgoing odd coefficient on the same branch

On the minimal isotropic outgoing `l=2` branch, the 2.5PN audit already fixed the odd coefficient algebraically as

`Gammabar_5 = 9 Kbar_2^(5/2) / Kbar_0^(3/2)`.

Substituting the actual one-pole conservative relations gives

`Gammabar_5 = 9 Kbar_0 / (32 Omega_Q^5).`

So the odd Burke–Thorne coefficient is not independent either.
It is already determined by the same two quantities `(Kbar_0, Omega_Q)`.

---

## 3. The GR target branch

The GR target on the same isotropic outgoing branch is

`Kbar_0^target = 64 G Omega_Q^5 / (45 c^5),`

equivalently, after using the already-carried geometric pole

`Omega_Q = 3 c_s / (2 a)`,

it becomes

`Kbar_0^target = 54 G c_s^5 / (5 a^5 c^5).`

Then automatically

`Kbar_2^target = Kbar_0^target / (4 Omega_Q^2),`

`Kbar_4^target = Kbar_0^target / (4 Omega_Q^4),`

`Gammabar_5^target = 2 G / (5 c^5).`

So the actual branch and the GR target have exactly the same algebraic structure.
The only possible mismatch is the overall normalization of `Kbar_0`.

---

## 4. Single normalization defect

Define the actual branch normalization defect by

`N_Q := Kbar_0 / Kbar_0^target`.

Then all three low-frequency targets scale by the same factor:

`Kbar_2 = N_Q Kbar_2^target,`

`Kbar_4 = N_Q Kbar_4^target,`

`Gammabar_5 = N_Q Gammabar_5^target = N_Q * 2 G / (5 c^5).`

Therefore the actual isotropic passive/outgoing branch satisfies the full reduced GR-like point-particle 2.5PN theorem **iff**

`N_Q = 1`.

Equivalently, the following defects are all the same number:

`R_0 := Kbar_0 / Kbar_0^target - 1,`

`R_2 := Kbar_2 / Kbar_2^target - 1,`

`R_4 := Kbar_4 / Kbar_4^target - 1,`

`R_5 := Gammabar_5 / (2G/(5c^5)) - 1,`

with

`R_0 = R_2 = R_4 = R_5 = N_Q - 1`.

So the final reduced theorem gap is now a **single scalar normalization defect**.

---

## 5. What this means physically

At this point, the reduced program has separated into two sharply different questions:

1. **support/source sufficiency** of the explicit branch,
2. **radiative normalization** of the actual outgoing quadrupole module.

The present stage shows that the second question is no longer a multi-parameter branch-selection problem.
Once the actual isotropic grouped-`P2` one-pole structure is accepted, it is only the one-number defect `N_Q`.

That is as narrow a reduced theorem gate as one can reasonably ask for before solving the full moving-throat PDE.

=== moving_throat_pde_stage081_family1_support_is_automatic.md ===

# Moving-Throat PDE — Stage 81: The Explicit Family-1 Support Test Is Automatic on the Actual Isotropic Branch

## Purpose

Stage 80 reduced the actual isotropic passive/outgoing quadrupole branch to one scalar normalization defect `N_Q`.

But the explicit branch still has a support/source theorem sitting beside that radiative theorem.
The next honest question is therefore:

> does the explicit Family-1 support/source side still need to be checked separately once the actual branch loading ratio is fixed?

This stage shows that the answer is effectively **no**.

On the actual isotropic branch,

`rho_alpha = 4/3`,

so even with finite blocking the exact support demand is so small that any explicit branch with support ceiling `zeta_max > 1` already passes it throughout the admissible blocked regime.

The Family-1 branch has

`zeta_max^(F1) ≈ 2.46752922945601 > 1`,

so its support/source side is automatic.

---

## 1. Exact blocked demand on the actual isotropic branch

Stage 68 gave the exact support-ratio demand

`zeta_req = (rho_alpha - 1) / (1 - eps_blk (2 - rho_alpha)).`

Inserting the actual isotropic branch value

`rho_alpha = 4/3`

gives

`zeta_req^(act)(eps_blk) = (1/3)/(1 - (2/3) eps_blk) = 1/(3 - 2 eps_blk).`

So the entire blocked support demand is one monotone increasing function of the blocking parameter.

---

## 2. Automatic support theorem for any branch with `zeta_max > 1`

Assume a branch has support ceiling `zeta_max > 1`, and use the same admissible blocking window already frozen earlier,

`0 <= eps_blk < 1/zeta_max`.

Because `zeta_req^(act)` increases with `eps_blk`, its worst value on the admissible window is bounded by

`zeta_req^(act) < 1 / (3 - 2/zeta_max).`

Now compare that to `zeta_max` itself. The inequality

`1 / (3 - 2/zeta_max) < zeta_max`

is equivalent to

`3 zeta_max (zeta_max - 1)/(3 zeta_max - 2) > 0,`

which is true for every `zeta_max > 1`.

Therefore:

> for the actual isotropic branch with `rho_alpha = 4/3`, any explicit support/source family whose constructive ceiling satisfies `zeta_max > 1` already passes the support test throughout the full admissible blocked regime.

So the support/source side is automatic once the actual isotropic branch is fixed.

---

## 3. Family-1 specialization

For the explicit Family-1 branch,

`zeta_max^(F1) ≈ 2.46752922945601 > 1`.

Therefore the automatic theorem above applies immediately.

At zero blocking,

`zeta_req^(act)(0) = 1/3,`

and even at the edge of the admissible Family-1 blocked window,

`eps_blk -> 1 / zeta_max^(F1)`,

the demand remains bounded by

`zeta_req^(act) < zeta_max^(F1)`.

Numerically,

`zeta_req^(act) < 0.456730991107963 < 2.46752922945601`.

So the explicit Family-1 support/source branch no longer carries an independent reduced theorem burden.
It is already safe once the actual isotropic outgoing branch is accepted.

---

## 4. What remains after this step

After Stages 80–81, the reduced theorem split is now fully sharp:

- **support/source theorem:** automatic on the actual isotropic branch,
- **radiative theorem:** controlled by the single normalization defect `N_Q`.

So the explicit Family-1 branch has now dropped out of the active uncertainty ledger.
The only remaining reduced theorem question is whether the completed moving-throat PDE gives

`N_Q = 1`.

=== moving_throat_pde_stage082_reduced_finish_line.md ===

# Moving-Throat PDE — Stage 82: The Reduced Finish Line After the Geometry-Lane Check

## Purpose

This note records the sharpest reduced-program verdict reached so far.

After the grouped-`P2` conservative closure, the geometry-lane check, the explicit Family-1 support/source analysis, and the passive/outgoing quadrupole reduction, the remaining theorem gap is no longer diffuse.

It is one scalar normalization datum.

---

## 1. What is now fixed inside the reduced hierarchy

On the actual natural isotropic branch we now have:

1. `eps_2 = eps_4 = 0`, so the geometry lane is dynamically inert through `O(omega^4)` with respect to the grouped real quadrupole carrier.
2. The conservative quadrupole module is exactly

   `Yhat_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)`.

3. Therefore

   `rho_alpha = 4/3`,

   `zeta_req = 1/3` in the unblocked limit.

4. On the explicit Family-1 branch the support/source side is automatic even in the admissible blocked regime.

So the reduced program is no longer carrying separate uncertainties in

- geometry contamination,
- contact/pole splitting,
- support/source sufficiency,
- or explicit Family-1 viability.

---

## 2. Single remaining reduced theorem gate

Write the actual isotropic passive/outgoing branch in canonical invariant form as

`Kbar_Q^cons(omega) = Kbar_0 [ 3/4 + (1/4)/(1 - omega^2/Omega_Q^2) ]`.

Then define

`N_Q := Kbar_0 / Kbar_0^target,`

with

`Kbar_0^target = 64 G Omega_Q^5 / (45 c^5)`

or equivalently

`Kbar_0^target = 54 G c_s^5 / (5 a^5 c^5)`

after using `Omega_Q = 3 c_s/(2a)`.

The full reduced GR-like point-particle 2.5PN closure on the actual isotropic branch is now equivalent to

`N_Q = 1`.

Everything else has already been reduced away.

---

## 3. Practical meaning

The completed moving-throat PDE no longer needs to answer a large family of loosely related questions in order to close the reduced theorem.

It only needs to determine the actual radiative normalization of the passive/outgoing grouped-`P2` quadrupole branch.

In other words, the open PDE task is now exactly one of the following equivalent statements:

- compute `Kbar_0` on the actual passive/outgoing branch,
- compute `Gammabar_5` on that branch,
- or compute the scalar defect `N_Q - 1`.

The three questions are equivalent on the actual isotropic branch.

---

## 4. Best current theorem statement

Inside the present reduced hierarchy, the moving-throat PDE program has reached the following finish line:

> the actual isotropic grouped-`P2` one-pole branch is conservatively clean, the explicit Family-1 support/source side is already sufficient, and the only remaining reduced theorem gap is the single passive/outgoing quadrupole normalization defect `N_Q - 1`.

That is the narrowest and strongest honest carry-forward statement available before solving the full moving-throat PDE normalization problem.

=== moving_throat_pde_stage083_outgoing_normalization_factorization.md ===

# Moving-Throat PDE — Stage 83: Exact Factorization of the Last 2.5PN Defect into Conservative and Outgoing Pieces

## Purpose

Stage 82 reduced the full isotropic point-particle 2.5PN problem to one scalar defect

`N_Q := Kbar_0 / Kbar_0^target`.

That statement was made on the *minimal passive/outgoing grouped-`P2` branch*, where the outgoing `l=2` odd coefficient already matches the canonical compact value. The next honest step is to relax that one assumption slightly and ask:

> if the actual passive/outgoing moving-throat branch is still a one-pole grouped-`P2` branch, but its leading outgoing odd coefficient is renormalized, how does the final defect factorize?

This stage shows that the answer is exact.

Introduce one dimensionless outgoing-normalization factor `chi_Q` by writing the actual normalized retarded branch as

`Yhat_Q^ret(omega)`
`= 3/4 + (1/4)/(1 - omega^2/Omega_Q^2 - i chi_Q sigma_Q^can omega^5) + O(omega^6),`

where the canonical compact outgoing coefficient is

`sigma_Q^can = 9/(8 Omega_Q^5) = 4 a^5/(27 c_s^5)`.

Then the invariant low-frequency tuple is

`Kbar_2 = Kbar_0/(4 Omega_Q^2),`

`Kbar_4 = Kbar_0/(4 Omega_Q^4),`

`Gammabar_5 = chi_Q * 9 Kbar_0/(32 Omega_Q^5).`

So the even conservative defect and the odd retarded defect separate cleanly:

- the even moments depend only on `N_Q`,
- the odd Burke–Thorne coefficient depends on `chi_Q N_Q`.

That is the exact factorization of the last 2.5PN obstruction.

---

## 1. Conservative target ratios remain scalar

Define the canonical invariant GR targets on the isotropic branch by

`Kbar_0^target = 64 G Omega_Q^5/(45 c^5),`

`Kbar_2^target = Kbar_0^target/(4 Omega_Q^2),`

`Kbar_4^target = Kbar_0^target/(4 Omega_Q^4).`

Then with

`N_Q := Kbar_0/Kbar_0^target`,

the exact even ratios are still

`Kbar_2/Kbar_2^target = N_Q,`

`Kbar_4/Kbar_4^target = N_Q.`

So the conservative/even side remains one-number clean.

---

## 2. Odd branch ratio picks up exactly one extra factor

The 2.5PN target odd coefficient is

`Gammabar_5^target = 2 G/(5 c^5)`.

With the renormalized outgoing branch,

`Gammabar_5 = chi_Q * 9 Kbar_0/(32 Omega_Q^5)`

and therefore

`Gammabar_5 / Gammabar_5^target = chi_Q N_Q.`

So all deviations from the exact compact outgoing `l=2` fingerprint are captured by the single multiplier `chi_Q`.

---

## 3. Observable point-particle normalization condition

The actual point-particle observable branch includes the source-map factor `mhat_0`, so the full odd normalization condition is

`mhat_0^2 Gammabar_5 = 2 G/(5 c^5)`.

Substituting the renormalized one-pole branch gives

`mhat_0^2 chi_Q N_Q = 1.`

Equivalently,

`N_Q = 1/(mhat_0^2 chi_Q).`

This is the exact factorized form of the last reduced 2.5PN defect.

---

## 4. Meaning

Stage 82 said the remaining reduced theorem gap was one scalar normalization defect `N_Q - 1`.

This stage refines that statement:

- if the branch is only *conservatively* known, the even defect is `N_Q - 1`;
- if the branch is *retarded* but not yet proven canonical, the full 2.5PN odd defect is `mhat_0^2 chi_Q N_Q - 1`;
- and all genuinely new retarded uncertainty sits in one number only:

`chi_Q`.

So the problem is no longer “some unknown outgoing structure.” It is exactly the leading outgoing-normalization factor of the actual grouped-`P2` one-pole branch.

=== moving_throat_pde_stage084_natural_source_map_reduction.md ===

# Moving-Throat PDE — Stage 84: On the Natural Source-Map Branch the Last Reduced 2.5PN Obstruction is Purely Outgoing

## Purpose

Stage 83 showed that the full odd normalization condition factorizes exactly as

`mhat_0^2 chi_Q N_Q = 1`.

The next honest step is to combine that with the already-carried source-map result of the natural orbital/worldtube STF branch:

`mhat_0 = 1 + O(a^2/r^2)`.

This stage shows that in the strict point-particle limit the remaining reduced 2.5PN obstruction is no longer a mixed source/outgoing product. It is purely the outgoing-normalization factor `chi_Q`.

---

## 1. Point-particle natural branch

On the natural source-map branch isolated by the 2.5PN audit, one has

`mhat_0 = 1 + O(a^2/r^2)`.

So in the strict point-particle limit,

`mhat_0 -> 1`.

Then the exact Stage-83 factorization becomes

`N_Q = 1/chi_Q`.

So:

- if `chi_Q = 1`, then `N_Q = 1`,
- if `chi_Q != 1`, then the entire reduced 2.5PN mismatch is just the inverse outgoing renormalization.

---

## 2. Canonical compact outgoing branch

Stage 80 already fixed the unique minimal passive/outgoing grouped-`P2` one-pole completion that matches the exact compact outgoing `l=2` fingerprint:

`Yhat_Q^ret(omega)`
`= 3/4 + (1/4)/(1 - omega^2/Omega_Q^2 - i sigma_Q^can omega^5) + O(omega^6),`

with

`sigma_Q^can = 9/(8 Omega_Q^5) = 4 a^5/(27 c_s^5)`.

By definition, that canonical compact outgoing branch has

`chi_Q = 1`.

Therefore on the natural source-map branch,

`N_Q = 1`.

So the reduced 2.5PN theorem closes automatically **if** the actual passive/outgoing moving-throat branch is exactly the canonical compact one-pole `l=2` completion.

---

## 3. Small deviation parameter

Define the outgoing-normalization defect by

`Delta_Q := chi_Q - 1`.

Then on the natural source-map branch,

`N_Q = 1/(1 + Delta_Q)`

and therefore

`N_Q - 1 = -Delta_Q/(1 + Delta_Q)`.

For a small outgoing deviation,

`N_Q - 1 = -Delta_Q + O(Delta_Q^2).`

So the last reduced theorem gap is linearly controlled by the first outgoing-normalization defect.

---

## 4. Meaning

At this point the reduced hierarchy has done everything it can without solving the actual passive/outgoing DtN problem.

What remains is no longer a broad question about conservative structure, support sufficiency, or source-map ambiguity. It is just this:

> does the actual passive/outgoing moving-throat quadrupole branch realize `chi_Q = 1`, or does it carry a nontrivial outgoing-normalization defect `Delta_Q`?

That is the cleanest reduced 2.5PN finish line reached so far.

=== moving_throat_pde_stage085_higher_odd_irrelevance.md ===

# Moving-Throat PDE — Stage 85: At 2.5PN the Only Live Retarded Obstruction is the Leading `omega^5` Outgoing Normalization

## Purpose

Stage 84 reduced the last finite reduced 2.5PN question to one outgoing-normalization factor `chi_Q`.

A natural objection remains:

> what if the true moving-throat PDE contains extra retarded structure beyond the canonical one-pole `l=2` denominator?

This stage shows that, at the level of the 2.5PN theorem, any extra retarded structure that first enters at `O(omega^7)` or higher is irrelevant.

So the only live retarded obstruction at 2.5PN is the leading `omega^5` outgoing-normalization factor `chi_Q`.

---

## 1. Generalized one-pole denominator with a higher odd tail

Take the most general one-pole grouped-`P2` denominator consistent with the already-fixed conservative branch but allowing one extra higher odd term:

`Yhat_Q^ret(omega)`
`= 3/4 + (1/4)/(1 - omega^2/Omega_Q^2 - i chi_Q sigma_Q^can omega^5 - i tau_Q omega^7) + O(omega^8).`

Expanding through `O(omega^5)` gives

`Yhat_Q^ret(omega)`
`= 1 + omega^2/(4 Omega_Q^2) + omega^4/(4 Omega_Q^4)`
`  + i chi_Q 9/(32 Omega_Q^5) omega^5 + O(omega^6).`

The extra `tau_Q` term does not appear until `O(omega^7)`.

So at 2.5PN order, the only retarded parameter that matters is `chi_Q`.

---

## 2. Consequence for the theorem ledger

Combining Stages 83–84 with the `O(omega^7)` irrelevance statement gives the exact reduced 2.5PN hierarchy:

1. `eps_2 = eps_4 = 0` kills the conservative geometry contamination,
2. `rho_alpha = 4/3` and the Family-1 analysis make support/source automatic,
3. `mhat_0 = 1 + O(a^2/r^2)` removes the source-map ambiguity on the natural branch,
4. every higher odd retarded coefficient starting at `omega^7` is invisible to the 2.5PN theorem,
5. so the only remaining reduced obstruction is

`Delta_Q = chi_Q - 1`.

If `Delta_Q = 0`, the reduced point-particle 2.5PN theorem is closed.

---

## 3. Best current finish-line statement

Inside the present reduced hierarchy, the moving-throat PDE program is now reduced to one sharply isolated PDE-facing question:

> does the actual passive/outgoing grouped-`P2` quadrupole branch have the canonical compact outgoing `omega^5` coefficient, i.e. `chi_Q = 1`?

Everything else that could have obstructed the reduced 2.5PN theorem has either been fixed or pushed above the relevant order.

=== moving_throat_pde_stage086_reduced_25pn_conditional_closure.md ===

# Moving-Throat PDE — Stage 86: Conditional Reduced 2.5PN Closure

## Statement

Inside the present reduced hierarchy, the following are now fixed:

- the actual isotropic grouped-`P2` conservative branch is geometry-clean through `O(omega^4)`,
- the conservative quadrupole split is exactly `3/4 + 1/4`,
- the selected support/source ratio is `rho_alpha = 4/3`,
- the explicit Family-1 support/source branch is automatically sufficient,
- the natural source-map branch gives `mhat_0 = 1 + O(a^2/r^2)`,
- and all higher odd retarded data beginning at `O(omega^7)` are irrelevant to the point-particle 2.5PN theorem.

Therefore the reduced 2.5PN theorem is conditionally closed by one and only one remaining branch datum:

`chi_Q`.

More precisely:

- if the actual passive/outgoing grouped-`P2` quadrupole branch satisfies `chi_Q = 1`,
  then the reduced GR-like point-particle 2.5PN theorem is closed;
- if `chi_Q != 1`, then the entire remaining reduced failure is measured by

`Delta_Q := chi_Q - 1`.

So the remaining PDE-facing problem is no longer “derive 2.5PN somehow.”
It is:

> compute the leading outgoing `omega^5` normalization of the actual passive/outgoing grouped-`P2` quadrupole branch and determine whether it equals the canonical compact outgoing value.

=== moving_throat_pde_stage087_outgoing_dtn_fingerprint.md ===

# Moving-Throat PDE — Stage 87: Exact Outgoing `l=2` DtN Fingerprint

## Goal

Compute the canonical compact passive/outgoing quadrupole normalization directly from an explicit outgoing `l=2` Dirichlet-to-Neumann model instead of carrying the coefficient `chi_Q` symbolically.

## Exact outgoing DtN model

Let
\[
z := \frac{a\omega}{c_s},
\]
and let the outgoing partial wave be the spherical Hankel mode
\[
h_2^{(1)}(z)=j_2(z)+i\,y_2(z).
\]

The exact `l=2` outgoing DtN operator is
\[
\Lambda_2^{\rm out}(z)
=
z\,\frac{d}{dz}\ln h_2^{(1)}(z)
=
z\,\frac{h_2^{(1)\prime}(z)}{h_2^{(1)}(z)}.
\]

Its small-\(z\) expansion is
\[
\boxed{
\Lambda_2^{\rm out}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}
+i\,\frac{z^5}{9}
-\frac{2z^6}{27}
-i\,\frac{z^7}{27}
+O(z^8).
}
\]

So the exact static slot is \(\Lambda_2^{\rm out}(0)=-3\).

## Normalized outgoing quadrupole admittance

Define the normalized outgoing branch by
\[
\widehat Y_2^{\rm out}(z)
:=
\frac{\Lambda_2^{\rm out}(0)}{\Lambda_2^{\rm out}(z)}
=
-\frac{3}{\Lambda_2^{\rm out}(z)}.
\]

Then
\[
\boxed{
\widehat Y_2^{\rm out}(z)
=
1+\frac{z^2}{9}
+\frac{4z^4}{81}
+i\,\frac{z^5}{27}
-\frac{11z^6}{729}
-i\,\frac{z^7}{243}
+O(z^8).
}
\]

Restoring \(\omega\),
\[
z=\frac{a\omega}{c_s},
\]
this becomes
\[
\boxed{
\widehat Y_2^{\rm out}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\,\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6).
}
\]

This is exactly the outgoing compact `l=2` fingerprint used earlier in the reduced 2.5PN program, but now derived directly from the explicit DtN model.

## Consequence

The canonical outgoing coefficient is therefore not free. On the exact spherical outgoing `l=2` DtN branch, the leading odd quadrupole coefficient is fixed to
\[
\Gamma_{5,\rm can}^{\rm DtN}=\frac{a^5}{27c_s^5}.
\]

So any later normalization mismatch must come from branch selection or source normalization, not from ambiguity in the canonical outgoing `l=2` DtN model itself.

=== moving_throat_pde_stage088_chiQ_fix_from_outgoing_dtn.md ===

# Moving-Throat PDE — Stage 88: Exact Fixing of `chi_Q`

## Goal

Use the explicit DtN fingerprint to determine the last reduced 2.5PN scalar
\[
\chi_Q
\]
in the canonical retarded grouped-`P2` one-pole-plus-contact module.

## Retarded minimal isotropic grouped-`P2` module

Write the retarded normalized quadrupole module as
\[
\widehat Y_Q^{\rm ret}(\omega)
=
\frac34
+
\frac14\,
\frac{1}{
1-\omega^2/\Omega_Q^2
-i\,\chi_Q\,\sigma_Q^{\rm can}\,\omega^5
}
+O(\omega^6).
\]

The conservative moment match already fixed
\[
\Omega_Q=\frac{3c_s}{2a}.
\]

Define
\[
\sigma_Q^{\rm can}:=\frac{9}{8\Omega_Q^5}.
\]
Using \(\Omega_Q=3c_s/(2a)\),
\[
\boxed{
\sigma_Q^{\rm can}
=
\frac{4a^5}{27c_s^5}.
}
\]

## Low-frequency expansion

Expanding through \(O(\omega^5)\) gives
\[
\widehat Y_Q^{\rm ret}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\,\chi_Q\,\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6).
\]

From Stage 87, the explicit outgoing DtN branch gives
\[
\widehat Y_2^{\rm out}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\,\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6).
\]

Matching the \(O(\omega^5)\) coefficient yields
\[
\boxed{\chi_Q=1.}
\]

So on the canonical compact passive/outgoing grouped-`P2` DtN branch, the last reduced normalization scalar is fixed exactly.

## General deformed DtN obstruction

A useful way to parametrize the only remaining PDE-facing freedom is to deform the outgoing DtN operator by
\[
\Lambda_2^{\rm def}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}
+i\,\xi_Q\,\frac{z^5}{9}+O(z^6).
\]

Then
\[
\widehat Y_2^{\rm def}(z)
=
1+\frac{z^2}{9}
+\frac{4z^4}{81}
+i\,\xi_Q\,\frac{z^5}{27}
+O(z^6),
\]
so
\[
\boxed{\chi_Q=\xi_Q.}
\]

That means the only reduced 2.5PN obstruction left after the present calculation is a deviation of the actual moving-throat DtN branch from the canonical outgoing `l=2` coefficient \(\xi_Q=1\).

=== moving_throat_pde_stage089_canonical_outgoing_reduced_closure.md ===

# Moving-Throat PDE — Stage 89: Reduced 2.5PN Closure on the Canonical Outgoing DtN Branch

## Goal

Insert the explicit DtN result
\[
\chi_Q=1
\]
into the reduced 2.5PN normalization stack and state the strongest honest closure now available.

## Reduced normalization stack

The current reduced program had already shown that
\[
\hat m_0^{\,2}\,\chi_Q\,N_Q=1,
\]
where:

- \(\hat m_0\) is the canonical source-map overlap,
- \(\chi_Q\) is the outgoing quadrupole normalization factor,
- \(N_Q\) is the conservative quadrupole normalization defect.

On the natural point-particle source-map branch,
\[
\hat m_0=1+O(a^2/r^2).
\]

From Stage 88, the canonical compact passive/outgoing grouped-`P2` DtN model gives
\[
\chi_Q=1.
\]

Therefore, in the strict point-particle limit,
\[
\boxed{N_Q=1.}
\]

## Canonical invariant coefficients

With \(N_Q=1\), the canonical invariant low-frequency quadrupole coefficients are fixed to their target values:
\[
\boxed{
\overline K_0=\frac{54Gc_s^5}{5a^5c^5},
\qquad
\overline K_2=\frac{6Gc_s^3}{5a^3c^5},
\qquad
\overline K_4=\frac{8Gc_s}{15ac^5},
\qquad
\overline \Gamma_5=\frac{2G}{5c^5}.
}
\]

Equivalently, the normalized odd coefficient is exactly the GR/Burke–Thorne target:
\[
\boxed{
\gamma_{\rm quad}^{\rm eff}
=
\hat m_0^{\,2}\Gamma_5
=
\frac{2G}{5c^5}
}
\]
on the canonical branch.

## Strongest honest theorem statement

Inside the present reduced hierarchy:

1. the isotropic grouped-`P2` conservative branch is geometry-clean through \(O(\omega^4)\),
2. the conservative split is exactly `3/4 + 1/4`,
3. the explicit Family-1 support/source branch is automatically sufficient,
4. the natural source map gives \(\hat m_0=1+O(a^2/r^2)\),
5. and the explicit compact passive/outgoing grouped-`P2` DtN model gives \(\chi_Q=1\).

Therefore the reduced nonspinning point-particle 2.5PN theorem is **closed on the canonical outgoing DtN branch**.

## What is still genuinely open

What remains open is no longer a reduced PN bookkeeping problem. It is the deeper PDE-side branch-selection theorem:

> Does the completed moving-throat PDE actually realize the canonical compact passive/outgoing grouped-`P2` DtN branch, rather than a deformed branch with \(\xi_Q\neq1\)?

So the reduced theorem is now finished **conditional only on actual branch realization**, not on any remaining reduced-sector normalization ambiguity.

=== moving_throat_pde_stage090_general_dtn_deformation.md ===

# Moving-Throat PDE — Stage 90: General Isotropic `l=2` DtN Deformation Algebra

## Goal

Replace the symbolic branch parameter `xi_Q` by an explicit low-frequency moving-throat DtN deformation model and derive the exact map from deformation data to the retarded quadrupole normalization factor `chi_Q`.

## Canonical outgoing branch

From the exact outgoing spherical `l=2` DtN model,
\[
\Lambda_2^{\rm out}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}+i\frac{z^5}{9}+O(z^6),
\qquad z:=\frac{a\omega}{c_s}.
\]

The corresponding normalized branch is
\[
\widehat Y_2^{\rm out}(z)
=1+\frac{z^2}{9}+\frac{4z^4}{81}+i\frac{z^5}{27}+O(z^6),
\]
so the canonical outgoing normalization is `chi_Q = 1`.

## General isotropic deformation model

Take the first explicit isotropic moving-throat deformation family
\[
\boxed{
\Lambda_2^{\rm def}(z)
=
S\,\Lambda_2^{\rm out}(\beta z)
+\Sigma_0+\Sigma_2 z^2+\Sigma_4 z^4+i\Sigma_5 z^5+O(z^6).
}
\]
Here:

- `S` is an overall mouth normalization,
- `beta` rescales the effective outgoing argument,
- `Sigma_0, Sigma_2, Sigma_4` are isotropic throat-core even self-energy data,
- `Sigma_5` is the first extra isotropic odd `l=2` core outlet.

Expanding,
\[
\Lambda_2^{\rm def}(z)=L_0+L_2 z^2+L_4 z^4+iL_5 z^5+O(z^6),
\]
with
\[
L_0=-3S+\Sigma_0,
\qquad
L_2=\frac{S\beta^2}{3}+\Sigma_2,
\qquad
L_4=\frac{S\beta^4}{9}+\Sigma_4,
\qquad
L_5=\frac{S\beta^5}{9}+\Sigma_5.
\]

## Normalized deformation law

Normalize by the actual static slot,
\[
\widehat Y_2^{\rm def}(z):=\frac{L_0}{L_0+L_2 z^2+L_4 z^4+iL_5 z^5}+O(z^6).
\]
Then the exact low-frequency coefficients are
\[
\boxed{
\widehat Y_2^{\rm def}(z)
=
1-\frac{L_2}{L_0}z^2
+\left(\frac{L_2^2}{L_0^2}-\frac{L_4}{L_0}\right)z^4
-i\frac{L_5}{L_0}z^5
+O(z^6).
}
\]

So the deformation-normalized quadrupole factor is
\[
\chi_Q=\frac{-L_5/L_0}{1/27}.
\]

## Exact canonical-even matching conditions

Demand that the deformed branch preserve the canonical conservative even fingerprint,
\[
\frac{z^2}{9},\qquad \frac{4z^4}{81}.
\]
Then
\[
-\frac{L_2}{L_0}=\frac19,
\qquad
\frac{L_2^2}{L_0^2}-\frac{L_4}{L_0}=\frac{4}{81}.
\]
Solving for the even throat-core coefficients gives
\[
\boxed{
\Sigma_2=-\frac{3S\beta^2-3S+\Sigma_0}{9},
\qquad
\Sigma_4=-\frac{3S\beta^4-3S+\Sigma_0}{27}.
}
\]

With those imposed, the exact odd normalization becomes
\[
\boxed{
\chi_Q=
\frac{3\big(S\beta^5+9\Sigma_5\big)}{3S-\Sigma_0}.
}
\]
Equivalently,
\[
\boxed{
\chi_Q-1=
\frac{3S(\beta^5-1)+\Sigma_0+27\Sigma_5}{3S-\Sigma_0}.
}
\]

## Consequence

This is the first explicit moving-throat DtN deformation model for the last reduced 2.5PN scalar. It shows exactly which isotropic branch data can move the canonical value `chi_Q = 1`:

- argument deformation `beta`,
- static additive throat-core shift `Sigma_0`,
- odd `l=2` throat-core outlet `Sigma_5`.

Overall scale `S` is not itself an independent obstruction; it only enters through the ratios above.

=== moving_throat_pde_stage091_robustness_classes.md ===

# Moving-Throat PDE — Stage 91: Robustness Classes for `chi_Q`

## Goal

Classify which explicit isotropic DtN perturbations leave the canonical outgoing normalization untouched and which ones genuinely shift it.

## Class A — pure scale deformation

Take
\[
\Lambda_2^{\rm def}(z)=S\Lambda_2^{\rm out}(z).
\]
Then
\[
\widehat Y_2^{\rm def}(z)=\widehat Y_2^{\rm out}(z),
\]
exactly, so
\[
\boxed{\chi_Q=1.}
\]
Thus pure mouth normalization is invisible to the normalized outgoing quadrupole fingerprint.

## Class B — pure scale+argument deformation

Take
\[
\Lambda_2^{\rm def}(z)=S\Lambda_2^{\rm out}(\beta z).
\]
Then
\[
\widehat Y_2^{\rm def}(z)
=
1+\frac{\beta^2 z^2}{9}+\frac{4\beta^4 z^4}{81}+i\frac{\beta^5 z^5}{27}+O(z^6).
\]
So if the conservative even fingerprint is kept canonical,
\[
\beta^2=1,
\qquad
\beta^4=1,
\]
which on the natural positive branch forces
\[
\boxed{\beta=1.}
\]
Hence
\[
\boxed{\chi_Q=1.}
\]
So a pure effective radius/sound-speed rescaling cannot move the canonical outgoing normalization without simultaneously spoiling the already fixed even moments.

## Class C — additive isotropic throat-core channel

Set `beta = 1` but allow a genuine throat-core self-energy,
\[
\Lambda_2^{\rm def}(z)
=
S\Lambda_2^{\rm out}(z)
+\Sigma_0+\Sigma_2 z^2+\Sigma_4 z^4+i\Sigma_5 z^5.
\]
If the even moments are held canonical, then
\[
\Sigma_2=-\frac{\Sigma_0}{9},
\qquad
\Sigma_4=-\frac{\Sigma_0}{27},
\]
and the odd normalization is
\[
\boxed{
\chi_Q=\frac{3(S+9\Sigma_5)}{3S-\Sigma_0}.
}
\]
So an additive throat-core channel can move `chi_Q` even while leaving the lower even moments unchanged.

A particularly clean special case is a purely even additive core with `Sigma_5 = 0`:
\[
\boxed{
\chi_Q=\frac{3S}{3S-\Sigma_0}.
}
\]
So a static additive core shift is a genuine candidate branch-selection effect.

## Exact preservation submanifold

The deformation preserves the canonical outgoing normalization iff
\[
\chi_Q=1.
\]
The exact condition is
\[
\boxed{
\Sigma_5=\frac{S(1-\beta^5)}{9}-\frac{\Sigma_0}{27}.
}
\]
This includes as special cases:

- pure scale deformation,
- pure scale+argument deformation with `beta = 1`,
- additive core deformations whose odd slot is locked to the static shift.

## Consequence

The canonical value `chi_Q = 1` is robust against:

1. overall mouth normalization,
2. pure effective radius/sound-speed rescaling once the conservative even fingerprint is fixed.

It is shifted only by a genuine isotropic throat-core self-energy that is not on the exact preservation submanifold above.

=== moving_throat_pde_stage092_linearized_branch_selection.md ===

# Moving-Throat PDE — Stage 92: Linearized Branch-Selection Law Near the Canonical Outgoing Branch

## Goal

Extract the first-order deformation law around the canonical outgoing branch and isolate the minimal PDE-facing branch-selection data.

## Linearized deformation ansatz

Expand around the canonical branch by writing
\[
S=1+\varepsilon s,
\qquad
\beta=1+\varepsilon b,
\qquad
\Sigma_0=\varepsilon a_0,
\qquad
\Sigma_5=\varepsilon a_5,
\]
with the even slots adjusted to preserve the canonical conservative fingerprint.

Then the exact Stage-90 formula
\[
\chi_Q=\frac{3(S\beta^5+9\Sigma_5)}{3S-\Sigma_0}
\]
expands to
\[
\boxed{
\chi_Q
=
1+
\varepsilon\left(5b+\frac{a_0}{3}+9a_5\right)
+O(\varepsilon^2).
}
\]

## Immediate implications

### 1. Pure overall scaling cancels
The coefficient `s` drops out to first order. This matches the exact Class-A invariance: overall mouth renormalization does not move the normalized outgoing quadrupole fingerprint.

### 2. Effective argument deformation matters only through `b`
A small outgoing-argument deformation contributes as
\[
\delta\chi_Q^{(\beta)}=5b.
\]
So a genuine effective radius/sound-speed shift that is not removed by the even-moment matching will move the outgoing normalization linearly.

### 3. Static additive slot and odd core slot are the remaining direct branch data
The first-order additive contributions are
\[
\delta\chi_Q^{(0)}=\frac{a_0}{3},
\qquad
\delta\chi_Q^{(5)}=9a_5.
\]
So the minimal isotropic branch-selection triple is
\[
\boxed{(b,\,a_0,\,a_5).}
\]

## Linearized preservation condition

To keep the canonical outgoing normalization at first order, the deformation must satisfy
\[
\boxed{
5b+\frac{a_0}{3}+9a_5=0.
}
\]
Equivalently,
\[
\boxed{
a_5=-\frac{5b}{9}-\frac{a_0}{27}.
}
\]

## Consequence for the reduced theorem program

After the grouped-`P2` conservative split, geometry cleaning, and Family-1 support sufficiency results, the final reduced 2.5PN branch-selection problem can now be stated very sharply:

> compute the isotropic moving-throat DtN branch data `(b, a_0, a_5)` and test whether they satisfy the exact nonlinear condition of Stage 91, or at least the linearized condition above.

So the remaining PDE-facing ambiguity is no longer an open-ended “deformed branch somehow.” It is a small explicit set of outgoing-branch deformation scalars.

=== moving_throat_pde_stage093_robin_outlet_model.md ===

# Moving-Throat PDE — Stage 93: Explicit Isotropic Robin Outlet Model

## Goal

Replace the abstract static deformation slot `a_0` by the first explicit isotropic geometric outlet model and compute its exact effect on the outgoing quadrupole normalization.

## Raw Robin DtN deformation

Take the canonical outgoing `l=2` branch
\[
\Lambda_2^{\rm out}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}+i\frac{z^5}{9}+O(z^6),
\qquad z:=\frac{a\omega}{c_s},
\]
and add a dimensionless isotropic Robin throat-core shift
\[
\boxed{\Lambda_2^{\rm R}(z)=\Lambda_2^{\rm out}(z)+\rho_R,}
\qquad \rho_R:=a h_R.
\]
This is the natural reduced form of an isotropic Robin-type mouth law.

Normalizing by the actual static slot,
\[
\widehat Y_2^{\rm R}(z)
:=
\frac{-3+\rho_R}{\Lambda_2^{\rm R}(z)},
\]
one finds the exact low-frequency expansion
\[
\boxed{
\widehat Y_2^{\rm R}(z)
=
1+
\frac{z^2}{9-3\rho_R}
+
\frac{4-\rho_R}{9(3-\rho_R)^2}z^4
+
 i\frac{z^5}{27-9\rho_R}
+O(z^6).
}
\]

So the raw isotropic Robin outlet changes both the even branch and the odd outgoing normalization.

## Effective quadrupole normalization factor

Relative to the canonical outgoing `l=2` fingerprint,
\[
\widehat Y_2^{\rm out}(z)=1+\frac{z^2}{9}+\frac{4z^4}{81}+i\frac{z^5}{27}+O(z^6),
\]
the raw Robin outlet carries the exact normalization factor
\[
\boxed{
\chi_Q^{\rm R}=rac{3}{3-\rho_R}.
}
\]
Expanding around the canonical branch,
\[
\chi_Q^{\rm R}=1+\frac{\rho_R}{3}+\frac{\rho_R^2}{9}+O(\rho_R^3).
\]
So a pure isotropic Robin core generically pushes the branch away from `chi_Q=1`.

## Branch-selection triple

In the Stage-92 linearized notation, the raw Robin outlet is the pure static deformation
\[
\boxed{(b,a_0,a_5)=(0,\rho_R,0).}
\]
Therefore the linearized outgoing-normalization shift is
\[
\delta\chi_Q^{\rm R}=\frac{\rho_R}{3}+O(\rho_R^2).
\]

## Consequence

A pure isotropic geometric Robin outlet is **not** automatically harmless. By itself it deforms both the canonical even fingerprint and the odd normalization. So if the already-fixed conservative grouped-`P_2` branch is to survive, a Robin core must either be negligible or be compensated by additional outlet structure.

=== moving_throat_pde_stage094_mixed_sidechannel_pole.md ===

# Moving-Throat PDE — Stage 94: Explicit Mixed `A_w/F_{\mu w}`-Type Side-Channel Pole

## Goal

Build the first explicit isotropic hidden mixed-sector side-channel model and check whether it can preserve the already-fixed conservative even `l=2` fingerprint.

## Hidden pole model

Take the canonical outgoing branch and add the first isotropic Schur-complement-style mixed side-channel
\[
\boxed{
\Lambda_2^{\rm mix}(z)
=
\Lambda_2^{\rm out}(z)
-
\frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5}
+O(z^6),
}
\]
with
- `sigma_W > 0` the static mixed loading,
- `kappa_W > 0` the hidden even dispersion scale,
- `gamma_W > 0` the hidden odd outgoing coefficient.

Expanding,
\[
\Lambda_2^{\rm mix}(z)
=
-(3+\sigma_W)
+
\left(\frac13-\sigma_W\kappa_W\right)z^2
+
\left(\frac19-\sigma_W\kappa_W^2\right)z^4
+
 i\left(\frac19-\sigma_W\gamma_W\right)z^5
+O(z^6).
\]

## Exact even-branch no-go

Demand that this branch preserve the canonical conservative even fingerprint. The `z^2` condition gives
\[
-\frac{L_2}{L_0}=\frac19
\quad\Longrightarrow\quad
\boxed{\kappa_W=-\frac19.}
\]
This is already incompatible with a standard positive passive pole parameter `\kappa_W>0`.

Even if one formally inserts that value, the `z^4` condition becomes
\[
\frac{L_2^2}{L_0^2}-\frac{L_4}{L_0}=\frac{4}{81}
\quad\Longrightarrow\quad
\boxed{\sigma_W=0.}
\]
So a standalone isotropic hidden pole of this type cannot sit on the already-fixed canonical even branch unless it is absent.

## Outgoing-normalization shift

The raw branch normalization factor is
\[
\boxed{
\chi_Q^{\rm mix}
=
\frac{3(1-9\sigma_W\gamma_W)}{3+\sigma_W}.
}
\]
For small loading,
\[
\chi_Q^{\rm mix}
=
1-
\sigma_W\left(\frac13+9\gamma_W\right)
+O(\sigma_W^2).
\]
So the linearized branch-selection triple is
\[
\boxed{(b,a_0,a_5)=(0,-\sigma_W,-\sigma_W\gamma_W).}
\]

## Consequence

A naive passive mixed `A_w/F_{\mu w}` side-channel pole is **too rigid**. It generically shifts the outgoing normalization and, more importantly, it cannot preserve the already-fixed canonical even `l=2` branch. If a mixed sector survives on the actual branch, it must appear in a more structured, compensated outlet law.

=== moving_throat_pde_stage095_hybrid_robin_mixed_compensation.md ===

# Moving-Throat PDE — Stage 95: Exact Robin–Mixed Compensation Law

## Goal

Combine the geometric Robin core and the hidden mixed side-channel and determine whether an explicit compensated moving-throat outlet can preserve the canonical outgoing quadrupole branch.

## Hybrid outlet model

Take
\[
\boxed{
\Lambda_2^{\rm hyb}(z)
=
\Lambda_2^{\rm out}(z)
+ho_R
-
\frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5}
+O(z^6).
}
\]
Expanding,
\[
L_0=-3+\rho_R-\sigma_W,
\qquad
L_2=\frac13-\sigma_W\kappa_W,
\qquad
L_4=\frac19-\sigma_W\kappa_W^2,
\qquad
L_5=\frac19-\sigma_W\gamma_W.
\]

## Exact canonical-even solutions

Imposing the canonical conservative even fingerprint,
\[
-\frac{L_2}{L_0}=\frac19,
\qquad
\frac{L_2^2}{L_0^2}-\frac{L_4}{L_0}=\frac{4}{81},
\]
yields exactly two branches:
\[
\boxed{\rho_R=\sigma_W,\qquad \kappa_W=0,}
\]
or
\[
\boxed{\rho_R=4\sigma_W,\qquad \kappa_W=\frac13.}
\]

The first is a trivial cancellation branch: the static Robin shift simply cancels the static mixed loading.
The second is the nontrivial compensated branch.

## Nontrivial compensated branch

On
\[
\rho_R=4\sigma_W,
\qquad
\kappa_W=\frac13,
\]
one finds
\[
\Lambda_2^{\rm hyb}(z)
=
-3(1-\sigma_W)
+
\frac{1-\sigma_W}{3}z^2
+
\frac{1-\sigma_W}{9}z^4
+
 i\left(\frac19-\sigma_W\gamma_W\right)z^5
+O(z^6).
\]
So the exact outgoing-normalization factor is
\[
\boxed{
\chi_Q^{\rm hyb}
=
\frac{1-9\sigma_W\gamma_W}{1-\sigma_W}.
}
\]
Canonical outgoing normalization is preserved iff
\[
\boxed{\gamma_W=\frac19.}
\]
With that value,
\[
\boxed{
\Lambda_2^{\rm hyb}(z)=(1-\sigma_W)\Lambda_2^{\rm out}(z)+O(z^6),
}
\]
so the whole hybrid outlet collapses to the harmless pure-scale deformation class.

## Branch-selection data

On the nontrivial compensated branch, the net Stage-92 deformation data are
\[
\boxed{(b,a_0,a_5)=\left(0,\,3\sigma_W,\,-\sigma_W\gamma_W\right).}
\]
Hence the linearized preservation condition becomes
\[
\frac{a_0}{3}+9a_5=\sigma_W(1-9\gamma_W)=0,
\]
again giving
\[
\boxed{\gamma_W=\frac19.}
\]

## Consequence

This is the first explicit compensated moving-throat outlet model that preserves the canonical outgoing quadrupole branch. It shows that neither a pure geometric Robin core nor a naive hidden mixed pole is sufficient by itself, but a specific Robin–mixed balance law can reduce the whole deformation to a pure mouth renormalization, which is exactly the robust class already identified in Stages 90–92.

=== moving_throat_pde_stage096_outlet_model_status.md ===

# Moving-Throat PDE — Stage 96: Outlet-Model Status Update

## What is now explicit

The abstract isotropic branch-selection triple `(b,a_0,a_5)` is no longer just symbolic. Three explicit outlet classes have now been checked:

1. **Pure geometric Robin core**
   \[
   \Lambda_2^{\rm R}=\Lambda_2^{\rm out}+\rho_R,
   \]
   which shifts the canonical outgoing normalization as
   \[
   \chi_Q^{\rm R}=\frac{3}{3-\rho_R}.
   \]

2. **Standalone mixed `A_w/F_{\mu w}`-type hidden pole**
   \[
   \Lambda_2^{\rm mix}=\Lambda_2^{\rm out}-\frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5},
   \]
   which cannot preserve the already-fixed even `l=2` branch unless it vanishes.

3. **Hybrid Robin–mixed outlet**
   \[
   \Lambda_2^{\rm hyb}=\Lambda_2^{\rm out}+\rho_R-
   \frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5},
   \]
   which admits one nontrivial compensated canonical-even branch
   \[
   \rho_R=4\sigma_W,
   \qquad
   \kappa_W=\frac13.
   \]
   On that branch,
   \[
   \chi_Q^{\rm hyb}=\frac{1-9\sigma_W\gamma_W}{1-\sigma_W},
   \]
   and canonical outgoing normalization is preserved exactly when
   \[
   \gamma_W=\frac19.
   \]

## PDE-facing consequence

The remaining isotropic branch-selection question is now much sharper than before.
The completed moving-throat PDE does **not** need to decide among an unlimited space of deformations. At low frequency, the first explicit outlet audit says the actual branch must land in one of two categories:

- either a nearly pure scale/argument deformation class, which is harmless,
- or a compensated Robin–mixed outlet of the exact type above.

Pure Robin loading alone and a naive standalone hidden mixed pole are not enough.

=== moving_throat_pde_stage097_concrete_core_schur.md ===


# Moving-Throat PDE — Stage 97: Concrete Two-Channel Core Outlet Model

## Goal

Replace the reduced outlet coefficients `(\rho_R,\sigma_W,\kappa_W,\gamma_W)` by a concrete throat-core response model with explicit internal variables.

## Core variables

Let `u(\omega)` be the isotropic `l=2` mouth amplitude seen by the exterior branch. Introduce two internal core variables:

- `s(\omega)`: a static geometric compliance coordinate,
- `q(\omega)`: a mixed `A_w/F_{\mu w}`-type side-channel coordinate.

Take the linear core system
\[
\begin{pmatrix}
K_s & \lambda\\[2pt]
\lambda & -K_q D_W^{\rm bare}(z)
\end{pmatrix}
\binom{s}{q}
=
u\binom{g_s}{g_q},
\qquad
z:=\frac{a\omega}{c_s},
\]
with bare mixed denominator
\[
D_W^{\rm bare}(z)=1-\kappa_0 z^2-i\gamma_0 z^5+O(z^6).
\]
The mouth feedback is defined by
\[
\delta\Lambda_{\rm core}(z)\,u = g_s s + g_q q.
\]

This is the simplest concrete isotropic core model that contains:
- a static geometric compliance,
- a dynamic mixed side-channel,
- and a nontrivial static/mixed hybridization `\lambda`.

## Exact Schur-complement outlet

Eliminating `(s,q)` gives the exact core correction
\[
\boxed{
\delta\Lambda_{\rm core}(z)
=
\frac{g_s^2}{K_s}
-
\frac{(K_s g_q-\lambda g_s)^2}
{K_s\big(K_sK_q D_W^{\rm bare}(z)+\lambda^2\big)}.
}
\]

Define the dimensionless hybridization ratio
\[
r_c:=\frac{\lambda^2}{K_sK_q}.
\]
Then the outlet takes the reduced Stage-95 form
\[
\boxed{
\delta\Lambda_{\rm core}(z)
=
\rho_c
-
\frac{\sigma_c}{1-\kappa_c z^2-i\gamma_c z^5}
+O(z^6),
}
\]
with exact core-level identifications
\[
\boxed{
\rho_c=\frac{g_s^2}{K_s},
}
\]
\[
\boxed{
\sigma_c=
\frac{(K_s g_q-\lambda g_s)^2}
{K_s^2K_q(1+r_c)},
}
\]
\[
\boxed{
\kappa_c=\frac{\kappa_0}{1+r_c},
\qquad
\gamma_c=\frac{\gamma_0}{1+r_c}.
}
\]

So the concrete core model reproduces the full reduced Robin–mixed hybrid outlet with no extra assumptions.

## Interpretation

This is already a significant narrowing. The outlet is no longer described by four unrelated reduced coefficients. It is controlled by:

- one static shell stiffness `K_s`,
- one mixed-channel stiffness `K_q`,
- one static/mixed hybridization `\lambda`,
- and two mouth couplings `(g_s,g_q)`,
- plus the bare mixed low-frequency pair `(\kappa_0,\gamma_0)`.

The next question is whether this concrete core model can *naturally* land on the compensated canonical branch found algebraically in Stage 95.

=== moving_throat_pde_stage098_core_balance_compensation.md ===


# Moving-Throat PDE — Stage 98: Exact Core-Balance Compensation Theorem

## Goal

Determine when the concrete two-channel core model of Stage 97 lands exactly on the nontrivial compensated canonical branch found algebraically in Stage 95.

## Canonical branch conditions

The Stage-95 nontrivial compensated outlet is
\[
\delta\Lambda_{\rm can}(z)
=
4\sigma_*-\frac{\sigma_*}{1-z^2/3-i z^5/9},
\]
which is exactly the branch that preserves the canonical outgoing quadrupole fingerprint.

For the concrete core model, this requires
\[
\boxed{
\rho_c=4\sigma_c,
\qquad
\kappa_c=\frac13,
\qquad
\gamma_c=\frac19.
}
\]

## Exact coupling-balance law

Using the Stage-97 identifications,
\[
\rho_c=\frac{g_s^2}{K_s},
\qquad
\sigma_c=
\frac{(K_s g_q-\lambda g_s)^2}
{K_s^2K_q(1+r_c)},
\qquad
r_c=\frac{\lambda^2}{K_sK_q},
\]
the nontrivial compensation condition becomes
\[
\boxed{
g_s^2\bigl(K_sK_q+\lambda^2\bigr)
=
4\bigl(K_s g_q-\lambda g_s\bigr)^2.
}
\]
Solving for the mixed coupling gives the exact two-branch law
\[
\boxed{
g_q=
\frac{g_s}{2K_s}
\left(
2\lambda \pm \sqrt{K_sK_q+\lambda^2}
\right).
}
\]

So the required Robin–mixed balance is not a mystery coefficient fit. It is an explicit codimension-one surface in the concrete core-coupling space.

## Bare mixed-channel normalization conditions

The even/odd preservation conditions become
\[
\boxed{
\kappa_0=\frac{1+r_c}{3},
\qquad
\gamma_0=\frac{1+r_c}{9}.
}
\]
So the bare mixed side-channel must itself be a scale-deformed copy of the canonical compact outgoing branch.

## Exact collapse to the Stage-95 branch

On the coupling-balance surface together with the bare-channel conditions above, the concrete core outlet reduces identically to
\[
\boxed{
\delta\Lambda_{\rm core}(z)
=
4\sigma_*-\frac{\sigma_*}{1-z^2/3-i z^5/9},
\qquad
\sigma_*=\frac{g_s^2}{4K_s}.
}
\]
Adding this to the canonical exterior branch
\[
\Lambda_2^{\rm out}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}+i\frac{z^5}{9}+O(z^6)
\]
gives a normalized response with exactly the same outgoing fingerprint:
\[
\widehat Y_2(z)=1+\frac{z^2}{9}+\frac{4z^4}{81}+i\frac{z^5}{27}+O(z^6).
\]

So the concrete core model does not merely mimic the reduced algebra. It reproduces the exact nontrivial compensated branch.

## Interpretation

This is the first explicit throat-core theorem in the outlet program:

- the Stage-95 compensation law is not just an algebraic accident,
- it is realized by a concrete two-channel core model,
- and the surviving free data are sharply reduced to one coupling-balance surface plus one scale-deformed bare mixed outlet.

=== moving_throat_pde_stage099_dn_mixed_tube_realization.md ===


# Moving-Throat PDE — Stage 99: Finite D/N Mixed-Tube Realization

## Goal

Give the concrete core-balance theorem a geometric throat-core realization rather than leaving the bare mixed-channel data `(\kappa_0,\gamma_0)` abstract.

## Bare mixed D/N tube

Take the mixed side-channel to live on a finite auxiliary tube of length `L_W` with the first D/N half-wave:
\[
q''+k^2 q=0,
\qquad
q(0)=0,
\qquad
q'(L_W)=0.
\]
Then
\[
\boxed{
k_W=\frac{\pi}{2L_W},
\qquad
\Omega_W=\frac{\pi c_s}{2L_W}.
}
\]
Writing the side-channel pole denominator in the outlet variable
\[
z=\frac{a\omega}{c_s},
\]
the bare even coefficient is
\[
\boxed{
\kappa_0=\frac{\omega^2}{\Omega_W^2}\Big/\! z^2
=\frac{4L_W^2}{\pi^2 a^2}.
}
\]

## Compensation-selected tube length

Stage 98 requires
\[
\kappa_0=\frac{1+r_c}{3},
\qquad
r_c=\frac{\lambda^2}{K_sK_q}.
\]
So the auxiliary mixed-tube length is fixed to
\[
\boxed{
L_W
=
\frac{\pi a}{2}\sqrt{\frac{1+r_c}{3}}.
}
\]
Thus the static/mixed hybridization `r_c` has a direct geometric meaning: it sets the ratio between the auxiliary mixed-tube length and the exterior mouth radius.

## Bare outgoing normalization

The same compensation theorem requires
\[
\gamma_0=\frac{1+r_c}{9}.
\]
A simple concrete realization is to take the bare mixed outlet to be a pure-scale deformation of the canonical compact outgoing `l=2` branch:
\[
\boxed{
D_W^{\rm bare}(z)
=
(1+r_c)\left(1-\frac{z^2}{3}-i\frac{z^5}{9}\right)+O(z^6).
}
\]
Then the hybridization factor `(1+r_c)` is removed exactly by the Stage-97 denominator renormalization, leaving the canonical final coefficients
\[
\kappa_c=\frac13,
\qquad
\gamma_c=\frac19.
\]

## Consequence

The outlet program is now geometrically concrete.

A fully compensated canonical branch exists whenever:

1. the core couplings lie on the exact balance surface
   \[
   g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2,
   \]
2. the auxiliary mixed side-channel is the first D/N half-wave of length
   \[
   L_W=\frac{\pi a}{2}\sqrt{\frac{1+r_c}{3}},
   \]
3. and its bare outgoing port is a pure-scale deformation of the canonical compact outgoing `l=2` branch.

So the remaining PDE-side question is no longer “some unknown outlet deformation.”
It is whether the actual moving-throat core realizes this specific D/N-tube + coupling-balance structure.

=== moving_throat_pde_stage100_outlet_core_status.md ===


# Moving-Throat PDE — Stage 100: Concrete Outlet-Core Status

The outlet-deformation problem is now much narrower.

## What is now explicit

The isotropic outlet is no longer described by four reduced coefficients
\[
(\rho_R,\sigma_W,\kappa_W,\gamma_W).
\]
In the first concrete core model it is generated by:

- a static shell compliance `K_s`,
- a mixed-channel stiffness `K_q`,
- a static/mixed hybridization `\lambda`,
- mouth couplings `(g_s,g_q)`,
- and an auxiliary mixed D/N tube of length `L_W`.

The exact eliminated outlet is
\[
\delta\Lambda_{\rm core}(z)
=
\rho_c-\frac{\sigma_c}{1-\kappa_c z^2-i\gamma_c z^5},
\]
with explicit formulas derived in Stages 97–99.

## What preserves the canonical outgoing branch

The concrete outlet lands exactly on the nontrivial compensated branch iff

\[
g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2,
\]
\[
L_W=\frac{\pi a}{2}\sqrt{\frac{1+r_c}{3}},
\qquad
r_c=\frac{\lambda^2}{K_sK_q},
\]
and the bare mixed outlet is a pure-scale deformation of the canonical compact outgoing `l=2` branch.

On that surface the full outlet reduces to
\[
\Lambda_{\rm eff}(z)
=
\Lambda_2^{\rm out}(z)
+
4\sigma_*-\frac{\sigma_*}{1-z^2/3-i z^5/9},
\qquad
\sigma_*=\frac{g_s^2}{4K_s},
\]
which preserves
\[
\widehat Y_2(z)=1+\frac{z^2}{9}+\frac{4z^4}{81}+i\frac{z^5}{27}+O(z^6).
\]

## What remains open

So the surviving PDE-facing question is no longer “is there some deformed outlet?”
It is much sharper:

> Does the actual moving-throat core realize the concrete coupling-balance surface and the auxiliary D/N-tube normalization identified above?

That is a genuinely microscopic throat-core question, not more reduced outlet algebra.

=== moving_throat_pde_stage101_parent_core_extraction.md ===


# Moving-Throat PDE — Stage 101: Parent-Action Extraction of Core Parameters

## Goal

Replace the reduced core parameters
\[
(K_s,\;K_q,\;\lambda,\;g_s,\;g_q)
\]
by explicit overlap formulas from one concrete **GNLS + localized-Maxwell throat-core ansatz**.

The point is not yet to prove the full moving-throat PDE. The point is to reduce the surviving outlet-core ambiguity to parent-action quantities that can be computed on an explicit throat branch.

## Parent sectors used

The parent 4D stack already gives:

- a gauged 4D GNLS matter sector with confinement \(V_{\rm conf}(\mathbf X;a,L)\) and frozen \(n=5\) EOS, so the bulk matter Hamiltonian has both gradient and compressional terms, fileciteturn33file4
- a localized \(4+1\) Maxwell sector with real mixed channels \((A_w,F_{\mu w},J^w)\) outside the strict far-field zero-mode reduction, fileciteturn33file1turn33file3
- and explicit mixed-sector EM energy channels such as \(E_w\), \(C_a=F_{aw}\), and \(S_{\rm EM}^w\), so a concrete throat-core mixed outlet is on-model rather than ad hoc. fileciteturn34file7turn34file12

## Core ansatz

Take one grouped-real \(P_2\) lane and introduce two core amplitudes:

- \(s(t)\): shell / compliance mode,
- \(q(t)\): mixed \(A_w/F_{\mu w}\)-type side-channel mode.

Use the factorized parent ansatz
\[
\delta\rho(\mathbf X,t)=s(t)\,\varrho_s(\mathbf X),
\qquad
\delta A_w(\mathbf X,t)=q(t)\,\mathcal A_q(\mathbf X),
\]
with all other fluctuating fields omitted at this first closure level.

For the shell mode, use the carried thin-wall GNLS profile
\[
\varrho_s(\mathbf X)=\rho_w\,\chi_s(y)\,Y_{2m}(\Omega),
\qquad
y:=\frac{r-a}{\ell},
\qquad
\chi_s(y)=f'(y),
\]
and on the canonical tanh wall,
\[
f(y)=\frac{1+\tanh y}{2},
\qquad
\chi_s(y)=\frac12\operatorname{sech}^2 y.
\]

For the mixed side-channel, use the first finite D/N half-wave on an auxiliary tube \(z\in[0,L_W]\):
\[
\mathcal A_q(\mathbf X)=\Xi_q(\Omega,w)\,\chi_{1/2}(z),
\qquad
\chi_{1/2}(z)=\sqrt{\frac{2}{L_W}}\sin\!\left(\frac{\pi z}{2L_W}\right).
\]

## 1. Shell stiffness \(K_s\)

Expanding the GNLS matter Hamiltonian to quadratic order in the shell mode gives
\[
\delta^2 H_s=\frac12 K_s s^2,
\]
with
\[
\boxed{
K_s
=
\int d^4X\left[
H_w\,\chi_s^2
+\frac{\hbar^2}{4m_\psi\rho_w}\,|\nabla\chi_s|^2
\right],
}
\]
where
\[
H_w=h'(\rho_w)=\frac{m_\psi c_{s,w}^2}{\rho_w}
\]
is the local compressional stiffness.

On the canonical thin wall, using the already carried moments
\[
I_f=\int_{-\infty}^{\infty}\chi_s^2\,dy=\frac13,
\qquad
I_g=\int_{-\infty}^{\infty}(\chi_s')^2\,dy=\frac4{15},
\]
this becomes
\[
\boxed{
K_s
=
4\pi a^2\left(
\frac{H_w\ell}{3}
+
\frac{\hbar^2}{15m_\psi\rho_w\,\ell}
\right).
}
\]

If the wall thickness is locked by the static GNLS healing relation
\[
\ell=\frac{\hbar}{2m_\psi c_{s,w}},
\]
then
\[
\boxed{
K_s
=
\frac{3\pi a^2\hbar^2}{5m_\psi\rho_w\,\ell}.
}
\]

## 2. Mixed side-channel stiffness \(K_q\)

Reducing the localized Maxwell action on the D/N half-wave gives
\[
\delta^2 H_q=\frac12 K_q q^2
\]
with
\[
\boxed{
K_q
=
\frac{\mathcal Z_q}{\mu_0}\int_0^{L_W} (\chi'_{1/2})^2\,dz,
}
\]
where
\[
\mathcal Z_q:=\int d\Omega\,dw\,Z(w)\,\Xi_q(\Omega,w)^2
\]
is the transverse localization norm of the mixed channel.

Since
\[
\int_0^{L_W}\chi_{1/2}^2\,dz=1,
\qquad
\int_0^{L_W}(\chi'_{1/2})^2\,dz=\frac{\pi^2}{4L_W^2},
\]
one gets
\[
\boxed{
K_q
=
\frac{\mathcal Z_q}{\mu_0}\,\frac{\pi^2 c_s^2}{4L_W^2}.
}
\]

So \(K_q\) is not arbitrary: it is fixed by one localization norm and the D/N tube length.

## 3. Static/mixed hybridization \(\lambda\)

The GNLS kinetic energy in Madelung form contains
\[
\mathcal H_{\rm kin}=\frac{m_\psi}{2}\rho\,v_Av_A.
\]
Take a static background mixed flow \(v_{w0}\neq0\), hold the phase fixed in the \(q\)-channel so that
\[
\delta v_w = -\frac{q_*}{m_\psi}\,\delta A_w,
\]
and expand
\[
\rho=\rho_0+s\,\varrho_s,
\qquad
A_w=A_{w0}+q\,\mathcal A_q.
\]
The bilinear \(sq\) term is exactly
\[
\delta^2 H_{sq}
=
- q_* v_{w0}\int d^4X\,\varrho_s\,\mathcal A_q\;s q.
\]
So the concrete hybridization is
\[
\boxed{
\lambda
=
- q_* v_{w0}\,\mathcal I_{sq},
\qquad
\mathcal I_{sq}:=\int d^4X\,\varrho_s\,\mathcal A_q.
}
\]

In the simplest uniform-core overlap closure,
\[
\mathcal I_{sq}=J_s I_q,
\qquad
J_s:=4\pi a^2\ell I_f=\frac{4\pi a^2\ell}{3},
\qquad
I_q:=\int_0^{L_W}\chi_{1/2}(z)\,dz=\frac{2\sqrt{2L_W}}{\pi},
\]
so
\[
\boxed{
\lambda
=
-\frac{8\sqrt2}{3}\,q_* v_{w0}\,a^2\ell\,\sqrt{L_W}.
}
\]

This is the first explicit parent-level source of the Stage-97 bilinear shell/mixed coupling.

## 4. Mouth couplings \(g_s\) and \(g_q\)

Use one external mouth amplitude \(u\) and let it couple to:

- shell traction at the mouth,
- mixed flux at the mouth.

### Shell coupling

If the mouth forcing is a lip traction \(\mathcal T_m u\), then the shell overlap gives
\[
\delta H_{\rm mouth}^{(s)}=-u g_s s,
\qquad
\boxed{
g_s=\mathcal T_m J_s=\mathcal T_m \frac{4\pi a^2\ell}{3}.
}
\]

### Mixed coupling

For the D/N mixed tube, the natural mouth observable is the derivative at the Dirichlet mouth:
\[
\delta H_{\rm mouth}^{(q)}=-u g_q q,
\qquad
g_q=\frac{\mathcal Z_q}{\mu_0}\chi'_{1/2}(0).
\]
Since
\[
\chi'_{1/2}(0)=\frac{\pi}{\sqrt2\,L_W^{3/2}},
\]
one gets
\[
\boxed{
g_q=\frac{\mathcal Z_q}{\mu_0}\,\frac{\pi}{\sqrt2\,L_W^{3/2}}.
}
\]

So every parameter in the Stage-97 core matrix now has a concrete parent-action interpretation.

## Result

The effective core matrix
\[
\begin{pmatrix}
K_s & \lambda\\
\lambda & -K_q D_W^{\rm bare}(z)
\end{pmatrix}
\binom{s}{q}
=
u\binom{g_s}{g_q}
\]
is no longer just a reduced model.

In the first explicit throat-core closure its entries are:

\[
\boxed{
K_s
=
4\pi a^2\left(
\frac{H_w\ell}{3}
+
\frac{\hbar^2}{15m_\psi\rho_w\,\ell}
\right),
}
\]
\[
\boxed{
K_q
=
\frac{\mathcal Z_q}{\mu_0}\,\frac{\pi^2 c_s^2}{4L_W^2},
}
\]
\[
\boxed{
\lambda
=
- q_* v_{w0}\,\mathcal I_{sq},
}
\]
\[
\boxed{
g_s=\mathcal T_m \frac{4\pi a^2\ell}{3},
\qquad
g_q=\frac{\mathcal Z_q}{\mu_0}\,\frac{\pi}{\sqrt2\,L_W^{3/2}}.
}
\]

The next step is to rewrite the compensation surface entirely in terms of the parent overlap ratios these formulas define.

=== moving_throat_pde_stage102_parent_balance_family.md ===


# Moving-Throat PDE — Stage 102: One-Parameter Parent Compensation Family

## Goal

Use the explicit parent-action formulas from Stage 101 to rewrite the compensated canonical outlet branch as a small set of dimensionless microscopic balance conditions.

## 1. Dimensionless parent ratios

Define the two dimensionless parent ratios
\[
\boxed{
\mathfrak r := \frac{\lambda}{\sqrt{K_sK_q}},
\qquad
\mathfrak g := \frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}}.
}
\]

Then the exact Stage-98 core-balance theorem
\[
g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2
\]
collapses to
\[
\boxed{
1+\mathfrak r^2 = 4(\mathfrak g-\mathfrak r)^2.
}
\]

So the compensated canonical branch is not a five-parameter tuning problem. It is a relation between just two normalized parent ratios.

## 2. Exact mouth-coupling balance law

Solving for \(\mathfrak g\) gives the exact two-branch law
\[
\boxed{
\mathfrak g
=
\mathfrak r
\pm
\frac12\sqrt{1+\mathfrak r^2}.
}
\]

This is the parent-action version of the Stage-98 coupling-balance surface.

Interpretation:

- \(\mathfrak r\) is the normalized static/mixed hybridization,
- \(\mathfrak g\) is the normalized relative mouth coupling of the mixed tube to the shell mode.

Once \(\mathfrak r\) is fixed by the core background, the required mouth-coupling ratio is fixed exactly.

## 3. D/N tube selection in parent variables

The Stage-98 D/N-tube condition
\[
\kappa_c=\frac13,
\qquad
r_c=\frac{\lambda^2}{K_sK_q}
\]
is simply
\[
r_c=\mathfrak r^2.
\]
So the auxiliary mixed-tube length becomes
\[
\boxed{
L_W
=
\frac{\pi a}{2}\sqrt{\frac{1+\mathfrak r^2}{3}}.
}
\]

This is a strong reduction of the geometric freedom: the D/N mixed-tube length is controlled by the **same** normalized hybridization that enters the coupling-balance law.

## 4. Parent formulas for \(\mathfrak r\) and \(\mathfrak g\)

From Stage 101,
\[
\mathfrak r
=
\frac{-q_* v_{w0}\mathcal I_{sq}}{\sqrt{K_sK_q}}.
\]

Under the simplest uniform-core overlap closure,
\[
\mathcal I_{sq}
=
\frac{8\sqrt2}{3}a^2\ell\sqrt{L_W},
\]
so
\[
\boxed{
\mathfrak r
=
-\frac{8\sqrt2}{3}\,
\frac{q_* v_{w0}\,a^2\ell\sqrt{L_W}}{\sqrt{K_sK_q}}.
}
\]

Likewise,
\[
g_s=\mathcal T_m J_s,
\qquad
J_s=\frac{4\pi a^2\ell}{3},
\qquad
g_q=\frac{\mathcal Z_q}{\mu_0}\frac{\pi}{\sqrt2\,L_W^{3/2}}.
\]
Using
\[
K_q=\frac{\mathcal Z_q}{\mu_0}\frac{\pi^2 c_s^2}{4L_W^2},
\]
the normalized mouth ratio becomes
\[
\boxed{
\mathfrak g
=
\frac{\sqrt{2\mathcal Z_q K_s}}
{\mathcal T_m J_s\,c_s\sqrt{\mu_0 L_W}}.
}
\]

So the compensated canonical outlet is selected by just two physical background controls:

- the normalized mixed background flow \(\mathfrak r\),
- the normalized mouth-traction strength \(\mathfrak g\).

## 5. Exact traction law

Solving the balance relation for the required traction amplitude gives
\[
\boxed{
\mathcal T_m
=
\frac{\sqrt{2\mathcal Z_q K_s}}
{J_s\,c_s\sqrt{\mu_0 L_W}}
\,
\frac{1}{\mathfrak r \pm \frac12\sqrt{1+\mathfrak r^2}}.
}
\]

So once the background mixed flow fixes \(\mathfrak r\), the traction amplitude required to preserve the canonical outgoing quadrupole fingerprint is fixed exactly.

## 6. Family-1 / healing-lock shell simplification

If the shell is on the carried healing-locked GNLS branch,
\[
\ell=\frac{\hbar}{2m_\psi c_{s,w}},
\]
then Stage 101 gave
\[
K_s=\frac{3\pi a^2\hbar^2}{5m_\psi\rho_w\,\ell}.
\]
So the only shell-side microscopic inputs left are \((a,\ell,\rho_w)\). All other shell data collapse into fixed prefactors.

## Result

The compensated canonical outlet has now been reduced to a **one-parameter parent family**:

1. choose the normalized mixed background flow \(\mathfrak r\),
2. then the D/N mixed-tube length is fixed by
   \[
   L_W=\frac{\pi a}{2}\sqrt{\frac{1+\mathfrak r^2}{3}},
   \]
3. and the required normalized mouth coupling is fixed by
   \[
   \mathfrak g=\mathfrak r\pm \frac12\sqrt{1+\mathfrak r^2}.
   \]

So the remaining PDE-facing question is no longer a large outlet-coefficient search. It is whether the actual GNLS + localized-Maxwell throat core picks one of these parent-balance branches.

=== moving_throat_pde_stage103_core_parameter_status.md ===


# Moving-Throat PDE — Stage 103: Core-Parameter Extraction Status

## Result

The concrete throat-core outlet is no longer described by free reduced coefficients.

After explicit GNLS + localized-Maxwell reduction, the surviving microscopic controls are:

- shell stiffness \(K_s\),
- mixed-tube stiffness \(K_q\),
- normalized shell/mixed hybridization
  \[
  \mathfrak r=\frac{\lambda}{\sqrt{K_sK_q}},
  \]
- normalized mouth-coupling ratio
  \[
  \mathfrak g=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}},
  \]
- and the D/N mixed-tube length \(L_W\).

## The sharp surviving theorem gate

The compensated canonical outgoing quadrupole branch exists **iff**
\[
1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2,
\]
with
\[
L_W=\frac{\pi a}{2}\sqrt{\frac{1+\mathfrak r^2}{3}}.
\]

So the surviving microscopic question is no longer:

> “What are the four reduced outlet coefficients?”

It is now:

> “What branch values of \((\mathfrak r,\mathfrak g)\) does the actual GNLS + localized-Maxwell throat core select?”

That is a substantially smaller target for the next derivation.

=== moving_throat_pde_stage108_positive_source_theorem.md ===


# Moving-Throat PDE — Stage 108: Positive Local Mouth-Source Theorem

## Goal

Replace the still-free normalized mouth-coupling ratio
\[
\mathfrak g=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}}
\]
by an explicit **localized mouth-source law** and determine which compensation branch is physically admissible.

## Localized source law

Take one scalar mouth amplitude \(u\) and one nonnegative normalized axial source profile
\[
\sigma(z)\ge 0,
\qquad
\int_0^L \sigma(z)\,dz = 1,
\qquad
z\in[0,L].
\]
Assume the shell mouth load uses the same positive source density throughout the mouth layer, while the mixed channel couples through the first D/N half-wave derivative on the same throat interval,
\[
\chi_{1/2}(z)=\sqrt{\frac{2}{L}}\sin\!\left(\frac{\pi z}{2L}\right),
\qquad
\chi_{1/2}'(z)=\sqrt{\frac{2}{L}}\frac{\pi}{2L}\cos\!\left(\frac{\pi z}{2L}\right).
\]

After normalizing to the point-source branch \(\sigma(z)=\delta(z)\), the mouth-bias factor becomes
\[
\boxed{
\mathfrak g[\sigma]
=
\int_0^L \sigma(z)\cos\!\left(\frac{\pi z}{2L}\right)\,dz.
}
\]

So the mouth-source bias is the first cosine moment of the positive axial source profile.

## Exact positivity theorem

Because
\[
0\le \cos\!\left(\frac{\pi z}{2L}\right)\le 1
\qquad \text{for } z\in[0,L],
\]
every positive normalized source law satisfies
\[
\boxed{
0\le \mathfrak g[\sigma]\le 1.
}
\]

This immediately selects the physically admissible compensation branch.

From Stage 105, on the explicit Family-1 branch,
\[
\mathfrak g_-^{F1}\approx 0.758035078944663,
\qquad
\mathfrak g_+^{F1}\approx 2.79795199200529.
\]
Therefore
\[
\boxed{
\mathfrak g_+^{F1}>1
\quad\Rightarrow\quad
\text{the upper compensated branch is impossible for any positive localized source law.}
}
\]

By contrast,
\[
\boxed{
0<\mathfrak g_-^{F1}<1
\quad\Rightarrow\quad
\text{the lower compensated branch is the unique physically admissible canonical branch under positive mouth sourcing.}
}
\]

## Interpretation

This is a strong narrowing.

The outlet-core ambiguity is no longer “which branch do we choose?” In any localized positive source model on the first D/N throat interval:

- the upper compensated branch is ruled out,
- the lower compensated branch is the only canonical candidate.

So the remaining question is not branch sign, but the **shape** of the positive source profile \(\sigma(z)\).

=== moving_throat_pde_stage109_positive_source_families.md ===


# Moving-Throat PDE — Stage 109: Explicit Positive Source Families and the Family-1 Compensation Point

## Goal

Evaluate explicit positive localized source laws and determine whether the physical lower compensated branch can be reached without exotic sign-changing mouth profiles.

## 1. Self-matched derivative profile

The most natural positive axial source on the first D/N interval is the normalized derivative profile itself,
\[
\boxed{
\sigma_{\rm match}(z)=k\cos(kz),
\qquad
k=\frac{\pi}{2L},
}
\]
since
\[
\int_0^L k\cos(kz)\,dz = 1.
\]

Its mouth-bias factor is exact:
\[
\mathfrak g_{\rm match}
=
\int_0^L \sigma_{\rm match}(z)\cos(kz)\,dz
=
k\int_0^L \cos^2(kz)\,dz
=
\frac{\pi}{4}.
\]
So
\[
\boxed{
\mathfrak g_{\rm match}=\frac{\pi}{4}\approx 0.785398163397448.
}
\]

This already lies much closer to the physical lower compensated branch than the naive point-source value \(\mathfrak g=1\).

## 2. Family-1 comparison

On the explicit geometric branch from Stages 104–105,
\[
\mathfrak g_-^{F1}\approx 0.758035078944663.
\]
Therefore the self-matched source overshoots the exact compensated lower branch by
\[
\Delta\mathfrak g_{\rm match}
=
\frac{\pi}{4}-\mathfrak g_-^{F1}
\approx 0.0273630844527852.
\]

Since traction scales as \(\mathfrak g^{-1}\),
\[
\boxed{
\frac{\mathcal T_m^{(-)}}{\mathcal T_m^{\rm match}}
=
\frac{\pi/4}{\mathfrak g_-^{F1}}
\approx 1.036097385480999.
}
\]
So the exact Family-1 compensated branch is only a **3.61% traction enhancement** away from the self-matched derivative profile.

## 3. Exact positive family interpolating through the compensated point

Introduce the convex positive family
\[
\boxed{
\sigma_\xi(z)
=
(1-\xi)\,k\cos(kz)+\xi\,\frac1L,
\qquad 0\le \xi\le 1.
}
\]
This is normalized and nonnegative on \([0,L]\). Its mouth-bias factor is
\[
\mathfrak g_\xi
=
(1-\xi)\frac{\pi}{4}+\xi\frac{2}{\pi}.
\]

Because
\[
\frac{2}{\pi}<\mathfrak g_-^{F1}<\frac{\pi}{4},
\]
there is a unique \(\xi_*\in(0,1)\) such that
\[
\mathfrak g_{\xi_*}=\mathfrak g_-^{F1}.
\]
The exact solution is
\[
\boxed{
\xi_*
=
\frac{\frac{\pi}{4}-\mathfrak g_-^{F1}}{\frac{\pi}{4}-\frac{2}{\pi}}
=
\frac{-37\sqrt3-5\pi^2+2\sqrt{4107-100\pi^2}}{5(8-\pi^2)}
\approx 0.183918405511538.
}
\]

So only an **18.4% admixture** of the fully washed positive profile \(1/L\) into the self-matched derivative profile reaches the exact lower compensated branch.

## Result

The explicit mouth-source bias is now much narrower than before:

- the point-source branch \(\mathfrak g=1\) is too high,
- the exact self-matched D/N derivative source gives \(\mathfrak g=\pi/4\),
- and the true Family-1 compensated branch is reached by a small \(18.4\%\) broadening of that already-natural positive source law.

So the canonical branch no longer looks like a delicate coefficient fit. It sits inside a simple exact family of positive localized mouth sources.

=== moving_throat_pde_stage110_penetration_families.md ===


# Moving-Throat PDE — Stage 110: Geometric Mouth-Penetration Families

## Goal

Translate the lower compensated mouth-source bias into simple geometric penetration scales.

## 1. Uniform slab source

Take the positive localized source profile
\[
\boxed{
\sigma_x^{\rm slab}(z)=\frac{1}{xL},
\qquad
0\le z\le xL,
\qquad
0<x\le 1,
}
\]
and zero elsewhere on \([0,L]\).

Then
\[
\mathfrak g_{\rm slab}(x)
=
\int_0^{xL}\frac{dz}{xL}\cos\!\left(\frac{\pi z}{2L}\right)
=
\boxed{
\frac{2}{\pi x}\sin\!\left(\frac{\pi x}{2}\right).
}
\]

Solving
\[
\mathfrak g_{\rm slab}(x)=\mathfrak g_-^{F1}
\]
gives the unique Family-1 slab depth
\[
\boxed{
x_*^{\rm slab}\approx 0.797839360904564.
}
\]

So a positive uniform mouth source extending over about \(80\%\) of the throat span reaches the lower compensated branch exactly.

## 2. Truncated exponential source

Take the normalized positive exponential family
\[
\boxed{
\sigma_x^{\exp}(z)
=
\frac{e^{-z/(xL)}}{xL\left(1-e^{-1/x}\right)},
\qquad
z\in[0,L],
\qquad
x>0.
}
\]

Its mouth-bias factor is
\[
\boxed{
\mathfrak g_{\exp}(x)
=
\frac{2\left(2+\pi x\,e^{-1/x}\right)}
{\left(4+\pi^2x^2\right)\left(1-e^{-1/x}\right)}.
}
\]

Solving
\[
\mathfrak g_{\exp}(x)=\mathfrak g_-^{F1}
\]
gives
\[
\boxed{
x_*^{\exp}\approx 0.662765402623160.
}
\]

So an exponentially localized mouth source with penetration depth about \(0.66L\) reaches the exact lower compensated branch.

## Result

These two geometric positive-source families both hit the same lower compensated branch with **moderate** mouth penetration depths:
\[
x_*^{\exp}\approx 0.66,
\qquad
x_*^{\rm slab}\approx 0.80.
\]

Combined with Stages 108–109, this means:

- the upper compensated branch is excluded by positivity,
- the lower compensated branch is the unique physical candidate,
- and simple positive localized source laws reach it without requiring sign-changing or finely oscillatory mouth forcing.

=== moving_throat_pde_stage111_mouth_source_bias_status.md ===


# Moving-Throat PDE — Stage 111: Mouth-Source Bias Status

## What is now fixed

The outlet-core branch-selection problem has narrowed further.

1. On the explicit Family-1 throat geometry,
\[
\mathfrak r_{F1}\approx 1.77799353547498
\]
is already fixed by geometry.

2. For any **positive normalized localized mouth source** on the first D/N interval,
\[
\mathfrak g[\sigma]
=
\int_0^L \sigma(z)\cos\!\left(\frac{\pi z}{2L}\right)\,dz
\]
satisfies
\[
0\le \mathfrak g\le 1.
\]

3. Therefore the upper compensated branch
\[
\mathfrak g_+^{F1}\approx 2.79795199200529
\]
is ruled out by positivity alone.

4. The lower compensated branch
\[
\mathfrak g_-^{F1}\approx 0.758035078944663
\]
is the unique physically admissible canonical branch under positive mouth sourcing.

## What the explicit source laws show

- Point source:
  \[
  \mathfrak g=1
  \]
  (the old naive equal-normalized branch).

- Self-matched D/N derivative source:
  \[
  \mathfrak g=\frac{\pi}{4}\approx 0.785398
  \]
  which is already only \(3.61\%\) away in traction from exact lower-branch compensation.

- Exact convex positive family:
  \[
  \sigma_\xi(z)=(1-\xi)\,k\cos(kz)+\xi/L
  \]
  hits the exact lower compensated branch at
  \[
  \xi_*\approx 0.183918.
  \]

- Uniform slab and truncated exponential penetration families also hit the exact lower branch at moderate penetration depths:
  \[
  x_*^{\rm slab}\approx 0.797839,
  \qquad
  x_*^{\exp}\approx 0.662765.
  \]

## Meaning

The branch-selection ambiguity is no longer “does some unknown source law maybe rescue the canonical outlet?”

It is now much sharper:

\[
\boxed{
\text{Under positive localized mouth sourcing, the lower compensated branch is uniquely selected and is easily reachable.}
}
\]

So the remaining PDE-facing question is not branch sign. It is the detailed shape of the actual mouth source profile \(\sigma(z)\), or equivalently the exact amount of positive mouth broadening away from the point-source limit.

=== moving_throat_pde_stage112_mouth_boundary_layer.md ===


# Moving-Throat PDE — Stage 112: Explicit GNLS + Localized-Maxwell Mouth Boundary Layer

## Goal

Replace the remaining abstract positive mouth-source family by the first explicit
boundary-layer law derived from a GNLS + localized-Maxwell electrochemical balance
at the mouth.

The point is not to solve the full moving-throat core. The point is to derive an
honest **local mouth source law** from the same parent structures already frozen
in the 4D stack:
- entropic/compressible source thermodynamics from the GNLS matter sector,
- a localized Maxwell scalar potential \(A_0\),
- and the mechanical mouth-traction/confinement channel already carried by
  \(\delta V_{\rm conf}\).

---

## 1. Local mouth free energy and electrochemical potential

Let \(z\in[0,L]\) denote the axial coordinate measured inward from the mouth.

Use the minimal positive source-density free energy
\[
\boxed{
F_{\rm mouth}[\sigma]
=
\int_0^L dz\,
\Big[
\Theta_\sigma\,\sigma\!\left(\ln\frac{\sigma}{\sigma_*}-1\right)
+
V_m(z)\,\sigma
\Big],
}
\]
with
\[
\Theta_\sigma>0
\]
the effective source compressibility/entropic scale.

Near the mouth, linearize the parent potential load as
\[
\boxed{
V_m(z)\approx V_1 z,
\qquad
V_1:=\left.\partial_z\delta V_{\rm conf}\right|_{\rm m}
-
q_*\left.\partial_z A_0\right|_{\rm m}.
}
\]
So \(V_1\) is the **net mouth-localizing slope** coming from the mechanical mouth
traction and the localized-Maxwell scalar-potential drop.

The corresponding electrochemical potential is
\[
\mu_\sigma^{\rm chem}(z)
=
\Theta_\sigma\ln\!\frac{\sigma}{\sigma_*}
+
V_1 z.
\]

---

## 2. Onsager current and stationary zero-flux branch

Take the standard positive Onsager current
\[
J_\sigma
=
-M_\sigma\,\sigma\,\partial_z\mu_\sigma^{\rm chem},
\qquad
M_\sigma>0.
\]
Then
\[
J_\sigma
=
-M_\sigma\left(\Theta_\sigma\,\sigma'(z)+V_1\sigma(z)\right).
\]

On the stationary recirculating mouth branch,
\[
\partial_t \sigma + \partial_z J_\sigma = 0,
\qquad
J_\sigma=0,
\]
so the source law is exactly
\[
\Theta_\sigma\,\sigma'(z)+V_1\sigma(z)=0.
\]

Therefore
\[
\sigma(z)=C\,e^{-V_1 z/\Theta_\sigma}.
\]

Normalizing to one total mouth source on \([0,L]\),
\[
\int_0^L dz\,\sigma(z)=1,
\]
gives the exact positive family
\[
\boxed{
\sigma_{\Pi}(z)
=
\frac{\Pi\,e^{-\Pi z/L}}{L\left(1-e^{-\Pi}\right)},
\qquad
\Pi:=\frac{V_1L}{\Theta_\sigma}>0.
}
\]

So the previously ad hoc truncated exponential family is now the exact zero-flux
equilibrium branch of a GNLS + localized-Maxwell mouth boundary layer.

---

## 3. Physical meaning of \(\Pi\)

The single remaining source-shape parameter is now

\[
\boxed{
\Pi_m=\frac{V_1L}{\Theta_\sigma}
=
\frac{L}{\Theta_\sigma}
\left(
\left.\partial_z\delta V_{\rm conf}\right|_{\rm m}
-
q_*\left.\partial_z A_0\right|_{\rm m}
\right).
}
\]

Interpretation:

- \(\Pi_m\to 0\): almost uniform source along the throat;
- \(\Pi_m\to\infty\): point-like source concentrated at the mouth;
- finite \(\Pi_m\): positive localized penetration into the throat interior.

So the mouth-source ambiguity has collapsed from a free profile \(\sigma(z)\) to a
single **dimensionless electrochemical bias** \(\Pi_m\).

---

## 4. What this changes

The branch-selection problem is now sharper:

1. positivity already killed the upper compensated branch;
2. the lower compensated branch was already known to be reachable by truncated
   exponentials;
3. now that truncated exponential is no longer just a convenient family — it is
   the exact equilibrium law of a concrete mouth boundary-layer model.

So the next question is no longer “what source family should we try?” It is simply:

\[
\boxed{
\text{what value of }\Pi_m\text{ does the actual moving-throat mouth layer select?}
}

=== moving_throat_pde_stage113_mouth_bias_map.md ===


# Moving-Throat PDE — Stage 113: Exact Mouth-Bias Map and Family-1 Compensation Point

## Goal

Compute the exact Family-1 mouth-bias factor for the explicit boundary-layer
profile
\[
\sigma_{\Pi}(z)=\frac{\Pi e^{-\Pi z/L}}{L(1-e^{-\Pi})},
\]
and solve the canonical lower-branch compensation condition on the explicit
Family-1 geometry.

---

## 1. Exact mouth-bias factor

The positive mouth-source branch is measured against the first D/N derivative
shape
\[
\cos\!\left(\frac{\pi z}{2L}\right).
\]
So the explicit mouth-bias factor is
\[
\mathfrak g_{\Pi}
=
\int_0^L dz\,
\sigma_{\Pi}(z)\,
\cos\!\left(\frac{\pi z}{2L}\right).
\]

This evaluates exactly to
\[
\boxed{
\mathfrak g_{\Pi}
=
\frac{2\Pi\bigl(2\Pi e^{\Pi}+\pi\bigr)}
{\left(4\Pi^2+\pi^2\right)\left(e^{\Pi}-1\right)}.
}
\]

Equivalent \(x=1/\Pi\) form:
\[
\mathfrak g_{\Pi}
=
\frac{2\left(2+\pi x\,e^{-1/x}\right)}
{\left(4+\pi^2x^2\right)\left(1-e^{-1/x}\right)},
\qquad x=\frac1\Pi,
\]
which matches the earlier truncated-exponential penetration family exactly.

---

## 2. Exact monotonicity law

Because \(\sigma_\Pi(z)\) is an exponential family,
\[
\partial_\Pi \sigma_\Pi
=
-\frac{1}{L}\,\sigma_\Pi
\Bigl(z-\langle z\rangle_\Pi\Bigr),
\]
so
\[
\boxed{
\frac{d\mathfrak g_\Pi}{d\Pi}
=
-\frac{1}{L}\,
\mathrm{Cov}_\Pi\!\left(
\cos\!\frac{\pi z}{2L},
\,z
\right).
}
\]

Since \(\cos(\pi z/2L)\) is strictly decreasing on \([0,L]\), the covariance is
negative, hence
\[
\boxed{
\frac{d\mathfrak g_\Pi}{d\Pi}>0.
}
\]

So the explicit GNLS + localized-Maxwell boundary-layer family is strictly
monotone:
\[
\mathfrak g_{\Pi}: \frac{2}{\pi}\longrightarrow 1
\qquad (\Pi:0^+\to+\infty).
\]

---

## 3. Family-1 compensation point

The explicit Family-1 lower compensated branch is already fixed from the earlier
core balance analysis:
\[
\mathfrak g_-^{F1}\approx 0.758035078944663.
\]

Solving
\[
\mathfrak g_{\Pi}=\mathfrak g_-^{F1}
\]
gives the unique explicit compensation point
\[
\boxed{
\Pi_* \approx 1.50882951349316.
}
\]

Equivalently,
\[
\boxed{
x_*=\frac{1}{\Pi_*}\approx 0.662765402623160,
}
\]
which is exactly the penetration-depth value already found earlier for the
truncated exponential family.

So the compensated Family-1 mouth source is selected by a **moderate**
electrochemical bias:
\[
\Pi_*\sim O(1).
\]

---

## 4. Interpretation

This is a sharper statement than the earlier positive-family result.

It now says:

- the lower compensated branch is not only positive and reachable;
- it is the unique target of a concrete monotone GNLS + localized-Maxwell
  boundary-layer family;
- and the required bias is moderate rather than extreme.

So the remaining branch-selection ambiguity has collapsed to one explicit number:

\[
\boxed{
\Pi_m \stackrel{?}{\approx} 1.50882951349.
}
\]

=== moving_throat_pde_stage114_parent_mouth_threshold.md ===


# Moving-Throat PDE — Stage 114: Parent Micro-Threshold for Canonical Mouth Compensation

## Goal

Translate the explicit Family-1 compensation point
\[
\Pi_* \approx 1.50882951349316
\]
into a direct parent-level threshold on the localized Maxwell + confinement
mouth data.

---

## 1. The parent bias parameter

From Stage 112,
\[
\Pi_m=\frac{V_1L}{\Theta_\sigma},
\qquad
V_1=
\left.\partial_z\delta V_{\rm conf}\right|_{\rm m}
-
q_*\left.\partial_zA_0\right|_{\rm m}.
\]

So the explicit canonical branch condition is simply
\[
\boxed{
\Pi_m=\Pi_*.
}
\]

Equivalently,
\[
\boxed{
V_1
=
\frac{\Pi_*\,\Theta_\sigma}{L}
\approx
1.50882951349316\,\frac{\Theta_\sigma}{L}.
}
\]

This is the first direct parent-level threshold for the mouth-source law.

---

## 2. Localized-Maxwell / mechanical split

Writing
\[
V_1=T_m-q_*A_0',
\qquad
T_m:=\left.\partial_z\delta V_{\rm conf}\right|_{\rm m},
\qquad
A_0':=\left.\partial_zA_0\right|_{\rm m},
\]
the canonical compensated branch requires
\[
\boxed{
T_m-q_*A_0'
=
1.50882951349316\,\frac{\Theta_\sigma}{L}.
}
\]

So the explicit outlet-core problem has collapsed again:

- \(T_m\) and \(A_0'\) may each vary,
- but only the **single linear combination**
  \[
  T_m-q_*A_0'
  \]
  matters for the mouth-source bias on this boundary-layer branch.

---

## 3. Linearized deviation law

Near the compensated point,
\[
\mathfrak g_{\Pi}
=
\mathfrak g_-^{F1}
+
\mathfrak g'(\Pi_*)(\Pi-\Pi_*)
+O((\Pi-\Pi_*)^2),
\]
with
\[
\mathfrak g'(\Pi_*)\approx 0.0714453558083195.
\]

Therefore
\[
\boxed{
\mathfrak g_{\Pi}-\mathfrak g_-^{F1}
\approx
0.0714453558083195\,(\Pi_m-\Pi_*).
}
\]

So the remaining source-shape uncertainty is now first-order equivalent to the
single parent bias mismatch
\[
\Pi_m-\Pi_*=
\frac{L}{\Theta_\sigma}\left(T_m-q_*A_0'\right)-\Pi_*.
\]

---

## 4. Meaning

The mouth-source problem is no longer:

- which sign branch,
- which family,
- or whether compensation is even reachable.

It is now just:

\[
\boxed{
\text{does the real mouth layer satisfy }
T_m-q_*A_0'
\approx
1.509\,\Theta_\sigma/L\ ?
}
\]

That is a concrete microscopic threshold, not a broad qualitative criterion.

=== moving_throat_pde_stage115_mouth_boundary_layer_status.md ===


# Moving-Throat PDE — Stage 115: Mouth Boundary-Layer Status After Explicit Source-Law Extraction

## What is now fixed

After Stages 112–114, the mouth-source side is no longer an abstract profile problem.

1. The actual positive source family is fixed by an explicit GNLS + localized-Maxwell
   boundary-layer law:
   \[
   \sigma_{\Pi}(z)=\frac{\Pi e^{-\Pi z/L}}{L(1-e^{-\Pi})}.
   \]

2. The corresponding Family-1 mouth-bias factor is exact:
   \[
   \mathfrak g_{\Pi}
   =
   \frac{2\Pi(2\Pi e^{\Pi}+\pi)}
        {(4\Pi^2+\pi^2)(e^{\Pi}-1)}.
   \]

3. This family is strictly monotone:
   \[
   \mathfrak g_{\Pi}: \frac{2}{\pi}\to 1
   \qquad(\Pi:0^+\to+\infty).
   \]

4. The unique canonical Family-1 compensation point is
   \[
   \Pi_* \approx 1.50882951349316.
   \]

5. The exact parent threshold is
   \[
   \left.\partial_z\delta V_{\rm conf}\right|_{\rm m}
   -
   q_*\left.\partial_zA_0\right|_{\rm m}
   =
   1.50882951349316\,\frac{\Theta_\sigma}{L}.
   \]

## Meaning

The outlet-core ambiguity has narrowed again.

The open question is no longer:

- which compensation branch is physical,
- whether positive localized sources can hit it,
- or what source family to use.

It is now simply:

\[
\boxed{
\text{what value of }
\Pi_m=\frac{L}{\Theta_\sigma}\left(
\left.\partial_z\delta V_{\rm conf}\right|_{\rm m}
-q_*\left.\partial_zA_0\right|_{\rm m}
\right)
\text{ does the real mouth layer select?}
}
\]

And the target value is moderate:
\[
\Pi_m\approx1.51.
\]

So the remaining gap is now one clean microscopic bias law, not a diffuse branch-selection problem.

=== moving_throat_pde_stage116_coupled_mouth_fixedpoint.md ===


# Moving-Throat PDE — Stage 116: Full Coupled Mouth-Layer Fixed-Point Law

## Goal

Replace the effective parent datum
\[
\Pi_m=\frac{L}{\Theta_\sigma}
\left(
\left.\partial_z\delta V_{{\rm conf}}\right|_{\rm m}
-
q_*\left.\partial_zA_0\right|_{\rm m}
\right)
\]
by the first explicit **coupled mouth-layer solve** in which the mechanical and
localized-Maxwell throat channels are solved together and the source law
\(\sigma_\Pi\) closes self-consistently.

The point is not yet to solve the full moving-throat PDE. The point is to reduce
the last mouth-bias ambiguity to a small set of explicit dimensionless gains.

---

## 1. Dimensionless coupled mouth-layer system

Let
\[
x:=\frac{z}{L}\in[0,1],
\qquad
\Sigma_\Pi(x)=\frac{\Pi e^{-\Pi x}}{1-e^{-\Pi}}
\]
be the normalized mouth-source law from Stage 112.

Introduce the two-component boundary-layer field
\[
\mathbf U(x)=
\begin{pmatrix}
U_s(x)\\
U_q(x)
\end{pmatrix},
\]
where:

- \(U_s\) is the shell / confinement-response field,
- \(U_q\) is the localized-Maxwell / mixed-channel field.

The first coupled linear mouth-layer model is
\[
\boxed{
\left[-\partial_x^2\,\mathbf I+\mathbf K\right]\mathbf U(x)
=
\Sigma_\Pi(x)\,\mathbf G,
}
\]
with D/N-type mouth/bottom conditions
\[
\boxed{
\mathbf U(0)=0,
\qquad
\mathbf U'(1)=0.
}
\]

Here \(\mathbf K\) is a constant positive symmetric \(2\times2\) stiffness matrix
and \(\mathbf G\) is the source-coupling vector.

The parent mouth bias is read from a linear projection of the mouth derivative,
\[
\boxed{
V_1 = \mathbf H^{\!T}\mathbf U'(0),
\qquad
\Pi=\frac{L}{\Theta_\sigma}V_1,
}
\]
for some channel-weight vector \(\mathbf H\).

So the source parameter \(\Pi\) is no longer an input. It is the fixed point of a
coupled mechanical/electromagnetic boundary-layer problem.

---

## 2. Exact diagonalization

Diagonalize
\[
\mathbf K = \mathbf R
\begin{pmatrix}
\kappa_+^2 & 0\\
0 & \kappa_-^2
\end{pmatrix}
\mathbf R^{T},
\qquad
\mathbf R\in O(2).
\]
Define eigenbasis couplings
\[
\mathbf G_\pm := \mathbf R^T\mathbf G,
\qquad
\mathbf H_\pm := \mathbf R^T\mathbf H.
\]

Then each eigenchannel \(u_\alpha(x)\) satisfies
\[
\left(-\partial_x^2+\kappa_\alpha^2\right)u_\alpha(x)
=
G_\alpha\,\Sigma_\Pi(x),
\qquad
u_\alpha(0)=0,\quad
u_\alpha'(1)=0,
\qquad
\alpha\in\{+,-\}.
\]

So the fully coupled two-channel mouth-layer problem reduces exactly to two scalar
D/N response problems.

---

## 3. Exact scalar D/N response to the exponential source

For one channel with stiffness \(\kappa\), solve
\[
\left(-\partial_x^2+\kappa^2\right)u(x)=G\,\Sigma_\Pi(x),
\qquad
u(0)=0,\quad u'(1)=0.
\]

A direct elementary solve gives
\[
u(x)=A\sinh(\kappa x)-C\cosh(\kappa x)+C e^{-\Pi x},
\]
with
\[
C=\frac{G\Pi}{(1-e^{-\Pi})(\kappa^2-\Pi^2)},
\qquad
A=
\frac{C\left(\kappa\sinh\kappa+\Pi e^{-\Pi}\right)}{\kappa\cosh\kappa}.
\]

Therefore the exact mouth derivative is
\[
u'(0)=G\,\mathcal S(\Pi,\kappa),
\]
where
\[
\boxed{
\mathcal S(\Pi,\kappa)
=
\frac{\Pi\left[\kappa\tanh\kappa+\Pi\left(e^{-\Pi}\operatorname{sech}\kappa-1\right)\right]}
{(1-e^{-\Pi})(\kappa^2-\Pi^2)}.
}
\]

This is the exact D/N mouth-response kernel for the exponential source branch.

A useful exact limit is
\[
\boxed{
\mathcal S(\Pi,0)=1,
}
\]
so a purely static shell-compliance lane contributes a constant unit slope factor.

---

## 4. Full fixed-point law

Returning to the coupled system,
\[
V_1
=
\sum_{\alpha=\pm} H_\alpha u_\alpha'(0)
=
\sum_{\alpha=\pm} H_\alpha G_\alpha\,\mathcal S(\Pi,\kappa_\alpha).
\]

So the explicit self-consistency law for the mouth-source bias is
\[
\boxed{
\Pi
=
\sum_{\alpha=\pm}
M_\alpha\,\mathcal S(\Pi,\kappa_\alpha),
\qquad
M_\alpha:=\frac{L\,H_\alpha G_\alpha}{\Theta_\sigma}.
}
\]

This is the first honest coupled replacement for the Stage-112 effective parent
datum.

The open mouth-source problem has now shrunk to:

1. the two D/N eigen-stiffnesses \(\kappa_\pm\),
2. the two dimensionless channel gains \(M_\pm\).

Everything else is fixed by the exact source law and the exact D/N response.

---

## Result

The actual mouth-layer selection problem is no longer “what profile \(\sigma(z)\)?”
and no longer “what effective slope \(V_1\)?”

It is now the explicit coupled fixed-point equation
\[
\boxed{
\Pi
=
M_+\mathcal S(\Pi,\kappa_+)+M_-\mathcal S(\Pi,\kappa_-).
}
\]

So the remaining ambiguity is now just a **small dimensionless gain/stiffness
quadruple** \((M_+,M_-,\kappa_+,\kappa_-)\).

=== moving_throat_pde_stage117_family1_mouth_fixedpoint.md ===


# Moving-Throat PDE — Stage 117: Family-1 Shell + First Mixed D/N Tube Reduction

## Goal

Specialize the general Stage-116 fixed-point law to the first explicit
moving-throat mouth-layer branch compatible with the carried Family-1 geometry:

- one **static shell-compliance** lane,
- one **first mixed D/N half-wave** lane.

This is the cleanest reduction of the coupled mouth-layer solve that still knows
about the localized-Maxwell mixed tube.

---

## 1. Channel identifications

Take the two diagonal D/N channels to be

- shell/compliance lane:
  \[
  \boxed{\kappa_s=0,}
  \]
  so it contributes the exact static limit \(\mathcal S(\Pi,0)=1\);

- mixed localized-Maxwell lane:
  \[
  \boxed{\kappa_q=\frac{\pi}{2},}
  \]
  the first D/N half-wave on the actual throat span.

Define the corresponding dimensionless gains
\[
M_s:=\frac{L\,H_sG_s}{\Theta_\sigma},
\qquad
M_q:=\frac{L\,H_qG_q}{\Theta_\sigma}.
\]

Then the coupled mouth-layer law collapses to
\[
\boxed{
\Pi = M_s + M_q\,\mathcal S_q(\Pi),
\qquad
\mathcal S_q(\Pi):=\mathcal S\!\left(\Pi,\frac{\pi}{2}\right).
}
\]

Using the exact Stage-116 kernel,
\[
\boxed{
\mathcal S_q(\Pi)
=
\frac{
\Pi\left[\frac{\pi}{2}\tanh\!\frac{\pi}{2}
+\Pi\left(e^{-\Pi}\operatorname{sech}\!\frac{\pi}{2}-1\right)\right]
}{
(1-e^{-\Pi})\left(\frac{\pi^2}{4}-\Pi^2\right)
}.
}
\]

So the Family-1 mouth-layer branch is now an explicit one-dimensional fixed-point
problem in the two gains \(M_s\) and \(M_q\).

---

## 2. Canonical compensation line

Stage 114 already fixed the exact source-shape compensation point
\[
\Pi_* \approx 1.50882951349316.
\]

Inserting this into the Family-1 fixed-point law gives the exact gain-line
condition
\[
\boxed{
M_s = \Pi_* - M_q\,\mathcal S_q(\Pi_*).
}
\]

Numerically,
\[
\boxed{
\mathcal S_q(\Pi_*) \approx 0.658075937605429.
}
\]
So the canonical Family-1 compensation line is
\[
\boxed{
M_s \approx 1.50882951349316 - 0.658075937605429\,M_q.
}
\]

Interpretation:

- if the mixed localized-Maxwell lane is mouth-localizing (\(M_q>0\)),
  the shell traction demand is reduced;
- if the mixed lane is de-localizing (\(M_q<0\)),
  the shell traction must be larger than the Stage-114 pure-mechanical threshold.

---

## 3. What has improved

Compared to Stage 114, the threshold is no longer written in terms of the effective
combination
\[
\left.\partial_z\delta V_{\rm conf}\right|_{\rm m}
-
q_*\left.\partial_zA_0\right|_{\rm m}
\]
alone.

It is now resolved into an explicit two-channel boundary-layer law with:

- a static shell lane,
- a first mixed D/N lane,
- and one exact nontrivial response factor
  \(\mathcal S_q(\Pi)\).

So the “actual \(\Pi_m\)” problem has narrowed from an effective slope datum to a
simple two-gain fixed-point problem.

---

## Result

On the first explicit Family-1 coupled mouth-layer branch, the actual source-shape
bias is determined by
\[
\boxed{
\Pi = M_s + M_q\,\mathcal S_q(\Pi),
}
\]
and the canonical compensated branch sits on the straight gain line
\[
\boxed{
M_s \approx 1.50882951349316 - 0.658075937605429\,M_q.
}
\]

The remaining ambiguity is therefore no longer profile selection and no longer a
free mouth slope: it is just the signed gain pair \((M_s,M_q)\).

=== moving_throat_pde_stage118_outlet_consistent_mouth_closure.md ===


# Moving-Throat PDE — Stage 118: Outlet-Consistent Mouth Closure

## Goal

Use the concrete two-channel throat-core compensation law from Stages 97–99 as a
first explicit closure for the Family-1 mouth-layer gains.

The point is to see what the coupled mouth-layer solve predicts when the
mouth-bias channels inherit the same shell-versus-mixed weighting that preserved
the canonical outgoing quadrupole branch.

---

## 1. Nontrivial compensated outlet ratio

The exact compensated core outlet found in Stage 98 was
\[
\delta\Lambda_{\rm core}(z)
=
4\sigma_*-\frac{\sigma_*}{1-z^2/3-i z^5/9}.
\]

So the nontrivial compensated outlet weights the static shell lane and the mixed
pole lane in the ratio
\[
\boxed{4 : -1.}
\]

As a first **outlet-consistent mouth closure**, impose the same ratio on the
dimensionless mouth-layer gains:
\[
\boxed{
M_s = 4\Sigma_m,
\qquad
M_q = -\Sigma_m,
}
\]
with \(\Sigma_m>0\) a single residual mouth-gain amplitude.

This is not yet a theorem of the full PDE. It is the cleanest way to make the
mouth-layer solve and the outlet-core compensation law talk to each other.

---

## 2. Reduced one-parameter fixed-point law

With
\[
M_s = 4\Sigma_m,
\qquad
M_q = -\Sigma_m,
\]
the Family-1 mouth-layer equation becomes
\[
\boxed{
\Pi
=
\Sigma_m\left[4-\mathcal S_q(\Pi)\right].
}
\]

So the actual source-shape bias is now selected by a **single** dimensionless
mouth gain \(\Sigma_m\).

Because \(\mathcal S_q(\Pi)\) is positive and bounded, the right-hand side stays
between \(3\Sigma_m\) and \(4\Sigma_m\), so this branch automatically lives at
moderate \(\Pi\) rather than extreme mouth-localization.

---

## 3. Canonical Family-1 value

Requiring the canonical compensation point
\[
\Pi=\Pi_*
\approx 1.50882951349316
\]
gives the exact outlet-consistent gain
\[
\boxed{
\Sigma_m^*
=
\frac{\Pi_*}{4-\mathcal S_q(\Pi_*)}
\approx 0.451485277739090.
}
\]

Therefore the corresponding shell and mixed gains are
\[
\boxed{
M_s^* = 4\Sigma_m^*
\approx 1.80594111095636,
\qquad
M_q^* = -\Sigma_m^*
\approx -0.451485277739090.
}
\]

So the canonical mouth-bias point is realized by a **moderate** one-parameter
gain rather than a delicate large-cancellation limit.

---

## 4. Direct comparison with the Stage-114 pure-mechanical threshold

At Stage 114 the pure-mechanical threshold was
\[
\Pi_* \approx 1.5088
\]
with no resolved mixed-lane correction.

The outlet-consistent coupled solve refines that picture:

- the shell lane must supply
  \[
  M_s^*\approx 1.80594111095636;
  \]
- but the mixed lane contributes back with
  \[
  M_q^*\mathcal S_q(\Pi_*)
  \approx -0.297111597463199151;
  \]
- leaving the exact net bias
  \[
  \Pi_*.
  \]

So the mixed D/N lane is not a tiny perturbation. It supplies an \(O(0.3)\)
correction to the shell demand while keeping the total bias in the moderate
canonical range.

---

## Result

Under the first outlet-consistent mouth closure, the full coupled mouth-layer
selection law has collapsed to
\[
\boxed{
\Pi
=
\Sigma_m\left[4-\mathcal S_q(\Pi)\right],
}
\]
and the canonical Family-1 branch is selected at
\[
\boxed{
\Sigma_m^*\approx 0.451485277739090.
}
\]

So the remaining mouth-layer ambiguity is now one dimensionless gain amplitude,
not an arbitrary source profile and not an arbitrary parent slope combination.

=== moving_throat_pde_stage119_coupled_mouth_status.md ===


# Moving-Throat PDE — Stage 119: Mouth-Layer Fixed-Point Status After the Coupled Solve

## What is now fixed

After Stages 116–118, the mouth-source selection problem has narrowed again.

1. The actual GNLS + localized-Maxwell mouth layer is no longer described by an
   effective slope alone. It is governed by the exact coupled fixed-point law
   \[
   \Pi
   =
   M_+\mathcal S(\Pi,\kappa_+)+M_-\mathcal S(\Pi,\kappa_-).
   \]

2. On the first explicit Family-1 branch with one static shell lane and one mixed
   D/N half-wave lane,
   \[
   \Pi = M_s + M_q\,\mathcal S_q(\Pi),
   \qquad
   \kappa_s=0,
   \qquad
   \kappa_q=\frac{\pi}{2}.
   \]

3. The canonical source-shape compensation point is now an exact gain line:
   \[
   M_s \approx 1.50882951349316 - 0.658075937605429\,M_q.
   \]

4. If the mouth-layer gains inherit the same \(4:-1\) shell-to-mixed weighting as
   the nontrivial compensated throat-core outlet, the selection law becomes
   \[
   \Pi
   =
   \Sigma_m\left[4-\mathcal S_q(\Pi)\right],
   \]
   with canonical Family-1 value
   \[
   \Sigma_m^* \approx 0.451485277739090.
   \]

## Meaning

The source-shape problem is no longer:

- which positive family to use,
- which sign branch is physical,
- or what effective combination
  \( \partial_z\delta V_{{\rm conf}}-q_*\partial_zA_0 \)
  should be guessed.

It is now simply:

\[
\boxed{
\text{what dimensionless mouth gain pair }(M_s,M_q)\text{ — or, under the
outlet-consistent closure, what single gain }\Sigma_m\text{ — does the real
moving-throat mouth layer select?}
}
\]

So the open microscopic bias problem has collapsed from a free profile question to
a small, explicit fixed-point law.

=== moving_throat_pde_stage120_core_to_mouth_gain_map.md ===

# Moving-Throat PDE — Stage 120: Explicit Core-to-Mouth Gain Map

## Goal

Derive the actual coupled mouth-layer gains `(M_s,M_q)` from an explicit
GNLS + localized-Maxwell throat-core ansatz instead of leaving them as abstract
fixed-point coefficients.

The key idea is to use the already derived concrete two-channel core response from
Stages 97–99 and embed it directly into the mouth electrochemical free energy.

---

## 1. Explicit electrochemical mouth free energy

Take the positive normalized mouth-source density `\sigma(z)` and let the brane-facing
mouth electrochemical potential be sourced by two scalar channels:

- a shell/compliance channel `U_s(z)`,
- a mixed localized-Maxwell channel `U_q(z)`.

Use the explicit free-energy density
\[
F_{\rm mouth}[\sigma,U_s,U_q]
=
\int_0^L dz\,
\Big[
\Theta_\sigma\,\sigma\big(\ln(\sigma/\sigma_*)-1\big)
+\sigma\big(\rho_c U_s(z)-\sigma_c U_q(z)\big)
\Big].
\]

Here the shell contribution enters with positive sign, while the mixed
localized-Maxwell contribution enters with the opposite sign, exactly as expected
for the parent electrochemical combination
\[
\partial_z\delta V_{\rm conf}-q_*\partial_z A_0.
\]

Varying with respect to `\sigma` gives the same exponential source law as before,
but with a self-consistent slope parameter
\[
\Pi
=
\frac{L}{\Theta_\sigma}
\Big[
\rho_c U_s'(0)-\sigma_c U_q'(0)
\Big].
\]

So the gains are now set by the *same* static and mixed coefficients that appear in
the explicit throat-core outlet.

---

## 2. Core coefficients from the exact Schur complement

Stage 97 already gave the concrete two-channel core outlet in the form
\[
\delta\Lambda_{\rm core}(z)
=
\rho_c-
\frac{\sigma_c}{1-\kappa_c z^2-i\gamma_c z^5},
\]
with exact coefficients
\[
\rho_c=\frac{g_s^2}{K_s},
\qquad
\sigma_c=
\frac{(K_sg_q-\lambda g_s)^2}{K_s(K_sK_q+\lambda^2)}.
\]

So the actual mouth-layer gains are
\[
\boxed{
M_s=
\frac{L}{\Theta_\sigma}\,\rho_c
=
\frac{L g_s^2}{K_s\Theta_\sigma},
}
\]
\[
\boxed{
M_q=
-\frac{L}{\Theta_\sigma}\,\sigma_c
=
-\frac{L (K_sg_q-\lambda g_s)^2}{K_s(K_sK_q+\lambda^2)\Theta_\sigma}.
}
\]

This is the first direct derivation of `(M_s,M_q)` from the explicit parent throat-core
ansatz.

The Stage-117 Family-1 fixed-point law therefore becomes
\[
\boxed{
\Pi=M_s+M_q\,\mathcal S_q(\Pi),
\qquad
\mathcal S_q(\Pi)=\mathcal S\!\left(\Pi,\frac{\pi}{2}\right).
}
\]

---

## 3. Immediate consequence

The abstract outlet-consistent closure of Stage 118 is no longer an extra guess.
It is the direct consequence of using the same Schur-complement core weights in the
mouth electrochemical source law.

The only remaining mouth-core data are now:

1. the shell stiffness `K_s`,
2. the mixed stiffness `K_q`,
3. the hybridization `\lambda`,
4. the mouth couplings `(g_s,g_q)`,
5. the source susceptibility `\Theta_\sigma`.

Once these are known on a branch, `(M_s,M_q)` are fixed.

---

## Result

The actual coupled mouth-layer gains are
\[
\boxed{
M_s=
\frac{L g_s^2}{K_s\Theta_\sigma},
\qquad
M_q=
-\frac{L (K_sg_q-\lambda g_s)^2}{K_s(K_sK_q+\lambda^2)\Theta_\sigma}.
}
\]

So the mouth fixed-point ambiguity has now collapsed from an abstract gain pair to
one explicit set of parent core quantities.

=== moving_throat_pde_stage121_normalized_mouth_gain_family.md ===

# Moving-Throat PDE — Stage 121: Normalized Mouth-Gain Family and Compensation Ratio

## Goal

Rewrite the explicit gain map of Stage 120 in the normalized parent variables already
used for the throat-core compensation family.

This turns the actual mouth gains into one overall amplitude times one exact ratio.

---

## 1. Normalized parent variables

Use the same core-normalized quantities as in Stages 102–106:
\[
\mathfrak r:=\frac{\lambda}{\sqrt{K_sK_q}},
\qquad
\mathfrak g_c:=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}},
\qquad
\Sigma_0:=\frac{L g_s^2}{K_s\Theta_\sigma}.
\]

Then Stage 120 gives immediately
\[
\boxed{M_s=\Sigma_0.}
\]

For the mixed gain,
\[
K_sg_q-\lambda g_s = g_s\sqrt{K_sK_q}(\mathfrak g_c-\mathfrak r),
\]
so
\[
\boxed{
M_q
=
-\Sigma_0\,
\frac{(\mathfrak g_c-\mathfrak r)^2}{1+\mathfrak r^2}.
}
\]

Define the exact mixed-to-shell gain ratio
\[
\boxed{
R_q:=-\frac{M_q}{M_s}=rac{(\mathfrak g_c-\mathfrak r)^2}{1+\mathfrak r^2}.
}
\]

So the full coupled Family-1 mouth law is
\[
\boxed{
\Pi = \Sigma_0\Big[1-R_q\,\mathcal S_q(\Pi)\Big].
}
\]

---

## 2. Exact compensation family

Stage 98 already fixed the core-balance family
\[
\mathfrak g_c
=
\mathfrak r\pm \frac12\sqrt{1+\mathfrak r^2}.
\]

Inserting this into the exact ratio gives
\[
\boxed{R_q=\frac14.}
\]

Therefore on the exact compensated branch,
\[
\boxed{
M_s=\Sigma_0,
\qquad
M_q=-\frac{\Sigma_0}{4}.
}
\]
If one defines
\[
\Sigma_m:=\frac{\Sigma_0}{4},
\]
then this is exactly the Stage-118 one-parameter closure
\[
\boxed{M_s=4\Sigma_m,\qquad M_q=-\Sigma_m.}
\]

So the Stage-118 closure is not independent bookkeeping. It is the normalized image of
Stage-120 plus the exact core-balance compensation surface.

---

## Result

The actual mouth gains are completely controlled by one overall amplitude `\Sigma_0`
and one exact ratio
\[
\boxed{R_q=\frac{(\mathfrak g_c-\mathfrak r)^2}{1+\mathfrak r^2}.}
\]

On the exact compensated branch this collapses to
\[
\boxed{R_q=\frac14,}
\]
so the outlet-consistent mouth closure is derived rather than assumed.

=== moving_throat_pde_stage122_family1_actual_mouth_gains.md ===

# Moving-Throat PDE — Stage 122: Actual Family-1 Mouth Gains

## Goal

Evaluate the explicit gain formulas of Stages 120–121 on the concrete Family-1 branch.

This answers the practical question:

> once the geometric core branch is fixed, what mouth gains does the actual throat-core
> ask for at the canonical compensation point `\Pi_*`?

---

## 1. Family-1 geometric input

From Stage 104,
\[
\mathfrak r_{F1}
=
\sqrt{\frac{12}{\pi^2}\left(\frac{37}{20}\right)^2-1}
\approx 1.77799353547498.
\]

From Stage 117,
\[
\Pi_* \approx 1.50882951349316,
\qquad
\mathcal S_q(\Pi_*)\approx 0.658075937605429.
\]

So the actual Family-1 fixed-point law is
\[
\Pi_* = M_s + M_q\,\mathcal S_q(\Pi_*).
\]
Using Stage 121,
\[
M_q=-R_q M_s,
\qquad
R_q=rac{(\mathfrak g_c-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2},
\]
so the exact shell gain required at `\Pi_*` is
\[
\boxed{
M_s^*(R_q)=rac{\Pi_*}{1-R_q\mathcal S_q(\Pi_*)},
\qquad
M_q^*(R_q)=-R_q M_s^*(R_q).
}
\]

---

## 2. Natural equal-normalized mouth-source branch

On the simplest natural branch,
\[
\mathfrak g_c=1.
\]
Then
\[
\boxed{
R_q^{\rm nat}
=
\frac{(1-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2}
\approx 0.145454452260421.
}
\]

So the actual Family-1 mouth gains selected at the canonical compensation point are
\[
\boxed{
M_s^{\rm nat,*}
=
\frac{\Pi_*}{1-R_q^{\rm nat}\mathcal S_q(\Pi_*)}
\approx 1.66854252965624,
}
\]
\[
\boxed{
M_q^{\rm nat,*}
=-R_q^{\rm nat}M_s^{\rm nat,*}
\approx -0.242696939724365.
}
\]

So the natural equal-normalized core source law does **not** exactly reproduce the
outlet-consistent canonical ratio, but it is already on the correct sign branch and not
far away from it.

---

## 3. Exact compensated branch

On the exact compensated family,
\[
\mathfrak g_c
=
\mathfrak r_{F1}-\frac12\sqrt{1+\mathfrak r_{F1}^2}
\approx 0.758035078944663,
\qquad
R_q=\frac14.
\]

Then the actual gains are
\[
\boxed{
M_s^{\rm comp,*}
=
\frac{\Pi_*}{1-\mathcal S_q(\Pi_*)/4}

\approx 1.80594111095636,
}
\]
\[
\boxed{
M_q^{\rm comp,*}=-\frac{M_s^{\rm comp,*}}{4}
\approx -0.451485277739090.
}
\]

This is exactly the Stage-118 one-parameter canonical branch.

---

## 4. Quantitative comparison

The natural and canonical gains differ only moderately:
\[
\frac{M_s^{\rm comp,*}}{M_s^{\rm nat,*}}-1
\approx 0.0823464663669,
\qquad
\frac{|M_q^{\rm comp,*}|}{|M_q^{\rm nat,*}|}
\approx 1.86097385480.
\]

So the shell gain changes by about `8.23%`, while the mixed gain magnitude must increase
by about a factor `1.86` to land exactly on the canonical compensated branch.

This is consistent with the earlier source-family result that the Family-1 branch was already
much closer to the lower compensated branch than to the forbidden upper one.

---

## Result

The actual Family-1 mouth gains are now explicit.

- natural equal-normalized core source:
  \[
  M_s\approx 1.66854,
  \qquad
  M_q\approx -0.24270;
  \]
- exact compensated canonical branch:
  \[
  M_s\approx 1.80594,
  \qquad
  M_q\approx -0.45149.
  \]

So the remaining ambiguity is no longer “what are the gains?” It is only whether the real
mouth core stays on the natural branch or shifts modestly toward the lower compensated one.

=== moving_throat_pde_stage123_selfmatched_mouth_susceptibility.md ===

# Moving-Throat PDE — Stage 123: Self-Matched Mouth Susceptibility Closure

## Goal

Push the explicit Family-1 gain formulas one step further by replacing the remaining
source susceptibility `\Theta_\sigma` with the first self-matched mouth-layer closure.

This does **not** solve the full PDE. It gives the first explicit parent formula for the
overall shell gain scale `\Sigma_0=M_s`.

---

## 1. Self-matched mouth susceptibility

Take the mouth source to live on the same active shell layer as the shell/compliance mode.
Then the entropic/source susceptibility is the Stage-43 type quantity
\[
\Theta_\sigma = h'(\rho_w) N_{\sigma\sigma}^{(m)}.
\]
On the self-matched mouth layer,
\[
N_{\sigma\sigma}^{(m)} = J_s = \frac{4\pi a^2\ell}{3},
\qquad
h'(\rho_w)=H_w=\frac{m_\psi c_{s,w}^2}{\rho_w}.
\]
So
\[
\boxed{
\Theta_\sigma = H_w J_s.
}
\]

Using the Stage-101 shell coupling
\[
g_s=\mathcal T_m J_s,
\]
and the healing-locked shell stiffness
\[
K_s=\frac{3\pi a^2\hbar^2}{5m_\psi\rho_w\ell},
\]
Stage 121 gives
\[
\Sigma_0=M_s=rac{L g_s^2}{K_s\Theta_\sigma}.
\]
A direct simplification yields
\[
\boxed{
\Sigma_0
=
\frac{20L\ell^2\rho_w^2\mathcal T_m^2}{9\hbar^2 c_{s,w}^2}.
}
\]

So the overall mouth shell gain is no longer abstract: it is fixed by one explicit
traction scale.

---

## 2. Convenient normalized traction parameter

Define the self-matched normalized mouth-traction amplitude
\[
\boxed{
\widehat T_m:=\frac{\rho_w\ell\sqrt{L}\,\mathcal T_m}{\hbar c_{s,w}}.
}
\]
Then
\[
\boxed{
\Sigma_0 = \frac{20}{9}\widehat T_m^2.
}
\]

So the full actual Family-1 gains are
\[
M_s = \frac{20}{9}\widehat T_m^2,
\qquad
M_q = -R_q\frac{20}{9}\widehat T_m^2.
\]

---

## 3. Required traction on the Family-1 branch

From Stage 122:

- natural equal-normalized branch requires
  \[
  M_s^{\rm nat,*}\approx 1.66854252965624,
  \qquad
  M_q^{\rm nat,*}\approx -0.242696939724365;
  \]
- exact compensated branch requires
  \[
  M_s^{\rm comp,*}\approx 1.80594111095636,
  \qquad
  M_q^{\rm comp,*}\approx -0.451485277739090.
  \]

Therefore the required normalized traction amplitudes are
\[
\boxed{
\widehat T_m^{\rm nat,*}
=\sqrt{\frac{9M_s^{\rm nat,*}}{20}}
\approx 0.866512630228382,
}
\]
\[
\boxed{
\widehat T_m^{\rm comp,*}
=\sqrt{\frac{9M_s^{\rm comp,*}}{20}}
\approx 0.901484054174206.
}
\]

So the exact outlet-consistent compensated branch requires only
\[
\boxed{
\frac{\widehat T_m^{\rm comp,*}}{\widehat T_m^{\rm nat,*}}-1
\approx 0.0403588161624,
}
\]
that is, only about a `4.04%` enhancement in the normalized mouth traction relative
to the natural equal-normalized branch.

---

## Result

Under the first self-matched mouth susceptibility closure, the overall shell gain is
fully explicit:
\[
\boxed{
M_s=\Sigma_0=\frac{20}{9}\widehat T_m^2.
}
\]

On the explicit Family-1 branch, the natural and exact-compensated mouth closures differ
by only about `4%` in the normalized traction amplitude.

That is the cleanest explicit parent-level narrowing of the mouth-gain problem so far.

=== moving_throat_pde_stage124_mouth_gain_status.md ===

# Moving-Throat PDE — Stage 124: Mouth-Gain Status Update

The coupled mouth-layer problem is now much tighter than it was at Stage 116.

## What is now explicit

1. The actual mouth gains are derived from the explicit throat-core ansatz:
   \[
   M_s=\frac{L g_s^2}{K_s\Theta_\sigma},
   \qquad
   M_q=-\frac{L (K_sg_q-\lambda g_s)^2}{K_s(K_sK_q+\lambda^2)\Theta_\sigma}.
   \]

2. In normalized core variables,
   \[
   M_q=-M_s\frac{(\mathfrak g_c-\mathfrak r)^2}{1+\mathfrak r^2}.
   \]

3. On the exact core-balance compensation family,
   \[
   M_q=-\frac{M_s}{4},
   \qquad
   M_s=4\Sigma_m,
   \qquad
   M_q=-\Sigma_m.
   \]

4. On the explicit Family-1 branch,
   the natural equal-normalized source law gives
   \[
   M_s\approx 1.66854,
   \qquad
   M_q\approx -0.24270,
   \]
   while the exact outlet-consistent canonical branch gives
   \[
   M_s\approx 1.80594,
   \qquad
   M_q\approx -0.45149.
   \]

5. Under the self-matched mouth susceptibility closure,
   \[
   M_s=\frac{20}{9}\widehat T_m^2,
   \]
   so the natural and canonical branches differ by only about `4%` in the normalized
   traction amplitude.

## What remains open

The mouth-layer ambiguity is no longer a profile ambiguity and no longer a free gain pair.
It has shrunk to two very concrete microscopic questions:

1. does the real moving-throat mouth source choose the natural equal-normalized branch
   \(\mathfrak g_c\approx1\) or the nearby lower compensated branch,
2. and does the actual mouth traction land at the corresponding `\widehat T_m` value?

That is a much smaller target than the original abstract `\Pi_m` problem.

=== moving_throat_pde_stage125_selfconsistent_mouth_branch.md ===

# Moving-Throat PDE — Stage 125: Self-Consistent Mouth-Branch Law

## Goal

Close the loop between:

1. the explicit positive mouth-layer source family
   \[
   \mathfrak g_\Pi
   =
   \frac{2\Pi(2\Pi e^\Pi+\pi)}
   {(4\Pi^2+\pi^2)(e^\Pi-1)},
   \]
2. the explicit core-to-mouth gain map
   \[
   M_q=-M_s\,\frac{(\mathfrak g_c-\mathfrak r)^2}{1+\mathfrak r^2},
   \]
3. and the coupled Family-1 mouth fixed-point law
   \[
   \Pi=M_s+M_q\,\mathcal S_q(\Pi).
   \]

This removes the last free branch label \(\mathfrak g_c\) from the mouth problem.

---

## 1. Self-consistent Family-1 gain law

On the explicit positive mouth-layer branch we must identify
\[
\mathfrak g_c=\mathfrak g_\Pi.
\]
So the Family-1 mixed-to-shell gain ratio becomes
\[
\boxed{
R_q(\Pi):=-\frac{M_q}{M_s}
=
\frac{(\mathfrak g_\Pi-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2}.
}
\]

The coupled Family-1 mouth equation then closes to
\[
\boxed{
\Pi
=
\Sigma_0\Big[1-R_q(\Pi)\,\mathcal S_q(\Pi)\Big],
}
\]
where
\[
\Sigma_0:=M_s,
\qquad
\mathcal S_q(\Pi)=\mathcal S\!\left(\Pi,\frac{\pi}{2}\right).
\]

Equivalently, the shell gain required to realize a given bias \(\Pi\) is
\[
\boxed{
\Sigma_0(\Pi)
=
\frac{\Pi}{1-R_q(\Pi)\mathcal S_q(\Pi)}.
}
\]

So the explicit mouth branch is now one scalar function \(\Sigma_0(\Pi)\), not a free gain pair.

---

## 2. Self-matched mouth susceptibility closure

Under the self-matched closure from Stage 123,
\[
\Sigma_0=\frac{20}{9}\,\widehat T_m^2.
\]
Therefore the normalized mouth-traction branch is
\[
\boxed{
\widehat T_m(\Pi)
=
\sqrt{\frac{9\Pi}{20\left[1-R_q(\Pi)\mathcal S_q(\Pi)\right]}}.
}
\]

This formula is the first explicit self-consistent map from the parent mouth bias \(\Pi\)
to the physical normalized traction \(\widehat T_m\).

---

## 3. Canonical Family-1 point

The already-frozen lower compensated Family-1 mouth point is
\[
\mathfrak g_-^{F1}
=
\mathfrak r_{F1}-\frac12\sqrt{1+\mathfrak r_{F1}^2}
\approx 0.758035078944663,
\]
and the unique mouth-layer bias solving
\[
\mathfrak g_\Pi=\mathfrak g_-^{F1}
\]
is
\[
\boxed{
\Pi_*\approx 1.50882951349316.
}
\]

At that point
\[
R_q(\Pi_*)=\frac14,
\qquad
\mathcal S_q(\Pi_*)\approx 0.658075937605429,
\]
so
\[
\boxed{
\Sigma_0(\Pi_*)\approx 1.80594111095636,
\qquad
\widehat T_m(\Pi_*)\approx 0.901484054174205.
}
\]

This reproduces the earlier outlet-consistent canonical branch exactly, but now from the
fully self-consistent positive mouth-layer + core map.

---

## 4. Meaning

The mouth problem is now much sharper:

- the source shape is no longer free,
- the gain ratio is no longer free,
- and the shell traction is no longer free.

Once the positive mouth layer is assumed, the branch is governed by one scalar variable
\(\Pi\), with a unique canonical compensation point at moderate finite bias and moderate
finite traction.

=== moving_throat_pde_stage126_equal_normalized_singular_limit.md ===

# Moving-Throat PDE — Stage 126: Equal-Normalized Branch Is a Singular Limit

## Goal

Determine whether the simple equal-normalized mouth-source branch
\[
\mathfrak g_c=1
\]
can occur at any **finite** positive mouth-layer bias \(\Pi\), and whether it corresponds
to a regular finite-traction mouth state.

---

## 1. Exact strict inequality for the exponential mouth layer

For the explicit positive mouth-layer family,
\[
\mathfrak g_\Pi
=
\frac{2\Pi(2\Pi e^\Pi+\pi)}
{(4\Pi^2+\pi^2)(e^\Pi-1)}.
\]
A direct rearrangement gives
\[
1-\mathfrak g_\Pi
=
\frac{\pi^2(e^\Pi-1)-2\pi\Pi-4\Pi^2}
{(4\Pi^2+\pi^2)(e^\Pi-1)}.
\]

The numerator splits exactly as
\[
\pi^2\!\left(e^\Pi-1-\Pi-\frac{\Pi^2}{2}\right)
+\Pi(\pi^2-2\pi)
+\Pi^2\!\left(\frac{\pi^2}{2}-4\right).
\]
Every term is positive for \(\Pi>0\), since
\[
e^\Pi-1-\Pi-\frac{\Pi^2}{2}>0,
\qquad
\pi^2-2\pi>0,
\qquad
\frac{\pi^2}{2}-4>0.
\]

Therefore
\[
\boxed{
0<\mathfrak g_\Pi<1
\qquad \text{for every finite } \Pi>0.
}
\]

So the equal-normalized branch \(\mathfrak g_c=1\) does **not** occur at finite positive bias.

---

## 2. Singular point-source limit

The same family obeys
\[
\lim_{\Pi\to\infty}\mathfrak g_\Pi=1.
\]
So the equal-normalized branch is recovered only in the singular point-source limit
\[
\boxed{
\Pi\to+\infty.
}
\]

In source-profile language, this is the delta-function mouth limit.

---

## 3. Traction divergence on the equal-normalized branch

From Stage 125,
\[
\widehat T_m(\Pi)
=
\sqrt{\frac{9\Pi}{20\left[1-R_q(\Pi)\mathcal S_q(\Pi)\right]}},
\qquad
R_q(\Pi)=\frac{(\mathfrak g_\Pi-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2}.
\]

As \(\Pi\to\infty\),
\[
\mathfrak g_\Pi\to 1,
\qquad
\mathcal S_q(\Pi)\to 1,
\qquad
R_q(\Pi)\to R_\infty
=
\frac{(1-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2}
\approx 0.145454452260420.
\]
Hence
\[
\Sigma_0(\Pi)
=
\frac{\Pi}{1-R_q(\Pi)\mathcal S_q(\Pi)}
\sim
\frac{\Pi}{1-R_\infty},
\]
and therefore
\[
\boxed{
\widehat T_m(\Pi)\sim
\sqrt{\frac{9}{20(1-R_\infty)}}\,\Pi^{1/2}
\approx 0.725669130700713\,\Pi^{1/2}.
}
\]

So the equal-normalized branch requires
\[
\boxed{
\widehat T_m\to\infty.
}
\]

---

## Result

The explicit positive exponential mouth-layer family proves that:

1. \(\mathfrak g_c=1\) is not a finite-bias branch,
2. it is reached only as \(\Pi\to\infty\),
3. and in that same limit the normalized mouth traction diverges.

So the naive equal-normalized mouth-source branch is **not** a regular finite branch of the
explicit mouth-layer dynamics. It is a singular point-source limit.

=== moving_throat_pde_stage127_unique_regular_canonical_branch.md ===

# Moving-Throat PDE — Stage 127: Unique Regular Canonical Mouth Branch

## Goal

Combine the positive-source theorem, the explicit mouth-layer bias family, and the
self-consistent core-to-mouth gain law to decide which compensated branch survives as a
**regular finite-bias / finite-traction** Family-1 solution.

---

## 1. Upper branch remains impossible

The explicit Family-1 compensated branches are
\[
\mathfrak g_\pm^{F1}
=
\mathfrak r_{F1}\pm\frac12\sqrt{1+\mathfrak r_{F1}^2},
\]
with
\[
\mathfrak g_-^{F1}\approx 0.758035078944663,
\qquad
\mathfrak g_+^{F1}\approx 2.79795199200529.
\]

But every positive normalized source profile on the first D/N interval satisfies
\[
0\le \mathfrak g_c\le 1.
\]
So
\[
\boxed{
\mathfrak g_+^{F1}>1
\quad\Rightarrow\quad
\text{the upper compensated branch is impossible.}
}
\]

---

## 2. Lower branch is uniquely reachable at finite bias

From Stage 113, \(\mathfrak g_\Pi\) is strictly monotone increasing in \(\Pi\), with range
\[
\frac{2}{\pi}<\mathfrak g_\Pi<1
\qquad (\Pi>0).
\]
Since
\[
\frac{2}{\pi}\approx 0.636619772367581
<
\mathfrak g_-^{F1}\approx 0.758035078944663
<
1,
\]
there exists a unique positive finite \(\Pi_*\) such that
\[
\mathfrak g_{\Pi_*}=\mathfrak g_-^{F1}.
\]
Numerically,
\[
\boxed{
\Pi_*\approx 1.50882951349316.
}
\]

So the lower compensated branch is not merely allowed. It is the **unique finite-bias**
compensated branch.

---

## 3. Finite regular traction at the canonical point

At \(\Pi=\Pi_*\), Stage 125 gives
\[
\Sigma_0(\Pi_*)\approx 1.80594111095636,
\qquad
\widehat T_m(\Pi_*)\approx 0.901484054174205.
\]
These are finite and moderate.

For comparison, the self-matched derivative-profile point \(\mathfrak g=\pi/4\) occurs at
\[
\Pi_{\rm match}\approx 1.90848600654854,
\qquad
\widehat T_m(\Pi_{\rm match})\approx 1.01132972803599.
\]
So the canonical compensated branch is reached **before** the self-matched derivative point
and far before the singular point-source limit.

---

## 4. Final branch-selection theorem

Within the explicit Family-1 core + positive exponential mouth-layer closure:

- the upper compensated branch is impossible,
- the equal-normalized branch is a singular limit,
- and the lower compensated branch is the unique regular finite-bias / finite-traction
  branch that preserves the canonical outgoing quadrupole fingerprint.

Equivalently,
\[
\boxed{
\text{the canonical Family-1 mouth branch is uniquely selected as the regular positive-source branch.}
}
\]

---

## Meaning

This is a real narrowing.

The remaining mouth ambiguity is no longer a branch-choice ambiguity at all.
On the explicit Family-1 positive mouth-layer closure, the canonical outgoing-preserving
branch is the only regular finite one.

=== moving_throat_pde_stage128_mouth_branch_selection_status.md ===

# Moving-Throat PDE — Stage 128: Mouth-Branch Selection Status

The mouth-source side is now much tighter than it was at Stage 124.

## What is now explicit

1. The positive mouth-layer family selects a unique bias factor
   \[
   \mathfrak g_\Pi
   =
   \frac{2\Pi(2\Pi e^\Pi+\pi)}
   {(4\Pi^2+\pi^2)(e^\Pi-1)},
   \]
   with range
   \[
   \frac{2}{\pi}<\mathfrak g_\Pi<1
   \qquad (\Pi>0).
   \]

2. The self-consistent Family-1 gain law is
   \[
   \Pi=\Sigma_0\bigl[1-R_q(\Pi)\mathcal S_q(\Pi)\bigr],
   \qquad
   R_q(\Pi)=\frac{(\mathfrak g_\Pi-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2},
   \]
   and under the self-matched closure
   \[
   \Sigma_0=\frac{20}{9}\widehat T_m^2.
   \]

3. The upper compensated branch is impossible because
   \[
   \mathfrak g_+^{F1}>1.
   \]

4. The naive equal-normalized branch \(\mathfrak g_c=1\) is not a finite branch:
   it is reached only as
   \[
   \Pi\to\infty,
   \qquad
   \widehat T_m\to\infty.
   \]

5. The lower compensated branch is uniquely reachable at finite values,
   \[
   \Pi_*\approx 1.50882951349316,
   \qquad
   \widehat T_m(\Pi_*)\approx 0.901484054174205.
   \]

## What remains open

The remaining PDE-facing ambiguity is now much smaller.

It is no longer:
- which compensation branch,
- or whether equal-normalized sourcing survives.

It is now only this:

1. how accurately the actual moving-throat mouth layer follows the explicit exponential
   electrochemical family,
2. and how large the finite corrections are to the unique regular canonical point
   \[
   (\Pi_*,\widehat T_{m,*}).
   \]

That is a much smaller target than the earlier abstract branch ambiguity.

=== moving_throat_pde_stage129_positive_deformation_expansion.md ===


# Moving-Throat PDE — Stage 129: Finite-Correction Expansion for Positive Mouth-Layer Deformations

## Goal

Test how rigid the unique regular Family-1 canonical mouth branch is under **non-exponential but still positive localized** mouth-source deformations.

The canonical exponential branch already fixed
\[
\Pi_* \approx 1.50882951349316,
\qquad
\widehat T_{m,*}\approx 0.901484054174205,
\qquad
\mathfrak g_*=\mathfrak g_-^{F1}\approx 0.758035078944663,
\qquad
\mathcal S_*=\mathcal S_q(\Pi_*)\approx 0.658075937605429.
\]

The question now is: if the actual mouth layer is **not exactly exponential**, how much must the canonical point move?

---

## 1. Positive normalized deformation family

Work on the dimensionless mouth interval
\[
x=\frac{z}{L}\in[0,1].
\]
Let the canonical exponential source be
\[
\Sigma_*(x)=\frac{\Pi_* e^{-\Pi_* x}}{1-e^{-\Pi_*}}.
\]

Take any other positive normalized mouth profile
\[
\varsigma(x)\ge 0,
\qquad
\int_0^1 \varsigma(x)\,dx=1,
\]
and define the exact convex positive deformation family
\[
\boxed{
\Sigma_\epsilon(x)=(1-\epsilon)\Sigma_*(x)+\epsilon\,\varsigma(x),
\qquad
0\le \epsilon\ll 1.
}
\]

This preserves positivity and normalization exactly.

---

## 2. The only two source moments that matter at first order

Define the Family-1 overlap and mixed-channel kernels
\[
c(x)=\cos\!\left(\frac{\pi x}{2}\right),
\qquad
K_q(x)=\frac{\cosh\!\left(\frac{\pi}{2}(1-x)\right)}{\cosh(\pi/2)}.
\]

Then the source shape enters the mouth-core system only through the two averages
\[
\bar g[\sigma]=\int_0^1 \sigma(x)\,c(x)\,dx,
\qquad
\bar S[\sigma]=\int_0^1 \sigma(x)\,K_q(x)\,dx.
\]

For the positive convex family this is exact:
\[
\boxed{
\bar g_\epsilon
=
\mathfrak g_*+\epsilon(\bar g_\varsigma-\mathfrak g_*),
\qquad
\bar S_\epsilon
=
\mathcal S_*+\epsilon(\bar S_\varsigma-\mathcal S_*),
}
\]
where
\[
\bar g_\varsigma=\int_0^1\varsigma(x)c(x)\,dx,
\qquad
\bar S_\varsigma=\int_0^1\varsigma(x)K_q(x)\,dx.
\]

So the first-order mouth-shape problem is **two-dimensional**:
only \((\bar g_\varsigma,\bar S_\varsigma)\) matter.

---

## 3. Retuning the electrochemical bias to stay on the canonical lower branch

To remain on the canonical outgoing-preserving lower branch, the overlap must stay fixed at
\[
\bar g = \mathfrak g_*.
\]
So if the source is deformed, the exponential bias must shift
\[
\Pi=\Pi_*+\delta\Pi.
\]

Let
\[
\mathfrak g'_*=
\left.\frac{d\mathfrak g_\Pi}{d\Pi}\right|_{\Pi_*},
\qquad
\mathcal S'_*
=
\left.\frac{d\mathcal S_q}{d\Pi}\right|_{\Pi_*}.
\]
Then the exact first-order compensation law is
\[
\boxed{
\delta\Pi
=
-\epsilon\,\frac{\bar g_\varsigma-\mathfrak g_*}{\mathfrak g'_*}.
}
\]

Numerically,
\[
\boxed{
\mathfrak g'_* \approx 0.0714453558083195,
\qquad
\mathcal S'_* \approx 0.0483709542125041.
}
\]

So the compensated mixed-channel response changes by
\[
\boxed{
\delta \mathcal S_q
=
\epsilon\left[
(\bar S_\varsigma-\mathcal S_*)
-
\frac{\mathcal S'_*}{\mathfrak g'_*}\,(\bar g_\varsigma-\mathfrak g_*)
\right].
}
\]

This is the first exact non-exponential correction formula.

=== moving_throat_pde_stage130_first_order_rigidity_kernel.md ===


# Moving-Throat PDE — Stage 130: First-Order Rigidity Kernel at the Canonical Family-1 Point

## Goal

Turn the general first-order correction formulas from Stage 129 into one explicit
rigidity law for the physical normalized mouth traction
\[
\widehat T_m.
\]

---

## 1. Canonical branch traction law

On the lower compensated branch one has
\[
R_q=\frac14,
\qquad
\Sigma_0=\frac{\Pi}{1-\mathcal S_q/4},
\qquad
\widehat T_m=\sqrt{\frac{9\Sigma_0}{20}}.
\]

So at fixed canonical branch structure,
\[
\delta \Sigma_0
=
\frac{\delta\Pi}{1-\mathcal S_*/4}
+
\frac{\Pi_*}{4(1-\mathcal S_*/4)^2}\,\delta\mathcal S_q.
\]

Substituting the Stage-129 retuning formulas gives
\[
\boxed{
\delta \widehat T_m
=
\epsilon\Big[
A_T\,(\bar g_\varsigma-\mathfrak g_*)
+
B_T\,(\bar S_\varsigma-\mathcal S_*)
\Big],
}
\]
with
\[
A_T
=
-\frac{9}{40\widehat T_{m,*}}
\left[
\frac{1}{\mathfrak g'_*(1-\mathcal S_*/4)}
+
\frac{\Pi_*\,\mathcal S'_*}{4\mathfrak g'_*(1-\mathcal S_*/4)^2}
\right],
\]
\[
B_T
=
\frac{9}{40\widehat T_{m,*}}
\frac{\Pi_*}{4(1-\mathcal S_*/4)^2}.
\]

Numerically,
\[
\boxed{
A_T \approx -4.27263956256927,
\qquad
B_T \approx 0.134875005736706.
}
\]

So the canonical branch is vastly more sensitive to overlap changes than to
mixed-kernel changes:
\[
\boxed{
\frac{|A_T|}{B_T}\approx 31.6785.
}
\]

---

## 2. One effective rigidity kernel

Because \(\Sigma_\epsilon-\Sigma_*\) integrates to zero, the traction shift can be written
as a single weighted overlap:
\[
\boxed{
\delta \widehat T_m
=
\epsilon
\int_0^1
\mathcal W_*(x)\,
\bigl[\varsigma(x)-\Sigma_*(x)\bigr]\,dx,
}
\]
with centered weight
\[
\boxed{
\mathcal W_*(x)
=
A_T\Big(c(x)-\mathfrak g_*\Big)
+
B_T\Big(K_q(x)-\mathcal S_*\Big).
}
\]

So after retuning the electrochemical bias, the entire first-order non-exponential
shape sensitivity collapses to **one scalar kernel** \(\mathcal W_*(x)\).

This is the right rigidity statement:

> once the canonical lower branch is fixed, a positive mouth-layer deformation can move the
> Family-1 point only through its overlap with one known rigidity kernel.

---

## 3. Meaning

The branch is not infinitely rigid. But it is much more rigid than the earlier
branch-choice ambiguity suggested:

- the deformation space collapses to two source moments,
- the canonical retuning removes one of them as an independent datum,
- and the remaining first-order traction shift is controlled by one explicit kernel.

That is a strong reduction of the mouth-side uncertainty.

=== moving_throat_pde_stage131_representative_positive_families.md ===


# Moving-Throat PDE — Stage 131: Representative Non-Exponential Positive Mouth Families

## Goal

Evaluate the first-order Family-1 correction formulas on two explicit positive
non-exponential mouth families and on their convex interpolation.

This converts the abstract rigidity kernel of Stage 130 into concrete scale estimates.

---

## 1. Uniform broadening family

Take the broadest positive normalized profile on the mouth interval,
\[
\varsigma_u(x)=1.
\]

Then
\[
\bar g_u=\int_0^1 \cos\!\left(\frac{\pi x}{2}\right)\,dx=\frac{2}{\pi},
\]
and
\[
\bar S_u
=
\int_0^1 K_q(x)\,dx
=
\frac{2\tanh(\pi/2)}{\pi}.
\]

Numerically,
\[
\bar g_u\approx 0.636619772367581,
\qquad
\bar S_u\approx 0.583877311158896.
\]

The canonical retuning shifts are therefore
\[
\boxed{
\frac{\delta\Pi_u}{\epsilon}\approx +1.69941496131430,
\qquad
\frac{\delta\widehat T_{m,u}}{\epsilon}\approx +0.508756302215085.
}
\]

So broadening the source toward uniformity forces the canonical branch to move to
**larger** bias and **larger** traction.

---

## 2. Self-matched derivative family

Take the positive self-matched derivative profile,
\[
\varsigma_d(x)=\frac{\pi}{2}\cos\!\left(\frac{\pi x}{2}\right).
\]

Then
\[
\bar g_d=\frac{\pi}{4},
\]
and
\[
\bar S_d
=
\frac{1+\sinh(\pi/2)}{2\cosh(\pi/2)}
\approx 0.657844575502830.
\]

The canonical retuning shifts are
\[
\boxed{
\frac{\delta\Pi_d}{\epsilon}\approx -0.382993186095921,
\qquad
\frac{\delta\widehat T_{m,d}}{\epsilon}\approx -0.116943802151809.
}
\]

So sharpening the source toward the self-matched derivative profile moves the
canonical branch to **smaller** bias and **smaller** traction.

---

## 3. Convex interpolation between the two positive non-exponential families

Now interpolate the two explicit positive families:
\[
\varsigma_\lambda(x)=(1-\lambda)\varsigma_u(x)+\lambda \varsigma_d(x),
\qquad 0\le \lambda\le 1.
\]

Because the first-order formulas are affine in \((\bar g_\varsigma,\bar S_\varsigma)\),
the canonical shifts are exactly
\[
\boxed{
\frac{\delta\Pi_\lambda}{\epsilon}
=
1.69941496131430
-
2.08240814741023\,\lambda,
}
\]
\[
\boxed{
\frac{\delta\widehat T_{m,\lambda}}{\epsilon}
=
0.508756302215085
-
0.625700104366894\,\lambda.
}
\]

The **bias-neutral** interpolation point is
\[
\boxed{
\lambda_{\Pi,0}\approx 0.816081594488463.
}
\]
Equivalently,
\[
1-\lambda_{\Pi,0}\approx 0.183918405511537,
\]
which exactly matches the earlier Stage-110 broadening fraction
\[
\xi_* \approx 0.183918405511538
\]
for the positive family that hit the Family-1 lower compensated branch.

The **traction-neutral** interpolation point is
\[
\boxed{
\lambda_{T,0}\approx 0.813099276577336.
}
\]

So the first-order finite-correction theory is perfectly consistent with the
earlier exact positive-family compensation result.

---

## 4. Meaning

For the two most natural explicit positive non-exponential mouth families:

- broadening toward uniformity changes the canonical point upward,
- sharpening toward the self-matched derivative profile changes it downward,
- and the zero-shift mixture lies very close to the earlier exact compensation
  broadening fraction.

This is the best evidence so far that the canonical Family-1 mouth branch is
**rigid but not brittle**: moderate positive non-exponential corrections move it in a
controlled, almost one-parameter way.

=== moving_throat_pde_stage132_mouth_rigidity_status.md ===


# Moving-Throat PDE — Stage 132: Non-Exponential Mouth-Rigidity Status

The mouth-side branch problem is now significantly tighter than it was at Stage 128.

## What is now explicit

1. Any positive normalized mouth deformation near the canonical exponential branch can be
   represented as
   \[
   \Sigma_\epsilon=(1-\epsilon)\Sigma_*+\epsilon\varsigma,
   \qquad
   \varsigma\ge0,
   \qquad
   \int_0^1\varsigma=1.
   \]

2. At first order, the entire deformation enters only through the two source moments
   \[
   \bar g_\varsigma=\int_0^1\varsigma\,c,
   \qquad
   \bar S_\varsigma=\int_0^1\varsigma\,K_q.
   \]

3. Retuning the electrochemical bias to stay on the canonical lower branch gives
   \[
   \delta\Pi
   =
   -\epsilon\frac{\bar g_\varsigma-\mathfrak g_*}{\mathfrak g'_*},
   \]
   and the physical traction shift reduces to
   \[
   \delta \widehat T_m
   =
   \epsilon\Big[
   A_T(\bar g_\varsigma-\mathfrak g_*)
   +
   B_T(\bar S_\varsigma-\mathcal S_*)
   \Big],
   \]
   with
   \[
   A_T\approx -4.27263956256927,
   \qquad
   B_T\approx 0.134875005736706.
   \]

4. So the canonical branch is controlled primarily by the overlap channel, not by the
   mixed-kernel channel:
   \[
   |A_T|/B_T\approx 31.68.
   \]

5. For representative positive non-exponential families:
   - uniform broadening raises the canonical point,
   - derivative matching lowers it,
   - and the zero-shift mixture coincides almost exactly with the earlier Stage-110
     exact compensation broadening fraction.

## Updated interpretation

Inside the explicit Family-1 positive mouth-layer closure, the mouth-side ambiguity is now
no longer a branch-selection problem and no longer a large shape-space uncertainty.

It has been reduced to:

1. a single canonical branch,
2. one explicit rigidity kernel \(\mathcal W_*(x)\),
3. and modest finite shifts under smooth positive non-exponential deformations.

So the remaining PDE-facing question is not “which mouth branch?” but rather:

\[
\boxed{
\text{what small non-exponential correction does the real moving-throat mouth layer induce around }(\Pi_*,\widehat T_{m,*})?
}
\]

That is a much smaller target than the earlier mouth-source ambiguity.

=== moving_throat_pde_stage133_full_profile_residual.md ===

# Moving-Throat PDE — Stage 133: Exact Full-Profile Mouth Potential and Curvature Residual

## Goal

Derive the **actual non-exponential correction** selected by the full coupled
GNLS + localized-Maxwell Family-1 mouth layer, instead of keeping only the linear
electrochemical slope that produced the exponential source law.

The canonical exponential branch came from replacing the full mouth potential by
its tangent
\[
\Phi_{\rm lin}(x)=\Pi_* x.
\]
The next honest step is to keep the full channel profiles.

---

## 1. Exact Family-1 channel profiles

On the Family-1 branch the two mouth channels are

- a static shell lane \(\kappa_s=0\),
- the first mixed D/N half-wave \(\kappa_q=\pi/2\).

For the exponential source branch
\[
\Sigma_*(x)=\frac{\Pi_* e^{-\Pi_* x}}{1-e^{-\Pi_*}},
\qquad x\in[0,1],
\]
the exact shell profile is
\[
\boxed{
T_s(x;\Pi_*)
=
\frac{1-e^{-\Pi_* x}}{\Pi_*(1-e^{-\Pi_*})}
-
\frac{x e^{-\Pi_*}}{1-e^{-\Pi_*}},
}
\]
while the exact mixed D/N profile is
\[
\boxed{
T_q(x;\Pi_*)
=
A_q\sinh\!\left(\frac{\pi x}{2}\right)
-
C_q\cosh\!\left(\frac{\pi x}{2}\right)
+
C_q e^{-\Pi_* x},
}
\]
with
\[
C_q=
\frac{\Pi_*}{(1-e^{-\Pi_*})(\pi^2/4-\Pi_*^2)},
\]
\[
A_q=
\frac{C_q\left(\frac{\pi}{2}\sinh(\pi/2)+\Pi_*e^{-\Pi_*}\right)}
{\frac{\pi}{2}\cosh(\pi/2)}.
\]

The compensated core outlet from Stage 118 fixes
\[
M_s^*=4\Sigma_m^*,
\qquad
M_q^*=-\Sigma_m^*,
\qquad
\Sigma_m^*\approx 0.451485277739090,
\qquad
\Pi_*=\Sigma_m^*\left[4-\mathcal S_q(\Pi_*)\right].
\]

So the exact full mouth potential generated by the canonical source is
\[
\boxed{
\Phi_*(x)=4\Sigma_m^*\,T_s(x;\Pi_*)-\Sigma_m^*\,T_q(x;\Pi_*).
}
\]

---

## 2. Curvature residual relative to the exponential tangent

Define the full-profile residual by subtracting the tangent exponential potential:
\[
\boxed{
R_*(x):=\Phi_*(x)-\Pi_* x.
}
\]

Using the fixed-point relation \(\Pi_*=\Sigma_m^*(4-\mathcal S_q(\Pi_*))\), this can
be rewritten exactly as
\[
\boxed{
R_*(x)=
\Sigma_m^*\Big[
4T_s(x;\Pi_*)-T_q(x;\Pi_*)-(4-\mathcal S_q(\Pi_*))x
\Big].
}
\]

The first exact mouth-shape theorem is

\[
\boxed{
R_*(0)=0,
\qquad
R_*'(0)=0.
}
\]

So the exponential source law is tangent to the full coupled mouth potential at the
mouth.

The second exact theorem is the curvature law
\[
\boxed{
R_*''(0)
=
-\,3\Sigma_m^*\frac{\Pi_*}{1-e^{-\Pi_*}}<0.
}
\]

So the full coupled Family-1 mouth potential is **sublinear** relative to the
exponential tangent right at the mouth. That means the actual self-consistent source
is broader than the tangent exponential branch.

---

## 3. Meaning

This is the first explicit non-exponential correction selected by the full coupled
mouth layer.

The important point is not just that a correction exists. The important point is:

- it is generated automatically by the exact shell + mixed D/N channel profiles,
- it is tangent-matched at the mouth,
- and its leading nontrivial effect is a negative curvature residual.

So the remaining mouth uncertainty is no longer a branch-choice issue.
It is a finite profile correction around the already-selected lower compensated
branch.
=== moving_throat_pde_stage134_first_order_selected_correction.md ===

# Moving-Throat PDE — Stage 134: First-Order Source Correction Selected by the Full Mouth Profile

## Goal

Project the exact full-profile mouth residual \(R_*(x)\) onto the previously derived
Family-1 rigidity formulas.

This produces the **actual** first-order non-exponential correction selected by the
full coupled GNLS + localized-Maxwell mouth layer.

---

## 1. First-order self-consistent source law

The exact self-consistent positive mouth source is
\[
\Sigma_{\rm full}(x)\propto e^{-\Phi_*(x)}
=
e^{-\Pi_* x - R_*(x)}.
\]

Expanding about the canonical exponential source
\[
\Sigma_*(x)=\frac{\Pi_*e^{-\Pi_*x}}{1-e^{-\Pi_*}}
\]
gives the normalized first-order correction
\[
\boxed{
\Sigma_{\rm act}(x)
=
\Sigma_*(x)\Big[1-\widetilde R_*(x)\Big]
+O(R_*^2),
\qquad
\widetilde R_*(x):=R_*(x)-\langle R_*\rangle_*,
}
\]
where
\[
\langle f\rangle_*:=\int_0^1 \Sigma_*(x)f(x)\,dx.
\]

So the actual selected deformation is not free:
it is exactly the centered residual of the full mouth potential.

---

## 2. Only two moment shifts matter

As in Stages 129–130, define
\[
c(x)=\cos\!\left(\frac{\pi x}{2}\right),
\qquad
K_q(x)=\frac{\cosh\!\left(\frac{\pi}{2}(1-x)\right)}{\cosh(\pi/2)}.
\]

Then the actual first-order shifts are
\[
\boxed{
\delta \mathfrak g_{\rm act}
=
\int_0^1 c(x)\,\delta\Sigma_*(x)\,dx
=
-\operatorname{Cov}_*(c,R_*),
}
\]
\[
\boxed{
\delta \mathcal S_{\rm act}
=
\int_0^1 K_q(x)\,\delta\Sigma_*(x)\,dx
=
-\operatorname{Cov}_*(K_q,R_*),
}
\]
with
\[
\operatorname{Cov}_*(f,h)
=
\langle fh\rangle_*-\langle f\rangle_*\langle h\rangle_*.
\]

So the full coupled mouth layer selects a unique pair of moment shifts; no new
branch ambiguity is introduced.

---

## 3. Projection onto the rigidity kernel

The canonical lower branch remains defined by
\[
\mathfrak g=\mathfrak g_*,
\]
so the electrochemical bias retunes by
\[
\boxed{
\delta\Pi_{\rm act}
=
-\frac{\delta\mathfrak g_{\rm act}}{\mathfrak g_*'}
=
\frac{\operatorname{Cov}_*(c,R_*)}{\mathfrak g_*'}.
}
\]

The normalized mouth traction shift is
\[
\boxed{
\delta \widehat T_{m,{\rm act}}
=
A_T\,\delta\mathfrak g_{\rm act}
+
B_T\,\delta\mathcal S_{\rm act}
=
- A_T\,\operatorname{Cov}_*(c,R_*)
- B_T\,\operatorname{Cov}_*(K_q,R_*).
}
\]

So the full-mouth selected correction is now completely explicit:
it is the projection of the exact residual \(R_*(x)\) against the same rigidity
kernel already derived in Stage 130.

---

## 4. Meaning

The mouth-side ambiguity is now reduced to one actual physical object:

\[
\boxed{
R_*(x)
}
\]

Once \(R_*(x)\) is known from the channel solve, the induced source correction,
bias retuning, and traction retuning follow automatically.

So the next question is no longer “which positive family?”
It is simply: what are the actual numerical covariances of \(R_*(x)\) on the
explicit Family-1 branch?
=== moving_throat_pde_stage135_family1_actual_correction.md ===

# Moving-Throat PDE — Stage 135: Actual Family-1 Mouth Correction and One-Step Nonlinear Check

## Goal

Evaluate the exact full-profile residual \(R_*(x)\) on the explicit Family-1
canonical branch and project it against the rigidity formulas.

---

## 1. Actual first-order Family-1 correction

On the compensated Family-1 branch
\[
\Pi_* \approx 1.50882951349316,
\qquad
\Sigma_m^*\approx 0.451485277739090,
\qquad
\mathfrak g_* \approx 0.758035078944663,
\qquad
\mathcal S_* \approx 0.658075937605429,
\]
the exact weighted residual covariances are
\[
\boxed{
\operatorname{Cov}_*(c,R_*)\approx 0.0648069687666328,
}
\]
\[
\boxed{
\operatorname{Cov}_*(K_q,R_*)\approx 0.0388718368650403.
}
\]

Therefore the actual first-order moment shifts are
\[
\boxed{
\delta\mathfrak g_{\rm act}\approx -0.0648069687666328,
\qquad
\delta\mathcal S_{\rm act}\approx -0.0388718368650403.
}
\]

So the full mouth profile broadens the source relative to the tangent exponential:
the overlap factor and the mixed D/N response factor both move downward.

---

## 2. Retuned canonical point

Using the previously frozen canonical derivatives and rigidity coefficients,
\[
\mathfrak g_*' \approx 0.0714453558083195,
\qquad
A_T\approx -4.27263956256927,
\qquad
B_T\approx 0.134875005736706,
\]
the actual retuning is
\[
\boxed{
\delta\Pi_{\rm act}
\approx 0.907084414842908,
}
\]
\[
\boxed{
\delta\widehat T_{m,{\rm act}}
\approx 0.271653979462338.
}
\]

So the corrected canonical Family-1 point is
\[
\boxed{
\Pi_{\rm corr}
=
\Pi_*+\delta\Pi_{\rm act}
\approx 2.41591392833607,
}
\]
\[
\boxed{
\widehat T_{m,\rm corr}
=
\widehat T_{m,*}+\delta\widehat T_{m,\rm act}
\approx 1.17313803363654.
}
\]

This is the first actual non-exponential correction selected by the full explicit
mouth-layer model.

---

## 3. One-step nonlinear Picard check

A useful nonlinear check is to replace the exponential source by the one-step
full-profile iterate
\[
\boxed{
\Sigma_1(x)=
\frac{e^{-\Pi_*x-R_*(x)}}{\int_0^1 e^{-\Pi_*y-R_*(y)}\,dy}.
}
\]

This is not yet the full nonlinear fixed point, but it keeps the entire finite
residual \(R_*(x)\) rather than only its linearized projection.

Its actual moments are
\[
\boxed{
\mathfrak g_1\approx 0.684423574065325,
\qquad
\mathcal S_1\approx 0.616333130570251.
}
\]

Retuning the canonical branch with those exact one-step moments gives
\[
\boxed{
\Pi_1\approx 2.53914847609768,
\qquad
\widehat T_{m,1}\approx 1.21036942084359.
}
\]

So the one-step nonlinear correction shifts the canonical point slightly more than
the linearized estimate, but in the same direction and by the same overall scale.

---

## 4. Effective positive-family interpretation

Comparing the actual first-order correction to the explicit positive-family
interpolation from Stage 131 gives
\[
\lambda_{\rm eff}^{(\Pi)}\approx 0.380487632771110,
\qquad
\lambda_{\rm eff}^{(T)}\approx 0.378939241176339.
\]

So the full coupled mouth-layer correction is well approximated by a point about
\[
\boxed{
\lambda_{\rm eff}\approx 0.38
}
\]
on the positive interpolation line between the uniform family and the self-matched
derivative family.

Equivalently, the actual selected correction corresponds to a **broadening fraction**
\[
1-\lambda_{\rm eff}\approx 0.62.
\]

So the full mouth profile behaves much more like a moderate broadening toward
uniformity than like a sharpening toward the self-matched derivative branch.
=== moving_throat_pde_stage136_full_mouth_correction_status.md ===

# Moving-Throat PDE — Stage 136 / v58 Status

## New closure result

The full explicit GNLS + localized-Maxwell Family-1 mouth boundary layer now selects
a **definite non-exponential correction** around the unique regular lower compensated
branch.

The exact compensated mouth potential
\[
\Phi_*(x)=4\Sigma_m^*T_s(x;\Pi_*)-\Sigma_m^*T_q(x;\Pi_*)
\]
defines a residual
\[
R_*(x)=\Phi_*(x)-\Pi_*x
\]
with
\[
R_*(0)=R_*'(0)=0,
\qquad
R_*''(0)=-3\Sigma_m^*\frac{\Pi_*}{1-e^{-\Pi_*}}<0.
\]
So the full mouth profile is tangent-matched but sublinear at the mouth, and therefore
broadens the actual source relative to the tangent exponential branch.

## Actual first-order Family-1 correction

Projecting the exact residual against the Stage-130 rigidity kernel gives
\[
\delta\mathfrak g_{\rm act}\approx -0.0648069688,
\qquad
\delta\mathcal S_{\rm act}\approx -0.0388718369,
\]
\[
\delta\Pi_{\rm act}\approx 0.9070844148,
\qquad
\delta\widehat T_{m,\rm act}\approx 0.2716539795.
\]

So the corrected canonical point is
\[
\Pi_{\rm corr}\approx 2.4159139283,
\qquad
\widehat T_{m,\rm corr}\approx 1.1731380336.
\]

## One-step nonlinear check

Using the one-step fully exponentiated profile
\[
\Sigma_1(x)\propto e^{-\Pi_*x-R_*(x)}
\]
gives
\[
\mathfrak g_1\approx 0.6844235741,
\qquad
\mathcal S_1\approx 0.6163331306,
\]
and a retuned point
\[
\Pi_1\approx 2.5391484761,
\qquad
\widehat T_{m,1}\approx 1.2103694208.
\]

So the nonlinear correction is slightly stronger than the linear estimate but follows
the same direction and scale.

## Interpretation

Inside the explicit Family-1 closure:

- branch selection is finished;
- the lower compensated branch remains the unique regular branch;
- the full coupled mouth profile does **not** destroy it;
- but it shifts the preferred point upward in both bias and normalized traction.

The mouth-side problem is therefore no longer one of branch ambiguity. It is now a
finite correction problem around a unique regular branch.

## Next serious step

The next PDE-facing step is to let the **core outlet coefficients and the mouth profile
co-evolve** self-consistently, instead of holding the compensated core branch fixed while
correcting only the mouth source law.
=== moving_throat_pde_stage137_coevolving_core_mouth_map.md ===

# Moving-Throat PDE — Stage 137: Exact Co-Evolving Core–Mouth Fixed-Point Map

## Goal

Replace the “fixed compensated core + corrected mouth source” closure by a **single
self-consistent co-evolving core–mouth map**.

The point is to let the actual positive mouth source profile modify the shell/mixed
loading ratio, and then feed that modified ratio back into the mouth potential.

---

## 1. The exact co-evolving Family-1 map

For any normalized positive mouth source
\[
\Sigma(x)\ge 0,
\qquad
\int_0^1 \Sigma(x)\,dx = 1,
\]
define the two mouth moments
\[
\mathfrak g[\Sigma]
=
\int_0^1 \Sigma(x)\cos\!\left(\frac{\pi x}{2}\right)dx,
\qquad
\mathcal S[\Sigma]
=
\int_0^1 \Sigma(x)\,
\frac{\cosh\!\left(\frac{\pi}{2}(1-x)\right)}{\cosh(\pi/2)}\,dx.
\]

On the explicit Family-1 core branch, the shell/mixed loading ratio is not free.
It is the exact core overlap law
\[
\boxed{
\mathcal R[\Sigma]
=
\frac{\bigl(\mathfrak g[\Sigma]-\mathfrak r_{F1}\bigr)^2}{1+\mathfrak r_{F1}^2},
\qquad
\mathfrak r_{F1}\approx 1.77799353547498.
}
\]

For the static shell channel,
\[
\mathcal T_s[\Sigma](x)=\int_0^1 \min(x,y)\,\Sigma(y)\,dy,
\]
and for the first mixed D/N half-wave,
\[
\mathcal T_q[\Sigma](x)
=
\int_0^1
\frac{\sinh\!\left(\frac{\pi}{2}\min(x,y)\right)
\cosh\!\left(\frac{\pi}{2}(1-\max(x,y))\right)}{(\pi/2)\cosh(\pi/2)}
\Sigma(y)\,dy.
\]

So the full co-evolving mouth potential is
\[
\boxed{
\Phi_{\Sigma_0}[\Sigma](x)
=
\Sigma_0
\Big[
\mathcal T_s[\Sigma](x)
-
\mathcal R[\Sigma]\,\mathcal T_q[\Sigma](x)
\Big].
}
\]

The self-consistent source law is therefore the nonlinear fixed-point equation
\[
\boxed{
\Sigma(x)
=
\frac{e^{-\Phi_{\Sigma_0}[\Sigma](x)}}
{\int_0^1 e^{-\Phi_{\Sigma_0}[\Sigma](y)}dy}.
}
\]

This is the first explicit Family-1 mouth problem in which the **source profile and
the core loading ratio co-evolve together**.

---

## 2. Canonical compensation inside the co-evolving map

The lower compensated throat-core branch remains the condition
\[
\mathfrak g[\Sigma]=\mathfrak g_*,
\qquad
\mathfrak g_*\approx 0.758035078944663.
\]

Because
\[
\mathfrak g_*=
\mathfrak r_{F1}-\frac12\sqrt{1+\mathfrak r_{F1}^2},
\]
this is exactly equivalent to
\[
\boxed{
\mathcal R[\Sigma]=\frac14.
}
\]

So the canonical outgoing quadrupole branch survives in the co-evolving theory iff
the self-consistent source profile lands back on the same lower compensation moment.

---

## 3. Exact first-order defect transport

Write
\[
\mathfrak g=\mathfrak g_*+\delta\mathfrak g.
\]
Then the exact Family-1 ratio becomes
\[
\boxed{
\mathcal R
=
\frac14
-
\frac{\delta\mathfrak g}{\sqrt{1+\mathfrak r_{F1}^2}}
+
\frac{(\delta\mathfrak g)^2}{1+\mathfrak r_{F1}^2}.
}
\]

So to first order,
\[
\boxed{
\delta\mathcal R
=
-\frac{\delta\mathfrak g}{\sqrt{1+\mathfrak r_{F1}^2}}
+O(\delta\mathfrak g^2).
}
\]

Numerically,
\[
\sqrt{1+\mathfrak r_{F1}^2}
\approx 2.039916913060632,
\qquad
\delta\mathcal R
\approx -0.490215\,\delta\mathfrak g.
\]

So broadening the source (negative \(\delta\mathfrak g\)) automatically drives the
mixed loading ratio **above** its compensated value \(1/4\).

---

## 4. Local slope / bias identity

For any self-consistent source on this branch, the actual mouth slope is
\[
\boxed{
\Pi[\Sigma]
=
\Phi'_{\Sigma_0}[\Sigma](0)
=
\Sigma_0\Bigl[1-\mathcal R[\Sigma]\,\mathcal S[\Sigma]\Bigr].
}
\]

Hence the co-evolving canonical branch is determined by the coupled conditions

\[
\Sigma
=
\frac{e^{-\Phi_{\Sigma_0}[\Sigma]}}
{\int e^{-\Phi_{\Sigma_0}[\Sigma]}},
\qquad
\mathfrak g[\Sigma]=\mathfrak g_*,
\qquad
\Pi=\Sigma_0\left(1-\frac14\mathcal S[\Sigma]\right).
\]

Under the self-matched susceptibility closure from Stage 123,
\[
\boxed{
\Sigma_0=\frac{20}{9}\widehat T_m^2.
}
\]

So the remaining ambiguity is now only the required **normalized mouth traction**
\(\widehat T_m\) that makes the co-evolving fixed point land on
\(\mathfrak g_*\).

=== moving_throat_pde_stage138_frozen_traction_fixedpoint.md ===

# Moving-Throat PDE — Stage 138: Family-1 Co-Evolving Fixed Point at Frozen Canonical Traction

## Goal

Solve the exact co-evolving core–mouth map on the explicit Family-1 branch while
holding the previously selected canonical traction fixed:
\[
\Sigma_0=\Sigma_0^*
\approx 1.80594111095636
\qquad
\left(
\widehat T_m=\widehat T_{m,*}
\approx 0.901484054174204
\right).
\]

This answers a very specific question:

> If the physical mouth traction is left at the old canonical value, does the
> fully co-evolving Family-1 fixed point still land on the lower compensated
> branch \(\mathcal R=1/4\)?

---

## 1. Unique positive fixed point

Iterating the exact nonlinear map
\[
\Sigma\mapsto
\frac{e^{-\Phi_{\Sigma_0^*}[\Sigma]}}
{\int_0^1 e^{-\Phi_{\Sigma_0^*}[\Sigma](y)}dy}
\]
from the canonical exponential seed converges to a unique positive fixed point
\(\Sigma_{\rm fp}\) on the explicit Family-1 branch.

Its selected moments are
\[
\boxed{
\mathfrak g_{\rm fp}
\approx 0.693352419668063,
\qquad
\mathcal S_{\rm fp}
\approx 0.6216013167514007.
}
\]

So relative to the exact lower compensated value
\[
\mathfrak g_*\approx 0.758035078944663,
\]
the fixed profile is broader:
\[
\boxed{
\delta\mathfrak g_{\rm fp}
=
\mathfrak g_{\rm fp}-\mathfrak g_*
\approx -0.0646826592766000<0.
}
\]

---

## 2. Co-evolving core ratio and mouth slope

Feeding that back into the exact Family-1 core law gives
\[
\boxed{
\mathcal R_{\rm fp}
=
\mathcal R[\Sigma_{\rm fp}]
\approx 0.2827139049082381.
}
\]

So the co-evolving fixed point **does not** stay on the compensated value
\(1/4\); it shifts upward by
\[
\boxed{
\delta\mathcal R_{\rm fp}
=
\mathcal R_{\rm fp}-\frac14
\approx 0.0327139049082381.
}
\]

The exact slope identity then gives
\[
\boxed{
\Pi_{\rm fp}
=
\Sigma_0^*
\Bigl[1-\mathcal R_{\rm fp}\mathcal S_{\rm fp}\Bigr]
\approx 1.4885734438300713.
}
\]

So at fixed canonical traction the co-evolving mouth profile actually stays very
close to the old canonical bias:
\[
\boxed{
\delta\Pi_{\rm fp}
=
\Pi_{\rm fp}-\Pi_*
\approx -0.0202560696630887.
}
\]

---

## 3. Exact first-order transport check

Using the exact transport law from Stage 137,
\[
\delta\mathcal R
=
-\frac{\delta\mathfrak g}{\sqrt{1+\mathfrak r_{F1}^2}}
+
\frac{(\delta\mathfrak g)^2}{1+\mathfrak r_{F1}^2},
\]
the observed fixed-point broadening predicts
\[
\delta\mathcal R_{\rm pred}
=
-\frac{\delta\mathfrak g_{\rm fp}}{\sqrt{1+\mathfrak r_{F1}^2}}
+
\frac{(\delta\mathfrak g_{\rm fp})^2}{1+\mathfrak r_{F1}^2}
pprox 0.0327139049082381,
\]
which matches the direct fixed-point value.

So the fixed-point drift is exactly the one implied by the co-evolving core law.

---

## 4. Meaning

The full Family-1 core–mouth solve gives a more nuanced answer than the earlier
fixed-core mouth correction:

- the source profile really does broaden,
- but the induced increase in the mixed loading ratio pushes back on the mouth bias,
- so at **fixed** canonical traction the self-consistent branch remains close in
  \(\Pi\), even though it is no longer exactly compensated in \(\mathfrak g\).

So the co-evolution does **not** destroy the regular Family-1 branch.
What it does is convert the old exact compensation point into a nearby
traction-dependent fixed point with
\[
\mathcal R_{\rm fp}>1/4.
\]

That means exact preservation of the canonical outgoing quadrupole fingerprint now
requires a **retuned traction**, not just the old canonical value.

=== moving_throat_pde_stage139_renormalized_canonical_branch.md ===

# Moving-Throat PDE — Stage 139: Renormalized Canonical Branch Under Full Core–Mouth Co-Evolution

## Goal

Determine the actual Family-1 traction required to restore the exact lower
compensated branch once the mouth profile and the core loading ratio are allowed
to co-evolve together.

Equivalently, solve the exact self-consistent condition
\[
\mathfrak g[\Sigma_{\Sigma_0}]=\mathfrak g_*,
\]
rather than imposing the old canonical traction by hand.

---

## 1. Unique renormalized canonical gain

Scanning the exact self-consistent Family-1 map shows that the fixed-point moment
\(\mathfrak g_{\rm fp}(\Sigma_0)\) increases monotonically over the physical
interval used in the solve. Therefore the compensation-restoration equation
\[
\mathfrak g_{\rm fp}(\Sigma_0)=\mathfrak g_*
\]
has a unique positive root.

That root is
\[
\boxed{
\Sigma_0^{\rm can}
\approx 4.651033550168876.
}
\]

Using the self-matched susceptibility closure,
\[
\Sigma_0=\frac{20}{9}\widehat T_m^2,
\]
this corresponds to the renormalized canonical traction
\[
\boxed{
\widehat T_{m,\rm can}
\approx 1.446708366456762.
}
\]

So the full co-evolving theory does preserve the lower compensated branch, but only
after a substantial upward retuning of the mouth traction.

---

## 2. The restored canonical fixed point

At that renormalized traction, the exact self-consistent fixed point satisfies
\[
\boxed{
\mathfrak g_{\rm can}=\mathfrak g_*
\approx 0.758035078944663,
\qquad
\mathcal R_{\rm can}=\frac14.
}
\]

Its selected mixed response and mouth bias are
\[
\boxed{
\mathcal S_{\rm can}
\approx 0.6703621156734617,
\qquad
\Pi_{\rm can}
\approx 3.8715643774790087.
}
\]

So the co-evolving compensated Family-1 branch is still perfectly regular and still
finite — it just sits at a higher traction/bias point than the earlier
fixed-core correction suggested.

---

## 3. Comparison with earlier mouth-only corrections

The earlier fixed-core mouth analysis gave approximately
\[
(\Pi_{\rm corr},\widehat T_{m,\rm corr})
\approx (2.4159,1.1731),
\]
with a one-step nonlinear check closer to
\[
(\Pi_1,\widehat T_{m,1})
\approx (2.5391,1.2104).
\]

The full co-evolving canonical solve instead gives
\[
\boxed{
(\Pi_{\rm can},\widehat T_{m,\rm can})
\approx
(3.8716,1.4467).
}
\]

Relative to the original canonical point
\[
(\Pi_*,\widehat T_{m,*})
\approx (1.5088,0.9015),
\]
the exact co-evolving compensation costs

\[
\boxed{
\frac{\Sigma_0^{\rm can}}{\Sigma_0^*}-1
\approx 1.5754070949223031,
}
\]

\[
\boxed{
\frac{\Pi_{\rm can}}{\Pi_*}-1
\approx 1.5659389234213572,
}
\]

\[
\boxed{
\frac{\widehat T_{m,\rm can}}{\widehat T_{m,*}}-1
\approx 0.6048074946616844.
}
\]

So exact preservation of the canonical outgoing quadrupole branch under full
core–mouth co-evolution requires roughly a

- \(60.48\%\) traction increase,
- and a \(156.59\%\) mouth-bias increase,

relative to the original lower compensated point.

---

## 4. Meaning

This is the sharpest Family-1 mouth/core result so far.

Inside the explicit Family-1 closure:

1. the lower compensated branch **survives** full co-evolution,
2. it remains the unique physically admissible compensated branch,
3. but the old canonical point is not the final self-consistent point,
4. and exact compensation is recovered only on a **renormalized** finite-bias,
   finite-traction branch.

So the mouth-side problem is no longer branch selection and no longer a broad
profile ambiguity. It is now a quantitative renormalization problem for the
unique regular canonical branch.

=== moving_throat_pde_stage140_core_mouth_coevolution_status.md ===

# Moving-Throat PDE — Stage 140 / v59 Status

## New closure result

The explicit Family-1 GNLS + localized-Maxwell program now has a **fully co-evolving
core–mouth closure**.

For a positive normalized source \(\Sigma\), the source moment
\[
\mathfrak g[\Sigma]
=
\int_0^1 \Sigma(x)\cos\!\left(\frac{\pi x}{2}\right)dx
\]
feeds the Family-1 core loading ratio
\[
\mathcal R[\Sigma]
=
\frac{(\mathfrak g[\Sigma]-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2},
\qquad
\mathfrak r_{F1}\approx 1.77799353547498,
\]
which then re-enters the full mouth potential
\[
\Phi_{\Sigma_0}[\Sigma](x)
=
\Sigma_0\Bigl[
\mathcal T_s[\Sigma](x)-\mathcal R[\Sigma]\mathcal T_q[\Sigma](x)
\Bigr].
\]

So the explicit Family-1 problem is now an honest nonlinear fixed point, not a
frozen-core correction.

## What the full solve says

At the old canonical traction
\[
\Sigma_0^*\approx 1.80594111095636
\quad
\left(
\widehat T_{m,*}\approx 0.901484054174204
\right),
\]
the co-evolving fixed point lands at
\[
\mathfrak g_{\rm fp}\approx 0.6933513438014128,
\qquad
\mathcal R_{\rm fp}\approx 0.2827144657621617,
\qquad
\Pi_{\rm fp}\approx 1.4885731675167777.
\]
So the branch survives and stays close in bias, but it is no longer exactly
compensated.

Restoring the exact lower compensated branch requires the unique renormalized
traction
\[
\Sigma_0^{\rm can}\approx 4.6510849445826352,
\qquad
\widehat T_{m,\rm can}\approx 1.4467163595750847,
\qquad
\Pi_{\rm can}\approx 3.8716072210898429.
\]

## Interpretation

The mouth/core ambiguity is now essentially gone inside the explicit Family-1
closure:

- the upper compensated branch is impossible,
- the equal-normalized branch is singular,
- the lower compensated branch remains the only regular canonical branch,
- and full self-consistency promotes it to a renormalized finite-traction,
  finite-bias fixed point.

## Next serious step

The next PDE-facing task is no longer branch selection. It is to derive the actual
deviation of the moving-throat mouth/core system from this explicit Family-1
co-evolving fixed point, and then translate that deviation into the remaining
outgoing quadrupole-normalization defect.
# Moving-Throat PDE Engine
## Raw 4D GNLS + localized Maxwell + moving-throat geometry, with the grouped real `P2` bridge

## 0. Why this document exists

This document is the **engine** handoff. It is meant to be sufficient for a fresh session to reconstruct the actual moving-throat PDE program without needing the conversational history.

It does four things:

1. states the exact **4D parent action** and the exact bulk equations already fixed by the program,
2. defines the **topological defect / throat** as a moving brane–bulk surface,
3. records the **boundary data and mode structure** needed for the finite-throat internal branch,
4. and states the reduced grouped real `P2` output that the completed PDE must eventually supply.

The document is intentionally explicit about **claim status**.

- **Exact** means it follows directly from the declared action, exact definitions, or exact algebra.
- **Controlled reduction** means it follows only after a stated ansatz or low-frequency/small-body reduction.
- **Protocol closure** means it is fixed only within the declared response hierarchy, not yet from a fully solved moving-throat PDE.
- **Open** means the full PDE still has to decide it.

The strongest current reading is:

- the parent 4D field theory is exact,
- the moving-throat geometry lift and linearized coupled bulk/interface skeleton are explicit,
- the grouped real `P2` conservative/output bridge is explicit,
- and the last serious unresolved theorem gate is still the **actual branch selection and normalization** of the outgoing quadrupole channel on the true moving-throat solution.

---

## 1. Core ontology and kinematic stage

### 1.1 Bulk spacetime and indices

The fundamental arena is a `(4+1)`-dimensional spacetime with coordinates
\[
x^M=(t,x,y,z,w),
\qquad
M,N\in\{0,1,2,3,4\}.
\]

The bulk spatial coordinates are
\[
\mathbf X=(x,y,z,w)\in\mathbb R^4,
\]
and the brane spatial coordinates are
\[
\mathbf x=(x,y,z)\in\mathbb R^3.
\]

Brane indices are
\[
\mu,\nu\in\{0,1,2,3\},
\]
bulk spatial indices are
\[
i,j\in\{x,y,z,w\},
\]
and brane spatial indices are
\[
a,b\in\{x,y,z\}.
\]

The flat bulk metric is
\[
\eta_{MN}=\mathrm{diag}(-1,+1,+1,+1,+1),
\]
with d’Alembertian
\[
\Box_5=-\partial_t^2+\nabla_3^2+\partial_w^2.
\]

### 1.2 Fields

The core dynamical variables are

\[
\psi(\mathbf X,t),\qquad A_M(x),\qquad \Sigma(\mathbf X,t),
\]
or equivalently the shape field
\[
R(\Omega,w,t)\quad\text{with}\quad \Sigma=r-R(\Omega,w,t),
\]
where
\[
r=\sqrt{x^2+y^2+z^2},\qquad \Omega=\mathbf x/r\in S^2.
\]

The bulk density is
\[
\rho=|\psi|^2.
\]

The localized gauge field is
\[
A_M=(A_0,A_i),\qquad
F_{MN}=\partial_M A_N-\partial_N A_M.
\]

The old finite-dimensional throat variables survive only as collective moments:
\[
a(t),\qquad L(t).
\]

### 1.3 Corrected charge ontology

The corrected electric-charge bookkeeping is
\[
\eta_Q=\pm 1,\qquad
q_\star=\eta_Q e_\star,\qquad
e_\star>0.
\]

After canonical zero-mode brane normalization,
\[
q_{\rm eff}=\frac{q_\star}{\sqrt{Z_{\rm int}}},
\qquad
e_{\rm eff}=\frac{e_\star}{\sqrt{Z_{\rm int}}},
\qquad
Z_{\rm int}=\int_{-\infty}^{+\infty} Z(w)\,dw.
\]

Important firewall:

- electric-charge sign is carried by \(\eta_Q\),
- circulation belongs to the magnetic/vortical sector,
- the historical gravity-side `q=1` is the mass-dressing coefficient \(\kappa_\rho=1\), not electric charge.

---

## 2. Exact 4D parent action

### 2.1 Matter sector: gauged 4D GNLS

The exact matter Lagrangian density is
\[
\boxed{
\mathcal L_\psi
=
\frac{i\hbar}{2}\left(\psi^*D_t\psi-\psi D_t\psi^*\right)
-\frac{\hbar^2}{2m}(D_i\psi)^*(D_i\psi)
- V_{\rm conf}(\mathbf X;\Sigma)\,\rho
- U(\rho).
}
\]

The covariant derivatives are
\[
D_t\psi=\partial_t\psi+\frac{i q_\star}{\hbar}A_0\psi,
\qquad
D_i\psi=\partial_i\psi-\frac{i q_\star}{\hbar}A_i\psi.
\]

The frozen stiff-polytropic equation of state is
\[
P(\rho)=K\rho^5,
\qquad
U(\rho)=\frac{K}{4}\rho^5,
\qquad
h(\rho)=\frac{dU}{d\rho}=\frac{5K}{4}\rho^4,
\]
so the bulk sound speed is
\[
c_s^2(\rho)=\frac{1}{m}\frac{dP}{d\rho}=\frac{5K}{m}\rho^4.
\]

### 2.2 Localized Maxwell sector

The exact localized Maxwell Lagrangian density is
\[
\boxed{
\mathcal L_{\rm EM}
=
-\frac{Z(w)}{4\mu_0}F_{MN}F^{MN}
-\frac{1}{2\xi\mu_0}(\partial\!\cdot\!A)^2
- A_M J_{\rm ext}^M.
}
\]

The total source entering Maxwell’s equation is
\[
J_{\rm tot}^M = J_\psi^M + J_{\rm ext}^M.
\]

The key bookkeeping rule is that varying the covariant matter action already generates the dynamical matter current \(J_\psi^M\), so one must not double-count it in the explicit EM source term.

### 2.3 Legacy finite-dimensional geometry closure

The original reduced geometry sector was
\[
\mathcal L_{\rm geom}
=
\frac12 M_a \dot a^2 + \frac12 M_L \dot L^2
- E_{\rm geom}(a,L) - E_{\rm ratio}(a,L),
\]
with optional penalty
\[
E_{\rm ratio}(a,L)=\frac{\kappa}{2}(L-\alpha a)^2,
\]
and damped laws
\[
M_a\ddot a+\Gamma_a\dot a = -\frac{\partial H_{\rm tot}}{\partial a},
\qquad
M_L\ddot L+\Gamma_L\dot L = -\frac{\partial H_{\rm tot}}{\partial L}.
\]

This survives only as the **lowest collective truncation** of the moving-throat field.

---

## 3. Exact bulk equations already fixed by the parent theory

### 3.1 Gauged 4D GNLS equation

Variation with respect to \(\psi^\ast\) gives
\[
\boxed{
i\hbar D_t\psi
=
\left[
-\frac{\hbar^2}{2m}D_iD_i
+V_{\rm conf}(\mathbf X;\Sigma)
+h(\rho)
\right]\psi.
}
\]

### 3.2 Exact current and continuity

The bulk number current is
\[
\boxed{
j^i=\frac{\hbar}{m}\,\Im(\psi^\ast D_i\psi),
}
\]
and exact continuity is
\[
\boxed{
\partial_t\rho+\partial_i j^i=0.
}
\]

Where \(\rho>0\), define the bulk velocity by
\[
j^i=\rho v^i.
\]

### 3.3 Exact localized Maxwell equation

Variation with respect to \(A_N\) gives
\[
\boxed{
\partial_M\!\left(Z(w)F^{MN}\right)
+\frac{1}{\xi}\partial^N(\partial\!\cdot\!A)
=
\mu_0 J_{\rm tot}^N.
}
\]

The Bianchi identities are
\[
\partial_{[L}F_{MN]}=0.
\]

Gauge invariance requires current conservation:
\[
\partial_M J_{\rm tot}^M = 0.
\]

### 3.4 Madelung rewrite and Euler-like form

Write
\[
\psi=\sqrt{\rho}\,e^{i\theta}.
\]

The gauge-invariant velocity is
\[
\boxed{
v_i=\frac{\hbar}{m}\partial_i\theta-\frac{q_\star}{m}A_i.
}
\]

The quantum potential is
\[
\boxed{
Q(\rho)=-\frac{\hbar^2}{2m}\frac{\nabla_4^2\sqrt{\rho}}{\sqrt{\rho}}.
}
\]

The exact Euler-like bulk equation is
\[
\boxed{
m(\partial_t+v_j\partial_j)v_i
=
q_\star(E_i+v_j B_{ij})
-\partial_i\!\left(V_{\rm conf}+h(\rho)+Q(\rho)\right),
}
\]
with
\[
E_i=-\partial_tA_i-\partial_iA_0,\qquad
B_{ij}=F_{ij}.
\]

### 3.5 Exact vorticity–gauge identity

Define the bulk vorticity 2-form
\[
\Omega_{ij}=\partial_i v_j-\partial_j v_i.
\]

Away from phase singularities,
\[
\boxed{
\Omega_{ij}=-\frac{q_\star}{m}F_{ij}.
}
\]

This is why circulation belongs to the magnetic/vortical sector rather than the electric-charge dictionary.

### 3.6 Exact mixed-sector gauge invariants

The mixed fields
\[
E_w = F_{w0} = -\partial_t A_w - \partial_w A_0,
\]
\[
C_a = F_{aw} = \partial_a A_w - \partial_w A_a,
\]
are exact gauge invariants under
\[
A_0\to A_0-\partial_t\chi,\qquad
A_a\to A_a+\partial_a\chi,\qquad
A_w\to A_w+\partial_w\chi.
\]

These mixed channels are **suppressed only in the strict far-field zero-mode brane reduction**. They remain part of the microscopic ontology and are essential for the honest outgoing bridge.

---

## 4. Projection versus reduction

### 4.1 Projection: exact brane observables

For a normalized brane weight \(W(w)\),
\[
\int W(w)\,dw=1,
\]
the projected brane observables are
\[
\rho_{\rm brane}(\mathbf x,t)=\int W(w)\rho(\mathbf x,w,t)\,dw,
\]
\[
\mathbf j_{\rm brane}(\mathbf x,t)=\int W(w)\mathbf j_{xyz}(\mathbf x,w,t)\,dw,
\]
\[
\mathbf v_{\rm brane}=\mathbf j_{\rm brane}/\rho_{\rm brane}.
\]

Projected continuity is exact:
\[
\boxed{
\partial_t\rho_{\rm brane}+\nabla_3\cdot \mathbf j_{\rm brane}=S_{\rm leak},
}
\]
with
\[
\boxed{
S_{\rm leak}
=
-\left[Wj^w\right]_{-\infty}^{+\infty}
+\int W'(w)j^w\,dw.
}
\]

### 4.2 Poisson hook

With the Helmholtz decomposition
\[
\mathbf v_{\rm brane}=\nabla_3\varphi+\mathbf v_T,\qquad \nabla_3\cdot \mathbf v_T=0,
\]
one obtains the exact longitudinal identity
\[
\boxed{
\rho_{\rm brane}\,\nabla_3^2\varphi
=
S_{\rm leak}
-\partial_t\rho_{\rm brane}
-(\nabla_3\rho_{\rm brane})\cdot(\nabla_3\varphi+\mathbf v_T).
}
\]

In the quasi-static longitudinal regime this reduces to a Poisson equation on the brane.

### 4.3 Zero-mode Maxwell reduction

Under the controlled far-field zero-mode assumptions
\[
A_w\approx 0,\qquad
\partial_w A_\mu\approx 0,\qquad
J^w\approx 0,\qquad
F_{\mu w}\approx 0,
\]
integration over \(w\) gives an effective brane Maxwell sector
\[
\partial_\mu F^{\mu\nu}=\mu_0^{\rm eff}J_{\rm eff}^\nu,
\qquad
\mu_0^{\rm eff}=\frac{\mu_0}{Z_{\rm int}}.
\]

This is a **controlled reduction**, not a denial of the mixed-core structure.

---

## 5. Topological defect / throat definition

### 5.1 Moving throat as a level-set surface

The smallest honest moving-throat lift is
\[
\boxed{
\Sigma(\mathbf X,t)=r-R(\Omega,w,t).
}
\]

The finite throat surface is
\[
\Sigma(\mathbf X,t)=0.
\]

The sign convention is

- exterior: \(\Sigma>0\),
- interior/support region: \(\Sigma<0\).

The reference stationary throat is
\[
\Sigma_0(\mathbf X)=r-R_0(w).
\]

### 5.2 Reference stationary throat

The reference branch is a finite throat with

- mouth at \(w=0\),
- finite depth \(0\le w\le L_0\),
- mouth radius
  \[
  a_0=R_0(0),
  \]
- bottom closure through either
  \[
  R_0(L_0)=0
  \]
  or an equivalent regular bottom condition.

The old collective variables are recovered as averages:

\[
a(t)=\frac{1}{4\pi}\int_{S^2}R(\Omega,0,t)\,d\Omega,
\]
and if \(W_b(\Omega,t)\) is the bottom defined by \(R(\Omega,W_b(\Omega,t),t)=0\),
\[
L(t)=\frac{1}{4\pi}\int_{S^2}W_b(\Omega,t)\,d\Omega.
\]

### 5.3 Moving-wall confinement coupling

Promote the confinement potential to the level-set form
\[
\boxed{
V_{\rm conf}(\mathbf X;\Sigma)=V_{\rm wall}\!\left(\frac{\Sigma(\mathbf X,t)}{\ell_c}\right).
}
\]

Linearizing around the reference throat,
\[
R(\Omega,w,t)=R_0(w)+\eta(\Omega,w,t),
\]
gives
\[
\delta\Sigma=-\eta,
\]
hence
\[
\boxed{
\delta V_{\rm conf}
=
-\frac{V'_{\rm wall}(\Sigma_0/\ell_c)}{\ell_c}\,\eta.
}
\]

This is the direct wall–bulk coupling that drives the linearized matter and gauge sectors.

---

## 6. Throat mode decomposition and harmonic content

### 6.1 Real spherical-harmonic decomposition on the mouth sphere

Write the wall displacement as
\[
\eta(\Omega,w,t)
=
\eta_0(w,t)Y_{00}(\Omega)
+
\sum_{m\in P_2({\rm real})} q_{2m}(w,t)\,Y_{2m}^{\rm real}(\Omega)
+
\eta_{\ge 3}(\Omega,w,t).
\]

The grouped real `P2` set is
\[
\{20,\ 21c,\ 21s,\ 22c,\ 22s\}.
\]

### 6.2 Monopole normalization bridge

With
\[
Y_{00}=\frac{1}{2\sqrt\pi},
\]
the physical mouth-average shift \(\delta a\) is related to the normalized monopole coefficient by
\[
\boxed{
q_{00}(0,t)=2\sqrt\pi\,\delta a(t).
}
\]

A useful axisymmetric split is
\[
\eta_0(w,t)=\alpha_a(w)\,\delta a(t)+\alpha_L(w)\,\delta L(t)+g(w,t),
\]
where \(g(w,t)\) is the remaining axisymmetric geometry lane orthogonal to the collective \(a\) and \(L\) moments.

### 6.3 Why the grouped real `P2` bundle matters

The grouped real `P2` sector is the first nontrivial harmonic family beyond the monopole. It is the literal moving-throat realization of the conservative grouped real quadrupole payload already exposed by the 2PN/3PN/2.5PN hierarchy.

---

## 7. Finite-throat support branch and boundary data

### 7.1 Exact D/N finite-throat support problem

On the finite throat interval
\[
z\in[0,L],
\]
the minimal internal support equation is
\[
\psi''+k^2\psi=0.
\]

The finite-throat D/N branch imposes

- **Dirichlet** at the mouth:
  \[
  \boxed{\psi(0)=0,}
  \]
- **Neumann** at the bottom:
  \[
  \boxed{\psi'(L)=0.}
  \]

The eigenmodes are
\[
\psi_j(z)=A_j\sin(k_j z),
\qquad
k_j=\frac{\pi}{L}\left(j+\frac12\right),
\qquad
j=0,1,2,\dots
\]
with frequencies
\[
\boxed{
\omega_j=c_s k_j
=
\frac{\pi c_s}{L}\left(j+\frac12\right).
}
\]

### 7.2 Mouth DtN operator

For prescribed mouth datum \(\psi_m\), the finite-throat D/N branch gives the exact mouth derivative
\[
\boxed{
Z_{00}(\omega)
=
-\frac{\omega}{c_s}\tan\!\left(\frac{\omega L}{c_s}\right).
}
\]

Its poles are exactly the D/N ladder above.

### 7.3 Trapped-branch round-trip closure

The scalar round-trip factor is
\[
R_{\rm rt}=r_0 r_L e^{2ikL}.
\]

For the D/N branch,
\[
r_D=-1,\qquad r_N=+1,
\]
so on the exact D/N ladder
\[
R_{\rm rt}=1,
\qquad
\phi_0\equiv 0 \pmod{2\pi}.
\]

This is the exact trapped-support closure used downstream.

---

## 8. Distributed wall action and reduced wall sectors

### 8.1 Minimal distributed wall action

The smallest passive distributed geometry action used in the moving-throat program is
\[
\boxed{
S_\eta^{(2)}
=
\frac12\int dt\,dw\,d\Omega\,\sqrt{\gamma_0}
\left[
\mu_\eta(w)(\partial_t\eta)^2
-
T_w(w)(\partial_w\eta)^2
-
T_\Omega(w)\,\eta(-\Delta_{S^2})\eta
-
K_\eta(w)\eta^2
\right].
}
\]

This yields the modal operator
\[
\mu_\eta q_{lm,tt}
-
\partial_w(T_w \partial_w q_{lm})
+
\bigl[K_\eta+l(l+1)T_\Omega\bigr]q_{lm}
=
S_{lm}^{(\psi,A)}+f_{lm}^{\rm ext}.
\]

### 8.2 Axisymmetric reduction back to \((a,L)\)

Using the two-mode axisymmetric truncation
\[
\eta_0(w,t)=2\sqrt\pi\,[\alpha_a(w)\delta a(t)+\alpha_L(w)\delta L(t)],
\]
one recovers the old reduced matrix system
\[
L_{\rm red}^{(0)}
=
\frac12 M_{AB}\dot Q^A\dot Q^B
-
\frac12 K_{AB}Q^A Q^B,
\qquad Q^A=(\delta a,\delta L),
\]
with
\[
M_{AB}=4\pi\int dw\,\mu_\eta\,\alpha_A\alpha_B,
\]
\[
K_{AB}=4\pi\int dw\,[T_w \alpha_A'\alpha_B' + K_0 \alpha_A\alpha_B].
\]

### 8.3 Grouped real `P2` reduction

For one grouped real quadrupole component,
\[
\eta_{2m}(\Omega,w,t)=\beta_2(w) q_{2m}(t) Y_{2m}^{\rm real}(\Omega),
\]
one obtains
\[
L_{2m}=\frac12 M_2 \dot q_{2m}^2 - \frac12 K_2 q_{2m}^2,
\]
with
\[
M_2=\int dw\,\mu_\eta \beta_2^2,
\]
\[
K_2=\int dw\,[T_w(\beta_2')^2 + (K_\eta+6T_\Omega)\beta_2^2].
\]

On an isotropic reference throat the grouped real `P2` channels are degenerate before additional matter/gauge anisotropy is turned on.

---

## 9. Linearized coupled bulk/interface skeleton

### 9.1 Stationary background

Take a stationary reference solution
\[
\psi(\mathbf X,t)=e^{-i\mu_0 t/\hbar}\psi_0(\mathbf X),
\qquad
A_M(\mathbf X,t)=A_{M0}(\mathbf X),
\qquad
R(\Omega,w,t)=R_0(w).
\]

### 9.2 Linearized BdG-type matter sector

Write
\[
\psi=e^{-i\mu_0 t/\hbar}(\psi_0+\delta\psi),
\qquad
\delta\rho=\psi_0^\ast\delta\psi+\psi_0\delta\psi^\ast.
\]

Schematically, the linearized matter sector is
\[
i\hbar \partial_t
\begin{bmatrix}
\delta\psi\\ \delta\psi^\ast
\end{bmatrix}
=
L_{\rm BdG}
\begin{bmatrix}
\delta\psi\\ \delta\psi^\ast
\end{bmatrix}
+
C_A[\delta A_M]
+
C_\eta[\eta],
\]
with wall source entering through
\[
\delta V_{\rm conf}
=
-\frac{V'_{\rm wall}(\Sigma_0/\ell_c)}{\ell_c}\,\eta.
\]

### 9.3 Linearized Maxwell sector

Write
\[
A_M=A_{M0}+\delta A_M.
\]

Then
\[
\boxed{
\partial_M\!\left(Z(w)\delta F^{MN}\right)
+\frac{1}{\xi}\partial^N(\partial\!\cdot\!\delta A)
=
\mu_0\,\delta J^N.
}
\]

The important ontology point is that the linearized mixed channels
\[
\delta A_w,\qquad \delta J^w,\qquad \delta F_{\mu w}
\]
all remain active.

### 9.4 Linearized geometry PDE

Varying the distributed wall action gives
\[
\boxed{
\mu_\eta \partial_t^2\eta
-
\partial_w(T_w\partial_w\eta)
-
T_\Omega \Delta_{S^2}\eta
+
K_\eta \eta
=
S_\eta^{(\psi)}[\delta\psi,\delta\psi^\ast]
+
S_\eta^{(A)}[\delta A_M]
+
f_{\rm ext}.
}
\]

Projecting onto harmonics gives separate modal equations for `l=0`, grouped real `l=2`, and higher multipoles.

---

## 10. Reduced conservative kernels produced by the PDE

The moving-throat PDE is expected to reduce, lane by lane, to the grouped real `P2` reduced bundle.

For each grouped lane \(A\in\{20,21,22\}\), define:

- wall/worldtube amplitude \(q_A\),
- stable BdG support modes \(X_{A\alpha}\) with frequencies \(\varpi_{A\alpha}\),
- localized brane-like gauge coordinates \(U_{A,r}\) with frequencies \(\Omega_{U,A,r}\),
- mixed `A_w/F_{\mu w}/J^w` coordinates \(W_{A,r}\) with frequencies \(\Omega_{W,A,r}\),
- internal mixed-sector couplings \(R_{A,r}\).

### 10.1 Conservative BdG moments

\[
B_{A,0}=\sum_\alpha \frac{c_{A\alpha}^2}{\varpi_{A\alpha}^2},
\qquad
B_{A,2}=\sum_\alpha \frac{c_{A\alpha}^2}{\varpi_{A\alpha}^4},
\qquad
B_{A,4}=\sum_\alpha \frac{c_{A\alpha}^2}{\varpi_{A\alpha}^6}.
\]

### 10.2 Conservative Maxwell/mixed moments

For each port \(r\),
\[
\Delta_{A,r}=\Omega_{U,A,r}^2\Omega_{W,A,r}^2-R_{A,r}^2,
\]
\[
S_{A,r}=\Omega_{U,A,r}^2+\Omega_{W,A,r}^2,
\]
\[
Q_{A,r}=g_{U,A,r}^2\Omega_{W,A,r}^2+2g_{U,A,r}g_{W,A,r}R_{A,r}+g_{W,A,r}^2\Omega_{U,A,r}^2,
\]
\[
G_{A,r}=g_{U,A,r}^2+g_{W,A,r}^2.
\]

Then
\[
Z_{A,0}^{(r)}=\frac{Q_{A,r}}{\Delta_{A,r}},
\]
\[
Z_{A,2}^{(r)}=\frac{Q_{A,r}S_{A,r}-G_{A,r}\Delta_{A,r}}{\Delta_{A,r}^2},
\]
\[
Z_{A,4}^{(r)}=\frac{Q_{A,r}(S_{A,r}^2-\Delta_{A,r})-S_{A,r}G_{A,r}\Delta_{A,r}}{\Delta_{A,r}^3}.
\]

Summing over ports gives \(Z_{A,n}\).

### 10.3 Outgoing-transfer moments

Define
\[
P_{A,r}=\Omega_{U,A,r}^2 g_{W,A,r}+R_{A,r}g_{U,A,r}.
\]

Then the outgoing-transfer moments begin with
\[
N_{A,0}^{(r)}=\frac{P_{A,r}^2}{\Delta_{A,r}^2},
\]
\[
N_{A,2}^{(r)}=\frac{2P_{A,r}\big(P_{A,r}S_{A,r}-\Delta_{A,r}g_{W,A,r}\big)}{\Delta_{A,r}^3},
\]
\[
N_{A,4}^{(r)}
=
\frac{
\Delta_{A,r}^2 g_{W,A,r}^2
-2\Delta_{A,r}P_{A,r}^2
-4\Delta_{A,r}P_{A,r}S_{A,r}g_{W,A,r}
+3P_{A,r}^2 S_{A,r}^2
}{\Delta_{A,r}^4}.
\]

Summing over ports gives \(N_{A,n}\).

### 10.4 Conservative grouped-lane operator

The full conservative grouped-lane operator is
\[
D_A^{\rm(cons)}(\omega)=D_{A,0}+D_{A,2}\omega^2+D_{A,4}\omega^4+O(\omega^6),
\]
with
\[
\boxed{
D_{A,0}=K_A-B_{A,0}-Z_{A,0},
}
\]
\[
\boxed{
D_{A,2}=-\big(M_A+B_{A,2}+Z_{A,2}\big),
}
\]
\[
\boxed{
D_{A,4}=-\big(B_{A,4}+Z_{A,4}\big).
}
\]

---

## 11. Grouped real `P2` normalized response and isotropy

Define the normalized grouped response
\[
Y_A(\omega)=\frac{D_{A,0}}{D_A^{\rm(cons)}(\omega)}
=
1+u_2^{(A)}\omega^2+u_4^{(A)}\omega^4+O(\omega^6).
\]

Then
\[
\boxed{
u_2^{(A)}=-\frac{D_{A,2}}{D_{A,0}},
}
\qquad
\boxed{
u_4^{(A)}=\frac{D_{A,2}^2-D_{A,0}D_{A,4}}{D_{A,0}^2}.
}
\]

For any grouped triple \(x_A\), define the weighted trace/anomaly variables
\[
x_{\rm bar}=\frac{x_{20}+2x_{21}+2x_{22}}{5},
\]
\[
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\qquad
b_x=\frac{x_{21}-x_{22}}{2}.
\]

Applied to \(u_2^{(A)}\) and \(u_4^{(A)}\), the isotropy gate is
\[
\boxed{
a_2=b_2=a_4=b_4=0.
}
\]

On the isotropic branch the three grouped lanes collapse to common coefficients
\[
D_{20,n}=D_{21,n}=D_{22,n}=D_n,
\qquad
N_{20,n}=N_{21,n}=N_{22,n}=N_n.
\]

---

## 12. Outgoing `l=2` bridge and the remaining theorem gap

### 12.1 Compact outgoing `l=2` fingerprint

The normalized outgoing compact `l=2` branch has the low-frequency fingerprint
\[
\boxed{
\widehat Y_2^{\rm(out)}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\,\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6).
}
\]

### 12.2 Outgoing prefactor

The grouped-lane outgoing prefactor is
\[
P_0=\frac{N_0}{D_0},
\]
and more generally
\[
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\]
\[
P_4=
\frac{
D_0^2N_4
-2D_0(D_2N_2+D_4N_0)
+3D_2^2N_0
}{D_0^3}.
\]

On the isotropic branch, the leading odd coefficient is controlled only by the static outgoing prefactor \(P_0\).

### 12.3 Universal quadrupole-normalization condition

The invariant normalization condition isolated by the 2.5PN bridge is
\[
\boxed{
m_{\hat 0}^{\,2} P_0
=
\frac{54\,G\,c_s^5}{5\,a^5 c^5}.
}
\]

Equivalently, at leading point-particle order on the natural source-map branch \(m_{\hat 0}=1+O(a^2/r^2)\),
\[
P_0=\frac{54\,G\,c_s^5}{5\,a^5 c^5}.
\]

This is the **single sharp normalization target** still left for the completed moving-throat PDE.

---

## 13. What is already solved, and what is still open

### 13.1 Solved at the reduced-theory level

- exact 4D GNLS + localized Maxwell parent equations,
- exact projection identities,
- exact level-set throat lift,
- exact finite-throat D/N support branch,
- exact distributed wall reduction,
- exact reduced grouped-lane conservative bundle formulas,
- exact grouped-response conversion formulas,
- exact statement of the outgoing quadrupole normalization target.

### 13.2 Still open

- the fully solved nonlinear moving-throat PDE,
- the true branch data \((K_A,M_A,B_{A,n},Z_{A,n},N_{A,n})\) of the completed solution,
- the exact selection of the isotropic/passive/outgoing quadrupole branch on that solution,
- and the final hit on
  \[
  m_{\hat 0}^{\,2}P_0=\frac{54Gc_s^5}{5a^5c^5}.
  \]

---

## 14. Minimal source anchors

This engine document was distilled from the exact/reduced stack anchored by

- `4d_summary.md`
- `4d_em_fields_summary.md`
- `4d_plasma_summary.md`
- `4d_1pn_bridge_summary.md`
- `moving_throat_pde_full.md`
- `moving_throat_pde_stage147_microscopic_log_channels.md`
- `moving_throat_pde_stage148_exact_branch_drifts.md`
- `moving_throat_pde_stage168_microscopic_monomials.md`
- `moving_throat_pde_stage169_similarity_orbit_closure.md`
- `moving_throat_pde_stage170_orbit_quotient_closure.md`

The intended use is: a fresh session should be able to start from **this document alone**, then only drill into the stage files if it wants the longer derivation path.
# Moving-Throat Translation Dictionary
## From 4D fields to reduced macroscopic variables, microscopic variables, and the three monomial invariants

## 0. Why this document exists

This document is the **ledger** handoff. It is the compact dictionary that turns the raw 4D/moving-throat PDE variables into the reduced response variables actually used in the later derivations.

It is meant to let a fresh session answer questions of the form:

- what exactly is a macroscopic variable here?
- what are the microscopic kernel variables?
- how do the grouped `P2` response coefficients arise?
- what are the direct coherent-branch defect coordinates?
- and what are the **three exact monomial invariants** that now carry the final closure theorem?

This document is structured from coarse to fine:

1. core 4D variables,
2. moving-throat and grouped-`P2` macroscopic/reduced variables,
3. coherent-branch microscopic variables,
4. microscopic slippages,
5. branch-adapted defect coordinates,
6. the three monomial invariants,
7. the similarity orbit / quotient theorem.

The intended use is simple: if a new session can read **this ledger + the PDE engine document**, it should be able to reconstruct the current theorem status without needing the conversational history.

---

## 1. Notation firewall and reading rules

### 1.1 Exact vs reduced vs open

- **Exact**: follows directly from the declared action, exact definitions, or exact algebra.
- **Controlled reduction**: follows only after a stated ansatz or reduction.
- **Protocol closure**: fixed only inside the declared response hierarchy.
- **Open**: still depends on the completed moving-throat PDE.

### 1.2 Notation firewall

The following notational separations are non-negotiable.

1. Electric charge is carried by
   \[
   \eta_Q,\ q_\star,\ q_{\rm eff},
   \]
   not by circulation.
2. The historical gravity-side bare `q=1` is really
   \[
   \kappa_\rho=1,
   \]
   not electric charge.
3. Grouped labels `20/21/22` refer to grouped real `P2` lanes, not spacetime indices.
4. The mixed channels
   \[
   A_w,\ J^w,\ F_{\mu w},\ E_w,\ C_a
   \]
   are suppressed only in the strict far-field brane reduction. They remain microscopic degrees of freedom.

---

## 2. Core 4D variables

### 2.1 Coordinates and fields

\[
x^M=(t,x,y,z,w),\qquad
\mathbf X=(x,y,z,w),\qquad
\mathbf x=(x,y,z).
\]

\[
\psi(\mathbf X,t),\qquad \rho=|\psi|^2,\qquad
A_M=(A_0,A_i),\qquad F_{MN}=\partial_MA_N-\partial_NA_M.
\]

### 2.2 Charge variables

\[
\eta_Q=\pm 1,\qquad
q_\star=\eta_Q e_\star,\qquad
e_\star>0.
\]

After zero-mode canonical normalization,
\[
q_{\rm eff}=\frac{q_\star}{\sqrt{Z_{\rm int}}},
\qquad
e_{\rm eff}=\frac{e_\star}{\sqrt{Z_{\rm int}}},
\qquad
Z_{\rm int}=\int Z(w)\,dw.
\]

### 2.3 Old collective throat variables

\[
a(t),\qquad L(t).
\]

These are **collective moments** of the moving-throat shape field, not the fundamental geometry variables.

---

## 3. Moving-throat geometry variables

### 3.1 Level-set and shape-field representation

\[
\Sigma(\mathbf X,t)=r-R(\Omega,w,t),
\qquad
r=\sqrt{x^2+y^2+z^2},
\qquad
\Omega=\mathbf x/r\in S^2.
\]

The throat surface is \(\Sigma=0\).

The stationary reference throat is
\[
\Sigma_0(\mathbf X)=r-R_0(w).
\]

### 3.2 Wall displacement

\[
R(\Omega,w,t)=R_0(w)+\eta(\Omega,w,t).
\]

### 3.3 Harmonic decomposition

\[
\eta(\Omega,w,t)
=
\eta_0(w,t)Y_{00}(\Omega)
+\sum_{m\in P_2({\rm real})}q_{2m}(w,t)\,Y_{2m}^{\rm real}(\Omega)
+\eta_{\ge 3}.
\]

The grouped real `P2` set is
\[
\{20,\ 21c,\ 21s,\ 22c,\ 22s\}.
\]

### 3.4 Monopole normalization bridge

With
\[
Y_{00}=\frac{1}{2\sqrt\pi},
\]
the physical mouth-average shift \(\delta a\) and the normalized monopole coefficient satisfy
\[
\boxed{
q_{00}(0,t)=2\sqrt\pi\,\delta a(t).
}
\]

---

## 4. Brane/macroscopic variables

### 4.1 Projected brane observables

\[
\rho_{\rm brane}=\int W(w)\rho(\mathbf x,w,t)\,dw,
\qquad
\mathbf j_{\rm brane}=\int W(w)\mathbf j_{xyz}(\mathbf x,w,t)\,dw,
\]
\[
\mathbf v_{\rm brane}=\mathbf j_{\rm brane}/\rho_{\rm brane}.
\]

### 4.2 Leakage

\[
S_{\rm leak}
=
-\left[Wj^w\right]_{-\infty}^{+\infty}
+\int W'(w)j^w\,dw.
\]

### 4.3 Brane velocity potential and Poisson hook

\[
\mathbf v_{\rm brane}=\nabla_3\varphi+\mathbf v_T,\qquad \nabla_3\cdot\mathbf v_T=0.
\]

In the quasi-static regime this gives the Poisson hook for \(\varphi\).

---

## 5. Family-1 / core–mouth macroscopic variables

These variables belong to the earlier core–mouth / compensated-branch part of the moving-throat program. They are not the final monomial invariants, but they remain part of the current dictionary.

### 5.1 Support/mouth source and stiffness variables

\[
K_s,\qquad K_q,\qquad \lambda,\qquad g_s,\qquad g_q.
\]

On the explicit throat-core branch,
\[
K_s=4\pi a^2\!\left(\frac{H_w\ell}{3}+\frac{\hbar^2}{15m_\psi\rho_w\,\ell}\right),
\]
\[
K_q=\frac{\mathcal Z_q}{\mu_0}\frac{\pi^2 c_s^2}{4L_W^2},
\]
\[
\lambda=-q_\star v_{w0}\mathcal I_{sq},
\]
\[
g_s=\mathcal T_m\frac{4\pi a^2\ell}{3},
\qquad
g_q=\frac{\mathcal Z_q}{\mu_0}\frac{\pi}{\sqrt2\,L_W^{3/2}}.
\]

### 5.2 Dimensionless Family-1 ratios

\[
\boxed{
\mathfrak r=\frac{\lambda}{\sqrt{K_sK_q}},
\qquad
\mathfrak g=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}}.
}
\]

### 5.3 Traction and source moments

The canonical mouth-traction variable is
\[
\boxed{
\Sigma_0=\frac{20}{9}\,\widehat T_m^2.
}
\]

The mixed loading ratio is
\[
\boxed{
\mathcal R=\frac{(\mathfrak g-\mathfrak r)^2}{1+\mathfrak r^2}.
}
\]

The support-shape functional is
\[
\mathcal S[\Sigma]=\int_0^1 K_q(x)\,\Sigma(x)\,dx,
\]
with the reference quadrupole kernel used in the fixed-point audits
\[
K_q(x)=\frac{\cosh\!\left(\frac{\pi}{2}(1-x)\right)}{\cosh(\pi/2)}.
\]

The mouth slope variable is
\[
\boxed{
\Pi=\Sigma_0\,[1-\mathcal R\,\mathcal S].
}
\]

The mouth source moments are
\[
\boxed{
M_s=\Sigma_0,
\qquad
M_q=-\Sigma_0\,\mathcal R.
}
\]

### 5.4 Useful lower-branch exact identities

The parent compensation condition is
\[
1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2.
\]

On the lower compensated branch,
\[
\mathfrak g_-(\mathfrak r)=\mathfrak r-\frac12\sqrt{1+\mathfrak r^2}.
\]

The first off-family normal coordinate is
\[
\delta_\perp
=
\delta\mathfrak g
-
\mathfrak g'_-(\mathfrak r_\ast)\,\delta\mathfrak r.
\]

This later collapses to explicit logarithmic microscopic imbalance channels.

---

## 6. Microscopic variables behind the core–mouth branch

### 6.1 Throat-core microscopic variables

\[
a,\qquad L_W,\qquad \rho_w,\qquad c_{s,w},\qquad c_s,\qquad \mathcal Z_q,\qquad \mathcal T_m,\qquad v_{w0}.
\]

### 6.2 Healing-lock shell variables

\[
\ell=\frac{\hbar}{2m_\psi c_{s,w}},
\qquad
K_s=\frac{3\pi a^2\hbar^2}{5m_\psi\rho_w\,\ell}
\]
on the carried healing-locked shell branch.

### 6.3 Exact lower-branch drift laws

On the exact lower compensated branch,
\[
\delta\ln L_W=\delta\ln a,
\]
\[
\delta\ln v_{w0}
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q}{\rho_w}\right)
+\frac32\,\delta\ln c_{s,w}
+\delta\ln c_s
-\frac52\,\delta\ln a,
\]
\[
\delta\ln \mathcal T_m
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q}{\rho_w}\right)
+\frac32\,\delta\ln c_{s,w}
-\delta\ln c_s
-\frac32\,\delta\ln a.
\]

These are exact reduced lower-branch transport identities.

---

## 7. Grouped real `P2` macroscopic response variables

For each grouped lane \(A\in\{20,21,22\}\), define the conservative operator
\[
D_A^{\rm(cons)}(\omega)=D_{A,0}+D_{A,2}\omega^2+D_{A,4}\omega^4+O(\omega^6).
\]

### 7.1 Normalized grouped response moments

\[
Y_A(\omega)=\frac{D_{A,0}}{D_A^{\rm(cons)}(\omega)}
=
1+u_2^{(A)}\omega^2+u_4^{(A)}\omega^4+O(\omega^6).
\]

Then
\[
u_2^{(A)}=-\frac{D_{A,2}}{D_{A,0}},
\qquad
u_4^{(A)}=\frac{D_{A,2}^2-D_{A,0}D_{A,4}}{D_{A,0}^2}.
\]

### 7.2 Grouped weighted trace/anomaly variables

For any grouped triple \(x_A\),
\[
x_{\rm bar}=\frac{x_{20}+2x_{21}+2x_{22}}{5},
\]
\[
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\qquad
b_x=\frac{x_{21}-x_{22}}{2}.
\]

Applied to the grouped response moments,
\[
u_{\rm bar,2},\ a_2,\ b_2,\qquad
u_{\rm bar,4},\ a_4,\ b_4
\]
are the isotropic trace and the two anisotropy defects.

The grouped isotropy gate is
\[
a_2=b_2=a_4=b_4=0.
\]

### 7.3 Outgoing prefactor data

On the isotropic branch,
\[
P_0=\frac{N_0}{D_0},
\]
\[
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\]
\[
P_4=
\frac{
D_0^2N_4
-2D_0(D_2N_2+D_4N_0)
+3D_2^2N_0
}{D_0^3}.
\]

The universal 2.5PN/4PN normalization target uses only \(P_0\) at leading order.

---

## 8. Coherent local-kernel branch variables

These are the central macroscopic variables of the later coherent-branch and invariant program.

### 8.1 Effective stiffnesses

\[
K_{U1}=K_U(1+\delta_U),
\]
\[
K_\eta^{(\mathrm{eff})}=K_\eta+6T_\Omega,
\]
\[
K_W^{(\mathrm{eff})}=K_W+\frac{\pi^2 T_W}{4L^2},
\]
\[
K_\phi^{(\mathrm{eff})}=K_\phi+\frac{\pi^2 T_\phi}{4L^2}.
\]

### 8.2 Dimensionless coherent ratios

\[
\boxed{
\epsilon_\eta=\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}},
}
\]
\[
\boxed{
\epsilon_W=\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}},
}
\]
\[
\boxed{
Z_W=\frac{\lambda_W^2}{K_\eta^{(\mathrm{eff})}K_W^{(\mathrm{eff})}},
}
\]
\[
\boxed{
\delta_U=\frac{\pi^2 T_U}{L^2 K_U},
}
\]
\[
\boxed{
\chi_0=\frac{\gamma c_{\eta U}}{K_U},
}
\]
\[
\boxed{
\zeta=\frac{\lambda_\phi^2 K_W^{(\mathrm{eff})}}{\lambda_W^2 K_\phi^{(\mathrm{eff})}},
}
\]
\[
\boxed{
\Lambda=\frac{27\pi^2 G c_s^5 K_W^{(\mathrm{eff})}}{20 a^5 c^5 \mu_W}.
}
\]

Useful coherent identities:
\[
\rho_0=\sigma_0=\chi_0,
\qquad
\epsilon_\phi=\zeta\epsilon_W,
\qquad
Z_\phi=\zeta Z_W.
\]

### 8.3 Split blocking ratio and tracking factor

The split mixed blocking ratio is
\[
\boxed{
\epsilon
=
\epsilon_W\!\left[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right].
}
\]

The exact tracking factor is
\[
\boxed{
R_{\rm tr}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0}
=
\frac{1+\chi_0+\delta_U}{(1+\chi_0)(1+\delta_U)}.
}
\]

On the constructive coherent branch,
\[
\frac{1}{1+\delta_U}<R_{\rm tr}<1.
\]

### 8.4 Mixed, support, and total baselines

\[
\boxed{
M_{\rm mix}
=
\frac{8 Z_W (1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\epsilon)}.
}
\]

Using \(Z_\phi=\zeta Z_W\) and \(\epsilon_\phi^{\rm(split)}=\zeta\epsilon\),
\[
\boxed{
M_{\rm supp}
=
\frac{8\zeta Z_W(1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\zeta\epsilon)}.
}
\]

The total baseline is
\[
\boxed{
M_{\rm tr}=M_{\rm mix}+M_{\rm supp}=M_{\rm mix}S(\zeta;\epsilon),
}
\]
with support-enhancement factor
\[
\boxed{
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon}.
}
\]

### 8.5 Transfer shape and normalization demand ratio

The coherent transfer shape is
\[
\boxed{
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2}.
}
\]

The exact selected-branch demand ratio is
\[
\boxed{
R_{\rm target}
=
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2}.
}
\]

Equivalent identity:
\[
\boxed{
\mathcal T^2
=
\frac{27\pi^2 G c_s^5}{20 a^5 c^5}\,
\frac{1-\epsilon_\eta}{R_{\rm target}}.
}
\]

A key coherent-branch fact is that \(R_{\rm target}\) is independent of the support ratio \(\zeta\).

---

## 9. Microscopic variables and their grouped weak-axisymmetric drifts

The key positive microscopic state vector is
\[
x=
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})},\ K_W^{(\mathrm{eff})},\ \mu_W,\ T_U).
\]

Its grouped weak-axisymmetric log-drift vector is
\[
\delta\mathbf x
=
\begin{pmatrix}
\lambda_1\\
c_1\\
\gamma_1\\
\kappa_U\\
\kappa_\eta\\
\kappa_W\\
\mu_1\\
\tau_1
\end{pmatrix}
=
\begin{pmatrix}
\delta\ln\lambda_W\\
\delta\ln c_{\eta U}\\
\delta\ln\gamma\\
\delta\ln K_U\\
\delta\ln K_\eta^{(\mathrm{eff})}\\
\delta\ln K_W^{(\mathrm{eff})}\\
\delta\ln\mu_W\\
\delta\ln T_U
\end{pmatrix}_{\rm grp}.
\]

---

## 10. Microscopic slippage variables

The exact coherent-kernel slippages are

\[
\boxed{
\Sigma_Z
=
2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W
=
\delta\ln\!\left(\frac{Z_W}{\Omega_W^2}\right),
}
\]
\[
\boxed{
\Sigma_\chi
=
\gamma_1+c_1-\kappa_U
=
\delta\ln\chi_0,
}
\]
\[
\boxed{
\Sigma_\eta
=
2c_1-\kappa_U-\kappa_\eta
=
\delta\ln\epsilon_\eta,
}
\]
\[
\boxed{
\Sigma_\epsilon
=
2\gamma_1+2\lambda_1-\kappa_U-\kappa_W
=
\delta\ln\epsilon_W,
}
\]
\[
\boxed{
\Sigma_\delta
=
\tau_1-\kappa_U
=
\delta\ln\delta_U.
}
\]

The tracking combination is
\[
\boxed{
\Sigma_{\rm tr}
=
(1+\chi_0)\Sigma_\delta + (1+\delta_U)\Sigma_\chi.
}
\]

The genuine nontracking transfer-shape slippage is
\[
\boxed{
\Sigma_{\rm nt}
=
\Sigma_Z
+\frac{2\epsilon_W}{1-\epsilon}\frac{11+9\delta_U}{11(1+\delta_U)}\,\Sigma_\epsilon
-\left[
\frac{2\chi_0}{1+\delta_U}
+\frac{4\epsilon_W\delta_U}{11(1-\epsilon)(1+\delta_U)^2}
\right]\Sigma_\delta.
}
\]

---

## 11. Observable defect variables

The first grouped weak-axisymmetric observables are

- \(\Theta_1\): tracking-factor drift,
- \(\Xi_1\): grouped transfer-shape drift,
- \(\mathcal R_1\): selected-branch demand-ratio drift.

The exact triangular normal form is

\[
\boxed{
\Theta_1
=
-\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}\,\Sigma_{\rm tr},
}
\]
\[
\boxed{
\Xi_1
=
\frac{2\chi_0}{(1+\chi_0)(1+\delta_U)}\,\Sigma_{\rm tr}
+\Sigma_{\rm nt},
}
\]
\[
\boxed{
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_\eta}{1-\epsilon_\eta}\,\Sigma_\eta.
}
\]

So:

- \(\Sigma_{\rm tr}\) carries tracking failure,
- \(\Sigma_{\rm nt}\) carries nontracking transfer-shape failure,
- \(\Sigma_\eta\) carries dressing failure.

---

## 12. Intermediate exact branch composites

Stage 167 packaged the same three directions into exact branch composites.

Define
\[
B_*=\frac{2(1+\chi_{0,*}+\delta_{U,*})}{\delta_{U,*}},
\qquad
C_*=\frac{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}{\chi_{0,*}\delta_{U,*}}.
\]

Then
\[
\mathfrak T_* = R_{\rm tr}^{-C_*},
\qquad
\mathfrak N_* = \mathcal T^2 R_{\rm tr}^{B_*},
\qquad
\mathfrak D = \epsilon_\eta,
\]
with
\[
\delta\ln\mathfrak T_*=\Sigma_{\rm tr},
\qquad
\delta\ln\mathfrak N_*=\Sigma_{\rm nt},
\qquad
\delta\ln\mathfrak D=\Sigma_\eta.
\]

These are useful, but the final closure is even sharper in the direct microscopic monomials below.

---

## 13. The three final monomial invariants

These are the final reduced invariants that matter.

### 13.1 Tracking monomial

\[
\boxed{
\mathfrak C_{{\rm tr},*}
=
\chi_0^{\,1+\delta_{U,*}}\,
\delta_U^{\,1+\chi_{0,*}}.
}
\]

It satisfies
\[
\boxed{
\delta\ln \mathfrak C_{{\rm tr},*}=\Sigma_{\rm tr}.
}
\]

### 13.2 Nontracking monomial

Define
\[
\boxed{
E_*
=
\frac{2\epsilon_{W,*}}{1-\epsilon_*}\,
\frac{11+9\delta_{U,*}}{11(1+\delta_{U,*})},
}
\]
\[
\boxed{
F_*
=
\frac{2\chi_{0,*}}{1+\delta_{U,*}}
+
\frac{4\epsilon_{W,*}\delta_{U,*}}{11(1-\epsilon_*)(1+\delta_{U,*})^2}.
}
\]

Then the nontracking monomial is
\[
\boxed{
\mathfrak C_{{\rm nt},*}
=
\frac{Z_W}{\Omega_W^2}\,
\epsilon_W^{E_*}\,
\delta_U^{-F_*}.
}
\]

It satisfies
\[
\boxed{
\delta\ln \mathfrak C_{{\rm nt},*}=\Sigma_{\rm nt}.
}
\]

### 13.3 Dressing invariant

The third invariant is simply
\[
\boxed{
\epsilon_\eta=\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}},
}
\]
with
\[
\boxed{
\delta\ln\epsilon_\eta=\Sigma_\eta.
}
\]

---

## 14. Similarity orbit and exact quotient theorem

### 14.1 Invariant map

Define the invariant map
\[
\boxed{
\mathcal I(x)=\bigl(\mathfrak C_{{\rm tr},*}(x),\ \mathfrak C_{{\rm nt},*}(x),\ \epsilon_\eta(x)\bigr).
}
\]

### 14.2 Exact monomial-drift matrix

The exact finite/infinite log-drift map is
\[
\begin{pmatrix}
\delta\ln\mathfrak C_{{\rm tr},*}\\
\delta\ln\mathfrak C_{{\rm nt},*}\\
\delta\ln\epsilon_\eta
\end{pmatrix}
=
M_*\,
\delta\mathbf x,
\]
with
\[
\boxed{
M_*=
\begin{pmatrix}
0 & 1+\delta_{U,*} & 1+\delta_{U,*} & -(2+\chi_{0,*}+\delta_{U,*}) & 0 & 0 & 0 & 1+\chi_{0,*}\\
2(1+E_*) & 0 & 2E_* & F_*-E_* & -1 & -(2+E_*) & 1 & -F_*\\
0 & 2 & 0 & -1 & -1 & 0 & 0 & 0
\end{pmatrix}.
}
\]

A useful minor gives
\[
\det M_*^{(\tau,\kappa_\eta,\mu_1)}=1+\chi_{0,*}>0,
\]
so the map has rank \(3\) and kernel dimension \(5\).

### 14.3 Exact five-parameter similarity orbit

Choose free co-scalings for
\[
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_W^{(\mathrm{eff})})
\]
and determine the remaining three by monomial preservation:
\[
K_\eta^{(\mathrm{eff})}\mapsto e^{\,2C-U}K_\eta^{(\mathrm{eff})},
\]
\[
T_U\mapsto
e^{\,U-\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Gamma+C-U)}T_U,
\]
\[
\mu_W\mapsto
e^{\,M(\Lambda,C,\Gamma,U,W)}\mu_W,
\]
where
\[
M(\Lambda,C,\Gamma,U,W)
=
2C-U+2W-2\Lambda
-
E_*(2\Gamma+2\Lambda-U-W)
-
F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Gamma+C-U).
\]

This defines the exact five-parameter similarity family \(\mathcal G_*\).

### 14.4 Exact finite quotient theorem

On the positive-coupling state space
\[
\mathcal M_+
=
\{(\lambda_W,c_{\eta U},\gamma,K_U,K_\eta^{(\mathrm{eff})},K_W^{(\mathrm{eff})},\mu_W,T_U)>0\},
\]
the fibres of \(\mathcal I\) are exactly the \(\mathcal G_*\)-orbits:
\[
\boxed{
\mathcal I(\widetilde x)=\mathcal I(x)
\iff
\widetilde x\in \mathcal G_*\cdot x.
}
\]

So
\[
\boxed{
\mathcal M_+/\mathcal G_*
\cong (\mathbb R_{>0})^3
\quad\text{with quotient coordinates}\quad
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta).
}
\]

This is the cleanest final reduced closure.

---

## 15. Final theorem ledger

### 15.1 Infinitesimal weak-axisymmetric closure

\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\ln\mathfrak C_{{\rm tr},*}
=
\delta\ln\mathfrak C_{{\rm nt},*}
=
\delta\ln\epsilon_\eta
=0.
}
\]

Equivalently,
\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\mathbf x\in T_{\rm id}\mathcal G_*.
}
\]

### 15.2 Finite closure

\[
\boxed{
\text{The coherent weak-axisymmetric defect vanishes exactly when the actual microscopic branch stays on one }
\mathcal G_*\text{-orbit,}
}
\]
that is, exactly when the three quotient coordinates
\[
(\mathfrak C_{{\rm tr},*},\ \mathfrak C_{{\rm nt},*},\ \epsilon_\eta)
\]
are preserved.

### 15.3 What is still open

The remaining theorem gap is now purely **branch-selective**:

> Does the actual completed moving-throat branch preserve the three exact quotient coordinates?

All algebraic compression is finished. The only remaining unknown is the true branch dynamics of the full PDE.

---

## 16. Minimal source anchors

This dictionary was distilled from the final moving-throat stages anchored by

- `moving_throat_pde_full.md`
- `moving_throat_pde_stage147_microscopic_log_channels.md`
- `moving_throat_pde_stage148_exact_branch_drifts.md`
- `moving_throat_pde_stage165_microscopic_coherent_slippage.md`
- `moving_throat_pde_stage166_triangular_normal_form.md`
- `moving_throat_pde_stage167_branch_invariant_coordinates.md`
- `moving_throat_pde_stage168_microscopic_monomials.md`
- `moving_throat_pde_stage169_similarity_orbit_closure.md`
- `moving_throat_pde_stage170_orbit_quotient_closure.md`

The intended rule is:

- use the **PDE engine** document for the equations and branch construction,
- use **this dictionary** for the reduced variable ledger and the final invariant theorem.
