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

From this point onward, the modal equations are written in a densitized
convention: after integrating over the reference sphere, the remaining
background surface weight is absorbed into the effective one-dimensional
coefficients `mu_eta, T_w, T_Omega, K_eta` and the wall amplitudes. Equivalently,
the axial operator is written with respect to the flat measure `dw`. If one kept
an explicit varying weight `g(w) = sqrt(gamma_0)`, the Euler-Lagrange equation
would instead contain the extra first-derivative term `T_w g' q_w / g`.

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
