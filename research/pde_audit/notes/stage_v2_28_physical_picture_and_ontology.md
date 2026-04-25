# Stage V2-28 - Physical Picture and Ontology Checklist

## Purpose

This note records the physical picture that the papers should state explicitly.
It exists because the model can be misread if the defect is described only by
reduced coefficients or mouth response variables.

The most important point:

```text
The particle is modeled as a finite brane-bulk throat defect: a puncture/open
conduit through the brane into the bulk, not a mere surface depression,
indentation, dimple, or capped pocket.
```

The code/audit variables `D0`, `C`, `P0`, `N2`, and `N4` are low-order
projected response readouts of this object.  They are not the full physical
description of the object.

---

## 1. Core ontology

The model treats a particle as a brane-bulk throat defect.

The brane sees a localized mouth field.  The interior carries throat/cavity
structure extending into the bulk.  Charge, mass, support moments, and outgoing
moments are reduction-layer quantities, not primitive labels attached to a
point particle.

The papers should say this early:

```text
The defect is a finite-radius opening of the brane into a bulk throat.  The
mouth is the brane-side cross-section of that throat.  The interior throat
supports bulk/superfluid, wall, BdG, Maxwell/mixed, and outgoing-port degrees
of freedom.  The reduced particle observables are extracted from this branch.
```

Do not describe the object as:

- a depression in the brane;
- a surface dimple;
- a shallow pocket;
- a closed bubble attached to the brane;
- a hard-capped tube, except when explicitly discussing obsolete toy models or
  negative controls.

---

## 2. Geometry

The effective moving-throat geometry is currently represented as a level-set or
shape-field lift:

```text
Sigma(X,t) = r - R(Omega,w,t),
```

where `w` is the transverse/bulk direction, `Omega` is the angular direction on
the brane-side sphere, and the finite throat surface is

```text
Sigma(X,t) = 0.
```

The convention carried by the notes is:

```text
exterior region: Sigma > 0,
interior/support region: Sigma < 0.
```

This geometry is an effective closure unless/until the parent throat action
`S_Sigma` is promoted and frozen.

The reference stationary branch has:

```text
mouth at w = 0,
finite depth 0 <= w <= L,
mouth radius R(0) = a,
open finite exit R(L) > 0.
```

The open-endpoint correction is load-bearing:

```text
R(0) = a,
R(L) > 0.
```

This is an open finite-radius conduit.  The old hard cap

```text
R(L) = 0
```

is an obsolete toy idealization unless explicitly declared as a negative
control or simplification.

---

## 3. Mouth versus interior

The mouth is the brane-side cross-section of the throat.  It is not the whole
defect.

Important distinctions:

- mouth geometry: radius, ellipticity, headless `P22` axis, area;
- interior geometry: axial throat profile, depth, support/cavity structure;
- open exit: finite-radius outlet into the bulk/reservoir side;
- projected branch observables: reduced coefficients extracted after solving
  or freezing the branch.

The papers should avoid language that collapses these into one scalar "size" or
"depth" unless the reduction being used has explicitly done that.

In particular, the same mouth radius can coexist with different internal
support, throughput, outgoing, and orbit-lock data.  A rigid-mouth reading is a
statement about the brane entrance geometry; it is not automatically a theorem
about every internal transfer or orbit-lock variable.

---

## 4. Open system and superfluid intake/output

The throat is an open-system object in the physical picture.

There are several flux-like quantities that must not be identified by fiat:

1. internal radial throughput entering the self-flow energy ledger;
2. projected brane-bulk exchange source built from transverse bulk current;
3. brane-side mouth-output flux measured through the mouth surface;
4. leakage/work terms through the open finite-radius junction;
5. outgoing/radiative port flux.

The superfluid material state is also not determined by these flux labels
alone.  A completed branch must separately specify or derive the density field
`rho`, EOS/internal energy, local sound speed `c_s(rho)`, and any effective
light-speed relation if the model makes the light cone density dependent.  See
`stage_v2_29_superfluid_material_closure_gap.md`.

The exact D/N trapped support mode is not itself a trans-brane current injector.
On that branch the support field can vanish at the mouth while its mouth
gradient is nonzero.  The first nontrivial mouth datum is therefore a quadratic
normal stress, not the mouth value.

The carried hammer-stress theorem gives a cycle-averaged normal mouth stress
after support-action normalization:

```text
T_nn_bar = pi*hbar*c_s/(2*L^2).
```

But the constitutive response coefficient that turns mouth stress into actual
mouth flux is not yet fully derived from the completed PDE.

A nonzero DC mouth source in the current notes appears only on an explicitly
open/radiative branch.  Closed conservative standing-wave support does not by
itself produce the required DC mouth flux.

---

## 5. Charge and circulation

The charge sign is tied to puncture orientation in the carried ontology:

```text
eta_Q = +/- 1.
```

The microscopic charge branch is

```text
q_* = eta_Q e_*.
```

The observable brane charge is thickness/localization dressed.

The fluxoid/circulation law is a topological law for loops surrounding the
mouth.  It quantizes the tangential vortical class.  It does not determine the
radial feed amplitude or the mouth-output strength by itself.

The papers should keep these separate:

- puncture orientation;
- microscopic charge branch;
- observable dressed charge;
- circulation/fluxoid sector;
- radial feed or throughput amplitude.

For the electromagnetic paper-facing status and claim boundary, see
`stage_v2_30_electromagnetic_ontology_and_status.md`.

---

## 6. Finite-mouth shape physics

Finite throat size is not a cosmetic correction.  It changes the first
shape-sensitive response problem.

For a finite throat, the first shape-sensitive external load is the partner
field Hessian across the mouth, not the raw scalar potential at a point.  This
load decomposes into:

```text
P0 trace/scalar channel,
P2 traceless quadrupole channel.
```

Under centering and area preservation, the first driven non-axisymmetric mouth
deformation is a real headless `P22` mouth ellipse.

The exact constant-area ellipse has semiaxes

```text
R1 = R_th exp(epsilon_ell),
R2 = R_th exp(-epsilon_ell),
R1*R2 = R_th^2,
Area = pi R_th^2.
```

Its first-order boundary perturbation is a real `P22` pair, not a dipole or
free monopole:

```text
delta r_m/R_th = epsilon_ell cos(2(phi - theta_mouth)).
```

The mouth quadrupole tensor is headless.  The axis satisfies

```text
theta_mouth ~ theta_mouth + pi.
```

This `P22` sector is essential for the atom/same-charge bridge notes, but it is
still reduced/conditional when continued to isolated same-charge structure.

---

## 7. Parent status and closure status

The physical picture must preserve the status firewall.

The strict parent action currently contains the GNLS/matter and localized
Maxwell sectors.  The moving throat is parent-level as a confinement-coupling
argument, but it is not yet parent-level as an autonomous dynamical field unless
`S_Sigma` is added to the total action.

Therefore the current wall/throat PDE should be described as:

```text
effective wall/throat closure
```

unless a paper explicitly derives and freezes the promoted parent action.

The reduced wall action and distributed wall equations are useful and
mathematically consistent inside the closure.  They should not be advertised as
already deriving the full parent moving-throat PDE.

---

## 8. What the response coefficients do not describe

The audit readouts

```text
D0, C, P0, N2, N4
```

are not physical ontology variables.  They are low-order projected response
moments extracted from a branch.

They do not directly describe:

- the throat surface `Sigma = 0`;
- the axial profile `R0(w)` or `R0(s)`;
- the mouth geometry;
- the open exit;
- superfluid intake;
- DC leakage;
- outgoing radiation;
- thermal/export partition;
- superfluid density, EOS, sound speed, or density-dependent effective
  light-speed behavior;
- internal support profiles;
- full BdG/Krein signatures;
- Maxwell/mixed mode shapes;
- finite-mouth `P22` orientation;
- branch selection;
- nonlinear saturation.

They are useful because they are the compressed packet that the GR-like audit
can test.  They are not a substitute for the physical mechanism that produces
the packet.

---

## 9. Paper checklist

Every paper in this bundle should make the following physical points explicit
when relevant:

1. The defect is a finite brane-bulk throat/puncture, not a surface depression.
2. The mouth is the brane-side cross-section, not the entire defect.
3. The throat is open in the branch-realization picture: `R(L)>0`.
4. Hard caps are obsolete toy idealizations or negative controls.
5. The wall/throat dynamics are currently an effective closure unless
   `S_Sigma` is promoted.
6. Charge sign, circulation, radial throughput, and mouth output are distinct
   physical quantities.
7. Electromagnetic language should preserve the charge/circulation firewall:
   puncture orientation carries `eta_Q`, while circulation belongs to the
   magnetic/vortical sector.
8. The exact D/N support mode gives a mouth gradient/stress, not automatically
   a mouth current.
9. Nonzero DC mouth output requires open/radiative or other dynamic
   rectification physics.
10. Finite-mouth response starts with Hessian/tidal loading, not point scalar
   depth.
11. The first non-axisymmetric finite-mouth deformation is headless real `P22`.
12. Rigid mouth is a statement about entrance geometry, not automatically about
   all internal transfer variables.
13. Reduced/effective/conditional/strict-parent claims must be labeled.

---

## 10. Recommended paper paragraph

A compact version suitable for an introduction or model section:

```text
In this program a particle is modeled as a finite brane-bulk throat defect: a
finite-radius puncture of the brane into an open bulk conduit, not a surface
depression or capped cavity.  The brane-side mouth is the entrance cross-section
of the throat, while the interior carries support, wall, Maxwell/mixed, and
outgoing-port degrees of freedom.  The current moving-throat equations are an
effective wall/throat closure unless the promoted parent action S_Sigma is
explicitly added.  Observable particle data are extracted only after a branch is
frozen; the audit coefficients D0, C, P0, N2, and N4 are therefore projected
response moments, not the full physical ontology of the defect.
```

---

## 11. Common misreadings to avoid

Misreading:

```text
The defect is just a depression in the brane.
```

Correction:

```text
The defect is a finite-radius throat/puncture through the brane into the bulk.
```

Misreading:

```text
The throat is a capped cavity.
```

Correction:

```text
The branch-realization geometry is an open finite-radius conduit with R(L)>0.
```

Misreading:

```text
The mouth value is the mouth source.
```

Correction:

```text
On the D/N support branch the mouth value can vanish; the first nontrivial
mouth datum is the normal stress from the boundary gradient.
```

Misreading:

```text
P0 scalar output can directly generate P22 shape control.
```

Correction:

```text
The scalar P0 hammer does not linearly drive the area-preserving P22 mouth
mode.  P22 requires genuine finite-mouth shape/bracing physics.
```

Misreading:

```text
D0, C, P0, N2, and N4 are the physical model.
```

Correction:

```text
They are the compressed audit readouts of a branch.  The physical model is the
throat, mouth, support, wall, mixed, flux, and outgoing-port system that must
produce those readouts.
```
