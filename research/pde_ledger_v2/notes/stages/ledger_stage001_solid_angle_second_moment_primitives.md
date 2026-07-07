# ledger_stage001_solid_angle_second_moment_primitives

## Status

EXACT -- closed geometry.

## Purpose

This stage isolates the solid-angle and normalized second-moment primitives used
by the pathA_21c force calculation. It proves the geometric inputs without
using the Noether-stress force algebra itself.

## Provenance

Consuming context: `software/stage1_solver/reports/pathA_21c_force_from_noether_stress_tensor.md`
uses these standard sphere integrals in the reduced three-dimensional and bulk
four-dimensional force lanes. The same geometry is the input expected by
stage002's force assembly.

Script-backed status: the exact residuals are checked by
`scripts/ledger_stage001_solid_angle_second_moment_primitives_sympy_audit.py`
and independently by
`mathematica/ledger_stage001_solid_angle_second_moment_primitives_mathematica_audit.wl`.
Both scripts derive the quantities by explicit surface integration and include
able-to-fail mutation probes.

## 0. Why needed

In the pathA_21c matter-stress integral, the exterior field of one defect is
approximately constant across a small sphere surrounding the second defect.
The angular part of the surface traction then reduces to two classes of
geometry:

- the solid angle of the unit sphere controlling the Gauss-law normalization;
- the normalized second moments of the outward unit normal controlling the
  convective and pressure angular averages.

For an ambient dimension `d`, isotropy predicts
`<n_i n_j> = delta_ij/d`. This stage verifies the `d=3` and `d=4` instances
directly from parametrized integrals.

## 1. The S2 solid angle

Use the unit sphere in `R^3` with

```text
n(theta, phi) = (cos theta, sin theta cos phi, sin theta sin phi),
0 <= theta <= pi, 0 <= phi <= 2 pi.
```

The induced surface element is

```text
dS = sin theta dphi dtheta.
```

Therefore

```text
Omega_2
= int_0^pi int_0^(2 pi) sin theta dphi dtheta
= (2 pi) [-cos theta]_0^pi
= 4 pi.
```

## 2. The S2 normalized second moments

The first two normal components are

```text
n_1 = cos theta,
n_2 = sin theta cos phi.
```

Normalize all moments by `Omega_2 = 4 pi`. For the diagonal moment,

```text
<n_1^2>_S2
= (1/Omega_2) int_0^pi int_0^(2 pi) cos^2 theta sin theta dphi dtheta
= (2 pi / 4 pi) int_0^pi cos^2 theta sin theta dtheta.
```

With `u = cos theta`,

```text
int_0^pi cos^2 theta sin theta dtheta
= int_-1^1 u^2 du
= 2/3,
```

so

```text
<n_1^2>_S2 = (1/2)(2/3) = 1/3.
```

For the cross moment,

```text
<n_1 n_2>_S2
= (1/Omega_2) int_0^pi int_0^(2 pi)
   cos theta sin theta cos phi sin theta dphi dtheta.
```

The `phi` factor vanishes:

```text
int_0^(2 pi) cos phi dphi = 0,
```

hence

```text
<n_1 n_2>_S2 = 0.
```

Thus the verified `S^2` second-moment rule is
`<n_i n_j>_S2 = delta_ij/3`.

## 3. The S3 solid angle

Use the unit 3-sphere in `R^4` with

```text
n(chi, theta, phi)
= (cos chi,
   sin chi cos theta,
   sin chi sin theta cos phi,
   sin chi sin theta sin phi),
0 <= chi <= pi, 0 <= theta <= pi, 0 <= phi <= 2 pi.
```

The induced surface element is

```text
dS = sin^2 chi sin theta dphi dtheta dchi.
```

Therefore

```text
Omega_3
= int_0^pi int_0^pi int_0^(2 pi)
   sin^2 chi sin theta dphi dtheta dchi
= (2 pi)(2) int_0^pi sin^2 chi dchi
= 4 pi (pi/2)
= 2 pi^2.
```

## 4. The S3 normalized second moments

The first two normal components are

```text
n_1 = cos chi,
n_2 = sin chi cos theta.
```

Normalize all moments by `Omega_3 = 2 pi^2`. For the diagonal moment,

```text
<n_1^2>_S3
= (1/Omega_3) int_0^pi int_0^pi int_0^(2 pi)
   cos^2 chi sin^2 chi sin theta dphi dtheta dchi
= (4 pi / 2 pi^2) int_0^pi cos^2 chi sin^2 chi dchi.
```

Using `sin^2 chi cos^2 chi = (1/8)(1 - cos 4 chi)`,

```text
int_0^pi cos^2 chi sin^2 chi dchi = pi/8.
```

Hence

```text
<n_1^2>_S3 = (2/pi)(pi/8) = 1/4.
```

For the cross moment,

```text
<n_1 n_2>_S3
= (1/Omega_3) int_0^pi int_0^pi int_0^(2 pi)
   cos chi sin chi cos theta sin^2 chi sin theta dphi dtheta dchi.
```

The `theta` factor vanishes:

```text
int_0^pi cos theta sin theta dtheta = 0,
```

so

```text
<n_1 n_2>_S3 = 0.
```

Thus the verified `S^3` second-moment rule is
`<n_i n_j>_S3 = delta_ij/4`.

## What this achieves physically

The `Omega_2 = 4 pi` and `Omega_3 = 2 pi^2` results set the Gauss-law
normalizations that produce the reduced `r^-2` and bulk `R^-3` force powers in
pathA_21c. The second moments supply the angular averages behind the
convective factor `-(1 + 1/d)` and the pressure contribution `1/d`, with
`d = 3` for the reduced lane and `d = 4` for the bulk lane.

## What is still missing

Nothing for this stage. The geometry primitives are closed.

## Next step

Stage002 assembles these primitives into the pathA_21c matter-stress force
calculation.
