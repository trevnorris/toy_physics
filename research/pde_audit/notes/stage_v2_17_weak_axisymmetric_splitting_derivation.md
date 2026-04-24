# Stage V2-17 — Weak-axisymmetric grouped-\(P_2\) splitting audit

## 0. Purpose and status

This stage audits the weak-axisymmetric grouped-\(P_2\) law that has been carried through the moving-throat program:

\[
(20,21,22)\sim \left(1,\frac12,-1\right),
\qquad
b=3a.
\]

The stage is **exact algebra inside the weak-axisymmetric grouped real \(P_2\) closure**. It does not solve the nonlinear moving-throat PDE. Its job is narrower:

1. derive the angular origin of the lane signature;
2. propagate that signature through the grouped trace/anomaly variables;
3. derive the first-order conservative response slope \(u_2^{(1)}\);
4. derive the outgoing-prefactor slope
   \[
   \Xi_1=\frac{P_1}{P_0};
   \]
5. identify the compensated branch where conservative even splitting is removed and only the prefactor-loading scalar remains.

The accompanying script is:

```text
stage_v2_17_weak_axisymmetric_splitting_sympy_audit.py
```

and its run output is:

```text
stage_v2_17_weak_axisymmetric_splitting_output.txt
```

The script reports:

```text
total_checks: 54
passed_checks: 54
failed_checks: 0
```

---

## 1. Angular source of the weak-axisymmetric signature

Use the normalized real \(l=2\) harmonics on the unit sphere:

\[
Y_{20}=\sqrt{\frac{5}{16\pi}}\,(2z^2-x^2-y^2),
\]

\[
Y_{21c}=\sqrt{\frac{15}{4\pi}}xz,
\qquad
Y_{21s}=\sqrt{\frac{15}{4\pi}}yz,
\]

\[
Y_{22c}=\sqrt{\frac{15}{16\pi}}(x^2-y^2),
\qquad
Y_{22s}=\sqrt{\frac{15}{4\pi}}xy.
\]

The script verifies orthonormality:

\[
\int_{S^2}Y_A Y_B\,d\Omega=\delta_{AB}.
\]

A weak axisymmetric quadrupole perturbation is proportional to \(Y_{20}\), so the lane-splitting matrix is controlled by

\[
M_{AA}^{(20)}=\int_{S^2}Y_A\,Y_{20}\,Y_A\,d\Omega.
\]

Using exact sphere monomial moments, the script obtains

\[
\int Y_{20}Y_{20}Y_{20}\,d\Omega
=
\frac{\sqrt5}{7\sqrt\pi},
\]

\[
\int Y_{21c}Y_{20}Y_{21c}\,d\Omega
=
\int Y_{21s}Y_{20}Y_{21s}\,d\Omega
=
\frac{\sqrt5}{14\sqrt\pi},
\]

\[
\int Y_{22c}Y_{20}Y_{22c}\,d\Omega
=
\int Y_{22s}Y_{20}Y_{22s}\,d\Omega
=
-\frac{\sqrt5}{7\sqrt\pi}.
\]

Therefore

\[
M^{(20)}
=
\kappa_\ast
\operatorname{diag}
\left(
1,\frac12,\frac12,-1,-1
\right),
\qquad
\kappa_\ast=\frac{\sqrt5}{7\sqrt\pi}.
\]

After grouping \(21c,21s\) into the \(21\) lane and \(22c,22s\) into the \(22\) lane, the grouped signature is

\[
\lambda_{20}=1,
\qquad
\lambda_{21}=\frac12,
\qquad
\lambda_{22}=-1.
\]

This is the angular origin of the weak-axisymmetric law.

---

## 2. Grouped trace/anomaly variables

The grouped metric is

\[
G_{\rm grp}=\operatorname{diag}(1,2,2).
\]

For any grouped vector

\[
x=(x_{20},x_{21},x_{22}),
\]

define

\[
\bar x=\frac{x_{20}+2x_{21}+2x_{22}}5,
\]

\[
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\]

\[
b_x=\frac{x_{21}-x_{22}}2.
\]

The inverse map is

\[
x_{20}=\bar x+4a_x,
\]

\[
x_{21}=\bar x-a_x+b_x,
\]

\[
x_{22}=\bar x-a_x-b_x.
\]

For

\[
\lambda=\left(1,\frac12,-1\right),
\]

the script verifies

\[
\bar\lambda=0,
\qquad
a_\lambda=\frac14,
\qquad
b_\lambda=\frac34.
\]

Thus

\[
\boxed{b_\lambda=3a_\lambda.}
\]

So for any weak coefficient split

\[
x_A=x_0+\epsilon\lambda_Ax_1,
\]

one gets

\[
\bar x=x_0,
\]

\[
a_x=\frac{\epsilon x_1}{4},
\]

\[
b_x=\frac{3\epsilon x_1}{4},
\]

and therefore

\[
\boxed{b_x=3a_x.}
\]

This applies to \(D\)-moments, \(N\)-moments, response moments, prefactor moments, and any other grouped scalar whose first-order lane dependence is sourced by a pure axisymmetric \(l=2\) perturbation.

---

## 3. Conservative response transport

Let the conservative grouped operator moments split as

\[
D_{A0}=D_0+\epsilon\lambda_A D_{01},
\]

\[
D_{A2}=D_2+\epsilon\lambda_A D_{21},
\]

\[
D_{A4}=D_4+\epsilon\lambda_A D_{41}.
\]

The normalized conservative response is

\[
Y_A(\omega)
=
\frac{D_{A0}}{D_{A0}+D_{A2}\omega^2+D_{A4}\omega^4+\cdots}
=
1+u_{2,A}\omega^2+u_{4,A}\omega^4+\cdots.
\]

At zeroth order,

\[
u_2=-\frac{D_2}{D_0}.
\]

At first order,

\[
u_{2,A}
=
u_2+\epsilon\lambda_Au_2^{(1)}+O(\epsilon^2),
\]

with

\[
\boxed{
u_2^{(1)}
=
-\frac{D_{21}+u_2D_{01}}{D_0}.
}
\]

The script verifies the equivalent expanded form

\[
u_2^{(1)}
=
-\frac{D_0D_{21}-D_{01}D_2}{D_0^2}.
\]

Therefore the first conservative even-preserving condition is

\[
u_2^{(1)}=0
\quad\Longleftrightarrow\quad
D_{21}+u_2D_{01}=0.
\]

On the canonical outgoing branch where

\[
u_2=\frac19,
\]

this becomes

\[
\boxed{
D_{21}=-\frac{D_{01}}9.
}
\]

---

## 4. Hidden-even relation at \(O(\omega^4)\)

The next normalized response coefficient is

\[
u_4=\frac{D_2^2-D_0D_4}{D_0^2}.
\]

The canonical compact outgoing inverse-DtN branch has

\[
u_2=\frac19,
\qquad
u_4=\frac4{81}.
\]

Equivalently, in operator moments,

\[
D_2=-\frac{D_0}{9},
\qquad
D_4=-\frac{D_0}{27}.
\]

The one-pole/hidden-even first-order relation is the linearized form of

\[
u_4=4u_2^2.
\]

Since \(u_2=1/9\), the linearization is

\[
u_4^{(1)}=\frac89u_2^{(1)}.
\]

The script derives

\[
u_2^{(1)}
=
-\frac{D_{01}+9D_{21}}{9D_0},
\]

\[
u_4^{(1)}
=
-\frac{5D_{01}+18D_{21}+81D_{41}}{81D_0},
\]

and

\[
u_4^{(1)}-\frac89u_2^{(1)}
=
\frac{D_{01}+18D_{21}-27D_{41}}{27D_0}.
\]

Therefore the hidden-even gate is

\[
D_{01}+18D_{21}-27D_{41}=0,
\]

or

\[
\boxed{
D_{41}=\frac23D_{21}+\frac{D_{01}}{27}.
}
\]

On the even-preserving branch

\[
D_{21}=-\frac{D_{01}}9,
\]

the hidden-even condition reduces to

\[
\boxed{
D_{41}=-\frac{D_{01}}{27}.
}
\]

Thus the compensated conservative branch is

\[
\boxed{
D_{21}=-\frac{D_{01}}9,
\qquad
D_{41}=-\frac{D_{01}}{27}.
}
\]

---

## 5. Outgoing-prefactor transport and \(\Xi_1\)

Let the outgoing transfer numerator split as

\[
N_{A0}=N_0+\epsilon\lambda_A N_{01}.
\]

The static outgoing prefactor is

\[
P_A=\frac{N_{A0}}{D_{A0}}.
\]

Expanding,

\[
P_A=P_0+\epsilon\lambda_A P_1+O(\epsilon^2),
\]

with

\[
P_0=\frac{N_0}{D_0},
\]

\[
\boxed{
P_1=
\frac{D_0N_{01}-D_{01}N_0}{D_0^2}.
}
\]

The dimensionless transported prefactor slope is

\[
\boxed{
\Xi_1=\frac{P_1}{P_0}.
}
\]

Substituting \(P_0=N_0/D_0\),

\[
\boxed{
\Xi_1
=
\frac{N_{01}}{N_0}
-
\frac{D_{01}}{D_0}.
}
\]

This is the same load-mismatch form used later in the monomial/similarity-orbit package.

The prefactor lanes obey the same grouped weak-axisymmetric law:

\[
P_{20}=P_0+\epsilon P_1,
\]

\[
P_{21}=P_0+\frac{\epsilon}{2}P_1,
\]

\[
P_{22}=P_0-\epsilon P_1.
\]

Therefore

\[
\bar P=P_0,
\]

\[
a_P=\frac{\epsilon P_1}{4},
\]

\[
b_P=\frac{3\epsilon P_1}{4},
\]

and

\[
\boxed{
b_P=3a_P.
}
\]

The prefactor-isotropy condition is

\[
\boxed{
\Xi_1=0
\quad\Longleftrightarrow\quad
\frac{N_{01}}{N_0}=\frac{D_{01}}{D_0}.
}
\]

---

## 6. Final V2-17 gate packet

The weak-axisymmetric first-order gate packet is:

\[
\lambda=(1,\tfrac12,-1),
\]

\[
b=3a,
\]

\[
u_2^{(1)}
=
-\frac{D_{21}+u_2D_{01}}{D_0},
\]

\[
\Xi_1
=
\frac{P_1}{P_0}
=
\frac{N_{01}}{N_0}
-
\frac{D_{01}}{D_0}.
\]

On the canonical compensated branch,

\[
D_{21}=-\frac{D_{01}}9,
\]

\[
D_{41}=-\frac{D_{01}}{27},
\]

and the only remaining first-order weak-axisymmetric normalization defect is

\[
\boxed{
\Xi_1
=
\frac{N_{01}}{N_0}
-
\frac{D_{01}}{D_0}.
}
\]

This is the scalar that future actual-branch calculations must compute.

---

## 7. Interpretation

V2-17 passes.

A pure weak axisymmetric \(l=2\) perturbation does not produce an arbitrary grouped anisotropy. It produces one fixed lane signature:

\[
(20,21,22)\sim(1,\tfrac12,-1).
\]

Therefore every first-order grouped defect produced by that angular mechanism must satisfy

\[
b=3a.
\]

After compensating the conservative even response, the remaining weak-axisymmetric first-order normalization problem is no longer a three-lane problem. It is one scalar loading mismatch:

\[
\Xi_1=\frac{P_1}{P_0}.
\]

That is the continuation target for the quotient/orbit stages.
