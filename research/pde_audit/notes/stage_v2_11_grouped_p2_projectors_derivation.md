# Stage V2-11 — Grouped Real \(P_2\) Projector Calculus

## 0. Purpose

This stage freezes the algebra of the grouped real \(P_2\) bundle before the next normalization and source-map stages.

The goal is not to solve a radial/axial moving-throat branch. The goal is to verify the exact bookkeeping used whenever the five real quadrupole lanes

\[
(20,\;21c,\;21s,\;22c,\;22s)
\]

are compressed into the grouped triple

\[
x=(x_{20},x_{21},x_{22}).
\]

The \(21\) and \(22\) grouped lanes each contain two real components, so the grouped triple is not Euclidean with equal weights. Its natural metric is

\[
G_{\rm grp}=\operatorname{diag}(1,2,2).
\]

This stage verifies:

1. the grouped metric;
2. the trace/anomaly basis;
3. the exact projectors;
4. the inverse grouped-coordinate map;
5. the anisotropy norm;
6. the weak-axisymmetric splitting law;
7. first-order transport formulas for \(u_2=-D_2/D_0\) and \(P_0=N_0/D_0\).

The SymPy audit passes all checks.

---

## 1. Grouped metric from the five real lanes

Let

\[
y=(y_{20},y_{21c},y_{21s},y_{22c},y_{22s})^T
\]

and impose the grouped embedding

\[
y_{20}=x_{20},\qquad
y_{21c}=y_{21s}=x_{21},\qquad
y_{22c}=y_{22s}=x_{22}.
\]

The embedding matrix is

\[
E=
\begin{pmatrix}
1&0&0\\
0&1&0\\
0&1&0\\
0&0&1\\
0&0&1
\end{pmatrix}.
\]

Therefore the grouped inner product inherited from the five-lane Euclidean inner product is

\[
G_{\rm grp}=E^T E
=
\begin{pmatrix}
1&0&0\\
0&2&0\\
0&0&2
\end{pmatrix}.
\]

So every grouped trace, orthogonality condition, and norm must use \(G_{\rm grp}\), not the naive identity metric.

---

## 2. Trace/anomaly basis

The useful grouped basis is

\[
e_{\rm tr}=(1,1,1)^T,
\]

\[
e_a=(4,-1,-1)^T,
\]

\[
e_b=(0,1,-1)^T.
\]

The script verifies exact \(G_{\rm grp}\)-orthogonality:

\[
e_{\rm tr}^TG_{\rm grp}e_a=0,
\qquad
e_{\rm tr}^TG_{\rm grp}e_b=0,
\qquad
e_a^TG_{\rm grp}e_b=0.
\]

The exact squared norms are

\[
e_{\rm tr}^TG_{\rm grp}e_{\rm tr}=5,
\]

\[
e_a^TG_{\rm grp}e_a=20,
\]

\[
e_b^TG_{\rm grp}e_b=4.
\]

Thus every grouped triple decomposes uniquely as

\[
x=\bar x\,e_{\rm tr}+a_x\,e_a+b_x\,e_b.
\]

The exact coordinates are

\[
\boxed{
\bar x=\frac{x_{20}+2x_{21}+2x_{22}}5,
}
\]

\[
\boxed{
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
}
\]

\[
\boxed{
b_x=\frac{x_{21}-x_{22}}2.
}
\]

The inverse map is

\[
\boxed{
x_{20}=\bar x+4a_x,
}
\]

\[
\boxed{
x_{21}=\bar x-a_x+b_x,
}
\]

\[
\boxed{
x_{22}=\bar x-a_x-b_x.
}
\]

---

## 3. Exact projectors

The \(G_{\rm grp}\)-orthogonal projector onto a basis vector \(e\) is

\[
P_e=\frac{e\,e^TG_{\rm grp}}{e^TG_{\rm grp}e}.
\]

The script verifies the three exact projectors:

\[
P_{\rm tr}
=
\frac15
\begin{pmatrix}
1&2&2\\
1&2&2\\
1&2&2
\end{pmatrix},
\]

\[
P_a
=
\frac1{20}
\begin{pmatrix}
16&-8&-8\\
-4&2&2\\
-4&2&2
\end{pmatrix},
\]

\[
P_b
=
\frac14
\begin{pmatrix}
0&0&0\\
0&2&-2\\
0&-2&2
\end{pmatrix}.
\]

They obey

\[
P_{\rm tr}+P_a+P_b=I_3,
\]

\[
P_i^2=P_i,
\]

\[
P_iP_j=0\qquad(i\neq j),
\]

and the \(G_{\rm grp}\)-self-adjointness condition

\[
P_i^TG_{\rm grp}=G_{\rm grp}P_i.
\]

So these are true \(G_{\rm grp}\)-orthogonal projectors.

---

## 4. Isotropy and anisotropy norm

The isotropic grouped branch is exactly

\[
a_x=0,\qquad b_x=0,
\]

which is equivalent to

\[
x_{20}=x_{21}=x_{22}.
\]

The grouped norm decomposes as

\[
x^TG_{\rm grp}x
=
5\bar x^2+20a_x^2+4b_x^2.
\]

The normalized anisotropy norm used downstream is the anisotropic part divided by the total grouped multiplicity \(5\):

\[
\boxed{
A_x^2
=
\frac{20a_x^2+4b_x^2}{5}
=
4a_x^2+\frac45 b_x^2.
}
\]

This is the exact norm behind the later grouped anisotropy measures.

---

## 5. Weak-axisymmetric splitting fingerprint

For a weak axisymmetric quadrupolar perturbation, the grouped splitting vector is

\[
\lambda=(1,\tfrac12,-1).
\]

The script verifies

\[
\bar\lambda=0,
\]

\[
a_\lambda=\frac14,
\]

\[
b_\lambda=\frac34.
\]

Therefore

\[
\boxed{
b_\lambda=3a_\lambda.
}
\]

So any first-order grouped coefficient of the form

\[
x_A=x^{(0)}+\epsilon\lambda_Ax^{(1)}
\]

has

\[
\bar x=x^{(0)},
\]

\[
a_x=\frac{\epsilon}{4}x^{(1)},
\]

\[
b_x=\frac{3\epsilon}{4}x^{(1)},
\]

and hence

\[
\boxed{b_x=3a_x.}
\]

This is the diagnostic line for pure weak-axisymmetric \(l=2\) symmetry breaking.

---

## 6. First-order transport for the normalized response \(u_2\)

Let

\[
u_2=-\frac{D_2}{D_0}.
\]

Perturb a grouped lane by

\[
D_{A0}=D_0+\epsilon\lambda_A\delta D_0,
\]

\[
D_{A2}=D_2+\epsilon\lambda_A\delta D_2.
\]

Then

\[
u_2^{(A)}
=
-\frac{D_{A2}}{D_{A0}}
=
u_2+\epsilon\lambda_A u_2^{(1)}+O(\epsilon^2),
\]

with

\[
\boxed{
u_2^{(1)}
=
-\frac{\delta D_2+u_2\,\delta D_0}{D_0}.
}
\]

Equivalently, using \(u_2=-D_2/D_0\),

\[
u_2^{(1)}
=
-\frac{D_0\delta D_2-D_2\delta D_0}{D_0^2}.
\]

On the weak-axisymmetric line, the \(u_2\) anisotropy obeys

\[
b_{u_2}=3a_{u_2}.
\]

---

## 7. First-order transport for the outgoing prefactor \(P_0\)

Let

\[
P_0=\frac{N_0}{D_0}.
\]

Perturb a grouped lane by

\[
N_{A0}=N_0+\epsilon\lambda_A\delta N_0,
\]

\[
D_{A0}=D_0+\epsilon\lambda_A\delta D_0.
\]

Then

\[
P_0^{(A)}
=
\frac{N_{A0}}{D_{A0}}
=
P_0+\epsilon\lambda_A P_0^{(1)}+O(\epsilon^2),
\]

with

\[
\boxed{
P_0^{(1)}
=
\frac{\delta N_0-P_0\delta D_0}{D_0}.
}
\]

Equivalently,

\[
P_0^{(1)}
=
\frac{D_0\delta N_0-N_0\delta D_0}{D_0^2}.
\]

On the weak-axisymmetric line, the outgoing-prefactor anisotropy also obeys

\[
b_{P_0}=3a_{P_0}.
\]

This is the exact \(P_1/P_0\) bookkeeping used later when \(\Xi_1\) is interpreted as the weak-axisymmetric outgoing-prefactor slope.

---

## 8. Result and carry-forward status

The script reports:

```text
checks_total: 41
checks_passed: 41
checks_failed: 0
```

So V2-11 gives a clean pass.

The frozen carry-forward packet is:

\[
G_{\rm grp}=\operatorname{diag}(1,2,2),
\]

\[
x=\bar x(1,1,1)+a_x(4,-1,-1)+b_x(0,1,-1),
\]

\[
A_x^2=4a_x^2+\frac45b_x^2,
\]

\[
\lambda_{\rm ax}=(1,\tfrac12,-1),
\qquad
b=3a,
\]

\[
u_2^{(1)}
=
-\frac{\delta D_2+u_2\delta D_0}{D_0},
\]

\[
P_0^{(1)}
=
\frac{\delta N_0-P_0\delta D_0}{D_0}.
\]

This stage is independent of the radial/axial open-throat boundary correction from V2-04. The open organ-pipe patch changes the radial/axial branch data inside \(D_n\) and \(N_n\), but not the grouped angular projector calculus.
