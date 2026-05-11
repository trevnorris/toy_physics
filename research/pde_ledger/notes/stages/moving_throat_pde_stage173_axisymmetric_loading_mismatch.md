# Moving-Throat PDE — Stage 173: Weak-Axisymmetric Transport of the Physical Slopes and Collapse to One Static Loading Mismatch

## Purpose

Stage 240 reduced the linear grouped outlet problem to the weak-axisymmetric physical slope pair
\[
u_2^{(1)},
\qquad
\frac{P_1}{P_0},
\]
and, on the even-preserving branch, to the single logarithmic prefactor slope \(P_1/P_0\).

That was already a strong collapse, but it still left the physical slopes looking like abstract first derivatives of the grouped response.

The next honest step is therefore:

> compute those slopes directly from the actual weak-axisymmetric grouped moving-throat response.

This stage does exactly that.

Using the exact Stage-7 axisymmetric grouped signature
\[
\lambda_{20}=1,\qquad \lambda_{21}=\frac12,\qquad \lambda_{22}=-1,
\]
I expand the actual grouped conservative operator moments and the static outgoing-transfer moment themselves,
\[
D_{A,0}=D_0+\epsilon\,\lambda_A D_{01}+O(\epsilon^2),
\]
\[
D_{A,2}=D_2+\epsilon\,\lambda_A D_{21}+O(\epsilon^2),
\qquad
D_{A,4}=D_4+\epsilon\,\lambda_A D_{41}+O(\epsilon^2),
\]
\[
N_{A,0}=N_0+\epsilon\,\lambda_A N_{01}+O(\epsilon^2).
\]

From those actual grouped response coefficients one gets the exact physical slope laws
\[
\boxed{
u_2^{(1)}=-\frac{D_{21}+u_2\,D_{01}}{D_0},
}
\]
\[
\boxed{
u_4^{(1)}=-\frac{5D_{01}+18D_{21}+81D_{41}}{81\,D_0}
\qquad
\text{on the canonical branch }(u_2,u_4)=\left(\frac19,\frac{4}{81}\right),
}
\]
and
\[
\boxed{
\frac{P_1}{P_0}
=
\frac{N_{01}}{N_0}
-
\frac{D_{01}}{D_0}.
}
\]

So the physical grouped slopes are no longer vague “measure the response somehow” objects.
They are the exact weak-axisymmetric logarithmic transport coefficients of the grouped conservative operator and the grouped outgoing transfer.

The sharpest consequence is on the even-preserving branch.
If the canonical conservative even fingerprint is preserved at first order,
\[
u_2^{(1)}=0,
\]
then
\[
\boxed{
D_{21}=-\frac{D_{01}}{9},
}
\]
and the hidden-even consistency relation
\[
u_4^{(1)}=\frac89\,u_2^{(1)}
\]
forces
\[
\boxed{
D_{41}=-\frac{D_{01}}{27}.
}
\]

So on the even-preserving branch the entire conservative grouped response is transported by **one** static operator slope \(D_{01}\), while the remaining linear grouped `2.5`PN defect collapses to one scalar mismatch
\[
\boxed{
\Xi_{\rm load}:=
\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
}
\]

Then
\[
\boxed{
\Delta_Q^{(20)}=\epsilon\,\Xi_{\rm load},
\qquad
\Delta_Q^{(21)}=\frac{\epsilon}{2}\,\Xi_{\rm load},
\qquad
\Delta_Q^{(22)}=-\epsilon\,\Xi_{\rm load}.
}
\]

So the remaining linear grouped normalization problem is now smaller than Stage 240 suggested.
It is no longer “compute \(u_2^{(1)}\) and \(P_1/P_0\) separately.”
On the even-preserving branch it is just:

> compute the single weak-axisymmetric **static loading mismatch** \(\Xi_{\rm load}\).

---

## 1. Exact weak-axisymmetric transport of the grouped conservative response

For each grouped lane \(A\in\{20,21,22\}\), the physical grouped response data are
\[
u_2^{(A)}=-\frac{D_{A,2}}{D_{A,0}},
\qquad
u_4^{(A)}=\frac{D_{A,2}^2-D_{A,0}D_{A,4}}{D_{A,0}^2},
\qquad
P_0^{(A)}=\frac{N_{A,0}}{D_{A,0}}.
\]

On the weak axisymmetric \(Y_{20}\) branch, write
\[
D_{A,0}=D_0+\epsilon\,\lambda_A D_{01}+O(\epsilon^2),
\]
\[
D_{A,2}=D_2+\epsilon\,\lambda_A D_{21}+O(\epsilon^2),
\]
\[
D_{A,4}=D_4+\epsilon\,\lambda_A D_{41}+O(\epsilon^2),
\]
\[
N_{A,0}=N_0+\epsilon\,\lambda_A N_{01}+O(\epsilon^2),
\]
with
\[
\lambda_{20}=1,\qquad \lambda_{21}=\frac12,\qquad \lambda_{22}=-1.
\]

Then the grouped physical slopes are defined by
\[
u_2^{(A)}=u_2+\epsilon\,\lambda_A u_2^{(1)}+O(\epsilon^2),
\]
\[
u_4^{(A)}=u_4+\epsilon\,\lambda_A u_4^{(1)}+O(\epsilon^2),
\]
\[
P_0^{(A)}=P_0+\epsilon\,\lambda_A P_1+O(\epsilon^2).
\]

Expanding the exact definitions gives
\[
\boxed{
u_2^{(1)}=-\frac{D_{21}+u_2 D_{01}}{D_0},
}
\]
\[
\boxed{
P_1
=
\frac{D_0 N_{01}-N_0 D_{01}}{D_0^2},
\qquad
\frac{P_1}{P_0}=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
}
\]

So the Stage 240 physical slope pair is now computed directly from the actual grouped moving-throat response coefficients.

---

## 2. Exact canonical-even transport law

On the canonical compensated branch one has
\[
u_2=\frac19,
\qquad
u_4=\frac{4}{81}.
\]

Expanding the exact grouped \(u_4\) definition gives
\[
\boxed{
u_4^{(1)}=
-\frac{5D_{01}+18D_{21}+81D_{41}}{81\,D_0}.
}
\]

Stage 240 already translated the one-parameter hidden-even relation into
\[
u_4^{(1)}=\frac89\,u_2^{(1)}.
\]
Substituting the explicit \(u_2^{(1)}\) and \(u_4^{(1)}\) formulas gives the equivalent operator law
\[
\boxed{
D_{41}=\frac23 D_{21}+\frac1{27}D_{01}.
}
\]

So the hidden-even consistency relation is now written directly in the weak-axisymmetric operator-slope variables of the actual grouped response.

---

## 3. Conservative collapse on the even-preserving branch

Now impose preservation of the canonical conservative even fingerprint:
\[
u_2^{(1)}=0.
\]
Then
\[
\boxed{
D_{21}=-\frac{D_{01}}{9}.
}
\]

Substituting into the hidden-even operator law above yields
\[
\boxed{
D_{41}=-\frac{D_{01}}{27}.
}
\]

So on the even-preserving branch, all three conservative grouped operator slopes are fixed by a single static slope \(D_{01}\):
\[
D_{A,0}=D_0\bigl[1+\epsilon\,\lambda_A \delta_D\bigr]+O(\epsilon^2),
\]
\[
D_{A,2}=D_2\bigl[1+\epsilon\,\lambda_A \delta_D\bigr]+O(\epsilon^2),
\]
\[
D_{A,4}=D_4\bigl[1+\epsilon\,\lambda_A \delta_D\bigr]+O(\epsilon^2),
\]
with
\[
\delta_D:=\frac{D_{01}}{D_0},
\qquad
D_2=-\frac{D_0}{9},
\qquad
D_4=-\frac{D_0}{27}.
\]

Equivalently,
\[
u_2^{(A)}=\frac19+O(\epsilon^2),
\qquad
u_4^{(A)}=\frac{4}{81}+O(\epsilon^2).
\]

So the weak-axisymmetric conservative branch is frozen to canonical order once \(u_2^{(1)}=0\).
The remaining first-order physics is only in the static loading.

---

## 4. The remaining linear grouped defect is one static loading mismatch

Define the logarithmic outgoing-transfer slope
\[
\delta_N:=\frac{N_{01}}{N_0},
\]
and the logarithmic static-operator slope
\[
\delta_D:=\frac{D_{01}}{D_0}.
\]

Then the exact prefactor slope becomes
\[
\boxed{
\frac{P_1}{P_0}=\delta_N-\delta_D.
}
\]

So on the even-preserving branch the only surviving linear grouped quadrupole-normalization defect is
\[
\boxed{
\Xi_{\rm load}:=\delta_N-\delta_D.
}
\]

In lane form,
\[
P_0^{(A)}=P_0\bigl[1+\epsilon\,\lambda_A\,\Xi_{\rm load}\bigr]+O(\epsilon^2),
\]
and therefore
\[
\boxed{
\Delta_Q^{(20)}=\epsilon\,\Xi_{\rm load},
\qquad
\Delta_Q^{(21)}=\frac{\epsilon}{2}\,\Xi_{\rm load},
\qquad
\Delta_Q^{(22)}=-\epsilon\,\Xi_{\rm load}.
}
\]

This is the main theorem of the stage.

It says that the remaining linear grouped `2.5`PN problem is a **single scalar mismatch between two static logarithmic slopes**.

---

## 5. Microscopic decomposition of the two static slopes

From the full grouped bundle,
\[
D_{A,0}=K_A-B_{A,0}-Z_{A,0},
\qquad
D_{A,2}=-(M_A+B_{A,2}+Z_{A,2}),
\qquad
N_{A,0}=N_A.
\]

So on the weak axisymmetric branch
\[
D_{01}=K_1-B_0^{(1)}-Z_0^{(1)},
\]
\[
D_{21}=-(M_1+B_2^{(1)}+Z_2^{(1)}),
\]
\[
N_{01}=N_1.
\]

Therefore the physical slope pair becomes
\[
\boxed{
u_2^{(1)}
=
-\frac{1}{D_0}
\left[
\frac19 K_1 - M_1
-\left(B_2^{(1)}+\frac19 B_0^{(1)}\right)
-\left(Z_2^{(1)}+\frac19 Z_0^{(1)}\right)
\right],
}
\]
\[
\boxed{
\frac{P_1}{P_0}
=
\frac{N_1}{N_0}
-
\frac{K_1-B_0^{(1)}-Z_0^{(1)}}{D_0}.
}
\]

So the Stage 239 obstruction pair has now been rewritten in the sharpest static form yet:

- one weak-axisymmetric conservative static slope,
- one weak-axisymmetric outgoing-transfer static slope.

On the even-preserving branch, only their difference matters.

---

## 6. Portwise logarithmic form of the outgoing-transfer slope

For each Maxwell/mixed port \(r\), write
\[
N_{0}^{(r)}=\frac{P_r^2}{\Delta_r^2},
\qquad
N_{1}^{(r)}=
\frac{2P_r P_{1r}}{\Delta_r^2}
-
\frac{2P_r^2 \Delta_{1r}}{\Delta_r^3}.
\]

Then
\[
\boxed{
\frac{N_1}{N_0}
=
\sum_r w_r
\left(
2\frac{P_{1r}}{P_r}
-
2\frac{\Delta_{1r}}{\Delta_r}
\right),
\qquad
w_r:=\frac{N_0^{(r)}}{\sum_s N_0^{(s)}}.
}
\]

So the outgoing-transfer slope is a weighted average of the portwise logarithmic deformation
\[
2\,\delta\ln P_r - 2\,\delta\ln\Delta_r.
\]

Likewise, for the conservative static Maxwell/mixed sector
\[
Z_{0}^{(r)}=\frac{Q_r}{\Delta_r},
\qquad
Z_{1}^{(r)}=
\frac{\Delta_r Q_{1r}-Q_r \Delta_{1r}}{\Delta_r^2},
\]
so
\[
\boxed{
Z_0^{(1)}
=
\sum_r Z_0^{(r)}
\left(
\frac{Q_{1r}}{Q_r}
-
\frac{\Delta_{1r}}{\Delta_r}
\right).
}
\]

This is the first reduced microscopic form in which the actual moving-throat branch can be interrogated directly:
the whole linear grouped defect is controlled by weighted logarithmic slopes of the static conservative and outgoing-transfer port data.

---

## 7. What Stage 241 changes

Before this stage, the next theorem gate was still phrased as:

> compute the physical grouped slopes \(u_2^{(1)}\) and \(P_1/P_0\).

After this stage, that is no longer the sharpest statement.

The new theorem status is:

1. the weak-axisymmetric physical slopes are exact transport coefficients of the actual grouped response moments,
   \[
   u_2^{(1)}=-\frac{D_{21}+u_2 D_{01}}{D_0},
   \qquad
   \frac{P_1}{P_0}=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0};
   \]
2. on the even-preserving branch,
   \[
   D_{21}=-\frac{D_{01}}{9},
   \qquad
   D_{41}=-\frac{D_{01}}{27},
   \]
   so the conservative grouped response is fully transported by one static slope \(D_{01}\);
3. the remaining linear grouped `2.5`PN defect is one scalar
   \[
   \Xi_{\rm load}=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
   \]

So the next honest theorem gate is now smaller than Stage 240 suggested:

> compute the weak-axisymmetric static operator slope \(D_{01}/D_0\) and the weak-axisymmetric static outgoing-transfer slope \(N_{01}/N_0\) on the actual moving-throat branch — and, on the even-preserving branch, only their difference \(\Xi_{\rm load}\).

That is the direct continuation point.
