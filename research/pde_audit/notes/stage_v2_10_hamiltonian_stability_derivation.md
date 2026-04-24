# Stage V2-10 — Hamiltonian / Stability Audit

## 0. Purpose

This stage combines the conservative wall, BdG, and localized Maxwell/mixed
reduced sectors into one stability ledger.

Earlier stages verified the individual Schur complements:

- stable BdG support modes soften the static wall stiffness and raise the
  effective inertia;
- the localized Maxwell/mixed block gives a conservative self-energy and a
  nonnegative outgoing-transfer factor;
- the outgoing `l=2` port supplies the odd `i omega^5` fingerprint.

The present stage asks the next necessary question:

> When all conservative pieces are active at once, what exact positivity gates
> keep the reduced branch from becoming a ghost, tachyon, or unstable
> near-softened denominator?

The result is a compact pass condition:

\[
M>0,\qquad \varpi^2>0,\qquad
\Omega_U^2>0,\qquad \Omega_W^2>0,
\]
\[
\Delta=\Omega_U^2\Omega_W^2-R^2>0,
\]
\[
D_0
=
K-\frac{c_B^2}{\varpi^2}
-
\frac{
g_U^2\Omega_W^2+2g_Ug_WR+g_W^2\Omega_U^2
}{
\Omega_U^2\Omega_W^2-R^2
}
>0.
\]

If these hold, the conservative reduced Hamiltonian is positive definite.

---

## 1. Conservative one-lane bundle

Use one representative grouped lane.  The reduced variables are

\[
q,\qquad X,\qquad U,\qquad W,
\]

where

- \(q\) is the wall/worldtube amplitude,
- \(X\) is a positive-energy BdG support coordinate,
- \(U\) is a brane-like localized Maxwell coordinate,
- \(W\) is the mixed \(A_w/F_{\mu w}/J^w\)-active coordinate.

The reduced Lagrangian is

\[
L
=
\frac12M\dot q^2-\frac12Kq^2
+
\frac12\dot X^2-\frac12\varpi^2X^2
+
c_BqX
\]
\[
+
\frac12\dot U^2-\frac12\Omega_U^2U^2
+
\frac12\dot W^2-\frac12\Omega_W^2W^2
+
RUW+g_UqU+g_WqW.
\]

The conservative potential Hessian for \((q,X,U,W)\) is

\[
P=
\begin{pmatrix}
K&-c_B&-g_U&-g_W\\
-c_B&\varpi^2&0&0\\
-g_U&0&\Omega_U^2&-R\\
-g_W&0&-R&\Omega_W^2
\end{pmatrix}.
\]

The kinetic matrix is

\[
T=\operatorname{diag}(M,1,1,1).
\]

So positive Hamiltonian kinetic energy requires \(M>0\).

---

## 2. Exact Schur-complement stability gate

The internal conservative block is

\[
A_{\rm int}
=
\begin{pmatrix}
\varpi^2&0&0\\
0&\Omega_U^2&-R\\
0&-R&\Omega_W^2
\end{pmatrix}.
\]

The mixed Maxwell sub-block is positive iff

\[
\Delta=\Omega_U^2\Omega_W^2-R^2>0
\]

with \(\Omega_U^2,\Omega_W^2>0\).

The wall-coupling vector is

\[
g_{\rm int}=(c_B,g_U,g_W)^T.
\]

The Schur complement of the internal block is

\[
D_0
=
K-g_{\rm int}^TA_{\rm int}^{-1}g_{\rm int}.
\]

The script verifies exactly that

\[
D_0
=
K-\frac{c_B^2}{\varpi^2}
-
\frac{
g_U^2\Omega_W^2+2g_Ug_WR+g_W^2\Omega_U^2
}{\Delta}.
\]

It also verifies

\[
\det P
=
\varpi^2\Delta D_0.
\]

Therefore the positive-energy conservative gate is

\[
\boxed{
M>0,\quad \varpi^2>0,\quad \Delta>0,\quad D_0>0.
}
\]

This is the precise version of “near-softening is allowed only below the
instability boundary.”

---

## 3. Low-frequency conservative moments

Define

\[
Q=
g_U^2\Omega_W^2+2g_Ug_WR+g_W^2\Omega_U^2,
\qquad
S=\Omega_U^2+\Omega_W^2,
\]
\[
H=g_U^2+g_W^2.
\]

The BdG moments are

\[
B_0=\frac{c_B^2}{\varpi^2},
\qquad
B_2=\frac{c_B^2}{\varpi^4},
\qquad
B_4=\frac{c_B^2}{\varpi^6}.
\]

The Maxwell/mixed conservative moments are

\[
Z_0=\frac{Q}{\Delta},
\]
\[
Z_2=\frac{QS-H\Delta}{\Delta^2},
\]
\[
Z_4=
\frac{Q(S^2-\Delta)-SH\Delta}{\Delta^3}.
\]

The script verifies the matrix identities

\[
Z_0=g^TA^{-1}g,\qquad
Z_2=g^TA^{-2}g,\qquad
Z_4=g^TA^{-3}g,
\]

where

\[
A=
\begin{pmatrix}
\Omega_U^2&-R\\
-R&\Omega_W^2
\end{pmatrix},
\qquad
g=(g_U,g_W)^T.
\]

Because \(A>0\), all these \(Z\)-moments are nonnegative quadratic forms.

The conservative wall operator is

\[
D(\omega)
=
D_0+D_2\omega^2+D_4\omega^4+O(\omega^6),
\]

with

\[
D_0=K-B_0-Z_0,
\]
\[
D_2=-(M+B_2+Z_2),
\]
\[
D_4=-(B_4+Z_4).
\]

The script verifies that this matches the direct low-frequency expansion of the
full conservative self-energy through \(O(\omega^4)\).

---

## 4. Outgoing-port passivity gate

For the Maxwell/mixed outgoing transfer factor,

\[
N_0
=
\frac{(\Omega_U^2g_W+Rg_U)^2}{\Delta^2}.
\]

This is a perfect square over \(\Delta^2\), so

\[
N_0\ge0
\]

provided the internal block is nonsingular.

For the compact outgoing \(l=2\) port,

\[
\Gamma_2^{\rm port}
=
\frac{a^5}{27c_s^5}.
\]

The wall-level odd coefficient is therefore

\[
\gamma_{\rm wall}
=
N_0\Gamma_2^{\rm port}
=
N_0\frac{a^5}{27c_s^5}
\ge0.
\]

With the convention \(e^{-i\omega t}\), the wall operator contribution

\[
\delta D_{\rm odd}(\omega)
=
-i\gamma_{\rm wall}\omega^5
\]

corresponds in time domain to

\[
+\gamma_{\rm wall} q^{(5)}.
\]

The script verifies the identity

\[
\dot q\,q^{(5)}
=
\frac{d}{dt}
\left(
\dot q\,q^{(4)}-\ddot q\,q^{(3)}
\right)
+
(q^{(3)})^2.
\]

Thus, after the usual Schott-energy reshuffling, the outgoing \(l=2\) term
dissipates

\[
\gamma_{\rm wall}(q^{(3)})^2\ge0.
\]

So the outgoing branch is passive whenever \(N_0\ge0\) and
\(\Gamma_2^{\rm port}>0\).

---

## 5. Failure modes

The audit isolates four concrete failure modes.

### 5.1 Internal mixed-block instability

If

\[
\Delta\le0,
\]

then the \(U/W\) internal block is indefinite.  The Maxwell/mixed support pair is
not a positive-energy conservative subsystem.

### 5.2 Static wall softening instability

If

\[
D_0\le0,
\]

then the Schur-complement stiffness of the wall is nonpositive.  The apparent
normalization enhancement from \(P_0=N_0/D_0\) has crossed into instability.

### 5.3 Ghost/Krein contamination

If the reduced BdG coordinates are not projected onto positive-energy stable
normal modes, the kinetic or symplectic signature assumption behind this
Hamiltonian ledger fails.  Such modes require a separate full BdG signature
audit.

### 5.4 Dark outgoing port

If

\[
\Omega_U^2g_W+Rg_U=0,
\]

then

\[
N_0=0.
\]

The branch can still be conservative and stable, but it cannot transfer the
leading outgoing \(l=2\) port to the wall at this order.

---

## 6. Numerical sanity checks

The script includes one stable rational example with

\[
D_0>0,\qquad \Delta>0,
\]

and all generalized eigenvalues \(\omega^2\) positive.

It then lowers \(K\) while holding the internal data fixed.  The resulting
example has

\[
D_0<0,
\]

and one generalized \(\omega^2\) eigenvalue becomes negative, confirming that
the Schur-complement gate is the static instability boundary.

It also includes a mixed-block example with \(\Delta<0\), showing that the
localized Maxwell/mixed internal block itself becomes indefinite.

---

## 7. Carry-forward statement

The combined reduced stability theorem is:

\[
\boxed{
T>0,\quad A_{\rm int}>0,\quad
D_0=K-g_{\rm int}^TA_{\rm int}^{-1}g_{\rm int}>0
\quad\Longrightarrow\quad
H_{\rm cons}>0.
}
\]

In scalar one-lane form,

\[
\boxed{
M>0,\quad
\varpi^2>0,\quad
\Omega_U^2\Omega_W^2-R^2>0,\quad
K-\frac{c_B^2}{\varpi^2}
-
\frac{
g_U^2\Omega_W^2+2g_Ug_WR+g_W^2\Omega_U^2
}{
\Omega_U^2\Omega_W^2-R^2
}>0.
}
\]

The outgoing \(l=2\) port is passive if

\[
\boxed{
\gamma_{\rm wall}
=
\frac{a^5}{27c_s^5}
\frac{(\Omega_U^2g_W+Rg_U)^2}
{(\Omega_U^2\Omega_W^2-R^2)^2}
\ge0.
}
\]

So the normalization strategy may use near-softening, but only with

\[
D_0>0
\]

kept as a hard branch gate.
