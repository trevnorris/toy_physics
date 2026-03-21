# 2PN P0/P2 mouth-port operator rebuild: current result

## 1. What was shown

The solved **added** conservative 2PN comparable-mass cross block is not just a target-matched polynomial.
It admits an exact **body-local mouth-port decomposition**.

The reconstructed block is built from:

- the already-frozen **1PN dipole wake**, and
- a genuinely new **2PN port layer** consisting of one monopole-like scalar port plus the full real quadrupole multiplet.

So the minimal new conservative 2PN port content is

a carried-forward dipole sector plus a new
\[
P_0 \oplus P_2
\]
response layer.

---

## 2. Frozen 1PN wake is already a dipole-kernel overlap

Write the separation direction as \(\mathbf n = \hat z\), and decompose each velocity into

\[
\mathbf v_A = \mathbf u_A + d_A \mathbf n,
\qquad
\mathbf v_B = \mathbf u_B + d_B \mathbf n,
\]
with \(\mathbf u_{A,B}\cdot \mathbf n = 0\).

Then the already-frozen 1PN wake
\[
L^{\rm wake}_{1PN} = -\frac72 \,\mathbf v_A\!\cdot\!\mathbf v_B - \frac12 d_A d_B
\]
can be written exactly as
\[
L^{\rm wake}_{1PN}
=
-\Big[
\Pi^{(1\pm)}_A\!\cdot\!\Pi^{(1\pm)}_B
+
\Pi^{(10)}_A\Pi^{(10)}_B
\Big],
\]
with body-local dipole ports
\[
\Pi^{(1\pm)}_A = \sqrt{\frac72}\,\mathbf u_A,
\qquad
\Pi^{(10)}_A = 2 d_A,
\]
and similarly for body \(B\).

So the frozen 1PN wake is already a body-local dipole kernel: transverse \(m=\pm1\) plus longitudinal \(m=0\).

---

## 3. After removing universal leg dressing, the new quartic residual is exactly \(P_0\oplus P_2\)

Subtract the universal carried-forward piece
\[
\frac12 (v_A^2+v_B^2)L^{\rm wake}_{1PN}
\]
from the solved quartic 2PN cross target.
The remaining quartic residual is reproduced **exactly** by the following local ports.

### Monopole-like scalar port
\[
\Pi^{(0)}_A = \frac{\sqrt5}{2}\,v_A^2.
\]

### Quadrupole \(m=0\) scalar port
\[
\Pi^{(20)}_A = \frac12\bigl(3d_A^2-v_A^2\bigr).
\]

### Quadrupole \(m=\pm1\) real pair
\[
\Pi^{(21)}_A = \sqrt2\,d_A\,\mathbf u_A.
\]

### Quadrupole \(m=\pm2\) real pair
Using \(\mathbf u_A=(u_{Ax},u_{Ay})\), define
\[
\Pi^{(22)}_A = \frac{1}{2\sqrt2}
\begin{pmatrix}
 u_{Ax}^2-u_{Ay}^2 \\
 2u_{Ax}u_{Ay}
\end{pmatrix}.
\]

Then the quartic residual is exactly
\[
Q^{\rm new}_{2PN}
=
\Pi^{(0)}_A\Pi^{(0)}_B
+
\Pi^{(20)}_A\Pi^{(20)}_B
+
\Pi^{(21)}_A\!\cdot\!\Pi^{(21)}_B
+
\Pi^{(22)}_A\!\cdot\!\Pi^{(22)}_B.
\]

So the genuinely new quartic 2PN tensor sector is a **positive Gram overlap** of one real \(P_0\) channel plus the full real \(P_2\) multiplet.

That is the first direct algebraic bridge between the solved 2PN comparable-mass cross block and the inner-throat notes’ mouth-port language.

---

## 4. The old scalar \((T,L)\) tensor block diagonalizes exactly into \(P_0\) and \(P_{20}\)

Earlier, the quartic scalar tensor sector appeared in the \((T,L)\) basis with kernel
\[
K_{TL}=
\begin{pmatrix}
\frac32 & \frac34 \\
\frac34 & \frac94
\end{pmatrix}.
\]

With
\[
P_0 = T+L = v^2,
\qquad
P_{20} = -T + 2L = 3(v\cdot n)^2 - v^2,
\]
that kernel diagonalizes **exactly** to
\[
K_{P_0P_{20}}=
\begin{pmatrix}
\frac54 & 0 \\
0 & \frac14
\end{pmatrix}.
\]

So the old scalar tensor sector was already hiding the \(P_0\) and axisymmetric \(P_2\) port structure.

---

## 5. The quadratic \(G^2/r^2\) velocity block excites only the scalar ports \(P_0\) and \(P_{20}\)

The solved quadratic-velocity block can be written as
\[
-\frac{15}{4}(m_A+m_B)(\mathbf v_A\!\cdot\!\mathbf v_B)
+
 m_A\Bigl(2v_A^2 + \frac58(3d_A^2-v_A^2)\Bigr)
+
 m_B\Bigl(2v_B^2 + \frac58(3d_B^2-v_B^2)\Bigr).
\]

Equivalently, after factoring the overall pair prefactor and using
\[
U_A=\frac{Gm_B}{r},
\qquad
U_B=\frac{Gm_A}{r},
\]
this is
\[
-\frac{15}{4}(U_A+U_B)(\mathbf v_A\!\cdot\!\mathbf v_B)
+
U_A\left(\frac{4}{\sqrt5}\Pi^{(0)}_B + \frac54\Pi^{(20)}_B\right)
+
U_B\left(\frac{4}{\sqrt5}\Pi^{(0)}_A + \frac54\Pi^{(20)}_A\right).
\]

So the scalar-potential dressing drives **only** the axisymmetric scalar ports \(P_0\) and \(P_{20}\).

---

## 6. The static cross term is a pure monopole-potential overlap

The added static cross piece is simply
\[
\frac54 U_A U_B
\]
in the same overall cross normalization.
So it is naturally interpreted as a pure monopole-potential overlap.

---

## 7. Full constructive form of the added conservative 2PN cross block

Putting everything together, the solved added conservative 2PN cross block is exactly
\[
L^{\rm add}_{2PN,\,cross}
=
\frac{Gm_A m_B}{c^4 r}
\Bigg[
\frac12(v_A^2+v_B^2)L^{\rm wake}_{1PN}
+
\sum_{\lambda\in\{0,20,21c,21s,22c,22s\}}
\Pi^{(\lambda)}_A\Pi^{(\lambda)}_B
\]
\[
\qquad
-
\frac{15}{4}(U_A+U_B)(\mathbf v_A\!\cdot\!\mathbf v_B)
+
U_A\left(\frac{4}{\sqrt5}\Pi^{(0)}_B + \frac54\Pi^{(20)}_B\right)
+
U_B\left(\frac{4}{\sqrt5}\Pi^{(0)}_A + \frac54\Pi^{(20)}_A\right)
+
\frac54 U_AU_B
\Bigg].
\]

This matches the full solved added ADM-lifted 2PN cross target **exactly**.

---

## 8. What this means physically

The constructive picture is now much sharper.

1. **1PN:** the wake sector is already a dipole mouth-response kernel.
2. **2PN:** the genuinely new conservative cross sector does **not** require arbitrarily high channels.
3. It closes exactly with:
   - the carried-forward dipole wake,
   - a new **monopole + full quadrupole** response layer,
   - scalar-potential source terms on the scalar \(P_0/P_{20}\) channels,
   - and the static monopole-potential overlap.

So the next inner-throat target is no longer vague. The constructive PDE / DtN job is to reproduce this small port hierarchy:
\[
\text{dipole }(\ell=1)\quad\text{carried from 1PN},
\qquad
P_0\oplus P_2\quad\text{new at 2PN}.
\]

That is a much more concrete finish-line for the full 2PN derivation.
