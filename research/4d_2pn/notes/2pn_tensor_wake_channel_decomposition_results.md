
# 2PN tensor-wake channel decomposition: current result

## 1. What was solved

The solved added 2PN comparable-mass cross block from the ADM lift was

\[
L_{2,\text{cross}}=
\frac{Gm_A m_B}{c^4 r}
\Big[
-\frac74\,v_{AB}(v_A^2+v_B^2)
-\frac14\,v_{An}v_{Bn}(v_A^2+v_B^2)
+\frac{11}{8}\,v_A^2v_B^2
+\frac14\,v_{AB}^2
-\frac58\,(v_A^2v_{Bn}^2+v_B^2v_{An}^2)
+\frac32\,v_{AB}v_{An}v_{Bn}
+\frac38\,v_{An}^2v_{Bn}^2
\Big]
\]
\[
\qquad
+\frac{G^2m_A m_B}{c^4 r^2}
\Big[
\frac{11}{8}(m_A v_A^2+m_B v_B^2)
-\frac{15}{4}(m_A+m_B)v_{AB}
+\frac{15}{8}(m_A v_{An}^2+m_B v_{Bn}^2)
\Big]
+\frac{5G^3m_A^2m_B^2}{4c^4 r^3}.
\]

The new question was whether this block is just a target-matched polynomial, or whether it admits a clean **constructive wake/tensor decomposition**.

It does.

---

## 2. Universal vector-wake dressing already fixes two quartic coefficients

Using the frozen 1PN wake coefficients
\[
C_\parallel=-\frac72,
\qquad
C_L=-\frac12,
\]
the quartic pair terms
\[
-\frac74\,v_{AB}(v_A^2+v_B^2)
-\frac14\,v_{An}v_{Bn}(v_A^2+v_B^2)
\]
are exactly

\[
\frac12\,(v_A^2+v_B^2)\,
\bigl(C_\parallel v_{AB}+C_L v_{An}v_{Bn}\bigr).
\]

So the first two 2PN quartic coefficients are not new fits at all: they are a universal
\[
\sigma=\frac12
\]
leg-dressing of the frozen 1PN vector wake.

---

## 3. The remaining quartic residual is a rank-2 tensor-wake overlap

After subtracting that vector-dressed piece, the quartic residual is

\[
\frac{Gm_A m_B}{c^4 r}
\left[
\frac{11}{8}v_A^2v_B^2
+\frac14 v_{AB}^2
-\frac58(v_A^2v_{Bn}^2+v_B^2v_{An}^2)
+\frac32 v_{AB}v_{An}v_{Bn}
+\frac38 v_{An}^2v_{Bn}^2
\right].
\]

Define the pair-\(AB\) projector-channel invariants relative to \(\mathbf n_{AB}\):

\[
T_A = v_A^2-v_{An}^2,
\qquad
L_A = v_{An}^2,
\]
\[
T_B = v_B^2-v_{Bn}^2,
\qquad
L_B = v_{Bn}^2,
\]
\[
S_{AB}=(v_{AB}-v_{An}v_{Bn})^2-\frac12 T_A T_B,
\]
\[
M_{AB}=2(v_{AB}-v_{An}v_{Bn})v_{An}v_{Bn}.
\]

These have the interpretation:

- \(T\): transverse-trace scalar channel,
- \(S\): transverse-shear channel,
- \(M\): mixed transverse-longitudinal tensor channel,
- \(L\): pure longitudinal scalar channel.

Then the quartic residual is **exactly**

\[
\frac{Gm_A m_B}{c^4 r}
\Big[
\frac32\,T_A T_B
+\frac14\,S_{AB}
+1\cdot M_{AB}
+\frac34\,(T_A L_B + L_A T_B)
+\frac94\,L_A L_B
\Big].
\]

So the residual quartic sector is not arbitrary; it is a **minimal 5-channel rank-2 tensor-wake overlap** with coefficients

\[
k_{TT}=\frac32,\qquad
k_S=\frac14,\qquad
k_M=1,\qquad
k_{TL}=\frac34,\qquad
k_{LL}=\frac94.
\]

The chosen tensor basis has rank \(5\), so this is an exact minimal solve in that channel space.

---

## 4. Positive scalar-sector response matrix

The scalar \((T,L)\) subsector is governed by the matrix

\[
K_{\rm scalar}=
\begin{pmatrix}
\frac32 & \frac34\\[4pt]
\frac34 & \frac94
\end{pmatrix}.
\]

This is positive:
\[
\det K_{\rm scalar}=\frac{45}{16}>0,
\qquad
\operatorname{tr}K_{\rm scalar}=\frac{15}{4}>0,
\]
with eigenvalues
\[
\lambda_\pm=\frac{15\pm 3\sqrt5}{8}.
\]

So the scalar tensor sector is constructive/stable, not a sign-indefinite fit.

---

## 5. The \(G^2/r^2\) velocity block also decomposes cleanly

The solved quadratic-velocity block is

\[
\frac{G^2m_A m_B}{c^4 r^2}
\Big[
\frac{11}{8}(m_A v_A^2+m_B v_B^2)
-\frac{15}{4}(m_A+m_B)v_{AB}
+\frac{15}{8}(m_A v_{An}^2+m_B v_{Bn}^2)
\Big].
\]

This is exactly reproduced by

1. a **purely parallel local-potential dressing** of the frozen 1PN vector wake,
2. plus diagonal **transverse / longitudinal tensor-potential** channels.

Using
\[
U_A=\frac{Gm_B}{r},\qquad U_B=\frac{Gm_A}{r},
\]
the constructive form is

\[
\frac{Gm_A m_B}{c^4 r}
\Big[
\tau_\parallel (U_A+U_B) C_\parallel v_{AB}
+\beta_T (U_A T_B + U_B T_A)
+\beta_L (U_A L_B + U_B L_A)
\Big],
\]
with the exact solution

\[
\tau_\parallel=\frac{15}{14},
\qquad
\beta_T=\frac{11}{8},
\qquad
\beta_L=\frac{13}{4}.
\]

Notably, there is **no** local-potential longitudinal vector dressing: the solved block needs only the parallel one.

---

## 6. Full constructive decomposition of the added 2PN cross block

Putting everything together, the solved added cross block is exactly

\[
L_{2,\text{cross}}
=
L_{2,\text{vec,kin}}
+
L_{2,\text{tensor,quartic}}
+
L_{2,\text{quad,pot}}
+
L_{2,\text{static,cross}},
\]
where

\[
L_{2,\text{vec,kin}}
=
\frac{Gm_A m_B}{c^4 r}\,
\frac12\,(v_A^2+v_B^2)\,
\bigl(C_\parallel v_{AB}+C_L v_{An}v_{Bn}\bigr),
\]

\[
L_{2,\text{tensor,quartic}}
=
\frac{Gm_A m_B}{c^4 r}
\Big[
\frac32\,T_A T_B
+\frac14\,S_{AB}
+M_{AB}
+\frac34\,(T_A L_B+L_A T_B)
+\frac94\,L_A L_B
\Big],
\]

\[
L_{2,\text{quad,pot}}
=
\frac{Gm_A m_B}{c^4 r}
\Big[
\frac{15}{14}(U_A+U_B) C_\parallel v_{AB}
+\frac{11}{8}(U_A T_B + U_B T_A)
+\frac{13}{4}(U_A L_B + U_B L_A)
\Big],
\]

\[
L_{2,\text{static,cross}}
=
\frac{5G^3m_A^2m_B^2}{4c^4 r^3}.
\]

This reconstructs the solved ADM-lift cross block **exactly**.

---

## 7. First 3-body predictions from the local-potential pieces

For the pair \(AB\) inside a 3-body environment \(A,B,C\), with
\[
U_A^{(3)}=G\left(\frac{m_B}{r_{AB}}+\frac{m_C}{r_{AC}}\right),
\qquad
U_B^{(3)}=G\left(\frac{m_A}{r_{AB}}+\frac{m_C}{r_{BC}}\right),
\]
the constructive local-potential module already predicts the following \(3\)-body coefficients for the pair-\(AB\) kinematics:

\[
-\frac{15}{4}
\quad\text{on}\quad
\frac{G^2 m_A m_B m_C}{c^4}
\left[
\frac{v_{AB}}{r_{AB}r_{AC}}
+
\frac{v_{AB}}{r_{AB}r_{BC}}
\right],
\]

\[
\frac{11}{8}
\quad\text{on}\quad
\frac{G^2 m_A m_B m_C}{c^4}
\left[
\frac{v_B^2}{r_{AB}r_{AC}}
+
\frac{v_A^2}{r_{AB}r_{BC}}
\right],
\]

\[
\frac{15}{8}
\quad\text{on}\quad
\frac{G^2 m_A m_B m_C}{c^4}
\left[
\frac{v_{Bn}^2}{r_{AB}r_{AC}}
+
\frac{v_{An}^2}{r_{AB}r_{BC}}
\right].
\]

So the constructive module is already making nontrivial \(N\)-body predictions beyond the matched 2-body sector.

---

## 8. Interpretation

This is a real step forward.

The 2PN cross sector is no longer just “some polynomial found by ADM matching.” It now has a concrete channel decomposition:

- **frozen 1PN vector wake** with universal kinematic dressing,
- **new rank-2 tensor wake** with a small positive projector-channel kernel,
- **local-potential dressing** of the vector/tensor channels,
- **scalar static cross term**.

That is exactly the kind of structure a genuine `tensor_wake_2pn_rebuild.wl` should aim to derive from the inner-throat / DtN side.

The hardest remaining question is now much sharper:

> can the inner-throat PDE / DtN machinery reproduce these channel strengths
> \[
> \sigma=\frac12,\quad
> k_{TT}=\frac32,\quad
> k_S=\frac14,\quad
> k_M=1,\quad
> k_{TL}=\frac34,\quad
> k_{LL}=\frac94,\quad
> \tau_\parallel=\frac{15}{14},\quad
> \beta_T=\frac{11}{8},\quad
> \beta_L=\frac{13}{4},
> \]
> rather than having them inserted algebraically?

That is a much better-defined nut than the one we started with.
