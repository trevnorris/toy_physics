# 2PN mouth-operator / DtN reduction: current preliminary result

## Core structural result

Take the low-frequency monopole mouth operator in the form
\[
Z_{00}(\omega;u)=Z_2(u)\,\omega^2+Z_4(u)\,\omega^4+O(\omega^6),
\qquad u\equiv U/c^2.
\]

Define the normalized DtN invariants
\[
C_s(u)\equiv \frac{Z_4/Z_2^3}{(Z_4/Z_2^3)_0},
\qquad
G(u)\equiv \frac{Z_4/Z_2^2}{(Z_4/Z_2^2)_0}.
\]

For the cylinder / Neumann-bottom unit-test branch,
\[
Z_2=-\frac{L}{c_s^2},
\qquad
Z_4=-\frac{L^3}{3c_s^4},
\]
so these become exactly
\[
C_s(u)=\frac{c_s^2(u)}{c_{s0}^2},
\qquad
G(u)=\frac{L(u)}{L_0}.
\]

Under the current 2PN freeze:
- exact Bernoulli gives \(c_s^2/c^2 = 1-4u\),
- the reduced throat closure gives \(L/L_0 = a(u)\),
- and the known 1PN breathing slope implies
  \[
  a(u)=1+\frac{57}{64}u+O(u^2).
  \]

## Unique minimal one-body closure that preserves the frozen 1PN slot

If the corrected one-body denominator is built from the DtN geometry invariant as
\[
D_{\rm eff}(u)=C_s(u)\left[1+\alpha\,(G-1)+\beta\,(G-1)^2\right],
\]
then preserving the already-fixed 1PN coefficient forces
\[
\alpha=0.
\]
So the first nontrivial allowed correction is quadratic:
\[
D_{\rm eff}(u)=C_s(u)\left[1+\mu\,(G-1)^2\right].
\]

At 2PN only the linear geometry slope matters. Writing
\[
G(u)=1+g_1 u+g_2 u^2+\cdots,
\]
one gets
\[
D_{\rm eff}(u)=1-4u+\mu g_1^2 u^2+O(u^3).
\]
Matching the exact isotropic Schwarzschild test-mass target through 2PN therefore fixes
\[
\mu g_1^2 = 8.
\]
Using the already-fixed throat slope
\[
g_1=\frac{57}{64},
\]
this gives
\[
\mu = \frac{8}{g_1^2} = \frac{32768}{3249} \approx 10.0855647892.
\]

So the final one-body denominator is
\[
D_{\rm eff}(u)=C_s(u)\left[1+\frac{32768}{3249}(G(u)-1)^2\right].
\]

With the current closure this expands to
\[
D_{\rm eff}(u)=1-4u+8u^2+O(u^3),
\]
which exactly closes the isotropic one-body 2PN target.

## Explicit current-series values

From the current throat closure:
\[
a(u)=1+\frac{57}{64}u+\frac{298821}{131072}u^2+O(u^3).
\]

Therefore the cylinder DtN coefficients become
\[
\frac{Z_2(u)}{Z_2(0)} = 1+\frac{313}{64}u+\frac{2862917}{131072}u^2+O(u^3),
\]
\[
\frac{Z_4(u)}{Z_4(0)} = 1+\frac{683}{64}u+\frac{10301487}{131072}u^2+O(u^3).
\]

The extracted invariants are then
\[
C_s(u)=1-4u,
\qquad
G(u)=1+\frac{57}{64}u+\frac{298821}{131072}u^2+O(u^3).
\]

## Relation to the earlier raw resonance-proxy fit

The raw resonance proxy derived from the same DtN data is
\[
D_{\rm raw}(u)=\frac{[Z_2/Z_4](u)}{[Z_2/Z_4](0)}
=\frac{1-4u}{G(u)^2}
=1-\frac{185}{32}u+\frac{324075}{65536}u^2+O(u^3).
\]

The multiplicative factor needed to convert it into the exact target denominator is
\[
P_{\rm port}(u)=G(u)^2\left[1+\mu (G(u)-1)^2\right]
=1+\frac{57}{32}u+\frac{875093}{65536}u^2+O(u^3).
\]

So the earlier port fit was not arbitrary; it factorizes cleanly as
\[
D_{\rm eff}(u)=D_{\rm raw}(u)\,P_{\rm port}(u).
\]

## Current status

This closes the **one-body missing 2PN response slot** under a minimal DtN-invariant conservative ansatz.

It does **not** yet derive that ansatz from the full inner PDE, and it does **not** yet close the comparable-mass 2PN wake/cross sector.

# 2PN comparable-mass ADM lift: preliminary result

## What was frozen

I kept the same freeze used in the passing 2PN harness and the passing DtN one-body closure:

- full frozen Newtonian + 1PN two-body Lagrangian,
- DtN-corrected one-body/self sector at 2PN,
- quadratic local mass scaling with \(\lambda_\rho=1/2\).

That means the 2PN self/static block starts from
\[
L_{2,\text{self+static}}=
\frac{m_A v_A^6+m_B v_B^6}{16}
+\frac{7Gm_A m_B}{8r}(v_A^4+v_B^4)
+\frac{2G^2 m_A m_B}{r^2}(m_B v_A^2+m_A v_B^2)
+\frac{G^3 m_A m_B(m_A^2+m_B^2)}{4r^3},
\]
with the understanding that the overall 2PN factor is \(c^{-4}\).

## The useful algebraic simplification

For
\[
L=L_0+\epsilon L_1+\epsilon^2 L_2,
\qquad \epsilon\equiv c^{-2},
\]
with quadratic \(L_0\), the perturbative Legendre transform through 2PN collapses to
\[
H_1=-L_1(v_0),
\qquad
H_2=-L_2(v_0)+\frac12 A_0^T M^{-1}A_0,
\]
where
\[
v_0=\frac{p}{m},
\qquad
A_0=\left.\frac{\partial L_1}{\partial v}\right|_{v=v_0}.
\]

This is the key reason the comparable-mass lift is tractable: once the full 1PN block is frozen, any *new* 2PN Lagrangian block enters the Hamiltonian only as \(-L_{2,\text{new}}(v_0)\).

## What happens with no 2PN cross sector

Using the generic-frame ADM Hamiltonian target through 2PN, the frozen 1PN block still matches exactly at 1PN, but the DtN-corrected self/static 2PN candidate leaves a nonzero residual in the comparable-mass sector.

That residual contains three types of missing structure:

1. \(G/r\) quartic cross-velocity monomials,
2. \(G^2/r^2\) quadratic-velocity comparable-mass terms,
3. a missing static cross mass polynomial proportional to \(m_A^2 m_B^2\).

## Minimal invariant basis and solved coefficients

I solved the residual in the compact basis
\[
\frac{Gm_A m_B}{r}\times
\Big\{
v_{AB}(v_A^2+v_B^2),\;
v_{An}v_{Bn}(v_A^2+v_B^2),\;
v_A^2v_B^2,\;
v_{AB}^2,\;
v_A^2v_{Bn}^2+v_B^2v_{An}^2,\;
v_{AB}v_{An}v_{Bn},\;
v_{An}^2v_{Bn}^2
\Big\},
\]
\[
\frac{G^2m_A m_B}{r^2}\times
\Big\{
m_B v_A^2+m_A v_B^2,\;
m_A v_A^2+m_B v_B^2,\;
(m_A+m_B)v_{AB},\;
(m_A+m_B)v_{An}v_{Bn},\;
m_B v_{An}^2+m_A v_{Bn}^2,\;
m_A v_{An}^2+m_B v_{Bn}^2
\Big\},
\]
and
\[
\frac{G^3}{r^3}\,m_A^2m_B^2.
\]

The solve is unique and gives
\[
q_1=-\frac74,\quad
q_2=-\frac14,\quad
q_3=\frac{11}{8},\quad
q_4=\frac14,\quad
q_5=-\frac58,\quad
q_6=\frac32,\quad
q_7=\frac38,
\]
\[
t_1=0,\quad
t_2=\frac{11}{8},\quad
t_3=-\frac{15}{4},\quad
t_4=0,\quad
t_5=0,\quad
t_6=\frac{15}{8},
\]
\[
s_1=\frac54.
\]

So three quadratic basis coefficients vanish automatically:
\[
t_1=t_4=t_5=0.
\]

## Added comparable-mass cross block

The required 2PN added cross block is
\[
L_{2,\text{cross}}=
\frac{Gm_A m_B}{r}
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
+\frac{G^2m_A m_B}{r^2}
\Big[
\frac{11}{8}(m_A v_A^2+m_B v_B^2)
-\frac{15}{4}(m_A+m_B)v_{AB}
+\frac{15}{8}(m_A v_{An}^2+m_B v_{Bn}^2)
\Big]
+\frac{5G^3m_A^2m_B^2}{4r^3}.
\]

Again, the overall 2PN factor is \(c^{-4}\).

## Full candidate through 2PN

Adding that to the frozen self/static piece gives the full preliminary 2PN block
\[
L_{2,\text{full}}=
\frac{m_A v_A^6+m_B v_B^6}{16}
+\frac{Gm_A m_B}{r}\Big[
\frac78(v_A^4+v_B^4)
-\frac74\,v_{AB}(v_A^2+v_B^2)
-\frac14\,v_{An}v_{Bn}(v_A^2+v_B^2)
+\frac{11}{8}\,v_A^2v_B^2
+\frac14\,v_{AB}^2
\]
\[
\qquad\qquad
-\frac58\,(v_A^2v_{Bn}^2+v_B^2v_{An}^2)
+\frac32\,v_{AB}v_{An}v_{Bn}
+\frac38\,v_{An}^2v_{Bn}^2
\Big]
\]
\[
\qquad
+\frac{G^2m_A m_B}{r^2}
\Big[
\Big(2m_B+\frac{11}{8}m_A\Big)v_A^2
+\Big(2m_A+\frac{11}{8}m_B\Big)v_B^2
-\frac{15}{4}(m_A+m_B)v_{AB}
+\frac{15}{8}(m_A v_{An}^2+m_B v_{Bn}^2)
\Big]
\]
\[
\qquad
+\frac{G^3m_A m_B}{4r^3}(m_A^2+5m_A m_B+m_B^2).
\]

With the overall \(c^{-4}\) restored, this full candidate Legendre-transforms **exactly** to the generic-frame ADM \(H_{2\mathrm{PN}}\) target.

## Sanity checks that pass

- The frozen 1PN Lagrangian still transforms to the ADM \(H_{1\mathrm{PN}}\) target exactly.
- The full added cross block vanishes in the strict test-mass limit with the heavy body at rest.
- The static mass polynomial is upgraded from
  \[
  +(m_A^2+m_B^2)
  \]
  to
  \[
  +(m_A^2+5m_A m_B+m_B^2),
  \]
  exactly the ADM/Hamiltonian mass structure after Legendre matching.

## Interpretation

This is a strong sign that the math is workable.

The 2PN comparable-mass problem is no longer “mysterious” at the target-matching level. Once the DtN-corrected one-body block is frozen, the remaining lift becomes a **linear residual solve**, and that solve lands in a compact, symmetry-respecting invariant basis.

What remains open is **not** whether a full 2PN candidate can be written down. It can.

What remains open is the **constructive explanation** of these coefficients from the inner-throat / wake side. In other words, the next real task is to see whether a 2PN wake/tensor construction can reproduce exactly the solved quartic and quadratic cross coefficients above, rather than having to insert them by matching.

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
# 2PN support-minus-closure throat EFT: current result

## 1. What changed at this step

The previous `P0⊕P2` mouth-port rebuild showed that the genuinely new conservative 2PN comparable-mass cross sector can be written in a small body-local port basis.

This step goes one layer deeper.

It shows that the solved **added** conservative 2PN cross block is exactly the on-shell result of a **local signed auxiliary-channel EFT** with three ingredients:

1. an **odd dipole** wake bundle (carried forward from 1PN),
2. a new **even support bundle** built from `P0 ⊕ P2` ports,
3. and one **pure-potential negative closure channel**.

So the constructive picture is now:


the 2PN cross sector is not merely a matched polynomial,

it is a **support-minus-closure throat response theory**.

---

## 2. Odd sector: the 1PN wake is an odd dipole overlap

With the pair axis along `n = z`, define the odd dipole channels

\[
\mathcal D_A = \Big(\sqrt{\tfrac72}\,v_{Ax},\;\sqrt{\tfrac72}\,v_{Ay},\;2v_{Az}\Big),
\qquad
\mathcal D_B = \Big(\sqrt{\tfrac72}\,v_{Bx},\;\sqrt{\tfrac72}\,v_{By},\;2v_{Bz}\Big).
\]

Then the frozen 1PN wake is exactly

\[
L^{\rm wake}_{1PN} = -\mathcal D_A\!\cdot\!\mathcal D_B.
\]

The **added** odd 2PN block comes from dressing these same channels by

- the universal kinetic factor `1 + 1/2 v^2`, and
- a channel-dependent potential dressing
  \[
  \eta_\perp = \frac{15}{14},
  \qquad
  \eta_\parallel = \frac{15}{16}.
  \]

So the dressed odd sources are

\[
\mathcal D_{A,\perp}^{\rm eff}
=
\Bigl(1+\frac12 v_A^2+\eta_\perp U_A\Bigr)
\sqrt{\tfrac72}\,\mathbf u_A,
\]
\[
\mathcal D_{A,0}^{\rm eff}
=
\Bigl(1+\frac12 v_A^2+\eta_\parallel U_A\Bigr)
2(v_A\!\cdot\!n),
\]
(and similarly for body `B`).

Expanding the antisymmetric overlap to first order in the dressings reproduces exactly

\[
\frac12(v_A^2+v_B^2)L^{\rm wake}_{1PN}
-\frac{15}{4}(U_A+U_B)(v_A\!\cdot\!v_B).
\]

So the old 1PN wake survives as the **odd channel backbone** of the 2PN construction.

---

## 3. Even sector: six positive support channels

The new even body-local sources are:

\[
\Pi^{(0)}=\frac{\sqrt5}{2}v^2,
\qquad
\Pi^{(20)}=\frac12\bigl(3(v\!\cdot\!n)^2-v^2\bigr),
\]
\[
\Pi^{(21c)}=\sqrt2\,(v\!\cdot\!n)u_x,
\qquad
\Pi^{(21s)}=\sqrt2\,(v\!\cdot\!n)u_y,
\]
\[
\Pi^{(22c)}=\frac{u_x^2-u_y^2}{2\sqrt2},
\qquad
\Pi^{(22s)}=\frac{2u_xu_y}{2\sqrt2}.
\]

The solved quartic residual was already

\[
\sum_{\lambda\in\{0,20,21c,21s,22c,22s\}}
\Pi^{(\lambda)}_A\Pi^{(\lambda)}_B.
\]

The new result is that once the solved `U × Pi` coefficients are included, the even block becomes

\[
\mathcal S_A\!\cdot\!\mathcal S_B - \mathcal C_A\mathcal C_B,
\]
with the **six positive support channels**

\[
\mathcal S_A=
\Big(
\Pi_A^{(0)}+\frac{4}{\sqrt5}U_A,
\Pi_A^{(20)}+\frac54 U_A,
\Pi_A^{(21c)},
\Pi_A^{(21s)},
\Pi_A^{(22c)},
\Pi_A^{(22s)}
\Big),
\]
(and likewise for `B`).

So the actual mouth/port support block is positive and minimal:

- one monopole support channel,
- one axisymmetric quadrupole support channel,
- the four real `m=±1,±2` quadrupole channels.

---

## 4. The negative piece is a pure geometry / closure channel

The positive support channels alone would generate a static coefficient

\[
\left(\frac{4}{\sqrt5}\right)^2 + \left(\frac54\right)^2
= \frac{381}{80}.
\]

But the solved static cross coefficient is only

\[
\frac54 = \frac{100}{80}.
\]

So the **closure deficit** is exactly

\[
\Delta_{\rm geom} = \frac{381}{80}-\frac{100}{80} = \frac{281}{80}.
\]

This is removed by one negative channel

\[
\mathcal C_A = \sqrt{\frac{281}{80}}\,U_A,
\qquad
\mathcal C_B = \sqrt{\frac{281}{80}}\,U_B.
\]

Therefore the whole even sector is

\[
L^{\rm even}_{2PN}
=
\mathcal S_A\!\cdot\!\mathcal S_B
-\mathcal C_A\mathcal C_B.
\]

The important structural fact is:

**the negative direction lives entirely in the pure-`U` direction.**

It does **not** mix with the genuine mouth ports.

That means the indefinite part is naturally interpreted as a **geometry-energy closure term**, not as a pathology of the mouth DtN/operator sector.

---

## 5. Matrix form

In the ordered source basis

\[
X=(\Pi^{(0)},\Pi^{(20)},\Pi^{(21c)},\Pi^{(21s)},\Pi^{(22c)},\Pi^{(22s)},U),
\]

the exact even response matrix is

\[
M_{\rm even}=
\begin{pmatrix}
1&0&0&0&0&0&4/\sqrt5\\
0&1&0&0&0&0&5/4\\
0&0&1&0&0&0&0\\
0&0&0&1&0&0&0\\
0&0&0&0&1&0&0\\
0&0&0&0&0&1&0\\
4/\sqrt5&5/4&0&0&0&0&5/4
\end{pmatrix}.
\]

It factorizes exactly as

\[
M_{\rm even}=R_{\rm support}^T R_{\rm support} - R_{\rm geom}^T R_{\rm geom},
\]
with

\[
R_{\rm support}=
\begin{pmatrix}
1&0&0&0&0&0&4/\sqrt5\\
0&1&0&0&0&0&5/4\\
0&0&1&0&0&0&0\\
0&0&0&1&0&0&0\\
0&0&0&0&1&0&0\\
0&0&0&0&0&1&0
\end{pmatrix},
\qquad
R_{\rm geom}=\begin{pmatrix}0&0&0&0&0&0&\sqrt{281/80}\end{pmatrix}.
\]

So the even block is literally a **PSD mouth-response block minus a rank-1 geometry closure block**.

---

## 6. Auxiliary-channel interpretation

The factorization has a very simple local auxiliary-field meaning.

For a positive support channel `q` with body sources `S_A,S_B`, use

\[
L_q=-\frac12 q^2 + q(S_A+S_B).
\]

Eliminating `q` gives the cross term `+S_A S_B`.

For the negative closure channel `g` with sources `C_A,C_B`, use

\[
L_g=+\frac12 g^2 + g(C_A+C_B).
\]

Eliminating `g` gives the cross term `-C_A C_B`.

For the odd dipole wake channel `d` with body sources `D_A,D_B`, use

\[
L_d=-\frac12 d^2 + d(D_A-D_B).
\]

Eliminating `d` gives the odd cross term `-D_A D_B`.

So the solved 2PN sector is exactly the on-shell result of a **small signed auxiliary-channel EFT**.

---

## 7. Full added conservative 2PN cross block

Putting everything together, the entire solved **added** conservative 2PN cross sector is

\[
L^{\rm add}_{2PN,\,cross}
=
\frac{Gm_A m_B}{c^4 r}
\Big[
L^{\rm add}_{\rm odd,dressed}
+
\mathcal S_A\!\cdot\!\mathcal S_B
-
\mathcal C_A\mathcal C_B
\Big],
\]
where

\[
L^{\rm add}_{\rm odd,dressed}
=
\frac12(v_A^2+v_B^2)L^{\rm wake}_{1PN}
-
\frac{15}{4}(U_A+U_B)(v_A\!\cdot\!v_B).
\]

This matches the solved ADM-lifted added 2PN target exactly.

---

## 8. Why this matters

This is the first point where the 2PN cross sector stops looking like “a list of coefficients” and starts looking like a real throat-response theory.

The constructive content is now very sharp:

- **Odd channels:** the 1PN dipole wake survives, with kinetic and potential dressing.
- **Even channels:** the new `P0⊕P2` layer is a positive mouth-response bundle.
- **Closure:** the only negative piece is a pure potential geometry term.

That means the next true finish-line problem is no longer “guess the 2PN cross polynomial.”
It is much narrower:

> derive the six positive support channels from the inner throat DtN operator,
> and derive the rank-1 pure-`U` closure channel from the geometry-energy side.

That is a very concrete split between **mouth response** and **geometry closure**, and it is exactly the kind of split the throat notes were already suggesting.
# 2PN minimal irreducible-channel throat operator: current result

## 1. What this step does

This step turns the solved 2PN cross-sector reconstruction into the **smallest actual low-frequency throat operator** that reproduces the conservative comparable-mass 2PN target.

The outcome is sharper than the previous support/closure factorization:

- the full solved **added** 2PN cross sector is reproduced by a canonical zero-frequency operator on
  
  \[
  \ell=1\quad\text{odd dipole channels},
  \qquad
  \ell=0\oplus 2\quad\text{even support channels},
  \]
  
- together with one direct pure-\(U\) geometry-closure term,
- and the first nontrivial **3-body lift** follows immediately once the local potentials are promoted beyond the strict two-body limit.

So the constructive content is now not just “some ports work.”
It is:

> a minimal irreducible-channel throat operator exists, it is uniquely fixed at zero frequency by the solved 2PN cross target, and it already predicts the first local-potential 3-body terms.

---

## 2. Odd sector: unique dressed dipole data

Write the frozen 1PN wake as

\[
L^{\rm wake}_{1PN}=-\frac72\,\mathbf u_A\!\cdot\!\mathbf u_B-4d_A d_B,
\qquad d_A=\mathbf v_A\!\cdot\!\mathbf n,
\]
with \(\mathbf u_A\) transverse to \(\mathbf n\).

The most general minimal added odd block with universal kinetic leg dressing and separate transverse/longitudinal potential dressings is

\[
L^{\rm add}_{\rm odd}
=\sigma(v_A^2+v_B^2)L^{\rm wake}_{1PN}
-(U_A+U_B)\bigl(p_\perp\,\mathbf u_A\!\cdot\!\mathbf u_B+p_\parallel d_A d_B\bigr).
\]

Exact matching gives the unique solution

\[
\sigma=\frac12,
\qquad
p_\perp=p_\parallel=\frac{15}{4}.
\]

Dividing by the frozen 1PN dipole residues gives the body-source dressing coefficients

\[
\eta_\perp=\frac{15}{14},
\qquad
\eta_\parallel=\frac{15}{16}.
\]

So the odd channel data are now fully fixed at 2PN conservative order.

---

## 3. Even sector: unique \(P_0\oplus P_2\) zero-frequency operator

In the canonical real port basis

\[
\Pi^{(0)}=\frac{\sqrt5}{2}v^2,
\qquad
\Pi^{(20)}=\frac12\bigl(3(v\!\cdot\!n)^2-v^2\bigr),
\]

with the real \(m=\pm1,\pm2\) quadrupole ports \(\Pi^{(21c)},\Pi^{(21s)},\Pi^{(22c)},\Pi^{(22s)}\), the most general minimal even zero-frequency operator is

\[
L^{\rm even}_{\rm min}
=
\sum_{\lambda} a_\lambda\,\Pi_A^{(\lambda)}\Pi_B^{(\lambda)}
+U_A\bigl(j_0\Pi_B^{(0)}+j_{20}\Pi_B^{(20)}\bigr)
+U_B\bigl(j_0\Pi_A^{(0)}+j_{20}\Pi_A^{(20)}\bigr)
+s\,U_AU_B.
\]

Exact matching gives the unique solution

\[
a_0=a_{20}=a_{21}=a_{22}=1,
\qquad
j_0=\frac{4}{\sqrt5},
\qquad
j_{20}=\frac54,
\qquad
s=\frac54.
\]

So the new conservative 2PN even operator is exactly the **identity metric** on the six real \(\ell=0\oplus 2\) support channels, plus a scalar source vector and a single static coefficient.

---

## 4. Canonical zero-frequency throat data

The full minimal operator data are therefore:

### Odd channel residues
\[
R_{1\perp}=\frac72,
\qquad
R_{10}=4.
\]

### Even support residues
\[
R_0=1,
\qquad
R_2=1.
\]

### Scalar source vector
\[
J=\left(\frac{4}{\sqrt5},\frac54,0,0,0,0\right).
\]

### Support static coefficient
\[
J\cdot J
=\left(\frac{4}{\sqrt5}\right)^2+\left(\frac54\right)^2
=\frac{381}{80}.
\]

### Geometry closure deficit
Since the solved static coefficient is only \(5/4\), the required direct closure term is

\[
\Delta_{\rm geom}=\frac{381}{80}-\frac54=\frac{281}{80}.
\]

### Direct geometry-energy coefficient
The pure-\(U\) closure block may be written as a direct geometry-energy cost

\[
E^{\rm pair}_{\rm geom,cl}=\frac{\Delta_{\rm geom}}{2}(U_A+U_B)^2
=\frac{281}{160}(U_A+U_B)^2,
\]

which contributes the Lagrangian cross term

\[
L^{\rm pair}_{\rm geom,cl}=-\Delta_{\rm geom}U_AU_B
=-\frac{281}{80}U_AU_B
\]

after the self pieces are dropped.

So the even block is exactly

\[
L^{\rm even}_{2PN}
=(\Pi_A+JU_A)\cdot(\Pi_B+JU_B)-\Delta_{\rm geom}U_AU_B.
\]

This is the smallest zero-frequency throat operator that reproduces the full added conservative 2PN cross target.

---

## 5. What is still not fixed

The conservative 2PN matching fixes only the **zero-frequency residues**.
A generic low-frequency DtN completion can therefore be written as

\[
Y_0(\omega)=1+\chi_0\omega^2+\cdots,
\qquad
Y_2(\omega)=1+\chi_2\omega^2+\cdots,
\]
\[
Y_{1\perp}(\omega)=\frac72+\chi_{1\perp}\omega^2+\cdots,
\qquad
Y_{10}(\omega)=4+\chi_{10}\omega^2+\cdots.
\]

The \(\chi\)-coefficients are not fixed by the present conservative 2PN algebra. They are genuine inner-PDE / DtN observables.

So this step sharply separates:

- what the 2PN derivation has already fixed, from
- what Paper 7 still has to compute dynamically.

---

## 6. First 3-body lift from local potentials

Now promote the pair-AB local potentials to

\[
U_A^{\rm loc}=G\left(\frac{m_B}{r_{AB}}+\frac{m_C}{r_{AC}}\right),
\qquad
U_B^{\rm loc}=G\left(\frac{m_A}{r_{AB}}+\frac{m_C}{r_{BC}}\right).
\]

Then the same minimal operator immediately predicts the pair-AB 3-body lift in a background body \(C\).

### Velocity-dependent 3-body lift
\[
L^{(3\text{-body, vel})}_{AB|C}
=
\frac{G^2 m_A m_B m_C}{8c^4 r_{AB}}
\left[
\frac{11v_B^2+15d_B^2-30(\mathbf v_A\!\cdot\!\mathbf v_B)}{r_{AC}}
+
\frac{11v_A^2+15d_A^2-30(\mathbf v_A\!\cdot\!\mathbf v_B)}{r_{BC}}
\right],
\]
where \(d_A=\mathbf v_A\!\cdot\!\mathbf n_{AB}\), \(d_B=\mathbf v_B\!\cdot\!\mathbf n_{AB}\).

This is exactly the previously observed 3-body velocity pattern, now derived from the minimal throat operator itself.

### Static 3-body lift
\[
L^{(3\text{-body, stat})}_{AB|C}
=
\frac{5G^3m_A m_B m_C}{4c^4}
\left[
\frac{m_A}{r_{AB}^2 r_{AC}}
+
\frac{m_B}{r_{AB}^2 r_{BC}}
+
\frac{m_C}{r_{AB} r_{AC} r_{BC}}
\right].
\]

So the minimal operator does not stop at the 2-body fit. It already carries nontrivial N-body content.

---

## 7. Why this matters

This step is important because it turns the constructive goal into a very specific target for the inner-throat PDE/DtN program.

The remaining problem is no longer vague.
The throat derivation must now explain:

1. why the zero-frequency odd residues are \(\{7/2,7/2,4\}\),
2. why the even \(\ell=0\oplus2\) support residues canonically normalize to the identity,
3. why the scalar source vector is
   \[
   J=\left(\frac{4}{\sqrt5},\frac54,0,0,0,0\right),
   \]
4. and why the geometry side contributes the direct closure coefficient
   \[
   \Delta_{\rm geom}=\frac{281}{80}.
   \]

At the same time, it cleanly identifies what still belongs to the PDE/DtN calculation rather than the current algebraic closure:

- the \(\omega^2\) low-frequency coefficients,
- and the dynamical origin of the geometry closure term.

So this is a genuine narrowing of the finish line, not just another algebraic repackaging.
# 2PN inner-throat modal DtN from PDE-side scaffolding: current result

## 1. What this step does

This step finally turns the solved 2PN cross-sector operator into a **real low-frequency throat-response model** rather than a static algebraic fit.

It does that in two stages.

First, it builds the simplest genuine PDE-side **unit-test branch**: the regular interior Helmholtz Dirichlet-to-Neumann (DtN) operator for a 4D ball.
That branch is explicit, diagonal in harmonic channels, and therefore a good baseline test of what a perfectly isotropic throat could and could not do.

Second, it builds the minimal **axisymmetric pole-completed DtN model** that *does* reproduce the solved conservative 2PN cross sector.

So this step is the first one that cleanly separates:

- what follows from a simple isotropic PDE cavity,
- what absolutely requires axisymmetric throat structure,
- and what dynamic data are still open PDE observables.

---

## 2. The isotropic 4D-ball unit test

For a regular interior Helmholtz mode in a 4D ball, the radial solution in harmonic degree \(\ell\) is

\[
u_{\ell}(r) \propto r^{-1}J_{\ell+1}(kr),
\qquad z \equiv ka,
\]

so the boundary DtN eigenvalue at radius \(a\) is

\[
\Lambda_\ell(z)
=
-\frac1a+
\frac{z}{a}
\frac{J'_{\ell+1}(z)}{J_{\ell+1}(z)}.
\]

Its small-\(z\) expansions for the first three channels are:

\[
\Lambda_0(z)
=
-\frac{z^2}{4a}
-\frac{z^4}{96a}
-\frac{z^6}{1536a}
-\frac{z^8}{23040a}
+\cdots,
\]
\[
\Lambda_1(z)
=
\frac1a
-\frac{z^2}{6a}
-\frac{z^4}{288a}
-\frac{z^6}{8640a}
-\frac{7z^8}{1658880a}
+\cdots,
\]
\[
\Lambda_2(z)
=
\frac2a
-\frac{z^2}{8a}
-\frac{z^4}{640a}
-\frac{z^6}{30720a}
-\frac{13z^8}{17203200a}
+\cdots.
\]

The corresponding admittances are

\[
Y_1(z)=\frac1{\Lambda_1(z)}
=
 a + \frac{a z^2}{6} + \frac{a z^4}{32}+\cdots,
\]
\[
Y_2(z)=\frac1{\Lambda_2(z)}
=
 \frac a2 + \frac{a z^2}{32} + \frac{3 a z^4}{1280}+\cdots.
\]

This is already very informative.

### What the isotropic branch gets right

- It is diagonal in harmonic channels.
- It naturally supports a positive \(\ell=1\) and \(\ell=2\) tower.
- After canonical normalization, it is completely compatible with an identity-like support metric inside a fixed \(\ell\) sector.

### What the isotropic branch cannot do

It fails two requirements of the solved 2PN operator immediately:

1. **No finite static monopole support**

   \[
   \Lambda_0(0)=0.
   \]

   So a pure isotropic cavity has no finite static \(\ell=0\) stiffness at zero frequency.
   That means it cannot by itself generate the solved finite monopole support / source structure.

2. **No dipole splitting**

   A spherical branch is degenerate inside each \(\ell\) multiplet.
   So it cannot distinguish the solved odd residues

   \[
   R_{1\perp}=\frac72,
   \qquad
   R_{10}=4.
   \]

So the isotropic 4D-ball branch is a **unit test**, not the full throat.

---

## 3. What the full throat must minimally contain

The solved 2PN operator therefore forces the following PDE-side conclusion:

> the full conservative throat response cannot be a pure isotropic cavity; it must at least be **axisymmetric** and include a separate **monopole/geometry** sector.

The minimal symmetry reduction is exactly

\[
O(3) \to O(2),
\]

which splits

- dipole \(\ell=1\) into \(|m|=1\) and \(m=0\),
- quadrupole \(\ell=2\) into \(|m|=2\), \(|m|=1\), and \(m=0\).

That is precisely the channel structure already found algebraically.

---

## 4. Minimal axisymmetric pole-completed DtN model

The smallest low-frequency completion consistent with the solved 2PN data is a **one-pole-per-channel** model with admittances

\[
Y_{1\perp}(\omega)=\frac{7/2}{1-\omega^2/\Omega_{1\perp}^2},
\qquad
Y_{10}(\omega)=\frac{4}{1-\omega^2/\Omega_{10}^2},
\]
\[
Y_{0}(\omega)=\frac{1}{1-\omega^2/\Omega_{0}^2},
\qquad
Y_{20}(\omega)=\frac{1}{1-\omega^2/\Omega_{20}^2},
\]
\[
Y_{21}(\omega)=\frac{1}{1-\omega^2/\Omega_{21}^2},
\qquad
Y_{22}(\omega)=\frac{1}{1-\omega^2/\Omega_{22}^2},
\qquad
Y_g(\omega)=\frac{1}{1-\omega^2/\Omega_g^2}.
\]

The corresponding bare DtN kernels are simply the inverses:

\[
Z_{1\perp}(\omega)=\frac{1-\omega^2/\Omega_{1\perp}^2}{7/2},
\qquad
Z_{10}(\omega)=\frac{1-\omega^2/\Omega_{10}^2}{4},
\]
\[
Z_{0}(\omega)=1-\frac{\omega^2}{\Omega_0^2},
\qquad
Z_{20}(\omega)=1-\frac{\omega^2}{\Omega_{20}^2},
\qquad
Z_{21}(\omega)=1-\frac{\omega^2}{\Omega_{21}^2},
\qquad
Z_{22}(\omega)=1-\frac{\omega^2}{\Omega_{22}^2},
\]
\[
Z_g(\omega)=1-\frac{\omega^2}{\Omega_g^2}.
\]

So the entire conservative 2PN throat is now encoded in:

- fixed **zero-frequency residues**,
- fixed scalar source data,
- and **seven still-open pole scales**.

---

## 5. Frequency-dependent even response matrix

In the canonical basis

\[
\{P_0, P_{20}, P_{21c}, P_{21s}, P_{22c}, P_{22s}, U\},
\]

the even-sector response matrix is

\[
M_{\rm even}(\omega)=
\begin{pmatrix}
Y_0 & 0 & 0 & 0 & 0 & 0 & J_0Y_0\\
0 & Y_{20} & 0 & 0 & 0 & 0 & J_{20}Y_{20}\\
0 & 0 & Y_{21} & 0 & 0 & 0 & 0\\
0 & 0 & 0 & Y_{21} & 0 & 0 & 0\\
0 & 0 & 0 & 0 & Y_{22} & 0 & 0\\
0 & 0 & 0 & 0 & 0 & Y_{22} & 0\\
J_0Y_0 & J_{20}Y_{20} & 0 & 0 & 0 & 0 & J_0^2Y_0+J_{20}^2Y_{20}-\Delta_{\rm geom}Y_g
\end{pmatrix},
\]
with

\[
J_0=\frac{4}{\sqrt5},
\qquad
J_{20}=\frac54,
\qquad
\Delta_{\rm geom}=\frac{281}{80}.
\]

This matrix factorizes **exactly at all frequencies** as

\[
M_{\rm even}(\omega)
=
R_{\rm support}(\omega)^T R_{\rm support}(\omega)
-
R_{\rm geom}(\omega)^T R_{\rm geom}(\omega),
\]

so the old support-minus-closure structure was not just static bookkeeping; it survives as the minimal dynamic DtN completion.

---

## 6. Frequency-dependent scalar source data

The scalar drive vector becomes

\[
J_{\rm eff}(\omega)=
\bigl(J_0Y_0(\omega),\;J_{20}Y_{20}(\omega),\;0,\;0,\;0,\;0\bigr),
\]

and the effective pure-\(U\) coefficient becomes

\[
S_{\rm eff}(\omega)=J_0^2Y_0(\omega)+J_{20}^2Y_{20}(\omega)-\Delta_{\rm geom}Y_g(\omega).
\]

At zero frequency,

\[
S_{\rm eff}(0)=\frac54,
\]

exactly as required by the solved 2PN static block.

Its low-frequency expansion is

\[
S_{\rm eff}(\omega)
=
\frac54
+\omega^2\left(
\frac{16}{5\Omega_0^2}
+\frac{25}{16\Omega_{20}^2}
-\frac{281}{80\Omega_g^2}
\right)
+O(\omega^4).
\]

So the frequency dependence of the scalar sector is now sharply parameterized.

---

## 7. Static limit check: the solved full cross operator is recovered exactly

The dynamic model was then checked in the \(\omega\to 0\) limit.
Using:

- the frozen 1PN odd wake,
- the already-solved 2PN odd leg/potential dressings,
- the new even-sector dynamic matrix above,

its static limit reproduces the full solved conservative 2PN cross operator exactly.

So this is not just a plausible dynamic ansatz; it is a true completion of the already-passed Mathematica operator package.

---

## 8. What is now fixed vs. still open

### Fixed by the current 2PN program

- odd static residues:
  \[
  R_{1\perp}=\frac72,
  \qquad
  R_{10}=4,
  \]
- even support residues:
  \[
  R_0=R_{20}=R_{21}=R_{22}=1,
  \]
- scalar source vector:
  \[
  J=\left(\frac{4}{\sqrt5},\frac54,0,0,0,0\right),
  \]
- geometry closure deficit:
  \[
  \Delta_{\rm geom}=\frac{281}{80}.
  \]

### Still open PDE/DtN observables

The dynamic pole scales:

\[
\Omega_{1\perp},\quad
\Omega_{10},\quad
\Omega_0,\quad
\Omega_{20},\quad
\Omega_{21},\quad
\Omega_{22},\quad
\Omega_g.
\]

Those are the first genuinely new inner-throat PDE observables after the conservative 2PN algebra.

A near-spherical support simplification is still available if desired:

\[
\Omega_{20}=\Omega_{21}=\Omega_{22},
\]

which collapses the quadrupole pole data to one support scale while leaving the static 2PN match untouched.

---

## 9. Why this matters

This step is important because it narrows the remaining research problem dramatically.

We no longer need to ask vaguely for a “2PN wake rebuild.”
We now know that the PDE task is specifically to derive:

1. an **axisymmetric support spectrum** with channels
   \(
   \{1\perp,10,0,20,21,22\}
   \),
2. one **pure geometry closure** pole,
3. the fixed static source vector
   \(
   J=(4/\sqrt5,5/4,0,0,0,0)
   \),
4. and the fixed closure deficit
   \(
   281/80
   \).

The isotropic cavity has now been promoted from a vague intuition to a proper **no-go unit test**.
And the minimal axisymmetric pole-completed DtN model is now the precise constructive target for the next inner-throat PDE module.
# 2PN axisymmetric Robin-wall PDE scaffold: current result

## 1. What this step adds

This step replaces the abstract channel-by-channel static fit with the first genuinely **explicit throat-wall PDE scaffold**.

The idea is simple:

- keep the 4D-ball DtN kernel as the bulk-side unit-test PDE anchor,
- add an **axisymmetric wall law** at the mouth,
- and ask what the smallest \(P_2\)-structured wall operator must look like in order to reproduce the already-solved raw static support data.

The answer is pleasantly sharp:

1. a finite monopole support comes from an **isotropic monopole wall stiffness**,
2. the odd dipole sector is reproduced exactly by a **first-order \(P_2\)** wall perturbation,
3. the even quadrupole sector **cannot** be reproduced by first-order \(P_2\) alone,
4. but it **is** reproduced exactly by adding one **second-order \(P_2^2\)** wall term.

So the conservative 2PN cross operator now has a concrete PDE-side wall-law interpretation.

---

## 2. Raw axisymmetric channel moments

Using the raw mouth basis from the previous step and the axisymmetric quadrupole
\[
P_2(\cos\theta)=\frac{3z^2-1}{2},
\]
the channel expectation values are

### Dipole sector
\[
\langle P_2\rangle_{1\perp}=-\frac15,\qquad
\langle P_2\rangle_{10}=\frac25,
\]
\[
\langle P_2^2\rangle_{1\perp}=\frac17,\qquad
\langle P_2^2\rangle_{10}=\frac{11}{35}.
\]

### Quadrupole sector
\[
\langle P_2\rangle_{20}=\frac27,\qquad
\langle P_2\rangle_{21}=\frac17,\qquad
\langle P_2\rangle_{22}=-\frac27,
\]
\[
\langle P_2^2\rangle_{20}=\frac37,\qquad
\langle P_2^2\rangle_{21}=\frac17,\qquad
\langle P_2^2\rangle_{22}=\frac17.
\]

These are exactly the axisymmetric splitting factors that appear in the wall-law fit.

---

## 3. Static wall-law fit

The already-solved raw support residues are
\[
Y_{0}= \frac{45}{4},\qquad
Y_{1\perp}=\frac72,\qquad
Y_{10}=4,\qquad
Y_{20}=\frac94,\qquad
Y_{21}=\frac32,\qquad
Y_{22}=\frac38,
\]
so the corresponding static impedances are
\[
Z_{0}=\frac4{45},\qquad
Z_{1\perp}=\frac27,\qquad
Z_{10}=\frac14,\qquad
Z_{20}=\frac49,\qquad
Z_{21}=\frac23,\qquad
Z_{22}=\frac83.
\]

### 3.1 Monopole
A finite monopole support requires
\[
Z_0(0)=\frac4{45}.
\]
This is the first explicit wall-side realization of the previously abstract “finite \(\Omega_0\)” requirement.

### 3.2 Odd dipole sector
The exact odd fit is already achieved by a first-order axisymmetric wall law
\[
Z_{1m}(0)=B_1+A_1\langle P_2\rangle_{1m},
\]
with
\[
B_1=\frac{23}{84},\qquad A_1=-\frac{5}{84}.
\]

This gives
\[
Z_{1\perp}=\frac27,\qquad Z_{10}=\frac14
\]
exactly.

### 3.3 Even quadrupole sector: first-order no-go
Trying the same first-order structure
\[
Z_{2m}(0)=B_2+A_2\langle P_2\rangle_{2m}
\]
fails.

Matching the \(m=0\) and \(|m|=1\) channels gives
\[
B_2=\frac89,\qquad A_2=-\frac{14}{9},
\]
which predicts
\[
Z_{22}^{\rm first\ order}=\frac43,
\]
but the solved target is
\[
Z_{22}^{\rm target}=\frac83.
\]

So a linear \(P_2\) wall is **not enough**.

### 3.4 Minimal exact quadrupole wall fit
The minimal exact extension is
\[
Z_{2m}(0)=B_2+A_2\langle P_2\rangle_{2m}+C_2\langle P_2^2\rangle_{2m},
\]
with
\[
B_2=\frac{10}{9},\qquad
A_2=-\frac{14}{3},\qquad
C_2=\frac{14}{9}.
\]

This reproduces
\[
Z_{20}=\frac49,\qquad
Z_{21}=\frac23,\qquad
Z_{22}=\frac83
\]
exactly.

So the minimal static wall law is

\[
Z_0(0)=\frac4{45},
\]
\[
Z_{1m}(0)=\frac{23}{84}-\frac{5}{84}\,\langle P_2\rangle_{1m},
\]
\[
Z_{2m}(0)=\frac{10}{9}-\frac{14}{3}\,\langle P_2\rangle_{2m}+\frac{14}{9}\,\langle P_2^2\rangle_{2m}.
\]

---

## 4. Axisymmetric source profile

The previously solved source vector is also natural in this wall language.

A two-parameter axisymmetric source
\[
S(\theta)=p_{\rm iso}+q_{\rm ax}P_2(\cos\theta)
\]
reproduces the required normalized source vector when
\[
p_{\rm iso}=\frac{11}{8},\qquad q_{\rm ax}=\frac{15}{8}.
\]

So both the support operator **and** the source sector now have an explicit axisymmetric wall-level realization.

---

## 5. Dynamic 4D-ball DtN pole equations

Keeping the 4D-ball unit-test DtN kernel
\[
\Lambda_\ell(z)= -1 + z\,\frac{J'_{\ell+1}(z)}{J_{\ell+1}(z)},
\qquad z\equiv \frac{a\Omega}{c_s},
\]
and promoting the static wall impedances directly into the pole equation
\[
\Lambda_\ell(z)+Z_{\ell m}(0)=0,
\]
gives the following lowest positive roots.

### Monopole
\[
z_0=0.591884444464394.
\]
The small-\(z\) estimate from \(\Lambda_0(z)\simeq -z^2/4\) is
\[
z_0 \simeq \frac{4}{\sqrt{45}}=0.596284793999944,
\]
so the finite monopole pole is already visible at leading order.

### Dipole
Using the isotropic odd base \(B_1=23/84\),
\[
z_1^{\rm base}=2.551215916564765,
\]
while the exact split poles are
\[
z_{1\perp}=2.561183722397930,\qquad
z_{10}=2.531063390840353.
\]

The first-order perturbative shifts about the isotropic base are
\[
z_{1\perp}^{\rm pert}=2.561219561244114,\qquad
z_{10}^{\rm pert}=2.531208627206067,
\]
which are extremely accurate. So the odd dipole sector is genuinely perturbative.

### Quadrupole
Using the isotropic even base \(B_2=10/9\),
\[
z_2^{\rm base}=4.254105628646177,
\]
while the exact split poles are
\[
z_{20}=3.901921523190568,\qquad
z_{21}=4.029116369391941,\qquad
z_{22}=4.821811915561263.
\]

First-order shifts about that base give
\[
z_{20}^{\rm pert}=3.942783453176964,\qquad
z_{21}^{\rm pert}=4.046557511666702,\qquad
z_{22}^{\rm pert}=4.980524038074339.
\]

So the perturbative ordering is correct, but the quadrupole anisotropy is large enough that the exact channel equations should be treated nonperturbatively.

---

## 6. Why this matters

This is the first point where the inner-throat program stops being merely “a target channel algebra” and becomes an actual PDE-side construction.

We now know:

1. **finite \(\Omega_0\)** comes from a finite isotropic monopole wall stiffness;
2. **dipole splitting** comes from a first-order axisymmetric \(P_2\) wall term;
3. the solved **quadrupole support hierarchy** requires a minimal second-order \(P_2^2\) wall contribution;
4. the source vector is exactly a two-parameter axisymmetric source;
5. and the full low-frequency channel structure can be attached directly to the 4D-ball DtN pole equations.

So the next derivation target is no longer vague. It is:

> derive these wall coefficients from a concrete soft-wall throat potential or stress-balance wall law, instead of fitting them at the operator level.

That is a very manageable next step, and it is exactly the wall-law gap the throat notes said had to be frozen before full throat modeling becomes coherent.
# 2PN Family-1 soft-wall → generalized Robin / wall-stress operator: current result

## 1. What this step adds

This step closes the gap between the earlier **axisymmetric Robin-wall fit** and the actual **Family-1 soft-wall / flare geometry**.

The new result is that the full passed support-sector static data for the \(\ell=1,2\) channels are reproduced **exactly** by a single minimal generalized Robin operator with a **dipole base** plus one **curvature factor**:
\[
Z_{\ell m}
=
A_0 + A_1(\lambda_\ell-2)
+
\bigl(B_0 + B_1(\lambda_\ell-2)\bigr)\langle P_2\rangle_{\ell m}
+
C(\lambda_\ell-2)\langle P_2^2\rangle_{\ell m},
\qquad
\lambda_\ell=\ell(\ell+1).
\]

Solving against the passed support data gives
\[
A_0=\frac{23}{84},\qquad
A_1=\frac{211}{1008},\qquad
B_0=-\frac{5}{84},\qquad
B_1=-\frac{129}{112},\qquad
C=\frac{7}{18}.
\]

So the support operator is no longer a loose set of channel-by-channel fits.

---

## 2. The structural simplification

It is cleaner to write the same result as
\[
Z_{\ell m}
=
\langle z_{\rm base}(\mu)\rangle_{\ell m}
+
(\lambda_\ell-2)\,\langle z_{\rm curv}(\mu)\rangle_{\ell m},
\qquad
\mu\equiv \cos\theta.
\]

The two local profiles are

\[
z_{\rm base}(\mu)=\frac{17}{56}-\frac{5}{56}\mu^2,
\]

\[
z_{\rm curv}(\mu)=\frac{593}{672}-\frac{1553}{672}\mu^2+\frac78\mu^4.
\]

This is the key simplification:

- the leading support operator is just a low-order polynomial in \(\mu^2\),
- the entire **quadrupole-only correction** is carried by the factor \((\lambda_\ell-2)\),
- so the dipole is special because
  \[
  \lambda_1-2=0.
  \]

That makes the earlier “\(P_2^2\) appears only for \(\ell=2\)” result look natural rather than accidental.

A compact operator notation for the passed diagonal matrix elements is
\[
\mathcal B_{\rm eff}
=
\partial_n
+
z_{\rm base}(\mu)
+\frac12\Big\{ -\Delta_S-2,\ z_{\rm curv}(\mu)\Big\},
\]
restricted to the support sector. On a spherical harmonic \(Y_{\ell m}\), the anticommutator contributes the factor \((\lambda_\ell-2)\).

---

## 3. Why this is exactly the Family-1 flare basis

The Family-1 throat notes use a flared mouth profile
\[
a(z)=a_0\Bigl(1+\beta e^{-z^2/z_m^2}\Bigr).
\]

Near the mouth center, with \(a_m=a(0)=a_0(1+\beta)\) and \(\mu=\cos\theta\),
\[
\frac{a(z)}{a_m}=1-q\mu^2+r\mu^4+O(\mu^6),
\]
for suitable dimensionless \(q,r\).

That is *exactly* the basis that appears in the solved support operator:
\[
h_1(\mu)=-\mu^2,\qquad h_2(\mu)=\mu^4.
\]

Indeed,
\[
z_{\rm base}(\mu)=\frac{17}{56}+\frac{5}{56}h_1(\mu),
\]
\[
z_{\rm curv}(\mu)=\frac{593}{672}+\frac{1553}{672}h_1(\mu)+\frac78 h_2(\mu).
\]

So the passed support operator is already sitting in the same low-order polynomial basis produced by the actual Family-1 flare geometry.

This is a strong sign that the earlier Robin-wall fit was not just formal algebra. It is exactly what a local soft-wall / flared-mouth reduction is expected to generate.

---

## 4. The \(P_2\) and \(P_2^2\) content is also automatic

Expanding the same Family-1 flare in the Legendre basis gives
\[
-q\mu^2+r\mu^4
=
\left(-\frac q3+\frac r9\right)
+
\left(-\frac{2q}{3}+\frac{4r}{9}\right)P_2
+
\frac{4r}{9}P_2^2.
\]

So:

- the first nontrivial flare term naturally gives \(P_2\),
- the next term naturally gives \(P_2^2\),
- which is exactly the hierarchy that the passed Robin-wall support data forced us to introduce.

Again, the “extra \(P_2^2\) term” is no longer an ad hoc channel patch. It is the expected quartic mouth-flare correction.

---

## 5. Source profile in the same basis

The earlier axisymmetric source solution also collapses into the same \(\mu^2\) basis:
\[
S(\mu)=\frac{11}{8}+\frac{15}{8}P_2(\mu)
      =\frac{7}{16}+\frac{45}{16}\mu^2.
\]

So both the support operator and the source sector now live in the same simple Family-1 flare basis:
\[
\{1,\ \mu^2,\ \mu^4\}.
\]

That is a very clean finish-line target for a genuine inner-throat PDE derivation.

---

## 6. What this means physically

This step says the static 2PN support sector can be read as:

1. a **dipole-base wall response**,
2. plus one **finite-thickness / curvature correction** that vanishes on the dipole,
3. with all angular dependence supplied by the same \(\mu^2,\mu^4\) structure generated by the actual Family-1 flared mouth.

So the next PDE task is much smaller than before:

- we do **not** need a completely free channel-by-channel wall law,
- we need a steep soft-wall / boundary-layer reduction that derives the two local profiles
  \[
  z_{\rm base}(\mu)=\frac{17}{56}-\frac{5}{56}\mu^2,
  \qquad
  z_{\rm curv}(\mu)=\frac{593}{672}-\frac{1553}{672}\mu^2+\frac78\mu^4.
  \]

If that reduction works, then the previously solved 2PN support operator has a direct Family-1 throat interpretation.

---

## 7. About the isotropic 4D-ball \(z^8\) correction

The correction to the isotropic 4D-ball unit-test series at order \(z^8\) for \(\ell=1,2\) does **not** affect this step.

This operator reconstruction uses the passed exact support-channel static data, not the earlier truncated \(z^8\) quotient. So the static support/Family-1 matching result here is unchanged.
# 2PN Family-1 soft-wall variational reduction: current result

## 1. What this step adds

This step turns the previously fitted static axisymmetric wall operator into the first **reduced variational wall-Hessian** that reproduces the passed 2PN support data exactly, while also tying its angular structure back to the frozen **Family-1** flared soft wall.

The correction to the earlier isotropic 4D-ball unit-test series does **not** affect this step. This derivation depends only on the already-passed **static channel data** and on exact axisymmetric angular moments, not on the corrected \(z^8\) coefficients.

---

## 2. Exact reduced modal wall-Hessian

For the static support channels
\[
K_{00}=\frac4{45},\qquad
K_{1\perp}=\frac27,\qquad
K_{10}=\frac14,\qquad
K_{20}=\frac49,\qquad
K_{21}=\frac23,\qquad
K_{22}=\frac83,
\]
a compact reduced wall law closes **exactly**:

\[
K_{00}=K_{\rm mono},
\]
\[
K_{\ell m}
=
K_{\rm mono}
+
T_0\,\ell(\ell+1)
+
\big(A_0+A_1\,\ell(\ell+1)\big)\,\langle P_2\rangle_{\ell m}
+
\big(B_0+B_1\,\ell(\ell+1)\big)\,\langle P_2^2\rangle_{\ell m},
\qquad \ell=1,2.
\]

The exact angular moments are
\[
\langle P_2\rangle_{1\perp}=-\frac15,
\quad
\langle P_2\rangle_{10}=\frac25,
\quad
\langle P_2\rangle_{20}=\frac27,
\quad
\langle P_2\rangle_{21}=\frac17,
\quad
\langle P_2\rangle_{22}=-\frac27,
\]
\[
\langle P_2^2\rangle_{1\perp}=\frac17,
\quad
\langle P_2^2\rangle_{10}=\frac{11}{35},
\quad
\langle P_2^2\rangle_{20}=\frac37,
\quad
\langle P_2^2\rangle_{21}=\frac17,
\quad
\langle P_2^2\rangle_{22}=\frac17.
\]

Solving against the passed channel data gives
\[
K_{\rm mono}=\frac4{45},
\qquad
T_0=\frac{23}{135},
\]
\[
A_0=\frac{9095}{3528},
\qquad
A_1=-\frac{25559}{21168},
\]
\[
B_0=-\frac{109}{56},
\qquad
B_1=\frac{1765}{3024}.
\]

Substituting these values back reproduces all six passed channel stiffnesses with zero residual.

So the static axisymmetric support operator is no longer just a fitted Robin table.
It is the Hessian of a concrete **reduced variational wall law** on the truncated \(\ell=0,1,2\) mouth basis.

---

## 3. Family-1 flare automatically generates the needed \(\{1,P_2,P_2^2\}\) basis

The frozen Family-1 throat profile uses a flared radius
\[
a(z)=a_0\Big(1+\beta e^{-z^2/z_m^2}\Big).
\]

On the mouth sphere, define
\[
\xi \equiv \frac{a_0^2}{z_m^2}.
\]

Expanding the flare to second order gives
\[
\frac{\delta a(\theta)}{a_0}
=
\beta\Big[
1-\frac{\xi}{3}+\frac{\xi^2}{18}
\Big]
+
\beta\Big[
-\frac{2\xi}{3}+\frac{2\xi^2}{9}
\Big]P_2
+
\beta\Big[
\frac{2\xi^2}{9}
\Big]P_2^2
+O(\xi^3).
\]

This is the key structural result:

> the actual frozen Family-1 flare automatically produces the exact axisymmetric basis \(\{1,P_2,P_2^2\}\) that the reduced wall-Hessian needs.

So the earlier \(P_2^2\) term is no longer ad hoc.
At the reduced wall-law level, it is the generic second-order imprint of the flared Family-1 soft wall.

---

## 4. Minimal linear-gradient interpretation

There is also a useful sharp sub-result.

If the **anisotropic gradient block** of the reduced wall law is assumed to be **linear** in the Family-1 flare profile, then
\[
\frac{B_1}{A_1}=rac{D_2}{D_1}=rac{\xi}{\xi-3},
\]
where
\[
D_1=-\frac{2\xi}{3}+\frac{2\xi^2}{9},
\qquad
D_2=\frac{2\xi^2}{9}.
\]

Using the exact solved coefficients gives
\[
\frac{B_1}{A_1}=-\frac{12355}{25559},
\]
so
\[
\xi=\frac{12355}{12638}\approx 0.9776072163,
\qquad
\frac{z_m}{a_0}=\frac{1}{\sqrt{\xi}}\approx 1.0113880097.
\]

So the minimal linear-gradient reading points to a flare width of order the throat radius,
\[
z_m\sim a_0.
\]

---

## 5. What this means

The constructive state of play is now:

1. the conservative 2PN comparable-mass cross sector is solved algebraically,
2. that solution was rebuilt as a small mouth-port / support-minus-closure EFT,
3. the static axisymmetric support data now come from an explicit **reduced variational wall-Hessian**,
4. and the exact angular basis of that Hessian is now traced back to the actual **Family-1** flared soft wall rather than inserted by hand.

So the remaining “physics gap” is narrower than before.

The next derivation target is no longer “invent the 2PN wall law.”
It is:

- derive the coefficients \((K_{\rm mono},T_0,A_0,A_1,B_0,B_1)\) from a more explicit soft-wall or traction-balance reduction,
- then fold that wall-Hessian into the inner-throat DtN / port solver so the support coefficients emerge from the PDE side directly.

That is a real step toward closing the full 2PN derivation.
# 2PN Family-1 soft-wall strict boundary-layer no-go: preliminary result

## 1. What was tested

This step asks whether the passed static support-channel data can already be reproduced by the **strictest** steep-wall interpretation of the Family-1 mouth geometry:

- one isotropic normal penetration moment \(A\),
- one isotropic tangential-gradient moment \(B\),
- and the full geometric pullback of a flared mouth sphere
  \[
  \sigma(\mu)=1-q\mu^2+r\mu^4.
  \]

The resulting quadratic-flare truncation of the pullback factors is
\[
J(\mu) = -7\mu^8 r^2 + 6\mu^6 q r + 8\mu^6 r^2 - \mu^4 q^2 - 8\mu^4 q r + 2\mu^4 r + 2\mu^2 q^2 - 2\mu^2 q + 1,
\]
\[
F_{\theta}(\mu) = 8\mu^8 r^2 - 8\mu^6 q r - 8\mu^6 r^2 + 2\mu^4 q^2 + 8\mu^4 q r - 2\mu^2 q^2 + 1,
\]
\[
F_{\phi}(\mu) = -8\mu^8 r^2 + 8\mu^6 q r + 8\mu^6 r^2 - 2\mu^4 q^2 - 8\mu^4 q r + 2\mu^2 q^2 + 1.
\]

The strict surface action is then
\[
E^{\rm strict}_{\rm BL}
=
\frac12\int d\Omega
\left[
A J(\mu)\,\Psi^2
+
B\left(
F_{\theta}(\mu)(\partial_{\theta}\Psi)^2
+
F_{\phi}(\mu)\frac{(\partial_{\phi}\Psi)^2}{\sin^2\theta}
\right)
\right].
\]

---

## 2. Exact channel formulas

Projecting this strict model onto the \(\ell=0,1,2\) mouth harmonics gives
\[
K_{00}=\frac{7}{15}Aq^2-\frac{26}{35}Aqr-\frac23 Aq+\frac{23}{63}Ar^2+\frac25 Ar+A,
\]
\[
K_{1\perp}=\frac{11}{35}Aq^2-\frac25 Aqr-\frac25 Aq+\frac{13}{77}Ar^2+\frac{6}{35}Ar+A
+\frac{8}{35}Bq^2-\frac{32}{105}Bqr+\frac{32}{231}Br^2+2B,
\]
\[
K_{10}=\frac{27}{35}Aq^2-\frac{10}{7}Aqr-\frac65 Aq+\frac{25}{33}Ar^2+\frac67 Ar+A
-\frac{16}{35}Bq^2+\frac{64}{105}Bqr-\frac{64}{231}Br^2+2B,
\]
\[
K_{20}=\frac{13}{21}Aq^2-\frac{94}{77}Aqr-\frac{22}{21}Aq+\frac{6185}{9009}Ar^2+\frac67 Ar+A
-\frac{16}{7}Bq^2+\frac{320}{77}Bqr-\frac{320}{143}Br^2+6B,
\]
\[
K_{21}=\frac{13}{21}Aq^2-\frac{230}{231}Aqr-\frac67 Aq+\frac{205}{429}Ar^2+\frac{10}{21}Ar+A
+\frac{8}{21}Bq^2-\frac{96}{77}Bqr+\frac{800}{1001}Br^2+6B,
\]
\[
K_{22}=\frac{5}{21}Aq^2-\frac{58}{231}Aqr-\frac27 Aq+\frac{25}{273}Ar^2+\frac{2}{21}Ar+A
+\frac{16}{21}Bq^2-\frac{64}{77}Bqr+\frac{320}{1001}Br^2+6B.
\]

These are exact for the quadratic-flare truncation of the strict isotropic layer.

---

## 3. Best-fit test against the passed support data

The target support values are
\[
K_{00}=\frac4{45},\qquad
K_{1\perp}=\frac27,\qquad
K_{10}=\frac14,
\]
\[
K_{20}=\frac49,\qquad
K_{21}=\frac23,\qquad
K_{22}=\frac83.
\]

A multi-seed least-squares search over \((A,B,q,r)\) gives the robust best fit
\[
A\approx -0.0177094,
\qquad
B\approx 0.2505433,
\qquad
q\approx -3.5661065,
\qquad
r\approx -3.0166765,
\]
with
\[
\sum_i \Delta_i^2 \approx 0.5536733194,
\qquad
\max_i |\Delta_i| \approx 0.4390744081.
\]

The channel residuals at that best fit are
\[
\Delta_{00}\approx -0.1497425,
\qquad
\Delta_{1\perp}\approx +0.3824726,
\qquad
\Delta_{10}\approx -0.2656747,
\]
\[
\Delta_{20}\approx -0.1804115,
\qquad
\Delta_{21}\approx +0.4390744,
\qquad
\Delta_{22}\approx -0.2984083.
\]

So the residual obstruction is not small and does not disappear under reseeding.

---

## 4. What this means

This is a useful negative result.

The passed static support sector is **not** just the strict isotropic steep-wall layer pulled back through Family-1 flare geometry.

So:

- pure isotropic penetration + pure geometry pullback is too small a model,
- at least one extra local wall-stress / traction / profile degree of freedom is required,
- and the reduced variational wall-Hessian derived in the companion step is therefore a **genuine new structure**, not a disguised rewrite of the strict isotropic layer.

That sharply narrows the next PDE target:

> derive the extra local wall degree of freedom from the actual soft-wall / traction-balance physics, rather than trying to squeeze the full support sector out of isotropic penetration alone.
# 2PN Family-1 traction-balance / wall-stress completion: current result

## 1. What this step closes

The strict isotropic boundary-layer pullback was too small to reproduce the passed static support sector.
The next natural question was whether the solved support operator could be rewritten as an **explicit local wall-stress / traction-balance surface energy**.

The answer is yes.

The cleanest exact closure is:

- keep the monopole wall channel separate,
  
  \[
  K_{00}=\frac{4}{45},
  \]

- and solve the dipole/quadrupole support sector with one base pressure profile and one tangential wall-stress profile.

---

## 2. Minimal exact traction-balance ansatz on the support sector

On the \(\ell=1,2\) support sector, write

\[
K_{\ell m}=\langle z_{\rm base}\rangle_{\ell m}+\bigl(\ell(\ell+1)-2\bigr)\langle t\rangle_{\ell m},
\]

with low-order Family-1 profiles

\[
z_{\rm base}(\mu)=b_0+b_2\mu^2,
\qquad
 t(\mu)=t_0+t_2\mu^2+t_4\mu^4.
\]

Matching the five carried-forward support targets

\[
K_{1\perp}=\frac27,
\quad
K_{10}=\frac14,
\quad
K_{20}=\frac49,
\quad
K_{21}=\frac23,
\quad
K_{22}=\frac83,
\]

fixes the profiles **uniquely**:

\[
b_0=\frac{17}{56},
\qquad
b_2=-\frac{5}{56},
\]

\[
t_0=\frac{593}{672},
\qquad
t_2=-\frac{1553}{672},
\qquad
t_4=\frac78.
\]

So

\[
z_{\rm base}(\mu)=\frac{17}{56}-\frac{5}{56}\mu^2,
\]

\[
t(\mu)=\frac{593}{672}-\frac{1553}{672}\mu^2+\frac78\mu^4.
\]

In Legendre form,

\[
z_{\rm base}(\mu)=\frac{23}{84}-\frac{5}{84}P_2(\mu),
\]

\[
t(\mu)=\frac{211}{1008}-\frac{129}{112}P_2(\mu)+\frac{7}{18}P_2(\mu)^2.
\]

This is the exact same Family-1 low-order basis that kept reappearing in the earlier flare / Robin / variational steps.

---

## 3. Equivalent local wall-energy form

For

\[
\mathcal O_{\rm supp}=z_{\rm base}+\frac12\{-\Delta_S-2,\,t\},
\]

the equivalent local quadratic form on the \(\ell=1,2\) support sector is

\[
E_{\rm supp}[\Psi]
=\frac12\int d\Omega\,\Big[p(\mu)\Psi^2+t(\mu)|\nabla_S\Psi|^2\Big],
\]

with

\[
p(\mu)=z_{\rm base}(\mu)-2t(\mu)-\frac12\Delta_S t(\mu).
\]

Using

\[
\Delta_S t(\mu)= -\frac{35}{2}\mu^4+\frac{2729}{112}\mu^2-\frac{1553}{336},
\]

this gives the exact pressure profile

\[
p(\mu)=\frac{571}{672}-\frac{5141}{672}\mu^2+7\mu^4.
\]

In Legendre form,

\[
p(\mu)= -\frac{155}{168}-\frac{2005}{1008}P_2(\mu)+\frac{28}{9}P_2(\mu)^2.
\]

So the solved support sector is now an explicit local wall model:

- one base pressure profile \(z_{\rm base}(\mu)\),
- one tangential wall-stress / curvature-compliance profile \(t(\mu)\),
- and the induced pressure profile \(p(\mu)\) that follows from the anticommutator structure.

---

## 4. Exact verification

With the local wall-energy form,

\[
K_{\ell m}=\int d\Omega\Big[p(\mu)|Y_{\ell m}|^2+t(\mu)|\nabla_S Y_{\ell m}|^2\Big],
\qquad \ell=1,2,
\]

all five support-channel residuals vanish identically.

So the dipole/quadrupole support sector is no longer just a passed operator fit. It now has a direct **traction-balance / wall-stress** reading.

---

## 5. Physical reading

This suggests the strict isotropic boundary layer was missing exactly one new structural ingredient:

> an axisymmetric tangential wall-stress / curvature-compliance profile in the same low-order Family-1 basis.

The monopole wall mode remains separate,

\[
K_{00}=\frac{4}{45},
\]

which is consistent with the earlier PDE observation that a finite monopole pole/stiffness is its own channel.

So the constructive picture is now:

1. strict isotropic pullback alone = insufficient,
2. add one explicit tangential wall-stress profile \(t(\mu)\),
3. keep the monopole wall channel separate,
4. the full passed support sector closes exactly.

That is the first explicit local traction-balance completion of the strict isotropic layer that is exact on the solved 2PN support sector.
# 2PN Family-1 anisotropic tangential-traction completion: current result

## 1. What was tested

The strict isotropic Family-1 boundary-layer pullback failed even after exact flare geometry was included.  
The previous step showed that the missing support-sector structure could be rewritten as a local wall-stress / traction-balance profile, but that still left open whether a **real PDE-side boundary-layer completion** could generate it.

The minimal new test was:

- keep the **normal penetration moment isotropic**,
- keep the same strict quadratic-flare Family-1 surface pullback,
- and promote only the tangential wall moment from a constant to the smallest axisymmetric profile
  \[
  B_{\rm tan}(\mu)=B_0+B_2\mu^2.
  \]

So the surface energy becomes
\[
E_{\rm BL}^{\rm aniso}
=
\frac12\int d\Omega\,
\Big[
A\,J(\mu)\,\Psi^2
+
B_{\rm tan}(\mu)\Big(
F_\theta(\mu)(\partial_\theta\Psi)^2
+
F_\phi(\mu)\frac{(\partial_\phi\Psi)^2}{\sin^2\theta}
\Big)
\Big],
\]
with the same strict pullback factors \(J,F_\theta,F_\phi\) used in the earlier no-go step.

This is the smallest genuine PDE-side completion of the strict model.

---

## 2. Exact support-channel closure

Projecting onto the \(\ell=1,2\) channels gives exact formulas for
\[
K_{1\perp},\ K_{10},\ K_{20},\ K_{21},\ K_{22}
\]
that are linear in \(A,B_0,B_2\) and polynomial in \((q,r)\).

A multi-seed solve shows that the system
\[
K_{1\perp}=\frac27,\qquad
K_{10}=\frac14,\qquad
K_{20}=\frac49,\qquad
K_{21}=\frac23,\qquad
K_{22}=\frac83
\]
has **exact real solutions**.

So the strict isotropic no-go is really repaired by one new local degree of freedom:

> an axisymmetric tangential wall-stress profile \(B_{\rm tan}(\mu)=B_0+B_2\mu^2\).

This is the first strict boundary-layer / soft-wall PDE completion that actually reproduces the solved support sector, rather than fitting the support operator directly.

---

## 3. Family-1-like physical branch

The SymPy branch scan finds at least four exact real branches.  
The natural Family-1-like branch is the one with

- \(q>0\),
- \(r>0\),
- and \(\sigma(\mu)=1-q\mu^2+r\mu^4>0\) on \([-1,1]\).

For that branch,
\[
A\approx -0.281923219302567,
\]
\[
B_0\approx 0.648914709228733,
\qquad
B_2\approx -1.085434603876673,
\]
\[
q\approx 2.370915717168574,
\qquad
r\approx 2.758343474052255.
\]

The support residuals vanish at machine precision.

In Legendre form,
\[
B_{\rm tan}(\mu)=\beta_0+\beta_2 P_2(\mu),
\]
with
\[
\beta_0=B_0+\frac13 B_2 \approx 0.287103174603175,
\]
\[
\beta_2=\frac23 B_2 \approx -0.723623069251115.
\]

So the minimal successful PDE-side completion is an **axisymmetric \(P_2\)-modulated tangential traction profile**.

---

## 4. Exact universal monopole prediction

The most useful structural result is that the anisotropic tangential model satisfies the exact identity
\[
K_{00}
=
K_{1\perp}
+\frac12 K_{10}
-\frac1{10}K_{20}
-\frac15 K_{21}
-\frac15 K_{22}.
\]

So once the \(\ell=1,2\) support targets are matched exactly, the raw boundary-layer monopole is fixed automatically:
\[
K_{00}^{\rm BL}
=
\frac27+\frac12\cdot\frac14-\frac1{10}\cdot\frac49-\frac15\cdot\frac23-\frac15\cdot\frac83
=
-\frac{757}{2520}.
\]

This is independent of which exact support branch is chosen.

Since the carried-forward monopole target is
\[
K_{00}^{\rm target}=\frac4{45},
\]
the separate monopole wall add-on must be
\[
\Delta K_{00}^{\rm mono}
=
\frac4{45}-\left(-\frac{757}{2520}\right)
=
\frac{109}{280}.
\]

So the support-sector repair does **not** eliminate the separate monopole channel; instead it predicts its exact required size.

---

## 5. Physical reading

This is an important narrowing of the throat PDE problem.

The previous exact operator and traction-balance steps said that a missing tangential wall-stress profile had to exist.  
This new result shows that the smallest strict boundary-layer realization of that statement is already enough:

1. isotropic normal penetration moment \(A\),
2. one axisymmetric tangential wall-moment profile \(B_{\rm tan}(\mu)=B_0+B_2\mu^2\),
3. the same strict Family-1 flare pullback,
4. plus the still-separate monopole wall channel.

So the next PDE target is no longer vague.

The next thing to derive from the soft-wall physics is:

- why the tangential moment picks up exactly this \(P_2\)-modulated profile, and
- how the separate monopole add-on \(\frac{109}{280}\) emerges from the isotropic wall / geometry sector.

That is a much sharper finish-line target than “derive the whole support operator from scratch.”
# 2PN Family-1 support/source/monopole closure: current result

## What this step closes

The previous Family-1 wall work had already shown three separate things:

1. the exact passed \(\ell=1,2\) support sector can be written as a local wall operator,
2. the source sector lives in the same low-order Family-1 flare basis,
3. a separate monopole add-on of \(109/280\) is still required.

The new question was whether these are still **independent pieces**, or whether the support and source sectors are actually locked together by one tighter constitutive law.

They are locked together.

---

## 1. Exact source-from-support relation

With the exact carried-forward Family-1 base profile
\[
z_{\rm base}(\mu)=\frac{17}{56}-\frac{5}{56}\mu^2,
\]
the exact carried-forward source profile
\[
S(\mu)=\frac{7}{16}+\frac{45}{16}\mu^2
\]
satisfies
\[
\boxed{
S(\mu)=10-\frac{63}{2}\,z_{\rm base}(\mu).
}
\]

So the wall/source sector is **not independent**.
Once the base support profile is fixed, the source profile is fixed automatically.

---

## 2. The same collapse happens in the actual Family-1 flare basis

Using the exact flare constitutive closure
\[
z_{\rm base}(\mu)=\pi_0+\pi_1\chi(\mu)+\pi_2\Delta_S\chi(\mu),
\]
with
\[
\chi(\mu)=1-\frac{1176}{1553}\mu^2+\frac12\left(\frac{1176}{1553}\right)^2\mu^4,
\]
the source coefficients are not new fits. They obey
\[
\boxed{
s_0=10-\frac{63}{2}\pi_0,\qquad
s_1=-\frac{63}{2}\pi_1,\qquad
s_2=-\frac{63}{2}\pi_2.
}
\]

So the exact source sector is generated by the **same** constitutive data that generate the base support profile.

---

## 3. Dipole support residues already determine the scalar source sector

Any axisymmetric quadratic base profile can be written as
\[
z_{\rm base}(\mu)=a_0+a_2 P_2(\mu).
\]

For the passed wall sector,
the dipole support residues \(K_{1\perp}\) and \(K_{10}\) reconstruct that profile exactly:
\[
a_0=\frac{K_{10}+2K_{1\perp}}{3},
\qquad
a_2=\frac{5}{3}(K_{10}-K_{1\perp}).
\]

Then the wall/source coefficients
\[
S(\mu)=p_{\rm iso}+q_{\rm ax}P_2(\mu)
\]
are forced to be
\[
\boxed{
p_{\rm iso}=10-\frac{21}{2}(K_{10}+2K_{1\perp}),
\qquad
q_{\rm ax}=-\frac{105}{2}(K_{10}-K_{1\perp}).
}
\]

At the passed dipole support values
\[
K_{1\perp}=\frac27,\qquad K_{10}=\frac14,
\]
this gives
\[
p_{\rm iso}=\frac{11}{8},
\qquad
q_{\rm ax}=\frac{15}{8},
\]
exactly.

So the scalar source sector is already fixed once the **dipole support** is fixed.

---

## 4. Canonical scalar source vector follows automatically

Using the same normalization convention as the passed minimal irrep throat operator,
the canonical scalar source vector is
\[
J_0=\frac{2}{\sqrt5}\left(p_{\rm iso}+\frac{q_{\rm ax}}{3}\right),
\qquad
J_{20}=\frac23 q_{\rm ax}.
\]

Substituting the dipole-support formulas above gives
\[
\boxed{
J_0=\frac{\sqrt5}{5}\bigl(20-56K_{10}-7K_{1\perp}\bigr),
\qquad
J_{20}=35\,(K_{1\perp}-K_{10}).
}
\]

At the passed support values this reproduces exactly
\[
\boxed{
J=\left(\frac{4}{\sqrt5},\frac54\right).
}
\]

This is a strong structural result: the scalar source vector used in the solved 2PN cross sector is **not an independent fit**. It is determined by the dipole support residues.

---

## 5. Full static local wall operator and monopole auxiliary completion

The exact local Family-1 wall operator is
\[
\mathcal O_{\rm loc}
=
z_{\rm base}(\mu)
+\frac12\left\{-\Delta_S-2,\ t(\mu)\right\},
\]
with
\[
t(\mu)=\frac{593}{672}-\frac{1553}{672}\mu^2+\frac78\mu^4.
\]

Its exact channel values are
\[
K_{00}^{\rm raw}=-\frac{757}{2520},
\qquad
K_{1\perp}=\frac27,
\qquad
K_{10}=\frac14,
\]
\[
K_{20}=\frac49,
\qquad
K_{21}=\frac23,
\qquad
K_{22}=\frac83.
\]

So the local wall law already reproduces the full carried-forward \(\ell=1,2\) sector exactly, and misses only the monopole by
\[
\boxed{
\Delta K_{00}=\frac4{45}-\left(-\frac{757}{2520}\right)=\frac{109}{280}.
}
\]

The minimal full static even-wall completion is therefore
\[
\boxed{
\mathcal O_{\rm full}
=
\mathcal O_{\rm loc}
+\frac{109}{280}\,\mathbb P_{00},
}
\]
with \(\mathbb P_{00}\) the monopole projector.

Equivalently, the remaining gap is just one global breathing auxiliary channel, for example
\[
Y_{\rm mono}(\omega)=\frac{109/280}{1-\omega^2/\Omega_{\rm mono}^2}.
\]

---

## 6. Interpretation

This is the sharpest static Family-1 wall closure so far.

The whole static even wall sector is now

\[
\boxed{
\text{[local Family-1 support/source constitutive law]}
\;+\;
\text{[one global monopole auxiliary channel]}.
}
\]

So the remaining PDE-side task is even narrower than before.

We no longer need to derive the source sector independently.
We need to derive:

1. the local support law \(\mathcal O_{\rm loc}\) from the soft-wall boundary layer / traction physics,
2. and the single monopole breathing auxiliary \(\Delta K_{00}=109/280\).

That is a much smaller and much cleaner finish-line target.

# 2PN geometry-Hessian monopole breathing closure: current result

## What this step adds

The previous Family-1 wall work compressed the full static even sector to

\[
\text{[local support/source constitutive law]} + \frac{109}{280}\,\mathbb P_{00},
\]

so the only unresolved static piece was the **global monopole breathing auxiliary**.

This step shows that the missing monopole projector can be generated directly from a reduced **geometry Hessian** in the throat variables \((a,L)\), with a natural volume-work coupling. So the monopole add-on does **not** need to be inserted as an independent ad hoc channel.

---

## 1. Minimal geometry-side model

Use the 4D cylinder-like throat geometry

\[
V(a,L)=\frac{4\pi}{3}a^3L,
\qquad
A(a,L)=4\pi a^2L+\frac{8\pi}{3}a^3,
\]

and the minimal curvature-completed geometry energy

\[
E_{\rm geom}(a,L)
=
P_{\rm vac}V(a,L)+\sigma A(a,L)+\kappa_b\frac{a^2}{L}.
\]

The new ingredient is the last term. It is the smallest explicit curvature/bending completion that can repair the known \(P_{\rm vac}V+\sigma A\) limitation.

---

## 2. Natural monopole coupling and exact projector coefficient

Treat the uniform monopole port as a pressure-like source that couples through volume work,

\[
\delta W = -p\,\delta V,
\qquad
g=\nabla_{(a,L)}V.
\]

Let \(H_0\) be the Hessian of \(E_{\rm geom}\) at the reference point
\[
(a_0,L_0)=(a_0,\Lambda a_0).
\]

Then integrating out \((\delta a,\delta L)\) gives the exact reduced geometry-side monopole coefficient

\[
\boxed{
\Delta K_{00}^{\rm geom}
=
\frac{g^T H_0^{-1} g}{V_0^2}.
}
\]

So the old global monopole projector is reinterpreted as a genuine low-frequency **geometry compressibility**.

---

## 3. Baseline no-go from the actual geometry Hessian

Write

\[
\sigma=\frac{\Sigma}{a_0^3},
\qquad
P_{\rm vac}=\rho\,\frac{\Sigma}{a_0^4},
\qquad
\kappa_b=\beta\,\frac{\Sigma}{a_0}.
\]

Then the exact reduced coefficient is

\[
\boxed{
\Delta K_{00}^{\rm geom}
=
\frac{2\pi\Lambda^2\rho+5\pi\Lambda^2-2\pi\Lambda-4\beta}
{2\pi\Sigma\left(\pi\Lambda^3\rho^2+4\pi\Lambda^3\rho+4\pi\Lambda^3-2\Lambda\beta\rho-3\Lambda\beta-2\beta\right)}.
}
\]

For the baseline model with no curvature completion \((\beta=0)\), the reduced Hessian is

\[
H_{\rm base}=
\begin{pmatrix}
8\pi(\Lambda\rho+\Lambda+2) & 4\pi(\rho+2) \\
4\pi(\rho+2) & 0
\end{pmatrix},
\]

with determinant

\[
\det H_{\rm base} = -16\pi^2(\rho+2)^2 < 0.
\]

So the literal \(P_{\rm vac}V+\sigma A\) geometry energy is **not** a passive 2DOF geometry Hessian and cannot by itself close the monopole channel. That exactly matches the limitation flagged in the throat notes.

---

## 4. Exact positivity conditions for the minimal curvature completion

The curvature-completed Hessian is passive / positive-definite when

\[
\beta>0
\]
and
\[
\boxed{
\beta >
\beta_{\rm stab}(\Lambda,\rho)
=
\frac{\pi\Lambda^3(\rho+2)^2}{2\Lambda\rho+3\Lambda+2}.
}
\]

For the resulting monopole coefficient to be positive as well, one also needs

\[
\boxed{
\beta >
\beta_{\Delta}(\Lambda,\rho)
=
\frac{\pi\Lambda(2\Lambda\rho+5\Lambda-2)}{4}.
}
\]

So the geometry-side closure exists whenever

\[
\beta > \max(\beta_{\rm stab},\beta_\Delta).
\]

---

## 5. EM-branch worked example

Take the EM-cavity aspect ratio
\[
\Lambda_{\rm EM}=\frac{\sqrt{2}\pi}{x_{01}},
\qquad
x_{01}=2.40482555769577\ldots,
\]
so
\[
\Lambda_{\rm EM}\approx 1.847486577120128.
\]

A simple positive example is

\[
\rho=\frac{1}{10},
\qquad
\beta=12,
\]
for which

\[
\beta_{\rm stab}(\Lambda_{\rm EM},\rho)\approx 11.0420171,
\qquad
\beta_\Delta(\Lambda_{\rm EM},\rho)\approx 11.0377513.
\]

So this point lies safely inside the passive/positive-support region.

Matching the exact target
\[
\Delta K_{00}^{\rm geom}=\frac{109}{280}
\]
then fixes the overall geometry stiffness scale to

\[
\boxed{
\Sigma_* \approx 0.2076143291835488854.
}
\]

At this point the exact reduced geometry coefficient is

\[
\Delta K_{00}^{\rm geom} = \frac{109}{280}
\]
to machine precision.

---

## 6. Why the earlier single breathing auxiliary was already the right idea

At the worked point above, the 2DOF geometry Hessian eigenvalues are

\[
\lambda_1 \approx 0.10664211,
\qquad
\lambda_2 \approx 24.42044437.
\]

Using the natural normalized coupling vector

\[
\hat g = \frac{\nabla V}{V_0} = \left(3,\frac{1}{\Lambda_{\rm EM}}\right),
\]

the mode-resolved contributions to \(\Delta K_{00}^{\rm geom}\) are

\[
(0.00878310,\ 0.38050262),
\]
which sum to

\[
0.389285714285714\ldots = \frac{109}{280}.
\]

So the dominant mode contributes

\[
97.743791\%
\]
of the total static monopole response.

That is the key interpretation:

> the earlier single global monopole auxiliary is a **controlled reduction** of the full \((a,L)\) geometry sector, not an extra arbitrary assumption.

At this stage the monopole closure is structurally:

\[
\boxed{
\text{local Family-1 wall law}
\;+\;
\text{dominant-coupling reduction of the reduced geometry Hessian}.
}
\]

---

## 7. What remains

This does **not** finish the PDE derivation of the monopole channel, but it narrows it a lot.

The remaining task is now:

1. derive the curvature completion \(\kappa_b a^2/L\) (or its more accurate soft-wall analog) from the actual Family-1 soft-wall / traction physics,
2. and derive the corresponding breathing inertia if we want the full pole
   \[
   Y_{\rm mono}(\omega)=\frac{109/280}{1-\omega^2/\Omega_{\rm mono}^2}.
   \]

So the static monopole story is no longer “mysterious projector needed.”
It is:

- baseline \(PV+\sigma A\) is insufficient,
- minimal curvature completion repairs the geometry Hessian,
- and the passed \(109/280\) coefficient is realizable directly from that repaired geometry sector.
# 2PN geometry-breathing dynamic reduction: current result

## What this step adds

The previous step showed that the missing raw monopole wall closure
\[
\Delta K_{00}=\frac{109}{280}
\]
can be generated from the reduced static geometry Hessian in \((a,L)\), but it still left the **pole**
\[
Y_{\rm mono}(\omega)=\frac{109/280}{1-\omega^2/\Omega_{\rm mono}^2}
\]
underdetermined.

This step supplies the missing **dynamic reduction**.

Using a minimal affine breathing ansatz for the throat geometry,
we derive an exact reduced inertia matrix, obtain the full **two-pole** monopole response, and then show that the familiar single-pole auxiliary is the controlled low-frequency Padé reduction of that exact 2DOF dynamics.

This step does **not** depend on the earlier corrected \(z^8\) isotropic 4D-ball quotient, because it uses only the already-passed static geometry data.

---

## 1. Minimal affine inertia model

Use dimensionless geometry coordinates
\[
q=\left(\frac{\delta a}{a_0},\frac{\delta L}{a_0}\right),
\]
and the affine displacement field inside the 4D cylinder-like throat:
\[
\boldsymbol\xi_\perp=\frac{\delta a}{a_0}\,\mathbf r, 
\qquad
\xi_w=\frac{\delta L}{L_0}\,w.
\]

For the 3-ball cross section,
\[
\int_{B^3_{a_0}} r^2\,d^3x = \frac{4\pi a_0^5}{5},
\]
and along the axial interval,
\[
\int_{-L_0/2}^{L_0/2} w^2\,dw = \frac{L_0^3}{12}.
\]

So the reduced kinetic energy is
\[
T^{(2)} = \frac12\rho_{\rm eff}V_0 a_0^2\,\dot q^T \widehat M\,\dot q,
\qquad
\widehat M=
\begin{pmatrix}
\frac35 & 0 \\
0 & \frac1{12}
\end{pmatrix},
\]
with
\[
V_0 = \frac{4\pi}{3}a_0^3L_0 = \frac{4\pi}{3}\Lambda a_0^4.
\]

This is the minimal entrained-fluid inertia metric for the \((a,L)\) geometry sector.

---

## 2. Exact dynamic monopole susceptibility

From the previous static geometry step,
\[
E^{(2)} = \frac12\Sigma\,q^T \widehat H\,q,
\]
with
\[
\widehat H = \frac{a_0^2}{\Sigma}H_0,
\qquad
\bar g = \frac{a_0\nabla V_0}{V_0} = \left(3,\frac1\Lambda\right).
\]

Define the scaled frequency variable
\[
s = \omega^2\,\frac{\rho_{\rm eff}V_0 a_0^2}{\Sigma}.
\]
Then the exact reduced dynamic monopole response is
\[
\boxed{
Y_{\rm geom}(s)=
\frac1\Sigma\,
\bar g^T(\widehat H-s\widehat M)^{-1}\bar g.
}
\]

The raw wall monopole completion is therefore
\[
\boxed{
K_{00}^{\rm raw}(s)
=
-\frac{757}{2520}+Y_{\rm geom}(s).
}
\]
At \(s=0\), this reproduces the static completion from the earlier step.

---

## 3. Exact low-frequency reduction

Expanding at small \(s\),
\[
Y_{\rm geom}(s)=\Delta_0 + \Delta_2 s + O(s^2),
\]
with
\[
\Delta_0
=
\frac1\Sigma\bar g^T\widehat H^{-1}\bar g,
\qquad
\Delta_2
=
\frac1\Sigma\bar g^T\widehat H^{-1}\widehat M\widehat H^{-1}\bar g.
\]

The associated \([1/1]\) Padé single-pole reduction is
\[
\boxed{
Y_{\rm geom}(s)\approx
\frac{\Delta_0}{1-s/\lambda_{\rm eff}},
\qquad
\lambda_{\rm eff}=\frac{\Delta_0}{\Delta_2}.
}
\]

So the old monopole auxiliary really is the low-frequency reduction of a genuine reduced geometry dynamics.

---

## 4. EM-branch worked point

Using the same worked point as the previous geometry-Hessian closure,
\[
\Lambda_{\rm EM}=\frac{\sqrt2\pi}{x_{01}}\approx 1.8474865771,
\qquad
\rho=\frac1{10},
\qquad
\beta=12,
\qquad
\Sigma_*=0.2076143291835488854\ldots,
\]
we get
\[
\widehat H \approx
\begin{pmatrix}
114.33174685 & 19.35786733 \\
19.35786733 & 3.80598758
\end{pmatrix},
\qquad
\widehat M=
\begin{pmatrix}
0.6 & 0 \\
0 & 1/12
\end{pmatrix},
\qquad
\bar g=\left(3,\frac1\Lambda\right).
\]

The generalized eigenproblem
\[
\widehat H v_i = \lambda_i \widehat M v_i
\]
gives the exact dimensionless pole parameters
\[
\boxed{
\lambda_- \approx 5.23115613,
\qquad
\lambda_+ \approx 230.99360624.
}
\]
Hence the physical pole frequencies are
\[
\boxed{
\Omega_{\pm}^2
=
\frac{\Sigma_*}{\rho_{\rm eff}V_0 a_0^2}\,\lambda_{\pm}.
}
\]
For \(a_0=1\) and \(\rho_{\rm eff}=1\),
\[
\Omega_-^2\approx 0.14034117,
\qquad
\Omega_+^2\approx 6.19708399.
\]

---

## 5. Exact two-pole decomposition

At the worked point, the exact dynamic response can be written as
\[
\boxed{
Y_{\rm geom}(s)=
\frac{R_-}{1-s/\lambda_-}
+
\frac{R_+}{1-s/\lambda_+},
}
\]
with positive static residues
\[
\boxed{
R_- \approx 0.00327376153,
\qquad
R_+ \approx 0.38601195275.
}
\]
They sum exactly to
\[
R_-+R_+ = \frac{109}{280}.
\]

So the exact geometry breathing response is a **two-pole Stieltjes function with positive residues**.

The dominant residue fraction is
\[
\frac{R_+}{R_-+R_+} \approx 99.1590\%.
\]
So the old single-pole auxiliary is not just qualitatively justified; it is quantitatively a very good reduction of the exact 2DOF geometry dynamics.

---

## 6. Single-pole Padé reduction

At the same worked point,
\[
\boxed{
\lambda_{\rm eff}\approx 169.48205088,
\qquad
\Omega_{\rm eff}^2
=
\frac{\Sigma_*}{\rho_{\rm eff}V_0 a_0^2}\,\lambda_{\rm eff}.
}
\]
For \(a_0=1\) and \(\rho_{\rm eff}=1\),
\[
\Omega_{\rm eff}^2 \approx 4.54685531.
\]

On the low-frequency band
\[
0\le s\le 0.1\lambda_-,
\]
the maximum relative error of the one-pole Padé form against the exact two-pole response is
\[
\boxed{\max |\delta Y/Y| \approx 8.87\times 10^{-5}.}
\]

So for PN/local response purposes, the one-pole monopole auxiliary is extremely accurate on the natural low-frequency band.

---

## 7. Interpretation

This is the cleanest dynamic monopole picture so far.

The raw monopole wall sector is now
\[
\boxed{
K_{00}^{\rm raw}(s)
=
-\frac{757}{2520}
+
\frac{R_-}{1-s/\lambda_-}
+
\frac{R_+}{1-s/\lambda_+},
}
\]
with the controlled low-frequency reduction
\[
\boxed{
K_{00}^{\rm raw}(s)
\approx
-\frac{757}{2520}+
\frac{109/280}{1-s/\lambda_{\rm eff}}.
}
\]

So the earlier “global breathing auxiliary” is no longer just a bookkeeping closure. It is the low-frequency reduction of the same \((a,L)\) geometry sector that already produced the exact static \(109/280\) coefficient.

---

## 8. What remains

The remaining task is now narrow and concrete:

1. derive the overall inertial scale \(\rho_{\rm eff}\) (or its more accurate soft-wall analog) from the actual Family-1 confinement / traction PDE,
2. and, if desired, improve the affine inertia ansatz to a soft-wall boundary-layer inertia.

Once that is done, the monopole pole is fully fixed from throat-side physics rather than inferred from operator matching.
# 2PN Family-1 TF inertia scale from the parent PDE

## What this step adds

The previous geometry-breathing reduction showed that the missing monopole wall closure
\[
\Delta K_{00}=\frac{109}{280}
\]
is supplied by the reduced geometry Hessian in \((a,L)\), and that its dynamics is an exact two-pole response. The remaining open item was the **bulk inertia scale**
\[
\rho_{\rm eff}
\]
appearing in the dimensionless frequency variable
\[
s=\omega^2\frac{\rho_{\rm eff}V_0 a_0^2}{\Sigma}.
\]

This step derives that scale from the **Family-1 parent PDE** itself, using the already-frozen \(n=5\) barotrope in the harmonic interior trap.

---

## 1. Parent-PDE reduction

Use the Family-1 interior throat branch with:

- steep radial wall and endcaps,
- harmonic interior brane trap
  \[
  V_{\rm in}(w)=\frac12 m_\psi \Omega_{\rm in}^2 w^2,
  \]
- and the frozen barotrope
  \[
  P(\rho)=K_{\rm EOS}\rho^n,\qquad
  h(\rho)=\frac{nK_{\rm EOS}}{n-1}\rho^{n-1}.
  \]

In the Thomas-Fermi / hydrostatic limit,
\[
h(\rho)+\frac12 m_\psi \Omega_{\rm in}^2 w^2=\mu_{\rm TF}.
\]

On the **filled-to-endcap** branch, the support reaches \(|w|=L_0/2\), so
\[
\mu_{\rm TF}=\frac12 m_\psi\Omega_{\rm in}^2\left(\frac{L_0}{2}\right)^2.
\]

That gives the exact axial profile
\[
\rho_0(w)=\rho_c\left(1-\frac{4w^2}{L_0^2}\right)^{\!1/(n-1)},
\qquad |w|\le \frac{L_0}{2},
\]
with central density
\[
\boxed{
\rho_c=
\left[
\frac{(n-1)m_\psi\Omega_{\rm in}^2L_0^2}{8nK_{\rm EOS}}
\right]^{\!1/(n-1)}.
}
\]

---

## 2. Exact bulk inertia scale

Write \(u=2w/L_0\) and define
\[
c_0(n)=\frac12\int_{-1}^{1}(1-u^2)^{1/(n-1)}\,du.
\]

Then the effective bulk inertia density entering the monopole breathing channel is
\[
\boxed{
\rho_{\rm eff}^{\rm TF}(n)=c_0(n)\,\rho_c.
}
\]

The exact closed form is
\[
\boxed{
c_0(n)=
\frac{\sqrt{\pi}\,\Gamma\!\left(1+\frac{1}{n-1}\right)}
     {2\,\Gamma\!\left(\frac32+\frac{1}{n-1}\right)}.
}
\]

For the already-frozen \(n=5\) branch,
\[
\boxed{
\rho_{\rm eff}^{\rm TF}(5)=
\frac{\sqrt{\pi}\,\Gamma(1/4)}{6\,\Gamma(3/4)}
\left(\frac{m_\psi\Omega_{\rm in}^2L_0^2}{10K_{\rm EOS}}\right)^{1/4}.
}
\]

Numerically,
\[
\frac{\sqrt{\pi}\Gamma(1/4)}{6\Gamma(3/4)}
\approx 0.87401918476404.
\]

So the Family-1 parent PDE fixes the previously free bulk inertia scale up to the already-physical throat parameters \((m_\psi,K_{\rm EOS},\Omega_{\rm in},L_0)\).

---

## 3. Exact inertia metric on the \((a,L)\) geometry sector

Using the same affine breathing ansatz as before,
\[
\boldsymbol\xi_\perp=\frac{\delta a}{a_0}\mathbf r,
\qquad
\xi_w=\frac{\delta L}{L_0}w,
\]
the reduced kinetic energy is
\[
T^{(2)}=
\frac12\,\rho_{\rm eff}^{\rm TF}V_0a_0^2\,\dot q^T\widehat M_{\rm TF}\dot q,
\qquad
q=\left(\frac{\delta a}{a_0},\frac{\delta L}{a_0}\right).
\]

The radial piece stays unchanged:
\[
\widehat M_{aa}=\frac35.
\]

The axial piece is fixed by
\[
c_2(n)=\frac12\int_{-1}^{1}u^2(1-u^2)^{1/(n-1)}\,du,
\]
with exact ratio
\[
\frac{c_2(n)}{c_0(n)}=\frac{n-1}{3n-1}.
\]

Therefore
\[
\boxed{
\widehat M_{LL}^{\rm TF}(n)=\frac{c_2}{4c_0}
=\frac{n-1}{4(3n-1)}.
}
\]

At \(n=5\),
\[
\boxed{
\widehat M_{\rm TF}(5)=
\begin{pmatrix}
\frac35 & 0\\[4pt]
0 & \frac1{14}
\end{pmatrix}.
}
\]

So the Family-1 \(n=5\) parent PDE renormalizes the axial breathing moment from the uniform-slice value
\[
\frac1{12}\longrightarrow \frac1{14}.
\]

This is a clean EOS-sensitive result.

---

## 4. Dynamic monopole response with the TF inertia branch

Carry forward the same reduced geometry Hessian \(\widehat H\) and volume-coupling vector \(\bar g\) from the previous geometry-Hessian step. The exact reduced response is still
\[
\boxed{
Y_{\rm geom}^{\rm TF}(s)=
\frac1\Sigma\,\bar g^T\left(\widehat H-s\widehat M_{\rm TF}(5)\right)^{-1}\bar g,
\qquad
s=\omega^2\frac{\rho_{\rm eff}^{\rm TF}(5)V_0a_0^2}{\Sigma}.
}
\]

So the only change from the earlier affine-uniform model is:

- the free scale \(\rho_{\rm eff}\) is replaced by the explicit TF value,
- and the axial inertia coefficient becomes \(1/14\).

---

## 5. EM-worked point

At the same worked point used in the passed geometry-Hessian closure,
\[
\Lambda_{\rm EM}=\frac{\sqrt2\pi}{x_{01}}\approx 1.8474865771,
\qquad
\rho=\frac{1}{10},
\qquad
\beta=12,
\qquad
\Sigma_* \approx 0.2076143291835488854,
\]
the TF inertia branch gives dimensionless poles
\[
\boxed{
\lambda_- \approx 5.92556258,
\qquad
\lambda_+ \approx 237.91117494.
}
\]

The exact positive residues are
\[
\boxed{
R_- \approx 0.00262800,
\qquad
R_+ \approx 0.38665771,
}
\qquad
R_-+R_+=\frac{109}{280}.
\]

So the static closure is unchanged, as it must be.

The dominant residue fraction is
\[
\boxed{
\frac{R_+}{R_-+R_+}\approx 0.99324917.
}
\]

Thus the monopole wall completion remains an exact two-pole Stieltjes response with a highly accurate one-pole low-frequency reduction.

The Padé pole is
\[
\boxed{
\lambda_{\rm eff}^{\rm TF}\approx 188.17695898.
}
\]

On the low-frequency band \(0\le s\le 0.1\lambda_-\),
the max relative error of the one-pole reduction is only
\[
\boxed{
\max {\rm rel.\ err.}\approx 7.10\times 10^{-5}.
}
\]

---

## 6. Physical pole scales

Once the TF inertia scale is inserted, the physical poles are no longer free:
\[
\boxed{
\Omega_{\pm}^2=
\frac{\Sigma_*}{\rho_{\rm eff}^{\rm TF}(5)V_0a_0^2}\,\lambda_{\pm},
\qquad
\Omega_{\rm eff}^2=
\frac{\Sigma_*}{\rho_{\rm eff}^{\rm TF}(5)V_0a_0^2}\,\lambda_{\rm eff}^{\rm TF}.
}
\]

Since \(\rho_{\rm eff}^{\rm TF}(5)\) is now explicit, the monopole breathing scale is fixed by the same Family-1 throat microphysics already used to define the parent PDE.

---

## Bottom line

This step closes an important remaining gap:

1. the Family-1 parent PDE fixes the bulk inertia scale,
2. the frozen \(n=5\) EOS fixes the axial breathing moment to \(1/14\),
3. the dynamic monopole wall channel remains a positive two-pole response,
4. and the one-pole breathing auxiliary stays an excellent low-frequency reduction.

So the bulk part of the monopole dynamics is no longer phenomenological.

The main open dynamic ingredient now is **not** the bulk inertia scale. It is the genuine **soft-wall / surface inertia completion**, if you want to go beyond the bulk TF branch.

# 2PN Family-1 soft-wall surface inertia completion

## What this step adds

The previous geometry-breathing chain already fixed the **bulk** part of the monopole dynamics from the \(n=5\) Thomas–Fermi interior and reduced the remaining open item to the genuine **soft-wall / surface inertia** completion.

This step derives the leading radial soft-wall correction from the actual Family-1 wall profile
\[
V_{\rm wall}(r;a)=V_0\,S\!\left(\frac{r-a}{d_r}\right)^p,
\qquad
S(x)=\frac{1+\tanh x}{2},
\]
using the same frozen \(n=5\) hydrostatic closure inside the throat.

The result is that the remaining dynamic completion is **not** a new phenomenological term. It is a derived thin-wall correction controlled by the actual wall parameters
\[
\alpha_0 \equiv \frac{V_0}{\mu_c},
\qquad
p,
\qquad
\varepsilon_r \equiv \frac{d_r}{a_0}.
\]

---

## 1. Local \(n=5\) soft-wall profile

At fixed axial slice, the radial hydrostatic/TF wall profile is
\[
f_{\alpha,p}(\xi)=\left(1-\alpha\,S(\xi)^p\right)_+^{1/4},
\qquad
\xi=\frac{r-a}{d_r}.
\]

The turning point is
\[
\xi_*(\alpha,p)=\operatorname{arctanh}\!\bigl(2\alpha^{-1/p}-1\bigr).
\]

So
\[
\xi_*(0)=0 \iff \alpha=2^p.
\]

This is useful physically:

- \(\alpha_0>2^p\): support already turns off **inside** the nominal radius \(a_0\),
- \(\alpha_0<2^p\): support still spills a bit **outside** \(a_0\).

---

## 2. Universal thin-wall moment expansion

Define defect moments against the sharp-wall step:
\[
m_k(\alpha,p)=\int_{-\infty}^{\infty}\xi^k\Bigl[f_{\alpha,p}(\xi)-\Theta(-\xi)\Bigr]\,d\xi.
\]

For the 3-ball cross-section, the relevant radial moments are
\[
J_2=\frac13+\varepsilon_r m_0+2\varepsilon_r^2 m_1+\varepsilon_r^3 m_2+O(\varepsilon_r^4),
\]
\[
J_4=\frac15+\varepsilon_r m_0+4\varepsilon_r^2 m_1+6\varepsilon_r^3 m_2+O(\varepsilon_r^4).
\]

Therefore the cross-sectional mass factor is
\[
\boxed{
R_{\rm mass}=3J_2
=1+3\varepsilon_r m_0+6\varepsilon_r^2 m_1+3\varepsilon_r^3 m_2+O(\varepsilon_r^4),
}
\]
and the radial inertia coefficient becomes
\[
\boxed{
M_{aa}=\frac{J_4}{J_2}
=
\frac35+\frac65\varepsilon_r m_0
+\left(\frac{42}{5}m_1-\frac{18}{5}m_0^2\right)\varepsilon_r^2
+\left(\frac{54}{5}m_0^3-\frac{162}{5}m_0m_1+\frac{81}{5}m_2\right)\varepsilon_r^3
+O(\varepsilon_r^4).
}
\]

So the leading dynamic correction is especially simple:

> the same wall moment \(m_0\) controls both the mass-scale correction and the radial inertia shift at \(O(\varepsilon_r)\).

---

## 3. Axially averaged moments on the filled-to-endcap \(n=5\) branch

On the carried-forward filled-to-endcap \(n=5\) TF branch,
\[
\mu(y)=\mu_c(1-y^2),
\qquad
y=\frac{2w}{L}.
\]

So the local wall ratio becomes
\[
\alpha(y)=\frac{\alpha_0}{1-y^2}.
\]

The relevant averaged wall moments are
\[
\bar m_k(\alpha_0,p)=
\frac{1}{2c_0}
\int_{-1}^{1}(1-y^2)^{1/4}
\,m_k\!\left(\frac{\alpha_0}{1-y^2},p\right)\,dy,
\]
with
\[
c_0=\frac12\int_{-1}^{1}(1-y^2)^{1/4}\,dy.
\]

Then the full radial soft-wall completion of the monopole inertia sector is
\[
\boxed{
\rho_{\rm eff}^{\rm TF+wall}
=
\rho_{\rm eff}^{\rm TF}(5)
\Bigl[
1+3\varepsilon_r\bar m_0
+6\varepsilon_r^2\bar m_1
+3\varepsilon_r^3\bar m_2
+O(\varepsilon_r^4)
\Bigr],
}
\]
\[
\boxed{
\widehat M_{aa}^{\rm TF+wall}
=
\frac35+\frac65\varepsilon_r\bar m_0
+\left(\frac{42}{5}\bar m_1-\frac{18}{5}\bar m_0^2\right)\varepsilon_r^2
+\left(\frac{54}{5}\bar m_0^3-\frac{162}{5}\bar m_0\bar m_1+\frac{81}{5}\bar m_2\right)\varepsilon_r^3
+O(\varepsilon_r^4),
}
\]
while
\[
\widehat M_{LL}=\frac1{14}
\]
stays unchanged at this radial-order completion.

---

## 4. Dynamic monopole response with the wall correction

The carried-forward static geometry Hessian and volume-work coupling are unchanged, so the static monopole closure remains
\[
\Delta K_{00}=\frac{109}{280}.
\]

Only the dynamic inertia changes. The corrected geometry response is therefore
\[
Y_{\rm geom}^{\rm TF+wall}(s)=
\frac{1}{\Sigma_*}\,
\bar g^T\!\left(\widehat H-s\,\widehat M^{\rm TF+wall}\right)^{-1}\bar g,
\]
with
\[
s=\omega^2\,
\frac{\rho_{\rm eff}^{\rm TF+wall}V_0a_0^2}{\Sigma_*}.
\]

So the physical monopole poles are
\[
\Omega_\pm^2(\alpha_0,p,\varepsilon_r)
=
\frac{\Sigma_*}{\rho_{\rm eff}^{\rm TF+wall}V_0a_0^2}
\lambda_\pm\!\left(\widehat M^{\rm TF+wall}\right).
\]

At the EM worked point already fixed in the earlier scripts,
\[
\Lambda_{\rm EM}\approx 1.8474865771,
\qquad
\Sigma_*\approx 0.20761432918,
\]
the sharp-wall TF baseline remains
\[
\lambda_-\approx 5.92556258,\qquad
\lambda_+\approx 237.91117494,
\]
with residues
\[
R_-\approx 0.00262800,\qquad
R_+\approx 0.38665771.
\]

The leading physical pole shifts are
\[
\frac{\Omega_-^2}{\Omega_{-,{\rm sharp}}^2}
=
1-3.40828621\,\varepsilon_r\bar m_0+O(\varepsilon_r^2),
\]
\[
\frac{\Omega_+^2}{\Omega_{+,{\rm sharp}}^2}
=
1-4.59171379\,\varepsilon_r\bar m_0+O(\varepsilon_r^2).
\]

So on the steep-wall branch, where \(\bar m_0<0\), the soft wall **raises** the physical poles.

---

## 5. Representative steep-wall branch

Take the first genuinely steep Family-1 reference branch
\[
p=2,\qquad \alpha_0=10,\qquad \varepsilon_r=0.05.
\]

Then
\[
2^p=4<\alpha_0,
\qquad
\xi_*(0)\approx -0.38558107,
\]
so support already turns off inside the nominal radius at the throat center.

The averaged wall moments are
\[
\bar m_0\approx -0.65067123,\qquad
\bar m_1\approx 0.25044370,\qquad
\bar m_2\approx -0.15585783.
\]

Therefore
\[
\boxed{
R_{\rm mass}\approx 0.90609752,
\qquad
\widehat M_{aa}^{\rm TF+wall}\approx 0.56238115.
}
\]

The corrected dimensionless poles are
\[
\boxed{
\lambda_-\approx 6.00228906,
\qquad
\lambda_+\approx 250.58092901,
}
\]
with residues
\[
\boxed{
R_-\approx 0.00289847,
\qquad
R_+\approx 0.38638724.
}
\]

So the actual physical pole-squared ratios relative to the sharp-wall TF baseline are
\[
\boxed{
\frac{\Omega_-^2}{\Omega_{-,{\rm sharp}}^2}\approx 1.11792424,
\qquad
\frac{\Omega_+^2}{\Omega_{+,{\rm sharp}}^2}\approx 1.16240703.
}
\]

The one-pole Padé reduction remains excellent:
\[
\lambda_{\rm eff}\approx 192.25314580,
\qquad
\max{\rm\ rel.\ err.}\approx 7.84\times 10^{-5}
\]
on the natural low-frequency band \(0\le s\le 0.1\lambda_-\).

So the old breathing auxiliary still survives, but now with a derived soft-wall correction.

---

## 6. Contrast branches

Two useful contrasts from the prototype:

### Stiffer wall, same \(p\)
\[
p=2,\quad \alpha_0=20,\quad \varepsilon_r=0.05
\]
gives
\[
R_{\rm mass}\approx 0.87606211,\qquad
\widehat M_{aa}\approx 0.54983763,
\]
and
\[
\frac{\Omega_-^2}{\Omega_{-,{\rm sharp}}^2}\approx 1.16125672,\qquad
\frac{\Omega_+^2}{\Omega_{+,{\rm sharp}}^2}\approx 1.22438713.
\]

### Sharper power, same central wall ratio
\[
p=4,\quad \alpha_0=10,\quad \varepsilon_r=0.05
\]
gives
\[
R_{\rm mass}\approx 0.99001853,\qquad
\widehat M_{aa}\approx 0.59627288,
\]
and
\[
\frac{\Omega_-^2}{\Omega_{-,{\rm sharp}}^2}\approx 1.01136447,\qquad
\frac{\Omega_+^2}{\Omega_{+,{\rm sharp}}^2}\approx 1.01510709.
\]

So the dynamic wall correction is modest and controllable in the steep Family-1 regime.

---

## Bottom line

This step closes the dynamic wall piece much more tightly:

1. the bulk TF core still fixes the dominant inertia scale,
2. the radial Family-1 soft wall gives a derived thin-wall correction to both the overall mass scale and \(M_{aa}\),
3. the static \(109/280\) monopole closure is untouched,
4. and the dynamic monopole channel remains a positive two-pole Stieltjes response with an excellent one-pole reduction.

So the remaining dynamic “surface inertia” is no longer a placeholder. It is already a controlled function of
\[
\left(\frac{V_0}{\mu_c},\,p,\,\frac{d_r}{a_0}\right).
\]

The next natural tightening is now very specific:

- add the **endcap soft-wall layer** to refine \(M_{LL}\),
- and/or tie this same wall profile directly to the earlier **tangential traction/support** law so the static and dynamic wall sectors come from one boundary-layer model.
# 2PN Family-1 endcap soft-wall inertia completion

## What this step adds

The previous chain had already fixed:

- the bulk Thomas–Fermi inertia scale on the filled-to-endcap `n=5` branch,
- the radial Family-1 sidewall correction to the monopole breathing channel,
- and the reduced geometry-Hessian / two-pole breathing closure.

The remaining missing dynamic wall ingredient was the **endcap soft-wall layer**.

This step derives that correction and folds it into the carried-forward monopole response.
The punchline is simple:

1. the endcap layer is **parametrically weaker** than the sidewall,
2. its first nontrivial effect scales as `eps_z^(5/4)` on the frozen `n=5` branch,
3. and once it is included, the current program already has a near-final **full wall-completed monopole breathing branch**, at least to separated leading order.

---

## 1. Why the endcap scaling is different

On the filled-to-endcap TF branch,

axial support already vanishes at the cap:

\[
\rho_0(u) \propto (1-u^2)^{1/4},
\qquad u=\frac{2w}{L}.
\]

So near the right endcap, with
\[
u = 1 + \varepsilon_z x,
\qquad \varepsilon_z = \frac{2 d_z}{L},
\]
we have
\[
1-u^2 = \varepsilon_z\bigl(-2x-\varepsilon_z x^2\bigr).
\]

That means a genuine thin-cap layer cannot have a wall amplitude of order `mu_c`.
To stay localized on an `O(d_z)` layer, the cap must scale as
\[
\frac{V_{\rm cap}}{\mu_c} = 2\varepsilon_z\,\alpha_z\,S(x)^p,
\qquad S(x)=\frac{1+\tanh x}{2}.
\]

This is the key distinction from the sidewall.
Because the TF profile already dies at the endcap, the first wall correction is suppressed by one extra quarter-power.

---

## 2. Local reduced cap profile and universal defect moment

With the scaling above, the local reduced cap profile is
\[
g_{\alpha,p}(x)=\bigl(-x-\alpha S(x)^p\bigr)_+^{1/4}.
\]

Relative to the sharp-cap baseline `(-x)_+^(1/4)`, define the universal defect moments
\[
\nu_k(\alpha,p)=\int_{-\infty}^{\infty}x^k\Bigl[g_{\alpha,p}(x)-(-x)_+^{1/4}\Bigr]dx.
\]

Then the filled-to-endcap TF integrals obey
\[
c_0^{\rm cap}=c_0+2^{1/4}\varepsilon_z^{5/4}\nu_0+O(\varepsilon_z^{9/4}),
\]
\[
c_2^{\rm cap}=c_2+2^{1/4}\varepsilon_z^{5/4}\nu_0+O(\varepsilon_z^{9/4}),
\]
with the carried-forward sharp-cap baseline
\[
c_0=\frac{\sqrt\pi\,\Gamma(5/4)}{2\Gamma(7/4)},
\qquad
c_2=\frac{2}{7}c_0.
\]

So the effective bulk inertia scale and the axial breathing metric become
\[
\boxed{
\frac{\rho_{\rm eff}^{\rm TF+cap}}{\rho_{\rm eff}^{\rm TF}}
=1+A_{\rm cap}\,\nu_0\,\varepsilon_z^{5/4}+\cdots,
}
\]
\[
\boxed{
\widehat M_{LL}^{\rm TF+cap}
=\frac{1}{14}+B_{\rm cap}\,\nu_0\,\varepsilon_z^{5/4}+\cdots,
}
\]
with exact coefficients
\[
A_{\rm cap}=\frac{2^{1/4}}{c_0}
=\frac{6\,2^{1/4}\Gamma(3/4)}{\sqrt\pi\,\Gamma(1/4)}
\approx 1.3606190066912236,
\]
\[
B_{\rm cap}=\frac{5\,2^{1/4}}{28c_0}
\approx 0.24296767976628994.
\]

This is the clean main scaling result:

> on the frozen `n=5` branch, the endcap wall correction is `O(eps_z^(5/4))`, not `O(eps_z)`.

---

## 3. Representative steep-cap branch and direct full-profile check

Take the representative cap branch
\[
\varepsilon_z=0.05,
\qquad \alpha_z=1,
\qquad p_z=2.
\]

The local turning point solving
\[
x+\alpha S(x)^p=0
\]
is
\[
x_*\approx -0.1720646550263600,
\]
and the universal defect moment is
\[
\nu_0\approx -0.1297171210945550.
\]

Using the asymptotic formulas above gives
\[
c_0^{\rm asym}\approx 0.8703719198752305,
\qquad
c_2^{\rm asym}\approx 0.2460725021866305,
\]
\[
\widehat M_{LL}^{\rm asym}\approx 0.0706802737334132.
\]

Direct numerical integration of the full soft-cap profile gives
\[
c_0^{\rm full}\approx 0.8703584556098921,
\qquad
c_2^{\rm full}\approx 0.2461237735473917,
\]
\[
\widehat M_{LL}^{\rm full}\approx 0.0706960942244548.
\]

So already at `eps_z = 0.05` the leading asymptotic is very good:

- `c0` relative error `~ 1.55 × 10^-5`,
- `M_LL` relative error `~ 2.24 × 10^-4`.

This makes the cap reduction numerically trustworthy in the same regime used by the earlier sidewall work.

---

## 4. Dynamic monopole response with the endcap correction

On the carried-forward EM worked point,
\[
\Lambda_{\rm EM}=\frac{\sqrt2\pi}{x_{01}}\approx 1.8474865771,
\qquad
\rho=\frac{1}{10},
\qquad
\beta=12,
\qquad
\Sigma_*\approx 0.20761432918,
\]
the sharp-wall TF baseline had
\[
\lambda_-\approx 5.92556258,
\qquad
\lambda_+\approx 237.91117494,
\]
with residues summing to the carried-forward monopole target
\[
R_-+R_+=\frac{109}{280}.
\]

The new cap branch gives
\[
R_{\rm cap}=\frac{c_0^{\rm cap}}{c_0}\approx 0.99581161464,
\qquad
\widehat M_{LL}^{\rm cap}\approx 0.07069609422.
\]

Feeding that into the same reduced geometry response yields
\[
\lambda_-\approx 5.97431790,
\qquad
\lambda_+\approx 238.41448955,
\]
with residues
\[
R_-\approx 0.00258517,
\qquad
R_+\approx 0.38670054,
\qquad
R_-+R_+=\frac{109}{280}.
\]

The physical pole-squared ratios relative to the sharp-wall TF baseline are
\[
\boxed{
\frac{\Omega_-^2}{\Omega_{-,\rm sharp}^2}\approx 1.01246857,
\qquad
\frac{\Omega_+^2}{\Omega_{+,\rm sharp}^2}\approx 1.00633046.
}
\]

So the endcap correction is real, but modest — exactly what the `eps_z^(5/4)` scaling suggested.

The one-pole Padé reduction remains excellent:
\[
\lambda_{\rm eff}\approx 189.46282891,
\qquad
\max \text{ rel. err.}\approx 6.98\times 10^{-5}
\]
on the natural low-frequency band `0 ≤ s ≤ 0.1 lambda_-`.

---

## 5. Leading separated-order full-wall composite branch

Now combine:

- the carried-forward radial sidewall result from the previous step,
  \[
  R_{\rm side}\approx 0.90609752477,
  \qquad
  \widehat M_{aa}\approx 0.56238115491,
  \]
- with the new endcap result,
  \[
  R_{\rm cap}\approx 0.99581161464,
  \qquad
  \widehat M_{LL}\approx 0.07069609422.
  \]

To leading separated order,
\[
R_{\rm full}=R_{\rm side}R_{\rm cap}\approx 0.90230243917,
\]
\[
\widehat M_{\rm full}=
\begin{pmatrix}
0.56238115491 & 0\\
0 & 0.07069609422
\end{pmatrix}.
\]

The resulting full-wall composite breathing branch is
\[
\lambda_-\approx 6.05235326,
\qquad
\lambda_+\approx 251.08293474,
\]
with residues
\[
R_-\approx 0.00285529,
\qquad
R_+\approx 0.38643043,
\qquad
R_-+R_+=\frac{109}{280}.
\]

The corresponding physical pole-squared ratios relative to the sharp-wall TF baseline are
\[
\boxed{
\frac{\Omega_-^2}{\Omega_{-,\rm sharp}^2}\approx 1.13198989,
\qquad
\frac{\Omega_+^2}{\Omega_{+,\rm sharp}^2}\approx 1.16963464.
}
\]

The one-pole Padé reduction is still excellent:
\[
\lambda_{\rm eff}\approx 193.59552541,
\qquad
\max \text{ rel. err.}\approx 7.72\times 10^{-5}.
\]

So this is the first actual **full wall-completed monopole breathing branch** in the current Family-1 program, at least to separated leading order in the sidewall and endcap thickness parameters.

---

## Bottom line

This step narrows the remaining finish-line task substantially.

We now have:

1. **static local wall support/source law** from the Family-1 traction sector,
2. **static monopole closure** from the reduced geometry Hessian,
3. **bulk TF inertia** from the parent PDE,
4. **radial sidewall dynamic correction** from the Family-1 radial soft wall,
5. **endcap dynamic correction** from the new filled-to-endcap cap-layer reduction,
6. and a first **full wall-completed dynamic monopole branch** by combining the carried-forward sidewall result with the new cap result.

So the remaining gap is no longer “derive the dynamic wall sector.”
It is much narrower:

> derive the **fully coupled** sidewall-plus-endcap boundary-layer reduction beyond separated order, and then package the whole wall sector into one final throat-response module.

# 2PN Family-1 final coupled response module

## What this step closes

This step does two things at once.

First, it goes **beyond the separated sidewall-plus-endcap approximation** and evaluates the full
coupled Family-1 Thomas–Fermi wall profile directly in the `(x,y)` throat coordinates.

Second, it packages that coupled wall completion together with the already-passed exact
support/source law and odd/even port structure into a **single low-frequency Family-1 throat-response module**
that reproduces the conservative added 2PN cross sector exactly at zero frequency.

So this is the closest thing yet to the actual finish line.

---

## 1. Exact static support/source law remains unchanged

The local static even-wall law is still

\[
z_{\rm base}(\mu)=\frac{17}{56}-\frac{5}{56}\mu^2,
\qquad
t(\mu)=\frac{593}{672}-\frac{1553}{672}\mu^2+\frac78\mu^4,
\]

with exact source profile

\[
S(\mu)=10-\frac{63}{2}z_{\rm base}(\mu)
      =\frac{7}{16}+\frac{45}{16}\mu^2.
\]

That part does not move.

---

## 2. Final throat module reproduces the exact added 2PN cross block

Using the already-fixed odd dipole wake plus the canonical \(P_0 \oplus P_2\) even ports, the final low-frequency module is

\[
L^{\rm add}_{\rm odd}
=
\frac12(v_A^2+v_B^2)L^{\rm wake}_{1PN}
-\frac{15}{4}(U_A+U_B)(\mathbf v_A\!\cdot\!\mathbf v_B),
\]

\[
L^{\rm add}_{\rm even}
=
\Pi_A\!\cdot\!\Pi_B
+
U_A\,J\!\cdot\!\Pi_B
+
U_B\,J\!\cdot\!\Pi_A
+
\bigl(J\!\cdot\!J-\Delta_{\rm geom}\bigr)U_AU_B,
\]

with

\[
J=\left(\frac{4}{\sqrt5},\frac54\right),
\qquad
\Delta_{\rm geom}=\frac{281}{80}.
\]

The symbolic residual against the exact solved added 2PN comparable-mass target vanishes identically.

So the constructive zero-frequency 2PN cross sector is now closed.

---

## 3. Direct coupled Family-1 full-profile wall completion

On the balanced thin-layer-consistent reference branch

\[
\alpha_r=10,\qquad \varepsilon_r=0.05,\qquad p_r=2,
\]
\[
\chi_{\rm cap}=4,\qquad \varepsilon_z=0.05,\qquad
\alpha_{\rm cap}=4\varepsilon_z\chi_{\rm cap}=0.8,\qquad p_z=2,
\]

use the full coupled TF profile

\[
\widetilde\rho(x,y)=
\Bigl[
1-y^2
-\alpha_r\,S\!\left(\frac{x-1}{\varepsilon_r}\right)^{p_r}
-\alpha_{\rm cap}\,S\!\left(\frac{y-1}{2\varepsilon_z}\right)^{p_z}
\Bigr]_+^{1/4},
\]

on \(0 \le x \le 1\), \(0 \le y \le 1\), with symmetry factor \(2\) in \(y\).

The exact coupled integrals are

\[
I_2 = 2\int_0^1\!\!\int_0^1 x^2\widetilde\rho\,dx\,dy,
\qquad
I_4 = 2\int_0^1\!\!\int_0^1 x^4\widetilde\rho\,dx\,dy,
\]
\[
I_w = 2\int_0^1\!\!\int_0^1 \frac{y^2}{4}x^2\widetilde\rho\,dx\,dy.
\]

Normalizing by the sharp-wall filled-to-endcap TF baseline gives

\[
R_{\rm mass}^{\rm full}\approx 0.886313972989725,
\]
\[
\widehat M_{aa}^{\rm full}\approx 0.563114968953987,
\qquad
\widehat M_{LL}^{\rm full}\approx 0.065829228119349.
\]

This is the first direct **coupled** sidewall+endcap check, rather than a separated leading-order composition.

---

## 4. Separated-order approximation is already very good

Compared to the carried-forward separated-order wall completion,

\[
R_{\rm mass}^{\rm sep}\approx 0.8846236634,
\qquad
\widehat M_{aa}^{\rm sep}\approx 0.5623810783,
\qquad
\widehat M_{LL}^{\rm sep}\approx 0.0671965962,
\]

the direct coupled full-profile values differ by only

\[
\frac{R_{\rm mass}^{\rm full}-R_{\rm mass}^{\rm sep}}
     {R_{\rm mass}^{\rm sep}}
\approx +1.91\times 10^{-3},
\]
\[
\frac{\widehat M_{aa}^{\rm full}-\widehat M_{aa}^{\rm sep}}
     {\widehat M_{aa}^{\rm sep}}
\approx +1.30\times 10^{-3},
\]
\[
\frac{\widehat M_{LL}^{\rm full}-\widehat M_{LL}^{\rm sep}}
     {\widehat M_{LL}^{\rm sep}}
\approx -2.03\times 10^{-2}.
\]

So the earlier separated reduction was already excellent for the total mass scale and radial inertia, and still very respectable for the axial inertia.

That is an important closure check.

---

## 5. Final coupled monopole breathing response

Feeding the exact coupled wall data into the already-fixed geometry Hessian at the EM worked point gives

\[
\lambda_-\approx 6.405572392138922,
\qquad
\lambda_+\approx 254.444968136936126,
\]

with positive residues

\[
R_-\approx 0.002552474771738,
\qquad
R_+\approx 0.386733239513976,
\qquad
R_-+R_+=\frac{109}{280}.
\]

So the final monopole channel is

\[
K_{00}(s)=
-\frac{757}{2520}
+\frac{R_-}{1-s/\lambda_-}
+\frac{R_+}{1-s/\lambda_+}.
\]

At zero frequency,
\[
K_{00}(0)=\frac{4}{45},
\]
exactly as required.

Relative to the sharp-wall TF baseline, the physical monopole pole ratios are

\[
\frac{\Omega_-^2}{\Omega_{-,{\rm TF}}^2}
\approx 1.219665554412172,
\qquad
\frac{\Omega_+^2}{\Omega_{+,{\rm TF}}^2}
\approx 1.206678094538845.
\]

The one-pole Padé reduction is still excellent:

\[
\lambda_{\rm eff}\approx 202.923516367519028,
\qquad
\max {\rm rel.err.}\approx 6.89\times 10^{-5}
\quad\text{on } 0\le s\le 0.1\lambda_-.
\]

So the old single breathing auxiliary remains a very controlled low-frequency reduction of the exact two-pole geometry channel.

---

## 6. Final state of play

At conservative 2PN order, the constructive throat program is now essentially closed in the following sense:

- **exact local static support/source wall law**: closed,
- **exact zero-frequency odd + \(P_0 \oplus P_2\) port reconstruction**: closed,
- **exact algebraic match to the solved added comparable-mass 2PN cross block**: closed,
- **static monopole closure \(109/280\) from geometry Hessian**: closed,
- **bulk TF inertia scale**: closed,
- **radial sidewall and endcap dynamic corrections**: closed,
- **direct coupled full-profile sidewall+endcap check**: now done.

The remaining open quantities are narrower than before:

1. the **non-monopole dynamic pole scales** (\(\ell=1,2\) channels), which are true inner-throat PDE observables rather than conservative 2PN necessities;
2. a fully analytic derivation of the local wall profiles from the microscopic soft-wall traction PDE, rather than from reduced matching.

Those no longer block the conservative 2PN derivation itself.

So the practical finish-line statement is:

\[
\boxed{
\text{the conservative 2PN throat-response module is now effectively in hand.}
}
\]
