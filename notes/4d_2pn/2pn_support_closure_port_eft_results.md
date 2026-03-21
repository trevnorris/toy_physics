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
