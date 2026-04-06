
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
