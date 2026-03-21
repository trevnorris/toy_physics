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
