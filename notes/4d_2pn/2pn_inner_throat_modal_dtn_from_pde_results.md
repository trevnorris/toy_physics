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
