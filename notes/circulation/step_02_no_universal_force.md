# Step 2 — 3D finite-mouth circulation does not by itself give a universal radial force

## Purpose

This step is the no-go part of the derivation. It asks what follows from a fixed 3D mouth circulation label \(n_A\) before imposing a current-loop or mixed-sector plumbing law.

## 3D difference from the 2D simulation picture

In 2D, a point vortex has a logarithmic Green function and it is tempting to read \(n_1n_2\) directly as a pair-potential sign. That intuition is not safe in the full 3D finite-mouth problem.

A localized 3D mouth/rim swirl is not a scalar vortex charge in the 2D sense. At large separation it is represented by a compact axial multipole: a circulation/current loop, vortex ring, or mixed-sector port with an orientation vector. Therefore the leading pair interaction is geometric and tensorial, not just a scalar function of \(n_1n_2\).

## Current-like force form, if a closure is supplied

If a later closure maps circulation into an effective localized current or dipole moment, write the effective global current sign as
\[
I_A = \alpha_A\,\sigma_A\,n_A.
\]
Here:

- \(n_A\) is the fluxoid/circulation integer;
- \(\alpha_A\) is a closure-dependent current/throughput coefficient;
- \(\sigma_A=\pm1\) encodes how the local mouth orientation maps into a common global axis.

For a coaxial fixed-current loop closure, the leading mutual inductance is
\[
M(d)=\frac{B}{d^3},\qquad B>0.
\]
The fixed-current mechanical potential is
\[
U_I(d)=-I_1I_2M(d),
\]
so the leading radial force is derived as
\[
F_d
=-\frac{\partial U_I}{\partial d}
= -\frac{3B}{d^4}\,\alpha_1\alpha_2\sigma_1\sigma_2 n_1n_2,
\]
where \(F_d<0\) means attraction.

This formula is useful, but it is **not** a consequence of the fluxoid constraint alone. It is a consequence of adding both \(M(d)\) and a current-like closure.

## Why the fluxoid alone is insufficient

The fluxoid constraint fixes
\[
\Gamma_A=\frac{2\pi\hbar}{m}n_A.
\]
It does **not** fix:

\[
\alpha_A,\qquad \sigma_A,\qquad I_A,\qquad J_A^w,
\qquad A_w/F_{\mu w}\ \text{transport},
\qquad \text{fixed-current versus fixed-flux ensemble}.
\]

Therefore the same pair of integers \((n_1,n_2)\) can be attractive or repulsive under different legitimate closures. For example, with \(n_1=n_2=+1\), \(\alpha_1\alpha_2>0\) attracts in a fixed-current coaxial closure, while \(\alpha_1\alpha_2<0\) repels.

There is also an ensemble issue. In ordinary magnetic systems, the field energy in current coordinates is
\[
E_B(I,d)=\frac12 I^T L(d)I.
\]
The fixed-current mechanical potential is the Legendre-transformed quantity
\[
G_I=E_B-I\cdot\Phi=-\frac12 I^TL(d)I,
\qquad
\Phi=L(d)I.
\]
Thus the mutual cross coefficient changes from \(+M\) in \(E_B\) to \(-M\) in \(G_I\). The mechanical force law must state which quantity is held fixed and which potential is being differentiated.

For fixed fluxes,
\[
U_\Phi=\frac12\Phi^T L(d)^{-1}\Phi.
\]
This is a separate ensemble expression, not a hand-flipped copy of the fixed-current potential.

At weak coupling with \(M(d)=B/d^3\), the fixed-flux force expands as
\[
F_\Phi
=-\frac{3B\Phi_1\Phi_2}{L_1L_2d^4}
+\frac{3B^2(L_1\Phi_2^2+L_2\Phi_1^2)}{L_1^2L_2^2d^7}
+O(B^3).
\]
The leading fixed-flux sign matches the fixed-current sign after the weak-coupling identification \(I\approx L^{-1}\Phi\). The point of the ensemble audit is not a universal leading sign flip; it is that the potential and held-fixed variables must be stated before differentiating.

## Step verdict

The fluxoid law is a holonomy law. It is not yet a radial-force law.

So the theorem status after Step 2 is:
\[
\boxed{\text{No universal attraction/repulsion sign follows from }n_1n_2\text{ alone in 3D.}}
\]

The next step supplies a named closure: the fixed-current Maxwell/current-loop closure.

## SymPy audit

The script derives the force from \(M(d)=B/d^3\), verifies that the same \((n_1,n_2)\) can flip force sign when the closure coefficient flips, performs the fixed-current Legendre transform from the inductance matrix, and expands the fixed-flux force through \(B^2\).

Run:

```bash
python step_02_no_universal_force_sympy.py
```
