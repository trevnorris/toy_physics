# Step 5 — Full 3D orientation law and finite-mouth validity

## Purpose

Steps 3–4 used the coaxial loop case because it is the cleanest model for two facing mouths. This step records the full 3D dipole-orientation law and clarifies what the finite-mouth correction does and does not prove.

## 3D dipole orientation law

For a current-like dipole closure, define vector dipole moments
\[
\mathbf m_A = I_A \pi R_A^2\,\hat s_A.
\]
Let \(\hat d\) point from mouth 1 to mouth 2, and define
\[
a=\hat m_1\cdot\hat d,
\qquad
b=\hat m_2\cdot\hat d,
\qquad
c=\hat m_1\cdot\hat m_2.
\]
The fixed-current dipole potential is
\[
U(d)
=-\frac{\mu_0}{4\pi d^3}
\left[3(\mathbf m_1\cdot\hat d)(\mathbf m_2\cdot\hat d)-\mathbf m_1\cdot\mathbf m_2\right].
\]
The audit builds \(\mathbf m_1\) and \(\mathbf m_2\) as actual 3-vectors for each case, then evaluates the dot products from components. The scalar notation
\[
S=3ab-c
\]
is only shorthand for the component calculation, not an independent input.

Differentiating the fixed-current potential gives
\[
F_d=-\frac{\partial U}{\partial d}
=-\frac{3\mu_0}{4\pi d^4}
\left[3(\mathbf m_1\cdot\hat d)(\mathbf m_2\cdot\hat d)-\mathbf m_1\cdot\mathbf m_2\right].
\]
So attraction requires
\[
\boxed{S=3ab-c>0.}
\]

## Orientation examples

| Geometry | \(S\) | Result |
|---|---:|---|
| coaxial aligned | \(+2\) | attraction |
| coaxial anti-aligned | \(-2\) | repulsion |
| side-by-side parallel | \(-1\) | repulsion |
| side-by-side anti-parallel | \(+1\) | attraction |
| one axial, one transverse | \(0\) | no leading dipole radial force |

This is why a scalar phrase like “same circulation attracts” is incomplete in 3D. The axis/orientation geometry matters.

## Finite-mouth correction in the coaxial far-field case

For coaxial loops, Step 3 gave
\[
F_d
= -\frac{3\mu_0\pi I_1I_2R_1^2R_2^2}{2d^4}
\left[1-\frac{5}{2}\frac{R_1^2+R_2^2}{d^2}+\cdots\right].
\]
The correction changes the magnitude and eventually the far-field expansion breaks down near contact. The audit also computes the next asymptotic term in the ratio,
\[
\frac{35}{8}\frac{R_1^4+3R_1^2R_2^2+R_2^4}{d^4}.
\]
This is a far-field asymptotic expansion, not a proof of convergence near contact. In the controlled regime
\[
d^2\gg R_1^2+R_2^2,
\]
the correction does not change the leading sign.

Self-energy terms depending only on \(R_1\) or \(R_2\) do not produce a pair force because they have no \(d\)-dependence.

## Step verdict

The 3D sign rule is geometric:
\[
\boxed{F_d\propto -(3ab-c).}
\]
For the coaxial facing-mouth branch, this reduces to the sign table in Step 4.

## SymPy audit

The script builds vector dipoles for the main 3D orientation cases, recovers the Step-3 coaxial leading force, and prints the finite-size asymptotic correction ratio.

Run:

```bash
python step_05_3d_orientation_and_finite_size_sympy.py
```
