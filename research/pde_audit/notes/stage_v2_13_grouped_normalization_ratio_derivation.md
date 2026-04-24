# Stage V2-13 — Grouped Normalization Ratio Audit

## Purpose

This stage freezes the exact algebraic bridge from grouped real `P2` conservative operator moments and outgoing-transfer moments to the normalized quantities used by the 2.5PN/4PN quadrupole-normalization package.

The stage treats one grouped lane first, then specializes to the isotropic full-bundle branch. It does **not** solve the moving-throat PDE. It identifies the exact target surface that the actual branch must hit.

The retained low-frequency inputs are

\[
D(\omega)=D_0+D_2\omega^2+D_4\omega^4+O(\omega^6),
\]

and

\[
N(\omega)=N_0+N_2\omega^2+N_4\omega^4+O(\omega^6).
\]

Here `D` is the conservative wall/worldtube operator and `N` is the outgoing-transfer numerator inherited from the Maxwell/mixed port.

A port-normalization factor \(\mathcal S_{\rm port}\) is carried explicitly. If `N` has already been gravitationally/port-normalized, set

\[
\mathcal S_{\rm port}=1.
\]

If `N` is the raw canonical mechanical transfer numerator, keep \(\mathcal S_{\rm port}\) until the comparison with the GR target.

---

## 1. Normalized conservative response

Define the normalized conservative response

\[
Y(\omega)=\frac{D_0}{D(\omega)}.
\]

Expanding,

\[
Y(\omega)=1+u_2\omega^2+u_4\omega^4+O(\omega^6).
\]

The exact conversion formulas are

\[
\boxed{u_2=-\frac{D_2}{D_0}},
\]

\[
\boxed{u_4=\frac{D_2^2-D_0D_4}{D_0^2}}.
\]

These formulas are purely algebraic and are independent of the microscopic origin of the coefficients.

---

## 2. Outgoing prefactor moments

The outgoing transfer appears through the internal prefactor

\[
{\rm Pref}(\omega)=\frac{D_0N(\omega)}{D(\omega)^2}.
\]

Write

\[
{\rm Pref}(\omega)=P_0+P_2\omega^2+P_4\omega^4+O(\omega^6).
\]

Then

\[
\boxed{P_0=\frac{N_0}{D_0}},
\]

\[
\boxed{P_2=\frac{D_0N_2-2D_2N_0}{D_0^2}},
\]

\[
\boxed{
P_4=
\frac{
D_0^2N_4
-2D_0(D_2N_2+D_4N_0)
+3D_2^2N_0
}{D_0^3}.
}
\]

Thus the static outgoing prefactor is the simple ratio

\[
\boxed{P_0=N_0/D_0}.
\]

This is the scalar that controls the leading universal quadrupole channel.

---

## 3. Constant-prefactor branch

The constant-prefactor branch through \(O(\omega^4)\) is defined by

\[
P_2=0,
\qquad
P_4=0.
\]

The exact branch conditions are

\[
\boxed{N_2=\frac{2D_2N_0}{D_0}},
\]

and

\[
\boxed{N_4=\frac{N_0(D_2^2+2D_0D_4)}{D_0^2}}.
\]

Equivalently, if the first condition is kept symbolic, the second may be written as

\[
N_4=
\frac{2D_0(D_2N_2+D_4N_0)-3D_2^2N_0}{D_0^2}.
\]

The important point is that the higher transfer moments do not need to vanish. They must be correlated with the conservative moments.

---

## 4. Multiplication by the compact outgoing \(l=2\) fingerprint

The compact outgoing quadrupole branch has normalized low-frequency form

\[
\widehat Y_2^{\rm out}(\omega)
=1+A\omega^2+B\omega^4+iG_5\omega^5+O(\omega^6),
\]

where, on the canonical compact branch,

\[
A=\frac{a^2}{9c_s^2},
\qquad
B=\frac{4a^4}{81c_s^4},
\qquad
G_5=\frac{a^5}{27c_s^5}.
\]

Multiplying by the internal prefactor gives

\[
\delta Y^{\rm out}(\omega)
={\rm Pref}(\omega)\widehat Y_2^{\rm out}(\omega).
\]

Thus

\[
K_0=P_0,
\]

\[
K_2=P_2+AP_0,
\]

\[
K_4=P_4+AP_2+BP_0,
\]

and

\[
\boxed{\Gamma_5=G_5P_0}.
\]

Therefore the leading odd 2.5PN coefficient depends only on the static prefactor \(P_0\), not on \(P_2\) or \(P_4\).

---

## 5. Isotropic full-bundle surface

On the isotropic grouped branch, the full coupled conservative operator is

\[
D_0=K-B_0-Z_0,
\]

\[
D_2=-(M+B_2+Z_2),
\]

\[
D_4=-(B_4+Z_4).
\]

Here:

- \(K,M\) are wall stiffness and inertia data,
- \(B_n\) are stable BdG support moments,
- \(Z_n\) are conservative Maxwell/mixed moments.

The one-pole conservative condition

\[
u_4=4u_2^2
\]

is equivalent to

\[
\boxed{
D_0(B_4+Z_4)=3(M+B_2+Z_2)^2.
}
\]

So the one-pole surface fixes the static denominator as

\[
D_0^{\rm one\ pole}
=
\frac{3(M+B_2+Z_2)^2}{B_4+Z_4}.
\]

The universal quadrupole target is

\[
\boxed{
\widehat m_0^{\,2}\,\mathcal S_{\rm port}\,P_0
=
\frac{54Gc_s^5}{5a^5c^5}.
}
\]

Since \(P_0=N_0/D_0\), this gives

\[
D_0^{\rm norm}
=
\frac{
\widehat m_0^{\,2}\mathcal S_{\rm port}N_0
}{T_{\rm GR}},
\]

with

\[
T_{\rm GR}=\frac{54Gc_s^5}{5a^5c^5}.
\]

The simultaneous one-pole and normalization compatibility surface is therefore

\[
\boxed{
\frac{
\widehat m_0^{\,2}\mathcal S_{\rm port}N_0
}{T_{\rm GR}}
=
\frac{3(M+B_2+Z_2)^2}{B_4+Z_4}.
}
\]

This is the exact isotropic scalar target surface for the reduced moving-throat branch.

---

## 6. GR coefficient check

Using

\[
\widehat m_0^{\,2}\mathcal S_{\rm port}P_0
=
\frac{54Gc_s^5}{5a^5c^5},
\]

and

\[
G_5=\frac{a^5}{27c_s^5},
\]

one obtains

\[
\widehat m_0^{\,2}\mathcal S_{\rm port}\Gamma_5
=
\frac{54Gc_s^5}{5a^5c^5}\frac{a^5}{27c_s^5}
=
\boxed{\frac{2G}{5c^5}}.
\]

So the normalization target is exactly equivalent to the standard quadrupole coefficient.

---

## 7. Weak-axisymmetric grouped transport

For a weak axisymmetric grouped perturbation with signature

\[
(20,21,22)\sim(1,\tfrac12,-1),
\]

let

\[
P_A=P_{\rm base}+\epsilon\lambda_A P_1.
\]

Then the grouped trace and anisotropy coordinates are

\[
\bar P=P_{\rm base},
\]

\[
a_P=\frac{\epsilon P_1}{4},
\]

\[
b_P=\frac{3\epsilon P_1}{4}.
\]

Therefore

\[
\boxed{b_P=3a_P}.
\]

This confirms that the weak-axisymmetric prefactor slope remains a one-scalar defect, not a generic three-lane deformation.

---

## 8. Stage result

The audit verifies 15 symbolic checks:

- normalized response formulas,
- outgoing prefactor formulas,
- constant-prefactor branch conditions,
- outgoing \(l=2\) multiplication,
- isotropic one-pole target surface,
- universal normalization to GR quadrupole coefficient,
- weak-axisymmetric grouped trace and \(b=3a\) law.

Final status:

\[
\boxed{\text{V2-13 passes as exact algebra.}}
\]

The remaining issue is not algebraic. It is branch realization: the actual moving-throat PDE must supply stable isotropic data satisfying the compatibility surface above.
