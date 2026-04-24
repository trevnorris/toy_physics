# Stage V2-19 — Isotropic Full-Bundle Target Surface

## Purpose

This stage consolidates the isotropic moving-throat full-bundle target into a finite pass/fail packet for an actual branch extraction.

Earlier Volume 2 stages verified the ingredients separately:

- the conservative grouped response bridge,
- the outgoing `l=2` fingerprint,
- the 2.5PN / 4PN normalization interface,
- the branch-freeze protocol,
- the weak-axisymmetric `b=3a` line,
- and the quotient / similarity-orbit structure.

V2-19 turns those into one scalar target surface for the isotropic branch.

The stage is deliberately algebraic. It does **not** solve the nonlinear moving-throat PDE. It states what the actual PDE branch must output after the branch data have been frozen.

---

## 1. Isotropic full-bundle definitions

On the isotropic grouped real `P2` branch, all three grouped lanes share the same conservative operator moments and outgoing-transfer moments.

Define

\[
D(\omega)=D_0+D_2\omega^2+D_4\omega^4+O(\omega^6),
\]

with

\[
D_0=K-B_0-Z_0,
\]

\[
D_2=-(M+B_2+Z_2),
\]

\[
D_4=-(B_4+Z_4).
\]

It is useful to abbreviate

\[
A\equiv M+B_2+Z_2,
\]

\[
C\equiv B_4+Z_4.
\]

Then

\[
D_2=-A,
\qquad
D_4=-C.
\]

Here:

- `K,M` are the wall/worldtube static stiffness and inertia data;
- `B_0,B_2,B_4` are the stable BdG support moments;
- `Z_0,Z_2,Z_4` are the conservative Maxwell/mixed moments;
- `N_0,N_2,N_4` are the outgoing-transfer moments.

The outgoing prefactor is defined by

\[
P(\omega)=\frac{D_0N(\omega)}{D(\omega)^2},
\]

where

\[
N(\omega)=N_0+N_2\omega^2+N_4\omega^4+O(\omega^6).
\]

---

## 2. Conservative one-pole surface

The normalized conservative response is

\[
Y(\omega)=\frac{D_0}{D(\omega)}.
\]

Expanding gives

\[
Y(\omega)=1+u_2\omega^2+u_4\omega^4+O(\omega^6),
\]

with

\[
u_2=\frac{A}{D_0},
\]

\[
u_4=\frac{A^2+D_0C}{D_0^2}.
\]

The one-pole condition is

\[
u_4=4u_2^2.
\]

Substitution gives

\[
\frac{A^2+D_0C}{D_0^2}=4\frac{A^2}{D_0^2},
\]

so the exact one-pole surface is

\[
\boxed{
D_0C=3A^2.
}
\]

Equivalently,

\[
\boxed{
D_0=\frac{3A^2}{C}.
}
\]

Since

\[
D_0=K-B_0-Z_0,
\]

the required wall stiffness on the one-pole surface is

\[
\boxed{
K=B_0+Z_0+\frac{3(M+B_2+Z_2)^2}{B_4+Z_4}.
}
\]

This is a branch-output condition. It is not a license to refit `K` after target evaluation.

---

## 3. Universal quadrupole normalization

The universal 2.5PN / 4PN quadrupole target is

\[
\widehat m_0^2\mathcal S_{\rm port}P_0
=
\frac{54Gc_s^5}{5a^5c^5}.
\]

Define

\[
T_{\rm GR}
\equiv
\frac{54Gc_s^5}{5a^5c^5}.
\]

Since

\[
P_0=\frac{N_0}{D_0},
\]

the normalization gate is

\[
\boxed{
\widehat m_0^2\mathcal S_{\rm port}\frac{N_0}{D_0}=T_{\rm GR}.
}
\]

Equivalently, as a polynomial residual,

\[
\boxed{
R_{\rm norm}
=
\widehat m_0^2\mathcal S_{\rm port}N_0
-T_{\rm GR}D_0=0.
}
\]

Solving for `N_0` gives

\[
\boxed{
N_0=\frac{T_{\rm GR}D_0}{\widehat m_0^2\mathcal S_{\rm port}}.
}
\]

On the one-pole surface this becomes

\[
\boxed{
N_0
=
\frac{3T_{\rm GR}A^2}{\widehat m_0^2\mathcal S_{\rm port}C}.
}
\]

Writing out `T_GR`,

\[
\boxed{
N_0
=
\frac{162Gc_s^5(M+B_2+Z_2)^2}
{5\mathcal S_{\rm port}a^5c^5\widehat m_0^2(B_4+Z_4)}.
}
\]

The SymPy audit also verifies that this is exactly equivalent to

\[
\gamma_{\rm quad}^{\rm eff}
=
\frac{2G}{5c^5},
\]

where

\[
\gamma_{\rm quad}^{\rm eff}
=
\widehat m_0^2\mathcal S_{\rm port}P_0\frac{a^5}{27c_s^5}.
\]

---

## 4. Constant-prefactor outgoing branch

The prefactor expansion is

\[
P(\omega)=P_0+P_2\omega^2+P_4\omega^4+O(\omega^6).
\]

The exact coefficients are

\[
P_0=\frac{N_0}{D_0},
\]

\[
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\]

\[
P_4=
\frac{
D_0^2N_4
-2D_0(D_2N_2+D_4N_0)
+3D_2^2N_0
}{D_0^3}.
\]

The constant-prefactor branch requires

\[
P_2=0,
\qquad
P_4=0.
\]

The exact moment constraints are therefore

\[
\boxed{
N_2=\frac{2D_2N_0}{D_0}.
}
\]

Since `D_2=-A`,

\[
\boxed{
N_2=-\frac{2AN_0}{D_0}.
}
\]

The next condition is

\[
\boxed{
N_4=\frac{N_0(D_2^2+2D_0D_4)}{D_0^2}.
}
\]

Since `D_2=-A` and `D_4=-C`,

\[
\boxed{
N_4=\frac{N_0(A^2-2D_0C)}{D_0^2}.
}
\]

On the one-pole surface `D_0C=3A^2`, this simplifies to

\[
\boxed{
N_4=-\frac{5A^2N_0}{D_0^2}.
}
\]

Using the normalization surface as well, the two outgoing-transfer target moments become

\[
\boxed{
N_2
=
-\frac{2AT_{\rm GR}}{\widehat m_0^2\mathcal S_{\rm port}},
}
\]

and

\[
\boxed{
N_4
=
-\frac{5CT_{\rm GR}}{3\widehat m_0^2\mathcal S_{\rm port}}.
}
\]

Equivalently, after writing out `T_GR`,

\[
\boxed{
N_2
=
-\frac{108Gc_s^5(M+B_2+Z_2)}
{5\mathcal S_{\rm port}a^5c^5\widehat m_0^2},
}
\]

\[
\boxed{
N_4
=
-\frac{18Gc_s^5(B_4+Z_4)}
{\mathcal S_{\rm port}a^5c^5\widehat m_0^2}.
}
\]

---

## 5. Optional tail-transport gate

If the 4PN tail transport is not already derived, carry the scalar gate

\[
\boxed{
R_{\rm tail}
=
\Theta_{\rm tail}\left(\frac{c}{c_s}\right)^3-1=0.
}
\]

Thus

\[
\boxed{
\Theta_{\rm tail}=\left(\frac{c_s}{c}\right)^3.
}
\]

On the branch `c_s=c`, this reduces to

\[
\boxed{
\Theta_{\rm tail}=1.
}
\]

This is separate from the grouped quadrupole normalization. It is not a new `P2` normalization slot.

---

## 6. Final V2-19 target packet

Let

\[
A=M+B_2+Z_2,
\qquad
C=B_4+Z_4.
\]

The isotropic full-bundle branch passes V2-19 only if the frozen branch output satisfies:

\[
\boxed{
D_0=K-B_0-Z_0=\frac{3A^2}{C},
}
\]

\[
\boxed{
\widehat m_0^2\mathcal S_{\rm port}\frac{N_0}{D_0}
=
\frac{54Gc_s^5}{5a^5c^5},
}
\]

\[
\boxed{
N_2=\frac{2D_2N_0}{D_0},
}
\]

\[
\boxed{
N_4=\frac{N_0(D_2^2+2D_0D_4)}{D_0^2},
}
\]

and, if the tail transport is not otherwise closed,

\[
\boxed{
\Theta_{\rm tail}\left(\frac{c}{c_s}\right)^3=1.
}
\]

---

## 7. Feasibility and stability implications

The one-pole surface gives

\[
D_0=\frac{3A^2}{C}.
\]

For a stable branch we require

\[
D_0>0.
\]

Therefore, for a nondegenerate branch with `A != 0`, we need

\[
\boxed{
C=B_4+Z_4>0.
}
\]

If `C <= 0`, the one-pole surface is incompatible with `D_0>0`, except for degenerate cases that would also collapse the target packet.

The normalized outgoing transfer is

\[
N_0=\frac{T_{\rm GR}D_0}{\widehat m_0^2\mathcal S_{\rm port}}.
\]

So if

\[
G>0,
\quad c_s>0,
\quad a>0,
\quad c>0,
\quad \widehat m_0^2>0,
\quad \mathcal S_{\rm port}>0,
\quad D_0>0,
\]

then the target branch has

\[
N_0>0.
\]

Thus it is not a dark port.

---

## 8. Local codimension check

The SymPy audit uses the residual vector

\[
R=
\begin{pmatrix}
R_{\rm pole}\\
R_{\rm norm}\\
R_{P2}\\
R_{P4}\\
R_{\rm tail}
\end{pmatrix}
\]

with output variables

\[
(K,N_0,N_2,N_4,\Theta_{\rm tail}).
\]

The Jacobian determinant is

\[
\boxed{
\det J
=
-\mathcal S_{\rm port}\,m_0^2
\left(\frac{c}{c_s}\right)^3
(B_4+Z_4)
(B_0-K+Z_0)^3.
}
\]

Since

\[
B_0-K+Z_0=-D_0,
\]

this determinant is nonzero whenever

\[
\mathcal S_{\rm port}\ne0,
\quad
m_0\ne0,
\quad
c_s\ne0,
\quad
B_4+Z_4\ne0,
\quad
D_0\ne0.
\]

On the target surface, the determinant becomes

\[
\boxed{
\det J\big|_{\rm surf}
=
\frac{
27\mathcal S_{\rm port}m_0^2c^3A^6
}{
c_s^3C^2
}.
}
\]

So the target packet is locally codimension five in those algebraic output slots. This is exactly why the branch-freeze rule matters: without freezing, the same slots could be tuned to hit the residuals.

---

## 9. SymPy audit result

The script performed 18 symbolic checks. All passed.

```text
checks_total: 18
checks_passed: 18
checks_failed: 0
```

The most important output formulas are:

```text
D0_surface = 3*(B2 + M + Z2)^2/(B4 + Z4)

K_surface = B0 + Z0 + 3*(B2 + M + Z2)^2/(B4 + Z4)

P0_surface = 54*G*c_s^5/(5*S_port*a_th^5*c^5*mhat0^2)

N0_surface = 162*G*c_s^5*(B2 + M + Z2)^2/
             (5*S_port*a_th^5*c^5*mhat0^2*(B4 + Z4))

N2_surface = -108*G*c_s^5*(B2 + M + Z2)/
             (5*S_port*a_th^5*c^5*mhat0^2)

N4_surface = -18*G*c_s^5*(B4 + Z4)/
             (S_port*a_th^5*c^5*mhat0^2)
```

---

## 10. Carry-forward statement

V2-19 gives the isotropic target sheet for the actual moving-throat branch:

> A frozen isotropic branch must output stable `D0`, non-dark `N0`, one-pole conservative moments, constant-prefactor outgoing moments, and the universal quadrupole normalization using one common `G`. If the tail scalar has not already been derived, the same branch must also satisfy the separate transport gate `Theta_tail(c/c_s)^3=1`.

This is the exact target packet that V2-20 should feed into a weak-form / numerical extraction plan.
