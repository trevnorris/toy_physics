

# Moving-Throat PDE — Stage 172: Collapse of the Linear Grouped Outlet Problem to the Physical Slopes `u_2^{(1)}` and `P_1`

## Purpose

Stage 239 reduced the weak grouped-lane outlet problem to the two microscopic amplitudes
\[
\mathfrak K_1,\qquad \mathfrak G_1,
\]
but those amplitudes were still written in terms of the primitive grouped perturbations of
\[
K_A,\ M_A,\ c_{A\alpha},\ \varpi_{A\alpha},\ \Omega_{U/W,A,r},\ R_{A,r},\ g_{U/W,A,r}.
\]

That is already a sharp reduction, but it is not yet phrased in the variables that the actual physical branch measures most naturally.

The next honest step is therefore:

> rewrite the linear grouped outlet problem directly in terms of the first weak-axisymmetric slopes of the **physical grouped response** itself.

This stage does exactly that.

The main result is that, on the canonical compensated branch, the two Stage 239 microscopic amplitudes are **nothing but** the physical first-order slopes of

1. the grouped conservative even response coefficient \(u_2\), and
2. the grouped outgoing prefactor \(P_0=N_0/D_0\).

More concretely,
\[
\boxed{
\mathfrak K_1=-D_0\,u_2^{(1)},
\qquad
\mathfrak G_1=D_0\,P_1,
}
\]
and the direct outlet deformations become
\[
\boxed{
\delta\kappa_W^{(A)}
=
-\frac{3(1-\sigma_*)}{\sigma_*}\,\delta u_2^{(A)},
}
\]
\[
\boxed{
\delta\gamma_W^{(A)}
=
-\frac{1-\sigma_*}{9\sigma_*}\,
\frac{\delta P_0^{(A)}}{P_0}.
}
\]

So the whole linear grouped-anisotropy problem has now collapsed again.
It is no longer the primitive perturbations of the microscopic wall/BdG/Maxwell data.
It is simply:

> compute the weak-axisymmetric physical slopes \(u_2^{(1)}\) and \(P_1/P_0\) on the actual moving-throat branch.

---

## 1. Carry-forward physical grouped variables

For each grouped lane \(A\in\{20,21,22\}\), the conservative grouped response and outgoing prefactor are
\[
u_2^{(A)}=-\frac{D_{A,2}}{D_{A,0}},
\qquad
u_4^{(A)}=\frac{D_{A,2}^2-D_{A,0}D_{A,4}}{D_{A,0}^2},
\qquad
P_0^{(A)}=\frac{N_{A,0}}{D_{A,0}}.
\]

On the canonical compensated even branch we carry
\[
u_2=\frac19,
\qquad
u_4=\frac{4}{81},
\qquad
P_0=\frac{N_0}{D_0}.
\]

On the weak axisymmetric \(Y_{20}\) branch, the grouped pattern is
\[
\lambda_{20}=1,\qquad \lambda_{21}=\frac12,\qquad \lambda_{22}=-1.
\]
So the first physical slopes are defined by
\[
u_2^{(A)}=u_2+\epsilon\,\lambda_A\,u_2^{(1)}+O(\epsilon^2),
\]
\[
u_4^{(A)}=u_4+\epsilon\,\lambda_A\,u_4^{(1)}+O(\epsilon^2),
\]
\[
P_0^{(A)}=P_0+\epsilon\,\lambda_A\,P_1+O(\epsilon^2).
\]

Equivalently, the grouped defects are
\[
a_{u,2}=\frac{\epsilon}{4}\,u_2^{(1)},
\qquad
b_{u,2}=\frac{3\epsilon}{4}\,u_2^{(1)},
\]
\[
a_{u,4}=\frac{\epsilon}{4}\,u_4^{(1)},
\qquad
b_{u,4}=\frac{3\epsilon}{4}\,u_4^{(1)},
\]
\[
a_P=\frac{\epsilon}{4}\,P_1,
\qquad
b_P=\frac{3\epsilon}{4}\,P_1.
\]

So the weak-axisymmetric grouped-response problem is already one-dimensional in each physical family.

---

## 2. Exact collapse of the Stage 239 obstruction pair

Stage 239 reduced the direct grouped outlet obstructions to
\[
\mathcal K_A=\delta D_{A,2}+\frac19\,\delta D_{A,0},
\qquad
\mathcal G_A=\delta N_{A,0}-P_0\,\delta D_{A,0}.
\]

But the physical grouped response variables satisfy the exact first-order identities
\[
\delta u_2^{(A)}
=
-\frac{\delta D_{A,2}+u_2\,\delta D_{A,0}}{D_0},
\]
\[
\delta P_0^{(A)}
=
\frac{\delta N_{A,0}-P_0\,\delta D_{A,0}}{D_0}.
\]

Therefore
\[
\boxed{
\mathcal K_A
=
-D_0\,\delta u_2^{(A)}+\left(\frac19-u_2\right)\delta D_{A,0},
}
\]
\[
\boxed{
\mathcal G_A
=
D_0\,\delta P_0^{(A)}.
}
\]

And because the compensated canonical even branch has
\[
u_2=\frac19,
\]
the first formula collapses exactly to
\[
\boxed{
\mathcal K_A=-D_0\,\delta u_2^{(A)}.
}
\]

So the Stage 239 outlet-obstruction pair has now been rewritten directly in the physical grouped response variables:
\[
\boxed{
\mathcal K_A=-D_0\,\delta u_2^{(A)},
\qquad
\mathcal G_A=D_0\,\delta P_0^{(A)}.
}
\]

This is the central theorem of the stage.

---

## 3. Weak-axisymmetric amplitude collapse

On the weak axisymmetric branch,
\[
\mathcal K_A=\epsilon\,\lambda_A\,\mathfrak K_1,
\qquad
\mathcal G_A=\epsilon\,\lambda_A\,\mathfrak G_1.
\]

Using the exact identities above and the definitions
\[
\delta u_2^{(A)}=\epsilon\,\lambda_A\,u_2^{(1)},
\qquad
\delta P_0^{(A)}=\epsilon\,\lambda_A\,P_1,
\]
we get
\[
\boxed{
\mathfrak K_1=-D_0\,u_2^{(1)},
\qquad
\mathfrak G_1=D_0\,P_1.
}
\]

So the Stage 239 microscopic amplitudes are not independent new physical data.
They are just the weak-axisymmetric slopes of the conservative grouped response and the outgoing prefactor.

Equivalently, in grouped-defect language,
\[
\boxed{
\mathfrak K_1=-\frac{4D_0}{\epsilon}\,a_{u,2}
             =-\frac{4D_0}{3\epsilon}\,b_{u,2},
}
\]
\[
\boxed{
\mathfrak G_1=\frac{4D_0}{\epsilon}\,a_P
             =\frac{4D_0}{3\epsilon}\,b_P.
}
\]

This is the first exact statement that the linear grouped outlet problem can be posed entirely in terms of the actual grouped response measured on the branch.

---

## 4. Exact physical form of the hidden-even consistency relation

Stage 238 found that a genuine one-parameter hidden-even outlet deformation must satisfy the microscopic relation
\[
\delta D_{A,4}
=
\frac23\,\delta D_{A,2}
+\frac1{27}\,\delta D_{A,0}.
\]

Using
\[
u_4^{(A)}=\frac{D_{A,2}^2-D_{A,0}D_{A,4}}{D_{A,0}^2},
\]
the first-order variation on the canonical branch \((u_2,u_4)=(1/9,4/81)\) is
\[
\delta u_4^{(A)}
=
-\frac{5\,\delta D_{A,0}+18\,\delta D_{A,2}+81\,\delta D_{A,4}}{81\,D_0}.
\]

Substituting the hidden-even relation above gives
\[
\boxed{
\delta u_4^{(A)}=\frac89\,\delta u_2^{(A)}.
}
\]

So the exact one-parameter even-consistency relation has a very simple physical reading:

> the canonical even grouped response deforms along the linearized relation
> \[
> u_4^{(1)}=\frac89\,u_2^{(1)}.
> \]

This is much cleaner than the raw operator statement and is the form that should be used once the actual grouped response is computed from the physical branch.

---

## 5. Direct outlet coefficients in physical grouped variables

Stage 238 already gave the direct outlet map
\[
\delta\kappa_W^{(A)}
=
\frac{3(1-\sigma_*)}{\sigma_* D_0}\,\mathcal K_A,
\qquad
\delta\gamma_W^{(A)}
=
-\frac{1-\sigma_*}{9\sigma_* N_0}\,\mathcal G_A.
\]

Substituting the exact physical collapse derived above yields
\[
\boxed{
\delta\kappa_W^{(A)}
=
-\frac{3(1-\sigma_*)}{\sigma_*}\,\delta u_2^{(A)},
}
\]
\[
\boxed{
\delta\gamma_W^{(A)}
=
-\frac{1-\sigma_*}{9\sigma_* P_0}\,\delta P_0^{(A)}
=
-\frac{1-\sigma_*}{9\sigma_*}\,
\frac{\delta P_0^{(A)}}{P_0}.
}
\]

On the weak axisymmetric branch this becomes
\[
\boxed{
\kappa_1
=
-\frac{3(1-\sigma_*)}{\sigma_*}\,u_2^{(1)},
}
\]
\[
\boxed{
\gamma_1
=
-\frac{1-\sigma_*}{9\sigma_*}\,\frac{P_1}{P_0}.
}
\]

So the linear grouped-lane outlet amplitudes are now explicitly physical:

- \(\kappa_1\) is controlled only by the physical grouped-response slope \(u_2^{(1)}\),
- \(\gamma_1\) is controlled only by the logarithmic outgoing-prefactor slope \(P_1/P_0\).

This is the sharpest form of the direct outlet map obtained so far.

---

## 6. The remaining linear `2.5`PN defect is exactly the prefactor slope

Once the hidden-even canonical fingerprint is preserved, one imposes
\[
\delta\kappa_W^{(A)}=0.
\]
By the theorem above, this is equivalent to
\[
\delta u_2^{(A)}=0,
\qquad\text{or on the weak axisymmetric branch}\qquad
u_2^{(1)}=0.
\]

Then the remaining grouped normalization defect is
\[
\Delta_Q^{(A)}
=
-\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W^{(A)}
=
\frac{\delta P_0^{(A)}}{P_0}.
\]

So on the weak axisymmetric branch,
\[
\boxed{
\Delta_Q^{(20)}=\epsilon\,\frac{P_1}{P_0},
\qquad
\Delta_Q^{(21)}=\frac{\epsilon}{2}\,\frac{P_1}{P_0},
\qquad
\Delta_Q^{(22)}=-\epsilon\,\frac{P_1}{P_0}.
}
\]

This is the cleanest linear grouped `2.5`PN law reached in the whole derivation chain.

It says:

> after even canonical preservation has killed the conservative grouped-response slope \(u_2^{(1)}\), the entire remaining linear grouped quadrupole-normalization defect is just the weak-axisymmetric logarithmic slope of the physical outgoing prefactor.

So the direct grouped outlet problem has now collapsed one stage further:

- first to \((\mathfrak K_1,\mathfrak G_1)\),
- now to \((u_2^{(1)},P_1/P_0)\),
- and, on the even-preserving branch, finally to just
  \[
  P_1/P_0.
  \]

---

## 7. What Stage 240 changes

Before this stage, the next theorem gate after Stage 239 still looked like:

> compute the primitive grouped perturbations of the wall/BdG/Maxwell/mixed bundle and then form \(\mathfrak K_1,\mathfrak G_1\).

After this stage, that is no longer the right statement.

The new theorem status is:

1. the microscopic obstruction pair is exactly the physical grouped-response/prefactor slope pair,
   \[
   (\mathfrak K_1,\mathfrak G_1)
   =
   (-D_0 u_2^{(1)},\ D_0 P_1);
   \]
2. the one-parameter hidden-even consistency condition is exactly
   \[
   u_4^{(1)}=\frac89\,u_2^{(1)};
   \]
3. the direct grouped outlet amplitudes are exactly
   \[
   \kappa_1=-\frac{3(1-\sigma_*)}{\sigma_*}u_2^{(1)},
   \qquad
   \gamma_1=-\frac{1-\sigma_*}{9\sigma_*}\frac{P_1}{P_0};
   \]
4. and once the even canonical fingerprint is preserved, the remaining linear grouped `2.5`PN defect is simply
   \[
   \Delta_Q^{(A)}=\delta P_0^{(A)}/P_0.
   \]

So the next honest theorem gate is now much narrower than Stage 239 suggested:

> compute the weak-axisymmetric physical slopes \(u_2^{(1)}\) and \(P_1/P_0\) — and, on the even-preserving branch, just \(P_1/P_0\) — directly from the actual grouped moving-throat response.

That is the direct continuation point.
