# Moving-Throat PDE — Stage 153: Exact Linear Grouped-`P2` Map into the Direct Outlet Coefficients `\delta\kappa_W` and `\delta\gamma_W`

## Purpose

Stage 152 proved that a pure weak grouped real `P2` anisotropy cannot linearly feed the scalar
off-bundle slippages
\[
\varepsilon_L,\qquad \varepsilon_v,\qquad \varepsilon_T,
\]
and therefore cannot linearly source the scalar normal defect `\delta_\perp`.

So the remaining **linear** grouped-anisotropy problem is now exactly what Stage 152 said it had to be:

> the direct outlet coefficients
> \[
> \delta\kappa_W,\qquad \delta\gamma_W.
> \]

This stage solves that map.

The main result is that the full linear grouped-`P2` outlet problem collapses to **two exact microscopic combinations** of the grouped bundle data:

1. one even combination
   \[
   \mathcal K_A:=\delta D_{A,2}+\frac19\,\delta D_{A,0},
   \]
   which feeds the hidden even pole deformation `\delta\kappa_W`;

2. one odd combination
   \[
   \mathcal G_A:=\delta N_{A,0}-P_0\,\delta D_{A,0},
   \]
   which feeds the hidden odd normalization deformation `\delta\gamma_W`.

Equivalently,
\[
\boxed{
\delta\kappa_W^{(A)}
=
\frac{3(1-\sigma_*)}{\sigma_* D_0}\,\mathcal K_A,
}
\qquad
\boxed{
\delta\gamma_W^{(A)}
=
-\frac{1-\sigma_*}{9\sigma_* N_0}\,\mathcal G_A.
}
\]

So the remaining linear grouped-anisotropy bottleneck is not the whole bundle
\[
(\delta D_{A,0},\delta D_{A,2},\delta D_{A,4},\delta N_{A,0}).
\]
It is just the pair
\[
(\mathcal K_A,\mathcal G_A).
\]

That is a substantial narrowing of the theorem gate.
It is a linear-response statement on the compensated isotropic branch; it does not
yet say which grouped anisotropies the full PDE dynamically realizes.

---

## 1. Carry-forward grouped-lane variables

For each grouped real quadrupole lane
\[
A\in\{20,21,22\},
\]
Stage 5/6 already gave the normalized conservative response moments
\[
u_2^{(A)}=-\frac{D_{A,2}}{D_{A,0}},
\qquad
u_4^{(A)}=\frac{D_{A,2}^2-D_{A,0}D_{A,4}}{D_{A,0}^2},
\]
and the static outgoing prefactor
\[
P_0^{(A)}=\frac{N_{A,0}}{D_{A,0}}.
\]

Now linearize around an isotropic compensated canonical branch:
\[
D_{A,n}=D_n+\delta D_{A,n},
\qquad
N_{A,0}=N_0+\delta N_{A,0}.
\]
On that branch the canonical values are
\[
u_2=\frac19,
\qquad
u_4=\frac{4}{81},
\qquad
P_0=\frac{N_0}{D_0}.
\]

Stage 6 already implies the first-order grouped formulas
\[
\boxed{
\delta u_2^{(A)}
=
-\frac{\delta D_{A,2}+\frac19\,\delta D_{A,0}}{D_0},
}
\]
\[
\boxed{
\delta P_0^{(A)}
=
\frac{\delta N_{A,0}-P_0\,\delta D_{A,0}}{D_0}.
}
\]

For later use, the exact first-order `u_4` transport is
\[
\boxed{
\delta u_4^{(A)}
=
-\frac{\delta D_{A,4}+\frac29\,\delta D_{A,2}+\frac{5}{81}\,\delta D_{A,0}}{D_0}.
}
\]

---

## 2. Direct hybrid-outlet defects on the pure grouped-anisotropy branch

Stage 142 gave the exact compensated-hybrid outlet formulas
\[
\delta E_2
=
\frac{\delta\mathcal C-9\sigma_*\,\delta\kappa_W}{27(1-\sigma_*)},
\]
\[
\delta E_4
=
\frac{5\delta\mathcal C-72\sigma_*\,\delta\kappa_W}{243(1-\sigma_*)},
\]
\[
\Delta_Q
=
\frac{\delta\mathcal C-27\sigma_*\,\delta\gamma_W}{3(1-\sigma_*)}.
\]

Stage 152 then proved that for a **pure linear grouped real `P2` anisotropy**
\[
\delta\mathcal C^{(1,P_2)}=0.
\]

So lane by lane the direct outlet defects reduce to
\[
\boxed{
\delta E_2^{(A)}
=
-\frac{\sigma_*}{3(1-\sigma_*)}\,\delta\kappa_W^{(A)},
}
\]
\[
\boxed{
\delta E_4^{(A)}
=
-\frac{8\sigma_*}{27(1-\sigma_*)}\,\delta\kappa_W^{(A)},
}
\]
\[
\boxed{
\Delta_Q^{(A)}
=
-\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W^{(A)}.
}
\]

Now identify the grouped conservative and outgoing variables:

- the even defect is just the grouped response defect
  \[
  \delta E_2^{(A)}=\delta u_2^{(A)};
  \]
- the odd normalization defect is the relative outgoing-prefactor defect
  \[
  \chi_Q^{(A)}=\frac{P_0^{(A)}}{P_0},
  \qquad
  \Delta_Q^{(A)}=\chi_Q^{(A)}-1=\frac{\delta P_0^{(A)}}{P_0}.
  \]

Therefore
\[
\boxed{
\delta u_2^{(A)}
=
-\frac{\sigma_*}{3(1-\sigma_*)}\,\delta\kappa_W^{(A)},
}
\qquad
\boxed{
\frac{\delta P_0^{(A)}}{P_0}
=
-\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W^{(A)}.
}
\]

Solving gives the exact linear map
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
-\frac{1-\sigma_*}{9\sigma_*}\,\frac{\delta P_0^{(A)}}{P_0}.
}
\]

Substituting the Stage-6 grouped-bundle transport laws then yields
\[
\boxed{
\delta\kappa_W^{(A)}
=
\frac{3(1-\sigma_*)}{\sigma_* D_0}
\left(
\delta D_{A,2}
+\frac19\,\delta D_{A,0}
\right),
}
\]
\[
\boxed{
\delta\gamma_W^{(A)}
=
-\frac{1-\sigma_*}{9\sigma_* N_0}
\left(
\delta N_{A,0}
-
P_0\,\delta D_{A,0}
\right).
}
\]

These are the two main formulas of the stage.

---

## 3. Exact one-parameter even-consistency relation

A direct hidden even pole deformation is a **one-parameter** object: once
`\delta\kappa_W^{(A)}` is fixed, both even grouped defects are fixed.

Indeed the Stage-142 ratio gives
\[
\delta E_4^{(A)}=\frac89\,\delta E_2^{(A)}.
\]
Since
\[
\delta E_2^{(A)}=\delta u_2^{(A)},
\qquad
\delta E_4^{(A)}=\delta u_4^{(A)},
\]
the grouped conservative moments must satisfy
\[
\boxed{
\delta u_4^{(A)}=\frac89\,\delta u_2^{(A)}.
}
\]

Using the explicit transport laws above, this is equivalent to
\[
-\frac{\delta D_{A,4}+\frac29\,\delta D_{A,2}+\frac{5}{81}\,\delta D_{A,0}}{D_0}
=
\frac89
\left[
-\frac{\delta D_{A,2}+\frac19\,\delta D_{A,0}}{D_0}
\right].
\]
So the exact microscopic even-consistency relation is
\[
\boxed{
\delta D_{A,4}
=
\frac23\,\delta D_{A,2}
+\frac{1}{27}\,\delta D_{A,0}.
}
\]

This is the direct grouped-lane criterion that says a given linear conservative anisotropy can really be interpreted as a hidden even-pole deformation `\delta\kappa_W^{(A)}` and not as some more general even outlet distortion.

So the linear grouped-even problem has collapsed twice:

1. first to the single combination
   \[
   \delta D_{A,2}+\frac19\,\delta D_{A,0},
   \]
2. and second to the exact `u_4` compatibility relation above.

---

## 4. Grouped trace/anomaly form

Let
\[
a_{D,0},\ b_{D,0},
\qquad
a_{D,2},\ b_{D,2},
\qquad
a_{D,4},\ b_{D,4},
\qquad
a_{N,0},\ b_{N,0}
\]
be the grouped anisotropy defects of the lane families `D_0,D_2,D_4,N_0`.

Then the exact Stage-6 projector formulas immediately give
\[
\boxed{
a_\kappa
=
\frac{3(1-\sigma_*)}{\sigma_* D_0}
\left(
a_{D,2}+\frac19\,a_{D,0}
\right),
\qquad
b_\kappa
=
\frac{3(1-\sigma_*)}{\sigma_* D_0}
\left(
b_{D,2}+\frac19\,b_{D,0}
\right),
}
\]
where `a_\kappa,b_\kappa` are the grouped anisotropy defects of the hidden even pole.

Similarly,
\[
\boxed{
a_\gamma
=
-\frac{1-\sigma_*}{9\sigma_* N_0}
\left(
a_{N,0}-P_0 a_{D,0}
\right),
\qquad
b_\gamma
=
-\frac{1-\sigma_*}{9\sigma_* N_0}
\left(
b_{N,0}-P_0 b_{D,0}
\right),
}
\]
where `a_\gamma,b_\gamma` are the grouped anisotropy defects of the hidden odd normalization.

And the exact even-consistency relation becomes
\[
\boxed{
a_{D,4}
=
\frac23\,a_{D,2}
+\frac{1}{27}\,a_{D,0},
\qquad
b_{D,4}
=
\frac23\,b_{D,2}
+\frac{1}{27}\,b_{D,0}.
}
\]

So the remaining linear grouped-anisotropy problem is completely transparent in the grouped projector language.

---

## 5. Weak axisymmetric branch

On the weak axisymmetric `Y_{20}` branch write
\[
\delta D_{A,n}=\epsilon\,\lambda_A\,D_n^{(1)},
\qquad
\delta N_{A,0}=\epsilon\,\lambda_A\,N_0^{(1)},
\]
with the standard grouped pattern
\[
(\lambda_{20},\lambda_{21},\lambda_{22})=\left(1,\frac12,-1\right).
\]

Then the direct outlet deformations inherit the same grouped signature:
\[
\boxed{
\delta\kappa_W^{(20)}=\epsilon\,\kappa_1,
\qquad
\delta\kappa_W^{(21)}=\frac\epsilon2\,\kappa_1,
\qquad
\delta\kappa_W^{(22)}=-\epsilon\,\kappa_1,
}
\]
with
\[
\boxed{
\kappa_1
=
\frac{3(1-\sigma_*)}{\sigma_* D_0}
\left(
D_2^{(1)}+\frac19\,D_0^{(1)}
\right),
}
\]
and
\[
\boxed{
\delta\gamma_W^{(20)}=\epsilon\,\gamma_1,
\qquad
\delta\gamma_W^{(21)}=\frac\epsilon2\,\gamma_1,
\qquad
\delta\gamma_W^{(22)}=-\epsilon\,\gamma_1,
}
\]
with
\[
\boxed{
\gamma_1
=
-\frac{1-\sigma_*}{9\sigma_* N_0}
\left(
N_0^{(1)}-P_0 D_0^{(1)}
\right).
}
\]

So the weak axisymmetric grouped-lane outlet problem collapses to **two scalar amplitudes only**:

- one even amplitude `\kappa_1`,
- one odd amplitude `\gamma_1`.

That is much smaller than the original grouped bundle.

---

## 6. What Stage 153 changes

Before this step, Stage 152 had shown only that the remaining **linear**
grouped-anisotropy problem sits in the direct outlet coefficients
\[
\delta\kappa_W,\qquad \delta\gamma_W,
\]
rather than in the scalar slippages.

After this step, that statement is now explicit.

The theorem status is:

1. the linear grouped-even hidden-pole deformation is controlled only by
   \[
   \mathcal K_A
   :=
   \delta D_{A,2}+\frac19\,\delta D_{A,0},
   \]
   together with the exact compatibility condition
   \[
   \delta D_{A,4}
   =
   \frac23\,\delta D_{A,2}
   +\frac1{27}\,\delta D_{A,0};
   \]

2. the linear grouped-odd normalization deformation is controlled only by
   \[
   \mathcal G_A
   :=
   \delta N_{A,0}-P_0\,\delta D_{A,0};
   \]

3. therefore the full linear grouped-anisotropy outlet bottleneck has collapsed to the pair
   \[
   (\mathcal K_A,\mathcal G_A).
   \]

So the next honest theorem gate is no longer

> “derive the linear grouped-anisotropy problem.”

That part is now finished.

The next gate is:

> compute the actual moving-throat bundle anisotropies
> \[
> (\delta D_{A,0},\delta D_{A,2},\delta D_{A,4},\delta N_{A,0})
> \]
> and test whether they satisfy
> \[
> \mathcal K_A=0,\qquad \mathcal G_A=0
> \]
> on the physical branch.

That is the direct continuation point.
