# 2PN Family-1 anisotropic tangential-traction completion: current result

## 1. What was tested

The strict isotropic Family-1 boundary-layer pullback failed even after exact flare geometry was included.  
The previous step showed that the missing support-sector structure could be rewritten as a local wall-stress / traction-balance profile, but that still left open whether a **real PDE-side boundary-layer completion** could generate it.

The minimal new test was:

- keep the **normal penetration moment isotropic**,
- keep the same strict quadratic-flare Family-1 surface pullback,
- and promote only the tangential wall moment from a constant to the smallest axisymmetric profile
  \[
  B_{\rm tan}(\mu)=B_0+B_2\mu^2.
  \]

So the surface energy becomes
\[
E_{\rm BL}^{\rm aniso}
=
\frac12\int d\Omega\,
\Big[
A\,J(\mu)\,\Psi^2
+
B_{\rm tan}(\mu)\Big(
F_\theta(\mu)(\partial_\theta\Psi)^2
+
F_\phi(\mu)\frac{(\partial_\phi\Psi)^2}{\sin^2\theta}
\Big)
\Big],
\]
with the same strict pullback factors \(J,F_\theta,F_\phi\) used in the earlier no-go step.

This is the smallest genuine PDE-side completion of the strict model.

---

## 2. Exact support-channel closure

Projecting onto the \(\ell=1,2\) channels gives exact formulas for
\[
K_{1\perp},\ K_{10},\ K_{20},\ K_{21},\ K_{22}
\]
that are linear in \(A,B_0,B_2\) and polynomial in \((q,r)\).

A multi-seed solve shows that the system
\[
K_{1\perp}=\frac27,\qquad
K_{10}=\frac14,\qquad
K_{20}=\frac49,\qquad
K_{21}=\frac23,\qquad
K_{22}=\frac83
\]
has **exact real solutions**.

So the strict isotropic no-go is really repaired by one new local degree of freedom:

> an axisymmetric tangential wall-stress profile \(B_{\rm tan}(\mu)=B_0+B_2\mu^2\).

This is the first strict boundary-layer / soft-wall PDE completion that actually reproduces the solved support sector, rather than fitting the support operator directly.

---

## 3. Family-1-like physical branch

The SymPy branch scan finds at least four exact real branches.  
The natural Family-1-like branch is the one with

- \(q>0\),
- \(r>0\),
- and \(\sigma(\mu)=1-q\mu^2+r\mu^4>0\) on \([-1,1]\).

For that branch,
\[
A\approx -0.281923219302567,
\]
\[
B_0\approx 0.648914709228733,
\qquad
B_2\approx -1.085434603876673,
\]
\[
q\approx 2.370915717168574,
\qquad
r\approx 2.758343474052255.
\]

The support residuals vanish at machine precision.

In Legendre form,
\[
B_{\rm tan}(\mu)=\beta_0+\beta_2 P_2(\mu),
\]
with
\[
\beta_0=B_0+\frac13 B_2 \approx 0.287103174603175,
\]
\[
\beta_2=\frac23 B_2 \approx -0.723623069251115.
\]

So the minimal successful PDE-side completion is an **axisymmetric \(P_2\)-modulated tangential traction profile**.

---

## 4. Exact universal monopole prediction

The most useful structural result is that the anisotropic tangential model satisfies the exact identity
\[
K_{00}
=
K_{1\perp}
+\frac12 K_{10}
-\frac1{10}K_{20}
-\frac15 K_{21}
-\frac15 K_{22}.
\]

So once the \(\ell=1,2\) support targets are matched exactly, the raw boundary-layer monopole is fixed automatically:
\[
K_{00}^{\rm BL}
=
\frac27+\frac12\cdot\frac14-\frac1{10}\cdot\frac49-\frac15\cdot\frac23-\frac15\cdot\frac83
=
-\frac{757}{2520}.
\]

This is independent of which exact support branch is chosen.

Since the carried-forward monopole target is
\[
K_{00}^{\rm target}=\frac4{45},
\]
the separate monopole wall add-on must be
\[
\Delta K_{00}^{\rm mono}
=
\frac4{45}-\left(-\frac{757}{2520}\right)
=
\frac{109}{280}.
\]

So the support-sector repair does **not** eliminate the separate monopole channel; instead it predicts its exact required size.

---

## 5. Physical reading

This is an important narrowing of the throat PDE problem.

The previous exact operator and traction-balance steps said that a missing tangential wall-stress profile had to exist.  
This new result shows that the smallest strict boundary-layer realization of that statement is already enough:

1. isotropic normal penetration moment \(A\),
2. one axisymmetric tangential wall-moment profile \(B_{\rm tan}(\mu)=B_0+B_2\mu^2\),
3. the same strict Family-1 flare pullback,
4. plus the still-separate monopole wall channel.

So the next PDE target is no longer vague.

The next thing to derive from the soft-wall physics is:

- why the tangential moment picks up exactly this \(P_2\)-modulated profile, and
- how the separate monopole add-on \(\frac{109}{280}\) emerges from the isotropic wall / geometry sector.

That is a much sharper finish-line target than “derive the whole support operator from scratch.”
