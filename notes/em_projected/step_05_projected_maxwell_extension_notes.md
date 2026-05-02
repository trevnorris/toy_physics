# Projected Maxwell Extension Notes

## What this note covers

This note extends the earlier projected-Maxwell derivation by including a **weighted gauge-fixing profile** \(H(w)\) in the localized \(4+1\) Maxwell sector.

The goal was to answer a precise question:

> If we keep the theory in **projection-first** form rather than going directly to the usual zero-mode reduction, how does the gauge-fixing sector look, and what changes when we compare \(H=1\) against \(H=Z\)?

It also sharpens a second issue:

> Under what conditions does the projection-first zero-mode equation actually collapse to the same effective Maxwell law as the reduction-first derivation?

---

## Starting point

The generalized localized bulk Maxwell equation is taken to be

\[
\partial_M\!\bigl(Z(w)F^{MN}\bigr)
+\frac{1}{\xi}\,\partial^N\!\bigl(H(w)B\bigr)
=\mu_0 J^N,
\qquad
B\equiv \partial\!\cdot\!A.
\]

Here:

- \(Z(w)\) localizes the Maxwell kinetic term.
- \(H(w)\) localizes the gauge-fixing term.
- The old unweighted choice is \(H=1\).
- The cleaned-up localized choice is \(H=Z\).

---

## Exact projected brane equation

Project with a normalized observation kernel \(W(w)\):

\[
\operatorname{Proj}_W[Q](x)=\int_{-\infty}^{\infty}W(w)Q(x,w)\,dw,
\qquad
\int W(w)\,dw=1.
\]

For a brane index \(\nu\in\{0,1,2,3\}\), the projected inhomogeneous law is

\[
\partial_\mu\operatorname{Proj}_W[Z F^{\mu\nu}]
+\operatorname{Boundary}[WZF^{w\nu}]
-\operatorname{Proj}_{W'}[Z F^{w\nu}]
+\frac{1}{\xi}\partial^\nu \operatorname{Proj}_W[H B]
=\mu_0\operatorname{Proj}_W[J^\nu].
\]

Under boundary decay in \(w\), this simplifies to

\[
\partial_\mu\operatorname{Proj}_W[Z F^{\mu\nu}]
-\operatorname{Proj}_{W'}[Z F^{w\nu}]
+\frac{1}{\xi}\partial^\nu \operatorname{Proj}_W[H B]
=\mu_0\operatorname{Proj}_W[J^\nu].
\]

### Immediate interpretation

The new weight \(H(w)\) only changes the **projected gauge-driver term**

\[
\frac{1}{\xi}\partial^\nu \operatorname{Proj}_W[H B].
\]

It does **not** change:

- the homogeneous Maxwell/Bianchi structure,
- the leakage term from \(ZF^{w\nu}\),
- or the basic projection/open-system character of the brane theory.

---

## Zero-mode projection-first effective equation

Now impose the usual zero-mode/far-field ansatz

\[
A_\mu(x,w)=a_\mu(x),
\qquad
A_w=0,
\qquad
\partial_w A_\mu=0,
\qquad
F^{w\nu}=0,
\qquad
J^\nu(x,w)=j^\nu(x)S(w).
\]

Then \(B=\partial\!\cdot\!a\) is independent of \(w\), and the projected equation becomes

\[
I_{WZ}\,\partial_\mu f^{\mu\nu}
+\frac{I_{WH}}{\xi}\,\partial^\nu(\partial\!\cdot\!a)
=\mu_0 I_{WS} j^\nu,
\]

with

\[
I_{WZ}=\int WZ\,dw,
\qquad
I_{WH}=\int WH\,dw,
\qquad
I_{WS}=\int WS\,dw.
\]

Dividing by \(I_{WZ}\) gives the projection-first effective Maxwell form

\[
\partial_\mu f^{\mu\nu}
+\frac{1}{\xi_{\rm eff}^{\rm proj}}\partial^\nu(\partial\!\cdot\!a)
=\mu_0^{\rm eff,proj} j^\nu,
\]

with

\[
\boxed{
\mu_0^{\rm eff,proj}=\mu_0\frac{I_{WS}}{I_{WZ}}
}
\]

and

\[
\boxed{
\frac{1}{\xi_{\rm eff}^{\rm proj}}=\frac{I_{WH}}{\xi I_{WZ}},
\qquad
\xi_{\rm eff}^{\rm proj}=\xi\frac{I_{WZ}}{I_{WH}}.
}
\]

This is the cleanest new result of the extension:

> In the projection-first zero-mode theory, there is not only an effective coupling \(\mu_0^{\rm eff}\); there is also an effective gauge parameter \(\xi_{\rm eff}\).

---

## Clean comparison: \(H=1\) versus \(H=Z\)

Because \(W\) is normalized,

\[
\int W(w)\,dw=1.
\]

### Case 1: unweighted gauge fixing \(H=1\)

Then

\[
I_{WH}=1,
\]

so

\[
\boxed{
\xi_{\rm eff}^{\rm proj}(H=1)=\xi I_{WZ}.
}
\]

So the projected gauge parameter depends on the observer kernel through \(I_{WZ}=\int WZ\,dw\).

### Case 2: localized gauge fixing \(H=Z\)

Then

\[
I_{WH}=I_{WZ},
\]

so

\[
\boxed{
\xi_{\rm eff}^{\rm proj}(H=Z)=\xi.
}
\]

This is exact for **any normalized** projection kernel \(W\).
The updated script proves this by starting from the symbolic projected integral
\(I_{WH}=\int WH\,dw\) and performing the substitution \(H\mapsto Z\) inside
SymPy, both in the integral itself and in the projected gauge-driver term,
rather than by typing \(I_{WH}=I_{WZ}\) from the start.

### Main gauge-sector conclusion

> The choice \(H=Z\) is structurally cleaner in projection-first form because it preserves the same gauge parameter \(\xi\) after projection.

By contrast, \(H=1\) leaves the projected gauge-driver term observer-kernel dependent.

---

## When projection-first matches reduction-first exactly

The effective coupling is

\[
\mu_0^{\rm eff,proj}=\mu_0\frac{I_{WS}}{I_{WZ}},
\]

so matching depends on the **source profile** \(S(w)\).

### Matching source profile

If the projected source follows the normalized localization profile

\[
S(w)=\frac{Z(w)}{Z_{\rm int}},
\qquad
Z_{\rm int}=\int Z(w)\,dw,
\]

then

\[
I_{WS}=\frac{I_{WZ}}{Z_{\rm int}},
\]

and therefore

\[
\boxed{
\mu_0^{\rm eff,proj}=\frac{\mu_0}{Z_{\rm int}}.
}
\]

This matches the reduction-first Maxwell coupling exactly, for **any** normalized \(W\).

### Full matching channel

If we also choose localized gauge fixing

\[
H=Z,
\]

then simultaneously

\[
\boxed{
\mu_0^{\rm eff,proj}=\frac{\mu_0}{Z_{\rm int}},
\qquad
\xi_{\rm eff}^{\rm proj}=\xi.
}
\]

So the projection-first zero-mode equation becomes

\[
\partial_\mu f^{\mu\nu}
+\frac{1}{\xi}\partial^\nu(\partial\!\cdot\!a)
=\frac{\mu_0}{Z_{\rm int}}j^\nu,
\]

which is exactly the same gauge-fixed Maxwell form as the localized-gauge reduction.

### Meaning

> Projection-first and reduction-first are not automatically the same, but there is an exact matching channel if both the source profile and gauge-fixing profile are aligned with the localized zero mode.

---

## Gaussian worked example

Take

\[
Z(w)=e^{-w^2/\lambda^2},
\qquad
Z_{\rm int}=\sqrt{\pi}\,\lambda,
\qquad
Z_{2,\rm int}=\int Z^2dw=\frac{\sqrt{2\pi}\,\lambda}{2}.
\]

Choose the matched observer kernel

\[
W(w)=\frac{Z(w)}{Z_{\rm int}}.
\]

Then

\[
I_{WZ}=\frac{Z_{2,\rm int}}{Z_{\rm int}}=\frac{1}{\sqrt{2}}.
\]

### Gauge sector

- For \(H=1\):
  \[
  \xi_{\rm eff}^{\rm proj}=\frac{\xi}{\sqrt{2}}.
  \]
- For \(H=Z\):
  \[
  \xi_{\rm eff}^{\rm proj}=\xi.
  \]

### Source sector A: delta-localized source

If the source is sharp on the brane,

\[
S(w)=\delta(w),
\]

then

\[
I_{WS}=W(0)=\frac{1}{\sqrt{\pi}\,\lambda},
\]

so

\[
\mu_0^{\rm eff,proj}
=\mu_0\frac{I_{WS}}{I_{WZ}}
=\frac{\sqrt{2}\,\mu_0}{\sqrt{\pi}\,\lambda}
=\frac{\mu_0}{Z_{2,\rm int}}.
\]

This is larger than the reduction-first value by a factor of \(\sqrt{2}\).

### Source sector B: matched distributed source

If instead

\[
S(w)=\frac{Z(w)}{Z_{\rm int}},
\]

then

\[
\mu_0^{\rm eff,proj}=\frac{\mu_0}{\sqrt{\pi}\,\lambda}=\frac{\mu_0}{Z_{\rm int}},
\]

which exactly matches the reduction-first result.

---

## Comparison with reduction-first gauge fixing

For localized action reduction with a finite gauge-fixing weight \(H(w)\), the reduced 4D gauge parameter is

\[
\xi_4=\xi\frac{Z_{\rm int}}{H_{\rm int}},
\qquad
H_{\rm int}=\int H(w)\,dw.
\]

### If \(H=Z\)

Then

\[
\xi_4=\xi.
\]

### If \(H=1\) on a regulator box \([-R,R]\)

Then

\[
\xi_4(R)=\xi\frac{Z_{\rm int}}{2R}
\longrightarrow 0
\qquad (R\to\infty).
\]

So the unweighted gauge-fixing term is unsafe as a noncompact zero-mode gauge-fixed action unless Lorenz gauge is imposed before reduction.

### Conceptual difference

- **Reduction-first** with \(H=1\): gauge-fixing is a bulk device; do not naively reduce it as a noncompact zero-mode action.
- **Projection-first** with \(H=1\): the projection stays finite because \(W\) is normalized, but the resulting \(\xi_{\rm eff}^{\rm proj}\) depends on \(W\).
- **Both pictures** become much cleaner for \(H=Z\).

---

## Final takeaways

### 1. Weighted gauge fixing can be carried through projection exactly

The exact projected law is

\[
\partial_\mu\operatorname{Proj}_W[ZF^{\mu\nu}]
-\operatorname{Proj}_{W'}[ZF^{w\nu}]
+\frac{1}{\xi}\partial^\nu\operatorname{Proj}_W[HB]
=\mu_0\operatorname{Proj}_W[J^\nu].
\]

### 2. Projection-first has two effective parameters

\[
\mu_0^{\rm eff,proj}=\mu_0\frac{\int WS}{\int WZ},
\qquad
\xi_{\rm eff}^{\rm proj}=\xi\frac{\int WZ}{\int WH}.
\]

### 3. \(H=Z\) is the clean gauge-sector choice

It gives

\[
\xi_{\rm eff}^{\rm proj}=\xi
\]

for any normalized \(W\).

### 4. Matching the coupling is a source-profile issue

If

\[
S=Z/Z_{\rm int},
\]

then

\[
\mu_0^{\rm eff,proj}=\mu_0/Z_{\rm int}
\]

for any normalized \(W\).

### 5. Exact projection/reduction agreement is possible

With

\[
H=Z,
\qquad
S=Z/Z_{\rm int},
\]

the projection-first zero-mode equation reproduces the same gauge-fixed Maxwell law as the localized reduction-first derivation.

### 6. Delta-localized sources still produce a mismatch

If the source is sharp in \(w\), the projection-first coupling generally differs from the reduction-first coupling even when \(H=Z\).

---

## Audit script

- `step_05_projected_maxwell_extension_sympy.py`

Run this script to make the derivation auditable and easy to extend later.
