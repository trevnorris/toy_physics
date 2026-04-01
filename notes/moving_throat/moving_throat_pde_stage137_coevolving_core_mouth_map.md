# Moving-Throat PDE — Stage 137: Exact Co-Evolving Core–Mouth Fixed-Point Map

## Goal

Replace the “fixed compensated core + corrected mouth source” closure by a **single
self-consistent co-evolving core–mouth map**.

The point is to let the actual positive mouth source profile modify the shell/mixed
loading ratio, and then feed that modified ratio back into the mouth potential.

---

## 1. The exact co-evolving Family-1 map

For any normalized positive mouth source
\[
\Sigma(x)\ge 0,
\qquad
\int_0^1 \Sigma(x)\,dx = 1,
\]
define the two mouth moments
\[
\mathfrak g[\Sigma]
=
\int_0^1 \Sigma(x)\cos\!\left(\frac{\pi x}{2}\right)dx,
\qquad
\mathcal S[\Sigma]
=
\int_0^1 \Sigma(x)\,
\frac{\cosh\!\left(\frac{\pi}{2}(1-x)\right)}{\cosh(\pi/2)}\,dx.
\]

On the explicit Family-1 core branch, the shell/mixed loading ratio is not free.
It is the exact core overlap law
\[
\boxed{
\mathcal R[\Sigma]
=
\frac{\bigl(\mathfrak g[\Sigma]-\mathfrak r_{F1}\bigr)^2}{1+\mathfrak r_{F1}^2},
\qquad
\mathfrak r_{F1}\approx 1.77799353547498.
}
\]

For the static shell channel,
\[
\mathcal T_s[\Sigma](x)=\int_0^1 \min(x,y)\,\Sigma(y)\,dy,
\]
and for the first mixed D/N half-wave,
\[
\mathcal T_q[\Sigma](x)
=
\int_0^1
\frac{\sinh\!\left(\frac{\pi}{2}\min(x,y)\right)
\cosh\!\left(\frac{\pi}{2}(1-\max(x,y))\right)}{(\pi/2)\cosh(\pi/2)}
\Sigma(y)\,dy.
\]

So the full co-evolving mouth potential is
\[
\boxed{
\Phi_{\Sigma_0}[\Sigma](x)
=
\Sigma_0
\Big[
\mathcal T_s[\Sigma](x)
-
\mathcal R[\Sigma]\,\mathcal T_q[\Sigma](x)
\Big].
}
\]

The self-consistent source law is therefore the nonlinear fixed-point equation
\[
\boxed{
\Sigma(x)
=
\frac{e^{-\Phi_{\Sigma_0}[\Sigma](x)}}
{\int_0^1 e^{-\Phi_{\Sigma_0}[\Sigma](y)}dy}.
}
\]

This is the first explicit Family-1 mouth problem in which the **source profile and
the core loading ratio co-evolve together**.

---

## 2. Canonical compensation inside the co-evolving map

The lower compensated throat-core branch remains the condition
\[
\mathfrak g[\Sigma]=\mathfrak g_*,
\qquad
\mathfrak g_*\approx 0.758035078944663.
\]

Because
\[
\mathfrak g_*=
\mathfrak r_{F1}-\frac12\sqrt{1+\mathfrak r_{F1}^2},
\]
this is exactly equivalent to
\[
\boxed{
\mathcal R[\Sigma]=\frac14.
}
\]

So the canonical outgoing quadrupole branch survives in the co-evolving theory iff
the self-consistent source profile lands back on the same lower compensation moment.

---

## 3. Exact first-order defect transport

Write
\[
\mathfrak g=\mathfrak g_*+\delta\mathfrak g.
\]
Then the exact Family-1 ratio becomes
\[
\boxed{
\mathcal R
=
\frac14
-
\frac{\delta\mathfrak g}{\sqrt{1+\mathfrak r_{F1}^2}}
+
\frac{(\delta\mathfrak g)^2}{1+\mathfrak r_{F1}^2}.
}
\]

So to first order,
\[
\boxed{
\delta\mathcal R
=
-\frac{\delta\mathfrak g}{\sqrt{1+\mathfrak r_{F1}^2}}
+O(\delta\mathfrak g^2).
}
\]

Numerically,
\[
\sqrt{1+\mathfrak r_{F1}^2}
\approx 2.039916913060632,
\qquad
\delta\mathcal R
\approx -0.490215\,\delta\mathfrak g.
\]

So broadening the source (negative \(\delta\mathfrak g\)) automatically drives the
mixed loading ratio **above** its compensated value \(1/4\).

---

## 4. Local slope / bias identity

For any self-consistent source on this branch, the actual mouth slope is
\[
\boxed{
\Pi[\Sigma]
=
\Phi'_{\Sigma_0}[\Sigma](0)
=
\Sigma_0\Bigl[1-\mathcal R[\Sigma]\,\mathcal S[\Sigma]\Bigr].
}
\]

Hence the co-evolving canonical branch is determined by the coupled conditions

\[
\Sigma
=
\frac{e^{-\Phi_{\Sigma_0}[\Sigma]}}
{\int e^{-\Phi_{\Sigma_0}[\Sigma]}},
\qquad
\mathfrak g[\Sigma]=\mathfrak g_*,
\qquad
\Pi=\Sigma_0\left(1-\frac14\mathcal S[\Sigma]\right).
\]

Under the self-matched susceptibility closure from Stage 123,
\[
\boxed{
\Sigma_0=\frac{20}{9}\widehat T_m^2.
}
\]

So the remaining ambiguity is now only the required **normalized mouth traction**
\(\widehat T_m\) that makes the co-evolving fixed point land on
\(\mathfrak g_*\).
