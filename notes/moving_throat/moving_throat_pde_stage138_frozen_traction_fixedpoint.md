# Moving-Throat PDE — Stage 138: Family-1 Co-Evolving Fixed Point at Frozen Canonical Traction

## Goal

Solve the exact co-evolving core–mouth map on the explicit Family-1 branch while
holding the previously selected canonical traction fixed:
\[
\Sigma_0=\Sigma_0^*
\approx 1.80594111095636
\qquad
\left(
\widehat T_m=\widehat T_{m,*}
\approx 0.901484054174204
\right).
\]

This answers a very specific question:

> If the physical mouth traction is left at the old canonical value, does the
> fully co-evolving Family-1 fixed point still land on the lower compensated
> branch \(\mathcal R=1/4\)?

---

## 1. Positive fixed point on the analyzed branch window

Iterating the exact nonlinear map
\[
\Sigma\mapsto
\frac{e^{-\Phi_{\Sigma_0^*}[\Sigma]}}
{\int_0^1 e^{-\Phi_{\Sigma_0^*}[\Sigma](y)}dy}
\]
from the canonical exponential seed converges to a positive fixed point
\(\Sigma_{\rm fp}\) on the explicit Family-1 branch. On the analyzed branch window,
this fixed point is unique.

Its selected moments are
\[
\boxed{
\mathfrak g_{\rm fp}
\approx 0.693352419668063,
\qquad
\mathcal S_{\rm fp}
\approx 0.6216013167514007.
}
\]

So relative to the exact lower compensated value
\[
\mathfrak g_*\approx 0.758035078944663,
\]
the fixed profile is broader:
\[
\boxed{
\delta\mathfrak g_{\rm fp}
=
\mathfrak g_{\rm fp}-\mathfrak g_*
\approx -0.0646826592766000<0.
}
\]

---

## 2. Co-evolving core ratio and mouth slope

Feeding that back into the exact Family-1 core law gives
\[
\boxed{
\mathcal R_{\rm fp}
=
\mathcal R[\Sigma_{\rm fp}]
\approx 0.2827139049082381.
}
\]

So the co-evolving fixed point **does not** stay on the compensated value
\(1/4\); it shifts upward by
\[
\boxed{
\delta\mathcal R_{\rm fp}
=
\mathcal R_{\rm fp}-\frac14
\approx 0.0327139049082381.
}
\]

The exact slope identity then gives
\[
\boxed{
\Pi_{\rm fp}
=
\Sigma_0^*
\Bigl[1-\mathcal R_{\rm fp}\mathcal S_{\rm fp}\Bigr]
\approx 1.4885734438300713.
}
\]

So at fixed canonical traction the co-evolving mouth profile actually stays very
close to the old canonical bias:
\[
\boxed{
\delta\Pi_{\rm fp}
=
\Pi_{\rm fp}-\Pi_*
\approx -0.0202560696630887.
}
\]

---

## 3. Exact first-order transport check

Using the exact transport law from Stage 137,
\[
\delta\mathcal R
=
-\frac{\delta\mathfrak g}{\sqrt{1+\mathfrak r_{F1}^2}}
+
\frac{(\delta\mathfrak g)^2}{1+\mathfrak r_{F1}^2},
\]
the observed fixed-point broadening predicts
\[
\delta\mathcal R_{\rm pred}
=
-\frac{\delta\mathfrak g_{\rm fp}}{\sqrt{1+\mathfrak r_{F1}^2}}
+
\frac{(\delta\mathfrak g_{\rm fp})^2}{1+\mathfrak r_{F1}^2}
\approx 0.0327139049082381,
\]
which matches the direct fixed-point value.

So the fixed-point drift is exactly the one implied by the co-evolving core law.

---

## 4. Meaning

The full Family-1 core–mouth solve gives a more nuanced answer than the earlier
fixed-core mouth correction:

- the source profile really does broaden,
- but the induced increase in the mixed loading ratio pushes back on the mouth bias,
- so at **fixed** canonical traction the self-consistent branch remains close in
  \(\Pi\), even though it is no longer exactly compensated in \(\mathfrak g\).

So the co-evolution does **not** destroy the regular Family-1 branch.
What it does is convert the old exact compensation point into a nearby
traction-dependent fixed point with
\[
\mathcal R_{\rm fp}>1/4.
\]

That means exact preservation of the canonical outgoing quadrupole fingerprint now
requires a **retuned traction**, not just the old canonical value.
