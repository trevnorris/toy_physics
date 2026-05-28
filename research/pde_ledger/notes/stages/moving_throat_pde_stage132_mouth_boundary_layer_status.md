
# Moving-Throat PDE — Stage 132: Mouth Boundary-Layer Status After Explicit Source-Law Extraction

## What is now fixed

After Stages 129–131, the mouth-source side is no longer an abstract profile problem.

1. The actual positive source family is fixed by an explicit GNLS + localized-Maxwell
   boundary-layer law:
   \[
   \sigma_{\Pi}(z)=\frac{\Pi e^{-\Pi z/L}}{L(1-e^{-\Pi})}.
   \]

2. The corresponding Family-1 mouth-bias factor is exact:
   \[
   \mathfrak g_{\Pi}
   =
   \frac{2\Pi(2\Pi e^{\Pi}+\pi)}
        {(4\Pi^2+\pi^2)(e^{\Pi}-1)}.
   \]

3. This family is strictly monotone:
   \[
   \mathfrak g_{\Pi}: \frac{2}{\pi}\to 1
   \qquad(\Pi:0^+\to+\infty).
   \]

4. The unique canonical Family-1 compensation point is
   \[
   \Pi_* \approx 1.50882951349316.
   \]

5. The exact parent threshold is
   \[
   \left.\partial_z\delta V_{\rm conf}\right|_{\rm m}
   -
   q_*\left.\partial_zA_0\right|_{\rm m}
   =
   1.50882951349316\,\frac{\Theta_\sigma}{L}.
   \]

## Meaning

The outlet-core ambiguity has narrowed again.

The open question is no longer:

- which compensation branch is physical,
- whether positive localized sources can hit it,
- or what source family to use.

It is now simply:

\[
\boxed{
\text{what value of }
\Pi_m=\frac{L}{\Theta_\sigma}\left(
\left.\partial_z\delta V_{\rm conf}\right|_{\rm m}
-q_*\left.\partial_zA_0\right|_{\rm m}
\right)
\text{ does the real mouth layer select?}
}
\]

And the target value is moderate:
\[
\Pi_m\approx1.51.
\]

So the remaining gap is now one clean microscopic bias law, not a diffuse branch-selection problem.
