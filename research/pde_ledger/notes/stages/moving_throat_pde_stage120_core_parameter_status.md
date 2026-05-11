
# Moving-Throat PDE — Stage 222: Core-Parameter Extraction Status

## Result

The concrete throat-core outlet is no longer described by free reduced coefficients.

After explicit GNLS + localized-Maxwell reduction, the surviving microscopic controls are:

- shell stiffness \(K_s\),
- mixed-tube stiffness \(K_q\),
- normalized shell/mixed hybridization
  \[
  \mathfrak r=\frac{\lambda}{\sqrt{K_sK_q}},
  \]
- normalized mouth-coupling ratio
  \[
  \mathfrak g=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}},
  \]
- and the D/N mixed-tube length \(L_W\).

## The sharp surviving theorem gate

The compensated canonical outgoing quadrupole branch exists **iff**
\[
1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2,
\]
with
\[
L_W=\frac{\pi a}{2}\sqrt{\frac{1+\mathfrak r^2}{3}}.
\]

So the surviving microscopic question is no longer:

> “What are the four reduced outlet coefficients?”

It is now:

> “What branch values of \((\mathfrak r,\mathfrak g)\) does the actual GNLS + localized-Maxwell throat core select?”

That is a substantially smaller target for the next derivation.
