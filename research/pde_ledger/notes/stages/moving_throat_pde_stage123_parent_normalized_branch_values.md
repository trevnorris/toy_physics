# Moving-Throat PDE — Stage 123: Parent-Normalized Branch Values on the Explicit Core Solve

## Goal

Translate the explicit branch data from Stages 121–122 into actual parent-level flow and traction requirements.

## 1. Parent-normalized background-flow number

Stage 119 gives
\[
\mathfrak r=\frac{\lambda}{\sqrt{K_sK_q}},
\qquad
\lambda=-q_*v_{w0}\mathcal I_{sq}.
\]
Using the healing-locked shell formulas from Stage 118 and identifying the mixed D/N tube with the actual throat length \(L_W=L\), define the parent-normalized mixed-flow number
\[
\boxed{
\Xi_v
:=
\frac{q_*\,\sqrt{\mu_0 m_\psi\rho_w}\;aL^{3/2}\ell^{3/2}}
{\hbar\sqrt{\mathcal Z_q}}\,
\frac{v_{w0}}{c_s}.
}
\]
Then the exact branch law is
\[
\boxed{
\Xi_v
=
-\frac{3\sqrt{30}\,\pi^{3/2}}{160}\,\mathfrak r.
}
\]

On the explicit Family-1 geometric branch,
\[
\boxed{
\Xi_v^{F1}
=
-\frac{3\sqrt{30}\,\pi^{3/2}}{160}\,\mathfrak r_{F1}
\approx -1.01675633282526.
}
\]

So the actual mixed background flow required by the compensated throat-core branch is already `O(1)` in normalized parent variables.

## 2. Parent-normalized traction number

Stage 119 also gives
\[
\mathfrak g
=
\frac{\sqrt{2\mathcal Z_qK_s}}
{\mathcal T_m J_s\,c_s\sqrt{\mu_0 L}}.
\]
Using the same healing-locked shell formulas, define the parent-normalized mouth-traction number
\[
\boxed{
\Xi_T
:=
\mathcal T_m\,
\frac{a\sqrt{\mu_0\rho_wL\ell}}
{\sqrt{\mathcal Z_q m_\psi}}.
}
\]
Then the exact branch law is
\[
\boxed{
\Xi_T
=
\frac{3\sqrt{30}}{10\sqrt{\pi}}\,
\frac{1}{\mathfrak g}.
}
\]

So the natural equal-normalized mouth-source branch \(\mathfrak g_{\rm nat}=1\) gives
\[
\boxed{
\Xi_T^{\rm nat}
=
\frac{3\sqrt{30}}{10\sqrt{\pi}}
\approx 0.927058084855655.
}
\]

On the two compensated branches,
\[
\boxed{
\Xi_T^{(\pm)}
=
\frac{3\sqrt{30}}{10\sqrt{\pi}}\,
\frac{1}{\mathfrak g_\pm^{F1}}.
}
\]
Numerically,
\[
\boxed{
\Xi_T^{(-)}\approx 1.22297517701464,
\qquad
\Xi_T^{(+)}\approx 0.331334521644609.
}
\]

## 3. Interpretation

These normalized branch values sharpen the outlet-core question a lot.

The explicit moving-throat branch now says:

- the geometry-selected normalized hybridization is
  \[
  \Xi_v^{F1}\approx -1.017;
  \]
- the naive equal-normalized mouth source would give
  \[
  \Xi_T^{\rm nat}\approx 0.927;
  \]
- the nearest compensated branch instead needs
  \[
  \Xi_T^{(-)}\approx 1.223
  \]

So the concrete throat-core problem is no longer “find some outlet coefficients.”
It is the much smaller question:

> Does the actual moving-throat core produce the normalized pair
> \[
> (\Xi_v,\Xi_T)\approx(-1.017,\;1.223)
> \]
> on the lower compensated branch, or does it remain on the naive equal-normalized source branch \((\Xi_T\approx0.927)\)?

## 4. Exact parent formulas for the corresponding raw parameters

If desired, these normalized values translate back to the raw parent parameters as
\[
\boxed{
\frac{v_{w0}}{c_s}
=
\frac{\hbar\sqrt{\mathcal Z_q}}{q_*\sqrt{\mu_0 m_\psi\rho_w}\;aL^{3/2}\ell^{3/2}}
\;\Xi_v,
}
\]
\[
\boxed{
\mathcal T_m
=
\frac{\sqrt{\mathcal Z_q m_\psi}}{a\sqrt{\mu_0\rho_wL\ell}}
\;\Xi_T.
}
\]

So once a concrete GNLS + localized-Maxwell throat-core solution supplies
\[
(\rho_w,\ell,\mathcal Z_q,q_*),
\]
the required background mixed flow and mouth traction are determined immediately.

## Result

The explicit core problem is now reduced to a very small finite target:

- one normalized background-flow number \(\Xi_v\),
- one normalized mouth-traction number \(\Xi_T\),
- and the choice between the naive equal-normalized source branch and the nearby lower compensated branch.

That is the smallest parent-level outlet-core target reached so far.
