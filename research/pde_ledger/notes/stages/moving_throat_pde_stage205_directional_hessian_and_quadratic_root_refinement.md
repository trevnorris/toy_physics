# Moving-Throat PDE — Stage 205: Directional Hessian, Exact Quadratic Log-Ray Refinement, and the Turning-Point Closure Test

## Status

**Exact within the carried Stage 203 scalar graph-slice theorem and the Stage 204 explicit free-quintuple log-ray family compiler.**

This stage does **not** introduce a new constitutive law.
It upgrades the Stage 204 first-order scalarized ray search to the exact **second-order** level, so the reduced branch search remains controlled even when the directional slope is small, vanishes, or is strongly curvature-corrected.

---

## Purpose

Stage 204 reduced the graph-lifted home-stretch problem to the scalar root equation
\[
\Phi_{\mathbf s}(\tau)=1,
\qquad
\Phi_{\mathbf s}(\tau):=\widehat\chi_Q(\mathbf y_{\mathbf s}(\tau)),
\]
and supplied two exact first-order local predictors:

- the affine predictor on `\(\widehat\chi_Q\)`, and
- the log-linear predictor on `\(\ln\widehat\chi_Q\)`.

That already solves the generic monotone case locally, but it leaves one practical gap:

> what is the exact second-order structure of the scalarized log-ray closure problem, and how do we treat rays whose first slope is weak or vanishes?

This stage answers that exactly.

The main outputs are:

1. the exact **directional Hessian operators** on free-quintuple log space,
2. the exact first/second derivative scalars
   \[
   \Phi_0,\ \Phi_1,\ \Phi_2,
   \qquad
   L_0,\ L_1,
   \]
   together with the exact identity
   \[
   \Phi_2=\Phi_0\,(L_1+L_0^2),
   \]
3. the exact **quadratic affine predictor** and **quadratic log predictor** for the closure point,
4. the exact discriminant conditions under which those predictors are real,
5. the exact **turning-point / tangency theorem** for rays with `\(\Phi_1=0\)`,
6. the curvature-corrected local expansions of the Stage 204 predictors,
7. and the exact theorem that the ordinary and logarithmic quadratic predictors agree through second order in the closure defect.

So Stage 205 turns the Stage 204 scalarized ray sweep into a true second-order spectral-placement tool.

---

## 1. Carry-forward scalarized graph-ray data

Keep the Stage 204 positive free quintuple
\[
\mathbf y=(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_W^{(\mathrm{eff})}),
\]
with logarithmic ray
\[
\mathbf y_{\mathbf s}(\tau)=\mathbf y_\circ\odot e^{\tau\mathbf s},
\qquad
\mathbf s=(s_\lambda,s_c,s_\gamma,s_U,s_W).
\]
The graph-lifted microscopic state remains
\[
\mathbf x_{\mathbf s}^{\rm graph}(\tau)=\mathbf x_*^{\rm graph}(\mathbf y_{\mathbf s}(\tau)),
\]
and Stage 204 already proved that this graph lift lies on the target orbit `\(\mathcal O_*\)` for **all** `\(\tau\)`.

So Packet B still vanishes identically:
\[
M_*\dot{\Delta\mathbf x}^{\rm graph}_{\mathbf s}=0.
\]
The only remaining reduced closure question is the Packet-A scalar slice
\[
\Phi_{\mathbf s}(\tau)=1,
\qquad
\Delta_{\mathbf s}(\tau):=\Phi_{\mathbf s}(\tau)-1.
\]

---

## 2. Exact directional Hessian on free log space

Introduce free log coordinates
\[
\boxed{
\boldsymbol\ell:=
(\ln\lambda_W,\ \ln c_{\eta U},\ \ln\gamma,\ \ln K_U,\ \ln K_W^{(\mathrm{eff})}).
}
\]
The Stage 204 directional derivative operator is
\[
\boxed{
\mathcal D_{\mathbf s}
:=
\sum_{i=1}^5 s_i\,\partial_{\ell_i}.
}
\]
The exact second-order directional operator is therefore
\[
\boxed{
\mathcal H_{\mathbf s}:=\mathcal D_{\mathbf s}^2
=\sum_{i,j=1}^5 s_i s_j\,\partial_{\ell_i}\partial_{\ell_j}.
}
\]
Equivalently, if
\[
H_{\widehat\chi_Q}(\mathbf y_\circ)
:=\bigl(\partial_{\ell_i}\partial_{\ell_j}\widehat\chi_Q\bigr)_{ij}\Big|_{\mathbf y_\circ},
\]
then the ordinary directional curvature is the quadratic form
\[
\boxed{
\Phi_2
:=
(\mathcal H_{\mathbf s}\widehat\chi_Q)(\mathbf y_\circ)
=
\mathbf s^T H_{\widehat\chi_Q}(\mathbf y_\circ)\,\mathbf s.
}
\]
Likewise, the logarithmic directional Hessian is
\[
H_{\ln\widehat\chi_Q}(\mathbf y_\circ)
:=\bigl(\partial_{\ell_i}\partial_{\ell_j}\ln\widehat\chi_Q\bigr)_{ij}\Big|_{\mathbf y_\circ},
\]
and its directional curvature is
\[
\boxed{
L_1
:=
(\mathcal H_{\mathbf s}\ln\widehat\chi_Q)(\mathbf y_\circ)
=
\mathbf s^T H_{\ln\widehat\chi_Q}(\mathbf y_\circ)\,\mathbf s.
}
\]

So once the Hessian matrices are known on the chosen free-quintuple base point, every ray curvature is obtained by a single quadratic form evaluation.

---

## 3. Exact first/second derivative scalars and ordinary/log identities

Define the base scalar values
\[
\boxed{
\Phi_0:=\Phi_{\mathbf s}(0)=\widehat\chi_Q(\mathbf y_\circ),
\qquad
\Phi_1:=\Phi_{\mathbf s}'(0)=(\mathcal D_{\mathbf s}\widehat\chi_Q)(\mathbf y_\circ),
\qquad
\Phi_2:=\Phi_{\mathbf s}''(0)=(\mathcal H_{\mathbf s}\widehat\chi_Q)(\mathbf y_\circ).
}
\]
Assume `\(\Phi_0>0\)`. Then the Stage 204 logarithmic slope is
\[
\boxed{
L_0:=\left.\frac{d}{d\tau}\ln\Phi_{\mathbf s}(\tau)\right|_{\tau=0}
=
\frac{\Phi_1}{\Phi_0}.
}
\]
The logarithmic directional curvature is
\[
\boxed{
L_1:=\left.\frac{d^2}{d\tau^2}\ln\Phi_{\mathbf s}(\tau)\right|_{\tau=0}
=
\frac{\Phi_2}{\Phi_0}-\frac{\Phi_1^2}{\Phi_0^2}.
}
\]
So the ordinary and logarithmic curvatures are related exactly by
\[
\boxed{
\Phi_2
=
\Phi_0\,(L_1+L_0^2).
}
\]
This is the central algebraic bridge of the stage.
It says the second-order Taylor data of `\(\widehat\chi_Q\)` and of `\(\ln\widehat\chi_Q\)` are not independent.

---

## 4. Exact quadratic scalar models on the ray

The ordinary second-order Taylor model of the closure residual is
\[
\boxed{
\Delta_{\mathbf s}^{\rm quad}(\tau)
:=
(\Phi_0-1)+\Phi_1\tau+\frac12\Phi_2\tau^2.
}
\]
The logarithmic second-order Taylor model is
\[
\boxed{
\mathcal L_{\mathbf s}^{\rm quad}(\tau)
:=
\ln\Phi_0+L_0\tau+\frac12 L_1\tau^2.
}
\]
The exact discriminants are therefore
\[
\boxed{
\Delta_{\rm aff}
:=
\Phi_1^2-2\Phi_2(\Phi_0-1),
}
\]
\[
\boxed{
\Delta_{\log}
:=
L_0^2-2L_1\ln\Phi_0.
}
\]
Whenever these are nonnegative, the corresponding quadratic predictor is real.

---

## 5. Exact quadratic affine predictor

Assume `\(\Phi_1\neq0\)` and choose the square-root branch continuously from the Stage 204 affine predictor. Then the exact quadratic predictor for the scalar slice is
\[
\boxed{
\tau_{\rm quad}
:=
\frac{2(1-\Phi_0)}{\Phi_1+\operatorname{sgn}(\Phi_1)\sqrt{\Delta_{\rm aff}}}.
}
\]
Equivalently, this is the root of
\[
(\Phi_0-1)+\Phi_1\tau+\frac12\Phi_2\tau^2=0
\]
that reduces continuously to the Stage 204 affine predictor when `\(\Phi_2\to0\)`.

### 5.1 Exact zero-curvature limit

As `\(\Phi_2\to0\)`,
\[
\boxed{
\tau_{\rm quad}\to\tau_{\rm aff}:=\frac{1-\Phi_0}{\Phi_1}.
}
\]
So Stage 204 is recovered exactly as the zero-curvature limit of the present construction.

---

## 6. Exact quadratic logarithmic predictor

Assume `\(L_0\neq0\)` and again choose the square-root branch continuously from the Stage 204 log-linear predictor. Then the exact logarithmic quadratic predictor is
\[
\boxed{
\tau_{\log,2}
:=
-\frac{2\ln\Phi_0}{L_0+\operatorname{sgn}(L_0)\sqrt{\Delta_{\log}}}.
}
\]
Equivalently, this is the root of
\[
\ln\Phi_0+L_0\tau+\frac12L_1\tau^2=0
\]
that reduces continuously to the Stage 204 log-linear predictor when `\(L_1\to0\)`.

### 6.1 Exact zero-curvature limit

As `\(L_1\to0\)`,
\[
\boxed{
\tau_{\log,2}\to\tau_{\log}:=-\frac{\ln\Phi_0}{L_0}.
}
\]
So the Stage 204 log-linear predictor is exactly the zero-log-curvature limit of the present construction.

---

## 7. Turning-point and tangency theorem

The new content of Stage 205 appears most sharply when the first directional slope vanishes.

### 7.1 Turning-point case: `\(\Phi_1=0\)` but `\(\Phi_0\neq1\)`

Then the quadratic scalar model becomes
\[
(\Phi_0-1)+\frac12\Phi_2\tau^2=0.
\]
So:

\[
\boxed{\textbf{Theorem (Stage 205 turning-point closure theorem).}}
\]

If `\(\Phi_1=0\)` and `\(\Phi_0\neq1\)`, then the quadratic model predicts real nearby closure points iff
\[
\boxed{
(1-\Phi_0)\Phi_2>0.
}
\]
In that case the two symmetric quadratic predictors are
\[
\boxed{
\tau_\pm
=
\pm\sqrt{\frac{2(1-\Phi_0)}{\Phi_2}}.
}
\]
If instead
\[
(1-\Phi_0)\Phi_2<0,
\]
then the quadratic model predicts **no nearby closure point** on that turning ray.

So Stage 205 provides the exact local replacement for the Stage 204 monotone theorem when the first slope vanishes.

### 7.2 Tangency case: `\(\Phi_0=1\)` and `\(\Phi_1=0\)`

If the base point already lies on the closure slice and the first derivative vanishes, then
\[
\boxed{
\Delta_{\mathbf s}(\tau)=\frac12\Phi_2\tau^2+O(\tau^3).
}
\]
So the ray is quadratically tangent to the closure slice at `\(\tau=0\)` unless `\(\Phi_2=0\)` as well, in which case higher-order data are needed.

This is the exact second-order tangency test for graph-lifted free-quintuple rays.

---

## 8. Curvature-corrected local expansions

Write the local closure defect as
\[
\Phi_0=1+\varepsilon,
\qquad |\varepsilon|\ll1.
\]
Then the Stage 204 and Stage 205 predictors are related by exact second-order expansions.

### 8.1 Ordinary quadratic correction

The affine predictor is
\[
\tau_{\rm aff}=\frac{1-\Phi_0}{\Phi_1}.
\]
The quadratic predictor satisfies
\[
\boxed{
\tau_{\rm quad}
=
\tau_{\rm aff}-\frac{\Phi_2}{2\Phi_1}\,\tau_{\rm aff}^2+O(\tau_{\rm aff}^3).
}
\]
So the ordinary second-order correction is controlled by the dimensionless ratio `\(\Phi_2/\Phi_1\)`.

### 8.2 Logarithmic quadratic correction

The Stage 204 log-linear predictor is
\[
\tau_{\log}=-\frac{\ln\Phi_0}{L_0}.
\]
The quadratic log predictor satisfies
\[
\boxed{
\tau_{\log,2}
=
\tau_{\log}-\frac{L_1}{2L_0}\,\tau_{\log}^2+O(\tau_{\log}^3).
}
\]
So the logarithmic second-order correction is controlled by the dimensionless ratio `\(L_1/L_0\)`.

### 8.3 Agreement theorem for the two quadratic predictors

The ordinary and logarithmic quadratic predictors agree more strongly than the Stage 204 first-order pair. In fact,
\[
\boxed{
\tau_{\log,2}-\tau_{\rm quad}=O((\Phi_0-1)^3).
}
\]
More explicitly, with `\(\Phi_0=1+\varepsilon\)`,
\[
\boxed{
\tau_{\log,2}-\tau_{\rm quad}
=
\frac{(L_0^2+3L_1)}{6L_0^3}\,\varepsilon^3+O(\varepsilon^4).
}
\]
So the second-order ordinary and logarithmic predictors coincide through quadratic order in the closure defect. Only at cubic order do they begin to separate.

This is the sharpest local predictor-equivalence statement reached so far in the free-quintuple ray program.

---

## 9. Best current reading after Stage 205

Stage 204 reduced the reduced closure search to a scalarized free-quintuple log-ray root problem.
Stage 205 now gives the exact second-order refinement of that problem.

The Packet-B orbit defect still vanishes identically on every graph-lifted ray, so the only remaining reduced question is the scalar Packet-A slice. But that scalar search is now controlled by:

1. the Stage 204 first derivative `\(L_0\)` or `\(\Phi_1\)`,
2. the Stage 205 directional Hessian `\(L_1\)` or `\(\Phi_2\)`,
3. the quadratic discriminants `\(\Delta_{\rm aff},\Delta_{\log}\)`,
4. and, when needed, the exact turning-point criterion.

So after Stage 205 the next honest continuation is no longer to invent a better local predictor by hand.
It is to evaluate the actual Hessian data of `\(\widehat\chi_Q\)` on candidate base points and free directions, and then use the quadratic crossing / tangency theorem to rank which free-quintuple rays are genuinely promising for the final branch search.
