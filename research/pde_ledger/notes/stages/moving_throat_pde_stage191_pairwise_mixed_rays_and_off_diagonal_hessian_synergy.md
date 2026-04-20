# Moving-Throat PDE — Stage 191: Pairwise Mixed Rays, Exact Off-Diagonal Hessian Synergy, and the Canonical Two-Ray Audit

## Status

**Exact within the carried Stage-189 scalarized search framework and the Stage-190 primitive certified sieve, once local pairwise diagonal/off-diagonal Hessian envelopes are supplied on the chosen mixed-ray intervals.**

This stage does **not** introduce a new constitutive law.
It is the first exact continuation beyond the primitive-ray table: the first place where the off-diagonal Hessian entries of
\[
\ln\widehat\chi_Q
\]
can actually enter the certified local search.

---

## Purpose

Stage 190 proved that every primitive ray uses only one diagonal Hessian restriction
\[
\partial_{\ell_i\ell_i}\ln\widehat\chi_Q,
\]
so no off-diagonal Hessian entry can affect the primitive certified brackets.
That already closed the primitive local sieve, but it also made the next gate completely sharp:

> what is the first exact mixed-ray family, how do the off-diagonal Hessian terms enter it, and what is the smallest certified pairwise audit that can possibly beat the primitive winners?

This stage answers that.

The main outputs are:

1. the exact **pairwise mixed-ray cone** built from two canonically oriented primitive descent directions,
2. the exact **gradient-synergy theorem** for the initial mixed slope,
3. the exact **off-diagonal curvature-synergy law** and the proof that the equal-mix ray maximizes cross-Hessian leverage,
4. the exact **mixed-ray curvature-envelope theorem** and certified bracket theorem for a fixed mixing ratio,
5. two exact canonical screen rays for each primitive pair:
   - the **gradient-optimal** ray,
   - the **equal-mix synergy** ray,
6. and the exact rule that the full one-parameter pairwise optimizer is now the only missing step beyond this first mixed-ray screen.

So Stage 191 is the first theorem-level audit of genuine two-coordinate mixed rays.
It does **not** yet solve the full pairwise optimization problem, but it reduces that problem to a one-parameter family with two exact canonical screening rays.

---

## 1. Carry-forward primitive data and the pairwise mixed-ray cone

Keep the free log coordinates
\[
\boldsymbol\ell=(\ell_\lambda,\ell_c,\ell_\gamma,\ell_U,\ell_W)
\]
and the scalar closure function
\[
h(\boldsymbol\ell)=\ln\widehat\chi_Q(\boldsymbol\ell).
\]
At the base point `\(\boldsymbol\ell_\circ\)`, Stage 190 already defined
\[
H_0=|h(\boldsymbol\ell_\circ)|>0,
\qquad
\Gamma_i=\sigma_0\,\partial_{\ell_i}h(\boldsymbol\ell_\circ),
\qquad
\sigma_0=\operatorname{sgn}(h(\boldsymbol\ell_\circ)).
\]
For every primitive axis with `\(\Gamma_i\neq0\)`, the canonical oriented primitive direction is
\[
\widehat{\mathbf e}_i:=\varepsilon_i\mathbf e_i,
\qquad
\varepsilon_i:=-\operatorname{sgn}(\Gamma_i),
\qquad
k_i:=|\Gamma_i|>0.
\]
So the forward primitive slope on `\(\widehat{\mathbf e}_i\)` is exactly
\[
K_i=-k_i<0.
\]

Now fix a **monotone–monotone primitive pair** `\(i\neq j\)` with
\[
\Gamma_i\neq0,
\qquad
\Gamma_j\neq0.
\]
The first genuine mixed-ray cone is the one-parameter family
\[
\boxed{
\widehat{\mathbf s}_{ij}(r)
:=
\frac{\widehat{\mathbf e}_i+r\widehat{\mathbf e}_j}{\sqrt{1+r^2}},
\qquad r\ge0.
}
\]
If all five primitive gradients are nonzero, this produces the full set of
\[
\binom{5}{2}=10
\]
monotone pairwise cones:
\[
(\lambda,c),\ (\lambda,\gamma),\ (\lambda,U),\ (\lambda,W),\ (c,\gamma),\ (c,U),\ (c,W),\ (\gamma,U),\ (\gamma,W),\ (U,W).
\]
Mixed rays involving a turning primitive axis are deferred to the next stage.

---

## 2. Exact gradient-synergy theorem

The oriented initial slope on the pairwise mixed ray is
\[
\boxed{
K_{ij}(r)
:=
\sigma_0\,\nabla_\ell h(\boldsymbol\ell_\circ)\cdot \widehat{\mathbf s}_{ij}(r)
=
-\frac{k_i+r k_j}{\sqrt{1+r^2}}.
}
\]
So the positive slope magnitude carried into the Stage-189 root map is
\[
\boxed{
k_{ij}(r):=\frac{k_i+r k_j}{\sqrt{1+r^2}}>0.
}
\]

Differentiating gives the exact first-order mixed-slope law
\[
\boxed{
\frac{d k_{ij}}{dr}
=
\frac{k_j-k_i r}{(1+r^2)^{3/2}}.
}
\]
Therefore the mixed-slope magnitude has a unique maximizer
\[
\boxed{
r_{ij}^{\rm grad}=\frac{k_j}{k_i},
}
\]
with exact maximum value
\[
\boxed{
k_{ij}^{\rm grad}=\sqrt{k_i^2+k_j^2}.
}
\]
So the first mixed-ray theorem is:

### Exact gradient-synergy theorem

For every monotone–monotone primitive pair,
\[
\boxed{
\max_{r\ge0} k_{ij}(r)=\sqrt{k_i^2+k_j^2}>\max(k_i,k_j).
}
\]
Thus a mixed ray can improve the **first-order** descent rate even when the off-diagonal Hessian entry vanishes.
The off-diagonal Hessian is therefore the first genuinely new **curvature** datum, not the first source of mixed-ray advantage.

---

## 3. Exact off-diagonal curvature-synergy law

Along the pairwise mixed ray define the oriented Hessian entries
\[
h_{ii}^{(ij)}(r;\tau)
:=
\sigma_0\,\partial_{\ell_i\ell_i}h\bigl(\boldsymbol\ell_\circ+\tau\widehat{\mathbf s}_{ij}(r)\bigr),
\]
\[
h_{jj}^{(ij)}(r;\tau)
:=
\sigma_0\,\partial_{\ell_j\ell_j}h\bigl(\boldsymbol\ell_\circ+\tau\widehat{\mathbf s}_{ij}(r)\bigr),
\]
\[
h_{ij}^{(ij)}(r;\tau)
:=
\sigma_0\,\varepsilon_i\varepsilon_j\,
\partial_{\ell_i\ell_j}h\bigl(\boldsymbol\ell_\circ+\tau\widehat{\mathbf s}_{ij}(r)\bigr).
\]
Then the exact mixed second directional derivative is
\[
\boxed{
H_{1,ij}(r;\tau)
=
\frac{h_{ii}^{(ij)}(r;\tau)+2r\,h_{ij}^{(ij)}(r;\tau)+r^2 h_{jj}^{(ij)}(r;\tau)}{1+r^2}.
}
\]
Equivalently,
\[
\boxed{
H_{1,ij}(r;\tau)
=
\frac{h_{ii}^{(ij)}+r^2 h_{jj}^{(ij)}}{1+r^2}
+
\frac{2r}{1+r^2}\,h_{ij}^{(ij)}.
}
\]
So the off-diagonal Hessian contributes with the exact cross weight
\[
\boxed{
w_{\times}(r):=\frac{2r}{1+r^2}.
}
\]
It satisfies
\[
0\le w_{\times}(r)\le1,
\qquad
\frac{d w_{\times}}{dr}=\frac{2(1-r^2)}{(1+r^2)^2},
\]
hence
\[
\boxed{
\max_{r\ge0} w_{\times}(r)=1
\quad\text{at}\quad r=1.
}
\]

This gives the second main theorem.

### Exact off-diagonal curvature-synergy theorem

The **equal-mix ray**
\[
\boxed{
\widehat{\mathbf s}_{ij}^{\rm eq}:=\frac{\widehat{\mathbf e}_i+\widehat{\mathbf e}_j}{\sqrt2}
}
\]
maximizes the absolute leverage of the off-diagonal Hessian entry inside the pairwise cone.

So there are two exact but different canonical ratios already at Stage 191:

- `\(r_{ij}^{\rm grad}=k_j/k_i\)` maximizes first-order descent,
- `\(r_{ij}^{\rm eq}=1\)` maximizes cross-Hessian leverage.

### Exact diagonal-neutrality law

If
\[
h_{ij}^{(ij)}(r;\tau)=0,
\]
then the pairwise curvature is the convex interpolation of the two diagonal restrictions:
\[
\boxed{
H_{1,ij}(r;\tau)=\frac{h_{ii}^{(ij)}(r;\tau)+r^2 h_{jj}^{(ij)}(r;\tau)}{1+r^2}.
}
\]
So the off-diagonal Hessian is exactly the first datum that can lower or raise the pairwise curvature beyond diagonal interpolation.

---

## 4. Exact mixed-ray curvature envelopes and certified brackets

Fix a pair `\((i,j)\)` and a chosen mixing ratio `\(r\ge0\)`.
Assume the pairwise mixed ray is valid on an interval
\[
\boxed{0\le \tau\le T_{ij}(r)}
\]
and that on this interval the three oriented Hessian entries admit bounds
\[
\underline h_{ii}\le h_{ii}^{(ij)}(r;\tau)\le \overline h_{ii},
\qquad
\underline h_{ij}\le h_{ij}^{(ij)}(r;\tau)\le \overline h_{ij},
\qquad
\underline h_{jj}\le h_{jj}^{(ij)}(r;\tau)\le \overline h_{jj}.
\]
Because the weights are nonnegative for `\(r\ge0\)`, the exact mixed curvature bounds are
\[
\boxed{
\underline\kappa_{ij}(r)
:=
\frac{\underline h_{ii}+2r\underline h_{ij}+r^2\underline h_{jj}}{1+r^2},
}
\]
\[
\boxed{
\overline\kappa_{ij}(r)
:=
\frac{\overline h_{ii}+2r\overline h_{ij}+r^2\overline h_{jj}}{1+r^2}.
}
\]
So the pairwise mixed ray is again inside the exact Stage-189 certified framework, now with slope magnitude `\(k_{ij}(r)\)` and curvature envelope `\([\underline\kappa_{ij}(r),\overline\kappa_{ij}(r)]\)`.

Using the carried monotone-ray root map
\[
\boxed{
\mathcal T(H_0,k;c):=
\frac{2H_0}{k+\sqrt{k^2-2cH_0}},
\qquad k>0,
}
\]
define the mixed discriminants
\[
\boxed{
\Delta_{ij}^{\rm lo}(r):=k_{ij}(r)^2-2\underline\kappa_{ij}(r)H_0,
\qquad
\Delta_{ij}^{\rm hi}(r):=k_{ij}(r)^2-2\overline\kappa_{ij}(r)H_0.
}
\]
Whenever both are nonnegative, the exact pairwise certified bracket is
\[
\boxed{
\tau_{ij,\rm lo}(r):=\mathcal T\bigl(H_0,k_{ij}(r);\underline\kappa_{ij}(r)\bigr),
\qquad
\tau_{ij,\rm hi}(r):=\mathcal T\bigl(H_0,k_{ij}(r);\overline\kappa_{ij}(r)\bigr).
}
\]
If also
\[
\boxed{
\tau_{ij,\rm hi}(r)\le T_{ij}(r),
}
\]
then there exists one unique true mixed-ray closure point on that ray, lying in
\[
\boxed{
\tau_*^{(ij)}(r)\in[\tau_{ij,\rm lo}(r),\tau_{ij,\rm hi}(r)].
}
\]
So Stage 191 reduces every fixed pairwise ratio to the same certified bracket theorem already used for primitive rows.

---

## 5. The two canonical pairwise audit rays

Although the full pairwise cone is one-parameter, two exact rays already isolate the two independent advantages a mixed ray can exploit.

### 5.1 Gradient-optimal ray

The exact gradient-optimal ratio is
\[
\boxed{
r_{ij}^{\rm grad}=\frac{k_j}{k_i}.}
\]
So the normalized gradient-optimal ray is
\[
\boxed{
\widehat{\mathbf s}_{ij}^{\rm grad}
=
\frac{k_i\widehat{\mathbf e}_i+k_j\widehat{\mathbf e}_j}{\sqrt{k_i^2+k_j^2}}.
}
\]
Its exact slope magnitude is
\[
\boxed{
k_{ij}^{\rm grad}=\sqrt{k_i^2+k_j^2}.
}
\]
Its exact curvature is the weighted Rayleigh quotient
\[
\boxed{
H_{1,ij}^{\rm grad}(\tau)
=
\frac{k_i^2 h_{ii}^{\rm grad}(\tau)+2k_i k_j h_{ij}^{\rm grad}(\tau)+k_j^2 h_{jj}^{\rm grad}(\tau)}{k_i^2+k_j^2}.
}
\]
So this ray maximizes first-order descent inside the pairwise cone.

### 5.2 Equal-mix synergy ray

The exact equal-mix ray is
\[
\boxed{
\widehat{\mathbf s}_{ij}^{\rm eq}
=
\frac{\widehat{\mathbf e}_i+\widehat{\mathbf e}_j}{\sqrt2}.
}
\]
Its exact slope magnitude is
\[
\boxed{
k_{ij}^{\rm eq}=\frac{k_i+k_j}{\sqrt2}.
}
\]
Its exact curvature is
\[
\boxed{
H_{1,ij}^{\rm eq}(\tau)
=
\frac{h_{ii}^{\rm eq}(\tau)+2h_{ij}^{\rm eq}(\tau)+h_{jj}^{\rm eq}(\tau)}{2}.
}
\]
So this ray maximizes the cross-Hessian leverage inside the pairwise cone.

### 5.3 Exact comparison of the two canonical rays

The gradient-optimal and equal-mix rays coincide iff
\[
\boxed{k_i=k_j.}
\]
Otherwise they are distinct.
So Stage 191 exposes the first genuine tradeoff in the mixed-ray search:

- one ray optimizes **descent magnitude**,
- the other optimizes **off-diagonal curvature leverage**.

That is why both belong in the canonical pairwise audit.

---

## 6. Exact canonical pairwise screen and the promotion rule

For every monotone primitive pair `\((i,j)\)`, define the two canonical mixed-ray rows
\[
\boxed{
\mathcal R_{ij}^{\rm grad}
:=
(H_0,k_{ij}^{\rm grad},\underline\kappa_{ij}^{\rm grad},\overline\kappa_{ij}^{\rm grad},T_{ij}^{\rm grad}),
}
\]
\[
\boxed{
\mathcal R_{ij}^{\rm eq}
:=
(H_0,k_{ij}^{\rm eq},\underline\kappa_{ij}^{\rm eq},\overline\kappa_{ij}^{\rm eq},T_{ij}^{\rm eq}).
}
\]
Here the envelope data are obtained by evaluating the Stage-191 mixed curvature-envelope formulas at the corresponding canonical ratio.

This gives the first exact pairwise promotion rule.

### Canonical pairwise promotion theorem

If one of the two canonical rows is admissible and its certified upper bracket lies strictly to the left of every surviving primitive lower bracket, i.e.
\[
\boxed{
\tau_{ij,\rm hi}^{\rm grad}<\min_m \tau_{m,\rm lo}^{\rm prim}
\quad\text{or}\quad
\tau_{ij,\rm hi}^{\rm eq}<\min_m \tau_{m,\rm lo}^{\rm prim},
}
\]
then that pair is already certified to beat the primitive sieve.

This theorem is exact because the corresponding canonical mixed ray is itself a valid candidate ray, not a bound.

### Canonical pairwise deferral rule

If both canonical rows fail admissibility or fail to beat the primitive table, the pair is **not eliminated**. It is merely deferred to the full one-parameter pairwise optimizer.

So Stage 191 is a true first mixed-ray audit, but it is not yet the final pairwise search theorem.

---

## 7. Minimal data packet for the full pairwise optimizer

The full pairwise optimizer of the next stage will need, for each monotone pair `\((i,j)\)`, the exact finite packet
\[
\boxed{
\mathcal P_{ij}^{\rm mix}
:=
\bigl(k_i,k_j,
[\underline h_{ii},\overline h_{ii}],
[\underline h_{ij},\overline h_{ij}],
[\underline h_{jj},\overline h_{jj}],
T_{ij}(r)\bigr),
}
\]
with the first five entries evaluated at the base point and the interval information supplied on the chosen pairwise cone.

Everything else is downstream algebra:

1. slope magnitude `\(k_{ij}(r)\)`,
2. curvature envelopes `\(\underline\kappa_{ij}(r),\overline\kappa_{ij}(r)\)`,
3. discriminants,
4. certified brackets,
5. and the pairwise winner comparison.

So Stage 191 turns the mixed-ray problem into a one-parameter certified search family, no longer an unconstrained five-dimensional speculation.

---

## 8. Best current reading after Stage 191

Stage 190 closed the primitive certified ray table and proved that off-diagonal Hessian data are invisible there.
Stage 191 now gives the first exact mixed-ray continuation:

1. a monotone–monotone primitive pair generates a one-parameter mixed-ray cone,
2. the first-order descent rate is optimized at
   \[
   r_{ij}^{\rm grad}=k_j/k_i,
   \]
3. the off-diagonal Hessian leverage is optimized at
   \[
   r_{ij}^{\rm eq}=1,
   \]
4. every fixed ratio still lies inside the exact Stage-189 certified bracket framework,
5. and the first honest mixed-ray audit is therefore the two-row canonical screen
   \[
   \mathcal R_{ij}^{\rm grad},\qquad \mathcal R_{ij}^{\rm eq}.
   \]

So the next honest continuation is now very sharply defined:

> carry out the full one-parameter optimization over `\(r\)` on each surviving pairwise cone, and decide whether any mixed ray can beat both the primitive and canonical-screen rows.

That is the natural Stage-192 theorem gate.
