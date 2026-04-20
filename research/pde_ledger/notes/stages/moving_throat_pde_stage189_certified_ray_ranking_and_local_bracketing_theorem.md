# Moving-Throat PDE — Stage 189: Curvature-Enveloped Ray Ranking, Certified Local Brackets, and the Search Sieve Theorem

## Status

**Exact within the carried Stage-188 scalarized log-ray / directional-Hessian framework, once a local curvature envelope is supplied on the chosen oriented ray.**

This stage does **not** introduce a new constitutive law.
It upgrades the Stage-188 quadratic predictor from a local approximation into a **certified local search tool** by adding exact curvature majorants/minorants for the oriented logarithmic closure residual.

---

## Purpose

Stage 188 gave the exact first/second derivative data
\[
\Phi_0,\ \Phi_1,\ \Phi_2,
\qquad
L_0,\ L_1,
\]
and the exact quadratic affine/log predictors for the scalarized free-quintuple ray problem
\[
\Phi_{\mathbf s}(\tau)=1.
\]
That already solved the local predictor problem, but it still left one practical gap:

> how do we turn those second-order predictors into a **certified** local root bracket, and how do we rank competing rays without running the full scalar solve on all of them?

This stage answers that exactly.

The main outputs are:

1. the exact **oriented logarithmic residual** and its local curvature envelope,
2. exact quadratic **minorant / majorant** functions on a ray,
3. the exact monotone-branch root map
   \[
   \mathcal T(H_0,K_0;c),
   \]
4. the exact **certified bracket theorem** for monotone rays,
5. the exact **turning-ray bracket theorem** when the first slope vanishes,
6. the exact bracket-width law and its small-envelope expansion,
7. the exact **pairwise ray-ordering theorem** for disjoint certified brackets,
8. and the final Stage-189 conclusion that the free-quintuple search has become a true **search sieve** rather than only a local predictor library.

So Stage 189 is the first place where the Stage-187/188 scalarized ray program becomes a rigorous branch-ranking tool.

---

## 1. Carry-forward scalarized ray and the oriented logarithmic residual

Keep the Stage-187 graph-lifted free-quintuple log ray
\[
\mathbf y_{\mathbf s}(\tau)=\mathbf y_\circ\odot e^{\tau\mathbf s},
\qquad
\Phi_{\mathbf s}(\tau):=\widehat\chi_Q(\mathbf y_{\mathbf s}(\tau)).
\]
Stage 188 already defined the logarithmic scalar
\[
h_{\mathbf s}(\tau):=\ln \Phi_{\mathbf s}(\tau).
\]
To make the local crossing problem one-sided, define the **oriented** logarithmic residual
\[
\boxed{
H_{\mathbf s}(\tau)
:=
\operatorname{sgn}\!\bigl(h_{\mathbf s}(0)\bigr)\,h_{\mathbf s}(\tau).
}
\]
If necessary, replace `\(\mathbf s\)` by `\(-\mathbf s\)` so that the forward search direction is toward the closure slice. Then the oriented local data are
\[
\boxed{
H_0:=H_{\mathbf s}(0)=\bigl|\ln\Phi_{\mathbf s}(0)\bigr|>0,
\qquad
K_0:=H_{\mathbf s}'(0)<0.
}
\]
The remaining closure condition is simply
\[
\boxed{
H_{\mathbf s}(\tau)=0.
}
\]
So every admissible local branch search can be normalized to the same canonical shape:

- positive initial oriented defect `\(H_0>0\)`,
- negative initial oriented slope `\(K_0<0\)`,
- and forward search parameter `\(\tau\ge0\)`.

This is the natural Stage-188 continuation because it removes the sign bookkeeping while keeping the exact logarithmic scalarization.

---

## 2. Exact local curvature envelope and quadratic comparison functions

Assume there is a controlled interval
\[
\boxed{I_T:=[0,T]}
\]
on which the oriented log curvature is bounded by exact constants
\[
\boxed{
\underline K_1
\le
H_{\mathbf s}''(\tau)
\le
\overline K_1
\qquad
(0\le \tau\le T).
}
\]
Integrating twice from `\(0\)` gives the exact quadratic comparison functions
\[
\boxed{
H_-(\tau):=H_0+K_0\tau+\frac12\underline K_1\tau^2,
\qquad
H_+(\tau):=H_0+K_0\tau+\frac12\overline K_1\tau^2,
}
\]
with the exact sandwich theorem
\[
\boxed{
H_-(\tau)
\le
H_{\mathbf s}(\tau)
\le
H_+(\tau)
\qquad
(0\le\tau\le T).
}
\]
So the Stage-188 second-order data become rigorous as soon as the directional Hessian of `\(\ln\widehat\chi_Q\)` is bounded on the local interval.

### 2.1 Exact discriminants

Define the quadratic discriminants
\[
\boxed{
\Delta(c):=K_0^2-2cH_0,
\qquad
\Delta_-:=\Delta(\underline K_1),
\qquad
\Delta_+:=\Delta(\overline K_1).
}
\]
Whenever `\(\Delta(c)\ge0\)`, the quadratic comparison curve with curvature `\(c\)` has a real forward root.

---

## 3. Exact monotone-branch root map

For any real curvature parameter `\(c\)` with `\(\Delta(c)\ge0\)`, define the forward quadratic root map
\[
\boxed{
\mathcal T(H_0,K_0;c)
:=
-\frac{2H_0}{K_0+\operatorname{sgn}(K_0)\sqrt{\Delta(c)}}.
}
\]
On the oriented search branch `\(K_0<0\)`, this becomes
\[
\boxed{
\mathcal T(H_0,K_0;c)
=
-\frac{2H_0}{K_0-\sqrt{K_0^2-2cH_0}}.
}
\]
It is the exact forward root of
\[
H_0+K_0\tau+\frac12 c\tau^2=0
\]
that reduces continuously to the Stage-187/188 zero-curvature predictor.

### 3.1 Zero-curvature limit

As `\(c\to0\)`,
\[
\boxed{
\mathcal T(H_0,K_0;0)=-\frac{H_0}{K_0}.
}
\]
So the Stage-187 log-linear predictor is the zero-curvature limit of the Stage-189 certified bracket map.

### 3.2 Exact monotonicity in the curvature parameter

Differentiate the defining quadratic relation implicitly. One gets
\[
\boxed{
\frac{\partial\mathcal T}{\partial c}
=
\frac{\mathcal T(H_0,K_0;c)^2}{2\sqrt{K_0^2-2cH_0}}
>0.
}
\]
So the forward root grows strictly with the curvature parameter.
This is the key exact ordering fact behind the certified bracket theorem below.

### 3.3 Collapse to the Stage-188 quadratic predictor

If the curvature envelope collapses to a point,
\[
\underline K_1=\overline K_1=L_1,
\]
then
\[
\boxed{
\tau_{\rm lo}=\tau_{\rm hi}=\mathcal T(H_0,K_0;L_1),
}
\]
which is exactly the Stage-188 quadratic logarithmic predictor written in the oriented variables.

So Stage 188 is recovered as the zero-width Stage-189 bracket.

---

## 4. Exact certified bracket theorem for monotone rays

Define the two quadratic comparison roots
\[
\boxed{
\tau_{\rm lo}:=\mathcal T(H_0,K_0;\underline K_1),
\qquad
\tau_{\rm hi}:=\mathcal T(H_0,K_0;\overline K_1).
}
\]
By the exact monotonicity of `\(\mathcal T\)`,
\[
\boxed{
0<\tau_{\rm lo}\le \tau_{\rm hi}
}
\]
whenever `\(\Delta_\pm\ge0\)`.

### 4.1 Exact descent sign on the whole bracket interval

On `\([0,\tau_{\rm hi}]\)` the derivative obeys
\[
H_{\mathbf s}'(\tau)
\le
K_0+\overline K_1\tau.
\]
If `\(\overline K_1\ge0\)`, the right-hand side is increasing in `\(\tau\)` and therefore
\[
K_0+\overline K_1\tau
\le
K_0+\overline K_1\tau_{\rm hi}
=
-\sqrt{\Delta_+}<0.
\]
If `\(\overline K_1<0\)`, the right-hand side is decreasing and already satisfies
\[
K_0+\overline K_1\tau \le K_0 <0.
\]
So in both cases
\[
\boxed{
H_{\mathbf s}'(\tau)<0
\qquad(0\le\tau\le\tau_{\rm hi}).
}
\]
Thus the true oriented residual is strictly decreasing throughout the certified bracket interval.

### 4.2 Certified local bracket theorem

\[
\boxed{\textbf{Theorem (Stage 189 certified monotone bracket theorem).}}
\]

Assume:

1. `\(H_0>0\)` and `\(K_0<0\)` on the oriented ray,
2. `\(\underline K_1\le H_{\mathbf s}''(\tau)\le\overline K_1\)` on `\([0,T]\)`,
3. `\(\Delta_-\ge0\)` and `\(\Delta_+\ge0\)`,
4. `\(\tau_{\rm hi}\le T\)`.

Then:

- for every `\(0\le\tau<\tau_{\rm lo}\)`, one has
  \[
  H_{\mathbf s}(\tau)\ge H_-(\tau)>0,
  \]
  so the true root cannot occur before `\(\tau_{\rm lo}\)`;
- at the upper comparison root,
  \[
  H_{\mathbf s}(\tau_{\rm hi})\le H_+(\tau_{\rm hi})=0,
  \]
  so the true root has occurred by `\(\tau_{\rm hi}\)`;
- and because `\(H_{\mathbf s}'(\tau)<0\)` on `\([0,\tau_{\rm hi}]\)`, the true root is unique.

Therefore there exists a unique actual closure point
\[
\boxed{
\tau_*\in[\tau_{\rm lo},\tau_{\rm hi}].
}
\]

So Stage 189 upgrades the Stage-188 quadratic predictor into a rigorous local bracket.

---

## 5. Exact turning-ray bracket theorem

Stage 188 already identified the turning-point criterion for rays with vanishing first slope.
Stage 189 promotes that criterion to a certified local bracket.

Assume the oriented initial slope vanishes,
\[
\boxed{K_0=0,\qquad H_0>0,}
\]
and the local curvature envelope is strictly negative:
\[
\boxed{
\underline K_1
\le
H_{\mathbf s}''(\tau)
\le
\overline K_1<0
\qquad(0\le\tau\le T).
}
\]
Then the quadratic comparison functions are
\[
H_-(\tau)=H_0+\frac12\underline K_1\tau^2,
\qquad
H_+(\tau)=H_0+\frac12\overline K_1\tau^2.
\]
Define the turning-ray comparison roots
\[
\boxed{
\tau_{\rm lo}^{\rm(tp)}:=\sqrt{-\frac{2H_0}{\underline K_1}},
\qquad
\tau_{\rm hi}^{\rm(tp)}:=\sqrt{-\frac{2H_0}{\overline K_1}}.
}
\]
Because `\(\underline K_1\le\overline K_1<0\)`, one has
\[
\boxed{
0<\tau_{\rm lo}^{\rm(tp)}\le \tau_{\rm hi}^{\rm(tp)}.
}
\]
Also
\[
H_{\mathbf s}'(\tau)
\le
\overline K_1\tau<0
\qquad(0<\tau\le\tau_{\rm hi}^{\rm(tp)}),
\]
so the turning ray becomes strictly decreasing immediately after the base point.

Hence:
\[
\boxed{\textbf{Theorem (Stage 189 certified turning-ray bracket theorem).}}
\]

If `\(\tau_{\rm hi}^{\rm(tp)}\le T\)`, then the true turning-ray closure point exists uniquely and obeys
\[
\boxed{
\tau_*\in[\tau_{\rm lo}^{\rm(tp)},\tau_{\rm hi}^{\rm(tp)}].
}
\]

So even the Stage-188 tangency/turning case can be turned into a certified bracket whenever the oriented curvature is strictly negative on a controlled interval.

---

## 6. Exact bracket width law and the small-envelope expansion

For a monotone admissible ray define the exact bracket width
\[
\boxed{
W_{\mathbf s}:=\tau_{\rm hi}-\tau_{\rm lo}.
}
\]
Because the root map is strictly increasing in the curvature parameter,
\[
\boxed{
W_{\mathbf s}>0
\iff
\overline K_1>\underline K_1.
}
\]
So the bracket width is a direct measure of unresolved curvature uncertainty on the chosen ray.

### 6.1 Symmetric-envelope expansion

Write the envelope as
\[
\underline K_1=\bar K_1-\frac12\Delta K_1,
\qquad
\overline K_1=\bar K_1+\frac12\Delta K_1.
\]
Then, by odd-order Taylor cancellation,
\[
\boxed{
W_{\mathbf s}
=
\frac{\mathcal T(H_0,K_0;\bar K_1)^2}{2\sqrt{K_0^2-2\bar K_1H_0}}
\,\Delta K_1
+O(\Delta K_1^3).
}
\]
So the width is linear in the curvature-envelope size at leading order.

### 6.2 Zero-curvature simplification

At `\(\bar K_1=0\)`, write
\[
\tau_0:=-\frac{H_0}{K_0},
\]
which is the oriented Stage-187 log-linear predictor. Then
\[
\boxed{
W_{\mathbf s}
=
\frac{\tau_0^2}{2|K_0|}
(\overline K_1-\underline K_1)
+O\!\bigl((\overline K_1-\underline K_1)^3\bigr).
}
\]
So, near the Stage-187 limit, rays with larger initial descent `\(|K_0|\)` and smaller curvature-envelope size are exactly the rays with the tightest certified brackets.

This is the first exact Stage-187/188 search-quality law.

---

## 7. Exact pairwise ray-ordering theorem

Let two admissible rays `\(\mathbf s_a\)` and `\(\mathbf s_b\)` have certified brackets
\[
[\tau_{\rm lo}^{(a)},\tau_{\rm hi}^{(a)}],
\qquad
[\tau_{\rm lo}^{(b)},\tau_{\rm hi}^{(b)}].
\]
Then:
\[
\boxed{\textbf{Theorem (Stage 189 certified pairwise ray-ordering theorem).}}
\]

If
\[
\boxed{
\tau_{\rm hi}^{(a)}<\tau_{\rm lo}^{(b)},
}
\]
then the actual closure points satisfy
\[
\boxed{
\tau_*^{(a)}<\tau_*^{(b)}.
}
\]
So disjoint certified brackets produce a strict and exact ordering of rays.

This is the first rigorous ranking theorem in the free-quintuple search program.

### 7.1 Canonical local priority pair

When certified brackets overlap, strict theorem-level ordering is not yet available. The natural audit ordering is then the exact priority pair
\[
\boxed{
\mathcal Q_{\mathbf s}:=(\tau_{\rm hi},\ W_{\mathbf s}).
}
\]
Smaller `\(\tau_{\rm hi}\)` means the ray reaches closure no later than a smaller guaranteed time, while smaller `\(W_{\mathbf s}\)` means the local placement is more sharply certified.

So the Stage-189 search sieve is:

1. reject rays with no admissible local bracket,
2. certify strict pairwise order whenever brackets are disjoint,
3. among overlapping candidates, sort lexicographically by `\(\mathcal Q_{\mathbf s}\)`.

This last step is an audit convention rather than a theorem, but it is now backed by exact bracket data rather than by heuristic slope comparisons alone.

---

## 8. Exact admissibility test for the local search sieve

A free-quintuple ray enters the Stage-189 search sieve only if one of the following exact local conditions is met.

### 8.1 Monotone admissible ray

A ray is locally admissible in the monotone class if:

1. the oriented base data satisfy `\(H_0>0\)` and `\(K_0<0\)`,
2. a finite curvature envelope `\([\underline K_1,\overline K_1]\)` is known on `\([0,T]\)`,
3. both quadratic discriminants are nonnegative,
4. the certified upper bracket satisfies `\(\tau_{\rm hi}\le T\)`.

### 8.2 Turning admissible ray

A ray is locally admissible in the turning class if:

1. `\(H_0>0\)` and `\(K_0=0\)`,
2. a strictly negative curvature envelope is known on `\([0,T]\)`,
3. the upper turning bracket satisfies `\(\tau_{\rm hi}^{\rm(tp)}\le T\)`.

Every other local ray is uncertified at this stage and should be deprioritized in the scalarized search.

---

## 9. Best current reading after Stage 189

Stage 187 turned the reduced home stretch into an explicit one-parameter free-quintuple log-ray search.
Stage 188 turned that search into an exact second-order predictor problem.
Stage 189 now turns it into a **certified search sieve**.

The final reduced search is no longer organized by vague statements like
“follow a promising ray.”
It is now organized by exact local data on each candidate ray:

1. the oriented defect `\(H_0\)`,
2. the oriented initial slope `\(K_0\)`,
3. the local curvature envelope `\([\underline K_1,\overline K_1]\)`.

Those three ingredients determine:

- whether the ray is locally admissible,
- an exact local root bracket,
- the bracket width,
- and, when brackets separate, an exact pairwise ray ordering.

So after Stage 189 the next honest continuation is no longer to invent another predictor.
It is to evaluate the actual directional Hessian envelope of `\(\ln\widehat\chi_Q\)` on candidate free-quintuple rays and let the certified bracket/ranking theorem decide which rays are genuinely worth following in the final scalar branch search.
