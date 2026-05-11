# Moving-Throat PDE — Stage 209: Exact Pairwise Ratio Optimizer, Finite Candidate Set, and the Mixed-Ray Winner Theorem

## Status

**Exact within the carried Stage 242 mixed-ray cone / curvature-envelope framework, once a compact pairwise ratio window and the corresponding entrywise Hessian envelopes are supplied on that window.**

This stage does **not** introduce a new constitutive law.
It closes the Stage 242 deferral step: the full one-parameter ratio search on each monotone pairwise cone is reduced to a finite algebraic candidate set.

---

## Purpose

Stage 242 proved that every monotone primitive pair `\((i,j)\)` generates a genuine one-parameter mixed-ray cone
\[
\widehat{\mathbf s}_{ij}(r)=\frac{\widehat{\mathbf e}_i+r\widehat{\mathbf e}_j}{\sqrt{1+r^2}},
\qquad r\ge 0,
\]
and that two canonical rays already isolate the two independent sources of mixed-ray advantage:

- the **gradient-optimal** ray `\(r_{ij}^{\rm grad}=k_j/k_i\)`,
- the **equal-mix synergy** ray `\(r_{ij}^{\rm eq}=1\)`.

That was enough for the first exact mixed-ray screen, but it also made the next theorem gate completely sharp:

> how does one optimize the full certified pairwise bracket over the ratio `\(r\)`, and how can one compare different mixed pairs without running a free-form continuum search on every cone?

This stage answers that.

The main outputs are:

1. the exact **algebraic form** of the pairwise certified upper/lower objective functions,
2. the exact **stationary numerator theorem** for the one-parameter ratio search,
3. the exact **quartic elimination theorem** for interior optimizer candidates,
4. the exact **finite candidate-set theorem** on a compact pairwise ratio window,
5. the exact **optimized pairwise bracket** for the best mixed ray inside one pair,
6. two exact special reductions:
   - the **diagonal-neutral curvature reduction** back to the gradient-optimal ray,
   - the **pair-symmetry reduction** forcing the equal-mix ray to be critical,
7. and the exact **mixed-ray winner theorem** that promotes one pair above all primitive and competing mixed rows.

So Stage 243 is the first place where the pairwise mixed-ray search becomes a finite certified optimization problem rather than a one-parameter continuum deferment.

---

## 1. Carry-forward pairwise cone and the compact ratio window

Fix one monotone primitive pair `\((i,j)\)` with positive primitive slope magnitudes
\[
k_i>0,
\qquad
k_j>0,
\]
and keep the Stage 242 pairwise cone
\[
\widehat{\mathbf s}_{ij}(r)=\frac{\widehat{\mathbf e}_i+r\widehat{\mathbf e}_j}{\sqrt{1+r^2}},
\qquad r\ge0.
\]
The exact positive initial slope magnitude on that cone is
\[
\boxed{
k_{ij}(r)=\frac{k_i+r k_j}{\sqrt{1+r^2}}.
}
\]

Stage 242 allowed arbitrary `\(r\ge0\)`. For the certified optimizer, the natural local object is a **compact pairwise ratio window**
\[
\boxed{
\mathcal R_{ij}:=[0,R_{ij}],
\qquad 0<R_{ij}<\infty,
}
\]
on which the pairwise curvature-envelope data and the local validity radius
\(T_{ij}(r)\) are available.

This does not weaken the search theorem. It is the natural local replacement of the infinite cone for the Stage 240/242 certified program, and it has one major technical payoff:

> every optimized pairwise certified bracket will now come from a **finite** candidate set.

---

## 2. Exact algebraic form of the upper/lower certified objectives

For either envelope label
\[
\star\in\{{\rm lo},{\rm hi}\},
\]
write the three oriented pairwise Hessian envelope entries as
\[
(u_\star,v_\star,w_\star)
=
\begin{cases}
(\underline h_{ii},\underline h_{ij},\underline h_{jj}), & \star={\rm lo},\\[4pt]
(\overline h_{ii},\overline h_{ij},\overline h_{jj}), & \star={\rm hi}.
\end{cases}
\]
Define the exact quadratic numerator coefficients
\[
\boxed{
A_\star:=k_i^2-2H_0 u_\star,
\qquad
B_\star:=2k_i k_j-4H_0 v_\star,
\qquad
C_\star:=k_j^2-2H_0 w_\star.
}
\]
Then the Stage 242 discriminant numerator becomes
\[
\boxed{
\Delta_{ij,\star}^{\sharp}(r)
:=
A_\star+B_\star r+C_\star r^2.
}
\]
Because
\[
k_{ij}(r)^2-2H_0\kappa_{ij,\star}(r)
=
\frac{\Delta_{ij,\star}^{\sharp}(r)}{1+r^2},
\]
the exact certified comparison roots take the explicit algebraic form
\[
\boxed{
\tau_{ij,\star}(r)
=
\frac{2H_0\sqrt{1+r^2}}
{k_i+r k_j+\sqrt{A_\star+B_\star r+C_\star r^2}}.
}
\]
So the full one-parameter pairwise certified search is already reduced to a single square-root rational function of `\(r\)`.

---

## 3. Exact admissible ratio set on a compact cone window

For each envelope label `\(\star\)`, define the admissible ratio set
\[
\boxed{
\mathcal A_{ij,\star}
:=
\Bigl\{
 r\in[0,R_{ij}]:
 \Delta_{ij,\star}^{\sharp}(r)\ge0,
 \tau_{ij,\star}(r)\le T_{ij}(r)
\Bigr\}.
}
\]
So the admissible set keeps exactly the two Stage 242 conditions:

1. the quadratic comparison discriminant is real,
2. the certified upper/lower root stays inside the available local validity interval.

For the rest of this stage, optimization always means optimization **on** `\(\mathcal A_{ij,\star}\)`.

---

## 4. Exact stationary numerator theorem

The exact pairwise certified search is easier to analyze through the denominator functional
\[
\boxed{
\Phi_{ij,\star}(r)
:=
\frac{k_i+r k_j+\sqrt{A_\star+B_\star r+C_\star r^2}}{\sqrt{1+r^2}},
}
\]
so that
\[
\boxed{
\tau_{ij,\star}(r)=\frac{2H_0}{\Phi_{ij,\star}(r)}.
}
\]
Hence minimizing `\(\tau_{ij,\star}\)` is equivalent to maximizing `\(\Phi_{ij,\star}\)`.

Differentiating gives the exact derivative law
\[
\boxed{
\frac{d\Phi_{ij,\star}}{dr}
=
\frac{\mathcal N_{ij,\star}(r)}
{2(1+r^2)^{3/2}\sqrt{A_\star+B_\star r+C_\star r^2}},
}
\]
where the exact stationary numerator is
\[
\boxed{
\mathcal N_{ij,\star}(r)
:=
2(k_j-k_i r)\sqrt{A_\star+B_\star r+C_\star r^2}
+
B_\star+2(C_\star-A_\star)r-B_\star r^2.
}
\]
Therefore every interior optimizer satisfies
\[
\boxed{
\mathcal N_{ij,\star}(r)=0.
}
\]

This is the first exact optimizer theorem of Stage 243:

### Exact stationary numerator theorem

For either the upper or lower certified pairwise objective, every interior ratio optimizer on `\(\mathcal A_{ij,\star}\)` is a zero of the explicit scalar numerator `\(\mathcal N_{ij,\star}\)`.

---

## 5. Exact quartic elimination theorem

The stationary numerator still contains the square root
\(
\sqrt{A_\star+B_\star r+C_\star r^2}
\).
Eliminating it gives an exact quartic candidate polynomial.

Square the stationary equation carefully and define
\[
\boxed{
\mathcal Q_{ij,\star}(r)
:=
\Bigl[B_\star+2(C_\star-A_\star)r-B_\star r^2\Bigr]^2
-
4(k_j-k_i r)^2\bigl(A_\star+B_\star r+C_\star r^2\bigr).
}
\]
Then every interior optimizer satisfies
\[
\boxed{
\mathcal Q_{ij,\star}(r)=0.
}
\]

So the full ratio search has become algebraic in the strongest possible local sense:

### Exact quartic elimination theorem

For each envelope label `\(\star\in\{{\rm lo},{\rm hi}\}\)`, every interior optimizer of the pairwise certified search on the compact ratio window belongs to the positive real root set of the quartic `\(\mathcal Q_{ij,\star}\)`.

The quartic equation may contain extraneous roots from squaring, so the genuine interior optimizer candidates are exactly the roots that also satisfy the unsquared stationary equation `\(\mathcal N_{ij,\star}=0\)`.

---

## 6. Exact finite candidate-set theorem

Because `\(\mathcal Q_{ij,\star}\)`` is quartic, it has at most four real roots.
Therefore each envelope search reduces to a finite set.

For `\(\star\in\{{\rm lo},{\rm hi}\}\)`, define the candidate set
\[
\boxed{
\mathcal C_{ij,\star}
:=
\partial\mathcal A_{ij,\star}
\cup
\Bigl\{r\in \operatorname{int}(\mathcal A_{ij,\star}):
\mathcal N_{ij,\star}(r)=0\Bigr\}.
}
\]
Here `\(\partial\mathcal A_{ij,\star}\)` means the admissible endpoints inherited from the compact ratio window and any internal admissibility boundaries.

Then:

### Exact finite candidate-set theorem

On a compact pairwise ratio window, every optimizer of `\(\tau_{ij,\star}(r)\)` occurs in the finite set `\(\mathcal C_{ij,\star}\)`.

In particular, if the admissible set is a single closed interval inside `\([0,R_{ij}]\)`, then each envelope optimizer requires at most

- the two interval endpoints, and
- up to four admissible quartic roots,

so each of the lower/upper searches needs at most **six** evaluations.

That means one mixed pair requires at most **twelve** exact candidate evaluations to obtain its optimized certified bracket.

---

## 7. Exact optimized pairwise bracket

Define the optimized lower and upper certified pairwise values by
\[
\boxed{
\tau_{ij,\min}^{\rm lo}:=
\min_{r\in\mathcal C_{ij,{\rm lo}}}
\tau_{ij,{\rm lo}}(r),
\qquad
\tau_{ij,\min}^{\rm hi}:=
\min_{r\in\mathcal C_{ij,{\rm hi}}}
\tau_{ij,{\rm hi}}(r).
}
\]
Now let
\[
\tau_{ij,*}^{\rm best}
:=
\min_{r\in\mathcal A_{ij}} \tau_{ij,*}(r)
\]
be the unknown exact best closure time on the pairwise cone, where `\(\tau_{ij,*}(r)\)` is the true mixed-ray closure time.
Because
\[
\tau_{ij,{\rm lo}}(r)
\le
\tau_{ij,*}(r)
\le
\tau_{ij,{\rm hi}}(r)
\qquad (r\in\mathcal A_{ij}),
\]
one gets the exact optimized bracket
\[
\boxed{
\tau_{ij,\min}^{\rm lo}
\le
\tau_{ij,*}^{\rm best}
\le
\tau_{ij,\min}^{\rm hi}.
}
\]
So Stage 243 turns each mixed pair into one finite certified interval for the **best** ray in that pairwise cone.

---

## 8. Two exact special reductions

The full quartic optimizer is the generic case.
But two important exact reductions show how Stage 243 collapses back to the canonical Stage 242 rays when symmetry removes the true one-parameter competition.

### 8.1 Diagonal-neutral curvature reduction

Suppose the chosen envelope is diagonal-neutral and isotropic in the pair:
\[
\boxed{
u_\star=w_\star=\kappa_\star,
\qquad
v_\star=0.}
\]
Then
\[
\kappa_{ij,\star}(r)=\kappa_\star
\]
is constant in `\(r\)`, so minimizing `\(\tau_{ij,\star}(r)\)` is equivalent to maximizing `\(k_{ij}(r)\)`.
By Stage 242, that gives the exact optimizer
\[
\boxed{
r_{ij,\star}^{\rm opt}=r_{ij}^{\rm grad}=\frac{k_j}{k_i}.
}
\]
So when the upper/lower curvature data contain no pairwise directional preference, the full Stage 243 optimizer collapses back to the gradient-optimal ray.

### 8.2 Pair-symmetry reduction

Suppose the pair is symmetric in the sense that
\[
\boxed{
k_i=k_j,
\qquad
u_\star=w_\star.
}
\]
Then the certified objective is invariant under
\[
 r\longmapsto 1/r,
\]
namely
\[
\boxed{
\tau_{ij,\star}(r)=\tau_{ij,\star}(1/r).
}
\]
Therefore
\[
\boxed{
\frac{d\tau_{ij,\star}}{dr}\Big|_{r=1}=0.
}
\]
So in a pair-symmetric mixed cone the equal-mix ray `\(r=1\)` is an exact critical ray for both the lower and upper certified searches.
If the admissible set is logarithmically symmetric and the optimizer is unique, then the Stage 243 optimizer collapses to the Stage 242 equal-mix screen.

These two reductions show that the Stage 242 canonical rays were not arbitrary screens. They are the exact optimizer whenever the pairwise data lose the asymmetry needed to create a true interior mixed-ratio competition.

---

## 9. Exact pairwise promotion and mixed-ray winner theorems

Let `\(\mathfrak P_{\rm mix}\)` be the set of monotone primitive pairs carried into Stage 243, and let the primitive sieve from Stage 241 still carry its certified rows
\[
[\tau_{m,\rm lo}^{\rm prim},\tau_{m,\rm hi}^{\rm prim}].
\]
For each mixed pair, Stage 243 now provides the optimized bracket
\[
[\tau_{ij,\min}^{\rm lo},\tau_{ij,\min}^{\rm hi}].
\]

### Exact pairwise promotion theorem

If, for one pair `\((i,j)\)`,
\[
\boxed{
\tau_{ij,\min}^{\rm hi}
<
\min_m \tau_{m,\rm lo}^{\rm prim},
}
\]
then that pair is already certified to beat the entire primitive sieve.

### Exact mixed-pair winner theorem

If one pair `\((i,j)\)` satisfies
\[
\boxed{
\tau_{ij,\min}^{\rm hi}
<
\min\Bigl(
\min_m \tau_{m,\rm lo}^{\rm prim},
\min_{(p,q)\in\mathfrak P_{\rm mix}\setminus\{(i,j)\}}\tau_{pq,\min}^{\rm lo}
\Bigr),
}
\]
then `\((i,j)\)` is the unique certified mixed-ray winner at the Stage 243 level.

So the full mixed-ray screen is now finite and ordered:

1. primitive certified rows from Stage 241,
2. canonical mixed-ray screens from Stage 242,
3. full pairwise optimized brackets from Stage 243.

No continuum search survives beyond the exact quartic candidate reduction.

---

## 10. Minimal packet for the next stage

After Stage 243, each monotone pair needs only the exact finite packet
\[
\boxed{
\mathcal P_{ij}^{\rm opt}
:=
\bigl(k_i,k_j,H_0,
(u_{\rm lo},v_{\rm lo},w_{\rm lo}),
(u_{\rm hi},v_{\rm hi},w_{\rm hi}),
R_{ij},T_{ij}(r)\bigr),
}
\]
from which the entire optimized bracket is downstream algebra.

So the next honest continuation is no longer to optimize over `\(r\)` again.
The natural continuation is to go **beyond pairwise cones** and ask whether any genuine three-coordinate mixed simplex can beat the Stage 243 pairwise winner.

---

## 11. Best current reading after Stage 243

Stage 242 reduced the first honest mixed-ray search to a one-parameter family plus two exact screen rays.
Stage 243 now finishes the pairwise part of that search:

1. every certified pairwise objective is an explicit square-root rational function of `\(r\)`,
2. every interior optimizer satisfies one exact stationary numerator equation,
3. that stationary equation reduces to one exact quartic,
4. each envelope search therefore collapses to a finite candidate set,
5. every monotone pair now carries one optimized certified bracket,
6. and the mixed-ray winner can be promoted exactly by comparing optimized upper and lower brackets.

So the next theorem gate is no longer “what is the best pairwise mixed ray?”
That question is now finite and exact.
The next real question is whether any **three-coordinate** mixed branch can still beat the Stage 243 pairwise winner.
