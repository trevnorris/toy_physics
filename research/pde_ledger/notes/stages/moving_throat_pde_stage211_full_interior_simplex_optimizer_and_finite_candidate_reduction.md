# Moving-Throat PDE — Stage 211: Full Interior Triple-Simplex Optimizer, Exact Stationary Elimination, and the Finite Algebraic Candidate Set

## Status

**Exact within the carried Stage 244 interior-simplex framework, once a compact interior ratio window and the corresponding validity map are supplied, and provided the optimizer is not sitting on an artificial chart boundary of that window.**

The true simplex boundaries are already closed by Stage 243 and imported into Stage 244, so this stage addresses only the genuinely new content: the **interior** optimizer on a three-coordinate mixed simplex.

---

## Purpose

Stage 244 proved three things at once.

1. Every boundary face of the positive three-coordinate simplex is already one of the Stage 243 pairwise cones.
2. A genuine interior triple ray always beats the pairwise edges at the purely first-order gradient level.
3. The smallest honest interior audit is the five-row screen made from three optimized pairwise boundary rows plus two canonical interior rows.

That makes the next theorem gate completely sharp:

> how do we optimize the full certified interior objective on the two-parameter simplex interior, and can that search be reduced to a finite algebraic candidate set rather than a free-form continuum scan?

This stage answers that.

The main outputs are:

1. the exact two-component stationary system for the interior ratio patch,
2. the exact elimination of the square root into one **quartic** and one **sextic** polynomial condition,
3. the exact finite algebraic pre-candidate set for all interior stationary rays,
4. two exact special reductions explaining why the Stage 244 canonical interior screens were the right first screens,
5. and the exact full interior winner theorem against the already-optimized pairwise boundaries.

So Stage 245 is the three-coordinate analogue of Stage 243: the interior optimizer is no longer a deferred continuum search; it is a finite algebraic candidate problem.

---

## 1. Carry-forward interior ratio patch and exact certified objective

Fix a monotone primitive triple `\((i,j,k)\)` with positive oriented primitive slope magnitudes
\[
 k_i>0,
 \qquad
 k_j>0,
 \qquad
 k_k>0.
\]
On the interior patch `\(a_i>0\)`, use the positive ratio coordinates
\[
 r:=\frac{a_j}{a_i}>0,
 \qquad
 s:=\frac{a_k}{a_i}>0.
\]
Then the positive spherical-simplex point is
\[
\boxed{
\mathbf a(r,s)
=
\frac{(1,r,s)}{\sqrt{1+r^2+s^2}}.
}
\]
The oriented initial slope magnitude is
\[
\boxed{
 k_{ijk}(r,s)
 =
 \frac{k_i+r k_j+s k_k}{\sqrt{1+r^2+s^2}}.
}
\]
For either envelope label
\[
\star\in\{{\rm lo},{\rm hi}\},
\]
write the oriented `\(3\times 3\)` Hessian-envelope block entries as
\[
 u_{ii,\star},\ u_{ij,\star},\ u_{ik,\star},\ u_{jj,\star},\ u_{jk,\star},\ u_{kk,\star}.
\]
Define the exact quadratic discriminant coefficients
\[
\boxed{
A_\star:=k_i^2-2H_0 u_{ii,\star},
\qquad
B_\star:=2k_i k_j-4H_0 u_{ij,\star},
\qquad
C_\star:=2k_i k_k-4H_0 u_{ik,\star},
}
\]
\[
\boxed{
D_\star:=k_j^2-2H_0 u_{jj,\star},
\qquad
E_\star:=2k_j k_k-4H_0 u_{jk,\star},
\qquad
F_\star:=k_k^2-2H_0 u_{kk,\star}.
}
\]
The exact interior discriminant numerator is then
\[
\boxed{
\Delta^{\sharp}_{ijk,\star}(r,s)
:=
A_\star+B_\star r+C_\star s+D_\star r^2+E_\star r s+F_\star s^2.
}
\]
The Stage 244 fixed-point comparison root becomes
\[
\boxed{
\tau_{ijk,\star}(r,s)
=
\frac{2H_0\sqrt{1+r^2+s^2}}
{k_i+r k_j+s k_k+\sqrt{\Delta^{\sharp}_{ijk,\star}(r,s)}}.
}
\]
Equivalently, define the denominator functional
\[
\boxed{
\Phi_{ijk,\star}(r,s)
:=
\frac{k_i+r k_j+s k_k+\sqrt{\Delta^{\sharp}_{ijk,\star}(r,s)}}{\sqrt{1+r^2+s^2}},
\qquad
\tau_{ijk,\star}(r,s)=\frac{2H_0}{\Phi_{ijk,\star}(r,s)}.
}
\]
So interior optimization again means maximizing `\(\Phi\)` or equivalently minimizing `\(\tau\)`.

Let
\[
\boxed{
\mathcal W_{ijk}:=[0,R_{ijk}]\times[0,S_{ijk}],
\qquad
0<R_{ijk},S_{ijk}<\infty,
}
\]
be a compact interior ratio window on which the envelope data and the local validity map `\(T_{ijk}(r,s)\)` are supplied.
The exact admissible interior set is
\[
\boxed{
\mathcal A^{\rm int}_{ijk,\star}
:=
\Bigl\{(r,s)\in(0,\infty)^2\cap\mathcal W_{ijk}:
\Delta^{\sharp}_{ijk,\star}(r,s)\ge0,
\ \tau_{ijk,\star}(r,s)\le T_{ijk}(r,s)
\Bigr\}.
}
\]

---

## 2. Exact two-component stationary numerator theorem

Differentiate `\(\Phi_{ijk,\star}\)` with respect to `\(r\)` and `\(s\)`.
Introduce the exact slope numerators
\[
\boxed{
M_r(r,s):=(1+r^2+s^2)k_j-r(k_i+r k_j+s k_k)=k_j(1+s^2)-r(k_i+s k_k),
}
\]
\[
\boxed{
M_s(r,s):=(1+r^2+s^2)k_k-s(k_i+r k_j+s k_k)=k_k(1+r^2)-s(k_i+r k_j).
}
\]
Also define the exact discriminant-transport numerators
\[
\boxed{
L_{r,\star}(r,s)
:=
(1+r^2+s^2)\,\partial_r\Delta^{\sharp}_{ijk,\star}(r,s)-2r\,\Delta^{\sharp}_{ijk,\star}(r,s),
}
\]
\[
\boxed{
L_{s,\star}(r,s)
:=
(1+r^2+s^2)\,\partial_s\Delta^{\sharp}_{ijk,\star}(r,s)-2s\,\Delta^{\sharp}_{ijk,\star}(r,s).
}
\]
Then the exact derivative laws are
\[
\boxed{
\frac{\partial\Phi_{ijk,\star}}{\partial r}
=
\frac{\mathcal N_{r,\star}(r,s)}
{2(1+r^2+s^2)^{3/2}\sqrt{\Delta^{\sharp}_{ijk,\star}(r,s)}},
}
\]
\[
\boxed{
\frac{\partial\Phi_{ijk,\star}}{\partial s}
=
\frac{\mathcal N_{s,\star}(r,s)}
{2(1+r^2+s^2)^{3/2}\sqrt{\Delta^{\sharp}_{ijk,\star}(r,s)}},
}
\]
with exact stationary numerators
\[
\boxed{
\mathcal N_{r,\star}(r,s)=2M_r(r,s)\sqrt{\Delta^{\sharp}_{ijk,\star}(r,s)}+L_{r,\star}(r,s),
}
\]
\[
\boxed{
\mathcal N_{s,\star}(r,s)=2M_s(r,s)\sqrt{\Delta^{\sharp}_{ijk,\star}(r,s)}+L_{s,\star}(r,s).
}
\]
So every interior stationary point satisfies
\[
\boxed{
\mathcal N_{r,\star}(r,s)=0,
\qquad
\mathcal N_{s,\star}(r,s)=0.
}
\]

### Exact interior stationary numerator theorem

For either the upper or lower certified interior objective, every interior stationary mixed ray is a common zero of the two explicit square-root numerators `\(\mathcal N_{r,\star}\)` and `\(\mathcal N_{s,\star}\)`.

This is the exact two-parameter analogue of the Stage 243 pairwise stationary-numerator theorem.

---

## 3. Exact quartic-sextic elimination theorem

The stationary equations still contain the square root
\(
\sqrt{\Delta^{\sharp}_{ijk,\star}(r,s)}
\).
Eliminate it in two steps.

### 3.1 Exact quartic cross-consistency polynomial

Multiply the stationary numerators crosswise to remove the square root.
Define
\[
\boxed{
\mathcal C_{ijk,\star}(r,s)
:=
M_s(r,s)L_{r,\star}(r,s)-M_r(r,s)L_{s,\star}(r,s).
}
\]
Then every interior stationary point satisfies
\[
\boxed{
\mathcal C_{ijk,\star}(r,s)=0.
}
\]
This polynomial is **quartic** in `\((r,s)\)`.

### 3.2 Exact sextic square conditions

From
\[
2M_r\sqrt{\Delta^{\sharp}}=-L_{r,\star},
\qquad
2M_s\sqrt{\Delta^{\sharp}}=-L_{s,\star},
\]
squaring gives the exact sextic eliminants
\[
\boxed{
\mathcal S_{r,ijk,\star}(r,s)
:=
L_{r,\star}(r,s)^2-4M_r(r,s)^2\,\Delta^{\sharp}_{ijk,\star}(r,s),
}
\]
\[
\boxed{
\mathcal S_{s,ijk,\star}(r,s)
:=
L_{s,\star}(r,s)^2-4M_s(r,s)^2\,\Delta^{\sharp}_{ijk,\star}(r,s).
}
\]
Every interior stationary point satisfies
\[
\boxed{
\mathcal S_{r,ijk,\star}(r,s)=0,
\qquad
\mathcal S_{s,ijk,\star}(r,s)=0.
}
\]
Both are **sextic** in `\((r,s)\)`.

### Exact elimination theorem

Every interior stationary point of the certified three-coordinate simplex objective is a common zero of
\[
\boxed{
\mathcal C_{ijk,\star}(r,s)=0,
\qquad
\mathcal S_{r,ijk,\star}(r,s)=0,
}
\]
and also of the equivalent pair
\[
\boxed{
\mathcal C_{ijk,\star}(r,s)=0,
\qquad
\mathcal S_{s,ijk,\star}(r,s)=0.
}
\]

Because the first polynomial is quartic and the second is sextic, the algebraic pre-candidate set is finite whenever the intersection is zero-dimensional. By Bézout, it has at most
\[
\boxed{4\cdot 6 = 24}
\]
complex roots counting multiplicity.

So the full interior stationary problem has been reduced from a two-parameter continuum search to a **finite algebraic pre-candidate set**.

---

## 4. Exact finite candidate-set theorem on a compact interior window

Define the algebraic pre-candidate set
\[
\boxed{
\widetilde{\mathcal C}^{\rm int}_{ijk,\star}
:=
\Bigl\{(r,s)\in\mathcal W_{ijk}:
\mathcal C_{ijk,\star}(r,s)=0,
\ \mathcal S_{r,ijk,\star}(r,s)=0
\Bigr\}.
}
\]
Because `\(\mathcal C\)` is quartic and `\(\mathcal S_r\)` is sextic, this set is finite whenever the common zeros are isolated.

Now filter by the original square-root and validity conditions. Define the exact admissible stationary candidate set
\[
\boxed{
\mathcal C^{\rm int}_{ijk,\star}
:=
\Bigl\{(r,s)\in\widetilde{\mathcal C}^{\rm int}_{ijk,\star}:
\Delta^{\sharp}_{ijk,\star}(r,s)\ge0,
\ \tau_{ijk,\star}(r,s)\le T_{ijk}(r,s),
\ \mathcal N_{r,\star}(r,s)=0,
\ \mathcal N_{s,\star}(r,s)=0
\Bigr\}.
}
\]
This final filtering removes any spurious roots introduced by squaring.

### Exact finite interior candidate-set theorem

Assume the optimizer does not lie on an artificial outer boundary of the chosen ratio window `\(\mathcal W_{ijk}\)`. Then every interior optimizer of `\(\tau_{ijk,\star}\)` belongs to the finite admissible set `\(\mathcal C^{\rm int}_{ijk,\star}\)`.

So the interior optimized bracket is obtained by finite evaluation:
\[
\boxed{
\tau_{ijk,\min}^{\star,\rm int}
=
\min_{(r,s)\in\mathcal C^{\rm int}_{ijk,\star}} \tau_{ijk,\star}(r,s),
\qquad \star\in\{{\rm lo},{\rm hi}\}.
}
\]

The true simplex boundaries are already solved by Stage 243 and imported in Stage 244, so this is the only genuinely new interior algebraic candidate problem.

---

## 5. Two exact special reductions

These two reductions justify the Stage 244 canonical interior screens.

### 5.1 Diagonal-isotropic curvature reduction

Suppose the interior Hessian envelope is diagonal-isotropic:
\[
 u_{ij,\star}=u_{ik,\star}=u_{jk,\star}=0,
 \qquad
 u_{ii,\star}=u_{jj,\star}=u_{kk,\star}=u_\star.
\]
Then
\[
\Delta^{\sharp}_{ijk,\star}(r,s)
=
(k_i+r k_j+s k_k)^2-2H_0 u_\star(1+r^2+s^2),
\]
so
\[
\tau_{ijk,\star}(r,s)
=
\frac{2H_0}
{k_{ijk}(r,s)+\sqrt{k_{ijk}(r,s)^2-2H_0 u_\star}}.
\]
The right-hand side depends only on the first-order slope magnitude `\(k_{ijk}(r,s)\)`, and is strictly decreasing in that quantity on the admissible branch.

Therefore the exact optimizer is the Stage 244 gradient-optimal interior ray
\[
\boxed{
\mathbf a_{ijk}^{\rm grad}
=
\frac{(k_i,k_j,k_k)}{\sqrt{k_i^2+k_j^2+k_k^2}}.
}
\]
So the Stage 244 gradient screen is not merely heuristic: it becomes the exact interior optimizer whenever the curvature envelope is diagonal-isotropic.

### 5.2 Full triple-symmetry reduction

Suppose the triple is fully symmetric:
\[
 k_i=k_j=k_k=:k,
\]
and the Hessian envelope is permutation-symmetric:
\[
 u_{ii,\star}=u_{jj,\star}=u_{kk,\star}=:u_{d,\star},
 \qquad
 u_{ij,\star}=u_{ik,\star}=u_{jk,\star}=:u_{x,\star}.
\]
Then the exact stationary numerators satisfy
\[
\boxed{
\mathcal N_{r,\star}(1,1)=0,
\qquad
\mathcal N_{s,\star}(1,1)=0.
}
\]
So the equal-mix barycenter
\[
\boxed{
\mathbf a_{ijk}^{\rm eq}=rac{(1,1,1)}{\sqrt 3}
}
\]
is an exact interior stationary ray. On a symmetric admissible window it is therefore the exact optimizer.

So the Stage 244 equal-mix screen also has an exact theorem regime behind it.

---

## 6. Exact interior winner theorem against the optimized pairwise boundaries

Keep the Stage 243 optimized pairwise boundary brackets
\[
\tau_{ij,\min}^{\rm lo},\qquad \tau_{ik,\min}^{\rm lo},\qquad \tau_{jk,\min}^{\rm lo},
\]
and similarly the corresponding upper certified values
\[
\tau_{ij,\min}^{\rm hi},\qquad \tau_{ik,\min}^{\rm hi},\qquad \tau_{jk,\min}^{\rm hi}.
\]
Once the finite interior candidate sets are evaluated, the full three-coordinate comparison is immediate.

### Exact interior winner theorem

If
\[
\boxed{
\tau_{ijk,\min}^{\rm hi,int}
<
\min\bigl(
\tau_{ij,\min}^{\rm lo},
\tau_{ik,\min}^{\rm lo},
\tau_{jk,\min}^{\rm lo}
\bigr),
}
\]
then there exists a genuine interior three-coordinate mixed ray whose certified root lies strictly below every optimized pairwise boundary winner.

### Exact interior non-improvement theorem

If
\[
\boxed{
\tau_{ijk,\min}^{\rm lo,int}
>
\min\bigl(
\tau_{ij,\min}^{\rm hi},
\tau_{ik,\min}^{\rm hi},
\tau_{jk,\min}^{\rm hi}
\bigr),
}
\]
then no admissible interior stationary candidate beats the best pairwise boundary winner.

So after Stage 245 the three-coordinate simplex no longer needs a free-form interior search to be compared with the already-optimized pairwise cones. The comparison is again finite and certified.

---

## 7. Minimal packet for the next stage

After Stage 245, the genuinely new three-coordinate content is compressed to the finite interior candidate problem plus the already-imported pairwise boundaries.

The smallest exact packet for the next stage is
\[
\boxed{
\mathcal P_{ijk}^{\rm int}
:=
\Bigl(
H_0,
(k_i,k_j,k_k),
H_{ijk,\rm lo},
H_{ijk,\rm hi},
T_{ijk}(r,s),
\mathcal C_{ijk,\star},
\mathcal S_{r,ijk,\star}
\Bigr).
}
\]
The natural continuation is now completely sharp:

> rank the full set of primitive triples by their finite certified interior winners, splice those winners together with the optimized boundary ledger, and determine whether any three-coordinate simplex truly beats all pairwise branches in the carried local search.

That is the natural content of Stage 246.

---

## 8. Best current reading after Stage 245

Stage 243 turned every pairwise cone into a finite certified candidate problem.
Stage 244 proved that the three-coordinate simplex has genuinely new interior content and isolated the first canonical interior screens.
Stage 245 now closes the missing algebraic step:

1. the interior stationary system is explicit,
2. the square root is eliminated exactly into one quartic and one sextic condition,
3. the admissible interior optimizer sits in a finite algebraic candidate set,
4. the gradient-optimal and equal-mix screens are justified by exact special-reduction theorems,
5. and the full interior-versus-boundary comparison is again a finite certified ranking problem.

So the next theorem gate is no longer “how do we search the interior at all?”
That problem is finished.
The next question is which primitive triples, if any, actually produce certified interior winners once the real PDE-derived local packet is supplied.
