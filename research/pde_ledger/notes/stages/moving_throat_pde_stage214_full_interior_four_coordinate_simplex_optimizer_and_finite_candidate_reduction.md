# Moving-Throat PDE — Stage 214: Full Interior Four-Coordinate Simplex Optimizer, Exact Lifted Stationary System, and the Finite Algebraic Candidate Set

## Status

**Exact within the carried Stage 213 four-coordinate mixed-simplex framework, once a compact interior ratio window and the corresponding validity map are supplied, and provided the optimizer does not sit on an artificial outer boundary of that chosen window.**

The whole codimension-one boundary of the four-coordinate simplex is already closed by the Stage 212 / Stage 213 ledger, so this stage addresses only the genuinely new content: the **three-parameter interior** optimizer on a primitive four-coordinate simplex.

---

## Purpose

Stage 213 did three things at once.

1. It proved that every codimension-one face of a primitive four-coordinate simplex is already one of the finished Stage 212 primitive triples.
2. It showed that the first honest support-cardinality-`4` audit is the four imported triple faces plus two canonical interior screens.
3. It isolated the exact next theorem gate:

> how do we optimize the **full** certified interior objective on the three-parameter simplex interior, and can that search be reduced to a finite algebraic candidate set rather than a free-form continuum scan?

This stage answers that.

The main outputs are:

1. the exact three-component stationary numerator system for the interior ratio patch,
2. the exact **lifted polynomial stationary system** with one auxiliary square-root variable,
3. the exact degree ledger `(3,3,3,2)` and the corresponding Bézout candidate bound `54`,
4. the direct square-root-free elimination into three quintic cross-consistency polynomials and three sextic square conditions,
5. the exact finite candidate-set theorem on a compact interior window,
6. two exact special reductions explaining why the Stage 213 canonical screens were the right first interior screens,
7. and the exact interior winner / non-improvement theorems against the already-closed Stage 213 four-face boundary ledger.

So Stage 214 is the four-coordinate analogue of Stage 211: the full interior optimizer is no longer a deferred continuum search. It is a finite algebraic candidate problem.

---

## 1. Carry-forward interior ratio patch and exact certified objective

Fix a monotone primitive quadruple
\[
Q=(i,j,k,l)\in\mathfrak Q_4,
\qquad
k_i,k_j,k_k,k_l>0.
\]
On the interior chart `a_i>0`, use the positive ratio coordinates
\[
r:=\frac{a_j}{a_i}>0,
\qquad
s:=\frac{a_k}{a_i}>0,
\qquad
t:=\frac{a_l}{a_i}>0.
\]
Then the interior simplex point is
\[
\boxed{
\mathbf a(r,s,t)
=
\frac{(1,r,s,t)}{\sqrt{1+r^2+s^2+t^2}}.
}
\]
The oriented initial slope magnitude is
\[
\boxed{
k_{ijkl}(r,s,t)
=
\frac{k_i+r k_j+s k_k+t k_l}{\sqrt{1+r^2+s^2+t^2}}.
}
\]

For either envelope label
\[
\star\in\{{\rm lo},{\rm hi}\},
\]
write the oriented `4 \times 4` Hessian-envelope block entries as
\[
\nu_{ii,\star},\nu_{ij,\star},\nu_{ik,\star},\nu_{il,\star},
\nu_{jj,\star},\nu_{jk,\star},\nu_{jl,\star},
\nu_{kk,\star},\nu_{kl,\star},\nu_{ll,\star}.
\]
Define the exact discriminant coefficients
\[
\boxed{
A_\star:=k_i^2-2H_0\nu_{ii,\star},
\qquad
B_\star:=2k_i k_j-4H_0\nu_{ij,\star},
\qquad
C_\star:=2k_i k_k-4H_0\nu_{ik,\star},
\qquad
D_\star:=2k_i k_l-4H_0\nu_{il,\star},
}
\]
\[
\boxed{
E_\star:=k_j^2-2H_0\nu_{jj,\star},
\qquad
F_\star:=2k_j k_k-4H_0\nu_{jk,\star},
\qquad
G_\star:=2k_j k_l-4H_0\nu_{jl,\star},
}
\]
\[
\boxed{
H_\star:=k_k^2-2H_0\nu_{kk,\star},
\qquad
I_\star:=2k_k k_l-4H_0\nu_{kl,\star},
\qquad
J_\star:=k_l^2-2H_0\nu_{ll,\star}.
}
\]
Then the exact interior discriminant numerator is
\[
\boxed{
\Delta^{\sharp}_{Q,\star}(r,s,t)
:=
A_\star+B_\star r+C_\star s+D_\star t
+E_\star r^2+F_\star rs+G_\star rt+H_\star s^2+I_\star st+J_\star t^2.
}
\]
The exact certified interior root function is
\[
\boxed{
\tau_{Q,\star}(r,s,t)
=
\frac{2H_0\sqrt{1+r^2+s^2+t^2}}
{k_i+r k_j+s k_k+t k_l+\sqrt{\Delta^{\sharp}_{Q,\star}(r,s,t)}}.
}
\]
Equivalently define the denominator functional
\[
\boxed{
\Phi_{Q,\star}(r,s,t)
:=
\frac{k_i+r k_j+s k_k+t k_l+\sqrt{\Delta^{\sharp}_{Q,\star}(r,s,t)}}{\sqrt{1+r^2+s^2+t^2}},
\qquad
\tau_{Q,\star}(r,s,t)=\frac{2H_0}{\Phi_{Q,\star}(r,s,t)}.
}
\]
So interior optimization again means maximizing `\(\Phi\)` or equivalently minimizing `\(\tau\)`.

Let the compact interior ratio window be
\[
\boxed{
\mathcal W_Q:=[0,R_Q]\times[0,S_Q]\times[0,T_Q],
\qquad
0<R_Q,S_Q,T_Q<\infty,
}
\]
and let the local validity map be denoted by
\[
\mathcal T_Q(r,s,t).
\]
The exact admissible interior set is
\[
\boxed{
\mathcal A^{\rm int}_{Q,\star}
:=
\Bigl\{(r,s,t)\in(0,\infty)^3\cap\mathcal W_Q:
\Delta^{\sharp}_{Q,\star}(r,s,t)\ge0,
\ \tau_{Q,\star}(r,s,t)\le \mathcal T_Q(r,s,t)
\Bigr\}.
}
\]

---

## 2. Exact three-component stationary numerator theorem

Differentiate `\(\Phi_{Q,\star}\)` with respect to `\(r\)`, `\(s\)`, and `\(t\)`.
Introduce the exact slope numerators
\[
\boxed{
M_r(r,s,t)
:=(1+r^2+s^2+t^2)k_j-r(k_i+r k_j+s k_k+t k_l)
=
k_j(1+s^2+t^2)-r(k_i+s k_k+t k_l),
}
\]
\[
\boxed{
M_s(r,s,t)
:=(1+r^2+s^2+t^2)k_k-s(k_i+r k_j+s k_k+t k_l)
=
k_k(1+r^2+t^2)-s(k_i+r k_j+t k_l),
}
\]
\[
\boxed{
M_t(r,s,t)
:=(1+r^2+s^2+t^2)k_l-t(k_i+r k_j+s k_k+t k_l)
=
k_l(1+r^2+s^2)-t(k_i+r k_j+s k_k).
}
\]
Also define the exact discriminant-transport numerators
\[
\boxed{
L_{r,Q,\star}
:=
(1+r^2+s^2+t^2)\partial_r\Delta^{\sharp}_{Q,\star}-2r\Delta^{\sharp}_{Q,\star},
}
\]
\[
\boxed{
L_{s,Q,\star}
:=
(1+r^2+s^2+t^2)\partial_s\Delta^{\sharp}_{Q,\star}-2s\Delta^{\sharp}_{Q,\star},
}
\]
\[
\boxed{
L_{t,Q,\star}
:=
(1+r^2+s^2+t^2)\partial_t\Delta^{\sharp}_{Q,\star}-2t\Delta^{\sharp}_{Q,\star}.
}
\]
Then define the exact stationary numerators
\[
\boxed{
\mathcal N_{r,Q,\star}=2M_r\sqrt{\Delta^{\sharp}_{Q,\star}}+L_{r,Q,\star},
}
\]
\[
\boxed{
\mathcal N_{s,Q,\star}=2M_s\sqrt{\Delta^{\sharp}_{Q,\star}}+L_{s,Q,\star},
}
\]
\[
\boxed{
\mathcal N_{t,Q,\star}=2M_t\sqrt{\Delta^{\sharp}_{Q,\star}}+L_{t,Q,\star}.
}
\]
The exact derivative laws are
\[
\boxed{
\frac{\partial\Phi_{Q,\star}}{\partial r}
=
\frac{\mathcal N_{r,Q,\star}}
{2(1+r^2+s^2+t^2)^{3/2}\sqrt{\Delta^{\sharp}_{Q,\star}}},
}
\]
\[
\boxed{
\frac{\partial\Phi_{Q,\star}}{\partial s}
=
\frac{\mathcal N_{s,Q,\star}}
{2(1+r^2+s^2+t^2)^{3/2}\sqrt{\Delta^{\sharp}_{Q,\star}}},
}
\]
\[
\boxed{
\frac{\partial\Phi_{Q,\star}}{\partial t}
=
\frac{\mathcal N_{t,Q,\star}}
{2(1+r^2+s^2+t^2)^{3/2}\sqrt{\Delta^{\sharp}_{Q,\star}}}.
}
\]

### Exact interior stationary numerator theorem

Every interior stationary point of the certified four-coordinate simplex objective satisfies
\[
\boxed{
\mathcal N_{r,Q,\star}(r,s,t)=0,
\qquad
\mathcal N_{s,Q,\star}(r,s,t)=0,
\qquad
\mathcal N_{t,Q,\star}(r,s,t)=0.
}
\]
This is the exact three-parameter analogue of the Stage 211 two-numerator interior theorem.

---

## 3. Exact lifted algebraic stationary system and the `54`-candidate bound

The stationary equations still contain the square root
\(
\sqrt{\Delta^{\sharp}_{Q,\star}(r,s,t)}
\).
The cleanest exact elimination is to introduce one auxiliary variable
\[
\boxed{y:=\sqrt{\Delta^{\sharp}_{Q,\star}(r,s,t)}\ge 0.}
\]
Then the interior stationary system becomes the **lifted polynomial system**
\[
\boxed{
\mathcal F_{r,Q,\star}(r,s,t,y):=2M_r(r,s,t)\,y+L_{r,Q,\star}(r,s,t)=0,
}
\]
\[
\boxed{
\mathcal F_{s,Q,\star}(r,s,t,y):=2M_s(r,s,t)\,y+L_{s,Q,\star}(r,s,t)=0,
}
\]
\[
\boxed{
\mathcal F_{t,Q,\star}(r,s,t,y):=2M_t(r,s,t)\,y+L_{t,Q,\star}(r,s,t)=0,
}
\]
\[
\boxed{
\mathcal F_{\Delta,Q,\star}(r,s,t,y):=y^2-\Delta^{\sharp}_{Q,\star}(r,s,t)=0.
}
\]
The exact degree ledger is
\[
\deg \mathcal F_r=3,
\qquad
\deg \mathcal F_s=3,
\qquad
\deg \mathcal F_t=3,
\qquad
\deg \mathcal F_\Delta=2.
\]
Therefore, whenever the lifted stationary set is zero-dimensional, Bézout gives the exact finite candidate bound
\[
\boxed{3\cdot 3\cdot 3\cdot 2 = 54.}
\]

### Exact lifted pre-candidate set

Define the lifted algebraic pre-candidate set
\[
\boxed{
\widetilde{\mathcal C}^{\rm int,lift}_{Q,\star}
:=
\Bigl\{(r,s,t,y)\in\mathcal W_Q\times\mathbb R_{\ge0}:
\mathcal F_{r,Q,\star}=0,
\ \mathcal F_{s,Q,\star}=0,
\ \mathcal F_{t,Q,\star}=0,
\ \mathcal F_{\Delta,Q,\star}=0
\Bigr\}.
}
\]
Now filter by the original admissibility conditions. Define the exact admissible lifted candidate set
\[
\boxed{
\mathcal C^{\rm int,lift}_{Q,\star}
:=
\Bigl\{(r,s,t,y)\in\widetilde{\mathcal C}^{\rm int,lift}_{Q,\star}:
\tau_{Q,\star}(r,s,t)\le\mathcal T_Q(r,s,t)
\Bigr\}.
}
\]
Because `\(\mathcal F_\Delta=0\)` already enforces `\(y^2=\Delta^{\sharp}\)` and the lifted set is restricted to `\(y\ge0\)`, the original square-root relation is automatic on this admissible set.

### Exact lifted finite candidate-set theorem

Assume the optimizer does not lie on an artificial outer boundary of the chosen ratio window `\(\mathcal W_Q\)` and that the lifted stationary set is zero-dimensional on that window. Then every interior optimizer of `\(\tau_{Q,\star}\)` belongs to the finite admissible set `\(\mathcal C^{\rm int,lift}_{Q,\star}\)`.

So the full interior optimized bracket is obtained by finite evaluation:
\[
\boxed{
\tau_{Q,\min}^{\star,\rm int}
=
\min_{(r,s,t,y)\in\mathcal C^{\rm int,lift}_{Q,\star}}\tau_{Q,\star}(r,s,t),
\qquad \star\in\{{\rm lo},{\rm hi}\}.
}
\]

This is the preferred candidate compiler for the four-coordinate interior problem.

---

## 4. Direct square-root-free elimination and the projected algebraic pre-candidate set

The lifted system is the cleanest exact compiler, but the direct square-root-free elimination is still useful.

### 4.1 Exact quintic cross-consistency polynomials

From
\[
2M_r y + L_r=0,
\qquad
2M_s y + L_s=0,
\qquad
2M_t y + L_t=0,
\]
eliminate `\(y\)` pairwise and define
\[
\boxed{
\mathcal C_{rs,Q,\star}:=M_s L_{r,Q,\star}-M_r L_{s,Q,\star},
}
\]
\[
\boxed{
\mathcal C_{rt,Q,\star}:=M_t L_{r,Q,\star}-M_r L_{t,Q,\star},
}
\]
\[
\boxed{
\mathcal C_{st,Q,\star}:=M_t L_{s,Q,\star}-M_s L_{t,Q,\star}.
}
\]
Each is **quintic** in `(r,s,t)`.

### 4.2 Exact sextic square conditions

Squaring any stationary numerator gives
\[
\boxed{
\mathcal S_{r,Q,\star}:=L_{r,Q,\star}^2-4M_r^2\Delta^{\sharp}_{Q,\star}=0,
}
\]
\[
\boxed{
\mathcal S_{s,Q,\star}:=L_{s,Q,\star}^2-4M_s^2\Delta^{\sharp}_{Q,\star}=0,
}
\]
\[
\boxed{
\mathcal S_{t,Q,\star}:=L_{t,Q,\star}^2-4M_t^2\Delta^{\sharp}_{Q,\star}=0.
}
\]
Each is **sextic** in `(r,s,t)`.

So every interior stationary point satisfies all six projected elimination equations.

### 4.3 One-chart projected candidate bound

If one works on a nondegenerate chart where the chosen `r`-stationary relation is not sitting on an auxiliary elimination degeneracy, then every interior stationary point lies in the projected algebraic pre-candidate set
\[
\boxed{
\widetilde{\mathcal C}^{\rm int}_{Q,\star;r}
:=
\Bigl\{(r,s,t)\in\mathcal W_Q:
\mathcal C_{rs,Q,\star}=0,
\ \mathcal C_{rt,Q,\star}=0,
\ \mathcal S_{r,Q,\star}=0
\Bigr\}.
}
\]
Because the two cross-consistency polynomials are quintic and the square condition is sextic, Bézout gives the projected one-chart bound
\[
\boxed{5\cdot 5\cdot 6 = 150.}
\]
This is a useful projected algebraic upper bound, but the lifted `54`-point system above is the preferred exact candidate compiler because it avoids the chart-dependent square-root degeneracy.

---

## 5. Two exact special reductions

These two reductions justify the Stage 213 canonical interior screens.

### 5.1 Diagonal-isotropic curvature reduction

Suppose the interior Hessian envelope is diagonal-isotropic:
\[
\nu_{ij,\star}=\nu_{ik,\star}=\nu_{il,\star}=\nu_{jk,\star}=\nu_{jl,\star}=\nu_{kl,\star}=0,
\]
\[
\nu_{ii,\star}=\nu_{jj,\star}=\nu_{kk,\star}=\nu_{ll,\star}=:\nu_\star.
\]
Then
\[
\Delta^{\sharp}_{Q,\star}(r,s,t)
=
(k_i+r k_j+s k_k+t k_l)^2-2H_0\nu_\star(1+r^2+s^2+t^2),
\]
so
\[
\tau_{Q,\star}(r,s,t)
=
\frac{2H_0}
{k_{ijkl}(r,s,t)+\sqrt{k_{ijkl}(r,s,t)^2-2H_0\nu_\star}}.
\]
The right-hand side depends only on the first-order slope magnitude `\(k_{ijkl}(r,s,t)\)` and is strictly decreasing in that quantity on the admissible branch.

Therefore the exact optimizer is the Stage 213 gradient-optimal interior ray
\[
\boxed{
\mathbf a_{Q}^{\rm grad}
=
\frac{(k_i,k_j,k_k,k_l)}{\sqrt{k_i^2+k_j^2+k_k^2+k_l^2}}.
}
\]
So the Stage 213 gradient screen is not merely heuristic: it becomes the exact interior optimizer whenever the curvature envelope is diagonal-isotropic.

### 5.2 Full quadruple-symmetry reduction

Suppose the quadruple is fully symmetric:
\[
 k_i=k_j=k_k=k_l=:k,
\]
and the Hessian envelope is permutation-symmetric:
\[
\nu_{ii,\star}=\nu_{jj,\star}=\nu_{kk,\star}=\nu_{ll,\star}=:\nu_{d,\star},
\]
\[
\nu_{ij,\star}=\nu_{ik,\star}=\nu_{il,\star}=\nu_{jk,\star}=\nu_{jl,\star}=\nu_{kl,\star}=:\nu_{x,\star}.
\]
Then the exact stationary numerators satisfy
\[
\boxed{
\mathcal N_{r,Q,\star}(1,1,1)=0,
\qquad
\mathcal N_{s,Q,\star}(1,1,1)=0,
\qquad
\mathcal N_{t,Q,\star}(1,1,1)=0.
}
\]
So the equal-mix barycenter
\[
\boxed{
\mathbf a_Q^{\rm eq}=\frac{(1,1,1,1)}{2}
}
\]
is an exact interior stationary ray. On a symmetric admissible window it is therefore the exact optimizer.

So the Stage 213 equal-mix screen also has an exact theorem regime behind it.

---

## 6. Exact interior winner theorem against the Stage 213 boundary ledger

Keep the imported Stage 213 four-face boundary bracket
\[
\boxed{
\beta_Q^{\rm lo}
\le
\tau_{Q,*}^{\rm best,\partial}
\le
\beta_Q^{\rm hi},
}
\]
where `\(\tau_{Q,*}^{\rm best,\partial}\)` is the unknown true best root on the full codimension-one boundary of the four-coordinate simplex.

Once the finite interior candidate sets are evaluated, the full four-coordinate comparison is immediate.

### Exact interior winner theorem

If
\[
\boxed{
\tau_{Q,\min}^{\rm hi,int} < \beta_Q^{\rm lo},
}
\]
then there exists a genuine interior four-coordinate mixed ray whose certified root lies strictly below the best full-boundary winner.

### Exact interior non-improvement theorem

If
\[
\boxed{
\tau_{Q,\min}^{\rm lo,int} > \beta_Q^{\rm hi},
}
\]
then no admissible interior stationary candidate beats the best boundary winner on that primitive four-coordinate simplex.

So after Stage 214 the full question on a primitive quadruple is no longer a continuum search. It is a finite comparison between:

- one imported exact boundary interval from Stage 213, and
- one finite admissible interior candidate set from the lifted algebraic system above.

The global support-cardinality-`4` ranking across the five primitive quadruples is deferred to the next stage.

---

## 7. What Stage 214 changes in the theorem problem

Stage 213 reduced the first support-cardinality-`4` audit to four imported triple faces plus two canonical interior screens.

Stage 214 now proves that the **full** interior search is still finite and algebraic.

The problem is therefore no longer

> search a three-parameter simplex interior by free-form scanning.

It is now

1. enumerate the finite lifted stationary candidate set,
2. filter it by admissibility,
3. evaluate its lower and upper certified interior minima,
4. compare those against the already-closed Stage 213 boundary ledger.

That is a major reduction in complexity.

---

## 8. Best current summary after Stage 214

The support-cardinality-`4` program is now organized in the same way as the support-`<=3` program.

- Stage 213 closed the boundary and canonical-screen audit.
- Stage 214 closes the **full interior optimizer** as a finite algebraic candidate problem.

So the next theorem gate is now completely sharp:

> rank the five primitive four-coordinate simplices by their full certified interior minima and compare the winning support-cardinality-`4` candidate against the already-finished support-`<=3` ledger.

That is the natural continuation point for Stage 215.
