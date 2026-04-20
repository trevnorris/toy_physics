# Moving-Throat PDE — Stage 193: Three-Coordinate Mixed-Simplex, Exact Interior Gradient/Curvature Synergy, and the Canonical Triple-Screen Audit

## Status

**Exact within the carried Stage-189 scalarized search framework and the fully optimized Stage-192 pairwise certified sieve, once a local oriented `3 x 3` Hessian-envelope block and a triple-simplex validity map are supplied on the chosen three-coordinate patch.**

This stage does **not** yet solve the full two-parameter interior optimizer on a three-coordinate branch.
It is the first exact continuation *beyond* the optimized pairwise cones: the first place where a genuinely interior three-coordinate mixed ray can appear and be audited without redoing the pairwise boundary problem.

---

## Purpose

Stage 192 closed the pairwise problem exactly. Every monotone pair now carries a finite optimized certified bracket, and every boundary pairwise winner can be ranked without a free-form continuum search.

That immediately makes the next theorem gate sharp:

> can a **genuine interior three-coordinate** mixed branch beat the fully optimized pairwise edges, and what is the smallest exact audit that can test that before solving the full two-parameter optimizer?

This stage answers that.

The main outputs are:

1. the exact **positive spherical simplex** for a monotone primitive triple,
2. the exact **boundary-reduction theorem** showing that all three boundary edges are already Stage-192 pairwise cones,
3. the exact **three-coordinate gradient-synergy theorem** and the unique interior gradient-optimal ray,
4. the exact **three-coordinate curvature law** and the theorem that the equal-mix barycenter maximizes total off-diagonal leverage,
5. the exact **fixed-simplex certified bracket** for any admissible interior point,
6. the exact **canonical triple-screen audit** consisting of three optimized boundary rows plus two interior canonical rows,
7. and the exact **interior-screen dominance criterion** that certifies a genuine three-coordinate improvement over all pairwise winners.

So Stage 193 is the three-coordinate analogue of Stage 191, but now with the Stage-192 boundary optimizer already finished and imported.

---

## 1. Carry-forward monotone primitive data and the positive spherical simplex

Keep the free log coordinates
\[
\boldsymbol\ell=(\ell_\lambda,\ell_c,\ell_\gamma,\ell_U,\ell_W)
\]
and the oriented scalar closure residual from Stages 189–192,
\[
H_0>0,
\qquad
\Gamma_i=\sigma_0\,\partial_{\ell_i}h(\boldsymbol\ell_\circ),
\qquad
\sigma_0=\operatorname{sgn}(h(\boldsymbol\ell_\circ)).
\]
For a monotone primitive axis,
\[
\widehat{\mathbf e}_i:=\varepsilon_i\mathbf e_i,
\qquad
\varepsilon_i:=-\operatorname{sgn}(\Gamma_i),
\qquad
k_i:=|\Gamma_i|>0.
\]
So the forward oriented primitive slope is
\[
K_i=-k_i<0.
\]

Now fix a monotone primitive triple `\((i,j,k)\)` with
\[
k_i>0,
\qquad
k_j>0,
\qquad
k_k>0.
\]
The first genuinely three-coordinate search set is the **positive spherical simplex**
\[
\boxed{
\Delta_{ijk}^{+}
:=
\Bigl\{
\mathbf a=(a_i,a_j,a_k)\in\mathbb R_{\ge0}^3:
 a_i^2+a_j^2+a_k^2=1
\Bigr\}.
}
\]
Each point generates the oriented mixed ray
\[
\boxed{
\widehat{\mathbf s}_{ijk}(\mathbf a)
:=
 a_i\widehat{\mathbf e}_i+a_j\widehat{\mathbf e}_j+a_k\widehat{\mathbf e}_k.
}
\]
Interior points have all three coordinates positive. Boundary points lie on one of the three pairwise edges.

---

## 2. Exact boundary-reduction theorem

The three simplex edges are
\[
E_{ij}:=\{a_k=0\},
\qquad
E_{ik}:=\{a_j=0\},
\qquad
E_{jk}:=\{a_i=0\}.
\]
Each edge is exactly one of the Stage-192 pairwise cones.

### 2.1 Edge `\(E_{ij}\)`
With the ratio coordinate
\[
r=\frac{a_j}{a_i}\ge0,
\]
we have
\[
\boxed{
\mathbf a_{ij}(r)=\frac{(1,r,0)}{\sqrt{1+r^2}},
\qquad
\widehat{\mathbf s}_{ijk}(\mathbf a_{ij}(r))
=
\frac{\widehat{\mathbf e}_i+r\widehat{\mathbf e}_j}{\sqrt{1+r^2}}.
}
\]
So `\(E_{ij}\)` is exactly the Stage-192 pairwise cone for `\((i,j)\)`.

### 2.2 Edge `\(E_{ik}\)`
With the ratio coordinate
\[
s=\frac{a_k}{a_i}\ge0,
\]
we have
\[
\boxed{
\mathbf a_{ik}(s)=\frac{(1,0,s)}{\sqrt{1+s^2}},
\qquad
\widehat{\mathbf s}_{ijk}(\mathbf a_{ik}(s))
=
\frac{\widehat{\mathbf e}_i+s\widehat{\mathbf e}_k}{\sqrt{1+s^2}}.
}
\]
So `\(E_{ik}\)` is exactly the Stage-192 pairwise cone for `\((i,k)\)`.

### 2.3 Edge `\(E_{jk}\)`
With the ratio coordinate
\[
u=\frac{a_k}{a_j}\ge0,
\]
we have
\[
\boxed{
\mathbf a_{jk}(\nu)=\frac{(0,1,\nu)}{\sqrt{1+\nu^2}},
\qquad
\widehat{\mathbf s}_{ijk}(\mathbf a_{jk}(\nu))
=
\frac{\widehat{\mathbf e}_j+\nu\widehat{\mathbf e}_k}{\sqrt{1+\nu^2}}.
}
\]
So `\(E_{jk}\)` is exactly the Stage-192 pairwise cone for `\((j,k)\)`.

This is the first theorem-level simplification of Stage 193:

> nothing on the boundary of the three-coordinate simplex is new. The only new content is the **interior** of `\(\Delta_{ijk}^{+}\)`.

---

## 3. Exact three-coordinate gradient-synergy theorem

The positive oriented initial slope magnitude on the simplex ray is
\[
\boxed{
k_{ijk}(\mathbf a)
:=
-\sigma_0\,\nabla_\ell h(\boldsymbol\ell_\circ)\cdot \widehat{\mathbf s}_{ijk}(\mathbf a)
=
a_i k_i+a_j k_j+a_k k_k.
}
\]
So the oriented forward slope is
\[
K_{ijk}(\mathbf a)=-k_{ijk}(\mathbf a)<0.
\]

Because `\(\mathbf a\)` is constrained only by the unit condition `\(a_i^2+a_j^2+a_k^2=1\)`, the exact Lagrange-multiplier solution is immediate:
\[
\boxed{
\mathbf a_{ijk}^{\rm grad}
=
\frac{(k_i,k_j,k_k)}{\sqrt{k_i^2+k_j^2+k_k^2}}.
}
\]
The exact maximum slope magnitude is
\[
\boxed{
\max_{\mathbf a\in\Delta_{ijk}^{+}} k_{ijk}(\mathbf a)
=
\sqrt{k_i^2+k_j^2+k_k^2}.
}
\]
In the interior ratio coordinates on `\(a_i>0\)`,
\[
\boxed{
r_{ijk}^{\rm grad}=\frac{k_j}{k_i},
\qquad
s_{ijk}^{\rm grad}=\frac{k_k}{k_i}.
}
\]

Because all three primitive slopes are positive, the interior gradient-optimal value is strictly larger than every pairwise edge maximum:
\[
\sqrt{k_i^2+k_j^2+k_k^2}
>
\max\!\bigl(\sqrt{k_i^2+k_j^2},\sqrt{k_i^2+k_k^2},\sqrt{k_j^2+k_k^2}\bigr).
\]
So Stage 193 yields the exact first-order theorem:

### Exact three-coordinate gradient-synergy theorem

If all three primitive slopes are nonzero and monotone, a genuine interior three-coordinate mixed ray always beats every pairwise edge at the **first-order** descent level.

That does **not** yet prove it wins after certified curvature enters, but it proves that the interior simplex is carrying genuinely new descent information already at linear order.

---

## 4. Exact three-coordinate curvature law and the total cross-leverage theorem

Let the oriented symmetric Hessian block along the triple ray be
\[
H_{ijk}(\mathbf a;\tau)
=
\begin{pmatrix}
 h_{ii} & h_{ij} & h_{ik}\\
 h_{ij} & h_{jj} & h_{jk}\\
 h_{ik} & h_{jk} & h_{kk}
\end{pmatrix}_{(ijk)}
\]
where, exactly as in Stage 191, each entry is evaluated on the oriented ray
\(
\boldsymbol\ell_\circ+\tau\widehat{\mathbf s}_{ijk}(\mathbf a)
\)
and already includes the sign/orientation factors.

Then the exact second directional derivative is
\[
\boxed{
H_{1,ijk}(\mathbf a;\tau)
=
\mathbf a^T H_{ijk}(\mathbf a;\tau)\mathbf a.
}
\]
Expanding,
\[
\boxed{
H_{1,ijk}
=
a_i^2 h_{ii}+a_j^2 h_{jj}+a_k^2 h_{kk}
+2a_i a_j h_{ij}+2a_i a_k h_{ik}+2a_j a_k h_{jk}.
}
\]
So the total off-diagonal leverage is carried by the exact weight
\[
\boxed{
w_{\Sigma}(\mathbf a)
:=2(a_i a_j+a_i a_k+a_j a_k).
}
\]
This weight admits the identity
\[
\boxed{
w_{\Sigma}(\mathbf a)
=
(a_i+a_j+a_k)^2-1,
}
\]
because `\(a_i^2+a_j^2+a_k^2=1\)` on the simplex.

Now use the exact inequality
\[
3(a_i^2+a_j^2+a_k^2)-(a_i+a_j+a_k)^2
=
(a_i-a_j)^2+(a_i-a_k)^2+(a_j-a_k)^2\ge0.
\]
Since `\(a_i^2+a_j^2+a_k^2=1\)`, it follows that
\[
(a_i+a_j+a_k)^2\le 3,
\qquad
w_{\Sigma}(\mathbf a)\le 2.
\]
Equality holds iff
\[
a_i=a_j=a_k=\frac{1}{\sqrt3}.
\]
So the **equal-mix barycenter**
\[
\boxed{
\mathbf a_{ijk}^{\rm eq}=
\frac{(1,1,1)}{\sqrt3}
}
\]
maximizes the total three-way off-diagonal leverage, with
\[
\boxed{
w_{\Sigma}(\mathbf a_{ijk}^{\rm eq})=2.
}
\]
For comparison, the equal-mix point on any pairwise edge has only
\[
w_{\Sigma}=1.
\]
So the interior barycenter doubles the maximum total cross leverage available on a pairwise equal-mix edge.

### Exact diagonal-neutral reduction

If all off-diagonal Hessian entries vanish,
\[
h_{ij}=h_{ik}=h_{jk}=0,
\]
then
\[
\boxed{
H_{1,ijk}(\mathbf a;\tau)=a_i^2 h_{ii}+a_j^2 h_{jj}+a_k^2 h_{kk},
}
\]
so there is no genuine three-way curvature synergy. In that case the interior simplex carries only the first-order gradient gain from Section 3.

---

## 5. Exact fixed-simplex curvature envelopes and certified local brackets

For either envelope label
\[
\star\in\{{\rm lo},{\rm hi}\},
\]
let the oriented symmetric envelope block be
\[
H_{ijk,\star}=
\begin{pmatrix}
 u_{ii,\star} & u_{ij,\star} & u_{ik,\star}\\
 u_{ij,\star} & u_{jj,\star} & u_{jk,\star}\\
 u_{ik,\star} & u_{jk,\star} & u_{kk,\star}
\end{pmatrix}.
\]
For any fixed simplex point `\(\mathbf a\)`, define the exact envelope curvature scalar
\[
\boxed{
\kappa_{ijk,\star}(\mathbf a)
:=
\mathbf a^T H_{ijk,\star}\mathbf a.
}
\]
Then the exact certified local comparison root at that simplex point is
\[
\boxed{
\tau_{ijk,\star}(\mathbf a)
:=
\mathcal T\bigl(H_0,-k_{ijk}(\mathbf a);\kappa_{ijk,\star}(\mathbf a)\bigr)
=
\frac{2H_0}{k_{ijk}(\mathbf a)+\sqrt{k_{ijk}(\mathbf a)^2-2H_0\kappa_{ijk,\star}(\mathbf a)}}.
}
\]
So the three-coordinate simplex introduces **no new root algebra** at a fixed point. The only new difficulty is the two-parameter interior search over `\(\mathbf a\)`.

### 5.1 Interior ratio coordinates

On the interior patch `\(a_i>0\)`, use the exact ratio coordinates
\[
r=\frac{a_j}{a_i},
\qquad
s=\frac{a_k}{a_i},
\qquad
(a_i,a_j,a_k)=\frac{(1,r,s)}{\sqrt{1+r^2+s^2}}.
\]
Define the six exact quadratic numerator coefficients
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
Then the exact discriminant numerator is
\[
\boxed{
\Delta_{ijk,\star}^{\sharp}(r,s)
:=
A_\star+B_\star r+C_\star s+D_\star r^2+E_\star rs+F_\star s^2.
}
\]
The certified bracket becomes
\[
\boxed{
\tau_{ijk,\star}(r,s)
=
\frac{2H_0\sqrt{1+r^2+s^2}}
{k_i+r k_j+s k_k+
\sqrt{A_\star+B_\star r+C_\star s+D_\star r^2+E_\star rs+F_\star s^2}}.
}
\]
Restricting to any edge recovers the exact Stage-192 pairwise formulas.

---

## 6. The canonical triple-screen audit

Stage 193 does **not** yet solve the full interior optimizer.
The smallest exact screen set is instead:

1. the three exact optimized boundary brackets imported from Stage 192,
   \[
   \mathcal B_{ij}^{\rm opt},
   \qquad
   \mathcal B_{ik}^{\rm opt},
   \qquad
   \mathcal B_{jk}^{\rm opt},
   \]
2. the unique interior gradient-optimal screen point
   \[
   \mathbf a_{ijk}^{\rm grad},
   \]
3. the unique interior equal-mix barycentric screen point
   \[
   \mathbf a_{ijk}^{\rm eq}.
   \]

So the exact Stage-193 triple-screen packet is
\[
\boxed{
\mathcal S_{ijk}^{\rm tri}
:=
\Bigl(
\mathcal B_{ij}^{\rm opt},
\mathcal B_{ik}^{\rm opt},
\mathcal B_{jk}^{\rm opt},
\mathbf a_{ijk}^{\rm grad},
\mathbf a_{ijk}^{\rm eq}
\Bigr).
}
\]
This is the exact three-coordinate analogue of the Stage-191 canonical pairwise screen, but with the three pairwise boundaries already fully optimized and imported rather than re-audited.

---

## 7. Exact interior-screen dominance criterion

Let
\[
\tau_{ijk,\rm hi}^{\rm grad}:=\tau_{ijk,\rm hi}(\mathbf a_{ijk}^{\rm grad}),
\qquad
\tau_{ijk,\rm hi}^{\rm eq}:=\tau_{ijk,\rm hi}(\mathbf a_{ijk}^{\rm eq}),
\]
and keep the three Stage-192 optimized boundary lower brackets
\[
\tau_{ij,\min}^{\rm lo},
\qquad
\tau_{ik,\min}^{\rm lo},
\qquad
\tau_{jk,\min}^{\rm lo}.
\]
Then the first exact three-coordinate certification rule is immediate:

### Exact interior-screen dominance theorem

If either canonical interior screen satisfies
\[
\boxed{
\tau_{ijk,\rm hi}^{\rm can}
<
\min\bigl(
\tau_{ij,\min}^{\rm lo},
\tau_{ik,\min}^{\rm lo},
\tau_{jk,\min}^{\rm lo}
\bigr),
\qquad
\text{can}\in\{{\rm grad},{\rm eq}\},
}
\]
then there exists a **genuine interior three-coordinate mixed ray** whose certified root lies strictly below every pairwise boundary winner.

The proof is immediate:

- the actual root at the chosen interior screen point is bounded above by its certified upper bracket,
- the actual pairwise boundary winners are bounded below by their certified lower brackets,
- so the interior screen point already certifies a strict interior improvement over all pairwise edges.

### Exact canonical non-improvement filter

Conversely, if
\[
\boxed{
\min\bigl(
\tau_{ijk,\rm lo}^{\rm grad},
\tau_{ijk,\rm lo}^{\rm eq}
\bigr)
>
\min\bigl(
\tau_{ij,\min}^{\rm hi},
\tau_{ik,\min}^{\rm hi},
\tau_{jk,\min}^{\rm hi}
\bigr),
}
\]
then neither canonical interior screen beats the best pairwise edge winner.

That is **not** a full no-go theorem for the interior simplex, because the full two-parameter optimizer has not yet been solved. But it is the first exact filter that can rule out the two canonical interior mechanisms before the full interior search is attempted.

---

## 8. Minimal packet for the next stage

After Stage 193, the boundary problem is already solved and the interior problem has only one unresolved part: the genuine two-parameter optimizer on the simplex interior.

The smallest exact packet for that next stage is
\[
\boxed{
\mathcal P_{ijk}^{\rm tri}
:=
\Bigl(
H_0,
(k_i,k_j,k_k),
H_{ijk,\rm lo},
H_{ijk,\rm hi},
T_{ijk}(\mathbf a)
\Bigr),
}
\]
where

- `\(H_0\)` is the oriented logarithmic defect,
- `\((k_i,k_j,k_k)\)` are the three primitive slope magnitudes,
- `\(H_{ijk,\rm lo},H_{ijk,\rm hi}\)` are the lower/upper oriented `3 x 3` Hessian-envelope blocks,
- `\(T_{ijk}(\mathbf a)\)` is the local validity-radius map on the simplex patch.

The natural continuation is now completely sharp:

> solve the full interior simplex optimizer by reducing the two-parameter stationary system to a finite algebraic candidate set, just as Stage 192 did for the one-parameter pairwise cones.

That is the natural content of Stage 194.

---

## 9. Best current reading after Stage 193

Stage 192 already proved that the pairwise boundary problem is finite and exact.
Stage 193 now shows what is genuinely new at the first three-coordinate level:

1. every boundary edge of the three-coordinate simplex is already one of the Stage-192 pairwise cones,
2. a genuine interior triple ray always beats the pairwise edges at the first-order gradient level,
3. the equal-mix barycenter uniquely maximizes the total three-way off-diagonal Hessian leverage,
4. every fixed interior simplex point already carries an exact certified local bracket,
5. and the smallest honest interior audit is now a five-row screen:
   - three optimized boundary rows,
   - one interior gradient-optimal row,
   - one interior equal-mix row.

So the next theorem gate is no longer “is there any three-coordinate effect at all?”
There is.
The next real question is whether the full two-parameter interior simplex optimizer can produce a certified winner that beats both canonical interior screens and all optimized pairwise edges.
