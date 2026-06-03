# Moving-Throat PDE — Stage 213: Four-Coordinate Mixed Simplex, Exact Face Reduction, and the Support-Cardinality-4 Theorem Gate

## Status

**Exact within the carried Stage 212 up-to-three-coordinate certified sieve, once a local oriented `4 x 4` Hessian-envelope block and the corresponding validity map are supplied on the chosen four-coordinate patch.**

This stage does **not** yet solve the full three-parameter interior optimizer on a primitive four-coordinate simplex.
It is the first exact continuation *beyond* the fully ranked support-`<=3` ledger: the first place where a genuinely interior support-cardinality-`4` mixed ray can appear and be audited without reopening the solved pairwise and triple boundary problems.

---

## Purpose

Stage 212 closed the full support-`<=3` local search.
Every primitive pair is already a finite certified candidate problem, every primitive triple carries a closed-simplex certified interval, and the full local sieve through three active primitive coordinates has been reduced to one finite certified ledger.

That makes the next theorem gate completely sharp:

> can a **genuine interior four-coordinate** mixed branch beat the already-finished support-`<=3` ledger, and what is the smallest exact audit that can test that before solving the full three-parameter interior optimizer?

This stage answers that.

The main outputs are:

1. the exact combinatorial ledger for the five primitive four-coordinate simplices,
2. the exact **positive spherical four-simplex** and its full face reduction back to the Stage 212 triple packets,
3. the exact **four-coordinate gradient-synergy theorem** and the unique interior gradient-optimal ray,
4. the exact **four-coordinate curvature law** and the theorem that the equal-mix barycenter maximizes the total six-way off-diagonal leverage,
5. the exact **fixed-simplex certified bracket** for any admissible interior point,
6. the exact **canonical quadruple-screen audit** consisting of four imported triple faces plus two interior canonical rows,
7. and the exact **support-cardinality-4 theorem gate** that certifies a genuine four-coordinate improvement over the entire up-to-three-coordinate ledger whenever one canonical interior screen already wins.

So Stage 213 is the four-coordinate analogue of Stage 210, but now with the whole support-`<=3` boundary ledger already finished and imported from Stage 212.

---

## 1. Carry-forward primitive quadruple ledger and the positive spherical four-simplex

Keep the free log coordinates
\[
\boldsymbol\ell=(\ell_\lambda,\ell_c,\ell_\gamma,\ell_U,\ell_W)
\]
and the oriented primitive slope magnitudes from the carried scalarized closure chain,
\[
H_0>0,
\qquad
\widehat{\mathbf e}_i:=\varepsilon_i\mathbf e_i,
\qquad
\varepsilon_i:=-\operatorname{sgn}(\Gamma_i),
\qquad
k_i:=|\Gamma_i|>0.
\]
So every monotone primitive axis has forward oriented slope
\[
K_i=-k_i<0.
\]

Let the primitive free-axis set be
\[
\boxed{
\mathfrak I_5:=\{\lambda,c,\gamma,U,W\}.
}
\]
The primitive quadruple set is then
\[
\boxed{
\mathfrak Q_4
:=
\bigl\{\{i,j,k,l\}\subset\mathfrak I_5:i<j<k<l\bigr\},
\qquad
\#\mathfrak Q_4=\binom54=5.
}
\]
Every primitive triple belongs to exactly two primitive quadruples, because once three axes are fixed there are exactly two remaining primitive directions that can complete the quadruple.

Now fix a monotone primitive quadruple
\[
Q=(i,j,k,l)\in\mathfrak Q_4,
\qquad
k_i,k_j,k_k,k_l>0.
\]
The first genuinely four-coordinate search set is the **positive spherical four-simplex**
\[
\boxed{
\Delta_{ijkl}^{+}
:=
\Bigl\{
\mathbf a=(a_i,a_j,a_k,a_l)\in\mathbb R_{\ge0}^4:
 a_i^2+a_j^2+a_k^2+a_l^2=1
\Bigr\}.
}
\]
Each point generates the oriented mixed ray
\[
\boxed{
\widehat{\mathbf s}_{ijkl}(\mathbf a)
:=
 a_i\widehat{\mathbf e}_i+a_j\widehat{\mathbf e}_j+a_k\widehat{\mathbf e}_k+a_l\widehat{\mathbf e}_l.
}
\]
Interior points have all four coordinates positive. Boundary points lie on one of the four three-coordinate faces.

---

## 2. Exact face-reduction theorem

The four codimension-one faces are
\[
F_{ijk}:=\{a_l=0\},
\qquad
F_{ijl}:=\{a_k=0\},
\qquad
F_{ikl}:=\{a_j=0\},
\qquad
F_{jkl}:=\{a_i=0\}.
\]
Each face is exactly one of the Stage 212 primitive triple **closed** simplices.

### 2.1 Face `\(F_{ijk}\)`
With positive ratio coordinates
\[
r=\frac{a_j}{a_i}\ge0,
\qquad
s=\frac{a_k}{a_i}\ge0,
\]
we have
\[
\boxed{
\mathbf a_{ijk}(r,s)=\frac{(1,r,s,0)}{\sqrt{1+r^2+s^2}},
\qquad
\widehat{\mathbf s}_{ijkl}(\mathbf a_{ijk}(r,s))
=
\frac{\widehat{\mathbf e}_i+r\widehat{\mathbf e}_j+s\widehat{\mathbf e}_k}{\sqrt{1+r^2+s^2}}.
}
\]
So `\(F_{ijk}\)` is exactly the Stage 212 closed simplex attached to the primitive triple `\((i,j,k)\)`.

### 2.2 Face `\(F_{ijl}\)`
With positive ratio coordinates
\[
r=\frac{a_j}{a_i}\ge0,
\qquad
t=\frac{a_l}{a_i}\ge0,
\]
we have
\[
\boxed{
\mathbf a_{ijl}(r,t)=\frac{(1,r,0,t)}{\sqrt{1+r^2+t^2}},
\qquad
\widehat{\mathbf s}_{ijkl}(\mathbf a_{ijl}(r,t))
=
\frac{\widehat{\mathbf e}_i+r\widehat{\mathbf e}_j+t\widehat{\mathbf e}_l}{\sqrt{1+r^2+t^2}}.
}
\]
So `\(F_{ijl}\)` is exactly the Stage 212 closed simplex for the primitive triple `\((i,j,l)\)`.

### 2.3 Face `\(F_{ikl}\)`
With positive ratio coordinates
\[
s=\frac{a_k}{a_i}\ge0,
\qquad
t=\frac{a_l}{a_i}\ge0,
\]
we have
\[
\boxed{
\mathbf a_{ikl}(s,t)=\frac{(1,0,s,t)}{\sqrt{1+s^2+t^2}},
\qquad
\widehat{\mathbf s}_{ijkl}(\mathbf a_{ikl}(s,t))
=
\frac{\widehat{\mathbf e}_i+s\widehat{\mathbf e}_k+t\widehat{\mathbf e}_l}{\sqrt{1+s^2+t^2}}.
}
\]
So `\(F_{ikl}\)` is exactly the Stage 212 closed simplex for the primitive triple `\((i,k,l)\)`.

### 2.4 Face `\(F_{jkl}\)`
With positive ratio coordinates
\[
u=\frac{a_k}{a_j}\ge0,
\qquad
w=\frac{a_l}{a_j}\ge0,
\]
we have
\[
\boxed{
\mathbf a_{jkl}(\nu,w)=\frac{(0,1,\nu,w)}{\sqrt{1+\nu^2+w^2}},
\qquad
\widehat{\mathbf s}_{ijkl}(\mathbf a_{jkl}(\nu,w))
=
\frac{\widehat{\mathbf e}_j+\nu\widehat{\mathbf e}_k+w\widehat{\mathbf e}_l}{\sqrt{1+\nu^2+w^2}}.
}
\]
So `\(F_{jkl}\)` is exactly the Stage 212 closed simplex for the primitive triple `\((j,k,l)\)`.

This gives the first exact four-coordinate simplification:

> nothing on the boundary of the four-coordinate simplex is new. The entire codimension-one boundary is already closed by the Stage 212 primitive-triple ledger.

Define the imported Stage 212 full-face certified intervals
\[
\mathcal I_{ijk}^{\triangle},
\qquad
\mathcal I_{ijl}^{\triangle},
\qquad
\mathcal I_{ikl}^{\triangle},
\qquad
\mathcal I_{jkl}^{\triangle},
\]
with lower and upper face minima
\[
\boxed{
\beta_Q^{\rm lo}
:=
\min\bigl(
\tau_{ijk,\min}^{\rm lo,\triangle},
\tau_{ijl,\min}^{\rm lo,\triangle},
\tau_{ikl,\min}^{\rm lo,\triangle},
\tau_{jkl,\min}^{\rm lo,\triangle}
\bigr),
}
\]
\[
\boxed{
\beta_Q^{\rm hi}
:=
\min\bigl(
\tau_{ijk,\min}^{\rm hi,\triangle},
\tau_{ijl,\min}^{\rm hi,\triangle},
\tau_{ikl,\min}^{\rm hi,\triangle},
\tau_{jkl,\min}^{\rm hi,\triangle}
\bigr).
}
\]
So the entire four-simplex boundary already carries one exact certified interval,
\[
\boxed{
\beta_Q^{\rm lo}
\le
\tau_{Q,*}^{\rm best,\partial}
\le
\beta_Q^{\rm hi},
}
\]
where `\(\tau_{Q,*}^{\rm best,\partial}\)` is the unknown true best closure time on the full boundary of the four-coordinate simplex.

---

## 3. Exact four-coordinate gradient-synergy theorem

The positive oriented initial slope magnitude on the four-simplex ray is
\[
\boxed{
k_{ijkl}(\mathbf a)
:=
-\sigma_0\,\nabla_\ell h(\boldsymbol\ell_\circ)\cdot \widehat{\mathbf s}_{ijkl}(\mathbf a)
=
a_i k_i+a_j k_j+a_k k_k+a_l k_l.
}
\]
So the oriented forward slope is
\[
K_{ijkl}(\mathbf a)=-k_{ijkl}(\mathbf a)<0.
\]

By the same exact Lagrange-multiplier argument as in Stage 210, the unique gradient-optimal interior point is
\[
\boxed{
\mathbf a_{ijkl}^{\rm grad}
=
\frac{(k_i,k_j,k_k,k_l)}{\sqrt{k_i^2+k_j^2+k_k^2+k_l^2}}.
}
\]
The exact maximum first-order slope magnitude is
\[
\boxed{
\max_{\mathbf a\in\Delta_{ijkl}^{+}} k_{ijkl}(\mathbf a)
=
\sqrt{k_i^2+k_j^2+k_k^2+k_l^2}.
}
\]
On the interior ratio patch `\(a_i>0\)`, the gradient-optimal ratios are
\[
\boxed{
r_{ijkl}^{\rm grad}=\frac{k_j}{k_i},
\qquad
s_{ijkl}^{\rm grad}=\frac{k_k}{k_i},
\qquad
t_{ijkl}^{\rm grad}=\frac{k_l}{k_i}.
}
\]

Because all four primitive slopes are positive, the interior four-coordinate gradient winner strictly beats every triple-face first-order maximum. For example,
\[
\bigl(k_{ijkl}^{\rm grad}\bigr)^2-igl(k_{ijk}^{\rm grad}\bigr)^2 = k_l^2>0,
\]
and similarly for the other three faces. Therefore:

### Exact four-coordinate gradient-synergy theorem

If all four primitive slopes are nonzero and monotone, a genuine interior four-coordinate mixed ray always beats every three-coordinate face at the **first-order** descent level.

That does **not** yet prove a four-coordinate winner after curvature enters, but it proves that the interior four-simplex carries genuinely new descent information already at linear order.

---

## 4. Exact four-coordinate curvature law and the total six-way cross-leverage theorem

Let the oriented symmetric Hessian block along the quadruple ray be
\[
H_{ijkl}(\mathbf a;\tau)
=
\begin{pmatrix}
 h_{ii} & h_{ij} & h_{ik} & h_{il}\\
 h_{ij} & h_{jj} & h_{jk} & h_{jl}\\
 h_{ik} & h_{jk} & h_{kk} & h_{kl}\\
 h_{il} & h_{jl} & h_{kl} & h_{ll}
\end{pmatrix}_{(ijkl)}.
\]
Then the exact second directional derivative is
\[
\boxed{
H_{1,ijkl}(\mathbf a;\tau)=\mathbf a^T H_{ijkl}(\mathbf a;\tau)\mathbf a.
}
\]
Expanding,
\[
\boxed{
\begin{aligned}
H_{1,ijkl}
={}&a_i^2 h_{ii}+a_j^2 h_{jj}+a_k^2 h_{kk}+a_l^2 h_{ll}\\
&+2a_i a_j h_{ij}+2a_i a_k h_{ik}+2a_i a_l h_{il}
+2a_j a_k h_{jk}+2a_j a_l h_{jl}+2a_k a_l h_{kl}.
\end{aligned}
}
\]
So the total off-diagonal leverage is carried by the exact weight
\[
\boxed{
w_{\Sigma}(\mathbf a)
:=
2\sum_{p<q} a_p a_q
=2(a_i a_j+a_i a_k+a_i a_l+a_j a_k+a_j a_l+a_k a_l).
}
\]
This weight admits the identity
\[
\boxed{
w_{\Sigma}(\mathbf a)=(a_i+a_j+a_k+a_l)^2-1,
}
\]
because `\(a_i^2+a_j^2+a_k^2+a_l^2=1\)` on the simplex.

Now use the exact Cauchy slack identity
\[
4(a_i^2+a_j^2+a_k^2+a_l^2)-(a_i+a_j+a_k+a_l)^2
=
\sum_{p<q}(a_p-a_q)^2\ge0.
\]
Since `\(a_i^2+a_j^2+a_k^2+a_l^2=1\)`, it follows that
\[
(a_i+a_j+a_k+a_l)^2\le 4,
\qquad
w_{\Sigma}(\mathbf a)\le 3.
\]
Equality holds iff
\[
a_i=a_j=a_k=a_l=\frac12.
\]
So the **equal-mix barycenter**
\[
\boxed{
\mathbf a_{ijkl}^{\rm eq}=\frac{(1,1,1,1)}{2}
}
\]
maximizes the total six-way off-diagonal leverage, with
\[
\boxed{
w_{\Sigma}(\mathbf a_{ijkl}^{\rm eq})=3.
}
\]
For comparison,

- the equal-mix point on any triple face has `\(w_{\Sigma}=2\)`,
- the equal-mix point on any pairwise edge has `\(w_{\Sigma}=1\)`.

So the four-coordinate interior barycenter adds one full new unit of total cross leverage beyond the triple-simplex barycenter.

### Exact diagonal-neutral reduction

If all off-diagonal Hessian entries vanish,
\[
h_{ij}=h_{ik}=h_{il}=h_{jk}=h_{jl}=h_{kl}=0,
\]
then
\[
\boxed{
H_{1,ijkl}(\mathbf a;\tau)
=
a_i^2 h_{ii}+a_j^2 h_{jj}+a_k^2 h_{kk}+a_l^2 h_{ll},
}
\]
so there is no genuine four-way curvature synergy. In that case the interior four-simplex carries only the first-order gradient gain from Section 3.

---

## 5. Exact fixed-simplex curvature envelopes and the certified local bracket

For either envelope label
\[
\star\in\{{\rm lo},{\rm hi}\},
\]
let the oriented symmetric `\(4\times 4\)` envelope block be
\[
H_{ijkl,\star}=
\begin{pmatrix}
 u_{ii,\star} & u_{ij,\star} & u_{ik,\star} & u_{il,\star}\\
 u_{ij,\star} & u_{jj,\star} & u_{jk,\star} & u_{jl,\star}\\
 u_{ik,\star} & u_{jk,\star} & u_{kk,\star} & u_{kl,\star}\\
 u_{il,\star} & u_{jl,\star} & u_{kl,\star} & u_{ll,\star}
\end{pmatrix}.
\]
For any fixed simplex point `\(\mathbf a\)`, define the exact envelope curvature scalar
\[
\boxed{
\kappa_{ijkl,\star}(\mathbf a):=
\mathbf a^T H_{ijkl,\star}\mathbf a.
}
\]
Then the exact certified local comparison root at that simplex point is
\[
\boxed{
\tau_{ijkl,\star}(\mathbf a)
:=
\mathcal T\bigl(H_0,-k_{ijkl}(\mathbf a);\kappa_{ijkl,\star}(\mathbf a)\bigr)
=
\frac{2H_0}{k_{ijkl}(\mathbf a)+\sqrt{k_{ijkl}(\mathbf a)^2-2H_0\kappa_{ijkl,\star}(\mathbf a)}}.
}
\]
So the four-coordinate simplex introduces **no new root algebra** at a fixed interior point. The only new difficulty is the genuine three-parameter interior search over `\(\mathbf a\)`.

### 5.1 Interior ratio patch `\(a_i>0\)`

Use the positive ratio coordinates
\[
r:=\frac{a_j}{a_i}>0,
\qquad
s:=\frac{a_k}{a_i}>0,
\qquad
t:=\frac{a_l}{a_i}>0,
\]
so that
\[
\boxed{
\mathbf a(r,s,t)
=
\frac{(1,r,s,t)}{\sqrt{1+r^2+s^2+t^2}}.
}
\]
Define the ten exact quadratic discriminant coefficients
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
D_\star:=2k_i k_l-4H_0 u_{il,\star},
\qquad
E_\star:=k_j^2-2H_0 u_{jj,\star},
\qquad
F_\star:=2k_j k_k-4H_0 u_{jk,\star},
}
\]
\[
\boxed{
G_\star:=2k_j k_l-4H_0 u_{jl,\star},
\qquad
H_\star:=k_k^2-2H_0 u_{kk,\star},
\qquad
I_\star:=2k_k k_l-4H_0 u_{kl,\star},
}
\]
\[
\boxed{
J_\star:=k_l^2-2H_0 u_{ll,\star}.
}
\]
Then the exact discriminant numerator is
\[
\boxed{
\Delta^{\sharp}_{ijkl,\star}(r,s,t)
:=
A_\star+B_\star r+C_\star s+D_\star t
+E_\star r^2+F_\star rs+G_\star rt+H_\star s^2+I_\star st+J_\star t^2.
}
\]
The certified interior bracket becomes
\[
\boxed{
\tau_{ijkl,\star}(r,s,t)
=
\frac{2H_0\sqrt{1+r^2+s^2+t^2}}
{k_i+r k_j+s k_k+t k_l+\sqrt{\Delta^{\sharp}_{ijkl,\star}(r,s,t)}}.
}
\]
Restricting to any face recovers the exact Stage 211 triple formulas.

---

## 6. The canonical quadruple-screen audit

Stage 213 does **not** yet solve the full three-parameter interior optimizer.
The smallest exact screen set is instead:

1. the four exact imported full-face intervals from Stage 212,
   \[
   \mathcal I_{ijk}^{\triangle},
   \qquad
   \mathcal I_{ijl}^{\triangle},
   \qquad
   \mathcal I_{ikl}^{\triangle},
   \qquad
   \mathcal I_{jkl}^{\triangle},
   \]
2. the unique interior gradient-optimal screen point
   \[
   \mathbf a_{ijkl}^{\rm grad},
   \]
3. the unique interior equal-mix barycentric screen point
   \[
   \mathbf a_{ijkl}^{\rm eq}.
   \]

So the exact Stage 213 quadruple-screen packet is
\[
\boxed{
\mathcal S_{ijkl}^{\rm quad}
:=
\Bigl(
\mathcal I_{ijk}^{\triangle},
\mathcal I_{ijl}^{\triangle},
\mathcal I_{ikl}^{\triangle},
\mathcal I_{jkl}^{\triangle},
\mathbf a_{ijkl}^{\rm grad},
\mathbf a_{ijkl}^{\rm eq}
\Bigr).
}
\]
This is the exact four-coordinate analogue of the Stage 210 canonical triple-screen packet, but now the whole codimension-one boundary is already fully optimized and imported rather than re-audited.

---

## 7. Exact support-cardinality-4 theorem gate

Let
\[
\tau_{ijkl,\rm hi}^{\rm grad}:=\tau_{ijkl,\rm hi}(\mathbf a_{ijkl}^{\rm grad}),
\qquad
\tau_{ijkl,\rm hi}^{\rm eq}:=\tau_{ijkl,\rm hi}(\mathbf a_{ijkl}^{\rm eq}),
\]
and similarly for the lower certified screen values.
Then the first exact four-coordinate certification rule is immediate.

### Exact interior-screen dominance theorem

If either canonical interior screen satisfies
\[
\boxed{
\tau_{ijkl,\rm hi}^{\rm can}
<
\beta_Q^{\rm lo},
\qquad
\text{can}\in\{{\rm grad},{\rm eq}\},
}
\]
then there exists a **genuine interior four-coordinate mixed ray** whose certified closure time lies strictly below every already-solved triple-face boundary winner.

The proof is the same as in Stage 210:

- the actual root at the chosen interior screen point is bounded above by its certified upper bracket,
- the actual face winners are bounded below by the imported certified lower face brackets,
- so the interior screen point already certifies a strict interior four-coordinate improvement over the full Stage 212 boundary ledger.

### Exact canonical non-improvement filter

Conversely, if
\[
\boxed{
\min\bigl(
\tau_{ijkl,\rm lo}^{\rm grad},
\tau_{ijkl,\rm lo}^{\rm eq}
\bigr)
>
\beta_Q^{\rm hi},
}
\]
then neither canonical four-coordinate screen beats the best Stage 212 boundary winner on that simplex.

That is **not** a full no-go theorem for the four-simplex interior, because the genuine three-parameter interior optimizer has not yet been solved. But it is the first exact filter that can rule out the two canonical four-way interior mechanisms before the full interior search is attempted.

---

## 8. Global support-cardinality-4 gate against the Stage 212 up-to-three-coordinate ledger

Carry forward the exact Stage 212 certified interval for the whole local support-`<=3` search,
\[
\boxed{
\tau_{\le 3,\min}^{\rm lo}
\le
\tau_{\le 3,*}^{\rm best}
\le
\tau_{\le 3,\min}^{\rm hi}.
}
\]
Now compare any primitive quadruple canonical interior screen against this already-finished ledger.

### Exact support-cardinality-4 improvement gate

If for some primitive quadruple `\(Q=(i,j,k,l)\)` and some canonical interior screen `\(\text{can}\in\{{\rm grad},{\rm eq}\}\)` one has
\[
\boxed{
\tau_{ijkl,\rm hi}^{\rm can}
<
\tau_{\le 3,\min}^{\rm lo},
}
\]
then there exists a genuine support-cardinality-`4` mixed ray whose certified closure time lies strictly below **every** already-ranked support-`<=3` ray.

### Exact canonical support-cardinality-4 non-improvement filter

If for **every** primitive quadruple `\(Q\)` one has
\[
\boxed{
\min\bigl(
\tau_{ijkl,\rm lo}^{\rm grad},
\tau_{ijkl,\rm lo}^{\rm eq}
\bigr)
>
\tau_{\le 3,\min}^{\rm hi},
}
\]
then no canonical four-way screen beats the current support-`<=3` winner.

Again, this is a filter on the canonical four-way screens, not yet a full no-go theorem for the entire four-coordinate interior.

So after Stage 213 the support-cardinality-`4` question is no longer opaque.
The first honest gate is now a direct interval comparison against the already-finished Stage 212 sieve.

---

## 9. Minimal packet for the next stage

After Stage 213, the full boundary of every primitive four-coordinate simplex is already solved, and the only unresolved part is the genuinely new three-parameter interior optimizer.

The smallest exact packet for that next stage is
\[
\boxed{
\mathcal P_{ijkl}^{(4)}
:=
\Bigl(
H_0,
(k_i,k_j,k_k,k_l),
H_{ijkl,\rm lo},
H_{ijkl,\rm hi},
T_{ijkl}(\mathbf a)
\Bigr),
}
\]
where

- `\(H_0\)` is the oriented logarithmic defect,
- `\((k_i,k_j,k_k,k_l)\)` are the four primitive slope magnitudes,
- `\(H_{ijkl,\rm lo},H_{ijkl,\rm hi}\)` are the lower/upper oriented `4 x 4` Hessian-envelope blocks,
- `\(T_{ijkl}(\mathbf a)\)` is the local validity-radius map on the four-simplex patch.

The natural continuation is now completely sharp:

> Stage 214 should solve the full interior four-simplex optimizer by reducing the three-parameter stationary system to a finite algebraic candidate set, just as Stage 211 did for the three-coordinate simplex interior.

---

## 10. Best current reading after Stage 213

Stage 212 already proved that the support-`<=3` search is finite and exact.
Stage 213 now shows what is genuinely new at the first four-coordinate level:

1. every codimension-one face of the four-coordinate simplex is already one of the Stage 212 primitive triples,
2. a genuine interior four-coordinate ray always beats every triple face at the first-order gradient level,
3. the equal-mix four-way barycenter uniquely maximizes the total six-way off-diagonal Hessian leverage,
4. every fixed interior four-simplex point already carries an exact certified local bracket,
5. and the smallest honest support-cardinality-`4` audit is now a six-row screen:
   - four imported triple-face rows,
   - one interior gradient-optimal row,
   - one interior equal-mix row.

So the next theorem gate is no longer “is there any four-coordinate effect at all?”
There is.
The next real question is whether the full three-parameter interior optimizer on a primitive four-simplex can produce a certified winner that beats both canonical interior screens and the now-finished support-`<=3` ledger.
