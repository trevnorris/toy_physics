# Moving-Throat PDE — Stage 216: Unique Five-Coordinate Positive Simplex, Exact Face Reduction, and the Support-Cardinality-5 Theorem Gate

## Status

**Exact within the carried Stage 215 support-`<=4` certified ledger, once a local oriented `5 x 5` Hessian-envelope block and the corresponding validity map are supplied on the chosen free-quintuple patch.**

This stage does **not** yet solve the full four-parameter interior optimizer on the unique five-coordinate simplex.
It is the exact continuation *between* the Stage 215 global support-`<=4` splice and the Stage 217 interior five-coordinate candidate reduction: the first place where a genuinely interior support-cardinality-`5` mixed ray can appear and be audited without reopening the already-finished quadruple boundary problem.

---

## Purpose

Stage 215 closed the full support-`<=4` local search.
Every primitive quadruple already carries a closed-simplex certified interval, and the whole local sieve through four active primitive coordinates has been reduced to one finite certified ledger.

That makes the next theorem gate completely sharp:

> can a **genuine interior five-coordinate** mixed branch beat the already-finished support-`<=4` ledger, and what is the smallest exact audit that can test that before solving the full four-parameter interior optimizer?

This stage answers that.

The main outputs are:

1. the exact **unique positive spherical five-simplex** and its full face reduction back to the Stage 215 primitive-quadruple packets,
2. the exact **five-coordinate gradient-synergy theorem** and the unique interior gradient-optimal ray,
3. the exact **five-coordinate curvature law** and the theorem that the equal-mix barycenter maximizes the total ten-way off-diagonal leverage,
4. the exact **fixed-simplex certified bracket** for any admissible interior point,
5. the exact **canonical five-way screen audit** consisting of five imported quadruple faces plus two interior canonical rows,
6. the exact **support-cardinality-`5` theorem gate** that certifies a genuine five-coordinate improvement over the entire support-`<=4` ledger whenever one canonical interior screen already wins,
7. and the exact statement that there are **no higher-support local mixed rays** beyond this stage because the free-quintuple has only five primitive axes.

So Stage 216 is the five-coordinate analogue of Stage 213, but now the whole codimension-one boundary is already fully optimized and imported from Stage 215, and the only genuinely new local content is the interior of the unique five-coordinate simplex.

---

## 1. The unique positive spherical five-simplex and exact face reduction

Keep the oriented primitive free axes
\[
\widehat{\mathbf e}_\lambda,\quad
\widehat{\mathbf e}_c,\quad
\widehat{\mathbf e}_\gamma,\quad
\widehat{\mathbf e}_U,\quad
\widehat{\mathbf e}_W,
\]
with positive oriented slope magnitudes
\[
k_\lambda,\quad k_c,\quad k_\gamma,\quad k_U,\quad k_W>0.
\]

The unique support-cardinality-`5` search set is the positive spherical five-simplex
\[
\boxed{
\Delta_5^+
:=
\Bigl\{
\mathbf a=(a_\lambda,a_c,a_\gamma,a_U,a_W)\in\mathbb R_{\ge 0}^5:
 a_\lambda^2+a_c^2+a_\gamma^2+a_U^2+a_W^2=1
\Bigr\}.
}
\]
Each simplex point generates the oriented mixed ray
\[
\boxed{
\widehat{\mathbf s}_5(\mathbf a)
=
a_\lambda\widehat{\mathbf e}_\lambda
+a_c\widehat{\mathbf e}_c
+a_\gamma\widehat{\mathbf e}_\gamma
+a_U\widehat{\mathbf e}_U
+a_W\widehat{\mathbf e}_W.
}
\]

Its five codimension-one faces are exactly the five Stage 215 primitive quadruple simplices:
\[
Q_{\widehat\lambda}=\{c,\gamma,U,W\},
\qquad
Q_{\widehat c}=\{\lambda,\gamma,U,W\},
\]
\[
Q_{\widehat\gamma}=\{\lambda,c,U,W\},
\qquad
Q_{\widehat U}=\{\lambda,c,\gamma,W\},
\qquad
Q_{\widehat W}=\{\lambda,c,\gamma,U\}.
\]
So the entire codimension-one boundary is already closed by the Stage 215 ledger.

Let the imported Stage 215 closed-simplex intervals be
\[
\mathcal I_{\widehat\lambda}^{\square},\quad
\mathcal I_{\widehat c}^{\square},\quad
\mathcal I_{\widehat\gamma}^{\square},\quad
\mathcal I_{\widehat U}^{\square},\quad
\mathcal I_{\widehat W}^{\square},
\]
with lower and upper face minima
\[
\boxed{
\beta_5^{\rm lo}
:=
\min\bigl(
\tau_{\widehat\lambda,\min}^{\rm lo,\square},
\tau_{\widehat c,\min}^{\rm lo,\square},
\tau_{\widehat\gamma,\min}^{\rm lo,\square},
\tau_{\widehat U,\min}^{\rm lo,\square},
\tau_{\widehat W,\min}^{\rm lo,\square}
\bigr),
}
\]
\[
\boxed{
\beta_5^{\rm hi}
:=
\min\bigl(
\tau_{\widehat\lambda,\min}^{\rm hi,\square},
\tau_{\widehat c,\min}^{\rm hi,\square},
\tau_{\widehat\gamma,\min}^{\rm hi,\square},
\tau_{\widehat U,\min}^{\rm hi,\square},
\tau_{\widehat W,\min}^{\rm hi,\square}
\bigr).
}
\]
Therefore the full boundary of the unique five-coordinate simplex already carries one exact certified interval,
\[
\boxed{
\beta_5^{\rm lo}
\le
\tau_{5,*}^{\rm best,\partial}
\le
\beta_5^{\rm hi},
}
\]
where `\(\tau_{5,*}^{\rm best,\partial}\)` is the unknown true best closure time on the full codimension-one boundary of `\(\Delta_5^+\)`.

So from this point onward, the only genuinely new content is the **interior** of the unique five-coordinate simplex.

---

## 2. Exact five-coordinate gradient-synergy theorem

The positive oriented initial slope magnitude on the five-simplex ray is
\[
\boxed{
k_5(\mathbf a)
:=
-\sigma_0\,\nabla_\ell h(\boldsymbol\ell_\circ)\cdot \widehat{\mathbf s}_5(\mathbf a)
=
a_\lambda k_\lambda+a_c k_c+a_\gamma k_\gamma+a_U k_U+a_W k_W.
}
\]
So the oriented forward slope is
\[
K_5(\mathbf a)=-k_5(\mathbf a)<0.
\]

By the exact Cauchy/Lagrange-multiplier argument used in the earlier support-cardinality stages, the unique gradient-optimal interior point is
\[
\boxed{
\mathbf a_5^{\rm grad}
=
\frac{(k_\lambda,k_c,k_\gamma,k_U,k_W)}{\sqrt{k_\lambda^2+k_c^2+k_\gamma^2+k_U^2+k_W^2}}.
}
\]
The exact maximum first-order slope magnitude is
\[
\boxed{
\max_{\mathbf a\in\Delta_5^+} k_5(\mathbf a)
=
\sqrt{k_\lambda^2+k_c^2+k_\gamma^2+k_U^2+k_W^2}.
}
\]
On the interior ratio patch `\(a_\lambda>0\)`, the gradient-optimal ratios are
\[
\boxed{
r_5^{\rm grad}=\frac{k_c}{k_\lambda},
\qquad
s_5^{\rm grad}=\frac{k_\gamma}{k_\lambda},
\qquad
t_5^{\rm grad}=\frac{k_U}{k_\lambda},
\qquad
u_5^{\rm grad}=\frac{k_W}{k_\lambda}.
}
\]

Because all five primitive slopes are positive, the interior five-coordinate gradient winner strictly beats every quadruple-face first-order maximum. For example,
\[
\bigl(k_5^{\rm grad}\bigr)^2-
\bigl(k_{Q_{\widehat\lambda}}^{\rm grad}\bigr)^2
=
k_\lambda^2>0,
\]
and similarly for the other four faces.
Therefore:

### Exact five-coordinate gradient-synergy theorem

If all five primitive slopes are nonzero and monotone, a genuine interior five-coordinate mixed ray always beats every four-coordinate face at the **first-order** descent level.

That does **not** yet prove a five-coordinate winner after curvature enters, but it proves that the interior five-simplex carries genuinely new descent information already at linear order.

---

## 3. Exact five-coordinate curvature law and the total ten-way cross-leverage theorem

Let the oriented symmetric Hessian block along the five-coordinate ray be
\[
H_5(\mathbf a;\tau)=\bigl(h_{pq}(\tau)\bigr)_{p,q\in\{\lambda,c,\gamma,U,W\}}.
\]
Then the exact second directional derivative is
\[
\boxed{
H_{1,5}(\mathbf a;\tau)=\mathbf a^T H_5(\mathbf a;\tau)\mathbf a.
}
\]
Expanding,
\[
\boxed{
H_{1,5}
=
\sum_p a_p^2 h_{pp}
+
2\sum_{p<q} a_p a_q h_{pq}.
}
\]
So the total off-diagonal leverage is carried by the exact weight
\[
\boxed{
w_\Sigma(\mathbf a)
:=
2\sum_{p<q} a_p a_q
=
2(a_\lambda a_c+a_\lambda a_\gamma+a_\lambda a_U+a_\lambda a_W+a_c a_\gamma+a_c a_U+a_c a_W+a_\gamma a_U+a_\gamma a_W+a_U a_W).
}
\]
This weight admits the identity
\[
\boxed{
w_\Sigma(\mathbf a)=\Bigl(\sum_p a_p\Bigr)^2-1,
}
\]
because `\(\sum_p a_p^2=1\)` on the simplex.

Now use the exact Cauchy slack identity
\[
5\sum_p a_p^2-\Bigl(\sum_p a_p\Bigr)^2
=
\sum_{p<q}(a_p-a_q)^2\ge 0.
\]
Since `\(\sum_p a_p^2=1\)`, it follows that
\[
\Bigl(\sum_p a_p\Bigr)^2\le 5,
\qquad
w_\Sigma(\mathbf a)\le 4.
\]
Equality holds iff
\[
a_\lambda=a_c=a_\gamma=a_U=a_W=\frac{1}{\sqrt5}.
\]
So the **equal-mix barycenter**
\[
\boxed{
\mathbf a_5^{\rm eq}=\frac{(1,1,1,1,1)}{\sqrt5}
}
\]
maximizes the total ten-way off-diagonal leverage, with
\[
\boxed{
w_\Sigma(\mathbf a_5^{\rm eq})=4.
}
\]
For comparison,

- the equal-mix point on any quadruple face has `\(w_\Sigma=3\)`,
- the equal-mix point on any triple face has `\(w_\Sigma=2\)`,
- the equal-mix point on any pairwise edge has `\(w_\Sigma=1\)`.

So the five-coordinate interior barycenter adds one full new unit of total cross leverage beyond the quadruple-simplex barycenter.

### Exact diagonal-neutral reduction

If all off-diagonal Hessian entries vanish,
\[
h_{pq}=0\qquad(p\neq q),
\]
then
\[
\boxed{
H_{1,5}(\mathbf a;\tau)=\sum_p a_p^2 h_{pp},
}
\]
so there is no genuine five-way curvature synergy. In that case the interior five-simplex carries only the first-order gradient gain from Section 2.

---

## 4. Exact fixed-simplex curvature envelopes and the certified local bracket

For either envelope label
\[
\star\in\{{\rm lo},{\rm hi}\},
\]
let the oriented symmetric `\(5\times 5\)` envelope block be
\[
H_{5,\star}=\bigl(u_{pq,\star}\bigr)_{p,q\in\{\lambda,c,\gamma,U,W\}}.
\]
For any fixed simplex point `\(\mathbf a\)`, define the exact envelope curvature scalar
\[
\boxed{
\kappa_{5,\star}(\mathbf a):=\mathbf a^T H_{5,\star}\mathbf a.
}
\]
Then the exact certified local comparison root at that simplex point is
\[
\boxed{
\tau_{5,\star}(\mathbf a)
:=
\mathcal T\bigl(H_0,-k_5(\mathbf a);\kappa_{5,\star}(\mathbf a)\bigr)
=
\frac{2H_0}{k_5(\mathbf a)+\sqrt{k_5(\mathbf a)^2-2H_0\kappa_{5,\star}(\mathbf a)}}.
}
\]
So the five-coordinate simplex introduces **no new root algebra** at a fixed interior point. The only new difficulty is the genuine four-parameter interior search over `\(\mathbf a\)`.

That full four-parameter interior optimizer is deferred deliberately to Stage 217.
Stage 216 keeps only the first exact gate needed before that solve: the full boundary splice plus the two canonical interior screen points.

---

## 5. The canonical five-way screen audit

Stage 216 does **not** yet solve the full four-parameter interior optimizer.
The smallest exact screen set is instead:

1. the five exact imported full-face intervals from Stage 215,
   \[
   \mathcal I_{\widehat\lambda}^{\square},\quad
   \mathcal I_{\widehat c}^{\square},\quad
   \mathcal I_{\widehat\gamma}^{\square},\quad
   \mathcal I_{\widehat U}^{\square},\quad
   \mathcal I_{\widehat W}^{\square},
   \]
2. the unique interior gradient-optimal screen point
   \[
   \mathbf a_5^{\rm grad},
   \]
3. the unique interior equal-mix barycentric screen point
   \[
   \mathbf a_5^{\rm eq}.
   \]

So the exact Stage 216 five-way screen packet is
\[
\boxed{
\mathcal S_5^{\rm can}
:=
\Bigl(
\mathcal I_{\widehat\lambda}^{\square},
\mathcal I_{\widehat c}^{\square},
\mathcal I_{\widehat\gamma}^{\square},
\mathcal I_{\widehat U}^{\square},
\mathcal I_{\widehat W}^{\square},
\mathbf a_5^{\rm grad},
\mathbf a_5^{\rm eq}
\Bigr).
}
\]
This is the exact five-coordinate analogue of the Stage 213 canonical quadruple-screen packet, but now there is only **one** support-cardinality-`5` simplex, and its whole codimension-one boundary is already fully optimized and imported rather than re-audited.

---

## 6. Exact support-cardinality-`5` theorem gate

Let
\[
\tau_{5,\rm hi}^{\rm grad}:=\tau_{5,\rm hi}(\mathbf a_5^{\rm grad}),
\qquad
\tau_{5,\rm hi}^{\rm eq}:=\tau_{5,\rm hi}(\mathbf a_5^{\rm eq}),
\]
and similarly for the lower certified screen values.
Then the first exact five-coordinate certification rule is immediate.

### Exact interior-screen dominance theorem

If either canonical interior screen satisfies
\[
\boxed{
\tau_{5,\rm hi}^{\rm can}
<
\beta_5^{\rm lo},
\qquad
\text{can}\in\{{\rm grad},{\rm eq}\},
}
\]
then there exists a **genuine interior five-coordinate mixed ray** whose certified closure time lies strictly below every already-solved four-coordinate boundary winner on the unique five-simplex.

The proof is the same as in Stage 213:

- the actual root at the chosen interior screen point is bounded above by its certified upper bracket,
- the actual face winners are bounded below by the imported certified lower face brackets,
- so the interior screen point already certifies a strict interior five-coordinate improvement over the full support-`<=4` boundary ledger.

### Exact canonical non-improvement filter

Conversely, if
\[
\boxed{
\min\bigl(
\tau_{5,\rm lo}^{\rm grad},
\tau_{5,\rm lo}^{\rm eq}
\bigr)
>
\beta_5^{\rm hi},
}
\]
then neither canonical five-coordinate screen beats the best support-`<=4` boundary winner on the closed five-simplex.

That is **not** a full no-go theorem for the five-simplex interior, because the genuine four-parameter interior optimizer has not yet been solved. But it is the first exact filter that can rule out the two canonical five-way interior mechanisms before the full Stage 217 interior search is attempted.

---

## 7. Global support-cardinality-`5` gate against the Stage 215 support-`<=4` ledger

Carry forward the exact Stage 215 certified interval for the whole local support-`<=4` search,
\[
\boxed{
\tau_{\le 4,\min}^{\rm lo}
\le
\tau_{\le 4,*}^{\rm best}
\le
\tau_{\le 4,\min}^{\rm hi}.
}
\]
Now compare the unique five-coordinate canonical interior screens against this already-finished ledger.

### Exact support-cardinality-`5` improvement gate

If for some canonical interior screen `\(\text{can}\in\{{\rm grad},{\rm eq}\}\)` one has
\[
\boxed{
\tau_{5,\rm hi}^{\rm can}
<
\tau_{\le 4,\min}^{\rm lo},
}
\]
then there exists a genuine support-cardinality-`5` mixed ray whose certified closure time lies strictly below **every** already-ranked support-`<=4` ray.

### Exact canonical support-cardinality-`5` non-improvement filter

If
\[
\boxed{
\min\bigl(
\tau_{5,\rm lo}^{\rm grad},
\tau_{5,\rm lo}^{\rm eq}
\bigr)
>
\tau_{\le 4,\min}^{\rm hi},
}
\]
then neither canonical five-coordinate screen beats the already-finished global support-`<=4` winner.

Again, that is not yet a complete no-go theorem for the five-coordinate interior, because the full four-parameter interior optimizer has not yet been solved. But it is the exact final **support-cardinality gate** before the unique interior optimizer of Stage 217.

---

## 8. Exact support-cardinality ceiling theorem

Because the free log space itself is the free quintuple
\[
(\lambda,c,\gamma,U,W),
\]
no local mixed ray can activate more than these five primitive coordinates at once.
Therefore:

\[
\boxed{
\text{there are no support-cardinality-`>5` local mixed rays in the free-quintuple search.}
}
\]

So Stage 216 is the **last combinatorial support-cardinality gate** in the local mixed-ray ladder.
What remains after it is not another support class, but only the unique four-parameter interior optimizer on `\(\Delta_5^+\)`.

---

## 9. What Stage 216 changes in the program

This stage does three things.

### 9.1 It identifies the exact missing bridge between Stage 215 and Stage 217

Stage 215 closed the support-`<=4` ledger.
Stage 217 solves the unique interior five-coordinate optimizer.
Stage 216 is the missing exact bridge between them: it defines the unique five-coordinate simplex, proves that its full codimension-one boundary is already closed, and isolates the two canonical five-way interior screens before the interior algebra is attacked.

### 9.2 It shows that the last genuinely new local content is interior

Nothing on the codimension-one boundary is new anymore.
Every boundary face is already a Stage 215 quadruple packet.
So the only new support-cardinality-`5` question is whether the unique four-parameter interior optimizer can beat the imported support-`<=4` ledger.

### 9.3 It makes the final search ladder conceptually complete

After this note there is no higher-support combinatorial stage left to build.
The whole local mixed-ray ladder is now:

- support-`<=3` closed by Stage 212,
- support-`4` globally ranked by Stage 215,
- support-`5` gated here at Stage 216,
- and the unique support-`5` interior optimizer solved next in Stage 217.

---

## 10. Immediate next derivation step

The next exact move is now completely sharp:

1. keep the imported Stage 215 five-face boundary ledger from this note,
2. move to the interior ratio chart on `\(\Delta_5^+\)`,
3. derive the full four-parameter stationary equations,
4. and reduce the unique interior five-coordinate optimizer to a finite algebraic candidate set.

That is exactly the Stage 217 task.
