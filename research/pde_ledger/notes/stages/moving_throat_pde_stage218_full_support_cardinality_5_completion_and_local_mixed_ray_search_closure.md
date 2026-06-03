# Moving-Throat PDE — Stage 218: Full Support-`<=5` Completion, Exact Boundary Identification, and the Final Local Mixed-Ray Closure Theorem

## Status

**Exact within the carried Stage 215 support-`<=4` certified ledger and the Stage 217 unique five-coordinate interior candidate reduction, once the admissible interior five-coordinate packets are supplied on the chosen free-quintuple patch.**

This stage does **not** introduce a new constitutive law or a new optimizer.
It splices the unique support-cardinality-`5` interior packet back into the already-finished support-`<=4` ledger and closes the **entire local mixed-ray search** on the positive free-quintuple simplex.

---

## Purpose

Stage 215 already finished the global support-`<=4` search: every primitive quadruple carries a closed-simplex certified interval, and the whole local search through four active primitive coordinates is finite and certified.

Stage 217 then solved the only genuinely new support-cardinality-`5` problem: the unique interior five-coordinate mixed simplex is reduced to a finite algebraic candidate set, with an exact certified interval for the best admissible interior winner.

That leaves one last exact continuation point:

> how do we splice the unique support-`5` interior packet to the already-finished support-`<=4` ledger so that the **entire** local mixed-ray search over the free quintuple is closed once and for all?

This stage answers that.

The main outputs are:

1. the exact **boundary-identification theorem** saying that the full boundary of the unique five-coordinate simplex is precisely the already-solved support-`<=4` search set,
2. the exact **support-cardinality ceiling theorem** showing that there are no higher-support local branches beyond `5`,
3. the exact **support-`<=5` splice theorem** for the global best local closure time,
4. the exact **support-`5` improvement / non-improvement theorems** against the imported support-`<=4` ledger,
5. the exact statement that any remaining ambiguity is **finite-candidate ambiguity**, not a new continuum scan,
6. and the exact **final evaluation-budget theorem** for the full local mixed-ray search.

So Stage 218 is the closure of the whole local mixed-ray sieve: after this note there are no more combinatorial support-cardinality stages left to build.

Script-backed status:
- `scripts/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_sympy_audit.py`
  checks the boundary-stratification counts, support-ceiling logic, final splice
  interval theorems, and carried budget formulas.
- `mathematica/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_mathematica_audit.wl`
  mirrors the same closure algebra in the second CAS and keeps the carried
  Stage 215 and Stage 217 budget constants explicit rather than hidden.

---

## 1. Carry-forward global packets

Let the free primitive index set be
\[
\boxed{
\mathfrak I_5:=\{\lambda,c,\gamma,U,W\}.
}
\]
There is only one positive spherical five-simplex on this free quintuple,
\[
\boxed{
\Delta_5^+
:=
\Bigl\{
\mathbf a=(a_\lambda,a_c,a_\gamma,a_U,a_W)\in\mathbb R_{\ge 0}^5:
\ a_\lambda^2+a_c^2+a_\gamma^2+a_U^2+a_W^2=1
\Bigr\}.
}
\]

### 1.1 Imported Stage 215 support-`<=4` global ledger

Carry forward the already-closed support-`<=4` certified interval
\[
\boxed{
\tau_{\le 4,\min}^{\rm lo}
\le
\tau_{\le 4,*}^{\rm best}
\le
\tau_{\le 4,\min}^{\rm hi}.
}
\]
This packet already includes

- every primitive ray,
- every optimized pairwise cone,
- every closed primitive triple simplex,
- and every closed primitive quadruple simplex.

Equivalently, using the Stage 217 notation for the five quadruple-face splice,
\[
\boxed{
\tau_{\le 4,\min}^{\rm lo}=\beta_5^{\rm lo},
\qquad
\tau_{\le 4,\min}^{\rm hi}=\beta_5^{\rm hi}.
}
\]

### 1.2 Imported Stage 217 support-`5` interior packet

Carry forward the exact interior five-coordinate certified interval
\[
\boxed{
\tau_{5,\min}^{\rm lo,int}
\le
\tau_{5,*}^{\rm best,int}
\le
\tau_{5,\min}^{\rm hi,int}.
}
\]
Here `\(\tau_{5,*}^{\rm best,int}\)` is the unknown true best admissible closure time over the **interior** of `\(\Delta_5^+\)`, while its whole boundary is already absorbed into the support-`<=4` ledger above.

So the full remaining local search content is now only the pair
\[
\boxed{
\bigl(
[\tau_{\le 4,\min}^{\rm lo},\tau_{\le 4,\min}^{\rm hi}],
[\tau_{5,\min}^{\rm lo,int},\tau_{5,\min}^{\rm hi,int}]
\bigr).
}
\]

---

## 2. Exact boundary identification theorem

The unique five-coordinate simplex `\(\Delta_5^+\)` has the usual positive-support face stratification.
Its codimension-one faces are the five primitive quadruple simplices
\[
Q_{\widehat\lambda}=\{c,\gamma,U,W\},
\quad
Q_{\widehat c}=\{\lambda,\gamma,U,W\},
\quad
Q_{\widehat\gamma}=\{\lambda,c,U,W\},
\]
\[
Q_{\widehat U}=\{\lambda,c,\gamma,W\},
\quad
Q_{\widehat W}=\{\lambda,c,\gamma,U\}.
\]
Their intersections generate all lower-support strata:

- `5` support-`4` facets,
- `10` support-`3` ridges,
- `10` support-`2` edges,
- `5` support-`1` vertices.

So the number of nonempty proper support strata is
\[
\boxed{
5+10+10+5=30=2^5-2.
}
\]
Every nonempty proper support set in `\mathfrak I_5` is contained in at least one codimension-one quadruple face, and in fact a support subset of size `k` belongs to exactly `5-k` quadruple faces.

### Exact boundary-identification theorem

Let `\(\mathcal S_{\le 4}^{\rm loc}\)` denote the already-solved local search set with support-cardinality at most `4` on the free quintuple. Then
\[
\boxed{
\partial\Delta_5^+ = \mathcal S_{\le 4}^{\rm loc}.
}
\]
Equivalently, the whole boundary of the unique support-`5` simplex is exactly the already-closed support-`<=4` search domain.

This theorem is stronger than a mere combinatorial count: because every quadruple face is carried as a **closed** simplex from Stage 215, all of its triple, pair, and ray subfaces are already included automatically.

---

## 3. Exact support-cardinality ceiling theorem

There are only five primitive free axes in `\mathfrak I_5`. Therefore every local mixed ray on the free quintuple has support-cardinality
\[
1\le \#\operatorname{supp}(\mathbf a) \le 5.
\]
There is no support-cardinality-`6` or higher local branch to consider.

### Exact support ceiling theorem

Let `\(\tau_{\le 5,*}^{\rm best}\)` be the true best local closure time on the full closed simplex `\(\Delta_5^+\)`. Then
\[
\boxed{
\tau_{\le 5,*}^{\rm best}
=
\min\bigl(
\tau_{\le 4,*}^{\rm best},
\tau_{5,*}^{\rm best,int}
\bigr).
}
\]
So after Stage 217 there is literally no further support-cardinality family left to add: the full local search is the support-`<=4` boundary winner versus the unique support-`5` interior winner.

---

## 4. Exact support-`<=5` splice theorem

Define the final support-`<=5` certified interval by
\[
\boxed{
\tau_{\le 5,\min}^{\rm lo}
:=
\min\bigl(
\tau_{\le 4,\min}^{\rm lo},
\tau_{5,\min}^{\rm lo,int}
\bigr),
}
\]
\[
\boxed{
\tau_{\le 5,\min}^{\rm hi}
:=
\min\bigl(
\tau_{\le 4,\min}^{\rm hi},
\tau_{5,\min}^{\rm hi,int}
\bigr).
}
\]

### Exact support-`<=5` splice theorem

The true best local closure time on the full free-quintuple simplex obeys
\[
\boxed{
\tau_{\le 5,\min}^{\rm lo}
\le
\tau_{\le 5,*}^{\rm best}
\le
\tau_{\le 5,\min}^{\rm hi}.
}
\]

The proof is now immediate from the previous section:
\[
\tau_{\le 5,*}^{\rm best}
=
\min\bigl(
\tau_{\le 4,*}^{\rm best},
\tau_{5,*}^{\rm best,int}
\bigr),
\]
and each of the two entries on the right-hand side already carries an exact certified interval.

So the entire local mixed-ray search is reduced to **one last interval splice**.

---

## 5. Exact support-`5` classification and winner theorems

Because there is only one support-`5` simplex, its contribution to the full search falls into three exact classes.

### 5.1 Genuine support-`5` improvement theorem

If the support-`5` interior upper bound already lies below the support-`<=4` lower bound,
\[
\boxed{
\tau_{5,\min}^{\rm hi,int}<\tau_{\le 4,\min}^{\rm lo},
}
\]
then every admissible support-`5` interior winner beats every admissible support-`<=4` winner. Therefore
\[
\boxed{
\tau_{\le 5,*}^{\rm best}=\tau_{5,*}^{\rm best,int},
}
\]
and the full local search has a genuine five-coordinate interior winner.

### 5.2 No genuine support-`5` improvement theorem

If the support-`5` interior lower bound already lies above the support-`<=4` upper bound,
\[
\boxed{
\tau_{5,\min}^{\rm lo,int}>\tau_{\le 4,\min}^{\rm hi},
}
\]
then no admissible support-`5` interior candidate beats the already-finished support-`<=4` ledger. Therefore
\[
\boxed{
\tau_{\le 5,*}^{\rm best}=\tau_{\le 4,*}^{\rm best},
}
\]
and the full local search closes without any genuine five-coordinate improvement.

### 5.3 Ambiguous but finite theorem gate

If instead the intervals overlap,
\[
\tau_{5,\min}^{\rm lo,int}\le \tau_{\le 4,\min}^{\rm hi}
\qquad\text{and}\qquad
\tau_{5,\min}^{\rm hi,int}\ge \tau_{\le 4,\min}^{\rm lo},
\]
then the support-`5` interior remains **certifiedly ambiguous at the interval level**.

But this is no longer a continuum ambiguity. Stage 217 already reduced the support-`5` interior to a finite admissible algebraic candidate set. So even in the ambiguous case, the only unresolved work is finite candidate comparison, not a new free multi-parameter scan.

This is the exact sense in which Stage 218 closes the search sieve even before the actual PDE data are inserted.

---

## 6. Exact final local mixed-ray closure theorem

Collect the carried finite ledgers:

- Stage 209: optimized primitive pair cones,
- Stage 212: full primitive triple closed-simplex packets,
- Stage 215: full primitive quadruple closed-simplex packets and support-`<=4` global ledger,
- Stage 217: unique support-`5` interior candidate packet.

Then the full local mixed-ray search on the free quintuple is completely exhausted.

### Exact final local mixed-ray closure theorem

The best local closure time on the free primitive quintuple is fully determined, up to the already-certified intervals of the carried packets, by the finite ledger
\[
\boxed{
\Bigl(
\tau_{\le 4,\min}^{\rm lo},\tau_{\le 4,\min}^{\rm hi},
\tau_{5,\min}^{\rm lo,int},\tau_{5,\min}^{\rm hi,int}
\Bigr),
}
\]
and there are no additional support-cardinality families beyond this stage.

Equivalently, the full local search on `\(\Delta_5^+\)` is no longer a continuum optimization problem. It is a finite certified ledger with one final splice.

So after Stage 218 there is no more generic local mixed-ray ranking work left to do.
The next honest move is no longer “rank larger support families.” It is to insert the actual PDE-derived branch data into this completed ledger.

---

## 7. Exact evaluation-budget theorem

The preferred lifted candidate compiler from Stage 217 contributes at most
\[
162
\]
interior five-coordinate stationary candidates **per** envelope, hence
\[
\boxed{324}
\]
across the `{lo,hi}` envelopes.

The already-finished support-`<=4` search budget from Stage 215 is
\[
\boxed{1140}.
\]
So the exact preferred total budget for the full support-`<=5` local search is
\[
\boxed{
1140+324=1464.
}
\]

If one instead falls back to the projected one-chart quintuple elimination bound from Stage 217, then the support-`5` interior contributes at most
\[
2\times 750 = 1500
\]
candidate evaluations, giving the fallback total
\[
\boxed{
1140+1500=2640.
}
\]

So even the fallback full local search is finite and explicit.

---

## 8. Best current summary after Stage 218

The local mixed-ray sieve is now finished.

- The whole support-`<=4` boundary is already closed by Stage 215.
- The only new support-`5` content was the unique interior five-coordinate optimizer, and Stage 217 reduced that to a finite algebraic candidate set.
- Stage 218 proves that the boundary of the unique five-simplex is exactly the already-solved support-`<=4` search set.
- Therefore the full local search over the positive free-quintuple simplex is reduced to one exact splice between the support-`<=4` global ledger and the unique support-`5` interior packet.
- The preferred exact total evaluation budget is `1464`, with fallback chart budget `2640`.

So the natural continuation after Stage 218 is **not** another support-cardinality theorem.
It is to begin inserting the actual PDE-derived Hessian-envelope and branch data into the completed search ledger — or, equivalently, to port this whole audited completion back into the compact PDE program master.
