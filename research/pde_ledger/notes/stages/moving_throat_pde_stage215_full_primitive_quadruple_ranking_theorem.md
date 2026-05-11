# Moving-Throat PDE — Stage 215: Full Primitive-Quadruple Ranking, Exact Boundary Splice, and the Up-to-Four-Coordinate Search Sieve

## Status

**Exact within the carried Stage 246 up-to-three-coordinate certified ledger and the Stage 248 finite interior four-coordinate candidate reduction, once the admissible interior quadruple packets are supplied on the chosen four-coordinate patches.**

This stage does **not** introduce a new constitutive law or a new optimizer.
It upgrades the Stage 248 interior four-coordinate candidate problem into a **global certified ranking ledger** for all primitive four-coordinate simplices, and it splices that ledger to the already-finished support-`<=3` search.

---

## Purpose

Stage 248 solved the genuinely new interior problem on a fixed primitive four-coordinate simplex: every admissible interior optimizer is contained in a finite algebraic candidate set, and the interior-versus-boundary comparison is again a certified interval problem.

That leaves one exact continuation point:

> how do we rank the **full set** of primitive four-coordinate simplices against one another, and how do we splice those interior winners together with the already-finished support-`<=3` ledger so the full search up through four active primitive coordinates becomes finite and certified?

This stage answers that.

The main outputs are:

1. the exact combinatorial ledger for the five primitive quadruples,
2. the exact **boundary-splice theorem** for every primitive quadruple,
3. the exact **local quadruple classification theorem** into interior-certified, boundary-certified, and unresolved simplices,
4. the exact **primitive-quadruple ranking theorem** by disjoint certified intervals,
5. the exact **global support-cardinality-4 improvement / non-improvement theorems** against the already-finished support-`<=3` ledger,
6. the exact **up-to-four-coordinate search sieve theorem**,
7. and the exact **finite evaluation budget theorem** for the whole support-`<=4` local search.

So Stage 249 is the global ranking complement of Stage 248: the four-coordinate interior candidate problem is no longer a deferred continuum search, and the full search up through four active primitive coordinates is now reduced to one finite certified ledger.

---

## 1. Carry-forward primitive triple and primitive quadruple ledgers

Let the free primitive index set be
\[
\boxed{
\mathfrak I_5:=\{\lambda,c,\gamma,U,W\}.
}
\]
The primitive triple set and primitive quadruple set are
\[
\boxed{
\mathfrak T_3
:=
\bigl\{\{i,j,k\}\subset\mathfrak I_5: i<j<k\bigr\},
\qquad
\#\mathfrak T_3=\binom53=10,
}
\]
\[
\boxed{
\mathfrak Q_4
:=
\bigl\{\{i,j,k,l\}\subset\mathfrak I_5: i<j<k<l\bigr\},
\qquad
\#\mathfrak Q_4=\binom54=5.
}
\]
Every primitive quadruple has exactly four codimension-one faces, and each face is one of the primitive triples. Dually, every primitive triple belongs to exactly two primitive quadruples, because after fixing three primitive directions there are exactly two remaining choices for the fourth.

### 1.1 Imported Stage 246 boundary packets

For each primitive triple
\[
T=(i,j,k)\in\mathfrak T_3,
\]
carry forward the exact closed-simplex certified interval
\[
\boxed{
\tau_{T,\min}^{\rm lo,\triangle}
\le
\tau_{T,*}^{\rm best,\triangle}
\le
\tau_{T,\min}^{\rm hi,\triangle}.
}
\]
These packets already include

- the full pairwise boundary splice, and
- the exact best certified closure time on the **closed** primitive triple simplex.

So by the time Stage 249 starts, every codimension-one face of a primitive quadruple is already solved.

### 1.2 Imported Stage 248 interior quadruple packets

For each primitive quadruple
\[
Q=(i,j,k,l)\in\mathfrak Q_4,
\]
carry forward the exact interior certified interval
\[
\boxed{
\tau_{Q,\min}^{\rm lo,int}
\le
\tau_{Q,*}^{\rm best,int}
\le
\tau_{Q,\min}^{\rm hi,int}.
}
\]
Here `\(\tau_{Q,*}^{\rm best,int}\)` is the unknown true best closure time over the **interior** of the positive four-coordinate simplex, while the entire codimension-one boundary is already the Stage 246 triple ledger.

---

## 2. Exact boundary-splice theorem for a primitive quadruple

Fix one primitive quadruple
\[
Q=(i,j,k,l)\in\mathfrak Q_4.
\]
Its four codimension-one faces are the primitive triple simplices
\[
(i,j,k),\qquad (i,j,l),\qquad (i,k,l),\qquad (j,k,l).
\]
Define the exact boundary minima
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
Let
\[
\tau_{Q,*}^{\rm best,\square}
\]
be the unknown true best closure time on the **closed** four-coordinate simplex, i.e. its interior plus all four triple faces.
Define the full quadruple certified interval
\[
\boxed{
\tau_{Q,\min}^{\rm lo,\square}
:=
\min\bigl(\beta_Q^{\rm lo},\tau_{Q,\min}^{\rm lo,int}\bigr),
}
\]
\[
\boxed{
\tau_{Q,\min}^{\rm hi,\square}
:=
\min\bigl(\beta_Q^{\rm hi},\tau_{Q,\min}^{\rm hi,int}\bigr).
}
\]

### Exact boundary-splice theorem

For every primitive quadruple `\(Q=(i,j,k,l)\)`, the best closure time on the closed simplex obeys
\[
\boxed{
\tau_{Q,\min}^{\rm lo,\square}
\le
\tau_{Q,*}^{\rm best,\square}
\le
\tau_{Q,\min}^{\rm hi,\square}.
}
\]

The proof is the same as in Stage 246: the true closed-simplex winner is the minimum of

- the true best boundary winner on the union of the four triple faces, and
- the true best interior winner.

Each of those pieces already has an exact certified interval, so the best over their union is bounded by the minima of the lower and upper certified bounds.

So after Stage 249, every primitive quadruple carries **one exact closed-simplex certified interval**.

---

## 3. Exact local quadruple classification theorem

The boundary splice immediately splits the primitive quadruples into three exact classes.

### 3.1 Interior-certified quadruples

If the interior upper bound already lies below the boundary lower bound,
\[
\boxed{
\tau_{Q,\min}^{\rm hi,int}<\beta_Q^{\rm lo},
}
\]
then every admissible interior winner beats every admissible boundary winner on that simplex.
So `\(Q\)` is **interior-certified**, and the closed-simplex winner is genuinely four-coordinate.

### 3.2 Boundary-certified quadruples

If the interior lower bound already lies above the boundary upper bound,
\[
\boxed{
\tau_{Q,\min}^{\rm lo,int}>\beta_Q^{\rm hi},
}
\]
then no admissible interior stationary candidate beats the best triple-face winner on that simplex.
So `\(Q\)` is **boundary-certified**, and that four-coordinate simplex contributes nothing beyond the already-ranked support-`<=3` ledger.

### 3.3 Ambiguous quadruples

Otherwise the certified interior and boundary intervals overlap, and `\(Q\)` remains **unresolved** at the certified level, even though its admissible interior candidate set is finite.

Define the three exact quadruple classes
\[
\boxed{
\mathfrak Q_4^{\rm int}
:=
\Bigl\{Q\in\mathfrak Q_4:
\tau_{Q,\min}^{\rm hi,int}<\beta_Q^{\rm lo}
\Bigr\},
}
\]
\[
\boxed{
\mathfrak Q_4^{\rm bdry}
:=
\Bigl\{Q\in\mathfrak Q_4:
\tau_{Q,\min}^{\rm lo,int}>\beta_Q^{\rm hi}
\Bigr\},
}
\]
\[
\boxed{
\mathfrak Q_4^{\rm amb}
:=
\mathfrak Q_4\setminus\bigl(\mathfrak Q_4^{\rm int}\cup\mathfrak Q_4^{\rm bdry}\bigr).
}
\]

So Stage 249 gives the exact local answer to the question

> is a given primitive quadruple truly a new support-cardinality-`4` simplex, or is it already exhausted by its triple-face boundary?

---

## 4. Exact primitive-quadruple ranking theorem

Take two primitive quadruples
\[
Q_1,Q_2\in\mathfrak Q_4.
\]
If their certified closed-simplex intervals are disjoint and ordered,
\[
\boxed{
\tau_{Q_1,\min}^{\rm hi,\square}
<
\tau_{Q_2,\min}^{\rm lo,\square},
}
\]
then every admissible winner on `\(Q_1\)` lies strictly below every admissible winner on `\(Q_2\)`.
So `\(Q_1\)` is certifiably better than `\(Q_2\)`.

### Exact unique certified primitive-quadruple winner theorem

If for some primitive quadruple `\(Q_\star\in\mathfrak Q_4\)` one has
\[
\boxed{
\tau_{Q_\star,\min}^{\rm hi,\square}
<
\min_{Q\in\mathfrak Q_4\setminus\{Q_\star\}}
\tau_{Q,\min}^{\rm lo,\square},
}
\]
then `\(Q_\star\)` is the **unique certified primitive-quadruple winner**.
Its exact closed-simplex best closure time is guaranteed to lie strictly below the best closure time on every other primitive quadruple simplex.

Define the global primitive-quadruple interval ledger
\[
\boxed{
\tau_{4,\min}^{\rm lo,\square}
:=
\min_{Q\in\mathfrak Q_4}
\tau_{Q,\min}^{\rm lo,\square},
\qquad
\tau_{4,\min}^{\rm hi,\square}
:=
\min_{Q\in\mathfrak Q_4}
\tau_{Q,\min}^{\rm hi,\square}.
}
\]
Then the best closure time over **all** primitive quadruple simplices obeys
\[
\boxed{
\tau_{4,\min}^{\rm lo,\square}
\le
\tau_{4,*}^{\rm best,\square}
\le
\tau_{4,\min}^{\rm hi,\square}.
}
\]

---

## 5. Exact global support-cardinality-4 improvement theorem against the support-`<=3` ledger

Carry forward the exact Stage 246 certified interval for the whole support-`<=3` search,
\[
\boxed{
\tau_{\le 3,\min}^{\rm lo}
\le
\tau_{\le 3,*}^{\rm best}
\le
\tau_{\le 3,\min}^{\rm hi}.
}
\]
Because every codimension-one face of a primitive quadruple is already a primitive triple, the only genuinely new support-cardinality-`4` content is the quadruple **interior**.

Define the global interior quadruple interval
\[
\boxed{
\tau_{4,\min}^{\rm lo,int}
:=
\min_{Q\in\mathfrak Q_4}
\tau_{Q,\min}^{\rm lo,int},
\qquad
\tau_{4,\min}^{\rm hi,int}
:=
\min_{Q\in\mathfrak Q_4}
\tau_{Q,\min}^{\rm hi,int}.
}
\]

### Exact global support-cardinality-4 improvement theorem

If
\[
\boxed{
\tau_{4,\min}^{\rm hi,int}<\tau_{\le 3,\min}^{\rm lo},
}
\]
then there exists a genuine four-coordinate interior mixed ray whose certified closure time lies strictly below every already-ranked support-`<=3` ray.

### Exact global support-cardinality-4 no-improvement theorem

If
\[
\boxed{
\tau_{4,\min}^{\rm lo,int}>\tau_{\le 3,\min}^{\rm hi},
}
\]
then no admissible four-coordinate interior stationary ray can beat the current support-`<=3` winner.

So after Stage 249, the full question

> do genuinely interior support-cardinality-`4` rays matter?

has become a single certified interval comparison.

---

## 6. Exact up-to-four-coordinate search sieve theorem

The support-`<=3` ledger already contains

- every primitive axis,
- every primitive pair,
- and every primitive triple,

while the quadruple interior ledger contains the genuinely new four-coordinate content.
Therefore the best closure time among **all** rays with support on at most four primitive coordinates obeys
\[
\boxed{
\tau_{\le 4,\min}^{\rm lo}
:=
\min\bigl(
\tau_{\le 3,\min}^{\rm lo},
\tau_{4,\min}^{\rm lo,int}
\bigr),
}
\]
\[
\boxed{
\tau_{\le 4,\min}^{\rm hi}
:=
\min\bigl(
\tau_{\le 3,\min}^{\rm hi},
\tau_{4,\min}^{\rm hi,int}
\bigr).
}
\]
Then the exact global support-`<=4` sieve theorem is
\[
\boxed{
\tau_{\le 4,\min}^{\rm lo}
\le
\tau_{\le 4,*}^{\rm best}
\le
\tau_{\le 4,\min}^{\rm hi}.
}
\]

So after Stage 249, the full local search up through four active primitive coordinates has been reduced to one finite certified interval comparison between

- the already-closed support-`<=3` ledger, and
- the finite interior quadruple ledger.

---

## 7. Exact finite evaluation budget theorem

Stage 248 already reduced every primitive quadruple interior to a finite algebraic candidate set with at most
\[
\boxed{54}
\]
admissible lifted stationary candidates per envelope.
There are exactly five primitive quadruples.
So the exact four-coordinate interior budget per envelope is
\[
\boxed{5\times 54 = 270}
\]
exact candidate evaluations.
Across the `\({\rm lo}/{\rm hi}\)` envelopes, the full support-cardinality-`4` interior budget is therefore
\[
\boxed{5\times 2\times 54 = 540}
\]
exact candidate evaluations.

Stage 246 already closed the support-`<=3` search with the exact finite budget
\[
\boxed{600}
\]
exact candidate evaluations.
Therefore the entire local search up through four active primitive coordinates has the exact finite budget
\[
\boxed{600+540=1140}
\]
exact candidate evaluations once the local slope data, envelope blocks, and admissible validity maps are supplied.

If the Stage 246 support-`<=3` ledger is already imported, then Stage 249 itself adds only the four-coordinate interior work:
\[
\boxed{540}
\]
new candidate evaluations plus exact interval comparisons.

This is the exact support-`<=4` cost theorem.

---

## 8. Minimal packet for the next stage

After Stage 249, the full search up through four active free coordinates is already finite and certified.
What remains beyond it is the unique five-coordinate positive simplex.

The natural next packet is therefore the five-coordinate simplex ledger
\[
\boxed{
\mathcal P_{\lambda c\gamma UW}^{(5)}
:=
\Bigl(
H_0,
(k_\lambda,k_c,k_\gamma,k_U,k_W),
H_{\lambda c\gamma UW,\rm lo},
H_{\lambda c\gamma UW,\rm hi},
T_{\lambda c\gamma UW}(\cdot)
\Bigr),
}
\]
with the exact boundary reduction now taken against the Stage 249 primitive-quadruple packets rather than against the Stage 246 triple ledger directly.

So the natural continuation is completely sharp:

> Stage 250 should build the exact five-coordinate positive simplex, reduce its full codimension-one boundary to the Stage 249 quadruple packets, and isolate the first canonical five-way interior screens before solving the full four-parameter interior optimizer.

---

## 9. Best current reading after Stage 249

Stage 246 turned every support-`<=3` search into a finite certified ledger.
Stage 248 did the same for the genuinely new interior of each four-coordinate simplex.
Stage 249 now closes the global support-cardinality-`4` ranking step.

1. There are exactly five primitive quadruples.
2. Every primitive quadruple now carries one exact closed-simplex certified interval obtained by splicing its finite interior packet to the imported triple-face boundary ledger.
3. The primitive quadruples split exactly into interior-certified, boundary-certified, and unresolved classes.
4. The full four-coordinate-versus-support-`<=3` question collapses to the single interval comparison
   \(
   \tau_{4,\min}^{\rm hi,int}\ ?\ \tau_{\le 3,\min}^{\rm lo}
   \).
5. The entire local search up through four active primitive coordinates now has a finite certified budget of at most `1140` exact candidate evaluations.

So the next theorem gate is no longer “is the four-coordinate search tractable?”
That problem is finished.
The next honest question is whether the unique five-coordinate simplex can produce a certified interior winner that beats the now-finite support-`<=4` ledger.
