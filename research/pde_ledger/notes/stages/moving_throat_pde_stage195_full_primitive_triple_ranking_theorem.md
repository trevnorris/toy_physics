# Moving-Throat PDE — Stage 195: Full Primitive-Triple Ranking, Exact Boundary Splice, and the Up-to-Three-Coordinate Search Sieve

## Status

**Exact within the carried Stage-192 pairwise certified optimizer and the Stage-194 finite interior triple-candidate reduction, once the local pairwise ledger and the admissible interior triple packets are supplied on the chosen free-quintuple patch.**

This stage does **not** introduce a new constitutive law or a new optimizer.
It upgrades the Stage-194 interior algebraic reduction into a **global certified ranking ledger** for all primitive three-coordinate simplices, and it splices that ledger to the already-finished pairwise cones.

---

## Purpose

Stage 194 solved the genuinely new algebraic problem on a fixed primitive triple: every interior optimizer on a three-coordinate mixed simplex is contained in a finite algebraic candidate set, and the interior-versus-boundary comparison is again a certified interval problem.

That leaves one exact continuation point:

> how do we rank the **full set** of primitive triples against one another, and how do we splice those triple-interior winners together with the already-optimized pairwise boundary ledger so the full search up through three-coordinate mixed rays becomes finite and certified?

This stage answers that.

The main outputs are:

1. the exact combinatorial ledger for the five primitive free directions,
2. the exact **boundary-splice theorem** for every primitive triple,
3. the exact **local triple classification theorem** into interior-certified, boundary-certified, and unresolved simplices,
4. the exact **primitive-triple ranking theorem** by disjoint certified intervals,
5. the exact **global three-coordinate improvement / non-improvement theorems** against the pairwise ledger,
6. and the exact **finite evaluation budget theorem** for the whole up-to-three-coordinate local sieve.

So Stage 195 is the global ranking complement of Stage 194: the interior candidate problem is no longer a continuum search, and the full search up through three active free coordinates is now reduced to a finite certified ledger.

---

## 1. Carry-forward primitive pair and primitive triple ledgers

Let the free primitive index set be
\[
\boxed{
\mathfrak I_5:=\{\lambda,c,\gamma,U,W\}.
}
\]
Then the primitive pair set and primitive triple set are
\[
\boxed{
\mathfrak P_2
:=
\bigl\{\{i,j\}\subset\mathfrak I_5: i<j\bigr\},
\qquad
\#\mathfrak P_2=\binom52=10,
}
\]
\[
\boxed{
\mathfrak T_3
:=
\bigl\{\{i,j,k\}\subset\mathfrak I_5: i<j<k\bigr\},
\qquad
\#\mathfrak T_3=\binom53=10.
}
\]
Every primitive pair belongs to exactly three primitive triples, because after fixing `\((i,j)\)` there are exactly three remaining primitive directions from which the third coordinate can be chosen.

### 1.1 Imported Stage-192 pairwise packets

For each primitive pair `\((i,j)\in\mathfrak P_2\)`, carry forward the exact optimized certified interval
\[
\boxed{
\tau_{ij,\min}^{\rm lo}
\le
\tau_{ij,*}^{\rm best}
\le
\tau_{ij,\min}^{\rm hi}.
}
\]
Because the Stage-192 admissible pairwise candidate set already contains the admissible ratio endpoints, this pairwise packet already subsumes the primitive one-axis rays on that edge.
So Stage 195 does **not** need a separate primitive-ray ledger.

### 1.2 Imported Stage-194 interior triple packets

For each primitive triple `\((i,j,k)\in\mathfrak T_3\)`, carry forward the exact interior certified interval
\[
\boxed{
\tau_{ijk,\min}^{\rm lo,int}
\le
\tau_{ijk,*}^{\rm best,int}
\le
\tau_{ijk,\min}^{\rm hi,int}.
}
\]
Here `\(\tau_{ijk,*}^{\rm best,int}\)` is the unknown true best closure time over the **interior** of the positive three-coordinate simplex, while the simplex boundaries themselves are already the Stage-192 pairwise cones.

---

## 2. Exact boundary splice theorem for a primitive triple

Fix one primitive triple
\[
T=(i,j,k)\in\mathfrak T_3.
\]
Its three simplex edges are the Stage-192 pairwise cones
\[
(i,j),\qquad (i,k),\qquad (j,k).
\]
Define the exact boundary minima
\[
\boxed{
\beta_T^{\rm lo}
:=
\min\bigl(
\tau_{ij,\min}^{\rm lo},
\tau_{ik,\min}^{\rm lo},
\tau_{jk,\min}^{\rm lo}
\bigr),
}
\]
\[
\boxed{
\beta_T^{\rm hi}
:=
\min\bigl(
\tau_{ij,\min}^{\rm hi},
\tau_{ik,\min}^{\rm hi},
\tau_{jk,\min}^{\rm hi}
\bigr).
}
\]
Let
\[
\tau_{T,*}^{\rm best,\triangle}
\]
be the unknown true best closure time on the **closed** three-coordinate simplex, i.e. interior plus all three edges.
Then define the full triple certified interval
\[
\boxed{
\tau_{T,\min}^{\rm lo,\triangle}
:=
\min\bigl(\beta_T^{\rm lo},\tau_{ijk,\min}^{\rm lo,int}\bigr),
}
\]
\[
\boxed{
\tau_{T,\min}^{\rm hi,\triangle}
:=
\min\bigl(\beta_T^{\rm hi},\tau_{ijk,\min}^{\rm hi,int}\bigr).
}
\]

### Exact boundary-splice theorem

For every primitive triple `\(T=(i,j,k)\)`, the best closure time on the closed simplex obeys
\[
\boxed{
\tau_{T,\min}^{\rm lo,\triangle}
\le
\tau_{T,*}^{\rm best,\triangle}
\le
\tau_{T,\min}^{\rm hi,\triangle}.
}
\]

The proof is immediate: the true simplex winner is the minimum of

- the true boundary winner on the union of the three edges, and
- the true interior winner.

Each of those pieces already has an exact certified interval, so the best over their union is bounded by the minima of the lower and upper certified bounds.

So after Stage 195, every primitive triple has been reduced to **one finite certified interval** on its full closed simplex.

---

## 3. Exact local triple classification theorem

The boundary splice immediately splits the primitive triples into three exact classes.

### 3.1 Interior-certified triples

If the interior upper bound already lies below the boundary lower bound,
\[
\boxed{
\tau_{ijk,\min}^{\rm hi,int}<\beta_T^{\rm lo},
}
\]
then every admissible interior winner beats every admissible boundary winner on that simplex.
So the triple is **interior-certified**, and the closed-simplex winner is genuinely three-coordinate.

### 3.2 Boundary-certified triples

If the interior lower bound already lies above the boundary upper bound,
\[
\boxed{
\tau_{ijk,\min}^{\rm lo,int}>\beta_T^{\rm hi},
}
\]
then no admissible interior stationary candidate beats the boundary winner.
So the triple is **boundary-certified**, and the three-coordinate simplex contributes nothing beyond its already-imported Stage-192 edges.

### 3.3 Ambiguous triples

Otherwise,
\[
\beta_T^{\rm lo}\le \tau_{ijk,\min}^{\rm hi,int}
\qquad\text{and}\qquad
\tau_{ijk,\min}^{\rm lo,int}\le \beta_T^{\rm hi},
\]
one has only an overlapping certified interval picture. The triple is then **unresolved** at the certified level, even though its candidate set is finite.

Define the three exact triple classes
\[
\boxed{
\mathfrak T_3^{\rm int}
:=
\Bigl\{T\in\mathfrak T_3:
\tau_{T,\min}^{\rm hi,int}<\beta_T^{\rm lo}
\Bigr\},
}
\]
\[
\boxed{
\mathfrak T_3^{\rm bdry}
:=
\Bigl\{T\in\mathfrak T_3:
\tau_{T,\min}^{\rm lo,int}>\beta_T^{\rm hi}
\Bigr\},
}
\]
\[
\boxed{
\mathfrak T_3^{\rm amb}
:=
\mathfrak T_3\setminus
\bigl(
\mathfrak T_3^{\rm int}\cup\mathfrak T_3^{\rm bdry}
\bigr).
}
\]

### Exact local triple classification theorem

Every primitive triple belongs to exactly one of these three classes. The interior-certified triples are the only simplices that are already proven to contain a genuine three-coordinate winner.

---

## 4. Exact primitive-triple ranking theorem

For each triple `\(T\in\mathfrak T_3\)`, keep its full certified interval
\[
\mathcal I_T^{\triangle}:=
\bigl[
\tau_{T,\min}^{\rm lo,\triangle},
\tau_{T,\min}^{\rm hi,\triangle}
\bigr].
\]
Now compare any two primitive triples `\(T_a,T_b\in\mathfrak T_3\)`.

### Exact interval-ordering relation

If
\[
\boxed{
\tau_{T_a,\min}^{\rm hi,\triangle}
<
\tau_{T_b,\min}^{\rm lo,\triangle},
}
\]
then every admissible simplex winner from `\(T_a\)` beats every admissible simplex winner from `\(T_b\)`.

So the certified intervals define an exact partial order on the primitive triples.

### Exact primitive-triple winner theorem

If there exists a triple `\(T_\star\in\mathfrak T_3\)` such that
\[
\boxed{
\tau_{T_\star,\min}^{\rm hi,\triangle}
<
\min_{T\in\mathfrak T_3\setminus\{T_\star\}}
\tau_{T,\min}^{\rm lo,\triangle},
}
\]
then `\(T_\star\)` is the **unique certified primitive-triple winner**.

Its exact closed-simplex best closure time is guaranteed to lie strictly below the best closure time on every other primitive triple simplex.

Define the global primitive-triple interval ledger
\[
\boxed{
\tau_{3,\min}^{\rm lo,\triangle}
:=
\min_{T\in\mathfrak T_3}
\tau_{T,\min}^{\rm lo,\triangle},
\qquad
\tau_{3,\min}^{\rm hi,\triangle}
:=
\min_{T\in\mathfrak T_3}
\tau_{T,\min}^{\rm hi,\triangle}.
}
\]
Then the best closure time over **all** primitive triple simplices obeys
\[
\boxed{
\tau_{3,\min}^{\rm lo,\triangle}
\le
\tau_{3,*}^{\rm best,\triangle}
\le
\tau_{3,\min}^{\rm hi,\triangle}.
}
\]

---

## 5. Exact global three-coordinate improvement theorem against the pairwise ledger

Define the global pairwise best interval
\[
\boxed{
\tau_{2,\min}^{\rm lo}
:=
\min_{(i,j)\in\mathfrak P_2}
\tau_{ij,\min}^{\rm lo},
\qquad
\tau_{2,\min}^{\rm hi}
:=
\min_{(i,j)\in\mathfrak P_2}
\tau_{ij,\min}^{\rm hi}.
}
\]
Because the primitive one-axis rays are already embedded in the Stage-192 pairwise edges, this interval is the full certified ledger for all searches with at most two active primitive coordinates.

Now define the global **interior** triple interval
\[
\boxed{
\tau_{3,\min}^{\rm lo,int}
:=
\min_{T\in\mathfrak T_3}
\tau_{T,\min}^{\rm lo,int},
\qquad
\tau_{3,\min}^{\rm hi,int}
:=
\min_{T\in\mathfrak T_3}
\tau_{T,\min}^{\rm hi,int}.
}
\]
This is the only genuinely new three-coordinate content, because triple boundaries are already pairwise.

### Exact global three-coordinate improvement theorem

If
\[
\boxed{
\tau_{3,\min}^{\rm hi,int}<\tau_{2,\min}^{\rm lo},
}
\]
then there exists a genuine three-coordinate interior mixed ray whose certified closure time lies strictly below every pairwise cone winner.

### Exact global three-coordinate no-improvement theorem

If
\[
\boxed{
\tau_{3,\min}^{\rm lo,int}>\tau_{2,\min}^{\rm hi},
}
\]
then no admissible three-coordinate interior stationary ray can beat the global pairwise winner.

So after Stage 195 the full question “do genuinely interior three-coordinate rays matter?” has become a single certified interval comparison.

---

## 6. Exact up-to-three-coordinate search sieve theorem

The pairwise ledger already contains all one-axis primitive rays, and the triple interior ledger contains the genuinely new three-coordinate content. Therefore the best closure time among **all** rays with support on at most three primitive coordinates obeys
\[
\boxed{
\tau_{\le 3,\min}^{\rm lo}
:=
\min\bigl(
\tau_{2,\min}^{\rm lo},
\tau_{3,\min}^{\rm lo,int}
\bigr),
}
\]
\[
\boxed{
\tau_{\le 3,\min}^{\rm hi}
:=
\min\bigl(
\tau_{2,\min}^{\rm hi},
\tau_{3,\min}^{\rm hi,int}
\bigr).
}
\]

### Exact up-to-three-coordinate search sieve theorem

Let `\(\tau_{\le3,*}^{\rm best}\)` be the unknown true best closure time among all primitive rays, all pairwise mixed cones, and all genuine interior three-coordinate mixed simplices on the supplied local patch. Then
\[
\boxed{
\tau_{\le 3,\min}^{\rm lo}
\le
\tau_{\le3,*}^{\rm best}
\le
\tau_{\le 3,\min}^{\rm hi}.
}
\]
So the local search up through three active primitive coordinates is now completely reduced to one finite certified interval.

---

## 7. Exact finite evaluation budget theorem

The combinatorics are now explicit.

### 7.1 Imported Stage-192 pairwise budget

There are
\[
\#\mathfrak P_2=10
\]
primitive pairs.
Each Stage-192 pair needs at most
\[
12
\]
exact candidate evaluations (six for the lower envelope and six for the upper envelope) on a compact admissible ratio interval.
So the full pairwise ledger costs at most
\[
\boxed{10\times 12 = 120}
\]
exact candidate evaluations.

### 7.2 New Stage-194 interior triple budget

There are
\[
\#\mathfrak T_3=10
\]
primitive triples.
Each Stage-194 interior simplex needs at most
\[
24
\]
admissible algebraic candidates per envelope, hence at most
\[
48
\]
interior candidate evaluations per triple.
So the full interior triple ledger costs at most
\[
\boxed{10\times 48 = 480}
\]
exact candidate evaluations.

### 7.3 Full up-to-three-coordinate budget

Therefore the entire local search up through three active primitive coordinates has the exact finite budget
\[
\boxed{
120+480=600
}
\]
exact candidate evaluations once the local slope data, envelope blocks, and admissible validity maps are supplied.

If the Stage-192 pairwise ledger is already imported, then Stage 195 itself adds only the triple-interior work:
\[
\boxed{480}
\]
new candidate evaluations plus exact interval comparisons.

This is the first exact global cost theorem for the primitive local sieve.

---

## 8. Minimal packet for the next stage

After Stage 195, the full search up through three active free coordinates is already finite and certified.
What remains beyond it is the first genuinely four-coordinate interior search.

The natural next packet is therefore the positive four-coordinate simplex ledger
\[
\boxed{
\mathcal P_{ijkl}^{(4)}
:=
\Bigl(
H_0,
(k_i,k_j,k_k,k_l),
H_{ijkl,\rm lo},
H_{ijkl,\rm hi},
T_{ijkl}(\cdot)
\Bigr),
}
\]
with the exact boundary reduction now taken against the Stage-195 triple ledger rather than against the pairwise ledger directly.

So the natural continuation is completely sharp:

> Stage 196 should build the exact four-coordinate positive simplex, reduce its full boundary to the Stage-195 triple packets, and isolate the first canonical interior four-way screens before solving the full three-parameter interior optimizer.

---

## 9. Best current reading after Stage 195

Stage 192 turned every pairwise cone into a finite certified candidate problem.
Stage 194 did the same for the genuinely new interior of each three-coordinate simplex.
Stage 195 now closes the global ranking step.

1. There are exactly ten primitive pairs and ten primitive triples.
2. Every primitive triple now carries one exact closed-simplex certified interval obtained by splicing its finite interior packet to the imported pairwise boundary ledger.
3. The primitive triples split exactly into interior-certified, boundary-certified, and unresolved classes.
4. The full three-coordinate-versus-pairwise question collapses to the single interval comparison
   \(
   \tau_{3,\min}^{\rm hi,int}\ ?\ \tau_{2,\min}^{\rm lo}
   \).
5. The entire local search up through three active primitive coordinates now has a finite certified budget of at most `600` exact candidate evaluations.

So the next theorem gate is no longer “is the three-coordinate search tractable?”
That problem is finished.
The next honest question is whether the first four-coordinate simplex can produce a certified interior winner that beats the now-finite up-to-three-coordinate ledger.
