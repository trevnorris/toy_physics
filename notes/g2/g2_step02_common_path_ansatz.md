# Step 2 — Minimal common quotient path and quartic matching

## Goal

Step 1 did two things:

1. it rebuilt the **exact current local anomaly law** as the baseline,
2. and it replaced the old staggered bookkeeping by the exact moving-throat quotient variables
   \[
   q(f)=\bigl(q_{\rm tr}(f),\,q_{\rm nt}(f),\,q_\eta(f)\bigr)
   =
   \bigl(\delta\ln \mathfrak C_{{\rm tr},*},\,
          \delta\ln \mathfrak C_{{\rm nt},*},\,
          \delta\ln \epsilon_\eta\bigr).
   \]

The next question is therefore no longer “what quartic coefficient do we want?”
It is:

> **How can a genuinely common charge–inertia quotient path modify the exact current local law without reopening the already-frozen lower orders?**

This step answers that question.

---

## Inputs carried forward

### From Step 1 / the anomaly write-up

The current local closure is
\[
g_{\rm loc}(f)
=
2\Bigl[Q_{\rm loc}(f)-\eta_1(f)f^2\Bigr],
\]
with the remaining discrepancy naturally sitting at
\[
O(f^4),
\]
which is exactly where the anomaly write-up says a **common charge–inertia layer** should first enter.

At the series level,
\[
\frac{g_{\rm loc}}{2}
=
1+f-\frac{47}{36}f^2
+
\left(c_{3,q}+c_{3,i}\right)f^3
+
a_{4,\rm staggered}f^4
+O(f^5),
\]
where
\[
c_{3,q}=\frac{4-\pi}{\pi^2\kappa},
\qquad
c_{3,i}=\frac{11}{6}\kappa.
\]

So the exact current local law already contains a definite staggered quartic coefficient.
The missing common layer has to be an **incremental** correction on top of that.

### From the moving-throat quotient bridge

The actual common deformation must be organized by the exact quotient variables
\[
q_{\rm tr},\qquad q_{\rm nt},\qquad q_\eta,
\]
because the five remaining microscopic directions are only similarity-orbit motion.

So at this stage the only honest new datum is the **quotient path**
\[
q(f).
\]

---

## Step 2A — Split the current local law into frozen lower orders plus transport residue

The first clean move is to separate the already-frozen lower-order ledger from the local transport residue.

### 1. Frozen sharp-boundary part

Define
\[
\frac{g_{\rm sharp}}{2}
=
1+f-\frac{47}{36}f^2.
\]

This contains the tree, Schwinger, and sharp self-loop terms that we do **not** want to disturb.

### 2. Charge-side local transport residue

Define
\[
T_q(f):=
Q_{\rm loc}(f)-\bigl(1+f-f^2\bigr).
\]
Using the Step-1 reconstruction,
\[
T_q(f)
=
\frac{4-\pi}{\pi^2\kappa}f^3
+
\frac{4(\pi-3)}{\pi^3\kappa}f^4
+
O(f^5).
\]

### 3. Inertia-side local transport residue

Write
\[
\eta_1(f)=\frac{11}{36}+\Delta\eta_1(f),
\]
and define the part that actually enters `g/2` by
\[
T_i(f):=
-\Delta\eta_1(f)f^2.
\]
Then
\[
T_i(f)
=
\frac{11}{6}\kappa f^3
-
\frac{55}{6}\kappa^2 f^4
+
O(f^5).
\]

### 4. Total local transport residue

So the current exact local law factorizes as
\[
\frac{g_{\rm loc}}{2}
=
\frac{g_{\rm sharp}}{2}
+
T_q(f)+T_i(f).
\]

The crucial structural point is:

\[
T_q(f)=O(f^3),
\qquad
T_i(f)=O(f^3),
\qquad
T_{\rm loc}(f):=T_q(f)+T_i(f)=O(f^3).
\]

This is the right object for the common PDE layer to act on.

---

## Step 2B — Why a naive whole-law dressing is wrong

Suppose one tries to introduce the common quotient path by multiplying the **whole** current law by a common factor
\[
e^{\Lambda(f)},
\qquad
\Lambda(f)=\Lambda_1 f+O(f^2).
\]

Then
\[
\bigl(e^{\Lambda(f)}-1\bigr)\frac{g_{\rm loc}(f)}{2}
=
\Lambda_1 f + O(f^2).
\]

So a linear-in-`f` quotient path would immediately reopen the already-frozen tree and one-loop sectors.
The same problem occurs if one multiplies the whole charge factor `Q_loc(f)` rather than only its transport-generated part.

That means:

- the quotient path **cannot** be allowed to dress the whole anomaly law,
- and the common layer must act only on the already-existing **transport residue**.

This is the main conceptual result of the step.

---

## Step 2C — The general transport-residue dressing ansatz

Now let the quotient path be analytic at small `f`:
\[
q_{\rm tr}(f)=s_{\rm tr}f+O(f^2),\qquad
q_{\rm nt}(f)=s_{\rm nt}f+O(f^2),\qquad
q_\eta(f)=s_\eta f+O(f^2).
\]

At first order, the most general linear functionals of the quotient path are
\[
\Lambda_Q(f)=
\bigl(a_{\rm tr}s_{\rm tr}+a_{\rm nt}s_{\rm nt}+a_\eta s_\eta\bigr)f+O(f^2)
=: \lambda_Q f+O(f^2),
\]
\[
\Lambda_I(f)=
\bigl(b_{\rm tr}s_{\rm tr}+b_{\rm nt}s_{\rm nt}+b_\eta s_\eta\bigr)f+O(f^2)
=: \lambda_I f+O(f^2).
\]

Use them to dress only the transport-generated pieces:
\[
\frac{g_{\rm common}}{2}
=
\frac{g_{\rm sharp}}{2}
+
e^{\Lambda_Q(f)}T_q(f)
+
e^{\Lambda_I(f)}T_i(f).
\]

The SymPy expansion gives the exact leading common correction
\[
\Delta\!\left(\frac g2\right)_{\rm common}
=
\Bigl(c_{3,q}\lambda_Q+c_{3,i}\lambda_I\Bigr)f^4
+
O(f^5).
\]

So at quartic order the common layer depends only on **one scalar combination**
of the charge-side and inertia-side tangent projections.

### Order theorem

Because the transport residue starts at `O(f^3)`, a quotient path that starts as
\[
q(f)=O(f^n)
\]
produces its first common correction at
\[
O(f^{3+n}).
\]

Therefore:

- if the common layer is supposed to close the observed quartic gap,
- then the quotient path must begin linearly,
  \[
  q_i(f)=s_i f+O(f^2).
  \]

This is an important sharpening.
A quotient path that starts only at `f^2` is already too small; it would first appear at `O(f^5)` in the transport-residue framework.

---

## Step 2D — The minimal single-common-scalar specialization

The most economical realization of a genuinely **common** layer is not to dress the charge and inertia residues separately, but to multiply the **entire local transport residue** by a single scalar built from the quotient path:
\[
\Lambda_{\rm common}(f)
=
\bigl(w_{\rm tr}s_{\rm tr}+w_{\rm nt}s_{\rm nt}+w_\eta s_\eta\bigr)f+O(f^2)
=: \Lambda_1 f+O(f^2).
\]

Then
\[
\frac{g_{\rm common}}{2}
=
\frac{g_{\rm sharp}}{2}
+
e^{\Lambda_{\rm common}(f)}\,T_{\rm loc}(f),
\]
and the quartic correction collapses to
\[
\Delta\!\left(\frac g2\right)_{\rm common}
=
c_{3,\rm total}\,\Lambda_1\,f^4
+
O(f^5),
\qquad
c_{3,\rm total}=c_{3,q}+c_{3,i}.
\]

So the entire quartic common layer is controlled, at leading order, by one number:
\[
\Lambda_1 = w\cdot s.
\]

That is the smallest common deformation compatible with both:

1. the exact quotient structure from the moving-throat paper, and
2. the requirement that no lower-order anomaly terms get disturbed.

---

## Step 2E — Numerical quartic matching

Using the frozen Step-1 benchmark
\[
\kappa=1.177746578880,
\qquad
f=0.001161409732093,
\qquad
\Delta g = g_e-g_{\rm loc}\approx 2.27204390584705\times 10^{-12},
\]
the residual quartic coefficient is
\[
a_{4,\rm resid}
=
\frac{\Delta g}{2f^4}
\approx 0.624374101073809.
\]

At the same `\kappa`,
\[
c_{3,q}\approx 0.0738485256041,
\qquad
c_{3,i}\approx 2.15920206128,
\qquad
c_{3,\rm total}\approx 2.23305058688410.
\]

### General quartic matching line

The two-functional transport ansatz therefore obeys
\[
c_{3,q}\lambda_Q+c_{3,i}\lambda_I=a_{4,\rm resid}.
\]

Numerically,
\[
\lambda_I
=
0.289168907473010
-
0.0342017669065772\,\lambda_Q.
\]

So the quartic gap fixes one line in the two-dimensional
\[
(\lambda_Q,\lambda_I)
\]
plane.

### Useful special points

If the common correction is carried only by the inertia-side dressing,
\[
\lambda_Q=0,
\qquad
\lambda_I\approx 0.289168907473010.
\]

If it is carried only by the charge-side dressing,
\[
\lambda_I=0,
\qquad
\lambda_Q\approx 8.45479440471249.
\]

This is already suggestive:
the quartic gap is much more naturally carried by an inertia-linked common deformation than by a pure charge-side one, because the inertia cubic coefficient is so much larger.

### Minimal single-common-scalar match

For the genuinely common specialization
\[
\lambda_Q=\lambda_I=\Lambda_1,
\]
one gets
\[
\boxed{
\Lambda_1
=
\frac{a_{4,\rm resid}}{c_{3,\rm total}}
\approx 0.279605891931464.
}
\]

So the exact Step-1 residual fixes the first quotient-path tangent projection required by the minimal common ansatz.

---

## What Step 2 establishes

### 1. The common layer must act on the transport residue, not on the whole anomaly law

That is the only way to let a linearly switched-on quotient path
\[
q_i(f)=s_i f+O(f^2)
\]
produce a quartic correction without contaminating the already-frozen lower orders.

### 2. The quartic problem is now reduced to one scalar tangent constraint

At the most general first nontrivial level,
\[
c_{3,q}\lambda_Q+c_{3,i}\lambda_I=a_{4,\rm resid}.
\]

In the minimal truly common specialization, this collapses to
\[
w\cdot s=\Lambda_1\approx 0.279605891931464.
\]

So the next derivation step does **not** need the whole nonlinear PDE.
It only needs the first quotient-path tangent on the actual moving-throat branch, projected onto the common scalar selected above.

### 3. The residual is naturally inertia-weighted

Because
\[
c_{3,i}\gg c_{3,q},
\]
the quartic closure strongly prefers a correction with substantial inertia-linked support unless the charge-side quotient projection is anomalously large.

That is a real physical hint, not just an algebraic curiosity.

---

## Immediate continuation point

The next clean step is now:

1. choose an explicit microscopic map from the quotient monomials
   \[
   \mathfrak C_{{\rm tr},*},\quad
   \mathfrak C_{{\rm nt},*},\quad
   \epsilon_\eta
   \]
   into the common scalar dressing functional `\Lambda_{\rm common}(f)`,
2. compute its first tangent
   \[
   \Lambda_1=w\cdot s,
   \]
   on the actual moving-throat branch,
3. and check whether it lands near
   \[
   0.279605891931464.
   \]

That is the first point where the new PDE paper can directly improve the anomaly calculation instead of only reorganizing the bookkeeping.
