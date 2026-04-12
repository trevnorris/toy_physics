# Step 1 — Exact moving-throat quotient bridge and staggered anomaly benchmark

## Goal

Start the rebuild from the **new moving-throat PDE language**, not from the old sequential collar story.
The purpose of this first step is:

1. to replace the old “charge update, then inertia update, then charge update again” bookkeeping by the exact moving-throat quotient coordinates, and
2. to rebuild the **current** staggered anomaly law exactly so that the later common charge–inertia correction is compared against the right baseline.

This is the correct beginning because the old write-up already says the remaining miss is the first genuinely **common** charge–inertia layer at `O(f^4)`, while the new PDE handoff compresses the microscopic coupled branch into three exact quotient coordinates.

---

## Inputs carried forward

### From the anomaly write-up

The current exact reduced closure has

- `f = alpha_fs / (2 pi)`,
- `g_loc ≈ 2.00231930435865`,
- a remaining miss of about `2.27e-12`,
- and the explicit warning that the unresolved piece is the first **common charge–inertia transport layer** at `O(f^4)`.

The anomaly text is also explicit that the previous successful layers were still partially **staggered**:

1. charge support moved into the collar,
2. inertia backreacted through the `P_22` self-loop and blur,
3. charge then acquired its own local collar mode.

So the next correction is expected to come from a genuinely coupled local law rather than one more one-sided update.

### From the moving-throat PDE handoff

The new exact microscopic state vector is

```math
x=(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})},\ K_W^{(\mathrm{eff})},\ \mu_W,\ T_U),
```

with grouped weak-axisymmetric log-drift vector

```math
\delta\mathbf x=
\begin{pmatrix}
\lambda_1\\ c_1\\ \gamma_1\\ \kappa_U\\ \kappa_\eta\\ \kappa_W\\ \mu_1\\ \tau_1
\end{pmatrix}.
```

The three exact quotient coordinates are the monomials

```math
\mathfrak C_{{\rm tr},*},\qquad
\mathfrak C_{{\rm nt},*},\qquad
\epsilon_\eta,
```

with exact invariant map

```math
\mathcal I(x)=\bigl(\mathfrak C_{{\rm tr},*}(x),\ \mathfrak C_{{\rm nt},*}(x),\ \epsilon_\eta(x)\bigr),
```

and exact log-drift law

```math
\begin{pmatrix}
\delta\ln \mathfrak C_{{\rm tr},*}\\
\delta\ln \mathfrak C_{{\rm nt},*}\\
\delta\ln \epsilon_\eta
\end{pmatrix}
=
M_*\,\delta\mathbf x.
```

The handoff further states that `rank(M_*) = 3`, so the true defect motion lives in a **3-dimensional quotient**, while the remaining 5 directions are pure similarity-orbit motion.

---

## Step 1A — Exact quotient coordinates and similarity-orbit decomposition

I define the quotient drift vector

```math
\mathbf q=
\begin{pmatrix}
q_{\rm tr}\\ q_{\rm nt}\\ q_\eta
\end{pmatrix}
:=
\begin{pmatrix}
\delta\ln \mathfrak C_{{\rm tr},*}\\
\delta\ln \mathfrak C_{{\rm nt},*}\\
\delta\ln \epsilon_\eta
\end{pmatrix}
=M_*\,\delta\mathbf x.
```

The SymPy audit verifies all of the following exactly.

### 1. `M_*` has rank 3 and nullity 5

So any microscopic grouped weak-axisymmetric drift splits into:

- **5 similarity directions** that preserve the exact monomials, and
- **3 quotient directions** that move the actual branch in
  `((mathfrak C_tr,*),(mathfrak C_nt,*),epsilon_eta)`.

### 2. There is an explicit right inverse `R`

Choosing the quotient lift in which the five free co-scalings

```math
(\lambda_1,c_1,\gamma_1,\kappa_U,\kappa_W)
```

are set to zero gives a concrete right inverse `R` with

```math
M_*R=I_3.
```

So the quotient motion is not abstract.  It can be injected back into the microscopic drift space exactly.

### 3. The dependent microscopic drifts solve exactly as

```math
\kappa_\eta
=
2c_1-\kappa_U-q_\eta,
```

```math
\tau_1
=
\kappa_U
-
\frac{1+\delta_{U,*}}{1+\chi_{0,*}}\,(\gamma_1+c_1-\kappa_U)
+
\frac{q_{\rm tr}}{1+\chi_{0,*}},
```

```math
\mu_1
=
2c_1-\kappa_U+2\kappa_W-2\lambda_1
-
E_*\,(2\gamma_1+2\lambda_1-\kappa_U-\kappa_W)
-
F_*\,\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\gamma_1+c_1-\kappa_U)
+
q_{\rm nt}-q_\eta+
\frac{F_*}{1+\chi_{0,*}}q_{\rm tr}.
```

Setting `q = 0` recovers the exact tangent **similarity orbit**.  So the new PDE picture makes the conceptual split precise:

- the five free co-scalings are not physical defect motion,
- the physical defect motion is exactly the 3-vector `q`.

### 4. The immediate consequence for the anomaly program

Any genuine **common** charge–inertia correction should be organized as a scalar built from

```math
q_{\rm tr},\qquad q_{\rm nt},\qquad q_\eta,
```

not from arbitrary sequential changes of charge and inertia support separately.

This is the key structural result of Step 1.

---

## Step 1B — Rebuild the current staggered anomaly law exactly

The common layer has to be added to the **current** local closure, not to a looser schematic formula.
So I also reconstructed the present staggered anomaly law directly from the Appendix F/G definitions.

### 1. Exact collar-mode integrals

Using the collar coordinate

```math
\tau(f)=1-\sqrt{1-f},
```

and the rescaled collar variable `s in [0,1]`, the charge-side mean subtraction and second moment are exact closed forms.

The exact mean subtraction is

```math
\bar c(\tau)
=
\frac{4\bigl(-2\tau+\pi(\tau-1)\bigr)}{\pi^2(\tau-2)}.
```

The exact second moment is

```math
\Xi(\tau)
=
\frac{2\tau\Bigl(-48\pi\tau^2+2\pi^2\tau^2+\pi^3\tau^2+96\tau^2
-3\pi^3\tau-4\pi^2\tau+48\pi\tau-8\pi^2+2\pi^3\Bigr)}{\pi^4(\tau-2)}.
```

This is useful because it means the current charge-side local mode can be benchmarked exactly, not only to cubic order.

### 2. Exact `Q_loc(f)` series through quartic order

The reconstructed charge second-moment ratio is

```math
Q_{\rm loc}(f)
=
1+f-f^2+
\frac{4-\pi}{\pi^2\kappa}f^3
+
\frac{4(\pi-3)}{\pi^3\kappa}f^4
+
O(f^5).
```

So the current **charge-side** quartic contribution is already fixed by the existing staggered closure.

### 3. Exact current staggered `g_loc(f)` series through quartic order

Combining the charge-side series with the exponential-blur inertia factor gives

```math
g_{\rm loc}(f)/2
=
1+f-
\frac{47}{36}f^2
+
\left(
\frac{11}{6}\kappa+
\frac{4-\pi}{\pi^2\kappa}
\right)f^3
+
a_{4,\rm staggered}f^4
+
O(f^5),
```

with

```math
a_{4,\rm staggered}
=
-
\frac{55}{6}\kappa^2
+
\frac{4(\pi-3)}{\pi^3\kappa}.
```

For the frozen benchmark value

```math
\kappa=1.177746578880,
```

this evaluates to

```math
a_{4,\rm staggered}\approx -12.6994546522869.
```

This is important:

> the current exact staggered law already carries a definite quartic coefficient.

So the unresolved common layer should **not** be identified with “the quartic coefficient of the full exact current local law.”
It is an **incremental** coupled correction on top of that law.

### 4. Exact physical benchmark with the current numbers

Using the frozen benchmark inputs, the reconstructed exact local closure gives

```math
g_{\rm loc}
=
2.002319304358647956\ldots
```

against the adopted target

```math
g_e = 2.00231930436092.
```

So the exact residual is

```math
\Delta g
=
g_e-g_{\rm loc}
=
2.27204390584705\times 10^{-12}.
```

That agrees with the `2.27e-12` scale reported in the write-up, but now we have the exact number used by the script.

---

## Step 1C — The quartic sign issue

If one compares the measured residual only against the **cubic truncation**

```math
g/2 = 1+f-c_2f^2+c_3f^3,
```

then the raw added quartic coefficient would be

```math
a_{4,\rm bench}
=
\frac{\Delta g}{2f^4}
\approx 0.624374101073809.
```

But Part VIII writes the next term as

```math
-c_4 f^4.
```

So in **that** sign convention,

```math
c_4 = -a_{4,\rm bench}
\approx -0.624374101073809.
```

That sign matters.
The positive number is the raw coefficient of `+f^4`; the Part VIII `c_4` is the negative of that.

---

## What Step 1 establishes

### 1. The new PDE base naturally replaces staggered bookkeeping by exact quotient motion

The moving-throat handoff gives a clean split:

- `5` similarity directions = zero-cost co-scalings,
- `3` quotient directions = true defect motion.

So the common charge–inertia layer should be written as a coupled evolution for the quotient vector

```math
\mathbf q=(q_{\rm tr},q_{\rm nt},q_\eta),
```

not as one more sequential update of the old charge and inertia ledgers.

### 2. The common layer must be benchmarked against the exact current local law

The exact present staggered closure already contains a real quartic series coefficient,
so the missing common term is **incremental**:

```math
\Delta g_{\rm common}(f)=O(f^4)
```

on top of the exact current `g_loc(f)`.

### 3. The Step-2 problem is now sharp

The next step should be:

1. choose a one-parameter coupled transport path
   ```math
   \mathbf q(f),
   ```
2. impose minimal conditions so that the new common correction begins at `O(f^4)`,
3. expand the lowest-order scalar correction built from
   `q_tr`, `q_nt`, and `q_eta`,
4. and match that correction against the exact residual **after** subtracting the current staggered baseline.

That is the cleanest beginning I know how to make from the new PDE paper.
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
# Step 3 — Exact coherent-defect normal form and the first PDE-fixed quartic constraint

## Goal

Step 2 showed that the first genuinely common correction can only enter through a **single scalar tangent**

a common scalar dressing functional `Lambda_common(f)` with first tangent
```math
Lambda_1.
```

But Step 2 still left that scalar abstract:
```math
Lambda_common = w\cdot q,
\qquad
q=(q_{\rm tr},q_{\rm nt},q_\eta).
```

The purpose of this step is to remove that arbitrariness.
Using the exact coherent-branch quotient theorem from the moving-throat notes, we can identify the first PDE-native scalar combinations explicitly.

---

## Inputs carried forward

### From Step 2

The quartic residual matching problem reduced to
```math
Lambda_1\approx 0.279605891931464,
```
with the common correction entering at `O(f^4)`.

### From the coherent moving-throat quotient theorem

On the constructive coherent branch, the exact quotient coordinates are the three monomial drifts
```math
q_{\rm tr}=\delta\ln \mathfrak C_{{\rm tr},*},
\qquad
q_{\rm nt}=\delta\ln \mathfrak C_{{\rm nt},*},
\qquad
q_\eta=\delta\ln \epsilon_\eta.
```

The same notes then give an exact triangular normal form for the three observable weak-axisymmetric defect channels:

```math
\Theta_1=-C_{{\rm tr},*}\,q_{\rm tr},
```

```math
\Xi_1=A_{{\rm tr},*}\,q_{\rm tr}+q_{\rm nt},
```

```math
\mathcal R_1+\Xi_1=-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,q_\eta,
```

with exact branch coefficients

```math
C_{{\rm tr},*}
=
\frac{\chi_{0,*}\,\delta_{U,*}}
{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})},
```

```math
A_{{\rm tr},*}
=
\frac{2\chi_{0,*}}
{(1+\chi_{0,*})(1+\delta_{U,*})}.
```

So the coherent branch already tells us something much stronger than Step 2 knew:
there are **three distinct scalar channels**, and only one of them is the direct grouped transport defect `Xi_1`.

---

## Step 3A — The arbitrary common scalar collapses to an exact coherent defect scalar

The direct observable coherent grouped defect is not a generic linear functional.
It is forced to be

```math
\boxed{
\Xi_1
=
A_{{\rm tr},*}\,q_{\rm tr}+q_{\rm nt}.
}
```

That is the first big simplification.
The Step-2 abstract weight vector is now fixed to

```math
w_{\Xi}=(A_{{\rm tr},*},\ 1,\ 0).
```

So if the missing quartic anomaly layer is identified with the **direct coherent transport defect**, then the common scalar is

```math
\Lambda_{\rm common}(f)=\Xi_1(f).
```

Its first tangent is therefore

```math
\boxed{
\Lambda_1
=
A_{{\rm tr},*}\,s_{\rm tr}+s_{\rm nt},
}
```

where
```math
q_i(f)=s_i f+O(f^2).
```

Using the Step-2 benchmark,

```math
\boxed{
A_{{\rm tr},*}\,s_{\rm tr}+s_{\rm nt}
\approx
0.279605891931464.
}
```

This is the first exact quartic matching condition tied directly to the new PDE notes rather than to an arbitrary ansatz.

---

## Step 3B — The dressing monomial is *not* part of the direct defect at this order

The third quotient coordinate appears only in the selected-branch residual:

```math
\boxed{
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,q_\eta.
}
```

So at leading common order there is a clean separation:

- `q_tr` = tracking-feedthrough channel,
- `q_nt` = nontracking transfer-shape channel,
- `q_eta` = dressing / selected-branch demand channel.

This matters because it says the first quartic closure of the **direct** anomaly defect can be formulated entirely in the
```math
(q_{\rm tr},q_{\rm nt})
```
plane unless one explicitly decides to dress the selected-branch residual as well.

So the Step-2 single-common-scalar problem is already substantially simplified by the PDE base:

```math
q_\eta
```
is not needed for the direct quartic closure.

---

## Step 3C — The stripped nontracking option

The same normal form also shows that `Xi_1` contains a universal tracking feed-through piece.
Because

```math
\Theta_1=-C_{{\rm tr},*}q_{\rm tr},
```

the exact inverse reconstruction gives

```math
q_{\rm nt}
=
\Xi_1+
\frac{A_{{\rm tr},*}}{C_{{\rm tr},*}}\Theta_1.
```

So if the old staggered baseline is interpreted as already absorbing the universal tracking feed-through, then the genuinely new common scalar is simply the nontracking monomial drift itself:

```math
\Lambda_{\rm common}^{\rm strip}(f)=q_{\rm nt}(f).
```

Its quartic tangent condition is then even simpler:

```math
\boxed{
s_{\rm nt}
\approx
0.279605891931464.
}
```

This is not yet a theorem that the stripped choice is the right anomaly scalar, but it is now a sharply defined **PDE-native alternative** to the raw observable choice.

---

## Step 3D — Canonical tangent paths

Once the direct observable condition

```math
A_{{\rm tr},*} s_{\rm tr}+s_{\rm nt}=\Lambda_1
```

is fixed, there is still a one-parameter family of solutions in the `(q_tr,q_nt)` plane.
A natural canonical choice is the **minimum-norm quotient tangent** in the full three-dimensional quotient space.
Because the direct condition does not constrain `s_eta`, the minimum-norm solution sets

```math
s_\eta=0.
```

Then the Euclidean minimum-norm direct-matching path is

```math
\boxed{
\begin{pmatrix}
s_{\rm tr}\\
s_{\rm nt}\\
s_\eta
\end{pmatrix}_{\rm min}
=
\frac{\Lambda_1}{1+A_{{\rm tr},*}^2}
\begin{pmatrix}
A_{{\rm tr},*}\\
1\\
0
\end{pmatrix}.
}
```

So the canonical direct observable path is the ray parallel to the exact defect vector itself.

Two special degenerate rays are also useful benchmarks:

### Pure tracking closure
```math
s_{\rm nt}=s_\eta=0,
\qquad
s_{\rm tr}=\frac{\Lambda_1}{A_{{\rm tr},*}}.
```

### Pure nontracking closure
```math
s_{\rm tr}=s_\eta=0,
\qquad
s_{\rm nt}=\Lambda_1.
```

The second one coincides with the stripped nontracking choice.

---

## What Step 3 establishes

### 1. The PDE notes fix the first common scalar much more strongly than Step 2 alone

The direct coherent grouped defect is exactly

```math
\Xi_1=A_{{\rm tr},*}q_{\rm tr}+q_{\rm nt},
```

so the Step-2 arbitrary weight vector `w` is replaced by a PDE-forced one.

### 2. The dressing channel decouples from the direct quartic closure

At leading common order,

```math
q_\eta
```

enters only the selected-branch demand residual, not the direct coherent grouped defect.

### 3. The quartic anomaly match now reduces to one exact coherent-branch tangent condition

For the raw observable defect,

```math
A_{{\rm tr},*} s_{\rm tr}+s_{\rm nt}
\approx 0.279605891931464.
```

If tracking feed-through is regarded as already frozen into the baseline, this reduces further to

```math
s_{\rm nt}\approx 0.279605891931464.
```

So the next derivation no longer needs a generic scalar ansatz.
It needs the actual coherent-branch tangent of the quotient monomials.

---

## Immediate continuation point

The next clean move is now:

1. determine whether the anomaly residual should be identified with the raw coherent defect `Xi_1` or the stripped nontracking scalar `q_nt`,
2. insert the actual moving-throat branch data for
   ```math
   \chi_{0,*},\ \delta_{U,*},\ \epsilon_{\eta,*},
   ```
3. compute the corresponding tangent combination of
   ```math
   (q_{\rm tr},q_{\rm nt},q_\eta),
   ```
4. and compare it directly to the Step-2 quartic target.

That is the first point where the new PDE paper stops reorganizing the anomaly problem and starts actually fixing the remaining coefficient.
# Step 4 — Exact transfer-shape / outgoing-prefactor bridge

## Goal

Step 3 fixed the coherent weak-axisymmetric defect to the exact branch-normal-form scalar
```math
\Xi_1=A_{{\rm tr},*}q_{\rm tr}+q_{\rm nt},
```
and separated the corrected nontracking coordinate
```math
q_{\rm nt}=\delta\ln \mathfrak N_*.
```

The next question is the obvious one:

> **what microscopic moving-throat quantity is `\Xi_1` actually measuring?**

This step answers that. It shows that `\Xi_1` is simultaneously

1. the logarithmic slope of the raw transfer shape `\mathcal T^2`,
2. the logarithmic slope of the effective outgoing prefactor `P_A`,
3. and therefore the first direct microscopic slope that can carry the missing quartic anomaly layer.

So this is the first point where the anomaly remainder is tied directly to actual branch data rather than to an abstract quotient tangent.

---

## Inputs carried forward

### From Step 3

On the coherent local branch,
```math
\Theta_1=-C_{{\rm tr},*}q_{\rm tr},
```
```math
\Xi_1=A_{{\rm tr},*}q_{\rm tr}+q_{\rm nt},
```
with
```math
A_{{\rm tr},*}=
\frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})},
\qquad
C_{{\rm tr},*}=
\frac{\chi_{0,*}\delta_{U,*}}
{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}.
```

The exact corrected nontracking branch composite was
```math
\mathfrak N_*:=\mathcal T^2 R_{\rm tr}^{B_*},
\qquad
B_*=
\frac{2(1+\chi_{0,*}+\delta_{U,*})}{\delta_{U,*}}.
```

### From Step 2 / Step 3

The quartic target carried forward from the already-frozen local anomaly law is
```math
\Lambda_1\approx 0.279605891931464.
```

At Step 3 there were still two logically distinct ways to use it:

- **raw coherent-defect closure**
  ```math
  \Xi_1=\Lambda_1,
  ```
- **corrected nontracking closure**
  ```math
  q_{\rm nt}=\Lambda_1.
  ```

Step 4 shows exactly what those two choices mean microscopically.

---

## Step 4A — The corrected branch composite removes the tracking feed-through exactly

Using
```math
\delta\ln R_{\rm tr}=-C_{{\rm tr},*}q_{\rm tr},
```
and
```math
\delta\ln\mathfrak N_*
=
\delta\ln\mathcal T^2 + B_*\,\delta\ln R_{\rm tr},
```
we get
```math
\delta\ln\mathfrak N_*
=
\delta\ln\mathcal T^2 - B_* C_{{\rm tr},*} q_{\rm tr}.
```
But the branch coefficients satisfy the exact identity
```math
A_{{\rm tr},*}=B_* C_{{\rm tr},*}.
```
So
```math
\delta\ln\mathcal T^2
=
A_{{\rm tr},*}q_{\rm tr}+q_{\rm nt}
=
\Xi_1.
```
Therefore

```math
\boxed{\delta\ln\mathcal T^2=\Xi_1,}
```

while the corrected composite gives

```math
\boxed{\delta\ln\mathfrak N_*=q_{\rm nt}.}
```

So the two Step-3 scalars are now interpreted very sharply:

- `\Xi_1` = **raw transfer-shape slope**,
- `q_nt` = **corrected nontracking transfer-shape slope**.

---

## Step 4B — `\Xi_1` is also the outgoing-prefactor slope

The later moving-throat grouped-`P2` notes show that the weak-axisymmetric grouped defect can be written as the logarithmic slope of the effective transfer shape,
```math
\Xi_1
=
\frac{\delta\ln \mathcal T_{{\rm eff},A}^2}{\epsilon\lambda_A},
```
and equally as the logarithmic slope of the outgoing prefactor,
```math
\Xi_1=\frac{P_1}{P_0},
```
if
```math
P_A=P_0+\epsilon\lambda_A P_1+O(\epsilon^2).
```

So the coherent defect is not just a bookkeeping scalar. It is literally the **first normalized slope of the outgoing quadrupole bridge**.

On a one-port branch,
```math
\mathcal T_A=\mathcal T_{A,0}\,e^{\epsilon\lambda_A\tau_{\rm eff}}+O(\epsilon^2),
```
so
```math
\boxed{\Xi_1=2\tau_{\rm eff}.}
```

This is the cleanest bridge yet between the conservative quotient story and the outgoing branch story.

---

## Step 4C — Exact microscopic one-port slope formulas

On the actual minimal one-port continuum branch,
```math
\mathcal T_A^2
=
\frac{Z_{W,A}(1+\rho_A)^2}
{\Omega_{W,A}^2(1-\epsilon_{W,A})^2}
=
\frac{27\pi^2Gc_s^5}{20a^5c^5}
\frac{1-\epsilon_{\eta,A}}{R_{{\rm target},A}}.
```

Taking the weak-axisymmetric logarithmic slope gives two exact ledgers for the same scalar `\Xi_1`.

### Port-side mixed-sector ledger

Define the dimensionless slopes
```math
\sigma_Z:=\delta\ln Z_W,
\qquad
\sigma_\Omega:=\delta\ln\Omega_W,
```
```math
\sigma_\rho:=\frac{\delta\rho}{1+\rho},
\qquad
\sigma_{\epsilon_W}:=\frac{\delta\epsilon_W}{1-\epsilon_W}.
```
Then
```math
\boxed{
\Xi_1
=
\sigma_Z + 2\sigma_\rho - 2\sigma_\Omega + 2\sigma_{\epsilon_W}.
}
```

### Geometric / selected-branch ledger

Define
```math
\sigma_{c_s}:=\delta\ln c_s,
\qquad
\sigma_a:=\delta\ln a,
\qquad
\sigma_R:=\delta\ln R_{\rm target},
```
```math
\sigma_{\epsilon_\eta}:=\frac{\delta\epsilon_\eta}{1-\epsilon_\eta}.
```
Then
```math
\boxed{
\Xi_1
=
5\sigma_{c_s} - 5\sigma_a - \sigma_{\epsilon_\eta} - \sigma_R.
}
```

So the same coherent defect can be read in two complementary microscopic ways:

- as a mixed-port overlap / blocking / frequency slope,
- or as a sound-speed / throat-size / selected-target / dressing slope.

This is exactly the kind of simplification the old staggered derivation did not have.

---

## Step 4D — Quartic anomaly target in direct microscopic variables

Now combine the outgoing-prefactor identity
```math
\Xi_1=\frac{P_1}{P_0}
```
with the Step-3 quartic target.

### 1. Direct coherent-defect closure

If the missing quartic anomaly layer is identified with the raw coherent defect, then
```math
\boxed{\frac{P_1}{P_0}=\Lambda_1.}
```
So the required prefactor slope is directly fixed:
```math
\boxed{\frac{P_1}{P_0}\approx 0.279605891931464.}
```

### 2. Corrected nontracking closure

If the physically correct common scalar is instead the corrected nontracking branch invariant,
```math
q_{\rm nt}=\Lambda_1,
```
then because
```math
q_{\rm nt}=\Xi_1-A_{{\rm tr},*}q_{\rm tr},
```
we get
```math
\boxed{
\frac{P_1}{P_0} - A_{{\rm tr},*} s_{\rm tr} = \Lambda_1.
}
```
So the raw prefactor slope and the corrected nontracking slope differ only by the universal tracking feed-through.

### 3. Tracking-rigid branch

If the coherent branch is tracking-rigid at this order,
```math
s_{\rm tr}=0,
```
then the two closures coincide:
```math
\boxed{
\frac{P_1}{P_0}
=
\Xi_1
=
q_{{\rm nt},1}
=
\Lambda_1.
}
```

That is the cleanest first microscopic closure of the missing quartic anomaly layer.

---

## Step 4E — Final microscopic matching equations

Using the port-side and geometric slope ledgers, the corrected quartic matching law becomes

```math
\boxed{
\sigma_Z + 2\sigma_\rho - 2\sigma_\Omega + 2\sigma_{\epsilon_W}
=
\Lambda_1 + A_{{\rm tr},*} s_{\rm tr},
}
```

and equivalently

```math
\boxed{
5\sigma_{c_s} - 5\sigma_a - \sigma_{\epsilon_\eta} - \sigma_R
=
\Lambda_1 + A_{{\rm tr},*} s_{\rm tr}.
}
```

On the tracking-rigid branch these reduce to

```math
\boxed{
\sigma_Z + 2\sigma_\rho - 2\sigma_\Omega + 2\sigma_{\epsilon_W}
=
5\sigma_{c_s} - 5\sigma_a - \sigma_{\epsilon_\eta} - \sigma_R
=
0.279605891931464.
}
```

This is the strongest reduced g-2 statement reached so far in the new PDE language:

> the missing quartic anomaly layer is the required weak-axisymmetric logarithmic slope of the effective outgoing transfer shape, equivalently of the outgoing prefactor, and it can now be written directly in microscopic moving-throat variables.

---

## Why this step matters

The old anomaly derivation stalled because the final correction still looked like “one more transport term.”
Step 4 removes that vagueness.

The remaining problem is no longer:

> invent the quartic coefficient.

It is now:

> compute the actual branch slopes of
> `Z_W`, `\rho`, `\Omega_W`, `\epsilon_W`, `c_s`, `a`, `\epsilon_\eta`, and `R_{\rm target}`
> on the moving-throat branch, and see whether their exact combination lands on
> `0.279605891931464`.

That is the first point where the new PDE framework genuinely improves the anomaly calculation rather than merely rephrasing it.
# Step 5 — Isotropic one-port loading mismatch and minimal drift families

## Goal

Step 4 rewrote the missing quartic anomaly layer as the weak-axisymmetric outgoing-prefactor slope
```math
\Xi_1=\frac{P_1}{P_0}.
```
The next question is then:

> **on the simplest isotropic one-port branch, what microscopic quantity is `P_1/P_0` actually measuring, and what is the smallest set of reduced drifts that can realize the required quartic target?**

This step answers that.

The main result is that on the weak-axisymmetric grouped branch the anomaly scalar is just the **static loading mismatch**
```math
\Xi_{\rm load}
:=
\frac{P_1}{P_0}
=
\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0},
```
and on the canonical even-preserving branch it is the **only** remaining first-order defect.

So the quartic anomaly problem is no longer “solve all grouped data at once.”
At this level it becomes:

> find a weak static outgoing-transfer strengthening `N_{01}/N_0` and/or a weak static conservative softening `D_{01}/D_0` whose difference equals the required anomaly scalar.

---

## Inputs carried forward

### From the anomaly write-up

The current exact local anomaly closure already reaches
```math
g_{\rm loc}=2.00231930435865,
```
with the remaining miss naturally sitting at
```math
O(f^4),
```
which is exactly where the missing common charge–inertia layer should first enter.
The carried quartic target from Steps 1–4 is
```math
\Lambda_1\approx 0.279605891931464.
```

### From Step 4

The missing quartic layer was rewritten as the weak-axisymmetric outgoing-prefactor slope
```math
\Xi_1=\frac{P_1}{P_0},
```
with the equivalent microscopic slope ledgers
```math
\sigma_Z+2\sigma_\rho-2\sigma_\Omega+2\sigma_{\epsilon_W} = \Xi_1,
```
```math
5\sigma_{c_s}-5\sigma_a-\sigma_{\epsilon_\eta}-\sigma_R = \Xi_1,
```
on the tracking-rigid branch.

### From the later moving-throat grouped notes

On the weak-axisymmetric grouped branch,
```math
P_0^{(A)}=\frac{N_{A,0}}{D_{A,0}},
```
and the first physical prefactor slope is defined by
```math
P_0^{(A)}=P_0+\epsilon\lambda_A P_1+O(\epsilon^2).
```
The same notes give the exact transport identity
```math
\frac{P_1}{P_0} = \frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
```
That is the key simplification used here.

---

## Step 5A — Exact weak-axisymmetric physical slopes

Start from the first-order grouped expansions
```math
D_{A,0}=D_0+\epsilon\lambda_A D_{01},
```
```math
D_{A,2}=D_2+\epsilon\lambda_A D_{21},
```
```math
D_{A,4}=D_4+\epsilon\lambda_A D_{41},
```
```math
N_{A,0}=N_0+\epsilon\lambda_A N_{01}.
```
Then the grouped conservative response and outgoing prefactor are
```math
u_2^{(A)}=-\frac{D_{A,2}}{D_{A,0}},
```
```math
u_4^{(A)}=\frac{D_{A,2}^2-D_{A,0}D_{A,4}}{D_{A,0}^2},
```
```math
P_0^{(A)}=\frac{N_{A,0}}{D_{A,0}}.
```
Expanding to first order gives the exact physical slopes
```math
u_2^{(1)}=-\frac{D_0D_{21}-D_{01}D_2}{D_0^2},
```
```math
u_4^{(1)}=
\frac{D_0(-D_0D_{41}-D_{01}D_4+2D_2D_{21})+2D_{01}(D_0D_4-D_2^2)}{D_0^3},
```
```math
P_1=\frac{D_0N_{01}-N_0D_{01}}{D_0^2},
```
and therefore
```math
\boxed{
\frac{P_1}{P_0}
=
\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
}
```
So the outgoing-prefactor slope is exactly the difference between

- the weak static **outgoing-transfer** slope `N_{01}/N_0`, and
- the weak static **conservative operator** slope `D_{01}/D_0`.

---

## Step 5B — Canonical compensated/even-preserving branch

On the canonical compensated branch,
```math
u_2=\frac19,
\qquad
u_4=\frac{4}{81},
```
so the isotropic conservative moments are
```math
D_2=-\frac{D_0}{9},
\qquad
D_4=-\frac{D_0}{27}.
```
Requiring the grouped conservative response to remain fixed to first order,
```math
u_2^{(1)}=0,
\qquad
u_4^{(1)}=0,
```
forces
```math
\boxed{D_{21}=-\frac{D_{01}}{9},}
```
```math
\boxed{D_{41}=-\frac{D_{01}}{27}.}
```
So on the even-preserving branch the conservative grouped response is fully transported by **one** static slope `D_{01}`.
The even coefficients themselves do not carry an independent first-order defect.

At that point the only remaining first-order grouped outlet scalar is
```math
\boxed{
\Xi_{\rm load}
:=
\frac{P_1}{P_0}
=
\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
}
```
This is the first major simplification of the new PDE-based anomaly rebuild.

---

## Step 5C — Quartic anomaly target as a normalized loading equation

Match the anomaly scalar to the carried quartic target:
```math
\Xi_{\rm load}=\Lambda_1.
```
Define the normalized static drift variables
```math
n:=\frac{N_{01}}{N_0},
\qquad
d:=\frac{D_{01}}{D_0}.
```
Then the anomaly matching law is simply
```math
\boxed{n-d=\Lambda_1.}
```
That is the smallest honest one-port isotropic continuation of the Step-4 quartic constraint.

Three benchmark realizations are immediate.

### 1. Pure outgoing-transfer realization
```math
n=\Lambda_1,
\qquad
d=0.
```
All of the missing quartic loading is placed in the outgoing-transfer strengthening.

### 2. Pure conservative-softening realization
```math
n=0,
\qquad
d=-\Lambda_1.
```
All of the missing quartic loading is placed in the static conservative softening.

### 3. Balanced minimum-norm realization
Minimizing
```math
n^2+d^2
```
subject to
```math
n-d=\Lambda_1
```
gives
```math
\boxed{n=\frac{\Lambda_1}{2},
\qquad
d=-\frac{\Lambda_1}{2}.}
```
So the canonical balanced branch splits the needed loading equally between

- transfer strengthening, and
- conservative softening.

Numerically,
```math
n\approx 0.139802945965732,
\qquad
d\approx -0.139802945965732.
```

---

## Step 5D — Minimal microscopic conservative split

The static conservative operator on the isotropic branch is
```math
D_0=K-B_0-Z_0,
```
so its first normalized weak slope can be written as
```math
d = k-b-z,
```
with
```math
k:=\frac{K_{01}}{D_0},
\qquad
b:=\frac{B_{01}}{D_0},
\qquad
z:=\frac{Z_{01}}{D_0}.
```
Here:
- `k` = wall / geometry stiffness drift,
- `b` = BdG support dressing drift,
- `z` = conservative Maxwell/mixed dressing drift.

Minimizing
```math
k^2+b^2+z^2
```
subject to
```math
k-b-z=d
```
gives the unique minimum-norm conservative split
```math
\boxed{
(k,b,z)=\left(\frac d3,
-\frac d3,
-\frac d3\right).
}
```
On the balanced branch `d=-\Lambda_1/2`, this becomes
```math
\boxed{
(k,b,z)=\left(-\frac{\Lambda_1}{6},
\frac{\Lambda_1}{6},
\frac{\Lambda_1}{6}\right).
}
```
Numerically,
```math
k\approx -0.0466009819885773,
\qquad
b\approx 0.0466009819885773,
\qquad
z\approx 0.0466009819885773.
```
So the most economical conservative-softening realization lowers the bare wall stiffness and raises the support and Maxwell/mixed dressings by equal reduced amounts.

---

## Step 5E — Minimal one-port outgoing-transfer deformation

For a single port-active outgoing block, the later grouped notes reduce the normalized transfer slope to a portwise logarithmic deformation of the form
```math
n = 2\pi_1 - 2\delta_1,
```
where
- `\pi_1` is the port-amplitude logarithmic slope,
- `\delta_1` is the port-denominator logarithmic slope.

Minimizing
```math
\pi_1^2+\delta_1^2
```
subject to
```math
2\pi_1-2\delta_1=n
```
gives
```math
\boxed{
\pi_1=\frac n4,
\qquad
\delta_1=-\frac n4.
}
```
So the most economical one-port transfer deformation splits the outgoing slope symmetrically between

- amplitude strengthening, and
- denominator weakening.

On the balanced anomaly branch `n=\Lambda_1/2`, that becomes
```math
\boxed{
\pi_1=\frac{\Lambda_1}{8},
\qquad
\delta_1=-\frac{\Lambda_1}{8}.
}
```
Numerically,
```math
\pi_1\approx 0.0349507364914330,
\qquad
\delta_1\approx -0.0349507364914330.
```

---

## Step 5F — Final reduced verdict

The old derivation asked for a new common charge–inertia transport layer without a sharp microscopic split.
The new PDE language sharpens that dramatically.

On the simplest isotropic one-port branch,
```math
\boxed{
\Xi_1
=
\frac{P_1}{P_0}
=
\Xi_{\rm load}
=
\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
}
```
And on the canonical even-preserving branch, that is the **only** remaining first-order defect.

So the quartic anomaly problem is no longer

> derive every grouped coefficient.

It is now

> derive one weak static loading mismatch between outgoing-transfer strengthening and conservative static softening.

That is a real simplification, and it gives a much cleaner target for the next moving-throat calculation.
# Step 6 — Static self-similarity and the outgoing-load theorem

## Goal

Step 5 showed that on the canonical even-preserving isotropic one-port branch the missing quartic anomaly layer is just one scalar,
```math
\Xi_{\rm load}
:=
\frac{P_1}{P_0}
=
\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
```
That was already a major simplification, but it still left the two slopes
```math
D_{01}/D_0,
\qquad
N_{01}/N_0
```
looking like unrelated microscopic inputs.

The next honest question is therefore:

> **is the remaining grouped loading defect a generic static mismatch of the whole moving-throat bundle, or does it collapse to a much narrower wall-referenced theorem?**

This step answers that.

The main result is that the entire remaining linear grouped defect is exactly a **weighted failure of static self-similarity relative to the wall baseline**. After the natural wall-normalized factorization, the conservative bundles reduce to pure shape variables, while the outgoing bundle reduces to a wall-loading law. On conservative-shape-preserving branches, the anomaly target becomes a pure **outgoing-load theorem**.

---

## Inputs carried forward

### From Step 5

On the canonical even-preserving branch,
```math
u_2^{(1)}=0,
\qquad
u_4^{(1)}=0
```
forces
```math
D_{21}=-\frac{D_{01}}{9},
\qquad
D_{41}=-\frac{D_{01}}{27},
```
so the grouped conservative response is transported by one static slope `D_{01}`.
The only remaining first-order outlet scalar is then
```math
\Xi_{\rm load}
=
\frac{P_1}{P_0}
=
\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
```
Matching the anomaly requires
```math
\Xi_{\rm load}=\Lambda_1,
\qquad
\Lambda_1\approx 0.279605891931464.
```

### From the later moving-throat grouped notes

The grouped weak-axisymmetric notes sharpen that scalar further: the static operator and outgoing-transfer slopes can themselves be rewritten in wall-referenced logarithmic form, and the whole remaining grouped loading defect can be expressed as a weighted failure of static self-similarity between

- the wall baseline,
- the BdG support bundle,
- the conservative Maxwell/mixed bundle,
- and the outgoing-transfer bundle.

That is the language adopted here.

---

## Step 6A — Exact wall-referenced decomposition of the loading scalar

Start from the static conservative operator
```math
D_0=K-B_0-Z_0,
\qquad
D_{01}=K_1-B_0^{(1)}-Z_0^{(1)}.
```
Define the logarithmic slopes
```math
\delta_K:=\frac{K_1}{K},
\qquad
\delta_B:=\frac{B_0^{(1)}}{B_0},
\qquad
\delta_Z:=\frac{Z_0^{(1)}}{Z_0},
\qquad
\delta_N:=\frac{N_1}{N_0},
```
and the static weights
```math
\omega_K:=\frac{K}{D_0},
\qquad
\omega_B:=\frac{B_0}{D_0},
\qquad
\omega_Z:=\frac{Z_0}{D_0}.
```
Then
```math
\delta_D:=\frac{D_{01}}{D_0}
=
\omega_K\delta_K-\omega_B\delta_B-\omega_Z\delta_Z,
```
with the exact identity
```math
\omega_K-\omega_B-\omega_Z=1.
```
Therefore the Step-5 scalar becomes
```math
\Xi_{\rm load}
=
\delta_N-\delta_D
=
(\delta_N-\delta_K)
+
\omega_B(\delta_B-\delta_K)
+
\omega_Z(\delta_Z-\delta_K).
```
So the natural reference slope is not arbitrary. It is the wall-baseline slope `\delta_K`.

This is the first sharpened result of Step 6:

> the remaining quartic anomaly layer is the weighted failure of the three support/transfer sectors to co-load with the wall baseline.

---

## Step 6B — Microscopic weighted-log forms of the three sectors

The next step is to rewrite the three sector slopes as explicit weighted logarithmic drifts.

### BdG support bundle

For each support mode `\alpha`,
```math
B_{0,\alpha}=\frac{c_\alpha^2}{\varpi_\alpha^2},
\qquad
B_0=\sum_\alpha B_{0,\alpha},
```
so
```math
\delta_B
=
\sum_\alpha \rho_\alpha^{(B)}
\,2\,\delta\ln\!\left(\frac{c_\alpha}{\varpi_\alpha}\right),
\qquad
\rho_\alpha^{(B)}:=\frac{B_{0,\alpha}}{B_0}.
```

### Conservative Maxwell/mixed static bundle

For each port `r`,
```math
Z_0^{(r)}=\frac{Q_r}{\Delta_r},
\qquad
Z_0=\sum_r Z_0^{(r)},
```
so
```math
\delta_Z
=
\sum_r \rho_r^{(Z)}
\,\delta\ln\!\left(\frac{Q_r}{\Delta_r}\right),
\qquad
\rho_r^{(Z)}:=\frac{Z_0^{(r)}}{Z_0}.
```

### Outgoing-transfer bundle

For each outgoing port,
```math
N_0^{(r)}=\frac{P_r^2}{\Delta_r^2},
\qquad
N_0=\sum_r N_0^{(r)},
```
so
```math
\delta_N
=
\sum_r \rho_r^{(N)}
\,2\,\delta\ln\!\left(\frac{P_r}{\Delta_r}\right),
\qquad
\rho_r^{(N)}:=\frac{N_0^{(r)}}{N_0}.
```

So the Step-5 loading scalar is already controlled by three weighted microscopic log drifts.

---

## Step 6C — Exact self-similarity defect fields

Because `\delta_K` is the natural wall-baseline slope, define the wall-referenced defect fields
```math
\Sigma_\alpha^{(B)}
:=
2\,\delta\ln\!\left(\frac{c_\alpha}{\varpi_\alpha}\right)-\delta_K,
```
```math
\Sigma_r^{(Z)}
:=
\delta\ln\!\left(\frac{Q_r}{\Delta_r}\right)-\delta_K,
```
```math
\Sigma_r^{(N)}
:=
2\,\delta\ln\!\left(\frac{P_r}{\Delta_r}\right)-\delta_K.
```
Then the whole grouped defect becomes
```math
\boxed{
\Xi_{\rm load}
=
\sum_r \rho_r^{(N)}\Sigma_r^{(N)}
+
\omega_B\sum_\alpha \rho_\alpha^{(B)}\Sigma_\alpha^{(B)}
+
\omega_Z\sum_r \rho_r^{(Z)}\Sigma_r^{(Z)}.
}
```

This is the sharpest exact Step-6 formula.

It says the remaining linear grouped anomaly defect is not sourced by the whole microscopic bundle independently. It is sourced only by the weighted failure of

1. the BdG support bundle,
2. the conservative Maxwell/mixed bundle,
3. and the outgoing-transfer bundle,

to co-load with the wall baseline.

---

## Step 6D — Wall-normalized factorization into shape and load variables

The previous defect fields sharpen further after factoring by the wall baseline.

### BdG support shape

Define
```math
\chi_\alpha:=\frac{c_\alpha}{\sqrt K\,\varpi_\alpha}.
```
Then
```math
B_{0,\alpha}=K\,\chi_\alpha^2,
```
and the defect field becomes
```math
\Sigma_\alpha^{(B)}=\delta\ln\chi_\alpha^2.
```
So the BdG sector contributes only through the drift of the wall-normalized support shape `\chi_\alpha`.

### Conservative port shape and outgoing load

Define the wall-normalized port variables
```math
\Upsilon_r:=\frac{Q_r}{K\Delta_r},
\qquad
\Lambda_r:=\frac{P_r}{\Delta_r}.
```
Then
```math
Z_0^{(r)}=K\,\Upsilon_r,
\qquad
N_0^{(r)}=\Lambda_r^2,
```
so
```math
\Sigma_r^{(Z)}=\delta\ln\Upsilon_r,
```
```math
\Sigma_r^{(N)}=\delta\ln\!\left(\frac{\Lambda_r^2}{K}\right)=2\,\delta\ln\Lambda_r-\delta_K.
```

This is the second major Step-6 sharpening:

- the conservative bundles reduce to pure **shape** drifts,
- the outgoing bundle reduces to a **wall-loading law** for `\Lambda_r`.

Define the weighted defect measures
```math
\Theta_B:=\sum_\alpha \rho_\alpha^{(B)}\,\delta\ln\chi_\alpha^2,
```
```math
\Theta_Z:=\sum_r \rho_r^{(Z)}\,\delta\ln\Upsilon_r,
```
```math
\Theta_N:=\sum_r \rho_r^{(N)}\,\delta\ln\!\left(\frac{\Lambda_r^2}{K}\right).
```
Then the whole grouped defect collapses to the compact formula
```math
\boxed{
\Xi_{\rm load}=\Theta_N+\omega_B\Theta_B+\omega_Z\Theta_Z.
}
```

---

## Step 6E — Conservative-shape theorem and outgoing-load theorem

Now specialize to the branch on which the conservative shapes are preserved:
```math
\delta\ln\chi_\alpha^2=0
\quad\text{for all }\alpha,
```
```math
\delta\ln\Upsilon_r=0
\quad\text{for all }r.
```
Then
```math
\Theta_B=0,
\qquad
\Theta_Z=0,
```
and the remaining grouped defect collapses to
```math
\boxed{
\Xi_{\rm load}
=
\sum_r \rho_r^{(N)}\left(2\,\delta\ln\Lambda_r-\delta_K\right).
}
```

So on conservative-shape-preserving branches, the full remaining linear grouped `2.5`PN defect is carried **only** by the outgoing load factor `\Lambda_r`.

That is the main theorem of the step.

---

## Step 6F — Two immediate consequences

### 1. Naive common-self-similarity no-go

If one also freezes the outgoing load factor itself,
```math
\delta\ln\Lambda_r=0
\quad\text{for all }r,
```
then
```math
\Xi_{\rm load}=-\delta_K.
```
So a naive common self-similarity branch does **not** kill the remaining grouped defect unless the wall baseline itself does not load.

This is a real no-go result.

### 2. Exact outgoing-load law

The defect vanishes if and only if
```math
\boxed{
2\sum_r \rho_r^{(N)}\,\delta\ln\Lambda_r
=
\delta_K.
}
```
A stronger sufficient condition is the portwise law
```math
2\,\delta\ln\Lambda_r=\delta_K
\qquad\text{for every outgoing port }r.
```
So the old vague requirement “the outgoing bundle must be self-similar” is too weak.
The actual requirement is a precise wall-loading law.

---

## Step 6G — Direct quartic anomaly law

Now feed the Step-5 anomaly target back in:
```math
\Xi_{\rm load}=\Lambda_1.
```
On a conservative-shape-preserving branch this becomes the exact quartic theorem gate
```math
\boxed{
2\sum_r \rho_r^{(N)}\,\delta\ln\Lambda_r
-
\delta_K
=
\Lambda_1.
}
```
So the anomaly does **not** ask for an arbitrary new static bundle deformation.
It asks for a specific outgoing-load mismatch relative to the wall baseline.

If the outgoing ports all share a common weak logarithmic slope,
```math
\delta\ln\Lambda_r=\ell,
```
then
```math
\Xi_{\rm load}=2\ell-\delta_K,
```
so
```math
\ell_{\rm kill}=\frac{\delta_K}{2}
```
would kill the grouped defect,
and anomaly matching instead requires
```math
\boxed{
\ell_{\rm anom}
=
\frac{\delta_K+\Lambda_1}{2}.
}
```
Therefore the anomaly target requires an **extra outgoing-load slippage above the defect-killing wall-tracking law** of
```math
\boxed{
\ell_{\rm anom}-\ell_{\rm kill}
=
\frac{\Lambda_1}{2}.
}
```
Numerically,
```math
\frac{\Lambda_1}{2}
\approx 0.139802945965732.
```

This is the cleanest continuation of Step 5 reached so far.

---

## Step 6H — Final reduced verdict

Step 5 showed that the quartic anomaly correction is one scalar loading mismatch.
Step 6 sharpens that again:

1. `\Xi_{\rm load}` is exactly a wall-referenced self-similarity defect.
2. After wall-normalization, the conservative bundles become pure shape variables `\chi_\alpha` and `\Upsilon_r`.
3. On conservative-shape-preserving branches, the whole remaining linear grouped defect is carried only by the outgoing load factor `\Lambda_r=P_r/\Delta_r`.
4. So the quartic anomaly target is not “derive every grouped coefficient.”
5. It is:
   > determine whether the actual moving-throat branch produces the required outgoing-load slippage relative to the wall baseline.

That is a much smaller and much cleaner theorem gate than we had before.
# Step 7 — Direct outgoing-port co-loading theorem

## Goal

Step 6 reduced the remaining quartic anomaly closure to the outgoing-load law on the conservative-shape-preserving branch:
```math
\Xi_1
:=
\frac{P_1}{P_0}
=
2\sum_r \rho_r^{(N)}\,\delta\ln\Lambda_r-\kappa_1,
\qquad
\Lambda_r:=\frac{P_r}{\Delta_r}.
```
That was already a strong simplification, but it still phrased the theorem gate in terms of the abstract load factor `\Lambda_r`.

The next honest question is therefore:

> **can the remaining defect be rewritten directly in terms of the weak-axisymmetric slopes of the actual outgoing-port data**
> ```math
> P_r,
> \qquad
> \Delta_r,
> \qquad
> N_{0}^{(r)}=\frac{P_r^2}{\Delta_r^2},
> ```
> **so that the g-2 closure becomes a direct port co-loading law?**

This step answers that.

---

## Important notation warning

There are now **two different `P` symbols** in play.

1. `P_1/P_0` is the **grouped outgoing-prefactor slope** carried from Steps 4–6.
2. `P_r` is the **actual static numerator** of outgoing port `r`.

This step keeps both because the whole point is to prove that the former collapses to a weighted slope law for the latter.

---

## Inputs carried forward

### From Step 6

On the conservative-shape-preserving branch,
```math
\Theta_B=0,
\qquad
\Theta_Z=0,
```
so the surviving grouped defect is
```math
\Xi_1
=
2\sum_r \rho_r^{(N)}\,\delta\ln\Lambda_r-\kappa_1,
\qquad
\sum_r \rho_r^{(N)}=1.
```
Here `\kappa_1` is the weak-axisymmetric wall-baseline slope,
```math
\delta\ln K_A=\epsilon\lambda_A\kappa_1,
```
with grouped branch weights
```math
\lambda_{20}=1,
\qquad
\lambda_{21}=\frac12,
\qquad
\lambda_{22}=-1.
```

### From the moving-throat port dictionary

The actual static outgoing-port data are
```math
P_r=\Omega_{U,r}^2 G_{W,r}+R_r G_{U,r},
```
```math
\Delta_r=\Omega_{U,r}^2\Omega_{W,r}^2-R_r^2,
```
```math
N_0^{(r)}=\frac{P_r^2}{\Delta_r^2}.
```
So the Step-6 load factor is simply
```math
\Lambda_r^2=N_0^{(r)}.
```
That is the bridge used below.

---

## Step 7A — Weak-axisymmetric slopes of the actual outgoing-port data

Define the weak-axisymmetric slopes of the actual numerator and actual detuning by
```math
\delta\ln P_{A,r}=\epsilon\lambda_A\,\mathfrak p_r,
\qquad
\delta\ln \Delta_{A,r}=\epsilon\lambda_A\,\mathfrak d_r.
```
Then the actual static outgoing-transfer coefficient satisfies
```math
N_{A,0}^{(r)}=\frac{P_{A,r}^2}{\Delta_{A,r}^2},
```
so its weak-axisymmetric logarithmic slope is
```math
\boxed{
\delta\ln N_{A,0}^{(r)}=\epsilon\lambda_A\,\nu_r,
\qquad
\nu_r:=2(\mathfrak p_r-\mathfrak d_r).
}
```
This is the first exact identity of the step.

So the wall-referenced outgoing-load defect of one port is just
```math
\delta\ln\!\left(\frac{\Lambda_{A,r}^2}{K_A}\right)
=
\delta\ln N_{A,0}^{(r)}-\delta\ln K_A
=
\epsilon\lambda_A(\nu_r-\kappa_1).
```

---

## Step 7B — Exact collapse of the full remaining grouped defect

Because `\Lambda_{A,r}^2=N_{A,0}^{(r)}`, the Step-6 theorem becomes
```math
\Xi_{\rm load}^{(A)}
=
\sum_r \rho_r^{(N)}
\,\delta\ln\!\left(\frac{\Lambda_{A,r}^2}{K_A}\right)
=
\epsilon\lambda_A
\sum_r \rho_r^{(N)}(\nu_r-\kappa_1).
```
So the scalar carried by the grouped branch is
```math
\boxed{
\Xi_1
=
\sum_r \rho_r^{(N)}(\nu_r-\kappa_1).
}
```
Define the outgoing-weighted static transfer slope
```math
\bar\nu_N:=\sum_r \rho_r^{(N)}\nu_r.
```
Since the weights sum to one, the whole defect collapses to
```math
\boxed{
\Xi_1=\bar\nu_N-\kappa_1.
}
```
This is the main theorem of the step.

It says that the whole remaining linear grouped `2.5`PN / quartic-anomaly defect is exactly the mismatch between

- the outgoing-weighted static transfer slope of the **actual** moving-throat ports,
- and the wall-baseline slope.

So the remaining g-2 theorem problem is no longer “compute every microscopic slippage variable.”
It is just:

> **do the actual outgoing ports co-load with the wall baseline?**

---

## Step 7C — Exact formula for the actual numerator slope `\mathfrak p_r`

The actual static numerator is
```math
P_r=\Omega_{U,r}^2 G_{W,r}+R_r G_{U,r}.
```
Define the positive numerator weights
```math
\boxed{
\alpha_r:=\frac{\Omega_{U,r}^2 G_{W,r}}{P_r},
\qquad
\beta_r:=\frac{R_r G_{U,r}}{P_r},
\qquad
\alpha_r+\beta_r=1.
}
```
Then the weak-axisymmetric numerator slope is exactly
```math
\boxed{
\mathfrak p_r
=
\alpha_r(\mathfrak o_{U,r}+\mathfrak g_{W,r})
+
\beta_r(\mathfrak r_r+\mathfrak g_{U,r}).
}
```
So `\mathfrak p_r` is the convex average of the two static numerator legs:

- the brane-like leg `\Omega_{U,r}^2 G_{W,r}`,
- the mixed leg `R_r G_{U,r}`.

---

## Step 7D — Exact formula for the actual detuning slope `\mathfrak d_r`

The actual static detuning is
```math
\Delta_r=\Omega_{U,r}^2\Omega_{W,r}^2-R_r^2.
```
Define the detuning weights
```math
\boxed{
\chi_r:=\frac{\Omega_{U,r}^2\Omega_{W,r}^2}{\Delta_r},
\qquad
\zeta_r:=\frac{R_r^2}{\Delta_r}.
}
```
Equivalently, with
```math
\mathcal H_r:=\frac{R_r^2}{\Omega_{U,r}^2\Omega_{W,r}^2},
```
one has
```math
\chi_r=\frac{1}{1-\mathcal H_r},
\qquad
\zeta_r=\frac{\mathcal H_r}{1-\mathcal H_r}.
```
Then the weak-axisymmetric detuning slope is exactly
```math
\boxed{
\mathfrak d_r
=
\chi_r(\mathfrak o_{U,r}+\mathfrak o_{W,r})
-
2\zeta_r\,\mathfrak r_r.
}
```
So `\mathfrak d_r` measures how the combined internal frequencies and the coupling `R_r` change the static port detuning.

---

## Step 7E — Static outgoing-transfer slope in actual port variables

Combining the last two results with `\nu_r=2(\mathfrak p_r-\mathfrak d_r)` gives
```math
\boxed{
\nu_r
=
2\alpha_r(\mathfrak o_{U,r}+\mathfrak g_{W,r})
+
2\beta_r(\mathfrak r_r+\mathfrak g_{U,r})
-
2\chi_r(\mathfrak o_{U,r}+\mathfrak o_{W,r})
+
4\zeta_r\,\mathfrak r_r.
}
```
This is the first direct formula for the static outgoing-transfer slope of an actual moving-throat port.

So the exact port-level g-2 closure is now fully explicit in the real port data.

---

## Step 7F — Exact equivalence to the earlier slippage language

Earlier steps used the port slippages
```math
\mathfrak m_r:=\mathfrak g_{W,r}-\mathfrak o_{W,r}-\frac12\kappa_1,
```
```math
\mathfrak i_r:=\mathfrak r_r+\mathfrak g_{U,r}-\mathfrak o_{U,r}-\mathfrak g_{W,r},
```
```math
\mathfrak h_r:=2\mathfrak r_r-\mathfrak o_{U,r}-\mathfrak o_{W,r},
```
with
```math
\mathcal I_r:=\frac{R_r G_{U,r}}{\Omega_{U,r}^2 G_{W,r}},
\qquad
\mathcal H_r:=\frac{R_r^2}{\Omega_{U,r}^2\Omega_{W,r}^2}.
```
Because
```math
\alpha_r=\frac{1}{1+\mathcal I_r},
\qquad
\beta_r=\frac{\mathcal I_r}{1+\mathcal I_r},
\qquad
\chi_r=\frac{1}{1-\mathcal H_r},
\qquad
\zeta_r=\frac{\mathcal H_r}{1-\mathcal H_r},
```
one gets the exact identity
```math
\boxed{
\nu_r
=
\kappa_1
+
2\mathfrak m_r
+
\frac{2\mathcal I_r}{1+\mathcal I_r}\,\mathfrak i_r
+
\frac{2\mathcal H_r}{1-\mathcal H_r}\,\mathfrak h_r.
}
```
So the earlier port amplitude `\sigma_r` is simply
```math
\boxed{
\sigma_r=\nu_r-\kappa_1.
}
```
This proves the present theorem is not a new branch of the algebra. It is the direct outgoing-port rewrite of the earlier slippage theorem.

---

## Step 7G — Outgoing-port co-loading theorem

### Exact zero-defect condition

The remaining linear grouped defect vanishes if and only if
```math
\boxed{
\bar\nu_N=\kappa_1.
}
```
Equivalently,
```math
\boxed{
\frac{P_1}{P_0}=0
\quad\Longleftrightarrow\quad
\bar\nu_N=\kappa_1.
}
```
So the outgoing-weighted static transfer slope must match the wall-baseline slope.

### Strong per-port sufficient condition

A stronger sufficient condition is
```math
\boxed{
\nu_r=\kappa_1
\qquad
\text{for every active outgoing port }r.
}
```
Then every static outgoing-transfer coefficient co-loads with the wall baseline lane by lane, and therefore `\Xi_1=0`.

### Dominant-port limit

If one outgoing port dominates,
```math
\rho_{r_*}^{(N)}\approx1,
```
then
```math
\boxed{
\Xi_1\approx \nu_{r_*}-\kappa_1.
}
```
So the last linear grouped defect is just the mismatch between the dominant port slope and the wall-baseline slope.

### Naive rigidity no-go

If the actual ports are rigid,
```math
\mathfrak p_r=\mathfrak d_r=0
\qquad\Longrightarrow\qquad
\nu_r=0,
```
then
```math
\boxed{
\Xi_1=-\kappa_1.
}
```
So freezing the outgoing ports is **not** enough. They must actively co-load with the wall baseline.

---

## Step 7H — Direct quartic anomaly law

The carried quartic target is
```math
\Xi_1=\Lambda_1,
\qquad
\Lambda_1\approx 0.279605891931464.
```
So the exact port-level anomaly gate is now
```math
\boxed{
\bar\nu_N=\kappa_1+\Lambda_1.
}
```
In the dominant-port regime this becomes
```math
\boxed{
\nu_{r_*}=\kappa_1+\Lambda_1.
}
```
This is the cleanest g-2 continuation reached so far:

> the remaining quartic anomaly correction is exactly the amount by which the outgoing-weighted static transfer slope must exceed the wall-baseline slope.

---

## Step 7I — Reduced verdict

Step 6 said the whole remaining defect is an outgoing-load law for `\Lambda_r=P_r/\Delta_r`.
Step 7 sharpens that again.

It proves that:

1. the actual static outgoing-transfer slope of each port is
   ```math
   \nu_r=2(\mathfrak p_r-\mathfrak d_r),
   ```
2. the full remaining grouped defect is exactly
   ```math
   \Xi_1=\bar\nu_N-\kappa_1,
   ```
3. the earlier slippage amplitude is just
   ```math
   \sigma_r=\nu_r-\kappa_1,
   ```
4. and anomaly matching requires
   ```math
   \bar\nu_N=\kappa_1+\Lambda_1.
   ```

So the next theorem gate is smaller again.

It is no longer
> “compute the whole microscopic slippage bundle.”

It is now simply
> **compute the actual outgoing-port slopes `\nu_r` on the moving-throat branch and see whether their outgoing-weighted average lands at `\kappa_1+\Lambda_1`.**
# Step 8 — Wall-normalized transfer-shape theorem

## Goal

Step 7 reduced the remaining quartic anomaly closure to the actual outgoing-port co-loading law
```math
\Xi_1
=
\bar\nu_N-\kappa_1,
\qquad
\bar\nu_N:=\sum_r \rho_r^{(N)}\nu_r,
```
where
```math
\nu_r
=
\frac{\delta\ln N_{A,0}^{(r)}}{\epsilon\lambda_A},
\qquad
N_{A,0}^{(r)}=\frac{P_{A,r}^2}{\Delta_{A,r}^2}.
```
That was already sharp, but it still phrased the theorem gate in terms of the raw actual-port pair
```math
P_r,
\qquad
\Delta_r.
```

The next honest question is therefore:

> **does each actual outgoing port admit a direct wall-normalized factorization, so that the remaining grouped defect becomes nothing more than the weak-axisymmetric drift of one dimensionless transfer shape?**

This step answers that.

---

## Inputs carried forward

### From Step 7

For each actual outgoing port,
```math
P_r=\Omega_{U,r}^2 G_{W,r}+R_r G_{U,r},
```
```math
\Delta_r=\Omega_{U,r}^2\Omega_{W,r}^2-R_r^2,
```
```math
N_0^{(r)}=\frac{P_r^2}{\Delta_r^2}.
```
On the weak-axisymmetric grouped branch,
```math
\delta\ln K_A=\epsilon\lambda_A\kappa_1,
\qquad
\lambda_{20}=1,
\quad
\lambda_{21}=\frac12,
\quad
\lambda_{22}=-1,
```
and Step 7 proved
```math
\Xi_1
=
\sum_r \rho_r^{(N)}(\nu_r-\kappa_1),
\qquad
\nu_r
=
\frac{\delta\ln N_{A,0}^{(r)}}{\epsilon\lambda_A}.
```

### From the moving-throat port dictionary

The wall baseline is the static grouped prefactor carrier `K_A`, so the natural next move is to normalize every actual outgoing quantity by `K_A`, `\Omega_{U,r}`, and `\Omega_{W,r}` before taking grouped weak-axisymmetric slopes.

---

## Step 8A — Exact wall-normalized factorization of the actual outgoing-transfer coefficient

Introduce the wall-normalized dimensionless port variables
```math
\boxed{
\widehat G_{W,r}:=\frac{G_{W,r}}{\Omega_{W,r}^2\sqrt K},
\qquad
\widehat G_{U,r}:=\frac{G_{U,r}}{\Omega_{U,r}\Omega_{W,r}\sqrt K},
\qquad
\widehat R_r:=\frac{R_r}{\Omega_{U,r}\Omega_{W,r}}.
}
```
Then the actual numerator and detuning factor exactly as
```math
P_r
=
\sqrt K\,\Omega_{U,r}^2\Omega_{W,r}^2
\bigl(\widehat G_{W,r}+\widehat R_r\widehat G_{U,r}\bigr),
```
```math
\Delta_r
=
\Omega_{U,r}^2\Omega_{W,r}^2(1-\widehat R_r^2).
```
Therefore
```math
\boxed{
\frac{N_0^{(r)}}{K}
=
\left[
\frac{\widehat G_{W,r}+\widehat R_r\widehat G_{U,r}}{1-\widehat R_r^2}
\right]^2.
}
```
Define the wall-normalized transfer shape
```math
\boxed{
\mathcal T_r
:=
\frac{\widehat G_{W,r}+\widehat R_r\widehat G_{U,r}}{1-\widehat R_r^2},
}
```
so that the actual static outgoing-transfer coefficient factors exactly as
```math
\boxed{
N_0^{(r)}=K\,\mathcal T_r^2.
}
```

This is the first main result of the step.
The actual outgoing-transfer coefficient is the wall baseline times one dimensionless transfer shape squared.

---

## Step 8B — Weak-axisymmetric transport of the wall-normalized port variables

Define the weak-axisymmetric grouped slopes of the wall-normalized port variables by
```math
\boxed{
\delta\ln\widehat G_{W,A,r}=\epsilon\lambda_A\,\mathfrak w_r,
}
```
```math
\boxed{
\delta\ln\widehat G_{U,A,r}=\epsilon\lambda_A\,\mathfrak u_r,
}
```
```math
\boxed{
\delta\ln\widehat R_{A,r}=\epsilon\lambda_A\,\mathfrak c_r.
}
```
In terms of the primitive Step-7 slopes,
```math
\boxed{
\mathfrak w_r
=
\mathfrak g_{W,r}-\mathfrak o_{W,r}-\frac12\kappa_1,
}
```
```math
\boxed{
\mathfrak u_r
=
\mathfrak g_{U,r}-\frac12\mathfrak o_{U,r}-\frac12\mathfrak o_{W,r}-\frac12\kappa_1,
}
```
```math
\boxed{
\mathfrak c_r
=
\mathfrak r_r-\frac12\mathfrak o_{U,r}-\frac12\mathfrak o_{W,r}.
}
```
So the wall-normalized variables already strip out the trivial wall-baseline scaling and leave only the genuine dimensionless port-shape drifts.

---

## Step 8C — Exact transfer-shape slope and the identity `\nu_r=\kappa_1+2\tau_r`

Because
```math
\mathcal T_r
=
\frac{\widehat G_{W,r}+\widehat R_r\widehat G_{U,r}}{1-\widehat R_r^2},
```
its weak-axisymmetric logarithmic slope is
```math
\boxed{
\delta\ln\mathcal T_{A,r}=\epsilon\lambda_A\,\tau_r,
}
```
with
```math
\boxed{
\tau_r
=
\widehat\alpha_r\,\mathfrak w_r
+
\widehat\beta_r\,(\mathfrak u_r+\mathfrak c_r)
+
\frac{2\widehat R_r^2}{1-\widehat R_r^2}\,\mathfrak c_r,
}
```
where
```math
\boxed{
\widehat\alpha_r:=
\frac{\widehat G_{W,r}}{\widehat G_{W,r}+\widehat R_r\widehat G_{U,r}},
\qquad
\widehat\beta_r:=
\frac{\widehat R_r\widehat G_{U,r}}{\widehat G_{W,r}+\widehat R_r\widehat G_{U,r}},
\qquad
\widehat\alpha_r+\widehat\beta_r=1.
}
```

Now use the exact factorization
```math
N_{A,0}^{(r)}=K_A\mathcal T_{A,r}^2.
```
Taking the weak-axisymmetric logarithmic slope gives the central identity
```math
\boxed{
\nu_r
=
\frac{\delta\ln N_{A,0}^{(r)}}{\epsilon\lambda_A}
=
\kappa_1+2\tau_r.
}
```
So the actual static outgoing-transfer slope is just the wall-baseline slope plus twice the transfer-shape slope.

---

## Step 8D — Exact collapse of the remaining grouped defect

Step 7 already gave
```math
\Xi_1
=
\sum_r \rho_r^{(N)}(\nu_r-\kappa_1),
\qquad
\sum_r \rho_r^{(N)}=1.
```
Substituting the transfer-shape identity above yields
```math
\boxed{
\Xi_1
=
2\sum_r \rho_r^{(N)}\tau_r.
}
```
Equivalently, the exact zero-defect condition becomes
```math
\boxed{
\sum_r \rho_r^{(N)}\tau_r=0.
}
```
This is the sharpest reduced theorem gate reached so far.

A stronger per-port sufficient condition is
```math
\boxed{
\tau_r=0
\qquad\text{for every active outgoing port }r,
}
```
which is equivalent to
```math
\boxed{
\delta\ln\mathcal T_{A,r}=0
\qquad\text{for every active outgoing port }r.
}
```
So the exact reduced meaning of port co-loading is now:

> **each wall-normalized transfer shape must be weak-axisymmetrically rigid.**

---

## Step 8E — Exact equivalence to the earlier slippage languages

Earlier steps used the slippage variables
```math
\mathfrak m_r,
\qquad
\mathfrak i_r,
\qquad
\mathfrak h_r,
```
with the port-shape composites
```math
\mathcal M_r:=\frac{G_{W,r}}{\Omega_{W,r}^2\sqrt K},
\qquad
\mathcal I_r:=\frac{R_r G_{U,r}}{\Omega_{U,r}^2 G_{W,r}},
\qquad
\mathcal H_r:=\frac{R_r^2}{\Omega_{U,r}^2\Omega_{W,r}^2}.
```
These are related to the present variables by
```math
\mathcal M_r=\widehat G_{W,r},
\qquad
\mathcal I_r=\frac{\widehat R_r\widehat G_{U,r}}{\widehat G_{W,r}},
\qquad
\mathcal H_r=\widehat R_r^2,
```
and
```math
\mathfrak m_r=\mathfrak w_r,
\qquad
\mathfrak i_r=(\mathfrak u_r+\mathfrak c_r)-\mathfrak w_r,
\qquad
\mathfrak h_r=2\mathfrak c_r.
```
With these substitutions the transfer-shape slope becomes exactly
```math
\boxed{
\tau_r
=
\mathfrak m_r
+
\frac{\mathcal I_r}{1+\mathcal I_r}\,\mathfrak i_r
+
\frac{\mathcal H_r}{1-\mathcal H_r}\,\mathfrak h_r.
}
```
So the earlier Step-7 port amplitude is simply
```math
\boxed{
\sigma_r=2\tau_r.
}
```
This proves the present theorem is not a different branch of the algebra. It is the exact compressed form of the Stage-159/160/161 slippage structure.

---

## Step 8F — Direct quartic anomaly law in transfer-shape language

The carried quartic target is still
```math
\Xi_1=\Lambda_1,
\qquad
\Lambda_1\approx 0.279605891931464.
```
Since
```math
\Xi_1=2\sum_r \rho_r^{(N)}\tau_r,
```
the exact transfer-shape anomaly gate is now
```math
\boxed{
\sum_r \rho_r^{(N)}\tau_r
=
\frac{\Lambda_1}{2}
\approx 0.139802945965732.
}
```
In the dominant-port limit,
```math
\boxed{
\tau_{r_*}
=
\frac{\Lambda_1}{2}
\approx 0.139802945965732.
}
```
So the remaining quartic anomaly correction is exactly the amount by which the outgoing-weighted transfer shape must drift away from weak-axisymmetric rigidity.

---

## Step 8G — Reduced verdict

Step 7 said the remaining quartic anomaly layer is a mismatch between the outgoing-weighted actual transfer slope and the wall-baseline slope.
Step 8 sharpens that one step further.

It proves that:

1. every actual outgoing-transfer coefficient factors exactly as
   ```math
   N_0^{(r)}=K\,\mathcal T_r^2,
   ```
2. the actual port slope is
   ```math
   \nu_r=\kappa_1+2\tau_r,
   ```
3. the full remaining grouped defect is exactly
   ```math
   \Xi_1=2\sum_r \rho_r^{(N)}\tau_r,
   ```
4. and anomaly matching requires
   ```math
   \sum_r \rho_r^{(N)}\tau_r=\Lambda_1/2.
   ```

So the next theorem gate is smaller again.

It is no longer
> “compute the raw outgoing-port slopes `\nu_r`.”

It is now simply
> **compute the weak-axisymmetric wall-normalized transfer-shape slopes `\tau_r` on the actual moving-throat branch, and test whether their outgoing-weighted average lands at `\Lambda_1/2`.**

That is the direct continuation point.
# Step 9 — Effective transfer shape collapse and the one-port continuum branch

## Goal

Step 8 reduced the remaining quartic anomaly layer to the outgoing-weighted transfer-shape slope
```math
\Xi_1 = 2\sum_r \rho_r^{(N)}\tau_r,
\qquad
N_{A,0}^{(r)} = K_A\,\mathcal T_{A,r}^2.
```
That was already sharp, but it still carried the theorem gate as a weighted many-port average.

The next honest question is therefore:

> **does that weighted average collapse to the slope of one single effective transfer shape, and if so, what is that shape on the actual minimal continuum branch?**

This step answers that.

---

## Inputs carried forward

### From Step 8

For each active outgoing port,
```math
N_{A,0}^{(r)} = K_A\,\mathcal T_{A,r}^2,
```
with weak-axisymmetric slope
```math
\delta\ln \mathcal T_{A,r} = \epsilon\lambda_A\,\tau_r.
```
So the remaining grouped defect is already
```math
\Xi_1 = 2\sum_r \rho_r^{(N)}\tau_r,
\qquad
\rho_r^{(N)} = \frac{N_{A,0}^{(r)}}{\sum_s N_{A,0}^{(s)}}.
```

### From the moving-throat continuum branch

On the minimal isotropic continuum branch there is only one active outgoing port and the actual static transfer coefficient is written as
```math
N_{A,0}=\beta_{0,A},
\qquad
K_A=K_{0,A}.
```
The continuum formulas carried forward are
```math
K_{0,A} = \frac{K_{\eta,A}^{(\mathrm{eff})}}{\mu_{\eta,A}},
```
```math
\beta_{0,A}
=
\frac{\mu_{W,A}}{\mu_{\eta,A}}
\frac{K_{\eta,A}^{(\mathrm{eff})}}{K_{W,A}^{(\mathrm{eff})}}
\frac{Z_{W,A}(1+\rho_A)^2}{(1-\epsilon_{W,A})^2}.
```
These are the only branch inputs needed here.

---

## Step 9A — Exact collapse from many ports to one effective transfer shape

Summing the portwise factorization gives
```math
N_{A,0}
=
\sum_r N_{A,0}^{(r)}
=
K_A\sum_r \mathcal T_{A,r}^2.
```
So define the **effective wall-normalized transfer shape** by
```math
\boxed{
\mathcal T_{\mathrm{eff},A}^2
:=
\sum_r \mathcal T_{A,r}^2
=
\frac{N_{A,0}}{K_A}.
}
```
Now perturb weak-axisymmetrically:
```math
\delta\ln \mathcal T_{A,r} = \epsilon\lambda_A\,\tau_r.
```
Then
```math
\frac{\delta\ln \mathcal T_{\mathrm{eff},A}^2}{\epsilon\lambda_A}
=
\frac{2\sum_r \mathcal T_r^2\tau_r}{\sum_s \mathcal T_s^2}
=
2\sum_r \rho_r^{(N)}\tau_r,
```
because
```math
\rho_r^{(N)}
=
\frac{N_0^{(r)}}{N_0}
=
\frac{\mathcal T_r^2}{\sum_s \mathcal T_s^2}.
```
Therefore
```math
\boxed{
\Xi_1
=
\frac{\delta\ln \mathcal T_{\mathrm{eff},A}^2}{\epsilon\lambda_A}.
}
```
Equivalently, if
```math
\delta\ln \mathcal T_{\mathrm{eff},A} = \epsilon\lambda_A\,\tau_{\mathrm{eff}},
```
then
```math
\boxed{
\tau_{\mathrm{eff}}=\sum_r \rho_r^{(N)}\tau_r,
\qquad
\Xi_1 = 2\tau_{\mathrm{eff}}.
}
```

So the whole remaining grouped defect is already the slope of **one** effective transfer shape.

---

## Step 9B — Exact one-port continuum transfer shape

On the actual minimal isotropic continuum branch there is only one active outgoing port, so
```math
\mathcal T_{\mathrm{eff},A}=\mathcal T_A,
\qquad
N_{A,0}=\beta_{0,A},
\qquad
K_A=K_{0,A}.
```
Hence
```math
\mathcal T_A^2 = \frac{\beta_{0,A}}{K_{0,A}}.
```
Substituting the continuum formulas gives immediately
```math
\boxed{
\mathcal T_A^2
=
\frac{\mu_{W,A}}{K_{W,A}^{(\mathrm{eff})}}
\frac{Z_{W,A}(1+\rho_A)^2}{(1-\epsilon_{W,A})^2}.
}
```
Using
```math
\Omega_{W,A}^2 = \frac{K_{W,A}^{(\mathrm{eff})}}{\mu_{W,A}},
```
this becomes
```math
\boxed{
\mathcal T_A^2
=
\frac{Z_{W,A}(1+\rho_A)^2}{\Omega_{W,A}^2(1-\epsilon_{W,A})^2}.
}
```

This is the first exact actual-port formula for the effective transfer shape on the real continuum branch.

A useful physical consequence is immediate:

**the wall–U dressing ratio `\epsilon_\eta` drops out of the direct continuum transfer-shape formula.**

At the direct port level the defect only sees

- the wall-to-mixed overlap ratio `Z_W`,
- the interference ratio `\rho`,
- the mixed-sector blocking ratio `\epsilon_W`,
- the mixed frequency scale `\Omega_W`.

So the wall–U dressing lane affects the port transfer only indirectly.

---

## Step 9C — Direct weak-axisymmetric slope law in continuum variables

Perturb the one-port continuum formula weak-axisymmetrically by
```math
\delta\ln Z_{W,A} = \epsilon\lambda_A\,\zeta_W,
```
```math
\delta\ln \Omega_{W,A}^2 = \epsilon\lambda_A\,\omega_W,
```
```math
\delta \rho_A = \epsilon\lambda_A\,\rho_1,
```
```math
\delta\epsilon_{W,A} = \epsilon\lambda_A\,\varepsilon_W.
```
Then from
```math
\mathcal T_A^2
=
\frac{Z_{W,A}(1+\rho_A)^2}{\Omega_{W,A}^2(1-\epsilon_{W,A})^2},
```
we obtain the exact one-port weak-axisymmetric defect law
```math
\boxed{
\Xi_1
=
\frac{\delta\ln \mathcal T_A^2}{\epsilon\lambda_A}
=
\zeta_W - \omega_W + \frac{2\rho_1}{1+\rho} + \frac{2\varepsilon_W}{1-\epsilon_W}.
}
```
Equivalently,
```math
\boxed{
\tau
=
\frac12\,\Xi_1
=
\frac12\left(
\zeta_W - \omega_W + \frac{2\rho_1}{1+\rho} + \frac{2\varepsilon_W}{1-\epsilon_W}
\right).
}
```

So on the actual minimal continuum port, the remaining grouped defect is controlled by **four** direct microscopic drift channels:

1. wall-to-mixed overlap drift `\zeta_W`,
2. mixed-frequency drift `\omega_W`,
3. interference-ratio drift `\rho_1`,
4. mixed blocking drift `\varepsilon_W`.

This is the first exact actual-port answer to the Step-8 “compute `\tau_r`” question.

---

## Step 9D — Selected-branch reformulation

The same one-port transfer shape can be rewritten in terms of the selected-branch placement variables.

The carried demand-ratio formula is
```math
R_{\mathrm{target},A}
=
\Lambda_A
\frac{(1-\epsilon_{\eta,A})(1-\epsilon_{W,A})^2}{Z_{W,A}(1+\rho_A)^2},
```
with
```math
\Lambda_A
=
\frac{27\pi^2 G c_s^5 K_{W,A}^{(\mathrm{eff})}}{20 a^5 c^5 \mu_{W,A}}
=
\frac{27\pi^2 G c_s^5}{20 a^5 c^5}\,\Omega_{W,A}^2.
```
So the transfer shape becomes
```math
\boxed{
\mathcal T_A^2
=
\frac{27\pi^2 G c_s^5}{20 a^5 c^5}
\frac{1-\epsilon_{\eta,A}}{R_{\mathrm{target},A}}.
}
```

At linear grouped weak-axisymmetric order, scalar observables such as `a` and `c_s` do not shift on the isotropic branch, so the front factor is inert at this order. If we define
```math
\delta\epsilon_{\eta,A} = \epsilon\lambda_A\,\eta_1,
\qquad
\delta\ln R_{\mathrm{target},A} = \epsilon\lambda_A\,\mathcal R_1,
```
then the same defect becomes
```math
\boxed{
\Xi_1
=
-\frac{\eta_1}{1-\epsilon_\eta} - \mathcal R_1.
}
```

So the one-port defect has a second exact interpretation:

it is the mismatch between

- the wall–U dressing drift, and
- the selected-branch demand-ratio drift.

---

## Step 9E — Exact one-port zero-defect theorem

The direct zero-defect condition is now explicit:
```math
\boxed{
\delta\ln\!\left[
\frac{Z_W(1+\rho)^2}{\Omega_W^2(1-\epsilon_W)^2}
\right] = 0.
}
```
Equivalently, in selected-branch variables,
```math
\boxed{
\delta\ln R_{\mathrm{target}}
=
-\frac{\delta\epsilon_\eta}{1-\epsilon_\eta}.
}
```
So the one-port continuum theorem is genuinely bidirectional.

Two useful corollaries follow immediately.

### Target-rigid branch

If the selected-branch demand ratio is weak-axisymmetrically rigid,
```math
\delta\ln R_{\mathrm{target}}=0,
```
then zero grouped defect requires
```math
\boxed{\delta\epsilon_\eta = 0.}
```
So on a target-rigid branch, the wall–U dressing ratio itself must be lane-rigid.

### Wall–U-rigid branch

If the wall–U dressing ratio is weak-axisymmetrically rigid,
```math
\delta\epsilon_\eta = 0,
```
then zero grouped defect requires
```math
\boxed{\delta\ln R_{\mathrm{target}}=0.}
```
So on an `\epsilon_\eta`-rigid branch, the selected demand ratio must be lane-rigid.

---

## Step 9F — Direct quartic anomaly gate on the actual one-port branch

The carried quartic anomaly target is still
```math
\Xi_1 = \Lambda_1,
\qquad
\Lambda_1 \approx 0.279605891931464.
```
So the exact one-port direct-continuum anomaly condition is
```math
\boxed{
\zeta_W - \omega_W + \frac{2\rho_1}{1+\rho} + \frac{2\varepsilon_W}{1-\epsilon_W}
=
\Lambda_1.
}
```
Equivalently,
```math
\boxed{
\tau = \frac{\Lambda_1}{2} \approx 0.139802945965732.
}
```
The selected-branch form is just as sharp:
```math
\boxed{
-\frac{\eta_1}{1-\epsilon_\eta} - \mathcal R_1 = \Lambda_1.
}
```

A useful **reference balanced split** in the direct continuum variables is obtained by sharing the defect equally among the four rescaled drift channels. That gives
```math
\boxed{
\zeta_W = \frac{\Lambda_1}{4},
\qquad
\omega_W = -\frac{\Lambda_1}{4},
\qquad
\rho_1 = \frac{1+\rho}{8}\,\Lambda_1,
\qquad
\varepsilon_W = \frac{1-\epsilon_W}{8}\,\Lambda_1.
}
```
Numerically,
```math
\frac{\Lambda_1}{4} \approx 0.0699014729828660.
```
A corresponding balanced split in the selected-branch variables is
```math
\boxed{
\eta_1 = -\frac{1-\epsilon_\eta}{2}\,\Lambda_1,
\qquad
\mathcal R_1 = -\frac{\Lambda_1}{2}.
}
```
Numerically,
```math
-\frac{\Lambda_1}{2} \approx -0.139802945965732.
```
These are not unique solutions; they are just the cleanest reference branch targets.

---

## Step 9G — Consequence for the grouped weak-axisymmetric quadrupole defect

The grouped weak-axisymmetric pattern is still
```math
\Delta_Q^{(20)} = \epsilon\,\Xi_1,
\qquad
\Delta_Q^{(21)} = \frac{\epsilon}{2}\,\Xi_1,
\qquad
\Delta_Q^{(22)} = -\epsilon\,\Xi_1.
```
Because `\Xi_1 = 2\tau` on the one-port branch, this becomes
```math
\boxed{
\Delta_Q^{(20)} = 2\epsilon\,\tau,
\qquad
\Delta_Q^{(21)} = \epsilon\,\tau,
\qquad
\Delta_Q^{(22)} = -2\epsilon\,\tau.
}
```
So once the actual one-port transfer-shape slope is known, the full weak-axisymmetric grouped normalization pattern is fixed immediately.

---

## Reduced verdict

Step 8 said the remaining quartic anomaly layer is the outgoing-weighted transfer-shape slope.
Step 9 sharpens that result in two decisive ways.

First, the many-port weighted average is itself exactly the slope of a **single effective transfer shape**,
```math
\mathcal T_{\mathrm{eff}}^2 = \frac{N_0}{K}.
```

Second, on the actual minimal one-port continuum branch, that transfer shape is no longer abstract. It is explicitly
```math
\mathcal T^2
=
\frac{Z_W(1+\rho)^2}{\Omega_W^2(1-\epsilon_W)^2}
=
\frac{27\pi^2 G c_s^5}{20 a^5 c^5}
\frac{1-\epsilon_\eta}{R_{\mathrm{target}}}.
```

So the remaining theorem gate is no longer
> “compute all raw outgoing-port slopes.”

It is now
> **determine whether the actual weak-axisymmetric moving-throat branch keeps this single continuum transfer shape rigid, and if not, which of its four direct microscopic drift channels carries the required `\Lambda_1`-sized defect.**

That is the direct continuation point.
# Step 10 — Coherent tracking-branch substitution and the support-blindness theorem

## Goal

Step 9 reduced the remaining quartic anomaly layer to the weak-axisymmetric slope of a **single** one-port continuum transfer shape,
```math
\mathcal T^2
=
\frac{Z_W(1+\rho)^2}{\Omega_W^2(1-\epsilon_W)^2}
=
\frac{27\pi^2 G c_s^5}{20 a^5 c^5}\,\frac{1-\epsilon_\eta}{R_{\mathrm{target}}},
\qquad
\Xi_1 = \delta\ln \mathcal T^2.
```
The next honest move is therefore not to invent a new transport layer by hand. It is to substitute the **actual coherent local D/N tracking branch** into that one-port law and see what survives.

That is what this step does.

The decisive outcome is:

> **coherent support can raise the steady normalization baseline, but it drops out identically from the first grouped weak-axisymmetric defect.**

So the next live g-2 calculation is no longer “how much coherent support is present?” It is “how the mixed/outgoing placement variables drift across the grouped weak-axisymmetric branch.”

---

## Inputs carried forward

### From Step 9

The one-port continuum transfer-shape law is
```math
\mathcal T^2
=
\frac{Z_W(1+\rho)^2}{\Omega_W^2(1-\epsilon_W)^2}
=
\frac{27\pi^2 G c_s^5}{20 a^5 c^5}\,\frac{1-\epsilon_\eta}{R_{\mathrm{target}}}.
```
Its grouped defect is
```math
\Xi_1 = \delta\ln \mathcal T^2,
```
and the quartic anomaly target remains
```math
\Xi_1 = \Lambda_1,
\qquad
\Lambda_1 \approx 0.279605891931464,
\qquad
\tau = \frac{\Lambda_1}{2} \approx 0.139802945965732.
```

### From the coherent local D/N tracking branch

The exact coherent branch data are
```math
R_{\mathrm{tr}}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0},
```
```math
\epsilon
=
\epsilon_W^{(\mathrm{split})}
=
\epsilon_W\!\left[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right],
```
```math
M_{\mathrm{tr}} = M_{\mathrm{mix}}\,S(\zeta;\epsilon),
\qquad
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
```
```math
R_{\mathrm{target}}
=
\Lambda\,
\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
\qquad
\Lambda = \frac{27\pi^2 G c_s^5}{20 a^5 c^5}\,\Omega_W^2.
```

---

## Step 10A — Exact coherent tracking-branch substitution

Substituting the coherent-branch demand ratio into the selected-branch transfer-shape law gives immediately
```math
\mathcal T^2
=
\frac{27\pi^2 G c_s^5}{20 a^5 c^5}\,\frac{1-\epsilon_\eta}{R_{\mathrm{target}}}
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2}.
```
So on the actual coherent tracking branch, the direct transfer shape depends only on

- the wall-to-mixed overlap ratio `Z_W`,
- the mixed frequency scale `\Omega_W`,
- the common interference ratio `\chi_0`,
- the split mixed blocking ratio `\epsilon`.

The support-enhancement factor `S(\zeta;\epsilon)` does **not** appear.

---

## Step 10B — Exact support-blindness theorem

The coherent support lane changes the total baseline through
```math
M_{\mathrm{tr}} = M_{\mathrm{mix}}\,S(\zeta;\epsilon),
\qquad
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
```
with
```math
\frac{\partial S}{\partial \zeta} = \frac{1-\epsilon}{(1-\zeta\epsilon)^2} > 0
```
on the physical branch. So support is a strict baseline enhancer.

But the direct transfer shape itself is
```math
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2},
```
and the selected-branch demand ratio is
```math
R_{\mathrm{target}}
=
\Lambda\,
\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2}.
```
Neither contains `\zeta`. Therefore
```math
\boxed{
\frac{\partial \ln \mathcal T^2}{\partial \zeta}=0,
\qquad
\frac{\partial \ln R_{\mathrm{target}}}{\partial \zeta}=0,
\qquad
\frac{\partial \Xi_1}{\partial \zeta}=0.
}
```
This is the exact **support-blindness theorem**.

A useful contrast is that the steady normalization product still feels support:
```math
R_{\mathrm{target}} M_{\mathrm{tr}}
=
\frac{8\Lambda(1-\epsilon)}{\pi^2}
S(\zeta;\epsilon).
```
So support can help the steady normalization test, but it cannot repair or spoil the first grouped weak-axisymmetric defect.

---

## Step 10C — Exact weak-axisymmetric defect law in physical branch variables

To avoid confusion with the split blocking ratio `\epsilon`, write the grouped perturbation parameter as `s\lambda_A`.
Take
```math
\delta\ln Z_{W,A}=s\lambda_A\,\zeta_Z,
\qquad
\delta\ln \Omega_{W,A}^2=s\lambda_A\,\omega_W,
```
```math
\delta\chi_{0,A}=s\lambda_A\,\chi_1,
\qquad
\delta\epsilon_{W,A}=s\lambda_A\,\varepsilon_W,
\qquad
\delta\delta_{U,A}=s\lambda_A\,\delta_{U,1},
```
```math
\delta\epsilon_{\eta,A}=s\lambda_A\,\eta_1.
```
Because
```math
\epsilon
=
\epsilon_W\!\left[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right],
```
its weak-axisymmetric drift is
```math
\boxed{
\epsilon_1
=
\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)\varepsilon_W
-
\frac{2\epsilon_W}{11(1+\delta_U)^2}\,\delta_{U,1}.
}
```
Then
```math
\ln \mathcal T^2
=
\ln Z_W + 2\ln(1+\chi_0) - \ln \Omega_W^2 - 2\ln(1-\epsilon),
```
so the grouped defect becomes
```math
\boxed{
\Xi_1
=
\zeta_Z - \omega_W + \frac{2\chi_1}{1+\chi_0} + \frac{2\epsilon_1}{1-\epsilon}.
}
```
Expanding `\epsilon_1` gives the full physical-branch law
```math
\boxed{
\Xi_1
=
\zeta_Z - \omega_W + \frac{2\chi_1}{1+\chi_0}
+
\frac{2}{1-\epsilon}
\left[
\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)\varepsilon_W
-
\frac{2\epsilon_W}{11(1+\delta_U)^2}\,\delta_{U,1}
\right].
}
```
So the first grouped defect is carried only by the mixed/outgoing placement drifts

1. `\zeta_Z` — wall-to-mixed overlap drift,
2. `\omega_W` — mixed-frequency drift,
3. `\chi_1` — common interference-ratio drift,
4. `\varepsilon_W` — bare mixed-blocking drift,
5. `\delta_{U,1}` — split-`U` axial drift.

The coherent support ratio `\zeta` is absent.

---

## Step 10D — Selected-branch reformulation

The Step-9 selected-branch identity remains exact on the coherent branch:
```math
\Xi_1 = -\frac{\eta_1}{1-\epsilon_\eta} - \mathcal R_1,
\qquad
\mathcal R_1 := \delta\ln R_{\mathrm{target},A}.
```
Since
```math
R_{\mathrm{target}}
=
\Lambda\,
\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
```
its weak-axisymmetric drift is
```math
\boxed{
\mathcal R_1
=
\omega_W - \frac{\eta_1}{1-\epsilon_\eta}
- \zeta_Z - \frac{2\chi_1}{1+\chi_0} - \frac{2\epsilon_1}{1-\epsilon}.
}
```
So the direct-port and selected-branch descriptions still agree exactly.

---

## Step 10E — Tracking-factor drift is not sufficient

The coherent tracking factor is
```math
R_{\mathrm{tr}}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0}
=
\frac{1+\chi_0+\delta_U}{(1+\chi_0)(1+\delta_U)}.
```
Its weak-axisymmetric logarithmic drift is
```math
\boxed{
\Theta_1
=
-\frac{\chi_0(1+\chi_0)\,\delta_{U,1} + \delta_U(1+\delta_U)\,\chi_1}
{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}.
}
```
So `\Theta_1` depends only on `(\chi_1,\delta_{U,1})`, whereas `\Xi_1` still depends on `(\zeta_Z,\omega_W,\varepsilon_W)` as well.

A simple explicit slice is
```math
\chi_1 = 0,
\qquad
\delta_{U,1}=0.
```
Then
```math
\Theta_1 = 0,
```
but
```math
\Xi_1
=
\zeta_Z - \omega_W
+
\frac{2}{1-\epsilon}
\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)\varepsilon_W,
```
which need not vanish.

So exact tracking-factor rigidity is **not** sufficient to kill the grouped defect.

---

## Step 10F — Coherent-branch quartic anomaly gate

The carried quartic target is still
```math
\Xi_1 = \Lambda_1,
\qquad
\Lambda_1 \approx 0.279605891931464,
\qquad
\tau = \frac{\Lambda_1}{2} \approx 0.139802945965732.
```
So the coherent tracking branch must satisfy
```math
\boxed{
\zeta_Z - \omega_W + \frac{2\chi_1}{1+\chi_0}
+
\frac{2}{1-\epsilon}
\left[
\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)\varepsilon_W
-
\frac{2\epsilon_W}{11(1+\delta_U)^2}\,\delta_{U,1}
\right]
=
\Lambda_1.
}
```
This is the first exact coherent-branch quartic anomaly gate written entirely in physical branch variables.

---

## Reduced verdict

Step 9 said the remaining quartic anomaly layer lives in one actual one-port transfer shape.
Step 10 sharpens that result in the most useful possible way.

First, on the actual coherent local D/N tracking branch,
```math
\mathcal T^2 = \frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2},
```
with `\epsilon = \epsilon_W[1 - (2/11)\delta_U/(1+\delta_U)]`.

Second, the coherent support lane drops out identically from
```math
\mathcal T^2,
\qquad
R_{\mathrm{target}},
\qquad
\Xi_1.
```
So support can raise the steady baseline but cannot carry the first grouped defect.

Third, the full coherent weak-axisymmetric defect is now carried only by the mixed/outgoing placement drifts
```math
(\zeta_Z,\omega_W,\chi_1,\varepsilon_W,\delta_{U,1}),
```
and exact tracking rigidity by itself is not enough to kill it.

So the next clean move is no longer “analyze coherent support.”
It is:

> **push the coherent branch one layer deeper, down to the microscopic coherent-kernel slippages, and separate the true tracking and nontracking defect coordinates.**

That is the direct continuation point.
# Step 11 — Microscopic slippage decomposition and exact triangular normal form

## Goal

Step 10 proved that the coherent local D/N tracking branch is **support-blind** at the level of the first grouped weak-axisymmetric defect:
```math
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2},
\qquad
\Xi_1=\delta\ln \mathcal T^2,
```
with
```math
\epsilon
=
\epsilon_W\!\left[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right].
```
So the next honest move is not to revisit coherent support. It is to push the defect one layer deeper, down to the actual **microscopic coherent-kernel slippages**.

That is what this step does.

The main result is that the whole coherent weak-axisymmetric problem collapses to three branch-adapted scalars,
```math
\Sigma_{\rm tr},\qquad \Sigma_{\rm nt},\qquad \Sigma_\eta,
```
rather than a broad raw-placement drift ledger.

---

## Inputs carried forward

From the coherent branch,
```math
\chi_0=\frac{\gamma c_{\eta U}}{K_U},
\qquad
\delta_U=\frac{\pi^2 T_U}{L^2K_U},
\qquad
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}},
```
```math
\epsilon_W=\frac{\gamma^2\lambda_W^2\sigma}{K_UK_W^{(\mathrm{eff})}},
\qquad
\frac{Z_W}{\Omega_W^2}
=
\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}.
```
The grouped defect from Step 10 was
```math
\Xi_1
=
\Sigma_Z
+\frac{2\chi_0}{1+\chi_0}\,\Sigma_\chi
+\frac{2\epsilon_W}{1-\epsilon}
\left[
\frac{11+9\delta_U}{11(1+\delta_U)}\,\Sigma_\epsilon
-\frac{2\delta_U}{11(1+\delta_U)^2}\,\Sigma_\delta
\right],
```
and the tracking-factor drift came from
```math
R_{\rm tr}=\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0}.
```

---

## Step 11A — Direct microscopic log coordinates

Perturb the microscopic coherent-kernel couplings multiplicatively:
```math
\gamma_A=\gamma e^{s\lambda_A\gamma_1},
\quad
c_{\eta U,A}=c_{\eta U}e^{s\lambda_A c_1},
\quad
\lambda_{W,A}=\lambda_W e^{s\lambda_A\lambda_1},
```
```math
K_{U,A}=K_U e^{s\lambda_A\kappa_U},
\quad
K_{\eta,A}=K_\eta e^{s\lambda_A\kappa_\eta},
\quad
K_{W,A}=K_W e^{s\lambda_A\kappa_W},
```
```math
\mu_{W,A}=\mu_W e^{s\lambda_A\mu_1},
\qquad
T_{U,A}=T_U e^{s\lambda_A\tau_1}.
```
Then the direct logarithmic slippages are
```math
\boxed{\Sigma_\chi=\delta\ln\chi_0=\gamma_1+c_1-\kappa_U,}
```
```math
\boxed{\Sigma_\delta=\delta\ln\delta_U=\tau_1-\kappa_U,}
```
```math
\boxed{\Sigma_\eta=\delta\ln\epsilon_\eta=2c_1-\kappa_U-\kappa_\eta,}
```
```math
\boxed{\Sigma_\epsilon=\delta\ln\epsilon_W=2\gamma_1+2\lambda_1-\kappa_U-\kappa_W,}
```
```math
\boxed{\Sigma_Z=\delta\ln\!\left(\frac{Z_W}{\Omega_W^2}\right)=2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W.}
```
So the first grouped defect already depends only on five direct microscopic slippages.

---

## Step 11B — Exact tracking combination

The tracking-factor drift depends only on the single combination
```math
\boxed{
\Sigma_{\rm tr}:=(1+\chi_0)\Sigma_\delta+(1+\delta_U)\Sigma_\chi.
}
```
A direct logarithmic expansion of
```math
R_{\rm tr}=\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0}
```
gives
```math
\boxed{
\Theta_1
=
-\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}\,\Sigma_{\rm tr}.
}
```
So exact tracking rigidity is simply
```math
\Theta_1=0 \iff \Sigma_{\rm tr}=0
```
on the physical branch.

---

## Step 11C — Genuine nontracking transfer-shape slippage

The grouped defect still contains one genuinely nontracking scalar. Define
```math
\boxed{
\Sigma_{\rm nt}
:=
\Sigma_Z
+\frac{2\epsilon_W}{1-\epsilon}\frac{11+9\delta_U}{11(1+\delta_U)}\,\Sigma_\epsilon
-\left[
\frac{2\chi_0}{1+\delta_U}
+\frac{4\epsilon_W\delta_U}{11(1-\epsilon)(1+\delta_U)^2}
\right]\Sigma_\delta.
}
```
Then the coherent grouped defect separates exactly into tracking and nontracking pieces:
```math
\boxed{
\Xi_1
=
A_{\rm tr}\,\Sigma_{\rm tr}+\Sigma_{\rm nt},
\qquad
A_{\rm tr}:=\frac{2\chi_0}{(1+\chi_0)(1+\delta_U)}.
}
```
So once the universal tracking feed-through is subtracted off, the grouped defect is controlled by `\Sigma_{\rm nt}` alone.

---

## Step 11D — Dressing slippage and exact triangular normal form

The selected-branch relation
```math
R_{\rm target}\,\mathcal T^2=\Lambda_0(1-\epsilon_\eta)
```
with scalar front factor `\Lambda_0` inert at grouped weak-axisymmetric order gives
```math
\boxed{
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_\eta}{1-\epsilon_\eta}\,\Sigma_\eta,
\qquad
\Sigma_\eta=\delta\ln\epsilon_\eta.
}
```
Putting everything together, the reduced observable drifts take the exact triangular form
```math
\boxed{
\Theta_1=-C_{\rm tr}\Sigma_{\rm tr},
}
```
```math
\boxed{
\Xi_1=A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt},
}
```
```math
\boxed{
\mathcal R_1+\Xi_1=-\frac{\epsilon_\eta}{1-\epsilon_\eta}\Sigma_\eta,
}
```
with
```math
C_{\rm tr}
=
\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}.
```
This is the strongest exact compression so far.

The whole coherent weak-axisymmetric problem is no longer a five-slippage bookkeeping problem. It is a **three-scalar normal form**.

---

## Step 11E — Exact inverse reconstruction formulas

Because the normal form is triangular, it can be inverted exactly.

### Tracking slippage from `\Theta_1`
```math
\boxed{
\Sigma_{\rm tr}
=
-\frac{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}{\chi_0\delta_U}\,\Theta_1.
}
```

### Nontracking slippage from `(\Theta_1,\Xi_1)`
Using `A_{\rm tr}/C_{\rm tr}=2(1+\chi_0+\delta_U)/\delta_U`,
```math
\boxed{
\Sigma_{\rm nt}
=
\Xi_1+\frac{2(1+\chi_0+\delta_U)}{\delta_U}\,\Theta_1.
}
```
So `\Sigma_{\rm nt}` is the grouped defect with the universal tracking feed-through removed.

### Dressing slippage from `(\mathcal R_1,\Xi_1)`
```math
\boxed{
\Sigma_\eta
=
-\frac{1-\epsilon_\eta}{\epsilon_\eta}\,(\mathcal R_1+\Xi_1).
}
```
So the selected-branch residual directly measures the dressing mismatch after the direct transfer-shape defect is added back.

---

## Step 11F — Quartic anomaly gate in the branch-adapted coordinates

The carried quartic anomaly target remains
```math
\Lambda_1\approx 0.279605891931464.
```
So the exact coherent-branch quartic gate is now
```math
\boxed{
A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt}=\Lambda_1.
}
```
This has two immediate specializations.

### If exact tracking rigidity holds
```math
\Sigma_{\rm tr}=0
\quad\Longrightarrow\quad
\boxed{\Sigma_{\rm nt}=\Lambda_1.}
```
So on a tracking-rigid branch the missing `O(f^4)` anomaly layer is purely the **nontracking transfer-shape slippage**.

### If tracking rigidity and selected-branch rigidity both hold
If in addition `\mathcal R_1=0`, then
```math
\boxed{
\Sigma_\eta
=
-\frac{1-\epsilon_\eta}{\epsilon_\eta}\,\Lambda_1.
}
```
So a target-rigid selected branch would require a compensating dressing slippage of precisely that size.

---

## What this changes

This is a real simplification.

Before this step, the remaining g-2 problem still looked like a loose collection of microscopic placement drifts.
After this step, it is exact and minimal:

1. `\Sigma_{\rm tr}` — tracking slippage,
2. `\Sigma_{\rm nt}` — nontracking transfer-shape slippage,
3. `\Sigma_\eta` — dressing slippage.

That is the correct reduced coordinate system for the next moving-throat test.

---

## Continuation point

The next clean move is now very sharp:

> determine whether the actual moving-throat branch preserves or drives the three branch-adapted scalars `(\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta)`, and then push them one step further down to exact branch composites / microscopic monomials.

# Step 12 — Direct microscopic monomials and the exact compatibility ledger

## Goal

Step 11 compressed the coherent weak-axisymmetric problem to three exact
branch-adapted slippages,
```math
\Sigma_{\rm tr},\qquad \Sigma_{\rm nt},\qquad \Sigma_\eta.
```
The next clean move is to remove even that intermediate layer and identify the
**direct microscopic quantities** whose logarithmic drifts are those three scalars.

That is what this step does.

The main result is that the coherent branch is now controlled by three direct
microscopic monomials:
```math
\mathfrak C_{{\rm tr},*},\qquad \mathfrak C_{{\rm nt},*},\qquad \epsilon_\eta,
```
and the full zero-defect condition becomes an explicit three-equation
compatibility ledger for the microscopic grouped drifts
```math
(\lambda_1,c_1,\gamma_1,\kappa_U,\kappa_\eta,\kappa_W,\mu_1,\tau_1).
```

So the continuation point is now even smaller than Step 11 suggested:
it is no longer “track three abstract slippages,” but rather

> determine whether the actual moving-throat branch preserves three direct microscopic kernel monomials.

---

## Inputs carried forward

From Step 11, the coherent-kernel ratios are
```math
\chi_0=\frac{\gamma c_{\eta U}}{K_U},
\qquad
\delta_U=\frac{\pi^2T_U}{L^2K_U},
\qquad
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}},
```
```math
\epsilon_W=\frac{\gamma^2\lambda_W^2\sigma}{K_UK_W^{(\mathrm{eff})}},
\qquad
\frac{Z_W}{\Omega_W^2}
=
\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}.
```
The direct microscopic logarithmic drifts are
```math
\Sigma_\chi=\gamma_1+c_1-\kappa_U,
\qquad
\Sigma_\delta=\tau_1-\kappa_U,
```
```math
\Sigma_Z=2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W,
\qquad
\Sigma_\epsilon=2\gamma_1+2\lambda_1-\kappa_U-\kappa_W,
```
```math
\Sigma_\eta=2c_1-\kappa_U-\kappa_\eta.
```
And Step 11 already defined
```math
\Sigma_{\rm tr}
=
(1+\chi_{0,*})\Sigma_\delta+(1+\delta_{U,*})\Sigma_\chi,
```
```math
\Sigma_{\rm nt}
=
\Sigma_Z
+
E_*\Sigma_\epsilon
-
F_*\Sigma_\delta,
```
with
```math
\epsilon_* = \epsilon_{W,*}\!\left(1-\frac{2}{11}\frac{\delta_{U,*}}{1+\delta_{U,*}}\right).
```
And
```math
E_*=
\frac{2\epsilon_{W,*}}{1-\epsilon_*}\,
\frac{11+9\delta_{U,*}}{11(1+\delta_{U,*})},
```
```math
F_*=
\frac{2\chi_{0,*}}{1+\delta_{U,*}}
+
\frac{4\epsilon_{W,*}\delta_{U,*}}{11(1-\epsilon_*)(1+\delta_{U,*})^2}.
```

---

## Step 12A — Direct microscopic tracking monomial

The tracking scalar is already linear in the two logarithmic drifts
```math
\Sigma_\chi,\Sigma_\delta.
```
So freeze the reference-branch coefficients and define
```math
\boxed{
\mathfrak C_{{\rm tr},*}
:=
\chi_0^{\,1+\delta_{U,*}}
\delta_U^{\,1+\chi_{0,*}}.
}
```
Then
```math
\delta\ln\mathfrak C_{{\rm tr},*}
=
(1+\delta_{U,*})\,\delta\ln\chi_0
+
(1+\chi_{0,*})\,\delta\ln\delta_U
=
\Sigma_{\rm tr}.
```
So the tracking coordinate is not an abstract reduced object anymore:
```math
\boxed{
\delta\ln\mathfrak C_{{\rm tr},*}=\Sigma_{\rm tr}.
}
```

### Explicit microscopic form
```math
\boxed{
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{1+\chi_{0,*}}.
}
```

---

## Step 12B — Direct microscopic nontracking monomial

The genuine nontracking transfer-shape scalar is also logarithmic in direct
microscopic ratios. Define
```math
\boxed{
\mathfrak C_{{\rm nt},*}
:=
\frac{Z_W}{\Omega_W^2}\,
\epsilon_W^{E_*}\,
\delta_U^{-F_*}.
}
```
Then
```math
\delta\ln\mathfrak C_{{\rm nt},*}
=
\delta\ln\!\left(\frac{Z_W}{\Omega_W^2}\right)
+
E_*\,\delta\ln\epsilon_W
-
F_*\,\delta\ln\delta_U
=
\Sigma_{\rm nt}.
```
So
```math
\boxed{
\delta\ln\mathfrak C_{{\rm nt},*}=\Sigma_{\rm nt}.
}
```

### Explicit microscopic form
```math
\boxed{
\mathfrak C_{{\rm nt},*}
=
\frac{\lambda_W^2\mu_W}
{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}
\left(
\frac{\gamma^2\lambda_W^2\sigma}
{K_UK_W^{(\mathrm{eff})}}
\right)^{E_*}
\left(
\frac{\pi^2T_U}{L^2K_U}
\right)^{-F_*}.
}
```

---

## Step 12C — Dressing invariant

The third coordinate is already direct:
```math
\boxed{
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}},
\qquad
\delta\ln\epsilon_\eta=\Sigma_\eta.
}
```

So the full coherent weak-axisymmetric zero-defect theorem now reads
```math
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\ln\mathfrak C_{{\rm tr},*}=0,
\quad
\delta\ln\mathfrak C_{{\rm nt},*}=0,
\quad
\delta\ln\epsilon_\eta=0.
}
```

That is a real sharpening. The direct weak-axisymmetric branch equations are now
microscopic monomial equations.

---

## Step 12D — Exact monomial-drift matrix

Now collect the microscopic grouped drift vector as
```math
\delta\mathbf x
=
(\lambda_1,\ c_1,\ \gamma_1,\ \kappa_U,\ \kappa_\eta,\ \kappa_W,\ \mu_1,\ \tau_1)^T.
```
Then the direct monomial drifts satisfy
```math
\begin{pmatrix}
\delta\ln\mathfrak C_{{\rm tr},*}\\
\delta\ln\mathfrak C_{{\rm nt},*}\\
\delta\ln\epsilon_\eta
\end{pmatrix}
=
M_*\,
\delta\mathbf x,
```
with
```math
\boxed{
M_*=
\begin{pmatrix}
0 & 1+\delta_{U,*} & 1+\delta_{U,*} & -(2+\chi_{0,*}+\delta_{U,*}) & 0 & 0 & 0 & 1+\chi_{0,*}\\
2(1+E_*) & 0 & 2E_* & F_*-E_* & -1 & -(2+E_*) & 1 & -F_*\\
0 & 2 & 0 & -1 & -1 & 0 & 0 & 0
\end{pmatrix}.
}
```

The useful minor built from the columns
```math
(\tau_1,\kappa_\eta,\mu_1)
```
has determinant
```math
\boxed{
\det M_*^{(\tau_1,\kappa_\eta,\mu_1)}=1+\chi_{0,*}>0.
}
```
So on the physical branch
```math
\operatorname{rank}(M_*)=3,
\qquad
\dim\ker M_*=5.
```

That is the first exact sign that the monomial-rigid branch should be a
five-parameter similarity family rather than a fine-tuned isolated locus.

---

## Step 12E — Exact microscopic compatibility ledger

Setting the three monomial drifts to zero gives an explicit three-equation
compatibility system.

### Tracking compatibility
```math
\boxed{
(1+\chi_{0,*})(\tau_1-\kappa_U)
+
(1+\delta_{U,*})(\gamma_1+c_1-\kappa_U)
=0.
}
```

### Dressing compatibility
```math
\boxed{
2c_1-\kappa_U-\kappa_\eta=0.
}
```

### Nontracking compatibility
```math
\boxed{
2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W
+
E_*(2\gamma_1+2\lambda_1-\kappa_U-\kappa_W)
-
F_*(\tau_1-\kappa_U)
=0.
}
```

The script solves this system exactly for the three dependent drifts
```math
(\tau_1,\kappa_\eta,\mu_1).
```

### Solved form
```math
\boxed{
\tau_1
=
\kappa_U
-
\frac{1+\delta_{U,*}}{1+\chi_{0,*}}
(\gamma_1+c_1-\kappa_U),
}
```
```math
\boxed{
\kappa_\eta=2c_1-\kappa_U,
}
```
```math
\boxed{
\mu_1
=
2c_1-\kappa_U+2\kappa_W-2\lambda_1
-
E_*(2\gamma_1+2\lambda_1-\kappa_U-\kappa_W)
-
F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}}
(\gamma_1+c_1-\kappa_U).
}
```

So the zero-defect branch is now an explicit microscopic rigidity ledger rather
than a broad slippage statement.

---

## Step 12F — Quartic anomaly gate in the monomial coordinates

For the anomaly problem, Step 11 gave
```math
\Xi_1=A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt},
\qquad
A_{\rm tr}=\frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})}.
```
In the new direct coordinates that becomes
```math
\boxed{
\Xi_1
=
A_{\rm tr}\,\delta\ln\mathfrak C_{{\rm tr},*}
+
\delta\ln\mathfrak C_{{\rm nt},*}.
}
```
So the carried quartic anomaly target
```math
\Xi_1=\Lambda_1
```
is now
```math
\boxed{
A_{\rm tr}\,\delta\ln\mathfrak C_{{\rm tr},*}
+
\delta\ln\mathfrak C_{{\rm nt},*}
=
\Lambda_1.
}
```
with the carried numerical target
```math
\Lambda_1\approx 0.279605891931464.
```

### Tracking-rigid specialization
If the branch preserves the tracking monomial exactly,
```math
\delta\ln\mathfrak C_{{\rm tr},*}=0,
```
then the entire quartic anomaly target collapses to
```math
\boxed{
\delta\ln\mathfrak C_{{\rm nt},*}=\Lambda_1.
}
```

That is the cleanest direct microscopic statement of the missing anomaly layer so far.

---

## What this changes

This step is important because it removes the last layer of abstraction before
the similarity-orbit geometry.

Before Step 12, the continuation point was still phrased in terms of three exact
branch-adapted slippages.

After Step 12, the continuation point is:

1. one tracking monomial,
2. one nontracking transfer-shape monomial,
3. one dressing ratio,
4. and one explicit `3 x 8` microscopic compatibility matrix.

That is the right base for the exact similarity-orbit closure.

---

## Continuation point

The next clean move is now immediate:

> identify the exact five-parameter monomial-preserving similarity orbit whose tangent space is cut out by the Step-12 compatibility matrix \(M_*\).

That will turn the present linear rigidity ledger into a full geometric closure statement.
# Step 13 — Exact similarity orbit and quotient closure

## Goal

Step 12 identified the three **direct microscopic monomials**
```math
\mathfrak C_{{\rm tr},*},\qquad \mathfrak C_{{\rm nt},*},\qquad \epsilon_\eta,
```
and proved that their logarithmic drifts are exactly
```math
\Sigma_{\rm tr},\qquad \Sigma_{\rm nt},\qquad \Sigma_\eta.
```

The remaining question was then unavoidable:

> what is the exact microscopic redundancy that leaves those three monomials unchanged?

This step answers that question completely.

It proves that the coherent weak-axisymmetric zero-defect branch is the tangent
space of an exact **five-parameter multiplicative similarity orbit**, and that the
three monomials above are the exact **quotient coordinates**.

For the g-2 derivation, that means something very concrete:

- the anomaly problem is no longer an 8-variable microscopic drift problem,
- it is exactly a 3-coordinate quotient problem,
- and the direct quartic anomaly gate itself lives only in a 2-dimensional
  quotient plane.

---

## Inputs carried forward

From Step 12, define the direct microscopic state vector
```math
x=(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})},\ K_W^{(\mathrm{eff})},\ \mu_W,\ T_U),
```
with all entries positive on the coherent branch.

The three exact direct monomials are
```math
\boxed{
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{1+\chi_{0,*}},
}
```
```math
\boxed{
\mathfrak C_{{\rm nt},*}
=
\left(\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}\right)
\left(\frac{\gamma^2\lambda_W^2\sigma}{K_UK_W^{(\mathrm{eff})}}\right)^{E_*}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{-F_*},
}
```
```math
\boxed{
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}}.
}
```
And the exact monomial-drift matrix from Step 12 is
```math
\boxed{
M_*=
\begin{pmatrix}
0 & 1+\delta_{U,*} & 1+\delta_{U,*} & -(2+\chi_{0,*}+\delta_{U,*}) & 0 & 0 & 0 & 1+\chi_{0,*}\\
2(1+E_*) & 0 & 2E_* & F_*-E_* & -1 & -(2+E_*) & 1 & -F_*\\
0 & 2 & 0 & -1 & -1 & 0 & 0 & 0
\end{pmatrix}.
}
```

---

## Step 13A — The finite fibre equations are controlled by the same matrix

Take two positive microscopic states `x` and `\widetilde x`, and define their
finite logarithmic ratio vector
```math
\Delta\mathbf x=
\begin{pmatrix}
\Delta_\lambda,
\Delta_c,
\Delta_\gamma,
\Delta_U,
\Delta_\eta,
\Delta_W,
\Delta_\mu,
\Delta_T
\end{pmatrix}^T
=
\begin{pmatrix}
\ln(\widetilde\lambda_W/\lambda_W),
\ln(\widetilde c_{\eta U}/c_{\eta U}),
\ln(\widetilde\gamma/\gamma),
\ln(\widetilde K_U/K_U),
\ln(\widetilde K_\eta/K_\eta),
\ln(\widetilde K_W/K_W),
\ln(\widetilde\mu_W/\mu_W),
\ln(\widetilde T_U/T_U)
\end{pmatrix}^T.
```
Then the exact invariant log-ratios are
```math
q_{\rm tr}:=\ln\frac{\widetilde{\mathfrak C}_{{\rm tr},*}}{\mathfrak C_{{\rm tr},*}},
\qquad
q_{\rm nt}:=\ln\frac{\widetilde{\mathfrak C}_{{\rm nt},*}}{\mathfrak C_{{\rm nt},*}},
\qquad
q_\eta:=\ln\frac{\widetilde\epsilon_\eta}{\epsilon_\eta}.
```
Because the three invariants are monomials, their exact finite log-ratios are
still linear in `\Delta\mathbf x`:
```math
\boxed{
\begin{pmatrix}
q_{\rm tr}\\
q_{\rm nt}\\
q_\eta
\end{pmatrix}
=
M_*\,\Delta\mathbf x.
}
```

This is the first key result of the step.
The matrix that governed the infinitesimal compatibility ledger in Step 12 also
controls the **exact finite invariant fibres**.

---

## Step 13B — Exact five-parameter similarity orbit

Choose the five free finite logarithmic co-scalings
```math
(\Delta_\lambda,\Delta_c,\Delta_\gamma,\Delta_U,\Delta_W)
```
for
```math
(\lambda_W,c_{\eta U},\gamma,K_U,K_W^{(\mathrm{eff})}).
```
Setting the invariant log-ratios to zero,
```math
q_{\rm tr}=q_{\rm nt}=q_\eta=0,
```
gives the exact solved dependent shifts
```math
\boxed{
\Delta_\eta=2\Delta_c-\Delta_U,
}
```
```math
\boxed{
\Delta_T
=
\Delta_U
-
\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Delta_\gamma+\Delta_c-\Delta_U),
}
```
```math
\boxed{
\Delta_\mu
=
2\Delta_c-\Delta_U+2\Delta_W-2\Delta_\lambda
-
E_*(2\Delta_\gamma+2\Delta_\lambda-\Delta_U-\Delta_W)
-
F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Delta_\gamma+\Delta_c-\Delta_U).
}
```
Exponentiating gives the exact five-parameter multiplicative similarity orbit
`\mathcal G_*`.

So the zero-defect microscopic branch is not an isolated compatibility locus.
It is an exact codimension-3 similarity family.

---

## Step 13C — Tangent generator matrix

Writing the free orbit coordinates as
```math
g=(g_1,g_2,g_3,g_4,g_5)^T
```
for
```math
(\Delta_\lambda,\Delta_c,\Delta_\gamma,\Delta_U,\Delta_W),
```
the exact tangent generator matrix is
```math
\boxed{
G=
\begin{pmatrix}
1 & 0 & 0 & 0 & 0\\
0 & 1 & 0 & 0 & 0\\
0 & 0 & 1 & 0 & 0\\
0 & 0 & 0 & 1 & 0\\
0 & 2 & 0 & -1 & 0\\
0 & 0 & 0 & 0 & 1\\
-2(1+E_*) & 2-F_*\alpha & -2E_*-F_*\alpha & -1+E_*+F_*\alpha & 2+E_*\\
0 & -\alpha & -\alpha & 1+\alpha & 0
\end{pmatrix},
\qquad
\alpha=\frac{1+\delta_{U,*}}{1+\chi_{0,*}}.
}
```
The script verifies
```math
\boxed{M_*G=0.}
```
So the tangent space of the similarity orbit is exactly the kernel of the
monomial-drift map.

That is the infinitesimal version of the orbit theorem.

---

## Step 13D — Canonical exact quotient section

Now solve the general finite fibre equation
```math
M_*\Delta\mathbf x=q,
\qquad
q=(q_{\rm tr},q_{\rm nt},q_\eta)^T,
```
using the same five orbit coordinates as free variables.
The dependent coordinates become
```math
\boxed{
\Delta_\eta=2\Delta_c-\Delta_U-q_\eta,
}
```
```math
\boxed{
\Delta_T
=
\Delta_U
-
\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Delta_\gamma+\Delta_c-\Delta_U)
+
\frac{q_{\rm tr}}{1+\chi_{0,*}},
}
```
```math
\boxed{
\Delta_\mu
=
\Delta_\mu^{\rm orbit}
+
q_{\rm nt}-q_\eta+
\frac{F_*}{1+\chi_{0,*}}q_{\rm tr}.
}
```
So the exact finite microscopic state splits into
```math
\boxed{
\Delta\mathbf x = Gg + Sq,
}
```
where the chosen right-inverse section is
```math
\boxed{
S=
\begin{pmatrix}
0 & 0 & 0\\
0 & 0 & 0\\
0 & 0 & 0\\
0 & 0 & 0\\
0 & 0 & -1\\
0 & 0 & 0\\
\dfrac{F_*}{1+\chi_{0,*}} & 1 & -1\\
\dfrac{1}{1+\chi_{0,*}} & 0 & 0
\end{pmatrix},
\qquad
M_*S=I_3.
}
```

This is the exact quotient closure in the form most useful for g-2.
It says that the eight microscopic log-ratios split into

- **five similarity-orbit directions** `g`, and
- **three exact quotient coordinates** `q`.

---

## Step 13E — Canonical finite quotient representative

The chosen gauge slice keeps the five similarity coordinates fixed and moves only
along the three quotient directions. In finite multiplicative form this means
```math
\boxed{
K_\eta^{(\mathrm{eff})}\mapsto e^{-q_\eta}K_\eta^{(\mathrm{eff})},
}
```
```math
\boxed{
T_U\mapsto e^{q_{\rm tr}/(1+\chi_{0,*})}T_U,
}
```
```math
\boxed{
\mu_W
\mapsto
\exp\!\left(q_{\rm nt}-q_\eta+\frac{F_*}{1+\chi_{0,*}}q_{\rm tr}\right)\mu_W,
}
```
with the other five microscopic variables held fixed.

The script verifies directly that this representative satisfies
```math
\Delta\ln\mathfrak C_{{\rm tr},*}=q_{\rm tr},
\qquad
\Delta\ln\mathfrak C_{{\rm nt},*}=q_{\rm nt},
\qquad
\Delta\ln\epsilon_\eta=q_\eta.
```

So the quotient coordinates are not abstract. They have an explicit finite
microscopic representative.

---

## Step 13F — Observables in exact quotient coordinates

Now the coherent weak-axisymmetric observables become completely transparent.
Using the carried triangular normal form,
```math
\boxed{
\Theta_1
=
-
\frac{\chi_{0,*}\delta_{U,*}}
{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}
q_{\rm tr},
}
```
```math
\boxed{
\Xi_1
=
A_{\rm tr}\,q_{\rm tr}+q_{\rm nt},
\qquad
A_{\rm tr}=\frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})},
}
```
```math
\boxed{
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}q_\eta.
}
```

This is one of the most useful simplifications so far.
The direct quartic anomaly gate depends only on
```math
(q_{\rm tr},q_{\rm nt}),
```
while the third exact quotient coordinate `q_\eta` is orthogonal to the direct
`\Xi_1` match at this order.

So the g-2 problem has collapsed from

- 8 microscopic grouped drifts,

to

- a 3-dimensional exact quotient,

and then further, for the direct quartic match itself, to

- a 2-dimensional quotient plane.

---

## Step 13G — Canonical quartic anomaly representatives

The carried quartic anomaly target is
```math
\Xi_1=\Lambda_1.
```
In exact quotient coordinates this is simply
```math
\boxed{
A_{\rm tr}q_{\rm tr}+q_{\rm nt}=\Lambda_1.
}
```
So on the chosen quotient section the general canonical representative is
```math
\Delta\mathbf x
=
S\begin{pmatrix}
q_{\rm tr}\\
\Lambda_1-A_{\rm tr}q_{\rm tr}\\
q_\eta
\end{pmatrix}.
```
That is the exact gauge-fixed form of the whole quartic matching family.

### Simplest direct slice
Take the tracking-rigid and dressing-rigid slice
```math
q_{\rm tr}=0,
\qquad
q_\eta=0.
```
Then
```math
q_{\rm nt}=\Lambda_1,
```
and the canonical representative reduces to
```math
\boxed{
\Delta\mathbf x
=
(0,0,0,0,0,0,\Lambda_1,0)^T.
}
```
So in the chosen quotient gauge, the carried quartic anomaly target is
represented entirely by a **pure `\mu_W` drift**.

That does **not** mean the real moving-throat branch must physically vary only
`\mu_W`. It means that after modding out the exact five-parameter similarity
redundancy, every tracking-rigid, dressing-rigid quartic match can be represented
in that minimal way.

This is the cleanest microscopic simplification of the anomaly problem so far.

---

## What this changes

This step is the point where the moving-throat quotient geometry becomes directly
useful for g-2.

Before Step 13, the direct anomaly statement still lived inside an 8-variable
microscopic drift ledger.

After Step 13:

1. the exact microscopic redundancy is known,
2. the exact quotient coordinates are explicit,
3. the observables are linear in those quotient coordinates,
4. the direct quartic anomaly gate uses only two of the three,
5. and the chosen canonical slice gives an explicit minimal microscopic
   representative of the whole matching family.

So the next step no longer needs to “simplify the similarity-orbit algebra.”
That is done.
The next step should instead use the quotient closure to decide which exact
quotient trajectory the actual moving-throat branch is most naturally taking.

---

## Continuation point

The clean next move is now:

> use the exact quotient section to choose a physically meaningful branch condition — for example tracking-rigid, dressing-rigid, or minimum-norm in the quotient plane — and translate that into the simplest microscopic law for the missing quartic anomaly layer.

At this point the direct g-2 problem is genuinely small enough to test specific closures rather than reorganizing the algebra again.
# Step 14 — Exact branch-selection closures in the quotient

## Goal

Step 13 finished the exact quotient split:

- the direct quartic anomaly gate is
  ```math
  
  \Xi_1=A_{\rm tr}q_{\rm tr}+q_{\rm nt},
  ```
- the tracking drift is
  ```math
  \Theta_1=-C_{\rm tr}q_{\rm tr},
  ```
- and the selected-branch residual is
  ```math
  \mathcal R_1=-\Xi_1-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}q_\eta.
  ```

So the algebra itself is no longer the bottleneck.
The live question is now purely physical:

> which exact quotient trajectory should the missing `O(f^4)` common layer take?

This step answers that by writing the natural closure families explicitly and pushing each one through the canonical quotient section.

---

## Inputs carried forward

Define
```math
A_{\rm tr}:=\frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})},
\qquad
\alpha:=\frac{1}{1+\chi_{0,*}},
\qquad
\beta:=\frac{F_*}{1+\chi_{0,*}}.
```

The carried quartic anomaly target is
```math
\Xi_1=\Lambda_1.
```
So the exact matching family in quotient coordinates is
```math
\boxed{
q_{\rm nt}=\Lambda_1-A_{\rm tr}q_{\rm tr}.
}
```

On the canonical Step-13 quotient section, only three microscopic variables move:
```math
\boxed{
\Delta\ln K_\eta^{(\mathrm{eff})}=-q_\eta,
\qquad
\Delta\ln \mu_W=\beta q_{\rm tr}+q_{\rm nt}-q_\eta,
\qquad
\Delta\ln T_U=\alpha q_{\rm tr}.
}
```
So once a quotient closure is chosen, the microscopic law is immediate.

---

## Step 14A — Tracking-rigid + dressing-rigid closure

The simplest exact closure is
```math
q_{\rm tr}=0,
\qquad
q_\eta=0.
```
Then the quartic gate forces
```math
\boxed{q_{\rm nt}=\Lambda_1.}
```
So
```math
\Theta_1=0,
\qquad
\Xi_1=\Lambda_1,
\qquad
\mathcal R_1=-\Lambda_1.
```

The canonical microscopic representative collapses to
```math
\boxed{
\Delta\ln K_\eta^{(\mathrm{eff})}=0,
\qquad
\Delta\ln \mu_W=\Lambda_1,
\qquad
\Delta\ln T_U=0.
}
```
So in the chosen quotient gauge the whole quartic layer is a **pure `\mu_W` drift**.

In finite monomial form,
```math
\boxed{
\mathfrak C_{{\rm tr},*}\mapsto \mathfrak C_{{\rm tr},*},
\qquad
\mathfrak C_{{\rm nt},*}\mapsto e^{\Lambda_1}\mathfrak C_{{\rm nt},*},
\qquad
\epsilon_\eta\mapsto \epsilon_\eta.
}
```
With the carried Step-2 value
```math
\Lambda_1\approx 0.279605891931464,
```
this finite nontracking amplification is
```math
\boxed{e^{\Lambda_1}\approx 1.322608458944212.}
```

So the direct anomaly layer can already be read as:

> keep tracking and dressing fixed, and multiply the nontracking composite by about `1.3226`.

---

## Step 14B — Tracking-rigid + selected-branch coherent closure

If one also demands that the selected-branch residual vanish at the same order,
```math
\mathcal R_1=0,
```
while keeping the tracking channel rigid,
```math
q_{\rm tr}=0,
```
then
```math
\Xi_1=\Lambda_1
\quad\Longrightarrow\quad
q_{\rm nt}=\Lambda_1,
```
and the selected-branch condition fixes the dressing quotient coordinate uniquely:
```math
\boxed{
q_\eta=-\frac{1-\epsilon_{\eta,*}}{\epsilon_{\eta,*}}\Lambda_1.
}
```

So this branch has
```math
\Theta_1=0,
\qquad
\Xi_1=\Lambda_1,
\qquad
\mathcal R_1=0.
```

The canonical microscopic representative becomes
```math
\boxed{
\Delta\ln K_\eta^{(\mathrm{eff})}=\frac{1-\epsilon_{\eta,*}}{\epsilon_{\eta,*}}\Lambda_1,
\qquad
\Delta\ln \mu_W=\frac{\Lambda_1}{\epsilon_{\eta,*}},
\qquad
\Delta\ln T_U=0.
}
```
So the selected branch stays coherent only if the quartic nontracking drift is accompanied by a locked dressing co-drift.

In finite monomial form,
```math
\boxed{
\mathfrak C_{{\rm tr},*}\mapsto \mathfrak C_{{\rm tr},*},
\qquad
\mathfrak C_{{\rm nt},*}\mapsto e^{\Lambda_1}\mathfrak C_{{\rm nt},*},
\qquad
\epsilon_\eta\mapsto
\exp\!\left[-\frac{1-\epsilon_{\eta,*}}{\epsilon_{\eta,*}}\Lambda_1\right]\epsilon_\eta.
}
```

This is the first really useful physical closure of the quotient problem:

> the missing quartic layer is not an arbitrary three-variable drift; it is a one-parameter nontracking update, with dressing either held fixed or slaved to selected-branch coherence.

---

## Step 14C — Dressing-rigid minimum norm in the quotient plane

A different exact criterion is purely geometric:
keep the dressing rigid,
```math
q_\eta=0,
```
and minimize the quotient length
```math
q_{\rm tr}^2+q_{\rm nt}^2
```
subject to
```math
A_{\rm tr}q_{\rm tr}+q_{\rm nt}=\Lambda_1.
```

This gives
```math
\boxed{
q_{\rm tr}=\frac{A_{\rm tr}}{1+A_{\rm tr}^2}\Lambda_1,
\qquad
q_{\rm nt}=\frac{1}{1+A_{\rm tr}^2}\Lambda_1.
}
```
So this branch deliberately spreads the quartic load between tracking and nontracking motion.

The corresponding canonical microscopic law is
```math
\boxed{
\Delta\ln K_\eta^{(\mathrm{eff})}=0,
\qquad
\Delta\ln T_U=\frac{\alpha A_{\rm tr}}{1+A_{\rm tr}^2}\Lambda_1,
\qquad
\Delta\ln \mu_W=
\left(
\frac{1}{1+A_{\rm tr}^2}
+
\frac{\beta A_{\rm tr}}{1+A_{\rm tr}^2}
\right)\Lambda_1.
}
```

In finite monomial form,
```math
\boxed{
\mathfrak C_{{\rm tr},*}\mapsto
\exp\!\left(\frac{A_{\rm tr}}{1+A_{\rm tr}^2}\Lambda_1\right)
\mathfrak C_{{\rm tr},*},
\qquad
\mathfrak C_{{\rm nt},*}\mapsto
\exp\!\left(\frac{1}{1+A_{\rm tr}^2}\Lambda_1\right)
\mathfrak C_{{\rm nt},*}.
}
```

So the quotient-minimum branch is mathematically clean, but it achieves that by reopening the tracking invariant.

---

## Step 14D — Dressing-rigid minimum norm in the canonical microscopic section

The quotient metric is not the only possible notion of “smallest” deformation.
A more microscopic criterion is to minimize the actual canonical section norm,
```math
(\Delta\ln \mu_W)^2+(\Delta\ln T_U)^2,
```
with `q_eta=0` and `\Xi_1=\Lambda_1`.

This gives the exact branch
```math
\boxed{
q_{\rm tr}=\frac{A_{\rm tr}-\beta}{(A_{\rm tr}-\beta)^2+\alpha^2}\Lambda_1,
\qquad
q_{\rm nt}=\Lambda_1-A_{\rm tr}q_{\rm tr}.
}
```
So the microscopic optimum only collapses back to the pure nontracking law when
```math
\boxed{\beta=A_{\rm tr}}
\qquad\Longleftrightarrow\qquad
\boxed{F_* = \frac{2\chi_{0,*}}{1+\delta_{U,*}}.}
```

That is a useful diagnostic:

- if the true moving-throat branch happens to satisfy this relation, then the microscopic optimum is exactly the pure `\mu_W` closure;
- otherwise the microscopic optimum generically carries a small tracking admixture even when dressing is frozen.

---

## Main result of the step

The quotient problem is now genuinely explicit.
The full quartic matching family is
```math
q_{\rm nt}=\Lambda_1-A_{\rm tr}q_{\rm tr},
```
but the physically meaningful exact closures collapse to a very short list:

### 1. Direct quartic closure
```math
q_{\rm tr}=0,
\qquad q_\eta=0,
\qquad q_{\rm nt}=\Lambda_1.
```
This is the **pure nontracking monomial law**.

### 2. Selected-branch coherent quartic closure
```math
q_{\rm tr}=0,
\qquad q_{\rm nt}=\Lambda_1,
\qquad q_\eta=-\frac{1-\epsilon_{\eta,*}}{\epsilon_{\eta,*}}\Lambda_1.
```
This is the **locked nontracking+dressing law**.

### 3. Purely geometric minimum-norm closures
These exist, but they do so by reintroducing tracking drift.

So the algebraic situation is now much cleaner than at the start of the anomaly repair:

> to improve the quartic layer, we do not need a vague “common charge-inertia transport PDE” anymore.
> We need the actual branch to tell us whether the first omitted common layer is
> - a pure nontracking update,
> - or a nontracking update with a slaved dressing correction.

That is a much smaller theorem target.

---

## Continuation point

The next clean move is now:

> translate these exact quotient closures back into the actual moving-throat branch composites
> ```math
> R_{\rm tr},\qquad \mathfrak N_*,\qquad \epsilon_\eta,
> ```
> and determine which of the two tracking-rigid laws the real branch is taking.

At that point the missing quartic anomaly layer will be expressed directly as a finite law for the exact branch composites rather than as an abstract quotient displacement.
# Step 15 — Branch-composite laws and the universal transfer-shape drift

## Goal

Step 14 solved the quartic anomaly layer in the **exact quotient coordinates**
```math
(q_{\rm tr},q_{\rm nt},q_\eta),
```
but the continuation point it left was still physical rather than algebraic:

> how do those quotient closures look when translated back into the **actual moving-throat branch composites**
> ```math
> R_{\rm tr},\qquad \mathfrak N_*,\qquad \epsilon_\eta,
> ```
> and, more importantly, what part of the quartic update is already fixed before the branch decides how to move?

This step answers that.

The main result is unexpectedly strong:

```math
\boxed{\delta\ln \mathcal T^2 = \Lambda_1}
```

at the carried first omitted common order, **for the full Step-14 matching family**, independent of the tracking choice `q_tr` and independent of the dressing choice `q_eta`.

So the real branch-selection question is no longer “what is the quartic transfer-shape update?”
That part is already fixed.
The only remaining question is whether that universal transfer-shape drift is accompanied by

- no dressing response, or
- the locked dressing co-drift that keeps the selected branch coherent.

---

## Inputs carried forward

From Step 14, the quartic matching family is
```math
q_{\rm nt}=\Lambda_1-A_{\rm tr}q_{\rm tr},
```
with
```math
A_{\rm tr}=\frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})}.
```

From the moving-throat branch-composite dictionary,
```math
\delta\ln R_{\rm tr}=-\frac{1}{C_*}q_{\rm tr},
\qquad
\delta\ln \mathfrak N_*=q_{\rm nt},
\qquad
\delta\ln \epsilon_\eta=q_\eta,
```
where
```math
\mathfrak N_*:=\mathcal T^2 R_{\rm tr}^{B_*},
```
with
```math
B_*:=\frac{2(1+\chi_{0,*}+\delta_{U,*})}{\delta_{U,*}},
\qquad
C_*:=\frac{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}{\chi_{0,*}\delta_{U,*}}.
```

A crucial identity is
```math
\boxed{A_{\rm tr}=\frac{B_*}{C_*}.}
```
The script verifies this exactly.

---

## Step 15A — The direct transfer shape is universal

Because
```math
\mathfrak N_* = \mathcal T^2 R_{\rm tr}^{B_*},
```
we have at the carried weak-axisymmetric order
```math
\delta\ln \mathcal T^2
=
\delta\ln \mathfrak N_* - B_*\,\delta\ln R_{\rm tr}.
```
Substitute the matching family and the branch-composite drift laws:
```math
\delta\ln \mathcal T^2
=
(\Lambda_1-A_{\rm tr}q_{\rm tr})
-
B_*\left(-\frac{q_{\rm tr}}{C_*}\right).
```
Using `A_tr = B_*/C_*`, the `q_tr` dependence cancels identically and one gets
```math
\boxed{\delta\ln \mathcal T^2 = \Lambda_1.}
```

So the quartic anomaly gate fixes the direct transfer shape completely.
Tracking admixture can move `R_tr` and therefore reshuffle `\mathfrak N_*`, but it cannot change the direct `\mathcal T^2` update.

This is the sharpest simplification reached so far.

---

## Step 15B — The selected-branch composite isolates the real remaining choice

Define the selected-branch dressing composite
```math
\mathfrak E := 1-\epsilon_\eta.
```
At the carried order,
```math
\delta\ln \mathfrak E
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}q_\eta,
```
and therefore
```math
\mathcal R_1 = \delta\ln \mathfrak E - \Lambda_1.
```
So the residual after the universal `\mathcal T^2` update is controlled entirely by the dressing response.

That means the physical branch question has collapsed to:

> does the branch let `\mathfrak E = 1-\epsilon_\eta` stay rigid while `\mathcal T^2` increases, or does it co-drift so that the selected branch remains coherent?

---

## Step 15C — The two tracking-rigid laws in actual branch variables

### 1. Direct tracking-rigid closure

Step 14’s simplest closure is
```math
q_{\rm tr}=0,
\qquad
q_\eta=0,
\qquad
q_{\rm nt}=\Lambda_1.
```
Translated back to the actual branch composites, this gives
```math
\boxed{\delta\ln R_{\rm tr}=0,}
```
```math
\boxed{\delta\ln \mathfrak N_*=\Lambda_1,}
```
```math
\boxed{\delta\ln \mathcal T^2=\Lambda_1,}
```
```math
\boxed{\delta\ln \epsilon_\eta=0,}
```
```math
\boxed{\delta\ln(1-\epsilon_\eta)=0,}
```
and therefore
```math
\boxed{\mathcal R_1=-\Lambda_1.}
```

So this branch says:

- keep the tracking factor fixed,
- keep the dressing ratio fixed,
- but increase the direct transfer shape by the universal quartic amount.

The selected branch then lags by exactly `-\Lambda_1`.

### 2. Tracking-rigid + selected-branch coherent closure

Step 14’s coherent tracking-rigid closure is
```math
q_{\rm tr}=0,
\qquad
q_{\rm nt}=\Lambda_1,
\qquad
q_\eta=-\frac{1-\epsilon_{\eta,*}}{\epsilon_{\eta,*}}\Lambda_1.
```
In actual branch variables this becomes
```math
\boxed{\delta\ln R_{\rm tr}=0,}
```
```math
\boxed{\delta\ln \mathfrak N_*=\Lambda_1,}
```
```math
\boxed{\delta\ln \mathcal T^2=\Lambda_1,}
```
```math
\boxed{\delta\ln \epsilon_\eta=-\frac{1-\epsilon_{\eta,*}}{\epsilon_{\eta,*}}\Lambda_1,}
```
```math
\boxed{\delta\ln(1-\epsilon_\eta)=\Lambda_1,}
```
and therefore
```math
\boxed{\mathcal R_1=0.}
```

So this branch says:

- keep the tracking factor fixed,
- impose the same universal transfer-shape update,
- and add the one dressing co-drift needed to keep the selected branch coherent.

### The crucial comparison

The two tracking-rigid closures have the **same** direct law
```math
\delta\ln \mathcal T^2=\Lambda_1.
```
They differ only in whether `\epsilon_\eta` / `(1-\epsilon_\eta)` responds.

That is the real reduction achieved by this step.

---

## Step 15D — Tracking admixture only repartitions `R_tr` and `\mathfrak N_*`

For the full dressing-rigid family,
```math
q_\eta=0,
\qquad
q_{\rm nt}=\Lambda_1-A_{\rm tr}q_{\rm tr},
```
we get
```math
\delta\ln R_{\rm tr}=-\frac{q_{\rm tr}}{C_*},
```
```math
\delta\ln \mathfrak N_*=
\Lambda_1-A_{\rm tr}q_{\rm tr},
```
but still
```math
\boxed{\delta\ln \mathcal T^2=\Lambda_1.}
```

So any nonzero `q_tr` only changes the partition between

- the tracking factor `R_tr`, and
- the corrected nontracking composite `\mathfrak N_*`,

while the direct transfer shape remains fixed.

This applies to both of Step 14’s minimum-norm closures:

- the quotient-plane minimum,
- and the minimum in the canonical microscopic section.

Neither changes the direct `\mathcal T^2` update.

---

## Main result of the step

After translating the Step-14 quotient closures back into the actual moving-throat branch composites, the quartic anomaly layer splits into:

### Universal part
```math
\boxed{\delta\ln \mathcal T^2 = \Lambda_1.}
```
This is fixed as soon as `\Xi_1=\Lambda_1` is imposed.

### Optional tracking part
```math
\boxed{\delta\ln R_{\rm tr}=-\frac{1}{C_*}q_{\rm tr}.}
```
This only repartitions the update between `R_tr` and `\mathfrak N_*`.

### Optional dressing part
```math
\boxed{\delta\ln \epsilon_\eta=q_\eta,}
```
which equivalently controls
```math
\delta\ln(1-\epsilon_\eta)
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}q_\eta.
```
This is the only part that decides whether the selected branch stays coherent.

So the continuation point is now much smaller than Step 14 left it:

> the real moving-throat branch does **not** need to tell us the direct transfer-shape quartic update anymore.
> That part is already fixed.
> It only needs to tell us whether the universal `\mathcal T^2` drift comes with
> - no dressing response, or
> - the locked dressing co-drift required for selected-branch coherence.

---

## Continuation point

The next clean move is now:

> insert the actual coherent local tracking-branch formulas for
> ```math
> R_{\rm tr}(\chi_0,\delta_U),
> \qquad
> 1-\epsilon_\eta,
> \qquad
> \mathcal T^2,
> ```
> and test whether the moving-throat branch leaves `(1-\epsilon_\eta)` rigid or makes it co-move with the universal transfer-shape update.

At that point the branch-selection problem will have collapsed to a single dressing-response test.
# Step 16 — Selected-branch demand-ratio law

## Goal

Step 15 translated the Step-14 quotient closures into the actual moving-throat branch composites and proved the strongest new fact so far:

```math
\boxed{\delta\ln \mathcal T^2 = \Lambda_1}
```

at the carried first omitted common order, for the **entire** quartic matching family.

So the direct transfer-shape update is no longer the live ambiguity.
The natural next question is now smaller:

> how does the same quartic layer look from the **selected-branch** side?
> Does it change the selected-branch demand ratio `R_{\rm target}`, or does it leave that ratio fixed and instead force the dressing sector to co-move?

This step answers that exactly.

---

## Inputs carried forward

The exact selected-branch identity is
```math
R_{\rm target}\,\mathcal T^2 = \Lambda_0(1-\epsilon_\eta),
```
with `\Lambda_0` inert at the carried grouped weak-axisymmetric order.

From Step 15,
```math
\delta\ln \mathcal T^2 = \Lambda_1.
```
And from the branch-composite dressing law,
```math
\delta\ln(1-\epsilon_\eta)
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}q_\eta.
```

---

## Step 16A — General demand-ratio law

Take the logarithmic drift of
```math
R_{\rm target}\,\mathcal T^2 = \Lambda_0(1-\epsilon_\eta).
```
Since `\Lambda_0` is fixed at this order,
```math
\delta\ln R_{\rm target} + \delta\ln \mathcal T^2
=
\delta\ln(1-\epsilon_\eta).
```
Substitute the Step-15 universal transfer-shape law:
```math
\delta\ln R_{\rm target}
=
\delta\ln(1-\epsilon_\eta)-\Lambda_1.
```
So the carried exact dressing-coordinate form is
```math
\boxed{
\delta\ln R_{\rm target}
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}q_\eta - \Lambda_1.
}
```

This already collapses the remaining branch-selection question to a single scalar.

---

## Step 16B — The selected-branch residual is literally the demand-ratio drift

From Step 15,
```math
\mathcal R_1 = \delta\ln(1-\epsilon_\eta)-\Lambda_1.
```
Comparing with the formula above gives
```math
\boxed{\mathcal R_1 = \delta\ln R_{\rm target}.}
```

That is the main interpretive result of the step.

The selected-branch residual is not an extra abstract defect anymore.
It is exactly the logarithmic drift of the selected-branch demand ratio.

So the branch-selection problem can now be asked as:

> does the omitted common quartic layer retarget `R_{\rm target}`, or does it preserve `R_{\rm target}` and instead force the dressing sector to absorb the universal `\mathcal T^2` update?

---

## Step 16C — The two tracking-rigid laws become a yes/no test on `R_target`

### 1. Direct tracking-rigid closure

The direct tracking-rigid branch has
```math
q_\eta = 0.
```
So
```math
\delta\ln(1-\epsilon_\eta)=0,
```
and therefore
```math
\boxed{\delta\ln R_{\rm target}=-\Lambda_1.}
```

This means:

- the direct transfer shape still increases by `\Lambda_1`,
- the dressing composite stays rigid,
- and the selected-branch demand ratio itself must shift downward by the same amount.

Equivalently,
```math
\boxed{\mathcal R_1=-\Lambda_1.}
```

### 2. Tracking-rigid + selected-branch coherent closure

The coherent tracking-rigid branch has
```math
q_\eta=-\frac{1-\epsilon_{\eta,*}}{\epsilon_{\eta,*}}\Lambda_1.
```
Then
```math
\delta\ln(1-\epsilon_\eta)=\Lambda_1,
```
so
```math
\boxed{\delta\ln R_{\rm target}=0.}
```

This means:

- the direct transfer shape still increases by `\Lambda_1`,
- the dressing composite co-moves by exactly the same logarithmic amount,
- and the selected-branch demand ratio stays fixed.

Equivalently,
```math
\boxed{\mathcal R_1=0.}
```

---

## Step 16D — The real branch-selection problem is now one scalar criterion

After Step 16, the two tracking-rigid laws are separated by a single condition:

### Direct law
```math
\boxed{\delta\ln R_{\rm target}=-\Lambda_1.}
```
The omitted common layer is absorbed by retargeting the selected-branch demand ratio.

### Coherent law
```math
\boxed{\delta\ln R_{\rm target}=0.}
```
The selected-branch demand ratio stays fixed, and the dressing sector co-moves to preserve coherence.

So the quartic branch-selection problem has become:

> does the actual moving-throat branch treat `R_{\rm target}` as fixed spectral data, or does it let the first omitted common layer shift that target itself?

That is much sharper than the original “common charge-inertia transport” bottleneck.

---

## Main result of the step

Step 15 fixed the universal direct law
```math
\delta\ln \mathcal T^2 = \Lambda_1.
```
Step 16 now fixes the selected-branch interpretation of the same update:

```math
\boxed{\mathcal R_1 = \delta\ln R_{\rm target}.}
```

So the remaining ambiguity in the quartic anomaly layer is no longer about the transfer shape at all.
It is only about whether the selected-branch demand ratio drifts.

That is the smallest clean formulation of the branch-selection question reached so far.

---

## Continuation point

The next clean move is now:

> insert the actual coherent local tracking-branch formulas for
> ```math
> R_{\rm target},\qquad R_{\rm tr},\qquad \epsilon_\eta,
> ```
> from the moving-throat kernel notes and test whether the physical branch naturally prefers
> ```math
> \delta\ln R_{\rm target}=0
> ```
> or
> ```math
> \delta\ln R_{\rm target}=-\Lambda_1.
> ```

At that point the quartic anomaly repair will have been reduced to one direct demand-ratio check on the real branch.
# Step 17 — Coherent-kernel microscopic demand-ratio ledger

## Goal

Step 16 reduced the quartic branch-selection question to one scalar criterion:

```math

delta\ln R_{\rm target}=0
\qquad\text{or}\qquad
\delta\ln R_{\rm target}=-\Lambda_1.
```

But that was still an abstract branch-composite statement.
The natural next move is to insert the **actual coherent local D/N kernel formulas** from the moving-throat notes and ask what the microscopic variables really control.

This step does that.

The main result is stronger than expected:

```math
\boxed{
\Lambda_1
=
-\delta\ln\Lambda
-2\,\delta\ln(1-\epsilon)
+\delta\ln Z_W
+2\,\delta\ln(1+\chi_0)
}
```

so the quartic universal transfer-shape drift is controlled entirely by the
**mixed/outgoing microscopic variables**.
The wall–`U` dressing coordinate `\epsilon_\eta` cancels identically from the inferred `\Lambda_1` law.

At the same time, the coherent support ratio `\zeta` does not enter either
`R_{\rm tr}` or `R_{\rm target}` at all. It only moves the total baseline `M_{\rm tr}`.

So the next branch-selection step is no longer “what does the support lane do to the target ratio?”
It does nothing to the target ratio directly.
The live question is whether the actual support lane compensates the exact tracking deficit by increasing the baseline at fixed target, or whether the completed PDE forces a genuine retargeting drift in the mixed/outgoing microdata.

---

## Inputs from the coherent local D/N kernel

From the coherent-kernel placement map in the moving-throat notes,

```math
R_{\rm tr}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0},
```

```math
\epsilon
=\epsilon_W^{(\mathrm{split})}
=\epsilon_W\Bigl[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\Bigr],
```

```math
M_{\rm mix}
=
\frac{8 Z_W (1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\epsilon)},
```

```math
M_{\rm supp}
=
\frac{8 \zeta Z_W (1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\zeta\epsilon)},
```

```math
M_{\rm tr}=M_{\rm mix}+M_{\rm supp}=M_{\rm mix}S(\zeta;\epsilon),
```
with
```math
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
```
and
```math
R_{\rm target}
=
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2}.
```

These are exactly the formulas that Stage 30 compresses the coherent local kernel into.

---

## Step 17A — The support ratio drops out of both `R_tr` and `R_target`

The script verifies

```math
\frac{\partial R_{\rm tr}}{\partial \zeta}=0,
\qquad
\frac{\partial R_{\rm target}}{\partial \zeta}=0.
```

So on the actual coherent local-kernel branch:

- `\zeta` does **not** move the tracking factor,
- `\zeta` does **not** move the demand ratio,
- `\zeta` only moves the total baseline through
  ```math
  M_{\rm tr}=M_{\rm mix}S(\zeta;\epsilon).
  ```

That already tells us something important:
whatever the support lane is doing, it is **not** directly retargeting the selected-branch demand ratio.

---

## Step 17B — Exact microscopic drift of the tracking factor

Use the natural logarithmic drift variables

```math
q_\chi:=\delta\ln(1+\chi_0),
\qquad
q_U:=\delta\ln(1+\delta_U).
```

Then the exact logarithmic drift of the tracking factor is

```math
\boxed{
\delta\ln R_{\rm tr}
=
-\frac{\delta_U\,q_\chi+\chi_0\,q_U}{1+\chi_0+\delta_U}.
}
```

So the constructive branch lowers `R_{\rm tr}` only through motion in the two microscopic placement variables `(1+\chi_0)` and `(1+\delta_U)`.
The support ratio `\zeta` is absent.

---

## Step 17C — Exact microscopic drift of the split blocking ratio

If

```math
q_W:=\delta\ln\epsilon_W,
```
then the split blocking ratio obeys

```math
\boxed{
q_\epsilon:=\delta\ln\epsilon
=
q_W-\frac{2}{11+9\delta_U}\,q_U.
}
```

So the split-blocking drift is fixed by the primitive mixed blocking drift plus the axial split drift.

---

## Step 17D — Exact microscopic drift of the demand ratio

Using the natural drifts

```math
q_\Lambda:=\delta\ln\Lambda,
\qquad
q_Z:=\delta\ln Z_W,
\qquad
q_\eta:=\delta\ln\epsilon_\eta,
\qquad
q_\epsilon:=\delta\ln\epsilon,
```

we get the exact microscopic demand-ratio law

```math
\boxed{
\delta\ln R_{\rm target}
=
q_\Lambda-q_Z-2q_\chi
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}q_\eta
-2\frac{\epsilon_*}{1-\epsilon_*}q_\epsilon.
}
```

Equivalently,

```math
\delta\ln R_{\rm target}
=
\delta\ln(1-\epsilon_\eta)
+\delta\ln\Lambda
+2\,\delta\ln(1-\epsilon)
-\delta\ln Z_W
-2\,\delta\ln(1+\chi_0).
```

So the actual coherent-kernel demand ratio is controlled by

- the common wall–`U` dressing factor `(1-\epsilon_\eta)`,
- the mixed/outgoing scale `\Lambda`,
- the split blocking `(1-\epsilon)`,
- the wall-to-mixed overlap `Z_W`,
- and the common interference ratio `(1+\chi_0)`.

Again, `\zeta` is absent.

---

## Step 17E — Exact cancellation theorem with the Step-16 law

Step 16 proved abstractly that

```math
\delta\ln R_{\rm target}
=
\delta\ln(1-\epsilon_\eta)-\Lambda_1.
```

Now insert the exact microscopic Stage-30 formula above. The `\delta\ln(1-\epsilon_\eta)` term cancels identically, leaving

```math
\boxed{
\Lambda_1
=
-\delta\ln\Lambda
-2\,\delta\ln(1-\epsilon)
+\delta\ln Z_W
+2\,\delta\ln(1+\chi_0).
}
```

Or, in the logarithmic drift variables,

```math
\boxed{
\Lambda_1
=
-q_\Lambda+q_Z+2q_\chi+2\frac{\epsilon_*}{1-\epsilon_*}q_\epsilon.
}
```

This is the main result of the step.

The quartic universal transfer-shape drift is **not** controlled by

- the support ratio `\zeta`, or
- the wall–`U` dressing coordinate `\epsilon_\eta`.

It is carried entirely by the mixed/outgoing microscopic variables.

---

## Step 17F — What the cancellation means physically

This exact cancellation changes the interpretation of the remaining anomaly problem.

The branch-selection ambiguity is no longer “support versus dressing.”
The support lane does not move `R_{\rm target}` at all, and the wall–`U` dressing coordinate cancels from the inferred `\Lambda_1` law.

So the next honest question is now much sharper:

> does the physical coherent support lane compensate the exact tracking-branch deficit by increasing the available baseline at fixed demand ratio, or does the completed PDE force a real retargeting drift in the mixed/outgoing microdata `(\Lambda,Z_W,\chi_0,\epsilon)`?

That is the direct continuation point.

---

## Continuation point

The next clean move is to use the exact support-compensation theorem from the coherent local D/N branch and see whether the physical support lane is naturally a **fixed-target / load-compensation** mechanism.

If it is, then the real PDE-side continuation will favor the Step-16 coherent side

```math
\delta\ln R_{\rm target}=0,
```

rather than the direct retargeting law.
# Step 18 — Support-compensation selection law

## Goal

Step 17 inserted the actual coherent local D/N kernel formulas and showed two exact facts:

1. the support ratio `\zeta` does **not** enter `R_{\rm tr}` or `R_{\rm target}`,
2. the quartic universal transfer-shape drift `\Lambda_1` is fixed entirely by the mixed/outgoing microscopic variables.

So the next clean question is no longer abstract:

> when the real coherent support lane turns on, does it behave like direct retargeting, or does it act by increasing the available baseline at fixed selected-branch demand ratio?

This step answers that using the exact support-compensation theorem from the moving-throat notes.

The main result is:

```math
\boxed{
\text{the coherent local D/N support lane is structurally a fixed-target / load-compensation mechanism.}
}
```

That does **not** prove that the whole quartic anomaly correction is literally nothing but `\zeta`-motion.
But it does show that the natural PDE-side continuation of the physical coherent branch lies on the **coherent** side of Step 16,

```math
\delta\ln R_{\rm target}=0,
```

not on the direct-retargeting side.

---

## Inputs from the coherent local D/N branch

The support-enhancement factor is

```math
M_{\rm tr}=M_{\rm mix}S(\zeta;\epsilon),
```
with
```math
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
```
for
```math
0<\epsilon<1,
\qquad
0\le \zeta<1/\epsilon.
```

The physical tracking branch obeys

```math
M_{\rm tr}=G_{\rm tr}(\xi,\delta;R_{\rm tr}),
```
where

```math
G_{\rm tr}(\xi,\delta;R)
=
\frac{9\xi(\xi+\delta)}{9\delta+(9+2R^2)\xi}.
```

Its exact critical load is

```math
M_{\rm crit}(\delta,R)
=
G_{\rm tr}(1,\delta;R)
=
\frac{9(1+\delta)}{9\delta+9+2R^2}.
```

And the normalization target is fixed by the selected-branch law

```math
R_{\rm target}=F_{\rm tr}(\xi,\delta;R_{\rm tr}).
```

---

## Step 18A — Exact inverse of the support-enhancement factor

The support factor is strictly increasing:

```math
\frac{dS}{d\zeta}=
\frac{1-\epsilon}{(1-\zeta\epsilon)^2}>0.
```

It obeys

```math
S(0;\epsilon)=1,
\qquad
\lim_{\zeta\to 1/\epsilon^-} S(\zeta;\epsilon)=+\infty.
```

So every finite required enhancement `S_{\rm req}>1` has a unique support ratio

```math
\boxed{
\zeta_{\rm req}
=
\frac{S_{\rm req}-1}{1+\epsilon(S_{\rm req}-2)}.
}
```

The exact stability margin below the blocking pole is

```math
\boxed{
\frac{1}{\epsilon}-\zeta_{\rm req}
=
\frac{1-\epsilon}{\epsilon[1+\epsilon(S_{\rm req}-2)]}>0.
}
```

So any finite required support enhancement sits strictly inside the stable window.

---

## Step 18B — Exact support-compensation theorem

For every finite target ratio `R_{\rm target}>1`, there exists a stable-side branch point

```math
\xi_{\rm req}\in(0,1)
```

with

```math
F_{\rm tr}(\xi_{\rm req},\delta;R_{\rm tr})=R_{\rm target}.
```

Define the corresponding required load

```math
M_{\rm req}:=G_{\rm tr}(\xi_{\rm req},\delta;R_{\rm tr}).
```

Because `G_{\rm tr}` is strictly increasing on the stable branch,

```math
0<M_{\rm req}<M_{\rm crit}(\delta,R_{\rm tr}).
```

If the mixed-only coherent branch is already strong enough,

```math
M_{\rm mix}\ge M_{\rm req},
```

then the target is reached with

```math
\zeta_{\rm req}=0.
```

If instead

```math
M_{\rm mix}<M_{\rm req},
```

then the exact required enhancement is

```math
S_{\rm req}=\frac{M_{\rm req}}{M_{\rm mix}}>1,
```

and the unique coherent support ratio that hits the target is

```math
\boxed{
\zeta_{\rm req}
=
\frac{S_{\rm req}-1}{1+\epsilon(S_{\rm req}-2)}.
}
```

Moreover,

```math
S_{\rm req}<S_{\rm crit}:=\frac{M_{\rm crit}}{M_{\rm mix}},
```
so
```math
\boxed{
\zeta_{\rm req}<\zeta_{\rm crit}<1/\epsilon.
}
```

The exact gap is

```math
\zeta_{\rm crit}-\zeta_{\rm req}
=
\frac{(S_{\rm crit}-S_{\rm req})(1-\epsilon)}{[1+\epsilon(S_{\rm crit}-2)][1+\epsilon(S_{\rm req}-2)]}>0.
```

So the target is reached **before** the selected branch softens out.

That is the exact reduced support-feasibility theorem.

---

## Step 18C — Support moves the branch deeper into the same tracking family

Combining the coherent support law and the tracking branch gives

```math
M_{\rm mix}S(\zeta;\epsilon)=G_{\rm tr}(\xi_{\rm phys},\delta;R_{\rm tr}).
```

Differentiate implicitly:

```math
\boxed{
\frac{d\xi_{\rm phys}}{d\zeta}
=
\frac{M_{\rm mix}\,dS/d\zeta}{(dG_{\rm tr}/d\xi)_{\xi_{\rm phys}}}>0.
}
```

So coherent support enhancement always drives the physical branch to larger softening depth.

This is the exact reduced statement behind the compensation theorem:

- support does **not** alter `R_{\rm target}`,
- it only increases the available baseline,
- and that increase moves the physical branch monotonically deeper into the same tracking family.

---

## Step 18D — Branch-selection implication for the quartic anomaly layer

Now compare with Step 16.

Step 16 said the two tracking-rigid options are separated by

```math
\delta\ln R_{\rm target}=0
\qquad\text{vs.}\qquad
\delta\ln R_{\rm target}=-\Lambda_1.
```

Step 17 already proved that `\zeta` does not enter `R_{\rm target}` at all.
Step 18 now proves that `\zeta` only increases `M_{\rm tr}` and drives the branch deeper into the same tracking family.

So the coherent local D/N support lane is structurally a

```math
\boxed{\text{fixed-target / load-compensation mechanism}.}
```

That means the natural PDE-side continuation of the physical branch favors the
**coherent** side of Step 16,

```math
\boxed{\delta\ln R_{\rm target}=0,}
```

rather than the direct retargeting law.

This is not yet a proof that the whole quartic anomaly correction is exhausted by support enhancement alone.
Any genuine retargeting would still have to come from the mixed/outgoing microscopic variables isolated in Step 17.

But it is the first exact reason, taken directly from the coherent local-kernel map, to say that the actual branch prefers the coherent interpretation.

---

## Main result of the step

The old branch-selection ambiguity has now become asymmetric.

- The mixed/outgoing microscopic variables control the universal transfer-shape drift `\Lambda_1`.
- The coherent support lane does not retarget `R_{\rm target}` at all.
- Instead it compensates the exact tracking deficit by increasing `M_{\rm tr}` and moving the physical branch deeper into the tracking family.

So the completed PDE does not need to invent a new branch-selection principle here.
It already has one: the support lane acts at fixed demand ratio.

---

## Continuation point

The next clean move is now very specific:

> combine the Step-17 mixed/outgoing microledger with the Step-18 fixed-target support-compensation law and solve for the minimal mixed/outgoing drift pattern that supplies the observed `\Lambda_1`, while the coherent support lane provides the needed baseline compensation.

That is the first point where the quartic g-2 repair can be turned into a concrete microscopic balance law instead of a branch-choice ambiguity.
# Step 19 — Coherent mixed/outgoing balance law

## Goal

Step 17 proved that the quartic transport-shape correction is carried entirely by
mixed/outgoing microscopic variables, while Step 18 proved that the coherent
support lane works at **fixed** demand ratio and only raises the available
baseline.

So the next honest question is now very sharp:

> on the coherent branch, what are the exact mixed/outgoing balance laws that
> carry the observed `\Lambda_1`, and how do they split from the separate
> baseline-compensation job of the support factor `S(\zeta;\epsilon)`?

This step answers that.

The main result is that the quartic repair splits into an exact two-equation
microledger:

```math
\boxed{
\Lambda_1
=
-q_\Lambda + q_Z + 2 q_\chi + \sigma q_\epsilon,
\qquad
\Delta_{\rm mix}
=
q_\Lambda - \frac{\sigma}{2} q_\epsilon,
}
```
with
```math
\sigma := \frac{2\epsilon_*}{1-\epsilon_*},
\qquad
\Delta_{\rm mix}:=\delta\ln M_{\rm mix}.
```

So the coherent branch does **not** leave one vague “missing quartic mechanism.”
It leaves two sharply separated mixed/outgoing tasks:

1. one linear combination produces the universal transfer-shape correction
   `\Lambda_1`,
2. another linear combination sets the mixed-only baseline drift
   `\Delta_{\rm mix}`,
3. and the support lane then carries any remaining baseline compensation through
   `M_{\rm tr}=M_{\rm mix}S(\zeta;\epsilon)`.

---

## Inputs from the coherent local D/N kernel

From the coherent-kernel placement map,

```math
M_{\rm mix}
=
\frac{8 Z_W (1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\epsilon)},
```

```math
R_{\rm target}
=
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
```

so the exact mixed-only product law is

```math
\boxed{
R_{\rm target} M_{\rm mix} = \frac{8\Lambda(1-\epsilon)}{\pi^2}.
}
```

That is the coherent-kernel version of the old product law: the mixed lane sets
one bare scale, and the support lane later multiplies it by the enhancement
factor `S(\zeta;\epsilon)`.

---

## Step 19A — Exact drift laws for `R_target`, `M_mix`, and their product

Use the logarithmic drift variables

```math
q_\Lambda:=\delta\ln\Lambda,
\qquad
q_Z:=\delta\ln Z_W,
\qquad
q_\chi:=\delta\ln(1+\chi_0),
\qquad
q_\eta:=\delta\ln\epsilon_\eta,
\qquad
q_\epsilon:=\delta\ln\epsilon.
```

Then the exact coherent-kernel drifts are

```math
\delta\ln R_{\rm target}
=
q_\Lambda-q_Z-2q_\chi
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}q_\eta
-2\frac{\epsilon_*}{1-\epsilon_*}q_\epsilon,
```

```math
\delta\ln M_{\rm mix}
=
q_Z+2q_\chi
+\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}q_\eta
+\frac{\epsilon_*}{1-\epsilon_*}q_\epsilon.
```

Adding them gives the exact mixed-only product drift

```math
\boxed{
\delta\ln(R_{\rm target}M_{\rm mix})
=
q_\Lambda - \frac{\epsilon_*}{1-\epsilon_*}q_\epsilon.
}
```

So the coherent branch already shows something useful:

- `q_Z` and `q_\chi` redistribute the product between `R_target` and `M_mix`,
- but the bare product scale itself is controlled only by `q_\Lambda` and
  `q_\epsilon`.

---

## Step 19B — Coherent fixed-target closure and the slaved dressing law

Step 18 picked the coherent branch-selection side

```math
\delta\ln R_{\rm target}=0,
```

and Step 14 / Step 16 already showed that the coherent tracking-rigid branch
slaves the wall–`U` dressing drift to the quartic correction:

```math
\boxed{
q_\eta
=
-\frac{1-\epsilon_{\eta,*}}{\epsilon_{\eta,*}}\,\Lambda_1.
}
```

Insert that into the exact `\delta\ln R_{\rm target}` law. The result is

```math
0
=
\Lambda_1
+q_\Lambda-q_Z-2q_\chi
-2\frac{\epsilon_*}{1-\epsilon_*}q_\epsilon.
```

Equivalently,

```math
\boxed{
\Lambda_1
=
-q_\Lambda + q_Z + 2 q_\chi + \sigma q_\epsilon,
\qquad
\sigma:=\frac{2\epsilon_*}{1-\epsilon_*}.
}
```

This is the exact coherent fixed-target transport-shape law.

---

## Step 19C — The mixed-only baseline drift

Now compute the coherent mixed-baseline drift by inserting the same slaved
`q_\eta` into `\delta\ln M_{\rm mix}` and then eliminating `\Lambda_1` with the
fixed-target relation from Step 19B. One gets

```math
\boxed{
\Delta_{\rm mix}
:=
\delta\ln M_{\rm mix}
=
q_\Lambda - \frac{\sigma}{2} q_\epsilon.
}
```

That is the second half of the balance law.

So on the coherent branch:

- the pair `(q_\Lambda,q_\epsilon)` controls the mixed-only product/baseline,
- the pair `(q_Z,q_\chi)` redistributes scale between target and baseline,
- and the wall–`U` dressing drift is no longer free but slaved to `\Lambda_1`.

---

## Step 19D — Exact solution for `(q_Λ, q_ε)`

Treating `q_Z`, `q_\chi`, `\Lambda_1`, and `\Delta_{\rm mix}` as the visible
independent data, the coherent balance pair solves exactly to

```math
\boxed{
q_\Lambda
=
\Lambda_1 + 2\Delta_{\rm mix} - q_Z - 2 q_\chi,
}
```

```math
\boxed{
q_\epsilon
=
\frac{1-\epsilon_*}{\epsilon_*}
\bigl(\Lambda_1 + \Delta_{\rm mix} - q_Z - 2 q_\chi\bigr).
}
```

This is the cleanest exact parameterization of the coherent mixed/outgoing
quartic layer so far.

A few useful special sub-branches follow immediately.

### Support-carried baseline closure

If the mixed/outgoing lane is chosen to supply **only** the transfer-shape
repair, while the support factor carries all remaining baseline compensation,
then set

```math
\Delta_{\rm mix}=0.
```

The exact family reduces to

```math
q_\Lambda = \Lambda_1 - q_Z - 2 q_\chi,
```

```math
q_\epsilon = \frac{1-\epsilon_*}{\epsilon_*}(\Lambda_1 - q_Z - 2 q_\chi).
```

### Frozen overlap/interference branch

If instead one freezes the wall-overlap and common-interference drifts,

```math
q_Z=q_\chi=0,
```

then

```math
q_\Lambda = \Lambda_1 + 2\Delta_{\rm mix},
```

```math
q_\epsilon = \frac{1-\epsilon_*}{\epsilon_*}(\Lambda_1 + \Delta_{\rm mix}).
```

So in that branch the quartic layer lives entirely in the outgoing scale and
split-blocking drifts.

---

## Step 19E — Minimum-norm reduced representative

Now treat the reduced mixed/outgoing variable set as

```math
x=(q_\Lambda,q_Z,q_\chi,q_\epsilon)^T,
```

with two exact constraints

```math
(-1,1,2,\sigma)\cdot x = \Lambda_1,
```

```math
(1,0,0,-\sigma/2)\cdot x = \Delta_{\rm mix}.
```

The minimum-Euclidean-norm solution is

```math
\boxed{
x_{\min}=A^T(AA^T)^{-1}b,
}
```
with

```math
A=
\begin{pmatrix}
-1 & 1 & 2 & \sigma \\
 1 & 0 & 0 & -\sigma/2
\end{pmatrix},
\qquad
b=
\begin{pmatrix}
\Lambda_1 \\
\Delta_{\rm mix}
\end{pmatrix}.
```

Explicitly,

```math
q_\Lambda^{\min}
=
\frac{2\Delta_{\rm mix}(\sigma^2+10)+\Lambda_1\sigma^2}{2(3\sigma^2+10)},
```

```math
q_Z^{\min}
=
\frac{2\Delta_{\rm mix}(\sigma^2+2)+\Lambda_1(\sigma^2+4)}{2(3\sigma^2+10)},
```

```math
q_\chi^{\min}
=
\frac{2\Delta_{\rm mix}(\sigma^2+2)+\Lambda_1(\sigma^2+4)}{3\sigma^2+10},
```

```math
q_\epsilon^{\min}
=
\frac{\sigma(\Lambda_1-4\Delta_{\rm mix})}{3\sigma^2+10}.
```

On the support-carried-baseline closure `\Delta_{\rm mix}=0`, this becomes

```math
q_\Lambda^{\min}
=
\frac{\Lambda_1\sigma^2}{2(3\sigma^2+10)},
```

```math
q_Z^{\min}
=
\frac{\Lambda_1(\sigma^2+4)}{2(3\sigma^2+10)},
qquad
q_\chi^{\min}
=
\frac{\Lambda_1(\sigma^2+4)}{3\sigma^2+10},
```

```math
q_\epsilon^{\min}
=
\frac{\Lambda_1\sigma}{3\sigma^2+10}.
```

So the smallest reduced coherent repair distributes the quartic correction across
all four mixed/outgoing coordinates rather than placing it in one slot.

---

## Main result of the step

The quartic branch problem is now much sharper than before.

The coherent branch is organized by an exact **mixed/outgoing balance pair**:

```math
\Lambda_1
=
-q_\Lambda + q_Z + 2 q_\chi + \sigma q_\epsilon,
```

```math
\Delta_{\rm mix}
=
q_\Lambda - \frac{\sigma}{2}q_\epsilon.
```

That means:

- `\Lambda_1` is the universal transfer-shape repair,
- `\Delta_{\rm mix}` is the mixed-only baseline drift,
- the support factor `S(\zeta;\epsilon)` handles whatever further baseline
  compensation is needed,
- and the wall–`U` dressing drift is no longer a free ambiguity but a slaved
  coherent response.

So after Step 19, the anomaly repair is no longer “choose a branch and hope.”
It is an exact microscopic balance law.

---

## Continuation point

The next clean move is now obvious:

> resolve the split-blocking drift
> ```math
> q_\epsilon
> = q_W - \frac{2}{11+9\delta_{U,*}} q_U
> ```
> back into the primitive coherent-kernel variables and derive the minimum-norm
> **primitive** microdrift pattern.

That will tell us whether the quartic repair is dominated by wall blocking,
axial splitting, overlap, interference, or outgoing-scale motion at the actual
microscopic ledger level.
# Step 20 — Split-`U` primitive resolution of the quartic layer

## Goal

Step 19 reduced the coherent quartic repair to an exact two-equation law in the
reduced mixed/outgoing variables

```math
(q_\Lambda,q_Z,q_\chi,q_\epsilon).
```

But `q_\epsilon` is not itself primitive. Step 17 already showed that the
split-blocking drift resolves as

```math
q_\epsilon
=
q_W - \frac{2}{11+9\delta_{U,*}}q_U.
```

So the next clean move is obvious:

> resolve the coherent quartic balance law all the way back to the primitive
> microscopic variables and identify the smallest primitive drift pattern that
> realizes the observed `\Lambda_1`.

That is exactly what this step does.

The main result is simple and useful:

```math
\boxed{
q_U^{\min} = -\beta\,q_W^{\min},
\qquad
\beta:=\frac{2}{11+9\delta_{U,*}}<\frac{2}{11}.
}
```

So on the minimum-norm primitive coherent closure, the split-`U` motion is
always a **suppressed counter-drift** of the wall-blocking drift. The quartic
repair is therefore dominated by `q_W`, not by large axial split motion.

---

## Inputs from Step 17 and Step 19

From Step 17C,

```math
q_\epsilon
=
q_W - \beta q_U,
\qquad
\beta:=\frac{2}{11+9\delta_{U,*}}.
```

From Step 19,

```math
\Lambda_1
=
-q_\Lambda + q_Z + 2 q_\chi + \sigma q_\epsilon,
```

```math
\Delta_{\rm mix}
=
q_\Lambda - \frac{\sigma}{2} q_\epsilon,
```

with

```math
\sigma:=\frac{2\epsilon_*}{1-\epsilon_*}>0.
```

Substituting the split-`U` resolution of `q_\epsilon` gives the exact primitive
coherent law.

---

## Step 20A — Exact primitive coherent balance law

The quartic transfer-shape correction becomes

```math
\boxed{
\Lambda_1
=
-q_\Lambda + q_Z + 2 q_\chi + \sigma q_W - \sigma\beta q_U,
}
```

and the mixed-only baseline drift becomes

```math
\boxed{
\Delta_{\rm mix}
=
q_\Lambda - \frac{\sigma}{2}q_W + \frac{\sigma\beta}{2}q_U.
}
```

So the primitive coherent microledger has five visible variables,

```math
(q_\Lambda,q_Z,q_\chi,q_W,q_U),
```

but only two exact balance constraints at this order.

---

## Step 20B — Exact solution for `(q_Λ, q_W)`

Solving the primitive balance pair for `(q_\Lambda,q_W)` gives

```math
\boxed{
q_\Lambda
=
\Lambda_1 + 2\Delta_{\rm mix} - q_Z - 2 q_\chi,
}
```

```math
\boxed{
q_W
=
\beta q_U + \frac{2}{\sigma}
\bigl(\Lambda_1 + \Delta_{\rm mix} - q_Z - 2 q_\chi\bigr).
}
```

So once `(q_Z,q_\chi,q_U)` are chosen, the coherent primitive repair is fully
fixed.

A few useful special branches follow immediately.

### Pure wall-blocking realization

If the split-`U` drift is set to zero,

```math
q_U=0,
```

then

```math
q_W
=
\frac{2}{\sigma}
\bigl(\Lambda_1 + \Delta_{\rm mix} - q_Z - 2 q_\chi\bigr).
```

### Pure split-`U` realization

If one instead forces `q_W=0`, then the entire split-blocking correction must be
carried by axial splitting:

```math
q_U
=
-\frac{2}{\beta\sigma}
\bigl(\Lambda_1 + \Delta_{\rm mix} - q_Z - 2 q_\chi\bigr).
```

Since `\beta<2/11`, this is the more expensive primitive realization. The same
net `q_\epsilon` requires a much larger `q_U` than `q_W`.

---

## Step 20C — Minimum-norm primitive representative

Now define the primitive variable vector

```math
x=(q_\Lambda,q_Z,q_\chi,q_W,q_U)^T,
```

and the two exact coherent constraints

```math
(-1,1,2,\sigma,-\sigma\beta)\cdot x = \Lambda_1,
```

```math
(1,0,0,-\sigma/2,\sigma\beta/2)\cdot x = \Delta_{\rm mix}.
```

The minimum-Euclidean-norm primitive solution is again the pseudoinverse
representative

```math
\boxed{
x_{\min}=A^T(AA^T)^{-1}b,
}
```
with

```math
A=
\begin{pmatrix}
-1 & 1 & 2 & \sigma & -\sigma\beta \\
 1 & 0 & 0 & -\sigma/2 & \sigma\beta/2
\end{pmatrix},
\qquad
b=
\begin{pmatrix}
\Lambda_1 \\
\Delta_{\rm mix}
\end{pmatrix}.
```

Explicitly,

```math
q_\Lambda^{\min}
=
\frac{2\Delta_{\rm mix}(\beta^2\sigma^2+\sigma^2+10)+\Lambda_1\sigma^2(1+\beta^2)}{2\,[3\sigma^2(1+\beta^2)+10]},
```

```math
q_Z^{\min}
=
\frac{2\Delta_{\rm mix}(\beta^2\sigma^2+\sigma^2+2)+\Lambda_1(\beta^2\sigma^2+\sigma^2+4)}{2\,[3\sigma^2(1+\beta^2)+10]},
```

```math
q_\chi^{\min}
=
\frac{2\Delta_{\rm mix}(\beta^2\sigma^2+\sigma^2+2)+\Lambda_1(\beta^2\sigma^2+\sigma^2+4)}{3\sigma^2(1+\beta^2)+10},
```

```math
q_W^{\min}
=
\frac{\sigma(\Lambda_1-4\Delta_{\rm mix})}{3\sigma^2(1+\beta^2)+10},
```

```math
q_U^{\min}
=
\frac{\beta\sigma(4\Delta_{\rm mix}-\Lambda_1)}{3\sigma^2(1+\beta^2)+10}.
```

Now the key cancellation appears immediately:

```math
\boxed{
q_U^{\min} = -\beta\,q_W^{\min}.
}
```

So the minimum-norm primitive solution automatically chooses the smallest split-`U`
excursion compatible with the needed reduced `q_\epsilon`.

On the support-carried-baseline closure `\Delta_{\rm mix}=0`, this simplifies to

```math
q_W^{\min}
=
\frac{\Lambda_1\sigma}{3\sigma^2(1+\beta^2)+10},
```

```math
q_U^{\min}
=
-\frac{\Lambda_1\beta\sigma}{3\sigma^2(1+\beta^2)+10}.
```

---

## Step 20D — Primitive dominance consequence

Because the constructive coherent branch has `\delta_{U,*}>0`, we have the exact
bound

```math
0<\beta=rac{2}{11+9\delta_{U,*}}<\frac{2}{11}.
```

Combine that with the minimum-norm identity `q_U^{\min}=-\beta q_W^{\min}`.
Then

```math
\boxed{
|q_U^{\min}| = \beta |q_W^{\min}| < \frac{2}{11}|q_W^{\min}|.
}
```

So the primitive quartic repair is always dominated by the wall-blocking drift
`q_W`. The split-`U` drift is a smaller opposite-sign companion.

This is the first primitive-microscopic hierarchy statement in the g-2 rebuild.
It says the coherent quartic repair naturally prefers

- wall blocking / mixed support placement,
- overlap / interference drift,
- outgoing-scale motion,

rather than a large axial-splitting correction.

---

## Main result of the step

After resolving `q_\epsilon` into primitive coherent-kernel variables, the
quartic repair is no longer abstract at all.

The primitive coherent balance law is

```math
\Lambda_1
=
-q_\Lambda + q_Z + 2 q_\chi + \sigma q_W - \sigma\beta q_U,
```

```math
\Delta_{\rm mix}
=
q_\Lambda - \frac{\sigma}{2}q_W + \frac{\sigma\beta}{2}q_U,
```

and the minimum-norm primitive realization satisfies the exact suppression law

```math
q_U^{\min} = -\beta q_W^{\min},
qquad
0<\beta<\frac{2}{11}.
```

So the primitive quartic repair is structurally a **wall-blocking-dominated**
coherent mixed/outgoing correction, with only a smaller split-`U`
counter-drift.

---

## Continuation point

The next clean move is now to insert actual coherent-branch values for

```math
\epsilon_*,\qquad \delta_{U,*},
```

so that

```math
\sigma = \frac{2\epsilon_*}{1-\epsilon_*},
\qquad
\beta = \frac{2}{11+9\delta_{U,*}},
```

become numerical, and then compare the sizes of the primitive weights
multiplying

```math
q_\Lambda,\ q_Z,\ q_\chi,\ q_W,\ q_U.
```

That will tell us which microscopic drift channel the actual coherent branch most
naturally uses to close the final `10^{-12}` anomaly sliver.
# Step 21 — Coherent-branch value audit and primitive weight ranking

## Goal

At the end of Step 20 the next obvious move was:

> insert the actual coherent-branch values of `\epsilon_*` and `\delta_{U,*}` and
> rank the primitive quartic weights multiplying
> `(q_\Lambda,q_Z,q_\chi,q_W,q_U)`.

After checking the newer moving-throat notes carefully, the honest situation is:

- the coherent branch **does** fix the exact functional dependence of the quartic
  repair on `\epsilon_*` and `\delta_{U,*}`,
- but the current notes still leave those two branch numbers **parametric** rather
  than pinning them to one unique numerical pair.

So Step 21 does the right thing for the current file stack:

1. it performs the exact branch-value audit,
2. derives the primitive weight formulas on the whole constructive branch,
3. proves the exact ordering identities,
4. and splits the allowed branch into the three ranking regimes that remain.

That is already enough to answer the physically important question:

> **which primitive microscopic drifts are generically the dominant carriers of the
> missing quartic anomaly layer?**

The answer is sharp:

```math
\boxed{q_\chi \text{ is always largest, } q_Z \text{ is always second, and } q_U
\text{ is never dominant.}}
```

The only nontrivial crossover left is whether `q_\Lambda` sits above or below `q_W`.

---

## What the coherent notes actually fix

The coherent branch gives the exact reduced parameters

```math
\sigma = \frac{2\epsilon_*}{1-\epsilon_*},
\qquad
\beta = \frac{2}{11+9\delta_{U,*}},
```

with constructive-branch inequalities

```math
0 < \epsilon_* < 1,
\qquad
\delta_{U,*} > 0,
\qquad
0 < \beta < \frac{2}{11}.
```

So the branch is already sharply constrained, but not yet numerically frozen.

---

## Input from Step 20

On the support-carried-baseline closure `\Delta_{\rm mix}=0`, the Step-20 minimum-norm
primitive representative has the form

```math
q_\Lambda = \Lambda_1 w_\Lambda,
\qquad
q_Z = \Lambda_1 w_Z,
\qquad
q_\chi = \Lambda_1 w_\chi,
\qquad
q_W = \Lambda_1 w_W,
\qquad
q_U = -\Lambda_1 |w_U|.
```

If one keeps `(\sigma,\beta)` as the primitive branch variables, the weights are

```math
w_\Lambda
=
\frac{\sigma^2(1+\beta^2)}{2\,[3\sigma^2(1+\beta^2)+10]},
```

```math
w_Z
=
\frac{\sigma^2(1+\beta^2)+4}{2\,[3\sigma^2(1+\beta^2)+10]},
```

```math
w_\chi
=
\frac{\sigma^2(1+\beta^2)+4}{3\sigma^2(1+\beta^2)+10},
```

```math
w_W
=
\frac{\sigma}{3\sigma^2(1+\beta^2)+10},
```

```math
|w_U|
=
\frac{\beta\sigma}{3\sigma^2(1+\beta^2)+10}.
```

But once the coherent branch relation

```math
\sigma = \frac{2\epsilon_*}{1-\epsilon_*}
```

is substituted, these collapse to a much simpler direct rational form.

---

## Step 21A — Exact primitive weights in `(\epsilon_*,\beta)`

Define the positive denominator

```math
N(\epsilon_*,\beta)
:=
5 - 10\epsilon_* + (11+6\beta^2)\epsilon_*^2.
```

Equivalently,

```math
N(\epsilon_*,\beta)
=
5(1-\epsilon_*)^2 + 6\epsilon_*^2(1+\beta^2) > 0.
```

Then the primitive weights are

```math
\boxed{
 w_\Lambda
 =
 \frac{\epsilon_*^2(1+\beta^2)}{N}
}
```

```math
\boxed{
 w_Z
 =
 \frac{1-2\epsilon_*+(2+\beta^2)\epsilon_*^2}{N}
}
```

```math
\boxed{
 w_\chi
 =
 \frac{2\bigl[1-2\epsilon_*+(2+\beta^2)\epsilon_*^2\bigr]}{N}
}
```

```math
\boxed{
 w_W
 =
 \frac{\epsilon_*(1-\epsilon_*)}{N}
}
```

```math
\boxed{
 |w_U|
 =
 \frac{\beta\epsilon_*(1-\epsilon_*)}{N}
}
```

with

```math
q_U = -\Lambda_1 |w_U|.
```

So Step 20’s primitive closure can now be read directly as a function of the two
actual coherent branch parameters.

---

## Step 21B — Exact identities and always-true orderings

The first exact identity is

```math
\boxed{w_\chi = 2w_Z.}
```

So `q_\chi` is always exactly twice the `q_Z` weight on the support-carried
minimum-norm closure.

Next,

```math
w_Z - w_\Lambda
=
\frac{(1-\epsilon_*)^2}{N} > 0.
```

So `q_Z` is always larger than `q_\Lambda`.

Also,

```math
w_Z - w_W
=
\frac{\beta^2\epsilon_*^2 + 3(\epsilon_* - 1/2)^2 + 1/4}{N} > 0,
```

so `q_Z` is always larger than `q_W`.

And

```math
w_W - |w_U|
=
\frac{\epsilon_*(1-\epsilon_*)(1-\beta)}{N} > 0,
```

because `0 < \beta < 2/11 < 1`. Hence `q_W` is always larger than `|q_U|`.

Finally,

```math
w_\chi - w_W
=
\frac{2\beta^2\epsilon_*^2 + 5(\epsilon_* - 1/2)^2 + 3/4}{N} > 0,
```

so `q_\chi` is always the largest primitive weight.

Putting these together gives the branch-independent part of the ranking:

```math
\boxed{
q_\chi > q_Z > q_W > |q_U|,
\qquad
q_Z > q_\Lambda.
}
```

So the only unresolved comparison is:

```math
q_W \ \,\text{vs.}\ \, q_\Lambda.
```

---

## Step 21C — The two exact threshold surfaces

### Threshold 1: `q_W` versus `q_\Lambda`

The crossover is

```math
w_W = w_\Lambda
\iff
\epsilon_* = \frac{1}{2+\beta^2}.
```

So:

- if
  ```math
  \epsilon_* < \frac{1}{2+\beta^2},
  ```
  then `q_W > q_\Lambda`;
- if
  ```math
  \epsilon_* > \frac{1}{2+\beta^2},
  ```
  then `q_\Lambda > q_W`.

Because `0<\beta<2/11`, this crossover always sits extremely close to half-blocking:

```math
\frac{121}{246} < \frac{1}{2+\beta^2} < \frac12.
```

So `q_\Lambda` can overtake `q_W` only when `\epsilon_*` is already quite near `1/2`.

### Threshold 2: `|q_U|` versus `q_\Lambda`

The second crossover is

```math
|w_U| = w_\Lambda
\iff
\epsilon_* = \frac{\beta}{1+\beta+\beta^2}.
```

So:

- if
  ```math
  \epsilon_* < \frac{\beta}{1+\beta+\beta^2},
  ```
  then `|q_U| > q_\Lambda`;
- if
  ```math
  \epsilon_* > \frac{\beta}{1+\beta+\beta^2},
  ```
  then `q_\Lambda > |q_U|`.

Since `0<\beta<2/11`, this threshold is always small:

```math
0 < \frac{\beta}{1+\beta+\beta^2} < \frac{22}{147} \approx 0.14966.
```

So `|q_U|` can beat `q_\Lambda` only in the very weak-blocking corner of the
constructive branch.

---

## Step 21D — The three exact ranking regimes

These two thresholds divide the whole constructive branch into three clean regions.

### Region I — very weak blocking

If

```math
0 < \epsilon_* < \frac{\beta}{1+\beta+\beta^2},
```

then

```math
\boxed{q_\chi > q_Z > q_W > |q_U| > q_\Lambda.}
```

### Region II — intermediate blocking

If

```math
\frac{\beta}{1+\beta+\beta^2}
< \epsilon_* <
\frac{1}{2+\beta^2},
```

then

```math
\boxed{q_\chi > q_Z > q_W > q_\Lambda > |q_U|.}
```

### Region III — strong blocking

If

```math
\epsilon_* > \frac{1}{2+\beta^2},
```

then

```math
\boxed{q_\chi > q_Z > q_\Lambda > q_W > |q_U|.}
```

So the branch picture is fully explicit now.

---

## Step 21E — Limiting behavior

### Weak-blocking limit `\epsilon_* \to 0^+`

One finds

```math
w_\Lambda \to 0,
\qquad
w_Z \to \frac15,
\qquad
w_\chi \to \frac25,
\qquad
w_W \to 0,
\qquad
|w_U| \to 0.
```

So in the extreme weak-blocking corner, the quartic sliver is dominated by the
interference / overlap pair `(q_\chi,q_Z)`, with the split-blocking lanes becoming
subleading.

### Strong-blocking limit `\epsilon_* \to 1^-`

One gets

```math
w_\Lambda \to \frac16,
\qquad
w_Z \to \frac16,
\qquad
w_\chi \to \frac13,
\qquad
w_W \to 0,
\qquad
|w_U| \to 0.
```

So at the strong-blocking end, the quartic repair collapses to the purely
interference/outgoing trio

```math
q_\chi : q_Z : q_\Lambda = 2 : 1 : 1,
```

while the split-blocking drifts turn off.

---

## Physical reading

This step is the cleanest answer the current coherent branch can honestly support.

The quartic anomaly repair is **not** generically driven by large axial split
motion.
Instead:

1. `q_\chi` is always the dominant primitive carrier;
2. `q_Z` is always the second-largest one;
3. `q_U` is always subleading to `q_W`;
4. `q_\Lambda` only overtakes `q_W` when `\epsilon_*` is already near the upper,
   strong-blocking part of the constructive branch.

So the missing quartic layer is generically an

```math
\text{interference / outgoing-scale / wall-blocking correction,}
```

not a large `U`-splitting correction.

That is already a strong improvement over the older staggered picture.

---

## What the next step should be

The honest next move is now sharper than before:

> either derive the actual coherent branch point `(\epsilon_*,\delta_{U,*})` from the
> moving-throat PDE / branch-selection side,
> or derive an independent normalization / stability law that pins one of the two
> coordinates.

Until that is done, Step 21 is the correct stopping point:

- the current notes do **not** yet give a unique numerical pair,
- but they do already give the **exact phase diagram** of primitive quartic weight
  dominance.
# Step 22 — Support-demand selector and the exact `\epsilon_*` window

## Goal

Step 21 left one real microscopic ambiguity in the quartic layer:

```math
q_W \text{ versus } q_\Lambda.
```

The threshold is exact:

```math
q_W = q_\Lambda
\iff
\epsilon_* = \frac{1}{2+\beta^2},
\qquad
\beta = \frac{2}{11+9\delta_{U,*}}.
```

So the obvious next question is no longer “what are the primitive weights?”
It is:

> **what does the moving-throat support-selection side allow the coherent branch
> to do in `\epsilon_*` at all?**

That question is now answerable, because the later moving-throat support notes
already sharpen the coherent branch to an exact support-demand law:

- the coherent support enhancement is
  ```math
  S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
  ```
- the mixed-only product scale is
  ```math
  C_{\rm mix}=\frac{8\Lambda(1-\epsilon_*)}{\pi^2},
  ```
- and the branch demand is carried by the selected-branch product `\Pi_{\rm tr}`. 

So Step 22 turns the support theorem into an exact selector on `\epsilon_*`, and
then intersects that selector with the Step-21 `q_W/q_\Lambda` crossover.

The main result is:

```math
\boxed{
\varrho := \frac{\pi^2\Pi_{\rm tr}}{16\Lambda}
}
```

collapses the whole support-selection problem to one scalar demand parameter, and

```math
\boxed{
\begin{aligned}
&\text{mixed-only enough} &&\iff&& 0<\epsilon_* \le 1-2\varrho,\\[1mm]
&\text{symmetric lowest twin enough} &&\iff&& 1-2\varrho < \epsilon_* \le 1-\varrho,\\[1mm]
&\text{non-twin asymmetry required} &&\iff&& \epsilon_* > 1-\varrho.
\end{aligned}
}
```

So the support-selection side already constrains the last live quartic dominance
ambiguity much more strongly than Step 21 by itself.

---

## Inputs from the moving-throat support theorem

From the coherent support-compensation stages,

```math
S_{\rm req} = \frac{\Pi_{\rm tr}}{C_{\rm mix}},
\qquad
C_{\rm mix}=\frac{8\Lambda(1-\epsilon_*)}{\pi^2}.
```

The exact regime split is

```math
\Pi_{\rm tr} \le C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{mixed-only branch already enough,}
```

```math
C_{\rm mix}<\Pi_{\rm tr}\le 2C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{symmetric lowest twin support already enough,}
```

```math
\Pi_{\rm tr} > 2C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{non-twin asymmetry is required.}
```

Those are the exact reduced-level support-selection facts carried forward from the
moving-throat notes.

---

## Step 22A — One scalar support-demand selector

Define

```math
\boxed{
\varrho := \frac{\pi^2\Pi_{\rm tr}}{16\Lambda}.
}
```

Then

```math
\Pi_{\rm tr} = \frac{16\Lambda\varrho}{\pi^2},
```

and the exact required support enhancement becomes

```math
S_{\rm req}
=\frac{\Pi_{\rm tr}}{C_{\rm mix}}
=\frac{2\varrho}{1-\epsilon_*}.
```

So the support theorem is now phrased in the simplest possible way: one scalar
`\varrho` competes against the split-blocking variable `\epsilon_*`.

---

## Step 22B — Exact `\epsilon_*` windows for the three support regimes

### Mixed-only enough

The condition

```math
\Pi_{\rm tr} \le C_{\rm mix}
```

becomes

```math
\frac{16\Lambda\varrho}{\pi^2}
\le
\frac{8\Lambda(1-\epsilon_*)}{\pi^2},
```

hence

```math
\boxed{
\epsilon_* \le 1 - 2\varrho.
}
```

So the mixed lane alone is sufficient only when the branch demand is low enough
that `\epsilon_*` lies below that ceiling.

### Symmetric lowest twin enough

The condition

```math
\Pi_{\rm tr} \le 2 C_{\rm mix}
```

becomes

```math
\boxed{
\epsilon_* \le 1 - \varrho.
}
```

Therefore the exact twin-sufficient window is

```math
\boxed{
1-2\varrho < \epsilon_* \le 1-\varrho.
}
```

### Non-twin asymmetry required

As soon as

```math
\Pi_{\rm tr} > 2 C_{\rm mix},
```

one has

```math
\boxed{
\epsilon_* > 1-\varrho,
}
```

and the symmetric lowest-twin support branch is no longer enough.

So the exact support selector is

```math
\boxed{
\begin{aligned}
&\text{mixed-only enough} &&\iff&& 0<\epsilon_* \le 1-2\varrho,\\[1mm]
&\text{symmetric lowest twin enough} &&\iff&& 1-2\varrho < \epsilon_* \le 1-\varrho,\\[1mm]
&\text{non-twin asymmetry required} &&\iff&& \epsilon_* > 1-\varrho.
\end{aligned}
}
```

### Existence conditions

Two immediate corollaries are exact.

- The mixed-only window is nonempty only if
  ```math
  \varrho < \frac12.
  ```
- The twin-sufficient window is nonempty only if
  ```math
  \varrho < 1.
  ```
- If
  ```math
  \varrho \ge 1,
  ```
  then every constructive branch point already lies in the non-twin-required
  regime.

So `\varrho` is already a branch classifier.

---

## Step 22C — Exact `\sigma` windows

The coherent quartic microledger does not use `\epsilon_*` directly. It uses

```math
\sigma = \frac{2\epsilon_*}{1-\epsilon_*}.
```

So the support selector is most useful after conversion to `\sigma`.

### Mixed-only ceiling

Insert `\epsilon_* = 1-2\varrho`:

```math
\sigma_{\rm mix,max}
=\frac{2(1-2\varrho)}{2\varrho}
=\boxed{\frac{1}{\varrho}-2}.
```

### Twin-sufficient ceiling

Insert `\epsilon_* = 1-\varrho`:

```math
\sigma_{\rm twin,max}
=\frac{2(1-\varrho)}{\varrho}
=\boxed{\frac{2}{\varrho}-2}.
```

So the exact `\sigma` windows are

```math
\boxed{
\begin{aligned}
&\text{mixed-only enough} &&\iff&& 0<\sigma \le \frac{1}{\varrho}-2,\\[1mm]
&\text{symmetric lowest twin enough} &&\iff&& \frac{1}{\varrho}-2 < \sigma \le \frac{2}{\varrho}-2,\\[1mm]
&\text{non-twin asymmetry required} &&\iff&& \sigma > \frac{2}{\varrho}-2.
\end{aligned}
}
```

That is the cleanest Step-22 bridge back into the quartic anomaly ledger.

---

## Step 22D — Intersecting the support selector with the Step-21 `q_W/q_\Lambda` crossover

Step 21 proved

```math
q_W = q_\Lambda
\iff
\epsilon_* = \epsilon_\times := \frac{1}{2+\beta^2},
```

equivalently

```math
\sigma_\times = \frac{2}{1+\beta^2}.
```

So the support-selection question is now simple:

> does the support-compatible `\epsilon_*` interval ever rise above
> `\epsilon_\times`?

### Mixed-only branch

The mixed-only branch can reach `q_\Lambda > q_W` only if its upper ceiling lies
above the crossover:

```math
1-2\varrho > \frac{1}{2+\beta^2}.
```

Solving gives the exact mixed-only threshold

```math
\boxed{
\varrho < \varrho_{\rm mix}^\times
:= \frac{1+\beta^2}{2(2+\beta^2)}.
}
```

So if

```math
\varrho \ge \varrho_{\rm mix}^\times,
```

then every mixed-only-compatible branch point is forced into

```math
q_W > q_\Lambda.
```

### Symmetric-twin-compatible branch

The twin-sufficient branch can reach `q_\Lambda > q_W` only if its upper ceiling
lies above the same crossover:

```math
1-\varrho > \frac{1}{2+\beta^2}.
```

So the exact twin threshold is

```math
\boxed{
\varrho < \varrho_{\rm twin}^\times
:= \frac{1+\beta^2}{2+\beta^2}.
}
```

Hence if

```math
\varrho \ge \varrho_{\rm twin}^\times,
```

then even the whole mixed-only-plus-twin-compatible support sector lies below the
Step-21 crossover, and the quartic branch is forced into

```math
q_W > q_\Lambda.
```

---

## Step 22E — Exact `(\varrho,\beta)` phase diagram for the last dominance ambiguity

The support selector and the Step-21 crossover now combine into three exact
phases.

### Phase A — strong enough support demand that `q_W` always wins on every support-compatible branch

If

```math
\boxed{
\varrho \ge \varrho_{\rm twin}^\times,
}
```

then even the twin ceiling

```math
\epsilon_* = 1-\varrho
```

sits below `\epsilon_\times`, so every mixed-only or symmetric-twin-compatible
branch obeys

```math
\boxed{q_W > q_\Lambda.}
```

### Phase B — mixed-only still wall-dominated, but the twin window can cross into `q_Λ` dominance

If

```math
\boxed{
\varrho_{\rm mix}^\times \le \varrho < \varrho_{\rm twin}^\times,
}
```

then

- the mixed-only branch is forced into `q_W > q_\Lambda`, but
- the upper part of the symmetric-twin window can cross into `q_\Lambda > q_W`.

So in this band, outgoing-scale dominance is possible only with support help or
on a still-harder non-twin branch.

### Phase C — low enough support demand that even the mixed-only branch can reach `q_Λ > q_W`

If

```math
\boxed{
\varrho < \varrho_{\rm mix}^\times,
}
```

then even the mixed-only branch can enter the Step-21 strong-blocking regime
where

```math
q_\Lambda > q_W.
```

So the last live quartic dominance ambiguity is now sharply tied to one support
selector.

---

## Step 22F — Numerical size of the thresholds on the constructive coherent branch

Because

```math
0 < \beta < \frac{2}{11},
```

these support-demand thresholds live in very narrow windows:

```math
\boxed{
\frac14 < \varrho_{\rm mix}^\times < \frac{125}{492} \approx 0.25407,
}
```

```math
\boxed{
\frac12 < \varrho_{\rm twin}^\times < \frac{125}{246} \approx 0.50813.
}
```

So the physics is clean.

- Once the normalized support demand is even moderately large,
  ```math
  \varrho \gtrsim 0.5,
  ```
  all mixed-only and symmetric-twin-compatible branches are forced into the
  `q_W`-dominant phase.
- Only at relatively low demand,
  ```math
  \varrho \lesssim 0.25,
  ```
  can the mixed-only branch itself reach `q_\Lambda > q_W`.

That is much sharper than the purely parametric Step-21 phase diagram.

---

## Main result of the step

The branch-selection side now contributes a real selector, not just more symbols.

The exact support-demand variable

```math
\varrho := \frac{\pi^2\Pi_{\rm tr}}{16\Lambda}
```

fixes the coherent support regimes by

```math
\begin{aligned}
&\text{mixed-only enough} &&\iff&& 0<\epsilon_* \le 1-2\varrho,\\
&\text{symmetric lowest twin enough} &&\iff&& 1-2\varrho < \epsilon_* \le 1-\varrho,\\
&\text{non-twin asymmetry required} &&\iff&& \epsilon_* > 1-\varrho,
\end{aligned}
```

or equivalently

```math
\begin{aligned}
&\text{mixed-only enough} &&\iff&& 0<\sigma \le \frac{1}{\varrho}-2,\\
&\text{symmetric lowest twin enough} &&\iff&& \frac{1}{\varrho}-2 < \sigma \le \frac{2}{\varrho}-2,\\
&\text{non-twin asymmetry required} &&\iff&& \sigma > \frac{2}{\varrho}-2.
\end{aligned}
```

Intersecting that with the Step-21 crossover shows:

```math
\boxed{
\begin{aligned}
&\text{mixed-only can reach } q_\Lambda > q_W
&&\iff&&
\varrho < \frac{1+\beta^2}{2(2+\beta^2)},\\[1mm]
&\text{twin-compatible branch can reach } q_\Lambda > q_W
&&\iff&&
\varrho < \frac{1+\beta^2}{2+\beta^2}.
\end{aligned}
}
```

So the branch-selection side has now turned the last live quartic dominance
ambiguity into one exact support-demand phase diagram.

For moderate or large support demand, every support-compatible coherent branch is
forced into the

```math
q_W > q_\Lambda
```

phase. The `q_\Lambda > q_W` regime is then a low-demand corner or a genuinely
non-twin effect.

---

## Continuation point

The next clean move is now very focused:

> derive or estimate the actual selected-branch demand parameter
> ```math
> \varrho_* = \frac{\pi^2\Pi_{{\rm tr},*}}{16\Lambda_*}
> ```
> from the moving-throat normalization side.

Once that is known, the support regime and the `q_W/q_\Lambda` dominance phase
are no longer qualitative at all; they become a direct branch verdict.
# Step 23 — Selected-branch loading ratio from the minimal isotropic quadrupole precursor

## Goal

Step 22 reduced the last support-selection ambiguity to the scalar demand parameter

```math
\varrho := \frac{\pi^2\Pi_{\rm tr}}{16\Lambda}.
```

But that still left one question open:

> **what does the selected passive/outgoing normalization side actually choose for
> `\Pi_{\rm tr}`?**

The later moving-throat notes answer that in two stages.

1. First, the selected-branch quadrupole-demand product **cancels** all separate
   outgoing-normalization amplitudes and depends only on the loading ratio
   ```math
   \rho_\alpha := \frac{\alpha_{\rm req}}{\alpha_{\rm mix}}.
   ```
2. Second, the natural minimal isotropic conservative quadrupole precursor fixes
   that loading ratio exactly through its contact/pole split.

So Step 23 is the point where the support selector `\varrho` is finally tied back
into the selected-branch normalization side.

The main result is:

```math
\boxed{
\rho_\alpha = \frac43,
\qquad
\zeta_{\rm req}=\frac13,
\qquad
\Pi_{\rm tr}=\frac43\,C_{\rm mix}.
}
```

Equivalently,

```math
\boxed{
\varrho = \frac{2(1-\epsilon_*)}{3},
\qquad
S_{\rm req}=\frac{\Pi_{\rm tr}}{C_{\rm mix}}=\frac43.
}
```

So the natural minimal isotropic passive/outgoing branch is **not** mixed-only and
**not** non-twin. It sits exactly on the symmetric-lowest-twin support slice.

---

## Inputs from the selected-branch normalization side

The later moving-throat notes give the exact product identities

```math
\Pi_{\rm tr}
=
\frac{N_Q^{(\rm target)}}{\beta_0}\,\alpha_{\rm req},
\qquad
C_{\rm mix}
=
\frac{N_Q^{(\rm target)}}{\beta_0}\,\alpha_{\rm mix}.
```

So immediately

```math
\boxed{
\frac{\Pi_{\rm tr}}{C_{\rm mix}}
=
\frac{\alpha_{\rm req}}{\alpha_{\rm mix}}
=:\rho_\alpha.
}
```

This is the key simplification: once the outgoing quadrupole branch is normalized,
all the separate selected-mode amplitudes drop out of the support test.

In the spectral notation of the selected branch,

```math
N_Q^{(\rm target)} = \hat m_-^2\,\beta_0\,\frac{s_-}{\lambda_-},
```

so the same identities can also be written as

```math
\Pi_{\rm tr} = \hat m_-^2\frac{s_-}{\lambda_-}\alpha_{\rm req},
\qquad
C_{\rm mix} = \hat m_-^2\frac{s_-}{\lambda_-}\alpha_{\rm mix}.
```

Again the ratio is just `\rho_\alpha`.

---

## Step 23A — Exact contact-plus-pole inverse formulas

The natural minimal conservative quadrupole precursor is written as

```math
Y_Q^{\rm cons}(\omega)
=
 c_0 + \frac{c_1}{1-\omega^2/\Omega_Q^2},
```

with normalized static limit

```math
c_0 + c_1 = 1.
```

On the explicit support/source branch, the natural reading is:

- the mixed baseline carries the static contact fraction,
- the extra support lane carries the finite conservative pole.

So the same precursor can be written as

```math
Y_Q^{\rm cons}(\omega)
=
\frac{\alpha_{\rm mix}}{\alpha_{\rm req}}
+
\frac{\alpha_{\rm req}-\alpha_{\rm mix}}{\alpha_{\rm req}}
\frac{1}{1-\omega^2/\Omega_Q^2}.
```

Introducing

```math
\rho_\alpha := \frac{\alpha_{\rm req}}{\alpha_{\rm mix}},
```

this becomes

```math
Y_Q^{\rm cons}(\omega)
=
\frac{1}{\rho_\alpha}
+
\frac{\rho_\alpha-1}{\rho_\alpha}
\frac{1}{1-\omega^2/\Omega_Q^2}.
```

Therefore the exact contact/pole data are

```math
\boxed{c_0 = \frac{1}{\rho_\alpha}},
\qquad
\boxed{c_1 = \frac{\rho_\alpha-1}{\rho_\alpha}},
```

and the inverse formulas are

```math
\boxed{\rho_\alpha = \frac{1}{c_0} = \frac{1}{1-c_1}},
```

```math
\boxed{\zeta_{\rm req} := \rho_\alpha-1 = \frac{c_1}{c_0}.}
```

So the support/source loading ratio is directly encoded in the static contact /
pole split of the conservative quadrupole precursor.

---

## Step 23B — Matching to the minimal isotropic quadrupole module

The 2.5PN quadrupole audit already fixed the smallest viable isotropic
conservative precursor to

```math
c_0 = \frac34,
\qquad
c_1 = \frac14,
\qquad
\Omega_Q = \frac{3c_s}{2a}.
```

Inserting those values into the exact inverse formulas gives immediately

```math
\boxed{\rho_\alpha = \frac{1}{3/4} = \frac43,}
```

```math
\boxed{\zeta_{\rm req} = \frac{1/4}{3/4} = \frac13.}
```

Then the selected demand product is

```math
\boxed{
\Pi_{\rm tr}
=
\rho_\alpha C_{\rm mix}
=
\frac43 C_{\rm mix}.
}
```

So the selected branch is no longer carrying an arbitrary support demand. The
natural minimal isotropic passive/outgoing branch fixes it exactly.

---

## Step 23C — Exact support-selector form of the selected branch

Step 22 defined

```math
\varrho := \frac{\pi^2\Pi_{\rm tr}}{16\Lambda},
\qquad
C_{\rm mix} = \frac{8\Lambda(1-\epsilon_*)}{\pi^2}.
```

Substituting

```math
\Pi_{\rm tr} = \frac43 C_{\rm mix}
```

gives

```math
\varrho
=
\frac{\pi^2}{16\Lambda}\cdot\frac43\cdot\frac{8\Lambda(1-\epsilon_*)}{\pi^2}
=
\boxed{\frac{2(1-\epsilon_*)}{3}.}
```

And the required support enhancement is simply

```math
S_{\rm req}
=
\frac{\Pi_{\rm tr}}{C_{\rm mix}}
=
\boxed{\frac43.}
```

So the selected branch is no longer scanning all support-demand sectors. It is
locked to one exact support ratio.

---

## Step 23D — Regime meaning

Stage 22 already split the support regimes by

```math
\Pi_{\rm tr} \le C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{mixed-only enough},
```

```math
C_{\rm mix} < \Pi_{\rm tr} \le 2C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{symmetric lowest twin enough},
```

```math
\Pi_{\rm tr} > 2C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{non-twin asymmetry required}.
```

Because the selected branch gives

```math
\Pi_{\rm tr} = \frac43 C_{\rm mix},
```

it follows exactly that

```math
\boxed{
C_{\rm mix} < \Pi_{\rm tr} < 2C_{\rm mix}.
}
```

So:

- mixed-only is **not** enough,
- the symmetric lowest twin **is** enough,
- and non-twin asymmetry is **not** required.

This is already a real simplification of the anomaly bridge.

---

## Main result of the step

The selected-branch normalization side has now fixed the support ratio carried by
the natural minimal isotropic passive/outgoing quadrupole branch:

```math
\boxed{
\rho_\alpha = \frac43,
\qquad
\zeta_{\rm req}=\frac13,
\qquad
\Pi_{\rm tr}=\frac43 C_{\rm mix}.
}
```

Equivalently,

```math
\boxed{
\varrho = \frac{2(1-\epsilon_*)}{3},
\qquad
S_{\rm req}=\frac43.
}
```

So the last support ambiguity has collapsed from three sectors

- mixed-only,
- symmetric lowest twin,
- non-twin asymmetry,

to exactly **one** selected support slice:

```math
\text{symmetric lowest twin, with demand ratio } \Pi_{\rm tr}/C_{\rm mix}=4/3.
```

---

## What the next step should be

The next honest move is now very sharp:

> restrict the Step-21 primitive quartic ranking problem to this selected
> twin-support branch, and ask how much of the remaining `q_W` versus
> `q_\Lambda` ambiguity survives there.

That is the smallest next derivation that still genuinely pushes the anomaly
closure forward.
# Step 24 — Exact primitive ranking on the selected twin-support branch

## Goal

Step 23 made the main support-side simplification:

```math
\frac{\Pi_{\rm tr}}{C_{\rm mix}} = \frac43.
```

So the natural minimal isotropic passive/outgoing quadrupole branch is **not**
allowed to roam over all three support sectors anymore. It lives on exactly one
selected support slice:

```math
\text{symmetric lowest twin, with } \Pi_{\rm tr}/C_{\rm mix}=4/3.
```

That means the old Step-21/22 phase diagram is now overkill. The real next
question is narrower:

> **once we restrict to that selected twin-support curve, how much of the
> primitive quartic ranking ambiguity survives?**

This step answers that exactly.

The main result is that the whole selected branch is one exact curve

```math
\boxed{\epsilon_* = 1 - \frac{3\varrho}{2}},
\qquad
\boxed{\sigma = \frac{4}{3\varrho}-2},
\qquad
0<\varrho<\frac23,
```

and along that curve only **two** ranking thresholds remain:

```math
\boxed{
\varrho_{W\Lambda}
=
\frac{2(1+\beta^2)}{3(2+\beta^2)},
}
```

```math
\boxed{
\varrho_{U\Lambda}
=
\frac{2(1+\beta^2)}{3(1+\beta+\beta^2)}.
}
```

So the full selected-branch primitive ranking is:

```math
\boxed{
\begin{aligned}
&0<\varrho<\varrho_{W\Lambda}
&&\Longrightarrow&&
q_\chi > q_Z > q_\Lambda > q_W > |q_U|,\\[1mm]
&\varrho_{W\Lambda}<\varrho<\varrho_{U\Lambda}
&&\Longrightarrow&&
q_\chi > q_Z > q_W > q_\Lambda > |q_U|,\\[1mm]
&\varrho_{U\Lambda}<\varrho<\frac23
&&\Longrightarrow&&
q_\chi > q_Z > q_W > |q_U| > q_\Lambda.
\end{aligned}
}
```

That is the cleanest anomaly ranking statement reached so far.

---

## Step 24A — The selected branch is an exact one-parameter twin-support curve

Step 22 defined

```math
\varrho := \frac{\pi^2\Pi_{\rm tr}}{16\Lambda},
\qquad
\sigma = \frac{2\epsilon_*}{1-\epsilon_*}.
```

Step 23 then fixed the selected support ratio to

```math
\frac{\Pi_{\rm tr}}{C_{\rm mix}} = \frac43,
\qquad
C_{\rm mix}=\frac{8\Lambda(1-\epsilon_*)}{\pi^2}.
```

So

```math
\varrho
=
\frac{\pi^2\Pi_{\rm tr}}{16\Lambda}
=
\frac{\pi^2}{16\Lambda}\cdot\frac43\cdot\frac{8\Lambda(1-\epsilon_*)}{\pi^2}
=
\frac{2(1-\epsilon_*)}{3}.
```

Hence

```math
\boxed{\epsilon_* = 1 - \frac{3\varrho}{2}.}
```

Since `0<\epsilon_*<1`, this gives the exact selected-branch range

```math
\boxed{0<\varrho<\frac23.}
```

Now convert to `\sigma`:

```math
\sigma
=
\frac{2\epsilon_*}{1-\epsilon_*}
=
\frac{2\bigl(1-3\varrho/2\bigr)}{3\varrho/2}
=
\boxed{\frac{4}{3\varrho}-2.}
```

So the selected branch is not a 2D region in `(epsilon_*,\varrho)` at all. It
is a single exact curve.

---

## Step 24B — The selected curve sits strictly inside the twin window

Step 22 gave the exact support windows in `\sigma`:

```math
0<\sigma\le \frac1\varrho-2
\quad\Longleftrightarrow\quad
\text{mixed-only enough},
```

```math
\frac1\varrho-2 < \sigma \le \frac2\varrho-2
\quad\Longleftrightarrow\quad
\text{symmetric lowest twin enough},
```

```math
\sigma > \frac2\varrho-2
\quad\Longleftrightarrow\quad
\text{non-twin asymmetry required}.
```

On the selected branch,

```math
\sigma_{\rm sel} = \frac{4}{3\varrho}-2.
```

Then

```math
\sigma_{\rm sel} - \left(\frac1\varrho-2\right)
=
\frac{1}{3\varrho} > 0,
```

and

```math
\left(\frac2\varrho-2\right) - \sigma_{\rm sel}
=
\frac{2}{3\varrho} > 0.
```

So for every allowed point on the selected branch,

```math
\boxed{
\frac1\varrho-2 < \sigma_{\rm sel} < \frac2\varrho-2.
}
```

That is the exact proof that the selected branch lies **strictly inside** the
symmetric-lowest-twin regime.

So mixed-only and non-twin branches are gone from the live anomaly closure.

---

## Step 24C — Surviving threshold 1: `q_W` versus `q_\Lambda`

Step 21 gave the exact crossover

```math
q_W = q_\Lambda
\iff
\epsilon_* = \frac{1}{2+\beta^2}.
```

Insert the selected-branch law

```math
\epsilon_* = 1 - \frac{3\varrho}{2}.
```

Then the threshold on the selected branch is

```math
1 - \frac{3\varrho}{2} = \frac{1}{2+\beta^2},
```

hence

```math
\boxed{
\varrho_{W\Lambda}
=
\frac{2(1+\beta^2)}{3(2+\beta^2)}.
}
```

Therefore:

- if
  ```math
  0<\varrho<\varrho_{W\Lambda},
  ```
  then `q_\Lambda > q_W`;
- if
  ```math
  \varrho>\varrho_{W\Lambda},
  ```
  then `q_W > q_\Lambda`.

So the outgoing-scale lane overtakes the wall-blocking lane only in the **low-`\varrho` / high-blocking** corner of the selected curve.

---

## Step 24D — Surviving threshold 2: `|q_U|` versus `q_\Lambda`

Step 21 also gave the exact crossover

```math
|q_U| = q_\Lambda
\iff
\epsilon_* = \frac{\beta}{1+\beta+\beta^2}.
```

Again insert

```math
\epsilon_* = 1 - \frac{3\varrho}{2},
```

and solve for `\varrho`:

```math
\boxed{
\varrho_{U\Lambda}
=
\frac{2(1+\beta^2)}{3(1+\beta+\beta^2)}.
}
```

So:

- if
  ```math
  \varrho<\varrho_{U\Lambda},
  ```
  then `q_\Lambda > |q_U|`;
- if
  ```math
  \varrho>\varrho_{U\Lambda},
  ```
  then `|q_U| > q_\Lambda`.

This is the selected-branch version of Step 21’s “very weak blocking” corner.

---

## Step 24E — Ordering of the two thresholds

The two thresholds are not independent. Their difference is

```math
\varrho_{U\Lambda} - \varrho_{W\Lambda}
=
\frac{2(1+\beta^2)(1-\beta)}{3(1+\beta+\beta^2)(2+\beta^2)} > 0
```

because `0<\beta<2/11<1`.

And

```math
\frac23 - \varrho_{U\Lambda}
=
\frac{2\beta}{3(1+\beta+\beta^2)} > 0.
```

So the exact threshold ordering on the selected branch is

```math
\boxed{0 < \varrho_{W\Lambda} < \varrho_{U\Lambda} < \frac23.}
```

That means the selected twin-support curve always splits into **three** ranking
regions and never fewer.

---

## Step 24F — Full primitive ranking on the selected twin-support branch

Step 21 already proved the branch-independent ordering facts

```math
q_\chi > q_Z,
\qquad
q_Z > q_W,
\qquad
q_W > |q_U|.
```

So only `q_\Lambda` moves relative to `q_W` and `|q_U|`.
Using the two selected-branch thresholds above, the complete ranking is now exact.

### Region I — low `\varrho`, strong blocking

If

```math
0<\varrho<\varrho_{W\Lambda},
```

then

```math
\boxed{q_\chi > q_Z > q_\Lambda > q_W > |q_U|.}
```

### Region II — intermediate `\varrho`

If

```math
\varrho_{W\Lambda}<\varrho<\varrho_{U\Lambda},
```

then

```math
\boxed{q_\chi > q_Z > q_W > q_\Lambda > |q_U|.}
```

### Region III — large `\varrho`, very weak blocking

If

```math
\varrho_{U\Lambda}<\varrho<\frac23,
```

then

```math
\boxed{q_\chi > q_Z > q_W > |q_U| > q_\Lambda.}
```

So the selected anomaly branch now has a completely explicit primitive ranking
phase diagram.

---

## Step 24G — Numerical size of the surviving thresholds

Using the constructive coherent bound

```math
0<\beta<\frac{2}{11},
```

one finds

```math
\boxed{
\frac13 < \varrho_{W\Lambda} < \frac{125}{369} \approx 0.338753,
}
```

and

```math
\boxed{
\frac{250}{441} \approx 0.566893 < \varrho_{U\Lambda} < \frac23.
}
```

So the selected twin-support curve has a very clean geometry.

- Only the **low-`\varrho`** end allows `q_\Lambda` to beat `q_W`.
- Across the middle of the selected curve, `q_W` beats `q_\Lambda` but
  `q_\Lambda` still beats `|q_U|`.
- Only near the **large-`\varrho` / very weak-blocking** end does `|q_U|`
  overtake `q_\Lambda`.

That is already much sharper than the full constructive-branch picture from Step 21.

---

## Main result of the step

The natural minimal isotropic passive/outgoing branch has collapsed the old
support-selection problem to one exact twin-support curve:

```math
\boxed{\epsilon_* = 1 - \frac{3\varrho}{2}},
\qquad
\boxed{\sigma = \frac{4}{3\varrho}-2},
\qquad
0<\varrho<\frac23.
```

On that curve, the primitive quartic hierarchy is controlled by only two exact
thresholds:

```math
\boxed{
\varrho_{W\Lambda}
=
\frac{2(1+\beta^2)}{3(2+\beta^2)},
\qquad
\varrho_{U\Lambda}
=
\frac{2(1+\beta^2)}{3(1+\beta+\beta^2)}.
}
```

So the complete selected-branch ranking is

```math
\boxed{
\begin{aligned}
&0<\varrho<\varrho_{W\Lambda}
&&\Longrightarrow&&
q_\chi > q_Z > q_\Lambda > q_W > |q_U|,\\[1mm]
&\varrho_{W\Lambda}<\varrho<\varrho_{U\Lambda}
&&\Longrightarrow&&
q_\chi > q_Z > q_W > q_\Lambda > |q_U|,\\[1mm]
&\varrho_{U\Lambda}<\varrho<\frac23
&&\Longrightarrow&&
q_\chi > q_Z > q_W > |q_U| > q_\Lambda.
\end{aligned}
}
```

This is the strongest quartic anomaly ranking statement reached so far from the
selected moving-throat branch side.

---

## What the next step should be

The next honest move is now very narrow:

> derive the actual physical position of the moving-throat branch on this
> selected twin-support curve — equivalently, pin `\epsilon_*` or `\varrho`
> rather than leaving it parametric.

Once that single coordinate is known, the quartic carrier hierarchy stops being a
phase diagram and becomes one definite anomaly prediction.
# Step 25 — The grouped-`P2` + static-geometry split forces the `3/4 + 1/4` conservative module

## Goal

Step 24 left the quartic anomaly hierarchy on the selected twin-support curve
parameterized by the single selector `\varrho`. That was the right support-side
picture, but it still treated the conservative passive/outgoing quadrupole branch
as if its contact/pole split were an external choice.

The moving-throat notes give a sharper route.

The 3PN conservative split already says that the higher conservative payload is
organized as

```math
\text{grouped real }P_2 \text{ middle block} + \text{static geometry completion},
```

and the 2.5PN quadrupole audit already fixed the minimal isotropic branch
identity

```math
K_0 K_4 = 4 K_2^2.
```

So the next honest question is:

> **if the actual conservative branch is one isotropic grouped-`P2` pole plus a
> purely static geometry completion, what contact/pole split is forced?**

This step answers that exactly.

The main result is

```math
\boxed{K_{\rm geom} = 3 K_{\rm pole}},
```

so

```math
\boxed{K_{\rm pole}=K_0/4,
\qquad
K_{\rm geom}=3K_0/4,}
```

and therefore the normalized conservative quadrupole module is forced to be

```math
\boxed{
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+
\frac{1}{4}
\frac{1}{1-\omega^2/\Omega_Q^2}.
}
```

Once that is true, the support/source loading ratio is no longer free:

```math
\boxed{\rho_\alpha = \frac43,
\qquad
\zeta_{\rm req}=\frac13,
\qquad
\Pi_{\rm tr}=\frac43 C_{\rm mix}.}
```

So this is the step where the old support-side selector stops being the real
bottleneck. The support ratio `4/3` follows as a **corollary** of the grouped
conservative branch organization.

---

## Step 25A — Minimal isotropic grouped-`P2` + geometry realization

Take the conservative isotropic quadrupole module in the minimal form

```math
K_Q^{\rm cons}(\omega)
=
K_{\rm geom}
+
\frac{K_{\rm pole}}{1-\omega^2/\Omega_Q^2}.
```

Here:

- `K_geom` is the static geometry completion,
- `K_pole` is the isotropic grouped-`P2` pole residue,
- `\Omega_Q` is the effective isotropic grouped-`P2` pole.

This is the smallest realization compatible with the already-frozen 3PN split if

- the grouped-`P2` side is the only dynamic conservative quadrupole lane,
- and geometry contributes only the static completion.

Expand through `O(\omega^4)`:

```math
K_Q^{\rm cons}(\omega)
=
K_0 + K_2 \omega^2 + K_4 \omega^4 + O(\omega^6),
```

with exact coefficients

```math
K_0 = K_{\rm geom}+K_{\rm pole},
```

```math
K_2 = \frac{K_{\rm pole}}{\Omega_Q^2},
```

```math
K_4 = \frac{K_{\rm pole}}{\Omega_Q^4}.
```

So the whole conservative branch is now parameterized by two amplitudes and one
pole scale.

---

## Step 25B — The minimal isotropic branch identity fixes the contact/pole split

The 2.5PN audit already froze the minimal isotropic conservative branch identity

```math
K_0 K_4 = 4 K_2^2.
```

Insert the one-pole + static-geometry coefficients:

```math
(K_{\rm geom}+K_{\rm pole})\frac{K_{\rm pole}}{\Omega_Q^4}
=
4\left(\frac{K_{\rm pole}}{\Omega_Q^2}\right)^2.
```

Assuming the branch is nontrivial (`K_pole\neq0`), the common factor cancels and
we get

```math
K_{\rm geom}+K_{\rm pole}=4K_{\rm pole},
```

hence

```math
\boxed{K_{\rm geom}=3K_{\rm pole}.}
```

So equivalently

```math
\boxed{K_{\rm pole}=\frac{K_0}{4},
\qquad
K_{\rm geom}=\frac{3K_0}{4}.}
```

This is the central algebraic result of the step.

The normalized conservative response therefore becomes

```math
\widehat Y_Q^{\rm cons}(\omega)
:=
\frac{K_Q^{\rm cons}(\omega)}{K_0}
=
\frac34+
\frac14\frac{1}{1-\omega^2/\Omega_Q^2}.
```

So the familiar `3/4 + 1/4` split is not an imported guess anymore. Under the
minimal grouped-`P2` + static-geometry realization, it is forced.

---

## Step 25C — Exact contact/pole map to the support/source loading ratio

On the explicit support/source branch, the natural contact-plus-pole reading is

```math
Y_Q^{\rm cons}(\omega)
=
\frac{\alpha_{\rm mix}}{\alpha_{\rm req}}
+
\frac{\alpha_{\rm req}-\alpha_{\rm mix}}{\alpha_{\rm req}}
\frac{1}{1-\omega^2/\Omega_Q^2}.
```

Introduce the selected loading ratio

```math
\rho_\alpha := \frac{\alpha_{\rm req}}{\alpha_{\rm mix}}.
```

Then

```math
Y_Q^{\rm cons}(\omega)
=
\frac{1}{\rho_\alpha}
+
\frac{\rho_\alpha-1}{\rho_\alpha}
\frac{1}{1-\omega^2/\Omega_Q^2}.
```

So the exact contact and pole fractions are

```math
c_0 = \frac{1}{\rho_\alpha},
\qquad
c_1 = \frac{\rho_\alpha-1}{\rho_\alpha},
```

with inverse formulas

```math
\rho_\alpha = \frac{1}{c_0} = \frac{1}{1-c_1},
```

```math
\zeta_{\rm req} := \rho_\alpha-1 = \frac{c_1}{c_0}.
```

Now insert the forced grouped-branch values

```math
c_0 = \frac34,
\qquad
c_1 = \frac14.
```

Then immediately

```math
\boxed{\rho_\alpha = \frac{1}{3/4} = \frac43,}
```

```math
\boxed{\zeta_{\rm req} = \frac{1/4}{3/4} = \frac13.}
```

And since

```math
\frac{\Pi_{\rm tr}}{C_{\rm mix}} = \rho_\alpha,
```

we also get

```math
\boxed{\Pi_{\rm tr} = \frac43 C_{\rm mix}.}
```

So the Stage-23 support ratio is now recovered directly from the grouped
conservative branch organization.

---

## Step 25D — What actually changes

Before this step, the selected twin-support curve looked like the main organizing
object. After this step, the logic reverses:

1. the higher conservative payload is read as
   ```math
   \text{grouped }P_2 + \text{static geometry},
   ```
2. the minimal isotropic branch identity then forces
   ```math
   3/4 + 1/4,
   ```
3. and only **then** does the support/source side inherit
   ```math
   \rho_\alpha = 4/3.
   ```

So the old support-phase picture is no longer the deepest one. The deeper theorem
gate is now:

> **does the actual moving-throat branch really realize one isotropic grouped-`P2`
> pole plus a purely static geometry completion?**

If yes, the support/source verdict follows immediately.

---

## Main result of the step

The minimal isotropic conservative quadrupole module is forced by the grouped
real `P2` + static-geometry split:

```math
\boxed{
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+
\frac14\frac{1}{1-\omega^2/\Omega_Q^2}.
}
```

And the support/source loading ratio follows as a corollary:

```math
\boxed{\rho_\alpha = \frac43,
\qquad
\zeta_{\rm req}=\frac13,
\qquad
\Pi_{\rm tr}=\frac43 C_{\rm mix}.}
```

So the next honest move is no longer to keep scanning support-sector phase
regions. It is to test whether the actual isotropic moving-throat branch really
has that minimal grouped-`P2` + static-geometry conservative structure.
# Step 26 — The actual isotropic branch collapses to one normalization defect

## Goal

Step 25 showed that if the conservative higher-order branch is

```math
\text{one isotropic grouped-}P_2\text{ pole} + \text{static geometry completion},
```

then the conservative quadrupole module is forced to be

```math
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+
\frac14\frac{1}{1-\omega^2/\Omega_Q^2},
```

and the support/source ratio follows as

```math
\rho_\alpha=\frac43.
```

But one real reduced ambiguity was still visible:

> **could the geometry lane contaminate the grouped quadrupole module at
> `O(\omega^2)` or `O(\omega^4)` and therefore spoil the forced `3/4 + 1/4`
> split?**

This step answers that and then pushes the whole branch to its sharpest reduced
form.

The main results are:

1. on the isotropic quadratic wall operator, the `l=0` geometry lane and the
   grouped real `l=2` bundle are exactly orthogonal,
2. therefore the dynamic geometry contamination numbers vanish,
   ```math
   \boxed{\epsilon_2 = \epsilon_4 = 0,}
   ```
3. so the actual isotropic branch really does realize the Step-25
   `3/4 + 1/4` conservative module,
4. and once that is true, the whole remaining reduced 2.5PN/4PN mismatch is just
   one scalar normalization defect
   ```math
   \boxed{N_Q := \bar K_0 / \bar K_0^{\rm target}.}
   ```

At the same time, the support/source side becomes automatic on that branch.

So this is the step where the old selected twin-support phase picture drops out
of the active theorem ledger.

---

## Step 26A — Exact isotropic geometry-decoupling theorem

Take the isotropic quadratic wall operator on the throat sphere. Because it is
`O(3)` invariant, its angular part depends only on the sphere Laplacian:

```math
L_{\rm ang} = a + b(-\Delta_{S^2}).
```

Expand the wall field into

- the scalar geometry lane `Y_{00}`,
- the grouped real quadrupole bundle `Y_{2A}` with
  ```math
  A\in\{20,21c,21s,22c,22s\}.
  ```

Using normalized real harmonics, the exact orthogonality relations are

```math
\int Y_{00} Y_{2A}\, d\Omega = 0.
```

And because

```math
-\Delta_{S^2}Y_{2A}=6Y_{2A},
```

we also have

```math
\langle Y_{00}, (a+b(-\Delta))Y_{2A} \rangle
=
(a+6b)\langle Y_{00},Y_{2A}\rangle
=0.
```

So on the isotropic branch the scalar/geometry `l=0` lane and the grouped real
`l=2` bundle are exactly block diagonal at linear order.

That is the key selection rule.

It means the scalar geometry lane cannot feed dynamic `\omega^2` or `\omega^4`
moments into the grouped quadrupole carrier on the actual isotropic branch.

---

## Step 26B — Dynamic geometry contamination vanishes

The obstruction formula for a geometry lane that carries its own dynamic even
moments is

```math
c_{\rm pole}
=
\frac{1+\epsilon_4}{4(1+\epsilon_2)^2},
```

with

```math
\epsilon_2 = \Omega_Q^2 K_{(g,2)}/K_{\rm pole},
\qquad
\epsilon_4 = \Omega_Q^4 K_{(g,4)}/K_{\rm pole}.
```

But the isotropic decoupling theorem above implies

```math
K_{(g,2)} = K_{(g,4)} = 0,
```

so

```math
\boxed{\epsilon_2 = \epsilon_4 = 0.}
```

Then the obstruction formula collapses to

```math
c_{\rm pole} = \frac14,
\qquad
c_{\rm geom} = 1-c_{\rm pole} = \frac34.
```

So the Step-25 contact/pole split is not merely a candidate branch value. It is
actually the isotropic reduced branch value once the geometry lane is checked.

---

## Step 26C — The full reduced mismatch is one scalar defect

Now write the actual isotropic passive/outgoing one-pole branch in invariant form
as

```math
\bar K_Q^{\rm cons}(\omega)
=
\bar K_0\left[
\frac34+
\frac14\frac{1}{1-\omega^2/\Omega_Q^2}
\right].
```

Define the GR target normalization

```math
\bar K_0^{\rm target} = \frac{64 G\Omega_Q^5}{45 c^5}.
```

Using the already-carried geometric pole

```math
\Omega_Q = \frac{3c_s}{2a},
```

this is equivalently

```math
\bar K_0^{\rm target}
=
\frac{54Gc_s^5}{5a^5c^5}.
```

Now define the actual branch normalization defect

```math
\boxed{N_Q := \frac{\bar K_0}{\bar K_0^{\rm target}}.}
```

Then every low-frequency invariant scales by the same factor:

```math
\bar K_2 = N_Q\,\bar K_2^{\rm target},
```

```math
\bar K_4 = N_Q\,\bar K_4^{\rm target},
```

```math
\bar\Gamma_5 = N_Q\,\frac{2G}{5c^5}.
```

So the full reduced GR-like point-particle 2.5PN closure on the actual isotropic
branch is now equivalent to the one equation

```math
\boxed{N_Q = 1.}
```

That is the narrowest reduced theorem gate reached so far.

---

## Step 26D — The support/source side becomes automatic

From Step 25 the actual isotropic branch has

```math
\rho_\alpha = \frac43.
```

The exact blocked support demand therefore becomes

```math
\zeta_{\rm req}^{\rm(act)}(\epsilon_{\rm blk})
=
\frac{\rho_\alpha-1}{1-\epsilon_{\rm blk}(2-\rho_\alpha)}
=
\frac{1}{3-2\epsilon_{\rm blk}}.
```

Now take any explicit support/source family with constructive ceiling

```math
\zeta_{\max} > 1
```

and admissible blocking window

```math
0 \le \epsilon_{\rm blk} < 1/\zeta_{\max}.
```

Then the worst blocked demand obeys

```math
\zeta_{\rm req}^{\rm(act)}
<
\frac{1}{3-2/\zeta_{\max}},
```

and

```math
\zeta_{\max} - \frac{1}{3-2/\zeta_{\max}}
=
\frac{3\zeta_{\max}(\zeta_{\max}-1)}{3\zeta_{\max}-2} > 0.
```

So any explicit family with `\zeta_{\max}>1` already passes the actual isotropic
support test throughout its admissible blocked regime.

In particular, the explicit Family-1 branch has

```math
\zeta_{\max}^{\rm(F1)} \approx 2.4675 > 1,
```

so its support/source side is automatic once the actual isotropic branch is
adopted.

---

## Step 26E — What actually changes

After this step, the active reduced theorem ledger is much smaller:

- **no** geometry-lane contamination on the isotropic branch,
- **no** extra contact/pole ambiguity,
- **no** remaining explicit support/source ambiguity,
- only the single passive/outgoing normalization defect `N_Q`.

So the old selected twin-support curve is no longer the live reduced object. It
was useful while the support/source side still looked like an independent
selection problem, but once the actual isotropic grouped-`P2` one-pole branch is
used, the support/source theorem becomes automatic.

The real remaining question is now simply:

> **does the completed moving-throat PDE give `N_Q = 1` on the actual passive/
> outgoing isotropic grouped-`P2` branch?**

---

## Main result of the step

On the actual isotropic grouped-`P2` one-pole branch,

```math
\boxed{\epsilon_2 = \epsilon_4 = 0,}
```

so the conservative quadrupole module is exactly

```math
\boxed{
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+
\frac14\frac{1}{1-\omega^2/\Omega_Q^2}.
}
```

The support/source side is then automatic, and the full remaining reduced
mismatch collapses to one scalar normalization defect,

```math
\boxed{N_Q = \bar K_0 / \bar K_0^{\rm target}.}
```

So the next honest move is no longer to keep refining support-side phase
selection. It is to compute the actual passive/outgoing quadrupole normalization
of the moving-throat branch itself.
# Step 27 — The actual isotropic passive/outgoing branch has one normalization defect

## Goal

Step 26 showed that on the actual isotropic branch the geometry lane is dynamically
inert through `O(ω^4)`, so the conservative grouped-`P2` carrier is exactly the
forced

```math
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2}.
```

The next honest move is to combine that conservative branch with the already-carried
minimal outgoing grouped-`P2` law and ask:

> **how many independent normalization defects are left on the actual isotropic
> passive/outgoing branch?**

This step shows the answer is **one**.

The main result is that, once the actual branch is isotropic and one-pole,

```math
\boxed{N_Q := \bar K_0/\bar K_0^{\rm target}}
```

controls **all** low-frequency invariant mismatches:

```math
\boxed{
\frac{\bar K_2}{\bar K_2^{\rm target}}
=
\frac{\bar K_4}{\bar K_4^{\rm target}}
=
\frac{\bar\Gamma_5}{\bar\Gamma_5^{\rm target}}
=
N_Q.
}
```

So the full reduced 2.5PN/4PN theorem gap on that branch is equivalent to the one
scalar condition

```math
\boxed{N_Q=1.}
```

---

## Step 27A — Exact low-frequency coefficients of the actual conservative module

Start from the actual isotropic grouped-`P2` conservative module

```math
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2}.
```

Expanding at low frequency gives

```math
\widehat Y_Q^{\rm cons}(\omega)
=
1+
\frac{\omega^2}{4\Omega_Q^2}
+
\frac{\omega^4}{4\Omega_Q^4}
+O(\omega^6).
```

So the exact conservative moments are

```math
\boxed{\bar K_2 = \bar K_0/(4\Omega_Q^2),}
```

```math
\boxed{\bar K_4 = \bar K_0/(4\Omega_Q^4).}
```

There is no extra even ambiguity once `\bar K_0` and `\Omega_Q` are fixed.

---

## Step 27B — The odd Burke–Thorne coefficient is not independent either

On the same minimal isotropic outgoing branch, the 2.5PN audit already fixed the
odd coefficient algebraically as

```math
\bar\Gamma_5
=
9\,\bar K_2^{5/2}/\bar K_0^{3/2}.
```

Substituting the conservative one-pole relation for `\bar K_2` gives

```math
\boxed{\bar\Gamma_5 = \frac{9\bar K_0}{32\Omega_Q^5}.}
```

So the odd Burke–Thorne coefficient is tied to the same two quantities
`(\bar K_0,\Omega_Q)` and is not a separate datum.

---

## Step 27C — Exact target branch

On the GR target branch,

```math
\bar K_0^{\rm target} = \frac{64G\Omega_Q^5}{45c^5}.
```

Then automatically

```math
\bar K_2^{\rm target} = \bar K_0^{\rm target}/(4\Omega_Q^2),
```

```math
\bar K_4^{\rm target} = \bar K_0^{\rm target}/(4\Omega_Q^4),
```

```math
\bar\Gamma_5^{\rm target} = \frac{2G}{5c^5}.
```

Using the already-carried geometric pole

```math
\Omega_Q = \frac{3c_s}{2a},
```

this becomes

```math
\boxed{
\bar K_0^{\rm target}
=
\frac{54Gc_s^5}{5a^5c^5}.
}
```

---

## Step 27D — One scalar defect controls the whole branch

Define the actual-branch normalization defect by

```math
\boxed{N_Q := \bar K_0/\bar K_0^{\rm target}.}
```

Then the actual branch can be written as

```math
\bar K_0 = N_Q\,\bar K_0^{\rm target},
```

and therefore

```math
\bar K_2 = N_Q\,\bar K_2^{\rm target},
```

```math
\bar K_4 = N_Q\,\bar K_4^{\rm target},
```

```math
\bar\Gamma_5 = N_Q\,\bar\Gamma_5^{\rm target} = N_Q\,\frac{2G}{5c^5}.
```

So the four natural defect measures

```math
R_0 := \bar K_0/\bar K_0^{\rm target}-1,
```

```math
R_2 := \bar K_2/\bar K_2^{\rm target}-1,
```

```math
R_4 := \bar K_4/\bar K_4^{\rm target}-1,
```

```math
R_5 := \bar\Gamma_5/\bar\Gamma_5^{\rm target}-1,
```

all collapse to the same number:

```math
\boxed{R_0=R_2=R_4=R_5=N_Q-1.}
```

---

## Main result of the step

On the actual isotropic grouped-`P2` one-pole passive/outgoing branch,

```math
\boxed{N_Q := \bar K_0/\bar K_0^{\rm target}}
```

is the **only** reduced normalization defect.

Everything else follows from it:

```math
\boxed{
\frac{\bar K_2}{\bar K_2^{\rm target}}
=
\frac{\bar K_4}{\bar K_4^{\rm target}}
=
\frac{\bar\Gamma_5}{\bar\Gamma_5^{\rm target}}
=
N_Q.
}
```

So the full reduced 2.5PN/4PN theorem gap is now equivalent to

```math
\boxed{N_Q=1.}
```

That is the cleanest conservative normalization statement reached so far.
# Step 28 — Explicit outgoing `l=2` DtN fix of `\chi_Q` and the branch-selection law

## Goal

Step 27 reduced the actual isotropic passive/outgoing grouped-`P2` branch to one
scalar conservative normalization defect,

```math
N_Q := \bar K_0/\bar K_0^{\rm target}.
```

The next honest question is the retarded one:

> **does the actual passive/outgoing branch carry the canonical compact outgoing
> `l=2` odd normalization, or is there still one live outgoing normalization factor?**

This step answers that in two parts.

1. The explicit compact outgoing `l=2` DtN model fixes the last reduced outgoing
   scalar exactly:
   ```math
   \boxed{\chi_Q = 1.}
   ```
2. If the true moving-throat branch is a deformed isotropic DtN branch, the only
   remaining first-order obstruction is the isotropic branch-selection triple
   `(b,a_0,a_5)` through
   ```math
   \boxed{\chi_Q = 1 + 5b + a_0/3 + 9a_5 + O(2).}
   ```

So after this step the remaining reduced uncertainty is no longer vague “outgoing
structure.” It is one explicit DtN branch-selection law.

---

## Step 28A — Exact compact outgoing `l=2` DtN fingerprint

Let

```math
z := \frac{a\omega}{c_s},
```

and take the outgoing partial wave

```math
h_2^{(1)}(z)=j_2(z)+i y_2(z).
```

The exact outgoing `l=2` Dirichlet-to-Neumann operator is

```math
\Lambda_2^{\rm out}(z)
=
z\,\frac{d}{dz}\ln h_2^{(1)}(z)
=
z\,\frac{h_2^{(1)\prime}(z)}{h_2^{(1)}(z)}.
```

Its small-`z` expansion is

```math
\boxed{
\Lambda_2^{\rm out}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}+i\frac{z^5}{9}
-\frac{2z^6}{27}-i\frac{z^7}{27}+O(z^8).
}
```

Normalizing by the static slot gives

```math
\widehat Y_2^{\rm out}(z)
:=
-\frac{3}{\Lambda_2^{\rm out}(z)}
=
1+\frac{z^2}{9}+\frac{4z^4}{81}+i\frac{z^5}{27}
-\frac{11z^6}{729}-i\frac{z^7}{243}+O(z^8).
```

So the canonical compact outgoing odd coefficient is fixed exactly.

---

## Step 28B — Matching to the retarded grouped-`P2` one-pole module

Write the retarded normalized grouped-`P2` one-pole branch as

```math
\widehat Y_Q^{\rm ret}(\omega)
=
\frac34+
\frac14\frac{1}{1-\omega^2/\Omega_Q^2-i\chi_Q\sigma_Q^{\rm can}\omega^5}
+O(\omega^6),
```

with

```math
\sigma_Q^{\rm can}:=\frac{9}{8\Omega_Q^5}.
```

Using the already-carried geometric pole

```math
\Omega_Q = \frac{3c_s}{2a},
```

this becomes

```math
\boxed{\sigma_Q^{\rm can} = \frac{4a^5}{27c_s^5}.}
```

Expanding through `O(ω^5)` gives

```math
\widehat Y_Q^{\rm ret}(\omega)
=
1+rac{a^2\omega^2}{9c_s^2}
+rac{4a^4\omega^4}{81c_s^4}
+i\chi_Q\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6).
```

Comparing with the explicit DtN fingerprint from Step 28A,

```math
\widehat Y_2^{\rm out}(\omega)
=
1+rac{a^2\omega^2}{9c_s^2}
+rac{4a^4\omega^4}{81c_s^4}
+i\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6),
```

forces

```math
\boxed{\chi_Q = 1.}
```

So on the canonical compact outgoing branch the last reduced outgoing normalization
scalar is fixed exactly.

---

## Step 28C — Exact factorization of the last reduced 2.5PN defect

The reduced normalization stack already had the exact factorized condition

```math
\hat m_0^{\,2}\,\chi_Q\,N_Q=1.
```

On the natural point-particle source-map branch,

```math
\hat m_0 = 1 + O(a^2/r^2),
```

so in the strict point-particle limit,

```math
N_Q = 1/\chi_Q.
```

Therefore on the canonical compact outgoing DtN branch,

```math
\boxed{N_Q = 1.}
```

So the reduced nonspinning point-particle 2.5PN theorem is closed on that
canonical outgoing branch.

---

## Step 28D — Exact isotropic DtN deformation law

To isolate the remaining PDE-facing freedom, deform the outgoing DtN operator by

```math
\Lambda_2^{\rm def}(z)
=
S\,\Lambda_2^{\rm out}(\beta z)
+\Sigma_0+\Sigma_2 z^2+\Sigma_4 z^4+i\Sigma_5 z^5+O(z^6).
```

Then the low-frequency coefficients are

```math
L_0 = -3S + \Sigma_0,
```

```math
L_2 = S\beta^2/3 + \Sigma_2,
```

```math
L_4 = S\beta^4/9 + \Sigma_4,
```

```math
L_5 = S\beta^5/9 + \Sigma_5.
```

Normalizing by the actual static slot,

```math
\widehat Y_2^{\rm def}(z)=\frac{L_0}{L_0+L_2 z^2+L_4 z^4+iL_5 z^5}+O(z^6),
```

and demanding that the even fingerprint stay canonical,

```math
\frac{z^2}{9},\qquad \frac{4z^4}{81},
```

fixes

```math
\Sigma_2=-\frac{3S\beta^2-3S+\Sigma_0}{9},
```

```math
\Sigma_4=-\frac{3S\beta^4-3S+\Sigma_0}{27}.
```

With those imposed, the outgoing normalization becomes

```math
\boxed{
\chi_Q=
\frac{3\bigl(S\beta^5+9\Sigma_5\bigr)}{3S-\Sigma_0}.
}
```

So the only isotropic branch data that can move the canonical outgoing value are:

- an argument deformation `β`,
- a static additive slot `Σ_0`,
- and an odd core outlet `Σ_5`.

---

## Step 28E — Linearized branch-selection law

Write a small deformation of the canonical branch as

```math
S = 1 + \varepsilon s,
\qquad
\beta = 1 + \varepsilon b,
\qquad
\Sigma_0 = \varepsilon a_0,
\qquad
\Sigma_5 = \varepsilon a_5.
```

Then the exact deformation law expands to

```math
\boxed{
\chi_Q
=
1+
\varepsilon\left(5b+\frac{a_0}{3}+9a_5\right)
+O(\varepsilon^2).
}
```

So at first order:

- pure overall mouth rescaling drops out,
- effective argument deformation contributes as `5b`,
- the static additive slot contributes as `a_0/3`,
- the odd core outlet contributes as `9a_5`.

Equivalently, the linearized preservation condition is

```math
\boxed{5b+\frac{a_0}{3}+9a_5=0.}
```

---

## Main result of the step

The explicit compact outgoing `l=2` DtN branch fixes the last reduced outgoing
scalar exactly:

```math
\boxed{\chi_Q=1.}
```

So on the natural point-particle source-map branch,

```math
\boxed{N_Q=1.}
```

Any remaining deviation from the canonical outgoing branch is an actual isotropic
DtN branch-selection effect, and at first order it can only enter through the
triple

```math
\boxed{(b,a_0,a_5),}
```

via

```math
\boxed{
\chi_Q
=
1+
5b+\frac{a_0}{3}+9a_5
+O(2).
}
```

That is the sharpest retarded closure statement reached so far.
# Step 29 — The quartic g-2 sliver is a tiny finite-`f` outgoing-normalization drift

## Goal

Step 28 fixed the canonical compact outgoing grouped-`P2` branch and showed that,
on the natural source-map branch,

```math
N_Q-1 = -\frac{\Delta_Q}{1+\Delta_Q},
\qquad
\Delta_Q := \chi_Q - 1.
```

The next honest question is:

> **can the remaining quartic electron-anomaly sliver be carried by a tiny
> finite-`f` drift of that same outgoing normalization, without spoiling the
> strict point-particle gravitational branch?**

This step shows the answer is yes.

The main result is that the missing quartic layer can be rewritten as an
`O(f)` outgoing-normalization defect:

```math
\boxed{
\Delta_Q(f)=\delta_1 f + O(f^2),
\qquad
\delta_1 = -\Lambda_1,
}
```

where `\Lambda_1` is the Step-2 common-path tangent fixed by the carried benchmark
residual,

```math
\Lambda_1 \approx 0.279605891931464.
```

So the exact tangent defect is

```math
\boxed{\delta_1 \approx -0.279605891931464.}
```

At the physical electron point this becomes a very small actual outgoing defect,

```math
\boxed{
\Delta_Q^{(e)}
=
-\frac{\Lambda_1 f}{1+\Lambda_1 f}
\approx -3.24631584151692\times 10^{-4},
}
```

with

```math
\chi_Q^{(e)}\approx 0.999675368415848,
\qquad
N_Q^{(e)}\approx 1.00032473700404.
```

So the quartic closure does **not** require a large renormalization of the
canonical outgoing branch. It requires only a tiny finite-`f` drift away from it.

---

## Step 29A — Exact bridge from outgoing defect to common conservative defect

On the natural source-map branch, Step 28 gave

```math
N_Q = \frac{1}{\chi_Q} = \frac{1}{1+\Delta_Q}.
```

Therefore

```math
\boxed{
N_Q - 1
=
-\frac{\Delta_Q}{1+\Delta_Q}.
}
```

This is the exact scalar bridge between the retarded outgoing defect `\Delta_Q`
and the conservative normalization defect `N_Q-1`.

---

## Step 29B — Why the outgoing defect must start at `O(f)`

Step 2 already showed that the missing common anomaly layer acts only on the
already-existing local transport residue, not on the whole anomaly law. The
general leading form was

```math
\Delta\!\left(\frac g2\right)
=
c_{3,\rm total}\,\Lambda_{\rm common}(f)\,f^3
+
O(f^5),
```

and to close a quartic gap one needs

```math
\Lambda_{\rm common}(f)=O(f).
```

If we now identify the common scalar with the conservative branch defect,

```math
\Lambda_{\rm common}(f)=N_Q(f)-1,
```

then the outgoing defect must also begin linearly:

```math
\boxed{
\Delta_Q(f)=\delta_1 f + O(f^2).
}
```

This is the key compatibility point.

- In the strict point-particle / weak-coupling limit `f\to0`, the outgoing defect
  vanishes and the canonical compact outgoing branch is recovered.
- At finite electron `f`, the same branch can still carry a tiny deformation that
  shows up only in the quartic anomaly layer.

So the gravitational reduced theorem and the electron-anomaly closure are not in
conflict.

---

## Step 29C — Quartic anomaly law in terms of `\delta_1`

Insert

```math
\Delta_Q(f)=\delta_1 f + O(f^2)
```

into

```math
N_Q(f)-1=-\frac{\Delta_Q(f)}{1+\Delta_Q(f)}.
```

Then

```math
N_Q(f)-1 = -\delta_1 f + O(f^2).
```

So the common anomaly correction becomes

```math
\Delta\!\left(\frac g2\right)
=
-c_{3,\rm total}\,\delta_1\,f^4
+
O(f^5).
```

Matching this to the Step-2 carried benchmark

```math
\Delta\!\left(\frac g2\right)
=
c_{3,\rm total}\,\Lambda_1\,f^4
+
O(f^5)
```

forces

```math
\boxed{
\delta_1 = -\Lambda_1.
}
```

Numerically,

```math
\boxed{
\delta_1 \approx -0.279605891931464.
}
```

So the quartic anomaly sliver is exactly the first tangent of the outgoing
normalization defect.

---

## Step 29D — Actual electron-point defect

The linear tangent coefficient `\delta_1` is not the physical deformation itself.
The actual electron-point defect is of order `f` smaller:

```math
\Delta_Q^{(e)}
=
-\frac{\Lambda_1 f}{1+\Lambda_1 f}.
```

Using the carried benchmark value

```math
f = 0.001161409732093,
```

gives

```math
\boxed{
\Lambda_1 f \approx 3.24737004039746\times 10^{-4},
}
```

and therefore

```math
\boxed{
\Delta_Q^{(e)}
\approx
-3.24631584151692\times 10^{-4}.
}
```

The corresponding branch values are

```math
\boxed{
\chi_Q^{(e)} = 1+\Delta_Q^{(e)} \approx 0.999675368415848,
}
```

```math
\boxed{
N_Q^{(e)} = \frac{1}{1+\Delta_Q^{(e)}} \approx 1.00032473700404.
}
```

So the needed defect is only a few parts in `10^4`.

---

## Step 29E — Exact improved anomaly law

The improved reduced anomaly law can now be written directly in terms of the
outgoing defect:

```math
\boxed{
\frac{g_{\rm imp}(f)}{2}
=
\frac{g_{\rm loc}(f)}{2}
-
c_{3,\rm total}\,
\frac{\Delta_Q(f)}{1+\Delta_Q(f)}\,f^3
+
O(f^5).
}
```

Choosing the electron-point defect above reproduces the carried exact residual:

```math
g_e-g_{\rm loc}\approx 2.27204390584705\times 10^{-12}.
```

So the quartic gap is exactly saturated by this tiny finite-`f` outgoing
normalization drift.

---

## Main result of the step

The remaining quartic g-2 sliver can be encoded as a tiny outgoing-normalization
defect that vanishes in the strict point-particle limit:

```math
\boxed{
\Delta_Q(f)= -\Lambda_1 f + O(f^2)
}
```

at tangent level, or more explicitly at the electron point,

```math
\boxed{
\Delta_Q^{(e)}
=
-\frac{\Lambda_1 f}{1+\Lambda_1 f}
\approx -3.24631584151692\times10^{-4}.
}
```

This is a strong simplification:

- the gravitational reduced theorem still wants the canonical compact outgoing
  branch as `f\to0`,
- while the finite electron anomaly needs only a very small `O(f)` departure from
  that branch.

So the next clean derivation is to translate this outgoing-defect tangent into the
explicit isotropic DtN branch-selection variables from Step 28.
# Step 30 — The quartic g-2 sliver is an explicit isotropic DtN branch-selection target

## Goal

Step 29 rewrote the remaining quartic anomaly sliver as a tiny outgoing-normalization
drift,

```math
\Delta_Q(f)= -\Lambda_1 f + O(f^2),
```

at tangent level, with

```math
\Lambda_1 \approx 0.279605891931464.
```

Step 28 had already shown that the most general first-order isotropic DtN branch
selection is

```math
\chi_Q
=
1+\varepsilon\left(5b+\frac{a_0}{3}+9a_5\right)+O(\varepsilon^2).
```

The next honest move is therefore obvious:

> **identify the anomaly small parameter with the DtN branch-selection parameter and solve for the required isotropic deformation combination.**

That is what this step does.

The main result is the exact tangent constraint

```math
\boxed{
5b+\frac{a_0}{3}+9a_5 = -\Lambda_1
\approx -0.279605891931464.
}
```

So the quartic electron-anomaly closure is no longer a vague “next common PDE
layer.” It is one explicit isotropic DtN branch-selection target.

---

## Step 30A — Identifying the small parameter

From Step 28,

```math
\chi_Q = 1+\varepsilon\left(5b+\frac{a_0}{3}+9a_5\right)+O(\varepsilon^2).
```

For the anomaly problem the only small reduced parameter already in play is

```math
f=\frac{\alpha_{\rm fs}}{2\pi}.
```

So the natural identification is

```math
\varepsilon=f.
```

Then

```math
\Delta_Q(f):=\chi_Q-1
=
f\left(5b+\frac{a_0}{3}+9a_5\right)+O(f^2).
```

But Step 29 fixed the same outgoing-defect tangent to be

```math
\Delta_Q(f)=-\Lambda_1 f+O(f^2).
```

Matching the two immediately gives

```math
\boxed{
5b+\frac{a_0}{3}+9a_5 = -\Lambda_1.
}
```

This is the exact isotropic DtN branch-selection law required by the quartic
electron-anomaly sliver.

---

## Step 30B — Direct quartic anomaly law in DtN variables

Since Step 29 also showed

```math
\Delta\!\left(\frac g2\right)
=
-c_{3,\rm total}\,\delta_1\,f^4+O(f^5),
\qquad
\delta_1=-\Lambda_1,
```

the quartic anomaly correction can now be written directly as

```math
\boxed{
\Delta\!\left(\frac g2\right)
=
-c_{3,\rm total}
\left(5b+\frac{a_0}{3}+9a_5\right)f^4
+
O(f^5).
}
```

So the missing `O(f^4)` layer is exactly the first isotropic DtN deformation of
the outgoing grouped-`P2` branch.

---

## Step 30C — Three pure isotropic realizations

The constraint

```math
5b+\frac{a_0}{3}+9a_5 = -\Lambda_1
```

does not yet tell us which microscopic branch the PDE picks.
But it already gives three clean pure realizations.

### Pure argument deformation

Set `a_0=a_5=0`. Then

```math
\boxed{
b = -\frac{\Lambda_1}{5}
\approx -0.0559211783862928.
}
```

### Pure static additive slot

Set `b=a_5=0`. Then

```math
\boxed{
a_0 = -3\Lambda_1
\approx -0.838817675794392.
}
```

### Pure odd core outlet

Set `b=a_0=0`. Then

```math
\boxed{
a_5 = -\frac{\Lambda_1}{9}
\approx -0.0310673213257182.
}
```

These are not all equally plausible physically, but they are exact first-order
bookkeeping realizations of the same anomaly target.

---

## Step 30D — Minimum-norm bookkeeping branch

If one ignores the physical dimensions of `(b,a_0,a_5)` and asks only for the
smallest Euclidean-norm tangent triple satisfying the constraint, the answer is

```math
\boxed{
(b,a_0,a_5)_{\min}
=
-\Lambda_1
\frac{(5,\;1/3,\;9)}{5^2+(1/3)^2+9^2}.
}
```

Numerically,

```math
\boxed{
(b,a_0,a_5)_{\min}
\approx
(-0.0131751467402261,\;
 -0.000878343116015070,\;
 -0.0237152641324069).
}
```

This is best read as a bookkeeping baseline, not as a real physical preference,
because the three deformation coordinates carry different meanings in the DtN
problem.

---

## Step 30E — Actual electron-point size of the deformation

The tangent combination itself is `O(1)`, but the actual electron-point
deformation is smaller by one factor of `f`:

```math
f\left(5b+\frac{a_0}{3}+9a_5\right)
=
-\Lambda_1 f
\approx -3.24737004039746\times 10^{-4}.
```

That is the linearized electron-point outgoing defect.
The exact Step-29 branch value was

```math
\Delta_Q^{(e)}
=
-\frac{\Lambda_1 f}{1+\Lambda_1 f}
\approx
-3.24631584151692\times 10^{-4},
```

so the linearized DtN law is already accurate to the expected next neglected order.

---

## Main result of the step

The quartic g-2 sliver is now pinned to one explicit isotropic DtN branch-selection
constraint:

```math
\boxed{
5b+\frac{a_0}{3}+9a_5 = -\Lambda_1
\approx -0.279605891931464.
}
```

Equivalently,

```math
\boxed{
\Delta\!\left(\frac g2\right)
=
-c_{3,\rm total}
\left(5b+\frac{a_0}{3}+9a_5\right)f^4
+O(f^5).
}
```

So the next PDE-facing question is no longer broad at all:

> **which actual moving-throat isotropic DtN branch produces the required combination `5b + a_0/3 + 9a_5`?**

That is the cleanest continuation point after the outgoing-defect bridge.
# Step 31 — The quartic g-2 sliver becomes an exact finite-`f` isotropic DtN surface

## Goal

Step 30 fixed only the **tangent** isotropic DtN branch-selection law

```math
5b+\frac{a_0}{3}+9a_5=-\Lambda_1.
```

The moving-throat notes already give the corresponding **exact** isotropic DtN
deformation law, after the canonical even `l=2` fingerprint is enforced,

```math
\chi_Q
=
\frac{3\big(S\beta^5+9\Sigma_5\big)}{3S-\Sigma_0}.
```

So the next honest move is to stop working only at tangent level and insert the
exact Step-29 electron-point target

```math
\chi_e=\frac{1}{1+\Lambda_1 f}.
```

That is what this step does.

The main result is the exact finite-`f` isotropic anomaly surface

```math
\boxed{
3\big(S\beta^5+9\Sigma_5\big)(1+\Lambda_1 f)=3S-\Sigma_0.
}
```

Equivalently,

```math
\boxed{
3S\!\left[(1+\Lambda_1 f)\beta^5-1\right]
+\Sigma_0
+27(1+\Lambda_1 f)\Sigma_5
=0.
}
```

So the quartic g-2 sliver is no longer only a first-order branch-selection
constraint. It is one exact finite-`f` isotropic DtN surface.

---

## Step 31A — Exact electron-point target

From Step 29,

```math
\Delta_Q^{(e)}
=
-\frac{\Lambda_1 f}{1+\Lambda_1 f},
\qquad
\chi_e = 1+\Delta_Q^{(e)} = \frac{1}{1+\Lambda_1 f}.
```

Using the carried benchmark numbers,

```math
\Lambda_1 \approx 0.279605891931464,
\qquad
f \approx 0.001161409732093,
```

one has

```math
\boxed{
\Lambda_1 f \approx 3.24737004039746\times 10^{-4},
}
```

and therefore

```math
\boxed{
\chi_e \approx 0.999675368415848.
}
```

So the electron-point branch is only a very small finite deformation away from
the canonical outgoing value `\chi_Q=1`.

---

## Step 31B — Exact isotropic DtN anomaly surface

The Stage-90/91 deformation law is

```math
\chi_Q
=
\frac{3\big(S\beta^5+9\Sigma_5\big)}{3S-\Sigma_0},
```

where

- `S` is the overall mouth normalization,
- `\beta` is the outgoing-argument deformation,
- `\Sigma_0` is the static isotropic throat-core shift,
- `\Sigma_5` is the extra odd isotropic `l=2` core outlet.

Matching `\chi_Q` to the exact electron target `\chi_e` gives

```math
\boxed{
3\big(S\beta^5+9\Sigma_5\big)(1+\Lambda_1 f)=3S-\Sigma_0.
}
```

Solving for the odd core slot gives the exact branch-selection surface

```math
\boxed{
\Sigma_5
=
\frac{3S-\Sigma_0}{27(1+\Lambda_1 f)}
-\frac{S\beta^5}{9}.
}
```

This is the finite-`f` replacement for the linear Step-30 tangent relation.

---

## Step 31C — Canonical-preservation surface as the zero-anomaly limit

If one instead asks for exact preservation of the canonical outgoing branch,
`\chi_Q=1`, then the same algebra gives

```math
\boxed{
\Sigma_5
=
\frac{S(1-\beta^5)}{9}
-\frac{\Sigma_0}{27}.
}
```

So the exact anomaly surface is a direct deformation of the exact
canonical-preservation surface, obtained simply by replacing `1` by
`1/(1+\Lambda_1 f)`.

That is a good consistency check: at `f\to 0`, the g-2 branch-selection law
collapses back to the canonical outgoing-preservation law.

---

## Step 31D — Three exact pure isotropic realizations

The exact finite-`f` surface already admits three clean pure realizations.

### Pure compensated argument deformation

Set `\Sigma_0=\Sigma_5=0`. Then

```math
\boxed{
\beta = \chi_e^{1/5} = (1+\Lambda_1 f)^{-1/5}.
}
```

Numerically,

```math
\boxed{
\beta_e \approx 0.999935065250674.
}
```

So the required argument deformation is tiny.

### Pure static additive core

Set `\beta=1` and `\Sigma_5=0`. Then

```math
\boxed{
\Sigma_0 = -3S\,\Lambda_1 f.
}
```

At unit scale `S=1`,

```math
\boxed{
\Sigma_0^{(e)} \approx -9.74211012119238\times 10^{-4}.
}
```

### Pure odd `l=2` core outlet

Set `\beta=1` and `\Sigma_0=0`. Then

```math
\boxed{
\Sigma_5 = -\frac{S\,\Lambda_1 f}{9(1+\Lambda_1 f)}.
}
```

At unit scale `S=1`,

```math
\boxed{
\Sigma_5^{(e)} \approx -3.60701760168546\times 10^{-5}.
}
```

So the exact electron-point outgoing defect can be carried by a very small
odd-core outlet as well.

---

## Step 31E — Exact surface reduces to the Step-30 tangent law

Write

```math
S=1+\varepsilon s,
\qquad
\beta=1+\varepsilon b,
\qquad
\Sigma_0=\varepsilon a_0,
\qquad
\Sigma_5=\varepsilon a_5,
```

and expand the exact surface together with

```math
\chi_e=\frac{1}{1+\varepsilon\Lambda_1}.
```

The first-order condition is exactly

```math
\boxed{
5b+\frac{a_0}{3}+9a_5+\Lambda_1=0,
}
```

i.e.

```math
\boxed{
5b+\frac{a_0}{3}+9a_5=-\Lambda_1.
}
```

So Step 30 is the strict tangent reduction of the present exact finite-`f`
branch-selection law.

---

## Main result of the step

The quartic g-2 sliver is now encoded by the exact isotropic DtN surface

```math
\boxed{
3\big(S\beta^5+9\Sigma_5\big)(1+\Lambda_1 f)=3S-\Sigma_0.
}
```

with exact pure-branch representatives

```math
\boxed{
\beta = (1+\Lambda_1 f)^{-1/5},
\qquad
\Sigma_0=-3S\Lambda_1 f,
\qquad
\Sigma_5=-\frac{S\Lambda_1 f}{9(1+\Lambda_1 f)}.
}
```

So the next PDE-facing question becomes sharper again:

> **which explicit moving-throat outlet model can realize this exact finite-`f` surface while preserving the already-fixed conservative even `l=2` branch?**
# Step 32 — The first explicit outlet audit leaves only a compensated Robin–mixed branch alive

## Goal

Step 31 turned the quartic g-2 sliver into an exact finite-`f` isotropic DtN
surface. The next honest question is no longer abstract:

> **which explicit low-frequency moving-throat outlet classes can realize that
> surface without spoiling the already-fixed conservative even `l=2`
> fingerprint?**

The moving-throat outlet notes already tested three explicit isotropic classes:

1. a pure geometric Robin core,
2. a standalone mixed `A_w/F_{\mu w}` side-channel pole,
3. a hybrid Robin–mixed outlet.

This step translates those outlet classes directly into the g-2 branch-selection
language.

The main result is sharp:

- a pure Robin core can match the electron outgoing defect exactly, but it
  distorts the even branch and therefore cannot be the whole story;
- a standalone mixed pole cannot preserve the canonical even branch unless it is
  absent;
- the first serious surviving candidate is the **compensated Robin–mixed
  outlet**.

---

## Step 32A — Pure geometric Robin core

The raw isotropic Robin outlet has the exact normalization factor

```math
\boxed{
\chi_Q^{\rm R}=\frac{3}{3-\rho_R},
}
```

where `\rho_R` is the dimensionless isotropic Robin loading.

Matching this directly to the exact electron target

```math
\chi_e=\frac{1}{1+\Lambda_1 f}
```

gives

```math
\boxed{
\rho_R = -3\Lambda_1 f.
}
```

Numerically,

```math
\boxed{
\rho_R^{(e)}
\approx
-9.74211012119238\times 10^{-4}.
}
```

So a pure Robin core can indeed carry the required *size* of the outgoing defect.

But the moving-throat notes also show that the raw Robin outlet changes the
canonical even `l=2` fingerprint at the same time. So it cannot be the full
answer if the conservative grouped-`P2` even branch is already fixed.

That makes the pure Robin result useful only as a scale diagnostic, not as a
finished outlet law.

---

## Step 32B — Standalone mixed side-channel pole

The raw mixed side-channel model is

```math
\Lambda_2^{\rm mix}(z)
=
\Lambda_2^{\rm out}(z)
-\frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5}.
```

For this branch, the canonical even conditions reduce to

```math
-\frac{L_2}{L_0}=\frac19,
\qquad
\frac{L_2^2}{L_0^2}-\frac{L_4}{L_0}=\frac{4}{81}.
```

Solving them gives first

```math
\kappa_W=-\frac19,
```

and then, after substitution,

```math
\boxed{
\sigma_W=0.
}
```

So a standalone passive mixed pole of this type cannot sit on the already-fixed
canonical even branch unless it vanishes.

That is a strong exclusion result:

> the mixed sector may still matter, but not as a naive isolated Schur-complement
> outlet of this simplest type.

---

## Step 32C — Hybrid Robin–mixed outlet

The hybrid outlet is

```math
\Lambda_2^{\rm hyb}(z)
=
\Lambda_2^{\rm out}(z)
+\rho_R
-\frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5}.
```

Imposing the canonical even `l=2` fingerprint yields exactly two branches:

```math
\boxed{
\rho_R=\sigma_W,\qquad \kappa_W=0,
}
```

or

```math
\boxed{
\rho_R=4\sigma_W,\qquad \kappa_W=\frac13.
}
```

The first is only a trivial cancellation branch.

The second is the first **nontrivial compensated branch**.
On that branch the exact normalization factor is

```math
\boxed{
\chi_Q^{\rm hyb}
=
\frac{1-9\sigma_W\gamma_W}{1-\sigma_W}.
}
```

For the canonical outgoing branch itself one recovers

```math
\boxed{
\gamma_W=\frac19
\quad\Longrightarrow\quad
\chi_Q^{\rm hyb}=1.
}
```

So the compensated hybrid outlet reproduces the same harmless pure-scale class
already identified in the Stage-90/91 robustness audit when its odd coefficient
is set to the canonical value.

---

## Step 32D — Exact electron-anomaly family on the compensated branch

Now replace the canonical target `\chi_Q=1` by the exact electron target
`\chi_e=1/(1+\Lambda_1 f)`.

Solving

```math
\frac{1-9\sigma_W\gamma_W}{1-\sigma_W}
=
\frac{1}{1+\Lambda_1 f}
```

gives the exact anomaly family

```math
\boxed{
\gamma_W
=
\frac{\sigma_W+\Lambda_1 f}
{9\sigma_W(1+\Lambda_1 f)}.
}
```

Equivalently,

```math
\boxed{
\gamma_W
=
\frac19
+
\frac{\Lambda_1 f\,(1-\sigma_W)}
{9\sigma_W(1+\Lambda_1 f)}.
}
```

So on the compensated hybrid branch, the electron anomaly is carried by a small
shift of the odd mixed outlet above its canonical value `1/9`.

That is exactly the kind of result we wanted:

- the even branch remains fixed,
- the outlet deformation is not broad or arbitrary,
- and the anomaly rides only in one controlled odd coefficient.

For a representative moderate loading `\sigma_W=1/2`, this gives

```math
\boxed{
\gamma_W \approx 0.111147181287128,
}
```

which is only slightly above

```math
\frac19 \approx 0.111111111111111.
```

So the required odd detuning is tiny, just as the anomaly itself is tiny.

---

## Step 32E — Interpretation of the outlet audit

This is the first really constraining answer after Step 31.

### What is ruled out

- **Pure Robin loading alone** is too blunt: it can carry the defect size but it
  distorts the even branch.
- **A naive standalone mixed pole** is too rigid: it cannot preserve the
  canonical even branch unless it disappears.

### What survives

- The first serious surviving outlet class is the **compensated Robin–mixed
  branch**
  ```math
  \rho_R=4\sigma_W,\qquad \kappa_W=\frac13.
  ```
- On that branch the anomaly is carried entirely by the odd side-channel
  coefficient
  ```math
  \gamma_W=\frac{\sigma_W+\Lambda_1 f}{9\sigma_W(1+\Lambda_1 f)}.
  ```

So the outlet problem is now much narrower than it was before.

It is no longer

> “some deformed isotropic DtN branch.”

It is now

> **an explicit compensated Robin–mixed outlet family with one tiny odd
> coefficient shift left to determine.**

---

## Main result of the step

The first explicit outlet audit leaves only one serious exact candidate class
alive for the g-2 quartic sliver:

```math
\boxed{
\rho_R=4\sigma_W,\qquad \kappa_W=\frac13,\qquad
\gamma_W=\frac{\sigma_W+\Lambda_1 f}{9\sigma_W(1+\Lambda_1 f)}.
}
```

So the next natural move is even sharper:

> **push this compensated outlet down into an actual throat-core model and see
> whether the required odd shift can be written directly in core couplings rather
> than as reduced DtN coefficients.**
# Step 33 — The surviving outlet class becomes a concrete balanced throat core with one tiny odd detuning

## Goal

Step 32 reduced the outlet problem to one surviving explicit class:

```math
\rho_R=4\sigma_W,\qquad
\kappa_W=\frac13,\qquad
\gamma_W=\frac{\sigma_W+\Lambda_1 f}{9\sigma_W(1+\Lambda_1 f)}.
```

That is still a reduced DtN description. The next natural move is to push it into
the first concrete throat-core model already worked out in the moving-throat
notes.

Those notes introduce a two-channel isotropic core with

- a static shell compliance coordinate,
- a mixed `A_w/F_{\mu w}` side-channel,
- and a static/mixed hybridization.

This step shows that the quartic g-2 sliver can be written directly in that core
language.

The main result is:

- the **same exact core-balance surface** that preserves the canonical even
  branch still survives,
- the **same auxiliary D/N mixed-tube length relation** still survives,
- and the anomaly is carried only by a tiny **odd detuning** of the balanced
  mixed side-channel.

---

## Step 33A — Concrete two-channel core data

For the explicit isotropic core, the reduced outlet coefficients are

```math
\boxed{
\rho_c=\frac{g_s^2}{K_s},
}
```

```math
\boxed{
\sigma_c=
\frac{(K_s g_q-\lambda g_s)^2}
{K_s^2K_q(1+r_c)},
\qquad
r_c=\frac{\lambda^2}{K_sK_q},
}
```

```math
\boxed{
\kappa_c=\frac{\kappa_0}{1+r_c},
\qquad
\gamma_c=\frac{\gamma_0}{1+r_c}.
}
```

So the reduced Robin–mixed outlet is no longer described by four unrelated
numbers. It is controlled by the concrete throat-core parameters

- `K_s`,
- `K_q`,
- `\lambda`,
- `g_s`,
- `g_q`,
- and the bare mixed low-frequency pair `(\kappa_0,\gamma_0)`.

That is already a major narrowing.

---

## Step 33B — Exact balance surface preserving the canonical even branch

The moving-throat notes already proved that the compensated canonical branch is
not a numerical accident. In the concrete core variables it is the exact
codimension-one balance law

```math
\boxed{
g_s^2\bigl(K_sK_q+\lambda^2\bigr)
=
4\bigl(K_s g_q-\lambda g_s\bigr)^2.
}
```

Solving for the mixed mouth coupling gives the two exact branches

```math
\boxed{
g_q=
\frac{g_s}{2K_s}
\left(
2\lambda \pm \sqrt{K_sK_q+\lambda^2}
\right).
}
```

On either branch one gets

```math
\boxed{
\sigma_c=\frac{g_s^2}{4K_s}.
}
```

So the same balance condition that preserved the canonical outgoing branch in the
reduced DtN language is realized exactly by a concrete two-channel throat core.

The even-preserving condition for the bare mixed channel is still

```math
\boxed{
\kappa_0=\frac{1+r_c}{3},
}
```

which means the effective mixed coefficient remains

```math
\boxed{
\kappa_c=\frac13.
}
```

So the g-2 anomaly does **not** ask us to reopen the even branch. The same
even-preserving core geometry survives intact.

---

## Step 33C — Exact electron-anomaly law on the balanced core

On the balanced-even core, the effective normalization law is still the hybrid
one,

```math
\chi_Q^{\rm core}
=
\frac{1-9\sigma_c\gamma_c}{1-\sigma_c}.
```

Matching this to the exact electron target

```math
\chi_e=\frac{1}{1+\Lambda_1 f}
```

gives the required **effective** odd coefficient

```math
\boxed{
\gamma_c
=
\frac{\sigma_c+\Lambda_1 f}
{9\sigma_c(1+\Lambda_1 f)}.
}
```

Equivalently,

```math
\boxed{
\gamma_c
=
\frac19
+
\frac{\Lambda_1 f\,(1-\sigma_c)}
{9\sigma_c(1+\Lambda_1 f)}.
}
```

So the anomaly is carried only by a small shift above the canonical value
`1/9`, provided `0<\sigma_c<1`.

Because

```math
\gamma_c=\frac{\gamma_0}{1+r_c},
```

the corresponding **bare** mixed odd coefficient is

```math
\boxed{
\gamma_0
=
\frac{1+r_c}{9}
+
\frac{(1+r_c)\Lambda_1 f\,(1-\sigma_c)}
{9\sigma_c(1+\Lambda_1 f)}.
}
```

That is the exact concrete-core version of the g-2 outlet target.

So the quartic anomaly sliver is no longer a generic branch deformation. It is a
tiny odd detuning of an otherwise canonical balanced core.

---

## Step 33D — The auxiliary D/N mixed-tube geometry is unchanged

If the bare mixed side-channel is realized by the first D/N half-wave on an
auxiliary tube of length `L_W`, then

```math
\boxed{
\kappa_0=\frac{4L_W^2}{\pi^2 a^2}.
}
```

Imposing the same even-preserving condition

```math
\kappa_0=\frac{1+r_c}{3}
```

gives exactly the same auxiliary-tube length as on the canonical compensated
branch:

```math
\boxed{
L_W
=
\frac{\pi a}{2}\sqrt{\frac{1+r_c}{3}}.
}
```

So the anomaly does **not** force a new even geometry for the side-channel tube.
It leaves the D/N half-wave length relation untouched and only detunes the odd
side-channel coefficient.

That is a particularly clean result.

---

## Step 33E — What changed conceptually

After Step 32, the best surviving outlet class was still described by reduced DtN
coefficients. After this step, it becomes a genuine microscopic core statement:

- balance the shell and mixed core on the exact surface
  ```math
  g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2,
  ```
- keep the same even-preserving auxiliary D/N geometry
  ```math
  L_W=\frac{\pi a}{2}\sqrt{\frac{1+r_c}{3}},
  ```
- and shift only the odd mixed outlet from its canonical value by the small
  amount above.

So the anomaly no longer looks like a request for a new structural mechanism.
It looks like a very small odd detuning of a core structure that the moving-throat
notes already know how to realize exactly on the canonical branch.

---

## Main result of the step

The surviving compensated outlet class has now been pushed into a concrete
throat-core model, and the electron anomaly requires only

```math
\boxed{
\gamma_c
=
\frac{\sigma_c+\Lambda_1 f}
{9\sigma_c(1+\Lambda_1 f)},
\qquad
\gamma_0=(1+r_c)\gamma_c,
}
```

while preserving the exact core-balance and D/N side-channel geometry:

```math
\boxed{
g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2,
\qquad
L_W=\frac{\pi a}{2}\sqrt{\frac{1+r_c}{3}}.
}
```

So the next PDE-facing question is now almost as sharp as it can be:

> **does the actual moving-throat core land on this balanced-even structure, and
> if so, does its odd mixed side-channel coefficient acquire precisely the tiny
> detuning above?**
# Step 34 — The balanced outlet collapses to a one-parameter parent compensation family

## Goal

Step 33 pushed the surviving compensated outlet into a concrete two-channel
throat core. The next natural move is to eliminate the remaining reduced-core
bookkeeping and rewrite the whole balanced branch directly in terms of the
**parent overlap ratios** already isolated in the moving-throat notes.

That matters because it tells us whether the anomaly is asking for a new core
mechanism, or only for a small detuning of a branch that is already naturally
present in the parent-action variables.

---

## Step 34A — The exact core-balance theorem in parent ratios

Introduce the two normalized parent ratios

```math
\boxed{
\mathfrak r:=\frac{\lambda}{\sqrt{K_sK_q}},
\qquad
\mathfrak g:=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}}.
}
```

Then the exact concrete core-balance law from Step 33,

```math
g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2,
```

collapses to the one-line parent relation

```math
\boxed{
1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2.
}
```

So the compensated-even outlet is no longer a five-parameter tuning problem.
It is a **one-parameter family**: once the normalized static/mixed hybridization
`\mathfrak r` is fixed by the throat background, the required normalized mouth
coupling `\mathfrak g` is fixed exactly.

The two exact branches are

```math
\boxed{
\mathfrak g
=
\mathfrak r
\pm
\frac12\sqrt{1+\mathfrak r^2}.
}
```

That is the parent-action version of the Stage-98 compensation surface.

---

## Step 34B — The side-tube geometry is fixed by the same ratio

On the balanced-even branch the bare mixed coefficient obeys

```math
\kappa_0=\frac{1+\mathfrak r^2}{3}.
```

If the mixed side-channel is realized by the first D/N half-wave on an auxiliary
 tube of length `L_W`, then

```math
\kappa_0=\frac{4L_W^2}{\pi^2 a^2}.
```

So the exact D/N tube-selection law becomes

```math
\boxed{
L_W=
\frac{\pi a}{2}\sqrt{\frac{1+\mathfrak r^2}{3}}.
}
```

So the same normalized hybridization `\mathfrak r` fixes both:

- the required mouth-coupling balance, and
- the auxiliary D/N geometry.

That is a very strong collapse of the surviving outlet freedom.

---

## Step 34C — The loading share simplifies to `\rho_c/4`

Recall the concrete core quantities

```math
\rho_c=\frac{g_s^2}{K_s},
\qquad
\sigma_c=
\frac{(K_s g_q-\lambda g_s)^2}{K_s^2K_q(1+r_c)},
\qquad
r_c=\frac{\lambda^2}{K_sK_q}=\mathfrak r^2.
```

In parent-ratio language,

```math
\boxed{
\sigma_c
=
\rho_c\,
\frac{(\mathfrak g-\mathfrak r)^2}{1+\mathfrak r^2}.
}
```

So on either exact compensation branch,

```math
\boxed{
\sigma_c=\frac{\rho_c}{4}.
}
```

That means the balanced branch does not leave both `\rho_c` and `\sigma_c`
independent. The static loading share is fixed to one quarter of the shell-side
normalization `\rho_c`.

---

## Step 34D — The electron anomaly becomes a tiny odd detuning of the parent family

Let

```math
x:=\Lambda_1 f.
```

The exact electron target from the earlier steps is

```math
\chi_e=\frac{1}{1+x}.
```

On the balanced hybrid outlet,

```math
\chi_Q=\frac{1-9\sigma_c\gamma_c}{1-\sigma_c}.
```

Substituting `\sigma_c=\rho_c/4` and matching to `\chi_e` gives

```math
\boxed{
\gamma_c
=
\frac{\rho_c+4x}{9\rho_c(1+x)}.
}
```

Because `\gamma_0=(1+\mathfrak r^2)\gamma_c`, the corresponding bare mixed odd
coefficient is

```math
\boxed{
\gamma_0
=
\frac{(1+\mathfrak r^2)(\rho_c+4x)}{9\rho_c(1+x)}.
}
```

So the anomaly is now even more sharply localized than in Step 33:

- the compensated-even branch itself is fixed by `\mathfrak r`,
- the static loading share is `\rho_c/4`,
- and the anomaly only asks for a small odd detuning above the canonical values
  ```math
  \gamma_c=\frac19,
  \qquad
  \gamma_0=\frac{1+\mathfrak r^2}{9}.
  ```

Indeed,

```math
\boxed{
\gamma_c-\frac19
=
\frac{x(4-\rho_c)}{9\rho_c(1+x)}.
}
```

So at fixed `\rho_c=O(1)` and tiny `x`, the required shift is necessarily tiny.

---

## Main result of the step

The surviving compensated outlet is no longer a free reduced core at all. It has
collapsed to a one-parameter parent family:

```math
\boxed{
1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2,
\qquad
L_W=\frac{\pi a}{2}\sqrt{\frac{1+\mathfrak r^2}{3}},
\qquad
\sigma_c=\frac{\rho_c}{4}.
}
```

And on that family the electron anomaly requires only

```math
\boxed{
\gamma_c=
\frac{\rho_c+4x}{9\rho_c(1+x)},
\qquad
\gamma_0=
\frac{(1+\mathfrak r^2)(\rho_c+4x)}{9\rho_c(1+x)}.
}
```

So the next real question is no longer “is there a parent family that can host the
outlet?” There is. The sharper question is whether the **actual mouth/source law**
selects a natural point on that family without fine tuning.
# Step 35 — Positive mouth sourcing selects the lower branch without fine tuning

## Goal

Step 34 reduced the compensated outlet to a one-parameter parent family labeled by
`\mathfrak r`, with the required mouth-coupling ratio fixed by

```math
\mathfrak g=\mathfrak r\pm\frac12\sqrt{1+\mathfrak r^2}.
```

The next real question is no longer algebraic. It is physical:

> does any natural localized mouth-source law actually land on that family, and if
> so, does it do it without fine tuning?

This step answers that question on the explicit Family-1 branch.

---

## Step 35A — The explicit Family-1 parent ratio and its two compensated branches

For the carried explicit Family-1 geometry,

```math
\boxed{
\mathfrak r_{F1}
=
\sqrt{\frac{12}{\pi^2}\left(\frac{37}{20}\right)^2-1}
\approx 1.77799353547498.
}
```

So the two compensated parent branches are

```math
\boxed{
\mathfrak g_-^{F1}
=
\mathfrak r_{F1}-\frac12\sqrt{1+\mathfrak r_{F1}^2}
\approx 0.758035078944663,
}
```

```math
\boxed{
\mathfrak g_+^{F1}
=
\mathfrak r_{F1}+\frac12\sqrt{1+\mathfrak r_{F1}^2}
\approx 2.79795199200529.
}
```

So the purely algebraic family still contains two branches. The physical mouth
source law has to decide which one is actually admissible.

---

## Step 35B — Positive mouth-source theorem

Let `\sigma(z)` be any nonnegative normalized axial source profile on the first
D/N throat interval,

```math
\sigma(z)\ge 0,
\qquad
\int_0^L \sigma(z)\,dz = 1,
\qquad
z\in[0,L].
```

The normalized mouth-bias factor is then

```math
\boxed{
\mathfrak g[\sigma]
=
\int_0^L \sigma(z)
\cos\!\left(\frac{\pi z}{2L}\right)
\,dz.
}
```

Because

```math
0\le \cos\!\left(\frac{\pi z}{2L}\right)\le 1
\qquad \text{on } [0,L],
```

any positive normalized source law satisfies

```math
\boxed{
0\le \mathfrak g[\sigma]\le 1.
}
```

That immediately kills the upper Family-1 branch, because

```math
\mathfrak g_+^{F1}>1.
```

So the upper compensated branch is physically impossible for any positive localized
mouth source.

By contrast,

```math
0<\mathfrak g_-^{F1}<1,
```

so the lower compensated branch is the **unique physically admissible canonical
candidate** under positive mouth sourcing.

---

## Step 35C — The natural self-matched derivative profile is already very close

The most natural positive axial source on the first D/N interval is the normalized
first-derivative profile itself,

```math
\boxed{
\sigma_{\rm match}(z)=k\cos(kz),
\qquad
k=\frac{\pi}{2L}.
}
```

It is normalized, and its mouth-bias factor is exactly

```math
\boxed{
\mathfrak g_{\rm match}=\frac{\pi}{4}\approx 0.785398163397448.
}
```

Comparing this to the exact Family-1 lower branch,

```math
\mathfrak g_-^{F1}\approx 0.758035078944663,
```

gives only a small traction mismatch. Since traction scales as
`\mathfrak g^{-1}`,

```math
\boxed{
\frac{\mathcal T_m^{(-)}}{\mathcal T_m^{\rm match}}
=
\frac{\pi/4}{\mathfrak g_-^{F1}}
\approx 1.036097385480999.
}
```

So the exact lower compensated Family-1 branch is only a **3.61% traction
enhancement** away from the most natural self-matched derivative source law.

That is already a strong “not fine tuned” signal.

---

## Step 35D — An explicit convex positive family hits the lower branch with only an 18.4% admixture

Introduce the exact convex positive family

```math
\boxed{
\sigma_\xi(z)
=
(1-\xi)\,k\cos(kz)+\xi\,\frac1L,
\qquad 0\le \xi\le 1.
}
```

It stays positive and normalized, and its bias factor is

```math
\boxed{
\mathfrak g_\xi
=
(1-\xi)\frac{\pi}{4}
+
\xi\frac{2}{\pi}.
}
```

Solving `\mathfrak g_\xi=\mathfrak g_-^{F1}` gives the exact interpolation point

```math
\boxed{
\xi_*
=
\frac{\frac{\pi}{4}-\mathfrak g_-^{F1}}{\frac{\pi}{4}-\frac{2}{\pi}}
\approx 0.183918405511538.
}
```

So only an **18.4% admixture** of the fully washed positive profile `1/L` into the
natural derivative profile reaches the exact lower compensated branch.

Again, that is a small deformation of a very natural source family, not a delicate
coefficient fit.

---

## Step 35E — The regular canonical point occurs at finite bias

For the explicit exponential mouth-layer family, the bias law is

```math
\boxed{
\mathfrak g_\Pi
=
\frac{2\Pi(2\Pi e^\Pi+\pi)}{(4\Pi^2+\pi^2)(e^\Pi-1)}.
}
```

A direct rearrangement gives an exact positivity decomposition for `1-\mathfrak g_\Pi`,
so in particular

```math
0<\mathfrak g_\Pi<1
\qquad \text{for every finite } \Pi>0,
```

and

```math
\lim_{\Pi\to\infty}\mathfrak g_\Pi = 1.
```

So the naive equal-normalized branch `\mathfrak g=1` is only reached in the
**singular point-source limit**, not at any finite regular mouth bias.

By contrast, the lower compensated Family-1 point occurs already at the finite bias

```math
\boxed{
\Pi_*\approx 1.50882951349316,
\qquad
\mathfrak g_\Pi(\Pi_*)=\mathfrak g_-^{F1}.
}
```

So the physically selected lower branch is regular and nearby, not a singular limit.

---

## Main result of the step

Positive mouth sourcing makes the branch selection problem almost trivial.

For the explicit Family-1 branch,

```math
\boxed{
\mathfrak g_-^{F1}\approx 0.758035078944663,
\qquad
\mathfrak g_+^{F1}\approx 2.79795199200529.
}
```

But every positive normalized source law satisfies

```math
0\le \mathfrak g[\sigma]\le 1,
```

so the upper branch is ruled out and the lower branch is the unique admissible
canonical candidate.

Even better, the lower branch is **close to natural source laws**:

- the self-matched derivative profile sits only `3.61%` away in traction,
- the exact lower branch is reached by only an `18.4%` admixture of the uniform
  positive profile into that derivative source,
- and the regular canonical point occurs already at the finite bias
  `\Pi_*\approx 1.50883`.

So the next question is no longer whether a natural source family can reach the
balanced outlet. It can. The sharper question is how the **actual coupled mouth/core
branch** deforms around that regular lower compensated point.
# Step 36 — Exact lower-branch co-transport and the collapse to four irreducible drifts

## Goal

Step 35 showed that the compensated Family-1 outlet is not a delicate fit: the
lower branch is the unique regular positive-source branch, and it sits close to
several very natural mouth-source laws.

The next natural question is then sharper:

> once we sit on that exact lower compensated branch, how many microscopic drifts
> are still genuinely free?

This step answers that by solving the exact lower-branch transport constraints.
The main result is that the branch already co-transports most of the apparent
microscopic freedom, and the old first-order off-family defect vanishes
identically on the exact lower branch.

---

## Step 36A — The two exact lower-branch channel constraints

On the explicit lower compensated branch, the parent ratios

```math
\mathfrak g=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}},
\qquad
\mathfrak r=\frac{\lambda}{\sqrt{K_sK_q}}
```

are fixed to first order. So the two logarithmic imbalance channels must vanish:

```math
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)=0,
\qquad
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right)=0.
```

Using the explicit throat-core branch formulas from the moving-throat notes, the
corresponding drift equations are

```math
0=
\delta\ln\mathcal Z_q
+3\,\delta\ln c_{s,w}
-\delta\ln\rho_w
-\delta\ln\mathcal T_m
-\delta\ln v_{w0}
-2\,\delta\ln a
-2\,\delta\ln L_W,
```

```math
0=
\delta\ln\mathcal Z_q
+2\,\delta\ln c_s
+3\,\delta\ln c_{s,w}
-\delta\ln\rho_w
-2\,\delta\ln v_{w0}
-2\,\delta\ln a
-3\,\delta\ln L_W.
```

So the lower branch is already imposing two exact co-transport conditions before
we solve any further PDE detail.

---

## Step 36B — The D/N geometry removes `L_W` as an independent drift

The exact D/N side-tube law is

```math
L_W=
\frac{\pi a}{2}\sqrt{\frac{1+\mathfrak r_*^2}{3}}.
```

Since `\mathfrak r_*` is fixed on the exact lower branch, differentiation gives

```math
\boxed{
\delta\ln L_W = \delta\ln a.
}
```

So the side-tube length is not an independent linearized variable on the lower
branch. It co-transports with the mouth radius.

---

## Step 36C — Exact drift laws for `v_{w0}` and `\mathcal T_m`

Substituting `\delta\ln L_W=\delta\ln a` into the two branch constraints and
solving gives

```math
\boxed{
\delta\ln v_{w0}
=
\frac12
\Big(
\delta\ln\mathcal Z_q
-\delta\ln\rho_w
+3\,\delta\ln c_{s,w}
+2\,\delta\ln c_s
-5\,\delta\ln a
\Big),
}
```

```math
\boxed{
\delta\ln\mathcal T_m
=
\frac12
\Big(
\delta\ln\mathcal Z_q
-\delta\ln\rho_w
+3\,\delta\ln c_{s,w}
-2\,\delta\ln c_s
-3\,\delta\ln a
\Big).
}
```

So the background mixed flow and the mouth traction are not independent branch
variables at first order. Once the shell/gauge/background drifts are fixed, the
branch already fixes both of them.

---

## Step 36D — Exact product/ratio factorization

The two drift laws split into a particularly clean product/ratio form.
Subtracting and adding them gives

```math
\boxed{
\delta\ln\!\left(\frac{v_{w0}}{\mathcal T_m}\right)
=
2\,\delta\ln c_s - \delta\ln a,
}
```

```math
\boxed{
\delta\ln\!(v_{w0}\mathcal T_m)
=
\delta\ln\mathcal Z_q
+3\,\delta\ln c_{s,w}
-\delta\ln\rho_w
-4\,\delta\ln a.
}
```

This is very informative.

- The **ratio** `v_{w0}/\mathcal T_m` only remembers the wave-speed/geometry side.
- The **product** `v_{w0}\mathcal T_m` only remembers the shell/localization side.

So the branch is already separating the two physical channels for us.

---

## Step 36E — Frozen `n=5` wall-EOS reduction

On the frozen wall GNLS branch,

```math
c_{s,w}^2\propto \rho_w^4,
```

so

```math
\delta\ln c_{s,w}=2\,\delta\ln\rho_w.
```

Then the exact lower-branch drift laws reduce to

```math
\boxed{
\delta\ln v_{w0}
=
\frac12\Big(
\delta\ln\mathcal Z_q
+5\,\delta\ln\rho_w
+2\,\delta\ln c_s
-5\,\delta\ln a
\Big),
}
```

```math
\boxed{
\delta\ln\mathcal T_m
=
\frac12\Big(
\delta\ln\mathcal Z_q
+5\,\delta\ln\rho_w
-2\,\delta\ln c_s
-3\,\delta\ln a
\Big).
}
```

And the product/ratio channels become

```math
\boxed{
\delta\ln\!\left(\frac{v_{w0}}{\mathcal T_m}\right)
=
2\,\delta\ln c_s-\delta\ln a,
}
```

```math
\boxed{
\delta\ln\!(v_{w0}\mathcal T_m)
=
\delta\ln\mathcal Z_q+5\,\delta\ln\rho_w-4\,\delta\ln a.
}
```

So after the exact branch constraints and the frozen `n=5` wall-EOS law are both
used, the apparently large drift space collapses to only four irreducible
microscopic directions:

```math
\boxed{
(\delta\ln\mathcal Z_q,\;\delta\ln\rho_w,\;\delta\ln c_s,\;\delta\ln a).
}
```

Everything else is already co-transported by the exact lower branch.

---

## Step 36F — The old off-family scalar vanishes identically on the exact lower branch

The earlier moving-throat notes reduced the first off-family defect to the scalar

```math
\delta_\perp
=
\mathfrak g_*\,\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)
+
\frac{1}{4\sqrt{1+\mathfrak r_*^2}}
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right).
```

But on the exact lower compensated branch, those two channels are precisely the
constraints we solved above, so

```math
\boxed{
\delta_\perp=0
}
```

identically.

That is the honest reason the old first-order defect disappears on the exact
lower branch: the lower-branch transport laws already enforce its cancellation.

---

## Main result of the step

The compensated Family-1 outlet is now much sharper than it looked even one step
ago.

The exact lower branch already forces

```math
\boxed{
\delta\ln L_W=\delta\ln a,
}
```

and fixes the mixed-flow and traction drifts by

```math
\boxed{
\delta\ln v_{w0}
=
\frac12(
\delta\ln\mathcal Z_q
-\delta\ln\rho_w
+3\,\delta\ln c_{s,w}
+2\,\delta\ln c_s
-5\,\delta\ln a),
}
```

```math
\boxed{
\delta\ln\mathcal T_m
=
\frac12(
\delta\ln\mathcal Z_q
-\delta\ln\rho_w
+3\,\delta\ln c_{s,w}
-2\,\delta\ln c_s
-3\,\delta\ln a).
}
```

After the frozen `n=5` wall-EOS reduction, the entire actual branch-drift problem
collapses to only

```math
\boxed{
(\delta\ln\mathcal Z_q,\;\delta\ln\rho_w,\;\delta\ln c_s,\;\delta\ln a),
}
```

and the old first-order off-family scalar vanishes identically.

So the next question is no longer “is there still a first-order branch defect?”
There is not. The sharper PDE-facing question is:

> what does the true moving-throat solution do to those **four** residual lower-branch variables?
# Step 37 — Exact bundle inversion of the last four irreducible branch drifts

## Goal

Step 36 reduced the exact lower compensated Family-1 branch to only four truly
irreducible microscopic drifts,

```math
\delta\ln\mathcal Z_q,
\qquad
\delta\ln\rho_w,
\qquad
\delta\ln c_s,
\qquad
\delta\ln a.
```

The right next move is not to guess those drifts one by one. The grouped
wall/BdG/Maxwell bundle already exposes four natural observables that determine
all of them algebraically:

```math
\Theta_w,
\qquad
K_s,
\qquad
K_q,
\qquad
P_0=\frac{N_0}{D_0}.
```

This step performs that inversion exactly.

---

## Step 37A — The four exact branch laws

On the explicit Family-1 wall branch,

```math
\Theta_w=25\lambda_\mu^2\rho_w^2,
```

so at fixed `\lambda_\mu`,

```math
\boxed{
\delta\ln\Theta_w = 2\,\delta\ln\rho_w.
}
```

On the healing-locked shell branch,

```math
K_s\propto a^2\rho_w,
```

so

```math
\boxed{
\delta\ln K_s = 2\,\delta\ln a + \delta\ln\rho_w.
}
```

On the exact lower compensated D/N branch,

```math
K_q\propto \mathcal Z_q\,\frac{c_s^2}{L_W^2},
qquad
\delta\ln L_W=\delta\ln a,
```

so

```math
\boxed{
\delta\ln K_q
=
\delta\ln\mathcal Z_q + 2\,\delta\ln c_s - 2\,\delta\ln a.
}
```

And on the isotropic outgoing quadrupole branch,

```math
P_0\propto \frac{c_s^5}{a^5},
```

so

```math
\boxed{
\delta\ln P_0 = 5\bigl(\delta\ln c_s - \delta\ln a\bigr).
}
```

These four equations are the entire inversion system.

---

## Step 37B — Exact inversion

Solving the system gives

```math
\boxed{
\delta\ln\rho_w
=
\frac12\,\delta\ln\Theta_w,
}
```

```math
\boxed{
\delta\ln a
=
\frac12\,\delta\ln K_s
-
\frac14\,\delta\ln\Theta_w,
}
```

```math
\boxed{
\delta\ln c_s
=
\frac12\,\delta\ln K_s
-
\frac14\,\delta\ln\Theta_w
+
\frac15\,\delta\ln P_0,
}
```

```math
\boxed{
\delta\ln\mathcal Z_q
=
\delta\ln K_q
-
\frac25\,\delta\ln P_0.
}
```

So the last four branch drifts are not open in any diffuse sense anymore. They are
exact algebraic images of the bundle observables `\Theta_w, K_s, K_q, P_0`.

---

## Step 37C — Full-bundle form using `P_0=N_0/D_0`

Because the grouped-bundle normalization is

```math
P_0=\frac{N_0}{D_0},
qquad
\delta\ln P_0 = \delta\ln N_0 - \delta\ln D_0,
```

we may rewrite the last two drifts directly in full-bundle language:

```math
\boxed{
\delta\ln c_s
=
\frac12\,\delta\ln K_s
-
\frac14\,\delta\ln\Theta_w
+
\frac15\bigl(\delta\ln N_0-\delta\ln D_0\bigr),
}
```

```math
\boxed{
\delta\ln\mathcal Z_q
=
\delta\ln K_q
-
\frac25\bigl(\delta\ln N_0-\delta\ln D_0\bigr).
}
```

So the full grouped wall/BdG/Maxwell bundle enters only through the isotropic
static normalization ratio `N_0/D_0`, exactly as the 2.5PN normalization package
had been suggesting all along.

---

## Step 37D — Frozen-wall corollary

If the explicit Family-1 wall datum is held fixed at first order,

```math
\delta\ln\Theta_w=0,
```

then the inversion simplifies to

```math
\boxed{
\delta\ln\rho_w = 0,
}
```

```math
\boxed{
\delta\ln a = \frac12\,\delta\ln K_s,
}
```

```math
\boxed{
\delta\ln c_s = \frac12\,\delta\ln K_s + \frac15\,\delta\ln P_0,
}
```

```math
\boxed{
\delta\ln\mathcal Z_q = \delta\ln K_q - \frac25\,\delta\ln P_0.
}
```

So even this explicit wall-preserving restriction already removes one of the four
would-be free drifts.

---

## Main result of the step

The lower-branch problem is now reduced as far as it can go without solving the
remaining moving-throat PDE.

The four residual microscopic drifts are exact algebraic images of the bundle
observables:

```math
\boxed{
\delta\ln\rho_w
=
\frac12\,\delta\ln\Theta_w,
qquad
\delta\ln a
=
\frac12\,\delta\ln K_s-rac14\,\delta\ln\Theta_w,
}
```

```math
\boxed{
\delta\ln c_s
=
\frac12\,\delta\ln K_s-rac14\,\delta\ln\Theta_w+rac15\,\delta\ln P_0,

\qquad
\delta\ln\mathcal Z_q
=
\delta\ln K_q-rac25\,\delta\ln P_0.
}
```

And because

```math
P_0=\frac{N_0}{D_0},
```

the only grouped-bundle quantity that still matters is the isotropic static
normalization ratio `N_0/D_0`.

So the next PDE-facing question is no longer “what are the four drifts?” It is:

> what are the first-order drifts of the four bundle observables `\Theta_w, K_s, K_q, P_0` on the actual isotropic moving-throat branch?
# Step 38 — Adiabatic-wall ground-state closure preserves the lower compensated branch

## Goal

You specified a new physical direction for the isolated-particle electron track:

- the wall is **adiabatic**,
- the wall’s thermal/entropy weighting is frozen,
- the wall density is frozen,
- the defect absorbs stress through coherent elastic squish rather than heating,
- the **core** remains compressible and can still change its stiffness.

The right next question is therefore not a numerical fit. It is structural:

> what does that closure do to the exact lower compensated Family-1 branch?

This step shows something stronger than I expected:

```math
\boxed{\text{the adiabatic-wall ground-state closure preserves the exact lower compensated parent surface identically.}}
```

So the electron-track branch does **not** fight the lower canonical outlet family at first order. It stays tangent automatically.

---

## Step 38A — Adiabatic-wall ground-state closure

The explicit Family-1 wall datum satisfies

```math
\Theta_w = 25\lambda_\mu^2\rho_w^2.
```

Your new physical choice is naturally encoded as

```math
\delta\ln\lambda_\mu = 0,
\qquad
\delta\ln\rho_w = 0.
```

Therefore

```math
\boxed{\delta\ln\Theta_w=0.}
```

So this is a **stronger** restriction than the earlier generic frozen-wall corollary: it is not just the composite datum `\Theta_w` that is frozen, but the wall-density lane itself.

---

## Step 38B — Exact lower-branch drift reduction

Using the exact Step-37 inversion with `\delta\ln\Theta_w=0`, the remaining microscopic drifts become

```math
\boxed{\delta\ln\rho_w = 0,}
```

```math
\boxed{\delta\ln a = \frac12\,\delta\ln K_s,}
```

```math
\boxed{
\delta\ln c_s
=
\frac12\,\delta\ln K_s
+\frac15\,\delta\ln P_0,
}
```

```math
\boxed{
\delta\ln\mathcal Z_q
=
\delta\ln K_q
-\frac25\,\delta\ln P_0.
}
```

So under the adiabatic-wall ground-state closure, the full lower branch already collapses from four irreducible drifts to only **three** bundle observables:

```math
\boxed{(K_s,\ K_q,\ P_0).}
```

Interpretation:

- `K_s` is the coherent elastic wall-squish lane,
- `K_q` is the mixed side-channel stiffness lane,
- `P_0` is the isotropic outgoing-normalization lane.

The wall-density lane is completely removed.

---

## Step 38C — Exact transport of `v_{w0}` and `\mathcal T_m`

From the exact lower-branch sum/difference laws,

```math
\delta\ln\!\left(\frac{v_{w0}}{\mathcal T_m}\right)
=
\frac12\,\delta\ln K_s
+\frac25\,\delta\ln P_0,
```

```math
\delta\ln(v_{w0}\mathcal T_m)
=
-2\,\delta\ln K_s
+\delta\ln K_q
-\frac25\,\delta\ln P_0.
```

Solving gives

```math
\boxed{
\delta\ln v_{w0}
=
-\frac34\,\delta\ln K_s
+\frac12\,\delta\ln K_q,
}
```

```math
\boxed{
\delta\ln\mathcal T_m
=
-\frac54\,\delta\ln K_s
+\frac12\,\delta\ln K_q
-\frac25\,\delta\ln P_0.
}
```

This is already suggestive physically:

- the mixed background flow `v_{w0}` is **blind** to the outgoing-normalization drift `P_0`,
- the mouth traction `\mathcal T_m` is where the `P_0` drift still enters directly.

So the adiabatic wall makes the transport split cleaner rather than messier.

---

## Step 38D — Parent-action transport of `(g_s,g_q,\lambda)`

Using the exact bundle transport laws,

```math
\boxed{
\delta\ln g_s
=
-\frac14\,\delta\ln K_s
+\frac12\,\delta\ln K_q
-\frac25\,\delta\ln P_0,
}
```

```math
\boxed{
\delta\ln g_q
=
-\frac34\,\delta\ln K_s
+\delta\ln K_q
-\frac25\,\delta\ln P_0,
}
```

```math
\boxed{
\delta\ln \lambda
=
\frac12\bigl(\delta\ln K_s+\delta\ln K_q\bigr).
}
```

So the bilinear shell/mixed coupling `\lambda` still ignores both the wall-depth drift and the outgoing-normalization drift, exactly as in the broader bundle theorem. Only the two stiffness lanes matter to it.

---

## Step 38E — Exact preservation of the lower compensated parent family

Now test the two exact parent-surface imbalance channels:

```math
\delta\ln\!\left(\frac{g_q K_s}{g_s\lambda}\right),
\qquad
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right).
```

Substituting the Step-38 transport laws gives

```math
\boxed{
\delta\ln\!\left(\frac{g_q K_s}{g_s\lambda}\right)=0,
}
```

```math
\boxed{
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right)=0.
}
```

Equivalently,

```math
\boxed{\delta\ln\mathfrak g = 0,}
```

```math
\boxed{\delta\ln\mathfrak r = 0,}
```

```math
\boxed{\delta\ln r_c = 0.}
```

So the exact lower compensated parent surface survives untouched.

This is the main result of the step.

The adiabatic-wall ground-state electron track does **not** generate a first-order off-family scalar. In the notation of the earlier parent-balance program,

```math
\boxed{\delta_\perp = 0.}
```

That is a strong structural simplification.

---

## Main result of the step

The adiabatic-wall ground-state choice does three useful things at once:

1. it removes the wall-density drift entirely,
2. it reduces the remaining bundle freedom to
   ```math
   (K_s, K_q, P_0),
   ```
3. and it preserves the exact lower compensated parent surface automatically.

So the electron-track branch is now much cleaner conceptually:

- the wall behaves as a coherent elastic shell,
- the core/outlet block carries the remaining response,
- and the branch does not drift off the canonical lower compensated family at first order.

That is exactly the kind of simplification we were hoping to get from a genuine thermodynamic ground-state closure.
# Step 39 — Minimal core-yield closure and the retained upper-sheet diagnostic

## Goal

Step 38 used your adiabatic-wall ground-state choice to reduce the electron-track branch to the three reduced observables

```math
(K_s,
K_q,
P_0),
```

while preserving the exact lower compensated parent surface.

The next natural question is then:

> what is the smallest additional closure that treats the wall as frozen/adiabatic while still letting the **core** yield?

The cleanest minimal answer is to freeze the transverse mixed-channel localization norm,

```math
\delta\ln\mathcal Z_q=0,
```

and let the remaining nontrivial response be carried only by

- coherent elastic wall squish, and
- outgoing/core compressibility.

This step shows that under that choice the whole reduced electron track becomes **rank 2**.

It also keeps your side note about the algebraic upper branch alive: `\mathfrak g_+` is not deleted, but it is reinterpreted as a sign-indefinite / pumped sheet rather than the passive positive-source electron branch.

---

## Step 39A — Minimal core-yield closure

Impose, in addition to Step 38,

```math
\boxed{\delta\ln\mathcal Z_q=0.}
```

This is not an exact theorem of the current file stack. It is a **minimal closure choice** consistent with your new physical picture:

- wall thermodynamics frozen,
- wall density frozen,
- coherent elastic squish allowed,
- core compressibility allowed,
- no extra side-channel relocalization introduced unless it is forced later.

Then the exact lower-branch law

```math
\delta\ln\mathcal Z_q
=
\delta\ln K_q-
\frac25\delta\ln P_0
```

immediately gives

```math
\boxed{
\delta\ln K_q
=
\frac25\,\delta\ln P_0.
}
```

So the mixed side-channel stiffness is no longer independent either.

---

## Step 39B — Exact two-amplitude reduced electron track

Using Step 38 plus `\delta\ln\mathcal Z_q=0`, the entire reduced branch becomes

```math
\boxed{
\delta\ln a = \frac12\,\delta\ln K_s,
}
```

```math
\boxed{
\delta\ln c_s
=
\frac12\,\delta\ln K_s
+\frac15\,\delta\ln P_0,
}
```

```math
\boxed{
\delta\ln K_q
=
\frac25\,\delta\ln P_0.
}
```

The transported outlet/background variables become

```math
\boxed{
\delta\ln v_{w0}
=
-\frac34\,\delta\ln K_s
+\frac15\,\delta\ln P_0,
}
```

```math
\boxed{
\delta\ln\mathcal T_m
=
-\frac54\,\delta\ln K_s
-\frac15\,\delta\ln P_0,
}
```

and the sum/difference channels are

```math
\boxed{
\delta\ln\!\left(\frac{v_{w0}}{\mathcal T_m}\right)
=
\frac12\,\delta\ln K_s
+\frac25\,\delta\ln P_0,
}
```

```math
\boxed{
\delta\ln(v_{w0}\mathcal T_m)
=
-2\,\delta\ln K_s.
}
```

So the product channel is now completely blind to `P_0`; only the elastic squish moves it.

That is a very sharp structural statement.

---

## Step 39C — Parent couplings on the two-amplitude branch

The exact parent-action couplings reduce to

```math
\boxed{
\delta\ln g_s
=
-\frac14\,\delta\ln K_s
-\frac15\,\delta\ln P_0,
}
```

```math
\boxed{
\delta\ln g_q
=
-\frac34\,\delta\ln K_s,
}
```

```math
\boxed{
\delta\ln\lambda
=
\frac12\,\delta\ln K_s
+\frac15\,\delta\ln P_0.
}
```

A useful consequence is that

```math
\delta\ln g_q
```

is completely blind to the outgoing-normalization drift `P_0` on this minimal branch. The mixed mouth coupling tracks only the elastic shell squish.

So the two surviving amplitudes now have a clean interpretation:

- `\delta\ln K_s` = **elastic squish amplitude**,
- `\delta\ln P_0` = **outgoing/core-yield amplitude**.

---

## Step 39D — Exact rank-2 statement

The reduced drift vector

```math
(
\delta\ln K_q,
\delta\ln a,
\delta\ln c_s,
\delta\ln v_{w0},
\delta\ln\mathcal T_m,
\delta\ln g_s,
\delta\ln g_q,
\delta\ln\lambda
)
```

depends only on the two amplitudes

```math
(\delta\ln K_s,
\delta\ln P_0).
```

The SymPy Jacobian of that map has rank exactly `2`.

So the adiabatic electron-track branch is now reduced to the smallest plausible parameter space we have reached so far:

```math
\boxed{
\text{electron track} 
\sim
(\text{elastic squish},\ \text{outgoing/core-yield}).
}
```

---

## Step 39E — The lower compensated surface stays exact

Even on this stronger two-amplitude branch, the exact compensation surface remains untouched:

```math
\boxed{
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)=0,
}
```

```math
\boxed{
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right)=0,
}
```

so

```math
\boxed{\delta\ln\mathfrak g = 0,}
```

```math
\boxed{\delta\ln\mathfrak r = 0,}
```

```math
\boxed{\delta\ln r_c = 0.}
```

So the minimal core-yield closure does **not** re-open the off-family scalar drift. The branch stays exactly on the lower compensated sheet.

---

## Step 39F — Keeping `\mathfrak g_+` alive as a system-wide diagnostic sheet

You asked not to throw away the algebraic upper value

```math
\mathfrak g_+^{F1}
\approx 2.79795199200529.
```

I agree that we should keep it.

The positive-source theorem only tells us that for any nonnegative normalized mouth source profile,

```math
0\le \mathfrak g[\sigma]\le 1.
```

So `\mathfrak g_+>1` means only that the upper branch is **not** the passive positive-source electron branch.

It does **not** mean the branch is mathematically meaningless.

To make that precise, write a sign-indefinite normalized source law as

```math
\sigma = \sigma_+ - \sigma_-,
\qquad
W_+ = \int \sigma_+,
\qquad
W_- = \int \sigma_-,
\qquad
W_+ - W_- = 1.
```

Since the cosine weight on the first D/N interval lies in `[0,1]`, one has

```math
\mathfrak g[\sigma]
\le W_+
=1+W_-.
```

Therefore any realization of `\mathfrak g>1` requires

```math
\boxed{W_- \ge \mathfrak g-1.}
```

Applied to the explicit Family-1 upper sheet,

```math
\boxed{
W_-^{\min}
=
\mathfrak g_+^{F1}-1
\approx 1.79795199200529,
}
```

with corresponding positive weight

```math
\boxed{
W_+^{\min}
=
\mathfrak g_+^{F1}
\approx 2.79795199200529.
}
```

So the upper sheet can be retained as a very specific diagnostic branch:

- not passive,
- not positive-source,
- requiring strong sign-indefinite / source-sink / pumped support,
- and therefore potentially telling us about a broader system-level sector rather than the isolated electron ground state.

That seems exactly in line with your side note.

---

## Main result of the step

The stronger adiabatic electron-track closure

```math
\delta\ln\mathcal Z_q=0
```

reduces the whole lower-branch problem to only two amplitudes:

```math
\boxed{
\delta\ln K_s
\quad\text{and}\quad
\delta\ln P_0.
}
```

So the g-2/electron-track search is no longer spread across many microscopic directions.
It has collapsed to

1. a coherent elastic wall-squish lane,
2. and a core/outgoing normalization lane.

At the same time, the algebraic upper sheet `\mathfrak g_+` is preserved as a flagged non-passive branch rather than deleted.

That is a good place to stand:

- the physical electron branch is cleaner,
- the non-electron algebraic sheet is still remembered,
- and the next step can finally ask how the quartic anomaly target picks out a particular combination of the two surviving amplitudes.
# Step 40 — On the strong adiabatic-wall branch the quartic anomaly becomes a tiny relative core-stiffening law

## Goal

The visible Step 39 closure already reduced the adiabatic electron branch to two amplitudes:

- `\sigma := \delta\ln K_s` — coherent elastic wall squish,
- `\ell := \delta\ln P_0` — isotropic static load / anomaly-carrying normalization drift.

So the next honest question is now very sharp:

> once the quartic g-2 sliver is inserted, what does it do to the adiabatic-wall branch microscopically?

This step answers that.

---

## Step 40A — Exact strong-closure kinematics

From the strong adiabatic-wall closure,

```math
\delta\ln\Theta_w=0,
\qquad
\delta\ln\mathcal Z_q=0,
```

the reduced branch laws are

```math
\delta\ln a = \frac12\sigma,
```

```math
\delta\ln c_s = \frac12\sigma + \frac15\ell,
```

```math
\delta\ln K_q = \frac25\ell.
```

So the relative core stiffening is

```math
\boxed{
\delta\ln\!\left(\frac{c_s}{a}\right)
=
\frac15\ell
=
\frac12\,\delta\ln K_q.
}
```

This is already strong physically:

- the adiabatic wall direction `\sigma` drops out completely from `c_s/a`,
- the only thing that changes the core stiffness **relative to wall size** is the anomaly-carrying load drift `\ell`.

---

## Step 40B — Exact outgoing-normalization law from the quartic sliver

The earlier quartic bridge already rewrote the anomaly miss as a tiny outgoing-normalization defect,

```math
\Delta_Q(f)= -\frac{\Lambda_1 f}{1+\Lambda_1 f},
```

with

```math
\Lambda_1\approx 0.279605891931464.
```

On the natural source-map branch,

```math
N_Q = \frac{1}{1+\Delta_Q}.
```

Substituting the exact `\Delta_Q(f)` gives

```math
\boxed{N_Q(f)=1+\Lambda_1 f.}
```

So on the strong adiabatic-wall branch,

```math
\boxed{
\ell
=
\delta\ln P_0
=
\ln N_Q
=
\ln(1+\Lambda_1 f).
}
```

At linear order,

```math
\ell \approx \Lambda_1 f,
```

but the exact logarithmic form is the cleaner finite-`f` law once the outgoing defect is already known.

---

## Step 40C — Exact microscopic anomaly laws on the strong adiabatic-wall branch

Substituting

```math
\ell = \ln(1+\Lambda_1 f)
```

into the strong closure gives

```math
\boxed{
\delta\ln K_q
=
\frac25\ln(1+\Lambda_1 f),
}
```

```math
\boxed{
\delta\ln\!\left(\frac{c_s}{a}\right)
=
\frac15\ln(1+\Lambda_1 f).
}
```

So the quartic anomaly is not asking for a broad thermodynamic rearrangement.
On this branch it asks only for a very small **relative core stiffening**.

That is the sharpest physical statement reached so far:

> with an adiabatic wall, the electron quartic sliver is carried by a tiny compressibility shift of the core relative to the wall size, while the coherent wall squish remains orthogonal.

---

## Step 40D — Electron-point numbers

Using the carried electron value

```math
f \approx 0.001161409732093,
```

one gets

```math
\Lambda_1 f \approx 3.24737004039746\times 10^{-4},
```

```math
\ell
=
\ln(1+\Lambda_1 f)
\approx
3.24684288391064\times 10^{-4},
```

```math
\delta\ln K_q
\approx
1.29873715356426\times 10^{-4},
```

```math
\delta\ln(c_s/a)
\approx
6.49368576782128\times 10^{-5}.
```

So the required correction is genuinely tiny.
The exact-vs-linear difference in `\ell` is only about

```math
5.27\times 10^{-8},
```

which means the linearized quartic bookkeeping is already numerically excellent.

---

## What remains free

This is also the first point where the branch splits cleanly into two jobs:

1. the thermodynamic ground state / coherent wall shape chooses
   ```math
   \sigma = \delta\ln K_s,
   ```
2. the quartic anomaly closure chooses
   ```math
   \ell = \delta\ln P_0.
   ```

So the anomaly does **not** determine the whole branch.
It only determines the very small compressibility/load correction.
The coherent elastic squish remains an orthogonal ground-state datum.

That is exactly what one would hope if the wall is truly adiabatic.
# Step 41 — Once the adiabatic ground state is fixed, the quartic anomaly is a pure core/outgoing retuning

## Goal

Step 40 showed that on the strong adiabatic-wall branch the quartic anomaly fixes only

```math
\ell := \delta\ln P_0 = \ln(1+\Lambda_1 f),
```

while the coherent elastic squish

```math
\sigma := \delta\ln K_s
```

remains orthogonal.

The natural next move is therefore to add one more physical closure:

> the thermodynamic ground state has already chosen `\sigma`, so the **incremental quartic anomaly correction** should be taken with
> ```math
> \delta\sigma = 0.
> ```

This step shows that, on that frozen ground state, the anomaly becomes a **pure core/outgoing retuning**.

---

## Step 41A — Exact anomaly-only drift vector

Setting

```math
\sigma=0
```

in the strong adiabatic-wall branch laws gives

```math
\boxed{\delta\ln a = 0,}
```

```math
\boxed{\delta\ln c_s = \frac15\ell,}
```

```math
\boxed{\delta\ln K_q = \frac25\ell,}
```

```math
\boxed{\delta\ln v_{w0} = \frac15\ell,}
```

```math
\boxed{\delta\ln\mathcal T_m = -\frac15\ell,}
```

```math
\boxed{\delta\ln g_s = -\frac15\ell,}
```

```math
\boxed{\delta\ln g_q = 0,}
```

```math
\boxed{\delta\ln\lambda = \frac15\ell.}
```

So the pure anomaly increment has a very clean structure:

- it does **not** move the mouth radius,
- it does **not** change the mixed coupling `g_q`,
- it increases the core sound scale,
- and it retunes the outgoing / traction side in equal-and-opposite logarithmic amounts.

---

## Step 41B — The lower compensated branch remains exact

Even after freezing the ground-state squish, the lower compensated invariants remain

```math
\boxed{\delta\ln\mathfrak g = 0,}
```

```math
\boxed{\delta\ln\mathfrak r = 0,}
```

```math
\boxed{\delta\ln r_c = 0.}
```

So the anomaly-only increment is still tangent to the same lower compensated sheet.
It does **not** push the system off the lower Family-1 branch.

That is a useful structural fact:

> the adiabatic-ground-state anomaly correction is an **internal retuning along the already selected electron branch**, not a branch-jump.

---

## Step 41C — Electron-point numbers

Using

```math
\ell = \ln(1+\Lambda_1 f)
approx 3.24684288391064\times 10^{-4},
```

the anomaly-only drifts are

```math
\delta\ln c_s
approx 6.49368576782128\times 10^{-5},
```

```math
\delta\ln K_q
approx 1.29873715356426\times 10^{-4},
```

```math
\delta\ln v_{w0}
approx 6.49368576782128\times 10^{-5},
```

```math
\delta\ln\mathcal T_m
approx -6.49368576782128\times 10^{-5},
```

```math
\delta\ln g_s
approx -6.49368576782128\times 10^{-5},
```

```math
\delta\ln g_q = 0,
```

```math
\delta\ln\lambda
approx 6.49368576782128\times 10^{-5}.
```

So once the adiabatic ground state is fixed, the anomaly correction is extremely localized in parameter space.
It is only an `O(10^{-5})` to `O(10^{-4})` retuning of the core/outgoing sector.

---

## Step 41D — The upper `g_+` sheet stays on the ledger, but the electron anomaly does not need it

The algebraic upper branch is still

```math
\mathfrak g_+^{F1} \approx 2.797951992,
```

and, as before, any realization of it requires a sign-indefinite or pumped source law with at least

```math
W_- \ge \mathfrak g_+ - 1 \approx 1.797951992.
```

So it still does **not** belong to the passive positive-source isolated-electron branch.

But it is worth keeping for later because it can still represent a deferred system-level sheet:

- a sign-changing mouth source,
- a pumped / plumbed configuration,
- or a non-electron excitation branch.

The important point for the present g-2 problem is narrower:

> the pure anomaly increment derived above does **not** move the system toward that upper sheet. It stays exactly tangent to the lower compensated electron branch.

---

## Main result of the step

With the adiabatic wall fixed to its ground state, the quartic electron anomaly reduces to a **pure core/outgoing retuning**:

```math
\boxed{
\delta\ln a = 0,
\qquad
\delta\ln g_q = 0,
\qquad
\delta\ln c_s = \delta\ln v_{w0} = \delta\ln\lambda = \frac15\ell,
\qquad
\delta\ln K_q = \frac25\ell,
\qquad
\delta\ln\mathcal T_m = \delta\ln g_s = -\frac15\ell.
}
```

with

```math
\ell = \ln(1+\Lambda_1 f).
```

So the improved g-2 picture is now very clean:

- the wall stays adiabatic,
- the ground-state elastic squish stays fixed,
- the anomaly rides only on a tiny core/outgoing correction,
- and the lower compensated branch remains exact while the upper `g_+` sheet stays deferred for later system-level interpretation.
# Step 42 — The adiabatic anomaly increment freezes the parent compensation sheet and detunes only the odd outlet

## Goal

Step 41 showed that once the adiabatic wall is fixed to its thermodynamic ground state, the quartic anomaly is a **pure core/outgoing retuning**:

```math

delta\ln K_s = 0,
\qquad
\ell := \delta\ln P_0,
```

with

```math
\delta\ln a = 0,
\qquad
\delta\ln c_s = \frac15\ell,
\qquad
\delta\ln K_q = \frac25\ell,
```

```math
\delta\ln v_{w0} = \frac15\ell,
\qquad
\delta\ln\mathcal T_m = -\frac15\ell,
\qquad
\delta\ln g_s = -\frac15\ell,
\qquad
\delta\ln g_q = 0,
\qquad
\delta\ln\lambda = \frac15\ell.
```

The next honest question is therefore sharper than before:

> does this anomaly-only increment move the system around on the parent compensation family, or does it stay on the same branch and only retune the outlet?

This step shows the stronger result:

```math
\boxed{
\text{the anomaly-only adiabatic branch leaves }\mathfrak g,\ \mathfrak r,\ r_c,\ L_W\ \text{all fixed.}
}
```

So the quartic sliver does **not** come from motion in the parent compensation ratios. It comes only from a small odd-outlet detuning on top of a slightly softened loading share.

---

## Step 42A — Exact parent-sheet invariants are frozen

The lower compensated family is parameterized by the parent ratios

```math
\mathfrak r := \frac{\lambda}{\sqrt{K_s K_q}},
\qquad
\mathfrak g := \frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}},
\qquad
r_c := \frac{\lambda^2}{K_s K_q} = \mathfrak r^2.
```

Substituting the Step-41 anomaly-only drift vector gives

```math
\boxed{
\delta\ln\mathfrak r
=
\delta\ln\lambda-
\frac12(\delta\ln K_s+\delta\ln K_q)
=0,
}
```

```math
\boxed{
\delta\ln\mathfrak g
=
\delta\ln g_q+
\frac12\delta\ln K_s-
\delta\ln g_s-
\frac12\delta\ln K_q
=0,
}
```

```math
\boxed{
\delta\ln r_c
=
2\,\delta\ln\lambda-
\delta\ln K_s-
\delta\ln K_q
=0.
}
```

So the anomaly-only adiabatic increment stays **exactly tangent to the same lower compensated parent sheet**.

That means the quartic g-2 sliver is **not** telling the electron branch to move to a new `\mathfrak g`/`\mathfrak r` balance point.

---

## Step 42B — The auxiliary D/N geometry stays fixed

The balanced auxiliary side-tube geometry obeys

```math
L_W = \frac{\pi a}{2}\sqrt{\frac{1+r_c}{3}}.
```

Since Step 41 already gave

```math
\delta\ln a = 0,
\qquad
\delta\ln r_c = 0,
```

we get immediately

```math
\boxed{\delta\ln L_W = 0.}
```

So the quartic anomaly does **not** ask for a new even D/N geometry either.
The balanced auxiliary tube length remains frozen.

This is important because it means the anomaly is not reopening the even outlet side of the problem.

---

## Step 42C — The balanced loading share drifts downward

Although the parent ratios stay fixed, the shell-side loading share does not.
On the balanced core,

```math
\rho_c := \frac{g_s^2}{K_s},
\qquad
\sigma_c = \frac{\rho_c}{4}.
```

Using

```math
\delta\ln g_s = -\frac15\ell,
\qquad
\delta\ln K_s = 0,
```

gives

```math
\boxed{
\delta\ln\rho_c = 2\,\delta\ln g_s - \delta\ln K_s = -\frac25\ell,
}
```

```math
\boxed{
\delta\ln\sigma_c = -\frac25\ell.
}
```

So the exact finite-`\ell` laws are

```math
\boxed{
\rho_c(\ell)=\rho_{c,*} e^{-2\ell/5},
\qquad
\sigma_c(\ell)=\sigma_{c,*} e^{-2\ell/5}.
}
```

Meanwhile the even-preserving mixed-side coefficient stays frozen:

```math
\boxed{\kappa_c = \frac13.}
```

So the anomaly-only adiabatic branch keeps the same **parent balance** and the same **even side-channel geometry**, but it reduces the loading share slightly because `g_s` softens while `K_s` is held fixed.

---

## Step 42D — Exact odd detuning law on the fixed parent sheet

The electron anomaly target on the outgoing branch is

```math
\chi_Q = e^{-\ell} = \frac{1}{1+x},
\qquad
x := e^{\ell}-1 = \Lambda_1 f.
```

On the balanced core / compensated outlet,

```math
\chi_Q = \frac{1-9\sigma_c\gamma_c}{1-\sigma_c}.
```

Solving for the odd coefficient gives

```math
\boxed{
\gamma_c(\ell)
=
\frac{1-e^{-\ell}(1-\sigma_c(\ell))}{9\sigma_c(\ell)}.
}
```

Using `x = e^{\ell}-1` this can be rewritten more compactly as

```math
\boxed{
\gamma_c
=
\frac{\sigma_c + x}{9\sigma_c(1+x)},
\qquad
\sigma_c = \sigma_{c,*}(1+x)^{-2/5}.
}
```

The corresponding **bare** odd coefficient is just

```math
\boxed{
\gamma_0 = (1+r_{c,*})\,\gamma_c,
}
```

because the parent hybridization ratio `r_c` is fixed.

So the anomaly does not move the branch in `\mathfrak g` or `\mathfrak r`; it forces only the odd outlet coefficient to detune so that the softened loading share still lands on the electron target.

At first order,

```math
\gamma_c
=
\frac19
+
\frac{1-\sigma_{c,*}}{9\sigma_{c,*}}\,\ell
+
O(\ell^2).
```

So the required shift is indeed small when `\ell` is small.

---

## Step 42E — Same result in compensated Robin–mixed outlet variables

On the compensated hybrid outlet,

```math
\rho_R = 4\sigma_W,
\qquad
\kappa_W = \frac13,
\qquad
\chi_Q = \frac{1-9\sigma_W\gamma_W}{1-\sigma_W}.
```

Identifying `\sigma_W=\sigma_c`, the adiabatic anomaly-only branch gives

```math
\boxed{
\sigma_W(\ell)=\sigma_{W,*}e^{-2\ell/5},
\qquad
\rho_R(\ell)=\rho_{R,*}e^{-2\ell/5},
\qquad
\kappa_W=\frac13,
}
```

and the exact odd detuning law

```math
\boxed{
\gamma_W
=
\frac{\sigma_W + x}{9\sigma_W(1+x)},
\qquad
x=\Lambda_1 f.
}
```

So the hybrid outlet statement is the same as the balanced-core statement:

- even branch fixed,
- parent compensation ratios fixed,
- loading share slightly reduced,
- odd outlet coefficient slightly increased above its canonical value.

---

## Main result of the step

The adiabatic anomaly-only branch has a very sharp parent-side meaning:

```math
\boxed{
\delta\ln\mathfrak g = 0,
\qquad
\delta\ln\mathfrak r = 0,
\qquad
\delta\ln r_c = 0,
\qquad
\delta\ln L_W = 0.
}
```

So the quartic electron sliver is **not** a motion of the parent compensation family or the even D/N geometry.

Instead it is:

```math
\boxed{
\sigma_c(\ell)=\sigma_{c,*}e^{-2\ell/5},
\qquad
\gamma_c
=
\frac{\sigma_c + \Lambda_1 f}{9\sigma_c(1+\Lambda_1 f)},
}
```

or equivalently in compensated Robin–mixed notation,

```math
\boxed{
\rho_R(\ell)=\rho_{R,*}e^{-2\ell/5},
\qquad
\kappa_W=\frac13,
\qquad
\gamma_W
=
\frac{\sigma_W + \Lambda_1 f}{9\sigma_W(1+\Lambda_1 f)}.
}
```

So the improved g-2 picture is now cleaner again:

- the parent compensation sheet stays fixed,
- the even side-channel stays fixed,
- the load share softens slightly,
- and the anomaly is carried only by a tiny odd outlet detuning.
# Step 43 — The adiabatic anomaly admits a pure odd normalized DtN representative

## Goal

Step 42 showed that on the anomaly-only adiabatic branch

- the parent compensation ratios stay fixed,
- the auxiliary D/N geometry stays fixed,
- and only the odd outlet coefficient detunes.

That still lives in the **microscopic** compensated-core / hybrid-outlet variables.
The next honest question is therefore:

> what does the same anomaly look like after we reduce it all the way back to the **normalized isotropic DtN branch-selection surface**?

This step shows an unexpectedly clean answer:

```math
\boxed{
\text{the entire adiabatic anomaly has an exact pure-odd normalized DtN representative.}
}
```

In that representative the even branch is frozen exactly and the whole observable branch motion sits in one isotropic odd slot.

---

## Step 43A — Exact isotropic DtN deformation surface

From the earlier isotropic DtN reduction,

```math
\chi_Q
=
\frac{3\bigl(S\beta^5+9\Sigma_5\bigr)}{3S-\Sigma_0},
```

with

- `S` = overall isotropic mouth normalization,
- `\beta` = outgoing-argument deformation,
- `\Sigma_0` = static isotropic additive deformation,
- `\Sigma_5` = isotropic odd `l=2` core outlet.

The adiabatic anomaly target is still

```math
\chi_Q = e^{-\ell} = \frac{1}{1+x},
\qquad x:=e^{\ell}-1=\Lambda_1 f.
```

So the problem is now very narrow: how is this exact target represented on the DtN surface once the adiabatic anomaly-only closure freezes the even branch?

---

## Step 43B — Fixed-even slice gives a pure odd exact law

Because Step 42 already proved that the anomaly-only branch does **not** move the parent compensation ratios or the even D/N geometry, the natural normalized DtN slice is

```math
\boxed{\beta = 1,\qquad \Sigma_0 = 0.}
```

Then the isotropic DtN surface collapses to

```math
\chi_Q = 1 + \frac{9\Sigma_5}{S}.
```

Matching this to `\chi_Q=e^{-\ell}` gives the exact pure-odd law

```math
\boxed{
\Sigma_5(\ell)
=
\frac{S}{9}\bigl(e^{-\ell}-1\bigr)
=
-\frac{Sx}{9(1+x)}.
}
```

So once the even branch is frozen, the exact electron anomaly is represented by **one** pure odd isotropic outlet coefficient and nothing else.

---

## Step 43C — This pure odd representative is independent of the microscopic loading share

On the compensated hybrid outlet we also have

```math
\chi_Q = \frac{1-9\sigma_W\gamma_W}{1-\sigma_W}.
```

Step 42 showed that `\sigma_W` drifts and `\gamma_W` retunes so that this still equals the electron target.

Now form the normalized DtN representative directly from `\chi_Q`:

```math
\Sigma_5 = \frac{S}{9}(\chi_Q-1).
```

After substituting the exact hybrid-outlet solution for `\gamma_W`, all explicit dependence on `\sigma_W` cancels, leaving

```math
\boxed{
\Sigma_5(\ell)=\frac{S}{9}(e^{-\ell}-1).
}
```

So the microscopic loading-share drift is real on the compensated-core side, but once the branch is repackaged as a **normalized outgoing DtN deformation**, the observable anomaly is carried entirely by the pure odd slot.

That is the cleanest reduced statement reached so far.

---

## Step 43D — Canonical normalized gauge and the tangent law

If we choose the canonical normalized gauge

```math
S=1,
```

then the exact finite-`f` odd slot is

```math
\boxed{
\Sigma_5
=
-\frac{\Lambda_1 f}{9(1+\Lambda_1 f)}.
}
```

Numerically at the electron point,

```math
\boxed{
\Sigma_5^{(e)}
\approx
-3.60701760168546\times 10^{-5}.
}
```

Now expand the full isotropic DtN surface with

```math
S=1,
\qquad
\beta = 1+\varepsilon b,
\qquad
\Sigma_0 = \varepsilon a_0,
\qquad
\Sigma_5 = \varepsilon a_5,
\qquad
\chi_Q = 1 - \varepsilon\Lambda_1 + O(\varepsilon^2).
```

The exact first-order constraint is the already-known tangent law

```math
5b + \frac{a_0}{3} + 9a_5 = -\Lambda_1.
```

On the fixed-even adiabatic slice `b=0`, `a_0=0`, this reduces to

```math
\boxed{a_5 = -\frac{\Lambda_1}{9}.}
```

So the adiabatic anomaly-only branch picks the pure odd tangent representative

```math
\boxed{b=0,\qquad a_0=0,\qquad a_5=-\Lambda_1/9.}
```

---

## Step 43E — What changed conceptually

There are now two equivalent descriptions of the same quartic anomaly layer.

### Microscopic compensated-core / hybrid-outlet language

The anomaly slightly softens the loading share and slightly retunes the odd outlet coefficient, while leaving the parent compensation ratios and even geometry fixed.

### Normalized isotropic DtN language

The exact same branch motion can be represented by

```math
\boxed{
\beta = 1,
\qquad
\Sigma_0 = 0,
\qquad
\Sigma_5 = \frac{S}{9}(e^{-\ell}-1).
}
```

So at the normalized outgoing-DtN level, the anomaly is **pure odd**.

That is a very useful simplification for the next moving-throat step, because it means the PDE-facing target is no longer “find a mix of even and odd isotropic deformations.” It is just:

> compute the physical odd quadrupole outlet coefficient on the true branch and see whether it lands on the pure-odd target above.

---

## Main result of the step

The adiabatic anomaly-only branch admits an exact pure-odd normalized DtN representative:

```math
\boxed{
\Sigma_5(\ell)=\frac{S}{9}(e^{-\ell}-1)=-\frac{S\Lambda_1 f}{9(1+\Lambda_1 f)},
\qquad
\beta=1,
\qquad
\Sigma_0=0.
}
```

In canonical normalized gauge,

```math
\boxed{
\Sigma_5
=
-\frac{\Lambda_1 f}{9(1+\Lambda_1 f)},
}
```

and at tangent level,

```math
\boxed{b=0,\qquad a_0=0,\qquad a_5=-\Lambda_1/9.}
```

So the improved g-2 picture is now extremely compressed:

- the adiabatic ground state freezes the wall and the parent compensation sheet,
- the microscopic core/hybrid outlet retunes only its odd slot,
- and the normalized outgoing DtN branch sees the whole anomaly as one exact pure-odd isotropic coefficient.
# Step 44 — The adiabatic pure-odd DtN deformation rescales the grouped-`P2` normalization product by `e^\ell`

## Goal

Step 43 compressed the quartic anomaly into a very clean normalized DtN statement:

```math
\chi_Q = e^{-\ell},
\qquad
\beta = 1,
\qquad
\Sigma_0 = 0,
\qquad
\Sigma_5 = \frac{S}{9}(e^{-\ell}-1).
```

The next honest move is to connect that result back to the moving-throat grouped-`P2`
normalization product. The right bridge already existed in reduced form:

```math
\hat m_0^{\,2}\,\chi_Q\,N_Q = 1.
```

This step shows that on the natural point-particle source-map branch the pure-odd
DtN deformation forces the grouped normalization defect to move in the **opposite**
direction,

```math
\boxed{N_Q=e^{+\ell},}
```

and therefore the moving-throat grouped-`P2` prefactor obeys

```math
\boxed{
P_0 = e^\ell\,P_0^{\rm target},
\qquad
P_0^{\rm target}=\frac{54Gc_s^5}{5a^5c^5}.
}
```

So the pure-odd outlet and the grouped-`P2` normalization product are reciprocal
descriptions of the same quartic branch motion.

---

## Step 44A — Reciprocal normalization law on the natural source-map branch

From the reduced factorization,

```math
\hat m_0^{\,2}\,\chi_Q\,N_Q = 1,
```

together with the Step-43 pure-odd law

```math
\chi_Q=e^{-\ell},
```

we get

```math
N_Q = \frac{e^\ell}{\hat m_0^{\,2}}.
```

On the natural point-particle source-map branch used throughout the reduced
electron closure,

```math
\hat m_0 = 1,
```

so this collapses to

```math
\boxed{N_Q=e^\ell.}
```

That already tells us something important: the same quartic anomaly that *reduces*
the normalized outlet factor `\chi_Q` *increases* the grouped normalization defect
`N_Q` by the inverse amount.

---

## Step 44B — The grouped-`P2` normalization product

From the grouped-`P2` normalization side, Step 27 had already shown that the single
reduced defect rescales the odd grouped coefficient:

```math
\bar\Gamma_5 = N_Q\,\bar\Gamma_5^{\rm target}
= N_Q\,\frac{2G}{5c^5}.
```

From the moving-throat grouped-`P2` bridge, the same coefficient is written as

```math
\bar\Gamma_5 = \frac{P_0 a^5}{27c_s^5}.
```

So the static prefactor is

```math
P_0 = N_Q\,\frac{54Gc_s^5}{5a^5c^5}.
```

On the natural source-map branch this becomes

```math
\boxed{
P_0=e^\ell\,P_0^{\rm target},
\qquad
P_0^{\rm target}:=\frac{54Gc_s^5}{5a^5c^5}.
}
```

So the Step-43 pure-odd branch selection is *exactly* the same branch motion as a
static grouped-`P2` prefactor rescaling by `e^\ell`.

---

## Step 44C — Exact reciprocal law

Multiplying the two descriptions together gives

```math
\boxed{
P_0\,\chi_Q = P_0^{\rm target}.
}
```

That is the cleanest statement of the bridge:

- in normalized DtN language the anomaly is a **pure odd outlet deformation**
  `\chi_Q=e^{-\ell}`,
- in grouped-`P2` normalization language the same anomaly is a **static prefactor
  enhancement** `P_0=e^\ell P_0^{\rm target}`.

Those are exact reciprocals on the natural source-map branch.

---

## Step 44D — Electron-point size

Using the carried electron value

```math
\ell = \ln(1+\Lambda_1 f),
```

we have

```math
e^\ell = 1+\Lambda_1 f.
```

Numerically,

```math
\ell \approx 3.24684288391064\times 10^{-4},
```

```math
\boxed{
e^\ell \approx 1.00032473700404.
}
```

So the grouped-`P2` normalization product is enhanced by about

```math
\boxed{324.737\ {\rm ppm}.}
```

That is the exact finite-`f` bridge from the quartic anomaly sliver to the
moving-throat normalization product.

---

## Main result of the step

The pure-odd adiabatic DtN deformation is not just an outlet-side statement.
On the natural point-particle source-map branch it forces the grouped-`P2`
normalization product to rescale inversely:

```math
\boxed{
\chi_Q=e^{-\ell}
\quad\Longleftrightarrow\quad
N_Q=e^\ell
\quad\Longleftrightarrow\quad
P_0=e^\ell\,P_0^{\rm target}.
}
```

with

```math
\boxed{
P_0^{\rm target}=\frac{54Gc_s^5}{5a^5c^5}.
}
```

So the quartic anomaly is now directly attached to the exact moving-throat
normalization product rather than only to an abstract DtN branch-selection scalar.
# Step 45 — On the adiabatic anomaly track the moving-throat normalization law survives as a retuned target curve

## Goal

Step 44 showed that the pure-odd DtN deformation from Step 43 does not stay only on
the outlet side. On the natural point-particle source-map branch it forces the
grouped-`P2` prefactor to rescale as

```math
P_0 = e^\ell P_0^{\rm target},
\qquad
P_0^{\rm target}=\frac{54Gc_s^5}{5a^5c^5}.
```

But Steps 40–41 had already shown that the same adiabatic anomaly-only branch retunes
the core sound scale by

```math
\delta\ln a = 0,
\qquad
\delta\ln c_s = \frac{\ell}{5}.
```

So the next honest question is:

> is Step 44 telling us that the anomaly moves **off** the universal moving-throat
> normalization law, or does it merely move the branch **along the same law** by
> retuning `c_s`?

This step shows the stronger and cleaner result:

```math
\boxed{
P_0(a,c_s e^{\ell/5}) = e^\ell P_0(a,c_s).
}
```

So the adiabatic anomaly rides the **same universal target curve**. It does not
break it.

---

## Step 45A — Exact retuned-target law

Take the universal moving-throat grouped-`P2` normalization curve

```math
\boxed{
P_0^{\rm target}(a,c_s)=\frac{54Gc_s^5}{5a^5c^5}.
}
```

Now insert the anomaly-only adiabatic retuning from Step 41,

```math
a \mapsto a,
\qquad
c_s \mapsto c_s e^{\ell/5}.
```

Then

```math
P_0^{\rm target}(a,c_s e^{\ell/5})
=
\frac{54G(c_s e^{\ell/5})^5}{5a^5c^5}
=
e^\ell P_0^{\rm target}(a,c_s).
```

So the Step-44 prefactor enhancement is **exactly** the same as evaluating the
same universal target formula at the retuned adiabatic branch values.

That is the cleanest physical reading so far:

> the quartic anomaly does not ask for a new normalization law; it asks for a very
> small retuning of the core sound scale along the same normalization law.

---

## Step 45B — Minimal constant-prefactor grouped-`P2` continuation

The pure-odd DtN representative from Step 43 does **not** by itself force a unique
choice of the higher grouped-`P2` prefactor moments `(P_2,P_4)`. But the moving-throat
notes already isolate one especially simple continuation:

```math
\boxed{P_2 = 0,\qquad P_4 = 0.}
```

On that minimal isotropic branch,

```math
\boxed{
K_2 = \frac{P_0 a^2}{9c_s^2},
\qquad
K_4 = \frac{4P_0 a^4}{81c_s^4}.
}
```

Substituting the adiabatic retuning `c_s\mapsto c_s e^{\ell/5}` with fixed `a`
gives the finite scaling laws

```math
\boxed{
P_0 \mapsto e^\ell P_0,
\qquad
K_2 \mapsto e^{3\ell/5} K_2,
\qquad
K_4 \mapsto e^{\ell/5} K_4.
}
```

So on the minimal grouped-`P2` continuation the quartic anomaly induces:

- a `324.7` ppm upward shift in the static prefactor,
- a `194.8` ppm upward shift in the quadratic even coefficient,
- and a `64.94` ppm upward shift in the quartic even coefficient.

At tangent level this is just

```math
\boxed{
\delta\ln P_0 = \ell,
\qquad
\delta\ln K_2 = \frac{3}{5}\ell,
\qquad
\delta\ln K_4 = \frac{1}{5}\ell.
}
```

---

## Step 45C — Exact constant-prefactor bundle conditions

The moving-throat grouped-`P2` bundle also gives exact formulas for the prefactor
moments in terms of the conservative bundle data:

```math
P_0 = \frac{N_0}{D_0},
```

```math
P_2 = \frac{D_0N_2 - 2D_2N_0}{D_0^2},
```

```math
P_4 =
\frac{D_0^2N_4 - 2D_0(D_2N_2+D_4N_0)+3D_2^2N_0}{D_0^3}.
```

The conservative normalized response moments are

```math
u_2 = -\frac{D_2}{D_0},
\qquad
u_4 = \frac{D_2^2 - D_0 D_4}{D_0^2}.
```

If we now impose the minimal constant-prefactor branch

```math
P_2=P_4=0,
```

the outgoing transfer moments are **not** free. They are slaved to the conservative
grouped response:

```math
\boxed{
N_2 = -2u_2 N_0,
}
```

```math
\boxed{
N_4 = (3u_2^2 - 2u_4)\,N_0.
}
```

So if the true moving-throat branch really selects this minimal continuation, the
PDE no longer gets to choose `N_2` and `N_4` independently once the conservative
bundle moments `(u_2,u_4)` and the static transfer weight `N_0` are known.

That is a sharp next theorem target.

---

## Step 45D — What is established, and what is still conditional

### Established now

The adiabatic anomaly branch is fully compatible with the universal moving-throat
normalization curve:

```math
\boxed{
P_0^{\rm actual}
=
\frac{54G c_s^5}{5a^5c^5}
}
```

provided `c_s` and `a` are understood as the **retuned adiabatic branch values**,
with

```math
\delta\ln a=0,
\qquad
\delta\ln c_s=\ell/5.
```

So the anomaly does **not** ask for a new normalization law.

### Still conditional

The extra step

```math
P_2=P_4=0
```

is the **minimal grouped-`P2` completion**, not yet a theorem forced by Step 43
alone. More complicated isotropic branches with `P_2,P_4\neq0` are still logically
available until the full moving-throat PDE fixes them.

So the honest status is:

- the target curve is now exact on the adiabatic anomaly branch,
- the constant-prefactor continuation is the cleanest next completion,
- but the PDE still has to decide whether that minimal branch is the actual one.

---

## Main result of the step

The adiabatic quartic anomaly does **not** push the moving-throat quadrupole
normalization off the universal target curve. It simply retunes the core sound
scale so that

```math
\boxed{
P_0(a,c_s e^{\ell/5}) = e^\ell P_0(a,c_s).
}
```

If we also choose the minimal constant-prefactor isotropic grouped-`P2` continuation,

```math
\boxed{
P_2=P_4=0,
}
```

then

```math
\boxed{
P_0 \mapsto e^\ell P_0,
\qquad
K_2 \mapsto e^{3\ell/5}K_2,
\qquad
K_4 \mapsto e^{\ell/5}K_4,
}
```

and the outgoing transfer moments are forced to obey

```math
\boxed{
N_2 = -2u_2 N_0,
\qquad
N_4 = (3u_2^2 - 2u_4)N_0.
}
```

So the next moving-throat question is no longer about the existence of a
normalization curve. It is whether the true PDE selects this minimal constant-
prefactor isotropic completion of the already-fixed adiabatic anomaly track.
# Step 46 — The adiabatic constant-prefactor outgoing bundle lies on a one-parameter algebraic target surface

## Goal

Step 45 showed two things at once:

1. on the moving-throat side, the natural isotropic outgoing branch is especially simple when the prefactor is constant,
   \[
   P_2=P_4=0,
   \]
   so that
   \[
   K_0=P_0,\qquad
   K_2=\frac{P_0 a^2}{9c_s^2},\qquad
   K_4=\frac{4P_0 a^4}{81c_s^4},\qquad
   \Gamma_5=\frac{P_0 a^5}{27c_s^5},
   \]
   with the universal normalization condition
   \[
   P_0=\frac{54Gc_s^5}{5a^5c^5}
   \]
   on the natural source-map branch; and

2. on the adiabatic anomaly track, that same normalization law survives as a
   retuned target curve rather than being replaced by a new rule.

The next clean move is therefore to eliminate the intermediate variables
`(P_0,a,c_s)` and ask:

> what **purely algebraic surface** in the observable outgoing bundle
> \((K_0,K_2,K_4,\Gamma_5)\) is selected by the adiabatic constant-prefactor branch?

This step shows that the answer is extremely rigid.

---

## Step 46A — Exact bundle values on the universal target curve

Start from the constant-prefactor branch

\[
K_0=P_0,
\qquad
K_2=\frac{P_0 a^2}{9c_s^2},
\qquad
K_4=\frac{4P_0 a^4}{81c_s^4},
\qquad
\Gamma_5=\frac{P_0 a^5}{27c_s^5},
\]
together with the universal normalization target
\[
P_0=\frac{54Gc_s^5}{5a^5c^5}.
\]

Substituting immediately gives

\[
\boxed{
K_0=\frac{54Gc_s^5}{5a^5c^5},
}
\]

\[
\boxed{
K_2=\frac{6Gc_s^3}{5a^3c^5},
}
\]

\[
\boxed{
K_4=\frac{8Gc_s}{15ac^5},
}
\]

\[
\boxed{
\Gamma_5=\frac{2G}{5c^5}.
}
\]

So the whole outgoing bundle collapses to a single scale ratio
\[
s:=\frac{c_s}{a}.
\]

In terms of \(s\),

\[
\boxed{
K_0=\frac{54G}{5c^5}s^5,\qquad
K_2=\frac{6G}{5c^5}s^3,\qquad
K_4=\frac{8G}{15c^5}s,\qquad
\Gamma_5=\frac{2G}{5c^5}.
}
\]

That is already a strong result: the adiabatic constant-prefactor bundle is not a
generic 4-parameter family. It is a **one-parameter leaf** with a universal odd slot.

---

## Step 46B — Exact algebraic surface in observable bundle space

Because the bundle depends only on the single ratio \(s=c_s/a\), the observable
coefficients satisfy exact algebraic relations.

The simplest is

\[
\boxed{
K_2^2=\frac14\,K_0K_4.
}
\]

Equivalently,

\[
\boxed{
K_0=\frac{4K_2^2}{K_4}.
}
\]

Using the universal odd coefficient \(\Gamma_5=2G/(5c^5)\), the same one-parameter
leaf can also be written as

\[
\boxed{
K_2=\frac{81}{64}\frac{K_4^3}{\Gamma_5^2},
}
\]

\[
\boxed{
K_0=\frac{6561}{1024}\frac{K_4^5}{\Gamma_5^4}.
}
\]

So the adiabatic constant-prefactor branch is a **codimension-2 algebraic target
surface** inside the 4-dimensional coefficient space \((K_0,K_2,K_4,\Gamma_5)\):

- one condition fixes the odd slot universally,
- one condition ties the three even slots together.

That is exactly the kind of theorem gate a future PDE output can be tested against
without ever referring back to \((P_0,a,c_s)\) explicitly.

---

## Step 46C — Physical reading

This outgoing bundle has a very clean interpretation.

### The odd coefficient is universal

\[
\boxed{\Gamma_5=\frac{2G}{5c^5}}
\]

is independent of the adiabatic retuning parameter.
So the universal Burke–Thorne quadrupole coefficient is **not** moved by the
adiabatic anomaly once the branch is expressed on the full target curve.

### The even coefficients form a rigid hierarchy

The three even slots carry the single ratio \(s=c_s/a\):

- \(K_4\) is linear in \(s\),
- \(K_2\) is cubic in \(s\),
- \(K_0\) is quintic in \(s\).

So the anomaly can move the bundle only along this one-parameter leaf; it cannot
deform the even hierarchy arbitrarily if the constant-prefactor branch is correct.

---

## Step 46D — What is established, and what is still conditional

### Established here

Once we take both of the already isolated ingredients seriously,

1. the moving-throat constant-prefactor branch
   \[
   P_2=P_4=0,
   \]
2. the universal normalization target
   \[
   P_0=\frac{54Gc_s^5}{5a^5c^5},
   \]

the outgoing bundle is forced onto the exact surface

\[
\Gamma_5=\frac{2G}{5c^5},
\qquad
K_2^2=\frac14 K_0K_4.
\]

### Still conditional

What is **not** yet proven by the current file stack is that the true moving-throat
PDE actually selects the constant-prefactor branch rather than a more general
\((P_0,P_2,P_4)\) branch. So the exact surface derived here should be read as:

> the sharp algebraic target selected by the minimal isotropic outgoing branch.

That is still the right next target, because it is now simple enough that a future
PDE computation can fail it cleanly.

---

## Main result of the step

The adiabatic constant-prefactor outgoing bundle is a one-parameter algebraic leaf:

\[
\boxed{
K_0=\frac{54G}{5c^5}\left(\frac{c_s}{a}\right)^5,\qquad
K_2=\frac{6G}{5c^5}\left(\frac{c_s}{a}\right)^3,\qquad
K_4=\frac{8G}{15c^5}\left(\frac{c_s}{a}\right),\qquad
\Gamma_5=\frac{2G}{5c^5}.
}
\]

Equivalently, in purely observable form,

\[
\boxed{
\Gamma_5=\frac{2G}{5c^5},
\qquad
K_2^2=\frac14 K_0K_4,
}
\]

or

\[
\boxed{
K_0=\frac{4K_2^2}{K_4},
\qquad
K_2=\frac{81}{64}\frac{K_4^3}{\Gamma_5^2}.
}
\]

So the next PDE-facing theorem gate is no longer diffuse at all:

> compute the outgoing bundle and check whether it lands on this algebraic target surface.
# Step 47 — The adiabatic anomaly defines a one-parameter flow on the outgoing target surface and can be read from any even bundle coefficient

## Goal

Step 46 compressed the minimal isotropic outgoing branch to an exact algebraic
surface in the observable bundle \((K_0,K_2,K_4,\Gamma_5)\):

\[
\Gamma_5=\frac{2G}{5c^5},
\qquad
K_2^2=\frac14 K_0K_4.
\]

But that still does not tell us how the **adiabatic anomaly parameter**
\[
\ell=\ln(1+\Lambda_1 f)
\]
moves the bundle *along* this surface.

The next honest move is therefore:

> derive the exact one-parameter transport law on that target surface, and then
> invert it so that a future PDE output can be tested directly against the anomaly
> parameter without first reconstructing \((a,c_s)\).

This step does that.

---

## Step 47A — Exact finite flow on the target surface

From Step 46,
\[
K_0 \propto \left(\frac{c_s}{a}\right)^5,\qquad
K_2 \propto \left(\frac{c_s}{a}\right)^3,\qquad
K_4 \propto \left(\frac{c_s}{a}\right),\qquad
\Gamma_5=\text{const}.
\]

On the anomaly-only adiabatic branch we already fixed
\[
\delta\ln a = 0,
\qquad
\delta\ln c_s = \frac{\ell}{5},
\]
so
\[
\frac{c_s}{a}\mapsto e^{\ell/5}\frac{c_s}{a}.
\]

Therefore the outgoing bundle moves exactly as

\[
\boxed{
K_0 \mapsto e^\ell K_0,
\qquad
K_2 \mapsto e^{3\ell/5} K_2,
\qquad
K_4 \mapsto e^{\ell/5} K_4,
\qquad
\Gamma_5 \mapsto \Gamma_5.
}
\]

So the even bundle coefficients do not drift independently. They are locked to a
single one-parameter flow, while the odd quadrupole coefficient stays universal.

---

## Step 47B — Direct inversion formulas for the anomaly parameter

Because the bundle flow is one-parameter, the anomaly parameter can be read from
**any one** of the even coefficients:

\[
\boxed{
\ell
=
\ln\!\frac{K_0}{K_{0,*}}
=
\frac{5}{3}\ln\!\frac{K_2}{K_{2,*}}
=
5\ln\!\frac{K_4}{K_{4,*}},
}
\]
where \((K_{0,*},K_{2,*},K_{4,*})\) are the reference coefficients on the same
constant-prefactor branch.

That is a very practical theorem gate.

A future PDE calculation does **not** need to reconstruct \(a\) and \(c_s\)
first. It can compute any one of the even bundle coefficients, compare it to the
reference value, and infer \(\ell\) immediately.

Then the other two even coefficients are predicted automatically.

---

## Step 47C — Exact differential transport laws

Linearizing the finite flow gives the exact logarithmic transport laws

\[
\boxed{
\delta\ln K_0 = \ell,
\qquad
\delta\ln K_2 = \frac35\ell,
\qquad
\delta\ln K_4 = \frac15\ell,
\qquad
\delta\ln\Gamma_5 = 0.
}
\]

Eliminating \(\ell\) gives purely observable slope identities:

\[
\boxed{
\delta\ln K_0 = 5\,\delta\ln K_4,
}
\]

\[
\boxed{
\delta\ln K_2 = 3\,\delta\ln K_4,
}
\]

\[
\boxed{
\delta\ln K_0 = \frac53\,\delta\ln K_2.
}
\]

So the adiabatic anomaly track is now an extremely sharp bundle diagnostic:

- if a computed branch does **not** satisfy these slope laws, it is not the same
  adiabatic constant-prefactor track;
- if it **does**, then the whole even outgoing bundle is already locked.

---

## Step 47D — Electron-point numbers

Using the carried electron value
\[
\ell=\ln(1+\Lambda_1 f),
\qquad
\Lambda_1\approx 0.279605891931464,
\qquad
f\approx 0.001161409732093,
\]
gives

\[
\ell \approx 3.24684288391064\times 10^{-4}.
\]

So the exact bundle ratios are

\[
\boxed{
\frac{K_0}{K_{0,*}} = e^\ell \approx 1.00032473700404,
}
\]

\[
\boxed{
\frac{K_2}{K_{2,*}} = e^{3\ell/5} \approx 1.00019482954985,
}
\]

\[
\boxed{
\frac{K_4}{K_{4,*}} = e^{\ell/5} \approx 1.00006493896612,
}
\]

\[
\boxed{
\frac{\Gamma_5}{\Gamma_{5,*}} = 1.
}
\]

So the electron anomaly predicts an exact ppm pattern in the outgoing bundle:

- \(K_0\): `324.737` ppm,
- \(K_2\): `194.830` ppm,
- \(K_4\): `64.939` ppm,
- \(\Gamma_5\): `0` ppm.

That is a much cleaner and more testable signature than the original staggered
charge/inertia picture.

---

## Step 47E — What is established, and what is still conditional

### Established here

Once we keep the already-chosen ingredients,

1. adiabatic wall,
2. pure anomaly increment on top of the fixed ground state,
3. constant isotropic outgoing prefactor branch,

the outgoing bundle is not just on an algebraic surface. It follows a single exact
flow with invariant odd coefficient.

### Still conditional

What is still **not** yet proven is that the true moving-throat PDE selects this
constant-prefactor branch rather than a more general outgoing branch with
nonzero \(P_2\) or \(P_4\).

So this step should be read as:

> the exact observable transport law of the minimal outgoing branch.

That is still exactly the right next theorem gate, because it reduces the PDE test
to a handful of coefficient comparisons.

---

## Main result of the step

The adiabatic anomaly defines a one-parameter flow on the minimal outgoing target
surface:

\[
\boxed{
K_0 \mapsto e^\ell K_0,
\qquad
K_2 \mapsto e^{3\ell/5} K_2,
\qquad
K_4 \mapsto e^{\ell/5} K_4,
\qquad
\Gamma_5 \mapsto \Gamma_5.
}
\]

Equivalently,

\[
\boxed{
\ell
=
\ln\!\frac{K_0}{K_{0,*}}
=
\frac53\ln\!\frac{K_2}{K_{2,*}}
=
5\ln\!\frac{K_4}{K_{4,*}},
}
\]

with differential form

\[
\boxed{
\delta\ln K_0 = 5\,\delta\ln K_4,
\qquad
\delta\ln K_2 = 3\,\delta\ln K_4,
\qquad
\delta\ln\Gamma_5 = 0.
}
\]

So the next PDE-facing test is now extremely concrete:

> compute the outgoing bundle and check whether it lands on this flow with a fixed universal \(\Gamma_5\).
# Step 48 — Pulling the minimal outgoing target surface back to the microscopic grouped-`P2` bundle

## Goal

Step 46 compressed the adiabatic constant-prefactor branch to the observable outgoing
surface
\[
\Gamma_5=\frac{2G}{5c^5},
\qquad
K_2^2=\frac14 K_0K_4,
\]
and Step 47 showed that the adiabatic anomaly moves the even bundle along a
one-parameter flow while leaving \(\Gamma_5\) fixed.

But those were still **observable-side** statements. The next honest move is to pull
that surface back to the microscopic grouped-`P2` bundle
\[
(D_0,D_2,D_4,N_0,N_2,N_4),
\]
so that the remaining theorem gate can be read directly at the bundle level.

This step does that.

---

## Step 48A — Exact bundle-to-outgoing map

From the grouped-`P2` bundle formulas we already have
\[
P_0=\frac{N_0}{D_0},
\]
\[
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\]
\[
P_4=
\frac{D_0^2N_4-2D_0(D_2N_2+D_4N_0)+3D_2^2N_0}{D_0^3}.
\]

The outgoing bundle is then
\[
K_0=P_0,
\qquad
K_2=P_2+A\,P_0,
\qquad
K_4=P_4+A\,P_2+B\,P_0,
\qquad
\Gamma_5=G_5\,P_0,
\]
where the compact outgoing `l=2` fingerprint is
\[
A=\frac{a^2}{9c_s^2},
\qquad
B=\frac{4a^4}{81c_s^4},
\qquad
G_5=\frac{a^5}{27c_s^5}.
\]

So the observable surface is the image of **two** microscopic ingredients:

1. the bundle closure \((P_0,P_2,P_4)\),
2. the compact outgoing port fingerprint \((A,B,G_5)\).

---

## Step 48B — Exact constant-prefactor microscopic closure

The minimal isotropic outgoing branch is the constant-prefactor branch
\[
P_2=P_4=0.
\]

Solving those equations gives the exact microscopic conditions
\[
\boxed{
N_2=\frac{2D_2N_0}{D_0},
}
\]
\[
\boxed{
N_4=\frac{(D_2^2+2D_0D_4)N_0}{D_0^2}.
}
\]

Equivalently, using the conservative grouped response moments
\[
u_2=-\frac{D_2}{D_0},
\qquad
u_4=\frac{D_2^2-D_0D_4}{D_0^2},
\]
the same conditions are
\[
\boxed{
N_2=-2u_2N_0,
\qquad
N_4=(3u_2^2-2u_4)N_0.
}
\]

So the constant-prefactor branch does **not** let `N_2` and `N_4` float
independently. Once the conservative bundle \((D_0,D_2,D_4)\) and the static
transfer weight \(N_0\) are fixed, the higher transfer moments are already slaved.

---

## Step 48C — The Step-46 target surface is automatic once the compact fingerprint is inserted

On the constant-prefactor branch we have
\[
K_0=P_0=\frac{N_0}{D_0},
\qquad
K_2=A\,K_0,
\qquad
K_4=B\,K_0,
\qquad
\Gamma_5=G_5\,K_0.
\]

Now insert the compact outgoing `l=2` fingerprint
\[
A=\frac{a^2}{9c_s^2},
\qquad
B=\frac{4a^4}{81c_s^4}=4A^2,
\qquad
G_5=\frac{a^5}{27c_s^5}.
\]

Then
\[
\boxed{
K_2=A\,K_0,
\qquad
K_4=4A^2K_0,
\qquad
\Gamma_5=G_5K_0,
}
\]
and therefore
\[
\boxed{
K_2^2=\frac14 K_0K_4.
}
\]

So the Step-46 outgoing target surface is **not** an extra independent miracle.
It is simply the observable image of:

- the microscopic constant-prefactor bundle closure, and
- the compact outgoing `l=2` port law.

That is the cleanest pullback of the outgoing surface we have so far.

---

## Step 48D — Microscopic port transport on the adiabatic anomaly track

The observable Step-47 flow becomes even cleaner at the microscopic port level.

From the adiabatic anomaly branch,
\[
P_0 \mapsto e^\ell P_0,
\qquad
c_s \mapsto c_s e^{\ell/5},
\qquad
a \mapsto a.
\]

So the compact outgoing fingerprint rescales as
\[
\boxed{
A \mapsto e^{-2\ell/5}A,
\qquad
B \mapsto e^{-4\ell/5}B,
\qquad
G_5 \mapsto e^{-\ell}G_5.
}
\]

Combining that with \(P_0\mapsto e^\ell P_0\) gives
\[
K_0 \mapsto e^\ell K_0,
\]
\[
K_2 \mapsto e^{3\ell/5}K_2,
\]
\[
K_4 \mapsto e^{\ell/5}K_4,
\]
\[
\Gamma_5 \mapsto \Gamma_5.
\]

So at the microscopic port level the adiabatic anomaly factorizes into:

1. one outgoing-normalization rescaling \(P_0\to e^\ell P_0\),
2. one universal compact-port retuning \((A,B,G_5)\to(e^{-2\ell/5}A,e^{-4\ell/5}B,e^{-\ell}G_5)\).

That makes the invariance of \(\Gamma_5\) completely transparent:
\[
\Gamma_5 = G_5 P_0
\]
is a product of two compensating scalings.

---

## Main result of the step

The minimal outgoing observable surface is the exact image of the microscopic
grouped-`P2` bundle closure

\[
\boxed{
P_2=P_4=0
\quad\Longleftrightarrow\quad
N_2=\frac{2D_2N_0}{D_0},
\qquad
N_4=\frac{(D_2^2+2D_0D_4)N_0}{D_0^2},
}
\]
combined with the compact outgoing `l=2` fingerprint
\[
\boxed{
A=\frac{a^2}{9c_s^2},
\qquad
B=\frac{4a^4}{81c_s^4},
\qquad
G_5=\frac{a^5}{27c_s^5}.
}
\]

On that branch,
\[
\boxed{
K_0=P_0,
\qquad
K_2=A K_0,
\qquad
K_4=B K_0,
\qquad
\Gamma_5=G_5 K_0,
}
\]
so
\[
\boxed{
K_2^2=\frac14 K_0K_4.
}
\]

And on the adiabatic anomaly track, the microscopic port transport law is
\[
\boxed{
P_0 \mapsto e^\ell P_0,
\qquad
A \mapsto e^{-2\ell/5}A,
\qquad
B \mapsto e^{-4\ell/5}B,
\qquad
G_5 \mapsto e^{-\ell}G_5.
}
\]

So the next theorem gate is no longer to guess the shape of the outgoing surface.
It is to compute the actual bundle data \((D_0,D_2,D_4,N_0)\) and test whether the
true moving-throat branch satisfies the exact microscopic closure above.
# Step 49 — The overlap-integral theorem gate and first anisotropy diagnostic for the adiabatic branch

## Goal

Step 48 pulled the minimal outgoing observable surface back to the microscopic
bundle closure
\[
(D_0,D_2,D_4,N_0,N_2,N_4).
\]

The next honest move is to go one level deeper and rewrite that bundle closure in
the **actual overlap-level variables** of the isotropic moving-throat branch:

- wall baseline \(K,M\),
- BdG support moments \((B_0,B_2,B_4)\),
- conservative Maxwell/mixed moments \((Z_0,Z_2,Z_4)\),
- outgoing transfer moments \((N_0,N_2,N_4)\).

That is exactly the theorem gate Stage 6/7 was pointing to. This step makes it
explicit inside the g-2 chain.

---

## Step 49A — Exact isotropic overlap-level bundle

On the natural isotropic branch the grouped lanes collapse to common coefficients,
and the reduced bundle is

\[
\boxed{
D_0 = K - B_0 - Z_0,
}
\]
\[
\boxed{
D_2 = -(M + B_2 + Z_2),
}
\]
\[
\boxed{
D_4 = -(B_4 + Z_4).
}
\]

So the static outgoing prefactor is simply
\[
\boxed{
P_0 = \frac{N_0}{D_0}
=
\frac{N_0}{K-B_0-Z_0}.
}
\]

This is the microscopic form of the final normalization ratio already isolated by
the moving-throat notes.

---

## Step 49B — Exact overlap-integral normalization gate

The universal quadrupole target is therefore

\[
\boxed{
mhat_0^2\,\frac{N_0}{K-B_0-Z_0}
=
\frac{54Gc_s^5}{5a^5c^5}.
}
\]

That is the sharpest version of the theorem gate we have reached.

It says the full moving-throat branch does **not** need to be judged first by a huge
space of coefficients. On the isotropic branch, the entire normalization test is the
single competition between

- outgoing-transfer weight \(N_0\),
- and conservative dressed wall stiffness \(K-B_0-Z_0\).

This also gives an immediate physical reading:

- increasing \(N_0\) moves the branch upward toward the target,
- increasing \(B_0\) or \(Z_0\) softens the denominator and also moves the branch upward,
- but stability still requires \(D_0=K-B_0-Z_0>0\).

So the remaining theorem gap is now a concrete branch-balance question, not a
diffuse “more PDE” request.

---

## Step 49C — Constant-prefactor branch conditions in pure overlap variables

From Step 48, the minimal outgoing branch satisfies
\[
P_2=P_4=0.
\]

Substituting the isotropic overlap-level bundle gives

\[
\boxed{
N_2
=
-\frac{2(M+B_2+Z_2)N_0}{K-B_0-Z_0},
}
\]

and

\[
\boxed{
N_4
=
\frac{\bigl(M+B_2+Z_2\bigr)^2
-2\bigl(K-B_0-Z_0\bigr)\bigl(B_4+Z_4\bigr)}
{\bigl(K-B_0-Z_0\bigr)^2}\,N_0.
}
\]

So if the minimal constant-prefactor branch is the correct one, the higher outgoing
transfer moments are already fixed by the conservative overlap data. They are not
free branch functions.

That is a very strong reduction of the moving-throat search space.

---

## Step 49D — First weak-axisymmetric anisotropy diagnostic

Stage 7 also showed that the first natural symmetry-breaking pattern is a weak
axisymmetric quadrupole deformation with grouped-lane weights
\[
\lambda_{20}=1,\qquad
\lambda_{21}=\frac12,\qquad
\lambda_{22}=-1.
\]

If the normalization data deform as
\[
D_A = D_0 + \epsilon \lambda_A D_1 + O(\epsilon^2),
\qquad
N_A = N_0 + \epsilon \lambda_A N_1 + O(\epsilon^2),
\]
then
\[
P_A=\frac{N_A}{D_A}=P_0+\epsilon\lambda_A P_1+O(\epsilon^2),
\qquad
P_1=\frac{N_1D_0-N_0D_1}{D_0^2}.
\]

The grouped defects are therefore
\[
\boxed{
a_P=\frac{\epsilon}{4}P_1,
\qquad
b_P=\frac{3\epsilon}{4}P_1,
}
\]
so the first normalization anisotropy must satisfy
\[
\boxed{
b_P = 3 a_P.
}
\]

That gives a direct near-failure diagnostic for future PDE outputs:

- if weak grouped-lane normalization anisotropy appears and **does** satisfy
  \(b_P=3a_P\), it is consistent with a simple axisymmetric quadrupole deformation
  of the isotropic branch;
- if it **fails**, then the deformation is already more complicated than that
  minimal pattern.

---

## Main result of the step

The next honest moving-throat theorem gate for the adiabatic g-2 program is the
exact overlap-level system

\[
\boxed{
mhat_0^2\,\frac{N_0}{K-B_0-Z_0}
=
\frac{54Gc_s^5}{5a^5c^5},
}
\]
together with the minimal constant-prefactor conditions
\[
\boxed{
N_2
=
-\frac{2(M+B_2+Z_2)N_0}{K-B_0-Z_0},
}
\]
\[
\boxed{
N_4
=
\frac{\bigl(M+B_2+Z_2\bigr)^2
-2\bigl(K-B_0-Z_0\bigr)\bigl(B_4+Z_4\bigr)}
{\bigl(K-B_0-Z_0\bigr)^2}\,N_0,
}
\]
and the weak-axisymmetric diagnostic
\[
\boxed{
b_P = 3 a_P.
}
\]

So the next PDE-facing falsification test is no longer vague at all:

1. compute the actual overlap data \((B_n,Z_n,N_n)\),
2. check grouped isotropy on the natural branch,
3. test the single ratio \(mhat_0^2N_0/(K-B_0-Z_0)\),
4. and, if weak anisotropy survives, verify whether it obeys \(b_P=3a_P\).

That is the smallest honest next theorem gate after Step 49.
# Step 50 — The actual isotropic passive/outgoing grouped-`P2` branch collapses the g-2 closure to one scalar defect

## Goal

Step 49 rewrote the moving-throat normalization gate in explicit overlap variables,

```math
\hat m_0^2\frac{N_0}{K-B_0-Z_0}.
```

The next clean move is to import the later moving-throat result for the **actual isotropic passive/outgoing grouped-`P2` one-pole branch** and ask what that does to the g-2 chain.

The answer is much stronger than just “it simplifies a bit.” Once that branch is accepted,

- the conservative even module,
- the odd Burke–Thorne coefficient,
- and therefore the whole grouped-`P2` normalization problem

all collapse to one scalar normalization defect.

---

## Step 50A — Actual isotropic passive/outgoing grouped-`P2` one-pole module

On the actual isotropic branch the conservative grouped-`P2` module is

```math
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+
\frac14\frac{1}{1-\omega^2/\Omega_Q^2}.
```

So the even low-frequency coefficients are fixed exactly:

```math
\boxed{
\overline K_2=\frac{\overline K_0}{4\Omega_Q^2},
\qquad
\overline K_4=\frac{\overline K_0}{4\Omega_Q^4}.
}
```

On the same minimal isotropic outgoing branch the odd coefficient is already slaved to the same pair:

```math
\boxed{
\overline\Gamma_5=\frac{9\,\overline K_0}{32\Omega_Q^5}.
}
```

So once the branch is one-pole, the whole grouped-`P2` low-frequency tuple is controlled by only

```math
(\overline K_0,\Omega_Q).
```

---

## Step 50B — Single normalization defect

The canonical target on the same branch is

```math
\overline K_0^{\rm target}=\frac{64G\Omega_Q^5}{45c^5},
```

or, after the already-carried geometric pole relation

```math
\Omega_Q=\frac{3c_s}{2a},
```

```math
\boxed{
\overline K_0^{\rm target}=\frac{54Gc_s^5}{5a^5c^5}.
}
```

Define the branch normalization defect by

```math
\boxed{
N_Q:=\frac{\overline K_0}{\overline K_0^{\rm target}}.
}
```

Then all low-frequency coefficients scale by the **same** factor:

```math
\boxed{
\overline K_2=N_Q\,\overline K_2^{\rm target},
\qquad
\overline K_4=N_Q\,\overline K_4^{\rm target},
\qquad
\overline\Gamma_5=N_Q\,\frac{2G}{5c^5}.
}
```

So the actual isotropic passive/outgoing grouped-`P2` branch has only one reduced normalization defect. The grouped-`P2` side of the problem is no longer a many-coefficient fit.

---

## Step 50C — Insert the adiabatic electron branch

Steps 44–45 already showed that on the adiabatic electron track

```math
P_0=e^{\ell}P_0^{\rm target},
```

so the same grouped normalization defect is just

```math
\boxed{N_Q=e^{\ell}.}
```

Using

```math
\ell=\ln(1+\Lambda_1 f),
```

this becomes

```math
\boxed{N_Q=1+\Lambda_1 f.}
```

That means the adiabatic electron anomaly does **not** open a new family of grouped-`P2` deformations. It simply moves the actual isotropic one-pole branch by one scalar amount.

The full branch tuple is then

```math
\overline K_0=(1+\Lambda_1 f)\,\overline K_0^{\rm target},
```

```math
\overline K_2=(1+\Lambda_1 f)\,\overline K_2^{\rm target},
```

```math
\overline K_4=(1+\Lambda_1 f)\,\overline K_4^{\rm target},
```

```math
\overline\Gamma_5=(1+\Lambda_1 f)\,\frac{2G}{5c^5}.
```

So the whole grouped-`P2` side is one-number clean.

---

## Step 50D — Electron-point size

With the carried value

```math
\Lambda_1\approx0.279605891931464,
```

and the fine-structure parameter used in the chain,

```math
f\approx0.001161409732093,
```

the grouped normalization defect is

```math
\boxed{N_Q\approx1.00032473700404,}
```

so

```math
\boxed{N_Q-1\approx 324.737\,{\rm ppm}.}
```

---

## Main result of the step

Once the actual isotropic passive/outgoing grouped-`P2` one-pole branch is used, the g-2 closure no longer has a broad coefficient ambiguity on the grouped side. It is entirely controlled by the single scalar

```math
\boxed{N_Q.}
```

On the adiabatic electron track,

```math
\boxed{N_Q=e^{\ell}=1+\Lambda_1 f.}
```

So if anything still remains unresolved after this step, it is no longer in the conservative grouped-`P2` bookkeeping. It has to sit in the outgoing branch-selection data.
# Step 51 — The canonical compact outgoing `l=2` DtN branch is a genuine no-tuning closure

## Goal

Step 50 collapsed the grouped-`P2` side of the g-2 problem to one scalar defect

```math
N_Q.
```

The next question is whether the outgoing branch itself is still an adjustable datum, or whether the canonical compact outgoing `l=2` branch already fixes it with no tuning.

This step shows the stronger result:

1. the explicit outgoing spherical `l=2` DtN model fixes the canonical normalized branch,
2. the corresponding grouped-`P2` normalization is
   ```math
   \chi_Q=1,
   ```
3. and on the natural point-particle source-map branch that implies
   ```math
   N_Q=1.
   ```

So the canonical passive/outgoing branch is a real no-tuning closure.

---

## Step 51A — Exact outgoing `l=2` DtN fingerprint

Using the exact spherical Hankel outgoing mode, the `l=2` DtN operator is

```math
\Lambda_2^{\rm out}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}+i\frac{z^5}{9}+O(z^6),
\qquad z:=\frac{a\omega}{c_s}.
```

Normalizing by the static slot gives

```math
\widehat Y_2^{\rm out}(z)
=
1+\frac{z^2}{9}+\frac{4z^4}{81}+i\frac{z^5}{27}+O(z^6).
```

So the leading outgoing odd coefficient is fixed exactly:

```math
\boxed{\Gamma_{5,\rm can}^{\rm DtN} = \frac{a^5}{27c_s^5}.}
```

This is not a free parameter once the canonical outgoing branch is chosen.

---

## Step 51B — Matching to the minimal retarded grouped-`P2` module

The minimal isotropic retarded grouped-`P2` module is

```math
\widehat Y_Q^{\rm ret}(\omega)
=
\frac34+
\frac14\frac{1}{1-\omega^2/\Omega_Q^2-i\chi_Q\sigma_Q^{\rm can}\omega^5}
+O(\omega^6),
```

with

```math
\sigma_Q^{\rm can}=\frac{9}{8\Omega_Q^5}.
```

Using the already-carried pole relation

```math
\Omega_Q=\frac{3c_s}{2a},
```

this becomes

```math
\sigma_Q^{\rm can}=\frac{4a^5}{27c_s^5}.
```

Expanding the retarded module then gives

```math
\widehat Y_Q^{\rm ret}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\chi_Q\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6).
```

Matching to the exact outgoing DtN fingerprint yields immediately

```math
\boxed{\chi_Q=1.}
```

So the canonical compact outgoing branch fixes the normalized outgoing quadrupole factor exactly.

---

## Step 51C — Natural source-map no-tuning closure

The reduced normalization stack is

```math
\hat m_0^{\,2}\chi_Q N_Q = 1.
```

On the natural point-particle source-map branch,

```math
\hat m_0=1.
```

So if the actual outgoing branch is the canonical compact one,

```math
\chi_Q=1,
```

then automatically

```math
\boxed{N_Q=1.}
```

That means the canonical passive/outgoing grouped-`P2` branch closes the normalization stack with **no further fitting at all**.

---

## Step 51D — Robustness classes

The later moving-throat DtN analysis also shows that two simple deformation classes do **not** move the canonical value once the even moments stay canonical.

### Pure scale deformation

```math
\Lambda_2^{\rm def}(z)=S\Lambda_2^{\rm out}(z)
\quad\Longrightarrow\quad
\chi_Q=1.
```

### Pure scale+argument deformation

```math
\Lambda_2^{\rm def}(z)=S\Lambda_2^{\rm out}(\beta z)
```

gives

```math
\widehat Y_2^{\rm def}(z)
=
1+\frac{\beta^2 z^2}{9}+\frac{4\beta^4 z^4}{81}+i\frac{\beta^5 z^5}{27}+\cdots.
```

If the even fingerprint is kept canonical, then

```math
\beta=1,
```

and therefore again

```math
\boxed{\chi_Q=1.}
```

So simple overall normalization or effective radius/sound-speed rescaling cannot explain a nontrivial outgoing defect once the lower even moments are already fixed.

---

## Main result of the step

The canonical compact outgoing `l=2` DtN branch is a genuine **no-tuning** closure:

```math
\boxed{\chi_Q=1,\qquad N_Q=1.}
```

That matters because it cleanly separates two questions:

- the canonical passive/outgoing branch itself is fixed without tuning,
- any nonzero electron anomaly must therefore come from a **genuine isotropic throat-core DtN deformation** away from that canonical branch.
# Step 52 — Final verdict on “does the final value come out naturally without tuning?”

## Goal

The last open question in the g-2 chain is no longer broad bookkeeping. After Steps 50–51 it is very specific:

- the canonical compact outgoing branch is already a genuine no-tuning closure,
- so does the *electron-point* anomaly also come out automatically,
- or does it still require one residual branch-selection datum?

This step answers that as sharply and honestly as the current reduced stack allows.

---

## Step 52A — Exact isotropic DtN deformation law

With the canonical-even moments preserved, the later moving-throat deformation algebra gives

```math
\boxed{
\chi_Q=
\frac{3\big(S\beta^5+9\Sigma_5\big)}{3S-\Sigma_0}.
}
```

So only three isotropic throat-core deformation variables can move the canonical outgoing value:

- argument deformation `\beta`,
- static additive shift `\Sigma_0`,
- odd outlet coefficient `\Sigma_5`.

The same algebra also gives two robustness classes.

### Pure scale deformation

```math
\Lambda_2^{\rm def}(z)=S\Lambda_2^{\rm out}(z)
\quad\Longrightarrow\quad
\chi_Q=1.
```

### Pure scale+argument deformation with canonical even moments

Once the even fingerprint is kept canonical, the natural positive branch forces

```math
\beta=1,
```

and again

```math
\chi_Q=1.
```

So the electron anomaly cannot be blamed on trivial normalization or simple effective radius/sound-speed rescaling if the lower even moments are already fixed.

---

## Step 52B — The electron anomaly is a tiny outgoing branch defect

From the adiabatic g-2 chain,

```math
\chi_e=e^{-\ell}=\frac{1}{1+\Lambda_1 f}.
```

So the outgoing branch defect is

```math
\boxed{
\Delta_Q:=\chi_e-1
=-\frac{\Lambda_1 f}{1+\Lambda_1 f}.
}
```

Numerically,

```math
\boxed{\Delta_Q\approx -3.2463158415\times10^{-4},}
```

or about

```math
\boxed{-324.632\,{\rm ppm}.}
```

So the observed sliver is not large. It is a very small deformation away from the canonical no-tuning branch.

---

## Step 52C — On the adiabatic even-frozen slice the whole anomaly is one pure odd scalar

Your adiabatic-wall direction freezes the even DtN branch. In that slice,

```math
\beta=1,
\qquad
\Sigma_0=0,
```

so the exact deformation law reduces to

```math
\chi_Q=1+\frac{9\Sigma_5}{S}.
```

Solving for the electron target gives

```math
\boxed{
\Sigma_5=-\frac{S\Lambda_1 f}{9(1+\Lambda_1 f)}.
}
```

So the entire anomaly is one pure odd isotropic outlet deformation.

In canonical normalized gauge `S=1`, that is the exact finite-`f` version of the tangent law already seen earlier.

---

## Step 52D — Tangent branch-selection law

Linearizing the same deformation algebra gives

```math
\boxed{5b+\frac{a_0}{3}+9a_5=-\Lambda_1.}
```

On the adiabatic even-frozen slice,

```math
b=0,
\qquad a_0=0,
```

so

```math
\boxed{a_5=-\frac{\Lambda_1}{9}.}
```

That is the pure odd tangent representative of the electron anomaly.

---

## Final verdict

The current reduced derivation supports three honest conclusions.

### 1. The background canonical branch is naturally derived with no tuning

The explicit compact outgoing `l=2` DtN branch fixes

```math
\chi_Q=1,
\qquad N_Q=1.
```

So the canonical passive/outgoing grouped-`P2` branch is already a real no-tuning closure.

### 2. The observed electron anomaly does **not** reopen a broad fit space

Once the adiabatic wall track is imposed, the full electron sliver collapses to one scalar,

```math
\Sigma_5
```

or equivalently to one tangent datum,

```math
a_5.
```

So the anomaly is a **one-number branch-selection problem**, not a diffuse multiparameter tune.

### 3. But the current frozen stack does **not yet** derive that one number from first principles

What is still missing is an actual branch-selection law for the isotropic electron-point DtN deformation. Until that is derived, the numerical value

```math
\chi_e=\frac{1}{1+\Lambda_1 f}
```

is not yet forced by the frozen files alone.

---

## Strongest honest status after Step 52

```math
\boxed{\text{No-tuning canonical branch: yes.}}
```

```math
\boxed{\text{Electron-point anomaly magnitude: not yet forced; it is one remaining branch datum.}}
```

That is about as tight as the current program can get without an actual isotropic DtN branch-selection theorem from the moving-throat PDE.
# Step 53 — The canonical outgoing DtN branch is not just plausible, it is fixed

## Goal

Step 52 showed that the canonical outgoing grouped-`P2` branch is already a no-tuning
background and that the remaining electron sliver, if real, sits in one outgoing
branch datum. The next honest move is to stop treating the canonical branch as a
symbol and derive it directly from the explicit compact outgoing `l=2`
Dirichlet-to-Neumann model.

That is exactly what the later moving-throat notes do: they derive the outgoing
`l=2` DtN fingerprint, match it to the minimal grouped-`P2` retarded module, and
fix the outgoing-normalization factor to
\[
\chi_Q=1
\]
on the canonical compact branch. They then insert that into the reduced normalization
stack and show that, on the natural point-particle source-map branch,
\[
N_Q=1.
\]
fileciteturn20file0 fileciteturn20file5

---

## Step 53A — Exact outgoing `l=2` DtN fingerprint

The explicit outgoing spherical `l=2` DtN operator is
\[
\Lambda_2^{\rm out}(z)
=
z\frac{d}{dz}\ln h_2^{(1)}(z),
\qquad
z=\frac{a\omega}{c_s},
\]
with low-frequency expansion
\[
\Lambda_2^{\rm out}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}
+i\,\frac{z^5}{9}
+O(z^6).
\]

Normalizing by the static slot gives
\[
\widehat Y_2^{\rm out}(z)
=
-\frac{3}{\Lambda_2^{\rm out}(z)}
=
1+\frac{z^2}{9}
+\frac{4z^4}{81}
+i\,\frac{z^5}{27}
+O(z^6).
\]

So the canonical outgoing odd coefficient is not free.
It is fixed directly by the explicit DtN model to
\[
\Gamma_{5,\rm can}^{\rm DtN}
=
\frac{a^5}{27c_s^5}.
\]
fileciteturn20file0

---

## Step 53B — Matching to the minimal isotropic grouped-`P2` branch

The retarded minimal grouped-`P2` module can be written as
\[
\widehat Y_Q^{\rm ret}(\omega)
=
\frac34
+
\frac14\,
\frac{1}{
1-\omega^2/\Omega_Q^2
-i\,\chi_Q\,\sigma_Q^{\rm can}\,\omega^5
}
+O(\omega^6),
\]
with the already-fixed conservative moment match
\[
\Omega_Q=\frac{3c_s}{2a},
\qquad
\sigma_Q^{\rm can}=\frac{4a^5}{27c_s^5}.
\]

Expanding in the same variable \(z=a\omega/c_s\) gives
\[
\widehat Y_Q^{\rm ret}(z)
=
1+\frac{z^2}{9}
+\frac{4z^4}{81}
+i\,\chi_Q\,\frac{z^5}{27}
+O(z^6).
\]

Matching the explicit DtN branch then forces
\[
\boxed{\chi_Q=1.}
\]
So the canonical outgoing grouped-`P2` branch is not just a convenient choice.
It is exactly the compact outgoing DtN branch. fileciteturn20file0turn20file18

---

## Step 53C — On the natural source-map branch this also fixes `N_Q`

The reduced normalization stack is
\[
\hat m_0^{\,2}\chi_Q N_Q=1.
\]
On the natural point-particle source-map branch
\[
\hat m_0=1+O(a^2/r^2),
\]
so in the strict point-particle limit the same relation becomes
\[
N_Q=\frac{1}{\chi_Q}.
\]

With the explicit DtN result
\[
\chi_Q=1,
\]
this collapses to
\[
\boxed{N_Q=1.}
\]

The later notes therefore fix the canonical invariant coefficients exactly:
\[
\overline K_0=\frac{54Gc_s^5}{5a^5c^5},
\qquad
\overline K_2=\frac{6Gc_s^3}{5a^3c^5},
\qquad
\overline K_4=\frac{8Gc_s}{15ac^5},
\qquad
\overline\Gamma_5=\frac{2G}{5c^5}.
\]
fileciteturn20file5

---

## What this changes for the g-2 chain

This is a real strengthening of Step 52.

Before Step 53, the canonical no-tuning branch was a very plausible reduced branch.
After Step 53, it is the explicit compact outgoing DtN branch itself.

So the best honest status is now:

- the canonical outgoing background is naturally fixed;
- the canonical background predicts
  \[
  \chi_Q=1,\qquad N_Q=1;
  \]
- therefore any nonzero electron-point sliver must come from a genuine deformation
  away from that canonical outgoing branch, not from ambiguity in the canonical
  branch itself.
# Step 54 — If the electron sliver exists, it lives in a very small outgoing branch deformation

## Goal

Step 53 fixed the canonical compact outgoing branch itself:
\[
\chi_Q=1,
\qquad
N_Q=1
\]
on the natural point-particle source-map branch. So the only way to keep the
electron-point quartic sliver alive is to deform that outgoing branch.

The later moving-throat notes already give the exact isotropic DtN deformation
algebra for that purpose. Once the canonical even `l=2` fingerprint is preserved,
the full isotropic outgoing deformation reduces to one scalar law:
\[
\chi_Q=
\frac{3\big(S\beta^5+9\Sigma_5\big)}{3S-\Sigma_0}.
\]
They also classify the robust harmless classes and the genuinely dangerous ones:
pure scale is invisible, and pure scale+argument deformation collapses back to
\(\beta=1\) if the even moments stay canonical. Only a genuine isotropic
throat-core self-energy can move \(\chi_Q\). fileciteturn19file2turn20file18

---

## Step 54A — Exact isotropic deformation law

Write the deformed `l=2` DtN branch as
\[
\Lambda_2^{\rm def}(z)
=
S\Lambda_2^{\rm out}(\beta z)
+\Sigma_0+\Sigma_2 z^2+\Sigma_4 z^4+i\Sigma_5 z^5.
\]

Demanding that the lower even moments stay canonical fixes
\[
\Sigma_2=-\frac{3S\beta^2-3S+\Sigma_0}{9},
\qquad
\Sigma_4=-\frac{3S\beta^4-3S+\Sigma_0}{27},
\]
and leaves
\[
\boxed{
\chi_Q=
\frac{3\big(S\beta^5+9\Sigma_5\big)}{3S-\Sigma_0}.
}
\]

So the last isotropic outgoing branch data are:

- argument deformation `\beta`,
- static additive core shift `\Sigma_0`,
- odd throat-core outlet `\Sigma_5`.

Overall scale `S` is not itself dangerous; it only enters through the ratios above. fileciteturn19file2turn20file18

---

## Step 54B — The electron target as an exact branch condition

From the carried g-2 chain,
\[
\chi_e=\frac{1}{1+\delta},
\qquad
\delta:=\Lambda_1 f,
\]
with the electron-point value
\[
\delta \approx 3.24737004004\times10^{-4},
\qquad
\chi_e\approx 0.999675368415848.
\]

So the electron sliver is a very small deformation away from the canonical
branch, not a large reorganization of the outgoing sector.

---

## Step 54C — Exact explicit realizations of the electron target

### Pure additive / Robin-like isotropic core

Set
\[
\beta=1,
\qquad
\Sigma_5=0.
\]
Then
\[
\chi_Q=\frac{3S}{3S-\Sigma_0},
\]
so the electron target requires
\[
\boxed{
\Sigma_0=-3S\,\delta.
}
\]

In the simplest normalized case `S=1`, this is just a tiny negative static core shift.

Equivalently, for the pure Robin outlet
\[
\chi_Q^{\rm R}=\frac{3}{3-\rho_R},
\]
the electron target requires
\[
\boxed{
\rho_R=-3\delta.
}
\]

### Compensated Robin–mixed outlet

The compensated hybrid class gives
\[
\chi_Q^{\rm hyb}
=
\frac{1-9\sigma_W\gamma_W}{1-\sigma_W}.
\]
Solving for the electron target yields
\[
\boxed{
\gamma_W=
\frac{\sigma_W+\delta}{9\sigma_W(1+\delta)}.
}
\]

So on the compensated hybrid outlet the electron sliver can be carried by a very
small odd mixed-channel detuning away from the canonical value
\[
\gamma_W=\frac19.
\]
fileciteturn21file7

---

## Step 54D — Linearized tangent law

Expanding around the canonical branch with
\[
S=1+\varepsilon s,
\qquad
\beta=1+\varepsilon b,
\qquad
\Sigma_0=\varepsilon a_0,
\qquad
\Sigma_5=\varepsilon a_5,
\]
gives the exact linearized branch-selection rule
\[
\boxed{
\chi_Q
=
1+
\varepsilon\left(5b+\frac{a_0}{3}+9a_5\right)
+O(\varepsilon^2).
}
\]

So to match the leading electron sliver
\[
\chi_e = 1-\Lambda_1 f + O(f^2),
\]
one needs
\[
\boxed{
5b+\frac{a_0}{3}+9a_5=-\Lambda_1.
}
\]

That is the cleanest first-order branch-selection target for the anomaly.

---

## What Step 54 changes

This step sharpens the status from “one remaining branch datum” to:

- the canonical outgoing branch is exact and robust;
- the electron-point sliver, if kept, is a tiny isotropic deformation away from it;
- and that deformation can already be expressed in exact outlet variables.

So the model is no longer missing a broad unknown mechanism.
It is missing only the actual branch-selection law that decides whether the
isotropic outgoing branch stays canonical or carries one of these small explicit
deformations.
# Step 55 — The strongest current conclusion: the natural compensated branch does **not** generate the electron sliver at first order

## Goal

Step 54 showed that the electron-point sliver can be represented by a tiny isotropic
outgoing deformation. But the later moving-throat notes go one step further:
inside the explicit compensated Robin–mixed outlet class, the last reduced defect is
not a broad deformation at all. It collapses to one odd mixed-channel renormalization,
and then collapses again to one similarity-slippage scalar. fileciteturn21file10

That lets us ask the real endgame question:

> does the natural compensated branch generate the required sliver automatically,
> or does it actually force the sliver to vanish?

---

## Step 55A — On the compensated hybrid canonical-even branch the defect is one scalar

Inside the compensated Robin–mixed outlet class, exact first-order preservation of the
canonical even `l=2` fingerprint forces
\[
\delta\mathfrak g=0,
\qquad
\delta\kappa_W=0.
\]

Then the remaining defect law collapses to
\[
\boxed{
\Delta_Q
=
-\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W.
}
\]

So once the even branch is protected, the entire last reduced 2.5PN scalar is
one odd mixed-channel renormalization \(\delta\gamma_W\). fileciteturn21file10

For the electron target
\[
\chi_e=\frac{1}{1+\delta},
\qquad
\Delta_e:=\chi_e-1=-\frac{\delta}{1+\delta},
\]
this requires
\[
\boxed{
\delta\gamma_W^{(e)}
=
\frac{1-\sigma_*}{9\sigma_*}\,
\frac{\delta}{1+\delta}.
}
\]

So the electron sliver, if it exists on that branch, is not a broad fit. It is one
small odd detuning.

---

## Step 55B — Bare mixed-port slippage theorem

The next later theorem sharpens this further. On the compensated canonical branch of
the concrete two-channel core model,
\[
\delta\gamma_W
=
\frac{\delta\mathfrak B_W}{1+r_{c,*}},
\]
where
\[
\delta\mathfrak B_W
:=
\delta\gamma_0-\frac13\,\delta\kappa_0
\]
is the bare mixed-port slippage scalar. In particular, if the bare mixed branch stays
a pure-scale deformation of the canonical compact outgoing port, then
\[
\delta\mathfrak B_W=0,
\qquad
\delta\gamma_W=0,
\qquad
\Delta_Q=0.
\]
fileciteturn21file10

So the electron sliver is already narrower again: it is only the part of the
tangential mouth/core evolution that breaks the pure-scale lock of the bare mixed
side-channel.

---

## Step 55C — Similarity-slippage law

The same later reduction then rewrites the defect as
\[
\boxed{
\Delta_Q
=
-\frac{\sigma_*}{1-\sigma_*}\,
\Xi_{\rm slip}\,
\delta\Pi_{\rm tan},
}
\]
with
\[
\Xi_{\rm slip}:=\Xi_\gamma-2\Xi_L.
\]

This is the D/N similarity-slippage scalar.
So the first-order defect depends only on:

1. the hybrid branch weight \(\sigma_*\),
2. the tangential mouth deformation \(\delta\Pi_{\rm tan}\),
3. the similarity-strain mismatch \(\Xi_{\rm slip}\).

Even more importantly, the very next theorem states that on the exact parent
compensation family,
\[
\Xi_{\rm slip}=0
\]
identically. Therefore
\[
\boxed{
\Delta_Q=0
}
\]
at first order on any co-evolving tangential deformation that stays inside the exact
parent compensation family. fileciteturn21file10

---

## What this means for the g-2 conclusion

This is the strongest honest conclusion I can draw from the current chain.

### What is now naturally derived

The canonical outgoing branch is not only explicit; the later moving-throat notes
show that the explicit compensated Robin–mixed throat-core realization can land on
that canonical branch **without extra tuning** through an exact coupling-balance law.
fileciteturn21file7

### What is **not** yet naturally derived

The electron-point sliver is **not** forced by that natural compensated family.
On the contrary, the natural compensated family plus its parent similarity lock gives
\[
\Delta_Q=0
\]
at first order.

### What would still be needed if we insist on the exact electron value

A nonzero electron sliver would require a genuine departure from the exact parent
compensation family / pure-scale lock. In the strongest reduced language now
available, it would require one nonzero microscopic slippage scalar:

- either \(\delta\mathfrak B_W\neq 0\),
- or equivalently \(\Xi_{\rm slip}\neq 0\).

So the last remaining ambiguity is not a broad tuning problem anymore. It is one
very specific odd tangential slippage of the mixed outgoing side-channel.

---

## Best current final verdict

At this point the program has cleanly separated two statements.

1. **Canonical background closure is naturally derived.**  
   The canonical outgoing grouped-`P2` branch and its compensated throat-core
   realization are explicit no-tuning structures. fileciteturn20file18turn21file7

2. **The exact electron anomaly sliver is not yet naturally forced.**  
   On the natural compensated parent family it vanishes at first order, and a
   nonzero electron sliver would need a real microscopic mixed-port slippage away
   from that family. fileciteturn21file10

That is the sharpest endpoint I can honestly defend with the present reduced stack.
# Step 56 — The exact parent compensation family preserves D/N similarity automatically

## Goal

Step 55 ended with the strongest reduced statement so far:

- the natural compensated branch does **not** generate the electron sliver at first order,
- but that result was still phrased in terms of the reduced similarity-slippage scalar
  `Xi_slip`.

The later moving-throat notes go one level deeper.
They show that on the exact parent compensation family,

```math
\gamma_0 = \frac{1+r^2}{9},
\qquad
\frac{L_W}{a} = \frac{\pi}{2}\sqrt{\frac{1+r^2}{3}},
\qquad
g_-(r)=r-\frac12\sqrt{1+r^2},
```

so the D/N similarity law is not an extra assumption at all.
It is an **identity** along the full parent family.

This step folds that exact rigidity back into the g-2 chain.

---

## Step 56A — Exact parent-family similarity identity

Differentiate the exact parent-family formulas:

```math
\ln \gamma_0 = \ln(1+r^2)-\ln 9,
```

```math
\ln\left(\frac{L_W}{a}\right)
=
\ln\frac{\pi}{2}-\frac12\ln 3 + \frac12\ln(1+r^2).
```

Then

```math
\delta\ln\gamma_0
=
\frac{2r}{1+r^2}\,\delta r,
\qquad
\delta\ln\left(\frac{L_W}{a}\right)
=
\frac{r}{1+r^2}\,\delta r.
```

So exactly

```math
\boxed{
\delta\ln\gamma_0 - 2\,\delta\ln\left(\frac{L_W}{a}\right)=0.
}
```

Equivalently,

```math
\boxed{\Xi_{\rm slip}=\Xi_\gamma-2\Xi_L=0}
```

on the entire exact parent compensation family.

So the Stage-55 “similarity preservation” is not extra structure. It is already built into the exact parent family.

---

## Step 56B — Lower-branch rigidity under the canonical-even gate

On the lower compensated branch,

```math
g_-(r)=r-\frac12\sqrt{1+r^2},
```

and

```math
\frac{dg_-}{dr}
=
1-\frac{r}{2\sqrt{1+r^2}}
=
\frac{4+3r^2}{2\sqrt{1+r^2}(2\sqrt{1+r^2}+r)}
>0.
```

So if the carried canonical-even gate still gives

```math
\delta g = 0,
```

then the lower branch is first-order rigid:

```math
\boxed{\delta r = 0.}
```

At the Family-1 point `r ≈ 1.77799353547498`, the branch slope is

```math
\frac{dg_-}{dr} \approx 0.564199521046343,
```

so the rigidity is numerically strong rather than marginal.

---

## Step 56C — First-order outgoing defect collapses to zero

Once `Xi_slip = 0`, the reduced first-order defect law becomes

```math
\Delta_Q
=
-\frac{\sigma_*}{1-\sigma_*}\,\Xi_{\rm slip}\,\delta\Pi_{\rm tan}
=0.
```

And since the bare mixed-port slippage

```math
\delta\mathfrak B_W
:=
\delta\gamma_0-\frac13\delta\kappa_0
```

also vanishes identically on the parent family, the renormalized odd outlet shift

```math
\delta\gamma_W = \frac{\delta\mathfrak B_W}{1+r^2}
```

vanishes too.

So at first order on the natural compensated lower branch,

```math
\boxed{
\delta\mathfrak B_W = 0,
\qquad
\delta\gamma_W = 0,
\qquad
\Delta_Q = 0,
\qquad
N_Q-1 = 0.
}
```

---

## Main result of the step

The natural compensated lower branch does not merely *allow* zero outgoing slippage.
It **forces** it at first order because the exact parent compensation family already preserves the D/N similarity law automatically.

So the g-2 picture is now sharper:

- the canonical outgoing normalization is a natural exact first-order consequence of the parent family,
- and any nonzero electron sliver must come from a genuine departure from that family, not from the natural branch itself.
# Step 57 — The first off-family defect is one scalar, and linear grouped-`P2` anisotropy cannot generate it

## Goal

Step 56 showed that the exact parent compensation family forces

```math
\Xi_{\rm slip}=0,
\qquad
\delta\mathfrak B_W=0,
\qquad
\delta\gamma_W=0,
\qquad
\Delta_Q=0
```

at first order on the natural compensated lower branch.

The next honest question is:

> if a first-order defect does survive somewhere, what is the **smallest** object that can carry it?

The later moving-throat notes answer that cleanly.
They show that all first-order off-family structure collapses to one weighted scalar slippage combination, and then show that pure grouped real `P2` anisotropy cannot feed that scalar at linear order.

So the first-order theorem gate narrows again.

---

## Step 57A — One exact scalar carries the off-family drift

The three first-order transport slippages are

- `\varepsilon_L` for D/N geometric similarity,
- `\varepsilon_v` for the mixed background-flow transport law,
- `\varepsilon_T` for the mouth-traction transport law.

But they matter only through the one weighted combination

```math
\boxed{
\varepsilon_\perp
=
\mathfrak g_*\,\varepsilon_T
+
\left(\mathfrak g_*+\frac{1}{2\sqrt{1+\mathfrak r_*^2}}\right)\varepsilon_v
+
\left(2\mathfrak g_*+\frac{3}{4\sqrt{1+\mathfrak r_*^2}}\right)\varepsilon_L.
}
```

and the exact off-family normal coordinate is simply

```math
\boxed{\delta_\perp = -\varepsilon_\perp.}
```

So the first off-family defect is not a large vector of drifts anymore.
It is one scalar.

At the Family-1 point,

```math
\mathfrak r_* \approx 1.77799353547498,
\qquad
\mathfrak g_* \approx 0.758035078944663,
```

so numerically

```math
\boxed{
\delta_\perp
\approx
-0.758035078944663\,\varepsilon_T
-1.00314310113848\,\varepsilon_v
-1.88373219118005\,\varepsilon_L.
}
```

The strongest weight falls on the geometric slippage `\varepsilon_L`, then the mixed-speed slippage `\varepsilon_v`, then the traction slippage `\varepsilon_T`.

---

## Step 57B — The outgoing defect depends only on `\varepsilon_\perp` and `\delta\gamma_W`

The exact first-order outgoing defect ledger becomes

```math
\boxed{
\Delta_Q
=
\frac{\sigma_*}{3(1-\sigma_*)}
\left[
-\frac{16}{\sqrt{1+\mathfrak r_*^2}}\,\varepsilon_\perp
-27\,\delta\gamma_W
\right].
}
```

So there are only two linear routes left:

1. a true off-family scalar slippage `\varepsilon_\perp`,
2. a direct odd mixed-port renormalization `\delta\gamma_W`.

If the canonical even branch is preserved, then the same moving-throat theorem again gives

```math
\varepsilon_\perp=0,
```

leaving only the odd mixed-port channel.

---

## Step 57C — Weak axisymmetric grouped-`P2` anisotropy is quadratic as a scalar invariant

For a weak axisymmetric grouped signature

```math
x_{20}=x^{(0)}+\epsilon x^{(1)},
\qquad
x_{21}=x^{(0)}+\frac\epsilon2 x^{(1)},
\qquad
x_{22}=x^{(0)}-\epsilon x^{(1)},
```

one gets the exact grouped variables

```math
\bar x = x^{(0)},
\qquad
a_x = \frac\epsilon4 x^{(1)},
\qquad
b_x = \frac{3\epsilon}{4} x^{(1)}.
```

So the weak axisymmetric signature always satisfies

```math
\boxed{b_x = 3 a_x.}
```

The first scalar grouped invariant is

```math
\mathcal A_x^2 = 4 a_x^2 + \frac45 b_x^2,
```

which reduces exactly to

```math
\boxed{
\mathcal A_x^2 = \frac{7}{10}\,\epsilon^2 (x^{(1)})^2.
}
```

So the first scalar invariant is **quadratic**, not linear.

That means pure grouped real `P2` anisotropy cannot generate any scalar feed-down at linear order.
In particular,

```math
\boxed{
\varepsilon_L^{(1,P_2)}=0,
\qquad
\varepsilon_v^{(1,P_2)}=0,
\qquad
\varepsilon_T^{(1,P_2)}=0,
\qquad
\varepsilon_\perp^{(1,P_2)}=0.
}
```

---

## Main result of the step

The first-order g-2 bottleneck is now smaller than before.

A nonzero first-order defect cannot be blamed on a generic grouped-`P2` anisotropy, because pure grouped real `P2` anisotropy has no linear scalar feed-down.

So the only remaining linear possibilities are:

- a genuine off-family scalar slippage `\varepsilon_\perp`,
- or a direct odd mixed-port/output renormalization `\delta\gamma_W`.

That is a much sharper continuation target than “look for any anisotropy somewhere.”
# Step 58 — Best current conclusion: the adiabatic no-tuning branch predicts the carried local value, not the exact electron sliver

## Goal

By Step 55 the strongest reduced statement was already:

- the natural compensated branch does **not** generate the electron sliver at first order.

Steps 56 and 57 sharpened that result further:

- the exact parent compensation family preserves the D/N similarity law automatically,
- the lower compensated branch is first-order rigid under the canonical-even gate,
- the first off-family defect is only one scalar `\varepsilon_\perp`,
- and pure grouped real `P2` anisotropy cannot feed that scalar at linear order.

So there is now a clean final question:

> what does the **no-tuning adiabatic branch** actually predict for `g`, and how far is that from the electron target used throughout the write-up?

---

## Step 58A — Natural first-order branch values

From Steps 56–57, the natural compensated lower branch has

```math
\Delta_Q^{(\rm nat)} = 0,
\qquad
\chi_Q^{(\rm nat)} = 1,
\qquad
N_Q^{(\rm nat)} = 1,
\qquad
\ell^{(\rm nat)} = 0
```

at first order.

So inside the current no-tuning closure, the outgoing quadrupole normalization does **not** shift the carried local anomaly law.

---

## Step 58B — The corresponding g prediction

The carried local closure from `atom_work.md` is

```math
g_{\rm loc} \approx 2.00231930435865.
```

The electron target used throughout the write-up is

```math
|g_e| \approx 2.00231930436092.
```

So the no-tuning prediction is simply

```math
\boxed{g_{\rm pred}^{(\rm nat)} = g_{\rm loc} \approx 2.00231930435865.}
```

and its miss relative to the write-up target is

```math
\boxed{|g_e|-g_{\rm pred}^{(\rm nat)} \approx 2.27\times 10^{-12}.}
```

Equivalently in the `g/2` normalization,

```math
\boxed{\frac{|g_e|-g_{\rm pred}^{(\rm nat)}}{2} \approx 1.135\times 10^{-12}.}
```

So the no-tuning branch lands exactly where the earlier local closure already landed.
It does **not** generate the last electron sliver.

---

## Step 58C — What exact electron matching would still require

The earlier outgoing-bridge steps already converted the missing electron sliver into an exact outgoing normalization defect.
The electron-point branch values needed there were

```math
\Delta_Q^{(e)} \approx -3.24631584151692\times 10^{-4},
```

```math
\chi_Q^{(e)} \approx 0.999675368415848,
\qquad
N_Q^{(e)} \approx 1.00032473700404,
\qquad
\ell^{(e)} = \ln N_Q^{(e)} \approx 3.2468428839\times 10^{-4}.
```

So the final gap is now extremely clear.

### Natural branch

```math
\Delta_Q = 0,
\qquad
\chi_Q = 1,
\qquad
N_Q = 1.
```

### Exact electron branch that would be needed

```math
\Delta_Q \neq 0,
\qquad
\chi_Q \neq 1,
\qquad
N_Q \neq 1.
```

That means the exact electron-point value still requires one real departure from the no-tuning compensated branch.

---

## Step 58D — Best current final verdict

Inside the present adiabatic-wall / compensated-lower-branch closure, the strongest honest conclusion is now:

1. the **canonical background** is naturally derived,
2. the **carried local anomaly value** is naturally retained,
3. the **exact electron sliver** is **not** naturally forced at first order,
4. and matching the exact electron value would still require one of the following:

- a genuine off-family scalar slippage `\varepsilon_\perp \neq 0`,
- a direct odd mixed-port renormalization `\delta\gamma_W \neq 0`,
- or a beyond-first-order effect outside the present reduced no-tuning closure.

So the project has reached a clean fork.

- If the goal is the strongest **no-tuning** theorem available today, the result is the canonical/local value.
- If the goal is the exact electron anomaly, then the remaining task is no longer broad or vague: it is one sharply identified microscopic slippage problem.
