
# Moving-Throat PDE — Stage 166: Exact Triangular Normal Form of the Coherent Weak-Axisymmetric Defect

## Purpose

Stage 165 reduced the coherent weak-axisymmetric grouped defect to five microscopic slippages
\[
(\Sigma_Z,\Sigma_\chi,\Sigma_\epsilon,\Sigma_\delta,\Sigma_\eta),
\]
with the tracking factor already carried by the single combination
\[
\Sigma_{\rm tr}:=(1+\chi_0)\Sigma_\delta+(1+\delta_U)\Sigma_\chi.
\]

That was already a strong narrowing, but the bookkeeping was still not in its smallest useful form.
The next honest step is now very sharp:

> reorganize the Stage-165 slippages into the smallest branch-adapted coordinates that directly control the three reduced observable drifts
> \[
> \Theta_1,\qquad \Xi_1,\qquad \mathcal R_1,
> \]
> and determine the exact inverse map.

That is what this stage does.

The main result is that the full coherent weak-axisymmetric problem collapses to three branch-adapted scalars:
\[
\Sigma_{\rm tr},
\qquad
\Sigma_{\rm nt},
\qquad
\Sigma_\eta,
\]
where \(\Sigma_{\rm nt}\) is the genuine **nontracking transfer-shape slippage**
\[
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
\]

In terms of these three scalars, the reduced observable drifts take the exact triangular form
\[
\boxed{
\Theta_1
=
-\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}\,\Sigma_{\rm tr},
}
\]
\[
\boxed{
\Xi_1
=
\frac{2\chi_0}{(1+\chi_0)(1+\delta_U)}\,\Sigma_{\rm tr}
+\Sigma_{\rm nt},
}
\]
\[
\boxed{
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_\eta}{1-\epsilon_\eta}\,\Sigma_\eta.
}
\]

So the grouped defect problem is no longer a five-slippage bookkeeping problem.
It is an exact three-scalar normal form:

- \(\Sigma_{\rm tr}\) controls tracking drift,
- \(\Sigma_{\rm nt}\) controls the nontracking transfer-shape defect,
- \(\Sigma_\eta\) controls the residual dressing mismatch between the direct transfer shape and the selected-branch demand ratio.

That is the strongest exact compression reached so far.

---

## 1. Branch-adapted defect coordinates

Stage 165 already gave the exact tracking/nontracking split
\[
\Xi_1
=
\Sigma_Z
+\frac{2\chi_0}{(1+\chi_0)(1+\delta_U)}\,\Sigma_{\rm tr}
+\frac{2\epsilon_W}{1-\epsilon}\frac{11+9\delta_U}{11(1+\delta_U)}\,\Sigma_\epsilon
-\left[
\frac{2\chi_0}{1+\delta_U}
+\frac{4\epsilon_W\delta_U}{11(1-\epsilon)(1+\delta_U)^2}
\right]\Sigma_\delta.
\]
This suggests defining the genuine nontracking transfer-shape slippage
\[
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
\]
Then the tracking piece and the nontracking piece are separated exactly:
\[
\boxed{
\Xi_1
=
A_{\rm tr}\,\Sigma_{\rm tr}
+\Sigma_{\rm nt},
\qquad
A_{\rm tr}:=\frac{2\chi_0}{(1+\chi_0)(1+\delta_U)}.
}
\]

So \(\Sigma_{\rm nt}\) is the first branch-adapted scalar beyond \(\Sigma_{\rm tr}\).

The third and final branch-adapted scalar is simply the selected-branch dressing slippage already identified in Stage 165:
\[
\boxed{
\Sigma_\eta:=2c_1-\kappa_U-\kappa_\eta=\frac{\eta_1}{\epsilon_\eta}.
}
\]

These are therefore the natural branch-adapted coordinates:
\[
\boxed{
(\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta).
}
\]

---

## 2. Exact triangular normal form

Stage 165 already gave
\[
\Theta_1
=
-\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}\,\Sigma_{\rm tr}.
\]
Define
\[
\boxed{
C_{\rm tr}:=
\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}.
}
\]
Then
\[
\boxed{
\Theta_1=-C_{\rm tr}\Sigma_{\rm tr}.
}
\]

Likewise, by construction of \(\Sigma_{\rm nt}\),
\[
\boxed{
\Xi_1=A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt}.
}
\]

Finally, the exact selected-branch identity from Stage 165 is
\[
\mathcal R_1
=
-\frac{\epsilon_\eta}{1-\epsilon_\eta}\Sigma_\eta
-\Xi_1.
\]
So
\[
\boxed{
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_\eta}{1-\epsilon_\eta}\Sigma_\eta.
}
\]

This is the exact triangular normal form of the coherent weak-axisymmetric defect.

### Interpretation

The observable hierarchy is now strictly ordered:

1. the tracking-factor drift \(\Theta_1\) measures only \(\Sigma_{\rm tr}\);
2. once the universal tracking feed-through is removed, the grouped defect \(\Xi_1\) measures only \(\Sigma_{\rm nt}\);
3. once the direct transfer-shape defect is added back, the selected-branch residual \(\mathcal R_1+\Xi_1\) measures only \(\Sigma_\eta\).

So the coherent branch has decomposed into a tracking sector, a nontracking transfer-shape sector, and a dressing sector.

---

## 3. Exact inverse reconstruction formulas

Because the normal form is triangular, it can be inverted exactly on the physical branch.

### 3.1 Tracking slippage from \(\Theta_1\)

On the constructive branch \(\chi_0>0\), \(\delta_U>0\), the prefactor \(C_{\rm tr}\) is strictly positive.
So
\[
\boxed{
\Sigma_{\rm tr}
=
-\frac{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}{\chi_0\delta_U}\,\Theta_1.
}
\]

### 3.2 Nontracking slippage from \((\Theta_1,\Xi_1)\)

Using
\[
\Xi_1=A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt},
\]
together with the previous inversion, gives
\[
\Sigma_{\rm nt}
=
\Xi_1+\frac{A_{\rm tr}}{C_{\rm tr}}\Theta_1.
\]
But
\[
\frac{A_{\rm tr}}{C_{\rm tr}}
=
\frac{2(1+\chi_0+\delta_U)}{\delta_U}.
\]
So the exact reconstructed nontracking slippage is
\[
\boxed{
\Sigma_{\rm nt}
=
\Xi_1+\frac{2(1+\chi_0+\delta_U)}{\delta_U}\,\Theta_1.
}
\]

This is a particularly useful formula: it subtracts off the universal tracking feed-through from \(\Xi_1\) and leaves the genuinely nontracking defect.

### 3.3 Dressing slippage from \((\mathcal R_1,\Xi_1)\)

Provided
\[
0<\epsilon_\eta<1,
\]
the dressing prefactor is nonzero, so
\[
\boxed{
\Sigma_\eta
=
-\frac{1-\epsilon_\eta}{\epsilon_\eta}\,(\mathcal R_1+\Xi_1).
}
\]

So the selected-branch dressing mismatch is reconstructed directly from the sum of the selected-branch demand drift and the grouped defect.

---

## 4. Exact zero-defect and rigidity theorems

The triangular normal form makes the rigidity conditions completely explicit.

### 4.1 Exact tracking rigidity

On the physical branch,
\[
\boxed{
\Theta_1=0
\iff
\Sigma_{\rm tr}=0.
}
\]

### 4.2 Exact grouped-defect cancellation

Without assuming tracking rigidity,
\[
\boxed{
\Xi_1=0
\iff
\Sigma_{\rm nt}=-A_{\rm tr}\Sigma_{\rm tr}.
}
\]
So a vanishing grouped defect can still hide a compensating tracking slippage unless \(\Theta_1\) is checked separately.

On the tracking-rigid branch, however,
\[
\boxed{
\Theta_1=0\ \text{and}\ \Xi_1=0
\iff
\Sigma_{\rm tr}=0,\ \Sigma_{\rm nt}=0.
}
\]

### 4.3 Exact selected-branch rigidity

Since
\[
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_\eta}{1-\epsilon_\eta}\Sigma_\eta,
\]
we have
\[
\boxed{
\mathcal R_1+\Xi_1=0
\iff
\Sigma_\eta=0
}
\]
on the physical branch \(0<\epsilon_\eta<1\).

Equivalently, if the direct transfer-shape defect already vanishes,
\[
\Xi_1=0,
\]
then
\[
\boxed{
\mathcal R_1=0
\iff
\Sigma_\eta=0.
}
\]

### 4.4 Full triple-rigidity theorem

Combining all three gives the cleanest statement of the stage:
\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\Sigma_{\rm tr}=\Sigma_{\rm nt}=\Sigma_\eta=0
}
\]
on the constructive coherent branch
\[
\chi_0>0,\qquad \delta_U>0,\qquad 0<\epsilon_\eta<1.
\]

So the full weak-axisymmetric normalization problem has collapsed to the vanishing of exactly three branch-adapted microscopic scalars.

---

## 5. Best current theorem statement after Stage 166

On the actual coherent local D/N tracking branch:

1. the five Stage-165 microscopic slippages compress exactly to the three branch-adapted defect coordinates
   \[
   (\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta);
   \]
2. the observable drift ledger
   \[
   (\Theta_1,\Xi_1,\mathcal R_1)
   \]
   has an exact triangular normal form in these three variables;
3. the inverse reconstruction is explicit:
   \[
   \Sigma_{\rm tr}
   =
   -\frac{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}{\chi_0\delta_U}\Theta_1,
   \]
   \[
   \Sigma_{\rm nt}
   =
   \Xi_1+\frac{2(1+\chi_0+\delta_U)}{\delta_U}\Theta_1,
   \]
   \[
   \Sigma_\eta
   =
   -\frac{1-\epsilon_\eta}{\epsilon_\eta}(\mathcal R_1+\Xi_1);
   \]
4. and exact simultaneous rigidity of tracking, grouped defect, and selected-branch demand is equivalent to
   \[
   \Sigma_{\rm tr}=\Sigma_{\rm nt}=\Sigma_\eta=0.
   \]

So the continuation point is now smaller than Stage 165 left it.

It is no longer

> compute the five microscopic coherent-kernel slippages
> \((\Sigma_Z,\Sigma_\chi,\Sigma_\epsilon,\Sigma_\delta,\Sigma_\eta)\).

It is now

> compute only the three branch-adapted defect coordinates
> \((\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta)\)
> on the actual moving-throat branch.

That is the direct next theorem gate.
