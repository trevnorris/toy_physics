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
