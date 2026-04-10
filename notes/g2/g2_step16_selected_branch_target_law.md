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
