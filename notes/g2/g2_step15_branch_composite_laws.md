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
