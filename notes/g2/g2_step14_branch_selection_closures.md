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
