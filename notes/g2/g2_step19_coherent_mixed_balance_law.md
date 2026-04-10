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
