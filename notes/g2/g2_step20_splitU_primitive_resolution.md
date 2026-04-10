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
