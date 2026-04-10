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
