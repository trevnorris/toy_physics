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
