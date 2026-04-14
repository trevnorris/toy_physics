# Same-Charge Barrier Audit — Stage 024: Exact Primitive Ranking on the Selected Twin-Support Branch

## Provenance

This standalone working note is extracted from the corresponding step in `g2_full_output.md`
so the stage-numbered same-charge barrier audit artifacts match the actual calculations that were done.



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