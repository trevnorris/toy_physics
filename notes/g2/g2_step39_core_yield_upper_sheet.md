# Step 39 — Minimal core-yield closure and the retained upper-sheet diagnostic

## Goal

Step 38 used your adiabatic-wall ground-state choice to reduce the electron-track branch to the three reduced observables

```math
(K_s,
K_q,
P_0),
```

while preserving the exact lower compensated parent surface.

The next natural question is then:

> what is the smallest additional closure that treats the wall as frozen/adiabatic while still letting the **core** yield?

The cleanest minimal answer is to freeze the transverse mixed-channel localization norm,

```math
\delta\ln\mathcal Z_q=0,
```

and let the remaining nontrivial response be carried only by

- coherent elastic wall squish, and
- outgoing/core compressibility.

This step shows that under that choice the whole reduced electron track becomes **rank 2**.

It also keeps your side note about the algebraic upper branch alive: `\mathfrak g_+` is not deleted, but it is reinterpreted as a sign-indefinite / pumped sheet rather than the passive positive-source electron branch.

---

## Step 39A — Minimal core-yield closure

Impose, in addition to Step 38,

```math
\boxed{\delta\ln\mathcal Z_q=0.}
```

This is not an exact theorem of the current file stack. It is a **minimal closure choice** consistent with your new physical picture:

- wall thermodynamics frozen,
- wall density frozen,
- coherent elastic squish allowed,
- core compressibility allowed,
- no extra side-channel relocalization introduced unless it is forced later.

Then the exact lower-branch law

```math
\delta\ln\mathcal Z_q
=
\delta\ln K_q-
\frac25\delta\ln P_0
```

immediately gives

```math
\boxed{
\delta\ln K_q
=
\frac25\,\delta\ln P_0.
}
```

So the mixed side-channel stiffness is no longer independent either.

---

## Step 39B — Exact two-amplitude reduced electron track

Using Step 38 plus `\delta\ln\mathcal Z_q=0`, the entire reduced branch becomes

```math
\boxed{
\delta\ln a = \frac12\,\delta\ln K_s,
}
```

```math
\boxed{
\delta\ln c_s
=
\frac12\,\delta\ln K_s
+\frac15\,\delta\ln P_0,
}
```

```math
\boxed{
\delta\ln K_q
=
\frac25\,\delta\ln P_0.
}
```

The transported outlet/background variables become

```math
\boxed{
\delta\ln v_{w0}
=
-\frac34\,\delta\ln K_s
+\frac15\,\delta\ln P_0,
}
```

```math
\boxed{
\delta\ln\mathcal T_m
=
-\frac54\,\delta\ln K_s
-\frac15\,\delta\ln P_0,
}
```

and the sum/difference channels are

```math
\boxed{
\delta\ln\!\left(\frac{v_{w0}}{\mathcal T_m}\right)
=
\frac12\,\delta\ln K_s
+\frac25\,\delta\ln P_0,
}
```

```math
\boxed{
\delta\ln(v_{w0}\mathcal T_m)
=
-2\,\delta\ln K_s.
}
```

So the product channel is now completely blind to `P_0`; only the elastic squish moves it.

That is a very sharp structural statement.

---

## Step 39C — Parent couplings on the two-amplitude branch

The exact parent-action couplings reduce to

```math
\boxed{
\delta\ln g_s
=
-\frac14\,\delta\ln K_s
-\frac15\,\delta\ln P_0,
}
```

```math
\boxed{
\delta\ln g_q
=
-\frac34\,\delta\ln K_s,
}
```

```math
\boxed{
\delta\ln\lambda
=
\frac12\,\delta\ln K_s
+\frac15\,\delta\ln P_0.
}
```

A useful consequence is that

```math
\delta\ln g_q
```

is completely blind to the outgoing-normalization drift `P_0` on this minimal branch. The mixed mouth coupling tracks only the elastic shell squish.

So the two surviving amplitudes now have a clean interpretation:

- `\delta\ln K_s` = **elastic squish amplitude**,
- `\delta\ln P_0` = **outgoing/core-yield amplitude**.

---

## Step 39D — Exact rank-2 statement

The reduced drift vector

```math
(
\delta\ln K_q,
\delta\ln a,
\delta\ln c_s,
\delta\ln v_{w0},
\delta\ln\mathcal T_m,
\delta\ln g_s,
\delta\ln g_q,
\delta\ln\lambda
)
```

depends only on the two amplitudes

```math
(\delta\ln K_s,
\delta\ln P_0).
```

The SymPy Jacobian of that map has rank exactly `2`.

So the adiabatic electron-track branch is now reduced to the smallest plausible parameter space we have reached so far:

```math
\boxed{
\text{electron track} 
\sim
(\text{elastic squish},\ \text{outgoing/core-yield}).
}
```

---

## Step 39E — The lower compensated surface stays exact

Even on this stronger two-amplitude branch, the exact compensation surface remains untouched:

```math
\boxed{
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)=0,
}
```

```math
\boxed{
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right)=0,
}
```

so

```math
\boxed{\delta\ln\mathfrak g = 0,}
```

```math
\boxed{\delta\ln\mathfrak r = 0,}
```

```math
\boxed{\delta\ln r_c = 0.}
```

So the minimal core-yield closure does **not** re-open the off-family scalar drift. The branch stays exactly on the lower compensated sheet.

---

## Step 39F — Keeping `\mathfrak g_+` alive as a system-wide diagnostic sheet

You asked not to throw away the algebraic upper value

```math
\mathfrak g_+^{F1}
\approx 2.79795199200529.
```

I agree that we should keep it.

The positive-source theorem only tells us that for any nonnegative normalized mouth source profile,

```math
0\le \mathfrak g[\sigma]\le 1.
```

So `\mathfrak g_+>1` means only that the upper branch is **not** the passive positive-source electron branch.

It does **not** mean the branch is mathematically meaningless.

To make that precise, write a sign-indefinite normalized source law as

```math
\sigma = \sigma_+ - \sigma_-,
\qquad
W_+ = \int \sigma_+,
\qquad
W_- = \int \sigma_-,
\qquad
W_+ - W_- = 1.
```

Since the cosine weight on the first D/N interval lies in `[0,1]`, one has

```math
\mathfrak g[\sigma]
\le W_+
=1+W_-.
```

Therefore any realization of `\mathfrak g>1` requires

```math
\boxed{W_- \ge \mathfrak g-1.}
```

Applied to the explicit Family-1 upper sheet,

```math
\boxed{
W_-^{\min}
=
\mathfrak g_+^{F1}-1
\approx 1.79795199200529,
}
```

with corresponding positive weight

```math
\boxed{
W_+^{\min}
=
\mathfrak g_+^{F1}
\approx 2.79795199200529.
}
```

So the upper sheet can be retained as a very specific diagnostic branch:

- not passive,
- not positive-source,
- requiring strong sign-indefinite / source-sink / pumped support,
- and therefore potentially telling us about a broader system-level sector rather than the isolated electron ground state.

That seems exactly in line with your side note.

---

## Main result of the step

The stronger adiabatic electron-track closure

```math
\delta\ln\mathcal Z_q=0
```

reduces the whole lower-branch problem to only two amplitudes:

```math
\boxed{
\delta\ln K_s
\quad\text{and}\quad
\delta\ln P_0.
}
```

So the g-2/electron-track search is no longer spread across many microscopic directions.
It has collapsed to

1. a coherent elastic wall-squish lane,
2. and a core/outgoing normalization lane.

At the same time, the algebraic upper sheet `\mathfrak g_+` is preserved as a flagged non-passive branch rather than deleted.

That is a good place to stand:

- the physical electron branch is cleaner,
- the non-electron algebraic sheet is still remembered,
- and the next step can finally ask how the quartic anomaly target picks out a particular combination of the two surviving amplitudes.
