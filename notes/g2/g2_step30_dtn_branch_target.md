# Step 30 — The quartic g-2 sliver is an explicit isotropic DtN branch-selection target

## Goal

Step 29 rewrote the remaining quartic anomaly sliver as a tiny outgoing-normalization
drift,

```math
\Delta_Q(f)= -\Lambda_1 f + O(f^2),
```

at tangent level, with

```math
\Lambda_1 \approx 0.279605891931464.
```

Step 28 had already shown that the most general first-order isotropic DtN branch
selection is

```math
\chi_Q
=
1+\varepsilon\left(5b+\frac{a_0}{3}+9a_5\right)+O(\varepsilon^2).
```

The next honest move is therefore obvious:

> **identify the anomaly small parameter with the DtN branch-selection parameter and solve for the required isotropic deformation combination.**

That is what this step does.

The main result is the exact tangent constraint

```math
\boxed{
5b+\frac{a_0}{3}+9a_5 = -\Lambda_1
\approx -0.279605891931464.
}
```

So the quartic electron-anomaly closure is no longer a vague “next common PDE
layer.” It is one explicit isotropic DtN branch-selection target.

---

## Step 30A — Identifying the small parameter

From Step 28,

```math
\chi_Q = 1+\varepsilon\left(5b+\frac{a_0}{3}+9a_5\right)+O(\varepsilon^2).
```

For the anomaly problem the only small reduced parameter already in play is

```math
f=\frac{\alpha_{\rm fs}}{2\pi}.
```

So the natural identification is

```math
\varepsilon=f.
```

Then

```math
\Delta_Q(f):=\chi_Q-1
=
f\left(5b+\frac{a_0}{3}+9a_5\right)+O(f^2).
```

But Step 29 fixed the same outgoing-defect tangent to be

```math
\Delta_Q(f)=-\Lambda_1 f+O(f^2).
```

Matching the two immediately gives

```math
\boxed{
5b+\frac{a_0}{3}+9a_5 = -\Lambda_1.
}
```

This is the exact isotropic DtN branch-selection law required by the quartic
electron-anomaly sliver.

---

## Step 30B — Direct quartic anomaly law in DtN variables

Since Step 29 also showed

```math
\Delta\!\left(\frac g2\right)
=
-c_{3,\rm total}\,\delta_1\,f^4+O(f^5),
\qquad
\delta_1=-\Lambda_1,
```

the quartic anomaly correction can now be written directly as

```math
\boxed{
\Delta\!\left(\frac g2\right)
=
-c_{3,\rm total}
\left(5b+\frac{a_0}{3}+9a_5\right)f^4
+
O(f^5).
}
```

So the missing `O(f^4)` layer is exactly the first isotropic DtN deformation of
the outgoing grouped-`P2` branch.

---

## Step 30C — Three pure isotropic realizations

The constraint

```math
5b+\frac{a_0}{3}+9a_5 = -\Lambda_1
```

does not yet tell us which microscopic branch the PDE picks.
But it already gives three clean pure realizations.

### Pure argument deformation

Set `a_0=a_5=0`. Then

```math
\boxed{
b = -\frac{\Lambda_1}{5}
\approx -0.0559211783862928.
}
```

### Pure static additive slot

Set `b=a_5=0`. Then

```math
\boxed{
a_0 = -3\Lambda_1
\approx -0.838817675794392.
}
```

### Pure odd core outlet

Set `b=a_0=0`. Then

```math
\boxed{
a_5 = -\frac{\Lambda_1}{9}
\approx -0.0310673213257182.
}
```

These are not all equally plausible physically, but they are exact first-order
bookkeeping realizations of the same anomaly target.

---

## Step 30D — Minimum-norm bookkeeping branch

If one ignores the physical dimensions of `(b,a_0,a_5)` and asks only for the
smallest Euclidean-norm tangent triple satisfying the constraint, the answer is

```math
\boxed{
(b,a_0,a_5)_{\min}
=
-\Lambda_1
\frac{(5,\;1/3,\;9)}{5^2+(1/3)^2+9^2}.
}
```

Numerically,

```math
\boxed{
(b,a_0,a_5)_{\min}
\approx
(-0.0131751467402261,\;
 -0.000878343116015070,\;
 -0.0237152641324069).
}
```

This is best read as a bookkeeping baseline, not as a real physical preference,
because the three deformation coordinates carry different meanings in the DtN
problem.

---

## Step 30E — Actual electron-point size of the deformation

The tangent combination itself is `O(1)`, but the actual electron-point
deformation is smaller by one factor of `f`:

```math
f\left(5b+\frac{a_0}{3}+9a_5\right)
=
-\Lambda_1 f
\approx -3.24737004039746\times 10^{-4}.
```

That is the linearized electron-point outgoing defect.
The exact Step-29 branch value was

```math
\Delta_Q^{(e)}
=
-\frac{\Lambda_1 f}{1+\Lambda_1 f}
\approx
-3.24631584151692\times 10^{-4},
```

so the linearized DtN law is already accurate to the expected next neglected order.

---

## Main result of the step

The quartic g-2 sliver is now pinned to one explicit isotropic DtN branch-selection
constraint:

```math
\boxed{
5b+\frac{a_0}{3}+9a_5 = -\Lambda_1
\approx -0.279605891931464.
}
```

Equivalently,

```math
\boxed{
\Delta\!\left(\frac g2\right)
=
-c_{3,\rm total}
\left(5b+\frac{a_0}{3}+9a_5\right)f^4
+O(f^5).
}
```

So the next PDE-facing question is no longer broad at all:

> **which actual moving-throat isotropic DtN branch produces the required combination `5b + a_0/3 + 9a_5`?**

That is the cleanest continuation point after the outgoing-defect bridge.
