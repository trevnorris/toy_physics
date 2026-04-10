# Step 29 — The quartic g-2 sliver is a tiny finite-`f` outgoing-normalization drift

## Goal

Step 28 fixed the canonical compact outgoing grouped-`P2` branch and showed that,
on the natural source-map branch,

```math
N_Q-1 = -\frac{\Delta_Q}{1+\Delta_Q},
\qquad
\Delta_Q := \chi_Q - 1.
```

The next honest question is:

> **can the remaining quartic electron-anomaly sliver be carried by a tiny
> finite-`f` drift of that same outgoing normalization, without spoiling the
> strict point-particle gravitational branch?**

This step shows the answer is yes.

The main result is that the missing quartic layer can be rewritten as an
`O(f)` outgoing-normalization defect:

```math
\boxed{
\Delta_Q(f)=\delta_1 f + O(f^2),
\qquad
\delta_1 = -\Lambda_1,
}
```

where `\Lambda_1` is the Step-2 common-path tangent fixed by the carried benchmark
residual,

```math
\Lambda_1 \approx 0.279605891931464.
```

So the exact tangent defect is

```math
\boxed{\delta_1 \approx -0.279605891931464.}
```

At the physical electron point this becomes a very small actual outgoing defect,

```math
\boxed{
\Delta_Q^{(e)}
=
-\frac{\Lambda_1 f}{1+\Lambda_1 f}
\approx -3.24631584151692\times 10^{-4},
}
```

with

```math
\chi_Q^{(e)}\approx 0.999675368415848,
\qquad
N_Q^{(e)}\approx 1.00032473700404.
```

So the quartic closure does **not** require a large renormalization of the
canonical outgoing branch. It requires only a tiny finite-`f` drift away from it.

---

## Step 29A — Exact bridge from outgoing defect to common conservative defect

On the natural source-map branch, Step 28 gave

```math
N_Q = \frac{1}{\chi_Q} = \frac{1}{1+\Delta_Q}.
```

Therefore

```math
\boxed{
N_Q - 1
=
-\frac{\Delta_Q}{1+\Delta_Q}.
}
```

This is the exact scalar bridge between the retarded outgoing defect `\Delta_Q`
and the conservative normalization defect `N_Q-1`.

---

## Step 29B — Why the outgoing defect must start at `O(f)`

Step 2 already showed that the missing common anomaly layer acts only on the
already-existing local transport residue, not on the whole anomaly law. The
general leading form was

```math
\Delta\!\left(\frac g2\right)
=
c_{3,\rm total}\,\Lambda_{\rm common}(f)\,f^3
+
O(f^5),
```

and to close a quartic gap one needs

```math
\Lambda_{\rm common}(f)=O(f).
```

If we now identify the common scalar with the conservative branch defect,

```math
\Lambda_{\rm common}(f)=N_Q(f)-1,
```

then the outgoing defect must also begin linearly:

```math
\boxed{
\Delta_Q(f)=\delta_1 f + O(f^2).
}
```

This is the key compatibility point.

- In the strict point-particle / weak-coupling limit `f\to0`, the outgoing defect
  vanishes and the canonical compact outgoing branch is recovered.
- At finite electron `f`, the same branch can still carry a tiny deformation that
  shows up only in the quartic anomaly layer.

So the gravitational reduced theorem and the electron-anomaly closure are not in
conflict.

---

## Step 29C — Quartic anomaly law in terms of `\delta_1`

Insert

```math
\Delta_Q(f)=\delta_1 f + O(f^2)
```

into

```math
N_Q(f)-1=-\frac{\Delta_Q(f)}{1+\Delta_Q(f)}.
```

Then

```math
N_Q(f)-1 = -\delta_1 f + O(f^2).
```

So the common anomaly correction becomes

```math
\Delta\!\left(\frac g2\right)
=
-c_{3,\rm total}\,\delta_1\,f^4
+
O(f^5).
```

Matching this to the Step-2 carried benchmark

```math
\Delta\!\left(\frac g2\right)
=
c_{3,\rm total}\,\Lambda_1\,f^4
+
O(f^5)
```

forces

```math
\boxed{
\delta_1 = -\Lambda_1.
}
```

Numerically,

```math
\boxed{
\delta_1 \approx -0.279605891931464.
}
```

So the quartic anomaly sliver is exactly the first tangent of the outgoing
normalization defect.

---

## Step 29D — Actual electron-point defect

The linear tangent coefficient `\delta_1` is not the physical deformation itself.
The actual electron-point defect is of order `f` smaller:

```math
\Delta_Q^{(e)}
=
-\frac{\Lambda_1 f}{1+\Lambda_1 f}.
```

Using the carried benchmark value

```math
f = 0.001161409732093,
```

gives

```math
\boxed{
\Lambda_1 f \approx 3.24737004039746\times 10^{-4},
}
```

and therefore

```math
\boxed{
\Delta_Q^{(e)}
\approx
-3.24631584151692\times 10^{-4}.
}
```

The corresponding branch values are

```math
\boxed{
\chi_Q^{(e)} = 1+\Delta_Q^{(e)} \approx 0.999675368415848,
}
```

```math
\boxed{
N_Q^{(e)} = \frac{1}{1+\Delta_Q^{(e)}} \approx 1.00032473700404.
}
```

So the needed defect is only a few parts in `10^4`.

---

## Step 29E — Exact improved anomaly law

The improved reduced anomaly law can now be written directly in terms of the
outgoing defect:

```math
\boxed{
\frac{g_{\rm imp}(f)}{2}
=
\frac{g_{\rm loc}(f)}{2}
-
c_{3,\rm total}\,
\frac{\Delta_Q(f)}{1+\Delta_Q(f)}\,f^3
+
O(f^5).
}
```

Choosing the electron-point defect above reproduces the carried exact residual:

```math
g_e-g_{\rm loc}\approx 2.27204390584705\times 10^{-12}.
```

So the quartic gap is exactly saturated by this tiny finite-`f` outgoing
normalization drift.

---

## Main result of the step

The remaining quartic g-2 sliver can be encoded as a tiny outgoing-normalization
defect that vanishes in the strict point-particle limit:

```math
\boxed{
\Delta_Q(f)= -\Lambda_1 f + O(f^2)
}
```

at tangent level, or more explicitly at the electron point,

```math
\boxed{
\Delta_Q^{(e)}
=
-\frac{\Lambda_1 f}{1+\Lambda_1 f}
\approx -3.24631584151692\times10^{-4}.
}
```

This is a strong simplification:

- the gravitational reduced theorem still wants the canonical compact outgoing
  branch as `f\to0`,
- while the finite electron anomaly needs only a very small `O(f)` departure from
  that branch.

So the next clean derivation is to translate this outgoing-defect tangent into the
explicit isotropic DtN branch-selection variables from Step 28.
