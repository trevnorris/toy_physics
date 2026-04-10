# Step 51 — The canonical compact outgoing `l=2` DtN branch is a genuine no-tuning closure

## Goal

Step 50 collapsed the grouped-`P2` side of the g-2 problem to one scalar defect

```math
N_Q.
```

The next question is whether the outgoing branch itself is still an adjustable datum, or whether the canonical compact outgoing `l=2` branch already fixes it with no tuning.

This step shows the stronger result:

1. the explicit outgoing spherical `l=2` DtN model fixes the canonical normalized branch,
2. the corresponding grouped-`P2` normalization is
   ```math
   \chi_Q=1,
   ```
3. and on the natural point-particle source-map branch that implies
   ```math
   N_Q=1.
   ```

So the canonical passive/outgoing branch is a real no-tuning closure.

---

## Step 51A — Exact outgoing `l=2` DtN fingerprint

Using the exact spherical Hankel outgoing mode, the `l=2` DtN operator is

```math
\Lambda_2^{\rm out}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}+i\frac{z^5}{9}+O(z^6),
\qquad z:=\frac{a\omega}{c_s}.
```

Normalizing by the static slot gives

```math
\widehat Y_2^{\rm out}(z)
=
1+\frac{z^2}{9}+\frac{4z^4}{81}+i\frac{z^5}{27}+O(z^6).
```

So the leading outgoing odd coefficient is fixed exactly:

```math
\boxed{\Gamma_{5,\rm can}^{\rm DtN} = \frac{a^5}{27c_s^5}.}
```

This is not a free parameter once the canonical outgoing branch is chosen.

---

## Step 51B — Matching to the minimal retarded grouped-`P2` module

The minimal isotropic retarded grouped-`P2` module is

```math
\widehat Y_Q^{\rm ret}(\omega)
=
\frac34+
\frac14\frac{1}{1-\omega^2/\Omega_Q^2-i\chi_Q\sigma_Q^{\rm can}\omega^5}
+O(\omega^6),
```

with

```math
\sigma_Q^{\rm can}=\frac{9}{8\Omega_Q^5}.
```

Using the already-carried pole relation

```math
\Omega_Q=\frac{3c_s}{2a},
```

this becomes

```math
\sigma_Q^{\rm can}=\frac{4a^5}{27c_s^5}.
```

Expanding the retarded module then gives

```math
\widehat Y_Q^{\rm ret}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\chi_Q\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6).
```

Matching to the exact outgoing DtN fingerprint yields immediately

```math
\boxed{\chi_Q=1.}
```

So the canonical compact outgoing branch fixes the normalized outgoing quadrupole factor exactly.

---

## Step 51C — Natural source-map no-tuning closure

The reduced normalization stack is

```math
\hat m_0^{\,2}\chi_Q N_Q = 1.
```

On the natural point-particle source-map branch,

```math
\hat m_0=1.
```

So if the actual outgoing branch is the canonical compact one,

```math
\chi_Q=1,
```

then automatically

```math
\boxed{N_Q=1.}
```

That means the canonical passive/outgoing grouped-`P2` branch closes the normalization stack with **no further fitting at all**.

---

## Step 51D — Robustness classes

The later moving-throat DtN analysis also shows that two simple deformation classes do **not** move the canonical value once the even moments stay canonical.

### Pure scale deformation

```math
\Lambda_2^{\rm def}(z)=S\Lambda_2^{\rm out}(z)
\quad\Longrightarrow\quad
\chi_Q=1.
```

### Pure scale+argument deformation

```math
\Lambda_2^{\rm def}(z)=S\Lambda_2^{\rm out}(\beta z)
```

gives

```math
\widehat Y_2^{\rm def}(z)
=
1+\frac{\beta^2 z^2}{9}+\frac{4\beta^4 z^4}{81}+i\frac{\beta^5 z^5}{27}+\cdots.
```

If the even fingerprint is kept canonical, then

```math
\beta=1,
```

and therefore again

```math
\boxed{\chi_Q=1.}
```

So simple overall normalization or effective radius/sound-speed rescaling cannot explain a nontrivial outgoing defect once the lower even moments are already fixed.

---

## Main result of the step

The canonical compact outgoing `l=2` DtN branch is a genuine **no-tuning** closure:

```math
\boxed{\chi_Q=1,\qquad N_Q=1.}
```

That matters because it cleanly separates two questions:

- the canonical passive/outgoing branch itself is fixed without tuning,
- any nonzero electron anomaly must therefore come from a **genuine isotropic throat-core DtN deformation** away from that canonical branch.
