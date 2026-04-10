# Step 28 — Explicit outgoing `l=2` DtN fix of `\chi_Q` and the branch-selection law

## Goal

Step 27 reduced the actual isotropic passive/outgoing grouped-`P2` branch to one
scalar conservative normalization defect,

```math
N_Q := \bar K_0/\bar K_0^{\rm target}.
```

The next honest question is the retarded one:

> **does the actual passive/outgoing branch carry the canonical compact outgoing
> `l=2` odd normalization, or is there still one live outgoing normalization factor?**

This step answers that in two parts.

1. The explicit compact outgoing `l=2` DtN model fixes the last reduced outgoing
   scalar exactly:
   ```math
   \boxed{\chi_Q = 1.}
   ```
2. If the true moving-throat branch is a deformed isotropic DtN branch, the only
   remaining first-order obstruction is the isotropic branch-selection triple
   `(b,a_0,a_5)` through
   ```math
   \boxed{\chi_Q = 1 + 5b + a_0/3 + 9a_5 + O(2).}
   ```

So after this step the remaining reduced uncertainty is no longer vague “outgoing
structure.” It is one explicit DtN branch-selection law.

---

## Step 28A — Exact compact outgoing `l=2` DtN fingerprint

Let

```math
z := \frac{a\omega}{c_s},
```

and take the outgoing partial wave

```math
h_2^{(1)}(z)=j_2(z)+i y_2(z).
```

The exact outgoing `l=2` Dirichlet-to-Neumann operator is

```math
\Lambda_2^{\rm out}(z)
=
z\,\frac{d}{dz}\ln h_2^{(1)}(z)
=
z\,\frac{h_2^{(1)\prime}(z)}{h_2^{(1)}(z)}.
```

Its small-`z` expansion is

```math
\boxed{
\Lambda_2^{\rm out}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}+i\frac{z^5}{9}
-\frac{2z^6}{27}-i\frac{z^7}{27}+O(z^8).
}
```

Normalizing by the static slot gives

```math
\widehat Y_2^{\rm out}(z)
:=
-\frac{3}{\Lambda_2^{\rm out}(z)}
=
1+\frac{z^2}{9}+\frac{4z^4}{81}+i\frac{z^5}{27}
-\frac{11z^6}{729}-i\frac{z^7}{243}+O(z^8).
```

So the canonical compact outgoing odd coefficient is fixed exactly.

---

## Step 28B — Matching to the retarded grouped-`P2` one-pole module

Write the retarded normalized grouped-`P2` one-pole branch as

```math
\widehat Y_Q^{\rm ret}(\omega)
=
\frac34+
\frac14\frac{1}{1-\omega^2/\Omega_Q^2-i\chi_Q\sigma_Q^{\rm can}\omega^5}
+O(\omega^6),
```

with

```math
\sigma_Q^{\rm can}:=\frac{9}{8\Omega_Q^5}.
```

Using the already-carried geometric pole

```math
\Omega_Q = \frac{3c_s}{2a},
```

this becomes

```math
\boxed{\sigma_Q^{\rm can} = \frac{4a^5}{27c_s^5}.}
```

Expanding through `O(ω^5)` gives

```math
\widehat Y_Q^{\rm ret}(\omega)
=
1+rac{a^2\omega^2}{9c_s^2}
+rac{4a^4\omega^4}{81c_s^4}
+i\chi_Q\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6).
```

Comparing with the explicit DtN fingerprint from Step 28A,

```math
\widehat Y_2^{\rm out}(\omega)
=
1+rac{a^2\omega^2}{9c_s^2}
+rac{4a^4\omega^4}{81c_s^4}
+i\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6),
```

forces

```math
\boxed{\chi_Q = 1.}
```

So on the canonical compact outgoing branch the last reduced outgoing normalization
scalar is fixed exactly.

---

## Step 28C — Exact factorization of the last reduced 2.5PN defect

The reduced normalization stack already had the exact factorized condition

```math
\hat m_0^{\,2}\,\chi_Q\,N_Q=1.
```

On the natural point-particle source-map branch,

```math
\hat m_0 = 1 + O(a^2/r^2),
```

so in the strict point-particle limit,

```math
N_Q = 1/\chi_Q.
```

Therefore on the canonical compact outgoing DtN branch,

```math
\boxed{N_Q = 1.}
```

So the reduced nonspinning point-particle 2.5PN theorem is closed on that
canonical outgoing branch.

---

## Step 28D — Exact isotropic DtN deformation law

To isolate the remaining PDE-facing freedom, deform the outgoing DtN operator by

```math
\Lambda_2^{\rm def}(z)
=
S\,\Lambda_2^{\rm out}(\beta z)
+\Sigma_0+\Sigma_2 z^2+\Sigma_4 z^4+i\Sigma_5 z^5+O(z^6).
```

Then the low-frequency coefficients are

```math
L_0 = -3S + \Sigma_0,
```

```math
L_2 = S\beta^2/3 + \Sigma_2,
```

```math
L_4 = S\beta^4/9 + \Sigma_4,
```

```math
L_5 = S\beta^5/9 + \Sigma_5.
```

Normalizing by the actual static slot,

```math
\widehat Y_2^{\rm def}(z)=\frac{L_0}{L_0+L_2 z^2+L_4 z^4+iL_5 z^5}+O(z^6),
```

and demanding that the even fingerprint stay canonical,

```math
\frac{z^2}{9},\qquad \frac{4z^4}{81},
```

fixes

```math
\Sigma_2=-\frac{3S\beta^2-3S+\Sigma_0}{9},
```

```math
\Sigma_4=-\frac{3S\beta^4-3S+\Sigma_0}{27}.
```

With those imposed, the outgoing normalization becomes

```math
\boxed{
\chi_Q=
\frac{3\bigl(S\beta^5+9\Sigma_5\bigr)}{3S-\Sigma_0}.
}
```

So the only isotropic branch data that can move the canonical outgoing value are:

- an argument deformation `β`,
- a static additive slot `Σ_0`,
- and an odd core outlet `Σ_5`.

---

## Step 28E — Linearized branch-selection law

Write a small deformation of the canonical branch as

```math
S = 1 + \varepsilon s,
\qquad
\beta = 1 + \varepsilon b,
\qquad
\Sigma_0 = \varepsilon a_0,
\qquad
\Sigma_5 = \varepsilon a_5.
```

Then the exact deformation law expands to

```math
\boxed{
\chi_Q
=
1+
\varepsilon\left(5b+\frac{a_0}{3}+9a_5\right)
+O(\varepsilon^2).
}
```

So at first order:

- pure overall mouth rescaling drops out,
- effective argument deformation contributes as `5b`,
- the static additive slot contributes as `a_0/3`,
- the odd core outlet contributes as `9a_5`.

Equivalently, the linearized preservation condition is

```math
\boxed{5b+\frac{a_0}{3}+9a_5=0.}
```

---

## Main result of the step

The explicit compact outgoing `l=2` DtN branch fixes the last reduced outgoing
scalar exactly:

```math
\boxed{\chi_Q=1.}
```

So on the natural point-particle source-map branch,

```math
\boxed{N_Q=1.}
```

Any remaining deviation from the canonical outgoing branch is an actual isotropic
DtN branch-selection effect, and at first order it can only enter through the
triple

```math
\boxed{(b,a_0,a_5),}
```

via

```math
\boxed{
\chi_Q
=
1+
5b+\frac{a_0}{3}+9a_5
+O(2).
}
```

That is the sharpest retarded closure statement reached so far.
