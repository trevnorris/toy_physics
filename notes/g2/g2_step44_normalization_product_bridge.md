# Step 44 — The adiabatic pure-odd DtN deformation rescales the grouped-`P2` normalization product by `e^\ell`

## Goal

Step 43 compressed the quartic anomaly into a very clean normalized DtN statement:

```math
\chi_Q = e^{-\ell},
\qquad
\beta = 1,
\qquad
\Sigma_0 = 0,
\qquad
\Sigma_5 = \frac{S}{9}(e^{-\ell}-1).
```

The next honest move is to connect that result back to the moving-throat grouped-`P2`
normalization product. The right bridge already existed in reduced form:

```math
\hat m_0^{\,2}\,\chi_Q\,N_Q = 1.
```

This step shows that on the natural point-particle source-map branch the pure-odd
DtN deformation forces the grouped normalization defect to move in the **opposite**
direction,

```math
\boxed{N_Q=e^{+\ell},}
```

and therefore the moving-throat grouped-`P2` prefactor obeys

```math
\boxed{
P_0 = e^\ell\,P_0^{\rm target},
\qquad
P_0^{\rm target}=\frac{54Gc_s^5}{5a^5c^5}.
}
```

So the pure-odd outlet and the grouped-`P2` normalization product are reciprocal
descriptions of the same quartic branch motion.

---

## Step 44A — Reciprocal normalization law on the natural source-map branch

From the reduced factorization,

```math
\hat m_0^{\,2}\,\chi_Q\,N_Q = 1,
```

together with the Step-43 pure-odd law

```math
\chi_Q=e^{-\ell},
```

we get

```math
N_Q = \frac{e^\ell}{\hat m_0^{\,2}}.
```

On the natural point-particle source-map branch used throughout the reduced
electron closure,

```math
\hat m_0 = 1,
```

so this collapses to

```math
\boxed{N_Q=e^\ell.}
```

That already tells us something important: the same quartic anomaly that *reduces*
the normalized outlet factor `\chi_Q` *increases* the grouped normalization defect
`N_Q` by the inverse amount.

---

## Step 44B — The grouped-`P2` normalization product

From the grouped-`P2` normalization side, Step 27 had already shown that the single
reduced defect rescales the odd grouped coefficient:

```math
\bar\Gamma_5 = N_Q\,\bar\Gamma_5^{\rm target}
= N_Q\,\frac{2G}{5c^5}.
```

From the moving-throat grouped-`P2` bridge, the same coefficient is written as

```math
\bar\Gamma_5 = \frac{P_0 a^5}{27c_s^5}.
```

So the static prefactor is

```math
P_0 = N_Q\,\frac{54Gc_s^5}{5a^5c^5}.
```

On the natural source-map branch this becomes

```math
\boxed{
P_0=e^\ell\,P_0^{\rm target},
\qquad
P_0^{\rm target}:=\frac{54Gc_s^5}{5a^5c^5}.
}
```

So the Step-43 pure-odd branch selection is *exactly* the same branch motion as a
static grouped-`P2` prefactor rescaling by `e^\ell`.

---

## Step 44C — Exact reciprocal law

Multiplying the two descriptions together gives

```math
\boxed{
P_0\,\chi_Q = P_0^{\rm target}.
}
```

That is the cleanest statement of the bridge:

- in normalized DtN language the anomaly is a **pure odd outlet deformation**
  `\chi_Q=e^{-\ell}`,
- in grouped-`P2` normalization language the same anomaly is a **static prefactor
  enhancement** `P_0=e^\ell P_0^{\rm target}`.

Those are exact reciprocals on the natural source-map branch.

---

## Step 44D — Electron-point size

Using the carried electron value

```math
\ell = \ln(1+\Lambda_1 f),
```

we have

```math
e^\ell = 1+\Lambda_1 f.
```

Numerically,

```math
\ell \approx 3.24684288391064\times 10^{-4},
```

```math
\boxed{
e^\ell \approx 1.00032473700404.
}
```

So the grouped-`P2` normalization product is enhanced by about

```math
\boxed{324.737\ {\rm ppm}.}
```

That is the exact finite-`f` bridge from the quartic anomaly sliver to the
moving-throat normalization product.

---

## Main result of the step

The pure-odd adiabatic DtN deformation is not just an outlet-side statement.
On the natural point-particle source-map branch it forces the grouped-`P2`
normalization product to rescale inversely:

```math
\boxed{
\chi_Q=e^{-\ell}
\quad\Longleftrightarrow\quad
N_Q=e^\ell
\quad\Longleftrightarrow\quad
P_0=e^\ell\,P_0^{\rm target}.
}
```

with

```math
\boxed{
P_0^{\rm target}=\frac{54Gc_s^5}{5a^5c^5}.
}
```

So the quartic anomaly is now directly attached to the exact moving-throat
normalization product rather than only to an abstract DtN branch-selection scalar.
