# Step 40 — On the strong adiabatic-wall branch the quartic anomaly becomes a tiny relative core-stiffening law

## Goal

The visible Step 39 closure already reduced the adiabatic electron branch to two amplitudes:

- `\sigma := \delta\ln K_s` — coherent elastic wall squish,
- `\ell := \delta\ln P_0` — isotropic static load / anomaly-carrying normalization drift.

So the next honest question is now very sharp:

> once the quartic g-2 sliver is inserted, what does it do to the adiabatic-wall branch microscopically?

This step answers that.

---

## Step 40A — Exact strong-closure kinematics

From the strong adiabatic-wall closure,

```math
\delta\ln\Theta_w=0,
\qquad
\delta\ln\mathcal Z_q=0,
```

the reduced branch laws are

```math
\delta\ln a = \frac12\sigma,
```

```math
\delta\ln c_s = \frac12\sigma + \frac15\ell,
```

```math
\delta\ln K_q = \frac25\ell.
```

So the relative core stiffening is

```math
\boxed{
\delta\ln\!\left(\frac{c_s}{a}\right)
=
\frac15\ell
=
\frac12\,\delta\ln K_q.
}
```

This is already strong physically:

- the adiabatic wall direction `\sigma` drops out completely from `c_s/a`,
- the only thing that changes the core stiffness **relative to wall size** is the anomaly-carrying load drift `\ell`.

---

## Step 40B — Exact outgoing-normalization law from the quartic sliver

The earlier quartic bridge already rewrote the anomaly miss as a tiny outgoing-normalization defect,

```math
\Delta_Q(f)= -\frac{\Lambda_1 f}{1+\Lambda_1 f},
```

with

```math
\Lambda_1\approx 0.279605891931464.
```

On the natural source-map branch,

```math
N_Q = \frac{1}{1+\Delta_Q}.
```

Substituting the exact `\Delta_Q(f)` gives

```math
\boxed{N_Q(f)=1+\Lambda_1 f.}
```

So on the strong adiabatic-wall branch,

```math
\boxed{
\ell
=
\delta\ln P_0
=
\ln N_Q
=
\ln(1+\Lambda_1 f).
}
```

At linear order,

```math
\ell \approx \Lambda_1 f,
```

but the exact logarithmic form is the cleaner finite-`f` law once the outgoing defect is already known.

---

## Step 40C — Exact microscopic anomaly laws on the strong adiabatic-wall branch

Substituting

```math
\ell = \ln(1+\Lambda_1 f)
```

into the strong closure gives

```math
\boxed{
\delta\ln K_q
=
\frac25\ln(1+\Lambda_1 f),
}
```

```math
\boxed{
\delta\ln\!\left(\frac{c_s}{a}\right)
=
\frac15\ln(1+\Lambda_1 f).
}
```

So the quartic anomaly is not asking for a broad thermodynamic rearrangement.
On this branch it asks only for a very small **relative core stiffening**.

That is the sharpest physical statement reached so far:

> with an adiabatic wall, the electron quartic sliver is carried by a tiny compressibility shift of the core relative to the wall size, while the coherent wall squish remains orthogonal.

---

## Step 40D — Electron-point numbers

Using the carried electron value

```math
f \approx 0.001161409732093,
```

one gets

```math
\Lambda_1 f \approx 3.24737004039746\times 10^{-4},
```

```math
\ell
=
\ln(1+\Lambda_1 f)
\approx
3.24684288391064\times 10^{-4},
```

```math
\delta\ln K_q
\approx
1.29873715356426\times 10^{-4},
```

```math
\delta\ln(c_s/a)
\approx
6.49368576782128\times 10^{-5}.
```

So the required correction is genuinely tiny.
The exact-vs-linear difference in `\ell` is only about

```math
5.27\times 10^{-8},
```

which means the linearized quartic bookkeeping is already numerically excellent.

---

## What remains free

This is also the first point where the branch splits cleanly into two jobs:

1. the thermodynamic ground state / coherent wall shape chooses
   ```math
   \sigma = \delta\ln K_s,
   ```
2. the quartic anomaly closure chooses
   ```math
   \ell = \delta\ln P_0.
   ```

So the anomaly does **not** determine the whole branch.
It only determines the very small compressibility/load correction.
The coherent elastic squish remains an orthogonal ground-state datum.

That is exactly what one would hope if the wall is truly adiabatic.
