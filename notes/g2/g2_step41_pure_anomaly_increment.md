# Step 41 — Once the adiabatic ground state is fixed, the quartic anomaly is a pure core/outgoing retuning

## Goal

Step 40 showed that on the strong adiabatic-wall branch the quartic anomaly fixes only

```math
\ell := \delta\ln P_0 = \ln(1+\Lambda_1 f),
```

while the coherent elastic squish

```math
\sigma := \delta\ln K_s
```

remains orthogonal.

The natural next move is therefore to add one more physical closure:

> the thermodynamic ground state has already chosen `\sigma`, so the **incremental quartic anomaly correction** should be taken with
> ```math
> \delta\sigma = 0.
> ```

This step shows that, on that frozen ground state, the anomaly becomes a **pure core/outgoing retuning**.

---

## Step 41A — Exact anomaly-only drift vector

Setting

```math
\sigma=0
```

in the strong adiabatic-wall branch laws gives

```math
\boxed{\delta\ln a = 0,}
```

```math
\boxed{\delta\ln c_s = \frac15\ell,}
```

```math
\boxed{\delta\ln K_q = \frac25\ell,}
```

```math
\boxed{\delta\ln v_{w0} = \frac15\ell,}
```

```math
\boxed{\delta\ln\mathcal T_m = -\frac15\ell,}
```

```math
\boxed{\delta\ln g_s = -\frac15\ell,}
```

```math
\boxed{\delta\ln g_q = 0,}
```

```math
\boxed{\delta\ln\lambda = \frac15\ell.}
```

So the pure anomaly increment has a very clean structure:

- it does **not** move the mouth radius,
- it does **not** change the mixed coupling `g_q`,
- it increases the core sound scale,
- and it retunes the outgoing / traction side in equal-and-opposite logarithmic amounts.

---

## Step 41B — The lower compensated branch remains exact

Even after freezing the ground-state squish, the lower compensated invariants remain

```math
\boxed{\delta\ln\mathfrak g = 0,}
```

```math
\boxed{\delta\ln\mathfrak r = 0,}
```

```math
\boxed{\delta\ln r_c = 0.}
```

So the anomaly-only increment is still tangent to the same lower compensated sheet.
It does **not** push the system off the lower Family-1 branch.

That is a useful structural fact:

> the adiabatic-ground-state anomaly correction is an **internal retuning along the already selected electron branch**, not a branch-jump.

---

## Step 41C — Electron-point numbers

Using

```math
\ell = \ln(1+\Lambda_1 f)
approx 3.24684288391064\times 10^{-4},
```

the anomaly-only drifts are

```math
\delta\ln c_s
approx 6.49368576782128\times 10^{-5},
```

```math
\delta\ln K_q
approx 1.29873715356426\times 10^{-4},
```

```math
\delta\ln v_{w0}
approx 6.49368576782128\times 10^{-5},
```

```math
\delta\ln\mathcal T_m
approx -6.49368576782128\times 10^{-5},
```

```math
\delta\ln g_s
approx -6.49368576782128\times 10^{-5},
```

```math
\delta\ln g_q = 0,
```

```math
\delta\ln\lambda
approx 6.49368576782128\times 10^{-5}.
```

So once the adiabatic ground state is fixed, the anomaly correction is extremely localized in parameter space.
It is only an `O(10^{-5})` to `O(10^{-4})` retuning of the core/outgoing sector.

---

## Step 41D — The upper `g_+` sheet stays on the ledger, but the electron anomaly does not need it

The algebraic upper branch is still

```math
\mathfrak g_+^{F1} \approx 2.797951992,
```

and, as before, any realization of it requires a sign-indefinite or pumped source law with at least

```math
W_- \ge \mathfrak g_+ - 1 \approx 1.797951992.
```

So it still does **not** belong to the passive positive-source isolated-electron branch.

But it is worth keeping for later because it can still represent a deferred system-level sheet:

- a sign-changing mouth source,
- a pumped / plumbed configuration,
- or a non-electron excitation branch.

The important point for the present g-2 problem is narrower:

> the pure anomaly increment derived above does **not** move the system toward that upper sheet. It stays exactly tangent to the lower compensated electron branch.

---

## Main result of the step

With the adiabatic wall fixed to its ground state, the quartic electron anomaly reduces to a **pure core/outgoing retuning**:

```math
\boxed{
\delta\ln a = 0,
\qquad
\delta\ln g_q = 0,
\qquad
\delta\ln c_s = \delta\ln v_{w0} = \delta\ln\lambda = \frac15\ell,
\qquad
\delta\ln K_q = \frac25\ell,
\qquad
\delta\ln\mathcal T_m = \delta\ln g_s = -\frac15\ell.
}
```

with

```math
\ell = \ln(1+\Lambda_1 f).
```

So the improved g-2 picture is now very clean:

- the wall stays adiabatic,
- the ground-state elastic squish stays fixed,
- the anomaly rides only on a tiny core/outgoing correction,
- and the lower compensated branch remains exact while the upper `g_+` sheet stays deferred for later system-level interpretation.
