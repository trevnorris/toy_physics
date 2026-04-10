# Step 42 — The adiabatic anomaly increment freezes the parent compensation sheet and detunes only the odd outlet

## Goal

Step 41 showed that once the adiabatic wall is fixed to its thermodynamic ground state, the quartic anomaly is a **pure core/outgoing retuning**:

```math

delta\ln K_s = 0,
\qquad
\ell := \delta\ln P_0,
```

with

```math
\delta\ln a = 0,
\qquad
\delta\ln c_s = \frac15\ell,
\qquad
\delta\ln K_q = \frac25\ell,
```

```math
\delta\ln v_{w0} = \frac15\ell,
\qquad
\delta\ln\mathcal T_m = -\frac15\ell,
\qquad
\delta\ln g_s = -\frac15\ell,
\qquad
\delta\ln g_q = 0,
\qquad
\delta\ln\lambda = \frac15\ell.
```

The next honest question is therefore sharper than before:

> does this anomaly-only increment move the system around on the parent compensation family, or does it stay on the same branch and only retune the outlet?

This step shows the stronger result:

```math
\boxed{
\text{the anomaly-only adiabatic branch leaves }\mathfrak g,\ \mathfrak r,\ r_c,\ L_W\ \text{all fixed.}
}
```

So the quartic sliver does **not** come from motion in the parent compensation ratios. It comes only from a small odd-outlet detuning on top of a slightly softened loading share.

---

## Step 42A — Exact parent-sheet invariants are frozen

The lower compensated family is parameterized by the parent ratios

```math
\mathfrak r := \frac{\lambda}{\sqrt{K_s K_q}},
\qquad
\mathfrak g := \frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}},
\qquad
r_c := \frac{\lambda^2}{K_s K_q} = \mathfrak r^2.
```

Substituting the Step-41 anomaly-only drift vector gives

```math
\boxed{
\delta\ln\mathfrak r
=
\delta\ln\lambda-
\frac12(\delta\ln K_s+\delta\ln K_q)
=0,
}
```

```math
\boxed{
\delta\ln\mathfrak g
=
\delta\ln g_q+
\frac12\delta\ln K_s-
\delta\ln g_s-
\frac12\delta\ln K_q
=0,
}
```

```math
\boxed{
\delta\ln r_c
=
2\,\delta\ln\lambda-
\delta\ln K_s-
\delta\ln K_q
=0.
}
```

So the anomaly-only adiabatic increment stays **exactly tangent to the same lower compensated parent sheet**.

That means the quartic g-2 sliver is **not** telling the electron branch to move to a new `\mathfrak g`/`\mathfrak r` balance point.

---

## Step 42B — The auxiliary D/N geometry stays fixed

The balanced auxiliary side-tube geometry obeys

```math
L_W = \frac{\pi a}{2}\sqrt{\frac{1+r_c}{3}}.
```

Since Step 41 already gave

```math
\delta\ln a = 0,
\qquad
\delta\ln r_c = 0,
```

we get immediately

```math
\boxed{\delta\ln L_W = 0.}
```

So the quartic anomaly does **not** ask for a new even D/N geometry either.
The balanced auxiliary tube length remains frozen.

This is important because it means the anomaly is not reopening the even outlet side of the problem.

---

## Step 42C — The balanced loading share drifts downward

Although the parent ratios stay fixed, the shell-side loading share does not.
On the balanced core,

```math
\rho_c := \frac{g_s^2}{K_s},
\qquad
\sigma_c = \frac{\rho_c}{4}.
```

Using

```math
\delta\ln g_s = -\frac15\ell,
\qquad
\delta\ln K_s = 0,
```

gives

```math
\boxed{
\delta\ln\rho_c = 2\,\delta\ln g_s - \delta\ln K_s = -\frac25\ell,
}
```

```math
\boxed{
\delta\ln\sigma_c = -\frac25\ell.
}
```

So the exact finite-`\ell` laws are

```math
\boxed{
\rho_c(\ell)=\rho_{c,*} e^{-2\ell/5},
\qquad
\sigma_c(\ell)=\sigma_{c,*} e^{-2\ell/5}.
}
```

Meanwhile the even-preserving mixed-side coefficient stays frozen:

```math
\boxed{\kappa_c = \frac13.}
```

So the anomaly-only adiabatic branch keeps the same **parent balance** and the same **even side-channel geometry**, but it reduces the loading share slightly because `g_s` softens while `K_s` is held fixed.

---

## Step 42D — Exact odd detuning law on the fixed parent sheet

The electron anomaly target on the outgoing branch is

```math
\chi_Q = e^{-\ell} = \frac{1}{1+x},
\qquad
x := e^{\ell}-1 = \Lambda_1 f.
```

On the balanced core / compensated outlet,

```math
\chi_Q = \frac{1-9\sigma_c\gamma_c}{1-\sigma_c}.
```

Solving for the odd coefficient gives

```math
\boxed{
\gamma_c(\ell)
=
\frac{1-e^{-\ell}(1-\sigma_c(\ell))}{9\sigma_c(\ell)}.
}
```

Using `x = e^{\ell}-1` this can be rewritten more compactly as

```math
\boxed{
\gamma_c
=
\frac{\sigma_c + x}{9\sigma_c(1+x)},
\qquad
\sigma_c = \sigma_{c,*}(1+x)^{-2/5}.
}
```

The corresponding **bare** odd coefficient is just

```math
\boxed{
\gamma_0 = (1+r_{c,*})\,\gamma_c,
}
```

because the parent hybridization ratio `r_c` is fixed.

So the anomaly does not move the branch in `\mathfrak g` or `\mathfrak r`; it forces only the odd outlet coefficient to detune so that the softened loading share still lands on the electron target.

At first order,

```math
\gamma_c
=
\frac19
+
\frac{1-\sigma_{c,*}}{9\sigma_{c,*}}\,\ell
+
O(\ell^2).
```

So the required shift is indeed small when `\ell` is small.

---

## Step 42E — Same result in compensated Robin–mixed outlet variables

On the compensated hybrid outlet,

```math
\rho_R = 4\sigma_W,
\qquad
\kappa_W = \frac13,
\qquad
\chi_Q = \frac{1-9\sigma_W\gamma_W}{1-\sigma_W}.
```

Identifying `\sigma_W=\sigma_c`, the adiabatic anomaly-only branch gives

```math
\boxed{
\sigma_W(\ell)=\sigma_{W,*}e^{-2\ell/5},
\qquad
\rho_R(\ell)=\rho_{R,*}e^{-2\ell/5},
\qquad
\kappa_W=\frac13,
}
```

and the exact odd detuning law

```math
\boxed{
\gamma_W
=
\frac{\sigma_W + x}{9\sigma_W(1+x)},
\qquad
x=\Lambda_1 f.
}
```

So the hybrid outlet statement is the same as the balanced-core statement:

- even branch fixed,
- parent compensation ratios fixed,
- loading share slightly reduced,
- odd outlet coefficient slightly increased above its canonical value.

---

## Main result of the step

The adiabatic anomaly-only branch has a very sharp parent-side meaning:

```math
\boxed{
\delta\ln\mathfrak g = 0,
\qquad
\delta\ln\mathfrak r = 0,
\qquad
\delta\ln r_c = 0,
\qquad
\delta\ln L_W = 0.
}
```

So the quartic electron sliver is **not** a motion of the parent compensation family or the even D/N geometry.

Instead it is:

```math
\boxed{
\sigma_c(\ell)=\sigma_{c,*}e^{-2\ell/5},
\qquad
\gamma_c
=
\frac{\sigma_c + \Lambda_1 f}{9\sigma_c(1+\Lambda_1 f)},
}
```

or equivalently in compensated Robin–mixed notation,

```math
\boxed{
\rho_R(\ell)=\rho_{R,*}e^{-2\ell/5},
\qquad
\kappa_W=\frac13,
\qquad
\gamma_W
=
\frac{\sigma_W + \Lambda_1 f}{9\sigma_W(1+\Lambda_1 f)}.
}
```

So the improved g-2 picture is now cleaner again:

- the parent compensation sheet stays fixed,
- the even side-channel stays fixed,
- the load share softens slightly,
- and the anomaly is carried only by a tiny odd outlet detuning.
