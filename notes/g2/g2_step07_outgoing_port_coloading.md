# Step 7 — Direct outgoing-port co-loading theorem

## Goal

Step 6 reduced the remaining quartic anomaly closure to the outgoing-load law on the conservative-shape-preserving branch:
```math
\Xi_1
:=
\frac{P_1}{P_0}
=
2\sum_r \rho_r^{(N)}\,\delta\ln\Lambda_r-\kappa_1,
\qquad
\Lambda_r:=\frac{P_r}{\Delta_r}.
```
That was already a strong simplification, but it still phrased the theorem gate in terms of the abstract load factor `\Lambda_r`.

The next honest question is therefore:

> **can the remaining defect be rewritten directly in terms of the weak-axisymmetric slopes of the actual outgoing-port data**
> ```math
> P_r,
> \qquad
> \Delta_r,
> \qquad
> N_{0}^{(r)}=\frac{P_r^2}{\Delta_r^2},
> ```
> **so that the g-2 closure becomes a direct port co-loading law?**

This step answers that.

---

## Important notation warning

There are now **two different `P` symbols** in play.

1. `P_1/P_0` is the **grouped outgoing-prefactor slope** carried from Steps 4–6.
2. `P_r` is the **actual static numerator** of outgoing port `r`.

This step keeps both because the whole point is to prove that the former collapses to a weighted slope law for the latter.

---

## Inputs carried forward

### From Step 6

On the conservative-shape-preserving branch,
```math
\Theta_B=0,
\qquad
\Theta_Z=0,
```
so the surviving grouped defect is
```math
\Xi_1
=
2\sum_r \rho_r^{(N)}\,\delta\ln\Lambda_r-\kappa_1,
\qquad
\sum_r \rho_r^{(N)}=1.
```
Here `\kappa_1` is the weak-axisymmetric wall-baseline slope,
```math
\delta\ln K_A=\epsilon\lambda_A\kappa_1,
```
with grouped branch weights
```math
\lambda_{20}=1,
\qquad
\lambda_{21}=\frac12,
\qquad
\lambda_{22}=-1.
```

### From the moving-throat port dictionary

The actual static outgoing-port data are
```math
P_r=\Omega_{U,r}^2 G_{W,r}+R_r G_{U,r},
```
```math
\Delta_r=\Omega_{U,r}^2\Omega_{W,r}^2-R_r^2,
```
```math
N_0^{(r)}=\frac{P_r^2}{\Delta_r^2}.
```
So the Step-6 load factor is simply
```math
\Lambda_r^2=N_0^{(r)}.
```
That is the bridge used below.

---

## Step 7A — Weak-axisymmetric slopes of the actual outgoing-port data

Define the weak-axisymmetric slopes of the actual numerator and actual detuning by
```math
\delta\ln P_{A,r}=\epsilon\lambda_A\,\mathfrak p_r,
\qquad
\delta\ln \Delta_{A,r}=\epsilon\lambda_A\,\mathfrak d_r.
```
Then the actual static outgoing-transfer coefficient satisfies
```math
N_{A,0}^{(r)}=\frac{P_{A,r}^2}{\Delta_{A,r}^2},
```
so its weak-axisymmetric logarithmic slope is
```math
\boxed{
\delta\ln N_{A,0}^{(r)}=\epsilon\lambda_A\,\nu_r,
\qquad
\nu_r:=2(\mathfrak p_r-\mathfrak d_r).
}
```
This is the first exact identity of the step.

So the wall-referenced outgoing-load defect of one port is just
```math
\delta\ln\!\left(\frac{\Lambda_{A,r}^2}{K_A}\right)
=
\delta\ln N_{A,0}^{(r)}-\delta\ln K_A
=
\epsilon\lambda_A(\nu_r-\kappa_1).
```

---

## Step 7B — Exact collapse of the full remaining grouped defect

Because `\Lambda_{A,r}^2=N_{A,0}^{(r)}`, the Step-6 theorem becomes
```math
\Xi_{\rm load}^{(A)}
=
\sum_r \rho_r^{(N)}
\,\delta\ln\!\left(\frac{\Lambda_{A,r}^2}{K_A}\right)
=
\epsilon\lambda_A
\sum_r \rho_r^{(N)}(\nu_r-\kappa_1).
```
So the scalar carried by the grouped branch is
```math
\boxed{
\Xi_1
=
\sum_r \rho_r^{(N)}(\nu_r-\kappa_1).
}
```
Define the outgoing-weighted static transfer slope
```math
\bar\nu_N:=\sum_r \rho_r^{(N)}\nu_r.
```
Since the weights sum to one, the whole defect collapses to
```math
\boxed{
\Xi_1=\bar\nu_N-\kappa_1.
}
```
This is the main theorem of the step.

It says that the whole remaining linear grouped `2.5`PN / quartic-anomaly defect is exactly the mismatch between

- the outgoing-weighted static transfer slope of the **actual** moving-throat ports,
- and the wall-baseline slope.

So the remaining g-2 theorem problem is no longer “compute every microscopic slippage variable.”
It is just:

> **do the actual outgoing ports co-load with the wall baseline?**

---

## Step 7C — Exact formula for the actual numerator slope `\mathfrak p_r`

The actual static numerator is
```math
P_r=\Omega_{U,r}^2 G_{W,r}+R_r G_{U,r}.
```
Define the positive numerator weights
```math
\boxed{
\alpha_r:=\frac{\Omega_{U,r}^2 G_{W,r}}{P_r},
\qquad
\beta_r:=\frac{R_r G_{U,r}}{P_r},
\qquad
\alpha_r+\beta_r=1.
}
```
Then the weak-axisymmetric numerator slope is exactly
```math
\boxed{
\mathfrak p_r
=
\alpha_r(\mathfrak o_{U,r}+\mathfrak g_{W,r})
+
\beta_r(\mathfrak r_r+\mathfrak g_{U,r}).
}
```
So `\mathfrak p_r` is the convex average of the two static numerator legs:

- the brane-like leg `\Omega_{U,r}^2 G_{W,r}`,
- the mixed leg `R_r G_{U,r}`.

---

## Step 7D — Exact formula for the actual detuning slope `\mathfrak d_r`

The actual static detuning is
```math
\Delta_r=\Omega_{U,r}^2\Omega_{W,r}^2-R_r^2.
```
Define the detuning weights
```math
\boxed{
\chi_r:=\frac{\Omega_{U,r}^2\Omega_{W,r}^2}{\Delta_r},
\qquad
\zeta_r:=\frac{R_r^2}{\Delta_r}.
}
```
Equivalently, with
```math
\mathcal H_r:=\frac{R_r^2}{\Omega_{U,r}^2\Omega_{W,r}^2},
```
one has
```math
\chi_r=\frac{1}{1-\mathcal H_r},
\qquad
\zeta_r=\frac{\mathcal H_r}{1-\mathcal H_r}.
```
Then the weak-axisymmetric detuning slope is exactly
```math
\boxed{
\mathfrak d_r
=
\chi_r(\mathfrak o_{U,r}+\mathfrak o_{W,r})
-
2\zeta_r\,\mathfrak r_r.
}
```
So `\mathfrak d_r` measures how the combined internal frequencies and the coupling `R_r` change the static port detuning.

---

## Step 7E — Static outgoing-transfer slope in actual port variables

Combining the last two results with `\nu_r=2(\mathfrak p_r-\mathfrak d_r)` gives
```math
\boxed{
\nu_r
=
2\alpha_r(\mathfrak o_{U,r}+\mathfrak g_{W,r})
+
2\beta_r(\mathfrak r_r+\mathfrak g_{U,r})
-
2\chi_r(\mathfrak o_{U,r}+\mathfrak o_{W,r})
+
4\zeta_r\,\mathfrak r_r.
}
```
This is the first direct formula for the static outgoing-transfer slope of an actual moving-throat port.

So the exact port-level g-2 closure is now fully explicit in the real port data.

---

## Step 7F — Exact equivalence to the earlier slippage language

Earlier steps used the port slippages
```math
\mathfrak m_r:=\mathfrak g_{W,r}-\mathfrak o_{W,r}-\frac12\kappa_1,
```
```math
\mathfrak i_r:=\mathfrak r_r+\mathfrak g_{U,r}-\mathfrak o_{U,r}-\mathfrak g_{W,r},
```
```math
\mathfrak h_r:=2\mathfrak r_r-\mathfrak o_{U,r}-\mathfrak o_{W,r},
```
with
```math
\mathcal I_r:=\frac{R_r G_{U,r}}{\Omega_{U,r}^2 G_{W,r}},
\qquad
\mathcal H_r:=\frac{R_r^2}{\Omega_{U,r}^2\Omega_{W,r}^2}.
```
Because
```math
\alpha_r=\frac{1}{1+\mathcal I_r},
\qquad
\beta_r=\frac{\mathcal I_r}{1+\mathcal I_r},
\qquad
\chi_r=\frac{1}{1-\mathcal H_r},
\qquad
\zeta_r=\frac{\mathcal H_r}{1-\mathcal H_r},
```
one gets the exact identity
```math
\boxed{
\nu_r
=
\kappa_1
+
2\mathfrak m_r
+
\frac{2\mathcal I_r}{1+\mathcal I_r}\,\mathfrak i_r
+
\frac{2\mathcal H_r}{1-\mathcal H_r}\,\mathfrak h_r.
}
```
So the earlier port amplitude `\sigma_r` is simply
```math
\boxed{
\sigma_r=\nu_r-\kappa_1.
}
```
This proves the present theorem is not a new branch of the algebra. It is the direct outgoing-port rewrite of the earlier slippage theorem.

---

## Step 7G — Outgoing-port co-loading theorem

### Exact zero-defect condition

The remaining linear grouped defect vanishes if and only if
```math
\boxed{
\bar\nu_N=\kappa_1.
}
```
Equivalently,
```math
\boxed{
\frac{P_1}{P_0}=0
\quad\Longleftrightarrow\quad
\bar\nu_N=\kappa_1.
}
```
So the outgoing-weighted static transfer slope must match the wall-baseline slope.

### Strong per-port sufficient condition

A stronger sufficient condition is
```math
\boxed{
\nu_r=\kappa_1
\qquad
\text{for every active outgoing port }r.
}
```
Then every static outgoing-transfer coefficient co-loads with the wall baseline lane by lane, and therefore `\Xi_1=0`.

### Dominant-port limit

If one outgoing port dominates,
```math
\rho_{r_*}^{(N)}\approx1,
```
then
```math
\boxed{
\Xi_1\approx \nu_{r_*}-\kappa_1.
}
```
So the last linear grouped defect is just the mismatch between the dominant port slope and the wall-baseline slope.

### Naive rigidity no-go

If the actual ports are rigid,
```math
\mathfrak p_r=\mathfrak d_r=0
\qquad\Longrightarrow\qquad
\nu_r=0,
```
then
```math
\boxed{
\Xi_1=-\kappa_1.
}
```
So freezing the outgoing ports is **not** enough. They must actively co-load with the wall baseline.

---

## Step 7H — Direct quartic anomaly law

The carried quartic target is
```math
\Xi_1=\Lambda_1,
\qquad
\Lambda_1\approx 0.279605891931464.
```
So the exact port-level anomaly gate is now
```math
\boxed{
\bar\nu_N=\kappa_1+\Lambda_1.
}
```
In the dominant-port regime this becomes
```math
\boxed{
\nu_{r_*}=\kappa_1+\Lambda_1.
}
```
This is the cleanest g-2 continuation reached so far:

> the remaining quartic anomaly correction is exactly the amount by which the outgoing-weighted static transfer slope must exceed the wall-baseline slope.

---

## Step 7I — Reduced verdict

Step 6 said the whole remaining defect is an outgoing-load law for `\Lambda_r=P_r/\Delta_r`.
Step 7 sharpens that again.

It proves that:

1. the actual static outgoing-transfer slope of each port is
   ```math
   \nu_r=2(\mathfrak p_r-\mathfrak d_r),
   ```
2. the full remaining grouped defect is exactly
   ```math
   \Xi_1=\bar\nu_N-\kappa_1,
   ```
3. the earlier slippage amplitude is just
   ```math
   \sigma_r=\nu_r-\kappa_1,
   ```
4. and anomaly matching requires
   ```math
   \bar\nu_N=\kappa_1+\Lambda_1.
   ```

So the next theorem gate is smaller again.

It is no longer
> “compute the whole microscopic slippage bundle.”

It is now simply
> **compute the actual outgoing-port slopes `\nu_r` on the moving-throat branch and see whether their outgoing-weighted average lands at `\kappa_1+\Lambda_1`.**
