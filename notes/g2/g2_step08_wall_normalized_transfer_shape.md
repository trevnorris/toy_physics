# Step 8 — Wall-normalized transfer-shape theorem

## Goal

Step 7 reduced the remaining quartic anomaly closure to the actual outgoing-port co-loading law
```math
\Xi_1
=
\bar\nu_N-\kappa_1,
\qquad
\bar\nu_N:=\sum_r \rho_r^{(N)}\nu_r,
```
where
```math
\nu_r
=
\frac{\delta\ln N_{A,0}^{(r)}}{\epsilon\lambda_A},
\qquad
N_{A,0}^{(r)}=\frac{P_{A,r}^2}{\Delta_{A,r}^2}.
```
That was already sharp, but it still phrased the theorem gate in terms of the raw actual-port pair
```math
P_r,
\qquad
\Delta_r.
```

The next honest question is therefore:

> **does each actual outgoing port admit a direct wall-normalized factorization, so that the remaining grouped defect becomes nothing more than the weak-axisymmetric drift of one dimensionless transfer shape?**

This step answers that.

---

## Inputs carried forward

### From Step 7

For each actual outgoing port,
```math
P_r=\Omega_{U,r}^2 G_{W,r}+R_r G_{U,r},
```
```math
\Delta_r=\Omega_{U,r}^2\Omega_{W,r}^2-R_r^2,
```
```math
N_0^{(r)}=\frac{P_r^2}{\Delta_r^2}.
```
On the weak-axisymmetric grouped branch,
```math
\delta\ln K_A=\epsilon\lambda_A\kappa_1,
\qquad
\lambda_{20}=1,
\quad
\lambda_{21}=\frac12,
\quad
\lambda_{22}=-1,
```
and Step 7 proved
```math
\Xi_1
=
\sum_r \rho_r^{(N)}(\nu_r-\kappa_1),
\qquad
\nu_r
=
\frac{\delta\ln N_{A,0}^{(r)}}{\epsilon\lambda_A}.
```

### From the moving-throat port dictionary

The wall baseline is the static grouped prefactor carrier `K_A`, so the natural next move is to normalize every actual outgoing quantity by `K_A`, `\Omega_{U,r}`, and `\Omega_{W,r}` before taking grouped weak-axisymmetric slopes.

---

## Step 8A — Exact wall-normalized factorization of the actual outgoing-transfer coefficient

Introduce the wall-normalized dimensionless port variables
```math
\boxed{
\widehat G_{W,r}:=\frac{G_{W,r}}{\Omega_{W,r}^2\sqrt K},
\qquad
\widehat G_{U,r}:=\frac{G_{U,r}}{\Omega_{U,r}\Omega_{W,r}\sqrt K},
\qquad
\widehat R_r:=\frac{R_r}{\Omega_{U,r}\Omega_{W,r}}.
}
```
Then the actual numerator and detuning factor exactly as
```math
P_r
=
\sqrt K\,\Omega_{U,r}^2\Omega_{W,r}^2
\bigl(\widehat G_{W,r}+\widehat R_r\widehat G_{U,r}\bigr),
```
```math
\Delta_r
=
\Omega_{U,r}^2\Omega_{W,r}^2(1-\widehat R_r^2).
```
Therefore
```math
\boxed{
\frac{N_0^{(r)}}{K}
=
\left[
\frac{\widehat G_{W,r}+\widehat R_r\widehat G_{U,r}}{1-\widehat R_r^2}
\right]^2.
}
```
Define the wall-normalized transfer shape
```math
\boxed{
\mathcal T_r
:=
\frac{\widehat G_{W,r}+\widehat R_r\widehat G_{U,r}}{1-\widehat R_r^2},
}
```
so that the actual static outgoing-transfer coefficient factors exactly as
```math
\boxed{
N_0^{(r)}=K\,\mathcal T_r^2.
}
```

This is the first main result of the step.
The actual outgoing-transfer coefficient is the wall baseline times one dimensionless transfer shape squared.

---

## Step 8B — Weak-axisymmetric transport of the wall-normalized port variables

Define the weak-axisymmetric grouped slopes of the wall-normalized port variables by
```math
\boxed{
\delta\ln\widehat G_{W,A,r}=\epsilon\lambda_A\,\mathfrak w_r,
}
```
```math
\boxed{
\delta\ln\widehat G_{U,A,r}=\epsilon\lambda_A\,\mathfrak u_r,
}
```
```math
\boxed{
\delta\ln\widehat R_{A,r}=\epsilon\lambda_A\,\mathfrak c_r.
}
```
In terms of the primitive Step-7 slopes,
```math
\boxed{
\mathfrak w_r
=
\mathfrak g_{W,r}-\mathfrak o_{W,r}-\frac12\kappa_1,
}
```
```math
\boxed{
\mathfrak u_r
=
\mathfrak g_{U,r}-\frac12\mathfrak o_{U,r}-\frac12\mathfrak o_{W,r}-\frac12\kappa_1,
}
```
```math
\boxed{
\mathfrak c_r
=
\mathfrak r_r-\frac12\mathfrak o_{U,r}-\frac12\mathfrak o_{W,r}.
}
```
So the wall-normalized variables already strip out the trivial wall-baseline scaling and leave only the genuine dimensionless port-shape drifts.

---

## Step 8C — Exact transfer-shape slope and the identity `\nu_r=\kappa_1+2\tau_r`

Because
```math
\mathcal T_r
=
\frac{\widehat G_{W,r}+\widehat R_r\widehat G_{U,r}}{1-\widehat R_r^2},
```
its weak-axisymmetric logarithmic slope is
```math
\boxed{
\delta\ln\mathcal T_{A,r}=\epsilon\lambda_A\,\tau_r,
}
```
with
```math
\boxed{
\tau_r
=
\widehat\alpha_r\,\mathfrak w_r
+
\widehat\beta_r\,(\mathfrak u_r+\mathfrak c_r)
+
\frac{2\widehat R_r^2}{1-\widehat R_r^2}\,\mathfrak c_r,
}
```
where
```math
\boxed{
\widehat\alpha_r:=
\frac{\widehat G_{W,r}}{\widehat G_{W,r}+\widehat R_r\widehat G_{U,r}},
\qquad
\widehat\beta_r:=
\frac{\widehat R_r\widehat G_{U,r}}{\widehat G_{W,r}+\widehat R_r\widehat G_{U,r}},
\qquad
\widehat\alpha_r+\widehat\beta_r=1.
}
```

Now use the exact factorization
```math
N_{A,0}^{(r)}=K_A\mathcal T_{A,r}^2.
```
Taking the weak-axisymmetric logarithmic slope gives the central identity
```math
\boxed{
\nu_r
=
\frac{\delta\ln N_{A,0}^{(r)}}{\epsilon\lambda_A}
=
\kappa_1+2\tau_r.
}
```
So the actual static outgoing-transfer slope is just the wall-baseline slope plus twice the transfer-shape slope.

---

## Step 8D — Exact collapse of the remaining grouped defect

Step 7 already gave
```math
\Xi_1
=
\sum_r \rho_r^{(N)}(\nu_r-\kappa_1),
\qquad
\sum_r \rho_r^{(N)}=1.
```
Substituting the transfer-shape identity above yields
```math
\boxed{
\Xi_1
=
2\sum_r \rho_r^{(N)}\tau_r.
}
```
Equivalently, the exact zero-defect condition becomes
```math
\boxed{
\sum_r \rho_r^{(N)}\tau_r=0.
}
```
This is the sharpest reduced theorem gate reached so far.

A stronger per-port sufficient condition is
```math
\boxed{
\tau_r=0
\qquad\text{for every active outgoing port }r,
}
```
which is equivalent to
```math
\boxed{
\delta\ln\mathcal T_{A,r}=0
\qquad\text{for every active outgoing port }r.
}
```
So the exact reduced meaning of port co-loading is now:

> **each wall-normalized transfer shape must be weak-axisymmetrically rigid.**

---

## Step 8E — Exact equivalence to the earlier slippage languages

Earlier steps used the slippage variables
```math
\mathfrak m_r,
\qquad
\mathfrak i_r,
\qquad
\mathfrak h_r,
```
with the port-shape composites
```math
\mathcal M_r:=\frac{G_{W,r}}{\Omega_{W,r}^2\sqrt K},
\qquad
\mathcal I_r:=\frac{R_r G_{U,r}}{\Omega_{U,r}^2 G_{W,r}},
\qquad
\mathcal H_r:=\frac{R_r^2}{\Omega_{U,r}^2\Omega_{W,r}^2}.
```
These are related to the present variables by
```math
\mathcal M_r=\widehat G_{W,r},
\qquad
\mathcal I_r=\frac{\widehat R_r\widehat G_{U,r}}{\widehat G_{W,r}},
\qquad
\mathcal H_r=\widehat R_r^2,
```
and
```math
\mathfrak m_r=\mathfrak w_r,
\qquad
\mathfrak i_r=(\mathfrak u_r+\mathfrak c_r)-\mathfrak w_r,
\qquad
\mathfrak h_r=2\mathfrak c_r.
```
With these substitutions the transfer-shape slope becomes exactly
```math
\boxed{
\tau_r
=
\mathfrak m_r
+
\frac{\mathcal I_r}{1+\mathcal I_r}\,\mathfrak i_r
+
\frac{\mathcal H_r}{1-\mathcal H_r}\,\mathfrak h_r.
}
```
So the earlier Step-7 port amplitude is simply
```math
\boxed{
\sigma_r=2\tau_r.
}
```
This proves the present theorem is not a different branch of the algebra. It is the exact compressed form of the Stage-159/160/161 slippage structure.

---

## Step 8F — Direct quartic anomaly law in transfer-shape language

The carried quartic target is still
```math
\Xi_1=\Lambda_1,
\qquad
\Lambda_1\approx 0.279605891931464.
```
Since
```math
\Xi_1=2\sum_r \rho_r^{(N)}\tau_r,
```
the exact transfer-shape anomaly gate is now
```math
\boxed{
\sum_r \rho_r^{(N)}\tau_r
=
\frac{\Lambda_1}{2}
\approx 0.139802945965732.
}
```
In the dominant-port limit,
```math
\boxed{
\tau_{r_*}
=
\frac{\Lambda_1}{2}
\approx 0.139802945965732.
}
```
So the remaining quartic anomaly correction is exactly the amount by which the outgoing-weighted transfer shape must drift away from weak-axisymmetric rigidity.

---

## Step 8G — Reduced verdict

Step 7 said the remaining quartic anomaly layer is a mismatch between the outgoing-weighted actual transfer slope and the wall-baseline slope.
Step 8 sharpens that one step further.

It proves that:

1. every actual outgoing-transfer coefficient factors exactly as
   ```math
   N_0^{(r)}=K\,\mathcal T_r^2,
   ```
2. the actual port slope is
   ```math
   \nu_r=\kappa_1+2\tau_r,
   ```
3. the full remaining grouped defect is exactly
   ```math
   \Xi_1=2\sum_r \rho_r^{(N)}\tau_r,
   ```
4. and anomaly matching requires
   ```math
   \sum_r \rho_r^{(N)}\tau_r=\Lambda_1/2.
   ```

So the next theorem gate is smaller again.

It is no longer
> “compute the raw outgoing-port slopes `\nu_r`.”

It is now simply
> **compute the weak-axisymmetric wall-normalized transfer-shape slopes `\tau_r` on the actual moving-throat branch, and test whether their outgoing-weighted average lands at `\Lambda_1/2`.**

That is the direct continuation point.
