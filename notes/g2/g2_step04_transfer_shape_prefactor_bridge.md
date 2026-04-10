# Step 4 — Exact transfer-shape / outgoing-prefactor bridge

## Goal

Step 3 fixed the coherent weak-axisymmetric defect to the exact branch-normal-form scalar
```math
\Xi_1=A_{{\rm tr},*}q_{\rm tr}+q_{\rm nt},
```
and separated the corrected nontracking coordinate
```math
q_{\rm nt}=\delta\ln \mathfrak N_*.
```

The next question is the obvious one:

> **what microscopic moving-throat quantity is `\Xi_1` actually measuring?**

This step answers that. It shows that `\Xi_1` is simultaneously

1. the logarithmic slope of the raw transfer shape `\mathcal T^2`,
2. the logarithmic slope of the effective outgoing prefactor `P_A`,
3. and therefore the first direct microscopic slope that can carry the missing quartic anomaly layer.

So this is the first point where the anomaly remainder is tied directly to actual branch data rather than to an abstract quotient tangent.

---

## Inputs carried forward

### From Step 3

On the coherent local branch,
```math
\Theta_1=-C_{{\rm tr},*}q_{\rm tr},
```
```math
\Xi_1=A_{{\rm tr},*}q_{\rm tr}+q_{\rm nt},
```
with
```math
A_{{\rm tr},*}=
\frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})},
\qquad
C_{{\rm tr},*}=
\frac{\chi_{0,*}\delta_{U,*}}
{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}.
```

The exact corrected nontracking branch composite was
```math
\mathfrak N_*:=\mathcal T^2 R_{\rm tr}^{B_*},
\qquad
B_*=
\frac{2(1+\chi_{0,*}+\delta_{U,*})}{\delta_{U,*}}.
```

### From Step 2 / Step 3

The quartic target carried forward from the already-frozen local anomaly law is
```math
\Lambda_1\approx 0.279605891931464.
```

At Step 3 there were still two logically distinct ways to use it:

- **raw coherent-defect closure**
  ```math
  \Xi_1=\Lambda_1,
  ```
- **corrected nontracking closure**
  ```math
  q_{\rm nt}=\Lambda_1.
  ```

Step 4 shows exactly what those two choices mean microscopically.

---

## Step 4A — The corrected branch composite removes the tracking feed-through exactly

Using
```math
\delta\ln R_{\rm tr}=-C_{{\rm tr},*}q_{\rm tr},
```
and
```math
\delta\ln\mathfrak N_*
=
\delta\ln\mathcal T^2 + B_*\,\delta\ln R_{\rm tr},
```
we get
```math
\delta\ln\mathfrak N_*
=
\delta\ln\mathcal T^2 - B_* C_{{\rm tr},*} q_{\rm tr}.
```
But the branch coefficients satisfy the exact identity
```math
A_{{\rm tr},*}=B_* C_{{\rm tr},*}.
```
So
```math
\delta\ln\mathcal T^2
=
A_{{\rm tr},*}q_{\rm tr}+q_{\rm nt}
=
\Xi_1.
```
Therefore

```math
\boxed{\delta\ln\mathcal T^2=\Xi_1,}
```

while the corrected composite gives

```math
\boxed{\delta\ln\mathfrak N_*=q_{\rm nt}.}
```

So the two Step-3 scalars are now interpreted very sharply:

- `\Xi_1` = **raw transfer-shape slope**,
- `q_nt` = **corrected nontracking transfer-shape slope**.

---

## Step 4B — `\Xi_1` is also the outgoing-prefactor slope

The later moving-throat grouped-`P2` notes show that the weak-axisymmetric grouped defect can be written as the logarithmic slope of the effective transfer shape,
```math
\Xi_1
=
\frac{\delta\ln \mathcal T_{{\rm eff},A}^2}{\epsilon\lambda_A},
```
and equally as the logarithmic slope of the outgoing prefactor,
```math
\Xi_1=\frac{P_1}{P_0},
```
if
```math
P_A=P_0+\epsilon\lambda_A P_1+O(\epsilon^2).
```

So the coherent defect is not just a bookkeeping scalar. It is literally the **first normalized slope of the outgoing quadrupole bridge**.

On a one-port branch,
```math
\mathcal T_A=\mathcal T_{A,0}\,e^{\epsilon\lambda_A\tau_{\rm eff}}+O(\epsilon^2),
```
so
```math
\boxed{\Xi_1=2\tau_{\rm eff}.}
```

This is the cleanest bridge yet between the conservative quotient story and the outgoing branch story.

---

## Step 4C — Exact microscopic one-port slope formulas

On the actual minimal one-port continuum branch,
```math
\mathcal T_A^2
=
\frac{Z_{W,A}(1+\rho_A)^2}
{\Omega_{W,A}^2(1-\epsilon_{W,A})^2}
=
\frac{27\pi^2Gc_s^5}{20a^5c^5}
\frac{1-\epsilon_{\eta,A}}{R_{{\rm target},A}}.
```

Taking the weak-axisymmetric logarithmic slope gives two exact ledgers for the same scalar `\Xi_1`.

### Port-side mixed-sector ledger

Define the dimensionless slopes
```math
\sigma_Z:=\delta\ln Z_W,
\qquad
\sigma_\Omega:=\delta\ln\Omega_W,
```
```math
\sigma_\rho:=\frac{\delta\rho}{1+\rho},
\qquad
\sigma_{\epsilon_W}:=\frac{\delta\epsilon_W}{1-\epsilon_W}.
```
Then
```math
\boxed{
\Xi_1
=
\sigma_Z + 2\sigma_\rho - 2\sigma_\Omega + 2\sigma_{\epsilon_W}.
}
```

### Geometric / selected-branch ledger

Define
```math
\sigma_{c_s}:=\delta\ln c_s,
\qquad
\sigma_a:=\delta\ln a,
\qquad
\sigma_R:=\delta\ln R_{\rm target},
```
```math
\sigma_{\epsilon_\eta}:=\frac{\delta\epsilon_\eta}{1-\epsilon_\eta}.
```
Then
```math
\boxed{
\Xi_1
=
5\sigma_{c_s} - 5\sigma_a - \sigma_{\epsilon_\eta} - \sigma_R.
}
```

So the same coherent defect can be read in two complementary microscopic ways:

- as a mixed-port overlap / blocking / frequency slope,
- or as a sound-speed / throat-size / selected-target / dressing slope.

This is exactly the kind of simplification the old staggered derivation did not have.

---

## Step 4D — Quartic anomaly target in direct microscopic variables

Now combine the outgoing-prefactor identity
```math
\Xi_1=\frac{P_1}{P_0}
```
with the Step-3 quartic target.

### 1. Direct coherent-defect closure

If the missing quartic anomaly layer is identified with the raw coherent defect, then
```math
\boxed{\frac{P_1}{P_0}=\Lambda_1.}
```
So the required prefactor slope is directly fixed:
```math
\boxed{\frac{P_1}{P_0}\approx 0.279605891931464.}
```

### 2. Corrected nontracking closure

If the physically correct common scalar is instead the corrected nontracking branch invariant,
```math
q_{\rm nt}=\Lambda_1,
```
then because
```math
q_{\rm nt}=\Xi_1-A_{{\rm tr},*}q_{\rm tr},
```
we get
```math
\boxed{
\frac{P_1}{P_0} - A_{{\rm tr},*} s_{\rm tr} = \Lambda_1.
}
```
So the raw prefactor slope and the corrected nontracking slope differ only by the universal tracking feed-through.

### 3. Tracking-rigid branch

If the coherent branch is tracking-rigid at this order,
```math
s_{\rm tr}=0,
```
then the two closures coincide:
```math
\boxed{
\frac{P_1}{P_0}
=
\Xi_1
=
q_{{\rm nt},1}
=
\Lambda_1.
}
```

That is the cleanest first microscopic closure of the missing quartic anomaly layer.

---

## Step 4E — Final microscopic matching equations

Using the port-side and geometric slope ledgers, the corrected quartic matching law becomes

```math
\boxed{
\sigma_Z + 2\sigma_\rho - 2\sigma_\Omega + 2\sigma_{\epsilon_W}
=
\Lambda_1 + A_{{\rm tr},*} s_{\rm tr},
}
```

and equivalently

```math
\boxed{
5\sigma_{c_s} - 5\sigma_a - \sigma_{\epsilon_\eta} - \sigma_R
=
\Lambda_1 + A_{{\rm tr},*} s_{\rm tr}.
}
```

On the tracking-rigid branch these reduce to

```math
\boxed{
\sigma_Z + 2\sigma_\rho - 2\sigma_\Omega + 2\sigma_{\epsilon_W}
=
5\sigma_{c_s} - 5\sigma_a - \sigma_{\epsilon_\eta} - \sigma_R
=
0.279605891931464.
}
```

This is the strongest reduced g-2 statement reached so far in the new PDE language:

> the missing quartic anomaly layer is the required weak-axisymmetric logarithmic slope of the effective outgoing transfer shape, equivalently of the outgoing prefactor, and it can now be written directly in microscopic moving-throat variables.

---

## Why this step matters

The old anomaly derivation stalled because the final correction still looked like “one more transport term.”
Step 4 removes that vagueness.

The remaining problem is no longer:

> invent the quartic coefficient.

It is now:

> compute the actual branch slopes of
> `Z_W`, `\rho`, `\Omega_W`, `\epsilon_W`, `c_s`, `a`, `\epsilon_\eta`, and `R_{\rm target}`
> on the moving-throat branch, and see whether their exact combination lands on
> `0.279605891931464`.

That is the first point where the new PDE framework genuinely improves the anomaly calculation rather than merely rephrasing it.
