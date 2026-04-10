# Step 5 — Isotropic one-port loading mismatch and minimal drift families

## Goal

Step 4 rewrote the missing quartic anomaly layer as the weak-axisymmetric outgoing-prefactor slope
```math
\Xi_1=\frac{P_1}{P_0}.
```
The next question is then:

> **on the simplest isotropic one-port branch, what microscopic quantity is `P_1/P_0` actually measuring, and what is the smallest set of reduced drifts that can realize the required quartic target?**

This step answers that.

The main result is that on the weak-axisymmetric grouped branch the anomaly scalar is just the **static loading mismatch**
```math
\Xi_{\rm load}
:=
\frac{P_1}{P_0}
=
\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0},
```
and on the canonical even-preserving branch it is the **only** remaining first-order defect.

So the quartic anomaly problem is no longer “solve all grouped data at once.”
At this level it becomes:

> find a weak static outgoing-transfer strengthening `N_{01}/N_0` and/or a weak static conservative softening `D_{01}/D_0` whose difference equals the required anomaly scalar.

---

## Inputs carried forward

### From the anomaly write-up

The current exact local anomaly closure already reaches
```math
g_{\rm loc}=2.00231930435865,
```
with the remaining miss naturally sitting at
```math
O(f^4),
```
which is exactly where the missing common charge–inertia layer should first enter.
The carried quartic target from Steps 1–4 is
```math
\Lambda_1\approx 0.279605891931464.
```

### From Step 4

The missing quartic layer was rewritten as the weak-axisymmetric outgoing-prefactor slope
```math
\Xi_1=\frac{P_1}{P_0},
```
with the equivalent microscopic slope ledgers
```math
\sigma_Z+2\sigma_\rho-2\sigma_\Omega+2\sigma_{\epsilon_W} = \Xi_1,
```
```math
5\sigma_{c_s}-5\sigma_a-\sigma_{\epsilon_\eta}-\sigma_R = \Xi_1,
```
on the tracking-rigid branch.

### From the later moving-throat grouped notes

On the weak-axisymmetric grouped branch,
```math
P_0^{(A)}=\frac{N_{A,0}}{D_{A,0}},
```
and the first physical prefactor slope is defined by
```math
P_0^{(A)}=P_0+\epsilon\lambda_A P_1+O(\epsilon^2).
```
The same notes give the exact transport identity
```math
\frac{P_1}{P_0} = \frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
```
That is the key simplification used here.

---

## Step 5A — Exact weak-axisymmetric physical slopes

Start from the first-order grouped expansions
```math
D_{A,0}=D_0+\epsilon\lambda_A D_{01},
```
```math
D_{A,2}=D_2+\epsilon\lambda_A D_{21},
```
```math
D_{A,4}=D_4+\epsilon\lambda_A D_{41},
```
```math
N_{A,0}=N_0+\epsilon\lambda_A N_{01}.
```
Then the grouped conservative response and outgoing prefactor are
```math
u_2^{(A)}=-\frac{D_{A,2}}{D_{A,0}},
```
```math
u_4^{(A)}=\frac{D_{A,2}^2-D_{A,0}D_{A,4}}{D_{A,0}^2},
```
```math
P_0^{(A)}=\frac{N_{A,0}}{D_{A,0}}.
```
Expanding to first order gives the exact physical slopes
```math
u_2^{(1)}=-\frac{D_0D_{21}-D_{01}D_2}{D_0^2},
```
```math
u_4^{(1)}=
\frac{D_0(-D_0D_{41}-D_{01}D_4+2D_2D_{21})+2D_{01}(D_0D_4-D_2^2)}{D_0^3},
```
```math
P_1=\frac{D_0N_{01}-N_0D_{01}}{D_0^2},
```
and therefore
```math
\boxed{
\frac{P_1}{P_0}
=
\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
}
```
So the outgoing-prefactor slope is exactly the difference between

- the weak static **outgoing-transfer** slope `N_{01}/N_0`, and
- the weak static **conservative operator** slope `D_{01}/D_0`.

---

## Step 5B — Canonical compensated/even-preserving branch

On the canonical compensated branch,
```math
u_2=\frac19,
\qquad
u_4=\frac{4}{81},
```
so the isotropic conservative moments are
```math
D_2=-\frac{D_0}{9},
\qquad
D_4=-\frac{D_0}{27}.
```
Requiring the grouped conservative response to remain fixed to first order,
```math
u_2^{(1)}=0,
\qquad
u_4^{(1)}=0,
```
forces
```math
\boxed{D_{21}=-\frac{D_{01}}{9},}
```
```math
\boxed{D_{41}=-\frac{D_{01}}{27}.}
```
So on the even-preserving branch the conservative grouped response is fully transported by **one** static slope `D_{01}`.
The even coefficients themselves do not carry an independent first-order defect.

At that point the only remaining first-order grouped outlet scalar is
```math
\boxed{
\Xi_{\rm load}
:=
\frac{P_1}{P_0}
=
\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
}
```
This is the first major simplification of the new PDE-based anomaly rebuild.

---

## Step 5C — Quartic anomaly target as a normalized loading equation

Match the anomaly scalar to the carried quartic target:
```math
\Xi_{\rm load}=\Lambda_1.
```
Define the normalized static drift variables
```math
n:=\frac{N_{01}}{N_0},
\qquad
d:=\frac{D_{01}}{D_0}.
```
Then the anomaly matching law is simply
```math
\boxed{n-d=\Lambda_1.}
```
That is the smallest honest one-port isotropic continuation of the Step-4 quartic constraint.

Three benchmark realizations are immediate.

### 1. Pure outgoing-transfer realization
```math
n=\Lambda_1,
\qquad
d=0.
```
All of the missing quartic loading is placed in the outgoing-transfer strengthening.

### 2. Pure conservative-softening realization
```math
n=0,
\qquad
d=-\Lambda_1.
```
All of the missing quartic loading is placed in the static conservative softening.

### 3. Balanced minimum-norm realization
Minimizing
```math
n^2+d^2
```
subject to
```math
n-d=\Lambda_1
```
gives
```math
\boxed{n=\frac{\Lambda_1}{2},
\qquad
d=-\frac{\Lambda_1}{2}.}
```
So the canonical balanced branch splits the needed loading equally between

- transfer strengthening, and
- conservative softening.

Numerically,
```math
n\approx 0.139802945965732,
\qquad
d\approx -0.139802945965732.
```

---

## Step 5D — Minimal microscopic conservative split

The static conservative operator on the isotropic branch is
```math
D_0=K-B_0-Z_0,
```
so its first normalized weak slope can be written as
```math
d = k-b-z,
```
with
```math
k:=\frac{K_{01}}{D_0},
\qquad
b:=\frac{B_{01}}{D_0},
\qquad
z:=\frac{Z_{01}}{D_0}.
```
Here:
- `k` = wall / geometry stiffness drift,
- `b` = BdG support dressing drift,
- `z` = conservative Maxwell/mixed dressing drift.

Minimizing
```math
k^2+b^2+z^2
```
subject to
```math
k-b-z=d
```
gives the unique minimum-norm conservative split
```math
\boxed{
(k,b,z)=\left(\frac d3,
-\frac d3,
-\frac d3\right).
}
```
On the balanced branch `d=-\Lambda_1/2`, this becomes
```math
\boxed{
(k,b,z)=\left(-\frac{\Lambda_1}{6},
\frac{\Lambda_1}{6},
\frac{\Lambda_1}{6}\right).
}
```
Numerically,
```math
k\approx -0.0466009819885773,
\qquad
b\approx 0.0466009819885773,
\qquad
z\approx 0.0466009819885773.
```
So the most economical conservative-softening realization lowers the bare wall stiffness and raises the support and Maxwell/mixed dressings by equal reduced amounts.

---

## Step 5E — Minimal one-port outgoing-transfer deformation

For a single port-active outgoing block, the later grouped notes reduce the normalized transfer slope to a portwise logarithmic deformation of the form
```math
n = 2\pi_1 - 2\delta_1,
```
where
- `\pi_1` is the port-amplitude logarithmic slope,
- `\delta_1` is the port-denominator logarithmic slope.

Minimizing
```math
\pi_1^2+\delta_1^2
```
subject to
```math
2\pi_1-2\delta_1=n
```
gives
```math
\boxed{
\pi_1=\frac n4,
\qquad
\delta_1=-\frac n4.
}
```
So the most economical one-port transfer deformation splits the outgoing slope symmetrically between

- amplitude strengthening, and
- denominator weakening.

On the balanced anomaly branch `n=\Lambda_1/2`, that becomes
```math
\boxed{
\pi_1=\frac{\Lambda_1}{8},
\qquad
\delta_1=-\frac{\Lambda_1}{8}.
}
```
Numerically,
```math
\pi_1\approx 0.0349507364914330,
\qquad
\delta_1\approx -0.0349507364914330.
```

---

## Step 5F — Final reduced verdict

The old derivation asked for a new common charge–inertia transport layer without a sharp microscopic split.
The new PDE language sharpens that dramatically.

On the simplest isotropic one-port branch,
```math
\boxed{
\Xi_1
=
\frac{P_1}{P_0}
=
\Xi_{\rm load}
=
\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
}
```
And on the canonical even-preserving branch, that is the **only** remaining first-order defect.

So the quartic anomaly problem is no longer

> derive every grouped coefficient.

It is now

> derive one weak static loading mismatch between outgoing-transfer strengthening and conservative static softening.

That is a real simplification, and it gives a much cleaner target for the next moving-throat calculation.
