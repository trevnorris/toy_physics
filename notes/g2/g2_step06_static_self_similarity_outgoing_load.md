# Step 6 — Static self-similarity and the outgoing-load theorem

## Goal

Step 5 showed that on the canonical even-preserving isotropic one-port branch the missing quartic anomaly layer is just one scalar,
```math
\Xi_{\rm load}
:=
\frac{P_1}{P_0}
=
\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
```
That was already a major simplification, but it still left the two slopes
```math
D_{01}/D_0,
\qquad
N_{01}/N_0
```
looking like unrelated microscopic inputs.

The next honest question is therefore:

> **is the remaining grouped loading defect a generic static mismatch of the whole moving-throat bundle, or does it collapse to a much narrower wall-referenced theorem?**

This step answers that.

The main result is that the entire remaining linear grouped defect is exactly a **weighted failure of static self-similarity relative to the wall baseline**. After the natural wall-normalized factorization, the conservative bundles reduce to pure shape variables, while the outgoing bundle reduces to a wall-loading law. On conservative-shape-preserving branches, the anomaly target becomes a pure **outgoing-load theorem**.

---

## Inputs carried forward

### From Step 5

On the canonical even-preserving branch,
```math
u_2^{(1)}=0,
\qquad
u_4^{(1)}=0
```
forces
```math
D_{21}=-\frac{D_{01}}{9},
\qquad
D_{41}=-\frac{D_{01}}{27},
```
so the grouped conservative response is transported by one static slope `D_{01}`.
The only remaining first-order outlet scalar is then
```math
\Xi_{\rm load}
=
\frac{P_1}{P_0}
=
\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
```
Matching the anomaly requires
```math
\Xi_{\rm load}=\Lambda_1,
\qquad
\Lambda_1\approx 0.279605891931464.
```

### From the later moving-throat grouped notes

The grouped weak-axisymmetric notes sharpen that scalar further: the static operator and outgoing-transfer slopes can themselves be rewritten in wall-referenced logarithmic form, and the whole remaining grouped loading defect can be expressed as a weighted failure of static self-similarity between

- the wall baseline,
- the BdG support bundle,
- the conservative Maxwell/mixed bundle,
- and the outgoing-transfer bundle.

That is the language adopted here.

---

## Step 6A — Exact wall-referenced decomposition of the loading scalar

Start from the static conservative operator
```math
D_0=K-B_0-Z_0,
\qquad
D_{01}=K_1-B_0^{(1)}-Z_0^{(1)}.
```
Define the logarithmic slopes
```math
\delta_K:=\frac{K_1}{K},
\qquad
\delta_B:=\frac{B_0^{(1)}}{B_0},
\qquad
\delta_Z:=\frac{Z_0^{(1)}}{Z_0},
\qquad
\delta_N:=\frac{N_1}{N_0},
```
and the static weights
```math
\omega_K:=\frac{K}{D_0},
\qquad
\omega_B:=\frac{B_0}{D_0},
\qquad
\omega_Z:=\frac{Z_0}{D_0}.
```
Then
```math
\delta_D:=\frac{D_{01}}{D_0}
=
\omega_K\delta_K-\omega_B\delta_B-\omega_Z\delta_Z,
```
with the exact identity
```math
\omega_K-\omega_B-\omega_Z=1.
```
Therefore the Step-5 scalar becomes
```math
\Xi_{\rm load}
=
\delta_N-\delta_D
=
(\delta_N-\delta_K)
+
\omega_B(\delta_B-\delta_K)
+
\omega_Z(\delta_Z-\delta_K).
```
So the natural reference slope is not arbitrary. It is the wall-baseline slope `\delta_K`.

This is the first sharpened result of Step 6:

> the remaining quartic anomaly layer is the weighted failure of the three support/transfer sectors to co-load with the wall baseline.

---

## Step 6B — Microscopic weighted-log forms of the three sectors

The next step is to rewrite the three sector slopes as explicit weighted logarithmic drifts.

### BdG support bundle

For each support mode `\alpha`,
```math
B_{0,\alpha}=\frac{c_\alpha^2}{\varpi_\alpha^2},
\qquad
B_0=\sum_\alpha B_{0,\alpha},
```
so
```math
\delta_B
=
\sum_\alpha \rho_\alpha^{(B)}
\,2\,\delta\ln\!\left(\frac{c_\alpha}{\varpi_\alpha}\right),
\qquad
\rho_\alpha^{(B)}:=\frac{B_{0,\alpha}}{B_0}.
```

### Conservative Maxwell/mixed static bundle

For each port `r`,
```math
Z_0^{(r)}=\frac{Q_r}{\Delta_r},
\qquad
Z_0=\sum_r Z_0^{(r)},
```
so
```math
\delta_Z
=
\sum_r \rho_r^{(Z)}
\,\delta\ln\!\left(\frac{Q_r}{\Delta_r}\right),
\qquad
\rho_r^{(Z)}:=\frac{Z_0^{(r)}}{Z_0}.
```

### Outgoing-transfer bundle

For each outgoing port,
```math
N_0^{(r)}=\frac{P_r^2}{\Delta_r^2},
\qquad
N_0=\sum_r N_0^{(r)},
```
so
```math
\delta_N
=
\sum_r \rho_r^{(N)}
\,2\,\delta\ln\!\left(\frac{P_r}{\Delta_r}\right),
\qquad
\rho_r^{(N)}:=\frac{N_0^{(r)}}{N_0}.
```

So the Step-5 loading scalar is already controlled by three weighted microscopic log drifts.

---

## Step 6C — Exact self-similarity defect fields

Because `\delta_K` is the natural wall-baseline slope, define the wall-referenced defect fields
```math
\Sigma_\alpha^{(B)}
:=
2\,\delta\ln\!\left(\frac{c_\alpha}{\varpi_\alpha}\right)-\delta_K,
```
```math
\Sigma_r^{(Z)}
:=
\delta\ln\!\left(\frac{Q_r}{\Delta_r}\right)-\delta_K,
```
```math
\Sigma_r^{(N)}
:=
2\,\delta\ln\!\left(\frac{P_r}{\Delta_r}\right)-\delta_K.
```
Then the whole grouped defect becomes
```math
\boxed{
\Xi_{\rm load}
=
\sum_r \rho_r^{(N)}\Sigma_r^{(N)}
+
\omega_B\sum_\alpha \rho_\alpha^{(B)}\Sigma_\alpha^{(B)}
+
\omega_Z\sum_r \rho_r^{(Z)}\Sigma_r^{(Z)}.
}
```

This is the sharpest exact Step-6 formula.

It says the remaining linear grouped anomaly defect is not sourced by the whole microscopic bundle independently. It is sourced only by the weighted failure of

1. the BdG support bundle,
2. the conservative Maxwell/mixed bundle,
3. and the outgoing-transfer bundle,

to co-load with the wall baseline.

---

## Step 6D — Wall-normalized factorization into shape and load variables

The previous defect fields sharpen further after factoring by the wall baseline.

### BdG support shape

Define
```math
\chi_\alpha:=\frac{c_\alpha}{\sqrt K\,\varpi_\alpha}.
```
Then
```math
B_{0,\alpha}=K\,\chi_\alpha^2,
```
and the defect field becomes
```math
\Sigma_\alpha^{(B)}=\delta\ln\chi_\alpha^2.
```
So the BdG sector contributes only through the drift of the wall-normalized support shape `\chi_\alpha`.

### Conservative port shape and outgoing load

Define the wall-normalized port variables
```math
\Upsilon_r:=\frac{Q_r}{K\Delta_r},
\qquad
\Lambda_r:=\frac{P_r}{\Delta_r}.
```
Then
```math
Z_0^{(r)}=K\,\Upsilon_r,
\qquad
N_0^{(r)}=\Lambda_r^2,
```
so
```math
\Sigma_r^{(Z)}=\delta\ln\Upsilon_r,
```
```math
\Sigma_r^{(N)}=\delta\ln\!\left(\frac{\Lambda_r^2}{K}\right)=2\,\delta\ln\Lambda_r-\delta_K.
```

This is the second major Step-6 sharpening:

- the conservative bundles reduce to pure **shape** drifts,
- the outgoing bundle reduces to a **wall-loading law** for `\Lambda_r`.

Define the weighted defect measures
```math
\Theta_B:=\sum_\alpha \rho_\alpha^{(B)}\,\delta\ln\chi_\alpha^2,
```
```math
\Theta_Z:=\sum_r \rho_r^{(Z)}\,\delta\ln\Upsilon_r,
```
```math
\Theta_N:=\sum_r \rho_r^{(N)}\,\delta\ln\!\left(\frac{\Lambda_r^2}{K}\right).
```
Then the whole grouped defect collapses to the compact formula
```math
\boxed{
\Xi_{\rm load}=\Theta_N+\omega_B\Theta_B+\omega_Z\Theta_Z.
}
```

---

## Step 6E — Conservative-shape theorem and outgoing-load theorem

Now specialize to the branch on which the conservative shapes are preserved:
```math
\delta\ln\chi_\alpha^2=0
\quad\text{for all }\alpha,
```
```math
\delta\ln\Upsilon_r=0
\quad\text{for all }r.
```
Then
```math
\Theta_B=0,
\qquad
\Theta_Z=0,
```
and the remaining grouped defect collapses to
```math
\boxed{
\Xi_{\rm load}
=
\sum_r \rho_r^{(N)}\left(2\,\delta\ln\Lambda_r-\delta_K\right).
}
```

So on conservative-shape-preserving branches, the full remaining linear grouped `2.5`PN defect is carried **only** by the outgoing load factor `\Lambda_r`.

That is the main theorem of the step.

---

## Step 6F — Two immediate consequences

### 1. Naive common-self-similarity no-go

If one also freezes the outgoing load factor itself,
```math
\delta\ln\Lambda_r=0
\quad\text{for all }r,
```
then
```math
\Xi_{\rm load}=-\delta_K.
```
So a naive common self-similarity branch does **not** kill the remaining grouped defect unless the wall baseline itself does not load.

This is a real no-go result.

### 2. Exact outgoing-load law

The defect vanishes if and only if
```math
\boxed{
2\sum_r \rho_r^{(N)}\,\delta\ln\Lambda_r
=
\delta_K.
}
```
A stronger sufficient condition is the portwise law
```math
2\,\delta\ln\Lambda_r=\delta_K
\qquad\text{for every outgoing port }r.
```
So the old vague requirement “the outgoing bundle must be self-similar” is too weak.
The actual requirement is a precise wall-loading law.

---

## Step 6G — Direct quartic anomaly law

Now feed the Step-5 anomaly target back in:
```math
\Xi_{\rm load}=\Lambda_1.
```
On a conservative-shape-preserving branch this becomes the exact quartic theorem gate
```math
\boxed{
2\sum_r \rho_r^{(N)}\,\delta\ln\Lambda_r
-
\delta_K
=
\Lambda_1.
}
```
So the anomaly does **not** ask for an arbitrary new static bundle deformation.
It asks for a specific outgoing-load mismatch relative to the wall baseline.

If the outgoing ports all share a common weak logarithmic slope,
```math
\delta\ln\Lambda_r=\ell,
```
then
```math
\Xi_{\rm load}=2\ell-\delta_K,
```
so
```math
\ell_{\rm kill}=\frac{\delta_K}{2}
```
would kill the grouped defect,
and anomaly matching instead requires
```math
\boxed{
\ell_{\rm anom}
=
\frac{\delta_K+\Lambda_1}{2}.
}
```
Therefore the anomaly target requires an **extra outgoing-load slippage above the defect-killing wall-tracking law** of
```math
\boxed{
\ell_{\rm anom}-\ell_{\rm kill}
=
\frac{\Lambda_1}{2}.
}
```
Numerically,
```math
\frac{\Lambda_1}{2}
\approx 0.139802945965732.
```

This is the cleanest continuation of Step 5 reached so far.

---

## Step 6H — Final reduced verdict

Step 5 showed that the quartic anomaly correction is one scalar loading mismatch.
Step 6 sharpens that again:

1. `\Xi_{\rm load}` is exactly a wall-referenced self-similarity defect.
2. After wall-normalization, the conservative bundles become pure shape variables `\chi_\alpha` and `\Upsilon_r`.
3. On conservative-shape-preserving branches, the whole remaining linear grouped defect is carried only by the outgoing load factor `\Lambda_r=P_r/\Delta_r`.
4. So the quartic anomaly target is not “derive every grouped coefficient.”
5. It is:
   > determine whether the actual moving-throat branch produces the required outgoing-load slippage relative to the wall baseline.

That is a much smaller and much cleaner theorem gate than we had before.
