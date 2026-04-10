# Step 1 — Exact moving-throat quotient bridge and staggered anomaly benchmark

## Goal

Start the rebuild from the **new moving-throat PDE language**, not from the old sequential collar story.
The purpose of this first step is:

1. to replace the old “charge update, then inertia update, then charge update again” bookkeeping by the exact moving-throat quotient coordinates, and
2. to rebuild the **current** staggered anomaly law exactly so that the later common charge–inertia correction is compared against the right baseline.

This is the correct beginning because the old write-up already says the remaining miss is the first genuinely **common** charge–inertia layer at `O(f^4)`, while the new PDE handoff compresses the microscopic coupled branch into three exact quotient coordinates.

---

## Inputs carried forward

### From the anomaly write-up

The current exact reduced closure has

- `f = alpha_fs / (2 pi)`,
- `g_loc ≈ 2.00231930435865`,
- a remaining miss of about `2.27e-12`,
- and the explicit warning that the unresolved piece is the first **common charge–inertia transport layer** at `O(f^4)`.

The anomaly text is also explicit that the previous successful layers were still partially **staggered**:

1. charge support moved into the collar,
2. inertia backreacted through the `P_22` self-loop and blur,
3. charge then acquired its own local collar mode.

So the next correction is expected to come from a genuinely coupled local law rather than one more one-sided update.

### From the moving-throat PDE handoff

The new exact microscopic state vector is

```math
x=(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})},\ K_W^{(\mathrm{eff})},\ \mu_W,\ T_U),
```

with grouped weak-axisymmetric log-drift vector

```math
\delta\mathbf x=
\begin{pmatrix}
\lambda_1\\ c_1\\ \gamma_1\\ \kappa_U\\ \kappa_\eta\\ \kappa_W\\ \mu_1\\ \tau_1
\end{pmatrix}.
```

The three exact quotient coordinates are the monomials

```math
\mathfrak C_{{\rm tr},*},\qquad
\mathfrak C_{{\rm nt},*},\qquad
\epsilon_\eta,
```

with exact invariant map

```math
\mathcal I(x)=\bigl(\mathfrak C_{{\rm tr},*}(x),\ \mathfrak C_{{\rm nt},*}(x),\ \epsilon_\eta(x)\bigr),
```

and exact log-drift law

```math
\begin{pmatrix}
\delta\ln \mathfrak C_{{\rm tr},*}\\
\delta\ln \mathfrak C_{{\rm nt},*}\\
\delta\ln \epsilon_\eta
\end{pmatrix}
=
M_*\,\delta\mathbf x.
```

The handoff further states that `rank(M_*) = 3`, so the true defect motion lives in a **3-dimensional quotient**, while the remaining 5 directions are pure similarity-orbit motion.

---

## Step 1A — Exact quotient coordinates and similarity-orbit decomposition

I define the quotient drift vector

```math
\mathbf q=
\begin{pmatrix}
q_{\rm tr}\\ q_{\rm nt}\\ q_\eta
\end{pmatrix}
:=
\begin{pmatrix}
\delta\ln \mathfrak C_{{\rm tr},*}\\
\delta\ln \mathfrak C_{{\rm nt},*}\\
\delta\ln \epsilon_\eta
\end{pmatrix}
=M_*\,\delta\mathbf x.
```

The SymPy audit verifies all of the following exactly.

### 1. `M_*` has rank 3 and nullity 5

So any microscopic grouped weak-axisymmetric drift splits into:

- **5 similarity directions** that preserve the exact monomials, and
- **3 quotient directions** that move the actual branch in
  `((mathfrak C_tr,*),(mathfrak C_nt,*),epsilon_eta)`.

### 2. There is an explicit right inverse `R`

Choosing the quotient lift in which the five free co-scalings

```math
(\lambda_1,c_1,\gamma_1,\kappa_U,\kappa_W)
```

are set to zero gives a concrete right inverse `R` with

```math
M_*R=I_3.
```

So the quotient motion is not abstract.  It can be injected back into the microscopic drift space exactly.

### 3. The dependent microscopic drifts solve exactly as

```math
\kappa_\eta
=
2c_1-\kappa_U-q_\eta,
```

```math
\tau_1
=
\kappa_U
-
\frac{1+\delta_{U,*}}{1+\chi_{0,*}}\,(\gamma_1+c_1-\kappa_U)
+
\frac{q_{\rm tr}}{1+\chi_{0,*}},
```

```math
\mu_1
=
2c_1-\kappa_U+2\kappa_W-2\lambda_1
-
E_*\,(2\gamma_1+2\lambda_1-\kappa_U-\kappa_W)
-
F_*\,\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\gamma_1+c_1-\kappa_U)
+
q_{\rm nt}-q_\eta+
\frac{F_*}{1+\chi_{0,*}}q_{\rm tr}.
```

Setting `q = 0` recovers the exact tangent **similarity orbit**.  So the new PDE picture makes the conceptual split precise:

- the five free co-scalings are not physical defect motion,
- the physical defect motion is exactly the 3-vector `q`.

### 4. The immediate consequence for the anomaly program

Any genuine **common** charge–inertia correction should be organized as a scalar built from

```math
q_{\rm tr},\qquad q_{\rm nt},\qquad q_\eta,
```

not from arbitrary sequential changes of charge and inertia support separately.

This is the key structural result of Step 1.

---

## Step 1B — Rebuild the current staggered anomaly law exactly

The common layer has to be added to the **current** local closure, not to a looser schematic formula.
So I also reconstructed the present staggered anomaly law directly from the Appendix F/G definitions.

### 1. Exact collar-mode integrals

Using the collar coordinate

```math
\tau(f)=1-\sqrt{1-f},
```

and the rescaled collar variable `s in [0,1]`, the charge-side mean subtraction and second moment are exact closed forms.

The exact mean subtraction is

```math
\bar c(\tau)
=
\frac{4\bigl(-2\tau+\pi(\tau-1)\bigr)}{\pi^2(\tau-2)}.
```

The exact second moment is

```math
\Xi(\tau)
=
\frac{2\tau\Bigl(-48\pi\tau^2+2\pi^2\tau^2+\pi^3\tau^2+96\tau^2
-3\pi^3\tau-4\pi^2\tau+48\pi\tau-8\pi^2+2\pi^3\Bigr)}{\pi^4(\tau-2)}.
```

This is useful because it means the current charge-side local mode can be benchmarked exactly, not only to cubic order.

### 2. Exact `Q_loc(f)` series through quartic order

The reconstructed charge second-moment ratio is

```math
Q_{\rm loc}(f)
=
1+f-f^2+
\frac{4-\pi}{\pi^2\kappa}f^3
+
\frac{4(\pi-3)}{\pi^3\kappa}f^4
+
O(f^5).
```

So the current **charge-side** quartic contribution is already fixed by the existing staggered closure.

### 3. Exact current staggered `g_loc(f)` series through quartic order

Combining the charge-side series with the exponential-blur inertia factor gives

```math
g_{\rm loc}(f)/2
=
1+f-
\frac{47}{36}f^2
+
\left(
\frac{11}{6}\kappa+
\frac{4-\pi}{\pi^2\kappa}
\right)f^3
+
a_{4,\rm staggered}f^4
+
O(f^5),
```

with

```math
a_{4,\rm staggered}
=
-
\frac{55}{6}\kappa^2
+
\frac{4(\pi-3)}{\pi^3\kappa}.
```

For the frozen benchmark value

```math
\kappa=1.177746578880,
```

this evaluates to

```math
a_{4,\rm staggered}\approx -12.6994546522869.
```

This is important:

> the current exact staggered law already carries a definite quartic coefficient.

So the unresolved common layer should **not** be identified with “the quartic coefficient of the full exact current local law.”
It is an **incremental** coupled correction on top of that law.

### 4. Exact physical benchmark with the current numbers

Using the frozen benchmark inputs, the reconstructed exact local closure gives

```math
g_{\rm loc}
=
2.002319304358647956\ldots
```

against the adopted target

```math
g_e = 2.00231930436092.
```

So the exact residual is

```math
\Delta g
=
g_e-g_{\rm loc}
=
2.27204390584705\times 10^{-12}.
```

That agrees with the `2.27e-12` scale reported in the write-up, but now we have the exact number used by the script.

---

## Step 1C — The quartic sign issue

If one compares the measured residual only against the **cubic truncation**

```math
g/2 = 1+f-c_2f^2+c_3f^3,
```

then the raw added quartic coefficient would be

```math
a_{4,\rm bench}
=
\frac{\Delta g}{2f^4}
\approx 0.624374101073809.
```

But Part VIII writes the next term as

```math
-c_4 f^4.
```

So in **that** sign convention,

```math
c_4 = -a_{4,\rm bench}
\approx -0.624374101073809.
```

That sign matters.
The positive number is the raw coefficient of `+f^4`; the Part VIII `c_4` is the negative of that.

---

## What Step 1 establishes

### 1. The new PDE base naturally replaces staggered bookkeeping by exact quotient motion

The moving-throat handoff gives a clean split:

- `5` similarity directions = zero-cost co-scalings,
- `3` quotient directions = true defect motion.

So the common charge–inertia layer should be written as a coupled evolution for the quotient vector

```math
\mathbf q=(q_{\rm tr},q_{\rm nt},q_\eta),
```

not as one more sequential update of the old charge and inertia ledgers.

### 2. The common layer must be benchmarked against the exact current local law

The exact present staggered closure already contains a real quartic series coefficient,
so the missing common term is **incremental**:

```math
\Delta g_{\rm common}(f)=O(f^4)
```

on top of the exact current `g_loc(f)`.

### 3. The Step-2 problem is now sharp

The next step should be:

1. choose a one-parameter coupled transport path
   ```math
   \mathbf q(f),
   ```
2. impose minimal conditions so that the new common correction begins at `O(f^4)`,
3. expand the lowest-order scalar correction built from
   `q_tr`, `q_nt`, and `q_eta`,
4. and match that correction against the exact residual **after** subtracting the current staggered baseline.

That is the cleanest beginning I know how to make from the new PDE paper.
