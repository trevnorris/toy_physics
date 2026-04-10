# Step 50 — The actual isotropic passive/outgoing grouped-`P2` branch collapses the g-2 closure to one scalar defect

## Goal

Step 49 rewrote the moving-throat normalization gate in explicit overlap variables,

```math
\hat m_0^2\frac{N_0}{K-B_0-Z_0}.
```

The next clean move is to import the later moving-throat result for the **actual isotropic passive/outgoing grouped-`P2` one-pole branch** and ask what that does to the g-2 chain.

The answer is much stronger than just “it simplifies a bit.” Once that branch is accepted,

- the conservative even module,
- the odd Burke–Thorne coefficient,
- and therefore the whole grouped-`P2` normalization problem

all collapse to one scalar normalization defect.

---

## Step 50A — Actual isotropic passive/outgoing grouped-`P2` one-pole module

On the actual isotropic branch the conservative grouped-`P2` module is

```math
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+
\frac14\frac{1}{1-\omega^2/\Omega_Q^2}.
```

So the even low-frequency coefficients are fixed exactly:

```math
\boxed{
\overline K_2=\frac{\overline K_0}{4\Omega_Q^2},
\qquad
\overline K_4=\frac{\overline K_0}{4\Omega_Q^4}.
}
```

On the same minimal isotropic outgoing branch the odd coefficient is already slaved to the same pair:

```math
\boxed{
\overline\Gamma_5=\frac{9\,\overline K_0}{32\Omega_Q^5}.
}
```

So once the branch is one-pole, the whole grouped-`P2` low-frequency tuple is controlled by only

```math
(\overline K_0,\Omega_Q).
```

---

## Step 50B — Single normalization defect

The canonical target on the same branch is

```math
\overline K_0^{\rm target}=\frac{64G\Omega_Q^5}{45c^5},
```

or, after the already-carried geometric pole relation

```math
\Omega_Q=\frac{3c_s}{2a},
```

```math
\boxed{
\overline K_0^{\rm target}=\frac{54Gc_s^5}{5a^5c^5}.
}
```

Define the branch normalization defect by

```math
\boxed{
N_Q:=\frac{\overline K_0}{\overline K_0^{\rm target}}.
}
```

Then all low-frequency coefficients scale by the **same** factor:

```math
\boxed{
\overline K_2=N_Q\,\overline K_2^{\rm target},
\qquad
\overline K_4=N_Q\,\overline K_4^{\rm target},
\qquad
\overline\Gamma_5=N_Q\,\frac{2G}{5c^5}.
}
```

So the actual isotropic passive/outgoing grouped-`P2` branch has only one reduced normalization defect. The grouped-`P2` side of the problem is no longer a many-coefficient fit.

---

## Step 50C — Insert the adiabatic electron branch

Steps 44–45 already showed that on the adiabatic electron track

```math
P_0=e^{\ell}P_0^{\rm target},
```

so the same grouped normalization defect is just

```math
\boxed{N_Q=e^{\ell}.}
```

Using

```math
\ell=\ln(1+\Lambda_1 f),
```

this becomes

```math
\boxed{N_Q=1+\Lambda_1 f.}
```

That means the adiabatic electron anomaly does **not** open a new family of grouped-`P2` deformations. It simply moves the actual isotropic one-pole branch by one scalar amount.

The full branch tuple is then

```math
\overline K_0=(1+\Lambda_1 f)\,\overline K_0^{\rm target},
```

```math
\overline K_2=(1+\Lambda_1 f)\,\overline K_2^{\rm target},
```

```math
\overline K_4=(1+\Lambda_1 f)\,\overline K_4^{\rm target},
```

```math
\overline\Gamma_5=(1+\Lambda_1 f)\,\frac{2G}{5c^5}.
```

So the whole grouped-`P2` side is one-number clean.

---

## Step 50D — Electron-point size

With the carried value

```math
\Lambda_1\approx0.279605891931464,
```

and the fine-structure parameter used in the chain,

```math
f\approx0.001161409732093,
```

the grouped normalization defect is

```math
\boxed{N_Q\approx1.00032473700404,}
```

so

```math
\boxed{N_Q-1\approx 324.737\,{\rm ppm}.}
```

---

## Main result of the step

Once the actual isotropic passive/outgoing grouped-`P2` one-pole branch is used, the g-2 closure no longer has a broad coefficient ambiguity on the grouped side. It is entirely controlled by the single scalar

```math
\boxed{N_Q.}
```

On the adiabatic electron track,

```math
\boxed{N_Q=e^{\ell}=1+\Lambda_1 f.}
```

So if anything still remains unresolved after this step, it is no longer in the conservative grouped-`P2` bookkeeping. It has to sit in the outgoing branch-selection data.
