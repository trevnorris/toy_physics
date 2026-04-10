# Step 27 — The actual isotropic passive/outgoing branch has one normalization defect

## Goal

Step 26 showed that on the actual isotropic branch the geometry lane is dynamically
inert through `O(ω^4)`, so the conservative grouped-`P2` carrier is exactly the
forced

```math
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2}.
```

The next honest move is to combine that conservative branch with the already-carried
minimal outgoing grouped-`P2` law and ask:

> **how many independent normalization defects are left on the actual isotropic
> passive/outgoing branch?**

This step shows the answer is **one**.

The main result is that, once the actual branch is isotropic and one-pole,

```math
\boxed{N_Q := \bar K_0/\bar K_0^{\rm target}}
```

controls **all** low-frequency invariant mismatches:

```math
\boxed{
\frac{\bar K_2}{\bar K_2^{\rm target}}
=
\frac{\bar K_4}{\bar K_4^{\rm target}}
=
\frac{\bar\Gamma_5}{\bar\Gamma_5^{\rm target}}
=
N_Q.
}
```

So the full reduced 2.5PN/4PN theorem gap on that branch is equivalent to the one
scalar condition

```math
\boxed{N_Q=1.}
```

---

## Step 27A — Exact low-frequency coefficients of the actual conservative module

Start from the actual isotropic grouped-`P2` conservative module

```math
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2}.
```

Expanding at low frequency gives

```math
\widehat Y_Q^{\rm cons}(\omega)
=
1+
\frac{\omega^2}{4\Omega_Q^2}
+
\frac{\omega^4}{4\Omega_Q^4}
+O(\omega^6).
```

So the exact conservative moments are

```math
\boxed{\bar K_2 = \bar K_0/(4\Omega_Q^2),}
```

```math
\boxed{\bar K_4 = \bar K_0/(4\Omega_Q^4).}
```

There is no extra even ambiguity once `\bar K_0` and `\Omega_Q` are fixed.

---

## Step 27B — The odd Burke–Thorne coefficient is not independent either

On the same minimal isotropic outgoing branch, the 2.5PN audit already fixed the
odd coefficient algebraically as

```math
\bar\Gamma_5
=
9\,\bar K_2^{5/2}/\bar K_0^{3/2}.
```

Substituting the conservative one-pole relation for `\bar K_2` gives

```math
\boxed{\bar\Gamma_5 = \frac{9\bar K_0}{32\Omega_Q^5}.}
```

So the odd Burke–Thorne coefficient is tied to the same two quantities
`(\bar K_0,\Omega_Q)` and is not a separate datum.

---

## Step 27C — Exact target branch

On the GR target branch,

```math
\bar K_0^{\rm target} = \frac{64G\Omega_Q^5}{45c^5}.
```

Then automatically

```math
\bar K_2^{\rm target} = \bar K_0^{\rm target}/(4\Omega_Q^2),
```

```math
\bar K_4^{\rm target} = \bar K_0^{\rm target}/(4\Omega_Q^4),
```

```math
\bar\Gamma_5^{\rm target} = \frac{2G}{5c^5}.
```

Using the already-carried geometric pole

```math
\Omega_Q = \frac{3c_s}{2a},
```

this becomes

```math
\boxed{
\bar K_0^{\rm target}
=
\frac{54Gc_s^5}{5a^5c^5}.
}
```

---

## Step 27D — One scalar defect controls the whole branch

Define the actual-branch normalization defect by

```math
\boxed{N_Q := \bar K_0/\bar K_0^{\rm target}.}
```

Then the actual branch can be written as

```math
\bar K_0 = N_Q\,\bar K_0^{\rm target},
```

and therefore

```math
\bar K_2 = N_Q\,\bar K_2^{\rm target},
```

```math
\bar K_4 = N_Q\,\bar K_4^{\rm target},
```

```math
\bar\Gamma_5 = N_Q\,\bar\Gamma_5^{\rm target} = N_Q\,\frac{2G}{5c^5}.
```

So the four natural defect measures

```math
R_0 := \bar K_0/\bar K_0^{\rm target}-1,
```

```math
R_2 := \bar K_2/\bar K_2^{\rm target}-1,
```

```math
R_4 := \bar K_4/\bar K_4^{\rm target}-1,
```

```math
R_5 := \bar\Gamma_5/\bar\Gamma_5^{\rm target}-1,
```

all collapse to the same number:

```math
\boxed{R_0=R_2=R_4=R_5=N_Q-1.}
```

---

## Main result of the step

On the actual isotropic grouped-`P2` one-pole passive/outgoing branch,

```math
\boxed{N_Q := \bar K_0/\bar K_0^{\rm target}}
```

is the **only** reduced normalization defect.

Everything else follows from it:

```math
\boxed{
\frac{\bar K_2}{\bar K_2^{\rm target}}
=
\frac{\bar K_4}{\bar K_4^{\rm target}}
=
\frac{\bar\Gamma_5}{\bar\Gamma_5^{\rm target}}
=
N_Q.
}
```

So the full reduced 2.5PN/4PN theorem gap is now equivalent to

```math
\boxed{N_Q=1.}
```

That is the cleanest conservative normalization statement reached so far.
