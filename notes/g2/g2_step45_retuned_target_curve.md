# Step 45 — On the adiabatic anomaly track the moving-throat normalization law survives as a retuned target curve

## Goal

Step 44 showed that the pure-odd DtN deformation from Step 43 does not stay only on
the outlet side. On the natural point-particle source-map branch it forces the
grouped-`P2` prefactor to rescale as

```math
P_0 = e^\ell P_0^{\rm target},
\qquad
P_0^{\rm target}=\frac{54Gc_s^5}{5a^5c^5}.
```

But Steps 40–41 had already shown that the same adiabatic anomaly-only branch retunes
the core sound scale by

```math
\delta\ln a = 0,
\qquad
\delta\ln c_s = \frac{\ell}{5}.
```

So the next honest question is:

> is Step 44 telling us that the anomaly moves **off** the universal moving-throat
> normalization law, or does it merely move the branch **along the same law** by
> retuning `c_s`?

This step shows the stronger and cleaner result:

```math
\boxed{
P_0(a,c_s e^{\ell/5}) = e^\ell P_0(a,c_s).
}
```

So the adiabatic anomaly rides the **same universal target curve**. It does not
break it.

---

## Step 45A — Exact retuned-target law

Take the universal moving-throat grouped-`P2` normalization curve

```math
\boxed{
P_0^{\rm target}(a,c_s)=\frac{54Gc_s^5}{5a^5c^5}.
}
```

Now insert the anomaly-only adiabatic retuning from Step 41,

```math
a \mapsto a,
\qquad
c_s \mapsto c_s e^{\ell/5}.
```

Then

```math
P_0^{\rm target}(a,c_s e^{\ell/5})
=
\frac{54G(c_s e^{\ell/5})^5}{5a^5c^5}
=
e^\ell P_0^{\rm target}(a,c_s).
```

So the Step-44 prefactor enhancement is **exactly** the same as evaluating the
same universal target formula at the retuned adiabatic branch values.

That is the cleanest physical reading so far:

> the quartic anomaly does not ask for a new normalization law; it asks for a very
> small retuning of the core sound scale along the same normalization law.

---

## Step 45B — Minimal constant-prefactor grouped-`P2` continuation

The pure-odd DtN representative from Step 43 does **not** by itself force a unique
choice of the higher grouped-`P2` prefactor moments `(P_2,P_4)`. But the moving-throat
notes already isolate one especially simple continuation:

```math
\boxed{P_2 = 0,\qquad P_4 = 0.}
```

On that minimal isotropic branch,

```math
\boxed{
K_2 = \frac{P_0 a^2}{9c_s^2},
\qquad
K_4 = \frac{4P_0 a^4}{81c_s^4}.
}
```

Substituting the adiabatic retuning `c_s\mapsto c_s e^{\ell/5}` with fixed `a`
gives the finite scaling laws

```math
\boxed{
P_0 \mapsto e^\ell P_0,
\qquad
K_2 \mapsto e^{3\ell/5} K_2,
\qquad
K_4 \mapsto e^{\ell/5} K_4.
}
```

So on the minimal grouped-`P2` continuation the quartic anomaly induces:

- a `324.7` ppm upward shift in the static prefactor,
- a `194.8` ppm upward shift in the quadratic even coefficient,
- and a `64.94` ppm upward shift in the quartic even coefficient.

At tangent level this is just

```math
\boxed{
\delta\ln P_0 = \ell,
\qquad
\delta\ln K_2 = \frac{3}{5}\ell,
\qquad
\delta\ln K_4 = \frac{1}{5}\ell.
}
```

---

## Step 45C — Exact constant-prefactor bundle conditions

The moving-throat grouped-`P2` bundle also gives exact formulas for the prefactor
moments in terms of the conservative bundle data:

```math
P_0 = \frac{N_0}{D_0},
```

```math
P_2 = \frac{D_0N_2 - 2D_2N_0}{D_0^2},
```

```math
P_4 =
\frac{D_0^2N_4 - 2D_0(D_2N_2+D_4N_0)+3D_2^2N_0}{D_0^3}.
```

The conservative normalized response moments are

```math
u_2 = -\frac{D_2}{D_0},
\qquad
u_4 = \frac{D_2^2 - D_0 D_4}{D_0^2}.
```

If we now impose the minimal constant-prefactor branch

```math
P_2=P_4=0,
```

the outgoing transfer moments are **not** free. They are slaved to the conservative
grouped response:

```math
\boxed{
N_2 = -2u_2 N_0,
}
```

```math
\boxed{
N_4 = (3u_2^2 - 2u_4)\,N_0.
}
```

So if the true moving-throat branch really selects this minimal continuation, the
PDE no longer gets to choose `N_2` and `N_4` independently once the conservative
bundle moments `(u_2,u_4)` and the static transfer weight `N_0` are known.

That is a sharp next theorem target.

---

## Step 45D — What is established, and what is still conditional

### Established now

The adiabatic anomaly branch is fully compatible with the universal moving-throat
normalization curve:

```math
\boxed{
P_0^{\rm actual}
=
\frac{54G c_s^5}{5a^5c^5}
}
```

provided `c_s` and `a` are understood as the **retuned adiabatic branch values**,
with

```math
\delta\ln a=0,
\qquad
\delta\ln c_s=\ell/5.
```

So the anomaly does **not** ask for a new normalization law.

### Still conditional

The extra step

```math
P_2=P_4=0
```

is the **minimal grouped-`P2` completion**, not yet a theorem forced by Step 43
alone. More complicated isotropic branches with `P_2,P_4\neq0` are still logically
available until the full moving-throat PDE fixes them.

So the honest status is:

- the target curve is now exact on the adiabatic anomaly branch,
- the constant-prefactor continuation is the cleanest next completion,
- but the PDE still has to decide whether that minimal branch is the actual one.

---

## Main result of the step

The adiabatic quartic anomaly does **not** push the moving-throat quadrupole
normalization off the universal target curve. It simply retunes the core sound
scale so that

```math
\boxed{
P_0(a,c_s e^{\ell/5}) = e^\ell P_0(a,c_s).
}
```

If we also choose the minimal constant-prefactor isotropic grouped-`P2` continuation,

```math
\boxed{
P_2=P_4=0,
}
```

then

```math
\boxed{
P_0 \mapsto e^\ell P_0,
\qquad
K_2 \mapsto e^{3\ell/5}K_2,
\qquad
K_4 \mapsto e^{\ell/5}K_4,
}
```

and the outgoing transfer moments are forced to obey

```math
\boxed{
N_2 = -2u_2 N_0,
\qquad
N_4 = (3u_2^2 - 2u_4)N_0.
}
```

So the next moving-throat question is no longer about the existence of a
normalization curve. It is whether the true PDE selects this minimal constant-
prefactor isotropic completion of the already-fixed adiabatic anomaly track.
