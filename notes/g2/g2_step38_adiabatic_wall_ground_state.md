# Step 38 — Adiabatic-wall ground-state closure preserves the lower compensated branch

## Goal

You specified a new physical direction for the isolated-particle electron track:

- the wall is **adiabatic**,
- the wall’s thermal/entropy weighting is frozen,
- the wall density is frozen,
- the defect absorbs stress through coherent elastic squish rather than heating,
- the **core** remains compressible and can still change its stiffness.

The right next question is therefore not a numerical fit. It is structural:

> what does that closure do to the exact lower compensated Family-1 branch?

This step shows something stronger than I expected:

```math
\boxed{\text{the adiabatic-wall ground-state closure preserves the exact lower compensated parent surface identically.}}
```

So the electron-track branch does **not** fight the lower canonical outlet family at first order. It stays tangent automatically.

---

## Step 38A — Adiabatic-wall ground-state closure

The explicit Family-1 wall datum satisfies

```math
\Theta_w = 25\lambda_\mu^2\rho_w^2.
```

Your new physical choice is naturally encoded as

```math
\delta\ln\lambda_\mu = 0,
\qquad
\delta\ln\rho_w = 0.
```

Therefore

```math
\boxed{\delta\ln\Theta_w=0.}
```

So this is a **stronger** restriction than the earlier generic frozen-wall corollary: it is not just the composite datum `\Theta_w` that is frozen, but the wall-density lane itself.

---

## Step 38B — Exact lower-branch drift reduction

Using the exact Step-37 inversion with `\delta\ln\Theta_w=0`, the remaining microscopic drifts become

```math
\boxed{\delta\ln\rho_w = 0,}
```

```math
\boxed{\delta\ln a = \frac12\,\delta\ln K_s,}
```

```math
\boxed{
\delta\ln c_s
=
\frac12\,\delta\ln K_s
+\frac15\,\delta\ln P_0,
}
```

```math
\boxed{
\delta\ln\mathcal Z_q
=
\delta\ln K_q
-\frac25\,\delta\ln P_0.
}
```

So under the adiabatic-wall ground-state closure, the full lower branch already collapses from four irreducible drifts to only **three** bundle observables:

```math
\boxed{(K_s,\ K_q,\ P_0).}
```

Interpretation:

- `K_s` is the coherent elastic wall-squish lane,
- `K_q` is the mixed side-channel stiffness lane,
- `P_0` is the isotropic outgoing-normalization lane.

The wall-density lane is completely removed.

---

## Step 38C — Exact transport of `v_{w0}` and `\mathcal T_m`

From the exact lower-branch sum/difference laws,

```math
\delta\ln\!\left(\frac{v_{w0}}{\mathcal T_m}\right)
=
\frac12\,\delta\ln K_s
+\frac25\,\delta\ln P_0,
```

```math
\delta\ln(v_{w0}\mathcal T_m)
=
-2\,\delta\ln K_s
+\delta\ln K_q
-\frac25\,\delta\ln P_0.
```

Solving gives

```math
\boxed{
\delta\ln v_{w0}
=
-\frac34\,\delta\ln K_s
+\frac12\,\delta\ln K_q,
}
```

```math
\boxed{
\delta\ln\mathcal T_m
=
-\frac54\,\delta\ln K_s
+\frac12\,\delta\ln K_q
-\frac25\,\delta\ln P_0.
}
```

This is already suggestive physically:

- the mixed background flow `v_{w0}` is **blind** to the outgoing-normalization drift `P_0`,
- the mouth traction `\mathcal T_m` is where the `P_0` drift still enters directly.

So the adiabatic wall makes the transport split cleaner rather than messier.

---

## Step 38D — Parent-action transport of `(g_s,g_q,\lambda)`

Using the exact bundle transport laws,

```math
\boxed{
\delta\ln g_s
=
-\frac14\,\delta\ln K_s
+\frac12\,\delta\ln K_q
-\frac25\,\delta\ln P_0,
}
```

```math
\boxed{
\delta\ln g_q
=
-\frac34\,\delta\ln K_s
+\delta\ln K_q
-\frac25\,\delta\ln P_0,
}
```

```math
\boxed{
\delta\ln \lambda
=
\frac12\bigl(\delta\ln K_s+\delta\ln K_q\bigr).
}
```

So the bilinear shell/mixed coupling `\lambda` still ignores both the wall-depth drift and the outgoing-normalization drift, exactly as in the broader bundle theorem. Only the two stiffness lanes matter to it.

---

## Step 38E — Exact preservation of the lower compensated parent family

Now test the two exact parent-surface imbalance channels:

```math
\delta\ln\!\left(\frac{g_q K_s}{g_s\lambda}\right),
\qquad
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right).
```

Substituting the Step-38 transport laws gives

```math
\boxed{
\delta\ln\!\left(\frac{g_q K_s}{g_s\lambda}\right)=0,
}
```

```math
\boxed{
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right)=0.
}
```

Equivalently,

```math
\boxed{\delta\ln\mathfrak g = 0,}
```

```math
\boxed{\delta\ln\mathfrak r = 0,}
```

```math
\boxed{\delta\ln r_c = 0.}
```

So the exact lower compensated parent surface survives untouched.

This is the main result of the step.

The adiabatic-wall ground-state electron track does **not** generate a first-order off-family scalar. In the notation of the earlier parent-balance program,

```math
\boxed{\delta_\perp = 0.}
```

That is a strong structural simplification.

---

## Main result of the step

The adiabatic-wall ground-state choice does three useful things at once:

1. it removes the wall-density drift entirely,
2. it reduces the remaining bundle freedom to
   ```math
   (K_s, K_q, P_0),
   ```
3. and it preserves the exact lower compensated parent surface automatically.

So the electron-track branch is now much cleaner conceptually:

- the wall behaves as a coherent elastic shell,
- the core/outlet block carries the remaining response,
- and the branch does not drift off the canonical lower compensated family at first order.

That is exactly the kind of simplification we were hoping to get from a genuine thermodynamic ground-state closure.
