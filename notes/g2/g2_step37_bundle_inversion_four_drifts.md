# Step 37 — Exact bundle inversion of the last four irreducible branch drifts

## Goal

Step 36 reduced the exact lower compensated Family-1 branch to only four truly
irreducible microscopic drifts,

```math
\delta\ln\mathcal Z_q,
\qquad
\delta\ln\rho_w,
\qquad
\delta\ln c_s,
\qquad
\delta\ln a.
```

The right next move is not to guess those drifts one by one. The grouped
wall/BdG/Maxwell bundle already exposes four natural observables that determine
all of them algebraically:

```math
\Theta_w,
\qquad
K_s,
\qquad
K_q,
\qquad
P_0=\frac{N_0}{D_0}.
```

This step performs that inversion exactly.

---

## Step 37A — The four exact branch laws

On the explicit Family-1 wall branch,

```math
\Theta_w=25\lambda_\mu^2\rho_w^2,
```

so at fixed `\lambda_\mu`,

```math
\boxed{
\delta\ln\Theta_w = 2\,\delta\ln\rho_w.
}
```

On the healing-locked shell branch,

```math
K_s\propto a^2\rho_w,
```

so

```math
\boxed{
\delta\ln K_s = 2\,\delta\ln a + \delta\ln\rho_w.
}
```

On the exact lower compensated D/N branch,

```math
K_q\propto \mathcal Z_q\,\frac{c_s^2}{L_W^2},
qquad
\delta\ln L_W=\delta\ln a,
```

so

```math
\boxed{
\delta\ln K_q
=
\delta\ln\mathcal Z_q + 2\,\delta\ln c_s - 2\,\delta\ln a.
}
```

And on the isotropic outgoing quadrupole branch,

```math
P_0\propto \frac{c_s^5}{a^5},
```

so

```math
\boxed{
\delta\ln P_0 = 5\bigl(\delta\ln c_s - \delta\ln a\bigr).
}
```

These four equations are the entire inversion system.

---

## Step 37B — Exact inversion

Solving the system gives

```math
\boxed{
\delta\ln\rho_w
=
\frac12\,\delta\ln\Theta_w,
}
```

```math
\boxed{
\delta\ln a
=
\frac12\,\delta\ln K_s
-
\frac14\,\delta\ln\Theta_w,
}
```

```math
\boxed{
\delta\ln c_s
=
\frac12\,\delta\ln K_s
-
\frac14\,\delta\ln\Theta_w
+
\frac15\,\delta\ln P_0,
}
```

```math
\boxed{
\delta\ln\mathcal Z_q
=
\delta\ln K_q
-
\frac25\,\delta\ln P_0.
}
```

So the last four branch drifts are not open in any diffuse sense anymore. They are
exact algebraic images of the bundle observables `\Theta_w, K_s, K_q, P_0`.

---

## Step 37C — Full-bundle form using `P_0=N_0/D_0`

Because the grouped-bundle normalization is

```math
P_0=\frac{N_0}{D_0},
qquad
\delta\ln P_0 = \delta\ln N_0 - \delta\ln D_0,
```

we may rewrite the last two drifts directly in full-bundle language:

```math
\boxed{
\delta\ln c_s
=
\frac12\,\delta\ln K_s
-
\frac14\,\delta\ln\Theta_w
+
\frac15\bigl(\delta\ln N_0-\delta\ln D_0\bigr),
}
```

```math
\boxed{
\delta\ln\mathcal Z_q
=
\delta\ln K_q
-
\frac25\bigl(\delta\ln N_0-\delta\ln D_0\bigr).
}
```

So the full grouped wall/BdG/Maxwell bundle enters only through the isotropic
static normalization ratio `N_0/D_0`, exactly as the 2.5PN normalization package
had been suggesting all along.

---

## Step 37D — Frozen-wall corollary

If the explicit Family-1 wall datum is held fixed at first order,

```math
\delta\ln\Theta_w=0,
```

then the inversion simplifies to

```math
\boxed{
\delta\ln\rho_w = 0,
}
```

```math
\boxed{
\delta\ln a = \frac12\,\delta\ln K_s,
}
```

```math
\boxed{
\delta\ln c_s = \frac12\,\delta\ln K_s + \frac15\,\delta\ln P_0,
}
```

```math
\boxed{
\delta\ln\mathcal Z_q = \delta\ln K_q - \frac25\,\delta\ln P_0.
}
```

So even this explicit wall-preserving restriction already removes one of the four
would-be free drifts.

---

## Main result of the step

The lower-branch problem is now reduced as far as it can go without solving the
remaining moving-throat PDE.

The four residual microscopic drifts are exact algebraic images of the bundle
observables:

```math
\boxed{
\delta\ln\rho_w
=
\frac12\,\delta\ln\Theta_w,
qquad
\delta\ln a
=
\frac12\,\delta\ln K_s-rac14\,\delta\ln\Theta_w,
}
```

```math
\boxed{
\delta\ln c_s
=
\frac12\,\delta\ln K_s-rac14\,\delta\ln\Theta_w+rac15\,\delta\ln P_0,

\qquad
\delta\ln\mathcal Z_q
=
\delta\ln K_q-rac25\,\delta\ln P_0.
}
```

And because

```math
P_0=\frac{N_0}{D_0},
```

the only grouped-bundle quantity that still matters is the isotropic static
normalization ratio `N_0/D_0`.

So the next PDE-facing question is no longer “what are the four drifts?” It is:

> what are the first-order drifts of the four bundle observables `\Theta_w, K_s, K_q, P_0` on the actual isotropic moving-throat branch?
