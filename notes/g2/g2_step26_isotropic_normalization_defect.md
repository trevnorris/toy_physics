# Step 26 — The actual isotropic branch collapses to one normalization defect

## Goal

Step 25 showed that if the conservative higher-order branch is

```math
\text{one isotropic grouped-}P_2\text{ pole} + \text{static geometry completion},
```

then the conservative quadrupole module is forced to be

```math
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+
\frac14\frac{1}{1-\omega^2/\Omega_Q^2},
```

and the support/source ratio follows as

```math
\rho_\alpha=\frac43.
```

But one real reduced ambiguity was still visible:

> **could the geometry lane contaminate the grouped quadrupole module at
> `O(\omega^2)` or `O(\omega^4)` and therefore spoil the forced `3/4 + 1/4`
> split?**

This step answers that and then pushes the whole branch to its sharpest reduced
form.

The main results are:

1. on the isotropic quadratic wall operator, the `l=0` geometry lane and the
   grouped real `l=2` bundle are exactly orthogonal,
2. therefore the dynamic geometry contamination numbers vanish,
   ```math
   \boxed{\epsilon_2 = \epsilon_4 = 0,}
   ```
3. so the actual isotropic branch really does realize the Step-25
   `3/4 + 1/4` conservative module,
4. and once that is true, the whole remaining reduced 2.5PN/4PN mismatch is just
   one scalar normalization defect
   ```math
   \boxed{N_Q := \bar K_0 / \bar K_0^{\rm target}.}
   ```

At the same time, the support/source side becomes automatic on that branch.

So this is the step where the old selected twin-support phase picture drops out
of the active theorem ledger.

---

## Step 26A — Exact isotropic geometry-decoupling theorem

Take the isotropic quadratic wall operator on the throat sphere. Because it is
`O(3)` invariant, its angular part depends only on the sphere Laplacian:

```math
L_{\rm ang} = a + b(-\Delta_{S^2}).
```

Expand the wall field into

- the scalar geometry lane `Y_{00}`,
- the grouped real quadrupole bundle `Y_{2A}` with
  ```math
  A\in\{20,21c,21s,22c,22s\}.
  ```

Using normalized real harmonics, the exact orthogonality relations are

```math
\int Y_{00} Y_{2A}\, d\Omega = 0.
```

And because

```math
-\Delta_{S^2}Y_{2A}=6Y_{2A},
```

we also have

```math
\langle Y_{00}, (a+b(-\Delta))Y_{2A} \rangle
=
(a+6b)\langle Y_{00},Y_{2A}\rangle
=0.
```

So on the isotropic branch the scalar/geometry `l=0` lane and the grouped real
`l=2` bundle are exactly block diagonal at linear order.

That is the key selection rule.

It means the scalar geometry lane cannot feed dynamic `\omega^2` or `\omega^4`
moments into the grouped quadrupole carrier on the actual isotropic branch.

---

## Step 26B — Dynamic geometry contamination vanishes

The obstruction formula for a geometry lane that carries its own dynamic even
moments is

```math
c_{\rm pole}
=
\frac{1+\epsilon_4}{4(1+\epsilon_2)^2},
```

with

```math
\epsilon_2 = \Omega_Q^2 K_{(g,2)}/K_{\rm pole},
\qquad
\epsilon_4 = \Omega_Q^4 K_{(g,4)}/K_{\rm pole}.
```

But the isotropic decoupling theorem above implies

```math
K_{(g,2)} = K_{(g,4)} = 0,
```

so

```math
\boxed{\epsilon_2 = \epsilon_4 = 0.}
```

Then the obstruction formula collapses to

```math
c_{\rm pole} = \frac14,
\qquad
c_{\rm geom} = 1-c_{\rm pole} = \frac34.
```

So the Step-25 contact/pole split is not merely a candidate branch value. It is
actually the isotropic reduced branch value once the geometry lane is checked.

---

## Step 26C — The full reduced mismatch is one scalar defect

Now write the actual isotropic passive/outgoing one-pole branch in invariant form
as

```math
\bar K_Q^{\rm cons}(\omega)
=
\bar K_0\left[
\frac34+
\frac14\frac{1}{1-\omega^2/\Omega_Q^2}
\right].
```

Define the GR target normalization

```math
\bar K_0^{\rm target} = \frac{64 G\Omega_Q^5}{45 c^5}.
```

Using the already-carried geometric pole

```math
\Omega_Q = \frac{3c_s}{2a},
```

this is equivalently

```math
\bar K_0^{\rm target}
=
\frac{54Gc_s^5}{5a^5c^5}.
```

Now define the actual branch normalization defect

```math
\boxed{N_Q := \frac{\bar K_0}{\bar K_0^{\rm target}}.}
```

Then every low-frequency invariant scales by the same factor:

```math
\bar K_2 = N_Q\,\bar K_2^{\rm target},
```

```math
\bar K_4 = N_Q\,\bar K_4^{\rm target},
```

```math
\bar\Gamma_5 = N_Q\,\frac{2G}{5c^5}.
```

So the full reduced GR-like point-particle 2.5PN closure on the actual isotropic
branch is now equivalent to the one equation

```math
\boxed{N_Q = 1.}
```

That is the narrowest reduced theorem gate reached so far.

---

## Step 26D — The support/source side becomes automatic

From Step 25 the actual isotropic branch has

```math
\rho_\alpha = \frac43.
```

The exact blocked support demand therefore becomes

```math
\zeta_{\rm req}^{\rm(act)}(\epsilon_{\rm blk})
=
\frac{\rho_\alpha-1}{1-\epsilon_{\rm blk}(2-\rho_\alpha)}
=
\frac{1}{3-2\epsilon_{\rm blk}}.
```

Now take any explicit support/source family with constructive ceiling

```math
\zeta_{\max} > 1
```

and admissible blocking window

```math
0 \le \epsilon_{\rm blk} < 1/\zeta_{\max}.
```

Then the worst blocked demand obeys

```math
\zeta_{\rm req}^{\rm(act)}
<
\frac{1}{3-2/\zeta_{\max}},
```

and

```math
\zeta_{\max} - \frac{1}{3-2/\zeta_{\max}}
=
\frac{3\zeta_{\max}(\zeta_{\max}-1)}{3\zeta_{\max}-2} > 0.
```

So any explicit family with `\zeta_{\max}>1` already passes the actual isotropic
support test throughout its admissible blocked regime.

In particular, the explicit Family-1 branch has

```math
\zeta_{\max}^{\rm(F1)} \approx 2.4675 > 1,
```

so its support/source side is automatic once the actual isotropic branch is
adopted.

---

## Step 26E — What actually changes

After this step, the active reduced theorem ledger is much smaller:

- **no** geometry-lane contamination on the isotropic branch,
- **no** extra contact/pole ambiguity,
- **no** remaining explicit support/source ambiguity,
- only the single passive/outgoing normalization defect `N_Q`.

So the old selected twin-support curve is no longer the live reduced object. It
was useful while the support/source side still looked like an independent
selection problem, but once the actual isotropic grouped-`P2` one-pole branch is
used, the support/source theorem becomes automatic.

The real remaining question is now simply:

> **does the completed moving-throat PDE give `N_Q = 1` on the actual passive/
> outgoing isotropic grouped-`P2` branch?**

---

## Main result of the step

On the actual isotropic grouped-`P2` one-pole branch,

```math
\boxed{\epsilon_2 = \epsilon_4 = 0,}
```

so the conservative quadrupole module is exactly

```math
\boxed{
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+
\frac14\frac{1}{1-\omega^2/\Omega_Q^2}.
}
```

The support/source side is then automatic, and the full remaining reduced
mismatch collapses to one scalar normalization defect,

```math
\boxed{N_Q = \bar K_0 / \bar K_0^{\rm target}.}
```

So the next honest move is no longer to keep refining support-side phase
selection. It is to compute the actual passive/outgoing quadrupole normalization
of the moving-throat branch itself.
