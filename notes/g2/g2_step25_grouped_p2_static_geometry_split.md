# Step 25 — The grouped-`P2` + static-geometry split forces the `3/4 + 1/4` conservative module

## Goal

Step 24 left the quartic anomaly hierarchy on the selected twin-support curve
parameterized by the single selector `\varrho`. That was the right support-side
picture, but it still treated the conservative passive/outgoing quadrupole branch
as if its contact/pole split were an external choice.

The moving-throat notes give a sharper route.

The 3PN conservative split already says that the higher conservative payload is
organized as

```math
\text{grouped real }P_2 \text{ middle block} + \text{static geometry completion},
```

and the 2.5PN quadrupole audit already fixed the minimal isotropic branch
identity

```math
K_0 K_4 = 4 K_2^2.
```

So the next honest question is:

> **if the actual conservative branch is one isotropic grouped-`P2` pole plus a
> purely static geometry completion, what contact/pole split is forced?**

This step answers that exactly.

The main result is

```math
\boxed{K_{\rm geom} = 3 K_{\rm pole}},
```

so

```math
\boxed{K_{\rm pole}=K_0/4,
\qquad
K_{\rm geom}=3K_0/4,}
```

and therefore the normalized conservative quadrupole module is forced to be

```math
\boxed{
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+
\frac{1}{4}
\frac{1}{1-\omega^2/\Omega_Q^2}.
}
```

Once that is true, the support/source loading ratio is no longer free:

```math
\boxed{\rho_\alpha = \frac43,
\qquad
\zeta_{\rm req}=\frac13,
\qquad
\Pi_{\rm tr}=\frac43 C_{\rm mix}.}
```

So this is the step where the old support-side selector stops being the real
bottleneck. The support ratio `4/3` follows as a **corollary** of the grouped
conservative branch organization.

---

## Step 25A — Minimal isotropic grouped-`P2` + geometry realization

Take the conservative isotropic quadrupole module in the minimal form

```math
K_Q^{\rm cons}(\omega)
=
K_{\rm geom}
+
\frac{K_{\rm pole}}{1-\omega^2/\Omega_Q^2}.
```

Here:

- `K_geom` is the static geometry completion,
- `K_pole` is the isotropic grouped-`P2` pole residue,
- `\Omega_Q` is the effective isotropic grouped-`P2` pole.

This is the smallest realization compatible with the already-frozen 3PN split if

- the grouped-`P2` side is the only dynamic conservative quadrupole lane,
- and geometry contributes only the static completion.

Expand through `O(\omega^4)`:

```math
K_Q^{\rm cons}(\omega)
=
K_0 + K_2 \omega^2 + K_4 \omega^4 + O(\omega^6),
```

with exact coefficients

```math
K_0 = K_{\rm geom}+K_{\rm pole},
```

```math
K_2 = \frac{K_{\rm pole}}{\Omega_Q^2},
```

```math
K_4 = \frac{K_{\rm pole}}{\Omega_Q^4}.
```

So the whole conservative branch is now parameterized by two amplitudes and one
pole scale.

---

## Step 25B — The minimal isotropic branch identity fixes the contact/pole split

The 2.5PN audit already froze the minimal isotropic conservative branch identity

```math
K_0 K_4 = 4 K_2^2.
```

Insert the one-pole + static-geometry coefficients:

```math
(K_{\rm geom}+K_{\rm pole})\frac{K_{\rm pole}}{\Omega_Q^4}
=
4\left(\frac{K_{\rm pole}}{\Omega_Q^2}\right)^2.
```

Assuming the branch is nontrivial (`K_pole\neq0`), the common factor cancels and
we get

```math
K_{\rm geom}+K_{\rm pole}=4K_{\rm pole},
```

hence

```math
\boxed{K_{\rm geom}=3K_{\rm pole}.}
```

So equivalently

```math
\boxed{K_{\rm pole}=\frac{K_0}{4},
\qquad
K_{\rm geom}=\frac{3K_0}{4}.}
```

This is the central algebraic result of the step.

The normalized conservative response therefore becomes

```math
\widehat Y_Q^{\rm cons}(\omega)
:=
\frac{K_Q^{\rm cons}(\omega)}{K_0}
=
\frac34+
\frac14\frac{1}{1-\omega^2/\Omega_Q^2}.
```

So the familiar `3/4 + 1/4` split is not an imported guess anymore. Under the
minimal grouped-`P2` + static-geometry realization, it is forced.

---

## Step 25C — Exact contact/pole map to the support/source loading ratio

On the explicit support/source branch, the natural contact-plus-pole reading is

```math
Y_Q^{\rm cons}(\omega)
=
\frac{\alpha_{\rm mix}}{\alpha_{\rm req}}
+
\frac{\alpha_{\rm req}-\alpha_{\rm mix}}{\alpha_{\rm req}}
\frac{1}{1-\omega^2/\Omega_Q^2}.
```

Introduce the selected loading ratio

```math
\rho_\alpha := \frac{\alpha_{\rm req}}{\alpha_{\rm mix}}.
```

Then

```math
Y_Q^{\rm cons}(\omega)
=
\frac{1}{\rho_\alpha}
+
\frac{\rho_\alpha-1}{\rho_\alpha}
\frac{1}{1-\omega^2/\Omega_Q^2}.
```

So the exact contact and pole fractions are

```math
c_0 = \frac{1}{\rho_\alpha},
\qquad
c_1 = \frac{\rho_\alpha-1}{\rho_\alpha},
```

with inverse formulas

```math
\rho_\alpha = \frac{1}{c_0} = \frac{1}{1-c_1},
```

```math
\zeta_{\rm req} := \rho_\alpha-1 = \frac{c_1}{c_0}.
```

Now insert the forced grouped-branch values

```math
c_0 = \frac34,
\qquad
c_1 = \frac14.
```

Then immediately

```math
\boxed{\rho_\alpha = \frac{1}{3/4} = \frac43,}
```

```math
\boxed{\zeta_{\rm req} = \frac{1/4}{3/4} = \frac13.}
```

And since

```math
\frac{\Pi_{\rm tr}}{C_{\rm mix}} = \rho_\alpha,
```

we also get

```math
\boxed{\Pi_{\rm tr} = \frac43 C_{\rm mix}.}
```

So the Stage-23 support ratio is now recovered directly from the grouped
conservative branch organization.

---

## Step 25D — What actually changes

Before this step, the selected twin-support curve looked like the main organizing
object. After this step, the logic reverses:

1. the higher conservative payload is read as
   ```math
   \text{grouped }P_2 + \text{static geometry},
   ```
2. the minimal isotropic branch identity then forces
   ```math
   3/4 + 1/4,
   ```
3. and only **then** does the support/source side inherit
   ```math
   \rho_\alpha = 4/3.
   ```

So the old support-phase picture is no longer the deepest one. The deeper theorem
gate is now:

> **does the actual moving-throat branch really realize one isotropic grouped-`P2`
> pole plus a purely static geometry completion?**

If yes, the support/source verdict follows immediately.

---

## Main result of the step

The minimal isotropic conservative quadrupole module is forced by the grouped
real `P2` + static-geometry split:

```math
\boxed{
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+
\frac14\frac{1}{1-\omega^2/\Omega_Q^2}.
}
```

And the support/source loading ratio follows as a corollary:

```math
\boxed{\rho_\alpha = \frac43,
\qquad
\zeta_{\rm req}=\frac13,
\qquad
\Pi_{\rm tr}=\frac43 C_{\rm mix}.}
```

So the next honest move is no longer to keep scanning support-sector phase
regions. It is to test whether the actual isotropic moving-throat branch really
has that minimal grouped-`P2` + static-geometry conservative structure.
