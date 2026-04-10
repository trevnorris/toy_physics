# Step 23 — Selected-branch loading ratio from the minimal isotropic quadrupole precursor

## Goal

Step 22 reduced the last support-selection ambiguity to the scalar demand parameter

```math
\varrho := \frac{\pi^2\Pi_{\rm tr}}{16\Lambda}.
```

But that still left one question open:

> **what does the selected passive/outgoing normalization side actually choose for
> `\Pi_{\rm tr}`?**

The later moving-throat notes answer that in two stages.

1. First, the selected-branch quadrupole-demand product **cancels** all separate
   outgoing-normalization amplitudes and depends only on the loading ratio
   ```math
   \rho_\alpha := \frac{\alpha_{\rm req}}{\alpha_{\rm mix}}.
   ```
2. Second, the natural minimal isotropic conservative quadrupole precursor fixes
   that loading ratio exactly through its contact/pole split.

So Step 23 is the point where the support selector `\varrho` is finally tied back
into the selected-branch normalization side.

The main result is:

```math
\boxed{
\rho_\alpha = \frac43,
\qquad
\zeta_{\rm req}=\frac13,
\qquad
\Pi_{\rm tr}=\frac43\,C_{\rm mix}.
}
```

Equivalently,

```math
\boxed{
\varrho = \frac{2(1-\epsilon_*)}{3},
\qquad
S_{\rm req}=\frac{\Pi_{\rm tr}}{C_{\rm mix}}=\frac43.
}
```

So the natural minimal isotropic passive/outgoing branch is **not** mixed-only and
**not** non-twin. It sits exactly on the symmetric-lowest-twin support slice.

---

## Inputs from the selected-branch normalization side

The later moving-throat notes give the exact product identities

```math
\Pi_{\rm tr}
=
\frac{N_Q^{(\rm target)}}{\beta_0}\,\alpha_{\rm req},
\qquad
C_{\rm mix}
=
\frac{N_Q^{(\rm target)}}{\beta_0}\,\alpha_{\rm mix}.
```

So immediately

```math
\boxed{
\frac{\Pi_{\rm tr}}{C_{\rm mix}}
=
\frac{\alpha_{\rm req}}{\alpha_{\rm mix}}
=:\rho_\alpha.
}
```

This is the key simplification: once the outgoing quadrupole branch is normalized,
all the separate selected-mode amplitudes drop out of the support test.

In the spectral notation of the selected branch,

```math
N_Q^{(\rm target)} = \hat m_-^2\,\beta_0\,\frac{s_-}{\lambda_-},
```

so the same identities can also be written as

```math
\Pi_{\rm tr} = \hat m_-^2\frac{s_-}{\lambda_-}\alpha_{\rm req},
\qquad
C_{\rm mix} = \hat m_-^2\frac{s_-}{\lambda_-}\alpha_{\rm mix}.
```

Again the ratio is just `\rho_\alpha`.

---

## Step 23A — Exact contact-plus-pole inverse formulas

The natural minimal conservative quadrupole precursor is written as

```math
Y_Q^{\rm cons}(\omega)
=
 c_0 + \frac{c_1}{1-\omega^2/\Omega_Q^2},
```

with normalized static limit

```math
c_0 + c_1 = 1.
```

On the explicit support/source branch, the natural reading is:

- the mixed baseline carries the static contact fraction,
- the extra support lane carries the finite conservative pole.

So the same precursor can be written as

```math
Y_Q^{\rm cons}(\omega)
=
\frac{\alpha_{\rm mix}}{\alpha_{\rm req}}
+
\frac{\alpha_{\rm req}-\alpha_{\rm mix}}{\alpha_{\rm req}}
\frac{1}{1-\omega^2/\Omega_Q^2}.
```

Introducing

```math
\rho_\alpha := \frac{\alpha_{\rm req}}{\alpha_{\rm mix}},
```

this becomes

```math
Y_Q^{\rm cons}(\omega)
=
\frac{1}{\rho_\alpha}
+
\frac{\rho_\alpha-1}{\rho_\alpha}
\frac{1}{1-\omega^2/\Omega_Q^2}.
```

Therefore the exact contact/pole data are

```math
\boxed{c_0 = \frac{1}{\rho_\alpha}},
\qquad
\boxed{c_1 = \frac{\rho_\alpha-1}{\rho_\alpha}},
```

and the inverse formulas are

```math
\boxed{\rho_\alpha = \frac{1}{c_0} = \frac{1}{1-c_1}},
```

```math
\boxed{\zeta_{\rm req} := \rho_\alpha-1 = \frac{c_1}{c_0}.}
```

So the support/source loading ratio is directly encoded in the static contact /
pole split of the conservative quadrupole precursor.

---

## Step 23B — Matching to the minimal isotropic quadrupole module

The 2.5PN quadrupole audit already fixed the smallest viable isotropic
conservative precursor to

```math
c_0 = \frac34,
\qquad
c_1 = \frac14,
\qquad
\Omega_Q = \frac{3c_s}{2a}.
```

Inserting those values into the exact inverse formulas gives immediately

```math
\boxed{\rho_\alpha = \frac{1}{3/4} = \frac43,}
```

```math
\boxed{\zeta_{\rm req} = \frac{1/4}{3/4} = \frac13.}
```

Then the selected demand product is

```math
\boxed{
\Pi_{\rm tr}
=
\rho_\alpha C_{\rm mix}
=
\frac43 C_{\rm mix}.
}
```

So the selected branch is no longer carrying an arbitrary support demand. The
natural minimal isotropic passive/outgoing branch fixes it exactly.

---

## Step 23C — Exact support-selector form of the selected branch

Step 22 defined

```math
\varrho := \frac{\pi^2\Pi_{\rm tr}}{16\Lambda},
\qquad
C_{\rm mix} = \frac{8\Lambda(1-\epsilon_*)}{\pi^2}.
```

Substituting

```math
\Pi_{\rm tr} = \frac43 C_{\rm mix}
```

gives

```math
\varrho
=
\frac{\pi^2}{16\Lambda}\cdot\frac43\cdot\frac{8\Lambda(1-\epsilon_*)}{\pi^2}
=
\boxed{\frac{2(1-\epsilon_*)}{3}.}
```

And the required support enhancement is simply

```math
S_{\rm req}
=
\frac{\Pi_{\rm tr}}{C_{\rm mix}}
=
\boxed{\frac43.}
```

So the selected branch is no longer scanning all support-demand sectors. It is
locked to one exact support ratio.

---

## Step 23D — Regime meaning

Stage 22 already split the support regimes by

```math
\Pi_{\rm tr} \le C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{mixed-only enough},
```

```math
C_{\rm mix} < \Pi_{\rm tr} \le 2C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{symmetric lowest twin enough},
```

```math
\Pi_{\rm tr} > 2C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{non-twin asymmetry required}.
```

Because the selected branch gives

```math
\Pi_{\rm tr} = \frac43 C_{\rm mix},
```

it follows exactly that

```math
\boxed{
C_{\rm mix} < \Pi_{\rm tr} < 2C_{\rm mix}.
}
```

So:

- mixed-only is **not** enough,
- the symmetric lowest twin **is** enough,
- and non-twin asymmetry is **not** required.

This is already a real simplification of the anomaly bridge.

---

## Main result of the step

The selected-branch normalization side has now fixed the support ratio carried by
the natural minimal isotropic passive/outgoing quadrupole branch:

```math
\boxed{
\rho_\alpha = \frac43,
\qquad
\zeta_{\rm req}=\frac13,
\qquad
\Pi_{\rm tr}=\frac43 C_{\rm mix}.
}
```

Equivalently,

```math
\boxed{
\varrho = \frac{2(1-\epsilon_*)}{3},
\qquad
S_{\rm req}=\frac43.
}
```

So the last support ambiguity has collapsed from three sectors

- mixed-only,
- symmetric lowest twin,
- non-twin asymmetry,

to exactly **one** selected support slice:

```math
\text{symmetric lowest twin, with demand ratio } \Pi_{\rm tr}/C_{\rm mix}=4/3.
```

---

## What the next step should be

The next honest move is now very sharp:

> restrict the Step-21 primitive quartic ranking problem to this selected
> twin-support branch, and ask how much of the remaining `q_W` versus
> `q_\Lambda` ambiguity survives there.

That is the smallest next derivation that still genuinely pushes the anomaly
closure forward.
