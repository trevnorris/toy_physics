# Step 34 — The balanced outlet collapses to a one-parameter parent compensation family

## Goal

Step 33 pushed the surviving compensated outlet into a concrete two-channel
throat core. The next natural move is to eliminate the remaining reduced-core
bookkeeping and rewrite the whole balanced branch directly in terms of the
**parent overlap ratios** already isolated in the moving-throat notes.

That matters because it tells us whether the anomaly is asking for a new core
mechanism, or only for a small detuning of a branch that is already naturally
present in the parent-action variables.

---

## Step 34A — The exact core-balance theorem in parent ratios

Introduce the two normalized parent ratios

```math
\boxed{
\mathfrak r:=\frac{\lambda}{\sqrt{K_sK_q}},
\qquad
\mathfrak g:=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}}.
}
```

Then the exact concrete core-balance law from Step 33,

```math
g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2,
```

collapses to the one-line parent relation

```math
\boxed{
1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2.
}
```

So the compensated-even outlet is no longer a five-parameter tuning problem.
It is a **one-parameter family**: once the normalized static/mixed hybridization
`\mathfrak r` is fixed by the throat background, the required normalized mouth
coupling `\mathfrak g` is fixed exactly.

The two exact branches are

```math
\boxed{
\mathfrak g
=
\mathfrak r
\pm
\frac12\sqrt{1+\mathfrak r^2}.
}
```

That is the parent-action version of the Stage-98 compensation surface.

---

## Step 34B — The side-tube geometry is fixed by the same ratio

On the balanced-even branch the bare mixed coefficient obeys

```math
\kappa_0=\frac{1+\mathfrak r^2}{3}.
```

If the mixed side-channel is realized by the first D/N half-wave on an auxiliary
 tube of length `L_W`, then

```math
\kappa_0=\frac{4L_W^2}{\pi^2 a^2}.
```

So the exact D/N tube-selection law becomes

```math
\boxed{
L_W=
\frac{\pi a}{2}\sqrt{\frac{1+\mathfrak r^2}{3}}.
}
```

So the same normalized hybridization `\mathfrak r` fixes both:

- the required mouth-coupling balance, and
- the auxiliary D/N geometry.

That is a very strong collapse of the surviving outlet freedom.

---

## Step 34C — The loading share simplifies to `\rho_c/4`

Recall the concrete core quantities

```math
\rho_c=\frac{g_s^2}{K_s},
\qquad
\sigma_c=
\frac{(K_s g_q-\lambda g_s)^2}{K_s^2K_q(1+r_c)},
\qquad
r_c=\frac{\lambda^2}{K_sK_q}=\mathfrak r^2.
```

In parent-ratio language,

```math
\boxed{
\sigma_c
=
\rho_c\,
\frac{(\mathfrak g-\mathfrak r)^2}{1+\mathfrak r^2}.
}
```

So on either exact compensation branch,

```math
\boxed{
\sigma_c=\frac{\rho_c}{4}.
}
```

That means the balanced branch does not leave both `\rho_c` and `\sigma_c`
independent. The static loading share is fixed to one quarter of the shell-side
normalization `\rho_c`.

---

## Step 34D — The electron anomaly becomes a tiny odd detuning of the parent family

Let

```math
x:=\Lambda_1 f.
```

The exact electron target from the earlier steps is

```math
\chi_e=\frac{1}{1+x}.
```

On the balanced hybrid outlet,

```math
\chi_Q=\frac{1-9\sigma_c\gamma_c}{1-\sigma_c}.
```

Substituting `\sigma_c=\rho_c/4` and matching to `\chi_e` gives

```math
\boxed{
\gamma_c
=
\frac{\rho_c+4x}{9\rho_c(1+x)}.
}
```

Because `\gamma_0=(1+\mathfrak r^2)\gamma_c`, the corresponding bare mixed odd
coefficient is

```math
\boxed{
\gamma_0
=
\frac{(1+\mathfrak r^2)(\rho_c+4x)}{9\rho_c(1+x)}.
}
```

So the anomaly is now even more sharply localized than in Step 33:

- the compensated-even branch itself is fixed by `\mathfrak r`,
- the static loading share is `\rho_c/4`,
- and the anomaly only asks for a small odd detuning above the canonical values
  ```math
  \gamma_c=\frac19,
  \qquad
  \gamma_0=\frac{1+\mathfrak r^2}{9}.
  ```

Indeed,

```math
\boxed{
\gamma_c-\frac19
=
\frac{x(4-\rho_c)}{9\rho_c(1+x)}.
}
```

So at fixed `\rho_c=O(1)` and tiny `x`, the required shift is necessarily tiny.

---

## Main result of the step

The surviving compensated outlet is no longer a free reduced core at all. It has
collapsed to a one-parameter parent family:

```math
\boxed{
1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2,
\qquad
L_W=\frac{\pi a}{2}\sqrt{\frac{1+\mathfrak r^2}{3}},
\qquad
\sigma_c=\frac{\rho_c}{4}.
}
```

And on that family the electron anomaly requires only

```math
\boxed{
\gamma_c=
\frac{\rho_c+4x}{9\rho_c(1+x)},
\qquad
\gamma_0=
\frac{(1+\mathfrak r^2)(\rho_c+4x)}{9\rho_c(1+x)}.
}
```

So the next real question is no longer “is there a parent family that can host the
outlet?” There is. The sharper question is whether the **actual mouth/source law**
selects a natural point on that family without fine tuning.
