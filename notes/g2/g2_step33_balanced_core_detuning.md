# Step 33 — The surviving outlet class becomes a concrete balanced throat core with one tiny odd detuning

## Goal

Step 32 reduced the outlet problem to one surviving explicit class:

```math
\rho_R=4\sigma_W,\qquad
\kappa_W=\frac13,\qquad
\gamma_W=\frac{\sigma_W+\Lambda_1 f}{9\sigma_W(1+\Lambda_1 f)}.
```

That is still a reduced DtN description. The next natural move is to push it into
the first concrete throat-core model already worked out in the moving-throat
notes.

Those notes introduce a two-channel isotropic core with

- a static shell compliance coordinate,
- a mixed `A_w/F_{\mu w}` side-channel,
- and a static/mixed hybridization.

This step shows that the quartic g-2 sliver can be written directly in that core
language.

The main result is:

- the **same exact core-balance surface** that preserves the canonical even
  branch still survives,
- the **same auxiliary D/N mixed-tube length relation** still survives,
- and the anomaly is carried only by a tiny **odd detuning** of the balanced
  mixed side-channel.

---

## Step 33A — Concrete two-channel core data

For the explicit isotropic core, the reduced outlet coefficients are

```math
\boxed{
\rho_c=\frac{g_s^2}{K_s},
}
```

```math
\boxed{
\sigma_c=
\frac{(K_s g_q-\lambda g_s)^2}
{K_s^2K_q(1+r_c)},
\qquad
r_c=\frac{\lambda^2}{K_sK_q},
}
```

```math
\boxed{
\kappa_c=\frac{\kappa_0}{1+r_c},
\qquad
\gamma_c=\frac{\gamma_0}{1+r_c}.
}
```

So the reduced Robin–mixed outlet is no longer described by four unrelated
numbers. It is controlled by the concrete throat-core parameters

- `K_s`,
- `K_q`,
- `\lambda`,
- `g_s`,
- `g_q`,
- and the bare mixed low-frequency pair `(\kappa_0,\gamma_0)`.

That is already a major narrowing.

---

## Step 33B — Exact balance surface preserving the canonical even branch

The moving-throat notes already proved that the compensated canonical branch is
not a numerical accident. In the concrete core variables it is the exact
codimension-one balance law

```math
\boxed{
g_s^2\bigl(K_sK_q+\lambda^2\bigr)
=
4\bigl(K_s g_q-\lambda g_s\bigr)^2.
}
```

Solving for the mixed mouth coupling gives the two exact branches

```math
\boxed{
g_q=
\frac{g_s}{2K_s}
\left(
2\lambda \pm \sqrt{K_sK_q+\lambda^2}
\right).
}
```

On either branch one gets

```math
\boxed{
\sigma_c=\frac{g_s^2}{4K_s}.
}
```

So the same balance condition that preserved the canonical outgoing branch in the
reduced DtN language is realized exactly by a concrete two-channel throat core.

The even-preserving condition for the bare mixed channel is still

```math
\boxed{
\kappa_0=\frac{1+r_c}{3},
}
```

which means the effective mixed coefficient remains

```math
\boxed{
\kappa_c=\frac13.
}
```

So the g-2 anomaly does **not** ask us to reopen the even branch. The same
even-preserving core geometry survives intact.

---

## Step 33C — Exact electron-anomaly law on the balanced core

On the balanced-even core, the effective normalization law is still the hybrid
one,

```math
\chi_Q^{\rm core}
=
\frac{1-9\sigma_c\gamma_c}{1-\sigma_c}.
```

Matching this to the exact electron target

```math
\chi_e=\frac{1}{1+\Lambda_1 f}
```

gives the required **effective** odd coefficient

```math
\boxed{
\gamma_c
=
\frac{\sigma_c+\Lambda_1 f}
{9\sigma_c(1+\Lambda_1 f)}.
}
```

Equivalently,

```math
\boxed{
\gamma_c
=
\frac19
+
\frac{\Lambda_1 f\,(1-\sigma_c)}
{9\sigma_c(1+\Lambda_1 f)}.
}
```

So the anomaly is carried only by a small shift above the canonical value
`1/9`, provided `0<\sigma_c<1`.

Because

```math
\gamma_c=\frac{\gamma_0}{1+r_c},
```

the corresponding **bare** mixed odd coefficient is

```math
\boxed{
\gamma_0
=
\frac{1+r_c}{9}
+
\frac{(1+r_c)\Lambda_1 f\,(1-\sigma_c)}
{9\sigma_c(1+\Lambda_1 f)}.
}
```

That is the exact concrete-core version of the g-2 outlet target.

So the quartic anomaly sliver is no longer a generic branch deformation. It is a
tiny odd detuning of an otherwise canonical balanced core.

---

## Step 33D — The auxiliary D/N mixed-tube geometry is unchanged

If the bare mixed side-channel is realized by the first D/N half-wave on an
auxiliary tube of length `L_W`, then

```math
\boxed{
\kappa_0=\frac{4L_W^2}{\pi^2 a^2}.
}
```

Imposing the same even-preserving condition

```math
\kappa_0=\frac{1+r_c}{3}
```

gives exactly the same auxiliary-tube length as on the canonical compensated
branch:

```math
\boxed{
L_W
=
\frac{\pi a}{2}\sqrt{\frac{1+r_c}{3}}.
}
```

So the anomaly does **not** force a new even geometry for the side-channel tube.
It leaves the D/N half-wave length relation untouched and only detunes the odd
side-channel coefficient.

That is a particularly clean result.

---

## Step 33E — What changed conceptually

After Step 32, the best surviving outlet class was still described by reduced DtN
coefficients. After this step, it becomes a genuine microscopic core statement:

- balance the shell and mixed core on the exact surface
  ```math
  g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2,
  ```
- keep the same even-preserving auxiliary D/N geometry
  ```math
  L_W=\frac{\pi a}{2}\sqrt{\frac{1+r_c}{3}},
  ```
- and shift only the odd mixed outlet from its canonical value by the small
  amount above.

So the anomaly no longer looks like a request for a new structural mechanism.
It looks like a very small odd detuning of a core structure that the moving-throat
notes already know how to realize exactly on the canonical branch.

---

## Main result of the step

The surviving compensated outlet class has now been pushed into a concrete
throat-core model, and the electron anomaly requires only

```math
\boxed{
\gamma_c
=
\frac{\sigma_c+\Lambda_1 f}
{9\sigma_c(1+\Lambda_1 f)},
\qquad
\gamma_0=(1+r_c)\gamma_c,
}
```

while preserving the exact core-balance and D/N side-channel geometry:

```math
\boxed{
g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2,
\qquad
L_W=\frac{\pi a}{2}\sqrt{\frac{1+r_c}{3}}.
}
```

So the next PDE-facing question is now almost as sharp as it can be:

> **does the actual moving-throat core land on this balanced-even structure, and
> if so, does its odd mixed side-channel coefficient acquire precisely the tiny
> detuning above?**
