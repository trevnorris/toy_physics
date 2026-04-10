# Step 32 — The first explicit outlet audit leaves only a compensated Robin–mixed branch alive

## Goal

Step 31 turned the quartic g-2 sliver into an exact finite-`f` isotropic DtN
surface. The next honest question is no longer abstract:

> **which explicit low-frequency moving-throat outlet classes can realize that
> surface without spoiling the already-fixed conservative even `l=2`
> fingerprint?**

The moving-throat outlet notes already tested three explicit isotropic classes:

1. a pure geometric Robin core,
2. a standalone mixed `A_w/F_{\mu w}` side-channel pole,
3. a hybrid Robin–mixed outlet.

This step translates those outlet classes directly into the g-2 branch-selection
language.

The main result is sharp:

- a pure Robin core can match the electron outgoing defect exactly, but it
  distorts the even branch and therefore cannot be the whole story;
- a standalone mixed pole cannot preserve the canonical even branch unless it is
  absent;
- the first serious surviving candidate is the **compensated Robin–mixed
  outlet**.

---

## Step 32A — Pure geometric Robin core

The raw isotropic Robin outlet has the exact normalization factor

```math
\boxed{
\chi_Q^{\rm R}=\frac{3}{3-\rho_R},
}
```

where `\rho_R` is the dimensionless isotropic Robin loading.

Matching this directly to the exact electron target

```math
\chi_e=\frac{1}{1+\Lambda_1 f}
```

gives

```math
\boxed{
\rho_R = -3\Lambda_1 f.
}
```

Numerically,

```math
\boxed{
\rho_R^{(e)}
\approx
-9.74211012119238\times 10^{-4}.
}
```

So a pure Robin core can indeed carry the required *size* of the outgoing defect.

But the moving-throat notes also show that the raw Robin outlet changes the
canonical even `l=2` fingerprint at the same time. So it cannot be the full
answer if the conservative grouped-`P2` even branch is already fixed.

That makes the pure Robin result useful only as a scale diagnostic, not as a
finished outlet law.

---

## Step 32B — Standalone mixed side-channel pole

The raw mixed side-channel model is

```math
\Lambda_2^{\rm mix}(z)
=
\Lambda_2^{\rm out}(z)
-\frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5}.
```

For this branch, the canonical even conditions reduce to

```math
-\frac{L_2}{L_0}=\frac19,
\qquad
\frac{L_2^2}{L_0^2}-\frac{L_4}{L_0}=\frac{4}{81}.
```

Solving them gives first

```math
\kappa_W=-\frac19,
```

and then, after substitution,

```math
\boxed{
\sigma_W=0.
}
```

So a standalone passive mixed pole of this type cannot sit on the already-fixed
canonical even branch unless it vanishes.

That is a strong exclusion result:

> the mixed sector may still matter, but not as a naive isolated Schur-complement
> outlet of this simplest type.

---

## Step 32C — Hybrid Robin–mixed outlet

The hybrid outlet is

```math
\Lambda_2^{\rm hyb}(z)
=
\Lambda_2^{\rm out}(z)
+\rho_R
-\frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5}.
```

Imposing the canonical even `l=2` fingerprint yields exactly two branches:

```math
\boxed{
\rho_R=\sigma_W,\qquad \kappa_W=0,
}
```

or

```math
\boxed{
\rho_R=4\sigma_W,\qquad \kappa_W=\frac13.
}
```

The first is only a trivial cancellation branch.

The second is the first **nontrivial compensated branch**.
On that branch the exact normalization factor is

```math
\boxed{
\chi_Q^{\rm hyb}
=
\frac{1-9\sigma_W\gamma_W}{1-\sigma_W}.
}
```

For the canonical outgoing branch itself one recovers

```math
\boxed{
\gamma_W=\frac19
\quad\Longrightarrow\quad
\chi_Q^{\rm hyb}=1.
}
```

So the compensated hybrid outlet reproduces the same harmless pure-scale class
already identified in the Stage-90/91 robustness audit when its odd coefficient
is set to the canonical value.

---

## Step 32D — Exact electron-anomaly family on the compensated branch

Now replace the canonical target `\chi_Q=1` by the exact electron target
`\chi_e=1/(1+\Lambda_1 f)`.

Solving

```math
\frac{1-9\sigma_W\gamma_W}{1-\sigma_W}
=
\frac{1}{1+\Lambda_1 f}
```

gives the exact anomaly family

```math
\boxed{
\gamma_W
=
\frac{\sigma_W+\Lambda_1 f}
{9\sigma_W(1+\Lambda_1 f)}.
}
```

Equivalently,

```math
\boxed{
\gamma_W
=
\frac19
+
\frac{\Lambda_1 f\,(1-\sigma_W)}
{9\sigma_W(1+\Lambda_1 f)}.
}
```

So on the compensated hybrid branch, the electron anomaly is carried by a small
shift of the odd mixed outlet above its canonical value `1/9`.

That is exactly the kind of result we wanted:

- the even branch remains fixed,
- the outlet deformation is not broad or arbitrary,
- and the anomaly rides only in one controlled odd coefficient.

For a representative moderate loading `\sigma_W=1/2`, this gives

```math
\boxed{
\gamma_W \approx 0.111147181287128,
}
```

which is only slightly above

```math
\frac19 \approx 0.111111111111111.
```

So the required odd detuning is tiny, just as the anomaly itself is tiny.

---

## Step 32E — Interpretation of the outlet audit

This is the first really constraining answer after Step 31.

### What is ruled out

- **Pure Robin loading alone** is too blunt: it can carry the defect size but it
  distorts the even branch.
- **A naive standalone mixed pole** is too rigid: it cannot preserve the
  canonical even branch unless it disappears.

### What survives

- The first serious surviving outlet class is the **compensated Robin–mixed
  branch**
  ```math
  \rho_R=4\sigma_W,\qquad \kappa_W=\frac13.
  ```
- On that branch the anomaly is carried entirely by the odd side-channel
  coefficient
  ```math
  \gamma_W=\frac{\sigma_W+\Lambda_1 f}{9\sigma_W(1+\Lambda_1 f)}.
  ```

So the outlet problem is now much narrower than it was before.

It is no longer

> “some deformed isotropic DtN branch.”

It is now

> **an explicit compensated Robin–mixed outlet family with one tiny odd
> coefficient shift left to determine.**

---

## Main result of the step

The first explicit outlet audit leaves only one serious exact candidate class
alive for the g-2 quartic sliver:

```math
\boxed{
\rho_R=4\sigma_W,\qquad \kappa_W=\frac13,\qquad
\gamma_W=\frac{\sigma_W+\Lambda_1 f}{9\sigma_W(1+\Lambda_1 f)}.
}
```

So the next natural move is even sharper:

> **push this compensated outlet down into an actual throat-core model and see
> whether the required odd shift can be written directly in core couplings rather
> than as reduced DtN coefficients.**
